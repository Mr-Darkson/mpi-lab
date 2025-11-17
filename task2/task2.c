#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

double random_range(double min, double max) {
    return min + ((double)rand() / RAND_MAX) * (max - min);
}

void generate_matrix_vector(double *A, double *x, int N) {
    for (int i = 0; i < N; i++) {
        x[i] = random_range(-5.0, 5.0);
        for (int j = 0; j < N; j++) {
            A[i * N + j] = random_range(-5.0, 5.0);
        }
    }
}

void multiply_local_rows(double *A, double *x, double *y, int local_rows, int N) {
    for (int i = 0; i < local_rows; i++) {
        y[i] = 0.0;
        for (int j = 0; j < N; j++)
            y[i] += A[i * N + j] * x[j];
    }
}

int check_result(double *A, double *x, double *y, int N, int my_rank) {
    if (my_rank != 0) return 1;
    
    double tolerance = 1e-10;
    int correct = 1;
    
    for (int i = 0; i < N; i++) {
        double correct_val = 0.0;
        for (int j = 0; j < N; j++) {
            correct_val += A[i * N + j] * x[j];
        }
        
        double diff = fabs(correct_val - y[i]);
        if (diff > tolerance) {
            printf("ERROR at y[%d]: expected %.10f, got %.10f, diff = %.2e\n", 
                   i, correct_val, y[i], diff);
            correct = 0;
            break;
        }
    }
    
    if (!correct) {
        printf("✗ RESULT INCORRECT\n");
    }
    
    return correct;
}

double multiply_rows_mpi(double *A, double *x, double *y, int N, int my_rank, int comm_sz) {
    int local_rows = N / comm_sz;
    double *local_A = (double *)malloc(local_rows * N * sizeof(double));
    double *local_y = (double *)malloc(local_rows * sizeof(double));

    double start_time, end_time;
    
    MPI_Barrier(MPI_COMM_WORLD);
    start_time = MPI_Wtime();

    MPI_Scatter(A, local_rows * N, MPI_DOUBLE, local_A, local_rows * N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(x, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    multiply_local_rows(local_A, x, local_y, local_rows, N);

    MPI_Gather(local_y, local_rows, MPI_DOUBLE, y, local_rows, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    end_time = MPI_Wtime();

    double global_start_time, global_end_time;
    MPI_Reduce(&start_time, &global_start_time, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(&end_time, &global_end_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    free(local_A);
    free(local_y);

    if (my_rank == 0) {
        check_result(A, x, y, N, my_rank);
    }
    return global_end_time - global_start_time;
}

double multiply_cols_mpi(double *A, double *x, double *y, int N, int my_rank, int comm_sz) {
    if (N % comm_sz != 0) {
        if (my_rank == 0) {
            printf("Error: N must be divisible by comm_sz\n");
        }
        return -1.0;
    }
    
    int local_cols = N / comm_sz;
    double *local_x = (double *)malloc(local_cols * sizeof(double));
    double *local_y = (double *)calloc(N, sizeof(double));
    double *local_A = (double *)malloc(N * local_cols * sizeof(double));

    if (local_x == NULL || local_y == NULL || local_A == NULL) {
        printf("Memory allocation failed on process %d\n", my_rank);
        return -1.0;
    }

    double start_time, end_time;
    MPI_Barrier(MPI_COMM_WORLD);
    start_time = MPI_Wtime();

    MPI_Scatter(x, local_cols, MPI_DOUBLE, local_x, local_cols, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (my_rank == 0) {
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < local_cols; j++) {
                local_A[i * local_cols + j] = A[i * N + j];
            }
        }
        
        for (int proc = 1; proc < comm_sz; proc++) {
            double *send_buffer = (double *)malloc(N * local_cols * sizeof(double));
            for (int i = 0; i < N; i++) {
                for (int j = 0; j < local_cols; j++) {
                    send_buffer[i * local_cols + j] = A[i * N + proc * local_cols + j];
                }
            }
            MPI_Send(send_buffer, N * local_cols, MPI_DOUBLE, proc, 0, MPI_COMM_WORLD);
            free(send_buffer);
        }
    } else {
        MPI_Recv(local_A, N * local_cols, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < local_cols; j++) {
            local_y[i] += local_A[i * local_cols + j] * local_x[j];
        }
    }

    MPI_Reduce(local_y, y, N, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    end_time = MPI_Wtime();

    double global_time;
    double local_time = end_time - start_time;
    double global_start_time, global_end_time;
    MPI_Reduce(&start_time, &global_start_time, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(&end_time, &global_end_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    free(local_A);
    free(local_x);
    free(local_y);

    if (my_rank == 0) {
        check_result(A, x, y, N, my_rank);
    }
    
    return global_end_time - global_start_time;
}
    

double multiply_blocks_mpi(double *A, double *x, double *y, int N, int my_rank, int comm_sz) {
    int q = (int)sqrt(comm_sz);
    if (q * q != comm_sz) {
        if (my_rank == 0) printf("Ошибка: число процессов должно быть квадратом для блочного разбиения!\n");
        return -1.0;
    }

    int block_size = N / q;
    int row_block = my_rank / q;
    int col_block = my_rank % q;

    double *local_A = (double *)malloc(block_size * block_size * sizeof(double));
    double *local_y = (double *)calloc(block_size, sizeof(double));

    double start_time, end_time;
    
    MPI_Barrier(MPI_COMM_WORLD);
    start_time = MPI_Wtime();

    MPI_Bcast(x, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(A, N * N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    for (int i = 0; i < block_size; i++) {
        for (int j = 0; j < block_size; j++) {
            int global_i = row_block * block_size + i;
            int global_j = col_block * block_size + j;
            local_A[i * block_size + j] = A[global_i * N + global_j];
        }
    }

    for (int i = 0; i < block_size; i++) {
        for (int j = 0; j < block_size; j++) {
            int global_j = col_block * block_size + j;
            local_y[i] += local_A[i * block_size + j] * x[global_j];
        }
    }

    double *all_blocks = NULL;
    if (my_rank == 0) {
        all_blocks = (double *)malloc(comm_sz * block_size * sizeof(double));
    }

    MPI_Gather(local_y, block_size, MPI_DOUBLE, 
               all_blocks, block_size, MPI_DOUBLE, 
               0, MPI_COMM_WORLD);

    if (my_rank == 0) {
        for (int i = 0; i < N; i++) {
            y[i] = 0.0;
        }
        
        for (int proc = 0; proc < comm_sz; proc++) {
            int proc_row = proc / q;
            int proc_col = proc % q;
            
            for (int i = 0; i < block_size; i++) {
                int global_i = proc_row * block_size + i;
                y[global_i] += all_blocks[proc * block_size + i];
            }
        }
    
        
        free(all_blocks);
    }

    end_time = MPI_Wtime();

    double global_start_time, global_end_time;
    MPI_Reduce(&start_time, &global_start_time, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(&end_time, &global_end_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    free(local_A);
    free(local_y);
    
    if (my_rank == 0) {
        check_result(A, x, y, N, my_rank);
    }
    return global_end_time - global_start_time;
}

int main(int argc, char **argv) {
    int my_rank, comm_sz;
    int N = 576;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
    srand(my_rank);

    double *A = NULL, *x = NULL, *y = NULL;

    
    A = (double *)malloc(N * N * sizeof(double));
    x = (double *)malloc(N * sizeof(double));
    y = (double *)malloc(N * sizeof(double));

    if (my_rank == 0) {
        generate_matrix_vector(A, x, N);
    } 

    double time_rows = multiply_rows_mpi(A, x, y, N, my_rank, comm_sz);
    double time_cols = multiply_cols_mpi(A, x, y, N, my_rank, comm_sz);
    double time_blocks = multiply_blocks_mpi(A, x, y, N, my_rank, comm_sz);

    if (my_rank == 0) {
        printf("%d,%.6f,%.6f,%.6f\n", comm_sz, time_rows, time_cols, time_blocks);
    
        FILE *csv_file = fopen("raw_times.csv", "a");
        fprintf(csv_file, "%d,%.6f,%.6f,%.6f\n", comm_sz, time_rows, time_cols, time_blocks);
        fclose(csv_file);
    }

    free(x);
    free(y);
    if (my_rank == 0) free(A);

    MPI_Finalize();
    return 0;
}
