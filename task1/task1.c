#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

double random_range(double min, double max) {
    return min + ((double)rand() / RAND_MAX) * (max - min);
}

long long calc_num_points_in_circle(long long local_num_points, int my_rank) {
    unsigned int seed = time(NULL) + my_rank;
    long long local_points_in_circle = 0;
    
    for (long long i = 0; i < local_num_points; i++) {
        double x = random_range(-1.0, 1.0);
        double y = random_range(-1.0, 1.0);
        
        if (x*x + y*y <= 1.0) {
            local_points_in_circle++;
        }
    }
    
    return local_points_in_circle;
}

int main() {
    int comm_sz, my_rank;
    long long total_num_points = 100000000;
    double start_time, end_time;
    
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    
    long long num_points_per_process = total_num_points / comm_sz;
    
    // Синхронизация перед началом измерения
    MPI_Barrier(MPI_COMM_WORLD);
    start_time = MPI_Wtime();
    
    long long points_in_circle_per_process = calc_num_points_in_circle(num_points_per_process, my_rank);
    
    long long total_points_in_circle;
    MPI_Reduce(&points_in_circle_per_process, &total_points_in_circle, 1, MPI_LONG_LONG, 
               MPI_SUM, 0, MPI_COMM_WORLD);
    
    end_time = MPI_Wtime();
    double local_execution_time = end_time - start_time;
    
    double total_execution_time;
    MPI_Reduce(&local_execution_time, &total_execution_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    
    if (my_rank == 0) {
        printf("%.6f\n", total_execution_time);
        fflush(stdout);  // Принудительно сбрасываем буфер вывода
    }
    
    MPI_Finalize();
    return 0;
}
