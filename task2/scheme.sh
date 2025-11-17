#!/bin/bash

# Компиляция программы
echo "Компиляция MPI программы..."
mpicc -o task2 task2.c -lm

# Создание файла с заголовками
echo "p,time_rows,time_cols,time_blocks" > raw_times.csv

# Запуск на разных количествах процессов
for processes in 1 4 9 16; do
    echo "Запуск на $processes процессах..."
    mpiexec -n $processes --oversubscribe ./task2
done

# Анализ результатов без создания временных файлов
echo "p,time_rows,speedup_rows,eff_rows,time_cols,speedup_cols,eff_cols,time_blocks,speedup_blocks,eff_blocks" > results.csv

# Чтение данных и расчет метрик
{
    read header # пропускаем заголовок
    first_line=true
    base_rows=""
    base_cols=""
    base_blocks=""
    
    while IFS=',' read -r p time_r time_c time_b; do
        if [ "$first_line" = true ]; then
            base_rows=$time_r
            base_cols=$time_c
            base_blocks=$time_b
            first_line=false
            
            # Для p=1 ускорение и эффективность всегда 1 и 100%
            echo "1,$time_r,1.0000,100.00,$time_c,1.0000,100.00,$time_b,1.0000,100.00" >> results.csv
        else
            # Расчет ускорения для строк
            speedup_r=$(echo "scale=4; $base_rows / $time_r" | bc 2>/dev/null || echo "0.0000")
            
            # Расчет эффективности для строк
            efficiency_r=$(echo "scale=2; $speedup_r / $p * 100" | bc 2>/dev/null || echo "0.00")
            
            # Расчет ускорения для столбцов
            speedup_c=$(echo "scale=4; $base_cols / $time_c" | bc 2>/dev/null || echo "0.0000")
            
            # Расчет эффективности для столбцов
            efficiency_c=$(echo "scale=2; $speedup_c / $p * 100" | bc 2>/dev/null || echo "0.00")
            
            # Расчет ускорения для блоков
            speedup_b=$(echo "scale=4; $base_blocks / $time_b" | bc 2>/dev/null || echo "0.0000")
            
            # Расчет эффективности для блоков
            efficiency_b=$(echo "scale=2; $speedup_b / $p * 100" | bc 2>/dev/null || echo "0.00")
            
            # Запись результатов
            echo "$p,$time_r,$speedup_r,$efficiency_r,$time_c,$speedup_c,$efficiency_c,$time_b,$speedup_b,$efficiency_b" >> results.csv
        fi
    done
} < raw_times.csv

echo "Анализ завершен. Результаты сохранены в results.csv"

# Вывод сводки
echo ""
echo "Сводка результатов:"
echo "==================="
cat results.csv
