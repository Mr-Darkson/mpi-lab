#!/bin/bash

PROGRAM="task1"
RESULTS_FILE="task1_results.csv"

# Проверяем наличие bc
if ! command -v bc &> /dev/null; then
    echo "Ошибка: bc не установлен. Установите: sudo apt-get install bc"
    exit 1
fi

mpicc task1.c -o $PROGRAM -lm

echo "Processes,Time,Speedup,Efficiency" > $RESULTS_FILE

echo "=== Запуск на 1 процессе ==="
time1=$(mpiexec -n 1 ./$PROGRAM | tr -d '\r' | tr -d '[:cntrl:]')
echo "Время: $time1 секунд"
echo "1,$time1,1.00,100.00" >> $RESULTS_FILE
echo ""

# Запускаем только до количества доступных ядер
max_processes=8
for ((n=2; n<=max_processes; n++)); do
    echo "=== Запуск на $n процессах ==="
    time_n=$(mpiexec --oversubscribe -n $n ./$PROGRAM | tr -d '\r' | tr -d '[:cntrl:]')
    
    # Проверяем, что time_n - корректное число
    if [[ $time_n =~ ^[0-9]+\.?[0-9]*$ ]] && [ $(echo "$time_n > 0" | bc) -eq 1 ]; then
        speedup=$(echo "scale=2; $time1 / $time_n" | bc -l)
        efficiency=$(echo "scale=2; $speedup / $n * 100" | bc -l)
        echo "$n,$time_n,$speedup,$efficiency" >> $RESULTS_FILE
        echo "Время: $time_n секунд"
        echo "Ускорение: $speedup"
        echo "Эффективность: $efficiency%"
    else
        echo "Ошибка: некорректное время выполнения: $time_n"
        echo "$n,ERROR,0,0" >> $RESULTS_FILE
    fi
    echo ""
done
