#!/bin/sh

: <<'DOC'
Questo script esegue un programma sfruttando OpenMP con
un numero di core variabile; ogni esecuzione considera
sempre la stessa dimensione dell'input, quindi i tempi misurati
possono essere usati per calcolare speedup e strong scaling efficiency.
DOC

if [ $# != 2 ]; then
echo "Usage:
      perf-eval [executable] [nbr-tries]"
else

CORES=$(grep -c processor < /proc/cpuinfo)
SCRIPT_NAME=$1
TRIES=$2
SCRIPT_DATA=$SCRIPT_NAME".csv"

printf '"Threads", "Time in seconds"\n' >> "${SCRIPT_DATA}"

for p in $(seq "$CORES"); do
    for _ in $(seq "$TRIES"); do
        EXEC_TIME="$( OMP_NUM_THREADS=$p ./omp-skyline < test6-N100000-D50.in )"
        EXEC_TIME=$(echo "$EXEC_TIME" | sed 's/Execution time: //')
        printf "%i, " "$p" >> "$SCRIPT_DATA"
        printf "%f\n" "${EXEC_TIME}" >> "$SCRIPT_DATA"
    done
done

fi
