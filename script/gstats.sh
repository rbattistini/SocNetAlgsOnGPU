#!/bin/sh

: <<'DOC'
Questo script misura lo speedup e il throughput degli algoritmi paralleli per
il calcolo della Betweenness Centrality su GPU.
DOC

if [ $# != 2 ]; then
echo "Usage:
      gstats [executable] [dataset_directory]"
else

EXE=$1
DATASET_DIR=$2
DEVICE=0
FILES=$(find "${DATASET_DIR}" -type f -name "*mtx")
NTECHNIQUES=3

# per ogni file .mtx nel dataset
#   per ogni tecnica disponibile

for FILE in ${FILES}; do
    INPUT_FILE=$FILE
    FILENAME="$(basename "$FILE" | cut -f 1 -d '.')"
    SCORES_FILE="scores_""${FILENAME}""$(basename "$FILE" | cut -f 1 -d '.')"".csv"
    STATS_FILE="stats_""${FILENAME}"".csv"

    printf '"Technique", "Total Time", "Load Time", "Unload Time", "BC Computation Time", "MTEPS"\n' >> "${STATS_FILE}"

    for TECH in $(seq "$((NTECHNIQUES - 1))"); do
        ${EXE} -i "${INPUT_FILE}" -t "${TECH}" -q -l -d ${DEVICE} -s "${STATS_FILE}"
    done
    ${EXE} -i "${INPUT_FILE}" -t ${NTECHNIQUES} -q -l -d ${DEVICE} -s "${STATS_FILE}" -b "${SCORES_FILE}"
done
fi
