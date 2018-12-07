#!/bin/bash
set -e

OMP_NUM_THREADS=4

# MPICXX_FLAGS="-O3 -fopenmp"
MPICXX_FLAGS=""
MPIRUN_FLAGS=""

MPICXX=mpicxx

SRC_PATH="./src"
OBJ_PATH="./obj"
BIN_PATH="./bin"

echo "Clean"

rm -rf ${OBJ_PATH}/*.o
rm -rf ${BIN_PATH}/*

echo "Compile"

for file in $(ls ${SRC_PATH})
do
    echo "Source \"${SRC_PATH}/${file}\""
    ${MPICXX} ${MPICXX_FLAGS} -c -o ${OBJ_PATH}/${file//.cpp/.o} ${SRC_PATH}/$file
done

${MPICXX} ${MPICXX_FLAGS} -o ${BIN_PATH}/mpiwaves ${OBJ_PATH}/*.o

echo "Run"

time mpirun -n $1 ${BIN_PATH}/mpiwaves $2