#!/bin/bash
set -e

MPICXX_FLAGS="-O3 -qsmp=omp"
MPIRUN_FLAGS=""

MPICXX=mpixlcxx_r

SRC_PATH="./src"
OBJ_PATH="./obj"
BIN_PATH="./bin"
LOG_PATH="./logs"

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