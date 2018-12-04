#!/bin/bash
set -e

export OMP_NUM_THREADS=4

FLAGS="-O3 -fopenmp"

CXX=g++
CC=gcc

INCLUDE_PATH="./include"
SRC_PATH="./src"
OBJ_PATH="./obj"
BIN_PATH="./bin"

echo "Compile"

for file in $(ls ${SRC_PATH})
do
    echo "Source \"${SRC_PATH}/${file}\""
    ${CXX} ${FLAGS} -I ${INCLUDE_PATH} -c -o ${OBJ_PATH}/${file//.cpp/.o} ${SRC_PATH}/$file
done

${CXX} ${FLAGS} -o ${BIN_PATH}/waves ${OBJ_PATH}/*.o

echo "Run"

time ${BIN_PATH}/waves