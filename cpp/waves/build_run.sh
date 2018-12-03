#!/bin/bash
set -e

export OMP_NUM_THREADS=4

FLAGS="-O3 -fopenmp"

CXX=g++
CC=gcc

LIB_PATH="/home/ruslixag/libs"
SRC_PATH="./src"
OBJ_PATH="./obj"
BIN_PATH="./bin"

echo "Compile"

for file in $(ls ${SRC_PATH})
do
    echo "Source \"${SRC_PATH}/${file}\""
    ${CXX} ${FLAGS} --std=c++11 -I "${LIB_PATH}/tqdm.cpp/include" -c -o ${OBJ_PATH}/*.o ${SRC_PATH}/*.cpp
done

${CXX} ${FLAGS} -o ${BIN_PATH}/waves ${OBJ_PATH}/*.o

echo "Run"

time ${BIN_PATH}/waves