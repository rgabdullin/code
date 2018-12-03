#!/bin/bash
set -e

CXX=g++
NVCC=nvcc

SRC_PATH="./src"
OBJ_PATH="./obj"
BIN_PATH="./bin"

echo "Compile"

${NVCC} --std=c++11 -c -o ${OBJ_PATH}/$1.o ${SRC_PATH}/$1.cpp
${CXX} -o ${BIN_PATH}/$1 ${OBJ_PATH}/$1.o

echo "Run"

${BIN_PATH}/$1