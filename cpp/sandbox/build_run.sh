#!/bin/bash
set -e

# MPIRUN_FLAGS="--tag-output"
MPIRUN_FLAGS=""

mpicxx -o $1 ./$1.cpp

mpirun $MPIRUN_FLAGS -n $2 ./$1