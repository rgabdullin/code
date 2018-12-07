#!/bin/bash
set -e

# MPIRUN_FLAGS="--tag-output"
MPIRUN_FLAGS=""

mpicxx -o mpi1 ./mpi1.cpp

mpirun $MPIRUN_FLAGS -n $1 ./mpi1