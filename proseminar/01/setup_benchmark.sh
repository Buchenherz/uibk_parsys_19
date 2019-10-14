#!/usr/bin/env bash
mkdir benchmark
cd benchmark
wget http://mvapich.cse.ohio-state.edu/download/mvapich/osu-micro-benchmarks-5.6.2.tar.gz
tar -zxvf osu-micro-benchmarks-5.6.2.tar.gz
module load openmpi/4.0.1
./configure CC=mpicc CXX=mpic++
make