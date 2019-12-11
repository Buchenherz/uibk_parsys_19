#!/bin/bash

if [$# == 0] 
then 
    echo "Usage: ./sr <job_name> <mode> <slots> <make>"

job_name=$1

# openmp / openmpi
mode=$2
slots=$3

# what should the makefile make
make=$4

