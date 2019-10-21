MPICC=mpic++
FLAGS=-lm
SHELL=/bin/bash
OBJECTS=pi_mpi.o pi_seq.o 

all: mpi seq

seq: pi_seq.o
	gcc pi_seq.c -o pi_seq $(FLAGS)

.PHONY: clean
clean:
	rm $(OBJECTS)