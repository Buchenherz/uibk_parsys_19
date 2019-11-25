CC=gcc
MPICC=mpicc
CC_FLAGS=-O2 -std=c99 -Wall -Wextra -pedantic -lm -fopenmp

.PHONEY: all
all: seq

seq: naive_2D_n-body_simulation_seq.c
	@$(CC) $(CC_FLAGS) $< -o 2D_n-body_simulation_seq_omp

.PHONEY: clean
clean:
	@rm 2D_n-body_simulation_seq
	

.PHONEY: run
run: all
	@echo "Sequential OMP:"
	@./2D_n-body_simulation_seq
	