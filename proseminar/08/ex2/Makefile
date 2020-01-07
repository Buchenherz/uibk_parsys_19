CC=gcc
MPICC=mpicc
CC_FLAGS=-O3 -std=c99 -Wall -Wextra -pedantic -lm 
OMP_FLAG=-fopenmp

.PHONEY: all
all: seq omp
seq: seq-mm
omp: omp-mm

seq-mm: matrix_mul_seq.c
	@$(CC) $(CC_FLAGS) $< -o matrix_mul_seq

omp-mm: matrix_mul_omp.c
	@$(CC) $(CC_FLAGS) $(OMP_FLAG) $< -o matrix_mul_omp

.PHONEY: clean
clean:
	@rm matrix_mul_seq
	@rm matrix_mul_omp

.PHONEY: run
run: all
	@echo "Sequential:"
	@./matrix_mul_seq 100 10000 1000
	@echo "OMP:"
	@./matrix_mul_omp 100 10000 1000