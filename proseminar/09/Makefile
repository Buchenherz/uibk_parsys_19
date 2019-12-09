CC=gcc-9
MPICC=mpicc
CC_FLAGS=-O3 -std=c99 -Wall -Wextra -pedantic -lm 
OMP_FLAGS=-lm -fopenmp -O3 -std=c99 -Wall -Wextra -pedantic 

.PHONEY: all
all: seq omp
seq: seq-nqueens
omp: omp-nqueens

seq-nqueens: naive-nqueens.c
	@$(CC) $(OMP_FLAGS) $< -o naive-nqueens.out


omp-nqueens: naive-nqueens-omp.c
	@$(CC) $(OMP_FLAGS) $< -o naive-nqueens-omp.out


.PHONEY: clean
clean:
	@rm naive-nqueens
	@rm naive-nqueens-omp
	

.PHONEY: run
run: all
	@echo "Sequential:"
	@./naive-nqueens 4
	@echo "OMP:"
	@./naive-nqueens-omp 4
	