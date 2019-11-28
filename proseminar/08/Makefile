CC=gcc-9
MPICC=mpicc
CC_FLAGS=-O3 -std=c99 -Wall -Wextra -pedantic -lm 
OMP_FLAGS=-lm -fopenmp -O3 -std=c99 -Wall -Wextra -pedantic 

.PHONEY: all
all: seq omp
seq: seq-pi seq-hs
omp: omp-pi omp-hs

seq-pi: pi_seq.c
	@$(CC) $(OMP_FLAGS) $< -o pi_seq

seq-hs: heat_stencil_2D_seq.c
	@$(CC) $(OMP_FLAGS) $< -o heat_stencil_2D_seq

omp-pi: pi_omp.c
	@$(CC) $(OMP_FLAGS) $< -o pi_omp

omp-hs: heat_stencil_2D_omp.c
	@$(CC) $(OMP_FLAGS) $< -o heat_stencil_2D_omp

.PHONEY: clean
clean:
	@rm heat_stencil_2D_seq
	@rm pi_seq
	@rm heat_stencil_2D_omp
	@rm pi_omp
	

.PHONEY: run
run: all
	@echo "Sequential OMP:"
	@./2D_n-body_simulation_seq
	