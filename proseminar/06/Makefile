CC=gcc
MPICC=mpicc
CC_FLAGS=-O2 -std=c99 -Wall -Wextra -pedantic -lm

.PHONEY: all
all: seq mpi

seq: 2D_n-body_simulation_seq.c
	@$(CC) $(CC_FLAGS) $< -o 2D_n-body_simulation_seq

mpi: 2D_n-body_simulation_mpi.c
	@$(MPICC) 2D_n-body_simulation_mpi.c -o 2D_n-body_simulation_mpi -lm

.PHONEY: clean
clean:
	@rm 2D_n-body_simulation_seq
	@rm 2D_n-body_simulation_mpi

.PHONEY: run
run: all
	@echo "Sequential:"
	@./2D_n-body_simulation_seq
	@echo "MPI:"
	@./2D_n-body_simulation_mpi