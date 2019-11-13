CC=gcc
MPICC=mpic++
CC_FLAGS=-O2 -std=c99 -Wall -Wextra -pedantic -lm

.PHONEY: all
all: 2D_n-body_simulation_seq 

seq: 2D_n-body_simulation_seq.c
	@$(CC) $(CC_FLAGS) $< -o 2D_n-body_simulation_seq

.PHONEY: clean
clean:
	@rm 2D_n-body_simulation_seq

.PHONEY: run
run: all
	@echo "Sequential:"
	@./2D_n-body_simulation_seq
