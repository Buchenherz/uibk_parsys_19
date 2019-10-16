EXECS=pi_mpi
MPICC=mpicc
FLAGS=lm

all: $(EXECS)

$(EXECS): pi_mpi.o
	module load openmpi/4.0.1
	$(MPICC) pi_mpi.c -o $(EXECS) $(FLAGS)

.PHONY: clean
clean:
	rm $(EXECS)