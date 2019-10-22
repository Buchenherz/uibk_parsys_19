# Exercise 1

## Tasks

### Write a sequential application `pi_seq` in C or C++ that computes π for a given number of samples (command line argument). Test your application for various, large sample sizes to verify the correctness of your implementation.

You can find the sequential implementation of our program in the file `pi_seq.c`. Within the comments are links to resources we used for information gathering. You can use `make seq` to compile it.

### Consider a parallelization strategy using MPI. Which communication pattern(s) would you choose and why?

We split the number of generated elements evenly across all ranks and then reduce the sum of these local hits to rank 0, where the average is build from all gathered hits multiplied by 4 (as specified in the definition of monte carlo)

### Implement your chosen parallelization strategy as a second application `pi_mpi`. Run it with varying numbers of ranks and sample sizes and verify its correctness by comparing the output to `pi_seq`.

The parallel program can be found in the file `pi_mpi.c`. We had problems with the Makefile which is why the program is compiled when submitting the job to the cluster. The job file is `pi_mpi.script`. Excel, Tableau and Graphics files are within this folder as well.

### Discuss the effects and implications of your parallelization.

Each rank but rank 0 has the same workload, e.g. number of times hits are calculated. Only after calculation, a single unsigned long int per rank is sent to rank 0, where
the rest of the calculation takes place. One should make sure that the input is divisible by the amount of ranks used to avoid splitting errors.
