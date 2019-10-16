#include<mpi.h>
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>

int main(int argc, char **argv) {
    // region OPENMPI INIT
    MPI_Init(&argc, &argv);  //initialize the MPI environment
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size); //get the number of ranks
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); //get the rankof the caller
    printf("Rank %d of %d online\n", rank, size);
    // endregion OPENMPI INIT

    double random_value;
    srand(time(NULL));
    unsigned long int gen, local_hits = 0;
    unsigned long int MAX_GEN_ELEMENTS = 1000000000;
    int num_gen_elements_per_proc = MAX_GEN_ELEMENTS / size;

    for (gen = 0; gen < num_gen_elements_per_proc; gen++) {
        // https://itp.tugraz.at/MML/MonteCarlo/MCIntro.pdf
        double random_x = (double) rand() / RAND_MAX; //float in range 0 to 1
        double random_y = (double) rand() / RAND_MAX;
        double random_point = pow(random_x, 2) + pow(random_y, 2);
        local_hits = (random_point < 1) ? ++local_hits : local_hits;
    }

    // Print the local_hits on each process
    printf("Local local_hits for process %d - local_hits %d - gens %d\n", rank, local_hits, gen);

    // Reduce all of the local local_hits into the global local_hits
    // This is, afaik, the longest possible long int. Just to be sure we have enough values to store our sum.
    unsigned long long int global_hits;
    MPI_Reduce(&local_hits, &global_hits, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        long int all_generated = size * num_gen_elements_per_proc;
        printf("PI = %f\n (All hits: %d | Genterated Points: %d)", ((double) global_hits * 4.0) / (double) all_generated,
               global_hits, all_generated);
    }

    MPI_Finalize(); //cleanup

    return EXIT_SUCCESS;
}