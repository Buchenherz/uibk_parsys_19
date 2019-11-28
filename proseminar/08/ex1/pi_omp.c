#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

// #define VERBOSE

int main(int argc, char *argv[]) {
    double start_time = omp_get_wtime();
    srand(time(NULL));

    long int MAX_GEN_ELEMENTS;
    int number_of_chunks = 8;

    // https://www.tutorialspoint.com/cprogramming/c_command_line_arguments.htm
    if (argc == 3) {
        MAX_GEN_ELEMENTS = atoi(argv[1]);
        number_of_chunks = atoi(argv[2]);
    } else if (argc == 2) {
        MAX_GEN_ELEMENTS = atoi(argv[1]);
    } else {
        printf("Usage: ./pi_mpi <Number of gen elements> <Number of chunks>\n");
        return EXIT_FAILURE;
    }

    int global_hits = 0;

    // calculate amount of work to be done by each thread
    int num_elements_per_thread = MAX_GEN_ELEMENTS / number_of_chunks;
#ifdef VERBOSE
    printf("num_elements_per_thread %d\n", num_elements_per_thread);
#endif

// Split the parallel workload up
#pragma omp parallel for reduction(+ : global_hits)
    for (int chunks = 0; chunks < number_of_chunks; chunks++) {
        // Set up local hits of chunk
        int local_hits = 0;
        for (int chunksize = 0; chunksize < num_elements_per_thread;
             chunksize++) {
            // https://itp.tugraz.at/MML/MonteCarlo/MCIntro.pdf
            // floats in range 0 to 1
            double random_x = (double)rand() / RAND_MAX;
            double random_y = (double)rand() / RAND_MAX;
            double random_point = pow(random_x, 2) + pow(random_y, 2);
            if (random_point < 1) {
                local_hits++;
            }
        }
        // Reduce local hits to global hits
        global_hits += local_hits;
#ifdef VERBOSE
        printf("Local hits: %d\n", local_hits);
#endif
    }
    double pi = (4.0 * (double)global_hits) / (double)MAX_GEN_ELEMENTS;
    double end_time = omp_get_wtime();
#ifdef VERBOSE
    printf("Total hits: %d ", global_hits);
#endif
    printf("PI, Walltime\n");
    printf("%f, %f\n", pi, (end_time - start_time));
}