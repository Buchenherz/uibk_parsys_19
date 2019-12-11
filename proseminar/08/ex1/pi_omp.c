#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

// #define VERBOSE

int main(int argc, char *argv[]) {
    double start_time = omp_get_wtime();
    srand(time(NULL));

    unsigned long long int MAX_GEN_ELEMENTS;

    // https://www.tutorialspoint.com/cprogramming/c_command_line_arguments.htm

    if (argc == 2) {
        MAX_GEN_ELEMENTS = atoi(argv[1]);
    } else {
        printf("Usage: ./pi_mpi <Number of gen elements> <Number of chunks>\n");
        return EXIT_FAILURE;
    }

    unsigned long long int global_hits = 0;

#ifdef VERBOSE
    printf("num_elements_per_thread %d\n", num_elements_per_thread);
#endif

// Split the parallel workload up
#pragma omp parallel for reduction(+ : global_hits)
    for (unsigned long long int i = 0; i < MAX_GEN_ELEMENTS; i++) {
        // https://itp.tugraz.at/MML/MonteCarlo/MCIntro.pdf
        // floats in range 0 to 1 

        // This version of rand is blocking when used by multiple threads, which is why there is a slowdown.
        // There is also a version of rand, called `rand_r()` that needs to be used. One could use the thread number
        // as the seed pointer.
        
        double random_x = (double)rand() / RAND_MAX;
        double random_y = (double)rand() / RAND_MAX;
        double random_point = pow(random_x, 2) + pow(random_y, 2);
        if (random_point < 1) {
            global_hits++;
        }
    }

#ifdef VERBOSE
    printf("Local hits: %d\n", local_hits);
#endif
    double pi = (4.0 * (double)global_hits) / (double)MAX_GEN_ELEMENTS;
    double end_time = omp_get_wtime();
#ifdef VERBOSE
    printf("Total hits: %d ", global_hits);
#endif
#ifdef PRINT_CSV_HEADER
    printf("PI, Walltime, Number_of_elements\n");
#endif
    printf("%f, %f, %llu\n", pi, end_time - start_time, MAX_GEN_ELEMENTS);
}