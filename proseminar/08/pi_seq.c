#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int main(int argc, char *argv[]) {
    double start_time = omp_get_wtime();
    srand(time(NULL));

    long int MAX_GEN_ELEMENTS;

    // region user input
    // https://www.tutorialspoint.com/cprogramming/c_command_line_arguments.htm
    if (argc == 2) {
        MAX_GEN_ELEMENTS = atoi(argv[1]);
    } else {
        printf("Usage: ./pi_mpi <Number of gen elements>\n");
        return EXIT_FAILURE;
    }
    // endregion
    int hits = 0;
    int gen = 0;
    for (gen; gen < MAX_GEN_ELEMENTS; gen++) {
        // https://itp.tugraz.at/MML/MonteCarlo/MCIntro.pdf
        double random_x = (double)rand() / RAND_MAX;  // float in range 0 to 1
        double random_y = (double)rand() / RAND_MAX;
        double random_point = pow(random_x, 2) + pow(random_y, 2);
        hits = (random_point < 1) ? ++hits : hits;
    }
    double pi = ((double)hits / (double)gen) * 4;
    double end_time = omp_get_wtime();
    printf("Gen =  %d\n", gen);
    printf("PI =  %f\n", pi);
    printf("Walltime =  %f\n", end_time - start_time);
}