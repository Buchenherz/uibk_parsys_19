#include <math.h>
#include <omp.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

typedef int value_t;

#define RESOLUTION_X 100
#define RESOLUTION_Y 20

// IFDEFS
// #define DEBUG
// #define PRINT
#define CSV

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// -- Matrix utilities --

typedef struct particle {
    struct position {
        double x;
        double y;
    } pos;
    struct velocity {
        double x;
        double y;
    } vel;
    double mass;
} particle;

typedef value_t **Matrix;

Matrix createMatrix(int Nx, int Ny);
void releaseMatrix(Matrix m);

void initParticles(particle *P, int particle_count, long long Nx, long long Ny,
                   int max_Mass);

/**
 * Calculates the forces of particle `i` in an array of particles `P` where P
 * contains `particle_count` total particles.
 */
void calculateParticleForces(particle *P, int i, int particle_count);

/**
 * Updates the position of particle `i` in an array of particles `P` containing
 * `particle_count` particles.
 */
void updateParticlePositions(particle *P, int i, long long *min_x,
                             long long *max_x, long long *min_y,
                             long long *max_y);

void printParticles(int particle_count, struct particle *P, long long *Mx,
                    long long *My, long long *Nx, long long *Ny);

long long MIN(long long x, long long y);
long long MAX(long long x, long long y);
double rand_gen();
// use it: sigma * normalRandom + µ
double normalRandom();
// -- simulation code ---

int main(int argc, char **argv) {
    // 'parsing' optional input parameter = problem size

    // Initial set of the scale of the universe
    long long Nx = 2000;        // columns
    long long Ny = 2000;        // rows
    int particle_count = 2000;  // particles

    int max_Mass = 1;

#ifdef DEBUG
    // Fixed random generator for debug purposes
    srand(111);
#else
    srand(time(NULL));
#endif

    int number_of_threads = 4;
    if (argc == 5) {
        Nx = atoi(argv[1]);
        Ny = atoi(argv[2]);
        particle_count = atoi(argv[3]);
        number_of_threads = atoi(argv[4]);
    } else if (argc == 4) {
        Nx = atoi(argv[1]);
        Ny = atoi(argv[2]);
        particle_count = atoi(argv[3]);
    } else if (argc == 3) {
        Nx = atoi(argv[1]);
        Ny = Nx;
        particle_count = atoi(argv[2]);
    } else if (argc < 4) {
        printf(
            "Usage: ./2D_n-body_simulation_seq Nx Ny particles num_threads\n");
        return EXIT_FAILURE;
    }

    int T = 50;  // (Nx < Ny ? Ny : Nx) * 500;
    double iteration_times_array[T];
#ifdef DEBUG
    printf(
        "Starting n-body simulation for Nx = %lld, Ny = %lld, "
        "no_of_threads = "
        "%d, timesteps = %d and "
        "particle_count = "
        "%d\n",
        Nx, Ny, number_of_threads, T, particle_count);
#endif

    double omp_start_time = omp_get_wtime();

    // ---------- setup ----------

    particle *P = malloc(particle_count * sizeof(*P));

    // set up initial conditions in P
    initParticles(P, particle_count, Nx, Ny, max_Mass);

    // allocate min / max values for future print scale
    long long *min_x = malloc(sizeof(long long)),
              *max_x = malloc(sizeof(long long));
    long long *min_y = malloc(sizeof(long long)),
              *max_y = malloc(sizeof(long long));

    // Init min / max
    *min_x = 0, *max_x = 0;
    *min_y = 0, *max_y = 0;

#ifdef PRINT
    printf("Initial:\n");
    printParticles(particle_count, P, min_x, min_y, &Nx, &Ny);
    printf("\n");
#endif

    // ---------- compute ----------

    // for each time step ..
    for (int t = 0; t < T; t++) {
        // Set min / max to 0 to minimize the print size - if the outest
        // particles move closer to the center
        *min_x = 0, *max_x = 0;
        *min_y = 0, *max_y = 0;

        // .. we propagate the positions
        double iteration_start_time = omp_get_wtime();

        omp_set_num_threads(number_of_threads);
#pragma omp parallel for schedule(static, 1)
        // https://software.intel.com/en-us/articles/openmp-loop-scheduling
        // http://ppc.cs.aalto.fi/ch3/schedule/
        for (int i = 0; i < particle_count; i++) {
            calculateParticleForces(P, i, particle_count);
        }

        for (int i = 0; i < particle_count; i++) {
            updateParticlePositions(P, i, min_x, max_x, min_y, max_y);
        }
        double iteration_end_time = omp_get_wtime();
        iteration_times_array[t] = iteration_end_time - iteration_start_time;

#ifdef DEBUG
        printf("Total time for iteration %d: %.16g\n", t, end - start);
#endif
        // show intermediate step
        if (!(t % 5)) {
#ifdef PRINT
            printf("Step t=%d:\n", t);
            printParticles(particle_count, P, min_x, min_y, max_x, max_y);
            printf("\n");
#endif
        }
    }

    // ---------- check ----------
#ifdef PRINT
    printf("Final:\n");
    printParticles(particle_count, P, min_x, min_y, max_x, max_y);
    printf("minx: %lld, miny: %lld\n", *min_x, *min_y);
    printf("maxx: %lld, maxy: %lld\n", *max_x, *max_y);
    printf("\n");
#endif
#ifdef DEBUG
    printf("\n--FINAL POSITIONS OF PARTICLES--\n");
    for (int i = 0; i < particle_count; i++) {
        printf("%f %f ", P[i].pos.x, P[i].pos.y);
    }
    printf("\n----\n");
#endif

    double omp_end_time = omp_get_wtime();
    double iteration_avg = 0.0;
    for (int t = 0; t < T; t++) {
        iteration_avg += iteration_times_array[t];
    }
    iteration_avg /= (double)T;

#ifdef PRINT
    printf("Total time for program execution: %.16g\n",
           omp_end_time - omp_start_time);
#endif
#ifdef CSV
    printf("%lld, %lld, %d, %d, %f, %f, %d, %d\n", Nx, Ny, particle_count, T,
           omp_end_time - omp_start_time, iteration_avg, max_Mass,
           number_of_threads);
#endif

    // ---------- cleanup ----------
    // free the min and max values
    free(min_x);
    free(max_x);
    free(min_y);
    free(max_y);
    // Free the Particle array
    free(P);

    // done
    return EXIT_SUCCESS;
}

Matrix createMatrix(int Nx, int Ny) {
    // create data and index Matrix
    value_t **mat = calloc(Ny, sizeof(int *));
    for (int i = 0; i < Ny; i++) {
        mat[i] = calloc(Nx, sizeof(int));
    }
    return mat;
}

void releaseMatrix(Matrix m) { free(m); }

// TODO refactor - extract to methode
// particleComputation()

void printParticles(int particle_count, struct particle *P, long long *Mx,
                    long long *My, long long *Nx, long long *Ny) {
    const char *colors = " .-:=+*^X#%@";
    const int numColors = 12;

    // set the 'render' resolution
    int W = RESOLUTION_X;
    int H = RESOLUTION_Y;

    // step size in each dimension
    int xW = (*Nx - *Mx) / W;
    int yW = (*Ny - *My) / H;

    Matrix A = createMatrix(W, H);

    for (int i = 0; i < particle_count; i++) {
        int y = MAX(MIN((P[i].pos.y - *My) / yW, H - 1), 0);
        int x = MAX(MIN((P[i].pos.x - *Mx) / xW, W - 1), 0);
        A[y][x]++;
    }

    // room
    for (int i = 0; i < H; i++) {
        if (i == 0) {
            printf("┌");
            for (int k = 0; k < W; k++) {
                printf("─");
            }
            printf("┐\n");
        }
        // left wall
        printf("│");
        // actual room
        for (int j = 0; j < W; j++) {
            int c = A[i][j] > numColors ? 11 : A[i][j];
            printf("%c", colors[c]);
        }
        // right wall
        printf("│\n");
        if (i == H - 1) {
            printf("└");
            for (int k = 0; k < W; k++) {
                printf("─");
            }
            printf("┘\n");
        }
    }

    releaseMatrix(A);
}

long long MIN(long long x, long long y) { return x > y ? y : x; }

long long MAX(long long x, long long y) { return x < y ? y : x; }

// return a uniformly distributed random value
double rand_gen() { return (double)rand() / (double)(RAND_MAX); }

// return a normally distributed random value
// (https://en.wikipedia.org/wiki/Box%E2%80%93Muller_transform)
double normalRandom() {
    double v1 = rand_gen();
    double v2 = rand_gen();
    return cos(2. * M_PI * v2) * sqrt(-2. * log(v1));
}

void initParticles(particle *P, int particle_count, long long Nx, long long Ny,
                   int max_Mass) {
    // #pragma omp parallel for schedule(dynamic) num_threads(1)
    // Even with a fixed seed changes may occur here because the access to a
    // rand variable is different on each run
    for (int i = 0; i < particle_count; i++) {
        P[i].pos.x =
            Nx / 10 * normalRandom() + Nx / 2;  // float in range 0 to Nx;
        P[i].pos.y =
            Ny / 10 * normalRandom() + Ny / 2;  // float in range 0 to Ny;
        P[i].mass = max_Mass * rand_gen() + 1;  // float in range 1 to max_Mass;
        P[i].vel.x = 0;
        P[i].vel.y = 0;
    }
}

void calculateParticleForces(particle *P, int i, int particle_count) {
    for (int j = 0; j < particle_count; j++) {
        if (i == j) {
            continue;
        }

        double deltaX, deltaY, force, old_vel_x, old_vel_y, j_pos_x, j_pos_y,
            i_pos_x, i_pos_y;

        j_pos_x = P[j].pos.x;

        j_pos_y = P[j].pos.y;

        i_pos_x = P[i].pos.x;

        i_pos_y = P[i].pos.y;

        deltaX = j_pos_x - i_pos_x;
        deltaY = j_pos_y - i_pos_y;

        // c² = a² + b² , no abs needed - square makes it positive
        double square_radius = pow(deltaX, 2) + pow(deltaY, 2);
        force = (P[i].mass * P[j].mass) / square_radius;

        // https://stackoverflow.com/questions/39818833/moving-an-object-from-one-point-to-another
        // calculate angle
        float angle = atan2(deltaY, deltaX);

        old_vel_x = P[i].vel.x;
        old_vel_y = P[i].vel.y;

        double new_vel_x = old_vel_x + (force / P[i].mass) * cos(angle);
        double new_vel_y = old_vel_y + (force / P[i].mass) * sin(angle);

        // When using atomics, this gets waaaay faster (from .5 per iteration
        // down to .1)
#pragma omp atomic write
        P[i].vel.x = new_vel_x;
#pragma omp atomic write
        P[i].vel.y = new_vel_y;

        // P[j].vel.x += (force / P[j].mass) * cos(angle + M_PI);
        // P[j].vel.y += (force / P[j].mass) * sin(angle + M_PI);
    }
}

void updateParticlePositions(particle *P, int i, long long *min_x,
                             long long *max_x, long long *min_y,
                             long long *max_y) {
    P[i].pos.x += P[i].vel.x;
    P[i].pos.y += P[i].vel.y;

    // Check and set new minimum and maximum values for future
    // print scaling
    *min_x = MIN((long long)P[i].pos.x, *min_x);
    *max_x = MAX((long long)P[i].pos.x, *max_x);
    *min_y = MIN((long long)P[i].pos.y, *min_y);
    *max_y = MAX((long long)P[i].pos.y, *max_y);
}
