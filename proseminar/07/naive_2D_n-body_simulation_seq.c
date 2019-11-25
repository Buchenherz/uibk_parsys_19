#include <math.h>
#include <omp.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

typedef int value_t;

#define RESOLUTION_X 100
#define RESOLUTION_Y 20

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
    int max_Mass = 10;
    bool print_csv = false;
    clock_t clock_time;
    bool print = true;
    srand(111);

    if (argc == 4) {
        Nx = atoi(argv[1]);
        Ny = atoi(argv[2]);
        particle_count = atoi(argv[3]);
    } else if (argc == 3) {
        Nx = atoi(argv[1]);
        Ny = Nx;
        particle_count = atoi(argv[2]);
    } else if (argc < 4) {
        printf("Usage: ./2D_n-body_simulation_seq Nx Ny particles\n");
        return EXIT_FAILURE;
    }

    int T = 50;  // (Nx < Ny ? Ny : Nx) * 500;

    if (print) {
        printf(
            "Starting n-body simulation for Nx = %lld, Ny = %lld and "
            "particle_count = "
            "%d\n",
            Nx, Ny, particle_count);
    }

    clock_time = clock();
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

    if (print) {
        printf("Initial:\n");
        printParticles(particle_count, P, min_x, min_y, &Nx, &Ny);
        printf("\n");
    }

    // ---------- compute ----------

    // for each time step ..
    for (int t = 0; t < T; t++) {
        // Set min / max to 0 to minimize the print size - if the outest
        // particles move closer to the center
        *min_x = 0, *max_x = 0;
        *min_y = 0, *max_y = 0;

        // .. we propagate the positions
        double start = omp_get_wtime();
        omp_set_num_threads(8);
        for (int i = 0; i < particle_count; i++) {
            // printf("%d\n ", omp_get_thread_num());
            calculateParticleForces(P, i, particle_count);
        }

#pragma omp parallel for schedule(dynamic)
        for (int i = 0; i < particle_count; i++) {
            updateParticlePositions(P, i, min_x, max_x, min_y, max_y);
        }
        double end = omp_get_wtime();
        printf("Total time for iteration %d: %.16g\n", t, end - start);

        // show intermediate step
        if (!(t % 5) && print) {
            printf("Step t=%d:\n", t);
            printParticles(particle_count, P, min_x, min_y, max_x, max_y);
            printf("\n");
        }
    }

    // ---------- check ----------
    if (print) {
        printf("Final:\n");
        printParticles(particle_count, P, min_x, min_y, max_x, max_y);
        printf("minx: %lld, miny: %lld\n", *min_x, *min_y);
        printf("maxx: %lld, maxy: %lld\n", *max_x, *max_y);
        printf("\n");
    }
    printf("\n--FINAL POSITIONS OF PARTICLES--\n");
    for (size_t i = 0; i < particle_count; i++) {
        printf("%d %d", P[i].pos.x, P[i].pos.y);
    }
    printf("\n----\n");

    int success = true;
    // TODO validation ???
    // for (long long i = 0; i < particle_count; i++)
    // {
    //   if (?)
    //     continue;
    //   success = 0;
    //   break;
    // }

    // printf("Verification: %s\n", (success) ? "OK" : "FAILED");

    // ---------- cleanup ----------

    clock_time = clock() - clock_time;
    double time_taken = ((double)clock_time) / CLOCKS_PER_SEC;  // in seconds
    if (print) {
        printf("2D_n-body_simulation_seq took %f seconds to execute \n",
               time_taken);
    } else if (print_csv) {
        // printf( "Nx, Ny, particle_count,T, walltime, max_mass\n", Nx, Ny,
        // particle_count, T, time_taken, max_Mass);
        printf("%lld, %lld, %d, %d, %f, %d\n", Nx, Ny, particle_count, T,
               time_taken, max_Mass);
    }

    // free the min and max values
    free(min_x);
    free(max_x);
    free(min_y);
    free(max_y);
    // Free the Particle array
    free(P);

    // done
    return (success) ? EXIT_SUCCESS : EXIT_FAILURE;
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
        P[i].mass =
            1;  // max_Mass * rand_gen() + 1;  // float in range 1 to max_Mass;
        P[i].vel.x = 0;
        P[i].vel.y = 0;
    }
}

void calculateParticleForces(particle *P, int i, int particle_count) {
    for (int j = 0; j < particle_count; j++) {
        if (i == j) {
            continue;
        }
        float deltaX = P[j].pos.x - P[i].pos.x;
        float deltaY = P[j].pos.y - P[i].pos.y;

        // c² = a² + b² , no abs needed - square makes it positive
        double square_radius = pow(deltaX, 2) + pow(deltaY, 2);
        double force = (P[i].mass * P[j].mass) / square_radius;

        // https://stackoverflow.com/questions/39818833/moving-an-object-from-one-point-to-another
        // calculate angle
        float angle = atan2(deltaY, deltaX);

        P[i].vel.x += (force / P[i].mass) * cos(angle);
        P[i].vel.y += (force / P[i].mass) * sin(angle);

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
