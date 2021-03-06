#include <math.h>
#include <mpi.h>
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
    int max_Mass = 10000;
    bool print_csv = false;
    clock_t clock_time;
    bool print = true;
    srand(time(NULL));

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

    int T = 5000;  // (Nx < Ny ? Ny : Nx) * 500;

    // #region MPI INIT
    MPI_Init(&argc, &argv);  // initialize the MPI environment
    int number_of_ranks;
    int rank;
    MPI_Comm_size(MPI_COMM_WORLD, &number_of_ranks);  // get the number of ranks
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);  // get the rankof the caller

    // #endregion

    // #region Set up MPI Type Structs
    MPI_Datatype mpi_pos;
    int pos_blocklengths[2] = {1, 1};
    MPI_Aint pos_displacements[2] = {offsetof(struct position, x),
                                     offsetof(struct position, y)};
    MPI_Datatype pos_datatypes[2] = {MPI_DOUBLE, MPI_DOUBLE};
    MPI_Type_create_struct(2, pos_blocklengths, pos_displacements,
                           pos_datatypes, &mpi_pos);
    MPI_Type_commit(&mpi_pos);

    MPI_Datatype mpi_vel;
    int vel_blocklengths[2] = {1, 1};
    MPI_Aint vel_displacements[2] = {offsetof(struct velocity, x),
                                     offsetof(struct velocity, y)};
    MPI_Datatype vel_datatypes[2] = {MPI_DOUBLE, MPI_DOUBLE};
    MPI_Type_create_struct(2, vel_blocklengths, vel_displacements,
                           vel_datatypes, &mpi_vel);
    MPI_Type_commit(&mpi_vel);

    MPI_Datatype mpi_particle;
    int part_blocklengths[3] = {1, 1, 1};
    MPI_Aint part_displacements[3] = {offsetof(particle, pos),
                                      offsetof(particle, vel),
                                      offsetof(particle, mass)};
    MPI_Datatype part_datatypes[3] = {mpi_pos, mpi_vel, MPI_DOUBLE};
    MPI_Type_create_struct(3, part_blocklengths, part_displacements,
                           part_datatypes, &mpi_particle);
    MPI_Type_commit(&mpi_particle);
    // #endregion
    clock_time = clock();
    // ---------- setup ----------

    particle *P = malloc(particle_count * sizeof(*P));

    if (print && rank == 0) {
        printf(
            "Starting n-body simulation for Nx = %lld, Ny = %lld and "
            "particle_count = "
            "%d\n",
            Nx, Ny, particle_count);
    }

    if (rank == 0) {
        // set up initial conditions in P
        initParticles(P, particle_count, Nx, Ny, max_Mass);
    }

    MPI_Bcast(&P[0], particle_count, mpi_particle, 0, MPI_COMM_WORLD);

    // allocate min / max values for future print scale
    long long *min_x = malloc(sizeof(long long)),
              *max_x = malloc(sizeof(long long)),
              *min_y = malloc(sizeof(long long)),
              *max_y = malloc(sizeof(long long));

    // Init min / max
    *min_x = 0, *max_x = 0;
    *min_y = 0, *max_y = 0;

    if (print && rank == 0) {
        printf("Initial:\n");
        printParticles(particle_count, P, min_x, min_y, &Nx, &Ny);
        printf("\n");
    }
    // ---------- compute ----------

    int local_particle_count = particle_count / number_of_ranks;
    // for each time step ..
    for (int t = 0; t < T; t++) {
        // printf("lpc: %d, rank: %d\n", local_particle_count, rank);
        // Set min / max to 0 to minimize the print size - if the outest
        // particles move closer to the center
        *min_x = 0, *max_x = 0;
        *min_y = 0, *max_y = 0;
        // .. we propagate the positions
        particle sbuf[local_particle_count];
        int j = 0;
        for (int i = local_particle_count * rank;
             i < local_particle_count * (rank + 1); i++) {
            calculateParticleForces(P, i, local_particle_count);
            updateParticlePositions(P, i, min_x, max_x, min_y, max_y);
            sbuf[j] = P[i];

            j++;
        }

        particle rbuf[local_particle_count];
        MPI_Allgather(&sbuf, local_particle_count, mpi_particle, &P[0],
                      local_particle_count, mpi_particle, MPI_COMM_WORLD);

        // show intermediate step
        if (!(t % 100) && print) {
            // printf("min_x = %lld, max_x = %lld, min_y = %lld, max_y =
            // %lld",
            //        *min_x, *max_x, *min_y, *max_y);
            if (rank == 0) {
                // https://stackoverflow.com/questions/17741574/in-place-mpi-reduce-crashes-with-openmpi
                MPI_Reduce(MPI_IN_PLACE, min_x, 1, MPI_LONG_LONG, MPI_MIN, 0,
                           MPI_COMM_WORLD);
                MPI_Reduce(MPI_IN_PLACE, max_x, 1, MPI_LONG_LONG, MPI_MAX, 0,
                           MPI_COMM_WORLD);
                MPI_Reduce(MPI_IN_PLACE, min_y, 1, MPI_LONG_LONG, MPI_MIN, 0,
                           MPI_COMM_WORLD);
                MPI_Reduce(MPI_IN_PLACE, max_y, 1, MPI_LONG_LONG, MPI_MAX, 0,
                           MPI_COMM_WORLD);
            } else {
                MPI_Reduce(min_x, min_x, 1, MPI_LONG_LONG, MPI_MIN, 0,
                           MPI_COMM_WORLD);
                MPI_Reduce(max_x, max_x, 1, MPI_LONG_LONG, MPI_MAX, 0,
                           MPI_COMM_WORLD);
                MPI_Reduce(min_y, min_y, 1, MPI_LONG_LONG, MPI_MIN, 0,
                           MPI_COMM_WORLD);
                MPI_Reduce(max_y, max_y, 1, MPI_LONG_LONG, MPI_MAX, 0,
                           MPI_COMM_WORLD);
            }
            if (rank == 0) {
                printf("Step t=%d:\n", t);
                printParticles(particle_count, P, min_x, min_y, max_x, max_y);
                printf("\n");
            }
        }
    }

    // ---------- check ----------
    if (print && rank == 0) {
        printf("Final:\n");
        printParticles(particle_count, P, min_x, min_y, max_x, max_y);
        printf("minx: %lld, miny: %lld\n", *min_x, *min_y);
        printf("maxx: %lld, maxy: %lld\n", *max_x, *max_y);
        printf(
            "END - No.Ranks: %d, Timesteps: %d, particle_count: %d, x: %lld, "
            "y: "
            "%lld\n",
            number_of_ranks, T, particle_count, Nx, Ny);
    }

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
    if (print && rank == 0) {
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

    // Finalise MPI
    MPI_Type_free(&mpi_particle);
    MPI_Type_free(&mpi_pos);
    MPI_Type_free(&mpi_vel);
    MPI_Finalize();

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
