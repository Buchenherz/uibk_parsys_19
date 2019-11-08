#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

typedef int value_t;

#define RESOLUTION_X 100
#define RESOLUTION_Y 20

// -- Matrix utilities --

struct particle {
    struct pos {
        double x;
        double y;
    } pos;
    struct velocity {
        double x;
        double y;
    } velocity;
    double mass;
};

typedef value_t **Matrix;

Matrix createMatrix(int Nx, int Ny);

void releaseMatrix(Matrix m);

void printParticles(int particle_count, struct particle *P, int Mx, int My,
                    int Nx, int Ny);

long long MIN(long long x, long long y);
long long MAX(long long x, long long y);
double rand_gen();
// use it: sigma * normalRandom + µ
double normalRandom();
// -- simulation code ---

int main(int argc, char **argv) {
    // 'parsing' optional input parameter = problem size
    int Nx = 2000;              // columns
    int Ny = 2000;              // rows
    int particle_count = 1000;  // particles
    int max_Mass = 10000000;
    clock_t clock_time;
    bool print = true;
    srand(time(NULL));

    if (argc > 1) {
        Nx = atoi(argv[1]);
        Ny = atoi(argv[2]);
        particle_count = atoi(argv[3]);
    }
    int T = 50;  // (Nx < Ny ? Ny : Nx) * 500;

    clock_time = clock();
    // ---------- setup ----------

    struct particle P[particle_count];

    // set up initial conditions in P
    for (int i = 0; i < particle_count; i++) {
        P[i].pos.x =
            Nx / 10 * normalRandom() + Nx / 2;  // float in range 0 to Nx;
        P[i].pos.y =
            Ny / 10 * normalRandom() + Ny / 2;  // float in range 0 to Ny;
        P[i].mass = max_Mass * rand_gen() + 1;  // float in range 1 to max_Mass;
        P[i].velocity.x = 0;
        P[i].velocity.y = 0;
    }

    if (print) {
        printf("Initial:\n");
        printParticles(particle_count, P, 0, 0, Nx, Ny);
        printf("\n");
    }

    // ---------- compute ----------

    long long min_x = 0, max_x = 0;
    long long min_y = 0, max_y = 0;

    // for each time step ..
    for (int t = 0; t < T; t++) {
        min_x = 0, max_x = 0;
        min_y = 0, max_y = 0;
        // .. we propagate the positions
        for (int i = 0; i < particle_count; i++) {
            for (int j = i + 1; j < particle_count; j++) {
                float deltaX = P[j].pos.x - P[i].pos.x;
                float deltaY = P[j].pos.y - P[i].pos.y;

                // c² = a² + b² , no abs needed - square makes it positive
                double square_radius = pow(deltaX, 2) + pow(deltaY, 2);
                double force = (P[i].mass * P[j].mass) / square_radius;

                // https://stackoverflow.com/questions/39818833/moving-an-object-from-one-point-to-another
                // calculate angle
                float angle = atan2(deltaY, deltaX);
                P[i].velocity.x += (force / P[i].mass) * cos(angle);
                P[i].velocity.y += (force / P[i].mass) * sin(angle);
                P[j].velocity.x += (force / P[j].mass) * cos(angle + M_PI);
                P[j].velocity.y += (force / P[j].mass) * sin(angle + M_PI);
            }
            P[i].pos.x += P[i].velocity.x;
            P[i].pos.y += P[i].velocity.y;

            min_x = MIN(P[i].pos.x, min_x);
            max_x = MAX(P[i].pos.x, max_x);
            min_y = MIN(P[i].pos.y, min_y);
            max_y = MAX(P[i].pos.y, max_y);
        }
        Nx = max_x;
        Ny = max_y;

        // show intermediate step
        if (!(t % 5) && print) {
            printf("Step t=%d:\n", t);
            printParticles(particle_count, P, min_x, min_y, Nx, Ny);
            printf("\n");
        }
    }

    // ---------- check ----------
    if (print) {
        printf("Final:\n");
        printParticles(particle_count, P, min_x, min_y, Nx, Ny);
        printf("minx: %llu, miny: %llu\n", min_x, min_y);
        printf("maxx: %llu, maxy: %llu\n", max_x, max_y);
        printf("\n");
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

    printf("Verification: %s\n", (success) ? "OK" : "FAILED");

    // ---------- cleanup ----------

    clock_time = clock() - clock_time;
    double time_taken = ((double)clock_time) / CLOCKS_PER_SEC;  // in seconds
    printf("2D_n-body_simulation_seq took %f seconds to execute \n",
           time_taken);

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

void printParticles(int particle_count, struct particle *P, int Mx, int My,
                    int Nx, int Ny) {
    const char *colors = " .-:=+*^X#%@";
    const int numColors = 12;

    // set the 'render' resolution
    int W = RESOLUTION_X;
    int H = RESOLUTION_Y;

    // step size in each dimension
    int xW = (Nx - Mx) / W;
    int yW = (Ny - My) / H;

    Matrix A = createMatrix(W, H);

    for (int i = 0; i < particle_count; i++) {
        int y = MAX(MIN(((int)P[i].pos.y - My) / yW, H - 1), 0);
        int x = MAX(MIN(((int)P[i].pos.x - Mx) / xW, W - 1), 0);
        A[y][x]++;
    }

    // room
    for (int i = 0; i < H; i++) {
        // left wall
        printf("|");
        // actual room
        for (int j = 0; j < W; j++) {
            int c = A[i][j] > numColors ? 11 : A[i][j];
            printf("%c", colors[c]);
        }
        // right wall
        printf("|\n");
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