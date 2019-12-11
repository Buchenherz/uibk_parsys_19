#include <omp.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

typedef double value_t;

#define RESOLUTION 100
// #define VERBOSE

// -- Matrix utilities --

typedef value_t **Matrix;

Matrix createMatrix(int N, int M);

void releaseMatrix(Matrix m);

void printTemperature(Matrix m, int N, int M);

// -- simulation code ---

int main(int argc, char **argv) {
    // 'parsing' optional input parameter = problem size
    int N = 200;  // columns
    int M = 200;  // rows

    if (argc > 1) {
        N = atoi(argv[1]);
        M = atoi(argv[2]);
    }
    int T = 10000;
#ifdef VERBOSE
    printf(
        "Computing heat-distribution for room size N=%d, M=%d for T=%d "
        "timesteps\n",
        N, M, T);
#endif
    double start_time = omp_get_wtime();
    // ---------- setup ----------

    // create a buffer for storing temperature fields
    Matrix A = createMatrix(N, M);

    // set up initial conditions in A
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
            A[i][j] = 273;  // temperature is 0° C everywhere (273 K)
        }
    }

    // and there is a heat source in one corner
    int source_x = N / 4;
    int source_y = M / 8;
    A[source_y][source_x] = 273 + 60;

#ifdef VERBOSE
    printf("Initial:\t");
    printTemperature(A, N, M);
    printf("\n");
#endif

    // ---------- compute ----------

    // create a second buffer for the computation
    Matrix B = createMatrix(N, M);

    // for each time step ..
    for (int t = 0; t < T; t++) {
        // .. we propagate the temperature
#pragma omp parallel for schedule(static)
        for (long long i = 0; i < M; i++) {
            for (long long j = 0; j < N; j++) {
                // center stays constant (the heat is still on)
                // Loop tiling and vectorisation are options one has for maybe optimizing this
                // program
                if (i == source_y && j == source_x) {
                    B[i][j] = A[i][j];
                    continue;
                }

                // get temperature at current position
                value_t tc = A[i][j];

                // get temperatures of adjacent cells
                value_t tl = (j != 0) ? A[i][j - 1] : tc;
                value_t tr = (j != N - 1) ? A[i][j + 1] : tc;
                value_t tu = (i != 0) ? A[i - 1][j] : tc;
                value_t td = (i != M - 1) ? A[i + 1][j] : tc;

                // compute new temperature at current position
                B[i][j] = tc + 0.2 * (tl + tr + tu + td + (-4 * tc));
            }
        }

        // swap matrices (just pointers, not content)
        Matrix H = A;
        A = B;
        B = H;

#ifdef VERBOSE
        // show intermediate step
        if (!(t % 1000)) {
            printf("Step t=%d:\n", t);
            printTemperature(A, N, M);
            printf("\n");
        }
#endif
    }

    releaseMatrix(B);

    // ---------- check ----------
#ifdef VERBOSE
    printf("Final:\n");
    printTemperature(A, N, M);
    printf("\n");
#endif

    int success = 1;
    for (long long i = 0; i < M; i++) {
        for (long long j = 0; j < N; j++) {
            value_t temp = A[i][j];
            if (273 <= temp && temp <= 273 + 60) continue;
            success = 0;
            break;
        }
    }
#ifdef VERBOSE
    printf("Verification: %s\n", (success) ? "OK" : "FAILED");
#endif
    // ---------- cleanup ----------

    releaseMatrix(A);
    double end_time = omp_get_wtime();
    double time_taken = end_time - start_time;
#ifdef PRINT_CSV_HEADER
    printf("Walltime, M, N\n");
#endif
    printf("%f, %d, %d\n", time_taken, M, N);

    // done
    return (success) ? EXIT_SUCCESS : EXIT_FAILURE;
}

Matrix createMatrix(int N, int M) {
    // create data and index Matrix
    value_t **mat = malloc(sizeof(value_t) * M);
    for (int i = 0; i < M; i++) {
        mat[i] = malloc(sizeof(value_t) * N);
    }
    return mat;
}

void releaseMatrix(Matrix m) { free(m); }

void printTemperature(Matrix m, int N, int M) {
    const char *colors = " .-:=+*^X#%@";
    const int numColors = 12;

    // boundaries for temperature (for simplicity hard-coded)
    const value_t max = 273 + 30;
    const value_t min = 273 + 0;

    // set the 'render' resolution
    int W = RESOLUTION;
    int H = RESOLUTION / 4;

    // step size in each dimension
    int xW = N / W;
    int yW = M / H;

    // room
    for (int i = 0; i < H; i++) {
        // left wall
        printf("X");
        // actual room
        for (int j = 0; j < W; j++) {
            // get max temperature in this tile
            value_t max_t = 0;
            for (int y = yW * i; y < yW * i + yW; y++) {
                for (int x = xW * j; x < xW * j + xW; x++) {
                    max_t = (max_t < m[y][x]) ? m[y][x] : max_t;
                }
            }
            value_t temp = max_t;

            // pick the 'color'
            int c = ((temp - min) / (max - min)) * numColors;
            c = (c >= numColors) ? numColors - 1 : ((c < 0) ? 0 : c);

            // print the average temperature
            printf("%c", colors[c]);
        }
        // right wall
        printf("X\n");
    }
}
