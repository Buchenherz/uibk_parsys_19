#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

typedef double value_t;

#define RESOLUTION 120

// -- vector utilities --

typedef value_t *Vector;

Vector createVector(int N);

void releaseVector(Vector m);

void printTemperature(Vector m, int N);

// -- simulation code ---

int main(int argc, char **argv) {
  printf("Hello");
  // region OPENMPI INIT
  MPI_Init(&argc, &argv);  //initialize the MPI environment
  int size;
  MPI_Comm_size(MPI_COMM_WORLD, &size); //get the number of ranks
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank); //get the rankof the caller
  printf("Rank %d of %d online\n", rank, size);
  // endregion OPENMPI INIT

  // 'parsing' optional input parameter = problem size
  int N = 2000;
  if (argc > 1) {
    N = atoi(argv[1]);
  }
  int local_N = N/size;
  int T = N * 500;
  // create a buffer for storing temperature fields
  Vector A;
  
  if (rank == 0) {
    printf("Computing heat-distribution for room size N=%d for T=%d timesteps\n", N, T);

    // ---------- setup ----------
    A = createVector(N);
    // set up initial conditions in A
    for (int i = 0; i < N; i++) {
      A[i] = 273; // temperature is 0° C everywhere (273 K)
    }

    // and there is a heat source in one corner
    int source_x = N / 4;
    A[source_x] = 273 + 60;


    printf("Initial:\t");
    printTemperature(A, N);
    printf("\n");
  }
  // ---------- compute ----------

  // create a second buffer for the computation
  Vector local_A = createVector(local_N+1);
  Vector local_B = createVector(local_N+1);

  // for each time step ..
  for (int t = 0; t < T; t++) {
    MPI_Scatter(A, local_N, MPI_DOUBLE, local_A, local_N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    // .. we propagate the temperature
    for (long long i = 0; i < local_N; i++) {
      // center stays constant (the heat is still on)
      if (local_A[i] == (273 + 60)) {
        local_B[i] = local_A[i];
        continue;
      }

      // get temperature at current position
      value_t tc = local_A[i];

      // get temperatures of adjacent cells
      value_t tl = (i != 0) ? local_A[i - 1] : tc;
      value_t tr = (i != N - 1) ? local_A[i + 1] : tc;

      // compute new temperature at current position
      local_B[i] = tc + 0.2 * (tl + tr + (-2 * tc));
    }

    MPI_Gather(local_B, local_N, MPI_DOUBLE, A, local_N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    // swap matrices (just pointers, not content)
    // Vector H = A;
    // A = local_B;
    // local_B = H;

    // show intermediate step
    if ((rank == 0) && !(t % 1000)) {
      printf("Step t=%d:\t", t);
      printTemperature(A, N);
      printf("\n");
    }
  }

  releaseVector(local_A);
  releaseVector(local_B);

  if (rank == 0) {
    // ---------- check ----------

    printf("Final:\t\t");
    printTemperature(A, N);
    printf("\n");

    int success = 1;
    for (long long i = 0; i < N; i++) {
      value_t temp = A[i];
      if (273 <= temp && temp <= 273 + 60)
        continue;
      success = 0;
      break;
    }

    printf("Verification: %s\n", (success) ? "OK" : "FAILED");
  
    // ---------- cleanup ----------

    releaseVector(A);

    // done
    return (success) ? EXIT_SUCCESS : EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}

Vector createVector(int N) {
  // create data and index vector
  return malloc(sizeof(value_t) * N);
}

void releaseVector(Vector m) { free(m); }

void printTemperature(Vector m, int N) {
  const char *colors = " .-:=+*^X#%@";
  const int numColors = 12;

  // boundaries for temperature (for simplicity hard-coded)
  const value_t max = 273 + 30;
  const value_t min = 273 + 0;

  // set the 'render' resolution
  int W = RESOLUTION;

  // step size in each dimension
  int sW = N / W;

  // room
  // left wall
  printf("X");
  // actual room
  for (int i = 0; i < W; i++) {
    // get max temperature in this tile
    value_t max_t = 0;
    for (int x = sW * i; x < sW * i + sW; x++) {
      max_t = (max_t < m[x]) ? m[x] : max_t;
    }
    value_t temp = max_t;

    // pick the 'color'
    int c = ((temp - min) / (max - min)) * numColors;
    c = (c >= numColors) ? numColors - 1 : ((c < 0) ? 0 : c);

    // print the average temperature
    printf("%c", colors[c]);
  }
  // right wall
  printf("X");
}
