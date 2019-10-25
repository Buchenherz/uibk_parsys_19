#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <mpi.h>

typedef double value_t;

/***
GENERAL CONCEPT
The general concept of this program is to split the work onto different ranks.
For that, we initialise each rank with the full matrix. For the calculation,
we give each rank a specific number of rows to work with. At the beginning of each
round of calculations, the first and last row of each section is shared with
the top and bottom ranks. After that, the sections are combined and send to rank 0
for printing.
***/

#define RESOLUTION 100

// -- Matrix utilities --

typedef value_t **Matrix;

Matrix createMatrix(int N, int M);

void releaseMatrix(Matrix m);

void printTemperature(Matrix m, int N, int M);

// Returns an int array that stores the rows that will be processed per rank
int *getRowSizePerRank(int M, int number_of_ranks);

// -- simulation code ---

int main(int argc, char **argv)
{
  // 'parsing' optional input parameter = problem size
  int N = 300; // columns
  int M = 100; // rows
  if (argc > 1)
  {
    N = atoi(argv[1]);
    M = atoi(argv[2]);
  }
  int T = (N < M ? M : N) * 500;

  // region OPENMPI INIT
  int number_of_ranks;
  int rank;
  MPI_Init(&argc, &argv);                          //initialize the MPI environment
  MPI_Comm_size(MPI_COMM_WORLD, &number_of_ranks); //get the number of ranks
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);            //get the rankof the caller
  printf("Rank %d of %d online\n", rank, number_of_ranks);
  // endregion

  if (rank == 0)
  {
    printf("Computing heat-distribution for room size N=%d, M=%d for T=%d timesteps\n", N, M, T);
  }

  int *rank_row_sizes = getRowSizePerRank(M, number_of_ranks);

  // ---------- setup ----------

  // create a buffer for storing temperature fields
  Matrix A = createMatrix(N, M);

  // set up initial conditions in A
  for (int i = 0; i < M; i++)
  {
    for (int j = 0; j < N; j++)
    {
      A[i][j] = 273; // temperature is 0° C everywhere (273 K)
    }
  }

  // and there is a heat source in one corner
  int source_x = N / 4;
  int source_y = M / 8;
  A[source_y][source_x] = 273 + 60;

  printf("Initial:\t");
  printTemperature(A, N, M);
  printf("\n");

  // ---------- compute ----------

  // create a second buffer for the computation
  Matrix B = createMatrix(N, M);

  // Set start and end row for each rank
  // rank 0 starts from 0 and goes to rrs[1] - 1
  int local_m_start = 0;
  if (rank != 0)
  {
    local_m_start = rank_row_sizes[rank];
  }
  int sum = 0;
  for (int i = 0; i < rank; i++)
  {
    sum = sum + rank_row_sizes[i];
  }
  int local_m_end = sum + rank_row_sizes[rank];

  // for each time step ..
  for (int t = 0; t < T; t++)
  {
    // .. we propagate the temperature
    for (long long i = local_m_start; i < local_m_end; i++)
    {
      // send rows here
      for (long long j = 0; j < N; j++)
      {
        // center stays constant (the heat is still on)
        if (i == source_y && j == source_x)
        {
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

    // show intermediate step
    if (!(t % 1000))
    {
      printf("Step t=%d:\n", t);
      printTemperature(A, N, M);
      printf("\n");
    }
  }

  releaseMatrix(B);

  // ---------- check ----------

  printf("Final:\n");
  printTemperature(A, N, M);
  printf("\n");

  int success = 1;
  for (long long i = 0; i < M; i++)
  {
    for (long long j = 0; j < N; j++)
    {
      value_t temp = A[i][j];
      if (273 <= temp && temp <= 273 + 60)
        continue;
      success = 0;
      break;
    }
  }

  printf("Verification: %s\n", (success) ? "OK" : "FAILED");

  // ---------- cleanup ----------

  releaseMatrix(A);
  free(rank_row_sizes);

  MPI_Finalize(); //cleanup

  // done
  return (success) ? EXIT_SUCCESS : EXIT_FAILURE;
}

Matrix createMatrix(int N, int M)
{
  // create data and index Matrix
  value_t **mat = malloc(sizeof(value_t) * M);
  for (int i = 0; i < M; i++)
  {
    mat[i] = malloc(sizeof(value_t) * N);
  }
  return mat;
}

void releaseMatrix(Matrix m) { free(m); }

int *getRowSizePerRank(int M, int number_of_ranks)
{
  int whole_div = (int)M / number_of_ranks;
  int *rows_per_rank;
  rows_per_rank = (int *)malloc(number_of_ranks * sizeof(int));

  for (int loop = number_of_ranks; loop >= 0; loop--)
  {
    rows_per_rank[loop] = 0;
  }
  int array_sum = 0;
  int i = 0;
  while (true)
  {
    array_sum = 0;
    if (i > number_of_ranks - 1)
    {
      i = 0;
    }
    for (int loop = number_of_ranks; loop >= 0; loop--)
    {
      array_sum = array_sum + rows_per_rank[loop];
    }
    if (array_sum >= M)
    {
      break;
    }
    rows_per_rank[i] += 1;

    i += 1;
  }
  return rows_per_rank;
}

void printTemperature(Matrix m, int N, int M)
{
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
  for (int i = 0; i < H; i++)
  {
    // left wall
    printf("X");
    // actual room
    for (int j = 0; j < W; j++)
    {
      // get max temperature in this tile
      value_t max_t = 0;
      for (int y = yW * i; y < yW * i + yW; y++)
      {
        for (int x = xW * j; x < xW * j + xW; x++)
        {
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
