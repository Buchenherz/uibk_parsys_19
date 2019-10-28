#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <mpi.h>
#include <unistd.h>

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
  int N = 10; // columns
  int M = 10; // rows
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

  MPI_Datatype row;
  MPI_Type_vector(N, 1, 1, MPI_DOUBLE, &row);
  MPI_Type_commit(&row);

  MPI_Datatype myType;
  int size[2] = {M, N};
  int subSize[2] = {M, N};
  int start[2] = {0, 0};
  MPI_Type_create_subarray(2, size, subSize, start, MPI_ORDER_C, MPI_DOUBLE, &myType);
  MPI_Type_commit(&myType);
  // Cartesian communicator
  MPI_Comm comm_cart;
  int ndims = 2;
  int dims[2] = {1, number_of_ranks};
  int periods[2] = {0, 0};
  int reorder = 0;

  MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, periods, reorder, &comm_cart);

  // calc upper / lower rank
  // If value is negative, rank is first or last in topology
  int upper_rank;
  int lower_rank;

  MPI_Cart_shift(comm_cart, 1, 1, &upper_rank, &lower_rank);

  printf("rank %d, %d, %d\n", rank, upper_rank, lower_rank);

  printf("Rank %d of %d online\n", rank, number_of_ranks - 1);

  // endregion

  if (rank == 0)
  {
    printf("Computing heat-distribution for room size N=%d, M=%d for T=%d timesteps\n", N, M, T);
  }

  int *rank_row_sizes = getRowSizePerRank(M, number_of_ranks);

  for (size_t i = 0; i < number_of_ranks; i++)
  {
    printf("%d ", rank_row_sizes[i]);
  }

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

  // printf("Initial:\t");
  // printTemperature(A, N, M);
  // printf("\n");

  // ---------- compute ----------

  // create a second buffer for the computation
  Matrix B = createMatrix(N, M);

  // Set start and end row for each rank
  // rank 0 starts from 0 and goes to rrs[1] - 1
  int local_m_start;
  int local_m_end;
  if (rank == 0)
  {
    local_m_start = 0;
    local_m_end = rank_row_sizes[0];
  }
  else
  {
    int sum = 0;
    for (int i = 0; i < rank; i++)
    {
      sum = sum + rank_row_sizes[i];
    }
    local_m_start = sum;
    local_m_end = sum + rank_row_sizes[rank];
  }

  int local_first_row_index = local_m_start;
  int local_last_row_index = local_m_end - 1;

  printf("LFRI %d, LLRI %d\n", local_first_row_index, local_last_row_index);
  printf("LMS %d, LME %d\n", local_m_start, local_m_end);
  // for each time step ..
  for (int t = 0; t < T; t++)
  {
    // send / receive rows here
    // rows that will be added on top / below local rows
    double received_upper_row[N];
    double received_lower_row[N];

    if (upper_rank < 0)
    {
      /* Upper rank does not exist, only receive one from below and send one to below */
      MPI_Send(A[local_last_row_index], 1, row, lower_rank, 0, comm_cart);
      MPI_Recv(received_lower_row, 1, row, lower_rank, 0, comm_cart, MPI_STATUS_IGNORE);

      /* Received the last local row from the lower rank, appending it as the new last local row */
      A[local_last_row_index + 1] = received_lower_row;
      // for (int index = 0; index < N; index++)
      // {
      //   printf("Recv: %f, ", received_lower_row[index]);
      // }
      // printf("\n");
    }
    else if (lower_rank < 0)
    {
      /* Lower rank does not exist, only receive one from above and send one to above */
      MPI_Recv(received_upper_row, 1, row, upper_rank, 0, comm_cart, MPI_STATUS_IGNORE);
      MPI_Send(A[local_first_row_index], 1, row, upper_rank, 0, comm_cart);

      /* Received the first local row from the upper rank, appending it as the new last local row */
      A[local_first_row_index - 1] = received_upper_row;

      // for (int index = 0; index < N; index++)
      // {
      //   printf("Recv: %f, ", received_upper_row[index]);
      // }
      // printf("\n");
    }
    else
    {
      MPI_Request request;
      /* Upper and lower ranks exist, send and receive both */
      MPI_Isend(A[local_first_row_index], 1, row, upper_rank, 0, comm_cart, &request);
      MPI_Isend(A[local_last_row_index], 1, row, lower_rank, 0, comm_cart, &request);
      MPI_Recv(received_upper_row, 1, row, upper_rank, 0, comm_cart, MPI_STATUS_IGNORE);
      MPI_Recv(received_lower_row, 1, row, lower_rank, 0, comm_cart, MPI_STATUS_IGNORE);
      A[local_first_row_index - 1] = received_upper_row;
      A[local_last_row_index + 1] = received_lower_row;
      // for (int index = 0; index < N; index++)
      // {
      //   printf("Recv: %f, ", received_upper_row[index]);
      // }
      // printf("\n");
      // for (int index = 0; index < N; index++)
      // {
      //   printf("Recv: %f, ", received_lower_row[index]);
      // }
      // printf("\n");
    }

    // .. we propagate the temperature
    for (long long i = local_m_start; i < local_m_end; i++)
    {

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
    if (!(t % 1000) && false)
    {
      if (rank == 0)
      {

        Matrix comb = createMatrix(N, M);

        for (int i = 0; i < local_m_end; i++)
        {
          comb[i] = A[i];
        }
        for (int orank = 1; orank < number_of_ranks; orank++)
        {
          Matrix recv = createMatrix(N, M);
          printf("\nRecv array %d\n", orank);
          MPI_Recv(&(recv[0][0]), N * M, MPI_DOUBLE, orank, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
          int sum = 0;
          for (int j = 0; j < orank; j++)
          {
            sum = sum + rank_row_sizes[j];
          }
          int o_m_start = sum;
          int o_m_end = sum + rank_row_sizes[orank];
          printf("\n omstart %d, omend %d\n", o_m_start, o_m_end);
          for (int k = o_m_start; k < o_m_end; k++)
          {
            comb[k] = recv[k];
          }
          free(recv);
        }

        printf("Step t=%d:\n", t);
        printTemperature(comb, N, M);
        printf("\n");
        free(comb);
      }
      else
      {
        // printf("\nSend array %d\n", rank);
        // MPI_Request req;
        MPI_Send(&(A[0][0]), N * M, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
        //MPI_Wait(&req, MPI_STATUS_IGNORE);
      }
    }
  }
  releaseMatrix(B);
  // ---------- check ----------
  int success = 1;
  if (rank == 0)
  {

    Matrix comb = createMatrix(N, M);
    // printf("Rank 0 rows:\n");
    for (size_t i = 0; i < local_m_end; i++)
    {
      comb[i] = A[i];
    }
    // sleep(1);
    for (int orank = 1; orank < number_of_ranks; orank++)
    {
      printf("\n");
      Matrix recv = createMatrix(N, M);
      printf("\nRecv array %d\n", orank);
      MPI_Recv(&recv[0][0], 1, myType, orank, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

      int sum = 0;
      for (int j = 0; j < orank; j++)
      {
        sum = sum + rank_row_sizes[j];
      }

      int o_m_start = sum;
      int o_m_end = sum + rank_row_sizes[orank];
      printf("omstart %d, omend %d\n", o_m_start, o_m_end);
      for (int k = o_m_start; k < o_m_end; k++)
      {
        comb[k] = recv[k];
        for (size_t g = 0; g < N; g++)
        {
          printf("%f ", recv[k][g]);
          /* code */
        }
        printf("\n");
      }
      releaseMatrix(recv);
      printf("Freed matrix recv\n");
    }

    printf("Final\n");
    //printTemperature(comb, N, M);
    printf("\n");

    for (long long i = 0; i < M; i++)
    {
      for (long long j = 0; j < N; j++)
      {
        value_t temp = comb[i][j];
        printf("%f ", temp);
        if (273 <= temp && temp <= 273 + 60)
          continue;
        success = 0;
        // printf("temp: %f, i: %lld, j: %lld \n", temp, i, j);
        break;
      }
      printf("\n");
    }

    printf("Verification: %s\n", (success) ? "OK" : "FAILED");
    releaseMatrix(comb);
  }
  else
  {
    // printf("\nSend array %d\n", rank);
    //MPI_Request req;
    // for (size_t k = 0; k < M; k++)
    // {
    //   for (size_t a = 0; a < N; a++)
    //   {
    //     printf("%f ", A[k][a]);
    //     /* code */
    //   }
    //   printf("\n");
    //   /* code */
    // }
    // printf("\n");
    // Matrix send = createMatrix(N, M);
    // printf("Rank 0 rows:\n");
    // for (size_t i = local_m_start; i < local_m_end; i++)
    // {
    //   send[i] = A[i];
    //   for (size_t k = 0; k < N; k++)
    //   {
    //     printf("%f ", send[i][k]);
    //     /* code */
    //   }
    //   printf("\n");
    // }

    value_t array[M][N];
    for (size_t i = 0; i < M; i++)
    {
      for (size_t k = 0; k < N; k++)
      {
        array[i][k] = A[i][k];
      }
    }

    MPI_Send(array, 1, myType, 0, 1, MPI_COMM_WORLD);
    printf("Freed matrix send\n");

    //MPI_Wait(&req, MPI_STATUS_IGNORE);
  }

  // ---------- cleanup ----------

  releaseMatrix(A);
  free(rank_row_sizes);
  MPI_Type_free(&row);
  MPI_Type_free(&myType);
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

void releaseMatrix(Matrix m)
{
  free(m);
}

int *getRowSizePerRank(int M, int number_of_ranks)
{
  int whole_div = (int)M / number_of_ranks;
  int *rows_per_rank;
  rows_per_rank = malloc(number_of_ranks * sizeof(int));

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
