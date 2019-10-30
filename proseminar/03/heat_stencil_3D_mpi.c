#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

typedef double value_t;

#define RESOLUTION 60

// -- Matrix utilities --

typedef value_t **Matrix;
typedef value_t ***Room;

Matrix createMatrix(int N, int M);
Room createRoom(int N, int M, int H);

void releaseMatrix(Matrix m);
void releaseRoom(Room m);

// Returns an int array that stores the planes that will be processed per rank
int *getPlaneSizePerRank(int H, int number_of_ranks);

void saveTemperature(Room m, int N, int M, int H);

void printTemperature(Room m, int N, int M, int H);

// -- simulation code ---

int main(int argc, char **argv)
{
  // 'parsing' optional input parameter = problem size
  int N = 96; // columns
  int M = 48; // planes
  int H = 12; // height
  if (argc > 1)
  {
    N = atoi(argv[1]);
    M = atoi(argv[2]);
    H = atoi(argv[3]);
  }
  int max = (N < M ? M : N);
  int T = (max < H ? H : max) * 500;

  // region OPENMPI INIT
  int number_of_ranks;
  int rank;

  MPI_Init(&argc, &argv);                          //initialize the MPI environment
  MPI_Comm_size(MPI_COMM_WORLD, &number_of_ranks); //get the number of ranks
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);            //get the rankof the caller

  MPI_Datatype plane;
  MPI_Type_vector(M, N, 1, MPI_DOUBLE, &plane);
  MPI_Type_commit(&plane);

  MPI_Datatype myType;
  int size[3] = {H, M, N};
  int subSize[3] = {H, M, N};
  int start[3] = {0, 0, 0};
  MPI_Type_create_subarray(3, size, subSize, start, MPI_ORDER_C, MPI_DOUBLE, &myType);
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

  if (rank == 0)
  {
    printf("Computing heat-distribution for room size N=%d, M=%d for T=%d timesteps\n", N, M, T);
  }

  int *rank_plane_sizes = getPlaneSizePerRank(H, number_of_ranks);

  for (int i = 0; i < number_of_ranks; i++)
  {
    printf("%d ", rank_plane_sizes[i]);
  }

  // ---------- setup ----------

  // create a buffer for storing temperature fields
  Room A = createRoom(N, M, H);

  // set up initial conditions in A
  for (int i = 0; i < H; i++)
  {
    for (int j = 0; j < M; j++)
    {
      for (int k = 0; k < N; k++)
      {
        A[i][j][k] = 273; // temperature is 0° C everywhere (273 K)
      }
    }
  }

  // and there is a heat source in one corner
  int source_x = N / 4;
  int source_y = M / 8;
  int source_z = H / 2;
  A[source_z][source_y][source_x] = 273 + 60;

  printf("Initial:\t");
  // saveTemperature(A, N, M);
  printTemperature(A, N, M, H);
  printf("\n");

  // ---------- compute ----------
  // create a second buffer for the computation
  Room B = createRoom(N, M, H);

  // Set start and end plane for each rank
  // rank 0 starts from 0 and goes to rrs[1] - 1
  int local_m_start;
  int local_m_end;
  if (rank == 0)
  {
    local_m_start = 0;
    local_m_end = rank_plane_sizes[0];
  }
  else
  {
    int sum = 0;
    for (int i = 0; i < rank; i++)
    {
      sum = sum + rank_plane_sizes[i];
    }
    local_m_start = sum;
    local_m_end = sum + rank_plane_sizes[rank];
  }

  int local_first_plane_index = local_m_start;
  int local_last_plane_index = local_m_end - 1;

  printf("LFRI %d, LLRI %d\n", local_first_plane_index, local_last_plane_index);
  printf("LMS %d, LME %d\n", local_m_start, local_m_end);
  // for each time step ..
  for (int t = 0; t < T; t++)
  {
    // send / receive planes here
    // planes that will be added on top / below local planes
    Matrix received_upper_plane = createMatrix(N,M);
    Matrix received_lower_plane = createMatrix(N,M);

    if (upper_rank < 0)
    {
      /* Upper rank does not exist, only receive one from below and send one to below */
      MPI_Send(A[local_last_plane_index], 1, plane, lower_rank, 0, comm_cart);
      MPI_Recv(received_lower_plane, 1, plane, lower_rank, 0, comm_cart, MPI_STATUS_IGNORE);

      /* Received the last local plane from the lower rank, appending it as the new last local plane */
      A[local_last_plane_index + 1] = received_lower_plane;
      for (int i = 0; i < M; i++)
      {
        for (int j = 0; j < N; j++)
        {
          printf("Recv: %lf, ", received_lower_plane[i][j]);
        }
        printf("\n");
      }
      printf("\n");
    }
    else if (lower_rank < 0)
    {
      /* Lower rank does not exist, only receive one from above and send one to above */
      MPI_Recv(received_upper_plane, 1, plane, upper_rank, 0, comm_cart, MPI_STATUS_IGNORE);
      MPI_Send(A[local_first_plane_index], 1, plane, upper_rank, 0, comm_cart);

      /* Received the first local plane from the upper rank, appending it as the new last local plane */
      A[local_first_plane_index - 1] = received_upper_plane;
      for (int i = 0; i < M; i++)
      {
        for (int j = 0; j < N; j++)
        {
          printf("Recv: %lf, ", received_upper_plane[i][j]);
        }
        printf("\n");
      }
      printf("\n");
    }
    else
    {
      MPI_Request request;
      /* Upper and lower ranks exist, send and receive both */
      MPI_Isend(A[local_first_plane_index], 1, plane, upper_rank, 0, comm_cart, &request);
      MPI_Isend(A[local_last_plane_index], 1, plane, lower_rank, 0, comm_cart, &request);
      MPI_Recv(received_upper_plane, 1, plane, upper_rank, 0, comm_cart, MPI_STATUS_IGNORE);
      MPI_Recv(received_lower_plane, 1, plane, lower_rank, 0, comm_cart, MPI_STATUS_IGNORE);
      A[local_first_plane_index - 1] = received_upper_plane;
      A[local_last_plane_index + 1] = received_lower_plane;
      for (int i = 0; i < M; i++)
      {
        for (int j = 0; j < N; j++)
        {
          printf("Recv: %lf, ", received_lower_plane[i][j]);
        }
        printf("\n");
      }
      printf("\n");
      for (int i = 0; i < M; i++)
      {
        for (int j = 0; j < N; j++)
        {
          printf("Recv: %lf, ", received_upper_plane[i][j]);
        }
        printf("\n");
      }
      printf("\n");
    }
    releaseMatrix(received_upper_plane);
    releaseMatrix(received_lower_plane);
    // .. we propagate the temperature
    for (long long i = 0; i < H; i++)
    {
      for (long long j = 0; j < M; j++)
      {
        for (long long k = 0; k < N; k++)
        {
          // center stays constant (the heat is still on)
          if (i == source_z && j == source_y && k == source_x)
          {
            B[i][j][k] = A[i][j][k];
            continue;
          }

          // get temperature at current position
          value_t tc = A[i][j][k];

          // get temperatures of adjacent cells
          value_t tl = (k != 0) ? A[i][j][k - 1] : tc;      // left
          value_t tr = (k != N - 1) ? A[i][j][k + 1] : tc;  // right
          value_t tu = (j != 0) ? A[i][j - 1][k] : tc;      // up
          value_t td = (j != M - 1) ? A[i][j + 1][k] : tc;  // down
          value_t th = (i != 0) ? A[i - 1][j][k] : tc;      // higher
          value_t tlo = (i != H - 1) ? A[i + 1][j][k] : tc;  // lower

          // compute new temperature at current position
          B[i][j][k] = tc + 0.1 * (tl + tr + tu + td +  th +  tlo + (-6 * tc));
        }
      }
    }

    // swap matrices (just pointers, not content)
    Room Help = A;
    A = B;
    B = Help;

    // // show intermediate step
    // if (!(t % 1000))
    // {
    //   if (rank == 0) {
    //     Room comb = createRoom(N, M, H);

    //     for (int i = 0, i_rank = -1; i < H; i++)
    //     {
    //       if (i % local_m_end) {
    //         i_rank++;
    //       }
    //       if (i < local_m_end) {
    //         comb[i] = A[i];
    //       }
    //       else {
    //         MPI_Recv(comb[i], 1, plane, i_rank , 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    //       }
    //     }
    //     printTemperature(comb, N, M, H);
    //     releaseRoom(comb);
    //   }
    //   else
    //   {
    //     // printf("\nSend array %d\n", rank);
    //     //MPI_Request req;
    //     for (int i = local_m_start; i < local_m_end; i++) {
    //       MPI_Send(A[i], 1, plane, 0, 1, MPI_COMM_WORLD);
    //     //MPI_Wait(&req, MPI_STATUS_IGNORE);
    //     }
    //   }
    // }
  }

  releaseRoom(B);

  // ---------- check ----------
  int success = 1;
  if (rank == 0) {
    // Room comb = createRoom(N, M, H);

    // for (int i = 0, i_rank = -1; i < H; i++)
    // {
    //   if (i % local_m_end) {
    //     i_rank++;
    //   }
    //   if (i < local_m_end) {
    //     comb[i] = A[i];
    //   }
    //   else {
    //     MPI_Recv(comb[i], 1, plane, i_rank , 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    //   }
    // }

    printf("Final:\n");
    printTemperature(A, N, M, H);
    printf("\n");

    for (long i = 0; i < H; i++)
    {
      for (long j = 0; j < M; j++)
      {
        for (long k = 0; k < N; k++)
        {
          value_t temp = A[i][j][k];
          if (273 <= temp && temp <= 273 + 60)
            continue;
          success = 0;
          printf("[%ld,%ld,%ld]",i,j,k);
          break;
        }
      }
    }

    printf("Verification: %s\n", (success) ? "OK" : "FAILED");

    // releaseRoom(comb);
  }
//   else
//   {
//     // printf("\nSend array %d\n", rank);
//     //MPI_Request req;
//     for (int i = local_m_start; i < local_m_end; i++)
//       MPI_Send(A[i], 1, plane, 0, 1, MPI_COMM_WORLD);
//     //MPI_Wait(&req, MPI_STATUS_IGNORE);
//   }

  // ---------- cleanup ----------

  releaseRoom(A);

  // done
  return (success) ? EXIT_SUCCESS : EXIT_FAILURE;
}

Room createRoom(int N, int M, int H)
{
  // create data and index Matrix
  value_t ***room = malloc(sizeof(value_t) * H);
  for (int i = 0; i < H ; i++)
  {
    room[i] = malloc(sizeof(value_t) * M);
    for (int j = 0; j < M ; j++)
    {
      room[i][j] = malloc(sizeof(value_t) * N);
    }
  }
  return room;
}

Matrix createMatrix(int N, int M)
{
  // create data and index Matrix
  value_t **mat = malloc(sizeof(value_t) * M);
  for (int i = 0; i < M ; i++)
  {
    mat[i] = malloc(sizeof(value_t) * N);
  }
  return mat;
}

void releaseMatrix(Matrix m) { free(m); }
void releaseRoom(Room m) { free(m); }

int *getPlaneSizePerRank(int H, int number_of_ranks)
{
  int whole_div = (int)H / number_of_ranks;
  int *planes_per_rank;
  planes_per_rank = (int *)malloc(number_of_ranks * sizeof(int));

  for (int loop = 0; loop < number_of_ranks; loop++)
  {
    planes_per_rank[loop] = whole_div;
    if (loop == 0) 
      planes_per_rank[loop] += H % number_of_ranks;
  }
  return planes_per_rank;
}

void printTemperature(Room m, int N, int M, int H)
{
  const char *colors = " .-:=+*^X#%@";
  const int numColors = 12;

  // boundaries for temperature (for simplicity hard-coded)
  const value_t max = 273 + 30;
  const value_t min = 273 + 0;

  // set the 'render' resolution
  int Xres = RESOLUTION;
  int Yres = RESOLUTION/4;
  int Zres = RESOLUTION/20;
  // printf("Resolution = %d x %d x %d\n",Xres,Yres,Zres);
  // step size in each dimension
  int xW = N / Xres;
  int yW = M / Yres;
  int zW = H / Zres;

  // room
  for (int i = 0; i < Zres; i++)
  {
    // top wall
    for (int j = 0; j < Xres; j++)
    {
      if (j == 0 || j == Xres-1)
        printf("+");
      else
        printf("-");
    }
    printf("\n");
    for (int j = 0; j < Yres; j++)
    {
      // left wall
      printf("|");
      // actual room
      for (int k = 0; k < Xres; k++)
      {
        // get max temperature in this tile
        value_t max_t = 0;
        for (int z = zW * i; z < zW * i + zW; z++)
        {
          for (int y = yW * j; y < yW * j + yW; y++)
          {
            for (int x = xW * k; x < xW * k + xW; x++)
            {
              max_t = (max_t < m[z][y][x]) ? m[z][y][x] : max_t;
              // printf("%d,%d,%d", z,y,x);
            }
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
      printf("|\n");
    }
    // bottom wall
    for (int j = 0; j < Xres; j++)
    {
      if (j == 0 || j == Xres-1)
        printf("+");
      else
        printf("-");
    }
    printf("\n");
  }
}
