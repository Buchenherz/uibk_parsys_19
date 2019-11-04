#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <stdbool.h>
#include <time.h>

typedef double value_t;

#define RESOLUTION 60

// -- Matrix utilities --

typedef value_t **Matrix;
typedef value_t ***Room;

clock_t walltime;

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
  int Nx = 10; // columns
  int Ny = 10; // planes
  int Nz = 10; // height
  if (argc > 1)
  {
    Nx = atoi(argv[1]);
    Ny = atoi(argv[2]);
    Nz = atoi(argv[3]);
  }
  int max = (Nx < Ny ? Ny : Nx);
  int T = (max < Nz ? Nz : max) * 500;

  // region OPENMPI INIT
  int number_of_ranks;
  int rank;

  MPI_Init(&argc, &argv);                          //initialize the MPI environment
  MPI_Comm_size(MPI_COMM_WORLD, &number_of_ranks); //get the number of ranks
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);            //get the rankof the caller

  MPI_Datatype plane;
  MPI_Type_vector(Nx, Nz, 1, MPI_DOUBLE, &plane);
  MPI_Type_commit(&plane);

  MPI_Datatype subarrayType;
  int size[3] = {Nz, Ny, Nx};
  int subSize[3] = {Nz, Ny, Nx};
  int start[3] = {0, 0, 0};
  MPI_Type_create_subarray(3, size, subSize, start, MPI_ORDER_C, MPI_DOUBLE, &subarrayType);
  MPI_Type_commit(&subarrayType);
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
    walltime = clock();
    printf("Computing heat-distribution for room size Nx=%d, Ny=%d, Nz=%d for T=%d timesteps\n", Nx, Ny, Nz, T);
  }

  int *rank_plane_sizes = getPlaneSizePerRank(Nx, number_of_ranks);

  for (int i = 0; i < number_of_ranks; i++)
  {
    printf("%d ", rank_plane_sizes[i]);
  }

  // ---------- setup ----------

  // create a buffer for storing temperature fields
  Room A = createRoom(Nx, Ny, Nz);

  // set up initial conditions in A
  for (int i = 0; i < Nx; i++)
  {
    for (int j = 0; j < Ny; j++)
    {
      for (int k = 0; k < Nz; k++)
      {
        A[i][j][k] = 273; // temperature is 0° C everywhere (273 K)
      }
    }
  }

  // and there is a heat source in one corner
  int source_x = Nx / 4;
  int source_y = Ny / 8;
  int source_z = Nz / 2;
  A[source_z][source_y][source_x] = 273 + 60;

  printf("Initial:\t");
  // saveTemperature(A, N, M);
  // printTemperature(A, N, M, H);
  printf("\n");

  // ---------- compute ----------
  // create a second buffer for the computation
  Room B = createRoom(Nx, Ny, Nz);

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
    Matrix received_upper_plane = createMatrix(Nz, Ny);
    Matrix received_lower_plane = createMatrix(Nz, Ny);

    if (upper_rank < 0)
    {
      /* Upper rank does not exist, only receive one from below and send one to below */
      MPI_Send(&(A[local_last_plane_index][0][0]), Nz * Ny, MPI_DOUBLE, lower_rank, 0, comm_cart);
      MPI_Recv(&(received_lower_plane[0][0]), Nz * Ny, MPI_DOUBLE, lower_rank, 0, comm_cart, MPI_STATUS_IGNORE);

      /* Received the last local plane from the lower rank, appending it as the new last local plane */
      A[local_last_plane_index + 1] = received_lower_plane;
    }
    else if (lower_rank < 0)
    {
      /* Lower rank does not exist, only receive one from above and send one to above */
      MPI_Recv(&(received_upper_plane[0][0]), Nz * Ny, MPI_DOUBLE, upper_rank, 0, comm_cart, MPI_STATUS_IGNORE);
      MPI_Send(&(A[local_first_plane_index][0][0]), Nz * Ny, MPI_DOUBLE, upper_rank, 0, comm_cart);

      /* Received the first local plane from the upper rank, appending it as the new last local plane */
      A[local_first_plane_index - 1] = received_upper_plane;
    }
    else
    {
      MPI_Request request;
      /* Upper and lower ranks exist, send and receive both */
      MPI_Isend(&(A[local_first_plane_index][0][0]), Nz * Ny, MPI_DOUBLE, upper_rank, 0, comm_cart, &request);
      MPI_Isend(&(A[local_last_plane_index][0][0]), Nz * Ny, MPI_DOUBLE, lower_rank, 0, comm_cart, &request);
      MPI_Recv(&(received_upper_plane[0][0]), Nz * Ny, MPI_DOUBLE, upper_rank, 0, comm_cart, MPI_STATUS_IGNORE);
      MPI_Recv(&(received_lower_plane[0][0]), Nz * Ny, MPI_DOUBLE, lower_rank, 0, comm_cart, MPI_STATUS_IGNORE);
      A[local_first_plane_index - 1] = received_upper_plane;
      A[local_last_plane_index + 1] = received_lower_plane;
    }

    // .. we propagate the temperature
    for (long long i = 0; i < Nx; i++)
    {
      for (long long j = 0; j < Ny; j++)
      {
        for (long long k = 0; k < Nz; k++)
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
          value_t tl = (k != 0) ? A[i][j][k - 1] : tc;       // left
          value_t tr = (k != Nx - 1) ? A[i][j][k + 1] : tc;  // right
          value_t tu = (j != 0) ? A[i][j - 1][k] : tc;       // up
          value_t td = (j != Ny - 1) ? A[i][j + 1][k] : tc;  // down
          value_t th = (i != 0) ? A[i - 1][j][k] : tc;       // higher
          value_t tlo = (i != Nz - 1) ? A[i + 1][j][k] : tc; // lower

          // compute new temperature at current position
          B[i][j][k] = tc + 0.1 * (tl + tr + tu + td + th + tlo + (-6 * tc));
        }
      }
    }

    // swap matrices (just pointers, not content)
    Room Help = A;
    A = B;
    B = Help;

    // show intermediate step
    if (!(t % 1000))
    {
      if (rank == 0)
      {
        Room comb = createRoom(Nx, Ny, Nz);

        for (size_t i = 0; i < local_m_end; i++)
        {
          comb[i] = A[i];
        }
        // sleep(1);
        for (int orank = 1; orank < number_of_ranks; orank++)
        {
          // printf("\n");
          Room recv = createRoom(Nx, Ny, Nz);
          // printf("\nRecv array %d\n", orank);
          MPI_Recv(&(recv[0][0][0]), Nx * Nz * Ny, MPI_DOUBLE, orank, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

          int sum = 0;
          for (int j = 0; j < orank; j++)
          {
            sum = sum + rank_plane_sizes[j];
          }
          int o_m_start = sum;
          int o_m_end = sum + rank_plane_sizes[orank];
          // printf("omstart %d, omend %d\n", o_m_start, o_m_end);
          for (int k = o_m_start; k < o_m_end; k++)
          {
            comb[k] = recv[k];
            // printf("\n");
          }
        }

        printTemperature(comb, Nx, Ny, Nz);
        releaseRoom(comb);
      }
      else
      {
        // printf("\nSend array %d\n", rank);
        //MPI_Request req;
        MPI_Send(&(A[0][0][0]), Nx * Nz * Ny, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
      }
    }
  }

  releaseRoom(B);

  // ---------- check ----------
  int success = true;
  if (rank == 0)
  {
    Room comb = createRoom(Nx, Ny, Nz);

    for (size_t i = 0; i < local_m_end; i++)
    {
      comb[i] = A[i];
    }
    // sleep(1);
    for (int orank = 1; orank < number_of_ranks; orank++)
    {
      // printf("\n");
      Room recv = createRoom(Nx, Ny, Nz);
      // printf("\nRecv array %d\n", orank);
      MPI_Recv(&(recv[0][0][0]), Nx * Ny * Nz, MPI_DOUBLE, orank, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

      int sum = 0;
      for (int j = 0; j < orank; j++)
      {
        sum = sum + rank_plane_sizes[j];
      }
      int o_m_start = sum;
      int o_m_end = sum + rank_plane_sizes[orank];
      printf("omstart %d, omend %d\n", o_m_start, o_m_end);
      for (int k = o_m_start; k < o_m_end; k++)
      {
        comb[k] = recv[k];
        // printf("\n");
      }
    }

    A = comb;

    printf("Final:\n");
    printTemperature(comb, Nx, Ny, Nz);
    printf("\n");

    for (long i = 0; i < Nx; i++)
    {
      for (long j = 0; j < Ny; j++)
      {
        for (long k = 0; k < Nz; k++)
        {
          value_t temp = comb[i][j][k];
          // printf("%lf ", A[i][j][k]);
          if (273 <= temp && temp <= 273 + 60)
            continue;
          // printf("[%ld,%ld,%ld]", i, j, k);
          success = 0;
          break;
        }
      }
    }
    walltime = clock() - walltime;
    double time_taken = ((double)walltime) / CLOCKS_PER_SEC; // in seconds
    printf("3d plane MPI heat stencil took %f seconds to execute \n", time_taken);

    printf("Verification: %s\n", (success) ? "OK" : "FAILED");

    // releaseRoom(comb);
  }
  else
  {
    // printf("\nSend array %d\n", rank);
    //MPI_Request req;
    MPI_Send(&(A[0][0][0]), Nx * Ny * Nz, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
  }

  // ---------- cleanup ----------

  releaseRoom(A);
  free(rank_plane_sizes);
  MPI_Type_free(&plane);
  MPI_Type_free(&subarrayType);
  MPI_Finalize();

  // done
  return (success) ? EXIT_SUCCESS : EXIT_FAILURE;
}

Room createRoom(int N, int M, int H)
{
  // https: //stackoverflow.com/questions/39108092/allocating-contiguous-memory-for-a-3d-array-in-c
  // create data and index Matrix
  value_t ***room = (value_t ***)malloc(sizeof(value_t **) * H);
  value_t *data = (value_t *)malloc(sizeof(value_t) * N * M * H);
  for (int i = 0; i < H; i++)
  {
    room[i] = (value_t **)malloc(sizeof(value_t *) * M);
    for (int j = 0; j < M; j++)
    {
      int idx = N * j + i * N * M;
      room[i][j] = &data[idx];
    }
  }
  return room;
}

Matrix createMatrix(int N, int M)
{
  // https://stackoverflow.com/questions/5901476/sending-and-receiving-2d-array-over-mpi
  // create data and index Matrix
  value_t **mat = (value_t **)malloc(sizeof(value_t *) * M);
  value_t *data = (value_t *)malloc(N * M * sizeof(value_t));
  for (int i = 0; i < M; i++)
  {
    mat[i] = &(data[N * i]);
  }
  return mat;
}

void releaseMatrix(Matrix m)
{
  free(m[0]);
}
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
  int Yres = RESOLUTION / 4;
  int Zres = RESOLUTION / 20;
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
      if (j == 0 || j == Xres - 1)
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
      if (j == 0 || j == Xres - 1)
        printf("+");
      else
        printf("-");
    }
    printf("\n");
  }
}
