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

void saveTemperature(Room m, int N, int M, int H);

void printTemperature(Room m, int N, int M, int H);

// -- simulation code ---

int main(int argc, char **argv)
{
  // 'parsing' optional input parameter = problem size
  int Nx = 9;  // columns
  int Ny = 11; // rows
  int Nz = 13; // height
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

  MPI_Datatype subarrayType;
  int size[3] = {Nz, Ny, Nx};
  int subSize[3] = {Nz, Ny, Nx};
  int start[3] = {0, 0, 0};
  MPI_Type_create_subarray(3, size, subSize, start, MPI_ORDER_C, MPI_DOUBLE, &subarrayType);
  MPI_Type_commit(&subarrayType);
  // Cartesian communicator
  MPI_Comm comm_cart;
  int ndims = 2;
  int dims[2] = {0, 0};
  int periods[2] = {0, 0};
  int reorder = 0;
  MPI_Dims_create(number_of_ranks, ndims, dims);
  // printf("Dims: %d and %d\n", dims[0],dims[1]);

  MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, periods, reorder, &comm_cart);

  // calc upper / lower rank
  // If value is negative, rank is first or last in topology
  int upper_rank;
  int lower_rank;
  int left_rank;
  int right_rank;

  MPI_Cart_shift(comm_cart, 1, 1, &upper_rank, &lower_rank);
  MPI_Cart_shift(comm_cart, 0, 1, &left_rank, &right_rank);
  printf("Rank %d: Up: %d, Low: %d, Left: %d, Right %d\n", rank, upper_rank, lower_rank, left_rank, right_rank);

  if (rank == 0)
  {
    walltime = clock();
    printf("Computing heat-distribution for room size N=%d, M=%d for T=%d timesteps\n", Nx, Ny, T);
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
  // printTemperature(A, N, M, H);
  printf("\n");

  int local_y_start;
  int local_z_start;
  int local_y_length;
  int local_z_length;
  int local_y_end;
  int local_z_end;

  int coord[2] = {0, 0};
  MPI_Cart_coords(comm_cart, rank, ndims, coord);

  local_y_start = coord[0] * Ny / dims[0];
  local_y_length = Ny / dims[0];
  if (right_rank == -2)
    local_y_length += Ny % dims[0];
  local_y_end = local_y_start + local_y_length;

  local_z_start = coord[1] * Nz / dims[1];
  local_z_length = Nz / dims[1];
  if (lower_rank == -2)
    local_z_length += Nz % dims[1];
  local_z_end = local_z_start + local_z_length;

  printf("Rank %d: #y = %d - %d, #z = %d - %d\n", rank, local_y_start, local_y_end, local_z_start, local_z_end);

  // ---------- compute ----------
  // create a second buffer for the computation
  Room B = createRoom(Nx, Ny, Nz);

  // Set start and end plane for each rank
  // rank 0 starts from 0 and goes to rrs[1] - 1

  Matrix received_upper_plane = createMatrix(Nx, local_y_length);
  Matrix received_lower_plane = createMatrix(Nx, local_y_length);
  Matrix upper_plane = createMatrix(Nx, local_y_length);
  Matrix lower_plane = createMatrix(Nx, local_y_length);
  Matrix received_left_plane = createMatrix(Nx, local_z_length);
  Matrix received_right_plane = createMatrix(Nx, local_z_length);
  Matrix left_plane = createMatrix(Nx, local_z_length);
  Matrix right_plane = createMatrix(Nx, local_z_length);

  // for each time step ..
  for (int t = 0; t < T; t++)
  {
    // send / receive planes here
    // ##########################################
    // lower upper ghostcell exchange
    // ##########################################

    for (int i = 0; i < local_y_length; i++)
    {
      lower_plane[i] = A[local_z_end][local_y_start + i];
      upper_plane[i] = A[local_z_start][local_y_start + i];
    }
    if (upper_rank == lower_rank)
    {
      /* There is no upper / lower rank, both are negative */
    }
    else if (upper_rank < 0)
    {
      // Upper rank does not exist, only receive one from below and send one to below
      MPI_Send(&(lower_plane[0][0]), Nx * local_y_length, MPI_DOUBLE, lower_rank, 0, comm_cart);
      MPI_Recv(&(received_lower_plane[0][0]), Nx * local_y_length, MPI_DOUBLE, lower_rank, 0, comm_cart, MPI_STATUS_IGNORE);

      // Received the last local plane from the lower rank, appending it as the new last local plane
      for (int i = 0; i < local_y_length; i++)
      {
        A[local_z_end + 1][local_y_start + i] = received_lower_plane[i];
      }
    }
    else if (lower_rank < 0)
    {

      // Lower rank does not exist, only receive one from above and send one to above
      MPI_Recv(&(received_upper_plane[0][0]), Nx * local_y_length, MPI_DOUBLE, upper_rank, 0, comm_cart, MPI_STATUS_IGNORE);
      MPI_Send(&(upper_plane[0][0]), Nx * local_y_length, MPI_DOUBLE, upper_rank, 0, comm_cart);

      // Received the first local plane from the upper rank, appending it as the new last local plane
      for (int i = 0; i < local_y_length; i++)
      {
        A[local_z_start - 1][local_y_start + i] = received_upper_plane[i];
      }
    }
    else
    {
      MPI_Request request;
      // Upper and lower ranks exist, send and receive both
      MPI_Isend(&(upper_plane[0][0]), Nx * local_y_length, MPI_DOUBLE, upper_rank, 0, comm_cart, &request);
      MPI_Isend(&(lower_plane[0][0]), Nx * local_y_length, MPI_DOUBLE, lower_rank, 0, comm_cart, &request);
      MPI_Recv(&(received_upper_plane[0][0]), Nx * local_y_length, MPI_DOUBLE, upper_rank, 0, comm_cart, MPI_STATUS_IGNORE);
      MPI_Recv(&(received_lower_plane[0][0]), Nx * local_y_length, MPI_DOUBLE, lower_rank, 0, comm_cart, MPI_STATUS_IGNORE);
      for (int i = 0; i < local_y_length; i++)
      {
        A[local_z_start - 1][local_y_start + i] = received_upper_plane[i];
        A[local_z_end + 1][local_y_start + i] = received_lower_plane[i];
      }
    }

    // ##########################################
    // left right ghostcell exchange
    // ##########################################

    for (int i = 0; i < local_z_length; i++)
    {
      left_plane[i] = A[local_z_start + i][local_y_end];
      right_plane[i] = A[local_z_start + i][local_y_start];
    }
    if (right_rank == left_rank)
    {
    }
    if (left_rank < 0)
    {
      // Left rank does not exist, only receive one from right and send one to right
      MPI_Send(&(right_plane[0][0]), Nx * local_z_length, MPI_DOUBLE, right_rank, 0, comm_cart);
      MPI_Recv(&(received_right_plane[0][0]), Nx * local_z_length, MPI_DOUBLE, right_rank, 0, comm_cart, MPI_STATUS_IGNORE);

      // Received the last local plane from the right rank, appending it as the new last local plane
      for (int i = 0; i < local_z_length; i++)
      {
        A[local_z_start + i][local_y_end + 1] = received_right_plane[i];
      }
    }
    else if (right_rank < 0)
    {
      // Right rank does not exist, only receive one from left and send one to left
      MPI_Recv(&(received_left_plane[0][0]), Nx * local_y_length, MPI_DOUBLE, left_rank, 0, comm_cart, MPI_STATUS_IGNORE);
      MPI_Send(&(left_plane[0][0]), Nx * local_y_length, MPI_DOUBLE, left_rank, 0, comm_cart);

      // Received the first local plane from the upper rank, appending it as the new last local plane
      for (int i = 0; i < local_z_length; i++)
      {
        A[local_z_start + i][local_y_start - 1] = received_left_plane[i];
      }
    }
    else
    {
      MPI_Request request;
      // Upper and lower ranks exist, send and receive both
      MPI_Isend(&(left_plane[0][0]), Nx * local_y_length, MPI_DOUBLE, left_rank, 0, comm_cart, &request);
      MPI_Isend(&(right_plane[0][0]), Nx * local_y_length, MPI_DOUBLE, right_rank, 0, comm_cart, &request);
      MPI_Recv(&(received_left_plane[0][0]), Nx * local_y_length, MPI_DOUBLE, left_rank, 0, comm_cart, MPI_STATUS_IGNORE);
      MPI_Recv(&(received_right_plane[0][0]), Nx * local_y_length, MPI_DOUBLE, right_rank, 0, comm_cart, MPI_STATUS_IGNORE);
      for (int i = 0; i < local_z_length; i++)
      {
        A[local_z_start + i][local_y_start - 1] = received_left_plane[i];
        A[local_z_start + i][local_y_end + 1] = received_right_plane[i];
      }
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

        for (int i = local_z_start; i < local_z_end; i++)
        {
          for (int j = local_y_start; j < local_y_end; j++)
          {
            comb[i][j] = A[i][j];
          }
        }
        // sleep(1);
        for (int orank = 1; orank < number_of_ranks; orank++)
        {

          int orank_cord[2] = {0, 0};
          MPI_Cart_coords(comm_cart, orank, ndims, orank_cord);

          int temp_y_start = orank_cord[0] * Ny / dims[0];
          int temp_y_length = Ny / dims[0];
          if (right_rank == -2)
            temp_y_length += Ny % dims[0];
          int temp_y_end = temp_y_start + temp_y_length;

          int temp_z_start = orank_cord[1] * Nz / dims[1];
          int temp_z_length = Nz / dims[1];
          if (lower_rank == -2)
            temp_z_length += Nz % dims[1];
          int temp_z_end = temp_z_start + temp_z_length;

          // printf("\n");
          Room recv = createRoom(Nx, temp_y_length, temp_z_length);
          // printf("\nRecv array %d\n", orank);
          MPI_Recv(&(recv[0][0][0]), Nx * temp_y_length * temp_z_length, MPI_DOUBLE, orank, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

          // printf("omstart %d, omend %d\n", o_m_start, o_m_end);
          for (int k = temp_z_start; k < temp_z_end; k++)
          {
            for (int g = temp_y_start; g < temp_y_end; g++)
            {
              //printf("r: %f ", recv[k][g]);
              comb[k][g] = recv[k][g];
              // code
            }
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
        Room send_subroom = createRoom(Nx, local_y_length, local_z_length);

        for (int i = 0; i < local_z_length; i++)
        {
          for (int j = 0; j < local_y_length; j++)
          {
            send_subroom[i][j] = A[local_z_start + i][local_y_start + j];
          }
        }
        MPI_Send(&(send_subroom[0][0][0]), Nx * local_y_length * local_z_length, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
        releaseRoom(send_subroom);
      }
    }
  }

  releaseRoom(B);

  // ---------- check ----------
  int success = true;
  if (rank == 0)
  {
    Room comb = createRoom(Nx, Ny, Nz);

    for (int i = local_z_start; i < local_z_end; i++)
    {
      for (int j = local_y_start; j < local_y_end; j++)
      {
        comb[i][j] = A[i][j];
      }
    }
    // sleep(1);
    for (int orank = 1; orank < number_of_ranks; orank++)
    {

      int orank_cord[2] = {0, 0};
      MPI_Cart_coords(comm_cart, orank, ndims, orank_cord);

      int temp_y_start = orank_cord[0] * Ny / dims[0];
      int temp_y_length = Ny / dims[0];
      if (right_rank == -2)
        temp_y_length += Ny % dims[0];
      int temp_y_end = temp_y_start + temp_y_length;

      int temp_z_start = orank_cord[1] * Nz / dims[1];
      int temp_z_length = Nz / dims[1];
      if (lower_rank == -2)
        temp_z_length += Nz % dims[1];
      int temp_z_end = temp_z_start + temp_z_length;

      // printf("\n");
      Room recv = createRoom(Nx, temp_y_length, temp_z_length);
      // printf("\nRecv array %d\n", orank);
      MPI_Recv(&(recv[0][0][0]), Nx * temp_y_length * temp_z_length, MPI_DOUBLE, orank, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

      // printf("omstart %d, omend %d\n", o_m_start, o_m_end);
      for (int k = temp_z_start; k < temp_z_end; k++)
      {
        for (int g = temp_y_start; g < temp_y_end; g++)
        {
          //printf("r: %f ", recv[k][g]);
          comb[k][g] = recv[k][g];
          // code
        }
        // printf("\n");
      }
    }

    A[0][0][0] = comb[0][0][0];
    releaseRoom(comb);

    printf("Final:\n");
    printTemperature(A, Nx, Ny, Nz);
    printf("\n");

    for (long i = 0; i < Nz; i++)
    {
      for (long j = 0; j < Ny; j++)
      {
        for (long k = 0; k < Nx; k++)
        {
          value_t temp = A[i][j][k];
          // printf("%lf ", A[i][j][k]);
          if (273 <= temp && temp <= 273 + 60)
            continue;
          success = 0;
          break;
        }
      }
    }
    walltime = clock() - walltime;
    double time_taken = ((double)walltime) / CLOCKS_PER_SEC; // in seconds
    printf("3d pole MPI heat stencil took %f seconds to execute \n", time_taken);

    printf("Verification: %s\n", (success) ? "OK" : "FAILED");
  }
  else
  {
    // printf("\nSend array %d\n", rank);
    //MPI_Request req;
    Room send_subroom = createRoom(Nx, local_y_length, local_z_length);

    for (int i = 0; i < local_z_length; i++)
    {
      for (int j = 0; j < local_y_length; j++)
      {
        send_subroom[i][j] = A[local_z_start + i][local_y_start + j];
      }
    }
    MPI_Send(&(send_subroom[0][0][0]), Nx * local_y_length * local_z_length, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
    releaseRoom(send_subroom);
  }

  // ---------- cleanup ----------
  releaseMatrix(received_upper_plane);
  releaseMatrix(received_lower_plane);
  releaseMatrix(upper_plane);
  releaseMatrix(lower_plane);
  releaseMatrix(received_left_plane);
  releaseMatrix(received_right_plane);
  releaseMatrix(left_plane);
  releaseMatrix(right_plane);

  releaseRoom(A);
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
