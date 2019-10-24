#include <stdio.h>
#include <stdlib.h>

typedef double value_t;

#define RESOLUTION 60

// -- Matrix utilities --

typedef value_t ***Matrix;

Matrix createMatrix(int N, int M, int H);

void releaseMatrix(Matrix m);

void saveTemperature(Matrix m, int N, int M, int H);

void printTemperature(Matrix m, int N, int M, int H);

// -- simulation code ---

int main(int argc, char **argv)
{
  // 'parsing' optional input parameter = problem size
  int N = 96; // columns
  int M = 48; // rows
  int H = 12; // height
  if (argc > 1)
  {
    N = atoi(argv[1]);
    M = atoi(argv[2]);
    H = atoi(argv[3]);
  }
  int max = (N < M ? M : N);
  int T = (max < H ? H : max) * 500;
  printf("Computing heat-distribution for room size N=%d, M=%d for T=%d timesteps\n", N, M, T);

  // ---------- setup ----------

  // create a buffer for storing temperature fields
  Matrix A = createMatrix(N, M, H);

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
  Matrix B = createMatrix(N, M, H);

  // for each time step ..
  for (int t = 0; t < T; t++)
  {
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
    Matrix Help = A;
    A = B;
    B = Help;

    // show intermediate step
    if (!(t % 1000))
    {
      printf("Step t=%d:\n", t);
      // saveTemperature(A, N, M);
      printTemperature(A, N, M, H);
      printf("\n");
    }
  }

  releaseMatrix(B);

  // ---------- check ----------

  printf("Final:\n");
  // saveTemperature(A, N, M, H);
  printTemperature(A, N, M, H);
  printf("\n");

  int success = 1;
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

  // ---------- cleanup ----------

  releaseMatrix(A);

  // done
  return (success) ? EXIT_SUCCESS : EXIT_FAILURE;
}

Matrix createMatrix(int N, int M, int H)
{
  // create data and index Matrix
  value_t ***mat = malloc(sizeof(value_t) * H);
  for (int i = 0; i < H ; i++)
  {
    mat[i] = malloc(sizeof(value_t) * M);
    for (int j = 0; j < M ; j++)
    {
      mat[i][j] = malloc(sizeof(value_t) * N);
    }
  }
  return mat;
}

void releaseMatrix(Matrix m) { free(m); }

void saveTemperature(Matrix m, int N, int M, int H) {
}

void printTemperature(Matrix m, int N, int M, int H)
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
