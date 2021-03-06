#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdbool.h>
#include <string.h>

typedef double value_t;

#define RESOLUTION 100

// -- Matrix utilities --

typedef value_t **Matrix;

Matrix createMatrix(int rows, int columns);

Matrix randomlyFillMatrix(Matrix m, int rows, int columns, bool print);

void printMatrix(Matrix m, int rows, int columns);

void releaseMatrix(Matrix m);

int rand_int(int min, int max);
// -- simulation code ---

int main(int argc, char **argv)
{
  // 'parsing' optional input parameter = problem size
  int N = 10; // columns
  int M = 100; // rows
  int L = 10; // rows
  clock_t start_time;
  bool print = false;

  if (argc > 3) {
    N = atoi(argv[1]);
    M = atoi(argv[2]);
    L = atoi(argv[3]);
    if (argc > 4)
      print = strcmp("print",argv[4]) == 0;
  }else {
    printf("usage: matrix_mul_seq <int Matrix1 rows> <int Matrix1 columns && Matrix2 rows> <int Matrix2 columns> <string print>\n");
  }
  
  start_time = clock();
  // ---------- setup ----------
  srand(1234);   // Initialization, should only be called once.

  // create matrixes
  Matrix A = createMatrix(N, M);
  Matrix B = createMatrix(M, L);
  Matrix C = createMatrix(N, L);

  A = randomlyFillMatrix(A, N, M, print);
  B = randomlyFillMatrix(B, M, L, print);

  // ---------- compute ----------
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < L; j++) {
      for (int k = 0; k < M; k++) {
        C[i][j] += A[i][k]*B[k][j];
      }
    }
  }
  if (print) {
    printMatrix(C,N,L);
  }

  // ---------- cleanup ----------
  releaseMatrix(A);
  releaseMatrix(B);
  releaseMatrix(C);

  start_time = clock() - start_time;
  double time_taken = ((double)start_time) / CLOCKS_PER_SEC; // in seconds
  printf("matrix multiplication took %f seconds to execute \n", time_taken);

  // done
  return EXIT_SUCCESS;
}

Matrix createMatrix(int rows, int columns)
{
  // create data and index Matrix
  value_t **mat = calloc(rows, sizeof(value_t));
  for (int i = 0; i < rows; i++) {
    mat[i] = calloc(columns, sizeof(value_t));
  }
  return mat;
}

Matrix randomlyFillMatrix(Matrix m, int rows, int columns, bool print) {
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < columns; j++) {
      m[i][j] = rand_int(0,100);
    }
  }
  if (print) {
    printMatrix(m,rows,columns);
  }
  return m;
}

void printMatrix(Matrix m, int rows, int columns)
{
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < columns; j++) {
      printf("%.2lf\t", m[i][j]);
    }
    printf("\n");
  }
}

void releaseMatrix(Matrix m) { free(m); }

// return a uniformly distributed random value
int rand_int(int min, int max) { return (rand() % (max - min)) + min; }

