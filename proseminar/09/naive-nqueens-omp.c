/* C/C++ program to solve N Queen Problem using
backtracking Original code by
https://www.geeksforgeeks.org/n-queen-problem-backtracking-3/.
Modified to allow multiple solution lookup.
There is also a trick with matrix as function parameter passing that works quite
well here:
https://stackoverflow.com/questions/18661702/passing-matrix-as-a-parameter-in-function
*/

#include <omp.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// #define DEBUG true

int min(int a, int b) {
    return a < b ? a : b;
}

/* A utility function to print solution */
void printSolution(int N, int board[N][N]) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) printf(" %d ", board[i][j]);
        printf("\n");
    }
}

/* A utility function to check if a queen can
be placed on board[row][col]. Note that this
function is called when "col" queens are
already placed in columns from 0 to col -1.
So we need to check only left side for
attacking queens */
/* there are no good ways to parllize this -
it takes much more time if you parllize all loops with a reduction of a boolean -
there would be too much unnecessary computation*/
bool isSafe(int row, int col, int N, int board[N][N]) {
    int i, j;

    /* Check this row on left side */
    for (i = 0; i < col; i++)
        if (board[row][i]) return false;

    /* Check upper diagonal on left side */
    for (i = row, j = col; i >= 0 && j >= 0; i--, j--)
        if (board[i][j]) return false;

    /* Check lower diagonal on left side */
    for (i = row, j = col; j >= 0 && i < N; i++, j--)
        if (board[i][j]) return false;

    return true;
}

/* A recursive utility function to solve N
Queen problem */
void solveNQUtil(int col, int N, int board[N][N], int* size) {
    /* base case: If all queens are placed
    then return true */
    if (col >= N) {
        #pragma omp atomic
        *size += 1;
#ifdef DEBUG
        if (size == 1) {
            // This will print the first possible solution
            printSolution(N, board);
        }
#endif
        return;
    }
    /* Consider this column and try placing
    this queen in all rows one by one */
    for (int i = 0; i < N; i++) {
        /* Check if the queen can be placed on
        board[i][col] */
        if (isSafe(i, col, N, board)) {
            /* Place this queen in board[i][col] */
            board[i][col] = 1;

            /* recur to place rest of the queens */
            solveNQUtil(col + 1, N, board, size);

            /* If placing queen in board[i][col]
            doesn't lead to a solution, then
            remove queen from board[i][col] */
            board[i][col] = 0;  // BACKTRACK
        }
    }
}

void baseSolveNQUtil(int N, int* size) {
    int board[N][N];
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            board[i][j] = 0;
        }
    }
    /* Consider this column and try placing
    this queen in all rows one by one */
    for (int i = 0; i < N; i++) {
        board[i][0] = 1;
        
        /* recur to place rest of the queens */
        #pragma omp task
        solveNQUtil(1, N, board, size);

        /* If placing queen in board[i][col]
        doesn't lead to a solution, then
        remove queen from board[i][col] */
        board[i][0] = 0;  // BACKTRACK
    }
}

/* This function solves the N Queen problem using
Backtracking. It mainly uses solveNQUtil() to
solve the problem. It returns false if queens
cannot be placed, otherwise, return true and
prints placement of queens in the form of 1s.
Please note that there may be more than one
solutions, this function prints one of the
feasible solutions.*/
bool solveNQ(int N) {
    // Total amount of possible solutions
    int* size = malloc(sizeof(int));
    *size = 0;
    double start_time = omp_get_wtime();
    #pragma omp parallel
    #pragma omp single
    baseSolveNQUtil(N, size);
    double end_time = omp_get_wtime();
    if (size == 0) {
        printf("%d, %f\n", *size, end_time - start_time);
        return false;
    } else {
        printf("%d, %f\n", *size, end_time - start_time);
        return true;
    }
}

// driver program to test above function
int main(int argc, char *argv[]) {
    int N;
    if (argc == 2) {
        N = atoi(argv[1]);
    } else {
        printf("Usage: ./program <Nqueens>\n");
        return EXIT_FAILURE;
    }
    solveNQ(N);
    return EXIT_SUCCESS;
}
