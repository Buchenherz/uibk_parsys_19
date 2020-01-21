use Random;
use Time;
use CyclicDist;

// 2552 number of rows and columns
config const rows_and_columns : int = 2552;
// debug mode that enables some useful print statements
config const debug: bool = false;

// To control the number of threads, set CHPL_RT_NUM_THREADS_PER_LOCALE=num_tasks
// (https://chapel-lang.org/docs/usingchapel/tasks.html#controlling-the-number-of-threads)

proc main {
    var timer: Timer;
    timer.start();
    
    var matrixA: [0..#rows_and_columns, 0..#rows_and_columns] real; // first matrix 
    var matrixB: [0..#rows_and_columns, 0..#rows_and_columns] real; // second matrix
    var matrixC: [0..#rows_and_columns, 0..#rows_and_columns] real; // matrix for results
    // use different seeds otherwise matrixA equals matrixB
    var seed1 = 13;
    var seed2 = 103;
    fillRandom(matrixA, seed1);
    fillRandom(matrixB, seed2);
    matrixC = 0.0;
    if debug {
        writeln(matrixA);
        writeln("\n---------------");
        writeln(matrixB);
        writeln("\n---------------");
    }
    const IndexSpace = {0..rows_and_columns-1} dmapped Cyclic(startIdx=0);

    // computation
    forall i in IndexSpace do 
        for (j,k) in {0..rows_and_columns-1,0..rows_and_columns-1} {
            matrixC(i,j) += matrixA(i,k)*matrixB(k,j);
        }
    if debug {
        writeln("\n---------------");
        writeln(matrixC);
        writeln("\n---------------");
    }
    writeln("Walltime (s): ",timer.elapsed());
}
