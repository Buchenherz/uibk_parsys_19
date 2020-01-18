use Random;
use Time;

// 10^9 number of points
config const rows_and_columns : uint = 2552;
// number of tasks / threads that are started
config const num_tasks: uint = 4;
// debug mode that enables some useful print statements
config const debug: bool = false;

// To control the number of threads, set CHPL_RT_NUM_THREADS_PER_LOCALE=num_tasks
// (https://chapel-lang.org/docs/usingchapel/tasks.html#controlling-the-number-of-threads)

proc main {
    var timer: Timer;
    timer.start();
    
    // Total hits, local hits will be reduced to this
    var total_hits: atomic uint;

    var matrixA = Matrix(#rows_and_columns);
    var matrixB = Matrix(#rows_and_columns);
    var matrixC = Matrix(#rows_and_columns);

    fillRandom(matrixA, seed);
    fillRandom(matrixB, seed);
    if debug {
        writeln(matrixA);
        writeln("\n---------------");
    }

    // coforall starts a new task each iteration 
    // coforall loops, the execution of the
    // parent task will not continue until all the children sync up.
    // for more info, visit https://chapel-lang.org/docs/users-guide/taskpar/coforall.html?highlight=coforall
    coforall taskID in 0..#num_tasks {
        if debug {
            writeln("Starting task ", taskID);
        }
        
        for i in rows_and_columns/num_tasks*taskID..rows_and_columns/num_tasks*(taskID+1) {
            for j in rows_and_columns/num_tasks*taskID..rows_and_columns/num_tasks*(taskID+1) {
                for k in rows_and_columns/num_tasks*taskID..rows_and_columns/num_tasks*(taskID+1) {
                    matrixC[i][j] += matrixA[i][k]*matrixB[k][j];
                }
            }
        }
    }
    if debug {
        writeln(matrixC);
        writeln("\n---------------");
        writeln("Time (s)");
    }
    writeln(timer.elapsed());
}
