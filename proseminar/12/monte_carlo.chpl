use Random;
use Time;

// 10^9 number of points
config const number_of_points : uint = 1000000000;
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
    // Chapels version of initialising randomness, setting parSafe to false to 
    // avoid locking overhead.
    var rs = new RandomStream(eltType = real, parSafe = false);

    // coforall starts a new task each iteration 
    // coforall loops, the execution of the
    // parent task will not continue until all the children sync up.
    // for more info, visit https://chapel-lang.org/docs/users-guide/taskpar/coforall.html?highlight=coforall
    coforall taskID in 1..#num_tasks {
        if debug {
            writeln("Starting task ", taskID);
        }
        var local_hits: uint = 0;
        for i in 1..number_of_points/num_tasks {
            var random_x = rs.getNext(0,1);
            var random_y = rs.getNext(0,1);
            if (random_x ** 2 + random_y ** 2) < 1.0 {
                local_hits+=1;
            }
        }
        if debug {
            writeln("Task ", taskID, " got ", local_hits, " local hits");
        }
        // add local hits atomically to total hits
        total_hits.add(local_hits);
    }
    var pi:real  = (total_hits.read():real / number_of_points) * 4;
    if debug {
        writeln("\n---------------");
        writeln("Pi, Time (s)");
    }
    writeln(pi, ", ", timer.elapsed());
}