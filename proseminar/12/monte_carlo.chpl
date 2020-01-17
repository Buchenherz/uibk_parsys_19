use Random;

config const number_of_points : int = 10000000;
config const num_tasks: int = 4;
config const debug: bool = false;

proc main {

    var total_hits: atomic int;
    var rs = new RandomStream(eltType = real, parSafe = true);
    coforall taskID in 1..#num_tasks {
        if debug {
            writeln("Hello from task ", taskID);
        }
        var local_hits: int = 0;
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
        total_hits.add(local_hits);
    }
    var pi:real  = (total_hits.read():real / number_of_points) * 4;
    writeln("Pi: ", pi);
}