/**
* @file  fast_monte_carlo.chpl
*
* @desription Calculation of number PI using Monte carlo method with
*             data parallelism model using reduction
*
* @source https://github.com/stefkeve/chapel_cray_imi_examples/blob/3f8a69dc183b4722f99eb1117eb4c667e394e1b4/03_montecarlo/pi_monte_carlo_dp.chpl
*
* @arguments --numOfPoints - max number of points
* @arguments --numOfTasks - max number of tasks to be used for reduction
*/

use Random;
use Time;

config const numOfTasks  : int = here.maxTaskPar;
config const numOfPoints : int = 1000000000; // 10^9
config const debug : bool = false;

const PI25DT = 3.141592653589793238462643;

/*
* main procedure
*/
proc main() {
    var timer : Timer;
    timer.start();
    var total_hits: atomic int = 0;
    if debug {
        writeln("Global variable total_hits is stored on ", total_hits.locale.id);
        writeln("Global variable numOfPoints is stored on ", numOfPoints.locale.id);
    }

    for loc in Locales do
        on loc do {
            if debug {
                writeln("Locale ", here.id, " starting with ", numOfTasks, " available tasks!");
                writeln("Locale ", here.id, " generating ", numOfPoints/numLocales, " points.");
            }

            // domain over the number of random points to generate
            var D = {1..numOfPoints/numLocales};

            var rs = new RandomStream(eltType = real, parSafe = false);

            /* Run the Monte Carlo method using data parallel reduction to compute count.
            * The reduction is over parallel loop iterating in zippered manner over
            * a tuple x,y generated randomly by RandomStream object.
            * https://chapel-lang.org/docs/modules/standard/Random/PCGRandom.html#PCGRandom.RandomStream
            * Esentially, rs.iterate(D, real) generates an array with D random variables in it. We generate two of these
            * arrays and zip them using zip. It is well explained using python: https://www.programiz.com/python-programming/methods/built-in/zip
            */
            var local_hits = + reduce [(x,y) in zip(rs.iterate(D, real), rs.iterate(D, real))] ((x**2 + y**2) <= 1.0);
            if debug {
                writeln("Locale ", here.id, " adding ", local_hits, " local hits to global hits");
            }
            total_hits.add(local_hits);
        }

    if debug {
        writeln("Total hits: ", total_hits.read():real);
        writeln("Num of points: ", numOfPoints);
    }

    var pi = (4.0 * total_hits.read():real) / numOfPoints;
    var wallTime = timer.elapsed();

    writeln("Pi, error, walltime, number of tasks, number of locales");
    writeln(pi, ", ", abs(pi - PI25DT), ", ", wallTime, ", ", numOfTasks, ", ", numLocales);
}