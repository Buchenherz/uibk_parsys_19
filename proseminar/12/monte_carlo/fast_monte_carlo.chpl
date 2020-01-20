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

const PI25DT = 3.141592653589793238462643;

/*
* main procedure
*/
proc main() {
    var timer : Timer;
    // domain over the number of random points to generate
    var D = {1..numOfPoints};

    var rs = new RandomStream(eltType = real, parSafe = false);

    timer.start();
    /* Run the Monte Carlo method using data parallel reduction to compute count.
     * The reduction is over parallel loop iterating in zippered manner over
     * a tuple x,y generated randomly by RandomStream object.
     * https://chapel-lang.org/docs/modules/standard/Random/PCGRandom.html#PCGRandom.RandomStream
     * Esentially, rs.iterate(D, real) generates an array with D random variables in it. We generate two of these
     * arrays and zip them using zip. It is well explained using python: https://www.programiz.com/python-programming/methods/built-in/zip
     */
    var globalCount = + reduce [(x,y) in zip(rs.iterate(D, real), rs.iterate(D, real))] ((x**2 + y**2) <= 1.0);

    var pi = 4.0 * globalCount / numOfPoints;
    var wallTime = timer.elapsed();

    writeln("Pi, error, walltime, number of tasks");
    writeln(pi, ", ", abs(pi - PI25DT), ", ", wallTime, ", ", numOfTasks);
}