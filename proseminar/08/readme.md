# Assignment 8, due December 4th 2019

The goal of this assignment is to parallelize several applications with OpenMP.

## Exercise 1

### Tasks

- Use OpenMP to parallelize the Monte Carlo computation of π of Assignment 2 and the 2D heat stencil simulation of Assignment 3.
- Measure the execution time of your OpenMP programs for several problem sizes and for 1 to 8 threads.
- Illustrate the data in appropriate speedup/efficiency figures and discuss them. What can you observe?
- Try to maximize their performance by considering all sequential and parallelism-related optimizations we discussed so far. Which did you choose and why?

#### Monte Carlo Obersvations

Interestingly, both the manual splitting and the automatic openmp splitting are slower with multiple threads. Up until 10000000 elements, every version takes arount the same time, with the sequential one being the fastest still. The 6 core version was the slowest when increasing the number of elements. Looking at the relative speedup it seems that s
ome strange things are going on during the calculation. Every single multithreaded execution is slower. Executions using 2-4 threads are massively varying in speedup / slowdown with increasing problem size, suggesting that running the jobs at the same time lowers their performance siginificantly.

#### Monte Carlo Optimizations

On the first parallel run of the program, we restarted several threads because we always divided the workload manually into 8 chunks, like in this solution: http://jakascorner.com/blog/2016/05/omp-monte-carlo-pi.html. On the second run we simplified the program to a single loop and let openmp do the splitting for us.

#### 2D Heat Stencil Obersvations

In this exercise we can clearly see the improvements of using multiple threads. Geneerally, the more threads we use, the better absolute time we get. In the instance of a 5000x5000 space, we gained a speedup of 4x when using 8 threads in parallel compared to the sequential version. Interestingly, we get a spke at 1000x1000 space when using more than 4 threads.

## Exercise 2

### Tasks

- Use OpenMP to develop a parallel matrix multiplication program.
- Measure the execution time of your OpenMP programs for several matrix sizes and for 1 to 8 threads.
- Illustrate the data in appropriate speedup/efficiency figures and discuss them. What can you observe?
- Try to maximize the performance by considering all sequential and parallelism-related optimizations we discussed so far. Which did you choose and why?

## General Notes

All the material required by the tasks above (e.g. code, figures, etc...) must be part of the solution that is handed in. Your experiments should be reproducible and comparable to your own measurements using the solution materials that you hand in. For source code, please provide a makefile or other, intuitive means of compiling with the required flags and settings.

**Every** member of your group must be able to explain the given problem, your solution, and possible findings. You may also need to answer detailed questions about any of these aspects.

**Please run any benchmarks or heavy CPU loads only on the compute nodes, not on the login node.**
If you want to do some interactive experimentation, use an _interactive job_ as outlined in the tutorial. Make sure to stop any interactive jobs once you are done.
