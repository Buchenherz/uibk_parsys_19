# Assignment 11, due January 15th 2020

The goal of this assignment is to parallelize an unknown application using profiling and OpenMP.

## Exercise 1

### Description

The file [real.tar.gz](real.tar.gz) contains a realistic implementation of a (simple) numerical algorithm. Imagine you are tasked with making this implementation faster by parallelizing it with OpenMP, without any further information.

### Tasks

- Familiarize yourself with the code. You are not required to look at every source line, but rather profile the code using the means discussed in the lecture and get a grasp on its computational hotspots and performance characteristics (computation-heavy, memory-heavy, etc.).
- Investigate any loops that carry larger workloads and determine if and how they can be parallelized. Parallelize them with OpenMP. Ensure that any code modification does not violate program correctness with respect to its output.
- Benchmark the original, sequential program and your parallelized version for 1, 2, 4 and 8 threads on LCC2 and enter your results in [this table](https://docs.google.com/spreadsheets/d/1hLTIc-VlzBOBrlZY2cSt1RIKc376UYyOLge2QcnJ7sQ/edit?usp=sharing).

## Task 1 
After decompressing the tar.gz file, we build the project using the delivered Makefile. Our first step was to analyse hotspots using the `google-perftools` profiler (since neither gprof nor prof are available for macOS (it just works)). We followed the steps taught in the lecture to gather the following results: 
 
 | No. profiling samples |  % profiling samples in function |  % profiling samples in functions printed so far |  No. profiling samples in this function and its callees |  % profiling samples in function and its callees |  Function name  | 
 |-----------------------|----------------------------------|--------------|---------------------|--------------|-----------------| 
 | 1662                  |  49.2%                           |  49.2%       |  1675               |  49.6%       |  resid          | 
 | 943                   |  27.9%                           |  77.2%       |  951                |  28.2%       |  psinv          | 
 | 287                   |  8.5%                            |  85.7%       |  287                |  8.5%        |  interp         | 
 | 199                   |  5.9%                            |  91.6%       |  199                |  5.9%        |  rprj3          | 
 | 135                   |  4.0%                            |  95.6%       |  135                |  4.0%        |  vranlc         | 
 | 40                    |  1.2%                            |  96.7%       |  40                 |  1.2%        |  zero3          | 
 | 33                    |  1.0%                            |  97.7%       |  66                 |  2.0%        |  norm2u3        | 
 | 21                    |  0.6%                            |  98.3%       |  21                 |  0.6%        |  __log10_finite | 
 | 21                    |  0.6%                            |  99.0%       |  21                 |  0.6%        |  comm3          | 
 | 19                    |  0.6%                            |  99.5%       |  166                |  4.9%        |  zran3          | 
 
 This clearly demonstrates that most of the execution time is spend in both `resid` and `psinv`, so that is where are parallelization focus needs to be.


## General Notes

All the material required by the tasks above (e.g. code, figures, etc...) must be part of the solution that is handed in. Your experiments should be reproducible and comparable to your own measurements using the solution materials that you hand in. For source code, please provide a makefile or other, intuitive means of compiling with the required flags and settings.

**Every** member of your group must be able to explain the given problem, your solution, and possible findings. You may also need to answer detailed questions about any of these aspects.

**Please run any benchmarks or heavy CPU loads only on the compute nodes, not on the login node.**
If you want to do some interactive experimentation, use an *interactive job* as outlined in the tutorial. Make sure to stop any interactive jobs once you are done.
