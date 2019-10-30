# Assignment 3, due October 30th 2019

The goal of this assignment is to extend the heat stencil application and measure its performance.

## Exercise 1

This exercise consists in extending the heat stencil application of Assignment 2 to two and three dimensions.

### Tasks

- Extend the heat stencil application to the two- and three-dimensional cases and name them `heat_stencil_2D` and `heat_stencil_3D`.
- Provide sequential and MPI implementations and run them with multiple problem and machine sizes.
  We provide seqential versions of both 2D and 3D as well as the MPI version of 2D. Comparison with sequential and MPI version can be found within the folder `2d heat stencil`.
- How can you verify the correctness of your applications?
  We verify our min max temperatures at the end of the program, where the matrix is recollected. We also compared the sequential and mpi versions and they are the same.

## Exercise 2

This exercise consists in measuring all heat stencil variants (1D, 2D and 3D) to get a grasp of their performance behavior.

### Tasks

- Measure the speedup and efficiency of all three stencil codes for varying problem and machine sizes/mappings. Consider using strong scalability, weak scalability, or both. Justify your choice.
- Illustrate the data in appropriate figures and discuss them. What can you observe?
- Bonus question: Measure and illustrate an application throughput metric. What can you observe?

## General Notes

All the material required by the tasks above (e.g. code, figures, etc...) must be part of the solution that is handed in. Your experiments should be reproducible and comparable to your own measurements using the solution materials that you hand in. For source code, please provide a makefile or other, intuitive means of compiling with the required flags and settings.

**Every** member of your group must be able to explain the given problem, your solution, and possible findings. You may also need to answer detailed questions about any of these aspects.

**Please run any benchmarks or heavy CPU loads only on the compute nodes, not on the login node.**
If you want to do some interactive experimentation, use an _interactive job_ as outlined in the tutorial. Make sure to stop any interactive jobs once you are done.
