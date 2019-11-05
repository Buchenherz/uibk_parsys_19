# Assignment 4, due November 6th 2019

The goal of this assignment is to extend the 3D heat stencil application of Assignment 3 to include multiple domain decomposition schemes and investigate the effect of non-blocking communication.

## Exercise 1 - 3d Planes

You can find the comparison to the sequential implementation as well as different amount of machines in the `plane_comparison.pdf` and the different `.csv` documents. If you have tableau installed, the notebook is provided as well.
It becomes clear, that with the parallel implemenation, only half of the time is needed compared to the sequential implementation. Unfortunately, sizes beyond 50x50x50 cubes were not possible because of memory errors from lcc2 (see below). The 4 and 2 System versions are basically equal within low ranges.

### Occured problems

Getting the program to run required to send matrixes of data, which was a struggle to do. Initially, we wanted to use custom data types for that (e.g. a subarray like we had with the 2d heat stencil from week 3); but it just would not work. This is why we switched from custom datatypes, to a different method, for which allocation of matrixes and rooms had to be changed. For a yet unkown reason, the program may crash when the input is not in cube format. It also does not work on 8 machines in parallel.

On lcc2, the 3D planes heat stencil sometimes exits without any errors or with `out of dynamic memory in opal_show_help_yylex()`. This occurs on cube sizes of 60 and above.
