This is the solution to the Laplace Equation using Jacobi Iteration written in C
The serial code is present in the file serial.c
The parallel code using OpenMP is present in the file parallel.c
The parallel implementaion of the scheme using OpenMPI is present in the file mpiRun.c

Solution strategy(Jacobi iteration):

1. Apply a square lattice with uniform spacing - label the points (i; j).
2. Apply the xed boundary condition values.
3. Make an initial guess for the interior points.
4. Iterate until convergence, using the "cross" scheme.

