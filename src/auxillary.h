#include <mpi.h>

double penta_solve(int size, double *data, int rank);
int prepare(double* data, int size, int rank);
double transpose_com(double value, int size, int rank);

