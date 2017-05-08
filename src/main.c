#include "auxillary.h"
#include <mpi.h>
#include <stdio.h>

int main(void) {
    double a[6], value;
    int size, rank;

    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    prepare(a, 10, rank);
    value = penta_solve(10, a, rank); 
    printf("I'm here with Id = %d and value = %f\n", rank, value);
    MPI_Finalize();
}
