#include "auxillary.h"
#include <mpi.h>
#include <stdio.h>

int main(void) {
    double a[6], value;
    int size, rank;

    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
/*    if (rank == 0) {
        a[0] = 0;
        a[1] = 0;
        a[2] = 92;
        a[3] = 26;
        a[4] = 1;
        a[5] = 120;
    }test
    else if (rank == 1) {
        a[0] = 0;
        a[1] = 27;
        a[2] = 66;
        a[3] = 26;
        a[4] = 1;
        a[5] = 120;
    }
    else if (rank == 6) {
        a[0] = 1;
        a[1] = 26;
        a[2] = 66;
        a[3] = 27;
        a[4] = 0;
        a[5] = 120;
    }
    else if (rank == 7) {
        a[0] = 1;
        a[1] = 26;
        a[2] = 92;
        a[3] = 0;
        a[4] = 0;
        a[5] = 120;
    }
    else {
        a[0] = 1;
        a[1] = 26;
        a[2] = 66;
        a[3] = 26;
        a[4] = 1;
        a[5] = 120;
    } */
    prepare(a, 10, rank);
    value = penta_solve(10, a, rank); 
    a[5] = 120 * transpose_com(value, 10, rank);
    value = penta_solve(10, a, rank); 
    printf("I'm here with Id = %d and value = %f\n", rank, value);
    MPI_Finalize();
}
     

