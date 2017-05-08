#include "auxillary.h"
#include <math.h>
#include <stdio.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>

double penta_solve(int size, double *data, int rank) {

    double t1, t2, t3, s1, s2, p, x1, x2, x3, x4;
    int first, second, third, fourth, iteration, i, j,
        max_iter, exponent, place, beginning, n;
    double first_data[6], second_data[6], third_data[6], 
        fourth_data[6], temp[6];
    gsl_matrix *m = gsl_matrix_calloc(4,4);
    gsl_vector *x = gsl_vector_calloc(4);
    gsl_vector *b = gsl_vector_calloc(4);
    gsl_permutation *k = gsl_permutation_alloc(4);
    MPI_Status status;
   
    iteration = 0;
    exponent = 1;
    beginning = rank / size;
    place = rank - beginning * size;
    first = place - 2;
    second = place - 1;
    third = place + 1;
    fourth = place + 2;
    max_iter = ((int) log2(size)) - 1;
    for (i=0; i<= max_iter; i++) {
        if (place - exponent < 0) {
            data[1] = 1;
            data[0] = 0;
        }
        else if (place - 2 * exponent < 0) {
            data[0] = 0;
        }
        if (place + exponent > size - 1) {
            data[3] = 1;
            data[4] = 0;
        }
        else if (place + 2 * exponent > size - 1) {
            data[4] = 0;
        }
        if (first >= 0) {
           MPI_Send(&data[0], 6, MPI_DOUBLE, 
                   beginning * size + first, 0, MPI_COMM_WORLD); 
        }
        if (fourth > size - 1) {
            fourth_data[0] = 0;
            fourth_data[1] = 0;
            fourth_data[2] = 1;
            fourth_data[3] = 1;
            fourth_data[4] = 0;
            fourth_data[5] = 0;
        }
        else {
           MPI_Recv(&fourth_data, 6, MPI_DOUBLE,
                   beginning * size + fourth, MPI_ANY_TAG, MPI_COMM_WORLD, &status); 
        }
        if (second >= 0) {
           MPI_Send(&data[0], 6, MPI_DOUBLE,
                   beginning * size + second, 0, MPI_COMM_WORLD); 
        }
        if (third > size - 1) {
            third_data[0] = 0;
            third_data[1] = 0;
            third_data[2] = 1;
            third_data[3] = 1;
            third_data[4] = 0;
            third_data[5] = 0;
        }
        else {
           MPI_Recv(&third_data, 6, MPI_DOUBLE, 
                   beginning * size + third, MPI_ANY_TAG, MPI_COMM_WORLD, &status); 
        }
        if (fourth <= size - 1) {
           MPI_Send(&data[0], 6, MPI_DOUBLE, 
                   beginning * size + fourth, 0, MPI_COMM_WORLD); 
        }
        if (first < 0) {
            first_data[0] = 0;
            first_data[1] = 1;
            first_data[2] = 1;
            first_data[3] = 0;
            first_data[4] = 0;
            first_data[5] = 0;
        }
        else {
           MPI_Recv(&first_data, 6, MPI_DOUBLE, 
                   beginning * size + first, MPI_ANY_TAG, MPI_COMM_WORLD, &status); 
        }
        if (third <= size - 1) {
           MPI_Send(&data[0], 6, MPI_DOUBLE, 
                   beginning * size + third, 0, MPI_COMM_WORLD); 
        }
        if (second < 0) {
            second_data[0] = 0;
            second_data[1] = 1;
            second_data[2] = 1;
            second_data[3] = 0;
            second_data[4] = 0;
            second_data[5] = 0;
        }
        else {
           MPI_Recv(&second_data, 6, MPI_DOUBLE, 
                   beginning * size + second, MPI_ANY_TAG, MPI_COMM_WORLD, &status); 
        }
        /*
        t1 = third_data[4] * fourth_data[1] - 
            third_data[2] * fourth_data[3];
        t2 = third_data[0] * fourth_data[3];
        t3 = first_data[3] * second_data[0];
        s1 = data[1] * t1 + data[3] * t2;
        s2 = first_data[1] * (second_data[4] * 
            data[1] - second_data[2] * data[3]) + 
            t3 * data[3];
        p = t3 * t1 - first_data[1] * (second_data[2] * 
            t1 + second_data[4] * t2); */
        gsl_matrix_set (m, 0, 0, first_data[1]);
        gsl_matrix_set (m, 1, 0, first_data[3]);
        gsl_matrix_set (m, 0, 1, second_data[0]);
        gsl_matrix_set (m, 1, 1, second_data[2]);
        gsl_matrix_set (m, 2, 1, second_data[4]);
        gsl_matrix_set (m, 1, 2, third_data[0]);
        gsl_matrix_set (m, 2, 2, third_data[2]);
        gsl_matrix_set (m, 3, 2, third_data[4]);
        gsl_matrix_set (m, 2, 3, fourth_data[1]);
        gsl_matrix_set (m, 3, 3, fourth_data[3]);
        gsl_vector_set (b, 1, -data[1]);
        gsl_vector_set (b, 2, -data[3]);
        gsl_linalg_LU_decomp(m, k, &j);
        gsl_linalg_LU_solve(m, k, b, x);
/*      x1 = -(second_data[0] * s1) / p;
        x2 = (first_data[1] * s1) / p;
        x3 = (fourth_data[3] * s2) / p;
        x4 = -(third_data[4] * s2) / p; */
        x1 = gsl_vector_get(x, 0);
        x2 = gsl_vector_get(x, 1);
        x3 = gsl_vector_get(x, 2);
        x4 = gsl_vector_get(x, 3);
        temp[0] = x1 * first_data[0];
        temp[1] = x1 * first_data[2] + x2 * second_data[1]
                  + data[0];
        temp[2] = x1 * first_data[4] + x2 * second_data[3]
                  + data[2] + x3 * third_data[1] + 
                  x4 * fourth_data[0];
        temp[3] = data[4] + x3 * third_data[3] + 
                  x4 * fourth_data[2];
        temp[4] = x4 * fourth_data[4];
        temp[5] = x1 * first_data[5] + x2 * second_data[5]
                  + data[5] + x3 * third_data[5] + x4 * 
                  fourth_data[5];
        for (n=0; n<6; n++) {
            data[n] = temp[n];
        }
/*        data = temp; */
        exponent = exponent * 2;
        first = place - 2 * exponent;
        second = place - exponent;
        third = place + exponent;
        fourth = place + 2 * exponent;
    }
    return (data[5]/data[2]);
}


int prepare(double* data, int size, int rank) {

    int beginning, place;

    beginning = rank/size;
    place = rank - beginning * size;
    if (place > 1 & place < size - 2) {
        data[0] = 1;
        data[1] = 26;
        data[2] = 66;
        data[3] = 26;
        data[4] = 1;
    }
    else if (place == 0) {
        data[0] = 0;
        data[1] = 0;
        data[2] = 93;
        data[3] = 26;
        data[4] = 1;
    }
    else if (place == 1) {
        data[0] = 0;
        data[1] = 27;
        data[2] = 66;
        data[3] = 26;
        data[4] = 1;
    }
    else if (place == size - 1) {
        data[0] = 1;
        data[1] = 26;
        data[2] = 93;
        data[3] = 0;
        data[4] = 0;
    }
    else if (place == size - 2) {
        data[0] = 1;
        data[1] = 26;
        data[2] = 66;
        data[3] = 27;
        data[4] = 0;
    }
    data[5] = 120;
    return 0;
}


double transpose_com(double value, int size, int rank) {

    int column, row, target;
    double temp;
    MPI_Status status;
    
    column = rank/size;
    row = rank - column * size;
    target = row * size + column;
    MPI_Sendrecv(&value, 1, MPI_DOUBLE, target, 0, &temp, 1,
            MPI_DOUBLE, target, 0, MPI_COMM_WORLD, &status);
    return temp;
}


    


    

    
        
        
