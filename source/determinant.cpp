#include "determinant.h"

/*
 * very simple matrix determinant calculation functions
 * for the algorythm I only need 3*3 and 2*2 matrix determinants and I cannot use libraries since the code is for microcontroller
 */
float simple_2_by_2_matrix_determinant(float matrix[4]) {

    return matrix[0] * matrix[3] - matrix[1] * matrix[2];
}

float simple_3_by_3_matrix_determinant(float matrix[9]) {
    return matrix[0] * (matrix[4] * matrix[8] - matrix[7] * matrix[5]) -
           matrix[1] * (matrix[3] * matrix[8] - matrix[5] * matrix[6]) +
           matrix[2] * (matrix[3] * matrix[7] - matrix[6] * matrix[4]);
}

