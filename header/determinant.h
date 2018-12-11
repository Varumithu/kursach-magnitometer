#ifndef DETERMINANT_H
#define DETERMINANT_H

/*
 * very simple matrix determinant calculation functions
 * for the algorythm I only need 3*3 and 2*2 matrix determinants and I cannot use libraries since the code is for microcontroller
 * I do not cheeck for size; if it is too little, std::vector will throw, if too large, i do not care
 */
double simple_2_by_2_matrix_determinant(double matrix[4]);
double simple_3_by_3_matrix_determinant(double matrix[9]);

#endif // DETERMINANT_H
