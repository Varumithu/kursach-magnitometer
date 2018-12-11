//
// Created by valmit on 23/11/18.
//

#include "magnitometer.h"
#include "determinant.h"

#include <cmath>
#include <algorithm>


// The constructor. Calls calculate_coeffs to calculate x_0, ... , beta
transform_calculator::transform_calculator(const float array[12]) {
    for (size_t i = 0; i < 4; ++i) {
        base_points[i] = point(array[i * 3], array[i * 3 + 1], array[i * 3 + 2]);
    }
    calculate_coeffs();
}

// A helper function to declare matrix and calculate determinant. Only a separate function to avoid cluttering calculate_coeffs
float transform_calculator::calculate_det_A_1 () const {
    float A_1[9]{base_points[0].x * base_points[0].x - base_points[1].x * base_points[1].x,
           base_points[0].y - base_points[1].y,
           base_points[0].z - base_points[1].z,
           base_points[0].x * base_points[0].x - base_points[2].x * base_points[2].x,
           base_points[0].y - base_points[2].y,
           base_points[0].z - base_points[2].z,
           base_points[0].x * base_points[0].x - base_points[3].x * base_points[3].x,
           base_points[0].y - base_points[3].y,
           base_points[0].z - base_points[3].z};
    return simple_3_by_3_matrix_determinant(A_1);
}

// A helper function to declare matrix and calculate determinant. Only a separate function to avoid cluttering calculate_coeffs
float transform_calculator::calculate_det_B_1 () const {
    float B_1[9] {base_points[1].y * (base_points[0].y - base_points[1].y),
           base_points[0].y - base_points[1].y,
           base_points[0].z - base_points[1].z,
           base_points[2].y * (base_points[0].y - base_points[2].y),
           base_points[0].y - base_points[2].y,
           base_points[0].z - base_points[2].z,
           base_points[3].y * (base_points[0].y - base_points[3].y),
           base_points[0].y - base_points[3].y,
           base_points[0].z - base_points[3].z};
    return simple_3_by_3_matrix_determinant(B_1);
}

// A helper function to declare matrix and calculate determinant. Only a separate function to avoid cluttering calculate_coeffs
float transform_calculator::calculate_det_C_1 () const {
    float C_1[9] =  {base_points[1].z * (base_points[0].z - base_points[1].z),
           base_points[0].y - base_points[1].y,
           base_points[0].z - base_points[1].z,
           base_points[2].z * (base_points[0].z - base_points[2].z),
           base_points[0].y - base_points[2].y,
           base_points[0].z - base_points[2].z,
           base_points[3].z * (base_points[0].z - base_points[3].z),
           base_points[0].y - base_points[3].y,
           base_points[0].z - base_points[3].z};
    return simple_3_by_3_matrix_determinant(C_1);
}

// A helper function to declare matrix and calculate determinant. Only a separate function to avoid cluttering calculate_coeffs
float transform_calculator::calculate_det_A_2 () const {
    float A_2[9] {base_points[0].x - base_points[1].x,
           base_points[1].x * (base_points[0].x - base_points[1].x),
           base_points[0].z - base_points[1].z,
           base_points[0].x - base_points[2].x,
           base_points[2].x * (base_points[0].x - base_points[2].x),
           base_points[0].z - base_points[2].z,
           base_points[0].x - base_points[3].x,
           base_points[3].x * (base_points[0].x - base_points[3].x),
           base_points[0].z - base_points[3].z};
    return simple_3_by_3_matrix_determinant(A_2);
}

// A helper function to declare matrix and calculate determinant. Only a separate function to avoid cluttering calculate_coeffs
float transform_calculator::calculate_det_B_2 () const {
    float B_2[9] {base_points[0].x - base_points[1].x,
           base_points[0].y * base_points[0].y - base_points[1].y * base_points[1].y,
           base_points[0].z - base_points[1].z,
           base_points[0].x - base_points[2].x,
           base_points[0].y * base_points[0].y - base_points[2].y * base_points[2].y,
           base_points[0].z - base_points[2].z,
           base_points[0].x - base_points[3].x,
           base_points[0].y * base_points[0].y - base_points[3].y * base_points[3].y,
           base_points[0].z - base_points[3].z};
    return simple_3_by_3_matrix_determinant(B_2);
}

// A helper function to declare matrix and calculate determinant. Only a separate function to avoid cluttering calculate_coeffs
float transform_calculator::calculate_det_C_2 () const {
    float C_2[9] {base_points[0].x - base_points[1].x,
           base_points[1].z * (base_points[0].z - base_points[1].z),
           base_points[0].z - base_points[1].z,
           base_points[0].x - base_points[2].x,
           base_points[2].x * (base_points[0].x - base_points[2].x),
           base_points[0].z - base_points[2].z,
           base_points[0].x - base_points[3].x,
           base_points[3].x * (base_points[0].x - base_points[3].x),
           base_points[0].z - base_points[3].z};
    return simple_3_by_3_matrix_determinant(C_2);
}

// A helper function to declare matrix and calculate determinant. Only a separate function to avoid cluttering calculate_coeffs
float transform_calculator::calculate_det_delta() const {
    float m_delta[9] {(base_points[0].x - base_points[1].x),
               (base_points[0].y - base_points[1].y),
               (base_points[0].z - base_points[1].z),
               (base_points[0].x - base_points[2].x),
               (base_points[0].y - base_points[2].y),
               (base_points[0].z - base_points[2].z),
               (base_points[0].x - base_points[3].x),
               (base_points[0].y - base_points[3].y),
               (base_points[0].z - base_points[3].z)};
    return simple_3_by_3_matrix_determinant(m_delta);
}

// A helper function to declare matrix and calculate determinant. Only a separate function to avoid cluttering calculate_coeffs
float transform_calculator::calculate_delta_z (const float a, const float b, const float c) const {
    float delta_z[9] {a,
               (base_points[0].y - base_points[1].y),
               (base_points[0].z - base_points[1].z),
               b,
               (base_points[0].y - base_points[2].y),
               (base_points[0].z - base_points[2].z),
               c,
               (base_points[0].y - base_points[3].y),
               (base_points[0].z - base_points[3].z)};
    return simple_3_by_3_matrix_determinant(delta_z);
}

// The a, b and c values are different only by which point is used in calculation which makes it easy to
// write one function to calculate any of them depending on select value
// select is simply the index of the point used
float transform_calculator::calculate_abc (const size_t select) const {
    if (select == 0 || select > 3) {
        return -1;
    }
    return alpha * 2 * (base_points[0].x * base_points[0].x - base_points[select].x * base_points[select].x) +
           beta * 2 * (base_points[0].y * base_points[0].y - base_points[select].y * base_points[select].y) +
           gamma * 2 * (base_points[0].z * base_points[0].z - base_points[select].z * base_points[select].z);
}

//function that does everything basically. After it's run the x_0, ..., alpha and beta variables should have their correct values
void transform_calculator::calculate_coeffs() {
    float delta = calculate_det_delta();
    float A_1 = calculate_det_A_1();
    float C_1 = calculate_det_C_1();
    float C_2 = calculate_det_C_2();
    float B_2 = calculate_det_B_2();
    float B_1 = calculate_det_B_1();
    float A_2 = calculate_det_A_2();
    float A = C_1 - (base_points[0].z * base_points[0].z - base_points[1].z * base_points[1].z);
    float B = C_2 - (base_points[0].z * base_points[0].z - base_points[2].z * base_points[2].z);
    float D_A_matrix[4] {A, (base_points[0].y * base_points[0].y - base_points[1].y * base_points[1].y) - B_1,
           B, (base_points[0].y * base_points[0].y - base_points[2].y * base_points[2].y) - B_2};
    float D_B_matrix[4] {(base_points[0].x * base_points[0].x - base_points[1].x * base_points[1].x) - A_1, A,
                  (base_points[0].x * base_points[0].x - base_points[2].x * base_points[2].x) - A_2, B};
    float D_matrix[4] {(base_points[0].x * base_points[0].x - base_points[1].x * base_points[1].x) - A_1,
                (base_points[0].y * base_points[0].y - base_points[1].y * base_points[1].y) - B_1,
                (base_points[0].x * base_points[0].x - base_points[2].x * base_points[2].x) - A_2,
                (base_points[0].y * base_points[0].y - base_points[2].y * base_points[2].y) - B_2};
    alpha = simple_2_by_2_matrix_determinant(D_A_matrix) / simple_2_by_2_matrix_determinant(D_matrix);
    beta = simple_2_by_2_matrix_determinant(D_B_matrix) / simple_2_by_2_matrix_determinant(D_matrix);
    gamma = 1.0;
    float a = calculate_abc(1);
    float b = calculate_abc(2);
    float c = calculate_abc(3);

    float delta_z = calculate_delta_z(a, b, c);
    float delta_x = alpha * A_1 + beta * B_1 + gamma * C_1;
    float delta_y = alpha * A_2 + beta * B_2 + gamma * C_2;
    x_0 = delta_x / (2 * alpha * delta);
    y_0 = delta_y / (2 * beta * delta);
    z_0 = delta_z / (2 * gamma * delta);

}

