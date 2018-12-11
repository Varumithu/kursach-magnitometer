//
// Created by valmit on 23/11/18.
//

#include "magnitometer.h"
#include "determinant.h"

#include <cmath>
#include <algorithm>


transform_calculator::transform_calculator(const double array[12]) { // I do not check length here, if it is too short vector will throw somwhere
    for (size_t i = 0; i < 4; ++i) {
        base_points[i] = point(array[i * 3], array[i * 3 + 1], array[i * 3 + 2]);
    }
    calculate_coeffs();
}

//transform_calculator::transform_calculator(std::vector<point>& base_points) : base_points(base_points) {
//    calculate_coeffs();
//}


double transform_calculator::calculate_det_A_1 () const {
    double A_1[9]{base_points[0].x * base_points[0].x - base_points[1].x * base_points[1].x,
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
double transform_calculator::calculate_det_B_1 () const {
    double B_1[9] {base_points[1].y * (base_points[0].y - base_points[1].y),
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
double transform_calculator::calculate_det_C_1 () const {
    double C_1[9] =  {base_points[1].z * (base_points[0].z - base_points[1].z),
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
double transform_calculator::calculate_det_A_2 () const {
    double A_2[9] {base_points[0].x - base_points[1].x,
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
double transform_calculator::calculate_det_B_2 () const {
    double B_2[9] {base_points[0].x - base_points[1].x,
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
double transform_calculator::calculate_det_C_2 () const {
    double C_2[9] {base_points[0].x - base_points[1].x,
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

double transform_calculator::calculate_det_delta() const {
    double m_delta[9] {(base_points[0].x - base_points[1].x),
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

double transform_calculator::calculate_delta_z (const double a, const double b, const double c) const {
    double delta_z[9] {a,
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

double transform_calculator::calculate_abc (const size_t select) const {
    if (select == 0 || select > 3) {
        return -1;
    }
    return alpha * 2 * (base_points[0].x * base_points[0].x - base_points[select].x * base_points[select].x) +
           beta * 2 * (base_points[0].y * base_points[0].y - base_points[select].y * base_points[select].y) +
           gamma * 2 * (base_points[0].z * base_points[0].z - base_points[select].z * base_points[select].z);
}

void transform_calculator::calculate_coeffs() {
    double delta = calculate_det_delta();
    double A_1 = calculate_det_A_1();
    double C_1 = calculate_det_C_1();
    double C_2 = calculate_det_C_2();
    double B_2 = calculate_det_B_2();
    double B_1 = calculate_det_B_1();
    double A_2 = calculate_det_A_2();
    double A = C_1 - (base_points[0].z * base_points[0].z - base_points[1].z * base_points[1].z);
    double B = C_2 - (base_points[0].z * base_points[0].z - base_points[2].z * base_points[2].z);
    double D_A_matrix[4] {A, (base_points[0].y * base_points[0].y - base_points[1].y * base_points[1].y) - B_1,
           B, (base_points[0].y * base_points[0].y - base_points[2].y * base_points[2].y) - B_2};
    double D_B_matrix[4] {(base_points[0].x * base_points[0].x - base_points[1].x * base_points[1].x) - A_1, A,
                  (base_points[0].x * base_points[0].x - base_points[2].x * base_points[2].x) - A_2, B};
    double D_matrix[4] {(base_points[0].x * base_points[0].x - base_points[1].x * base_points[1].x) - A_1,
                (base_points[0].y * base_points[0].y - base_points[1].y * base_points[1].y) - B_1,
                (base_points[0].x * base_points[0].x - base_points[2].x * base_points[2].x) - A_2,
                (base_points[0].y * base_points[0].y - base_points[2].y * base_points[2].y) - B_2};
    alpha = simple_2_by_2_matrix_determinant(D_A_matrix) / simple_2_by_2_matrix_determinant(D_matrix);
    beta = simple_2_by_2_matrix_determinant(D_B_matrix) / simple_2_by_2_matrix_determinant(D_matrix);
    gamma = 1.0;
    double a = calculate_abc(1);
    double b = calculate_abc(2);
    double c = calculate_abc(3);

    double delta_z = calculate_delta_z(a, b, c);
    double delta_x = alpha * A_1 + beta * B_1 + gamma * C_1;
    double delta_y = alpha * A_2 + beta * B_2 + gamma * C_2;
    x_0 = delta_x / (2 * alpha * delta);
    y_0 = delta_y / (2 * beta * delta);
    z_0 = delta_z / (2 * gamma * delta);

}

