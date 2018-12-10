//
// Created by valmit on 23/11/18.
//

#include "magnitometer.h"

#include <cmath>
#include <algorithm>

transform_calculator::transform_calculator(std::vector<double>& array) { // I do not check length here, if it is too short vector will throw somwhere
    for (size_t i = 0; i < 4; ++i) {
        base_points.emplace_back(array[i * 3], array[i * 3 + 1], array[i * 3 + 2]);
    }
    calculate_coeffs();
}

transform_calculator::transform_calculator(std::vector<point>& base_points) : base_points(base_points) {
    calculate_coeffs();
}


double transform_calculator::calculate_det_A_1 () const {
    Eigen::Matrix<double, 3, 3> A_1;
    A_1 << base_points[0].x * base_points[0].x - base_points[1].x * base_points[1].x,
           base_points[0].y - base_points[1].y,
           base_points[0].z - base_points[1].z,
           base_points[0].x * base_points[0].x - base_points[2].x * base_points[2].x,
           base_points[0].y - base_points[2].y,
           base_points[0].z - base_points[2].z,
           base_points[0].x * base_points[0].x - base_points[3].x * base_points[3].x,
           base_points[0].y - base_points[3].y,
           base_points[0].z - base_points[3].z;
    return A_1.determinant();
}
double transform_calculator::calculate_det_B_1 () const {
    Eigen::Matrix<double, 3, 3> B_1;
    B_1 << base_points[1].y * (base_points[0].y - base_points[1].y),
           base_points[0].y - base_points[1].y,
           base_points[0].z - base_points[1].z,
           base_points[2].y * (base_points[0].y - base_points[2].y),
           base_points[0].y - base_points[2].y,
           base_points[0].z - base_points[2].z,
           base_points[3].y * (base_points[0].y - base_points[3].y),
           base_points[0].y - base_points[3].y,
           base_points[0].z - base_points[3].z;
    return B_1.determinant();
}
double transform_calculator::calculate_det_C_1 () const {
    Eigen::Matrix<double, 3, 3> C_1;
    C_1 << base_points[1].z * (base_points[0].z - base_points[1].z),
           base_points[0].y - base_points[1].y,
           base_points[0].z - base_points[1].z,
           base_points[2].z * (base_points[0].z - base_points[2].z),
           base_points[0].y - base_points[2].y,
           base_points[0].z - base_points[2].z,
           base_points[3].z * (base_points[0].z - base_points[3].z),
           base_points[0].y - base_points[3].y,
           base_points[0].z - base_points[3].z;
    return C_1.determinant();
}
double transform_calculator::calculate_det_A_2 () const {
    Eigen::Matrix<double, 3, 3> A_2;
    A_2 << base_points[0].x - base_points[1].x,
           base_points[1].x * (base_points[0].x - base_points[1].x),
           base_points[0].z - base_points[1].z,
           base_points[0].x - base_points[2].x,
           base_points[2].x * (base_points[0].x - base_points[2].x),
           base_points[0].z - base_points[2].z,
           base_points[0].x - base_points[3].x,
           base_points[3].x * (base_points[0].x - base_points[3].x),
           base_points[0].z - base_points[3].z;
    return A_2.determinant();
}
double transform_calculator::calculate_det_B_2 () const {
    Eigen::Matrix<double, 3, 3> B_2;
    B_2 << base_points[0].x - base_points[1].x,
           base_points[0].y * base_points[0].y - base_points[1].y * base_points[1].y,
           base_points[0].z - base_points[1].z,
           base_points[0].x - base_points[2].x,
           base_points[0].y * base_points[0].y - base_points[2].y * base_points[2].y,
           base_points[0].z - base_points[2].z,
           base_points[0].x - base_points[3].x,
           base_points[0].y * base_points[0].y - base_points[3].y * base_points[3].y,
           base_points[0].z - base_points[3].z;
    return B_2.determinant();
}
double transform_calculator::calculate_det_C_2 () const {
    Eigen::Matrix<double, 3, 3> C_2;
    C_2 << base_points[0].x - base_points[1].x,
           base_points[1].z * (base_points[0].z - base_points[1].z),
           base_points[0].z - base_points[1].z,
           base_points[0].x - base_points[2].x,
           base_points[2].x * (base_points[0].x - base_points[2].x),
           base_points[0].z - base_points[2].z,
           base_points[0].x - base_points[3].x,
           base_points[3].x * (base_points[0].x - base_points[3].x),
           base_points[0].z - base_points[3].z;
    return C_2.determinant();
}

double transform_calculator::calculate_det_delta() const {
    Eigen::Matrix<double, 3, 3> m_delta;
    m_delta << (base_points[0].x - base_points[1].x),
               (base_points[0].y - base_points[1].y),
               (base_points[0].z - base_points[1].z),
               (base_points[0].x - base_points[2].x),
               (base_points[0].y - base_points[2].y),
               (base_points[0].z - base_points[2].z),
               (base_points[0].x - base_points[3].x),
               (base_points[0].y - base_points[3].y),
               (base_points[0].z - base_points[3].z);
    return m_delta.determinant();
}

double transform_calculator::calculate_delta_z (const double a, const double b, const double c) const {
    Eigen::Matrix<double, 3, 3> delta_z_m;
    delta_z_m << a,
               (base_points[0].y - base_points[1].y),
               (base_points[0].z - base_points[1].z),
               b,
               (base_points[0].y - base_points[2].y),
               (base_points[0].z - base_points[2].z),
               c,
               (base_points[0].y - base_points[3].y),
               (base_points[0].z - base_points[3].z);
    return delta_z_m.determinant();
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
    Eigen::Matrix<double, 2, 2> D_A_matrix;
    D_A_matrix << A, (base_points[0].y * base_points[0].y - base_points[1].y * base_points[1].y) - B_1,
           B, (base_points[0].y * base_points[0].y - base_points[2].y * base_points[2].y) - B_2;
    Eigen::Matrix<double, 2, 2> D_B_matrix;
    D_B_matrix << (base_points[0].x * base_points[0].x - base_points[1].x * base_points[1].x) - A_1, A,
                  (base_points[0].x * base_points[0].x - base_points[2].x * base_points[2].x) - A_2, B;
    Eigen::Matrix<double, 2, 2> D_matrix;
    D_matrix << (base_points[0].x * base_points[0].x - base_points[1].x * base_points[1].x) - A_1,
                (base_points[0].y * base_points[0].y - base_points[1].y * base_points[1].y) - B_1,
                (base_points[0].x * base_points[0].x - base_points[2].x * base_points[2].x) - A_2,
                (base_points[0].y * base_points[0].y - base_points[2].y * base_points[2].y) - B_2;
    alpha = D_A_matrix.determinant() / D_matrix.determinant();
    beta = D_B_matrix.determinant() / D_matrix.determinant();
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

