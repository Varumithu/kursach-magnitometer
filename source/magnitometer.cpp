//
// Created by valmit on 23/11/18.
//

#include "magnitometer.h"

#include <cmath>
#include <algorithm>


double calculate_squared_range(const point& first, const point& second) {//
    return (first.x - second.x) * (first.x - second.x) + (first.y - second.y) * (first.y - second.y) + (first.z - second.z) * (first.z - second.z);
}

std::vector<size_t> magnitometer_data::find_base_points() {
    std::vector<size_t> res{0, 0, 0, 0};
    //TODO

}

double magnitometer_data::calculate_det_A_1 (const std::vector<point>& base_points) const {
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
double magnitometer_data::calculate_det_B_1 (const std::vector<point>& base_points) const {
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
double magnitometer_data::calculate_det_C_1 (const std::vector<point>& base_points) const {
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
double magnitometer_data::calculate_det_A_2 (const std::vector<point>& base_points) const {
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
double magnitometer_data::calculate_det_B_2 (const std::vector<point>& base_points) const {
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
double magnitometer_data::calculate_det_C_2 (const std::vector<point>& base_points) const {
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

double magnitometer_data::calculate_det_delta(const std::vector<point>& base_points) const {
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

void magnitometer_data::calculate_coeffs(std::vector<point>& base_points) {
    double delta = calculate_det_delta(base_points);
    double A_1 = calculate_det_A_1(base_points);
    double C_1 = calculate_det_C_1(base_points);
    double C_2 = calculate_det_C_2(base_points);
    double B_2 = calculate_det_B_2(base_points);
    double B_1 = calculate_det_B_1(base_points);
    double A_2 = calculate_det_A_2(base_points);
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
    double alpha = D_A_matrix.determinant() / D_matrix.determinant();
    double beta = D_B_matrix.determinant() / D_matrix.determinant();

}


size_t magnitometer_data::find_nearest(const point &&that) const { // the first search implementation that comes to mind, maybe should change it later
    double range = calculate_squared_range(data[0], that);
    size_t res = 0;
    for (size_t i = 1; i < data.size(); ++i) {
        double cur_range = calculate_squared_range(data[i], that);
        if (cur_range < range) {
            range = cur_range;
            res = i;
        }
    }
    return  res;
}

point magnitometer_data::calculate_offset() const {
    double x_max = std::max_element(data.begin(), data.end(), [](point first, point second){ return (first.x < second.x); })->x;
    double y_max = std::max_element(data.begin(), data.end(), [](point first, point second){ return (first.y < second.y); })->y;
    double z_max = std::max_element(data.begin(), data.end(), [](point first, point second){ return (first.z < second.z); })->z;
    double x_min = std::min_element(data.begin(), data.end(), [](point first, point second){ return (first.x < second.x); })->x;
    double y_min = std::min_element(data.begin(), data.end(), [](point first, point second){ return (first.y < second.y); })->y;
    double z_min = std::min_element(data.begin(), data.end(), [](point first, point second){ return (first.z < second.z); })->z;
    return {x_max / 2 + x_min / 2, y_max / 2 + y_min / 2, z_max / 2 + z_min / 2};
}

void magnitometer_data::apply_offset() {
    point offset = calculate_offset();
    for (size_t i = 0; i < data.size(); ++i) {
        data[i].x = data[i].x - offset.x;
        data[i].y = data[i].y - offset.y;
        data[i].z = data[i].z - offset.z;
    }
}
