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

void magnitometer_data::calculate_coeffs(std::vector<size_t>& base_points) {
    Eigen::Matrix<double, 3, 3> m_delta;
    m_delta << (data[base_points[0]].x - data[base_points[1]].x),
               (data[base_points[0]].y - data[base_points[1]].y),
               (data[base_points[0]].z - data[base_points[1]].z),
               (data[base_points[0]].x - data[base_points[2]].x),
               (data[base_points[0]].y - data[base_points[2]].y),
               (data[base_points[0]].z - data[base_points[2]].z),
               (data[base_points[0]].x - data[base_points[3]].x),
               (data[base_points[0]].y - data[base_points[3]].y),
               (data[base_points[0]].z - data[base_points[3]].z);
    double delta = m_delta.determinant();
    Eigen::Matrix<double, 3, 3> A_1, A_2, B_1, B_2, C_1, C_2;
    A_1 << data[base_points[0]].x * data[base_points[0]].x - data[base_points[1]].x * data[base_points[1]].x,
           data[base_points[0]].y - data[base_points[1]].y,
           data[base_points[0]].z - data[base_points[1]].z,
           data[base_points[0]].x * data[base_points[0]].x - data[base_points[2]].x * data[base_points[2]].x,
           data[base_points[0]].y - data[base_points[2]].y,
           data[base_points[0]].z - data[base_points[2]].z,
           data[base_points[0]].x * data[base_points[0]].x - data[base_points[3]].x * data[base_points[3]].x,
           data[base_points[0]].y - data[base_points[3]].y,
           data[base_points[0]].z - data[base_points[3]].z;
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