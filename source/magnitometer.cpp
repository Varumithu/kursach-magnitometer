//
// Created by valmit on 23/11/18.
//

#include "magnitometer.h"

#include <cmath>
#include <algorithm>


double calculate_squared_range(const point& first, const point& second) {//
    return (first.x - second.x) * (first.x - second.x) + (first.y - second.y) * (first.y - second.y) + (first.z - second.z) * (first.z - second.z);
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
    double x_offset = std::max_element(data.begin(), data.end(), [](point first, point second){ return (first.x < second.x); })->x;
    double y_offset = std::max_element(data.begin(), data.end(), [](point first, point second){ return (first.y < second.y); })->y;
    double z_offset = std::max_element(data.begin(), data.end(), [](point first, point second){ return (first.z < second.z); })->z;
    return {x_offset, y_offset, z_offset};
}

void magnitometer_data::apply_offset() {
    point offset = calculate_offset();
    for (size_t i = 0; i < data.size(); ++i) {
        data[i].x = data[i].x - offset.x;
        data[i].y = data[i].y - offset.y;
        data[i].z = data[i].z - offset.z;
    }
}