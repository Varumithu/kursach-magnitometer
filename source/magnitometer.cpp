//
// Created by valmit on 23/11/18.
//

#include "magnitometer.h"
#include <cmath>

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