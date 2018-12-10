#include <iostream>
#include <vector>

#include "magnitometer.h"

int main() {
    point p1 (1, 2, 3), p2(3, 4, 5), p3(-2, 4, -4), p4(3, 2, 5);
    std::vector<point> base_points;
    base_points.push_back(p1);
    base_points.push_back(p2);
    base_points.push_back(p3);
    base_points.push_back(p4);
    transform_calculator test(base_points);
    std::cout << test.alpha << ' ' << test.beta << ' ' << test.gamma << ' ' << test.x_0 << ' ' << test.y_0 << ' ' << test.z_0;
}
