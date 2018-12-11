#include <iostream>
#include <vector>

#include "magnitometer.h"

int main() {
    float points[12]{1, 2, 3, 3, 4, 5, -2, 4, -4, 3, 2, 5};
    transform_calculator test1(points);
    std::cout << test1.alpha << ' ' << test1.beta << ' ' << test1.gamma << ' ' << test1.x_0 << ' ' << test1.y_0 << ' ' << test1.z_0 << std::endl;
}
