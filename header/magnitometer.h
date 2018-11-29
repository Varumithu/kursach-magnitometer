//
// Created by valmit on 23/11/18.
//

#ifndef KURSACH_MAGNITOMETER_MAGNITOMETER_H
#define KURSACH_MAGNITOMETER_MAGNITOMETER_H

#include <vector>
#include <Eigen/Dense>

class point final {
public:
    double x, y, z;
    point(double x, double y, double z) : x(x), y(y), z(z) {}
    point() = default;
    point(const point& other) noexcept : x(other.x), y(other.y), z(other.z) {}

};



class magnitometer_data final {
public:

    std::vector<point> data;

    size_t find_nearest(const point&& that) const; // returns position of the nearest element in the data vector

    point calculate_offset() const; // returns a point just for convenience, three numbers that must be substracted from coordintaes
    void apply_offset();
    std::vector<size_t> find_base_points();
    void calculate_coeffs(std::vector<size_t>& base_points);
};

double calculate_squared_range(point& first, point& second);












#endif //KURSACH_MAGNITOMETER_MAGNITOMETER_H
