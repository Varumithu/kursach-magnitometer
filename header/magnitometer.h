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



class transform_calculator final {
private:
    double calculate_abc (const size_t select) const;
    double calculate_det_delta () const;
    double calculate_det_A_1 () const;
    double calculate_det_B_1 () const;
    double calculate_det_C_1 () const;
    double calculate_det_A_2 () const;
    double calculate_det_B_2 () const;
    double calculate_det_C_2 () const;
    double calculate_delta_z (const double a, const double b, const double c) const;
public:
    double x_0, y_0, z_0, alpha, beta, gamma = 1.0;
    transform_calculator(std::vector<point>& base_points);

    std::vector<point> base_points;

    void calculate_coeffs();
};











#endif //KURSACH_MAGNITOMETER_MAGNITOMETER_H
