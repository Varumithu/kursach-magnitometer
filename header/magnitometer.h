//
// Created by valmit on 23/11/18.
//

#ifndef KURSACH_MAGNITOMETER_MAGNITOMETER_H
#define KURSACH_MAGNITOMETER_MAGNITOMETER_H


// for size_t defenition
#include <cstddef>


//helper class to make working with points more intuitive
class point final {
public:
    double x, y, z;
    point(double x, double y, double z) : x(x), y(y), z(z) {}
    point() = default;
    point(const point& other) noexcept : x(other.x), y(other.y), z(other.z) {}

};


/* main class that does the calculation of desired coefficients according to the algorythm
 * The constructor accepts an array of 12 doubles being 4 sets of three coordinates aka 4 points,
 * constructs base_points array and calls the calculate_coeffs methods which sets the x_0, ..., alpha and beta data members
 * They are public because i did not want to write getters, but they should not be changed
*/
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

    void calculate_coeffs();
public:
    double x_0, y_0, z_0, alpha, beta, gamma = 1.0;
    transform_calculator(const double array[12]);
    point base_points[4];


};











#endif //KURSACH_MAGNITOMETER_MAGNITOMETER_H
