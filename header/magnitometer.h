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
    float x, y, z;
    point(float x, float y, float z) : x(x), y(y), z(z) {}
    point() = default;
    point(const point& other) noexcept : x(other.x), y(other.y), z(other.z) {}

};


/* main class that does the calculation of desired coefficients according to the algorythm
 * The constructor accepts an array of 12 floats being 4 sets of three coordinates aka 4 points,
 * constructs base_points array and calls the calculate_coeffs methods which sets the x_0, ..., alpha and beta data members
 * They are public because i did not want to write getters, but they should not be changed
*/
class transform_calculator final {
private:
    float calculate_abc (const size_t select) const;
    float calculate_det_delta () const;
    float calculate_det_A_1 () const;
    float calculate_det_B_1 () const;
    float calculate_det_C_1 () const;
    float calculate_det_A_2 () const;
    float calculate_det_B_2 () const;
    float calculate_det_C_2 () const;
    float calculate_delta_z (const float a, const float b, const float c) const;

    void calculate_coeffs();
public:
    float x_0, y_0, z_0, alpha, beta, gamma = 1.0;
    transform_calculator(const float array[12]);
    point base_points[4];


};











#endif //KURSACH_MAGNITOMETER_MAGNITOMETER_H
