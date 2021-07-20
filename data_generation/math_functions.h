#ifndef __MATH_F_H_INCLUDED__
#define __MATH_F_H_INCLUDED__


#include <math.h>
#include <limits>

#include "math_functions.h"

bool const INF_CHECK = false;
double const INF_WELL_RADIUS = 15;
double const OMEGA = 0.00001;
double const H = 6.62607015*pow(10,-34);
double const H_BAR = 1.0545*pow(10,-34);
double const H_BAR_SQUARED = pow(H_BAR,2);
double const PI = 3.14159265358979323846;
double const INF = std::numeric_limits<double>::infinity();
double const Z = 1.0;
double const E = 1.60217*pow(10,-19);//1.0;
double const ME =  9.109*pow(10,-31);//1.0
double const PERM = 8.854*pow(10,-12);
double const A0 = 4*PI*PERM*H_BAR_SQUARED/(ME*E*E);
double const B1 = 4*Z*ME*E*E/H_BAR_SQUARED;
double const B2 = B1/3;
double const LAMBDA = 0.01;
double const MASS_OVER_HBAR = (9.109/1.0545)*pow(10,3);

struct complex {
  double re;
  double im;
};


complex hermite_contour(double x, double y, double n);
complex f(double t, double x, double n);
complex f2(double t, double x, double n);

double j_d(double x, double y, double n);

double get_theta(double x, double y);


double factorial_double(double n);



/**
    using a 4th degree runge kutta method to do the trajectory evolution.

    @param (*f) the function used in the method
    @param timestep the time different used in the runge kutta method
*/
double runge_kutta(double (*f)(int*, double, int, int), double time_step);



/**
    fitting a polunomial to x and y

    @param x the axis locations
    @param y the values to be fit
    @param coeffs a pointer to the coefficients which start off 0 and get filled at the end of the function
    @param degree the degree of the polynomial we're fitting
*/
void fit_polynomial(double *x, double *y, double *coeffs, int degree);



/**
    taking a nth degree derivative of a polynomial

    @param coeffs a pointer holding the coefficients of the polynomial we're taking the derivative of
    @param x_val the point to be evaluated
    @param degree_of_fit the number of coefficients - 1
    @param n the degree of the derivative
*/
double nth_derivative_polynomial(double *coeffs, double x_val, int degree_of_fit, int n);



/**
    converting the coordinates to an indexing

    @param coordiantes the coordinates of the point
    @param config_dimension the spatial_dimension*number of particles, the dimension of our grid
    @param grid_length the length of the grid in each dimension
*/
int coordinates_to_index(int* coordiantes, int config_dimension, int grid_length);



/**
    converting the index to coordinates

    @param coordiantes the coordinates of the point
    @param index the index that will be converted
    @param config_dimension the spatial_dimension*number of particles, the dimension of our grid
    @param grid_length the length of the grid in each dimension
*/
void index_to_coordinates(int* coordinates, int index, int config_dimension, int grid_length);



/**
    calculating the distance between two points, used by the potential_function

    @param position_1
    @param position_2
    @param spatial_dimension
*/
double calc_distance(double *position_1, double *position_2, int spatial_dimension);



/**
    a power function for integers, x^p, because normal pow had funny behavior when converting ints to doubles

    @param x base
    @param p power
*/
int int_pow(int x, int p);



/**
    printing progress during point initialization

    @param progress the percentage of points that have been initialized
*/
void print_load_bar(double progress);



/**
    a description you could pick up from the title

    @param array yep
    @param size definitely not the size of the array
*/
void print_double_array(double *array, int size);


/**
    the first kind bessel function for n = 1, up to a certain term in the sum

    @param x the point to evaluate at
*/
double j_1(double x, int grid_length);

double j_n(double x, int n);
double j_n_non_int(double x, double n);
int factorial(int n);
double non_int_hermite(double x, double n, int mult);
#endif
