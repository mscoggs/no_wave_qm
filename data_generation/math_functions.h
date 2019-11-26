#ifndef __MATH_F_H_INCLUDED__
#define __MATH_F_H_INCLUDED__


#include <math.h>
#include "math_functions.h"






int coordinates_to_index(int* coordiantes, int config_dimension, int grid_length);
void index_to_coordinates(int* coordinates, int index, int config_dimension, int grid_length);
void fit_polynomial(double *x, double *y, double *coeffs, int degree);
double derivative_polynomial(double *coeffs, double x_val);
double second_derivative_polynomial(double *coeffs, double x_val);
int int_pow(int x, int p);

double calc_rho_initial(int *coordinates, double time, int config_dimension, int grid_length);
void calc_velocities_initial(int *coordinates, double time, double *velocities, int config_dimension);
double calc_distance(double *position_1, double *position_2, int spatial_dimension);
double calc_potential(int *coordinates, int num_particles, int spatial_dimension);
double runge_kutta(double (*f)(int*, double, int, int), double time_step);

// double divergence(double (*f)(int*, double, int, int), int * coordinates, double time, int config_dimension);
// double derivative_i_potential(double (*f)(int*, int, int), int *coordinates, int config_dimension, int num_particles, int spatial_dimension, double epsilon, int index);
// double derivative_i_quantum_potential(double (*f)(double*, int*, double, int), double* mass, int *coordinates, int config_dimension, double time , double epsilon, int index);
// double derivative_i_velocity(double (*f)(int*, double, int, int), int *coordinates, double time, int config_dimension, double epsilon, int index, int j);
// double derivative_i(double (*f)(int*, double, int, int), int *coordinates, double time, int config_dimension, double epsilon, int index);



double const H_BAR = 6.62607004*pow(10,-34);
double const H_BAR_SQUARED = pow(H_BAR,2);
double const EPSILON = 2*pow(10,-2);
int    const DEGREE_OF_FIT = 2;
double const POINT_DISTANCE = 10;
double const OFFSET = 2*POINT_DISTANCE;

#endif
