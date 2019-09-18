#ifndef __MATH_F_H_INCLUDED__
#define __MATH_F_H_INCLUDED__


#include <math.h>
#include "math_functions.h"



double const H_BAR = 6.62607004*pow(10,-34);
double const EPSILON = 2*pow(10,-2);




double calc_distance(double *position_1, double *position_2);


double calc_velocity_i(int *coordinates, double time, int i);


double calc_rho();


double calc_potential();


double calc_quantum_potential();


double gradient_i(double (*f)(double), double x, int i);


double divergence(double (*f)(double), double x);


double laplacian(double (*f)(double), double x);


double derivative_i(double (*f)(double), double x, int i);

#endif
