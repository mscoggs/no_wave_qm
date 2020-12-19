#include <math.h>
#include <fstream>


#include "input_functions.h"
#include "math_functions.h"



Psi calc_psi(double x, double y, int psi_function, double coord_to_distance){
  x = x*coord_to_distance, y = y*coord_to_distance;

  double re, im, r = sqrt(x*x+y*y), theta= get_theta(x,y);

  switch(psi_function){
    case 12: //perturbed inf-well
      re = j_d(x,y,1)*cos(theta);
      im = j_d(x,y,1)*sin(theta);
    case 13: //perturbed inf-well
      re = j_d(x,y,1)*cos(theta) + LAMBDA*j_d(x,y,1.5)*cos(1.5*theta);
      im = j_d(x,y,1)*sin(theta) + LAMBDA*j_d(x,y,1.5)*sin(1.5*theta);
  }

  Psi psi = {re, im};
  return psi;
}



double calc_rho_vel_initial(int *coordinates,double *velocities, int config_dimension, int grid_length, double *mass, int psi_function, double coord_to_distance){
  int a = psi_function;
  double x, y, center, A,epsilon = pow(10,-3), d0,d1;
  double *x3,*y3, *abc;
  x3   = new double[3]();
  y3   = new double[3]();
  abc = new double[3]();

  center = (grid_length-1.0)/2.0;
  x = static_cast<double>(coordinates[0])-center;
  y = static_cast<double>(coordinates[1])-center;
  Psi psi = calc_psi(x,y,a,coord_to_distance), psi0, psi2;
  A = H_BAR/(mass[0]*(1 + pow(psi.imaginary/psi.real,2)));


  x3[0] = 2-epsilon, x3[1] = 2, x3[2] = 2+epsilon;
  psi0 = calc_psi(x-epsilon,y,a,coord_to_distance), psi2 = calc_psi(x+epsilon,y,a,coord_to_distance);
  y3[0] = psi0.imaginary/psi0.real, y3[1] = psi.imaginary/psi.real, y3[2] = psi2.imaginary/psi2.real;
  fit_polynomial(x3, y3, abc, 2);
  d0 = nth_derivative_polynomial(abc, 2, 2, 1);

  //x3[0] = y-epsilon, x3[1] = y, x3[2] = 2+epsilon;
  psi0 = calc_psi(x,y-epsilon,a,coord_to_distance), psi2 = calc_psi(x,y+epsilon,a,coord_to_distance);
  y3[0] = psi0.imaginary/psi0.real, y3[1] = psi.imaginary/psi.real, y3[2] = psi2.imaginary/psi2.real;
  fit_polynomial(x3, y3, abc, 2);
  d1 = nth_derivative_polynomial(abc, 2, 2, 1);


  velocities[0] =  A*d0;
  velocities[1] =  A*d1;


  return pow(psi.real,2)+pow(psi.imaginary,2);
}



double calc_potential(int *coordinates, int num_particles, int spatial_dimension, int v_function, int grid_length, double *mass){
  double x, y, center, r;
  double e = 0.5, V=0, *positioni, *positionj, distance =0;

  center = (grid_length-1.0)/2.0;
  x = static_cast<double>(coordinates[0])-center;
  y = static_cast<double>(coordinates[1])-center;

  switch(v_function){
    case 1: //plane_wave
      return 0;
    case 2: //holland_ex4.11_pg169
      return V;
    case 3:
      return 0;
    case 4: //2d_harmonic_oscillator
      return 0.5*mass[0]*OMEGA*OMEGA*(x*x + y*y);
    case 5: //inf_circular_well

      if(sqrt(x*x + y*y) < INF_WELL_RADIUS) return 0;
      else return INF;
    case 6: //2D Hydrogen
      r = sqrt(x*x + y*y);
      return -(Z*E*E)/r;
    case 9:
      return 0;
    case 10:
      return 0;
  }
}
