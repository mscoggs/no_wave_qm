#include <math.h>
#include <fstream>


#include "input_functions.h"
#include "math_functions.h"






Psi calc_psi(double x, double y, int psi_function, double coord_to_distance, double *mass){
  x = x*coord_to_distance, y = y*coord_to_distance;

  double re, im, r = sqrt(x*x+y*y), theta= get_theta(x,y);
  double constant, k;
  switch(psi_function){
    case 12: //perturbed inf-well
      re = j_d(x,y,1)*cos(theta);
      im = j_d(x,y,1)*sin(theta);
      break;
    case 13: //perturbed inf-well
      re = j_d(x,y,1)*cos(theta) + LAMBDA*j_d(x,y,1.5)*cos(1.5*theta);
      im = j_d(x,y,1)*sin(theta) + LAMBDA*j_d(x,y,1.5)*sin(1.5*theta);
      break;
    case 15: //2d ho 10 01
      k = mass[0]*OMEGA/H_BAR;
      constant = (k/sqrt(2*PI))*exp(-k*(r*r)/2);
      re = (constant/sqrt(pow(2,1)*factorial(1)*factorial(0)))*x;
      im = (constant/sqrt(pow(2,1)*factorial(0)*factorial(1)))*y;
      break;
    case 16: //2d ho 30 01
      k = mass[0]*OMEGA/H_BAR;
      constant = sqrt(k/PI)*exp(-k*(r*r)/2);
      x = sqrt(k)*x;
      y = sqrt(k)*y;
      re = (constant/sqrt(pow(2,3)*factorial(3)*factorial(0)))*(pow(x,3) - 3*x);
      im = (constant/sqrt(pow(2,1)*factorial(0)*factorial(1)))*y;
      break;
  }

  Psi psi = {re, im};
  return psi;
}


double calc_potential(double time, int *coordinates, int num_particles, int spatial_dimension, int potential_function, int grid_length, double *mass){
  double x, y, center, r,  v, x_0,y_0, r_x, r_y, width;

  center = (grid_length-1.0)/2.0;
  x = static_cast<double>(coordinates[0])-center;
  y = static_cast<double>(coordinates[1])-center;
  r = sqrt(x*x + y*y);

  switch(potential_function){
    case 4: //2d_harmonic_oscillator stationary
      return 0.5*mass[0]*OMEGA*OMEGA*(x*x + y*y);
    case 5: //inf_circular_well
      if(sqrt(x*x + y*y) < INF_WELL_RADIUS) return 0;
      else return INF;
    case 6: //2D Hydrogen
      return -(Z*E*E)/r;
    case 9: //wave 1
      v = static_cast<double>(grid_length)/20000;
      x_0 = -grid_length/2.0;
      r_x = x-v*time - x_0;
      width = grid_length/10.0;
      return 50*mass[0]*OMEGA*OMEGA*exp(-r_x*r_x/(2*width*width)) + 0.5*mass[0]*OMEGA*OMEGA*(x*x + y*y);
    case 10: // passing off-center node
      v = static_cast<double>(grid_length)/20000;
      x_0 = -grid_length/2.0;
      y_0 = -grid_length/4.0;
      r_x = x-v*time - x_0;
      r_y = y - y_0;
      width = grid_length/10.0;
      return 50*mass[0]*OMEGA*OMEGA*exp(-(r_x*r_x + r_y*r_y)/(2*width*width)) + 0.5*mass[0]*OMEGA*OMEGA*(x*x + y*y);
    case 11: // passing off-center node which vanishes halfway
      v = static_cast<double>(grid_length)/20000;
      if(time > 10000) return 0.5*mass[0]*OMEGA*OMEGA*(x*x + y*y);
      x_0 = -grid_length/2.0;
      y_0 = -grid_length/4.0;
      r_x = x-v*time - x_0;
      r_y = y - y_0;
      width = grid_length/10.0;
      return 50*mass[0]*OMEGA*OMEGA*exp(-(r_x*r_x + r_y*r_y)/(2*width*width)) + 0.5*mass[0]*OMEGA*OMEGA*(x*x + y*y);
  }
}

double calc_rho_vel_initial(int *coordinates,double *velocities, int config_dimension, int grid_length, double *mass, int psi_function, double coord_to_distance, int velocity_perturbation, int density_perturbation){
  int a = psi_function;
  double x, y, center = (grid_length-1.0)/2.0, A,epsilon = EPSILON_INIT_CALCS, d0,d1, rho_pert = 0, rho_scale=1, *x3,*y3, *abc;


  x3   = new double[3]();
  y3   = new double[3]();
  abc = new double[3]();
  x = static_cast<double>(coordinates[0])-center;
  y = static_cast<double>(coordinates[1])-center;
  Psi psi = calc_psi(x,y,a,coord_to_distance,mass), psi0, psi2;


  if(abs(psi.real) == 0){
    velocities[0] =  0;
    velocities[1] =  0;
  }else{
    x3[0] = 2-epsilon, x3[1] = 2, x3[2] = 2+epsilon;
    psi0 = calc_psi(x-epsilon,y,a,coord_to_distance,mass), psi2 = calc_psi(x+epsilon,y,a,coord_to_distance,mass);
    y3[0] = psi0.imaginary/psi0.real, y3[1] = psi.imaginary/psi.real, y3[2] = psi2.imaginary/psi2.real;
    fit_polynomial(x3, y3, abc, 2);
    d0 = nth_derivative_polynomial(abc, 2, 2, 1);

    x3[0] = 2-epsilon, x3[1] = 2, x3[2] = 2+epsilon;
    psi0 = calc_psi(x,y-epsilon,a,coord_to_distance,mass), psi2 = calc_psi(x,y+epsilon,a,coord_to_distance,mass);
    y3[0] = psi0.imaginary/psi0.real, y3[1] = psi.imaginary/psi.real, y3[2] = psi2.imaginary/psi2.real;
    fit_polynomial(x3, y3, abc, 2);
    d1 = nth_derivative_polynomial(abc, 2, 2, 1);

    A = H_BAR/(mass[0]*(1 + pow(psi.imaginary/psi.real,2)));
    velocities[0] =  A*d0;
    velocities[1] =  A*d1;
  }

  switch(velocity_perturbation){
    case 0:
      break;
    case 1:
      velocities[0] += -0.01*H_BAR*y/(mass[0]*(x*x+y*y));
      velocities[1] += 0.01*H_BAR*x/(mass[0]*(x*x+y*y));
      break;
  }

  switch(density_perturbation){
    case 0:
      break;
    case 1:
      rho_pert = 0;
      break;
    case 2:
      rho_scale = get_random_double(0.8, 1.2);
      break;
  }

  delete[] x3, delete[] y3, delete[] abc;
  return rho_scale*pow(psi.real,2)+pow(psi.imaginary,2) + rho_pert;
}
