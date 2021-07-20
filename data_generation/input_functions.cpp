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
    case 13: //perturbed inf-well
      re = j_d(x,y,1)*cos(theta) + LAMBDA*j_d(x,y,1.5)*cos(1.5*theta);
      im = j_d(x,y,1)*sin(theta) + LAMBDA*j_d(x,y,1.5)*sin(1.5*theta);
    case 15: //2d ho 10 01
      k = mass[0]*OMEGA/H_BAR;
      constant = (k/sqrt(2*PI))*exp(-k*(r*r)/2);
      //x = sqrt(k)*x;
      //y = sqrt(k)*y;
      re = (constant/sqrt(pow(2,1)*factorial(1)*factorial(0)))*x;
      im = (constant/sqrt(pow(2,1)*factorial(0)*factorial(1)))*y;
    case 16: //2d ho 30 01
      k = mass[0]*OMEGA/H_BAR;
      constant = sqrt(k/PI)*exp(-k*(r*r)/2);
      x = sqrt(k)*x;
      y = sqrt(k)*y;
      re = (constant/sqrt(pow(2,3)*factorial(3)*factorial(0)))*(pow(x,3) - 3*x);
      im = (constant/sqrt(pow(2,1)*factorial(0)*factorial(1)))*y;
  }
  //
  // printf("re: %f, x: %f\n",re, x );
  // printf("mass[0] %f\n", mass[0]);
  // printf("HBAR[0] %f\n", H_BAR);
  // printf("k %f\n", k);
  // printf("r: %f\n", r);
  // printf("exp(-k*(r*r)/2) %f\n", exp(-k*(r*r)/2));
  Psi psi = {re, im};
  return psi;
}


double calc_potential(double time, int *coordinates, int num_particles, int spatial_dimension, int potential_function, int grid_length, double *mass){
  double x, y, center, r;
  double e = 0.5, V=0, *positioni, *positionj, distance =0;

  center = (grid_length-1.0)/2.0;
  x = static_cast<double>(coordinates[0])-center;
  y = static_cast<double>(coordinates[1])-center;

  r = sqrt(x*x + y*y);
  double v = grid_length/200;
  double x_0 = grid_length/2.0;
  double z = x-v*time + x_0;
  double width = grid_length/10.0;

  switch(potential_function){
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
      return -(Z*E*E)/r;
    case 9:
      return 0.5*mass[0]*OMEGA*OMEGA*exp(-z*z/(2*width*width));
    case 10:
      return 0;
  }
}

double calc_rho_vel_initial(int *coordinates,double *velocities, int config_dimension, int grid_length, double *mass, int psi_function, double coord_to_distance, bool velocity_perturbation, double v_perturb){
  int a = psi_function;
  double x, y, center, A,epsilon = pow(10,-1), d0,d1;
  double *x3,*y3, *abc;
  x3   = new double[3]();
  y3   = new double[3]();
  abc = new double[3]();

  center = (grid_length-1.0)/2.0;
  x = static_cast<double>(coordinates[0])-center;
  y = static_cast<double>(coordinates[1])-center;
  Psi psi = calc_psi(x,y,a,coord_to_distance,mass), psi0, psi2;


  //THIS IS A TEMPORARY FIX -- y3 becomes nan in the case of psi.real being 0. In this case, the velocity at this point should be 0.
  //Hardset 0 and continue.
  if(abs(psi.real) == 0){
    printf("SKIPPING\n");
    velocities[0] =  0;
    velocities[1] =  0;
    return pow(psi.real,2)+pow(psi.imaginary,2);
  }


  A = H_BAR/(mass[0]*(1 + pow(psi.imaginary/psi.real,2)));


  x3[0] = 2-epsilon, x3[1] = 2, x3[2] = 2+epsilon;
  psi0 = calc_psi(x-epsilon,y,a,coord_to_distance,mass), psi2 = calc_psi(x+epsilon,y,a,coord_to_distance,mass);
  y3[0] = psi0.imaginary/psi0.real, y3[1] = psi.imaginary/psi.real, y3[2] = psi2.imaginary/psi2.real;
  fit_polynomial(x3, y3, abc, 2);
  d0 = nth_derivative_polynomial(abc, 2, 2, 1);
  if(d0 != d0){
    printf("coordinates %f, %f", x,y);
    printf("x %f %f %f", x3[0],x3[1], x3[2]);
    printf("y %f %f %f\n\n", y3[0],y3[1], y3[2]);
    printf("psi re, psi im: %f, %f", psi.real, psi.imaginary);
    exit(0);
  }

  x3[0] = 2-epsilon, x3[1] = 2, x3[2] = 2+epsilon;
  psi0 = calc_psi(x,y-epsilon,a,coord_to_distance,mass), psi2 = calc_psi(x,y+epsilon,a,coord_to_distance,mass);
  y3[0] = psi0.imaginary/psi0.real, y3[1] = psi.imaginary/psi.real, y3[2] = psi2.imaginary/psi2.real;
  fit_polynomial(x3, y3, abc, 2);
  d1 = nth_derivative_polynomial(abc, 2, 2, 1);


  velocities[0] =  A*d0;
  velocities[1] =  A*d1;

  if(velocity_perturbation){
    velocities[0] += -v_perturb*H_BAR*y/(mass[0]*(x*x+y*y));
    velocities[1] += v_perturb*H_BAR*x/(mass[0]*(x*x+y*y));

  }


  // printf("v0: %f\n", velocities[0] );
  // printf("rho: %f\n", pow(psi.real,2)+pow(psi.imaginary,2));
  delete[] x3, delete[] y3, delete[] abc;
  return pow(psi.real,2)+pow(psi.imaginary,2);
}
