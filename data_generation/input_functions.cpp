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

  //
  //
  // exit(0);



  return pow(psi.real,2)+pow(psi.imaginary,2);
}



//
// void calc_velocities_initial(int *coordinates, double *velocities, int config_dimension, int vel_function,  int grid_length, double *mass){
//   int i, j;
//   double x,y, center, k, constant,m,h,w, exp1,r, x0,x1,y0,y1;
//   center = (grid_length-1.0)/2.0;
//   x = static_cast<double>(coordinates[0])-center;
//   y = static_cast<double>(coordinates[1])-center;
//
//
//   switch(vel_function){
//     case 1: //plane_wave
//         for(i=0; i<config_dimension; i++) velocities[i] = 0.1*coordinates[0] + 0.3*coordinates[1]*coordinates[0];
//         break;
//     case 2:  // holland_ex4.11_pg169
//         k = sqrt(2*OMEGA*mass[0]/H_BAR);
//         velocities[0] = -(H_BAR/mass[0])*(y/(x*x + y*y));
//         velocities[1] =  (H_BAR/mass[1])*(x/(x*x + y*y) + k);
//         break;
//     case 3: //2d_harmonic_oscillator
//         velocities[0] = 0;
//         velocities[1] =  (pow(H_BAR/mass[0], 1.5)/pow(2*OMEGA, 0.5)) / (y*y + H_BAR/2*OMEGA*mass[0]);
//         break;
//     case 4: //2d_harmonic_oscillator #2
//         m = mass[0], w = OMEGA, h = H_BAR;
//
//         velocities[0] = (15*h*y*(m*w*x*x-h))/(8*m*m*w*w*pow(x,6) - 48*h*m*w*pow(x,4) + 72*h*h*x*x + 3*h*h*y*y);
//         velocities[1] = 1.22/((2*m*w*pow(x,3)/h - 6*x)*(3*y*y/(2*pow(2*m*w*pow(x,3)/h - 6*x,2)) + 1));
//         break;
//     case 5:  //inf_circular_well
//         if(sqrt(x*x + y*y) < INF_WELL_RADIUS){
//           velocities[0] = -(H_BAR/mass[0])*(y/(x*x + y*y));
//           velocities[1] =  (H_BAR/mass[1])*(x/(x*x + y*y));
//           break;
//         }else{
//           velocities[0] = 0;
//           velocities[1] = 0;
//           break;
//         }
//       case 6: //2D Hydrogen # 1
//         velocities[0] = 0;
//         velocities[1] = 0;
//         break;
//       case 7: //2D Hydrogen #2
//         velocities[0] = -(H_BAR/mass[0])*(y/(x*x + y*y));
//         velocities[1] =  (H_BAR/mass[1])*(x/(x*x + y*y));
//         break;
//       case 8: //2D Hydrogen #3
//         r = sqrt(x*x + y*y);
//         exp1 = exp(-2*B1*r/3);
//         velocities[0] = -(H_BAR/mass[0])*(1/(1+ 486*exp1/(pow(r,2)*pow(B1,2))))*(sqrt(486)/B1)*(x*exp1*(-2*B1*r/3 - 1)/pow(r,3));
//         velocities[1] = -(H_BAR/mass[0])*(1/(1+ 486*exp1/(pow(r,2)*pow(B1,2))))*(y*exp1*(-2*B1*r/3 - 1)/pow(r,3));
//         break;
//       case 9:
//         velocities[0] = -(H_BAR/mass[0])*(y/(x*x + y*y));
//         velocities[1] =  (H_BAR/mass[1])*(x/(x*x + y*y));
//         break;
//       case 10:
//         velocities[0] = (H_BAR/mass[0])*(x/(x*x + y*y));
//         velocities[1] =  (H_BAR/mass[1])*(y/(x*x + y*y));
//         break;
//       case 11:
//         x = x+center;
//         x0 = x-center/2;
//         x1 = x-1.5*center;
//
//         velocities[1] = (H_BAR/mass[0])*((x0/(x0*x0 + y*y)) + (x1/(x1*x1 + y*y)));
//         velocities[0] =  -(H_BAR/mass[1])*((y/(x0*x0 + y*y)) + (y/(x1*x1 + y*y)));
//         break;
//       case 12:
//
//         velocities[1] = 0;
//         velocities[0] =  0;
//         break;
//       case 13: //numerical approach to the  infinite well perturbation
//         r = sqrt(x*x + y*y);
//         double *x3,*y3, *abc, A, epsilon = 1.0545*pow(10,-5), derivative;
//         A = (H_BAR/mass[0])/(1 + pow((j_d(x,y,1)*(y/r) + LAMBDA * j_d(x,y,1.5)*sin(1.5*get_theta(x,y)))/(j_d(x,y,1)*(x/r) + LAMBDA * j_d(x,y,1.5)*cos(1.5*get_theta(x,y))),2));
//         x3   = new double[3]();
//         y3   = new double[3]();
//         abc = new double[3]();
//
//         x3[0] = x-epsilon, x3[1] = x, x3[2] = x+epsilon;
//         for(i=0; i<3;i++) y3[i] = (j_d(x3[i],y,1)*(y/r) + LAMBDA * j_d(x3[i],y,1.5)*sin(1.5*get_theta(x3[i],y)))/(j_d(x3[i],y,1)*(x3[i]/r) + LAMBDA * j_d(x3[i],y,1.5)*cos(1.5*get_theta(x3[i],y)));
//         fit_polynomial(x3, y3, abc, 2);
//         derivative = nth_derivative_polynomial(abc, x, 2, 1);
//         velocities[0] =  A*derivative;
//
//
//         x3[0] = y-epsilon, x3[1] = y, x3[2] = y+epsilon;
//         for(i=0; i<3;i++) y3[i] = (j_d(x,x3[i],1)*(x3[i]/r) + LAMBDA * j_d(x ,x3[i],1.5)*sin(1.5*get_theta(x,x3[i])))/(j_d(x,x3[i],1)*(x/r) + LAMBDA * j_d(x,x3[i],1.5)*cos(1.5*get_theta(x,x3[i])));
//         fit_polynomial(x3, y3, abc, 2);
//         derivative = nth_derivative_polynomial(abc, y, 2, 1);
//         velocities[1] =  A*derivative;
//         break;
//
//
//
//
//
//
//   }
// }




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
