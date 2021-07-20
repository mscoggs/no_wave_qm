#include <math.h>
#include <fstream>


#include "input_functions.h"
#include "math_functions.h"


double calc_rho_initial(int *coordinates, int config_dimension, int grid_length, int rho_function, double *mass){
  int i;
  double rho=0, x, y, norm, center, A,r,moh=0,x0,x1,r0,r1;
  center = (grid_length-1.0)/2.0;
  x = static_cast<double>(coordinates[0])-center;
  y = static_cast<double>(coordinates[1])-center;

  switch(rho_function){
    case 1: //plane_wave
      return 1.0;
    case 2: //holland_ex4.11_pg169

      norm = 1.5/(static_cast<double>(grid_length*grid_length*grid_length*grid_length));
      return (x*x+y*y)*norm;
      case 3: //2d_harmonic_oscillator
        norm = 0.5;
        A = (2/3.14)*(mass[0]*OMEGA/H_BAR)*(mass[0]*OMEGA/H_BAR)*norm;
        return A*x*x*((2*mass[0]*OMEGA/H_BAR)*y*y + 1)*exp((-mass[0]*OMEGA/(H_BAR)) * (x*x + y*y));
    case 4: //2d_harmonic_oscillator #2
          norm = 0.5;
          moh = mass[0]*OMEGA/H_BAR;
          A = ((2/3.14)*pow(moh, 2))*norm*exp((-mass[0]*OMEGA/H_BAR) * (x*x + y*y));
          return A*(0.66*pow(moh, 2)*pow(x,6) - 2*moh*pow(x,4) + 1.5*pow(x,2) + pow(y,2));
    case 5:  //inf_circular_well
        if(sqrt(x*x + y*y) < INF_WELL_RADIUS){
          return pow(j_n((10.173/INF_WELL_RADIUS)*sqrt(x*x + y*y), 1),2);
        }else{
          return 0;
        }
    case 6: //2D Hydrogen #1
      r = sqrt(x*x + y*y);
      return (pow(B1,2)/(2*3.14)) * exp(-B1*r);
    case 7: //2D Hydrogen #2
      r = sqrt(x*x + y*y);
      return (r*r*pow(B2,4)/(12*3.14)) * exp(-B2*r);
    case 8: ///2D Hydrogen #3
      r = sqrt(x*x + y*y);
      return (r*r*pow(B2,4)/(12*3.14)) * exp(-B2*r) + (pow(B1,2)/(2*3.14)) * exp(-B1*r);
    case 9:
      return 1.0/(grid_length*grid_length);
    case 11:
      x = x+center;
      x0 = x-center/2;
      x1 = x-1.5*center;
      r1 = sqrt(x1*x1 + y*y);
      r0 = sqrt(x0*x0 + y*y);
      return (pow(B1,2)/(2*3.14)) * (exp(-r1)+exp(-r0));
    case 12:
      r = sqrt(x*x + y*y);
      return pow(j_n_non_int((10.904/INF_WELL_RADIUS)*r,1.5),2);
    case 13:
      r = sqrt(x*x + y*y);
      double j1 = j_d(x,y,1);
      double j1p5 = j_d(x,y,1.5);
      return pow(j1,2) + 2*LAMBDA*j1*j1p5*cos(0.5*get_theta(x,y)) + pow(LAMBDA*j1p5,2);
  }
}




void calc_velocities_initial(int *coordinates, double *velocities, int config_dimension, int vel_function,  int grid_length, double *mass){
  int i, j;
  double x,y, center, k, constant,m,h,w, exp1,r, x0,x1,y0,y1;
  center = (grid_length-1.0)/2.0;
  x = static_cast<double>(coordinates[0])-center;
  y = static_cast<double>(coordinates[1])-center;


  switch(vel_function){
    case 1: //plane_wave
        for(i=0; i<config_dimension; i++) velocities[i] = 0.1*coordinates[0] + 0.3*coordinates[1]*coordinates[0];
        break;
    case 2:  // holland_ex4.11_pg169
        k = sqrt(2*OMEGA*mass[0]/H_BAR);
        velocities[0] = -(H_BAR/mass[0])*(y/(x*x + y*y));
        velocities[1] =  (H_BAR/mass[1])*(x/(x*x + y*y) + k);
        break;
    case 3: //2d_harmonic_oscillator
        velocities[0] = 0;
        velocities[1] =  (pow(H_BAR/mass[0], 1.5)/pow(2*OMEGA, 0.5)) / (y*y + H_BAR/2*OMEGA*mass[0]);
        break;
    case 4: //2d_harmonic_oscillator #2
        m = mass[0], w = OMEGA, h = H_BAR;

        velocities[0] = (15*h*y*(m*w*x*x-h))/(8*m*m*w*w*pow(x,6) - 48*h*m*w*pow(x,4) + 72*h*h*x*x + 3*h*h*y*y);
        velocities[1] = 1.22/((2*m*w*pow(x,3)/h - 6*x)*(3*y*y/(2*pow(2*m*w*pow(x,3)/h - 6*x,2)) + 1));
        break;
    case 5:  //inf_circular_well
        if(sqrt(x*x + y*y) < INF_WELL_RADIUS){
          velocities[0] = -(H_BAR/mass[0])*(y/(x*x + y*y));
          velocities[1] =  (H_BAR/mass[1])*(x/(x*x + y*y));
          break;
        }else{
          velocities[0] = 0;
          velocities[1] = 0;
          break;
        }
      case 6: //2D Hydrogen # 1
        velocities[0] = 0;
        velocities[1] = 0;
        break;
      case 7: //2D Hydrogen #2
        velocities[0] = -(H_BAR/mass[0])*(y/(x*x + y*y));
        velocities[1] =  (H_BAR/mass[1])*(x/(x*x + y*y));
        break;
      case 8: //2D Hydrogen #3
        r = sqrt(x*x + y*y);
        exp1 = exp(-2*B1*r/3);
        velocities[0] = -(H_BAR/mass[0])*(1/(1+ 486*exp1/(pow(r,2)*pow(B1,2))))*(sqrt(486)/B1)*(x*exp1*(-2*B1*r/3 - 1)/pow(r,3));
        velocities[1] = -(H_BAR/mass[0])*(1/(1+ 486*exp1/(pow(r,2)*pow(B1,2))))*(y*exp1*(-2*B1*r/3 - 1)/pow(r,3));
        break;
      case 9:
        velocities[0] = -(H_BAR/mass[0])*(y/(x*x + y*y));
        velocities[1] =  (H_BAR/mass[1])*(x/(x*x + y*y));
        break;
      case 10:
        velocities[0] = (H_BAR/mass[0])*(x/(x*x + y*y));
        velocities[1] =  (H_BAR/mass[1])*(y/(x*x + y*y));
        break;
      case 11:
        x = x+center;
        x0 = x-center/2;
        x1 = x-1.5*center;

        velocities[1] = (H_BAR/mass[0])*((x0/(x0*x0 + y*y)) + (x1/(x1*x1 + y*y)));
        velocities[0] =  -(H_BAR/mass[1])*((y/(x0*x0 + y*y)) + (y/(x1*x1 + y*y)));
        break;
      case 12:

        velocities[1] = 0;
        velocities[0] =  0;
        break;
      case 13: //numerical approach to the  infinite well perturbation
        r = sqrt(x*x + y*y);
        double *x3,*y3, *abc, A, epsilon = 1.0545*pow(10,-5), derivative;
        A = (H_BAR/mass[0])/(1 + pow((j_d(x,y,1)*(y/r) + LAMBDA * j_d(x,y,1.5)*sin(1.5*get_theta(x,y)))/(j_d(x,y,1)*(x/r) + LAMBDA * j_d(x,y,1.5)*cos(1.5*get_theta(x,y))),2));
        x3   = new double[3]();
        y3   = new double[3]();
        abc = new double[3]();

        x3[0] = x-epsilon, x3[1] = x, x3[2] = x+epsilon;
        for(i=0; i<3;i++) y3[i] = (j_d(x3[i],y,1)*(y/r) + LAMBDA * j_d(x3[i],y,1.5)*sin(1.5*get_theta(x3[i],y)))/(j_d(x3[i],y,1)*(x3[i]/r) + LAMBDA * j_d(x3[i],y,1.5)*cos(1.5*get_theta(x3[i],y)));
        fit_polynomial(x3, y3, abc, 2);
        derivative = nth_derivative_polynomial(abc, x, 2, 1);
        velocities[0] =  A*derivative;


        x3[0] = y-epsilon, x3[1] = y, x3[2] = y+epsilon;
        for(i=0; i<3;i++) y3[i] = (j_d(x,x3[i],1)*(x3[i]/r) + LAMBDA * j_d(x ,x3[i],1.5)*sin(1.5*get_theta(x,x3[i])))/(j_d(x,x3[i],1)*(x/r) + LAMBDA * j_d(x,x3[i],1.5)*cos(1.5*get_theta(x,x3[i])));
        fit_polynomial(x3, y3, abc, 2);
        derivative = nth_derivative_polynomial(abc, y, 2, 1);
        velocities[1] =  A*derivative;
        break;






  }
}




double calc_potential(double time, int *coordinates, int num_particles, int spatial_dimension, int potential_function, int grid_length, double *mass){
  double x, y, center, r;
  double e = 0.5, V=0, *positioni, *positionj, distance =0;

  center = (grid_length-1.0)/2.0;
  x = static_cast<double>(coordinates[0])-center;
  y = static_cast<double>(coordinates[1])-center;

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
      r = sqrt(x*x + y*y);
      return -(Z*E*E)/r;
    case 9:
      return 0;
    case 10:
      return 0;
  }
}
