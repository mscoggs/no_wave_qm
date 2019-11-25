#include <math.h>
#include <iostream>
#include <string.h>
#include <cstring>

#include "math_functions.h"



void index_to_coordinates(int* coordinates, int index, int config_dimension, int grid_length){
  int i=0,j, r;
  while(index != 0){
    r = index%grid_length;
    coordinates[i] = r;
    index = (index-r)/grid_length;
    i++;
  }
  for(j=i;j<config_dimension;j++) coordinates[j] = 0;
}

int coordinates_to_index(int* coordinates, int config_dimension, int grid_length){
  int i, index;
  for(i=0;i<config_dimension; i++) index+=coordinates[i]*pow(grid_length, i);
  return  index;
}




void fit_polynomial(int *x, double *y, double *coeffs, int degree){
  int i,j,k, n = degree, N = sizeof(x)/sizeof(int);
  double *sigma_x, *sigma_y, temp, t;
  sigma_x = new double[2*n+1]();
  sigma_y = new double[n+1]();
  double B[n+1][n+2];

  for (i=0;i<2*n+1;i++) for (j=0;j<N;j++) sigma_x[i]=sigma_x[i]+pow(x[j],i);
  for (i=0;i<=n;i++) for (j=0;j<=n;j++) B[i][j]=sigma_x[i+j];
  for (i=0;i<n+1;i++) for (j=0;j<N;j++) sigma_y[i]=sigma_y[i]+pow(x[j],i)*y[j];

  for (i=0;i<=n;i++) B[i][n+1]=sigma_y[i];
  n=n+1;

  for (i=0;i<n;i++){
      for (k=i+1;k<n;k++){
          if (B[i][i]<B[k][i]){
              for (j=0;j<=n;j++){
                  temp=B[i][j];
                  B[i][j]=B[k][j];
                  B[k][j]=temp;
              }
          }
      }
  }
  for (i=0;i<n-1;i++){
      for (k=i+1;k<n;k++){
              t=B[k][i]/B[i][i];
              for (j=0;j<=n;j++) B[k][j]=B[k][j]-t*B[i][j];
      }
  }
  for (i=n-1;i>=0;i--){
      coeffs[i]=B[i][n];
      for (j=0;j<n;j++) if (j!=i) coeffs[i]=coeffs[i]-B[i][j]*coeffs[j];
      coeffs[i]=coeffs[i]/B[i][i];
  }
}


double derivative_polynomial(double *coeffs, double x_val){
  int i;
  double value = coeffs[1];
  for(i=2; i<sizeof(coeffs)/sizeof(double); i++){
    value += i*coeffs[i]*pow(x_val,i-1);
  }
  return value;
}


double second_derivative_polynomial(double *coeffs, double x_val){
  int i;
  double value = coeffs[2];
  for(i=3; i<sizeof(coeffs)/sizeof(double); i++){
    value += (i-1)*i*coeffs[i]*pow(x_val, i-2);
  }
  return value;
}






//SPECIFY THE DENSITY EQUATION HERE
double calc_rho(int *coordinates, double time, int config_dimension){
  double spatial_dependence = 1, time_dependence =time/7;
  int i;
  for(i=0; i<config_dimension; i++) spatial_dependence += double(coordinates[i]*coordinates[i])/(i+1);
  return time_dependence+spatial_dependence;
}


//SPECIFY THE VELOCITY EQUATION HERE
double calc_velocity_i(int *coordinates, double time, int config_dimension, int i){
  double vi = 0;
  int j;

  for(j=0; j<config_dimension; j++){
    vi += (coordinates[j]+time)/(j+1);
  }
  return vi/(i+7);
}


//SPECIFY THE POTENTIAL EQUATION HERE
double calc_potential(int *coordinates, int num_particles, int spatial_dimension){
  // double e = 0.5, V=0, *positioni, *positionj, distance =0;
  // int i,j,k;
  // positioni = new double[spatial_dimension];
  // positionj = new double[spatial_dimension];
  //
  // for(i = 0; i<num_particles; i++){
  //
  //   for(k =0; k<spatial_dimension; k++) positioni[k] = coordinates[i*spatial_dimension + k];
  //
  //   for(j = i+1; j < num_particles; j++){
  //
  //     for(k =0; k<spatial_dimension; k++) positionj[k] = coordinates[j*spatial_dimension + k];
  //     distance = calc_distance(positioni, positionj, spatial_dimension);
  //     if(distance > 0) V += -pow(e,2)/distance;
  //   }
  // }
  // return V;
  return 0;
}











double calc_sqrt_rho(int *coordinates,double time,  int config_dimension, int i){
  return sqrt(calc_rho(coordinates, time, config_dimension));
}



double calc_rho_times_v(int *coordinates, double time, int config_dimension, int i){
  return calc_rho(coordinates, time, config_dimension)*calc_velocity_i(coordinates, time, config_dimension, i);
}





double calc_quantum_potential(double *mass, int *coordinates, double time, int config_dimension){
  double Q = 0, Q_i = 0;
  int i, j;
  for(i = 0; i<config_dimension; i++){

    Q_i = derivative_i(calc_sqrt_rho, coordinates, time, config_dimension, EPSILON, i) / mass[i];
    Q += Q_i;
  }

  //return (-pow(H_BAR,2) * Q) / (2 * sqrt(calc_rho(coordinates, time, velocities, config_dimension)));
  return -Q / (2 * sqrt(calc_rho(coordinates, time, config_dimension)));

}




void calc_velocities(int *coordinates, double time, double *velocities, int config_dimension){
  int i;
  for(i=0; i<config_dimension; i++){
    velocities[i] = calc_velocity_i(coordinates, time, config_dimension, i);
  }
}







double calc_distance(double *position_1, double *position_2, int spatial_dimension){
  int i;
  double distance = 0;
  for(i=0; i<spatial_dimension; i++){
    distance += pow(position_2[i]-position_1[i],2);
  }
  return sqrt(distance);
}





double divergence(double (*f)(int*, double, int, int), int * coordinates, double time, int config_dimension){
  int i;
  double divergence =0;
  for(i=0; i<config_dimension; i++){
    divergence += derivative_i(calc_rho_times_v, coordinates, time, config_dimension, EPSILON, i);
  }
  return divergence;
}




double derivative_i_potential(double (*f)(int*, int, int), int *coordinates, int config_dimension, int num_particles, int spatial_dimension, double epsilon, int index){
  double rise=0, run = 2*epsilon;
  int *coord_upper, *coord_lower;
  coord_upper = new int[config_dimension]();
  coord_lower = new int[config_dimension]();
  for(int i = 0; i< config_dimension; i++){
    coord_upper[i] = coordinates[i];
    coord_lower[i] = coordinates[i];
  }
  coord_upper[index] += epsilon;
  coord_lower[index] -= epsilon;

  rise = f(coord_upper, num_particles, spatial_dimension) - f(coord_lower, num_particles, spatial_dimension);
  return rise/run;
}





double derivative_i_quantum_potential(double (*f)(double*, int*, double, int), double* mass, int *coordinates, int config_dimension, double time , double epsilon, int index){
  double rise=0, run = 2*epsilon;
  int *coord_upper, *coord_lower;
  coord_upper = new int[config_dimension]();
  coord_lower = new int[config_dimension]();
  for(int i = 0; i< config_dimension; i++){
    coord_upper[i] = coordinates[i];
    coord_lower[i] = coordinates[i];
  }
  coord_upper[index] += epsilon;
  coord_lower[index] -= epsilon;

  rise = f(mass, coord_upper, time, config_dimension) - f(mass, coord_lower, time, config_dimension);
  return rise/run;
}




double derivative_i_velocity(double (*f)(int*, double, int, int), int *coordinates, double time, int config_dimension, double epsilon, int index, int j){
  double rise=0, run = 2*epsilon;
  int *coord_upper, *coord_lower;
  coord_upper = new int[config_dimension]();
  coord_lower = new int[config_dimension]();
  for(int i = 0; i< config_dimension; i++){
    coord_upper[i] = coordinates[i];
    coord_lower[i] = coordinates[i];
  }
  coord_upper[index] += epsilon;
  coord_lower[index] -= epsilon;

  rise = f(coord_upper, time, config_dimension, j) - f(coord_lower, time, config_dimension, j);
  return rise/run;
}



double derivative_i(double (*f)(int*, double, int, int), int *coordinates, double time, int config_dimension, double epsilon, int index){
  double rise=0, run = 2*epsilon;
  int *coord_upper, *coord_lower;
  coord_upper = new int[config_dimension]();
  coord_lower = new int[config_dimension]();
  for(int i = 0; i< config_dimension; i++){
    coord_upper[i] = coordinates[i];
    coord_lower[i] = coordinates[i];
  }
  coord_upper[index] += epsilon;
  coord_lower[index] -= epsilon;

  rise = f(coord_upper, time, config_dimension, index) - f(coord_lower, time, config_dimension, index);
  return rise/run;
}






double runge_kutta(double (*f)(int*, double, int, int), double time_step){
  // double t4, k1,k2,k3,k4, h, h2;
  // h = time_step;
  // h2 = h/2;
  // t4 = (k1+2*k2+2*k3+k$)/6
  // k1 = f(x   , y);
  // k2 = f(x+h2, y + h2*k1);
  // k3 = f(x+h2, y + h2*k2);
  // k4 = f(x+h , y + h *k3);
  //
  // return change;
}
