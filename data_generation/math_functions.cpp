#include <math.h>

#include "math_functions.h"


double calc_distance(double *position_1, double *position_2){
  int i;
  double distance = 0;
  for(i=0; i<sizeof(position_1)/sizeof(double); i++){
    distance += pow(position_2[i]-position_1[i],2);
  }
  return sqrt(distance);
}





double calc_velocity_i(int *coordinates, double time, int i){
  //Specify the velocity equation, eventually this will get replaced with a function that comes from the input file
  double v = 0;
  int j;
  for(j=0; j<sizeof(coordinates)/sizeof(double); j++) v+= (coordinates[j]+time)/(j+1) + i/time;
  return v;
}



double calc_rho(){
  //Specify the rho equation, eventually this will get replaced with a function that comes from the input file
  return 0;
}



double calc_potential(){
  //Specify the V equation, eventually this will get replaced with a function that comes from the input file
  return 0;
}


double calc_quantum_potential(){
  //Specify the Q equation
  return 0;
}

double gradient_i(double (*f)(double), double x, int i){
  return 0;
}


double divergence(double (*f)(double), double x){
  return 0;
}


double laplacian(double (*f)(double), double x){
  return 0;
}


double derivative_i(double (*f)(double), double x, int i){
  return 0;
}
