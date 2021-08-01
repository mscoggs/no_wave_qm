#include <math.h>
#include <iostream>
#include <string.h>
#include <cstring>
#include <random>


#include "math_functions.h"


double factorial_double(double n){
  double step = 0.001;
  double sum = 0, lhs=0, rhs = step;
  double int_limit = 1000;

  while(rhs < int_limit){
    sum += ((pow(lhs,n)/exp(lhs)) + (pow(rhs,n)/exp(rhs)))/2;
    lhs = rhs;
    rhs += step;
  }
  return sum/int_limit;
}




void fit_polynomial(double *x, double *y, double *coeffs, int degree){
  int i,j,k, n = degree, N = 3;
  double *sigma_x, *sigma_y, temp, t;
  bool skip = true;



  //If all of the values are equal, we skip the fit and return a constant
  for (i=0;i<N-1;i++) if(y[i] != y[(i+1)]) skip = false;
  if(skip){
    coeffs[0] = y[0];
    return;
  }
  if(x[0] == x[1]){
    x[1] = x[2];
    y[1] = y[2];
    N = 2;
    n = 1;
  }
  if(x[1] == x[2]){
    N = 2;
    n = 1;
  }

  double norm_aug[n+1][n+2];
  sigma_x = new double[2*n+1]();
  sigma_y = new double[n+1]();



  for (i=0;i<2*n+1;i++) for (j=0;j<N;j++) sigma_x[i]=sigma_x[i]+pow(x[j],i);
  for (i=0;i<=n;i++) for (j=0;j<=n;j++) norm_aug[i][j]=sigma_x[i+j];
  for (i=0;i<n+1;i++) for (j=0;j<N;j++) sigma_y[i]=sigma_y[i]+pow(x[j],i)*y[j];
  for (i=0;i<=n;i++) norm_aug[i][n+1]=sigma_y[i];
  n=n+1;

  for (i=0;i<n;i++){
      for (k=i+1;k<n;k++){
          if (norm_aug[i][i]<norm_aug[k][i]){
              for (j=0;j<=n;j++){
                  temp=norm_aug[i][j];
                  norm_aug[i][j]=norm_aug[k][j];
                  norm_aug[k][j]=temp;
              }
          }
      }
  }
  for (i=0;i<n-1;i++){
      for (k=i+1;k<n;k++){
              t=norm_aug[k][i]/norm_aug[i][i];
              for (j=0;j<=n;j++) norm_aug[k][j]=norm_aug[k][j]-t*norm_aug[i][j];
      }
  }
  for (i=n-1;i>=0;i--){
      coeffs[i]=norm_aug[i][n];
      for (j=0;j<n;j++) if (j!=i) coeffs[i]=coeffs[i]-norm_aug[i][j]*coeffs[j];
      coeffs[i]=coeffs[i]/norm_aug[i][i];
  }
  delete[] sigma_x, delete[] sigma_y;
}




double nth_derivative_polynomial(double *coeffs, double x_val, int degree_of_fit, int n){
  int i, j;
  double value, value_per_degree;

  if(n == 1) return 2*coeffs[2]*x_val + coeffs[1];
  if(n == 2) return 2*coeffs[2];
  // for(i=n; i<degree_of_fit+1; i++){
  //   value_per_degree = coeffs[i]*pow(x_val,i-n);
  //   for(j=n; j>0; j--) value_per_degree = value_per_degree*n;
  //   value += value_per_degree;
  // }
  // return value;
}




int coordinates_to_index(int* coordinates, int config_dimension, int grid_length){
  int i, index=0;


  for(i=0;i<config_dimension; i++) index+=(coordinates[i])*int_pow(grid_length, i);
  return  index;
}




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




double calc_distance(double *position_1, double *position_2, int spatial_dimension){
  int i;
  double distance = 0;


  for(i=0; i<spatial_dimension; i++){
    distance += int_pow(position_2[i]-position_1[i],2);
  }
  return sqrt(distance);
}




int int_pow(int x, int p) {


  if (p == 0) return 1;
  if (p == 1) return x;
  return x * int_pow(x, p-1);
}




void print_load_bar(double progress){
  int barWidth = 70, pos, i;


  std::cout << "      [";
  pos = barWidth * progress;
  for (i=0;i<barWidth; ++i) {
      if (i < pos) std::cout << "=";
      else if (i == pos) std::cout << ">";
      else std::cout << " ";
  }
  std::cout << "] " << int(progress * 100.0) << " %\r";
  std::cout.flush();
}



void print_double_array(double *array, int size){
   int i = 0;

   printf("\n\nPrinting the array of size %i\n[", size);
   for(i=0;i<size;i++){
	   printf("%f, ",array[i]);
   }
   printf("]\n");
}



double j_1(double x, int grid_length){
  double sum = 0, term_i;
  int N = grid_length*2, i,j;
  for(i=0;i<N;i++){
    term_i = x/2.0;
    for(j=1;j<=i;j++) term_i = term_i*x*x/(4.0*j*(j+1));
    if(i%2 == 1) term_i = -term_i;
    sum += term_i;
  }
  return sum;
}


double j_n(double x, int n){
  double o=0, sum =0, left, right, step = 0.0001;
  while(o<PI){
    left = cos(x*sin(o) - n*o);
    right = cos(x*sin(o+step) - n*(o+step));
    sum += (right+left)*step/2;
    o += step;
  }
  return sum/PI;
}


double j_n_non_int(double x, double n){
  double sum1 =0, sum2 =0, left, right, step = 0.001;
  double o=step;
  left = 1;
  while(o<=PI+step/2){
    right = cos(x*sin(o) - n*(o));
    sum1 += (right+left)*step/2;
    left = right;
    o += step;
  }


  o=step;
  left = 1;
  while(o<=20){
    right = exp(-n*o - x*sinh(o));
    sum2 += (right+left)*step/2;
    left = right;
    o += step;
  }

  sum2 = sum2 * sin(n*PI);
  return (sum1-sum2)/PI;
}


complex f(double t, double x, double n){
  complex en, ef, power, result;
  power.re = -cos(4*PI*t)+2*cos(2*PI*t)*x;
  power.im = -sin(4*PI*t)+2*sin(2*PI*t)*x;


  ef.re = exp(power.re) * cos(power.im );
  ef.im = exp(power.re) * sin(power.im );

  en.re = cos(-2*PI*n*t);
  en.im = sin(-2*PI*n*t);


  result.re = ef.re*en.re - ef.im*en.im;
  result.im = ef.im*en.re + ef.re*en.im;
  return result;
}




complex f2(double t, double x, double n){
  complex en, ef, power, result;




  power.re = -cos(4*PI*t)+2*cos(2*PI*t)*x;
  power.im = -sin(4*PI*t)+2*sin(2*PI*t)*x;


  result.re = cos(-2*PI*n*t + power.im);
  result.im = sin(-2*PI*n*t + power.im);

  result.re = result.re *exp(power.re);
  result.im = result.im *exp(power.re);

  return result;
}


complex hermite_contour(double x, double y, double n){


    double r = sqrt(x*x+y*y);
    double t = 0, stop = 1, dt = 0.01;
    complex sum, result;

    sum.re = 0, sum.im = 0;

    while(t < stop){
      result = f(t,r,n);
      //printf("result: %f, %f\n", result.re, result.im);
      sum.re += result.re*dt;
      sum.im += result.im*dt;
      //printf("sum: %f, %f\n\n\n", sum.re, sum.im);

      t+=dt;
    }
    sum.re = sum.re * factorial_double(n);
    sum.im = sum.im * factorial_double(n);

    //printf("r, n : %f, %f \n", r, n);
    // printf("sum.re : %f \n", sum.re);
    // printf("sum.im : %f \n", sum.im);
    return sum;
}





double j_d(double x, double y, double n){
  double sum1 =0, sum2 =0, left, right, step = 0.001;
  double o=step;
  double z = 0;

  if(n == 1.5) z = 10.904;
  else if(n == 1) z = 10.173;
  else{
    printf("No Z_nm for the given n\n\n");
    exit(0);
  }
  x = (z/INF_WELL_RADIUS)*sqrt(x*x + y*y);
  left = 1;
  while(o<=PI+step/2){
    right = cos(x*sin(o) - n*(o));
    sum1 += (right+left)*step/2;
    left = right;
    o += step;
  }



  o=step;
  left = 1;
  while(o<=20){
    right = exp(-n*o - x*sinh(o));
    sum2 += (right+left)*step/2;
    left = right;
    o += step;
  }

  sum2 = sum2 * sin(n*PI);
  return (sum1-sum2)/PI;
}

int factorial(int n){
  int total = 1;
  int i = 1;
  while(i <= n){
    total = total*i;
    i+=1;
  }
  return total;
}

double get_random_double(double lower_bound, double upper_bound){
  std::uniform_real_distribution<double> unif(lower_bound,upper_bound);
  std::default_random_engine re;
  double a_random_double = unif(re);
  return a_random_double;
}

double get_theta(double x, double y){
  double r = sqrt(x*x+y*y);
  double theta = asin(y/r);
  if(x < 0) theta = PI -theta;
  return theta;
}

// double sin_ao(double x, double y, double a){
//   return sin(a*get_theta(x,y));
// }
//
// double cos_ao(double x, double y, double a){
//   double r = sqrt(x*x+y*y);
//
//   double theta = arccos(y/r);
//   if(x < 0) theta += PI;
//   return sin(a*theta);
// }
double non_int_hermite(double x, double n, int mult){
  printf("WARNING, THIS FUNCTION IS NOT WORKING");
  double total =  0;
  double term_total = 0;
  //int mult = 10;
  int limit = (int)ceil(n*mult);
  int pow_2 = 1;
  int i,j;
  double n_prod;
  for(i=0; i<limit;i++){


    if(i%2 == 0) pow_2 = i/2;
    else pow_2 = (i-1)/2;



    n_prod=1;
    for(j=0; j<i-1;j++){
      if((i-j)%2 == 0) n_prod = n_prod*(j-n);
    }


    term_total = n_prod*pow(x,i)*pow(2,pow_2)/factorial(i);
    // printf("n_prod: %f  ",n_prod);
    // printf("pow(x,i),pow_2%f %i",pow(x,i),pow_2);
    // //printf("pow_2: %f   ",pow_2);
    // printf("pow(2,pow_2): %f   ",pow(2,pow_2));
    // printf("factorial(i): %i   \n\n",factorial(i));

    total += term_total;
  }
  //exit(0);
  return total;
}


// double divergence(double (*f)(int*, double, int, int), int * coordinates, double time, int config_dimension){
//   int i;
//   double divergence =0;
//   for(i=0; i<config_dimension; i++){
//     divergence += derivative_i(calc_rho_initial_times_v, coordinates, time, config_dimension, EPSILON, i);
//   }
//   return divergence;
// }
//
//
//
//
// double derivative_i_potential(double (*f)(int*, int, int), int *coordinates, int config_dimension, int num_particles, int spatial_dimension, double epsilon, int index){
//   double rise=0, run = 2*epsilon;
//   int *coord_upper, *coord_lower;
//   coord_upper = new int[config_dimension]();
//   coord_lower = new int[config_dimension]();
//   for(int i = 0; i< config_dimension; i++){
//     coord_upper[i] = coordinates[i];
//     coord_lower[i] = coordinates[i];
//   }
//   coord_upper[index] += epsilon;
//   coord_lower[index] -= epsilon;
//
//   rise = f(coord_upper, num_particles, spatial_dimension) - f(coord_lower, num_particles, spatial_dimension);
//   return rise/run;
// }
//
//
//
//
// double derivative_i_quantum_potential(double (*f)(double*, int*, double, int), double* mass, int *coordinates, int config_dimension, double time , double epsilon, int index){
//   double rise=0, run = 2*epsilon;
//   int *coord_upper, *coord_lower;
//   coord_upper = new int[config_dimension]();
//   coord_lower = new int[config_dimension]();
//   for(int i = 0; i< config_dimension; i++){
//     coord_upper[i] = coordinates[i];
//     coord_lower[i] = coordinates[i];
//   }
//   coord_upper[index] += epsilon;
//   coord_lower[index] -= epsilon;
//
//   rise = f(mass, coord_upper, time, config_dimension) - f(mass, coord_lower, time, config_dimension);
//   return rise/run;
// }
//
//
//
//
// double derivative_i_velocity(double (*f)(int*, double, int, int), int *coordinates, double time, int config_dimension, double epsilon, int index, int j){
//   double rise=0, run = 2*epsilon;
//   int *coord_upper, *coord_lower;
//   coord_upper = new int[config_dimension]();
//   coord_lower = new int[config_dimension]();
//   for(int i = 0; i< config_dimension; i++){
//     coord_upper[i] = coordinates[i];
//     coord_lower[i] = coordinates[i];
//   }
//   coord_upper[index] += epsilon;
//   coord_lower[index] -= epsilon;
//
//   rise = f(coord_upper, time, config_dimension, j) - f(coord_lower, time, config_dimension, j);
//   return rise/run;
// }
//
//
//
//
// double derivative_i(double (*f)(int*, double, int, int), int *coordinates, double time, int config_dimension, double epsilon, int index){
//   double rise=0, run = 2*epsilon;
//   int *coord_upper, *coord_lower;
//   coord_upper = new int[config_dimension]();
//   coord_lower = new int[config_dimension]();
//   for(int i = 0; i< config_dimension; i++){
//     coord_upper[i] = coordinates[i];
//     coord_lower[i] = coordinates[i];
//   }
//   coord_upper[index] += epsilon;
//   coord_lower[index] -= epsilon;
//
//   rise = f(coord_upper, time, config_dimension, index) - f(coord_lower, time, config_dimension, index);
//   return rise/run;
// }
