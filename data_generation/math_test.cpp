#include <fstream>
#include <iostream>
#include <string>
#include <algorithm>
#include <complex.h>
#include <math.h>


#include "math_functions.h"
#include "input_functions.h"




main(int argc, char *argv[]){

  printf("non_int_hermite(2,4) = %f\n", non_int_hermite(0.5,4,5));
  printf("non_int_hermite(2,4) = %f\n", non_int_hermite(0.5,4,10));
  printf("non_int_hermite(2,4) = %f\n", non_int_hermite(0.5,4,20));
  printf("non_int_hermite(2,4) = %f\n", non_int_hermite(0.5,4,100));
  // printf("non_int_hermite(5.7,4) = %f\n", non_int_hermite(5.7,4));
  // printf("non_int_hermite(1,3.1) = %f\n", non_int_hermite(1,3.1));
  // printf("non_int_hermite(1,3) = %f\n", non_int_hermite(1,3));
  // printf("non_int_hermite(6,3.1) = %f\n", non_int_hermite(6,3.1));
  // printf("non_int_hermite(6,3) = %f\n", non_int_hermite(6,3));
}
