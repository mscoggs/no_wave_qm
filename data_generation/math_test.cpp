// #include <fstream>
// #include <iostream>
// #include <string>
// #include <algorithm>
// #include <complex.h>
// #include <math.h>
//
//
// #include "math_functions.h"
// #include "input_functions.h"
//
//
//
// 
// main(int argc, char *argv[]){
//   std::ofstream file;
//   complex result;
//   double n, x;
//
//   int upper_n, lower_n, upper_x, lower_x;
//   double step_n, step_x;
//
//
//   upper_n = 5;
//   lower_n = 0;
//   upper_x = 5;
//   lower_x = 0;
//   step_n = 0.1;
//   step_x = 0.1;
//
//   file.open("hermite_test.txt");
//   file << "n      x      real       imaginary";
//   for(n=lower_n; n<upper_n; n+=step_n){
//     for(x=lower_x; x<upper_x; x+=step_x){
//       printf("n, x: %f %f\n",n,x);
//       result = hermite_contour(x,0,n);
//       file << n << "     " << x << "     " << result.re << "     " << result.im << "\n";
//     }
//   }
// }


#include <random>

int main()
{
   double lower_bound = 0;
   double upper_bound = 10000;
   std::uniform_real_distribution<double> unif(lower_bound,upper_bound);
   std::default_random_engine re;
   double a_random_double = unif(re);

   return 0;
}
