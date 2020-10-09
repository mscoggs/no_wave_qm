#include <fstream>
#include <iostream>
#include <string>
#include <algorithm>
#include <complex.h>
#include <math.h>


#include "point.h"
#include "grid.h"
#include "math_functions.h"
#include "input_functions.h"


/*
g++ -o main main.cpp grid.cpp point.cpp math_functions.cpp -std=c++11
./main

OR

make run
*/

main(int argc, char *argv[]){
  int i, j;
	Grid grid;


  printf("\n...Initializing Grid...\n\n");
  grid.init_grid();
  printf("\n...Evolving System...\n\n");
  grid.evolve();
  printf("\n...Done Evolving...\n\n");
  exit(0);
}
