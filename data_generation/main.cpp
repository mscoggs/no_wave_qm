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




main(int argc, char *argv[]){
  int i, j;
	Grid grid;


  printf("\n...Initializing Grid...\n\n");
  grid.init_grid();
  printf("\n...Evolving System...\n\n\n");
  grid.evolve();
  printf("\n...Done Evolving...\n\n");
  exit(0);
}
