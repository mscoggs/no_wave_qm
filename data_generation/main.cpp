#include <fstream>
#include <iostream>
#include <string>
#include <algorithm>
#include <complex.h>
#include <math.h>

#include "point.h"
#include "grid.h"
#include "math_functions.h"

//g++ -o main main.cpp grid.cpp point.cpp math_functions.cpp -std=c++11
//./main





//Structure of this program:
//A file holding all of the mathematical functions (as a way to bypass inputting a function through a  text file)
//main
//a file holding the grids class
//the grid object holds a point class
//each Point holds the things we wish to evolve over time




main(int argc, char *argv[])
{
  int i, j;
	Grid grid;
  grid.init_grid();
  grid.print_grid_info();
  grid.evolve();

}
