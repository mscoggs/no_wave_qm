#include <fstream>
#include <iostream>
#include <string>
#include <algorithm>
#include <complex.h>
#include <math.h>

#include "grid.h"
#include "point.h"
#include "math_functions.h"



void Grid::evolve(){
    while(time < total_time){
      //evolve n stuff, update vs and rhos for each point!

      time += time_step;
    }
}




void Grid::init_grid(){
  get_data();
  set_points();
}



void Grid::get_data(){
  std::ifstream input("../inputs/input.txt");
	if(input.is_open()){
		std::string line, delimiter = " ", token;
		size_t pos = 0;
		int i =0, j=0;
		while(getline(input, line)){

      if(line=="GRID_LENGTH"){

				getline(input, line);
        grid_length = stoi(line);
				continue;
			}
			if(line=="SPATIAL_DIMENSION"){
				getline(input, line);
				spatial_dimension = stoi(line);
				continue;
			}
			if(line=="NUMBER_OF_PARTICLES"){
				getline(input, line);
				num_particles = stoi(line);
        config_dimension = num_particles*spatial_dimension;
        mass = new double[config_dimension]();
				continue;
			}

			if(line=="MASS"){
				getline(input, line);
				for(i=0; i<num_particles; i++){
						pos = line.find(delimiter);
						token = line.substr(0, pos);
            for(j=0; j<spatial_dimension; j++) mass[i*spatial_dimension + j] = stod(token);
						line.erase(0, pos + delimiter.length());
				}
        if(line=="TOTAL_TIME"){
          getline(input, line);
          total_time = stod(line);
          time = 0;
          continue;
        }
        if(line=="TIME_STEP"){
          getline(input, line);
          time_step = stod(line);
          continue;
        }
			}}}
	input.close();
}



void Grid::set_points(){
  int i;
  total_points = pow(grid_length, config_dimension);
  points = new Point[total_points]();

  for(i = 0; i<total_points; i++){
    points[i].init_point(config_dimension, grid_length, i);
  }
}



void Grid::print_grid_info(){
  printf("\n\n");
  printf("|| spatial_dimension         %4i || num_particles     %4i ||\n", spatial_dimension, num_particles);
  printf("|| config_space_dimension    %4i || grid_length       %4i ||\n", config_dimension, grid_length);
  printf("|| mass: [");
  for(int i = 0; i<config_dimension; i++) printf("%3.2f, ", mass[i]);
  printf("]\n\n\n");
}
