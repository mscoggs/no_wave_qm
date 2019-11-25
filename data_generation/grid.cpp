#include <fstream>
#include <iostream>
#include <string>
#include <algorithm>
#include <complex.h>
#include <math.h>
#include <string.h>

#include "grid.h"
#include "point.h"
#include "math_functions.h"
#include <algorithm>


void Grid::evolve(){
  int i;
  while(time < total_time){
      print_evolution_info();
      step_rho();
      //step_velocities();
      //update_Q();

      //old_to_current();

      time += time_step;
    }
}


void Grid::step_rho(){
  double dpdt, *y, *coeffs;
  int i, j, *coords_of_point, *coords_above, *coords_below, *x;
  Point point_to_evolve, other_point;

  x   = new int[3]();
  y   = new double[3]();
  coeffs = new double[3]();
  coords_of_point = new int[config_dimension]();
  coords_above    = new int[config_dimension]();
  coords_below    = new int[config_dimension]();


  for(i=0; i<total_points; i++){

    dpdt = 0;
    point_to_evolve = points[i];
    point_to_evolve.get_coordinates(coords_of_point, config_dimension);
    for(j=0; j<config_dimension; j++) printf("%i, ", coords_of_point[j]);

    for(j=0; j<config_dimension; j++){

      std::fill_n(x, 3, 0);
      std::fill_n(y, 3, point_to_evolve.get_rho_old()*point_to_evolve.get_velocity_old_i(j));
      std::fill_n(coeffs, 3, 0);
      memcpy(coords_above, coords_of_point, sizeof(coords_of_point));
      memcpy(coords_below, coords_of_point, sizeof(coords_of_point));

      if(coords_of_point[j] < grid_length - 1){
        x[2] =  1;
        coords_above[j]+=1;
        int index = coordinates_to_index(coords_above, config_dimension,grid_length);
        printf("index: %i", index);

        other_point =  points[index];
        printf("ITERATION j:%i\n", j);

        y[2] =other_point.get_rho_old()*other_point.get_velocity_old_i(j);
      }
      printf("ITERATION j:%i\n", j);

      if(coords_of_point[j] > 0) {
        x[0] =  -1;
        coords_below[j]-=1;
        other_point =  points[coordinates_to_index(coords_below, config_dimension,grid_length)];
        y[0] =other_point.get_rho_old()*other_point.get_velocity_old_i(j);
      }

      fit_polynomial(x, y, coeffs, DEGREE_OF_FIT);
      dpdt-= derivative_polynomial(coeffs,0);
    }
    point_to_evolve.change_rho(dpdt*time_step);
  }
}



void Grid::step_velocities(){
  double dvdt, dkdx, dvdx, dqdx, *yk, *yv, *yq, *coeffsk, *coeffsv, *coeffsq,*vel_of_point, *vel_below, *vel_above;
  int i, j,k, *coords_of_point, *coords_above, *coords_below, *x;
  Point point_to_evolve, other_point;

  x       = new int[3]();
  yk      = new double[3]();
  yv      = new double[3]();
  yq      = new double[3]();
  coeffsk = new double[3]();
  coeffsv = new double[3]();
  coeffsq = new double[3]();
  coords_of_point = new int[config_dimension]();
  coords_above    = new int[config_dimension]();
  coords_below    = new int[config_dimension]();
  vel_above       = new double[config_dimension]();
  vel_below       = new double[config_dimension]();

  for(i=0; i<total_points; i++){
        point_to_evolve = points[j];
        point_to_evolve.get_coordinates(coords_of_point, config_dimension);
        point_to_evolve.get_velocities(vel_of_point, config_dimension);

        for(j=0; j<config_dimension; j++){

          memcpy(coords_above, coords_of_point, sizeof(coords_of_point));
          memcpy(coords_below, coords_of_point, sizeof(coords_of_point));
          memcpy(vel_above, vel_of_point, sizeof(vel_of_point));
          memcpy(vel_below, vel_of_point, sizeof(vel_of_point));
          dvdt=0, dkdx=0, dqdx=0, dvdx=0;

          std::fill_n(x, 3, 0);
          std::fill_n(yv, 3, point_to_evolve.get_V_old());
          std::fill_n(yq, 3, point_to_evolve.get_Q_old());
          std::fill_n(coeffsk, 3, 0);
          std::fill_n(coeffsq, 3, 0);
          std::fill_n(coeffsv, 3, 0);

          if(coords_of_point[j] < grid_length - 1){
            x[2] =  1;
            coords_above[j]+=1;
            other_point =  points[coordinates_to_index(coords_above, config_dimension,grid_length)];
            other_point.get_velocities(vel_above, config_dimension);
            yv[2] = other_point.get_V_old();
            yq[2] = other_point.get_Q_old();
            coords_above[j]-=1;
          }
          if(coords_of_point[j] > 0){
            x[0] =  -1;
            coords_below[j]-=1;
            other_point =  points[coordinates_to_index(coords_below, config_dimension,grid_length)];
            other_point.get_velocities(vel_below, config_dimension);
            yv[0] = other_point.get_V_old();
            yq[0] = other_point.get_Q_old();
            coords_below[j]+=1;
          }

          fit_polynomial(x, yv, coeffsv, DEGREE_OF_FIT);
          fit_polynomial(x, yq, coeffsq, DEGREE_OF_FIT);
          dvdx = derivative_polynomial(coeffsv, 0);
          dqdx = derivative_polynomial(coeffsq, 0);

          for(k=0; k<config_dimension; k++){
            yk[0] = pow(vel_below[k],2);
            yk[1] = pow(vel_of_point[k],2);
            yk[2] = pow(vel_above[k],2);
            fit_polynomial(x, yk, coeffsk, DEGREE_OF_FIT);
            dkdx += (mass[k]/2) * derivative_polynomial(coeffsk, 0);
          }

          dvdt = -(dvdx + dqdx + dkdx);
          point_to_evolve.change_velocity_i(j, dvdt*time_step);
        }
      }
}



void Grid::update_Q(){
  double dvdt, Q_new,Q_j, *y, *coeffs;
  int i, j,k, *coords_of_point, *coords_above, *coords_below, *x;
  Point point_to_evolve, other_point;

  x       = new int[3]();
  y       = new double[3]();
  coeffs  = new double[3]();
  coords_of_point = new int[config_dimension]();
  coords_above    = new int[config_dimension]();
  coords_below    = new int[config_dimension]();

  for(i=0; i<total_points; i++){
        point_to_evolve = points[j];
        point_to_evolve.get_coordinates(coords_of_point, config_dimension);

        for(j=0; j<config_dimension; j++){

          memcpy(coords_above, coords_of_point, sizeof(coords_of_point));
          memcpy(coords_below, coords_of_point, sizeof(coords_of_point));
          Q_j=0;

          std::fill_n(x, 3, 0);
          std::fill_n(y, 3, sqrt(point_to_evolve.get_rho_old()));
          std::fill_n(coeffs, 3, 0);

          if(coords_of_point[j] < grid_length - 1){
            x[2] =  1;
            coords_above[j]+=1;
            other_point =  points[coordinates_to_index(coords_above, config_dimension,grid_length)];
            y[2] = sqrt(other_point.get_rho_old());
            coords_above[j]-=1;
          }
          if(coords_of_point[j] > 0){
            x[0] =  -1;
            coords_below[j]-=1;
            other_point =  points[coordinates_to_index(coords_below, config_dimension,grid_length)];
            y[0] = sqrt(other_point.get_rho_old());
            coords_below[j]+=1;
          }

          fit_polynomial(x, y, coeffs, DEGREE_OF_FIT);
          Q_j = second_derivative_polynomial(coeffs, 0);
          Q_new += Q_j/mass[j];
        }
       point_to_evolve.change_Q(Q_new*H_BAR_SQUARED/(2*sqrt(point_to_evolve.get_rho_old())));
       }
}








void Grid::init_grid(){
  int i,j;
  printf("      ...Loading Parameters...\n");
  get_data();
  printf("      ...Parameters Successfully Loaded...\n      ...Displaying Simulation Parameters...\n      .............................................................................\n");
  printf("      || spatial_dimension  %2i || num_particles %2i || configuration_dimension %3i ||\n", spatial_dimension, num_particles, config_dimension);
  printf("      || grid_length      %4i || total_time  %4.2f || time_step              %4.2f ||\n", grid_length, total_time, time_step);
  printf("      || total_points     %4i || mass [", (grid_length*grid_length));
  for(i = 0; i<config_dimension; i+= spatial_dimension) printf("%3.2f, ", mass[i]);
  printf("]\n      .............................................................................\n      ...Initializing Points...");

  set_points();
  printf("      ...Points Successfully Initialized...\n");
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
			}}
	input.close();
}



void Grid::set_points(){
  int i;
  total_points = pow(grid_length, config_dimension);
  points = new Point[total_points]();
  printf("\n");
  for(i = 0; i<total_points; i++){
    points[i].init_point(config_dimension, mass, grid_length, num_particles, spatial_dimension, i);
    print_load_bar(i*1.0/total_points);
  }
  print_load_bar(1.0);
  printf("\n");
}



void Grid::print_evolution_info(){
  Point point_to_print = points[total_points/2+20];
  int i;
  double *vel;
  vel = new double[config_dimension]();
  point_to_print.get_velocities(vel, config_dimension);

  printf("     Time  %7.3f || P %10.3f || Q %10.3f || V %7.3f || ", time, point_to_print.get_rho(),point_to_print.get_Q() ,point_to_print.get_V());
  printf("velocity: [");
  for(i = 0; i<config_dimension; i++) printf("%6.2f, ", vel[i]);
  printf("]\n");

}


void Grid::print_load_bar(double progress){
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
