#include <fstream>
#include <string.h>


#include "point.h"
#include "math_functions.h"
#include "input_functions.h"




void Point::init_point(int config_dimension, double *mass, int grid_length, int num_particles, int spatial_dimension, int index, int psi_function, int v_function){
  int i, remainder;


  coordinates = new int[config_dimension]();
  velocities = new double[config_dimension]();
  velocities_old = new double[config_dimension]();
  velocities_oldest = new double[config_dimension]();
  index_to_coordinates(coordinates, index,config_dimension, grid_length);

  rho = calc_rho_vel_initial(coordinates, velocities, config_dimension, grid_length, mass,psi_function);
  memcpy(velocities_old, velocities, sizeof(double)*config_dimension);
  memcpy(velocities_oldest, velocities, sizeof(double)*config_dimension);
  V = calc_potential(coordinates, num_particles, spatial_dimension, v_function, grid_length, mass);
  rho_oldest = rho;
  rho_old = rho;
}




void Point::update_old_values(int config_dimension){

  rho_oldest = rho_old;
  rho_old = rho;
  memcpy(velocities_oldest, velocities_old, sizeof(double)*config_dimension);
  memcpy(velocities_old, velocities, sizeof(double)*config_dimension);
}
