#include "point.h"
#include "math_functions.h"
#include <fstream>
#include <string.h>







void Point::init_point(int config_dimension, double *mass, int grid_length, int num_particles, int spatial_dimension, int index){
  int i, remainder;
  double time = 0;
  coordinates = new int[config_dimension]();
  velocities = new double[config_dimension]();
  velocities_old = new double[config_dimension]();
  index_to_coordinates(coordinates, index,config_dimension, grid_length);
  calc_velocities(coordinates, 0.0, velocities, config_dimension);
  memcpy(velocities_old, velocities, sizeof(velocities));
  Q = calc_quantum_potential(mass, coordinates, 0.0, config_dimension);
  Q_old = Q;
  V = calc_potential(coordinates, num_particles, spatial_dimension);
  V_old = V;
  rho = calc_rho(coordinates, 0.0, config_dimension);
  rho_old = rho;
}



void Point::update_old_values(){
  rho_old = rho;
  Q_old =Q;
  V_old =V;
  memcpy(velocities_old, velocities, sizeof(velocities));
}
