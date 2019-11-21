#include "point.h"
#include "math_functions.h"
#include <fstream>







void Point::init_point(int config_dimension, double *mass, int grid_length, int num_particles, int spatial_dimension, int index){
  int i, remainder;
  double time = 0;
  coordinates = new int[config_dimension]();
  velocities = new double[config_dimension]();

  index_to_coordinates(coordinates, index, grid_length);
  calc_velocities(coordinates, 0.0, velocities, config_dimension);
  calc_velocities(coordinates, 0.0, velocities_old, config_dimension);
  Q = calc_quantum_potential(mass, coordinates, 0.0, config_dimension);
  Q_old = Q;
  V = calc_potential(coordinates, num_particles, spatial_dimension);
  V_old = V;
  rho = calc_rho(coordinates, 0.0, config_dimension);
  rho_old = rho;
}
