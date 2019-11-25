#ifndef __GRID_H_INCLUDED__
#define __GRID_H_INCLUDED__

#include "point.h"

class Grid{

public:
  void init_grid();
  void print_evolution_info();
  void print_load_bar(double progress);

  void evolve();
  void step_rho();
  void step_velocities();
  void update_Q();

  int get_spatial_dimension(){return spatial_dimension;}
  int get_grid_length(){return grid_length;}
  int get_num_particles(){return num_particles;}
  int get_config_dimension(){return config_dimension;}
  int get_total_points(){return total_points;}
  double get_mass_i(int i){return mass[i];}
  Point get_point_i(int i){return points[i];}

private:
  void get_data();
  void set_points();
  int spatial_dimension, num_particles, config_dimension, grid_length, total_points;
  double *mass, time, time_step, total_time;
  Point *points;
};

#endif
