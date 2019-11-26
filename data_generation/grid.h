#ifndef __GRID_H_INCLUDED__
#define __GRID_H_INCLUDED__

#include "point.h"

class Grid{

public:
  void init_grid();
  void print_load_bar(double progress);
  void evolve();
  void step_rho();
  void step_velocities();
  void calc_Q();

private:
  void print_evolution_info();
  void get_data();
  void set_points();
  int spatial_dimension, num_particles, config_dimension, grid_length, total_points;
  double *mass, time, time_step, total_time;
  Point *points;
};

#endif
