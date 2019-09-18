#ifndef __POINT_H_INCLUDED__
#define __POINT_H_INCLUDED__

class Point{

public:
  void init_point(int config_dimension, int grid_length, int index);
  void update_velocities(int config_dimension, int num_particles, double *mass, double time_step, double time);
  void update_rho();

  int get_coordinates_i(int i){return coordinates[i];}
  int get_rho(){return rho;}
private:
  double rho, Q, V, *velocities;
  int *coordinates;
};

#endif
