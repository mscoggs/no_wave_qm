#ifndef __POINT_H_INCLUDED__
#define __POINT_H_INCLUDED__

class Point{

public:
  void init_point(int config_dimension,double *mass, int grid_length, int num_particles, int spatial_dimension, int index);

  int* get_coordinates(){return coordinates;}
  int get_coordinate_i(int i){return coordinates[i];}
  double get_velocity_i(int i){return velocities[i];}
  double get_rho(){return rho;}
  double get_V(){return V;}
  double get_Q(){return Q;}
private:
  double rho, rho_old, Q, Q_old, V, V_old, *velocities, *velocities_old;
  int *coordinates;
};

#endif
