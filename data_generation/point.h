#ifndef __POINT_H_INCLUDED__
#define __POINT_H_INCLUDED__

class Point{

public:
  void init_point(int config_dimension,double *mass, int grid_length, int num_particles, int spatial_dimension, int index);

  double get_rho(){return rho;}
  double get_rho_old(){return rho_old;}
  double get_V(){return V;}
  double get_V_old(){return V_old;}
  double get_Q(){return Q;}
  double get_Q_old(){return Q_old;}
  double get_velocity_old_i(int i){return velocities[i];}
  void get_coordinates(int *coords, int config_dimension){
    int i;
    for(i=0;i<config_dimension;i++) coords[i] = coordinates[i];}
  void get_velocities(double *vel, int config_dimension){
    int i;
    for(i=0;i<config_dimension;i++) vel[i] = velocities[i];}
  void change_velocity_i(int i, double change){velocities[i] += change;}
  void change_rho(double change){rho += change;}
  void change_Q(double change){Q += change;}

private:
  double rho, rho_old, Q, Q_old, V, V_old, *velocities, *velocities_old;
  int *coordinates;
};

#endif
