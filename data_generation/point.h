#ifndef __POINT_H_INCLUDED__
#define __POINT_H_INCLUDED__


#include "input_functions.h"
class Point{


public:


  /**
      initializing the points

      @param config_dimension the dimension of the grid, num_particles*spatial_dimension
      @param mass the mass of each particle
      @param grid_length the length of the grid in each dimension
      @param num_particles the number of particles in our system
      @param spatial_dimension the spatial config_dimension
      @param index the index-th point in the grid, used to get the coordinates of the point

  */
  void init_point(int config_dimension,double *mass, int grid_length, int num_particles, int spatial_dimension, int index, int psi_function, int v_function);



  double get_rho(){
    return rho;
    // if(rho>=0)return rho;
    // else return 0;
  }
  double get_rho_old(){
    return rho_old;
    // if(rho_old>=0)return rho_old;
    // else return 0;
  }
  double get_rho_oldest(){
    return rho_oldest;
    // if(rho_old>=0)return rho_old;
    // else return 0;
  }
  void get_velocities(double *vel, int config_dimension){
    int i;
    for(i=0;i<config_dimension;i++) vel[i] = velocities[i];
  }
  void get_velocities_old(double *vel, int config_dimension){
    int i;
    for(i=0;i<config_dimension;i++) vel[i] = velocities_old[i];
  }
  void get_velocities_oldest(double *vel, int config_dimension){
    int i;
    for(i=0;i<config_dimension;i++) vel[i] = velocities_oldest[i];
  }
  double get_velocity_i(int i){return velocities[i];}
  double get_velocity_old_i(int i){return velocities_old[i];}
  double get_velocity_oldest_i(int i){return velocities_oldest[i];}

  double get_V(){return V;}
  double get_Q(){return Q;}


  void get_coordinates(int *coords, int config_dimension){
    int i;
    for(i=0;i<config_dimension;i++) coords[i] = coordinates[i];
  }
  int get_coordinates_i(int i){return coordinates[i];}


  void change_velocity_i(int i, double change){velocities[i] += change;}


  void change_rho(double change){
    rho += change;
    //if(rho < 0) rho -=change;
  }


  void set_Q(double val){Q = val;}


  void update_old_values(int config_dimension);


private:

  double rho, rho_old, rho_oldest, Q,V, *velocities, *velocities_old, *velocities_oldest;
  int *coordinates;
};

#endif
