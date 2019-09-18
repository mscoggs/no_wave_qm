#include "point.h"
#include "math_functions.h"

void Point::init_point(int config_dimension, int grid_length, int index){
  int i, remainder;
  double time = 0;
  coordinates = new int[config_dimension]();
  velocities = new double[config_dimension]();

  for(i = config_dimension-1; i>=0;i--){
    remainder = index%grid_length;
    coordinates[i] = remainder;
    index = (index - remainder)/grid_length;
  }


  Q = calc_quantum_potential();
  V = calc_potential();
  rho = calc_rho();


  for(i = 0; i<config_dimension;i++) velocities[i] = calc_velocity_i(coordinates, 0.0, i);

}




void Point::update_velocities(int config_dimension, int num_particles, double *mass, double time_step, double time){
  double dvdt = 0, *velocities_temp;
  int i,j;
  velocities_temp = new double[config_dimension]();

  for(i=0; i<config_dimension; i++){
      for(j=0; j<config_dimension; j++){
        //dvdt += mass[j]*velocities[j]*derivative_i(calc_velocity_i, 0, i);
      //dvdt = -(dvdt + derivative_i(calc_potential,0, i) + derivative_i(calc_quantum_potential, 0,i));
      velocities_temp[i] = velocities[i] + dvdt * time_step;
      }
    }
  for(i=0; i<config_dimension; i++) velocities[i] = velocities_temp[i];
}

void Point::update_rho(){
  //p_new = p_old +dpdt
}
