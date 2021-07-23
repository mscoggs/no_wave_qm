#ifndef __GRID_H_INCLUDED__
#define __GRID_H_INCLUDED__


#include <string>


#include "point.h"




class Grid{


public:
  /**
      initializing the grid, printing the grid parameters, and initializing all of the points (set_points())
  */
  void init_grid();



  /**
      Evolving the system, iteratively stepping rho/velocities/Q for each point
  */
  void evolve();



  /**
      Stepping the rho for each point
  */
  void step_rho();



  /**
      Stepping the velocity for each point
  */
  void step_velocities();



  /**
      calculating the 'quantum potential'
  */
  void calc_Q();



  void step_V();


  /**
      update the values after each step in the evolution
  */
  void update_grid();


private:
  /**
      Initializing each point's values according the input functions
  */
  void set_points();



  /**
      printing all of the parameters for the evolution
  */
  void print_evolution_info();



  /**
      grabbing all the data from input.txt file located in ../inputs/
  */
  void get_data();



  /**
      writing the data for one time step to a file titled the current time which is locoated in ../data/dir_  where dir is named by the inputs in input.txt
  */
  void write_data();


  /**
      calculating the loop integral [int_0^2\pi  v \cdot dl] for a fixed radius and center, hardset in the function.
  */
  void calc_loop();
  void get_loop_velocity(double *position, double *vel);


  /**
      Stepping the velocity by grabbing the velocity of the 'nearest neighbor' -- implementing a runge-kutta method soon
  */
  void step_trajectories();
  void step_trajectory_runge_kutta(double *position);
  void get_runge_kutta_velocity(double t, double *position_temp, double *velocity_of_trajectory);




  int spatial_dimension, num_particles, config_dimension, grid_length, total_points, degree_of_fit, num_trajectories, num_neighbors, num_loops, vel_function, rho_function, potential_function, psi_function, velocity_perturbation, density_perturbation;
  double *mass, *trajectories, time, time_step, total_time, coord_to_distance, point_offset, time_step_save,v_perturb, *k1, *k2, *k3, *k4, *loop_vals, *loop_radii;
  Point *points;
  bool save;
  std::string dir, write_folder;
};

#endif
