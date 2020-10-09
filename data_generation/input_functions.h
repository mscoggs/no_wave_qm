#ifndef __INPUT_FUNCTIONS_H_INCLUDED__
#define __INPUT_FUNCTIONS_H_INCLUDED__



/**
    setting the density during init_grid

    @param coordinates the coordinates of the point
    @param config_dimension the dimension of the grid, num_particles*spatial_dimension
    @param grid_length the length of the grid in each dimension
    @param rho_function specifying which function to use in a switch statement
*/
double calc_rho_initial(int *coordinates, int config_dimension, int grid_length, int rho_function, double *mass);



/**
    setting the velocity during init_grid

    @param coordinates the coordinates of the point
    @param velocities an array acting as the velocity vector for this point
    @param config_dimension the dimension of the grid, num_particles*spatial_dimension
    @param vel_function specifying which function to use in a switch statement
*/
void calc_velocities_initial(int *coordinates, double *velocities, int config_dimension, int vel_function,  int grid_length, double *mass);


/**
    calculating the potential function

    @param coordinates the coordinates of the point
    @param num_particles number of particles in our system
    @param spatial_dimension
    @param v_function specifying which function to use in a switch statement
*/
double calc_potential(int *coordinates, int num_particles, int spatial_dimension, int v_function, int grid_length, double *mass);



#endif
