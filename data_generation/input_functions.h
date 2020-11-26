#ifndef __INPUT_FUNCTIONS_H_INCLUDED__
#define __INPUT_FUNCTIONS_H_INCLUDED__


struct Psi
{
     double real;
     double imaginary;
};

Psi calc_psi(double x, double y, int psi_function);

double calc_rho_vel_initial(int *coordinates,double *velocities, int config_dimension, int grid_length, double *mass, int psi_function);
/**
    calculating the potential function

    @param coordinates the coordinates of the point
    @param num_particles number of particles in our system
    @param spatial_dimension
    @param v_function specifying which function to use in a switch statement
*/
double calc_potential(int *coordinates, int num_particles, int spatial_dimension, int v_function, int grid_length, double *mass);



#endif
