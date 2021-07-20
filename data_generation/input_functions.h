#ifndef __INPUT_FUNCTIONS_H_INCLUDED__
#define __INPUT_FUNCTIONS_H_INCLUDED__


struct Psi
{
     double real;
     double imaginary;
};

Psi calc_psi(double x, double y, int psi_function,double coord_to_distance, double *mass, bool velocity_perturbation, double v_perturb);

double calc_rho_vel_initial(int *coordinates,double *velocities, int config_dimension, int grid_length, double *mass, int psi_function, double coord_to_distance, bool velocity_perturbation, double v_perturb);
/**
    calculating the potential function

    @param coordinates the coordinates of the point
    @param num_particles number of particles in our system
    @param spatial_dimension
    @param potential_function specifying which function to use in a switch statement
*/
double calc_potential(double time, int *coordinates, int num_particles, int spatial_dimension, int potential_function, int grid_length, double *mass);



#endif
