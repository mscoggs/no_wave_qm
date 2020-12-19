#include <fstream>
#include <iostream>
#include <string>
#include <algorithm>
#include <complex.h>
#include <math.h>
#include <string.h>
#include <algorithm>
#include <sys/stat.h>
#include <sys/types.h>


#include "grid.h"
#include "point.h"
#include "math_functions.h"
#include "input_functions.h"




void Grid::init_grid(){
  int i,j;


  printf("      ...Loading Parameters...\n");
  get_data();
  printf("      ...Parameters Successfully Loaded...\n      ...Displaying Simulation Parameters...\n      .............................................................................\n");
  printf("      || spatial_dimension  %2i || num_particles %2i || configuration_dimension %3i ||\n", spatial_dimension, num_particles, config_dimension);
  printf("      || grid_length    %3i || total_time  %4.2f || time_step              %4.2f ||\n", grid_length, total_time, time_step);
  printf("      || total_points     %4i ||\n      || mass [", total_points);
  for(i = 0; i<config_dimension; i+= spatial_dimension) printf("%3.2f, ", mass[i]);
  printf("]\n      || initial_trajectories  [");
  for(i = 0; i<num_trajectories; i++){
    printf("[");
    for(j = 0; j<config_dimension; j++) printf("%3.2f, ", trajectories[i*config_dimension + j]);
    printf("],");
  }
  printf("]\n      .............................................................................\n      ...Initializing Points...\n");
  set_points();
  printf("      ...Points Successfully Initialized...\n");
}




void Grid::evolve(){
  int i = 0;

  while(time < total_time){
    print_evolution_info();
    if(i > 1 && i%2==0) step_trajectories();
    calc_Q();
    if(fmod(time+0.5*time_step, time_step_save) < time_step) write_data();
    //if(abs(time-7885)<0.5) write_data();
    step_rho();
    step_velocities();
    update_grid();
    time += time_step;
    i+=1;
    }
  write_data();
}




void Grid::step_rho(){
  double dpdt, *y, *x, *coeffs,offset,center = (grid_length-1.0)/2.0, x_, y_;
  int i, j, *coords_of_point, *coords_above, *coords_below;
  Point other_point;


  x   = new double[3]();
  y   = new double[3]();
  coeffs = new double[3]();
  coords_of_point = new int[config_dimension]();
  coords_above    = new int[config_dimension]();
  coords_below    = new int[config_dimension]();

  for(i=0; i<total_points; i++){

    dpdt = 0;
    points[i].get_coordinates(coords_of_point, config_dimension);

    x_ = static_cast<double>(coords_of_point[0])-center;
    y_ = static_cast<double>(coords_of_point[1])-center;
    if(sqrt(x_*x_ + y_*y_)  >= INF_WELL_RADIUS && write_folder=="inf_circular_well") continue;

    for(j=0; j<config_dimension; j++){

      std::fill_n(x, 3, point_offset);
      std::fill_n(y, 3, points[i].get_rho_old()*points[i].get_velocity_old_i(j));
      std::fill_n(coeffs, 3, 0.0);
      memcpy(coords_above, coords_of_point, config_dimension*sizeof(int));
      memcpy(coords_below, coords_of_point, config_dimension*sizeof(int));

      if(coords_above[j] < grid_length - 1){
        x[2] +=  coord_to_distance;
        coords_above[j]+=1;
        other_point =  points[coordinates_to_index(coords_above, config_dimension,grid_length)];
        y[2] =other_point.get_rho_old()*other_point.get_velocity_old_i(j);
      }

      if(coords_below[j] > 0) {
        x[0] -=  coord_to_distance;
        coords_below[j]-=1;
        other_point =  points[coordinates_to_index(coords_below, config_dimension,grid_length)];
        y[0] =other_point.get_rho_old()*other_point.get_velocity_old_i(j);
      }

      fit_polynomial(x, y, coeffs, degree_of_fit);
      dpdt-= nth_derivative_polynomial(coeffs, point_offset, degree_of_fit, 1);

    }
    points[i].change_rho(dpdt*time_step);
    if(INF_CHECK){
      if(isnan(points[i].get_rho())) printf("dpdt: %f\n\n\n", dpdt), exit(0);
    }
  }
  delete[] x, delete[] y, delete[] coeffs, delete[] coords_of_point, delete[] coords_below, delete[] coords_above;
}




void Grid::step_velocities(){
  double dvdt, dkdx, dvdx, dqdx, *x, *yk, *yv, *yq, *coeffsk, *coeffsv, *coeffsq,*vel_of_point, *vel_below, *vel_above, center = (grid_length-1.0)/2.0, x_,y_;
  int i, j,k, *coords_of_point, *coords_above, *coords_below;
  Point other_point;


  x       = new double[3]();
  yk      = new double[3]();
  yv      = new double[3]();
  yq      = new double[3]();
  coeffsk = new double[3]();
  coeffsv = new double[3]();
  coeffsq = new double[3]();
  coords_of_point = new int[config_dimension]();
  coords_above    = new int[config_dimension]();
  coords_below    = new int[config_dimension]();
  vel_of_point       = new double[config_dimension]();
  vel_above       = new double[config_dimension]();
  vel_below       = new double[config_dimension]();

  for(i=0; i<total_points; i++){

    points[i].get_coordinates(coords_of_point, config_dimension);
    points[i].get_velocities(vel_of_point, config_dimension);


    for(j=0; j<config_dimension; j++){

      memcpy(coords_above, coords_of_point, config_dimension*sizeof(int));
      memcpy(coords_below, coords_of_point, config_dimension*sizeof(int));
      memcpy(vel_above, vel_of_point, config_dimension*sizeof(double));
      memcpy(vel_below, vel_of_point, config_dimension*sizeof(double));
      dvdt=0, dkdx=0, dqdx=0, dvdx=0;
      std::fill_n(x, 3, point_offset);
      std::fill_n(yv, 3, points[i].get_V());
      std::fill_n(yq, 3, points[i].get_Q());
      std::fill_n(coeffsk, 3, 0.0);
      std::fill_n(coeffsq, 3, 0.0);
      std::fill_n(coeffsv, 3, 0.0);

      if(coords_of_point[j] < grid_length-1){
        x[2] +=  coord_to_distance;
        coords_above[j]+=1;
        other_point =  points[coordinates_to_index(coords_above, config_dimension,grid_length)];
        other_point.get_velocities_old(vel_above, config_dimension);
        yv[2] = other_point.get_V();
        yq[2] = other_point.get_Q();
        coords_above[j]-=1;
      }

      if(coords_of_point[j] > 0){
        x[0] -= coord_to_distance;
        coords_below[j]-=1;
        other_point =  points[coordinates_to_index(coords_below, config_dimension,grid_length)];
        other_point.get_velocities_old(vel_below, config_dimension);
        yv[0] = other_point.get_V();
        yq[0] = other_point.get_Q();
        coords_below[j]+=1;
      }

      fit_polynomial(x, yv, coeffsv, degree_of_fit);
      fit_polynomial(x, yq, coeffsq, degree_of_fit);

      dvdx = nth_derivative_polynomial(coeffsv, point_offset, degree_of_fit, 1);
      dqdx = nth_derivative_polynomial(coeffsq, point_offset, degree_of_fit, 1);

       for(k=0; k<config_dimension; k++){
          yk[0] = pow(abs(vel_below[k]),2), yk[1] = pow(abs(vel_of_point[k]),2), yk[2] = pow(abs(vel_above[k]),2);
          fit_polynomial(x, yk, coeffsk, degree_of_fit);
          dkdx += (mass[k]/2) * nth_derivative_polynomial(coeffsk, point_offset, degree_of_fit, 1);
        }
        dvdt = -(dvdx + dqdx + dkdx)/mass[j];

        points[i].change_velocity_i(j, dvdt*time_step);
      }
    }
    delete[] x, delete[] yk, delete[] yv, delete[] yq, delete[] coeffsk, delete[] coeffsv, delete[] coeffsq, delete[] coords_of_point, delete[] coords_below, delete[] coords_above, delete[] vel_of_point, delete[] vel_above, delete[] vel_below;
}




void Grid::calc_Q(){
  double dvdt, Q_new,Q_j, *y, *x, *coeffs, center = (grid_length-1.0)/2.0, x_, y_;
  int i, j,k, *coords_of_point, *coords_above, *coords_below, edge_sign;
  Point point_above, point_below;


  x       = new double[3]();
  y       = new double[3]();
  coeffs  = new double[3]();
  coords_of_point = new int[config_dimension]();
  coords_above    = new int[config_dimension]();
  coords_below    = new int[config_dimension]();

  for(i=0; i<total_points; i++){

        points[i].get_coordinates(coords_of_point, config_dimension);

        Q_new = 0;
        for(j=0; j<config_dimension; j++){

          memcpy(coords_above, coords_of_point, config_dimension*sizeof(int));
          memcpy(coords_below, coords_of_point, config_dimension*sizeof(int));
          Q_j=0;
          std::fill_n(x, 3, point_offset);
          std::fill_n(y, 3, sqrt(points[i].get_rho_old()));
          std::fill_n(coeffs, 3, 0.0);

          if((coords_of_point[j] == 0) || coords_of_point[j] == grid_length-1){ //this is an edge case, fit a polynomial centered around the point interior to this.
            if(coords_of_point[j] == 0){
              x[2] += 2*coord_to_distance;
              coords_above[j] += 2;
              point_above =  points[coordinates_to_index(coords_above, config_dimension,grid_length)];
              y[2] = sqrt(point_above.get_rho_old());

              x[1] += coord_to_distance;
              coords_below[j] += 1;
              point_below =  points[coordinates_to_index(coords_below, config_dimension,grid_length)];
              y[1] = sqrt(point_below.get_rho_old());

              edge_sign = 1;
            }
            else{
              x[1] -= coord_to_distance;
              coords_above[j] -= 1;
              point_above =  points[coordinates_to_index(coords_above, config_dimension,grid_length)];
              y[1] = sqrt(point_above.get_rho_old());

              x[0] -= 2*coord_to_distance;
              coords_below[j] -=2;
              point_below =  points[coordinates_to_index(coords_below, config_dimension,grid_length)];
              y[0] = sqrt(point_below.get_rho_old());


              edge_sign = -1;
            }

            fit_polynomial(x, y, coeffs, degree_of_fit);
            Q_j = nth_derivative_polynomial(coeffs, point_offset-edge_sign*coord_to_distance, degree_of_fit, 2);

          }
          else{
            x[2] +=  coord_to_distance;
            coords_above[j]+=1;
            point_above =  points[coordinates_to_index(coords_above, config_dimension,grid_length)];
            y[2] = sqrt(point_above.get_rho_old());

            x[0] -= coord_to_distance;
            coords_below[j]-=1;
            point_below =  points[coordinates_to_index(coords_below, config_dimension,grid_length)];
            y[0] = sqrt(point_below.get_rho_old());

            fit_polynomial(x, y, coeffs, degree_of_fit);

            Q_j = nth_derivative_polynomial(coeffs, point_offset, degree_of_fit, 2);

          }

          Q_new += Q_j/mass[j];
        }
        if(points[i].get_rho_old() == 0) points[i].set_Q(0);
        else points[i].set_Q(-1*Q_new*H_BAR_SQUARED/(2*sqrt(points[i].get_rho_old())));
       }
       delete[] x, delete[] y, delete[] coeffs, delete[] coords_of_point, delete[] coords_below, delete[] coords_above;
}





void Grid::update_grid(){
  int i;
  for(i = 0; i<total_points; i++) points[i].update_old_values(config_dimension);
}




void Grid::set_points(){
  int i;
  total_points = int_pow(grid_length, config_dimension);
  points = new Point[total_points]();
  for(i = 0; i<total_points; i++){
    points[i].init_point(config_dimension, mass, grid_length, num_particles, spatial_dimension, i, psi_function, v_function, coord_to_distance);
    print_load_bar(i*1.0/total_points);
  }
  calc_Q();
  print_load_bar(1.0);
  printf("\n");
}




void Grid::print_evolution_info(){
  Point point_to_print = points[total_points/2+2];
  int i;
  double *vel;
  vel = new double[config_dimension]();
  point_to_print.get_velocities(vel, config_dimension);

  printf("     Time %6.3f || P %9.8f || Q %9.8f || V %9.8f || ", time, point_to_print.get_rho(),point_to_print.get_Q() ,point_to_print.get_V());
  printf("velocity: [");
  for(i = 0; i<config_dimension; i++) printf("%9.8f, ", vel[i]);
  printf("]\n");
  delete[] vel;
}




void Grid::get_data(){


  std::ifstream input("../inputs/input.txt");
	if(input.is_open()){
		std::string line, delimiter = " ", token;
		size_t pos = 0;
		int i =0, j=0;
		while(getline(input, line)){

      if(line=="GRID_LENGTH"){
				getline(input, line);
        grid_length = stoi(line);
				continue;
			}
			if(line=="SPATIAL_DIMENSION"){
				getline(input, line);
				spatial_dimension = stoi(line);
				continue;
			}
			if(line=="NUMBER_OF_PARTICLES"){
				getline(input, line);
				num_particles = stoi(line);
        config_dimension = num_particles*spatial_dimension;
        mass = new double[config_dimension]();
				continue;
			}
			if(line=="MASS"){
				getline(input, line);
				for(i=0; i<num_particles; i++){
						pos = line.find(delimiter);
						token = line.substr(0, pos);
            for(j=0; j<spatial_dimension; j++) mass[i*spatial_dimension + j] = stod(token);
						line.erase(0, pos + delimiter.length());
				}
      }
      if(line=="TOTAL_TIME"){
        getline(input, line);
        total_time = stod(line);
        time = 0;
        continue;
      }
      if(line=="TIME_STEP"){
        getline(input, line);
        time_step = stod(line);
        continue;
      }
      if(line=="COORD_TO_DISTANCE"){
        getline(input, line);
        coord_to_distance = stod(line);
        point_offset = 0;//coord_to_distance+ 0.01;
        continue;
      }
      if(line=="DEGREE_OF_FIT"){
        getline(input, line);
        degree_of_fit = stoi(line);
        continue;
      }
      if(line=="DEGREE_OF_FIT"){
        getline(input, line);
        degree_of_fit = stoi(line);
        continue;
      }
      if(line=="DEGREE_OF_FIT"){
        getline(input, line);
        degree_of_fit = stoi(line);
        continue;
      }
      if(line=="NUMBER_OF_TRAJECTORIES"){
        getline(input, line);
        num_trajectories = stoi(line);
        k1 = new double[3*num_trajectories*config_dimension]();
        k2 = new double[3*num_trajectories*config_dimension]();
        k3 = new double[3*num_trajectories*config_dimension]();
        k4 = new double[3*num_trajectories*config_dimension]();
        num_neighbors = pow(2, config_dimension);
        continue;
      }
      if(line=="INITIAL_TRAJECTORIES"){
        getline(input, line);
        trajectories = new double[config_dimension*num_trajectories]();
        std::stringstream ss(line);
        for(i=0; i<config_dimension*num_trajectories; i++){
            pos = line.find(delimiter);
            token = line.substr(0, pos);
            trajectories[i] = stod(token);
            line.erase(0, pos + delimiter.length());
        }
        continue;
      }
      if(line=="INITIAL_VELOCITY_FUNCTION"){
        getline(input, line);
        vel_function = stoi(line);
        continue;
      }
      if(line=="INITIAL_DENSITY_FUNCTION"){
        getline(input, line);
        rho_function = stoi(line);
        continue;
      }
      if(line=="PSI_FUNCTION"){
        getline(input, line);
        psi_function = stoi(line);
        continue;
      }
      if(line=="POTENTIAL_FUNCTION"){
        getline(input, line);
        v_function = stoi(line);
        continue;
      }
      if(line=="WRITE_FOLDER"){
        getline(input, line);
        write_folder = line;
        continue;
      }
      if(line=="TIME_STEP_SAVE"){
        getline(input, line);
        time_step_save = stod(line);
        continue;
      }
      		}
  }
	input.close();
  total_points = int_pow(grid_length, config_dimension);
}




void Grid::write_data(){
  int i,j;
  std::string coord_string = "[", velocities_string = "[";
  std::string file_name_evolve = "time_" +std::to_string(time)+ ".txt";;
  std::string file_name_trajectories = "trajectories.txt";
  std::ofstream file;


  //creating the directory/wiping the trajectory file for the first call to this function
  if(time == 0){

    dir = "../data/" + write_folder;

    std::string cmd1 = "rm -rf ";
    std::string cmd2 = "mkdir ";
    int status = system((cmd1 + dir).c_str());
    status = system((cmd2 + dir).c_str());
    std::ifstream src("../inputs/input.txt", std::ios::binary);
    std::ofstream dest(dir+"/input.txt", std::ios::binary);
    dest << src.rdbuf();
    file.open(dir + "/" + file_name_trajectories);
    file << "time    ";
    for(i = 0; i<num_trajectories; i++) file << "position" <<  i << "   ";
    file << "\n";
    file.close();
  }

  //writing the evolution information for every point to a single file
  file.open(dir + "/" + file_name_evolve);
  file << "coordinates  index   rho      Q      V      velocities\n";
  for(i = 0; i<total_points; i++){
    for(j = 0; j<config_dimension; j++){
      coord_string += std::to_string(points[i].get_coordinates_i(j)) + ",";
      velocities_string += std::to_string(points[i].get_velocity_i(j)) + ",";
    }
    coord_string.append("]");
    velocities_string.append("]");
    file  << coord_string << "   " << i << "   " << points[i].get_rho() << "   " << points[i].get_Q() << "   " << points[i].get_V() << "   " << velocities_string << "\n";
    coord_string = "[", velocities_string = "[";
  }
  file.close();

  //writing the trajectories
  file.open(dir + "/" + file_name_trajectories, std::ios::app);
  file << time << "  ";
  for(i = 0; i<num_trajectories; i++){
    file << "[";
    for(j = 0; j<config_dimension; j++) file << trajectories[i*config_dimension + j] << ",";
    file << "]      ";
  }
  file << "\n";
  file.close();
}






void Grid::step_trajectories(){
  int i,j;
  double  *coords_of_trajectory;

  coords_of_trajectory = new double[config_dimension]();

  for(i=0; i<num_trajectories; i++){
    //making sure the coordinates aren't outside the boundaries of the grid, otherwise we skip them
    start:
    for(j=0; j<config_dimension; j++){
      if(trajectories[i*config_dimension + j] < 0 || trajectories[i*config_dimension + j] > grid_length-1){
        i++;
        goto start;
      }
    }

    for(j=0; j<config_dimension; j++) coords_of_trajectory[j] = trajectories[i*config_dimension + j];
    step_trajectory_runge_kutta(coords_of_trajectory);
    for(j=0; j<config_dimension; j++) trajectories[i*config_dimension + j] = coords_of_trajectory[j];


  }
  delete[] coords_of_trajectory;
}






void Grid::step_trajectory_runge_kutta(double *position){
  int i;
  double *dr1, *dr2, *dr3, *dr4, *velocities_at_trajectory, *position_temp;

  dr1 = new double[config_dimension]();
  dr2 = new double[config_dimension]();
  dr3 = new double[config_dimension]();
  dr4 = new double[config_dimension]();
  position_temp = new double[config_dimension]();
  velocities_at_trajectory = new double[config_dimension]();

  //k1 = hf(x_n, y_n)
  memcpy(position_temp, position, config_dimension*sizeof(double));
  get_runge_kutta_velocity(time-2*time_step, position_temp, velocities_at_trajectory);
  for(i=0; i<config_dimension; i++) dr1[i] = 2*time_step*velocities_at_trajectory[i];

  //k2 = hf(x_n+1, y_n+k1/2)
  memcpy(position_temp, position, config_dimension*sizeof(double));
  for(i=0; i<config_dimension; i++) position_temp[i] += dr1[i]/2;
  get_runge_kutta_velocity(time-time_step, position_temp, velocities_at_trajectory);
  for(i=0; i<config_dimension; i++) dr2[i] = 2*time_step*velocities_at_trajectory[i];

  //k3 = hf(x_n+1, y_n+k2/2)
  memcpy(position_temp, position, config_dimension*sizeof(double));
  for(i=0; i<config_dimension; i++) position_temp[i] += dr2[i]/2;
  get_runge_kutta_velocity(time-time_step, position_temp, velocities_at_trajectory);
  for(i=0; i<config_dimension; i++) dr3[i] = 2*time_step*velocities_at_trajectory[i];

  //k4 = hf(x_n+2, y_n+k3)
  memcpy(position_temp, position, config_dimension*sizeof(double));
  for(i=0; i<config_dimension; i++) position_temp[i] += dr3[i];
  get_runge_kutta_velocity(time, position_temp, velocities_at_trajectory);
  for(i=0; i<config_dimension; i++) dr3[i] = 2*time_step*velocities_at_trajectory[i];


  for(i=0; i<config_dimension; i++) position[i] += (dr1[i]+dr4[i])/6 + (dr2[i]+dr3[i])/3;

  delete[] dr1, delete[] dr2, delete[] dr3, delete[] dr4, delete[] position_temp, delete[] velocities_at_trajectory;
}




void Grid::get_runge_kutta_velocity(double t, double *position_temp, double *velocity_of_trajectory){
  int i,j,k, *neighbor_coords, n;
  Point point;
  double *neighbor_vel, *delta_n;


  neighbor_coords = new int[config_dimension*num_neighbors]();
  neighbor_vel = new double[config_dimension*num_neighbors]();
  delta_n = new double[spatial_dimension]();


  //getting all neighbor combinations using bit a bit operation
  for(i=0; i<num_neighbors; i++){
    for(j=0; j<config_dimension; j++){
      if(i&(1<<j)) neighbor_coords[i*config_dimension + j] = ceil(position_temp[j]);
      else neighbor_coords[i*config_dimension + j] = floor(position_temp[j]);
    }
  }

  //grabbing the velocities of all the neighbors
  for(j=0; j<num_neighbors; j++){
    point = points[coordinates_to_index(&neighbor_coords[j*config_dimension], config_dimension,grid_length)];
    if(t == time) point.get_velocities(&neighbor_vel[j*config_dimension], config_dimension);
    else if(t == time-time_step) point.get_velocities_old(&neighbor_vel[j*config_dimension], config_dimension);
    else if(t == time-2*time_step) point.get_velocities_oldest(&neighbor_vel[j*config_dimension], config_dimension);
    else printf("\n\n\n\nAKLSJDALKJSDKLAJSDKLAJSDLKJASLKDHALKSHJAS\n\n\n\n\n");
  }


  for(i=0; i<spatial_dimension; i++) delta_n[i] = position_temp[i] - floor(position_temp[i]);


  //doing a simple linear interpolation, first interpolating in the 'x' direction, then y, then z ... should work for an abitrary configuration_dimension
  n = num_neighbors;
  for(i=0; i<config_dimension; i++){
    for(j=0; j<n; j+=2){
      for(k=0; k<config_dimension;k++) neighbor_vel[j*config_dimension +k] += (neighbor_vel[(j+1)*config_dimension +k]-neighbor_vel[(j)*config_dimension +k]) *delta_n[i];
    }
    n = n/2;
    for(j=0; j<n; j++){
      for(k=0; k<config_dimension;k++) neighbor_vel[j*config_dimension +k] = neighbor_vel[j*config_dimension*2 +k];
    }
  }
  for(i=0; i<config_dimension; i++) velocity_of_trajectory[i] = neighbor_vel[i];


  delete[] neighbor_coords, delete[] neighbor_vel, delete[] delta_n;
}
