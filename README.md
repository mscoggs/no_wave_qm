# no_wave_qm
Simulating the evolution according to a hamilton-jacobi formulation of QM which replaces the wave with a configuration space density and equations of motion. Trajectory tracking using a 4-th order Runge Kutta technique.



# data_generation

## REQUIREMENTS
- c++
- g++ compiler

## USAGE
Open up input_functions.cpp and add your case to calc_psi(), which shoud be a quantum state which determines initial conditions. The program will numerically calculate all of the intitial conditions given this state. Add your desired potential function to calc_potential().  Go to the inputs file and specifiy the simulation parameters, along with trajectories that you wish to track. After that:
```bash  
make run
```


# data_analysis

## REQUIREMENTS
 - anaconda
     - jupyter notebook
     - python
     - glob 
     - matplotlib 
     - pandas 
 
 ## USAGE
 ```bash  
jupyter notebook analyze_data.ipynb
```
then follow the markdown in the notebook



# Paper (in prep)
