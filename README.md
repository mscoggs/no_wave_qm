# Hamilton-Jacobi No-Wave Quantum Simulator
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



# Paper 
```
@ARTICLE{2023arXiv230304959R,
       author = {{Roser}, Philipp and {Scoggins}, Matthew T.},
        title = "{Non-Quantum Behaviors of Configuration-Space Density Formulations of quantum mechanics}",
      journal = {arXiv e-prints},
     keywords = {Quantum Physics},
         year = 2023,
        month = mar,
          eid = {arXiv:2303.04959},
        pages = {arXiv:2303.04959},
          doi = {10.48550/arXiv.2303.04959},
archivePrefix = {arXiv},
       eprint = {2303.04959},
 primaryClass = {quant-ph},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2023arXiv230304959R},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
```




