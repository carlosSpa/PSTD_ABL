# PSTD_ABL

Introduction
------------
Source code files and documentation of Calibration of Absorbing Boundary Layers (ABL) for Geoacoustic Wave Modeling in Pseudo-Spectral Time-Domain Methods.


Code compilation: 
----------------
g++ -lfftw3 -lfftw3f -fopenmp -O3 *.c -o executable

./executable

Dependencies:
-------------
* fftw3
* openmp


**Test 1 (CALIBRATION)**
---------

 AUTHOR: Carlos Spa
 
 DATE: 22-06-2022
 
 FOLDERS: Damped, PML, Sponge

 COMMENTS: Each folder contains some C++ codes comprising a main code and some libraries:
  
 MAIN CODE: main.c 
 
 INTERNAL LIBRARIES: "analisis.h", "pstd_opt.h"
 
 EXTERNAL LIBRARIES: <fftw3.h>, <omp.h>
 
 COMPILATION: g++ -lfftw3 -fopenmp -O3 *.c -o executable
 
 EXECUTION: ./executable
 
 OUTPUT: binary file "Energy.dat"
 
COMMENTS: This code calculates the energy of 600 different simulations defined on an homogeneous cube of 500x500x500 m³ and a propagation velocity of c=2000m/s.
A point source corresponding to a ricker wavelet of f=10Hz located at the center of the cube.
In all the cases, we also fix the uniform spatial sampling Dx=40m, the temporal step Dt=0.002s and the total simulation time to 2s.
We solve this problem by using use the numerical method provided in this folder. The code performs 600 simulations by varying the tuple (Nabl,Absorbing parameter) in 30x20 different cases according to Table 1 in the paper. The output of these simulations are stored in "Energy.dat".

DATA ANALYSIS CODE: energy.m

EXECUTION: octave or matlab

OUTPUT: Fig 2 and Table 2. 

COMMENTS: By running the script "energy.m", the results of Fig. 2 and Table 2 can be reproduced. 



**Test 2 (GEOPHYSICAL IMAGING)**
--------- 

 AUTHOR: Carlos Spa

 DATE: 22-06-2022
 
 FOLDERS: damped (D), sponge (S), euler (E)
 
 COMMENTS: A single C main code is available at each folder along with some libraries: 
 MAIN CODE: main.c 
 
 INTERNAL LIBRARIES: "analisis.h", "pstd_opt.h", "gestor_sim.h"
 
 EXTERNAL LIBRARIES: <fftw3.h>, <omp.h>
 
 COMPILATION: g++ -lfftw3 -lfftw3f -fopenmp -O3 *.c -o executable
 
 EXECUTION: ./executable
 
 OUTPUT: binary file 2x7x3 Files: "P_type_Acc_abl_method.dat" where type=F,B  abl=1,2,3,4,5,6,7 and method=D,S,E.
 The results are stored at each corresponding folder.
 
COMMENTS: This code calculates the forward/backward  (x2) simulations for each different method (D,S,E) for all cases listed in table 2. The simulation domain is an homogeneous cube of 4000x4000x4000 m³ and a propagation velocity of c=2000m/s.
In all the cases, we also fix the spatial sampling Dx=40m, the temporal step Dt=0.002s and the total simulation time to 4s.
For the Forward simulation, a point source ricker wavelet of f=10Hz is located at the grid position (Nx/2,Ny/2,4).
For the Backward simulation, three point source ricker wavelets of f=10Hz are located at the grid position (Nx/2,Ny/2,4).
The output data of each ABL method is stored at the corresponding folder. 
Each method performs 7x2 simulations by varying the tuple (Nabl,Absorbing parameter) according to Table 2 in the paper. 

DATA ANALYSIS CODE: analisis_METHOD.m where METHOD=D,S,E

LIBRARIES: PE_misfit.m, ricker_wavelet_tis0.m

EXECUTION: octave or matlab

OUTPUT: Table 3. 

COMMENTS: Each folder contains a matlab code that generates the results of Table 3
 

**Test 3 (EAGE SEG-SALT)**
-------------

 AUTHOR: Carlos Spa

 DATE: 22-06-2022
 
 FOLDERS: the folder "Mesh" should be downloaded from: https://doi.org/10.5281/zenodo.7821703
 
 COMMENTS: It should be unzip at the same folder where the main C code is stored. 
 The "Mesh" folder contains three sub-folders (damped,sponge,euler) with different binary files used as input files by the main code
 
 MAIN CODE: main.c 
 
 INTERNAL LIBRARIES: "analisis.h", "pstd_opt.h", "gestor_sim.h"
 
 EXTERNAL LIBRARIES: <fftw3.h>, <omp.h>
 
 COMPILATION: g++ -lfftw3 -lfftw3f -fopenmp -O3 *.c -o executable
 
 EXECUTION: ./executable
 
 INPUT: A binary mesh is required and is located at ./Mesh/METHOD/3d.outAcc_NUMBER.data, where METHOD=damped,sponge,euler and NUMBER=1,2,3,4,5,6,7
 
 OUTPUT: 7x3 binary files EnergyAcc_NUMBER_METHOD where NUMBER=1,2,3,4,5,6,7 and METHOD=D,S,E
 
COMMENTS: This code performs the simulation for a unit impulse source signal in the SEG SALT EAGE model for each different method (D,S,E) and for all cases listed in table 2. The simulation domain is a cube of 251x251x122 nodes, where the SEG SALT EAGE model is embedded in a homogenous medium of propagation velocity of c=2000m/s.
In all the cases, we also fix the spatial sampling Dx=40m, the temporal step Dt=0.002s and the total simulation time to 4s.
For the Forward simulation, a ricker point source wavelet of f=10Hz is emitted at the grid position (Nx/2,Ny/2,4).
The output data of each method is stored at the corresponding folder. 
Each method performs 7 simulations by varying the tuple (Nabl,Absorbing parameter) according to Table 2 in the paper, in adition to one reference simulation perfomed by the PML method.

DATA ANALYSIS CODE: energy.m

EXECUTION: octave or matlab

OUTPUT: Table 4. 


 
 
 
