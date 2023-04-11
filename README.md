# PSTD_ABL

Introduction
------------
Source code files and documentation of Absorbing Boundary Layer (ABL) techniques in Pseudo-Spectral Time-Domain methods


Code compilation: 
----------------
g++ -lfftw3 -fopenmp -O3 *.c -o executable

./executable

Dependences:
-------------
* fftw3
* openmp


**Test 1**
---------

 AUTHOR: Carlos Spa
 
 DATE: 22-06-2022
 
 FOLDERS: Damped, PML, Sponge

 COMMENTS: Each foldier contains a code in c++ that is composed by a main code and some libraries:
  
 MAIN CODE: main.c 
 
 INTERNAL LIBRARIES: "analisis.h", "pstd_opt.h"
 
 EXTERNAL LIBRARIES: <fftw3.h>, <omp.h>
 
 COMPILATION: g++ -lfftw3 -fopenmp -O3 *.c -o executable
 
 EXECUTION: ./executable
 
 OUTPUT: binary file "Energy.dat"
 
COMMENTS: This code calculates the energy of 600 different simulations defined in an homogeneous cube of 500x500x500m and a propagation velocity of c=2000m/s.
A ricker wavelet of f=10Hz is emitted at the center of the cube.
In all the cases, we also  fix the spatial sampling Dx=40m, the temporal step Dt=0.002s and the total simulation time 2s.
We use the numerical method corresponding to the foldier to solve this problem.
The 600 simulations are defined varying the tuple (Nabl,Absorbing parameter) in 30x20 different cases according to Table 1 in the paper.

DATA ANALYSIS CODE: energy.m

EXECUTION: octave or matlab

OUTPUT: Fig 2 and Table 2. 

COMMENTS: Each folder contains a code in matlab that generates the results of Fig. 2 and Table 2.



**Test 2**
--------- 

AUTHOR: Carlos Spa

 DATE: 22-06-2022
 
 
 FOLDERS: damped, sponge, euler  (D,S,E)
 
 COMMENTS: Each folder contains a code in c++ that is composed by a main code and some libraries:
 
 
 MAIN CODE: main.c 
 
 INTERNAL LIBRARIES: "analisis.h", "pstd_opt.h", "gestor_sim.h"
 
 EXTERNAL LIBRARIES: <fftw3.h>, <omp.h>
 
 COMPILATION: g++ -lfftw3 -fopenmp -O3 *.c -o executable
 
 EXECUTION: ./executable
 
 OUTPUT: binary file 2x7x3 Files: "P_type_Acc_abl_method.dat" where type=F,B  abl=1,2,3,4,5,6,7 and method=D,S,E
 
COMMENTS: This code calculates the forward/backward  (x2) simulations for each different method (D,S,E) and with the cases of table 2 defined in an homogeneous 
cube of 4000x4000x4000m and a propagation velocity of c=2000m/s.
In all the cases, we also  fix the spatial sampling Dx=40m, the temporal step Dt=0.002s and the total simulation time 4s.
For the Forward simulation, A ricker wavelet of f=10Hz is emitted at the position (Nx/2,Ny/2,4) of the cube.
For the Backward dsimulation, Three ricker wavelets of f=10Hz are emitted at the position (Nx/2,Ny/2,4) of the cube.
The output data of each method is recorded at the corresponding foldiers. 
The 21x2 simulations are defined varying the tuple (Nabl,Absorbing parameter)  according to Table 2 in the paper.

DATA ANALYSIS CODE: energy.m

EXECUTION: octave or matlab

OUTPUT: Table 3. 

 COMMENTS: Each folder contains a code in matlab that generates the results of Table 3
 

**Test 3**
-------------

 AUTHOR: Carlos Spa

DATE: 22-06-2022
 
FOLDERS: malla2, malla2/damped,malla2/sponge, malla2/euler 
 
COMMENTS: These folders contain different binary files used as input files in the code
 
 MAIN CODE: main.c 
 
 INTERNAL LIBRARIES: "analisis.h", "pstd_opt.h", "gestor_sim.h"
 
 EXTERNAL LIBRARIES: <fftw3.h>, <omp.h>
 
 COMPILATION:g++ -lfftw3 -fopenmp -O3 *.c -o executable
 
 EXECUTION: ./executable
 
 INPUT: A binary mesh is required and is located ./mallas2/METHOD/3d.outAcc_NUMBER.data where METHOD=damped,sponge,euler and NUMBER=1,2,3,4,5,6,7
 
 OUTPUT: 7x3 binary files EnergyAcc_NUMBER_METHOD where NUMBER=1,2,3,4,5,6,7 and METHOD:D,S,E
 
COMMENTS: This code calculates the  simulation of a unit impulse signal in a SEG SALT EAGE model for each different method (D,S,E) and with the cases of table 2 defined in anhomogeneous cube of 4000x4000x4000m and a propagation velocity of c=2000m/s.
In all the cases, we also  fix the spatial sampling Dx=40m, the temporal step Dt=0.002s and the total simulation time 4s.
For the Forward simulation, A ricker wavelet of f=10Hz is emitted at the position (Nx/2,Ny/2,4) of the cube.
The output data of each method is recorded at the corresponding . 
The 21x2 simulations are defined varying the tuple (Nabl,Absorbing parameter)  according to Table 2 in de paper.

DATA ANALYSIS CODE: energy.m

EXECUTION: octave or matlab

OUTPUT: Table 4. 


 
 
 
