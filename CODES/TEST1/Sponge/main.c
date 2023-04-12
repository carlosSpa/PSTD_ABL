/*
 ====================================================================================================================================================================
 AUTHOR: Carlos Spa
 DATE: 22-06-2020
 NAME: main.c 
 INTERNAL LIBRARIES: "analisis.h", "pstd_opt.h"
 EXTERNAL LIBRARIES: <fftw3.h>, <omp.h>
 COMPILATION:g++ -lfftw3 -fopenmp -O3 *.c -o executable
 EXECUTION: ./executable
 OUTPUT: binary file "Energy.dat"
 =====================================================================================================================================================================
COMMENTS: This code calculates the energy of 600 different simulations defined in an homogeneous cube of 500x500x500m and a propagation velocity of c=2000m/s.
A ricker wavelet of f=10Hz is emitted at the center of the cube.
In all the cases, we also  fix the spatial sampling Dx=40m, the temporal step Dt=0.002s and the total simulation time 2s.
We use the SPONGE METHOD formulation to solve this problem.
The 600 simulations are defined varying the tuple (Nabl,cerjan) in 30x20 different cases.
======================================================================================================================================================================
*/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "pstd_opt.h"
#include "analisis.h"
#include <fftw3.h>
#include <omp.h>
int main(){

	domain D;
	results R;
	int t,nx,ny,nz,nt,l,Snode,typeBC,Nb,Nparam;
	Nb=30;//Number of  cases of absorbing boundary layers Nabl
	Nparam=20;//Number of cases of absorbing parameters
	
	double S,Dx,Dt,Tiempo,c,Lx,Ly,Lz,t0,f0,sigma,cerjan;
	Lx=500;//Length in x dimension (meters)
	Ly=500;//Lenght in y dimension (meters)
	Lz=500;//Lenght in z dimension (meters)
	c=2000;//Propagation Velocity (meter/second)
	
	Tiempo=2;//Total time of simulation (seconds)
	Dx=40;//spatial discretization
	sigma=0.0;//Absorbing parameter set to 0
	cerjan=0.0;//Absorbing parameter set to 0
	Dt=0.002;//Temporal Discretization
	int auxPar=0;
	
	nt=(int)(Tiempo/Dt);//number of temporal nodes
	nx=(int)(Lx/Dx);//number of nodes at x direction
	ny=(int)(Ly/Dx);//number of nodes at y direction
	nz=(int)(Lz/Dx);//number of nodes at z direction
    t0=0.1;//Delay on the ricker expressed in seconds
    f0=10;//Central freq of the ricker wavelet
	int nPML;
    
	printf("Cube of 500x500x500 meters with numerical values: \n");
    printf("nx=%d, ny=%d, nz=%d with Dx=%lf and Dt=%lf. \n",nx,ny,nz,Dx,Dt);
    printf("The total time of  %lf s., or nt=%d iterations\n",Tiempo,nt);
    
    //Start Firt Kernel: Loop number of nodes
	for(int m=0;m<Nb;m++){

		//Number of absorbing bounday nodes from 4 to 34
		nPML=2*(4+m);//Number of absorbing bounday nodes from 4 to 34
		
		Snode=(int)((nx+nPML)/2-1)+(int)(((ny+nPML)/2-1)*(nx+nPML)+((nz+nPML)/2-1)*(ny+nPML)*(nx+nPML));//Nodal Location of the source.(In this case, in the middle)
		initParam(&D,Lx,Ly,Lz,Dx,Dt,S,Tiempo,nPML,sigma,cerjan,Nb,Nparam);
		initArrays(&D);
		defineDomain(&D);
		ricker_wavelet_tis0(&D,Dt,D.nt,t0,f0,0);
		initDerivatives(&D);
        //Start Second Loop: Absorbing parameters
		for(int k=0;k<Nparam;k++){
			t=0;
			D.cerjan=0.001+0.0404*k/19;
			tick(&D);tick(&D);//Data set to 0
			Energy(&D,&R,t,m,k);
            printf("Simulation m=%d de Nb=%d y k=%d de Nparam= %d\n",m+1,Nb,k+1,Nparam);
            //Start the simulation Kernel
			for(t=1;t<D.nt;t++){
#pragma omp parallel for 
				for(auxPar=0;auxPar<3;auxPar++){
					if (auxPar==0)
						updatePx(&D);
					if(auxPar==1)
						updatePy(&D);
					if(auxPar==2)
						updatePz(&D);
				}
#pragma end 
				source(&D,0,Snode,t);
				tick(&D);
				
				Energy(&D,&R,t,m,k);
		
	
			}
			//End simulation Kernel
		}
		freeVectors(&D);
		//End Second Loop
		
	}
    //End First Loop

}		
		

