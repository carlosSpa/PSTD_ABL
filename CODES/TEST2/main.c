 /*
 ====================================================================================================================================================================
 AUTHOR: Carlos Spa
 DATE: 22-06-2020
 NAME: main.c 
 INTERNAL LIBRARIES: "analisis.h", "pstd_opt.h", "gestor_sim.h"
 EXTERNAL LIBRARIES: <fftw3.h>, <omp.h>
 COMPILATION:g++ -lfftw3 -fopenmp -O3 *.c -o executable
 EXECUTION: ./executable
 OUTPUT: 7x3x2 binary files: "PBAcc_NUMBER_METHOD.dat" and "PFAcc_NUMBER_METHOD" where NUMBER=1,2,3,4,5,6,7 and METHOD:D (Damped), S (Sponge) and E (Euler PML).
 =====================================================================================================================================================================
COMMENTS: This code calculates the forward/backward simulations defined in an homogeneous cube of 4000x4000x4000m and a propagation velocity of c=2000m/s. 
This information is essential for obtaining a gepohysical imaging that is built with only one source/receiver located at z=0 and x=y=50 nodes of distance. 
For the forward simulation,  a ricker wavelet of f=10Hz is emitted whereas for the backward simulation three ricker wavelet of 10Hz are emitted at 1, 2 and 3 seconds.
In all the cases, we also  fix the spatial sampling Dx=40m, the temporal step Dt=0.002s and the total simulation time 4s.
We use the DAMPED WAVE EQUATION, the SPONGE BOUNDARY LAYERS and the EULER PML  formulations to solve this problem.
The  simulations are defined varying the duple (Nabl,sigma) in 7 different cases according to Table 2.
======================================================================================================================================================================
*/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "pstd_opt.h"
#include "analisis.h"
#include "gestor_sim.h"
#include <fftw3.h>
#include <omp.h>
int main(){

	domain D;
	results R;
	
  
	int t,nx,ny,nz,nt,l,Snode,type;
	l=0;
	float S,Dx,Dt,Tiempo,c,rho,t0,f0,sigma,sponge;
    float sigmaE[7],sigmaD[7],sigmaS[7];
    int NablE[7],NablD[7],NablS[7];
    
//The  7x3 cases presented in Table 2 are considered taking the different pairs of (Nabl,Absorption) for the simulations.
    
    NablE[0]=4*2;NablE[1]=4*2;NablE[2]=5*2;NablE[3]=6*2;NablE[4]=9*2;NablE[5]=12*2;NablE[6]=16*2;
    sigmaE[0]=0.065;sigmaE[1]=0.097;sigmaE[2]=0.097;sigmaE[3]=0.097;sigmaE[4]=0.097;sigmaE[5]=0.065;sigmaE[6]=0.065;
    
    NablD[0]=2*5;NablD[1]=2*7;NablD[2]=2*10;NablD[3]=2*14;NablD[4]=2*18;NablD[5]=2*25;NablD[6]=2*32;
    sigmaD[0]=0.086;sigmaD[1]=0.071;sigmaD[2]=0.056;sigmaD[3]=0.041;sigmaD[4]=0.041;sigmaD[5]=0.025;sigmaD[6]=0.025;
    
    NablS[0]=2*7;NablS[1]=2*8;NablS[2]=2*11;NablS[3]=2*14;NablS[4]=2*17;NablS[5]=2*23;NablS[6]=2*30;
    sigmaS[0]=0.031;sigmaS[1]=0.03;sigmaS[2]=0.02;sigmaS[3]=0.016;sigmaS[4]=0.012;sigmaS[5]=0.007;sigmaS[6]=0.005;
//=======================================================================================================================0
    
    Tiempo=4;//Total simulation time
	Dx=40;//Spatial Discretization
	Dt=0.002;//Temporal Discretization
	nt=(int)(Tiempo/Dt);//Temporal nodes
    nx=100;//x-spatial nodes
	ny=100;//y-spatial nodes
	nz=100;//z-spatial nodes

//ricker wavelet of 10 Hz and t0=0.1 seconds
    t0=0.1;
    f0=10;
//=================================
    int malla;
    malla=0;
    int nPML;
    
    
// MAIN LOOP: We simulate the backward/forward simulations of the 7 different cases for each different method.
    
    for(malla=0;malla<7;malla++){
        
        //EULER PML SIMULATIONS
        nPML=NablE[malla];
        sigma=sigmaE[malla];
        Snode=(int)((nx+nPML)/2-1)+(int)(((ny+nPML)/2-1)*(nx+nPML)+((nPML)/2+4)*(ny+nPML)*(nx+nPML));
        D.SN=Snode;
        type=0;//FORWARD SIMULATION 
        initParamEuler(&D,nx,ny,nz,Dx,Dt,S,Tiempo,nPML,sigma,0);
        ricker_wavelet_tis0(&D,Dt,D.nt,t0,f0,type);
        simulacionEuler(&D,&R,Snode,type,malla);//This function is at gestor_sim.h library
        
        type=1;//BACKWARD SIMULATION
        initParamEuler(&D,nx,ny,nz,Dx,Dt,S,Tiempo,nPML,sigma,0);
        ricker_wavelet_tis0(&D,Dt,D.nt,t0,f0,type);
        simulacionEuler(&D,&R,Snode,type,malla);//This function is at gestor_sim.h library
        
        //DAMPED WAVE EQUATION SIMULATIONS
        nPML=NablD[malla];
        sigma=sigmaD[malla];
        Snode=(int)((nx+nPML)/2-1)+(int)(((ny+nPML)/2-1)*(nx+nPML)+((nPML)/2+4)*(ny+nPML)*(nx+nPML));
        D.SN=Snode;
        type=0;//FORWARD SIMULATION
        initParamDamped(&D,nx,ny,nz,Dx,Dt,S,Tiempo,nPML,sigma,0);
        ricker_wavelet_tis0(&D,Dt,D.nt,t0,f0,type);
        simulacionDamped(&D,&R,Snode,type,malla);//This function is at gestor_sim.h library
      
        type=1;//BACKWARD SIMULATION
        initParamDamped(&D,nx,ny,nz,Dx,Dt,S,Tiempo,nPML,sigma,0);
        ricker_wavelet_tis0(&D,Dt,D.nt,t0,f0,type);
        simulacionDamped(&D,&R,Snode,type,malla);//This function is at gestor_sim.h library
        
        //SPONGE ABSORBING LAYERS SIMULATIONS
        nPML=NablS[malla];
        sigma=0.0;
        sponge=sigmaS[malla];
        Snode=(int)((nx+nPML)/2-1)+(int)(((ny+nPML)/2-1)*(nx+nPML)+((nPML)/2+4)*(ny+nPML)*(nx+nPML));
        D.SN=Snode;
        
        type=0;//FORWARD SIMULATION
        initParamSponge(&D,nx,ny,nz,Dx,Dt,S,Tiempo,nPML,sigma,sponge);
        ricker_wavelet_tis0(&D,Dt,D.nt,t0,f0,type);
        simulacionSponge(&D,&R,Snode,type,malla);//This function is at gestor_sim.h library
        
        type=1;//BACKWARD SIMULATION
        initParamSponge(&D,nx,ny,nz,Dx,Dt,S,Tiempo,nPML,sigma,sponge);
        ricker_wavelet_tis0(&D,Dt,D.nt,t0,f0,type);
        simulacionSponge(&D,&R,Snode,type,malla);//This function is at gestor_sim.h library
            
        
    }
    
	
}
		
		

