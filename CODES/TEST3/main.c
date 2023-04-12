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
	
  
	int t,nx,ny,nz,nt,l,Snode;
	l=0;
	float S,Dx,Dt,Tiempo,c,rho,Lx,Ly,Lz,t0,f0,sigma,sponge;
	//Defino parametros fisicos 
    
	Tiempo=4;
	Dx=30;
    float sigmaE[8],sigmaD[7],sigmaS[7];
    int NablE[8],NablD[7],NablS[7];
    NablE[0]=120*2;NablE[1]=4*2;NablE[2]=4*2;NablE[3]=5*2;NablE[4]=6*2;NablE[5]=9*2;NablE[6]=12*2;NablE[7]=16*2;
    sigmaE[0]=0.065;sigmaE[1]=0.065;sigmaE[2]=0.097;sigmaE[3]=0.097;sigmaE[4]=0.097;sigmaE[5]=0.097;sigmaE[6]=0.065;sigmaE[7]=0.065;
    
    NablD[0]=2*5;NablD[1]=2*7;NablD[2]=2*10;NablD[3]=2*14;NablD[4]=2*18;NablD[5]=2*25;NablD[6]=2*32;
    sigmaD[0]=0.086;sigmaD[1]=0.071;sigmaD[2]=0.056;sigmaD[3]=0.041;sigmaD[4]=0.041;sigmaD[5]=0.025;sigmaD[6]=0.025;
    
    NablS[0]=2*7;NablS[1]=2*8;NablS[2]=2*11;NablS[3]=2*14;NablS[4]=2*17;NablS[5]=2*23;NablS[6]=2*30;
    sigmaS[0]=0.031;sigmaS[1]=0.03;sigmaS[2]=0.02;sigmaS[3]=0.016;sigmaS[4]=0.012;sigmaS[5]=0.007;sigmaS[6]=0.005;
    
    
	Dt=0.002;
	nt=(int)(Tiempo/Dt);
    nx=251;
	ny=251;
	nz=122;
    t0=0.1;f0=10;
    int malla;
    malla=0;
    int nPML;
   
     
    nPML=NablE[malla];
	sigma=sigmaE[malla];
    Snode=(int)((nx+nPML)/2-1)+(int)(((ny+nPML)/2-1)*(nx+nPML)+((nPML)/2+4)*(ny+nPML)*(nx+nPML));
    D.SN=Snode;
   
    initParamEuler(&D,Lx,Ly,Lz,Dx,Dt,S,Tiempo,nPML,sigma,0);
	ricker_wavelet_tis0(&D,Dt,D.nt,t0,f0,0);
   
	simulacionEuler(&D,&R,Snode,0,malla);
    
    
     
    
    for(malla=1;malla<8;malla++){
        
        
        nPML=NablE[malla];
        sigma=sigmaE[malla];
        Snode=(int)((nx+nPML)/2-1)+(int)(((ny+nPML)/2-1)*(nx+nPML)+((nPML)/2+4)*(ny+nPML)*(nx+nPML));
        D.SN=Snode;

        initParamEuler(&D,Lx,Ly,Lz,Dx,Dt,S,Tiempo,nPML,sigma,0);
        ricker_wavelet_tis0(&D,Dt,D.nt,t0,f0,0);
        simulacionEuler(&D,&R,Snode,0,malla);
     
        
        nPML=NablD[malla-1];
        sigma=sigmaD[malla-1];
        Snode=(int)((nx+nPML)/2-1)+(int)(((ny+nPML)/2-1)*(nx+nPML)+((nPML)/2+4)*(ny+nPML)*(nx+nPML));
        D.SN=Snode;
    
        initParamDamped(&D,Lx,Ly,Lz,Dx,Dt,S,Tiempo,nPML,sigma,0);
        ricker_wavelet_tis0(&D,Dt,D.nt,t0,f0,0);
        simulacionDamped(&D,&R,Snode,0,malla);
      
       
        
        
        
        nPML=NablS[malla-1];
        sigma=0.0;
        sponge=sigmaS[malla-1];
        Snode=(int)((nx+nPML)/2-1)+(int)(((ny+nPML)/2-1)*(nx+nPML)+((nPML)/2+4)*(ny+nPML)*(nx+nPML));
        D.SN=Snode;
        
        initParamSponge(&D,Lx,Ly,Lz,Dx,Dt,S,Tiempo,nPML,sigma,sponge);
        ricker_wavelet_tis0(&D,Dt,D.nt,t0,f0,0);
        simulacionSponge(&D,&R,Snode,0,malla);
            
        
    }
    
	
}
		
		

