
#ifndef PSTD_OPT_H
#define PSTD_OPT_H
#include <malloc.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <omp.h>
#include <fftw3.h>




typedef struct
{	
	
	int *EUx,*EUy,*EUz,*EVy,*EVx,*EVz;
	int PML;
	int nLTx,nLTy,nLTz;

}scenario;

typedef struct
{
	int nFFT,N,N2,SN;
	scenario s;
	float CteDer,CteDer2,sigma,y,nu,tau;
	int nx,ny,nz,nx2,ny2,nz2,nPML,nt;
	fftwf_complex *Der, *Der2, *Der3, *Der4,*Der5,*Der6,*DerP,*DerDP;
    float *DerR, *Der2R, *Der3R, *Der4R,*Der5R,*Der6R,*DerPR,*DerDPR;
	fftwf_plan planfordDer,planbackDer,planfordDer2,planbackDer2,planfordDer3,planbackDer3,planfordDer4,planbackDer4,planfordDer5,planbackDer5,planfordDer6,planbackDer6,planfordDer3d,planbackDer3d,planfordDer3d2,planbackDer3d2;	
	float *vx,*vy,*vz,*pxn,*pyn,*pzn,*pxa,*pya,*pza,*px,*py,*pz,a1,a2,cerjan,*PXn,*DtPXn,*PYn,*DtPYn,*PZn,*DtPZn,*PX,*PY,*p,*dpxn,*dpyn,*dpzn,*dpxa,*dpya,*dpza;
	float Dx,S,Lx,Ly,Lz,Dt,*source,*sourceB,c,*vel,*rho,*rhoVx,*rhoVy,*rhoVz;
} domain;


void initDerivatives(domain* D);

void source3D(domain* D,int type,int nS,int t);
void initParamEuler(domain* D,int Lx,int Ly,int Lz,float Dx,float Dt, float S,float tiempo, int nPML,float sigma,float cerjan);
void defineDomain3D(domain* D,int malla);
void initArrays3D(domain* D);
void initDerivatives(domain* D);


void executeDerShift(domain* D,fftwf_complex* p,fftwf_plan planfordDer,fftwf_plan planbackdDer,int N);
void executeDerShift2(domain* D,fftwf_complex* p,fftwf_plan planfordDer,fftwf_plan planbackDer,int N);


void updatePx(domain* D);
void updateVx(domain* D);
void updatePy(domain* D);
void updateVy(domain* D);
void updatePz(domain* D);
void updateVz(domain* D);
void tickEuler(domain* D);

void ricker_wavelet_tis0(domain* D,float dt,int nt,float t0,float f0,int l);

void sourceDamped(domain* D,int type,int nS,int t);
void initParamDamped(domain* D,int Lx,int Ly,int Lz,float Dx,float Dt, float S, float tiempo,int nPML,float sigma,float cerjan);
void initArraysDamped(domain* D);
void defineDomainDamped(domain* D,int malla);
void initDerivativesDamped(domain* D);
void executeDerDamped(domain* D,fftwf_complex* p,fftwf_plan planfordDer,fftwf_plan planbackDer,int N);
void updatePxDamped(domain* D);
void updatePyDamped(domain* D);
void updatePzDamped(domain* D);
void tickDamped(domain* D);
void freeEuler(domain* D);
void freeDamped(domain* D);
void freeSponge(domain* D);
void sourceSponge(domain* D,int type,int nS,int t);
void initParamSponge(domain* D,int Lx,int Ly,int Lz,float Dx,float Dt, float S, float tiempo,int nPML,float sigma,float cerjan);
void initArraysSponge(domain* D);
void defineDomainSponge(domain* D,int malla);
void initDerivativesSponge(domain* D);
void executeDerSponge(domain* D,fftwf_complex* p,fftwf_plan planfordDer,fftwf_plan planbackDer,int N);
void updatePxSponge(domain* D);
void updatePySponge(domain* D);
void updatePzSponge(domain* D);
void tickSponge(domain* D);


#endif
