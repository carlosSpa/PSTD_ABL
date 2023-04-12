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
	
	int *EUx,*EUy,*EUz;
	int PML;
	int nLTx,nLTy,nLTz;

}scenario;

typedef struct
{
	int nFFT;
	scenario s;
	double CteDer,sigma;
	int nx,ny,nz,nx2,ny2,nz2,nPML,nt,Nb,Nparam;
	fftw_complex *Der, *Der2, *Der3;
	fftw_plan planfordDer,planbackDer,planfordDer2,planbackDer2,planfordDer3,planbackDer3;	
	double *pxa,*pya,*pza,*pxn,*pyn,*pzn,*px,*py,*pz,a1,a2,cerjan;
	double Dx,S,Lx,Ly,Lz,Dt,*source,c,*vel;
} domain;

void source(domain* D,int type,int nS,int t);
void initParam(domain* D,double Lx,double Ly,double Lz,double Dx,double Dt, double S,double tiempo, int nPML,double sigma,double cerjan,int Nb,int Nparam);
void defineDomain(domain* D);
void initArrays(domain* D);

void initDerivatives(domain* D);
void executeDer(domain* D,fftw_complex* p,fftw_plan planfordDer,fftw_plan planbackdDer,int N);
void updatePx(domain* D);
void updatePy(domain* D);
void updatePz(domain* D);

void ricker_wavelet_tis0(domain* D,double dt,int nt,double t0,double f0,int l);
void tick(domain* D);
void freeVectors(domain* D);
#endif
