#ifndef ANALISIS_H
#define ANALISIS_H
#include <malloc.h>
#include "pstd_opt.h"

typedef struct
{
	int mNodes;
	float Recep1,Recep2,Recep3;
	FILE *r1,*r2,*r3,*r4;
	
} results;

//Declare Functions


void Energy(domain *D,results *R,int t,int malla);
void EnergyDamped(domain *D,results *R,int t,int malla);
void EnergySponge(domain *D,results *R,int t,int malla);
#endif
