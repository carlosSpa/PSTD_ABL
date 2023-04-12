#ifndef GESTOR_SIM_H
#define GESTOR_SIM_H

#include "pstd_opt.h"
#include "analisis.h"


void simulacionEuler(domain* D,results* R,int Snode,int type,int malla);
void simulacionDamped(domain* D,results* R,int Snode,int type,int malla);
void simulacionSponge(domain* D,results* R,int Snode,int type,int malla);
#endif
