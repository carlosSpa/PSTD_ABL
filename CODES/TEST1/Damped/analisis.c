#include "analisis.h"
#include <omp.h> 




void Energy(domain *D,results *R,int t,int m,int k){
	int nData,N,node;
	double energy;
	
	energy=0.0;
	if ((t==0)&&(m==0)&&(k==0))
		R->r1=fopen("./Energy.dat","wb");
	for(int i=0;i<D->nx;i++){
		for(int j=0;j<D->ny;j++){
			for(int k=0;k<D->nz;k++){
				node=i+j*D->nx+k*D->nx*D->ny;
				if((i>D->nPML-1)&&(i<D->nx-D->nPML)&&(j>D->nPML-1)&&(j<D->ny-D->nPML)&&(k>D->nPML-1)&&(k<D->nz-D->nPML)){
					energy=energy+(D->pxa[node]+D->pya[node]+D->pza[node])*(D->pxa[node]+D->pya[node]+D->pza[node]);
					
				}			
			}
		}
	}
	nData=fwrite(&energy,sizeof(double),1,R->r1);
	if((t==D->nt-1)&&(m==D->Nb-1)&&(k==D->Nparam-1))
		fclose(R->r1);

}

