#include "analisis.h"
#include <omp.h> 



void Energy(domain *D,results *R,int t,int malla){
	int nData,N,node;
	float energy;
	
	energy=0.0;
	if (t==0){
     
        if(malla==0)
            R->r1=fopen("./EnergyRef.dat","wb");
        if(malla==1)
            R->r1=fopen("./EnergyAcc1E.dat","wb");
        if(malla==2)
            R->r1=fopen("./EnergyAcc2E.dat","wb");
        if(malla==3)
            R->r1=fopen("./EnergyAcc3E.dat","wb");
        if(malla==4)
            R->r1=fopen("./EnergyAcc4E.dat","wb");
        if(malla==5)
            R->r1=fopen("./EnergyAcc5E.dat","wb");
        if(malla==6)
            R->r1=fopen("./EnergyAcc6E.dat","wb");
        if(malla==7)
            R->r1=fopen("./EnergyAcc7E.dat","wb");   
        
    }
	for(int i=0;i<D->nx;i++){
		for(int j=0;j<D->ny;j++){
			for(int k=0;k<D->nz;k++){
				node=i+j*D->nx+k*D->nx*D->ny;
				if((i>D->nPML-1)&&(i<D->nx-D->nPML)&&(j>D->nPML-1)&&(j<D->ny-D->nPML)&&(k>D->nPML-1)&&(k<D->nz-D->nPML)){
					energy=energy+(D->px[node]+D->py[node]+D->pz[node])*(D->px[node]+D->py[node]+D->pz[node]);
					
				}			
			}
		}
	}
	nData=fwrite(&energy,sizeof(float),1,R->r1);
	if(t==D->nt-1)
		fclose(R->r1);

}

void EnergyDamped(domain *D,results *R,int t,int malla){
	int nData,N,node;
	float energy;
	
	energy=0.0;
	if (t==0){
     
        if(malla==1)
            R->r1=fopen("./EnergyAcc1D.dat","wb");
        if(malla==2)
            R->r1=fopen("./EnergyAcc2D.dat","wb");
        if(malla==3)
            R->r1=fopen("./EnergyAcc3D.dat","wb");
        if(malla==4)
            R->r1=fopen("./EnergyAcc4D.dat","wb");
        if(malla==5)
            R->r1=fopen("./EnergyAcc5D.dat","wb");
        if(malla==6)
            R->r1=fopen("./EnergyAcc6D.dat","wb");
        if(malla==7)
            R->r1=fopen("./EnergyAcc7D.dat","wb");   
        
    }
	for(int i=0;i<D->nx;i++){
		for(int j=0;j<D->ny;j++){
			for(int k=0;k<D->nz;k++){
				node=i+j*D->nx+k*D->nx*D->ny;
				if((i>D->nPML-1)&&(i<D->nx-D->nPML)&&(j>D->nPML-1)&&(j<D->ny-D->nPML)&&(k>D->nPML-1)&&(k<D->nz-D->nPML)){
					energy=energy+(D->px[node]+D->py[node]+D->pz[node])*(D->px[node]+D->py[node]+D->pz[node]);
					
				}			
			}
		}
	}
	nData=fwrite(&energy,sizeof(float),1,R->r1);
	if(t==D->nt-1)
		fclose(R->r1);

}

void EnergySponge(domain *D,results *R,int t,int malla){
	int nData,N,node;
	float energy;
	
	energy=0.0;
	if (t==0){
     
        if(malla==1)
            R->r1=fopen("./EnergyAcc1S.dat","wb");
        if(malla==2)
            R->r1=fopen("./EnergyAcc2S.dat","wb");
        if(malla==3)
            R->r1=fopen("./EnergyAcc3S.dat","wb");
        if(malla==4)
            R->r1=fopen("./EnergyAcc4S.dat","wb");
        if(malla==5)
            R->r1=fopen("./EnergyAcc5S.dat","wb");
        if(malla==6)
            R->r1=fopen("./EnergyAcc6S.dat","wb");
        if(malla==7)
            R->r1=fopen("./EnergyAcc7S.dat","wb");   
        
    }
	for(int i=0;i<D->nx;i++){
		for(int j=0;j<D->ny;j++){
			for(int k=0;k<D->nz;k++){
				node=i+j*D->nx+k*D->nx*D->ny;
				if((i>D->nPML-1)&&(i<D->nx-D->nPML)&&(j>D->nPML-1)&&(j<D->ny-D->nPML)&&(k>D->nPML-1)&&(k<D->nz-D->nPML)){
					energy=energy+(D->pxn[node]+D->pyn[node]+D->pzn[node])*(D->pxn[node]+D->pyn[node]+D->pzn[node]);
					
				}			
			}
		}
	}
	nData=fwrite(&energy,sizeof(float),1,R->r1);
	if(t==D->nt-1)
		fclose(R->r1);

}


