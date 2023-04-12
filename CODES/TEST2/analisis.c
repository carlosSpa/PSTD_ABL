#include "analisis.h"
#include <omp.h> 



void Pressure1DEuler(domain *D,results *R,int t,int malla,int type){
	int nData,N,node;
	float pressure;
	
	pressure=0.0;
	if (t==0){
     
        if(type==0){
            if(malla==0)
                R->r1=fopen("./euler/PFAcc1E.dat","wb");
            if(malla==1)
                R->r1=fopen("./euler/PFAcc2E.dat","wb");
            if(malla==2)
                R->r1=fopen("./euler/PFAcc3E.dat","wb");
            if(malla==3)
                R->r1=fopen("./euler/PFAcc4E.dat","wb");
            if(malla==4)
                R->r1=fopen("./euler/PFAcc5E.dat","wb");
            if(malla==5)
                R->r1=fopen("./euler/PFAcc6E.dat","wb");
            if(malla==6)
                R->r1=fopen("./euler/PFAcc7E.dat","wb");
        }
        else{
            
            if(malla==0)
                R->r1=fopen("./euler/PBAcc1E.dat","wb");
            if(malla==1)
                R->r1=fopen("./euler/PBAcc2E.dat","wb");
            if(malla==2)
                R->r1=fopen("./euler/PBAcc3E.dat","wb");
            if(malla==3)
                R->r1=fopen("./euler/PBAcc4E.dat","wb");
            if(malla==4)
                R->r1=fopen("./euler/PBAcc5E.dat","wb");
            if(malla==5)
                R->r1=fopen("./euler/PBAcc6E.dat","wb");
            if(malla==6)
                R->r1=fopen("./euler/PBAcc7E.dat","wb");
            
            
            
        }
        
    }
	
    for(int k=0;k<D->nz;k++){
        node=D->nx/2-1+(D->ny/2-1)*D->nx+k*D->nx*D->ny;
        if((k>D->nPML-1)&&(k<D->nz-D->nPML)){
            pressure=(D->pxn[node]+D->pyn[node]+D->pzn[node]);
            nData=fwrite(&pressure,sizeof(float),1,R->r1);
        }			
			
        
	}
	
	if(t==D->nt-1)
		fclose(R->r1);

}

void Pressure1DDamped(domain *D,results *R,int t,int malla,int type){
	int nData,N,node;
	float pressure;
	
	pressure=0.0;
	if (t==0){
     
        if(type==0){
            if(malla==0)
                R->r1=fopen("./damped/PFAcc1D.dat","wb");
            if(malla==1)
                R->r1=fopen("./damped/PFAcc2D.dat","wb");
            if(malla==2)
                R->r1=fopen("./damped/PFAcc3D.dat","wb");
            if(malla==3)
                R->r1=fopen("./damped/PFAcc4D.dat","wb");
            if(malla==4)
                R->r1=fopen("./damped/PFAcc5D.dat","wb");
            if(malla==5)
                R->r1=fopen("./damped/PFAcc6D.dat","wb");
            if(malla==6)
                R->r1=fopen("./damped/PFAcc7D.dat","wb");
        }
        else{
            
            if(malla==0)
                R->r1=fopen("./damped/PBAcc1D.dat","wb");
            if(malla==1)
                R->r1=fopen("./damped/PBAcc2D.dat","wb");
            if(malla==2)
                R->r1=fopen("./damped/PBAcc3D.dat","wb");
            if(malla==3)
                R->r1=fopen("./damped/PBAcc4D.dat","wb");
            if(malla==4)
                R->r1=fopen("./damped/PBAcc5D.dat","wb");
            if(malla==5)
                R->r1=fopen("./damped/PBAcc6D.dat","wb");
            if(malla==6)
                R->r1=fopen("./damped/PBAcc7D.dat","wb");
            
            
            
        }
        
    }
	
    for(int k=0;k<D->nz;k++){
        node=D->nx/2-1+(D->ny/2-1)*D->nx+k*D->nx*D->ny;
        if((k>D->nPML-1)&&(k<D->nz-D->nPML)){
            pressure=(D->pxn[node]+D->pyn[node]+D->pzn[node]);
            nData=fwrite(&pressure,sizeof(float),1,R->r1);
        }			
			
        
	}
	
	if(t==D->nt-1)
		fclose(R->r1);

}

void Pressure1DSponge(domain *D,results *R,int t,int malla,int type){
	int nData,N,node;
	float pressure;
	
	pressure=0.0;
	if (t==0){
     
        if(type==0){
            if(malla==0)
                R->r1=fopen("./sponge/PFAcc1S.dat","wb");
            if(malla==1)
                R->r1=fopen("./sponge/PFAcc2S.dat","wb");
            if(malla==2)
                R->r1=fopen("./sponge/PFAcc3S.dat","wb");
            if(malla==3)
                R->r1=fopen("./sponge/PFAcc4S.dat","wb");
            if(malla==4)
                R->r1=fopen("./sponge/PFAcc5S.dat","wb");
            if(malla==5)
                R->r1=fopen("./sponge/PFAcc6S.dat","wb");
            if(malla==6)
                R->r1=fopen("./sponge/PFAcc7S.dat","wb");
        }
        else{
            
            if(malla==0)
                R->r1=fopen("./sponge/PBAcc1S.dat","wb");
            if(malla==1)
                R->r1=fopen("./sponge/PBAcc2S.dat","wb");
            if(malla==2)
                R->r1=fopen("./sponge/PBAcc3S.dat","wb");
            if(malla==3)
                R->r1=fopen("./sponge/PBAcc4S.dat","wb");
            if(malla==4)
                R->r1=fopen("./sponge/PBAcc5S.dat","wb");
            if(malla==5)
                R->r1=fopen("./sponge/PBAcc6S.dat","wb");
            if(malla==6)
                R->r1=fopen("./sponge/PBAcc7S.dat","wb");
            
            
            
        }
        
    }
	
    for(int k=0;k<D->nz;k++){
        node=D->nx/2-1+(D->ny/2-1)*D->nx+k*D->nx*D->ny;
        if((k>D->nPML-1)&&(k<D->nz-D->nPML)){
            pressure=(D->pxn[node]+D->pyn[node]+D->pzn[node]);
            nData=fwrite(&pressure,sizeof(float),1,R->r1);
        }			
			
        
	}
	
	if(t==D->nt-1)
		fclose(R->r1);
}


