#include "gestor_sim.h"
 
void simulacionEuler(domain* D,results* R,int Snode,int type,int malla){
    int t;
    t=0;
    
    initArrays3D(D);
	
	
	defineDomain3D(D,malla);
	
	
	
	initDerivatives(D);
	
	
	printf("Simulacion Euler de un cubo de nx=%d ny=%d nz=%d y un Dx=%lf Dt=%lf\n",D->nx,D->ny,D->nz,D->Dx,D->Dt);
    Pressure1DEuler(D,R,t,malla,type);
	for(t=1;t<D->nt;t++){
		
		if(t%500==0)
                printf("iteracion %d de %d\n",t,D->nt);
		
#pragma omp parallel for simd        
        for (int l=0;l<3;l++){
            
            if(l==0){
                updateVx(D);
                updatePx(D);
            }
            if(l==1){
                updateVy(D);
                updatePy(D);
            } 
            if(l==2){
                updateVz(D);
                updatePz(D);
            }
        
        }
#pragma for end
        source3D(D,type,Snode,t);
         Pressure1DEuler(D,R,t,malla,type);
        tickEuler(D);
       
    }
	freeEuler(D);

    
}

void simulacionDamped(domain* D,results* R,int Snode,int type,int malla){
    int t;
   
    initArraysDamped(D);
    
    defineDomainDamped(D,malla);
    
    initDerivativesDamped(D);
    t=0;
    printf("Simulacion Damped de un cubo de nx=%d ny=%d nz=%d y un Dx=%lf Dt=%lf\n",D->nx,D->ny,D->nz,D->Dx,D->Dt);
    Pressure1DDamped(D,R,t,malla,type);
	for(t=1;t<D->nt;t++){
    	
		if(t%500==0)
			printf("iteracion %d de %d\n",t,D->nt);
		
#pragma omp parallel for simd        
        for (int l=0;l<3;l++){
            if(l==0)
                updatePxDamped(D); 
            if(l==1)
                updatePyDamped(D); 
            if(l==2)
                updatePzDamped(D);
        }
#pragma for end
   
        sourceDamped(D,type,Snode,t);
        Pressure1DDamped(D,R,t,malla,type);
        tickDamped(D);
    
    }
    freeDamped(D);
}


void simulacionSponge(domain* D,results* R,int Snode,int type,int malla){
    int t;
   
    initArraysSponge(D);
   
    defineDomainSponge(D,malla);

    initDerivativesSponge(D);
     
    t=0;
    printf("Simulacion Sponge de un cubo de nx=%d ny=%d nz=%d y un Dx=%f Dt=%f\n",D->nx,D->ny,D->nz,D->Dx,D->Dt);
    Pressure1DSponge(D,R,t,malla,type);
	for(t=1;t<D->nt;t++){
    	
		if(t%500==0)
			printf("iteracion %d de %d\n",t,D->nt);
		
#pragma omp parallel for simd        
        for (int l=0;l<3;l++){
            if(l==0)
                updatePxSponge(D); 
            if(l==1)
                updatePySponge(D); 
            if(l==2)
                updatePzSponge(D);
        }
#pragma for end
   
        sourceSponge(D,type,Snode,t);
        Pressure1DSponge(D,R,t,malla,type);
        tickSponge(D);
    
    }
    
    freeSponge(D);
    
}
