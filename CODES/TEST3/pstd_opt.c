#include "pstd_opt.h"



void ricker_wavelet_tis0(domain* D,float dt,int nt,float t0,float f0,int l){

	D->source=(float*)malloc((D->nt)*sizeof(float));
	float B=0;float Pi;  Pi=3.141592654;
	B=Pi*Pi*f0*f0;
	if(l==0){
		for(int t=0;t<D->nt;t++){
			if(t!=0)
				D->source[t]=100*((1-2*B*(t*D->Dt-t0)*(t*D->Dt-t0))*exp(-B*(t*D->Dt-t0)*(t*D->Dt-t0))); 
				
			else
				D->source[t]=0.0;
		}
	}
	

}

void source3D(domain* D,int type,int nS,int t){
    double aux;aux=0.0;
	if(type==0){
        for(int n=0;n<t;n++)
            aux=aux+D->source[n];
		for(int i=-1;i<2;i++){
			for(int j=-1;j<2;j++){
				for(int k=-1;k<2;k++){
					D->pxn[nS+i+j*D->nx+k*D->nx*D->ny]=aux*exp(-1*((i)*(i)+(j)*(j)+(k)*(k)))+D->pxn[nS+i+j*D->nx+k*D->nx*D->ny];
					D->pyn[nS+i+j*D->nx+k*D->nx*D->ny]=aux*exp(-1*((i)*(i)+(j)*(j)+(k)*(k)))+D->pyn[nS+i+j*D->nx+k*D->nx*D->ny];
					D->pzn[nS+i+j*D->nx+k*D->nx*D->ny]=aux*exp(-1*((i)*(i)+(j)*(j)+(k)*(k)))+D->pzn[nS+i+j*D->nx+k*D->nx*D->ny];

				}
			}
		}
		
	}
	
}

void initParamEuler(domain* D,float Lx,float Ly,float Lz,float Dx,float Dt, float S, float tiempo,int nPML,float sigma,float cerjan){
	D->c=2000;float rho=1;
	D->a1=rho*Dt;
	D->a2=Dt/rho;
	D->S=D->c*Dt/Dx;
	D->Dx=Dx;
	D->nx=251+nPML;
	D->ny=251+nPML;
	D->nz=122+nPML;
	D->sigma=sigma;
	D->nPML=nPML/2;
	D->nt=(int)(tiempo/Dt);
	D->Dt=Dt;	
	D->cerjan=cerjan;
}


void initArrays3D(domain* D){
	int nodes;
	nodes=D->nx*D->ny*D->nz;	
	D->px=(float*)malloc((nodes)*sizeof(float));
	D->py=(float*)malloc((nodes)*sizeof(float));
	D->pz=(float*)malloc((nodes)*sizeof(float));
	D->vx=(float*)malloc((nodes)*sizeof(float));
	D->vy=(float*)malloc((nodes)*sizeof(float));
	D->vz=(float*)malloc((nodes)*sizeof(float));
	D->pxn=(float*)malloc((nodes)*sizeof(float));
	D->pyn=(float*)malloc((nodes)*sizeof(float));
	D->pzn=(float*)malloc((nodes)*sizeof(float));
	D->vel=(float*)malloc((nodes)*sizeof(float));
	D->s.EUx=(int*)malloc((nodes)*sizeof(int));
	D->s.EUy=(int*)malloc((nodes)*sizeof(int));
	D->s.EUz=(int*)malloc((nodes)*sizeof(int));

	for(int i=0;i<nodes;i++){
		D->px[i]=D->py[i]=D->pz[i]=D->vx[i]=D->vy[i]=D->vz[i]=D->pxn[i]=D->pyn[i]=D->pzn[i]=0.0;
		D->s.EUx[i]=D->s.EUy[i]=D->s.EUz[i]=0;D->vel[i]=0.0;
	}
}



void defineDomain3D(domain* D,int malla){
 
	for(int i=0;i<D->nx;i++){
		for(int j=0;j<D->ny;j++){
			for(int k=0;k<D->nz;k++){
				int node;
				node=i+j*D->nx+k*D->nx*D->ny;
				
				if(i<D->nPML)
					D->s.EUx[node]=D->nPML-i;
				if(i>D->nx-1-D->nPML)
					D->s.EUx[node]=i-(D->nx-1-D->nPML);
				if(j<D->nPML)
					D->s.EUy[node]=D->nPML-j;
				if(j>D->ny-1-D->nPML)
					D->s.EUy[node]=j-(D->ny-1-D->nPML);
				if(k<D->nPML){
					D->s.EUz[node]=D->nPML-k;
                }
				if(k>D->nz-1-D->nPML)
					D->s.EUz[node]=k-(D->nz-1-D->nPML);
              

				if(i==0){//0
					D->s.EUx[node]=40;D->s.EUy[node]=40;D->s.EUz[node]=40;}
				if(i==D->nx-1){
					D->s.EUx[node]=40;D->s.EUy[node]=40;D->s.EUz[node]=40;}	
				if(j==0){
					D->s.EUx[node]=40;D->s.EUy[node]=40;D->s.EUz[node]=40;}
				if(j==D->ny-1){
					D->s.EUx[node]=40;D->s.EUy[node]=40;D->s.EUz[node]=40;}
				if(k==0){
					D->s.EUx[node]=40;D->s.EUy[node]=40;D->s.EUz[node]=40;}
				if(k==D->nz-1){
					D->s.EUx[node]=40;D->s.EUy[node]=40;D->s.EUz[node]=40;}
				
			}
		}
	}


	int nData;FILE *r1;double av;
    if(malla==0)
        r1=fopen("./Mesh/Euler/3d.outRef.bin","rb");
    if(malla==1)
        r1=fopen("./Mesh/Euler/3d.outAcc2_1.bin","rb");
    if(malla==2)
        r1=fopen("./Mesh/Euler/3d.outAcc2_1.bin","rb");
    if(malla==3)
        r1=fopen("./Mesh/Euler/3d.outAcc3.bin","rb");
    if(malla==4)
        r1=fopen("./Mesh/Euler/3d.outAcc4.bin","rb");
    if(malla==5)
        r1=fopen("./Mesh/Euler/3d.outAcc5.bin","rb");
    if(malla==6)
        r1=fopen("./Mesh/Euler/3d.outAcc6.bin","rb");
    if(malla==7)
        r1=fopen("./Mesh/Euler/3d.outAcc7.bin","rb");
    
    
	for(int i=0;i<D->nx;i++){
		for(int j=0;j<D->ny;j++){
			for(int k=0;k<D->nz;k++){
				int node;
				node=i+j*D->nx+k*D->nx*D->ny;
				
				nData=fread(&av,sizeof(double),1,r1);
                D->vel[node]=(float)(av);
                

			}
		}
	}
	fclose(r1);
	
	
}


void initDerivatives(domain* D){

	double pi=3.14159265358979323846;
	int N,p2,N2;
	N=D->nx;
	p2=2;
	D->nFFT=0;
   
	if(N<D->ny)
		N=D->ny;
	if(N<D->nz)
		N=D->nz;
	
	N=N+D->nFFT;
   
    if(N%2!=0)
        N=N+1;
    D->N=N;
    N2=N/2+1;
    if(N2%2!=0){
        N=N+2;
        D->N=N;
        N2=N/2+1;
    }
    
    D->N2=N2;
    
	D->CteDer=2*pi/(N*D->Dx*N);
	D->Der=(fftwf_complex *)fftw_malloc(N2*sizeof(fftwf_complex));
    D->DerR=(float *)malloc(N*sizeof(float));
	D->planfordDer = fftwf_plan_dft_r2c_1d(N, D->DerR, D->Der,  FFTW_ESTIMATE);
	D->planbackDer = fftwf_plan_dft_c2r_1d(N, D->Der, D->DerR,  FFTW_ESTIMATE);
	
    D->Der2=(fftwf_complex *)fftw_malloc(N2*sizeof(fftwf_complex));
    D->Der2R=(float *)malloc(N*sizeof(float));
	D->planfordDer2 = fftwf_plan_dft_r2c_1d(N, D->Der2R, D->Der2,  FFTW_ESTIMATE);
	D->planbackDer2 = fftwf_plan_dft_c2r_1d(N, D->Der2, D->Der2R,  FFTW_ESTIMATE);
	
    
    D->Der3=(fftwf_complex *)fftw_malloc(N2*sizeof(fftwf_complex));
    D->Der3R=(float *)malloc(N*sizeof(float));
	D->planfordDer3 = fftwf_plan_dft_r2c_1d(N, D->Der3R, D->Der3,  FFTW_ESTIMATE);
	D->planbackDer3 = fftwf_plan_dft_c2r_1d(N, D->Der3, D->Der3R,  FFTW_ESTIMATE);
    

	

	for(int i=0;i<N;i++){
        D->DerR[i]=D->Der2R[i]=D->Der3R[i]=0.0; 
        if(i<N2){
            D->Der[i][0]=D->Der[i][1]=0.0;D->Der2[i][0]=D->Der2[i][1]=0.0;D->Der3[i][0]=D->Der3[i][1]=0.0;
           
        }
	}
}






void executeDerShift(domain* D,fftwf_complex* p,fftwf_plan planfordDer,fftwf_plan planbackDer,int N){
	float pi=3.14159265358979323846;
	fftwf_execute(planfordDer);
	for(int i=0;i<N/2+1;i++){
		float aux1,aux2;
			aux1=(-D->CteDer*i*p[i][1]*cos(D->CteDer*i*N*D->Dx/2)-D->CteDer*i*p[i][0]*sin(D->CteDer*i*N*D->Dx/2));
			aux2=(-D->CteDer*i*p[i][1]*sin(D->CteDer*i*N*D->Dx/2)+D->CteDer*i*p[i][0]*cos(D->CteDer*i*N*D->Dx/2));
		
		p[i][0]=aux1;
		p[i][1]=aux2;
	}

	fftwf_execute(planbackDer);
}

void executeDerShift2(domain* D,fftwf_complex* p,fftwf_plan planfordDer,fftwf_plan planbackDer,int N){
	float pi=3.14159265358979323846;
	fftwf_execute(planfordDer);
	for(int i=0;i<N/2+1;i++){
		float aux1,aux2;
		
			aux1=(-D->CteDer*i*p[i][1]*cos(D->CteDer*i*N*D->Dx/2)+D->CteDer*i*p[i][0]*sin(D->CteDer*i*N*D->Dx/2));
			aux2=(D->CteDer*i*p[i][1]*sin(D->CteDer*i*N*D->Dx/2)+D->CteDer*i*p[i][0]*cos(D->CteDer*i*N*D->Dx/2));
            
		p[i][0]=aux1;
		p[i][1]=aux2;
	}
   
	fftwf_execute(planbackDer);
   
}







void updatePx(domain* D){
	
	
	int N;float sigma;
	N=D->N;
	for(int j=0;j<D->ny;j++){
		for(int k=0;k<D->nz;k++){
			int node;
			node=j*D->nx+k*D->nx*D->ny;
			
			for(int i=0;i<D->nx;i++){
				D->DerR[i]=D->vx[node+i];
				
			}
			
			for(int i=0;i<N-D->nx;i++)
				D->DerR[i+D->nx]=0.0;

			
			executeDerShift(D,D->Der,D->planfordDer,D->planbackDer,N);
			
			
			for(int i=0;i<D->nx;i++){
				if(D->s.EUx[i+node]!=40){
					sigma=D->s.EUx[i+node]*D->sigma/(D->nPML);
					D->pxn[i+node]=(1-sigma)*D->px[i+node]-D->vel[i+node]*D->vel[i+node]*D->a1*D->DerR[i];
				
				}
				else
					D->pxn[i+node]=0.0;
			}
			
			
			
		}

	}
}


void updateVx(domain* D){
	
	
	int N;float sigma;
	N=D->N;
    
	for(int j=0;j<D->ny;j++){
		for(int k=0;k<D->nz;k++){
			int node;
			node=j*D->nx+k*D->nx*D->ny;
			
			for(int i=0;i<D->nx;i++){
                
				D->DerR[i]=D->px[node+i]+D->py[node+i]+D->pz[i+node];
				
			}
			
			for(int i=0;i<N-D->nx;i++){
                
				D->DerR[i+D->nx]=0.0;

            }
			
			executeDerShift2(D,D->Der,D->planfordDer,D->planbackDer,N);
       
			
			for(int i=0;i<D->nx;i++){
				if(D->s.EUx[i+node]!=40){
                        sigma=D->s.EUx[i+node]*D->sigma/(D->nPML);
                        
						
					D->vx[i+node]=(D->vx[i+node]-D->a2*D->DerR[i])/(1+sigma);	
				}
				else
					D->vx[i+node]=0.0;
			}
			
			
			
		}

	}
}


void updatePy(domain* D){
	
	
	int N;float sigma;
	N=D->N;
	for(int k=0;k<D->nz;k++){
		for(int i=0;i<D->nx;i++){
			int node;
			node=i+k*D->nx*D->ny;
			
			
			for(int j=0;j<D->ny;j++){
				D->Der2R[j]=D->vy[j*D->nx+node];
				
			}
			
			for(int j=0;j<N-D->ny;j++)
				D->Der2R[j+D->ny]=0.0;

			
			executeDerShift(D,D->Der2,D->planfordDer2,D->planbackDer2,N);
			
			
			for(int j=0;j<D->ny;j++){
				if(D->s.EUy[j*D->nx+node]!=40){
					sigma=D->s.EUy[j*D->nx+node]*D->sigma/(D->nPML);
					D->pyn[j*D->nx+node]=(1-sigma)*D->py[j*D->nx+node]-D->vel[j*D->nx+node]*D->vel[j*D->nx+node]*D->a1*D->Der2R[j];
				}
				else
					D->pyn[j*D->nx+node]=0.0;
			}
		}

	}
}

void updateVy(domain* D){
	
	
	int N;float sigma;
	N=D->N;
	for(int k=0;k<D->nz;k++){
		for(int i=0;i<D->nx;i++){
			int node;
			node=i+k*D->nx*D->ny;
			
			
			for(int j=0;j<D->ny;j++){
				D->Der2R[j]=D->px[j*D->nx+node]+D->py[j*D->nx+node]+D->pz[j*D->nx+node];
			}
			
			for(int j=0;j<N-D->ny;j++)
				D->Der2R[j+D->ny]=0.0;

			
			executeDerShift2(D,D->Der2,D->planfordDer2,D->planbackDer2,N);

			
			for(int j=0;j<D->ny;j++){
				if(D->s.EUy[j*D->nx+node]!=40){
                    
                    sigma=D->s.EUy[j*D->nx+node]*D->sigma/(D->nPML);
					D->vy[j*D->nx+node]=(D->vy[j*D->nx+node]-D->a2*D->Der2R[j])/(1+sigma);	
					
				}
				else
					D->vy[j*D->nx+node]=0.0;
			}
		}

	}
}


void updatePz(domain* D){
	
	
	int N;float sigma;
	N=D->N;
	for(int j=0;j<D->ny;j++){
		for(int i=0;i<D->nx;i++){
			int node;
			node=i+j*D->nx;
			
			for(int k=0;k<D->nz;k++){
				D->Der3R[k]=D->vz[k*D->nx*D->ny+node];
			}
			
			for(int k=0;k<N-D->nz;k++)
				D->Der3R[k+D->nz]=0.0;

			executeDerShift(D,D->Der3,D->planfordDer3,D->planbackDer3,N);
			
			
			for(int k=0;k<D->nz;k++){
				if(D->s.EUz[k*D->nx*D->ny+node]!=40){
					sigma=D->s.EUz[k*D->nx*D->ny+node]*D->sigma/(D->nPML);	
					D->pzn[k*D->nx*D->ny+node]=(1-sigma)*D->pz[k*D->nx*D->ny+node]-D->vel[k*D->nx*D->ny+node]*D->vel[k*D->nx*D->ny+node]*D->a1*D->Der3R[k];
								
				}
				else
					D->pzn[k*D->nx*D->ny+node]=0.0;
			}
			
			
			
		}

	}
}

void updateVz(domain* D){
	
	
	int N;float sigma;
	N=D->N;
	for(int j=0;j<D->ny;j++){
		for(int i=0;i<D->nx;i++){
			int node;
			node=i+j*D->nx;
			
			for(int k=0;k<D->nz;k++){
				D->Der3R[k]=D->px[k*D->nx*D->ny+node]+D->py[k*D->nx*D->ny+node]+D->pz[k*D->nx*D->ny+node];
			}
			
			for(int k=0;k<N-D->nz;k++)
				D->Der3R[k+D->nz]=0.0;

			executeDerShift2(D,D->Der3,D->planfordDer3,D->planbackDer3,N);

			
			for(int k=0;k<D->nz;k++){
				if(D->s.EUz[k*D->nx*D->ny+node]!=40){
                    
                    sigma=D->s.EUz[k*D->nx*D->ny+node]*D->sigma/(D->nPML);
                    
					D->vz[k*D->nx*D->ny+node]=(D->vz[k*D->nx*D->ny+node]-D->a2*D->Der3R[k])/(1+sigma);	
					
					
				}
				else
					D->vz[k*D->nx*D->ny+node]=0.0;
			}
			
			
			
		}

	}
}

void tickEuler(domain* D){

	for (int i=0;i<D->nx*D->ny*D->nz;i++){
		D->px[i]=D->pxn[i];
		D->py[i]=D->pyn[i];
		D->pz[i]=D->pzn[i];
		D->pxn[i]=0;
		D->pyn[i]=0;
		D->pzn[i]=0;
	}


}











void sourceDamped(domain* D,int type,int nS,int t){

	if(type==0){

		for(int i=-1;i<2;i++){
			for(int j=-1;j<2;j++){
				for(int k=-1;k<2;k++){
					D->pxn[nS+i+j*D->nx+k*D->nx*D->ny]=D->source[t]*exp(-1*((i)*(i)+(j)*(j)+(k)*(k)))+D->pxn[nS+i+j*D->nx+k*D->nx*D->ny];
					D->pyn[nS+i+j*D->nx+k*D->nx*D->ny]=D->source[t]*exp(-1*((i)*(i)+(j)*(j)+(k)*(k)))+D->pyn[nS+i+j*D->nx+k*D->nx*D->ny];
					D->pzn[nS+i+j*D->nx+k*D->nx*D->ny]=D->source[t]*exp(-1*((i)*(i)+(j)*(j)+(k)*(k)))+D->pzn[nS+i+j*D->nx+k*D->nx*D->ny];
					
				}
			}
		}
		
	}
	
}
void initParamDamped(domain* D,float Lx,float Ly,float Lz,float Dx,float Dt, float S, float tiempo,int nPML,float sigma,float cerjan){
	D->c=2000;float rho=1;
	D->a1=D->c*D->c*Dt*Dt;
	D->a2=0.0;
	D->S=D->c*Dt/Dx;
	D->Dx=Dx;
	D->nx=251+nPML;
	D->ny=251+nPML;
	D->nz=122+nPML;
	D->sigma=sigma/Dt;
	D->nPML=nPML/2;
	D->nt=(int)(tiempo/Dt);
	D->Dt=Dt;	
	D->cerjan=cerjan;
}

void initArraysDamped(domain* D){
	int nodes;
   
	nodes=D->nx*D->ny*D->nz;	
	D->px=(float*)malloc((nodes)*sizeof(float));
	D->py=(float*)malloc((nodes)*sizeof(float));
	D->pz=(float*)malloc((nodes)*sizeof(float));
	D->pxa=(float*)malloc((nodes)*sizeof(float));
	D->pya=(float*)malloc((nodes)*sizeof(float));
	D->pza=(float*)malloc((nodes)*sizeof(float));
	D->pxn=(float*)malloc((nodes)*sizeof(float));
	D->pyn=(float*)malloc((nodes)*sizeof(float));
	D->pzn=(float*)malloc((nodes)*sizeof(float));
	D->vel=(float*)malloc((nodes)*sizeof(float));
	D->s.EUx=(int*)malloc((nodes)*sizeof(int));
	D->s.EUy=(int*)malloc((nodes)*sizeof(int));
	D->s.EUz=(int*)malloc((nodes)*sizeof(int));
 
	for(int i=0;i<nodes;i++){
		D->px[i]=D->py[i]=D->pz[i]=D->pxa[i]=D->pya[i]=D->pza[i]=D->pxn[i]=D->pyn[i]=D->pzn[i]=0.0;
		D->s.EUx[i]=D->s.EUy[i]=D->s.EUz[i]=0;
	}
}


void defineDomainDamped(domain* D,int malla){
   
	for(int i=0;i<D->nx;i++){
		for(int j=0;j<D->ny;j++){
			for(int k=0;k<D->nz;k++){
				int node;
				node=i+j*D->nx+k*D->nx*D->ny;
				D->vel[node]=D->c;
				if(i<D->nPML)
					D->s.EUx[node]=D->nPML-i;
				if(i>D->nx-1-D->nPML)
					D->s.EUx[node]=i-(D->nx-1-D->nPML);
				if(j<D->nPML)
					D->s.EUy[node]=D->nPML-j;
				if(j>D->ny-1-D->nPML)
					D->s.EUy[node]=j-(D->ny-1-D->nPML);
				if(k<D->nPML){
					D->s.EUz[node]=D->nPML-k;}
				if(k>D->nz-1-D->nPML)
					D->s.EUz[node]=k-(D->nz-1-D->nPML);

				if(i==0){//0
					D->s.EUx[node]=40;D->s.EUy[node]=40;D->s.EUz[node]=40;}
				if(i==D->nx-1){
					D->s.EUx[node]=40;D->s.EUy[node]=40;D->s.EUz[node]=40;}	
				if(j==0){
					D->s.EUx[node]=40;D->s.EUy[node]=40;D->s.EUz[node]=40;}
				if(j==D->ny-1){
					D->s.EUx[node]=40;D->s.EUy[node]=40;D->s.EUz[node]=40;}
				if(k==0){
					D->s.EUx[node]=40;D->s.EUy[node]=40;D->s.EUz[node]=40;}
				if(k==D->nz-1){
					D->s.EUx[node]=40;D->s.EUy[node]=40;D->s.EUz[node]=40;}
				
			}
		}
	}




  int nData;FILE *r1;double av;
    
    if(malla==1)
        r1=fopen("./Mesh/Damped/3d.outAcc1.bin","rb");
    if(malla==2)
        r1=fopen("./Mesh/Damped/3d.outAcc2.bin","rb");
    if(malla==3)
        r1=fopen("./Mesh/Damped/3d.outAcc3.bin","rb");
    if(malla==4)
        r1=fopen("./Mesh/Damped/3d.outAcc4.bin","rb");
    if(malla==5)
        r1=fopen("./Mesh/Damped/3d.outAcc5.bin","rb");
    if(malla==6)
        r1=fopen("./Mesh/Damped/3d.outAcc6.bin","rb");
    if(malla==7)
        r1=fopen("./Mesh/Damped/3d.outAcc7.bin","rb");
    
    
	for(int i=0;i<D->nx;i++){
		for(int j=0;j<D->ny;j++){
			for(int k=0;k<D->nz;k++){
				int node;
				node=i+j*D->nx+k*D->nx*D->ny;
				nData=fread(&av,sizeof(double),1,r1);
                D->vel[node]=(float)(av);

			}
		}
	}
	fclose(r1);


}



void initDerivativesDamped(domain* D){
 
	float pi=3.14159265358979323846;
	int N,p2,N_2;
	N=D->nx;
	p2=2;
	D->nFFT=0;
	if(N<D->ny)
		N=D->ny;
	if(N<D->nz)
		N=D->nz;
	if(N%2!=0)
        N=N+1;
    N_2=N/2+1;
	D->CteDer=4*pi*pi/(N*N*D->Dx*D->Dx*N);
	D->Der=(fftwf_complex *)fftw_malloc(N_2*sizeof(fftwf_complex));
    D->DerR=(float *)malloc(N*sizeof(float));
	D->planfordDer = fftwf_plan_dft_r2c_1d(N, D->DerR, D->Der, FFTW_ESTIMATE);
	D->planbackDer = fftwf_plan_dft_c2r_1d(N, D->Der, D->DerR, FFTW_ESTIMATE);
	
	D->Der2=(fftwf_complex *)fftw_malloc(N_2*sizeof(fftwf_complex));
    D->Der2R=(float *)malloc(N*sizeof(float));
	D->planfordDer2 = fftwf_plan_dft_r2c_1d(N, D->Der2R, D->Der2, FFTW_ESTIMATE);
	D->planbackDer2 = fftwf_plan_dft_c2r_1d(N, D->Der2, D->Der2R, FFTW_ESTIMATE);
    
    D->Der3=(fftwf_complex *)fftw_malloc(N_2*sizeof(fftwf_complex));
    D->Der3R=(float *)malloc(N*sizeof(float));
	D->planfordDer3 = fftwf_plan_dft_r2c_1d(N, D->Der3R, D->Der3, FFTW_ESTIMATE);
	D->planbackDer3 = fftwf_plan_dft_c2r_1d(N, D->Der3, D->Der3R, FFTW_ESTIMATE);

	for(int i=0;i<N_2;i++){
		D->Der[i][0]=D->Der[i][1]=0.0;D->Der2[i][0]=D->Der2[i][1]=0.0;D->Der3[i][0]=D->Der3[i][1]=0.0;
	}

    
}

void executeDerDamped(domain* D,fftwf_complex* p,fftwf_plan planfordDer,fftwf_plan planbackDer,int N){

	fftwf_execute(planfordDer);
	for(int i=0;i<N/2+1;i++){
		float aux1,aux2;
		
        aux1=-D->CteDer*i*i*p[i][0];
        aux2=-D->CteDer*i*i*p[i][1];
    
		p[i][0]=aux1;
		p[i][1]=aux2;
	}

	fftwf_execute(planbackDer);
}

void updatePxDamped(domain* D){
	
	
	int N;float sigma;
	N=D->nx;
	if(N<D->ny)
		N=D->ny;
	else if(N<D->nz)
		N=D->nz;
    if(N%2!=0)
        N=N+1;
	for(int j=0;j<D->ny;j++){
		for(int k=0;k<D->nz;k++){
			int node;
			node=j*D->nx+k*D->nx*D->ny;
			
			for(int i=0;i<D->nx;i++){
				D->DerR[i]=D->px[node+i]+D->py[node+i]+D->pz[i+node];
				
			}
			
			for(int i=0;i<N-D->nx;i++)
				D->DerR[i+D->nx]=0.0;

			
			executeDerDamped(D,D->Der,D->planfordDer,D->planbackDer,N);

			
			for(int i=0;i<D->nx;i++){
				if(D->s.EUx[i+node]!=40){
					sigma=sqrt(D->s.EUx[i+node]*D->s.EUx[i+node]+D->s.EUy[i+node]*D->s.EUy[i+node]+D->s.EUz[i+node]*D->s.EUz[i+node])*D->sigma/D->nPML;
					D->pxn[i+node]=(4*D->px[i+node]/(2+sigma*D->Dt)+2*D->vel[i+node]*D->vel[i+node]*D->Dt*D->Dt*D->DerR[i]/(2+sigma*D->Dt)
							+(sigma*D->Dt-2)*D->pxa[i+node]/(2+sigma*D->Dt));
					
				}
				else
					D->pxn[i+node]=0.0;
			}
			
			
			
		}

	}
}


void updatePyDamped(domain* D){
	
	
	int N;float sigma;
	N=D->nx;
	if(N<D->ny)
		N=D->ny;
	else if(N<D->nz)
		N=D->nz;
    if(N%2!=0)
        N=N+1;
	for(int k=0;k<D->nz;k++){
		for(int i=0;i<D->nx;i++){
			int node;
			node=i+k*D->nx*D->ny;
			
			
			for(int j=0;j<D->ny;j++){
				D->Der2R[j]=D->px[j*D->nx+node]+D->py[j*D->nx+node]+D->pz[j*D->nx+node];
				
			}
			
			for(int j=0;j<N-D->ny;j++)
				D->Der2R[j+D->ny]=0.0;

			
			executeDerDamped(D,D->Der2,D->planfordDer2,D->planbackDer2,N);

			
			for(int j=0;j<D->ny;j++){
				if(D->s.EUy[j*D->nx+node]!=40){
					sigma=sqrt(D->s.EUx[j*D->nx+node]*D->s.EUx[j*D->nx+node]
						+D->s.EUy[j*D->nx+node]*D->s.EUy[j*D->nx+node]+D->s.EUz[j*D->nx+node]*D->s.EUz[j*D->nx+node])*D->sigma/D->nPML;
					D->pyn[j*D->nx+node]=(4*D->py[j*D->nx+node]/(2+sigma*D->Dt)+2*D->vel[j*D->nx+node]*D->vel[j*D->nx+node]*D->Dt*D->Dt*D->Der2R[j]/(2+sigma*D->Dt)
							+(sigma*D->Dt-2)*D->pya[j*D->nx+node]/(2+sigma*D->Dt));
				}
				else
					D->pyn[j*D->nx+node]=0.0;
			}
		}

	}
}


void updatePzDamped(domain* D){
	
	
	int N;float sigma;
	N=D->nx;
	if(N<D->ny)
		N=D->ny;
	else if(N<D->nz)
		N=D->nz;
    if(N%2!=0)
        N=N+1;
	for(int j=0;j<D->ny;j++){
		for(int i=0;i<D->nx;i++){
			int node;
			node=i+j*D->nx;
			
			for(int k=0;k<D->nz;k++){
				D->Der3R[k]=D->px[k*D->nx*D->ny+node]+D->py[k*D->nx*D->ny+node]+D->pz[k*D->nx*D->ny+node];
				
			}
			
			for(int k=0;k<N-D->nz;k++)
				D->Der3R[k+D->nz]=0.0;

			executeDerDamped(D,D->Der3,D->planfordDer3,D->planbackDer3,N);

			
			for(int k=0;k<D->nz;k++){
				if(D->s.EUz[k*D->nx*D->ny+node]!=40){
					sigma=sqrt(D->s.EUx[k*D->nx*D->ny+node]*D->s.EUx[k*D->nx*D->ny+node]
						+D->s.EUy[k*D->nx*D->ny+node]*D->s.EUy[k*D->nx*D->ny+node]+D->s.EUz[k*D->nx*D->ny+node]*D->s.EUz[k*D->nx*D->ny+node])*D->sigma/D->nPML;
					D->pzn[k*D->nx*D->ny+node]=(4*D->pz[k*D->nx*D->ny+node]/(2+sigma*D->Dt)
								+2*D->vel[k*D->nx*D->ny+node]*D->vel[k*D->nx*D->ny+node]*D->Dt*D->Dt*D->Der3R[k]/(2+sigma*D->Dt)
							+(sigma*D->Dt-2)*D->pza[k*D->nx*D->ny+node]/(2+sigma*D->Dt));
				}
				else
					D->pzn[k*D->nx*D->ny+node]=0.0;
			}
			
			
			
		}

	}
}



void tickDamped(domain* D){

	for (int i=0;i<D->nx*D->ny*D->nz;i++){
		D->pxa[i]=D->px[i];
		D->px[i]=D->pxn[i];
		D->pya[i]=D->py[i];
		D->py[i]=D->pyn[i];
		D->pza[i]=D->pz[i];
		D->pz[i]=D->pzn[i];
		D->pxn[i]=0;
		D->pyn[i]=0;
		D->pzn[i]=0;
	}


}
void freeDamped(domain* D){
    
 free(D->s.EUx);free(D->s.EUy);free(D->s.EUz);free(D->vel);
 free(D->pxn);free(D->pyn);free(D->pzn);
 free(D->px);free(D->py);free(D->pz);
 free(D->pxa);free(D->pya);free(D->pza);
 free(D->DerR);free(D->Der2R);free(D->Der3R);
 fftwf_free(D->Der);fftwf_free(D->Der2);fftwf_free(D->Der3);
}
void freeEuler(domain* D){
    
 free(D->s.EUx);free(D->s.EUy);free(D->s.EUz);free(D->vel);
 free(D->pxn);free(D->pyn);free(D->pzn);
 free(D->px);free(D->py);free(D->pz);
 free(D->vx);free(D->vy);free(D->vz);
 free(D->DerR);free(D->Der2R);free(D->Der3R);
 fftwf_free(D->Der);fftwf_free(D->Der2);fftwf_free(D->Der3);
}




void freeSponge(domain* D){
    
 free(D->s.EUx);free(D->s.EUy);free(D->s.EUz);free(D->vel);
 free(D->pxn);free(D->pyn);free(D->pzn);
 free(D->dpxn);free(D->dpyn);free(D->dpzn);
 free(D->pxa);free(D->pya);free(D->pza);
 free(D->dpxa);free(D->dpya);free(D->dpza);
 free(D->DerR);free(D->Der2R);free(D->Der3R);
 fftwf_free(D->Der);fftwf_free(D->Der2);fftwf_free(D->Der3);
}

void sourceSponge(domain* D,int type,int nS,int t){

        float aux=0.0;
        for( int l=0;l<t;l++)
            aux=aux+D->source[l];
		for(int i=-1;i<2;i++){
			for(int j=-1;j<2;j++){
				for(int k=-1;k<2;k++){
                    if(t!=0){
                        D->pxn[nS+i+j*D->nx+k*D->nx*D->ny]=aux*exp(-1*((i)*(i)+(j)*(j)+(k)*(k)))+D->pxn[nS+i+j*D->nx+k*D->nx*D->ny];
                        D->pyn[nS+i+j*D->nx+k*D->nx*D->ny]=aux*exp(-1*((i)*(i)+(j)*(j)+(k)*(k)))+D->pyn[nS+i+j*D->nx+k*D->nx*D->ny];
                        D->pzn[nS+i+j*D->nx+k*D->nx*D->ny]=aux*exp(-1*((i)*(i)+(j)*(j)+(k)*(k)))+D->pzn[nS+i+j*D->nx+k*D->nx*D->ny];
                    }
				}
			}
		}
		
	
	
}
void initParamSponge(domain* D,float Lx,float Ly,float Lz,float Dx,float Dt, float S, float tiempo,int nPML,float sigma,float cerjan){
	D->c=2000;float rho=1;
	D->a1=D->c*D->c*Dt;
	D->a2=Dt;
	D->S=S;
	D->Dx=Dx;
	D->nx=251+nPML;
	D->ny=251+nPML;
	D->nz=122+nPML;
	D->sigma=sigma;
	D->nPML=nPML/2;
	D->nt=(int)(tiempo/Dt);
	D->Dt=Dt;	
	D->cerjan=cerjan;
}

void initArraysSponge(domain* D){
	int nodes;
	nodes=D->nx*D->ny*D->nz;
	D->dpxn=(float*)malloc((nodes)*sizeof(float));
	D->dpyn=(float*)malloc((nodes)*sizeof(float));
	D->dpzn=(float*)malloc((nodes)*sizeof(float));	
	D->dpxa=(float*)malloc((nodes)*sizeof(float));
	D->dpya=(float*)malloc((nodes)*sizeof(float));
	D->dpza=(float*)malloc((nodes)*sizeof(float));
	D->pxa=(float*)malloc((nodes)*sizeof(float));
	D->pya=(float*)malloc((nodes)*sizeof(float));
	D->pza=(float*)malloc((nodes)*sizeof(float));
	D->pxn=(float*)malloc((nodes)*sizeof(float));
	D->pyn=(float*)malloc((nodes)*sizeof(float));
	D->pzn=(float*)malloc((nodes)*sizeof(float));
	D->vel=(float*)malloc((nodes)*sizeof(float));
	D->s.EUx=(int*)malloc((nodes)*sizeof(int));
	D->s.EUy=(int*)malloc((nodes)*sizeof(int));
	D->s.EUz=(int*)malloc((nodes)*sizeof(int));

	for(int i=0;i<nodes;i++){
		D->dpxn[i]=D->dpyn[i]=D->dpzn[i]=D->dpxa[i]=D->dpya[i]=D->dpza[i]=D->pxa[i]=D->pya[i]=D->pza[i]=D->pxn[i]=D->pyn[i]=D->pzn[i]=0.0;
		D->s.EUx[i]=D->s.EUy[i]=D->s.EUz[i]=0;
	}
}


void defineDomainSponge(domain* D,int malla){

	for(int i=0;i<D->nx;i++){
		for(int j=0;j<D->ny;j++){
			for(int k=0;k<D->nz;k++){
				int node;
				node=i+j*D->nx+k*D->nx*D->ny;
				D->vel[node]=D->c;
				if(i<D->nPML)
					D->s.EUx[node]=D->nPML-i;
				if(i>D->nx-1-D->nPML)
					D->s.EUx[node]=i-(D->nx-1-D->nPML);
				if(j<D->nPML)
					D->s.EUy[node]=D->nPML-j;
				if(j>D->ny-1-D->nPML)
					D->s.EUy[node]=j-(D->ny-1-D->nPML);
				if(k<D->nPML){
					D->s.EUz[node]=D->nPML-k;}
				if(k>D->nz-1-D->nPML)
					D->s.EUz[node]=k-(D->nz-1-D->nPML);

				if(i==0){//0
					D->s.EUx[node]=40;D->s.EUy[node]=40;D->s.EUz[node]=40;}
				if(i==D->nx-1){
					D->s.EUx[node]=40;D->s.EUy[node]=40;D->s.EUz[node]=40;}	
				if(j==0){
					D->s.EUx[node]=40;D->s.EUy[node]=40;D->s.EUz[node]=40;}
				if(j==D->ny-1){
					D->s.EUx[node]=40;D->s.EUy[node]=40;D->s.EUz[node]=40;}
				if(k==0){
					D->s.EUx[node]=40;D->s.EUy[node]=40;D->s.EUz[node]=40;}
				if(k==D->nz-1){
					D->s.EUx[node]=40;D->s.EUy[node]=40;D->s.EUz[node]=40;}
				
				
			}
		}
	}

	int nData;FILE *r1;double av;
    
    if(malla==1)
        r1=fopen("./Mesh/Sponge/3d.outAcc1.bin","rb");
    if(malla==2)
        r1=fopen("./Mesh/Sponge/3d.outAcc2.bin","rb");
    if(malla==3)
        r1=fopen("./Mesh/Sponge/3d.outAcc3.bin","rb");
    if(malla==4)
        r1=fopen("./Mesh/Sponge/3d.outAcc4.bin","rb");
    if(malla==5)
        r1=fopen("./Mesh/Sponge/3d.outAcc5.bin","rb");
    if(malla==6)
        r1=fopen("./Mesh/Sponge/3d.outAcc6.bin","rb");
    if(malla==7)
        r1=fopen("./Mesh/Sponge/3d.outAcc7.bin","rb");
    
    
	for(int i=0;i<D->nx;i++){
		for(int j=0;j<D->ny;j++){
			for(int k=0;k<D->nz;k++){
				int node;
               
				node=i+j*D->nx+k*D->nx*D->ny;
				
				nData=fread(&av,sizeof(double),1,r1);
                D->vel[node]=(float)(av);
                

			}
		}
	}
	 fclose(r1);
     
	
}



void initDerivativesSponge(domain* D){

	float pi=3.14159265358979323846;
	int N,p2,N_2;
	N=D->nx;
	p2=2;
	D->nFFT=0;
	if(N<D->ny)
		N=D->ny;
	if(N<D->nz)
		N=D->nz;
	if(N%2!=0)
        N=N+1;
    D->N=N;
	N=N+D->nFFT;
    N_2=N/2+1;
     
	D->CteDer=4*pi*pi/(N*N*D->Dx*D->Dx*N);
	D->Der=(fftwf_complex *)fftwf_malloc(N_2*sizeof(fftwf_complex));
    D->DerR=(float *)malloc(N*sizeof(float));
	D->planfordDer = fftwf_plan_dft_r2c_1d(N, D->DerR, D->Der, FFTW_ESTIMATE);
	D->planbackDer = fftwf_plan_dft_c2r_1d(N, D->Der, D->DerR, FFTW_ESTIMATE);

	D->Der2=(fftwf_complex *)fftwf_malloc(N_2*sizeof(fftwf_complex));
    D->Der2R=(float *)malloc(N*sizeof(float));
	D->planfordDer2 = fftwf_plan_dft_r2c_1d(N, D->Der2R, D->Der2, FFTW_ESTIMATE);
	D->planbackDer2 = fftwf_plan_dft_c2r_1d(N, D->Der2, D->Der2R, FFTW_ESTIMATE);

	D->Der3=(fftwf_complex *)fftwf_malloc(N_2*sizeof(fftwf_complex));
    D->Der3R=(float *)malloc(N*sizeof(float));
	D->planfordDer3 = fftwf_plan_dft_r2c_1d(N, D->Der3R, D->Der3, FFTW_ESTIMATE);
	D->planbackDer3 = fftwf_plan_dft_c2r_1d(N, D->Der3, D->Der3R, FFTW_ESTIMATE);
	
	for(int i=0;i<N;i++){
        D->DerR[i]=D->Der2R[i]=0.0;D->Der3R[i]=0.0;
        if(i<N_2){
            D->Der[i][0]=D->Der[i][1]=0.0;D->Der2[i][0]=D->Der2[i][1]=0.0;D->Der3[i][0]=D->Der3[i][1]=0.0;
        }
	}
}

void executeDerSponge(domain* D,fftwf_complex* p,fftwf_plan planfordDer,fftwf_plan planbackDer,int N){

	fftwf_execute(planfordDer);
	for(int i=0;i<N/2+1;i++){
		float aux1,aux2;
        aux1=-D->CteDer*i*i*p[i][0];
        aux2=-D->CteDer*i*i*p[i][1];
		
		p[i][0]=aux1;
		p[i][1]=aux2;
	}

	fftwf_execute(planbackDer);
}

void updatePxSponge(domain* D){
	
	
	int N;float c,a1;float dist;
	N=D->N;
	
	for(int j=0;j<D->ny;j++){
		for(int k=0;k<D->nz;k++){
			int node;
			node=j*D->nx+k*D->nx*D->ny;
			
			for(int i=0;i<N;i++){
                if(i<D->nx)
                    D->DerR[i]=D->pxa[node+i]+D->pya[node+i]+D->pza[i+node];
                else
                    D->DerR[i]=0.0;
			}
			
			
			
			executeDerSponge(D,D->Der,D->planfordDer,D->planbackDer,N);

			
			for(int i=0;i<D->nx;i++){
				if(D->s.EUx[i+node]!=40){
					dist=sqrt(D->s.EUx[i+node]*D->s.EUx[i+node]+D->s.EUy[i+node]*D->s.EUy[i+node]+D->s.EUz[i+node]*D->s.EUz[i+node]);
					a1=D->vel[i+node]*D->vel[i+node]*D->Dt;
					c=exp(-D->cerjan*D->cerjan*dist*dist);
					D->dpxn[i+node]=(D->dpxa[i+node]+a1*D->DerR[i]);
					D->pxn[i+node]=(D->pxa[i+node]+D->a2*D->dpxn[i+node]);
					D->dpxn[i+node]=c*(D->dpxn[i+node]);
					D->pxn[i+node]=c*(D->pxn[i+node]);
					
				}
				else
					D->pxn[i+node]=D->dpxn[i+node]=0.0;
			}
			
			
			
		}

	}
}


void updatePySponge(domain* D){
	
	
	int N;float c,a1;float dist;
	N=D->N;
    
	for(int k=0;k<D->nz;k++){
		for(int i=0;i<D->nx;i++){
			int node;
			node=i+k*D->nx*D->ny;
			
			
			for(int j=0;j<N;j++){
                if(i<D->ny)
                    D->Der2R[j]=D->pxa[j*D->nx+node]+D->pya[j*D->nx+node]+D->pza[j*D->nx+node];
                else
                    D->Der2R[j]=0.0;
			}
			
			executeDerSponge(D,D->Der2,D->planfordDer2,D->planbackDer2,N);

			
			for(int j=0;j<D->ny;j++){
				if(D->s.EUy[j*D->nx+node]!=40){
					dist=sqrt(D->s.EUx[j*D->nx+node]*D->s.EUx[j*D->nx+node]+D->s.EUy[j*D->nx+node]*D->s.EUy[j*D->nx+node]+D->s.EUz[j*D->nx+node]*D->s.EUz[j*D->nx+node]);
					a1=D->vel[j*D->nx+node]*D->vel[j*D->nx+node]*D->Dt;
					c=exp(-D->cerjan*D->cerjan*dist*dist);
					D->dpyn[j*D->nx+node]=(D->dpya[j*D->nx+node]+a1*D->Der2R[j]);
					D->pyn[j*D->nx+node]=(D->pya[j*D->nx+node]+D->a2*D->dpyn[j*D->nx+node]);
					D->dpyn[j*D->nx+node]=c*(D->dpyn[j*D->nx+node]);
					D->pyn[j*D->nx+node]=c*(D->pyn[j*D->nx+node]);
				}
				else
					D->dpyn[j*D->nx+node]=D->pyn[j*D->nx+node]=0.0;
			}
		}

	}
}


void updatePzSponge(domain* D){
	
	
	int N;float c,a1;float dist;
	N=D->N;

	for(int j=0;j<D->ny;j++){
		for(int i=0;i<D->nx;i++){
			int node;
			node=i+j*D->nx;
			
			for(int k=0;k<N;k++){
                if(k<D->nz)
                    D->Der3R[k]=D->pxa[k*D->nx*D->ny+node]+D->pya[k*D->nx*D->ny+node]+D->pza[k*D->nx*D->ny+node];
                else
                    D->Der3R[k]=0.0;
			}
			
			executeDerSponge(D,D->Der3,D->planfordDer3,D->planbackDer3,N);

			
			for(int k=0;k<D->nz;k++){
				if(D->s.EUz[k*D->nx*D->ny+node]!=40){
					dist=sqrt(D->s.EUx[k*D->nx*D->ny+node]*D->s.EUx[k*D->nx*D->ny+node]
					+D->s.EUy[k*D->nx*D->ny+node]*D->s.EUy[k*D->nx*D->ny+node]+D->s.EUz[k*D->nx*D->ny+node]*D->s.EUz[k*D->nx*D->ny+node]);
					a1=D->vel[k*D->nx*D->ny+node]*D->vel[k*D->nx*D->ny+node]*D->Dt;
					c=exp(-D->cerjan*D->cerjan*dist*dist);
					D->dpzn[k*D->nx*D->ny+node]=(D->dpza[k*D->nx*D->ny+node]+a1*D->Der3R[k]);
					D->pzn[k*D->nx*D->ny+node]=(D->pza[k*D->nx*D->ny+node]+D->a2*D->dpzn[k*D->nx*D->ny+node]);
					D->dpzn[k*D->nx*D->ny+node]=c*(D->dpzn[k*D->nx*D->ny+node]);
					D->pzn[k*D->nx*D->ny+node]=c*(D->pzn[k*D->nx*D->ny+node]);
				}
				else
					D->dpzn[k*D->nx*D->ny+node]=D->pzn[k*D->nx*D->ny+node]=0.0;
			}
			
			
			
		}

	}
}



void tickSponge(domain* D){

	for (int i=0;i<D->nx*D->ny*D->nz;i++){
		D->dpxa[i]=D->dpxn[i];
		D->pxa[i]=D->pxn[i];
		D->dpya[i]=D->dpyn[i];
		D->pya[i]=D->pyn[i];
		D->dpza[i]=D->dpzn[i];
		D->pza[i]=D->pzn[i];
		
		D->pxn[i]=0;
		D->pyn[i]=0;
		D->pzn[i]=0;
		D->dpxn[i]=0;
		D->dpyn[i]=0;
		D->dpzn[i]=0;
	}


}

