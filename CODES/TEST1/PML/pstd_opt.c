#include "pstd_opt.h"



void ricker_wavelet_tis0(domain* D,double dt,int nt,double t0,double f0,int l){

	D->source=(double*)malloc((D->nt)*sizeof(double));
	double B=0;double Pi;  Pi=3.141592654;
	B=Pi*Pi*f0*f0;
	if(l==0){
		for(int t=0;t<D->nt;t++){
				D->source[t]=100*((1-2*B*(t*D->Dt-t0)*(t*D->Dt-t0))*exp(-B*(t*D->Dt-t0)*(t*D->Dt-t0))); 
				
		}
	}

}

void source(domain* D,int type,int nS,int t){

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
void initParam(domain* D,double Lx,double Ly,double Lz,double Dx,double Dt, double S, double tiempo,int nPML,double sigma,double cerjan,int nB,int nParam){
	D->c=2000;double rho=1;//1.21;
	D->a1=rho*Dt;
	D->a2=Dt/rho;//Dt/(rho);
	D->S=D->c*Dt/Dx;
	D->Dx=Dx;//D->c*Dt/S;
	D->nx=(int)(Lx/D->Dx)+nPML;
	D->ny=(int)(Ly/D->Dx)+nPML;
	D->nz=(int)(Lz/D->Dx)+nPML;
	D->sigma=sigma;
	D->nPML=nPML/2;
	D->nt=(int)(tiempo/Dt);
	D->Dt=Dt;	
	D->cerjan=cerjan;
	D->Nb=nB;
	D->Nparam=nParam;
}

void initArrays(domain* D){
	int nodes;
	nodes=D->nx*D->ny*D->nz;	
	D->px=(double*)malloc((nodes)*sizeof(double));
	D->py=(double*)malloc((nodes)*sizeof(double));
	D->pz=(double*)malloc((nodes)*sizeof(double));
	D->vx=(double*)malloc((nodes)*sizeof(double));
	D->vy=(double*)malloc((nodes)*sizeof(double));
	D->vz=(double*)malloc((nodes)*sizeof(double));
	D->pxn=(double*)malloc((nodes)*sizeof(double));
	D->pyn=(double*)malloc((nodes)*sizeof(double));
	D->pzn=(double*)malloc((nodes)*sizeof(double));
	D->vel=(double*)malloc((nodes)*sizeof(double));
	D->s.EUx=(int*)malloc((nodes)*sizeof(int));
	D->s.EUy=(int*)malloc((nodes)*sizeof(int));
	D->s.EUz=(int*)malloc((nodes)*sizeof(int));

	for(int i=0;i<nodes;i++){
		D->px[i]=D->py[i]=D->pz[i]=D->vx[i]=D->vy[i]=D->vz[i]=D->pxn[i]=D->pyn[i]=D->pzn[i]=0.0;
		D->s.EUx[i]=D->s.EUy[i]=D->s.EUz[i]=0;
	}
}


void defineDomain(domain* D){

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
				if(k<D->nPML)
					D->s.EUz[node]=D->nPML-k;
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

	D->s.nLTx=D->s.nLTy=D->s.nLTz=2;
}



void initDerivatives(domain* D){

	double pi=3.14159265358979323846;
	int N,p2;
	N=D->nx;
	p2=2;
	D->nFFT=0;
	if(N<D->ny)
		N=D->ny;
	else if(N<D->nz)
		N=D->nz;
	
	N=N+D->nFFT;
	D->CteDer=2*pi/(N*D->Dx*N);
	D->Der=(fftw_complex *)fftw_malloc(N*sizeof(fftw_complex));
	D->planfordDer = fftw_plan_dft_1d(N, D->Der, D->Der, FFTW_FORWARD, FFTW_ESTIMATE);
	D->planbackDer = fftw_plan_dft_1d(N, D->Der, D->Der, FFTW_BACKWARD, FFTW_ESTIMATE);
	
	D->Der2=(fftw_complex *)fftw_malloc(N*sizeof(fftw_complex));
	D->planfordDer2 = fftw_plan_dft_1d(N, D->Der2, D->Der2, FFTW_FORWARD, FFTW_ESTIMATE);
	D->planbackDer2 = fftw_plan_dft_1d(N, D->Der2, D->Der2, FFTW_BACKWARD, FFTW_ESTIMATE);

	D->Der3=(fftw_complex *)fftw_malloc(N*sizeof(fftw_complex));
	D->planfordDer3 = fftw_plan_dft_1d(N, D->Der3, D->Der3, FFTW_FORWARD, FFTW_ESTIMATE);
	D->planbackDer3 = fftw_plan_dft_1d(N, D->Der3, D->Der3, FFTW_BACKWARD, FFTW_ESTIMATE);

	for(int i=0;i<N;i++){
		D->Der[i][0]=D->Der[i][1]=0.0;D->Der2[i][0]=D->Der2[i][1]=0.0;D->Der3[i][0]=D->Der3[i][1]=0.0;
	}
}


void executeDerShift(domain* D,fftw_complex* p,fftw_plan planfordDer,fftw_plan planbackDer,int N){
	double pi=3.14159265358979323846;
	fftw_execute(planfordDer);
	for(int i=0;i<N;i++){
		double aux1,aux2;
		if(i<N/2-1){
			aux1=(-D->CteDer*i*p[i][1]*cos(D->CteDer*i*N*D->Dx/2)-D->CteDer*i*p[i][0]*sin(D->CteDer*i*N*D->Dx/2));
			aux2=(-D->CteDer*i*p[i][1]*sin(D->CteDer*i*N*D->Dx/2)+D->CteDer*i*p[i][0]*cos(D->CteDer*i*N*D->Dx/2));
		}
		else if(i==N/2-1){
			aux1=aux2=0.0;
		}
		else{
			aux1=(-D->CteDer*(i-N)*p[i][1]*cos(D->CteDer*(i-N)*N*D->Dx/2)-D->CteDer*(i-N)*p[i][0]*sin(D->CteDer*(i-N)*N*D->Dx/2));
			aux2=(-D->CteDer*(i-N)*p[i][1]*sin(D->CteDer*(i-N)*N*D->Dx/2)+D->CteDer*(i-N)*p[i][0]*cos(D->CteDer*(i-N)*N*D->Dx/2));
			
		}
		p[i][0]=aux1;
		p[i][1]=aux2;
	}

	fftw_execute(planbackDer);
}


void executeDerShift2(domain* D,fftw_complex* p,fftw_plan planfordDer,fftw_plan planbackDer,int N){
	double pi=3.14159265358979323846;
	fftw_execute(planfordDer);
	for(int i=0;i<N;i++){
		double aux1,aux2;
		if(i<N/2-1){
			aux1=(-D->CteDer*i*p[i][1]*cos(D->CteDer*i*N*D->Dx/2)+D->CteDer*i*p[i][0]*sin(D->CteDer*i*N*D->Dx/2));
			aux2=(D->CteDer*i*p[i][1]*sin(D->CteDer*i*N*D->Dx/2)+D->CteDer*i*p[i][0]*cos(D->CteDer*i*N*D->Dx/2));
		}
		else if(i==N/2-1){
			aux1=aux2=0.0;
		}
		else{
			aux1=(-D->CteDer*(i-N)*p[i][1]*cos(D->CteDer*(i-N)*N*D->Dx/2)+D->CteDer*(i-N)*p[i][0]*sin(D->CteDer*(i-N)*N*D->Dx/2));
			aux2=(D->CteDer*(i-N)*p[i][1]*sin(D->CteDer*(i-N)*N*D->Dx/2)+D->CteDer*(i-N)*p[i][0]*cos(D->CteDer*(i-N)*N*D->Dx/2));
			
		}
		p[i][0]=aux1;
		p[i][1]=aux2;
	}

	fftw_execute(planbackDer);
}

void updatePx(domain* D){
	
	
	int N;double sigma;
	N=D->nx;
	if(N<D->ny)
		N=D->ny;
	else if(N<D->nz)
		N=D->nz;
	for(int j=0;j<D->ny;j++){
		for(int k=0;k<D->nz;k++){
			int node;
			node=j*D->nx+k*D->nx*D->ny;
			
			for(int i=0;i<D->nx;i++){
				D->Der[i][0]=D->vx[node+i];
				D->Der[i][1]=0.0;
			}
			
			for(int i=0;i<N-D->nx;i++)
				D->Der[i+D->nx][0]=D->Der[i+D->nx][1]=0.0;

			
			executeDerShift(D,D->Der,D->planfordDer,D->planbackDer,N);

			
			for(int i=0;i<D->nx;i++){
				if(D->s.EUx[i+node]!=40){
					sigma=D->s.EUx[i+node]*D->sigma/D->nPML;
					D->pxn[i+node]=(1-sigma)*D->px[i+node]-D->vel[i+node]*D->vel[i+node]*D->a1*D->Der[i][0];
				
				}
				else
					D->pxn[i+node]=0.0;
			}
			
			
			
		}

	}
}


void updateVx(domain* D){
	
	
	int N;double sigma;
	N=D->nx;
	if(N<D->ny)
		N=D->ny;
	else if(N<D->nz)
		N=D->nz;
	for(int j=0;j<D->ny;j++){
		for(int k=0;k<D->nz;k++){
			int node;
			node=j*D->nx+k*D->nx*D->ny;
			
			for(int i=0;i<D->nx;i++){
				D->Der[i][0]=D->px[node+i]+D->py[node+i]+D->pz[i+node];
				D->Der[i][1]=0.0;
			}
			
			for(int i=0;i<N-D->nx;i++)
				D->Der[i+D->nx][0]=D->Der[i+D->nx][1]=0.0;

			
			executeDerShift2(D,D->Der,D->planfordDer,D->planbackDer,N);

			
			for(int i=0;i<D->nx;i++){
				if(D->s.EUx[i+node]!=40){
					sigma=D->s.EUx[i+node]*D->sigma/D->nPML;
						
					D->vx[i+node]=(D->vx[i+node]-D->a2*D->Der[i][0])/(1+sigma);	
				}
				else
					D->vx[i+node]=0.0;
			}
			
			
			
		}

	}
}



void updatePy(domain* D){
	
	
	int N;double sigma;
	N=D->nx;
	if(N<D->ny)
		N=D->ny;
	else if(N<D->nz)
		N=D->nz;
	for(int k=0;k<D->nz;k++){
		for(int i=0;i<D->nx;i++){
			int node;
			node=i+k*D->nx*D->ny;
			
			
			for(int j=0;j<D->ny;j++){
				D->Der2[j][0]=D->vy[j*D->nx+node];
				D->Der2[j][1]=0.0;
			}
			
			for(int j=0;j<N-D->ny;j++)
				D->Der2[j+D->ny][0]=D->Der2[j+D->ny][1]=0.0;

			
			executeDerShift(D,D->Der2,D->planfordDer2,D->planbackDer2,N);

			
			for(int j=0;j<D->ny;j++){
				if(D->s.EUy[j*D->nx+node]!=40){
					sigma=D->s.EUy[j*D->nx+node]*D->sigma/D->nPML;
					D->pyn[j*D->nx+node]=(1-sigma)*D->py[j*D->nx+node]-D->vel[j*D->nx+node]*D->vel[j*D->nx+node]*D->a1*D->Der2[j][0];
				}
				else
					D->pyn[j*D->nx+node]=0.0;
			}
		}

	}
}

void updateVy(domain* D){
	
	
	int N;double sigma;
	N=D->nx;
	if(N<D->ny)
		N=D->ny;
	else if(N<D->nz)
		N=D->nz;
	for(int k=0;k<D->nz;k++){
		for(int i=0;i<D->nx;i++){
			int node;
			node=i+k*D->nx*D->ny;
			
			
			for(int j=0;j<D->ny;j++){
				D->Der2[j][0]=D->px[j*D->nx+node]+D->py[j*D->nx+node]+D->pz[j*D->nx+node];
				D->Der2[j][1]=0.0;
			}
			
			for(int j=0;j<N-D->ny;j++)
				D->Der2[j+D->ny][0]=D->Der2[j+D->ny][1]=0.0;

			
			executeDerShift2(D,D->Der2,D->planfordDer2,D->planbackDer2,N);

			
			for(int j=0;j<D->ny;j++){
				if(D->s.EUy[j*D->nx+node]!=40){
					sigma=D->s.EUy[j*D->nx+node]*D->sigma/D->nPML;
					D->vy[j*D->nx+node]=(D->vy[j*D->nx+node]-D->a2*D->Der2[j][0])/(1+sigma);	
					
				}
				else
					D->vy[j*D->nx+node]=0.0;
			}
		}

	}
}


void updatePz(domain* D){
	
	
	int N;double sigma;
	N=D->nx;
	if(N<D->ny)
		N=D->ny;
	else if(N<D->nz)
		N=D->nz;
	for(int j=0;j<D->ny;j++){
		for(int i=0;i<D->nx;i++){
			int node;
			node=i+j*D->nx;
			
			for(int k=0;k<D->nz;k++){
				D->Der3[k][0]=D->vz[k*D->nx*D->ny+node];
				D->Der3[k][1]=0.0;
			}
			
			for(int k=0;k<N-D->nz;k++)
				D->Der3[k+D->nz][0]=D->Der3[k+D->nz][1]=0.0;

			executeDerShift(D,D->Der3,D->planfordDer3,D->planbackDer3,N);

			
			for(int k=0;k<D->nz;k++){
				if(D->s.EUz[k*D->nx*D->ny+node]!=40){
					sigma=D->s.EUz[k*D->nx*D->ny+node]*D->sigma/D->nPML;	
					D->pzn[k*D->nx*D->ny+node]=(1-sigma)*D->pz[k*D->nx*D->ny+node]-D->vel[k*D->nx*D->ny+node]*D->vel[k*D->nx*D->ny+node]*D->a1*D->Der3[k][0];
								
				}
				else
					D->pzn[k*D->nx*D->ny+node]=0.0;
			}
			
			
			
		}

	}
}

void updateVz(domain* D){
	
	
	int N;double sigma;
	N=D->nx;
	if(N<D->ny)
		N=D->ny;
	else if(N<D->nz)
		N=D->nz;
	for(int j=0;j<D->ny;j++){
		for(int i=0;i<D->nx;i++){
			int node;
			node=i+j*D->nx;
			
			for(int k=0;k<D->nz;k++){
				D->Der3[k][0]=D->px[k*D->nx*D->ny+node]+D->py[k*D->nx*D->ny+node]+D->pz[k*D->nx*D->ny+node];
				D->Der3[k][1]=0.0;
			}
			
			for(int k=0;k<N-D->nz;k++)
				D->Der3[k+D->nz][0]=D->Der3[k+D->nz][1]=0.0;

			executeDerShift2(D,D->Der3,D->planfordDer3,D->planbackDer3,N);

			
			for(int k=0;k<D->nz;k++){
				if(D->s.EUz[k*D->nx*D->ny+node]!=40){
					sigma=D->s.EUz[k*D->nx*D->ny+node]*D->sigma/D->nPML;
					D->vz[k*D->nx*D->ny+node]=(D->vz[k*D->nx*D->ny+node]-D->a2*D->Der3[k][0])/(1+sigma);	
					
					
				}
				else
					D->vz[k*D->nx*D->ny+node]=0.0;
			}
			
			
			
		}

	}
}

void tick(domain* D){

	for (int i=0;i<D->nx*D->ny*D->nz;i++){
		D->px[i]=D->pxn[i];
		D->py[i]=D->pyn[i];
		D->pz[i]=D->pzn[i];
		D->pxn[i]=0;
		D->pyn[i]=0;
		D->pzn[i]=0;
	}


}
void SetV(domain* D){
 for (int i=0;i<D->nx*D->ny*D->nz;i++){
		D->vx[i]=0;
		D->vy[i]=0;
		D->vz[i]=0;
	}   
    
}

void freeVectors(domain* D){
    
  free(D->px);free(D->py);free(D->pz);free(D->pxn);free(D->pyn);free(D->pzn);free(D->vx);free(D->vy);free(D->vz);
  free(D->vel);free(D->s.EUx);free(D->s.EUy);free(D->s.EUz);
  fftw_free(D->Der);fftw_free(D->Der2);fftw_free(D->Der3);  
}
