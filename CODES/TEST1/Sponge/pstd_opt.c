#include "pstd_opt.h"



void ricker_wavelet_tis0(domain* D,double dt,int nt,double t0,double f0,int l){

	D->source=(double*)malloc((D->nt)*sizeof(double));
	double B=0;double Pi,aux;  Pi=3.141592654;
	B=Pi*Pi*f0*f0;
	if(l==0){
		for(int t=0;t<D->nt;t++){
			aux=0.0;
			if(t!=0){
				for(int j=0;j<=t;j++)
					aux=100*((1-2*B*(j*D->Dt-t0)*(j*D->Dt-t0))*exp(-B*(j*D->Dt-t0)*(j*D->Dt-t0)))+aux; 
				D->source[t]=aux*dt;
			}
			else
				D->source[t]=0.0;
		}
	}
	else{
		for(int t=D->nt-1;t>=0;t--){
			aux=0.0;
		
			for(int j=D->nt-1;j>=t;j--)
				aux=100*((1-2*B*(j*D->Dt-t0-1)*(j*D->Dt-t0-1))*exp(-B*(j*D->Dt-t0-1)*(j*D->Dt-t0-1))
					+(1-2*B*(j*D->Dt-t0-2)*(j*D->Dt-t0-2))*exp(-B*(j*D->Dt-t0-2)*(j*D->Dt-t0-2))
					+(1-2*B*(j*D->Dt-t0-3)*(j*D->Dt-t0-3))*exp(-B*(j*D->Dt-t0-3)*(j*D->Dt-t0-3)))+aux;
			D->source[t]=aux*dt;
			
			
		}
	}

}
void initSource(domain* D, double frec,int nf,double n0){

	D->source=(double*)malloc((D->nt)*sizeof(double));

	for(int t=0;t<D->nt;t++){
		if(t<nf)
			D->source[t]=100*sin(2*3.1415926*frec*(t-n0)*D->Dt)/(2*3.1415926*frec*(t-n0)*D->Dt);
		else
			D->source[t]=0.0;
	}

}
void source(domain* D,int type,int nS,int t){

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
	else{
		for(int i=-1;i<2;i++){
			for(int j=-1;j<2;j++){
				for(int k=-1;k<2;k++){
					D->pxn[nS+i+j*D->nx+k*D->nx*D->ny]=D->source[D->nt-1-t]*exp(-1*((i)*(i)+(j)*(j)+(k)*(k)))+D->pxn[nS+i+j*D->nx+k*D->nx*D->ny];
					D->pyn[nS+i+j*D->nx+k*D->nx*D->ny]=D->source[D->nt-1-t]*exp(-1*((i)*(i)+(j)*(j)+(k)*(k)))+D->pyn[nS+i+j*D->nx+k*D->nx*D->ny];
					D->pzn[nS+i+j*D->nx+k*D->nx*D->ny]=D->source[D->nt-1-t]*exp(-1*((i)*(i)+(j)*(j)+(k)*(k)))+D->pzn[nS+i+j*D->nx+k*D->nx*D->ny];
				}
			}
		}
		
	}
}
void initParam(domain* D,double Lx,double Ly,double Lz,double Dx,double Dt, double S, double tiempo,int nPML,double sigma,double cerjan,int Nb,int Nparam){
	D->c=2000;double rho=1;//1.21;
	D->a1=D->c*D->c*Dt;
	D->a2=Dt;//Dt/(rho);
	D->S=S;
	D->Dx=Dx;//D->c*Dt/S;
	D->nx=(int)(Lx/D->Dx)+nPML;
	D->ny=(int)(Ly/D->Dx)+nPML;
	D->nz=(int)(Lz/D->Dx)+nPML;
	D->sigma=sigma;
	D->nPML=nPML/2;
	D->nt=(int)(tiempo/Dt);
	D->Dt=Dt;	
	D->cerjan=cerjan;
	D->Nb=Nb;
	D->Nparam=Nparam;
}

void initArrays(domain* D){
	int nodes;
	nodes=D->nx*D->ny*D->nz;
	D->dpxn=(double*)malloc((nodes)*sizeof(double));
	D->dpyn=(double*)malloc((nodes)*sizeof(double));
	D->dpzn=(double*)malloc((nodes)*sizeof(double));	
	D->dpxa=(double*)malloc((nodes)*sizeof(double));
	D->dpya=(double*)malloc((nodes)*sizeof(double));
	D->dpza=(double*)malloc((nodes)*sizeof(double));
	D->pxa=(double*)malloc((nodes)*sizeof(double));
	D->pya=(double*)malloc((nodes)*sizeof(double));
	D->pza=(double*)malloc((nodes)*sizeof(double));
	D->pxn=(double*)malloc((nodes)*sizeof(double));
	D->pyn=(double*)malloc((nodes)*sizeof(double));
	D->pzn=(double*)malloc((nodes)*sizeof(double));
	D->vel=(double*)malloc((nodes)*sizeof(double));
	D->s.EUx=(int*)malloc((nodes)*sizeof(int));
	D->s.EUy=(int*)malloc((nodes)*sizeof(int));
	D->s.EUz=(int*)malloc((nodes)*sizeof(int));

	
	

	for(int i=0;i<nodes;i++){
		D->dpxn[i]=D->dpyn[i]=D->dpzn[i]=D->dpxa[i]=D->dpya[i]=D->dpza[i]=D->pxa[i]=D->pya[i]=D->pza[i]=D->pxn[i]=D->pyn[i]=D->pzn[i]=0.0;
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
				
				/*if((i<D->nPML)&&(j<D->nPML)){//&&(k<D->nPML))
					D->s.EUx[node]=40;D->s.EUy[node]=40;D->s.EUz[node]=40;
				}
				else if((i>D->nx-1-D->nPML)&&(j<D->nPML)){//&&(k<D->nPML))
					D->s.EUx[node]=40;D->s.EUy[node]=40;D->s.EUz[node]=40;}
				else if((i<D->nPML)&&(j>D->ny-1-D->nPML)){//&&(k<D->nPML))
					D->s.EUx[node]=40;D->s.EUy[node]=40;D->s.EUz[node]=40;}
				else if((i>D->nx-1-D->nPML)&&(j>D->ny-1-D->nPML)){//&&(k>D->nz-1+D->nPML))
					D->s.EUx[node]=40;D->s.EUy[node]=40;D->s.EUz[node]=40;}
				else if((k<D->nPML)&&(j<D->nPML)){//&&(k<D->nPML))
					D->s.EUx[node]=40;D->s.EUy[node]=40;D->s.EUz[node]=40;}
				else if((k>D->nz-1-D->nPML)&&(j<D->nPML)){//&&(k<D->nPML))
					D->s.EUx[node]=40;D->s.EUy[node]=40;D->s.EUz[node]=40;}
				else if((k<D->nPML)&&(j>D->ny-1-D->nPML)){//&&(k<D->nPML))
					D->s.EUx[node]=40;D->s.EUy[node]=40;D->s.EUz[node]=40;}
				else if((k>D->nz-1-D->nPML)&&(j>D->ny-1-D->nPML)){//&&(k>D->nz-1+D->nPML))
					D->s.EUx[node]=40;D->s.EUy[node]=40;D->s.EUz[node]=40;}
				else if((i<D->nPML)&&(k<D->nPML)){//&&(k<D->nPML))
					D->s.EUx[node]=40;D->s.EUy[node]=40;D->s.EUz[node]=40;}
				else if((i>D->nx-1-D->nPML)&&(k<D->nPML)){//&&(k<D->nPML))
					D->s.EUx[node]=40;D->s.EUy[node]=40;D->s.EUz[node]=40;}
				else if((i<D->nPML)&&(k>D->nz-1-D->nPML)){//&&(k<D->nPML))
					D->s.EUx[node]=40;D->s.EUy[node]=40;D->s.EUz[node]=40;}
				else if((i>D->nx-1-D->nPML)&&(k>D->nz-1-D->nPML)){//&&(k>D->nz-1+D->nPML))
					D->s.EUx[node]=40;D->s.EUy[node]=40;D->s.EUz[node]=40;}
				
				else{
					
				}*/
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
	D->CteDer=4*pi*pi/(N*N*D->Dx*D->Dx*N);
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

void executeDer(domain* D,fftw_complex* p,fftw_plan planfordDer,fftw_plan planbackDer,int N){

	fftw_execute(planfordDer);
	for(int i=0;i<N;i++){
		double aux1,aux2;
		if(i<N/2-1){
			aux1=-D->CteDer*i*i*p[i][0];
			aux2=-D->CteDer*i*i*p[i][1];
		}
		else if(i==N/2-1){
			aux1=aux2=0.0;
		}
		else{
			aux1=-D->CteDer*(i-N)*(i-N)*p[i][0];
			aux2=-D->CteDer*(i-N)*(i-N)*p[i][1];
		}
		p[i][0]=aux1;
		p[i][1]=aux2;
	}

	fftw_execute(planbackDer);
}

void updatePx(domain* D){
	
	
	int N;double c;double dist;
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
				D->Der[i][0]=D->pxa[node+i]+D->pya[node+i]+D->pza[i+node];
				D->Der[i][1]=0.0;
			}
			
			for(int i=0;i<N-D->nx;i++)
				D->Der[i+D->nx][0]=D->Der[i+D->nx][1]=0.0;

			
			executeDer(D,D->Der,D->planfordDer,D->planbackDer,N);

			
			for(int i=0;i<D->nx;i++){
				if(D->s.EUx[i+node]!=40){
					dist=sqrt(D->s.EUx[i+node]*D->s.EUx[i+node]+D->s.EUy[i+node]*D->s.EUy[i+node]+D->s.EUz[i+node]*D->s.EUz[i+node]);
					D->a1=D->vel[i+node]*D->vel[i+node]*D->Dt;
					c=exp(-D->cerjan*D->cerjan*dist*dist);
					D->dpxn[i+node]=(D->dpxa[i+node]+D->a1*D->Der[i][0]);
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


void updatePy(domain* D){
	
	
	int N;double c;double dist;
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
				D->Der2[j][0]=D->pxa[j*D->nx+node]+D->pya[j*D->nx+node]+D->pza[j*D->nx+node];
				D->Der2[j][1]=0.0;
			}
			
			for(int j=0;j<N-D->ny;j++)
				D->Der2[j+D->ny][0]=D->Der2[j+D->ny][1]=0.0;

			
			executeDer(D,D->Der2,D->planfordDer2,D->planbackDer2,N);

			
			for(int j=0;j<D->ny;j++){
				if(D->s.EUy[j*D->nx+node]!=40){
					dist=sqrt(D->s.EUx[j*D->nx+node]*D->s.EUx[j*D->nx+node]+D->s.EUy[j*D->nx+node]*D->s.EUy[j*D->nx+node]+D->s.EUz[j*D->nx+node]*D->s.EUz[j*D->nx+node]);
					D->a1=D->vel[j*D->nx+node]*D->vel[j*D->nx+node]*D->Dt;
					c=exp(-D->cerjan*D->cerjan*dist*dist);
					D->dpyn[j*D->nx+node]=(D->dpya[j*D->nx+node]+D->a1*D->Der2[j][0]);
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


void updatePz(domain* D){
	
	
	int N;double c;double dist;
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
				D->Der3[k][0]=D->pxa[k*D->nx*D->ny+node]+D->pya[k*D->nx*D->ny+node]+D->pza[k*D->nx*D->ny+node];
				D->Der3[k][1]=0.0;
			}
			
			for(int k=0;k<N-D->nz;k++)
				D->Der3[k+D->nz][0]=D->Der3[k+D->nz][1]=0.0;

			executeDer(D,D->Der3,D->planfordDer3,D->planbackDer3,N);

			
			for(int k=0;k<D->nz;k++){
				if(D->s.EUz[k*D->nx*D->ny+node]!=40){
					dist=sqrt(D->s.EUx[k*D->nx*D->ny+node]*D->s.EUx[k*D->nx*D->ny+node]
					+D->s.EUy[k*D->nx*D->ny+node]*D->s.EUy[k*D->nx*D->ny+node]+D->s.EUz[k*D->nx*D->ny+node]*D->s.EUz[k*D->nx*D->ny+node]);
					D->a1=D->vel[k*D->nx*D->ny+node]*D->vel[k*D->nx*D->ny+node]*D->Dt;
					c=exp(-D->cerjan*D->cerjan*dist*dist);
					D->dpzn[k*D->nx*D->ny+node]=(D->dpza[k*D->nx*D->ny+node]+D->a1*D->Der3[k][0]);
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



void tick(domain* D){

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
void freeVectors(domain* D){
    free(D->dpxa);free(D->dpya);free(D->dpza);free(D->dpxn);free(D->dpyn);free(D->dpzn);free(D->pxn);free(D->pyn);free(D->pzn);free(D->pxa);free(D->pya);free(D->pza);
    free(D->vel);free(D->s.EUx);free(D->s.EUy);free(D->s.EUz);
    fftw_free(D->Der);fftw_free(D->Der2);fftw_free(D->Der3);
}
