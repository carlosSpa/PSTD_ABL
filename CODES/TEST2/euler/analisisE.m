close all;
clear;


Nt=1999;
N=100;
ResultsB=zeros(N,Nt);
ResultsF=zeros(N,Nt);
Image=zeros(N,7);
PM=zeros(1,7);
EM=zeros(1,7);
for l=1:7
    if l==1
    fid1=fopen('./PBAcc1E.dat','rb');
    fid2=fopen('./PFAcc1E.dat','rb');
    end
    if l==2
    fid1=fopen('./PBAcc2E.dat','rb');
    fid2=fopen('./PFAcc2E.dat','rb');
    end
    if l==3
    fid1=fopen('./PBAcc3E.dat','rb');
    fid2=fopen('./PFAcc3E.dat','rb');
    end
    if l==4
    fid1=fopen('./PBAcc4E.dat','rb');
    fid2=fopen('./PFAcc4E.dat','rb');
    end
    if l==5
    fid1=fopen('./PBAcc5E.dat','rb');
    fid2=fopen('./PFAcc5E.dat','rb');
    end
    if l==6
    fid1=fopen('./PBAcc6E.dat','rb');
    fid2=fopen('./PFAcc6E.dat','rb');
    end
    if l==7
    fid1=fopen('./PBAcc7E.dat','rb');
    fid2=fopen('./PFAcc7E.dat','rb');
    end
for n=1:Nt
   for i=1:N
      ResultsB(i,Nt+1-n)=fread(fid1,1,'float'); 
      ResultsF(i,n)=fread(fid2,1,'float'); 
       
   end  
end


for i=1:N
for n=1:Nt
    
   Image(i,l)=Image(i,l)+(ResultsF(i,n)*ResultsB(i,n))/(max(ResultsF(i,:))*max(ResultsF(i,:))); 
end
end
Image2=zeros(1,N);
for i=5:N
   Image2(1,i-4)=Image(i,l); 
end
figure
plot(Image2/max(abs(Image2)));
dx=40;
nx=100;
c=2000;
f0=sqrt(3)*1/(c/10);
x0=1000;
sol=ricker_wavelet_tis0(dx,nx,x0,f0)+ricker_wavelet_tis0(dx,nx,x0+1000,f0)+ricker_wavelet_tis0(dx,nx,x0+2000,f0);
hold on;
plot(sol);

[PM(1,l), EM(1,l),a,b]=PE_misfit(N,sol,Image2/(max(abs(Image2))),'no');

end
