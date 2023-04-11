close all;clear;
fid=fopen('./R2dEuler.dat','rb');
%fid2=fopen('./R2dDamped.dat','rb');
%fid3=fopen('./R2dSponge.dat','rb');
Nx=fread(fid,1,'int');
Ny=fread(fid,1,'int');
nT=fread(fid,1,'int');
results=zeros(Ny,Nx);
%NxD=fread(fid2,1,'int');
%NyD=fread(fid2,1,'int');
%nTd=fread(fid2,1,'int');
%resultsD=zeros(Nx,Ny);

%NxS=fread(fid3,1,'int');
%NyS=fread(fid3,1,'int');
%nTs=fread(fid3,1,'int');
%resultsS=zeros(Nx,Ny);

source=zeros(1,nT);
%sourceD=zeros(1,nT);
%sourceS=zeros(1,nT);
tiempo=zeros(1,nT);
for i=1:nT
  tiempo(1,i)=(i-1)*0.002;
   for l=1:Nx
     for j=1:Ny
     results(Ny-j+1,l)=fread(fid,1,'float');
%     resultsD(l,j)=fread(fid2,1,'float');
%     resultsS(l,j)=fread(fid3,1,'float');
   end
end 

  source(1,i)=results(ceil(Nx/2),ceil(Ny/2));
 % sourceD(1,i)=resultsD(ceil(Nx/2),ceil(Ny/2));
 % sourceS(1,i)=resultsS(ceil(Nx/2),ceil(Ny/2));
  if(mod(i,5)==0)
i
  
   imagesc(results);
   colormap gray;cmap = colormap;cmap = flipud(cmap);colormap(cmap),set(gca,'ydir','normal');
   colorbar;
   %caxis([-1000 1500])
   %plot(results(ceil(Nx/2),:));
   %hold on;
   %plot(resultsD(ceil(Nx/2),:));ylim([-60 120])
   
  pause(0.0002);clc;
   end
end

figure 
plot(tiempo,source/(abs(max(source))));
%hold on;
%plot(tiempo,sourceD/(abs(max(sourceD))));
%hold on;
%plot(tiempo,sourceS/(abs(max(sourceS))));