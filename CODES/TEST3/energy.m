close all;clear;

fid=fopen('./EnergyAcc1E.dat','rb');
fid2=fopen('./EnergyAcc2E.dat','rb');
fid3=fopen('./EnergyAcc3E.dat','rb');
fid4=fopen('./EnergyAcc4E.dat','rb');
fid5=fopen('./EnergyAcc5E.dat','rb');
fid6=fopen('./EnergyAcc6E.dat','rb');
fid7=fopen('./EnergyAcc7E.dat','rb');
fid8=fopen('./EnergyRef.dat','rb');

fid11=fopen('./EnergyAcc1D.dat','rb');
fid12=fopen('./EnergyAcc2D.dat','rb');
fid13=fopen('./EnergyAcc3D.dat','rb');
fid14=fopen('./EnergyAcc4D.dat','rb');
fid15=fopen('./EnergyAcc5D.dat','rb');
fid16=fopen('./EnergyAcc6D.dat','rb');
fid17=fopen('./EnergyAcc7D.dat','rb');

fid21=fopen('./EnergyAcc1S.dat','rb');
fid22=fopen('./EnergyAcc2S.dat','rb');
fid23=fopen('./EnergyAcc3S.dat','rb');
fid24=fopen('./EnergyAcc4S.dat','rb');
fid25=fopen('./EnergyAcc5S.dat','rb');
fid26=fopen('./EnergyAcc6S.dat','rb');
fid27=fopen('./EnergyAcc7S.dat','rb');
nt=1999;
ener=zeros(1,nt);
ener2=zeros(1,nt);
ener3=zeros(1,nt);
ener4=zeros(1,nt);
ener5=zeros(1,nt);
ener6=zeros(1,nt);
ener7=zeros(1,nt);
enerRef=zeros(1,nt);

enerD=zeros(1,nt);
ener2D=zeros(1,nt);
ener3D=zeros(1,nt);
ener4D=zeros(1,nt);
ener5D=zeros(1,nt);
ener6D=zeros(1,nt);
ener7D=zeros(1,nt);

enerS=zeros(1,nt);
ener2S=zeros(1,nt);
ener3S=zeros(1,nt);
ener4S=zeros(1,nt);
ener5S=zeros(1,nt);
ener6S=zeros(1,nt);
ener7S=zeros(1,nt);

for t=1:nt
    ener(1,t)=fread(fid,1,'float');
    ener2(1,t)=fread(fid2,1,'float');
    ener3(1,t)=fread(fid3,1,'float');
    ener4(1,t)=fread(fid4,1,'float');
    ener5(1,t)=fread(fid5,1,'float');
    ener6(1,t)=fread(fid6,1,'float');
    ener7(1,t)=fread(fid7,1,'float');
    enerRef(1,t)=fread(fid8,1,'float');
    
    enerD(1,t)=fread(fid11,1,'float');
    ener2D(1,t)=fread(fid12,1,'float');
    ener3D(1,t)=fread(fid13,1,'float');
    ener4D(1,t)=fread(fid14,1,'float');
    ener5D(1,t)=fread(fid15,1,'float');
    ener6D(1,t)=fread(fid16,1,'float');
    ener7D(1,t)=fread(fid17,1,'float');
    
    enerS(1,t)=fread(fid21,1,'float');
    ener2S(1,t)=fread(fid22,1,'float');
    ener3S(1,t)=fread(fid23,1,'float');
    ener4S(1,t)=fread(fid24,1,'float');
    ener5S(1,t)=fread(fid25,1,'float');
    ener6S(1,t)=fread(fid26,1,'float');
    ener7S(1,t)=fread(fid27,1,'float');
    
end
ener=ener/max(ener);
ener2=ener2/max(ener2);
ener3=ener3/max(ener3);
ener4=ener4/max(ener4);
ener5=ener5/max(ener5);
ener6=ener6/max(ener6);
ener7=ener7/max(ener7);
enerRef=enerRef/max(enerRef);

enerD=enerD/max(enerD);
ener2D=ener2D/max(ener2D);
ener3D=ener3D/max(ener3D);
ener4D=ener4D/max(ener4D);
ener5D=ener5D/max(ener5D);
ener6D=ener6D/max(ener6D);
ener7D=ener7D/max(ener7D);

enerS=enerS/max(enerS);
ener2S=ener2S/max(ener2S);
ener3S=ener3S/max(ener3S);
ener4S=ener4S/max(ener4S);
ener5S=ener5S/max(ener5S);
ener6S=ener6S/max(ener6S);
ener7S=ener7S/max(ener7S);

%plot(ener3);hold on;
%plot(ener4/max(ener4));hold on;
%plot(ener);
%hold on;
%plot(enerRef);
%hold on;
%plot(enerRef);
%hold on;
%plot(ener4S,'k');

DeltaE=zeros(7,1);
DeltaD=zeros(7,1);
DeltaS=zeros(7,1);

for t=1:nt
   DeltaE(1,1)=DeltaE(1,1)+abs(enerRef(1,t)-ener(1,t));% *(enerRef(1,t)-ener(1,t)); 
    DeltaE(2,1)=DeltaE(2,1)+abs(enerRef(1,t)-ener2(1,t));% *(enerRef(1,t)-ener2(1,t)); 
    DeltaE(3,1)=DeltaE(3,1)+abs(enerRef(1,t)-ener3(1,t));% *(enerRef(1,t)-ener3(1,t)); 
    DeltaE(4,1)=DeltaE(4,1)+abs(enerRef(1,t)-ener4(1,t));% *(enerRef(1,t)-ener4(1,t)); 
    DeltaE(5,1)=DeltaE(5,1)+abs(enerRef(1,t)-ener5(1,t));% *(enerRef(1,t)-ener5(1,t)); 
    DeltaE(6,1)=DeltaE(6,1)+abs(enerRef(1,t)-ener6(1,t));% *(enerRef(1,t)-ener6(1,t)); 
    DeltaE(7,1)=DeltaE(7,1)+abs(enerRef(1,t)-ener7(1,t));% *(enerRef(1,t)-ener7(1,t)); 
    
    DeltaD(1,1)=DeltaD(1,1)+abs(enerRef(1,t)-enerD(1,t));
    DeltaD(2,1)=DeltaD(2,1)+abs(enerRef(1,t)-ener2D(1,t));% *(enerRef(1,t)-ener2(1,t)); 
    DeltaD(3,1)=DeltaD(3,1)+abs(enerRef(1,t)-ener3D(1,t));% *(enerRef(1,t)-ener3(1,t)); 
    DeltaD(4,1)=DeltaD(4,1)+abs(enerRef(1,t)-ener4D(1,t));% *(enerRef(1,t)-ener4(1,t)); 
    DeltaD(5,1)=DeltaD(5,1)+abs(enerRef(1,t)-ener5D(1,t));% *(enerRef(1,t)-ener5(1,t)); 
    DeltaD(6,1)=DeltaD(6,1)+abs(enerRef(1,t)-ener6D(1,t));
    DeltaD(7,1)=DeltaD(7,1)+abs(enerRef(1,t)-ener7D(1,t));
    
    DeltaS(1,1)=DeltaS(1,1)+abs(enerRef(1,t)-enerS(1,t));
    DeltaS(2,1)=DeltaS(2,1)+abs(enerRef(1,t)-ener2S(1,t));% *(enerRef(1,t)-ener2(1,t)); 
    DeltaS(3,1)=DeltaS(3,1)+abs(enerRef(1,t)-ener3S(1,t));% *(enerRef(1,t)-ener3(1,t)); 
    DeltaS(4,1)=DeltaS(4,1)+abs(enerRef(1,t)-ener4S(1,t));% *(enerRef(1,t)-ener4(1,t)); 
    DeltaS(5,1)=DeltaS(5,1)+abs(enerRef(1,t)-ener5S(1,t));
    DeltaS(6,1)=DeltaS(6,1)+abs(enerRef(1,t)-ener6S(1,t));
    DeltaS(7,1)=DeltaS(7,1)+abs(enerRef(1,t)-ener7S(1,t));
end

for i=1:7
   DeltaE(i,1)=(DeltaE(i,1))*0.002; 
   
   DeltaD(i,1)=(DeltaD(i,1))*0.002; 
   
   
   DeltaS(i,1)=(DeltaS(i,1))*0.002; 
     
end