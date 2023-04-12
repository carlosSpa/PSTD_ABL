close all;clear;

fid=fopen('./Energy.dat','rb');
nt=1000;
nB=30;
Dt=0.002;
nParam=20;
ener=zeros(1,nt);
epsilon=zeros(nB,nParam);
ejex=zeros(1,nB);
ejey=zeros(1,nParam);
for m=1:nB
  ejex(1,m)=(4+m-1);
    for k=1:nParam
      ejey(1,k)=0.001+0.0404*(k-1)/19;
        for t=1:nt
         ener(1,t)=fread(fid,1,'double');
    
        end
        ener=ener/max(ener);
        for t=208:nt
           epsilon(m,k)=epsilon(m,k)+ener(1,t)*Dt; 
        end



    end
end
colormap(gray);
colormap(flipud(gray));
contourf(ejey,ejex,log10(epsilon));caxis([-6.5 0]);colorbar;

auxN1=40;
auxK1=20;
auxN2=40;
auxK2=20;
auxN3=40;
auxK3=20;
auxN4=40;
auxK4=20;
auxN5=40;
auxK5=20;
auxN6=40;
auxK6=20;
auxN7=40;
auxK7=20;
for m=1:nB
  for k=1:nParam    
    if ((log10(epsilon(m,k))<-3)&&(log10(epsilon(m,k))>-3.3))
      
        if(auxN1>m)
        
          auxN1=m;
          auxK1=k;
        end   
    end
    
    if ((log10(epsilon(m,k))<-3.5)&&(log10(epsilon(m,k))>-3.7))
      
        if(auxN2>m)
       
          auxN2=m;
          auxK2=k;
        end   
    end
    if ((log10(epsilon(m,k))<-4)&&(log10(epsilon(m,k))>-4.3))
      
        if(auxN3>m)
        
          auxN3=m;
          auxK3=k;
        end   
    end
    if ((log10(epsilon(m,k))<-4.5)&&(log10(epsilon(m,k))>-4.7))
      
        if(auxN4>m)
       
          auxN4=m;
          auxK4=k;
        end   
    end
    if ((log10(epsilon(m,k))<-5)&&(log10(epsilon(m,k))>-5.3))
      
        if(auxN5>m)
       
          auxN5=m;
          auxK5=k;
        end   
    end
    if ((log10(epsilon(m,k))<-5.5)&&(log10(epsilon(m,k))>-5.7))
      
        if(auxN6>m)
       
          auxN6=m;
          auxK6=k;
        end   
    end
    if ((log10(epsilon(m,k))<-6)&&(log10(epsilon(m,k))>-6.3))
      
        if(auxN7>m)
       
          auxN7=m;
          auxK7=k;
        end   
    end
  endfor
endfor
log10(epsilon(auxN1,auxK1))
N1_acc=ejex(1,auxN1)
auxK1_acc=ejey(1,auxK1)
log10(epsilon(auxN2,auxK2))
N2_acc=ejex(1,auxN2)
auxK2_acc=ejey(1,auxK2)
log10(epsilon(auxN3,auxK3))
N3_acc=ejex(1,auxN3)
auxK3_acc=ejey(1,auxK3)
log10(epsilon(auxN4,auxK4))
N4_acc=ejex(1,auxN4)
auxK4_acc=ejey(1,auxK4)
log10(epsilon(auxN5,auxK5))
N5_acc=ejex(1,auxN5)
auxK5_acc=ejey(1,auxK5)
log10(epsilon(auxN6,auxK6))
N6_acc=ejex(1,auxN6)
auxK6_acc=ejey(1,auxK6)
log10(epsilon(auxN7,auxK7))
N7_acc=ejex(1,auxN7)
auxK7_acc=ejey(1,auxK7)