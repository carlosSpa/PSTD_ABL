% ricker_wavelet(dt,nt,t0,nu1);
%  This function generates a Ricker wavelet with the following parameters:
%    dt:  time discretization
%    nt:  number of time samples
%    t0:  time offset
%    f0:  peak frequency
%    

function [A]=ricker_wavelet_tis0(dt,nt,t0,f0);

t=[1:nt]*dt-dt;
B=pi*pi*f0*f0;
A=(1-2*B*(t-t0).*(t-t0)).*exp(-B*(t-t0).*(t-t0));
A(1)=0;     %Imposed zero value at t=0.
end