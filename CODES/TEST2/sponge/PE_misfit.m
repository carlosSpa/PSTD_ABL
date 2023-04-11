% [PM,EM,TPM,TEM]=PE_misfit(nt,SREF,S,plotme);
%  This function computes the phase and envelope misfit between two time
%  signals, as per Kristekova et al (2006)
%    nt:   number of time samples
%    SREF: reference time signal
%    S:    compared time signal
%    plotme: (optional) set to 'yes' if you want an image of the misfits (default:
%    'no')

function [PM,EM,TPM,TEM]=PE_misfit(nt,SREF,S,varargin);

%%% Default vales
plotme=0;
    
%%% Check argument list length
if nargin>=5
    error('Too many input arguments') 
end

%%% Set optional values
if nargin>=4
plotme = varargin{1};
end

if strcmp(char(varargin{1}),'yes')==1
   plotme=1;
elseif strcmp(char(varargin{1}),'no')==1
   plotme=0;
else
   error('plotme can only have values "yes" or "no"')
end

%%% Compute misfits
HREF=hilbert(SREF);
H=hilbert(S);

TEM=abs(HREF)-abs(H);
TPM=abs(HREF).*angle(HREF./H)/pi;

TEM=TEM/max(abs(HREF));
TPM=TPM/max(abs(HREF));

NP=0;D=0;NE=0;
for i=1:nt
    NP=NP+(abs(HREF(i))*angle(HREF(i)/H(i)))^2;
    D=D+abs(HREF(i))^2;
    NE=NE+(abs(HREF(i))-abs(H(i)))^2;
end

PM=sqrt(NP)/pi/sqrt(D);
EM=sqrt(NE)/sqrt(D);

%%% Plot misfits
if plotme==1
x1=[1:nt];
y1=SREF;
y1b=S;
y2=TPM;
y2b=TEM;

figure
%title('Windowed phase and envelope misfits')
subplot(211)
plot(x1,y1,'k'),hold on,plot(x1,y1b,'b'),legend('Wave 1','Wave 2'),set(gca,'XTickLabel',[]),grid on%,set(gca,'YGrid','off')
subplot(212)
plot(x1,y2,'r'),hold on,plot(x1,y2b,'g'),grid on,legend('PM','EM')
end

end
