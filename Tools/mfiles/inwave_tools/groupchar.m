function [t,eta_tot,env_eta,bwtot]=groupchar(ndur,Dt,E1,df1,h0)

% h0= 21;   % mean water depth at the boundary
% ndur=1;   % duration of each spectra in hours
% Dt=1;     % time interval for the envelope reconstruction (sec)
% E1=       % 1D wave spectra
% df1=      % frequency increment

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% CREATE TIME VECTOR %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

time=[0:Dt:(ndur*3600)-Dt]';  %sec
np=length(time);

freq=E1(:,1);
Nt = length(freq);
T = Nt/freq(end)/2;
dT = 1/(2*freq(end));
if mod(np,2),np=np+1;end % make sure it is even
df = 1/(np*dT);

Se=E1(:,2);
fs=freq;
fcutoff=1/30;
fmin=1/400;

% interpolate for freq.  [1:(N/2)-1]*df and create 2-sided, uncentered spectra
% ----------------------------------------------------------------------------
f = [1:(np/2)-1]'*df;

Si = interp1q(fs,Se,f);
%Si = interp1(fs,Se,f);
Si(isnan(Si)==1)=0;
dir_int=interp1(E1(:,1),E1(:,3),f);
dir_int(isnan(dir_int)==1)=0;

Su=[0; Si; 0; Si((np/2)-1:-1:1)];

Zi = [zeros(1,1); randn((np/2)-1,1).*2*pi; zeros(1,1)];
a1=((2*Si).*df).^0.5;
ttm=(a1/2).*exp(Zi(2:end-1).*(-1*sqrt(-1)));

% Compute the Fourier Coefficients
% ----------------------------------------------------------------------------

A                = zeros(np,1);
A(2:(np/2),:)    = ttm;
A((np/2+2):np,:) = conj(A(np/2:-1:2,:));

% Compute the free surface elevation time series (eta_tot)
% ----------------------------------------------------------------------------

eta_tot= real(ifft(A))/df;
t(:,1)     = linspace(0,(np-1)*dT,np)';

clear ttm A

%%%%%%%%%%%%%%%%%%%%%%%%% COMPUTE THE WAVE ENVELOPE   %%%%%%%%%%%%%%%%%%%%%%%%

env_eta=abs(hilbert(eta_tot));

%%%%%%%%%%%%%%%%%%%%%%%%%% COMPUTE THE BOUND WAVE   %%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize arrays
% ---------------------------------------------------------

%k=wavenum(f,h0);
g=9.81;
wa=2*pi*f;
for mm=1:length(f)
% k(mm)=wavenumber(g,wa(mm),0,h0,0,0);
  k(mm)=waven(1/f(mm),h0);
end
nf=length(f);
A3=zeros(nf,nf);
Z_bw=zeros(nf,nf);
DTOT=zeros(nf,nf);
bwtot=zeros(length(time),1);

% Compute the energy transfer for each pair of frequency bins
% -----------------------------------------------------------

dum=find(f<fcutoff);

  for f1=dum(end):nf-1
    for f2=f1+1:nf
        
        DDf=f(f2)-f(f1);%frequency difference     
        DDtheta=(dir_int(f2)-dir_int(f1)).*pi/180; %direction difference    
        k3=sqrt(k(f2).^2+k(f1).^2-2*k(f2)*k(f1)*cos(DDtheta)); %wave number      
        Z_bw(f1,f2)=Zi(f1)-Zi(f2)+pi;  % phase of the bound wave
        
        D1=9.81*k(f1)*k(f2)*cos(DDtheta+pi)/(8*pi^2*f(f1)*f(f2))*cosh(k3*h0)/(cosh(k(f1)*h0)*cosh(k(f2)*h0));
        D2=-9.81*(DDf)/((9.81*k3*tanh(k3*h0)-(2*pi).^2*DDf^2)*f(f1)*f(f2));
        D3=DDf*(((2*pi)^4*(f(f1)*f(f2))^2)/(9.81^2)-k(f1)*k(f2)*cos(DDtheta+pi));
        D4=-0.5*(-f(f1)*k(f2)^2/(cosh(k(f2)*h0)).^2+f(f2)*k(f1)^2/(cosh(k(f1)*h0)).^2);
        DTOT(f1,f2)=D1+D2*(D3+D4);        

        if DDf>fcutoff
            DTOT(f1,f2)=0;  % do not consider frequencies higher than fcutoff
        end   
        
        if DDf<fmin
            DTOT(f1,f2)=0;  % do not consider frequencies lower than fmin
        end   
        if (isnan(DTOT(f1,f2))==1) 
            DTOT(f1,f2)=0;
        end
        
        E3(f1,f2)=2.*DTOT(f1,f2).^2*Si(f1)*Si(f2)*df;  %energy of the bound wave
        A3(f1,f2)=sqrt(2*E3(f1,f2)*df);
        
        if DDf<=fcutoff
            bwtot(:,1)=bwtot(:,1)+A3(f1,f2)./2*cos(2*pi*DDf.*t+Z_bw(f1,f2));
        end      
       
    end
  end
end