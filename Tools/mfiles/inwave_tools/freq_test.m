Tp=13.;
fp=1/Tp;
fnyq=3.0*fp
fmax=fnyq
fmin=fp/20.0
%fmin=0.;
nfreq=1000
df=(fmax-fmin)/(nfreq-1)
%df=round(df,8)
%df=2.5e-4
fmin=df;
for i=1:nfreq
  f(i)=(i-1)*df+fmin;
end
dur=1./df
dt=.1
Insteps=round(dur/dt,0)
%Insteps=floor(dur/dt)
CompFn=zeros(1,Insteps*1);
for i=1:Insteps*1
  for j=1:nfreq
    amp(j)=j*.0001;
    phase(j)=rand(1);
    cff=(i-1)*dt;
    CompFn(i)=CompFn(i)+...
      amp(j)*cos(2*pi*f(j)*cff+phase(j)*2*pi);
  end
end

figure
hold on
plot([1:Insteps*1],CompFn,'r-+')
%plot([Insteps*3+1+10:Insteps*20+Insteps*3+10],CompFn,'k')
%A=0;
%plot([Insteps+1+A:Insteps*3+Insteps+A],CompFn,'g-+')
xhil=hilbert(CompFn);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% start here
[freq,dir_sort,E2]=read_swan_spc('offbreach_2012103100.spc2d');
E2=E2/1025;
figure
pcolorjw(dir_sort,freq,E2)

    ndir=length(dir_sort);  
    nfreq=length(freq);  % freq Hz=[1/sec]
    dir_sort=dir_sort;   % dir_sort [degrees]
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%% COMPUTE SPECTRAL CHARACTERISTICS  %%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    clear df
    for ii=1:nfreq-1
        df(ii)=freq(ii+1)-freq(ii);
    end
    df(nfreq)=df(nfreq-1);
    for dd=1:ndir-1
        dtheta(dd)=dir_sort(dd+1)-dir_sort(dd);
    end
    dtheta(ndir)=dtheta(dd-1);
    dtheta=-dtheta;
       
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%   WRITE 1-DIMENSIONAL SPECTRUM   %%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    E1=zeros(nfreq,3);
    
    for i=1:nfreq
        E1(i,1)=freq(i);
        E1(i,2)=sum(E2(i,1:ndir).*dtheta(ndir));
        E1(i,3)= (sum(E2(i,1:ndir).*dir_sort(1:ndir,1)'))./sum(E2(i,1:ndir));
    end


    [Hs,Tp]=specchar(E2,E1,df,dtheta,nfreq,ndir);
  %seems big hsig

Dt=0.1
ndur=1
time=[0:Dt:(ndur*3600+5*60)-Dt]';  %sec  % extra 10 minutes for time transitions
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

figure
plot(eta_tot,'b')
hold on
plot(env_eta,'g')


% do a filter on ampzeta
      zeta_filt=env_eta;
      env_eta_filt=env_eta;
       
      dt=Dt
        cff=0.
        numavg=floor(5./dt)
        numavg=numavg-1*(1-mod(numavg,2)) %!force odd
        for i=1:numavg
          cff=cff+zeta_filt(i);
        end
        numavgh=(numavg-1)/2
        for i=numavgh+1:length(zeta_filt)-numavgh
          env_eta_filt(i)=cff/numavg;
          cff=cff-zeta_filt(i-numavgh)+zeta_filt(i+numavgh);
        end
plot(env_eta_filt,'r')





env_eta_filt2=env_eta_filt;
for mm=150:length(zeta_filt)-150
  env_eta_filt2(mm)=max(env_eta_filt(mm-149:mm+149));
end
      env_eta_filt3=env_eta_filt2;
        for i=1:numavg
          cff=cff+env_eta_filt2(i);
        end
        for i=numavgh+1:length(zeta_filt)-numavgh
          env_eta_filt3(i)=cff/numavg;
          cff=cff-env_eta_filt2(i-numavgh)+env_eta_filt2(i+numavgh);
        end

























