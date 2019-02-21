% plot_bound_wave.m
%
% script to plot the computed bound wave based on 
% readinng a 2d SWAN spec file.
% Also needs user input params.
% Based on Maitane's bound_wave.
%: This Matlab file will compute mean wave parameters (Hs, Mwd          %
%                and Tp) for specific point and time step selected from %
%                a SWAN spectrum output file, and it will compute the   %
%                wave envelope and the associated bound wave (Hasselman,%
%                1962; Herbers et al., 1994; Van Dongeren et al., 2003).%
%                A double summation technique                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%  USER DEFINED PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h0= 18;                                   % mean water depth at the boundary
ndur=1;                                   % duration of each spectra in hours
Dt=1;                                     % time interval for the envelope reconstruction
spec_file='offbreach_2012103100.spc2d';   % name of SWAN 2d spec file

%%%%%%%%%%%%%%%  END of USER defined parameters.  %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%1) LOAD THE 2D SPECTRA
[freq,dir_sort,E2]=read_swan_spc(spec_file);
    
ndir=length(dir_sort);  
nfreq=length(freq);  % freq Hz=[1/sec]
dir_sort=dir_sort;   % dir_sort [degrees]
%    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%2) COMPUTE SPECTRAL CHARACTERISTICS
for ii=1:nfreq-1
    df(ii)=freq(ii+1)-freq(ii);
end
df(nfreq)=df(nfreq-1);
for dd=1:ndir-1
    dtheta(dd)=dir_sort(dd+1)-dir_sort(dd);
    if (dtheta(dd)<0); dtheta(dd)=dtheta(dd)*-1.; end
end
dtheta(ndir)=dtheta(dd-1);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%3) WRITE 1-DIMENSIONAL SPECTRUM
E1=zeros(nfreq,3);
for i=1:nfreq
    E1(i,1)=freq(i);
    E1(i,2)=sum(E2(i,1:ndir).*dtheta(ndir));
    E1(i,3)= (sum(E2(i,1:ndir).*dir_sort(1:ndir,1)'))./sum(E2(i,1:ndir)+1e-10);
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%4) COMPUTE MEAN WAVE PARAMETERS
[Hs,Tp]=specchar(E2,E1,df,dtheta,nfreq,ndir);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%5) COMPUTE FREE SURFACE TIME SERIES
%
[t,eta_tot,env_eta,bwtot]=groupchar(ndur,Dt,E1,df,h0);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%6) Make some plots
%
time_tot=[0:Dt:ndur*3600-Dt]';
figure
hold on
plot(time_tot,eta_tot)
plot(time_tot,env_eta,'r')
plot(time_tot,bwtot,'k')
xlabel('Time (s)')
ylabel('\eta, m')
legend('eta\_tot','env\_eta','bound wave')




