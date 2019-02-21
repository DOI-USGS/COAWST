%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Maitane Olabarrieta                                           %
% Date:   29/07/2014                                                    %
%: This Matlab file will compute mean wave parameters (Hs, Mwd          %
%                and Tp) for specific point and time step selected from %
%                a SWAN spectrum output file, and it will compute the   %
%                wave envelope and the associated bound wave (Hasselman,%
%                1962; Herbers et al., 1994; Van Dongeren et al., 2003).%
%                A double summation technique                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%  USER DEFINED PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h0= 21;  % mean water depth at the boundary
ndur=1;   % duration of each spectra in hours
Dt=1;     % time interval for the envelope reconstruction

%%% LOAD TIDAL SIGNAL

temp_tide=load('WLevel_Cascais_201401_1h.dat');
YY=temp_tide(:,1);
MM=temp_tide(:,2);
DD=temp_tide(:,3);
HH=temp_tide(:,4);
MINU=temp_tide(:,5);
SEC=temp_tide(:,6);

tide_time=datenum(YY,MM,DD,HH,MINU,SEC)
tide_zeta=temp_tide(:,7);

clear temp_tide

plot((tide_time-tide_time(1,1)).*24,tide_zeta)

%%% SPECIFY THE DIFFERENT SPECTRAL FILES, NUMBER OF FILES MUST BE EQUAL TO
%%% nfiles

pathfiles='H:\My_Projects\FAR_INFRAGRAVITY\CASCAIS_SPECTRA\Waves_Spectra_SWAN_tides\'
count=0;
for ndays=5:7
    for nhours=1:24
        count=count+1;
        if nhours<11
            filename(:,count)=['Cascais_Spec2D_2014010',num2str(ndays),'0',num2str(nhours-1),'0000.spc'];
        else
            filename(:,count)=['Cascais_Spec2D_2014010',num2str(ndays),num2str(nhours-1),'0000.spc'];
        end
    end
end

nfiles=count; % number of spectral files to consider

for file_num=1:nfiles
    
    file_num
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%% LOAD THE 2D SPECTRA %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [freq,dir_sort,E2]=read_swan_spc([pathfiles,filename(:,file_num)']);
    
    ndir=length(dir_sort);  
    nfreq=length(freq);  % freq Hz=[1/sec]
    dir_sort=dir_sort;   % dir_sort [degrees]
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%% COMPUTE SPECTRAL CHARACTERISTICS  %%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    for ii=1:nfreq-1
        df(ii)=freq(ii+1)-freq(ii);
    end
    df(nfreq)=df(nfreq-1);
    for dd=1:ndir-1
        dtheta(dd)=dir_sort(dd+1)-dir_sort(dd);
    end
    dtheta(ndir)=dtheta(dd-1);
       
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%   WRITE 1-DIMENSIONAL SPECTRUM   %%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    E1=zeros(nfreq,3);
    
    for i=1:nfreq
        E1(i,1)=freq(i);
        E1(i,2)=sum(E2(i,1:ndir).*dtheta(ndir));
        E1(i,3)= (sum(E2(i,1:ndir).*dir_sort(1:ndir,1)'))./sum(E2(i,1:ndir));
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%    COMPUTE MEAN WAVE PARAMETERS   %%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
    [Hs,Tp]=specchar(E2,E1,df,dtheta,nfreq,ndir);
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%% COMPUTE FREE SURFACE TIME SERIES %%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [t,eta_tot,env_eta,bwtot]=groupchar(ndur,Dt,E1,df,h0+tide_zeta(file_num));
               
    eta_dur(:,file_num)=eta_tot;
    env_dur(:,file_num)=env_eta;
    bwtot_dur(:,file_num)=bwtot;
    
    clearvars -except file_num eta_dur env_dur bwtot_dur filename h0 ndur Dt nfiles pathfiles tide_zeta
 
end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%  LINK BOUND WAVES %%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 

time_tot=[0:Dt:(nfiles)*ndur*3600-Dt]';
bw_ent=zeros(size(time_tot));
env_ent=zeros(size(time_tot));
eta_ent=zeros(size(time_tot));
npo=(ndur*3600)/Dt;
tpo=(5*60)/Dt;
pond1=[0:1/tpo:(tpo-1)*1/tpo]';
pond2=1-pond1;

for fn=1:nfiles
     if fn==1 
        bw_ent(1:tpo)=bwtot_dur(1:tpo,fn);
        env_ent(1:tpo)=env_dur(1:tpo,fn);
        eta_ent(1:tpo)=eta_dur(1:tpo,fn);
     else
        bw_ent((fn-1)*npo+1:(fn-1)*npo+tpo)=bwtot_dur(end-tpo+1:end,fn-1).*pond2+bwtot_dur(1:tpo,fn).*pond1;
        env_ent((fn-1)*npo+1:(fn-1)*npo+tpo)=env_dur(end-tpo+1:end,fn-1).*pond2+env_dur(1:tpo,fn).*pond1;
        eta_ent((fn-1)*npo+1:(fn-1)*npo+tpo)=eta_dur(end-tpo+1:end,fn-1).*pond2+eta_dur(1:tpo,fn).*pond1;
     end     
        bw_ent((fn-1)*npo+tpo+1:(fn)*npo)=bwtot_dur(tpo+1:npo,fn);    
        env_ent((fn-1)*npo+tpo+1:(fn)*npo)=env_dur(tpo+1:npo,fn); 
        eta_ent((fn-1)*npo+tpo+1:(fn)*npo)=eta_dur(tpo+1:npo,fn); 
     if fn==file_num
        bw_ent(end-tpo+1:end)=bwtot_dur(npo+1:end,fn);
        env_ent(end-tpo+1:end)=env_dur(npo+1:end,fn);
        eta_ent(end-tpo+1:end)=eta_dur(npo+1:end,fn);
     else
        bw_ent(fn*npo+1:fn*npo+tpo)=bwtot_dur(npo+tpo:end,fn).*pond2+bwtot_dur(1:tpo,fn+1).*pond1;
        env_ent(fn*npo+1:fn*npo+tpo)=env_dur(npo+tpo:end,fn).*pond2+env_dur(1:tpo,fn+1).*pond1;
        eta_ent(fn*npo+1:fn*npo+tpo)=eta_dur(npo+tpo:end,fn).*pond2+eta_dur(1:tpo,fn+1).*pond1;
     end
end

plot(time_tot./3600,eta_ent)
hold on
plot(time_tot./3600,env_ent,'r')
plot(time_tot./3600,bw_ent,'k')
xlabel('Time,min')
ylabel('\eta, m')
legend('Bound wave')

env_ent1=1/8*1025*9.81*(2.*env_ent).^2;

xx=[time_tot bw_ent env_ent1];

fid = fopen('H:\My_Projects\FAR_INFRAGRAVITY\CASCAIS_SPECTRA\Cascais_bw_r5_1.txt','w');
fprintf(fid,'%12.2f  %12.8f %12.8f\n',xx');
fclose(fid);

save('H:\My_Projects\FAR_INFRAGRAVITY\CASCAIS_SPECTRA\realization5_1.mat','bwtot_dur','time_tot')

