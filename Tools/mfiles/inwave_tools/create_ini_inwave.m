clear all
close all


inifile='F:\MAITANE\PROJECTS\INWAVE\InWave_ini.nc';

disp(' ')
disp([' Creating the file : ',inifile])
disp(' ')

Lp=101;   % xi_rho
Mp=101;   % eta_rho
Nangle=20;  % Number of desired directional bins; 

L=Lp-1;
M=Mp-1;

nbin=[1:1:Nangle];
dir_fin=(nbin-1).*360./(Nangle);
pd_fin=ones(size(dir_fin)).*360./(Nangle);

initime=0;

%
%  Create the initial file
%

type = 'INITIAL file' ; 
history = 'InWave' ;
nc = netcdf(inifile,'clobber');
result = redef(nc);


%
%  Create dimensions
%

nc('xi_u') = L;
nc('xi_rho') = Lp;
nc('eta_v') = M;
nc('eta_rho') = Mp;
nc('energy_time') = length(initime);
nc('energy_angle') = Nangle;
nc('TA_dim') = 1;

%
%  Create variables and attributes
%

nc{'energy_time'} = ncdouble('energy_time') ;
nc{'energy_time'}.long_name = ncchar('time for energy envelope');
nc{'energy_time'}.long_name = 'time for energy envelope';
nc{'energy_time'}.units = ncchar('seconds');
nc{'energy_time'}.units = 'seconds';

nc{'energy_angle'} = ncdouble('energy_angle') ;
nc{'energy_angle'}.long_name = ncchar('direction respect to the north of the bin');
nc{'energy_angle'}.long_name = 'direction respect to the north of the bin';
nc{'energy_angle'}.units = ncchar('degrees');
nc{'energy_angle'}.units = 'degrees';

nc{'TA_dim'} = ncdouble('TA_dim') ;
nc{'TA_dim'}.long_name = ncchar('representative absolute peak period');
nc{'TA_dim'}.long_name = 'representative absolute peak period';
nc{'TA_dim'}.units = ncchar('Seconds');
nc{'TA_dim'}.units = 'Seconds';

nc{'AC'} = ncdouble('energy_angle','energy_time','eta_rho','xi_rho') ;
nc{'AC'}.long_name = ncchar('wave energy envelope');
nc{'AC'}.long_name = 'wave energy envelope';
nc{'AC'}.units = ncchar('Joules');
nc{'AC'}.units = 'Joules';
  
nc{'TA'} = ncdouble('TA_dim') ;
nc{'TA'}.long_name = ncchar('Representative absolute peak period');
nc{'TA'}.long_name = 'Representative absolute peak period';
nc{'TA'}.units = ncchar('Seconds');
nc{'TA'}.units = 'Seconds';
  

%
% Create global attributes
%

nc.date = ncchar(date);
nc.date = date;
nc.clim_file = ncchar(inifile);
nc.clim_file = inifile;
nc.type = ncchar(type);
nc.type = type;

%
% Leave define mode
%

result = endef(nc);

%
% Write variables
%

nc{'energy_time'}(:) =  initime; 
nc{'energy_angle'}(:) = dir_fin; 
nc{'TA_dim'}(:) = 1; 


energy=zeros(length(initime),Nangle,Mp,Lp);
energy(:,16,:,:)=100;  % Including energy in the bin corresponding to 270 degrees

nc{'AC'}(:,:,:,:) = energy(:,:,:,:); 
nc{'TA'}(:) = 10;

close(nc)




