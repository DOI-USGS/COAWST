function create_inwave_ini(Lp,Mp,Nbins,Bindir,Bindir_c,pd,Ac,Cx,Cy,Ct,TA,inifile)

disp(' ')
disp(['## Creating the file : ',inifile])

initime=0;

%
%  Create the initial file
%

type = 'INITIAL file' ; 
history = 'InWave' ;
nc = netcdf(inifile,'clobber');
result = redef(nc);


L=Lp-1;
M=Mp-1;

%
%  Create dimensions
%

nc('xi_u') = L;
nc('xi_rho') = Lp;
nc('eta_v') = M;
nc('eta_rho') = Mp;
nc('energy_time') = length(initime);
nc('energy_angle') = Nbins+1;
nc('energy_angle_c') = Nbins;
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

nc{'energy_angle_c'} = ncdouble('energy_angle_c') ;
nc{'energy_angle_c'}.long_name = ncchar('direction respect to the north of the bin');
nc{'energy_angle_c'}.long_name = 'direction respect to the north of the bin';
nc{'energy_angle_c'}.units = ncchar('degrees');
nc{'energy_angle_c'}.units = 'degrees';

nc{'TA_dim'} = ncdouble('TA_dim') ;
nc{'TA_dim'}.long_name = ncchar('representative absolute peak period');
nc{'TA_dim'}.long_name = 'representative absolute peak period';
nc{'TA_dim'}.units = ncchar('Seconds');
nc{'TA_dim'}.units = 'Seconds';

nc{'AC'} = ncdouble('energy_angle_c','eta_rho','xi_rho') ;
nc{'AC'}.long_name = ncchar('wave energy envelope');
nc{'AC'}.long_name = 'wave energy envelope';
nc{'AC'}.units = ncchar('Joules');
nc{'AC'}.units = 'Joules';

nc{'cx'} = ncdouble('energy_angle_c','eta_rho','xi_u') ;
nc{'cx'}.long_name = ncchar('x component of the group celerity');
nc{'cx'}.long_name = 'x component of the group celerity';
nc{'cx'}.units = ncchar('m/s');
nc{'cx'}.units = 'm/s';

nc{'cy'} = ncdouble('energy_angle_c','eta_v','xi_rho') ;
nc{'cy'}.long_name = ncchar('y component of the group celerity');
nc{'cy'}.long_name = 'y component of the group celerity';
nc{'cy'}.units = ncchar('m/s');
nc{'cy'}.units = 'm/s';

nc{'ct'} = ncdouble('energy_angle','eta_rho','xi_rho') ;
nc{'ct'}.long_name = ncchar('directional component of the group celerity');
nc{'ct'}.long_name = 'directional component of the group celerity';
nc{'ct'}.units = ncchar('rad/s');
nc{'ct'}.units = 'rad/s';
  
nc{'Ta'} = ncdouble('TA_dim') ;
nc{'Ta'}.long_name = ncchar('Representative absolute peak period');
nc{'Ta'}.long_name = 'Representative absolute peak period';
nc{'Ta'}.units = ncchar('Seconds');
nc{'Ta'}.units = 'Seconds';

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
nc{'energy_angle'}(1:Nbins+1) = Bindir(1,1:Nbins+1); 
nc{'energy_angle_c'}(1:Nbins) = Bindir_c(1,1:Nbins); 
nc{'TA_dim'}(:) = 1; 


nc{'AC'}(:,:,:) = Ac(:,:,:); 
nc{'cx'}(:,:,:) = Cx(:,:,:); 
nc{'cy'}(:,:,:) = Cy(:,:,:); 
nc{'ct'}(:,:,:) = Ct(:,:,:); 
nc{'Ta'}(:) = TA;

close(nc)

disp(['Created initial file -->   ',inifile])

end




