function create_inwave_bnd(Lp, Mp, Nangle_bnd, Dir_bnd, obc, ...
    AcN, AcE, AcS, AcW, TA, time, bndfile)

disp(' ')
disp(['## Creating the file : ',bndfile])

L=Lp-1;
M=Mp-1;

%  Create the boundary file
%

type = 'BOUNDARY file for the InWave model' ; 
nc = netcdf(bndfile,'clobber');
result = redef(nc);

nbin=[1:1:Nangle_bnd];
dir=Dir_bnd;

%
%  Create dimensions
%

nc('xi_u') = L;
nc('xi_rho') = Lp;
nc('eta_v') = M;
nc('eta_rho') = Mp;
nc('energy_time') = length(time);
nc('energy_angle_c') = Nangle_bnd;
nc('TA_dim') = 1;

NT=length(time);

%
%  Create variables and attributes
%

nc{'energy_time'} = ncdouble('energy_time') ;
nc{'energy_time'}.long_name = ncchar('time for energy envelope');
nc{'energy_time'}.long_name = 'time for energy envelope';
nc{'energy_time'}.units = ncchar('seconds');
nc{'energy_time'}.units = 'seconds';

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

%

if obc(3)==1

%   Southern boundary

  nc{'AC_south'} = ncdouble('energy_time','energy_angle_c','xi_rho') ;
  nc{'AC_south'}.long_name = ncchar('southern boundary wave action envelope');
  nc{'AC_south'}.long_name = 'southern boundary wave action envelope';
  nc{'AC_south'}.units = ncchar('Joules');
  nc{'AC_south'}.units = 'Joules';
  
  nc{'TA_south'} = ncdouble('TA_dim') ;
  nc{'TA_south'}.long_name = ncchar('southern boundary representative absolute peak period');
  nc{'TA_south'}.long_name = 'southern boundary representative absolute peak period';
  nc{'TA_south'}.units = ncchar('Seconds');
  nc{'TA_south'}.units = 'Seconds';
  
end

if obc(2)==1

%   Eastern boundary

  nc{'AC_east'} = ncdouble('energy_time','energy_angle_c','eta_rho') ;
  nc{'AC_east'}.long_name = ncchar('eastern boundary wave action envelope');
  nc{'AC_east'}.long_name = 'eastern boundary wave action envelope';
  nc{'AC_east'}.units = ncchar('Joules');
  nc{'AC_east'}.units = 'Joules';
    
  nc{'TA_east'} = ncdouble('TA_dim') ;
  nc{'TA_east'}.long_name = ncchar('eastern boundary representative absolute peak period');
  nc{'TA_east'}.long_name = 'eastern boundary representative absolute peak period';
  nc{'TA_east'}.units = ncchar('Seconds');
  nc{'TA_east'}.units = 'Seconds';

end

if obc(1)==1
%
%   Northern boundary
%
  nc{'AC_north'} = ncdouble('energy_time','energy_angle_c','xi_rho') ;
  nc{'AC_north'}.long_name = ncchar('northern boundary wave action envelope');
  nc{'AC_north'}.long_name = 'northern boundary wave action envelope';
  nc{'AC_north'}.units = ncchar('Joules');
  nc{'AC_north'}.units = 'Joules';
  
  nc{'TA_north'} = ncdouble('TA_dim') ;
  nc{'TA_north'}.long_name = ncchar('northern boundary representative absolute peak period');
  nc{'TA_north'}.long_name = 'northern boundary representative absolute peak period';
  nc{'TA_north'}.units = ncchar('Seconds');
  nc{'TA_north'}.units = 'Seconds';
%
end
%
if obc(4)==1
%
%   Western boundary
%
  nc{'AC_west'} = ncdouble('energy_time','energy_angle_c','eta_rho') ;
  nc{'AC_west'}.long_name = ncchar('western boundary wave action envelope');
  nc{'AC_west'}.long_name = 'western boundary wave action envelope';
  nc{'AC_west'}.units = ncchar('Joules');
  nc{'AC_west'}.units = 'Joules';
  
  nc{'TA_west'} = ncdouble('TA_dim') ;
  nc{'TA_west'}.long_name = ncchar('westhern boundary representative absolute peak period');
  nc{'TA_west'}.long_name = 'westhern boundary representative absolute peak period';
  nc{'TA_west'}.units = ncchar('Seconds');
  nc{'TA_west'}.units = 'Seconds';
  

end
%
%
% Create global attributes
%

nc.date = ncchar(date);
nc.date = date;
nc.clim_file = ncchar(bndfile);
nc.clim_file = bndfile;
nc.type = ncchar(type);
nc.type = type;

%
% Leave define mode
%

result = endef(nc);

%
% Write variables
%

nc{'energy_time'}(1:NT) =  time(1,1:NT); 
nc{'energy_angle_c'}(:) =  dir; 
nc{'TA_dim'}(:) = 1; 

if obc(3)==1
   nc{'AC_south'}(:,:,:) =AcS(:,:,:);
   nc{'TA_south'}(:)=TA;
end

if obc(2)==1
   nc{'AC_east'}(:,:,:) =AcE(:,:,:); 
   nc{'TA_east'}(:)=TA;
end

if obc(1)==1
   nc{'AC_north'}(:,:,:) =AcN(:,:,:); 
   nc{'TA_north'}(:)=TA;
end

clear AcN AcE AcS Ac_east Ac Ac1

if obc(4)==1
   nc{'AC_west'}(:,:,:) =AcW(:,:,:); 
   nc{'TA_west'}(:)=TA;
end

close(nc)

disp(['Created boundary file -->   ',bndfile])
disp(' ')

return

