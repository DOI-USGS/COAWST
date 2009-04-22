function create_roms_netcdf_clm(fn,gn,tid1,tid2)

[eta xi]=size(gn.lon_rho);

nc=netcdf(fn,'clobber');
if isempty(nc), return, end

disp(' ## Defining Global Attributes...')
nc.history = ncchar(['Created by updatclim on ' datestr(now)]);
nc.type = ncchar('climate forcing file from http://hycom.coaps.fsu.edu:8080/thredds/dodsC/glb_analysis');

% Dimensions:

disp(' ## Defining Dimensions...')

LP=xi;
MP=eta;
L=LP-1;
M=MP-1;
N=0;

nc('xi_psi') = L;
nc('xi_rho') = LP;
nc('xi_u') = L;
nc('xi_v') = LP;

nc('eta_psi') = M;
nc('eta_rho') = MP;
nc('eta_u') = MP;
nc('eta_v') = M;
nc('s_rho') = gn.N;

nc('ocean_time')=length(tid1:tid2);
nc('zeta_time')=length(tid1:tid2);
nc('v2d_time')=length(tid1:tid2);
nc('v3d_time')=length(tid1:tid2);
nc('salt_time')=length(tid1:tid2);
nc('temp_time')=length(tid1:tid2);
nc('one') = 1;

% Variables and attributes:
disp(' ## Defining Dimensions, Variables, and Attributes...')

nc{'ocean_time'} = ncdouble('ocean_time'); %% 1 element.
nc{'ocean_time'}.long_name = ncchar('wind field time');
nc{'ocean_time'}.units = ncchar('days');
nc{'ocean_time'}.field = ncchar('wave_time, scalar, series');

nc{'zeta_time'} = ncdouble('zeta_time');
nc{'zeta_time'}.long_name = ncchar('zeta_time');
nc{'zeta_time'}.units = ncchar('days');
nc{'zeta_time'}.field = ncchar('zeta_time, scalar, series');

nc{'v2d_time'} = ncdouble('v2d_time');
nc{'v2d_time'}.long_name = ncchar('v2d_time');
nc{'v2d_time'}.units = ncchar('days');
nc{'v2d_time'}.field = ncchar('v2d_time, scalar, series');

nc{'v3d_time'} = ncdouble('v3d_time');
nc{'v3d_time'}.long_name = ncchar('v3d_time');
nc{'v3d_time'}.units = ncchar('days');
nc{'v3d_time'}.field = ncchar('v3d_time, scalar, series');

nc{'salt_time'} = ncdouble('salt_time');
nc{'salt_time'}.long_name = ncchar('salt_time');
nc{'salt_time'}.units = ncchar('days');
nc{'salt_time'}.field = ncchar('salt_time, scalar, series');

nc{'temp_time'} = ncdouble('temp_time');
nc{'temp_time'}.long_name = ncchar('temp_time');
nc{'temp_time'}.units = ncchar('days');
nc{'temp_time'}.field = ncchar('temp_time, scalar, series');

nc{'lon_rho'} = ncfloat('eta_rho', 'xi_rho');
nc{'lon_rho'}.long_name = ncchar('lon_rho');
nc{'lon_rho'}.units = ncchar('meter');
nc{'lon_rho'}.FillValue_ = ncfloat(100000.);
nc{'lon_rho'}.missing_value = ncfloat(100000.);
nc{'lon_rho'}.field = ncchar('xp, scalar, series');

nc{'lat_rho'} = ncfloat('eta_rho', 'xi_rho');
nc{'lat_rho'}.long_name = ncchar('lat_rho');
nc{'lat_rho'}.units = ncchar('meter');
nc{'lat_rho'}.missing_value = ncfloat(100000.);
nc{'lat_rho'}.FillValue_ = ncfloat(100000.);
nc{'lat_rho'}.field = ncchar('yp, scalar, series');

nc{'zeta'} = ncfloat('zeta_time','eta_rho','xi_rho');
nc{'zeta'}.long_name = ncchar('zeta');
nc{'zeta'}.units = ncchar('meter');
nc{'zeta'}.field = ncchar('zeta, scalar, series');

nc{'salt'} = ncfloat('salt_time','s_rho','eta_rho','xi_rho');
nc{'salt'}.long_name = ncchar('salt');
nc{'salt'}.units = ncchar('meter');
nc{'salt'}.field = ncchar('salt, scalar, series');

nc{'temp'} = ncfloat('temp_time','s_rho','eta_rho','xi_rho');
nc{'temp'}.long_name = ncchar('temp');
nc{'temp'}.units = ncchar('meter');
nc{'temp'}.field = ncchar('temp, scalar, series');

nc{'u'} = ncfloat('v3d_time','s_rho','eta_u','xi_u');
nc{'u'}.long_name = ncchar('velx');
nc{'u'}.units = ncchar('meter second-1');
nc{'u'}.field = ncchar('velx, scalar, series');

nc{'v'} = ncfloat('v3d_time','s_rho','eta_v','xi_v');
nc{'v'}.long_name = ncchar('vely');
nc{'v'}.units = ncchar('meter second-1');
nc{'v'}.field = ncchar('vely, scalar, series');

nc{'ubar'} = ncfloat('v2d_time','eta_u','xi_u');
nc{'ubar'}.long_name = ncchar('velx');
nc{'ubar'}.units = ncchar('meter second-1');
nc{'ubar'}.field = ncchar('velx, scalar, series');

nc{'vbar'} = ncfloat('v2d_time','eta_v','xi_v');
nc{'vbar'}.long_name = ncchar('vely');
nc{'vbar'}.units = ncchar('meter second-1');
nc{'vbar'}.field = ncchar('vely, scalar, series');

ncclose