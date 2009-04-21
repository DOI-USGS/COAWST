function create_roms_netcdf_bndry(fn,gn,t_clim)
%
% create the boundary file for US East
% so this is special and only has certain arrays.
%
% jcw April 18, 2009
%
%

% Get some grid info. 
  [MP,LP]=size(gn.lon_rho);
  L=LP-1;
  Lm=L-1;
  M=MP-1;
  Mm=M-1;
  L  = Lm+1;
  M  = Mm+1;
  xi_psi  = L;
  xi_rho  = LP;
  xi_u    = L;
  xi_v    = LP;
  eta_psi = M;
  eta_rho = MP;
  eta_u   = MP;
  eta_v   = M;
  s       = gn.N;

%create bndry file
nc_bndry=netcdf(fn,'clobber');
 
%% Global attributes:
disp(' ## Defining Global Attributes...')
nc_bndry.history = ncchar(['Created by "' mfilename '" on ' datestr(now)]);
nc_bndry.type = ncchar('Initialization file from create_roms_boundary_bndrying.m');

%% Dimensions:

disp(' ## Defining Dimensions...')
 
nc_bndry('xi_psi') = L;
nc_bndry('xi_rho') = LP;
nc_bndry('xi_u') = L;
nc_bndry('xi_v') = LP;

nc_bndry('eta_psi') = M;
nc_bndry('eta_rho') = MP;
nc_bndry('eta_u') = MP;
nc_bndry('eta_v') = M;

nc_bndry('s_rho') = s;

nc_bndry('zeta_time')=length(t_clim);
nc_bndry('v2d_time')=length(t_clim);
nc_bndry('v3d_time')=length(t_clim);
nc_bndry('salt_time')=length(t_clim);
nc_bndry('temp_time')=length(t_clim);
 
%% Variables and attributes:
disp(' ## Defining Dimensions, Variables, and Attributes...')
 
nc_bndry{'zeta_time'} = ncdouble('zeta_time');
nc_bndry{'zeta_time'}.long_name = ncchar('zeta_time');
nc_bndry{'zeta_time'}.units = ncchar('seconds');
nc_bndry{'zeta_time'}.field = ncchar('zeta_time, scalar, series');

nc_bndry{'v2d_time'} = ncdouble('v2d_time');
nc_bndry{'v2d_time'}.long_name = ncchar('v2d_time');
nc_bndry{'v2d_time'}.units = ncchar('seconds');
nc_bndry{'v2d_time'}.field = ncchar('v2d_time, scalar, series');

nc_bndry{'v3d_time'} = ncdouble('v3d_time');
nc_bndry{'v3d_time'}.long_name = ncchar('v3d_time');
nc_bndry{'v3d_time'}.units = ncchar('seconds');
nc_bndry{'v3d_time'}.field = ncchar('v3d_time, scalar, series');

nc_bndry{'salt_time'} = ncdouble('salt_time');
nc_bndry{'salt_time'}.long_name = ncchar('salt_time');
nc_bndry{'salt_time'}.units = ncchar('seconds');
nc_bndry{'salt_time'}.field = ncchar('salt_time, scalar, series');

nc_bndry{'temp_time'} = ncdouble('temp_time');
nc_bndry{'temp_time'}.long_name = ncchar('temp_time');
nc_bndry{'temp_time'}.units = ncchar('seconds');
nc_bndry{'temp_time'}.field = ncchar('temp_time, scalar, series');

nc_bndry{'zeta_south'} = ncdouble('zeta_time','xi_rho');
nc_bndry{'zeta_south'}.long_name = ncchar('free-surface southern boundary condition');
nc_bndry{'zeta_south'}.units = ncchar('meter');
nc_bndry{'zeta_south'}.field = ncchar('zeta_south, scalar, series');

nc_bndry{'zeta_east'} = ncdouble('zeta_time','eta_rho');
nc_bndry{'zeta_east'}.long_name = ncchar('free-surface eastern boundary condition');
nc_bndry{'zeta_east'}.units = ncchar('meter');
nc_bndry{'zeta_east'}.field = ncchar('zeta_east, scalar, series');
 
nc_bndry{'ubar_south'} = ncdouble('v2d_time','xi_u');
nc_bndry{'ubar_south'}.long_name = ncchar('2D u-momentum southern boundary condition');
nc_bndry{'ubar_south'}.units = ncchar('meter second-1');
nc_bndry{'ubar_south'}.field = ncchar('ubar_south, scalar, series');

nc_bndry{'ubar_east'} = ncdouble('v2d_time','eta_u');
nc_bndry{'ubar_east'}.long_name = ncchar('2D u-momentum eastern boundary condition');
nc_bndry{'ubar_east'}.units = ncchar('meter second-1');
nc_bndry{'ubar_east'}.field = ncchar('ubar_east, scalar, series');
 
nc_bndry{'vbar_south'} = ncdouble('v2d_time','xi_v');
nc_bndry{'vbar_south'}.long_name = ncchar('2D v-momentum southern boundary condition');
nc_bndry{'vbar_south'}.units = ncchar('meter second-1');
nc_bndry{'vbar_south'}.field = ncchar('vbar_south, scalar, series');

nc_bndry{'vbar_east'} = ncdouble('v2d_time','eta_v');
nc_bndry{'vbar_east'}.long_name = ncchar('2D v-momentum eastern boundary condition');
nc_bndry{'vbar_east'}.units = ncchar('meter second-1');
nc_bndry{'vbar_east'}.field = ncchar('vbar_east, scalar, series');

nc_bndry{'u_south'} = ncdouble('v3d_time','s_rho','xi_u');
nc_bndry{'u_south'}.long_name = ncchar('3D u-momentum southern boundary condition');
nc_bndry{'u_south'}.units = ncchar('meter second-1');
nc_bndry{'u_south'}.field = ncchar('u_south, scalar, series');

nc_bndry{'u_east'} = ncdouble('v3d_time','s_rho','eta_u');
nc_bndry{'u_east'}.long_name = ncchar('3D u-momentum eastern boundary condition');
nc_bndry{'u_east'}.units = ncchar('meter second-1');
nc_bndry{'u_east'}.field = ncchar('u_east, scalar, series');

nc_bndry{'v_south'} = ncdouble('v3d_time','s_rho','xi_v');
nc_bndry{'v_south'}.long_name = ncchar('3D v-momentum southern boundary condition');
nc_bndry{'v_south'}.units = ncchar('meter second-1');
nc_bndry{'v_south'}.field = ncchar('v_south, scalar, series');

nc_bndry{'v_east'} = ncdouble('v3d_time','s_rho','eta_v');
nc_bndry{'v_east'}.long_name = ncchar('3D v-momentum eastern boundary condition');
nc_bndry{'v_east'}.units = ncchar('meter second-1');
nc_bndry{'v_east'}.field = ncchar('v_east, scalar, series');

nc_bndry{'temp_south'} = ncdouble('temp_time','s_rho','xi_rho');
nc_bndry{'temp_south'}.long_name = ncchar('3D temperature southern boundary condition');
nc_bndry{'temp_south'}.units = ncchar('degrees C');
nc_bndry{'temp_south'}.field = ncchar('temp_south, scalar, series');

nc_bndry{'temp_east'} = ncdouble('temp_time','s_rho','eta_rho');
nc_bndry{'temp_east'}.long_name = ncchar('3D temperature eastern boundary condition');
nc_bndry{'temp_east'}.units = ncchar('degrees C');
nc_bndry{'temp_east'}.field = ncchar('temp_east, scalar, series');

nc_bndry{'salt_south'} = ncdouble('salt_time','s_rho','xi_rho');
nc_bndry{'salt_south'}.long_name = ncchar('3D salinity southern boundary condition');
nc_bndry{'salt_south'}.units = ncchar('--');
nc_bndry{'salt_south'}.field = ncchar('salt_south, scalar, series');

nc_bndry{'salt_east'} = ncdouble('salt_time','s_rho','eta_rho');
nc_bndry{'salt_east'}.long_name = ncchar('3D salinity eastern boundary condition');
nc_bndry{'salt_east'}.units = ncchar('--');
nc_bndry{'salt_east'}.field = ncchar('salt_east, scalar, series');

%close file
close(nc_bndry)

