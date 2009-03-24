% script create_griz_bay_forcing
%
% Create a netcdf file that contains forcing data for ROMS:
%
% 'zeta_east'  -   'free surface on east side at RHO-points'
% 'zeta_west'  -   'free surface on west side at RHO-points'
%
% jcw 2-18-2009
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Begin user input section                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%1) Enter name of netcdf forcing file to be created.
%   If it already exists it will be overwritten!!.
    forc_file='griz_bay_forc.nc';

%2) Enter name of grid file that forcing will be applied to.
    grid_file='C:\work\models\COAWST\Projects\Griz_Bay\griz_grid2.nc';

%3) create zeta_* time series
zeta_time=[0:900:10*24*3600];   %10 days
ramp=min(zeta_time./(3600*12),1);
zeta_east=-2*sin(2*pi.*(zeta_time-3600)/(12.42*3600)).*ramp;
zeta_west=-2*sin(2*pi.*zeta_time/(12.42*3600)).*ramp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  END of USER INPUT                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get some grid info. 
  ncload(grid_file);
  [MP,LP]=size(h);
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

%increase zeta_* to be size of grid
zeta_east=repmat(zeta_east,eta_rho,1)';
zeta_west=repmat(zeta_west,eta_rho,1)';

%create bndry file
nc_forc=netcdf(forc_file,'clobber');
 
%% Global attributes:
disp(' ## Defining Global Attributes...')
nc_forc.history = ncchar(['Created by "' mfilename '" on ' datestr(now)]);
nc_forc.type = ncchar('Initialization file from create_griz_bay_forcing.m');

%% Dimensions:

disp(' ## Defining Dimensions...')
 
nc_forc('xi_psi') = L;
nc_forc('xi_rho') = LP;
nc_forc('xi_u') = L;
nc_forc('xi_v') = LP;

nc_forc('eta_psi') = M;
nc_forc('eta_rho') = MP;
nc_forc('eta_u') = MP;
nc_forc('eta_v') = M;

nc_forc('zeta_time')=length(zeta_time);
 
%% Variables and attributes:
disp(' ## Defining Dimensions, Variables, and Attributes...')
 
nc_forc{'zeta_time'} = ncdouble('zeta_time');
nc_forc{'zeta_time'}.long_name = ncchar('zeta_time');
nc_forc{'zeta_time'}.units = ncchar('seconds');
nc_forc{'zeta_time'}.field = ncchar('zeta_time, scalar, series');

nc_forc{'zeta_east'} = ncdouble('zeta_time','eta_rho');
nc_forc{'zeta_east'}.long_name = ncchar('free-surface eastern boundary condition');
nc_forc{'zeta_east'}.units = ncchar('meter');
nc_forc{'zeta_east'}.field = ncchar('zeta_east, scalar, series');
 
nc_forc{'zeta_west'} = ncdouble('zeta_time','eta_rho');
nc_forc{'zeta_west'}.long_name = ncchar('free-surface western boundary condition');
nc_forc{'zeta_west'}.units = ncchar('meter');
nc_forc{'zeta_west'}.field = ncchar('zeta_west, scalar, series');
 
%now write the data from the arrays to the netcdf file
disp(' ## Filling Variables in netcdf file with data...')

nc_forc{'zeta_time'}(:) = zeta_time;
nc_forc{'zeta_west'}(:) = zeta_west;
nc_forc{'zeta_east'}(:) = zeta_east;

%close file
close(nc_forc)
