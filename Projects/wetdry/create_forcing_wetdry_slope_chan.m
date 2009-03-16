% script create_forcing_wetdry_slope_chan
%
% create forcing file for roms
%
% for the wetdry sloping channel test case
%
% jcw 15Mar09
%

%1) enter name of the grid file here
grid_file='wetdry_slope_chan_grd.nc';

%2) enter name of output file
forc_name='wetdry_slope_chan_forc.nc';

%3) enter time series of zeta_east
zeta_time=[0:10:43200*1];                                     %secs
zeta_east=10.0*sin(2.0*pi/1.0*zeta_time*3600*24)-10.0+0.10;   %m

%%%%%%%%%%%%%%%%%  END of USER SECTION %%%%%%%%%%%%%%%%%

%get some grid info here
nc=netcdf(grid_file);
h=nc{'h'}(:);
[MP,LP]=size(h);
ncclose

% Create the Roms NetCDF file.

nc = netcdf(forc_name, 'clobber');
if isempty(nc), return, end
 
%% Global attributes:
disp(' ## Defining Global Attributes...')
 
nc.type = ncchar('Gridpak file');
nc.gridid = 'forc file';
nc.history = ncchar(['Created by "' mfilename '" on ' datestr(now)]);

L = LP-1;       % The psi dimension.
M = MP-1;       % The psi dimension.

disp(' ## Defining Dimensions...')
 
nc('xi_psi') = L;
nc('xi_rho') = LP;
nc('xi_u') = L;
nc('xi_v') = LP;

nc('eta_psi') = M;
nc('eta_rho') = MP;
nc('eta_u') = MP;
nc('eta_v') = M;

nc('zeta_time') = length(zeta_time);

%% Variables and attributes:

disp(' ## Defining Variables and Attributes...')
 
nc{'zeta_time'} = ncdouble('zeta_time');
nc{'zeta_time'}.long_name = ncchar('zeta_time');
nc{'zeta_time'}.units = ncchar('s');
nc{'zeta_time'}.field = ncchar('scalar');

nc{'zeta_east'} = ncdouble('zeta_time','eta_rho');
nc{'zeta_east'}.long_name = ncchar('free surface elevation');
nc{'zeta_east'}.units = ncchar('meter');
nc{'zeta_east'}.field = ncchar('m, scalar');
 
close(nc)

%now fill the arrays with data
nc = netcdf(forc_name, 'write')

nc{'zeta_time'}(:) = zeta_time;
nc{'zeta_east'}(:) = zeta_east'*ones(MP,1)';

ncclose

