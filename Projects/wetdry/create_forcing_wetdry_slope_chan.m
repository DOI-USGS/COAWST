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
zeta_times=[0:10:43200*1]/(3600*24);            %days
zeta_east=10.0*sin(2*pi/1*zeta_times)-10.0;   %m

%%%%%%%%%%%%%%%%%  END of USER SECTION %%%%%%%%%%%%%%%%%

%get some grid info here
h=ncread(grid_file,'h');
[LP,MP]=size(h);

% Create the Roms NetCDF file.
nc = netcdf.create(forc_name, 'clobber');
if isempty(nc), return, end
 
%% Global attributes:
disp(' ## Defining Global Attributes...')
 
netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'type','ROMS file');
netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'gridid','forc file');
netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'history',['Created by ', mfilename ', on ', datestr(now)]);
netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'title','ROMS Application')

L = LP-1;       % The psi dimension.
M = MP-1;       % The psi dimension.

disp(' ## Defining Dimensions...')
 
xi_psi = netcdf.defDim(nc,'xi_psi',L);
xi_rho = netcdf.defDim(nc,'xi_rho',LP);
xi_u   = netcdf.defDim(nc,'xi_u',L);
xi_v   = netcdf.defDim(nc,'xi_v',LP);

eta_psi = netcdf.defDim(nc,'eta_psi',M);
eta_rho = netcdf.defDim(nc,'eta_rho',MP);
eta_u   = netcdf.defDim(nc,'eta_u',MP);
eta_v   = netcdf.defDim(nc,'eta_v',M);

timedimID   = netcdf.defDim(nc,'zeta_time',length(zeta_times));

%% Variables and attributes:

disp(' ## Defining Variables and Attributes...')

zeta_timeID = netcdf.defVar(nc,'zeta_time','double',timedimID);
netcdf.putAtt(nc,zeta_timeID,'long_name','zeta_time');
netcdf.putAtt(nc,zeta_timeID,'units','days');
netcdf.putAtt(nc,zeta_timeID,'field','zeta_time, scalar, series');

zetaID = netcdf.defVar(nc,'zeta_east','float',[eta_rho timedimID]);
netcdf.putAtt(nc,zetaID,'long_name','free-surface');
netcdf.putAtt(nc,zetaID,'units','meter');
netcdf.putAtt(nc,zetaID,'field','free-surface, scalar, series');

netcdf.close(nc)

%now fill the arrays with data
ncwrite(forc_name,'zeta_time',zeta_times);
zz=(zeta_east'*ones(1,MP))';
ncwrite(forc_name,'zeta_east',zz);

