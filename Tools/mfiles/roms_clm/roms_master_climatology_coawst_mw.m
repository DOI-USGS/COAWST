
%
% Driver to create 3 files for roms: clim, bndry, and init.
%
% This is currently set up to use opendap calls to acquire data
% from HYCOM and interp to US_East grid.
%
% based on efforts by:
% written by Mingkui Li, May 2008
% Modified by Brandy Armstrong March 2009
% jcwarner April 20, 2009
%

warning off

%1) Enter start and end dates to get climatology data here.
%(time step is 1 day).
T1=floor(now+1); %start date

%2) Enter working directory here.
wdr='d:/data/models/COAWST/Tools/mfiles';
eval(['cd ',wdr])

%3) Enter path\grid_name here.
modelgrid='D:\data\Carolinas\modeling\Grids\USeast_grd17.nc'
eval(['gridname=''',modelgrid,''';']);

%4) Set grid vertical coordinate params here.
theta_s=5;
theta_b=0.4;
Tcline=50;
N=16;

%***********  END USER INPUT *******************************
disp('getting roms grid dimensions ...');
gn=roms_get_grid_mw(gridname,[theta_s theta_b Tcline N]);

tic

% call to create clm file
disp('going to create clm file')
fn=updatclim_coawst_mw(T1,gn,'coawst_clm.nc',wdr)

% call to create boundary file
disp('going to create bndry file')
% fn = filename from updatclim, we should change this later as we see fit.
updatbdry_coawst_mw(fn,gn,'coawst_bdy.nc',wdr)

% call to create init file
disp('going to create init file')
updatinit_coawst_mw(fn,gn,'coawst_ini.nc',wdr)
    
toc


