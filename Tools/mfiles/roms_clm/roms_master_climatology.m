% create_roms_clim.m
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
tic

%1) Enter start and end dates to get climatology data here.
%(time step is 1 day).
T1=datenum(2003,12,1); %start date
T2=datenum(2003,12,1); %end date

%2) Enter working directory here.
wdr='g:\data2\Carolinas\modeling\bc_ic';

%3) Enter path/grid_name here.
gridname='g:\data2\Carolinas\modeling\Grids\USeast_grd13.nc';

%4) Set grid vertical coordinate params here.
theta_s=5;
theta_b=0.4;
Tcline=50;
N=16;

%***********  END USER INPUT *******************************
disp('getting roms grid dimensions ...');
gn=roms_get_grid(gridname,[theta_s theta_b Tcline N]);

% call to create climatology file
disp('going to create clm file')
updatclim

% call to create boundary file
disp('going to create bndry file')
% fn = filename from updatclim, we should change this later as we see fit.
updatbdry(fn,gn)

% call to create init file
disp('going to create init file')
updatinit(fn,gn)


