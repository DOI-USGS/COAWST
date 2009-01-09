%PLOT_HOLE07_BATHY read in and plot flow and wave grid for 1st test of
%       Nearshore ROMS on one of J.List's Delft3D morphological evolution 
%       test cases

% Some grid parameters below are given below from the program that made the
%  grids.  The coast is oriented N/S, with the water towards the East.
%
% flow grid:
% out_name='shelf01_exp_hole7';
% xorig=400000;            %x-origin (a N.C. UTM value)
% yorig=3997800;           %y-origin (a N.C. UTM value at about 36 deg N)
%                          %xlength is determined from depth array
% ylength=8400;            %y length 
% deltax=20;               %delta x
% deltay=80;               %delta y
% 
% wave grid:
% out_name='shelf02_exp_hole7';
% xorig=400000;            %x-origin (a N.C. UTM value)
% yorig=3994200;           %y-origin (a N.C. UTM value at about 36 deg N)
%                          %xlength is determined from depth array
% ylength=15600;            %y length
% deltax=20;               %delta x
% deltay=80;               %delta y

%% Load the flow grid:
load shelf01_exp_hole7 

% view the surface:
figure
surf(x_array,y_array,-depth_array_mod)
daspect([1 1 .01])
shading interp

% contour the surface
figure
[C,h] = contourf(x_array,y_array,-depth_array_mod);
axis equal
clabel(C,h);

%% Load the wave grid:
load shelf02_exp_hole7 

% view the surface:
figure
surf(x_array,y_array,-depth_array_mod)
daspect([1 1 .01])
shading interp

% contour the surface:
figure
[C,h] = contourf(x_array,y_array,-depth_array_mod);
axis equal
clabel(C,h);

