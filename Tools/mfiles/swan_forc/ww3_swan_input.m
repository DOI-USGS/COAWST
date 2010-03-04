% script ww3_swan_input.m
%
% This script is the main driver to 
% download, convert and interpolate WW3 data to TPAR input for SWAN.
%
% Written 4/23/09 by Brandy Armstrong
% some mods, jcwarner Arpil 27, 2009
%
% READ THE INSTRUCTIONS IN THE COAWST MANUAL 
% FOR SWAN BC's SECTION 11.
% If you are using this to run in hindcast mode, then
% it is assumed that you have already installed
% the njToolbox and the jar files from
% http://www2.msstate.edu/~skb12/nopp/njTbx/index.html
% and that you have acquired the grib files from
% ftp://polar.ncep.noaa.gov/pub/history/waves/
%
% For forecast mode, just plow ahead below.
%

% ************* BEGIN USER INPUT   ****************************

%1) Enter WORKING DIRECTORY.
% This is the location of ww3 grb files downloaded (if needed) and the
% location of output for the TPAR files to be created.
% ***WARNING***
% The TPAR files created are saved in the working directory and are named
% generically (numbered). Any existing TPAR files will be overwritten !!!!
%
working_drive='g:';
working_dir='data2\Carolinas\modeling\bc_ic\ww3_200312_sc56\';
eval(['cd ',working_drive,'\',working_dir,';']);

%2) Enter dates of data requested.
yearww3='2003';    %input year of data yyyy 
mmww3='12';        %input month of data mm

% For recent data, i.e. the last 7 days, input ddww3='dd'
% For HISTORICAL data, keep ddww3='00'
ddww3='00';        % dd

%3) Enter the ww3 grid area
% choices are: akw, enp, nah, nph, nww3, wna
ww3_area='wna';    %western north atlantic

%4) Enter path\name of SWAN netcdf grid. This is typically the same
% as the roms grid.
%modelgrid='G:\data2\Carolinas\modeling\Grids\USeast_grd15.nc';
modelgrid='G:\data2\Carolinas\modeling\Grids\sc_grid_56.nc';

%5) Enter the spacings of TPAR file locations around the perimeter
% of the grid. One TPAR file every 'specres' point.
% ww3_specpoints assumes masking of 0 for land and NaN for water
specres=20; % spec point resolution


% *************END OF USER INPUT****************************

% Call routine to compute TPAR files.
% ww3gb_2TPAR = for hindcast mode
% ww3nc_2TPAR = for forecast mode
if str2num(ddww3)==0; %assumes these are grb files
    ww3gb_2TPAR(modelgrid,yearww3,mmww3,working_dir,working_drive,ww3_area,ddww3,specres)
else
    ww3nc_2TPAR(modelgrid,yearww3,mmww3,working_dir,working_drive,ww3_area,ddww3,specres)
end

% After creating the TPAR files, tell the user what info is needed to 
% be added to INPUT file.
% Write out boundary file lines for INPUT file
bdry_com %script writes out file Bound_spec_command to working directory
display('BOUNDSPEC command lines can be found in the file Bound_spec_command');
display('You must copy the lines from that file into your SWAN INPUT file');


