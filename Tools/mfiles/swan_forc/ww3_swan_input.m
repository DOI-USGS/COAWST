% script ww3_swan_input.m
%
% This script is the main driver to 
% download, convert and interpolate WW3 data to TPAR input for SWAN.
%
% Written 4/23/09 by Brandy Armstrong
% some mods, jcwarner Arpil 27, 2009
%
% READ THE INSTRUCTIONS IN THE COAWST MANUAL 
% FOR SWAN BC's SECTION 10.
%
% First, acquire the necessary grib files from
% ftp://polar.ncep.noaa.gov/pub/history/waves/
%

% ************* BEGIN USER INPUT   ****************************

%1) Enter WORKING DIRECTORY.
% This is the location of ww3 grb files downloaded and the
% location of output for the TPAR files to be created.
% ***WARNING***
% The TPAR files created are saved in the working directory and are named
% generically (numbered). Any existing TPAR files will be overwritten !!!!
%
working_drive='d:';
working_dir='data\models\COAWST\Tools\mfiles\swan_forc\';
eval(['cd ',working_drive,'\',working_dir,';']);

%2) Enter dates of data requested.
yearww3='2005';    %input year of data yyyy 
mmww3='05';        %input month of data mm
ddww3='00';        %keep this as '00'

%3) Enter the ww3 grid area
ww3_area='multi_1.glo_30m';    %western north atlantic

%4) Enter path\name of SWAN netcdf grid. This is typically the same
% as the roms grid.
modelgrid='D:/data/Carolinas/modeling/Grids/USeast_grd17.nc';

%5) Enter the spacings of TPAR file locations around the perimeter
% of the grid. One TPAR file every 'specres' point.
% ww3_specpoints assumes masking of 0 for land and NaN for water
specres=20; % spec point resolution


% *************END OF USER INPUT****************************

% Call routine to compute TPAR files.
ww3gb_2TPAR(modelgrid,yearww3,mmww3,ww3_area,ddww3,specres)

% After creating the TPAR files, tell the user what info is needed to 
% be added to INPUT file.
% Write out boundary file lines for INPUT file
bdry_com %script writes out file Bound_spec_command to working directory
display('BOUNDSPEC command lines can be found in the file Bound_spec_command');
display('You must copy the lines from that file into your SWAN INPUT file');


