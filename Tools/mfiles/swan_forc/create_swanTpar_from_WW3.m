% script create_swanTpar_from_WW3.m
%
% This script is the main driver to create TPAR boundary files for SWAN
% from WW3 data.
%
%  !!!!   This data is only from Feb 2005 to May 2019  !!!!
%  https://www.ncei.noaa.gov/thredds-ocean/catalog.html
%

% ************* BEGIN USER INPUT   ****************************

% 1) Enter WORKING DIRECTORY.
% This is the location where the forcing files to be created.
%
working_dir='E:\data\models\COAWST'

% 2) Enter the year and month of the data requested.
yearww3='2012';    %input year of data yyyy 
mmww3='10';        %input month of data mm
%

% 3) Enter path\name of SWAN grid. This is set up to use the roms grid as the same for swan.
modelgrid='E:\data\models\COAWST\Projects\Sandy\Sandy_roms_grid.nc';

% 4) Enter the spacings of the forcing file locations around the perimeter
% of the grid. One forcings file spans between the 'specres' points.
specres=20; % spec point resolution

% 5)
% Enter the WW3 grid:
% wc_4m, wc_10m, glo_30m, ep_10m, at_4m, at_10m, ao_30m, ak_4m, ak_10m
% grid descriptions are here:  https://polar.ncep.noaa.gov/waves/implementations.php
%
ww3_grid='glo_30m';

% *************END OF USER INPUT****************************

%cd to user dir
eval(['cd ',working_dir,';']);

% call to get the spectral points
[specpts]=ww3_specpoints(modelgrid,specres);

% Call routine to compute TPAR files.
readww3_2TPAR(modelgrid,yearww3,mmww3,ww3_grid,specpts)

% After creating the TPAR files, tell the user what info is needed to 
% be added to INPUT file.
% Write out boundary file lines for INPUT file
bdry_com(specpts) %script writes out file Bound_spec_command to working directory
display('BOUNDSPEC command lines can be found in the file Bound_spec_command');
display('You must copy the lines from that file into your SWAN INPUT file');

