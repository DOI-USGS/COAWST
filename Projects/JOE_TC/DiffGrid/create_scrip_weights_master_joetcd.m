% function create_scrip_weights_master.m
%
% This m file is part of a series of routines driven by "create_scrip_weights_master.m"
%
% This routine is the Master driver m file to create sparse matrix interpolation 
% weights file using SCRIP.
%
% jcwarner 29June2014
%


%%%%%%%  Begin USER Input Section  %%%%%%%%%%%%%

%step 1:  Enter the number of grids for each model:
Ngrids_roms=1;
Ngrids_swan=1;
Ngrids_wrf=1;

%step 2a: Enter the grids for the ROMS model:
roms_grids{1}='joe_tc_coarse_grd.nc';

%step 2b: Enter the grids for the SWAN model:
swan_coord{1}='joe_tc_coarse_grid_coord.grd';
Numx{1}=100;  Numy{1}=75;   % for swan grid 1
cartesian{1}=1;  % this was not in lon/lat space
bath_file{1}='joe_tc_coarse_roms_bathy.bot';

%step 2c: Enter the grids for the WRF model:
wrf_grids{1}='wrfinput_d01';

%step 3: enter location of scrip.exe
scrip_exe='e:/data/models/COAWST/Lib/SCRIP/scrip.exe';

%step 4: enter working dir
wdir='e:\data\models\COAWST\Projects\JOE_TCd';

%step 5: Select the process steps to create the SCRIP files.

% step5a) The "create_scrip_files" converts roms, swan, and wrf grids to 
% a format that SCRIP likes.
create_scrip_masks=1;

% step5b) The "create_scrip_weights" calls scrip to compute the weights 
% and produces the netcdf weights files.
create_scrip_files=1;

%%%%%%%  END USER Input Section  %%%%%%%%%%%%%

eval(['cd ',wdir])

if (create_scrip_masks)
  create_scrip_masks_driver
end

if (create_scrip_files)
  create_scrip_files_driver
end

display('finished creating scrip weights')

