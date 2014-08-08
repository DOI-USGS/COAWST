% function create_scrip_weights_master.m
%
% This m file is part of a series of routines driven by "create_scrip_weights_master.m"
%
% This routine is the Master driver m file to create sparse matrix interpolation 
% weights file using SCRIP.
%
% jcwarner 29June2014
%
% specialized for Project/Sandy


%%%%%%%  Begin USER Input Section  %%%%%%%%%%%%%

%step 1:  Enter the number of grids for each model:
Ngrids_roms=2;
Ngrids_swan=2;
Ngrids_wrf=2;

%step 2a: Enter the grids for the ROMS model:
roms_grids{1}='Sandy_roms_grid.nc';
roms_grids{2}='Sandy_roms_grid_ref3.nc';

%step 2b: Enter the grids for the SWAN model:
%Sandy
swan_coord{1}='Sandy_swan_coord.grd';
swan_coord{2}='Sandy_swan_coord_ref3.grd';
Numx{1}=84;  Numy{1}=64;   % for swan grid 1
Numx{2}=122; Numy{2}=83;   % for swan grid 2
cartesian{1}=0;  % this was in lon/lat space
cartesian{2}=0;  % this was in lon/lat space
bath_file{1}='Sandy_swan_bathy.bot';
bath_file{2}='Sandy_swan_bathy_ref3.bot';

%step 2c: Enter the grids for the WRF model:
wrf_grids{1}='wrfinput_d01';
wrf_grids{2}='wrfinput_d02';

%step 3: enter location of scrip.exe
scrip_exe='e:/data/models/COAWST/Lib/SCRIP/scrip.exe';

%step 4: enter working dir
wdir='e:\data\models\COAWST\Projects\Sandy';

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

