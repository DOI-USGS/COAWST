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
Numx{2}=116; Numy{2}=86;   % for swan grid 2
cartesian{1}=0;  % this was in lon/lat space
cartesian{2}=0;  % this was in lon/lat space
bath_file{1}='Sandy_swan_bathy.bot';
bath_file{2}='Sandy_swan_bathy_ref3.bot';

%step 2c: Enter the grids for the WRF model:
wrf_grids{1}='wrfinput_d01';
wrf_grids{2}='wrfinput_d02';

%step 3: Select to use scrip as a compiled executable or as matlab m files.
use_scrip_exe=1;      % put this = 1 to use scrip exe 
use_scrip_matlab=0;   % put this = 1 to use scrip m files
if (use_scrip_exe)
  scrip_exe='e:/data/models/COAWST/Lib/SCRIP/scrip.exe';
end

%step 4: enter working dir
wdir='c:\work\models\COAWST\Projects\Sandy2';

%step 5: Select the process steps to create the SCRIP files. If you are not
% sure, then just leave these 3 options to =1 so that they all run.

% step5a) Setting this flag to =1 calls "create_scrip_masks" that computes 
%  roms, swan, and wrf masks used by SCRIP. For refinement, these masks
%  allow multple grids to act as a single source to provide a 
%  combined data set. This also prevents land areas from participating as 
%  a source to the ocean or wave.
create_scrip_masks=1;

% step5b) Setting this flag to =1 calls "create_scrip_files" that calls 
%  scrip to compute the weights.
create_scrip_files=1;

% step5c) Setting this flag to =1 calls "check_scrip_weights" that calls 
%  an mfile to check that all the destination cells are getting information
%  from a source grid. This is currently only used for
% atm2ocn and atm2wav weights.
check_scrip_weights=1;

%%%%%%%  END USER Input Section  %%%%%%%%%%%%%

eval(['cd ',wdir])

if (create_scrip_masks)
  create_scrip_masks_driver
end

if (create_scrip_files)
  create_scrip_files_driver
end

if (check_scrip_weights)
  check_scrip_weights_driver
end

display('finished creating scrip weights')

