% function create_scrip_weights_master.m
%
% This m file is part of a series of routines driven by "create_scrip_weights_master.m"
%
% This routine is the Master driver m file to create sparse matrix interpolation 
% weights file using SCRIP.
%
% jcwarner 29June2014
%
% set up for Projects/Inlet_test/Refined
%

%%%%%%%  Begin USER Input Section  %%%%%%%%%%%%%

%step 1:  Enter the number of grids for each model:
Ngrids_roms=2;
Ngrids_swan=2;
Ngrids_wrf=0;

%step 2a: Enter the grids for the ROMS model:
roms_grids{1}='inlet_test_grid.nc';
roms_grids{2}='inlet_test_grid_ref3.nc';

%step 2b: Enter the grids for the SWAN model:
%Inlet_test_refined
 swan_coord{1}='inlet_test_grid_coord.grd';
 swan_coord{2}='inlet_test_grid_coord_ref3.grd';
 Numx{1}=77; Numy{1}=72;   % this was for inlet_test_grid_coord
 Numx{2}=92; Numy{2}=50;   % this was for inlet_test_grid_ref3_coord
 cartesian{1}=1;  % this was in meters
 cartesian{2}=1;  % this was in meters
 bath_file{1}='inlet_test_bathy.bot';
 bath_file{2}='inlet_test_bathy_ref3.bot';

%step 2c: Enter the grids for the WRF model:
wrf_grids{1}='wrfinput_d01';

%step 3: enter location of scrip.exe
scrip_exe='c:/work/models/COAWST/Lib/SCRIP/scrip.exe';

%step 4: enter working dir
wdir='c:\work\models\COAWST\Projects\Inlet_test\Refined';

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
  create_scrip_masks_driver_alexandra
end

if (create_scrip_files)
  create_scrip_files_driver_alexandra
end

if (check_scrip_weights)
  check_scrip_weights_driver
end

display('finished creating scrip weights')

