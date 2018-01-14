% create_sandy_application
%
% jcwarner 19August2014
%          29July2016 - general updates for workshop.
%          10 Jan 2018 - add WAVEWATCH III coupling
%
% This is a list of procedures and notes used to create the 
% Projects/Sandy files for roms, swan, and ww3. The wrf files
% are assumed to already be created.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% cd Projects/Sandy
% This is a place to create all the necessary files.
%
% Step 0.  We are drawing heavily on the work of the scientists at
% NCAR who already created the WRF input files. These were provided
% by Dave Gill and Cindy Bruyere. They provided a link to a WRF
% tutorial here:
%
% If you are interested in the WRF-only processing (allowing you to define 
% the atmospheric domain), please take a look at the following web page:
% http://www2.mmm.ucar.edu/wrf/OnLineTutorial/compilation_tutorial.php
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 1:  ROMS grid
%
% As a first cut, you can use the wrf grid. 
% Load the wrf grid to get coastline data
netcdf_load('wrfinput_d01')
figure
pcolorjw(XLONG,XLAT,double(1-LANDMASK))
hold on
title('WRF LANDMASK grid, Lambert conformal proj')
xlabel('longitude'); ylabel('latitiude')
%
% pick 4 corners for roms parent grid
Istr=27; Iend=84;
Jstr=7 ; Jend=64;
plot(XLONG(Istr,Jstr),XLAT(Istr,Jstr),'g+')
plot(XLONG(Iend,Jstr),XLAT(Iend,Jstr),'g+')
plot(XLONG(Istr,Jend),XLAT(Istr,Jend),'g+')
plot(XLONG(Iend,Jend),XLAT(Iend,Jend),'g+')

% pick number of points in the grid
numx=84;
numy=64;

% Make a matrix of the lons and lats
% You should not need to change this.
% start in lower left and go clockwise
corner_lon=double([XLONG(Istr,Jstr) XLONG(Istr,Jend) XLONG(Iend,Jend) XLONG(Iend,Jstr) ]);
corner_lat=double([ XLAT(Istr,Jstr)  XLAT(Istr,Jend)  XLAT(Iend,Jend)  XLAT(Iend,Jstr) ]);
x=[1 numx; 1 numx]; y=[1 1; numy numy];
z=[corner_lon(1) corner_lon(2) corner_lon(4) corner_lon(3)];
F = TriScatteredInterp(x(:),y(:),z(:));
[X,Y]=meshgrid(1:numx,1:numy);
lon=F(X,Y).';
%
x=[1 numx; 1 numx]; y=[1 1; numy numy];
z=[corner_lat(1) corner_lat(2) corner_lat(4) corner_lat(3)];
F = TriScatteredInterp(x(:),y(:),z(:));
lat=F(X,Y).';
plot(lon,lat,'k-')
plot(lon',lat','k-')
text(-75,27,'- - - roms grid')
%
% Call generic grid creation
%
roms_grid='Sandy_roms_grid.nc';
rho.lat=lat;
rho.lon=lon;
rho.depth=zeros(size(rho.lon))+100;  % for now just make zeros
rho.mask=zeros(size(rho.lon));       % for now just make zeros
spherical='T';
projection='lambert conformal conic';
%projection='mercator';
save temp_jcw33.mat rho spherical projection
eval(['mat2roms_mw(''temp_jcw33.mat'',''',roms_grid,''');'])
!del temp_jcw33.mat
%User needs to edit roms variables
disp(['     '])
disp(['Created roms grid -->   ',roms_grid])
disp(['     '])
disp(['You need to edit depth in ',roms_grid])
disp(['     '])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 2: Masking
%
% base this on the WRF mask.
F = TriScatteredInterp(double(XLONG(:)),double(XLAT(:)), ...
                       double(1-LANDMASK(:)),'nearest');
roms_mask=F(lon,lat);
figure
pcolorjw(lon,lat,roms_mask)
title('ROMS 1-LANDMASK grid, Lambert conformal proj')
xlabel('longitude'); ylabel('latitiude')
%
% compute mask on rho, u, v, and psi points.
water = double(roms_mask);
u_mask = water(1:end-1,:) & water(2:end,:);
v_mask= water(:,1:end-1) & water(:,2:end);
psi_mask= water(1:end-1,1:end-1) & water(1:end-1,2:end) & water(2:end,1:end-1) & water(2:end,2:end);
ncwrite('Sandy_roms_grid.nc','mask_rho',roms_mask);
ncwrite('Sandy_roms_grid.nc','mask_u',double(u_mask));
ncwrite('Sandy_roms_grid.nc','mask_v',double(v_mask));
ncwrite('Sandy_roms_grid.nc','mask_psi',double(psi_mask));
%
% use coastline tool to get coast:
% use GEODAS software http://www.ngdc.noaa.gov/mgg/geodas/geodas.html
% select lat bounds of 45 / 25 lon of -84 / -60 
%
lon=coastline(:,1);
lat=coastline(:,2);
save coastline.mat lon lat
%
% To create better mask, use editmask with the 
% grid file Sandy_roms_grid.nc and 
% the coastline file coastline.mat
%
editmask
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 3: Bathymetry
%
% You can use a source like this:
% http://www.ngdc.noaa.gov/mgg/global/global.html
% i am using data from a local file
load USeast_bathy.mat
netcdf_load(roms_grid)
h=griddata(h_lon,h_lat,h_USeast,lon_rho,lat_rho);
h(isnan(h))=5;
%smooth h a little
h(2:end-1,2:end-1)=0.2*(h(1:end-2,2:end-1)+h(2:end-1,2:end-1)+h(3:end,2:end-1)+ ...
                        h(2:end-1,1:end-2)+h(2:end-1,3:end));
% Bathymetry can be smoothed using
% http://drobilica.irb.hr/~mathieu/Bathymetry/index.html 
figure
pcolorjw(lon_rho,lat_rho,h)
hold on
load coastline.mat
plot(lon,lat,'k')
caxis([5 2500]); colorbar
title('ROMS bathy')
xlabel('longitude'); ylabel('latitiude')
ncwrite(roms_grid,'h',h);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 4: roms child grid
%
% plot wrf d01 and d02 grids
netcdf_load('wrfinput_d01')
figure
pcolorjw(XLONG,XLAT,double(1-LANDMASK))
hold on
netcdf_load('wrfinput_d02')
pcolorjw(XLONG,XLAT,double(1-LANDMASK))
plot(XLONG(1,:),XLAT(1,:),'r'); plot(XLONG(end,:),XLAT(end,:),'r')
plot(XLONG(:,1),XLAT(:,1),'r'); plot(XLONG(:,end),XLAT(:,end),'r')
% plot roms parent grid
plot(lon_rho(1,:),lat_rho(1,:),'k'); plot(lon_rho(end,:),lat_rho(end,:),'k')
plot(lon_rho(:,1),lat_rho(:,1),'k'); plot(lon_rho(:,end),lat_rho(:,end),'k')
text(-75,29,'roms parent grid')
text(-77,27,'wrf parent grid')
text(-77.2,34,'wrf child grid')
title('LANDMASKS')
xlabel('longitude'); ylabel('latitiude')
%
% Select child indices
Istr=22; Iend=60; Jstr=26; Jend=54;
%now make some plots and call the coarse-> fine
plot(lon_rho(Istr,Jstr),lat_rho(Istr,Jstr),'m+')
plot(lon_rho(Istr,Jend),lat_rho(Istr,Jend),'m+')
plot(lon_rho(Iend,Jstr),lat_rho(Iend,Jstr),'m+')
plot(lon_rho(Iend,Jend),lat_rho(Iend,Jend),'m+')
ref_ratio=3;
roms_child_grid='Sandy_roms_grid_ref3.nc';
F=coarse2fine('Sandy_roms_grid.nc','Sandy_roms_grid_ref3.nc', ...
              ref_ratio,Istr,Iend,Jstr,Jend);
Gnames={'Sandy_roms_grid.nc','Sandy_roms_grid_ref3.nc'};
[S,G]=contact(Gnames,'Sandy_roms_contact.nc');
%
% Recompute h.
%
netcdf_load('Sandy_roms_grid_ref3.nc')
load USeast_bathy.mat
h=griddata(h_lon,h_lat,h_USeast,lon_rho,lat_rho);
h(isnan(h))=5;
h=h;
h(2:end-1,2:end-1)=0.2*(h(1:end-2,2:end-1)+h(2:end-1,2:end-1)+h(3:end,2:end-1)+ ...
                        h(2:end-1,1:end-2)+h(2:end-1,3:end));
figure
pcolorjw(lon_rho,lat_rho,h)
hold on
plot(lon,lat,'r')
caxis([5 2500]); colorbar
title('ROMS bathy')
xlabel('longitude'); ylabel('latitiude')
ncwrite('Sandy_roms_grid_ref3.nc','h',h);
%
% recompute child mask based on WRF mask
%
netcdf_load('wrfinput_d02');
F = TriScatteredInterp(double(XLONG(:)),double(XLAT(:)), ...
                       double(1-LANDMASK(:)),'nearest');
roms_mask=F(lon_rho,lat_rho);
figure
pcolorjw(lon_rho,lat_rho,roms_mask)
title('ROMS child mask')
xlabel('longitude'); ylabel('latitiude')
hold on
plot(lon,lat,'r')
water = double(roms_mask);
u_mask = water(1:end-1,:) & water(2:end,:);
v_mask= water(:,1:end-1) & water(:,2:end);
psi_mask= water(1:end-1,1:end-1) & water(1:end-1,2:end) & water(2:end,1:end-1) & water(2:end,2:end);
ncwrite('Sandy_roms_grid_ref3.nc','mask_rho',roms_mask);
ncwrite('Sandy_roms_grid_ref3.nc','mask_u',double(u_mask));
ncwrite('Sandy_roms_grid_ref3.nc','mask_v',double(v_mask));
ncwrite('Sandy_roms_grid_ref3.nc','mask_psi',double(psi_mask));
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 5: 3D BC's
%
% Need to download nctoolbox from
% http://code.google.com/p/nctoolbox/
setup_nctoolbox
%
% Create roms init conditions, bc's, and nudging files for the parent grid.
% I copied Tools/mfiles/roms_clm/roms_master_climatology_coawst_mw.m
% to Projects/Sandy/roms_master_climatology_sandy.m and ran it:
roms_master_climatology_sandy
% this created three files:
% merged_coawst_clm.nc
% merged_coawst_bdy.nc
% coawst.ini.nc
% I renamed those to Sandy_clm.nc, Sandy_bdy.nc, and Sandy_ini.nc.
%
% For child grid, you need to create init and clm file. You can use the same 
% procedure as above using the master_climatology file, or use:
create_roms_child_init('Sandy_roms_grid.nc','Sandy_roms_grid_ref3.nc', ...
                       'Sandy_ini.nc','Sandy_ini_ref3.nc')
%and
create_roms_child_clm('Sandy_roms_grid.nc','Sandy_roms_grid_ref3.nc', ...
                       'Sandy_clm.nc','Sandy_clm_ref3.nc')
% that will interpolate the data from the parent.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 6: 2D BC's
%
% Get the tidal data at
svn checkout https://coawstmodel.sourcerepo.com/coawstmodel/data .
% edit Tools/mfiles/tides/create_roms_tides to select tidal database
% and start time. then run that file
create_roms_tides
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 7: Surface forcings.
%
% To create ROMS surface forcing files, there are several options such as:
% - create_roms_forcings.m: converts data from matlab workspce to a netcdf forcing file.
% - narrnc2roms.m: you need to get files from ftp://ftp.cdc.noaa.gov/Datasets/NARR/monolevel/
%                  and then run the mfile to create netcdf forcing file.
% - nam_narr_2roms.m: same as narrnc2roms, but the data is read via THREDDS.
%                     Suggest you try this one first.
nam_narr_2roms
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 8: roms input files.
%
edit sandy.h
% some of model setup cppdefs:
#define ROMS_MODEL
#define NESTING
#undef  WRF_MODEL
#undef  SWAN_MODEL
#undef  MCT_LIB
#undef  MCT_INTERP_OC2AT
#undef  MCT_INTERP_WV2AT
#undef  MCT_INTERP_OC2WV

%and
edit ocean_sandy.in
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 9: build file
%
edit and run coawst.bash
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 10: run file
%
edit and run run_coawst
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  ROMS - WRF coupling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 11: create SCRIP weights
% cd to LIib/SCRIP_COAWST and run
./scrip_coawst scrip_coawst_sandy.in
% to create all the scrip interpolation weights
% You can use scrip_coawst_sandy.in file to create weights for a moving 
% sandy nest or a static nest.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 12: input files.
%
edit sandy.h
% some of model setup cppdefs:
#define ROMS_MODEL
#define NESTING
#define WRF_MODEL
#undef  SWAN_MODEL
#define MCT_LIB
#define MCT_INTERP_OC2AT
#undef  MCT_INTERP_WV2AT
#undef  MCT_INTERP_OC2WV

%and
edit ocean_sandy.in
%and
edit coupling_sandy.in
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 13: build file
%
edit and run coawst.bash 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 14: run file
%
copy wrf* and namelist.input from Projects/Sandy to the root directory.
edit run_coawst
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  SWAN only application
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 15:  create swan grids: use roms2swan
roms2swan('Sandy_roms_grid.nc')
% this created swan_coord.grd and swan_bathy.bot and i renamed them to be
% Sandy_swan_bathy.bot and Sandy_swan_coord.grd
roms2swan('Sandy_roms_grid_ref3.nc')
% this created swan_coord.grd and swan_bathy.bot and i renamed them to be
% Sandy_swan_bathy_ref3.bot and Sandy_swan_coord_ref3.grd
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 16 - create SWAN forcing files
% To create a surface wind forcing file for SWAN, similar to nam_narr_2roms, use
nam_narr_2swan
% This m-file uses THREDDS to access data via internet to create the forcing file.
% The output file is set to be swan_narr_Oct2012.dat.  Make a copy of this for
% the child swan grid.
copy swan_narr_Oct2012.dat swan_narr_Oct2012_ref3.dat
% In the swan INPUT file, you need to set the WIND Inpgrid commnand line
% to be consistent with the data file created. The INPGRID line would look like
%
&& KEYWORD TO CREATE WIND GRID &&
INPGRID WIND REGULAR -105 10 0 550 400 0.1 0.1 &
       NONSTATIONARY 20121028.000000 3 HR 20121031.000000
READINP WIND 1 'Projects/Sandy/swan_narr_Oct2012.dat' 4 0 FREE
%
% These values are from the grid size  and time stamps:
%  lon_rho=[255:0.1:310]-360;
%  lat_rho=[ 10:0.1:50 ];  % Create a 0.1 degree lon-lat grid
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 17 - create SWAN boundary TPAR files
%
% READ THE INSTRUCTIONS IN THE COAWST MANUAL 
% FOR SWAN BC's SECTION 10.
% Acquire the necessary grib [hs, tp, and dp] files from
% ftp://polar.ncep.noaa.gov/pub/history/waves/
%  I got these 3 files:
%  multi_1.at_10m.tp.201210.grb2
%  multi_1.at_10m.hs.201210.grb2
%  multi_1.at_10m.dp.201210.grb2
%
% then run the ww3_swan_input m-file to create the TPAR files
ww3_swan_input
% - this file prodcues TPAR files. Place the TPAR files in 
%  your project folder. 
% - Also, the m files creates some commands lines called
% "BOUNDSPEC." You must copy the command lines from the file 
% Bound_spec_command and place them into your SWAN INPUT file,
% such as swan_sandy.in. This is only done for the parent grid.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 18 - create SWAN init files.
% To create init file for swan, you can run SWAN in stationary mode
% and create a ‘hot start’ file(s). 
edit sandy.h
% some of model setup cppdefs:
#undef  ROMS_MODEL
#define NESTING
#undef  WRF_MODEL
#define SWAN_MODEL
#undef  MCT_LIB
#undef  MCT_INTERP_OC2AT
#undef  MCT_INTERP_WV2AT
#undef  MCT_INTERP_OC2WV
% and rerun coawst.bash to build a new executable.
%
% Edit the swan_sandy.in files to just have INIT 
%
& Restart name **********************************
INIT
%
% and set to compute stationary and create a hotfile.
%
COMPUTE STAT 20121028.120000
HOTFILE 'Sandy_init.hot'
%
% do the same for the sandy ref3 input file.
%
% To run SWAN and create the init files use:
mpirun -np 16 -machinefile $PBS_NODEFILE ./coawstM Projects/Sandy/swan_sandy.in Projects/Sandy/swan_sandy_ref3.in > cwstv3.out
% this step created the *.hot files. 
% Copy those files to the Projects/Sandy folder.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 19 - run SWAN by itself for multiday simulation.
% To run SWAN for the multi day simulation, edit the swan.in files and list the init files:
& Restart name **********************************
INITIAL HOTSTART SINGLE 'Projects/Sandy/Sandy_init.hot'
&INIT
%
% and change the run to be nonstationary
%
COMPUTE NONSTAT 20121028.120000 180 SEC 20121030.120000
%
% notice the command to create hourly restart files:
RESTART 'swan_rst.dat' FREE 1 HR
%
% then run the program
mpirun -np 16 -machinefile $PBS_NODEFILE ./coawstM Projects/Sandy/swan_sandy.in Projects/Sandy/swan_sandy_ref3.in > cwstv3.out
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  ROMS - WRF - SWAN coupling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Step 20 - to run ROMS + WRF + SWAN
edit sandy.h
% some of model setup cppdefs:
#define ROMS_MODEL
#define NESTING
#define WRF_MODEL
#define SWAN_MODEL
#define MCT_LIB
#define MCT_INTERP_OC2AT
#define MCT_INTERP_WV2AT
#define MCT_INTERP_OC2WV
%
% recompile and run
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  WAVEWATCH III only
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 21 - I added text to the COAWST USer manual Section 11 for WW3 by itself.
% please read that part of the manual.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  ROMS - WRF - WAVEWATCH III coupling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Step 22 - to run ROMS + WRF + WW3
edit sandy.h
% some of model setup cppdefs:
#define ROMS_MODEL
#define NESTING
#define WRF_MODEL
#define WW3_MODEL
#undef  SWAN_MODEL
#define MCT_LIB
#define MCT_INTERP_OC2AT
#define MCT_INTERP_WV2AT
#define MCT_INTERP_OC2WV
%
%Step 23 edit coawst.bash
% Set these 5 environment variables:
#1) COAWST_WW3_DIR is a pointer to root WW3 code, do not change.
export   COAWST_WW3_DIR=${MY_ROOT_DIR}/WW3

#2) WWATCH3_NETCDF can be NC3 or NC4. We need NC4 for COAWST. do not change.
export   WWATCH3_NETCDF=NC4

#3) WWATCH_ENV points to WW3 environment listing. do not change.
export   WWATCH_ENV=${COAWST_WW3_DIR}/wwatch.env

#4) NETCDF_CONFIG is needed by WW3. You need to set this:
export   NETCDF_CONFIG=/usr/bin/nc-config

#5) WW3_SWITCH_FILE is like cpp options for WW3. You need to create it and
#        list the name here. This is described below.
export  WW3_SWITCH_FILE=sandy_coupled

3) As with all applications, set the MPI and compiler information (this will depend on your system):
export           USE_MPI=on            # distributed-memory parallelism
export        USE_MPIF90=on            # compile with mpif90 script
export         which_MPI=openmpi       # compile with OpenMPI library
export              FORT=ifort

4) Set the location of the header and analytical files:
export     MY_HEADER_DIR=${MY_PROJECT_DIR}/Projects/Sandy
export MY_ANALYTICAL_DIR=${MY_PROJECT_DIR}/Projects/Sandy

5) You may need to edit some WW3 specific build files.  So if you selected ifort (for example), then we need to look at
COAWST/WW3/bin/comp.Intel
COAWST/WW3/bin/link.Intel
Please look at these fiels (or the files needed for your compiler) to see if the flags are correct.

6) Create a switch file.
The switch file is like the cppdefs file for roms. It tells WW3 what features to compile. We have a file COAWST/WW3/bin/switch_sandy_coupled. List part of the name in the coawst.bash file. The file needs to start with ‘switch_’ and list the last part in the coawst.bash file. Here is the switch_sandy_coupled:

F90 NOGRB COAWST LRB4 SCRIP SCRIPNC NC4 TRKNC DIST MPI PR3 UQ FLX0 LN1 ST4 STAB0 NL1 BT4 DB1 MLIM TR0 BS0 IC2 IS2 REF1 IG1 XX0 WNT2 WNX1 RWND CRT1 CRX1 TIDE O0 O1 O2 O2a O2b O2c O3 O4 O5 O6 O7

Please read the WW3 manual in the WW3 directory. Some important options are:
NOPA	(not listed above) -this is for stand alone. If couled to another model, use COAWST and not NOPA. Do not use PALM.
COAWST	list this if WW3 coupled to another model, and then do not list NOPA.
SCRIP, SCRIPNC	I don’t use the WW3 scrip but we may in the future. So these are not really needed.
DIST, MPI	this is to use mpi, needed.

7) compile the code by running ./coawst.bash at the command prompt.
If it compiles correctly, you will get a coawstM.

8) Now we need to create some WW3 grids files. Here is a way to do that. There are probalbly many other ways to do this.
You can use Tools/mfiles/mtools/create_ww3_grid_files to create an x, y, bath, and mask files for WW3. This m file uses the roms grid to create the WW3 files of:
ww3_sandy_xcoord.dat, ww3_sandy_ycoord.dat, ww3_sandy_bathy.bot, and ww3_sandy_mapsta.inp.

9) Forcinng files for WW3: again, there are probably many ways to do this, but here is one way.
This is only needed if you are not coupled to WRF. You can use Tools/mfiles/mtools/create_ww3_wind_forcing to create a wind forcing file ww3_sandy_winnd_forc.dat.

10) We then need to create some WW3 run files using 4 of their utilities:
- edit WW3/work/ww3_grid.inp. Here is where you enter number of bins, advection choices, time steps (we chose 180), and the grid settings 
     'CURV'  T  'NONE'
      84     64
      40 1. 0. 1 1 'FREE' 'NAME' 'ww3_sandy_xcoord.dat'
      40 1. 0. 1 1 'FREE' 'NAME' 'ww3_sandy_ycoord.dat'
      -1.0  -1.0  40 -1. 1 1 '(....)' 'NAME' 'ww3_sandy_bathy.bot'
Then cd to WW3/work and run ../exe/ww3_grid to create mask.ww3, mapsta.ww3, mod_def.ww3.

- edit WW3/work/ww3_strt.inp. Here is where you enter init information.
Then cd to WW3/work and run ../exe/ww3_strt to create restart.ww3.

- edit WW3/work/ww3_prep.inp. Here is where you enter forcing data. You do not need this if it is couled to WRF. We have:
$
   'WND' 'LL' T T
$ Additional input format type 'LL' ---------------------------------- $
$ Grid range (degr. or m) and number of points for axes, respectively.
$ Example for longitude-latitude grid.
    -105. -50. 551   10. 50. 401
This tells WW3 that the wind is on a grid with these lon/lat values.
$ Define data files -------------------------------------------------- $
$ The first input line identifies the file format with FROM, IDLA and
$ IDFM, the second (third) lines give the file unit number and name.
  'NAME' 1 1 '(....)' '(....)'
   40 'ww3_sandy_wind_forc.dat'
This tells WW3 the name of the wind forcing file.
Then cd to WW3/work and run ../exe/ww3_prep to create wind.ww3 (not needed if coupled to wrf).

- edit WW3/work/ww3_shel.inp. Here is where you enter run information.
Here we choose:
   C F     Water levels
   C F     Currents
   C F     Winds
The first letter is to set if the field is to be not used (F), used as read in from a file (T), or used and read in from a coupled model ( C). So those settings are to get coupled water levels, currents, and winds.

You need to set the run times here:
$ Type 1 : Fields of mean wave parameters
   20121028 120000   1800  20121030 120000
and set to have the fields of Hwave, Dwave, and Lwave written out:
$ HS  LM  T02 T0M1 T01 FP DIR SPR DP
  T   T   T   T   T   T   T   T   T

We do NOT need to run ../exe/ww3_shel because that would run the wave model. We will run the wave model as coawstM.

11) Now we are ready to run coawstM.  When the model is done you can visualize the output by using 
cd WW3/work and type
../exe/ww3_ounf
This will create a ww3.*.nc file

