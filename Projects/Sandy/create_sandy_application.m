% create_sandy_application
%
% jcwarner 19August2014
% a list of procedures and notes
% used to create the Sandy app files
% for roms and swan. The wrf files
% are assumed to already be created.
%
%
% cd Projects/Sandy
% This is a place to create all the necessary files.
%
% Step 1:  ROMS grid
%
% As a first cut, you can use the wrf grid. 
% Load the wrf grid to get coastline data
netcdf_load('wrfinput_d01')
figure
pcolorjw(XLONG,XLAT,double(1-LANDMASK))
hold on
%
% pick 4 corners for roms parent grid
Istr=27; Iend=84;
Jstr=7 ; Jend=64;
plot(XLONG(Istr,Jstr),XLAT(Istr,Jstr),'g+')
plot(XLONG(Iend,Jstr),XLAT(Iend,Jstr),'g+')
plot(XLONG(Istr,Jend),XLAT(Istr,Jend),'g+')
plot(XLONG(Iend,Jend),XLAT(Iend,Jend),'g+')

%step 3: pick number of points in the grid
numx=84;
numy=64;

% Make a matrix of the lons and lats
% You should not need to change this.
% start in lower left and go clockwise
corner_lon=double([XLONG(Istr,Jstr) XLONG(Istr,Jend) XLONG(Iend,Jend) XLONG(Iend,Jstr) ]);
corner_lat=double([XLAT(Istr,Jstr) XLAT(Istr,Jend) XLAT(Iend,Jend) XLAT(Iend,Jstr) ]);
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
%
% Call generic grid creation
%
roms_grid='Sandy_roms_grid.nc';
rho.lat=lat;
rho.lon=lon;
rho.depth=zeros(size(rho.lon))+100;  % for now just make zeros
rho.mask=zeros(size(rho.lon));       % for now just make zeros
spherical='T';
%projection='lambert conformal conic';
projection='mercator';
save temp_jcw33.mat rho spherical projection
eval(['mat2roms_mw(''temp_jcw33.mat'',''',roms_grid,''');'])
!del temp_jcw33.mat
%User needs to edit roms variables
disp(['     '])
disp(['Created roms grid -->   ',roms_grid])
disp(['     '])
disp(['You need to edit depth in ',roms_grid])
disp(['     '])
%
% Step 2: Masking
%
% base this on the WRF mask.
F = TriScatteredInterp(double(XLONG(:)),double(XLAT(:)), ...
                       double(1-LANDMASK(:)),'nearest');
roms_mask=F(lon,lat);
figure
pcolorjw(lon,lat,roms_mask)
%
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
% http://www.ngdc.noaa.gov/mgg/coast/
% select lat bounds of 45 / 25 lon of -84 / -60 
%
lon=coastline(:,1);
lat=coastline(:,2);
save coastline.mat lon lat
%
% use editmask to create better mask
%
editmask
%
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
hnew=h;
hnew(2:end-1,2:end-1)=0.2*(h(1:end-2,2:end-1)+h(2:end-1,2:end-1)+h(3:end,2:end-1)+ ...
                       h(2:end-1,1:end-2)+h(2:end-1,3:end));
% Bathymetry can be smoothed using
% http://drobilica.irb.hr/~mathieu/Bathymetry/index.html 
figure
pcolorjw(lon_rho,lat_rho,h)
hnew(isnan(hnew))=5;
ncwrite(roms_grid,'h',hnew);
%
% Step 4: child grid
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
hnew=h;
hnew(2:end-1,2:end-1)=0.2*(h(1:end-2,2:end-1)+h(2:end-1,2:end-1)+h(3:end,2:end-1)+ ...
                       h(2:end-1,1:end-2)+h(2:end-1,3:end));
% Bathymetry can be smoothed using
% http://drobilica.irb.hr/~mathieu/Bathymetry/index.html 
figure
pcolorjw(lon_rho,lat_rho,hnew)
ncwrite('Sandy_roms_grid_ref3.nc','h',hnew);
%
% recompute child mask based on WRF mask
%
netcdf_load('wrfinput_d02');
F = TriScatteredInterp(double(XLONG(:)),double(XLAT(:)), ...
                       double(1-LANDMASK(:)),'nearest');
roms_mask=F(lon_rho,lat_rho);
figure
pcolorjw(lon_rho,lat_rho,roms_mask)
water = double(roms_mask);
u_mask = water(1:end-1,:) & water(2:end,:);
v_mask= water(:,1:end-1) & water(:,2:end);
psi_mask= water(1:end-1,1:end-1) & water(1:end-1,2:end) & water(2:end,1:end-1) & water(2:end,2:end);
ncwrite('Sandy_roms_grid_ref3.nc','mask_rho',roms_mask);
ncwrite('Sandy_roms_grid_ref3.nc','mask_u',double(u_mask));
ncwrite('Sandy_roms_grid_ref3.nc','mask_v',double(v_mask));
ncwrite('Sandy_roms_grid_ref3.nc','mask_psi',double(psi_mask));
%
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
% Step 6: 2D BC's
%
% Get the tidal data at
svn checkout https://coawstmodel.sourcerepo.com/coawstmodel/data .
% edit Tools/mfiles/tides/create_roms_tides to select tidal database
% and start time. then run that file
create_roms_tides
%
% Step 7: Surface forcings
%
% To create a forcing file with just wind forcing
% you need to get wind data, such as from  http://nomads.ncdc.noaa.gov/
% I got uwnd.10m.2012.nc and vwnd.10m.2012.nc
% then run the m file:
narr2romsnc
% this also creates a SWAN wind forcing file.
% I did this rwice: once for the parent and once for the child.
%
% To create a forcing file that has winds and many other vars like
% Uwind, Vwind, Pair, Tair, Qair, rain, swrad, lwrad you can use:
% COAWST/Tools/mfiles/mtools/create_roms_forcings 
%
% Step 8: roms input files.
%
edit sandy.h
%and
edit ocean_sandy.in
%
% Step 9: build file
%
edit and run coawst.bash 
%
% Step 10: run file
%
edit run_nemo


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  coupling to WRF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Step 11: create SCRIP weights
create_scrip_weights_master_sandy
% to create all the scrip interpolation weights
%
% Step 12: input files.
%
edit sandy.h
%and
edit ocean_sandy.in
%and
edit coupling_sandy.in
%
% Step 13: build file
%
edit and run coawst.bash 
%
% Step 14: run file
%
copy wrf* and namelist.input from Projects/Sandy to the root directory.
edit run_nemo

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  coupling to SWAN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Step 15:  create swan grids: run roms2swan twice
roms2swan('Sandy_roms_grid.nc')
% this created swan_coord.grd and swan_bathy.bot and i renamed them to be
% Sandy_swan_bathy.bot and Sandy_swan_coord.grd
roms2swan('Sandy_roms_grid_ref3.nc')
% this created swan_coord.grd and swan_bathy.bot and i renamed them to be
% Sandy_swan_bathy_ref3.bot and Sandy_swan_coord_ref3.grd
%
% Step 16 - create SWAN boundary TPAR files
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
% then run the m file to create the TPAR files
COAWST\Tools\mfiles\swan_forc\ww3_swan_input.m
% - this file prodcues TPAR files. Place the TPAR files in 
%  your project folder. 
% - Also, the m files creates some commands lines called
% "BOUNDSPEC." You must copy the command lines from the file 
% Bound_spec_command and place them into your SWAN INPUT file.
% This is only done for the parent grid.
%
% Step 17 - create SWAN init files.
% To create an init file for swan, you can run SWAN in stationary mode
% and create a ‘hot start’ file(s). 
%
