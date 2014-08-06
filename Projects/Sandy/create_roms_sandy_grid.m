% create roms_sandy_grid
%
% There are many ways to make a ROMS grid.
% This is to creatre a basic roms grid, and we want to 
% make sure it fits in part of the WRF grid.
%

% step 1: load the wrf grid
%cd Projects/Sandy
netcdf_load('wrfinput_d01')
figure
pcolorjw(XLONG,XLAT,(1-LANDMASK))
%dasp(35);
hold on

% step 2:  pick 4 corners for roms parent grid
% start in lower left and go clockwise
corner_lon=[ -76.5209 -82.2818  -65.9208  -60.1599 ];
corner_lat=[  27.6300  35.4648   47.3553   39.6587];
plot(corner_lon,corner_lat,'g.','markersize',25)

%step 3: pick 3 of points in the grid
numx=84;
numy=64;

% step4: make a matrix of the lons and lats
% You should not need to change this.
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
plot(lon,lat,'g+')

%step5: make a land mask based on WRF mask
F = TriScatteredInterp(double(XLONG(:)),double(XLAT(:)), ...
                       double(1-LANDMASK(:)),'nearest');
roms_mask=F(lon,lat);
pcolorjw(lon,lat,roms_mask)

%step6: enter name of parent roms grid
roms_grid='Sandy_roms_grid.nc';

% step 7: call to make the grid file
rho.lat=lat;
rho.lon=lon;
rho.depth=zeros(size(rho.lon))+100;
rho.mask=roms_mask;
spherical='T';
%projection='lambert conformal conic';
projection='mercator';

%call generic grid creation
save temp_jcw33.mat
eval(['mat2roms_mw(''temp_jcw33.mat'',''',roms_grid,''');'])
!del temp_jcw33.mat

%User needs to edit roms variables
disp(['     '])
disp(['Created roms grid -->   ',roms_grid])
disp(['     '])
disp(['You need to edit depth in ',roms_grid])
disp(['     '])

%step 8: put bathy into ocean grid.
% you can use a source like this:
% http://www.ngdc.noaa.gov/mgg/global/global.html
% i am using data from a local file
load USeast_bathy.mat
eval(['netcdf_load(''',roms_grid,''');'])
h=griddata(h_lon,h_lat,h_USeast,lon_rho,lat_rho);
%smooth h a little
hnew=h;
hnew(2:end-1,2:end-1)=0.2*(h(1:end-2,2:end-1)+h(2:end-1,2:end-1)+h(3:end,2:end-1)+ ...
                       h(2:end-1,1:end-2)+h(2:end-1,3:end));
figure
pcolorjw(lon_rho,lat_rho,h)
h(isnan(h))=5;
ncwrite(roms_grid,'h',h);

%step 9: create ocean child grid.
% Select child indices
Istr=20; Iend=60; Jstr=35; Jend=62;
%now make some plots and call the coarse-> fine
figure
pcolorjw(lon_rho,lat_rho,mask_rho)
hold on
plot(lon_rho(Istr,Jstr),lat_rho(Istr,Jstr),'g+')
plot(lon_rho(Istr,Jend),lat_rho(Istr,Jend),'g+')
plot(lon_rho(Iend,Jstr),lat_rho(Iend,Jstr),'g+')
plot(lon_rho(Iend,Jend),lat_rho(Iend,Jend),'g+')
ref_ratio=3;
F=coarse2fine('Sandy_roms_grid.nc','Sandy_roms_grid_ref3.nc', ...
              ref_ratio,Istr,Iend,Jstr,Jend);
Gnames={'Sandy_roms_grid.nc','Sandy_roms_grid_ref3.nc'};
[S,G]=contact(Gnames,'Sandy_roms_contact.nc');

%step 10: create roms init conditions, bc's, and nudging files.
% I copied Tools/mfiles/roms_clm/roms_master_climatology_coawst_mw.m
% to Projects/Sandy/roms_master_climatology_sandy.m
% I then edited the roms_master_climatology_sandy.m m file and ran it once
% for the parent grid 'Sandy_roms_grid.nc'. That produced three important files:
% merged_coawst_clm.nc
% merged_coawst_bdy.nc
% coawst.ini.nc
% I renamed those to Sandy_clm.nc, Sandy_bdy.nc, and Sandy_ini.nc.
% Then i ran roms_master_climatology_sandy.m again for the child grid 'Sandy_roms_grid_ref3.nc'
% and renamed the merged files and ini file to create
% Sandy_clm_ref3.nc, Sandy_bdy_ref3.nc, and Sandy_ini_ref3.nc.

% step 11:  create swan grids: run roms2swan twice
roms2swan('Sandy_roms_grid.nc')
% this created swan_coord.grd and swan_bathy.bot and i renamed them to be
% Sandy_swan_bathy.bot and Sandy_swan_coord.grd
roms2swan('Sandy_roms_grid_ref3.nc')
% this created swan_coord.grd and swan_bathy.bot and i renamed them to be
% Sandy_swan_bathy_ref3.bot and Sandy_swan_coord_ref3.grd

%step 12: create SCRIP weights
% use  create_scrip_weights_master.m
% to create all the scrip interpolation weights



