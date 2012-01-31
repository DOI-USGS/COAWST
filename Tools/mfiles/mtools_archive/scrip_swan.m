% script scrip_swan.m
%
% convert swan grid to scrip input format
%
% jcw July 20, 2010
%

%%%%%%%%%% BEGIN user input   %%%%%%%%%%%%
%
%1) Enter name of swan .grd ascii grid file.
grid_file = 'Projects\Inlet_test\DiffGrid\inlet_test_grid_coord2.grd';

%2) Enter the grid size. This will be +1 from the values entered for the 
%   SWAN input file. For example, if you use:
%   CGRID CURVILINEAR 86 81 EXC 9.999000e+003 & ...
%   THEN enter Numx = 87;  Numy = 82;
Numx=87;
Numy=82;

%3) Are the grid coordinates in lat/lon or in meters?
%   Enter 0 for lat/lon,  1 for meters.
%   If no, then they should be 
grid_coords=1;  % this was in meters

%3) Enter name of the bathy file. This will be used to compute the
%   grid mask.
bath_file = 'D:\data\models\COAWST\Projects\Inlet_test\DiffGrid\inlet_test_bathy2.bot';

%4) Enter name of netcdf output file to use by scrip.
out_file = 'Projects\Inlet_test\DiffGrid\inlet_test_swan_scrip.nc';

%%%%%%%%%% END of user input   %%%%%%%%%%%%

%load the swan grid and bathy
bath=load(bath_file);
bath=reshape(bath,Numx,Numy)';
mask=ones(size(bath));
mask(bath==9999)=0;

grd=load(grid_file);
gridsize=length(grd)/2;
lon_rho=grd(1:gridsize);
lat_rho=grd(gridsize+1:end);

%If user has grid locations in meters, then convert to lat/lon.
if (grid_coords)
  projection='mercator';
  m_proj(projection);
  [lon_rho, lat_rho] = m_xy2ll(lon_rho/6371000, lat_rho/6371000);   % Degrees.
end

lon_rho=reshape(lon_rho,Numx,Numy)';
lat_rho=reshape(lat_rho,Numx,Numy)';

XLONGI=interp2(lon_rho,'bilinear');
XLATI=interp2(lat_rho,'bilinear');
lon_psi=XLONGI(2:2:end,2:2:end);
lat_psi=XLATI(2:2:end,2:2:end);

[MP,LP]=size(lon_rho);
gridsize=LP*MP;

%create a full set of psi points
[x_full_grid,y_full_grid]=create_extra_rho_grid(lon_psi,lat_psi);

%create the srcip netcdf grid file

nc = netcdf(out_file, 'clobber');

%% Global attributes:
disp(' ## Defining Global Attributes...')
 
nc.title = ncchar('Scrip file');
nc.history = ncchar(['Created by "' mfilename '" on ' datestr(now)]);

L = LP-1;       % The psi dimension.
M = MP-1;       % The psi dimension.

disp(' ## Defining Dimensions...')
 
nc('grid_size') = gridsize;
nc('grid_corners') = 4;
nc('grid_rank') = 2;

disp(' ## Defining Variables')
 
nc{'grid_dims'} = ncshort('grid_rank');

nc{'grid_imask'} = ncshort('grid_size');
nc{'grid_imask'}.units = ncchar('---');

nc{'grid_center_lat'} = ncdouble('grid_size');
nc{'grid_center_lat'}.units = ncchar('radians');

nc{'grid_center_lon'} = ncdouble('grid_size');
nc{'grid_center_lon'}.units = ncchar('radians');

nc{'grid_corner_lat'} = ncdouble('grid_size', 'grid_corners');
nc{'grid_corner_lat'}.units = ncchar('radians');

nc{'grid_corner_lon'} = ncdouble('grid_size', 'grid_corners');
nc{'grid_corner_lon'}.units = ncchar('radians');

ncclose

%now fill that netcdf file
eval(['nc=netcdf(''',out_file,''',''w'');'])

nc{'grid_dims'}(:) = [MP, LP];
scale=pi/180;

%get grid mask
ztemp=reshape(mask',1,MP*LP);
nc{'grid_imask'}(:) = ztemp;

%get grid centers
disp('step 1/4, filling grid lat centers')
ztemp=reshape(lat_rho',1,MP*LP);
nc{'grid_center_lat'}(:) = ztemp(:)*scale;
clear ztemp

disp('step 2/4, filling grid lon centers')
ztemp=reshape(lon_rho',1,MP*LP);
nc{'grid_center_lon'}(:) = ztemp(:)*scale;
clear ztemp

%get grid corners, counterclockwise
disp('step 3/4, filling grid lat corners')
c1=reshape(y_full_grid(1:MP,1:LP)',MP*LP,1);
c2=reshape(y_full_grid(1:MP,2:LP+1)',MP*LP,1);
c3=reshape(y_full_grid(2:MP+1,2:LP+1)',MP*LP,1);
c4=reshape(y_full_grid(2:MP+1,1:LP)',MP*LP,1);
ztemp(:,:)=[c1 c2 c3 c4];
nc{'grid_corner_lat'}(:) = ztemp(:)*scale;
clear ztemp

disp('step 4/4, filling grid lon corners')
c1=reshape(x_full_grid(1:MP,1:LP)',MP*LP,1);
c2=reshape(x_full_grid(1:MP,2:LP+1)',MP*LP,1);
c3=reshape(x_full_grid(2:MP+1,2:LP+1)',MP*LP,1);
c4=reshape(x_full_grid(2:MP+1,1:LP)',MP*LP,1);
ztemp(:,:)=[c1 c2 c3 c4];
nc{'grid_corner_lon'}(:) = ztemp(:)*scale;
clear ztemp

%close the file.
if ~isempty(close(nc))
	disp(' ## Unable too close the ROMS output file.')
end
ncclose('nc')

