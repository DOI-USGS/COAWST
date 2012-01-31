% script scrip_wrf.m (formerly wrf2scrip.m)
%
% convert wrf grid to scrip netcdf input grid.
%
% jcw 14Sep2008, 20Feb2009
%

%%%%%%%%%% BEGIN user input   %%%%%%%%%%%%
%
%1) Enter name of file with WRF grid info.
grid_file = 'wrfinput_d01';

%2) Enter name of netcdf output file to use by scrip.
out_file = 'joe_tc_wrf_scrip.nc';

%%%%%%%%%% END user input   %%%%%%%%%%%%

% get wrf grid stuff
nc=netcdf(grid_file);
lon_rho=nc{'XLONG'}(:);
lat_rho=nc{'XLAT'}(:);

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
nc{'grid_imask'}(:) = ones(1,MP*LP);
scale=pi/180;

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

