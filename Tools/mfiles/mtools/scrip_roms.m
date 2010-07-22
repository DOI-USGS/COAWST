% script scrip_roms.m (formerly roms2scrip.m)
%
% convert roms grid to scrip input format
%
% jcw Sept 7, 2008
%     Feb 20, 2009
%

%%%%%%%%%% BEGIN user input   %%%%%%%%%%%%
%
%1) Enter name of roms netcdf grid file.
%grid_file = 'joe_tc_coarse_grd.nc';
grid_file = 'Projects\Inlet_test\DiffGrid\inlet_test_grid.nc';

%2) Enter name of netcdf output file to use by scrip.
%out_file = 'joe_tc_coarse_roms_scrip.nc';
out_file = 'Projects\Inlet_test\DiffGrid\inlet_test_roms_scrip.nc';

%%%%%%%%%% END of user input   %%%%%%%%%%%%

%load the roms grid
ncload(grid_file)

[MP,LP]=size(h);
gridsize=LP*MP;

%if ((spherical=='F')||(spherical=='f'))
%  lon_rho=x_rho;
%  lat_rho=y_rho;
%  lon_psi=x_psi;
%  lat_psi=y_psi;
%end

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


disp('step 0/4, filling grid mask')
%get grid mask
count=0;
for jj=1:MP
  for ii=1:LP
    count=count+1;
    ztemp(count)=mask_rho(jj,ii);
  end
end
nc{'grid_imask'}(:) = ztemp(:);
clear ztemp
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

