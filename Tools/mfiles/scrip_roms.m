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
grid_file = '../../Projects/JOE_TC/joe_tc_coarse_grd.nc';

%2) Enter name of netcdf output file to use by scrip.
out_file = '../../Projects/JOE_TC/joe_tc_coarse_roms_scrip.nc';

%%%%%%%%%% END of user input   %%%%%%%%%%%%

%load the roms grid
ncload(grid_file)

[MP,LP]=size(h);
gridsize=LP*MP;

if ((spherical=='F')||(spherical=='f'))
  lon_rho=x_rho;
  lat_rho=y_rho;
end

%create a full set of psi points
[x_full_grid,y_full_grid]=create_extra_rho_grid(lon_rho,lat_rho);

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
nc{'grid_center_lat'}.units = ncchar('meters');

nc{'grid_center_lon'} = ncdouble('grid_size');
nc{'grid_center_lon'}.units = ncchar('meters');

nc{'grid_corner_lat'} = ncdouble('grid_size', 'grid_corners');
nc{'grid_corner_lat'}.units = ncchar('meters');

nc{'grid_corner_lon'} = ncdouble('grid_size', 'grid_corners');
nc{'grid_corner_lon'}.units = ncchar('meters');

ncclose

%now fill that netcdf file
eval(['nc=netcdf(''',out_file,''',''w'');'])

nc{'grid_dims'}(:) = [MP, LP];
nc{'grid_imask'}(:) = ones(1,MP*LP);

disp('step 1/4, filling grid lat centers')
%get grid centers
count=0;
for jj=1:MP
  for ii=1:LP
    count=count+1;
    ztemp(count)=lat_rho(jj,ii);
  end
end
nc{'grid_center_lat'}(:) = ztemp(:)*pi/180;
clear ztemp

disp('step 2/4, filling grid lon centers')
count=0;
for jj=1:MP
  for ii=1:LP
    count=count+1;
    ztemp(count)=lon_rho(jj,ii);
  end
end
nc{'grid_center_lon'}(:) = ztemp(:)*pi/180;
clear ztemp

disp('step 3/4, filling grid lat corners')
%get grid corners, counterclockwise
count=0;
for jj=2:MP+1
  for ii=2:LP+1
    count=count+1;
    ztemp(count,:)=[y_full_grid(jj-1,ii) y_full_grid(jj,ii+1) ...
                    y_full_grid(jj+1,ii) y_full_grid(jj,ii-1)];
  end
end
nc{'grid_corner_lat'}(:) = ztemp(:)*pi/180;

disp('step 4/4, filling grid lon corners')
count=0;
for jj=2:MP+1
  for ii=2:LP+1
    count=count+1;
    ztemp(count,:)=[x_full_grid(jj-1,ii) x_full_grid(jj,ii+1) ...
                    x_full_grid(jj+1,ii) x_full_grid(jj,ii-1)];
  end
end
nc{'grid_corner_lon'}(:) = ztemp(:)*pi/180;

%close the file.
if ~isempty(close(nc))
	disp(' ## Unable too close the ROMS output file.')
end
ncclose('nc')

