% script scrip_wrf.m (formerly wrf2scrip.m)
%
% convert wrf grid to scrip netcdf input grid.
%
% jcw 14Sep2008, 20Feb2009
%

%%%%%%%%%% BEGIN user input   %%%%%%%%%%%%
%
%1) Enter name of file with WRF grid info.
grid_file = '../../Projects/JOE_TC/wrfinput_d01';

%2) Enter name of netcdf output file to use by scrip.
out_file = '../../Projects/JOE_TC/joe_tc_wrf_scrip.nc';

%%%%%%%%%% END user input   %%%%%%%%%%%%

% get wrf grid stuff
nc=netcdf(grid_file);
lon_rho=nc{'XLONG'}(:);
lat_rho=nc{'XLAT'}(:);




LP=size(XLONG_U,2);
MP=size(XLAT_V,1);

gridsize=LP*MP;

lon_rho_full=repmat(XLONG_U(1,:),MP,1);
lat_rho_full=repmat(XLAT_V(:,1),1,LP);

%create a full set of psi points
[x_full_grid,y_full_grid]=create_extra_rho_grid(lon_rho_full,lat_rho_full);

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

%get grid centers
count=0;
for jj=1:MP
  for ii=1:LP
    count=count+1;
    ztemp(count)=lat_rho_full(jj,ii);
  end
end
nc{'grid_center_lat'}(:) = ztemp(:)*pi/180;
clear ztemp

count=0;
for jj=1:MP
  for ii=1:LP
    count=count+1;
    ztemp(count)=lon_rho_full(jj,ii);
  end
end
nc{'grid_center_lon'}(:) = ztemp(:)*pi/180;
clear ztemp

%get grid corners, counterclockwise
disp('getting corners')

count=0;
for jj=2:MP+1
  for ii=2:LP+1
    count=count+1;
    ztemp(count,:)=[y_full_grid(jj-1,ii) y_full_grid(jj,ii+1) ...
                    y_full_grid(jj+1,ii) y_full_grid(jj,ii-1)];
  end
end
nc{'grid_corner_lat'}(:) = ztemp(:)*pi/180;

count=0;
for jj=2:MP+1
  for ii=2:LP+1
    count=count+1;
%   ztemp(count,:)=[x_full_grid(jj,ii)     x_full_grid(jj,ii+1) ...
%                   x_full_grid(jj+1,ii+1) x_full_grid(jj+1,ii)];
    ztemp(count,:)=[x_full_grid(jj-1,ii) x_full_grid(jj,ii+1) ...
                    x_full_grid(jj+1,ii) x_full_grid(jj,ii-1)];
  end
end
nc{'grid_corner_lon'}(:) = ztemp(:)*pi/180;

%close the file.
if ~isempty(close(nc))
	disp(' ## Unable to close the netcdf file.')
end
ncclose('nc')

