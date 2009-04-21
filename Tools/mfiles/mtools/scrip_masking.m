% script scrip_masking.m
%
% Adjusts the masking in the wrf scrip grid.
% Determines what cells are not in the roms grid.
%
% jcw Sept 7, 2008
%     Feb 20, 2009
%

%%%%%%%%%% BEGIN user input   %%%%%%%%%%%%
%
%1) Enter name of roms netcdf grid file.
grid_filer = 'joe_tc_coarse_grd.nc';

%2) Enter name of file with WRF grid info.
grid_filew = 'wrfinput_d01';

%%%%%%%%%% END of user input   %%%%%%%%%%%%

%get roms grid
ncload(grid_filer)
if ((spherical=='F')||(spherical=='f'))
  lon_rho=x_rho;
  lat_rho=y_rho;
end
cornersx=[lon_rho(1,1) lon_rho(1,end) lon_rho(end,end) lon_rho(end,1)];
cornersy=[lat_rho(1,1) lat_rho(1,end) lat_rho(end,end) lat_rho(end,1)];

%get wrf grid
nc=netcdf(grid_filew);
lon_rho_wrf=nc{'XLONG'}(:);
lat_rho_wrf=nc{'XLAT'}(:);
[MP,LP]=size(lon_rho_wrf);
ncclose('nc')

%find wrf points inside roms grid
wrfinroms=inpolygon(lon_rho_wrf,lat_rho_wrf,cornersx,cornersy);

%now put this into the mask array
ncload ocn2atm_weights.nc
nc=netcdf('ocn2atm_weights.nc','w');
count=0;
ztemp=zeros(LP*MP,1);
for jj=1:MP
  for ii=1:LP
    count=count+1;
    ztemp(count)=wrfinroms(jj,ii);
  end
end
nc{'dst_grid_imask'}(:) = ztemp(:);
clear ztemp
ncclose('nc')
