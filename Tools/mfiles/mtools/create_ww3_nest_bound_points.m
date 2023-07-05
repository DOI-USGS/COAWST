clear
cd D:\NOPP_hurricanes_modeling\grids

roms_grid='L2_MEXI_GOMSAB_2km_CUDEM_v3.nc';
%roms_grid='useast_grd5_2_cnapsv2.nc';
%roms_grid='GOMSAB_2km_ext_smooth.nc';
netcdf_load(roms_grid)
[LP, MP]=size(h);
%
lon=[lon_rho(:,1)' lon_rho(end,:) lon_rho(end:-1:1,end)' lon_rho(1,end:-1:1)];
lat=[lat_rho(:,1)' lat_rho(end,:) lat_rho(end:-1:1,end)' lat_rho(1,end:-1:1)];
%
fid = fopen('bc_points_ww3gridinp.dat','w');
for i = 1:length(lon)
  fprintf(fid,'%10.6f %10.6f %3.2f %3.2f %2i\n',lon(i),lat(i),0.0,0.0,1);
end
fclose(fid);
%

