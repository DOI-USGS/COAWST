%plot_delilah

cd D:\models\COAWST_updates\COAWST_v3.9\COAWST_v3.9_test18_all_PROJECTS_RUN
netcdf_load('Projects/estuary_test2/estuary_test2_grid.nc')
netcdf_load('Output_estuary_test2/ocean_estuary2_his.nc')

figure
pcolorjw(lon_rho,lat_rho,squeeze(h(:,:)))
colorbar
xlabel('Longitude')
ylabel('Latitude')
title('Depth, m')
print -dpng 'D:\models\COAWST\Projects\Estuary_test2\depth.png'

figure
N=size(salt,3);
ot=60; %length(ocean_time)
igrid=1 %rho
z_rho=set_depth(Vtransform, Vstretching, ...
                       theta_s, theta_b, hc, N, ...
                       igrid, h, zeta(:,:,ot));
pcolorjw(repmat(lon_rho(:,31),1,N),squeeze(z_rho(:,31,:)),squeeze(salt(:,31,:,ot)))
zoom on
colorbar
xlabel('Longitude')
ylabel('Depth, m')
title('Along channel salinity at time 60')
print -dpng 'D:\models\COAWST\Projects\Estuary_test2\salt.png'

