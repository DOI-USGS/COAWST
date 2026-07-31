%plot_delilah

cd D:\models\COAWST_updates\COAWST_v3.9\COAWST_v3.9_test18_all_PROJECTS_RUN
netcdf_load('Output_delilah/ocean_delilah_his.nc')
netcdf_load('Projects/Delilah/InWave_delilah_grd5_per.nc')
figure
pcolorjw(x_rho,y_rho,squeeze(Hwave(:,:,7)))
colorbar
caxis([0 2.5])
xlabel('Cross-shore distance, m')
ylabel('Alongshore distance, m')
title('Wave Heights at time 328.00, m')

print -dpng 'D:\models\COAWST\Projects\Delilah\delilah_waves.png'