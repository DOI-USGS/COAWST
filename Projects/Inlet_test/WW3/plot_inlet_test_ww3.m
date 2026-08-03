%plot_delilah

cd D:\models\COAWST_updates\COAWST_v3.9\COAWST_v3.9_test18_all_PROJECTS_RUN
netcdf_load('Output_inlet_test_ww3/ocean_his_coupled.nc')

figure
subplot(121)
  pcolorjw(x_rho,y_rho,squeeze(Hwave(:,:,13)))
  colorbar
  xlabel('x-rho, m')
  ylabel('y-rho, m')
  caxis([0 1.2])
  title('Hsig at 12 hours, m')
subplot(122)
  pcolorjw(x_rho,y_rho,squeeze(vWave(:,:,13)))
  colorbar
  xlabel('x-rho, m')
  ylabel('y-rho, m')
  caxis([-0.1 1.6])
  title('vWave sent to wave model, freq bin13, at 12 hours, m/s')
print -dpng 'D:\models\COAWST\Projects\Inlet_test\WW3\hwave_vwave.png'
