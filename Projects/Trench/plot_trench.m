%plot_delilah

cd D:\models\COAWST_updates\COAWST_v3.9\COAWST_v3.9_test18_all_PROJECTS_RUN
netcdf_load('Output_trench/ocean_trench_his.nc')

figure
  plot(x_rho(:,4),-squeeze(bath(:,4,1)),'k')
  hold on
  plot(x_rho(:,4),-squeeze(bath(:,4,10)),'r')
  title('init and final along channel bathymetry')
  legend ('Initial','Final')
  xlabel('along channel distance, m')
  ylabel('bathymetry, m')
  
  print -dpng 'D:\models\COAWST\Projects\Trench\trench.png'
