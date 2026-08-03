% wetdry_dam_break
%
% plot time series of water level from
% dambreak test case.
%
cd D:\models\COAWST_updates\COAWST_v3.9\COAWST_v3.9_test18_all_PROJECTS_RUN
netcdf_load('Output_wetdry/wetdry_dam_break_his.nc')
% make 6 panel plan view plot - surf
tidx=[2 3 4 11 21 41];
frame=0;
figure
for mm=1:6
  colormap('gray')
  eval(['subplot(2,3,mm)'])
  zz=surf(x_rho,y_rho,double(zeta(:,:,tidx(mm))));
  caxis([-0.2 0.7])
  colorbar
  title(['time = ',num2str(ocean_time(tidx(mm))),' sec'])
  set(gca,'CLim',[-0.2 0.7], ...
          'View',[157.5 56], ...
          'XLim',[0 4], ...
          'YLim',[0 2], ...
          'ZLim',[0 0.6])
  set(gca,'ztick',[0:0.2:0.6])
end

print -dpng 'D:\models\COAWST\Projects\Wetdry\wetdry_dambreak.png'
