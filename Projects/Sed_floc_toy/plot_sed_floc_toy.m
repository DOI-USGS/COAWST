%plot_delilah

cd D:\models\COAWST_updates\COAWST_v3.9\COAWST_v3.9_test18_all_PROJECTS_RUN
netcdf_load('Output_sed_floc_toy/ocean_sed_floc_toy_his.nc')

figure
subplot(121)
  hold on
  for mm=1:15
    toadd=['00',num2str(mm)];toadd=toadd(end-1:end);
    eval(['plot(squeeze(mud_',toadd,'(2,2,:,73)))'])
  end
  title('final vertical mud profiles')
  legend
subplot(122)
  hold on
  for mm=1:4
    toadd=['00',num2str(mm)];toadd=toadd(end-1:end);
    eval(['plot(squeeze(sand_',toadd,'(2,2,:,73)))'])
  end
  legend
  title('final vertical sand profiles')

  print -dpng 'D:\models\COAWST\Projects\Sed_floc_toy\mud_sand.png'
