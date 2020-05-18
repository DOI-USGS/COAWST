% to create ww3 wind forcing file.
%
% jcwarner 13Nov2017
%

% This is very hacky for now. In future we will provide 
% a file that directly reads NCEP data and writes out a
% WW3 wind forcing ascii file.
%
netcdf_load('../romsforc_NARR_Sandy2012.nc')
[LP, MP, ntimes]=size(Uwind);
%
fid = fopen('ww3_sandy_wind_forc.dat','w');
for mm=1:ntimes
% write time stamp
  zz=datevec(wind_time(mm)+datenum(1858,11,17,0,0,0));
  zz1=num2str(zz(1));
  zz2=['00',num2str(zz(2))];zz2=zz2(end-1:end);
  zz3=['00',num2str(zz(3))];zz3=zz3(end-1:end);
  zz4=['00',num2str(zz(4))];zz4=zz4(end-1:end);
  zz5=['00',num2str(zz(5))];zz5=zz5(end-1:end);
  zz6=['00',num2str(zz(6))];zz6=zz6(end-1:end);
  dtstr=[zz1,zz2,zz3,' ',zz4,zz5,zz6]
  fprintf(fid,'%15s',dtstr);
  fprintf(fid,'\n');
  zzu=squeeze(Uwind(:,:,mm));
  zzv=squeeze(Vwind(:,:,mm));
  for index1 = 1:MP;
    fprintf(fid,'%12.4f',zzu(:,index1));
    fprintf(fid,'\n');
  end
  for index1 = 1:MP;
    fprintf(fid,'%12.4f',zzv(:,index1));
    fprintf(fid,'\n');
  end
end
fclose(fid);



