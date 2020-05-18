% plot_scrip_weights
%
% jcw 04Jan2018
%

%1) enter scrip weight file name
%fname='scrip_sandy_nowavenest.nc';
fname='scrip_florence_1romshydro.nc';

%2) enter number of wrf,roms,ww3,and swan grids.
NGRIDS_ROMS=1;
NGRIDS_SWAN=0;
NGRIDS_WW3=0;
NGRIDS_WRF=0;
NGRIDS_HYDRO=1;

%%%%%%%%%%%%%%%%%  END OF USER INPUT  %%%%%%%%%%%%%%%%%%%

count=0;
clear strnames;
for mo=1:NGRIDS_ROMS
  for mw=1:NGRIDS_SWAN
    count=count+1;
    strnames(count,:)=['ocn',num2str(mo),'_to_wav',num2str(mw)];
  end
  for mw=1:NGRIDS_WW3
    count=count+1;
    strnames(count,:)=['ocn',num2str(mo),'_to_wav',num2str(mw)];
  end
  for ma=1:NGRIDS_WRF
    count=count+1;
    strnames(count,:)=['ocn',num2str(mo),'_to_atm',num2str(ma)];
  end
  for mh=1:NGRIDS_HYDRO
    count=count+1;
    strnames(count,:)=['ocn',num2str(mo),'_to_hyd',num2str(mh)];
  end
end
for mw=1:NGRIDS_SWAN
  for mo=1:NGRIDS_ROMS
    count=count+1;
    strnames(count,:)=['wav',num2str(mo),'_to_ocn',num2str(mw)];
  end
  for ma=1:NGRIDS_WRF
    count=count+1;
    strnames(count,:)=['wav',num2str(mo),'_to_atm',num2str(ma)];
  end
end
for mw=1:NGRIDS_WW3
  for mo=1:NGRIDS_ROMS
    count=count+1;
    strnames(count,:)=['wav',num2str(mw),'_to_ocn',num2str(mo)];
  end
  for ma=1:NGRIDS_WRF
    count=count+1;
    strnames(count,:)=['wav',num2str(mw),'_to_atm',num2str(ma)];
  end
end
for ma=1:NGRIDS_WRF
  for mw=1:NGRIDS_SWAN
    count=count+1;
    strnames(count,:)=['atm',num2str(ma),'_to_wav',num2str(mw)];
  end
  for mw=1:NGRIDS_WW3
    count=count+1;
    strnames(count,:)=['atm',num2str(ma),'_to_wav',num2str(mw)];
  end
  for mo=1:NGRIDS_ROMS
    count=count+1;
    strnames(count,:)=['atm',num2str(ma),'_to_ocn',num2str(mo)];
  end
end
for mh=1:NGRIDS_HYDRO
  for mo=1:NGRIDS_ROMS
    count=count+1;
    strnames(count,:)=['hyd',num2str(mh),'_to_ocn',num2str(mo)];
  end
end

for mm=1:count
  str=strnames(mm,:);
  remap=ncread(fname,['/',str,'_weights.nc/remap_matrix']);
  src=ncread(fname,['/',str,'_weights.nc/src_address']);
  dst=ncread(fname,['/',str,'_weights.nc/dst_address']);
  dst_mask=ncread(fname,['/',str,'_weights.nc/dst_grid_imask']);
  dst_lat=ncread(fname,['/',str,'_weights.nc/dst_grid_center_lat']);
  dst_lon=ncread(fname,['/',str,'_weights.nc/dst_grid_center_lon']);
  dst_cor_lat=ncread(fname,['/',str,'_weights.nc/dst_grid_corner_lat']);
  dst_cor_lon=ncread(fname,['/',str,'_weights.nc/dst_grid_corner_lon']);
  src_lon=ncread(fname,['/',str,'_weights.nc/src_grid_center_lon']);
  src_lat=ncread(fname,['/',str,'_weights.nc/src_grid_center_lat']);
  src_cor_lat=ncread(fname,['/',str,'_weights.nc/src_grid_corner_lat']);
  src_cor_lon=ncread(fname,['/',str,'_weights.nc/src_grid_corner_lon']);
  src_size=ncread(fname,['/',str,'_weights.nc/src_grid_dims']);
  dst_size=ncread(fname,['/',str,'_weights.nc/dst_grid_dims']);
%
  zd=double([1:dst_size(1)*dst_size(2)]*0.);
  zs=double([1:src_size(1)*src_size(2)]*0.);
  for mm=1:length(dst)
    zd(dst(mm))=zd(dst(mm))+remap(1,mm);
    zs(src(mm))=zs(src(mm))+1;
  end
%
  figure
  subplot(211)
  zs2=reshape(zs,src_size(1),src_size(2));
  src_lon(src_lon>pi)=src_lon(src_lon>pi)-2*pi;
  src_lon=reshape(src_lon,src_size(1),src_size(2))*180/pi;
  src_lat=reshape(src_lat,src_size(1),src_size(2))*180/pi;
  pcolorjw(src_lon,src_lat,zs2); colorbar
  title([str(1:4),' to ',str(9:12),' number of times these src location is used'])
  subplot(212)
  zd2=reshape(zd,dst_size(1),dst_size(2));
  dst_lon(dst_lon>pi)=dst_lon(dst_lon>pi)-2*pi;
  dst_lon=reshape(dst_lon,dst_size(1),dst_size(2))*180/pi;
  dst_lat=reshape(dst_lat,dst_size(1),dst_size(2))*180/pi;
  dst_mask=double(reshape(dst_mask,dst_size(1),dst_size(2)));
  pcolorjw(dst_lon,dst_lat,zd2.*dst_mask); colorbar
  title([str(1:4),' to ',str(9:12),' sum of weigths to each location on dst grid'])
  hold on
  %src_cor_lon(src_cor_lon>pi)=src_cor_lon(src_cor_lon>pi)-2*pi;
  %plot(src_cor_lon(:)*180/pi+360,src_cor_lat(:)*180/pi,'r+')
end

% compute sum of weights on destination grids
count=0;
clear zstring;
for nn=1:NGRIDS_ROMS
  count=count+1;
  zstring(count,:)=['to_ocn',num2str(nn)];
end
for nn=1:NGRIDS_SWAN
  count=count+1;
  zstring(count,:)=['to_wav',num2str(nn)];
end
for nn=1:NGRIDS_WW3
  count=count+1;
  zstring(count,:)=['to_wav',num2str(nn)];
end
for nn=1:NGRIDS_WRF
  count=count+1;
  zstring(count,:)=['to_atm',num2str(nn)];
end
for nn=1:NGRIDS_HYDRO
  count=count+1;
  zstring(count,:)=['to_hyd',num2str(nn)];
end
%
% Plot sum of destination weights to each grid
%
for aa=1:count
  start=1;
  for ii=1:size(strnames,1)
    if strmatch(zstring(aa,:),strnames(ii,end-6:end))
      str=strnames(ii,:)
      remap=ncread(fname,['/',str,'_weights.nc/remap_matrix']);
      dst=ncread(fname,['/',str,'_weights.nc/dst_address']);
      dst_mask=double(ncread(fname,['/',str,'_weights.nc/dst_grid_imask']));
      dst_lat=ncread(fname,['/',str,'_weights.nc/dst_grid_center_lat']);
      dst_lon=ncread(fname,['/',str,'_weights.nc/dst_grid_center_lon']);
      dst_size=ncread(fname,['/',str,'_weights.nc/dst_grid_dims']);
      if start==1
        zd=double([1:dst_size(1)*dst_size(2)]*0.);
        start=0
      end
      for mm=1:length(dst)
        zd(dst(mm))=zd(dst(mm))+remap(1,mm)*dst_mask(dst(mm));
      end
    end
  end

  figure
  zd2=reshape(zd,dst_size(1),dst_size(2));
  dst_lon(dst_lon>pi)=dst_lon(dst_lon>pi)-2*pi;
  dst_lon=reshape(dst_lon,dst_size(1),dst_size(2))*180/pi;
  dst_lat=reshape(dst_lat,dst_size(1),dst_size(2))*180/pi;
  dst_mask=double(reshape(dst_mask,dst_size(1),dst_size(2)));
  pcolorjw(dst_lon,dst_lat,zd2); colorbar
  title(['Sum of all weights to dest',zstring(aa,4:end),' should be 0 to 1'])

end

