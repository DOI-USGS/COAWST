% plot_scrip_weights
%
% jcw 04Jan2018
%

%1) enter scrip weight file name
fname='scrip_sandy_nowavenest.nc';

%2) enter number of wrf,roms,ww3,and swan grids.
NGRIDS_ROMS=1;
NGRIDS_SWAN=0;
NGRIDS_WW3=1;
NGRIDS_WRF=1;

%%%%%%%%%%%%%%%%%  END OF USER INPUT  %%%%%%%%%%%%%%%%%%%

%num_str=(NGRIDS_ROMS+NGRIDS_SWAN+NGRIDS_WW3+NGRIDS_WRF);
%for mm=1:num_str
%  strnames(mm,:)='            ';
%end
count=0;
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

%if (count~=num_str)
%  disp('error : numstrings not = number of grid connections');
%  return
%end
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
% add_src=ncread(fname,['/',str,'_weights.nc/add_src_address']);
% add_dst=ncread(fname,['/',str,'_weights.nc/add_dst_address']);
% add_remap=ncread(fname,['/',str,'_weights.nc/add_remap_matrix']);

%
  zd=double([1:dst_size(1)*dst_size(2)]*0.);
  zs=double([1:src_size(1)*src_size(2)]*0.);
  for mm=1:length(dst)
    zd(dst(mm))=zd(dst(mm))+remap(1,mm);
    zs(src(mm))=zs(src(mm))+remap(1,mm);
  end
%  for mm=1:length(add_dst)
%    zd(add_dst(mm))=zd(add_dst(mm))+add_remap(mm);
%    zs(add_src(mm))=zs(add_src(mm))+add_remap(mm);
%  end
  
  if (0)
    figure
    plot(zd,'b+')
    hold on
    plot(zs,'r+')
  end
  figure
  subplot(211)
  zs2=reshape(zs,src_size(1),src_size(2));
  src_lon(src_lon>pi)=src_lon(src_lon>pi)-2*pi;
  src_lon=reshape(src_lon,src_size(1),src_size(2))*180/pi;
  src_lat=reshape(src_lat,src_size(1),src_size(2))*180/pi;
  pcolorjw(src_lon,src_lat,zs2); colorbar
  title(str)
  subplot(212)
  zd2=reshape(zd,dst_size(1),dst_size(2));
  dst_lon(dst_lon>pi)=dst_lon(dst_lon>pi)-2*pi;
  dst_lon=reshape(dst_lon,dst_size(1),dst_size(2))*180/pi;
  dst_lat=reshape(dst_lat,dst_size(1),dst_size(2))*180/pi;
  dst_mask=double(reshape(dst_mask,dst_size(1),dst_size(2)));
  pcolorjw(dst_lon,dst_lat,zd2.*dst_mask); colorbar
  hold on
  %src_cor_lon(src_cor_lon>pi)=src_cor_lon(src_cor_lon>pi)-2*pi;
  %plot(src_cor_lon(:)*180/pi+360,src_cor_lat(:)*180/pi,'r+')
end

