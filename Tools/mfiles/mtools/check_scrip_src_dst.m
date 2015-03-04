 function stat=check_scrip_src_dst
%
% compare scrip src and dst nc files
% and see if points are coincident. 
% The test looks for less than 2.0e-7 rads
% =  1.1459e-05 degs =~ 1.3 m
% If so, then offset the dst by 2.0e-7 rads
%
% jcwarner 04March2015
%
src_center_lon=ncread('src.nc','grid_center_lon');
src_center_lat=ncread('src.nc','grid_center_lat');
src_corner_lon=ncread('src.nc','grid_corner_lon');
src_corner_lat=ncread('src.nc','grid_corner_lat');

dst_center_lon=ncread('dst.nc','grid_center_lon');
dst_center_lat=ncread('dst.nc','grid_center_lat');
dst_corner_lon=ncread('dst.nc','grid_corner_lon');
dst_corner_lat=ncread('dst.nc','grid_corner_lat');

[Ls,Ms]=size(src_center_lon);
[Ld,Md]=size(dst_center_lon);
for ii=1:Ls
  for jj=1:Ms
    lons=src_center_lon(ii,jj);
    lats=src_center_lat(ii,jj);
    for mm=1:Ld
      for nn=1:Md
        lond=dst_center_lon(mm,nn);
        latd=dst_center_lat(mm,nn);
        zz=min(abs(lons-lond),abs(lats-latd));
        if (zz< 2.0e-7)
          dst_center_lon=dst_center_lon+2.0e-7;
          dst_center_lat=dst_center_lat+2.0e-7;
          dst_corner_lon=dst_corner_lon+2.0e-7;
          dst_corner_lat=dst_corner_lat+2.0e-7;
          ncwrite('dst.nc','grid_center_lon',dst_center_lon);
          ncwrite('dst.nc','grid_center_lat',dst_center_lat);
          ncwrite('dst.nc','grid_corner_lon',dst_corner_lon);
          ncwrite('dst.nc','grid_corner_lat',dst_corner_lat);
          return
        end
      end
    end
  end
end



