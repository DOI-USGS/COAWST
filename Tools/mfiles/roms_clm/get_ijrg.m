function get_ijrg(gn,url2)

if (0)
  % This is what i did to get the indices for the entire HYCOM grid.
  % url2='http://tds.hycom.org/thredds/dodsC/GLBa0.08/expt_90.9'; %/2011';
  display(url2);
  nc=ncgeodataset(url2);
  hycom_lon=zeros(3298,4500);
  hycom_lat=zeros(3298,4500);
  for jj=1220:3298
    jj
    hycom_lon(jj,:)=lon.data(jj,:);
    hycom_lat(jj,:)=lat.data(jj,:);
  end
  close(nc)
  %
  % I then saved every 4th point, so that the file was only 4 megs.
  % These values are only used to get a rough idea of what
  % chunk size of data to grab via opendap so that we dont
  % have to grab the whole thing.
   hycom_lon_small=hycom_lon(1:4:end,1:4:end);
   hycom_lat_small=hycom_lat(1:4:end,1:4:end);
   save hycom_lonlat_small_grid.mat hycom_lat_small hycom_lon_small
end

%  Load the hycom indices 
load hycom_lonlat_small_grid.mat   %  this is every 4th index

%get roms grid dims
disp('getting roms grid dimensions ...');
xl=min(min(gn.lon_rho));xr=max(max(gn.lon_rho));
yb=min(min(gn.lat_rho));yt=max(max(gn.lat_rho));

% optimize the chunk size to obtain from hycom
disp('optimizing grid dimensions ...');
xg=hycom_lon_small;
xg(xg>=180)=(xg(xg>=180)-360);
yg=hycom_lat_small;
clear hycom_lon_small hycom_lat_small

[ym xm]=size(xg);
%ii=find(xg>=xl & xg<=xr & yg>=yb & yg<=yt);jj=fix((ii-1)/ym)+1;ii=mod(ii,ym);
[ii,jj] = find(xg>=xl & xg<=xr & yg>=yb & yg<=yt);

ig0=(min(ii)-1); ig1=(max(ii)+1); jg0=(min(jj)-1); jg1=(max(jj)+1);
%multiply all these *4 since i saved every 4th point.
ig0=(ig0-1)*4;
ig1=(ig1+1)*4;
jg0=(jg0-1)*4;
jg1=(jg1+1)*4;
% Constrain indexes to lie within the full HYCOM grid.
ig0 = max(ig0, 1);
jg0 = max(jg0, 1);
ig1 = min(ig1, 3298);
jg1 = min(jg1, 4000);

%save irg.mat irg
%save jrg.mat jrg

irg2=[num2str(ig0) ':' num2str(ig1)];
jrg2=[num2str(jg0) ':' num2str(jg1)];
firstidx=['[',num2str(ig0),' ',num2str(jg0),']']; 
lastidx=['[',num2str(ig1),' ',num2str(jg1),']'];

%url2='http://tds.hycom.org/thredds/dodsC/GLBa0.08/expt_90.9'; %/2011';
nc=ncgeodataset(url2);
eval(['lon=nc.data(''Longitude'',',firstidx,',',lastidx,');']);
xg=lon;
xg(xg>=180)=(xg(xg>=180)-360);
clear lon
eval(['lat=nc.data(''Latitude'',',firstidx,',',lastidx,');']);
yg=lat;
clear lat
clm.lon=double(xg);
clm.lat=double(yg);
lev=nc.data('Depth');
clm.z=double(lev); clear lev;
save hycom_info.mat clm irg2 jrg2 firstidx lastidx

close(nc)
clear xg yg ii jj;

