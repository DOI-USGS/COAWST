function [gn, clm]=get_ijrg(url, modelgrid, theta_s, theta_b, Tcline, N, Vtransform, Vstretching)
%
% here we get the indices from the hycom grid,
% compare them to the roms grid,
% and then just determine a subset of the hycom grid to obtain data from.
%
% jcw, revised, Feb 10, 2019
%

%
% Read ROMS grid info
%
disp('getting roms grid dimensions ...');
Sinp.theta_s     =theta_s;      %surface control parameter
Sinp.theta_b     =theta_b;      %bottom  control parameter
Sinp.Tcline      =Tcline;       %surface/bottom stretching width
Sinp.N           =N;            %number of vertical levels
Sinp.Vtransform  =Vtransform;   %vertical transformation equation
Sinp.Vstretching =Vstretching;  %vertical stretching function
if (Vtransform==1)
  h=ncread(modelgrid,'h');
  hmin=min(h(:));
  hc=min(max(hmin,0),Tcline);
elseif (Vtransform==2)
  hc=Tcline;
end
Sinp.hc          =hc;           %stretching width used in ROMS
gn=get_roms_grid(modelgrid,Sinp);
gn.z_r=shiftdim(gn.z_r,2);
gn.z_u=shiftdim(gn.z_u,2);
gn.z_v=shiftdim(gn.z_v,2);
gn.z_w=shiftdim(gn.z_w,2);

%
% Read HYCOM lon lat depth
%
display(['getting HYCOM grid data from ', url])
lonlat_fullgrid=1;
try
  numX=ncread(url,'X');
  numY=ncread(url,'Y');
% hycom_lon=ncread(url,'Longitude',[1 1],[length(numX) 1]);
% hycom_lat=ncread(url,'Latitude',[1 1],[1 length(numY)]);
  hycom_lon=ncread(url,'Longitude');
  hycom_lat=ncread(url,'Latitude');
  hycom_depth=ncread(url,'Depth');
catch
  hycom_lon=ncread(url,'lon');
  hycom_lat=ncread(url,'lat');
  hycom_depth=ncread(url,'depth');
  lonlat_fullgrid=0;
  numX=length(hycom_lon);
  numY=length(hycom_lat);
end

%
% Get roms grid limits
%
disp('getting roms grid dimensions ...');
xl=min(min(gn.lon_rho));xr=max(max(gn.lon_rho));
yb=min(min(gn.lat_rho));yt=max(max(gn.lat_rho));
%
% optimize the chunk size to obtain from hycom
%
disp('optimizing grid dimensions ...');
%
% now use xg and yg becasue we are modifying the lon
%
xg=hycom_lon;
%for exp930 dont need this next line
xg(xg>=180)=xg(xg>=180)-360;
yg=hycom_lat;
%
% Find the indices of the roms grid (xl xr yb yt) that are inside the 
% hycom grid (xg yg)
%
%  lon and Longitude are the same in the lat direction.
%  lat changes in the longitude dir for older data
[ii] = find(xg(:,1)>=xl & xg(:,1)<=xr);
ig0=ii(1);   ig1=ii(end);
%
if (lonlat_fullgrid)
  [jj] = find(yg(1200,:)>=yb & yg(1200,:)<=yt);
  jg0=jj(1);   jg1=jj(end);
  if (yt>20)
    jg1=numY(end);   %just take all the way north
  end
else
  [jj] = find(yg>=yb & yg<=yt);
  jg0=jj(1);   jg1=jj(end);
end
%
% Constrain indexes to lie within the full HYCOM grid.
%
ig0 = max(ig0-1, 1);
jg0 = max(jg0-1, 1);
ig1 = min(ig1+1, numX(end));
jg1 = min(jg1+1, numY(end));
%
% also save indices as strings
%
irg2=[num2str(ig0) ':' num2str(ig1)];
jrg2=[num2str(jg0) ':' num2str(jg1)];
%
if (lonlat_fullgrid)
  clm.lon=double(xg(ig0:ig1,jg0:jg1));
  clm.lat=double(yg(ig0:ig1,jg0:jg1));
else
  %clm.lon=double(xg(ig0:ig1));
  %%for exp930 need this one
  %%clm.lat=double(yg(jg0:jg1)');
  clm.lon=double(repmat(xg(ig0:ig1),1,jg1-jg0+1));
  clm.lat=double(repmat(yg(jg0:jg1)',ig1-ig0+1,1));
end

clm.z=double(hycom_depth);
clm.irg2=irg2;
clm.jrg2=jrg2;
clm.ig0=ig0; % lon idx strt
clm.ig1=ig1; % lon idx end
clm.jg0=jg0; % lat idx strt
clm.jg1=jg1; % lat idx end
