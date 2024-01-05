% function create_roms_netcdf_bndry_file_thredds
%
% read USEast forecast from thredds and 
% create bc file.
%
% LIMITED TO ZETA, UBAR, and VBAR for now  !!!!!!!!!
%
% jcw 15Nov2022
%
clc;
echo off
clear

%1 enter your working dir here
cd /Users/shimakasaei/Desktop/ROMS/Gridding/Simulation/David

%2 enter the file name to be created
fn='dec22_bc_masked.nc';

%3 enter the grid to interpolate to (your roms grid)
grid.name='JB05_grd_masked.nc';

%4 enter number of grid layers you will use
grid.N=8;

%5 enter start and end dates, the forecast has 2 files for each day.
tstart=datenum(2022,12,20,0,0,0);%datenum(2021,08,28,23,0,0); %ida  % best to start at hour 0.
tend=datenum(2022,12,25,0,0,0);%datenum(2021,09,03,0,0,0); %ida    % best to end at hr 0 or hr 12.
%tend=datenum(2021,09,03,0,0,0);    % best to end at hr 0 or hr 12.

%forcing ramp time
ramp_hours=6;  %hours of ramp at start

%%%%%%%%%%%%%%%%%%%% end of user input here %%%%%%%%%%%%%%%%%%%%%%%%%%%

%get COAWST grid data
disp(['getting COAWST grid data at first time step'])
timey=datestr(tstart,'yyyymmdd');
eval(['filecoawst=''http://geoport.whoi.edu/thredds/dodsC/vortexfs1/usgs/Projects/COAWST/2021/coawst_us_',timey,'_01.nc''']);
zt=ncread(filecoawst,'ocean_time')/3600/24+datenum(1858,11,17);
%                                           you can set a smaller size if you want
lonr_use=ncread(filecoawst,'lon_rho');   Istrr=1; Iendr=size(lonr_use,1);
latr_use=ncread(filecoawst,'lat_rho');   Jstrr=1; Jendr=size(lonr_use,2);
lonu_use=ncread(filecoawst,'lon_u');     Istru=1; Iendu=size(lonu_use,1);
latu_use=ncread(filecoawst,'lat_u');     Jstru=1; Jendu=size(lonu_use,2);
lonv_use=ncread(filecoawst,'lon_v');     Istrv=1; Iendv=size(lonv_use,1);
latv_use=ncread(filecoawst,'lat_v');     Jstrv=1; Jendv=size(lonv_use,2);
%  resize the lons lats if user sets smaller size
lonr_use=lonr_use(Istrr:Iendr,Jstrr:Jendr);
latr_use=latr_use(Istrr:Iendr,Jstrr:Jendr);
lonu_use=lonu_use(Istru:Iendu,Jstru:Jendu);
latu_use=latu_use(Istru:Iendu,Jstru:Jendu);
lonv_use=lonv_use(Istrv:Iendv,Jstrv:Jendv);
latv_use=latv_use(Istrv:Iendv,Jstrv:Jendv);
%
h_use=ncread(filecoawst,'h');
maskr_use=ncread(filecoawst,'mask_rho');
masku_use=ncread(filecoawst,'mask_u');
maskv_use=ncread(filecoawst,'mask_v');
angler_use=ncread(filecoawst,'angle');
%
s_rho_use=ncread(filecoawst,'s_rho');
s_rho=ncread('JB05_grd_masked.nc','s_rho');
%N=8;   %%%%%%%%% Number of layers I guess (shima)

%get user grid
netcdf_load(grid.name)
%make fig of grids
figure
subplot(121)
  zz=h_use;zz(maskr_use==0)=nan;
  pcolor(lonr_use,latr_use,zz)
  hold on
  zz=h;zz(mask_rho==0)=nan;
  pcolor(lon_rho,lat_rho,h)
  colormap('jet')
  zoom on
  shading flat
subplot(122)
  zz=h_use;zz(maskr_use==0)=nan;
  pcolor(lonr_use,latr_use,zz)
  hold on
  zz=h;zz(mask_rho==0)=nan;
  pcolor(lon_rho,lat_rho,zz)
  colormap('jet')
  zoom on
  shading flat
  caxis([0 30])
  set(gca,'xlim',[-74.4 -73],'ylim',[39.8 41.2])

%first lets go get the data to find out how many time steps

timecount=0;
zcount=0;
ubcount=0;
vbcount=0;
numfiles=ceil(tend-tstart)*2
for mm=1:numfiles
  ct=mod(mm,2); suff=ct*1+(ct-1)*13;
  toadd=['00',num2str(suff)];toadd=toadd(end-1:end);
  timey=datestr(tstart+floor((mm-1)/2),'yyyymmdd');
  eval(['filecoawst=''http://geoport.whoi.edu/thredds/dodsC/vortexfs1/usgs/Projects/COAWST/2021/coawst_us_',timey,'_',toadd,'.nc'''])
  zt=ncread(filecoawst,'ocean_time')/3600/24+datenum(1858,11,17);
  existyes=length(zt>0);
  if (existyes)
  % time
    zt=ncread(filecoawst,'ocean_time');
    ntimes=length(zt);
    for tidx=1:ntimes
      timecount=timecount+1;
      ot(timecount)=zt(tidx);
    end
  % zeta
    disp(['getting COAWST zeta data at ', timey])
    zt=ncread(filecoawst,'zeta',[Istrr Jstrr 1],[Iendr-Istrr+1 Jendr-Jstrr+1 Inf]);
    for tidx=1:ntimes
      zcount=zcount+1;
      ztt=double(squeeze(zt(:,:,tidx))); ztt(maskr_use==0)=nan; ztt=maplev(ztt); ztt=ztt(:);
      if (zcount==1)
        Fr = scatteredInterpolant(lonr_use(:),latr_use(:),ztt);
      else
        Fr.Values = ztt;
      end
      zeta_south(:,zcount)=Fr(lon_rho(:,1),lat_rho(:,1));  %%%lon_rho is from input *.nc file
      zeta_north(:,zcount)=Fr(lon_rho(:,end),lat_rho(:,end));
      zeta_west(:,zcount)=Fr(lon_rho(1,:),lat_rho(1,:));
      zeta_east(:,zcount)=Fr(lon_rho(end,:),lat_rho(end,:));
    end
  % ubar and vbar
    disp(['getting COAWST ubar vbar data at ', timey])
    zu=ncread(filecoawst,'ubar',[Istru Jstru 1],[Iendu-Istru+1 Jendu-Jstru+1 Inf]);
    zv=ncread(filecoawst,'vbar',[Istrv Jstrv 1],[Iendv-Istrv+1 Jendv-Jstrv+1 Inf]);
    for tidx=1:ntimes
      ubcount=ubcount+1;
      vbcount=vbcount+1;
     %interpolate u v bar from coarse grid to local grid rho points
      ztu=double(squeeze(zu(:,:,tidx))); ztu(masku_use==0)=nan; ztu=maplev(ztu); ztu=ztu(:);
      if (ubcount==1)
        Fu = scatteredInterpolant(lonu_use(:),latu_use(:),ztu);
      else
        Fu.Values = ztu;
      end
      ztv=double(squeeze(zv(:,:,tidx))); ztv(maskv_use==0)=nan; ztv=maplev(ztv); ztv=ztv(:);
      if (vbcount==1)
        Fv = scatteredInterpolant(lonv_use(:),latv_use(:),ztv);
      else
        Fv.Values = ztv;
      end
      if (ubcount==1)
        Fa = scatteredInterpolant(lonr_use(:),latr_use(:),angler_use(:));
        angle_i=Fa(lon_rho,lat_rho);
      end
      ubar_i=Fu(lon_rho,lat_rho);
      vbar_i=Fv(lon_rho,lat_rho);
     %rotate to east and north on the local grid
      ubar_e=ubar_i.*cos(angle_i)-vbar_i.*sin(angle_i);
      vbar_n=ubar_i.*sin(angle_i)+vbar_i.*cos(angle_i);
     %rotate to new grid, still at rho points
      ubar_rot= ubar_e.*cos(angle)+vbar_n.*sin(angle);
      vbar_rot=-ubar_e.*sin(angle)+vbar_n.*cos(angle);
     %avg to u and v points on local grid
      ubar(:,:)=(ubar_rot(1:end-1,:)+ubar_rot(2:end,:))/2;
      vbar(:,:)=(vbar_rot(:,1:end-1)+vbar_rot(:,2:end))/2;
     %save to each 4 boundaries
      ubar_south(:,ubcount)=ubar(:,1);
      ubar_north(:,ubcount)=ubar(:,end);
      ubar_west(:,ubcount)=ubar(1,:);
      ubar_east(:,ubcount)=ubar(end,:);
      vbar_south(:,vbcount)=vbar(:,1);
      vbar_north(:,vbcount)=vbar(:,end);
      vbar_west(:,vbcount)=vbar(1,:);
      vbar_east(:,vbcount)=vbar(end,:);
    end
  end
end

save yourdata.mat

%now create the bndry file
t_clim=length(ot);
gn.N=grid.N;
gn.lon_rho=lon_rho;
create_roms_netcdf_bndry_mwUL(fn,gn,t_clim,t_clim)

%now fill the data
ot=ot/3600/24;
%datevec(ot(1)+datenum(1858,11,17))
ncwrite(fn,'zeta_time', ot);
ncwrite(fn,'v2d_time',  ot);
ncwrite(fn,'v3d_time',  ot);
ncwrite(fn,'salt_time', ot);
ncwrite(fn,'temp_time', ot);

%replace any nans from interpolation with 0's
zeta_south(isnan(zeta_south))=0;
zeta_north(isnan(zeta_north))=0;
zeta_east(isnan(zeta_east))=0;
zeta_west(isnan(zeta_west))=0;
%
ubar_south(isnan(ubar_south))=0;
ubar_north(isnan(ubar_north))=0;
ubar_east(isnan(ubar_east))=0;
ubar_west(isnan(ubar_west))=0;
%
vbar_south(isnan(vbar_south))=0;
vbar_north(isnan(vbar_north))=0;
vbar_east(isnan(vbar_east))=0;
vbar_west(isnan(vbar_west))=0;

%apply ramp and any other scaling
scale=[ot-ot(1)]*24;
scale=min(scale/ramp_hours,1);
for tidx=1:length(ot)
  zeta_north(:,tidx)=zeta_north(:,tidx).*scale(tidx);
  zeta_east(:,tidx)=zeta_east(:,tidx).*scale(tidx);
  zeta_south(:,tidx)=zeta_south(:,tidx).*scale(tidx);
  zeta_west(:,tidx)=zeta_west(:,tidx).*scale(tidx);

  ubar_north(:,tidx)=ubar_north(:,tidx).*scale(tidx);
  ubar_east(:,tidx)=ubar_east(:,tidx).*scale(tidx);
  ubar_south(:,tidx)=ubar_south(:,tidx).*scale(tidx);
  ubar_west(:,tidx)=ubar_west(:,tidx).*scale(tidx);

  vbar_north(:,tidx)=vbar_north(:,tidx).*scale(tidx);
  vbar_east(:,tidx)=vbar_east(:,tidx).*scale(tidx);
  vbar_south(:,tidx)=vbar_south(:,tidx).*scale(tidx);
  vbar_west(:,tidx)=vbar_west(:,tidx).*scale(tidx);
end

ncwrite(fn,'zeta_north',zeta_north);
ncwrite(fn,'zeta_east',zeta_east);
ncwrite(fn,'zeta_west',zeta_west);
ncwrite(fn,'zeta_south',zeta_south);
%
ncwrite(fn,'ubar_north',ubar_north);
ncwrite(fn,'ubar_east',ubar_east);
ncwrite(fn,'ubar_west',ubar_west);
ncwrite(fn,'ubar_south',ubar_south);
%
ncwrite(fn,'vbar_north',vbar_north);
ncwrite(fn,'vbar_east',vbar_east);
ncwrite(fn,'vbar_west',vbar_west);
ncwrite(fn,'vbar_south',vbar_south);

