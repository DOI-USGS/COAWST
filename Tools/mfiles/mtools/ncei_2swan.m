%   script ncei_2swan.m
%
%   inputs-     3hour NAM 12km data grib files
%               3hour NARR 32 km grib files
%               6 hour GFS  0.5 deg grib files
%               Can combine and interpolates to a common grid.
%               All data is read thru THREDDs from
%               https://www.ncei.noaa.gov/thredds/model/model.html
% 
%   user selects:   - time interval [nam_start:nam_end]
%                   - spatial interval [roms grid or generic grid]
%                   - variables to be inlucded
%                   - to get NAM, NARR, GFS.
%                   NAM seems to be from  2017 to present
%                   NARR and GFS seem to be for longer time periods. 
%                   Please check the data is available for your time period!
%
%   output-     SWAN ASCII wind forcing file
%   needs-      native matlab to read opendap
%
% 29Sept2014 - jcwarner
% 19Sept2016 - jcwarner add GFS read option
%                       add rotation of NAM and NARR from Lamb Conf to Lat Lon.
% 04Feb2019  - jcwarner update name from nomads_2swan to be ncei_2swan
%                       use https://www.ncei.noaa.gov ...
%
% 22Dec2023 - jcwarner:  allow use of GFS + NAM by adding time averaging of GFS to 3 hr.
%                        you can use GFS, NARR, NAM, NARR+NAM, or GFS+NAM
echo off
%%%%%%%%%%%%%%%%%%%%%   START OF USER INPUT  %%%%%%%%%%%%%%%%%%%%%%%%%%

%(1) Select which variables to include in this ascii forcing file.
%  put a '1' if you want to include it, '0' otherwise.
get_Wind=1;       % surface u- and v- winds (m/s)

%(2) Enter name of output SWAN forcing file
SWAN_forc_name='swan_GFS_Oct2012.dat';

%(3) Enter start and end dates - needed if using get_NARR or get_NAM
namnarr_start = datenum('28-Oct-2012');
namnarr_end   = datenum('01-Nov-2012');

%(4) Select which data to obtain: NAM, NARR, or GFS.
get_NARR=0;  %NARR-A grid 221 32km data, dt = 3 hr
get_NAM=0;   %NAM grid 218 12km data, dt = 3 hr
get_GFS=1;   %GFS 0.5 degree, dt = 6 hr
%
% -- or  --
%
% read in data from a netcdf file
get_netcdf_file=0;
if (get_netcdf_file)
  ncfile='roms_namnarr_30Sept10Oct2012.nc';
end

%(5) Select to interpolate to a swan grid or a user defined grid.
% Set one of these to a 1, the other to a 0.
interpto_swan_grid = 0;
interpto_user_grid = 1;
%
if (interpto_swan_grid)
  model_grid='Sandy_roms_grid.coord';
  nx=84;
  ny=24;
elseif (interpto_user_grid)
% Provide lon_rho, lat_rho, and angle_rho.
% NAM / NARR grids are centered near -100 deg lon; GFS = 0:360 lon.
% Always select global lon values (0-360) for lon_rho.
% These offsets will allow the mixed use of 
% NAM-or-NARR, or NAM-NARR + GFS to use 
% the west longitde convention (-100 for example).
  if (get_GFS);  offset=0;    end
  if (get_NARR); offset=-360; end
  if (get_NAM);  offset=-360; end
% Select grid resolution.
% I used 0.25 for Sandy case, but you should use 0.1.
  lon_rho = [255:0.1:310]+offset;      % always use global values in the [0:360]
  lat_rho = [ 10:0.1:50 ];
  lon_rho = repmat(lon_rho,size(lat_rho,2),1)';
  lat_rho = repmat(lat_rho',1,size(lon_rho,1))';
  angle_rho = lon_rho*0;
else
  disp('getting data from file')
end
%%%%%%%%%%%%%%%%%%%%%   END OF USER INPUT  %%%%%%%%%%%%%%%%%%%%%%%%%%

if (get_NARR + get_NAM + get_GFS) > 0
  %get some grid info
  if (interpto_swan_grid)
    grd=textread(model_grid);
    gridsize=length(grd)/2;
    zz=grd(1:gridsize);
    lon_rho=reshape(zz,nx,ny);
    zz=grd(gridsize+1:end);
    lat_rho=reshape(zz,nx,ny);
  end
  [Lp,Mp]=size(lon_rho);
  L=Lp-1;
  M=Mp-1;
end

% now figure out what year they want
if (get_NARR + get_NAM + get_GFS) > 0
  NAMNARR_time=[namnarr_start:3/24:namnarr_end];
  ntimes=length(NAMNARR_time);
  Time=NAMNARR_time-datenum(1858,11,17,0,0,0);
  nskip=2;
end
if ((get_GFS > 0) && (get_NARR + get_NAM < 1))
  disp('doing GFS only')
  NAMNARR_time=[namnarr_start:6/24:namnarr_end];
  ntimes=length(NAMNARR_time);
  Time=NAMNARR_time-datenum(1858,11,17,0,0,0);
  nskip=1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% pre allocate some arrays
%
if (get_NARR + get_NAM + get_GFS) > 0
  if (get_Wind)
    Uwind=zeros(size(lon_rho,1),size(lon_rho,2),ntimes);
    Vwind=zeros(size(lon_rho,1),size(lon_rho,2),ntimes); 
  end
end
%
if (get_GFS)
  disp('going to get GFS grid 4 0.5deg data');
  C=datevec(time_start);
  gotF=0;   %this is the F = scattered interp func, only need to call it once
  for mm = 1:nskip:ntimes
        
        dd = datestr(time(mm) + datenum(1858,11,17,0,0,0),'yyyymmddTHHMMSS');
%
%       this is old
%       url = ['https://www.ncei.noaa.gov/thredds/dodsC/gfs-g4-anl-files/',dd(1:6),'/',dd(1:8),'/gfsanl_4_',dd(1:8),'_',dd(10:11),'00_000.grb2'];
%
%       if you go here , you will see these options:
%       https://www.ncei.noaa.gov/thredds/catalog/model/gfs.html
%
        if (C(1)>=2021)
%         GFS grid 4  file units  2020-2023                   no rad
%         https://www.ncei.noaa.gov/thredds/dodsC/model-gfs-004-files/202312/20231217/gfs_4_20231217_0000_000.grb2
          url = ['https://www.ncei.noaa.gov/thredds/dodsC/model-gfs-004-files/',dd(1:6),'/',dd(1:8),'/gfs_4_',dd(1:8),'_',dd(10:11),'00_000.grb2'];
%
%         GFS grid 4  historical file units 2019-2020           yes rad   ??
%
%         GFS grid 4 Analysis Set file units 2020-2023      - no Radiations   dont use this one
%         https://www.ncei.noaa.gov/thredds/dodsC/model-gfs-g4-anl-files/202311/20231130/gfs_4_20231130_1800_006.grb2
%         url = ['https://www.ncei.noaa.gov/thredds/dodsC/model-gfs-g4-anl-files/',dd(1:6),'/',dd(1:8),'/gfs_4_',dd(1:8),'_',dd(10:11),'00_000.grb2'];
%
        else
%         GFS grid 4 Analysis Set historical file units 2004-2020 - yes radiations, use this one
%         https://www.ncei.noaa.gov/thredds/dodsC/model-gfs-g4-anl-files-old/202005/20200515/gfsanl_4_20200515_0600_006.grb2
          url = ['https://www.ncei.noaa.gov/thredds/dodsC/model-gfs-g4-anl-files-old/',dd(1:6),'/',dd(1:8),'/gfsanl_4_',dd(1:8),'_',dd(10:11),'00_000.grb2'];
%
        end

        if (mm==1)
          x=ncread(url,'lon');
          y=ncread(url,'lat');
          [nlon,nlat]=meshgrid(x,y);
        end
        nlon=double(nlon); nlat=double(nlat);
        if (get_NARR + get_NAM)>0 
          addoff = -offset;
        else
          addoff = 0;
        end
%
        if (get_Wind)
          var=squeeze(ncread(url,'u-component_of_wind_height_above_ground'));
          var=double(squeeze(var(:,:,1)));
          var=var.';
          if gotF == 0
            F = scatteredInterpolant(nlon(:),nlat(:),var(:));
            gotF=1;
          else
            F.Values=var(:);
          end
          cff=F(double(lon_rho),double(lat_rho));
          cff(isnan(cff))=0;
          Uwind_ll=cff;
      %
          var=squeeze(ncread(url,'v-component_of_wind_height_above_ground'));
          var=double(squeeze(var(:,:,1)));
          var=var.';
          if gotF == 0
            F = scatteredInterpolant(nlon(:),nlat(:),var(:));
            gotF=1;
          else
            F.Values=var(:);
          end
          cff=F(double(lon_rho),double(lat_rho));
          cff(isnan(cff))=0;
          Vwind_ll=cff;
      %
      %   Rotate winds to ROMS or user grid.
      %
          cffx=Uwind_ll.*cos(angle_rho)+Vwind_ll.*sin(angle_rho);
          cffy=Vwind_ll.*cos(angle_rho)-Uwind_ll.*sin(angle_rho);
          Uwind(:,:,mm)=cffx;
          Vwind(:,:,mm)=cffy;
        end
    end
  if (get_NARR + get_NAM) > 0
% here we need to interp to every other missing data, the 3 hr ones.  3 9 15 21
    for mm = 2:nskip:ntimes
      if get_Wind
        Uwind(:,:,mm) = 0.5*(Uwind(:,:,mm-1)+Uwind(:,:,mm+1));
        Vwind(:,:,mm) = 0.5*(Vwind(:,:,mm-1)+Vwind(:,:,mm+1));
      end
    end
  end
  save GFS_data.mat
end
if (get_NARR)
  disp('going to get NARR-A grid 221 32km data');
  gotF=0;   %this is the F = scattered interp func, only need to call it once
%
  for mm=1:ntimes
    dd=datestr(Time(mm)+datenum(1858,11,17,0,0,0),'yyyymmddTHHMMSS');
    disp(['getting NARR-A grid 221 32km data at ',dd]);
    url=['http://www.ncei.noaa.gov/thredds/dodsC/narr-a-files/',dd(1:6),'/',dd(1:8),'/narr-a_221_',dd(1:8),'_',dd(10:11),'00_000.grb'];
%
    if (mm==1)
      x=ncread(url,'x');
      y=ncread(url,'y');
      clo=-107.0;   clat=50.0;
      earth_rad=6371.2;
      [X,Y]=meshgrid(x,y);
      m_proj('lambert conformal conic','clongitude',clo,'lat',[clat clat]);
      [nlon,nlat]=m_xy2ll(X/earth_rad,Y/earth_rad);
      nlon=double(nlon); nlat=double(nlat);
    end
%
    if (get_Wind)
      var=squeeze(ncread(url,'u-component_of_wind_height_above_ground'));
      var=double(squeeze(var(:,:,1)));
      var=var.';
      if gotF == 0
        F = scatteredInterpolant(nlon(:),nlat(:),var(:));
        gotF=1;
      else
        F.Values=var(:);
      end
      cff=F(double(lon_rho),double(lat_rho));
      cff(isnan(cff))=0;
      Uwind_lamb=cff;
%
      var=squeeze(ncread(url,'v-component_of_wind_height_above_ground'));
      var=double(squeeze(var(:,:,1)));
      var=var.';
      if gotF == 0
        F = scatteredInterpolant(nlon(:),nlat(:),var(:));
        gotF=1;
      else
        F.Values=var(:);
      end
      cff=F(double(lon_rho),double(lat_rho));
      cff(isnan(cff))=0;
      Vwind_lamb=cff;
  %
  %   Rotate winds to earth lon lat based on http://ruc.noaa.gov/RUC.faq.html
  %
  %   ROTCON_P          R  WIND ROTATION CONSTANT, = 1 FOR POLAR STEREO
  %                         AND SIN(LAT_TAN_P) FOR LAMBERT CONFORMAL
  %   LON_XX_P          R  MERIDIAN ALIGNED WITH CARTESIAN X-AXIS(DEG)
  %   LAT_TAN_P         R  LATITUDE AT LAMBERT CONFORMAL PROJECTION
  %                         IS TRUE (DEG)
      lat_tan_p  =  clat;                    % 50.0 for NARR;
      lon_xx_p   =  clo;                     % -107.0 for NARR;
      rotcon_p   =  sin(lat_tan_p*pi/180);
      deg2rad=2*pi/360;
%
      angle2 = rotcon_p*(lon_rho-lon_xx_p)*deg2rad;
      sinx2 = sin(angle2);
      cosx2 = cos(angle2);
      Uwind_rot = cosx2.*Uwind_lamb+sinx2.*Vwind_lamb;
      Vwind_rot =-sinx2.*Uwind_lamb+cosx2.*Vwind_lamb;
  %
  %   Rotate winds to ROMS or user grid.
  %
      cffx=Uwind_rot.*cos(angle_rho)+Vwind_rot.*sin(angle_rho);
      cffy=Vwind_rot.*cos(angle_rho)-Uwind_rot.*sin(angle_rho);
      Uwind(:,:,mm)=cffx;
      Vwind(:,:,mm)=cffy;
    end
  end
  save NARR_data.mat
end
%
if (get_NAM)
  disp('going to get NAM grid 218 12km data');
  gotF=0;   %this is the F = scattered interp func, only need to call it once
%
  for mm=1:ntimes
    nstp=mod(mm,8);
    dd=datestr(Time(mm)+datenum(1858,11,17,0,0,0),'yyyymmddTHHMMSS');
    disp(['getting NAM grid 218 12km data at ',dd]);
    if ismember(nstp,[1 2])
      first='0000';
    elseif ismember(nstp,[3 4])
      first='0600';
    elseif ismember(nstp,[5 6])
      first='1200';
    elseif ismember(nstp,[7 0])
      first='1800';
    end
    if ismember(nstp,[1 3 5 7])
      second='000';
    else
      second='003';
    end
%   url=['http://nomads.ncdc.noaa.gov/thredds/dodsC/namanl/',dd(1:6),'/',dd(1:8),'/namanl_218_',dd(1:8),'_',first,'_',second,'.grb'];
    url = ['https://www.ncei.noaa.gov/thredds/dodsC/namanl/',dd(1:6),'/',dd(1:8),'/namanl_218_',dd(1:8),'_',first,'_',second,'.grb'];
    try
      if (mm==1)
        x=ncread(url,'x');
        y=ncread(url,'y');
        clo=-95.0;   clat=25.0;
        earth_rad=6371.2;
        [X,Y]=meshgrid(x,y);
        m_proj('lambert conformal conic','clongitude',clo,'lat',[clat clat]);
        [nlon,nlat]=m_xy2ll(X/earth_rad,Y/earth_rad);
        nlon=double(nlon); nlat=double(nlat);
    %
    % find the indices of the lon_rho lat_rho grid that are inside the NAM
    % data. we will just use these points from NAM and take the rest from NARR.
    %
        disp('computing mask to merge NARR and NAM')
        mask=zeros(size(lon_rho));
        X=[nlon(:,1); nlon(end,:)' ;nlon(end:-1:1,end); nlon(1,end:-1:1)'];
        Y=[nlat(:,1); nlat(end,:)' ;nlat(end:-1:1,end); nlat(1,end:-1:1)'];
        zz=inpolygon(lon_rho,lat_rho,X, Y);
        mask(zz==1)=1;
      end
    %
      if (get_Wind)
        var=squeeze(ncread(url,'u_wind_height_above_ground'));
        var=double(squeeze(var(:,:,1)));
        var=var.';
        if gotF == 0
          F = scatteredInterpolant(nlon(:),nlat(:),var(:));
          gotF=1;
        else
          F.Values=var(:);
        end
        cff=F(lon_rho,lat_rho);
        cff(isnan(cff))=0;
        Uwind_lamb=cff;
    %
        var=squeeze(ncread(url,'v_wind_height_above_ground'));
        var=double(squeeze(var(:,:,1)));
        var=var.';
        if gotF == 0
          F = scatteredInterpolant(nlon(:),nlat(:),var(:));
          gotF=1;
        else
          F.Values=var(:);
        end
        cff=F(lon_rho,lat_rho);
        cff(isnan(cff))=0;
        Vwind_lamb=cff;
  %
  %   Rotate winds to earth lon lat based on http://ruc.noaa.gov/RUC.faq.html
  %
  %   ROTCON_P          R  WIND ROTATION CONSTANT, = 1 FOR POLAR STEREO
  %                         AND SIN(LAT_TAN_P) FOR LAMBERT CONFORMAL
  %   LON_XX_P          R  MERIDIAN ALIGNED WITH CARTESIAN X-AXIS(DEG)
  %   LAT_TAN_P         R  LATITUDE AT LAMBERT CONFORMAL PROJECTION
  %                         IS TRUE (DEG)
        lat_tan_p  =  clat;                    % 25.0 for NAM;
        lon_xx_p   =  clo;                     % -95.0 for NAM;
        rotcon_p   =  sin(lat_tan_p*pi/180);
        deg2rad=2*pi/360;
%
        angle2 = rotcon_p*(lon_rho-lon_xx_p)*deg2rad;
        sinx2 = sin(angle2);
        cosx2 = cos(angle2);
        Uwind_rot = cosx2.*Uwind_lamb+sinx2.*Vwind_lamb;
        Vwind_rot =-sinx2.*Uwind_lamb+cosx2.*Vwind_lamb;
  %
  %   Rotate winds to ROMS or user grid and merge with previous data.
  %
        cffx=Uwind_rot.*cos(angle_rho)+Vwind_rot.*sin(angle_rho);
        cffy=Vwind_rot.*cos(angle_rho)-Uwind_rot.*sin(angle_rho);
        Uwind(:,:,mm)=squeeze(Uwind(:,:,mm)).*(1-mask)+cffx.*mask;
        Vwind(:,:,mm)=squeeze(Vwind(:,:,mm)).*(1-mask)+cffy.*mask;
      end
    catch ME
      disp(['cldnt get that data at ', url])
    end
  end
end
%
% write data to ascii file
%
if(get_netcdf_file)
  netcdf_load(ncfile);
  Time=wind_time;
  get_Wind=1;
end
%
fid = fopen(SWAN_forc_name,'w');
for i=1:length(Time)
  if (get_Wind)
    disp(['Writing winds for SWAN at ',datestr(Time(i)+datenum(1858,11,17,0,0,0))])
    uswan=squeeze(Uwind(:,:,i)');
    vswan=squeeze(Vwind(:,:,i)');
    fprintf(fid,'%10.2f\n',uswan');
    fprintf(fid,'%10.2f\n',vswan');
  end
end
fclose(fid);

%
disp(['------------ wrote ',SWAN_forc_name,' ------------']);


