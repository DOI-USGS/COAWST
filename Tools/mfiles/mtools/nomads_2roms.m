%   script nomads_2roms.m
%
%   input choices-  3hour NAM 12km data grib files
%                   3hour NARR 32 km grib files
%                   3hour GFS  0.5 deg grib files
%                   Can combine and interpolates to a common grid.
%                   All data is read thru THREDDs from
%                   https://www.ncdc.noaa.gov/nomads/data-products
% 
%   user selects:   - time interval [nam_start:nam_end]
%                   - spatial interval [roms grid or generic grid]
%                   - variables to be inlucded
%                   - to get 1) NAM and/or NARR, -- or -- 2) GFS.
%
%   output-     ROMS netcdf forcing file
%   needs-      - native matlab netcdf to create a forcing file for ROMS.
%               - native matlab to read opendap
%
%   Example of options to be used in project.h
%   # define BULK_FLUXES
%   # define ANA_BTFLUX
%   # define ANA_BSFLUX
%   # define ATM_PRESS
%   # define EMINUSP
%   # define SOLAR_SOURCE
%
% 24Sept2014 - jcwarner convert from Alfredo's + several others.
% 19Sept2016 - jcwarner add GFS read option
%                       add rotation of NAM and NARR from Lamb Conf to Lat Lon.
%

%%%%%%%%%%%%%%%%%%%%%   START OF USER INPUT  %%%%%%%%%%%%%%%%%%%%%%%%%%

%(1) Select which variables to include in this netcdf forcing file.
%  put a '1' if you want to include it, '0' otherwise.
%  1/0          Var description                   Units   File expected
get_lwrad=1;      % gets downward and upward longwave and computes net flux of down-up (W/m2)
                  % this will also store lwrad_down so you can use LONGWAVE option.
get_swrad=1;      % gets downward and updward shortwave and computes net down-up flux (W/m2)
get_rain=1;       % precipitation rate (kg/m2/s)
get_Tair=1;       % surface air temperature (C)
get_Pair=1;       % pressure at surface (Pa)
get_Qair=1;       % relative_humidity (percent)
get_Wind=1;       % surface u- and v- winds (m/s)

%(2) Enter name of output ROMS forcing file
ROMS_forc_name='roms_namnarr_Sandy2012_GFS.nc';

%(3) Enter start and end dates
namnarr_start = datenum('28-Oct-2012');
namnarr_end   = datenum('29-Oct-2012');

%(4) Select to interpolate to a roms grid or a user defined grid.
% Set one of these to a 1, the other to a 0.
interpto_roms_grid = 0;
interpto_user_grid = 1;
if (interpto_roms_grid)
  model_grid='Sandy_roms_grid.nc';
elseif (interpto_user_grid)
% Need to provide lon_rho, lat_rho, and angle_rho.
% NAM / NARR grids are centered at~ -100 deg lon; GFS = 0:360 lon
  lon_rho=[255:0.1:310];%-360;
  lat_rho=[ 10:0.1:50 ];  % Create a 0.1 degree lon-lat grid
  lon_rho=repmat(lon_rho,length(lat_rho),1)';
  lat_rho=repmat(lat_rho',1,size(lon_rho,1))';
  angle_rho = lon_rho*0;
else
  disp('pick a grid')
end

%5) Select which data to obtain: NAM, NARR, both NAM+NARR -- or -- GFS.
get_NARR=0;  %NARR-A grid 221 32km data
get_NAM=0;   %NAM grid 218 12km data
% --- or    ---
get_GFS=1;   %GFS 0.5 degree
% GFS is 6 hr and NAM/NARR is 3 hr. I dont have time interpolation
% added in so you have to pick NAM/NARR or GFS.
%
%%%%%%%%%%%%%%%%%%%%%   END OF USER INPUT  %%%%%%%%%%%%%%%%%%%%%%%%%%

%determine the grid to interpolate to
if (interpto_roms_grid)
  nc_model_grid= netcdf.open(model_grid);
  varid = netcdf.inqVarID(nc_model_grid,'lon_rho'); 
  lon_rho = netcdf.getVar(nc_model_grid,varid);
  varid = netcdf.inqVarID(nc_model_grid,'lat_rho'); 
  lat_rho = netcdf.getVar(nc_model_grid,varid);
  varid = netcdf.inqVarID(nc_model_grid,'angle'); 
  angle_rho = netcdf.getVar(nc_model_grid,varid);
end
[Lp,Mp]=size(lon_rho);
L=Lp-1;
M=Mp-1;

% now figure out what year they want
if ((get_NARR+get_NAM)>0)
  NAMNARR_time=[namnarr_start:3/24:namnarr_end];
  ntimes=length(NAMNARR_time);
  Time=NAMNARR_time-datenum(1858,11,17,0,0,0);
end
if (get_GFS>0)
  NAMNARR_time=[namnarr_start:6/24:namnarr_end];
  ntimes=length(NAMNARR_time);
  Time=NAMNARR_time-datenum(1858,11,17,0,0,0);
end
if ((get_NARR+get_NAM)>0 && (get_GFS)>0)
  disp('cant do GFS and NAM/NARR'); return;
end
% Creation of NetCDF file for NARR data
nc = netcdf.create(ROMS_forc_name,'nc_clobber');

% Global variables
netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'type', 'bulk fluxes forcing file');
netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'gridid','combined grid');
netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'history',['Created by "' mfilename '" on ' datestr(now)]);

% Dimensions
disp(' ## Defining Dimensions...')
dimid_time=netcdf.defDim(nc,'narr_time',length(Time));
lon_dimID=netcdf.defDim(nc,'xr',Lp);
lat_dimID=netcdf.defDim(nc,'er',Mp);
t=length(Time);
t_dimID = netcdf.defDim(nc,'time',t);

% Variables and attributes:
disp(' ## Defining Variables, and Attributes...')

tID = netcdf.defVar(nc,'ocean_time','double',t_dimID);
netcdf.putAtt(nc,tID,'long_name','atmospheric forcing time');
netcdf.putAtt(nc,tID,'units','days');
netcdf.putAtt(nc,tID,'field','time, scalar, series');

lonID = netcdf.defVar(nc,'lon','double',[lon_dimID lat_dimID]);
netcdf.putAtt(nc,lonID,'long_name','longitude');
netcdf.putAtt(nc,lonID,'units','degrees_east');
netcdf.putAtt(nc,lonID,'field','xp, scalar, series');

latID = netcdf.defVar(nc,'lat','double',[lon_dimID lat_dimID]);
netcdf.putAtt(nc,latID,'long_name','latitude');
netcdf.putAtt(nc,latID,'units','degrees_north');
netcdf.putAtt(nc,latID,'field','yp, scalar, series');

if (get_Wind)
  wt_dimID = netcdf.defDim(nc,'wind_time',t);
  wtID = netcdf.defVar(nc,'wind_time','double',wt_dimID);
  netcdf.putAtt(nc,wtID,'long_name','wind_time');
  netcdf.putAtt(nc,wtID,'units','days since 1858-11-17 00:00:00 UTC');
  netcdf.putAtt(nc,wtID,'units','days');
  netcdf.putAtt(nc,wtID,'field','wind_time, scalar, series');

  UwindID = netcdf.defVar(nc,'Uwind','double',[lon_dimID lat_dimID wt_dimID]);
  netcdf.putAtt(nc,UwindID,'long_name','surface u-wind component');
  netcdf.putAtt(nc,UwindID,'units','meter second-1');
  netcdf.putAtt(nc,UwindID,'field','Uwind, scalar, series');
  netcdf.putAtt(nc,UwindID,'coordinates','lon lat');
  netcdf.putAtt(nc,UwindID,'time','wind_time');

  VwindID = netcdf.defVar(nc,'Vwind','double',[lon_dimID lat_dimID wt_dimID]);
  netcdf.putAtt(nc,VwindID,'long_name','surface v-wind component');
  netcdf.putAtt(nc,VwindID,'units','meter second-1');
  netcdf.putAtt(nc,VwindID,'field','Vwind, scalar, series');
  netcdf.putAtt(nc,VwindID,'coordinates','lon lat');
  netcdf.putAtt(nc,VwindID,'time','wind_time');
end

if (get_Pair)
  Pat_dimID = netcdf.defDim(nc,'pair_time',t);
  PatID = netcdf.defVar(nc,'pair_time','double',Pat_dimID);
  netcdf.putAtt(nc,PatID,'long_name','pair_time');
  netcdf.putAtt(nc,PatID,'units','days since 1858-11-17 00:00:00 UTC');
  netcdf.putAtt(nc,PatID,'units','days');
  netcdf.putAtt(nc,PatID,'field','pair_time, scalar, series');

  PairID = netcdf.defVar(nc,'Pair','double',[lon_dimID lat_dimID Pat_dimID]);
  netcdf.putAtt(nc,PairID,'long_name','surface air pressure');
  netcdf.putAtt(nc,PairID,'units','millibar');
  netcdf.putAtt(nc,PairID,'field','Pair, scalar, series');
  netcdf.putAtt(nc,PairID,'coordinates','lon lat');
  netcdf.putAtt(nc,PairID,'time','pair_time');
end

if (get_Tair)
  Tat_dimID = netcdf.defDim(nc,'tair_time',t);
  TatID = netcdf.defVar(nc,'tair_time','double',Tat_dimID);
  netcdf.putAtt(nc,TatID,'long_name','tair_time');
  netcdf.putAtt(nc,TatID,'units','days since 1858-11-17 00:00:00 UTC');
  netcdf.putAtt(nc,TatID,'units','days');
  netcdf.putAtt(nc,TatID,'field','tair_time, scalar, series');

  TairID = netcdf.defVar(nc,'Tair','double',[lon_dimID lat_dimID Tat_dimID]);
  netcdf.putAtt(nc,TairID,'long_name','surface air temperature');
  netcdf.putAtt(nc,TairID,'units','Celsius');
  netcdf.putAtt(nc,TairID,'field','Tair, scalar, series');
  netcdf.putAtt(nc,TairID,'coordinates','lon lat');
  netcdf.putAtt(nc,TairID,'time','tair_time');
end

if (get_Qair)
  Qat_dimID = netcdf.defDim(nc,'qair_time',t);
  QatID = netcdf.defVar(nc,'qair_time','double',Qat_dimID);
  netcdf.putAtt(nc,QatID,'long_name','qair_time');
  netcdf.putAtt(nc,QatID,'units','days since 1858-11-17 00:00:00 UTC');
  netcdf.putAtt(nc,QatID,'units','days');
  netcdf.putAtt(nc,QatID,'field','qair_time, scalar, series');

  QairID = netcdf.defVar(nc,'Qair','double',[lon_dimID lat_dimID Qat_dimID]);
  netcdf.putAtt(nc,QairID,'long_name','surface air relative humidity');
  netcdf.putAtt(nc,QairID,'units','percentage');
  netcdf.putAtt(nc,QairID,'field','Qair, scalar, series');
  netcdf.putAtt(nc,QairID,'coordinates','lon lat');
  netcdf.putAtt(nc,QairID,'time','qair_time');
end

if (get_rain)
  rt_dimID = netcdf.defDim(nc,'rain_time',t);
  rtID = netcdf.defVar(nc,'rain_time','double',rt_dimID);
  netcdf.putAtt(nc,rtID,'long_name','rain_time');
  netcdf.putAtt(nc,rtID,'units','days since 1858-11-17 00:00:00 UTC');
  netcdf.putAtt(nc,rtID,'units','days');
  netcdf.putAtt(nc,rtID,'field','rain_time, scalar, series');

  rainID = netcdf.defVar(nc,'rain','double',[lon_dimID lat_dimID rt_dimID]);
  netcdf.putAtt(nc,rainID,'long_name','rain fall rate');
  netcdf.putAtt(nc,rainID,'units','kilogram meter-2 second-1');
  netcdf.putAtt(nc,rainID,'field','rain, scalar, series');
  netcdf.putAtt(nc,rainID,'coordinates','lon lat');
  netcdf.putAtt(nc,rainID,'time','rain_time');
end

if (get_swrad)
  swrt_dimID = netcdf.defDim(nc,'srf_time',t);
  swrtID = netcdf.defVar(nc,'srf_time','double',swrt_dimID);
  netcdf.putAtt(nc,swrtID,'long_name','srf_time');
  netcdf.putAtt(nc,swrtID,'units','days since 1858-11-17 00:00:00 UTC');
  netcdf.putAtt(nc,swrtID,'units','days');
  netcdf.putAtt(nc,swrtID,'field','srf_time, scalar, series');

  swradID = netcdf.defVar(nc,'swrad','double',[lon_dimID lat_dimID swrt_dimID]);
  netcdf.putAtt(nc,swradID,'long_name','net solar shortwave radiation');
  netcdf.putAtt(nc,swradID,'units','Watts meter-2');
  netcdf.putAtt(nc,swradID,'positive_value','downward flux, heating');
  netcdf.putAtt(nc,swradID,'negative_value','upward flux, cooling');
  netcdf.putAtt(nc,swradID,'field','swrad, scalar, series');
  netcdf.putAtt(nc,swradID,'coordinates','lon lat');
  netcdf.putAtt(nc,swradID,'time','srf_time');
end

if (get_lwrad)
  lwrt_dimID = netcdf.defDim(nc,'lrf_time',t);
  lwrtID = netcdf.defVar(nc,'lrf_time','double',lwrt_dimID);
  netcdf.putAtt(nc,lwrtID,'long_name','lrf_time');
  netcdf.putAtt(nc,lwrtID,'units','days since 1858-11-17 00:00:00 UTC');
  netcdf.putAtt(nc,lwrtID,'units','days');
  netcdf.putAtt(nc,lwrtID,'field','lrf_time, scalar, series');

  lwradID = netcdf.defVar(nc,'lwrad','double',[lon_dimID lat_dimID lwrt_dimID]);
  netcdf.putAtt(nc,lwradID,'long_name','net downward solar longwave radiation');
  netcdf.putAtt(nc,lwradID,'units','Watts meter-2');
  netcdf.putAtt(nc,lwradID,'positive_value','downward flux, heating');
  netcdf.putAtt(nc,lwradID,'negative_value','upward flux, cooling');
  netcdf.putAtt(nc,lwradID,'field','lwrad, scalar, series');
  netcdf.putAtt(nc,lwradID,'coordinates','lon lat');
  netcdf.putAtt(nc,lwradID,'time','lrf_time');

  lwradID = netcdf.defVar(nc,'lwrad_down','double',[lon_dimID lat_dimID lwrt_dimID]);
  netcdf.putAtt(nc,lwradID,'long_name','downward solar longwave radiation');
  netcdf.putAtt(nc,lwradID,'units','Watts meter-2');
  netcdf.putAtt(nc,lwradID,'positive_value','downward flux, heating');
  netcdf.putAtt(nc,lwradID,'negative_value','upward flux, cooling');
  netcdf.putAtt(nc,lwradID,'field','lwrad_down, scalar, series');
  netcdf.putAtt(nc,lwradID,'coordinates','lon lat');
  netcdf.putAtt(nc,lwradID,'time','lrf_time');
end
netcdf.close(nc)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% write lon and lat to netcdf file
%
nc=ROMS_forc_name;
ncwrite(nc,'lon',lon_rho);
ncwrite(nc,'lat',lat_rho);
ncwrite(nc,'ocean_time',Time);
%
% pre allocate some arrays
%
if (get_lwrad)
  lwrad=zeros(size(lon_rho,1),size(lon_rho,2),ntimes);
  lwrad_down=zeros(size(lon_rho,1),size(lon_rho,2),ntimes);
end
if (get_swrad)
  swrad=zeros(size(lon_rho,1),size(lon_rho,2),ntimes);
end
if(get_rain)
  rain=zeros(size(lon_rho,1),size(lon_rho,2),ntimes);
end
if (get_Tair)
  Tair=zeros(size(lon_rho,1),size(lon_rho,2),ntimes);
end
if (get_Pair)
  Pair=zeros(size(lon_rho,1),size(lon_rho,2),ntimes);
end
if (get_Qair)
  Qair=zeros(size(lon_rho,1),size(lon_rho,2),ntimes);
end
if (get_Wind)
  Uwind=zeros(size(lon_rho,1),size(lon_rho,2),ntimes);
  Vwind=zeros(size(lon_rho,1),size(lon_rho,2),ntimes);
end
%
if (get_GFS)
  disp('going to get GFS grid 4 0.5deg data');
%
  for mm=1:ntimes
    dd=datestr(Time(mm)+datenum(1858,11,17,0,0,0),'yyyymmddTHHMMSS');
    disp(['getting GFS grid 4 0.5deg data at ',dd]);
    url=['http://nomads.ncdc.noaa.gov/thredds/dodsC/gfs-004-anl/',dd(1:6),'/',dd(1:8),'/gfsanl_4_',dd(1:8),'_',dd(10:11),'00_000.grb2'];
    if (mm==1)
      x=ncread(url,'lon');
      y=ncread(url,'lat');
      [nlon,nlat]=meshgrid(x,y);
    end
%
    if (get_lwrad)
      down=squeeze(ncread(url,'Downward_Long-Wave_Rad_Flux'));
      down=down.';
      up=squeeze(ncread(url,'Upward_Long-Wave_Rad_Flux_surface'));
      up=up.';
      var=down-up;
      F=TriScatteredInterp(nlon(:),nlat(:),double(var(:)));
      cff=F(lon_rho,lat_rho);
      cff(isnan(cff))=0;
      lwrad(:,:,mm)=cff;
      F=TriScatteredInterp(nlon(:),nlat(:),double(down(:)));
      cff=F(lon_rho,lat_rho);
      cff(isnan(cff))=0;
      lwrad_down(:,:,mm)=cff;
    end
    if (get_swrad)
      down=squeeze(ncread(url,'Downward_Short-Wave_Rad_Flux'));
      down=down.';
      up=squeeze(ncread(url,'Upward_Short-Wave_Rad_Flux_surface'));
      up=up.';
      var=down-up;
      F=TriScatteredInterp(nlon(:),nlat(:),double(var(:)));
      cff=F(double(lon_rho),double(lat_rho));
      cff(isnan(cff))=0;
      swrad(:,:,mm)=cff;
    end
    if(get_rain)
      var=squeeze(ncread(url,'Precipitation_rate'));
      var=var.';
      F=TriScatteredInterp(nlon(:),nlat(:),double(var(:)));
      cff=F(double(lon_rho),double(lat_rho));
      cff(isnan(cff))=0;
      rain(:,:,mm)=cff;
    end
    if (get_Tair)
%     var=squeeze(ncread(url,'Temperature_surface'));
      var=squeeze(ncread(url,'Temperature_height_above_ground'));
      var=squeeze(var(:,:,1));
      var=var.';
      var=var-273.15; % Kelvin to Centigrades
      F=TriScatteredInterp(nlon(:),nlat(:),double(var(:)));
      cff=F(double(lon_rho),double(lat_rho));
      cff(isnan(cff))=0;
      Tair(:,:,mm)=cff;
    end
    if (get_Pair)
      var=squeeze(ncread(url,'Pressure_reduced_to_MSL'));
      var=var.';
      var=var*0.01; %Pa to db
      F=TriScatteredInterp(nlon(:),nlat(:),double(var(:)));
      cff=F(double(lon_rho),double(lat_rho));
      cff(isnan(cff))=0;
      Pair(:,:,mm)=cff;
    end
    if (get_Qair)
      var=squeeze(ncread(url,'Relative_humidity_height_above_ground'));
      var=var.';
      F=TriScatteredInterp(nlon(:),nlat(:),double(var(:)));
      cff=F(double(lon_rho),double(lat_rho));
      cff(isnan(cff))=0;
      Qair(:,:,mm)=cff;
    end
    if (get_Wind)
      var=squeeze(ncread(url,'U-component_of_wind_height_above_ground'));
      var=squeeze(var(:,:,1));
      var=var.';
      F=TriScatteredInterp(nlon(:),nlat(:),double(var(:)));
      cff=F(double(lon_rho),double(lat_rho));
      cff(isnan(cff))=0;
      Uwind_ll=cff;
  %
      var=squeeze(ncread(url,'V-component_of_wind_height_above_ground'));
      var=squeeze(var(:,:,1));
      var=var.';
      F=TriScatteredInterp(nlon(:),nlat(:),double(var(:)));
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
  save GFS_data.mat
end
if (get_NARR)
  disp('going to get NARR-A grid 221 32km data');
%
  for mm=1:ntimes
    dd=datestr(Time(mm)+datenum(1858,11,17,0,0,0),'yyyymmddTHHMMSS');
    disp(['getting NARR-A grid 221 32km data at ',dd]);
    url=['http://nomads.ncdc.noaa.gov/thredds/dodsC/narr/',dd(1:6),'/',dd(1:8),'/narr-a_221_',dd(1:8),'_',dd(10:11),'00_000.grb'];
    if (mm==1)
      x=ncread(url,'x');
      y=ncread(url,'y');
      clo=-107.0;   clat=50.0;
      earth_rad=6371.2;
      [X,Y]=meshgrid(x,y);
      m_proj('lambert conformal conic','clongitude',clo,'lat',[clat clat]);
      [nlon,nlat]=m_xy2ll(X/earth_rad,Y/earth_rad);
    end
%
    if (get_lwrad)
      down=squeeze(ncread(url,'Downward_longwave_radiation_flux'));
      down=down.';
      up=squeeze(ncread(url,'Upward_long_wave_radiation_flux_surface'));
      up=up.';
      var=down-up;
      F=TriScatteredInterp(nlon(:),nlat(:),double(var(:)));
      cff=F(lon_rho,lat_rho);
      cff(isnan(cff))=0;
      lwrad(:,:,mm)=cff;
      F=TriScatteredInterp(nlon(:),nlat(:),double(down(:)));
      cff=F(lon_rho,lat_rho);
      cff(isnan(cff))=0;
      lwrad_down(:,:,mm)=cff;
    end
    if (get_swrad)
      down=squeeze(ncread(url,'Downward_shortwave_radiation_flux'));
      down=down.';
      up=squeeze(ncread(url,'Upward_short_wave_radiation_flux_surface'));
      up=up.';
      var=down-up;
      F=TriScatteredInterp(nlon(:),nlat(:),double(var(:)));
      cff=F(double(lon_rho),double(lat_rho));
      cff(isnan(cff))=0;
      swrad(:,:,mm)=cff;
    end
    if(get_rain)
      var=squeeze(ncread(url,'Precipitation_rate'));
      var=var.';
      F=TriScatteredInterp(nlon(:),nlat(:),double(var(:)));
      cff=F(double(lon_rho),double(lat_rho));
      cff(isnan(cff))=0;
      rain(:,:,mm)=cff;
    end
    if (get_Tair)
%     var=squeeze(ncread(url,'Temperature_surface'));
      var=squeeze(ncread(url,'Temperature_height_above_ground'));
      var=squeeze(var(:,:,1));
      var=var.';
      var=var-273.15; % Kelvin to Centigrades
      F=TriScatteredInterp(nlon(:),nlat(:),double(var(:)));
      cff=F(double(lon_rho),double(lat_rho));
      cff(isnan(cff))=0;
      Tair(:,:,mm)=cff;
    end
    if (get_Pair)
      var=squeeze(ncread(url,'Pressure_reduced_to_MSL'));
      var=var.';
      var=var*0.01; %Pa to db
      F=TriScatteredInterp(nlon(:),nlat(:),double(var(:)));
      cff=F(double(lon_rho),double(lat_rho));
      cff(isnan(cff))=0;
      Pair(:,:,mm)=cff;
    end
    if (get_Qair)
      var=squeeze(ncread(url,'Relative_humidity'));
      var=var.';
      F=TriScatteredInterp(nlon(:),nlat(:),double(var(:)));
      cff=F(double(lon_rho),double(lat_rho));
      cff(isnan(cff))=0;
      Qair(:,:,mm)=cff;
    end
    if (get_Wind)
      var=squeeze(ncread(url,'u_wind_height_above_ground'));
      var=squeeze(var(:,:,1));
      var=var.';
      F=TriScatteredInterp(nlon(:),nlat(:),double(var(:)));
      cff=F(double(lon_rho),double(lat_rho));
      cff(isnan(cff))=0;
      Uwind_lamb=cff;
  %
      var=squeeze(ncread(url,'v_wind_height_above_ground'));
      var=squeeze(var(:,:,1));
      var=var.';
      F=TriScatteredInterp(nlon(:),nlat(:),double(var(:)));
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
if (get_NAM)
  disp('going to get NAM grid 218 12km data');
  %
  for mm=1:ntimes
    nstp=mod(mm,8);
    dd=datestr(Time(mm)+datenum(1858,11,17,0,0,0),'yyyymmddTHHMMSS');
    disp(['getting NAM grid 218 12km data at ',dd]);
          %http://nomads.ncdc.noaa.gov/thredds/dodsC/namanl/201210/20121028/namanl_218_20121028_0000_000.grb
    %url=['http://nomads.ncdc.noaa.gov/thredds/dodsC/namanl/',dd(1:6),'/',dd(1:8),'/namanl_218_',dd(1:8),'_',dd(10:11),'00_000.grb'];
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
    url=['http://nomads.ncdc.noaa.gov/thredds/dodsC/namanl/',dd(1:6),'/',dd(1:8),'/namanl_218_',dd(1:8),'_',first,'_',second,'.grb'];
    try
        if (mm==1)
          x=ncread(url,'x');
          y=ncread(url,'y');
          clo=-95.0;   clat=25.0;
          earth_rad=6371.2;
          [X,Y]=meshgrid(x,y);
          m_proj('lambert conformal conic','clongitude',clo,'lat',[clat clat]);
          [nlon,nlat]=m_xy2ll(X/earth_rad,Y/earth_rad);
      %
      % find the indices of the lon_rho lat_rho grid that are inside the NAM
      % data. we will just use these points from NAM and take the rest from NARR.
      %
      %   disp('computing mask to merge NARR and NAM')
          mask=zeros(size(lon_rho));
          X=[nlon(:,1); nlon(end,:)' ;nlon(end:-1:1,end); nlon(1,end:-1:1)'];
          Y=[nlat(:,1); nlat(end,:)' ;nlat(end:-1:1,end); nlat(1,end:-1:1)'];
          zz=inpolygon(lon_rho,lat_rho,X, Y);
          mask(zz==1)=1;
        end
      %
        if (get_lwrad)
          down=squeeze(ncread(url,'Downward_long_wave_flux'));
          down=down.';
          up=squeeze(ncread(url,'Upward_long_wave_flux'));
          up=up.';
          var=down-up;
          F=TriScatteredInterp(nlon(:),nlat(:),double(var(:)));
          zz=F(lon_rho,lat_rho);
          zz(isnan(zz))=0;
          cff=squeeze(lwrad(:,:,mm)).*(1-mask)+zz.*mask;
          lwrad(:,:,mm)=cff;
      %
          F=TriScatteredInterp(nlon(:),nlat(:),double(down(:)));
          zz=F(lon_rho,lat_rho);
          zz(isnan(zz))=0;
          cff=squeeze(lwrad_down(:,:,mm)).*(1-mask)+zz.*mask;
          lwrad_down(:,:,mm)=cff;
        end
        if (get_swrad)
          down=squeeze(ncread(url,'Downward_short_wave_flux'));
          down=down.';
          up=squeeze(ncread(url,'Upward_short_wave_flux'));
          up=up.';
          var=down-up;
          F=TriScatteredInterp(nlon(:),nlat(:),double(var(:)));
          zz=F(lon_rho,lat_rho);
          zz(isnan(zz))=0;
          cff=squeeze(swrad(:,:,mm)).*(1-mask)+zz.*mask;
          swrad(:,:,mm)=cff;
        end
        if(get_rain)
      %    var=squeeze(ncread(url,'Total_precipitation'));
      %    var=var.';
      %    F=TriScatteredInterp(nlon(:),nlat(:),double(var(:)));
      %    zz=F(lon_rho,lat_rho);
      %    zz(isnan(zz))=0;
      %    cff=squeeze(rain(:,:,mm)).*(1-mask)+zz.*mask;
      %    rain(:,:,mm)=cff;
        end
        if (get_Tair)
%         var=squeeze(ncread(url,'Temperature_surface'));
          var=squeeze(ncread(url,'Temperature_height_above_ground'));
          var=squeeze(var(:,:,1));
          var=var.';
          var=var-273.15; % Kelvin to Centigrades
          F=TriScatteredInterp(nlon(:),nlat(:),double(var(:)));
          zz=F(lon_rho,lat_rho);
          zz(isnan(zz))=0;
          cff=squeeze(Tair(:,:,mm)).*(1-mask)+zz.*mask;
          Tair(:,:,mm)=cff;
        end
        if (get_Pair)
          var=squeeze(ncread(url,'Pressure_reduced_to_MSL'));
          var=var.';
          var=var*0.01; %Pa to db
          F=TriScatteredInterp(nlon(:),nlat(:),double(var(:)));
          zz=F(lon_rho,lat_rho);
          zz(isnan(zz))=0;
          cff=squeeze(Pair(:,:,mm)).*(1-mask)+zz.*mask;
          Pair(:,:,mm)=cff;
        end
        if (get_Qair)
          var=squeeze(ncread(url,'Relative_humidity_height_above_ground'));
          var=var.';
          F=TriScatteredInterp(nlon(:),nlat(:),double(var(:)));
          zz=F(lon_rho,lat_rho);
          zz(isnan(zz))=0;
          cff=squeeze(Qair(:,:,mm)).*(1-mask)+zz.*mask;
          Qair(:,:,mm)=cff;
        end
        if (get_Wind)
          var=squeeze(ncread(url,'u_wind_height_above_ground'));
          var=squeeze(var(:,:,1));
          var=var.';
          F=TriScatteredInterp(nlon(:),nlat(:),double(var(:)));
          cff=F(lon_rho,lat_rho);
          cff(isnan(cff))=0;
          Uwind_lamb=cff;
      %
          var=squeeze(ncread(url,'v_wind_height_above_ground'));
          var=squeeze(var(:,:,1));
          var=var.';
          F=TriScatteredInterp(nlon(:),nlat(:),double(var(:)));
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
      disp(['could not get that data at ', url])
    end
  end
end
%
% write data to netcdf file
%
nc=ROMS_forc_name;
if (get_lwrad)
  ncwrite(nc,'lrf_time',Time);
  ncwrite(nc,'lwrad',lwrad);
  ncwrite(nc,'lwrad_down',lwrad_down);
end
if (get_swrad)
  ncwrite(nc,'srf_time',Time);
  ncwrite(nc,'swrad',swrad);
end
if(get_rain)
  ncwrite(nc,'rain_time',Time);
  rain(rain<0)=0;
  ncwrite(nc,'rain',rain);
end
if (get_Tair)
  ncwrite(nc,'tair_time',Time);
  Tair(Tair<-100)=0;
  ncwrite(nc,'Tair',Tair);
end
if (get_Pair)
  ncwrite(nc,'pair_time',Time);
  Pair(Pair<0)=0;
  ncwrite(nc,'Pair',Pair);
end
if (get_Qair)
  ncwrite(nc,'qair_time',Time);
  Qair(Qair<0)=0;
  ncwrite(nc,'Qair',Qair);
end
if (get_Wind)
  ncwrite(nc,'wind_time',Time);
  Uwind(Uwind<-100)=0;
  Vwind(Vwind<-100)=0;
  ncwrite(nc,'Uwind',Uwind);
  ncwrite(nc,'Vwind',Vwind);
end
%
disp(['------------ wrote ',ROMS_forc_name,' ------------']);
%datevec(wind_time(1)+datenum(1858,11,17,0,0,0))
%http://nomads.ncdc.noaa.gov/thredds/dodsC/narr/201210/20121028/narr-a_221_20121028_0000_000.grb.html


