%   script narrnc2roms.m
%
%   input-  reads 1-year-long 3-hour-averaged NARR data netcdf files
%           Go here-->  from ftp://ftp.cdc.noaa.gov/Datasets/NARR/monolevel/
%           and get these files:
%    downward longwave radiation flux     dlwrf.YYYY.nc
%    downward shortwave radiation flux    dswrf.YYYY.nc
%    surface_upward_latent_heat_flux      lhtfl.YYYY.nc   (maybe)
%    surface_upward_sensible_heat_flux    shtfl.YYYY.nc   (maybe)
%    accum total evaporation amount       evap.YYYY.nc    (maybe)
%    precipitation rate                   prate.YYYY.nc
%    surface air temperature              air.sfc.YYYY.nc
%    pressure at surface                  pres.sfc.YYYY.nc
%    relative_humidity                    rhum.2m.YYYY.nc
%    surface u-wind component             uwnd.10m.YYYY.nc
%    surface v-wind component             vwnd.10m.YYYY.nc
%    these variable definitions are decsribed here:
%    http://www.esrl.noaa.gov/psd/data/gridded/data.narr.monolevel.html
%
%   user selects:   - time interval [narr_start:narr_end]
%                   - spatial interval [roms grid or generic grid]
%                   - variables to be inlucded
%
%   output- uses native matlab netcdf to create a forcing file for ROMS.
%
%   Example of options to be used in project.h
%   # define BULK_FLUXES
%   # define ANA_BTFLUX
%   # define ANA_BSFLUX
%   # define ATM_PRESS
%   # define LONGWAVE_OUT
%   # define EMINUSP
%   # define SOLAR_SOURCE
%
% 2012 - 04: I. Safak converted from the non-native netcdf setup of N.Kumar
% 23May2012 - jcwarner, modify slightly to be user defined.
% 12Sept2014 - jcwarner combine narr2roms + create_roms_forcings
%

%%%%%%%%%%%%%%%%%%%%%   START OF USER INPUT  %%%%%%%%%%%%%%%%%%%%%%%%%%

%(1) Select which variables to include in this netcdf forcing file.
%  put a '1' if you want to include it, '0' otherwise.
%  1/0          Var description                   Units   File expected
lwrad_down=1; % downward longwave radiation flux  W/m2    dlwrf.YYYY.nc
swrad=1;      % downward shortwave radiation flux W/m2    dswrf.YYYY.nc
latent=0;     % surface_upward_latent_heat_flux   W/m2    lhtfl.YYYY.nc
sensible=0;   % surface_upward_sensible_heat_flux W/m2    shtfl.YYYY.nc
evap=0;       % accum total evaporation amount    kg/m2   evap.YYYY.nc
rain=1;       % precipitation rate                kg/m2/s prate.YYYY.nc
Tair=1;       % surface air temperature           K       air.sfc.YYYY.nc
Pair=1;       % pressure at surface               Pa      pres.sfc.YYYY.nc
Qair=1;       % relative_humidity                 %       rhum.2m.YYYY.nc
Winds=1;      % surface u-wind                    m/s     uwnd.10m.YYYY.nc
              % surface v-wind                    m/s     vwnd.10m.YYYY.nc
% selecting 'Winds' above creates both u- and v-, both files needed.
%
% The air temp is converted from K to C below.
% This gets both downward and upward long wave an computes the net as
% net = down - up.  You can change this if you want.


%(2) Enter name of output ROMS forcing file
ROMS_NARR_name='roms_narr_Sandy2012.nc';

%(3) Enter start and end dates
narr_start = datenum('27-Oct-2012');
narr_end   = datenum('2-Nov-2012');

%(4) Select to interpolate to a roms grid or a user defined grid.
% Set one of these to a 1, the other to a 0.
interpto_roms_grid = 1;
interpto_user_grid = 0;
if (interpto_roms_grid)
  model_grid='Sandy_roms_grid.nc';
else
  lon_rho=[255:0.25:310]-360;
  lat_rho=[ 10:0.25:50 ];  % Create a 1/4 degree lat-lon grid
end

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
else
  lon_rho=repmat(lon_rho',1,length(lat_rho));
  lat_rho=repmat(lat_rho,size(lon_rho,1),1);
end
[Lp,Mp]=size(lon_rho);
L=Lp-1;
M=Mp-1;

% now figure out what year they want
zz=datevec(narr_start);
YYYY=num2str(zz(:,1));

%now lets figure out what time indices to get.
if (lwrad_down)
  eval(['ncfile = netcdf.open(''dlwrf.',YYYY,'.nc'',''NC_NOWRITE'');'])
elseif (swrad)
  eval(['ncfile = netcdf.open(''dswrf.',YYYY,'.nc'',''NC_NOWRITE'');'])
elseif (latent)
  eval(['ncfile = netcdf.open(''lhtfl.',YYYY,'.nc'',''NC_NOWRITE'');'])
elseif (sensible)
  eval(['ncfile = netcdf.open(''shtfl.',YYYY,'.nc'',''NC_NOWRITE'');'])
elseif (evap)
  eval(['ncfile = netcdf.open(''evap.',YYYY,'.nc'',''NC_NOWRITE'');'])
elseif (rain)
  eval(['ncfile = netcdf.open(''prate.',YYYY,'.nc'',''NC_NOWRITE'');'])
elseif (Tair)
  eval(['ncfile = netcdf.open(''air.sfc.',YYYY,'.nc'',''NC_NOWRITE'');'])
elseif (Pair)
  eval(['ncfile = netcdf.open(''pres.sfc.',YYYY,'.nc'',''NC_NOWRITE'');'])
elseif (Qair)
  eval(['ncfile = netcdf.open(''rhum.2m.',YYYY,'.nc'',''NC_NOWRITE'');'])
elseif (Winds)
  eval(['ncfile = netcdf.open(''uwnd.10m.',YYYY,'.nc'',''NC_NOWRITE'');'])
% eval(['ncfile = netcdf.open(''vwnd.10m.',YYYY,'.nc'',''NC_NOWRITE'');'])
end
varid = netcdf.inqVarID(ncfile,'time'); 
NARR_time = netcdf.getVar(ncfile,varid);
NARR_time = NARR_time./24+julian(1800,1,1,0);
NARR_time = datenum(gregorian(NARR_time));

%compute start and end wind indices
narr_start_idx = floor(interp1(NARR_time,[1:1:length(NARR_time)],narr_start));
narr_end_idx   = floor(interp1(NARR_time,[1:1:length(NARR_time)],narr_end));
ntimes=narr_end_idx-narr_start_idx+1;
Time = NARR_time(narr_start_idx:narr_end_idx)-datenum(1858,11,17,0,0,0);

% Creation of NetCDF file for NARR data
nc = netcdf.create(ROMS_NARR_name,'nc_clobber');

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

if (Winds)
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

if (Pair)
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

if (Tair)
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

if (Qair)
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

if (rain)
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

if (evap)
  et_dimID = netcdf.defDim(nc,'evap_time',t);
  etID = netcdf.defVar(nc,'evap_time','double',et_dimID);
  netcdf.putAtt(nc,etID,'long_name','evap_time');
  netcdf.putAtt(nc,etID,'units','days since 1858-11-17 00:00:00 UTC');
  netcdf.putAtt(nc,etID,'units','days');
  netcdf.putAtt(nc,etID,'field','evap_time, scalar, series');

  evapID = netcdf.defVar(nc,'evaporation','double',[lon_dimID lat_dimID et_dimID]);
  netcdf.putAtt(nc,evapID,'long_name','evaporation rate');
  netcdf.putAtt(nc,evapID,'units','kilogram meter-2 second-1');
  netcdf.putAtt(nc,evapID,'field','evap, scalar, series');
  netcdf.putAtt(nc,evapID,'coordinates','lon lat');
  netcdf.putAtt(nc,evapID,'time','evap_time');
end

if (swrad)
  swrt_dimID = netcdf.defDim(nc,'srf_time',t);
  swrtID = netcdf.defVar(nc,'srf_time','double',swrt_dimID);
  netcdf.putAtt(nc,swrtID,'long_name','srf_time');
  netcdf.putAtt(nc,swrtID,'units','days since 1858-11-17 00:00:00 UTC');
  netcdf.putAtt(nc,swrtID,'units','days');
  netcdf.putAtt(nc,swrtID,'field','srf_time, scalar, series');

  swradID = netcdf.defVar(nc,'swrad','double',[lon_dimID lat_dimID swrt_dimID]);
  netcdf.putAtt(nc,swradID,'long_name','solar shortwave radiation');
  netcdf.putAtt(nc,swradID,'units','Watts meter-2');
  netcdf.putAtt(nc,swradID,'positive_value','downward flux, heating');
  netcdf.putAtt(nc,swradID,'negative_value','upward flux, cooling');
  netcdf.putAtt(nc,swradID,'field','swrad, scalar, series');
  netcdf.putAtt(nc,swradID,'coordinates','lon lat');
  netcdf.putAtt(nc,swradID,'time','srf_time');
end

if (lwrad_down)
  lwrt_dimID = netcdf.defDim(nc,'lrf_time',t);
  lwrtID = netcdf.defVar(nc,'lrf_time','double',lwrt_dimID);
  netcdf.putAtt(nc,lwrtID,'long_name','lrf_time');
  netcdf.putAtt(nc,lwrtID,'units','days since 1858-11-17 00:00:00 UTC');
  netcdf.putAtt(nc,lwrtID,'units','days');
  netcdf.putAtt(nc,lwrtID,'field','lrf_time, scalar, series');

  lwradID = netcdf.defVar(nc,'lwrad_down','double',[lon_dimID lat_dimID lwrt_dimID]);
  netcdf.putAtt(nc,lwradID,'long_name','downward solar longwave radiation');
  netcdf.putAtt(nc,lwradID,'units','Watts meter-2');
  netcdf.putAtt(nc,lwradID,'positive_value','downward flux, heating');
  netcdf.putAtt(nc,lwradID,'negative_value','upward flux, cooling');
  netcdf.putAtt(nc,lwradID,'field','lwrad, scalar, series');
  netcdf.putAtt(nc,lwradID,'coordinates','lon lat');
  netcdf.putAtt(nc,lwradID,'time','lrf_time');
end

if (latent)
  lhft_dimID = netcdf.defDim(nc,'lhf_time',t);
  lhftID = netcdf.defVar(nc,'lhf_time','double',lhft_dimID);
  netcdf.putAtt(nc,lhftID,'long_name','latent heat flux time');
  netcdf.putAtt(nc,lhftID,'units','days since 1858-11-17 00:00:00 UTC');
  netcdf.putAtt(nc,lhftID,'units','days');
  netcdf.putAtt(nc,lwrtID,'field','lhf_time, scalar, series');

  latentID = netcdf.defVar(nc,'latent','double',[lon_dimID lat_dimID lhft_dimID]);
  netcdf.putAtt(nc,latentID,'long_name','latent heat flux');
  netcdf.putAtt(nc,latentID,'units','Watts meter-2');
  netcdf.putAtt(nc,latentID,'field','latent, scalar, series');
  netcdf.putAtt(nc,latentID,'coordinates','lon lat');
  netcdf.putAtt(nc,latentID,'time','lhf_time');
end

if (sensible)
  shft_dimID = netcdf.defDim(nc,'shf_time',t);
  shftID = netcdf.defVar(nc,'shf_time','double',shft_dimID);
  netcdf.putAtt(nc,shftID,'long_name','sensible heat flux time');
  netcdf.putAtt(nc,shftID,'units','days since 1858-11-17 00:00:00 UTC');
  netcdf.putAtt(nc,shftID,'units','days');
  netcdf.putAtt(nc,swrtID,'field','shf_time, scalar, series');

  sensibleID = netcdf.defVar(nc,'sensible','double',[lon_dimID lat_dimID shft_dimID]);
  netcdf.putAtt(nc,sensibleID,'long_name','sensible heat flux');
  netcdf.putAtt(nc,sensibleID,'units','Watts meter-2');
  netcdf.putAtt(nc,sensibleID,'field','sensible, scalar, series');
  netcdf.putAtt(nc,sensibleID,'coordinates','lon lat');
  netcdf.putAtt(nc,sensibleID,'time','shf_time');
end
netcdf.close(nc)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% write lon and lat to netcdf file
%
nc=ROMS_NARR_name;
ncwrite(nc,'lon',lon_rho);
ncwrite(nc,'lat',lat_rho);
ncwrite(nc,'ocean_time',Time);
%
disp('going to get narr data');
%
num_loops=lwrad_down+swrad+latent+sensible+evap+rain+Tair+Pair+Qair+Winds;
FIELDS=[''];
if (lwrad_down);FIELDS=[FIELDS;'lwrad_down ']; end
if (swrad);     FIELDS=[FIELDS;'swrad      ']; end
if (latent);    FIELDS=[FIELDS;'latent     ']; end
if (sensible);  FIELDS=[FIELDS;'sensible   ']; end
if (evap);      FIELDS=[FIELDS;'evaporation']; end
if (rain);      FIELDS=[FIELDS;'rain       ']; end
if (Tair);      FIELDS=[FIELDS;'Tair       ']; end
if (Pair);      FIELDS=[FIELDS;'Pair       ']; end
if (Qair);      FIELDS=[FIELDS;'Qair       ']; end
if (Winds);     FIELDS=[FIELDS;'Winds      ']; end

for mm=1:num_loops
  expr=FIELDS(mm,:);
  vector=0;
  switch expr
    case('lwrad_down ')
      ncv=['dlwrf.',YYYY,'.nc']; varname='dlwrf'; time_name='lrf_time';  scale=1; offset=0;
    case('swrad      ')
      ncv=['dswrf.',YYYY,'.nc']; varname='dswrf'; time_name='srf_time';  scale=1; offset=0;
    case('latent     ')
      ncv=['lhtfl.',YYYY,'.nc']; varname='lhtfl'; time_name='lhf_time';  scale=-1; offset=0;
    case('sensible   ')
      ncv=['shtfl.',YYYY,'.nc']; varname='shtfl'; time_name='shf_time';  scale=-1; offset=0;
    case('evaporation')
      ncv=['evap.',YYYY,'.nc']; varname='evap'; time_name='evap_time';  scale=1; offset=0;
    case('rain       ')
      ncv=['prate.',YYYY,'.nc']; varname='prate'; time_name='rain_time';  scale=1; offset=0;
    case('Tair       ')
      ncv=['air.sfc.',YYYY,'.nc']; varname='air'; time_name='tair_time';  scale=1; offset=-273.15;
    case('Pair       ')
      ncv=['pres.sfc.',YYYY,'.nc']; varname='pres'; time_name='pair_time';  scale=0.01; offset=0;
    case('Qair       ')
      ncv=['rhum.2m.',YYYY,'.nc']; varname='rhum'; time_name='qair_time';  scale=1; offset=0;
    case('Winds      ')
      ncv =['uwnd.10m.',YYYY,'.nc']; varname='uwnd'; time_name='wind_time';  scale=1; offset=0;
      ncv2=['vwnd.10m.',YYYY,'.nc']; varname2='vwnd'; vector=1;
  end

  nlon = ncread(ncv,'lon');   nlat = ncread(ncv,'lat');
  numx=size(nlon,1);          numy=size(nlon,2);
  var = ncread(ncv,varname,[1 1 narr_start_idx],[numx numy ntimes]);
  var=var*scale+offset;
  if (vector)  %Winds
    var2=ncread(ncv2,varname2,[1 1 narr_start_idx],[numx numy ntimes]);
    var2=var2*scale+offset;
  end

% Pre-allocate variables that will contain wind data on the ROMS grid
  var_roms=zeros(Lp,Mp,ntimes);
  if (vector)  %Winds
    var2_roms=zeros(Lp,Mp,ntimes);
  end
  count = 0;
% Interpolate NARR data onto the ROMS grid
  for i = 1:ntimes
      count=count+1;
      disp(['NARR data for variable ', varname, ' at ',datestr(NARR_time(narr_start_idx-1+i))])

      un = var(:,:,i);
      un(abs(un)>99999)=0;
      F  = TriScatteredInterp(double(nlon(:)),double(nlat(:)),un(:));
      var_roms(:,:,count) = F(lon_rho,lat_rho);

      if (vector)  %Winds
        un = var2(:,:,i);
        un(abs(un)>99999)=0;
        F  = TriScatteredInterp(double(nlon(:)),double(nlat(:)),un(:));
        var2_roms(:,:,count) = F(lon_rho,lat_rho);
        if (interpto_roms_grid)
          cffx= var_roms(:,:,count).*cos(angle_rho)+var2_roms(:,:,count).*sin(angle_rho);
          cffy=var2_roms(:,:,count).*cos(angle_rho)-var_roms(:,:,count).*sin(angle_rho);
          var_roms(:,:,count)=cffx;
          var2_roms(:,:,count)=cffy;
        end
      end
  end
  ncwrite(nc,time_name,Time);
  if (vector)  %Winds
    ncwrite(nc,'Uwind',var_roms);
    ncwrite(nc,'Vwind',var2_roms);
    clear var2_roms
  else
    ncwrite(nc,strtrim(FIELDS(mm,:)),var_roms);
  end
end

disp(['------------ wrote ',ROMS_NARR_name,' ------------']);
%datevec(wind_time(1)+datenum(1858,11,17,0,0,0))



