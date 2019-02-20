%   script ncei_2roms.m
%
%   input choices-  3hour NAM 12km data grib files
%                   3hour NARR 32 km grib files
%                   3hour GFS  0.5 deg grib files
%                   can combine and interpolate to a common grid
%                   all data is read thru THREDDs from:
%                   https://www.ncei.noaa.gov/thredds/model/model.html
% 
%   user selects:   - time interval
%                   - spatial interval [roms grid or generic grid]
%                   - variables to include
%                   - to get 1) NAM and/or NARR, -- or -- 2) GFS
%                   NAM seems to be from  2017 to present
%                   NARR and GFS seem to be for longer time periods. 
%                   Please check the data is available for your time period!
%
%                   - use of 1) matlab native ncread -- or -- 2) nctoolbox
%
%   output-     ROMS netcdf forcing file
%
%   needs-      - native matlab netcdf to create a forcing file for ROMS.
%               - native matlab to read opendap
%               - optional: nctoolbox to read opendap
%                 see: https://github.com/nctoolbox/nctoolbox for
%                 installation and setup instructions
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
% 1Nov2017 - chegermiller update base url at lines to reflect use of NCEI server
%            instead of NOMADS server
% 2Nov2017 - chegermiller notes would be good to add try catch statement at urls to account for
%            the potential for .grb or .grb2 
%            For newer files, slight change of variable names, might fail
%            there too
% 9Nov2017 - chegermiller adds nctoolbox option based on Maitane's edits
%
% 2Nov2017 - chegermiller NOTE: for older matlab versions, this still might fail because
% ncread is looking for http: to indicate an OPeNDAP file. Here, the server
% is now https:

%%%%%%%%%%%%%%%%%%%%%%%%%%  START OF USER INPUT  %%%%%%%%%%%%%%%%%%%%%%%%%%

% (1) Select which variables to include in this netcdf forcing file.
% put a '1' if you want to include it, '0' otherwise.
% 1/0              Var description (Units)
get_lwrad = 1;    % gets downward and upward longwave and computes net flux of down-up (W/m2)
                  % this will also store lwrad_down so you can use LONGWAVE option.
get_swrad = 1;    % gets downward and updward shortwave and computes net down-up flux (W/m2)
get_rain  = 1;    % precipitation rate (kg/m2/s) at surface
get_Tair  = 1;    % surface air temperature (C) at 2m
get_Pair  = 1;    % pressure reduce to MSL (Pa)
get_Qair  = 1;    % relative_humidity (percent) at 2m
get_Wind  = 1;    % surface u- and v- winds (m/s) at 10m

% (2) Enter name of output ROMS forcing file
%ROMS_force_name = 'romsforc_GFS_Sandy2012.nc';
ROMS_force_name = 'romsforc_NARR_Sandy2012.nc';

% (3) Enter start and end dates
time_start = datenum('28-Oct-2012');
time_end   = datenum('31-Oct-2012');

% (4) Select which data to obtain: NAM, NARR, both NAM+NARR -- or -- GFS.
get_NARR = 1;  % NARR-A grid 221 32km data, available 1979-2014
get_NAM  = 0;  % NAM grid 218 12km data
% --- or ---
get_GFS  = 0;  % GFS 0.5 degree
% GFS is 6 hr and NAM/NARR is 3 hr. I dont have time interpolation
% added in so you have to pick NAM/NARR or GFS.

% (5) Select to interpolate to a roms grid or a user defined grid.
% Set one of these to a 1, the other to a 0.
interpto_roms_grid = 0;
interpto_user_grid = 1;

% The code has been written so that longitude is in -deg E. 
if interpto_roms_grid
  model_grid = 'C:\models\grids\USEast_grd31.nc';
elseif interpto_user_grid
% Provide lon_rho, lat_rho, and angle_rho.
% NAM / NARR grids are centered at~ -100 deg lon; GFS = 0:360 lon
  if (get_NARR); offset=-360; end
  if (get_NAM);  offset=-360; end
  if (get_GFS);  offset=0;    end
% You Probably want to make a finer resolution from 0.2 to 0.1.
  lon_rho = [255:0.2:310]+offset;
  lat_rho = [ 10:0.2:50 ];
  lon_rho = repmat(lon_rho,size(lat_rho,2),1)';
  lat_rho = repmat(lat_rho',1,size(lon_rho,1))';
  angle_rho = lon_rho*0;
else
  disp('pick a grid')
end

% (6) Select which method to use: ncread or nctoolbox
use_matlab = 1;
use_nctoolbox = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%  END OF USER INPUT  %%%%%%%%%%%%%%%%%%%%%%%%%%%

if use_nctoolbox
  setup_nctoolbox
end

% determine the grid to interpolate to
if interpto_roms_grid
  lon_rho = ncread(model_grid,'lon_rho');
  lat_rho = ncread(model_grid,'lat_rho');
  angle_rho = ncread(model_grid,'angle');
end
[Lp,Mp] = size(lon_rho);

% determine date range
if (get_NARR + get_NAM) > 0
  time = [time_start:3/24:time_end] - datenum(1858,11,17,0,0,0);
  ntimes = length(time);
end
if get_GFS > 0
  time = [time_start:6/24:time_end] - datenum(1858,11,17,0,0,0);
  ntimes = length(time);
end
if (get_NARR + get_NAM) > 0 && get_GFS > 0
  disp('cant do GFS and NAM/NARR')
  return
end

% create NetCDF file for forcing data
disp([' ## Creating ' ROMS_force_name '...'])
nc = netcdf.create(ROMS_force_name,'nc_clobber');

% global variables
netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'type', 'bulk fluxes forcing file');
netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'gridid','combined grid');
netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'history',['Created by "' mfilename '" on ' datestr(now)]);

% dimensions
disp(' ## Defining Dimensions...')
lon_dimID = netcdf.defDim(nc,'xr',Lp);
lat_dimID = netcdf.defDim(nc,'er',Mp);
t_dimID = netcdf.defDim(nc,'time',ntimes);

% variables and attributes
disp(' ## Defining Variables and Attributes...')

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

if get_Wind
  wt_dimID = netcdf.defDim(nc,'wind_time',ntimes);
  wtID = netcdf.defVar(nc,'wind_time','double',wt_dimID);
  netcdf.putAtt(nc,wtID,'long_name','wind_time');
  netcdf.putAtt(nc,wtID,'units','days since 1858-11-17 00:00:00 UTC');
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

if get_Pair
  Pat_dimID = netcdf.defDim(nc,'pair_time',ntimes);
  PatID = netcdf.defVar(nc,'pair_time','double',Pat_dimID);
  netcdf.putAtt(nc,PatID,'long_name','pair_time');
  netcdf.putAtt(nc,PatID,'units','days since 1858-11-17 00:00:00 UTC');
  netcdf.putAtt(nc,PatID,'field','pair_time, scalar, series');

  PairID = netcdf.defVar(nc,'Pair','double',[lon_dimID lat_dimID Pat_dimID]);
  netcdf.putAtt(nc,PairID,'long_name','surface air pressure');
  netcdf.putAtt(nc,PairID,'units','millibar');
  netcdf.putAtt(nc,PairID,'field','Pair, scalar, series');
  netcdf.putAtt(nc,PairID,'coordinates','lon lat');
  netcdf.putAtt(nc,PairID,'time','pair_time');
end

if get_Tair
  Tat_dimID = netcdf.defDim(nc,'tair_time',ntimes);
  TatID = netcdf.defVar(nc,'tair_time','double',Tat_dimID);
  netcdf.putAtt(nc,TatID,'long_name','tair_time');
  netcdf.putAtt(nc,TatID,'units','days since 1858-11-17 00:00:00 UTC');
  netcdf.putAtt(nc,TatID,'field','tair_time, scalar, series');

  TairID = netcdf.defVar(nc,'Tair','double',[lon_dimID lat_dimID Tat_dimID]);
  netcdf.putAtt(nc,TairID,'long_name','surface air temperature');
  netcdf.putAtt(nc,TairID,'units','Celsius');
  netcdf.putAtt(nc,TairID,'field','Tair, scalar, series');
  netcdf.putAtt(nc,TairID,'coordinates','lon lat');
  netcdf.putAtt(nc,TairID,'time','tair_time');
end

if get_Qair
  Qat_dimID = netcdf.defDim(nc,'qair_time',ntimes);
  QatID = netcdf.defVar(nc,'qair_time','double',Qat_dimID);
  netcdf.putAtt(nc,QatID,'long_name','qair_time');
  netcdf.putAtt(nc,QatID,'units','days since 1858-11-17 00:00:00 UTC');
  netcdf.putAtt(nc,QatID,'field','qair_time, scalar, series');

  QairID = netcdf.defVar(nc,'Qair','double',[lon_dimID lat_dimID Qat_dimID]);
  netcdf.putAtt(nc,QairID,'long_name','surface air relative humidity');
  netcdf.putAtt(nc,QairID,'units','percentage');
  netcdf.putAtt(nc,QairID,'field','Qair, scalar, series');
  netcdf.putAtt(nc,QairID,'coordinates','lon lat');
  netcdf.putAtt(nc,QairID,'time','qair_time');
end

if get_rain
  rt_dimID = netcdf.defDim(nc,'rain_time',ntimes);
  rtID = netcdf.defVar(nc,'rain_time','double',rt_dimID);
  netcdf.putAtt(nc,rtID,'long_name','rain_time');
  netcdf.putAtt(nc,rtID,'units','days since 1858-11-17 00:00:00 UTC');
  netcdf.putAtt(nc,rtID,'field','rain_time, scalar, series');

  rainID = netcdf.defVar(nc,'rain','double',[lon_dimID lat_dimID rt_dimID]);
  netcdf.putAtt(nc,rainID,'long_name','rain fall rate');
  netcdf.putAtt(nc,rainID,'units','kilogram meter-2 second-1');
  netcdf.putAtt(nc,rainID,'field','rain, scalar, series');
  netcdf.putAtt(nc,rainID,'coordinates','lon lat');
  netcdf.putAtt(nc,rainID,'time','rain_time');
end

if get_swrad
  swrt_dimID = netcdf.defDim(nc,'srf_time',ntimes);
  swrtID = netcdf.defVar(nc,'srf_time','double',swrt_dimID);
  netcdf.putAtt(nc,swrtID,'long_name','srf_time');
  netcdf.putAtt(nc,swrtID,'units','days since 1858-11-17 00:00:00 UTC');
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

if get_lwrad
  lwrt_dimID = netcdf.defDim(nc,'lrf_time',ntimes);
  lwrtID = netcdf.defVar(nc,'lrf_time','double',lwrt_dimID);
  netcdf.putAtt(nc,lwrtID,'long_name','lrf_time');
  netcdf.putAtt(nc,lwrtID,'units','days since 1858-11-17 00:00:00 UTC');
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

% write lon and lat
ncwrite(ROMS_force_name,'lon',lon_rho);
ncwrite(ROMS_force_name,'lat',lat_rho);
ncwrite(ROMS_force_name,'ocean_time',time);

% pre allocate some arrays
fill = zeros([size(lon_rho) ntimes]);

if get_lwrad
  lwrad = fill;
  lwrad_down = fill;
end
if get_swrad
  swrad = fill;
end
if get_rain
  rain = fill;
end
if get_Tair
  Tair = fill;
end
if get_Pair
  Pair = fill;
end
if get_Qair
  Qair = fill;
end
if get_Wind
  Uwind = fill;
  Vwind = fill;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if get_GFS
    disp('going to get GFS grid 4 0.5deg data')
    
    for mm = 1:ntimes
        
        dd = datestr(time(mm) + datenum(1858,11,17,0,0,0),'yyyymmddTHHMMSS');
        url = ['https://www.ncei.noaa.gov/thredds/dodsC/gfs-g4-anl-files/',dd(1:6),'/',dd(1:8),'/gfsanl_4_',dd(1:8),'_',dd(10:11),'00_000.grb2'];
        
        disp(['getting GFS grid 4 0.5 deg data at ',dd])
        if use_nctoolbox
            geo = ncgeodataset(url);
        end
        if mm == 1
            if use_matlab
                x = double(ncread(url,'lon'));
                y = double(ncread(url,'lat'));
            elseif use_nctoolbox
                x = double(geo{'lon'}(:))-360;
                y = double(geo{'lat'}(:));
            end
            [nlon,nlat] = meshgrid(x,y);
            nlon=nlon.'; nlat=nlat.';
        end
        if get_lwrad && mm == 1
            % down = double(ncread(url,'Downward_Long-Wave_Rad_Flux')).';
            % up   = double(ncread(url,'Upward_Long-Wave_Rad_Flux_surface')).';
            % var  = down-up;
            % F = scatteredInterpolant(nlon(:),nlat(:),var(:));
            % cff = F(lon_rho,lat_rho);
            % cff(isnan(cff)) = 0;
            % lwrad(:,:,mm) = cff;
            
            % F = scatteredInterpolant(nlon(:),nlat(:),down(:));
            % cff = F(lon_rho,lat_rho);
            % cff(isnan(cff)) = 0;
            % lwrad_down(:,:,mm) = cff;
            disp('lwrad is not available for GFS');
        end
        if get_swrad && mm == 1
            % down = double(ncread(url,'Downward_Short-Wave_Rad_Flux')).';
            % up   = double(ncread(url,'Upward_Short-Wave_Rad_Flux_surface')).';
            % var  = down-up;
            % F = scatteredInterpolant(nlon(:),nlat(:),var(:));
            % cff = F(lon_rho,lat_rho);
            % cff(isnan(cff)) = 0;
            % swrad(:,:,mm) = cff;
            disp('swrad is not available for GFS');
        end
        if get_rain && mm == 1
            % var = double(ncread(url,'Precipitation_rate')).';
            % F = scatteredInterpolant(nlon(:),nlat(:),var(:));
            % cff = F(lon_rho,lat_rho);
            % cff(isnan(cff)) = 0;
            % rain(:,:,mm) = cff;
            disp('rain is not available for GFS');
        end
        if get_Tair
            if use_matlab
                % var = double(ncread(url,'Temperature_surface'));
                var = double(ncread(url,'Temperature_height_above_ground'));
                var = squeeze(var(:,:,1)); % 1st index is 2 m above gound
            elseif use_nctoolbox
                var = double(squeeze(geo{'Temperature_height_above_ground'}(:)));
                var = squeeze(var(1,:,:));
            end 
            var = var - 273.15; % K to degC
            F = scatteredInterpolant(nlon(:),nlat(:),var(:));
            cff = F(lon_rho,lat_rho);
            cff(isnan(cff)) = 0;
            Tair(:,:,mm) = cff;
        end
        if get_Pair
            if use_matlab
                var = double(ncread(url,'Pressure_reduced_to_MSL_msl'));
            elseif use_nctoolbox
                var = double(squeeze(geo{'Pressure_reduced_to_MSL_msl'}(:)));
            end
            var = var*0.01; %Pa to db
            F = scatteredInterpolant(nlon(:),nlat(:),var(:));
            cff = F(lon_rho,lat_rho);
            cff(isnan(cff)) = 0;
            Pair(:,:,mm) = cff;
        end
        if get_Qair
            if use_matlab
                var = double(ncread(url,'Relative_humidity_height_above_ground'));
            elseif use_nctoolbox
                var = double(squeeze(geo{'Relative_humidity_height_above_ground'}(:)));
            end
            F = scatteredInterpolant(nlon(:),nlat(:),var(:));
            cff = F(lon_rho,lat_rho);
            cff(isnan(cff)) = 0;
            Qair(:,:,mm) = cff;
        end
        if get_Wind
            if use_matlab
                var = double(ncread(url,'u-component_of_wind_height_above_ground'));
                var = squeeze(var(:,:,1));
            elseif use_nctoolbox
                var = double(squeeze(geo{'u-component_of_wind_height_above_ground'}(:)));
                var = squeeze(var(1,:,:));
            end
            F = scatteredInterpolant(nlon(:),nlat(:),var(:));
            cff = F(lon_rho,lat_rho);
            cff(isnan(cff)) = 0;
            Uwind_ll = cff;

            if use_matlab
                var = double(ncread(url,'v-component_of_wind_height_above_ground'));
                var = squeeze(var(:,:,1));
            elseif use_nctoolbox
                var = double(squeeze(geo{'v-component_of_wind_height_above_ground'}(:)));
                var = squeeze(var(1,:,:));
            end
            F = scatteredInterpolant(nlon(:),nlat(:),var(:));
            cff = F(lon_rho,lat_rho);
            cff(isnan(cff)) = 0;
            Vwind_ll = cff;
            
            % rotate winds to ROMS or user grid
            cffx = Uwind_ll.*cos(angle_rho) + Vwind_ll.*sin(angle_rho);
            cffy = Vwind_ll.*cos(angle_rho) - Uwind_ll.*sin(angle_rho);
            Uwind(:,:,mm) = cffx;
            Vwind(:,:,mm) = cffy;
        end
        if use_nctoolbox
          close(geo);
        end
    end
    save GFS_data.mat
end

if get_NARR
    disp('going to get NARR-A grid 221 32km data')

    for mm = 1:ntimes
        
        dd = datestr(time(mm) + datenum(1858,11,17,0,0,0),'yyyymmddTHHMMSS');
%       url = ['https://nomads.ncdc.noaa.gov/thredds/dodsC/narr/',dd(1:6),'/',dd(1:8),'/narr-a_221_',dd(1:8),'_',dd(10:11),'00_000.grb'];
        url=['http://www.ncei.noaa.gov/thredds/dodsC/narr-a-files/',dd(1:6),'/',dd(1:8),'/narr-a_221_',dd(1:8),'_',dd(10:11),'00_000.grb'];

        disp(['getting NARR-A grid 221 32km data at ',dd])
        if use_nctoolbox
            geo = ncgeodataset(url);
        end
        if mm == 1
            if use_matlab
                x = double(ncread(url,'x'));
                y = double(ncread(url,'y'));
            elseif use_nctoolbox
                x = geo{'x'}(:);
                y = geo{'y'}(:);
            end
            clon = -107.0;
            clat = 50.0;
            earth_rad = 6371.2;
            [x,y] = meshgrid(x,y);
            m_proj('lambert conformal conic','clongitude',clon,'lat',[clat clat]);
            [nlon,nlat] = m_xy2ll(x/earth_rad,y/earth_rad);
        end
        if get_lwrad
            if use_matlab
                down = double(ncread(url,'Downward_longwave_radiation_flux')).';
                up   = double(ncread(url,'Upward_long_wave_radiation_flux_surface')).';
            elseif use_nctoolbox
                down = double(squeeze(geo{'Downward_longwave_radiation_flux'}(:)));
                up   = double(squeeze(geo{'Upward_long_wave_radiation_flux_surface'}(:)));
            end
            var = down-up;
            F = scatteredInterpolant(nlon(:),nlat(:),var(:));
            cff = F(lon_rho,lat_rho);
            cff(isnan(cff)) = 0;
            lwrad(:,:,mm) = cff;
            
            F = scatteredInterpolant(nlon(:),nlat(:),down(:));
            cff = F(lon_rho,lat_rho);
            cff(isnan(cff)) = 0;
            lwrad_down(:,:,mm) = cff;
        end
        if get_swrad
            if use_matlab
                down = double(ncread(url,'Downward_shortwave_radiation_flux')).';
                up   = double(ncread(url,'Upward_short_wave_radiation_flux_surface')).';
            elseif use_nctoolbox
                down = double(squeeze(geo{'Downward_shortwave_radiation_flux'}(:)));
                up   = double(squeeze(geo{'Upward_long_wave_radiation_flux_surface'}(:)));
            end 
            var = down-up;
            F = scatteredInterpolant(nlon(:),nlat(:),var(:));
            cff = F(lon_rho,lat_rho);
            cff(isnan(cff)) = 0;
            swrad(:,:,mm) = cff;
        end
        if get_rain 
            if use_matlab
                var = double(ncread(url,'Precipitation_rate')).';
            elseif use_nctoolbox
                var = double(squeeze(geo{'Precipitation_rate'}(:)));
            end 
            F = scatteredInterpolant(nlon(:),nlat(:),var(:));
            cff = F(lon_rho,lat_rho);
            cff(isnan(cff)) = 0;
            rain(:,:,mm) = cff;
        end
        if get_Tair
            if use_matlab
                % var = double(ncread(url,'Temperature_surface'));
                var = double(ncread(url,'Temperature_height_above_ground'));
                var = squeeze(var(:,:,1)).';
            elseif use_nctoolbox
                var = double(squeeze(geo{'Temperature_height_above_ground'}(:)));
                var = squeeze(var(1,:,:));
            end
            var = var - 273.15; % K to degC
            F = scatteredInterpolant(nlon(:),nlat(:),var(:));
            cff = F(lon_rho,lat_rho);
            cff(isnan(cff)) = 0;
            Tair(:,:,mm) = cff;
        end
        if get_Pair
            if use_matlab
                var = double(ncread(url,'Pressure_reduced_to_MSL')).';
            elseif use_nctoolbox
                var = double(squeeze(geo{'Pressure_reduced_to_MSL'}(:)));
            end
            var = var*0.01; %Pa to db
            F = scatteredInterpolant(nlon(:),nlat(:),var(:));
            cff = F(lon_rho,lat_rho);
            cff(isnan(cff)) = 0;
            Pair(:,:,mm) = cff;
        end
        if get_Qair
            if use_matlab
                var = double(ncread(url,'Relative_humidity')).';
            elseif use_nctoolbox
                var = double(squeeze(geo{'Relative_humidity'}(:)));
            end
            F = scatteredInterpolant(nlon(:),nlat(:),var(:));
            cff = F(lon_rho,lat_rho);
            cff(isnan(cff)) = 0;
            Qair(:,:,mm) = cff;
        end
        if get_Wind
            if use_matlab
                var = double(ncread(url,'u_wind_height_above_ground'));
                var = squeeze(var(:,:,1)).';
            elseif use_nctoolbox
                var = double(squeeze(geo{'u_wind_height_above_ground'}(:)));
                var = squeeze(var(1,:,:));
            end
            F = scatteredInterpolant(nlon(:),nlat(:),var(:));
            cff = F(lon_rho,lat_rho);
            cff(isnan(cff)) = 0;
            Uwind_lamb = cff;

            if use_matlab
                var = double(ncread(url,'v_wind_height_above_ground'));
                var = squeeze(var(:,:,1)).';
            elseif use_nctoolbox
                var = double(squeeze(geo{'v_wind_height_above_ground'}(:)));
                var = squeeze(var(1,:,:));
            end
            F = scatteredInterpolant(nlon(:),nlat(:),var(:));
            cff = F(lon_rho,lat_rho);
            cff(isnan(cff)) = 0;
            Vwind_lamb = cff;

            %   Rotate winds to earth lon lat based on http://ruc.noaa.gov/RUC.faq.html
            %
            %   ROTCON_P          R  WIND ROTATION CONSTANT, = 1 FOR POLAR STEREO
            %                         AND SIN(LAT_TAN_P) FOR LAMBERT CONFORMAL
            %   LON_XX_P          R  MERIDIAN ALIGNED WITH CARTESIAN X-AXIS(DEG)
            %   LAT_TAN_P         R  LATITUDE AT LAMBERT CONFORMAL PROJECTION
            %                         IS TRUE (DEG)
            lat_tan_p = clat;                    % 50.0 for NARR;
            lon_xx_p = clon;                     % -107.0 for NARR;
            rotcon_p = sin(lat_tan_p*pi/180);
            deg2rad = 2*pi/360;
            
            angle2 = rotcon_p*(lon_rho-lon_xx_p)*deg2rad;
            sinx2 = sin(angle2);
            cosx2 = cos(angle2);
            Uwind_rot =  cosx2.*Uwind_lamb+sinx2.*Vwind_lamb;
            Vwind_rot = -sinx2.*Uwind_lamb+cosx2.*Vwind_lamb;
            
            % rotate winds to ROMS or user grid
            cffx = Uwind_rot.*cos(angle_rho) + Vwind_rot.*sin(angle_rho);
            cffy = Vwind_rot.*cos(angle_rho) - Uwind_rot.*sin(angle_rho);
            Uwind(:,:,mm) = cffx;
            Vwind(:,:,mm) = cffy;
        end
        if use_nctoolbox
          close(geo);
        end
    end
    save NARR_data.mat
end

if get_NAM
    disp('going to get NAM grid 218 12km data')
    
    for mm = 1:ntimes
        
        dd = datestr(time(mm) + datenum(1858,11,17,0,0,0),'yyyymmddTHHMMSS');
        disp(['getting NAM grid 218 12km data at ',dd]);
        
        nstp = mod(mm,8);
        if ismember(nstp,[1 2])
            first = '0000';
        elseif ismember(nstp,[3 4])
            first = '0600';
        elseif ismember(nstp,[5 6])
            first = '1200';
        elseif ismember(nstp,[7 0])
            first = '1800';
        end
        if ismember(nstp,[1 3 5 7])
            second = '000';
        else
            second = '003';
        end
        
        %url = ['https://www.ncei.noaa.gov/thredds/dodsC/namanl/',dd(1:6),'/',dd(1:8),'/namanl_218_',dd(1:8),'_',first,'_',second,'.grb2'];
        url = ['https://www.ncei.noaa.gov/thredds/dodsC/namanl/',dd(1:6),'/',dd(1:8),'/namanl_218_',dd(1:8),'_',first,'_',second,'.grb'];
        
        try
            if use_nctoolbox
                geo = ncgeodataset(url);
            end
            if mm == 1
                if use_matlab
                    x = double(ncread(url,'x'));
                    y = double(ncread(url,'y'));
                elseif use_nctoolbox
                    x = double(geo{'x'}(:));
                    y = double(geo{'y'}(:));
                end
                clon = -95.0;
                % clon = 265.0;   
                clat = 25.0;
                earth_rad = 6371.2;
                [x,y] = meshgrid(x,y);
                m_proj('lambert conformal conic','clongitude',clon,'lat',[clat clat]);
                [nlon,nlat] = m_xy2ll(x/earth_rad,y/earth_rad);
                
                % find the indices of the lon_rho lat_rho grid that are inside the NAM
                % data. we will just use these points from NAM and take the rest from NARR
                
                disp('computing mask to merge NARR and NAM')
                mask = zeros(size(lon_rho));
                X = [nlon(:,1); nlon(end,:)'; nlon(end:-1:1,end); nlon(1,end:-1:1)'];
                Y = [nlat(:,1); nlat(end,:)'; nlat(end:-1:1,end); nlat(1,end:-1:1)'];
                zz = inpolygon(lon_rho,lat_rho,X, Y);
                mask(zz == 1) = 1;
            end
            if get_lwrad
                if use_matlab
                    down = double(ncread(url,'downward_long_wave_rad_flux_surface')).';
                    up = double(ncread(url,'upward_long_wave_rad_flux_surface')).';
                    % down = double(ncread(url,'Downward_Long-Wave_Radp_Flux_surface')).';
                    % up   = double(ncread(url,'Upward_Long-Wave_Radp_Flux_surface')).';
                elseif use_nctoolbox
                    down = double(squeeze(geo{'downward_long_wave_rad_flux_surface'}(:)));
                    up   = double(squeeze(geo{'upward_long_wave_rad_flux_surface'}(:)));
                    % down = double(squeeze(geo{'Downward_Long-Wave_Radp_Flux_surface'}(:)));
                    % up   = double(squeeze(geo{'Upward_Long-Wave_Radp_Flux_surface'}(:)));
                end
                var = down-up;
                F = scatteredInterpolant(nlon(:),nlat(:),var(:));
                zz = F(lon_rho,lat_rho);
                zz(isnan(zz)) = 0;
                cff = squeeze(lwrad(:,:,mm)).*(1-mask)+zz.*mask;
                lwrad(:,:,mm)=cff;
                
                F = scatteredInterpolant(nlon(:),nlat(:),down(:));
                zz = F(lon_rho,lat_rho);
                zz(isnan(zz)) = 0;
                cff = squeeze(lwrad_down(:,:,mm)).*(1-mask)+zz.*mask;
                lwrad_down(:,:,mm) = cff;
            end
            if get_swrad
                if use_matlab
                    down = squeeze(ncread(url,'downward_short_wave_rad_flux_surface'));
                    up = squeeze(ncread(url,'upward_short_wave_rad_flux_surface'));
                    % down = double(ncread(url,'Downward_Short-Wave_Radiation_Flux_surface')).';
                    % up   = double(ncread(url,'Upward_Short-Wave_Radiation_Flux_surface')).';
                elseif use_nctoolbox
                    down = double(squeeze(geo{'downward_short_wave_rad_flux_surface'}(:)));
                    up   = double(squeeze(geo{'upward_short_wave_rad_flux_surface'}(:)));
                    % down = double(squeeze(geo{'Downward_Short-Wave_Radiation_Flux_surface'}(:)));
                    % up   = double(squeeze(geo{'Upward_Short-Wave_Radiation_Flux_surface'}(:)));
                end
                var = down-up;
                F = scatteredInterpolant(nlon(:),nlat(:),var(:));
                zz = F(lon_rho,lat_rho);
                zz(isnan(zz)) = 0;
                cff = squeeze(swrad(:,:,mm)).*(1-mask)+zz.*mask;
                swrad(:,:,mm) = cff;
            end
            if get_rain && mm == 1
                % var = double(ncread(url,'Total_precipitation_surface_0_Hour_Accumulation')).';
                % F = scatteredInterpolant(nlon(:),nlat(:),var(:));
                % zz = F(lon_rho,lat_rho);
                % zz(isnan(zz)) = 0;
                % cff = squeeze(rain(:,:,mm)).*(1-mask)+zz.*mask;
                % rain(:,:,mm) = cff;
                disp('rain is not available for NAM');
            end
            if get_Tair
                if use_matlab
                    var = double(ncread(url,'Temperature_height_above_ground'));
                    var = squeeze(var(:,:,1)).';
                elseif use_nctoolbox
                    var = double(squeeze(geo{'Temperature_height_above_ground'}(:)));
                    var = squeeze(var(1,:,:));
                end
                var = var - 273.15; % K to degC
                F = scatteredInterpolant(nlon(:),nlat(:),var(:));
                zz = F(lon_rho,lat_rho);
                zz(isnan(zz)) = 0;
                cff = squeeze(Tair(:,:,mm)).*(1-mask)+zz.*mask;
                Tair(:,:,mm) = cff;
            end
            if get_Pair
                if use_matlab
                    var = double(ncread(url,'Pressure_reduced_to_MSL_msl')).';
                elseif use_nctoolbox
                    var = double(squeeze(geo{'Pressure_reduced_to_MSL_msl'}(:)));
                end
                var = var*0.01; % Pa to db
                F = scatteredInterpolant(nlon(:),nlat(:),var(:));
                zz = F(lon_rho,lat_rho);
                zz(isnan(zz)) = 0;
                cff = squeeze(Pair(:,:,mm)).*(1-mask)+zz.*mask;
                Pair(:,:,mm) = cff;
            end
            if get_Qair
                if use_matlab
                    var = double(ncread(url,'Relative_humidity_height_above_ground')).';
                elseif use_nctoolbox
                    var = double(squeeze(geo{'Relative_humidity_height_above_ground'}(:)));
                end
                F = scatteredInterpolant(nlon(:),nlat(:),var(:));
                zz = F(lon_rho,lat_rho);
                zz(isnan(zz)) = 0;
                cff = squeeze(Qair(:,:,mm)).*(1-mask)+zz.*mask;
                Qair(:,:,mm) = cff;
            end
            if get_Wind
                if use_matlab
                    var = double(ncread(url,'u-component_of_wind_height_above_ground'));
                    var = squeeze(var(:,:,1)).';
                elseif use_nctoolbox
                    var = double(squeeze(geo{'u-component_of_wind_height_above_ground'}(:)));
                    var = squeeze(var(1,:,:));
                end
                F = scatteredInterpolant(nlon(:),nlat(:),var(:));
                cff = F(lon_rho,lat_rho);
                cff(isnan(cff)) = 0;
                Uwind_lamb = cff;
                
                if use_matlab
                    var = double(ncread(url,'v-component_of_wind_height_above_ground'));
                    var = squeeze(var(:,:,1)).';
                elseif use_nctoolbox
                    var = double(squeeze(geo{'v-component_of_wind_height_above_ground'}(:)));
                    var = squeeze(var(1,:,:));
                end
                F = scatteredInterpolant(nlon(:),nlat(:),var(:));
                cff = F(lon_rho,lat_rho);
                cff(isnan(cff)) = 0;
                Vwind_lamb = cff;
                
                %   Rotate winds to earth lon lat based on http://ruc.noaa.gov/RUC.faq.html
                %
                %   ROTCON_P          R  WIND ROTATION CONSTANT, = 1 FOR POLAR STEREO
                %                         AND SIN(LAT_TAN_P) FOR LAMBERT CONFORMAL
                %   LON_XX_P          R  MERIDIAN ALIGNED WITH CARTESIAN X-AXIS(DEG)
                %   LAT_TAN_P         R  LATITUDE AT LAMBERT CONFORMAL PROJECTION
                %                         IS TRUE (DEG)
                lat_tan_p = clat;                   % 25.0 for NAM;
                lon_xx_p = clon;                    % -95.0 for NAM;
                rotcon_p = sin(lat_tan_p*pi/180);
                deg2rad = 2*pi/360;
                
                angle2 = rotcon_p*(lon_rho-lon_xx_p)*deg2rad;
                sinx2 = sin(angle2);
                cosx2 = cos(angle2);
                Uwind_rot =  cosx2.*Uwind_lamb+sinx2.*Vwind_lamb;
                Vwind_rot = -sinx2.*Uwind_lamb+cosx2.*Vwind_lamb;
                
                % rotate winds to ROMS or user grid and merge with previous data
                cffx = Uwind_rot.*cos(angle_rho) + Vwind_rot.*sin(angle_rho);
                cffy = Vwind_rot.*cos(angle_rho) - Uwind_rot.*sin(angle_rho);
                Uwind(:,:,mm) = squeeze(Uwind(:,:,mm)).*(1-mask)+cffx.*mask;
                Vwind(:,:,mm) = squeeze(Vwind(:,:,mm)).*(1-mask)+cffy.*mask;
            end
        catch ME
            disp(['could not get that data at ', url])
        end
        if use_nctoolbox
          close(geo);
        end
    end
    save NAM_data.mat
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% write data to netcdf file

nc = ROMS_force_name;
if get_lwrad
  ncwrite(nc,'lrf_time',time);
  ncwrite(nc,'lwrad',lwrad);
  ncwrite(nc,'lwrad_down',lwrad_down);
end
if get_swrad
  ncwrite(nc,'srf_time',time);
  ncwrite(nc,'swrad',swrad);
end
if get_rain
  ncwrite(nc,'rain_time',time);
  rain(rain<0) = 0;
  ncwrite(nc,'rain',rain);
end
if get_Tair
  ncwrite(nc,'tair_time',time);
  Tair(Tair<-100) = 0;
  ncwrite(nc,'Tair',Tair);
end
if get_Pair
  ncwrite(nc,'pair_time',time);
  Pair(Pair<0) = 0;
  ncwrite(nc,'Pair',Pair);
end
if get_Qair
  ncwrite(nc,'qair_time',time);
  Qair(Qair<0) = 0;
  ncwrite(nc,'Qair',Qair);
end
if get_Wind
  ncwrite(nc,'wind_time',time);
  Uwind(Uwind<-100) = 0;
  Vwind(Vwind<-100) = 0;
  ncwrite(nc,'Uwind',Uwind);
  ncwrite(nc,'Vwind',Vwind);
end

disp(['------------ wrote ',ROMS_force_name,' ------------']);
