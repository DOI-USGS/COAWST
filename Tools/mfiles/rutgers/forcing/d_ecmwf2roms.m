%
% D_ECMWF2ROMS:  Driver script to create ROMS forcing NetCDF from ECMWF
%                ERA-Interim Dataset
%
% This is a user modifiable script showing how to create ROMS forcing
% NetCDF files. The data source is the ECMWF's ERA-Interim Dataset
%
% http://apps.ecmwf.int/datasets/data/interim_full_daily/
%
% This global data is available from Jan 1, 1979 to the present. It is
% usually extracted to a regional grid from the ECMWF data server.  No
% attempt is made to interpolate the forcing fields to a particular
% ROMS grid.  ROMS has the capability to perform the spatial and
% temporal interpolation of 2D forcing fields.
%
% This is a template script showing how to create several NetCDF files
% (one per each atmospheric forcing field) using ROMS metadata structure,
% which follows the schema of "nc_inq" or native 'ncinfo" function.
%
% The following parameters are used to extract ERA-Interim fields:
%
% Select time:   00:00:00     12:00:00
%
% Select step:   0  3  6  9  12
%

% svn $Id: d_ecmwf2roms.m 996 2020-01-10 04:28:56Z arango $
%=========================================================================%
%  Copyright (c) 2002-2020 The ROMS/TOMS Group      Hernan G. Arango      %
%    Licensed under a MIT/X style license           John Wilkin           %
%    See License_ROMS.txt                                                 %
%=========================================================================%
%
% Set file names. Input ERA-Interim input data file names.  The data
% was extracted for the Gulf of Mexico (100.5W, 10.5N) to (70.5W, 40.5N) 
% every 3-hours. The dataset was extracted in several annual files:
%
%                                               ERA Variables
%
%   ecmwf_era_atms_2000.nc              time, msl, tcc, v10u, v10u
%   ecmwf_era_flux_2000.nc              time, ewss, nsss, e, tp
%   ecmwf_era_heat_2000.nc              time, sshf, slhf, ssr, str
%   ecmwf_era_temp_2000.nc              time, v2t, v2d, ssrd, par
%   ...
%   ecmwf_era_atms_2012.nc
%   ecmwf_era_flux_2012.nc
%   ecmwf_era_heat_2012.nc
%   ecmwf_era_temp_2012.nc
%
% ERA-Interim variable names:  @ denotes accumulated quantity that must
%                              be divided by the time interval into the
%                              cycle 3, 6, 9 or 12 hours
%
% @  sshf    W m-2 s        surface sensible heat flux
% @  slhf    W m-2 s        surface latent heat flux
% @  ssr     W m-2 s        surface net solar radiation (shortwave)
% @  str     W m-2 s        surface net thermal radiation (longwave)
% @  strd    W m-2 s        surface thermal radiation downwards
% @  ewss    N m-2 s        east-west surface stress
% @  nsss    N m-2 s        north-south surface stress
% @  e       m              evaporation
% @  ro      m              runoff
% @  tcc     nondimensional total cloud cover [0:1]
% @  tp      m              total precipitation
% @  par     W m-2 s        photosynthetically active radiation at surface
%    msl     Pa             mean sea level pressure
%    v10v    m s-1          10 metre U wind component
%    vl0u    m s-1          10 metre V wind component
%    v2t     K              2 metre temperature
%    v2d     K              2 metre dewpoint temperature
%
% This dataset is written in compact way (short numbers). We need to
% convert to floating-point data and scale to ROMS units:
%
%   Uwind       (m s-1)         v10u
%   Vwind       (m s-1)         v10v
%   sustr       (N m-2)         ewss / (3*3600);   3-hour step
%   svstr       (N m-2)         nsss / (3*3600);   3-hour step
%   shflux      (W m-2)         (ssr+str+sshf+slhf) / (3*3600)
%   swrad       (W m-2)         ssr  / (3*3600);   3-hour step
%   lwrad_down  (W m-2)         strd / (3*3600);   3-hour step
%   latent      (W m-2)         slhf / (3*3600);   3-hour step
%   sensible    (W m-2)         sshf / (3*3600):   3-hour step
%   rain        (kg m-2 s-1)    tp * Rho_w / (3*3600)
%   evaporation (kg m-2 s-1)    e  * Rho_w / (3*3600)
%   swflux      (cm day-1)      (e - tp) * 100 / (3/24);  0.125 day step
%   cloud       (nondimesional) tcc
%   Pair        (mb)            msl / 100;   (1 mb = 100 Pa)
%   Tair        (Celsius)       t2m - 273.15;   (1 C = 273.15 K)
%   Qair        (percentage)    100 * (E/Es)
%
% where
%
%   Rho_w = 1000 kg m-3  (water density)
%
%   E  = 6.11 * 10.0 ^ (7.5 * d2m / (237.7 + d2m))    vapor pressure (mb)
%                                                     d2m in Celsius
%
%   Es = 6.11 * 10.0 ^ (7.5 * t2m / (237.7 + t2m))    saturation vapor
%                                                     pressure (mb)
%                                                     t2m in Celsius
%
%--------------------------------------------------------------------------

% Path to downloaded ERA data files.

Dir = fullfile(getenv('HOME'),                                          ...
               'ocean/repository/Projects/gom/Data/Forcing2');

% Base date for ROMS forcing files: "days since 2000-01-01 00:00:00".

mybasedate = datenum(2000,1,1,0,0,0);

% Set ROMS and ECMWF-ERA field NetCDF variable names for each forcing
% variable, input file name prefix, output file name. Notice that output
% forcing variables follows ROMS metadata design.
%
% If the ROMS option BULK_FLUXES is NOT activatived, we just need to
% process elements 1:5 in the structure vector (F).  Otherwise, we need
% process almost all the fields below.

clear S F

F( 1).Vname  = {'sustr', 'ewss'};
F( 1).Tname  = {'sms_time', 'time'};
F( 1).input  = 'ecmwf_era_flux_';
F( 1).output = 'gom_sms_era.nc';
F( 1).scale  = -1.0/(3*3600.0);

F( 2).Vname  = {'svstr', 'nsss'};
F( 2).Tname  = {'sms_time', 'time'};
F( 2).input  = 'ecmwf_era_flux_';
F( 2).output = 'gom_sms_era.nc';
F( 2).scale  = -1.0/(3*3600.0);

F( 3).Vname  = {'shflux', 'null'};
F( 3).Tname  = {'shf_time', 'time'};
F( 3).input  = 'ecmwf_era_heat_';
F( 3).output = 'gom_shflux_era.nc';
F( 3).scale  = -1.0/(3*3600.0);

F( 4).Vname  = {'swflux', 'null'};
F( 4).Tname  = {'swf_time', 'time'};
F( 4).input  = 'ecmwf_era_flux_';
F( 4).output = 'gom_swflux_era.nc';
F( 4).scale  = -100.0/(3*3600.0)*(24*3600.0);

F( 5).Vname  = {'swrad', 'ssr'};
F( 5).Tname  = {'srf_time',  'time'};
F( 5).input  = 'ecmwf_era_heat_';
F( 5).output = 'gom_swrad_era.nc';
F( 5).scale  = -1.0/(3*3600.0);

F( 6).Vname  = {'Uwind', 'u10'};
F( 6).Tname  = {'wind_time', 'time'};
F( 6).input  = 'ecmwf_era_atms_';
F( 6).output = 'gom_wind_era.nc';
F( 6).scale  = 1.0;

F( 7).Vname  = {'Vwind', 'v10'};
F( 7).Tname  = {'wind_time', 'time'};
F( 7).input  = 'ecmwf_era_atms_';
F( 7).output = 'gom_wind_era.nc';
F( 7).scale  = 1.0;

F( 8).Vname  = {'lwrad', 'str'};
F( 8).Tname  = {'lrf_time',  'time'};
F( 8).input  = 'ecmwf_era_heat_';
F( 8).output = 'gom_lwrad_era.nc';
F( 8).scale  = -1.0/(3*3600.0);

F( 9).Vname  = {'lwrad_down', 'strd'};
F( 9).Tname  = {'lrf_time',  'time'};
F( 9).input  = 'ecmwf_era_temp_';
F( 9).output = 'gom_lwrad_era.nc';
F( 9).scale  = -1.0/(3*3600.0);

F(10).Vname  = {'latent', 'slhf'};
F(10).Tname  = {'lhf_time',  'time'};
F(10).input  = 'ecmwf_era_heat_';
F(10).output = 'gom_latent_era.nc';
F(10).scale  = -1.0/(3*3600.0);

F(11).Vname  = {'sensible', 'sshf'};
F(11).Tname  = {'sen_time',  'time'};
F(11).input  = 'ecmwf_era_heat_';
F(11).output = 'gom_sensible_era.nc';
F(11).scale  = -1.0/(3*3600.0);

F(12).Vname  = {'cloud', 'tcc'};
F(12).Tname  = {'cloud_time',  'time'};
F(12).input  = 'ecmwf_era_atms_';
F(12).output = 'gom_cloud_era.nc';
F(12).scale  = 1.0;

F(13).Vname  = {'rain', 'tp'};
F(13).Tname  = {'rain_time',  'time'};
F(13).input  = 'ecmwf_era_flux_';
F(13).output = 'gom_rain_era.nc';
F(13).scale  = -1000.0/(3*3600.0);

F(14).Vname  = {'Pair', 'msl'};
F(14).Tname  = {'pair_time',  'time'};
F(14).input  = 'ecmwf_era_atms_';
F(14).output = 'gom_Pair_era.nc';
F(14).scale  = 0.01;

F(15).Vname  = {'Tair', 'v2t'};
F(15).Tname  = {'tair_time',  'time'};
F(15).input  = 'ecmwf_era_temp_';
F(15).output = 'gom_Tair_era.nc';
F(15).scale  = 1.0;

F(16).Vname  = {'Qair', 'v2d'};         % Use temperature (t2m) and
F(16).Tname  = {'qair_time',  'time'};  % dewpoint temperature (d2m)
F(16).input  = 'ecmwf_era_temp_';       % to compute relative humidity
F(16).output = 'gom_Qair_era.nc';
F(16).scale  = 1.0;

F(17).Vname  = {'PAR', 'par'};
F(17).Tname  = {'srf_time',  'time'};
F(17).input  = 'ecmwf_era_temp_';
F(17).output = 'gom_PAR_era.nc';
F(17).scale  = -1.0/(3*3600.0);

% Set field element in structure F to process.

doFields = [1 2 3 4 5];

% Various parameters.

spherical = true;                       % Spherical switch

LonMin    = -110.5;                     % GOM left-bottom longitude
LonMax    = -70.5;                      % GOM right-top   longitude
LatMin    = 10.5;                       % GOM left-bottom latitude
LatMax    = 40.5;                       % GOM right-top   latitude

FlipLat   = true;                       % Flip data in ERA NetCDF files
                                        % (see below information)

StrDay    = datenum('01-Jan-2000');     % starting day to process
EndDay    = datenum('31-Dec-2011');     % ending   day to process

nctype    = 'nc_float';                 % Input data is in single precision
Unlimited = true;                       % time dimension is umlimited in
                                        % output files

%--------------------------------------------------------------------------
% Create surface forcing NetCDF files: build creation parameters
% structure, S.  We want to create a single file for each forcing
% field.  However, the wind components "Uwind" and "Vwind" are
% saved in the same NetCDF file.
%--------------------------------------------------------------------------

Tindex       = [];
ReplaceValue = NaN;
PreserveType = true;

mode = netcdf.getConstant('CLOBBER');                    % overwrite!!!
mode = bitor(mode,netcdf.getConstant('64BIT_OFFSET'));

% The strategy here is to build manually the NetCDF metadata structure to
% facilitate creating several NetCDF file in a generic and compact way.
% This structure is similar to that returned by "nc_inq" or native Matlab
% function "ncinfo".
%
% Notice that we call the "roms_metadata" function to create the fields
% in the structure.  Then, we call "check_metadata" for fill unassigned
% values and to check for consistency.

OutFiles = {};

for n = doFields,

  ncname = char(F(n).output);
  Vname  = char(F(n).Vname{1});
  Tname  = char(F(n).Tname{1});
  
  if (n == 1),
    create_file = true;
  else
    create_file = ~any(strcmp(OutFiles, ncname));
  end
  OutFiles = [OutFiles, ncname];
  InpFile  = fullfile(Dir, strcat(char(F(n).input),'2000.nc'));
  
  if (create_file),
    disp(blanks(1));
    disp(['** Creating ROMS NetCDF forcing file: ', ncname,' **']);
    disp(blanks(1))
  end
  
  S.Filename = ncname;

  lon = nc_read(InpFile, 'longitude',                                   ...
                Tindex, ReplaceValue, PreserveType);

  lat = nc_read(InpFile, 'latitude',                                    ...
                Tindex, ReplaceValue, PreserveType);

% We want the longitude in the range from -180 to 180 instead of 0 to 360.
% Notice that the latitude is flip so the origin is at southwest corner
% of the extracted grid (LonMin,LatMin)

  ROMS_lon = lon; 
  ind = find(ROMS_lon > 180);
  if (~isempty(ind)),
    ROMS_lon(ind) = ROMS_lon(ind) - 360;
  end

  if (FlipLat),
    ROMS_lat = flipud(lat);
  else
    ROMS_lat = lat;
  end

  Im = length(ROMS_lon);
  Jm = length(ROMS_lat);

  S.Attributes(1).Name      = 'type';
  S.Attributes(1).Value     = 'FORCING file';

  S.Attributes(2).Name      = 'title';
  S.Attributes(2).Value     = ['ECMWF ERA-Interim Dataset, ',           ...
                               'Gulf of Mexico Region'];

  S.Attributes(3).Name      = 'history';
  S.Attributes(3).Value     = ['Forcing file created with ',            ...
                               which(mfilename) ' on ' date_stamp];

  S.Dimensions(1).Name      = 'lon';
  S.Dimensions(1).Length    = Im;
  S.Dimensions(1).Unlimited = false;

  S.Dimensions(2).Name      = 'lat';
  S.Dimensions(2).Length    = Jm;
  S.Dimensions(2).Unlimited = false;

  S.Dimensions(3).Name      = Tname;
  S.Dimensions(3).Length    = nc_constant('nc_unlimited');
  S.Dimensions(3).Unlimited = true;

  S.Variables(1) = roms_metadata('spherical');
  S.Variables(2) = roms_metadata('lon');
  S.Variables(3) = roms_metadata('lat');
  S.Variables(4) = roms_metadata(Tname, [], [], Unlimited);
  S.Variables(5) = roms_metadata(Vname, spherical, nctype, Unlimited);

% Edit the time variable 'units" attribute for the correct reference
% time and add calendar attribute.

  natts = length(S.Variables(4).Attributes);
  iatt  = strcmp({S.Variables(4).Attributes.Name}, 'units');
  S.Variables(4).Attributes(iatt).Value = ['days since ' datestr(mybasedate,31)];

  S.Variables(4).Attributes(natts+1).Name  = 'calendar';
  S.Variables(4).Attributes(natts+1).Value = 'gregorian';

% Change the dimensions of "shflux", "swflux", "sustr" or "svstr" to
% (lon,lat) dimensions to allow interpolation from a coarser grid in
% ROMS. As default, 'roms_metadata' uses ROMS native coordinates.
%
% This illustrates also how to manipulate the default metadata values set
% in 'roms_metadata'.  It has to be done before calling 'check_metadata'.

  if (strcmp(Vname, 'shflux') ||                                        ...
      strcmp(Vname, 'swflux') ||                                        ...
      strcmp(Vname, 'sustr')  ||                                        ...
      strcmp(Vname, 'svstr')),
    S.Variables(5).Dimensions(1).Name = 'lon';
    S.Variables(5).Dimensions(2).Name = 'lat';
  end
  
% Check ROMS metadata structure.  Fill unassigned fields.

  S = check_metadata(S);
  
% Create forcing NetCDF files.  Write our coordinates. Notice that
% "nc_append" is used here since we wand both wind components "Uwind"
% and "Vwind" to be in the same output NetCDF file for CF conventions.

  if (create_file),
    ncid = nc_create(S.Filename, mode, S); % create a new NetCDF file

    lon = repmat(ROMS_lon,  [1 Jm]);
    lat = repmat(ROMS_lat', [Im 1]);
  
    status = nc_write(S.Filename, 'spherical', int32(spherical));
    status = nc_write(S.Filename, 'lon',       lon);
    status = nc_write(S.Filename, 'lat',       lat);
  else
    nc_append(S.Filename, S)               % append to existing NetCDF file
  end

end

%--------------------------------------------------------------------------
% Convert forcing data to floating-point and scale to ROMS units.
%--------------------------------------------------------------------------

% This particular data set has a time coordinate in days starting
% on 1-Jan-1900.

epoch = datenum('01-Jan-1900');

% Set reference time for this application using the specified base date.

ref_time = (mybasedate - epoch);

StrYear  = str2num(datestr(StrDay, 'yyyy'));
EndYear  = str2num(datestr(EndDay, 'yyyy'));

MyRec = zeros([1 length(F)]);

disp(blanks(1));

for n = doFields,

  year  = StrYear;

  while (year <= EndYear)
    suffix  = strcat(num2str(year), '.nc');
    InpFile = fullfile(Dir, strcat(char(F(n).input),suffix));
    OutFile = char(F(n).output);
    Troms   = char(F(n).Tname{1});
    Tecmwf  = char(F(n).Tname{2});
    Vroms   = char(F(n).Vname{1});
    Vecmwf  = char(F(n).Vname{2});
    scale   = abs(F(n).scale);

% Determine time records to process.

    time = nc_read(InpFile, Tecmwf);
    time = time ./ 24;                              % hours to day

    StrRec = 1;
    EndRec = length(time);

% Read and write forcing fields. The forcing field are scaled to ROMS
% units.

    for Rec=StrRec:EndRec,

      MyRec(n) = MyRec(n) + 1;

      mydate = datestr(epoch+time(Rec));

      disp(['** Processing: ', Vroms, '  for  ',mydate,' **']);  

      frc_time = time(Rec) - ref_time;

      switch Vroms
        case 'Tair'
          field = nc_read(InpFile, Vecmwf, Rec);
          field = field - 273.15;                   % Kelvin to Celsius
        case 'Qair'     
          tsur  = nc_read(InpFile, 'v2t', Rec);     % 2m temperature
          tdew  = nc_read(InpFile, 'v2d', Rec);     % 2m dewpoint
          tsur  = tsur - 273.15;
          tdew  = tdew - 273.15;
          E     = 6.11 .* 10.0 .^ (7.5 .* tdew ./ (237.7 + tdew));
          Es    = 6.11 .* 10.0 .^ (7.5 .* tsur ./ (237.7 + tsur));
          field = 100.0 .* (E ./ Es);
        case 'swflux'     
          evap  = nc_read(InpFile, 'e' , Rec);      % evaporation
          prec  = nc_read(InpFile, 'tp', Rec);      % precipitation
          field = (evap - prec) .* scale;          
        case 'shflux'   
          sensbl = nc_read(InpFile, 'sshf' , Rec);  % sensible
          latent = nc_read(InpFile, 'slhf', Rec);   % latent
          nlwrad = nc_read(InpFile, 'str', Rec);    % net longwave
          nswrad = nc_read(InpFile, 'ssr', Rec);    % net shortwave
          field = (sensbl+latent+nlwrad+nswrad) .* scale;
        otherwise
          field = nc_read(InpFile, Vecmwf, Rec);
          field = field.*scale;
      end

      fieldfinal = field;

% If the scale F(n).scale is set to negative, the input ECMWF data is a
% cummulative integral in forecast cycle from hour zero.
 % For steps at 6, 9 and 12 hours we must separate last 3 hours of 
% integration from previous accumulation.
% At 3 hour step don't change anything

      if (F(n).scale < 0),
        step = rem(frc_time,0.5)*24;
        if step == 3
          fieldfinal = field;
        else
          fieldfinal = field - field_previous;  % At other steps subtract
        end                                     % the previous accumulation

        frc_time = frc_time - 1.5/24;	        % Center forcing time on the
                                                % accumulation interval

        field_previous = field;                 % Save this accumulation
                                                % to on the next step
      end

% If appropriate, flip the data so the origin (1,1) corresponds to
% the lower-left corner (LonMin,LatMin) of the extracted region.
% ERA data is written into NetCDF files with the origin at (1,end).
% That is, the origin is the left-upper corner of the extracted grid
% (LonMin,LatMax).  Users need to check and plot the extracted data
% from the ECMWF server.

      if (FlipLat),
        fieldfinal  = fliplr(fieldfinal);
      end

% Write out data.

      status = nc_write(OutFile, Troms, frc_time, MyRec(n));
      status = nc_write(OutFile, Vroms, fieldfinal, MyRec(n));
    end

%  Advance year counter.  Recall all input files are split in annual
%  files.

    year = year + 1;

  end
end
