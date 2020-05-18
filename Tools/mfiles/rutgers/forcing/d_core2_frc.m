%
% D_CORE2_FRC:  Driver script to create ROMS forcing NetCDF from CORE2
%               Global Dataset.
%
% This is a user modifiable script showing how to create ROMS forcing
% NetCDF files. The data source is the CORE 2 Global Air-Sea Flux Dataset.
% Notice that the original data set is sampled for the Gulf of Mexico (GOM)
% from 1 Jan 2000 to 31 Dec 2007.  No attempt is made to interpolate the
% forcing fields to a particular ROMS grid.  ROMS has the capability to
% perform the spatial and temporal interpolation of 2D forcing fields.
%
% This is a template script showing how to create several NetCDF files
% (one per each atmospheric forcing field) using ROMS metadata structure,
% which follows the schema of "nc_inq" or native 'ncinfo" function.
%

% svn $Id$
%=========================================================================%
%  Copyright (c) 2002-2020 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%
%
% Set file names. Input CORE2 global dataset NetCDF files are courtesy from
% Enrique Curchitser (Rutgers University).

Dir = 'http://alnilam.marine.rutgers.edu:8080/thredds/dodsC/Forcing/CORE2';

InpData = {fullfile(Dir, 'Uwind.1948-2007.nc'),                         ...
           fullfile(Dir, 'Vwind.1948-2007.nc'),                         ...
           fullfile(Dir, 'swrad.1948-2007.nc'),                         ...
           fullfile(Dir, 'lwrad.1948-2007.nc'),                         ...
           fullfile(Dir, 'rain.1948-2007.nc'),                          ...
           fullfile(Dir, 'Tair.1948-2007.nc'),                          ...
           fullfile(Dir, 'Qair.1948-2007.nc'),                          ...
           fullfile(Dir, 'Pair.1948-2007.nc')};

OutData = {'gom_wind_core2.nc',                                         ...
           'gom_wind_core2.nc',                                         ...
           'gom_swrad_core2.nc',                                        ...
           'gom_lwrad_core2.nc',                                        ...
           'gom_rain_core2.nc',                                         ...
           'gom_Tair_core2.nc',                                         ...
           'gom_Qair_core2.nc',                                         ...
           'gom_Pair_core2.nc'};

Nfiles = length(InpData);

% Set input and output time NetCDF variable names for each forcing
% variable.

InpTimeVar  = {'wind_time',                                             ...
               'wind_time',                                             ...
               'srf_time',                                              ...
               'lrf_time',                                              ...
               'rain_time',                                             ...
               'tair_time',                                             ...
               'qair_time',                                             ...
               'pair_time'};

OutTimeVar  = {'wind_time',                                             ...
               'wind_time',                                             ...
               'srf_time',                                              ...
               'lrf_time',                                              ...
               'rain_time',                                             ...
               'tair_time',                                             ...
               'qair_time',                                             ...
               'pair_time'};

% Set input and output field NetCDF variable names for each forcing
% variable. Notice that output variable follows ROMS metadata design.

InpFieldVar = {'Uwind',                                                 ...
               'Vwind',                                                 ...
               'swrad',                                                 ...
               'lwrad_down',                                            ...
               'rain',                                                  ...
               'Tair',                                                  ...
               'Qair',                                                  ...
               'Pair'};

OutFieldVar = {'Uwind',                                                 ...
               'Vwind',                                                 ...
               'swrad',                                                 ...
               'lwrad_down',                                            ...
               'rain',                                                  ...
               'Tair',                                                  ...
               'Qair',                                                  ...
               'Pair'};

% Various parameters.

spherical = true;                       % Spherical switch

LonMin    = -110;                       % GOM left-bottom longitude
LonMax    = -70;                        % GOM right-top   longitude
LatMin    = 10;                         % GOM left-bottom latitude
LatMax    = 40;                         % GOM right-top   latitude

StrDay    = datenum('01-Jan-2000');     % starting day to process
EndDay    = datenum('01-Dec-2008');     % ending   day to process

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

for n=1:Nfiles,

  disp(blanks(1));
  disp(['** Creating NetCDF file: ',char(OutData(n)),' **']);
  disp(blanks(1))
  
  S.Filename = char(OutData(n));

  lon = nc_read(char(InpData(n)), 'lon',                                ...
                Tindex, ReplaceValue, PreserveType);

  lat = nc_read(char(InpData(n)), 'lat',                                ...
                Tindex, ReplaceValue, PreserveType);

% Get data indices covering the GOM.  We want the longitude in a range
% from -180 to 180 instead of 0 to 360.

  lon_ind = find(lon >= mod(LonMin,360) & lon <= mod(LonMax,360));
  lat_ind = find(lat >= LatMin & lat <= LatMax);

  GOM_lon = lon(lon_ind); 
  GOM_lat = lat(lat_ind);

  ind = find(GOM_lon > 180);
  if (~isempty(ind)),
    GOM_lon(ind) = GOM_lon(ind) - 360;
  end

  Im = length(GOM_lon);
  Jm = length(GOM_lat);

  S.Attributes(1).Name      = 'type';
  S.Attributes(1).Value     = 'FORCING file';

  S.Attributes(2).Name      = 'title';
  S.Attributes(2).Value     = ['CORE 2 Global Air-Sea Flux ',           ...
                               'Dataset, Gulf of Mexico Region'];

  S.Attributes(3).Name      = 'history';
  S.Attributes(3).Value     = ['Forcing file created from ',            ...
                               'gom_core2: ', date_stamp];

  S.Dimensions(1).Name      = 'lon';
  S.Dimensions(1).Length    = Im;
  S.Dimensions(1).Unlimited = false;

  S.Dimensions(2).Name      = 'lat';
  S.Dimensions(2).Length    = Jm;
  S.Dimensions(2).Unlimited = false;

  S.Dimensions(3).Name      = char(OutTimeVar(n));
  S.Dimensions(3).Length    = nc_constant('nc_unlimited');
  S.Dimensions(3).Unlimited = true;

  S.Variables(1) = roms_metadata('spherical');
  S.Variables(2) = roms_metadata('lon');
  S.Variables(3) = roms_metadata('lat');
  S.Variables(4) = roms_metadata(char(OutTimeVar(n)), [], [], Unlimited);
  S.Variables(5) = roms_metadata(char(OutFieldVar(n)),                  ...
                                 spherical, nctype, Unlimited);

% Edit the time variable 'units" attribute for the correct reference
% time and add calendar attribute.

  natts = length(S.Variables(4).Attributes);  
  iatt  = strcmp({S.Variables(4).Attributes.Name}, 'units');
  S.Variables(4).Attributes(iatt).Value = 'days since 2000-01-01 00:00:00';

  S.Variables(4).Attributes(natts+1).Name  = 'calendar';
  S.Variables(4).Attributes(natts+1).Value = 'gregorian';

% Check ROMS metadata structure.  Fill unassigned fields.

  S = check_metadata(S);
  
% Create forcing NetCDF files.  Write our coordinates. Notice that
% "nc_append" is used here since we wand both wind components "Uwind"
% and "Vwind" to be in the same output NetCDF file for CF conventions.

  if (strcmp(OutFieldVar(n), 'Vwind')),
    nc_append(S.Filename, S)               % append to existing NetCDF file
  else
    ncid = nc_create(S.Filename, mode, S); % create a new NetCDF file

    lon = repmat(GOM_lon,  [1 Jm]);
    lat = repmat(GOM_lat', [Im 1]); 
  
    status = nc_write(S.Filename, 'spherical', int32(spherical));
    status = nc_write(S.Filename, 'lon',       lon);
    status = nc_write(S.Filename, 'lat',       lat);
  end

end

%--------------------------------------------------------------------------
% Extract and write desired time records.
%--------------------------------------------------------------------------

% This particular data set has a time coordinate in days starting
% on 1-Jan-1900.

epoch = datenum('01-Jan-1900');

% Set reference time for this application such that the output time units
% are "days since 2000-01-01 00:00:00".

ref_time = (datenum('01-Jan-2000')-datenum('01-Jan-1900'));

for n=1:Nfiles,

  InpFile = char(InpData(n));
  OutFile = char(OutData(n));

% Determine longitude and latitude indices to extract.  Notice that
% not spatial or temporal interpolation is carried out. 

  lon = nc_read(char(InpData(n)), 'lon',                                ...
                Tindex, ReplaceValue, PreserveType);
  lat = nc_read(char(InpData(n)), 'lat',                                ...
                Tindex, ReplaceValue, PreserveType);

  lon_ind = find(lon >= mod(LonMin,360) & lon <= mod(LonMax,360));
  lat_ind = find(lat >= LatMin & lat <= LatMax);

  if (~isempty(lon_ind)),
    Imin=min(lon_ind);
    Imax=max(lon_ind);
  else
    error('GOM_CORE2: longitude range to extract is empty.');
  end

  if (~isempty(lat_ind)),
    Jmin=min(lat_ind);
    Jmax=max(lat_ind);
  else
    error('GOM_CORE2: latitude range to extract is empty.');
  end

% Determine dataset time record to process.

  time = nc_read(InpFile, char(InpTimeVar(n)));

  ind = find(time >= StrDay-epoch & time <= EndDay-epoch);
  if (~isempty(ind)),
    StrRec = ind(1)-1;     % substract 1 to bound initialization forcing
    EndRec = ind(end);
  else
    error('GOM_CORE2: ''StrDay'' and ''EndDay'' not available in dataset');
  end

% Extract and write forcing fields.  Notice that the final time string
% can be computed as:
%
% datestr(datenum('01-Jan-2000') + frc_time).

  MyRec = 0;

  for Rec=StrRec:EndRec,

    MyRec = MyRec + 1;

    mydate = datestr(epoch+time(Rec));

    disp(blanks(1))
    disp(['** Processing: ', mydate,' **']);  

    frc_time = time(Rec) - ref_time;

    field = nc_read(InpFile, char(InpFieldVar(n)), Rec);
    F     = field(Imin:Imax, Jmin:Jmax);

    status = nc_write(OutFile, char(OutTimeVar(n)), frc_time, MyRec);
    status = nc_write(OutFile, char(OutFieldVar(n)), F, MyRec);
  end

end
