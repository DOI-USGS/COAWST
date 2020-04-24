function [data]=load_sst_pfeg(GRDfile, StartDay, EndDay, varargin)

%
% LOAD_SST_PFEG:  Loads SST data for the region and time period
%
% [DATA]=load_sst_pfeg(GRDfile, StartDay, EndDay, sst_URL)
% 
%  Given a ROMS grid NetCDF, this function loads the SST data from the
%  extensive OpenDAP catalog maintained by NOAA PFEG OceanWatch. There
%  are several SST near real-time or composite blended products at
%  various spatial and temporal scales. The entire OceanWatch catalog
%  can be found at:  
%  
%  http://oceanwatch.pfeg.noaa.gov/thredds/catalog.html
%
% On Input:
%
%    GRDname       NetCDF grid file name (string)
%
%    StartDay      Starting period of interest (date number)
%
%                    Example:   startday=datenum(2004,1, 1)=731947
%
%    EndDay        Ending   period of interest (date number)
%
%                    Example:   startday=datenum(2004,1,15)=731961
%
%                    See Matlab intrinsic datenum(Y,Mo,D,H,Mi,S)
%                    for details.
%
%    sst_URL       SST dataset OpenDAP URL (OPTIONAL)
%
% On Output:
%
%    data          Composite SST data (structure array):
%
%                    Data.time     time of extracted data (date number)
%                    Data.lon      longitude of extracted data
%                    Data.lat      latitude  of extracted data
%                    Data.sst      sea surface temperatures
%
% You can use any of the following URL for "sst_URL":
%
% http://oceanwatch.pfeg.noaa.gov/thredds/dodsC/satellite/AA/ssta/1day
% http://oceanwatch.pfeg.noaa.gov/thredds/dodsC/satellite/AA/ssta/3day
% http://oceanwatch.pfeg.noaa.gov/thredds/dodsC/satellite/AA/ssta/5day
% http://oceanwatch.pfeg.noaa.gov/thredds/dodsC/satellite/AA/ssta/8day
% http://oceanwatch.pfeg.noaa.gov/thredds/dodsC/satellite/AA/ssta/14day
% http://oceanwatch.pfeg.noaa.gov/thredds/dodsC/satellite/AA/ssta/mday
%
% http://oceanwatch.pfeg.noaa.gov/thredds/dodsC/satellite/BA/ssta/5day
% http://oceanwatch.pfeg.noaa.gov/thredds/dodsC/satellite/BA/ssta/8day
% http://oceanwatch.pfeg.noaa.gov/thredds/dodsC/satellite/BA/ssta/mday
%
% Warning: This function either uses native Matlab NetCDF interface with
%          OpenDAP support (version 2012a or higher) or 'nc_varget' from
%          SNCTOOLS with OpenDAP to read data.
%

% svn $Id: load_sst_pfeg.m 996 2020-01-10 04:28:56Z arango $
%=========================================================================%
%  Copyright (c) 2002-2020 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license           Hernan G. Arango      %
%    See License_ROMS.txt                           John Wilkin           %
%=========================================================================%

%  Set optional arguments: Use the 5-day composite experimental product
%  as default.

sst_URL = 'http://oceanwatch.pfeg.noaa.gov/thredds/dodsC/satellite/BA/ssta/5day';

switch numel(varargin),
 case 1
   sst_URL = varargin{1};
end

%  Check arguments.

if (nargin < 3),
  error([' LOAD_SST_PFEG: You must specify a grid file along with',     ...
         ' starting and ending times']);
end

if (StartDay > EndDay),
  error([' LOAD_SST_PFEG: Your starting time must be greater than',     ...
         ' the ending time']);
end

data=[];

%--------------------------------------------------------------------------
%  Extract SST data for the period of interest.
%--------------------------------------------------------------------------

%  Choose NetCDF file interface (native or java from SNCTOOLS)

[method,~,~] = nc_interface(sst_URL);

%  Determine SST variables name;

Info = nc_vnames(sst_URL);

index = strfind({Info.Variables.Name}, 'sst');
index = ~cellfun(@isempty, index);

if (any(index)),
  sst_vname = Info.Variables(index).Name;
else
  error(' LOAD_SST_PFEG: unable to determine input SST variable name.');
end

index = strfind({Info.Variables.Name}, 'lon');
index = ~cellfun(@isempty, index);

if (any(index)),
  lon_vname = Info.Variables(index).Name;
else
  error(' LOAD_SST_PFEG: unable to determine input LON variable name.');
end

index = strfind({Info.Variables.Name}, 'lat');
index = ~cellfun(@isempty, index);

if (any(index)),
  lat_vname = Info.Variables(index).Name;
else
  error(' LOAD_SST_PFEG: unable to determine input LAT variable name.');
end

index = strfind({Info.Variables.Name}, 'time');
index = ~cellfun(@isempty, index);

if (any(index)),
  time_vname = Info.Variables(index).Name;
else
  error(' LOAD_SST_PFEG: unable to determine input TIME variable name.');
end

%  Find the time period of interest.

epoch = datenum([1970 1 1 0 0 0]);

switch(method),
  case {'native'}
    sst_time = epoch + ncread(sst_URL, time_vname)/86400;
  case {'java'}
    sst_time = epoch + nc_varget(sst_URL, time_vname)/86400;
end

T = find(sst_time >= StartDay & sst_time <= EndDay);

if (isempty(T)),
  disp(' LOAD_SST_PFEG: no data found for period of interest.');
  return;
end

switch(method),
  case {'java'}
    T = T - 1;            %  substract 1 because SNCTOOLS is 0-based.
end

%  Read SST longitudes and latitudes.

switch(method),
  case {'native'}
    sst_lon = ncread(sst_URL, lon_vname);
    sst_lat = ncread(sst_URL, lat_vname);
  case {'java'}
    sst_lon = nc_varget(sst_URL, lon_vname);
    sst_lat = nc_varget(sst_URL, lat_vname);
end

%  Read in application grid longitude and latitude.

switch(method),
  case {'native'}
    rlon = ncread(GRDfile, 'lon_rho');
    rlat = ncread(GRDfile, 'lat_rho');
  case {'java'}
    rlon = nc_varget(GRDfile, 'lon_rho');
    rlat = nc_varget(GRDfile, 'lat_rho');
end

MinLon = min(rlon(:))-0.5;
MaxLon = max(rlon(:))+0.5;
MinLat = min(rlat(:))-0.5;
MaxLat = max(rlat(:))+0.5;

%  Check how western longitudes are handled.

ind = find(sst_lon > MaxLon);
if (~isempty(ind)),
  sst_lon(ind) = sst_lon(ind) - 360;
end

ind = find(sst_lon < MinLon);
if (~isempty(ind)),
  sst_lon(ind) = sst_lon(ind) + 360;
end

%  Grab the indices for the application grid.

I = find(sst_lon >= MinLon & sst_lon <= MaxLon);
J = find(sst_lat >= MinLat & sst_lat <= MaxLat);

if (isempty(I) || isempty(J))
  disp(' LOAD_SST_PFEG: no data found for application grid.');
  return;
end

switch(method),
  case {'java'}
    I = I - 1;        %  substract 1 because SNCTOOLS is 0-based.
    J = J - 1;        %  substract 1 because SNCTOOLS is 0-based.
end

%  Read again coordinates for the selected region and time period to
%  be safe.

switch(method),
  case {'native'}
    data.time = ncread(sst_URL, time_vname, T(1), length(T));
    data.lon  = ncread(sst_URL, lon_vname,  I(1), length(I));
    data.lat  = ncread(sst_URL, lat_vname,  J(1), length(J));
  case {'java'}
    data.time = nc_varget(sst_URL, time_vname, T(1), length(T));
    data.lon  = nc_varget(sst_URL, lon_vname,  I(1), length(I));
    data.lat  = nc_varget(sst_URL, lat_vname,  J(1), length(J));
end

data.time = epoch + data.time./86400;

ind = find(data.lon > MaxLon);
if (~isempty(ind)),
  data.lon(ind) = data.lon(ind) - 360;
end

ind = find(data.lon < MinLon);
if (~isempty(ind)),
  data.lon(ind) = data.lon(ind) + 360;
end

%  Get the SST data (time,lat,lon). The data are actually 4D with second
%  coordinate being altitude.

switch(method),
  case {'native'}
    sst = ncread(sst_URL, sst_vname, [I(1) J(1) 1 T(1)],                ...
                [length(I) length(J) 1 length(T)]);
    data.sst = permute(squeeze(sst), [3 2 1]);
  case {'java'}
    sst = nc_varget(sst_URL, sst_vname, [I(1) J(1) 0 T(1)],             ...
                    [length(I) length(J) 1 length(T)]);   
    data.sst = permute(squeeze(sst), [3 2 1]);
end

return
