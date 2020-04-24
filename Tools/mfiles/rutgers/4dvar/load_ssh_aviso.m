function [D]=load_ssh_aviso(GRDfile, StartDay, EndDay, ssh_URL)

%
% LOAD_SSH_AVISO:  Loads AVISO data for the region and time period
%
% [D]=load_ssh_aviso(GRDfile, StartDay, EndDay, ssh_URL)
% 
%  Given a ROMS grid NetCDF, this function loads the AVISO sea level
%  anomaly covering the region for the time period specified (StartDay
%  to EndDay).
%   
% On Input:
%
%    GRDname      NetCDF grid file name (string)
%
%    StartDay     Starting period of interest (date number)
%
%                   Example:   StartDay=datenum(2004,1,1)=731947
%
%    EndDay       Ending   period of interest (date number)
%
%                   Example:   StartDay=datenum(2004,2,1)=731978
%
%                   See Matlab intrinsic datenum(Y,Mo,D,H,Mi,S)
%                   for details.
%
%    ssh_URL      SSH dataset OpenDAP URL (string). The AVISO dataset
%                   is password proctected and users need to request
%                   access by filling the following form from:
%                   http://www.aviso.oceanobs.com/en/data/registration-form
%                
% On Output:
%
%    D            AVISO sea level data (structure array):
%
%                   D.ssh      sea level anomalies, cm (time,lat,lon)
%                   D.time     time of analysis (date number)
%                   D.lon      longitude of extracted data (lat,lon)
%                   D.lat      latitude  of extracted data (lat,lon)
%
% There are various options for "ssh_URL". Check the following link for a
% summary of all sea level anomaly products:
%
%  http://www.aviso.oceanobs.com/en/data/data-access-services/opendap/opendap-sla-products.html
%  
% For example, to use Ssalto/Duacs merged near-real time sea level anomaly
% data (merged NRT-MSLA) set "ssh_URL" to:
%
% Dir = 'http://USERNAME:PASSWORD@opendap.aviso.oceanobs.com/thredds/dodsC'
%
% ssh_URL = strcat(Dir, '/dataset-duacs-nrt-over30d-global-merged-msla-h')
%
% Other possibilities are:
%
% ssh_URL = strcat(Dir, '/dataset-duacs-dt-ref-global-merged-msla-h-daily')
% ssh_URL = strcat(Dir, '/dataset-duacs-dt-upd-global-merged-msla-h-daily')
%
% Warning: This function either uses native Matlab NetCDF interface with
%          OpenDAP support (version 2012a or higher) or 'nc_varget' from
%          SNCTOOLS with OpenDAP to read data.
%

% svn $Id: load_ssh_aviso.m 996 2020-01-10 04:28:56Z arango $
%=========================================================================%
%  Copyright (c) 2002-2020 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license           Hernan G. Arango      %
%    See License_ROMS.txt                           Brian Powell          %
%=========================================================================%

%  Check arguments.

if (nargin < 4),
  error([' LOAD_SSH_DATA: You must specify a grid file along with',     ...
         ' starting and ending times and OpenDAP URL']);
end,

if (StartDay > EndDay),
  error([' LOAD_SSH_DATA: Your starting time must be greater than',     ...
         ' the ending time']);
end,

D=[];

%--------------------------------------------------------------------------
%  Extract AVISO data for the period of interest.
%--------------------------------------------------------------------------

%  Choose NetCDF file interface (native or java from SNCTOOLS)

[method,~,~] = nc_interface(ssh_URL);

%  Determine SST variables name;

Info = nc_vnames(ssh_URL);

index = strfind({Info.Variables.Name}, 'Grid');
index = ~cellfun(@isempty, index);

if (any(index)),
  ssh_vname = Info.Variables(index).Name;
else
  error(' LOAD_SSH_AVISO: unable to determine input SSH variable name.');
end

index = strfind({Info.Variables.Name}, 'Longitude');
index = ~cellfun(@isempty, index);

if (any(index)),
  lon_vname = Info.Variables(index).Name;
else
  error([' LOAD_SSH_AVISO: unable to determine input longitude',        ...
         ' variable name.']);
end

index = strfind({Info.Variables.Name}, 'Latitude');
index = ~cellfun(@isempty, index);

if (any(index)),
  lat_vname = Info.Variables(index).Name;
else
  error([' LOAD_SSH_AVISO: unable to determine input latitude',         ...
         ' variable name.']);
end

index = strfind({Info.Variables.Name}, 'time');
index = ~cellfun(@isempty, index);

if (any(index)),
  time_vname = Info.Variables(index).Name;
else
  error(' LOAD_SST_PFEG: unable to determine input time variable name.');
end

%  Find the time period of interest.  The AVISO data has a reference time
%  is hours 1-Jan-1950.

epoch = datenum(1950,1,1);

switch(method),
  case {'native'}
    ssh_time = epoch + double(ncread(ssh_URL, time_vname)) ./ 24;
  case {'java'}
    ssh_time = epoch + nc_varget(ssh_URL, time_vname) ./ 24;
end

T = find(ssh_time >= StartDay & ssh_time <= EndDay);

if (isempty(T)),
  disp(' LOAD_SSH_AVISO: no data found for rquested period.');
  disp([blanks(17), 'Available Times: ', datestr(ssh_time(1))]);
  disp([blanks(34), datestr(ssh_time(end))]);
  return;
end

D.time = ssh_time(T);

%  Find the region of interest.

switch(method),
  case {'native'}
    ssh_lon = ncread(ssh_URL, lon_vname);
    ssh_lat = ncread(ssh_URL, lat_vname);
  case {'java'}
    ssh_lon = nc_varget(ssh_URL, lon_vname);
    ssh_lat = nc_varget(ssh_URL, lat_vname);
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

ind = find(ssh_lon > MaxLon);
if (~isempty(ind)),
  ssh_lon(ind) = ssh_lon(ind) - 360;
end

ind = find(ssh_lon < MinLon);
if (~isempty(ind)),
  ssh_lon(ind) = ssh_lon(ind) + 360;
end

% Grab the indices for the application grid.

I = find( ssh_lon >= MinLon & ssh_lon <= MaxLon );
J = find( ssh_lat >= MinLat & ssh_lat <= MaxLat );

if (isempty(I) || isempty(J))
  disp(' LOAD_SST_AVISO: no data found for application grid.');
  return;
end

[D.lon, D.lat] = meshgrid(ssh_lon(I), ssh_lat(J));

% Get the SSH data.

switch(method),
  case {'native'}
    ssh = ncread(ssh_URL, ssh_vname, [J(1) I(1) T(1)],                  ...
                 [length(J) length(I) length(T)]);

    D.ssh = permute(ssh, [3 1 2]);
  case {'java'}
    T = T - 1;        %  substract 1 because SNCTOOLS is 0-based.
    I = I - 1;        %  substract 1 because SNCTOOLS is 0-based.
    J = J - 1;        %  substract 1 because SNCTOOLS is 0-based.

    ssh = nc_varget(ssh_URL, ssh_vname, [J(1) I(1) T(1)],               ...
                 [length(J) length(I) length(T)]);

    D.ssh = permute(ssh, [3 1 2]);
end

return
