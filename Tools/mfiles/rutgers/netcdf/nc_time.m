function T = nc_time(ncfile, Rdate, Tname, varargin)

%
% NC_TIME:  Reads time variable NetCDF referenced to given date
%
% T = nc_time(ncfile, Rdate, Tname, Tindex)
%
% This function reads the time variable from specified NetCDF file,
% and computes the elapsed time since the specified reference date.
% If the time record (Tindex) is omitted, it process all the
% available records.
%
% On Input:
%
%    ncfile     NetCDF file name or URL file name (string)
%    Rdate      Reference date (seconds; scalar)
%    Tname      NetCDF time variable name (string)
%    Tindex     time record to process  (scalar, opional)
%
% On Output:
%
%    T          Time since reference date (original units)
%
% Examples:
%
% T = nc_time ('myfile.nc', datenum(1900,1,1)*86400, 'time')
% T = nc_time ('myurl', datenum('1-Jan-2000 12:00:00')*86400, 'time', 1)
%
  
% svn $Id: nc_time.m 996 2020-01-10 04:28:56Z arango $
%=========================================================================%
%  Copyright (c) 2002-2020 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%

T = [];

switch numel(varargin)
  case 0
    Tindex=[];
  case 1
    Tindex=varargin{1};
end

% Read time variable.

time = nc_read(ncfile, Tname, Tindex);

% Get time units attribute.

A = nc_getatt(ncfile, 'units', Tname);

% Covert to elapsed time since reference date.

ind = strfind(lower(A), 'since');

if (~isempty(ind))
  units = A(1:ind-2);
  Adate = A(ind+6:end);
  Tdate = datenum(Adate)*86400;
  if (Tdate ~= Rdate)
    switch (lower(units))
      case {'second', 'seconds'}
        T = (Tdate + time) - Rdate;
      case {'hour', 'hours'}
        T = (Tdate/3600 + time) - Rdate/3600;
      case {'day', 'days'}
        T = (Tdate/86400 + time) - Rdate/86400;
    end
  else
    T = time;
  end
else
  units = strtrim(lower(A));
  if (contains(units, 'julian day'))
    time = time - 2440000;
  end
  if (contains(units, 'second'))
    T = time + Rdate;
  elseif (contains(units, 'hour'))
    T = time + Rdate/3600;
  elseif (contains(units, 'day'))
    T = time + Rdate/86400;
  end
end

return
