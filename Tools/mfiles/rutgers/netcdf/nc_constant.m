function value = nc_constant(parameter)

%
% NC_CONSTANT:  Gets numeric value of a named NetCDF library constant.
%
% value = nc_constant(parameter)
%
% This function returns the numeric value associated with the name of a
% NetCDF library parameter. Its is an umbrella on top the native function
% or older NetCDF interface in Matlab.  It needed to support older and
% newer version of Matlab.
%
% On Input:
%
%    parameter  NetCDF library parameter name (string)
%
% On Output:
%
%    value      NetCDF parameter numerical value
%
% Example:
%
%    value = nc_constant('NC_FLOAT')
% or
%    value = nc_constant('nc_float')
%

% svn $Id: nc_constant.m 711 2014-01-23 20:36:13Z arango $
%=========================================================================%
%  Copyright (c) 2002-2014 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%

% Determine NetCDF library constant numerical value from it name string.
  
if (~isempty(which('netcdf.open'))),
  value = netcdf.getConstant(parameter);             % native interface
else
  if (~isempty(which('mexnc'))),
    value = mexnc('parameter', parameter);           % MEXNC interface
  else
    switch(upper(parameter))
      case {'NC_GLOBAL'}
        value = -1;
      case {'NC_CLOBBER'}
        value = 0;
      case {'NC_UNLIMITED'}
        value = 0;
      case {'NC_FILL'}
        value = 0;
      case {'NC_NOWRITE'}
        value = 0;
      case {'NC_BYTE'}
        value = 1;
      case {'NC_NOWRITE'}
        value = 1;
      case {'NC_CHAR'}
        value = 2;
      case {'NC_SHORT'}
        value = 3;
      case {'NC_INT'}
        value = 4;
      case {'NC_NOCLOBBER'}
        value = 4;
      case {'NC_FORMAT_NETCDF4_CLASSIC'}
        value = 4;
      case {'NC_FLOAT'}
        value = 5;
      case {'NC_DOUBLE'}
        value = 6;
      case {'NC_NOFILL'}
        value = 256;
      case {'NC_64BIT_OFFSET'}
        value = 512;
      case {'NC_LOCK'}
        value = 1024;
      case {'NC_SHARE'}
        value = 2048;
      otherwise
        error('NC_CONSTANT: unable to determine NetCDF constant value');
    end
  end
end

return
