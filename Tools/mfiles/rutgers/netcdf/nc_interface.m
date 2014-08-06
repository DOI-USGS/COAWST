function [method,url,ftype] = nc_interface(varargin)

%
% NC_INTERFACE:  Deternine the interface to process NetCDF files
%
% [method,url,ftype] = nc_interface(ncfile)
%
% This function determines the interface to use when processing a
% NetCDF file. The strategy is to always use the 'native' method
% when posible.
%
% On Input:
%
%    ncfile     NetCDF file name or URL name (string)
%
% On Output:
%
%    method     NetCDF interface for Matlab (string):
%
%                 'native'         Native Matlab interface
%                 'mexnc'          MEXNC interface
%                 'java'           SNCTOOLS Java interface
%
%    url        Switch indicating an OpenDAP file (logical)
%
%    ftype      File type (string)
%
%                 'NetCDF'         Classic NetCDF-3 file
%                 'NetCDF4'        NetCDF-4/HDF5 file
%                 'OpenDAP'        OpenDAP file
%
% Notes:
%
%   * The 'native' interface was introduced in Matlab version '2008b' for
%     NetCDF-3 type files. The NetCDF-4 support was introduced in version
%     '2010b'. The support for HDF5 files was completed in version '2011a'.
%     The OpenDAP support started in version '2012a'.
%
%   * The 'snctools' have support for NetCDF-3, NetCDF-4, HDF5, and a Java
%     interface for OpenDAP files.
%
%   * The 'mexcdf' is the oldest interface. It is no longer developed and
%     only has support for NetCDF-3 and NetCDF-4 classic files. It does not
%     have Java support for OpenDAP files.
%

% svn $Id: nc_interface.m 711 2014-01-23 20:36:13Z arango $
%=========================================================================%
%  Copyright (c) 2002-2014 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%

% Get Matlab version and check if the current Matlab version has the
% native NetCDF interface.

Mversion = version('-release');
Vyear    = sscanf(Mversion, '%i');
native   = ~isempty(which('netcdf.open'));

% If NetCDF file not provided, set the information available.

if (numel(varargin) == 0),
 if (native),
   method = 'native';
 else
   method = 'mexnc';
 end
 url = false;
 ftype = [];
 return
else
  ncfile = varargin{1};
end

% Check NetCDF type.  Adapted from SNCTOOLS.

if (~isempty(regexp(ncfile,'https*://', 'once'))),
  url = true;
else
  url = false;
end

if (~url),
  fid = fopen(ncfile,'r');
  if (ne(fid,-1)),
    signature = fread(fid,8,'uint8=>uint8');
    if (strcmp(char(signature(1:3))', 'CDF')),
      if (signature(4) == 1)
        ftype = 'NetCDF';                           % NetCDF classic
      elseif (signature(4) == 2)
        ftype = 'NetCDF';                           % NetCDF 64bit offset
      end         
    elseif (strcmp(char(signature(2:4))', 'HDF'))
      ftype = 'NetCDF4';                            % NetCDF4/HDF5
    end
  end
  fclose(fid);
else
  ftype = 'OpenDAP';                                % OpenDAP file
end

% Determine NetCDF interface to use. The 'native' interface is the
% preferred method. The native function 'ncinfo' is not available in
% 2010 versions.

method = 'native';

if (url && Vyear < 2012),
  method = 'java';
else
  if (~native || Vyear == 2010),
    method = 'mexnc';
  end
end

return
