function D = nc_dinfo(fname)

%
% NC_DINFO:  Inquire about the dimensions in a NetCDF file
%
% D = nc_dinfo(fname)
%
% This function gets all dimension information about the requested NetCDF
% file.
%
% On Input:
%
%    fname      NetCDF file name or URL file name (string)
%
% On Output:
%
%    D          Dimensions information (Struct array):
%
%                 D(:).Name       dimension name (string)
%                 D(:).Length     dimension length (double)
%                 D(:).Unlimited  unlimited dimension (logical)
%                 D(:).dimid      NetCDF dimension id (double)
%

% svn $Id: nc_dinfo.m 711 2014-01-23 20:36:13Z arango $
%=========================================================================%
%  Copyright (c) 2002-2014 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%

% Choose NetCDF file interface.

[method,~,~] = nc_interface(fname);

switch(method),
  case {'native'}
    D = nc_dinfo_matlab(fname);                 % Matlab native interface
  case {'java'}
    D = nc_dinfo_java  (fname);                 % SNCTOOLS JAVA interface
  case {'mexnc'}
    D = nc_dinfo_mexnc (fname);                 % MEXNC inteface
  otherwise
    error(['NC_DINFO: unable to determine NetCDF processing interface']);
end

return

%--------------------------------------------------------------------------

function [D]=nc_dinfo_matlab(fname)

%
% NC_DINFO_MATLAB:  Inquire about the dimensions in a NetCDF file
%
% [dnames,dsizes,recdim]=nc_dinfo_matlab(fname)
%
% This function gets dimensions information about requested NetCDF
% file (regular or OpenDAP URL). It uses native matlab function
% "ncinfo" that it is available in version 2012a or higher.
%

% Initialize output structure.

D = struct('Name', '',                                                  ...
           'Length', [],                                                ...
           'Unlimited', [],                                             ...
	   'dimid', []);

% Inquire information from URL NetCDF file.

Info = ncinfo(fname); 

% Extract dimension information. Notice thar in NetCDF the ID count
% start from zero.

ndims=length(Info.Dimensions);

for n=1:ndims,
  D(n).Name      = Info.Dimensions(n).Name;
  D(n).Length    = Info.Dimensions(n).Length;
  D(n).Unlimited = Info.Dimensions(n).Unlimited;
  D(n).dimid     = n-1;
end

return

%--------------------------------------------------------------------------

function D = nc_dinfo_java(fname)

%
% NC_DINFO_JAVA:  Inquire about the dimensions in a NetCDF file
%
% D = nc_dinfo_java(fname)
%
% This function gets dimensions information about requested URL
% OpenDAP file(s). It uses SNCTOOLS function "nc_info".
%

% Initialize output structure.

D = struct('Name', '',                                                  ...
           'Length', [],                                                ...
           'Unlimited', [],                                             ...
	   'dimid', []);

% Inquire information from URL NetCDF file.

Info=nc_info(fname); 

% Extract dimension information. Notice thar in NetCDF the ID count
% start from zero.

ndims=length(Info.Dimension);

for n=1:ndims,
  D(n).Name      = Info.Dimension(n).Name;
  D(n).Length    = Info.Dimension(n).Length;
  D(n).Unlimited = Info.Dimension(n).Unlimited;
  D(n).dimid     = n-1;
end

return

%--------------------------------------------------------------------------

function D = nc_dinfo_mexnc(fname)

%
% NC_DINFO_MEXNC:  Inquire about the dimensions in a NetCDF file
%
% D = nc_dinfo_mexnc(fname)
%
% This function gets dimensions information about requested NetCDF
% file. It uses MEXNC functions. Therefore, it cannot process a URL
% OpenDAP file.
%

% Initialize output structure.

D = struct('Name', '',                                                  ...
           'Length', [],                                                ...
           'Unlimited', [],                                             ...
	   'dimid', []);

% Open NetCDF file.

[ncid,status]=mexnc('open',fname,'nc_nowrite');
if (status ~= 0),
  disp('  ');
  disp(mexnc('strerror',status));
  error(['NC_DINFO: ncopen - unable to open file: ', fname]);
end
 
% Inquire about contents.

[ndims,~,~,recdim,status]=mexnc('inq',ncid);
if (status ~= 0),
  disp('  ');
  disp(mexnc('strerror',status));
  error(['NC_DINFO: INQ - cannot inquire file: ',fname]);
end

% Extract dimension information. Notice thar in NetCDF the ID count
% start from zero.

D.Filename = fname;

for n=1:ndims;
  dimid=n-1;
  [name,size,status]=mexnc('inq_dim',ncid,dimid);
  if (status ~= 0),
    disp('  ');
    disp(mexnc('strerror',status));
    error(['NC_DINFO: INQ_DIM - unable to inquire about dimension ',  ...
           'ID: ',num2str(n)]);
  else
    D(n).Name   = name;
    D(n).Length = size;
    if (dimid == recdim),
      D(n).Unlimited = true;
    else
      D(n).Unlimited = false;
    end
    D(n).dimid = dimid;
  end
end

% Close NetCDF file.

[status]=mexnc('close',ncid);
if (status ~= 0),
  disp('  ');
  disp(mexnc('strerror',status));
  error(['NC_DINFO: CLOSE - unable to close file: ', fname]);
end

return
