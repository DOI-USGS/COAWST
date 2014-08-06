function [status]=nc_drename(ncfile, Dname_old, Dname_new)

%
% NC_DRENAME:  Changes the name of a dimension in a NetCDF file
%
% status = nc_drename(ncfile, Dname_old, Dname_new)
%
% This function renames requested dimension in a NetCDF file.
%
% On Input:
%
%    ncfile       NetCDF file name   (string)
%
%    Dname_old    Old dimension name (string)
%
%    Dname_new    New dimension name (string)
%
% On Output:
%
%    status       Error flag
%

% svn $Id: nc_drename.m 711 2014-01-23 20:36:13Z arango $
%=========================================================================%
%  Copyright (c) 2002-2014 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%

% Inquire about the contents of the NetCDF file.

Info = nc_inq(ncfile);

if (~any(strcmp({Info.Dimensions.Name}, Dname_old)))
  nc_inq(ncfile, true);
  disp(' ');
  error(['NC_DRENAME: cannot find NetCDF dimension: ',Vname_old]);
end

% Choose NetCDF file interface.

[method,~,~] = nc_interface(ncfile);

switch(method),
  case {'native'}
    status = nc_drename_matlab(ncfile, Dname_old, Dname_new, Info);
  case {'mexnc'}
    status = nc_drename_mexnc (ncfile, Dname_old, Dname_new, Info);
  case {'java'}
    error('NC_DRENAME: it is not possible to write into an OpenDAP file.');
  otherwise
    error('NC_DRENAME: unable to determine NetCDF processing interface.');
end

return

%--------------------------------------------------------------------------

function status = nc_drename_matlab(ncfile, Dname_old, Dname_new, Info)

%
% NC_DRENAME_MATLAB:  Changes the name of a dimension in a NetCDF file
%
% status = nc_drename_matlab(ncfile, Dname_old, Dname_new, Info)
%
% This function renames requested dimension in a NetCDF file using the
% native Matlab interface.
%
% On Input:
%
%    ncfile       NetCDF file name   (string)
%
%    Dname_old    Old dimension name (string)
%
%    Dname_new    New dimension name (string)
%
%    Info         NetCDF information structure (struct array)
%
% On Output:
%
%    status       Error flag
%

% Initialize.

status = 0;

% Make sure that the new variable name does not exist in the NetCDF
% file.

if (any(strcmp({Info.Dimensions.Name}, Dname_new)))
  nc_inq(ncfile, true);
  disp(' ');
  error(['NC_DRENAME_MATLAB: new name already exist: ', Dname_new]);
end

% Open NetCDF file and get variable ID.

ncid  = netcdf.open(ncfile, 'nc_write');
dimid = netcdf.inqDimID(ncid, Dname_old);

% Rename dimension.

netcdf.renameDim(ncid, dimid, Dname_new);

% Close NetCDF file.

netcdf.close(ncid);	
	
return	

%--------------------------------------------------------------------------

function status = nc_drename_mexnc(ncfile, Dname_old, Dname_new, Info)

%
% NC_DRENAME_MEXNC:  Changes the name of a dimension in a NetCDF file
%
% status = nc_drename(ncfile, Dname_old, Dname_new, Info)
%
% This function renames requested dimension in a NetCDF file using the
% MEXNC interface.
%
% On Input:
%
%    ncfile       NetCDF file name   (string)
%
%    Dname_old    Old dimension name (string)
%
%    Dname_new    New dimension name (string)
%
%    Info         NetCDF information structure (struct array)
%
% On Output:
%
%    status       Error flag
%

% Initialize error flag.

status = -1;

% Make sure that the new variable name does not exist in the NetCDF
% file.

if (any(strcmp({Info.Dimensions.Name}, Dname_new)))
  nc_inq(ncfile, true);
  disp(' ');
  error(['NC_DRENAME_MEXNC: new name already exist: ', Dname_new]);
end

% Open NetCDF file.

ncid = mexnc('open',ncfile,'nc_write');
if (ncid < 0),
  disp(' ');
  error(['NC_DRENAME_MEXNC: open - unable to open file: ', ncfile]);
end

% Get dimension ID.

dimid = mexnc('inq_dimid',ncid,Dname_old);
if (dimid < 0),
  [status]=mexnc('close',ncid);
  disp(' ');
  error(['NC_DRENAME_MEXNC: ncvarid - cannot find variable: ',Dname_old]);
end

% Put open file into define mode.

status = mexnc('redef',ncid);
if (status < 0),
  disp(' ');
  disp(mexnc('strerror',status));
  [status]=mexnc('close',ncid);
  error(['NC_DRENAME_MEXNC: redef - unable to put in definition ',      ...
         'mode: ',ncfile]);
end

% Change name of requested variable.

status = mexnc('rename_dim',ncid,dimid,Dname_new);
if (status < 0),
  disp(' ');
  disp(mexnc('strerror',status));
  error(['NC_DRENAME_MEXNC: rename_dim - unable to rename dimension: ', ...
         Dname_old]);
end

% Exit definition mode.

status = mexnc('enddef',ncid);
if (status < 0),
  disp(' ');
  disp(mexnc('strerror',status));
  error(['NC_DRENAME_MEXNC: enddef - unable to exit definition mode: ', ...
         ncfile]);
end

% Close NetCDF file.

cstatus = mexnc('ncclose',ncid);
if (cstatus < 0),
  disp(' ');
  disp(mexnc('strerror',status));
  error(['NC_DRENAME_MEXNC: ncclose - unable to close NetCDF file: ',   ...
         ncfile])
end

return
