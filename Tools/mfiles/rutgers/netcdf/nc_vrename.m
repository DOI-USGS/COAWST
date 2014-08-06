function status = nc_vrename(Fname, Vname_old, Vname_new)

%
% NC_VRENAME:  Changes the name of a variable in a NetCDF file
%
% status = nc_vrename(Fname, Vname_old, Vname_new)
%
% This function renames requested variable in a NetCDF file.
%
% On Input:
%
%    Fname        NetCDF file name  (string)
%
%    Vname_old    Old variable name (string)
%
%    Vname_new    New variable name (string)
%
% On Output:
%
%    status       Error flag
%

% svn $Id: nc_vrename.m 711 2014-01-23 20:36:13Z arango $
%=========================================================================%
%  Copyright (c) 2002-2014 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%

% Inquire about the contents of the NetCDF file.

Info = nc_inq(Fname);

if (~any(strcmp({Info.Variables.Name}, Vname_old)))
  nc_inq(Fname, true);
  disp(' ');
  error(['NC_VRENAME: cannot find NetCDF variable: ',Vname_old]);
end

% Choose NetCDF file interface.

[method,~,~] = nc_interface(Fname);

switch(method),
  case {'native'}
    status = nc_vrename_matlab(Fname, Vname_old, Vname_new, Info);
  case {'mexnc'}
    status = nc_vrename_mexnc (Fname, Vname_old, Vname_new, Info);
  case {'java'}
    error('NC_VRENAME: it is not possible to write into an OpenDAP file.');
  otherwise
    error('NC_VRENAME: unable to determine NetCDF processing interface.');
end

return

%--------------------------------------------------------------------------

function status = nc_vrename_matlab(Fname, Vname_old, Vname_new, Info)
  
%
% NC_VRENAME_MATLAB:  Changes the name of a variable in a NetCDF file
%
% status = nc_vrename_matlab(Fname, Vname_old, Vname_new, Info)
%
% This function renames requested variable in a NetCDF file using the
% native Matlab interface.
%
% On Input:
%
%    Fname        NetCDF file name  (string)
%
%    Vname_old    Old variable name (string)
%
%    Vname_new    New variable name (string)
%
%    Info         NetCDF information structure (struct array)
%
% On Output:
%
%    status       error flag
%

% Initialize.

status = 0;

% Make sure that the new variable name does not exist in the NetCDF
% file.

if (any(strcmp({Info.Variables.Name}, Vname_new)))
  nc_inq(Fname, true);
  disp(' ');
  error(['NC_VRENAME_MATLAB: new name already exist: ', Vname_new]);
end

% Open NetCDF file and get variable ID.

ncid  = netcdf.open(Fname, 'nc_write');
varid = netcdf.inqVarID(ncid, Vname_old);

% Rename variable.

netcdf.renameVar(ncid, varid, Vname_new);
	
% Close NetCDF file.

netcdf.close(ncid);	
	
return

%--------------------------------------------------------------------------

function status = nc_vrename_mexnc(Fname, Vname_old, Vname_new, Info)
  
%
% NC_VRENAME_MEXNC:  Changes the name of a variable in a NetCDF file
%
% status = nc_vrename_mexnc(Fname, Vname_old, Vname_new, Info)
%
% This function renames requested variable in a NetCDF file.
%
% On Input:
%
%    Fname        NetCDF file name  (string)
%
%    Vname_old    Old variable name (string)
%
%    Vname_new    New variable name (string)
%
%    Info         NetCDF information structure (struct array)
%
% On Output:
%
%    status     error flag
%

% Initialize error flag.

status = -1;

% Make sure that the new variable name does not exist in the NetCDF
% file.

if (any(strcmp({Info.Variables.Name}, Vname_new)))
  nc_inq(Fname, true);
  disp(' ');
  error(['NC_VRENAME_MEXNC: new name already exist: ', Vname_new]);
end

% Open NetCDF file.

ncid = mexnc('open',Fname,'nc_write');
if (ncid < 0),
  disp(' ');
  error(['NC_VRENAME_MEXNC: open - unable to open file: ', Fname]);
end

% Get variable ID.

varid = mexnc('ncvarid',ncid,Vname_old);
if (varid < 0),
  status = mexnc('close',ncid);
  disp(' ');
  error(['NC_VRENAME_MEXNC: ncvarid - cannot find variable: ',Vname_old]);
end

% Put open file into define mode.

status = mexnc('redef',ncid);
if (status < 0),
  disp(' ');
  disp(mexnc('strerror',status));
  status = mexnc('close',ncid);
  error(['NC_VRENAME_MEXNC: redef - unable to put in definition ',      ...
         'mode: ', Fname]);
end

% Change name of requested variable.

status = mexnc('rename_var',ncid,varid,Vname_new);
if (status < 0),
  disp(' ');
  disp(mexnc('strerror',status));
  error(['NC_VRENAME_MEXNC: rename_var - unable to rename ',            ...
         'variable: ',Vname_old]);
end

% Exit definition mode.

status = mexnc('enddef',ncid);
if (status < 0),
  disp(' ');
  disp(mexnc('strerror',status));
  error(['NC_VRENAME_MEXNC: enddef - unable to exit definition ',       ...
         'mode: ',Fname]);
end

% Close NetCDF file.

cstatus = mexnc('ncclose',ncid);
if (cstatus < 0),
  disp(' ');
  disp(mexnc('strerror',status));
  error(['NC_VRENAME_MEXNC: ncclose - unable to close NetCDF file: ',   ...
	 Fname]);
end

return
