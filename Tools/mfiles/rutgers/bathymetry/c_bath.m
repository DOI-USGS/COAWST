function [status]=c_bath(Im, Jm, Bname)

%
% C_BATH:  Create a topography NetCDF file
%
% [status]=c_bath(Im,Jm,Bname)
%
% This function creates a topography NetCDF file.
%
% On Input:
%
%    Im          Number of RHO-points in the X-direction.
%    Jm          Number of RHO-points in the Y-direction.
%    Bname       Topography NetCDF file name.
%
% On Output:
%
%    status      Error flag.
%

% svn $Id$
%=========================================================================%
%  Copyright (c) 2002-2025 The ROMS Group                                 %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.md                            Hernan G. Arango      %
%=========================================================================%

%--------------------------------------------------------------------------
%  Inquire dimensions from a existing NeTCDF file.
%--------------------------------------------------------------------------

Dname = {'lon', 'lat', 'bath'};
	 
Dsize.lon =Im;
Dsize.lat =Jm;
Dsize.bath=nc_constant('nc_unlimited');

Vname.lon ='lon';
Vname.lat ='lat';
Vname.bath='hraw';

%--------------------------------------------------------------------------
%  Create topography NetCDF file.
%--------------------------------------------------------------------------

mode = netcdf.getConstant('CLOBBER');
mode = bitor(mode, netcdf.getConstant('64BIT_OFFSET'));
ncid = netcdf.create(Bname, mode);

%--------------------------------------------------------------------------
%  Create global attribute(s).
%--------------------------------------------------------------------------

varid = netcdf.getConstant('nc_global');

type = 'GRID file';
netcdf.putAtt(ncid, varid, 'type', type);

history = ['GRID file using Matlab script: c_bath, ', date_stamp];
netcdf.putAtt(ncid, varid, 'history', history);

%--------------------------------------------------------------------------
%  Define dimensions.
%--------------------------------------------------------------------------

for value = Dname
  field = char(value);
  did.(field) = netcdf.defDim(ncid, field, Dsize.(field));
end

%--------------------------------------------------------------------------
%  Define variables.
%--------------------------------------------------------------------------

% Define spherical switch.

Var.name          = 'spherical';
Var.type          = nc_constant('nc_int');
Var.dimid         = [];
Var.long_name     = 'grid type logical switch';
Var.flag_values   = [0 1];
Var.flag_meanings = ['Cartesian', blanks(1), ...
                     'spherical'];
[~,status]=nc_vdef(ncid,Var);
if (status ~= 0), return, end,
clear Var

%  Longitude.

Var.name          = Vname.lon;
Var.type          = nc_constant('nc_double');
Var.dimid         = [did.lat did.lon];
Var.long_name     = 'longitude';
Var.units         = 'degree_east';
Var.standard_name = 'longitude';
[~,status]=nc_vdef(ncid,Var);
if (status ~= 0), return, end,
clear Var

%  Latitude.

Var.name          = Vname.lat;
Var.type          = nc_constant('nc_double');
Var.dimid         = [did.lat did.lon];
Var.long_name     = 'latitute';
Var.units         = 'degree_north';
Var.standard_name = 'latitude';
[~,status]=nc_vdef(ncid,Var);
if (status ~= 0), return, end,
clear Var

%  Topography.

Var.name          = Vname.bath;
Var.type          = nc_constant('nc_double');
Var.dimid         = [did.bath did.lat did.lon];
Var.long_name     = 'topography';
Var.units         = 'meter';
Var.coordinates   = strcat([Vname.lon,' ',Vname.lat]);
[~,status]=nc_vdef(ncid,Var);
if (status ~= 0), return, end,
clear Var

%---------------------------------------------------------------------------
%  Leave definition mode and close NetCDF file.
%---------------------------------------------------------------------------

netcdf.endDef(ncid);
netcdf.close(ncid);

return

