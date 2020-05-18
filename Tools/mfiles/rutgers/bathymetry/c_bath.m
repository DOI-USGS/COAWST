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

% svn $Id: c_bath.m 996 2020-01-10 04:28:56Z arango $
%=========================================================================%
%  Copyright (c) 2002-2020 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%

%--------------------------------------------------------------------------
%  Inquire dimensions from a existing NeTCDF file.
%--------------------------------------------------------------------------

Dname.lon ='lon';   Dsize.lon =Im;
Dname.lat ='lat';   Dsize.lat =Jm;
Dname.bath='bath';  Dsize.bath=nc_constant('nc_unlimited');

Vname.lon ='lon';
Vname.lat ='lat';
Vname.bath='hraw';

%--------------------------------------------------------------------------
%  Create topography NetCDF file.
%--------------------------------------------------------------------------

[ncid,status]=mexnc('nccreate',Bname,'nc_write');
if (ncid == -1),
  disp('  ');
  disp(mexnc('strerror',status));
  error(['C_BATH: nccreate - unable to create file: ', Bname]);
end

%--------------------------------------------------------------------------
%  Create global attribute(s).
%--------------------------------------------------------------------------

type='GRID file';
lstr=max(size(type));
[status]=mexnc('ncattput',ncid,nc_constant('nc_global'),'type',         ...
               nc_constant('nc_char'),lstr,type);
if (status == -1),
  disp('  ');
  disp(mexnc('strerror',status));
  error('C_BATH: ncattput - unable to global attribure: type.');
end

history=['GRID file using Matlab script: c_bath, ', date_stamp];
lstr=max(size(history));
[status]=mexnc('ncattput',ncid,nc_constant('nc_global'),'history',      ...
               nc_constant('nc_char'),lstr,history);
if (status == -1),
  disp('  ');
  disp(mexnc('strerror',status));
  error('C_BATH: ncattput - unable to global attribure: history.');
end

%--------------------------------------------------------------------------
%  Define dimensions.
%--------------------------------------------------------------------------

[did.lon]=mexnc('ncdimdef',ncid,Dname.lon,Dsize.lon);
if (did.lon == -1),
  disp('  ');
  disp(mexnc('strerror',status));
  error(['C_BATH: ncdimdef - unable to define dimension: ',Dname.lon]);
end,

[did.lat]=mexnc('ncdimdef',ncid,Dname.lat,Dsize.lat);
if (did.lat == -1),
  disp('  ');
  disp(mexnc('strerror',status));
  error(['C_BATH: ncdimdef - unable to define dimension: ',Dname.lat]);
end,

[did.bath]=mexnc('ncdimdef',ncid,Dname.bath,Dsize.bath);
if (did.bath == -1),
  disp('  ');
  disp(mexnc('strerror',status));
  error(['C_BATH: ncdimdef - unable to define dimension: ',Dname.bath]);
end,

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

[status]=mexnc('ncendef',ncid);
if (status == -1),
  disp('  ');
  disp(mexnc('strerror',status));
  error('C_BATH: ncendef - unable to leave definition mode.');
end

[status]=mexnc('ncclose',ncid);
if (status == -1),
  disp('  ');
  disp(mexnc('strerror',status));
  error(['C_BATH: ncclose - unable to close GRID NetCDF file: ', Bname]);
end

return

