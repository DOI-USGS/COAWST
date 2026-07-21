function [status]=c_initial(S)

%
% C_INITIAL:  Create a ROMS initial conditions NetCDF file
%
% [status]=c_initial(S)
%
% This function creates a ROMS initial conditions NetCDF file using
% specified parameters in structure array, S.
%
% On Input:
%
%    S           Initial condidions creation parameters (structure array):
%
%                  S.ncname           NetCDF file name
%                  S.spherical        Spherical grid switch
%                  S.Vtransform       Vertical transformation equation
%                  S.Lm               Number of interior RHO-points in X
%                  S.Mm               Number of interior RHO-points in Y
%                  S.N                Number of vertical levels
%                  S.NT               Number of active and passive tracers
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
%  Get initial condition creation parameters.
%--------------------------------------------------------------------------

if (isfield(S,'ncname')),
  ncname=S.ncname;
else
  error(['C_INITIAL - Cannot find dimension parameter: ncname, ',       ...
         'in structure array S']);
end

if (isfield(S,'spherical')),
  spherical=S.spherical;
else
  spherical=0;
end

if (isfield(S,'Vtransform')),
  Vtransform=S.Vtransform;
else
  error(['C_INITIAL - Cannot find dimension parameter: Vtransform, ',   ...
         'in structure array S']);
end

if (isfield(S,'Lm')),
  Lp=S.Lm+2;
else
  error(['C_INITIAL - Cannot find dimension parameter: Lm, ',           ...
         'in structure array S']);
end

if (isfield(S,'Mm')),
  Mp=S.Mm+2;
else
  error(['C_INITIAL - Cannot find dimension parameter: Mm, ',           ...
         'in structure array S']);
end

if (isfield(S,'N')),
  N=S.N;
else
  error(['C_INITIAL - Cannot find dimension parameter: N, ',            ...
         'in structure array S']);
end

if (isfield(S,'NT')),
  NT=S.NT;
else
  error(['C_INITIAL - Cannot find dimension parameter: NT, ',           ...
         'in structure S']);
end

%--------------------------------------------------------------------------
%  Create dimensions.
%--------------------------------------------------------------------------

Dname = {'xi_rho',  'xi_u',  'xi_v',  'xi_psi',                       ...
         'eta_rho', 'eta_u', 'eta_v', 'eta_psi',                      ...
         's_rho', 's_w', 'tracer', 'ocean_time'};

for value = Dname
  dim = char(value);
  switch (dim)
    case {'xi_rho', 'xi_v'}
      Dsize.(dim) = Lp;
    case {'xi_psi', 'xi_u'}
      Dsize.(dim) = Lp-1;
    case {'eta_rho', 'eta_u'}
      Dsize.(dim) = Mp;
    case {'eta_psi', 'eta_v'}
      Dsize.(dim) = Mp-1;
    case {'s_rho'}
      Dsize.(dim) = N;
    case {'s_w'}
      Dsize.(dim) = N+1;
    case {'tracer'}
      Dsize.(dim) = NT;
    case {'ocean_time'}
      Dsize.(dim) = netcdf.getConstant('UNLIMITED');
  end
end

%--------------------------------------------------------------------------
%  Set grid Variables.
%--------------------------------------------------------------------------

grd_vars = {'spherical', 'Vtransform', 'Vstretching',                   ...
            'theta_s', 'theta_b', 'Tcline', 'hc',                       ...
            's_rho', 's_w', 'Cs_r', 'Cs_w', 'h'};
	    
if (spherical)
  grd_vars = [grd_vars, 'lon_rho', 'lat_rho', 'lon_psi', 'lat_psi',     ...
                        'lon_u',   'lat_u',   'lon_v',   'lat_v'];
else
  grd_vars = [grd_vars, 'x_rho', 'y_rho', 'x_psi', 'y_psi',             ...
                        'x_u',   'y_u',   'x_v',   'y_v'];
end

grd_vars = [grd_vars, 'mask_rho', 'mask_psi', 'mask_u', 'mask_v'];

%--------------------------------------------------------------------------
%  Set initial conditions variables.
%--------------------------------------------------------------------------

ini_vars = {'ocean_time', 'zeta', 'ubar', 'vbar', 'u', 'v',             ...
            'temp', 'salt'};

%--------------------------------------------------------------------------
%  Create initial conditions NetCDF file.
%--------------------------------------------------------------------------

mode = netcdf.getConstant('CLOBBER');
mode = bitor(mode, netcdf.getConstant('64BIT_OFFSET'));

ncid = netcdf.create(ncname, mode);

%--------------------------------------------------------------------------
%  Define dimensions.
%--------------------------------------------------------------------------

for value = Dname,
  field = char(value);
  did.(field) = netcdf.defDim(ncid, field, Dsize.(field));
end

%--------------------------------------------------------------------------
%  Define NetCDF variables.
%--------------------------------------------------------------------------

vars = [grd_vars, ini_vars];

for value = vars
  field = char(value);
  Vname.(field) = field;
end

% Define spherical switch.

Var.name          = Vname.spherical;
Var.type          = netcdf.getConstant('nc_int');
Var.dimid         = [];
Var.long_name     = 'grid type logical switch';
Var.flag_values   = [0 1];
Var.flag_meanings = ['Cartesian', blanks(1),                            ...
                     'spherical'];
[~,status] = nc_vdef(ncid,Var);
if (status ~= 0), return, end
clear Var

% Define vertical coordinate variables.

Var.name            = Vname.Vtransform;
Var.type            = nc_constant('nc_int');
Var.dimid           = [];
Var.long_name       = 'vertical terrain-following transformation equation';
[~,status]=nc_vdef(ncid,Var);
if (status ~= 0), return, end
clear Var

Var.name            = Vname.Vstretching;
Var.type            = nc_constant('nc_int');
Var.dimid           = [];
Var.long_name       = 'vertical terrain-following stretching function';
[~,status]=nc_vdef(ncid,Var);
if (status ~= 0), return, end
clear Var

Var.name            = Vname.theta_s;
Var.type            = nc_constant('nc_double');
Var.dimid           = [];
Var.long_name       = 'S-coordinate surface control parameter';
[~,status]=nc_vdef(ncid,Var);
if (status ~= 0), return, end
clear Var

Var.name            = Vname.theta_b;
Var.type            = nc_constant('nc_double');
Var.dimid           = [];
Var.long_name       = 'S-coordinate bottom control parameter';
[~,status]=nc_vdef(ncid,Var);
if (status ~= 0), return, end
clear Var

Var.name            = Vname.Tcline;
Var.type            = nc_constant('nc_double');
Var.dimid           = [];
Var.long_name       = 'S-coordinate surface/bottom layer width';
Var.units           = 'meter';
[~,status]=nc_vdef(ncid,Var);
if (status ~= 0), return, end
clear Var

Var.name            = Vname.hc;
Var.type            = nc_constant('nc_double');
Var.dimid           = [];
Var.long_name       = 'S-coordinate parameter, critical depth';
Var.units           = 'meter';
[~,status]=nc_vdef(ncid,Var);
if (status ~= 0), return, end
clear Var

Var.name            = Vname.s_rho;
Var.type            = nc_constant('nc_double');
Var.dimid           = [did.s_rho];
Var.long_name       = 'S-coordinate at RHO-points';
Var.valid_min       = -1;
Var.valid_max       = 0;
Var.positive        = 'up';
if (Vtransform == 1),
  Var.standard_name = 'ocean_s_coordinate_g1';
elseif (Vtransform == 2),
  Var.standard_name = 'ocean_s_coordinate_g2';
end
Var.formula_terms   = 's: s_rho C: Cs_r eta: zeta depth: h depth_c: hc';
[~,status]=nc_vdef(ncid,Var);
if (status ~= 0), return, end
clear Var

Var.name            = Vname.s_w;
Var.type            = nc_constant('nc_double');
Var.dimid           = [did.s_w];
Var.long_name       = 'S-coordinate at W-points';
Var.valid_min       = -1;
Var.valid_max       = 0;
Var.positive        = 'up';
if (Vtransform == 1),
  Var.standard_name = 'ocean_s_coordinate_g1';
elseif (Vtransform == 2),
  Var.standard_name = 'ocean_s_coordinate_g2';
end
Var.formula_terms   = 's: s_w C: Cs_w eta: zeta depth: h depth_c: hc';
[~,status]=nc_vdef(ncid,Var);
if (status ~= 0), return, end
clear Var

Var.name            = Vname.Cs_r;
Var.type            = nc_constant('nc_double');
Var.dimid           = [did.s_rho];
Var.long_name       = 'S-coordinate stretching function at RHO-points';
Var.valid_min       = -1;
Var.valid_max       = 0;
[~,status]=nc_vdef(ncid,Var);
if (status ~= 0), return, end
clear Var

Var.name            = Vname.Cs_w;
Var.type            = nc_constant('nc_double');
Var.dimid           = [did.s_w];
Var.long_name       = 'S-coordinate stretching function at W-points';
Var.valid_min       = -1;
Var.valid_max       = 0;
[~,status]=nc_vdef(ncid,Var);
if (status ~= 0), return, end
clear Var

%  Define bathymetry.

Var.name            = Vname.h;
Var.type            = nc_constant('nc_double');
Var.dimid           = [did.eta_rho did.xi_rho];
Var.long_name       = 'bathymetry at RHO-points';
Var.units           = 'meter';
if (spherical)
  Var.coordinates   = 'lon_rho lat_rho';
else
  Var.coordinates   = 'x_rho y_rho';
end
[~,status]=nc_vdef(ncid,Var);
if (status ~= 0), return, end
clear Var

%  Define horizontal grid variables.

if (spherical)
  Var.name          = Vname.lon_rho;
  Var.type          = nc_constant('nc_double');
  Var.dimid         = [did.eta_rho did.xi_rho];
  Var.long_name     = 'longitude of RHO-points';
  Var.units         = 'degree_east';
  Var.standard_name = 'longitude';
  [~,status]=nc_vdef(ncid,Var);
  if (status ~= 0), return, end
  clear Var

  Var.name          = Vname.lat_rho;
  Var.type          = nc_constant('nc_double');
  Var.dimid         = [did.eta_rho did.xi_rho];
  Var.long_name     = 'latitute of RHO-points';
  Var.units         = 'degree_north';
  Var.standard_name = 'latitude';
  [~,status]=nc_vdef(ncid,Var);
  if (status ~= 0), return, end
  clear Var

  Var.name          = Vname.lon_u;
  Var.type          = nc_constant('nc_double');
  Var.dimid         = [did.eta_u did.xi_u];
  Var.long_name     = 'longitude of U-points';
  Var.units         = 'degree_east';
  Var.standard_name = 'longitude';
  [~,status]=nc_vdef(ncid,Var);
  if (status ~= 0), return, end
  clear Var

  Var.name          = Vname.lat_u;
  Var.type          = nc_constant('nc_double');
  Var.dimid         = [did.eta_u did.xi_u];
  Var.long_name     = 'latitute of U-points';
  Var.units         = 'degree_north';
  Var.standard_name = 'latitude';
  [~,status]=nc_vdef(ncid,Var);
  if (status ~= 0), return, end
  clear Var

  Var.name          = Vname.lon_v;
  Var.type          = nc_constant('nc_double');
  Var.dimid         = [did.eta_v did.xi_v];
  Var.long_name     = 'longitude of V-points';
  Var.units         = 'degree_east';
  Var.standard_name = 'longitude';
  [~,status]=nc_vdef(ncid,Var);
  if (status ~= 0), return, end
  clear Var

  Var.name          = Vname.lat_v;
  Var.type          = nc_constant('nc_double');
  Var.dimid         = [did.eta_v did.xi_v];
  Var.long_name     = 'latitute of V-points';
  Var.units         = 'degree_north';
  Var.standard_name = 'latitude';
  [~,status]=nc_vdef(ncid,Var);
  if (status ~= 0), return, end
  clear Var

else

  Var.name          = Vname.x_rho;
  Var.type          = nc_constant('nc_double');
  Var.dimid         = [did.eta_rho did.xi_rho];
  Var.long_name     = 'X-location of RHO-points';
  Var.units         = 'meter';
  [~,status]=nc_vdef(ncid,Var);
  if (status ~= 0), return, end
  clear Var

  Var.name          = Vname.y_rho;
  Var.type          = nc_constant('nc_double');
  Var.dimid         = [did.eta_rho did.xi_rho];
  Var.long_name     = 'Y-location of RHO-points';
  Var.units         = 'meter';
  [~,status]=nc_vdef(ncid,Var);
  if (status ~= 0), return, end
  clear Var

  Var.name          = Vname.x_u;
  Var.type          = nc_constant('nc_double');
  Var.dimid         = [did.eta_u did.xi_u];
  Var.long_name     = 'X-location of U-points';
  Var.units         = 'meter';
  [~,status]=nc_vdef(ncid,Var);
  if (status ~= 0), return, end
  clear Var

  Var.name          = Vname.y_u;
  Var.type          = nc_constant('nc_double');
  Var.dimid         = [did.eta_u did.xi_u];
  Var.long_name     = 'Y-location of U-points';
  Var.units         = 'meter';
  [~,status]=nc_vdef(ncid,Var);
  if (status ~= 0), return, end
  clear Var

  Var.name          = Vname.x_v;
  Var.type          = nc_constant('nc_double');
  Var.dimid         = [did.eta_v did.xi_v];
  Var.long_name     = 'X-location of V-points';
  Var.units         = 'meter';
  [~,status]=nc_vdef(ncid,Var);
  if (status ~= 0), return, end
  clear Var

  Var.name          = Vname.y_v;
  Var.type          = nc_constant('nc_double');
  Var.dimid         = [did.eta_v did.xi_v];
  Var.long_name     = 'Y-location of V-points';
  Var.units         = 'meter';
  [~,status]=nc_vdef(ncid,Var);
  if (status ~= 0), return, end
  clear Var
  
end

%--------------------------------------------------------------------------
%  Define initial conditions variables.
%--------------------------------------------------------------------------

Var.name            = Vname.ocean_time;
Var.type            = nc_constant('nc_double');
Var.dimid           = [did.ocean_time];
Var.long_name       = 'time since initialization';
Var.units           = 'seconds';
[~,status]=nc_vdef(ncid,Var);
if (status ~= 0), return, end
clear Var

Var.name            = Vname.zeta;
Var.type            = nc_constant('nc_double');
Var.dimid           = [did.ocean_time did.eta_rho did.xi_rho];
Var.long_name       = 'free-surface';
Var.units           = 'meter';
Var.time            = Vname.ocean_time;
if (spherical)
  Var.coordinates   = strcat([Vname.lon_rho,' ',Vname.lat_rho,' ',      ...
                              Vname.ocean_time]); 
else
  Var.coordinates   = strcat([Vname.x_rho,' ',Vname.y_rho,' ',          ...
                              Vname.ocean_time]); 
end
[~,status]=nc_vdef(ncid,Var);
if (status ~= 0), return, end
clear Var

Var.name            = Vname.ubar;
Var.type            = nc_constant('nc_double');
Var.dimid           = [did.ocean_time did.eta_u did.xi_u];
Var.long_name       = 'vertically integrated u-momentum component';
Var.units           = 'meter second-1';
Var.time            = Vname.ocean_time;
if (spherical)
  Var.coordinates   = strcat([Vname.lon_u,' ',Vname.lat_u,' ',          ...
                              Vname.ocean_time]); 
else
  Var.coordinates   = strcat([Vname.x_u,' ',Vname.y_u,' ',              ...
                              Vname.ocean_time]); 
end
[~,status]=nc_vdef(ncid,Var);
if (status ~= 0), return, end
clear Var

Var.name            = Vname.vbar;
Var.type            = nc_constant('nc_double');
Var.dimid           = [did.ocean_time did.eta_v did.xi_v];
Var.long_name       = 'vertically integrated v-momentum component';
Var.units           = 'meter second-1';
Var.time            = Vname.ocean_time;
if (spherical)
  Var.coordinates   = strcat([Vname.lon_v,' ',Vname.lat_v,' ',          ...
                              Vname.ocean_time]); 
else
  Var.coordinates   = strcat([Vname.x_v,' ',Vname.y_v,' ',              ...
                              Vname.ocean_time]); 
end
[~,status]=nc_vdef(ncid,Var);
if (status ~= 0), return, end
clear Var

Var.name            = Vname.u;
Var.type            = nc_constant('nc_double');
Var.dimid           = [did.ocean_time did.s_rho did.eta_u did.xi_u];
Var.long_name       = 'u-momentum component';
Var.units           = 'meter second-1';
Var.time            = Vname.ocean_time;
if (spherical)
  Var.coordinates   = strcat([Vname.lon_u,' ',Vname.lat_u,' ',          ...
                              Vname.s_rho,' ',Vname.ocean_time]); 
else
  Var.coordinates   = strcat([Vname.x_rho,' ',Vname.y_rho,' ',          ...
                              Vname.s_rho,' ',Vname.ocean_time]); 
end
[~,status]=nc_vdef(ncid,Var);
if (status ~= 0), return, end
clear Var

Var.name            = Vname.v;
Var.type            = nc_constant('nc_double');
Var.dimid           = [did.ocean_time did.s_rho did.eta_v did.xi_v];
Var.long_name       = 'v-momentum component';
Var.units           = 'meter second-1';
Var.time            = Vname.ocean_time;
if (spherical)
  Var.coordinates   = strcat([Vname.lon_v,' ',Vname.lat_v,' ',          ...
                              Vname.s_rho,' ',Vname.ocean_time]); 
else
  Var.coordinates   = strcat([Vname.x_v,' ',Vname.y_v,' ',              ...
                              Vname.s_rho,' ',Vname.ocean_time]); 
end
[~,status]=nc_vdef(ncid,Var);
if (status ~= 0), return, end
clear Var

Var.name            = Vname.temp;
Var.type            = nc_constant('nc_double');
Var.dimid           = [did.ocean_time did.s_rho did.eta_rho did.xi_rho];
Var.long_name       = 'potential temperature';
Var.units           = 'Celsius';
Var.time            = Vname.ocean_time;
if (spherical)
  Var.coordinates   = strcat([Vname.lon_rho,' ',Vname.lat_rho,' ',      ...
                              Vname.s_rho,' ',Vname.ocean_time]); 
else
  Var.coordinates   = strcat([Vname.x_rho,' ',Vname.y_rho,' ',          ...
                              Vname.s_rho,' ',Vname.ocean_time]); 
end
[~,status]=nc_vdef(ncid,Var);
if (status ~= 0), return, end
clear Var

Var.name            = Vname.salt;
Var.type            = nc_constant('nc_double');
Var.dimid           = [did.ocean_time did.s_rho did.eta_rho did.xi_rho];
Var.long_name       = 'salinity';
Var.time            = Vname.ocean_time;
if (spherical)
  Var.coordinates   = strcat([Vname.lon_rho,' ',Vname.lat_rho,' ',      ...
                              Vname.s_rho,' ',Vname.ocean_time]); 
else
  Var.coordinates   = strcat([Vname.x_rho,' ',Vname.y_rho,' ',          ...
                              Vname.s_rho,' ',Vname.ocean_time]); 
end
[~,status]=nc_vdef(ncid,Var);
if (status ~= 0), return, end
clear Var

%--------------------------------------------------------------------------
%  Create global attributes.
%--------------------------------------------------------------------------

varid  = netcdf.getConstant('nc_global');

type='ROMS INITIAL file';
netcdf.putAtt(ncid, varid, 'type', type);

history=['Initial file using Matlab script: c_initial, ',date_stamp];
netcdf.putAtt(ncid, varid, 'history', history);

%--------------------------------------------------------------------------
%  Leave definition mode and close NetCDF file.
%--------------------------------------------------------------------------

netcdf.endDef(ncid);
netcdf.close(ncid);

return

