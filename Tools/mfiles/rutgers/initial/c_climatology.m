function [status]=c_climatology(S)

%
% C_CLIMATOLOGY:  Create a ROMS climatology NetCDF file
%
% [status]=c_climatology(S)
%
% This function creates a ROMS climatology NetCDF file using specified
% parameters in structure array, S.
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
%                  S.def_zeta         Switch to define sea surfac height
%                  S.def_v2d          Switch to define 2D momentum
%                  S.def_v3d          Switch to define 3D momentum
%                  S.def_temp         Switch to define potential temperature
%                  S.def_salt         Switch to define salinity
%
%                  S.zeta_time        Number of sea surface height records
%                  S.v2d_time         Number of 2D momentum records
%                  S.v3d_time         Number of 3D momentum records
%                  S.temp_time        Number of potential temperature records
%                  S.salt_time        Number of salinity records
%
% On Output:
%
%    status      Error flag.
%

% svn $Id: c_climatology.m 996 2020-01-10 04:28:56Z arango $
%=========================================================================%
%  Copyright (c) 2002-2020 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%

%--------------------------------------------------------------------------
%  Get climatology creation parameters.
%--------------------------------------------------------------------------

if (isfield(S,'ncname')),
  ncname=S.ncname;
else
  error(['C_CLIMATOLOGY - Cannot find dimension parameter: ncname, ',   ...
         'in structure array S']);
end,

if (isfield(S,'spherical')),
  spherical=S.spherical;
else
  spherical=0;
end

if (isfield(S,'Vtransform')),
  Vtransform=S.Vtransform;
else
  error(['C_CLIMATOLOGY - Cannot find dimension parameter: ',           ...
         'Vtransform, in structure array S']);
end

if (isfield(S,'Lm')),
  Lp=S.Lm+2;
else
  error(['C_CLIMATOLOGY - Cannot find dimension parameter: Lm, ',       ...
          'in structure array S']);
end

if (isfield(S,'Mm')),
  Mp=S.Mm+2;
else
  error(['C_CLIMATOLOGY - Cannot find dimension parameter: Mm, ',       ...
         'in structure array S']);
end,

if (isfield(S,'N')),
  N=S.N;
else
  error(['C_CLIMATOLOGY - Cannot find dimension parameter: N, ',        ...
         'in structure array S']);
end

if (isfield(S,'NT')),
  NT=S.NT;
else
  error(['C_CLIMATOLOGY - Cannot find dimension parameter: NT, ',       ...
         'in structure array S']);
end

if (isfield(S,'def_zeta')),
  if (S.def_zeta && ~isfield(S,'zeta_time')),
    error(['C_CLIMATOLOGY - Cannot find dimension parameter: ',         ...
           'zeta_time, in structure array S']);
  end
else  
  S.def_zeta=0;
  S.zeta_time=0;
end

if (isfield(S,'def_v2d')),
  if (S.def_v2d && ~isfield(S,'v2d_time')),
    error(['C_CLIMATOLOGY - Cannot find dimension parameter: ',         ...
           'v2d_time, in structure array S']);
  end
else   
  S.def_v2d=0;
  S.v2d_time=0;
end,

if (isfield(S,'def_v3d')),
  if (S.def_v3d && ~isfield(S,'v3d_time')),
    error(['C_CLIMATOLOGY - Cannot find dimension parameter: ',         ...
           'v3d_time, in structure array S']);
  end
else   
  S.def_v3d=0;
  S.v3d_time=0;
end,

if (isfield(S,'def_temp')),
  if (S.def_temp && ~isfield(S,'temp_time')),
    error(['C_CLIMATOLOGY - Cannot find dimension parameter: ',         ...
           'temp_time, in structure array S']);
  end
else  
  S.def_temp=0;
  S.temp_time=0;
end,

if (isfield(S,'def_salt')),
  if (S.def_salt && ~isfield(S,'salt_time')),
    error(['C_CLIMATOLOGY - Cannot find dimension parameter: ',         ...
           'salt_time, in structure array S']);
  end
else  
  S.def_salt=0;
  S.salt_time=0;
end
  
%--------------------------------------------------------------------------
%  Set dimensions.
%--------------------------------------------------------------------------

Dname.xr = 'xi_rho';       Dsize.xr = Lp;
Dname.xu = 'xi_u';         Dsize.xu = Lp-1;
Dname.xv = 'xi_v';         Dsize.xv = Lp;
Dname.xp = 'xi_psi';       Dsize.xp = Lp-1;
Dname.yr = 'eta_rho';      Dsize.yr = Mp;
Dname.yu = 'eta_u';        Dsize.yu = Mp;
Dname.yv = 'eta_v';        Dsize.yv = Mp-1;
Dname.yp = 'eta_psi';      Dsize.yp = Mp-1;
Dname.Nr = 's_rho';        Dsize.Nr = N;
Dname.Nw = 's_w';          Dsize.Nw = N+1;
Dname.NT = 'tracer';       Dsize.NT = NT;

if (S.def_zeta),
 Dname.zeta_time = 'zeta_time';    Dsize.zeta_time = S.zeta_time;
end

if (S.def_v2d),
 Dname.v2d_time  = 'v2d_time';     Dsize.v2d_time  = S.v2d_time;
end

if (S.def_v3d),
 Dname.v2d_time  = 'v3d_time';     Dsize.v3d_time  = S.v3d_time;
end

if (S.def_temp),
 Dname.temp_time = 'temp_time';    Dsize.temp_time = S.temp_time;
end

if (S.def_salt),
 Dname.salt_time = 'salt_time';    Dsize.salt_time = S.temp_time;
end

%--------------------------------------------------------------------------
%  Set Variables.
%--------------------------------------------------------------------------

%  Vertical grid variables.

Vname.Vtransform  = 'Vtransform';
Vname.Vstretching = 'Vstretching';
Vname.theta_s     = 'theta_s';
Vname.theta_b     = 'theta_b';
Vname.Tcline      = 'Tcline';
Vname.hc          = 'hc';
Vname.s_rho       = 's_rho';
Vname.s_w         = 's_w';
Vname.Cs_r        = 'Cs_r';
Vname.Cs_w        = 'Cs_w';

%  Horizontal grid variables.

Vname.spherical   = 'spherical';
Vname.h           = 'h';

if (spherical),
  Vname.rlon      = 'lon_rho';
  Vname.rlat      = 'lat_rho';
  Vname.ulon      = 'lon_u';
  Vname.ulat      = 'lat_u';
  Vname.vlon      = 'lon_v';
  Vname.vlat      = 'lat_v';
else
  Vname.rx        = 'x_rho';
  Vname.ry        = 'y_rho';
  Vname.ux        = 'x_u';
  Vname.uy        = 'y_u';
  Vname.vx        = 'x_v';
  Vname.vy        = 'y_v';
end

%  climatology variables.

Vname.zeta_time   = 'zeta_time';
Vname.v2d_time    = 'v2d_time';
Vname.v3d_time    = 'v3d_time';
Vname.temp_time   = 'temp_time';
Vname.salt_time   = 'salt_time';

Vname.zeta        = 'zeta';
Vname.ubar        = 'ubar';
Vname.vbar        = 'vbar';
Vname.u           = 'u';
Vname.v           = 'v';
Vname.temp        = 'temp';
Vname.salt        = 'salt';

%--------------------------------------------------------------------------
%  Create initial conditions NetCDF file.
%--------------------------------------------------------------------------

[ncid,status]=mexnc('create',ncname,'clobber');
if (status ~= 0),
  disp('  ');
  disp(mexnc('strerror',status));
  error(['C_CLIMATOLOGY: CREATE - unable to create file: ', ncname]);
end

%--------------------------------------------------------------------------
%  Define dimensions.
%--------------------------------------------------------------------------

[did.xr,status]=mexnc('def_dim',ncid,Dname.xr,Dsize.xr); 
if (status ~= 0),
  disp('  ');
  disp(mexnc('strerror',status));
  error(['C_CLIMATOLOGY: DEF_DIM - unable to define dimension: ',Dname.xr]);
end

[did.xu,status]=mexnc('def_dim',ncid,Dname.xu,Dsize.xu);
if (status ~= 0),
  disp('  ');
  disp(mexnc('strerror',status));
  error(['C_CLIMATOLOGY: DEF_DIM - unable to define dimension: ',Dname.xu]);
end

[did.xv,status]=mexnc('def_dim',ncid,Dname.xv,Dsize.xv);
if (status ~= 0),
  disp('  ');
  disp(mexnc('strerror',status));
  error(['C_CLIMATOLOGY: DEF_DIM - unable to define dimension: ',Dname.xv]);
end

[did.yr,status]=mexnc('def_dim',ncid,Dname.yr,Dsize.yr);
if (status ~= 0),
  disp('  ');
  disp(mexnc('strerror',status));
  error(['C_CLIMATOLOGY: DEF_DIM - unable to define dimension: ',Dname.yr]);
end

[did.yu,status]=mexnc('def_dim',ncid,Dname.yu,Dsize.yu);
if (status ~= 0),
  disp('  ');
  disp(mexnc('strerror',status));
  error(['C_CLIMATOLOGY: DEF_DIM - unable to define dimension: ',Dname.yu]);
end

[did.yv,status]=mexnc('def_dim',ncid,Dname.yv,Dsize.yv);
if (status ~= 0),
  disp('  ');
  disp(mexnc('strerror',status));
  error(['C_CLIMATOLOGY: DEF_DIM - unable to define dimension: ',Dname.yv]);
end

[did.Nr,status]=mexnc('def_dim',ncid,Dname.Nr,Dsize.Nr);
if (status ~= 0),
  disp('  ');
  disp(mexnc('strerror',status));
  error(['C_CLIMATOLOGY: DEF_DIM - unable to define dimension: ',Dname.Nr]);
end

[did.Nw,status]=mexnc('def_dim',ncid,Dname.Nw,Dsize.Nw);
if (status ~= 0),
  disp('  ');
  disp(mexnc('strerror',status));
  error(['C_CLIMATOLOGY: DEF_DIM - unable to define dimension: ',Dname.Nw]);
end

[did.NT,status]=mexnc('def_dim',ncid,Dname.NT,Dsize.NT);
if (status ~= 0),
  disp('  ');
  disp(mexnc('strerror',status));
  error(['C_CLIMATOLOGY: DEF_DIM - unable to define dimension: ',Dname.NT]);
end

if (S.def_zeta),
  [did.zeta_time,status]=mexnc('def_dim',ncid,Dname.zeta_time,Dsize.zeta_time);
  if (status ~= 0),
    disp('  ');
    disp(mexnc('strerror',status));
    error(['C_CLIMATOLOGY: DEF_DIM - unable to define dimension: ',     ...
           Dname.zeta_time]);
  end
end

if (S.def_v2d),
  [did.v2d_time,status]=mexnc('def_dim',ncid,Dname.v2d_time,Dsize.v2d_time);
  if (status ~= 0),
    disp('  ');
    disp(mexnc('strerror',status));
    error([ 'C_CLIMATOLOGY: DEF_DIM - unable to define dimension: ',    ...
	    Dname.v2d_time]);
  end
end

if (S.def_v3d),
  [did.v3d_time,status]=mexnc('def_dim',ncid,Dname.v3d_time,Dsize.v3d_time);
  if (status ~= 0),
    disp('  ');
    disp(mexnc('strerror',status));
    error(['C_CLIMATOLOGY: DEF_DIM - unable to define dimension: ',     ...
           Dname.v3d_time]);
  end
end

if (S.def_temp),
  [did.temp_time,status]=mexnc('def_dim',ncid,Dname.temp_time,Dsize.temp_time);
  if (status ~= 0),
    disp('  ');
    disp(mexnc('strerror',status));
    error(['C_CLIMATOLOGY: DEF_DIM - unable to define dimension: ',     ...
           Dname.temp_time]);
  end
end

if (S.def_salt),
  [did.salt_time,status]=mexnc('def_dim',ncid,Dname.salt_time,Dsize.salt_time);
  if (status ~= 0),
    disp('  ');
    disp(mexnc('strerror',status));
    error(['C_CLIMATOLOGY: DEF_DIM - unable to define dimension: ',     ...
           Dname.salt_time]);
  end
end

%--------------------------------------------------------------------------
%  Create global attributes.
%--------------------------------------------------------------------------

type='CLIMATOLOGY file';
lstr=max(size(type));
[status]=mexnc('put_att_text',ncid,nc_constant('nc_global'),            ...
               'type',nc_constant('nc_char'),lstr,type);
if (status ~= 0),
  disp('  ');
  disp(mexnc('strerror',status));
  error('C_CLIMATOLOGY: PUT_ATT_TEXT - unable to global attribure: type.');
end

history=['Climatology file using Matlab script: c_climatology, ',date_stamp];
lstr=max(size(history));
[status]=mexnc('put_att_text',ncid,nc_constant('nc_global'),            ...
               'history',nc_constant('nc_char'),lstr,history);
if (status ~= 0),
  disp('  ');
  disp(mexnc('strerror',status));
  error('C_CLIMATOLOGY: PUT_ATT_TEXT - unable to global attribure: history.');
end

%--------------------------------------------------------------------------
%  Define configuration variables.
%--------------------------------------------------------------------------

% Define spherical switch.

Var.name            = Vname.spherical;
Var.type            = nc_constant('nc_int');
Var.dimid           = [];
Var.long_name       = 'grid type logical switch';
Var.flag_values     = [0 1];
Var.flag_meanings   = ['Cartesian', blanks(1), ...
                       'spherical'];
[~,status]=nc_vdef(ncid,Var);
if (status ~= 0), return, end,
clear Var

% Define vertical coordinate variables.

Var.name            = Vname.Vtransform;
Var.type            = nc_constant('nc_int');
Var.dimid           = [];
Var.long_name       = 'vertical terrain-following transformation equation';
[~,status]=nc_vdef(ncid,Var);
if (status ~= 0), return, end,
clear Var

Var.name            = Vname.Vstretching;
Var.type            = nc_constant('nc_int');
Var.dimid           = [];
Var.long_name       = 'vertical terrain-following stretching function';
[~,status]=nc_vdef(ncid,Var);
if (status ~= 0), return, end,
clear Var

Var.name            = Vname.theta_s;
Var.type            = nc_constant('nc_double');
Var.dimid           = [];
Var.long_name       = 'S-coordinate surface control parameter';
[~,status]=nc_vdef(ncid,Var);
if (status ~= 0), return, end,
clear Var

Var.name            = Vname.theta_b;
Var.type            = nc_constant('nc_double');
Var.dimid           = [];
Var.long_name       = 'S-coordinate bottom control parameter';
[~,status]=nc_vdef(ncid,Var);
if (status ~= 0), return, end,
clear Var

Var.name            = Vname.Tcline;
Var.type            = nc_constant('nc_double');
Var.dimid           = [];
Var.long_name       = 'S-coordinate surface/bottom layer width';
Var.units           = 'meter';
[~,status]=nc_vdef(ncid,Var);
if (status ~= 0), return, end,
clear Var

Var.name            = Vname.hc;
Var.type            = nc_constant('nc_double');
Var.dimid           = [];
Var.long_name       = 'S-coordinate parameter, critical depth';
Var.units           = 'meter';
[~,status]=nc_vdef(ncid,Var);
if (status ~= 0), return, end,
clear Var

Var.name            = Vname.s_rho;
Var.type            = nc_constant('nc_double');
Var.dimid           = [did.Nr];
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
if (status ~= 0), return, end,
clear Var

Var.name            = Vname.s_w;
Var.type            = nc_constant('nc_double');
Var.dimid           = [did.Nw];
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
if (status ~= 0), return, end,
clear Var

Var.name            = Vname.Cs_r;
Var.type            = nc_constant('nc_double');
Var.dimid           = [did.Nr];
Var.long_name       = 'S-coordinate stretching function at RHO-points';
Var.valid_min       = -1;
Var.valid_max       = 0;
[~,status]=nc_vdef(ncid,Var);
if (status ~= 0), return, end,
clear Var

Var.name            = Vname.Cs_w;
Var.type            = nc_constant('nc_double');
Var.dimid           = [did.Nw];
Var.long_name       = 'S-coordinate stretching function at W-points';
Var.valid_min       = -1;
Var.valid_max       = 0;
[~,status]=nc_vdef(ncid,Var);
if (status ~= 0), return, end,
clear Var

%  Define bathymetry.

Var.name            = Vname.h;
Var.type            = nc_constant('nc_double');
Var.dimid           = [did.yr did.xr];
Var.long_name       = 'bathymetry at RHO-points';
Var.units           = 'meter';
if (spherical),
  Var.coordinates   = 'lon_rho lat_rho';
else
  Var.coordinates   = 'x_rho y_rho';
end,
[~,status]=nc_vdef(ncid,Var);
if (status ~= 0), return, end,
clear Var

%  Define horizontal grid variables.

if (spherical),
  Var.name          = Vname.rlon;
  Var.type          = nc_constant('nc_double');
  Var.dimid         = [did.yr did.xr];
  Var.long_name     = 'longitude of RHO-points';
  Var.units         = 'degree_east';
  Var.standard_name = 'longitude';
  [~,status]=nc_vdef(ncid,Var);
  if (status ~= 0), return, end,
  clear Var

  Var.name          = Vname.rlat;
  Var.type          = nc_constant('nc_double');
  Var.dimid         = [did.yr did.xr];
  Var.long_name     = 'latitute of RHO-points';
  Var.units         = 'degree_north';
  Var.standard_name = 'latitude';
  [~,status]=nc_vdef(ncid,Var);
  if (status ~= 0), return, end,
  clear Var

  Var.name          = Vname.ulon;
  Var.type          = nc_constant('nc_double');
  Var.dimid         = [did.yu did.xu];
  Var.long_name     = 'longitude of U-points';
  Var.units         = 'degree_east';
  Var.standard_name = 'longitude';
  [~,status]=nc_vdef(ncid,Var);
  if (status ~= 0), return, end,
  clear Var

  Var.name          = Vname.ulat;
  Var.type          = nc_constant('nc_double');
  Var.dimid         = [did.yu did.xu];
  Var.long_name     = 'latitute of U-points';
  Var.units         = 'degree_north';
  Var.standard_name = 'latitude';
  [~,status]=nc_vdef(ncid,Var);
  if (status ~= 0), return, end,
  clear Var

  Var.name          = Vname.vlon;
  Var.type          = nc_constant('nc_double');
  Var.dimid         = [did.yv did.xv];
  Var.long_name     = 'longitude of V-points';
  Var.units         = 'degree_east';
  Var.standard_name = 'longitude';
  [~,status]=nc_vdef(ncid,Var);
  if (status ~= 0), return, end,
  clear Var

  Var.name          = Vname.vlat;
  Var.type          = nc_constant('nc_double');
  Var.dimid         = [did.yv did.xv];
  Var.long_name     = 'latitute of V-points';
  Var.units         = 'degree_north';
  Var.standard_name = 'latitude';
  [~,status]=nc_vdef(ncid,Var);
  if (status ~= 0), return, end,
  clear Var

else

  Var.name          = Vname.rx;
  Var.type          = nc_constant('nc_double');
  Var.dimid         = [did.yr did.xr];
  Var.long_name     = 'X-location of RHO-points';
  Var.units         = 'meter';
  [~,status]=nc_vdef(ncid,Var);
  if (status ~= 0), return, end,
  clear Var

  Var.name          = Vname.ry;
  Var.type          = nc_constant('nc_double');
  Var.dimid         = [did.yr did.xr];
  Var.long_name     = 'Y-location of RHO-points';
  Var.units         = 'meter';
  [~,status]=nc_vdef(ncid,Var);
  if (status ~= 0), return, end,
  clear Var

  Var.name          = Vname.ux;
  Var.type          = nc_constant('nc_double');
  Var.dimid         = [did.yu did.xu];
  Var.long_name     = 'X-location of U-points';
  Var.units         = 'meter';
  [~,status]=nc_vdef(ncid,Var);
  if (status ~= 0), return, end,
  clear Var

  Var.name          = Vname.uy;
  Var.type          = nc_constant('nc_double');
  Var.dimid         = [did.yu did.xu];
  Var.long_name     = 'Y-location of U-points';
  Var.units         = 'meter';
  [~,status]=nc_vdef(ncid,Var);
  if (status ~= 0), return, end,
  clear Var

  Var.name          = Vname.vx;
  Var.type          = nc_constant('nc_double');
  Var.dimid         = [did.yv did.xv];
  Var.long_name     = 'X-location of V-points';
  Var.units         = 'meter';
  [~,status]=nc_vdef(ncid,Var);
  if (status ~= 0), return, end,
  clear Var

  Var.name          = Vname.vy;
  Var.type          = nc_constant('nc_double');
  Var.dimid         = [did.yv did.xv];
  Var.long_name     = 'Y-location of V-points';
  Var.units         = 'meter';
  [~,status]=nc_vdef(ncid,Var);
  if (status ~= 0), return, end,
  clear Var
  
end

%--------------------------------------------------------------------------
%  Define climatology variables.
%--------------------------------------------------------------------------

%  Time variables.

if (S.def_zeta),
  Var.name          = Vname.zeta_time;
  Var.type          = nc_constant('nc_double');
  Var.dimid         = [did.zeta_time];
  Var.long_name     = 'time for sea surface height climatology';
  Var.units         = 'day';
  [~,status]=nc_vdef(ncid,Var);
  if (status ~= 0), return, end,
  clear Var
end

if (S.def_v2d),
  Var.name          = Vname.v2d_time;
  Var.type          = nc_constant('nc_double');
  Var.dimid         = [did.v2d_time];
  Var.long_name     = 'time for 2D momentum climatology';
  Var.units         = 'day';
  [~,status]=nc_vdef(ncid,Var);
  if (status ~= 0), return, end,
  clear Var
end

if (S.def_v3d),
  Var.name          = Vname.v3d_time;
  Var.type          = nc_constant('nc_double');
  Var.dimid         = [did.v3d_time];
  Var.long_name     = 'time for 3D momentum climatology';
  Var.units         = 'day';
  [~,status]=nc_vdef(ncid,Var);
  if (status ~= 0), return, end,
  clear Var
end

if (S.def_temp),
  Var.name          = Vname.temp_time;
  Var.type          = nc_constant('nc_double');
  Var.dimid         = [did.temp_time];
  Var.long_name     = 'time for potential temperature climatology';
  Var.units         = 'day';
  [~,status]=nc_vdef(ncid,Var);
  if (status ~= 0), return, end,
  clear Var
end,

if (S.def_salt),
  Var.name          = Vname.salt_time;
  Var.type          = nc_constant('nc_double');
  Var.dimid         = [did.salt_time];
  Var.long_name     = 'time for salinity climatology';
  Var.units         = 'day';
  [~,status]=nc_vdef(ncid,Var);
  if (status ~= 0), return, end,
  clear Var
end

%  Climatology fields.

if (S.def_zeta),
  Var.name          = Vname.zeta;
  Var.type          = nc_constant('nc_double');
  Var.dimid         = [did.zeta_time did.yr did.xr];
  Var.long_name     = 'sea surface height climatology';
  Var.units         = 'meter';
  Var.time          = Vname.zeta_time;
  if (spherical),
    Var.coordinates = strcat([Vname.rlon,' ',Vname.rlat,                ...
                              ' ',Vname.zeta_time]); 
  else
    Var.coordinates = strcat([Vname.rx,' ',Vname.ry,                    ...
                              ' ',Vname.zeta_time]); 
  end
  [~,status]=nc_vdef(ncid,Var);
  if (status ~= 0), return, end,
  clear Var
end

if (S.def_v2d),
  Var.name          = Vname.ubar;
  Var.type          = nc_constant('nc_double');
  Var.dimid         = [did.v2d_time did.yu did.xu];
  Var.long_name     = 'vertically integrated u-momentum component climatology';
  Var.units         = 'meter second-1';
  Var.time          = Vname.v2d_time;
  if (spherical),
    Var.coordinates = strcat([Vname.ulon,' ',Vname.ulat,                ...
                              ' ',Vname.v2d_time]);
  else
    Var.coordinates = strcat([Vname.ux,' ',Vname.uy,                    ...
                              ' ',Vname.v2d_time]);
  end
  [~,status]=nc_vdef(ncid,Var);
  if (status ~= 0), return, end,
  clear Var

  Var.name          = Vname.vbar;
  Var.type          = nc_constant('nc_double');
  Var.dimid         = [did.v2d_time did.yv did.xv];
  Var.long_name     = 'vertically integrated v-momentum component climatology';
  Var.units         = 'meter second-1';
  Var.time          = Vname.v2d_time;
  if (spherical),
    Var.coordinates = strcat([Vname.vlon,' ',Vname.vlat,                ...
                              ' ',Vname.v2d_time]); 
  else
    Var.coordinates = strcat([Vname.vx,' ',Vname.vy,                    ...
                              ' ',Vname.v2d_time]); 
  end
  [~,status]=nc_vdef(ncid,Var);
  if (status ~= 0), return, end,
  clear Var
end

if (S.def_v3d),
  Var.name          = Vname.u;
  Var.type          = nc_constant('nc_double');
  Var.dimid         = [did.v3d_time did.Nr did.yu did.xu];
  Var.long_name     = 'u-momentum component climatology';
  Var.units         = 'meter second-1';
  Var.time          = Vname.v3d_time;
  if (spherical),
    Var.coordinates = strcat([Vname.ulon,' ',Vname.ulat,' ',            ...
                              Vname.s_rho,' ',Vname.v3d_time]); 
  else
    Var.coordinates = strcat([Vname.ux,' ',Vname.uy,' ',                ...
                              Vname.s_rho,' ',Vname.v3d_time]); 
  end
  [~,status]=nc_vdef(ncid,Var);
  if (status ~= 0), return, end,
  clear Var

  Var.name          = Vname.v;
  Var.type          = nc_constant('nc_double');
  Var.dimid         = [did.v3d_time did.Nr did.yv did.xv];
  Var.long_name     = 'v-momentum component climatology';
  Var.units         = 'meter second-1';
  Var.time          = Vname.v3d_time;
  if (spherical),
    Var.coordinates = strcat([Vname.vlon,' ',Vname.vlat,' ',            ...
                              Vname.s_rho,' ',Vname.v3d_time]); 
  else
    Var.coordinates = strcat([Vname.vx,' ',Vname.vy,' ',                ...
                              Vname.s_rho,' ',Vname.v3d_time]); 
  end
  [~,status]=nc_vdef(ncid,Var);
  if (status ~= 0), return, end,
  clear Var
end

if (S.def_temp),
  Var.name          = Vname.temp;
  Var.type          = nc_constant('nc_double');
  Var.dimid         = [did.temp_time did.Nr did.yr did.xr];
  Var.long_name     = 'potential temperature climatology';
  Var.units         = 'Celsius';
  Var.time          = Vname.temp_time;
  if (spherical),
    Var.coordinates = strcat([Vname.rlon,' ',Vname.rlat,' ',            ...
                              Vname.s_rho,' ',Vname.temp_time]); 
  else
    Var.coordinates = strcat([Vname.rx,' ',Vname.ry,' ',                ...
                              Vname.s_rho,' ',Vname.temp_time]); 
  end
  [~,status]=nc_vdef(ncid,Var);
  if (status ~= 0), return, end,
  clear Var
end

if (S.def_salt),
  Var.name          = Vname.salt;
  Var.type          = nc_constant('nc_double');
  Var.dimid         = [did.salt_time did.Nr did.yr did.xr];
  Var.long_name     = 'salinity climatology';
  Var.time          = Vname.salt_time;
  if (spherical),
    Var.coordinates = strcat([Vname.rlon,' ',Vname.rlat,' ',            ...
                              Vname.s_rho,' ',Vname.salt_time]); 
  else
    Var.coordinates = strcat([Vname.rx,' ',Vname.ry,' ',                ...
                              Vname.s_rho,' ',Vname.salt_time]); 
  end
  [~,status]=nc_vdef(ncid,Var);
  if (status ~= 0), return, end,
  clear Var
end

%--------------------------------------------------------------------------
%  Leave definition mode and close NetCDF file.
%--------------------------------------------------------------------------

[status]=mexnc('enddef',ncid);
if (status == -1),
  disp('  ');
  disp(mexnc('strerror',status));
  error('C_CLIMATOLOGY: ENDDEF - unable to leave definition mode.');
end

[status]=mexnc('close',ncid);
if (status == -1),
  disp('  ');
  disp(mexnc('strerror',status));
  error([ 'C_CLIMATOLOGY: CLOSE - unable to close NetCDF file: ', ncname]);
end

return

