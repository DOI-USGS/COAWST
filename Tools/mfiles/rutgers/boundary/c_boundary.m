function [status]=c_boundary(S)

%
% C_BOUNDARY:  Create a ROMS lateral boundary conditions NetCDF file
%
% [status]=c_boundary(S)
%
% This function creates a ROMS boundary conditions NetCDF file using
% specified parameters in structure array, S.
%
% On Input:
%
%    S           Boundary condidions creation parameters (structure array):
%
%                  S.ncname           NetCDF file name
%                  S.boundary(:)      Boundary edge switch to process
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

% svn $Id: c_boundary.m 996 2020-01-10 04:28:56Z arango $
%=========================================================================%
%  Copyright (c) 2002-2020 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%

%--------------------------------------------------------------------------
%  Set some NetCDF parameters.
%--------------------------------------------------------------------------

 vartype=nc_constant('nc_float');         % single precision
%vartype=nc_constant('nc_double');        % double precision

%--------------------------------------------------------------------------
%  Get lateral boundary condition creation parameters.
%--------------------------------------------------------------------------

if (isfield(S,'ncname'))
  ncname=S.ncname;
else
  error(['C_BOUNDARY - Cannot find dimension parameter: ncname, ',      ...
         'in structure array S']);
end

if (isfield(S,'spherical'))
  spherical=S.spherical;
else
  spherical=0;
end

if (isfield(S,'Vtransform'))
  Vtransform=S.Vtransform;
else
  error(['C_BOUNDARY - Cannot find dimension parameter: Vtransform, ',  ...
         'in structure array S']);
end

if (isfield(S,'Lm'))
  Lp=S.Lm+2;
else
  error(['C_BOUNDARY - Cannot find dimension parameter: Lm, ',          ...
         'in structure array S']);
end

if (isfield(S,'Mm'))
  Mp=S.Mm+2;
else
  error(['C_BOUNDARY - Cannot find dimension parameter: Mm, ',          ...
         'in structure array S']);
end

if (isfield(S,'N'))
  N=S.N;
else
  error(['C_BOUNDARY - Cannot find dimension parameter: N, ',           ...
         'in structure array S']);
end

if (isfield(S,'NT'))
  NT=S.NT;
else
  error(['C_BOUNDARY - Cannot find dimension parameter: NT, ',          ...
         'in structure S']);
end

%--------------------------------------------------------------------------
%  Set dimensions.
%--------------------------------------------------------------------------

Dname.xr   = 'xi_rho';       Dsize.xr   = Lp;
Dname.xu   = 'xi_u';         Dsize.xu   = Lp-1;
Dname.xv   = 'xi_v';         Dsize.xv   = Lp;
Dname.xp   = 'xi_psi';       Dsize.xp   = Lp-1;
Dname.yr   = 'eta_rho';      Dsize.yr   = Mp;
Dname.yu   = 'eta_u';        Dsize.yu   = Mp;
Dname.yv   = 'eta_v';        Dsize.yv   = Mp-1;
Dname.yp   = 'eta_psi';      Dsize.yp   = Mp-1;
Dname.Nr   = 's_rho';        Dsize.Nr   = N;
Dname.Nw   = 's_w';          Dsize.Nw   = N+1;
Dname.NT   = 'tracer';       Dsize.NT   = NT;
Dname.time = 'bry_time';     Dsize.time = nc_constant('nc_unlimited');

%--------------------------------------------------------------------------
%  Set Variables.
%--------------------------------------------------------------------------

Vname.spherical   = 'spherical';

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

%  Open boundary conditions variables.

Vname.time        = 'bry_time';
Vname.zeta        = 'zeta';
Vname.ubar        = 'ubar';
Vname.vbar        = 'vbar';
Vname.u           = 'u';
Vname.v           = 'v';
Vname.temp        = 'temp';
Vname.salt        = 'salt';

%  Set boundary strings.

Bname = {'west','east','south','north'};

Blong = {'western boundary condition',                                  ...
         'eastern boundary condition',                                  ...
         'southern boundary condition',                                 ...
         'northern boundary condition'};

if (spherical)
  BcoorXr = {'lon_rho_west','lon_rho_east','lon_rho_south','lon_rho_north'};
  BcoorYr = {'lat_rho_west','lat_rho_east','lat_rho_south','lat_rho_north'};
  BcoorXu = {'lon_u_west','lon_u_east','lon_u_south','lon_u_north'};
  BcoorYu = {'lat_u_west','lat_u_east','lat_u_south','lat_u_north'};
  BcoorXv = {'lon_v_west','lon_v_east','lon_v_south','lon_v_north'};
  BcoorYv = {'lat_v_west','lat_v_east','lat_v_south','lat_v_north'};
else
  BcoorXr = {'x_rho_west','x_rho_east','x_rho_south','x_rho_north'};
  BcoorYr = {'y_rho_west','y_rho_east','y_rho_south','y_rho_north'};
  BcoorXu = {'x_u_west','x_u_east','x_u_south','x_u_north'};
  BcoorYu = {'y_u_west','y_u_east','y_u_south','y_u_north'};
  BcoorXv = {'x_v_west','x_v_east','x_v_south','x_v_north'};
  BcoorYv = {'y_v_west','y_v_east','y_v_south','y_v_north'};
end

%--------------------------------------------------------------------------
%  Create open boundary conditions NetCDF file.
%--------------------------------------------------------------------------

[ncid,status]=mexnc('create',ncname,'clobber');
if (status ~= 0)
  disp('  ');
  disp(mexnc('strerror',status));
  error(['C_BOUNDARY: CREATE - unable to create file: ', ncname]);
end

%--------------------------------------------------------------------------
%  Define dimensions.
%--------------------------------------------------------------------------

[did.xr,status]=mexnc('def_dim',ncid,Dname.xr,Dsize.xr); 
if (status ~= 0)
  disp('  ');
  disp(mexnc('strerror',status));
  error([ 'C_BOUNDARY: ncdimdef - unable to define dimension: ',Dname.xr]);
end

[did.xu,status]=mexnc('def_dim',ncid,Dname.xu,Dsize.xu);
if (status ~= 0)
  disp('  ');
  disp(mexnc('strerror',status));
  error(['C_BOUNDARY: DEF_DIM - unable to define dimension: ',Dname.xu]);
end

[did.xv,status]=mexnc('def_dim',ncid,Dname.xv,Dsize.xv);
if (status ~= 0)
  disp('  ');
  disp(mexnc('strerror',status));
  error(['C_BOUNDARY: DEF_DIM - unable to define dimension: ',Dname.xv]);
end

[did.yr,status]=mexnc('def_dim',ncid,Dname.yr,Dsize.yr);
if (status ~= 0)
  disp('  ');
  disp(mexnc('strerror',status));
  error(['C_BOUNDARY: DEF_DIM - unable to define dimension: ',Dname.yr]);
end

[did.yu,status]=mexnc('def_dim',ncid,Dname.yu,Dsize.yu);
if (status ~= 0)
  disp('  ');
  disp(mexnc('strerror',status));
  error(['C_BOUNDARY: DEF_DIM - unable to define dimension: ',Dname.yu]);
end

[did.yv,status]=mexnc('def_dim',ncid,Dname.yv,Dsize.yv);
if (status ~= 0)
  disp('  ');
  disp(mexnc('strerror',status));
  error(['C_BOUNDARY: DEF_DIM - unable to define dimension: ',Dname.yv]);
end

[did.Nr,status]=mexnc('def_dim',ncid,Dname.Nr,Dsize.Nr);
if (status ~= 0)
  disp('  ');
  disp(mexnc('strerror',status));
  error(['C_BOUNDARY: DEF_DIM - unable to define dimension: ',Dname.Nr]);
end

[did.Nw,status]=mexnc('def_dim',ncid,Dname.Nw,Dsize.Nw);
if (status ~= 0)
  disp('  ');
  disp(mexnc('strerror',status));
  error(['C_BOUNDARY: DEF_DIM - unable to define dimension: ',Dname.Nw]);
end

[did.NT,status]=mexnc('def_dim',ncid,Dname.NT,Dsize.NT);
if (status ~= 0)
  disp('  ');
  disp(mexnc('strerror',status));
  error(['C_BOUNDARY: DEF_DIM - unable to define dimension: ',Dname.NT]);
end

[did.time,status]=mexnc('def_dim',ncid,Dname.time,Dsize.time);
if (status ~= 0)
  disp('  ');
  disp(mexnc('strerror',status));
  error(['C_BOUNDARY: DEF_DIM - unable to define dimension: ',Dname.time]);
end

%  Set dimensions for each boundary segment.

did.br=[did.yr did.yr did.xr did.xr];
did.bu=[did.yu did.yu did.xu did.xu];
did.bv=[did.yv did.yv did.xv did.xv];

%--------------------------------------------------------------------------
%  Create global attributes.
%--------------------------------------------------------------------------

type='BOUNDARY file';
lstr=max(size(type));
[status]=mexnc('PUT_ATT_TEXT',ncid,nc_constant('nc_global'),            ...
               'type',nc_constant('nc_char'),lstr,type);
if (status ~= 0)
  disp('  ');
  disp(mexnc('strerror',status));
  error('C_BOUNDARY: PUT_ATT_TEXT - unable to global attribure: type.');
end

history=['Boundary file using Matlab script: c_boundary, ',date_stamp];
lstr=max(size(history));
[status]=mexnc('put_att_text',ncid,nc_constant('nc_global'),            ...
               'history',nc_constant('nc_char'),lstr,history);
if (status ~= 0)
  disp('  ');
  disp(mexnc('strerror',status));
  error('C_BOUNDARY: PUT_ATT_TEXT - unable to global attribure: history.');
end

%--------------------------------------------------------------------------
%  Define configuration variables.
%--------------------------------------------------------------------------

% Define spherical switch.

Var.name                = Vname.spherical;
Var.type                = nc_constant('nc_int');
Var.dimid               = [];
Var.long_name           = 'grid type logical switch';
Var.flag_values         = [0 1];
Var.flag_meanings       = ['Cartesian', blanks(1), ...
                           'spherical'];
[~,status]=nc_vdef(ncid,Var);
if (status ~= 0), return, end
clear Var

% Define vertical coordinate variables.

Var.name                = Vname.Vtransform;
Var.type                = nc_constant('nc_int');
Var.dimid               = [];
Var.long_name           = 'vertical terrain-following transformation equation';
[~,status]=nc_vdef(ncid,Var);
if (status ~= 0), return, end
clear Var

Var.name                = Vname.Vstretching;
Var.type                = nc_constant('nc_int');
Var.dimid               = [];
Var.long_name           = 'vertical terrain-following stretching function';
[~,status]=nc_vdef(ncid,Var);
if (status ~= 0), return, end
clear Var

Var.name                = Vname.theta_s;
Var.type                = nc_constant('nc_double');
Var.dimid               = [];
Var.long_name           = 'S-coordinate surface control parameter';
[~,status]=nc_vdef(ncid,Var);
if (status ~= 0), return, end
clear Var

Var.name                = Vname.theta_b;
Var.type                = nc_constant('nc_double');
Var.dimid               = [];
Var.long_name           = 'S-coordinate bottom control parameter';
[~,status]=nc_vdef(ncid,Var);
if (status ~= 0), return, end
clear Var

Var.name                = Vname.Tcline;
Var.type                = nc_constant('nc_double');
Var.dimid               = [];
Var.long_name           = 'S-coordinate surface/bottom layer width';
Var.units               = 'meter';
[~,status]=nc_vdef(ncid,Var);
if (status ~= 0), return, end
clear Var

Var.name                = Vname.hc;
Var.type                = nc_constant('nc_double');
Var.dimid               = [];
Var.long_name           = 'S-coordinate parameter, critical depth';
Var.units               = 'meter';
[~,status]=nc_vdef(ncid,Var);
if (status ~= 0), return, end
clear Var

Var.name                = Vname.s_rho;
Var.type                = nc_constant('nc_double');
Var.dimid               = [did.Nr];
Var.long_name           = 'S-coordinate at RHO-points';
Var.valid_min           = -1;
Var.valid_max           = 0;
Var.positive            = 'up';
if (Vtransform == 1)
  Var.standard_name     = 'ocean_s_coordinate_g1';
elseif (Vtransform == 2)
  Var.standard_name     = 'ocean_s_coordinate_g2';
end
Var.formula_terms       = 's: s_rho C: Cs_r eta: zeta depth: h depth_c: hc';
[~,status]=nc_vdef(ncid,Var);
if (status ~= 0), return, end
clear Var

Var.name                = Vname.s_w;
Var.type                = nc_constant('nc_double');
Var.dimid               = [did.Nw];
Var.long_name           = 'S-coordinate at W-points';
Var.valid_min           = -1;
Var.valid_max           = 0;
Var.positive            = 'up';
if (Vtransform == 1)
  Var.standard_name     = 'ocean_s_coordinate_g1';
elseif (Vtransform == 2)
  Var.standard_name     = 'ocean_s_coordinate_g2';
end
Var.formula_terms       = 's: s_w C: Cs_w eta: zeta depth: h depth_c: hc';
[~,status]=nc_vdef(ncid,Var);
if (status ~= 0), return, end
clear Var

Var.name                = Vname.Cs_r;
Var.type                = nc_constant('nc_double');
Var.dimid               = [did.Nr];
Var.long_name           = 'S-coordinate stretching function at RHO-points';
Var.valid_min           = -1;
Var.valid_max           = 0;
[~,status]=nc_vdef(ncid,Var);
if (status ~= 0), return, end
clear Var

Var.name                = Vname.Cs_w;
Var.type                = nc_constant('nc_double');
Var.dimid               = [did.Nw];
Var.long_name           = 'S-coordinate stretching function at W-points';
Var.valid_min           = -1;
Var.valid_max           = 0;
[~,status]=nc_vdef(ncid,Var);
if (status ~= 0), return, end
clear Var

%  Define horizontal grid variables.

if (spherical)
  for ib=1:4
    if (S.boundary(ib))
      Var.name          = char(BcoorXr(ib));
      Var.type          = vartype;
      Var.dimid         = did.br(ib);
      Var.long_name     = strcat(['longitude of RHO-points, ',          ...
                                  char(Blong(ib))]);
      Var.units         = 'degree_east';
      Var.standard_name = 'longitude';
      [~,status]=nc_vdef(ncid,Var);
      if (status ~= 0), return, end
      clear Var

      Var.name          = char(BcoorYr(ib));
      Var.type          = vartype;
      Var.dimid         = did.br(ib);
      Var.long_name     = strcat(['latitute of RHO-points, ',           ...
                                  char(Blong(ib))]);
      Var.units         = 'degree_north';
      Var.standard_name = 'latitude';
      [~,status]=nc_vdef(ncid,Var);
      if (status ~= 0), return, end
      clear Var
    end
  end
  
  for ib=1:4
    if (S.boundary(ib))
      Var.name          = char(BcoorXu(ib));
      Var.type          = vartype;
      Var.dimid         = did.bu(ib);
      Var.long_name     = strcat(['longitude of U-points, ',            ...
                                  char(Blong(ib))]);
      Var.units         = 'degree_east';
      Var.standard_name = 'longitude';
      [~,status]=nc_vdef(ncid,Var);
      if (status ~= 0), return, end
      clear Var

      Var.name          = char(BcoorYu(ib));
      Var.type          = vartype;
      Var.dimid         = did.bu(ib);
      Var.long_name     = strcat(['latitute of U-points, ',             ...
                                  char(Blong(ib))]);
      Var.units         = 'degree_north';
      Var.standard_name = 'latitude';
      [~,status]=nc_vdef(ncid,Var);
      if (status ~= 0), return, end
      clear Var
    end
  end
      
  for ib=1:4
    if (S.boundary(ib))
      Var.name          = char(BcoorXv(ib));
      Var.type          = vartype;
      Var.dimid         = did.bv(ib);
      Var.long_name     = strcat(['longitude of V-points, ',            ...
                                  char(Blong(ib))]);
      Var.units         = 'degree_east';
      Var.standard_name = 'longitude';
      [~,status]=nc_vdef(ncid,Var);
      if (status ~= 0), return, end
      clear Var

      Var.name          = char(BcoorYv(ib));
      Var.type          = vartype;
      Var.dimid         = did.bv(ib);
      Var.long_name     = strcat(['latitute of V-points, ',             ...
                                  char(Blong(ib))]);
      Var.units         = 'degree_north';
      Var.standard_name = 'latitude';
      [~,status]=nc_vdef(ncid,Var);
      if (status ~= 0), return, end
      clear Var
    end
  end

else

  for ib=1:4
    if (S.boundary(ib))
      Var.name          = char(BcoorXr(ib));
      Var.type          = vartype;
      Var.dimid         = did.br(ib);
      Var.long_name     = strcat(['X-location of RHO-points, ',         ...
                                  char(Blong(ib))]);
      Var.units         = 'meter';
      [~,status]=nc_vdef(ncid,Var);
      if (status ~= 0), return, end
      clear Var

      Var.name          = char(BcoorYr(ib));
      Var.type          = vartype;
      Var.dimid         = did.br(ib);
      Var.long_name     = strcat(['Y-location of RHO-points, ',         ...
                                 char(Blong(ib))]);
      Var.units         = 'meter';
      [~,status]=nc_vdef(ncid,Var);
      if (status ~= 0), return, end
      clear Var
    end
  end

  for ib=1:4
    if (S.boundary(ib))
      Var.name          = char(BcoorXu(ib));
      Var.type          = vartype;
      Var.dimid         = did.bu(ib);
      Var.long_name     = strcat(['X-location of U-points, ',           ...
                                  char(Blong(ib))]);
      Var.units         = 'meter';
      [~,status]=nc_vdef(ncid,Var);
      if (status ~= 0), return, end
      clear Var

      Var.name          = char(BcoorYu(ib));
      Var.type          = vartype;
      Var.dimid         = did.bu(ib);
      Var.long_name     = strcat(['Y-location of U-points, ',           ...
                                  char(Blong(ib))]);
      Var.units         = 'meter';
      [~,status]=nc_vdef(ncid,Var);
      if (status ~= 0), return, end
      clear Var
    end
  end
  
  for ib=1:4
    if (S.boundary(ib))
      Var.name          = char(BcoorXv(ib));
      Var.type          = vartype;
      Var.dimid         = did.bv(ib);
      Var.long_name     = strcat(['X-location of V-points, ',           ...
                                  char(Blong(ib))]);
      Var.units         = 'meter';
      [~,status]=nc_vdef(ncid,Var);
      if (status ~= 0), return, end
      clear Var

      Var.name          = char(BcoorYv(ib));
      Var.type          = vartype;
      Var.dimid         = did.bv(ib);
      Var.long_name     = strcat(['Y-location of V-points, ',           ...
                                  char(Blong(ib))]);
      Var.units         = 'meter';
      [~,status]=nc_vdef(ncid,Var);
      if (status ~= 0), return, end
      clear Var
    end
  end
  
end

%--------------------------------------------------------------------------
%  Define open boundary conditions variables.
%--------------------------------------------------------------------------

Var.name                = Vname.time;
Var.type                = nc_constant('nc_double');
Var.dimid               = [did.time];
Var.long_name           = 'time since initialization';
Var.units               = 'seconds';
[~,status]=nc_vdef(ncid,Var);
if (status ~= 0), return, end
clear Var

for ib=1:4
  if (S.boundary(ib))
    Var.name            = strcat([Vname.zeta,'_',char(Bname(ib))]);
    Var.type            = vartype;
    Var.dimid           = [did.time did.br(ib)];
    Var.long_name       = strcat(['free-surface, ',char(Blong(ib))]);
    Var.units           = 'meter';
    Var.time            = Vname.time;
    Var.coordinates     = strcat([char(BcoorXr(ib)),' ',                ...
                                  char(BcoorYr(ib)),' ',                ...
                                  Vname.time]);
    [~,status]=nc_vdef(ncid,Var);
    if (status ~= 0), return, end
    clear Var
  end
end

for ib=1:4
  if (S.boundary(ib))
    Var.name            = strcat([Vname.ubar,'_',char(Bname(ib))]);
    Var.type            = vartype;
    Var.dimid           = [did.time did.bu(ib)];
    Var.long_name       = strcat(['vertically integrated u-momentum ',  ...
                                  'component, ',char(Blong(ib))]);
    Var.units           = 'meter second-1';
    Var.time            = Vname.time;
    Var.coordinates     = strcat([char(BcoorXu(ib)),' ',                ...
                                  char(BcoorYu(ib)),' ',                ...
                                  Vname.time]); 
    [~,status]=nc_vdef(ncid,Var);
    if (status ~= 0), return, end
    clear Var
  end
end

for ib=1:4
  if (S.boundary(ib))
    Var.name            = strcat([Vname.vbar,'_',char(Bname(ib))]);
    Var.type            = vartype;
    Var.dimid           = [did.time did.bv(ib)];
    Var.long_name       = strcat(['vertically integrated v-momentum ',  ...
                                  'component, ',char(Blong(ib))]);
    Var.units           = 'meter second-1';
    Var.time            = Vname.time;
    Var.coordinates     = strcat([char(BcoorXv(ib)),' ',                ...
                                  char(BcoorYv(ib)),' ',                ...
                                  Vname.time]); 
    [~,status]=nc_vdef(ncid,Var);
    if (status ~= 0), return, end
    clear Var
  end
end

for ib=1:4
  if (S.boundary(ib))
    Var.name            = strcat([Vname.u,'_',char(Bname(ib))]);
    Var.type            = vartype;
    Var.dimid           = [did.time did.Nr did.bu(ib)];
    Var.long_name       = strcat(['u-momentum component, ',             ...
                                  char(Blong(ib))]);
    Var.units           = 'meter second-1';
    Var.time            = Vname.time;
    Var.coordinates     = strcat([char(BcoorXu(ib)),' ',                ...
                                  char(BcoorYu(ib)),' ',...
                                  Vname.s_rho,' ',Vname.time]); 
    [~,status]=nc_vdef(ncid,Var);
    if (status ~= 0), return, end
    clear Var
  end
end

for ib=1:4
  if (S.boundary(ib))
    Var.name            = strcat([Vname.v,'_',char(Bname(ib))]);
    Var.type            = vartype;
    Var.dimid           = [did.time did.Nr did.bv(ib)];
    Var.long_name       = strcat(['v-momentum component, ',             ...
                                  char(Blong(ib))]);
    Var.units           = 'meter second-1';
    Var.time            = Vname.time;
    Var.coordinates     = strcat([char(BcoorXv(ib)),' ',                ...
                                  char(BcoorYv(ib)),' ',                ...
                                  Vname.s_rho,' ',Vname.time]); 
    [~,status]=nc_vdef(ncid,Var);
    if (status ~= 0), return, end
    clear Var
  end
end

for ib=1:4
  if (S.boundary(ib))
    Var.name            = strcat([Vname.temp,'_',char(Bname(ib))]);
    Var.type            = vartype;
    Var.dimid           = [did.time did.Nr did.br(ib)];
    Var.long_name       = strcat(['potential temperature, ',            ...
                                  char(Blong(ib))]);
    Var.units           = 'Celsius';
    Var.time            = Vname.time;
    Var.coordinates     = strcat([char(BcoorXr(ib)),' ',                ...
                                  char(BcoorYr(ib)),' ',                ...
                                  Vname.s_rho,' ',Vname.time]); 
    [~,status]=nc_vdef(ncid,Var);
    if (status ~= 0), return, end
    clear Var
  end
end

for ib=1:4
  if (S.boundary(ib))
    Var.name            = strcat([Vname.salt,'_',char(Bname(ib))]);
    Var.type            = vartype;
    Var.dimid           = [did.time did.Nr did.br(ib)];
    Var.long_name       = strcat(['salinity, ',char(Blong(ib))]);
    Var.time            = Vname.time;
    Var.coordinates     = strcat([char(BcoorXr(ib)),' ',                ...
                                  char(BcoorYr(ib)),' ',                ...
                                  Vname.s_rho,' ',Vname.time]); 
    [~,status]=nc_vdef(ncid,Var);
    if (status ~= 0), return, end
    clear Var
  end
end

%--------------------------------------------------------------------------
%  Leave definition mode and close NetCDF file.
%--------------------------------------------------------------------------

[status]=mexnc('enddef',ncid);
if (status == -1)
  disp('  ');
  disp(mexnc('strerror',status));
  error('C_BOUNDARY: ENDDEF - unable to leave definition mode.');
end

[status]=mexnc('close',ncid);
if (status == -1)
  disp('  ');
  disp(mexnc('strerror',status));
  error(['C_BOUNDARY: CLOSE - unable to close NetCDF file: ', ncname]);
end

return

