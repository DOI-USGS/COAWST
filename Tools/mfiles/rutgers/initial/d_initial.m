%
%  D_INITIAL:  Driver script to create a ROMS initial conditions file.
%
%  This a user modifiable script that can be used to prepare ROMS initial
%  conditions NetCDF file.  It sets-up all the necessary parameters and
%  variables. USERS can use this as a prototype for their application.
%

% svn $Id: d_initial.m 996 2020-01-10 04:28:56Z arango $
%=========================================================================%
%  Copyright (c) 2002-2020 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%

%  Set input/output NetCDF files.

 my_root = '~/ocean/repository/Projects/damee';

 GRDname = fullfile(my_root, 'Data/netcdf3', 'damee4_grid_a.nc');
 OAname  = fullfile(my_root, 'Data/OA',      'oa4_lev94_feb.nc');
%INIname = fullfile(my_root, 'Data/netcdf3', 'damee4_levfeb_b.nc');

 INIname = 'damee4_levfeb_b.nc';

%  Set local variables.

 OA_INTERPOLATE = 1;               % switch to interpolate from OA fields
%method='linear';                  % linear interpolation
 method='spline';                  % spline interpolation

%--------------------------------------------------------------------------
%  Set application parameters in structure array, S.
%--------------------------------------------------------------------------

%  Initial conditions output file name.

S.ncname = INIname;

%  Spherical grid switch.
%
%            [0] Cartesian grid
%            [1] Spherical grid

S.spherical = 1;

%  Number of interior RHO-points (Lm, Mm) in the X- and Y-directions. These
%  are the values specified in ROMS input script.

S.Lm = 128;
S.Mm = 128;

%  Number of vertical levels (N) at RHO-points (center grid cell).

S.N = 20;

%  Number of total active and passive tracers. Usually, NT=2 for potential
%  temperature and salinity.

S.NT = 2;

%  Set vertical stretching parameters. These values are specified in ROMS
%  input script.

%  Vertical transfomation equation (Vtransform).
%    (Check  https://www.myroms.org/wiki/index.php/Vertical_S-coordinate)
%
%            [1] Original tranformation equation
%            [2] New transformation equation

S.Vtransform = 2;

%  Vertical stretching function (Vstretching).
%    (Check  https://www.myroms.org/wiki/index.php/Vertical_S-coordinate)
%
%            [1] Original function (Song and Haidvogel, 1994)
%            [2] A. Shchepetkin function
%            [3] R. Geyer function

S.Vstretching = 2;

%  Vertical stretching parameters:
%
%            tetha_s:  S-coordinate surface control parameter.
%            tetha_b:  S-coordinate bottom control parameter.
%            Tcline:   S-coordinate surface/bottom stretching width (m)

S.theta_s = 7.0;
S.theta_b = 0.1;
S.Tcline  = 200.0;

S.hc = S.Tcline;

%--------------------------------------------------------------------------
%  Create initial condition Netcdf file.
%--------------------------------------------------------------------------

[~]=c_initial(S);

%  Set attributes for "ocean_time".

avalue='seconds since 0001-01-01 00:00:00';
[~]=nc_attadd(INIname,'units',avalue,'ocean_time');
  
avalue='360.0 days in every year';
[~]=nc_attadd(INIname,'calendar',avalue,'ocean_time');

%--------------------------------------------------------------------------
%  Set grid variables.
%--------------------------------------------------------------------------

V=nc_vnames(GRDname);
nvars=length(V.Variables);

%  Horizontal grid variables. Read in for input GRID NetCDF file.

if (S.spherical),
  S.lon_rho = nc_read(GRDname, 'lon_rho');
  S.lat_rho = nc_read(GRDname, 'lat_rho');
  
  S.lon_u   = nc_read(GRDname, 'lon_u');
  S.lat_u   = nc_read(GRDname, 'lat_u');
  
  S.lon_v   = nc_read(GRDname, 'lon_v');
  S.lat_v   = nc_read(GRDname, 'lat_v');
else  
  S.x_rho   = nc_read(GRDname, 'x_rho');
  S.y_rho   = nc_read(GRDname, 'y_rho');
  
  S.x_u     = nc_read(GRDname, 'x_u');
  S.y_u     = nc_read(GRDname, 'y_u');
  
  S.x_v     = nc_read(GRDname, 'x_v');
  S.y_v     = nc_read(GRDname, 'y_v');  
end  

%  Read in Land/Sea mask, if appropriate.

for n=1:nvars,
  name=char(V.Variables(n).Name);
  switch (name),
    case 'mask_rho'
      S.mask_rho = nc_read(GRDname, 'mask_rho');
    case 'mask_u'
      S.mask_u   = nc_read(GRDname, 'mask_u');
    case 'mask_v'
      S.mask_v   = nc_read(GRDname, 'mask_v');
  end,
end,


%  Bathymetry.

S.h = nc_read(GRDname, 'h');

%  Set vertical grid variables.

[S.s_rho, S.Cs_r]=stretching(S.Vstretching, ...
                             S.theta_s, S.theta_b, S.hc, S.N,           ...
			     0, 1);

[S.s_w,   S.Cs_w]=stretching(S.Vstretching, ...
                             S.theta_s, S.theta_b, S.hc, S.N,           ...
			     1, 1);

%--------------------------------------------------------------------------
%  Set zero initial conditions.
%--------------------------------------------------------------------------

Lr = S.Lm+2;   Lu = Lr-1;   Lv = Lr;
Mr = S.Mm+2;   Mu = Mr;     Mv = Mr-1;

S.zeta = zeros([Lr Mr]);
S.ubar = zeros([Lu Mu]);
S.vbar = zeros([Lv Mv]);
S.u    = zeros([Lu Mu S.N]);
S.v    = zeros([Lv Mv S.N]);
S.temp = zeros([Lr Mr S.N]);
S.salt = zeros([Lr Mr S.N]);

%  If Land/Sea masking arrays are not found, initialize them to unity.

if (~isfield(S, 'mask_rho')),  S.mask_rho = ones([Lr Mr]);  end,
if (~isfield(S, 'mask_u'  )),  S.mask_u   = ones([Lu Mu]);  end,
if (~isfield(S, 'mask_v'  )),  S.mask_v   = ones([Lv Mv]);  end,

%--------------------------------------------------------------------------
%  Write out grid variables.
%--------------------------------------------------------------------------
			 
[~]=nc_write(INIname,   'spherical',   S.spherical);

[~]=nc_write(INIname,   'Vtransform',  S.Vtransform);
[~]=nc_write(INIname,   'Vstretching', S.Vstretching);
[~]=nc_write(INIname,   'theta_s',     S.theta_s);
[~]=nc_write(INIname,   'theta_b',     S.theta_b);
[~]=nc_write(INIname,   'Tcline',      S.Tcline);
[~]=nc_write(INIname,   'hc',          S.hc);

[~]=nc_write(INIname,   's_rho',       S.s_rho);
[~]=nc_write(INIname,   's_w',         S.s_w);
[~]=nc_write(INIname,   'Cs_r',        S.Cs_r);
[~]=nc_write(INIname,   'Cs_w',        S.Cs_w);

[~]=nc_write(INIname,   'h',           S.h);

if (S.spherical),
  [~]=nc_write(INIname, 'lon_rho',     S.lon_rho);
  [~]=nc_write(INIname, 'lat_rho',     S.lat_rho);
  [~]=nc_write(INIname, 'lon_u',       S.lon_u);
  [~]=nc_write(INIname, 'lat_u',       S.lat_u);
  [~]=nc_write(INIname, 'lon_v',       S.lon_v);
  [~]=nc_write(INIname, 'lat_v',       S.lat_v);
else
  [~]=nc_write(INIname, 'x_rho',       S.x_rho);
  [~]=nc_write(INIname, 'y_rho',       S.y_rho);
  [~]=nc_write(INIname, 'x_u',         S.x_u);
  [~]=nc_write(INIname, 'y_u',         S.y_u);
  [~]=nc_write(INIname, 'x_v',         S.x_v);
  [~]=nc_write(INIname, 'y_v',         S.y_v);
end

%--------------------------------------------------------------------------
%  Compute depths at horizontal and vertical RHO-points.
%--------------------------------------------------------------------------

igrid = 1;

[z_r] = set_depth(S.Vtransform, S.Vstretching,                          ...
                  S.theta_s, S.theta_b, S.hc, S.N,                      ...
                  igrid, S.h, S.zeta);

%--------------------------------------------------------------------------
%  Interpolate OA of temperature and salinity from standard levels to
%  model depths.
%--------------------------------------------------------------------------

if (OA_INTERPOLATE),

  disp(' ')
  disp('Interpolating from OA fields, please wait ...');
  
  InpRec = 1;

  Zoa=nc_read(OAname, 'zout');

  oa_temp=nc_read(OAname, 'temp', InpRec);
  oa_salt=nc_read(OAname, 'salt', InpRec);

  for j=1:Mr,
    for i=1:Lr,
      Zroms = squeeze(z_r(i,j,:));
      Toa   = squeeze(oa_temp(i,j,:));
      Soa   = squeeze(oa_salt(i,j,:));

      S.temp(i,j,:) = interp1(Zoa, Toa, Zroms, method);
      S.salt(i,j,:) = interp1(Zoa, Soa, Zroms, method);
    end,
  end,
  
end,

%---------------------------------------------------------------------------
%  Write out initial conditions.
%---------------------------------------------------------------------------

IniRec = 1;                               % NetCDF time record

S.ocean_time = 30.0*86400;                % initial conditions time (s)

[~]=nc_write(INIname, 'ocean_time', S.ocean_time, IniRec);

[~]=nc_write(INIname, 'zeta', S.zeta, IniRec);
[~]=nc_write(INIname, 'ubar', S.ubar, IniRec);
[~]=nc_write(INIname, 'vbar', S.vbar, IniRec);
[~]=nc_write(INIname, 'u',    S.u,    IniRec);
[~]=nc_write(INIname, 'v',    S.v,    IniRec);
[~]=nc_write(INIname, 'temp', S.temp, IniRec);
[~]=nc_write(INIname, 'salt', S.salt, IniRec);
