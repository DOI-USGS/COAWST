%
%  D_CLIMATOLOGY:  Driver script to create a ROMS climatology file.
%
%  This a user modifiable script that can be used to prepare ROMS
%  climatology NetCDF file.  It sets-up all the necessary parameters
%  and variables. USERS can use this as a prototype for their
%  application.
%

% svn $Id: d_climatology.m 996 2020-01-10 04:28:56Z arango $
%=========================================================================%
%  Copyright (c) 2002-2020 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%

%  Set input/output NetCDF files.

 my_root = '~/ocean/repository/Projects/damee';

 GRDname = fullfile(my_root, 'Data/netcdf3', 'damee4_grid_a.nc');
%CLMname = fullfile(my_root, 'Data/netcdf3', 'damee4_Lclm_b.nc');

 CLMname = 'damee4_Lclm_b.nc';

 OAname( 1,:)  = fullfile(my_root, 'Data/OA', 'oa4_lev94_jan.nc');
 OAname( 2,:)  = fullfile(my_root, 'Data/OA', 'oa4_lev94_feb.nc');
 OAname( 3,:)  = fullfile(my_root, 'Data/OA', 'oa4_lev94_mar.nc');
 OAname( 4,:)  = fullfile(my_root, 'Data/OA', 'oa4_lev94_apr.nc');
 OAname( 5,:)  = fullfile(my_root, 'Data/OA', 'oa4_lev94_may.nc');
 OAname( 6,:)  = fullfile(my_root, 'Data/OA', 'oa4_lev94_jun.nc');
 OAname( 7,:)  = fullfile(my_root, 'Data/OA', 'oa4_lev94_jul.nc');
 OAname( 8,:)  = fullfile(my_root, 'Data/OA', 'oa4_lev94_aug.nc');
 OAname( 9,:)  = fullfile(my_root, 'Data/OA', 'oa4_lev94_sep.nc');
 OAname(10,:)  = fullfile(my_root, 'Data/OA', 'oa4_lev94_oct.nc');
 OAname(11,:)  = fullfile(my_root, 'Data/OA', 'oa4_lev94_nov.nc');
 OAname(12,:)  = fullfile(my_root, 'Data/OA', 'oa4_lev94_dec.nc');

%  Set local variables.

 OA_INTERPOLATE = 1;               % switch to interpolate from OA fields
%method='linear';                  % linear interpolation
 method='spline';                  % spline interpolation

%--------------------------------------------------------------------------
%  Set application parameters in structure array, S.
%--------------------------------------------------------------------------

%  Climatology output file name.

S.ncname = CLMname;

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

%  Set climatology fields to process and the number of time records.
%    (S.def_* = 0,  do not process
%     S.def_* = 1,  process)

S.def_zeta = 0;       S.zeta_time = 12;
S.def_v2d  = 0;       S.v2d_time  = 12;
S.def_v3d  = 0;       S.v3d_time  = 12;
S.def_temp = 1;       S.temp_time = 12;
S.def_salt = 1;       S.salt_time = 12;

%--------------------------------------------------------------------------
%  Create climatology Netcdf file.
%--------------------------------------------------------------------------

[~]=c_climatology(S);

%  Update attributes for time variables.

if (S.def_zeta),
  avalue='days since 0001-01-01 00:00:00';
  [~]=nc_attadd(CLMname,'units',avalue,'zeta_time');
  
  avalue='360.0 days in every year';
  [~]=nc_attadd(CLMname,'calendar',avalue,'zeta_time');

  [~]=nc_attadd(CLMname,'cycle_length',360.0,'zeta_time');
end

if (S.def_v2d),
  avalue='days since 0001-01-01 00:00:00';
  [~]=nc_attadd(CLMname,'units',avalue,'v2d_time');
  
  avalue='360.0 days in every year';
  [~]=nc_attadd(CLMname,'calendar',avalue,'v2d_time');

  [~]=nc_attadd(CLMname,'cycle_length',360.0,'v2d_time');
end

if (S.def_v3d),
  avalue='days since 0001-01-01 00:00:00';
  [~]=nc_attadd(CLMname,'units',avalue,'v3d_time');
  
  avalue='360.0 days in every year';
  [~]=nc_attadd(CLMname,'calendar',avalue,'v3d_time');

  [~]=nc_attadd(CLMname,'cycle_length',360.0,'v3d_time');
end

if (S.def_temp),
  avalue='days since 0001-01-01 00:00:00';
  [~]=nc_attadd(CLMname,'units',avalue,'temp_time');
  
  avalue='360.0 days in every year';
  [~]=nc_attadd(CLMname,'calendar',avalue,'temp_time');

  [~]=nc_attadd(CLMname,'cycle_length',360.0,'temp_time');
end

if (S.def_salt),
  avalue='days since 0001-01-01 00:00:00';
  [~]=nc_attadd(CLMname,'units',avalue,'salt_time');
  
  avalue='360.0 days in every year';
  [~]=nc_attadd(CLMname,'calendar',avalue,'salt_time');

  [~]=nc_attadd(CLMname,'cycle_length',360.0,'salt_time');
end

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
  end
end


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
%  If Land/Sea masking arrays are not found, initialize them to unity.
%--------------------------------------------------------------------------

Lr = S.Lm+2;   Lu = Lr-1;   Lv = Lr;
Mr = S.Mm+2;   Mu = Mr;     Mv = Mr-1;

if (~isfield(S, 'mask_rho')),  S.mask_rho = ones([Lr Mr]);  end,
if (~isfield(S, 'mask_u'  )),  S.mask_u   = ones([Lu Mu]);  end,
if (~isfield(S, 'mask_v'  )),  S.mask_v   = ones([Lv Mv]);  end,

%--------------------------------------------------------------------------
%  Write out grid variables.
%--------------------------------------------------------------------------
			 
[~]=nc_write(CLMname,   'spherical',   S.spherical);

[~]=nc_write(CLMname,   'Vtransform',  S.Vtransform);
[~]=nc_write(CLMname,   'Vstretching', S.Vstretching);
[~]=nc_write(CLMname,   'theta_s',     S.theta_s);
[~]=nc_write(CLMname,   'theta_b',     S.theta_b);
[~]=nc_write(CLMname,   'Tcline',      S.Tcline);
[~]=nc_write(CLMname,   'hc',          S.hc);

[~]=nc_write(CLMname,   's_rho',       S.s_rho);
[~]=nc_write(CLMname,   's_w',         S.s_w);
[~]=nc_write(CLMname,   'Cs_r',        S.Cs_r);
[~]=nc_write(CLMname,   'Cs_w',        S.Cs_w);

[~]=nc_write(CLMname,   'h',           S.h);

if (S.spherical),
  [~]=nc_write(CLMname, 'lon_rho',     S.lon_rho);
  [~]=nc_write(CLMname, 'lat_rho',     S.lat_rho);
  [~]=nc_write(CLMname, 'lon_u',       S.lon_u);
  [~]=nc_write(CLMname, 'lat_u',       S.lat_u);
  [~]=nc_write(CLMname, 'lon_v',       S.lon_v);
  [~]=nc_write(CLMname, 'lat_v',       S.lat_v);
else
  [~]=nc_write(CLMname, 'x_rho',       S.x_rho);
  [~]=nc_write(CLMname, 'y_rho',       S.y_rho);
  [~]=nc_write(CLMname, 'x_u',         S.x_u);
  [~]=nc_write(CLMname, 'y_u',         S.y_u);
  [~]=nc_write(CLMname, 'x_v',         S.x_v);
  [~]=nc_write(CLMname, 'y_v',         S.y_v);
end

%--------------------------------------------------------------------------
%  Compute depths at horizontal and vertical RHO-points.
%--------------------------------------------------------------------------

igrid = 1;

S.zeta = zeros(size(S.h));

[z_r] = set_depth(S.Vtransform, S.Vstretching,                          ...
                  S.theta_s, S.theta_b, S.hc, S.N,                      ...
                 igrid, S.h, S.zeta);

%--------------------------------------------------------------------------
%  Interpolate OA of temperature and salinity from standard levels to
%  model depths.
%--------------------------------------------------------------------------

if (OA_INTERPOLATE),

  Nrec=size(OAname,1);
  
  disp(' ')
  disp('Interpolating from OA fields, please wait ...');
  disp(' ')
  
  for rec=1:Nrec,

    [dir,name]=fileparts(OAname(rec,:));
    
    disp([ '  Procesing OA file: ', name]);

    Zoa=nc_read(OAname(rec,:), 'zout');

    oa_temp=nc_read(OAname(rec,:), 'temp', 1);
    oa_salt=nc_read(OAname(rec,:), 'salt', 1);

    S.temp_time(rec)=nc_read(OAname(rec,:), 'time', 1);
    S.salt_time(rec)=nc_read(OAname(rec,:), 'time', 1);
    
    for j=1:Mr,
      for i=1:Lr,
        Zroms = squeeze(z_r(i,j,:));
        Toa   = squeeze(oa_temp(i,j,:));
        Soa   = squeeze(oa_salt(i,j,:));

        S.temp(i,j,:,rec) = interp1(Zoa, Toa, Zroms, method);
        S.salt(i,j,:,rec) = interp1(Zoa, Soa, Zroms, method);
      end
    end

  end
  
end

%--------------------------------------------------------------------------
%  Write out climatology.
%--------------------------------------------------------------------------

for rec=1:Nrec,
  
  if (S.def_zeta),
    [~]=nc_write(CLMname, 'zeta_time', S.zeta_time(rec), rec);
    [~]=nc_write(CLMname, 'zeta', S.zeta(:,:,rec), rec);
  end

  if (S.def_v2d),
    [~]=nc_write(CLMname, 'v2d_time', S.v2d_time(rec), rec);
    [~]=nc_write(CLMname, 'ubar', S.ubar(:,:,rec), rec);
    [~]=nc_write(CLMname, 'vbar', S.vbar(:,:,rec), rec);
  end

  if (S.def_v3d),
    [~]=nc_write(CLMname, 'v3d_time', S.v3d_time(rec), rec);
    [~]=nc_write(CLMname, 'u', S.u(:,:,:,rec), rec);
    [~]=nc_write(CLMname, 'v', S.v(:,:,:,rec), rec);
  end

  if (S.def_temp),
    [~]=nc_write(CLMname, 'temp_time', S.temp_time(rec), rec);
    [~]=nc_write(CLMname, 'temp', S.temp(:,:,:,rec), rec);
  end

  if (S.def_salt),
    [~]=nc_write(CLMname, 'salt_time', S.salt_time(rec), rec);
    [~]=nc_write(CLMname, 'salt', S.salt(:,:,:,rec), rec);
  end

end
