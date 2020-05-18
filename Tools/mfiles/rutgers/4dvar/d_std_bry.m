%
%  D_STD_BRY:  Driver script to compute and write 4D-Var open boundary
%              conditions standard deviations.
%
%  This a user modifiable script that can be used to prepare ROMS 4D-Var
%  open boundary standard deviations which are used to convert modeled
%  error correlations to error covariances.
%
%  There are many ways to compute the open boundary conditions standard
%  deviations.  It depends on the source of the lateral boundary
%  conditions data.
% 
%  As an example, we are extracting the standard deviations from an
%  existing initial conditions standard deviation file. That is, a
%  full grid NetCDF file.  This not a good way to do it but it is
%  done here to show the structure of such computation.
%

% svn $Id: d_std_bry.m 996 2020-01-10 04:28:56Z arango $
%=========================================================================%
%  Copyright (c) 2002-2020 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%

%--------------------------------------------------------------------------
%  User tunable parameters.
%--------------------------------------------------------------------------

%  Set output standard deviation NetCDF file.

STDfile = 'wc13_std_b.nc';

%  Set input grid file.

my_root = '/home/arango/ocean/toms/repository/test';

GRDfile = strcat(my_root, '/WC13/Data/wc13_grd.nc');

%  Set input initial conditions standard deviation.

my_root = '/home/arango/ocean/toms/repository/test';

INIfile = strcat(my_root, '/WC13/Data/wc13_std_i.nc');

%  Set state variables dynamical fields (cell array) to process.

field_list = {'zeta', 'ubar', 'vbar', 'u', 'v', 'temp', 'salt'};

%  Set grid variables dynamical fields (cell array) to write in
%  output file(s).

grid_list  = {'theta_s', 'theta_b', 'Tcline' , 'hc'     , 's_rho'  ,    ...
              's_w'    , 'Cs_r'   , 'Cs_w'   , 'h'      , 'lon_rho',    ...
              'lat_rho', 'lon_u'  , 'lat_u'  , 'lon_v'  , 'lat_v'};

%  Initialize working structure.

S.title = 'California Current System, 1/3 degree resolution (WC13)';

S.ncname = STDfile;                 % output file

S.grd_file = GRDfile;               % application grid

S.do_zeta = true;                   % free-surface
S.do_ubar = true;                   % vertically integrated u-momentum
S.do_vbar = true;                   % vertically integrated v-momentum
S.do_u    = true;                   % u-momentum
S.do_v    = true;                   % v-momentum
S.do_temp = true;                   % temperature
S.do_salt = true;                   % salinity

%--------------------------------------------------------------------------
%  Extract open boundary conditions standard deviations from initial
%  conditions standard deviation file.
%--------------------------------------------------------------------------

%  Set boundary edges indices.

iwest  = 1;
isouth = 2;
ieast  = 3;
inorth = 4;

%  Set NetCDF file creation parameters.

V = nc_vnames(INIfile);
nvars = length(V.Variables);
 
S.curvilinear = false;
S.masking     = false;
S.Vtransform  = 1;
S.Vstretching = 1;

for n=1:nvars,
  name = char(V.Variables(n).Name);
  switch name
    case 'spherical'
      S.spherical = nc_read(INIfile, 'spherical');
      if (ischar(S.spherical)),
        if (S.spherical == 'T' || S.spherical == 't');
          S.spherical = 1;
        else
          S.spherical = 0;
        end
      end
    case 'Vtransform'
      S.Vtransform  = nc_read(INIfile, 'Vtransform');
    case 'Vstretching'
      S.Vstretching = nc_read(INIfile, 'Vstretching');
    case 'angle'
      S.angle = nc_read(INIfile, 'angle');
      S.curvilinear = true;
    case 'mask_rho'
      S.rmask = nc_read(INIfile, 'mask_rho');
      S.umask = nc_read(INIfile, 'mask_u');
      S.vmask = nc_read(INIfile, 'mask_v');
      S.masking = true;
  end
end

if (S.curvilinear),
  grid_list = [grid_list, 'angle'];
end
  
if (S.masking),
  grid_list = [grid_list, 'mask_rho', 'mask_u', 'mask_v'];
end

%  Get grid size.

[Lr,Mr,Nr] = size(nc_read(INIfile, 'temp', 1));

S.Lm = Lr - 2;                      % number of interior RHO x-points
S.Mm = Mr - 2;                      % number of interior RHO y-points
S.N  = Nr;                          % number of vertical RHO levels

IorJ = max(Lr, Mr);                 % maximum number of lateral points

%  Read in grid.

S.rlon = nc_read(INIfile, 'lon_rho');
S.rlat = nc_read(INIfile, 'lat_rho');

S.ulon = nc_read(INIfile, 'lon_u');
S.ulat = nc_read(INIfile, 'lat_u');

S.vlon = nc_read(INIfile, 'lon_v');
S.vlat = nc_read(INIfile, 'lat_v');

%  Extract open boundary conditions standard deviations for the initial
%  time.  This is just a template for user to see how the boundary arrays
%  are built and not a good meassure of the open boundady conditions
%  standard deviation. Notice that all the boundary arrays are built
%  in a compact way:
%
%      bry2d(IorJ, 4)           2D field
%      bry3d(IorJ, N, 4)        3D field
%
%  The IorJ dimension is the maximum of the size in the X- or Y-directions
%  at RHO-points. The U-points variables need to be shifted by one in the
%  X-direction and the V-points variables need to be shifted by one in the
%  Y-direction (see "extract_bry.m").  This is because of ROMS staggered
%  C-grid variables. In Fortran,
%
%      RHO-variable:            0:Lm+1, 0:Mm+1
%        U-variable:            1:Lm+1, 0:Mm+1
%        V-variable:            0:Lm+1, 1:Mm+1
%
%  It is up the user to compute the appropriate value of the standard
%  deviation for the open boundary conditions. The example here is to
%  show how the boundary arrays are built and written.

rec = 1;                                  % record to process

compact = true;                           % extract compact boundary data

for fval = field_list,
  field     = char(fval);                 % convert cell to string
  field_std = [field, '_std'];            % standard deviation field

  S.(field_std) = extract_bry(INIfile, field, rec, compact);
end

%--------------------------------------------------------------------------
%  Write out standard deviation fields.
%--------------------------------------------------------------------------

disp(' ');
disp([ 'Writting out standard deviation, file = ', S.ncname]);
disp(' ');

%  Create standard deviations NetCDF file.

s = c_std_bry(S);

%  Write out grid data.

v = 'spherical';    s = nc_write(S.ncname, v, S.spherical);
v = 'Vtransform';   s = nc_write(S.ncname, v, S.Vtransform);
v = 'Vstretching';  s = nc_write(S.ncname, v, S.Vstretching);

for fval = grid_list,
  field = char(fval);                     % convert cell to string

  f = nc_read(INIfile, field);  s = nc_write(S.ncname, field, f);
end

% Write out standard deviation data.

rec = 1;
  
s = nc_write(S.ncname, 'ocean_time', 0, rec);

for fval = field_list,
  field     = char(fval);                 % convert cell to string
  field_obc = [field, '_obc'];            % open boundary condition field
  field_std = [field, '_std'];            % standard deviation field

  s = nc_write(S.ncname, field_obc, S.(field_std), rec);
end

disp(' ');
disp('Done.');
disp(' ');

