%
%  D_NUDGCOEF:  Driver script to create a ROMS nudging coefficients file.
%
%  This a user modifiable script that can be used to prepare ROMS
%  nudging inverse time scales NetCDF file.  It sets-up  all the
%  necessary parameters and variables. USERS can use this as a
%  prototype for their application.
%
%  Nudging to climatology can be used in ROMS for various purposes:
%
%  (1) Improve the behavior of open boundary conditions.
%  (2) Used in conjunction with sponges.
%  (3) Minimize numerical diapycnal mixing of tracers over steep
%      bathymetry (improve T-S properties in deep water masses).  For
%      example, we can nudge to T-S climatology is areas depeer than
%      1500 m.
%
%  The inverse nudging coefficients have units of 1/time.  The default
%  input units in ROMS is 1/day but 1/second is also possible. The
%  routine 'get_nudgcoef.F' will check the 'units' attribute to compute
%  the conversion factor for 1/second. Users need to be sure that the
%  'units' variable attribute is consistent with the written data.
%
%  The variable names for the nudging coefficients is as follows:
%
%     M2_NudgeCoef       for 2D momentum
%     M3_NudgeCoef       for 3D momentum
%     temp_NudgeCoef     for potential temperature
%     salt_NudgeCoef     for salinity
%     ...
%     NO3_NudgeCoef      for nitrate
%     ...
%     tracer_NudgeCoef   for any generic tracer
%
%  They are all defined at RHO-points. If the nudging coefficients for
%  a specific tracer are available in the NetCDF, ROMS will read that
%  NetCDF variable. If NOT and the generic coefficients 'tracer_NudgeCoef'
%  are available, ROMS will process those values instead.
%
%  Notice that the input swicth 'LnudgeTCLM(itrc,ng)' in ROMS input
%  script 'ocean.in' will control which tracer to nudge in the desired
%  grid.
%
%  Currently, the nudging coefficients are time invariant in ROMS.  The
%  same scales are used for the entire simulation.
%

% svn $Id: d_nudgcoef.m 996 2020-01-10 04:28:56Z arango $
%=========================================================================%
%  Copyright (c) 2002-2020 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%

% Set input/output NetCDF files.

 my_root = '~/ocean/repository/Projects/damee';

 GRDname = fullfile(my_root, 'Data/netcdf3', 'damee4_grid_b.nc');
 INIname = fullfile(my_root, 'Data/netcdf3', 'damee4_levfeb_b.nc');

 NUDname = 'damee4_nudgcoef_b.nc';
 
% Get grid structure.

G = get_roms_grid(GRDname, INIname);

spherical = G.spherical;

[Lr,Mr] = size(G.h);
Nr = length(G.s_rho);

% Set switches for state variables to nudge.

LnudgeM2CLM    = true;           % nudging 2D momentum
LnudgeM3CLM    = true;           % nudging 3D momentum
LnudgeTCLM     = true;           % nudging tracers (usually T-S)
LnudgeTgeneric = true;           % nudging generic tracers

% Set NetCDF variables to process.  Initialize inverse nudging
% coefficients with zeros.  Recall that division by zero is not
% defined and we will give "Inf".  Therefore, we just need to set
% only the values in the desired areas and the rest can be zero
% (for no nudging) because the nudging in ROMS is:
%
%      F(...,new) = F(...,new) +
%                   dt * F_nudgcoef * (Fclm - F(...,new))

VarNUD = [];
F = [];

if (spherical),
  VarNUD = {'lon_rho', 'lat_rho'};
else
  VarNUD = {'x_rho', 'y_rho'};
end

if (LnudgeM2CLM),
  VarNUD = [VarNUD, 'M2_NudgeCoef'];
  F.M2_NudgeCoef = zeros(Lr,Mr);                % RHO-points
end

if (LnudgeM3CLM),
  VarNUD = [VarNUD, 'M3_NudgeCoef'];
  F.M3_NudgeCoef = zeros(Lr,Mr,Nr);             % RHO-points
end

if (LnudgeTCLM),
  VarNUD = [VarNUD, 'temp_NudgeCoef', 'salt_NudgeCoef'];
  F.temp_NudgeCoef = zeros(Lr,Mr,Nr);
  F.salt_NudgeCoef = zeros(Lr,Mr,Nr);
end

if (LnudgeTgeneric),
  VarNUD = [VarNUD, 'tracer_NudgeCoef'];
  F.tracer_NudgeCoef = zeros(Lr,Mr,Nr);
end

%--------------------------------------------------------------------------
% Create Nudging coefficients NetCDF file: build creation parameters
% structure, S.
%--------------------------------------------------------------------------

S            = [];
Tindex       = [];
ReplaceValue = NaN;
PreserveType = true;
Unlimited    = false;                   % time dimension is umlimited
nctype       = 'nc_double';             % input data is in double precision

mode = netcdf.getConstant('CLOBBER');                    % overwrite!!!
mode = bitor(mode,netcdf.getConstant('64BIT_OFFSET'));

% The strategy here is to build manually the NetCDF metadata structure to
% facilitate creating several NetCDF file in a generic and compact way.
% This structure is similar to that returned by "nc_inq" or native Matlab
% function "ncinfo".
%
% Notice that we call the "roms_metadata" function to create the fields
% in the structure.  Then, we call "check_metadata" for fill unassigned
% values and to check for consistency.

ncname = NUDname;

disp(blanks(1));
disp(['** Creating NetCDF file: ', ncname,' **']);
disp(blanks(1))

S.Filename = ncname;

S.Attributes(1).Name      = 'type';
S.Attributes(1).Value     = 'Nudging Coeffcients file';

S.Attributes(2).Name      = 'title';
S.Attributes(2).Value     = ['North Atlantic Damee #4, 0.75 Resolution'];

S.Attributes(3).Name      = 'grd_file';
S.Attributes(3).Value     = GRDname;

S.Attributes(4).Name      = 'ini_file';
S.Attributes(4).Value     = INIname;

S.Attributes(5).Name      = 'history';
S.Attributes(5).Value     = ['Nudging coefficient file created from ',  ...
                             'd_nudgcoef.m: ', date_stamp];

S.Dimensions(1).Name      = 'xi_rho';
S.Dimensions(1).Length    = Lr;
S.Dimensions(1).Unlimited = false;

S.Dimensions(2).Name      = 'eta_rho';
S.Dimensions(2).Length    = Mr;
S.Dimensions(2).Unlimited = false;

S.Dimensions(3).Name        = 's_rho';
S.Dimensions(3).Length      = Nr;
S.Dimensions(3).Unlimited   = false;

S.Variables(1) = roms_metadata('spherical');
if (spherical),
  S.Variables(2) = roms_metadata('lon_rho');
  S.Variables(3) = roms_metadata('lat_rho');
else
  S.Variables(2) = roms_metadata('x_rho');
  S.Variables(3) = roms_metadata('y_rho');
end

% Process inverse nudging coefficient variables.

for n = 3:length(VarNUD),

  Vname  = char(VarNUD{n});
  
  i = n + 1;
  S.Variables(i) = roms_metadata(Vname, spherical, nctype, Unlimited);

end
  
% Check ROMS metadata structure.  Fill unassigned fields.

S = check_metadata(S);
  
% Create forcing NetCDF files.  Write grid coordinates.

ncid = nc_create(S.Filename, mode, S);    % create a new NetCDF file

status = nc_write(S.Filename, 'spherical', int32(spherical));
status = nc_write(S.Filename, 'lon_rho',   G.lon_rho);
status = nc_write(S.Filename, 'lat_rho',   G.lat_rho);

%--------------------------------------------------------------------------
% Set inverse time scales.
%--------------------------------------------------------------------------

IstrR = 0;
IendR = Lr-1;
JstrR = 0;
JendR = Mr-1;

% In the North Atlantic DAMEE_4 application the nudging is done in the
% southern and northern domain edges over a 8-point linearly tapered
% nudging scales of 5 to 60 days.

inner = 1/60;                        % 60 days at interior limit
outer = 1/5;                         %  5 days at boundary
width = 8;                           %  8 points

work  = zeros(Lr,Mr);

for j=JstrR:width,                   % Southern boundary
  for i=IstrR:IendR,
    work(i+1,j+1) = inner + (width - j) * (outer - inner) / width;
  end
end

for j=JendR-width:JendR,             % Northern boundary
  for i=IstrR:IendR,
    work(i+1,j+1) = outer + (JendR - j) * (inner - outer) / width;
  end
end

fac = (7 * outer - inner)/6;

for j=74:80,                         % Mediterranean outflow
  for i=102:106,
    cff = sqrt((i-109)^2 + (j-77)^2);
    work(i+1,j+1) = max(0, (fac + cff * (inner - outer) / 6));
  end
end

% Load values into structure

if (LnudgeM2CLM),
  F.M2_NudgeCoef = work;
end

if (LnudgeM3CLM),
  F.M3_NudgeCoef = repmat(work, [1 1 Nr]);
end

if (LnudgeTCLM),
  F.temp_NudgeCoef = repmat(work, [1 1 Nr]);
  F.salt_NudgeCoef = repmat(work, [1 1 Nr]);
end

if (LnudgeTgeneric),
  F.tracer_NudgeCoef = repmat(work, [1 1 Nr]);
end

%--------------------------------------------------------------------------
% Write out inverse nudging coefficients.
%--------------------------------------------------------------------------

Nfields = length(VarNUD);

for i=3:Nfields,
  field = char(VarNUD(i));
  status = nc_write(NUDname, field, F.(field));
end
