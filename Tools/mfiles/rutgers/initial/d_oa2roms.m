%
% D_OA2ROMS:  Driver script to create ROMS climatology NetCDF file
%               from OA NetCDF files.
%
% This is a user modifiable script showing how to create ROMS climatology
% NetCDF files. The data source are objective analyzed (OA) Annual and
% Monthly temperature and salinity fields from the Levitus (1998) dataset.
% These fields are generated using ROMS OA package.
%
% This is a template script showing how to create several NetCDF files
% (one per each atmospheric forcing field) using ROMS metadata structure,
% which follows the schema of "nc_inq" or native 'ncinfo" function.
%

% svn $Id: d_oa2roms.m 996 2020-01-10 04:28:56Z arango $
%=========================================================================%
%  Copyright (c) 2002-2020 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%
%
% Input ROMS Grid NetCDF file.

Dir = fullfile(getenv('HOME'), 'ocean/repository/Projects/gom/Data');

GrdFile =  fullfile(Dir, 'gom_grid_a.nc');

% Set input OA file names.

AnnFile =  fullfile(Dir, 'OA', 'HDAT', 'oa_lev98_ann.nc');

MonFile = {fullfile(Dir, 'OA', 'HDAT', 'oa_lev98_jan.nc'),              ...
           fullfile(Dir, 'OA', 'HDAT', 'oa_lev98_feb.nc'),              ...
           fullfile(Dir, 'OA', 'HDAT', 'oa_lev98_mar.nc'),              ...
           fullfile(Dir, 'OA', 'HDAT', 'oa_lev98_apr.nc'),              ...
           fullfile(Dir, 'OA', 'HDAT', 'oa_lev98_may.nc'),              ...
           fullfile(Dir, 'OA', 'HDAT', 'oa_lev98_jun.nc'),              ...
           fullfile(Dir, 'OA', 'HDAT', 'oa_lev98_jul.nc'),              ...
           fullfile(Dir, 'OA', 'HDAT', 'oa_lev98_aug.nc'),              ...
           fullfile(Dir, 'OA', 'HDAT', 'oa_lev98_sep.nc'),              ...
           fullfile(Dir, 'OA', 'HDAT', 'oa_lev98_oct.nc'),              ...
           fullfile(Dir, 'OA', 'HDAT', 'oa_lev98_nov.nc'),              ...
           fullfile(Dir, 'OA', 'HDAT', 'oa_lev98_dec.nc')};

% Set output concatenated files.  The monthly Levitus climatology is
% available from the surface to 1500m.  The other bottom levels are
% extracted from the annual climatology. The premise here is that the
% ocean belon 1500m does not change on monthly/seasonal time scales.

OutFile = {fullfile(Dir, 'OA', 'gom_lev98_jan.nc'),                     ...
           fullfile(Dir, 'OA', 'gom_lev98_feb.nc'),                     ...
           fullfile(Dir, 'OA', 'gom_lev98_mar.nc'),                     ...
           fullfile(Dir, 'OA', 'gom_lev98_apr.nc'),                     ...
           fullfile(Dir, 'OA', 'gom_lev98_may.nc'),                     ...
           fullfile(Dir, 'OA', 'gom_lev98_jun.nc'),                     ...
           fullfile(Dir, 'OA', 'gom_lev98_jul.nc'),                     ...
           fullfile(Dir, 'OA', 'gom_lev98_aug.nc'),                     ...
           fullfile(Dir, 'OA', 'gom_lev98_sep.nc'),                     ...
           fullfile(Dir, 'OA', 'gom_lev98_oct.nc'),                     ...
           fullfile(Dir, 'OA', 'gom_lev98_nov.nc'),                     ...
           fullfile(Dir, 'OA', 'gom_lev98_dec.nc')};

% Set output ROMS climatology files.

AnnCLM  =  fullfile(Dir, 'gom_lev98_ann.nc');
MonCLM  =  fullfile(Dir, 'gom_lev98_clm.nc');

%  Set intepolation parameters and switches.

method       = 'linear';       % vetical interpolation method
Concatenate  = true;           % concatenate monthly and annual OA files
nctype       = 'nc_float';     % Input data is in single precision
Unlimited    = true;           % time dimension is umlimited in
                               % output files
Tindex       = [];
ReplaceValue = NaN;
PreserveType = true;

mode = netcdf.getConstant('CLOBBER');                    % overwrite!!!
mode = bitor(mode,netcdf.getConstant('64BIT_OFFSET'));

%--------------------------------------------------------------------------
% Concatenate annual and monthly OA NetCDF files.
%--------------------------------------------------------------------------

Nfiles = length(MonFile);

if (Concatenate),
  for n=1:Nfiles,
    oa_cat(char(OutFile(n)), AnnFile, char(MonFile(n)));
  end
end

%--------------------------------------------------------------------------
%  Set application parameters in structure array, S.
%--------------------------------------------------------------------------

[Lr,Mr] = size(nc_read(GrdFile, 'h'));

Lu = Lr-1;
Lv = Lr;
Mu = Mr;
Mv = Mr-1;

S.ncname      = MonCLM;    % output NetCDF file

S.spherical   = 1;          % spherical grid

S.Lm          = Lr-2;       % number of interior RHO-points, X-direction
S.Mm          = Mr-2;       % number of interior RHO-points, Y-direction
S.N           = 40;         % number of vertical levels at RHO-points
S.NT          = 2;          % total number of tracers

S.Vtransform  = 2;          % vertical transfomation equation
S.Vstretching = 4;          % vertical stretching function

S.theta_s     = 7.0;        % S-coordinate surface control parameter
S.theta_b     = 2.0;        % S-coordinate bottom control parameter
S.Tcline      = 250.0;      % S-coordinate surface/bottom stretching width
S.hc          = S.Tcline;   % S-coordinate stretching width

%--------------------------------------------------------------------------
% Set variables to process.
%--------------------------------------------------------------------------

VarGrd  = {'spherical',                                                 ...
           'Vtransform', 'Vstretching',                                 ...
           'theta_s', 'theta_b', 'Tcline', 'hc',                        ...
           's_rho', 'Cs_r', 'h'};

if (S.spherical),
  VarGrd = [VarGrd, 'lon_rho', 'lat_rho', 'mask_rho'];
else
  VarGrd = [VarGrd, 'x_rho', 'y_rho', 'mask_rho'];
end

VarCLM = {'temp', 'salt', 'rho'};

ROMSvars = [VarGrd, 'clm_time', VarCLM];

%--------------------------------------------------------------------------
% Set ROMS Grid structure. The depths are for an unperturbed state
% (zeta = 0).
%--------------------------------------------------------------------------

G = get_roms_grid(GrdFile, S);

%--------------------------------------------------------------------------
% Create annual and montly climatology NetCDF files: build creation
% parameters structure, C.
%--------------------------------------------------------------------------
%
% The strategy here is to build manually the NetCDF metadata structure to
% facilitate creating several NetCDF file in a generic and compact way.
% This structure is similar to that returned by "nc_inq" or native Matlab
% function "ncinfo".
%
% Notice that we call the "roms_metadata" function to create the fields
% in the structure.  Then, we call "check_metadata" for fill unassigned
% values and to check for consistency.

C.Filename = AnnCLM;

C.Attributes(1).Name      = 'type';
C.Attributes(1).Value     = 'CLIMATOLOGY file';

C.Attributes(2).Name      = 'title';
C.Attributes(2).Value     = ['Annual Levitus 1998 Climatology ',        ...
                             'Dataset, Gulf of Mexico Region'];

C.Attributes(3).Name      = 'history';
C.Attributes(3).Value     = sprintf('%s:  Created by %s with %s.m',     ...
                                    datestr(now), getenv('USER'),       ...
                                    mfilename);

C.Dimensions(1).Name      = 'xi_rho';
C.Dimensions(1).Length    = Lr;
C.Dimensions(1).Unlimited = false;

C.Dimensions(2).Name      = 'eta_rho';
C.Dimensions(2).Length    = Mr;
C.Dimensions(2).Unlimited = false;

C.Dimensions(3).Name      = 's_rho';
C.Dimensions(3).Length    = S.N;
C.Dimensions(3).Unlimited = false;

C.Dimensions(4).Name      = 'clm_time';
C.Dimensions(4).Length    = nc_constant('nc_unlimited');
C.Dimensions(4).Unlimited = true;

% Set ROMS variables metadata.

ic = 0;

for var = ROMSvars,
  vname = char(var);
  ic = ic + 1;
  C.Variables(ic) = roms_metadata(vname, S.spherical, nctype, Unlimited);
end

% Edit the time variable 'units" attribute for the correct reference
% time and add calendar attribute.

ivar   = strcmp({C.Variables.Name}, 'clm_time');
nvatts = length(C.Variables(ivar).Attributes);  
iatt   = strcmp({C.Variables(ivar).Attributes.Name}, 'units');

C.Variables(ivar).Attributes(iatt).Value = 'days since 2000-01-01 00:00:00';

C.Variables(ivar).Attributes(nvatts+1).Name  = 'calendar';
C.Variables(ivar).Attributes(nvatts+1).Value = '365_day';

C.Variables(ivar).Attributes(nvatts+2).Name  = 'cycle_length';
C.Variables(ivar).Attributes(nvatts+2).Value = double(365);

% Check ROMS metadata structure.  Fill unassigned fields.

C = check_metadata(C);
  
% Create annual climatology NetCDF file and write grid arrays.

ncid = nc_create(AnnCLM, mode, C);

for var = VarGrd,
  vname = char(var);
  status = nc_write(AnnCLM, vname, G.(vname));
  if (status ~= 0),
    error(['D_OA2ROMS: error while writing variable: ', vname]);
  end
end

% Create annual climatology NetCDF file.

C.Filename = MonCLM;

C.Attributes(2).Value = ['Monthly Levitus 1998 Climatology ',           ...
                         'Dataset, Gulf of Mexico Region'];

ncid = nc_create(MonCLM, mode, C);

for var = VarGrd,
  vname = char(var);
  status = nc_write(MonCLM, vname, G.(vname));
  if (status ~= 0),
    error(['D_OA2ROMS: error while writing variable: ', vname]);
  end
end

%--------------------------------------------------------------------------
% Vertically interpolate OA fields to ROMS terrain-following coordinates.
% Notice that ROMS in situ  density is computed from potential temperature
% and salinity.  This will be used for nudging density in the future.
%--------------------------------------------------------------------------

Tindex = 1;                    % The OA files have a single time record

% Annual climatology fields.

A.temp = oa2roms(G, AnnFile, 'temp', Tindex, method);
A.salt = oa2roms(G, AnnFile, 'salt', Tindex, method);
A.rho  = roms_eos(A.temp, A.salt, G.z_r);

for var = VarCLM,
  field = char(var);
  status = nc_write(AnnCLM, field, A.(field), Tindex);
  if (status ~= 0),
    error(['D_OA2ROMS: error while writing annual field: ', field]);
  end
end

% Monthly climatology fields.

for n=1:Nfiles,
  M.temp = oa2roms(G, char(OutFile(n)), 'temp', Tindex, method);
  M.salt = oa2roms(G, char(OutFile(n)), 'salt', Tindex, method);
  M.rho  = roms_eos(M.temp, M.salt, G.z_r);

  for var = VarCLM,
    field = char(var);
    status = nc_write(MonCLM, field, M.(field), n);
    if (status ~= 0),
      error(['D_OA2ROMS: error while monthly field: ', field]);
    end
  end
end
