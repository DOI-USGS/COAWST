function [status] = c_grid(Lp,Mp,Gname,varargin)

%
% C_GRID:  Creates or modifies a ROMS grid NetCDF file.
%
% [status] = c_grid(Lp,Mp,Gname,NewFile,spherical)
%
% This function creates or modifies a existing ROMS grid NetCDF file.
%
% On Input:
%
%    Lp          Number of RHO-points in the XI-direction
%    Mp          Number of RHO-points in the ETA-direction
%    Gname       Grid NetCDF file name (string)
%    NewFile     Switch to create a new file (logical, OPTIONAL)
%    spherical   Spherical grid switch (default true)
%
% On Output:
%
%    status      Error flag.

% svn $Id: c_grid.m 711 2014-01-23 20:36:13Z arango $
%=========================================================================%
%  Copyright (c) 2002-2014 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%

% Set optional arguments.

NewFile = false;
spherical = true;

switch numel(varargin),
  case 1
    NewFile = varargin{1};
  case 2
    NewFile   = varargin{1};
    spherical = varargin{2};
end

GOT_NCFILE = false;
Cartesian = true;

% Check if grid NetCDF already exist.

disp(' ');
if (~NewFile),
  if (exist(Gname,'file')),
    GOT_NCFILE = true;
    disp(['Appending to existing GRID NetCDF file: ',Gname]);
  else
    disp(['Creating a new GRID NetCDF file: ',Gname]);
  end
  disp(' ');
end

%--------------------------------------------------------------------------
%  Inquire dimensions from a existing NeTCDF file.
%--------------------------------------------------------------------------

Dname = {'xi_rho',  'xi_u',  'xi_v',  'xi_psi',                       ...
         'eta_rho', 'eta_u', 'eta_v', 'eta_psi',                      ...
         'bath'};

for value = Dname,
  field = char(value);
  gotDim.(field) = false;
  switch (field)
    case {'xi_rho', 'xi_v'}
      Dsize.(field) = Lp;
    case {'xi_psi', 'xi_u'}
      Dsize.(field) = Lp-1;
    case {'eta_rho', 'eta_u'}
      Dsize.(field) = Mp;
    case {'eta_psi', 'eta_v'}
      Dsize.(field) = Mp-1;
    case {'bath'}
      Dsize.(field) = netcdf.getConstant('UNLIMITED');
  end
end

if (GOT_NCFILE),

  D = nc_dinfo(Gname);
  ndims = length(D);
  for n=1:ndims,
    name = char(D(n).Name);
    if (sum(strcmp(Dname,name)) == 1),
      Dsize.(name) = D(n).Length;
      gotDim.(name) = true;
      did.(name) = D(n).dimid;
    end
  end

  if (gotDim.xi_rho),
    if (~gotDim.xi_u),   Dsize.xi_u   = Dsize.xi_rho - 1; end
    if (~gotDim.xi_v),   Dsize.xi_v   = Dsize.xi_rho;     end
    if (~gotDim.xi_eta), Dsize.xi_psi = Dsize.xi_rho - 1; end
  end

  if (gotDim.eta_rho),
    if (~gotDim.eta_u),   Dsize.eta_u   = Dsize.eta_rho;     end
    if (~gotDim.eta_v),   Dsize.eta_v   = Dsize.eta_rho - 1; end
    if (~gotDim.eta_psi), Dsize.eta_psi = Dsize.eta_rho - 1; end
  end
  if (~gotDim.bath), Dsize.bath = nc_constant('nc_unlimited'); end

end

%--------------------------------------------------------------------------
%  Inquire Variables from a existing NeTCDF file.
%--------------------------------------------------------------------------

grd_vars = {'spherical', 'xl', 'el',                                  ...
            'angle', 'pm', 'pn', 'dndx', 'dmde',                      ...
            'f', 'h', 'hraw'};

if (Cartesian),
  grd_vars = [grd_vars, 'x_rho', 'y_rho', 'x_psi', 'y_psi',           ...
                        'x_u',   'y_u',   'x_v',   'y_v'];
end

if (spherical),
  grd_vars = [grd_vars, 'lon_rho', 'lat_rho', 'lon_psi', 'lat_psi',   ...
                        'lon_u',   'lat_u',   'lon_v',   'lat_v'];
end
  
grd_vars = [grd_vars, 'mask_rho', 'mask_psi', 'mask_u', 'mask_v'];

for value = grd_vars,
  field = char(value);
  Vname.(field) = field;
  gotVar.(field) = false;
end
  
if (GOT_NCFILE),
  V = nc_vnames(Gname);
  nvars = length(V.Variables);
  for n = 1:nvars,
    name = char(V.Variables(n).Name);
    if (sum(strcmp(grd_vars,name)) == 1),
      gotVar.(name) = true;
    end
  end
end

%--------------------------------------------------------------------------
%  If applicable, create initial/climatology NetCDF file.
%--------------------------------------------------------------------------

if (~GOT_NCFILE),
  mode = netcdf.getConstant('CLOBBER');
  mode = bitor(mode, netcdf.getConstant('64BIT_OFFSET'));

  ncid = netcdf.create(Gname, mode);
end

%--------------------------------------------------------------------------
%  If applicable, open GRID NetCDF file and put in definition mode.
%--------------------------------------------------------------------------

if (GOT_NCFILE),
  ncid  = netcdf.open(Gname, 'nc_write');
  netcdf.reDef(ncid);
end

%--------------------------------------------------------------------------
%  Create global attribute(s).
%--------------------------------------------------------------------------

if (~GOT_NCFILE),
  varid  = netcdf.getConstant('nc_global');

  type = 'GRID file';
  netcdf.putAtt(ncid, varid, 'type', type);
  
  history = ['GRID file using Matlab script: c_grid, ', date_stamp];
  netcdf.putAtt(ncid, varid, 'history', history);
end

%--------------------------------------------------------------------------
%  If appropriate, define dimensions.
%--------------------------------------------------------------------------

for value = Dname,
  field = char(value);
  if (~gotDim.(field)),
    did.(field) = netcdf.defDim(ncid, field, Dsize.(field));
  end
end

%--------------------------------------------------------------------------
%  If appropriate, define variables.
%--------------------------------------------------------------------------

% Define spherical switch.

if (~gotVar.spherical),
  Var.name          = Vname.spherical;
  Var.type          = netcdf.getConstant('nc_int');
  Var.dimid         = [];
  Var.long_name     = 'grid type logical switch';
  Var.flag_values   = [0 1];
  Var.flag_meanings = ['Cartesian', blanks(1),                        ...
                       'spherical'];
  [~,status] = nc_vdef(ncid,Var);
  if (status ~= 0), return, end
  clear Var
end

% Define basin length.

if (~gotVar.xl),
  Var.name          = Vname.xl;
  Var.type          = netcdf.getConstant('nc_double');
  Var.dimid         = [];
  Var.long_name     = 'basin length in the XI-direction';
  Var.units         = 'meter';
  [~,status] = nc_vdef(ncid,Var);
  if (status ~= 0), return, end
  clear Var
end

if (~gotVar.el),
  Var.name          = Vname.el;
  Var.type          = netcdf.getConstant('nc_double');
  Var.dimid         = [];
  Var.long_name     = 'basin length in the ETA-direction';
  Var.units         = 'meter';
  [~,status] = nc_vdef(ncid,Var);
  if (status ~= 0), return, end
  clear Var
end

%  Curvilinear rotation angle on RHO-points.

if (~gotVar.angle),
  Var.name          = Vname.angle;
  Var.type          = netcdf.getConstant('nc_double');
  Var.dimid         = [did.eta_rho did.xi_rho];
  Var.long_name     = 'angle between XI-axis and EAST';
  Var.units         = 'radians';
  [~,status] = nc_vdef(ncid,Var);
  if (status ~= 0), return, end
  clear Var
end

%  Curvilinear coordinates metrics at RHO-points.

if (~gotVar.pm),
  Var.name          = Vname.pm;
  Var.type          = netcdf.getConstant('nc_double');
  Var.dimid         = [did.eta_rho did.xi_rho];
  Var.long_name     = 'curvilinear coordinate metric in XI';
  Var.units         = 'meter-1';
  [~,status] = nc_vdef(ncid,Var);
  if (status ~= 0), return, end
  clear Var
end

if (~gotVar.pn),
  Var.name          = Vname.pn;
  Var.type          = netcdf.getConstant('nc_double');
  Var.dimid         = [did.eta_rho did.xi_rho];
  Var.long_name     = 'curvilinear coordinate metric in ETA';
  Var.units         = 'meter-1';
  [~,status] = nc_vdef(ncid,Var);
  if (status ~= 0), return, end
  clear Var
end

%  Curvilinear coordinates inverse metric derivative.

if (~gotVar.dndx),
  Var.name          = Vname.dndx;
  Var.type          = netcdf.getConstant('nc_double');
  Var.dimid         = [did.eta_rho did.xi_rho];
  Var.long_name     = 'XI-derivative of inverse metric factor pn';
  Var.units         = 'meter';
  [~,status] = nc_vdef(ncid,Var);
  if (status ~= 0), return, end
  clear Var
end

if (~gotVar.dmde),
  Var.name          = Vname.dmde;
  Var.type          = netcdf.getConstant('nc_double');
  Var.dimid         = [did.eta_rho did.xi_rho];
  Var.long_name     = 'ETA-derivative of inverse metric factor pm';
  Var.units         = 'meter';
  [~,status] = nc_vdef(ncid,Var);
  if (status ~= 0), return, end
  clear Var
end

%  Coriolis Parameter at RHO-points.

if (~gotVar.f),
  Var.name          = Vname.f;
  Var.type          = netcdf.getConstant('nc_double');
  Var.dimid         = [did.eta_rho did.xi_rho];
  Var.long_name     = 'Coriolis parameter at RHO-points';
  Var.units         = 'second-1';
  [~,status] = nc_vdef(ncid,Var);
  if (status ~= 0), return, end
  clear Var
end

%  Raw bathymetry at RHO-points.

if (~gotVar.hraw),
  Var.name          = Vname.hraw;
  Var.type          = netcdf.getConstant('nc_double');
  Var.dimid         = [did.bath did.eta_rho did.xi_rho];
  Var.long_name     = 'Working bathymetry at RHO-points';
  Var.units         = 'meter';
  [~,status] = nc_vdef(ncid,Var);
  if (status ~= 0), return, end
  clear Var
end

%  Model bathymetry at RHO-points.

if (~gotVar.h),
  Var.name          = Vname.h;
  Var.type          = netcdf.getConstant('nc_double');
  Var.dimid         = [did.eta_rho did.xi_rho];
  Var.long_name     = 'model bathymetry at RHO-points';
  Var.units         = 'meter';
  [~,status] = nc_vdef(ncid,Var);
  if (status ~= 0), return, end
  clear Var
end

%  Cartesian locations of RHO-points.

if (Cartesian),
  if (~gotVar.x_rho),
    Var.name          = Vname.x_rho;
    Var.type          = netcdf.getConstant('nc_double');
    Var.dimid         = [did.eta_rho did.xi_rho];
    Var.long_name     = 'X-location of RHO-points';
    Var.units         = 'meter';
    [~,status] = nc_vdef(ncid,Var);
    if (status ~= 0), return, end
    clear Var
  end

  if (~gotVar.y_rho),
    Var.name          = Vname.y_rho;
    Var.type          = netcdf.getConstant('nc_double');
    Var.dimid         = [did.eta_rho did.xi_rho];
    Var.long_name     = 'Y-location of RHO-points';
    Var.units         = 'meter';
    [~,status] = nc_vdef(ncid,Var);
    if (status ~= 0), return, end
    clear Var
  end  
end

%  Cartesian locations of PSI-points.

if (Cartesian),
  if (~gotVar.x_psi),
    Var.name          = Vname.x_psi;
    Var.type          = netcdf.getConstant('nc_double');
    Var.dimid         = [did.eta_psi did.xi_psi];
    Var.long_name     = 'X-location of PSI-points';
    Var.units         = 'meter';
    [~,status] = nc_vdef(ncid,Var);
    if (status ~= 0), return, end
    clear Var
  end

  if (~gotVar.y_psi),
    Var.name          = Vname.y_psi;
    Var.type          = netcdf.getConstant('nc_double');
    Var.dimid         = [did.eta_psi did.xi_psi];
    Var.long_name     = 'Y-location of PSI-points';
    Var.units         = 'meter';
    [~,status]=nc_vdef(ncid,Var);
    if (status ~= 0), return, end
    clear Var
  end
end

%  Cartesian locations of U-points.

if (Cartesian),
  if (~gotVar.x_u),
    Var.name          = Vname.x_u;
    Var.type          = netcdf.getConstant('nc_double');
    Var.dimid         = [did.eta_u did.xi_u];
    Var.long_name     = 'X-location of U-points';
    Var.units         = 'meter';
    [~,status]=nc_vdef(ncid,Var);
    if (status ~= 0), return, end
    clear Var
  end

  if (~gotVar.y_u),
    Var.name          = Vname.y_u;
    Var.type          = netcdf.getConstant('nc_double');
    Var.dimid         = [did.eta_u did.xi_u];
    Var.long_name     = 'Y-location of U-points';
    Var.units         = 'meter';
    [~,status]=nc_vdef(ncid,Var);
    if (status ~= 0), return, end
    clear Var
  end
end

%  Cartesian locations of V-points.

if (Cartesian),
  if (~gotVar.x_v),
    Var.name          = Vname.x_v;
    Var.type          = netcdf.getConstant('nc_double');
    Var.dimid         = [did.eta_v did.xi_v];
    Var.long_name     = 'X-location of V-points';
    Var.units         = 'meter';
    [~,status] = nc_vdef(ncid,Var);
    if (status ~= 0), return, end
    clear Var
  end

  if (~gotVar.y_v),
    Var.name          = Vname.y_v;
    Var.type          = netcdf.getConstant('nc_double');
    Var.dimid         = [did.eta_v did.xi_v];
    Var.long_name     = 'Y-location of V-points';
    Var.units         = 'meter';
    [~,status] = nc_vdef(ncid,Var);
    if (status ~= 0), return, end
    clear Var
  end
end

%  Longitude/latitude of RHO-points.

if (spherical),
  if (~gotVar.lon_rho),
    Var.name          = Vname.lon_rho;
    Var.type          = netcdf.getConstant('nc_double');
    Var.dimid         = [did.eta_rho did.xi_rho];
    Var.long_name     = 'longitude of RHO-points';
    Var.units         = 'degree_east';
    Var.standard_name = 'longitude';
    [~,status] = nc_vdef(ncid,Var);
    if (status ~= 0), return, end
    clear Var
  end

  if (~gotVar.lat_rho),
    Var.name          = Vname.lat_rho;
    Var.type          = netcdf.getConstant('nc_double');
    Var.dimid         = [did.eta_rho did.xi_rho];
    Var.long_name     = 'latitute of RHO-points';
    Var.units         = 'degree_north';
    Var.standard_name = 'latitude';
    [~,status] = nc_vdef(ncid,Var);
    if (status ~= 0), return, end
    clear Var
  end
end

%  Longitude/latitude of PSI-points.

if (spherical)
  if (~gotVar.lon_psi),
    Var.name          = Vname.lon_psi;
    Var.type          = netcdf.getConstant('nc_double');
    Var.dimid         = [did.eta_psi did.xi_psi];
    Var.long_name     = 'longitude of PSI-points';
    Var.units         = 'degree_east';
    Var.standard_name = 'longitude';
    [~,status] = nc_vdef(ncid,Var);
    if (status ~= 0), return, end
    clear Var
  end

  if (~gotVar.lat_psi),
    Var.name          = Vname.lat_psi;
    Var.type          = netcdf.getConstant('nc_double');
    Var.dimid         = [did.eta_psi did.xi_psi];
    Var.long_name     = 'latitute of PSI-points';
    Var.units         = 'degree_north';
    Var.standard_name = 'latitude';
    [~,status] = nc_vdef(ncid,Var);
    if (status ~= 0), return, end
    clear Var
  end
end

%  Longitude/latitude of U-points.

if (spherical),
  if (~gotVar.lon_u),
    Var.name          = Vname.lon_u;
    Var.type          = netcdf.getConstant('nc_double');
    Var.dimid         = [did.eta_u did.xi_u];
    Var.long_name     = 'longitude of U-points';
    Var.units         = 'degree_east';
    Var.standard_name = 'longitude';
    [~,status] = nc_vdef(ncid,Var);
    if (status ~= 0), return, end
    clear Var
  end
  
  if (~gotVar.lat_u),
    Var.name          = Vname.lat_u;
    Var.type          = netcdf.getConstant('nc_double');
    Var.dimid         = [did.eta_u did.xi_u];
    Var.long_name     = 'latitute of U-points';
    Var.units         = 'degree_north';
    Var.standard_name = 'latitude';
    [~,status] = nc_vdef(ncid,Var);
    if (status ~= 0), return, end
    clear Var
  end
end

%  Longitude/latitude of V-points.

if (spherical),
  if (~gotVar.lon_v),
    Var.name          = Vname.lon_v;
    Var.type          = netcdf.getConstant('nc_double');
    Var.dimid         = [did.eta_v did.xi_v];
    Var.long_name     = 'longitude of V-points';
    Var.units         = 'degree_east';
    Var.standard_name = 'longitude';
    [~,status] = nc_vdef(ncid,Var);
    if (status ~= 0), return, end
    clear Var
  end

  if (~gotVar.lat_v),
    Var.name          = Vname.lat_v;
    Var.type          = netcdf.getConstant('nc_double');
    Var.dimid         = [did.eta_v did.xi_v];
    Var.long_name     = 'latitute of V-points';
    Var.units         = 'degree_north';
    Var.standard_name = 'latitude';
    [~,status] = nc_vdef(ncid,Var);
    if (status ~= 0), return, end
    clear Var
  end
end

%  Land/sea mask on RHO-points.

if (~gotVar.mask_rho),
  Var.name          = Vname.mask_rho;
  Var.type          = netcdf.getConstant('nc_double');
  Var.dimid         = [did.eta_rho did.xi_rho];
  Var.long_name     = 'mask on RHO-points';
  Var.flag_values   = [0.0 1.0];
  Var.flag_meanings = ['land', blanks(1),                             ...
                       'water'];
  [~,status] = nc_vdef(ncid,Var);
  if (status ~= 0), return, end
  clear Var
end

%  Land/sea mask on PSI-points.

if (~gotVar.mask_psi),
  Var.name          = Vname.mask_psi;
  Var.type          = netcdf.getConstant('nc_double');
  Var.dimid         = [did.eta_psi did.xi_psi];
  Var.long_name     = 'mask on PSI-points';
  Var.flag_values   = [0.0 1.0];
  Var.flag_meanings = ['land', blanks(1), ...
                       'water'];
  [~,status] = nc_vdef(ncid,Var);
  if (status ~= 0), return, end
  clear Var
end

%  Land/sea mask on U-points.

if (~gotVar.mask_u),
  Var.name          = Vname.mask_u;
  Var.type          = netcdf.getConstant('nc_double');
  Var.dimid         = [did.eta_u did.xi_u];
  Var.long_name     = 'mask on U-points';
  Var.flag_values   = [0.0 1.0];
  Var.flag_meanings = ['land', blanks(1),                             ...
                       'water'];
  [~,status] = nc_vdef(ncid,Var);
  if (status ~= 0), return, end
  clear Var
end

%  Land/sea mask on V-points.

if (~gotVar.mask_v),
  Var.name          = Vname.mask_v;
  Var.type          = netcdf.getConstant('nc_double');
  Var.dimid         = [did.eta_v did.xi_v];
  Var.long_name     = 'mask on V-points';
  Var.flag_values   = [0.0 1.0];
  Var.flag_meanings = ['land', blanks(1), ...
                       'water'];
  [~,status] = nc_vdef(ncid,Var);
  if (status ~= 0), return, end
  clear Var
end

%--------------------------------------------------------------------------
%  Leave definition mode and close NetCDF file.
%--------------------------------------------------------------------------

netcdf.endDef(ncid);
netcdf.close(ncid);

return
