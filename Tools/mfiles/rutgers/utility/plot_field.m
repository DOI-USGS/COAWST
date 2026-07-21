function F=plot_field(Gname, Hname, Vname, Tindex, varargin)

%
% PLOT_FIELD:  Plot requested ROMS variable from input NetCDF file
%
% F=plot_field(Gname, Hname, Vname, Tindex, Level, Caxis, Mmap, ptype,
%              wrtPNG)
%
% This function plots requested ROMS variable from input history
% NetCDF file. This function is very useful when debugging a ROMS
% application. The plotting is not that fancy but it provides enough
% information for browsing ROMS variables very quickly. The location
% of the minumum is marked with a filled magenta circle whereas the
% maximum is marked with a filled magenta square.
%
% It also draws the coastline if the data is available in the
% grid file or structure.
%
% On Input:
%
%    Gname         ROMS Grid NetCDF file/URL name (string)
%              or, an existing ROMS grid structure (struct array)
%
%    Hname         ROMS history NetCDF file/URL name (string)
%
%    Vname         ROMS NetCDF variable name to plot (string)
%
%    Tindex        Time record index used to process (scalar)
%                    (Use Inf or inf for last record)
%
%    Level         If 3D variable, vertical level to plot (scalar)
%                    (Optional, default: surface level)
%
%                     Level > 0,    terrain-following level
%                     Level < 0,    depth (field interpolation)
%
%    Caxis         Color axis range (optional; vector)
%                    (default: [-Inf Inf])
%
%    Mmap          Switch to use m_map utility (optional; integer)
%                    Mmap = 0,  no map projection (default)
%                    Mmap = 1,  'm_map' utility
%                    Mmap = 2,  native Matlab toolbox
%
%    ptype         Plot type (optional; integer)
%                    ptype < 0     use 'contourf' with abs(ptype) colors
%                    ptype = 1     use 'pcolor'
%                    ptype = 2     use 'pcolorjw'
%
%    wrtPNG        Switch to write out PNG file (optional; switch)
%                    (default: 0)
%
%                  if wrtPNG < 1, ommit figure title, doTitle = false
%
% On Output:
%
%    F             Requested 2D or 3D variable (structure)
%

% svn $Id$
%=======================================================================%
%  Copyright (c) 2002-2025 The ROMS Group                               %
%    Licensed under a MIT/X style license                               %
%    See License_ROMS.md                            Hernan G. Arango    %
%=======================================================================%

% Initialize.

F = struct('ncname'     , [], 'Vname'     , [],                       ...
           'Tindex'     , [], 'Tname'     , [], 'Tstring'   , [],     ...
           'Level'      , [], 'is3d'      , [],                       ...
           'X'          , [], 'Y'         , [],                       ...
           'value'      , [], 'min'       , [], 'max'       , [],     ...
           'Caxis'      , [], 'doMap'     , [], 'projection', [],     ...
           'ptype'      , [],                                         ...
           'gotCoast'   , [], 'lon_coast' , [], 'lat_coast' , [],     ...
           'shading'    , [], 'pltHandle' , [], 'wrtPNG'    , []);

F.projection = 'mercator';

F.ncname = Hname;
F.Vname = Vname;
F.gotCoast = false;

got.Mname = false;
got.Xname = false;
got.Yname = false;
got.Zname = false;

is3d  = false;
isr3d = false;
isw3d = false;
iszflat = false;

F.Tname = [];
Tsize = 0;
recordless = true;
Tstring = blanks(2);

D = nc_dinfo(Hname);
N = D(strcmp({D.Name}, 's_rho')).Length;

%  Optional arguments.

switch numel(varargin)
  case 0
    Level  = N;
    Caxis  = [-Inf Inf];
    Mmap   = false;
    ptype  = 0;
    wrtPNG = false;
  case 1
    if (~isinf(varargin{1}))
      Level = varargin{1};
    else
      Level = N;
    end
    Caxis  = [-Inf Inf];
    Mmap   = false;
    ptype  = 0;
    wrtPNG = false;
  case 2
    if (~isinf(varargin{1}))
      Level = varargin{1};
    else
      Level = N;
    end
    Caxis  = varargin{2};
    Mmap   = false;
    ptype  = 0;
    wrtPNG = false;
  case 3
    if (~isinf(varargin{1}))
      Level = varargin{1};
    else
      Level = N;
    end
    Caxis  = varargin{2};
    Mmap   = varargin{3};
    ptype  = 0;
    wrtPNG = false;
  case 4
    if (~isinf(varargin{1}))
      Level = varargin{1};
    else
      Level = N;
    end
    Caxis  = varargin{2};
    Mmap   = varargin{3};
    ptype  = varargin{4};
    wrtPNG = false;
  case 5
    if (~isinf(varargin{1}))
      Level = varargin{1};
    else
      Level = N;
    end
    Caxis  = varargin{2};
    Mmap   = varargin{3};
    ptype  = varargin{4};
    wrtPNG = varargin{5};
end

F.Level   = Level;
F.Caxis   = Caxis;
F.doMap   = Mmap;
F.ptype   = ptype;
F.shading = 'interp';
F.wrtPNG  = wrtPNG;

% Set ROMS Grid structure.

if (~isstruct(Gname))
  G = get_roms_grid(Gname);
else
  G = Gname;
end

%--------------------------------------------------------------------------
% Get Variable information.
%--------------------------------------------------------------------------

I = nc_vinfo(Hname, Vname);

%  Check variable dimensions and determine horizontal/vertical
%  coordinates and Land/Sea mask arrays.

nvdims = length(I.Dimensions);

if (nvdims > 0)
  for n=1:nvdims
    dimnam = char(I.Dimensions(n).Name);
    switch dimnam
      case {'s_rho'}
        isr3d = true;
      case 's_w'
        isw3d = true;
      case 'z_slice'
        iszflat = true;
        z_slice = nc_read(Hname, 'z_slice');
        F.Level = z_slice(Level);
      case {'xi_rho','eta_rho'}
        Mname = 'mask_rho';
        got.Mname = true;
        if (~(got.Xname || got.Yname))
          if (G.spherical)
            Xname = 'lon_rho';
            Yname = 'lat_rho';
          else
            Xname = 'x_rho';
            Yname = 'y_rho';
          end
          Zname = 'z_r';
          got.Xname = true;
          got.Yname = true;
          got.Zname = true;
        end
      case {'xi_psi','eta_psi'}
        Mname = 'mask_psi';
        got.Mname = true;
        if (~(got.Xname || got.Yname))
          if (G.spherical)
            Xname = 'lon_psi';
            Yname = 'lat_psi';
          else
            Xname = 'x_psi';
            Yname = 'y_psi';
          end
          Zname = 'z_p';
          got.Xname = true;
          got.Yname = true;
          got.Zname = true;
        end
      case {'xi_u','eta_u'}
        Mname = 'mask_u';
        got.Mname = true;
        if (~(got.Xname || got.Yname))
          if (G.spherical)
            Xname = 'lon_u';
            Yname = 'lat_u';
          else
            Xname = 'x_u';
            Yname = 'y_u';
          end
          Zname = 'z_u';
          got.Xname = true;
          got.Yname = true;
          got.Zname = true;
        end
      case {'xi_v','eta_v'}
        Mname = 'mask_v';
        got.Mname = true;
        if (~(got.Xname || got.Yname))
          if (G.spherical)
            Xname = 'lon_v';
            Yname = 'lat_v';
          else
            Xname = 'x_v';
            Yname = 'y_v';
          end
          Zname = 'z_v';
          got.Xname = true;
          got.Yname = true;
          got.Zname = true;
        end
      case {'ocean_time', 'time', ~contains(dimnam,'time')}
        recordless = false;
        Tsize = I.Dimensions(n).Length;
    end
  end
  if (isw3d)
    Zname = 'z_w';
  end
end

is3d = isr3d || isw3d || iszflat;

F.is3d = is3d;

% Check if 'time' attribute is available.

itime = strcmp({I.Attributes.Name}, 'time');
if (any(itime))
  F.Tname = I.Attributes(itime).Value;
end

%------------------------------------------------------------------------
% Get coordinates.
%------------------------------------------------------------------------

if (isfield(G,Xname))
  if (~isempty(G.(Xname)))
    X = G.(Xname);
  else
    error([' PLOT_FIELD - field '', Xname, ''',                       ...
           ' is empty in receiver grid structure: G']);
  end
else
  error([' PLOT_FIELD - unable to find field '', Xname, ''',          ...
         ' in receiver grid structure: G']);
end

if (isfield(G,Yname))
  if (~isempty(G.(Yname)))
    Y = G.(Yname);
  else
    error([' PLOT_FIELD - field '', Yname, ''',                       ...
           ' is empty in receiver grid structure: G']);
  end
else
  error([' PLOT_FIELD - unable to find field '', Yname, ''',          ...
         ' in receiver grid structure: G']);
end

if (is3d)
  if (isfield(G,Zname))
    if (~isempty(G.(Zname)))
      if (iszflat)
        Z = z_slice;
      else
        Z = G.(Zname);
      end
    else
      error([' PLOT_FIELD - field '', Zname, ''',                     ...
             ' is empty in receiver grid structure: G']);
    end
  else
    error([' PLOT_FIELD - unable to find field '', Zname, ''',        ...
           ' in receiver grid structure: G']);
  end
end

if (isfield(G,Mname))
  if (~isempty(G.(Mname)))
    mask = G.(Mname);
  else
    error([' PLOT_FIELD - field '', Mname, ''',                       ...
          ' is empty in receiver grid structure: G']);
  end
else
  error([' PLOT_FIELD - unable to find field '', Mname, ''',          ...
         ' in receiver grid structure: G']);
end

if (isfield(G,'lon_coast') && isfield(G,'lat_coast'))
  F.lon_coast = G.lon_coast;
  F.lat_coast = G.lat_coast;
  F.gotCoast = true;
end

if (~G.spherical)
  X = 0.001 .* X;                 % km
  Y = 0.001 .* Y;                 % km
end

F.X = X;
F.Y = Y;

%------------------------------------------------------------------------
% Read in requested variable from NetCDF file.
%------------------------------------------------------------------------

if (~recordless && Tindex > Tsize)
  Tindex = Tsize;                    % process last time record available
end
F.Tindex = Tindex;

if (~isempty(F.Tname))
  Tvalue = nc_read(Hname,F.Tname,Tindex);
  Tattr  = nc_getatt(Hname,'units',F.Tname);
  Tdays  = true;
  if (contains(Tattr, 'second'))
    Tvalue = Tvalue/86400;                    % seconds to days
    Tdays  = false;
  end
  iatt = strfind(Tattr, 'since');
  if (~isempty(iatt))
    Torigin = Tattr(iatt+6:end);
    epoch   = datenum(Torigin,31);            % 'yyyy-mm-dd HH:MM:SS'
    Tstring = datestr(epoch+Tvalue);
  else
    Tstring = num2str(Tvalue);
  end
end
F.Tstring = Tstring;

ReplaceValue = NaN;
PreserveType = false;

if (recordless)
  field = nc_read(Hname,Vname,[],ReplaceValue,PreserveType);
  F.shading = 'interp';
else
  field = nc_read(Hname,Vname,Tindex,ReplaceValue,PreserveType);
end

if (is3d || iszflat)
  if (Level > 0)
    value = squeeze(field(:,:,Level));
  else
    [~,~,value] = nc_slice(G,Hname,Vname,Level,Tindex);
  end
else
  value = field;
end

F.value = value;

%------------------------------------------------------------------------
% Plot requested field
%------------------------------------------------------------------------

F = hplot(G, F);

return
