function [x,y,f] = nc_slice(Ginp, Hfile, Vname, depth, varargin)

%
% NC_SLICE:  Interpolate requested slice from a 3D NetCDf variable
%
% [x,y,f] = nc_slice(Gfile, Hfile, Vname, depth, Tindex, zflag, method)
%
% This function computes a horizontal slice from a NetCDF file generated
% by ROMS. This function should only be used when the grid has variable
% bathymetry.  The field slice is interpolated at the requested depth
% from the input terrain-following coordinate data.
%
% On Input:
%
%    Ginp        Grid NetCDF file name (string): if no grid file, use
%                  field file name instead
%            or, an existing grid structure (struct array)
%
%    Hfile       Field history NetCDF file name (string)
%
%    Vname       NetCDF variable name to process (string)
%
%    depth       Slice depth (scalar; meters, negative)
%
%    Tindex      Time index to read (Optional, scalar):
%
%                  Tindex = Inf    Process last record (default)
%                  Tindex = 0      All time records are processed
%                  Tindex > 0      only "Tindex" record is processed
%
%    zflag       Variable depths computation switch (Optional, scalar):
%
%                  zflag  = 0      Use zero free-surface (default)
%                  zflag  = 1      Read in free-surface at "Tindex"
%
%    method      Interpolation method (string):
%
%                  'nearest'       Nearest neighbor
%                  'linear'        Linear (default)
%                  'spline'        Piecewise cubic spline
%                  'cubic'         Cubic
%
% On Output:
%
%    x           Slice X-positions (2D array)
%
%    y           Slice Y-positions (2D array)
%
%    f           Field slice (array)
%

% svn $Id: nc_slice.m 711 2014-01-23 20:36:13Z arango $
%=========================================================================%
%  Copyright (c) 2002-2014 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%

% Initialize.

f = [];
x = [];
y = [];
Z = [];

isr3d = false;
isw3d = false;

Nrecs = 0;
Tname = [];
Tsize = 0;

recordless = true;

% Set optional arguments.

Tindex = Inf;
zflag  = 0;
method = 'linear';

switch numel(varargin)
  case 1
    Tindex = varargin{1};
  case 2
    Tindex = varargin{1};
    zflag  = varargin{2};
  case 3
    Tindex = varargin{1};
    zflag  = varargin{2};
    method = varargin{3};
end

% If applicable, get ROMS grid structure.

if (~isstruct(Ginp)),
  G = get_roms_grid(Ginp);
else
  G = Ginp;
end

%--------------------------------------------------------------------------
% Determine horizontal positions and Land/Sea masking variables.
%--------------------------------------------------------------------------

I = nc_vinfo(Hfile,Vname);

%  Check variable dimensions and determine horizontal/vertical
%  coordinates and Land/Sea mask arrays.

nvdims = length(I.Dimensions);

if (nvdims > 0),
  for n=1:nvdims,
    Dname = char(I.Dimensions(n).Name);
    switch Dname
      case 's_rho'
        isr3d = true;
      case 's_w'
        isw3d = true;
      case 'xi_rho'
        Mask = G.mask_rho;
        if (G.spherical),
          x = G.lon_rho;
          y = G.lat_rho;
        else
          x = G.x_rho;
          y = G.y_rho;
        end
        if (isfield(G, 'z_r')),
	  Z = G.z_r;
        end
      case 'xi_psi'
        Mask = G.mask_psi;
        if (G.spherical),
          x = G.lon_psi;
          y = G.lat_psi;
        else
          x = G.x_psi;
          y = G.y_psi;
        end
        if (isfield(G, 'z_psi')),
	  Z = G.z_psi;
        end
      case 'xi_u'
        Mask = G.mask_u;
        if (G.spherical),
          x = G.lon_u;
          y = G.lat_u;
        else
          x = G.x_u;
          y = G.y_u;
        end 
        if (isfield(G, 'z_u')),
	  Z = G.z_u;
        end
      case 'xi_v'
        Mname = G.mask_v;
        if (G.spherical),
          x = G.lon_v;
          y = G.lat_v;
        else
          x = G.x_v;
          y = G.y_v;
        end
        if (isfield(G, 'z_v')),
	  Z = G.z_v;
        end
      case 'ocean_time'
        recordless = false;    
        Nrecs = I.Size(end);
        Tsize = I.Dimensions(n).Length;
    end
  end
  if (isw3d),
    Zname = 'z_w';
  end  
end

is3d = isr3d || isw3d;

if (~is3d),
  error(['nc_slice - cannot vertically interpolate: ',Vname]);
end

%  If appropriate, process last time record available.

if (~recordless && Tindex > Tsize),
  Tindex = Tsize;
end 

% Determine starting and ending record to process.

if (Tindex == 0),
  StrRec = 1;
  EndRec = Nrecs;
else
  StrRec = Tindex;
  EndRec = Tindex;
end

%--------------------------------------------------------------------------
% Compute grid depths.
%--------------------------------------------------------------------------

igrid = I.Cgridtype.Value;

if (zflag == 0),
  ZetaRec = 0;
else
  ZetaRec = StrRec;
end

if (isempty(Z) || zflag > 0),
  if (~isstruct(Ginp)),
    Z = depths(Hfile,Ginp,igrid,0,ZetaRec);
  else
    Z = depths(Hfile,Hfile,igrid,0,ZetaRec);
  end
end

[Im Jm Km] = size(Z);

%--------------------------------------------------------------------------
% Interpolate slice to requested homogeneous depth.
%--------------------------------------------------------------------------

f  = NaN([Im Jm EndRec-StrRec+1]);
ic = 0;

for Rec=StrRec:EndRec,
  ic = ic + 1;
  V  = nc_read(Hfile,Vname,Rec);

  for j=1:Jm,
    for i=1:Im,
      Zwrk = reshape(Z(i,j,:),1,Km);
      Vwrk = reshape(V(i,j,:),1,Km);
      f(i,j,ic) = interp1(Zwrk,Vwrk,depth,method);
    end
  end
end

%  Remove singlenton dimension if processing a single record.

f = squeeze(f);

return
