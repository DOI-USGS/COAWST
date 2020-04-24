function [lon,lat,u,v] = vector4stream(U,V,G,dx,dy,varargin)

%
% VECTOR4STREAM:  Interpolates vector to monotonic and plaid coordinates
%                 for ploting streamlines
%
% [lon,lat,u,v] = vector4stream(U,V,G,dx,dy,Llon,Rlon,Blat,Tlat)
%
% Given velocity components and grid information, it interpolates data to
% a monotonic and plaid grid (as if produced by MESHGRID) for the plotting
% of streamlines elsewhere using 'streamslice' or 'm_streamline'.  The user
% may specify the monotonic grid corners bounded by (Llon,Blat) and
% (Rlon,Tlat).
%
%                 ______ (Rlon, Tlat)
%                |      |
%                |      |
%                |______|
%    (Llon, Blat)
%
% or it can be determined internally to include the application domain.
% The interpolation will yield NaN values for unbouded grid points.
%
% If needed, it rotates the velocity components to geographical EAST and
% NORTH components.
%
% On Input:
%
%    U          Input U-velocity component (2D or 3D array) at staggered
%                 U-points or RHO-points.
%
%    V          Input V-velocity component (2D or 3D array) at staggered
%                 V-points or RHO-points.
%
%    G          Application grid NetCDF filename (string) 
%           or, an existing ROMS grid structure (struct array)
%
%    dx         Monotonic grid spacing in the longitude axis (positive)
%
%    dy         Monotonic grid spacing in the latitude  axis (positive)
%
%    Llon       Monotonic grid left-edge   longitude (degrees, OPTIONAL)
%
%    Rlon       Monotonic grid right-edge  longitude (degrees, OPTIONAL)
%
%    Blat       Monotonic grid bottom-edge latitude  (degress, OPTIONAL)
%
%    Tlat       Monotonic grid top-edge    latitude  (degress, OPTIONAL)
%
% On Output:
%
%    lon        Monotonic grid longitude (2D array)
%
%    lat        Monotonic grid latitude  (2D array)
%
%    u          Interpolated U-velocity data into monotonic grid
%
%    v          Interpolated V-velocity data into monotonic grid
%

% svn $Id: vector4stream.m 996 2020-01-10 04:28:56Z arango $
%=========================================================================%
%  Copyright (c) 2002-2020 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%

% Initialize.

Smethod = 'linear';                            % interpolation method

% Get ROMS application grid structure.

if (~isstruct(G)),
  ncfile = G;
  G = get_roms_grid(ncfile);
end

% Set monotonic and plaid grid.

switch numel(varargin)
  case 0
    Llon = min(G.lon_rho(:));
    Rlon = max(G.lon_rho(:));
    Blat = min(G.lat_rho(:));
    Tlat = max(G.lat_rho(:));
  case 1
    Llon = varargin{1};
    Rlon = max(G.lon_rho(:));
    Blat = min(G.lat_rho(:));
    Tlat = max(G.lat_rho(:));
  case 2
    Llon = varargin{1};
    Rlon = varargin{2};
    Blat = min(G.lat_rho(:));
    Tlat = max(G.lat_rho(:));
  case 3
    Llon = varargin{1};
    Rlon = varargin{2};
    Blat = varargin{3};
    Tlat = max(G.lat_rho(:));
  case 4
    Llon = varargin{1};
    Rlon = varargin{2};
    Blat = varargin{3};
    Tlat = varargin{4};
end

[lon,lat] = meshgrid(Llon:dx:Rlon, Blat:dy:Tlat);

%  Determine if 2D or 3D vector field and dimension size.

[Lp,Mp]=size(G.lon_rho);

L=Lp-1;  Lm=L-1;  Lm2=Lm-1;
M=Mp-1;  Mm=M-1;  Mm2=Mm-1;

if (length(size(U)) == 2),
  is3d=false;
  [Lu,Mu]=size(U);
  [Lv,Mv]=size(V);
elseif (length(size(U)) == 3),
  is3d=true;
  [Lu,Mu,N]=size(U);
  [Lv,Mv,N]=size(V);
end

%  Determine if input vector components are at RHO-points.


if ((Lu == Lp) && (Mv == Mp)),
  rho_points = true;
  ulon = G.lon_rho;
  ulat = G.lat_rho;
  vlon = G.lon_rho;
  vlat = G.lat_rho;
else
  rho_points = false;
  ulon = G.lon_u;
  ulat = G.lat_u;
  vlon = G.lon_v;
  vlat = G.lat_v;
end

if (is3d)
  lon = repmat(lon,[1,1,N]);
  lat = repmat(lat,[1,1,N]);

  ulon = repmat(ulon,[1,1,N]);
  ulat = repmat(ulat,[1,1,N]);

  vlon = repmat(vlon,[1,1,N]);
  vlat = repmat(vlat,[1,1,N]);
end

% Determine the points inside of the application domain perimenter.

[IN,ON]=inpolygon(lon(:), lat(:), G.lon_perimeter, G.lat_perimeter);

IN(ON) = true;  

%-------------------------------------------------------------------------
% If applicable, rotate input velocity components to geographical EAST and
% NORTH components.
%-------------------------------------------------------------------------

irotate = 1;                          % rotate from (XI, ETA) to (lon,lat)

if (G.vector_rotation)
  [U,V] = rotate_vec(U, V, G.angle, irotate);    
end

%-------------------------------------------------------------------------
% Interpolate input velocity components to a monotonic and plaid grid.
%-------------------------------------------------------------------------

Fu = scatteredInterpolant(ulon(:), ulat(:), U(:), Smethod);
u  = Fu(lon,lat);

Fv = scatteredInterpolant(vlon(:), vlat(:), V(:), Smethod);
v  = Fv(lon,lat);

% Flag out points outside of the application domain perimeter.

u(~IN) = NaN;
v(~IN) = NaN;

return
