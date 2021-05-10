function R = set_grid(rlon, rlat, ncname, varargin)

%
% SET_GRID:  Creates a ROMS grid NetCDF file.
%
% set_grid(rlon, rlat, Gname, method, database)
%
% Given the longitude and latitude (rlon, rlat) at RHO-points (cell center)
% this function creates a ROMS grid NetCDF file. 
%
% On Input:
%
%    rlon        Longitude of RHO-points (2D array; degree_east)
%
%    rlat        Latitude  of RHO-points (2D array, degree_north)
%
%    ncname      Grid NetCDF file name (string)
%
%    method      Interpolation method for griddedInterpolant (OPTIONAL)
%                  'linear'    (default)
%                  'spline'
%                  'cubic'
%
%    database    GSHHS coastline database for land/sea masking (OPTIONAL)
%                  'f'    full resolution
%                  'h'    high resolution
%                  'i'    intermediate resolution (default)
%                  'l'    load resolution
%                  'c'    crude resolution  
%
% On Output:
%
%    R           ROMS Grid structure
%

% svn $Id$
%=========================================================================%
%  Copyright (c) 2002-2020 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%

% Initialize.

switch numel(varargin)
  case 0
    method = 'linear';
    extract_coast = false;
  case 1
    method = varargin{1};
    extract_coast = false;
  case 2
    method = varargin{1};
    database = varargin{2};
    extract_coast = true;
end
    
G.name = ncname;
G.spherical = 1;
G.lon_rho = rlon;
G.lat_rho = rlat;

% Set fractional coordinates.

[Lp, Mp] = size(G.lon_rho);

L = Lp-1;
M = Mp-1;

[Yp, Xp] = meshgrid(1.0:1:M     , 1.0:1:L     );
[Yr, Xr] = meshgrid(0.5:1:Mp-0.5, 0.5:1:Lp-0.5);
[Yu, Xu] = meshgrid(0.5:1:Mp-0.5, 1.0:1:L     );
[Yv, Xv] = meshgrid(1.0:1:M     , 0.5:1:Lp-0.5);

% Interpolate coordinates at PSI-, U-, V-points.

Fr = griddedInterpolant(Xr, Yr, G.lon_rho, method);

                         G.lon_psi = Fr(Xp, Yp);
                         G.lon_u   = Fr(Xu, Yu);
                         G.lon_v   = Fr(Xv, Yv);

Fr.Values = G.lat_rho;   G.lat_psi = Fr(Xp, Yp);
                         G.lat_u   = Fr(Xu, Yu);
                         G.lat_v   = Fr(Xv, Yv);

% Compute ROMS grid metrics and load its values into G structure.

S = roms_metrics(G);

G.x_rho = S.x_rho;
G.y_rho = S.y_rho;
G.x_psi = S.x_psi;
G.y_psi = S.y_psi;
G.x_u   = S.x_u;
G.y_u   = S.y_u;
G.x_v   = S.x_v;
G.y_v   = S.y_v;

G.xl    = S.xl;
G.el    = S.el;

G.angle = S.angle;

G.pm    = S.pm;
G.pn    = S.pn;
G.dndx  = S.dndx;
G.dmde  = S.dmde;

G.f     = S.f

G.h     = zeros(size(G.lon_rho));

% Initialize land/sea masking arrays.

G.mask_psi = ones(size(G.lon_psi));
G.mask_rho = ones(size(G.lon_rho));
G.mask_u   = ones(size(G.lon_u));
G.mask_v   = ones(size(G.lon_v));

%-------------------------------------------------------------------------
% Create ROMS grid NetCDF file and write available variables.
%-------------------------------------------------------------------------

c_grid(Lp, Mp, G.name, true);

% Write out available variables.

grd_vars = {'spherical', 'xl', 'el',                                    ...
            'f', 'h', 'angle', 'pm', 'pn', 'dndx', 'dmde',              ...
            'x_rho',   'y_rho',   'x_psi',   'y_psi',                   ...
            'x_u',     'y_u',     'x_v',     'y_v',                     ...
            'lon_rho', 'lat_rho', 'lon_psi', 'lat_psi',                 ...
            'lon_u',   'lat_u',   'lon_v',   'lat_v',                   ...
	    'mask_rho', 'mask_psi', 'mask_u', 'mask_v'};

for var = grd_vars
  field = char(var);
  nc_write(G.name, field, G.(field));
end

%-------------------------------------------------------------------------
% Extract coastlines and compute/write land/sea masking arrays.
%-------------------------------------------------------------------------

if (extract_coast)

% Compute and write land/sea mask.
  
  F = landsea(G.name, database);

  G.mask_psi = F.mask_psi;
  G.mask_rho = F.mask_rho;
  G.mask_u   = F.mask_u;
  G.mask_v   = F.mask_v;
  
  G.lon_coast = F.clon;
  G.lat_coast = F.clat;

% Add coastline to grid NetCDF file.

  add_coastline(G.name, G.lon_coast, G.lat_coast);

end

%-------------------------------------------------------------------------
% Recompute ROMS full grid structure.
%-------------------------------------------------------------------------

R = get_roms_grid(G.name);

return
