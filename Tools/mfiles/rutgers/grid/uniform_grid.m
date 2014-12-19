function G = uniform_grid (dx,dy,L,M)

%
% UNIFORM_GRID:  Sets a ROMS grid uniform coordinates and metrics
%
% Given the grid spacing (dx, dy) and number of PSI-points (L, M) in
% the X- and Y-direction, this function sets the ROMS grid Cartesian
% coordinates and metrics.
%
% Notice that ROMS has a staggered Arakawa's C-grid:
%
%    L          M                Number of PSI-points
%    Lr = L+1   Mr = M+1         Number of RHO-points
%    Lu = L     Mu = M+1         Number of U-points
%    Lv = L+1   Mv = M           Number of V-points
%
% On Input:
%
%    dx           Grid spacing in the X-direction (meters)
%    dy           Grid spacing in the Y-direction (meters) 
%    L            Number of PSI grid-points in the X-direction
%    M            Number of PSI grid-points in the Y-direction
%
% On Output:
%
%    G            Uniform grid coordinates and metrics (struct array)
%
%                   G.x_psi    X-coordinates at PSI-points (meters)
%                   G.y_psi    Y-coordinates at PSI-points (meters)
%                   G.x_rho    X-coordinates at RHO-points (meters)
%                   G.y_rho    Y-coordinates at RHO-points (meters)
%                   G.x_u      X-coordinates at U-points (meters)
%                   G.y_u      Y-coordinates at U-points (meters)
%                   G.x_v      X-coordinates at V-points (meters)
%                   G.y_v      Y-coordinates at V-points (meters)
%                   G.pm       X-coordinates metric "m" (1/dx, 1/meters)
%                   G.pn       Y-coordinates metric "n" (1/dy, 1/meters)
%                   G.dndx     Inverse metric, d(1/pn)/d(x)
%                   G.dmde     Inverse metric, d(1/pm)/d(y) 
%

% svn $Id: uniform_grid.m 711 2014-01-23 20:36:13Z arango $
%=========================================================================%
%  Copyright (c) 2002-2014 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%

% Initialize.
  
Lr = L+1;   Mr = M+1;
Lu = L;     Mu = M+1;
Lv = L+1;   Mv = M;

% Compute uniform Cartesian coordinates in meters.

halfdx = 0.5*dx;
halfdy = 0.5*dy;

[G.y_psi, G.x_psi] = meshgrid(      0:dy:(M -1)*dy,                     ...
                                    0:dx:(L -1)*dx);
[G.y_rho, G.x_rho] = meshgrid(-halfdy:dy:(Mr-1)*dy,                     ...
                              -halfdx:dx:(Lr-1)*dx);
[G.y_u  , G.x_u  ] = meshgrid(-halfdy:dy:(Mu-1)*dy,                     ...
                                    0:dx:(Lu-1)*dx); 
[G.y_v  , G.x_v  ] = meshgrid(      0:dy:(Mv-1)*dy,                     ...
                              -halfdx:dx:(Lv-1)*dx); 

% Compute inverse grid spacing metrics.

G.pm = ones(size(G.x_rho))./dx;
G.pn = ones(size(G.y_rho))./dy;

% Compute inverse metric derivatives. It should be zero for uniform
% grid spacing:
%
% G.dndx(2:L,2:M) = 0.5.*(1.0./G.pn(3:Lr,2:M ) - 1.0./G.pn(1:L-1,2:M  ));
% G.dmde(2:L,2:M) = 0.5.*(1.0./G.pm(2:L ,3:Mr) - 1.0./G.pm(2:L  ,1:M-1));

G.dndx = zeros(size(G.x_rho));
G.dmde = zeros(size(G.y_rho));

return
