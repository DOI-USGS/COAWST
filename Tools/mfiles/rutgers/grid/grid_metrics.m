function [pm, pn, dndx, dmde] = grid_metrics(G, GreatCircle)

%
% GRID_METRICS:  Compute ROMS Grid horizontal metrics
%
% [pm, pn, dndx, dmde] = grid_metrics(G, GreatCircle)
%
% This function computes horizontal grid spacing metrics from
% Grid NetCDF file or Grid structure G.
%
% On Input:
%
%    G            Either Grid NetCDF file name (character string) or
%                   Grid structure (structure array)
%
%    GreatCircle  Switch indicating how to compute the grid distance:
%                   GreatCircle = true     Great-circle distance
%                   GreatCircle = false    Cartesian distance
%
% On Output:
%
%    pm         Curvilinear coordinate metric in the XI-direction
%                 (1/meters; dx = 1/pm)
%
%    pm         Curvilinear coordinate metric in the ETA-direction
%                 (1/meters; dy = 1/pn)
%
%    dndx       XI-derivative  of inverse metric factor pn (meters),
%                 d(pn)/d(XI)
%
%    dmde       ETA-derivative of inverse metric factor pm (meters),
%                 d(pm)/d(ETA)
%
% If G is a Grid structure and GreatCircle=true, the following values are
% needed to compute the great circle distances (G.spherical must be 1):
%
%    G.spherical     Grid spherical flag [0 | 1]
%    G.lon_rho       RHO-points longitude (decimal degrees)
%    G.lat_rho       RHO-points latitude  (decimal degrees)
%    G.lon_u         U-points   longitude (decimal degrees)
%    G.lat_u         U-points   latitude  (decimal degrees)
%    G.lon_v         V-points   longitude (decimal degrees)
%    G.lat_v         V-points   latitude  (decimal degrees)
%
%                    longitude:  positive East,   negative West
%                    latitude:   positive North,  negative South
%
% Otherwise, if G is a grid structure and GreatCircle=false, the following
% values are needed to compute Cartesian distances regardless the value of
% G.spherical:
%
%    G.spherical     Grid spherical flag [0 | 1]
%    G.x_rho         RHO-points X-location (meters)
%    G.y_rho         RHO-points Y-location (meters)
%    G.x_u           U-points   X-location (meters)
%    G.x_u           U-points   Y-location (meters)
%    G.x_v           V-points   X-location (meters)
%    G.x_v           V-points   Y-location (meters)

  
% svn $Id: grid_metrics.m 711 2014-01-23 20:36:13Z arango $
%=========================================================================%
%  Copyright (c) 2002-2014 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%

% Get Grid coordinates.

if (ischar(G)),
  G = get_roms_grid(G);
end

spherical = G.spherical;

if (GreatCircle && spherical),
  Xr = G.lon_rho;
  Yr = G.lat_rho;
  Xu = G.lon_u;
  Yu = G.lat_u;
  Xv = G.lon_v;
  Yv = G.lat_v;
else
  Xr = G.x_rho;
  Yr = G.y_rho;
  Xu = G.x_u;
  Yu = G.y_u;
  Xv = G.x_v;
  Yv = G.y_v;
end
  
%----------------------------------------------------------------------------  
%  Compute grid spacing (meters).
%----------------------------------------------------------------------------  

[Lp,Mp] = size(Xr);

L = Lp-1;   Lm = L-1;
M = Mp-1;   Mm = M-1;

dx = zeros(size(Xr));
dy = zeros(size(Xr));

% Compute grid spacing.

if (GreatCircle && spherical),            % great circle distances
  BEARING = false;
  
  dx(2:L,1:Mp) = gcircle(Xu(1:Lm,1:Mp), Yu(1:Lm,1:Mp),                ...
                         Xu(2:L ,1:Mp), Yu(2:L ,1:Mp), BEARING);
  dx(1  ,1:Mp) = gcircle(Xr(1   ,1:Mp), Yr(1   ,1:Mp),                ...
                         Xu(1   ,1:Mp), Yu(1   ,1:Mp), BEARING).*2.0;
  dx(Lp ,1:Mp) = gcircle(Xr(Lp  ,1:Mp), Yr(Lp  ,1:Mp),                ...
                         Xu(L   ,1:Mp), Yu(L   ,1:Mp), BEARING).*2.0;  

  dy(1:Lp,2:M) = gcircle(Xv(1:Lp,1:Mm), Yv(1:Lp,1:Mm),                ...
                         Xv(1:Lp,2:M ), Yv(1:Lp,2:M ), BEARING);
  dy(1:Lp,1  ) = gcircle(Xr(1:Lp,1   ), Yr(1:Lp,1   ),                ...
                         Xv(1:Lp,1   ), Yv(1:Lp,1   ), BEARING).*2.0;
  dy(1:Lp,Mp ) = gcircle(Xr(1:Lp,Mp  ), Yr(1:Lp,Mp  ),                ...
                         Xv(1:Lp,M   ), Yv(1:Lp,M   ), BEARING).*2.0;

  dx = dx .* 1000;       % great circle function computes
  dy = dy .* 1000;       % distances in kilometers
  
else                                      % Cartesian distances

  dx(2:L,1:Mp) = sqrt((Xu(2:L ,1:Mp) - Xu(1:Lm,1:Mp)).^2 +            ...
                      (Yu(2:L ,1:Mp) - Yu(1:Lm,1:Mp)).^2);
  dx(1  ,1:Mp) = sqrt((Xu(1   ,1:Mp) - Xr(1   ,1:Mp)).^2 +            ...
                      (Yu(1   ,1:Mp) - Yr(1   ,1:Mp)).^2).*2.0;
  dx(Lp ,1:Mp) = sqrt((Xr(Lp  ,1:Mp) - Xu(L   ,1:Mp)).^2 +            ...
                      (Yr(Lp  ,1:Mp) - Yu(L   ,1:Mp)).^2).*2.0;

  dy(1:Lp,2:M) = sqrt((Xv(1:Lp,2:M ) - Xv(1:Lp,1:Mm)).^2 +            ...
                      (Yv(1:Lp,2:M ) - Yv(1:Lp,1:Mm)).^2);
  dy(1:Lp,1  ) = sqrt((Xv(1:Lp,1   ) - Xr(1:Lp,1   )).^2 +            ...
                      (Yv(1:Lp,1   ) - Yr(1:Lp,1   )).^2).*2.0;
  dy(1:Lp,Mp ) = sqrt((Xr(1:Lp,Mp  ) - Xv(1:Lp,M   )).^2 +            ...
                      (Yr(1:Lp,Mp  ) - Yv(1:Lp,M   )).^2).*2.0;
end

% Compute inverse grid spacing metrics.
  
pm = 1.0./dx;
pn = 1.0./dy;

%----------------------------------------------------------------------------  
% Compute inverse metric derivatives.
%----------------------------------------------------------------------------  

dndx = zeros(size(Xr));
dmde = zeros(size(Xr));

if (~G.uniform),
  dndx(2:L,2:M) = 0.5.*(1.0./pn(3:Lp,2:M ) - 1.0./pn(1:Lm,2:M ));
  dmde(2:L,2:M) = 0.5.*(1.0./pm(2:L ,3:Mp) - 1.0./pm(2:L ,1:Mm));
% jcw mirror boundary data
  dndx(1,:)=dndx(2,:); dndx(end,:)=dndx(end-1,:);
  dndx(:,1)=dndx(:,2); dndx(:,end)=dndx(:,end-1);
  dndx(1,1)=    0.5*(dndx(1,2)      +dndx(2,1));
  dndx(1,end)=  0.5*(dndx(1,end-1)  +dndx(2,end));
  dndx(end,1)=  0.5*(dndx(end,2)    +dndx(end-1,1));
  dndx(end,end)=0.5*(dndx(end-1,end)+dndx(end,end-1));

  dmde(1,:)=dmde(2,:); dmde(end,:)=dmde(end-1,:);
  dmde(:,1)=dmde(:,2); dmde(:,end)=dmde(:,end-1);
  dmde(1,1)=    0.5*(dmde(1,2)      +dmde(2,1));
  dmde(1,end)=  0.5*(dmde(1,end-1)  +dmde(2,end));
  dmde(end,1)=  0.5*(dmde(end,2)    +dmde(end-1,1));
  dmde(end,end)=0.5*(dmde(end-1,end)+dmde(end,end-1));
end

return
