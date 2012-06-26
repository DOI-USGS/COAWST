function [h]=pltscope(Gname, ptype, C);

%
% PLTSCOPE:  Plots adjoint sensitivity scope masks
%
% [h]=pltscope(Gname, ptype, C)
%
% This function plots the adjoint sensitivity scope mask. If
% appropriate, It also ovelays the given coastline data.
%
% On Input:
%
%    Gname       NetCDF file name (character string).
%    ptype       plot type flag:
%                  ptype=0   => Cartesian axis.
%                  ptype=1   => Spherical axis.
%    C           Coastline indices (structure array):
%                  C.grid    => Grid file name.
%                  C.coast   => Coastline file name.
%                  C.indices => Coastline indices file name.
%                  C.lon     => Coastline longitudes.
%                  C.lat     => Coastline latitudes.
%                  C.Icst    => Coastline I-grid coordinates, (0:L).
%                  C.Jcst    => Coastline J-grid coordinates, (0:M).
%

% svn $Id: pltscope.m 436 2010-01-02 17:07:55Z arango $
%===========================================================================%
%  Copyright (c) 2002-2010 The ROMS/TOMS Group                              %
%    Licensed under a MIT/X style license                                   %
%    See License_ROMS.txt                           Hernan G. Arango        %
%===========================================================================%

%  Read spherical switch.

spherical=0;
string=nc_read(Gname,'spherical');
if (upper(string) == 'T'),
  spherical=1;
end,

%  Set coastline switch

ICOAST=0;
if (spherical & (nargin >2)),
  ICOAST=1;
end,

%---------------------------------------------------------------------------
% Plot scope mask on RHO-points.
%---------------------------------------------------------------------------

figure;
mask=nc_read(Gname,'scope_rho');
[Lp,Mp]=size(mask);
if (spherical),
  lon=nc_read(Gname,'lon_rho');
  lat=nc_read(Gname,'lat_rho');
end,
[x,y]=ndgrid(1:Lp,1:Mp);
cmask=mask*2+mod(mod(x,2)+mod(y,2),2);

if (ptype & spherical),
  h=surface(lon,lat,cmask); shading flat;
  colormap([0 1 0;.5 1 0;0 .7 1;.3 0 1]);
  set(gca,'layer','top');
  axis tight;
  grid on;
  hold on
  plot(C.lon,C.lat,'k');
  xlabel('Longitude');
  ylabel('Latitude');
else,
  x=0:1:Lp-1; x=x';
  y=0:1:Mp-1;
  h=image(x,y,cmask','cdatamapping','scaled');
  colormap([0 1 0;.5 1 0;0 .7 1;.3 0 1]);
  set(gca,'YDir','normal','drawmode','fast','layer','top');
  grid on;
  hold on
  plot(C.Icst,C.Jcst,'k');
  xlabel('I-grid');
  ylabel('J-grid');
end,
title('Adjoint Sensitivity Scope Mask on RHO-points');

%---------------------------------------------------------------------------
% Plot scope mask on U-points.
%---------------------------------------------------------------------------

figure;
mask=nc_read(Gname,'scope_u');
[L,Mp]=size(mask);
if (spherical),
  lon=nc_read(Gname,'lon_u');
  lat=nc_read(Gname,'lat_u');
end,
[x,y]=ndgrid(1:L,1:Mp);
cmask=mask*2+mod(mod(x,2)+mod(y,2),2);

if (ptype & spherical),
  h=surface(lon,lat,cmask); shading flat;
  colormap([0 1 0;.5 1 0;0 .7 1;.3 0 1]);
  set(gca,'layer','top');
  axis tight;
  grid on;
  hold on
  plot(C.lon,C.lat,'k');
  xlabel('Longitude');
  ylabel('Latitude');
else,
  x=1:1:L; x=x';
  y=0:1:Mp;
  h=image(x,y,cmask','cdatamapping','scaled');
  colormap([0 1 0;.5 1 0;0 .7 1;.3 0 1]);
  set(gca,'YDir','normal','drawmode','fast','layer','top');
  grid on;
  hold on
  plot(C.Icst,C.Jcst,'k');
  xlabel('I-grid');
  ylabel('J-grid');
end,
title('Adjoint Sensitivity Scope Mask on U-points');

%---------------------------------------------------------------------------
% Plot scope mask on V-points.
%---------------------------------------------------------------------------

figure;
mask=nc_read(Gname,'scope_v');
[Lp,M]=size(mask);
if (spherical),
  lon=nc_read(Gname,'lon_v');
  lat=nc_read(Gname,'lat_v');
end,
[x,y]=ndgrid(1:Lp,1:M);
cmask=mask*2+mod(mod(x,2)+mod(y,2),2);

if (ptype & spherical),
  h=surface(lon,lat,cmask); shading flat;
  colormap([0 1 0;.5 1 0;0 .7 1;.3 0 1]);
  set(gca,'layer','top');
  axis tight;
  grid on;
  hold on
  plot(C.lon,C.lat,'k');
  xlabel('Longitude');
  ylabel('Latitude');
else,
  x=0:1:Lp; x=x';
  y=1:1:M;
  h=image(x,y,cmask','cdatamapping','scaled');
  colormap([0 1 0;.5 1 0;0 .7 1;.3 0 1]);
  set(gca,'YDir','normal','drawmode','fast','layer','top');
  grid on;
  hold on
  plot(C.Icst,C.Jcst,'k');
  xlabel('I-grid');
  ylabel('J-grid');
end,
title('Adjoint Sensitivity Mask of V-points');

return
