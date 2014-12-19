function [x,y,dhdx,dhde,slope,r]=hslope(fname,iprint,iplot);
 
%
% HSLOPE:  Compute and plot ROMS bathymetry slope
%
% [x,y,dhdx,dhde,slope,r]=hslope(fname,iprint,iplot)
%
% This function computes the bathymetry slope from a ROMS NetCDF file.
%
% On Input:
%
%    fname       NetCDF file name (character string)
%    iprint      Optional, flag for printing out information (integer):
%                  iprint = 0, do not print information
%                  iprint = 1, print information
%    iplot       Optional, flag for plotting (integer):
%                  iplot  = 0, do not plot data
%                  iplot  = 1, plot data
%
% On Output:
%
%    x           X-positions (2D array)
%    y           Y-positions (2D array)
%    dhdx        Bathymetry gradient in the XI-direction (2D array)
%    dhde        Bathymetry gradient in the ETA-direction (2D array)
%    slope       Bathymetry slope (2D array)
%    r           R-value (2D array)
%

% svn $Id: hslope.m 738 2014-10-14 21:49:14Z arango $
%===========================================================================%
%  Copyright (c) 2002-2014 The ROMS/TOMS Group                              %
%    Licensed under a MIT/X style license                                   %
%    See License_ROMS.txt                           Hernan G. Arango        %
%===========================================================================%


%  If not provided, set-up printing and ploting data switch.

if (nargin < 3),
  iplot=1;
end
if (nargin < 2),
  iprint=1;
end

%  Incquire about variables in NetCDF file. Read in bathymetry.

S=nc_vnames(fname);
nvars=length(S.Variables);

got.angle=0;
got.dndx=0;
got.dmde=0;

for n=1:nvars,
  name=char(S.Variables(n).Name);
  switch name
    case 'h'
      h=nc_read(fname,'h',0);
    case 'bath'
      h=nc_read(fname,'bath',1);
    case 'angle'
      got.angle=1;
      angle=nc_read(fname,'angle');
    case 'dndx'
      got.dndx=1;
      dndx=nc_read(fname,'dndx');
    case 'dmde'
      got.dmde=1;
      dmde=nc_read(fname,'dmde');
  end
end

hmin=min(min(h));
hmax=max(max(h));
havg=mean(mean(h));
hmed=median(median(h));

[Lp,Mp]=size(h);
L=Lp-1;
M=Mp-1;

%----------------------------------------------------------------------------
%  Read in curvilinear grid metrics.
%----------------------------------------------------------------------------

pm=nc_read(fname,'pm',0);
pn=nc_read(fname,'pn',0);

dx=1./pm; dx=dx./1000;
dy=1./pn; dy=dy./1000;

dxmin=min(min(dx));
dxmax=max(max(dx));
dxavg=mean(mean(dx));
dxmed=median(median(dx));

dymin=min(min(dy));
dymax=max(max(dy));
dyavg=mean(mean(dy));
dymed=median(median(dy));

%----------------------------------------------------------------------------
%  Read in Land/Sea mask.
%----------------------------------------------------------------------------

got.rmask=false;
got.pmask=false;

for n=1:nvars,
  name=char(S.Variables(n).Name);
  switch (name),
    case 'mask_rho'
      got.rmask=true;
    case 'mask_psi'
      got.pmask=true;
  end
end

if (got.rmask),
 rmask=nc_read(fname,'mask_rho');
else
 rmask=ones([Lp Mp]);
end

if (got.pmask),
 pmask=nc_read(fname,'mask_psi');
else
 pmask=ones([L M]);
end

%----------------------------------------------------------------------------
%  Read in positions.
%----------------------------------------------------------------------------

spherical=nc_read(fname,'spherical');
if (ischar(spherical)),
  if (spherical == 'T' | spherical == 't'),
    spherical=true;
  else
    spherical=false;
  end
end

if (spherical),
  xr=nc_read(fname,'lon_rho',0);
  yr=nc_read(fname,'lat_rho',0);
  x(1:L,1:M)=0.25.*(xr(1:L,1:M)+xr(2:Lp,1:M)+xr(1:L,2:Mp)+xr(2:Lp,2:Mp));
  y(1:L,1:M)=0.25.*(yr(1:L,1:M)+yr(2:Lp,1:M)+yr(1:L,2:Mp)+yr(2:Lp,2:Mp));
else
  xr=nc_read(fname,'x_rho',0);
  yr=nc_read(fname,'y_rho',0);
  x(1:L,1:M)=0.25.*(xr(1:L,1:M)+xr(2:Lp,1:M)+xr(1:L,2:Mp)+xr(2:Lp,2:Mp));
  y(1:L,1:M)=0.25.*(yr(1:L,1:M)+yr(2:Lp,1:M)+yr(1:L,2:Mp)+yr(2:Lp,2:Mp));
end

%----------------------------------------------------------------------------
%  Compute bathymetry gradients at PSI-points.
%----------------------------------------------------------------------------

hx(1:L,1:Mp)=0.5.*(pm(1:L,1:Mp)+pm(2:Lp,1:Mp)).*(h(2:Lp,1:Mp)-h(1:L,1:Mp));
hy(1:Lp,1:M)=0.5.*(pm(1:Lp,1:M)+pm(1:Lp,2:Mp)).*(h(1:Lp,2:Mp)-h(1:Lp,1:M));

dhdx(1:L,1:M)=0.5.*(hx(1:L,1:M)+hx(1:L,2:Mp));
dhde(1:L,1:M)=0.5.*(hy(1:L,1:M)+hy(2:Lp,1:M));

%----------------------------------------------------------------------------
%  Compute R-factor.
%----------------------------------------------------------------------------

hx(1:L,1:Mp)=abs(h(2:Lp,1:Mp)-h(1:L,1:Mp))./(h(2:Lp,1:Mp)+h(1:L,1:Mp));
hy(1:Lp,1:M)=abs(h(1:Lp,2:Mp)-h(1:Lp,1:M))./(h(1:Lp,2:Mp)+h(1:Lp,1:M));

r(1:L,1:M)=max(max(hx(1:L,1:M),hx(1:L,2:Mp)), ...
               max(hy(1:L,1:M),hy(2:Lp,1:M)));

rmin=min(min(r));
rmax=max(max(r));
ravg=mean(mean(r));
rmed=median(median(r));

%----------------------------------------------------------------------------
%  Compute slope PSI-points.
%----------------------------------------------------------------------------

slope=sqrt(dhdx.*dhdx + dhde.*dhde);

smin=min(min(slope));
smax=max(max(slope));
savg=mean(mean(slope));
smed=median(median(slope));
srms=std(std(slope));

%----------------------------------------------------------------------------
%  Compute barotropic timestep (s) limit.
%----------------------------------------------------------------------------

dt=0.5.*1000.0.*sqrt(dx.*dx + dy.*dy)./sqrt(9.807.*h);

ind=find(rmask < 0.5);
if (~isempty(ind)),
  dt(ind)=NaN;
end

dtmin=min(min(dt));
dtmax=max(max(dt));

%----------------------------------------------------------------------------
%  Print information.
%----------------------------------------------------------------------------

if (iprint),
  disp('  ');
  disp(['Minimum dx (km) = ', num2str(dxmin)]);
  disp(['Maximum dx (km) = ', num2str(dxmax)]);
  disp(['Mean    dx (km) = ', num2str(dxavg)]);
  disp(['Median  dx (km) = ', num2str(dxmed)]);

  disp('  ');
  disp(['Minimum dy (km) = ', num2str(dymin)]);
  disp(['Maximum dy (km) = ', num2str(dymax)]);
  disp(['Mean    dy (km) = ', num2str(dyavg)]);
  disp(['Median  dy (km) = ', num2str(dymed)]);

  disp('  ');
  disp(['Minimum dt (s) = ', num2str(dtmin)]);
  disp(['Maximum dt (s) = ', num2str(dtmax)]);

  disp('  ');
  disp(['Minimum depth (m) = ', num2str(hmin)]);
  disp(['Maximum depth (m) = ', num2str(hmax)]);
  disp(['Mean    depth (m) = ', num2str(havg)]);
  disp(['Median  depth (m) = ', num2str(hmed)]);

  disp('  ');
  disp(['Minimum r-value = ', num2str(rmin)]);
  disp(['Maximum r-value = ', num2str(rmax)]);
  disp(['Mean    r-value = ', num2str(ravg)]);
  disp(['Median  r-value = ', num2str(rmed)]);

  disp('  ');
  disp(['Minimum Slope = ', num2str(smin)]);
  disp(['Maximum Slope = ', num2str(smax)]);
  disp(['Mean    Slope = ', num2str(savg)]);
  disp(['Median  Slope = ', num2str(smed)]);
  disp(['    RMS Slope = ', num2str(srms)]);
end,

%----------------------------------------------------------------------------
%  Check if CURVILINEAR flag needs to be activated.
%----------------------------------------------------------------------------

curvilinear=0;

if (got.angle),
  if (min(min(angle)) < 0 | max(max(angle)) > 0),
    curvilinear=1;
  end
end

if (got.dndx),
  if (min(min(dndx)) < 0 | max(max(dndx)) > 0),
    curvilinear=1;
  end
end

if (got.dmde),
  if (min(min(dmde)) < 0 | max(max(dmde)) > 0),
    curvilinear=1;
  end
end

if (iprint),
  if (curvilinear),
    disp('  ');
    disp(['You need to activate the CURVGRID flag in your header file!!!']);
    disp('  ');
  else
    disp('  ');
    disp(['There is not need to activate the CURVGRID flag in your ',...
          'header file...']);
    disp('  ');
  end
end  

%----------------------------------------------------------------------------
%  Plot data.
%----------------------------------------------------------------------------

if (iplot),

 figure;
 pcolor(xr,yr,h); shading interp; colorbar;
 title('Bathymetry (m)');
 xlabel(['Min  = ', num2str(hmin),'   Max  = ', num2str(hmax), ...
         '   Mean  = ', num2str(havg), '   Median  = ', num2str(hmed)]);
 grid;

 figure;
 pcolor(x,y,slope); shading interp; colorbar;
 title('Bathymetric Slope');
 xlabel(['Min  = ', num2str(smin),'   Max  = ', num2str(smax), ...
         '   Mean  = ', num2str(savg), '   Median  = ', num2str(smed)]);
 grid;

 figure;
 pcolor(x,y,r); shading interp; colorbar;
 title('Bathymetry r-value');
 xlabel(['Min  = ', num2str(rmin),'   Max  = ', num2str(rmax), ...
         '   Mean  = ', num2str(ravg), '   Median  = ', num2str(rmed)]);
 grid;

 figure;
 pcolor(xr,yr,dt); shading interp; colorbar;
 title('Barotropic Timestep Limit (seconds)');
 xlabel(['Min  = ', num2str(dtmin),'   Max  = ', num2str(dtmax)]);
 grid;
end

return


