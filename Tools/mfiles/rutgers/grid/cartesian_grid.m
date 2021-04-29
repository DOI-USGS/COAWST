function R = rotate_cartesian (D, Im, Jm, dx, dy, Xc, Yc, Theta, Lplt)

%
% ROTATE_CARTESIAN: Computes a rotated Cartesian grid
%
% R = rotate_cartesian (D, Im, Jm, dx, dy, Xc, Yc, theta, Lplt)
%
% Given a larger donor grid (D), this function extracts a rotated Cartesian
% grid centered at (Xc, Yc) and rotated by angle "Theta". The larger donor
% grid is used during plotting to display the location of rotated grid.
%
% On Input:
%
%    D          Donor ROMS grid NetCDF filename (string)
%           or, an existing ROMS grid structure (struct array)
%
%    Im         Number of PSI-points (Lm+1) in the XI-direction
%
%    Jm         Number of PSI-points (Mm+1) in the ETA-direction
%
%    dx         Rotated grid spacing in the XI-direction (meter)
%
%    dy         Rotated grid spacing in the ETA-direction (meter)
%
%    Xc         Rotated target grid center XI-coordinate (m)
%
%    Yc         Rotated target grid center ETA-coordinate (m)
%
%    Theta      Rotation angle between positive XI-axis relative to EAST
%                (degrees)
%
%    Lplt       Switch to plot computed rotated grid (OPTIONAL, logical,
%                 default = false)
%
% On Output:
%
%    R          Rotated grid structure (struct)
%  
%                 R.spherical          Spherical grid switch
%                 R.curvilinear        curvilinear switch
%                 R.vector_rotation    Vector rotation switch
%                 R.lon_psi            PSI-points longitude
%                 R.lat_psi            PSI-points latitude
%                 R.lon_rho            RHO-points longitude
%                 R.lat_rho            RHO-points latitude
%                 R.lon_u              U-points longitude
%                 R.lat_u              U-points latitude
%                 R.lon_v              V-points longitude
%                 R.lat_v              V-points latitude
%                 R.pm                 Inverse grid XI-spacing
%                 R.pn                 Inverse grid ETA-spacing
%                 R.angle              Rotation angle between XI and EAST
%
% For Example, the rotated grid in the LAKE_JERSEY application is
% extracted using the following command:
%
%    R = rotate_cartesian('lake_jersey_grd_a.nc', 164, 85, 100, 100, ...
%                         14000, 30000, -24.5, true);
%

% svn $Id$
%=========================================================================%
%  Copyright (c) 2002-2019 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%

% Set switch to plot rotated grid.
  
Lplot = false;
if (nargin == 9)
  Lplot = Lplt;
end

% Initialize.

pi = 3.14159265358979323846;
deg2rad = pi/180;                 % degrees to radians
theta = Theta * deg2rad;          % rotation angle in radians

% Compute half of the extend in the retagular domain in the east-west and
% north-south directions.

Xmin = Xc - (Im/2) * dx;
Xmax = Xc + (Im/2) * dx;
Ymin = Yc - (Jm/2) * dy;
Ymax = Yc + (Jm/2) * dy;

Xsize = Xmax - Xmin;                  % meter
Ysize = Ymax - Ymin;                  % meters

Xoff = 0.5 * dx;                      % half X-cell offset
Yoff = 0.5 * dy;                      % half Y-cell offset
  
% Compute rotated Cartisian coordinates (meter) centered on (Xc,Yc).
% Notice that the grid points in the Cartesian plane are shifted
% [X-Xc; Y-Yc] so the desired center of rotation is moved to the
% origin (0,0). After rotation, the grid points are  moved back to
% the desired center.

CosAngle = cos(theta);
SinAngle = sin(theta);

[Yp, Xp] = meshgrid(Ymin:dy:Ymax, Xmin:dx:Xmax);
x_psi = (Xp-Xc)*CosAngle + (Yp-Yc)*SinAngle + Xc;
y_psi = (Yp-Yc)*CosAngle - (Xp-Xc)*SinAngle + Yc;

[Yr, Xr] = meshgrid(Ymin-Yoff:dy:Ymax+Yoff, Xmin-Xoff:dx:Xmax+Xoff);
x_rho = (Xr-Xc)*CosAngle + (Yr-Yc)*SinAngle + Xc;
y_rho = (Yr-Yc)*CosAngle - (Xr-Xc)*SinAngle + Yc;

[Yu, Xu] = meshgrid(Ymin-Yoff:dy:Ymax+Yoff, Xmin:dx:Xmax);
x_u = (Xu-Xc)*CosAngle + (Yu-Yc)*SinAngle + Xc;
y_u = (Yu-Yc)*CosAngle - (Xu-Xc)*SinAngle + Yc;

[Yv, Xv] = meshgrid(Ymin:dy:Ymax, Xmin-Xoff:dx:Xmax+Xoff);
x_v = (Xv-Xc)*CosAngle + (Yv-Yc)*SinAngle + Xc;
y_v = (Yv-Yc)*CosAngle - (Xv-Xc)*SinAngle + Yc;

% Set inverse grid spacing and rotation angle.

pm = ones(size(x_rho)) * dx;          % 1/m
pn = ones(size(x_rho)) * dy;          % 1/m

angle = ones(size(x_rho)) * theta;    % radians

%--------------------------------------------------------------------------
% Load fields into output structure.
%--------------------------------------------------------------------------

R.spherical = true;

if (length(unique(pm(:))) > 1 || length(unique(pn(:))) > 1)
 R.curvilinear = true;
else
 R.curvilinear = false;
end

if (length(unique(angle(:))) > 1 ||  unique(angle(:))  > 0)
  R.vector_rotation = true;
else
  R.vector_rotation = true;
end

R.xl = Xsize;
R.el = Ysize;

R.x_psi = x_psi;
R.y_psi = y_psi;

R.x_rho = x_rho;
R.y_rho = y_rho;

R.x_u = x_u;
R.y_u = y_u;

R.x_v = x_v;
R.y_v = y_v;

R.pm = pm;
R.pn = pn;

R.angle = angle;

%--------------------------------------------------------------------------
% Plot rotated grid.
%--------------------------------------------------------------------------

if (Lplot)

  if (~isstruct(D)),
    ncname = D;
    D = get_roms_grid(ncname);
  end

  figure;
  plot(D.x_psi(:,1:end)/1000,  D.y_psi(:,1:end)/1000,  'k-',            ...
       D.x_psi(1:end,:)'/1000, D.y_psi(1:end,:)'/1000, 'k-',            ...
       R.x_psi(:)/1000, R.y_psi(:)/1000, 'r+');
  axis tight;
  axis equal;
  xlabel ('\bf\fontsize{18}X_{\psi}')
  ylabel ('\bf\fontsize{18}Y_{\psi}')
  title (['Xc = ', num2str(Xc/1000), ' km', blanks(4),                  ...
          'Yc = ', num2str(Yc/1000), ' km', blanks(4),                  ...
          'dx = ', num2str(dx/1000), ' km', blanks(4),                  ...
          'dy = ', num2str(dy/1000), ' km', blanks(4),                  ...
          'theta = ', num2str(Theta), blanks(4),                        ...
          'Im = ', num2str(Im), blanks(4),                              ...
          'Jm = ', num2str(Jm)]);

  figure;
  plot(R.x_psi(:)/1000, R.y_psi(:)/1000, 'r+');

  hold on;
  pcolorjw(D.x_psi/1000, D.y_psi/1000, D.mask_psi);
  shading flat; colorbar;
  alpha (0.5);
  axis tight;
  axis equal;
  grid on;
  xlabel ('\bf\fontsize{18}X_{\psi}')
  ylabel ('\bf\fontsize{18}X_{\psi}')
  title (['Xc = ', num2str(Xc/1000), ' km', blanks(4),                  ...
          'Yc = ', num2str(Yc/1000), ' km', blanks(4),                  ...
          'dx = ', num2str(dx/1000), ' km', blanks(4),                  ...
          'dy = ', num2str(dx/1000), ' km', blanks(4),                  ...
          'theta = ', num2str(Theta), blanks(4),                        ...
          'Im = ', num2str(Im), blanks(4),                              ...
          'Jm = ', num2str(Jm)]);
  hold off;
  
end

return
