function S = roms_metrics(Ginp, varargin)

%
% ROMS_METRICS: computes several ROMS grid metrics.
%
% S = roms_metrics (Ginp, Lplot, Ldiff)
%
% Computes the several ROMS grid metrics from the (lon,lat) coordinates.
%
% On Input:
%
%    G          ROMS grid NetCDF filename (string)
%           or, an existing ROMS grid structure (struct array)
%
%    Lplot      Switch to plot metrics (logical; optional)
%
%                 Lplot = false  (default)
%
%    Ldiff      Switch to plot metrics differences: G-S (logical; optional)
%
%                 Ldiff = false  (default)
%
% On Output:
%
%    S          Grid metric structure (struc):
%
%                 S.lon        longitude of RHO-points
%                 S.lat        latitude  of RHO-points
%                 S.x_rho      X-location of RHO-points
%                 S.y_rho      Y-location of RHO-points
%                 S.x_psi      X-location of PSI-points
%                 S.y_psi      Y-location of PSI-points
%                 S.x_u        X-location of U-points
%                 S.y_u        Y-location of U-points
%                 S.x_v        X-location of V-points
%                 S.y_v        Y-location of V-points
%                 S.angle      rotation angle between XI and EAST
%                 S.dx         grid spacing in the X-direction (m)
%                 S.dy         grid spacing in the Y-direction (m)
%                 S.pm         metric factor in the XI-direction (m-1)
%                 S.pn         metric factor in the ETA-direction (m-1)
%                 S.dmde       Inverse metric derivative, d(pm)/d(eta)
%                 S.dndx       Inverse metric derivative, d(pn)/d(xi)
%                 S.f          Coriolis parameter
  
% svn $Id$
%=========================================================================%
%  Copyright (c) 2002-2020 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                                Hernan G. Arango %
%=========================================================================%
  
switch numel(varargin)
  case 0 
    Lplot = false;
    Ldiff = false;
  case 1
    Lplot = varargin{1};
    Ldiff = false;
  case 2
    Lplot = varargin{1};
    Ldiff = varargin{2};
end

if (~isstruct(Ginp))
  G = get_roms_grid(Ginp);
else
  G = Ginp;
end

% Initialize.

pi = 3.14159265358979323846;
deg2rad = pi/180;                           % degrees to radians
rad2deg = 180/pi;                           % radians to degrees
Eradius = 6371315;                          % Earth Radius (m)


S = struct('lon'     , [], 'lat'     , [],                              ...
           'x_rho'   , [], 'y_rho'   , [],                              ...
           'x_psi'   , [], 'y_psi'   , [],                              ...
           'x_u'     , [], 'y_u'     , [],                              ...
           'x_v'     , [], 'y_v'     , [],                              ...
           'dx'      , [], 'dy'      , [],                              ...
           'pm'      , [], 'pn'      , [],                              ...
           'dmde'    , [], 'dndx'    , [],                              ...
           'angle'   , [], 'f'       , [],                              ...
           'xl'      , [], 'el'      , []);
	   
% Build a larger (lon,lat) arrays at PSI points and fill the extra
% edges values.

[Im,Jm]=size(G.lon_rho);

lonp = nan(Im+1,Jm+1);
latp = nan(Im+1,Jm+1);

lonp(2:Im,2:Jm) = G.lon_psi;
latp(2:Im,2:Jm) = G.lat_psi;

lonp(1,2:Jm) = G.lon_v(1,1:Jm-1);           % western edge logitude
latp(1,2:Jm) = G.lat_v(1,1:Jm-1);           % western edge latitude

lonp(2:Im,1) = G.lon_u(1:Im-1,1);           % southern edge logitude
latp(2:Im,1) = G.lat_u(1:Im-1,1);           % southern edge latitude

lonp(end,2:Jm) = G.lon_v(end,1:Jm-1);       % eastern edge longitude
latp(end,2:Jm) = G.lat_v(end,1:Jm-1);       % eastern edge latitude

lonp(2:Im,end) = G.lon_u(1:Im-1,end);       % southern edge logitude
latp(2:Im,end) = G.lat_u(1:Im-1,end);       % southern edge latitude

lonp(1,1) = G.lon_rho(1,1);                 % southwest corner longitude
latp(1,1) = G.lat_rho(1,1);                 % southwest corner latitude

lonp(end,1) = G.lon_rho(end,1);             % southeast corner longitude
latp(end,1) = G.lat_rho(end,1);             % southeast corner latitude

lonp(end,end) = G.lon_rho(end,end);         % northeast corner longitude
latp(end,end) = G.lat_rho(end,end);         % northeast corner latitude

lonp(1,end) = G.lon_rho(1,end);             % northwest corner longitude
latp(1,end) = G.lat_rho(1,end);             % northwest corner latitude

% Convert longitude and latitude to radians.

lonr = deg2rad .* G.lon_rho;
latr = deg2rad .* G.lat_rho;

lonp = deg2rad .* lonp;
latp = deg2rad .* latp;

% Compute inverse grid spacing and angle (From Shchepetkin grid GUI).

pm = zeros(Im,Jm);
pn = zeros(Im,Jm);
ang = zeros(Im,Jm);

for j=1:Jm
  for i=1:Im
    dLnX1=lonp(i+1,j+1)-lonp(i,j+1);        % Most of the time "lonp"
    if (dLnX1 > pi)                         % is expected to be a smooth
      dLnX1=dLnX1-2*pi;                     % function of its indices,
    elseif (dLnX1 < -pi)                    % however it may experience
      dLnX1=dLnX1+2*pi;                     % 2*pi jumps due to periodicity
    end                                     % resulting in large vales of
    dLnX=lonp(i+1,j)-lonp(i,j);             % differences dLnX,dLnY. If
    if (dLnX > pi)                          % this happens, 2*pi is added
      dLnX=dLnX-2*pi;                       % or subtracted as appropriate.
    elseif (dLnX < -pi)
      dLnX=dLnX+2*pi;
    end

    dLnY1=lonp(i+1,j+1)-lonp(i+1,j);
    if (dLnY1 > pi)
      dLnY1=dLnY1-2*pi;
    elseif (dLnY1 < -pi)
      dLnY1=dLnY1+2*pi;
    end
    dLnY=lonp(i,j+1)-lonp(i,j);
    if (dLnY > pi)
      dLnY=dLnY-2*pi;
    elseif (dLnY < -pi)
      dLnY=dLnY+2*pi;
    end

    cff=0.5*cos(latr(i,j));

    a11=cff*(dLnX+dLnX1);                   % cos(Lat)*dLon/dX

    a12=cff*(dLnY+dLnY1);                   % cos(Lat)*dLon/dY

    a21=0.5*(latp(i+1,j+1)-latp(i,j+1)+                                 ...
             latp(i+1,j  )-latp(i,j  )) ;   % dLat/dX

    a22=0.5*(latp(i,j+1)+latp(i+1,j+1)-                                 ... 
             latp(i,j  )-latp(i+1,j  ));    % dLat/dY

                                            % Note that this computation
    pm(i,j)=1./(Eradius*sqrt(a11^2+a21^2)); % of "pm" and "pn" overwrites
    pn(i,j)=1./(Eradius*sqrt(a12^2+a22^2)); % their previously computed
                                            % analytical counterparts.
    if (a21 < -abs(a11))
      ang1=-0.5*pi-atan(a11/a21);           % For a perfectly othogonal
    elseif (a21 > abs(a11))                 % curvilinear grid "ang1" and
      ang1= 0.5*pi-atan(a11/a21);           % "ang2" should be the same.
    else                                    % Here we chose to compute both
      ang1=atan(a21/a11);                   % of them, use their average as
      if (a11 < 0.)                         % the final accepted value for
        if (a21 < 0.)                       % east angle and also compute
          ang1=ang1-pi;                     % and save their difference as
        else                                % the measure for orthogonality
          ang1=ang1+pi;                     % error.
        end
      end
    end

    if (a12 < -abs(a22))
      ang2= 0.5*pi+atan(a22/a12);
    elseif (a12 > abs(a22))
      ang2=-0.5*pi+atan(a22/a12);
    else
      ang2=-atan(a12/a22);
      if (a22 < 0.)
        if (a12 < 0.)
          ang2=ang2+pi;
        else
          ang2=ang2-pi;
        end
      end
    end

    ang(i,j)=0.5*(ang1+ang2);        % the two should be same or very close

  end
end

% Compute grid spacing (m).

dx = 1 ./ pm;
dy = 1 ./ pn;

% Double the values at the edges because the fill (lon,lat) values at for
% the half locations.

dx(1,:)   = 2 .* dx(1,:);
dx(end,:) = 2 .* dx(end,:);

dy(:,1)   = 2 .* dy(:,1);
dy(:,end) = 2 .* dy(:,end);

% Load fields into structure.

S.lon = G.lon_rho;
S.lat = G.lat_rho;
S.angle = ang;
S.dx = dx;
S.dy = dy;
S.pm = 1 ./ dx;
S.pn = 1 ./ dy;

% Compute derivatives of inverse metric factors: d(pm)/d(eta) and
% d(pn)/d(xi).

[Lp,Mp]=size(lonr);
L = Lp-1;             Lm = L-1;
M = Mp-1;             Mm = M-1;

S.dndx(2:L,2:M) = 0.5 .* (1.0 ./ S.pn(3:Lp,2:M ) - 1.0 ./ S.pn(1:Lm,2:M ));
S.dmde(2:L,2:M) = 0.5 .* (1.0 ./ S.pm(2:L ,3:Mp) - 1.0 ./ S.pm(2:L ,1:Mm));

S.dndx(1 , :) = S.dndx(2,:);
S.dndx(Lp, :) = S.dndx(L,:);
S.dndx(: , 1) = S.dndx(:,2);
S.dndx(: ,Mp) = S.dndx(:,M);

S.dndx(1 , 1) = 0.5 * (S.dndx(1 , 2) + S.dndx(2 , 1));
S.dndx(1 ,Mp) = 0.5 * (S.dndx(1 , M) + S.dndx(2 ,Mp));
S.dndx(Lp, 1) = 0.5 * (S.dndx(Lp, 2) + S.dndx(L , 1));
S.dndx(Lp,Mp) = 0.5 * (S.dndx(L ,Mp) + S.dndx(Lp, M));

S.dmde(1 , :) = S.dmde(2,:);
S.dmde(Lp, :) = S.dmde(L,:);
S.dmde(: , 1) = S.dmde(:,2);
S.dmde(: ,Mp) = S.dmde(:,M);

S.dmde(1 , 1) = 0.5 * (S.dmde(1 , 2) + S.dmde(2 , 1));
S.dmde(1 ,Mp) = 0.5 * (S.dmde(1 , M) + S.dmde(2 ,Mp));
S.dmde(Lp, 1) = 0.5 * (S.dmde(Lp, 2) + S.dmde(L , 1));
S.dmde(Lp,Mp) = 0.5 * (S.dmde(L ,Mp) + S.dmde(Lp, M));

% Compute Cartesian coordinates.

x_rho = zeros(size(dx));
for j=1:Mp
  x_rho(1,j) = -dx(1,j);
  for i=1:L
    x_rho(i+1,j) = x_rho(i,j) + dx(i+1,j);
  end
  x_rho( 1,j) = x_rho( 1,j) + 0.5*dx( 1,j);
  x_rho(Lp,j) = x_rho(Lp,j) - 0.5*dx(Lp,j);
end

y_rho = zeros(size(dy));
for i=1:Lp
  y_rho(i,1) = -dy(i,1);
  for j=1:M
    y_rho(i,j+1) = y_rho(i,j) + dy(i,j+1);
  end
  y_rho(i,1 ) = y_rho(i,1 ) + 0.5*dy(i,1 );  
  y_rho(i,Mp) = y_rho(i,Mp) - 0.5*dy(i,Mp);
end
  
S.x_rho = x_rho;
S.y_rho = y_rho;

S.x_psi = 0.25.*(x_rho(1:L,1:M ) + x_rho(2:Lp,1:M ) +                   ...
                 x_rho(1:L,2:Mp) + x_rho(2:Lp,2:Mp));
S.y_psi = 0.25.*(y_rho(1:L,1:M ) + y_rho(2:Lp,1:M ) +                   ...
                 y_rho(1:L,2:Mp) + y_rho(2:Lp,2:Mp));

S.x_u   = 0.5.*(x_rho(1:L,1:Mp) + x_rho(2:Lp,1:Mp));
S.y_u   = 0.5.*(y_rho(1:L,1:Mp) + y_rho(2:Lp,1:Mp));

S.x_v   = 0.5.*(x_rho(1:Lp,1:M) + x_rho(1:Lp,2:Mp));
S.y_v   = 0.5.*(y_rho(1:Lp,1:M) + y_rho(1:Lp,2:Mp));

S.xl    = max(S.x_psi(:));
S.el    = max(S.y_psi(:));

% Compute Coriolis parameter.

omega = 7.2921D-5;              % Earth's rotation (radians/s)

S.f = 2.*omega.*sin(G.lat_rho.*deg2rad);

% If requested, plot metrics.

if (Lplot)
  figure;
  pcolor(S.lon, S.lat, 0.001.*S.dx);
  colorbar; title('DX (km)')

  figure;
  pcolor(S.lon, S.lat, 0.001.*S.dy);
  colorbar; title('DY (km)')
  
  figure;
  pcolor(S.lon, S.lat, S.pm);
  colorbar; title('Inverse Metric factor pm (m-1)')

  figure;
  pcolor(S.lon, S.lat, S.pn);
  colorbar; title('Inverse Metric factor pn (m-1)')

  figure;
  pcolor(S.lon, S.lat, S.dmde);
  colorbar; title('Metric derivative d(pm)/d(eta)')

  figure;
  pcolor(S.lon, S.lat, S.dndx);
  colorbar; title('Metric derivative d(pn)/d(xi)')
  
  figure;
  pcolor(S.lon, S.lat, rad2deg.*S.angle);
  colorbar; title('Rotation Angle between XI and EAST')
end  

if (Ldiff)
  figure;
  pcolor(S.lon, S.lat, 0.001.*((1./G.pm)-S.dx));
  colorbar; title('G.DX - S.DX (km)')

  figure;
  pcolor(S.lon, S.lat, 0.001.*((1./G.pn)-S.dy));
  colorbar; title('G.DY - S.DY (km)')
  
  figure;
  pcolor(S.lon, S.lat, (G.pm-S.pm));
  colorbar; title('Inverse Metric factor pm (m-1)')

  figure;
  pcolor(S.lon, S.lat, (G.pn-S.pn));
  colorbar; title('Inverse Metric factor pn (m-1)')

  figure;
  pcolor(S.lon, S.lat, (G.dmde-S.dmde));
  colorbar; title('Metric derivative G.dmde - S.dmde')

  figure;
  pcolor(S.lon, S.lat, (G.dndx-S.dndx));
  colorbar; title('Metric derivative G.dndx - S.dndx')
  
  figure;
  pcolor(S.lon, S.lat, rad2deg.*(G.angle-S.angle));
  colorbar; title('Rotation Angle G.angle - S.angle')
end

return
