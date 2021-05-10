function R = spherical_grid (D,nx,ny,DX,DY,Clat,Rlon,Rlat,Alpha,Lplt)

%
% SPHERICAL_GRID: Computes a rotated spherical from specified parameters.
%
% R = spherical_grid (D,nx,ny,DX,DY,Clat,Rlon,Rlat,Alpha,Lplt)
%
% Given a larger donor grid (D), this function extracts a rotated spherical
% grid centered at (Tlon, Tlat) and rotated by angle alpha.
%
% On Input:
%
%    D          Donor ROMS grid structure (struct)
%
%    nx         Number of interior RHO-points (Lm) in the XI-direction
%
%    ny         Number of interior RHO-points (Mm) in the ETA-direction
%
%    DX         Rotated grid spacing in the XI-direction (meter)
%
%    DY         Rotated grid spacing in the ETA-direction (meter)
%
%    Clat       Initial center latitude of the Cartesian grid on the
%               Mercator plane on the Greenwich meridian (Clon = 0)
%
%    Rlon       Rotated target grid center longitude (degrees_east)
%
%    Rlat       Rotated target grid center latitude (degrees_north)
%
%    Alpha      Azimuthal angle of grid direction (degrees)
%                 Alpha = 0   means that the positive XI-direction at the
%                             center of the grid points exactly to EAST
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
%                 R.orthogonality      Orthogonality error, (dx/dy)-1
%
% Adapted from A. F. Shchepetkin function "compute_grid.m" from his ROMS
% Grid generation GUI.
%

%
% It producea a logically rectangular curvilinear orthogonal grid with
% minimal distortion. The computation involves 5 stages:
%
% (1) Set up a patch of uniform-resolution Cartesian grid on Mercator plane
%     centered at the user-specified point Clat on the Greenwich meridian
%     (Clon=0), and with boundaries
%
%        xsize = nx * DX   (m)
%        ysize = ny * DY   (m)
%
%     corresponding to the east-west and north-south directions measured
%     along the longitudinal line passing through the center of the patch.
%
%     Specifying Clat=0 makes the patch centered around Equator; a positive
%     (negative) value moves it to the north (south). When projected back
%     to the sphere this makes the grid tapered toward one of the poles
%     with dx, dy ~cos(lat-Clat).  An infinitely small patch is almost
%     rectangular if placed on Equator, while behaving more and more like
%     a sector of polar coordinates if placed closer to one of the poles.
%
% (2) Transfer from Mercator plane to the sphere [x,y] --> [lon,lat].
%
% (3) move the grid toward Equator as a solid object, azimuthally rotate it
%     around it center by angle "Alpha", then moved to the north or south
%     to latitude "Rlat", and finally moved east or west to longitude
%     "Rlon".
%
% (4) Compute metric terms, angle between local positive XI-direction and
%     direction to the east, and orthogonality error. 
%

% svn $Id$
%=========================================================================%
%  Copyright (c) 2002-2019 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                        Alexander F. Shchepetkin %
%=========================================================================%

% Set switch to plot rotated grid.
  
Lplot = false;
if (nargin == 10)
  Lplot = Lplt;
end

% Initialize.

pi = 3.14159265358979323846;
deg2rad = pi/180;                 % degrees to radians
rad2deg = 180/pi;                 % radians to degrees
Eradius = 6371315;                % Earth Radius (m)

xsize = fix(nx) * DX;             % m
ysize = fix(ny) * DY;             % m

Lplot = true;

%--------------------------------------------------------------------------
% Stage 1: Create a patch of Cartesian grid on Mercator plane (x,y).
%          At first, find the extents xmin,xmax,ymin,ymax corresponding
%          to user-specified sizes and center latitude, then create the
%          patch, then extract its perimeter (if for plotting only). 
%--------------------------------------------------------------------------

theta=deg2rad*Clat ;              % latitude of grid center in radians
dtht=0.5*ysize/Eradius ;          % half of north-south extent in radians
ymax=lat_to_eta(theta+dtht);      % north and south bounds of the initial
ymin=lat_to_eta(theta-dtht);      % rectangular domain on Mercator plane

ycent=0.5*(ymin+ymax);            % center on Mercator plane
y_size=ymax-ymin;
dy=y_size/double(ny);             % Temporarily set north and south bounds
ymin=-0.5*y_size;                 % symmetrically about the Equator - their
ymax=+0.5*y_size;                 % original values can be reinstated by
                                  % adding ycent to both.

x_size=xsize/(Eradius*cos(theta));             % east-west extents of
xmin=-0.5*x_size;                              % the grid expressed in
xmax=+0.5*x_size ;                             % radians
dx=x_size/double(nx) ;                    

disp(['ISOTROPY ERROR (dx/dy)-1 = ', num2str(dx/dy -1.) ])

% Create uniform retangular Cartesian grid on Mercator plane at PSI-,
% RHO-, U-, and V-points. The PSI-points have additional points to
% facilitate the computation of the metrics "pm", "pn", and "angle".

[eta_psi, xi_psi]=meshgrid(-1.0:1.0:ny+1.0, -1.0:1.0:nx+1.0);
[eta_rho, xi_rho]=meshgrid(-0.5:1.0:ny+0.5, -0.5:1.0:nx+0.5);
[eta_u,   xi_u  ]=meshgrid(-0.5:1.0:ny+0.5,  0.0:1.0:nx    );
[eta_v,   xi_v  ]=meshgrid( 0.0:1.0:ny    , -0.5:1.0:nx+0.5);

Xp=xmin+dx*xi_psi;
Yp=ymin+dy*eta_psi;  Yp=Yp+ycent;

Xr=xmin+dx*xi_rho;
Yr=ymin+dy*eta_rho;  Yr=Yr+ycent;

Xu=xmin+dx*xi_u;
Yu=ymin+dy*eta_u;    Yu=Yu+ycent;

Xv=xmin+dx*xi_v;
Yv=ymin+dy*eta_v;    Yv=Yv+ycent;

clear eta_psi eta_rho eta_u eta_v xi_psi xi_rho xi_u xi_v

%--------------------------------------------------------------------------
% Stage 2: transfer from Mercator plane to the sphere [x,y] --> [lon,lat]
%          Longitude and latitude are in radians.
%--------------------------------------------------------------------------

latp=asin(tanh(Yp));                   % extended PSI-points
lonp=Xp;

latr=asin(tanh(Yr));                   % RHO-points
lonr=Xr;

latu=asin(tanh(Yu));                   % U-points
lonu=Xu;
					       
latv=asin(tanh(Yv));                   % V-points
lonv=Xv;

clear Xp Xr Xu Xv Yp Yr Yu Yv

%--------------------------------------------------------------------------
% Stage 3: move the grid to the desired place and rotate it azimuthally 
%--------------------------------------------------------------------------

[lonp,latp]=MoveAndTurn(Rlon, Rlat, Alpha, Clat, lonp, latp);
[lonr,latr]=MoveAndTurn(Rlon, Rlat, Alpha, Clat, lonr, latr);
[lonu,latu]=MoveAndTurn(Rlon, Rlat, Alpha, Clat, lonu, latu);
[lonv,latv]=MoveAndTurn(Rlon, Rlat, Alpha, Clat, lonv, latv);

%--------------------------------------------------------------------------
% Stage 4: Compute angles of local grid positive XI-axis relative to east,
%          metric coefficients "pm" and "pn", and check for othogonality.
%--------------------------------------------------------------------------
%
% a11: cos(Lat) * dLon / dX
% a12: cos(Lat) * dLon / dY
% a21: dLat/dX
% a22: dLat/dY
%
% ang = 0.5*cos(latr);                                      % temporally
%
% a11 = ang .* ( lonp(2:ii+1,2:jj+1) -lonp(1:ii,2:jj+1) ...
%               +lonp(2:ii+1,1:jj  ) -lonp(1:ii,1:jj  ) ) ;
%
% a12 = ang .* ( lonp(1:ii,2:jj+1) +lonp(2:ii+1,2:jj+1) ... 
%               -lonp(1:ii,1:jj  ) -lonp(2:ii+1,1:jj  ) ) ;
%
% a21 =  0.5 * ( latp(2:ii+1,2:jj+1) -latp(1:ii,2:jj+1) ... 
%               +latp(2:ii+1,1:jj  ) -latp(1:ii,1:jj  ) ) ;
%
% a22 =  0.5 * ( latp(1:ii,2:jj+1) +latp(2:ii+1,2:jj+1) ... 
%               -latp(1:ii,1:jj  ) -latp(2:ii+1,1:jj  ) ) ;

[Im,Jm]=size(latr);

pm=zeros(Im,Jm);
pn=zeros(Im,Jm);

ang=zeros(Im,Jm);
orterr=zeros(Im,Jm);

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

    dLnY1=lonp(i+1,j+1)-lonp(i+1,j);       % Because of this possibility,
    if (dLnY1 > pi)                        % the elegant Matlab-style array
      dLnY1=dLnY1-2*pi;                    % syntax code above is commented
    elseif (dLnY1 < -pi)                   % out permanently, but still
      dLnY1=dLnY1+2*pi;                    % left here for Matlab lovers
    end                                    % to admire, while quasi-Fortran
    dLnY=lonp(i,j+1)-lonp(i,j);            % code of the left does the job.
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
    orterr(i,j)=ang1-ang2;           % orthogonality error

  end
end

%--------------------------------------------------------------------------
% Load fields into output structure.
%--------------------------------------------------------------------------

R.spherical = true;

if (length(unique(pm(:))) > 1 || length(unique(pn(:))) > 1)
 R.curvilinear = true;
else
 R.curvilinear = false;
end

if (length(unique(ang(:))) > 1 ||  unique(ang(:))  > 0)
  R.vector_rotation = true;
else
  R.vector_rotation = true;
end

R.xl = xsize;
R.el = ysize;

R.lon_psi = rad2deg * lonp(2:end-1,2:end-1);
R.lat_psi = rad2deg * latp(2:end-1,2:end-1);

R.lon_rho = rad2deg * lonr;
R.lat_rho = rad2deg * latr;

R.lon_u = rad2deg * lonu;
R.lat_u = rad2deg * latu;

R.lon_v = rad2deg * lonv;
R.lat_v = rad2deg * latv;

R.pm = pm;
R.pn = pn;

R.angle = ang;
R.orthogonality = orterr;

%--------------------------------------------------------------------------
% Plot rotated grid.
%--------------------------------------------------------------------------

if (Lplot)
  
  figure;
  plot(D.lon_psi(:,1:end),  D.lat_psi(:,1:end),  'k-',                  ...
       D.lon_psi(1:end,:)', D.lat_psi(1:end,:)', 'k-',                  ...
       R.lon_psi(:), R.lat_psi(:), 'r+');
  axis tight;
  xlabel ('\bf\fontsize{18}Lon_{\psi}')
  ylabel ('\bf\fontsize{18}Lat_{\psi}')
  title (['Xcenter = ', num2str(Rlon), blanks(4),                       ...
          'Ycenter = ', num2str(Rlat), blanks(4),                       ...
          'dx = ', num2str(DX/1000), ' km', blanks(4),                  ...
          'dy = ', num2str(DY/1000), ' km', blanks(4),                  ...
          'Alpha = ', num2str(Alpha), blanks(4),                        ...
          'nx = ', num2str(nx), blanks(4),                              ...
          'nx = ', num2str(ny)]);

  figure;
  plot(R.lon_psi(:), R.lat_psi(:), 'r+');

  hold on;
  pcolorjw(D.lon_psi, D.lat_psi, D.mask_psi);
  shading flat; colorbar;
  alpha (0.5);
  axis tight;
  grid on;
  xlabel ('\bf\fontsize{18}Lon_{\psi}')
  ylabel ('\bf\fontsize{18}Lat_{\psi}')
  title (['Xcenter = ', num2str(Rlon), blanks(4),                       ...
          'Ycenter = ', num2str(Rlat), blanks(4),                       ...
          'dx = ', num2str(DX/1000), ' km', blanks(4),                  ...
          'dy = ', num2str(DY/1000), ' km', blanks(4),                  ...
          'Alpha = ', num2str(Alpha), blanks(4),                        ...
          'nx = ', num2str(nx), blanks(4),                              ...
          'ny = ', num2str(ny)]);
  hold off;
  
end

return

%==========================================================================
function eta = lat_to_eta(theta)
%==========================================================================
%
% It converts latitude (radians) into y-coordinate of Mercator pojection.

pi=3.14159265358979323846;

if (-0.5*pi < theta && theta < 0.5*pi)
  cff=sin(theta);
  eta=0.5*log((1.+cff)/(1.-cff));
else
  disp(['lat_to_eta :: theta=', num2str(theta)])
  error('ETA: Latitude range exception.')
end

return

%==========================================================================
function [lon1,lat1] = MoveAndTurn(psid, thetad, alphad, lambdad, lon, lat)
%==========================================================================
%
% The incoming arguments are as follows:
%
%      psid,thetad are geographical lon,lat coordinates of the
%                      desired location of the center of the grid;
%           alphad is azimuthal rotation angle of the grid;
%          lambdad is the initial latitude of the grid center;
%
% Hence the center of the grid travels as follows
%
%         (lon,lat) = (0,lambda) --> (0,0) --> (psi,theta)
%
% Incoming "lon, lat" are expected to be in radians, however all
% four rotation angles are in degrees, so convert them into radians.

pi=3.14159265358979323846;
deg2rad=pi/180.;

psi0=deg2rad*psid;
theta0=deg2rad*thetad;

alpha=deg2rad*alphad;
lambda=deg2rad*lambdad;

% The system of X,Y,Z Cartesian coordinates adopted here is as follows:
%
%   X-axis originates at the center of the Earth and exits through
%          Greenwich meridian at Equator
%
%   Y-axis exists at 90 degrees East (thus, X-Y plane is Equatorial plane)
%
%   Z-axis exists through North pole

[a11,a12,a13,a21,a22,a23,a31,a32,a33]=mtrx_YXY_rotate(theta0,alpha,lambda);

[nxi,neta]=size(lon);

lon1=zeros(size(lon));
lat1=zeros(size(lat));

for j=1:neta
  for i=1:nxi
    csT = cos(lat(i,j));

    xr = csT*cos(lon(i,j));
    yr = csT*sin(lon(i,j));

    zr = sin(lat(i,j));

    x = a11*xr + a12*yr + a13*zr;
    y = a21*xr + a22*yr + a23*zr;
    z = a31*xr + a32*yr + a33*zr;

    if (y < -abs(x))
      psi = -0.5*pi-atan(x/y);   % Here "psi" is defined
    elseif (y > abs(x))          % to be within the range
      psi =  0.5*pi-atan(x/y);   % of
    else                         %    -pi < psi <= pi
      psi = atan(y/x);           %
      if (x < 0.)                % this is the most natural
        if (y < 0.)              % way as the initial grid
          psi = psi-pi;          % is centered around 0E in
        else                     % longitudinal direction.
          psi = psi+pi;
        end
      end
    end

    rd=sqrt(x*x+y*y);
    if (z < -rd)                 % The resultant "tht"
      tht =-0.5*pi-atan(rd/z);   % is within the range
    elseif (z > rd)              % -pi/2 <= psi <= pi/2
      tht = 0.5*pi-atan(rd/z);   % Note that technically
    else                         % speaking, code on the
      tht=atan(z/rd);            % left can handle the
    end                          % situation when rd=0.

    lon1(i,j) = psi + psi0;      % add longitudinal shift
    lat1(i,j) = tht;
  end
end

return

%==========================================================================
function [a11,a12,a13, a21,a22,a23, a31,a32,a33] =                      ...
                                    mtrx_YXY_rotate(theta,alpha,lambda)
%==========================================================================
%
% The following matrix applies three successive turns to vector [x,y,z]:
%
%  1. around Y-axis by angle "lambda", counterclockwise as seen from the
%     top of the arrow looking toward the center of the Earth (i.e, from
%     positive Y-direction to negative);
%
%  2. azimuthal rotation around X-axis by angle "alpha" counterclockwise
%     as seen from the top of the arrow;
%
%  3. around Y-axis again, this time clockwise by angle "theta";
%
% Note that in the case when alpha=0 (hence csA=1 and snA=0 below) this
% matrix degenerates into rotation by angle theta-lambda, as all surviving
% in this case csTh, csLm, snTh, and snLm appear in combinations which add
% up to either cos or sin of difference of the two angles.

csTh = cos(theta);
csA  = cos(alpha);
csLm = cos(lambda);

snTh = sin(theta);
snA  = sin(alpha);
snLm = sin(lambda);

a11 = csTh*csLm + snTh*csA*snLm;
a12 =-snTh*snA;
a13 = csTh*snLm - snTh*csA*csLm;

a21 = snA*snLm;
a22 = csA;
a23 =-snA*csLm;

a31 = snTh*csLm - csTh*csA*snLm;
a32 = csTh*snA;
a33 = snTh*snLm + csTh*csA*csLm;

return
