function [z]=set_depth(Vtransform, Vstretching, ...
                       theta_s, theta_b, hc, N, ...
                       igrid, h, zeta, report)
%
% SET_DEPTH:  Compute ROMS grid depth from vertical stretched variables
%
% [z]=set_depth(Vtransform, Vstretching, theta_s, theta_b, hc, N, ...
%               igrid, h, zeta)
%
% Given a batymetry (h), free-surface (zeta) and terrain-following
% parameters, this function computes the 3D depths for the requested
% C-grid location. If the free-surface is not provided, a zero value
% is assumed resulting in unperturb depths.  This function can be
% used when generating initial conditions or climatology data for
% an application. Check the following link for details:
%
%    https://www.myroms.org/wiki/index.php/Vertical_S-coordinate
%
% On Input:
%
%    Vtransform    Vertical transformation equation:
%
%                    Vtransform = 1,   original transformation
%
%                    z(x,y,s,t)=Zo(x,y,s)+zeta(x,y,t)*[1+Zo(x,y,s)/h(x,y)]
%
%                    Zo(x,y,s)=hc*s+[h(x,y)-hc]*C(s)
%
%                    Vtransform = 2,   new transformation
%
%                    z(x,y,s,t)=zeta(x,y,t)+[zeta(x,y,t)+h(x,y)]*Zo(x,y,s)
%
%                    Zo(x,y,s)=[hc*s(k)+h(x,y)*C(k)]/[hc+h(x,y)]
%
%    Vstretching   Vertical stretching function:
%                    Vstretching = 1,  original (Song and Haidvogel, 1994)
%                    Vstretching = 2,  A. Shchepetkin (UCLA-ROMS, 2005)
%                    Vstretching = 3,  R. Geyer BBL refinement
%                    Vstretching = 4,  A. Shchepetkin (UCLA-ROMS, 2010)
%
%    theta_s       S-coordinate surface control parameter (scalar)
%
%    theta_b       S-coordinate bottom control parameter (scalar)
%
%    hc            Width (m) of surface or bottom boundary layer in which
%                    higher vertical resolution is required during
%                    stretching (scalar)
%
%    N             Number of vertical levels (scalar)
%
%    igrid         Staggered grid C-type (integer):
%                    igrid=1  => density points
%                    igrid=2  => streamfunction points
%                    igrid=3  => u-velocity points
%                    igrid=4  => v-velocity points
%                    igrid=5  => w-velocity points
%
%    h             Bottom depth, 2D array at RHO-points (m, positive),
%                    h(1:Lp+1,1:Mp+1)
%
%    zeta          Free-surface, 2D array at RHO-points (m), OPTIONAL,
%                    zeta(1:Lp+1,1:Mp+1)
%
%    report        Flag to report detailed information (OPTIONAL):
%                    report = 0,       do not report
%                    report = 1,       report information
%
% On Output:
%
%    z             Depths (m, negative), 3D array
%

% svn $Id: set_depth.m 711 2014-01-23 20:36:13Z arango $
%=========================================================================%
%  Copyright (c) 2002-2014 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%

z=[];

%--------------------------------------------------------------------------
%  Set several parameters.
%--------------------------------------------------------------------------

if (nargin < 8),
  disp(' ');
  disp('*** Error:  SET_DEPTH - too few arguments.');
  disp(['                     number of supplied arguments: ',          ...
       num2str(nargin)]);
  disp('                     number of required arguments: 8');
  disp(' ');
  return
end

if (Vtransform < 1 || Vtransform > 2),
  disp(' ');
  disp(['*** Error:  SET_DEPTH - Illegal parameter Vtransform = '       ...
	    num2str(Vtransfrom)]);
  disp(' ');
  return
end

if (Vstretching < 1 || Vstretching > 4),
  disp(' ');
  disp(['*** Error:  SET_DEPTH - Illegal parameter Vstretching = '      ...
        num2str(Vstretching)]);
  disp(' ');
  return
end

if (hc > min(min(h)) && Vtransform == 1),
  disp(' ');
  disp(['*** Error:  SET_DEPTH - critical depth exceeds minimum'        ...
        ' bathymetry value.']);
  disp(['                        Vtransform = ', num2str(Vtransform)]);
  disp(['                        hc         = ', num2str(hc)]);
  disp(['                        hmax       = ', num2str(min(min(h)))]);
  disp(' ');
  return
end

if (nargin < 9),
  zeta=zeros(size(h));
end

if (nargin < 10),
  report=1;
end

Np=N+1;
[Lp Mp]=size(h);
L=Lp-1;
M=Mp-1;

hmin=min(min(h));
hmax=max(max(h));

%--------------------------------------------------------------------------
% Compute vertical stretching function, C(k):
%--------------------------------------------------------------------------

if (report),
  disp(' ');
  if (Vtransform == 1),
    disp(['Vtransform  = ',num2str(Vtransform), '   original ROMS']);
  elseif (Vtransform == 2),
    disp(['Vtransform  = ',num2str(Vtransform), '   ROMS-UCLA']);
  end

  switch ( igrid ),
    case 1
      disp(['   igrid    = ',num2str(igrid),                            ...
            '   at horizontal RHO-points']);
    case 2
      disp(['   igrid    = ',num2str(igrid),                            ...
            '   at horizontal PSI-points']);
    case 3
      disp(['   igrid    = ',num2str(igrid),                            ...
            '   at horizontal U-points']);
    case 4
      disp(['   igrid    = ',num2str(igrid),                            ...
            '   at horizontal V-points']);
    case 5
      disp(['   igrid    = ',num2str(igrid),                            ...
            '   at horizontal RHO-points']);
  end
end

if (igrid == 5),
  kgrid=1;
else
  kgrid=0;
end

[s,C]=stretching(Vstretching, theta_s, theta_b, hc, N, kgrid, report);

%--------------------------------------------------------------------------
%  Average bathymetry and free-surface at requested C-grid type.
%--------------------------------------------------------------------------

switch ( igrid ),
  case 1
    hr=h;
    zetar=zeta;
  case 2
    hp=0.25.*(h(1:L,1:M)+h(2:Lp,1:M)+h(1:L,2:Mp)+h(2:Lp,2:Mp));
    zetap=0.25.*(zeta(1:L,1:M )+zeta(2:Lp,1:M )+                        ...
                 zeta(1:L,2:Mp)+zeta(2:Lp,2:Mp));
  case 3
    hu=0.5.*(h(1:L,1:Mp)+h(2:Lp,1:Mp));
    zetau=0.5.*(zeta(1:L,1:Mp)+zeta(2:Lp,1:Mp));
  case 4
    hv=0.5.*(h(1:Lp,1:M)+h(1:Lp,2:Mp));
    zetav=0.5.*(zeta(1:Lp,1:M)+zeta(1:Lp,2:Mp));
  case 5
    hr=h;
    zetar=zeta;
end,

%--------------------------------------------------------------------------
% Compute depths (m) at requested C-grid location.
%--------------------------------------------------------------------------

if (Vtransform == 1),
  switch ( igrid ),
    case 1
      for k=1:N,
	z0=(s(k)-C(k))*hc + C(k).*hr;
        z(:,:,k)=z0 + zetar.*(1.0 + z0./hr);
      end
    case 2
      for k=1:N,
        z0=(s(k)-C(k))*hc + C(k).*hp;
        z(:,:,k)=z0 + zetap.*(1.0 + z0./hp);
      end
    case 3
      for k=1:N,
        z0=(s(k)-C(k))*hc + C(k).*hu;
        z(:,:,k)=z0 + zetau.*(1.0 + z0./hu);
      end
    case 4
      for k=1:N,
        z0=(s(k)-C(k))*hc + C(k).*hv;
        z(:,:,k)=z0 + zetav.*(1.0 + z0./hv);
      end
    case 5
      z(:,:,1)=-hr;
      for k=2:Np,
        z0=(s(k)-C(k))*hc + C(k).*hr;
        z(:,:,k)=z0 + zetar.*(1.0 + z0./hr);
      end
  end,
elseif (Vtransform == 2),
  switch ( igrid ),
    case 1
      for k=1:N,
        z0=(hc.*s(k)+C(k).*hr)./(hc+hr);
        z(:,:,k)=zetar+(zeta+hr).*z0;
      end
    case 2
      for k=1:N,
        z0=(hc.*s(k)+C(k).*hp)./(hc+hp);
        z(:,:,k)=zetap+(zetap+hp).*z0;
      end
    case 3
      for k=1:N,
        z0=(hc.*s(k)+C(k).*hu)./(hc+hu);
        z(:,:,k)=zetau+(zetau+hu).*z0;
      end
    case 4
      for k=1:N,
        z0=(hc.*s(k)+C(k).*hv)./(hc+hv);
        z(:,:,k)=zetav+(zetav+hv).*z0;
      end
    case 5
      for k=1:Np,
        z0=(hc.*s(k)+C(k).*hr)./(hc+hr);
        z(:,:,k)=zetar+(zetar+hr).*z0;
      end
  end
end

return
