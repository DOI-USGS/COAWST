function [lon,lat,z]=ijk_lonlatz(Hinp, I, J, K)
  
%
% IJK_LONLATZ:  Converts (I,J,K) fractional coordinates to (lon,lat,z)
%
% [lon,lat,z]=ijk_lonlatz(Hinp, I, J, K)
%
% Converts observation locations from ROMS (I,J,K) fractional coordinates
% to geographical coordinates (lon,lat,z).  The observation locations are
% assumed to be at RHO-points.  The depths are time-independent (zeta=0).
%
% On Input:
%
%    Hinp              ROMS History NetCDF filename (string)
%                  or, an existing grid structure (struct array)
%
%    I                 Observation I-fractional coordinate (vector)
%
%    J                 Observation J-fractional coordinate (vector)
%
%    K                 Observation K-fractional coordinate (vector)
%
% On Output:
%
%    lon               Observation longitude (vector; degrees_east)
%
%    lat               Observation latitude (vector; degrees_north)
%
%    z                 Observation depth (vector; meters, negative)
%

% svn $Id$
%=========================================================================%
%  Copyright (c) 2002-2025 The ROMS Group                                 %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.md                            Hernan G. Arango      %
%=========================================================================%

% Initialize.

if (~isstruct(Hinp))
  G = get_roms_grid(Hinp);
else
  G = Hinp;
end

[Lp,Mp,Np] = size(G.z_w);          % depths at W-points (RHO-cell) 
method = 'linear';                 % interpolation method

Xmin = 0;
Xmax = Lp-1;
Ymin = 0;
Ymax = Mp-1;

lon = NaN(size(I));
lat = NaN(size(J));
z   = NaN(size(K));

%--------------------------------------------------------------------------
% Interpolate from fractional observation (I,J) coordinates to (lon,lat).
% Notice that in 4D-Var, the arbitrary origin of the fractional coordinates
% is (0.0, 0.0).
%--------------------------------------------------------------------------

[Jgrid, Igrid]=meshgrid(0:Mp-1, 0:Lp-1);


F = griddedInterpolant(Igrid, Jgrid, G.lon_rho, method);

                       lon = F(I, J);
F.Values = G.lat_rho;  lat = F(I, J);

%--------------------------------------------------------------------------
% Interpolate vertical coordinates to depths. Same way as ROMS operator.
%--------------------------------------------------------------------------

for n=1:length(K)
  if (((Xmin <= I(n)) && (I(n) < Xmax)) &&                              ...
      ((Ymin <= J(n)) && (J(n) < Ymax)))
    i1=fix(I(n));
    j1=fix(J(n));
    i2=i1+1;
    j2=j1+1;
    if (i2 > Lp)
      i2=i1;                       % observation at the eastern boundary
    end
    if (j2 > Mp)
      i2=i1;                       % observation at the eastern boundary
    end
    p2=(i2-i1)*(I(n)-i1);
    q2=(j2-j1)*(J(n)-j1);
    p1=1-p2;
    q1=1-q2;
    w11=p1*q1;
    w21=p2*q1;
    w22=p2*q2;
    w12=p1*q2;
    depth=K(n);
    if (depth > 0)
      if (abs(depth) <  1)
        Linterpolate=false;                   % surface observation
        z(n)=0;
      else
        Linterpolate=true;                    % fractional level
        k1=max(0,fix(depth-0.5));             % W-point
        k2=min(k1+1,Np);
        r2=(k2-k1)*(depth-k1);
        r1=1-r2;
      end
    else
      Linterpolate=false;                     % already depths in meters
      z(n)=depth;
    end
    if (Linterpolate)
      if ((r1+r2) > 0)
        Hmat(1)=w11*r1;                       % every spatial index in z_w
        Hmat(2)=w21*r1;                       % needs to be shifted by 1
        Hmat(3)=w22*r1;                       % because of the zero lower
        Hmat(4)=w12*r1;                       % bound in all dimensions
        Hmat(5)=w11*r2;
        Hmat(6)=w21*r2;
        Hmat(7)=w22*r2;
        Hmat(8)=w12*r2;
        z(n)=Hmat(1)*G.z_w(i1+1,j1+1,k1+1)+                             ...
             Hmat(2)*G.z_w(i2+1,j1+1,k1+1)+                             ...
             Hmat(3)*G.z_w(i2+1,j2+1,k1+1)+                             ...
             Hmat(4)*G.z_w(i1+1,j2+1,k1+1)+                             ...
             Hmat(5)*G.z_w(i1+1,j1+1,k2+1)+                             ...
             Hmat(6)*G.z_w(i2+1,j1+1,k2+1)+                             ...
             Hmat(7)*G.z_w(i2+1,j2+1,k2+1)+                             ...
             Hmat(8)*G.z_w(i1+1,j2+1,k2+1);
      else
        z(n)=NaN;
      end
    end
  end
end

return
