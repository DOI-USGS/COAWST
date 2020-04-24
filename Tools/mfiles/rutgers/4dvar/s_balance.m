function [deltaS_b]=s_balance(K,T,S,deltaT,ml_depth,dTdz_min);

%
% S_BALANCE:  Computes the balanced salinity anomaly
%
% [deltaS_b]=s_balance(K,T,S,deltaT,ml_depth,dTdz_min)
%
% This function computes the balanced salinity anomaly for the
% error covariance balance operator:
%
%          deltaS_b = fac * dS/dz * dz/dT * deltaT
%
% where fac is a coefficient that depends on the mixed layer depth:
%
%               fac = 1.0 - EXP (z_r / ml_depth)
%
% On Input:
%
%    K           Balance operator structure array:
%
%                  K.ml_depth  Default mixed-layer depth (100 m)
%                  K.dTdz_min  Default minimun dT/dz allowed (0.001 Celsius/m)
%                  K.Norder    Shapiro filter order
%                  K.Fcoef     Shapiro filter coefficients
%                  K.Hz        Vertical level thicknesses (m)
%                  K.Zr        Depths at vertical RHO-points (m, negative)
%                  K.Zw        Depths at vertical W-points   (m, negative)
%                  
%    T           Temperature (3D array)
%    S           Salinity (3D array)
%    deltaT      Temperature anomaly (T-Tavg) from a time mean (3d array)
%    ml_depth    Mixed-layer depth (m, positive), scalar, OPTIONAL
%    dTdz_min    Minimum dT/dz value allowed (Celsius/m) scalar, OPTIONAL
%
% On Output:
%
%    deltaS_b    Balanced salinity anomaly.
%

% svn $Id: s_balance.m 996 2020-01-10 04:28:56Z arango $
%===========================================================================%
%  Copyright (c) 2002-2020 The ROMS/TOMS Group                              %
%    Licensed under a MIT/X style license                                   %
%    See License_ROMS.txt                           Hernan G. Arango        %
%===========================================================================%

% If applicable, set default parameters.

if (nargin < 6),
  dTdz_min=K.dTdz_min;
else,
  dTdz_min=abs(dTdz_min);
end,

if (nargin < 5),
  ml_depth=K.ml_depth;
else,
  ml_depth=abs(ml_depth);
end,

% Compute temperature (dTdz) and salinity (dSdz) shears at vertical
% W-points. Use parabolic splines.

[Lp,Mp,N]=size(K.Zr);

dTdz=zeros([Lp,Mp,N+1]);
dSdz=zeros([Lp,Mp,N+1]);
FC=zeros([Lp,Mp,N+1]);

for k=1:N-1,
  cff=1.0./(2.0.*K.Hz(:,:,k+1)+K.Hz(:,:,k).*(2.0-FC(:,:,k)));
  FC(:,:,k+1)=cff.*K.Hz(:,:,k+1);
  dTdz(:,:,k+1)=cff.*(6.0.*(T(:,:,k+1)-T(:,:,k))-K.Hz(:,:,k).*dTdz(:,:,k));
  dSdz(:,:,k+1)=cff.*(6.0.*(S(:,:,k+1)-S(:,:,k))-K.Hz(:,:,k).*dSdz(:,:,k));
end,

for k=N:-1:1,
  dTdz(:,:,k)=dTdz(:,:,k)-FC(:,:,k).*dTdz(:,:,k+1);
  dSdz(:,:,k)=dSdz(:,:,k)-FC(:,:,k).*dSdz(:,:,k+1);
end,

clear FC

% Compute dz/dT at vertical RHO-points.  If dT/dz is less than a threshold
% value (dTdz_min), set dz/dT = 0.

dTdz_r=0.5.*(dTdz(1:Lp,1:Mp,1:N)+dTdz(1:Lp,1:Mp,2:N+1));

ind=find(dTdz_r < dTdz_min);
dTdz_r(ind)=1;

dzdT=1.0./dTdz_r;
dzdT(ind)=0.0;

clear dTdz_r ind

% Compute dS/dT at RHO-points.

dSdT=0.5.*(dSdz(1:Lp,1:Mp,1:N)+dSdz(1:Lp,1:Mp,2:N+1)).*dzdT;

clear dTdz dSdz

% Shapiro filter dS/dT.

dSdT_filter=zeros(size(dSdT));

for order=1:K.Norder/2,
  if (order ~= K.Norder/2),
    dSdT_filter(:,:,1)=2.0*(dSdT(:,:,1)-dSdT(:,:,2));
    dSdT_filter(:,:,N)=2.0*(dSdT(:,:,N)-dSdT(:,:,N-1));
  else,
    dSdT_filter(:,:,1)=0.0;
    dSdT_filter(:,:,N)=0.0;
  end,

  dSdT_filter(1:Lp,1:Mp,2:N-1)=2.0.*dSdT(1:Lp,1:Mp,2:N-1)- ...
                               dSdT(1:Lp,1:Mp,1:N-2)- ...
			       dSdT(1:Lp,1:Mp,3:N);
  
  dSdT=dSdT-K.Fcoef(K.Norder/2).*dSdT_filter;
end,

clear dSdT_filter

% Compute balance salinity. Notice that multiplication to deltaT
% will be done elsewhere.

deltaS_b=(1.0-exp(K.Zr./ml_depth)).*dSdT.*deltaT;

return
