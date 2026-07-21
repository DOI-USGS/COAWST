function [Obs]=obs_k2z(ObsData, Hname)

%
% OBS_K2S:  Computes observation depths from fractional vertical coordinate
%
% [Obs]=obs_k2z(ObsData, Hname)
%
% This function computes the observation depths from fractional vertical
% coordinates. It uses ROMS algorithm.
%
% On Input:
%
%    ObsData     ROMS observation or 4D-Var MOD file name (string)
%            or, an existing observation structure (struct array)
%  
%    Hname      ROMS history NetCDF file name (string)
%    
% On Output:
%
%    Obs        Observation structure (struct array)
%
  
% svn $Id$
%=========================================================================%
%  Copyright (c) 2002-2025 The ROMS Group                                 %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.md                            Hernan G. Arango      %
%=========================================================================%
  
% Read in observation structure. 

if (~isstruct(ObsData))
  Obs=obs_read(Obs);
else
  Obs=ObsData;
end

% Set ROMS application depth at W-points.

Zw=depths(Hname,Hname,5);

[Im,Jm,Km]=size(Zw);

Xmin=0;
Xmax=Im-1;
Ymin=0;
Ymax=Jm-1;

%--------------------------------------------------------------------------
% Interpolate vertical fractional coordinates to depths.
%--------------------------------------------------------------------------

for i=1:Obs.Ndatum
  if (((Xmin <= Obs.Xgrid(i)) && (Obs.Xgrid(i) < Xmax)) &&              ...
      ((Ymin <= Obs.Ygrid(i)) && (Obs.Ygrid(i) < Ymax)))
    i1=fix(Obs.Xgrid(i));
    j1=fix(Obs.Ygrid(i));
    i2=i1+1;
    j2=j1+1;
    if (i2 > Im)
      i2=i1;                       % observation at the eastern boundary
    end
    if (j2 > Jm)
      i2=i1;                       % observation at the eastern boundary
    end
    p2=(i2-i1)*(Obs.Xgrid(i)-i1);
    q2=(j2-j1)*(Obs.Ygrid(i)-j1);
    p1=1-p2;
    q1=1-q2;
    w11=p1*q1;
    w21=p2*q1;
    w22=p2*q2;
    w12=p1*q2;
%   ObsDepth=Obs.depth(i);
    ObsDepth=Obs.Zgrid(i);
    if (ObsDepth > 0)
      if (abs(ObsDepth) <  1)
        Linterpolate=false;                   % surface observation
        Zobs(i)=0;
      else
        Linterpolate=true;                    % fractional level
        k1=max(0,fix(ObsDepth-0.5));          % W-point
        k2=min(k1+1,Km);
        r2=(k2-k1)*(ObsDepth-k1);
        r1=1-r2;
      end
    else
      Linterpolate=false;                     % already depths in meters
      Zobs(i)=ObsDepth;
    end
    if (Linterpolate)
      if ((r1+r2) > 0)
        Hmat(1)=w11*r1;                       % every spatial index in Zw
        Hmat(2)=w21*r1;                       % needs to be shifted by 1
        Hmat(3)=w22*r1;                       % because of the zero lower
        Hmat(4)=w12*r1;                       % bound in all dimensions
        Hmat(5)=w11*r2;
        Hmat(6)=w21*r2;
        Hmat(7)=w22*r2;
        Hmat(8)=w12*r2;
        Zobs(i)=Hmat(1)*Zw(i1+1,j1+1,k1+1)+                             ...
                Hmat(2)*Zw(i2+1,j1+1,k1+1)+                             ...
                Hmat(3)*Zw(i2+1,j2+1,k1+1)+                             ...
                Hmat(4)*Zw(i1+1,j2+1,k1+1)+                             ...
                Hmat(5)*Zw(i1+1,j1+1,k2+1)+                             ...
                Hmat(6)*Zw(i2+1,j1+1,k2+1)+                             ...
                Hmat(7)*Zw(i2+1,j2+1,k2+1)+                             ...
                Hmat(8)*Zw(i1+1,j2+1,k2+1);
      else
	Zobs(i)=NaN;
      end
    end
  end
end

Obs.depth=Zobs;

return
