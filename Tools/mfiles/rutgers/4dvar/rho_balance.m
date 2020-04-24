function [deltaR_b]=rho_balance(K,deltaT,deltaS_b);

%
% RHO_BALANCE:  Computes the balanced density anomaly
%
% [deltaR_b]=rho_balance(K,deltaT,deltaS_b)
%
% This function computes the balanced density anomaly using a linear
% equation of state for the error covariance balance operator:
%
%          deltaR_b = -rho0 * alpha * deltaT + rho0 * beta * deltaS_b
%
% where alpha is thermal expansion coefficient and beta is the saline
% contraction coefficient.
%
% On Input:
%
%    K           Balance operator structure array:
%
%                  K.rho0      Mean density (Kg/m3) used when the Boussinesq
%                                approximation is inferred
%                  K.alpha     Surface thermal expansion coefficient
%                                (1/Celsius)
%                  K.beta      Surface saline contraction coefficient
%                                (nondimensional)
%
%    deltaT      Temperature anomaly (T-Tavg) from a time mean (3d array)
%    deltaS_b    Balanced salinity anomaly (3D array), use the values
%                  computed from function "s_balance".
%
% On Output:
%
%    deltaR_b    Balanced density anomaly.
%

% svn $Id: rho_balance.m 996 2020-01-10 04:28:56Z arango $
%===========================================================================%
%  Copyright (c) 2002-2020 The ROMS/TOMS Group                              %
%    Licensed under a MIT/X style license                                   %
%    See License_ROMS.txt                           Hernan G. Arango        %
%===========================================================================%

[Lp,Mp,N]=size(deltaT);

deltaR_b=-K.rho0.*repmat(K.alpha,[1 1 N]).*deltaT +      ...
          K.rho0.*repmat(K.beta ,[1 1 N]).*deltaS_b;

return
