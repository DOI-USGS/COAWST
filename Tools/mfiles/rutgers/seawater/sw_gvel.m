function vel = sw_gvel(ga,lat,lon)

% SW_GVEL    Geostrophic velocity
%===================================================================
% GEOVEL   $Id: sw_gvel.m 330 2009-03-10 05:57:42Z arango $
%          Copyright (C) CSIRO, Phil Morgan 1992
%
% USAGE:  vel = sw_gvel(ga,lat,lon)
%
% DESCRIPTION:
%    Calculates geostrophic velocity given the geopotential anomaly
%    and position of each station.
%
% INPUT:
%    ga   = geopotential anomoly relative to the sea surface.
%           dim(mxnstations)
%    lat  = latitude  of each station (+ve = N, -ve = S) [ -90.. +90]
%    lon  = longitude of each station (+ve = E, -ve = W) [-180..+180]
%
% OUTPUT:
%    vel  = geostrophic velocity RELATIVE to the sea surface.
%           dim(m,nstations-1)
%
% AUTHOR:   Phil Morgan   1992/03/26  (morgan@ml.csiro.au)
%
% DISCLAIMER:
%   This software is provided "as is" without warranty of any kind.
%   See the file sw_copy.m for conditions of use and licence.
%
% REFERENCE: S. Pond & G.Pickard  2nd Edition 1986
%            Introductory Dynamical Oceanogrpahy
%            Pergamon Press Sydney.  ISBN 0-08-028728-X
%            Equation 8.9A p73  Pond & Pickard
%
% NOTE: This calls sw_dist.m.  You can replace the call to this
%       routine if you have a more appropraite distance routine.
%

% CALLER:   general purpose
% CALLEE:   sw_dist.m
%


DEG2RAD = pi/180;
RAD2DEG = 180/pi;
OMEGA   = 7.292e-5;  % Angular velocity of Earth  [radians/sec]

% You may replace the call to sw_dist if you have
% a more appropriate distance routine.
distm = 1000*sw_dist(lat,lon,'km');
[m,n] = size(ga);
f     = 2*OMEGA*sin( (lat(1:n-1)+lat(2:n))*DEG2RAD/2 );
lf    = f.*distm;

LF = lf(ones(m,1),:);

vel   = -( ga(:,2:n)-ga(:,1:n-1) ) ./ LF;

return
