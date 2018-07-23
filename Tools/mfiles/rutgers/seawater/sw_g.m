function g = sw_g(LAT,z)

% SW_G       Gravitational acceleration
%===========================================================================
% SW_G   $Id: sw_g.m 330 2009-03-10 05:57:42Z arango $
%        Copyright (C) CSIRO, Phil Morgan 1993.
%
% USAGE:  g = sw_g(lat,z)
%
% DESCRIPTION:
%    Calculates acceleration due to gravity as function of latitude.
%
% INPUT:  (all must have same dimensions)
%   lat = Latitude in decimal degress north [-90..+90]
%   z   = height in metres (+ve above sea surface, -ve below)
%
% OUTPUT:
%  g    = gravity [m/s^2]
%
% AUTHOR:  Phil Morgan 93-04-20  (morgan@ml.csiro.au)
%
% DISCLAIMER:
%   This software is provided "as is" without warranty of any kind.
%   See the file sw_copy.m for conditions of use and licence.
%
% REFERENCES:
%   Unesco 1983. Algorithms for computation of fundamental properties of
%   seawater, 1983. _Unesco Tech. Pap. in Mar. Sci._, No. 44, 53 pp.
%
%   A.E. Gill 1982. p.597
%   "Atmosphere-Ocean Dynamics"
%   Academic Press: New York.  ISBN: 0-12-283522-0
%

% CALLER:  general purpose
% CALLEE:  none

%-------------
% CHECK INPUTS
%-------------
if ~(nargin==1 | nargin==2)
   error('sw_g.m:  Requires one or two input arguments')
end %if
if nargin == 1
  z = zeros(size(LAT));
end %if

[mL,nL] = size(LAT);
[mz,nz] = size(z);
if ~(mL==mz | nL==nz)
   error('sw_g.m:  Input arguments should have same dimensions')
end %if

%-------------
% BEGIN
%-------------
% Eqn p27.  Unesco 1983.
a       = 6371000;    % mean radius of earth  A.E.Gill
DEG2RAD = pi/180;
LAT     = abs(LAT);
X       = sin(LAT*DEG2RAD);  % convert to radians
sin2    = X.*X;
g       = 9.780318*(1.0+(5.2788E-3+2.36E-5*sin2).*sin2);
if any(any(z))
   g    = g./((1+z/a).^2);    % from A.E.Gill p.597
end %if
return
