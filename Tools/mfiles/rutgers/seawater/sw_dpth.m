function DEPTHM = sw_dpth(P,LAT)

% SW_DPTH    Depth from pressure
%===========================================================================
% SW_DPTH   $Id: sw_dpth.m 330 2009-03-10 05:57:42Z arango $
%           Copyright (C) CSIRO, Phil Morgan 1992.
%
% USAGE:  dpth = sw_dpth(P,LAT)
%
% DESCRIPTION:
%    Calculates depth in metres from pressure in dbars.
%
% INPUT:  (all must have same dimensions)
%   P   = Pressure    [db]
%   LAT = Latitude in decimal degress north [-90..+90]
%         (lat may have dimensions 1x1 or 1xn where P(mxn).
%
% OUTPUT:
%  dpth = depth [metres]
%
% AUTHOR:  Phil Morgan 92-04-06  (morgan@ml.csiro.au)
%
% DISCLAIMER:
%   This software is provided "as is" without warranty of any kind.
%   See the file sw_copy.m for conditions of use and licence.
%
% REFERENCES:
%    Unesco 1983. Algorithms for computation of fundamental properties of
%    seawater, 1983. _Unesco Tech. Pap. in Mar. Sci._, No. 44, 53 pp.

% Modifications
% 99-06-25. Lindsay Pender, Fixed transpose of row vectors.

% CALLER:  general purpose
% CALLEE:  none

%-------------
% CHECK INPUTS
%-------------
[mP,nP] = size(P);
[mL,nL] = size(LAT);
if mL==1 & nL==1                    % LAT scalar - fill to size of P
  LAT = LAT*ones(size(P));

elseif nP == nL & mL == 1           % LAT is row vector
  LAT = LAT(ones(1, mP), :);        % Coppy down each column

elseif mP == mL & nL == 1           % LAT is column vector
  LAT = LAT(:, ones(1, nP));        % Copy across each row

elseif mP == mL & nP == nL
  % Ok

else
   error('sw_depth.m:  Inputs arguments have wrong dimensions')
end %if

%-------------
% BEGIN
%-------------
% Eqn 25, p26.  Unesco 1983.

DEG2RAD = pi/180;
c1 = +9.72659;
c2 = -2.2512E-5;
c3 = +2.279E-10;
c4 = -1.82E-15;
gam_dash = 2.184e-6;

LAT = abs(LAT);
X   = sin(LAT*DEG2RAD);  % convert to radians
X   = X.*X;
bot_line = 9.780318*(1.0+(5.2788E-3+2.36E-5*X).*X) + gam_dash*0.5*P;
top_line = (((c4*P+c3).*P+c2).*P+c1).*P;
DEPTHM   = top_line./bot_line;
return
