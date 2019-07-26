function f = sw_f(lat)

% SW_F       Coriolis factor "f"
%===========================================================================
% SW_F   $Id: sw_f.m 330 2009-03-10 05:57:42Z arango $
%        Copyright (C) CSIRO, Phil Morgan 1993.
%
% USAGE:  f = sw_f(lat)
%
% DESCRIPTION:
%    Calculates the Coriolis factor "f" defined by
%       f = 2*Omega*Sin(lat)  where Omega = 7.292e-5 radians/sec
%
% INPUT:
%   lat = Latitude in decimal degress north [-90..+90]
%
% OUTPUT:
%  f    = Coriolis Factor "f" [s-1]
%
% AUTHOR:  Phil Morgan 93-04-20  (morgan@ml.csiro.au)
%
% DISCLAIMER:
%   This software is provided "as is" without warranty of any kind.
%   See the file sw_copy.m for conditions of use and licence.
%
% REFERENCE:
%   S. Pond & G.Pickard  2nd Edition 1986
%   Introductory Dynamical Oceanogrpahy
%   Pergamon Press Sydney.  ISBN 0-08-028728-X
%
%   A.E. Gill 1982. p.597
%   "Atmosphere-Ocean Dynamics"
%   Academic Press: New York.  ISBN: 0-12-283522-0

% CALLER:  general purpose
% CALLEE:  none

%-------------
% CHECK INPUTS
%-------------
if nargin ~= 1
   error('sw_f.m:  Requires one input argument')
end %if

%-------------
% BEGIN
%-------------
% Eqn p27.  Unesco 1983.
DEG2RAD = pi/180;
OMEGA   = 7.292e-5;     %s-1   A.E.Gill p.597
f       = 2*OMEGA*sin(lat*DEG2RAD);

return


