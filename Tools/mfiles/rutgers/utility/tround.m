function [Y]=tround(X, CT)

% TROUND: Rounding function with a Fuzzy or Tolerant Floor function
%
% [x]=tround(X, CT)
%
% Floating point rounding function with a Fuzzy or Tolerant Floor
% function.
%  
% On Input:
%
%    X         Double precision argument to be operated on. It is
%              assumed that X is represented with m mantissa bits.
%
%    CT        Comparison Tolerance such that 0 < CT <= 3-SQRT(5)/2.
%              If the relative difference between X and a whole
%              number is less than CT, then TFLOOR is returned as
%              this whole number. By treating the floating-point
%              numbers as a finite ordered set, note that the
%              heuristic EPS=2.**(-(m-1)) and CT=3*eps causes
%              arguments of TFLOOR/TCEIL to be treated as whole
%              numbers if they are exactly whole numbers or are
%              immediately adjacent to whole number representations.
%              Since EPS, the  "distance"  between  floating-point
%              numbers on the unit interval, and m, the number of
%              bits in X mantissa, exist on every  floating-point
%              computer, TFLOOR/TCEIL are consistently definable
%              on every floating-point computer.
%
% On Input:
%
%    Y         Double precision round of X.
%
% Usage:
%
%    CT = 3 * eps(X)         That is, CT is about 1 bit on either
%                            side of X mantissa bits.
%     Y = round(X, CT)
%
% References:
%
%    P. E. Hagerty, 1978: More on Fuzzy Floor and Ceiling, APL QUOTE
%      QUAD 8(4):20-24. (The TFLOOR=FL5 took five years of refereed
%      evolution publication).
%
%    L. M. Breed, 1978: Definitions for Fuzzy Floor and Ceiling, APL
%      QUOTE QUAD 8(3):16-23.
%
% Adapted from H.D. Knoble code (Penn State University).

% svn $Id: tround.m 996 2020-01-10 04:28:56Z arango $
%=========================================================================%
%  Copyright (c) 2002-2020 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%

Y=tfloor(X+0.5d0, CT);

return

%--------------------------------------------------------------------------

function [Y]=tceil(X, CT)

% Double precision Tolerant ceiling function.

Y=-tfloor(-X, CT);

return

%--------------------------------------------------------------------------

function [Y]=tfloor(X, CT)

% Double precision Tolerant floor function: Hagerty FL5 function.

Q=1.0d0;
if (X < 0.0d0)
  Q=1.0d0-CT;
end
RMAX=Q./(2.0d0-CT);
EPS5=CT./Q;
Y=ufloor(X+max(CT, min(RMAX, EPS5.*abs(1.0d0+ufloor(X)))));
if ((X <= 0.0d0) | (Y-X) < RMAX)
else  
  Y=Y-1.0d0;
end

return

%--------------------------------------------------------------------------

function [Y]=ufloor(X, CT)

% Compute the largest integer algebraically less than or equal to X;
% that is, the unfuzzy Floor Function.

Y=X-mod(X,1.0d0)-mod(2.0d0+sign(X),3.0d0);

return
