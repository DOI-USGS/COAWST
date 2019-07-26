function S = sw_sals(Rt,T)

% SW_SALS    Salinity of sea water
%=========================================================================
% SW_SALS  $Id: sw_sals.m 330 2009-03-10 05:57:42Z arango $
%          Copyright (C) CSIRO, Phil Morgan 1993.
%
% USAGE:  S = sw_sals(Rt,T)
%
% DESCRIPTION:
%    Salinity of sea water as a function of Rt and T.
%    UNESCO 1983 polynomial.
%
% INPUT:
%   Rt = Rt(S,T) = C(S,T,0)/C(35,T(IPTS-68),0)
%   T  = temperature [degree C (ITS-90)]
%
% OUTPUT:
%   S  = salinity    [psu      (PSS-78)]
%
% AUTHOR:  Phil Morgan 93-04-17, Lindsay Pender (Lindsay.Pender@csiro.au)
%
% DISCLAIMER:
%   This software is provided "as is" without warranty of any kind.
%   See the file sw_copy.m for conditions of use and licence.
%
% REFERENCES:
%    Fofonoff, P. and Millard, R.C. Jr
%    Unesco 1983. Algorithms for computation of fundamental properties of
%    seawater, 1983. _Unesco Tech. Pap. in Mar. Sci._, No. 44, 53 pp.
%

% Modifications
% 03-12-12. Lindsay Pender, Converted to ITS-90.

% CALLER: sw_salt
% CALLEE: none

%--------------------------
% CHECK INPUTS
%--------------------------
if nargin~=2
  error('sw_sals.m: requires 2 input arguments')
end %if

[mrt,nrt] = size(Rt);
[mT,nT]   = size(T);
if ~(mrt==mT | nrt==nT)
   error('sw_sals.m: Rt and T must have the same shape')
end %if

%--------------------------
% eqn (1) & (2) p6,7 unesco
%--------------------------

del_T68 = T * 1.00024 - 15;

a0 =  0.0080;
a1 = -0.1692;
a2 = 25.3851;
a3 = 14.0941;
a4 = -7.0261;
a5 =  2.7081;

b0 =  0.0005;
b1 = -0.0056;
b2 = -0.0066;
b3 = -0.0375;
b4 =  0.0636;
b5 = -0.0144;

k  =  0.0162;

Rtx   = sqrt(Rt);
del_S = (del_T68 ./ (1+k*del_T68) ) .* ...
        ( b0 + (b1 + (b2+ (b3 + (b4 + b5.*Rtx).*Rtx).*Rtx).*Rtx).*Rtx);

S = a0 + (a1 + (a2 + (a3 + (a4 + a5.*Rtx).*Rtx).*Rtx).*Rtx).*Rtx;

S = S + del_S;

return
