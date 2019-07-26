function dS = sw_salds(Rtx,delT)

% SW_SALDS   Differiential dS/d(sqrt(Rt)) at constant T.
%=========================================================================
% SW_SALDS   $Id: sw_salds.m 330 2009-03-10 05:57:42Z arango $
%            Copyright (C) CSIRO, Phil Morgan 1993.
%
% USAGE:  dS = sw_salds(Rtx,delT)
%
% DESCRIPTION:
%   Calculates Salinity differential dS/d(sqrt(Rt)) at constant T.
%   UNESCO 1983 polynomial.
%
% INPUT: (all must have same dimensions)
%   Rtx   = sqrt(Rt) where Rt defined in sw_salt.m
%   delT  = T-15     [degree C (IPTS-68)]
%
% OUTPUT:
%   dS = S differential dS/d(sqrt(Rt)) at constant T.
%
% AUTHOR:  Phil Morgan 93-04-21, Lindsay Pender (Lindsay.Pender@csiro.au)
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

% CALLER: sw_cndr.m
% CALLEE: none

%-------------
% CHECK INPUTS
%-------------
if nargin~=2
   error('sw_salds.m: must have 2 input arguments')
end %if

[m1,n1] = size(Rtx);
[m2,n2] = size(delT);
if ~(m1==m2 | n1==n2)
  error('sw_salds.m: Rtx and delT must have the same shape')
end %if

%-------
% BEGIN
%-------

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

dS =  a1 + (2*a2 + (3*a3 + (4*a4 + 5*a5.*Rtx).*Rtx).*Rtx).*Rtx + ...
     (delT./(1+k*delT))* ...
     (b1 + (2*b2 + (3*b3 + (4*b4 + 5*b5.*Rtx).*Rtx).*Rtx).*Rtx);

return
