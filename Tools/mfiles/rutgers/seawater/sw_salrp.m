
function Rp = sw_salrp(R,T,P)

% SW_SALRP   Conductivity ratio   Rp(S,T,P) = C(S,T,P)/C(S,T,0)
%=========================================================================
% SW_SALRP   $Id: sw_salrp.m 330 2009-03-10 05:57:42Z arango $
%            Copyright (C) CSIRO, Phil Morgan 1993.
%
% USAGE:  Rp = sw_salrp(R,T,P)
%
% DESCRIPTION:
%    Equation Rp(S,T,P) = C(S,T,P)/C(S,T,0) used in calculating salinity.
%    UNESCO 1983 polynomial.
%
% INPUT: (All must have same shape)
%   R = Conductivity ratio  R =  C(S,T,P)/C(35,15(IPTS-68),0) [no units]
%   T = temperature [degree C (ITS-90)]
%   P = pressure    [db]
%
% OUTPUT:
%   Rp = conductivity ratio  Rp(S,T,P) = C(S,T,P)/C(S,T,0)  [no units]
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

%-------------------
% CHECK INPUTS
%-------------------
if nargin~=3
  error('sw_salrp.m: requires 3 input arguments')
end %if

[mr,nr] = size(R);
[mp,np] = size(P);
[mt,nt] = size(T);
if ~(mr==mp | mr==mt | nr==np | nr==nt)
   error('sw_salrp.m: R,T,P must all have the same shape')
end %if

%-------------------
% eqn (4) p.8 unesco.
%-------------------

T68 = T * 1.00024;

d1 =  3.426e-2;
d2 =  4.464e-4;
d3 =  4.215e-1;
d4 = -3.107e-3;

e1 =  2.070e-5;
e2 = -6.370e-10;
e3 =  3.989e-15;

Rp = 1 + ( P.*(e1 + e2.*P + e3.*P.^2) ) ...
     ./ (1 + d1.*T68 + d2.*T68.^2 +(d3 + d4.*T68).*R);

return

