
function fp = sw_fp(S,P)

% SW_FP      Freezing point of sea water
%=========================================================================
% SW_FP % $Id: sw_fp.m 330 2009-03-10 05:57:42Z arango $
%         Copyright (C) CSIRO, Phil Morgan 1993.
%
% USAGE:  fp = sw_fp(S,P)
%
% DESCRIPTION:
%    Freezing point of Sea Water using UNESCO 1983 polynomial.
%
% INPUT:  (all must have same dimensions)
%   S = salinity    [psu      (PSS-78)]
%   P = pressure    [db]
%       (P may have dims 1x1, mx1, 1xn or mxn for S(mxn) )
%
% OUTPUT:
%   fp = Freezing Point temperature [degree C (ITS-90)]
%
% AUTHOR:  Phil Morgan 93-04-20, Lindsay Pender (Lindsay.Pender@csiro.au)
%
% DISCLAIMER:
%   This software is provided "as is" without warranty of any kind.
%   See the file sw_copy.m for conditions of use and licence.
%
% REFERENCES:
%    Fofonff, P. and Millard, R.C. Jr
%    Unesco 1983. Algorithms for computation of fundamental properties of
%    seawater, 1983. _Unesco Tech. Pap. in Mar. Sci._, No. 44, 53 pp.
%

% Modifications
% 99-06-25. Lindsay Pender, Fixed transpose of row vectors.
% 03-12-12. Lindsay Pender, Converted to ITS-90.

% CALLER: general purpose
% CALLEE: none

%----------------------
% CHECK INPUT ARGUMENTS
%----------------------
if nargin ~=2
   error('Must pass 2 parameters')
end %if

[ms,ns] = size(S);
[mp,np] = size(P);

% CHECK OPTIONAL SHAPES FOR P
if     mp==1  & np==1      % P is a scalar.  Fill to size of S
   P = P(1)*ones(ms,ns);
elseif np==ns & mp==1      % P is row vector with same cols as S
   P = P( ones(1,ms), : ); %   Copy down each column.
elseif mp==ms & np==1      % P is column vector
   P = P( :, ones(1,ns) ); %   Copy across each row
elseif mp==ms & np==ns     % PR is a matrix size(S)
   % shape ok
else
  error('P has wrong dimensions')
end %if

%------
% BEGIN
%------
%P = P/10; % to convert db to Bar as used in Unesco routines

%------------
% eqn  p.29
%------------
a0 = -0.0575;
a1 = 1.710523e-3;
a2 = -2.154996e-4;
b  = -7.53e-4;

fp = (a0.*S + a1.*S.*sqrt(S) + a2.*S.^2 + b.*P) / 1.00024;

return
