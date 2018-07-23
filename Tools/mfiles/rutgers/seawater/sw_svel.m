function svel = sw_svel(S,T,P)

% SW_SVEL    Sound velocity of sea water
%=========================================================================
% SW_SVEL  $Id: sw_svel.m 330 2009-03-10 05:57:42Z arango $
%          Copyright (C) CSIRO, Phil Morgan 1993.
%
% USAGE:  svel = sw_svel(S,T,P)
%
% DESCRIPTION:
%    Sound Velocity in sea water using UNESCO 1983 polynomial.
%
% INPUT:  (all must have same dimensions)
%   S = salinity    [psu      (PSS-78)]
%   T = temperature [degree C (ITS-90)]
%   P = pressure    [db]
%       (P may have dims 1x1, mx1, 1xn or mxn for S(mxn) )
%
% OUTPUT:
%   svel = sound velocity  [m/s]
%
% AUTHOR:  Phil Morgan 93-04-20, Lindsay Pender (Lindsay.Pender@csiro.au)
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
% 99-06-25. Lindsay Pender, Fixed transpose of row vectors.
% 03-12-12. Lindsay Pender, Converted to ITS-90.

% CALLER: general purpose
% CALLEE: none

% UNESCO 1983. eqn.33  p.46

%----------------------
% CHECK INPUT ARGUMENTS
%----------------------
if nargin ~=3
   error('sw_svel.m: Must pass 3 parameters')
end %if

% CHECK S,T,P dimensions and verify consistent
[ms,ns] = size(S);
[mt,nt] = size(T);
[mp,np] = size(P);


% CHECK THAT S & T HAVE SAME SHAPE
if (ms~=mt) | (ns~=nt)
   error('check_stp: S & T must have same dimensions')
end %if

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
   error('check_stp: P has wrong dimensions')
end %if

%***check_stp

%---------
% BEGIN
%--------

P = P/10;  % convert db to bars as used in UNESCO routines
T68 = T * 1.00024;

%------------
% eqn 34 p.46
%------------
c00 = 1402.388;
c01 =    5.03711;
c02 =   -5.80852e-2;
c03 =    3.3420e-4;
c04 =   -1.47800e-6;
c05 =    3.1464e-9;

c10 =  0.153563;
c11 =  6.8982e-4;
c12 = -8.1788e-6;
c13 =  1.3621e-7;
c14 = -6.1185e-10;

c20 =  3.1260e-5;
c21 = -1.7107e-6;
c22 =  2.5974e-8;
c23 = -2.5335e-10;
c24 =  1.0405e-12;

c30 = -9.7729e-9;
c31 =  3.8504e-10;
c32 = -2.3643e-12;

Cw = ((((c32.*T68 + c31).*T68 + c30).*P + ...
       ((((c24.*T68 + c23).*T68 + c22).*T68 + c21).*T68 + c20)).*P + ...
       ((((c14.*T68 + c13).*T68 + c12).*T68 + c11).*T68 + c10)).*P + ...
       ((((c05.*T68 + c04).*T68 + c03).*T68 + c02).*T68 + c01).*T68 + c00;

%-------------
% eqn 35. p.47
%-------------
a00 =  1.389;
a01 = -1.262e-2;
a02 =  7.164e-5;
a03 =  2.006e-6;
a04 = -3.21e-8;

a10 =  9.4742e-5;
a11 = -1.2580e-5;
a12 = -6.4885e-8;
a13 =  1.0507e-8;
a14 = -2.0122e-10;

a20 = -3.9064e-7;
a21 =  9.1041e-9;
a22 = -1.6002e-10;
a23 =  7.988e-12;

a30 =  1.100e-10;
a31 =  6.649e-12;
a32 = -3.389e-13;

A = ((((a32.*T68 + a31).*T68 + a30).*P + ...
    (((a23.*T68 + a22).*T68 + a21).*T68 + a20)).*P + ...
    ((((a14.*T68 + a13).*T68 + a12).*T68 + a11).*T68 + a10)).*P + ...
    (((a04.*T68 + a03).*T68 + a02).*T68 + a01).*T68 + a00;

%------------
% eqn 36 p.47
%------------
b00 = -1.922e-2;
b01 = -4.42e-5;
b10 =  7.3637e-5;
b11 =  1.7945e-7;

B = b00 + b01.*T68 + (b10 + b11.*T68).*P;

%------------
% eqn 37 p.47
%------------
d00 =  1.727e-3;
d10 = -7.9836e-6;

D = d00 + d10.*P;

%------------
% eqn 33 p.46
%------------
svel = Cw + A.*S + B.*S.*sqrt(S) + D.*S.^2;

return

