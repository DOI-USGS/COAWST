function cp = sw_cp(S,T,P)

% SW_CP      Heat Capacity (Cp) of sea water
%=========================================================================
% SW_CP  $Id: sw_cp.m 330 2009-03-10 05:57:42Z arango $
%         Copyright (C) CSIRO, Phil Morgan 1993.
%
% USAGE: cp = sw_cp(S,T,P)
%
% DESCRIPTION:
%    Heat Capacity of Sea Water using UNESCO 1983 polynomial.
%
% INPUT:  (all must have same dimensions)
%   S = salinity    [psu      (PSS-78)]
%   T = temperature [degree C (ITS-90)]
%   P = pressure    [db]
%       (P may have dims 1x1, mx1, 1xn or mxn for S(mxn) )
%
% OUTPUT:
%   cp = Specific Heat Capacity  [J kg^-1 C^-1]
%
% AUTHOR:  Phil Morgan, Lindsay Pender (Lindsay.Pender@csiro.au)
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
if nargin ~=3
   error('Must pass 3 parameters')
end %if

% CHECK S,T,P dimensions and verify consistent
[ms,ns] = size(S);
[mt,nt] = size(T);
[mp,np] = size(P);


% CHECK THAT S & T HAVE SAME SHAPE
if (ms~=mt) | (ns~=nt)
   error('S & T must have same dimensions')
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
   error('P has wrong dimensions')
end %if

%***check_stp

%------
% BEGIN
%------
P = P/10; % to convert db to Bar as used in Unesco routines
T68 = T * 1.00024;

%------------
% eqn 26 p.32
%------------
c0 = 4217.4;
c1 =   -3.720283;
c2 =    0.1412855;
c3 =   -2.654387e-3;
c4 =    2.093236e-5;

a0 = -7.64357;
a1 =  0.1072763;
a2 = -1.38385e-3;

b0 =  0.1770383;
b1 = -4.07718e-3;
b2 =  5.148e-5;

Cpst0 = (((c4.*T68 + c3).*T68 + c2).*T68 + c1).*T68 + c0 + ...
        (a0 + a1.*T68 + a2.*T68.^2).*S + ...
    (b0 + b1.*T68 + b2.*T68.^2).*S.*sqrt(S);

%------------
% eqn 28 p.33
%------------
a0 = -4.9592e-1;
a1 =  1.45747e-2;
a2 = -3.13885e-4;
a3 =  2.0357e-6;
a4 =  1.7168e-8;

b0 =  2.4931e-4;
b1 = -1.08645e-5;
b2 =  2.87533e-7;
b3 = -4.0027e-9;
b4 =  2.2956e-11;

c0 = -5.422e-8;
c1 =  2.6380e-9;
c2 = -6.5637e-11;
c3 =  6.136e-13;

del_Cp0t0 =  (((((c3.*T68 + c2).*T68 + c1).*T68 + c0).*P + ...
             ((((b4.*T68 + b3).*T68 + b2).*T68 + b1).*T68 + b0)).*P + ...
             ((((a4.*T68 + a3).*T68 + a2).*T68 + a1).*T68 + a0)).*P;

%------------
% eqn 29 p.34
%------------
d0 =  4.9247e-3;
d1 = -1.28315e-4;
d2 =  9.802e-7;
d3 =  2.5941e-8;
d4 = -2.9179e-10;

e0 = -1.2331e-4;
e1 = -1.517e-6;
e2 =  3.122e-8;

f0 = -2.9558e-6;
f1 =  1.17054e-7;
f2 = -2.3905e-9;
f3 =  1.8448e-11;

g0 =  9.971e-8;

h0 =  5.540e-10;
h1 = -1.7682e-11;
h2 =  3.513e-13;

j1 = -1.4300e-12;
S3_2  = S.*sqrt(S);

del_Cpstp = [((((d4.*T68 + d3).*T68 + d2).*T68 + d1).*T68 + d0).*S + ...
             ((e2.*T68 + e1).*T68 + e0).*S3_2].*P                + ...
        [(((f3.*T68 + f2).*T68 + f1).*T68 + f0).*S            + ...
         g0.*S3_2].*P.^2                                  + ...
         [((h2.*T68 + h1).*T68 + h0).*S                      + ...
         j1.*T68.*S3_2].*P.^3;


cp = Cpst0 + del_Cp0t0 + del_Cpstp;

return
