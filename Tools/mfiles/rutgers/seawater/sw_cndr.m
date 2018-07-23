function R = sw_cndr(S,T,P)

% SW_CNDR    Conductivity ratio   R = C(S,T,P)/C(35,15(IPTS-68),0)
%=========================================================================
% SW_CNDR  $Id: sw_cndr.m 330 2009-03-10 05:57:42Z arango $
%          Copyright (C) CSIRO, Phil Morgan 1993.
%
% USAGE:  cndr = sw_cndr(S,T,P)
%
% DESCRIPTION:
%   Calculates conductivity ratio from S,T,P.
%
% INPUT:  (all must have same dimensions)
%   S = salinity    [psu      (PSS-78) ]
%   T = temperature [degree C (ITS-90)]
%   P = pressure    [db]
%       (P may have dims 1x1, mx1, 1xn or mxn for S(mxn) )
%
% OUTPUT:
%   cndr = Conductivity ratio     R =  C(S,T,P)/C(35,15(IPTS-68),0) [no units]
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

% Modifications
% 99-06-25. Lindsay Pender, Fixed transpose of row vectors.
% 03-12-12. Lindsay Pender, Converted to ITS-90.

% CALLER: general purpose
% CALLEE: sw_salds.m sw_sals.m sw_salrt.m

%--------------
% check inputs
%-------------
if nargin~=3
  error('sw_cndr.m: must have 3 input arguments')
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

%-------
% BEGIN
%-------

T68 = T * 1.00024;

for i = 1:ms
  for j = 1:ns
    %---------------------------------------------------------------------
    % DO A NEWTON-RAPHSON ITERATION FOR INVERSE INTERPOLATION OF Rt FROM S.
    %---------------------------------------------------------------------
    S_loop  = S(i,j);  % S in the loop
    T_loop  = T(i,j);  % T in the loop
    Rx_loop = sqrt(S_loop/35.0);                % first guess at Rx = sqrt(Rt)
    SInc    = sw_sals(Rx_loop.*Rx_loop,T_loop); % S INCrement (guess) from Rx
    iloop    = 0;
    end_loop = 0;
    while ~end_loop
       Rx_loop = Rx_loop + (S_loop - SInc)./sw_salds(Rx_loop,T_loop - 15);
       SInc    = sw_sals(Rx_loop.*Rx_loop,T_loop);
       iloop   = iloop + 1;
       dels    = abs(SInc-S_loop);
       if (dels>1.0e-4 & iloop<10)
          end_loop = 0;
       else
          end_loop = 1;
       end %if
    end %while

    Rx(i,j) = Rx_loop;

  end %for j
end %for i

%------------------------------------------------------
% ONCE Rt FOUND, CORRESPONDING TO EACH (S,T) EVALUATE R
%------------------------------------------------------
% eqn(4) p.8 Unesco 1983

d1 =  3.426e-2;
d2 =  4.464e-4;
d3 =  4.215e-1;
d4 = -3.107e-3;

e1 =  2.070e-5;
e2 = -6.370e-10;
e3 =  3.989e-15;

A  = (d3 + d4.*T68);
B  = 1 + d1.*T68 + d2.*T68.^2;
C  = P.*(e1 + e2.*P + e3.*P.^2);

% eqn(6) p.9 UNESCO 1983.
Rt    = Rx.*Rx;
rt    = sw_salrt(T);
Rtrt  = rt.*Rt;
D     = B - A.*rt.*Rt;
E     = rt.*Rt.*A.*(B+C);
R     = sqrt(abs(D.^2+4*E)) - D;
R     = 0.5*R./A;

return
