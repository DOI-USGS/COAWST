function [ga] = sw_gpan(S,T,P)

% SW_GPAN    Geopotential anomaly
%=========================================================================
% SW_GPAN  $Id: sw_gpan.m 330 2009-03-10 05:57:42Z arango $
%          Copyright (C) CSIRO, Phil Morgan 1992.
%
% USAGE:  [gpan]= sw_gpan(S,T,P)
%
% DESCRIPTION:
%   Geopotential Anomaly calculated as the integral of svan from the
%   the sea surface to the bottom.  Thus RELATIVE TO SEA SURFACE.
%
% INPUT:  (all must have same dimensions)
%   S = salinity    [psu      (PSS-78)]
%   T = temperature [degree C (ITS-90)]
%   P = Pressure    [db]
%       (P may have dims 1x1, mx1, 1xn or mxn for S(mxn) )
%
% OUTPUT:
%  gpan = Geopotential Anomaly  [m^3 kg^-1 Pa == m^2 s^-2 == J kg^-1]
%
% AUTHOR:  Phil Morgan 92-11-05, Lindsay Pender (Lindsay.Pender@csiro.au)
%
% DISCLAIMER:
%   This software is provided "as is" without warranty of any kind.
%   See the file sw_copy.m for conditions of use and licence.
%
% REFERENCE: S. Pond & G.Pickard  2nd Edition 1986
%            Introductory Dynamical Oceanogrpahy
%            Pergamon Press Sydney.  ISBN 0-08-028728-X
%
% Note that older literature may use units of "dynamic decimeter' for above.
%
% Adapted method from Pond and Pickard (p76) to calc gpan rel to sea
% surface whereas P&P calculated relative to the deepest common depth.
%

% Modifications
% 03-12-12. Lindsay Pender, Converted to ITS-90.

%
% CALLER: general purpose
% CALLEE: sw_svan.m

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
[mp,np] = size(P);



% IF ALL ROW VECTORS ARE PASSED THEN LET US PRESERVE SHAPE ON RETURN.
Transpose = 0;
if mp == 1  % row vector
   P       =  P(:);
   T       =  T(:);
   S       =  S(:);

   Transpose = 1;
end %if
%***check_stp

%------
% BEGIN
%------
db2Pascal  = 1e4;
[m,n]      = size(P);
svan       = sw_svan(S,T,P);
mean_svan  = 0.5*(svan(2:m,:) + svan(1:m-1,:) );

if n==1
   top = svan(1,1).*P(1,1)*db2Pascal;
else
   top = svan(1,:).*P(1,:)*db2Pascal;
end %if

%press_diff = diff(P);

delta_ga   = (mean_svan.*diff(P))*db2Pascal;
ga         = cumsum([ top; delta_ga]);

if Transpose
   ga = ga';
end %if

return
