function [BETA] = sw_beta(S, T, P, keyword)

% SW_BETA    Saline contraction coefficient (beta)
%========================================================================
% SW_BETA  $Id: sw_beta.m 330 2009-03-10 05:57:42Z arango $
%   %      Copyright (C) CSIRO, Nathan Bindoff 1993.
%
% USAGE:  [BETA] = sw_beta(S, T, P, {keyword} )
%
%         [BETA] = sw_beta(S, T,    P, 'temp')     %default
%         [BETA] = sw_beta(S, PTMP, P, 'ptmp')
%
% DESCRIPTION
%   The saline contraction coefficient as defined by T.J. McDougall.
%
% INPUT:  (all must have same dimensions)
%   S       = salinity              [psu      (PSS-78) ]
% * PTMP    = potential temperature [degree C (ITS-90)]
% * T       = temperature           [degree C (ITS-90)]
%   P       = pressure              [db]
%             (P may have dims 1x1, mx1, 1xn or mxn for S(mxn) )
%
%   keyword = optional string to identify if temp or ptmp passed.
%           = No argument defaults to 'temp'
%           = 'temp' assumes (S,T,P) passed.    Will execute slower
%                    as ptmp will be calculated internally.
%           = 'ptmp' assumes (S,PTMP,P) passed. Will execute faster.
%
% OUTPUT
%   BETA = Saline Contraction Coefficient  [psu.^-1]
%
% AUTHOR:   N.L. Bindoff  1993, Lindsay Pender (Lindsay.pender@csiro.au)
%
% DISCLAIMER:
%   This software is provided "as is" without warranty of any kind.
%   See the file sw_copy.m for conditions of use and licence.
%
% REFERENCE:
%    McDougall, T.J. 1987. "Neutral Surfaces"
%    Journal of Physical Oceanography vol 17 pages 1950-1964,
%
% CHECK VALUE:
%    beta=0.72088e-3 psu.^-1 at S=40.0 psu, ptmp = 10.0 C (ITS-68), p=4000 db
%

% Modifications
% 93-04-22. Phil Morgan,  Help display modified to suit library
% 93-04-23. Phil Morgan,  Input argument checking
% 94-10-15. Phil Morgan,  Pass S,T,P and keyword for 'ptmp'
% 99-06-25. Lindsay Pender, Fixed transpose of row vectors.
% 03-12-12. Lindsay Pender, Converted to ITS-90.

% CHECK INPUT ARGUMENTS
if ~(nargin==3 | nargin==4)
  error('sw_beta.m: requires 3 or 4 input arguments')
end %if
if nargin == 3
  keyword = 'temp';
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

% ENSURE WE USE PTMP IN CALCULATIONS
if ~strcmp(lower(keyword),'ptmp')
  T = sw_ptmp(S,T,P,0); % now have ptmp
end %if

T = T * 1.00024;

% BEGIN

        c1=fliplr([ 0.785567e-3, -0.301985e-5 ...
         0.555579e-7, -0.415613e-9]);
    c2=fliplr([ -0.356603e-6, 0.788212e-8]);
    c3=fliplr([0.0 0.408195e-10, -0.602281e-15]);
    c4=[0.515032e-8];
    c5=fliplr([-0.121555e-7, 0.192867e-9, -0.213127e-11]);
        c6=fliplr([0.176621e-12 -0.175379e-14]);
    c7=[0.121551e-17];
%
% Now calaculate the thermal expansion saline contraction ratio adb
%
    [m,n] = size(S);
    sm35  = S-35*ones(m,n);
    BETA  = polyval(c1,T) + sm35.*(polyval(c2,T) + ...
            polyval(c3,P)) + c4*(sm35.^2) + ...
            P.*polyval(c5,T) + (P.^2).*polyval(c6,T) ...
                +c7*( P.^3);

return
