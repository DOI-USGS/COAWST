function [ALPHA] = sw_alpha(S, T, P, keyword)

% SW_ALPHA   Thermal expansion coefficient (alpha)
%================================================================
% SW_ALPHA  $Id: sw_alpha.m 330 2009-03-10 05:57:42Z arango $
%           Copyright (C) CSIRO, Nathan Bindoff 1993.
%
% USAGE:  [ALPHA] = alpha(S, T, P, keyword)
%
%         [ALPHA] = alpha(S, T,    P, 'temp')    %default
%         [ALPHA] = alpha(S, PTMP, P, 'ptmp')
%
% DESCRIPTION:
%    A function to calculate the thermal expansion coefficient.
%
% INPUT:
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
% OUTPUT:
%   ALPHA = Thermal expansion coeff (alpha) [degree_C.^-1]
%
% AUTHOR:   N.L. Bindoff  1993, Lindsay Pender (Lindsay.Pender@csiro.au)
%
% DISCLAIMER:
%   This software is provided "as is" without warranty of any kind.
%   See the file sw_copy.m for conditions of use and licence.
%
% REFERENCE:
%    McDougall, T.J. 1987.  "Neutral Surfaces"
%    Journal of Physical Oceanography vol 17 pages 1950-1964,
%
% CHECK VALUE:
%    See sw_beta.m amd sw_aonb.m
%

% Modifications
% 93-04-22. Phil Morgan,  Help display modified to suit library
% 93-04-23. Phil Morgan,  Input argument checking
% 94-10-15. Phil Morgan,  Pass S,T,P and keyword for 'ptmp'
% 99-06-25. Lindsay Pender, Fixed transpose of row vectors.
% 03-12-12. Lindsay Pender, Converted to ITS-90.

% CHECK INPUT ARGUMENTS
if ~(nargin==3 | nargin==4)
  error('sw_alpha.m: requires 3 or 4 input arguments')
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

% BEGIN

ALPHA = sw_aonb(S,T,P,keyword).*sw_beta(S,T,P,keyword);

return
