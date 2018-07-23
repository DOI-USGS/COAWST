function speed = sw_swvel(len, depth)

% SW_SWVEL    Surface wave velocity
%===========================================================================
% SW_SWVEL   $Id: sw_swvel.m 330 2009-03-10 05:57:42Z arango $
%           Copyright (C) CSIRO, Phil Morgan 1993.
%
% USAGE:  speed = sw_swvel(len, depth)
%
% DESCRIPTION:
%    Calculates surface wave velocity.
%
% INPUT:  (all must have same dimensions)
%   len   = wave length
%   depth = water depth [metres]
%
% OUTPUT:
%  speed   = Surface wave speed (m/s)
%
% AUTHOR:  Lindsay Pender 2005
%
% DISCLAIMER:
%   This software is provided "as is" without warranty of any kind.
%   See the file sw_copy.m for conditions of use and licence.
%

% CALLER:  general purpose
% CALLEE:  none

%-------------
% CHECK INPUTS
%-------------
[mD,nD] = size(depth);
[mL,nL] = size(len);
if mD==1 & nD==1                    % depth scalar - fill to size of len
  depth = depth*ones(size(len));

elseif nL == nD & mD == 1           % depth is row vector
  depth = depth(ones(1, mL), :);    % Copy down each column

elseif mL == mD & nD == 1           % depth is column vector
  depth = depth(:, ones(1, nL));    % Copy across each row

elseif mD == mL & nD == nL
  % Ok

else
   error('sw_swvel.m:  Inputs arguments have wrong dimensions')
end %if

%-------------
% BEGIN
%-------------

g = 9.8;
k = 2.0 * pi ./ len;
speed = sqrt(g * tanh(k .* depth) ./ k);
return

