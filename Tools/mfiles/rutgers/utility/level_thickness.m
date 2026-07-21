function [Hz]=level_thickness(fname, gname, tindex);

%
% LEVEL_THICKNESS:  Compute the ROMS level thickness
%
% [Hz]=level_thickness(fname, gname, tindex)
%
% This function computes the level thickness (m) at RHO-points. If
% the time record (tindex) is not provided, a zero free-surface is
% assumed and the unperturbed depths are returned.
%
% On Input:
%
%    fname       NetCDF data file name (character string)
%    gname       NetCDF grid file name (character string)
%    tindex      Time index (integer)
%
% On Output:
%
%    Hz          Level thickness (3D array; meters)
%

% svn $Id$
%=========================================================================%
%  Copyright (c) 2002-2025 The ROMS Group                                 %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.md                            Hernan G. Arango      %
%=========================================================================%

% Check arguments and set default values.

if (nargin < 3)
  tindex=0;
end

%--------------------------------------------------------------------------
% Compute ROMS level thickness at RHO-points (m).
%--------------------------------------------------------------------------
%
% Compute depths (m) at RHO- and Wpoints.

Zw = depths(fname, gname, 5, 0, tindex);

N = size(Zw,3) - 1;

% Compute level thickness (m).

Hz = Zw(:,:,2:N+1) - Zw(:,:,1:N);

return
