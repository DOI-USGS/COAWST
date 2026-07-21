function [F]=merge_field(Slon, Slat, S, Dlon, Dlat, D)

% MERGE_FIELD:  Merges fields from two different datasets.
%
% [F]=merge_field(Slon, Slat, S, Dlon, Dlat, D)
%
% This function merges source and destination fields from two different
% datasets and grids. The mergin only takes in the common are between
% grids.
%
% On Input:
%
%    Slon        Source field longitude (2D array)
%    Slat        Source field latitude  (2D array)
%    S           Source field (2D or 3D array)
%    Dlon        Destination field longitude (2D array)
%    Dlat        Destination field latitude  (2D array)
%    D           Destination field (2D or 3D array)
%
% On Output:
%
%    F           Merged destination field (structure)
%
%                  F.index   common destination water points linear indices
%                  F.common  interpolated field at common water points
%                  F.merged  merged source and destination fields
%                  

% svn $Id$
%=========================================================================%
%  Copyright (c) 2002-2025 The ROMS Group                                 %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.md                            Hernan G. Arango      %
%=========================================================================%
  
% Find source and destination fields NaN values corresponding to land/sea
% mask.

Sland =  isnan(S);
Ssea  = ~isnan(S);

Dland =  isnan(D);
Dsea  = ~isnan(D);

% Set source field polygonal region.

[Im,Jm] = size(Slon);

Slon_perimeter = [squeeze(Slon(1:Im,1));                                ...
                  squeeze(Slon(Im,2:Jm))';                              ...
                  squeeze(flipud(Slon(1:Im-1,Jm)));                     ...
                  squeeze(fliplr(Slon(1,1:Jm-1)))'];

Slat_perimeter = [squeeze(Slat(1:Im,1));                                ...
                  squeeze(Slat(Im,2:Jm))';                              ...
                  squeeze(flipud(Slat(1:Im-1,Jm)));                     ...
                  squeeze(fliplr(Slat(1,1:Jm-1)))'];

% Find destination points inside or at the boundary of the source field
% perimeter.

[IN,ON] = inpolygon(Dlon(:), Dlat(:),                                   ...
                    Slon_perimeter, Slat_perimeter);

IN(ON) = true;                   % add points on perimeter
IN(Dland) = false;               % remove land points

% Interpolate source field points to destination points (only common
% water points.

Fi = scatteredInterpolant(Slon(Ssea), Slat(Ssea), S(Ssea));

Di = nan(size(D));
Di(IN) = Fi(Dlon(IN), Dlat(IN));

% Meld source and destination fields: replace values at source points.

D(IN) = 0.3 .* D(IN)+ 0.7 .* Di(IN);

% Set output melded field structure.

F=struct('index', [], 'common', [], 'merged', []);

F.index  = IN;
F.common = Di;
F.merged = D;

return

