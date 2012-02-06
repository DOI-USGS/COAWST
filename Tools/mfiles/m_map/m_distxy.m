function dist = m_distxy(varargin)

% M_DISTXY Distance between X,Y locations on a map.  Distances
% are given in meters.

% USAGE: M_DIST(X,Y), for example 
%
%        G = GINPUT(2); km_dist = M_DIST(G(:,1),G(:,2))/1000;

% Deirdre Byrne, dbyrne@umeoce.maine.edu 00/06/23
%
% This software is provided "as is" without warranty of any kind.

global MAP_PROJECTION MAP_VAR_LIST

if isempty(MAP_PROJECTION),
  disp('No Map Projection initialized - call M_PROJ first!');
  return;
end;

if nargin < 2;
  help m_distxy
  return
end

varargin = varargin(:);
s = size(varargin,1);

% get lon, lat coords
[lon,lat] = m_xy2ll(varargin{1},varargin{2});

% calculate distance in meters
dist = geodist(lat,lon,varargin{3:s});

if size(lon,2) == 1;
    dist = dist(:);
end

return
