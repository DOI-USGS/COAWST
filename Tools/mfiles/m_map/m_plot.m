function h = m_plot(varargin)

% M_PLOT Plot objects on an M_MAP plot.   All of the normal
% Matlab plot options are available.  

% USAGE: M_PLOT(LON,LAT,[OPTIONS]) 

% Deirdre Byrne, dbyrne@umeoce.maine.edu 00/06/23
%
% This software is provided "as is" without warranty of any kind.

global MAP_PROJECTION MAP_VAR_LIST

if isempty(MAP_PROJECTION),
  disp('No Map Projection initialized - call M_PROJ first!');
  return;
end;

if nargin < 2;
  help m_plot
  return
end

[x,y] = m_ll2xy(varargin{1},varargin{2});
varargin = varargin(:);
s = size(varargin,1);
h=plot(x,y,varargin{3:s});

if nargout == 0
  clear h
end

