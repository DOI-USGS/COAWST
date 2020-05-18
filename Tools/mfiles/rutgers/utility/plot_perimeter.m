function h=plot_perimeter(G, LineType)

%
% PLOT_PERIMETER:  Adds domain perimeter to an existing figure
%
% h=plot_perimeter(G)
%
% This function is used to add a grid perimeter outline to an
% existing figure plotted with 'plot_nesting'.  It is used for
% browsing quickly output ROMS variables during nesting.
%
% On Input:
%
%    G            An existing ROMS grid structure (struct array)
%
%    LineType     Line color and type (string)
%
% On Output:
%
%    h            Plot handler (vector)
%

% svn $Id: plot_perimeter.m 996 2020-01-10 04:28:56Z arango $
%=========================================================================%
%  Copyright (c) 2002-2020 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%

% Check if input variable G is a structure.
  
if (~isstruct(G))
  disp(blanks(1));
  error('Input grid variable is not a structure...');
end
  
% Initialize line color and type.

if (nargin < 2)
  LineType = 'r-';
end

% Draw perimeter to current figure.

hold on;

Ngrids = length(G);

h = zeros([1 Ngrids]);

for ng = 1:S.Ngrids,
  if (G(ng).spherical),
    h(ng) = plot(G(ng).lon_perimeter,                                   ...
                 G(ng).lat_perimeter,                                   ...
                 LineType, 'LineWidth', 2);
  else
    h(ng) = plot(G(ng).x_perimeter,                                     ...
                 G(ng).y_perimeter,                                     ...
                 LineType, 'LineWidth', 2);
  end
end

hold off

return

 
