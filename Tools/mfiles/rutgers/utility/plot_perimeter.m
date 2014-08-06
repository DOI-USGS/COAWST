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

% svn $Id: plot_perimeter.m 735 2014-04-28 23:15:37Z arango $
%=========================================================================%
%  Copyright (c) 2002-2014 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%

% Check if input variable G is a structure.
  
if (~isstruct(G)),
  disp(blanks(1));
  error('Input grid variable is not a structure...');
end
  
% Initialize line color and type.

if (nargin < 2),
  LineType = 'r-'
end

% Get grids perimeter structure.

S = grid_perimeter(G);

% Draw perimeter to current figure.

hold on;

h = zeros([1 S.Ngrids]);

for ng = 1:S.Ngrids,
  h(ng) = plot(S.grid(ng).perimeter.X_psi,                              ...
               S.grid(ng).perimeter.Y_psi,                              ...
               LineType, 'LineWidth', 2);
end

hold off

return

 
