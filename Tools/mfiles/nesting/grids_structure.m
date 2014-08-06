function G = grids_structure(Gnames)

%
% G = grids_structure(Gnames)
%
% This function builds ROMS nested grids structure array, G(:), containing
% all the variables associated with the application's horizontal and
% vertical grids.
%
% On Input:
%
%    Gnames        ROMS Grid/History NetCDF file/URL names containing
%                    all grid variables (cell array)
% On Output:
%
%    G(:)          Nested grids structure (1 x Ngrid struct array)
%
%
% Example:
%
%    G = grids_structure({'my_his_coarse.nc',  ...
%                         'my_his_fine1.nc',   ... 
%                         'my_his_fine2.nc'});
%
% Calls to External Functions:
%
%    get_roms_grid     Gets Information Grids Structure, G(ng) 
%

% svn $Id: grids_structure.m 711 2014-01-23 20:36:13Z arango $
%=========================================================================%
%  Copyright (c) 2002-2014 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%

%--------------------------------------------------------------------------
% Get nested grid structures.
%--------------------------------------------------------------------------  

% If the grid structure have the parent fields, remove them to have an
% array of similar structures.

parent = {'parent_grid',                                                ...
          'parent_Imin', 'parent_Imax',                                 ...
          'parent_Jmin', 'parent_Jmax'};

Ngrids = length(Gnames);

for n=1:Ngrids,
  g = get_roms_grid(char(Gnames(n)));
  if (isfield(g, 'parent_grid')),
    G(n) = rmfield(g, parent);
  else
    G(n) = g;
  end
end

return
