function G = grids_structure(Gnames, Hnames)

%
% G = grids_structure(Gnames);
% G = grids_structure(Hnames);
% G = grids_structure(Gnames, Hnames);
%
% This function builds ROMS nested grids structure array, G(:), containing
% all the variables associated with the application's horizontal and
% vertical grids.
%  
% On Input:
%
%    Gnames        ROMS Grid or History NetCDF file/URL names containing
%                    all grid variables (cell array)
%
%    Hnames        ROMS History NetCDF file/URL names containing
%                    vertical grid information (OPTIONAL; cell array)
%
% On Output:
%
%    G(:)          Nested grids structure (1 x Ngrid struct array)
%
%
% If the coastline information is needed in the output structure array,
% you need to provide the Grid NetCDF files in Gnames. These data is not
% available in histr=ory files. Of course, the coastline variables
% (lon_coast, lat_coast) need to be in the Grid NetCDF file(s).  It is
% very good idea to have the coastline data used when processing grids.
% The script "add_coastline.m" can be used to append such data.
%
% If the Grid NetCDF files are provided in Gnames, you will need to
% provide the history files in Hnames in order to process vertical
% grid arrays.  This information is not available in the Grid NetCDF
% files.
%
% Example:
%
%    G = grids_structure({'my_grd_coarse.nc',  ...
%                         'my_grd_fine1.nc',   ...
%                         'my_grd_fine2.nc'},  ...
%                        {'my_his_coarse.nc',  ...
%                         'my_his_fine1.nc',   ... 
%                         'my_his_fine2.nc'});
%
% Calls to External Functions:
%
%    get_roms_grid     Gets Information Grids Structure, G(ng) 
%

% svn $Id: grids_structure.m 738 2014-10-14 21:49:14Z arango $
%=========================================================================%
%  Copyright (c) 2002-2014 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%

% Initialize.
			  
if (nargin > 1),
  got_his = true;
else
  got_his = false;
end

%--------------------------------------------------------------------------
% Get nested grid structures.
%--------------------------------------------------------------------------  

% If the grid structure have the parent fields, remove them to have an
% array of similar structures.

parent = {'parent_grid',                                                ...
          'parent_Imin', 'parent_Imax',                                 ...
          'parent_Jmin', 'parent_Jmax'};

Ngrids = length(Gnames);

if (got_his),
  for n=1:Ngrids,
    g = get_roms_grid(char(Gnames(n)), char(Hnames(n)));
    if (isfield(g, 'parent_grid')),
      G(n) = rmfield(g, parent);
    else
      G(n) = g;
    end
  end
else
  for n=1:Ngrids,
    g = get_roms_grid(char(Gnames(n)));
    if (isfield(g, 'parent_grid')),
      G(n) = rmfield(g, parent);
    else
      G(n) = g;
    end
  end
end

return
