%
% ROMS Miscellaneous Matlab scripts
% =================================
%
% This utility contains several generic Matlab scripts to pre- and
% post-processing ROMS data.
%
% Terrain-Following Coordinates:
%
%   depths        - Computes ROMS depths associated with a 3D variable in
%                     a NetCDF file.
%   hslope        - Computes and plot ROMS grid bathymetry slope.
%   rfactor       - Computes bathymetry stiffness ratio, r-factor.
%   scoord        - Computes and plot ROMS vertical stretched coordinates.
%   set_depth     - Computes ROMS depths for a 3D variable during
%                     pre-processing, like initial conditions, climatology,
%                     etc.
%   smooth_bath   - Smooths bathymetry as function of the r-factor.
%   stretching    - Computes ROMS vertical coordinate stretching function.
%
% ROMS Data Processing:
%
%   nanland       - Mask Land points to NaN to facilitate plotting.
%   plot_field    - Plots requested ROMS variable from input history
%                     NetCDF file. This function is very usefull when
%                     debugging a ROMS application. The ploting is not
%                     that fancy but it provides enough information for
%                     browsing ROMS variables very quickly. The location
%                     of the minimum is marked with a filled magenta
%                     circle whereas the maximum is marked with a filled
%                     magenta square.
%   roms_vectors  - Processes vector data for either the full grid
%                     or boundary edges. The strategy is to get any
%                     horizontal vector field at RHO-points for the
%                     event that a rotation to ROMS curvilinear grid
%                     is needed.  Then, they are computed at the
%                     appropriate Arakawa C-grid location.
%   rotate_vec    - Rotates vector components from TRUE East and North
%                     to curvilinear coordinates (XI,ETA) orientation
%                     or viceversa. Input vector components may be
%                     located at the center of the cell (RHO-points)
%                     or at staggered Arakawa's C-grid locations.
%   sample_grid   - Gets Parent grid indices range of the polygon that
%                     tightly contains the Target Grid.
%
%   uv_barotropic - Computes vertically integrated velocity components
%                      for ROMS full grid or boundaries.
%
% Filters:
%
%   shapiro1      - 1D Shapiro filter.
%   shapiro2      - 2D Shapiro filter.
%
% Geophysical:
%
%   eos           - ROMS equation of state for seawater.
%   gcircle       - Great circle distance between two (lon,lat) points.
%   geodesic_dist - Geodesic distance between two (lon,lat) points.
%   roms_eos.m    - Computes 'in situ' density using ROMS nonlinear
%                     equation of state for seawater. It assumes no
%                     pressure variation along geopotential surfaces,
%                     that is, depth (meters; negative) and pressure
%                     (dbar; assumed negative here) are interchangeable.
%
% Time Management:
%
%   caldate      - Converts Julian day number to calendar date structure.
%   date_stamp   - Sets current date string.
%   day_code     - Computes day of the week for a given date.
%   gregorian    - Converts Julian day number to Gregorian calendar date.
%   greg2str     - Converts Gregorian date array to string.
%   hms2h        - Converts hours, minutes, and seconds to decimal hours.
%   julian       - Converts Gregorian calendar date to Julian day numbers.
%   s2hms        - Converts decimal seconds to integer hour, minute,
%                    seconds.
%
% Parallelism:
%
%   ptile        - Plot (overlay) ROMS horizontal tile partitions.
%   tile         - Compute ROMS horizontal tile partitions indices.
%
% Plotting:
%
%   pcolorjw       - modified version of "pcolor" the recomputes X,Y at the
%                      mid points of the input, and pads C so that the
%                      effect is to shift each square one half space to
%                      the lower left.  The perimeter squares are only half
%                      the width/height of all the others as can be seen if
%                      one sets shading('faceted').
%
%   plot_nesting   - plots requested ROMS nesting variable from input
%                      history NetCDF files. This function is very useful
%                      when debugging a ROMS nesting application. The
%                      plotting is not that fancy but it provides enough
%                      information for browsing ROMS variables very quickly.
%
%   plot_perimeter - Adds a grid perimeter outline to an existing figure
%                      plotted with 'plot_nesting'. 
%

% svn $Id: Contents.m 735 2014-04-28 23:15:37Z arango $
%=========================================================================%
%  Copyright (c) 2002-2014 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%
