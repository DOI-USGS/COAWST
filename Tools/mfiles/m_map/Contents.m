% M_Map - mapping toolbox (Author: rich@ocgy.ubc.ca)
% Version 1.3 6/Nov/2000
%
% You have collected your data, loaded it into Matlab, analyzed 
% everything to death, and now you want to make a simple map showing 
% how it relates to the world. 
%
% But you can't. 
%
% Instead you have to figure out how to save all your data, and 
% then read it into a mapping program, and then spend all that extra 
% time figuring out why the mapping program doesn't give you what
% you expected it would...
%
% No more! 
%
%                            Announcing M_Map v1.3f! 
%
% M_Map is a set of mapping tools written for Matlab v5. These include: 
%
%    1. Routines to project data in 18 different spherical 
%       projections (and determine inverse mappings) 
%    2. A grid generation routine to make nice axes with 
%       limits either in long/lat terms or in planar
%       X/Y terms. 
%    3. A coastline database (with 1/4 degree resolution) 
%    4. A global elevation database (1 degree resolution) 
%    5. Hooks into freely available high-resolution coastlines and
%       bathymetry/topography.
%
% M_Map v1.3 is available via the web at 
%
%       http://www.ocgy.ubc.ca/~rich/
%
%
% Toolbox contents
%
%    Contents.m    - This file
%    m_demo.m      - demonstrates a few different maps.
%
%  User-callable functions
%
%    m_proj.m      - initializes projections
%
%    m_grid.m      - draws grids 
%    m_scale       - forces map to a given scale.
%
%    m_ungrid.m    - erases map elements (if you want to change parameters)
%
%    m_coast.m     - draws a coastline
%    m_elev.m      - draws elevation data from 1 degree database
%
%    m_tbase.m     - draws elevation data from 5-minute TerrainBase database
%    m_gshhs_c.m   - draws coastline from GSHHS crude database
%    m_gshhs_l.m   - draws coastline from GSHHS low-resolution database
%    m_gshhs_i.m   - draws coastline from GSHHS intermediate-resolution database
%    m_gshhs_h.m   - draws coastline from GSHHS high-resolution database
%    m_gshhs_f.m   - draws coastline from GSHHS full database
%    m_plotbndry.m - draws a political boundary from the DCW 
%    m_usercoast.m - draws a coastline using a user-specified subset database.
%
%    m_plot.m      - draws line data in map coords
%    m_line.m      - draws line data in map coords
%    m_text.m      - adds text data in map coords
%    m_legend.m    - Draw a legend box
%    m_quiver      - draws arrows for vector data
%    m_contour     - draws contour lines for gridded data
%    m_contourf    - draws filled contours
%    m_patch       - draws patch data
%    m_track       - draws annotated tracklines
%    m_range_ring  - draws range rings
%
%    m_ll2xy.m     - converts from long/lat to map coordinates
%    m_xy2ll.m     - converts from map coordinates to long/lat
%
%    m_lldist      - distance between points (long/lat coordinates)
%    m_xydist      - distance between points (map projection coordinates)
%
%    m_tba2b.m     - used in installing high-resolution elevation database.
%
%    m_vec.m       - fancy arrows
%
%  Internal functions (not meant to be user-callable)
%
%    private/mp_azim.m   - azimuthal projections
%    private/mp_cyl.m    - cylindrical projections (equatorial)
%    private/mp_conic.m  - conic projections
%    private/mp_tmerc.m  - transverse cylindrical projections
%    private/mp_utm.m    - elliptical universal transverse cylindrical projections
%    private/mp_omerc.m  - oblique cylindrical projection
%
%    private/mu_util.m   - various utility routines
%    private/mu_coast.m  - routines to handle coastlines.
%
%    private/clabel.m    - patched version of clabel 
%                         (matlab v5.1 version does not contain
%                         capabilities for different text properties).    
%
%    private/m_coasts.mat - coastline data
%
%  HTML documentation
%
%    map.html           - Home page, examples
%    private/mapug.html - User's guide
%    private/*gif       - examples.
%  
%
% Questions or problems; email me - rich@ocgy.ubc.ca.
%
% Rich Pawlowicz
% Oceanography, Dept. of Earth and Ocean Sciences, Univ. of British Columbia, 
% 6270 University Blvd., Vancouver, B.C. CANADA V6T 1Z4
% email: rich@ocgy.ubc.ca 


    

