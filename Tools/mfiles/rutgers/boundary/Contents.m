%
% ROMS Open Boundary Conditions
% =============================
%
% These functions are used for preparing ROMS open boundary conditions
% NetCDF file.
%
%
%   c_boundary      - Creates ROMS open boundary conditions NetCDF
%                       file.
%
%   d_obc_mercator  - Driver to extract open boundary conditions from
%                       Mercator dataset. It creates open boundary
%                       conditions NetCDF file, interpolates data to
%                       application grid, and writes out data. It is
%                       a template showing how to use the 'obc_mercator'
%                       script.
%   d_obc_roms2roms - Driver to extract open boundary conditions from
%                       a ROMS dataset. It creates open boundary
%                       conditions NetCDF file, interpolates data to
%                       application grid, and writes out data.  It is
%                       a template showing how to use the 'obc_roms2roms'
%                       script.
%
%   interp_boundary - Interpolates lateral boundary conditions for a
%                       ROMS generic 2D or 3D state variable from a
%                       Donor to Receiver. If 3D interpolation, the
%                       Donor Grid data is interpolated first to the
%                       Receiver Grid horizontal locations using
%                       'TriScatteredInterp' at each of the Donor
%                       Grid vertical levels. Then, 'interp1' is used
%                       to interpolate to Receiver Grid vertical
%                       locations.
%
%   plot_boundary   - Plots requested ROMS variable from input lateral
%                       boundary NetCDF file. This function is very
%                       useful for browsing 3D lateral boundary
%                       conditions for ROMS in a vertical slab. The
%                       location of the minumum is marked with a filled
%                       magenta circle whereas the maximum is marked
%                       with a filled magenta square. 
%
%   obc_mercator    - Interpolates requested variable open boundary
%                       conditions from Mercator to ROMS grid boundary
%                       edges.
%   obc_roms2roms   - Interpolates requested variable open boundary
%                       conditions from ROMS to ROMS grids.
%
%   extract_bry     - Reads requested variable from a ROMS NetCDF file
%                       at the specified time record and extracts the
%                       lateral boundary edges. No interpolation is
%                       carried out.
%

% svn $Id: Contents.m 996 2020-01-10 04:28:56Z arango $
%=========================================================================%
%  Copyright (c) 2002-2020 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%
