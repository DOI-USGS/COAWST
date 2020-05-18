%
% ROMS Initial Conditions and Climatology
% =======================================
%
% These functions are used for preparing ROMS initial conditions
% and climatology NetCDF file.
%
%
%   c_biology       - Defines ROMS biology variables to an existing
%                       initial conditions or climatology NetCDF file.
%
%   c_climatology   - Creates ROMS climatology NetCDF file.
%
%   c_initial       - Creates ROMS initial conditions NetCDF file. 
%
%   d_climatology   - Driver to process generic ROMS climatology file.
%
%   d_initial       - Driver ro process generic ROMS initial conditions
%                       file.
%
%   d_mercator2roms - Driver template to process ROMS initial conditions
%                       using the 'mercator2roms' script.
%
%   d_nudgcoef      - User modifiable script that can be used to prepare
%                       ROMS nudging inverse time scales NetCDF file.
%
%   d_oa2roms.m     - User modifiable script showing how to create ROMS
%                       climatology NetCDF file(s). The data source are
%                       objective analyzed (OA) Annual and Monthly
%                       temperature and salinity fields from the Levitus
%                       (1998) dataset. These fields are generated using
%                       ROMS OA package. 
%
%   d_roms2roms     - Driver template to process ROMS initial conditions
%                       using the 'roms2roms' script.
%
%   interp_field    - Interpolates a generic ROMS 2D or 3D field variable
%                       from a Donor to Receiver Grid. If 3D interpolation,
%                       the Donor Grid data is interpolated first to the
%                       Receiver Grid horizontal locations using
%                       'TriScatteredInterp' at each of the Donor Grid
%                       vertical levels. Then, 'interp1' is used to
%                       interpolate to Receiver Grid vertical locations.
%
%   mercator2roms   - Interpolates requested variable data from
%                       Mercator to ROMS grids.
%
%   oa_cat.m        - Reads Annual and Monthly objective analysis (OA)
%                       files of the Levitus climatology and appends
%                       bottom Annual levels to Monthly OA fields. The
%                       Monthly Levitus Climatology is available from
%                       the surface to 1000m (1994 dataset) or 1500m
%                       (1998 dataset). The annual fields are used for
%                       the missing bottom levels. 
%
%   oa2roms.m       - Vertically interpolates requested OA variable
%                       to ROMS terrain-following vertical grid. The
%                       OA package generates qfields on a constant
%                       depth: standard depth levels. 
%
%   roms2roms       - Interpolates requested variable data from
%                       ROMS to ROMS grids.

% svn $Id: Contents.m 996 2020-01-10 04:28:56Z arango $
%=========================================================================%
%  Copyright (c) 2002-2020 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%
