%
% Coupling Melding Weight Factors:
% ===============================
%
% This directory contains functions to compute melding weights between
% DATA and ROMS coupling components for exporting fields to the atmosphere
% component.
%
% c_weithts.m       - Creates a new weight factor NetCDF file for melding
%                       DATA and OCEAN components during ESMF coupling.
% coamps_weights.m  - Computes and writes the COAMPS melding weights to
%                       combine fields from DATA and ROMS components after
%                       ESMF regridding to the atmosphere grid.
% smooth_weights.m  - Computes smooth melding weights to combine fields
%                       from DATA and OCEAN components. The merging
%                       factors change gradually in an area next to the
%                       OCEAN grid open boundary.
% wrf_weights.m     - Computes and writes the WRF melding weights to
%                       combine fields from DATA and ROMS components after
%                       ESMF regridding to the atmosphere grid.
%

% svn $Id: Contents.m 996 2020-01-10 04:28:56Z arango $
%=========================================================================%
%  Copyright (c) 2002-2020 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%
