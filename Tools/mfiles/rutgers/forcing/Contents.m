%
% ROMS Forcing Fields
% ===================
%
% These functions are used for preparing ROMS forcing NetCDF files.
%
%
%   add_heatflux  - Adds the heat flux components (shortwave, longwave,
%                     latent, and sensible) variables to an existent
%                     NetCdF file if not present.
%   d_core2_frc   - Driver template script showing how to create ROMS
%                     forcing NetCDF file(s) using ROMS metadata
%                     structure. The data source is the CORE 2 Global
%                     Air-Sea Flux Dataset. Notice that the original
%                     data set is sampled for the Gulf of Mexico (GOM).
%
%   d_ecmwf2roms  - Driver template script showing how to create ROMS
%                     forcing NetCDF file(s) using ROMS metadata
%                     structure. The data source is the ECMWF ERA-Interim
%                     (Jan 1979 - present).  The example is provided
%                     for the Gulf of Mexico (GOM).
%
%   otps2roms     - Generates a ROMS tidal forcing from the OTPS dataset.
%                     The 'base_date' is related to ROMS input parameter
%                     the TIDE_START. It is used in ROMS to compute the
%                     correct phase lag with respect initialization time.
%                     I still not happy how this is done and will revisit
%                     this function in the future to get the data directly
%                     without the awkward Fortran interface.
%
%  write_tides    - Creates ROMS tidal forcing NetCDF file and writes
%                     data extracted from either OTPS or ADCIRC and
%                     processed with the "t_tide" utility.
%

% svn $Id: Contents.m 996 2020-01-10 04:28:56Z arango $
%=========================================================================%
%  Copyright (c) 2002-2020 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%
