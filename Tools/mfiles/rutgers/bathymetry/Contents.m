%
% Bathymetry Extracting Functions
% ===============================
%
% These functions are used for extracting bathymetry from a dataset. 
%
%
%   etopo5.nc     - ETOPO5 Earth surface topography dataset
%                   (5 minute resolution).
%
%   c_bath        - Creates a bathymetry NetCDF file.
%   extract_bath  - Driver to extract bathymetry.  The extracted data is
%                     either written into a Matlab file that can be use in
%                     "seagrid" or a NetCDF file.
%   get_bath      _ Extract bathymetry an elevation from a ETOPO NetCDF
%                     file. It is similar to the "extract_bath" driver but
%                     in function format.
%   x_etopo5      - Extract requested bathymetry from ETOPO-5 dataset.
%

% svn $Id: Contents.m 996 2020-01-10 04:28:56Z arango $
%===========================================================================%
%  Copyright (c) 2002-2020 The ROMS/TOMS Group                              %
%    Licensed under a MIT/X style license                                   %
%    See License_ROMS.txt                           Hernan G. Arango        %
%===========================================================================%
