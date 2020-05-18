%
% Coastline Extracting Functions
% ==============================
%
% These functions are used for extracting coastline data from the GSHHS
% dataset.
%
%
%   boundary.dat  - Example file to load boundary edges to "seagrid".
%
%   extract_coast - Template to extract coastline.  The extracted data is
%                     either written into a Matlab file that can be use in
%                     "seagrid" or an ASCII file for ROMS plotting package.
%   get_coast     - Extracts coastline data from the Global Self-consistent
%                     Hierarchical, High-resolution Shoreline (GSHHS)
%                     database for the box corners bounded by (Llon,Blat)
%                     and (Rlon,Tlat).
%   grid_corners  - Template to set the grid corners for "seagrid".
%
%   r_gshhs       - Read and extract coastline data from GSHHS dataset.
%   x_gshhs       - Process extracted GSHHS coastline data.
%

% svn $Id: Contents.m 996 2020-01-10 04:28:56Z arango $
%=========================================================================%
%  Copyright (c) 2002-2020 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%
