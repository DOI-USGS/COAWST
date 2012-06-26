%
% ROMS Land/Sea Mask Utility
% ==========================
%
% This utility is a graphical User interface enabling the User to easily
% create and modify a land mask for ROMS. To accelerate the proccessing,
% the Land/Sea mask is edited in (I,J) grid coordinates. This requires a
% convertion of coastline data (used in SeaGRid) to (I,J) indices.  This
% utility calls the MEXNC interface for reading and writing to NetCDF
% files.
%
%
% Drivers:
%
%   editmask     - Interactive ROMS Land/Sea mask editing driver.
%   editscope    - Interactive ROMS adjoint sensitivity scope mask editor.
%   landsea      - Authomatic ROMS Land/Sea processing.
%
% Input/Output:  NetCDF I/O processing in ROMS GRID file
%
%   read_mask    - Reads in Land/Sea mask data.
%   read_scope   - Reads in adjoint sensitivity scope mask.
%   write_mask   - Writes out Land/Sea mask data.
%   write_scope   - Writes out adjoint sensitivity scope mask data.
%
% Land/Sea mask:
%
%   uvp_masks    - Computes the Land/Sea mask data on U-, V-, and PSI-points.
%
% Menu Interface:
%
%   axisscroll   - Draws horizontal or vertical scroll bars.
%   button       - Creates a menu button.
%   pointer      - Sets custon pointer.
%   radiobox     - Creates a group of radio bottons.
%   textbox      - Creates a textbox in a frame.
%
% Orthogonal Grid Coordinate Manipulation:
%
%   ijcoast      - Converts coastline (lon,lat) coordinates to fractional
%                  coordinates.
%
% Miscellaneous:
%
%   pltmask      - Plots Land/Sea masks.
%   pltscope     - Plots adjoint sensitivity scope masks.
%

% svn $Id: Contents.m 436 2010-01-02 17:07:55Z arango $
%===========================================================================%
%  Copyright (c) 2002-2010 The ROMS/TOMS Group                              %
%    Licensed under a MIT/X style license                                   %
%    See License_ROMS.txt                           Hernan G. Arango        %
%===========================================================================%
