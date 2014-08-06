%
% ROMS NetCDF Proccessing Scripts
% ===============================
%
% Starting Matlab Version 2012a, released Feb 9, 2012, the 'native' Matlab
% interface to NetCDF is the preferred method for processing NetCDF files.
% This version has OpenDAP support. Otherwise, the MEXNC toolbox or the
% SNCTOOLS java interface for OpenDAP files are used for older versions
% Matlab.  The 'native' NetCDF Matlab support started in version 2008b.
% 
% The NetCDF scripts for MEXNC toolbox that can be downloaded from:
%
%   svn checkout url ~/matlab/mexnc
%
% where
%
%   url = https://mexcdf.svn.sourceforge.net/svnroot/mexcdf/mexnc/trunk
%
% Similarly, the SNCTOOLS uses the NetCDF Java interface for OpenDAP
% that can be downloaded from:
%
%   svn checkout url ~/matlab/snctools
%
% where
%
%   url = https://mexcdf.svn.sourceforge.net/svnroot/mexcdf/snctools/trunk
%
% Both the MEXNC and SNCTOOLS intefaces are becoming obsolete and I doubt
% that they will further developed in the future.
%
% NetCDF I/O Processing:
%
%   nc_append      - Defines/appends new variables to an existing NetCDF
%                      file. The variable(s) metadata is provided in input
%                      structure, which includes information about variable
%                      type, dimensions and attributes.
%   nc_attadd      - Adds/modifies a global or variable NetCDF attribute.
%   nc_attdel      - Deletes requested global or variable NetCDF attribute.
%   nc_check       - Checks the information structure returned from calls
%                      to "nc_inq" or native "ncinfo" for compliance and
%                      changes variables types and attributes.
%   nc_constant    - Gets numeric value of a named NetCDF library contant.
%   nc_create      - Creates a new NetCDF file according to the file
%                      creation mode flag. If a structure S is provided,
%                      it defines the dimensions, global attributes,
%                      variables and attributes stored in the structure.
%                      This structure can be created using "nc_inq" or
%                      nativ function "ncinfo".
%   nc_dinfo       - Inquires about the dimensions in a NetCDF file.
%   nc_drename     - Renames a NetCDF dimension.
%   nc_getatt      - Gets a global or variable NetCDF attribute.
%   nc_inq         - Inquires about the contents of a NetCDF file.
%   nc_interface   - Sets what NetCDF interface to use: 'native' Matlab
%                      interface (default), 'java' SNCTOOLS interface for
%                      OpenDAP files (Matlab version < 2012a), or 'mexnc'
%                      interface (Matlab Verison < 2008b).
%   nc_test        - Creates a NetCDF using data from the peaks(40)
%                      function. Several datatype variables are created
%                      to test the NetCDF interface in Matlab.
%   nc_vdef        - Creates a ROMS variable in a NetCDF file.
%   nc_vinfo       -  Inquires information about requested NetCDF variable.
%   nc_vname       - Gets the names of all variables in a NetCDF file. This
%                      function is OBSOLETE but it is kept for backward
%                      compatibility. Use "nc_vnames" in the future.
%   nc_vnames      - Gets the names of all variables in a NetCDF file.
%   nc_vrename     - Renames a NetCDF variable.
%
%   nc_read        - Generic function to read  requested NetCDF variable.
%   nc_write       - Generic function to write requested NetCDF variable.
%
%   nc_slice       - Interpolates requested slice from a 3D NetCDF variable.
%
% ROMS NetCDF Metadata Structure:
%
%   check_metadata - Checks ROMS metadata structure for consistency and
%                      fills unassigned fields. This structure will be
%                      used elsewhere to create NetCDF files in a
%                      compact way.
%
%   roms_metadata  - Sets metadata structure for requested ROMS NetCDF
%                      variable. The structure contains the same fields
%                      that are returned by 'nc_inq' or native interface
%                      'ncinfo'. This can be used elsewhere to create
%                      NetCDF a file in a compact way.
%

% svn $Id: Contents.m 711 2014-01-23 20:36:13Z arango $
%=========================================================================%
%  Copyright (c) 2002-2014 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%
