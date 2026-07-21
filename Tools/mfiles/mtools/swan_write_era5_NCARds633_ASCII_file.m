% Script swan_write_era5_NCARds633_ASCII_file.m
%
% This is based on John L. Wilkin's code 
%
%     roms_write_era5_NCARds633_frcfile.m
%
% -------------------------------------------------------------------------
%
% Create a ROMS meteorology forcing file from ERA5 reanalysis extracted
% from NCAR dataset ds633.0 using roms_get_era5_ncar_ds633.m, i.e. ...
%
% First run e.g. E = roms_get_era5_ncar_ds633(year,month,bounding_box...)
% and then run this script to create the SWAN forcing ASCII file.
%
% See the help on roms_get_era5_ncar_ds633 regarding obtaining login
% credentials to access the ERA5 archive at NCAR Research Data Archive.
%
% -------------------------------------------------------------------------
%
% NOTE: This is an example script. The user needs to configure the ROMS
% time coordinate basedate (Time0) and the output filename (ncname) and
% title and any other pertinent metadata.
%
% This script uses routines from the myroms.org Matlab tools
% (roms_metadata, nc_create, nc_write, nc_constant etc.).
%
% John Wilkin - December 2020
%
% Copyright (c) - 2021 John L. Wilkin - jwilkin@rutgers.edu
% $Id: roms_write_era5_NCARds633_frcfile.m 600 2020-12-29 19:19:51Z wilkin $
%
% Obtain an up-to-date version of this code from 
% https://github.com/johnwilkin/roms_wilkin
%
% See also roms_get_era5_NCARds633_bulkflux
%
%
% -------------------------------------------------------------------------
%
%  Updated by Salme E Cook - February 2021 - secook@usgs.gov
%
% ------------------------------------------------------------------------
%

% USER SETS PARAMETERS IN THIS BLOCK

% shift time to the ROMS basedate you want to use
Time0 = E.time.data(1);
time = E.time.data - Time0;

% yyyy and mm are inherited from the call to roms_get_era5_ncar_ds633
YYYY = upper(int2str(E.yyyy));
MM = upper(sprintf('%02d',E.mm));

% yyyy and mm are inherited from the call to roms_get_era5_ncar_ds633
YYYY = upper(int2str(E.yyyy));
MM = upper(sprintf('%02d',E.mm));

% Set the output file name prefix
SWAN_forc_name=strcat('swan_ERA5_',YYYY,MM,'.dat');


% SWAN input file needs the following in an ASII file format
%
% For Example: A Regular Input Wind Grid and Nonstationary in time
%
% REGULAR WIND grid
%
% INPGRID WIND REGULAR options
% see: http://swanmodel.sourceforge.net/online_doc/swanuse/node26.html
%
%     [xpinp]	geographic location (x -coordinate) of the origin of the input grid in	 
%         problem coordinates (in m) if Cartesian coordinates are used or in degrees if	 
%         spherical coordinates are use (see command COORD).	 
%         Default: [xpinp] = 0. In case of spherical coordinates there is no default, the	 
%         user must give a value.	 
%     [ypinp]	geographic location (y -coordinate) of the origin of the input grid in	 
%         problem coordinates (in m) if Cartesian coordinates are used or in degrees if	 
%         spherical coordinates are use (see command COORD).	 
%         Default: [ypinp] = 0. In case of spherical coordinates there is no default, the	 
%         user must give a value.	 
%     [alpinp]	direction of the positive x -axis of the input grid (in degrees, Cartesian convention).	 
%         See command COORD.	 
%         Default: [alpinp] = 0.	 
%     [mxinp]	number of meshes in x -direction of the input grid (this number is one less	 
%         than the number of grid points in this direction!).	 
%     [myinp]	number of meshes in y -direction of the input grid (this number is one less	 
%         than the number of grid points in this direction!).	 
%         In 1D-mode, [myinp] should be 0.	 
%     [dxinp]	mesh size in x -direction of the input grid,	 
%         in m in case of Cartesian coordinates or	 
%         in degrees if spherical coordinates are used, see command COORD.	 
%     [dyinp]	mesh size in y -direction of the input grid,	 
%         in m in case of Cartesian coordinates or	 
%         in degrees if spherical coordinates are used, see command COORD.	 
%         In 1D-mode, [dyinp] may have any value.	 
%         Default: [dyinp] = [dxinp].
%
% NONSTATIONARY options
% 
%     [tbeginp]	begin time of the first field of the variable, the format is:	 
%         ISO-notation 19870530.153000	 
%
%     [deltinp]	time interval between fields, the unit is indicated in the next option:	 
%         SEC unit seconds	 
%         MIN unit minutes	 
%         HR unit hours	 
%         DAY unit days	 
%
%     [tendinp]	end time of the last field of the variable, the format is:	 
%           ISO-notation 19870530.153000	 

% ------------------------------------------------------------------------
%

% Coordinates and all the data necessary to write this file are in
% structure E loaded by roms_get_era5_ncar_ds633 for a requested year,
% month, and lon/lat bounding box.

lon = E.lon.data;
lat = E.lat.data;

xpinp = lon(1,1);
ypinp = lat(1,1);
alpinp = 0; % defailt direction of positive x - axis 
mxinp = length(lon)-1;
myinp = length(lat)-1;
dxinp = lon(2)-lon(1);
dyinp = lat(2)-lat(1);
 
% shift time to the ROMS/SWAN basedate you want to use
tbeginp = E.time.data(1);
tendinp = E.time.data(end);
deltatinp = (E.time.data(2)-E.time.data(1))*24; % HR

% write the data ----------------------------------------

Time = E.time.data;

fid = fopen(SWAN_forc_name,'w');
for i=1:length(Time)

    disp(['Writing winds for SWAN at ',datestr(Time(i))])
    uswan = squeeze(E.u10.data(:,:,i))';
    vswan = squeeze(E.v10.data(:,:,i));
    fprintf(fid,'%10.2f\n',uswan');
    fprintf(fid,'%10.2f\n',vswan');

end
fclose(fid);

    disp(['xinp for SWAN input file  = ',num2str(xpinp)])
    disp(['yinp for SWAN input file  = ',num2str(ypinp)])
    disp(['alpinp for SWAN input file  = ',num2str(alpinp)])
    disp(['mxinp for SWAN input file  = ',num2str(mxinp)])
    disp(['myinp for SWAN input file  = ',num2str(myinp)])
    disp(['dxinp for SWAN input file  = ',num2str(dxinp)])
    disp(['dyinp for SWAN input file  = ',num2str(dyinp)])
    disp(['tbeginp for SWAN input file  = ',datestr(tbeginp)])
    disp(['deltatinp for SWAN input file  = ', num2str(deltatinp)])
    disp(['tendinp for SWAN input file  = ',datestr(tendinp)])

