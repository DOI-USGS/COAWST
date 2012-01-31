function roms_grid_create_netcdf(theRomsFile,LP,MP)
 
% Copyright (C) 1999 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
 
% Version of 18-Jun-1999 17:11:59.
% Updated    23-Oct-2000 11:30:58.

RADIAN_CONVERSION_FACTOR = 180/pi;
EARTH_RADIUS_METERS = 6378*1000;   % Equatorial radius.

% With grid_x of size [m, n], the grid itself has
%  [m-1, n-1] cells.  The latter size corresponds
%  to the size of the mask and bathymetry.  These
%  cell-centers are called the "rho" points.

% Create the Roms NetCDF file.

%% ncdump('priapus:MATLAB:toolbox:Roms:tasman_grd.nc')   %% Generated 17-Jun-1999 17:37:30
 
nc = netcdf(theRomsFile, 'clobber');
if isempty(nc), return, end
 
%% Global attributes:

disp(' ## Defining Global Attributes...')
 
nc.type = ncchar('Gridpak file');
nc.gridid = 'combined grid';
nc.history = ncchar(['Created by "' mfilename '" on ' datestr(now)]);

nc.CPP_options = ncchar('DCOMPLEX, DBLEPREC, NCARG_32, PLOTS,');
name(nc.CPP_options, 'CPP-options')

%LP = (m-1)/2;   % The rho dimension.
L = LP-1;       % The psi dimension.

% The eta direction (up-down):

%MP = (n-1)/2;   % The rho dimension.
M = MP-1;       % The psi dimension.

disp(' ## Defining Dimensions...')
 
nc('xi_psi') = L;
nc('xi_rho') = LP;
nc('xi_u') = L;
nc('xi_v') = LP;

nc('eta_psi') = M;
nc('eta_rho') = MP;
nc('eta_u') = MP;
nc('eta_v') = M;

nc('two') = 2;
nc('bath') = 0; %% (record dimension)
 
%% Variables and attributes:

disp(' ## Defining Variables and Attributes...')
 
nc{'xl'} = ncdouble; %% 1 element.
nc{'xl'}.long_name = ncchar('domain length in the XI-direction');
nc{'xl'}.units = ncchar('meter');
 
nc{'el'} = ncdouble; %% 1 element.
nc{'el'}.long_name = ncchar('domain length in the ETA-direction');
nc{'el'}.units = ncchar('meter');
 
nc{'JPRJ'} = ncchar('two'); %% 2 elements.
nc{'JPRJ'}.long_name = ncchar('Map projection type');

nc{'JPRJ'}.option_ME_ = ncchar('Mercator');
nc{'JPRJ'}.option_ST_ = ncchar('Stereographic');
nc{'JPRJ'}.option_LC_ = ncchar('Lambert conformal conic');
name(nc{'JPRJ'}.option_ME_, 'option(ME)')
name(nc{'JPRJ'}.option_ST_, 'option(ST)')
name(nc{'JPRJ'}.option_LC_, 'option(LC)')
 
nc{'PLAT'} = ncfloat('two'); %% 2 elements.
nc{'PLAT'}.long_name = ncchar('Reference latitude(s) for map projection');
nc{'PLAT'}.units = ncchar('degree_north');
 
nc{'PLONG'} = ncfloat; %% 1 element.
nc{'PLONG'}.long_name = ncchar('Reference longitude for map projection');
nc{'PLONG'}.units = ncchar('degree_east');
 
nc{'ROTA'} = ncfloat; %% 1 element.
nc{'ROTA'}.long_name = ncchar('Rotation angle for map projection');
nc{'ROTA'}.units = ncchar('degree');
 
nc{'JLTS'} = ncchar('two'); %% 2 elements.
nc{'JLTS'}.long_name = ncchar('How limits of map are chosen');

nc{'JLTS'}.option_CO_ = ncchar('P1, .. P4 define two opposite corners ');
nc{'JLTS'}.option_MA_ = ncchar('Maximum (whole world)');
nc{'JLTS'}.option_AN_ = ncchar('Angles - P1..P4 define angles to edge of domain');
nc{'JLTS'}.option_LI_ = ncchar('Limits - P1..P4 define limits in u,v space');
name(nc{'JLTS'}.option_CO_, 'option(CO)')
name(nc{'JLTS'}.option_MA_, 'option(MA)')
name(nc{'JLTS'}.option_AN_, 'option(AN)')
name(nc{'JLTS'}.option_LI_, 'option(LI)')
 
nc{'P1'} = ncfloat; %% 1 element.
nc{'P1'}.long_name = ncchar('Map limit parameter number 1');
 
nc{'P2'} = ncfloat; %% 1 element.
nc{'P2'}.long_name = ncchar('Map limit parameter number 2');
 
nc{'P3'} = ncfloat; %% 1 element.
nc{'P3'}.long_name = ncchar('Map limit parameter number 3');
 
nc{'P4'} = ncfloat; %% 1 element.
nc{'P4'}.long_name = ncchar('Map limit parameter number 4');
 
nc{'XOFF'} = ncfloat; %% 1 element.
nc{'XOFF'}.long_name = ncchar('Offset in x direction');
nc{'XOFF'}.units = ncchar('meter');
 
nc{'YOFF'} = ncfloat; %% 1 element.
nc{'YOFF'}.long_name = ncchar('Offset in y direction');
nc{'YOFF'}.units = ncchar('meter');
 
nc{'depthmin'} = ncshort; %% 1 element.
nc{'depthmin'}.long_name = ncchar('Shallow bathymetry clipping depth');
nc{'depthmin'}.units = ncchar('meter');
 
nc{'depthmax'} = ncshort; %% 1 element.
nc{'depthmax'}.long_name = ncchar('Deep bathymetry clipping depth');
nc{'depthmax'}.units = ncchar('meter');
 
nc{'spherical'} = ncchar; %% 1 element.
nc{'spherical'}.long_name = ncchar('Grid type logical switch');
nc{'spherical'}.option_T_ = ncchar('spherical');
nc{'spherical'}.option_F_ = ncchar('Cartesian');
name(nc{'spherical'}.option_T_, 'option(T)')
name(nc{'spherical'}.option_F_, 'option(F)')
 
nc{'hraw'} = ncdouble('bath', 'eta_rho', 'xi_rho'); %% 0 elements.
nc{'hraw'}.long_name = ncchar('Working bathymetry at RHO-points');
nc{'hraw'}.units = ncchar('meter');
nc{'hraw'}.field = ncchar('bath, scalar');
 
nc{'h'} = ncdouble('eta_rho', 'xi_rho'); %% 16900 elements.
nc{'h'}.long_name = ncchar('Final bathymetry at RHO-points');
nc{'h'}.units = ncchar('meter');
nc{'h'}.field = ncchar('bath, scalar');
 
nc{'f'} = ncdouble('eta_rho', 'xi_rho'); %% 16900 elements.
nc{'f'}.long_name = ncchar('Coriolis parameter at RHO-points');
nc{'f'}.units = ncchar('second-1');
nc{'f'}.field = ncchar('Coriolis, scalar');
 
nc{'pm'} = ncdouble('eta_rho', 'xi_rho'); %% 16900 elements.
nc{'pm'}.long_name = ncchar('curvilinear coordinate metric in XI');
nc{'pm'}.units = ncchar('meter-1');
nc{'pm'}.field = ncchar('pm, scalar');
 
nc{'pn'} = ncdouble('eta_rho', 'xi_rho'); %% 16900 elements.
nc{'pn'}.long_name = ncchar('curvilinear coordinate metric in ETA');
nc{'pn'}.units = ncchar('meter-1');
nc{'pn'}.field = ncchar('pn, scalar');
 
nc{'dndx'} = ncdouble('eta_rho', 'xi_rho'); %% 16900 elements.
nc{'dndx'}.long_name = ncchar('xi derivative of inverse metric factor pn');
nc{'dndx'}.units = ncchar('meter');
nc{'dndx'}.field = ncchar('dndx, scalar');
 
nc{'dmde'} = ncdouble('eta_rho', 'xi_rho'); %% 16900 elements.
nc{'dmde'}.long_name = ncchar('eta derivative of inverse metric factor pm');
nc{'dmde'}.units = ncchar('meter');
nc{'dmde'}.field = ncchar('dmde, scalar');
 
nc{'x_rho'} = ncdouble('eta_rho', 'xi_rho'); %% 16900 elements.
nc{'x_rho'}.long_name = ncchar('x location of RHO-points');
nc{'x_rho'}.units = ncchar('meter');
 
nc{'y_rho'} = ncdouble('eta_rho', 'xi_rho'); %% 16900 elements.
nc{'y_rho'}.long_name = ncchar('y location of RHO-points');
nc{'y_rho'}.units = ncchar('meter');
 
nc{'x_psi'} = ncdouble('eta_psi', 'xi_psi'); %% 16641 elements.
nc{'x_psi'}.long_name = ncchar('x location of PSI-points');
nc{'x_psi'}.units = ncchar('meter');
 
nc{'y_psi'} = ncdouble('eta_psi', 'xi_psi'); %% 16641 elements.
nc{'y_psi'}.long_name = ncchar('y location of PSI-points');
nc{'y_psi'}.units = ncchar('meter');
 
nc{'x_u'} = ncdouble('eta_u', 'xi_u'); %% 16770 elements.
nc{'x_u'}.long_name = ncchar('x location of U-points');
nc{'x_u'}.units = ncchar('meter');
 
nc{'y_u'} = ncdouble('eta_u', 'xi_u'); %% 16770 elements.
nc{'y_u'}.long_name = ncchar('y location of U-points');
nc{'y_u'}.units = ncchar('meter');
 
nc{'x_v'} = ncdouble('eta_v', 'xi_v'); %% 16770 elements.
nc{'x_v'}.long_name = ncchar('x location of V-points');
nc{'x_v'}.units = ncchar('meter');
 
nc{'y_v'} = ncdouble('eta_v', 'xi_v'); %% 16770 elements.
nc{'y_v'}.long_name = ncchar('y location of V-points');
nc{'y_v'}.units = ncchar('meter');
 
nc{'lat_rho'} = ncdouble('eta_rho', 'xi_rho'); %% 16900 elements.
nc{'lat_rho'}.long_name = ncchar('latitude of RHO-points');
nc{'lat_rho'}.units = ncchar('degree_north');
 
nc{'lon_rho'} = ncdouble('eta_rho', 'xi_rho'); %% 16900 elements.
nc{'lon_rho'}.long_name = ncchar('longitude of RHO-points');
nc{'lon_rho'}.units = ncchar('degree_east');
 
nc{'lat_psi'} = ncdouble('eta_psi', 'xi_psi'); %% 16641 elements.
nc{'lat_psi'}.long_name = ncchar('latitude of PSI-points');
nc{'lat_psi'}.units = ncchar('degree_north');
 
nc{'lon_psi'} = ncdouble('eta_psi', 'xi_psi'); %% 16641 elements.
nc{'lon_psi'}.long_name = ncchar('longitude of PSI-points');
nc{'lon_psi'}.units = ncchar('degree_east');
 
nc{'lat_u'} = ncdouble('eta_u', 'xi_u'); %% 16770 elements.
nc{'lat_u'}.long_name = ncchar('latitude of U-points');
nc{'lat_u'}.units = ncchar('degree_north');
 
nc{'lon_u'} = ncdouble('eta_u', 'xi_u'); %% 16770 elements.
nc{'lon_u'}.long_name = ncchar('longitude of U-points');
nc{'lon_u'}.units = ncchar('degree_east');
 
nc{'lat_v'} = ncdouble('eta_v', 'xi_v'); %% 16770 elements.
nc{'lat_v'}.long_name = ncchar('latitude of V-points');
nc{'lat_v'}.units = ncchar('degree_north');
 
nc{'lon_v'} = ncdouble('eta_v', 'xi_v'); %% 16770 elements.
nc{'lon_v'}.long_name = ncchar('longitude of V-points');
nc{'lon_v'}.units = ncchar('degree_east');
 
nc{'mask_rho'} = ncdouble('eta_rho', 'xi_rho'); %% 16900 elements.
nc{'mask_rho'}.long_name = ncchar('mask on RHO-points');
nc{'mask_rho'}.option_0_ = ncchar('land');
nc{'mask_rho'}.option_1_ = ncchar('water');
name(nc{'mask_rho'}.option_0_, 'option(0)')
name(nc{'mask_rho'}.option_1_, 'option(1)')
 
nc{'mask_u'} = ncdouble('eta_u', 'xi_u'); %% 16770 elements.
nc{'mask_u'}.long_name = ncchar('mask on U-points');
nc{'mask_u'}.option_0_ = ncchar('land');
nc{'mask_u'}.option_1_ = ncchar('water');
%name(nc{'mask_u'}.option_0_, 'option(0)')
%name(nc{'mask_u'}.option_1_, 'option(1)')
%		 nc{'mask_u'}.FillValue_ = ncdouble(1);
 
nc{'mask_v'} = ncdouble('eta_v', 'xi_v'); %% 16770 elements.
nc{'mask_v'}.long_name = ncchar('mask on V-points');
nc{'mask_v'}.option_0_ = ncchar('land');
nc{'mask_v'}.option_1_ = ncchar('water');
%name(nc{'mask_v'}.option_0_, 'option(0)')
%name(nc{'mask_v'}.option_1_, 'option(1)')
%		 nc{'mask_v'}.FillValue_ = ncdouble(1);
 
nc{'mask_psi'} = ncdouble('eta_psi', 'xi_psi'); %% 16641 elements.
nc{'mask_psi'}.long_name = ncchar('mask on PSI-points');
nc{'mask_psi'}.option_0_ = ncchar('land');
nc{'mask_psi'}.option_1_ = ncchar('water');
%name(nc{'mask_psi'}.option_0_, 'option(0)')
%name(nc{'mask_psi'}.option_1_, 'option(1)')
%		 nc{'mask_psi'}.FillValue_ = ncdouble(1);

% Now, what about depths: "h" and "hraw".  <== DEPTHS.

% The following seems mistaken:
 
nc{'angle'} = ncdouble('eta_rho', 'xi_rho'); %% 16900 elements.
nc{'angle'}.long_name = ncchar('angle between xi axis and east');
nc{'angle'}.units = ncchar('degree');

close(nc)
