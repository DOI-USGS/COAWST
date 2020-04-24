%
% ROMS 4D-Var data assimilation Matlab scripts
% ===========================================
%
% This package contains several generic Matlab scripts to process
% data for ROMS 4D-Var data assimilation algorithms.
%
% 4D-Var Observations:
%
%   c_observations    - Creates 4D-Var observation NetCDF file.
%   d_observations    - Driver to process 4D-Var observation NetCDF file.
%
%
%   obs_extract       - Extracts and creates observation NetCDF files
%                       from input data assimilation observation NetCDF
%                       files at the requested interval.
%   obs_read          - Reads observation NetCDF file and load all data
%                         into a structure array.
%   obs_write         - Writes all observation data in structure array into
%                         an existing NetCDF file.
%
%   obs_depth         - Computes the fraction z-grid locations of the 
%                         observations.
%   obs_merge         - Merges data from several 4D-Var observations NetCDF
%                         files.
%
%   obs_ijpos         - Computes observations locations in ROMS fractional
%                         coordinates. It uses the 'inpolygon' intrinsic
%                         Matlab function to check if the observations are
%                         inside or at the egde of the application grid.
%
%   super_obs         - Checks the provided observation data (4D-Var NetCDF
%                         file or structure) and creates super observations if
%                         there is more than one datum of the same observation
%                         type per grid cell.
%   plot_super        - Sample script to compute and plot super observations
%                         from specified 4D-Var observations NetCDF file.
%
%   d_ssh_obs         - Driver template to extract SSH observations for
%                         AVISO. It creates and writes observation file.
%                         Then, it computes super observations and creates
%                         and writes super observations NetCDF file.
%   d_sst_merge       - Driver template to merge several SST ROMS 4D-Var
%                         observations NetCDF files. It allows the user
%                         to merge data from several satellites and
%                         process super observations.
%   d_sst_obs         - Driver template to extract SST observations from
%                         satellite data. It creates and writes observation
%                         file. Then, it computes super observations and
%                         creates and writes super observations NetCDF file.
%   d_ts_metoffice    - Driver template to extract hydrographic potential
%                         temperature and salinity profile data from the
%                         UK Met Office observation datasets.
%
%   load_ssh_aviso    - Extracts AVISO sea level anomaly for the period of
%                         interest and specified region from ROMS Grid file.
%   load_sst_AMSRE    - Extracts satellite sea surface temperature for the
%                         period of interest and specified region from
%                         ROMS Grid file. The SST data is from the OpenDAP
%                         catalog maintained by NOAA PFEG CoastWatch in
%                         California. The resolution is 0.025 degree global
%                         1-day average product.
%   load_sst_pfeg     - Extracts satellite sea surface temperature for the
%                         period of interest and specified region from
%                         ROMS Grid file. The SST data is from the OpenDAP
%                         catalog maintained by NOAA PFEG CoastWatch in
%                         California. The resolution is 0.1 degree global
%                         5-day average composite.
%   load_ts_metoffice - Extracts potential temperature and salinity profiles
%                         from the UK Met Office, quality controlled, EN3
%                         observation datasets. The data are available from
%                         1950 to the present and include XBT, CTD, buoys,
%                         thermistor chain, and ARGO floats.
%
% Error Covariance Matrix:
%
%   average:          - Computes the time average of requested NetCDF
%                         variable.
%   variance:         - Computes the variance of requested NetCDF variable
%                         from its specified time mean.
%
% Error Covariance Matrix Balance Operator:
%
%   balance_4dvar     - Computes 4D-Var balance operator.
%   biconj            - Biconjugate gradient solver for the SSH elliptic
%                         equation.
%   ini_balance       - Initializes balance operator structure array.
%                         It sets internal parameters, reads needed grid
%                         metrics and computes several quantities.
%   lateral_obc       - Sets lateral boundary conditions for a 2D or 3D
%                         field.
%   rho_balance       - Computes balanced density anomaly using a linear
%                         equation of state.
%   s_balance         - Given a temperature anomaly, deltaT=T-Tavg, it 
%                         computes balanced salinity anomaly using a T-S
%                         empirical formula.
%   ssh_reference     - Computes the balance operator reference sea surface
%                         height.
%   uv_balance        - Computes balanced, baroclinic U- and V-momentum
%                         anomalies (m/s) using the geostrophic balance.
%   zeta_balance      - Computes balanced, baroclinic free-surface anomaly
%                         by solving an elliptical equation OR integrating
%                         the hydrostatic equation from surface to bottom.
%
% Error Covariance standard deviation:
%
%   c_std             - Creates 4D-Var initial conditions or model error
%                         standard deviation NetCDF file.
%   c_std_bry         - Creates 4D-Var open boundary conditions standard
%                         deviations NetCDF file.
%   c_std_frc         - Creates 4D-Var surface forcing standard deviations
%                         NetCDF file.
%   d_std             - Driver template to compute and write 4D-Var
%                         standard deviations for initial conditions or
%                         model error.
%   d_std_bry         - Driver template to extract and write open boundary
%                         standard deviations. The data is extracted from
%                         the initial standard deviation file.
%   d_std_frc         - Driver template to compute and write 4D-Var
%                         standard deviations for surface forcing.
%   d_std_unbalanced  - Driver template to compute and write 4D-Var
%                         unbalanced standard deviations for initial
%                         conditions.
%

% svn $Id: Contents.m 996 2020-01-10 04:28:56Z arango $
%===========================================================================%
%  Copyright (c) 2002-2020 The ROMS/TOMS Group                              %
%    Licensed under a MIT/X style license                                   %
%    See License_ROMS.txt                           Hernan G. Arango        %
%===========================================================================%
