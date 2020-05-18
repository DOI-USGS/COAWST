%
%  D_TS_METOFFICE:  Driver script to create a 4D-Var T,S observations file.
%
%  This a user modifiable script that can be used to prepare ROMS 4D-Var
%  SST observations NetCDF file. The T,S data is extracted from the UK
%  Met office Hadley Centre, quality controlled, EN3 observation datasets.
%  It uses the script 'load_ts_metoffice.m' to extract profiles of potential
%  temperature and salinity for the application grid and requested times.
%
%  The hydrographic data are available from 1950 to the present and stored
%  in annual tar files (said, EN3_v2a_Profiles_2004.tar) containing a NetCDF
%  for each month (said, EN3_v2a_Profiles_200401.nc.gz). Since this dataset
%  is not available in an OpenDAP server, the user needs to get the tar
%  files and extract the montly NetCDF files in it. The URL is:
%
%    http://hadobs.metoffice.com/en3/data/EN3_v2a/download_EN3_v2a.html
%
%  The script 'load_ts_metoffice.m' can processed compressed NetCDF files.
%  If compressed files are passed, the script will uncompress them using
%  "gunzip" and compress back after processing with "gzip". The user can
%  pass a single monthly file or a cell string with several monthly files.
%
%  USERS can use this script as a prototype for their application.
%

% svn $Id: d_ts_metoffice.m 996 2020-01-10 04:28:56Z arango $
%=========================================================================%
%  Copyright (c) 2002-2020 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%

%  Set input application Grid and history NetCDF files. The history file
%  is used to compute the depth of model vertical levels.

 my_root = '~/ocean/repository/test';

 GRDfile = strcat(my_root, '/WC13/Data/wc13_grd.nc'); 
 HISfile = strcat(my_root, '/WC13/Data/wc13_ini.nc'); 
 
%  Set input UK Met Office NetCDF file name(s).

 MET_dir = '~/ocean/repository/test/WC13/OBS/EN3';

 METfile = strcat(MET_dir, '/EN3_v2a_Profiles_200401.nc.gz');

 METcell = {strcat(MET_dir, '/EN3_v2a_Profiles_200401.nc.gz'),          ...
            strcat(MET_dir, '/EN3_v2a_Profiles_200402.nc.gz'),          ...
            strcat(MET_dir, '/EN3_v2a_Profiles_200403.nc.gz'),          ...
            strcat(MET_dir, '/EN3_v2a_Profiles_200404.nc.gz'),          ...
            strcat(MET_dir, '/EN3_v2a_Profiles_200405.nc.gz'),          ...
            strcat(MET_dir, '/EN3_v2a_Profiles_200406.nc.gz'),          ...
            strcat(MET_dir, '/EN3_v2a_Profiles_200407.nc.gz'),          ...
            strcat(MET_dir, '/EN3_v2a_Profiles_200408.nc.gz'),          ...
            strcat(MET_dir, '/EN3_v2a_Profiles_200409.nc.gz'),          ...
            strcat(MET_dir, '/EN3_v2a_Profiles_200410.nc.gz'),          ...
            strcat(MET_dir, '/EN3_v2a_Profiles_200411.nc.gz'),          ...
            strcat(MET_dir, '/EN3_v2a_Profiles_200412.nc.gz')};

%  Set observations output NetCDF files.
 
 OBSfile = 'wc13_ts_obs.nc';
 SUPfile = 'wc13_ts_super_obs.nc';

%  Set ROMS state variable type classification.

Nstate=7;            % number of ROMS state variables

state.zeta = 1;      % free-surface
state.ubar = 2;      % vertically integrated u-momentum
state.vbar = 3;      % vertically integrated v-momentum
state.u    = 4;      % u-momentum
state.v    = 5;      % v-momentum
state.temp = 6;      % temperature
state.salt = 7;      % salinity

%  Set observations provenance.  This is an arbitrary classification.
%  You need to choose unique values for instruments, type of measurement,
%  or data survey.

provenance.ssh_aviso = 1;    % AVISO SSH
provenance.sst_blend = 2;    % blended SST
provenance.Txbt_MetO = 3;    % XBT temperature from Met Office
provenance.Tctd_MetO = 4;    % CTD temperature from Met Office
provenance.Sctd_MetO = 5;    % CTD salinity from Met Office
provenance.Targo     = 6;    % ARGO floats temperature from Met Office
provenance.Sargo     = 7;    % ARGO floats salinity from Met Office
provenance.Tctd_CalC = 8;    % CTD temperature from CalCOFI
provenance.Sctd_CalC = 9;    % CTD salinity from CalCOFI
provenance.Tctd_GLOB = 10;   % CTD temperature from GLOBEC
provenance.Sctd_GLOB = 11;   % CTD salinity from GLOBEC
provenance.Tbuy_MetO = 12;   % Buoys, thermistor temperature form MetOffice

%  Set observations error by instrument.  Square values since we
%  need variances.

error.Targo = 0.1;        error.Targo = error.Targo ^2;
error.Tbuoy = 0.1;        error.Tbuoy = error.Tbuoy ^2;
error.Tctd  = 0.1;        error.Tctd  = error.Tctd  ^2;
error.Txbt  = 0.1;        error.Txbt  = error.Txbt  ^2;

error.Sargo = 0.01;       error.Sargo = error.Sargo ^2;
error.Sctd  = 0.01;       error.Sctd  = error.Sctd  ^2;

%  Set number of vertical levels in application. The depth of the
%  satellite data is assign to the surface level.

Nsur = length(nc_read(HISfile,'s_rho'));

%  Set 'grid_Lm_Mm_N' global attribute. It needs to be integer. We need
%  to have this attribute for data quality control (like binning) in other
%  programs.

[Lr,Mr]=size(nc_read(GRDfile,'h'));

C.grid_Lm_Mm_N = int32([Lr-2 Mr-2 Nsur]);

%  Set switch to apply small correction due to spherical/curvilinear
%  grids (see "obs_ijpos.m").

Correction = true;

%  Set switch to include observations on the applicationopen boundary
%  edge (see "obs_ijpos.m").

obc_edge = false;

%  Set ROMS application I- and J-grid offset to process observation away
%  from boundaries, if so desired.

Ioffset(1)=0;     % I-grid offset on the edge where Istr=1
Ioffset(2)=0;     % I-grid offset on the edge where Iend=Lm

Joffset(1)=0;     % J-grid offset on the edge where Jstr=1
Joffset(2)=0;     % J-grid offset on the edge where Jend=Lm

%---------------------------------------------------------------------------
%  Extract SST observations and store them into structure array D.
%---------------------------------------------------------------------------

%  Set spherical switch.

obs.spherical = 1;

%  Set time to process.

  StartDay = datenum(2004,1, 1);
  EndDay   = datenum(2004,1,15);

% [T,S] = load_ts_metoffice(METcell, GRDfile);
  [T,S] = load_ts_metoffice(METfile, GRDfile, StartDay, EndDay);

%  The script 'load_ts_metoffice.m' loads the extracted potential 
%  temperature and salinity in two different structures (T, S).
%  Notice that salinity data in available in XBT or thermistor
%  chain profilers. Each individual sctructure has 1D array for
%  each datum and the time has been sorted in increasing time
%  order.  We get:
%
%       T.time     S.time       time (date number)
%       T.lon      S.lon        longitude (degree_east)
%       T.lat      S.lat        latitude
%       T.depth    S.depth      depth (meter, negative)
%       T.WMO      S.WMO        WMO instrument type
%       T.value    S.value      observation values

T.type  = ones(size(T.value)) .* state.temp;
S.type  = ones(size(S.value)) .* state.salt;

T.error = zeros(size(T.value));
S.error = zeros(size(S.value));

T.provenance = zeros(size(T.value));
S.provenance = zeros(size(S.value));

ind = find(T.WMO == 401);
if (~isempty(ind)),
  T.provenance(ind) = provenance.Txbt_MetO;
  T.error(ind)      = error.Txbt;
end
ind = find(T.WMO == 741);
if (~isempty(ind)),
  T.provenance(ind) = provenance.Tctd_MetO;
  T.error(ind)      = error.Tctd;
end
ind = find(T.WMO == 820);
if (~isempty(ind)),
  T.provenance(ind) = provenance.Tbuy_MetO;
  T.error(ind)      = error.Tbuoy;
end
ind = find(T.WMO == 831);
if (~isempty(ind)),
  T.provenance(ind) = provenance.Targo;
  T.error(ind)      = error.Targo;
end

ind = find(S.WMO == 741);
if (~isempty(ind)),
  S.provenance(ind) = provenance.Sctd_MetO;
  S.error(ind)      = error.Sctd;
end
ind = find(S.WMO == 831);
if (~isempty(ind)),
  S.provenance(ind) = provenance.Sargo;
  S.error(ind)      = error.Sargo;
end

%  Build observation sctructure. Sort observations in time
%  ascending order.

obs.type       = [T.type,       S.type      ];
obs.provenance = [T.provenance, S.provenance];
obs.time       = [T.time,       S.time      ];
obs.lon        = [T.lon,        S.lon       ];
obs.lat        = [T.lat,        S.lat       ];
obs.depth      = [T.depth,      S.depth     ];
obs.error      = [T.error,      S.error     ];
obs.value      = [T.value,      S.value     ];

[Y,I] = sort(obs.time, 'ascend');

obs.type       = obs.type(I);
obs.provenance = obs.provenance(I);
obs.time       = obs.time(I);
obs.lon        = obs.lon(I);
obs.lat        = obs.lat(I);
obs.depth      = obs.depth(I);
obs.error      = obs.error(I);
obs.value      = obs.value(I);

%  Compute observation fractional grid coordinates in term
%  of ROMS grid.

[obs.Xgrid, obs.Ygrid] = obs_ijpos(GRDfile, obs.lon, obs.lat, ...
                                   Correction, obc_edge, ...
                                   Ioffset, Joffset);

obs.Xgrid = transpose(obs.Xgrid);  % take the transpose of the locations     
obs.Ygrid = transpose(obs.Ygrid);  % to faciliate structure expansion

ind = find(isnan(obs.Xgrid) & isnan(obs.Ygrid));
if (~isempty(ind));                % remove NaN's from data
  obs.type      (ind) = [];
  obs.provenance(ind) = [];
  obs.time      (ind) = [];
  obs.lon       (ind) = [];
  obs.lat       (ind) = [];
  obs.depth     (ind) = [];
  obs.Xgrid     (ind) = [];
  obs.Ygrid     (ind) = [];
  obs.error     (ind) = [];
  obs.value     (ind) = [];
end

%  Determine number of unique surveys times and number of
%  observation per survey.  They are already sorted in
%  increased time order.

obs.survey_time = unique(obs.time);
obs.Nsurvey     = length(obs.survey_time);
obs.Ndatum      = length(obs.value);

for n=1:obs.Nsurvey,
  ind = find(obs.time == obs.survey_time(n));
  obs.Nobs(n) = length(ind);
end
obs.Nobs = obs.Nobs;

%  Set fractional z-grid coordinates for the observations.

zflag = 0;                          % use zero free-surface

[obs] = obs_depth(HISfile, obs, zflag);

%  Initialize global variance per state variables.

obs.variance = zeros([1 Nstate]);

obs.variance(state.temp) = error.Targo;
obs.variance(state.salt) = error.Sargo;

%--------------------------------------------------------------------------
%  Set observation file creation parameter in structure array, C.
%--------------------------------------------------------------------------

%  Observations output file name.

C.ncfile = OBSfile;

%  Application grid NetCDF file name.

C.grd_file = GRDfile;

%  Set application title.

C.title = 'California Current System, 1/3 degree resolution (WC13)';

%  Spherical grid switch.
%
%            [0] Cartesian grid
%            [1] Spherical grid

C.spherical = 1;

%  Set switches to include the 'obs_lon' and 'obs_lat' variables.
%  This are not used inside ROMS but are needed during pre- and
%  post-processing.

C.do_longitude = 1;
C.do_latitude  = 1;

%  Number of state variables. Usually, ROMS uses 7 variables in the
%  observation state vector (zeta, ubar, vbar, u, v, temperature,
%  and salinity). This number need to be increased when considering
%  additional tracer variables.  This is the value of the NetCDF
%  dimension 'state_variable'.

C.Nstate = Nstate;

%  Number of data surveys. That is, number of unique survey times
%  in the observational dataset. This is the value of the NetCDF
%  dimension 'survey'.

C.Nsurvey = obs.Nsurvey;

%  Total number of observations in space and time. This is the value
%  of the NetCDF dimension 'datum'.  If zero, an unlimited dimension
%  is used.

%C.Ndatum = 0;             % use unlimited datum record dimension
 C.Ndatum = obs.Ndatum;    % fixed datum dimension
 
%  Set attributes for 'obs_type' variable which assigns the model
%  state variable associated with the observation. Usually, ROMS
%  uses 7 variables in the observation state vector:
%
%      obs_type = 1    free-surface
%      obs_type = 2    vertically integrated u-momentum component
%      obs_type = 3    vertically integrated v-momentum component
%      obs_type = 4    u-momentum component
%      obs_type = 5    v-momentum component
%      obs_type = 6    potential temperature
%      obs_type = 7    salinity
%      obs_type = ...  other passive tracers NAT+1:NT
%  
%  NOTE: We are following the CF compliance rules for variable
%        attributes 'flag_values' and 'flag_meanings'.

C.state_flag_values = 1:1:7;

C.state_flag_meanings = ['zeta', blanks(1),                             ...
                        'ubar', blanks(1),                              ...
                        'vbar', blanks(1),                              ...
                        'u', blanks(1),                                 ...
                        'v', blanks(1),                                 ...
                        'temperature', blanks(1),                       ...
                        'salinity'];

%  Set attributes for 'obs_provenance' variable which assigns different
%  flags for each instrument or data source.  This information is used
%  in the observation impact and observation sensitivity analysis. The
%  user has a lot of latitute here.
%
%  NOTE: We are following the CF compliance rules for variable
%        attributes 'flag_values' and 'flag_meanings'. All the
%        blank spaces for each meaning needs to be separate with
%        underscores.


C.origin_flag_values = 1:1:12;

C.origin_flag_meanings = ['gridded_AVISO_SLA', blanks(1),               ...
                         'blended_SST', blanks(1),                      ...
                         'XBT_Met_Office', blanks(1),                   ...
                         'CTD_temperature_Met_Office', blanks(1),       ...
                         'CTD_salinity_Met_Office', blanks(1),          ...
                         'ARGO_temperature_Met_Office', blanks(1),      ...
                         'ARGO_salinity_Met_Office', blanks(1),         ...
                         'CTD_temperature_CalCOFI', blanks(1),          ...
                         'CTD_salinity_CalCOFI', blanks(1),             ...
                         'CTD_temperature_GLOBEC', blanks(1),           ...
                         'CTD_salinity_GLOBEC', blanks(1),              ...
                         'buoy_temperature_Met_Office'];

%  The attribute association between 'flag_values' and 'flag_meanings'
%  is difficult to read when a moderate number of flags are use. To
%  aid the decoding, two readable global attributes are added:
%  'state_variables' and 'obs_provenance' which are stored in
%  S.global_variables and S.global_provenance, respectively.
%
%  NOTE: the 'state_variables' attribute include the units.

newline=sprintf('\n');

C.global_variables=[newline, ...
       '1: free-surface (m) ', newline, ...
       '2: vertically integrated u-momentum component (m/s) ', newline, ...
       '3: vertically integrated v-momentum component (m/s) ', newline, ...
       '4: u-momentum component (m/s) ', newline,                       ...
       '5: v-momentum component (m/s) ', newline,                       ...
       '6: potential temperature (Celsius) ', newline,                  ...
       '7: salinity (nondimensional)'];

C.global_provenance=[newline, ...
       ' 1: gridded AVISO sea level anomaly ', newline,                 ...
       ' 2: blended satellite SST ', newline,                           ...
       ' 3: XBT temperature from Met Office ', newline,                 ...
       ' 4: CTD temperature from Met Office ', newline,                 ...
       ' 5: CTD salinity from Met Office ', newline,                    ...
       ' 6: ARGO floats temperature from Met Office ', newline,         ...
       ' 7: ARGO floats salinity from Met Office ', newline,            ...
       ' 8: CTD temperature from CalCOFI ', newline,                    ...
       ' 9: CTD salinity from CalCOFI ', newline,                       ...
       '10: CTD temperature from GLOBEC ', newline,                     ...
       '11: CTD salinity from GLOBEC ', newline,                        ...
       '12: buoy, thermistor temperature from Met Office'];

%  Set the observation data sources global attribute 'obs_sources'
%  which is stored in S.global_sources (OPTIONAL).

C.global_sources=[newline, ...
       'http://opendap.aviso.oceanobs.com/thredds/dodsC ', newline, ...
       'http://thredds1.pfeg.noaa.gov:8080/thredds/dodsC/satellite/BA/' ...
             'ssta/5day ', newline, ...
       'http://hadobs.metoffice.com/en3'];

%--------------------------------------------------------------------------
%  Create 4D-Var observations NetCDF file.
%--------------------------------------------------------------------------

[status]=c_observations(C);

%  Update 'units' attribute for time variables. Notice that the time
%  of the observations in the NetCDF file is in DAYS.

avalue='days since 1968-05-23 00:00:00 GMT';

[status]=nc_attadd(OBSfile,'units',avalue,'survey_time');
[status]=nc_attadd(OBSfile,'calendar','gregorian','survey_time');

[status]=nc_attadd(OBSfile,'units',avalue,'obs_time');
[status]=nc_attadd(OBSfile,'calendar','gregorian','obs_time');

%--------------------------------------------------------------------------
%  Write 4D-Var observations NetCDF file.
%--------------------------------------------------------------------------

[status]=obs_write(OBSfile, obs);

%--------------------------------------------------------------------------
%  Super observations.
%--------------------------------------------------------------------------

%  Add needed fields to "obs" structure for the case that it is used latter
%  when computing super observations. We can either use the structure of
%  the NetCDF file OBSfile created above. This allows flexibility.

obs.grid_Lm_Mm_N = C.grid_Lm_Mm_N;
obs.ncfile       = OBSfile;

%  It is possible that more that one observations associatiated to the
%  same ROMS state variable is available at the same time in a particular
%  grid cell. If this is the case, we need average the data within that
%  cell and create a super observation. The following script just do
%  that and saves everything in observation structure OBS.

[OBS]=super_obs(OBSfile);

%  Meld inital observation error with the values computed when binning
%  the data into super observations. Take the larger value.

OBS.error = sqrt(OBS.error.^2 + OBS.std.^2);

%  Write new structure to a new NetCDF.

[status]=c_observations(OBS, SUPfile);

avalue='days since 1968-05-23 00:00:00 GMT';

[status]=nc_attadd(SUPfile,'units',avalue,'survey_time');
[status]=nc_attadd(SUPfile,'calendar','gregorian','survey_time');

[status]=nc_attadd(SUPfile,'units',avalue,'obs_time');
[status]=nc_attadd(SUPfile,'calendar','gregorian','obs_time');

[status]=obs_write(SUPfile, OBS);

disp(' ');
disp('Done.');
disp(' ');
