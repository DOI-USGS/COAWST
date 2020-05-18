%
%  D_SST_OBS:  Driver script to create a 4D-Var SST observations file.
%
%  This a user modifiable script that can be used to prepare ROMS 4D-Var
%  SST observations NetCDF file. The SST data is extracted from the
%  extensive OpenDAP catalog maintained by NOAA PFEG Coastwatch Global 
%  using script 'load_sst_pfeg.m'.  USERS can use this as a prototype
%  for their application.
%

% svn $Id: d_sst_obs.m 996 2020-01-10 04:28:56Z arango $
%=========================================================================%
%  Copyright (c) 2002-2020 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%

%  Set SST product to use from CoastWatch OpenDAP server:

CoastWatch = 'http://oceanwatch.pfeg.noaa.gov/thredds/dodsC/satellite';

%  sst_URL = strcat(CoastWatch, '/AA/ssta/1day');    % SST, Aqua AMSR-E,
%  sst_URL = strcat(CoastWatch, '/AA/ssta/3day');    % Global
%  sst_URL = strcat(CoastWatch, '/AA/ssta/5day');
%  sst_URL = strcat(CoastWatch, '/AA/ssta/8day');
%  sst_URL = strcat(CoastWatch, '/AA/ssta/14day');
%  sst_URL = strcat(CoastWatch, '/AA/ssta/mday');

   sst_URL = strcat(CoastWatch, '/BA/ssta/5day');    % SST, Blended, Global
%  sst_URL = strcat(CoastWatch, '/BA/ssta/8day');    % Experimental
%  sst_URL = strcat(CoastWatch, '/BA/ssta/mday');    % Product
   
%  Set input/output NetCDF files.

 my_root = '~/ocean/repository/test';

 GRDfile = strcat(my_root, '/WC13/Data/wc13_grd.nc');
 OBSfile = 'wc13_sst_obs.nc';
 SUPfile = 'wc13_sst_super_obs.nc';

% Base date for ROMS observation files: "days since 1968-05-23 00:00:00".
% (The WC13 application has a modified Julian day number as reference
%  time: May-23-1968).

mybasedate = datenum(1968,05,23,0,0,0);

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
provenance.Targo     = 7;    % ARGO floats salinity from Met Office
provenance.Tctd_CalC = 8;    % CTD temperature from CalCOFI
provenance.Sctd_CalC = 9;    % CTD salinity from CalCOFI
provenance.Tctd_GLOB = 10;   % CTD temperature from GLOBEC
provenance.Sctd_GLOB = 11;   % CTD salinity from GLOBEC
provenance.Tbuy_MetO = 12;   % Buoys, thermistor temperature form MetOffice

%  Set data error to 0.4 degrees Celsius.  Square values since we
%  need variances.

error.zeta = 0;           error.zeta = error.zeta ^2;
error.ubar = 0;           error.ubar = error.ubar ^2;
error.vbar = 0;           error.vbar = error.ubar ^2;
error.u    = 0;           error.u    = error.u ^2;
error.v    = 0;           error.v    = error.u ^2;
error.temp = 0.4;         error.temp = error.temp ^2;
error.salt = 0;           error.salt = error.salt ^2;

%  Set number of vertical levels in application. The depth of the
%  satellite data is assign to the surface level.

Nsur=30;

%  Set 'grid_Lm_Mm_N' global attribute. It needs to be integer. We need
%  to have this attribute for data quality control (like binning) in other
%  programs.

[Lr,Mr]=size(nc_read(GRDfile,'h'));

S.grid_Lm_Mm_N = int32([Lr-2 Mr-2 Nsur]);

%  Set switch to apply small correction due to spherical/curvilinear
%  grids (see "obs_ijpos.m").

Correction = true;

%  Set switch to include observations on the application open boundary
%  edge (see "obs_ijpos.m").

obc_edge = false;

%  Set ROMS application I- and J-grid offset to process observation away
%  from boundaries, if so desired.

Ioffset(1)=0;     % I-grid offset on the edge where Istr=1
Ioffset(2)=0;     % I-grid offset on the edge where Iend=Lm

Joffset(1)=0;     % J-grid offset on the edge where Jstr=1
Joffset(2)=0;     % J-grid offset on the edge where Jend=Lm

%--------------------------------------------------------------------------
%  Extract SST observations and store them into structure array D.
%--------------------------------------------------------------------------

%  Initialize "obs" structure.

obs = [];

%  Set spherical switch.

obs.spherical = 1;

%  The 'load_sst_pfeg' stores SST data as:   D.sst, D.time, D.lon, D.lat.

StartDay = datenum(2004,1, 1);
EndDay   = datenum(2004,1,15);

D = load_sst_pfeg(GRDfile, StartDay, EndDay, sst_URL);

%  Convert time to specified base date. The SST data from 'load_sst_pfeg'
%  has a time in day numbers. 

D.time = D.time - mybasedate;

%  Convert data to one-dimenension array and replicate the data to
%  the same dimension of D.sst. Notice that we get the following
%  dimensions after running 'load_sst_pfeg':
%
%       D.time(time)           already sorted in increased time order
%       D.lon (lon)
%       D.lat (lat)
%       D.sst (time,lat,lon)

if (length(D.time) == 1),
  D.sst = reshape(D.sst, [1, size(D.sst)]);  % add singleton dimension
end                                          % for time

[it,Jm,Im] = size(D.sst);

obs.time = repmat(D.time,[1 Jm Im]);
obs.lon  = permute(repmat(repmat(D.lon',[Jm 1]),[1 1 it]),[3 1 2]);
obs.lat  = permute(repmat(repmat(D.lat ,[1 Im]),[1 1 it]),[3 1 2]);

obs.time  = obs.time(:);           % no time sorting is necessary
obs.lon   = obs.lon(:);
obs.lat   = obs.lat(:);
obs.value = D.sst(:);

ind = find(isnan(obs.value));
if (~isempty(ind));                % remove NaN's from data, if any
  obs.time (ind) = [];
  obs.lon  (ind) = [];
  obs.lat  (ind) = [];
  obs.value(ind) = [];
end

%  Compute observation fractional grid coordinates in term
%  of ROMS grid.

[obs.Xgrid, obs.Ygrid] = obs_ijpos(GRDfile, obs.lon, obs.lat,           ...
                                   Correction, obc_edge,                ...
                                   Ioffset, Joffset);

ind = find(isnan(obs.Xgrid) & isnan(obs.Ygrid));
if (~isempty(ind));                % remove NaN's from data
  obs.time (ind) = [];
  obs.lon  (ind) = [];
  obs.lat  (ind) = [];
  obs.Xgrid(ind) = [];
  obs.Ygrid(ind) = [];
  obs.value(ind) = [];
end

%  Assign ROMS associated state variable and provenance.

obs.type = ones(size(obs.value)) .* state.temp;

obs.provenance = ones(size(obs.value)) .* provenance.sst_blend;

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

%  Set depths and fractional z-grid coordinates for the observations.
%  The SST data is a the surface.  Therefore, we need to assign
%  'obs_depth' to surface level and 'obs_Zgrid' to zero.  ROMS will
%  compute the value of 'obs_zgrid' when relevant.

obs.depth = ones(size(obs.value)) .* Nsur;

obs.Zgrid = zeros(size(obs.value));

%  Set observation error covariance (squared units).

obs.error = ones(size(obs.value)) .* error.temp;

%  Initialize global variance per state variables.

obs.variance = zeros([1 Nstate]);

obs.variance(state.temp) = error.temp;

%--------------------------------------------------------------------------
%  Set observation file creation parameter in structure array, S.
%--------------------------------------------------------------------------

%  Observations output file name.

S.ncfile = OBSfile;

%  Application grid NetCDF file name.

S.grd_file = GRDfile;

%  Set application title.

S.title = 'California Current System, 1/3 degree resolution (WC13)';

%  Spherical grid switch.
%
%            [0] Cartesian grid
%            [1] Spherical grid

S.spherical = 1;

%  Set switches to include the 'obs_lon' and 'obs_lat' variables.
%  This are not used inside ROMS but are needed during pre- and
%  post-processing.

S.do_longitude = 1;
S.do_latitude  = 1;

%  Number of state variables. Usually, ROMS uses 7 variables in the
%  observation state vector (zeta, ubar, vbar, u, v, temperature,
%  and salinity). This number need to be increased when considering
%  additional tracer variables.  This is the value of the NetCDF
%  dimension 'state_variable'.

S.Nstate = Nstate;

%  Number of data surveys. That is, number of unique survey times
%  in the observational dataset. This is the value of the NetCDF
%  dimension 'survey'.

S.Nsurvey = obs.Nsurvey;

%  Total number of observations in space and time. This is the value
%  of the NetCDF dimension 'datum'.  If zero, an unlimited dimension
%  is used.

%S.Ndatum = 0;             % use unlimited datum record dimension
 S.Ndatum = obs.Ndatum;    % fixed datum dimension
 
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

S.state_flag_values  =(1:1:7);

S.state_flag_meanings=['zeta', blanks(1),                               ...
                       'ubar', blanks(1),                               ...
                       'vbar', blanks(1),                               ...
                       'u', blanks(1),                                  ...
                       'v', blanks(1),                                  ...
                       'temperature', blanks(1),                        ...
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


S.origin_flag_values = 1:1:12;

S.origin_flag_meanings = ['gridded_AVISO_SLA', blanks(1),               ...
                          'blended_SST', blanks(1),                     ...
                          'XBT_Met_Office', blanks(1),                  ...
                          'CTD_temperature_Met_Office', blanks(1),      ...
                          'CTD_salinity_Met_Office', blanks(1),         ...
                          'ARGO_temperature_Met_Office', blanks(1),     ...
                          'ARGO_salinity_Met_Office', blanks(1),        ...
                          'CTD_temperature_CalCOFI', blanks(1),         ...
                          'CTD_salinity_CalCOFI', blanks(1),            ...
                          'CTD_temperature_GLOBEC', blanks(1),          ...
                          'CTD_salinity_GLOBEC', blanks(1),             ...
                          'buoy_temperature_Met_Office'];

%  The attribute association between 'flag_values' and 'flag_meanings'
%  is difficult to read when a moderate number of flags are use. To
%  aid the decoding, two readable global attributes are added:
%  'state_variables' and 'obs_provenance' which are stored in
%  S.global_variables and S.global_provenance, respectively.
%
%  NOTE: the 'state_variables' attribute include the units.

newline=sprintf('\n');

S.global_variables=[newline,                                            ...
       '1: free-surface (m) ', newline,                                 ...
       '2: vertically integrated u-momentum component (m/s) ', newline, ...
       '3: vertically integrated v-momentum component (m/s) ', newline, ...
       '4: u-momentum component (m/s) ', newline,                       ...
       '5: v-momentum component (m/s) ', newline,                       ...
       '6: potential temperature (Celsius) ', newline,                  ...
       '7: salinity (nondimensional)'];

S.global_provenance=[newline,                                           ...
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

S.global_sources=[newline,                                              ...
       'http://opendap.aviso.oceanobs.com/thredds/dodsC ', newline,     ...
        sst_URL, newline,                                               ...
       'http://hadobs.metoffice.com/en3'];

%--------------------------------------------------------------------------
%  Create 4D-Var observations NetCDF file.
%--------------------------------------------------------------------------

[status]=c_observations(S);

%  Update 'units' attribute for time variables. Notice that the time
%  of the observations in the NetCDF file is in DAYS.

avalue=['days since ' datestr(mybasedate,31)];

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

obs.grid_Lm_Mm_N = S.grid_Lm_Mm_N;
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

[status]=c_observations(OBS,SUPfile);

avalue=['days since ' datestr(mybasedate,31)];

[status]=nc_attadd(SUPfile,'units',avalue,'survey_time');
[status]=nc_attadd(SUPfile,'calendar','gregorian','survey_time');

[status]=nc_attadd(SUPfile,'units',avalue,'obs_time');
[status]=nc_attadd(SUPfile,'calendar','gregorian','obs_time');

[status]=obs_write(SUPfile,OBS);

disp(' ');
disp('Done.');
disp(' ');

