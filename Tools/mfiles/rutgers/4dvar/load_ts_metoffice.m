function [T,S]=load_ts_metoffice(ncfile, GRDfile, StartDay, EndDay)

%
% LOAD_TS_metoffice:  Loads T,S data from the UK Met Office datasets
%
% [T,S]=load_ts_metoffice(ncfile, GRDfile, StartDay, EndDay)
% 
% This function loads temperature and salinity profile data from the UK Met
% Office Hadley Centre, quality controlled, EN3 observation datasets:
%
%    http://hadobs.metoffice.com/en3/data/EN3_v2a/download_EN3_v2a.html
%
% The hydrographic data are available from 1950 to the present and stored
% in annual tar files (said, EN3_v2a_Profiles_2004.tar) containing a NetCDF
% for each month (said, EN3_v2a_Profiles_200401.nc). Since this dataset is
% not available in an OpenDAP server, the user needs to get the tar files
% and extract the files in it.
%
% Only the following observations types (WMO code table 177) are processed:
%
%    401 - XBT
%    741 - TESAC, XCTD, CTD
%    820 - buoys, thermistor chain
%    831 - ARGO floats
%   
% On Input:
%
%    ncfile        UK Met Office dataset monthly NetCDF file name(s)
%                    (string or cell array). If the file is compressed,
%                    it will uncompress using "gunzip" and compress
%                    back after processing with "gzip". Original monthly
%                    files are "gzip" compressed.
%
%    GRDname       NetCDF grid file name (string, OPTIONAL).
%                    
%    StartDay      Starting period of interest (date number, OPTIONAL)
%
%                    Example:   startday=datenum(2004,1, 1)=731947
%
%    EndDay        Ending   period of interest (date number, OPTIONAL)
%
%                    Example:   startday=datenum(2004,1,15)=731961
%
%                    See Matlab intrinsic datenum(Y,Mo,D,H,Mi,S)
%                    for details.
%
%    If the OPTIONAL argument GRDname is provided, only the data inside
%    the application grid is extracted. Similarly, OPTIONAL arguments
%    StartDay and EndDay can be used to extract available data for the
%    specified period.
%
% On Output:
%
%    T             Potential temperature profiles (structure array):
%
%                    T.time         time (date number)
%                    T.lon          longitude (degree_east)
%                    T.lat          latitude
%                    T.depth        depth (meter, negative)
%                    T.WMO          WMO instrument type (numeric)
%                    T.value        potential temperature value (Celsius)
%
%                    T.Nsurvey      number of observation surveys
%                    T.Ndatum       total number of observation
%                    T.Nobs         number of observations per survey
%                    T.survey_time  time for each survey
%
%    S             Salinity profiles (structure array):
%
%                    S.time         time (date number)
%                    S.lon          longitude (degree_east)
%                    S.lat          latitude
%                    S.depth        depth (meter, negative)
%                    S.WMO          WMO instrument type (numeric)
%                    S.value        salinity value
%
%                    S.Nsurvey      number of observation surveys
%                    S.Ndatum       total number of observation
%                    S.Nobs         number of observations per survey
%                    S.survey_time  time for each survey
%

% svn $Id: load_ts_metoffice.m 996 2020-01-10 04:28:56Z arango $
%=========================================================================%
%  Copyright (c) 2002-2020 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%

debugging  = true;
got_grid   = false;
got_period = false;

%  Check arguments.

if (nargin > 1),
  got_grid = true;
end,

if (nargin > 2),
  got_period = true;
  if (StartDay > EndDay),
    error([' LOAD_TS_METOFFICE: Your starting time must be greater',    ...
           ' than the ending time']);
  end,
end,

%  Determine number of input monthly files.

if (iscell(ncfile)),
  Nfiles = length(ncfile);       % multiple monthly file names cell array
else
  Nfiles = 1;                    % single monthly file name string
end

%  Read in application grid longitude and latitude. Set grid application
%  polygon.

if (got_grid),
  rlon = nc_read(GRDfile,'lon_rho');
  rlat = nc_read(GRDfile,'lat_rho');

  [Im,Jm]=size(rlon);

  Xbox=[squeeze(rlon(:,1));                                             ...
        squeeze(rlon(Im,2:Jm))';                                        ...
        squeeze(flipud(rlon(1:Im-1,Jm)));                               ...
        squeeze(fliplr(rlon(1,1:Jm-1)))'];

  Ybox=[squeeze(rlat(:,1));                                             ...
        squeeze(rlat(Im,2:Jm))';                                        ...
        squeeze(flipud(rlat(1:Im-1,Jm)));                               ...
        squeeze(fliplr(rlat(1,1:Jm-1)))'];
end

%  Set observations dynamical fields (cell array) in output structures.

field_list = {'time', 'lon', 'lat', 'depth', 'WMO', 'QC', 'value'};

%  Initialize.

for value = field_list,
  field = char(value);
  T.(field) = [];
  S.(field) = [];
end

%--------------------------------------------------------------------------
%  Read in UK Met Office dataset(s).
%--------------------------------------------------------------------------

tindex    = [];                         % data do not have record dimension
FillValue = NaN;                        % replace fill values with NaN
epoch     = datenum([1950 1 1 0 0 0]);  % dataset time origin

N = 1;                                  % input data file counter

while (N <= Nfiles),

%  Get monthly file name.

  if (iscell(ncfile)),
    ncname = char(ncfile(N));
  else
    ncname = ncfile;
  end

%  Check if input file is compressed.  If so uncmpress before processing.
%  It slows down the processing a little but it saves disk space.

  lstr = length(ncname);
  uncompress = false;
  if (strfind(ncname,'.gz')),
    s = unix(['gunzip ', ncname]);
    ncname = ncname(1:lstr-3);
    uncompress = true;
  end

%  Read in time, spatial location, and WMO instrument type string.
%
%    JULD           day since midnight 1/1/1950           (1:Nprof)
%    LONGITUDE      observation longitude, -180 to 180    (1:Nprof)
%    LATITUDE       observation latitude, -90 to 90       (1:Nprof)
%    DEPH_CORRECTED observation depth, negative, meters   (1:Nlev, 1:Nprof)
%    WMO_INST_TYPE  WMO instrument type, code table 1770  (4xNprof)

  time  = nc_read(ncname, 'JULD'          , [], FillValue) + epoch;
  lon   = nc_read(ncname, 'LONGITUDE'     , [], FillValue);
  lat   = nc_read(ncname, 'LATITUDE'      , [], FillValue);
  depth = nc_read(ncname, 'DEPH_CORRECTED', [], FillValue);
  WMO   = nc_read(ncname, 'WMO_INST_TYPE');

  [Nlev, Nprof] = size(depth);

%  Make sure that depths are negative.

  depth = - abs(depth);

%  Convert WMO code to numbers and replicate for each level.

  WMO   = str2num(WMO');

%  Replicate data to each level.

  time  = repmat(time', [Nlev 1]);
  lon   = repmat(lon' , [Nlev 1]);
  lat   = repmat(lat' , [Nlev 1]);
  WMO   = repmat(WMO' , [Nlev 1]);

%  Read in corrected potential temperature (Celsius) and practical salinity
%  of each level of each profile (1:Nlev, 1:Nprof).  Replace fill values
%  with NaNs.

  temp  = nc_read(ncname, 'POTM_CORRECTED', tindex, FillValue);   
  salt  = nc_read(ncname, 'PSAL_CORRECTED', tindex, FillValue);

%  Read in quality control character (0 to 9) for potential temperature and
%  salinity (1:Nlev, 1:Nprof). The fill value is '0'. Values greater than 4
%  should be rejected. We are only processing data with a quality control
%  flag of 1.

  temp_qc = nc_read(ncname, 'POTM_CORRECTED_QC', tindex, '0');
  salt_qc = nc_read(ncname, 'PSAL_CORRECTED_QC', tindex, '0');

%  Convert data to a vector to facilitate screening. Need the transpose
%  here to allow expansion of the structure with more data from other
%  files.

  Temp.time  = time(:);       Salt.time  = time(:);
  Temp.lon   = lon(:);        Salt.lon   = lon(:);
  Temp.lat   = lat(:);        Salt.lat   = lat(:);
  Temp.depth = depth(:);      Salt.depth = depth(:);
  Temp.WMO   = WMO(:);        Salt.WMO   = WMO(:);
  Temp.QC    = temp_qc(:);    Salt.QC    = salt_qc(:);
  Temp.value = temp(:);       Salt.value = salt(:);

  NobsT(1) = length (Temp.time);
  NobsS(1) = length (Salt.time);

  clear WMO depth lat lon salt salt_qc temp temp_qc time

%  Screen and remove data for quality control and missing values.

  ind_T = find((Temp.QC ~= '1')  | isnan(Temp.value) |                  ...
               isnan(Temp.lon)   | isnan(Temp.lat)   |                  ...
               isnan(Temp.depth) | isnan(Temp.time));

  ind_S = find((Salt.QC ~= '1')  | isnan(Salt.value) |                  ...
               isnan(Salt.lon)   | isnan(Salt.lat)   |                  ...
               isnan(Salt.depth) | isnan(Salt.time));

  if (~isempty(ind_T)),
    for value = field_list,
      field = char(value);
      Temp.(field)(ind_T) = [];
    end
  end

  if (~isempty(ind_S)),
    for value = field_list,
      field = char(value);
      Salt.(field)(ind_S) = [];
    end
  end

  NobsT(2) = length (Temp.time);
  NobsS(2) = length (Salt.time);

  clear ind_S ind_T

%  Extract only requested WMO code data.

  ind_T = find((Temp.WMO ~= 401) & (Temp.WMO ~= 741) &                  ...
               (Temp.WMO ~= 820) & (Temp.WMO ~= 831));

  ind_S = find((Salt.WMO ~= 401) & (Salt.WMO ~= 741) &                  ...
               (Salt.WMO ~= 820) & (Salt.WMO ~= 831));
  
  if (~isempty(ind_T)),
    for value = field_list,
      field = char(value);
      Temp.(field)(ind_T) = [];
    end
  end

  if (~isempty(ind_S)),
    for value = field_list,
      field = char(value);
      Salt.(field)(ind_S) = [];
    end
  end

  NobsT(3) = length (Temp.time);
  NobsS(3) = length (Salt.time);

  clear ind_S ind_T

%  If application grid is provided, extract observations inside the grid
%  polygon.

  if (got_grid),
    bounded = false(size(Temp.lon));
    [IN,ON] = inpolygon(Temp.lon, Temp.lat, Xbox, Ybox);
    bounded(IN) = true;
    bounded(ON) = true;
   
    ind_T = find(~bounded);

    if (~isempty(ind_T)),
      for value = field_list,
        field = char(value);
        Temp.(field)(ind_T) = [];
      end
    end

    bounded = false(size(Salt.lon));
    [IN,ON] = inpolygon(Salt.lon, Salt.lat, Xbox, Ybox);
    bounded(IN) = true;
    bounded(ON) = true;
   
    ind_S = find(~bounded);

    if (~isempty(ind_S)),
      for value = field_list,
        field = char(value);
        Salt.(field)(ind_S) = [];
      end
    end

    NobsT(4) = length (Temp.time);
    NobsS(4) = length (Salt.time);

    clear IN ON bounded ind_S ind_T
  end

%  If data period is provided, extract observations for the specified
%  times.

  if (got_period),
    ind_T = find(Temp.time >= StartDay & EndDay <= Temp.time);
    ind_S = find(Salt.time >= StartDay & EndDay <= Salt.time);

    if (~isempty(ind_T)),
      for value = field_list,
        field = char(value);
        Temp.(field)(ind_T) = [];
      end
    end

    if (~isempty(ind_S)),
      for value = field_list,
        field = char(value);
        Salt.(field)(ind_S) = [];
      end,
    end,

    NobsT(5) = length (Temp.time);
    NobsS(5) = length (Salt.time);

    clear ind_S ind_T
  end

%  Load extracted data into structure.  Sort in ascending order of time.

  if (debugging),
    disp(' ');
    disp(['                            Observations file = ',           ...
          ncname]);
    disp(['                  Number of observations read = ',           ...
          num2str(NobsT(1),'%8.8i'), '  ', num2str(NobsS(1),'%8.8i')]);
    disp([' Number after QC screening and missing values = ',           ...
          num2str(NobsT(2),'%8.8i'), '  ', num2str(NobsS(2),'%8.8i')]);
    disp(['         Number after WMO intrument screening = ',           ...
          num2str(NobsT(3),'%8.8i'), '  ', num2str(NobsS(3),'%8.8i')]);
    if (got_grid),
      disp(['                Number inside aplication grid = ', ...
            num2str(NobsT(4),'%8.8i'), '  ', num2str(NobsS(4),'%8.8i')]);
    end,
    if (got_period),
      disp(['                  Number after time screening = ',         ...
            num2str(NobsT(5),'%8.8i'), '  ', num2str(NobsS(5),'%8.8i')]);
    end
  end
  
  if (isempty(Temp.time) && isempty(Salt.time)),

    disp([' LOAD_TS_METOFFICE: no data extracted from file: ', ncname]);
  
  else   % found observations in this period
    
    [~,I] = sort(Temp.time, 'ascend');       % sort in time ascending order
  
    for value = field_list,
      field = char(value);
      T.(field) = [T.(field), transpose(Temp.(field)(I))];
    end

    [~,I] = sort(Salt.time, 'ascend');       % sort in time ascending order
  
    for value = field_list,
      field = char(value);
      S.(field) = [S.(field), transpose(Salt.(field)(I))];
    end
  
    clear I Temp Salt Y
    
  end


%  Increase file name counter. If applicable, compress back the input
%  file.
  
  N = N + 1;

  if (uncompress),
    s = unix(['gzip ', ncname]);
  end

end

%  Add number of observations, survery time, and number of obserbations
%  per survey time.

Tsurvey   = unique(T.time);
Ssurvey   = unique(S.time);

T.Nsurvey = length(Tsurvey);
S.Nsurvey = length(Ssurvey);

T.Ndatum  = length(T.time);
S.Ndatum  = length(S.time);

for n=1:T.Nsurvey,
  ind = find(T.time == Tsurvey(n));
  T.Nobs(n) = length(ind);
end
T.Nobs = T.Nobs';

for n=1:S.Nsurvey,
  ind = find(S.time == Ssurvey(n));
  S.Nobs(n) = length(ind);
end
S.Nobs = S.Nobs';

T.survey_time = Tsurvey;
S.survey_time = Ssurvey;

T.WMO_types = unique(T.WMO);
S.WMO_types = unique(S.WMO);

%  Remove QC field from structure.

T = rmfield(T, 'QC');
S = rmfield(S, 'QC');

return
