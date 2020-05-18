function [Sout]=super_obs(Sinp);

%
% SUPER_OBS:  Creates super observations when necessary
%
% [Sout]=super_obs(Sinp)
%
% This function checks the provided observation data and creates
% super observations when there are more than one meassurement of
% the same state variable per grid cell. At input, Sinp is either a
% 4D-Var observation NetCDF file or data structure.
%
% An additional field (Sout.std) is added to the output observation
% structure containing the standard deviation of the binning which
% can be used as observation error. You may choose to assign this
% value as observation error before writing to NetCDF file:
%
%    Sout.error = max(Sout.error, Sout.std)
%
% That is, the larger observation variace is chosen. A zero value
% for "Sout.std" indicates that the observation did not required
% binning.  Use "c_observations.m" to create NetCDF file and
% "obs_write" to write out data.
%
% On Input:
%
%    Sinp    Observations data structure or NetCDF file name
%
%
% On Output:
%
%    Sout    Binned observations data (structure array):
%
%              Sout.ncfile       NetCDF file name (string)
%              Sout.Ndatum       total number of observations
%              Sout.spherical    spherical grid switch
%              Sout.Nobs         number of observations per survey
%              Sout.survey_time  time for each survey time
%              Sout.variance     global variance per state variable
%              Sout.type         state variable associated with observation
%              Sout.time         time for each observation
%              Sout.depth        depth of observation
%              Sout.Xgrid        observation fractional x-grid location
%              Sout.Ygrid        observation fractional y-grid location
%              Sout.Zgrid        observation fractional z-grid location
%              Sout.error        observation error
%              Sout.value        observation value
%
%              Sout.std          binning observation standard deviation
%
%            The following optional variables will be process if available:
%
%              Sout.provenance   observation origin
%              Sout.lon          observation longitude
%              Sout.lat          observation latitude
% 
  
% svn $Id: super_obs.m 996 2020-01-10 04:28:56Z arango $
%===========================================================================%
%  Copyright (c) 2002-2020 The ROMS/TOMS Group                              %
%    Licensed under a MIT/X style license           Hernan G. Arango        %
%    See License_ROMS.txt                           John L. Wilkin          %
%===========================================================================%

%  Read observations if 'Sinp' is a 4D-Var observation NetCDF file.

if (ischar(Sinp)),
  Sinp=obs_read(Sinp);
end,

Lm = Sinp.grid_Lm_Mm_N(1);
Mm = Sinp.grid_Lm_Mm_N(2);
N  = Sinp.grid_Lm_Mm_N(3);

%  Check if 'provenace', 'lon', and 'lat' fields are available.

has.lonlat=false;
if (isfield(Sinp,'lon') && isfield(Sinp,'lat')),
  has.lonlat=true;
end,

has.provenance=false;
if (isfield(Sinp,'provenance')),
  has.provenance=true;
end,

%  Get number of data surveys.

Nsurvey=length(Sinp.survey_time);

%  Set observations dynamical fields (cell array) in structure, S.

field_list = {'Xgrid', 'Ygrid', 'Zgrid', 'depth', 'error', 'value'};

if (has.lonlat),
  field_list = [field_list, 'lon', 'lat'];
end,

if (has.provenance),
  field_list = [field_list, 'provenance'];
end,

%  Insure that the vector fields in the input structure have the
%  singleton as the first dimension to allow vector concatenation.

for value = field_list,
  field = char(value);
  if (size(Sinp.(field),1) > 1),
    Sinp.(field) = transpose(Sinp.(field));
  end,
end,
if (size(Sinp.survey_time,1) > 1),
  Sinp.survey_time = transpose(Sinp.survey_time);
end,

%  Add binning standard deviation field to input structure. Initialize
%  to input observation error.

Sinp.std = zeros(size(Sinp.error));

%----------------------------------------------------------------------------
%  Find observations associated with the same state variable.
%----------------------------------------------------------------------------

state_vars=unique(Sinp.type);
Nstate=length(state_vars);

%  Initialize output structure.

Sout.ncfile       = Sinp.ncfile;
Sout.Nsurvey      = Nsurvey;
Sout.Nstate       = length(Sinp.variance);
Sout.Ndatum       = [];
Sout.spherical    = Sinp.spherical;
Sout.Nobs         = zeros(size(Sinp.Nobs));
Sout.survey_time  = Sinp.survey_time;
Sout.variance     = Sinp.variance;

Sout.type         = [];            % not included in dynamic fields
Sout.time         = [];            % for efficiency

for value = field_list,
  field = char(value);             % convert from cell to string
  Sout.(field) = [];               % initilize to empty
end,
Sout.std = [];                     % binning standard deviation

%----------------------------------------------------------------------------
%  Compute super observations when needed.
%----------------------------------------------------------------------------

for m=1:Nsurvey,   %%% SURVEY TIME LOOP %%%

%  Extract observations with the same survey time.

  ind_t=find(Sinp.time == Sinp.survey_time(m));

  T.type = Sinp.type(ind_t);                  % not included in dynamic
  T.time = Sinp.time(ind_t);                  % fields for efficiency
  
  for value = field_list,
    field = char(value);                      % convert from cell to string
    T.(field) = Sinp.(field)(ind_t);          % initilize to Sinp structure
  end,
  T.std = Sinp.std(ind_t);                    % binning standard deviation
  
  for n=1:Nstate,  %%% STATE VARIABLE LOOP %%%

%   Now extract observations associated with the same state variable.

    ind_v=find(T.type == state_vars(n));

    if isempty(ind_v),
      continue;                               % none found, exit n-loop
    end,

    V.type = T.type(ind_v);                   % not included in dynamic
    V.time = T.time(ind_v);                   % fields for efficiency
    
    for value = field_list,
      field = char(value);                    % convert from cell to string
      V.(field) = T.(field)(ind_v);           % initilize to T structure
    end,
    V.std = T.std(ind_v);                     % binning standard deviation

%  Set binning parameters. The processing is done in fractional (x,y,z) grid
%  locations

    Xmin = min(V.Xgrid);
    Xmax = max(V.Xgrid);

    Ymin = min(V.Ygrid);
    Ymax = max(V.Ygrid);

    Zmin = min(V.Zgrid);
    Zmax = max(V.Zgrid);

    dx = 1.0;
    dy = 1.0;
    dz = 1.0;
    
    minObs= 1;

%  Compute the index in each dimension of the grid cell in which the
%  observation is located.

    Xbin = 1.0 + floor((V.Xgrid - Xmin) ./ dx);
    Ybin = 1.0 + floor((V.Ygrid - Ymin) ./ dy);
    Zbin = 1.0 + floor((V.Zgrid - Zmin) ./ dz);

%  Similarly, compute the maximum averaging grid size.

    Xsize = 1.0 + floor((Xmax - Xmin) ./ dx);
    Ysize = 1.0 + floor((Ymax - Ymin) ./ dy);
    Zsize = 1.0 + floor((Zmax - Zmin) ./ dz);
    
    matsize = [Ysize, Xsize, Zsize];

%  Combine the indices in each dimension into one index. It is like stacking
%  all the matrix in one column vector.

    varInd    = transpose(sub2ind(matsize, Ybin, Xbin, Zbin));
    onesCol   = transpose(ones(size(varInd)));
    maxVarInd = max(varInd);

%  Accumulate values in bins using "accumarray" function. Count how many
%  observations fall in each bin.

    count = accumarray(varInd, onesCol, [], @sum, [], true);
    
%  Bins with no observations are not keept. Index vector "isdata" will be
%  used in the conversion from the sparse output from "accumarray" to the
%  corresponding full vector assigned to Sout.

    isdata = find(count ~= 0);
    Nsuper = length(isdata);
    
%  Loop through list of fields that require binning. The "accumarray"
%  function take the sum of all values having the same bin index. The last
%  argument activates sparse. It turns out that method below is (somewhat)
%  faster than asking directly for the mean, thus:
%
%     binned = accumarray(varInd, V.(field), [], @mean, [], true);

    for fval = field_list,
      field = char(fval);
      
      binned = accumarray(varInd, V.(field), [], @sum, [], true);
      binned = full(binned(isdata) ./ count(isdata));
      
      Sout.(field) = [Sout.(field) transpose(binned)];
    
      if (strcmp(field, 'value')),        % save binned observation value
        Vmean = binned;                   % to compute its associated
      end,                                % variance (biased estimate)
    end,

%  Compute the standard deviation of the binning. The user can use these
%  values to estimate observation error.

    binned = accumarray(varInd, V.value.^2, [], @sum, [], true);
    binned = full(binned(isdata) ./ count(isdata)) - Vmean .^ 2;

    Sout.std = [Sout.std transpose(sqrt(binned))];

%  Set "type" and "time" fields for binned observations.

    Sout.type  = [Sout.type, ones([1 Nsuper]).*state_vars(n)];
   
    Sout.time  = [Sout.time, ones([1 Nsuper]).*Sinp.survey_time(m)];

    Sout.Nobs(m) = Sout.Nobs(m) + Nsuper;

%  Make sure that "provenance", if available, is not a fractional number.
%  Binning provenance is weird here. Hopefully, this is always a full
%  number since the extrategy is to create one single observation file per
%  dataset or instrument. Then, the 4D-Var observations for ROMS are merged
%  together using "obs_merge.m".
    
    if (has.provenance),
      Sout.provenance = floor(Sout.provenance);
    end,
    
  end,             %%% end of STATE VARIABLE LOOP %%% 
end,               %%% end of SURVEY TIME LOOP    %%%

Sout.Ndatum = sum(Sout.Nobs);

%  Add other information fields from input structure. This will facilitates
%  writing of the new structure.

if (has.lonlat),
  Sout.do_longitude=true;
  Sout.do_latitude =true;
end,

if (isfield(Sinp,'title'));
  Sout.title = Sinp.title;
end,

if (isfield(Sinp,'state_flag_values'));
  Sout.state_flag_values = Sinp.state_flag_values;
end,

if (isfield(Sinp,'state_flag_meanings'));
  Sout.state_flag_meanings = Sinp.state_flag_meanings;
end,

if (isfield(Sinp,'origin_flag_values'));
  Sout.origin_flag_values = Sinp.origin_flag_values;
end,

if (isfield(Sinp,'origin_flag_meanings'));
  Sout.origin_flag_meanings = Sinp.origin_flag_meanings;
end,

if (isfield(Sinp,'grd_file'));
  Sout.grd_file = Sinp.grd_file;
end,

if (isfield(Sinp,'grid_Lm_Mm_N'));
  Sout.grid_Lm_Mm_N = Sinp.grid_Lm_Mm_N;
end,

if (isfield(Sinp,'global_variables'));
  Sout.global_variables = Sinp.global_variables;
end,

if (isfield(Sinp,'global_provenance'));
  Sout.global_provenance = Sinp.global_provenance;
end,

if (isfield(Sinp,'global_sources'));
  Sout.global_sources = Sinp.global_sources;
end,

return
