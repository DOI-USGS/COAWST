function D = desroziers(ncname)

%
% DESROZIERS:  Computes 4D-Var consistency diagnostics
%
% D = desroziers(ncname)
%
% This function computes Desroziers et al. (2005) consistency diagnostics
% for the innovations, background errors, observation errors, and analysis
% error from ROMS output 4D-Var (MODname) NetCDF file(s).  The diagnotics
% are separated by surface and sub-surface state variables. Usually, there
% are diagnostics for 13 state variables classes in the following order:
%
%    1   zeta      free-surface (altimetry SSH observations)
%    2   Usur      surface U-velocity (HFR observations)
%    3   Usub      sub-surface U-velocity (in situ observations)
%    4   Uvel      surface and sub-surface U-velocity
%    5   Vsur      surface V-velocity (HFR observations)
%    6   Vsub      sub-surface V-velocity (in situ observations)
%    7   Vvel      surface and sub-surface V-velocity
%    8   Tsur      surface temperature (satellite SST observations)
%    9   Tsub      sub-surface temperature (in situ observations)
%   10   Temp      surface and sub-surface temperature
%   11   Ssur      surface salinity (satellite SSS observations)
%   12   Ssub      sub-surface salinity (in situ observations)
%   13   Salt      surface and sub-surface salinity
%
% On Input:
%
%    ncname        ROMS 4D-Var NetCDF filename(s) (sring or cell array)
%
% On Output:
%
%    D             Desroziers consistency diagnostics (struct)
%                    (state variable units)
%
%                    D.ncname        4D-Var NetCDF filename(s)
%                    D.Vstate        state variables
%                    D.Ncycles       number of 4D-Var cycles (files)
%                    D.epoch         reference time (datenum)
%                    D.time          days since reference time
%                    D.date          calendar date (DD-MMM-YYYY)
%                    D.Ndatum        number of viable observations
%                    D.reject        number of rejected observations
%                    D.count         number of observations per variable
%                    D.EsigmaI       diagnostic on innovations
%                    D.EsigmaB       diagnostic on background error
%                    D.EsigmaO       diagnostic on observation error
%                    D.EsigmaA       diagnostic of analysis error
%                    D.sigmaB        background error standard deviation
%                    D.sigmaO        observation error standard deviation 
%
% where:
%
%   EsigmaI = E[d_{b}^{o} (d_{b}^{o})']   
%   EsigmaB = E[d_{b}^{a} (d_{b}^{o})']
%   EsigmaO = E[d_{a}^{o} (d_{b}^{o})']
%   EsigmaA = E[d_{b}^{a} (d_{a}^{o})']
%
%   d_{b}^{o} = O - B                   (innovation)
%   d_{b}^{a} = A - B                   (increment)
%   d_{a}^{o} = O - A                   (residual)
%
% NOTES/WARNINGS:
%
%   * The diagnostics only provide information about the diagonal of
%     B and R error covariances and not about the cross-covariaces.
%   * EsigmaI is always positive
%   * If EsigmaB or EsigmaO is negative (complex number since we are
%     taken the squared-root), there is an inconsistency in the error
%     model assumptions and the User needs to adjust the background
%     and observation error covariances in a substantial way.
%   * Optimally, EsigmaA is positive but often comes out negative. It
%     indicates that B and R are not correct. But we know that anyway.
%
%--------------------------------------------------------------------------
% There are various ways in Matlab to create a cell array of various
% NetCDF files across many directories of 4D-Var cycles:
%
%   S = dir('Run01/201*/*_mod*.nc);     % creates a structure of files
%   [~,K] = sort([S.datenum]);          % sort by date of creation
%   ncname = {S(K).name};               % create ordered cell array
%
% For example, we may have the following directory tree containing a full
% year of 4D-Var cycles:
%
%   Run01/2017.01.01/wc12_mod_20170101.nc
%   Run01/2017.01.06/wc12_mod_20170106.nc
%   ...
%   Run01/2017.12.27/wc12_mod_20171227.nc
%     
%--------------------------------------------------------------------------
%
% Reference:
%
%   Desroziers, G., L. Berre, B. Chapnik, and P. Poli, 2005: Diagnosis of
%     observation, background and analysis-error statisticsin observation
%     space, Q.J.R. Meteorol. Soc., 131, 3385-3396, doi: 10.1256/qj.05.108.
%  

% svn $Id$
%=========================================================================%
%  Copyright (c) 2002-2025 The ROMS Group                                 %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.md                            Hernan G. Arango      %
%=========================================================================%
  
% Initialize

Lprint = true;                  % print out diagnotics.

FlipMultiplication = true;      % dot vector product flip:  A * B = B * A 
                                % (cummulative property)

% Determine number of input NetCDF files to process.
  
if (iscell(ncname))
  nfiles = length(ncname);

  obs_type = nc_read(ncname{1}, 'obs_type');
  Avalue   = nc_getatt(ncname{1}, 'units', 'obs_time');
else
  nfiles = 1;

  obs_type = nc_read(ncname, 'obs_type');
  Avalue   = nc_getatt(ncname, 'units', 'obs_time');
end

% Get reference time ('days since ...');

ind   = strfind(lower(Avalue), 'since');
rtime = Avalue(ind+6:ind+24);
epoch = datenum(rtime);

% Determine state variables to process.

types = unique(obs_type);
is3d  = any(types > 3);

if (is3d)
  nvars  = 13;
  Vstate = {'zeta', 'Usur', 'Usub', 'Uvel',                             ...
                    'Vsur', 'Vsub', 'Vvel',                             ...
                    'Tsur', 'Tsub', 'Temp',                             ...
                    'Ssur', 'Ssub', 'Salt'};

  izeta  = 1;          % SSH observations
  iUsur  = 2;          % surface U-velocity (HF radar observations)
  iUsub  = 3;          % sub-surface U-velocity (in situ observations)
  iUvel  = 4;          % surface and sub-surface U-velocity
  iVsur  = 5;          % surface V-velocity (HF radar observations)
  iVsub  = 6;          % sub-surface V-velocity (in situ observations)
  iVvel  = 7;          % surface and sub-surface V-velocity
  iTsur  = 8;          % surface temperature (satellite SST observations)
  iTsub  = 9;          % sub-surface temperature (in situ observations)
  iTemp  = 10;         % surface and sub-surface temperature
  iSsur  = 11;         % surface salinity (satellite SSS observations)
  iSsub  = 12;         % sub-surface salinity (in situ observations)
  iSalt  = 13;         % surface and sub-surface salinity
else
  nvars  = 3;
  Vstate = {'zeta', 'Ubar', 'Ubar'};

  izeta  = 1;          % SSH observations
  iUbar  = 2;          % vertically-integraded V-velocity
  iVbar  = 3;          % vertically-integraded V-velocity
end

% Initialize output structure.

D = struct('ncname'     , [], 'Vstate'    , [], 'Ncycles'   , [],       ...
           'epoch'      , [], 'time'      , [], 'date'      , [],       ...
           'Ndatum'     , [], 'reject'    , [], 'count'     , [],       ...
           'EsigmaI'    , [], 'EsigmaB'   , [], 'EsigmaO'   , [],       ...
           'EsigmaA'    , [],                                           ...
           'sigmaB'     , [], 'sigmaO'    , []);

D.ncname  = ncname;
D.Vstate  = Vstate;
D.Ncycles = nfiles;
D.epoch   = epoch;
D.time    = nan(nfiles);
D.date    = cell(nfiles);

D.Ndatum  = nan(nfiles);
D.reject  = nan(nfiles);
D.count   = nan(nvars, nfiles);

D.EsigmaI = nan(nvars, nfiles);
D.EsigmaB = nan(nvars, nfiles);
D.EsigmaO = nan(nvars, nfiles);
D.EsigmaA = nan(nvars, nfiles);

D.sigmaB  = nan(nvars, nfiles);
D.sigmaO  = nan(nvars, nfiles);

% Set NetCDF variables to process.

VarList = {'obs_error', 'obs_provenance', 'obs_scale', 'obs_time',      ...
	   'obs_depth', 'obs_type',  'obs_value',                       ...
	   'innovation', 'increment', 'residual',                       ...
           'NLmodel_initial', 'NLmodel_final', 'NLmodel_value',         ...
	   'BgError_value'};

%--------------------------------------------------------------------------
% Process input ROMS 4D-Var NetCDF file(s).
%--------------------------------------------------------------------------

if (Lprint)
  disp(blanks(1));
end

for n = 1:nfiles
  if (iscell(ncname))
    fname = ncname{n};
  else
    fname = ncname;
  end

% Inquire about the contents of NetCDF file.

  I = nc_inq(fname);

  Nouter = I.Dimensions(strcmp({I.Dimensions.Name}, 'Nouter')).Length;

% Read in variables.

  for var = VarList
    field = char(var);
    vindx = strcmp({I.Variables.Name}, field);
    if (any(vindx))
      got.(field) = true;
      S.(field) = nc_read(fname, field);
    else
      got.(field) = false;
    end
  end
  
% Get 4D-Var cycle start time and date to facilitate plotting elsewhere.

  iatt = strcmp({I.Attributes.Name}, 'str_day');
  if (any(iatt))
    str_day = I.Attributes(iatt).Value;
  else
    ind = find(S.obs_scale > 0);
    str_day = floor(mean(S.obs_time(ind)));
  end
  str_date = datestr(epoch+str_day,1);
  
  D.time(n) = str_day;
  D.date(n) = {str_date};
  
% Check if 'innovation', 'increment', are 'residual' available. If not,
% we need to compute them.  These variables are not available in old 
% version of the MOD files. In some version of ROMS, we didnt' multiply
% by the screeing flag 'obs_scale'.

  if (~got.NLmodel_final)
    S.NLmodel_final = S.NLmodel_value(:,Nouter);
  end

  if (~got.innovation)                             % d_{b}^{o} = O-B
    S.innovation = S.obs_scale .*                                       ...
                  (S.obs_value - S.NLmodel_initial);
  else 
    S.innovation = S.innovation .* S.obs_scale;
  end
  
  if (~got.increment)                              % d_{b}^{a} = A-B
    S.increment  = S.obs_scale .*                                       ...
                   (S.NLmodel_final - S.NLmodel_initial);
  else
    S.increment  = S.increment .* S.obs_scale;
  end
  
  if (~got.residual)                               % d_{a}^{o} = O-A
    S.residual   = S.obs_scale .*                                       ...
                   (S.obs_value - S.NLmodel_final);
  else 
    S.residual   = S.residual  .* S.obs_scale;
  end

% Get unique observation types and provenances.

  types  = unique(S.obs_type);                     % state variable type
  origin = unique(S.obs_provenance);

% Determine the number of viable observations for current 4D-Var cycle.

  D.Ndatum(n) = sum(S.obs_scale > 0);
  D.reject(n) = sum(S.obs_scale == 0);

% Compute Desrozier consistency diagnostics for the innovations, background
% errors, observation errors, and analysis errors in terms of the state
% variable type. Notice that the squared root is taken so the diagnostics
% are the RMSE (state variable units).

  if (Lprint)
    disp(['** 4D-Var Cycle: ', sprintf('%4s',num2str(n)), ', ', str_date]);
    disp(blanks(1));
  end

  for m = 1:length(types)
    itype  = types(m);
    isur   = find(S.obs_type == itype & S.obs_depth == 0);  % surface
    isub   = find(S.obs_type == itype & S.obs_depth <  0);  % sub-surface
    ival   = find(S.obs_type == itype & itype > 3);         % 3D variable

    if (~isempty(isur))                      % surface observations
      ic = size(S.obs_error(isur) ~=0, 1);
      Count = sum(isur > 0);
      
      if (FlipMultiplication)
        EsigmaI_sur = sqrt(S.innovation(isur)' * S.innovation(isur) / ic);
        EsigmaB_sur = sqrt(S.innovation(isur)' * S.increment(isur)  / ic);
        EsigmaO_sur = sqrt(S.innovation(isur)' * S.residual(isur)   / ic);
        EsigmaA_sur = sqrt(S.residual(isur)'   * S.increment(isur)  / ic);
      else
        EsigmaI_sur = sqrt(S.innovation(isur)' * S.innovation(isur) / ic);
        EsigmaB_sur = sqrt(S.increment(isur)'  * S.innovation(isur) / ic);
        EsigmaO_sur = sqrt(S.residual(isur)'   * S.innovation(isur) / ic);
        EsigmaA_sur = sqrt(S.increment(isur)'  * S.residual(isur)   / ic);
      end
	
      sigmaB_sur = sqrt(nanmean(S.BgError_value(isur).*                 ...
                                S.BgError_value(isur)));
      sigmaO_sur = sqrt(nanmean(S.obs_error(isur)));

      switch (itype)
        case 1
	  index = izeta;                    % SSSH
        case 2
	  index = iUbar;                    % Ubar, is3d = false   
        case 3
	  index = iVbar;                    % Vbar, is3d = false
	case 4
          index = iUsur;                    % Usur, is3d = true    
        case 5
          index = iVsur;                    % Vsur, is3d = true
        case 6
	  index = iTsur;                    % SST,  is3d = true
        case 7
          index = iSsur;                    % SSS,  is3d = true
      end

      D.count  (index, n) = Count;
      D.EsigmaI(index, n) = EsigmaI_sur;
      D.EsigmaB(index, n) = EsigmaB_sur;
      D.EsigmaO(index, n) = EsigmaO_sur;
      D.EsigmaA(index, n) = EsigmaA_sur;
      
      D.sigmaB (index, n) = sigmaB_sur;
      D.sigmaO (index, n) = sigmaO_sur;
    
      if (Lprint)
        disp([sprintf('   %3s',char(Vstate(index))),' obs count   = ',  ...
              num2str(Count), ', (ic = ',num2str(ic),')']);
        disp([sprintf('   %3s',char(Vstate(index))),' EsigmaI_sur = ',  ...
              num2str(EsigmaI_sur)]);
        disp([sprintf('   %3s',char(Vstate(index))),' EsigmaB_sur = ',  ...
              num2str(EsigmaB_sur)]);
        disp([sprintf('   %3s',char(Vstate(index))),' EsigmaO_sur = ',  ...
              num2str(EsigmaO_sur)]);
        disp([sprintf('   %3s',char(Vstate(index))),' EsigmaA_sur = ',  ...
              num2str(EsigmaA_sur)]);
        disp([sprintf('   %3s',char(Vstate(index))),' sigmaB_sur  = ',  ...
              num2str(sigmaB_sur)]);
        disp([sprintf('   %3s',char(Vstate(index))),' sigmaO_sur  = ',  ...
              num2str(sigmaO_sur)]);
        disp(blanks(1));
      end
    end      
      
    if (~isempty(isub))                      % sub-surface observations
      ic = size(S.obs_error(isub) ~=0, 1);
      Count = sum(isub > 0);

      if (FlipMultiplication)
        EsigmaI_sub = sqrt(S.innovation(isub)' * S.innovation(isub) / ic);
        EsigmaB_sub = sqrt(S.innovation(isub)' * S.increment(isub)  / ic);
        EsigmaO_sub = sqrt(S.innovation(isub)' * S.residual(isub)   / ic);
        EsigmaA_sub = sqrt(S.residual(isub)'   * S.increment(isub)  / ic);
      else
        EsigmaI_sub = sqrt(S.innovation(isub)' * S.innovation(isub) / ic);
        EsigmaB_sub = sqrt(S.increment(isub)'  * S.innovation(isub) / ic);
        EsigmaO_sub = sqrt(S.residual(isub)'   * S.innovation(isub) / ic);
        EsigmaA_sub = sqrt(S.increment(isub)'  * S.residual(isub)   / ic);
      end

      sigmaB_sub = sqrt(nanmean(S.BgError_value(isub).*                 ...
                                S.BgError_value(isub)));
      sigmaO_sub = sqrt(nanmean(S.obs_error(isub)));
    
      switch (itype)
	case 4
          index = iUsub;                    % U,    is3d = true    
        case 5
          index = iVsub;                    % V,    is3d = true
        case 6
	  index = iTsub;                    % T,    is3d = true
        case 7
          index = iSsub;                    % S,    is3d = true
      end
    
      D.count  (index, n) = Count;
      D.EsigmaI(index, n) = EsigmaI_sub;
      D.EsigmaB(index, n) = EsigmaB_sub;
      D.EsigmaO(index, n) = EsigmaO_sub;
      D.EsigmaA(index, n) = EsigmaA_sub;
      
      D.sigmaB (index, n) = sigmaB_sub;
      D.sigmaO (index, n) = sigmaO_sub;

      if (Lprint)
        disp([sprintf('   %3s',char(Vstate(index))),' obs count   = ',  ...
              num2str(Count), ', (ic = ',num2str(ic),')']);
        disp([sprintf('   %3s',char(Vstate(index))),' EsigmaI_sub = ',  ...
              num2str(EsigmaI_sub)]);
        disp([sprintf('   %3s',char(Vstate(index))),' EsigmaB_sub = ',  ...
              num2str(EsigmaB_sub)]);
        disp([sprintf('   %3s',char(Vstate(index))),' EsigmaO_sub = ',  ...
              num2str(EsigmaO_sub)]);
        disp([sprintf('   %3s',char(Vstate(index))),' EsigmaA_sub = ',  ...
              num2str(EsigmaA_sub)]);
        disp([sprintf('   %3s',char(Vstate(index))),' sigmaB_sub  = ',  ...
              num2str(sigmaB_sub)]);
        disp([sprintf('   %3s',char(Vstate(index))),' sigmaO_sub  = ',  ...
              num2str(sigmaO_sub)]);
        disp(blanks(1));
      end
    end

    if (~isempty(ival))                      % sub-surface and sub-surface
      ic = size(S.obs_error(ival) ~=0, 1);   % observations
      Count = sum(ival > 0);

      if (FlipMultiplication)
        EsigmaI_val = sqrt(S.innovation(ival)' * S.innovation(ival) / ic);
        EsigmaB_val = sqrt(S.innovation(ival)' * S.increment(ival)  / ic);
        EsigmaO_val = sqrt(S.innovation(ival)' * S.residual(ival)   / ic);
        EsigmaA_val = sqrt(S.residual(ival)'   * S.increment(ival)  / ic);
      else
        EsigmaI_val = sqrt(S.innovation(ival)' * S.innovation(ival) / ic);
        EsigmaB_val = sqrt(S.increment(ival)'  * S.innovation(ival) / ic);
        EsigmaO_val = sqrt(S.residual(ival)'   * S.innovation(ival) / ic);
        EsigmaA_val = sqrt(S.increment(ival)'  * S.residual(ival)   / ic);
      end

      sigmaB_val = sqrt(nanmean(S.BgError_value(ival).*                 ...
                                S.BgError_value(ival)));
      sigmaO_val = sqrt(nanmean(S.obs_error(ival)));
    
      switch (itype)
	case 4
          index = iUvel;                    % U,    is3d = true    
        case 5
          index = iVvel;                    % V,    is3d = true
        case 6
	  index = iTemp;                    % T,    is3d = true
        case 7
          index = iSalt;                    % S,    is3d = true
      end
    
      D.count  (index, n) = Count;
      D.EsigmaI(index, n) = EsigmaI_val;
      D.EsigmaB(index, n) = EsigmaB_val;
      D.EsigmaO(index, n) = EsigmaO_val;
      D.EsigmaA(index, n) = EsigmaA_val;
      
      D.sigmaB (index, n) = sigmaB_val;
      D.sigmaO (index, n) = sigmaO_val;

      if (Lprint)
        disp([sprintf('   %3s',char(Vstate(index))),' obs count   = ',  ...
              num2str(Count), ', (ic = ',num2str(ic),')']);
        disp([sprintf('   %3s',char(Vstate(index))),' EsigmaI_val = ',  ...
              num2str(EsigmaI_val)]);
        disp([sprintf('   %3s',char(Vstate(index))),' EsigmaB_val = ',  ...
              num2str(EsigmaB_val)]);
        disp([sprintf('   %3s',char(Vstate(index))),' EsigmaO_val = ',  ...
              num2str(EsigmaO_val)]);
        disp([sprintf('   %3s',char(Vstate(index))),' EsigmaA_val = ',  ...
              num2str(EsigmaA_val)]);
        disp([sprintf('   %3s',char(Vstate(index))),' sigmaB_val  = ',  ...
              num2str(sigmaB_val)]);
        disp([sprintf('   %3s',char(Vstate(index))),' sigmaO_val  = ',  ...
              num2str(sigmaO_val)]);
        disp(blanks(1));
      end
    end
    
  end  % state variable types
end    % file per 4D-Var cycle

return
