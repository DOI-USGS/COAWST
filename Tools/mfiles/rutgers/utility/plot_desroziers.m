function D=plot_desroziers(MODname, Tstr, Tend, dt, varargin)

%
% PLOT_DESROZIERS: Computes and plots Desroziers 4D-Var diagnostics
%
% [D]=plot_desroziers(ncname, Tstr, Tend, Tdelta, Lread, wrtPNG)
%
% It computes and plots the Desroziers et al. (2005) diagnostics for a
% given 4D-Var data assimilation cycle using the innovation, increment,
% and residual vectors. It computes the diagnosed observation and
% background variances per observation type associated with the
% control vector so it can be compared with the specified values used
% in the minimization. It is a measure of how correct or incorrect
% are the specified error hypothesis.
%
% On Input:
%
%    MODname       4D-Var DAV/MOD output NetCDF filename (string or
%                    cell array for multiple files and DA cycles)
%
%                    MODname = {'dir1/cyle1', 'dir2/cycle2', ...}
%
%    Tstr          4DVar cycle starting time (datenum)
%
%    Tend          4DVar cycle ending   time (datenum)
%
%    dt            Time binning or bucketing interval (hours, say 12)
%
%    Lread         Switch to Read or compute innovation, increment,
%                    and residual vectors (OPTIONAL, default=true)
%
%    wrtPNG        Switch to write PNG file (OPTIONAL; default false)
%
% On Input:
%
%    D(:)          Desroziers et al. (2005) diagnostics (struct array)
%
%                    D(var).SigmaB(:)    diagnosed Background  error
%                    D(var).SigmaO(:)    diagnosed Observation error
%                    D(var).SigmaBe(:)   specified Background  error
%                    D(var).SigmaOe(:)   specified Observation error
%                    D(var).meanB(:)     mean diagnosed Background error
%                    D(var).meanO(:)     mean diagnosed Observation error
%                    D(var).ratioB       mean Bck diagnosed/specified
%                    D(var).ratioO       mean Obs diagnosed/specified
%                    D(var).time(:)      averaged binning times
%                    D(var).timeSince(:) time since Tstr
%                    D(var).Nobs(:)      number of observations processed
% Example:
%
%  D=plot_desroziers(MyFile, datenum(2008,1,3), datenum(2008,1,6), 12);
%

% git $Id$
%=======================================================================%
%  Copyright (c) 2002-2025 The ROMS Group                               %
%    Licensed under a MIT/X style license                               %
%    See License_ROMS.md                            Hernan G. Arango    %
%=======================================================================%

% Initialize

switch numel(varargin)
  case 0
    Lread  = true;
    wrtPNG = false;
  case 1
    Lread  = varargin{1};
    wrtPNG = false;
  case 2
    Lread  = varargin{1};
    wrtPNG = varargin{2};
end

D = struct('ncname',  [], 'obs_var', [],                              ...
           'SigmaB',  [], 'SigmaO',  [],                              ...
           'SigmaBe', [], 'SigmaOe', [],                              ...
           'meanB',   [], 'meanO',   [],                              ...
           'meanBe',  [], 'meanOe',  [],                              ...
           'ratioB',  [], 'ratioO',  [],                              ...
           'Nobs',    [], 'time',    [], 'timeSince',  [],            ...
           'date',    []);

if iscell(MODname)
  Nfiles = length(MODname);
else
  Nfiles = 1;
  MODname = {MODname};
end

ic = 0;      % time binning or bucketing counter for D-structure

%------------------------------------------------------------------------
% If applicable, process one file at the time.
%------------------------------------------------------------------------

% Initialize data for appending multiple DA cycles.

Innovation = [];
Increment  = [];
Residual   = [];
BgError    = [];
NL_final   = [];
NL_initial = [];
ObsError   = [];
ObsLon     = [];
ObsLat     = [];
ObsProv    = [];
ObsScale   = [];
ObsTime    = [];
ObsType    = [];
ObsValue   = [];
ObsZgrid   = [];

disp(blanks(1));

foundit(1:Nfiles) = false;

for nf = 1:Nfiles
  ncname = MODname{nf};

  if exist(ncname,'file')

    foundit(nf) = true;

% Inquire about the input NetCDF file.

    I = nc_inq(ncname);

% Read observation type, provenance, and time.

    myNL_final   = nc_read(ncname, 'NLmodel_final');
    myNL_initial = nc_read(ncname, 'NLmodel_initial');
    myNL_value   = nc_read(ncname, 'NLmodel_value');
    myObsType    = nc_read(ncname, 'obs_type');
    myObsProv    = nc_read(ncname, 'obs_provenance');
    myObsLon     = nc_read(ncname, 'obs_lon');
    myObsLat     = nc_read(ncname, 'obs_lat');
    myObsError   = nc_read(ncname, 'obs_error');
    myObsScale   = nc_read(ncname, 'obs_scale');
    myObsTime    = nc_read(ncname, 'obs_time');
    myObsValue   = nc_read(ncname, 'obs_value');
    myObsZgrid   = nc_read(ncname, 'obs_Zgrid');

    if (any(strcmp({I.Variables.Name}, 'BgError_value')))
      myBgError  = nc_read(ncname, 'BgError_value');
    else
      error(['Unable to find variable: BgError_value in ',ncname]);
    end

    N = max(myObsZgrid(:));                    % number of vertical levels

    str_date{nf} = nc_getatt(ncname, 'str_date');   % 4D-Var starting date
    end_date{nf} = nc_getatt(ncname, 'end_date');   % 4D-Var ending   date

    disp(['Processing File: ', ncname]);
    disp(['       DA Cycle: ', str_date{nf}, ' to ', end_date{nf}]);

    if (wrtPNG & nf == 1)
      Fsuffix = datestr(str_date{nf}, 'yyyymmdd');
    end

% Determine time reference.

    Tattr = nc_getatt(ncname, 'units', 'obs_time');
    if (contains(Tattr, 'second'))
      myObsTime = myObsTime/86400;            % seconds to days
    end

    if (Nfiles > 1)
      Tscale = 1;                             % Time in days
    else
      Tscale = 24;                            % Time in hours
    end

    iatt = strfind(Tattr, 'since');
    if (~isempty(iatt))
      Torigin = Tattr(iatt+6:end);
      epoch   = datenum(Torigin);             % 'yyyy-mm-dd HH:MM:SS'
    end
    myObsTime = epoch + myObsTime;            % obs in date number
    myObsTime = myObsTime * Tscale;           % obs in hours

% Get innovations (observations minus background), increment
% (analysis minus background), residual (observations minus analysis).

    if (Lread)
      myInnovation = nc_read(ncname, 'innovation');
      myIncrement  = nc_read(ncname, 'increment');
      myResidual   = nc_read(ncname, 'residual');
    else
      myInnovation = myObsScale .* (myObsValue - myNL_initial);
      myIncrement  = myObsScale .* (myNL_final - myNL_initial);
      myResidual   = myObsScale .* (myObsValue - myNL_final);
    end

% Remove rejected observations.

    ind = find(myObsScale == 0);
    if ~isempty(ind)
      myInnovation(ind) = [];
      myIncrement (ind) = [];
      myResidual  (ind) = [];
      myBgError   (ind) = [];
      myNL_final  (ind) = [];
      myNL_initial(ind) = [];
      myObsError  (ind) = [];
      myObsLon    (ind) = [];
      myObsLat    (ind) = [];
      myObsProv   (ind) = [];
      myObsScale  (ind) = [];
      myObsTime   (ind) = [];
      myObsType   (ind) = [];
      myObsValue  (ind) = [];
      myObsZgrid  (ind) = [];
    end

% Remove replicated SSH observations from analysis by inquiring about
% its unique location. The strategy here is to set the SSH provenance
% of the repetitive values to a negative value

    LremoveReplicateSSH = false;

    indSSH = find(myObsType == 1);
    if ~isempty(indSSH)
      [ll,IA,IC] = unique(complex(myObsLon(indSSH), myObsLon(indSSH)));
      if ~isempty(IA)
        SSHprov = -myObsProv(indSSH);      % Set negative SSH provenance
        for n = 1:length(IA)
          SSHprov(IA(n))= -SSHprov(IA(n)); % Turn positive unique values
        end
      end
      myObsProv(indSSH) = SSHprov;         % overwrite SSH provenance
      LremoveReplicateSSH = true;
    end

    if (LremoveReplicateSSH)
      ind = find(myObsType == 1 & myObsProv < 0);
      if ~isempty(ind)
        myInnovation(ind) = [];
        myIncrement (ind) = [];
        myResidual  (ind) = [];
        myBgError   (ind) = [];
        myNL_final  (ind) = [];
        myNL_initial(ind) = [];
        myObsError  (ind) = [];
        myObsLon    (ind) = [];
        myObsLat    (ind) = [];
        myObsProv   (ind) = [];
        myObsScale  (ind) = [];
        myObsTime   (ind) = [];
        myObsType   (ind) = [];
        myObsValue  (ind) = [];
        myObsZgrid  (ind) = [];
      end
    end

% Accumulate data for all DA cycles.

    Innovation = [Innovation; myInnovation];
    Increment  = [Increment;  myIncrement];
    Residual   = [Residual;   myResidual];
    BgError    = [BgError;    myBgError];
    NL_final   = [NL_final;   myNL_final];
    NL_initial = [NL_initial; myNL_initial];
    ObsError   = [ObsError;   myObsError];
    ObsLon     = [ObsLon;     myObsLon];
    ObsLat     = [ObsLat;     myObsLat];
    ObsProv    = [ObsProv;    myObsProv];
    ObsScale   = [ObsScale;   myObsScale];
    ObsTime    = [ObsTime;    myObsTime];
    ObsType    = [ObsType;    myObsType];
    ObsValue   = [ObsValue;   myObsValue];
    ObsZgrid   = [ObsZgrid;   myObsZgrid];
  end
end

if ~any(foundit)
  for n = 1:Nfiles
    disp(['Cannot find input files: ', MODname{n}]);
  end
  return
end

%------------------------------------------------------------------------
%  Compute Desroziers et al. (2005) diagnostics.
%------------------------------------------------------------------------

% Categorize observations by their type.

otypes = unique(ObsType);
ntypes = length(otypes);

Stype = zeros(7,1);
Svars = cell(7,1);
Svars = [];

varid = 0;

for n = 1:ntypes
  switch (otypes(n))
    case 1
      if ~isempty(find(ObsType == 1))
        varid = varid + 1;
        SSHid = varid;
        Stype(varid) = 1;
        Svars{SSHid} = 'SSH';
      end
    case 4
      if ~isempty(find(ObsType == 4))
        varid = varid + 1;
        Uid = varid;
        Stype(varid) = 4;
        Svars{Uid} = 'U-velocity';
      end
    case 5
      if ~isempty(find(ObsType == 5))
        varid = varid + 1;
        Vid = varid;
        Stype(varid) = 5;
        Svars{Vid} = 'V-velocity';
      end
    case 6
      Tindex = find(ObsType == 6);
      if ~isempty(Tindex)
        varid = varid + 1;
        Tid = varid;
        Stype(varid) = 6;
        Svars{Tid} = 'All Temperature';

        if ~isempty(find(ObsZgrid(Tindex) == N))
          varid = varid + 1;
          SSTid = varid;
          Stype(varid) = 6.1;
          Svars{SSTid} = 'SST';
        end
      end
    case 7
      if ~isempty(find(ObsType == 7))
        varid = varid + 1;
        Sid = varid;
        Stype(varid) = 7;
        Svars{Sid} = 'Salinity';
      end
  end
end

Nvars  = varid;

Tstr = Tstr * Tscale;
Tend = Tend * Tscale;

if (Nfiles > 1)
  dt = dt / 24;                                % interval: hours to days
end

if (Nfiles > 1)
  Tcycle_str = datenum(str_date).* Tscale - Tstr;
end

Ntimes = length(Tstr:dt:Tend-dt)+1;

for n = 1:Nvars
  D(n).ncname  = MODname;
  D(n).obs_var = Svars{n};
end

%  Binning or bucketing 4D-Var data into "dt" intervals for diagnostics
%  computations.

Trange = [0 Tend-Tstr];

for time = Tstr:dt:Tend-dt                     % hours or days
  ic = ic + 1;
  midtime   = time + 0.5 * dt;
  timeSince = midtime - Tstr;
  midtime   = midtime / Tscale;
  my_date   = datestr(midtime, 'yyyy-mm-dd HH:MM');
  ind = find(time <= ObsTime & ObsTime < time+dt);

  for n = 1:Nvars
    D(n).SigmaB(ic)    = NaN;
    D(n).SigmaO(ic)    = NaN;
    D(n).SigmaBe(ic)   = NaN;
    D(n).SigmaOe(ic)   = NaN;
    D(n).Nobs(ic)      = 0;
    D(n).time(ic)      = midtime;
    D(n).timeSince(ic) = timeSince;
    D(n).date{ic}      = my_date;
  end

  if ~isempty(ind)
    inn   = Innovation(ind);
    inc   = Increment(ind);
    res   = Residual(ind);
    Berr  = BgError(ind);
    Oerr  = ObsError(ind);
    type  = ObsType(ind);
    Zgrid = ObsZgrid(ind);

    my_otypes = unique(type);
    my_ntypes = length(my_otypes);

    for n = 1:Nvars
      switch (Stype(n))
        case 1
          vindex = find(type == 1);
          varid = SSHid;
        case 4
          vindex = find(type == 4);
          varid = Uid;
        case 5
          vindex = find(type == 5);
          varid = Vid;
        case 6
          vindex = find(type == 6 & Zgrid < N);
          varid = Tid;
        case 6.1
          vindex = find(type == 6 & Zgrid == N);
          varid = SSTid;
        case 7
          vindex = find(type == 7);
          varid = Sid;
      end
      if ~isempty(vindex)
        my_inn  = inc(vindex);
        my_inc  = inc(vindex);
        my_res  = res(vindex);
        my_Berr = Berr(vindex);            % standard deviation
        my_Oerr = sqrt(Oerr(vindex));      % covariance: squared units!
        Nobs    = length(my_inn);

        D(varid).SigmaB(ic)  = dot(my_inc',  my_inn);
        D(varid).SigmaO(ic)  = dot(my_res',  my_inn);
        D(varid).SigmaBe(ic) = dot(my_Berr', my_Berr);
        D(varid).SigmaOe(ic) = dot(my_Oerr', my_Oerr);
        D(varid).Nobs(ic)    = Nobs;
      end
    end
  end
end                                        % process next file, if any

Nrecs = ic;

%  Compute mean variance and mean ratios between diagnosed and specified.

for n = 1:Nvars
  Mobs = sum(D(n).Nobs(:));
  D(n).meanB  = sqrt(abs(sum(D(n).SigmaB (:),'omitnan') / Mobs));
  D(n).meanBe = sqrt(abs(sum(D(n).SigmaBe(:),'omitnan') / Mobs));
  D(n).meanO  = sqrt(abs(sum(D(n).SigmaO (:),'omitnan') / Mobs));
  D(n).meanOe = sqrt(abs(sum(D(n).SigmaOe(:),'omitnan') / Mobs));

  D(n).ratioB = D(n).meanB / D(n).meanBe;
  D(n).ratioO = D(n).meanO / D(n).meanOe;
end

%  Take the square root of time-binned variances.

for n = 1:Nvars
  for ic = 1:Nrecs
     D(n).SigmaB (ic) = sqrt(abs(D(n).SigmaB (ic)/D(n).Nobs(ic)));
     D(n).SigmaBe(ic) = sqrt(abs(D(n).SigmaBe(ic)/D(n).Nobs(ic)));
     D(n).SigmaO (ic) = sqrt(abs(D(n).SigmaO (ic)/D(n).Nobs(ic)));
     D(n).SigmaOe(ic) = sqrt(abs(D(n).SigmaOe(ic)/D(n).Nobs(ic)));
  end
end

%------------------------------------------------------------------------
%  Plot Desroziers et al. (2005) diagnostics.
%------------------------------------------------------------------------

Cred1    = [255  69   0]./255;        % orange red
Cred2    = [220  20  50]./255;        % crimson red
Cred3    = [199  21 133]./255;        % medium violet red
Cred4    = [255  20 147]./255;        % deep pink
Cred5    = [255  99  71]./255;        % tomato

Corange1 = [255 140   0]./255;        % dark orange
Corange2 = [236 177  32]./255;        % orange
Csienna  = [160  82  45]./255;        % Sienna brown

Cpurple1 = [236  47 142]./255;        % purple
Cpurple2 = [128   0 128]./255;        % purple
Cpurple3 = [139   0 139]./255;        % dark magneta
Corchid  = [186  85 211]./255;        % mediun orchid

Cgreen1  = [ 34 139  34]./255;        % forest green
Cgreen2  = [119 172  45]./255;        % olive green
Cgreen3  = [107 142  35]./255;        % olive drab
Cgreen4  = [154 205  50]./255;        % yellow green
Cgreen5  = [  0 139 139]./255;        % dark cyan

Cblue1   = [ 70 130 170]./255;        % steel blue
Cblue2   = [ 30 144 255]./255;        % dodger blue
Cblue3   = [ 77 190 238]./255;        % sky blue
Cblue4   = [138  43 226]./255;        % blue violet
Cblue5   = [  0 206 209]./255;        % dark turquoise

%  Plot function parameters observations/background variance/error.

LtypeOV = '-bs';                      % Obs       line type and marker
LtypeOE = '-r^';                      % Obs Error line type and marker
MEC_OV  = Cblue2;                     % Obs       MarkerEdgeColor
MEC_OE  = Corange1;                   % Obs Error MarkerEdgeColor
MFC_OV  = Cblue2;                     % Obs       MarkerFaceColor
MFC_OE  = Corange1;                   % Obs Error MarKerFaceColor

LW      = 1.5;                        % LineWidth
MS      = 6;                          % MarkerSize

LtypeBV = '-gs';                      % Backg       line type and marker
LtypeBE = '-mv';                      % Backg Error line type and marker
MEC_BV  = Cgreen3;                    % Backg       MarkerEdgeColor
MEC_BE  = Corchid;                    % Backg Error MarkerEdgeColor
MFC_BV  = Cgreen3;                    % Backg       MarKerFaceColor
MFC_BE  = Corchid;                    % Backg Error MarkerFaceColor

if (Nfiles > 1)
  myXlabel = ['Days since ', datestr(Tstr/Tscale, 'dd-mmm-yyyy')];
else
  myXlabel = ['Hours since ', datestr(Tstr/Tscale, 'dd-mmm-yyyy')];
end

figure;

set(gcf, 'Units', 'Normalized',                                       ...
         'Position', [0.2 0.1 0.6 0.8],                               ...
         'PaperOrientation', 'landscape',                             ...
         'PaperUnits', 'Normalized',                                  ...
         'PaperPosition', [0.2 0.1 0.6 0.8]);

m = 0;
for n = 1:Nvars
  switch (Stype(n))
    case 1
      m = m + 1;
      ax11 = subplot(Nvars,2,m);
      h11a = plot (D(SSHid).timeSince, D(SSHid).SigmaO, LtypeOV,      ...
                   'LineWidth', LW,                                   ...
                   'MarkerSize', MS,                                  ...
                   'MarkerEdgeColor',MEC_OV,                          ...
                   'MarkerFaceColor', MFC_OV);
      hold on;
      h11b = plot (D(SSHid).timeSince, D(SSHid).SigmaOe, LtypeOE,     ...
                   'LineWidth', LW,                                   ...
                   'MarkerSize', MS,                                  ...
                   'MarkerEdgeColor', MEC_OE,                         ...
                   'MarkerFaceColor', MFC_OE);
      if (Nfiles > 1 & dt < 1)
        for l=1:length(Tcycle_str)
          plot([Tcycle_str(l) Tcycle_str(l)] , get(gca, 'Ylim'), 'k:');
        end
      end
      xlim(Trange);
      if (n == Nvars)
        xlabel(myXlabel, 'fontweight', 'bold');
      end
      if (n == 1)
        hl1 = legend({'$\sigma^{O}$', '$\sigma^{O}_{error}$'},        ...
                     'Interpreter', 'latex', 'FontSize', 14);
        title(hl1, {'$\Leftarrow$ Observation'},                      ...
            'Interpreter', 'latex');
      end
      hold off;
      title([D(SSHid).obs_var, ':', blanks(2),                        ...
             'mean = ', num2str(D(SSHid).meanO),                      ...
             ', ratio = ', num2str(D(SSHid).ratioO)]);

      m = m + 1;
      ax12 = subplot(Nvars,2,m);
      h12a = plot (D(SSHid).timeSince, D(SSHid).SigmaB, LtypeBV,      ...
                   'LineWidth', LW,                                   ...
                   'MarkerSize', MS,                                  ...
                   'MarkerEdgeColor', MEC_BV,                         ...
                   'MarkerFaceColor', MFC_BV);
      hold on;
      h12b = plot (D(SSHid).timeSince, D(SSHid).SigmaBe, LtypeBE,     ...
                   'LineWidth', LW,                                   ...
                   'MarkerSize', MS,                                  ...
                   'MarkerEdgeColor', MEC_BE,                         ...
                   'MarkerFaceColor', MFC_BE);
      if (Nfiles > 1 & dt < 1)
        for l=1:length(Tcycle_str)
          plot([Tcycle_str(l) Tcycle_str(l)] , get(gca, 'Ylim'), 'k:');
        end
      end
      xlim(Trange);
      if (n == Nvars)
        xlabel(myXlabel, 'fontweight', 'bold');
      end
      if (n == 1)
        hl2 = legend({'$\sigma^{B}$', '$\sigma^{B}_{error}$'},        ...
                     'Interpreter', 'latex', 'FontSize', 14);
        title(hl2, {'Background $\Rightarrow$'},                      ...
              'Interpreter', 'latex');
      end
      hold off;
      title([D(SSHid).obs_var, ':', blanks(2),                        ...
             'mean = ', num2str(D(SSHid).meanB),                      ...
             ', ratio = ', num2str(D(SSHid).ratioB)]);
    case 4
      m = m + 1;
      ax41 = subplot(Nvars,2,m);
      h41a = plot (D(Uid).timeSince, D(Uid).SigmaO, LtypeOV,          ...
                   'LineWidth', LW,                                   ...
                   'MarkerSize', MS,                                  ...
                   'MarkerEdgeColor', MEC_OV,                         ...
                   'MarkerFaceColor', MFC_OV);
      hold on;
      h41b = plot (D(Uid).timeSince, D(Uid).SigmaOe, LtypeOE,         ...
                   'LineWidth', LW,                                   ...
                   'MarkerSize', MS,                                  ...
                   'MarkerEdgeColor', MEC_OE,                         ...
                   'MarkerFaceColor', MFC_OE);
      if (Nfiles > 1 & dt < 1)
        for l=1:length(Tcycle_str)
          plot([Tcycle_str(l) Tcycle_str(l)] , get(gca, 'Ylim'), 'k:');
        end
      end
      xlim(Trange);
      if (n == Nvars)
        xlabel(myXlabel, 'fontweight', 'bold');
      end
      if (n == 1)
        hl1 = legend({'$\sigma^{O}$', '$\sigma^{O}_{error}$'},        ...
                     'Interpreter', 'latex', 'FontSize', 14);
        title(hl1, {'$\Leftarrow$ Observation'},                      ...
            'Interpreter', 'latex');
      end
      hold off;
      title([D(Uid).obs_var, ':', blanks(2),                          ...
             'mean = ', num2str(D(Uid).meanO),                        ...
             ', ratio = ', num2str(D(Uid).ratioO)]);

      m = m + 1;
      ax42 = subplot(Nvars,2,m);
      h42a = plot (D(Uid).timeSince, D(Uid).SigmaB, LtypeBV,          ...
                   'LineWidth', LW,                                   ...
                   'MarkerSize', MS,                                  ...
                   'MarkerEdgeColor', MEC_BV,                         ...
                   'MarkerFaceColor', MFC_BV);
      hold on;
      h42b = plot (D(Uid).timeSince, D(Uid).SigmaBe, LtypeBE,         ...
                   'LineWidth', LW,                                   ...
                   'MarkerSize', MS,                                  ...
                   'MarkerEdgeColor', MEC_BE,                         ...
                   'MarkerFaceColor', MFC_BE);
      if (Nfiles > 1 & dt < 1)
        for l=1:length(Tcycle_str)
          plot([Tcycle_str(l) Tcycle_str(l)] , get(gca, 'Ylim'), 'k:');
        end
      end
      xlim(Trange);
      if (n == Nvars)
        xlabel(myXlabel, 'fontweight', 'bold');
      end
      if (n == 1)
        hl2 = legend({'$\sigma^{B}$', '$\sigma^{B}_{error}$'},        ...
                     'Interpreter', 'latex', 'FontSize', 14);
        title(hl2, {'Background $\Rightarrow$'},                      ...
              'Interpreter', 'latex');
      end
      hold off;
      title([D(Uid).obs_var, ':', blanks(2),                          ...
             'mean = ', num2str(D(Uid).meanB),                        ...
             ', ratio = ', num2str(D(Uid).ratioB)]);
    case 5
      m = m + 1;
      ax51 = subplot(Nvars,2,m);
      h51a = plot (D(Vid).timeSince, D(Vid).SigmaO, LtypeOV,          ...
                   'LineWidth', LW,                                   ...
                   'MarkerSize', MS,                                  ...
                   'MarkerEdgeColor', MEC_OV,                         ...
                   'MarkerFaceColor', MFC_OV);
      hold on;
      h51b = plot (D(Vid).timeSince, D(Vid).SigmaOe, LtypeOE,         ...
                   'LineWidth', LW,                                   ...
                   'MarkerSize', MS,                                  ...
                   'MarkerEdgeColor', MEC_OE,                         ...
                   'MarkerFaceColor', MFC_OE);
      if (Nfiles > 1 & dt < 1)
        for l=1:length(Tcycle_str)
          plot([Tcycle_str(l) Tcycle_str(l)] , get(gca, 'Ylim'), 'k:');
        end
      end
      xlim(Trange);
      if (n == Nvars)
        xlabel(myXlabel, 'fontweight', 'bold');
      end
      if (n == 1)
        hl1 = legend({'$\sigma^{O}$', '$\sigma^{O}_{error}$'},        ...
                     'Interpreter', 'latex', 'FontSize', 14);
        title(hl1, {'$\Leftarrow$ Observation'},                      ...
            'Interpreter', 'latex');
      end
      hold off;
      title([D(Vid).obs_var, ':', blanks(2),                          ...
             'mean = ', num2str(D(Vid).meanO),                        ...
             ', ratio = ', num2str(D(Vid).ratioO)]);

      m = m + 1;
      ax52 = subplot(Nvars,2,m);
      h52a = plot (D(Vid).timeSince, D(Vid).SigmaB, LtypeBV,          ...
                   'LineWidth', LW,                                   ...
                   'MarkerSize', MS,                                  ...
                   'MarkerEdgeColor', MEC_BV,                         ...
                   'MarkerFaceColor', MFC_BV);
      hold on;
      h52b = plot (D(Vid).timeSince, D(Vid).SigmaBe, LtypeBE,         ...
                   'LineWidth', LW,                                   ...
                   'MarkerSize', MS,                                  ...
                   'MarkerEdgeColor', MEC_BE,                         ...
                   'MarkerFaceColor', MFC_BE);
      if (Nfiles > 1 & dt < 1)
        for l=1:length(Tcycle_str)
          plot([Tcycle_str(l) Tcycle_str(l)] , get(gca, 'Ylim'), 'k:');
        end
      end
      xlim(Trange);
      if (n == Nvars)
        xlabel(myXlabel, 'fontweight', 'bold');
      end
      if (n == 1)
        hl2 = legend({'$\sigma^{B}$', '$\sigma^{B}_{error}$'},        ...
                     'Interpreter', 'latex', 'FontSize', 14);
        title(hl2, {'Background $\Rightarrow$'},                      ...
              'Interpreter', 'latex');
      end
      hold off;
      title([D(Vid).obs_var, ':', blanks(2),                          ...
             'mean = ', num2str(D(Vid).meanB),                        ...
             ', ratio = ', num2str(D(Vid).ratioB)]);
    case 6
      m = m + 1;
      ax61 = subplot(Nvars,2,m);
      h61a = plot (D(Tid).timeSince, D(Tid).SigmaO, LtypeOV,          ...
                   'LineWidth', LW,                                   ...
                   'MarkerSize', MS,                                  ...
                   'MarkerEdgeColor', MEC_OV,                         ...
                   'MarkerFaceColor', MFC_OV);
      hold on;
      h61b = plot (D(Tid).timeSince, D(Tid).SigmaOe, LtypeOE,         ...
                   'LineWidth', LW,                                   ...
                   'MarkerSize', MS,                                  ...
                   'MarkerEdgeColor', MEC_OE,                         ...
                   'MarkerFaceColor', MFC_OE);
      if (Nfiles > 1 & dt < 1)
        for l=1:length(Tcycle_str)
          plot([Tcycle_str(l) Tcycle_str(l)] , get(gca, 'Ylim'), 'k:');
        end
      end
      xlim(Trange);
      if (n == Nvars)
        xlabel(myXlabel, 'fontweight', 'bold');
      end
      if (n == 1)
        hl1 = legend({'$\sigma^{O}$', '$\sigma^{O}_{error}$'},        ...
                     'Interpreter', 'latex', 'FontSize', 14);
        title(hl1, {'$\Leftarrow$ Observation'},                      ...
            'Interpreter', 'latex');
      end
      hold off;
      title([D(Tid).obs_var, ':', blanks(2),                          ...
             'mean = ', num2str(D(Tid).meanO),                        ...
             ', ratio = ', num2str(D(Tid).ratioO)]);

      m = m + 1;
      ax62 = subplot(Nvars,2,m);
      h62a = plot (D(Tid).timeSince, D(Tid).SigmaB, LtypeBV,          ...
                   'LineWidth', LW,                                   ...
                   'MarkerSize', MS,                                  ...
                   'MarkerEdgeColor', MEC_BV,                         ...
                   'MarkerFaceColor', MFC_BV);
      hold on;
      h62b = plot (D(Tid).timeSince, D(Tid).SigmaBe, LtypeBE,         ...
                   'LineWidth', LW,                                   ...
                   'MarkerSize', MS,                                  ...
                   'MarkerEdgeColor', MEC_BE,                         ...
                   'MarkerFaceColor', MFC_BE);
      if (Nfiles > 1 & dt < 1)
        for l=1:length(Tcycle_str)
          plot([Tcycle_str(l) Tcycle_str(l)] , get(gca, 'Ylim'), 'k:');
        end
      end
      xlim(Trange);
      if (n == Nvars)
        xlabel(myXlabel, 'fontweight', 'bold');
      end
      if (n == 1)
        hl2 = legend({'$\sigma^{B}$', '$\sigma^{B}_{error}$'},        ...
                     'Interpreter', 'latex', 'FontSize', 14);
        title(hl2, {'Background $\Rightarrow$'},                      ...
              'Interpreter', 'latex');
      end
      hold off;
      title([D(Tid).obs_var, ':', blanks(2),                          ...
             'mean = ', num2str(D(Tid).meanB),                        ...
             ', ratio = ', num2str(D(Tid).ratioB)]);
    case 6.1
      m = m + 1;
      ax611 = subplot(Nvars,2,m);
      h611a = plot (D(SSTid).timeSince, D(SSTid).SigmaO, LtypeOV,     ...
                    'LineWidth', LW,                                  ...
                    'MarkerSize', MS,                                 ...
                    'MarkerEdgeColor', MEC_OV,                        ...
                    'MarkerFaceColor', MFC_OV);
      hold on;
      h611b = plot (D(SSTid).timeSince, D(SSTid).SigmaOe, LtypeOE,    ...
                    'LineWidth', LW,                                  ...
                    'MarkerSize', MS,                                 ...
                    'MarkerEdgeColor', MEC_OE,                        ...
                    'MarkerFaceColor', MFC_OE);
      if (Nfiles > 1 & dt < 1)
        for l=1:length(Tcycle_str)
          plot([Tcycle_str(l) Tcycle_str(l)] , get(gca, 'Ylim'), 'k:');
        end
      end
      xlim(Trange);
      if (n == Nvars)
        xlabel(myXlabel, 'fontweight', 'bold');
      end
      if (n == 1)
        hl1 = legend({'$\sigma^{O}$', '$\sigma^{O}_{error}$'},        ...
                     'Interpreter', 'latex', 'FontSize', 14);
        title(hl1, {'$\Leftarrow$ Observation'},                      ...
            'Interpreter', 'latex');
      end
      hold off;
      title([D(SSTid).obs_var, ':', blanks(2),                        ...
             'mean = ', num2str(D(SSTid).meanO),                      ...
             ', ratio = ', num2str(D(SSTid).ratioO)]);

      m = m + 1;
      ax612 = subplot(Nvars,2,m);
      h612a = plot (D(SSTid).timeSince, D(SSTid).SigmaB, LtypeBV,     ...
                    'LineWidth', LW,                                  ...
                    'MarkerSize', MS,                                 ...
                    'MarkerEdgeColor', MEC_BV,                        ...
                    'MarkerFaceColor', MFC_BV);
      hold on;
      h612b = plot (D(SSTid).timeSince, D(SSTid).SigmaBe, LtypeBE,    ...
                    'LineWidth', LW,                                  ...
                    'MarkerSize', MS,                                 ...
                    'MarkerEdgeColor', MEC_BE,                        ...
                    'MarkerFaceColor', MFC_BE);
      if (Nfiles > 1 & dt < 1)
        for l=1:length(Tcycle_str)
          plot([Tcycle_str(l) Tcycle_str(l)] , get(gca, 'Ylim'), 'k:');
        end
      end
      xlim(Trange);
      if (n == Nvars)
        xlabel(myXlabel, 'fontweight', 'bold');
      end
      if (n == 1)
        hl2 = legend({'$\sigma^{B}$', '$\sigma^{B}_{error}$'},        ...
                     'Interpreter', 'latex', 'FontSize', 14);
        title(hl2, {'Background $\Rightarrow$'},                      ...
              'Interpreter', 'latex');
      end
      hold off;
      title([D(SSTid).obs_var, ':', blanks(2),                        ...
             'mean = ', num2str(D(SSTid).meanB),                      ...
             ', ratio = ', num2str(D(SSTid).ratioB)]);
    case 7
      m = m + 1;
      ax71 = subplot(Nvars,2,m);
      h71a = plot (D(Sid).timeSince, D(Sid).SigmaO, LtypeOV,          ...
                   'LineWidth', LW,                                   ...
                   'MarkerSize', MS,                                  ...
                   'MarkerEdgeColor', MEC_OV,                         ...
                   'MarkerFaceColor', MFC_OV);
      hold on;
      h71b = plot (D(Sid).timeSince, D(Sid).SigmaOe, LtypeOE,         ...
                   'LineWidth', LW,                                   ...
                   'MarkerSize', MS,                                  ...
                   'MarkerEdgeColor', MEC_OE,                         ...
                   'MarkerFaceColor', MFC_OE);
      if (Nfiles > 1 & dt < 1)
        for l=1:length(Tcycle_str)
          plot([Tcycle_str(l) Tcycle_str(l)] , get(gca, 'Ylim'), 'k:');
        end
      end
      xlim(Trange);
      if (n == Nvars)
        xlabel(myXlabel, 'fontweight', 'bold');
      end
      if (n == 1)
        hl1 = legend({'$\sigma^{O}$', '$\sigma^{O}_{error}$'},        ...
                     'Interpreter', 'latex', 'FontSize', 14);
        title(hl1, {'$\Leftarrow$ Observation'},                      ...
            'Interpreter', 'latex');
      end
      hold off;
      title([D(Sid).obs_var, ':', blanks(2),                          ...
             'mean = ', num2str(D(Sid).meanO),                        ...
             ', ratio = ', num2str(D(Sid).ratioO)]);

      m = m + 1;
      ax72 = subplot(Nvars,2,m);
      h72a = plot (D(Sid).timeSince, D(Sid).SigmaB, LtypeBV,          ...
                   'LineWidth', LW,                                   ...
                   'MarkerSize', MS,                                  ...
                   'MarkerEdgeColor', MEC_BV,                         ...
                   'MarkerFaceColor', MFC_BV);
      hold on;
      h72b = plot (D(Sid).timeSince, D(Sid).SigmaBe, LtypeBE,         ...
                   'LineWidth', LW,                                   ...
                   'MarkerSize', MS,                                  ...
                   'MarkerEdgeColor', MEC_BE,                         ...
                   'MarkerFaceColor', MFC_BE);
      if (Nfiles > 1 & dt < 1)
        for l=1:length(Tcycle_str)
          plot([Tcycle_str(l) Tcycle_str(l)] , get(gca, 'Ylim'), 'k:');
        end
      end
      xlim(Trange);
      if (n == Nvars)
        xlabel(myXlabel, 'fontweight', 'bold');
      end
      if (n == 1)
        hl2 = legend({'$\sigma^{B}$', '$\sigma^{B}_{error}$'},        ...
                     'Interpreter', 'latex', 'FontSize', 14);
        title(hl2, {'Background $\Rightarrow$'},                      ...
              'Interpreter', 'latex');
      end
      hold off;
      title([D(Sid).obs_var, ':', blanks(2),                          ...
             'mean = ', num2str(D(Sid).meanB),                        ...
             ', ratio = ', num2str(D(Sid).ratioB)]);
  end
end

% Move "legend" labels between the two subplots columns;

hl1.Position(1) = hl1.Position(1) + 0.09;  % Observations column legend

hl2.Position(1) = hl1.Position(1);         % Background   column legend
hl2.Position(2) = hl2.Position(2) - 0.09;  % Background   column legend

if (wrtPNG)
  png_file = strcat('desroziers_', Fsuffix, '.png');
  exportgraphics(gcf, png_file, 'resolution', 300);
end

return
