function V=plot_4dvar_vectors(ncname, varargin)

%
% PLOT_4DVAR_STATS: Plots control vectors from 4D-Var data assimilation
%
% V=plot_4dvar_vectors(ncname, Lread, wrtPNG)
%
% It plots 4D-Var cycle control vectors: innovation. increment, residual,
% and background errors at the observation locations.
%
% On Input:
%
%    ncname        4D-Var DAV/MOD output NetCDF filename (string)
%
%    Lread         Switch to Read or compute innovation, increment, and
%                    residual vectors (OPTIONAL; default true, read)
%
%    wrtPNG        Switch to write out PNG files (OPTIONAL; default false)
%
% On Output:
%
%    V             4D-Var control vectors structure
%
%                    V.ncname         4D-Var DAV/MOD filename
%                    V.str_date       4D-Var starting date
%                    V.end_date       4D-Var ending   date
%                    V.observation    state vector observations
%                    V.background     prior at observation locations
%                    V.backgError     prior error at observation locations
%                    V.analysis       posterior at observation locations
%                    V.innovation     observations minus background
%                    V.increment      analysis minus background
%                    V.residual       observations minus analysis
%                    V.time           observation time (datenum)
%                    V.provenance     observation provenance
%                    V.type           observation state variable type
%

% git $Id$
%=========================================================================%
%  Copyright (c) 2002-2025 The ROMS Group                                 %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.md                            Hernan G. Arango      %
%=========================================================================%

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

V = struct('ncname',       [], 'str_date', [], 'end_date', [],         ...
           'observation',  [],                                         ...
           'background',   [],                                         ...
           'backgError',   [],                                         ...
           'analysis',     [],                                         ...
           'innovation',   [],                                         ...
           'increment',    [],                                         ...
           'residual',     [],                                         ...
           'time',         [],                                         ...
           'provenance',   [],                                         ...
           'type',         []);

% Inquire about input NetCDF file.

I = nc_inq(ncname);

% Read observation type, provenance, and time.

NL_final   = nc_read(ncname, 'NLmodel_final');
NL_initial = nc_read(ncname, 'NLmodel_initial');
NL_value   = nc_read(ncname, 'NLmodel_value');
ObsLon     = nc_read(ncname, 'obs_lon');
ObsLat     = nc_read(ncname, 'obs_lat');
ObsProv    = nc_read(ncname, 'obs_provenance');
ObsScale   = nc_read(ncname, 'obs_scale');
ObsTime    = nc_read(ncname, 'obs_time');
ObsType    = nc_read(ncname, 'obs_type');
ObsValue   = nc_read(ncname, 'obs_value');
ObsZgrid   = nc_read(ncname, 'obs_Zgrid');

if (any(strcmp({I.Variables.Name}, 'BgError_value')))
  BgError  = nc_read(ncname, 'BgError_value');
end

N = max(ObsZgrid(:));                       % number of vertical levels

str_date = nc_getatt(ncname, 'str_date');   % 4D-Var starting date
end_date = nc_getatt(ncname, 'end_date');   % 4D-Var ending   date

if (wrtPNG)
  Fsuffix = datestr(str_date, 'yyyymmdd');
end

% Determine time reference.

Tattr = nc_getatt(ncname, 'units', 'obs_time');
if (contains(Tattr, 'second'))
  ObsTime = ObsTime/86400;                  % seconds to days
end
iatt = strfind(Tattr, 'since');
if (~isempty(iatt))
  Torigin = Tattr(iatt+6:end);
  epoch   = datenum(Torigin);               % 'yyyy-mm-dd HH:MM:SS'
end
ObsTime = epoch + ObsTime;                  % obs in date number

% Get innovations (observations minus background), increment
% (analysis minus background), residual (observations minus analysis).

if (Lread)
  innovation = nc_read(ncname, 'innovation');
  increment  = nc_read(ncname, 'increment');
  residual   = nc_read(ncname, 'residual');
else
  innovation = ObsScale .* (ObsValue - NL_initial);
  increment  = ObsScale .* (NL_final - NL_initial);
  residual   = ObsScale .* (ObsValue - NL_final);
end

% Remove rejected observations.

ind = find(ObsScale == 0);
if ~isempty(ind)
  NL_initial(ind) = [];
  NL_final  (ind) = [];
  innovation(ind) = [];
  increment (ind) = [];
  residual  (ind) = [];
  ObsLon    (ind) = [];
  ObsLat    (ind) = [];
  ObsProv   (ind) = [];
  ObsScale  (ind) = [];
  ObsTime   (ind) = [];
  ObsType   (ind) = [];
  ObsValue  (ind) = [];
  ObsZgrid  (ind) = [];
  if (any(strcmp({I.Variables.Name}, 'BgError_value')))
    BgError (ind) = [];
  end
end

% Remove replicated SSH observations from analysis by inquiring its
% unique location. The strategy here is to set the SSH provenance of
% the repetitive values to a negative value

LremoveReplicateSSH = false;

indSSH = find(ObsType == 1);
if ~isempty(indSSH)
  [lonlat, IA, IC] = unique(complex(ObsLon(indSSH), ObsLon(indSSH)));
  if ~isempty(IA)
    SSHprov = -ObsProv(indSSH);         % All SSH provenance are negative
    for n = 1:length(IA)
      SSHprov(IA(n))= -SSHprov(IA(n));  % Turn positive unique SSH values
    end
  end
  ObsProv(indSSH) = SSHprov;            % overwrite SSH provenance
  LremoveReplicateSSH = true;
end

if (LremoveReplicateSSH)
  ind = find(ObsType == 1 & ObsProv < 0);
  if ~isempty(ind)
    NL_initial(ind) = [];
    NL_final  (ind) = [];
    innovation(ind) = [];
    increment (ind) = [];
    residual  (ind) = [];
    ObsLon    (ind) = [];
    ObsLat    (ind) = [];
    ObsProv   (ind) = [];
    ObsScale  (ind) = [];
    ObsTime   (ind) = [];
    ObsType   (ind) = [];
    ObsValue  (ind) = [];
    ObsZgrid  (ind) = [];
    if (any(strcmp({I.Variables.Name}, 'BgError_value')))
      BgError (ind) = [];
    end
  end
end

% Categorize observation by they type.

otypes = unique(ObsType);
ntypes = length(otypes);

Stype = zeros(7,1);
Svars = cell(7,1);
Svars = [];

varid = 0;

for n = 1:ntypes
  switch (otypes(n))
    case 1
      SSHindex = find(ObsType == 1);
      if ~isempty(SSHindex)
        varid = varid + 1;
        SSHid = varid;
        Stype(varid) = 1;
        Svars{SSHid} = 'SSH Data (x10^3)';
      end
    case 4
      Uindex = find(ObsType == 4);
      if ~isempty(Uindex)
        varid = varid + 1;
        Uid = varid;
        Stype(varid) = 4;
        Svars{Uid} = 'U-velocity Data (x10^3)';
      end
    case 5
      Vindex = find(ObsType == 5);
      if ~isempty(Vindex)
        varid = varid + 1;
        Vid = varid;
        Stype(varid) = 5;
        Svars{Vid} = 'V-velocity Data (x10^3)';
      end
    case 6
      Tindex = find(ObsType == 6);
      if ~isempty(Tindex)
        varid = varid + 1;
        Tid = varid;
        Stype(varid) = 6;
        Svars{Tid} = 'All Temperature Data (x10^3)';

        SSTindex = find(ObsZgrid(Tindex) == N);
        if ~isempty(SSTindex)
          varid = varid + 1;
          SSTid = varid;
          Stype(varid) = 6.1;
          Svars{SSTid} = 'SST Data (x10^3)';
        end
      end
    case 7
      Sindex = find(ObsType == 7);
      if ~isempty(Sindex)
        varid = varid + 1;
        Sid = varid;
        Stype(varid) = 7;
        Svars{Sid} = 'Salinity Data (x10^3)';
      end
  end
end
Nvars  = varid;

% Load output struture.

V.ncname      = ncname;
V.str_date    = str_date;
V.end_date    = end_date;
V.observation = ObsValue;
V.background  = NL_initial;
V.backgError  = BgError;
V.analysis    = NL_final;
V.innovation  = innovation;
V.increment   = increment;
V.residual    = residual;
V.time        = ObsTime;
V.provenance  = ObsProv;
V.type        = ObsType;

%--------------------------------------------------------------------------
% Plot 4D-Var cycle statistics.
%--------------------------------------------------------------------------

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

Xscale = true;                        % scale number of observation axis
scale  = 1/1000;                      % thousands of observations

% Plot innovations.

figure;

set(gcf, 'Units', 'Normalized',                                       ...
         'Position', [0.2 0.1 0.6 0.8],                               ...
         'PaperOrientation', 'landscape',                             ...
         'PaperUnits', 'Normalized',                                  ...
         'PaperPosition', [0.2 0.1 0.6 0.8]);

for n = 1:Nvars
  switch (Stype(n))
    case 1
      ax1 = subplot(Nvars,1,n);
      y = innovation(SSHindex);
      x = 1:1:length(y);
      if (Xscale), x = x.*scale; end
      bar(x, y, 'FaceColor', Corchid);
      if (Xscale), xlim([0 max(x)]); end
      xlabel([Svars{SSHid}, ':  ',                                    ...
              ' Min=',   sprintf('%12.5e',min(y(:))),                 ...
              ', Max=',  sprintf('%12.5e',max(y(:))),                 ...
              ', Mean=', sprintf('%12.5e',mean(y(:)))]);
    case 4
      ax4 = subplot(Nvars,1,n);
      y = innovation(Uindex);
      x = 1:1:length(y);
      if (Xscale), x = x.*scale; end
      bar(x, y, 'FaceColor', Cgreen3);
      if (Xscale), xlim([0 max(x)]); end
      xlabel([Svars{Uid}, ':  ',                                      ...
              ' Min=',   sprintf('%12.5e',min(y(:))),                 ...
              ', Max=',  sprintf('%12.5e',max(y(:))),                 ...
              ', Mean=', sprintf('%12.5e',mean(y(:)))]);
    case 5
      ax5 = subplot(Nvars,1,n);
      y = innovation(Vindex);
      x = 1:1:length(y);
      if (Xscale), x = x.*scale; end
      bar(x, y, 'FaceColor', Corange2);
      if (Xscale), xlim([0 max(x)]); end
      xlabel([Svars{Vid}, ':  ',                                      ...
              ' Min=',   sprintf('%12.5e',min(y(:))),                 ...
              ', Max=',  sprintf('%12.5e',max(y(:))),                 ...
              ', Mean=', sprintf('%12.5e',mean(y(:)))]);
    case 6
      ax6 = subplot(Nvars,1,n);
      y = innovation(Tindex);
      x = 1:1:length(y);
      if (Xscale), x = x.*scale; end
      bar(x, y, 'FaceColor', Cred2);
      if (Xscale), xlim([0 max(x)]); end
      xlabel([Svars{Tid}, ':  ',                                      ...
              ' Min=',   sprintf('%12.5e',min(y(:))),                 ...
              ', Max=',  sprintf('%12.5e',max(y(:))),                 ...
              ', Mean=', sprintf('%12.5e',mean(y(:)))]);
    case 6.1
      ax61 = subplot(Nvars,1,n);
      y = innovation(SSTindex);
      x = 1:1:length(y);
      if (Xscale), x = x.*scale; end
      bar(x, y, 'FaceColor', Cred4);
      if (Xscale), xlim([0 max(x)]); end
      xlabel([Svars{SSTid}, ':  ',                                    ...
              ' Min=',   sprintf('%12.5e',min(y(:))),                 ...
              ', Max=',  sprintf('%12.5e',max(y(:))),                 ...
              ', Mean=', sprintf('%12.5e',mean(y(:)))]);
    case 7
      ax7 = subplot(Nvars,1,n);
      y = innovation(Sindex);
      x = 1:1:length(y);
      if (Xscale), x = x.*scale; end
      bar(x, y, 'FaceColor', Cblue2);
      if (Xscale), xlim([0 max(x)]); end
      xlabel([Svars{Sid}, ':  ',                                      ...
              ' Min=',   sprintf('%12.5e',min(y(:))),                 ...
              ', Max=',  sprintf('%12.5e',max(y(:))),                 ...
              ', Mean=', sprintf('%12.5e',mean(y(:)))]);
  end
  if (n==1)
    title('Innovations: Observations minus Background');
  end
end

if (wrtPNG)
  png_file = strcat('inn_', Fsuffix, '.png');
  exportgraphics(gcf, png_file, 'resolution', 300);
end

% Plot increment.

figure;

set(gcf, 'Units', 'Normalized',                                       ...
         'Position', [0.2 0.1 0.6 0.8],                               ...
         'PaperOrientation', 'landscape',                             ...
         'PaperUnits', 'Normalized',                                  ...
         'PaperPosition', [0.2 0.1 0.6 0.8]);

for n = 1:Nvars
  switch (Stype(n))
    case 1
      ax1 = subplot(Nvars,1,n);
      y = increment(SSHindex);
      x = 1:1:length(y);
      if (Xscale), x = x.*scale; end
      bar(x, y, 'FaceColor', Corchid);
      if (Xscale), xlim([0 max(x)]); end
      xlabel([Svars{SSHid}, ':  ',                                    ...
              ' Min=',   sprintf('%12.5e',min(y(:))),                 ...
              ', Max=',  sprintf('%12.5e',max(y(:))),                 ...
              ', Mean=', sprintf('%12.5e',mean(y(:)))]);
    case 4
      ax4 = subplot(Nvars,1,n);
      y = increment(Uindex);
      x = 1:1:length(y);
      if (Xscale), x = x.*scale; end
      bar(x, y, 'FaceColor', Cgreen3 );
      xlabel([Svars{Uid}, ':  ',                                      ...
              ' Min=',   sprintf('%12.5e',min(y(:))),                 ...
              ', Max=',  sprintf('%12.5e',max(y(:))),                 ...
              ', Mean=', sprintf('%12.5e',mean(y(:)))]);
    case 5
      ax5 = subplot(Nvars,1,n);
      y = increment(Vindex);
      x = 1:1:length(y);
      if (Xscale), x = x.*scale; end
      bar(x, y, 'FaceColor', Corange2);
      if (Xscale), xlim([0 max(x)]); end
      xlabel([Svars{Vid}, ':  ',                                      ...
              ' Min=',   sprintf('%12.5e',min(y(:))),                 ...
              ', Max=',  sprintf('%12.5e',max(y(:))),                 ...
              ', Mean=', sprintf('%12.5e',mean(y(:)))]);
    case 6
      ax6 = subplot(Nvars,1,n);
      y = increment(Tindex);
      x = 1:1:length(y);
      if (Xscale), x = x.*scale; end
      bar(x, y, 'FaceColor', Cred2);
      xlabel([Svars{Tid}, ':  ',                                      ...
              ' Min=',   sprintf('%12.5e',min(y(:))),                 ...
              ', Max=',  sprintf('%12.5e',max(y(:))),                 ...
              ', Mean=', sprintf('%12.5e',mean(y(:)))]);
    case 6.1
      ax61 = subplot(Nvars,1,n);
      y = increment(SSTindex);
      x = 1:1:length(y);
      if (Xscale), x = x.*scale; end
      bar(x, y, 'FaceColor', Cred4);
      if (Xscale), xlim([0 max(x)]); end
      xlabel([Svars{SSTid}, ':  ',                                    ...
              ' Min=',   sprintf('%12.5e',min(y(:))),                 ...
              ', Max=',  sprintf('%12.5e',max(y(:))),                 ...
              ', Mean=', sprintf('%12.5e',mean(y(:)))]);
    case 7
      ax7 = subplot(Nvars,1,n);
      y = increment(Sindex);
      x = 1:1:length(y);
      if (Xscale), x = x.*scale; end
      bar(x, y, 'FaceColor', Cblue2);
      if (Xscale), xlim([0 max(x)]); end
      xlabel([Svars{Sid}, ':  ',                                      ...
              ' Min=',   sprintf('%12.5e',min(y(:))),                 ...
              ', Max=',  sprintf('%12.5e',max(y(:))),                 ...
              ', Mean=', sprintf('%12.5e',mean(y(:)))]);
  end
  if (n==1)
    title('Increment: Analysis minus Background');
  end
end

if (wrtPNG)
  png_file = strcat('inc_', Fsuffix, '.png');
  exportgraphics(gcf, png_file, 'resolution', 300);
end


% Plot residual.

figure;

set(gcf, 'Units', 'Normalized',                                       ...
         'Position', [0.2 0.1 0.6 0.8],                               ...
         'PaperOrientation', 'landscape',                             ...
         'PaperUnits', 'Normalized',                                  ...
         'PaperPosition', [0.2 0.1 0.6 0.8]);

for n = 1:Nvars
  switch (Stype(n))
    case 1
      ax1 = subplot(Nvars,1,n);
      y = residual(SSHindex);
      x = 1:1:length(y);
      if (Xscale), x = x.*scale; end
      bar(x, y, 'FaceColor', Corchid);
      if (Xscale), xlim([0 max(x)]); end
      xlabel([Svars{SSHid}, ':  ',                                    ...
              ' Min=',   sprintf('%12.5e',min(y(:))),                 ...
              ', Max=',  sprintf('%12.5e',max(y(:))),                 ...
              ', Mean=', sprintf('%12.5e',mean(y(:)))]);
    case 4
      ax4 = subplot(Nvars,1,n);
      y = residual(Uindex);
      x = 1:1:length(y);
      if (Xscale), x = x.*scale; end
      bar(x, y, 'FaceColor', Cgreen3);
      if (Xscale), xlim([0 max(x)]); end
      xlabel([Svars{Uid}, ':  ',                                      ...
              ' Min=',   sprintf('%12.5e',min(y(:))),                 ...
              ', Max=',  sprintf('%12.5e',max(y(:))),                 ...
              ', Mean=', sprintf('%12.5e',mean(y(:)))]);
    case 5
      ax5 = subplot(Nvars,1,n);
      y = residual(Vindex);
      x = 1:1:length(y);
      if (Xscale), x = x.*scale; end
      bar(x, y, 'FaceColor', Corange2);
      if (Xscale), xlim([0 max(x)]); end
      xlabel([Svars{Vid}, ':  ',                                      ...
              ' Min=',   sprintf('%12.5e',min(y(:))),                 ...
              ', Max=',  sprintf('%12.5e',max(y(:))),                 ...
              ', Mean=', sprintf('%12.5e',mean(y(:)))]);
    case 6
      ax6 = subplot(Nvars,1,n);
      y = residual(Tindex);
      x = 1:1:length(y);
      if (Xscale), x = x.*scale; end
      bar(x, y, 'FaceColor', Cred2);
      if (Xscale), xlim([0 max(x)]); end
      xlabel([Svars{Tid}, ':  ',                                      ...
              ' Min=',   sprintf('%12.5e',min(y(:))),                 ...
              ', Max=',  sprintf('%12.5e',max(y(:))),                 ...
              ', Mean=', sprintf('%12.5e',mean(y(:)))]);
    case 6.1
      ax61 = subplot(Nvars,1,n);
      y = residual(Tindex);
      x = 1:1:length(y);
      if (Xscale), x = x.*scale; end
      bar(x, y, 'FaceColor', Cred4);
      if (Xscale), xlim([0 max(x)]); end
      xlabel([Svars{SSTid}, ':  ',                                    ...
              ' Min=',   sprintf('%12.5e',min(y(:))),                 ...
              ', Max=',  sprintf('%12.5e',max(y(:))),                 ...
              ', Mean=', sprintf('%12.5e',mean(y(:)))]);
    case 7
      ax7 = subplot(Nvars,1,n);
      y = residual(Sindex);
      x = 1:1:length(y);
      if (Xscale), x = x.*scale; end
      bar(x, y, 'FaceColor', Cblue2);
      if (Xscale), xlim([0 max(x)]); end
      xlabel([Svars{Sid}, ':  ',                                      ...
              ' Min=',   sprintf('%12.5e',min(y(:))),                 ...
              ', Max=',  sprintf('%12.5e',max(y(:))),                 ...
              ', Mean=', sprintf('%12.5e',mean(y(:)))]);
  end
  if (n==1)
    title('Residual: Observations minus Analysis');
  end
end

if (wrtPNG)
  png_file = strcat('res_', Fsuffix, '.png');
  exportgraphics(gcf, png_file, 'resolution', 300);
end

% Plot background error standard deviation.

if (any(strcmp({I.Variables.Name}, 'BgError_value')))
  figure;

  set(gcf, 'Units', 'Normalized',                                     ...
           'Position', [0.2 0.1 0.6 0.8],                             ...
           'PaperOrientation', 'landscape',                           ...
           'PaperUnits', 'Normalized',                                ...
           'PaperPosition', [0.2 0.1 0.6 0.8]);

  ind = find(BgError == 0);
  if ~isempty(ind)
    BgError(ind) = NaN;
  end

  for n = 1:Nvars
    switch (Stype(n))
      case 1
        ax1 = subplot(Nvars,1,n);
        y = BgError(SSHindex);
        x = 1:1:length(y);
        if (Xscale), x = x.*scale; end
        bar(x, y, 'FaceColor', Corchid);
        if (Xscale), xlim([0 max(x)]); end
        xlabel([Svars{SSHid}, ':  ',                                  ...
                ' Min=',   sprintf('%12.5e',min(y(:))),               ...
                ', Max=',  sprintf('%12.5e',max(y(:))),               ...
                ', Mean=', sprintf('%12.5e',mean(y(:),'omitnan'))]);
      case 4
        ax4 = subplot(Nvars,1,n);
        y = BgError(Uindex);
        x = 1:1:length(y);
        if (Xscale), x = x.*scale; end
        bar(x, y, 'FaceColor', Cgreen3);
        if (Xscale), xlim([0 max(x)]); end
        xlabel([Svars{Uid}, ':  ',                                    ...
                ' Min=',   sprintf('%12.5e',min(y(:))),               ...
                ', Max=',  sprintf('%12.5e',max(y(:))),               ...
                ', Mean=', sprintf('%12.5e',mean(y(:),'omitnan'))]);
      case 5
        ax5 = subplot(Nvars,1,n);
        y = BgError(Vindex);
        x = 1:1:length(y);
        if (Xscale), x = x.*scale; end
        bar(x, y, 'FaceColor', Corange2);
        if (Xscale), xlim([0 max(x)]); end
        xlabel([Svars{Vid}, ':  ',                                    ...
                ' Min=',   sprintf('%12.5e',min(y(:))),               ...
                ', Max=',  sprintf('%12.5e',max(y(:))),               ...
                ', Mean=', sprintf('%12.5e',mean(y(:),'omitnan'))]);
      case 6
        ax6 = subplot(Nvars,1,n);
        y = BgError(Tindex);
        x = 1:1:length(y);
        if (Xscale), x = x.*scale; end
        bar(x, y, 'FaceColor', Cred2);
        if (Xscale), xlim([0 max(x)]); end
        xlabel([Svars{Tid}, ':   ',                                   ...
                ' Min=',   sprintf('%12.5e',min(y(:))),               ...
                ', Max=',  sprintf('%12.5e',max(y(:))),               ...
                ', Mean=', sprintf('%12.5e',mean(y(:),'omitnan'))]);
      case 6.1
        ax61 = subplot(Nvars,1,n);
        y = BgError(SSTindex);
        x = 1:1:length(y);
        if (Xscale), x = x.*scale; end
        bar(x, y, 'FaceColor', Cred4);
        if (Xscale), xlim([0 max(x)]); end
        xlabel([Svars{SSTid}, ':   ',                                 ...
                ' Min=',   sprintf('%12.5e',min(y(:))),               ...
                ', Max=',  sprintf('%12.5e',max(y(:))),               ...
                ', Mean=', sprintf('%12.5e',mean(y(:),'omitnan'))]);
      case 7
        ax7 = subplot(Nvars,1,n);
        y = BgError(Sindex);
        x = 1:1:length(y);
        if (Xscale), x = x.*scale; end
        bar(x, y, 'FaceColor', Cblue2);
        if (Xscale), xlim([0 max(x)]); end
        xlabel([Svars{Sid}, ':  ',                                    ...
                ' Min=',   sprintf('%12.5e',min(y(:))),               ...
                ', Max=',  sprintf('%12.5e',max(y(:))),               ...
                ', Mean=', sprintf('%12.5e',mean(y(:),'omitnan'))]);
    end
    if (n==1)
      title('Backround Error Standard Deviation');
    end
  end

  if (wrtPNG)
    png_file = strcat('err_', Fsuffix, '.png');
    exportgraphics(gcf, png_file, 'resolution', 300);
  end
end

return
