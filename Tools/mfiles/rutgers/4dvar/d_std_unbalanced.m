%
%  D_STD_UNBALANCED:  Driver script to compute and write 4D-Var unbalaced
%                     standard deviations.
%
%  This a user modifiable script that can be used to prepare ROMS 4D-Var
%  standard deviations when the balance operator is activated. They are
%  used to convert modeled error correlations to error covariances.
%
%  The first step is to run the model application for a period that is
%  long enough to compute meaningful circulation statistics, like mean
%  and standard deviations for all prognostic state variables: zeta, u,
%  v, temp, and salt.
%
%  As an example, we compute the 4D-Var standard deviations for the WC13
%  application from daily history files from 1/1/2000 to 12/25/2004.
%

% svn $Id: d_std_unbalanced.m 996 2020-01-10 04:28:56Z arango $
%=========================================================================%
%  Copyright (c) 2002-2020 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%

%--------------------------------------------------------------------------
%  User tunable parameters.
%--------------------------------------------------------------------------

%  Set standard deviation NetCDF file. The file name is edited and the
%  month will be appended as *iu_jan.nc:

STDfile = 'wc13_std_iu.nc';

mstr = {'jan', 'feb', 'mar', 'apr', 'may', 'jun', ...
        'jul', 'aug', 'sep', 'oct', 'nov', 'dec'};

%  Set input grif file.

my_root = '/home/arango/ocean/toms/repository/test';

GRDfile = strcat(my_root, '/WC13/Data/wc13_grd.nc');

%  Set input history files (string cell structure).

my_root = '/home/arango/ocean/toms/repository/test';

HISdir  = strcat(my_root, '/WC13/STD/Data');

HISfile = dir(fullfile(HISdir, 'wc*.nc'));

nfiles = length(HISfile);

%  Initialize working structure.

S.title = 'California Current System, 1/3 degree resolution (WC13)';

S.grd_file = GRDfile;               % aplication grid

S.do_zeta  = true;                  % free-surface
S.do_ubar  = true;                  % vertically integrated u-momentum
S.do_vbar  = true;                  % vertically integrated v-momentum
S.do_u     = true;                  % u-momentum
S.do_v     = true;                  % v-momentum
S.do_temp  = true;                  % temperature
S.do_salt  = true;                  % salinity

%  Initialize balance operator parameters (B structure):

B.Gname = GRDfile;      % application's grid

%  Set switch for the computation of balanced, baroclinic free-surface
%  using an elliptic equation (B.elliptic=true) or by integrating the
%  hydrostatic equation (B.elliptic=false). If false, the integration
%  is from bottom to surface (B.LNM_depth=0) or from level of no
%  motion to the surface (B.LNM_depth>0).

% B.elliptic = false;   % integrate hydrostatic equation
  B.elliptic = true;    % solve SSH elliptic equation

  if (B.elliptic),
%   B.Niter = 300;      % Number of elliptical solver iterations
  end                   % (default 200)
  
  if (~B.elliptic),
    B.LNM_depth = 0;    % integrate from bottom to surface
%   B.LNM_depth = 500;  % integrate from z=-500 to surface
  end    
    
%  Internal parameters for the computation of the balanced salinity
%  in terms of the temperature. These parameters are set in
%  "ini_balance.m"

% B.ml_depth=150;       % mixed-layer depth (m, positive)
%                       % (default 100)

% B.dTdz-min=0.0001;    % minimum dT/dz allowed (Celsius/m)
%                       % (default 0.001)

%--------------------------------------------------------------------------
%  Compute monthly averages and standard deviations.
%--------------------------------------------------------------------------

%  Set state variables dynamical fields (cell array) to process.

field_list = {'zeta', 'ubar', 'vbar', 'u', 'v', 'temp', 'salt'};

%  Set grid variables dynamical fields (cell array) to write in
%  output file(s).

grid_list  = {'theta_s', 'theta_b', 'Tcline' , 'hc'     , 's_rho'  , ...
              's_w'    , 'Cs_r'   , 'Cs_w'   , 'h'      , 'lon_rho', ...
              'lat_rho', 'lon_u'  , 'lat_u'  , 'lon_v'  , 'lat_v'};

%  Extract first history file name.

HisFile1 = fullfile(HISdir, HISfile(1).name);

%  Set NetCDF file creation parameters.

V = nc_vnames(HisFile1);
nvars = length(V.Variables);
 
S.curvilinear = false;
S.masking     = false;
S.Vtransform  = 1;
S.Vstretching = 1;

for n=1:nvars,
  name=char(V.Variables(n).Name);
  switch name
    case 'spherical'
      S.spherical = nc_read(HisFile1, 'spherical');
      if (ischar(S.spherical)),
        if (S.spherical == 'T' || S.spherical == 't');
          S.spherical = 1;
        else
          S.spherical = 0;
        end
      end
    case 'Vtransform'
      S.Vtransform  = nc_read(HisFile1, 'Vtransform');
    case 'Vstretching'
      S.Vstretching = nc_read(HisFile1, 'Vstretching');
    case 'angle'
      S.angle = nc_read(HisFile1, 'angle');
      S.curvilinear = true;
    case 'mask_rho'
      S.rmask = nc_read(HisFile1, 'mask_rho');
      S.umask = nc_read(HisFile1, 'mask_u');
      S.vmask = nc_read(HisFile1, 'mask_v');
      S.masking = true;
  end
end

if (S.curvilinear),
  grid_list = [grid_list, 'angle'];
end
  
if (S.masking),
  grid_list = [grid_list, 'mask_rho', 'mask_u', 'mask_v'];
end

%  Get grid size.

[Lr,Mr,Nr] = size(nc_read(HisFile1, 'temp', 1));

S.Lm = Lr - 2;                      % number of interior RHO x-points
S.Mm = Mr - 2;                      % number of interior RHO y-points
S.N  = Nr;                          % number of vertical RHO levels

%  Read in grid.

S.rlon = nc_read(HisFile1, 'lon_rho');
S.rlat = nc_read(HisFile1, 'lat_rho');

S.ulon = nc_read(HisFile1, 'lon_u');
S.ulat = nc_read(HisFile1, 'lat_u');

S.vlon = nc_read(HisFile1, 'lon_v');
S.vlat = nc_read(HisFile1, 'lat_v');

%  Process month-by-month.

for m=1:12,

%  Set standard deviation file name.

  lstr = length(STDfile)-3;
  S.ncname = char(strcat(STDfile(1:lstr), '_', mstr(m), '.nc'));
  
%  Initialize mean and variance arrays.

  Rcount = 0;                         % record counter
  rec    = 1;                         % initialization record

  S.month = m;

  for fval = field_list,
    field     = char(fval);
    field_avg = [field, '_avg'];

    try
      S.(field_avg) = zeros(size(nc_read(HisFile1, field, rec))); 
    catch
      disp([' D_STD: error while processing, rec = ', num2str(rec)]);
      disp(['        for variable : ', field]);
      disp(['        in file: ', HisFile1]);
      return
    end
  end
 
  disp(' ');
  disp([ 'Computing mean fields, month = ', num2str(m), ' ...']);
  disp(' ');

%  Accumulate montly fields.

  for n=1:nfiles,

    ncfile = fullfile(HISdir, HISfile(n).name);

    time = nc_read(ncfile,'ocean_time');
    Nrec = length(time);

    for rec=1:Nrec,

      [year,month,day]=datevec(datenum(1968,5,23) + time(rec)/86400);
      
      if (month == m),

        mydate=datestr(datenum(1968,5,23) + time(rec)/86400, 0);
        disp([ '*** Processing Averages: ', mydate]);

        Rcount = Rcount + 1;
	
        for fval = field_list,
          field     = char(fval);           % convert cell to string
          field_avg = [field, '_avg'];      % average field 

          try
            F = nc_read(ncfile, field, rec);

          catch
            disp([' D_STD: error while processing, rec = ', num2str(rec)]);
            disp(['        for variable : ', field]);
            disp(['        in file: ', HisFile1]);
            return
          end,

          S.(field_avg) = S.(field_avg) + F;
        end
    
      end
    
    end
  
  end

%  Compute monthly mean fields.

  for fval = field_list,
    field     = char(fval);                 % convert cell to string
    field_avg = [field, '_avg'];            % average field 

    S.(field_avg) = S.(field_avg) ./ Rcount;
  end

%--------------------------------------------------------------------------
%  Compute monthly unbalanced error covariance standard deviations.
%--------------------------------------------------------------------------

  disp(' ');
  disp(['   Computing unbalanced standard deviation, month = ',         ...
        num2str(m), ' ...']);
  disp(' ');

%  Initialize balance operator structure, A.  See 'balance_driver.m' for
%  discussion about some of these parameters. The default values for the
%  parameters are set in 'ini_balance.m', the user can overwrite such
%  values here.

  if (exist('A', 'var')),
    clear('A');         % clear structure A
  end

%  Load balance operator parameters (stored in structure B).

  A = B;

%  Process monthly fields.

  Nvar = 0;             % initialize variance counter

  for n=1:nfiles,

    ncfile = fullfile(HISdir, HISfile(n).name);

    time = nc_read(ncfile,'ocean_time');
    Nrec = length(time);

    for rec=1:Nrec,

      [year,month,day]=datevec(datenum(1968,5,23)+time(rec)/86400);
      
      if (month == m),

        mydate=datestr(datenum(1968,5,23) + time(rec)/86400, 0);
        disp([ '*** Processing Variance: ', mydate]);

%  Set time record to use in the computation of the thermal expansion
%  and saline contraction coefficients in "ini_balance.m", which are
%  used to compute the balanced salinity in "s_balance.m".

        A.Hname = ncfile;               % history file to process

        A.HisTimeRec = rec;             % needed in balance_4dvar

%  Get basic state to use in the balance operator. Read in selected record
%  from History NetCDF file, "A.Hname". Only temperature and salinity are
%  needed in "balance_4dbar.m".

        for fval = field_list,
          field = char(fval);           % convert cell to string

          try
            A.(field) = nc_read(A.Hname, field, A.HisTimeRec);
          catch
            disp([' D_STD_UNBALANCED: error while processing, rec = ',  ...
                  num2str(rec)]);
            disp(['        for variable : ', field]);
            disp(['        in file: ', HisFile1]);
            return
          end
        
        end

%  Compute state anomalies from computed time mean.  Only temperature
%  and salinity anomalies (A.temp_ano, A.salt_ano) are used in 
%  "balanced_4dvar.m". The temperature anomaly is used to stablish
%  the balance part of the other state variables.

        for fval = field_list,
          field     = char(fval);              % convert cell to string
          field_ano = [field, '_ano'];         % anomaly field
          field_avg = [field, '_avg'];         % average field

          A.(field_ano) = A.(field) - S.(field_avg);
        end
 
%  Set first guess free-surface for elliptic equation.  It is only used
%  when A.elliptic = 1. Use basic state values.

        A.zeta_guess = A.zeta;                 % needed in balance_4dvar

%  Compute balanced state anomalies.

        [K] = balance_4dvar(A);
 
%  Compute unbalanced state by substracting the balanced term to the
%  anomaly field. Notice that the balaced operator is K^(-1). The minus
%  sign below accounts for the inverse operator.
%
%  Also notice that the balanced temperature anomaly (K.temp_bal) is set
%  to zero in "balance_4dvar" since temperature is the starting point.
%  This facilitates the dynamical field use in the "A" structure below
%  and yields A.temp_unb = A.temp_ano.


        for fval = field_list,
          field     = char(fval);              % convert cell to string
          field_ano = [field, '_ano'];         % anomaly field
          field_bal = [field, '_bal'];         % balanced field
          field_unb = [field, '_unb'];         % unbalanced field

          A.(field_unb) = A.(field_ano) - K.(field_bal);
        end

%  Accumulate unbalanced error covariance matrix standard deviation.

        if (Nvar == 0),                        % initialize
          for fval = field_list,
            field     = char(fval);            % convert cell to string
            field_std = [field, '_std'];       % standard deviation field

            A.(field_std) = zeros(size(A.(field)));
          end
        end

        Nvar = Nvar + 1;

        for fval = field_list,
          field     = char(fval);              % convert cell to string
          field_std = [field, '_std'];         % standard deviation field
          field_unb = [field, '_unb'];         % unbalanced field

          A.(field_std) = A.(field_std) + A.(field_unb) .^ 2;
        end

      end

    end

  end

%  Compute unbalanced error covariance matrix standard deviations.
%  Notice that we are computing a unbiased estimate of the variance
%  by dividing by (Nvar-1).

  fac = 1 / max(1, Nvar-1);
  
  for fval = field_list,
    field     = char(fval);                    % convert cell to string
    field_std = [field, '_std'];               % standard deviation field

    A.(field_std) = sqrt(fac * A.(field_std));
  end

%--------------------------------------------------------------------------
%  Write out standard deviation fields.
%--------------------------------------------------------------------------

  disp(' ');
  disp(['Writting out unbalanced standard deviation, file = ', S.ncname]);
  disp(' ');

%  Create standard deviations NetCDF file.

  s = c_std(S);
  
%  Write out grid data.

  v = 'spherical';    s = nc_write(S.ncname, v, S.spherical);
  v = 'Vtransform';   s = nc_write(S.ncname, v, S.Vtransform);
  v = 'Vstretching';  s = nc_write(S.ncname, v, S.Vstretching);

  for fval = grid_list,
    field = char(fval);                     % convert cell to string

    f = nc_read(HisFile1, field);  s = nc_write(S.ncname, field, f);
  end

% Write out standard deviation data.

  rec = 1;
  
  s = nc_write(S.ncname, 'ocean_time', 0 , rec);

  for fval = field_list,
    field     = char(fval);                 % convert cell to string
    field_std = [field, '_std'];            % standard deviation field

    s = nc_write(S.ncname, field, A.(field_std), rec);
  end

% Process next month.

end

disp(' ');
disp('Done.');
disp(' ');
