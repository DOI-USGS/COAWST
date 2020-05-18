%
% BALANCE_DRIVER: This a template script showing how to compute the 4DVar
%                 balance operator. For a more realistic application see
%                 script "d_std_unbalanced".
%

% svn $Id: balance_driver.m 996 2020-01-10 04:28:56Z arango $
%===========================================================================%
%  Copyright (c) 2002-2020 The ROMS/TOMS Group                              %
%    Licensed under a MIT/X style license                                   %
%    See License_ROMS.txt                           Hernan G. Arango        %
%===========================================================================%

%---------------------------------------------------------------------------
%  User tunable parameters.
%---------------------------------------------------------------------------
%
% The user has a lot of latitude here.  In "ini_balance.m", the Grid
% NetCDF is used to read in bathymetry (h), Coriolis parameter (f),
% curvilinear metrics (pm, pn), and land/sea masking arrays, if any.
% The history file is used to compute the depths (Zr and Zw) and to
% compute the thermal expansion and saline contraction coefficients
% which are used in "s_balance.m".  When computing the depths, a
% free-surface value of zero is used since we need the rest state
% depths.

 my_root = '/home/arango/ocean/toms/repository/Projects/WC13';

 A.Gname = strcat(my_root, '/Data/wc13_grd.nc');
 A.Hname = strcat(my_root, '/Forward/wc13_his.nc');
 A.Aname = strcat(my_root, '/Forward/wc13_avg.nc');

% Set output error covariance standard deviation file.

 A.Sname = 'wc13_std_i.nc';

%  Initialize working structure.

 S.title = 'California Current System, 1/3 degree resolution (WC13)';

 S.grd_file = A.Gname;               % aplication grid

 S.ncname   = A.Sname;               % standard deviation file
 
 S.do_zeta  = true;                  % free-surface
 S.do_ubar  = true;                  % vertically integrated u-momentum
 S.do_vbar  = true;                  % vertically integrated v-momentum
 S.do_u     = true;                  % u-momentum
 S.do_v     = true;                  % v-momentum
 S.do_temp  = true;                  % temperature
 S.do_salt  = true;                  % salinity

% Initialize balance operator parameters:  It is possible to
% over-write some of the internal parameters in the computation of
% the balanced operator. These parameters are set in "ini_balance.m"
% and default values are shown in parenthesis.
%
% Set switch for the computation of balanced, baroclinic free-surface
% using an elliptic equation (B.elliptic=true) or by integrating the
% hydrostatic equation (B.elliptic=false). If false, the integration
% is from bottom to surface (B.LNM_depth=0) or from level of no
% motion to the surface (B.LNM_depth>0).

 A.elliptic = false;     % integrate hydrostatic equation
%A.elliptic = true;      % solve SSH elliptic equation

 if (A.elliptic),
%  A.Niter = 300;        % Number of elliptical solver iterations
 end,                    % (default 200)
  
 if (~A.elliptic),
%  A.LNM_depth = 0;      % integrate from bottom to surface
   A.LNM_depth = 1000;   % integrate from z=-1000 to surface
 end,   
    
% Internal parameters for the computation of the balanced salinity
% in terms of the temperature. These parameters are set in
% "ini_balance.m"

%A.ml_depth=150;         % mixed-layer depth (m, positive)
%                        % (default 100)

%A.dTdz-min=0.0001;      % minimum dT/dz allowed (Celsius/m)
%                        % (default 0.001)

%---------------------------------------------------------------------------
% Initialize several variables needed for processing.
%---------------------------------------------------------------------------

% Set state variables dynamical fields (cell array) to process.

 field_list = {'zeta', 'ubar', 'vbar', 'u', 'v', 'temp', 'salt'};

% Set grid variables dynamical fields (cell array) to write in
% output file(s).

 grid_list  = {'theta_s', 'theta_b', 'Tcline' , 'hc'     , 's_rho'  , ...
               's_w'    , 'Cs_r'   , 'Cs_w'   , 'h'      , 'lon_rho', ...
               'lat_rho', 'lon_u'  , 'lat_u'  , 'lon_v'  , 'lat_v'};

% Set NetCDF file creation parameters.

 V = nc_vnames(A.Hname);
 nvars = length(V.Variables);

 S.curvilinear = false;
 S.masking     = false;
 S.Vtransform  = 1;
 S.Vstretching = 1;

 for n=1:nvars,
   name = char(V.Variables(n).Name);
   switch name
     case 'spherical'
       S.spherical = nc_read(A.Hname, 'spherical');
       if (ischar(S.spherical)),
         if (S.spherical == 'T' | S.spherical == 't');
           S.spherical = 1;
         else,
           S.spherical = 0;
         end,
       end,
     case 'Vtransform'
       S.Vtransform  = nc_read(A.Hname, 'Vtransform');
     case 'Vstretching'
       S.Vstretching = nc_read(A.Hname, 'Vstretching');
     case 'angle'
       S.angle = nc_read(A.Hname, 'angle');
       S.curvilinear = true;
     case 'mask_rho'
       S.rmask = nc_read(A.Hname, 'mask_rho');
       S.umask = nc_read(A.Hname, 'mask_u');
       S.vmask = nc_read(A.Hname, 'mask_v');
       S.masking = true;
   end,
 end,

 if (S.curvilinear),
   grid_list = [grid_list, 'angle'];
 end,
  
 if (S.masking),
   grid_list = [grid_list, 'mask_rho', 'mask_u', 'mask_v'];
 end,

%  Get grid size.

 [Lr,Mr,Nr] = size(nc_read(A.Hname, 'temp', 1));

 S.Lm = Lr - 2;                      % number of interior RHO x-points
 S.Mm = Mr - 2;                      % number of interior RHO y-points
 S.N  = Nr;                          % number of vertical RHO levels

%  Read in grid.

 S.rlon = nc_read(A.Hname, 'lon_rho');
 S.rlat = nc_read(A.Hname, 'lat_rho');

 S.ulon = nc_read(A.Hname, 'lon_u');
 S.ulat = nc_read(A.Hname, 'lat_u');

 S.vlon = nc_read(A.Hname, 'lon_v');
 S.vlat = nc_read(A.Hname, 'lat_v');

%===========================================================================
% Compute time averaged fields.
%===========================================================================

% Get number of time records in average file. You have the choice to
% select the values of "Tstr" and "Tend", which are optional arguments
% to "average.m", to compute the desired time-average window.

 Nrec = length(nc_read(A.Aname,'ocean_time'));
 
 Tstr = 1;
 Tend = Nrec;

% We are computing the time average from a ROMS time-averaged NetCDF
% file with a lot of records.
%
% The user may want to compute this in a different way.  See template
% "d_std_unbalanced.m" for another example.

 for fval = field_list,
   field     = char(fval);                 % convert cell to string
   field_avg = [field, '_avg'];            % average field 

   A.(field_avg) = average(A.Aname, field, Tstr, Tend);
 end,

%===========================================================================
% Compute balanced error covariance standard deviation.
%===========================================================================

 disp(' ');
 disp(['   Computing unbalanced error covariance standard deviation ...']);
 disp(' ');

% Set time records to process from the History NetCDF file.

 Nrec = length(nc_read(A.Hname,'ocean_time'));

 Tstr = 1;                                % first record to process
 Tend = Nrec;                             % last  record to process

 Nvar = 0;                                % counter for variance

 for rec=Tstr:Tend,

   time = nc_read(A.Hname, 'ocean_time', rec);
   mydate=datestr(datenum(1968,5,23) + time/86400, 0);
   disp([ '** Processing: ', mydate]);
   
% Set time record to use in the computation of the thermal expansion
% and saline contraction coefficients in "ini_balance.m", which are
% used to compute the balanced salinity in "s_balance.m".

   A.HisTimeRec = rec;                    % needed in balance_4dvar
   
%---------------------------------------------------------------------------
% Get basic state to use in the balance operator. Read in selected record
% from History NetCDF file, "A.Hname".
%---------------------------------------------------------------------------

% Only temperature and salinity are needed in "balance_4dbar.m".

   for fval = field_list,
     field = char(fval);                  % convert cell to string

     try,
       A.(field) = nc_read(A.Hname, field, A.HisTimeRec);
     catch,
       disp([' BALANCE_DRIVER: error while processing, rec = ', ...
             num2str(rec)]);
       disp(['        for variable : ', field]);
       disp(['        in file: ', A.Hname]);
       return
     end,
   
   end,
 
%---------------------------------------------------------------------------
% Compute state anomalies from computed time mean. Only temperature
% and salinity anomalies (A.temp_ano, A.salt_ano) are used in 
% "balanced_4dvar.m". The temperature anomaly is used to stablish
% the balance part of the other state variables.
%---------------------------------------------------------------------------

   for fval = field_list,
     field     = char(fval);              % convert cell to string
     field_ano = [field, '_ano'];         % anomaly field
     field_avg = [field, '_avg'];         % average field

     A.(field_ano) = A.(field) - A.(field_avg);
   end,

%---------------------------------------------------------------------------
% Set first guess free-surface for elliptic equation.  It is only used
% when A.elliptic = 1.
%---------------------------------------------------------------------------
%
% Use basic state values.

   A.zeta_guess = A.zeta;                          % needed in balance_4dvar

%---------------------------------------------------------------------------
% Compute balanced state anomalies.
%---------------------------------------------------------------------------

   [K]=balance_4dvar(A);
 
%---------------------------------------------------------------------------
% Compute unbalanced state.
%---------------------------------------------------------------------------
%
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
   end,

%---------------------------------------------------------------------------
% Accumulate unbalance error covariance matrix variance.
%---------------------------------------------------------------------------

   if (rec == Tstr),                      % initialize
     for fval = field_list,
       field     = char(fval);            % convert cell to string
       field_std = [field, '_std'];       % standard deviation field

       A.(field_std) = zeros(size(A.(field)));
     end,
   end,

   Nvar = Nvar + 1;

   for fval = field_list,
     field     = char(fval);              % convert cell to string
     field_std = [field, '_std'];         % standard deviation field
     field_unb = [field, '_unb'];         % unbalanced field

     A.(field_std) = A.(field_std) + A.(field_unb) .^ 2;
   end,

 end,

%---------------------------------------------------------------------------
% Compute unbalance error covariance matrix standard deviations.
%---------------------------------------------------------------------------

%  Notice the when computing the variance we are dividing by "Nvar"
%  instead of "Nvar-1". This is fine since we don't a statistically
%  significant long record anyway.

 for fval = field_list,
   field     = char(fval);                % convert cell to string
   field_std = [field, '_std'];           % standard deviation field

   A.(field_std) = sqrt(A.(field_std) ./ Nvar);
 end,

%---------------------------------------------------------------------------
% Write out unbalanced error covariance standard devaitions to a NetCDF
% file
%---------------------------------------------------------------------------

 disp(' ');
 disp(['Writting unbalanced standard deviation, file = ',S.ncname]);
 disp(' ');

% Create standard deviation NetCDF file.

 [status]=c_std(S);  

%  Write out grid variables.

 v = 'spherical';    s = nc_write(S.ncname, v, S.spherical);
 v = 'Vtransform';   s = nc_write(S.ncname, v, S.Vtransform);
 v = 'Vstretching';  s = nc_write(S.ncname, v, S.Vstretching);

 for fval = grid_list,
   field = char(fval);                    % convert cell to string

   var = nc_read(A.Hname, field);  s = nc_write(S.ncname, field, var);
 end,

% Write out standard deviation data.

 rec = 1;
  
 s = nc_write(S.ncname, 'ocean_time', 0 , rec);

 for fval = field_list,
   field     = char(fval);                % convert cell to string
   field_std = [field, '_std'];           % standard deviation field

   s = nc_write(S.ncname, field, A.(field_std), rec);
 end,
 
 disp(' ');
 disp('Done.');
 disp(' ');
