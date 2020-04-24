function [K]=balance_4dvar(A);

%
% BALANCE:  Computes ROMS 4DVar balance operator fields
%
% [K]=balance(A)
%
% This function computes ROMS 4DVar balance operator fields that may
% be used to compute the error covariance matrix, B. The balance
% operator imposed a multivariate constraint on the error covariance
% such that the unobserved variables information is extracted from
% observed data.
%
% If follows the approach proposed by Weaver et al. (2005). The state
% vector is split into balanced and unbalanced components, except
% for temperature, which is used to stablish the balanced part of
% other state variables.
%
% The error covariance matrix is represented as:
%
%     B = K Bu transpose(K)
%
% where
%
%     B : background/model error covariance matrix.
%     Bu: unbalanced background/model error covariance matrix modeled
%         with the generalized diffusion operator.
%     K : balance matrix operator.
%
% The multivariate formulation is obtained by  establishing  balance
% relationships with the other state variables  using  T-S empirical
% formulas, hydrostatic balance, and geostrophic balance.
%
% The intent of this function is to compute balance state so the
% correct statistics of Bu are computed using the K^(-1) operator.
%
% Reference:
%
%   Weaver, A.T., C. Deltel, E. Machu, S. Ricci, and N. Daget, 2005:
%     A multivariate balance operator for variational data assimilation,
%     Q. J. R. Meteorol. Soc., 131, 3605-3625.
%
% On Input:
%
%    A         Application structure array: If the OPTIONAL variables
%                                           are not present, the default
%                                           values that are specified in
%                                           "ini_balance.m" will be used.
%
%              A.Gname        ROMS Grid NetCdF file name (string)
%              A.Hname        ROMS History NetCDF file name (string)
%              A.elliptic     Balanced, baroclinic free-surface computation
%                               switch:
%
%                               A.elliptic=0, integrate hydrostatic equation
%                               A.elliptic=1, solve elliptic equation
%
%              A.HisTimeRec   History NetCDF time record to use in the
%                               computation of thermal expansion and 
%                               saline contraction coefficients (integer)
%              A.temp         Basic state temperature (Celsius) used to
%                               compute dT/dz in "s_balance.m" (3D array)
%              A.salt         Basic state salinity used to compute
%                               dS/dz in "s_balance.m" (3D array)
%              A.temp_ano     Given temperature (Celsius) anomaly from a
%                               time mean (temp-temp_avg), 3D array
%              A.salt_ano     Given salinity anomaly from a
%                               time mean (salt-salt_avg), 3D array
%              A.zeta_guess   Free-surface first guess for elliptic
%                               equation, 2D array. Only needed if
%                               A.elliptic=1
%              A.ml_depth     Mixed-layer depth (m, positive) used in the
%                               computation of salinity anomaly from a T-S
%                               empirical formula, see "s_balance.m"
%                               (OPTIONAL)
%              A.dTdz_min     Minimum dT/dz (Celsius/m) allowed during the
%                               computation of salinity anomaly, see
%                               "s_balance.m" (OPTIONAL)
%              A.Niter        Number of elliptical solver iterations used
%                               in the computation of balanced, baroclinic
%                               free-surface, see "biconj.m" (OPTIONAL)
% On Output:
%
%    K         Balance operator structure array:
%
%              K.Gname        ROMS Grid NetCdF file name (string)
%              K.Hname        ROMS History NetCDF file name (string)
%              K.Lm           Number of interior points in XI-direction
%              K.Mm           Number of interior points in ETA-direction
%              K.N            Number of vertical level
%              K.Niter        Number of iterations in the biconjugate
%                               gradient solver.
%              K.EW_periodic  Switch for East-West periodic boundaries
%              K.NS_periodic  Switch for North-South periodic boundaries
%              K.g            Accelerarion due to gravity (m/s2)
%              K.rho0         Mean density (Kg/m3) used when the Boussinesq
%                               approximation is inferred
%              K.ml_depth     Mixed-layer depth (100 m)
%              K.dTdz_min     Minimum dT/dz allowed (0.001 Celsius/m)
%              K.Norder       Shapiro filter order to use
%              K.Fcoef        Shapiro filter coefficients
%              K.h            Bathymetry (m) at RHO-points
%              K.f            Coriolis parameter (1/s) at RHO-points
%              K.pm           Curvilinear X-coordinate metric (1/m)
%              K.pn           Curvilinear Y-coordinate metric (1/m)
%              K.pmon_u       Curvilinear metric pm/pn at U-points
%              K.pnom_v       Curvilinear metric pn/pm at V-points
%              K.Pmask        Land/Sea mask on PSI-points
%              K.rmask        Land/Sea mask on RHO-points
%              K.umask        Land/Sea mask on U-points
%              K.vmask        Land/Sea mask on V-points
%              K.alpha        Surface thermal expansion coefficient
%                               (1/Celsius)
%              K.beta         Surface saline contraction coefficient
%                               (nondimensional)
%              K.rho          Basic state in situ density (kg/m3)
%              K.Hz           Vertical level thicknesses (m)
%              K.Zr           Depths at vertical RHO-points (m, negative)
%              K.Zw           Depths at vertical W-points   (m, negative)
%
%              ............................................................
%
%              K.ssh_ref      Basic state reference sea-surface height (m)
%              K.rho_ano      Basic state linear density anomaly (kg/m3)
%              K.temp_bal     Balanced zero temperature anomaly (Celsius)
%              K.salt_bal     Balanced salinity anomaly
%              K.u_bal        Balanced, baroclinic U-momentum anomaly (m/s)
%              K.v_bal        Balanced, baroclinic V-momentum anomaly (m/s)
%              K.ubar_bal     Balanced, barotropic U-momentum anomaly (m/s)
%              K.vbar_bal     Balanced, barotropic V-momentum anomaly (m/s)
%              K.zeta_bal     Balanced, baroclinic free-surface anomaly (m)
%
% Calls:
%
%    ini_balance    Initializes balance operator structure array, K.
%                   It sets internal parameters, reads needed grid
%                   metrics and computes several quantities.
%
%    s_balance      Given a temperature anomaly, deltaT = temp-temp_avg, 
%                   it computes balanced salinity anomaly using a T-S
%                   empirical formula.
%
%    rho_balance    Computes balanced density anomaly using a linear
%                   equation of state.
%
%    ssh_reference  Computes the balance operartor reference sea surface
%                   height.
%
%    uv_balance     Computes balanced, baroclinic U- and V-momentum
%                   anomalies (m/s) using the geostrophic balance.
%
%    uv_barotropic  Computes vertically integrated U- and V-momentum
%                   components (m/s).
%
%    zeta_balance   Computes balanced, baroclinic free-surface anomaly
%                   by solving an elliptical equation OR integrating
%                   the hydrostatic equation from surface to bottom.
%

% svn $Id: balance_4dvar.m 996 2020-01-10 04:28:56Z arango $
%===========================================================================%
%  Copyright (c) 2002-2020 The ROMS/TOMS Group                              %
%    Licensed under a MIT/X style license           Hernan G. Arango        %
%    See License_ROMS.txt                           Andrew M. Moore         %
%===========================================================================%

% Check if required variables are present in input structure array.

if (~isfield(A,'Gname')),
  disp(' ');
  error(['BALANCE_4DVAR - unable to find ROMS Grid NetCDF file name ',   ...
         'variable: A.Gname']);
  return
end,

if (~isfield(A,'Hname')),
  disp(' ');
  error(['BALANCE_4DVAR - unable to find ROMS History NetCDF file ',     ...
         'name variable: A.Hname']);
  return
end,

if (~isfield(A,'elliptic')),
  disp(' ');
  error(['BALANCE_4DVAR - unable to find elliptic equation switch: ',    ...
         'A.elliptic']);
  return
end,

if (~A.elliptic),
  if (~isfield(A,'LNM_depth')),
    disp(' ');
    error(['BALANCE_4DVAR - unable to find level of no motion flag: ',   ...
           'A.LNM_depth']);
    return
  end,
end,
  
if (~isfield(A,'HisTimeRec')),
  disp(' ');
  error(['BALANCE_4DVAR - unable to find History time record to use: ',  ...
         'A.HisTimeRec']);
  return
end,

if (~isfield(A,'temp')),
  disp(' ');
  error(['BALANCE_4DVAR - unable to find basic state temperature ',      ...
         'variable: A.temp']);
  return
end,

if (~isfield(A,'salt')),
  disp(' ');
  error(['BALANCE_4DVAR - unable to find basic state salinity ',         ...
         'variable: A.salt']);
  return
end,

if (~isfield(A,'temp_ano')),
  disp(' ');
  error(['BALANCE_4DVAR - unable to find given temperature anomaly ',    ...
         'variable: A.temp_ano']);
  return
end,

if (~isfield(A,'salt_ano')),
  disp(' ');
  error(['BALANCE_4DVAR - unable to find given salinity anomaly ',       ...
         'variable: A.salt_ano']);
  return
end,

if (~isfield(A,'zeta_guess'));
  disp(' ');
  error(['BALANCE_4DVAR - unable to find first guess free-surface ',     ...
         'variable: A.zeta_guess']);
  return
end,  

%---------------------------------------------------------------------------
% Initialize balance operator structure array: set internal parameters,
% read needed grid metrics and compute several quantities.
%---------------------------------------------------------------------------

[K]=ini_balance(A.Gname,A.Hname,A.HisTimeRec);

% Check few input OPTIONAL parameters.

if (~isfield(A,'ml_depth')),
  ml_depth = K.ml_depth;              % use default
else,
  K.ml_depth = A.ml_depth;            % over-write default
end,

if (~isfield(A,'dTdz_min')),
  dTdz_min = K.dTdz_min;              % use default
else,
  K.dTdz_min = A.dTdz_min;            % over-write default
end,

if (~isfield(A,'Niter')),
  Niter = K.Niter;                    % use default
else,
  K.Niter = A.Niter;                  % over-write default
  Niter = A.Niter;
end,

%  Set level of no motion flag, if appropriate.

if (~A.elliptic),
  K.LNM_depth = abs(A.LNM_depth);
end,

%---------------------------------------------------------------------------
% Given temperature anomaly, temp_ano = temp - temp_avg, compute balanced
% salinity anomaly using T-S empirical formula.
%---------------------------------------------------------------------------

[K.salt_bal]=s_balance(K, A.temp, A.salt, A.temp_ano, ml_depth, dTdz_min);

%---------------------------------------------------------------------------
%  Compute density anomaly (kg/m3) using a linear equation of state.
%---------------------------------------------------------------------------

[K.rho_ano]=rho_balance(K, A.temp_ano, A.salt_ano);

%---------------------------------------------------------------------------
%  Compute balanced, baroclinic U- and V-momentum anomalies (m/s) using
%  the geostrophic balance. It uses a version of ROMS pressure gradient
%  "prsgrd31.h" algorithm.
%---------------------------------------------------------------------------

[K.u_bal, K.v_bal, K.zeta_rhs]=uv_balance(K, K.rho_ano);

%---------------------------------------------------------------------------
%  Compute basic state, reference sea surface height (m), solve elliptic
%  equation.
%---------------------------------------------------------------------------

if (A.elliptic),
  [K.ssh_ref, K.ssh_err]=ssh_reference(K, K.rho, A.zeta_guess, Niter);  
end,

%---------------------------------------------------------------------------
%  Compute balanced, baroclinic free-surface anomaly (m).
%---------------------------------------------------------------------------

if (A.elliptic),

%  Solve SSH elliptic equation as in Fukumori et al. (1998) and
%  Weaver et al. (2005).

  [K.zeta_bal, K.zeta_err]=zeta_balance(K, A.elliptic, K.zeta_rhs, ...
                                        A.zeta_guess, Niter);
else,
  
%  Integrate hydrostatic equation from surface to bottom

  [K.zeta_bal]=zeta_balance(K, A.elliptic, K.rho_ano);

end,

%---------------------------------------------------------------------------
%  Set the balanced temperature to zero since it is the same as the
%  starting temperature ananomaly (A.temp_ano).  This facilitates the
%  compact computation of the unbalanced contribution elsewhere with
%  dynamical fields in the structure.
%---------------------------------------------------------------------------

K.temp_bal = zeros(size(A.temp_ano));

%---------------------------------------------------------------------------
%  Vertically integrate balanced baroclinic momentum anomalies. This
%  quantites are not used in 4D-Var but are computed here for diagnostic
%  purposes.
%---------------------------------------------------------------------------

[K.ubar_bal, K.vbar_bal]=uv_barotropic(K.u_bal, K.v_bal, K.Hz);

return
