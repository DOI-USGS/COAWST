function [deltaZ_b, zeta_err]=zeta_balance(K, elliptic, r, zeta, Niter);

%
% ZETA_BALANCE:  Computes the balanced, baroclinic free-surface
%
% [deltaZ_b] = zeta_balance(K, elliptic, r, zeta, Niter);
%
% This function computes balanced, baroclinic free-surface.
%
% There are two algorithms to compute the balanced, baroclinic
% free_surface:
%
%  (1) Intergrate hydrostatic equation from bottom to surface
%      (elliptic=0, K.LNM_depth=0) or from level of no motion
%      K.LNM_depth to surface (elliptic=0, K.LNM_depth > 0).
%
%  (2) Solve an elliptic equation as in Fukumori et al. (1998)
%      and Weaver et al. (2005), elliptic=1
%
%      div (h grad(zeta)) = - div (int{int{grad(rho)z'/rho0) dz'} dz})
%
% On Input:
%
%    K           Balance operator structure array, as computed
%                  in function "ini_balance"
%
%    elliptic    Balanced, baroclinic free-surface computation switch:
%
%                  elliptic=0,  integrate hydrostatic equation
%                  elliptic=1,  solve elliptic equation
%
%    r           Either balanced density OR elliptic equation RHS:
%
%                  elliptic=0,  balance density (kg/m3) as computed
%                               from function "rho_balance" (3D array)
%
%                  elliptic=1,  RHS term for elliptic equation (m3/s2)
%                               as computed from the pressure gradient
%                               in function "uv_balance" (2D array):
%
%                                 int{int{grad(rho)z'/rho0) dz'} dz}
%
%    zeta        First guess free-surface (m), use backgound as
%                  starting value (2D array); only needed if elliptic=1
%
%    Niter       Number of biconjugate gradient iterations, OPTIONAL
%                  (default: K.Niter); only needed if elliptic=1
%
% On Output:
%
%    deltaZ_b    Balanced, baroclinic free-surface anomaly (m), 2D array
%    zeta_err    Error in elliptical solver after "K.Niter" iterations
%
% References:
%
%   Fukumori, I., R. Raghunath and L. Fu, 1998: Nature of global
%     large-scale sea level variability in relation to atmospheric
%     forcing: a modeling study, J. Geophys. Res., 103, 5493-5512.
%
%   Weaver, A.T., C. Deltel, E. Machu, S. Ricci, and N. Daget, 2005:
%     A multivariate balance operator for variational data assimilation,
%     Q. J. R. Meteorol. Soc., 131, 3605-3625.
%
% Calls:
%
%    biconj:     Biconjugate gradient algorithm
%

% svn $Id: zeta_balance.m 996 2020-01-10 04:28:56Z arango $
%===========================================================================%
%  Copyright (c) 2002-2020 The ROMS/TOMS Group                              %
%    Licensed under a MIT/X style license           Andrew M. Moore         %
%    See License_ROMS.txt                           Hernan G. Arango        %
%===========================================================================%

%  Initialize internal parameters.

Istr = 1;
Iend = K.Lm+2;
Jstr = 1;
Jend = K.Mm+2;
N    = K.N;

zeta_err = [];

%---------------------------------------------------------------------------
%  Integrate hydrostatic equation.
%---------------------------------------------------------------------------

if (~ elliptic),

%  Integrate from bottom to surface.

  if (K.LNM_depth == 0)

    cff = -1.0 / K.rho0;

    deltaZ_b(1:Iend,1:Jend) = 0;

    for k=1:N,
      deltaZ_b(Istr:Iend,Jstr:Jend) = deltaZ_b(Istr:Iend,Jstr:Jend) +    ...
                                      cff .* r(Istr:Iend,Jstr:Jend,k) .* ...
                                          K.Hz(Istr:Iend,Jstr:Jend,k) .* ...
	                               K.rmask(Istr:Iend,Jstr:Jend);
    end,

  else,
    
%  Integrate from level of no motion (K.LNM_depth > 0) to the surface
%  or integrate from local bottom if shallower than K.LNM_depth. Note
%  that the W-points depths (K.Zw) are negative.
%
%  The vertical integration is done by recomputing the level thickness
%  (Hz) from (K.Zw) and zero out its values below the level of no motion
%  (K.LNM_depth). A very elegant and efficient approach.  The density, r,
%  is assumed to be constant in grid cell containing the level of no
%  motion, as it should be in ROMS finite volume formulation.

    cff = -1.0 / K.rho0;

    Zw = K.Zw;                            % load W-points depths

    ind = find((Zw + K.LNM_depth) < 0);   % find values below LNM_depth
    if (~isempty(ind)),                   % and substitute its values
      Zw(ind) = - K.LNM_depth;            % with -LNM_depth, so we get
    end,                                  % zero vertical difference in Hz

    Hz = diff(Zw, 1, 3);                  % Recompute vertical thickness

    deltaZ_b(1:Iend,1:Jend) = 0;          % Then, integrate vertically
    
    for k=1:N,
      deltaZ_b(Istr:Iend,Jstr:Jend) = deltaZ_b(Istr:Iend,Jstr:Jend) +    ...
                                      cff .* r(Istr:Iend,Jstr:Jend,k) .* ...
                                            Hz(Istr:Iend,Jstr:Jend,k) .* ...
	                               K.rmask(Istr:Iend,Jstr:Jend);
    end,

  end,

%---------------------------------------------------------------------------
%  Solve elliptic equation: use biconjugate gradient algorithm.
%---------------------------------------------------------------------------

else,

  if (nargin < 4),
    Niter = K.Niter;
  end,
  
  [deltaZ_b, zeta_err] = biconj(K, r, zeta, Niter);

end,

return
