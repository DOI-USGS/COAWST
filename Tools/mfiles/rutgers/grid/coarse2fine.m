function F = coarse2fine(Ginp,Gout,Gfactor,varargin)

%
% COARSE2FINE:  Creates a finer resolution ROMS Grid NetCDF file
%
% F = coarse2fine(Ginp,Gout,Gfactor,Imin,Imax,Jmin,Jmax)
%
% Given a coarse resolution Grid NetCDF file (Ginp), this function creates
% a finer resolution grid in the region specified by the coarser grid
% coordinates (Imin,Jmin) and (Imax,Jmax). Notice that (Imin,Jmin), and
% (Imax,Jmax) indices are in terms of the PSI-points because it actually
% define the physical boundaries of the refined grid. The grid refinement
% coefficient is specified with Gfactor.
%
% On Input:
%
%    Ginp       Input  coaser Grid NetCDF file name (character string)
%    Gout       Output finer  Grid NetCDF file name (character string)
%    Gfactor    Grid refinement factor (1,3,5,7,9,11,13,15,...)
%    Imin       Coarse grid lower-left  I-coordinate (PSI-point)
%    Imax       Coarse grid upper-right I-coordinate (PSI-point)
%    Jmin       Coarse grid lower-left  J-coordinate (PSI-point)
%    Jmax       Coarse grid upper-right J-coordinate (PSI-point)
%
% On Output:
%
%    F          Fine resolution Grid structure
%

% svn $Id: coarse2fine.m 996 2020-01-10 04:28:56Z arango $
%=========================================================================%
%  Copyright (c) 2002-2020 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%

% Check refinement factor. The values below (1; 3:2:27) are somewhat
% rediculous for an application in ROMS. However, these values are
% possible if the strategy is extract the final grids to be used
% in ROMS from an intermediate very fine resolution grid. It Also
% allows Gfactor=1 to extract a small grid from a larger one.

if 0~=mod(Gfactor-1,2) || Gfactor >= 27
  error([' COARSE2FINE: illegal refinement factor, Gfactor = ',         ...
         num2str(Gfactor)]);
end

% Get coarse grid structure.

C = get_roms_grid(Ginp);

% Set coarse grid refinement PSI-indices.

[Lp,Mp] = size(C.h);

% Set offset "half" factor which is used to extract boundary points
% outside of the PSI-points perimeter defined by Imin, Imax, Jmin,
% and Jmax. This also takes into account the exterior contact points.

half = floor(Gfactor-1)/2;

Imin = half;                   % Default PSI-points range to apply
Imax = Lp-half;                % refinement to the full input grid
Jmin = half;
Jmax = Mp-half;

switch numel(varargin)
  case 1
    Imin = varargin{1};
  case 2
    Imin = varargin{1};
    Imax = varargin{2};
  case 3
    Imin = varargin{1};
    Imax = varargin{2};
    Jmin = varargin{3};
  case 4
    Imin = varargin{1};
    Imax = varargin{2};
    Jmin = varargin{3};
    Jmax = varargin{4};
end

F.refine_factor = Gfactor;   % These values are needed in 'refined_gridvar'
F.parent_Imin = Imin;        % for the 9-points quadratic conservative
F.parent_Imax = Imax;        % interpolation
F.parent_Jmin = Jmin;
F.parent_Jmax = Jmax;

% Set spherical switch.

spherical = C.spherical;

if spherical
  F.spherical = 1;
else
  F.spherical = 0;
end

% Set uniform grid switch.

if isfield(C, 'uniform')
  F.uniform = C.uniform;
else
  F.uniform = 0;
end

% Set curvilinear switch.

curvilinear = C.curvilinear;

if curvilinear
  F.curvilinear = 1;
else
  F.curvilinear = 0;
end

% Get grid lengths.

F.xl = 0;
F.el = 0;

% Check available grid variables. Get input file information
% structure, I.

I = nc_inq(Ginp);
vnames = {I.Variables.Name};

got_list = {'lon_rho'  , 'lat_rho'  , 'lon_psi'  , 'lat_psi'  ,         ...
            'lon_u'    , 'lat_u'    , 'lon_v'    , 'lat_v'    ,         ...
            'x_rho'    , 'y_rho'    , 'x_psi'    , 'y_psi'    ,         ...
            'x_u'      , 'y_u'      , 'x_v'      , 'y_v'      ,         ...
            'hraw'     , 'angle'    ,                                   ...
            'mask_rho' , 'mask_psi' , 'mask_u'   , 'mask_v'   ,         ...
            'lon_coast', 'lat_coast'};

for value = got_list
  field = char(value);
  got.(field) =  any(strcmp(vnames, field));
end

% Set fields to process.

field_list = {'f', 'h', 'pm', 'pn'};

if curvilinear
  field_list = [field_list, 'dmde', 'dndx'];
end

if got.angle
  field_list = [field_list, 'angle'];
end

if got.x_rho && got.y_rho
  field_list = [field_list, 'x_rho', 'y_rho'];
end

if got.x_psi && got.y_psi
  field_list = [field_list, 'x_psi', 'y_psi'];
end

if got.x_u && got.y_u
  field_list = [field_list, 'x_u', 'y_u'];
end

if got.x_v && got.y_v
  field_list = [field_list, 'x_v', 'y_v'];
end

if spherical
  field_list = [field_list, 'lon_rho', 'lat_rho',                       ...
                            'lon_psi', 'lat_psi',                       ...
                            'lon_u',   'lat_u',                         ...
                            'lon_v',   'lat_v'];
end

if got.mask_rho || got.mask_psi || got.mask_u || got.mask_v
  field_list = [field_list, 'mask_rho', 'mask_psi',                     ...
                            'mask_u',   'mask_v'];
end

% Set coarser grid fractional coordinates.

[Lp,Mp] = size(C.h);          % RHO-points

L = Lp-1;
M = Mp-1;

[YrC,XrC] = meshgrid(0.5:1:Mp-0.5, 0.5:1:Lp-0.5);
[YpC,XpC] = meshgrid(1.0:1:M     , 1.0:1:L     );
[YuC,XuC] = meshgrid(0.5:1:Mp-0.5, 1.0:1:L     );
[YvC,XvC] = meshgrid(1.0:1:M     , 0.5:1:Lp-0.5);

C.Lr = Lp;
C.Mr = Mp;

% Check refinement region.

if Imin >= Imax
  error([' COARSE2FINE: Imin >= Imax,    ',                             ...
         ' Imin = ', num2str(Imin), ',',                                ...
         ' Imax = ', num2str(Imax)]);
end

if Jmin >= Jmax
  error([' COARSE2FINE: Jmin >= Jmax,    ',                             ...
         ' Jmin = ', num2str(Jmin), ','                                 ...
         ' Jmax = ', num2str(Jmax)]);
end

if Imin < 1
  error([' COARSE2FINE: Imin < 1,   ',                                  ...
         ' Imin = ', num2str(Imax)]);
end

if Imax > L
  error([' COARSE2FINE: Imax > L,    ',                                 ...
         ' Imax = ', num2str(Imax), ',',                                ...
         ' L = ', num2str(Lp-half)]);
end

if Jmin < 1
  error([' COARSE2FINE: Jmin < 1,   ',                                  ...
         ' Jmin = ', num2str(Imax)]);
end

if Jmax > M
  error([' COARSE2FINE: Jmax > M,    ',                                 ...
         ' Jmax = ', num2str(Jmax), ',',                                ...
         ' M = ', num2str(Mp-half)]);
end

% Set finer grid fractional coordinates.

delta = 1.0/Gfactor;
halfdelta  = 0.5*delta;

IpF = (Imin:delta:Imax);                           % PSI-points
JpF = (Jmin:delta:Jmax);                           % PSI-points
IrF = [IpF(1)-halfdelta IpF+halfdelta];            % RHO-points
JrF = [JpF(1)-halfdelta JpF+halfdelta];            % RHO-points

[YrF,XrF] = meshgrid(JrF,IrF);                     % RHO-points
[YpF,XpF] = meshgrid(JpF,IpF);                     % PSI-points
[YuF,XuF] = meshgrid(JrF,IpF);                     % U-points
[YvF,XvF] = meshgrid(JpF,IrF);                     % V-points

F.Lr = 0;
F.Mr = 0;

%==========================================================================
% Interpolate grid variables.
%==========================================================================

Rutgers=1; COAWST=0;
if (Rutgers)
% Check Matlab version and use "griddedInterpolant" (release 2011b), if
% available. Otherwise, use "interp2".  The "griddedInterpolant" is
% faster the interpolant can be used several times for the same C-grid
% type.  Nowadays, "interp2" calls "griddedInterpolant" and generates
% identical results but is less efficient.

Mversion = version('-release');
Vyear    = sscanf(Mversion, '%i');

if (Vyear > 2011)
  UseGriddedInterpolant = true;
else
  UseGriddedInterpolant = false;
end

% There are a couple of options for processing special grid variables:
% 'h', 'f', 'angle', 'pm', 'pn', 'dmde', and 'dndx'.  We can either
% use quadratic conservative interpolation or regular interpolation.

Qconservative = false;        % quadratic conservative interpolation
                              % off by default
%COAWST - sometimes this helps for small grid spacings to ensure dx_child*Nrefined=dx_parent etc
%Qconservative = true;        %  Use Nrefined scaling to compute pm pn dmde dndx angle f h hraw

% Set method (linear or cubic spline) for regular interpolation of grid
% coordinates (x,y) and/or (lon,lat) coordinates. The method is 'linear'
% for idealized Cartesian coordinates applications or 'cubic' for
% spherical applications.

if spherical
  method = 'spline';          % C2 continuity: first and second derivatives
else
  method = 'linear';          % C0 continuity
end

% The interpolation is done in nondimensional (XI,ETA) coordinates:
% (XrC,YrC) and (XrF,YrF) for RHO-point variables, and so on for other
% C-grid type variables.

disp(' ');
disp('Interpolating from coarse to fine ...');

% Grid locations at RHO-points.

if (got.x_rho && got.y_rho) || (got.lon_rho && got.lat_rho)

  if got.x_rho && got.y_rho
    if UseGriddedInterpolant
      RCr = griddedInterpolant(XrC, YrC, C.x_rho, method);

                                F.x_rho = RCr(XrF, YrF);
      RCr.Values = C.y_rho;     F.y_rho = RCr(XrF, YrF);
    else
      F.x_rho = interp2(XrC', YrC', C.x_rho', XrF, YrF, method);
      F.y_rho = interp2(XrC', YrC', C.y_rho', XrF, YrF, method);
    end
  end

  if got.lon_rho && got.lat_rho
    if UseGriddedInterpolant
      RCr = griddedInterpolant(XrC, YrC, C.lon_rho, method);

                                F.lon_rho = RCr(XrF, YrF);
      RCr.Values = C.lat_rho;   F.lat_rho = RCr(XrF, YrF);
    else
      F.lon_rho = interp2(XrC', YrC', C.lon_rho', XrF, YrF, method);
      F.lat_rho = interp2(XrC', YrC', C.lat_rho', XrF, YrF, method);
    end
  end

else
  
  error(' COARSE2FINE: unable to find coordinates at RHO-points')
  
end

% Grid locations at PSI-points.

if (got.x_psi && got.y_psi) || (got.lon_psi && got.lat_psi)

  if got.x_psi && got.y_psi
    if UseGriddedInterpolant
      RCp = griddedInterpolant(XpC, YpC, C.x_psi, method);

                                 F.x_psi = RCp(XpF, YpF);
      RCp.Values = C.y_psi;      F.y_psi = RCp(XpF, YpF);
    else
      F.x_psi = interp2(XpC', YpC', C.x_psi', XpF, YpF, method);
      F.y_psi = interp2(XpC', YpC', C.y_psi', XpF, YpF, method);
    end
  end

  if got.lon_psi && got.lat_psi
    if UseGriddedInterpolant
      RCp = griddedInterpolant(XpC, YpC, C.lon_psi, method);

                                 F.lon_psi = RCp(XpF, YpF);
      RCp.Values = C.lat_psi;    F.lat_psi = RCp(XpF, YpF);
    else
      F.lon_psi = interp2(XpC', YpC', C.lon_psi', XpF, YpF, method);
      F.lat_psi = interp2(XpC', YpC', C.lat_psi', XpF, YpF, method);
    end
  end

else
  
  error(' COARSE2FINE: unable to find coordinates at PSI-points')

end

% Grid locations at U-points.

if (got.x_u && got.y_u) || (got.lon_u && got.lat_u)

  if got.x_u && got.y_u
    if UseGriddedInterpolant
      RCu = griddedInterpolant(XuC, YuC, C.x_u, method);

                               F.x_u = RCu(XuF, YuF);
      RCu.Values = C.y_u(:);   F.y_u = RCu(XuF, YuF);
    else
      F.x_u = interp2(XuC', YuC', C.x_u', XuF, YuF, method);
      F.y_u = interp2(XuC', YuC', C.y_u', XuF, YuF, method);
    end
  end

  if got.lon_u && got.lat_u
    if UseGriddedInterpolant
      RCu = griddedInterpolant(XuC, YuC, C.lon_u, method);

                               F.lon_u = RCu(XuF, YuF);
      RCu.Values = C.lat_u;    F.lat_u = RCu(XuF, YuF);
    else
      F.lon_u = interp2(XuC', YuC', C.lon_u', XuF, YuF, method);
      F.lat_u = interp2(XuC', YuC', C.lat_u', XuF, YuF, method);
    end
  end

else
  
  error(' COARSE2FINE: unable to find coordinates at U-points')

end

% Grid locations at V-points.

if (got.x_v && got.y_v) || (got.lon_v && got.lat_v)

  if got.x_v && got.y_v
    if UseGriddedInterpolant
      RCv = griddedInterpolant(XvC, YvC, C.x_v, method);

                             F.x_v = RCv(XvF, YvF);
      RCv.Values = C.y_v;    F.y_v = RCv(XvF, YvF);
    else
      F.x_v = interp2(XvC', YvC', C.x_v', XvF, YvF, method);
      F.y_v = interp2(XvC', YvC', C.y_v', XvF, YvF, method);
    end
  end

  if got.lon_v && got.lat_v
    if UseGriddedInterpolant
      RCv = griddedInterpolant(XvC, YvC, C.lon_v, method);

                             F.lon_v = RCv(XvF, YvF);
      RCv.Values = C.lat_v;  F.lat_v = RCv(XvF, YvF);
    else
      F.lon_v = interp2(XvC', YvC', C.lon_v', XvF, YvF, method);
      F.lat_v = interp2(XvC', YvC', C.lat_v', XvF, YvF, method);
    end
  end

else
  
  error(' COARSE2FINE: unable to find coordinates at V-points')
  
end

else %COAWST

% Grid locations at RHO-points.

if (got.x_rho && got.y_rho),
  F.x_rho=interp2(XrC',YrC',C.x_rho',XrF',YrF','spline')';
  F.y_rho=interp2(XrC',YrC',C.y_rho',XrF',YrF','spline')';
  if (got.lon_rho && got.lat_rho),
    F.lon_rho=interp2(XrC',YrC',C.lon_rho',XrF',YrF','spline')';
    F.lat_rho=interp2(XrC',YrC',C.lat_rho',XrF',YrF','spline')';
  end
elseif (got.lon_rho && got.lat_rho),
  F.lon_rho=interp2(XrC',YrC',C.lon_rho',XrF',YrF','spline')';
  F.lat_rho=interp2(XrC',YrC',C.lat_rho',XrF',YrF','spline')';
  if (got.x_rho && got.y_rho),
    F.x_rho=interp2(XrC',YrC',C.x_rho',XrF',YrF','spline')';
    F.y_rho=interp2(XrC',YrC',C.y_rho',XrF',YrF','spline')';
  end
else
  error(' COARSE2FINE: unable to find coordinates at RHO-points')
end

% Grid locations at PSI-points.

if (got.x_psi && got.y_psi),
  F.x_psi=interp2(XpC',YpC',C.x_psi',XpF',YpF','spline')';
  F.y_psi=interp2(XpC',YpC',C.y_psi',XpF',YpF','spline')';
  if (got.lon_psi && got.lat_psi),
    F.lon_psi=interp2(XpC',YpC',C.lon_psi',XpF',YpF','spline')';
    F.lat_psi=interp2(XpC',YpC',C.lat_psi',XpF',YpF','spline')';
  end
elseif (got.lon_psi && got.lat_psi),
  F.lon_psi=interp2(XpC',YpC',C.lon_psi',XpF',YpF','spline')';
  F.lat_psi=interp2(XpC',YpC',C.lat_psi',XpF',YpF','spline')';
  if (got.x_psi && got.y_psi),
    F.x_psi=interp2(XpC',YpC',C.x_psi',XpF',YpF','spline')';
    F.y_psi=interp2(XpC',YpC',C.y_psi',XpF',YpF','spline')';
  end
else
  error(' COARSE2FINE: unable to find coordinates at PSI-points')
end

% Grid locations at U-points.

if (got.x_u && got.y_u),
  F.x_u=interp2(XuC',YuC',C.x_u',XuF',YuF','spline')';
  F.y_u=interp2(XuC',YuC',C.y_u',XuF',YuF','spline')';
  if (got.lon_u && got.lat_u),
    F.lon_u=interp2(XuC',YuC',C.lon_u',XuF',YuF','spline')';
    F.lat_u=interp2(XuC',YuC',C.lat_u',XuF',YuF','spline')';
  end
elseif (got.lon_u && got.lat_u),
  F.lon_u=interp2(XuC',YuC',C.lon_u',XuF',YuF','spline')';
  F.lat_u=interp2(XuC',YuC',C.lat_u',XuF',YuF','spline')';
  if (got.x_u && got.y_u),
    F.x_u=interp2(XuC',YuC',C.x_u',XuF',YuF','spline')';
    F.y_u=interp2(XuC',YuC',C.y_u',XuF',YuF','spline')';
  end
else
  error(' COARSE2FINE: unable to find coordinates at U-points')
end

% Grid locations at V-points.

if (got.x_v && got.y_v),
  F.x_v=interp2(XvC',YvC',C.x_v',XvF',YvF','spline')';
  F.y_v=interp2(XvC',YvC',C.y_v',XvF',YvF','spline')';
  if (got.lon_v && got.lat_v),
    F.lon_v=interp2(XvC',YvC',C.lon_v',XvF',YvF','spline')';
    F.lat_v=interp2(XvC',YvC',C.lat_v',XvF',YvF','spline')';
  end
elseif (got.lon_v && got.lat_v),
  F.lon_v=interp2(XvC',YvC',C.lon_v',XvF',YvF','spline')';
  F.lat_v=interp2(XvC',YvC',C.lat_v',XvF',YvF','spline')';
  if (got.x_v && got.y_v),
    F.x_v=interp2(XvC',YvC',C.x_v',XvF',YvF','spline')';
    F.y_v=interp2(XvC',YvC',C.y_v',XvF',YvF','spline')';
  end
else
  error(' COARSE2FINE: unable to find coordinates at V-points')
end

end % Rutgers or COAWST


% Get grid lengths.

if got.x_psi && got.y_psi
  F.xl = max(F.x_psi(:)) - min(F.x_psi(:));
  F.el = max(F.y_psi(:)) - min(F.y_psi(:));
else
  F.xl = 0;
  F.el = 0;
end

%--------------------------------------------------------------------------
% Important grid variables ('h', 'f', 'angle', 'pm', 'pn', 'dmde', and
% 'dndx') may required special treatment. We can eiter use quadratic
% conservative interpolation or regular interpolation.
%
% The quadratic conservative interopolation (Clark and Farley,1984;
% equations 30-36) is done in function 'refined_gridvar'. The interpolation
% is reversible between coarse-to-fine and fine-to-coarse:
%
%       AVG[F.field(:,:)_i , i=1,Gfactor**2] = C.field(Idg,Jdg)
%
% where (Idg,Jdg) are the indices of the coarse grid donor cell contatining
% the finer grid.
%
% Clark, T.L. and R.D. Farley, 1984:  Severe Downslope Windstorm
%   Calculations in Two and Three Spatial Dimensions Using Anelastic
%   Interative Grid Nesting: A Possible Mechanism for Gustiness,
%   J. Atmos. Sci., 329-350.
%--------------------------------------------------------------------------

% Coriolis parameter.

if Qconservative
  F.f = refined_gridvar(C, F, 'f');
else
  if UseGriddedInterpolant
    RCr.Values = C.f;
    F.f = RCr(XrF, YrF);
  else
    F.f = interp2(XrC', YrC', C.f', XrF, YrF, method);
  end
end

% Curvilinear angle.

if got.angle
  if Qconservative
    F.angle = refined_gridvar(C, F, 'angle');
  else
    if UseGriddedInterpolant
      RCr.Values = C.angle;
      F.angle = RCr(XrF, YrF);
    else
      F.angle = interp2(XrC', YrC', C.angle', XrF, YrF, method);
    end
  end
end

% Bathymetry.

hmin = min(C.h(:));

if Qconservative
  F.h = refined_gridvar(C, F, 'h');
else
  if UseGriddedInterpolant
    RCr.Values = C.h;
    F.h = RCr(XrF, YrF);
  else
    F.h = interp2(XrC', YrC', C.h', XrF, YrF, method);
  end
end

% If appropriate, clip bathymetry to donor grid minimum value.

clip_bath = false;           % off by default

if (clip_bath)
  ind = find(F.h < hmin);
  if (~isempty(ind))
    F.h(ind) = hmin;
  end
end
  
% Raw bathymetry.

if got.hraw
  try
    C.hraw = nc_read(Ginp, 'hraw', 1);
    if Qconservative
      F.hraw = refined_gridvar(C, F, 'hraw');
    else
      if UseGriddedInterpolant
        RCr.Values = C.hraw;
        F.hraw = RCr(XrF, YrF);
      else
        F.hraw = interp2(XrC', YrC', C.hraw', XrF, YrF, method);
      end
    end
  catch
    got.hraw = false;
  end
end

% Inverse grid spacing ('pm', 'pn') and curvilinear metric terms
% d(n)/d(xi) and d(m)/d(eta).

if Qconservative
  DivideCoarse = true;      % The coarse value is divided by refinement
                            % factor to guaranttee exact area conservation
  
  F.pm = refined_gridvar(C, F, 'pm', DivideCoarse);
  F.pn = refined_gridvar(C, F, 'pn', DivideCoarse);

  F.dndx = refined_gridvar(C, F, 'dndx');
  F.dmde = refined_gridvar(C, F, 'dmde');
else
  if C.uniform              % Uniform grid distribution
    dx = 1/unique(C.pm(:));
    dy = 1/unique(C.pn(:));

    F.pm = ones(size(F.h)) .* (Gfactor/dx);
    F.pn = ones(size(F.h)) .* (Gfactor/dy);
  
    F.dndx = zeros(size(F.h));
    F.dmde = zeros(size(F.h));
  else
    disp(' ');
    if spherical
      GreatCircle = true;
      disp('Computing grid spacing: great circle distances');
    else
      GreatCircle = false;
      disp('Computing grid spacing: Cartesian distances');
    end
    [F.pm, F.pn, F.dndx, F.dmde]=grid_metrics(F, GreatCircle); 
  end
end

% Land/sea masking: use nondimensional fractional coordinates. Use
% nearest point interpolation.

if (got.mask_rho || got.mask_psi || got.mask_u || got.mask_v)
  F.mask_rho = interp2(XrC', YrC', C.mask_rho', XrF, YrF, 'nearest');

  [F.mask_u, F.mask_v, F.mask_psi]=uvp_masks(F.mask_rho);
end

%--------------------------------------------------------------------------
% Create finer resolution Grid NetCDF file and write out data.
%--------------------------------------------------------------------------

% Set number of grid points.

[LpF,MpF] = size(F.h);                  % RHO-points

disp(' ');
disp(['Number of points:',                                              ...
      ' Coarse = ', num2str(Lp) ,' x ', num2str(Mp),                    ...
      ',  fine = ', num2str(LpF),' x ', num2str(MpF)]);

% Check information structure for compliance. If applicable, change
% spherical variable to integer. Modify the Land/Sea masking attributes
% for compliance. If necessary, add the "coordinates" attribute.

I = nc_check(I);

% Create ROMS Grid NetCDF file. Rewrite dimension in input file
% information structure, I.

I.Dimensions(strcmp({I.Dimensions.Name},'xi_rho' )).Length = LpF;
I.Dimensions(strcmp({I.Dimensions.Name},'eta_rho')).Length = MpF;

I.Dimensions(strcmp({I.Dimensions.Name},'xi_psi' )).Length = LpF-1;
I.Dimensions(strcmp({I.Dimensions.Name},'eta_psi')).Length = MpF-1;

I.Dimensions(strcmp({I.Dimensions.Name},'xi_u'   )).Length = LpF-1;
I.Dimensions(strcmp({I.Dimensions.Name},'eta_u'  )).Length = MpF;

I.Dimensions(strcmp({I.Dimensions.Name},'xi_v'   )).Length = LpF;
I.Dimensions(strcmp({I.Dimensions.Name},'eta_v'  )).Length = MpF-1;

mode = netcdf.getConstant('CLOBBER');
mode = bitor(mode,netcdf.getConstant('64BIT_OFFSET'));

nc_create(Gout, mode, I);

% Set global attributes.  The extracting factors global attributes names
% depend on "Gfactor" since refinement applications (Gfactor > 1) they
% are used to process the contact points.  The Gfactor=1 extracted grid
% may be used as the coaser grid in a nested application or in a non-
% nested application. It is always good idea to keep this information
% in the extrated file global attributes.

if (Gfactor == 1)
  status = nc_attadd(Gout, 'donor_grid', Ginp);
  if (status ~= 0), return, end

  status = nc_attadd(Gout, 'donor_Imin', int32(Imin));
  if (status ~= 0), return, end

  status = nc_attadd(Gout, 'donor_Imax', int32(Imax));
  if (status ~= 0), return, end

  status = nc_attadd(Gout, 'donor_Jmin', int32(Jmin));
  if (status ~= 0), return, end

  status = nc_attadd(Gout, 'donor_Jmax', int32(Jmax));
  if (status ~= 0), return, end

  status = nc_attadd(Gout, 'sampling_factor', int32(Gfactor));
  if (status ~= 0), return, end
else
  status = nc_attadd(Gout, 'parent_grid', Ginp);
  if (status ~= 0), return, end

  status = nc_attadd(Gout, 'parent_Imin', int32(Imin));
  if (status ~= 0), return, end

  status = nc_attadd(Gout, 'parent_Imax', int32(Imax));
  if (status ~= 0), return, end

  status = nc_attadd(Gout, 'parent_Jmin', int32(Jmin));
  if (status ~= 0), return, end

  status = nc_attadd(Gout, 'parent_Jmax', int32(Jmax));
  if (status ~= 0), return, end

  status = nc_attadd(Gout, 'refine_factor', int32(Gfactor));
  if (status ~= 0), return, end
end

status = nc_attdel(Gout, 'history');
if (status ~= 0), return, end

history = ['GRID file created using Matlab script ',                    ...
           which(mfilename), blanks(1), date_stamp];
status = nc_attadd(Gout, 'history', history);
if (status ~= 0), return, end

% Write out fine resolution grid variables.

disp(['Writing finer grid variables into: ', Gout]);
disp(' ');

% Write out fine resolution grid variables.

status = nc_write (Gout, 'spherical', F.spherical);
if (status ~= 0), return, end

status = nc_write (Gout, 'xl', F.xl);
if (status ~= 0), return, end

status = nc_write (Gout, 'el', F.el);
if (status ~= 0), return, end

if got.hraw
  status = nc_write (Gout, 'hraw', F.hraw, 1);
  if (status ~= 0), return, end
end

for value = field_list
  field = char(value);
  status = nc_write (Gout, field, F.(field));
  if (status ~= 0), return, end
end

% If appropriate, write coastline data.

if spherical && got.lon_coast && got.lat_coast
  add_coastline (Gout, C.lon_coast, C.lat_coast);
  if (status ~= 0), return, end
end

%--------------------------------------------------------------------------
% Get full refinement grid structure.
%--------------------------------------------------------------------------

F = get_roms_grid(Gout);

return


