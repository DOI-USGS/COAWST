function F = interpolator(Ginp, field, Xpath, Ypath, ad_F, varargin)

% INTERPOLATOR: Spatially interpolate a ROMS field to horizontal
%               path/trajectory
%
% F = interpolator(Ginp, field, Xpath, Ypath, ad_F, [options])
%
% Given any 2D or 3D ROMS variable FIELD, linearly interpolate to the
% requested path (Xpath, Ypath) in spherical or Cartesian coordinates.
% If input variable is 3D, the horizontally interpolated values are
% returned at all levels of the s-coordinate.
%
% Output includes elemental cell thicknesses and path segment lengths to
% facilitate spatial integration of the interpolated field.
%
% On Input:
%
%    Ginp        ROMS NetCDF file name (string) containing coordinates
%            or, ROMS grid structure (struct array) with depth arrays so
%                3D variables can be processed (see get_roms_grid.m).
%
%    field       ROMS field to interpolate from (2D or 3D array).
%                  Avoid passing field time records. This function cannot
%                  be used for temporal interpolation.
%
%    Xpath       Longitude (degrees_east) or Cartesian (m) coordinates in
%                  the XI-direction (1D or 2D array).
%
%    Ypath       Latitude (degrees_north) or Cartesian (m) coordinates in
%                  the XI-direction (1D or 2D array).
%
%    [ ad_F ]    Adjoint array of the same size of X if 2D array. If
%                  processing 3D variables, it needs an additional
%                  dimension vertical levels. Not used here but added for
%                  compatability with associated adjoint function.
%
% On Output:
%
%    F           Interpolated field (struct array)
%
%                  F.field      variable to process (2D or 3D)
%                  F.Xgrd       variable X-locations (degrees_east or m)
%                  F.Ygrd       variable Y-locations (degrees_north or m)
%                  F.size       variable size
%                  F.Ctype      variable C-grid type
%                  F.method     interpolation method
%                  F.val        interpolated values at (X,Y) or (X,Y,k)
%                  F.Xval       interpolation X-locations (input Xpath)
%                  F.Yval       interpolation Y-locations (input Ypath)
%                  F.Lindex1    linear indices (I  ,J  ) containing (X,Y)
%                  F.Lindex2    linear indices (I+1,J  ) containing (X,Y)
%                  F.Lindex3    linear indices (I+1,J+1) containing (X,Y)
%                  F.Lindex4    linear indices (I  ,J+1) containing (X,Y)
%                  F.Jindex     lower-left field J-index containing (X,Y)
%                  F.h          bathymetry at requested (X,Y)
%                  F.dis        distance along the path (in km)
%                  F.dislen     path segment lengths for path integral
%                  F.z_w        cell center depths at requested (X,Y,N+1)
%                  F.z_r        cell face depths at requested (X,Y,N)
%                  F.dz         cell thicknesses at requested (X,Y,N)
%
%    dis, z_r, dz and dislen are returned as 2D variables so that the
%    interpolated output can be plotted with pcolor(F.dis,F.z_r,F.val)
%    and integrated with sum(F.val .* F.dislen .* F.dz)
%
% Notice that this function is adjointed, see ad_interpolator.m.  It can
% be used to set-up functionals in various adjoint based algorithms, like
% adjoint sensitivities.
%
% The output linear indices are as follows:
%
%                     Lindex4           Lindex3
%
%                        4________________3
%                        |                |
%                        |                |
%                        |<---p--->       |
%                        |                |
%                        |........o   :   |     field cell elements
%                        |       .    :   |     containing (Xpath,Ypath)
%                        |       .    q   |     as Matlab array
%                        |       .    :   |
%                        |____________:___|
%                        1                2
%
%                     Lindex1           Lindex2

% svn $Id: interpolator.m 996 2020-01-10 04:28:56Z arango $
%=========================================================================%
%  Copyright (c) 2002-2020 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license           John L. Wilkin        %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%

% The input variable ad_F is not used here but it is added
% for compatibility with the adjoint version of this function.
% If given, check it is empty so we know the user is not trying to so
% something with it

if nargin > 4
  if ~isempty(ad_F)
    warning('INTERPOLATOR does not use the ad_F input')
    error('Pass an empty [ ] argument if further inputs are used')
  end
end
if nargin > 5
  error(['There is no code in ' which(mfilename)                        ...
         ' to process additional inputs'])
end

% Interpolation method in scatteredInterpolant.

method = 'linear';

% Initialize output structure.

F = struct('field'   , [], 'Xgrd'    , [], 'Ygrd'    , [],              ...
           'size'    , [], 'Cgrid'   , [], 'method'  , [],              ...
           'val'     , [], 'Xval'    , [], 'Yval'    , [],              ...
           'Lindex1' , [], 'Lindex2' , [], 'Lindex3' , [],              ...
           'Lindex4' , [], 'h'       , [],                              ...
           'dis'     , [], 'dislen'  , [],                              ...
           'z_w'     , [], 'z_r'     , [], 'dz'      , []);

if ~isvector(Xpath)
  F = rmfield(F, 'dis');                % computed only for vector
  F = rmfield(F, 'dislen');             % trajectories
end

% Remove possible singleton (time) dimension of input variable.

field = squeeze(field);

% Set ROMS grid structure, G.

if ~isstruct(Ginp)
  G = get_roms_grid(Ginp);
else
  G = Ginp;
end

% Set ROMS grid size.

L = G.Lm+1;  Lp = L+1;
M = G.Mm+1;  Mp = M+1;

% Check input field dimensions - 2D or 3D.

my_rank = ndims(field);
is2d = false;
is3d = false;

switch my_rank
  case 2
    if isvector(field)
      error('Input FIELD is a 1D vector - that makes no sense');
    end
    is2d = true;
    [Im,Jm] = size(field);
  case 3
    is3d = true;
    [Im,Jm,Km] = size(field);
  otherwise
    disp('Input FIELD must be a 2D or 3D array at 1 time slice.')
    error([' Number of dimensions = ',num2str(my_rank)])
end

if ~is3d
  F = rmfield(F, 'z_w');             % remove structure fields associated
  F = rmfield(F, 'z_r');             % with a 3D input variable.
  F = rmfield(F, 'dz');
end

F.field = field;
F.size = size(field);
F.method = method;

% Extract the ROMS grid lon/lat or x/y coordinates for this field.

if (Im == L && Jm == M)
  
  F.Cgrid = 'psi_point';
  if G.spherical
    F.Xgrd = G.lon_psi;
    F.Ygrd = G.lat_psi;
  else
    F.Xgrd = G.x_psi;
    F.Ygrd = G.y_psi;
  end

  if isfield(G,'mask_psi')
    mask = G.mask_psi;
  else
    mask = ones(size(F.Xgrd));
  end

  h  = 0.25 .* (G.h(1:L,1:M ) + G.h(2:Lp,1:M ) +                        ...
                G.h(1:L,2:Mp) + G.h(2:Lp,2:Mp));

elseif (Im == Lp && Jm == Mp)
  
  F.Cgrid = 'rho_point';
  if (G.spherical)
    F.Xgrd = G.lon_rho;
    F.Ygrd = G.lat_rho;
  else
    F.Xgrd = G.x_rho;
    F.Ygrd = G.y_rho;
  end

  if isfield(G,'mask_rho')
    mask = G.mask_rho;
  else
    mask = ones(size(F.Xgrd));
  end

  h  = G.h;

elseif (Im == L && Jm == Mp)
  
  F.Cgrid = 'u_point';
  if (G.spherical)
    F.Xgrd = G.lon_u;
    F.Ygrd = G.lat_u;
  else
    F.Xgrd = G.x_u;
    F.Ygrd = G.y_u;
  end

  if isfield(G,'mask_u')
    mask = G.mask_u;
  else
    mask = ones(size(F.Xgrd));
  end

  h  = 0.5 .* (G.h(1:L,1:Mp) + G.h(2:Lp,1:Mp));

elseif (Im == Lp && Jm == M)

  F.Cgrid = 'v_point';
  if (G.spherical)
    F.Xgrd = G.lon_v;
    F.Ygrd = G.lat_v;
  else
    F.Xgrd = G.x_v;
    F.Ygrd = G.y_v;
  end

  if (isfield(G,'mask_v')),
    mask = G.mask_v;
  else
    mask = ones(size(F.Xgrd));
  end

  h  = 0.5 .* (G.h(1:Lp,1:M) + G.h(1:Lp,2:Mp));

else

  error([' INTERPOLATOR: unable to determine Arakawa C-grid type:',     ...
         ' Cgrid is inconsistent with ROMS grid.']);
end

%--------------------------------------------------------------------------
% Get bathymetry and z-coords along requested path.
%--------------------------------------------------------------------------

% Use scatteredInterpolant - ROMS horizontal coordinates are not plaid.

Fi = scatteredInterpolant(F.Xgrd(:), F.Ygrd(:), h(:), method);
F.h = Fi(Xpath(:), Ypath(:));

if is3d
  
  % Compute depths (for plotting) and layer thinknesses (for vertical
  % integration of metric/functionals)
  
  zeta = zeros(size(F.h));
  igrid = 1;
  z_r = set_depth(G.Vtransform, G.Vstretching, G.theta_s,               ...
                  G.theta_b, G.hc, G.N, igrid, F.h, zeta, false);
  F.z_r = squeeze(z_r);
  
  igrid = 5;
  z_w = set_depth(G.Vtransform, G.Vstretching, G.theta_s,               ...
                  G.theta_b, G.hc, G.N, igrid, F.h, zeta, false);
  F.z_w = squeeze(z_w);
  F.dz = diff(F.z_w,1,2);
  
  if ~isvector(Xpath)
    F.z_w = reshape(F.z_w, [size(Xpath) Km+1]);
    F.z_r = reshape(F.z_r, [size(Xpath) Km]);
    F.dz  = reshape(F.dz, [size(Xpath) Km]);
  end

end

if ~isvector(Xpath)
  F.h  = reshape(F.h, size(Xpath));
end

%--------------------------------------------------------------------------
% Path distance and segment lengths (for along-path integration).
%--------------------------------------------------------------------------

if isvector(Xpath)
  
  if G.spherical
    dis = cumsum([0; sw_dist(Ypath(:),Xpath(:),'km')]);
    
    if is3d
      F.dis = repmat(dis,[1 Km]);
    else
      F.dis = dis;
    end
    
    % Path segment lengths for along path integral.
    % The cumulative distance from s=0 is at the path points. Average to
    % between-path points so that we can use diff to calculate path
    % segment lengths centered on the path points. But first Pad with the
    % first and last values of the original dis. Then diff will give half
    % cell widths for the end points of the path

    dtmp = 0.5*(dis(1:(end-1))+dis(2:end));
    dtmp = [dis(1); dtmp; dis(end)];
    dislen = diff(dtmp);
    
    if is3d
      F.dislen = repmat(dislen,[1 Km]);
    else
      F.dislen = dislen;
    end
    
  else
    error('still need code for path in grid or normalized coordinates')
  end
  
else
   warning('No code yet to compute areas for input coordinates no a grid')
end

%--------------------------------------------------------------------------
% Set fractional coordinates and linear indices containing the points
% to interpolate. The fractional coordinates origin is irrelevant here
% since is only used locally.
%--------------------------------------------------------------------------

[Yfrac, Xfrac] = meshgrid(1:Jm, 1:Im);

Fi = scatteredInterpolant(F.Xgrd(:), F.Ygrd(:), Xfrac(:), method);
Ifrac = Fi(Xpath(:), Ypath(:));
Fi.Values = Yfrac(:);
Jfrac = Fi(Xpath(:), Ypath(:));

I   = fix(Ifrac);                                     % I
Ip1 = min(Im, I+1);                                   % I+1
J   = fix(Jfrac);                                     % J
Jp1 = min(Jm, J+1);                                   % J+1

% Set grid cell linear indices used in the interpolation.

index1 = sub2ind([Im,Jm], I  , J  );                  % (I  , J  )
index2 = sub2ind([Im,Jm], Ip1, J  );                  % (I+1, J  )
index3 = sub2ind([Im,Jm], Ip1, Jp1);                  % (I+1, J+1)
index4 = sub2ind([Im,Jm], I  , Jp1);                  % (I  , J+1)

F.Lindex1 = index1;
F.Lindex2 = index2;
F.Lindex3 = index3;
F.Lindex4 = index4;

p = Ifrac - Xfrac(index1);    % recall that we have unity spacing
q = Jfrac - Yfrac(index1);    % in the fractional coordinates

W1 = mask(index1) .* (1 - p) .* (1 - q);
W2 = mask(index2) .* p .* (1 - q);
W3 = mask(index3) .* p .* q;
W4 = mask(index4) .* (1 -p) .* q;

MaskSum = mask(index1) +                                                ...
          mask(index2) +                                                ...
          mask(index3) +                                                ...
          mask(index4);

if any(MaskSum < 4)                        % at least one point is land
  for n = 1:numel(Xpath)
    if MaskSum(n) < 4
      Wsum = W1(n) + W2(n) + W3(n) + W4(n);
      if (Wsum > 0)
        W1(n) = W1(n) ./ Wsum;             % using only water points
        W2(n) = W2(n) ./ Wsum;
        W3(n) = W3(n) ./ Wsum;
        W4(n) = W4(n) ./ Wsum;
      else
        W1(n) = 0;                         % all points are on liness
        W2(n) = 0;
        W3(n) = 0;
        W4(n) = 0;
      end
    end
  end
end

%--------------------------------------------------------------------------
% Interpolate. This is the only part that is adjointable.
%--------------------------------------------------------------------------

F.Xval = Xpath;
F.Yval = Ypath;

if is2d
  
  % 2D field

  F.val = W1 .* F.field(index1) +                                       ...
          W2 .* F.field(index2) +                                       ...
          W3 .* F.field(index3) +                                       ...
          W4 .* F.field(index4);

  if ~isvector(Xpath)
    F.val = reshape(F.val, size(Xpath));
  end

else

  % 3D field - do all k levels

  for k=1:Km
    vartmp = squeeze(F.field(:,:,k));
    F.val(:,k) = W1 .* vartmp(index1) +                                 ...
                 W2 .* vartmp(index2) +                                 ...
                 W3 .* vartmp(index3) +                                 ...
                 W4 .* vartmp(index4);
  end

  if ~isvector(Xpath)
    F.val = reshape(F.val, [size(Xpath) Km]);
  end
  
end

