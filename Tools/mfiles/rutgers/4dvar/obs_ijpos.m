function [Xgrid, Ygrid]=obs_ijpos(GRDname, obs_lon, obs_lat, ...
                                  Correction, obc_edge, ...
                                  Ioffset, Joffset);

%
% OBS_IJPOS:  Computes observation locations in ROMS fractional coordinates
%
% [Xgrid,Ygrid]=obs_ijpos(GRDname,obs_lon,obs_lat,Correction,obc_edge,
%                         Ioffset,Joffset)
%
% This function computes the observation locations (Xgrid,Ygrid) in terms
% of ROMS fractional (I,J) coordinates. This is done to facilitate the
% processing of the observation operators inside ROMS. All the observations
% are assumed to be located at RHO-points.  If the Ioffset and Joffset
% vectors are provided, the polygon defined by the application grid is
% smaller by the number of grid points provided in the offset.  This is
% done to avoid processing next to the applications boundary edges.
%
% On Input:
%
%    GRDname       NetCDF grid file name (string)
%    obs_lon       observation longitude (positive, degrees_east)
%    obs_lat       observation latitude  (positive, degrees_north)
%    Correction    switch to apply correction due to spherical/curvilinear
%                    grids (false, true)
%    obc_edge      switch to include observations on open boundary edges
%                    (false, true)
%    Ioffset       Application I-grid offset when defining polygon (vector):
%                    Ioffset(1):  I-grid offset on the edge where Istr=1
%                    Ioffset(2):  I-grid offset on the edge where Iend=Lm
%    Joffset       Application J-grid offset when defining polygon (vector):
%                    Joffset(1):  I-grid offset on the edge where Jstr=1
%                    Joffset(2):  I-grid offset on the edge where Jend=Mm
%
% On Ouput:
%
%    Xgrid         observation fractional x-grid location
%    Ygrid         observation fractional y-grid location
%
% Notice:  Outlier observations outside of ROMS grid has an NaN value in
%          Xgrid and Ygrid.
%
  
% svn $Id: obs_ijpos.m 996 2020-01-10 04:28:56Z arango $
%===========================================================================%
%  Copyright (c) 2002-2020 The ROMS/TOMS Group                              %
%    Licensed under a MIT/X style license                                   %
%    See License_ROMS.txt                           Hernan G. Arango        %
%===========================================================================%

IPLOT=0;          % switch for plotting during debugging

% If appropriate, transpose input vectors so they are of size: [length 1].

lsize=size(obs_lon);
if (lsize(1) < lsize(2)),
  obs_lon = obs_lon';
end,

lsize=size(obs_lat);
if (lsize(1) < lsize(2)),
  obs_lat = obs_lat';
end,

if (nargin < 4),
  Correction = false;
end,

if (nargin < 5),
  obc_edge = false;
end,

if (nargin < 6),
  Ioffset(1)=0;
  Ioffset(2)=0;
end,

if (nargin < 7),
  Joffset(1)=0;
  Joffset(2)=0;
end,

%----------------------------------------------------------------------------
% Read in model grid variables.
%----------------------------------------------------------------------------

V = nc_vnames(GRDname);
nvars = length(V.Variables);
 
got.spherical = false;
got.angle     = false;
got.mask_rho  = false;
got.lon_rho   = false;
got.lat_rho   = false;
got.coast     = false;

for n=1:nvars
  name = char(V.Variables(n).Name);
  switch name
    case 'spherical'
      got.spherical = true;
    case 'lon_rho'
      got.lon_rho   = true;
    case 'lat_rho'
      got.lat_rho   = true;
    case 'angle'
      got.angle     = true;
    case 'mask_rho'
      got.mask_rho  = true;
    case 'lon_coast'
      got.coast     = true;
      clon = nc_read(GRDname, 'lon_coast');
      clat = nc_read(GRDname, 'lat_coast');
  end,
end,

if (got.spherical),
  spherical=nc_read(GRDname,'spherical');
  if (ischar(spherical)),
    if (spherical == 'T' || spherical == 't'),
      spherical = true;
    else,
      spherical = false;
    end,
  end,
else,
  spherical = true;
end,

if (got.lon_rho),
  rlon = nc_read(GRDname, 'lon_rho');
else,
  error('OBS_IJPOS - cannot find variable: lon_rho');
end,

if (got.lat_rho),
  rlat = nc_read(GRDname, 'lat_rho');
else,
  error('OBS_IJPOS - cannot find variable: lat_rho');
end,

if (got.angle),
  angle = nc_read(GRDname,'angle');
else
  angle = zeros(size(rlon));
end,

if (got.mask_rho),
  rmask = nc_read(GRDname,'mask_rho');
else
  rmask = ones(size(rlon));
end,

%----------------------------------------------------------------------------
%  Extract polygon defining application grid box.
%----------------------------------------------------------------------------

%  Set grid application polygon.

[Im,Jm]=size(rlon);

Istr=1 +Ioffset(1);
Iend=Im-Ioffset(2);
Jstr=1 +Joffset(1);
Jend=Jm-Joffset(2);

Xbox=[squeeze(rlon(Istr:Iend,Jstr)); ...
      squeeze(rlon(Iend,Jstr+1:Jend))'; ...
      squeeze(flipud(rlon(Istr:Iend-1,Jend))); ...
      squeeze(fliplr(rlon(Istr,Jstr:Jend-1)))'];

Ybox=[squeeze(rlat(Istr:Iend,Jstr)); ...
      squeeze(rlat(Iend,Jstr+1:Jend))'; ...
      squeeze(flipud(rlat(Istr:Iend-1,Jend))); ...
      squeeze(fliplr(rlat(Istr,Jstr:Jend-1)))'];

%  Find observation inside (IN) or on the edge (ON) the polygon defined
%  by (Xbox,Ybox).

[IN ON]=inpolygon(obs_lon, obs_lat, Xbox, Ybox);

%  Flag outlier observations as bounded=false.  We are only considering
%  observations inside the polygon.

bounded=false(size(obs_lon));

bounded(IN) = true;
if (obc_edge),                 % process observations on boundary edges
  bounded(ON) = true;
end,

% Plot observations in the application grid.

if (IPLOT),
  pcolor(rlon,rlat,ones(size(rlon)));
  hold on;
  set(gca, 'fontsize', 14, 'fontweight', 'bold');
  if (got.coast),
    plot(clon, clat, 'k-');
  end,
  plot(Xbox, Ybox, 'b-', ...
       obs_lon(IN), obs_lat(IN), 'b.', ...
       obs_lon(ON), obs_lat(ON), 'r.', ...
       obs_lon(~IN), obs_lat(~IN), 'm.');
  title(['Blue points (inside),  ', ...
         'Red points (boundary),  ', ...
         'Magenta Points (outliers)'], ...
         'fontsize', 14, 'fontweight', 'bold');
  hold off;
end,

clear IN ON Xbox Ybox

%----------------------------------------------------------------------------
% Compute model grid fractional (I,J) locations at observation locations
% via interpolation.
%----------------------------------------------------------------------------

Igrid=repmat([0:1:Im-1]', [1 Jm]);
Jgrid=repmat([0:1:Jm-1] , [Im 1]);

if (got.mask_rho),
  ind = find(rmask < 1);
  if (~isempty(ind)),
    Igrid(ind) = NaN;
    Jgrid(ind) = NaN;
  end,
end,

Xgrid = ones(size(obs_lon)) .* NaN;    % initialize unbounded observations
Ygrid = ones(size(obs_lon)) .* NaN;    % to NaN

Xgrid(bounded) = griddata(rlon, rlat, Igrid, obs_lon(bounded), obs_lat(bounded));
Ygrid(bounded) = griddata(rlon, rlat, Jgrid, obs_lon(bounded), obs_lat(bounded));

clear Igrid Jgrid

%  If land/sea masking, find the observation in land (Xgrid=Ygrid=NaN);

ind = find(isnan(Xgrid) | isnan(Ygrid));
if (~isempty(ind)),
  bounded(ind) = false;
end,

clear ind

%----------------------------------------------------------------------------
%  Spherical/Curvilinear corrections.
%----------------------------------------------------------------------------

if (Correction),
  [Xgrid, Ygrid] = correction(rlon, rlat, angle, obs_lon, obs_lat, ...
                              bounded, Xgrid, Ygrid);
end,

return


function [X,Y] = correction(rlon, rlat, angle, obs_lon, obs_lat, bounded, X, Y);

%
% CORRECTION:  Apply curvilinear coordinates correction to observations
%              fractional (X,Y) locations.
%
% The maximum variable size allowed by Matlab can be exceeded very quickly
% if the observation vector is large. This function is coded either using  
% block temporary arrays or a simple loop where the corrections are computed
% one by one. This is one of the few instances in Matlab that actually is
% more efficient to do complex scalar operations than vector operations.
% Notice that the number of satellite observation can be large and often
% exceeds 1E5.
% 
% On Input:
%
%    rlon          Application grid longitude at RHO-points (matrix)
%    rlat          Application grid latitude  at RHO-points (matrix)
%    angle         curvilinear grid rotation (radians; matrix)
%    obs_lon       Observation longitude locations (vector)
%    obs_lat       Observation latitude  locations (vector)
%    bounded       Switch marking outlier (false) points (vector)
%    X             Observation fractional x-grid location (first guess)
%    Y             Observation fractional y-grid location (first guess)
%
% On Output:
%
%    X             Adjusted observation fractional x-grid location
%    Y             Adjusted observation fractional y-grid location
%

%===========================================================================%
%  Copyright (c) 2002-2020 The ROMS/TOMS Group                              %
%    Licensed under a MIT/X style license                                   %
%    See License_ROMS.txt                           Hernan G. Arango        %
%===========================================================================%

debugging = false;                  % debugging switch

Nobs = length(obs_lon);             % number of observations

block_length = 100000;              % size of block arrays

Eradius = 6371315.0;                % Earth radius (meters)
deg2rad = pi/180;                   % degrees to radians factor

[Lr, Mr] = size(rlon);              % Number of rho-points

Lm = Lr - 1;
Mm = Mr - 1;

%----------------------------------------------------------------------------
%  Set (I,J) coordinates of the grid cell containing the observation
%  need to add 1 because zero lower bound in ROMS "rlon"
%----------------------------------------------------------------------------

I = fix(X);
ind = find((0 <= I) & (I < Lm));
if (~isempty(ind)),
  I(ind) = I(ind) + 1;
end,

J = fix(Y);
ind = find((0 <= J) & (J < Mm));
if (~isempty(ind)),
  J(ind) = J(ind) + 1;
end,

if (debugging),
  disp(' ');
  disp(['  Xmin = ', num2str(min(X),'%6.2f'), ...
        '  Xmax = ', num2str(max(X),'%6.2f'), ...
        '  Imin = ', num2str(min(I),'%3.3i'), ...
        '  Imax = ', num2str(max(I),'%3.3i'), ...
        '  Lr = ',   num2str(Lr,'%3.3i')]);

  disp(['  Ymin = ', num2str(min(Y),'%6.2f'), ...
        '  Ymax = ', num2str(max(Y),'%6.2f'), ...
        '  Jmin = ', num2str(min(J),'%3.3i'), ...
        '  Jmax = ', num2str(max(J),'%3.3i'), ...
        '  Mr = ',   num2str(Mr,'%3.3i')]);
  disp(' ');
end,

clear ind

%----------------------------------------------------------------------------
%  It is possible that we are processing a large number of observations.
%  Therefore, the observation vector is processed by blocks to reduce
%  the memory requirements.
%----------------------------------------------------------------------------

N  = ceil(Nobs / block_length);
n1 = 0;
n2 = 0;

while (n2 < Nobs),

  n1 = n2 + 1;
  n2 = n1 + block_length;
  if (n2 > Nobs)
    n2 = Nobs;
  end,

  iobs = n1:n2;
  ind  = find(~bounded(n1:n2));
  if (~isempty(ind)),             % remove unbouded observations, if any
    iobs(ind) = [];
  end,

  i_j   = sub2ind(size(rlon), I(iobs)  , J(iobs)  );
  ip1_j = sub2ind(size(rlon), I(iobs)+1, J(iobs)  );
  i_jp1 = sub2ind(size(rlon), I(iobs)  , J(iobs)+1);

  if (debugging),
    disp(['  Processing observation vector, n1:n2 = ' ...
          num2str(n1,'%7.7i'), ' - ', num2str(n2,'%7.7i') ...
          '  size = ', num2str(length(iobs))]);
  end,

%  Convert all positions to meters first.  

  yfac = Eradius * deg2rad;
  xfac = yfac .* cos(obs_lat(iobs) .* deg2rad);

  xpp  = (obs_lon(iobs) - rlon(i_j)) .* xfac;
  ypp  = (obs_lat(iobs) - rlat(i_j)) .* yfac;

%  Use Law of Cosines to get cell parallelogram "shear" angle.

  diag2 = (rlon(ip1_j) - rlon(i_jp1)) .^ 2 + ...
          (rlat(ip1_j) - rlat(i_jp1)) .^ 2;

  aa2   = (rlon(i_j)   - rlon(ip1_j)) .^ 2 + ...
          (rlat(i_j)   - rlat(ip1_j)) .^ 2;

  bb2   = (rlon(i_j)   - rlon(i_jp1)) .^ 2 + ...
          (rlat(i_j)   - rlat(i_jp1)) .^ 2;

  phi   = asin((diag2 - aa2- bb2) ./ (2 .* sqrt(aa2 .* bb2)));

%  Transform fractional locations into curvilinear coordinates. Assume
%  the cell is rectanglar, for now.

  dx = xpp .* cos(angle(i_j)) + ypp .* sin(angle(i_j));

  dy = ypp .* cos(angle(i_j)) - xpp .* sin(angle(i_j));

%  Correct for parallelogram.

  dx = dx + dy .* tan(phi);
  dy = dy ./ cos(phi);

%  Scale with cell side lengths to translate into cell indexes.

  dx = dx ./ sqrt(aa2) ./ xfac;

  ind = find(dx < 0);
  if (~isempty(ind)),
    dx(ind) = 0;
  end,

  ind = find(dx > 1);
  if (~isempty(ind)),
    dx(ind) = 1;
  end,

  dy = dy ./ sqrt(bb2) ./ yfac;

  ind = find(dy < 0);
  if (~isempty(ind)),
    dy(ind) = 0;
  end,

  ind = find(dy > 1);
  if (~isempty(ind)),
    dy(ind) = 1;
  end,

  X(iobs) = fix(X(iobs)) + dx;
  Y(iobs) = fix(Y(iobs)) + dy;

  clear aa2 bb2 diag2 dx dy ind iobs i_j ip1_j i_jp1 phi xfac xpp ypp

end,

return
