function [dist,bearing]=gcircle(lon1,lat1,lon2,lat2,BEARING);

%
% GCIRCLE:  Coumputes great circle distance between to points
%
% [dist,bearing]=gcircle(lon1,lat1,lon2,lat2,BEARING);
%
%  This function computes great circle distance and initial bearing
%  between two longitude and latitude points. The Earth is assumed
%  to be a sphere. This approximation is valid for short distances.
%
%  On Input:  Longitude is positive to the east and negative to the
%             west.  Latitude is positive to the north and negative
%             to the south.
%
%     lon1     longitude point 1 (decimal degrees)
%     lat1     latitude  point 1 (decimal degrees)
%     lon2     longitude point 2 (decimal degrees)
%     lat2     latitude  point 2 (decimal degrees)
%     BEARING  Logical switch to compute bearing angle (0/1)
%
%  On Output:
%
%     bearing  Azimuth from point 1 to point 2 counterclockwise
%                from true EAST (decimal degrees)
%     dist     great circle distance between point 1 and point 2
%                (kilometers)
%
% Adapted from routine written by Pat J. Haley (Harvard University).
%

% svn $Id: gcircle.m 938 2019-01-28 06:35:10Z arango $
%===========================================================================%
%  Copyright (c) 2002-2019 The ROMS/TOMS Group                              %
%    Licensed under a MIT/X style license                                   %
%    See License_ROMS.txt                           Hernan G. Arango        %
%===========================================================================%

if (nargin < 5),
  BEARING=1;
end,

%----------------------------------------------------------------------------
% Set often used parameters.
%----------------------------------------------------------------------------

radius  = 6371.315;
deg2rad = pi/180.0;
rad2deg = 180.0/pi;
tng_cut = 2*eps;
small   = sqrt(realmin);

%----------------------------------------------------------------------------
% Convert to radians.
%----------------------------------------------------------------------------

slon = lon1.*deg2rad;
slat = lat1.*deg2rad;
elon = lon2.*deg2rad;
elat = lat2.*deg2rad;

%----------------------------------------------------------------------------
% Compute distance along great circle (kilometers).
%----------------------------------------------------------------------------

alpha = sin(slat).*sin(elat) + cos(slat).*cos(elat).*cos(elon-slon);

ind = find (abs(alpha)>1);
if (~isempty(ind)), alpha(ind) = sign(alpha(ind)); end;

alpha=acos(alpha);
dist=radius.*alpha;

%----------------------------------------------------------------------------
%  Compute bearing angle: angle from the first point to the second 
%  point clockwise from true NORTH.
%----------------------------------------------------------------------------

if (BEARING),

  phi = cos(elat).*sin(elon-slon);
  ang = zeros(size(alpha));

  ind1 = find (sin(alpha) ~= 0);
  if (~isempty(ind1))
    ang(ind1) = phi(ind1)./sin(alpha(ind1));
    ind2 = find (abs(ang(ind1)) > 1);
    if (~isempty(ind2)), ang(ind1(ind2)) = sign(ang(ind1(ind2))); end;
    ang(ind1) = asin(ang(ind1)).*rad2deg;
  end;

  ind1 = find (sin(alpha) == 0);
  if (~isempty(ind))
    ind2 = find (phi(ind1) == 0);
    if (~isempty(ind2)), ang(ind1(ind2)) = zeros(size(ind2)); end;

    ind2 = find (phi(ind1) > 0);
    if (~isempty(ind2)), ang(ind1(ind2)) = 90.*ones(size(ind2)); end;

    ind2 = find (phi(ind1) < 0);
    if (~isempty(ind2)), ang(ind1(ind2)) = -90.*ones(size(ind2)); end;
  end;

  ind1 = find (elat < slat);
% if (~isempty(ind1)), ang(ind1) = 180 - ang(ind1); end;

  bearing = rem (ang+360, 360);

%  Compute bearing from point 1 to point 2 counterclockwise from
%  true EAST.

  ind1 = find ( (0 <= bearing) & (bearing <= 90) );
  if (~isempty(ind1)),
    bearing(ind1) = 90 - bearing(ind1);
    ind2 = find ( alpha(ind1) == 0);
    if (~isempty(ind2)), bearing(ind1((ind2))) = zeros(size(ind2)); end;
  end;

  ind1 = find ( (0 > bearing) | (bearing > 90) );
  if (~isempty(ind1)),
    bearing(ind1) = rem ( 450-bearing(ind1), 360 );
  end;

else,

  bearing=ones(size(dist)).*NaN;

end,
  
return
