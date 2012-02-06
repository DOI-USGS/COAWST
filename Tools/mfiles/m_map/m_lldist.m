function dist = m_lldist(long,lat)
% M_LLDIST Spherical earth distance between points in long/lat coordinates. 
%   RANGE=M_LLDIST(LONG,LAT) gives the distance in meters between
%   successive points in the vectors LONG and LAT, computed
%   using the Haversine formula on a spherical earth of radius
%   6378.137km. Distances are probably good to better than 1% of the
%   "true" distance on the ellipsoidal earth
%
%   See also M_XYDIST

% Rich Pawlowicz (rich@ocgy.ubc.ca) 6/Nov/00
% This software is provided "as is" without warranty of any kind. But
% it's mine, so you can't sell it.

pi180=pi/180;
earth_radius=6378.137e3;


long1=long(1:end-1)*pi180;
long2=long(2:end)*pi180;
lat1=lat(1:end-1)*pi180;
lat2=lat(2:end)*pi180;

dlon = long2 - long1; 
dlat = lat2 - lat1; 
a = (sin(dlat/2)).^2 + cos(lat1) .* cos(lat2) .* (sin(dlon/2)).^2;
c = 2 * atan2( sqrt(a), sqrt(1-a) );
dist = earth_radius * c;


