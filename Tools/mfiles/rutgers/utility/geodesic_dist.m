function [gdist,galpha]=geodesic_dist(lon1,lat1,lon2,lat2,flag);

%
% GEODESIC_DIST:  Computes geodesic distance between two points
%
% [gdist,galpha]=geodesic_dist(lon1,lat1,lon2,lat2,flag);
%
% Inverse, non-iterative solutions for distance and geodesic azimuth
% between two points on the ellipsoid (The Earth) from the equations,
% second order in spheroidal flatttening, given by:
%
% Sodano , E.M., and T. Robinson, 1963: Direct and inverse solutions
%   of geodesics, Army Map Service Technical Report No. 7, AD 657591.
%
% On Input:    Longitude is positive to the east and negative to the
%              west.  Latitude is positive to the north and negative
%              to the south
%
%    lon1      Longitude point 1 (decimal degrees)
%    lat1      Latitude  point 1 (decimal degrees)
%    lon2      Longitude point 2 (decimal degrees)
%    lat2      Latitude  point 2 (decimal degrees)
%    flag      flag for distance units on output:
%                flag=1  => meters
%                flag=2  => nautical miles
%                flag=3  => feet
%                flag=4  => kilometers
%                flag=5  => statute miles
%
% On Output:
%
%    gdist     Geodesic distance between point 1 and point 2.
%    galpha    Geodesic azimuth from point 1 to point 2 clockwise from
%                North (decimal degrees)
%

% svn $Id: geodesic_dist.m 711 2014-01-23 20:36:13Z arango $
%===========================================================================%
%  Copyright (c) 2002-2014 The ROMS/TOMS Group                              %
%    Licensed under a MIT/X style license                                   %
%    See License_ROMS.txt                           Hernan G. Arango        %
%===========================================================================%

%  Determine output distance units.

if (nargin < 5),
  flag=1;
end,

switch (flag),
  case 1
    dscale=1.0;                 % meters
  case 2  
    dscale=5.396d-4;            % nautical miles
  case 3
    dscale=3.281;               % feet
  case 4
    dscale=0.001;               % kilometers
  case 5
    dscale=6.214d-4;            % statute mile
end,

deg2rad=pi/180;
rad2deg=180/pi;

%  Define parameters on first pass (SMIN: Ellipsoid semi-minor axis in
%  meters; SMAJ: Ellipsoid semi-major axis in meters; F: spheroidal
%  flattening).

smin=6356750.52;
smaj=6378135.0;
f=1-(smin/smaj);

%  Return if zero distance.

if ((lon1 == lon2) & (lat1 == lat2)),
  gdist=0;
  galpha=0;
  return
end,

%  Determine proper longitudinal shift.

delta=lon2-lon1;
l=abs(delta);
if (l >= 180),
  l=360-abs(lon1-lon2);
end,

%  Convert Decimal degrees to radians.

r_lat1=lat1*deg2rad;
r_lat2=lat2*deg2rad;
l=l*deg2rad;

%  Calculate S/Bo subformulas.

beta1=atan(tan(r_lat1)*(1-f));
beta2=atan(tan(r_lat2)*(1-f));
a=sin(beta1)*sin(beta2);
b=cos(beta1)*cos(beta2);
ct=a+b*cos(l);
st=sqrt(((sin(l)*cos(beta2))^2)+ ...
       (((sin(beta2)*cos(beta1))-(sin(beta1)*cos(beta2)*cos(l)))^2));
t=asin(st);
c=(b*sin(l))/st;
m=1-(c*c);

%  Calculate S/Bo term.

q=f+(f*f);
z=0.5*f*f;
x=0.0625*f*f;
y=0.125*f*f;
w=0.25*f*f;

sob=((1+q)*t)+(a*((q*st)-(z*(t*t)*(1/sin(t))))) + ...
    (m*(((-0.5*q)*t)-((0.5*q)*st*ct)+(z*(t*t)*(1/tan(t))))) + ...
    ((a*a)*(-z*st*ct)) + ...
    ((m*m)*((x*t)+(x*st*ct)-(z*(t*t)*(1/tan(t))) - ...
    (y*st*(ct*ct*ct)))) + ...
    ((a*m)*((z*(t*t)*(1/sin(t)))+(z*st*(ct*ct))));

gdist=dscale*sob*smin;

%  Compute geodesic azimuth from point 1 to point 2 clockwise from
%  North, alpha.

lambda=q*t+a*(-z*st-f*f*t*t/sin(t)) + ...
       m*(-5*w*t+w*st*cos(t)+f*f*t*t/tan(t));
lambda=c*lambda+l;

if (lambda == 0),
  if (lat1 < lat2), galpha=0; end,
  if (lat1 > lat2), galpha=180; end,
  return
end,

cott=(sin(beta2)*cos(beta1)-cos(lambda)*sin(beta1)*cos(beta2))/ ...
     (sin(lambda)*cos(beta2));
if (cott == 0),
  alpha=90;
else
  alpha=atan(1/cott)*rad2deg;
end,

%  Compute heading from point 1 to point 2 clockwise from north.
 
if (delta > 0),
  if (cott > 0),
    galpha=alpha;
  elseif (cott < 0),
    galpha=180+alpha;
  end,
end,

if (delta < 0),
  if (cott < 0),
    galpha=180-alpha;
  elseif (cott > 0),
    galpha=360-alpha;
  end,
end,

return
