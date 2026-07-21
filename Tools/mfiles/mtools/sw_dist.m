
function [dist,phaseangle] = distance(lat,lon,units)

% SW_DIST    Distance between two lat,lon coordinates
%===================================================================
% SW_DIST  $Revision: 1.4 $  $Date: 1994/10/10 04:55:23 $
%          Copyright (C) CSIRO, Phil Morgan & Steve Rintoul 1992. 
%
% USAGE:  [dist,phaseangle] = distance(lat,lon {,units} )
%
% DESCRIPTION:
%   Calculate distance between two positions on glode using the "Plane
%   Sailing" method.  Also uses simple geometry to calculate the bearing of
%   the path between position pairs.
% 
% INPUT:
%    lat      = decimal degrees (+ve N, -ve S) [- 90.. +90]
%    lon      = decimal degrees (+ve E, -ve W) [-180..+180]
%    units    = optional string specifing units of distance
%               'nm'  = nautical miles (default)
%               'km'  = kilometres
%
% OUTPUT:
%    dist        = distance between positions in units
%    phaseangle  = angle of line between stations with x axis (East).
%                  Range of values are -180..+180. (E=0, N=90, S=-90)
%
% AUTHOR:   Phil Morgan and Steve Rintoul 92-02-10
%
% DISCLAIMER:
%   This software is provided "as is" without warranty of any kind.  
%   See the file sw_copy.m for conditions of use and licence.
% 
% REFERENCE:
%    The PLANE SAILING method as descriibed in "CELESTIAL NAVIGATION" 1989 by
%    Dr. P. Gormley. The Australian Antartic Division.
%==================================================================

% CALLER:   general purpose
% CALLEE:   none

%----------------------
% CHECK INPUT ARGUMENTS
%----------------------
if nargin > 3
  error('sw_dist.m: No more than 3 arguments allowed')
elseif nargin==3 
  if ~isstr(units)
      error('sw_dist.m: units argument must be string')
  end %if
elseif nargin==2
  units = 'nm';  % default units
else
  error('sw_dist.m: wrong number of arguments')
end%if

[mlat,nlat] = size(lat);
if mlat~=1 & nlat~=1
   error('sw_dist.m: lat, lon must be vectors.  No matrices allowed')
else
  if mlat == 1
    Transpose = 1;  % row vector passed in
  else
    Transpose = 0;  % accept column vector
  end%if
end%if
lat = lat(:); %force to column vectors
lon = lon(:);
if length(lat)~=length(lon)
   error('sw_dist.m: lat and lon must have same number of elements')
end%if

%-----------------
% DEFINE CONSTANTS
%-----------------
DEG2RAD = (2*pi/360);
RAD2DEG = 1/DEG2RAD;
DEG2MIN = 60;
DEG2NM  = 60;
NM2KM   = 1.8520;    % Defined in Pond & Pickard p303.

% BEGIN
npositions = length(lat);
ind=1:npositions-1;     % index to first of position pairs

dlon = diff(lon);
if any(abs(dlon)>180)
   flag = find(abs(dlon)>180);
   for ii=1:length(flag)
     dlon(flag(ii))= -sign(dlon(flag(ii))) * (360 - abs(dlon(flag(ii))) );
   end %for
end %if
latrad = abs(lat*DEG2RAD);
dep    = cos( (latrad(ind+1)+latrad(ind))./2 ) .* dlon;
dlat   = diff(lat);
dist   = DEG2NM*sqrt(dlat.^2 + dep.^2);  % in n.miles

if strcmp(units,'km')   % defaults to n.miles
    dist = dist * NM2KM;
end %if

% CALCUALTE ANGLE TO X AXIS
phaseangle  = angle(dep+dlat*sqrt(-1))*RAD2DEG;
 
if Transpose
  dist = dist';
  phaseangle = phaseangle';
end %if

return
%--------------------------------------------------------------------

