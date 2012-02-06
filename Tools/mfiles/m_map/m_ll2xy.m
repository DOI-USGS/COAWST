function [X,Y]=m_ll2xy(varargin);
% M_LL2XY Converts long,lat to X,Y coordinates using the current projection
%         [X,Y]=m_ll2xy(LONGITUDE,LATITUDE);
%
%         Extra properties can be added after the latitude variable:
%          ...,'clip', ( 'on' | 'off' | 'patch'  | 'point' )
%         where normal clipping sets out-of-range values to NaN, and patch
%         clipping sets out-of-range values to border values (useful for
%         drawing patches). A border point is interpolated for line
%         segments crossing the border; line segments are assumed for
%         one-dimensional input arrays and for the columns of 2-dimensional
%         arrays. The interpolation can be disabled using the 'point'
%         option (i.e. data is composed of discrete unrelated points).
%

% Rich Pawlowicz (rich@ocgy.ubc.ca) 2/Apr/1997
%
% This software is provided "as is" without warranty of any kind. But
% it's mine, so you can't sell it.

% 6/Nov/00 - eliminate returned stuff if ';' neglected (thx to D Byrne)

global MAP_PROJECTION 

if nargin==0 | isstr(varargin{1}),
  disp(' Usage');
  disp(' [X,Y]=m_ll2xy(LONGITUDES,LATITUDES <,''clip'',( ''on''|''off''|''patch'' | ''point'' ) >)');
else
   % Sneaky way of making default clipping on (sneaky 'cause only the 4th
   % input parameter is checked for the clipping property)
  [X,Y]=feval(MAP_PROJECTION.routine,'ll2xy',varargin{:},'clip','on');
end;

if nargout==0,
 clear X Y
end;

