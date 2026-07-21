function [ALat, ALon, spN, spE, Y, X] = frfCoord(p1, p2)
% function [ALat, ALon, spN, spE, Y, X] = frfCoord(p1, p2)
%
%  15 Dec 2014
%  Kent Hathaway.
%  Uses new fit (angles and scales) Bill Birkemeier determined in Nov 2014
%  
%  This version will determine the input based on values, outputs FRF, lat/lon,
%  and state plane coordinates.  Uses NAD83-2011.
%
%  IO:
%  p1 = FRF X (m), or Longitude (deg + or -), or state plane Easting (m) 
%  p2 = FRF Y (m), or Latitude (deg), or state plane Northing (m)
%
%  X = FRF cross-shore (m)
%  Y = FRF longshore (m)
%  ALat = latitude (decimal degrees)
%  ALon = longitude (decimal degrees, positive, or W)
%  spN = state plane northing (m)
%  spE = state plane easting (m)

% NAD83-86	2014 
% Origin Latitude          36.1775975
% Origin Longitude         75.7496860
% m/degLat             110963.357 
% m/degLon              89953.364 
% GridAngle (deg)          18.1465
% Angle FRF to Lat/Lon     71.8535
% Angle FRF to State Grid  69.9747
% FRF Origin Northing  274093.1562 
% Easting              901951.6805 

%  Debugging values
% p1=566.93;  p2=515.11;  % south rail at 1860
% ALat = 36.1836000
% ALon = 75.7454804
% p2= 36.18359977;  p1=-75.74548109;
% SP:  p1 = 902307.92; 	p2 = 274771.22; 


%  ************************************************************************

r2d = 180.0 / pi;

Eom=901951.6805;               % E Origin State Plane
Nom=274093.1562;               % N Origin State Plane
%ALat0=10.65583950;            % Origin Lat minutes
%ALon0=44.9811435;             % Origin Lon minutes
ALat0=36.1775975;              % Origin Lat minutes
ALon0=75.7496860;              % Origin Lon minutes
DegLat = 110963.35726;         % m/deg Lat
DegLon = 89953.36413;          % m/deg long
GridAngle=18.1465 ./r2d;
spAngle = (90 - 69.974707831) ./ r2d;

% Determine data type
if (floor(abs(p1)) == 75 & floor(p2) == 36)        % lat/lon input
	% to FRF
	ALat=p1;  ALon=p2;
	if (p1 < 0); p1 = -p1; end;
	ALatLeng = (p2 - ALat0) * DegLat;
	ALonLeng = -(p1 - ALon0) * DegLon;
	R = sqrt(ALatLeng.^2 + ALonLeng.^2 );
	Ang1 = atan2(ALonLeng, ALatLeng);
	Ang2 = Ang1 + GridAngle;
	X = R .* sin(Ang2);
	Y = R .* cos(Ang2);
	
	% to state plane 
	Ang2 = Ang2 - spAngle;
	AspN = R .* cos(Ang2);
	AspE = R .* sin(Ang2);             
	spN = AspN + Nom;
	spE = AspE + Eom;
	
elseif (p1 > 800000 & p2 > 200000)                % stat plane input
	spE=p1;  spN=p2;
	% to FRF
	spLengE = p1 - Eom;
	spLengN = p2 - Nom;
	R = sqrt(spLengE.^2 + spLengN.^2 );
	Ang1 = atan2(spLengE, spLengN);
	Ang2 = Ang1 + spAngle;
	X = R .* sin(Ang2);
	Y = R .* cos(Ang2);
	
	% to lat/lon
	Ang2 = Ang1 - (GridAngle-spAngle);            % 
	ALatLeng = R .* cos(Ang2);
	ALonLeng = R .* sin(-Ang2);         % neg to go west
	ALat = ALatLeng./DegLat + ALat0;    % was 60 * ALatLeng./DegLat + ALat0;
	ALon = ALonLeng./DegLon + ALon0;

elseif (p1 > -10000 & p1 < 10000 & p2 > -10000 & p2 < 10000)        % FRF input (+/- 10km)
	% to lat/lon
	X=p1;  Y=p2;
	R = sqrt(p1.^2 + p2.^2);
	Ang1 = atan2(p1, p2);               % CW from Y
	Ang2 = Ang1 - GridAngle;            % 
	ALatLeng = R .* cos(Ang2);
	ALonLeng = R .* sin(-Ang2);         % neg to go west
	ALat = ALatLeng./DegLat + ALat0;    % was 60 * ALatLeng./DegLat + ALat0;
	ALon = ALonLeng./DegLon + ALon0;

	% to state plane 
	Ang2 = Ang1 - spAngle;
	AspN = R .* cos(Ang2);
	AspE = R .* sin(Ang2);             
	spN = AspN + Nom;
	spE = AspE + Eom;
	
else
	disp '<EE> Cound not determine input type, returning NaNs'
	ALat=NaN; ALon=NaN; spN=NaN; spE=NaN; Y=NaN; X=NaN;
end


return
