function xfac=dasp(lat)
% DASP sets data aspect ratio to 1
% or if optional lat is supplied, crude projection
% Usage: xfac=dasp(lat); 
% where: lat = latitude in degrees
% xfac : the ratio lon/lat
if(nargin==1),
   if (isnan(lat)),
      set(gca,'DataAspectRatioMode','auto');
   else
      xfac=cos(lat*pi/180);
      set (gca, 'DataAspectRatio', [1 xfac 1] );
   end
else
   set (gca, 'DataAspectRatio', [1 1 1] );
end
