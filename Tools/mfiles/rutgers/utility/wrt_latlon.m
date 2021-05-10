function wrt_latlon (lon, lat, Oname)

% WRT_LATLON:  Write (lat,lon) pairs into an ASCII file.
%
% wrt_latlon(lon, lat, Oname)  
%
% This function writes the (lat,lon) pairs into an ASCII so it may be used
% in the ROMS plotting package.
%
% On Input:
%
%    lon        Longitude coordinate (vector)
%    lat        Latitude  coordinate (vector)
%    Oname      Output ASCII filename (string)
%

% svn $Id: wrt_latlon.m 996 2020-01-10 04:28:56Z arango $
%=========================================================================%
%  Copyright (c) 2002-2020 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%  

% Write out (lat,lon) pairs.  Notice that in ROMS plotting file the
% latitude is written in the first column and longitude in the second
% column. The special value is a separator for pen up when drawing. 
  
spv=999;
fid=fopen(Oname,'w');
if (fid ~= -1)
  fprintf(fid, '%11.6f  %11.6f\n', spv, spv);
  for i=1:length(lon)
    fprintf(fid, '%11.6f  %11.6f\n', lat(i), lon(i));
  end
  fprintf(fid, '%11.6f  %11.6f\n', spv, spv);
  fclose(fid);
end

return
