function [hour,min,sec]=s2hms(secs)
% S2HMS:  converts seconds to integer hour,minute,seconds
%  
%  Usage: 
%        [hour,min,sec]=s2hms(secs)
%
%  Rich Signell rsignell@usgs.gov
 
hour=floor(secs/3600);
min=floor(rem(secs,3600)/60);
sec=round(rem(secs,60));
