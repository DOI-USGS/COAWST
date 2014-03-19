function wrf2roms_mw(theWRFFile, theROMSFile)

% wrf2roms_mw('theWRFFile', 'theROMSFile')
% wrf2roms_mw -- Convert grid data from a WRF file to a ROMS file
% using the native mathworks (mw) netcdf interface.
 
% Version 05-Feb-2012 John C. Warner

% open the wrf file

rho.lat=ncread(theWRFFile,'XLAT');
rho.lon=ncread(theWRFFile,'XLONG');
rho.depth=zeros(size(rho.lon))+100;
rho.mask=1.0-ncread(theWRFFile,'LANDMASK');
spherical='T';
%projection='lambert conformal conic';
projection='mercator';


%call generic grid creation
save temp_jcw33.mat
eval(['mat2roms_mw(''temp_jcw33.mat'',''',theROMSFile,''');'])
!del temp_jcw33.mat

%User needs to edit roms variables
disp(['     '])
disp(['Created roms grid -->   ',theROMSFile])
disp(['     '])
disp(['You need to edit depth and masking in ',theROMSFile])
disp(['     '])

end
