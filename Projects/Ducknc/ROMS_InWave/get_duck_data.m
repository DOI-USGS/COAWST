% get duck data
%
% main site is https://chldata.erdc.dren.mil
% https://chldata.erdc.dren.mil/thredds/catalog/frf/catalog.html
%
% all data
% https://chldata.erdc.dren.mil/thredds/dodsC/frf/geomorphology/DEMs/surveyDEM/surveyDEM.ncml
%
% jcw 26May2026

%cd to a wroking dir
cd D:\models\COAWST_updates\COAWST_v3.9\Ducknc_xshore_test\Projects\Ducknc\ROMS_InWave

% get bathy data from Nov 11 2023
url='https://chldata.erdc.dren.mil/thredds/dodsC/frf/geomorphology/DEMs/surveyDEM/data/FRF_geomorphology_DEMs_surveyDEM_20231109.nc'
time=ncread(url,'time');
xFRF=ncread(url,'xFRF');
yFRF=ncread(url,'yFRF');
elev=ncread(url,'elevation');
lat=ncread(url,'latitude');
lon=ncread(url,'longitude');
%save bath.mat
%%%%%%

%here we make Jonswap forcing files
clear f
f(1)=0.038;
for mm=2:33
  f(mm)=f(mm-1)*1.1;
end
Ohm = 2*pi./f;
Hs=2; Tp=12;
[ S, Amp, Phase ] = jonswap_lip( Ohm, Hs, Tp);
figure
hold on
plot(f,S,'r+-')

%
%  ok, now write the data out into a 2D spec file for InWave.
% 
fid=fopen('duck_xshore.spc2d','w');
HEAD=['SWAN   1                                Swan standard spectral file, version'; ...
      '$   Data converted to a SWAN version 40.91A                                 '; ...
      '$   Project: Delilah      ;  run number:                                    '; ...
      'TIME                                    time-dependent data                 '; ...
      '     1                                  time coding option                  '; ...
      'LOCATIONS                               locations in x-y-space              '; ...
      '     1                                  number of locations                 '; ...
      '     950.0000    800.0000                                                   '; ...
      'AFREQ                                   absolute frequencies in Hz          '; ...
      '    33                                  number of frequencies               '];
for mm=1:size(HEAD,1)
  fprintf(fid,HEAD(mm,:),'%76s \n');
  fprintf(fid,'\n');
end
for mm=1:length(f)
  fprintf(fid,num2str(f(mm)),'%s \n');
  fprintf(fid,'\n');
end
fprintf(fid,'NDIR','%s\n');
fprintf(fid,'\n');
fprintf(fid,'    36','%s\n');
fprintf(fid,'\n');
for mm=1:36
  fprintf(fid,num2str(0.0+(mm-1)*10.),'%s \n');
  fprintf(fid,'\n');
end
%
HEAD2=['QUANT                                                                   '; ...
       '     1                                  number of quantities in table   '; ...
       'EnDens                                  energy densities in J/m2/Hz/degr'; ...
       'J/m2/Hz/degr                            unit                            '; ...
       '   -0.9900E+02                          exception value                 '];
for mm=1:size(HEAD2,1)
  fprintf(fid,HEAD2(mm,:),'%72s \n');
  fprintf(fid,'\n');
end
for mm=1:2
  if mm==1; datestr='20231101.000000                         date and time'; end
  if mm==2; datestr='20231111.000000                         date and time'; end
  fprintf(fid,datestr,'%s\n');
  fprintf(fid,'\n');
  fprintf(fid,'FACTOR','%s\n');
  fprintf(fid,'\n');
  fprintf(fid,'    1.0E+0','%s\n');
  fprintf(fid,'\n');
  for ii=1:33
    for jj=1:36
      if (jj==10)
        fprintf(fid,num2str(S(ii)*1025.*9.81/10),' %s ');      % m^2/Hz kg/m^3 m/s^2 1/deg = J/(m^2 Hz deg)
      else
        fprintf(fid,num2str(S(ii)*1025.*9.81*0.),' %s ');
      end
      fprintf(fid,' ',' %s ');
    end
    fprintf(fid,'\n');
  end
end
fclose(fid);

