% create_swan_2dspec_from_ww3
%
% jcw 04April2025
%

%0) cd to working dir
cd D:\Projects\NOPP_hurricanes_modeling\Ian_2022\ian150

%1) load ww3 file
netcdf_load('ww3.202209_spec_20220928_12to13.nc')
%netcdf_load('ww3.202209_spec_20220928_14.nc')
%netcdf_load('ww3.202209_spec_20220928_15.nc')
%netcdf_load('ww3.202209_spec_20220928_16to22.nc')

%2) what is name of new pseudo swan spec file
spec_file='SWAN2D_sanible1.txt';
%spec_file='SWAN2D_sanible2.txt';
%spec_file='SWAN2D_sanible3.txt';
%spec_file='SWAN2D_sanible4.txt';

%what point number do you want to use
spec_pnt=1;

%%%%%%  end of user

%create text for SWAN spec file
fid=fopen(spec_file,'w');
nline=['SWAN   1                                Swan standard spectral file, version'];
  fprintf(fid,nline); 
  fprintf(fid,'\n');
nline=['$   Data produced by SWAN version 40.91A              '];
  fprintf(fid,nline); 
  fprintf(fid,'\n');
nline=['$   Project: Inlet Test      ;  run number:     '];
  fprintf(fid,nline); 
  fprintf(fid,'\n');
nline=['TIME                                    time-dependent data'];
  fprintf(fid,nline); 
  fprintf(fid,'\n');
nline=['     1                                  time coding option'];
  fprintf(fid,nline); 
  fprintf(fid,'\n');
nline=['LOCATIONS                               locations in x-y-space'];
  fprintf(fid,nline); 
  fprintf(fid,'\n');
nline=['     1                                  number of locations'];
  fprintf(fid,nline); 
  fprintf(fid,'\n');
nline=['     ',num2str(longitude(spec_pnt)),'     ',num2str(latitude(spec_pnt))];
  fprintf(fid,nline); 
  fprintf(fid,'\n');
nline=['AFREQ                                   absolute frequencies in Hz'];
  fprintf(fid,nline); 
  fprintf(fid,'\n');
nline=['    ',num2str(length(frequency)),'                                  number of frequencies'];
  fprintf(fid,nline); 
  fprintf(fid,'\n');
for mm=1:length(frequency)
  nline=['    ',num2str(frequency(mm))];
  fprintf(fid,nline); 
  fprintf(fid,'\n');
end
nline=['NDIR                                    spectral nautical directions in degr'];
  fprintf(fid,nline); 
  fprintf(fid,'\n');
  NDIR=num2str(length(direction));
nline=['    ',NDIR,'                                  number of directions'];
  fprintf(fid,nline); 
  fprintf(fid,'\n');
for mm=1:length(direction)
  nline=['    ',num2str(direction(mm))];
  fprintf(fid,nline); 
  fprintf(fid,'\n');
end
nline=['QUANT'];
  fprintf(fid,nline); 
  fprintf(fid,'\n');
nline=['     1                                  number of quantities in table'];
  fprintf(fid,nline); 
  fprintf(fid,'\n');
nline=['EnDens                                  energy densities in J/m2/Hz/degr'];
  fprintf(fid,nline); 
  fprintf(fid,'\n');
nline=['J/m2/Hz/degr                            unit'];
  fprintf(fid,nline); 
  fprintf(fid,'\n');
nline=['   -0.9900E+02                          exception value'];
  fprintf(fid,nline); 
  fprintf(fid,'\n');
%
wave_time=time+datenum(1990,01,01,0,0,0);
for tt=1:length(wave_time)
  Z=datevec(wave_time(tt));
  YY=Z(:,1);  YY=[num2str(YY)];
  MM=Z(:,2);  MM=['0',num2str(MM)]; MM=MM(end-1:end);
  DD=Z(:,3);  DD=['0',num2str(DD)]; DD=DD(end-1:end);
  hh=Z(:,4);  hh=['0',num2str(hh)]; hh=hh(end-1:end);
  mm=Z(:,5);  mm=['0',num2str(mm)]; mm=mm(end-1:end);
  ss=Z(:,6);  ss=['0',num2str(ss)]; ss=ss(end-1:end);
  nline=[YY,MM,DD,'.',hh,mm,ss,                        '                         date and time'];
    fprintf(fid,nline); 
    fprintf(fid,'\n');
  nline=['FACTOR'];
    fprintf(fid,nline); 
    fprintf(fid,'\n');
  nline=['    1.0'];
    fprintf(fid,nline); 
    fprintf(fid,'\n');
  % wirte out spec at this time at this point
  %  1 J = kg m2 /s2
  % efth is m2 s rad-1, convert to m2/hz/deg:  2pirad/360deg  1Hz 1s
  %
  %   m2 s  *  2 pi rad    1       * 1025 kg   9.81 m     = (2*pi/360)*1025*9.81    kg    * (m2)  = 2*pi*1025*9.81/360 J 
  %   rad      360 deg   1Hz 1s           m3       s2                           deg s2 Hz   (m2)                   deg Hz m2
  %
  fac=2*pi*1025*9.81/360;
  %efth(dir freq #sta #times)
  for ii=1:length(frequency)
    fprintf(fid, '%12.4f  ',squeeze(efth(:,ii,spec_pnt,tt))'*fac);
    fprintf(fid,'\n');
  end
end
fprintf(fid,'\n');
fclose(fid);
