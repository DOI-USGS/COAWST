function ww3gb_2TPAR(modelgrid,yearww3,mmww3,ww3_grid,specpts)

%set urls of the hsig, peak period, and dominant period
if (str2num(yearww3)<2017)
  eval(['hsurl=''https://www.ncei.noaa.gov/thredds-ocean/dodsC/ncep/nww3/',yearww3,'/',mmww3,'/',ww3_grid,'/multi_1.',ww3_grid,'.hs.',yearww3,mmww3,'.grb2'';'])
  eval(['tpurl=''https://www.ncei.noaa.gov/thredds-ocean/dodsC/ncep/nww3/',yearww3,'/',mmww3,'/',ww3_grid,'/multi_1.',ww3_grid,'.tp.',yearww3,mmww3,'.grb2'';'])
  eval(['dpurl=''https://www.ncei.noaa.gov/thredds-ocean/dodsC/ncep/nww3/',yearww3,'/',mmww3,'/',ww3_grid,'/multi_1.',ww3_grid,'.dp.',yearww3,mmww3,'.grb2'';'])
else
  eval(['hsurl=''https://www.ncei.noaa.gov/thredds-ocean/dodsC/ncep/nww3/',yearww3,'/',mmww3,'/gribs/multi_1.',ww3_grid,'.hs.',yearww3,mmww3,'.grb2'';'])
  eval(['tpurl=''https://www.ncei.noaa.gov/thredds-ocean/dodsC/ncep/nww3/',yearww3,'/',mmww3,'/gribs/multi_1.',ww3_grid,'.tp.',yearww3,mmww3,'.grb2'';'])
  eval(['dpurl=''https://www.ncei.noaa.gov/thredds-ocean/dodsC/ncep/nww3/',yearww3,'/',mmww3,'/gribs/multi_1.',ww3_grid,'.dp.',yearww3,mmww3,'.grb2'';'])
end
%
%first lets get the lon, lat, and time of the ww3 data
disp(' getting ww3 lon lat and time')
lonww3_1d=double(ncread(hsurl,'lon'));  zlon=length(lonww3_1d);
latww3_1d=double(ncread(hsurl,'lat'));  zlat=length(latww3_1d);
%latww3_1d=flipud(latww3_1d);
lonww3=repmat(lonww3_1d,1,zlat);
latww3=repmat(latww3_1d,1,zlon)';
%latww3=fliplr(latww3);
%
timeww3=ncread(hsurl,'time');   % Hour since the start of the month, this includes first hour of next month so
timeww3=timeww3(1:end-1);       % we will stop 1 short to not allow double count of end of this month and start of next.
%
% put time in a datenum format
for mm=1:length(timeww3)
  time(mm)=datenum(str2num(yearww3),str2num(mmww3),1,timeww3(mm),0,0);
end
%
% Now prep the swan/roms grid location arrays
%
offset=0;
if (min(specpts(:,1)<0)); offset=360; end
disp(['We are adding ',num2str(offset),' to the roms longitude - line 34 ww3gb_2TPAR.m'])
gx=specpts(:,1)+offset;
gy=specpts(:,2);
xl=min(gx(:))-3;xr=max(gx(:))+3;     %this adds a buffer to make sure we get enough data.
yb=min(gy(:))-3;yt=max(gy(:))+3;
%
% determine dimensions of ww3 data needed for interpolation
%
ii=find(lonww3_1d>=xl & lonww3_1d<=xr);
ig0=(min(ii)-1); ig1=(max(ii)+1);
jj=find(latww3_1d>=yb & latww3_1d<=yt);
jg0=(min(jj)-1); jg1=(max(jj)+1);
daplon=lonww3(ig0:ig1,jg0:jg1); daplat=latww3(ig0:ig1,jg0:jg1);
clear ii jj;
%
% Now for each moment in time, we read the hs, tp, dp and then interpolate
% those fields to the specpoints.
%
method='linear';
for tidx=1:length(timeww3)
%
  disp(['getting hs tp and dp for ',datestr(time(tidx),'yyyymmdd.HHMM')])
% get hs
  hs=double(squeeze(ncread(hsurl,'Significant_height_of_combined_wind_waves_and_swell_surface',[ig0 jg0 tidx],[ig1-ig0+1 jg1-jg0+1 1])));
  zz=hs>1000;
  hs(zz)=0; %make bad data 0, swan not like NaNs
  hs(isnan(hs))=0;
% get tp
  tp=double(squeeze(ncread(tpurl,'Primary_wave_mean_period_surface',[ig0 jg0 tidx],[ig1-ig0+1 jg1-jg0+1 1])));
  zz=tp>1000;
  tp(zz)=0; %make bad data 0, swan not like NaNs
  tp(isnan(tp))=0;
% get dp
  dp=double(squeeze(ncread(dpurl,'Primary_wave_direction_surface',[ig0 jg0 tidx],[ig1-ig0+1 jg1-jg0+1 1])));
  zz=dp>1000;
  dp(zz)=0; %make bad data 0, swan not like NaNs
  dp(isnan(dp))=0;
%
% interpolate hs
  FCr = griddedInterpolant(daplon,fliplr(daplat),fliplr(hs),method);
  HS(:,tidx)=FCr(gx,gy);
% interpolate tp
  FCr.Values=fliplr(tp);
  TP(:,tidx)=FCr(gx,gy); 
% interpolate dp
  Dwave_Ax=1*cos((90-double(dp))*pi/180);%grid of x
  Dwave_Ay=1*sin((90-double(dp))*pi/180);%grid of y
  FCr.Values=fliplr(Dwave_Ax);
  Zx(:)=FCr(gx,gy);
  FCr.Values=fliplr(Dwave_Ay);
  Zy(:)=FCr(gx,gy);
  Z=90. - atan2(Zy,Zx)*180./pi; %calculate direction
  if Z<0
    Z=Z+360;
  end
  DP(:,tidx)=Z; 
end
%
% now write out the spec files
%
% set some min values
HS(HS<0.01)=0.01;
TP(TP<0.01)=0.01;
DP(DP<0.0)=DP(DP<0.0)+360.;
%
for i=1:length(specpts)
    disp(['writing out file TPAR',num2str(i)])
    TPAR_time=str2num(datestr(time,'yyyymmdd.HHMM'));
    ofile=['TPAR',num2str(i),'.txt'];
    fid=fopen(ofile,'w');
    fprintf(fid,'TPAR \n');
    for tidx=1:length(time)
      fprintf(fid,'%8.4f',TPAR_time(tidx));
      fprintf(fid,'         %3.2f        %3.2f     %3.f.       %2.f\n',[HS(i,tidx) TP(i,tidx) DP(i,tidx) 20]);
    end
    fclose(fid);
end

