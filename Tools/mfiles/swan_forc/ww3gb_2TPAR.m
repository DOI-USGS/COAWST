function ww3gb_2TPAR(modelgrid,yearww3,mmww3,ww3_area,ddww3,specres)

%this assumes data is historical data and is already downloaded to working
%directory
eval(['dpname=''',ww3_area,'.dp.',yearww3,mmww3,'.grb2'';']);
eval(['hsname=''',ww3_area,'.hs.',yearww3,mmww3,'.grb2'';']);
eval(['tpname=''',ww3_area,'.tp.',yearww3,mmww3,'.grb2'';']);

hsnc=ncgeodataset(hsname);
tpnc=ncgeodataset(tpname);
dpnc=ncgeodataset(dpname);
xg=dpnc{'lon'}(:);%lat of ww3 data
xg(xg>0)=xg(xg>0)-360;
yg=dpnc{'lat'}(:);%lon of ww3 data
[xg,yg]=meshgrid(xg,yg);
time=dpnc{'time'}(:);%time interval 3 hours
time2=datenum(str2num(yearww3),str2num(mmww3),str2num(ddww3)+1):3/24:datenum(str2num(yearww3),str2num(mmww3),str2num(ddww3)+((length(time)-1)/(24/3))+1);
time=time2;
%  limit the length of time.  Could use better logic here
%  for user to set this in ww3_swan_input file.
tstart=1; tend=length(time);
time=time(tstart:tend);

%determine spec pts from grid
% specpoints assumes a masking of 0 for land and NaN for water
[specpts]=ww3_specpoints(modelgrid,specres);

for i=1:length(specpts)
    gx=specpts(i,1);
    gy=specpts(i,2);
    xl=gx-3;xr=gx+3;
    yb=gy-3;yt=gy+3;

    %determine dimensions of ww3 data needed for interpolation
    [ym xm]=size(xg);
    ii=find(xg>=xl & xg<=xr & yg>=yb & yg<=yt);jj=fix((ii-1)/ym)+1;ii=mod(ii,ym);
    ig0=(min(ii)-1); ig1=(max(ii)+1); jg0=(min(jj)-1); jg1=(max(jj)+1);
    irg=[num2str(ig0) ':' num2str(ig1) ]; % irg='[1671:2042]';
    jrg=[num2str(jg0) ':' num2str(jg1) ]; % jrg='[2318:2722]';
    daplon=xg(ig0:ig1,jg0:jg1); daplat=yg(ig0:ig1,jg0:jg1);
    clear ii jj;
    
    %Interpolate the data to each point and create/write TPAR file
    for wavet=1:length(time)
%       eval(['hst=squeeze(hsnc{''Significant_height_of_combined_wind_waves_and_swell''}(wavet,',irg,',',jrg,'));']);
        eval(['hst=squeeze(hsnc{''Significant_height_of_combined_wind_waves_and_swell_surface''}(wavet,',irg,',',jrg,'));']); 
        zz=hst>1000;
        hst(zz)=0; %make bad data 0, swan not like NaNs
        Z1=interp2(daplon,daplat,hst,gx,gy);
        Z1(isnan(Z1))=0;
        TPAR(wavet,2)=Z1;
    end
%
    for wavet=1:length(time)
%       eval(['tpt=squeeze(tpnc{''Primary_wave_mean_period''}(wavet,',irg,',',jrg,'));']);
        eval(['tpt=squeeze(tpnc{''Primary_wave_mean_period_surface''}(wavet,',irg,',',jrg,'));']); 
        zz=tpt>1000;
        tpt(zz)=0; %make bad data 0, swan not like NaNs
        Z1=interp2(daplon,daplat,tpt,gx,gy);
        Z1(isnan(Z1))=0;
        TPAR(wavet,3)=Z1;
    end
%
    for wavet=1:length(time)
%       eval(['dpt=squeeze(dpnc{''Primary_wave_direction''}(wavet,',irg,',',jrg,'));']);
        eval(['dpt=squeeze(dpnc{''Primary_wave_direction_degree_true_surface''}(wavet,',irg,',',jrg,'));']); 
        zz=dpt>1000;
        dpt(zz)=0; %make bad data 0, swan not like NaNs
        Z1=interp2(daplon,daplat,dpt,gx,gy);
        Z1(isnan(Z1))=0;
        TPAR(wavet,4)=Z1;
    end
%
    TPAR(1:length(time),5)=20;
    TPAR_time(1:length(time),1)=str2num(datestr(time,'yyyymmdd.HHMM'));
    l=num2str(i);
    ofile=['TPAR',l,'.txt'];
    fid=fopen(ofile,'w');
    fprintf(fid,'TPAR \n');
%   fprintf(fid,'%8.4f         %3.2f        %3.2f     %3.f.       %2.f\n',TPAR');
    for wavet=1:length(time)
      fprintf(fid,'%8.4f',TPAR_time(wavet));
      fprintf(fid,'         %3.2f        %3.2f     %3.f.       %2.f\n',TPAR(wavet,2:5)');
    end
    fclose(fid);
end
close(dpnc)
close(hsnc);
close(tpnc);

