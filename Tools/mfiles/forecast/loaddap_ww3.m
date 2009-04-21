function [lat,lon,time,dp,hs,tp,info,indicator]=loaddap_ww3;
cd D:\Temporary\COAWST

t=datestr(now-1,'yyyymmdd');
dirt=['wna',t];
loc=['http://nomads.ncep.noaa.gov:9090/dods/wave/wna/',dirt,'/wna',t,'_00z'];
%loc=['http://nomads.ncep.noaa.gov:9090/dods/wave/nww3/nww3',t,'/nww3',t,'_00z'];
%check to see if data is updated
[s,status]=urlread(loc);
if length(s)>1000
    % loc2=['http://nomads.ncep.noaa.gov:9090/dods/wave/wna/',dirt,'/wna',t,'_0 6z'];
    % to get header information:
    % data is downloaded same day at ~ 5 UTC/GMT
    %info=loaddap('-A',loc);
    %y=loaddap('-A',loc2]);

    lat=loaddap([loc,'?lat']);
    lon=loaddap([loc,'?lon']);
    time=loaddap([loc,'?time']);%time interval 3 hours, length 8 days from next day, not now

    %convert to datenum, data convention is days from 1 1 1 00
%     if str2num(dirt(1,8:9))<2
%         time=time.time+366;
        time=time.time+datenum(0000,12,30,0,0,0);
%     else
%         time=time.time+365;
%     end

    dp=loaddap([loc,'?dirpwsfc.dirpwsfc']);
    hs=loaddap([loc,'?htsgwsfc.htsgwsfc']);
    tp=loaddap([loc,'?perpwsfc.perpwsfc']);
    indicator=1;
else
    indicator=0;
    lat=[];
    lon=[];
    time=[];
    dp=[];
    hs=[];
    tp=[];
    info=[];
end

