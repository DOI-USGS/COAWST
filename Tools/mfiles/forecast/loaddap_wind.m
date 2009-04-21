function [windstart,windend,indicator]=loaddap_wind(modelgrid);
cd D:\Temporary\COAWST\SWAN\forcing
t=datestr(now-1,'yyyymmdd');
dirt=[t];
loc=['http://nomads.ncep.noaa.gov:9090/dods/nam/nam',dirt,'/nam_00z'];
%check to see if data is updated
s=urlread(loc);
if length(s)>1000
    %loc=['http://nomads.ncdc.noaa.gov:80/dods/NCEP_NAM/',dirt(1,1:6),'/',dirt,'/nam_218_',dirt,'_1200_fff'];
    %info=loaddap('-A',loc);
    loaddap([loc,'?ugrd10m.ugrd10m']);%[4:28][0:380][480:1037]');
    loaddap([loc,'?vgrd10m.vgrd10m']);%[4:28][0:380][480:1037]');
    
    %take out "bad data" so calculations will work correctly
    ugrd10m(ugrd10m>9999)=0;
    vgrd10m(vgrd10m>9999)=0;
    
    loaddap([loc,'?vgrd10m.time']);%[4:28]
    loaddap([loc,'?vgrd10m.lat']);%[0:380]
    loaddap([loc,'?vgrd10m.lon']);%[480:1037]
    tn=length(time);
%     if str2num(dirt(1,5:6))<2
%         time=time+366;
        time=time+datenum(0000,12,30,0,0,0);
%     else
%         time=time+365;%convert to datenum, data convention is days from 1 1 1 00
%     end
    windstart=datestr(time(1),'yyyymmdd.HHMMSS');
    windend=datestr(time(end),'yyyymmdd.HHMMSS');
    [xwind ywind]=meshgrid(lon,lat);

    %Grid to interpolate to
    ncg=netcdf(modelgrid);
    gx=ncg{'lon_rho'}(:);
    gy=ncg{'lat_rho'}(:);
    [eta xi]=size(gy);
    angle=ncg{'angle'}(:);
    ncclose

    %interpolate to grid
    ugrid=ones(eta,xi,tn);
    vgrid=ones(eta,xi,tn);
    for i=1:tn
        ugrid(:,:,i)=griddata(xwind,ywind,ugrd10m(:,:,i),gx,gy);
        vgrid(:,:,i)=griddata(xwind,ywind,vgrd10m(:,:,i),gx,gy);
    end

    %get rid of problematic NANs
    ugrid(isnan(ugrid)==1)=999;
    vgrid(isnan(vgrid)==1)=999;   
    ugrid(ugrid==0)=999;
    vgrid(vgrid==0)=999;   
   
    clear ugrd10m vgrd10m xwind ywind

    %rotate to grid
    [ugridrot, vgridrot]=rotation(ugrid, vgrid, angle);

    clear ugrid vgrid

    % write wind to ascii file

    fid = fopen('namwind_update.dat','w'); %name of data file to write
    for file=1:tn;
        for index = 1:eta;
            %for index2 = 1:xi;
                fprintf(fid,'%3.2f\n',ugridrot(index,:,file));
            %end
            %fprintf(fid,'\n');
        end
        for index = 1:eta;
            %for index2 = 1:xi;
                fprintf(fid,'%3.2f\n',vgridrot(index,:,file));
            %end
            %fprintf(fid,'\n');
        end
    end
    fclose(fid);
    indicator=1;
else
    %determine wind start and wind end
    f=textread('D:\Temporary\COAWST\SWAN\run\INPUT','%s');
    ind=find(strcmp(f,'NONSTATIONARY'),2);
     windstart=(char(f(ind(2)+1)));
     windend=(char(f(ind(2)+4)));
%     windstart=str2num(char(f(ind(2)+1)));
%     windend=str2num(char(f(ind(2)+4)));
    indicator=0;
end



%check
% fid = fopen('namwind_update.dat');
% for tidx=1:9999999
%   tidx
%   [A,count]=fscanf(fid,'%f',896*336);
%   if count==0; break; end
%   wind_temp=reshape(A,896,336);
%   wind_temp=wind_temp';
%   pcolorjw(wind_temp);
%   %write wind
%    %fprintf(fid2,'%f\n',winds');
%    pause(1);
% end
%fid =fopen('namwind_update.dat')
%x=fread(fid,336*896,'float')
%x=reshape(x,336,896);
%pcolorjw(x);







