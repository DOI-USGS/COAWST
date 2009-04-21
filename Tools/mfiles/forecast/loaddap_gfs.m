function [windstart,windend,indicator]=loaddap_gfs(modelgrid);
%0.5x0.5 degree http://nomads.ncep.noaa.gov:9090/dods/gfs_hd/gfs_hd20090323/gfs_hd_00z
%Longitude: 0.00000000000°E to 359.50000000000°E  (720 points, avg. res. 0.5°)  
%Latitude: -90.00000000000°N to 90.00000000000°N  (361 points, avg. res. 0.5°)  
%Time: 61 points, avg. res. 0.125 days
cd D:\Temporary\COAWST\SWAN\forcing
t=datestr(now-1,'yyyymmdd');
dirt=[t];
loc=['http://nomads.ncep.noaa.gov:9090/dods/gfs_hd/gfs_hd',dirt,'/gfs_hd_00z'];
%check to see if data is updated
s=urlread(loc);
if length(s)>1000 %data is updated
    clear s
    %info=loaddap('-A',loc);
    loaddap([loc,'?ugrd10m.ugrd10m[0:60][201:281][515:617]']);%[:][201:281][515:617]');
    loaddap([loc,'?vgrd10m.vgrd10m[0:60][201:281][515:617]']);%[:][201:281][515:617]');
    
    %take out "bad data" so calculations will work correctly
    ugrd10m(ugrd10m>9999)=0;
    vgrd10m(vgrd10m>9999)=0;
    
    loaddap([loc,'?vgrd10m.time']);%[4:28]
    loaddap([loc,'?vgrd10m.lat[201:281]']);%[0:380]
    loaddap([loc,'?vgrd10m.lon[515:617]']);%[480:1037]
    tn=length(time);
%     if str2num(dirt(1,5:6))<2
%         time_gfs=time+366;
        time_gfs=time+datenum(0000,12,30,0,0,0);
%     else
%         time_gfs=time+365;%convert to datenum, data convention is days from 1 1 1 00
%     end
    windstart=datestr(time_gfs(1),'yyyymmdd.HHMMSS');
    windend=datestr(time_gfs(end),'yyyymmdd.HHMMSS');
    [x_gfs y_gfs]=meshgrid(lon,lat);
    clear lon lat time

    %Grid to interpolate to
    ncg=netcdf(modelgrid);
    gx=ncg{'lon_rho'}(:);gx2=gx;gx2(gx<0) = gx2(gx<0)+360;clear gx;
    gy=ncg{'lat_rho'}(:);
    [eta xi]=size(ncg{'lon_rho'}(:));
    angle=ncg{'angle'}(:);
    ncclose

    %interpolate to grid
    ugrid=ones(eta,xi,tn);
    vgrid=ones(eta,xi,tn);
    for i=1:tn
        ugrid(:,:,i)=interp2(x_gfs,y_gfs,ugrd10m(:,:,i),gx2,gy);
        vgrid(:,:,i)=interp2(x_gfs,y_gfs,vgrd10m(:,:,i),gx2,gy);
    end
    clear ugrd10m vgrd10m x_gfs y_gfs gx2 gy
    
    %rotate to grid
    [ugridrot, vgridrot]=rotation(ugrid, vgrid, angle);

    clear ugrid vgrid

    %get rid of problematic NANs
    ugridrot(isnan(ugridrot)==1)=0;
    vgridrot(isnan(vgridrot)==1)=0;

    clear ugrid vgrid
    % write wind to ascii file

    fid = fopen('gfswind_update.dat','w'); %name of data file to write
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
else %data is not updated
    %determine wind start and wind end
    f=textread('D:\Temporary\COAWST\SWAN\run\INPUT','%s');
    ind=find(strcmp(f,'NONSTATIONARY'),2);
    windstart=char(f(ind(2)+1));
    windend=[num2str(str2num(windstart)+7.12) '0000'];
    indicator=0;
end
%1.0x1.0 degree http://nomads.ncep.noaa.gov:9090/dods/gfs/gfs20090323/gfs_00z
%2.5x2.5 degree
%http://nomads.ncep.noaa.gov:9090/dods/gfs_2p5/gfs_2p520090323/gfs2p5_00z



