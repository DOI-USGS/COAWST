function [windstart,windend]=loaddap_narr_hind(modelgrid,dirt);
%loaddap_narr_hind uses NARR data to create .dat wind input for SWAN
%dirt should be date as string in format 'yyyymmdd'
%modelgrid should be nc file with fields lon_rho and lat_rho 

loc=['http://nomads.ncdc.noaa.gov:80/dods/NCEP_NARR_DAILY/',dirt(1,1:6),'/',dirt,'/narr-a_221_',dirt,'_0000_000'];
%NARR data is available for 1979-present, -220°E to -0.625°E, 0°N to 89.625°N
%http://nomads.ncdc.noaa.gov/dods/NCEP_NARR_DAILY/yyyymm/yyyymmdd/narr-a_221_yyyymmdd_0000_000
%check to see if data exists
try 
    
    %info=loaddap('-A',loc);
    
    loaddap([loc,'?lon']);    
    xg=lon;
    loaddap([loc,'?lat']);
    yg=lat;
    [xg yg]=meshgrid(xg,yg);
    
    %get time
    loaddap([loc,'?time']);
    time=time+datenum(0000,12,30,0,0,0);
    tn=length(time);
    trg=['[0:',num2str(length(time)-1),']'];
    %trg2=['[8:',num2str(7+length(time)),']'];
    
    %get grid dimensions
    ncg=netcdf(modelgrid);
    xl=min(min(ncg{'lon_rho'}(:)));xr=max(max(ncg{'lon_rho'}(:)));
    yb=min(min(ncg{'lat_rho'}(:)));yt=max(max(ncg{'lat_rho'}(:)));     

    %determine dimensions of data needed from wind grid
    disp('optimizing grid dimensions ...');
    [ym xm]=size(xg);
    ii=find(xg>=xl & xg<=xr & yg>=yb & yg<=yt);jj=fix((ii-1)/ym)+1;ii=mod(ii,ym);
    ig0=(min(ii)-1); ig1=(max(ii)+1); jg0=(min(jj)-1); jg1=(max(jj)+1);
    irg=['[' num2str(ig0-1) ':' num2str(ig1-1) ']'];
    jrg=['[' num2str(jg0-1) ':' num2str(jg1-1) ']'];
    lat=lat(ig0:ig1);
    lon=lon(jg0:jg1);
    
    loaddap([loc,'?ugrd10m.ugrd10m',num2str(trg) num2str(irg) num2str(jrg)]);
    loaddap([loc,'?vgrd10m.vgrd10m',num2str(trg) num2str(irg) num2str(jrg)]);
    
    %take out "bad data" so calculations will work correctly
    ugrd10m(ugrd10m>9999)=0;
    vgrd10m(vgrd10m>9999)=0;

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

    fid = fopen(['narr_',dirt,'.dat'],'w'); %name of data file to write
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
catch
    disp('Error, unable to make wind file');    
end



% % check
% fid = fopen(['narr_',dirt,'.dat']);
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








