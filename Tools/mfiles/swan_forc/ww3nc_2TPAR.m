function ww3nc_2TPAR(modelgrid,yearww3,mmww3,working_dir,working_drive,ww3_area,ddww3)

t=[yearww3,mmww3,ddww3];
dirt=[ww3_area,t];
loc=['http://nomads.ncep.noaa.gov:9090/dods/wave/',ww3_area,'/',dirt,'/',dirt,'_00z'];
%loc=['http://nomads.ncep.noaa.gov:9090/dods/wave/nww3/nww3',t,'/nww3',t,'_00z'];
%check to see if data is updated
[s,status]=urlread(loc);

if length(s)>1000
    xg=loaddap([loc,'?lon']);%load ww3 lon
    yg=loaddap([loc,'?lat']);%load ww3 lat
    yg=yg.lat;
    xg=xg.lon-360;
    [xg,yg]=meshgrid(xg,yg);
    time=loaddap([loc,'?time']);%time interval 3 hours
    
    %determine spec pts from grid
    % specpoints assumes a masking of 0 for land and NaN for water
    [specpts]=ww3_specpoints(modelgrid,50);
    
    for i=1:length(specpts)
        gx=specpts(i,1);
        gy=specpts(i,2);
        xl=gx-3;xr=gx+3;
        yb=gy-3;yt=gy+3;

        %determine dimensions of ww3 data needed for interpolation
        [ym xm]=size(xg);
        ii=find(xg>=xl & xg<=xr & yg>=yb & yg<=yt);jj=fix((ii-1)/ym)+1;ii=mod(ii,ym);
        ig0=(min(ii)-1); ig1=(max(ii)+1); jg0=(min(jj)-1); jg1=(max(jj)+1);
        irg=['[' num2str(ig0) ':' num2str(ig1) ']']; % irg='[1671:2042]';
        jrg=['[' num2str(jg0) ':' num2str(jg1) ']']; % jrg='[2318:2722]';
        dap.lon=xg(ig0:ig1,jg0:jg1); dap.lat=yg(ig0:ig1,jg0:jg1);
        clear ii jj;

        %download the data using loaddap  
        dp=loaddap([loc '?dirpwsfc.dirpwsfc[0:' num2str(length(time.time)-1) ']' irg jrg ]);
        dp=dp.dirpwsfc;

        hs=loaddap('+v',[loc '?htsgwsfc.htsgwsfc[0:' num2str(length(time.time)-1) ']' irg jrg ]);
        hs=hs.htsgwsfc;

        tp=loaddap('+v',[loc '?perpwsfc.perpwsfc[0:' num2str(length(time.time)-1) ']' irg jrg ]);
        tp=tp.perpwsfc;
        indicator=1;
        
        %Interpolate the data to each point and create/write TPAR file
        for wavet=1:length(time.time)
                hst=squeeze(hs(:,:,wavet));
                zz=hst>1000;
                hst(zz)=0; %make bad data 0, swan not like NaNs
                %Z1=interp2(lon2,lat2,hs,specpts(r,1),specpts(r,2));
                Z1=griddata(dap.lon,dap.lat,hst,specpts(i,1),specpts(i,2));
                TPAR(wavet,2)=Z1;
            end
            for wavet=1:length(time.time)
                tpt=squeeze(tp(:,:,wavet));
                zz=tpt>1000;
                tpt(zz)=0; %make bad data 0, swan not like NaNs
                %Z1=interp2(lon2,lat2,tp,specpts(r,1),specpts(r,2));
                Z1=griddata(dap.lon,dap.lat,tpt,specpts(i,1),specpts(i,2));
                TPAR(wavet,3)=Z1;
            end
            for wavet=1:length(time.time)
                dpt=squeeze(dp(:,:,wavet));
                zz=dpt>1000;
                dpt(zz)=0; %make bad data 0, swan not like NaNs
                %Z1=interp2(lon2,lat2,dp,specpts(r,1),specpts(r,2));
                Z1=griddata(dap.lon,dap.lat,dpt,specpts(i,1),specpts(i,2));
                TPAR(wavet,4)=Z1;
            end
            TPAR(1:length(time.time),1)=str2num(datestr(time.time+datenum(0000,12,30,0,0,0),'yyyymmdd.HHMM'));
            TPAR(1:length(time.time),5)=20;
            l=num2str(i);
            ofile=['TPAR',l,'.txt'];
            fid=fopen(ofile,'w');
            fprintf(fid,'TPAR \n');
            fprintf(fid,'%8.4f         %3.2f        %3.2f     %3.f.       %2.f\n',TPAR');
            fclose(fid);
    end
end

