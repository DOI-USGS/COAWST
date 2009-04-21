function ww3totparnew(ww3_dpgr,ww3_hsgr,ww3_tpgr,timestr,lat,lon,modelgrid,gridname,time);
cd D:\Temporary\COAWST\SWAN\forcing
lon=lon-360;
[lon2,lat2]=meshgrid(lon,lat);
%interpolate ww3 variables to USEast_grd
ncg=netcdf(modelgrid);
gx=ncg{'lon_rho'}(:);gx2=gx;%gx2(gx<0) = gx2(gx<0)+360;
gy=ncg{'lat_rho'}(:);
[eta xi]=size(gy);
mask=ncg{'mask_rho'}(:);
%mask(mask==0)=nan;
mask(isnan(mask))=1;
gx3=gx2.*mask;
gy3=gy.*mask;
rr=1;
spec_res=50;
for ct=1:spec_res:eta
    if gx3(ct,1)~=0
        if ct+spec_res<eta
            specpts(rr,1)=gx2(ct,1)+0.001;
            specpts(rr,2)=gy(ct,1)+0.001;
            specpts(rr,3)=gx2(ct+spec_res,1)+0.001;
            specpts(rr,4)=gy(ct+spec_res,1)+0.001;
        else
            specpts(rr,1)=gx2(ct,1)+0.001;
            specpts(rr,2)=gy(ct,1)+0.001;
            specpts(rr,3)=gx2(end,1)+0.001;
            specpts(rr,4)=gy(end,1)+0.001;
        end
        rr=rr+1;
    end

end
for ct=1:spec_res:eta
    if gx3(ct,end)~=0
        if ct+spec_res<eta
            specpts(rr,1)=gx2(ct,end)-0.001;
            specpts(rr,2)=gy(ct,end)-0.001;
            specpts(rr,3)=gx2(ct+spec_res,end)-0.001;
            specpts(rr,4)=gy(ct+spec_res,end)-0.001;
        else
            specpts(rr,1)=gx2(ct,end)-0.001;
            specpts(rr,2)=gy(ct,end)-0.001;
            specpts(rr,3)=gx2(end,end)-0.001;
            specpts(rr,4)=gy(end,end)-0.001;
        end
        rr=rr+1;
    end

end
for ct=1:spec_res:xi
    if gx3(1,ct)~=0
        if ct+spec_res<xi
            specpts(rr,1)=gx2(1,ct);
            specpts(rr,2)=gy(1,ct)+0.001;
            specpts(rr,3)=gx2(1,ct+spec_res);
            specpts(rr,4)=gy(1,ct+spec_res)+0.001;
        else
            specpts(rr,1)=gx2(1,ct);
            specpts(rr,2)=gy(1,ct)+0.001;
            specpts(rr,3)=gx2(1,end)-0.001;
            specpts(rr,4)=gy(1,end)+0.001;
        end
        rr=rr+1;
    end

end
for ct=1:spec_res:xi
    if gx3(end,ct)~=0
        if ct+spec_res<xi
            specpts(rr,1)=gx2(end,ct)-0.001;
            specpts(rr,2)=gy(end,ct);
            specpts(rr,3)=gx2(end,ct+spec_res)-0.001;
            specpts(rr,4)=gy(end,ct+spec_res);
        else
            specpts(rr,1)=gx2(end,ct)-0.001;
            specpts(rr,2)=gy(end,ct);
            specpts(rr,3)=gx2(end,end)-0.001;
            specpts(rr,4)=gy(end,end);
        end
        rr=rr+1;
    end

end
clear gx gx2 gy gx3 gy3;
ncclose;
[m n]=size(time);
%m=24;%set this to only use partial data set, 25 is ~3 days (@ 3hr interval)
%time loop to create time series for each spec point
%input longitude and latitude in decimal degrees

[r c]=size(specpts);
%timestr=datestr(time(1:m,n),'yyyymmdd.HHMM');
%take out bad wave data so that calculations will work
ww3_dpgr(ww3_dpgr>9999)=0;
ww3_hsgr(ww3_hsgr>9999)=0;
ww3_tpgr(ww3_tpgr>9999)=0;
save specpts.mat specpts
for pt=1:r
    for wavet=1:m
        hs=squeeze(ww3_hsgr(:,:,wavet));
        zz=hs>1000;
        hs(zz)=0; %make bad data 0, swan not like NaNs
        %Z1=interp2(lon2,lat2,hs,specpts(pt,1),specpts(pt,2));
        Z1=griddata(lon2,lat2,hs,specpts(pt,1),specpts(pt,2));
        TPAR(wavet,2)=Z1;
    end
    for wavet=1:m
        tp=squeeze(ww3_tpgr(:,:,wavet));
        zz=tp>1000;
        tp(zz)=0; %make bad data 0, swan not like NaNs
        %Z1=interp2(lon2,lat2,tp,specpts(pt,1),specpts(pt,2));
        Z1=griddata(lon2,lat2,tp,specpts(pt,1),specpts(pt,2));
        TPAR(wavet,3)=Z1;
    end
    for wavet=1:m
        dp=squeeze(ww3_dpgr(:,:,wavet));
        zz=dp>1000;
        dp(zz)=0; %make bad data 0, swan not like NaNs
        %Z1=interp2(lon2,lat2,dp,specpts(pt,1),specpts(pt,2));
        Z1=griddata(lon2,lat2,dp,specpts(pt,1),specpts(pt,2));
        TPAR(wavet,4)=Z1;
    end
    TPAR(1:m,1)=str2num(timestr(1:m,:));
    TPAR(1:m,5)=20;
    l=num2str(pt);
    ofile=['TPAR',l,'.txt'];
    fid=fopen(ofile,'w');
    fprintf(fid,'TPAR \n');
    fprintf(fid,'%8.4f         %3.2f        %3.2f     %3.f.       %2.f\n',TPAR');
    fclose(fid);

end

