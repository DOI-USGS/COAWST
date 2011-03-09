% written by Mingkui Li, May 2008
% Modified by Brandy Armstrong March 2009
% Modified again by jcw, April 2009
% Url of the global climate data
% url='http://hycom.coaps.fsu.edu/thredds/dodsC/glb_analysis';
Time1=datestr(now+1,'yyyymmdd');
Time_before=datestr(now,'yyyymmdd');


%***********USER CHOOSE FSU or RTOFS ***************************
%url=['http://nomads.ncep.noaa.gov:9090/dods/ofs/ofs',Time1,'/daily/rtofs_native_000_atl']
%http://nomads.ncep.noaa.gov:9090/dods/ofs/ofs20090702/daily/rtofs_native_000_atl%experimental server
url=['http://nomads.ncep.noaa.gov:9090/dods/ofs/ofs',Time_before,'/daily/rtofs_native_024_atl'];
url1=['http://nomad1.ncep.noaa.gov:9090/dods/rtofs_native/daily/ofs.',Time_before,'/rtofs_DODS_atl'];
try
    url2='http://tds.hycom.org/thredds/dodsC/GLBa0.08/expt_90.9/2011';
    nc=mDataset(url2);
    lontest=(double(nc{'Longitude'}(1,:)));
    if size(lontest)==0
        clear lontest;
        close(nc);
        lever=0;
    else
        lever=1;
            clear lontest
        close(nc);
    end
catch
    lever=0;
    close(nc);
end
%***************************************************************

% get grid lat and lon
disp('getting roms grid dimensions ...');
xl=min(min(gn.lon_rho));xr=max(max(gn.lon_rho));
yb=min(min(gn.lat_rho));yt=max(max(gn.lat_rho));

%download hycom dimensions
disp('get hycom dimensions ...');
if lever==0;
    try
         nc=mDataset(url);
        lon=double(nc{'lon'}(:)); %NJ Toolbox way
%        lon=nj_varget(url,'lon'); %NJ Toolbox way
        %loaddap('+v',[url,'?lon']); %loaddap way
        %xg=Longitude.Longitude-360;clear Longitude;
        %lon=nc{'lon'}(:); %Mexnc way

        xg=lon;
        clear lon
        lat=double(nc{'lat'}(:)); %NJ Toolbox way
        
%         lat=nj_varget(url,'lat'); %NJ Toolbox way
        %loaddap('+v',[url,'?lat']); %loaddap way
        %     yg=Latitude.Latitude;clear Latitude;
        %lat=nc{'lat'}(:); %Mexnc way

        yg=lat;
        clear lat
        [xg yg]=meshgrid(xg,yg);
        eval(['tmp=double(nc{''wtmpcdsl''}(1,:,[547:1087],[24:391]));']);
        clear tmp;
    catch
        close(nc)
        display('using nomad1.ncep.noaa.gov:9090');
        nc=mDataset(url1);        
        lon=double(nc{'lon'}(:)); %NJ Toolbox way
%         lon=nj_varget(url1,'lon'); %NJ Toolbox way
        %loaddap('+v',[url,'?lon']); %loaddap way
        %xg=Longitude.Longitude-360;clear Longitude;

        %lon=nc{'lon'}(:); %Mexnc way

        xg=lon;
        clear lon
        lat=double(nc{'lat'}(:)); %NJ Toolbox way
%         lat=nj_varget(url1,'lat'); %NJ Toolbox way

        %loaddap('+v',[url,'?lat']); %loaddap way
        %     yg=Latitude.Latitude;clear Latitude;

        %lat=nc{'lat'}(:); %Mexnc way

        yg=lat;
        clear lat
        [xg yg]=meshgrid(xg,yg);
    end
end
if lever==1;
    %*****************************special fsu hycom
    load irg.mat %fsu data
    load jrg.mat %fsu data
    load hycom_info.mat

    irg2=irg;
    jrg2=jrg;

    irg=['[' num2str(str2num(irg(2:5))-1) ':' num2str(str2num(irg(7:10))-1) ']'];
    jrg=['[' num2str(str2num(jrg(2:5))-1) ':' num2str(str2num(jrg(7:10))-1) ']'];

    %******************************
%     clm.lon=xg(ig0:ig1,jg0:jg1); clm.lat=yg(ig0:ig1,jg0:jg1);
%     eval(['clm.lon=xg(',irg2,',',jrg2,');']);
%     eval(['clm.lat=yg(',irg2,',',jrg2,');']);

            display('using http://tds.hycom.org/thredds/dodsC/GLBa0.08/expt_90.8/2011');
        nc=mDataset(url2);        
%     lon=74.12003:((434.11993-74.12003)/4499):434.11993;
%     lon=double(nc{'Longitude'}(:));
     eval(['lon=double(nc{''Longitude''}(',irg2,',',jrg2,'));']);
    %lon=nj_varget(url,'Longitude'); %NJ Toolbox way,fsu data

    %loaddap('+v',[url,'?Longitude']); %fsu url var name
    %xg=Longitude.Longitude-360;clear Longitude;
    %lon=nc{'Longitude'}(:); %Mexnc way

    xg=lon;
%    xg(xg>=360)=-360+(xg(xg>=360)-360);
    xg(xg>=180)=(xg(xg>=180)-360);
%     xg(xg<=-180)=-(xg(xg<=-180)+180);
    clear lon
%    lat=-78.6400:((89.98425-(-78.6400))/3297):89.98425;
    %lat=double(nc{'Latitude'}(:));
     eval(['lat=double(nc{''Latitude''}(',irg2,',',jrg2,'));']);
    %lat=nj_varget(url,'Latitude'); %NJ Toolbox way, fsu data

    %    loaddap('+v',[url,'?Latitude']); %fsu url var name
    %     yg=Latitude.Latitude;clear Latitude;
    %lat=nc{'Latitude'}(:); %Mexnc way
    yg=lat;
    clear lat
%     [xg yg]=meshgrid(xg,yg);
    
end

%loaddap('+v',[url,'?lev']); %loaddap way
if lever==0;
    lev=double(nc{'lev'}(:));
    %lev=nj_varget(url,'lev'); %NJ Toolbox way, rtofs lev is in m
    %lev=nc{'lev'}(:); %Mexnc way
    %    loaddap('+v',[url,'?lev']); %fsu url var name
    clm.z=lev; clear lev;
end
if lever==1;
    lev=double(nc{'Depth'}(:));
    %lev=nj_varget(url,'Depth'); %NJ Toolbox way, fsu data depth is in m
    %lev=nc{'Depth'}(:); %Mexnc way
    %    loaddap('+v',[url,'?Depth']); %fsu url var name
    %    clm.z=Depth; clear Depth;
    clm.z=lev; clear lev;
end

%determine dimensions of hycom data needed for interpolation
disp('optimizing grid dimensions ...');
[ym xm]=size(xg);
ii=find(xg>=xl & xg<=xr & yg>=yb & yg<=yt);jj=fix((ii-1)/ym)+1;ii=mod(ii,ym);
ig0=(min(ii)-1); ig1=(max(ii)+1); jg0=(min(jj)-1); jg1=(max(jj)+1);




if lever==0;

    %********************************
    irg=['[' num2str(ig0-1) ':' num2str(ig1-1) ']'];
    jrg=['[' num2str(jg0-1) ':' num2str(jg1-1) ']'];


    irg2=[num2str(ig0) ':' num2str(ig1)];
    jrg2=[num2str(jg0) ':' num2str(jg1)];
    %*********************************
clm.lon=xg(ig0:ig1,jg0:jg1); clm.lat=yg(ig0:ig1,jg0:jg1);

end

clear xg yg ii jj;

% catch

% end


%determine indices for time period of interpolation
disp('getting the number of time records ...');
if lever==0;
    % get hycom times
    %loaddap('+v',[url,'?time']); %loaddap way
    % time=nc{'time'}(:); %Mexnc way
    time=double(nc{'time'}(:));
    %time=nj_varget(url,'time'); %NJ Toolbox way
    tg=time+datenum(0000,12,30,0,0,0);
    tg2=julian(str2num(datestr(tg,'yyyy')),str2num(datestr(tg,'mm')),str2num(datestr(tg,'dd')),str2num(datestr(tg,'HH')))-2400001; %% 1 element.
end
if lever==1;
    %get hycom time fsu
    %these are used to convert time stamp to datenum
    t0=datenum(1900,12,31); tr0=datenum(1858,11,17);
    time=double(nc{'MT'}(:));
    MT=time;
    %time=nj_varget(url,'MT');
    %loaddap('+v',[url '?MT']);
    tg=MT+t0;clear MT;
    % time=nc{'MT'}(:); %Mexnc way
        tg2=julian(str2num(datestr(tg,'yyyy')),str2num(datestr(tg,'mm')),str2num(datestr(tg,'dd')),str2num(datestr(tg,'HH')))-2400001; %% 1 element.
end

%NJToolbox dims [start] [count] [stride]
starts=['[1 ',num2str(ig0),' ',num2str(jg0),']'];
counts=['[',num2str(length(time)),' ',num2str(ig1-ig0+1),' ',num2str(jg1-jg0+1),']'];

starts3d=['[1 1 ',num2str(ig0),' ',num2str(jg0),']'];
counts3d=['[',num2str(length(time)),' inf ',num2str(ig1-ig0+1),' ',num2str(jg1-jg0+1),']'];

clear time;

% get user times
%tid1=T1-tg(1)+1; tid2=T2-tg(1)+1;
[junk,tid1,ib]=intersect(tg,T1);
if tid1
else
    tid1=length(tg);
end
[junk,tid2,ib]=intersect(tg,T2);
if tid2
else
    tid2=length(tg);
end
disp([num2str(tid2-tid1+1) ' available record(s): ' num2str([tid1 tid2])]);
disp(['interpolating from ' datestr(T1) ' to ' datestr(T2)]);

%Create netcdf clm file
tag=datestr(tg(tid1),'yyyymmdd');
disp([datestr(tg(tid1)) ' at ' datestr(now)]);
%fn=[wdr 'USE_' tag '_clm.nc'];
fn=[wdr 'USE_coawst_clm.nc'];
disp(['creating netcdf file ',fn]);
create_roms_netcdf_clm_2_3(fn,gn,tid1,tid2);

%fill grid dims
RN=netcdf(fn,'w');
RN{'lon_rho'}(:) =gn.lon_rho;
RN{'lat_rho'}(:) =gn.lat_rho;
close(RN)

%%
disp('interpolating u,v,ubar,vbar ...');
for t=tid1:tid2
    tid=t-tid1+1;
    trg=num2str(t-1); disp([datestr(tg(t)) ' at ' datestr(now)]);
    clm.u=zeros([length(clm.z) size(gn.lon_rho)]);
    ttu=1;
    while ttu==1;
        try
            if lever==0;
                %eval(['tmp=double(nj_varget(''',url,''',''uogrddsl'',',starts3d,',',counts3d,'));']);
                eval(['tmp=double(nc{''uogrddsl''}(',num2str(tid1),',:,',irg2,',',jrg2,'));']);
            end
            if lever==1;
                for k=1:length(clm.z)
                    %eval(['tmp=double(nj_varget(''',url,''',''u'',',starts3d,',',counts3d,'));']);%fsu data
                    eval(['tmp=double(nc{''u''}(',num2str(tid1),',k,',irg2,',',jrg2,'));']);
                    %                    eval(['tmp=double(nc{''u''}(',num2str(tid1),',:,',num2str(ig0),':',num2str(ig1),',',num2str(jg0),':',num2str(jg1),'));']);

                    disp(['doing griddata u for ' num2str(k) ' at ' datestr(now)]);
                    clm.u(k,:,:)=griddata(clm.lon,clm.lat,squeeze(tmp(:,:)),gn.lon_rho,gn.lat_rho);%,'spline');
                    clear tmp;
                    clm.u(k,:,:)=maplev(squeeze(clm.u(k,:,:)));
                end
            end
            ttu=0;
        catch
            disp(['catch u Unable to download HYCOM u data at' datestr(now)]);
            fid=fopen('coawstlog.txt','a');
            fprintf(fid,'Unable to download HYCOM u data at');
            fprintf(fid,datestr(now));
            fprintf(fid,'\n');
        end
    end
    %== Vertical interpolation (t,s,u,v) from standard z-level to s-level
    u=roms_from_stdlev(gn.lon_rho,gn.lat_rho,clm.z,clm.u,gn,'u',0);
    clm=rmfield(clm,'u');
    save u.mat
    clear u;

    clm.v=zeros([length(clm.z) size(gn.lon_rho)]);
    ttv=1;
    while ttv==1;
        try
            if lever==0;
                %eval(['tmp=double(nj_varget(''',url,''',''vogrddsl'',',starts3d,',',counts3d,'));']);
                eval(['tmp=double(nc{''vogrddsl''}(',num2str(tid1),',:,',irg2,',',jrg2,'));']);
            end
            if lever==1;
                for k=1:length(clm.z)
                    %eval(['tmp=double(nj_varget(''',url,''',''v'',',starts3d,',',counts3d,'));']);%fsu data
                    eval(['tmp=double(nc{''v''}(',num2str(tid1),',k,',irg2,',',jrg2,'));']);
                    %                        eval(['tmp=double(nc{''v''}(',num2str(tid1),',:,',num2str(ig0),':',num2str(ig1),',',num2str(jg0),':',num2str(jg1),'));']);

                    disp(['doing griddata v for ' num2str(k) ' at ' datestr(now)]);
                    clm.v(k,:,:)=griddata(clm.lon,clm.lat,squeeze(tmp(:,:)),gn.lon_rho,gn.lat_rho);%,'spline');
                    clear tmp;
                    clm.v(k,:,:)=maplev(squeeze(clm.v(k,:,:)));
                end
            end
            ttv=0;
        catch
            disp('catch v');
            fid=fopen('coawstlog.txt','a');
            fprintf(fid,'Unable to download HYCOM v data at');
            fprintf(fid,datestr(now));
            fprintf(fid,'\n');
        end
    end
    %== Vertical interpolation (t,s,u,v) from standard z-level to s-level
    v=roms_from_stdlev(gn.lon_rho,gn.lat_rho,clm.z,clm.v,gn,'v',0);
    clm=rmfield(clm,'v');
    save v.mat
    clear v;

    %== Rotate the velocity
    theta=exp(-sqrt(-1)*mean(mean(gn.angle)));
    RN=netcdf(fn,'w');
    load u.mat; load v.mat

    if lever==0
        load rtofs_phi_USeast.mat
        phi2(phi2==0)=nan;
        phi2=maplev(phi2);
        for k=1:16
            [u2,v2]=rtofs2geo(squeeze(u(k,2:end,:)),squeeze(v(k,:,2:end)),phi2(2:end,2:end));
            u(k,2:end,:)=u2;
            v(k,:,2:end)=v2;
        end
    end
    disp('doing rotation to grid for u and v');
    uv=(u2rho_3d(u)+sqrt(-1)*v2rho_3d(v)).*theta;
    u=rho2u_3d(real(uv)); v=rho2v_3d(imag(uv));

    clear uv 
    
  %% == output
    RN{'u'}(tid,:,:,:)=u;
    RN{'v'}(tid,:,:,:)=v;
    
    clear u; clear v;
    %    end
    RN{'ocean_time'}(tid)=tg2(t)*3600*24;
    RN{'zeta_time'}(tid)=tg2(t);
    RN{'v2d_time'}(tid)=tg2(t);
    RN{'v3d_time'}(tid)=tg2(t);
    RN{'salt_time'}(tid)=tg2(t);
    RN{'temp_time'}(tid)=tg2(t);
    close(RN);

    %== Depth averaging u, v to get Ubar
    load u.mat; load v.mat

    cc=roms_zint(u,gn);  ubar=rho2u_2d(u2rho_2d(cc)./gn.h);
    cc=roms_zint(v,gn);  vbar=rho2v_2d(v2rho_2d(cc)./gn.h);
    %== Rotate the velocity
%     uv=(u2rho_2d(ubar)+sqrt(-1)*v2rho_2d(vbar)).*theta;
%     ubar=rho2u_2d(real(uv)); vbar=rho2v_2d(imag(uv));
    clear u;
    clear v;
    RN=netcdf(fn,'w');
    RN{'ubar'}(tid,:,:)=ubar;
    RN{'vbar'}(tid,:,:)=vbar;
    close(RN);
    tt=0;
    clear ubar
    clear vbar
end
%% interpolate the data
if lever==1
    close(nc);
    nc=mDataset(url2);  
end

disp('interpolating zeta ...');
for t=tid1:tid2
    tid=t-tid1+1;
    trg=['[' num2str(t-1) ']'];
    disp([datestr(tg(t)) ' at ' datestr(now)]);
    if lever==0;
        %eval(['tmp=double(nj_varget(''',url,''',''sshgsfc'',',starts,',',counts,'));']);
        eval(['tmp=double(nc{''sshgsfc''}(',num2str(tid1),',',irg2,',',jrg2,'));']);
    end
    if lever==1;
        %eval(['tmp=double(nj_varget(''',url,''',''ssh'',',starts,',',counts,'));']);%fsu analysis
         eval(['tmp=double(nc{''ssh''}(',num2str(tid1),',',irg2,',',jrg2,'));']);
%        eval(['tmp=double(nc{''ssh''}(',num2str(tid1),',',num2str(ig0),':',num2str(ig1),',',num2str(jg0),':',num2str(jg1),'));']);
    end
    zeta=griddata(clm.lon,clm.lat,tmp,gn.lon_rho,gn.lat_rho);%,'spline');
    zeta=maplev(zeta);
    clear tmp
    %== output
    RN=netcdf(fn,'w');
    RN{'zeta'}(tid,:,:) =zeta;
    close(RN);clear zeta;
end

%%
disp('interpolating temp ...');
if lever==1
    close(nc);
    nc=mDataset(url2);  
end
for t=tid1:tid2
    tid=t-tid1+1;
    trg=num2str(t-1); disp([datestr(tg(t)) ' at ' datestr(now)]);
    clm.temp=zeros([length(clm.z) size(gn.lon_rho)]);
    if lever==0;
        %eval(['tmp=double(nj_varget(''',url,''',''wtmpcdsl'',',starts3d,',',counts3d,'));']);
        eval(['tmp=double(nc{''wtmpcdsl''}(',num2str(tid1),',:,',irg2,',',jrg2,'));']);
%        eval(['tmp=double(nc{''wtmpdsl''}(',num2str(tid1),',:,',irg2,',',jrg2,'));']);
    end
    if lever==1;
        for k=1:length(clm.z)
            %eval(['tmp=double(nj_varget(''',url,''',''temperature'',',starts3d,',',counts3d,'));']);%fsu analysis
                    eval(['tmp=double(nc{''temperature''}(',num2str(tid1),',k,',irg2,',',jrg2,'));']);
            %eval(['tmp=double(nc{''temperature''}(',num2str(tid1),',k,[',num2str(ig0),':',num2str(ig1),'],[',num2str(jg0),':',num2str(jg1),']));']);

            clm.temp(clm.temp<0)=nan;
            clm.temp(k,:,:)=griddata(clm.lon,clm.lat,tmp(:,:),gn.lon_rho,gn.lat_rho);%,'spline');
            clear tmp
            clm.temp(k,:,:)=maplev(squeeze(clm.temp(k,:,:)));
        end
end    %== Vertical interpolation (t,s,u,v) from standard z-level to s-level
    temp=roms_from_stdlev(gn.lon_rho,gn.lat_rho,clm.z,clm.temp,gn,'rho',0);
    clm=rmfield(clm,'temp');
    %== output
    RN=netcdf(fn,'w');
    RN{'temp'}(tid,:,:,:)=temp;
    close(RN);clear temp;
end

%%
disp('interpolating salt ...');
if lever==1
    close(nc);
    nc=mDataset(url2);  
end
for t=tid1:tid2
    tid=t-tid1+1;
    trg=['[' num2str(t-1) ']']; disp([datestr(tg(t)) ' at ' datestr(now)]);
    clm.salt=zeros([length(clm.z) size(gn.lon_rho)]);
    if lever==0;
        %eval(['tmp=double(nj_varget(''',url,''',''salindsl'',',starts3d,',',counts3d,'));']);
        eval(['tmp=double(nc{''salindsl''}(',num2str(tid1),',:,',irg2,',',jrg2,'));']);
    end
    if lever==1;
        for k=1:length(clm.z)
            %eval(['tmp=double(nj_varget(''',url,''',''salinity'',',starts3d,',',counts3d,'));']);%fsu analysis
            eval(['tmp=double(nc{''salinity''}(',num2str(tid1),',k,',irg2,',',jrg2,'));']);
            %         eval(['tmp=double(nc{''salinity''}(',num2str(tid1),',:,',num2str(ig0),':',num2str(ig1),',',num2str(jg0),':',num2str(jg1),'));']);

            clm.salt(k,:,:)=griddata(clm.lon,clm.lat,(tmp(:,:)),gn.lon_rho,gn.lat_rho);%,'spline');
            clear tmp;
            clm.salt(k,:,:)=maplev(squeeze(clm.salt(k,:,:)));
        end
    end
    %== Vertical interpolation (t,s,u,v) from standard z-level to s-level
    salt=roms_from_stdlev(gn.lon_rho,gn.lat_rho,clm.z,clm.salt,gn,'rho',0);
    clm=rmfield(clm,'salt');
    %== output
    RN=netcdf(fn,'w');
    RN{'salt'}(tid,:,:,:)=salt;
    close(RN);clear salt;
end

%%
close(nc)
ncclose