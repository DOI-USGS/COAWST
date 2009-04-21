% written by Mingkui Li, May 2008
% Modified by Brandy Armstrong March 2009
% Modified again by jcw, April 2009

% Url of the global climate data
% url='http://hycom.coaps.fsu.edu/thredds/dodsC/glb_analysis';
url='http://hycom.coaps.fsu.edu:8080/thredds/dodsC/glb_analysis';

%these are used to convert time stamp to datenum
t0=datenum(1900,12,31); tr0=datenum(1858,11,17);

% get grid lat and lon
%disp('getting roms grid dimensions ...');
%gn=roms_get_grid(gridname,[theta_s theta_b Tcline N]);
xl=min(min(gn.lon_rho));xr=max(max(gn.lon_rho));
yb=min(min(gn.lat_rho));yt=max(max(gn.lat_rho));

%download hycom dimensions
disp('get hycom dimensions ...');
loaddap('+v',[url '?Longitude']);
xg=Longitude.Longitude-360;clear Longitude;
loaddap('+v',[url '?Latitude']);
yg=Latitude.Latitude;clear Latitude;
loaddap('+v',[url '?Depth']);
clm.z=Depth; clear Depth;

%determine dimensions of hycom data needed for interpolation
disp('optimizing grid dimensions ...');
[ym xm]=size(xg);
ii=find(xg>=xl & xg<=xr & yg>=yb & yg<=yt);jj=fix((ii-1)/ym)+1;ii=mod(ii,ym);
ig0=(min(ii)-1); ig1=(max(ii)+1); jg0=(min(jj)-1); jg1=(max(jj)+1);
irg=['[' num2str(ig0) ':' num2str(ig1) ']'];
jrg=['[' num2str(jg0) ':' num2str(jg1) ']'];
clm.lon=xg(ig0:ig1,jg0:jg1); clm.lat=yg(ig0:ig1,jg0:jg1);
clear xg yg ii jj;

%determine indices for time period of interpolation
disp('getting the number of time records ...');
% get hycom times
loaddap('+v',[url '?MT']);
tg=MT+t0;clear MT;
% get user times
%tid1=T1-tg(1)+1; tid2=T2-tg(1)+1;
[junk,tid1,ib]=intersect(tg,T1);
[junk,tid2,ib]=intersect(tg,T2);
disp([num2str(tid2-tid1+1) ' available record(s): ' num2str([tid1 tid2])]);
disp(['interpolating from ' datestr(T1) ' to ' datestr(T2)]);

%Create netcdf clm file
tag=datestr(tg(tid1),'yyyymm');
disp([datestr(tg(tid1)) ' at ' datestr(now)]);
fn=[wdr 'ncoda_' tag '.nc'];
disp(['creating netcdf file ',fn]);
create_roms_netcdf_clm(fn,gn,tid1,tid2);

%fill grid dims
RN=netcdf(fn,'w');
RN{'lon_rho'}(:) =gn.lon_rho;
RN{'lat_rho'}(:) =gn.lat_rho;
close(RN)

%% interpolate the data     
disp('interpolating zeta ...');
for t=tid1:tid2
    tid=t-tid1+1;
%    tag=datestr(tg(t),'yyyymm'); nr0=datenum([tag '01'],'yyyymmdd');
    trg=['[' num2str(t-1) ']'];
    disp([datestr(tg(t)) ' at ' datestr(now)]);
    tmp=loaddap('+v',[url '?ssh' trg irg jrg ]);
    tmp=tmp.ssh.ssh;  ii=find(tmp>5 | tmp<-5); tmp(ii)=nan; tmp=maplev(tmp);
    %roms.zeta=griddata(clm.lon,clm.lat,tmp,roms.lon,roms.lat);%,'spline');
    zeta=griddata(clm.lon,clm.lat,tmp,gn.lon_rho,gn.lat_rho);%,'spline');
    clear tmp
    %== output
    RN=netcdf(fn,'w');
    RN{'zeta'}(tid,:,:) =zeta;
    close(RN);clear zeta;
end

disp('interpolating temp ...');
for t=tid1:tid2
    tid=t-tid1+1;
    trg=['[' num2str(t-1) ']']; disp([datestr(tg(t)) ' at ' datestr(now)]);
    clm.temp=zeros([length(clm.z) size(gn.lon_rho)]);
    for k=1:length(clm.z)
        
        krg=['[' num2str(k-1) ']'];
        tmp=loaddap('+v',[url '?temperature' trg krg irg jrg ]);
        tmp=tmp.temperature.temperature; tmp(abs(tmp)>50)=NaN; tmp=maplev(tmp);
        %        clm.temp(k,:,:)=griddata(clm.lon,clm.lat,tmp,roms.lon,roms.lat);%,'spline');
        clm.temp(k,:,:)=griddata(clm.lon,clm.lat,tmp,gn.lon_rho,gn.lat_rho);%,'spline');
        clear tmp
    end
    %== Vertical interpolation (t,s,u,v) from standard z-level to s-level
    temp=roms_from_stdlev(gn.lon_rho,gn.lat_rho,clm.z,clm.temp,gn,'rho',0);
    clm=rmfield(clm,'temp');
    %== output
    RN=netcdf(fn,'w');
    RN{'temp'}(tid,:,:,:)=temp;
    close(RN);clear temp;
end

disp('interpolating salt ...');
for t=tid1:tid2
    tid=t-tid1+1;
    trg=['[' num2str(t-1) ']']; disp([datestr(tg(t)) ' at ' datestr(now)]);
    clm.salt=zeros([length(clm.z) size(gn.lon_rho)]);
    for k=1:length(clm.z)
        krg=['[' num2str(k-1) ']'];
        tmp=loaddap('+v',[url '?salinity' trg krg irg jrg ]);
        tmp=tmp.salinity.salinity; tmp(abs(tmp)>50)=NaN; tmp=maplev(tmp);
        %clm.salt(k,:,:)=griddata(clm.lon,clm.lat,tmp,roms.lon,roms.lat);%,'spline');
        clm.salt(k,:,:)=griddata(clm.lon,clm.lat,tmp,gn.lon_rho,gn.lat_rho);%,'spline');
        clear tmp
    end
    %== Vertical interpolation (t,s,u,v) from standard z-level to s-level
    salt=roms_from_stdlev(gn.lon_rho,gn.lat_rho,clm.z,clm.salt,gn,'rho',0);
    clm=rmfield(clm,'salt');
    %== output
    RN=netcdf(fn,'w');
    RN{'salt'}(tid,:,:,:)=salt;
    close(RN);clear salt;
end

disp('interpolating u,v,ubar,vbar ...');
for t=tid1:tid2
    tid=t-tid1+1;
    trg=['[' num2str(t-1) ']']; disp([datestr(tg(t)) ' at ' datestr(now)]);

    clm.u=zeros([length(clm.z) size(gn.lon_rho)]);
    for k=1:length(clm.z)
        krg=['[' num2str(k-1) ']'];
        tmp=loaddap('+v',[url '?u' trg krg irg jrg ]);
        tmp=tmp.u.u; tmp(abs(tmp)>50)=NaN; tmp=maplev(tmp);
        %clm.u(k,:,:)=griddata(clm.lon,clm.lat,tmp,roms.lon,roms.lat);%,'spline');
        clm.u(k,:,:)=griddata(clm.lon,clm.lat,tmp,gn.lon_rho,gn.lat_rho);%,'spline');
        clear tmp
    end
    %== Vertical interpolation (t,s,u,v) from standard z-level to s-level
    u=roms_from_stdlev(gn.lon_rho,gn.lat_rho,clm.z,clm.u,gn,'u',0);
    clm=rmfield(clm,'u');
%    RN=netcdf(fn,'w');
%    RN{'u'}(tid,:,:,:)=u;  %still needs to be rotated
%    close(RN);
     save u.mat
     clear u;

    clm.v=zeros([length(clm.z) size(gn.lon_rho)]);
    for k=1:length(clm.z)
        krg=['[' num2str(k-1) ']'];
        tmp=loaddap('+v',[url '?v' trg krg irg jrg ]);
        tmp=tmp.v.v; tmp(abs(tmp)>50)=NaN; tmp=maplev(tmp);
        %clm.v(k,:,:)=griddata(clm.lon,clm.lat,tmp,roms.lon,roms.lat);%,'spline');
        clm.v(k,:,:)=griddata(clm.lon,clm.lat,tmp,gn.lon_rho,gn.lat_rho);%,'spline');
        clear tmp
    end
    %== Vertical interpolation (t,s,u,v) from standard z-level to s-level
    v=roms_from_stdlev(gn.lon_rho,gn.lat_rho,clm.z,clm.v,gn,'v',0);
    clm=rmfield(clm,'v');
%    RN=netcdf(fn,'w');
%    RN{'v'}(tid,:,:,:)=v;  %still needs to be rotated
%    close(RN);
     save v.mat
     clear v;
    
    %== Rotate the velocity
    theta=exp(-sqrt(-1)*mean(mean(gn.angle)));
    RN=netcdf(fn,'w');
%    for k=1:length(clm.z)
%      u=squeeze(RN{'u'}(tid,k,:,:));
%      v=squeeze(RN{'v'}(tid,k,:,:));
      load u.mat; load v.mat
      uv=(u2rho_3d(u)+sqrt(-1)*v2rho_3d(v)).*theta;
      u=rho2u_3d(real(uv)); v=rho2v_3d(imag(uv));
      clear uv
      %== output
      RN{'u'}(tid,:,:,:)=u;
      RN{'v'}(tid,:,:,:)=v;
%    end
    RN{'ocean_time'}(tid)=tg(t)-tr0;
    close(RN);

    %== Depth averaging u, v to get Ubar
%   load u.mat; load v.mat
    cc=roms_zint(u,gn);  ubar=rho2u_2d(u2rho_2d(cc)./gn.h);
    cc=roms_zint(v,gn);  vbar=rho2v_2d(v2rho_2d(cc)./gn.h);
    %== Rotate the velocity
    uv=(u2rho_2d(ubar)+sqrt(-1)*v2rho_2d(vbar)).*theta;
    ubar=rho2u_2d(real(uv)); vbar=rho2v_2d(imag(uv));
    clear uv
    RN=netcdf(fn,'w');
    RN{'ubar'}(tid,:,:)=ubar;
    RN{'vbar'}(tid,:,:)=vbar;
    close(RN);
    
end

toc

%% Call to create and fill boundary file
%udatbdry(tg,fn);


