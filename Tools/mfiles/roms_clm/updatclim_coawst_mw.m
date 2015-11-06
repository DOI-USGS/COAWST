function [fn]=updatclim_coawst_mw(T1,gn,clmname,wdr,url2)
% Modified by Brandy Armstrong January 2012 to use only NCTOOLBOX 
% and Matlab builtin functions to read and write netcdf files

%T1 = date for climatology file
%gn = data from grid
%wdr = the working directory
%gname_pre = grid name prefix for climatology filenames

Time1=datestr(T1,'yyyymmdd');
Time_before=datestr(T1-1,'yyyymmdd');

%%
lever=1;
%url for hycom data
%url2='http://tds.hycom.org/thredds/dodsC/GLBa0.08/expt_90.9'; %/2011';
nc=ncgeodataset(url2);
load hycom_info.mat
%% 

%determine indices for time period of interpolation
disp('getting the number of time records ...');
if lever==1;
    %get hycom time fsu
    t0=datenum(1900,12,31); tr0=datenum(1858,11,17);
    time=nc.data('MT');
    MT=time;
    tg=MT+t0;clear MT;
    tg2=julian(str2num(datestr(tg,'yyyy')),str2num(datestr(tg,'mm')),str2num(datestr(tg,'dd')),str2num(datestr(tg,'HH')))-2400001; %% 1 element.
end

clear time;

% get user times
%[junk,tid1,ib]=intersect(tg,T1);
[junk,tid1,ib]=intersect(tg,floor(T1)); %modify to be nearest jcw 23Aug2014
if tid1
else
    tid1=length(tg);
end

disp(['interpolating for ' datestr(T1)]);

%Create netcdf clm file
tag=datestr(tg(tid1),'yyyymmdd');
disp([datestr(tg(tid1)) ' at ' datestr(now)]);
fn=[clmname];
disp(['creating netcdf file ',fn]);
create_roms_netcdf_clm_mwUL(fn,gn,1);% converted to BI functions

%fill grid dims using builtin (BI) functions
RN=netcdf.open(fn,'NC_WRITE');
lonid=netcdf.inqVarID(RN,'lon_rho');
netcdf.putVar(RN,lonid,gn.lon_rho);
latid=netcdf.inqVarID(RN,'lat_rho');
netcdf.putVar(RN,latid,gn.lat_rho);
netcdf.close(RN)

%%
disp('interpolating u,v,ubar,vbar ...');

disp([datestr(tg(tid1)) ' at ' datestr(now)]);
tz_levs=length(clm.z);
ttu=1;
while ttu==1;
    try
        if lever==0;
            clm.u=zeros([(tz_levs/2)+1 size(gn.lon_rho)]);
            level=1;
                tmpt=nc.geovariable('uogrddsl');
            for k=1:2:tz_levs%had to many levels, ran out of MEMORY
                eval(['tmp=tmpt.data(',num2str(tid1),',',num2str(k),',',irg2,',',jrg2,');']);
                disp(['doing griddata u for level ' num2str(k)]);
                clm.u(level,:,:)=griddata(clm.lon,clm.lat,double(squeeze(tmp)),gn.lon_rho,gn.lat_rho);%,'spline');
                clear tmp;
                clm.u(level,:,:)=maplev(squeeze(clm.u(level,:,:)));
                level=level+1;
            end
            %make sure to get the last level so it has top and bottom
            eval(['tmp=tmpt.data(',num2str(tid1),',',num2str(length(clm.z)),',',irg2,',',jrg2,');']);
            disp(['doing griddata u for level ' num2str(tz_levs)]);
            clm.u(level,:,:)=griddata(clm.lon,clm.lat,double(squeeze(tmp)),gn.lon_rho,gn.lat_rho);%,'spline');
            clm.u(level,:,:)=maplev(squeeze(clm.u(level,:,:)));
            clm.z=[clm.z(1:2:end);clm.z(end)];
        end
        if lever==1;
            clm.u=zeros([length(clm.z) size(gn.lon_rho)]);
                tmpt=nc.geovariable('u');
            for k=1:tz_levs
                eval(['tmp=tmpt.data(',num2str(tid1),',',num2str(k),',',irg2,',',jrg2,');']);
                disp(['doing griddata u for HYCOM level ' num2str(k)]);
                clm.u(k,:,:)=griddata(clm.lon,clm.lat,double(squeeze(tmp)),gn.lon_rho,gn.lat_rho);%,'spline');
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
save MW_test.mat
%== Vertical interpolation (t,s,u,v) from standard z-level to s-level
u=roms_from_stdlev_mw(gn.lon_rho,gn.lat_rho,clm.z,clm.u,gn,'u',0);
clm=rmfield(clm,'u');
save u.mat u
clear u;

ttv=1;
clm.v=zeros([length(clm.z) size(gn.lon_rho)]);
while ttv==1;
    try
        if lever==0;
            level=1;
                tmpt=nc.geovariable('vogrddsl');
            for k=1:2:tz_levs%had to many levels, ran out of MEMORY
                eval(['tmp=tmpt.data(',num2str(tid1),',',num2str(k),',',irg2,',',jrg2,');']);
                disp(['doing griddata v for level ' num2str(k)]);
                clm.v(level,:,:)=griddata(clm.lon,clm.lat,double(squeeze(tmp)),gn.lon_rho,gn.lat_rho);%,'spline');
                clear tmp;
                clm.v(level,:,:)=maplev(squeeze(clm.v(level,:,:)));
                level=level+1;
            end
            %make sure to get the last level so it has top and bottom
            eval(['tmp=tmpt.data(',num2str(tid1),',',num2str(length(clm.z)),',',irg2,',',jrg2,');']);
            disp(['doing griddata v for level ' num2str(tz_levs)]);
            clm.v(level,:,:)=griddata(clm.lon,clm.lat,double(squeeze(tmp)),gn.lon_rho,gn.lat_rho);%,'spline');
            clm.v(level,:,:)=maplev(squeeze(clm.v(level,:,:)));
        end
        if lever==1;
               tmpt=nc.geovariable('v');
            for k=1:tz_levs
                eval(['tmp=tmpt.data(',num2str(tid1),',',num2str(k),',',irg2,',',jrg2,');']);
                disp(['doing griddata v for HYCOM level ' num2str(k)]);
                clm.v(k,:,:)=griddata(clm.lon,clm.lat,double(squeeze(tmp)),gn.lon_rho,gn.lat_rho);%,'spline');
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
v=roms_from_stdlev_mw(gn.lon_rho,gn.lat_rho,clm.z,clm.v,gn,'v',0);
clm=rmfield(clm,'v');
save v.mat v
clear v;

%== Rotate the velocity
theta=exp(-sqrt(-1)*mean(mean(gn.angle)));

load u.mat; load v.mat

if lever==0
    load rtofs_phi_USeast_MW.mat
    phi2(phi2==0)=nan;
    phi2=maplev(phi2);
    for k=1:16
        [u2,v2]=rtofs2geo(squeeze(u(k,:,2:end)),squeeze(v(k,2:end,:)),phi2(2:end,2:end));
        u(k,:,2:end)=u2;
        v(k,2:end,:)=v2;
    end
end
disp('doing rotation to grid for u and v');
uv=(u2rho_3d_mw(u)+sqrt(-1)*v2rho_3d_mw(v)).*theta;
u=rho2u_3d_mw(real(uv)); v=rho2v_3d_mw(imag(uv));

clear uv

%% == output
RN=netcdf.open(fn,'NC_WRITE');

tempid=netcdf.inqVarID(RN,'u');
netcdf.putVar(RN,tempid,shiftdim(u,1));

tempid=netcdf.inqVarID(RN,'v');
netcdf.putVar(RN,tempid,shiftdim(v,1));

clear u; clear v;
tempid=netcdf.inqVarID(RN,'ocean_time');
netcdf.putVar(RN,tempid,tg2(tid1));
tempid=netcdf.inqVarID(RN,'zeta_time');
netcdf.putVar(RN,tempid,tg2(tid1));
tempid=netcdf.inqVarID(RN,'v2d_time');
netcdf.putVar(RN,tempid,tg2(tid1));
tempid=netcdf.inqVarID(RN,'v3d_time');
netcdf.putVar(RN,tempid,tg2(tid1));
tempid=netcdf.inqVarID(RN,'salt_time');
netcdf.putVar(RN,tempid,tg2(tid1));
tempid=netcdf.inqVarID(RN,'temp_time');
netcdf.putVar(RN,tempid,tg2(tid1));
netcdf.close(RN);
%%
%== Depth averaging u, v to get Ubar
load u.mat; load v.mat
cc=roms_zint_mw(u,gn);  ubar=rho2u_2d_mw(u2rho_2d_mw(cc)./gn.h);
cc=roms_zint_mw(v,gn);  vbar=rho2v_2d_mw(v2rho_2d_mw(cc)./gn.h);
%== Rotate the velocity
uv=(u2rho_2d_mw(ubar)+sqrt(-1)*v2rho_2d_mw(vbar)).*theta;
ubar=rho2u_2d_mw(real(uv)); vbar=rho2v_2d_mw(imag(uv));
clear u
clear v

RN=netcdf.open(fn,'NC_WRITE');
tempid=netcdf.inqVarID(RN,'ubar');
netcdf.putVar(RN,tempid,ubar);
tempid=netcdf.inqVarID(RN,'vbar');
netcdf.putVar(RN,tempid,vbar);
netcdf.close(RN);

clear ubar
clear vbar
clear uv

%% interpolate the zeta data
disp('interpolating zeta ...');
disp([datestr(tg(tid1)) ' at ' datestr(now)]);
if lever==0;
    tmpt=nc.geovariable('sshgsfc');
    eval(['tmp=tmpt.data(',num2str(tid1),',',irg2,',',jrg2,');']);
end
if lever==1;
    tmpt=nc.geovariable('ssh');
    eval(['tmp=tmpt.data(',num2str(tid1),',',irg2,',',jrg2,');']);
end
zeta=griddata(clm.lon,clm.lat,double(squeeze(tmp)),gn.lon_rho,gn.lat_rho);%,'spline');
zeta=maplev(zeta);
clear tmp


%== output
RN=netcdf.open(fn,'NC_WRITE');
tempid=netcdf.inqVarID(RN,'zeta');
netcdf.putVar(RN,tempid,zeta);
netcdf.close(RN);
clear zeta;


%%
disp('interpolating temp ...');

disp([datestr(tg(tid1)) ' at ' datestr(now)]);
clm.temp=zeros([length(clm.z) size(gn.lon_rho)]);
if lever==0;
    level=1;
        tmpt=nc.geovariable('wtmpcdsl');
    for k=1:2:tz_levs
        disp(['doing griddata temp for level ' num2str(k)]);
        eval(['tmp=tmpt.data(',num2str(tid1),',',num2str(k),',',irg2,',',jrg2,');']);
        clm.temp(clm.temp<0)=nan;
        clm.temp(level,:,:)=griddata(clm.lon,clm.lat,double(squeeze(tmp)),gn.lon_rho,gn.lat_rho);%,'spline');
        clear tmp
        clm.temp(level,:,:)=maplev(squeeze(clm.temp(level,:,:)));
        level=level+1;
    end
    tmpt=nc.geovariable('wtmpcdsl');
    eval(['tmp=tmpt.data(',num2str(tid1),',',num2str(tz_levs),',',irg2,',',jrg2,');']);
    clm.temp(clm.temp<0)=nan;
    clm.temp(level,:,:)=griddata(clm.lon,clm.lat,double(squeeze(tmp)),gn.lon_rho,gn.lat_rho);%,'spline');
    clear tmp
    clm.temp(level,:,:)=maplev(squeeze(clm.temp(level,:,:)));
end
if lever==1;
        tmpt=nc.geovariable('temperature');
    for k=1:tz_levs
        disp(['doing griddata temp for HYCOM level ' num2str(k)]);
        eval(['tmp=tmpt.data(',num2str(tid1),',',num2str(k),',',irg2,',',jrg2,');']);
        clm.temp(clm.temp<0)=nan;
        clm.temp(k,:,:)=griddata(clm.lon,clm.lat,double(squeeze(tmp)),gn.lon_rho,gn.lat_rho);%,'spline');
        clear tmp
        clm.temp(k,:,:)=maplev(squeeze(clm.temp(k,:,:)));
    end
end    

%== Vertical interpolation (t,s,u,v) from standard z-level to s-level
temp=roms_from_stdlev_mw(gn.lon_rho,gn.lat_rho,clm.z,clm.temp,gn,'rho',0);
clm=rmfield(clm,'temp');

%== output
RN=netcdf.open(fn,'NC_WRITE');
tempid=netcdf.inqVarID(RN,'temp');
netcdf.putVar(RN,tempid,shiftdim(temp,1));
netcdf.close(RN);
clear temp;


%%
disp('interpolating salt ...');

disp([datestr(tg(tid1)) ' at ' datestr(now)]);
clm.salt=zeros([length(clm.z) size(gn.lon_rho)]);
if lever==0;
    level=1;
        tmpt=nc.geovariable('salindsl');
    for k=1:2:tz_levs
        disp(['doing griddata salt for level ' num2str(k)]);
        eval(['tmp=tmpt.data(',num2str(tid1),',',num2str(k),',',irg2,',',jrg2,');']);
        clm.salt(level,:,:)=griddata(clm.lon,clm.lat,double(squeeze(tmp)),gn.lon_rho,gn.lat_rho);
        clear tmp;
        clm.salt(level,:,:)=maplev(squeeze(clm.salt(level,:,:)));
        level-level+1;
    end
    eval(['tmp=tmpt.data(',num2str(tid1),',',num2str(k),',',irg2,',',jrg2,');']);
    clm.salt(level,:,:)=griddata(clm.lon,clm.lat,double(squeeze(tmp)),gn.lon_rho,gn.lat_rho);
    clear tmp;
    clm.salt(level,:,:)=maplev(squeeze(clm.salt(level,:,:)));
end
if lever==1;
        tmpt=nc.geovariable('salinity');
    for k=1:tz_levs
        disp(['doing griddata salt for HYCOM level ' num2str(k)]);
        eval(['tmp=tmpt.data(',num2str(tid1),',',num2str(k),',',irg2,',',jrg2,');']);
        clm.salt(k,:,:)=griddata(clm.lon,clm.lat,double(squeeze(tmp)),gn.lon_rho,gn.lat_rho);
        clear tmp;
        clm.salt(k,:,:)=maplev(squeeze(clm.salt(k,:,:)));
    end
end
%== Vertical interpolation (t,s,u,v) from standard z-level to s-level
salt=roms_from_stdlev_mw(gn.lon_rho,gn.lat_rho,clm.z,clm.salt,gn,'rho',0);
clm=rmfield(clm,'salt');

%== output
RN=netcdf.open(fn,'NC_WRITE');
tempid=netcdf.inqVarID(RN,'salt');
netcdf.putVar(RN,tempid,shiftdim(salt,1));
netcdf.close(RN);
clear salt;
close(nc)

disp(['Finished creating clim file at ' datestr(now)]);
%%
