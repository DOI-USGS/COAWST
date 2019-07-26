function [fn]=updatclim_coawst_mw(T1, gn, clm, clmname, wdr, url)
% Modified by Brandy Armstrong January 2012 to use only NCTOOLBOX 
% and Matlab builtin functions to read and write netcdf files
% jcw Feb 2019 - only use matalb BI
%
%T1 = date for climatology file
%gn = data from grid
%clm = data of hycom indices
%wdr = the working directory
%clmname = grid name prefix for climatology filenames
%url = where get data from

%
%determine indices for time period of interpolation
%
disp('getting the number of time records ...');
t0=datenum(1900,12,31); % tr0=datenum(1858,11,17);
time=ncread(url,'MT');
tg=time+t0;
tg2=julian(str2num(datestr(tg,'yyyy')),str2num(datestr(tg,'mm')),str2num(datestr(tg,'dd')),str2num(datestr(tg,'HH')))-2400001;
%
% get user times
%
[junk,tid1,ib]=intersect(tg,floor(T1)); %modify to be nearest jcw 23Aug2014
if isempty(tid1)
  tid1=length(tg);
end

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
tz_levs=length(clm.z);
X=repmat(clm.lon,1,length(clm.lat));
Y=repmat(clm.lat,length(clm.lon),1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['Interpolating u for ',datestr(tg(tid1))]);
ttu=1;
clm.u=zeros([length(clm.z) size(gn.lon_rho)]);
while ttu==1;
    try
        tmpt=ncread(url,'u',[clm.ig0 clm.jg0 1 tid1],[clm.ig1-clm.ig0+1 clm.jg1-clm.jg0+1 tz_levs 1 ] );
        for k=1:tz_levs
            disp(['doing griddata u for HYCOM level ' num2str(k)]);
            tmp=double(squeeze(tmpt(:,:,k)));
            F = scatteredInterpolant(X(:),Y(:),tmp(:));
            cff = F(gn.lon_rho,gn.lat_rho);
            clm.u(k,:,:)=maplev(cff);
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
u=roms_from_stdlev_mw(gn.lon_rho,gn.lat_rho,clm.z,clm.u,gn,'u',0);
clm=rmfield(clm,'u');
save u.mat u
clear u;

disp(['Interpolating v for ',datestr(tg(tid1))]);
ttv=1;
clm.v=zeros([length(clm.z) size(gn.lon_rho)]);
while ttv==1;
    try
        tmpt=ncread(url,'v',[clm.ig0 clm.jg0 1 tid1],[clm.ig1-clm.ig0+1 clm.jg1-clm.jg0+1 tz_levs 1 ] );
        for k=1:tz_levs
            disp(['doing griddata v for HYCOM level ' num2str(k)]);
            tmp=double(squeeze(tmpt(:,:,k)));
            F = scatteredInterpolant(X(:),Y(:),tmp(:));
            cff = F(gn.lon_rho,gn.lat_rho);
            clm.v(k,:,:)=maplev(cff);
        end
        ttv=0;
    catch
        disp(['catch v Unable to download HYCOM v data at' datestr(now)]);
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% interpolate the zeta data
disp(['Interpolating zeta for ',datestr(tg(tid1))]);
ttz=1;
while ttz==1;
    try
        tmpt=ncread(url,'ssh',[clm.ig0 clm.jg0 tid1],[clm.ig1-clm.ig0+1 clm.jg1-clm.jg0+1 1 ] );
        tmp=double(squeeze(tmpt(:,:)));
        disp(['doing griddata zeta for HYCOM ']);
        F = scatteredInterpolant(X(:),Y(:),tmp(:));
        cff = F(gn.lon_rho,gn.lat_rho);
        zeta=maplev(cff);
        ttz=0;
    catch
        disp(['catch z Unable to download HYCOM ssh data at' datestr(now)]);
        fid=fopen('coawstlog.txt','a');
        fprintf(fid,'Unable to download HYCOM ssh data at');
        fprintf(fid,datestr(now));
        fprintf(fid,'\n');
    end
end
clear tmp
%
%== output zeta
%
RN=netcdf.open(fn,'NC_WRITE');
tempid=netcdf.inqVarID(RN,'zeta');
netcdf.putVar(RN,tempid,zeta);
netcdf.close(RN);
clear zeta;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['Interpolating temp for ',datestr(tg(tid1))]);
ttt=1;
clm.temp=zeros([length(clm.z) size(gn.lon_rho)]);
while ttt==1;
    try
        tmpt=ncread(url,'temperature',[clm.ig0 clm.jg0 1 tid1],[clm.ig1-clm.ig0+1 clm.jg1-clm.jg0+1 tz_levs 1 ] );
        for k=1:tz_levs
            disp(['doing griddata temp for HYCOM level ' num2str(k)]);
            tmp=double(squeeze(tmpt(:,:,k)));
            F = scatteredInterpolant(X(:),Y(:),tmp(:));
            cff = F(gn.lon_rho,gn.lat_rho);
%           cff(cff<0)=nan;
            clm.temp(k,:,:)=maplev(cff);
        end
        ttt=0;
    catch
        disp(['catch temp Unable to download HYCOM temp data at' datestr(now)]);
        fid=fopen('coawstlog.txt','a');
        fprintf(fid,'Unable to download HYCOM temp data at');
        fprintf(fid,datestr(now));
        fprintf(fid,'\n');
    end
end
%
%== Vertical interpolation (t,s,u,v) from standard z-level to s-level
%
temp=roms_from_stdlev_mw(gn.lon_rho,gn.lat_rho,clm.z,clm.temp,gn,'rho',0);
clm=rmfield(clm,'temp');
%
%== output temp
%
RN=netcdf.open(fn,'NC_WRITE');
tempid=netcdf.inqVarID(RN,'temp');
netcdf.putVar(RN,tempid,shiftdim(temp,1));
netcdf.close(RN);
clear temp;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['Interpolating salt for ',datestr(tg(tid1))]);
tts=1;
clm.salt=zeros([length(clm.z) size(gn.lon_rho)]);
while tts==1;
    try
        tmpt=ncread(url,'salinity',[clm.ig0 clm.jg0 1 tid1],[clm.ig1-clm.ig0+1 clm.jg1-clm.jg0+1 tz_levs 1 ] );
        for k=1:tz_levs
            disp(['doing griddata salt for HYCOM level ' num2str(k)]);
            tmp=double(squeeze(tmpt(:,:,k)));
            F = scatteredInterpolant(X(:),Y(:),tmp(:));
            cff = F(gn.lon_rho,gn.lat_rho);
            cff(cff<0)=nan;
            clm.salt(k,:,:)=maplev(cff);
        end
        tts=0;
    catch
        disp(['catch temp Unable to download HYCOM temp data at' datestr(now)]);
        fid=fopen('coawstlog.txt','a');
        fprintf(fid,'Unable to download HYCOM temp data at');
        fprintf(fid,datestr(now));
        fprintf(fid,'\n');
    end
end
%
%== Vertical interpolation (t,s,u,v) from standard z-level to s-level
%
salt=roms_from_stdlev_mw(gn.lon_rho,gn.lat_rho,clm.z,clm.salt,gn,'rho',0);
clm=rmfield(clm,'salt');
%
%== output salt
%
RN=netcdf.open(fn,'NC_WRITE');
tempid=netcdf.inqVarID(RN,'salt');
netcdf.putVar(RN,tempid,shiftdim(salt,1));
netcdf.close(RN);
clear salt;

disp(['Finished creating clim file at ' datestr(now)]);
%%
