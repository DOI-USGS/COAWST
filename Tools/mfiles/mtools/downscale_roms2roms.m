function downscale_from_L0_new(grdpath,L0,L1,time_ini,time_end)

% This routine creates boundary and initial condition files for ROMS:
%    coawst_bdy.nc ; coawst_ini.nc
%    on a user-defined grid for a user-defined date.
%
% This is currently set up to use L0 grid outputs to force L1 grids for the
% NOPP-NHCI project, based on efforts by:
% Maitane Olabarrieta 09/26/2021
% modified by John Warner 08/31/2022
% remodified by Maitane 08/31/2022

% Input paths and names of the ini, clm and bry files to generate
init_file=['L0_',L1.name,'_ini.nc'];
clm_file=['L0_',L1.name,'_clm.nc'];
bry_file=['L0_',L1.name,'_bry.nc'];

%% Define the constant Z levels to which interpolate L0 3-dimensional variables
ZCL=[-2,2,5,10,20,35,50,75,100,150,200,250,375,550,750,1000,1250,1500,2000,2500,3000,3500,4000,4500,5000];
ZCL= -ZCL;
ZCL_N = length(ZCL);

%% Extract the characteristics of the child grids
eval(['L1.gridname=''',L1.modelgrid,''';']);
disp('getting roms grid dimensions ...');

if (L1.Vtransform==1)
    h=ncread(L1.gridname,'h');
    hmin=min(h(:));
    hc=min(max(hmin,0),L1.Tcline);
elseif (L1.Vtransform==2)
    h=ncread(L1.gridname,'h');
    hmin=max(0.1,min(h(:)));
    hc=L1.Tcline;
end

L1.hc=hc;           %stretching width used in ROMS
gn=get_roms_grid(L1.gridname,L1);
[nxr,nyr]=size(gn.lon_rho);
[nxu,nyu]=size(gn.lon_u);
[nxv,nyv]=size(gn.lon_v);
L1.maskr=gn.mask_rho;

%% Initialize the variables of interest
time=[time_ini:1/24:time_end];
time4d=[time_ini:1:time_end];

L1_u=zeros(nxu,nyu,L1.N,length(time4d));
L1_v=zeros(nxv,nyv,L1.N,length(time4d));
L1_temp=zeros(nxr,nyr,L1.N,length(time4d));
L1_salt=zeros(nxr,nyr,L1.N,length(time4d));
L1_ubar_4d=zeros(nxu,nyu,length(time4d));
L1_vbar_4d=zeros(nxv,nyv,length(time4d));
L1_zeta_4d=zeros(nxr,nyr,length(time4d));
L1_ubar=zeros(nxu,nyu,length(time));
L1_vbar=zeros(nxv,nyv,length(time));
L1_zeta=zeros(nxr,nyr,length(time));

L1_lon_max=max(max(gn.lon_rho))+0.1;
L1_lon_min=min(min(gn.lon_rho))-0.1;
L1_lat_max=max(max(gn.lat_rho))+0.1;
L1_lat_min=min(min(gn.lat_rho))-0.1;

%% READ L0 GRID INFORMATION AND EXTRACT VARIABLES
lonr=ncread(L0.grid,'lon_rho');
latr=ncread(L0.grid,'lat_rho');
[IIr,JJr]=find((lonr<=L1_lon_max & lonr>=L1_lon_min) & (latr<=L1_lat_max & latr>=L1_lat_min));
L0_xinir=IIr(1);
L0_xendr=IIr(end);
L0_yinir=JJr(1);
L0_yendr=JJr(end);

L0_lonr=lonr(IIr(1):IIr(end),JJr(1):JJr(end));
L0_latr=latr(IIr(1):IIr(end),JJr(1):JJr(end));
[L0_nxr,L0_nyr]=size(L0_lonr);

L0_xiniu=IIr(1);
L0_xendu=IIr(end)-1;
L0_yiniu=JJr(1);
L0_yendu=JJr(end);
L0_xiniv=IIr(1);
L0_xendv=IIr(end);
L0_yiniv=JJr(1);
L0_yendv=JJr(end)-1;

clear lonr latr IIr JJr

L0_lonu=ncread(L0.grid,'lon_u',[L0_xinir,L0_yinir],[L0_nxr-1,L0_nyr]);
L0_latu=ncread(L0.grid,'lat_u',[L0_xinir,L0_yinir],[L0_nxr-1,L0_nyr]);
L0_lonv=ncread(L0.grid,'lon_v',[L0_xinir,L0_yinir],[L0_nxr,L0_nyr-1]);
L0_latv=ncread(L0.grid,'lat_v',[L0_xinir,L0_yinir],[L0_nxr,L0_nyr-1]);
[L0_nxu,L0_nyu]=size(L0_lonu);
[L0_nxv,L0_nyv]=size(L0_lonv);

L0_maskr=ncread(L0.grid,'mask_rho',[L0_xinir,L0_yinir],[L0_nxr,L0_nyr]);
L0_angler=ncread(L0.grid,'angle',[L0_xinir,L0_yinir],[L0_nxr,L0_nyr]);
L0_h=ncread(L0.grid,'h',[L0_xinir,L0_yinir],[L0_nxr,L0_nyr]);
L0_masku=ncread(L0.grid,'mask_u',[L0_xinir,L0_yinir],[L0_nxu,L0_nyu]);
L0_maskv=ncread(L0.grid,'mask_v',[L0_xinir,L0_yinir],[L0_nxv,L0_nyv]);
L0_hc=ncread(L0.out,'hc');
[L0_N,pp] = size(ncread(L0.out,'s_rho'));
clear pp

for zz=1:L0_N
    L0_lonr_z(:,:,zz)=L0_lonr;
    L0_latr_z(:,:,zz)=L0_latr;
end

ocean_time=double(ncread(L0.qck,'ocean_time'))./(3600*24)+datenum(1858,11,17,00,00,00);
dum=find(ocean_time<=time_end & ocean_time>=time_ini);
L0_nt2d=length(dum);
L0_tini2d=dum(1);
L0_tend2d=dum(end);
L0_ocean_time2d=ocean_time(dum);

ocean_time=double(ncread(L0.out,'ocean_time'))./(3600*24)+datenum(1858,11,17,00,00,00);
dum=find(ocean_time<=time_end & ocean_time>=time_ini);
L0_nt3d=length(dum);
L0_tini3d=dum(1);
L0_tend3d=dum(end);
L0_ocean_time3d=ocean_time(dum);

for zz=1:L1.N
    lon_rho_romsz(:,:,zz)=gn.lon_rho;
    lat_rho_romsz(:,:,zz)=gn.lat_rho;
end
lon_rho_roms(:,:)=gn.lon_rho;
lat_rho_roms(:,:)=gn.lat_rho;

tic

for tt=1:L0_nt2d
    
    disp(tt)
    
    % Interpolate 3-dimensional variables
    
    % Free surface elevation
    L0_zeta=double(ncread(L0.qck,'zeta',[L0_xinir,L0_yinir,L0_tini2d+tt-1],[L0_nxr,L0_nyr,1]));
    aa(:,:)=squeeze(L0_zeta(:,:,1));
    aa(L0_maskr==0)=nan;
    aa(:,:)=maplev(aa);
    zz=griddata(L0_lonr,L0_latr,aa,gn.lon_rho,gn.lat_rho);
    L1_zeta(:,:,tt)=zz(:,:);
    clear aa zz
    display('zeta:',num2str(tt))
    
    %   Ubar
    L0_ubar=double(ncread(L0.qck,'ubar',[L0_xiniu,L0_yiniu,L0_tini2d+tt-1],[L0_nxu,L0_nyu,1]));
    au(:,:)=squeeze(L0_ubar(:,:,1));
    au(L0_masku==0)=nan;
    au(:,:)=maplev(au);
    aur=u2rho_2d_mw(au);
    display('ubar:',num2str(tt))
    
    %   Vbar
    L0_vbar=double(ncread(L0.qck,'vbar',[L0_xiniv,L0_yiniv,L0_tini2d+tt-1],[L0_nxv,L0_nyv,1]));
    av(:,:)=squeeze(L0_vbar(:,:,1));
    av(L0_maskv==0)=nan;
    av(:,:)=maplev(av);
    avr=v2rho_2d_mw(av);
    display('vbar:',num2str(tt))
    
    % Compute Northward and Eastward velocities, important!
    vel=aur + avr.*sqrt(-1);
    vel=vel .* exp ( sqrt(-1) * L0_angler);
    velu=real(vel);
    velv=imag(vel);
    
    velu1=griddata(L0_lonr,L0_latr,velu,gn.lon_rho,gn.lat_rho);
    velv1=griddata(L0_lonr,L0_latr,velv,gn.lon_rho,gn.lat_rho);
    
    % Rotate velocities to ROMS grid, important!
    ubar1(:,:)=velu1.*cos(gn.angle)+velv1.*sin(gn.angle);
    vbar1(:,:)=velv1.*cos(gn.angle)-velu1.*sin(gn.angle);
    L1_ubar(:,:,tt)=rho2u_2d_mw(ubar1);  % defined at u points
    L1_vbar(:,:,tt)=rho2v_2d_mw(vbar1);  % defined at v points
    
    clear vel au av aur avr velu velv velu1 velv1 ubar1 vbar1
    clear L0_ubar L0_vbar
    
end
toc

%save('barotropic')

tic
tt=0;
for tt1=1:8:L0_nt3d
    tt=tt+1;
    disp(tt)
    
    % Interpolate 3-dimensional variables
    
    % Free surface elevation
    L0_zeta=double(ncread(L0.qck,'zeta',[L0_xinir,L0_yinir,L0_tini2d+(tt-1)*24],[L0_nxr,L0_nyr,1]));
    aa(:,:)=squeeze(L0_zeta(:,:,1));
    aa(L0_maskr==0)=nan;
    aa(:,:)=maplev(aa);
    zz=griddata(L0_lonr,L0_latr,aa,gn.lon_rho,gn.lat_rho);
    L1_zeta_4d(:,:,tt)=zz(:,:);
    clear aa zz
    
    L0_u=double(ncread(L0.out,'u',[L0_xiniu,L0_yiniu,1,L0_tini3d+tt1-1],[L0_nxu,L0_nyu,inf,1]));
    L0_v=double(ncread(L0.out,'v',[L0_xiniv,L0_yiniv,1,L0_tini3d+tt1-1],[L0_nxv,L0_nyv,inf,1]));
    
    L0_z(:,:,:)=set_depth(L0.Vtransform,L0.Vstretching,L0.theta_s,L0.theta_b,L0_hc,L0_N, ...
        5,L0_h,squeeze(L0_zeta(:,:,1)));
    L0_z(isnan(L0_z)==1)=0;

    for zz=1:L0_N
        L0_Hz(:,:,zz)=abs(L0_z(:,:,zz+1)-L0_z(:,:,zz));
    end
    
    [L0_ubar(:,:),L0_vbar(:,:)]=uv_barotropic(squeeze(L0_u(:,:,:,1)),squeeze(L0_v(:,:,:,1)),squeeze(L0_Hz(:,:,:,1)));    
    clear L0_Hz
    
    au(:,:)=L0_ubar(:,:);
    au(L0_masku==0)=nan;
    au(:,:)=maplev(au);
    aur=u2rho_2d_mw(au);
    
    av(:,:)=L0_vbar(:,:);
    av(L0_maskv==0)=nan;
    av(:,:)=maplev(av);
    avr=v2rho_2d_mw(av);
    
    % Compute Northward and Eastward velocities, important!
    vel=aur + avr.*sqrt(-1);
    vel=vel .* exp ( sqrt(-1) * L0_angler);
    velu=real(vel);
    velv=imag(vel);
    velu1=griddata(L0_lonr,L0_latr,velu,gn.lon_rho,gn.lat_rho);
    velv1=griddata(L0_lonr,L0_latr,velv,gn.lon_rho,gn.lat_rho);
    
    % Rotate velocities to ROMS grid, important!
    ubar1(:,:)=velu1.*cos(gn.angle)+velv1.*sin(gn.angle);
    vbar1(:,:)=velv1.*cos(gn.angle)-velu1.*sin(gn.angle);
    L1_ubar_4d(:,:,tt)=rho2u_2d_mw(ubar1);  % defined at u points
    L1_vbar_4d(:,:,tt)=rho2v_2d_mw(vbar1);  % defined at v points
    
    clear vel au av aur avr velu velv velu1 velv1 ubar1 vbar1
    clear L0_ubar L0_vbar
    
    % Interpolate 4-dimensional variables
    L0_temp=double(ncread(L0.out,'temp',[L0_xinir,L0_yinir,1,L0_tini3d+tt1-1],[L0_nxr,L0_nyr,inf,1]));
    L0_salt=double(ncread(L0.out,'salt',[L0_xinir,L0_yinir,1,L0_tini3d+tt1-1],[L0_nxr,L0_nyr,inf,1]));
    
    temp1(:,:,:)=L0_temp(:,:,:);
    salt1(:,:,:)=L0_salt(:,:,:);
    u1(:,:,:)=L0_u(:,:,:);
    v1(:,:,:)=L0_v(:,:,:);
    clear L0_temp L0_salt L0_u L0_v
    
    aa(:,:)=squeeze(L0_zeta(:,:,1));
    aa(L0_maskr==0)=nan;
    aa(:,:)=maplev(aa);
    L0_zr(:,:,:)=set_depth(L0.Vtransform,L0.Vstretching,L0.theta_s,L0.theta_b,L0_hc,L0_N, ...
        1,L0_h,aa);
    clear aa
    
    % Interpolate to constant Z vertical levels (ZCL)
    [nxr_L0,nyr_L0,nz_L0]=size(temp1);
    [nxu_L0,nyu_L0,nz_L0]=size(u1);
    [nxv_L0,nyv_L0,nz_L0]=size(v1);
    temp2=nan(nxr_L0,nyr_L0,ZCL_N);
    salt2=nan(nxr_L0,nyr_L0,ZCL_N);
    
    for ii=1:nxr_L0
        for jj=1:nyr_L0
            if (L0_maskr(ii,jj)==1)
                temp2(ii,jj,:)=interp1(squeeze(L0_zr(ii,jj,:)),squeeze(temp1(ii,jj,:)),ZCL,'linear');
                salt2(ii,jj,:)=interp1(squeeze(L0_zr(ii,jj,:)),squeeze(salt1(ii,jj,:)),ZCL,'linear');
            end
        end
    end
    temp2(:,:,1)=temp2(:,:,3);
    temp2(:,:,2)=temp2(:,:,3);
    salt2(:,:,1)=salt2(:,:,3);
    salt2(:,:,2)=salt2(:,:,3);
    clear temp1 salt1
    
    for ii=1:nxu_L0
        for jj=1:nyu_L0
            u2(ii,jj,:)=interp1(squeeze(L0_zr(ii,jj,:)),squeeze(u1(ii,jj,:)),ZCL,'linear');
        end
    end
    u2(:,:,1)=u2(:,:,3);
    u2(:,:,2)=u2(:,:,3);
    clear u1
    
    for ii=1:nxv_L0
        for jj=1:nyv_L0
            v2(ii,jj,:)=interp1(squeeze(L0_zr(ii,jj,:)),squeeze(v1(ii,jj,:)),ZCL,'linear');
        end
    end
    v2(:,:,1)=v2(:,:,3);
    v2(:,:,2)=v2(:,:,3);
    clear v1
    
    for zz=1:ZCL_N
        ur(:,:,zz)=u2rho_2d_mw(u2(:,:,zz));
        vr(:,:,zz)=v2rho_2d_mw(v2(:,:,zz));
    end
    clear u2 v2
    
    % Compute Northward and Eastward velocities
    for zz=1:ZCL_N
        vel(:,:)=ur(:,:,zz)+vr(:,:,zz).*sqrt(-1);
        vel=vel .* exp ( sqrt(-1) * L0_angler);
        ur(:,:,zz)=real(vel);
        vr(:,:,zz)=imag(vel);
    end
    clear vel
    
    %   Interpolate each constant ZCL layer to the horizontal L1 grid
    
    for zz=1:ZCL_N
        temp3(:,:,zz)= griddata(L0_lonr,L0_latr,squeeze(temp2(:,:,zz)),lon_rho_roms,lat_rho_roms);
        salt3(:,:,zz)= griddata(L0_lonr,L0_latr,squeeze(salt2(:,:,zz)),lon_rho_roms,lat_rho_roms);
        ur3(:,:,zz)= griddata(L0_lonr,L0_latr,squeeze(ur(:,:,zz)),lon_rho_roms,lat_rho_roms);
        vr3(:,:,zz)= griddata(L0_lonr,L0_latr,squeeze(vr(:,:,zz)),lon_rho_roms,lat_rho_roms);
    end
    clear ur vr temp2 salt2
    
    % Rotate velocities to ROMS grid, important!
    for zz=1:ZCL_N
        ur(:,:,zz)=ur3(:,:,zz).*cos(gn.angle)+vr3(:,:,zz).*sin(gn.angle);
        vr(:,:,zz)=vr3(:,:,zz).*cos(gn.angle)-ur3(:,:,zz).*sin(gn.angle);
        % Pass velocities from rho to u and v points
        u(:,:,zz)=rho2u_2d_mw(ur(:,:,zz));  % defined at u points
        v(:,:,zz)=rho2v_2d_mw(vr(:,:,zz));  % defined at v points
    end
    clear ur vr
    
    %   Interpolate to the L1 S levels
    %   gn.z_r assumed zeta =0 for L1. now we can use the actual L1 zeta to get the z_r's.
    L1_z(:,:,:)=set_depth(gn.Vtransform,gn.Vstretching,gn.theta_s,gn.theta_b,gn.hc,gn.N, ...
        1,gn.h,squeeze(L1_zeta_4d(:,:,tt)));
      
    for ii=1:nxr
        for jj=1:nyr
            if (gn.mask_rho(ii,jj)==1)
                L1_temp(ii,jj,:,tt)=interp1(ZCL,squeeze(temp3(ii,jj,:)),squeeze(L1_z(ii,jj,:)));
                L1_salt(ii,jj,:,tt)=interp1(ZCL,squeeze(salt3(ii,jj,:)),squeeze(L1_z(ii,jj,:)));
            else
                L1_temp(ii,jj,:,tt)=nan;
                L1_salt(ii,jj,:,tt)=nan;
            end
        end
    end
    clear temp3 salt3 L1_z
    
    L1_zu(:,:,:)=set_depth(gn.Vtransform,gn.Vstretching,gn.theta_s,gn.theta_b,gn.hc,gn.N, ...
        3,gn.h,squeeze(L1_zeta_4d(:,:,tt)));
    
    for ii=1:nxu
        for jj=1:nyu
            if (gn.mask_u(ii,jj)==1)
                L1_u(ii,jj,:,tt)=interp1(ZCL,squeeze(u(ii,jj,:)),squeeze(L1_zu(ii,jj,:)));
            else
                L1_u(ii,jj,:,tt)=nan;
            end
        end
    end
    clear u L1_zu
    
    L1_zv(:,:,:)=set_depth(gn.Vtransform,gn.Vstretching,gn.theta_s,gn.theta_b,gn.hc,gn.N, ...
        4,gn.h,squeeze(L1_zeta_4d(:,:,tt)));
    
    for ii=1:nxv
        for jj=1:nyv
            if (gn.mask_v(ii,jj)==1)
                L1_v(ii,jj,:,tt)=interp1(ZCL,squeeze(v(ii,jj,:)),squeeze(L1_zv(ii,jj,:)));
            else
                L1_v(ii,jj,:,tt)=nan;
            end
        end
    end
    clear v L1_zv
    
    % Remove possible Nan values
    for zz=1:L1.N
        aa(:,:)=squeeze(L1_temp(:,:,zz,tt));
        aa=maplev(aa);
        L1_temp(:,:,zz,tt)=double(aa);
        clear aa
        
        aa(:,:)=squeeze(L1_salt(:,:,zz,tt));
        aa=maplev(aa);
        L1_salt(:,:,zz,tt)=double(aa);
        clear aa
        
        aa(:,:)=squeeze(L1_u(:,:,zz,tt));
        aa=maplev(aa);
        L1_u(:,:,zz,tt)=double(aa);
        clear aa
        
        aa(:,:)=squeeze(L1_v(:,:,zz,tt));
        aa=maplev(aa);
        L1_v(:,:,zz,tt)=double(aa);
        clear aa
    end
end

toc
% save('E:\FLORIDA_GRIDS\L1\L0_L1_GOMSAB.mat','time','time4d','L1_temp','L1_salt','L1_u','L1_v','L1_ubar_4d','L1_vbar_4d','L1_zeta_4d','L1_ubar','L1_vbar','L1_zeta','-v7.3');

%% CREATE INITIAL CONDITION
init_time=time4d(1,1)-datenum(1858,11,17);
create_roms_init_from_coawst(L1.modelgrid,init_file,init_time,...
    L1.theta_s,L1.theta_b,L1.Tcline,L1.Vtransform,L1.Vstretching,...
    L1.N,L1_u(:,:,:,1),L1_v(:,:,:,1),L1_ubar(:,:,1),L1_vbar(:,:,1),...
    L1_temp(:,:,:,1),L1_salt(:,:,:,1),L1_zeta(:,:,1));

ncwrite(init_file,'hc',L1.Tcline)

%
% clear zeta_coawst ubar_coawst vbar_coawst u_coawst v_coawst temp_coawst salt_coawst
%
%% CREATE CLIMATOLOGY
clm_time=time4d-datenum(1858,11,17);
nt_ts=length(init_time);

create_roms_clm_from_coawst(L1.modelgrid,clm_file,clm_time,...
    L1.theta_s,L1.theta_b,L1.Tcline,L1.Vtransform,L1.Vstretching,...
    L1.N,L1_u,L1_v,L1_ubar_4d,L1_vbar_4d,L1_temp,L1_salt,L1_zeta_4d);

%% CREATE BOUNDARY CONDITION
time_bry=time-datenum(1858,11,17);
time4d_bry=time4d-datenum(1858,11,17);
zeta_north=squeeze(L1_zeta(:,end,:));
ubar_north=squeeze(L1_ubar(:,end,:));
vbar_north=squeeze(L1_vbar(:,end,:));
u_north=squeeze(L1_u(:,end,:,:));
v_north=squeeze(L1_v(:,end,:,:));
salt_north=squeeze(L1_salt(:,end,:,:));
temp_north=squeeze(L1_temp(:,end,:,:));
%
zeta_south=squeeze(L1_zeta(:,1,:));
ubar_south=squeeze(L1_ubar(:,1,:));
vbar_south=squeeze(L1_vbar(:,1,:));
u_south=squeeze(L1_u(:,1,:,:));
v_south=squeeze(L1_v(:,1,:,:));
salt_south=squeeze(L1_salt(:,1,:,:));
temp_south=squeeze(L1_temp(:,1,:,:));
%
zeta_east=squeeze(L1_zeta(end,:,:));
ubar_east=squeeze(L1_ubar(end,:,:));
vbar_east=squeeze(L1_vbar(end,:,:));
u_east=squeeze(L1_u(end,:,:,:));
v_east=squeeze(L1_v(end,:,:,:));
salt_east=squeeze(L1_salt(end,:,:,:));
temp_east=squeeze(L1_temp(end,:,:,:));
%
zeta_west=squeeze(L1_zeta(1,:,:));
ubar_west=squeeze(L1_ubar(1,:,:));
vbar_west=squeeze(L1_vbar(1,:,:));
u_west=squeeze(L1_u(1,:,:,:));
v_west=squeeze(L1_v(1,:,:,:));
salt_west=squeeze(L1_salt(1,:,:,:));
temp_west=squeeze(L1_temp(1,:,:,:));

create_roms_bry_from_coawst(L1.modelgrid,bry_file,time_bry,time4d_bry,...
    L1.theta_s,L1.theta_b,L1.Tcline,L1.Vtransform,L1.Vstretching,L1.N,...
    zeta_north,ubar_north,vbar_north,u_north,v_north,salt_north,temp_north,...
    ubar_south,vbar_south,zeta_south,u_south,v_south,salt_south,temp_south,...
    zeta_east,ubar_east,vbar_east,u_east,v_east,salt_east,temp_east,...
    zeta_west,ubar_west,vbar_west,u_west,v_west,salt_west,temp_west);
end
