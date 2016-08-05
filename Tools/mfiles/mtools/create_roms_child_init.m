function create_roms_child_init(parent_grid,child_grid,parent_ini,child_ini)
%
% Create a netcdf file that contains initialization data for a ROMS 
% child grid. % Initializes temp, salt, u, v, ubar, vbar, and 
% all sediment parameters.
%
% parent_grid - input parent grid
% child_grid  - input child grid
% parent_ini  - input parent init file
% child_ini   - output child init file
%
% jcw 5-25-2005
% jcw 3-7-07 add L and M for get grid
% jcwarner adapt froma Brandy m file, 12Aug2014
%
%!         W-level  RHO-level                                           !
%!                                                                   
%!            N     _________                                           !
%!                 |         |                                          !
%!                 |    N    |                                          !
%!          N-1    |_________|  S                                       !
%!                 |         |  E                                       !
%!                 |   N-1   |  A                                       !
%!            2    |_________|  W                                       !
%!                 |         |  A                                       !
%!                 |    2    |  T                                       !
%!            1    |_________|  E                                       !
%!                 |         |  R                                       !
%!                 |    1    |                                          !
%!            0    |_________|_____ bathymetry                          !
%!                 |/////////|                                          !
%!                 |    1    |                                          !
%!            1    |_________|  S                                       !
%!                 |         |  E                                       !
%!                 |    2    |  D                                       !
%!            2    |_________|  I                                       !
%!                 |         |  M                                       !
%!                 |  Nbed-1 |  E                                       !
%!        Nbed-1   |_________|  N                                       !
%!                 |         |  T                                       !
%!                 |  Nbed   |                                          !
%!         Nbed    |_________|                                          !
%!                                                                      !

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  No NEED for user input section                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%1) Enter name of netcdf initial file to be created.
%   If it already exists it will be overwritten!!.
    init_file=child_ini;

%2) Enter start time of initial file, in seconds.
%   This time needs to be consistent with model time (ie dstart and time_ref).
%   See *.in files for more detail. 
   init_time=ncread(parent_ini,'ocean_time');
   init_time=init_time(1);

%3) Enter number of vertical sigma levels in model.
    Cs_r=ncread(parent_ini,'Cs_r');
    Cs_w=ncread(parent_ini,'Cs_w');
%   s_rho=ncread(parent_ini,'s_rho');
%   s_w=ncread(parent_ini,'s_w');
    N=length(Cs_r);

%4) Enter the values of theta_s, theta_b, and Tcline from your *.in file.
    theta_s=ncread(parent_ini,'theta_s');
    theta_b=ncread(parent_ini,'theta_b');
    Tcline=ncread(parent_ini,'Tcline');
    Vtransform=ncread(parent_ini,'Vtransform');
    Vstretching=ncread(parent_ini,'Vstretching');

%5) Enter value of h, Lm, and Mm.
    h=ncread(child_grid,'h');
    hmin=min(min(h(:)));
    hc=min([hmin,Tcline]);
    [LP,MP]=size(h);
    L  = LP-1;
    M  = MP-1;
    xi_psi  = L;
    xi_rho  = LP;
    xi_u    = L;
    xi_v    = LP;
    eta_psi = M;
    eta_rho = MP;
    eta_u   = MP;
    eta_v   = M;
%
% These are just copied from above, and then we call get_roms_grid.
%
   Sinp.N           =N;            %number of vertical levels
   Sinp.Vtransform  =Vtransform;   %vertical transformation equation
   Sinp.Vstretching =Vstretching;  %vertical stretching function
   Sinp.theta_s     =theta_s;      %surface control parameter
   Sinp.theta_b     =theta_b;      %bottom  control parameter
   Sinp.Tcline      =Tcline;       %surface/bottom stretching width
   Sinp.hc          =hc;           %stretching width used in ROMS
%
   Gout=get_roms_grid(child_grid,Sinp);
%
%6) Initialize zeta, salt, temp, u, v, ubar, vbar.
%
% Init values for zeta.
% display('Initializing zeta')
%
%get zeta from init conditions for us_east grd
lor=ncread(parent_grid,'lon_rho');
lar=ncread(parent_grid,'lat_rho');
lou=ncread(parent_grid,'lon_u');
lau=ncread(parent_grid,'lat_u');
lov=ncread(parent_grid,'lon_v');
lav=ncread(parent_grid,'lat_v');
mar=double(ncread(parent_grid,'mask_rho'));
mau=double(ncread(parent_grid,'mask_u'));
mav=double(ncread(parent_grid,'mask_v'));
mar(mar==0)=nan;
mau(mau==0)=nan;
mav(mav==0)=nan;

display('Initializing zeta')
zt=ncread(parent_ini,'zeta');
if (size(zt)>2)
  zt=squeeze(zt(:,:,1));
end
zt=double(zt.*mar);
zt(zt<-9999)=nan;
zt(zt>9999)=nan;
zeta=griddata(lor,lar,zt,Gout.lon_rho,Gout.lat_rho);
zeta=maplev(squeeze(zeta(:,:)));
clear zt

display('Initializing ubar and vbar')
ub=ncread(parent_ini,'ubar');
if (size(ub)>2)
  ub=squeeze(ub(:,:,1));
end
ub=double(ub.*mau);
ub(ub<-9999)=nan;
ub(ub>9999)=nan;
ubar=griddata(lou,lau,ub,Gout.lon_u,Gout.lat_u);
ubar=maplev(ubar);
clear ub

vb=ncread(parent_ini,'vbar');
if (size(vb)>2)
  vb=squeeze(vb(:,:,1));
end
vb=double(vb.*mav);
vb(vb<-9999)=nan;
vb(vb>9999)=nan;
vbar=griddata(lov,lav,vb,Gout.lon_v,Gout.lat_v);
vbar=maplev(vbar);
clear vb

display('Initializing u')
ub=ncread(parent_ini,'u');
if (size(ub)>3)
  ub=squeeze(ub(:,:,:,1));
end
ub=double(ub);
for k=1:N
    ub2=squeeze(ub(:,:,k).*mau);
    ub2(ub2<-9999)=nan;
    ub2(ub2>9999)=nan;
    u(:,:,k)=griddata(lou,lau,ub2,Gout.lon_u,Gout.lat_u);
end
for k=1:N
  u(:,:,k)=maplev(squeeze(u(:,:,k)));
end
clear ub; clear ub2;

display('Initializing v')
vb=ncread(parent_ini,'v');
if (size(vb)>3)
  vb=squeeze(vb(:,:,:,1));
end
vb=double(vb);
for k=1:N
    vb2=squeeze(vb(:,:,k).*mav);
    vb2(vb2<-9999)=nan;
    vb2(vb2>9999)=nan;
    v(:,:,k)=griddata(lov,lav,vb2,Gout.lon_v,Gout.lat_v);
end
for k=1:N
  v(:,:,k)=maplev(squeeze(v(:,:,k)));
end
clear vb

display('Initializing salt')
sa=ncread(parent_ini,'salt');
if (size(sa)>3)
  sa=squeeze(sa(:,:,:,1));
end
sa=double(sa);
for k=1:N
    sa2=squeeze(sa(:,:,k).*mar);
    sa2(sa2<0)=nan;
    sa2(sa2>9999)=nan;
    salt(:,:,k)=griddata(lor,lar,sa2,Gout.lon_rho,Gout.lat_rho);
end
for k=1:N
  salt(:,:,k)=maplev(squeeze(salt(:,:,k)));
end
clear sa; clear sa2;

display('Initializing temp')
te=ncread(parent_ini,'temp');
if (size(te)>3)
  te=squeeze(te(:,:,:,1));
end
te=double(te);
te(te<=0)=nan;
for k=1:N
    te2=squeeze(te(:,:,k).*mar);
    te2(te2<0)=nan;
    te2(te2>9999)=nan;
    temp(:,:,k)=griddata(lor,lar,te2,Gout.lon_rho,Gout.lat_rho);
end
for k=1:N
  temp(:,:,k)=maplev(squeeze(temp(:,:,k)));
end
clear te; clear te2;
%
%7) Get number of noncohesive and cohesive sediments.
finfo=ncinfo(parent_ini);
Numvars=length(finfo.Variables);
NCS=0;
for j=1:50 %pick some high number here
  zz=['000',num2str(j)];
  zz=zz(end-1:end);
  for i=1:Numvars
    if (strcmp(['mud_',zz],finfo.Variables(i).Name))
      NCS=NCS+1
    end
  end
end
NNS=0;
for j=1:50 %pick some high number here
  zz=['000',num2str(j)];
  zz=zz(end-1:end);
  for i=1:Numvars
    if (strcmp(['sand_',zz],finfo.Variables(i).Name))
      NNS=NNS+1
    end
  end
end
%
% calc sed parameters. Do not alter.
%
   NST = NCS + NNS;     % total number of sed tracers.
   NAT=2;               % Number of active tracers: usually 2 (temp + salt).
   NT = NAT+NST;        % total number of tracers.

%8) Enter number of bed layers
Numdims=length(finfo.Dimensions);
Nbed=0;
for i=1:Numdims
  if (strcmp(['Nbed'],finfo.Dimensions(i).Name))
      Nbed=finfo.Dimensions(i).Length
  end
end

%9) Provide initial sediment properties in water column.
  display('Initializing suspended sediments.')
%
% mud.
%
for idsed=1:NCS
  display('Initializing Mud')
  count=['0',num2str(idsed)];
  count=count(end-1:end);
  sa=ncread(parent_ini,['mud_',count]);
  if (size(sa)>3)
    sa=squeeze(sa(:,:,:,1));
  end
  sa=double(sa);
  for k=1:Nbed
    sa2=squeeze(sa(:,:,k).*mar);
    sa2(sa2<0)=nan;
    sa2(sa2>9999)=nan;
    eval(['mud_',count,'(:,:,k)=griddata(lor,lar,sa2,Gout.lon_rho,Gout.lat_rho);'])
  end
  for k=1:Nbed
    eval(['mud_',count,'(:,:,k)=maplev(squeeze(mud_',count,'(:,:,k)));'])
  end
  clear sa; clear sa2;
end
%
% sand.
%
for idsed=1:NNS
  display('Initializing Sand')
  count=['0',num2str(idsed)];
  count=count(end-1:end);
  sa=ncread(parent_ini,['sand_',count]);
  if (size(sa)>3)
    sa=squeeze(sa(:,:,:,1));
  end
  sa=double(sa);
  for k=1:Nbed
    sa2=squeeze(sa(:,:,k).*mar);
    sa2(sa2<0)=nan;
    sa2(sa2>9999)=nan;
    eval(['sand_',count,'(:,:,k)=griddata(lor,lar,sa2,Gout.lon_rho,Gout.lat_rho);'])
  end
  for k=1:Nbed
    eval(['sand_',count,'(:,:,k)=maplev(squeeze(sand_',count,'(:,:,k)));'])
  end
  clear sa; clear sa2;
end

%10) Provide initial sediment properties in bed.
%   bed properties
if (NST>0)
  display('Initializing bed thickness.')
  zt=ncread(parent_ini,'bed_thickness');
  if (size(zt)>3)
    zt=squeeze(zt(:,:,:,1));
  end
  zt=double(zt);
  for k=1:Nbed
    zt2=squeeze(zt(:,:,k).*mar);
    zt2(zt2<0)=nan;
    zt2(zt2>9999)=nan;
    bed_thickness(:,:,k)=griddata(lor,lar,zt2,Gout.lon_rho,Gout.lat_rho);
  end
  for k=1:Nbed
    bed_thickness(:,:,k)=maplev(squeeze(bed_thickness(:,:,k)));
  end
  clear zt; clear zt2;
%
  display('Initializing bed age.')
  zt=ncread(parent_ini,'bed_age');
  if (size(zt)>3)
    zt=squeeze(zt(:,:,:,1));
  end
  zt=double(zt);
  for k=1:Nbed
    zt2=squeeze(zt(:,:,k).*mar);
    zt2(zt2<0)=nan;
    zt2(zt2>9999999)=nan;
    bed_age(:,:,k)=griddata(lor,lar,zt2,Gout.lon_rho,Gout.lat_rho);
  end
  for k=1:Nbed
    bed_age(:,:,k)=maplev(squeeze(bed_age(:,:,k)));
  end
  clear zt; clear zt2;
%
  display('Initializing bed porosity.')
  zt=ncread(parent_ini,'bed_porosity');
  if (size(zt)>3)
    zt=squeeze(zt(:,:,:,1));
  end
  zt=double(zt);
  for k=1:Nbed
    zt2=squeeze(zt(:,:,k).*mar);
    zt2(zt2<0)=nan;
    zt2(zt2>9999)=nan;
    bed_porosity(:,:,k)=griddata(lor,lar,zt2,Gout.lon_rho,Gout.lat_rho);
  end
  for k=1:Nbed
    bed_porosity(:,:,k)=maplev(squeeze(bed_porosity(:,:,k)));
  end
  clear zt; clear zt2;
%
  display('Initializing bed biodiff.')
  zt=ncread(parent_ini,'bed_biodiff');
  if (size(zt)>3)
    zt=squeeze(zt(:,:,:,1));
  end
  zt=double(zt);
  for k=1:Nbed
    zt2=squeeze(zt(:,:,k).*mar);
    zt2(zt2<0)=nan;
    zt2(zt2>9999)=nan;
    bed_biodiff(:,:,k)=griddata(lor,lar,zt2,Gout.lon_rho,Gout.lat_rho);
  end
  for k=1:Nbed
    bed_biodiff(:,:,k)=maplev(squeeze(bed_biodiff(:,:,k)));
  end
  clear zt; clear zt2;
end
%
% for mud
%
for idsed=1:NCS
  display('Initializing Mudmass and Mudfrac')
  count=['0',num2str(idsed)];
  count=count(end-1:end);
  sa=ncread(parent_ini,['mudfrac_',count]);
  if (size(sa)>3)
    sa=squeeze(sa(:,:,:,1));
  end
  sa=double(sa);
  for k=1:Nbed
    sa2=squeeze(sa(:,:,k).*mar);
    sa2(sa2<0)=nan;
    sa2(sa2>9999)=nan;
    eval(['mudfrac_',count,'(:,:,k)=griddata(lor,lar,sa2,Gout.lon_rho,Gout.lat_rho);'])
  end
  for k=1:Nbed
    eval(['mudfrac_',count,'(:,:,k)=maplev(squeeze(mudfrac_',count,'(:,:,k)));'])
    eval(['mudmass_',count,'(:,:,k)=squeeze(mudfrac_',count,'(:,:,k)).*bed_thickness(:,:,k)*2650.*(1.0-bed_porosity(:,:,k));'])
  end
  clear sa; clear sa2;
end
%
% for sand
%
for idsed=1:NNS
  display('Initializing Sandmass and Sandfrac')
  count=['0',num2str(idsed)];
  count=count(end-1:end);
  sa=ncread(parent_ini,['sandfrac_',count]);
  if (size(sa)>3)
    sa=squeeze(sa(:,:,:,1));
  end
  sa=double(sa);
  for k=1:Nbed
    sa2=squeeze(sa(:,:,k).*mar);
    sa2(sa2<0)=nan;
    sa2(sa2>9999)=nan;
    eval(['sandfrac_',count,'(:,:,k)=griddata(lor,lar,sa2,Gout.lon_rho,Gout.lat_rho);'])
  end
  for k=1:Nbed
    eval(['sandfrac_',count,'(:,:,k)=maplev(squeeze(sandfrac_',count,'(:,:,k)));'])
    eval(['sandmass_',count,'(:,:,k)=squeeze(sandfrac_',count,'(:,:,k)).*bed_thickness(:,:,k)*2650.*(1.0-bed_porosity(:,:,k));'])
  end
  clear sa; clear sa2;
%
  eval(['bedload_Usand_',count,'(1:xi_u,1:eta_u,1:length(init_time)) = 0;'])              %bed load
  eval(['bedload_Vsand_',count,'(1:xi_v,1:eta_v,1:length(init_time)) = 0;'])              %bed load
end

%11)
% set some surface properties
if (NST>0)
  display('Initializing sediment surface properties.')
%
  zt=ncread(parent_ini,'grain_diameter');
  if (size(zt)>2)
    zt=squeeze(zt(:,:,1));
  end
  zt=double(zt.*mar);
  zt(zt<-9999)=nan;
  zt(zt>9999)=nan;
  grain_diameter=griddata(lor,lar,zt,Gout.lon_rho,Gout.lat_rho);
  grain_diameter=maplev(squeeze(grain_diameter(:,:)));
  clear zt
%
  zt=ncread(parent_ini,'grain_density');
  if (size(zt)>2)
    zt=squeeze(zt(:,:,1));
  end
  zt=double(zt.*mar);
  zt(zt<-9999)=nan;
  zt(zt>9999)=nan;
  grain_density=griddata(lor,lar,zt,Gout.lon_rho,Gout.lat_rho);
  grain_density=maplev(squeeze(grain_density(:,:)));
  clear zt
%
  zt=ncread(parent_ini,'settling_vel');
  if (size(zt)>2)
    zt=squeeze(zt(:,:,1));
  end
  zt=double(zt.*mar);
  zt(zt<-9999)=nan;
  zt(zt>9999)=nan;
  settling_vel=griddata(lor,lar,zt,Gout.lon_rho,Gout.lat_rho);
  settling_vel=maplev(squeeze(settling_vel(:,:)));
  clear zt
%
  zt=ncread(parent_ini,'erosion_stress');
  if (size(zt)>2)
    zt=squeeze(zt(:,:,1));
  end
  zt=double(zt.*mar);
  zt(zt<-9999)=nan;
  zt(zt>9999)=nan;
  erosion_stress=griddata(lor,lar,zt,Gout.lon_rho,Gout.lat_rho);
  erosion_stress=maplev(squeeze(erosion_stress(:,:)));
  clear zt
%
  zt=ncread(parent_ini,'ripple_length');
  if (size(zt)>2)
    zt=squeeze(zt(:,:,1));
  end
  zt=double(zt.*mar);
  zt(zt<-9999)=nan;
  zt(zt>9999)=nan;
  ripple_length=griddata(lor,lar,zt,Gout.lon_rho,Gout.lat_rho);
  ripple_length=maplev(squeeze(ripple_length(:,:)));
  clear zt
%
  zt=ncread(parent_ini,'ripple_height');
  if (size(zt)>2)
    zt=squeeze(zt(:,:,1));
  end
  zt=double(zt.*mar);
  zt(zt<-9999)=nan;
  zt(zt>9999)=nan;
  ripple_height=griddata(lor,lar,zt,Gout.lon_rho,Gout.lat_rho);
  ripple_height=maplev(squeeze(ripple_height(:,:)));
  clear zt
%
  zt=ncread(parent_ini,'dmix_offset');
  if (size(zt)>2)
    zt=squeeze(zt(:,:,1));
  end
  zt=double(zt.*mar);
  zt(zt<-9999)=nan;
  zt(zt>9999)=nan;
  dmix_offset=griddata(lor,lar,zt,Gout.lon_rho,Gout.lat_rho);
  dmix_offset=maplev(squeeze(dmix_offset(:,:)));
  clear zt
%
  zt=ncread(parent_ini,'dmix_slope');
  if (size(zt)>2)
    zt=squeeze(zt(:,:,1));
  end
  zt=double(zt.*mar);
  zt(zt<-9999)=nan;
  zt(zt>9999)=nan;
  dmix_slope=griddata(lor,lar,zt,Gout.lon_rho,Gout.lat_rho);
  dmix_slope=maplev(squeeze(dmix_slope(:,:)));
  clear zt
%
  zt=ncread(parent_ini,'dmix_time');
  if (size(zt)>2)
    zt=squeeze(zt(:,:,1));
  end
  zt=double(zt.*mar);
  zt(zt<-9999)=nan;
  zt(zt>9999)=nan;
  dmix_time=griddata(lor,lar,zt,Gout.lon_rho,Gout.lat_rho);
  dmix_time=maplev(squeeze(dmix_time(:,:)));
  clear zt
end

%12 Vegetation
Numdims=length(finfo.Dimensions);
NVEG=1;
for i=1:Numdims
  if (strcmp(['Nveg'],finfo.Dimensions(i).Name))
      NVEG=finfo.Dimensions(i).Length
  end
end

if (NVEG>0)
  display('Initializing vegetation properties.')
%
  zt=ncread(parent_ini,'plant_density');
  if (size(zt)>2)
    zt=squeeze(zt(:,:,1));
  end
  zt=double(zt.*mar);
  zt(zt<-9999)=nan;
  zt(zt>9999)=nan;
  plant_density=griddata(lor,lar,zt,Gout.lon_rho,Gout.lat_rho);
  plant_density=maplev(squeeze(plant_density(:,:)));
  clear zt
%
  zt=ncread(parent_ini,'plant_height');
  if (size(zt)>2)
    zt=squeeze(zt(:,:,1));
  end
  zt=double(zt.*mar);
  zt(zt<-9999)=nan;
  zt(zt>9999)=nan;
  plant_height=griddata(lor,lar,zt,Gout.lon_rho,Gout.lat_rho);
  plant_height=maplev(squeeze(plant_height(:,:)));
  clear zt
%
  zt=ncread(parent_ini,'plant_diameter');
  if (size(zt)>2)
    zt=squeeze(zt(:,:,1));
  end
  zt=double(zt.*mar);
  zt(zt<-9999)=nan;
  zt(zt>9999)=nan;
  plant_diameter=griddata(lor,lar,zt,Gout.lon_rho,Gout.lat_rho);
  plant_diameter=maplev(squeeze(plant_diameter(:,:)));
  clear zt
%
  zt=ncread(parent_ini,'plant_thickness');
  if (size(zt)>2)
    zt=squeeze(zt(:,:,1));
  end
  zt=double(zt.*mar);
  zt(zt<-9999)=nan;
  zt(zt>9999)=nan;
  plant_thickness=griddata(lor,lar,zt,Gout.lon_rho,Gout.lat_rho);
  plant_thickness=maplev(squeeze(plant_thickness(:,:)));
  clear zt
%
  zt=ncread(parent_ini,'marsh_mask');
  if (size(zt)>2)
    zt=squeeze(zt(:,:,1));
  end
  zt=double(zt.*mar);
  zt(zt<-9999)=nan;
  zt(zt>9999)=nan;
  marsh_mask=griddata(lor,lar,zt,Gout.lon_rho,Gout.lat_rho);
  marsh_mask=maplev(squeeze(marsh_mask(:,:)));
  clear zt
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%create init file
create_roms_netcdf_init_mw(init_file,Gout,Nbed,NNS,NCS,NVEG)

ncwrite(init_file,'theta_s',Gout.theta_s);
ncwrite(init_file,'theta_b',Gout.theta_b);
ncwrite(init_file,'Tcline',Gout.Tcline);
ncwrite(init_file,'Cs_r',Gout.Cs_r);
ncwrite(init_file,'Cs_w',Gout.Cs_w);
ncwrite(init_file,'sc_w',Gout.s_w);
ncwrite(init_file,'sc_r',Gout.s_rho);
ncwrite(init_file,'hc',Gout.hc);
ncwrite(init_file,'Vtransform',Gout.Vtransform);
ncwrite(init_file,'Vstretching',Gout.Vstretching);
ncwrite(init_file,'spherical',Gout.spherical);

ncwrite(init_file,'ocean_time',init_time);

ncwrite(init_file,'temp',temp);
ncwrite(init_file,'salt',salt);
ncwrite(init_file,'u',u);
ncwrite(init_file,'ubar',ubar);
ncwrite(init_file,'v',v);
ncwrite(init_file,'vbar',vbar);
ncwrite(init_file,'zeta',zeta);

if (NST>0)
  for mm=1:NCS
    count=['00',num2str(mm)];
    count=count(end-1:end);
    eval(['ncwrite(init_file,''mud_',count,''',mud_',count,');'])           %sed conc in water column
    eval(['ncwrite(init_file,''mudfrac_',count,''',mudfrac_',count,');'])           %sed conc in water column
    eval(['ncwrite(init_file,''mudmass_',count,''',mudmass_',count,');'])           %sed conc in water column
  end
  for mm=1:NNS
    count=['00',num2str(mm)];
    count=count(end-1:end);
    eval(['ncwrite(init_file,''sand_',count,''',sand_',count,');'])           %sed conc in water column
    eval(['ncwrite(init_file,''sandfrac_',count,''',sandfrac_',count,');'])           %sed conc in water column
    eval(['ncwrite(init_file,''sandmass_',count,''',sandmass_',count,');'])           %sed conc in water column
%
    eval(['ncwrite(init_file,''bedload_Usand_',count,''',bedload_Usand_',count,');'])   %bedload
    eval(['ncwrite(init_file,''bedload_Vsand_',count,''',bedload_Vsand_',count,');'])   %bedload
  end

  ncwrite(init_file,'bed_thickness',bed_thickness);
  ncwrite(init_file,'bed_age',bed_age);
  ncwrite(init_file,'bed_porosity',bed_porosity);
  ncwrite(init_file,'bed_biodiff',bed_biodiff);

  ncwrite(init_file,'grain_diameter',grain_diameter);
  ncwrite(init_file,'grain_density',grain_density);
  ncwrite(init_file,'settling_vel',settling_vel);
  ncwrite(init_file,'erosion_stress',erosion_stress); 
  ncwrite(init_file,'ripple_height',ripple_height);
  ncwrite(init_file,'ripple_length',ripple_length);
  ncwrite(init_file,'dmix_offset',dmix_offset);
  ncwrite(init_file,'dmix_slope',dmix_slope);
  ncwrite(init_file,'dmix_time',dmix_time);
end
if (NVEG>0)
  ncwrite(init_file,'plant_height',plant_height);
  ncwrite(init_file,'plant_diameter',plant_diameter);
  ncwrite(init_file,'plant_density',plant_density);
  ncwrite(init_file,'plant_thickness',plant_thickness);
  ncwrite(init_file,'marsh_mask',marsh_mask);
end
