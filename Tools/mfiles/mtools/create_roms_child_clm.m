function create_roms_child_clm(parent_grid,child_grid,parent_clm,child_clm)
%
% Create a netcdf file that contains climatology data for a ROMS 
% child grid. Reads temp, salt, u, v, ubar, vbar, and 
% all sediment parameters.
%
% parent_grid - input parent grid
% child_grid  - input child grid
% parent_clm  - input parent climatology file
% child_clm   - output child climatology file
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
    clm_file=child_clm;

%2) Enter start time of initial file, in seconds.
%   This time needs to be consistent with model time (ie dstart and time_ref).
%   See *.in files for more detail. 
    clm_time=ncread(parent_clm,'ocean_time');

%3) Enter number of vertical sigma levels in model.
     N=ncread(parent_clm,'temp',[1 1 1 1],[1 1 inf 1]);
     N=max(size(N)); gn.N=N;

%4) Enter the values of theta_s, theta_b, and Tcline from your *.in file.
%    theta_s=ncread(parent_clm,'theta_s');
%    theta_b=ncread(parent_clm,'theta_b');
%    Tcline=ncread(parent_clm,'Tcline');
%    Vtransform=ncread(parent_clm,'Vtransform');
%    Vstretching=ncread(parent_clm,'Vstretching');

%5) Enter value of h, Lm, and Mm.
      h=ncread(child_grid,'h');
%     hmin=min(min(h(:)));
%     hc=min([hmin,Tcline]);
      lon_rho=ncread(child_grid,'lon_rho');gn.lon_rho=lon_rho;
      lat_rho=ncread(child_grid,'lat_rho');gn.lat_rho=lat_rho;
      lon_u=ncread(child_grid,'lon_u');
      lat_u=ncread(child_grid,'lat_u');
      lon_v=ncread(child_grid,'lon_v');
      lat_v=ncread(child_grid,'lat_v');
      [LP,MP]=size(h);
      L  = LP-1;
      M  = MP-1;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calc some grid stuff here - do not change this.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   xi_psi  = L;
   xi_rho  = LP;
   xi_u    = L;
   xi_v    = LP;
   eta_psi = M;
   eta_rho = MP;
   eta_u   = MP;
   eta_v   = M;

%create clm file
create_roms_netcdf_clm_mwUL(clm_file,gn,length(clm_time))
%ncwrite(clm_file,'theta_s',theta_s);
%ncwrite(clm_file,'theta_b',theta_b);
%ncwrite(clm_file,'Tcline',Tcline);
%ncwrite(clm_file,'Cs_r',Cs_r);
%ncwrite(clm_file,'Cs_w',Cs_w);
%ncwrite(clm_file,'sc_w',sc_w);
%ncwrite(clm_file,'sc_r',sc_r);
%%ncwrite(clm_file,'hc',hc);
%ncwrite(clm_file,'Vtransform',Vtransform);
%ncwrite(clm_file,'Vstretching',Vstretching);
%ncwrite(clm_file,'ocean_time',clm_time);


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

display('Initializing zeta')
for mm=1:length(clm_time)
  zt=ncread(parent_clm,'zeta',[1 1 mm],[inf inf 1]);
  zt=double(zt);
  zt(zt<-9999)=nan;
  zt(zt>9999)=nan;
  zeta(:,:,mm)=griddata(lor,lar,zt,lon_rho,lat_rho);
  zeta(:,:,mm)=maplev(squeeze(zeta(:,:,mm)));
  clear zt
end
ncwrite(clm_file,'zeta',zeta);
 
display('Initializing ubar and vbar')
for mm=1:length(clm_time)
  ub=ncread(parent_clm,'ubar',[1 1 mm],[inf inf 1]);
  ub=double(ub);
  ub(ub<-9999)=nan;
  ub(ub>9999)=nan;
  ubar(:,:,mm)=griddata(lou,lau,ub,lon_u,lat_u);
  ubar(:,:,mm)=maplev(squeeze(ubar(:,:,mm)));
  clear ub
end
ncwrite(clm_file,'ubar',ubar);

for mm=1:length(clm_time)
  vb=ncread(parent_clm,'vbar',[1 1 mm],[inf inf 1]);
  vb=double(vb);
  vb(vb<-9999)=nan;
  vb(vb>9999)=nan;
  vbar(:,:,mm)=griddata(lov,lav,vb,lon_v,lat_v);
  vbar(:,:,mm)=maplev(squeeze(vbar(:,:,mm)));
  clear vb
end
ncwrite(clm_file,'vbar',vbar);

display('Initializing u')
for mm=1:length(clm_time)
  ub=ncread(parent_clm,'u',[1 1 1 mm],[inf inf inf 1]);
  ub=double(ub);
  for k=1:N
    ub2=squeeze(ub(:,:,k));
    ub2(ub2<-9999)=nan;
    ub2(ub2>9999)=nan;
    u(:,:,k,mm)=griddata(lou,lau,ub2,lon_u,lat_u);
    u(:,:,k,mm)=maplev(squeeze(u(:,:,k,mm)));
    clear ub2
  end
  clear ub
end
ncwrite(clm_file,'u',u);

display('Initializing v')
for mm=1:length(clm_time)
  vb=ncread(parent_clm,'v',[1 1 1 mm],[inf inf inf 1]);
  vb=double(vb);
  for k=1:N
    vb2=squeeze(vb(:,:,k));
    vb2(vb2<-9999)=nan;
    vb2(vb2>9999)=nan;
    v(:,:,k,mm)=griddata(lov,lav,vb2,lon_v,lat_v);
    v(:,:,k,mm)=maplev(squeeze(v(:,:,k,mm)));
    clear vb2
  end
  clear vb
end
ncwrite(clm_file,'v',v);

display('Initializing salt')
for mm=1:length(clm_time)
  sa=ncread(parent_clm,'salt',[1 1 1 mm],[inf inf inf 1]);
  sa=double(sa);
  for k=1:N
    sa2=squeeze(sa(:,:,k));
    sa2(sa2<-9999)=nan;
    sa2(sa2>9999)=nan;
    salt(:,:,k,mm)=griddata(lor,lar,sa2,lon_rho,lat_rho);
    salt(:,:,k,mm)=maplev(squeeze(salt(:,:,k,mm)));
    clear sa2
  end
  clear sa
end
ncwrite(clm_file,'salt',salt);

display('Initializing temp')
for mm=1:length(clm_time)
  te=ncread(parent_clm,'temp',[1 1 1 mm],[inf inf inf 1]);
  te=double(te);
  for k=1:N
    te2=squeeze(te(:,:,k));
    te2(te2<-9999)=nan;
    te2(te2>9999)=nan;
    temp(:,:,k,mm)=griddata(lor,lar,te2,lon_rho,lat_rho);
    temp(:,:,k,mm)=maplev(squeeze(temp(:,:,k,mm)));
    clear te2
  end
  clear te
end
ncwrite(clm_file,'temp',temp);

ncwrite(clm_file,'ocean_time',clm_time)
ncwrite(clm_file,'zeta_time',clm_time)
ncwrite(clm_file,'v2d_time',clm_time)
ncwrite(clm_file,'v3d_time',clm_time)
ncwrite(clm_file,'salt_time',clm_time)
ncwrite(clm_file,'temp_time',clm_time)

