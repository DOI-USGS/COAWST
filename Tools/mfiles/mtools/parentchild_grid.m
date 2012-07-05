% function create_fine_grid(ncfile_coarse,ncfile_fine,Istr,Iend,Jstr,Jend,scale)
%
% called from create_child_grid.m
%
% user provides: - netcdf file name of coarse grid
%                - netcdf file name for new grid
%                - scale factor of increased resolution
%
% output is a netcdf roms grid with finer resolution
%
% User will need to modify the masking and bathymetry of new grid!!!!!
%
% Version jcwarner Aug 3, 2007
% Updated clea denamiel Feb 5, 2012
% Updated jcwarner 26Jun2012, make this a main driver.

%get the coarser grid
netcdf_load(ncfile_coarse)

%set some arrays for interpolation here, psi points
% _c for coarse, _f for fine
offset=ceil(4/scale);  % need 4 points to the left, 3 to the right
numx_c=((Iend+1)-(Istr-offset))+1;
numy_c=((Jend+1)-(Jstr-offset))+1;
numx_f=scale*(numx_c-1)+1;
numy_f=scale*(numy_c-1)+1;

xc_c=[1:numx_c];
xc_c=repmat(xc_c',1,numy_c);
yc_c=[1:numy_c];
yc_c=repmat(yc_c,numx_c,1);

xc_f=[1:1/scale:numx_c];
xc_f=repmat(xc_f',1,numy_f);
yc_f=[1:1/scale:numy_c];
yc_f=repmat(yc_f,numx_f,1);

%%%%%%%%%  Lon Lat SPACE  %%%%%%%%%%%
%establish a full grid of psi points, including an extra row and column
x_full_grid=interp2(xc_c', yc_c', lon_psi(Istr-offset:Iend+1,Jstr-offset:Jend+1)', xc_f, yc_f);
y_full_grid=interp2(xc_c', yc_c', lat_psi(Istr-offset:Iend+1,Jstr-offset:Jend+1)', xc_f, yc_f);

%set the rho points
ztemp=interp2(x_full_grid(offset+2-1:end-(scale-3+1)+1,offset+2-1:end-(scale-3+1)+1),1);
fine.lon.rho=ztemp(2:2:end-1,2:2:end-1);
ztemp=interp2(y_full_grid(offset+2-1:end-(scale-3+1)+1,offset+2-1:end-(scale-3+1)+1),1);
fine.lat.rho=ztemp(2:2:end-1,2:2:end-1);

%%%%%%%%%  X Y SPACE  %%%%%%%%%%%
%establish a full grid of psi points, including an extra row and column
x_full_grid=interp2(xc_c', yc_c', x_psi(Istr-offset:Iend+1,Jstr-offset:Jend+1)', xc_f, yc_f);
y_full_grid=interp2(xc_c', yc_c', y_psi(Istr-offset:Iend+1,Jstr-offset:Jend+1)', xc_f, yc_f);

%set the rho points
ztemp=interp2(x_full_grid(offset+2-1:end-(scale-3+1)+1,offset+2-1:end-(scale-3+1)+1),1);
fine.x.rho=ztemp(2:2:end-1,2:2:end-1);
ztemp=interp2(y_full_grid(offset+2-1:end-(scale-3+1)+1,offset+2-1:end-(scale-3+1)+1),1);
fine.y.rho=ztemp(2:2:end-1,2:2:end-1);

%%%%%% GRID h and masking  %%%%%%%%%
[XI,YI]=meshgrid([min(x_psi(:)):min(1./pm(:))/2:max(x_psi(:))], ...
                 [min(y_psi(:)):min(1./pn(:))/2:max(y_psi(:))]);
ZI=griddata(x_rho,y_rho,h,XI,YI);
%ZI=interp2(x_rho,y_rho,h,XI,YI);
fine.h=interp2(XI,YI,ZI,fine.x.rho,fine.y.rho);

ZI=griddata(x_rho,y_rho,mask_rho,XI,YI);
%ZI=interp2(x_rho,y_rho,mask_rho,XI,YI);
fine.mask.rho=interp2(XI,YI,ZI,fine.x.rho,fine.y.rho,'nearest');
%force masking to = 0 if it is less than 1.
fine.mask.rho(fine.mask.rho>0.95)=1;
fine.mask.rho(fine.mask.rho<0.96)=0;

%%%%%%%%%%%%%%%%%%% finished creating grid %%%%%%%%%%%%%%%%

%call routine to create netcdf file

  projection='mercator';
  rho.x=fine.x.rho;
  rho.y=fine.y.rho;  
  rho.lon=fine.lon.rho;  
  rho.lat=fine.lat.rho;  
  rho.depth=fine.h;
  rho.mask = fine.mask.rho;
 
  save temp_jcw33.mat rho spherical projection
  eval(['mat2roms_mw(''temp_jcw33.mat'',''',ncfile_fine,''');'])
  !del temp_jcw33.mat
  disp(['Created roms grid -->   ',ncfile_fine])
  disp([' '])
  disp(['you really need to look at the masking and bathy in the new file '])
  

%to create swan grid then use:
  disp([' '])
  if (spherical)
    roms2swan(fine.lon.rho,fine.lat.rho,fine.h,fine.mask.rho);
  else
    roms2swan(fine.x.rho,fine.y.rho,fine.h,fine.mask.rho);
  end
  disp(['Created swan grid + bathy files: swan_bathy.bot and swan_coord.grd'])
  disp(['You should probably rename these so they dont get overwritten. '])

%%%%%%%%%%%%% END OF USER SECTION %%%%%%%%%%%%%%%%

% Plot new grid points inside old grid
for mm=1:2
  if (mm==1); netcdf_load(ncfile_coarse); end
  if (mm==2); netcdf_load(ncfile_fine); end
  if (spherical)
    xr=lon_rho;  yr=lat_rho;
    xp=lon_psi;  yp=lat_psi;
    xu=lon_u;    yu=lat_u;
    xv=lon_v;    yv=lat_v;
  else     
    xr=x_rho;  yr=y_rho;
    xp=x_psi;  yp=y_psi;
    xu=x_u;    yu=y_u;
    xv=x_v;    yv=y_v;
  end
  maskr=mask_rho;
  maskr(mask_rho==0)=nan;
  
% Grid lines
  figure(1)
  if (mm==1) 
    plot(xp,yp,'c')
    hold on
    plot(xp',yp','c')
    x_rhoc=xr; y_rhoc=yr;  %save these to plot later
  else
    plot(xp,yp,'k')
    plot(xp',yp','k')
    plot(xr,yr,'ro')
    plot(xu,yu,'kx')
    plot(xv,yv,'ms')
    plot(x_rhoc,y_rhoc,'bo') %replot coarse rhos points on top
  end
  title('grid lines and rho, u, v, points')
% Bathymetry
  figure(2)
  pcolorjw(xr,yr,h.*maskr)
  hold on
  shading flat
  colorbar
  if (mm==2)
    title('bathymetry (m)')
    plot(xp(:,1),yp(:,1),'r'); plot(xp(:,end),yp(:,end),'r'); 
    plot(xp(1,:),yp(1,:),'r'); plot(xp(end,:),yp(end,:),'r'); 
  end
  
% PM
  figure(3)
  pcolorjw(xr,yr,pm.*maskr)
  shading flat
  hold on
  colorbar
  if (mm==2)
    title('pm');
    plot(xp(:,1),yp(:,1),'r'); plot(xp(:,end),yp(:,end),'r'); 
    plot(xp(1,:),yp(1,:),'r'); plot(xp(end,:),yp(end,:),'r'); 
  end

% PN
  figure(4)
  pcolorjw(xr,yr,pn.*maskr)
  shading flat
  hold on
  colorbar
  if (mm==2)
    title('pn');
    plot(xp(:,1),yp(:,1),'r'); plot(xp(:,end),yp(:,end),'r'); 
    plot(xp(1,:),yp(1,:),'r'); plot(xp(end,:),yp(end,:),'r'); 
  end

% dndx
  figure(5)
  pcolorjw(xr,yr,dndx.*maskr)
  shading flat
  hold on
  colorbar
  if (mm==2)
    title('dndx');
    plot(xp(:,1),yp(:,1),'r'); plot(xp(:,end),yp(:,end),'r'); 
    plot(xp(1,:),yp(1,:),'r'); plot(xp(end,:),yp(end,:),'r'); 
  end

% dmde
  figure(6)
  pcolorjw(xr,yr,dmde.*maskr)
  shading flat
  colorbar
  hold on
  if (mm==2)
    title('dmde');
    plot(xp(:,1),yp(:,1),'r'); plot(xp(:,end),yp(:,end),'r'); 
    plot(xp(1,:),yp(1,:),'r'); plot(xp(end,:),yp(end,:),'r'); 
  end

% angle
  figure(7)
  pcolorjw(xr,yr,angle.*maskr)
  shading flat
  hold on
  colorbar
  if (mm==2);
    title('angle');
    plot(xp(:,1),yp(:,1),'r'); plot(xp(:,end),yp(:,end),'r'); 
    plot(xp(1,:),yp(1,:),'r'); plot(xp(end,:),yp(end,:),'r'); 
  end
end

