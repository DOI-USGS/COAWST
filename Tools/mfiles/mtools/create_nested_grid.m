% create_nested_grid.m
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

%%%%%%%%%%%%% START OF USER SECTION %%%%%%%%%%%%%%%%

%1) ENTER NAME OF EXISITNG COARSE GRID FILE
%ncfile_coarse='refined_chan_grid.nc';
ncfile_coarse='inlet_test_grid.nc';

%2) ENTER NAME OF NEW FINE GRID FILE
%ncfile_fine='refined_chan_grid_ref5.nc';
ncfile_fine='inlet_test_grid_ref5_test.nc';

%3) ENTER THE START AND END INDICES OF THE PSI POINTS OF THE COARSE GRID
%   THAT IDENTIFY THE OUTER BOUNDS OF THE FINE GRID
%   Mathworks netcdf is Fortran convention. To determine these locations, use:
%   netcdf_load (ncfile_coarse)
%   figure; plot(x_psi,y_psi,'k'); hold on; plot(x_psi', y_psi','k')
%   and select [Istr, Iend] as the lower left and lower right horizontal indices
%   and select [Jstr, Jend] as the lower left and upper left vertical indices.
% 
%Istr=30; Iend=50; Jstr=2; Jend=5;  % test_chan_refined
Istr=40; Iend=56; Jstr=24; Jend=54;  % inlet_test_refined

%4) ENTER SCALE FACTOR FOR INCREASED RESOLUTION (use 3 or 5)
scale=5;

%%%%%%%%%%%%% END OF USER SECTION %%%%%%%%%%%%%%%%

%set some random variables here
msize1=5;      %marker size for plots

%get the coarser grid
netcdf_load(ncfile_coarse)

%number of rho points in coarse grid section to have increased resolution
Ilen_c=Iend-Istr;
Jlen_c=Jend-Jstr;

%Now set dimensions for fine grid.
%xi direction (left-right)
LP = Ilen_c*scale+7;   %rho dimension
L = LP-1;              %psi dimension

%eta direction (up-down)
MP = Jlen_c*scale+7;   % rho dimension
M = MP-1;              % psi dimension

%sizes for fine grid
xi_psi = L;
xi_rho = LP;
xi_u = L;
xi_v = LP;

eta_psi = M;
eta_rho = MP;
eta_u = MP;
eta_v = M;

%set some arrays for interpolation here, psi points
% _c for coarse, _f for fine
offset=ceil(4/scale);  % need 4 points to the left, 3 to the right
numx_c=((Iend+1)-(Istr-offset))+1;
numy_c=((Jend+1)-(Jstr-offset))+1;
numx_f=scale*(numx_c-1)+1;
numy_f=scale*(numy_c-1)+1;

xc_c=[1:numx_c];
xc_c=repmat(xc_c,numy_c,1);
yc_c=[1:numy_c];
yc_c=repmat(yc_c',1,numx_c);

xc_f=[1:1/scale:numx_c];
xc_f=repmat(xc_f,numy_f,1);
yc_f=[1:1/scale:numy_c];
yc_f=repmat(yc_f',1,numx_f);

%%%%%%%%%  Lon Lat SPACE  %%%%%%%%%%%
%establish a full grid of psi points, including an extra row and column
x_full_grid=interp2(xc_c, yc_c, lon_psi(Jstr-offset:Jend+1,Istr-offset:Iend+1), xc_f, yc_f);
y_full_grid=interp2(xc_c, yc_c, lat_psi(Jstr-offset:Jend+1,Istr-offset:Iend+1), xc_f, yc_f);

%set the rho points
ztemp=interp2(x_full_grid(offset+2-1:end-(scale-3+1)+1,offset+2-1:end-(scale-3+1)+1),1);
fine.lon.rho=ztemp(2:2:end-1,2:2:end-1);
ztemp=interp2(y_full_grid(offset+2-1:end-(scale-3+1)+1,offset+2-1:end-(scale-3+1)+1),1);
fine.lat.rho=ztemp(2:2:end-1,2:2:end-1);


%%%%%%%%%  X Y SPACE  %%%%%%%%%%%
%establish a full grid of psi points, including an extra row and column
x_full_grid=interp2(xc_c, yc_c, x_psi(Jstr-offset:Jend+1,Istr-offset:Iend+1), xc_f, yc_f);
y_full_grid=interp2(xc_c, yc_c, y_psi(Jstr-offset:Jend+1,Istr-offset:Iend+1), xc_f, yc_f);

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
fine.mask.rho=interp2(XI,YI,ZI,fine.x.rho,fine.y.rho);
%force masking to = 0 if it is less than 1.
[js,is]=size(fine.mask.rho);
for i=1:is
  for j=1:js
    if (fine.mask.rho(j,i) < 1)
      fine.mask.rho(j,i)=0;
    end
  end
end
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
  disp(['Created swan grid + bathy files:'])
  roms2swan(fine.x.rho,fine.y.rho,fine.h,fine.mask.rho);

%%%%%%%%%%%%% END OF USER SECTION %%%%%%%%%%%%%%%%

% Plot new grid points inside old grid
for mm=1:2
  if (mm==1); netcdf_load(ncfile_coarse); end
  if (mm==2); netcdf_load(ncfile_fine); end

% Grid lines
  figure(1)
  if (mm==1) 
    plot(x_psi,y_psi,'c')
    hold on
    plot(x_psi',y_psi','c')
    x_rhoc=x_rho; y_rhoc=y_rho;
  else
    plot(x_psi,y_psi,'k')
    plot(x_psi',y_psi','k')
    plot(x_rho,y_rho,'ro')
    plot(x_u,y_u,'kx')
    plot(x_v,y_v,'ms')
    plot(x_rhoc,y_rhoc,'bo')
  end

% Bathymetry
  figure(2)
  subplot(2,1,mm);
  pcolor(x_rho,y_rho,h)
  shading flat
  colorbar
  if (mm==1); title('Parent grid - bathymetry (m)'); end
  if (mm==2); title(' Child grid - bathymetry (m)'); end

% PM
  figure(3)
  subplot(2,1,mm);
  pcolor(x_rho,y_rho,pm)
  shading flat
  colorbar
  if (mm==1); title('Parent grid - pm'); end
  if (mm==2); title(' Child grid - pm'); end


% PN
  figure(4)
  subplot(2,1,mm);
  pcolor(x_rho,y_rho,pn)
  shading flat
  colorbar
  if (mm==1); title('Parent grid - pn'); end
  if (mm==2); title(' Child grid - pn'); end


% dndx
  figure(5)
  subplot(2,1,mm);
  pcolor(x_rho,y_rho,dndx)
  shading flat
  colorbar
  if (mm==1); title('Parent grid - dndx'); end
  if (mm==2); title(' Child grid - dndx'); end

% dmde
  figure(6)
  subplot(2,1,mm);
  pcolor(x_rho,y_rho,dmde)
  shading flat
  colorbar
  if (mm==1); title('Parent grid - dmde'); end
  if (mm==2); title(' Child grid - dmde'); end

% angle
  figure(7)
  subplot(2,1,mm);
  pcolor(x_rho,y_rho,angle)
  shading flat
  colorbar
  if (mm==1); title('Parent grid - angle'); end
  if (mm==2); title(' Child grid - angle'); end
end

