% parentchild_bathy
%
% Loads in parent and child bathymetry and modifies them such that:
% - boundary points of the child grid are adjusted to match the parent,
% - overlapping region in the parent is repalced with average child
%
% Users need to still decide on land masking independently of this script.
% Need to specify parent and child grids.
%

%specify corners - same as used for create_nested_grid
LLi=Istr;
LLj=Jstr;
URi=Iend;
URj=Jend;
ULi=Istr;
ULj=Jend;
LRi=Iend;
LRj=Jstr;

%load coordinates and bathymetry
netcdf_load(ncfile_coarse)
lon_psi_coarse=lon_psi;
lat_psi_coarse=lat_psi;
lon_rho_coarse=lon_rho;
lat_rho_coarse=lat_rho;
h_coarse=h;
h_coarse_new=h_coarse;
%%
netcdf_load(ncfile_fine)
lon_psi_fine=lon_psi;
lat_psi_fine=lat_psi;
lon_rho_fine=lon_rho;
lat_rho_fine=lat_rho;
area_fine=(ones(size(h))./(pm.*pn));
h_fine=h;
h_fine_new=h_fine;
%%

% loop on parent cells: [X,Y] is the cell polygon determined by psi (corner) points
% and replace parent bathy with avg child bathy.
for j=LLj+1:ULj
  for i=LLi+1:LRi
    X=[lon_psi_coarse(i-1,j-1) lon_psi_coarse(i,j-1) lon_psi_coarse(i,j) lon_psi_coarse(i-1,j)];
    Y=[lat_psi_coarse(i-1,j-1) lat_psi_coarse(i,j-1) lat_psi_coarse(i,j) lat_psi_coarse(i-1,j)];
    %find child cells within parent cell
    ind=find(inpolygon(lon_rho_fine(:),lat_rho_fine(:),X,Y));
    if ~isempty(ind),
      h_coarse_new(i,j)=sum(h_fine(ind).*area_fine(ind))./sum(area_fine(ind));
    end
  end
end
disp('finished averaging loop')
%
if (1)
% loop on the perimeter edges of the child grid and replace with parent
% bathy.
    %left strip
    for j=LLj:ULj+1
      for i=LLi:LLi %+1
        X=[lon_psi_coarse(i-1,j-1) lon_psi_coarse(i,j-1) lon_psi_coarse(i,j) lon_psi_coarse(i-1,j)];
        Y=[lat_psi_coarse(i-1,j-1) lat_psi_coarse(i,j-1) lat_psi_coarse(i,j) lat_psi_coarse(i-1,j)];
        %find child cells within parent cell
        ind=find(inpolygon(lon_rho_fine(:),lat_rho_fine(:),X,Y));
        if ~isempty(ind),
          h_fine_new(ind)=h_coarse_new(i,j);  % assign child depth as parent depth
        end
      end
    end
    disp('finished left strip')
    %right strip
    for j=LRj:URj+1
      for i=LRi+1:LRi+1
        X=[lon_psi_coarse(i-1,j-1) lon_psi_coarse(i,j-1) lon_psi_coarse(i,j) lon_psi_coarse(i-1,j)];
        Y=[lat_psi_coarse(i-1,j-1) lat_psi_coarse(i,j-1) lat_psi_coarse(i,j) lat_psi_coarse(i-1,j)];
        %find child cells within parent cell
        ind=find(inpolygon(lon_rho_fine(:),lat_rho_fine(:),X,Y));
        if ~isempty(ind),
          h_fine_new(ind)=h_coarse_new(i,j);  % assign child depth as parent depth
        end
      end
    end
    disp('finished right strip')
    %top strip
    for j=URj+1:URj+1
      for i=ULi:URi+1
        X=[lon_psi_coarse(i-1,j-1) lon_psi_coarse(i,j-1) lon_psi_coarse(i,j) lon_psi_coarse(i-1,j)];
        Y=[lat_psi_coarse(i-1,j-1) lat_psi_coarse(i,j-1) lat_psi_coarse(i,j) lat_psi_coarse(i-1,j)];
        %find child cells within parent cell
        ind=find(inpolygon(lon_rho_fine(:),lat_rho_fine(:),X,Y));
        if ~isempty(ind),
          h_fine_new(ind)=h_coarse_new(i,j);  % assign child depth as parent depth
        end
      end
    end
    disp('finished top strip')
    %bottom strip
    for j=LLj:LLj %+1
      for i=LLi:LRi+1
        X=[lon_psi_coarse(i-1,j-1) lon_psi_coarse(i,j-1) lon_psi_coarse(i,j) lon_psi_coarse(i-1,j)];
        Y=[lat_psi_coarse(i-1,j-1) lat_psi_coarse(i,j-1) lat_psi_coarse(i,j) lat_psi_coarse(i-1,j)];
        %find child cells within parent cell
        ind=find(inpolygon(lon_rho_fine(:),lat_rho_fine(:),X,Y));
        if ~isempty(ind),
          h_fine_new(ind)=h_coarse_new(i,j);  % assign child depth as parent depth
        end
      end
    end
    disp('finished bottom strip')
%   write new depths to child
    ncwrite(ncfile_fine,'h',h_fine_new);
end

%%
% write new depths to parent
ncwrite(ncfile_coarse,'h',h_coarse_new)

%%
% plot the new bathy
figure
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
 
% Bathymetry
  pcolorjw(xr,yr,h.*maskr)
  hold on
  shading flat
  colorbar
  if (mm==2)
    title('bathymetry - modified (m)')
    plot(xp(:,1),yp(:,1),'r'); plot(xp(:,end),yp(:,end),'r'); 
    plot(xp(1,:),yp(1,:),'r'); plot(xp(end,:),yp(end,:),'r'); 
  end
end
caxis([min(h(:)) max(h(:))])
set (gca, 'DataAspectRatio', [1 cos(mean(lat_rho_fine(:))*pi/180) 1] );

