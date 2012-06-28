%fix_parentchild_mask

% This script loads in parent and child masks and 
% modifies the child mask to be consistent at the boundaries 
% with the parent.

% User would still need to go in and make sure the mask transition is
% sensible.

% Istr, Jstr, Iend, Jstr are from create_nested_grid.

%specify corners
LLi=Istr;
LLj=Jstr;
URi=Iend;
URj=Jend;
ULi=Istr;
ULj=Jend;
LRi=Iend;
LRj=Jstr;

%load coordinates and bathymetry
netcdf_load(ncfile_coarse);
lon_psi_coarse=lon_psi;
lat_psi_coarse=lat_psi;
lon_rho_coarse=lon_rho;
lat_rho_coarse=lat_rho;
mask_coarse=mask_rho;

%%
netcdf_load(ncfile_fine)
lon_psi_fine=lon_psi;
lat_psi_fine=lat_psi;
lon_rho_fine=lon_rho;
lat_rho_fine=lat_rho;
mask_fine=mask_rho;
mask_fine_new=mask_fine;
%%

%left strip
for j=LLj:ULj+1
  for i=LLi:LLi+1
    X=[lon_psi_coarse(i-1,j-1) lon_psi_coarse(i,j-1) lon_psi_coarse(i,j) lon_psi_coarse(i-1,j)];
    Y=[lat_psi_coarse(i-1,j-1) lat_psi_coarse(i,j-1) lat_psi_coarse(i,j) lat_psi_coarse(i-1,j)];
    %find child cells within parent cell
    ind=find(inpolygon(lon_rho_fine(:),lat_rho_fine(:),X,Y));
    if ~isempty(ind),
      mask_fine_new(ind)=mask_coarse(i,j);  % assign child depth as parent depth
    end
  end
end
disp('finished left strip')
%right strip
for j=LRj:URj+1
  for i=LRi:LRi+1
    X=[lon_psi_coarse(i-1,j-1) lon_psi_coarse(i,j-1) lon_psi_coarse(i,j) lon_psi_coarse(i-1,j)];
    Y=[lat_psi_coarse(i-1,j-1) lat_psi_coarse(i,j-1) lat_psi_coarse(i,j) lat_psi_coarse(i-1,j)];
    %find child cells within parent cell
    ind=find(inpolygon(lon_rho_fine(:),lat_rho_fine(:),X,Y));
    if ~isempty(ind),
      mask_fine_new(ind)=mask_coarse(i,j);  % assign child depth as parent depth
    end
  end
end
disp('finished right strip')
%top strip
for j=URj:URj+1
  for i=ULi:URi+1
    X=[lon_psi_coarse(i-1,j-1) lon_psi_coarse(i,j-1) lon_psi_coarse(i,j) lon_psi_coarse(i-1,j)];
    Y=[lat_psi_coarse(i-1,j-1) lat_psi_coarse(i,j-1) lat_psi_coarse(i,j) lat_psi_coarse(i-1,j)];
    %find child cells within parent cell
    ind=find(inpolygon(lon_rho_fine(:),lat_rho_fine(:),X,Y));
    if ~isempty(ind),
      mask_fine_new(ind)=mask_coarse(i,j);  % assign child depth as parent depth
    end
  end
end
disp('finished top strip')
%bottom strip
for j=LLj:LLj+1
  for i=LLi:LRi+1
    X=[lon_psi_coarse(i-1,j-1) lon_psi_coarse(i,j-1) lon_psi_coarse(i,j) lon_psi_coarse(i-1,j)];
    Y=[lat_psi_coarse(i-1,j-1) lat_psi_coarse(i,j-1) lat_psi_coarse(i,j) lat_psi_coarse(i-1,j)];
    %find child cells within parent cell
    ind=find(inpolygon(lon_rho_fine(:),lat_rho_fine(:),X,Y));
    if ~isempty(ind),
      mask_fine_new(ind)=mask_coarse(i,j);  % assign child depth as parent depth
    end
  end
end
disp('finished bottom strip')


%%
% write new mask to child
u_mask = double(mask_fine_new(1:end-1,:) & mask_fine_new(2:end,:));
v_mask = double(mask_fine_new(:,1:end-1) & mask_fine_new(:,2:end));
p_mask = double( mask_fine_new(1:end-1,1:end-1) & mask_fine_new(1:end-1,2:end) & ...
                 mask_fine_new(2:end,1:end-1) & mask_fine_new(2:end,2:end));

ncwrite(ncfile_fine,'mask_rho',mask_fine_new);
ncwrite(ncfile_fine,'mask_u',u_mask);
ncwrite(ncfile_fine,'mask_v',v_mask);
ncwrite(ncfile_fine,'mask_psi',p_mask);




