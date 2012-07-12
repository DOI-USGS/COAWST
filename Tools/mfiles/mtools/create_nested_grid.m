% create_nested_grid.m
%
% user provides: - netcdf file name of coarse grid
%                - netcdf file name for new grid
%                - indices for location of child grid
%                - scale factor of increased resolution
%
% output is a netcdf roms grid and ascii swan grid with finer resolution
%
% User will need to modify the masking and bathymetry of new grid!!!!!
%
% Version jcwarner Aug 3, 2007
% Updated clea denamiel Feb 5, 2012
% Updated jcwarner 26Jun2012, make this a main driver.

%%%%%%%%%%%%% START OF USER SECTION %%%%%%%%%%%%%%%%

%1) ENTER NAME OF EXISITNG COARSE GRID FILE
%   The bathy in this file will be modified below, so you might want to
%   make a copy first.
%ncfile_coarse='refined_chan_grid.nc';
%ncfile_coarse='inlet_test_grid.nc';
ncfile_coarse='USeast_grd19_train.nc';

%2) ENTER NAME OF NEW FINE GRID FILE
%ncfile_fine='refined_chan_grid_ref5.nc';
%ncfile_fine='inlet_test_grid_ref5_test.nc';
ncfile_fine='Carolinas_grd6_train.nc';

%3) ENTER THE START AND END INDICES OF THE PSI POINTS OF THE COARSE GRID
%   THAT IDENTIFY THE OUTER BOUNDS OF THE FINE GRID
%   Mathworks netcdf is Fortran convention. Determine locations for child grid.
%   Select (Istr,Jstr) to locate the lower left corner of the new grid and
%   and    (Iend,Jend) as the upper right corner uisng
%   netcdf_load (ncfile_coarse)
%   figure; plot(lon_psi,lat_psi,'k'); hold on; plot(lon_psi', lat_psi','k')
%   plot(lon_psi(Istr,Jstr),lat_psi(Istr,Jstr),'r+')
%   plot(lon_psi(Iend,Jend),lat_psi(Iend,Jend),'r+')
% 
%Istr=30; Iend=50; Jstr=2; Jend=5;        % test_chan_refined
%Istr=24; Iend=54; Jstr=40; Jend=56;       % inlet_test_refined
Istr=411; Iend=549; Jstr=210; Jend=259;    % car_grd6 from useast

%4) ENTER SCALE FACTOR FOR INCREASED RESOLUTION (use 3 or 5)
scale=5;

%5) ENTER THE PROCEDURE THAT YOU WANT TO DO (1 = YES, 0 = NO)
create_chid_grid=1
merge_par_child_bathy=0
merge_par_child_mask=0


%%%%%%%%%%%%% END OF USER SECTION %%%%%%%%%%%%%%%%

%Call to create the child grid
if (create_chid_grid)
  parentchild_grid
end

%call to correct the bathy in both the parent and child
if (merge_par_child_bathy)
  parentchild_bathy
end

%call to correct masking
if (merge_par_child_mask)
  parentchild_mask
end


