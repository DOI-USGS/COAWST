      function resize_remap_vars(nmap, increment)

%      use grids                      ! module with grid information
       global grid1_size grid2_size grid1_rank grid2_rank
       global grid1_corners grid2_corners grid1_dims grid2_dims
       global grid1_name grid2_name grid1_units grid2_units grid1_mask grid2_mask
       global grid1_center_lat grid1_center_lon grid2_center_lat grid2_center_lon
       global grid1_area grid2_area grid1_area_in grid2_area_in grid1_frac grid2_frac
       global grid1_corner_lat grid1_corner_lon grid2_corner_lat grid2_corner_lon
       global luse_grid_centers luse_grid1_area luse_grid2_area grid1_bound_box grid2_bound_box
       global restrict_type num_srch_bins bin_addr1 bin_addr2 bin_lats bin_lons

%      use remap_vars                 ! common remapping variables
       global norm_opt_none norm_opt_dstarea norm_opt_frcarea
       global map_type_conserv map_type_bilinear map_type_bicubic map_type_distwgt
       global max_links_map1 num_links_map1 max_links_map2 num_links_map2
       global num_maps num_wts map_type norm_opt  resize_increment
       global grid1_add_map1 grid2_add_map1 grid1_add_map2 grid2_add_map2
       global wts_map1 wts_map2

     
%!-----------------------------------------------------------------------
%!
%!     this routine resizes remapping arrays by increasing(decreasing)
%!     the max_links by increment
%!
%!-----------------------------------------------------------------------

%!-----------------------------------------------------------------------
%!
%!     input variables
%!
%!-----------------------------------------------------------------------

%      integer (kind=int_kind), intent(in) ::
%     &     nmap,      ! identifies which mapping array to resize
%     &     increment  ! the number of links to add(subtract) to arrays

%!-----------------------------------------------------------------------
%!
%!     local variables
%!
%!-----------------------------------------------------------------------

%      integer (kind=int_kind) ::
%     &   ierr,     ! error flag
%     &   mxlinks   ! size of link arrays

%      integer (kind=int_kind), dimension(:), allocatable ::
%     &   add1_tmp, ! temp array for resizing address arrays
%     &   add2_tmp  ! temp array for resizing address arrays

%      real (kind=dbl_kind), dimension(:,:), allocatable ::
%     &   wts_tmp   ! temp array for resizing weight arrays

%!-----------------------------------------------------------------------
%!
%!     resize map 1 arrays if required.
%!
%!-----------------------------------------------------------------------

      switch (nmap)
      case(1)

%        !***
%        !*** allocate temporaries to hold original values
%        !***

        mxlinks = length(grid1_add_map1);
%        allocate (add1_tmp(mxlinks), add2_tmp(mxlinks), 
%     &            wts_tmp(num_wts,mxlinks))

        add1_tmp = grid1_add_map1;
        add2_tmp = grid2_add_map1;
        wts_tmp  = wts_map1;

%       !***
%       !*** deallocate originals and increment max_links then
%       !*** reallocate arrays at new size
%       !***

%       deallocate (grid1_add_map1, grid2_add_map1, wts_map1)
        clear global grid1_add_map1 grid2_add_map1 wts_map1

        max_links_map1 = mxlinks + increment;
%        allocate (grid1_add_map1(max_links_map1),
%     &            grid2_add_map1(max_links_map1),
%     &            wts_map1(num_wts,max_links_map1))
        global grid1_add_map1 grid2_add_map1 wts_map1

%       !***
%       !*** restore original values from temp arrays and
%       !*** deallocate temps
%       !***

        mxlinks = min(mxlinks, max_links_map1);
        grid1_add_map1(1:mxlinks) = add1_tmp (1:mxlinks);
        grid2_add_map1(1:mxlinks) = add2_tmp (1:mxlinks);
        wts_map1    (:,1:mxlinks) = wts_tmp(:,1:mxlinks);
        clear add1_tmp add2_tmp wts_tmp

%!-----------------------------------------------------------------------
%!
%!     resize map 2 arrays if required.
%!
%!-----------------------------------------------------------------------

      case(2)

%        !***
%        !*** allocate temporaries to hold original values
%        !***

        mxlinks = size(grid1_add_map2);
%        allocate (add1_tmp(mxlinks), add2_tmp(mxlinks), 
%     &            wts_tmp(num_wts,mxlinks),stat=ierr)
        if (ierr ~= 0)
          disp('error allocating temps in resize: ')
%          stop
        end

        add1_tmp = grid1_add_map2;
        add2_tmp = grid2_add_map2;
        wts_tmp  = wts_map2;
        
%        !***
%        !*** deallocate originals and increment max_links then
%        !*** reallocate arrays at new size
%        !***

        clear grid1_add_map2 grid2_add_map2 wts_map2
        max_links_map2 = mxlinks + increment;
%        allocate (grid1_add_map2(max_links_map2),
%     &            grid2_add_map2(max_links_map2),
%     &            wts_map2(num_wts,max_links_map2),stat=ierr)
        if (ierr ~= 0)
          disp('error allocating new arrays in resize: ')
        end


%        !***
%        !*** restore original values from temp arrays and
%        !*** deallocate temps
%        !***

        mxlinks = min(mxlinks, max_links_map2);
        grid1_add_map2(1:mxlinks) = add1_tmp (1:mxlinks);
        grid2_add_map2(1:mxlinks) = add2_tmp (1:mxlinks);
        wts_map2    (:,1:mxlinks) = wts_tmp(:,1:mxlinks);
        clear add1_tmp add2_tmp wts_tmp

      end

%!-----------------------------------------------------------------------

      end %subroutine resize_remap_vars

