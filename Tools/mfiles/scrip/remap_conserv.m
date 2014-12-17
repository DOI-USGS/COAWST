%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%!
%!     this module contains necessary routines for computing addresses
%!     and weights for a conservative interpolation  between any two 
%!     grids on a sphere.  the weights are computed by performing line 
%!     integrals around all overlap regions of the two grids.  see 
%!     Dukowicz and Kodis, SIAM J. Sci. Stat. Comput. 8, 305 (1987) and
%!     Jones, P.W. Monthly Weather Review (submitted).
%!
%!-----------------------------------------------------------------------
%!
%!     CVS:$Id: remap_conserv.f,v 1.10 2001/08/21 21:05:13 pwjones Exp $
%!
%!     Copyright (c) 1997, 1998 the Regents of the University of 
%!       California.
%!
%!     This software and ancillary information (herein called software) 
%!     called SCRIP is made available under the terms described here.  
%!     The software has been approved for release with associated 
%!     LA-CC Number 98-45.
%!
%!     Unless otherwise indicated, this software has been authored
%!     by an employee or employees of the University of California,
%!     operator of the Los Alamos National Laboratory under Contract
%!     No. W-7405-ENG-36 with the U.S. Department of Energy.  The U.S.
%!     Government has rights to use, reproduce, and distribute this
%!     software.  The public may copy and use this software without
%!     charge, provided that this Notice and any statement of authorship
%!     are reproduced on all copies.  Neither the Government nor the
%!     University makes any warranty, express or implied, or assumes
%!     any liability or responsibility for the use of this software.
%!
%!     If software is modified to produce derivative works, such modified
%!     software should be clearly marked, so as not to confuse it with 
%!     the version available from Los Alamos National Laboratory.
%!
%!***********************************************************************

%     module remap_conservative

%!-----------------------------------------------------------------------

 %     use kinds_mod    ! defines common data types
 %     use constants    ! defines common constants
 %     use timers       ! module for timing
 %     use grids        ! module containing grid information
 %     use remap_vars   ! module containing remap information

%      implicit none

%!-----------------------------------------------------------------------
%!
%!     module variables
%!
%!-----------------------------------------------------------------------

%       integer (kind=int_kind), save :: 
%      &        num_srch_cells ! num cells in restricted search arrays
% 
%       integer (kind=int_kind), dimension(:), allocatable, save :: 
%      &        srch_add       ! global address of cells in srch arrays
% 
%       real (kind=dbl_kind), parameter :: 
%      &     north_thresh = 1.45_dbl_kind, ! threshold for coord transf.
%      &     south_thresh =-2.00_dbl_kind  ! threshold for coord transf.
% 
%       real (kind=dbl_kind), dimension(:,:), allocatable, save ::
%      &     srch_corner_lat,  ! lat of each corner of srch cells
%      &     srch_corner_lon   ! lon of each corner of srch cells

%!***********************************************************************

%      contains

%!***********************************************************************

      function remap_conserv

%use grids                      ! module with grid information
       global grid1_size grid2_size grid1_rank grid2_rank
       global grid1_corners grid2_corners grid1_dims grid2_dims
       global grid1_name grid2_name grid1_units grid2_units grid1_mask grid2_mask
       global grid1_center_lat grid1_center_lon grid2_center_lat grid2_center_lon
       global grid1_area grid2_area grid1_area_in grid2_area_in grid1_frac grid2_frac
       global grid1_corner_lat grid1_corner_lon grid2_corner_lat grid2_corner_lon
       global luse_grid_centers luse_grid1_area luse_grid2_area grid1_bound_box grid2_bound_box
       global restrict_type num_srch_bins bin_addr1 bin_addr2 bin_lats bin_lons

%use remap_vars                 ! common remapping variables
       global norm_opt_none norm_opt_dstarea norm_opt_frcarea
       global map_type_conserv map_type_bilinear map_type_bicubic map_type_distwgt
       global max_links_map1 num_links_map1 max_links_map2 num_links_map2
       global num_maps num_wts map_type norm_opt  resize_increment
       global grid1_add_map1 grid2_add_map1 grid1_add_map2 grid2_add_map2
       global wts_map1 wts_map2

%use remap_conservative         ! routines for conservative remap
      global num_srch_cells srch_add north_thresh south_thresh srch_corner_lat srch_corner_lon

%intersection
      global last_loc lthresh intrsct_lat_off intrsct_lon_off
      global intrsct_lat intrsct_lon location

%pole intersection
      global luse_last intrsct_x intrsct_y avoid_pole_count avoid_pole_offset

%store_link_cnsrv
      global first_call link_add1 link_add2

% line integral
      global weights

      tiny=1e-14;
      north_thresh = 1.45; %! threshold for coord transf.
      south_thresh =-2.00; %! threshold for coord transf.
      lthresh = false;     %! flags segments crossing threshold bndy jcw
      luse_last = false;
      avoid_pole_count = 0;
      avoid_pole_offset = tiny;
      first_call=1;
      intrsct_lat=0;  % init
      intrsct_lon=0;  % init
      lcoinc=0;       % init
      weights=[1:6]*0;       % init
%!-----------------------------------------------------------------------
%!
%!     this routine traces the perimeters of every grid cell on each
%!     grid checking for intersections with the other grid and computing
%!     line integrals for each subsegment.
%!
%!-----------------------------------------------------------------------

%!-----------------------------------------------------------------------
%!
%!     local variables
%!
%!-----------------------------------------------------------------------

%      integer (kind=int_kind), parameter :: 
       max_subseg = 10000;%! max number of subsegments per segment
                          %! to prevent infinite loop
%      integer (kind=int_kind) :: 
%     &        grid1_add,  ! current linear address for grid1 cell
%     &        grid2_add,  ! current linear address for grid2 cell
%     &        min_add,    ! addresses for restricting search of
%     &        max_add,    !   destination grid
%     &        n, nwgt,    ! generic counters
%     &        corner,     ! corner of cell that segment starts from
%     &        next_corn,  ! corner of cell that segment ends on
%     &        num_subseg  ! number of subsegments 

%       logical (kind=log_kind) :: 
%      &        lcoinc,  ! flag for coincident segments
%      &        lrevers, ! flag for reversing direction of segment
%      &        lbegin   ! flag for first integration of a segment
% 
%       logical (kind=log_kind), dimension(:), allocatable ::
%      &        srch_mask   ! mask for restricting searches
% 
%       real (kind=dbl_kind) ::
%      &     intrsct_lat, intrsct_lon,       ! lat/lon of next intersect
%      &     beglat, endlat, beglon, endlon, ! endpoints of current seg.
%      &     norm_factor                     ! factor for normalizing wts
% 
%       real (kind=dbl_kind), dimension(:), allocatable ::
%      &       grid2_centroid_lat, grid2_centroid_lon, ! centroid coords
%      &       grid1_centroid_lat, grid1_centroid_lon  ! on each grid
% 
%       real (kind=dbl_kind), dimension(2) :: begseg ! begin lat/lon for
%                                                    ! full segment
% 
%       real (kind=dbl_kind), dimension(6) :: weights ! local wgt array
% 
% !-----------------------------------------------------------------------
% !
% !     initialize centroid arrays
% !
% !-----------------------------------------------------------------------
% 
%       allocate( grid1_centroid_lat(grid1_size),
%      &          grid1_centroid_lon(grid1_size),
%      &          grid2_centroid_lat(grid2_size),
%      &          grid2_centroid_lon(grid2_size))
% 
      grid1_centroid_lat = zeros(grid1_size,1);
      grid1_centroid_lon = zeros(grid1_size,1);
      grid2_centroid_lat = zeros(grid2_size,1);
      grid2_centroid_lon = zeros(grid2_size,1);

% !-----------------------------------------------------------------------
% !
% !     integrate around each cell on grid1
% !
% !-----------------------------------------------------------------------
% 
%     allocate(srch_mask(grid2_size))
      srch_mask=zeros(grid2_size,1);

      disp('grid1 sweep ')
      for grid1_add = 1:grid1_size

%        !***
%        !*** restrict searches first using search bins
%        !***

%       call timer_start(1)
        min_add = grid2_size;
        max_add = 1;
        for n=1:num_srch_bins
          if (grid1_add >= bin_addr1(1,n) & ...
              grid1_add <= bin_addr1(2,n))
            min_add = min(min_add, bin_addr2(1,n));
            max_add = max(max_add, bin_addr2(2,n));
          end
        end

%         !***
%         !*** further restrict searches using bounding boxes
%         !***

        num_srch_cells = 0;
        for grid2_add = min_add:max_add
          srch_mask(grid2_add) = (grid2_bound_box(1,grid2_add) <=  ...
                                 grid1_bound_box(2,grid1_add)) &  ...
                                (grid2_bound_box(2,grid2_add) >=  ...
                                 grid1_bound_box(1,grid1_add)) &  ...
                                (grid2_bound_box(3,grid2_add) <=  ...
                                 grid1_bound_box(4,grid1_add)) &  ...
                                (grid2_bound_box(4,grid2_add) >=  ...
                                 grid1_bound_box(3,grid1_add));

          if (srch_mask(grid2_add)); num_srch_cells = num_srch_cells+1; end
        end

%         !***
%         !*** create search arrays
%         !***

%        allocate(srch_add(num_srch_cells),
%     &           srch_corner_lat(grid2_corners,num_srch_cells),
%     &           srch_corner_lon(grid2_corners,num_srch_cells))
        global srch_corner_lat srch_corner_lon srch_add
        srch_add=zeros(num_srch_cells,1);
        srch_corner_lat=zeros(grid2_corners,num_srch_cells);
        srch_corner_lon=zeros(grid2_corners,num_srch_cells);

        n = 0;
        for grid2_add = min_add:max_add
          if (srch_mask(grid2_add))
            n = n+1;
            srch_add(n) = grid2_add;
            srch_corner_lat(:,n) = grid2_corner_lat(:,grid2_add);
            srch_corner_lon(:,n) = grid2_corner_lon(:,grid2_add);
          end
        end
%        call timer_stop(1)
% 
%         !***
%         !*** integrate around this cell
%         !***

        for corner = 1:grid1_corners
          next_corn = mod(corner,grid1_corners) + 1;

%           !***
%           !*** define endpoints of the current segment
%           !***

          beglat = grid1_corner_lat(corner,grid1_add);
          beglon = grid1_corner_lon(corner,grid1_add);
          endlat = grid1_corner_lat(next_corn,grid1_add);
          endlon = grid1_corner_lon(next_corn,grid1_add);
          lrevers = false;

%           !***
%           !*** to ensure exact path taken during both
%           !*** sweeps, always integrate segments in the same 
%           !*** direction (SW to NE).
%           !***

          if ((endlat < beglat) |  ...
             (endlat == beglat & endlon < beglon))
            beglat = grid1_corner_lat(next_corn,grid1_add);
            beglon = grid1_corner_lon(next_corn,grid1_add);
            endlat = grid1_corner_lat(corner,grid1_add);
            endlon = grid1_corner_lon(corner,grid1_add);
            lrevers = true;
          end

          begseg(1) = beglat;
          begseg(2) = beglon;
          lbegin = true;
          num_subseg = 0;

%           !***
%           !*** if this is a constant-longitude segment, skip the rest 
%           !*** since the line integral contribution will be zero.
%           !***

          if (endlon ~= beglon)

%           !***
%           !*** integrate along this segment, detecting intersections 
%           !*** and computing the line integral for each sub-segment
%           !***

            while (beglat ~= endlat | beglon ~= endlon)

%             !***
%             !*** prevent infinite loops if integration gets stuck
%             !*** near cell or threshold boundary
%             !***

            num_subseg = num_subseg + 1;
            if (num_subseg > max_subseg)
              disp(['STOPPED  max subseg exceeded'])
              return %stop 'integration stalled: num_subseg exceeded limit'
            end

%             !***
%             !*** find next intersection of this segment with a grid
%             !*** line on grid 2.
%             !***

%           call timer_start(2)
            location=grid2_add;
            intersection(location,intrsct_lat,intrsct_lon,lcoinc, ...
                         beglat, beglon, endlat, endlon, begseg,  ...
                         lbegin, lrevers)
            grid2_add=location;

%            call timer_stop(2)
            lbegin = false;

%             !***
%             !*** compute line integral for this subsegment.
%             !***

%            call timer_start(3)
            if (grid2_add ~= 0)
              line_integral(weights, num_wts, ...
                              beglon, intrsct_lon, beglat, intrsct_lat, ...
                              grid1_center_lat(grid1_add),  ...
                              grid1_center_lon(grid1_add),  ...
                              grid2_center_lat(grid2_add),  ...
                              grid2_center_lon(grid2_add));
            else
              line_integral(weights, num_wts, ...
                              beglon, intrsct_lon, beglat, intrsct_lat, ...
                              grid1_center_lat(grid1_add),  ...
                              grid1_center_lon(grid1_add),  ...
                              grid1_center_lat(grid1_add),  ...
                              grid1_center_lon(grid1_add));
           end

%            call timer_stop(3)

%             !***
%             !*** if integrating in reverse order, change
%             !*** sign of weights
%             !***

            if (lrevers)
              weights = -weights;
            end

%             !***
%             !*** store the appropriate addresses and weights. 
%             !*** also add contributions to cell areas and centroids.
%             !***
% 
%             !if (grid1_add == 119247) then
%             !  print *,grid1_add,grid2_add,corner,weights(1)
%             !  print *,grid1_corner_lat(:,grid1_add)
%             !  print *,grid1_corner_lon(:,grid1_add)
%             !  print *,grid2_corner_lat(:,grid2_add)
%             !  print *,grid2_corner_lon(:,grid2_add)
%             !  print *,beglat,beglon,intrsct_lat,intrsct_lon
%             !endif

            if (grid2_add ~= 0)
              if (grid1_mask(grid1_add))
%               call timer_start(4)
                add1=grid1_add; add2=grid2_add;
                store_link_cnsrv(add1, add2, weights)
%               call timer_stop(4)
                grid1_frac(grid1_add) = grid1_frac(grid1_add) + ...
                                       weights(1);
                grid2_frac(grid2_add) = grid2_frac(grid2_add) + ...
                                       weights(num_wts+1);
              end
            end

            grid1_area(grid1_add) = grid1_area(grid1_add) + weights(1);
            grid1_centroid_lat(grid1_add) = grid1_centroid_lat(grid1_add) + weights(2);
            grid1_centroid_lon(grid1_add) = grid1_centroid_lon(grid1_add) + weights(3);

%             !***
%             !*** reset beglat and beglon for next subsegment.
%             !***

            beglat = intrsct_lat;
            beglon = intrsct_lon;
            end

          end

%           !***
%           !*** end of segment
%           !***

        end

%         !***
%         !*** finished with this cell: deallocate search array and
%         !*** start on next cell

        clear global srch_add srch_corner_lat srch_corner_lon

      end

      clear global srch_mask

% !-----------------------------------------------------------------------
% !
% !     integrate around each cell on grid2
% !
% !-----------------------------------------------------------------------

      global srch_mask
      srch_mask=zeros(grid1_size,1);

      disp('grid2 sweep ')
      for grid2_add = 1:grid2_size

%         !***
%         !*** restrict searches first using search bins
%         !***

%       call timer_start(5)
        min_add = grid1_size;
        max_add = 1;
        for n=1:num_srch_bins
          if (grid2_add >= bin_addr2(1,n) & ...
             grid2_add <= bin_addr2(2,n))
            min_add = min(min_add, bin_addr1(1,n));
            max_add = max(max_add, bin_addr1(2,n));
          end
        end

%         !***
%         !*** further restrict searches using bounding boxes
%         !***

        num_srch_cells = 0;
        for grid1_add = min_add:max_add
          srch_mask(grid1_add) = (grid1_bound_box(1,grid1_add) <= ...
                                 grid2_bound_box(2,grid2_add)) & ...
                                (grid1_bound_box(2,grid1_add) >=  ...
                                 grid2_bound_box(1,grid2_add)) & ...
                                (grid1_bound_box(3,grid1_add) <= ...
                                 grid2_bound_box(4,grid2_add)) & ...
                                (grid1_bound_box(4,grid1_add) >= ...
                                 grid2_bound_box(3,grid2_add));

          if (srch_mask(grid1_add)); num_srch_cells = num_srch_cells+1; end
        end

        global srch_corner_lat srch_corner_lon srch_add
        srch_add=zeros(num_srch_cells,1);
        srch_corner_lat=zeros(grid1_corners,num_srch_cells);
        srch_corner_lon=zeros(grid1_corners,num_srch_cells);

        n = 0;
        for grid1_add = min_add:max_add
          if (srch_mask(grid1_add))
            n = n+1;
            srch_add(n) = grid1_add;
            srch_corner_lat(:,n) = grid1_corner_lat(:,grid1_add);
            srch_corner_lon(:,n) = grid1_corner_lon(:,grid1_add);
          end
        end
%        call timer_stop(5)

%         !***
%         !*** integrate around this cell
%         !***

        for corner = 1:grid2_corners
          next_corn = mod(corner,grid2_corners) + 1;

          beglat = grid2_corner_lat(corner,grid2_add);
          beglon = grid2_corner_lon(corner,grid2_add);
          endlat = grid2_corner_lat(next_corn,grid2_add);
          endlon = grid2_corner_lon(next_corn,grid2_add);
          lrevers = false;

%           !***
%           !*** to ensure exact path taken during both
%           !*** sweeps, always integrate in the same direction
%           !***

          if ((endlat < beglat) | ...
              (endlat == beglat & endlon < beglon))
            beglat = grid2_corner_lat(next_corn,grid2_add);
            beglon = grid2_corner_lon(next_corn,grid2_add);
            endlat = grid2_corner_lat(corner,grid2_add);
            endlon = grid2_corner_lon(corner,grid2_add);
            lrevers = true;
          end

          begseg(1) = beglat;
          begseg(2) = beglon;
          lbegin = true;

%           !***
%           !*** if this is a constant-longitude segment, skip the rest 
%           !*** since the line integral contribution will be zero.
%           !***

          if (endlon ~= beglon)
            num_subseg = 0;

%           !***
%           !*** integrate along this segment, detecting intersections 
%           !*** and computing the line integral for each sub-segment
%           !***

            while (beglat ~= endlat || beglon ~= endlon)

%             !***
%             !*** prevent infinite loops if integration gets stuck
%             !*** near cell or threshold boundary
%             !***

              num_subseg = num_subseg + 1;
              if (num_subseg > max_subseg)
                disp('integration stalled: num_subseg exceeded limit');
                return;
              end

%             !***
%             !*** find next intersection of this segment with a line 
%             !*** on grid 2.
%             !***

%             call timer_start(6)
              location=grid1_add;
              intersection(location,intrsct_lat,intrsct_lon,lcoinc, ...
                               beglat, beglon, endlat, endlon, begseg, ...
                               lbegin, lrevers)
              grid1_add=location;
%              call timer_stop(6)
               lbegin = false;
 
%               !***
%               !*** compute line integral for this subsegment.
%               !***

%              call timer_start(7)
              if (grid1_add ~= 0)
                line_integral(weights, num_wts,   ...
                                beglon, intrsct_lon, beglat, intrsct_lat, ...
                                grid1_center_lat(grid1_add), ...
                                grid1_center_lon(grid1_add), ...
                                grid2_center_lat(grid2_add), ...
                                grid2_center_lon(grid2_add))
              else
                line_integral(weights, num_wts, ...
                                beglon, intrsct_lon, beglat, intrsct_lat, ...
                                grid2_center_lat(grid2_add), ...
                                grid2_center_lon(grid2_add), ...
                                grid2_center_lat(grid2_add), ...
                                grid2_center_lon(grid2_add))
              end
%             call timer_stop(7)

              if (lrevers)
                weights = -weights;
              end

%             !***
%             !*** store the appropriate addresses and weights. 
%             !*** also add contributions to cell areas and centroids.
%             !*** if there is a coincidence, do not store weights
%             !*** because they have been captured in the previous loop.
%             !*** the grid1 mask is the master mask
%             !***
% 
%             !if (grid1_add == 119247) then
%             !  print *,grid1_add,grid2_add,corner,weights(1)
%             !  print *,grid1_corner_lat(:,grid1_add)
%             !  print *,grid1_corner_lon(:,grid1_add)
%             !  print *,grid2_corner_lat(:,grid2_add)
%             !  print *,grid2_corner_lon(:,grid2_add)
%             !  print *,beglat,beglon,intrsct_lat,intrsct_lon
%             !endif

              if (~lcoinc & grid1_add ~= 0)
                if (grid1_mask(grid1_add))
%                 call timer_start(8)
                  add1=grid1_add; add2=grid2_add;
                  store_link_cnsrv(add1, add2, weights)
%                 call timer_stop(8)
                  grid1_frac(grid1_add) = grid1_frac(grid1_add) + ...
                                         weights(1);
                  grid2_frac(grid2_add) = grid2_frac(grid2_add) + ...
                                         weights(num_wts+1);
                end
              end

              grid2_area(grid2_add) = grid2_area(grid2_add) + ...
                                             weights(num_wts+1);
              grid2_centroid_lat(grid2_add) = ...
              grid2_centroid_lat(grid2_add) + weights(num_wts+2);
              grid2_centroid_lon(grid2_add) = ...
              grid2_centroid_lon(grid2_add) + weights(num_wts+3);

%             !***
%             !*** reset beglat and beglon for next subsegment.
%             !***

              beglat = intrsct_lat;
              beglon = intrsct_lon;
            end

          end

%           !***
%           !*** end of segment
%           !***

        end

%         !***
%         !*** finished with this cell: deallocate search array and
%         !*** start on next cell

        clear global srch_add srch_corner_lat srch_corner_lon

      end

      clear global srch_mask

% !-----------------------------------------------------------------------
% !
% !     correct for situations where N/S pole not explicitly included in
% !     grid (i.e. as a grid corner point). if pole is missing from only
% !     one grid, need to correct only the area and centroid of that 
% !     grid.  if missing from both, do complete weight calculation.
% !
% !-----------------------------------------------------------------------

%     !*** North Pole
      weights(1) =  2*pi;
      weights(2) =  pi*pi;
      weights(3) =  0;
      weights(4) =  2*pi;
      weights(5) =  pi*pi;
      weights(6) =  0;

      grid1_add = 0;
      for n=1:grid1_size
        if (grid1_area(n) < -3*pi/2 & ...
            grid1_center_lat(n) > 0)
          grid1_add = n;
          break; %      exit pole_loop1
        end
      end

      grid2_add = 0;
      for n=1:grid2_size
        if (grid2_area(n) < -3*pi/2 & ...
           grid2_center_lat(n) > 0)
          grid2_add = n;
          break; %exit pole_loop2
        end
      end

      if (grid1_add ~=0)
        grid1_area(grid1_add) = grid1_area(grid1_add) + weights(1);
        grid1_centroid_lat(grid1_add) = ...
        grid1_centroid_lat(grid1_add) + weights(2);
        grid1_centroid_lon(grid1_add) = ...
        grid1_centroid_lon(grid1_add) + weights(3);
      end

      if (grid2_add ~=0)
        grid2_area(grid2_add) = grid2_area(grid2_add) + ...
                                        weights(num_wts+1);
        grid2_centroid_lat(grid2_add) = ...
        grid2_centroid_lat(grid2_add) + weights(num_wts+2);
        grid2_centroid_lon(grid2_add) = ...
        grid2_centroid_lon(grid2_add) + weights(num_wts+3);
      end

      if (grid1_add ~= 0 & grid2_add ~=0)
        add1=grid1_add; add2=grid2_add;
        store_link_cnsrv(add1, add2, weights)
        grid1_frac(grid1_add) = grid1_frac(grid1_add) + ...
                                weights(1);
        grid2_frac(grid2_add) = grid2_frac(grid2_add) + ...
                                weights(num_wts+1);
      end

%     !*** South Pole
      weights(1) =  2*pi;
      weights(2) = -pi*pi;
      weights(3) =  0;
      weights(4) =  2*pi;
      weights(5) = -pi*pi;
      weights(6) =  0;

      grid1_add = 0;
      for n=1:grid1_size
        if (grid1_area(n) < -3*pi/2 & ...
            grid1_center_lat(n) < 0)
          grid1_add = n;
          break; %exit pole_loop3
        end
      end

      grid2_add = 0;
      for n=1:grid2_size
        if (grid2_area(n) < -3*pi/2 & ...
            grid2_center_lat(n) < 0)
          grid2_add = n;
          break; %exit pole_loop4
        end
      end
 
      if (grid1_add ~=0)
        grid1_area(grid1_add) = grid1_area(grid1_add) + weights(1);
        grid1_centroid_lat(grid1_add) = ...
        grid1_centroid_lat(grid1_add) + weights(2);
        grid1_centroid_lon(grid1_add) = ...
        grid1_centroid_lon(grid1_add) + weights(3);
      end

      if (grid2_add ~=0)
        grid2_area(grid2_add) = grid2_area(grid2_add) + ...
                                        weights(num_wts+1);
        grid2_centroid_lat(grid2_add) = ...
        grid2_centroid_lat(grid2_add) + weights(num_wts+2);
        grid2_centroid_lon(grid2_add) = ...
        grid2_centroid_lon(grid2_add) + weights(num_wts+3);
      end

      if (grid1_add ~= 0 & grid2_add ~=0)
        add1=grid1_add; add2=grid2_add;
        store_link_cnsrv(add1, add2, weights)

        grid1_frac(grid1_add) = grid1_frac(grid1_add) + ...
                                weights(1);
        grid2_frac(grid2_add) = grid2_frac(grid2_add) + ...
                                weights(num_wts+1);
      end

% !-----------------------------------------------------------------------
% !
% !     finish centroid computation
% !
% !-----------------------------------------------------------------------

%      where (grid1_area ~= zero)
%        grid1_centroid_lat = grid1_centroid_lat/grid1_area
%        grid1_centroid_lon = grid1_centroid_lon/grid1_area
%      end where
%      where (grid2_area /= zero)
%        grid2_centroid_lat = grid2_centroid_lat/grid2_area
%        grid2_centroid_lon = grid2_centroid_lon/grid2_area
%      end where
      zz=grid1_area~=0;
      grid1_centroid_lat(zz) = grid1_centroid_lat(zz)./grid1_area(zz);
      grid1_centroid_lon(zz) = grid1_centroid_lon(zz)./grid1_area(zz);
      zz=grid2_area~=0;
      grid2_centroid_lat(zz) = grid2_centroid_lat(zz)./grid2_area(zz);
      grid2_centroid_lon(zz) = grid2_centroid_lon(zz)./grid2_area(zz);

% !-----------------------------------------------------------------------
% !
% !     include centroids in weights and normalize using destination
% !     area if requested
% !
% !-----------------------------------------------------------------------

      for n=1:num_links_map1
        grid1_add = grid1_add_map1(n);
        grid2_add = grid2_add_map1(n);
        for nwgt=1:num_wts
          weights(        nwgt) = wts_map1(nwgt,n);
          if (num_maps > 1)
            weights(num_wts+nwgt) = wts_map2(nwgt,n);
          end
        end

        switch (norm_opt)
        case ('norm_opt_dstarea')
          if (grid2_area(grid2_add) ~= 0)
            if (luse_grid2_area)
              norm_factor = one/grid2_area_in(grid2_add);
            else
              norm_factor = one/grid2_area(grid2_add);
            end
          else
            norm_factor = 0;
          end
        case ('norm_opt_frcarea')
          if (grid2_frac(grid2_add) ~= 0)
            if (luse_grid2_area)
              norm_factor = grid2_area(grid2_add)/ ...
                           (grid2_frac(grid2_add)* ...
                            grid2_area_in(grid2_add));
            else
              norm_factor = 1/grid2_frac(grid2_add);
            end
          else
            norm_factor = 0;
          end
        case ('norm_opt_none')
          norm_factor = 1;
        end

        wts_map1(1,n) =  weights(1)*norm_factor;
        wts_map1(2,n) = (weights(2) - weights(1)* ...
                                    grid1_centroid_lat(grid1_add))* ...
                                    norm_factor;
        wts_map1(3,n) = (weights(3) - weights(1)* ...
                                    grid1_centroid_lon(grid1_add))* ...
                                    norm_factor;

        if (num_maps > 1)
          switch (norm_opt)
          case ('norm_opt_dstarea')
            if (grid1_area(grid1_add) ~= 0)
              if (luse_grid1_area)
                norm_factor = 1/grid1_area_in(grid1_add);
              else
                norm_factor = 1/grid1_area(grid1_add);
              end
            else
              norm_factor = 0;
            end
          case ('norm_opt_frcarea')
            if (grid1_frac(grid1_add) ~= zero)
              if (luse_grid1_area)
                norm_factor = grid1_area(grid1_add)/ ...
                             (grid1_frac(grid1_add)* ...
                              grid1_area_in(grid1_add));
              else
                norm_factor = 1/grid1_frac(grid1_add);
              end
            else
              norm_factor = 0;
            end
          case ('norm_opt_none')
            norm_factor = 1;
          end

          wts_map2(1,n) =  weights(num_wts+1)*norm_factor;
          wts_map2(2,n) = (weights(num_wts+2) - weights(num_wts+1)* ...
                                      grid2_centroid_lat(grid2_add))* ...
                                      norm_factor;
          wts_map2(3,n) = (weights(num_wts+3) - weights(num_wts+1)* ...
                                      grid2_centroid_lon(grid2_add))* ...
                                      norm_factor;
        end

      end

      disp(['Total number of links = ',num2str(num_links_map1)])

%     where (grid1_area /= zero) grid1_frac = grid1_frac/grid1_area
%     where (grid2_area /= zero) grid2_frac = grid2_frac/grid2_area
      zz=grid1_area~=0;  grid1_frac(zz) = grid1_frac(zz)./grid1_area(zz);
      zz=grid2_area~=0;  grid2_frac(zz) = grid2_frac(zz)./grid2_area(zz);

% !-----------------------------------------------------------------------
% !
% !     perform some error checking on final weights
% !
% !-----------------------------------------------------------------------

%     grid2_centroid_lat = 0;
%     grid2_centroid_lon = 0;
      grid2_centroid_lat = zeros(grid2_size,1);
      grid2_centroid_lon = zeros(grid2_size,1);

      for n=1:grid1_size
        if (grid1_area(n) < -.01)
          disp(['Grid 1 area error: ',num2str(n),grid1_area(n)])
        end
        if (grid1_centroid_lat(n) < -pi/2-.01 | ...
            grid1_centroid_lat(n) >  pi/2+.01)
          disp(['Grid 1 centroid lat error: ',num2str(n),grid1_centroid_lat(n)])
        end
        grid1_centroid_lat(n) = 0;
        grid1_centroid_lon(n) = 0;
      end

      for n=1:grid2_size
        if (grid2_area(n) < -.01)
          disp(['Grid 2 area error: ',num2str(n),grid2_area(n)])
        end
        if (grid2_centroid_lat(n) < -pi/2-.01 | ...
            grid2_centroid_lat(n) >  pi/2+.01)
          disp(['Grid 2 centroid lat error: ',num2str(n),grid2_centroid_lat(n)])
        end
        grid2_centroid_lat(n) = 0;
        grid2_centroid_lon(n) = 0;
      end

      for n=1:num_links_map1
        grid1_add = grid1_add_map1(n);
        grid2_add = grid2_add_map1(n);
        
        if (wts_map1(1,n) < -.01)
          disp(['Map 1 weight < 0 ',num2str(grid1_add),num2str(grid2_add),num2str(wts_map1(1,n))])
        end
        if (norm_opt ~= norm_opt_none & wts_map1(1,n) > 1.01)
          disp(['Map 1 weight > 1 ',num2str(grid1_add),num2str(grid2_add),num2str(wts_map1(1,n))])
        end
        grid2_centroid_lat(grid2_add) = ...
        grid2_centroid_lat(grid2_add) + wts_map1(1,n);

        if (num_maps > 1)
          if (wts_map2(1,n) < -.01)
            disp(['Map 2 weight < 0 ',num2str(grid1_add),num2str(grid2_add), ...
                                      num2str(wts_map2(1,n))])
          end
          if (norm_opt ~= norm_opt_none & wts_map2(1,n) > 1.01)
            disp(['Map 2 weight < 0 ',num2str(grid1_add),num2str(grid2_add), ...
                                      num2str(wts_map2(1,n))])
          end
          grid1_centroid_lat(grid1_add) =  ...
          grid1_centroid_lat(grid1_add) + wts_map2(1,n);
        end
      end

      for n=1:grid2_size
        switch (norm_opt)
        case ('norm_opt_dstarea')
          norm_factor = grid2_frac(grid2_add);
        case ('norm_opt_frcarea')
          norm_factor = 1;
        case ('norm_opt_none')
          if (luse_grid2_area)
            norm_factor = grid2_area_in(grid2_add);
          else
            norm_factor = grid2_area(grid2_add);
          end
        end
        if (abs(grid2_centroid_lat(grid2_add)-norm_factor) > .01)
          disp(['Error: sum of wts for map1 ',num2str(grid2_add), ...
                  num2str(grid2_centroid_lat(grid2_add)),num2str(norm_factor)])
        end
      end

      if (num_maps > 1)
        for n=1:grid1_size
          switch (norm_opt)
          case ('norm_opt_dstarea')
            norm_factor = grid1_frac(grid1_add);
          case ('norm_opt_frcarea')
            norm_factor = 1;
          case ('norm_opt_none')
            if (luse_grid1_area)
              norm_factor = grid1_area_in(grid1_add);
            else
              norm_factor = grid1_area(grid1_add);
            end
          end
          if (abs(grid1_centroid_lat(grid1_add)-norm_factor) > .01)
            disp(['Error: sum of wts for map2 ',num2str(grid1_add), ...
                  num2str(grid1_centroid_lat(grid1_add)),num2str(norm_factor)])
          end
        end
      end
%!-----------------------------------------------------------------------

      end %subroutine remap_conserv

%!***********************************************************************

      function intersection(location,intrsct_lat,intrsct_lon,lcoinc, ...
                            beglat, beglon, endlat, endlon, begseg, ...
                            lbegin, lrevers)

% !-----------------------------------------------------------------------
% !
% !     this routine finds the next intersection of a destination grid 
% !     line with the line segment given by beglon, endlon, etc.
% !     a coincidence flag is returned if the segment is entirely 
% !     coincident with an ocean grid line.  the cells in which to search
% !     for an intersection must have already been restricted in the
% !     calling routine.
% !
% !-----------------------------------------------------------------------
% 
% !-----------------------------------------------------------------------
% !
% !     intent(in): 
% !
% !-----------------------------------------------------------------------

%       logical (kind=log_kind), intent(in) ::
%      &     lbegin, ! flag for first integration along this segment
%      &     lrevers ! flag whether segment integrated in reverse
% 
%       real (kind=dbl_kind), intent(in) :: 
%      &     beglat, beglon,  ! beginning lat/lon endpoints for segment
%      &     endlat, endlon   ! ending    lat/lon endpoints for segment
% 
%       real (kind=dbl_kind), dimension(2), intent(inout) :: 
%      &     begseg ! begin lat/lon of full segment
% 
% !-----------------------------------------------------------------------
% !
% !     intent(out): 
% !
% !-----------------------------------------------------------------------
% 
%       integer (kind=int_kind), intent(out) ::
%      &        location  ! address in destination array containing this
%                         ! segment
% 
%       logical (kind=log_kind), intent(out) ::
%      &        lcoinc    ! flag segments which are entirely coincident
%                         ! with a grid line
% 
%       real (kind=dbl_kind), intent(out) ::
%      &     intrsct_lat, intrsct_lon ! lat/lon coords of next intersect.
% 
% !-----------------------------------------------------------------------
% !
% !     local variables
% !
% !-----------------------------------------------------------------------

%use grids                      ! module with grid information
       global grid1_size grid2_size grid1_rank grid2_rank
       global grid1_corners grid2_corners grid1_dims grid2_dims
       global grid1_name grid2_name grid1_units grid2_units grid1_mask grid2_mask
       global grid1_center_lat grid1_center_lon grid2_center_lat grid2_center_lon
       global grid1_area grid2_area grid1_area_in grid2_area_in grid1_frac grid2_frac
       global grid1_corner_lat grid1_corner_lon grid2_corner_lat grid2_corner_lon
       global luse_grid_centers luse_grid1_area luse_grid2_area grid1_bound_box grid2_bound_box
       global restrict_type num_srch_bins bin_addr1 bin_addr2 bin_lats bin_lons

%use remap_vars                 ! common remapping variables
       global norm_opt_none norm_opt_dstarea norm_opt_frcarea
       global map_type_conserv map_type_bilinear map_type_bicubic map_type_distwgt
       global max_links_map1 num_links_map1 max_links_map2 num_links_map2
       global num_maps num_wts map_type norm_opt  resize_increment
       global grid1_add_map1 grid2_add_map1 grid1_add_map2 grid2_add_map2
       global wts_map1 wts_map2

%use remap_conservative         ! routines for conservative remap
      global num_srch_cells srch_add north_thresh south_thresh srch_corner_lat srch_corner_lon

%intersection
      global last_loc lthresh intrsct_lat_off intrsct_lon_off
      global intrsct_lat intrsct_lon location

%pole intersection
      global luse_last intrsct_x intrsct_y avoid_pole_count avoid_pole_offset

%store_link_cnsrv
      global first_call

%       integer (kind=int_kind) :: n, next_n, cell, srch_corners, pole_loc
% 
%       integer (kind=int_kind), save :: 
%      &     last_loc  ! save location when crossing threshold
% 
%       logical (kind=log_kind) :: 
%      &     loutside  ! flags points outside grid
% 
%       logical (kind=log_kind), save :: 
%      &     lthresh = .false.  ! flags segments crossing threshold bndy
% 
%       real (kind=dbl_kind) ::
%      &     lon1, lon2,       ! local longitude variables for segment
%      &     lat1, lat2,       ! local latitude  variables for segment
%      &     grdlon1, grdlon2, ! local longitude variables for grid cell
%      &     grdlat1, grdlat2, ! local latitude  variables for grid cell
%      &     vec1_lat, vec1_lon, ! vectors and cross products used
%      &     vec2_lat, vec2_lon, ! during grid search
%      &     cross_product, 
%      &     eps, offset,        ! small offset away from intersect
%      &     s1, s2, determ,     ! variables used for linear solve to
%      &     mat1, mat2, mat3, mat4, rhs1, rhs2  ! find intersection
% 
%       real (kind=dbl_kind), save ::
%      &     intrsct_lat_off, intrsct_lon_off ! lat/lon coords offset 
%                                             ! for next search
% 
% !-----------------------------------------------------------------------
% !
% !     initialize defaults, flags, etc.
% !
% !-----------------------------------------------------------------------

      location = 0;
      lcoinc = false;
      intrsct_lat = endlat;
      intrsct_lon = endlon;
      tiny=1e-14;
      eps=1e-14;


      if (num_srch_cells == 0); return; end

      if (beglat > north_thresh || beglat < south_thresh)

        if (lthresh); location = last_loc; end
        pole_intersection(location, ...
                    intrsct_lat,intrsct_lon,lcoinc,lthresh, ...
                    beglat, beglon, endlat, endlon, begseg, lrevers)
        if (lthresh)
          last_loc = location;
          intrsct_lat_off = intrsct_lat;
          intrsct_lon_off = intrsct_lon;
        end
        return

      end

      loutside = false;
      if (lbegin)
        lat1 = beglat;
        lon1 = beglon;
      else
        lat1 = intrsct_lat_off;
        lon1 = intrsct_lon_off;
      end
      lat2 = endlat;
      lon2 = endlon;
      if ((lon2-lon1) > 3*pi/2)
        lon2 = lon2 - pi2;
      elseif ((lon2-lon1) < -3*pi/2)
        lon2 = lon2 + pi2;
      end
      s1 = 0;

% !-----------------------------------------------------------------------
% !
% !     search for location of this segment in ocean grid using cross
% !     product method to determine whether a point is enclosed by a cell
% !
% !-----------------------------------------------------------------------

%     call timer_start(12)
      srch_corners = size(srch_corner_lat,1);
%     srch_loop: do
      srch_loop=1;
      while (srch_loop)

%         !***
%         !*** if last segment crossed threshold, use that location
%         !***
% 
        if (lthresh)
          for cell=1:num_srch_cells
            if (srch_add(cell) == last_loc)
              location = last_loc;
              eps = tiny;
              srch_loop=0;
              break %exit srch_loop
            end
          end
        end

%         !***
%         !*** otherwise normal search algorithm
%         !***

        for cell=1:num_srch_cells
          src_corners=1;
          while (src_corners)
            jcw_count=0;
            for n=1:srch_corners
              next_n = mod(n,srch_corners) + 1;

%             !***
%             !*** here we take the cross product of the vector making 
%             !*** up each cell side with the vector formed by the vertex
%             !*** and search point.  if all the cross products are 
%             !*** positive, the point is contained in the cell.
%             !***

              vec1_lat = srch_corner_lat(next_n,cell) - ...
                         srch_corner_lat(n     ,cell);
              vec1_lon = srch_corner_lon(next_n,cell) - ...
                         srch_corner_lon(n     ,cell);
              vec2_lat = lat1 - srch_corner_lat(n,cell);
              vec2_lon = lon1 - srch_corner_lon(n,cell);

%             !***
%             !*** if endpoint coincident with vertex, offset
%             !*** the endpoint
%             !***

              if ((vec2_lat == 0) && (vec2_lon == 0))
                lat1 = lat1 + 1.d-10*(lat2-lat1);
                lon1 = lon1 + 1.d-10*(lon2-lon1);
                vec2_lat = lat1 - srch_corner_lat(n,cell);
                vec2_lon = lon1 - srch_corner_lon(n,cell);
              end

%             !***
%             !*** check for 0,2pi crossings
%             !***

              if (vec1_lon >  pi)
                vec1_lon = vec1_lon - pi2;
              elseif (vec1_lon < -pi)
                vec1_lon = vec1_lon + pi2;
              end
              if (vec2_lon >  pi)
                vec2_lon = vec2_lon - pi2;
              elseif (vec2_lon < -pi)
                vec2_lon = vec2_lon + pi2;
              end

              cross_product = vec1_lon*vec2_lat - vec2_lon*vec1_lat;

%             !***
%             !*** if the cross product for a side is zero, the point 
%             !***   lies exactly on the side or the side is degenerate
%             !***   (zero length).  if degenerate, set the cross 
%             !***   product to a positive number.  otherwise perform 
%             !***   another cross product between the side and the 
%             !***   segment itself. 
%             !*** if this cross product is also zero, the line is 
%             !***   coincident with the cell boundary - perform the 
%             !***   dot product and only choose the cell if the dot 
%             !***   product is positive (parallel vs anti-parallel).
%             !***

              if (cross_product == 0)
                if ((vec1_lat ~= 0) | (vec1_lon ~= 0))
                  vec2_lat = lat2 - lat1;
                  vec2_lon = lon2 - lon1;

                  if (vec2_lon >  pi)
                    vec2_lon = vec2_lon - pi2;
                  elseif (vec2_lon < -pi)
                    vec2_lon = vec2_lon + pi2;
                  end

                  cross_product = vec1_lon*vec2_lat - vec2_lon*vec1_lat;
                else
                  cross_product = 1;
                end

                if (cross_product == 0)
                  lcoinc = true;
                  cross_product = vec1_lon*vec2_lon + vec1_lat*vec2_lat;
                  if (lrevers); cross_product = -cross_product; end
                end
              end

%             !***
%             !*** if cross product is less than zero, this cell
%             !*** doesn't work
%             !***

              if (cross_product < 0); src_corners=0; break; end %exit corner_loop
              jcw_count=jcw_count+1;
            end %do corner_loop
            src_corners=0;
          end

%         !***
%         !*** if cross products all positive, we found the location
%         !***

%          if (n > srch_corners)
           if (jcw_count+1 > srch_corners)
            location = srch_add(cell);

%             !***
%             !*** if the beginning of this segment was outside the
%             !*** grid, invert the segment so the intersection found
%             !*** will be the first intersection with the grid
%             !***

            if (loutside)
              lat2 = beglat;
              lon2 = beglon;
              location = 0;
              eps  = -tiny;
            else
              eps  = tiny;
            end

            srch_loop=0;
            break
%           break %exit srch_loop
          end
          if (srch_loop==0); break; end

%           !***
%           !*** otherwise move on to next cell
%           !***

        end %do cell_loop

        if (srch_loop==0); break; end
%         !***
%         !*** if still no cell found, the point lies outside the grid.
%         !***   take some baby steps along the segment to see if any
%         !***   part of the segment lies inside the grid.  
%         !***

        loutside = true;
        s1 = s1 + 0.001;
        lat1 = beglat + s1*(endlat - beglat);
        lon1 = beglon + s1*(lon2   - beglon);

%         !***
%         !*** reached the end of the segment and still outside the grid
%         !*** return no intersection
%         !***

        if (s1 >= 1); return; end

%        srch_loop=0;
      end %do srch_loop
%     call timer_stop(12)

% !-----------------------------------------------------------------------
% !
% !     now that a cell is found, search for the next intersection.
% !     loop over sides of the cell to find intersection with side
% !     must check all sides for coincidences or intersections
% !
% !-----------------------------------------------------------------------

%     call timer_start(13)
      intrsct_loop=1;
      while(intrsct_loop)
        for n=1:srch_corners
          next_n = mod(n,srch_corners) + 1;

          grdlon1 = srch_corner_lon(n     ,cell);
          grdlon2 = srch_corner_lon(next_n,cell);
          grdlat1 = srch_corner_lat(n     ,cell);
          grdlat2 = srch_corner_lat(next_n,cell);

%         !***
%         !*** set up linear system to solve for intersection
%         !***

          mat1 = lat2 - lat1;
          mat2 = grdlat1 - grdlat2;
          mat3 = lon2 - lon1;
          mat4 = grdlon1 - grdlon2;
          rhs1 = grdlat1 - lat1;
          rhs2 = grdlon1 - lon1;

          if (mat3 >  pi)
            mat3 = mat3 - pi2;
          elseif (mat3 < -pi)
            mat3 = mat3 + pi2;
          end
          if (mat4 >  pi)
            mat4 = mat4 - pi2;
          elseif (mat4 < -pi)
            mat4 = mat4 + pi2;
          end
          if (rhs2 >  pi)
            rhs2 = rhs2 - pi2;
          elseif (rhs2 < -pi)
            rhs2 = rhs2 + pi2;
          end

          determ = mat1*mat4 - mat2*mat3;

%         !***
%         !*** if the determinant is zero, the segments are either 
%         !***   parallel or coincident.  coincidences were detected 
%         !***   above so do nothing.
%         !*** if the determinant is non-zero, solve for the linear 
%         !***   parameters s for the intersection point on each line 
%         !***   segment.
%         !*** if 0<s1,s2<1 then the segment intersects with this side.
%         !***   return the point of intersection (adding a small
%         !***   number so the intersection is off the grid line).
%         !***

          if (abs(determ) > 1.e-30)

            s1 = (rhs1*mat4 - mat2*rhs2)/determ;
            s2 = (mat1*rhs2 - rhs1*mat3)/determ;

            if (s2 >= 0 & s2 <= 1 & ...
                s1 >  0 & s1 <= 1)

%             !***
%             !*** recompute intersection based on full segment
%             !*** so intersections are consistent for both sweeps
%             !***

              if (~loutside)
                mat1 = lat2 - begseg(1);
                mat3 = lon2 - begseg(2);
                rhs1 = grdlat1 - begseg(1);
                rhs2 = grdlon1 - begseg(2);
              else
                mat1 = begseg(1) - endlat;
                mat3 = begseg(2) - endlon;
                rhs1 = grdlat1 - endlat;
                rhs2 = grdlon1 - endlon;
              end

              if (mat3 >  pi)
                mat3 = mat3 - pi2;
              elseif (mat3 < -pi)
                mat3 = mat3 + pi2;
              end
              if (rhs2 >  pi)
                rhs2 = rhs2 - pi2;
              elseif (rhs2 < -pi)
                rhs2 = rhs2 + pi2;
              end

              determ = mat1*mat4 - mat2*mat3;

%             !***
%             !*** sometimes due to roundoff, the previous 
%             !*** determinant is non-zero, but the lines
%             !*** are actually coincident.  if this is the
%             !*** case, skip the rest.
%             !***

              if (determ ~= 0)
                s1 = (rhs1*mat4 - mat2*rhs2)/determ;
                s2 = (mat1*rhs2 - rhs1*mat3)/determ;

                offset = s1 + eps/determ;
                if (offset > 1); offset = 1; end

                if (~loutside)
                  intrsct_lat = begseg(1) + mat1*s1;
                  intrsct_lon = begseg(2) + mat3*s1;
                  intrsct_lat_off = begseg(1) + mat1*offset;
                  intrsct_lon_off = begseg(2) + mat3*offset;
                else
                  intrsct_lat = endlat + mat1*s1;
                  intrsct_lon = endlon + mat3*s1;
                  intrsct_lat_off = endlat + mat1*offset;
                  intrsct_lon_off = endlon + mat3*offset;
                end
                intrsct_loop=0;
                break %exit intrsct_loop
              end
              if (intrsct_loop==0); break; end
            end
            if (intrsct_loop==0); break; end
          end

%         !***
%         !*** no intersection this side, move on to next side
%         !***
        end %do intrsct_loop
        intrsct_loop=0;
      end %while intrsct_loop
%     call timer_stop(13)

% !-----------------------------------------------------------------------
% !
% !     if the segment crosses a pole threshold, reset the intersection
% !     to be the threshold latitude.  only check if this was not a
% !     threshold segment since sometimes coordinate transform can end
% !     up on other side of threshold again.
% !
% !-----------------------------------------------------------------------

      if (lthresh)
        if (intrsct_lat < north_thresh || intrsct_lat > south_thresh); lthresh = false; end
      elseif (lat1 > 0 & intrsct_lat > north_thresh)
          intrsct_lat = north_thresh + tiny;
          intrsct_lat_off = north_thresh + eps*mat1;
          s1 = (intrsct_lat - begseg(1))/mat1;
          intrsct_lon     = begseg(2) + s1*mat3;
          intrsct_lon_off = begseg(2) + (s1+eps)*mat3;
          last_loc = location;
          lthresh = true;
      elseif (lat1 < 0 & intrsct_lat < south_thresh)
          intrsct_lat = south_thresh - tiny;
          intrsct_lat_off = south_thresh + eps*mat1;
          s1 = (intrsct_lat - begseg(1))/mat1;
          intrsct_lon     = begseg(2) + s1*mat3;
          intrsct_lon_off = begseg(2) + (s1+eps)*mat3;
          last_loc = location;
          lthresh = true;
      end

%!-----------------------------------------------------------------------

      end %subroutine intersection

%!***********************************************************************

      function pole_intersection(location, ...
                   intrsct_lat,intrsct_lon,lcoinc,lthresh, ...
                   beglat, beglon, endlat, endlon, begseg, lrevers)

% !-----------------------------------------------------------------------
% !
% !     this routine is identical to the intersection routine except
% !     that a coordinate transformation (using a Lambert azimuthal
% !     equivalent projection) is performed to treat polar cells more
% !     accurately.
% !
% !-----------------------------------------------------------------------
% 
% !-----------------------------------------------------------------------
% !
% !     intent(in): 
% !
% !-----------------------------------------------------------------------

%       real (kind=dbl_kind), intent(in) :: 
%      &     beglat, beglon,  ! beginning lat/lon endpoints for segment
%      &     endlat, endlon   ! ending    lat/lon endpoints for segment
% 
%       real (kind=dbl_kind), dimension(2), intent(inout) :: 
%      &     begseg ! begin lat/lon of full segment
% 
%       logical (kind=log_kind), intent(in) ::
%      &        lrevers   ! flag true if segment integrated in reverse
% 
% !-----------------------------------------------------------------------
% !
% !     intent(out): 
% !
% !-----------------------------------------------------------------------
% 
%       integer (kind=int_kind), intent(inout) ::
%      &        location  ! address in destination array containing this
%                         ! segment -- also may contain last location on
%                         ! entry
% 
%       logical (kind=log_kind), intent(out) ::
%      &        lcoinc    ! flag segment coincident with grid line
% 
%       logical (kind=log_kind), intent(inout) ::
%      &        lthresh   ! flag segment crossing threshold boundary
% 
%       real (kind=dbl_kind), intent(out) ::
%      &     intrsct_lat, intrsct_lon ! lat/lon coords of next intersect.
% 
% !-----------------------------------------------------------------------
% !
% !     local variables
% !
% !-----------------------------------------------------------------------

%       integer (kind=int_kind) :: n, next_n, cell, srch_corners, pole_loc
% 
%       logical (kind=log_kind) :: loutside ! flags points outside grid
% 
%       real (kind=dbl_kind) :: pi4, rns, ! north/south conversion
%      &     x1, x2,       ! local x variables for segment
%      &     y1, y2,       ! local y variables for segment
%      &     begx, begy,   ! beginning x,y variables for segment
%      &     endx, endy,   ! beginning x,y variables for segment
%      &     begsegx, begsegy,   ! beginning x,y variables for segment
%      &     grdx1, grdx2, ! local x variables for grid cell
%      &     grdy1, grdy2, ! local y variables for grid cell
%      &     vec1_y, vec1_x, ! vectors and cross products used
%      &     vec2_y, vec2_x, ! during grid search
%      &     cross_product, eps, ! eps=small offset away from intersect
%      &     s1, s2, determ,     ! variables used for linear solve to
%      &     mat1, mat2, mat3, mat4, rhs1, rhs2  ! find intersection
% 
%       real (kind=dbl_kind), dimension(:,:), allocatable ::
%      &     srch_corner_x,  ! x of each corner of srch cells
%      &     srch_corner_y   ! y of each corner of srch cells
% 
%       !***
%       !*** save last intersection to avoid roundoff during coord
%       !*** transformation
%       !***

%use grids                      ! module with grid information
       global grid1_size grid2_size grid1_rank grid2_rank
       global grid1_corners grid2_corners grid1_dims grid2_dims
       global grid1_name grid2_name grid1_units grid2_units grid1_mask grid2_mask
       global grid1_center_lat grid1_center_lon grid2_center_lat grid2_center_lon
       global grid1_area grid2_area grid1_area_in grid2_area_in grid1_frac grid2_frac
       global grid1_corner_lat grid1_corner_lon grid2_corner_lat grid2_corner_lon
       global luse_grid_centers luse_grid1_area luse_grid2_area grid1_bound_box grid2_bound_box
       global restrict_type num_srch_bins bin_addr1 bin_addr2 bin_lats bin_lons

%use remap_vars                 ! common remapping variables
       global norm_opt_none norm_opt_dstarea norm_opt_frcarea
       global map_type_conserv map_type_bilinear map_type_bicubic map_type_distwgt
       global max_links_map1 num_links_map1 max_links_map2 num_links_map2
       global num_maps num_wts map_type norm_opt  resize_increment
       global grid1_add_map1 grid2_add_map1 grid1_add_map2 grid2_add_map2
       global wts_map1 wts_map2

%use remap_conservative         ! routines for conservative remap
      global num_srch_cells srch_add north_thresh south_thresh srch_corner_lat srch_corner_lon

%intersection
      global last_loc lthresh intrsct_lat_off intrsct_lon_off

%pole intersection
      global luse_last intrsct_x intrsct_y avoid_pole_count avoid_pole_offset
      global intrsct_lat intrsct_lon location

%store_link_cnsrv
      global first_call
      
%      logical (kind=log_kind), save :: luse_last = .false.

%      real (kind=dbl_kind), save :: 
%     &     intrsct_x, intrsct_y  ! x,y for intersection

%      !***
%      !*** variables necessary if segment manages to hit pole
%      !***

%      integer (kind=int_kind), save :: 
%     &     avoid_pole_count = 0  ! count attempts to avoid pole

%      real (kind=dbl_kind), save :: 
%     &     avoid_pole_offset = tiny  ! endpoint offset to avoid pole

% !-----------------------------------------------------------------------
% !
% !     initialize defaults, flags, etc.
% !
% !-----------------------------------------------------------------------
% 
      if (~lthresh); location = 0; end
      lcoinc = false;
      intrsct_lat = endlat;
      intrsct_lon = endlon;

      loutside = false;
      s1 = 0;

% !-----------------------------------------------------------------------
% !
% !     convert coordinates
% !
% !-----------------------------------------------------------------------

%      allocate(srch_corner_x(size(srch_corner_lat,DIM=1),
%     &                       size(srch_corner_lat,DIM=2)),
%     &         srch_corner_y(size(srch_corner_lat,DIM=1),
%     &                       size(srch_corner_lat,DIM=2)))
      srch_corner_x=zeros(size(srch_corner_lat,1),size(srch_corner_lat,2));
      srch_corner_y=zeros(size(srch_corner_lat,1),size(srch_corner_lat,2));

      if (beglat > zero)
        pi4 = quart*pi;
        rns = 1;
      else
        pi4 = -quart*pi;
        rns = -1;
      end

      if (luse_last)
        x1 = intrsct_x;
        y1 = intrsct_y;
      else
        x1 = rns*two*sin(pi4 - half*beglat)*cos(beglon);
        y1 =     two*sin(pi4 - half*beglat)*sin(beglon);
        luse_last = true;
      end
      x2 = rns*two*sin(pi4 - half*endlat)*cos(endlon);
      y2 =     two*sin(pi4 - half*endlat)*sin(endlon);
      srch_corner_x = rns*two*sin(pi4 - half*srch_corner_lat)* ...
                              cos(srch_corner_lon);
      srch_corner_y =     two*sin(pi4 - half*srch_corner_lat)* ...
                              sin(srch_corner_lon);

      begx = x1;
      begy = y1;
      endx = x2;
      endy = y2;
      begsegx = rns*two*sin(pi4 - half*begseg(1))*cos(begseg(2));
      begsegy =     two*sin(pi4 - half*begseg(1))*sin(begseg(2));
      intrsct_x = endx;
      intrsct_y = endy;

% !-----------------------------------------------------------------------
% !
% !     search for location of this segment in ocean grid using cross
% !     product method to determine whether a point is enclosed by a cell
% !
% !-----------------------------------------------------------------------

%     call timer_start(12)
      srch_corners = size(srch_corner_lat,1);
      srch_loop=1;
      while(srch_loop)
%     srch_loop: do

%        !***
%        !*** if last segment crossed threshold, use that location
%        !***

        if (lthresh)
          for cell=1:num_srch_cells
            if (srch_add(cell) == location)
              eps = tiny;
              srch_loop=0;
              break %srch_loop
            end
          end
        end

%         !***
%         !*** otherwise normal search algorithm
%         !***

        for cell=1:num_srch_cells
          srch_corner=1;
          while(srch_corner)
            for n=1:srch_corners
              next_n = MOD(n,srch_corners) + 1;

%             !***
%             !*** here we take the cross product of the vector making 
%             !*** up each cell side with the vector formed by the vertex
%             !*** and search point.  if all the cross products are 
%             !*** positive, the point is contained in the cell.
%             !***

              vec1_x = srch_corner_x(next_n,cell) - ...
                       srch_corner_x(n     ,cell);
              vec1_y = srch_corner_y(next_n,cell) - ...
                       srch_corner_y(n     ,cell);
              vec2_x = x1 - srch_corner_x(n,cell);
              vec2_y = y1 - srch_corner_y(n,cell);

%             !***
%             !*** if endpoint coincident with vertex, offset
%             !*** the endpoint
%             !***

              if (vec2_x == 0 & vec2_y == 0)
                x1 = x1 + 1.d-10*(x2-x1);
                y1 = y1 + 1.d-10*(y2-y1);
                vec2_x = x1 - srch_corner_x(n,cell);
                vec2_y = y1 - srch_corner_y(n,cell);
              end

              cross_product = vec1_x*vec2_y - vec2_x*vec1_y;

%             !***
%             !*** if the cross product for a side is zero, the point 
%             !***   lies exactly on the side or the length of a side
%             !***   is zero.  if the length is zero set det > 0.
%             !***   otherwise, perform another cross 
%             !***   product between the side and the segment itself. 
%             !*** if this cross product is also zero, the line is 
%             !***   coincident with the cell boundary - perform the 
%             !***   dot product and only choose the cell if the dot 
%             !***   product is positive (parallel vs anti-parallel).
%             !***

              if (cross_product == 0)
                if (vec1_x ~= zero | vec1_y ~= 0)
                  vec2_x = x2 - x1;
                  vec2_y = y2 - y1;
                  cross_product = vec1_x*vec2_y - vec2_x*vec1_y;
                else
                  cross_product = 1;
                end

                if (cross_product == 0)
                  lcoinc = true;
                  cross_product = vec1_x*vec2_x + vec1_y*vec2_y;
                  if (lrevers); cross_product = -cross_product; end
                end
              end

%             !***
%             !*** if cross product is less than zero, this cell
%             !*** doesn't work
%             !***

              if (cross_product < zero); srch_corner=0; end%exit corner_loop

            end %do corner_loop
            srch_corner=0;
          end %while

%           !***
%           !*** if cross products all positive, we found the location
%           !***

          if (n > srch_corners)
            location = srch_add(cell);

%             !***
%             !*** if the beginning of this segment was outside the
%             !*** grid, invert the segment so the intersection found
%             !*** will be the first intersection with the grid
%             !***

            if (loutside)
              x2 = begx;
              y2 = begy;
              location = 0;
              eps  = -tiny;
            else
              eps  = tiny;
            end

            srch_loop=0;
            break %exit srch_loop
          end

%           !***
%           !*** otherwise move on to next cell
%           !***

        end %do cell_loop

%         !***
%         !*** if no cell found, the point lies outside the grid.
%         !***   take some baby steps along the segment to see if any
%         !***   part of the segment lies inside the grid.  
%         !***

        loutside = true;
        s1 = s1 + 0.001;
        x1 = begx + s1*(x2 - begx);
        y1 = begy + s1*(y2 - begy);

%         !***
%         !*** reached the end of the segment and still outside the grid
%         !*** return no intersection
%         !***

        if (s1 >= 1)
          clear srch_corner_x srch_corner_y
          luse_last = false;
          return
        end

        srch_loop=0;
      end %do srch_loop
%     call timer_stop(12)

% !-----------------------------------------------------------------------
% !
% !     now that a cell is found, search for the next intersection.
% !     loop over sides of the cell to find intersection with side
% !     must check all sides for coincidences or intersections
% !
% !-----------------------------------------------------------------------

%     call timer_start(13)
      intrsct_loop=1;
      while(intrsct_loop)
      for n=1:srch_corners
        next_n = mod(n,srch_corners) + 1;

        grdy1 = srch_corner_y(n     ,cell);
        grdy2 = srch_corner_y(next_n,cell);
        grdx1 = srch_corner_x(n     ,cell);
        grdx2 = srch_corner_x(next_n,cell);

%         !***
%         !*** set up linear system to solve for intersection
%         !***

        mat1 = x2 - x1;
        mat2 = grdx1 - grdx2;
        mat3 = y2 - y1;
        mat4 = grdy1 - grdy2;
        rhs1 = grdx1 - x1;
        rhs2 = grdy1 - y1;

        determ = mat1*mat4 - mat2*mat3;

%         !***
%         !*** if the determinant is zero, the segments are either 
%         !***   parallel or coincident or one segment has zero length.  
%         !***   coincidences were detected above so do nothing.
%         !*** if the determinant is non-zero, solve for the linear 
%         !***   parameters s for the intersection point on each line 
%         !***   segment.
%         !*** if 0<s1,s2<1 then the segment intersects with this side.
%         !***   return the point of intersection (adding a small
%         !***   number so the intersection is off the grid line).
%         !***

        if (abs(determ) > 1.e-30)

          s1 = (rhs1*mat4 - mat2*rhs2)/determ;
          s2 = (mat1*rhs2 - rhs1*mat3)/determ;

          if (s2 >= 0 & s2 <= 1 & ...
              s1 >  0 & s1 <= 1)

%             !***
%             !*** recompute intersection using entire segment
%             !*** for consistency between sweeps
%             !***

            if (~loutside)
              mat1 = x2 - begsegx;
              mat3 = y2 - begsegy;
              rhs1 = grdx1 - begsegx;
              rhs2 = grdy1 - begsegy;
            else
              mat1 = x2 - endx;
              mat3 = y2 - endy;
              rhs1 = grdx1 - endx;
              rhs2 = grdy1 - endy;
            end

            determ = mat1*mat4 - mat2*mat3;

%             !***
%             !*** sometimes due to roundoff, the previous 
%             !*** determinant is non-zero, but the lines
%             !*** are actually coincident.  if this is the
%             !*** case, skip the rest.
%             !***

            if (determ ~= 0)
              s1 = (rhs1*mat4 - mat2*rhs2)/determ;
              s2 = (mat1*rhs2 - rhs1*mat3)/determ;

              if (~loutside)
                intrsct_x = begsegx + s1*mat1;
                intrsct_y = begsegy + s1*mat3;
              else
                intrsct_x = endx + s1*mat1;
                intrsct_y = endy + s1*mat3;
              end

%               !***
%               !*** convert back to lat/lon coordinates
%               !***

              intrsct_lon = rns*atan2(intrsct_y,intrsct_x);
              if (intrsct_lon < 0); intrsct_lon = intrsct_lon + pi2; end
              if (abs(intrsct_x) > 1.d-10)
                intrsct_lat = (pi4 - ...
                  asin(rns*half*intrsct_x/cos(intrsct_lon)))*2 ;
              elseif (abs(intrsct_y) > 1.d-10)
                 intrsct_lat = (pi4 - ...
                 asin(half*intrsct_y/sin(intrsct_lon)))*2
              else
                intrsct_lat = two*pi4;
              end

%               !***
%               !*** add offset in transformed space for next pass.
%               !***

              if (s1 - eps/determ < 1)
                intrsct_x = intrsct_x - mat1*(eps/determ);
                intrsct_y = intrsct_y - mat3*(eps/determ);
              else
                if (~loutside)
                  intrsct_x = endx;
                  intrsct_y = endy;
                  intrsct_lat = endlat;
                  intrsct_lon = endlon;
                else
                  intrsct_x = begsegx;
                  intrsct_y = begsegy;
                  intrsct_lat = begseg(1);
                  intrsct_lon = begseg(2);
                end
              end

              intrsct_loop=0; %exit intrsct_loop
            end
          end
        end

%         !***
%         !*** no intersection this side, move on to next side
%         !***
 
      end %do intrsct_loop
      intrsct_loop=0;
      end
%     call timer_stop(13)

      clear srch_corner_x srch_corner_y

% !-----------------------------------------------------------------------
% !
% !     if segment manages to cross over pole, shift the beginning 
% !     endpoint in order to avoid hitting pole directly
% !     (it is ok for endpoint to be pole point)
% !
% !-----------------------------------------------------------------------

      if (abs(intrsct_x) < 1.e-10 & abs(intrsct_y) < 1.e-10 & ...
          (endx ~= 0 & endy ~=0))
        if (avoid_pole_count > 2)
           avoid_pole_count = 0;
           avoid_pole_offset = 10.*avoid_pole_offset;
        end

        cross_product = begsegx*(endy-begsegy) - begsegy*(endx-begsegx);
        intrsct_lat = begseg(1);
        if (cross_product*intrsct_lat > 0)
          intrsct_lon = beglon    + avoid_pole_offset;
          begseg(2)   = begseg(2) + avoid_pole_offset;
        else
          intrsct_lon = beglon    - avoid_pole_offset;
          begseg(2)   = begseg(2) - avoid_pole_offset;
        end

        avoid_pole_count = avoid_pole_count + 1;
        luse_last = false;
      else
        avoid_pole_count = 0;
        avoid_pole_offset = tiny;
      end

% !-----------------------------------------------------------------------
% !
% !     if the segment crosses a pole threshold, reset the intersection
% !     to be the threshold latitude and do not reuse x,y intersect
% !     on next entry.  only check if did not cross threshold last
% !     time - sometimes the coordinate transformation can place a
% !     segment on the other side of the threshold again
% !
% !-----------------------------------------------------------------------

      if (lthresh)
        if (intrsct_lat > north_thresh | intrsct_lat < south_thresh);lthresh = false; end
      elseif (beglat > 0 & intrsct_lat < north_thresh)
        mat4 = endlat - begseg(1);
        mat3 = endlon - begseg(2);
        if (mat3 >  pi);  mat3 = mat3 - pi2; end
        if (mat3 < -pi);  mat3 = mat3 + pi2; end
        intrsct_lat = north_thresh - tiny;
        s1 = (north_thresh - begseg(1))/mat4;
        intrsct_lon = begseg(2) + s1*mat3;
        luse_last = false;
        lthresh = true;
      elseif (beglat < 0 & intrsct_lat > south_thresh)
        mat4 = endlat - begseg(1);
        mat3 = endlon - begseg(2);
        if (mat3 >  pi); mat3 = mat3 - pi2; end
        if (mat3 < -pi); mat3 = mat3 + pi2; end
        intrsct_lat = south_thresh + tiny;
        s1 = (south_thresh - begseg(1))/mat4;
        intrsct_lon = begseg(2) + s1*mat3;
        luse_last = false;
        lthresh = true;
      end

%       !***
%       !*** if reached end of segment, do not use x,y intersect 
%       !*** on next entry
%       !***

      if (intrsct_lat == endlat & intrsct_lon == endlon)
        luse_last = false;
      end

%!-----------------------------------------------------------------------

      end %subroutine pole_intersection

%!***********************************************************************

      function line_integral(weights, num_wts, ...
                            in_phi1, in_phi2, theta1, theta2, ...
                            grid1_lat, grid1_lon, grid2_lat, grid2_lon)

%use grids                      ! module with grid information
       global grid1_size grid2_size grid1_rank grid2_rank
       global grid1_corners grid2_corners grid1_dims grid2_dims
       global grid1_name grid2_name grid1_units grid2_units grid1_mask grid2_mask
       global grid1_center_lat grid1_center_lon grid2_center_lat grid2_center_lon
       global grid1_area grid2_area grid1_area_in grid2_area_in grid1_frac grid2_frac
       global grid1_corner_lat grid1_corner_lon grid2_corner_lat grid2_corner_lon
       global luse_grid_centers luse_grid1_area luse_grid2_area grid1_bound_box grid2_bound_box
       global restrict_type num_srch_bins bin_addr1 bin_addr2 bin_lats bin_lons

%use remap_vars                 ! common remapping variables
       global norm_opt_none norm_opt_dstarea norm_opt_frcarea
       global map_type_conserv map_type_bilinear map_type_bicubic map_type_distwgt
       global max_links_map1 num_links_map1 max_links_map2 num_links_map2
       global num_maps num_wts map_type norm_opt  resize_increment
       global grid1_add_map1 grid2_add_map1 grid1_add_map2 grid2_add_map2
       global wts_map1 wts_map2

%use remap_conservative         ! routines for conservative remap
      global num_srch_cells srch_add north_thresh south_thresh srch_corner_lat srch_corner_lon

%intersection
      global last_loc lthresh intrsct_lat_off intrsct_lon_off
      global intrsct_lat intrsct_lon location

%pole intersection
      global luse_last intrsct_x intrsct_y avoid_pole_count avoid_pole_offset

%store_link_cnsrv
      global first_call

% line integral
     global weights

% !-----------------------------------------------------------------------
% !
% !     this routine computes the line integral of the flux function 
% !     that results in the interpolation weights.  the line is defined
% !     by the input lat/lon of the endpoints.
% !
% !-----------------------------------------------------------------------
% 
% !-----------------------------------------------------------------------
% !
% !     intent(in):
% !
% !-----------------------------------------------------------------------

%      integer (kind=int_kind), intent(in) ::
%     &        num_wts  ! number of weights to compute

%      real (kind=dbl_kind), intent(in) :: 
%     &     in_phi1, in_phi2,     ! longitude endpoints for the segment
%     &     theta1, theta2,       ! latitude  endpoints for the segment
%     &     grid1_lat, grid1_lon, ! reference coordinates for each
%     &     grid2_lat, grid2_lon  ! grid (to ensure correct 0,2pi interv.

%!-----------------------------------------------------------------------
%!
%!     intent(out):
%!
%!-----------------------------------------------------------------------

%      real (kind=dbl_kind), dimension(2*num_wts), intent(out) ::
%     &     weights   ! line integral contribution to weights

%!-----------------------------------------------------------------------
%!
%!     local variables
%!
%!-----------------------------------------------------------------------

%      real (kind=dbl_kind) :: dphi, sinth1, sinth2, costh1, costh2, fac,
%     &                        phi1, phi2, phidiff1, phidiff2, sinint
%      real (kind=dbl_kind) :: f1, f2, fint

%!-----------------------------------------------------------------------
%!
%!     weights for the general case based on a trapezoidal approx to
%!     the integrals.
%!
%!-----------------------------------------------------------------------

      sinth1 = sin(theta1);
      sinth2 = sin(theta2);
      costh1 = cos(theta1);
      costh2 = cos(theta2);

      dphi = in_phi1 - in_phi2;
      if (dphi >  pi)
        dphi = dphi - 2*pi;
      elseif (dphi < -pi)
        dphi = dphi + 2*pi;
      end
      dphi = 0.5*dphi;

%!-----------------------------------------------------------------------
%!
%!     the first weight is the area overlap integral. the second and
%!     fourth are second-order latitude gradient weights.
%!
%!-----------------------------------------------------------------------

      weights(        1) = dphi*(sinth1 + sinth2);
      weights(num_wts+1) = dphi*(sinth1 + sinth2);
      weights(        2) = dphi*(costh1 + costh2 + (theta1*sinth1 + ...
                                                    theta2*sinth2));
      weights(num_wts+2) = dphi*(costh1 + costh2 + (theta1*sinth1 + ...
                                                    theta2*sinth2));

% !-----------------------------------------------------------------------
% !
% !     the third and fifth weights are for the second-order phi gradient
% !     component.  must be careful of longitude range.
% !
% !-----------------------------------------------------------------------

      f1 = 0.5*(costh1*sinth1 + theta1);
      f2 = 0.5*(costh2*sinth2 + theta2);

      phi1 = in_phi1 - grid1_lon;
      if (phi1 >  pi)
        phi1 = phi1 - 2*pi;
      elseif (phi1 < -pi)
        phi1 = phi1 + 2*pi;
      end

      phi2 = in_phi2 - grid1_lon;
      if (phi2 >  pi)
        phi2 = phi2 - 2*pi;
      elseif (phi2 < -pi)
        phi2 = phi2 + 2*pi;
      end

      if ((phi2-phi1) <  pi & (phi2-phi1) > -pi)
        weights(3) = dphi*(phi1*f1 + phi2*f2);
      else
        if (phi1 > zero)
          fac = pi;
        else
          fac = -pi;
        end
        fint = f1 + (f2-f1)*(fac-phi1)/abs(dphi);
        weights(3) = 0.5*phi1*(phi1-fac)*f1 - ...
                    0.5*phi2*(phi2+fac)*f2 + ...
                    0.5*fac*(phi1+phi2)*fint;
      end

      phi1 = in_phi1 - grid2_lon;
      if (phi1 >  pi)
        phi1 = phi1 - 2*pi;
      elseif (phi1 < -pi)
        phi1 = phi1 + 2*pi;
      end

      phi2 = in_phi2 - grid2_lon;
      if (phi2 >  pi)
        phi2 = phi2 - 2*pi;
      elseif (phi2 < -pi)
        phi2 = phi2 + 2*pi;
      end

      if ((phi2-phi1) <  pi && (phi2-phi1) > -pi)
        weights(num_wts+3) = dphi*(phi1*f1 + phi2*f2);
      else
        if (phi1 > zero)
          fac = pi;
        else
          fac = -pi;
        end
        fint = f1 + (f2-f1)*(fac-phi1)/abs(dphi);
        weights(num_wts+3) = half*phi1*(phi1-fac)*f1 - ...
                             half*phi2*(phi2+fac)*f2 + ...
                             half*fac*(phi1+phi2)*fint;
      end

%!-----------------------------------------------------------------------

      end %subroutine line_integral

%!***********************************************************************

      function store_link_cnsrv(add1, add2, weights)

% !-----------------------------------------------------------------------
% !
% !     this routine stores the address and weight for this link in
% !     the appropriate address and weight arrays and resizes those
% !     arrays if necessary.
% !
% !-----------------------------------------------------------------------
% 
% !-----------------------------------------------------------------------
% !
% !     input variables
% !
% !-----------------------------------------------------------------------
% 
%       integer (kind=int_kind), intent(in) ::
%      &        add1,  ! address on grid1
%      &        add2   ! address on grid2
% 
%       real (kind=dbl_kind), dimension(:), intent(in) ::
%      &        weights ! array of remapping weights for this link
% 
% !-----------------------------------------------------------------------
% !
% !     local variables
% !
% !-----------------------------------------------------------------------
% 
%       integer (kind=int_kind) :: nlink, min_link, max_link ! link index
% 
%       integer (kind=int_kind), dimension(:,:), allocatable, save ::
%      &        link_add1,  ! min,max link add to restrict search
%      &        link_add2   ! min,max link add to restrict search

%     logical (kind=log_kind), save :: first_call = .true.

%use grids                      ! module with grid information
       global grid1_size grid2_size grid1_rank grid2_rank
       global grid1_corners grid2_corners grid1_dims grid2_dims
       global grid1_name grid2_name grid1_units grid2_units grid1_mask grid2_mask
       global grid1_center_lat grid1_center_lon grid2_center_lat grid2_center_lon
       global grid1_area grid2_area grid1_area_in grid2_area_in grid1_frac grid2_frac
       global grid1_corner_lat grid1_corner_lon grid2_corner_lat grid2_corner_lon
       global luse_grid_centers luse_grid1_area luse_grid2_area grid1_bound_box grid2_bound_box
       global restrict_type num_srch_bins bin_addr1 bin_addr2 bin_lats bin_lons

%use remap_vars                 ! common remapping variables
       global norm_opt_none norm_opt_dstarea norm_opt_frcarea
       global map_type_conserv map_type_bilinear map_type_bicubic map_type_distwgt
       global max_links_map1 num_links_map1 max_links_map2 num_links_map2
       global num_maps num_wts map_type norm_opt  resize_increment
       global grid1_add_map1 grid2_add_map1 grid1_add_map2 grid2_add_map2
       global wts_map1 wts_map2

%use remap_conservative         ! routines for conservative remap
      global num_srch_cells srch_add north_thresh south_thresh srch_corner_lat srch_corner_lon

%intersection
      global last_loc lthresh intrsct_lat_off intrsct_lon_off

%pole intersection
      global luse_last intrsct_x intrsct_y avoid_pole_count avoid_pole_offset
      global intrsct_lat intrsct_lon location

%store_link_cnsrv
      global first_call link_add1 link_add2

% line integral
      global weights

% !-----------------------------------------------------------------------
% !
% !     if all weights are zero, do not bother storing the link
% !
% !-----------------------------------------------------------------------

      if (all(weights == 0)); return; end

% !-----------------------------------------------------------------------
% !
% !     restrict the range of links to search for existing links
% !
% !-----------------------------------------------------------------------

      if (first_call)
        link_add1=zeros(2,grid1_size);
        link_add2=zeros(2,grid2_size);
        first_call = 0;
        min_link = 1;
        max_link = 0;
      else
        min_link = min(link_add1(1,add1),link_add2(1,add2));
        max_link = max(link_add1(2,add1),link_add2(2,add2));
        if (min_link == 0)
          min_link = 1;
          max_link = 0;
        end
      end

% !-----------------------------------------------------------------------
% !
% !     if the link already exists, add the weight to the current weight
% !     arrays
% !
% !-----------------------------------------------------------------------

      for nlink=min_link:max_link
        if (add1 == grid1_add_map1(nlink))
        if (add2 == grid2_add_map1(nlink))

          wts_map1(:,nlink) = wts_map1(:,nlink) + weights(1:num_wts)';
          if (num_maps == 2)
            wts_map2(:,nlink) = wts_map2(:,nlink) + ...
                                        weights(num_wts+1:2*num_wts);
          end
          return

        end
        end
      end

% !-----------------------------------------------------------------------
% !
% !     if the link does not yet exist, increment number of links and 
% !     check to see if remap arrays need to be increased to accomodate 
% !     the new link.  then store the link.
% !
% !-----------------------------------------------------------------------

      num_links_map1  = num_links_map1 + 1;
      if (num_links_map1 > max_links_map1) 
        resize_remap_vars(1,resize_increment)
      end

      grid1_add_map1(num_links_map1) = add1;
      grid2_add_map1(num_links_map1) = add2;
      wts_map1    (:,num_links_map1) = weights(1:num_wts);

      if (num_maps > 1)
        num_links_map2  = num_links_map2 + 1;
        if (num_links_map2 > max_links_map2)
          resize_remap_vars(2,resize_increment)
        end

        grid1_add_map2(num_links_map2) = add1;
        grid2_add_map2(num_links_map2) = add2;
        wts_map2    (:,num_links_map2) = weights(num_wts+1:2*num_wts);
      end

      if (link_add1(1,add1) == 0); link_add1(1,add1) = num_links_map1; end
      if (link_add2(1,add2) == 0); link_add2(1,add2) = num_links_map1; end
      link_add1(2,add1) = num_links_map1;
      link_add2(2,add2) = num_links_map1;

%!-----------------------------------------------------------------------

      end %subroutine store_link_cnsrv

%!***********************************************************************

%      end module remap_conservative


%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
