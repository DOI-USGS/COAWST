%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%!
%!     this module contains necessary variables for remapping between
%!     two grids.  also routines for resizing and initializing these
%!     variables.
%!
%!-----------------------------------------------------------------------
%!
%!     CVS:$Id: remap_vars.f,v 1.5 2000/04/19 21:56:26 pwjones Exp $
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
%
%      module remap_vars
%
%      use kinds_mod
%      use constants
%      use grids
%
%      implicit none
%
%!-----------------------------------------------------------------------
%!
%!     module variables
%!
%!-----------------------------------------------------------------------
%
%      integer (kind=int_kind), parameter ::
%       global norm_opt_none  norm_opt_dstarea norm_opt_frcarea
%       global map_type_conserv map_type_bilinear map_type_bicubic map_type_distwgt

%       norm_opt_none    = 1;
%       norm_opt_dstarea = 2;
%       norm_opt_frcarea = 3;

%      integer (kind=int_kind), parameter ::
%       map_type_conserv  = 1;
%       map_type_bilinear = 2;
%       map_type_bicubic  = 3;
%       map_type_distwgt  = 4;

%      integer (kind=int_kind), save :: 
%     &      max_links_map1  ! current size of link arrays
%     &,     num_links_map1  ! actual number of links for remapping
%     &,     max_links_map2  ! current size of link arrays
%     &,     num_links_map2  ! actual number of links for remapping
%     &,     num_maps        ! num of remappings for this grid pair
%     &,     num_wts         ! num of weights used in remapping
%     &,     map_type        ! identifier for remapping method
%     &,     norm_opt        ! option for normalization (conserv only)
%     &,     resize_increment ! default amount to increase array size

%      integer (kind=int_kind), dimension(:), allocatable, save ::
%     &      grid1_add_map1, ! grid1 address for each link in mapping 1
%     &      grid2_add_map1, ! grid2 address for each link in mapping 1
%     &      grid1_add_map2, ! grid1 address for each link in mapping 2
%     &      grid2_add_map2  ! grid2 address for each link in mapping 2

%      real (kind=dbl_kind), dimension(:,:), allocatable, save ::
%     &      wts_map1, ! map weights for each link (num_wts,max_links)
%     &      wts_map2  ! map weights for each link (num_wts,max_links)
%
%!***********************************************************************

%      contains

%!***********************************************************************

      function init_remap_vars(map_type)

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


       norm_opt_none    = 1;
       norm_opt_dstarea = 2;
       norm_opt_frcarea = 3;

       map_type_conserv  = 1;
       map_type_bilinear = 2;
       map_type_bicubic  = 3;
       map_type_distwgt  = 4;

%!-----------------------------------------------------------------------
%!
%!     this routine initializes some variables and provides an initial
%!     allocation of arrays (fairly large so frequent resizing 
%!     unnecessary).
%!
%!-----------------------------------------------------------------------

%!-----------------------------------------------------------------------
%!
%!     determine the number of weights
%!
%!-----------------------------------------------------------------------

      switch (map_type)
      case('map_type_conserv')
        num_wts = 3;
      case('map_type_bilinear')
        num_wts = 1;
      case('map_type_bicubic')
        num_wts = 4;
      case('map_type_distwgt')
        num_wts = 1;
      end

%!-----------------------------------------------------------------------
%!
%!     initialize num_links and set max_links to four times the largest 
%!     of the destination grid sizes initially (can be changed later).
%!     set a default resize increment to increase the size of link
%!     arrays if the number of links exceeds the initial size
%!   
%!-----------------------------------------------------------------------

      num_links_map1 = 0;
      max_links_map1 = 4*grid2_size;
      if (num_maps > 1)
        num_links_map2 = 0;
        max_links_map1 = max(4*grid1_size,4*grid2_size);
        max_links_map2 = max_links_map1;
      end

      resize_increment = 0.1*max(grid1_size,grid2_size);

%!-----------------------------------------------------------------------
%!
%!     allocate address and weight arrays for mapping 1
%!   
%!-----------------------------------------------------------------------

       grid1_add_map1=zeros(max_links_map1,1);
       grid2_add_map1=zeros(max_links_map1,1);
       wts_map1=zeros(num_wts, max_links_map1);

%!-----------------------------------------------------------------------
%!
%!     allocate address and weight arrays for mapping 2 if necessary 
%!   
%!-----------------------------------------------------------------------

      if (num_maps > 1)
       grid1_add_map2=zeros(max_links_map2,1);
       grid2_add_map2=zeros(max_links_map2,1);
       wts_map2=zeros(num_wts, max_links_map2);
      end

%!-----------------------------------------------------------------------

      end %subroutine init_remap_vars

%!***********************************************************************

%!***********************************************************************

%      end module remap_vars

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
