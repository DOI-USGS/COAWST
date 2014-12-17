function scrip
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%!
%!     This routine is the driver for computing the addresses and weights 
%!     for interpolating between two grids on a sphere.
%!
%!-----------------------------------------------------------------------
%!
%!     CVS:$Id: scrip.f,v 1.6 2001/08/21 21:06:44 pwjones Exp $
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
%      program scrip
%
%!-----------------------------------------------------------------------
%
%      use kinds_mod                  ! module defining data types
%      use constants                  ! module for common constants
%      use iounits                    ! I/O unit manager
%      use timers                     ! CPU timers
%      use grids                      ! module with grid information
%      use remap_vars                 ! common remapping variables
%      use remap_conservative         ! routines for conservative remap
%      use remap_distance_weight      ! routines for dist-weight remap
%      use remap_bilinear             ! routines for bilinear interp
%      use remap_bicubic              ! routines for bicubic  interp
%      use remap_write                ! routines for remap output
%

       warning('OFF', 'all')

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

%      use remap_conservative         ! routines for conservative remap
       global num_srch_cells srch_add north_thresh south_thresh srch_corner_lat srch_corner_lon

%      implicit none
%
%!-----------------------------------------------------------------------
%!
%!     input namelist variables
%!
%!-----------------------------------------------------------------------
%
%      character (char_len) :: 
%     &           grid1_file,   ! filename of grid file containing grid1
%     &           grid2_file,   ! filename of grid file containing grid2
%     &           interp_file1, ! filename for output remap data (map1)
%     &           interp_file2, ! filename for output remap data (map2)
%     &           map1_name,    ! name for mapping from grid1 to grid2
%     &           map2_name,    ! name for mapping from grid2 to grid1
%     &           map_method,   ! choice for mapping method
%     &           normalize_opt,! option for normalizing weights
%     &           output_opt    ! option for output conventions
%
%      integer (kind=int_kind) ::
%     &           nmap          ! number of mappings to compute (1 or 2)
%
%      namelist /remap_inputs/ grid1_file, grid2_file, 
%     &                        interp_file1, interp_file2,
%     &                        map1_name, map2_name, num_maps,
%     &                        luse_grid1_area, luse_grid2_area,
%     &                        map_method, normalize_opt, output_opt,
%     &                        restrict_type, num_srch_bins
%
%!-----------------------------------------------------------------------
%!
%!     local variables
%!
%!-----------------------------------------------------------------------
%
%      integer (kind=int_kind) :: n,     ! dummy counter
%     &                           iunit  ! unit number for namelist file
%
%!-----------------------------------------------------------------------
%!
%!     initialize timers
%!
%!-----------------------------------------------------------------------
%
%      call timers_init
%      do n=1,max_timers
%        call timer_clear(n)
%      end do
%
%!-----------------------------------------------------------------------
%!
%!     read input namelist
%!
%!-----------------------------------------------------------------------
%

      grid1_file    = 'unknown';
      grid2_file    = 'unknown';
      interp_file1  = 'unknown';
      interp_file2  = 'unknown';
      map1_name     = 'unknown';
      map2_name     = 'unknown';
      luse_grid1_area = false;
      luse_grid2_area = false;
      num_maps      = 2;
      map_type      = 1;
      normalize_opt = 'fracarea';
      output_opt    = 'scrip';
      restrict_type = 'latitude';
      num_srch_bins = 900;

%     call get_unit(iunit)
%     open(iunit, file='scrip_in', status='old', form='formatted')
%     read(iunit, nml=remap_inputs)
%     call release_unit(iunit)
      fid=fopen('scrip_in');
      for ii=1:500
        tline=fgetl(fid);
        if (tline==-1); break; end
        tline=[deblank(tline),'                '];
        zz=strfind(tline,'=')+1; if isempty(zz); zz=1;end
        if (strcmp(tline(1:8),'num_maps'));
          num_maps=str2num(tline(zz:end));
        elseif (strcmp(tline(1:10),'grid1_file')); 
          grid1_file=tline(zz:end);
          a=strfind(grid1_file,''''); grid1_file=grid1_file(a(1)+1:a(2)-1);
        elseif (strcmp(tline(1:10),'grid2_file')); 
          grid2_file=tline(zz:end);
          a=strfind(grid2_file,''''); grid2_file=grid2_file(a(1)+1:a(2)-1);
        elseif (strcmp(tline(1:12),'interp_file1')); 
          interp_file1=tline(zz:end);
          a=strfind(interp_file1,''''); interp_file1=interp_file1(a(1)+1:a(2)-1);
        elseif (strcmp(tline(1:12),'interp_file2')); 
          interp_file2=tline(zz:end);
          a=strfind(interp_file2,''''); interp_file2=interp_file2(a(1)+1:a(2)-1);
        elseif (strcmp(tline(1:9),'map1_name')); 
          map1_name=tline(zz:end);
          a=strfind(map1_name,''''); map1_name=map1_name(a(1)+1:a(2)-1);
        elseif (strcmp(tline(1:9),'map2_name')); 
          map2_name=tline(zz:end);
          a=strfind(map2_name,''''); map2_name=map2_name(a(1)+1:a(2)-1);
        elseif (strcmp(tline(1:10),'map_method')); 
          map_method=tline(zz:end);
          a=strfind(map_method,''''); map_method=map_method(a(1)+1:a(2)-1);
        elseif (strcmp(tline(1:13),'normalize_opt')); 
          normalize_opt=tline(zz:end);
          a=strfind(normalize_opt,''''); normalize_opt=normalize_opt(a(1)+1:a(2)-1);
        elseif (strcmp(tline(1:10),'output_opt')); 
          output_opt=tline(zz:end);
          a=strfind(output_opt,''''); output_opt=output_opt(a(1)+1:a(2)-1);
        elseif (strcmp(tline(1:13),'restrict_type')); 
          restrict_type=tline(zz:end);
          a=strfind(restrict_type,''''); restrict_type=restrict_type(a(1)+1:a(2)-1);
        elseif (strcmp(tline(1:13),'num_srch_bins')); 
          num_srch_bins=str2num(tline(zz:end));
        elseif (strcmp(tline(1:15),'luse_grid1_area')); 
          luse_grid1_area=deblank(tline(zz+1:end));
          luse_grid1_area=luse_grid1_area(2:end-1);
          if (strcmp('false',luse_grid1_area))
            luse_grid1_area=0;
          else
            luse_grid1_area=1;
          end
        elseif (strcmp(tline(1:15),'luse_grid2_area')); 
          luse_grid2_area=deblank(tline(zz+1:end));
          luse_grid2_area=luse_grid2_area(2:end-1);
          if (strcmp('false',luse_grid2_area))
            luse_grid2_area=0;
          else
            luse_grid2_area=1;
          end
        end
      end
      fclose(fid);

      switch (map_method)
      case ('conservative')
        map_type = 'map_type_conserv';
        luse_grid_centers = false;
      case ('bilinear')
        map_type = 'map_type_bilinear';
        luse_grid_centers = true;
      case ('bicubic')
        map_type = 'map_type_bicubic';
        luse_grid_centers = true;
      case ('distwgt')
        map_type = 'map_type_distwgt';
        luse_grid_centers = true;
      case default
        disp('unknown mapping method')
      end
 
      switch (normalize_opt(1:4))
      case ('none')
        norm_opt = 'norm_opt_none';
      case ('frac')
        norm_opt = 'norm_opt_frcarea';
      case ('dest')
        norm_opt = 'norm_opt_dstarea';
      case default
        disp('unknown normalization option');
      end

%!-----------------------------------------------------------------------
%!
%!     initialize grid information for both grids
%!
%!-----------------------------------------------------------------------

      scrip_grid_init(grid1_file, grid2_file)

      disp([' Computing remappings between: ',grid1_file])
      disp(['                          and  ',grid2_file])

%!-----------------------------------------------------------------------
%!
%!     initialize some remapping variables.
%!
%!-----------------------------------------------------------------------

      init_remap_vars(map_type)

%!-----------------------------------------------------------------------
%!
%!     call appropriate interpolation setup routine based on type of
%!     remapping requested.
%!
%!-----------------------------------------------------------------------

      switch (map_type)
      case('map_type_conserv')
        remap_conserv
      case('map_type_bilinear')
        remap_bilin
      case('map_type_distwgt')
        remap_distwgt
      case('map_type_bicubic')
        remap_bicub
      case default
        disp('Invalid Map Type')
      end
 
%!-----------------------------------------------------------------------
%!
%!     reduce size of remapping arrays and then write remapping info
%!     to a file.
%!
%!-----------------------------------------------------------------------

      if (num_links_map1 ~= max_links_map1)
        resize_remap_vars(1, num_links_map1-max_links_map1)
      end
      if ((num_maps > 1) & (num_links_map2 ~= max_links_map2))
        resize_remap_vars(2, num_links_map2-max_links_map2)
      end

     write_remap(map1_name, map2_name, ... 
                 interp_file1, interp_file2, output_opt)

%!-----------------------------------------------------------------------

      end %program scrip

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
