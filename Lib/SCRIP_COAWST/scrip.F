!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!     Note SCRIP_COAWST required the original file from SCRIP package 
!     to be converted into a module including a subroutine. 
!     -Arrays have been read in and deallocated in the end of file 
!        
!---- Written by John C. Warner-----------------------------------------
!-----         Tarandeep S. Kalra --------------------------------------
!--------------Date: 10/04/2015-----------------------------------------
!
!-------ORIGINAL SCRIP COMMENTS----------------------------------------- 
!     This routine is the driver for computing the addresses and weights 
!     for interpolating between two grids on a sphere.
!
!-----------------------------------------------------------------------
!
!     CVS:$Id: scrip.f,v 1.6 2001/08/21 21:06:44 pwjones Exp $
!
!     Copyright (c) 1997, 1998 the Regents of the University of 
!       California.
!
!     This software and ancillary information (herein called software) 
!     called SCRIP is made available under the terms described here.  
!     The software has been approved for release with associated 
!     LA-CC Number 98-45.
!
!     Unless otherwise indicated, this software has been authored
!     by an employee or employees of the University of California,
!     operator of the Los Alamos National Laboratory under Contract
!     No. W-7405-ENG-36 with the U.S. Department of Energy.  The U.S.
!     Government has rights to use, reproduce, and distribute this
!     software.  The public may copy and use this software without
!     charge, provided that this Notice and any statement of authorship
!     are reproduced on all copies.  Neither the Government nor the
!     University makes any warranty, express or implied, or assumes
!     any liability or responsibility for the use of this software.
!
!     If software is modified to produce derivative works, such modified
!     software should be clearly marked, so as not to confuse it with 
!     the version available from Los Alamos National Laboratory.
!
!***********************************************************************

      module scrip

!-----------------------------------------------------------------------

      use kinds_mod                  ! module defining data types
      use constants                  ! module for common constants
      use iounits                    ! I/O unit manager
      use timers                     ! CPU timers
      use scripwrap_mod              ! 
      use grids                      ! module with grid information
      use remap_vars                 ! common remapping variables
      use remap_conservative         ! routines for conservative remap
      use remap_distance_weight      ! routines for dist-weight remap
      use remap_bilinear             ! routines for bilinear interp
      use remap_bicubic              ! routines for bicubic  interp
      use remap_write                ! routines for remap output

      contains

      subroutine scrip_package( grid1_file, grid2_file,                 &
     &                          grid1_xdim, grid1_ydim,                 &
     &                          grid1_lon_rho, grid1_lat_rho,           &
     &                          grid1_lon_psi, grid1_lat_psi,           &
     &                          grid2_xdim, grid2_ydim,                 &
     &                          grid2_lon_rho, grid2_lat_rho,           &
     &                          grid2_lon_psi, grid2_lat_psi,           &
     &                          interp_file1, interp_file2,             &
     &                          map1_name, map2_name,                   &
     &                          src_mask, dst_mask,                     &
     &                          dst_mask_unlim,                         &
     &                          counter_grid,                           &
     &                          Ngrids_comb_total,                      &
     &                          output_ncfile, MyComm, My_map_type,     &
     &                          samegrid)

      implicit none 

      character (char_len), intent(in) ::                               &
     &           grid1_file,                                            &
                               ! filename of grid file containing grid1
     &           grid2_file,                                            &
                               ! filename of grid file containing grid2
     &           interp_file1,                                          &
                               ! filename for output remap data (map2)
     &           interp_file2,                                          &
                               ! filename for output remap data (map2)
     &           map1_name,                                             &
                               ! name for mapping from grid1 to grid2
     &           map2_name,                                             &
                               ! name for mapping from grid2 to grid1
     &           output_ncfile                                          
                               ! name for output netcdf file

      integer (kind=int_kind), intent(in) :: grid1_xdim, grid1_ydim
      integer (kind=int_kind), intent(in) :: grid2_xdim, grid2_ydim
      integer (kind=int_kind), intent(in) :: src_mask(:,:)
      integer (kind=int_kind), intent(inout) :: dst_mask(:,:)
      integer (kind=int_kind), intent(in) :: dst_mask_unlim(:,:)
      integer (kind=int_kind), intent(in) :: counter_grid
      integer (kind=int_kind), intent(in) :: Ngrids_comb_total
      integer (kind=int_kind), intent(in) :: MyComm, My_map_type
      integer (kind=int_kind), intent(in) :: samegrid
!
      real    (kind=dbl_kind), intent(in) :: grid1_lon_rho(:,:)
      real    (kind=dbl_kind), intent(in) :: grid1_lat_rho(:,:)
      real    (kind=dbl_kind), intent(in) :: grid1_lon_psi(:,:)
      real    (kind=dbl_kind), intent(in) :: grid1_lat_psi(:,:)
      real    (kind=dbl_kind), intent(in) :: grid2_lon_rho(:,:)
      real    (kind=dbl_kind), intent(in) :: grid2_lat_rho(:,:)
      real    (kind=dbl_kind), intent(in) :: grid2_lon_psi(:,:)
      real    (kind=dbl_kind), intent(in) :: grid2_lat_psi(:,:)
!-----------------------------------------------------------------------
!     local variables
!-----------------------------------------------------------------------
      character (char_len) ::                                           &
     &           map_method,                                            &
                               ! choice for mapping method
     &           normalize_opt,                                         &
                               ! option for normalizing weights
     &           output_opt    ! option for output conventions

      integer (kind=int_kind) ::                                        &
     &           nmap         ! number of mappings to compute (1 or 2)

!      namelist /remap_inputs/ grid1_file, grid2_file, 
!     &                        interp_file1, interp_file2,
!     &                        map1_name, map2_name, num_maps,
!     &                        luse_grid1_area, luse_grid2_area,
!     &                        map_method, normalize_opt, output_opt,
!     &                        restrict_type, num_srch_bins
      integer (kind=int_kind) :: n,                                     &
                                        ! dummy counter
     &                           iunit  ! unit number for namelist file
!#ifdef MPI
      integer (kind=int_kind) :: MyError, MyRank
!#endif

#ifdef MPI
      CALL mpi_comm_rank (MyComm, MyRank, MyError)
#else
      MyRank=0
#endif

!-----------------------------------------------------------------------
!
!     initialize timers
!
!-----------------------------------------------------------------------

      call timers_init
      do n=1,max_timers
        call timer_clear(n)
      end do
!
!  These are hardwired for SCRIP-COAWST Package 
!  Now this is selected for each grid pair 18May2020
      num_maps      = 1
      IF (My_map_type.eq.1) THEN
        map_type      = 1
        map_method = 'conservative'
      ELSE IF (My_map_type.eq.2) THEN
        map_type      = 2
        map_method = 'bilinear'
      ELSE IF (My_map_type.eq.3) THEN
        map_type      = 3
        map_method = 'bicubic'
      ELSE IF (My_map_type.eq.4) THEN
        map_type      = 4
        map_method = 'distwgt'
      END IF
      normalize_opt = 'fracarea'
      output_opt = 'scrip'
      restrict_type = 'latitude'
      num_srch_bins = 90
      luse_grid1_area = .false.
      luse_grid2_area = .false.

      select case(map_method)
      case ('conservative')
        map_type = map_type_conserv
        luse_grid_centers = .false.
      case ('bilinear')
        map_type = map_type_bilinear
        luse_grid_centers = .true.
      case ('bicubic')
        map_type = map_type_bicubic
        luse_grid_centers = .true.
      case ('distwgt')
        map_type = map_type_distwgt
        luse_grid_centers = .true.
      case default
        stop 'unknown mapping method'
      end select

      select case(normalize_opt(1:4))
      case ('none')
        norm_opt = norm_opt_none
      case ('frac')
        norm_opt = norm_opt_frcarea
      case ('dest')
        norm_opt = norm_opt_dstarea
      case default
        stop 'unknown normalization option'
      end select

!-----------------------------------------------------------------------
!
!     initialize grid information for both grids
!
!-----------------------------------------------------------------------
      call grid_init_coawst(grid1_file, grid2_file,                     &
     &                      grid1_xdim, grid1_ydim,                     &
     &                      grid1_lon_rho, grid1_lat_rho,               &
     &                      grid1_lon_psi, grid1_lat_psi,               &
     &                      grid2_xdim, grid2_ydim,                     &
     &                      grid2_lon_rho, grid2_lat_rho,               &
     &                      grid2_lon_psi, grid2_lat_psi,               &
     &                      src_mask, dst_mask, MyRank)

#ifdef MPI
!     CALL mpi_comm_rank (MyComm, MyRank, MyError)
      IF (MyRank.eq.0) THEN
#endif
      write(stdout, *) ' Computing remappings between: ',grid1_file
      write(stdout, *) '                          and  ',grid2_file
#ifdef MPI
      END IF
      CALL mpi_barrier (MyComm, MyError)
#endif

!-----------------------------------------------------------------------
!
!     initialize some remapping variables.
!
!-----------------------------------------------------------------------

      call init_remap_vars

!-----------------------------------------------------------------------
!
!     call appropriate interpolation setup routine based on type of
!     remapping requested.
!
!-----------------------------------------------------------------------

      select case(map_type)
      case(map_type_conserv)
        call remap_conserv (MyComm)
      case(map_type_bilinear)
        call remap_bilin (MyComm, samegrid)
      case(map_type_distwgt)
        call remap_distwgt
      case(map_type_bicubic)
        call remap_bicub
      case default
        stop 'Invalid Map Type'
      end select

#ifdef MPI
!     CALL mpi_comm_rank (MyComm, MyRank, MyError)
      IF (MyRank.eq.0) THEN
#endif
!-----------------------------------------------------------------------
!
!     reduce size of remapping arrays and then write remapping info
!     to a file.
!
!-----------------------------------------------------------------------
      if (num_links_map1.gt.0) then
        if (num_links_map1 /= max_links_map1) then
          call resize_remap_vars(1, num_links_map1-max_links_map1)
        endif
        if ((num_maps > 1).and.(num_links_map2 /= max_links_map2)) then
          call resize_remap_vars(2, num_links_map2-max_links_map2)
        endif
      endif

      if (map1_name(1:11).eq.'WRF to ROMS') then
        call check_weights (dst_mask, dst_mask_unlim,                   &
     &                      grid2_lon_rho, grid2_lat_rho, map_method)
      end if

      call write_remap(map1_name, map2_name,                            &
     &                 interp_file1, interp_file2, output_opt,          &
     &                 counter_grid, Ngrids_comb_total, output_ncfile)

#ifdef MPI
      END IF
      CALL mpi_barrier (MyComm, MyError)
#endif
!-----------------------------------------------------------------------
!     DEALLOCATE HERE for SCRIP_COAWST package
#ifdef MPI
!     CALL mpi_comm_rank (MyComm, MyRank, MyError)
      IF (MyRank.eq.0) THEN
#endif
      write(stdout,*) "-------------------------------------------"
      write(stdout,*) "Reached the end of mapping one set of grids"
#ifdef MPI
      END IF
      CALL mpi_barrier (MyComm, MyError)
#endif

!     deallocate arrays from grids.f
      deallocate ( grid1_dims, grid2_dims )
      deallocate ( grid1_area, grid2_area )
      deallocate ( grid1_frac, grid2_frac )
      deallocate ( grid1_mask, grid2_mask)
      deallocate ( grid1_center_lon, grid2_center_lon)
      deallocate ( grid1_center_lat, grid2_center_lat)
      deallocate ( grid1_corner_lon, grid2_corner_lon)
      deallocate ( grid1_corner_lat, grid2_corner_lat)
      deallocate ( grid1_bound_box, grid2_bound_box)
      deallocate( bin_addr1, bin_addr2, bin_lats, bin_lons)
      deallocate( grid1_add_map1, grid2_add_map1)
      deallocate( wts_map1)

      if (map_method.eq.'conservative') then
!       deallocate arrays from remap_conserv.f
        if (num_links_map1.gt.0) then
          deallocate(link_add1,link_add2)
        endif
      endif

#ifdef MPI
!      CALL mpi_comm_rank (MyComm, MyRank, MyError)
      IF (MyRank.eq.0) THEN
#endif
      if (counter_grid.eq.Ngrids_comb_total) then
        write(*,*) 'end of scrip package '
      end if
#ifdef MPI
      END IF
      CALL mpi_barrier (MyComm, MyError)
#endif
      end subroutine scrip_package

!======================================================================                                                                       
      SUBROUTINE check_weights (dst_mask, dst_mask_unlim,               &
     &                          grid2_lon_rho, grid2_lat_rho,           &
     &                          map_method)
!
      implicit none
!
!     impoprted variables.
!
      integer (kind=int_kind), intent(inout) :: dst_mask(:,:)
      integer (kind=int_kind), intent(in) :: dst_mask_unlim(:,:)
      real(dbl_kind), intent(in) :: grid2_lon_rho(:,:)
      real(dbl_kind), intent(in) :: grid2_lat_rho(:,:)
      character (char_len), intent(in) :: map_method
!
!     local variables
!
      integer (int_kind) :: grid2_xdim, grid2_ydim
      integer(int_kind)  :: i, j, mm, nx, ny, add_wts, mxlinks
      integer(int_kind)  :: indx, Ikeep, Jkeep
      real(dbl_kind)     :: dist1, dist_max, dlon
      real(dbl_kind)     :: latrad1, latrad2, dep, dlat
      real(dbl_kind)     :: xx1, yy1, xx2, yy2

      integer(int_kind), allocatable :: iloc(:)
      integer(int_kind), allocatable :: jloc(:)
      integer(int_kind), allocatable :: add_src_address(:)
      integer(int_kind), allocatable :: add_dst_address(:)
      real(dbl_kind), allocatable :: add_remap_matrix(:)
!
      integer(int_kind), allocatable :: add1_tmp(:)
      integer(int_kind), allocatable :: add2_tmp(:)
      real(dbl_kind), allocatable :: wts_tmp(:,:)
!
! set some grid dim sizes
!
      grid2_xdim=grid2_dims(1)
      grid2_ydim=grid2_dims(2)
!
! find out how many locations have dst mask unlimited is 1 
! but dst mask is a 0.
! (this means the atm is land but the ocean thinks it is water)
!
      add_wts=0
      do i=1,grid2_xdim
        do j=1,grid2_ydim
          if ((dst_mask_unlim(i,j)-dst_mask(i,j)).eq.1) then
            add_wts=add_wts+1
          end if
        end do
      end do
!
      if (add_wts.gt.0) then
        allocate (iloc(add_wts))
        allocate (jloc(add_wts))
        allocate (add_src_address(add_wts))
        allocate (add_dst_address(add_wts))
        allocate (add_remap_matrix(add_wts))
!
! now with allocated arrays, redo the same logic to fill
! the i, j, remap, and correct dst_mask
!
        add_wts=0
        do i=1,grid2_xdim
          do j=1,grid2_ydim
            if ((dst_mask_unlim(i,j)-dst_mask(i,j)).eq.1) then
              add_wts=add_wts+1
              iloc(add_wts)=i
              jloc(add_wts)=j
              add_dst_address(add_wts)=(j-1)*grid2_xdim+i
              add_src_address(add_wts)=-9999
              add_remap_matrix(add_wts)=1
            end if
          end do
        end do
!
!  for all the wts locations, need to set remap matrix to 0
!  so it does not take weigths from atm land points.
!
        do nx=1,add_wts
          i=iloc(nx)
          j=jloc(nx)
          indx=(j-1)*grid2_xdim+i
          do mm=1,size(grid2_add_map1)
            if (indx.eq.grid2_add_map1(mm)) then
              wts_map1(1,mm)=0
              wts_map1(2,mm)=0
              wts_map1(3,mm)=0
            end if
          end do
        end do
!
! for all the locations, find the closest dst locations 
! (grid2_add_map1) and use their corresponding src address
! (grid1_add_map1).
!
        do mm=1,add_wts
          nx=iloc(mm)
          ny=jloc(mm)
          dist_max=10e6
          Ikeep=1
          Jkeep=1
          xx2=grid2_lon_rho(nx,ny)
          yy2=grid2_lat_rho(nx,ny)
          do i=1,grid2_xdim
            do j=1,grid2_ydim
              xx1=grid2_lon_rho(i,j)
              yy1=grid2_lat_rho(i,j)
              dlon = xx1-xx2
              latrad1=abs(yy1*deg2rad)
              latrad2=abs(yy2*deg2rad)
              dep=cos(0.5*(latrad2+latrad1))*dlon
              dlat=yy2-yy1
              dist1=1852.0*60.0*sqrt(dlat**2+dep**2)
              if((dist1<=dist_max).and.(dst_mask(i,j).eq.1)) then
                dist_max=dist1
                Ikeep=i
                Jkeep=j
              endif
            enddo
          enddo
          indx=(Jkeep-1)*grid2_xdim+Ikeep
          do i=1,size(grid1_add_map1)
            if (grid2_add_map1(i).eq.indx) then
              add_src_address(mm)=grid1_add_map1(i)
            end if
          end do
        end do
!
! set all the new dst masks to = 1
!
        do mm=1,add_wts
          i=iloc(mm)
          j=jloc(mm)
          dst_mask(i,j)=1
        end do
!
! correct grid2_mask based on dst_mask
!
        mm=0
        do j=1,grid2_ydim
          do i=1,grid2_xdim
            mm=mm+1
            if (dst_mask(i,j).eq.1) then
              grid2_mask(mm)=.true.
            else
              grid2_mask(mm)=.false.
            end if
          end do 
        end do 
!
! here we add the src, dst, and wts to the global arrays
!
        !***
        !*** allocate temporaries to hold original values
        !***

        mxlinks = size(grid1_add_map1)

        allocate (add1_tmp(mxlinks), add2_tmp(mxlinks),                 &
     &            wts_tmp(num_wts,mxlinks))

        add1_tmp = grid1_add_map1
        add2_tmp = grid2_add_map1
        wts_tmp  = wts_map1
        
        !***
        !*** deallocate originals and increment max_links then
        !*** reallocate arrays at new size
        !***

        deallocate (grid1_add_map1, grid2_add_map1, wts_map1)

        num_links_map1 = mxlinks + add_wts
        allocate (grid1_add_map1(num_links_map1),                       &
     &            grid2_add_map1(num_links_map1),                       &
     &            wts_map1(num_wts,num_links_map1))

        !***
        !*** restore original values from temp arrays and
        !*** deallocate temps
        !***
        do mm=1,mxlinks
          grid1_add_map1(mm) = add1_tmp (mm)
          grid2_add_map1(mm) = add2_tmp (mm)
          wts_map1    (:,mm) = wts_tmp(:,mm)
        end do
        i=0
        do mm=mxlinks+1,num_links_map1
          i=i+1
          if ((add_src_address(i).le.0).or.                             &
     &        (add_dst_address(i).le.0)) then
            grid1_add_map1(mm) = 1
            grid2_add_map1(mm) = 1
            wts_map1    (1,mm) = 0
            if (map_method.eq.'conservative') then
              wts_map1    (2,mm) = 0
              wts_map1    (3,mm) = 0
            endif
          else
            grid1_add_map1(mm) = add_src_address(i)
            grid2_add_map1(mm) = add_dst_address(i)
            wts_map1    (1,mm) = add_remap_matrix(i)
            if (map_method.eq.'conservative') then
              wts_map1    (2,mm) = 0
              wts_map1    (3,mm) = 0
            endif
          endif
        end do

        deallocate(add1_tmp, add2_tmp, wts_tmp)
        deallocate (iloc, jloc)
        deallocate (add_src_address, add_dst_address, add_remap_matrix)
      end if

      RETURN
      END SUBROUTINE check_weights

      end module scrip

