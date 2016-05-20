!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Note SCRIP_COAWST required the original file from SCRIP package 
!     to be converted into a module including a subroutine. 
!     -Arrays come in through scrip.f and read here for both grids1&2. 
!        
!---- Written by John C. Warner-----------------------------------------
!-----         Tarandeep S. Kalra --------------------------------------
!--------------Date: 10/04/2015-----------------------------------------
! 
!-----ORIGINAL SCRIP COMMENTS-------------------------------------------
!
!     This module reads in and initializes two grids for remapping.
!     NOTE: grid1 must be the master grid -- the grid that determines
!           which cells participate (e.g. land mask) and the fractional
!           area of grid2 cells that participate in the remapping.
!
!-----------------------------------------------------------------------
!
!     CVS:$Id: grids.f,v 1.6 2001/08/21 21:06:41 pwjones Exp $
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

      module grids

!-----------------------------------------------------------------------

      use kinds_mod    ! defines data types
      use constants    ! common constants
      use iounits      ! I/O unit manager
      use netcdf_mod   ! netCDF stuff
      
      implicit none

!-----------------------------------------------------------------------
!
!     variables that describe each grid
!
!-----------------------------------------------------------------------

      integer (kind=int_kind), save ::
     &             grid1_size, grid2_size, ! total points on each grid
     &             grid1_rank, grid2_rank, ! rank of each grid
     &             grid1_corners, grid2_corners ! number of corners
                                                ! for each grid cell

      integer (kind=int_kind), dimension(:), allocatable, save ::
     &             grid1_dims, grid2_dims  ! size of each grid dimension

      character(char_len), save :: 
     &             grid1_name, grid2_name  ! name for each grid

      character (char_len), save :: 
     &             grid1_units, ! units for grid coords (degs/radians)
     &             grid2_units  ! units for grid coords

      real (kind=dbl_kind), parameter ::
     &      deg2rad = pi/180.   ! conversion for deg to rads

!-----------------------------------------------------------------------
!
!     grid coordinates and masks
!
!-----------------------------------------------------------------------

      logical (kind=log_kind), dimension(:), allocatable, save ::
     &             grid1_mask,        ! flag which cells participate
     &             grid2_mask         ! flag which cells participate

      real (kind=dbl_kind), dimension(:), allocatable, save ::
     &             grid1_center_lat,  ! lat/lon coordinates for
     &             grid1_center_lon,  ! each grid center in radians
     &             grid2_center_lat, 
     &             grid2_center_lon,
     &             grid1_area,        ! tot area of each grid1 cell
     &             grid2_area,        ! tot area of each grid2 cell
     &             grid1_area_in,     ! area of grid1 cell from file
     &             grid2_area_in,     ! area of grid2 cell from file
     &             grid1_frac,        ! fractional area of grid cells
     &             grid2_frac         ! participating in remapping

      real (kind=dbl_kind), dimension(:,:), allocatable, save ::
     &             grid1_corner_lat,  ! lat/lon coordinates for
     &             grid1_corner_lon,  ! each grid corner in radians
     &             grid2_corner_lat, 
     &             grid2_corner_lon

      logical (kind=log_kind), save ::
     &             luse_grid_centers ! use centers for bounding boxes
     &,            luse_grid1_area   ! use area from grid file
     &,            luse_grid2_area   ! use area from grid file

      real (kind=dbl_kind), dimension(:,:), allocatable, save ::
     &             grid1_bound_box,  ! lat/lon bounding box for use
     &             grid2_bound_box   ! in restricting grid searches

!-----------------------------------------------------------------------
!
!     bins for restricting searches
!
!-----------------------------------------------------------------------

      character (char_len), save ::
     &        restrict_type  ! type of bins to use

      integer (kind=int_kind), save ::
     &        num_srch_bins  ! num of bins for restricted srch

      integer (kind=int_kind), dimension(:,:), allocatable, save ::
     &        bin_addr1, ! min,max adds for grid1 cells in this lat bin
     &        bin_addr2  ! min,max adds for grid2 cells in this lat bin

      real(kind=dbl_kind), dimension(:,:), allocatable, save ::
     &        bin_lats   ! min,max latitude for each search bin
     &,       bin_lons   ! min,max longitude for each search bin

!***********************************************************************

      contains

!***********************************************************************
      subroutine grid_init_coawst(grid1_file, grid2_file,               &
     &                            grid1_xdim, grid1_ydim,               &
     &                            grid1_lon_rho, grid1_lat_rho,         &
     &                            grid1_lon_psi, grid1_lat_psi,         &
     &                            grid2_xdim, grid2_ydim,               &
     &                            grid2_lon_rho, grid2_lat_rho,         &
     &                            grid2_lon_psi, grid2_lat_psi,         &
     &                            src_mask, dst_mask)
!     necessary changes (e.g. for 0,2pi longitude range)
!
!-----------------------------------------------------------------------

      character(char_len), intent(in) :: grid1_file, grid2_file

      integer (kind=int_kind), intent(in) :: grid1_xdim, grid1_ydim
      integer (kind=int_kind), intent(in) :: grid2_xdim, grid2_ydim
      integer (kind=int_kind), intent(in) :: src_mask(:,:)
      integer (kind=int_kind), intent(in) :: dst_mask(:,:)
      real    (kind=dbl_kind), intent(in) :: grid1_lon_rho(:,:) 
      real    (kind=dbl_kind), intent(in) :: grid1_lat_rho(:,:) 
      real    (kind=dbl_kind), intent(in) :: grid1_lon_psi(:,:) 
      real    (kind=dbl_kind), intent(in) :: grid1_lat_psi(:,:) 
      real    (kind=dbl_kind), intent(in) :: grid2_lon_rho(:,:) 
      real    (kind=dbl_kind), intent(in) :: grid2_lat_rho(:,:) 
      real    (kind=dbl_kind), intent(in) :: grid2_lon_psi(:,:) 
      real    (kind=dbl_kind), intent(in) :: grid2_lat_psi(:,:) 
!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

      integer (kind=int_kind) :: 
     &  n      ! loop counter
     &, nele   ! element loop counter
     &, iunit  ! unit number for opening files
     &, i,j    ! logical 2d addresses
     &, ip1,jp1
     &, n_add, e_add, ne_add
     &, nx, ny
     &, counter
     &, MP,LP 

      integer (kind=int_kind), dimension(:), allocatable :: 
     &                            imask ! integer mask read from file

      real (kind=dbl_kind) ::
     &  dlat,dlon           !lat/lon intervals for search bins

      real (kind=dbl_kind), dimension(4) ::
     &  tmp_lats, tmp_lons  !temps for computing bounding boxes

!-----------------------------------------------------------------------
!
!     read grid size/name data
!
!-----------------------------------------------------------------------
      real(dbl_kind), allocatable :: lat_rho(:,:), lon_rho(:,:) 
      real(dbl_kind), allocatable :: lat_psi(:,:), lon_psi(:,:) 
 
!     These numbers are same for all 3 grid file types
      luse_grid1_area = .false. ! earlier an input 
      grid1_corners=4
      grid1_rank=2
      grid1_size=grid1_xdim*grid1_ydim 
      write(stdout,*)"----------------------------------------------"
      write(stdout,*)"grid1 dimensions=",grid1_xdim, grid1_ydim
      write(stdout,*)"----------------------------------------------"
!-----------------------------------------------------------------------
!
!     allocate grid coordinates/masks and read data
!
!-----------------------------------------------------------------------
      allocate( grid1_dims(grid1_rank))
      allocate( grid1_area(grid1_size)) 
      allocate( grid1_frac(grid1_size))
      allocate( grid1_mask(grid1_size))
      allocate( imask(grid1_size))
      allocate( grid1_center_lon(grid1_size))
      allocate( grid1_center_lat(grid1_size))
      allocate( grid1_corner_lon(grid1_corners, grid1_size))
      allocate( grid1_corner_lat(grid1_corners, grid1_size))
      allocate( grid1_bound_box(4, grid1_size))


      grid1_dims(1)=grid1_xdim
      grid1_dims(2)=grid1_ydim

      grid1_area = zero
      grid1_frac = zero

!-----------------------------------------------------------------------
!
!     initialize logical mask and convert lat/lon units if required
!
!-----------------------------------------------------------------------

! jcw only do deg2rad if spherical and grid in degreees !!!

      counter=0
      do j=1,grid1_dims(2)
        do i=1,grid1_dims(1)
          counter=counter+1
!
          grid1_center_lon(counter)=grid1_lon_rho(i,j)*deg2rad
          grid1_center_lat(counter)=grid1_lat_rho(i,j)*deg2rad
!
          grid1_corner_lon(1,counter)=grid1_lon_psi(i,j)*deg2rad
          grid1_corner_lat(1,counter)=grid1_lat_psi(i,j)*deg2rad
          grid1_corner_lon(2,counter)=grid1_lon_psi(i+1,j)*deg2rad
          grid1_corner_lat(2,counter)=grid1_lat_psi(i+1,j)*deg2rad
          grid1_corner_lon(3,counter)=grid1_lon_psi(i+1,j+1)*deg2rad
          grid1_corner_lat(3,counter)=grid1_lat_psi(i+1,j+1)*deg2rad
          grid1_corner_lon(4,counter)=grid1_lon_psi(i,j+1)*deg2rad
          grid1_corner_lat(4,counter)=grid1_lat_psi(i,j+1)*deg2rad
!
          imask(counter)=src_mask(i,j)
!          
        end do 
      end do
                             
      where (imask == 1)
        grid1_mask = .true.
      elsewhere
        grid1_mask = .false.
      endwhere
       
      deallocate(imask)

!-----------------------------------------------------------------------
!
!     Read the second grid 
!
!-----------------------------------------------------------------------
!     These numbers are same for all 3 grid file types
      luse_grid2_area = .false. ! earlier an input 
      grid2_corners=4
      grid2_rank=2
      grid2_size=grid2_xdim*grid2_ydim 
      
      write(stdout,*)"----------------------------------------------"
      write(stdout,*)"grid2 dimensions=",grid2_xdim, grid2_ydim
      write(stdout,*)"----------------------------------------------"
!-----------------------------------------------------------------------
!
!     allocate grid coordinates/masks and read data
!
!-----------------------------------------------------------------------
      allocate( grid2_dims(grid2_rank))
      allocate( grid2_area(grid2_size)) 
      allocate( grid2_frac(grid2_size))
      allocate( grid2_mask(grid2_size))
      allocate( imask(grid2_size))
      allocate( grid2_center_lon(grid2_size))
      allocate( grid2_center_lat(grid2_size))
      allocate( grid2_corner_lon(grid2_corners, grid2_size))
      allocate( grid2_corner_lat(grid2_corners, grid2_size))
      allocate( grid2_bound_box(4, grid2_size))

        
      grid2_dims(1)=grid2_xdim
      grid2_dims(2)=grid2_ydim

      grid2_area = zero
      grid2_frac = zero
!-----------------------------------------------------------------------
!
!     initialize logical mask and convert lat/lon units if required
!
!-----------------------------------------------------------------------
      counter=0
      do j=1,grid2_dims(2)
        do i=1,grid2_dims(1)
          counter=counter+1
          grid2_center_lon(counter)=grid2_lon_rho(i,j)*deg2rad
          grid2_center_lat(counter)=grid2_lat_rho(i,j)*deg2rad
!
          grid2_corner_lon(1,counter)=grid2_lon_psi(i,j)*deg2rad
          grid2_corner_lat(1,counter)=grid2_lat_psi(i,j)*deg2rad
          grid2_corner_lon(2,counter)=grid2_lon_psi(i+1,j)*deg2rad
          grid2_corner_lat(2,counter)=grid2_lat_psi(i+1,j)*deg2rad
          grid2_corner_lon(3,counter)=grid2_lon_psi(i+1,j+1)*deg2rad
          grid2_corner_lat(3,counter)=grid2_lat_psi(i+1,j+1)*deg2rad
          grid2_corner_lon(4,counter)=grid2_lon_psi(i,j+1)*deg2rad
          grid2_corner_lat(4,counter)=grid2_lat_psi(i,j+1)*deg2rad
!
          imask(counter)=dst_mask(i,j)
!         
        end do 
      end do 

      where (imask == 1)
        grid2_mask = .true.
      elsewhere
        grid2_mask = .false.
      endwhere

      deallocate(imask)

!-----------------------------------------------------------------------
!
!     convert longitudes to 0,2pi interval
!
!-----------------------------------------------------------------------

      where (grid1_center_lon .gt. pi2)  grid1_center_lon =
     &                                   grid1_center_lon - pi2
      where (grid1_center_lon .lt. zero) grid1_center_lon =
     &                                   grid1_center_lon + pi2
      where (grid2_center_lon .gt. pi2)  grid2_center_lon =
     &                                   grid2_center_lon - pi2
      where (grid2_center_lon .lt. zero) grid2_center_lon =
     &                                   grid2_center_lon + pi2
      where (grid1_corner_lon .gt. pi2)  grid1_corner_lon =
     &                                   grid1_corner_lon - pi2
      where (grid1_corner_lon .lt. zero) grid1_corner_lon =
     &                                   grid1_corner_lon + pi2
      where (grid2_corner_lon .gt. pi2)  grid2_corner_lon =
     &                                   grid2_corner_lon - pi2
      where (grid2_corner_lon .lt. zero) grid2_corner_lon =
     &                                   grid2_corner_lon + pi2

!-----------------------------------------------------------------------
!
!     make sure input latitude range is within the machine values
!     for +/- pi/2 
!
!-----------------------------------------------------------------------

      where (grid1_center_lat >  pih) grid1_center_lat =  pih
      where (grid1_corner_lat >  pih) grid1_corner_lat =  pih
      where (grid1_center_lat < -pih) grid1_center_lat = -pih
      where (grid1_corner_lat < -pih) grid1_corner_lat = -pih

      where (grid2_center_lat >  pih) grid2_center_lat =  pih
      where (grid2_corner_lat >  pih) grid2_corner_lat =  pih
      where (grid2_center_lat < -pih) grid2_center_lat = -pih
      where (grid2_corner_lat < -pih) grid2_corner_lat = -pih

!----------------------------------------------------------------------
!
!     compute bounding boxes for restricting future grid searches
!
!-----------------------------------------------------------------------

      if (.not. luse_grid_centers) then
        grid1_bound_box(1,:) = minval(grid1_corner_lat, DIM=1)
        grid1_bound_box(2,:) = maxval(grid1_corner_lat, DIM=1)
        grid1_bound_box(3,:) = minval(grid1_corner_lon, DIM=1)
        grid1_bound_box(4,:) = maxval(grid1_corner_lon, DIM=1)

        grid2_bound_box(1,:) = minval(grid2_corner_lat, DIM=1)
        grid2_bound_box(2,:) = maxval(grid2_corner_lat, DIM=1)
        grid2_bound_box(3,:) = minval(grid2_corner_lon, DIM=1)
        grid2_bound_box(4,:) = maxval(grid2_corner_lon, DIM=1)
      else

        nx = grid1_dims(1)
        ny = grid1_dims(2)

        do n=1,grid1_size

          !*** find N,S and NE points to this grid point

          j = (n - 1)/nx +1
          i = n - (j-1)*nx

          if (i < nx) then
            ip1 = i + 1
          else
            !*** assume cyclic
            ip1 = 1
            !*** but if it is not, correct
            e_add = (j - 1)*nx + ip1
            if (abs(grid1_center_lat(e_add) - 
     &              grid1_center_lat(n   )) > pih) then
              ip1 = i
            endif
          endif

          if (j < ny) then
            jp1 = j+1
          else
            !*** assume cyclic
            jp1 = 1
            !*** but if it is not, correct
            n_add = (jp1 - 1)*nx + i
            if (abs(grid1_center_lat(n_add) - 
     &              grid1_center_lat(n   )) > pih) then
              jp1 = j
            endif
          endif

          n_add = (jp1 - 1)*nx + i
          e_add = (j - 1)*nx + ip1
          ne_add = (jp1 - 1)*nx + ip1

          !*** find N,S and NE lat/lon coords and check bounding box

          tmp_lats(1) = grid1_center_lat(n)
          tmp_lats(2) = grid1_center_lat(e_add)
          tmp_lats(3) = grid1_center_lat(ne_add)
          tmp_lats(4) = grid1_center_lat(n_add)

          tmp_lons(1) = grid1_center_lon(n)
          tmp_lons(2) = grid1_center_lon(e_add)
          tmp_lons(3) = grid1_center_lon(ne_add)
          tmp_lons(4) = grid1_center_lon(n_add)

          grid1_bound_box(1,n) = minval(tmp_lats)
          grid1_bound_box(2,n) = maxval(tmp_lats)
          grid1_bound_box(3,n) = minval(tmp_lons)
          grid1_bound_box(4,n) = maxval(tmp_lons)
        end do

        nx = grid2_dims(1)
        ny = grid2_dims(2)

        do n=1,grid2_size

          !*** find N,S and NE points to this grid point

          j = (n - 1)/nx +1
          i = n - (j-1)*nx

          if (i < nx) then
            ip1 = i + 1
          else
            !*** assume cyclic
            ip1 = 1
            !*** but if it is not, correct
            e_add = (j - 1)*nx + ip1
            if (abs(grid2_center_lat(e_add) - 
     &              grid2_center_lat(n   )) > pih) then
              ip1 = i
            endif
          endif

          if (j < ny) then
            jp1 = j+1
          else
            !*** assume cyclic
            jp1 = 1
            !*** but if it is not, correct
            n_add = (jp1 - 1)*nx + i
            if (abs(grid2_center_lat(n_add) - 
     &              grid2_center_lat(n   )) > pih) then
              jp1 = j
            endif
          endif

          n_add = (jp1 - 1)*nx + i
          e_add = (j - 1)*nx + ip1
          ne_add = (jp1 - 1)*nx + ip1

          !*** find N,S and NE lat/lon coords and check bounding box

          tmp_lats(1) = grid2_center_lat(n)
          tmp_lats(2) = grid2_center_lat(e_add)
          tmp_lats(3) = grid2_center_lat(ne_add)
          tmp_lats(4) = grid2_center_lat(n_add)

          tmp_lons(1) = grid2_center_lon(n)
          tmp_lons(2) = grid2_center_lon(e_add)
          tmp_lons(3) = grid2_center_lon(ne_add)
          tmp_lons(4) = grid2_center_lon(n_add)

          grid2_bound_box(1,n) = minval(tmp_lats)
          grid2_bound_box(2,n) = maxval(tmp_lats)
          grid2_bound_box(3,n) = minval(tmp_lons)
          grid2_bound_box(4,n) = maxval(tmp_lons)
        end do

      endif

      where (abs(grid1_bound_box(4,:) - grid1_bound_box(3,:)) > pi)
        grid1_bound_box(3,:) = zero
        grid1_bound_box(4,:) = pi2
      end where

      where (abs(grid2_bound_box(4,:) - grid2_bound_box(3,:)) > pi)
        grid2_bound_box(3,:) = zero
        grid2_bound_box(4,:) = pi2
      end where

      !***
      !*** try to check for cells that overlap poles
      !***

      where (grid1_center_lat > grid1_bound_box(2,:))
     &  grid1_bound_box(2,:) = pih

      where (grid1_center_lat < grid1_bound_box(1,:))
     &  grid1_bound_box(1,:) = -pih

      where (grid2_center_lat > grid2_bound_box(2,:))
     &  grid2_bound_box(2,:) = pih

      where (grid2_center_lat < grid2_bound_box(1,:))
     &  grid2_bound_box(1,:) = -pih

!-----------------------------------------------------------------------
!
!     set up and assign address ranges to search bins in order to 
!     further restrict later searches
!
!-----------------------------------------------------------------------

      select case (restrict_type)

      case ('latitude')
        write(stdout,*) 'Using latitude bins to restrict search.'

        allocate(bin_addr1(2,num_srch_bins))
        allocate(bin_addr2(2,num_srch_bins))
        allocate(bin_lats (2,num_srch_bins))
        allocate(bin_lons (2,num_srch_bins))

        dlat = pi/num_srch_bins

        do n=1,num_srch_bins
          bin_lats(1,n) = (n-1)*dlat - pih
          bin_lats(2,n) =     n*dlat - pih
          bin_lons(1,n) = zero
          bin_lons(2,n) = pi2
          bin_addr1(1,n) = grid1_size + 1
          bin_addr1(2,n) = 0
          bin_addr2(1,n) = grid2_size + 1
          bin_addr2(2,n) = 0
        end do

        do nele=1,grid1_size
          do n=1,num_srch_bins
            if (grid1_bound_box(1,nele) <= bin_lats(2,n) .and.
     &          grid1_bound_box(2,nele) >= bin_lats(1,n)) then
              bin_addr1(1,n) = min(nele,bin_addr1(1,n))
              bin_addr1(2,n) = max(nele,bin_addr1(2,n))
            endif
          end do
        end do

        do nele=1,grid2_size
          do n=1,num_srch_bins
            if (grid2_bound_box(1,nele) <= bin_lats(2,n) .and.
     &          grid2_bound_box(2,nele) >= bin_lats(1,n)) then
              bin_addr2(1,n) = min(nele,bin_addr2(1,n))
              bin_addr2(2,n) = max(nele,bin_addr2(2,n))
            endif
          end do
        end do

      case ('latlon')
        write(stdout,*) 'Using lat/lon boxes to restrict search.'

        dlat = pi /num_srch_bins
        dlon = pi2/num_srch_bins

        allocate(bin_addr1(2,num_srch_bins*num_srch_bins))
        allocate(bin_addr2(2,num_srch_bins*num_srch_bins))
        allocate(bin_lats (2,num_srch_bins*num_srch_bins))
        allocate(bin_lons (2,num_srch_bins*num_srch_bins))

        n = 0
        do j=1,num_srch_bins
        do i=1,num_srch_bins
          n = n + 1

          bin_lats(1,n) = (j-1)*dlat - pih
          bin_lats(2,n) =     j*dlat - pih
          bin_lons(1,n) = (i-1)*dlon
          bin_lons(2,n) =     i*dlon
          bin_addr1(1,n) = grid1_size + 1
          bin_addr1(2,n) = 0
          bin_addr2(1,n) = grid2_size + 1
          bin_addr2(2,n) = 0
        end do
        end do

        num_srch_bins = num_srch_bins**2

        do nele=1,grid1_size
          do n=1,num_srch_bins
            if (grid1_bound_box(1,nele) <= bin_lats(2,n) .and.
     &          grid1_bound_box(2,nele) >= bin_lats(1,n) .and.
     &          grid1_bound_box(3,nele) <= bin_lons(2,n) .and.
     &          grid1_bound_box(4,nele) >= bin_lons(1,n)) then
              bin_addr1(1,n) = min(nele,bin_addr1(1,n))
              bin_addr1(2,n) = max(nele,bin_addr1(2,n))
            endif
          end do
        end do

        do nele=1,grid2_size
          do n=1,num_srch_bins
            if (grid2_bound_box(1,nele) <= bin_lats(2,n) .and.
     &          grid2_bound_box(2,nele) >= bin_lats(1,n) .and.
     &          grid2_bound_box(3,nele) <= bin_lons(2,n) .and.
     &          grid2_bound_box(4,nele) >= bin_lons(1,n)) then
              bin_addr2(1,n) = min(nele,bin_addr2(1,n))
              bin_addr2(2,n) = max(nele,bin_addr2(2,n))
            endif
          end do
        end do

      case default
        stop 'unknown search restriction method'
      end select

!-----------------------------------------------------------------------

      end subroutine grid_init_coawst

      end module grids
