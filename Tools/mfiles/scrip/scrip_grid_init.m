%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%!
%!     This module reads in and initializes two grids for remapping.
%!     NOTE: grid1 must be the master grid -- the grid that determines
%!           which cells participate (e.g. land mask) and the fractional
%!           area of grid2 cells that participate in the remapping.
%!
%!-----------------------------------------------------------------------
%!
%!     CVS:$Id: grids.f,v 1.6 2001/08/21 21:06:41 pwjones Exp $
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

%      module grids

%!-----------------------------------------------------------------------

%      use kinds_mod    ! defines data types
%      use constants    ! common constants
%      use iounits      ! I/O unit manager
%      use netcdf_mod   ! netCDF stuff

%      implicit none

%!-----------------------------------------------------------------------
%!
%!     variables that describe each grid
%!
%!-----------------------------------------------------------------------
%
%      integer (kind=int_kind), save ::
%     &             grid1_size, grid2_size, ! total points on each grid
%     &             grid1_rank, grid2_rank, ! rank of each grid
%     &             grid1_corners, grid2_corners ! number of corners
%                                                ! for each grid cell
%
%      integer (kind=int_kind), dimension(:), allocatable, save ::
%     &             grid1_dims, grid2_dims  ! size of each grid dimension
%
%      character(char_len), save :: 
%     &             grid1_name, grid2_name  ! name for each grid
%
%      character (char_len), save :: 
%     &             grid1_units, ! units for grid coords (degs/radians)
%     &             grid2_units  ! units for grid coords
%
%      real (kind=dbl_kind), parameter ::
%     &      deg2rad = pi/180.   ! conversion for deg to rads
%
%!-----------------------------------------------------------------------
%!
%!     grid coordinates and masks
%!
%!-----------------------------------------------------------------------
%
%      logical (kind=log_kind), dimension(:), allocatable, save ::
%     &             grid1_mask,        ! flag which cells participate
%     &             grid2_mask         ! flag which cells participate
%
%      real (kind=dbl_kind), dimension(:), allocatable, save ::
%     &             grid1_center_lat,  ! lat/lon coordinates for
%     &             grid1_center_lon,  ! each grid center in radians
%     &             grid2_center_lat, 
%     &             grid2_center_lon,
%     &             grid1_area,        ! tot area of each grid1 cell
%     &             grid2_area,        ! tot area of each grid2 cell
%     &             grid1_area_in,     ! area of grid1 cell from file
%     &             grid2_area_in,     ! area of grid2 cell from file
%     &             grid1_frac,        ! fractional area of grid cells
%     &             grid2_frac         ! participating in remapping

%      real (kind=dbl_kind), dimension(:,:), allocatable, save ::
%     &             grid1_corner_lat,  ! lat/lon coordinates for
%     &             grid1_corner_lon,  ! each grid corner in radians
%     &             grid2_corner_lat, 
%     &             grid2_corner_lon
%
%      logical (kind=log_kind), save ::
%     &             luse_grid_centers ! use centers for bounding boxes
%     &,            luse_grid1_area   ! use area from grid file
%     &,            luse_grid2_area   ! use area from grid file
%
%      real (kind=dbl_kind), dimension(:,:), allocatable, save ::
%     &             grid1_bound_box,  ! lat/lon bounding box for use
%     &             grid2_bound_box   ! in restricting grid searches
%
%!-----------------------------------------------------------------------
%!
%!     bins for restricting searches
%!
%!-----------------------------------------------------------------------
%%       character (char_len), save ::
%     &        restrict_type  ! type of bins to use
%
%      integer (kind=int_kind), save ::
%     &        num_srch_bins  ! num of bins for restricted srch
%
%      integer (kind=int_kind), dimension(:,:), allocatable, save ::
%     &        bin_addr1, ! min,max adds for grid1 cells in this lat bin
%     &        bin_addr2  ! min,max adds for grid2 cells in this lat bin
%
%      real(kind=dbl_kind), dimension(:,:), allocatable, save ::
%     &        bin_lats   ! min,max latitude for each search bin
%     &,       bin_lons   ! min,max longitude for each search bin
%
%!***********************************************************************
%
%      contains
%
%!***********************************************************************

%     subroutine grid_init(grid1_file, grid2_file)
      function scrip_grid_init(grid1_file, grid2_file)

%      use grids                      ! module with grid information
       global grid1_size grid2_size grid1_rank grid2_rank
       global grid1_corners grid2_corners grid1_dims grid2_dims
       global grid1_name grid2_name grid1_units grid2_units grid1_mask grid2_mask
       global grid1_center_lat grid1_center_lon grid2_center_lat grid2_center_lon
       global grid1_area grid2_area grid1_area_in grid2_area_in grid1_frac grid2_frac
       global grid1_corner_lat grid1_corner_lon grid2_corner_lat grid2_corner_lon
       global luse_grid_centers luse_grid1_area luse_grid2_area grid1_bound_box grid2_bound_box
       global restrict_type num_srch_bins bin_addr1 bin_addr2 bin_lats bin_lons

%!-----------------------------------------------------------------------
%!
%!     this routine reads grid info from grid files and makes any
%!     necessary changes (e.g. for 0,2pi longitude range)
%!
%!-----------------------------------------------------------------------

%!-----------------------------------------------------------------------
%!
%!     input variables
%!
%!-----------------------------------------------------------------------

 %     character(char_len), intent(in) :: 
 %    &             grid1_file, grid2_file  ! grid data files

%!-----------------------------------------------------------------------
%!
%!     local variables
%!
%!-----------------------------------------------------------------------

%       integer (kind=int_kind) :: 
%      &  n      ! loop counter
%      &, nele   ! element loop counter
%      &, iunit  ! unit number for opening files
%      &, i,j    ! logical 2d addresses
%      &, ip1,jp1
%      &, n_add, e_add, ne_add
%      &, nx, ny
% 
%       integer (kind=int_kind) :: 
%      &         ncstat,           ! netCDF status variable
%      &         nc_grid1_id,       ! netCDF grid file id
%      &         nc_grid2_id,       ! netCDF grid file id
%      &         nc_grid1size_id,   ! netCDF grid size dim id
%      &         nc_grid2size_id,   ! netCDF grid size dim id
%      &         nc_grid1corn_id,   ! netCDF grid corner dim id
%      &         nc_grid2corn_id,   ! netCDF grid corner dim id
%      &         nc_grid1rank_id,   ! netCDF grid rank dim id
%      &         nc_grid2rank_id,   ! netCDF grid rank dim id
%      &         nc_grid1area_id,   ! netCDF grid rank dim id
%      &         nc_grid2area_id,   ! netCDF grid rank dim id
%      &         nc_grid1dims_id,   ! netCDF grid dimension size id
%      &         nc_grid2dims_id,   ! netCDF grid dimension size id
%      &         nc_grd1imask_id,   ! netCDF grid imask var id
%      &         nc_grd2imask_id,   ! netCDF grid imask var id
%      &         nc_grd1crnrlat_id, ! netCDF grid corner lat var id
%      &         nc_grd2crnrlat_id, ! netCDF grid corner lat var id
%      &         nc_grd1crnrlon_id, ! netCDF grid corner lon var id
%      &         nc_grd2crnrlon_id, ! netCDF grid corner lon var id
%      &         nc_grd1cntrlat_id, ! netCDF grid center lat var id
%      &         nc_grd2cntrlat_id, ! netCDF grid center lat var id
%      &         nc_grd1cntrlon_id, ! netCDF grid center lon var id
%      &         nc_grd2cntrlon_id  ! netCDF grid center lon var id
% 
%       integer (kind=int_kind), dimension(:), allocatable :: 
%      &                            imask ! integer mask read from file
% 
%       real (kind=dbl_kind) :: 
%      &  dlat,dlon           ! lat/lon intervals for search bins
% 
%       real (kind=dbl_kind), dimension(4) ::
%      &  tmp_lats, tmp_lons  ! temps for computing bounding boxes
% 
% !-----------------------------------------------------------------------
% !
% !     open grid files and read grid size/name data
% !
% !-----------------------------------------------------------------------

%      ncstat = nf_open(grid1_file, NF_NOWRITE, nc_grid1_id)
%      call netcdf_error_handler(ncstat)
%      ncstat = nf_open(grid2_file, NF_NOWRITE, nc_grid2_id)
%      call netcdf_error_handler(ncstat)
      nc1 = netcdf.open(grid1_file,'NC_NOWRITE');
      nc2 = netcdf.open(grid2_file,'NC_NOWRITE');

%     ncstat = nf_inq_dimid(nc_grid1_id, 'grid_size', nc_grid1size_id)
%     call netcdf_error_handler(ncstat)
%     ncstat = nf_inq_dimlen(nc_grid1_id, nc_grid1size_id, grid1_size)
%     call netcdf_error_handler(ncstat)
      dimid = netcdf.inqDimID(nc1,'grid_size');
      [dimname, grid1_size] = netcdf.inqDim(nc1,dimid);

%     ncstat = nf_inq_dimid(nc_grid2_id, 'grid_size', nc_grid2size_id)
%     call netcdf_error_handler(ncstat)
%     ncstat = nf_inq_dimlen(nc_grid2_id, nc_grid2size_id, grid2_size)
%     call netcdf_error_handler(ncstat)
      dimid = netcdf.inqDimID(nc2,'grid_size');
      [dimname, grid2_size] = netcdf.inqDim(nc2,dimid);

%     ncstat = nf_inq_dimid(nc_grid1_id, 'grid_rank', nc_grid1rank_id)
%     call netcdf_error_handler(ncstat)
%     ncstat = nf_inq_dimlen(nc_grid1_id, nc_grid1rank_id, grid1_rank)
%     call netcdf_error_handler(ncstat)
      dimid = netcdf.inqDimID(nc1,'grid_rank');
      [dimname, grid1_rank] = netcdf.inqDim(nc1,dimid);

%     ncstat = nf_inq_dimid(nc_grid2_id, 'grid_rank', nc_grid2rank_id)
%     call netcdf_error_handler(ncstat)
%     ncstat = nf_inq_dimlen(nc_grid2_id, nc_grid2rank_id, grid2_rank)
%     call netcdf_error_handler(ncstat)
      dimid = netcdf.inqDimID(nc2,'grid_rank');
      [dimname, grid2_rank] = netcdf.inqDim(nc2,dimid);

%     ncstat = nf_inq_dimid(nc_grid1_id,'grid_corners',nc_grid1corn_id)
%     call netcdf_error_handler(ncstat)
%     ncstat = nf_inq_dimlen(nc_grid1_id,nc_grid1corn_id,grid1_corners)
%     call netcdf_error_handler(ncstat)
      dimid = netcdf.inqDimID(nc1,'grid_corners');
      [dimname, grid1_corners] = netcdf.inqDim(nc1,dimid);

%     ncstat = nf_inq_dimid(nc_grid2_id,'grid_corners',nc_grid2corn_id)
%     call netcdf_error_handler(ncstat)
%     ncstat = nf_inq_dimlen(nc_grid2_id,nc_grid2corn_id,grid2_corners)
%     call netcdf_error_handler(ncstat)
      dimid = netcdf.inqDimID(nc2,'grid_corners');
      [dimname, grid2_corners] = netcdf.inqDim(nc2,dimid);

%      allocate( grid1_dims(grid1_rank),
%     &          grid2_dims(grid2_rank))

%      ncstat = nf_get_att_text(nc_grid1_id, nf_global, 'title',
%     &                         grid1_name)
%      call netcdf_error_handler(ncstat)

%      ncstat = nf_get_att_text(nc_grid2_id, nf_global, 'title',
%     &                         grid2_name)
%      call netcdf_error_handler(ncstat)

%!-----------------------------------------------------------------------
%!
%!     allocate grid coordinates/masks and read data
%!
%!-----------------------------------------------------------------------

%       allocate( grid1_mask      (grid1_size),
%      &          grid2_mask      (grid2_size),
%      &          grid1_center_lat(grid1_size), 
%      &          grid1_center_lon(grid1_size),
%      &          grid2_center_lat(grid2_size), 
%      &          grid2_center_lon(grid2_size),
%      &          grid1_area      (grid1_size),
%      &          grid2_area      (grid2_size),
%      &          grid1_frac      (grid1_size),
%      &          grid2_frac      (grid2_size),
%      &          grid1_corner_lat(grid1_corners, grid1_size),
%      &          grid1_corner_lon(grid1_corners, grid1_size),
%      &          grid2_corner_lat(grid2_corners, grid2_size),
%      &          grid2_corner_lon(grid2_corners, grid2_size),
%      &          grid1_bound_box (4            , grid1_size),
%      &          grid2_bound_box (4            , grid2_size))
% 
%       allocate(imask(grid1_size))


        grid1_mask=zeros(grid1_size,1);
        grid2_mask=zeros(grid2_size,1);
        grid1_center_lat=zeros(grid1_size,1);
        grid1_center_lon=zeros(grid1_size,1);
        grid2_center_lat=zeros(grid2_size,1);
        grid2_center_lon=zeros(grid2_size,1);
        grid1_area=zeros(grid1_size,1);
        grid2_area=zeros(grid2_size,1);
        grid1_frac=zeros(grid1_size,1);
        grid2_frac=zeros(grid2_size,1);
        grid1_corner_lat=zeros(grid1_corners, grid1_size);
        grid1_corner_lon=zeros(grid1_corners, grid1_size);
        grid2_corner_lat=zeros(grid2_corners, grid2_size);
        grid2_corner_lon=zeros(grid2_corners, grid2_size);
        grid1_bound_box=zeros (4            , grid1_size);
        grid2_bound_box=zeros (4            , grid2_size);
        imask=zeros(grid1_size,1);

%       ncstat = nf_inq_varid(nc_grid1_id, 'grid_dims', nc_grid1dims_id)
%       call netcdf_error_handler(ncstat)
%       ncstat = nf_inq_varid(nc_grid1_id, 'grid_imask', nc_grd1imask_id)
%       call netcdf_error_handler(ncstat)
%       ncstat = nf_inq_varid(nc_grid1_id, 'grid_center_lat', 
%      &                                   nc_grd1cntrlat_id)
%       call netcdf_error_handler(ncstat)
%       ncstat = nf_inq_varid(nc_grid1_id, 'grid_center_lon', 
%      &                                   nc_grd1cntrlon_id)
%       call netcdf_error_handler(ncstat)
%       ncstat = nf_inq_varid(nc_grid1_id, 'grid_corner_lat', 
%      &                                   nc_grd1crnrlat_id)
%       call netcdf_error_handler(ncstat)
%       ncstat = nf_inq_varid(nc_grid1_id, 'grid_corner_lon', 
%      &                                   nc_grd1crnrlon_id)
%       call netcdf_error_handler(ncstat)

%     ncstat = nf_get_var_int(nc_grid1_id, nc_grid1dims_id, grid1_dims)
%     call netcdf_error_handler(ncstat)
      grid1_dims=ncread(grid1_file,'grid_dims');

%     ncstat = nf_get_var_int(nc_grid1_id, nc_grd1imask_id, imask)
%     call netcdf_error_handler(ncstat)
      imask=ncread(grid1_file,'grid_imask');

%      ncstat = nf_get_var_double(nc_grid1_id, nc_grd1cntrlat_id, 
%     &                                       grid1_center_lat)
%      call netcdf_error_handler(ncstat)
      grid1_center_lat=ncread(grid1_file,'grid_center_lat');

%      ncstat = nf_get_var_double(nc_grid1_id, nc_grd1cntrlon_id, 
%     &                                       grid1_center_lon)
%      call netcdf_error_handler(ncstat)
      grid1_center_lon=ncread(grid1_file,'grid_center_lon');

%      ncstat = nf_get_var_double(nc_grid1_id, nc_grd1crnrlat_id, 
%     &                                       grid1_corner_lat)
%      call netcdf_error_handler(ncstat)
      grid1_corner_lat=ncread(grid1_file,'grid_corner_lat');

%      ncstat = nf_get_var_double(nc_grid1_id, nc_grd1crnrlon_id, 
%     &                                       grid1_corner_lon)
%      call netcdf_error_handler(ncstat)
      grid1_corner_lon=ncread(grid1_file,'grid_corner_lon');

      if (luse_grid1_area)
%        allocate (grid1_area_in(grid1_size))
%        ncstat = nf_inq_varid(nc_grid1_id, 'grid_area', nc_grid1area_id)
%        call netcdf_error_handler(ncstat)
%        ncstat = nf_get_var_double(nc_grid1_id, nc_grid1area_id, 
%     &                                          grid1_area_in)
%        call netcdf_error_handler(ncstat)
        grid1_area_in=ncread(grid1_file,'grid_area');
      end

%      grid1_area = zero
%      grid1_frac = zero

%!-----------------------------------------------------------------------
%!
%!     initialize logical mask and convert lat/lon units if required
%!
%!-----------------------------------------------------------------------

%     where (imask == 1)
%       grid1_mask = .true.
%     elsewhere
%       grid1_mask = .false.
%     endwhere
%     deallocate(imask)
      grid1_mask=logical(imask);
      clear imask

%      grid1_units = ' ';
%      ncstat = nf_get_att_text(nc_grid1_id, nc_grd1cntrlat_id, 'units',
%     &                         grid1_units)
%      call netcdf_error_handler(ncstat)
      grid1_units=ncreadatt(grid1_file,'grid_center_lat','units');

      switch (grid1_units(1:7))
      case ('degrees')

        grid1_center_lat = grid1_center_lat*deg2rad
        grid1_center_lon = grid1_center_lon*deg2rad

      case ('radians')

%        !*** no conversion necessary

      case default

        disp('unknown units supplied for grid1 center lat/lon: ');
        disp('proceeding assuming radians');

      end

%      grid1_units = ' '
%      ncstat = nf_get_att_text(nc_grid1_id, nc_grd1crnrlat_id, 'units',
%     &                         grid1_units)
%      call netcdf_error_handler(ncstat)
      grid1_units=ncreadatt(grid1_file,'grid_corner_lat','units');

      switch (grid1_units(1:7))
      case ('degrees')

        grid1_corner_lat = grid1_corner_lat*deg2rad
        grid1_corner_lon = grid1_corner_lon*deg2rad

      case ('radians')

%        !*** no conversion necessary

      case default

        disp('unknown units supplied for grid1 corner lat/lon: ');
        disp('proceeding assuming radians');

      end

%     ncstat = nf_close(nc_grid1_id)
%     call netcdf_error_handler(ncstat)
      netcdf.close(nc1);
%!-----------------------------------------------------------------------
%!
%!     read data for grid 2
%!
%!-----------------------------------------------------------------------

%      allocate(imask(grid2_size))

%      ncstat = nf_inq_varid(nc_grid2_id, 'grid_dims', nc_grid2dims_id)
%      call netcdf_error_handler(ncstat)
%      ncstat = nf_inq_varid(nc_grid2_id, 'grid_imask', nc_grd2imask_id)
%      call netcdf_error_handler(ncstat)
%      ncstat = nf_inq_varid(nc_grid2_id, 'grid_center_lat', 
%     &                                   nc_grd2cntrlat_id)
%      call netcdf_error_handler(ncstat)
%      ncstat = nf_inq_varid(nc_grid2_id, 'grid_center_lon', 
%     &                                   nc_grd2cntrlon_id)
%      call netcdf_error_handler(ncstat)
%      ncstat = nf_inq_varid(nc_grid2_id, 'grid_corner_lat', 
%     &                                   nc_grd2crnrlat_id)
%      call netcdf_error_handler(ncstat)
%      ncstat = nf_inq_varid(nc_grid2_id, 'grid_corner_lon', 
%     &                                   nc_grd2crnrlon_id)
%      call netcdf_error_handler(ncstat)

%     ncstat = nf_get_var_int(nc_grid2_id, nc_grid2dims_id, grid2_dims)
%     call netcdf_error_handler(ncstat)
      grid2_dims=ncread(grid2_file,'grid_dims');

 %     ncstat = nf_get_var_int(nc_grid2_id, nc_grd2imask_id, imask)
 %     call netcdf_error_handler(ncstat)
      imask=ncread(grid2_file,'grid_imask');

%      ncstat = nf_get_var_double(nc_grid2_id, nc_grd2cntrlat_id, 
%     &                                       grid2_center_lat)
%      call netcdf_error_handler(ncstat)
      grid2_center_lat=ncread(grid2_file,'grid_center_lat');

%      ncstat = nf_get_var_double(nc_grid2_id, nc_grd2cntrlon_id, 
%     &                                       grid2_center_lon)
%      call netcdf_error_handler(ncstat)
      grid2_center_lon=ncread(grid2_file,'grid_center_lon');

%      ncstat = nf_get_var_double(nc_grid2_id, nc_grd2crnrlat_id, 
%     &                                       grid2_corner_lat)
%      call netcdf_error_handler(ncstat)
      grid2_corner_lat=ncread(grid2_file,'grid_corner_lat');

%      ncstat = nf_get_var_double(nc_grid2_id, nc_grd2crnrlon_id, 
%     &                                       grid2_corner_lon)
%      call netcdf_error_handler(ncstat)
      grid2_corner_lon=ncread(grid2_file,'grid_corner_lon');

      if (luse_grid2_area) then
%        allocate (grid2_area_in(grid2_size))
%        ncstat = nf_inq_varid(nc_grid2_id, 'grid_area', nc_grid2area_id)
%        call netcdf_error_handler(ncstat)
%        ncstat = nf_get_var_double(nc_grid2_id, nc_grid2area_id, 
%     &                                          grid2_area_in)
%        call netcdf_error_handler(ncstat)
        grid2_area_in=ncread(grid2_file,'grid_area');
      end

%      grid2_area = zero
%      grid2_frac = zero

%!-----------------------------------------------------------------------
%!
%!     initialize logical mask and convert lat/lon units if required
%!
%!-----------------------------------------------------------------------

%      where (imask == 1)
%        grid2_mask = .true.
%      elsewhere
%        grid2_mask = .false.
%      endwhere
%      deallocate(imask)
      grid2_mask=logical(imask);
      clear imask

%      grid2_units = ' '
%      ncstat = nf_get_att_text(nc_grid2_id, nc_grd2cntrlat_id, 'units',
%     &                         grid2_units)
%      call netcdf_error_handler(ncstat)
      grid2_units=ncreadatt(grid2_file,'grid_center_lat','units');

      switch (grid2_units(1:7))
      case ('degrees')

        grid2_center_lat = grid2_center_lat*deg2rad
        grid2_center_lon = grid2_center_lon*deg2rad

      case ('radians')

%        !*** no conversion necessary

      case default

        disp('unknown units supplied for grid2 center lat/lon: ');
        disp('proceeding assuming radians');

      end

%      grid2_units = ' '
%      ncstat = nf_get_att_text(nc_grid2_id, nc_grd2crnrlat_id, 'units',
%     &                         grid2_units)
%      call netcdf_error_handler(ncstat)
      grid2_units=ncreadatt(grid2_file,'grid_corner_lat','units');

      switch (grid2_units(1:7))
      case ('degrees')

        grid2_corner_lat = grid2_corner_lat*deg2rad
        grid2_corner_lon = grid2_corner_lon*deg2rad

      case ('radians')

%        !*** no conversion necessary

      case default

        disp('no units supplied for grid2 corner lat/lon: ');
        disp('proceeding assuming radians');

      end

%      ncstat = nf_close(nc_grid2_id)
%      call netcdf_error_handler(ncstat)
      netcdf.close(nc2);


%!-----------------------------------------------------------------------
%!
%!     convert longitudes to 0,2pi interval
%!
%!-----------------------------------------------------------------------

%      where (grid1_center_lon .gt. pi2)  grid1_center_lon =
%     &                                   grid1_center_lon - pi2
%      where (grid1_center_lon .lt. zero) grid1_center_lon =
%     &                                   grid1_center_lon + pi2
%      where (grid2_center_lon .gt. pi2)  grid2_center_lon =
%     &                                   grid2_center_lon - pi2
%      where (grid2_center_lon .lt. zero) grid2_center_lon =
%     &                                   grid2_center_lon + pi2
%      where (grid1_corner_lon .gt. pi2)  grid1_corner_lon =
%     &                                   grid1_corner_lon - pi2
%      where (grid1_corner_lon .lt. zero) grid1_corner_lon =
%     &                                   grid1_corner_lon + pi2
%      where (grid2_corner_lon .gt. pi2)  grid2_corner_lon =
%     &                                   grid2_corner_lon - pi2
%      where (grid2_corner_lon .lt. zero) grid2_corner_lon =
%     &                                   grid2_corner_lon + pi2
     grid1_center_lon(grid1_center_lon>2*pi)=grid1_center_lon(grid1_center_lon>2*pi)-2*pi;
     grid1_center_lon(grid1_center_lon<0)   =grid1_center_lon(grid1_center_lon<0)+2*pi;
     grid2_center_lon(grid2_center_lon>2*pi)=grid2_center_lon(grid2_center_lon>2*pi)-2*pi;
     grid2_center_lon(grid2_center_lon<0)   =grid2_center_lon(grid2_center_lon<0)+2*pi;

     grid1_corner_lon(grid1_corner_lon>2*pi)=grid1_corner_lon(grid1_corner_lon>2*pi)-2*pi;
     grid1_corner_lon(grid1_corner_lon<0)   =grid1_corner_lon(grid1_corner_lon<0)+2*pi;
     grid2_corner_lon(grid2_corner_lon>2*pi)=grid2_corner_lon(grid2_corner_lon>2*pi)-2*pi;
     grid2_corner_lon(grid2_corner_lon<0)   =grid2_corner_lon(grid2_corner_lon<0)+2*pi;


%!-----------------------------------------------------------------------
%!
%!     make sure input latitude range is within the machine values
%!     for +/- pi/2 
%!
%!-----------------------------------------------------------------------

%      where (grid1_center_lat >  pih) grid1_center_lat =  pih
%      where (grid1_corner_lat >  pih) grid1_corner_lat =  pih
%      where (grid1_center_lat < -pih) grid1_center_lat = -pih
%      where (grid1_corner_lat < -pih) grid1_corner_lat = -pih

%      where (grid2_center_lat >  pih) grid2_center_lat =  pih
%      where (grid2_corner_lat >  pih) grid2_corner_lat =  pih
%      where (grid2_center_lat < -pih) grid2_center_lat = -pih
%      where (grid2_corner_lat < -pih) grid2_corner_lat = -pih

     grid1_center_lat(grid1_center_lat>pi/2)=pi/2;
     grid1_corner_lat(grid1_corner_lat>pi/2)=pi/2;
     grid1_center_lat(grid1_center_lat<-pi/2)=-pi/2;
     grid1_corner_lat(grid1_corner_lat<-pi/2)=-pi/2;

     grid2_center_lat(grid2_center_lat>pi/2)=pi/2;
     grid2_corner_lat(grid2_corner_lat>pi/2)=pi/2;
     grid2_center_lat(grid2_center_lat<-pi/2)=-pi/2;
     grid2_corner_lat(grid2_corner_lat<-pi/2)=-pi/2;


%!-----------------------------------------------------------------------
%!
%!     compute bounding boxes for restricting future grid searches
%!
%!-----------------------------------------------------------------------

      if (~luse_grid_centers)
%        grid1_bound_box(1,:) = minval(grid1_corner_lat, DIM=1)
%        grid1_bound_box(2,:) = maxval(grid1_corner_lat, DIM=1)
%        grid1_bound_box(3,:) = minval(grid1_corner_lon, DIM=1)
%        grid1_bound_box(4,:) = maxval(grid1_corner_lon, DIM=1)
%        grid2_bound_box(1,:) = minval(grid2_corner_lat, DIM=1)
%        grid2_bound_box(2,:) = maxval(grid2_corner_lat, DIM=1)
%        grid2_bound_box(3,:) = minval(grid2_corner_lon, DIM=1)
%        grid2_bound_box(4,:) = maxval(grid2_corner_lon, DIM=1)

        grid1_bound_box(1,:) = min(grid1_corner_lat(1:4,:));
        grid1_bound_box(2,:) = max(grid1_corner_lat(1:4,:));
        grid1_bound_box(3,:) = min(grid1_corner_lon(1:4,:));
        grid1_bound_box(4,:) = max(grid1_corner_lon(1:4,:));

        grid2_bound_box(1,:) = min(grid2_corner_lat(1:4,:));
        grid2_bound_box(2,:) = max(grid2_corner_lat(1:4,:));
        grid2_bound_box(3,:) = min(grid2_corner_lon(1:4,:));
        grid2_bound_box(4,:) = max(grid2_corner_lon(1:4,:));

      else

        nx = grid1_dims(1);
        ny = grid1_dims(2);

        for n=1,grid1_size

%          !*** find N,S and NE points to this grid point

          j = (n - 1)/nx +1
          i = n - (j-1)*nx

          if (i < nx) then
            ip1 = i + 1
          else
%            !*** assume cyclic
            ip1 = 1
%            !*** but if it is not, correct
            e_add = (j - 1)*nx + ip1
            if (abs(grid1_center_lat(e_add) - ...
                   grid1_center_lat(n   )) > pih) then
              ip1 = i
            end
          end

          if (j < ny) then
            jp1 = j+1
          else
%            !*** assume cyclic
            jp1 = 1
%            !*** but if it is not, correct
            n_add = (jp1 - 1)*nx + i
            if (abs(grid1_center_lat(n_add) -  ...
                   grid1_center_lat(n   )) > pih) then
              jp1 = j
            end
          end

          n_add = (jp1 - 1)*nx + i
          e_add = (j - 1)*nx + ip1
          ne_add = (jp1 - 1)*nx + ip1

 %         !*** find N,S and NE lat/lon coords and check bounding box

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
        end

        nx = grid2_dims(1)
        ny = grid2_dims(2)

        for n=1,grid2_size

%          !*** find N,S and NE points to this grid point

          j = (n - 1)/nx +1
          i = n - (j-1)*nx

          if (i < nx) then
            ip1 = i + 1
          else
%            !*** assume cyclic
            ip1 = 1
%            !*** but if it is not, correct
            e_add = (j - 1)*nx + ip1
            if (abs(grid2_center_lat(e_add) - ...
                   grid2_center_lat(n   )) > pih) then
              ip1 = i
            end
          end

          if (j < ny) then
            jp1 = j+1
          else
%            !*** assume cyclic
            jp1 = 1
%            !*** but if it is not, correct
            n_add = (jp1 - 1)*nx + i
            if (abs(grid2_center_lat(n_add) - ...
                   grid2_center_lat(n   )) > pih) then
              jp1 = j
            end
          end

          n_add = (jp1 - 1)*nx + i
          e_add = (j - 1)*nx + ip1
          ne_add = (jp1 - 1)*nx + ip1

%          !*** find N,S and NE lat/lon coords and check bounding box

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
        end

      end

%      where (abs(grid1_bound_box(4,:) - grid1_bound_box(3,:)) > pi)
%        grid1_bound_box(3,:) = zero
%        grid1_bound_box(4,:) = pi2
%      end where
      zz=abs(grid1_bound_box(4,:) - grid1_bound_box(3,:));
      grid1_bound_box(3,zz>pi) = 0;
      grid1_bound_box(4,zz>pi) = 2*pi;
      
%      where (abs(grid2_bound_box(4,:) - grid2_bound_box(3,:)) > pi)
%        grid2_bound_box(3,:) = zero
%        grid2_bound_box(4,:) = pi2
%      end where
      zz=abs(grid2_bound_box(4,:) - grid2_bound_box(3,:));
      grid2_bound_box(3,zz>pi) = 0;
      grid2_bound_box(4,zz>pi) = 2*pi;

%      !***
%      !*** try to check for cells that overlap poles
%      !***

%      where (grid1_center_lat > grid1_bound_box(2,:))
%     &  grid1_bound_box(2,:) = pih
      zz=(grid1_center_lat > grid1_bound_box(2,:).');
      grid1_bound_box(2,zz) = pi/2;

%      where (grid1_center_lat < grid1_bound_box(1,:))
%     &  grid1_bound_box(1,:) = -pih
      zz=(grid1_center_lat < grid1_bound_box(1,:).');
      grid1_bound_box(1,zz) = -pi/2;

%      where (grid2_center_lat > grid2_bound_box(2,:))
%     &  grid2_bound_box(2,:) = pih
      zz=(grid2_center_lat > grid2_bound_box(2,:).');
      grid2_bound_box(2,zz) = pi/2;

%      where (grid2_center_lat < grid2_bound_box(1,:))
%     &  grid2_bound_box(1,:) = -pih
      zz=(grid2_center_lat < grid2_bound_box(1,:).');
      grid2_bound_box(1,zz) = -pi/2;

%!-----------------------------------------------------------------------
%!
%!     set up and assign address ranges to search bins in order to 
%!     further restrict later searches
%!
%!-----------------------------------------------------------------------

      switch (restrict_type)

      case ('latitude')
        disp('Using latitude bins to restrict search.');

        bin_addr1=zeros(2,num_srch_bins);
        bin_addr2=zeros(2,num_srch_bins);
        bin_lats=zeros(2,num_srch_bins);
        bin_lons=zeros(2,num_srch_bins);

        dlat = pi/num_srch_bins;

        for n=1:num_srch_bins
          bin_lats(1,n) = (n-1)*dlat - pi/2;
          bin_lats(2,n) =     n*dlat - pi/2;
          bin_lons(1,n) = 0;
          bin_lons(2,n) = 2*pi;
          bin_addr1(1,n) = grid1_size + 1;
          bin_addr1(2,n) = 0;
          bin_addr2(1,n) = grid2_size + 1;
          bin_addr2(2,n) = 0;
        end

        for nele=1:grid1_size
          for n=1:num_srch_bins
            if (grid1_bound_box(1,nele) <= bin_lats(2,n) & ...
               grid1_bound_box(2,nele) >= bin_lats(1,n))
              bin_addr1(1,n) = min(nele,bin_addr1(1,n));
              bin_addr1(2,n) = max(nele,bin_addr1(2,n));
            end
          end
        end

        for nele=1:grid2_size
          for n=1:num_srch_bins
            if (grid2_bound_box(1,nele) <= bin_lats(2,n) & ...
               grid2_bound_box(2,nele) >= bin_lats(1,n))
              bin_addr2(1,n) = min(nele,bin_addr2(1,n));
              bin_addr2(2,n) = max(nele,bin_addr2(2,n));
            end
          end
        end

      case ('latlon')
        disp('Using lat/lon boxes to restrict search.')

        dlat = pi /num_srch_bins
        dlon = pi2/num_srch_bins

        allocate(bin_addr1(2,num_srch_bins*num_srch_bins))
        allocate(bin_addr2(2,num_srch_bins*num_srch_bins))
        allocate(bin_lats (2,num_srch_bins*num_srch_bins))
        allocate(bin_lons (2,num_srch_bins*num_srch_bins))

        n = 0
        for j=1,num_srch_bins
        for i=1,num_srch_bins
          n = n + 1

          bin_lats(1,n) = (j-1)*dlat - pih
          bin_lats(2,n) =     j*dlat - pih
          bin_lons(1,n) = (i-1)*dlon
          bin_lons(2,n) =     i*dlon
          bin_addr1(1,n) = grid1_size + 1
          bin_addr1(2,n) = 0
          bin_addr2(1,n) = grid2_size + 1
          bin_addr2(2,n) = 0
        end
        end

        num_srch_bins = num_srch_bins^2

        for nele=1,grid1_size
          for n=1,num_srch_bins
            if (grid1_bound_box(1,nele) <= bin_lats(2,n) & ...
               grid1_bound_box(2,nele) >= bin_lats(1,n) & ...
               grid1_bound_box(3,nele) <= bin_lons(2,n) & ...
               grid1_bound_box(4,nele) >= bin_lons(1,n))
              bin_addr1(1,n) = min(nele,bin_addr1(1,n));
              bin_addr1(2,n) = max(nele,bin_addr1(2,n));
            end
          end
        end

        for nele=1,grid2_size
          for n=1,num_srch_bins
            if (grid2_bound_box(1,nele) <= bin_lats(2,n) & ...
               grid2_bound_box(2,nele) >= bin_lats(1,n) & ...
               grid2_bound_box(3,nele) <= bin_lons(2,n) & ...
               grid2_bound_box(4,nele) >= bin_lons(1,n))
              bin_addr2(1,n) = min(nele,bin_addr2(1,n));
              bin_addr2(2,n) = max(nele,bin_addr2(2,n));
            end
          end
        end

      case ('default')
        disp('unknown search restriction method')
      end

%!-----------------------------------------------------------------------

      end %subroutine grid_init

%!***********************************************************************

%      end %module grids

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

