!=======================================================================
!     This module defines the global variables for running 
!     scripwrap_mod.f.
!
!     This is a wrapper around SCRIP package developed at USGS, 
!     Woods Hole Science And Marine Center. This routine calls the 
!     SCRIP's package routine that is the driver for computing the 
!     addresses and weights for interpolating between two grids on 
!     a sphere.
!
!---- Written by John C. Warner----------------------------------------- 
!-----         Tarandeep S. Kalra --------------------------------------
!--------------Date: 08/04/2015-----------------------------------------
!=======================================================================

      module scripwrap_mod

      use kinds_mod          ! define common data types
      implicit none
      save

!-----------------------------------------------------------------------
      character(char_len), dimension(5) :: roms_grids
      character(char_len), dimension(5) :: swan_coord, swan_bath
      character(char_len), dimension(5) :: ww3_xcoord, ww3_ycoord
      character(char_len), dimension(5) :: ww3_bath
      character(char_len), dimension(5) :: wrf_grids, hydro_grids
      character(char_len) :: output_ncfile 

      integer(int_kind) :: iunit, counter_grid
      integer(int_kind) :: Ngrids_roms, Ngrids_swan, Ngrids_wrf
      integer(int_kind) :: Ngrids_ww3, Ngrids_hyd
      integer(int_kind) :: Ngrids_comb_total
      integer(int_kind) :: roms_swan_samegrid
      integer(int_kind),dimension(5) :: swan_numx, swan_numy
      integer(int_kind),dimension(5) :: ww3_numx, ww3_numy
      integer(int_kind),dimension(5) :: cartesian
      integer(int_kind),dimension(5) :: parent_grid_ratio
      integer(int_kind),dimension(5) :: parent_id

!-----------------------------------------------------------------------
      type get_swan_grid
        integer(int_kind) :: istr_w, jstr_w,                            &
     &                       iend_w, jend_w, ref_fac
        integer(int_kind) :: Numx_swan, Numy_swan 
        real(dbl_kind), allocatable :: xx(:,:)
! rho points
        real(dbl_kind), allocatable :: lon_rho_w(:,:), lat_rho_w(:,:)
! full grid psi/corner points
        real(dbl_kind), allocatable ::x_full_grid(:,:), y_full_grid(:,:)
! masking values
        integer(int_kind), allocatable :: mask_rho_w(:,:)
        integer(int_kind), allocatable :: src_mask(:,:), dst_mask(:,:)
      end type get_swan_grid

      type(get_swan_grid), allocatable :: ngrd_sw(:)
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      type get_ww3_grid
        integer(int_kind) :: istr_w, jstr_w,                            &
     &                       iend_w, jend_w, ref_fac
        integer(int_kind) :: Numx_ww3, Numy_ww3 
        real(dbl_kind), allocatable :: xx(:,:)
! rho points
        real(dbl_kind), allocatable :: lon_rho_w(:,:), lat_rho_w(:,:)
! full grid psi/corner points
        real(dbl_kind), allocatable ::x_full_grid(:,:), y_full_grid(:,:)
! masking values
        integer(int_kind), allocatable :: mask_rho_w(:,:)
        integer(int_kind), allocatable :: src_mask(:,:), dst_mask(:,:)
      end type get_ww3_grid

      type(get_ww3_grid), allocatable :: ngrd_w3(:)
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      type get_roms_grid
        integer(int_kind) :: istr_o, jstr_o,                            &
     &                       iend_o, jend_o, ref_fac
        integer(int_kind) :: xi_size, eta_size
! rho points
        real(dbl_kind), allocatable :: lon_rho_o(:,:), lat_rho_o(:,:)
        real(dbl_kind), allocatable :: mask_rho_or(:,:)
! full grid psi/corner points
        real(dbl_kind), allocatable ::x_full_grid(:,:), y_full_grid(:,:)
! masking values
        integer(int_kind), allocatable :: mask_rho_o(:,:)
        integer(int_kind), allocatable :: src_mask(:,:), dst_mask(:,:)
      end type get_roms_grid

      type(get_roms_grid), allocatable :: ngrd_rm(:)
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      type get_wrf_grid
        integer(int_kind) :: we_size, sn_size
        integer(int_kind) :: istr_a, jstr_a,                            &
     &                       iend_a, jend_a
! rho points
        real(dbl_kind), allocatable :: lon_rho_a(:,:),lat_rho_a(:,:)
! full grid psi/corner points
        real(dbl_kind), allocatable ::x_full_grid(:,:), y_full_grid(:,:)
! masking values
        integer(int_kind), allocatable :: mask_rho_a(:,:) 
        integer(int_kind), allocatable :: src_mask(:,:), dst_mask(:,:)
      end type get_wrf_grid

      type(get_wrf_grid), allocatable :: ngrd_wr(:)
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      type get_hydro_grid
        integer(int_kind) :: we_size, sn_size
        integer(int_kind) :: istr_h, jstr_h,                            &
     &                       iend_h, jend_h
! rho points
        real(dbl_kind), allocatable :: lon_rho_h(:,:),lat_rho_h(:,:)
! full grid psi/corner points
        real(dbl_kind), allocatable ::x_full_grid(:,:), y_full_grid(:,:)
! masking values
        integer(int_kind), allocatable :: mask_rho_h(:,:) 
        integer(int_kind), allocatable :: src_mask(:,:), dst_mask(:,:)
      end type get_hydro_grid

      type(get_hydro_grid), allocatable :: ngrd_hy(:)
!-----------------------------------------------------------------------

      end module scripwrap_mod

!=======================================================================
