!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     This module assigns the masking for each model and calls
!     the SCRIP Package routines
!
!-----------------------------------------------------------------------
!     This is a wrapper around SCRIP package developed at USGS,
!     Woods Hole Science And Marine Center. This routine calls the
!     SCRIP's package routine that is the driver for computing the
!     addresses and weights for interpolating between two grids on
!     a sphere.
!
!---- Written by John C. Warner-----------------------------------------
!-----           Tarandeep S. Kalra ------------------------------------
!--------------  Date: 08/20/2015---------------------------------------
!***********************************************************************

      module create_masks

!-----------------------------------------------------------------------
      use kinds_mod     ! defines common data types
      use scripwrap_mod ! defines input file variable names
      use scrip         ! calls the scrip package subroutine

      implicit none

!     Variables that are used in all the subroutines below
      character(char_len) :: grid1_file, grid2_file
      character(char_len) :: interp_file1, interp_file2
      character(char_len) :: map1_name, map2_name
      character(char_len) :: mo_string, mw_string, ma_string
      integer(int_kind)   :: ncstat, nc_file_id, offset
      integer(int_kind)   :: i, j, c, Ikeep, Jkeep
      integer(int_kind)   :: mo, ma, mw, nx, ny, grdsize
      integer(int_kind)   :: N, INOUT, count
      real(dbl_kind)      :: dist1, dist_max, dlon
      real(dbl_kind)      :: latrad1, latrad2, dep, dlat
      real(dbl_kind)      :: xx1, yy1, xx2, yy2, PX, PY
      real(dbl_kind)      :: lons, lats, lond, latd, tol, cff
      integer(int_kind), allocatable :: src_mask_unlim(:,:)
      integer(int_kind), allocatable :: dst_mask_unlim(:,:)
      integer(int_kind), allocatable :: do_adjust(:)

      contains
!=======================================================================
      subroutine ocn2wav_mask()

      implicit none

!     local variables
      real(dbl_kind), allocatable :: XX(:), YY(:)

!  Allocate the source mask
      do mo=1,Ngrids_roms
        nx=ngrd_rm(mo)%xi_size
        ny=ngrd_rm(mo)%eta_size
        allocate(ngrd_rm(mo)%src_mask(nx,ny))
      end do
!  Save the ocean grid mask to src_mask
      do mo=1,Ngrids_roms
        do j=1,ngrd_rm(mo)%eta_size
          do i=1,ngrd_rm(mo)%xi_size
            ngrd_rm(mo)%src_mask(i,j)=ngrd_rm(mo)%mask_rho_o(i,j)
          end do
        end do
      end do
!  Mask out child portion of the src mask.
      do mo=1,Ngrids_roms-1
          do j=ngrd_rm(mo)%jstr_o,ngrd_rm(mo)%jend_o
            do i=ngrd_rm(mo)%istr_o,ngrd_rm(mo)%iend_o
              ngrd_rm(mo)%src_mask(i,j)=0
            end do
          end do
      end do
!  If a child grid, then set perimeter src_mask=0 because
!  this region comes from its parent.
      if (Ngrids_roms.gt.1) then
        do mo=2,Ngrids_roms
          i=ngrd_rm(mo)%xi_size
          do j=1,ngrd_rm(mo)%eta_size
            ngrd_rm(mo)%src_mask(1,j)=0
            ngrd_rm(mo)%src_mask(i,j)=0
          end do
          j=ngrd_rm(mo)%eta_size
          do i=1,ngrd_rm(mo)%xi_size
            ngrd_rm(mo)%src_mask(i,1)=0
            ngrd_rm(mo)%src_mask(i,j)=0
          end do
        end do
      end if
!  NOTICE piped swan and roms loops to end
      allocate(do_adjust(Ngrids_swan))
      do mw=1,Ngrids_swan
        do_adjust(mw)=0
        do mo=1,Ngrids_roms
!
!  Allocate the destination mask
!
          nx=ngrd_sw(mw)%Numx_swan
          ny=ngrd_sw(mw)%Numy_swan
          allocate(ngrd_sw(mw)%dst_mask(nx,ny))
          allocate(dst_mask_unlim(nx,ny))
!
!  Compute swan grid dst mask. Start by using locations on dst grid
!  where src_mask=1.
!
          do nx=1,ngrd_sw(mw)%Numx_swan
            do ny=1,ngrd_sw(mw)%Numy_swan
              dist_max=10e6
              Ikeep=1
              Jkeep=1
              xx2=ngrd_sw(mw)%lon_rho_w(nx,ny)
              yy2=ngrd_sw(mw)%lat_rho_w(nx,ny)
              do j=1,ngrd_rm(mo)%eta_size
                do i=1,ngrd_rm(mo)%xi_size
                  dist1=sqrt((ngrd_rm(mo)%lon_rho_o(i,j)-xx2)**2+       &
     &                       (ngrd_rm(mo)%lat_rho_o(i,j)-yy2)**2)
!                 xx1=ngrd_rm(mo)%lon_rho_o(i,j)
!                 yy1=ngrd_rm(mo)%lat_rho_o(i,j)
!                 dlon = xx1-xx2
!                 latrad1=abs(yy1*deg2rad)
!                 latrad2=abs(yy2*deg2rad)
!                 dep=cos(0.5*(latrad2+latrad1))*dlon
!                 dlat=yy2-yy1
!                 dist1=1852.0*60.0*sqrt(dlat**2+dep**2)
                  if(dist1<=dist_max)then
                    dist_max=dist1
                    Ikeep=i
                    Jkeep=j
                  endif
                enddo
              enddo
              ngrd_sw(mw)%dst_mask(nx,ny)=                              &
     &                                 ngrd_rm(mo)%src_mask(Ikeep,Jkeep)
            enddo
          enddo
!  Apply the swan grid Land/Sea mask to dst_mask
          do j=1,ngrd_sw(mw)%Numy_swan
            do i=1,ngrd_sw(mw)%Numx_swan
              ngrd_sw(mw)%dst_mask(i,j)=ngrd_sw(mw)%dst_mask(i,j)*      &
     &                                  ngrd_sw(mw)%mask_rho_w(i,j)
            end do
          end do
!  Remove points on dst mask that are outside of the src grid.
!  Create a vector of bounding src grid psi points.
          N=(ngrd_rm(mo)%xi_size+ngrd_rm(mo)%eta_size+2)*2-4
          allocate(XX(N), YY(N))
          count=0
          do i=1,ngrd_rm(mo)%eta_size+1
            count=count+1
            XX(count)=ngrd_rm(mo)%x_full_grid(1,i)
            YY(count)=ngrd_rm(mo)%y_full_grid(1,i)
          end do
          do i=2,ngrd_rm(mo)%xi_size+1
            count=count+1
            XX(count)=ngrd_rm(mo)%x_full_grid(i,ngrd_rm(mo)%eta_size+1)
            YY(count)=ngrd_rm(mo)%y_full_grid(i,ngrd_rm(mo)%eta_size+1)
          end do
          do i=ngrd_rm(mo)%eta_size,1,-1
            count=count+1
            XX(count)=ngrd_rm(mo)%x_full_grid(ngrd_rm(mo)%xi_size+1,i)
            YY(count)=ngrd_rm(mo)%y_full_grid(ngrd_rm(mo)%xi_size+1,i)
          end do
          do i=ngrd_rm(mo)%xi_size,2,-1
            count=count+1
            XX(count)=ngrd_rm(mo)%x_full_grid(i,1)
            YY(count)=ngrd_rm(mo)%y_full_grid(i,1)
          end do
!
!  For each dst pt call the pnpoly to see if it is inside the src grd.
!
          do nx=1,ngrd_sw(mw)%Numx_swan
            do ny=1,ngrd_sw(mw)%Numy_swan
              PX=ngrd_sw(mw)%lon_rho_w(nx,ny)
              PY=ngrd_sw(mw)%lat_rho_w(nx,ny)
              INOUT=0
              call PNPOLY( PX, PY, XX, YY, N, INOUT )
              ngrd_sw(mw)%dst_mask(nx,ny)=                              &
     &          ngrd_sw(mw)%dst_mask(nx,ny)*                            &
     &          MIN(INOUT+1,1)
            enddo
          enddo
          deallocate(XX, YY)
!
!  Compare src and dst lon points and see if points are coincident. 
!  The test looks for less than 2.0e-7 rads =  1.1459e-05 degs =~ 1.3m
!  If so, then offset the dst by 2.0e-7 rads
!
            tol=1.1459e-05
            do i=1,ngrd_rm(mo)%xi_size
              if (do_adjust(mw).eq.1) exit
              do j=1,ngrd_rm(mo)%eta_size
                if (do_adjust(mw).eq.1) exit
                lons=ngrd_rm(mo)%lon_rho_o(i,j)
                lats=ngrd_rm(mo)%lat_rho_o(i,j)
                do nx=1,ngrd_sw(mw)%Numx_swan
                  if (do_adjust(mw).eq.1) exit
                  do ny=1,ngrd_sw(mw)%Numy_swan
                    lond=ngrd_sw(mw)%lon_rho_w(nx,ny)
                    latd=ngrd_sw(mw)%lat_rho_w(nx,ny)
                    cff=min(abs(lons-lond),abs(lats-latd))
                    if (cff.le.tol) then
                      do_adjust(mw)=1
                  write(*,*) 'Do adjust for coincident points grid ', mw
                      exit
                    end if
                  enddo
                enddo
              enddo
            enddo
            if (do_adjust(mw).eq.1) then
              do nx=1,ngrd_sw(mw)%Numx_swan
                do ny=1,ngrd_sw(mw)%Numy_swan
                  ngrd_sw(mw)%lon_rho_w(nx,ny)=                         &
     &                                  ngrd_sw(mw)%lon_rho_w(nx,ny)+tol
                  ngrd_sw(mw)%lat_rho_w(nx,ny)=                         &
     &                                  ngrd_sw(mw)%lat_rho_w(nx,ny)+tol
                enddo
              enddo
              do nx=1,ngrd_sw(mw)%Numx_swan+1
                do ny=1,ngrd_sw(mw)%Numy_swan+1
                  ngrd_sw(mw)%x_full_grid(nx,ny)=                       &
     &                                ngrd_sw(mw)%x_full_grid(nx,ny)+tol
                  ngrd_sw(mw)%y_full_grid(nx,ny)=                       &
     &                                ngrd_sw(mw)%y_full_grid(nx,ny)+tol
                enddo
              enddo
            end if
!
!  Send the data to SCRIP routines
!
          write(mw_string,'(i1)')mw
          write(mo_string,'(i1)')mo
          interp_file1='ocn'//trim(mo_string)//'_to_'//'wav'//          &
     &                       trim(mw_string)//'_weights.nc'
          interp_file2='unknown'
          map1_name='ROMS to SWAN distwgt Mapping'
          map2_name='unknown'

          write(stdout,*)"============================================="
          write(stdout,*)"Begin mapping between the two grids"
          write(stdout,*)"---------------------------------------------"
          write(stdout,*)"The interp file is: ", interp_file1
          write(stdout,*)"The src grid is: ", roms_grids(mo)
          write(stdout,*)"The dst grid is: ", swan_coord(mw)

          grid1_file=roms_grids(mo)
          grid2_file=swan_coord(mw)

          counter_grid=counter_grid+1

          call scrip_package(grid1_file, grid2_file,                    &
     &                       ngrd_rm(mo)%xi_size,                       &
     &                       ngrd_rm(mo)%eta_size,                      &
     &                       ngrd_rm(mo)%lon_rho_o,                     &
     &                       ngrd_rm(mo)%lat_rho_o,                     &
     &                       ngrd_rm(mo)%x_full_grid,                   &
     &                       ngrd_rm(mo)%y_full_grid,                   &
     &                       ngrd_sw(mw)%Numx_swan,                     &
     &                       ngrd_sw(mw)%Numy_swan,                     &
     &                       ngrd_sw(mw)%lon_rho_w,                     &
     &                       ngrd_sw(mw)%lat_rho_w,                     &
     &                       ngrd_sw(mw)%x_full_grid,                   &
     &                       ngrd_sw(mw)%y_full_grid,                   &
     &                       interp_file1, interp_file2,                &
     &                       map1_name, map2_name,                      &
     &                       ngrd_rm(mo)%mask_rho_o,                    &
     &                       ngrd_sw(mw)%dst_mask,                      &
     &                       dst_mask_unlim,                            &
     &                       counter_grid,                              &
     &                       Ngrids_comb_total,                         &
     &                       output_ncfile)
!
!  Undo the grid adjust.
          if ((do_adjust(mw).eq.1).and.(mo.lt.Ngrids_roms)) then
            do nx=1,ngrd_sw(mw)%Numx_swan
              do ny=1,ngrd_sw(mw)%Numy_swan
                ngrd_sw(mw)%lon_rho_w(nx,ny)=                           &
     &                                ngrd_sw(mw)%lon_rho_w(nx,ny)-tol
                ngrd_sw(mw)%lat_rho_w(nx,ny)=                           &
     &                                ngrd_sw(mw)%lat_rho_w(nx,ny)-tol
              enddo
            enddo
            do nx=1,ngrd_sw(mw)%Numx_swan+1
              do ny=1,ngrd_sw(mw)%Numy_swan+1
                ngrd_sw(mw)%x_full_grid(nx,ny)=                         &
     &                              ngrd_sw(mw)%x_full_grid(nx,ny)-tol
                ngrd_sw(mw)%y_full_grid(nx,ny)=                         &
     &                              ngrd_sw(mw)%y_full_grid(nx,ny)-tol
              enddo
            enddo
          end if
          do_adjust(mw)=0
!
          deallocate(ngrd_sw(mw)%dst_mask)
          deallocate(dst_mask_unlim)
        end do
      end do
      deallocate(do_adjust)
      do mo=1,Ngrids_roms
        deallocate(ngrd_rm(mo)%src_mask)
      end do

      end subroutine ocn2wav_mask

!======================================================================

      subroutine wav2ocn_mask()

      implicit none

!     local variables
      real(dbl_kind), allocatable :: XX(:), YY(:)

!  Allocate the source mask
      do mw=1,Ngrids_swan
        nx=ngrd_sw(mw)%Numx_swan
        ny=ngrd_sw(mw)%Numy_swan
        allocate(ngrd_sw(mw)%src_mask(nx,ny))
      end do
!  Save the wav grid mask to src_mask
      do mw=1,Ngrids_swan
        do j=1,ngrd_sw(mw)%Numy_swan
          do i=1,ngrd_sw(mw)%Numx_swan
            ngrd_sw(mw)%src_mask(i,j)=ngrd_sw(mw)%mask_rho_w(i,j)
          end do
        end do
      end do
!  Mask out child portion of the src mask.
      do mw=1,Ngrids_swan-1
        do j=ngrd_sw(mw)%jstr_w,ngrd_sw(mw)%jend_w
          do i=ngrd_sw(mw)%istr_w,ngrd_sw(mw)%iend_w
            ngrd_sw(mw)%src_mask(i,j)=0
           end do
        end do
      end do
!  If a child grid, then set perimeter src_mask=0 because
!  this region comes from its parent.
      if (Ngrids_swan.gt.1) then
        do mw=2,Ngrids_swan
          i=ngrd_sw(mw)%Numx_swan
          do j=1,ngrd_sw(mw)%Numy_swan
            ngrd_sw(mw)%src_mask(1,j)=0
            ngrd_sw(mw)%src_mask(i,j)=0
          end do
          j=ngrd_sw(mw)%Numy_swan
          do i=1,ngrd_sw(mw)%Numx_swan
            ngrd_sw(mw)%src_mask(i,1)=0
            ngrd_sw(mw)%src_mask(i,j)=0
          end do
        end do
      end if
!  NOTICE piped swan and roms loops to end
      do mo=1,Ngrids_roms
        do mw=1,Ngrids_swan
!
!  Allocate the destination mask
!
          nx=ngrd_rm(mo)%xi_size
          ny=ngrd_rm(mo)%eta_size
          allocate(ngrd_rm(mo)%dst_mask(nx,ny))
          allocate(dst_mask_unlim(nx,ny))
!
!  Compute roms grid dst mask. Start by using locations on dst grid
!  where src_mask=1.
!
          do nx=1,ngrd_rm(mo)%xi_size
            do ny=1,ngrd_rm(mo)%eta_size
              dist_max=10e6
              Ikeep=1
              Jkeep=1
              xx2=ngrd_rm(mo)%lon_rho_o(nx,ny)
              yy2=ngrd_rm(mo)%lat_rho_o(nx,ny)
              do j=1,ngrd_sw(mw)%Numy_swan
                do i=1,ngrd_sw(mw)%Numx_swan
                  dist1=sqrt((ngrd_sw(mw)%lon_rho_w(i,j)-xx2)**2+       &
     &                       (ngrd_sw(mw)%lat_rho_w(i,j)-yy2)**2)
!                 xx1=ngrd_sw(mw)%lon_rho_w(i,j)
!                 yy1=ngrd_sw(mw)%lat_rho_w(i,j)
!                 dlon = xx1-xx2
!                 latrad1=abs(yy1*deg2rad)
!                 latrad2=abs(yy2*deg2rad)
!                 dep=cos(0.5*(latrad2+latrad1))*dlon
!                 dlat=yy2-yy1
!                 dist1=1852.0*60.0*sqrt(dlat**2+dep**2)
                  if(dist1<=dist_max)then
                    dist_max=dist1
                    Ikeep=i
                    Jkeep=j
                  endif
                enddo
              enddo
              ngrd_rm(mo)%dst_mask(nx,ny)=                              &
     &                                 ngrd_sw(mw)%src_mask(Ikeep,Jkeep)
            enddo
          enddo
!  Apply the roms grid Land/Sea mask to dst_mask
          do j=1,ngrd_rm(mo)%eta_size
            do i=1,ngrd_rm(mo)%xi_size
              ngrd_rm(mo)%dst_mask(i,j)=ngrd_rm(mo)%dst_mask(i,j)*      &
     &                                  ngrd_rm(mo)%mask_rho_o(i,j)
            end do
          end do
!  Remove points on dst mask that are outside of the src grid.
!  Create a vector of bounding src grid psi points.
          N=(ngrd_sw(mw)%Numx_swan+ngrd_sw(mw)%Numy_swan+2)*2-4
          allocate(XX(N), YY(N))
          count=0
          do i=1,ngrd_sw(mw)%Numy_swan+1
            count=count+1
            XX(count)=ngrd_sw(mw)%x_full_grid(1,i)
            YY(count)=ngrd_sw(mw)%y_full_grid(1,i)
          end do
          do i=2,ngrd_sw(mw)%Numx_swan+1
            count=count+1
            XX(count)=ngrd_sw(mw)%x_full_grid(i,ngrd_sw(mw)%Numy_swan+1)
            YY(count)=ngrd_sw(mw)%y_full_grid(i,ngrd_sw(mw)%Numy_swan+1)
          end do
          do i=ngrd_sw(mw)%Numy_swan,1,-1
            count=count+1
            XX(count)=ngrd_sw(mw)%x_full_grid(ngrd_sw(mw)%Numx_swan+1,i)
            YY(count)=ngrd_sw(mw)%y_full_grid(ngrd_sw(mw)%Numx_swan+1,i)
          end do
          do i=ngrd_sw(mw)%Numx_swan,2,-1
            count=count+1
            XX(count)=ngrd_sw(mw)%x_full_grid(i,1)
            YY(count)=ngrd_sw(mw)%y_full_grid(i,1)
          end do
!
!  For each dst pt call the pnpoly to see if it is inside the src grd.
!
          do nx=1,ngrd_rm(mo)%xi_size
            do ny=1,ngrd_rm(mo)%eta_size
              PX=ngrd_rm(mo)%lon_rho_o(nx,ny)
              PY=ngrd_rm(mo)%lat_rho_o(nx,ny)
              INOUT=0
              call PNPOLY( PX, PY, XX, YY, N, INOUT )
              ngrd_rm(mo)%dst_mask(nx,ny)=                              &
     &          ngrd_rm(mo)%dst_mask(nx,ny)*                            &
     &          MIN(INOUT+1,1)
            enddo
          enddo
          deallocate(XX, YY)
!
!  Send the data to SCRIP routines
!
          write(mo_string,'(i1)')mo
          write(mw_string,'(i1)')mw
          interp_file1='wav'//trim(mw_string)//'_to_'//'ocn'//          &
     &                       trim(mo_string)//'_weights.nc'
          interp_file2='unknown'
          map1_name='SWAN to ROMS distwgt Mapping'
          map2_name='unknown'

          write(stdout,*)"============================================="
          write(stdout,*)"Begin mapping between the two grids"
          write(stdout,*)"---------------------------------------------"
          write(stdout,*)"The interp file is: ", interp_file1
          write(stdout,*)"The src grid is: ", swan_coord(mw)
          write(stdout,*)"The dst grid is: ", roms_grids(mo)

          grid1_file=swan_coord(mw)
          grid2_file=roms_grids(mo)

          counter_grid=counter_grid+1

          call scrip_package(grid1_file, grid2_file,                    &
     &                       ngrd_sw(mw)%Numx_swan,                     &
     &                       ngrd_sw(mw)%Numy_swan,                     &
     &                       ngrd_sw(mw)%lon_rho_w,                     &
     &                       ngrd_sw(mw)%lat_rho_w,                     &
     &                       ngrd_sw(mw)%x_full_grid,                   &
     &                       ngrd_sw(mw)%y_full_grid,                   &
     &                       ngrd_rm(mo)%xi_size,                       &
     &                       ngrd_rm(mo)%eta_size,                      &
     &                       ngrd_rm(mo)%lon_rho_o,                     &
     &                       ngrd_rm(mo)%lat_rho_o,                     &
     &                       ngrd_rm(mo)%x_full_grid,                   &
     &                       ngrd_rm(mo)%y_full_grid,                   &
     &                       interp_file1, interp_file2,                &
     &                       map1_name, map2_name,                      &
     &                       ngrd_sw(mw)%mask_rho_w,                    &
     &                       ngrd_rm(mo)%dst_mask,                      &
     &                       dst_mask_unlim,                            &
     &                       counter_grid,                              &
     &                       Ngrids_comb_total,                         &
     &                       output_ncfile)

          deallocate(ngrd_rm(mo)%dst_mask)
          deallocate(dst_mask_unlim)
        end do
      end do
      do mw=1,Ngrids_swan
        deallocate(ngrd_sw(mw)%src_mask)
      end do

      end subroutine wav2ocn_mask

!======================================================================

      subroutine ocn2atm_mask()

      implicit none

!     local variables
      real(dbl_kind), allocatable :: XX(:), YY(:)

!  Allocate the source mask
      do mo=1,Ngrids_roms
        nx=ngrd_rm(mo)%xi_size
        ny=ngrd_rm(mo)%eta_size
        allocate(ngrd_rm(mo)%src_mask(nx,ny))
      end do
!  Save the ocean grid mask to src_mask
      do mo=1,Ngrids_roms
        do j=1,ngrd_rm(mo)%eta_size
          do i=1,ngrd_rm(mo)%xi_size
            ngrd_rm(mo)%src_mask(i,j)=ngrd_rm(mo)%mask_rho_o(i,j)
          end do
        end do
      end do
!  Mask out child portion of the src mask.
      do mo=1,Ngrids_roms-1
          do j=ngrd_rm(mo)%jstr_o,ngrd_rm(mo)%jend_o
            do i=ngrd_rm(mo)%istr_o,ngrd_rm(mo)%iend_o
              ngrd_rm(mo)%src_mask(i,j)=0
            end do
          end do
      end do
!  If a child grid, then set perimeter src_mask=0 because
!  this region comes from its parent.
      if (Ngrids_roms.gt.1) then
        do mo=2,Ngrids_roms
          i=ngrd_rm(mo)%xi_size
          do j=1,ngrd_rm(mo)%eta_size
            ngrd_rm(mo)%src_mask(1,j)=0
            ngrd_rm(mo)%src_mask(i,j)=0
          end do
          j=ngrd_rm(mo)%eta_size
          do i=1,ngrd_rm(mo)%xi_size
            ngrd_rm(mo)%src_mask(i,1)=0
            ngrd_rm(mo)%src_mask(i,j)=0
          end do
        end do
      end if
!  NOTICE piped wrf and roms loops to end
      allocate(do_adjust(Ngrids_wrf))
      do ma=1,Ngrids_wrf
        do_adjust(ma)=0
        do mo=1,Ngrids_roms
!
!  Allocate the destination mask
!
          nx=ngrd_wr(ma)%we_size
          ny=ngrd_wr(ma)%sn_size
          allocate(ngrd_wr(ma)%dst_mask(nx,ny))
          allocate(dst_mask_unlim(nx,ny))
!
!  Compute wrf grid dst mask. Start by using locations on dst grid
!  where src_mask=1.
!
          offset=MIN(mo-1,1)
          do nx=1,ngrd_wr(ma)%we_size
            do ny=1,ngrd_wr(ma)%sn_size
              dist_max=10e6
              Ikeep=1
              Jkeep=1
              xx2=ngrd_wr(ma)%lon_rho_a(nx,ny)
              yy2=ngrd_wr(ma)%lat_rho_a(nx,ny)
              do j=1+offset,ngrd_rm(mo)%eta_size-offset
                do i=1+offset,ngrd_rm(mo)%xi_size-offset
                  xx1=ngrd_rm(mo)%lon_rho_o(i,j)
                  yy1=ngrd_rm(mo)%lat_rho_o(i,j)
                  dlon = xx1-xx2
                  latrad1=abs(yy1*deg2rad)
                  latrad2=abs(yy2*deg2rad)
                  dep=cos(0.5*(latrad2+latrad1))*dlon
                  dlat=yy2-yy1
                  dist1=1852.0*60.0*sqrt(dlat**2+dep**2)
                  if(dist1<=dist_max)then
                    dist_max=dist1
                    Ikeep=i
                    Jkeep=j
                  endif
                enddo
              enddo
              ngrd_wr(ma)%dst_mask(nx,ny)=                              &
     &                                 ngrd_rm(mo)%src_mask(Ikeep,Jkeep)
            enddo
          enddo
!  Apply the wrf grid Land/Sea mask to dst_mask
          do j=1,ngrd_wr(ma)%sn_size
            do i=1,ngrd_wr(ma)%we_size
              ngrd_wr(ma)%dst_mask(i,j)=ngrd_wr(ma)%dst_mask(i,j)*      &
     &                                  ngrd_wr(ma)%mask_rho_a(i,j)
            end do
          end do
!  Remove points on dst mask that are outside of the src grid.
!  Create a vector of bounding src grid psi points.
          N=(ngrd_rm(mo)%xi_size+ngrd_rm(mo)%eta_size+2)*2-4
          allocate(XX(N), YY(N))
          count=0
          do i=1,ngrd_rm(mo)%eta_size+1
            count=count+1
            XX(count)=ngrd_rm(mo)%x_full_grid(1,i)
            YY(count)=ngrd_rm(mo)%y_full_grid(1,i)
          end do
          do i=2,ngrd_rm(mo)%xi_size+1
            count=count+1
            XX(count)=ngrd_rm(mo)%x_full_grid(i,ngrd_rm(mo)%eta_size+1)
            YY(count)=ngrd_rm(mo)%y_full_grid(i,ngrd_rm(mo)%eta_size+1)
          end do
          do i=ngrd_rm(mo)%eta_size,1,-1
            count=count+1
            XX(count)=ngrd_rm(mo)%x_full_grid(ngrd_rm(mo)%xi_size+1,i)
            YY(count)=ngrd_rm(mo)%y_full_grid(ngrd_rm(mo)%xi_size+1,i)
          end do
          do i=ngrd_rm(mo)%xi_size,2,-1
            count=count+1
            XX(count)=ngrd_rm(mo)%x_full_grid(i,1)
            YY(count)=ngrd_rm(mo)%y_full_grid(i,1)
          end do
!
!  For each dst pt call the pnpoly to see if it is inside the src grd.
!
          do nx=1,ngrd_wr(ma)%we_size
            do ny=1,ngrd_wr(ma)%sn_size
              PX=ngrd_wr(ma)%lon_rho_a(nx,ny)
              PY=ngrd_wr(ma)%lat_rho_a(nx,ny)
              INOUT=0
              call PNPOLY( PX, PY, XX, YY, N, INOUT )
              ngrd_wr(ma)%dst_mask(nx,ny)=                              &
     &          ngrd_wr(ma)%dst_mask(nx,ny)*                            &
     &          MIN(INOUT+1,1)
            enddo
          enddo
          deallocate(XX, YY)
!
!  Compare src and dst lon points and see if points are coincident. 
!  The test looks for less than 2.0e-7 rads =  1.1459e-05 degs =~ 1.3m
!  If so, then offset the dst by 2.0e-7 rads
!
            tol=1.1459e-05
            do i=1,ngrd_rm(mo)%xi_size
              if (do_adjust(ma).eq.1) exit
              do j=1,ngrd_rm(mo)%eta_size
                if (do_adjust(ma).eq.1) exit
                lons=ngrd_rm(mo)%lon_rho_o(i,j)
                lats=ngrd_rm(mo)%lat_rho_o(i,j)
                do nx=1,ngrd_wr(ma)%we_size
                  if (do_adjust(ma).eq.1) exit
                  do ny=1,ngrd_wr(ma)%sn_size
                    lond=ngrd_wr(ma)%lon_rho_a(nx,ny)
                    latd=ngrd_wr(ma)%lat_rho_a(nx,ny)
                    cff=min(abs(lons-lond),abs(lats-latd))
                    if (cff.le.tol) then
                      do_adjust(ma)=1
                  write(*,*) 'Do adjust for coincident points grid ', ma
                      exit
                    end if
                  enddo
                enddo
              enddo
            enddo
            if (do_adjust(ma).eq.1) then
              do nx=1,ngrd_wr(ma)%we_size
                do ny=1,ngrd_wr(ma)%sn_size
                  ngrd_wr(ma)%lon_rho_a(nx,ny)=                         &
     &                                  ngrd_wr(ma)%lon_rho_a(nx,ny)+tol
                  ngrd_wr(ma)%lat_rho_a(nx,ny)=                         &
     &                                  ngrd_wr(ma)%lat_rho_a(nx,ny)+tol
                enddo
              enddo
              do nx=1,ngrd_wr(ma)%we_size+1
                do ny=1,ngrd_wr(ma)%sn_size+1
                  ngrd_wr(ma)%x_full_grid(nx,ny)=                       &
     &                                ngrd_wr(ma)%x_full_grid(nx,ny)+tol
                  ngrd_wr(ma)%y_full_grid(nx,ny)=                       &
     &                                ngrd_wr(ma)%y_full_grid(nx,ny)+tol
                enddo
              enddo
            end if
!
!  Send the data to SCRIP routines
!
          write(ma_string,'(i1)')ma
          write(mo_string,'(i1)')mo
          interp_file1='ocn'//trim(mo_string)//'_to_'//'atm'//          &
     &                       trim(ma_string)//'_weights.nc'
          interp_file2='unknown'
          map1_name='ROMS to WRF distwgt Mapping'
          map2_name='unknown'

          write(stdout,*)"============================================="
          write(stdout,*)"Begin mapping between the two grids"
          write(stdout,*)"---------------------------------------------"
          write(stdout,*)"The interp file is: ", interp_file1
          write(stdout,*)"The src grid is: ", roms_grids(mo)
          write(stdout,*)"The dst grid is: ", wrf_grids(ma)

          grid1_file=roms_grids(mo)
          grid2_file=wrf_grids(ma)

          counter_grid=counter_grid+1

          call scrip_package(grid1_file, grid2_file,                    &
     &                       ngrd_rm(mo)%xi_size,                       &
     &                       ngrd_rm(mo)%eta_size,                      &
     &                       ngrd_rm(mo)%lon_rho_o,                     &
     &                       ngrd_rm(mo)%lat_rho_o,                     &
     &                       ngrd_rm(mo)%x_full_grid,                   &
     &                       ngrd_rm(mo)%y_full_grid,                   &
     &                       ngrd_wr(ma)%we_size,                       &
     &                       ngrd_wr(ma)%sn_size,                       &
     &                       ngrd_wr(ma)%lon_rho_a,                     &
     &                       ngrd_wr(ma)%lat_rho_a,                     &
     &                       ngrd_wr(ma)%x_full_grid,                   &
     &                       ngrd_wr(ma)%y_full_grid,                   &
     &                       interp_file1, interp_file2,                &
     &                       map1_name, map2_name,                      &
     &                       ngrd_rm(mo)%mask_rho_o,                    &
     &                       ngrd_wr(ma)%dst_mask,                      &
     &                       dst_mask_unlim,                            &
     &                       counter_grid,                              &
     &                       Ngrids_comb_total,                         &
     &                       output_ncfile)
!
!  Undo the grid adjust.
          if ((do_adjust(ma).eq.1).and.(mo.lt.Ngrids_roms)) then
            do nx=1,ngrd_wr(ma)%we_size
              do ny=1,ngrd_wr(ma)%sn_size
                ngrd_wr(ma)%lon_rho_a(nx,ny)=                           &
     &                                ngrd_wr(ma)%lon_rho_a(nx,ny)-tol
                ngrd_wr(ma)%lat_rho_a(nx,ny)=                           &
     &                                ngrd_wr(ma)%lat_rho_a(nx,ny)-tol
              enddo
            enddo
            do nx=1,ngrd_wr(ma)%we_size+1
              do ny=1,ngrd_wr(ma)%sn_size+1
                ngrd_wr(ma)%x_full_grid(nx,ny)=                         &
     &                              ngrd_wr(ma)%x_full_grid(nx,ny)-tol
                ngrd_wr(ma)%y_full_grid(nx,ny)=                         &
     &                              ngrd_wr(ma)%y_full_grid(nx,ny)-tol
              enddo
            enddo
          end if
          do_adjust(ma)=0
!
          deallocate(ngrd_wr(ma)%dst_mask)
          deallocate(dst_mask_unlim)
        end do
      end do
      deallocate(do_adjust)
      do mo=1,Ngrids_roms
        deallocate(ngrd_rm(mo)%src_mask)
      end do

      end subroutine ocn2atm_mask

!======================================================================

      subroutine atm2ocn_mask()

      implicit none

!     local variables
      real(dbl_kind), allocatable :: XX(:), YY(:)

!  Allocate the source mask
      do ma=1,Ngrids_wrf
        nx=ngrd_wr(ma)%we_size
        ny=ngrd_wr(ma)%sn_size
        allocate(ngrd_wr(ma)%src_mask(nx,ny))
      end do
!  Save the atm grid mask to src_mask
      do ma=1,Ngrids_wrf
        do j=1,ngrd_wr(ma)%sn_size
          do i=1,ngrd_wr(ma)%we_size
            ngrd_wr(ma)%src_mask(i,j)=ngrd_wr(ma)%mask_rho_a(i,j)
          end do
        end do
      end do
!  Mask out child portion of the src mask.
      do ma=1,Ngrids_wrf-1
        if (wrf_grids(ma+1)/="moving") then
          do j=ngrd_wr(ma)%jstr_a,ngrd_wr(ma)%jend_a
            do i=ngrd_wr(ma)%istr_a,ngrd_wr(ma)%iend_a
              ngrd_wr(ma)%src_mask(i,j)=0
            end do
          end do
        end if
      end do
!  NOTICE piped wrf and roms loops to end
      do mo=1,Ngrids_roms
        do ma=1,Ngrids_wrf
!
!  create an  unlimited src mask to determine dst cells removed
!  due to src land masking.
!
          nx=ngrd_wr(ma)%we_size
          ny=ngrd_wr(ma)%sn_size
          allocate(src_mask_unlim(nx,ny))
          do j=1,ngrd_wr(ma)%sn_size
            do i=1,ngrd_wr(ma)%we_size
              src_mask_unlim(i,j)=1
            end do
          end do
!  Mask out child portion of the src mask.
          if ((ma.lt.Ngrids_wrf).and.(wrf_grids(ma+1)/="moving")) then
            do j=ngrd_wr(ma)%jstr_a,ngrd_wr(ma)%jend_a
              do i=ngrd_wr(ma)%istr_a,ngrd_wr(ma)%iend_a
                src_mask_unlim(i,j)=0
              end do
            end do
          end if
!
!  Allocate the destination mask
!
          nx=ngrd_rm(mo)%xi_size
          ny=ngrd_rm(mo)%eta_size
          allocate(ngrd_rm(mo)%dst_mask(nx,ny))
          allocate(dst_mask_unlim(nx,ny))
!
!  Compute roms grid dst mask. Start by using locations on dst grid
!  where src_mask=1.
!
          do nx=1,ngrd_rm(mo)%xi_size
            do ny=1,ngrd_rm(mo)%eta_size
              dist_max=10e6
              Ikeep=1
              Jkeep=1
              xx2=ngrd_rm(mo)%lon_rho_o(nx,ny)
              yy2=ngrd_rm(mo)%lat_rho_o(nx,ny)
              do j=1,ngrd_wr(ma)%sn_size
                do i=1,ngrd_wr(ma)%we_size
                  xx1=ngrd_wr(ma)%lon_rho_a(i,j)
                  yy1=ngrd_wr(ma)%lat_rho_a(i,j)
                  dlon = xx1-xx2
                  latrad1=abs(yy1*deg2rad)
                  latrad2=abs(yy2*deg2rad)
                  dep=cos(0.5*(latrad2+latrad1))*dlon
                  dlat=yy2-yy1
                  dist1=1852.0*60.0*sqrt(dlat**2+dep**2)
                  if(dist1<=dist_max)then
                    dist_max=dist1
                    Ikeep=i
                    Jkeep=j
                  endif
                enddo
              enddo
              count=ngrd_wr(ma)%src_mask(Ikeep,Jkeep)
              ngrd_rm(mo)%dst_mask(nx,ny)=count
              count=src_mask_unlim(Ikeep,Jkeep)
              dst_mask_unlim(nx,ny)=count
            enddo
          enddo
!  Apply the roms grid Land/Sea mask to dst_mask
!  Check to see if there are any dst points that have valid mask
!  but are not being picked up by the src grid
          do j=1,ngrd_rm(mo)%eta_size
            do i=1,ngrd_rm(mo)%xi_size
              ngrd_rm(mo)%dst_mask(i,j)=ngrd_rm(mo)%dst_mask(i,j)*      &
     &                                  ngrd_rm(mo)%mask_rho_o(i,j)
              dst_mask_unlim(i,j)=dst_mask_unlim(i,j)*                  &
     &                                  ngrd_rm(mo)%mask_rho_o(i,j)
            end do
          end do
!  Remove points on dst mask that are outside of the src grid.
!  Create a vector of bounding src grid psi points.
          N=(ngrd_wr(ma)%we_size+ngrd_wr(ma)%sn_size+2)*2-4
          allocate(XX(N), YY(N))
          count=0
          do i=1,ngrd_wr(ma)%sn_size+1
            count=count+1
            XX(count)=ngrd_wr(ma)%x_full_grid(1,i)
            YY(count)=ngrd_wr(ma)%y_full_grid(1,i)
          end do
          do i=2,ngrd_wr(ma)%we_size+1
            count=count+1
            XX(count)=ngrd_wr(ma)%x_full_grid(i,ngrd_wr(ma)%sn_size+1)
            YY(count)=ngrd_wr(ma)%y_full_grid(i,ngrd_wr(ma)%sn_size+1)
          end do
          do i=ngrd_wr(ma)%sn_size,1,-1
            count=count+1
            XX(count)=ngrd_wr(ma)%x_full_grid(ngrd_wr(ma)%we_size+1,i)
            YY(count)=ngrd_wr(ma)%y_full_grid(ngrd_wr(ma)%we_size+1,i)
          end do
          do i=ngrd_wr(ma)%we_size,2,-1
            count=count+1
            XX(count)=ngrd_wr(ma)%x_full_grid(i,1)
            YY(count)=ngrd_wr(ma)%y_full_grid(i,1)
          end do
!
!  For each dst pt call the pnpoly to see if it is inside the src grd.
!
          do nx=1,ngrd_rm(mo)%xi_size
            do ny=1,ngrd_rm(mo)%eta_size
              PX=ngrd_rm(mo)%lon_rho_o(nx,ny)
              PY=ngrd_rm(mo)%lat_rho_o(nx,ny)
              INOUT=0
              call PNPOLY( PX, PY, XX, YY, N, INOUT )
              ngrd_rm(mo)%dst_mask(nx,ny)=                              &
     &          ngrd_rm(mo)%dst_mask(nx,ny)*MIN(INOUT+1,1)
              dst_mask_unlim(nx,ny)=                                    &
     &          dst_mask_unlim(nx,ny)*MIN(INOUT+1,1)
            enddo
          enddo
          deallocate(XX, YY)
!
!  Send the data to SCRIP routines
!
          write(mo_string,'(i1)')mo
          write(ma_string,'(i1)')ma
          interp_file1='atm'//trim(ma_string)//'_to_'//'ocn'//          &
     &                       trim(mo_string)//'_weights.nc'
          interp_file2='unknown'
          map1_name='WRF to ROMS distwgt Mapping'
          map2_name='unknown'

          write(stdout,*)"============================================="
          write(stdout,*)"Begin mapping between the two grids"
          write(stdout,*)"---------------------------------------------"
          write(stdout,*)"The interp file is: ", interp_file1
          write(stdout,*)"The src grid is: ", wrf_grids(ma)
          write(stdout,*)"The dst grid is: ", roms_grids(mo)

          grid1_file=wrf_grids(ma)
          grid2_file=roms_grids(mo)

          counter_grid=counter_grid+1

          call scrip_package(grid1_file, grid2_file,                    &
     &                       ngrd_wr(ma)%we_size,                       &
     &                       ngrd_wr(ma)%sn_size,                       &
     &                       ngrd_wr(ma)%lon_rho_a,                     &
     &                       ngrd_wr(ma)%lat_rho_a,                     &
     &                       ngrd_wr(ma)%x_full_grid,                   &
     &                       ngrd_wr(ma)%y_full_grid,                   &
     &                       ngrd_rm(mo)%xi_size,                       &
     &                       ngrd_rm(mo)%eta_size,                      &
     &                       ngrd_rm(mo)%lon_rho_o,                     &
     &                       ngrd_rm(mo)%lat_rho_o,                     &
     &                       ngrd_rm(mo)%x_full_grid,                   &
     &                       ngrd_rm(mo)%y_full_grid,                   &
     &                       interp_file1, interp_file2,                &
     &                       map1_name, map2_name,                      &
     &                       ngrd_wr(ma)%mask_rho_a,                    &
     &                       ngrd_rm(mo)%dst_mask,                      &
     &                       dst_mask_unlim,                            &
     &                       counter_grid,                              &
     &                       Ngrids_comb_total,                         &
     &                       output_ncfile)
!
          deallocate(ngrd_rm(mo)%dst_mask)
          deallocate(src_mask_unlim, dst_mask_unlim)
        end do
      end do
      do ma=1,Ngrids_wrf
        deallocate(ngrd_wr(ma)%src_mask)
      end do

      end subroutine atm2ocn_mask

!======================================================================

      subroutine atm2wav_mask()

      implicit none

!     local variables
      real(dbl_kind), allocatable :: XX(:), YY(:)

!  Allocate the source mask
      do ma=1,Ngrids_wrf
        nx=ngrd_wr(ma)%we_size
        ny=ngrd_wr(ma)%sn_size
        allocate(ngrd_wr(ma)%src_mask(nx,ny))
      end do
!  Save the atm grid mask to src_mask
      do ma=1,Ngrids_wrf
        do j=1,ngrd_wr(ma)%sn_size
          do i=1,ngrd_wr(ma)%we_size
            ngrd_wr(ma)%src_mask(i,j)=ngrd_wr(ma)%mask_rho_a(i,j)
          end do
        end do
      end do
!  Mask out child portion of the src mask.
      do ma=1,Ngrids_wrf-1
        if (wrf_grids(ma+1)/="moving") then
          do j=ngrd_wr(ma)%jstr_a,ngrd_wr(ma)%jend_a
            do i=ngrd_wr(ma)%istr_a,ngrd_wr(ma)%iend_a
              ngrd_wr(ma)%src_mask(i,j)=0
            end do
          end do
        end if
      end do
!
!  NOTICE piped wrf and swan loops to end
      do mw=1,Ngrids_swan
        do ma=1,Ngrids_wrf
!
!  Allocate the destination mask
!
          nx=ngrd_sw(mw)%Numx_swan
          ny=ngrd_sw(mw)%Numy_swan
          allocate(ngrd_sw(mw)%dst_mask(nx,ny))
          allocate(dst_mask_unlim(nx,ny))
!
!  Compute swan grid dst mask. Start by using locations on dst grid
!  where src_mask=1.
!
          do nx=1,ngrd_sw(mw)%Numx_swan
            do ny=1,ngrd_sw(mw)%Numy_swan
              dist_max=10e6
              Ikeep=1
              Jkeep=1
              xx2=ngrd_sw(mw)%lon_rho_w(nx,ny)
              yy2=ngrd_sw(mw)%lat_rho_w(nx,ny)
              do j=1,ngrd_wr(ma)%sn_size
                do i=1,ngrd_wr(ma)%we_size
                  xx1=ngrd_wr(ma)%lon_rho_a(i,j)
                  yy1=ngrd_wr(ma)%lat_rho_a(i,j)
                  dlon = xx1-xx2
                  latrad1=abs(yy1*deg2rad)
                  latrad2=abs(yy2*deg2rad)
                  dep=cos(0.5*(latrad2+latrad1))*dlon
                  dlat=yy2-yy1
                  dist1=1852.0*60.0*sqrt(dlat**2+dep**2)
                  if(dist1<=dist_max)then
                    dist_max=dist1
                    Ikeep=i
                    Jkeep=j
                  endif
                enddo
              enddo
              ngrd_sw(mw)%dst_mask(nx,ny)=                              &
     &                                 ngrd_wr(ma)%src_mask(Ikeep,Jkeep)
            enddo
          enddo
!  Apply the swan grid Land/Sea mask to dst_mask
          do j=1,ngrd_sw(mw)%Numy_swan
            do i=1,ngrd_sw(mw)%Numx_swan
              ngrd_sw(mw)%dst_mask(i,j)=ngrd_sw(mw)%dst_mask(i,j)*      &
     &                                  ngrd_sw(mw)%mask_rho_w(i,j)
            end do
          end do
!  Remove points on dst mask that are outside of the src grid.
!  Create a vector of bounding src grid psi points.
          N=(ngrd_wr(ma)%we_size+ngrd_wr(ma)%sn_size+2)*2-4
          allocate(XX(N), YY(N))
          count=0
          do i=1,ngrd_wr(ma)%sn_size+1
            count=count+1
            XX(count)=ngrd_wr(ma)%x_full_grid(1,i)
            YY(count)=ngrd_wr(ma)%y_full_grid(1,i)
          end do
          do i=2,ngrd_wr(ma)%we_size+1
            count=count+1
            XX(count)=ngrd_wr(ma)%x_full_grid(i,ngrd_wr(ma)%sn_size+1)
            YY(count)=ngrd_wr(ma)%y_full_grid(i,ngrd_wr(ma)%sn_size+1)
          end do
          do i=ngrd_wr(ma)%sn_size,1,-1
            count=count+1
            XX(count)=ngrd_wr(ma)%x_full_grid(ngrd_wr(ma)%we_size+1,i)
            YY(count)=ngrd_wr(ma)%y_full_grid(ngrd_wr(ma)%we_size+1,i)
          end do
          do i=ngrd_wr(ma)%we_size,2,-1
            count=count+1
            XX(count)=ngrd_wr(ma)%x_full_grid(i,1)
            YY(count)=ngrd_wr(ma)%y_full_grid(i,1)
          end do
!
!  For each dst pt call the pnpoly to see if it is inside the src grd.
!
          do nx=1,ngrd_sw(mw)%Numx_swan
            do ny=1,ngrd_sw(mw)%Numy_swan
              PX=ngrd_sw(mw)%lon_rho_w(nx,ny)
              PY=ngrd_sw(mw)%lat_rho_w(nx,ny)
              INOUT=0
              call PNPOLY( PX, PY, XX, YY, N, INOUT )
              ngrd_sw(mw)%dst_mask(nx,ny)=                              &
     &          ngrd_sw(mw)%dst_mask(nx,ny)*                            &
     &          MIN(INOUT+1,1)
            enddo
          enddo
          deallocate(XX, YY)
!
!  Send the data to SCRIP routines
!
          write(mw_string,'(i1)')mw
          write(ma_string,'(i1)')ma
          interp_file1='atm'//trim(ma_string)//'_to_'//'wav'//          &
     &                       trim(mw_string)//'_weights.nc'
          interp_file2='unknown'
          map1_name='WRF to SWAN distwgt Mapping'
          map2_name='unknown'

          write(stdout,*)"============================================="
          write(stdout,*)"Begin mapping between the two grids"
          write(stdout,*)"---------------------------------------------"
          write(stdout,*)"The interp file is: ", interp_file1
          write(stdout,*)"The src grid is: ", wrf_grids(ma)
          write(stdout,*)"The dst grid is: ", swan_coord(mw)

          grid1_file=wrf_grids(ma)
          grid2_file=swan_coord(mw)

          counter_grid=counter_grid+1

          call scrip_package(grid1_file, grid2_file,                    &
     &                       ngrd_wr(ma)%we_size,                       &
     &                       ngrd_wr(ma)%sn_size,                       &
     &                       ngrd_wr(ma)%lon_rho_a,                     &
     &                       ngrd_wr(ma)%lat_rho_a,                     &
     &                       ngrd_wr(ma)%x_full_grid,                   &
     &                       ngrd_wr(ma)%y_full_grid,                   &
     &                       ngrd_sw(mw)%Numx_swan,                     &
     &                       ngrd_sw(mw)%Numy_swan,                     &
     &                       ngrd_sw(mw)%lon_rho_w,                     &
     &                       ngrd_sw(mw)%lat_rho_w,                     &
     &                       ngrd_sw(mw)%x_full_grid,                   &
     &                       ngrd_sw(mw)%y_full_grid,                   &
     &                       interp_file1, interp_file2,                &
     &                       map1_name, map2_name,                      &
     &                       ngrd_wr(ma)%mask_rho_a,                    &
     &                       ngrd_sw(mw)%dst_mask,                      &
     &                       dst_mask_unlim,                            &
     &                       counter_grid,                              &
     &                       Ngrids_comb_total,                         &
     &                       output_ncfile)

          deallocate(ngrd_sw(mw)%dst_mask)
          deallocate(dst_mask_unlim)
        end do
      end do
      do ma=1,Ngrids_wrf
        deallocate(ngrd_wr(ma)%src_mask)
      end do

      end subroutine atm2wav_mask

!======================================================================

      subroutine wav2atm_mask()

      implicit none

!     local variables
      real(dbl_kind), allocatable :: XX(:), YY(:)

!  Allocate the source mask
      do mw=1,Ngrids_swan
        nx=ngrd_sw(mw)%Numx_swan
        ny=ngrd_sw(mw)%Numy_swan
        allocate(ngrd_sw(mw)%src_mask(nx,ny))
      end do
!  Save the wav grid mask to src_mask
      do mw=1,Ngrids_swan
        do j=1,ngrd_sw(mw)%Numy_swan
          do i=1,ngrd_sw(mw)%Numx_swan
            ngrd_sw(mw)%src_mask(i,j)=ngrd_sw(mw)%mask_rho_w(i,j)
          end do
        end do
      end do
!  Mask out child portion of the src mask.
      do mw=1,Ngrids_swan-1
          do j=ngrd_sw(mw)%jstr_w,ngrd_sw(mw)%jend_w
            do i=ngrd_sw(mw)%istr_w,ngrd_sw(mw)%iend_w
              ngrd_sw(mw)%src_mask(i,j)=0
            end do
          end do
      end do
!  If a child grid, then set perimeter src_mask=0 because
!  this region comes from its parent.
      if (Ngrids_swan.gt.1) then
        do mw=2,Ngrids_swan
          i=ngrd_sw(mw)%Numx_swan
          do j=1,ngrd_sw(mw)%Numy_swan
            ngrd_sw(mw)%src_mask(1,j)=0
            ngrd_sw(mw)%src_mask(i,j)=0
          end do
          j=ngrd_sw(mw)%Numy_swan
          do i=1,ngrd_sw(mw)%Numx_swan
            ngrd_sw(mw)%src_mask(i,1)=0
            ngrd_sw(mw)%src_mask(i,j)=0
          end do
        end do
      end if
!  NOTICE piped swan and roms loops to end
      do ma=1,Ngrids_wrf
        do mw=1,Ngrids_swan
!
!  Allocate the destination mask
!
          nx=ngrd_wr(ma)%we_size
          ny=ngrd_wr(ma)%sn_size
          allocate(ngrd_wr(ma)%dst_mask(nx,ny))
          allocate(dst_mask_unlim(nx,ny))
!
!  Compute wrf grid dst mask. Start by using locations on dst grid
!  where src_mask=1.
!
          offset=MIN(mw-1,1)
          do nx=1,ngrd_wr(ma)%we_size
            do ny=1,ngrd_wr(ma)%sn_size
              dist_max=10e6
              Ikeep=1
              Jkeep=1
              xx2=ngrd_wr(ma)%lon_rho_a(nx,ny)
              yy2=ngrd_wr(ma)%lat_rho_a(nx,ny)
              do j=1+offset,ngrd_sw(mw)%Numy_swan-offset
                do i=1+offset,ngrd_sw(mw)%Numx_swan-offset
                  xx1=ngrd_sw(mw)%lon_rho_w(i,j)
                  yy1=ngrd_sw(mw)%lat_rho_w(i,j)
                  dlon = xx1-xx2
                  latrad1=abs(yy1*deg2rad)
                  latrad2=abs(yy2*deg2rad)
                  dep=cos(0.5*(latrad2+latrad1))*dlon
                  dlat=yy2-yy1
                  dist1=1852.0*60.0*sqrt(dlat**2+dep**2)
                  if(dist1<=dist_max)then
                    dist_max=dist1
                    Ikeep=i
                    Jkeep=j
                  endif
                enddo
              enddo
              ngrd_wr(ma)%dst_mask(nx,ny)=                              &
     &                                 ngrd_sw(mw)%src_mask(Ikeep,Jkeep)
            enddo
          enddo
!  Apply the wrf grid Land/Sea mask to dst_mask
          do j=1,ngrd_wr(ma)%sn_size
            do i=1,ngrd_wr(ma)%we_size
              ngrd_wr(ma)%dst_mask(i,j)=ngrd_wr(ma)%dst_mask(i,j)*      &
     &                                  ngrd_wr(ma)%mask_rho_a(i,j)
            end do
          end do
!  Remove points on dst mask that are outside of the src grid.
!  Create a vector of bounding src grid psi points.
          N=(ngrd_sw(mw)%Numx_swan+ngrd_sw(mw)%Numy_swan+2)*2-4
          allocate(XX(N), YY(N))
          count=0
          do i=1,ngrd_sw(mw)%Numy_swan+1
            count=count+1
            XX(count)=ngrd_sw(mw)%x_full_grid(1,i)
            YY(count)=ngrd_sw(mw)%y_full_grid(1,i)
          end do
          do i=2,ngrd_sw(mw)%Numx_swan+1
            count=count+1
            XX(count)=ngrd_sw(mw)%x_full_grid(i,ngrd_sw(mw)%Numy_swan+1)
            YY(count)=ngrd_sw(mw)%y_full_grid(i,ngrd_sw(mw)%Numy_swan+1)
          end do
          do i=ngrd_sw(mw)%Numy_swan,1,-1
            count=count+1
            XX(count)=ngrd_sw(mw)%x_full_grid(ngrd_sw(mw)%Numx_swan+1,i)
            YY(count)=ngrd_sw(mw)%y_full_grid(ngrd_sw(mw)%Numx_swan+1,i)
          end do
          do i=ngrd_sw(mw)%Numx_swan,2,-1
            count=count+1
            XX(count)=ngrd_sw(mw)%x_full_grid(i,1)
            YY(count)=ngrd_sw(mw)%y_full_grid(i,1)
          end do
!
!  For each dst pt call the pnpoly to see if it is inside the src grd.
!
          do nx=1,ngrd_wr(ma)%we_size
            do ny=1,ngrd_wr(ma)%sn_size
              PX=ngrd_wr(ma)%lon_rho_a(nx,ny)
              PY=ngrd_wr(ma)%lat_rho_a(nx,ny)
              INOUT=0
              call PNPOLY( PX, PY, XX, YY, N, INOUT )
              ngrd_wr(ma)%dst_mask(nx,ny)=                              &
     &          ngrd_wr(ma)%dst_mask(nx,ny)*                            &
     &          MIN(INOUT+1,1)
            enddo
          enddo
          deallocate(XX, YY)
!
!  Send the data to SCRIP routines
!
          write(ma_string,'(i1)')ma
          write(mw_string,'(i1)')mw
          interp_file1='wav'//trim(mw_string)//'_to_'//'atm'//          &
     &                       trim(ma_string)//'_weights.nc'
          interp_file2='unknown'
          map1_name='SWAN to WRF distwgt Mapping'
          map2_name='unknown'

          write(stdout,*)"============================================="
          write(stdout,*)"Begin mapping between the two grids"
          write(stdout,*)"---------------------------------------------"
          write(stdout,*)"The interp file is: ", interp_file1
          write(stdout,*)"The src grid is: ", swan_coord(mw)
          write(stdout,*)"The dst grid is: ", wrf_grids(ma)

          grid1_file=swan_coord(mw)
          grid2_file=wrf_grids(ma)

          counter_grid=counter_grid+1

          call scrip_package(grid1_file, grid2_file,                    &
     &                       ngrd_sw(mw)%Numx_swan,                     &
     &                       ngrd_sw(mw)%Numy_swan,                     &
     &                       ngrd_sw(mw)%lon_rho_w,                     &
     &                       ngrd_sw(mw)%lat_rho_w,                     &
     &                       ngrd_sw(mw)%x_full_grid,                   &
     &                       ngrd_sw(mw)%y_full_grid,                   &
     &                       ngrd_wr(ma)%we_size,                       &
     &                       ngrd_wr(ma)%sn_size,                       &
     &                       ngrd_wr(ma)%lon_rho_a,                     &
     &                       ngrd_wr(ma)%lat_rho_a,                     &
     &                       ngrd_wr(ma)%x_full_grid,                   &
     &                       ngrd_wr(ma)%y_full_grid,                   &
     &                       interp_file1, interp_file2,                &
     &                       map1_name, map2_name,                      &
     &                       ngrd_sw(mw)%mask_rho_w,                    &
     &                       ngrd_wr(ma)%dst_mask,                      &
     &                       dst_mask_unlim,                            &
     &                       counter_grid,                              &
     &                       Ngrids_comb_total,                         &
     &                       output_ncfile)

          deallocate(ngrd_wr(ma)%dst_mask)
          deallocate(dst_mask_unlim)
        end do
      end do
      do mw=1,Ngrids_swan
        deallocate(ngrd_sw(mw)%src_mask)
      end do

      end subroutine wav2atm_mask




!!!!!!!!NEWWWWWWWWWWw

!=======================================================================
      subroutine ocn2ww3_mask()

      implicit none

!     local variables
      real(dbl_kind), allocatable :: XX(:), YY(:)

!  Allocate the source mask
      do mo=1,Ngrids_roms
        nx=ngrd_rm(mo)%xi_size
        ny=ngrd_rm(mo)%eta_size
        allocate(ngrd_rm(mo)%src_mask(nx,ny))
      end do
!  Save the ocean grid mask to src_mask
      do mo=1,Ngrids_roms
        do j=1,ngrd_rm(mo)%eta_size
          do i=1,ngrd_rm(mo)%xi_size
            ngrd_rm(mo)%src_mask(i,j)=ngrd_rm(mo)%mask_rho_o(i,j)
          end do
        end do
      end do
!  Mask out child portion of the src mask.
      do mo=1,Ngrids_roms-1
          do j=ngrd_rm(mo)%jstr_o,ngrd_rm(mo)%jend_o
            do i=ngrd_rm(mo)%istr_o,ngrd_rm(mo)%iend_o
              ngrd_rm(mo)%src_mask(i,j)=0
            end do
          end do
      end do
!  If a child grid, then set perimeter src_mask=0 because
!  this region comes from its parent.
      if (Ngrids_roms.gt.1) then
        do mo=2,Ngrids_roms
          i=ngrd_rm(mo)%xi_size
          do j=1,ngrd_rm(mo)%eta_size
            ngrd_rm(mo)%src_mask(1,j)=0
            ngrd_rm(mo)%src_mask(i,j)=0
          end do
          j=ngrd_rm(mo)%eta_size
          do i=1,ngrd_rm(mo)%xi_size
            ngrd_rm(mo)%src_mask(i,1)=0
            ngrd_rm(mo)%src_mask(i,j)=0
          end do
        end do
      end if
!  NOTICE piped ww3 and roms loops to end
      allocate(do_adjust(Ngrids_ww3))
      do mw=1,Ngrids_ww3
        do_adjust(mw)=0
        do mo=1,Ngrids_roms
!
!  Allocate the destination mask
!
          nx=ngrd_w3(mw)%Numx_ww3
          ny=ngrd_w3(mw)%Numy_ww3
          allocate(ngrd_w3(mw)%dst_mask(nx,ny))
          allocate(dst_mask_unlim(nx,ny))
!
!  Compute ww3 grid dst mask. Start by using locations on dst grid
!  where src_mask=1.
!
          do nx=1,ngrd_w3(mw)%Numx_ww3
            do ny=1,ngrd_w3(mw)%Numy_ww3
              dist_max=10e6
              Ikeep=1
              Jkeep=1
              xx2=ngrd_w3(mw)%lon_rho_w(nx,ny)
              yy2=ngrd_w3(mw)%lat_rho_w(nx,ny)
              do j=1,ngrd_rm(mo)%eta_size
                do i=1,ngrd_rm(mo)%xi_size
                  dist1=sqrt((ngrd_rm(mo)%lon_rho_o(i,j)-xx2)**2+       &
     &                       (ngrd_rm(mo)%lat_rho_o(i,j)-yy2)**2)
!                 xx1=ngrd_rm(mo)%lon_rho_o(i,j)
!                 yy1=ngrd_rm(mo)%lat_rho_o(i,j)
!                 dlon = xx1-xx2
!                 latrad1=abs(yy1*deg2rad)
!                 latrad2=abs(yy2*deg2rad)
!                 dep=cos(0.5*(latrad2+latrad1))*dlon
!                 dlat=yy2-yy1
!                 dist1=1852.0*60.0*sqrt(dlat**2+dep**2)
                  if(dist1<=dist_max)then
                    dist_max=dist1
                    Ikeep=i
                    Jkeep=j
                  endif
                enddo
              enddo
              ngrd_w3(mw)%dst_mask(nx,ny)=                              &
     &                                 ngrd_rm(mo)%src_mask(Ikeep,Jkeep)
            enddo
          enddo
!  Apply the ww3 grid Land/Sea mask to dst_mask
          do j=1,ngrd_w3(mw)%Numy_ww3
            do i=1,ngrd_w3(mw)%Numx_ww3
              ngrd_w3(mw)%dst_mask(i,j)=ngrd_w3(mw)%dst_mask(i,j)*      &
     &                                  ngrd_w3(mw)%mask_rho_w(i,j)
            end do
          end do
!  Remove points on dst mask that are outside of the src grid.
!  Create a vector of bounding src grid psi points.
          N=(ngrd_rm(mo)%xi_size+ngrd_rm(mo)%eta_size+2)*2-4
          allocate(XX(N), YY(N))
          count=0
          do i=1,ngrd_rm(mo)%eta_size+1
            count=count+1
            XX(count)=ngrd_rm(mo)%x_full_grid(1,i)
            YY(count)=ngrd_rm(mo)%y_full_grid(1,i)
          end do
          do i=2,ngrd_rm(mo)%xi_size+1
            count=count+1
            XX(count)=ngrd_rm(mo)%x_full_grid(i,ngrd_rm(mo)%eta_size+1)
            YY(count)=ngrd_rm(mo)%y_full_grid(i,ngrd_rm(mo)%eta_size+1)
          end do
          do i=ngrd_rm(mo)%eta_size,1,-1
            count=count+1
            XX(count)=ngrd_rm(mo)%x_full_grid(ngrd_rm(mo)%xi_size+1,i)
            YY(count)=ngrd_rm(mo)%y_full_grid(ngrd_rm(mo)%xi_size+1,i)
          end do
          do i=ngrd_rm(mo)%xi_size,2,-1
            count=count+1
            XX(count)=ngrd_rm(mo)%x_full_grid(i,1)
            YY(count)=ngrd_rm(mo)%y_full_grid(i,1)
          end do
!
!  For each dst pt call the pnpoly to see if it is inside the src grd.
!
          do nx=1,ngrd_w3(mw)%Numx_ww3
            do ny=1,ngrd_w3(mw)%Numy_ww3
              PX=ngrd_w3(mw)%lon_rho_w(nx,ny)
              PY=ngrd_w3(mw)%lat_rho_w(nx,ny)
              INOUT=0
              call PNPOLY( PX, PY, XX, YY, N, INOUT )
              ngrd_w3(mw)%dst_mask(nx,ny)=                              &
     &          ngrd_w3(mw)%dst_mask(nx,ny)*                            &
     &          MIN(INOUT+1,1)
            enddo
          enddo
          deallocate(XX, YY)
!
!  Compare src and dst lon points and see if points are coincident. 
!  The test looks for less than 2.0e-7 rads =  1.1459e-05 degs =~ 1.3m
!  If so, then offset the dst by 2.0e-7 rads
!
            tol=1.1459e-05
            do i=1,ngrd_rm(mo)%xi_size
              if (do_adjust(mw).eq.1) exit
              do j=1,ngrd_rm(mo)%eta_size
                if (do_adjust(mw).eq.1) exit
                lons=ngrd_rm(mo)%lon_rho_o(i,j)
                lats=ngrd_rm(mo)%lat_rho_o(i,j)
                do nx=1,ngrd_w3(mw)%Numx_ww3
                  if (do_adjust(mw).eq.1) exit
                  do ny=1,ngrd_w3(mw)%Numy_ww3
                    lond=ngrd_w3(mw)%lon_rho_w(nx,ny)
                    latd=ngrd_w3(mw)%lat_rho_w(nx,ny)
                    cff=min(abs(lons-lond),abs(lats-latd))
                    if (cff.le.tol) then
                      do_adjust(mw)=1
                  write(*,*) 'Do adjust for coincident points grid ', mw
                      exit
                    end if
                  enddo
                enddo
              enddo
            enddo
            if (do_adjust(mw).eq.1) then
              do nx=1,ngrd_w3(mw)%Numx_ww3
                do ny=1,ngrd_w3(mw)%Numy_ww3
                  ngrd_w3(mw)%lon_rho_w(nx,ny)=                         &
     &                                  ngrd_w3(mw)%lon_rho_w(nx,ny)+tol
                  ngrd_w3(mw)%lat_rho_w(nx,ny)=                         &
     &                                  ngrd_w3(mw)%lat_rho_w(nx,ny)+tol
                enddo
              enddo
              do nx=1,ngrd_w3(mw)%Numx_ww3+1
                do ny=1,ngrd_w3(mw)%Numy_ww3+1
                  ngrd_w3(mw)%x_full_grid(nx,ny)=                       &
     &                                ngrd_w3(mw)%x_full_grid(nx,ny)+tol
                  ngrd_w3(mw)%y_full_grid(nx,ny)=                       &
     &                                ngrd_w3(mw)%y_full_grid(nx,ny)+tol
                enddo
              enddo
            end if
!
!  Send the data to SCRIP routines
!
          write(mw_string,'(i1)')mw
          write(mo_string,'(i1)')mo
          interp_file1='ocn'//trim(mo_string)//'_to_'//'wav'//          &
     &                       trim(mw_string)//'_weights.nc'
          interp_file2='unknown'
          map1_name='ROMS to WW3 distwgt Mapping'
          map2_name='unknown'

          write(stdout,*)"============================================="
          write(stdout,*)"Begin mapping between the two grids"
          write(stdout,*)"---------------------------------------------"
          write(stdout,*)"The interp file is: ", interp_file1
          write(stdout,*)"The src grid is: ", roms_grids(mo)
          write(stdout,*)"The dst grid is: ", ww3_xcoord(mw)

          grid1_file=roms_grids(mo)
          grid2_file=ww3_xcoord(mw)

          counter_grid=counter_grid+1

          call scrip_package(grid1_file, grid2_file,                    &
     &                       ngrd_rm(mo)%xi_size,                       &
     &                       ngrd_rm(mo)%eta_size,                      &
     &                       ngrd_rm(mo)%lon_rho_o,                     &
     &                       ngrd_rm(mo)%lat_rho_o,                     &
     &                       ngrd_rm(mo)%x_full_grid,                   &
     &                       ngrd_rm(mo)%y_full_grid,                   &
     &                       ngrd_w3(mw)%Numx_ww3,                      &
     &                       ngrd_w3(mw)%Numy_ww3,                      &
     &                       ngrd_w3(mw)%lon_rho_w,                     &
     &                       ngrd_w3(mw)%lat_rho_w,                     &
     &                       ngrd_w3(mw)%x_full_grid,                   &
     &                       ngrd_w3(mw)%y_full_grid,                   &
     &                       interp_file1, interp_file2,                &
     &                       map1_name, map2_name,                      &
     &                       ngrd_rm(mo)%mask_rho_o,                    &
     &                       ngrd_w3(mw)%dst_mask,                      &
     &                       dst_mask_unlim,                            &
     &                       counter_grid,                              &
     &                       Ngrids_comb_total,                         &
     &                       output_ncfile)
!
!  Undo the grid adjust.
          if ((do_adjust(mw).eq.1).and.(mo.lt.Ngrids_roms)) then
            do nx=1,ngrd_w3(mw)%Numx_ww3
              do ny=1,ngrd_w3(mw)%Numy_ww3
                ngrd_w3(mw)%lon_rho_w(nx,ny)=                           &
     &                                ngrd_w3(mw)%lon_rho_w(nx,ny)-tol
                ngrd_w3(mw)%lat_rho_w(nx,ny)=                           &
     &                                ngrd_w3(mw)%lat_rho_w(nx,ny)-tol
              enddo
            enddo
            do nx=1,ngrd_w3(mw)%Numx_ww3+1
              do ny=1,ngrd_w3(mw)%Numy_ww3+1
                ngrd_w3(mw)%x_full_grid(nx,ny)=                         &
     &                              ngrd_w3(mw)%x_full_grid(nx,ny)-tol
                ngrd_w3(mw)%y_full_grid(nx,ny)=                         &
     &                              ngrd_w3(mw)%y_full_grid(nx,ny)-tol
              enddo
            enddo
          end if
          do_adjust(mw)=0
!
          deallocate(ngrd_w3(mw)%dst_mask)
          deallocate(dst_mask_unlim)
        end do
      end do
      deallocate(do_adjust)
      do mo=1,Ngrids_roms
        deallocate(ngrd_rm(mo)%src_mask)
      end do

      end subroutine ocn2ww3_mask

!======================================================================
!
      subroutine ww32ocn_mask()

      implicit none

!     local variables
      real(dbl_kind), allocatable :: XX(:), YY(:)

!  Allocate the source mask
      do mw=1,Ngrids_ww3
        nx=ngrd_w3(mw)%Numx_ww3
        ny=ngrd_w3(mw)%Numy_ww3
        allocate(ngrd_w3(mw)%src_mask(nx,ny))
      end do
!  Save the wav grid mask to src_mask
      do mw=1,Ngrids_ww3
        do j=1,ngrd_w3(mw)%Numy_ww3
          do i=1,ngrd_w3(mw)%Numx_ww3
            ngrd_w3(mw)%src_mask(i,j)=ngrd_w3(mw)%mask_rho_w(i,j)
          end do
        end do
      end do
!  Mask out child portion of the src mask.
      do mw=1,Ngrids_ww3-1
        do j=ngrd_w3(mw)%jstr_w,ngrd_w3(mw)%jend_w
          do i=ngrd_w3(mw)%istr_w,ngrd_w3(mw)%iend_w
            ngrd_w3(mw)%src_mask(i,j)=0
           end do
        end do
      end do
!  If a child grid, then set perimeter src_mask=0 because
!  this region comes from its parent.
      if (Ngrids_ww3.gt.1) then
        do mw=2,Ngrids_ww3
          i=ngrd_w3(mw)%Numx_ww3
          do j=1,ngrd_w3(mw)%Numy_ww3
            ngrd_w3(mw)%src_mask(1,j)=0
            ngrd_w3(mw)%src_mask(i,j)=0
          end do
          j=ngrd_w3(mw)%Numy_ww3
          do i=1,ngrd_w3(mw)%Numx_ww3
            ngrd_w3(mw)%src_mask(i,1)=0
            ngrd_w3(mw)%src_mask(i,j)=0
          end do
        end do
      end if
!  NOTICE piped ww3 and roms loops to end
      do mo=1,Ngrids_roms
        do mw=1,Ngrids_ww3
!
!  Allocate the destination mask
!
          nx=ngrd_rm(mo)%xi_size
          ny=ngrd_rm(mo)%eta_size
          allocate(ngrd_rm(mo)%dst_mask(nx,ny))
          allocate(dst_mask_unlim(nx,ny))
!
!  Compute roms grid dst mask. Start by using locations on dst grid
!  where src_mask=1.
!
          do nx=1,ngrd_rm(mo)%xi_size
            do ny=1,ngrd_rm(mo)%eta_size
              dist_max=10e6
              Ikeep=1
              Jkeep=1
              xx2=ngrd_rm(mo)%lon_rho_o(nx,ny)
              yy2=ngrd_rm(mo)%lat_rho_o(nx,ny)
              do j=1,ngrd_w3(mw)%Numy_ww3
                do i=1,ngrd_w3(mw)%Numx_ww3
                  dist1=sqrt((ngrd_w3(mw)%lon_rho_w(i,j)-xx2)**2+       &
     &                       (ngrd_w3(mw)%lat_rho_w(i,j)-yy2)**2)
!                 xx1=ngrd_sw(mw)%lon_rho_w(i,j)
!                 yy1=ngrd_sw(mw)%lat_rho_w(i,j)
!                 dlon = xx1-xx2
!                 latrad1=abs(yy1*deg2rad)
!                 latrad2=abs(yy2*deg2rad)
!                 dep=cos(0.5*(latrad2+latrad1))*dlon
!                 dlat=yy2-yy1
!                 dist1=1852.0*60.0*sqrt(dlat**2+dep**2)
                  if(dist1<=dist_max)then
                    dist_max=dist1
                    Ikeep=i
                    Jkeep=j
                  endif
                enddo
              enddo
              ngrd_rm(mo)%dst_mask(nx,ny)=                              &
     &                                 ngrd_w3(mw)%src_mask(Ikeep,Jkeep)
            enddo
          enddo
!  Apply the roms grid Land/Sea mask to dst_mask
          do j=1,ngrd_rm(mo)%eta_size
            do i=1,ngrd_rm(mo)%xi_size
              ngrd_rm(mo)%dst_mask(i,j)=ngrd_rm(mo)%dst_mask(i,j)*      &
     &                                  ngrd_rm(mo)%mask_rho_o(i,j)
            end do
          end do
!  Remove points on dst mask that are outside of the src grid.
!  Create a vector of bounding src grid psi points.
          N=(ngrd_w3(mw)%Numx_ww3+ngrd_w3(mw)%Numy_ww3+2)*2-4
          allocate(XX(N), YY(N))
          count=0
          do i=1,ngrd_w3(mw)%Numy_ww3+1
            count=count+1
            XX(count)=ngrd_w3(mw)%x_full_grid(1,i)
            YY(count)=ngrd_w3(mw)%y_full_grid(1,i)
          end do
          do i=2,ngrd_w3(mw)%Numx_ww3+1
            count=count+1
            XX(count)=ngrd_w3(mw)%x_full_grid(i,ngrd_w3(mw)%Numy_ww3+1)
            YY(count)=ngrd_w3(mw)%y_full_grid(i,ngrd_w3(mw)%Numy_ww3+1)
          end do
          do i=ngrd_w3(mw)%Numy_ww3,1,-1
            count=count+1
            XX(count)=ngrd_w3(mw)%x_full_grid(ngrd_w3(mw)%Numx_ww3+1,i)
            YY(count)=ngrd_w3(mw)%y_full_grid(ngrd_w3(mw)%Numx_ww3+1,i)
          end do
          do i=ngrd_w3(mw)%Numx_ww3,2,-1
            count=count+1
            XX(count)=ngrd_w3(mw)%x_full_grid(i,1)
            YY(count)=ngrd_w3(mw)%y_full_grid(i,1)
          end do
!
!  For each dst pt call the pnpoly to see if it is inside the src grd.
!
          do nx=1,ngrd_rm(mo)%xi_size
            do ny=1,ngrd_rm(mo)%eta_size
              PX=ngrd_rm(mo)%lon_rho_o(nx,ny)
              PY=ngrd_rm(mo)%lat_rho_o(nx,ny)
              INOUT=0
              call PNPOLY( PX, PY, XX, YY, N, INOUT )
              ngrd_rm(mo)%dst_mask(nx,ny)=                              &
     &          ngrd_rm(mo)%dst_mask(nx,ny)*                            &
     &          MIN(INOUT+1,1)
            enddo
          enddo
          deallocate(XX, YY)
!
!  Send the data to SCRIP routines
!
          write(mo_string,'(i1)')mo
          write(mw_string,'(i1)')mw
          interp_file1='wav'//trim(mw_string)//'_to_'//'ocn'//          &
     &                       trim(mo_string)//'_weights.nc'
          interp_file2='unknown'
          map1_name='WW3 to ROMS distwgt Mapping'
          map2_name='unknown'

          write(stdout,*)"============================================="
          write(stdout,*)"Begin mapping between the two grids"
          write(stdout,*)"---------------------------------------------"
          write(stdout,*)"The interp file is: ", interp_file1
          write(stdout,*)"The src grid is: ", ww3_xcoord(mw)
          write(stdout,*)"The dst grid is: ", roms_grids(mo)

          grid1_file=ww3_xcoord(mw)
          grid2_file=roms_grids(mo)

          counter_grid=counter_grid+1

          call scrip_package(grid1_file, grid2_file,                    &
     &                       ngrd_w3(mw)%Numx_ww3,                      &
     &                       ngrd_w3(mw)%Numy_ww3,                      &
     &                       ngrd_w3(mw)%lon_rho_w,                     &
     &                       ngrd_w3(mw)%lat_rho_w,                     &
     &                       ngrd_w3(mw)%x_full_grid,                   &
     &                       ngrd_w3(mw)%y_full_grid,                   &
     &                       ngrd_rm(mo)%xi_size,                       &
     &                       ngrd_rm(mo)%eta_size,                      &
     &                       ngrd_rm(mo)%lon_rho_o,                     &
     &                       ngrd_rm(mo)%lat_rho_o,                     &
     &                       ngrd_rm(mo)%x_full_grid,                   &
     &                       ngrd_rm(mo)%y_full_grid,                   &
     &                       interp_file1, interp_file2,                &
     &                       map1_name, map2_name,                      &
     &                       ngrd_w3(mw)%mask_rho_w,                    &
     &                       ngrd_rm(mo)%dst_mask,                      &
     &                       dst_mask_unlim,                            &
     &                       counter_grid,                              &
     &                       Ngrids_comb_total,                         &
     &                       output_ncfile)

          deallocate(ngrd_rm(mo)%dst_mask)
          deallocate(dst_mask_unlim)
        end do
      end do
      do mw=1,Ngrids_ww3
        deallocate(ngrd_w3(mw)%src_mask)
      end do

      end subroutine ww32ocn_mask

      
!======================================================================

      subroutine atm2ww3_mask()

      implicit none

!     local variables
      real(dbl_kind), allocatable :: XX(:), YY(:)

!  Allocate the source mask
      do ma=1,Ngrids_wrf
        nx=ngrd_wr(ma)%we_size
        ny=ngrd_wr(ma)%sn_size
        allocate(ngrd_wr(ma)%src_mask(nx,ny))
      end do
!  Save the atm grid mask to src_mask
      do ma=1,Ngrids_wrf
        do j=1,ngrd_wr(ma)%sn_size
          do i=1,ngrd_wr(ma)%we_size
            ngrd_wr(ma)%src_mask(i,j)=ngrd_wr(ma)%mask_rho_a(i,j)
          end do
        end do
      end do
!  Mask out child portion of the src mask.
      do ma=1,Ngrids_wrf-1
        if (wrf_grids(ma+1)/="moving") then
          do j=ngrd_wr(ma)%jstr_a,ngrd_wr(ma)%jend_a
            do i=ngrd_wr(ma)%istr_a,ngrd_wr(ma)%iend_a
              ngrd_wr(ma)%src_mask(i,j)=0
            end do
          end do
        end if
      end do
!
!  NOTICE piped wrf and ww3 loops to end
      do mw=1,Ngrids_ww3
        do ma=1,Ngrids_wrf
!
!  Allocate the destination mask
!
          nx=ngrd_w3(mw)%Numx_ww3
          ny=ngrd_w3(mw)%Numy_ww3
          allocate(ngrd_w3(mw)%dst_mask(nx,ny))
          allocate(dst_mask_unlim(nx,ny))
!
!  Compute ww3 grid dst mask. Start by using locations on dst grid
!  where src_mask=1.
!
          do nx=1,ngrd_w3(mw)%Numx_ww3
            do ny=1,ngrd_w3(mw)%Numy_ww3
              dist_max=10e6
              Ikeep=1
              Jkeep=1
              xx2=ngrd_w3(mw)%lon_rho_w(nx,ny)
              yy2=ngrd_w3(mw)%lat_rho_w(nx,ny)
              do j=1,ngrd_wr(ma)%sn_size
                do i=1,ngrd_wr(ma)%we_size
                  xx1=ngrd_wr(ma)%lon_rho_a(i,j)
                  yy1=ngrd_wr(ma)%lat_rho_a(i,j)
                  dlon = xx1-xx2
                  latrad1=abs(yy1*deg2rad)
                  latrad2=abs(yy2*deg2rad)
                  dep=cos(0.5*(latrad2+latrad1))*dlon
                  dlat=yy2-yy1
                  dist1=1852.0*60.0*sqrt(dlat**2+dep**2)
                  if(dist1<=dist_max)then
                    dist_max=dist1
                    Ikeep=i
                    Jkeep=j
                  endif
                enddo
              enddo
              ngrd_w3(mw)%dst_mask(nx,ny)=                              &
     &                                 ngrd_wr(ma)%src_mask(Ikeep,Jkeep)
            enddo
          enddo
!  Apply the ww3 grid Land/Sea mask to dst_mask
          do j=1,ngrd_w3(mw)%Numy_ww3
            do i=1,ngrd_w3(mw)%Numx_ww3
              ngrd_w3(mw)%dst_mask(i,j)=ngrd_w3(mw)%dst_mask(i,j)*      &
     &                                  ngrd_w3(mw)%mask_rho_w(i,j)
            end do
          end do
!  Remove points on dst mask that are outside of the src grid.
!  Create a vector of bounding src grid psi points.
          N=(ngrd_wr(ma)%we_size+ngrd_wr(ma)%sn_size+2)*2-4
          allocate(XX(N), YY(N))
          count=0
          do i=1,ngrd_wr(ma)%sn_size+1
            count=count+1
            XX(count)=ngrd_wr(ma)%x_full_grid(1,i)
            YY(count)=ngrd_wr(ma)%y_full_grid(1,i)
          end do
          do i=2,ngrd_wr(ma)%we_size+1
            count=count+1
            XX(count)=ngrd_wr(ma)%x_full_grid(i,ngrd_wr(ma)%sn_size+1)
            YY(count)=ngrd_wr(ma)%y_full_grid(i,ngrd_wr(ma)%sn_size+1)
          end do
          do i=ngrd_wr(ma)%sn_size,1,-1
            count=count+1
            XX(count)=ngrd_wr(ma)%x_full_grid(ngrd_wr(ma)%we_size+1,i)
            YY(count)=ngrd_wr(ma)%y_full_grid(ngrd_wr(ma)%we_size+1,i)
          end do
          do i=ngrd_wr(ma)%we_size,2,-1
            count=count+1
            XX(count)=ngrd_wr(ma)%x_full_grid(i,1)
            YY(count)=ngrd_wr(ma)%y_full_grid(i,1)
          end do
!
!  For each dst pt call the pnpoly to see if it is inside the src grd.
!
          do nx=1,ngrd_w3(mw)%Numx_ww3
            do ny=1,ngrd_w3(mw)%Numy_ww3
              PX=ngrd_w3(mw)%lon_rho_w(nx,ny)
              PY=ngrd_w3(mw)%lat_rho_w(nx,ny)
              INOUT=0
              call PNPOLY( PX, PY, XX, YY, N, INOUT )
              ngrd_w3(mw)%dst_mask(nx,ny)=                              &
     &          ngrd_w3(mw)%dst_mask(nx,ny)*                            &
     &          MIN(INOUT+1,1)
            enddo
          enddo
          deallocate(XX, YY)
!
!  Send the data to SCRIP routines
!
          write(mw_string,'(i1)')mw
          write(ma_string,'(i1)')ma
          interp_file1='atm'//trim(ma_string)//'_to_'//'wav'//          &
     &                       trim(mw_string)//'_weights.nc'
          interp_file2='unknown'
          map1_name='WRF to WW3 distwgt Mapping'
          map2_name='unknown'

          write(stdout,*)"============================================="
          write(stdout,*)"Begin mapping between the two grids"
          write(stdout,*)"---------------------------------------------"
          write(stdout,*)"The interp file is: ", interp_file1
          write(stdout,*)"The src grid is: ", wrf_grids(ma)
          write(stdout,*)"The dst grid is: ", ww3_xcoord(mw)

          grid1_file=wrf_grids(ma)
          grid2_file=ww3_xcoord(mw)

          counter_grid=counter_grid+1

          call scrip_package(grid1_file, grid2_file,                    &
     &                       ngrd_wr(ma)%we_size,                       &
     &                       ngrd_wr(ma)%sn_size,                       &
     &                       ngrd_wr(ma)%lon_rho_a,                     &
     &                       ngrd_wr(ma)%lat_rho_a,                     &
     &                       ngrd_wr(ma)%x_full_grid,                   &
     &                       ngrd_wr(ma)%y_full_grid,                   &
     &                       ngrd_w3(mw)%Numx_ww3,                      &
     &                       ngrd_w3(mw)%Numy_ww3,                      &
     &                       ngrd_w3(mw)%lon_rho_w,                     &
     &                       ngrd_w3(mw)%lat_rho_w,                     &
     &                       ngrd_w3(mw)%x_full_grid,                   &
     &                       ngrd_w3(mw)%y_full_grid,                   &
     &                       interp_file1, interp_file2,                &
     &                       map1_name, map2_name,                      &
     &                       ngrd_wr(ma)%mask_rho_a,                    &
     &                       ngrd_w3(mw)%dst_mask,                      &
     &                       dst_mask_unlim,                            &
     &                       counter_grid,                              &
     &                       Ngrids_comb_total,                         &
     &                       output_ncfile)

          deallocate(ngrd_w3(mw)%dst_mask)
          deallocate(dst_mask_unlim)
        end do
      end do
      do ma=1,Ngrids_wrf
        deallocate(ngrd_wr(ma)%src_mask)
      end do

      end subroutine atm2ww3_mask

!======================================================================

      subroutine ww32atm_mask()

      implicit none

!     local variables
      real(dbl_kind), allocatable :: XX(:), YY(:)

!  Allocate the source mask
      do mw=1,Ngrids_ww3
        nx=ngrd_w3(mw)%Numx_ww3
        ny=ngrd_w3(mw)%Numy_ww3
        allocate(ngrd_w3(mw)%src_mask(nx,ny))
      end do
!  Save the wav grid mask to src_mask
      do mw=1,Ngrids_ww3
        do j=1,ngrd_w3(mw)%Numy_ww3
          do i=1,ngrd_w3(mw)%Numx_ww3
            ngrd_w3(mw)%src_mask(i,j)=ngrd_w3(mw)%mask_rho_w(i,j)
          end do
        end do
      end do
!  Mask out child portion of the src mask.
      do mw=1,Ngrids_ww3-1
          do j=ngrd_w3(mw)%jstr_w,ngrd_w3(mw)%jend_w
            do i=ngrd_w3(mw)%istr_w,ngrd_w3(mw)%iend_w
              ngrd_w3(mw)%src_mask(i,j)=0
            end do
          end do
      end do
!  If a child grid, then set perimeter src_mask=0 because
!  this region comes from its parent.
      if (Ngrids_ww3.gt.1) then
        do mw=2,Ngrids_ww3
          i=ngrd_w3(mw)%Numx_ww3
          do j=1,ngrd_w3(mw)%Numy_ww3
            ngrd_w3(mw)%src_mask(1,j)=0
            ngrd_w3(mw)%src_mask(i,j)=0
          end do
          j=ngrd_w3(mw)%Numy_ww3
          do i=1,ngrd_w3(mw)%Numx_ww3
            ngrd_w3(mw)%src_mask(i,1)=0
            ngrd_w3(mw)%src_mask(i,j)=0
          end do
        end do
      end if
!  NOTICE piped ww3 and roms loops to end
      do ma=1,Ngrids_wrf
        do mw=1,Ngrids_ww3
!
!  Allocate the destination mask
!
          nx=ngrd_wr(ma)%we_size
          ny=ngrd_wr(ma)%sn_size
          allocate(ngrd_wr(ma)%dst_mask(nx,ny))
          allocate(dst_mask_unlim(nx,ny))
!
!  Compute wrf grid dst mask. Start by using locations on dst grid
!  where src_mask=1.
!
          offset=MIN(mw-1,1)
          do nx=1,ngrd_wr(ma)%we_size
            do ny=1,ngrd_wr(ma)%sn_size
              dist_max=10e6
              Ikeep=1
              Jkeep=1
              xx2=ngrd_wr(ma)%lon_rho_a(nx,ny)
              yy2=ngrd_wr(ma)%lat_rho_a(nx,ny)
              do j=1+offset,ngrd_w3(mw)%Numy_ww3-offset
                do i=1+offset,ngrd_w3(mw)%Numx_ww3-offset
                  xx1=ngrd_w3(mw)%lon_rho_w(i,j)
                  yy1=ngrd_w3(mw)%lat_rho_w(i,j)
                  dlon = xx1-xx2
                  latrad1=abs(yy1*deg2rad)
                  latrad2=abs(yy2*deg2rad)
                  dep=cos(0.5*(latrad2+latrad1))*dlon
                  dlat=yy2-yy1
                  dist1=1852.0*60.0*sqrt(dlat**2+dep**2)
                  if(dist1<=dist_max)then
                    dist_max=dist1
                    Ikeep=i
                    Jkeep=j
                  endif
                enddo
              enddo
              ngrd_wr(ma)%dst_mask(nx,ny)=                              &
     &                                 ngrd_w3(mw)%src_mask(Ikeep,Jkeep)
            enddo
          enddo
!  Apply the wrf grid Land/Sea mask to dst_mask
          do j=1,ngrd_wr(ma)%sn_size
            do i=1,ngrd_wr(ma)%we_size
              ngrd_wr(ma)%dst_mask(i,j)=ngrd_wr(ma)%dst_mask(i,j)*      &
     &                                  ngrd_wr(ma)%mask_rho_a(i,j)
            end do
          end do
!  Remove points on dst mask that are outside of the src grid.
!  Create a vector of bounding src grid psi points.
          N=(ngrd_w3(mw)%Numx_ww3+ngrd_w3(mw)%Numy_ww3+2)*2-4
          allocate(XX(N), YY(N))
          count=0
          do i=1,ngrd_w3(mw)%Numy_ww3+1
            count=count+1
            XX(count)=ngrd_w3(mw)%x_full_grid(1,i)
            YY(count)=ngrd_w3(mw)%y_full_grid(1,i)
          end do
          do i=2,ngrd_w3(mw)%Numx_ww3+1
            count=count+1
            XX(count)=ngrd_w3(mw)%x_full_grid(i,ngrd_w3(mw)%Numy_ww3+1)
            YY(count)=ngrd_w3(mw)%y_full_grid(i,ngrd_w3(mw)%Numy_ww3+1)
          end do
          do i=ngrd_w3(mw)%Numy_ww3,1,-1
            count=count+1
            XX(count)=ngrd_w3(mw)%x_full_grid(ngrd_w3(mw)%Numx_ww3+1,i)
            YY(count)=ngrd_w3(mw)%y_full_grid(ngrd_w3(mw)%Numx_ww3+1,i)
          end do
          do i=ngrd_w3(mw)%Numx_ww3,2,-1
            count=count+1
            XX(count)=ngrd_w3(mw)%x_full_grid(i,1)
            YY(count)=ngrd_w3(mw)%y_full_grid(i,1)
          end do
!
!  For each dst pt call the pnpoly to see if it is inside the src grd.
!
          do nx=1,ngrd_wr(ma)%we_size
            do ny=1,ngrd_wr(ma)%sn_size
              PX=ngrd_wr(ma)%lon_rho_a(nx,ny)
              PY=ngrd_wr(ma)%lat_rho_a(nx,ny)
              INOUT=0
              call PNPOLY( PX, PY, XX, YY, N, INOUT )
              ngrd_wr(ma)%dst_mask(nx,ny)=                              &
     &          ngrd_wr(ma)%dst_mask(nx,ny)*                            &
     &          MIN(INOUT+1,1)
            enddo
          enddo
          deallocate(XX, YY)
!
!  Send the data to SCRIP routines
!
          write(ma_string,'(i1)')ma
          write(mw_string,'(i1)')mw
          interp_file1='wav'//trim(mw_string)//'_to_'//'atm'//          &
     &                       trim(ma_string)//'_weights.nc'
          interp_file2='unknown'
          map1_name='WW3 to WRF distwgt Mapping'
          map2_name='unknown'

          write(stdout,*)"============================================="
          write(stdout,*)"Begin mapping between the two grids"
          write(stdout,*)"---------------------------------------------"
          write(stdout,*)"The interp file is: ", interp_file1
          write(stdout,*)"The src grid is: ", ww3_xcoord(mw)
          write(stdout,*)"The dst grid is: ", wrf_grids(ma)

          grid1_file=ww3_xcoord(mw)
          grid2_file=wrf_grids(ma)

          counter_grid=counter_grid+1

          call scrip_package(grid1_file, grid2_file,                    &
     &                       ngrd_w3(mw)%Numx_ww3,                      &
     &                       ngrd_w3(mw)%Numy_ww3,                      &
     &                       ngrd_w3(mw)%lon_rho_w,                     &
     &                       ngrd_w3(mw)%lat_rho_w,                     &
     &                       ngrd_w3(mw)%x_full_grid,                   &
     &                       ngrd_w3(mw)%y_full_grid,                   &
     &                       ngrd_wr(ma)%we_size,                       &
     &                       ngrd_wr(ma)%sn_size,                       &
     &                       ngrd_wr(ma)%lon_rho_a,                     &
     &                       ngrd_wr(ma)%lat_rho_a,                     &
     &                       ngrd_wr(ma)%x_full_grid,                   &
     &                       ngrd_wr(ma)%y_full_grid,                   &
     &                       interp_file1, interp_file2,                &
     &                       map1_name, map2_name,                      &
     &                       ngrd_w3(mw)%mask_rho_w,                    &
     &                       ngrd_wr(ma)%dst_mask,                      &
     &                       dst_mask_unlim,                            &
     &                       counter_grid,                              &
     &                       Ngrids_comb_total,                         &
     &                       output_ncfile)

          deallocate(ngrd_wr(ma)%dst_mask)
          deallocate(dst_mask_unlim)
        end do
      end do
      do mw=1,Ngrids_ww3
        deallocate(ngrd_w3(mw)%src_mask)
      end do

      end subroutine ww32atm_mask


!======================================================================
!                                                                       
         SUBROUTINE PNPOLY( PX, PY, XX, YY, N, INOUT )                     
!
!        https://www.ecse.rpi.edu/~wrf/Research/Short_Notes/pnpoly.html
!        Copyright (c) 1970-2003, Wm. Randolph Franklin
!        Permission is hereby granted, free of charge, to any person 
!        obtaining a copy of this software and associated documentation
!        files (the "Software"), to deal in the Software without 
!        restriction, including without limitation the rights to use, 
!        copy, modify, merge, publish, distribute, sublicense, and/or 
!        sell copies of the Software, and to permit persons to whom the
!        Software is furnished to do so, subject to the following 
!        conditions: 
!        1.Redistributions of source code must retain the above 
!          copyright notice, this list of conditions and the following
!          disclaimers. 
!        2.Redistributions in binary form must reproduce the above
!          copyright notice in the documentation and/or other materials
!          provided with the distribution. 
!        3.The name of W. Randolph Franklin may not be used to endorse
!          or promote products derived from this Software without 
!          specific prior written permission. 
!
!        THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
!        EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
!        OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND 
!        NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT 
!        HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, 
!        WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING 
!        FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR 
!        OTHER DEALINGS IN THE SOFTWARE. 
!                                                                       
!        PURPOSE                                                        
!           TO DETERMINE WHETHER A POINT IS INSIDE A POLYGON            
!                                                                       
!        USAGE                                                          
!           CALL PNPOLY (PX, PY, XX, YY, N, INOUT )                     
!                                                                       
!        DESCRIPTION OF THE PARAMETERS                                  
!           PX      - X-COORDINATE OF POINT IN QUESTION.                
!           PY      - Y-COORDINATE OF POINT IN QUESTION.                
!           XX      - N LONG VECTOR CONTAINING X-COORDINATES OF         
!                     VERTICES OF POLYGON.                              
!           YY      - N LONG VECTOR CONTAING Y-COORDINATES OF           
!                     VERTICES OF POLYGON.                              
!           N       - NUMBER OF VERTICES IN THE POLYGON.                
!           INOUT   - THE SIGNAL RETURNED:                              
!                     -1 IF THE POINT IS OUTSIDE OF THE POLYGON,        
!                      0 IF THE POINT IS ON AN EDGE OR AT A VERTEX,     
!                      1 IF THE POINT IS INSIDE OF THE POLYGON.         
!                                                                       
!        REMARKS                                                        
!           THE VERTICES MAY BE LISTED CLOCKWISE OR ANTICLOCKWISE.      
!           THE FIRST MAY OPTIONALLY BE REPEATED, IF SO N MAY           
!           OPTIONALLY BE INCREASED BY 1.                               
!           THE INPUT POLYGON MAY BE A COMPOUND POLYGON CONSISTING      
!           OF SEVERAL SEPARATE SUBPOLYGONS. IF SO, THE FIRST VERTEX    
!           OF EACH SUBPOLYGON MUST BE REPEATED, AND WHEN CALCULATING   
!           N, THESE FIRST VERTICES MUST BE COUNTED TWICE.              
!           INOUT IS THE ONLY PARAMETER WHOSE VALUE IS CHANGED.         
!           THE SIZE OF THE ARRAYS MUST BE INCREASED IF N > MAXDIM      
!           WRITTEN BY RANDOLPH FRANKLIN, UNIVERSITY OF OTTAWA, 7/70.   
!                                                                       
!        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  
!           NONE                                                        
!                                                                       
!        METHOD                                                         
!           A VERTICAL LINE IS DRAWN THRU THE POINT IN QUESTION. IF IT  
!           CROSSES THE POLYGON AN ODD NUMBER OF TIMES, THEN THE        
!           POINT IS INSIDE OF THE POLYGON.                             
!                                                                       
!        Modifications:
!          - jcw 12Dec2015 I/O var declarations, 
!                          allow more than 200 maxdim
!     ..................................................................
!                                                                       
!     SUBROUTINE PNPOLY(PX,PY,XX,YY,N,INOUT)
!     REAL X(200),Y(200),XX(N),YY(N)
      real(dbl_kind), intent(in) :: PX, PY, XX(1:N), YY(1:N)
      integer, intent(in) :: N
      integer, intent(inout) :: INOUT
!
      real(dbl_kind) :: X(N),Y(N)
      LOGICAL MX,MY,NX,NY
      INTEGER O, MAXDIM
!     OUTPUT UNIT FOR PRINTED MESSAGES
      DATA O/6/
      MAXDIM=N+1
      IF(N.LE.MAXDIM)GO TO 6
      WRITE(O,7)
7     FORMAT('0WARNING:',I5,' TOO GREAT. RESULTS INVALID')
      RETURN
6     DO 1 I=1,N
        X(I)=XX(I)-PX
1       Y(I)=YY(I)-PY
      INOUT=-1
      DO 2 I=1,N
        J=1+MOD(I,N)
        MX=X(I).GE.0.0
        NX=X(J).GE.0.0
        MY=Y(I).GE.0.0
        NY=Y(J).GE.0.0
        IF(.NOT.((MY.OR.NY).AND.(MX.OR.NX)).OR.(MX.AND.NX)) GO TO 2
        IF(.NOT.(MY.AND.NY.AND.(MX.OR.NX).AND..NOT.(MX.AND.NX))) GO TO 3
        INOUT=-INOUT
        GO TO 2
3       IF((Y(I)*X(J)-X(I)*Y(J))/(X(J)-X(I))) 2,4,5
4       INOUT=0
        RETURN
5       INOUT=-INOUT
2     CONTINUE
      RETURN
      END SUBROUTINE PNPOLY

      end module create_masks 
