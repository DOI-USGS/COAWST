#include "cppdefs.h"
#if defined MCT_INTERP_OC2AT || defined MCT_INTERP_WV2AT
      SUBROUTINE get_wrf_moving_grids
!
!svn $Id: get_wrf_moving_grids.F 1336 2008-01-24 02:45:56Z jcwarner $
!================================================== John C. Warner =====
!  Copyright (c) 2002-2020 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine reads the sparse matrix wrf grid sizes to determine    !
!  if the user has these as moving grids.                              !
!                                                                      !
!=======================================================================
!
      USE mod_coupler_iounits
      USE netcdf

#ifdef ROMS_MODEL
      USE mod_iounits
#endif
!
      implicit none
!
!  Imported variable declarations.
!
!
!  Local variable declarations
!
      integer :: ia, status, varid
      integer :: ncid, ncSPSid, grpid
      integer :: exit_flag

      character (len=70) :: nc_name
      character (len=70) :: grp_name

      integer, dimension(2) :: src_dims
      integer, dimension(2) :: start, total

  35  FORMAT(a3,i1,a7,i1,a11)
!      IF (Myrank.eq.MyMaster) THEN
!
!  Set moving nest of all grids to be zero.
        DO ia=1,Natm_grids
          moving_nest(ia)=0
        END DO
!
!  Now do the other wrf grids.
!
        DO ia=2,Natm_grids
!
!  Open grid NetCDF file for reading.
!
# ifdef MCT_INTERP_OC2AT
          IF (scrip_opt.eq.1) THEN
            write(nc_name,35) 'atm',ia,'_to_ocn',1,'_weights.nc'
          ELSE
            nc_name=A2Oname(ia,1)
          END IF
# else
          IF (scrip_opt.eq.1) THEN
            write(nc_name,35) 'atm',ia,'_to_wav',1,'_weights.nc'
          ELSE
            nc_name=A2Wname(ia,1)
          END IF
# endif
          IF (scrip_opt.eq.1) THEN
            grp_name=TRIM(ADJUSTL(nc_name))
            nc_name=SCRIPname
          END IF
          ncSPSid=-1
          IF (ncSPSid.eq.-1) THEN
            status=nf90_open(TRIM(nc_name), nf90_nowrite, ncSPSid)
            IF (status.ne.nf90_noerr) THEN
              WRITE (stdout,5) TRIM(nc_name)
              exit_flag=2
              ioerror=status
              RETURN
            END IF
          END IF
!
! Get the group id.
!
          IF (scrip_opt.eq.1) THEN
            status=nf90_inq_ncid(ncSPSid,grp_name,grpid)
            ncid=grpid
          ELSE
            ncid=ncSPSid
          END IF
!
! Get the dst grid dimensions.
!
!
          status=nf90_inq_varid(ncid,'src_grid_dims', varid)
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,30) TRIM('src_grid_dims')
            exit_flag=4
            ioerror=status
          END IF
          start(1)=1
          total(1)=2
          status=nf90_get_var(ncid, varid, src_dims, start, total)
!
          IF ((wrf_e_we(ia)-1).ne.src_dims(1)) THEN
            moving_nest(ia)=1
          END IF
!
! Close GRID NetCDF file.
!
          status=nf90_close(ncSPSid)
          ncSPSid=-1
!
        END DO
!      END IF
!
!      CALL mpi_bcast(moving_nest, Natm_grids, MPI_INTEGER,              &
!     &               MyMaster, MPI_COMM_WORLD, MyError)
!
!
  5   FORMAT (/,' GET_WRF_MOV - error while opening file: ', a)
 10   FORMAT (/,' GET_WRF_MOV - error while reading attribute: ', a, &
     &          ' for variable: ', a)
 20   FORMAT (/,' GET_WRF_MOV - error while inquiring attribute: ',  &
     &        a,' for variable: ', a)
 30   FORMAT (/,' GET_WRF_MOV - cannot inquire ID for variable: ', a)
      RETURN
      END SUBROUTINE get_wrf_moving_grids
#else
      SUBROUTINE get_wrf_moving_grids
      RETURN
      END SUBROUTINE get_wrf_moving_grids
#endif
