#include "cppdefs.h"
#if defined MCT_INTERP_OC2AT || defined MCT_INTERP_OC2WV
      SUBROUTINE get_sparse_matrix (ng, nc_name, num_sparse_elems,      &
     &                              src_grid_dims, dst_grid_dims)
!
!svn $Id: get_sparse_matrix.F 1336 2008-01-24 02:45:56Z jcwarner $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine reads the sparse matrix weights.                       !
!                                                                      !
!=======================================================================
!
      USE mod_iounits
      USE mod_netcdf
      USE mod_scalars
      USE mod_ncparam
      USE mod_coupler
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng

      character (len=70), intent(in) :: nc_name

      integer, intent(out) :: num_sparse_elems
      integer, dimension(2), intent(inout) :: src_grid_dims,            &
     &                                        dst_grid_dims
!
!  Local variable declarations
!
      integer ::  ndims, status, varid, ncid, numelem, ncSPSid
      integer, dimension(2) :: start, total
      integer, dimension(nf90_max_var_dims) :: dimIDs

!  Open grid NetCDF file for reading.
!
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
      ncid=ncSPSid
!
! Determine number of weights.
!
      status=nf90_inq_varid(ncSPSid,'src_address', varid)
      status=nf90_inquire_variable(ncSPSid,varid,dimids = dimIDs)
      IF (status.ne.nf90_noerr) THEN
        WRITE (stdout,30) TRIM('num_links')
        exit_flag=4
        ioerror=status
      END IF
      status=nf90_inquire_dimension(ncid, dimIDs(1), len=src_grid_dims(1))
      num_sparse_elems=src_grid_dims(1)
!
! Allocate arrays.
!
!      IF (.not.allocated(sparse_rows)) THEN
        allocate ( sparse_rows(num_sparse_elems) )
!      END IF
!      IF (.not.allocated(sparse_cols)) THEN
        allocate ( sparse_cols(num_sparse_elems) )
!      END IF
!      IF (.not.allocated(sparse_weights)) THEN
        allocate ( sparse_weights(num_sparse_elems) )
!      END IF
!
! Get the src and dest grid dimensions.
!
      status=nf90_inq_varid(ncSPSid,'src_grid_dims', varid)
      IF (status.ne.nf90_noerr) THEN
        WRITE (stdout,30) TRIM('src_grid_dims')
        exit_flag=4
        ioerror=status
      END IF
      start(1)=1
      total(1)=2
      status=nf90_get_var(ncid, varid, src_grid_dims, start, total)
!
      status=nf90_inq_varid(ncSPSid,'dst_grid_dims', varid)
      IF (status.ne.nf90_noerr) THEN
        WRITE (stdout,30) TRIM('dst_grid_dims')
        exit_flag=4
        ioerror=status
      END IF
      start(1)=1
      total(1)=2
      status=nf90_get_var(ncid, varid, dst_grid_dims, start, total)
!
! Get the rows.
!
      status=nf90_inq_varid(ncSPSid,'dst_address', varid)
      IF (status.ne.nf90_noerr) THEN
        WRITE (stdout,30) TRIM('dst_address')
        exit_flag=4
        ioerror=status
      END IF
      start(1)=1
      total(1)=num_sparse_elems
      status=nf90_get_var(ncid, varid, sparse_rows, start, total)
!
! Get the cols.
!
      status=nf90_inq_varid(ncSPSid,'src_address', varid)
      IF (status.ne.nf90_noerr) THEN
        WRITE (stdout,30) TRIM('src_address')
        exit_flag=4
        ioerror=status
      END IF
      start(1)=1
      total(1)=num_sparse_elems
      status=nf90_get_var(ncid, varid, sparse_cols, start, total)
!
! Get the weights.
!
      status=nf90_inq_varid(ncSPSid,'remap_matrix', varid)
      IF (status.ne.nf90_noerr) THEN
        WRITE (stdout,30) TRIM('remap_matrix')
        exit_flag=4
        ioerror=status
      END IF
      start(2)=1
      total(2)=num_sparse_elems
      start(1)=1
      total(1)=1
      status=nf90_get_var(ncid, varid, sparse_weights, start, total)
!
! Get the destination grid mask
!
      status=nf90_inq_varid(ncSPSid,'dst_grid_imask', varid)
      IF (status.ne.nf90_noerr) THEN
        WRITE (stdout,30) TRIM('dst_grid_imask')
        exit_flag=4
        ioerror=status
      END IF
      numelem=dst_grid_dims(1)*dst_grid_dims(2)
      allocate ( dst_grid_imask(numelem) )
      start(1)=1
      total(1)=numelem
      status=nf90_get_var(ncid, varid, dst_grid_imask, start, total)
!
! Close GRID NetCDF file.
!
      status=nf90_close(ncSPSid)
      ncSPSid=-1
!

  5   FORMAT (/,' GET_SPARSE - error while opening file: ', a)
 10   FORMAT (/,' GET_SPARSE - error while reading attribute: ', a, &
     &          ' for variable: ', a)
 20   FORMAT (/,' GET_SPARSE - error while inquiring attribute: ',  &
     &        a,' for variable: ', a)
 30   FORMAT (/,' GET_SPARSE - cannot inquire ID for variable: ', a)
 40   FORMAT (/,' GET_SPARSE - error while inquiring dimensions',   &
     &          ' for variable: ', a)
 50   FORMAT (/,' GET_SPARSE - error while reading variable: ', a)
      RETURN
      END SUBROUTINE get_sparse_matrix
#else
      SUBROUTINE get_sparse_matrix (ng, tile)
      RETURN
      END SUBROUTINE get_sparse_matrix
#endif
