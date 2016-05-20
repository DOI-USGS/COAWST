#include "cppdefs.h"
#if defined MCT_INTERP_OC2AT || defined MCT_INTERP_OC2WV || \
    defined MCT_INTERP_WV2AT
      SUBROUTINE get_sparse_matrix (ng, nc_name, num_sparse_elems,      &
     &                              src_dims, dst_dims)
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
      integer, intent(in) :: ng
      integer :: exit_flag

      character (len=70), intent(inout) :: nc_name
      character (len=70) :: grp_name

      integer, intent(out) :: num_sparse_elems
      integer, dimension(2), intent(inout) :: src_dims,                 &
     &                                        dst_dims
!
!  Local variable declarations
!
      integer :: i, j, ndims, status, varid
      integer :: ncid, numelem, ncSPSid, grpid
      integer, dimension(2) :: start, total
      integer, dimension(nf90_max_var_dims) :: dimIDs

      integer :: num_elems_main, num_elems_add
! these 3 are used to read the main normal weigths from scrip.
      integer, dimension(:), pointer :: sparse_rows_main
      integer, dimension(:), pointer :: sparse_cols_main
      real(m8), dimension(:), pointer :: sparse_weights_main
! these 3 are added to read in the weights from, ie, cells
! that wrf thinks are land but roms thinks are ocean
      integer, dimension(:), pointer :: sparse_rows_add
      integer, dimension(:), pointer :: sparse_cols_add
      real(m8), dimension(:), pointer :: sparse_weights_add
!  Then all 6 are added to create the sparse_rows, sparse_cols,
!  and sparse_weights.
!
      num_elems_add=0
!  Open grid NetCDF file for reading.
!
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
      END IF
!
! Get the src and dst grid dimensions.
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
      status=nf90_inq_varid(ncid,'dst_grid_dims', varid)
      IF (status.ne.nf90_noerr) THEN
        WRITE (stdout,30) TRIM('dst_grid_dims')
        exit_flag=4
        ioerror=status
      END IF
      start(1)=1
      total(1)=2
      status=nf90_get_var(ncid, varid, dst_dims, start, total)
!
! Get the destination grid mask
!
      status=nf90_inq_varid(ncid,'dst_grid_imask', varid)
      IF (status.ne.nf90_noerr) THEN
        WRITE (stdout,30) TRIM('dst_grid_imask')
        exit_flag=4
        ioerror=status
      END IF
      numelem=dst_dims(1)*dst_dims(2)
      allocate ( dst_grid_imask(numelem) )
      start(1)=1
      total(1)=numelem
      status=nf90_get_var(ncid, varid, dst_grid_imask, start, total)
!
! Determine number of weights.
!
      status=nf90_inq_varid(ncid,'src_address', varid)
      status=nf90_inquire_variable(ncid,varid,dimids = dimIDs)
      IF (status.ne.nf90_noerr) THEN
        WRITE (stdout,30) TRIM('num_links')
        exit_flag=4
        ioerror=status
      END IF
      status=nf90_inquire_dimension(ncid, dimIDs(1), len=num_elems_main)
      num_sparse_elems=num_elems_main
!
!     These 3 are for the main normal weights.
!
      allocate ( sparse_rows_main(num_elems_main) )
      allocate ( sparse_cols_main(num_elems_main) )
      allocate ( sparse_weights_main(num_elems_main) )
!
! Get the rows.
!
      status=nf90_inq_varid(ncid,'dst_address', varid)
      IF (status.ne.nf90_noerr) THEN
        WRITE (stdout,30) TRIM('dst_address')
        exit_flag=4
        ioerror=status
      END IF
      start(1)=1
      total(1)=num_elems_main
      status=nf90_get_var(ncid, varid, sparse_rows_main, start, total)
!
! Get the cols.
!
      status=nf90_inq_varid(ncid,'src_address', varid)
      IF (status.ne.nf90_noerr) THEN
        WRITE (stdout,30) TRIM('src_address')
        exit_flag=4
        ioerror=status
      END IF
      start(1)=1
      total(1)=num_elems_main
      status=nf90_get_var(ncid, varid, sparse_cols_main, start, total)
!
! Get the weights.
!
      status=nf90_inq_varid(ncid,'remap_matrix', varid)
      IF (status.ne.nf90_noerr) THEN
        WRITE (stdout,30) TRIM('remap_matrix')
        exit_flag=4
        ioerror=status
      END IF
      start(2)=1
      total(2)=num_elems_main
      start(1)=1
      total(1)=1
      status=nf90_get_var(ncid, varid, sparse_weights_main, start, total)
!
! Check to see if the user computed additional weights to account
! for land areas of the src model driving water pts in dst model.
!
      status=nf90_inq_varid(ncid,'add_src_address', varid)
      IF (status.eq.nf90_noerr) THEN
        status=nf90_inquire_variable(ncid,varid,dimids = dimIDs)
        status=nf90_inquire_dimension(ncid, dimIDs(1), len=num_elems_add)
        num_sparse_elems=num_sparse_elems+num_elems_add
!
!       These 3 are the added arrays.
!
        allocate ( sparse_rows_add(num_elems_add) )
        allocate ( sparse_cols_add(num_elems_add) )
        allocate ( sparse_weights_add(num_elems_add) )
!
! Get the rows.
!
        status=nf90_inq_varid(ncid,'add_dst_address', varid)
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,30) TRIM('add_dst_address')
          exit_flag=4
          ioerror=status
        END IF
        start(1)=1
        total(1)=num_elems_add
        status=nf90_get_var(ncid, varid, sparse_rows_add, start, total)
!
! Get the cols.
!
        status=nf90_inq_varid(ncid,'add_src_address', varid)
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,30) TRIM('add_src_address')
          exit_flag=4
          ioerror=status
        END IF
        start(1)=1
        total(1)=num_elems_add
        status=nf90_get_var(ncid, varid, sparse_cols_add, start, total)
!
! Get the weights.
!
        status=nf90_inq_varid(ncid,'add_remap_matrix', varid)
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,30) TRIM('add_remap_matrix')
          exit_flag=4
          ioerror=status
        END IF
        start(1)=1
        total(1)=num_elems_add
        status=nf90_get_var(ncid, varid, sparse_weights_add, start, total)
      END IF
!
! Allocate total arrays.
!
!     These 3 are the total arrays.
      allocate ( sparse_rows(num_sparse_elems) )
      allocate ( sparse_cols(num_sparse_elems) )
      allocate ( sparse_weights(num_sparse_elems) )
!
!  Copy the main and add parts to the total.
!
      DO i=1,num_elems_main
        sparse_rows(i)=sparse_rows_main(i)
        sparse_cols(i)=sparse_cols_main(i)
        sparse_weights(i)=sparse_weights_main(i)
      END DO
      IF (num_elems_add.gt.0) THEN
        DO i=1,num_elems_add
          j=num_elems_main+i
          sparse_rows(j)=sparse_rows_add(i)
          sparse_cols(j)=sparse_cols_add(i)
          sparse_weights(j)=sparse_weights_add(i)
        END DO
        deallocate (sparse_rows_add)
        deallocate (sparse_cols_add)
        deallocate (sparse_weights_add)
      END IF
!
      deallocate (sparse_rows_main)
      deallocate (sparse_cols_main)
      deallocate (sparse_weights_main)
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
