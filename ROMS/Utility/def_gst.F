#include "cppdefs.h"
      MODULE def_gst_mod
#if defined PROPAGATOR && defined CHECKPOINTING
!
!git $Id$
!svn $Id: def_gst.F 1151 2023-02-09 03:08:53Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2023 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This module creates the checkpointing restart file for the GST      !
!  (Generalized Stability Theory) analysis drivers using either the    !
!  standard NetCDF library or the Parallel-IO (PIO) library.           !
!                                                                      !
!  The NetCDF file contains all the necessary eigenproblem fields      !
!  for restating ARPACK.                                               !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_ncparam
      USE mod_scalars
      USE mod_storage
!
      USE def_dim_mod,  ONLY : def_dim
      USE def_info_mod, ONLY : def_info
      USE def_var_mod,  ONLY : def_var
      USE strings_mod,  ONLY : FoundError
!
      implicit none
!
      PUBLIC  :: def_gst
      PRIVATE :: def_gst_nf90
# if defined PIO_LIB && defined DISTRIBUTE
      PRIVATE :: def_gst_pio
# endif
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE def_gst (ng, model)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
!
!  Local variable declarations.
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__
!
!-----------------------------------------------------------------------
!  Create a new history file according to IO type.
!-----------------------------------------------------------------------
!
      SELECT CASE (GST(ng)%IOtype)
        CASE (io_nf90)
          CALL def_gst_nf90 (ng, model)

# if defined PIO_LIB && defined DISTRIBUTE
        CASE (io_pio)
          CALL def_gst_pio (ng, model)
# endif
        CASE DEFAULT
          IF (Master) THEN
            WRITE (stdout,10) GST(ng)%IOtype
          END IF
          exit_flag=3
      END SELECT
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
  10  FORMAT (' DEF_GST - Illegal output file type, io_type = ',i0,     &
     &        /,15x,'Check KeyWord ''OUT_LIB'' in ''roms.in''.')
!
      RETURN
      END SUBROUTINE def_gst
!
!***********************************************************************
      SUBROUTINE def_gst_nf90 (ng, model)
!***********************************************************************
!
      USE mod_netcdf
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
!
!  Local variable declarations.
!
      integer, parameter :: Natt = 25
!
      integer :: Mstatedim, NCVdim, SworkDdim, SworkLdim
      integer :: char1dim, char2dim, iaitrdim, iaupddim, iaup2dim
      integer :: iparamdim, ipntrdim, laitrdim, laup2dim, raitrdim
      integer :: raup2dim
      integer :: i, j, status, varid
      integer :: DimIDs(nDimID), vardim(2)
!
      real(r8) :: Aval(6)
!
      character (len=256)    :: ncname
      character (len=MaxLen) :: Vinfo(Natt)

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", def_gst_nf90"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Set and report file name.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ncname=GST(ng)%name
!
      IF (Master) THEN
        IF (.not.LrstGST) THEN
          WRITE (stdout,10) ng, TRIM(ncname)
        ELSE
          WRITE (stdout,20) ng, TRIM(ncname)
        END IF
      END IF
!
!=======================================================================
!  Create a new averages NetCDF file.
!=======================================================================
!
      DEFINE : IF (.not.LrstGST) THEN
        CALL netcdf_create (ng, model, TRIM(ncname), GST(ng)%ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) THEN
          IF (Master) WRITE (stdout,30) TRIM(ncname)
          RETURN
        END IF
!
!-----------------------------------------------------------------------
!  Define the dimensions of staggered fields.
!-----------------------------------------------------------------------
!
        status=def_dim(ng, model, GST(ng)%ncid, ncname, 'Mstate',       &
     &                 Mstate(ng), Mstatedim)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, GST(ng)%ncid, ncname, 'NCV',          &
     &                 NCV, NCVdim)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, GST(ng)%ncid, ncname, 'LworkD',       &
     &                 3*Mstate(ng), SworkDdim)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, GST(ng)%ncid, ncname, 'LworkL',       &
     &                 LworkL, SworkLdim)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, GST(ng)%ncid, ncname, 'iparam',       &
     &                 SIZE(iparam), iparamdim)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, GST(ng)%ncid, ncname, 'ipntr',        &
     &                 SIZE(ipntr), ipntrdim)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, GST(ng)%ncid, ncname, 'iaupd',        &
     &                 SIZE(iaupd), iaupddim)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, GST(ng)%ncid, ncname, 'laitr',        &
     &                 SIZE(laitr), laitrdim)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, GST(ng)%ncid, ncname, 'iaitr',        &
     &                 SIZE(iaitr), iaitrdim)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, GST(ng)%ncid, ncname, 'raitr',        &
     &                 SIZE(raitr), raitrdim)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, GST(ng)%ncid, ncname, 'laup2',        &
     &                 SIZE(laup2), laup2dim)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, GST(ng)%ncid, ncname, 'iaup2',        &
     &                 SIZE(iaup2), iaup2dim)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, GST(ng)%ncid, ncname, 'raup2',        &
     &                 SIZE(raup2), raup2dim)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, GST(ng)%ncid, ncname, 'char2',        &
     &                 2, char2dim)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!-----------------------------------------------------------------------
!  Define global attributes.
!-----------------------------------------------------------------------
!
        CALL def_info (ng, model, GST(ng)%ncid, ncname, DimIDs)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!-----------------------------------------------------------------------
!  Define variables and their attributes.
!-----------------------------------------------------------------------
!
!  Initialize local information variable arrays.
!
        DO i=1,Natt
          DO j=1,LEN(Vinfo(1))
            Vinfo(i)(j:j)=' '
          END DO
        END DO
        DO i=1,6
          Aval(i)=0.0_r8
        END DO
!
!  Define number of eigenvalues to compute.
!
        Vinfo( 1)='NEV'
        Vinfo( 2)='number of eigenvalues to compute'
        status=def_var(ng, model, GST(ng)%ncid, varid, nf90_int,        &
     &                 1, (/0/), Aval, Vinfo, ncname,                   &
     &                 SetParAccess = .FALSE.)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Define number of Lanczos vectors to compute.
!
        Vinfo( 1)='NCV'
        Vinfo( 2)='number of Lanczos vectors to compute'
        status=def_var(ng, model, GST(ng)%ncid, varid, nf90_int,        &
     &                 1, (/0/), Aval, Vinfo, ncname,                   &
     &                 SetParAccess = .FALSE.)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Define size of eigenvalue problem.
!
        Vinfo( 1)='Mstate'
        Vinfo( 2)='total size of eigenvalue problem'
        status=def_var(ng, model, GST(ng)%ncid, varid, nf90_int,        &
     &                 1, (/0/), Aval, Vinfo, ncname,                   &
     &                 SetParAccess = .FALSE.)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Define iteration number.
!
        Vinfo( 1)='iter'
        Vinfo( 2)='iteration number'
        status=def_var(ng, model, GST(ng)%ncid, varid, nf90_int,        &
     &                 1, (/0/), Aval, Vinfo, ncname,                   &
     &                 SetParAccess = .FALSE.)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Define reverse communications flag.
!
        Vinfo( 1)='ido'
        Vinfo( 2)='reverse communications flag'
        status=def_var(ng, model, GST(ng)%ncid, varid, nf90_int,        &
     &                 1, (/0/), Aval, Vinfo, ncname,                   &
     &                 SetParAccess = .FALSE.)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Define information and error flag.
!
        Vinfo( 1)='info'
        Vinfo( 2)='information and error flag'
        status=def_var(ng, model, GST(ng)%ncid, varid, nf90_int,        &
     &                 1, (/0/), Aval, Vinfo, ncname,                   &
     &                 SetParAccess = .FALSE.)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Define eigenvalue problem type.
!
        Vinfo( 1)='bmat'
        Vinfo( 2)='eigenvalue problem type'
        status=def_var(ng, model, GST(ng)%ncid, varid, nf90_char,       &
     &                 1, (/0/), Aval, Vinfo, ncname,                   &
     &                 SetParAccess = .FALSE.)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Define Ritz eigenvalues to compute.
!
        Vinfo( 1)='which'
        Vinfo( 2)='Ritz eigenvalues to compute'
        status=def_var(ng, model, GST(ng)%ncid, varid, nf90_char,       &
     &                 1, (/char2dim/), Aval, Vinfo, ncname,            &
     &                 SetParAccess = .FALSE.)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Define form of basis function.
!
        Vinfo( 1)='howmany'
        Vinfo( 2)='form of basis function'
        status=def_var(ng, model, GST(ng)%ncid, varid, nf90_char,       &
     &                 1, (/0/), Aval, Vinfo, ncname,                   &
     &                 SetParAccess = .FALSE.)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Define relative accuracy of computed Ritz values.
!
        Vinfo( 1)='Ritz_tol'
        Vinfo( 2)='relative accuracy of computed Ritz values'
        status=def_var(ng, model, GST(ng)%ncid, varid, NF_FRST,         &
     &                 1, (/0/), Aval, Vinfo, ncname,                   &
     &                 SetParAccess = .FALSE.)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Define eigenproblem parameters.
!
        Vinfo( 1)='iparam'
        Vinfo( 2)='eigenproblem parameters'
        status=def_var(ng, model, GST(ng)%ncid, varid, nf90_int,        &
     &                 1, (/iparamdim/), Aval, Vinfo, ncname,           &
     &                 SetParAccess = .FALSE.)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Define pointers to mark starting location in work arrays.
!
        Vinfo( 1)='ipntr'
        Vinfo( 2)='pointers to mark starting location in work arrays'
        status=def_var(ng, model, GST(ng)%ncid,varid, nf90_int,         &
     &                 1, (/ipntrdim/), Aval, Vinfo, ncname,            &
     &                 SetParAccess = .FALSE.)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Define ARPACK internal integer parameters to _aupd routines.
!
        Vinfo( 1)='iaupd'
        Vinfo( 2)='ARPACK internal integer parameters to _aupd routines'
        status=def_var(ng, model, GST(ng)%ncid, varid, nf90_int,        &
     &                 1, (/iaupddim/), Aval, Vinfo, ncname,            &
     &                 SetParAccess = .FALSE.)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Define ARPACK internal integer parameters to _aitr routines.
!
        Vinfo( 1)='iaitr'
        Vinfo( 2)='ARPACK internal integer parameters to _aitr routines'
        status=def_var(ng, model, GST(ng)%ncid, varid, nf90_int,        &
     &                 1, (/iaitrdim/), Aval, Vinfo, ncname,            &
     &                 SetParAccess = .FALSE.)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Define ARPACK internal integer parameters to _aup2 routines.
!
        Vinfo( 1)='iaup2'
        Vinfo( 2)='ARPACK internal integer parameters to _aup2 routines'
        status=def_var(ng, model, GST(ng)%ncid, varid, nf90_int,        &
     &                 1, (/iaup2dim/), Aval, Vinfo, ncname,            &
     &                 SetParAccess = .FALSE.)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Define ARPACK internal logical parameters to _aitr routines.
!
        Vinfo( 1)='laitr'
        Vinfo( 2)='ARPACK internal logical parameters to _aitr routines'
        Vinfo( 9)='.FALSE.'
        Vinfo(10)='.TRUE.'
        status=def_var(ng, model, GST(ng)%ncid, varid, nf90_int,        &
     &                 1, (/laitrdim/), Aval, Vinfo, ncname,            &
     &                 SetParAccess = .FALSE.)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Define ARPACK internal logical parameters to _aup2 routines.
!
        Vinfo( 1)='laup2'
        Vinfo( 2)='ARPACK internal logical parameters to _aup2 routines'
        Vinfo( 9)='.FALSE.'
        Vinfo(10)='.TRUE.'
        status=def_var(ng, model, GST(ng)%ncid, varid, nf90_int,        &
     &                 1, (/laup2dim/), Aval, Vinfo, ncname,            &
     &                 SetParAccess = .FALSE.)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Define ARPACK internal real parameters to _aitr routines.
!
        Vinfo( 1)='raitr'
        Vinfo( 2)='ARPACK internal real parameters to _aitr routines'
        status=def_var(ng, model, GST(ng)%ncid, varid, NF_FRST,         &
     &                 1, (/raitrdim/), Aval, Vinfo, ncname,            &
     &                 SetParAccess = .FALSE.)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Define ARPACK internal real parameters to _aup2 routines.
!
        Vinfo( 1)='raup2'
        Vinfo( 2)='ARPACK internal real parameters to _aup2 routines'
        status=def_var(ng, model, GST(ng)%ncid, varid, NF_FRST,         &
     &                 1, (/raup2dim/), Aval, Vinfo, ncname,            &
     &                 SetParAccess = .FALSE.)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Define Lanczos/Arnoldi basis vectors.
!
        Vinfo( 1)='Bvec'
        Vinfo( 2)='Lanczos/Arnoldi basis vectors'
        vardim(1)=Mstatedim
        vardim(2)=NCVdim
        status=def_var(ng, model, GST(ng)%ncid, varid, NF_FRST,         &
     &                 2, vardim, Aval, Vinfo, ncname)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Define eigenproblem residual vector.
!
        Vinfo( 1)='resid'
        Vinfo( 2)='eigenproblem residual vector'
        status=def_var(ng, model, GST(ng)%ncid, varid, NF_FRST,         &
     &                 1, (/Mstatedim/), Aval, Vinfo, ncname)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Define state reverse communications work array.
!
        Vinfo( 1)='SworkD'
        Vinfo( 2)='reverse communications state array'
        status=def_var(ng, model, GST(ng)%ncid, varid, NF_FRST,         &
     &                 1, (/SworkDdim/), Aval, Vinfo, ncname)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Define eigenproblem work array.
!
        Vinfo( 1)='SworkL'
        Vinfo( 2)='eigenproblem work array'
        status=def_var(ng, model, GST(ng)%ncid, varid, NF_FRST,         &
     &                 1, (/SworkLdim/), Aval, Vinfo, ncname)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!-----------------------------------------------------------------------
!  Leave definition mode.
!-----------------------------------------------------------------------
!
        CALL netcdf_enddef (ng, model, ncname, GST(ng)%ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      END IF DEFINE
!
  10  FORMAT (/,2x,'DEF_GST_NF90     - creating checkpointing file,',   &
     &        t56,'Grid ',i2.2,': ',a)
  20  FORMAT (/,2x,'DEF_GST_NF90     - inquiring checkpointing file,',  &
     &        t56,'Grid ',i2.2,': ',a)
  30  FORMAT (/,' DEF_GST_NF90 - unable to create checkpointing',       &
     &        ' NetCDF file: ',a)
!
      RETURN
      END SUBROUTINE def_gst_nf90

# if defined PIO_LIB && defined DISTRIBUTE
!
!***********************************************************************
      SUBROUTINE def_gst_pio (ng, model)
!***********************************************************************
!
      USE mod_pio_netcdf
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
!
!  Local variable declarations.
!
      integer, parameter :: Natt = 25
!
      integer :: Mstatedim, NCVdim, SworkDdim, SworkLdim
      integer :: char1dim, char2dim, iaitrdim, iaupddim, iaup2dim
      integer :: iparamdim, ipntrdim, laitrdim, laup2dim, raitrdim
      integer :: raup2dim
      integer :: i, j, status
      integer :: DimIDs(nDimID), vardim(2)
!
      real(r8) :: Aval(6)
!
      character (len=256)    :: ncname
      character (len=MaxLen) :: Vinfo(Natt)

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", def_gst_pio"
!
      TYPE (Var_desc_t)  :: VarDesc
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Set and report file name.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ncname=GST(ng)%name
!
      IF (Master) THEN
        IF (.not.LrstGST) THEN
          WRITE (stdout,10) ng, TRIM(ncname)
        ELSE
          WRITE (stdout,20) ng, TRIM(ncname)
        END IF
      END IF
!
!=======================================================================
!  Create a new averages NetCDF file.
!=======================================================================
!
      DEFINE : IF (.not.LrstGST) THEN
        CALL pio_netcdf_create (ng, model, TRIM(ncname),                &
     &                          GST(ng)%pioFile)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) THEN
          IF (Master) WRITE (stdout,30) TRIM(ncname)
          RETURN
        END IF
!
!-----------------------------------------------------------------------
!  Define the dimensions of staggered fields.
!-----------------------------------------------------------------------
!
        status=def_dim(ng, model, GST(ng)%pioFile, ncname, 'Mstate',    &
     &                 Mstate(ng), Mstatedim)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, GST(ng)%pioFile, ncname, 'NCV',       &
     &                 NCV, NCVdim)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, GST(ng)%pioFile, ncname, 'LworkD',    &
     &                 3*Mstate(ng), SworkDdim)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, GST(ng)%pioFile, ncname, 'LworkL',    &
     &                 LworkL, SworkLdim)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, GST(ng)%pioFile, ncname, 'iparam',    &
     &                 SIZE(iparam), iparamdim)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, GST(ng)%pioFile, ncname, 'ipntr',     &
     &                 SIZE(ipntr), ipntrdim)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, GST(ng)%pioFile, ncname, 'iaupd',     &
     &                 SIZE(iaupd), iaupddim)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, GST(ng)%pioFile, ncname, 'laitr',     &
     &                 SIZE(laitr), laitrdim)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, GST(ng)%pioFile, ncname, 'iaitr',     &
     &                 SIZE(iaitr), iaitrdim)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, GST(ng)%pioFile, ncname, 'raitr',     &
     &                 SIZE(raitr), raitrdim)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, GST(ng)%pioFile, ncname, 'laup2',     &
     &                 SIZE(laup2), laup2dim)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, GST(ng)%pioFile, ncname, 'iaup2',     &
     &                 SIZE(iaup2), iaup2dim)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, GST(ng)%pioFile, ncname, 'raup2',     &
     &                 SIZE(raup2), raup2dim)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, GST(ng)%pioFile, ncname, 'char2',     &
     &                 2, char2dim)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!-----------------------------------------------------------------------
!  Define global attributes.
!-----------------------------------------------------------------------
!
        CALL def_info (ng, model, GST(ng)%pioFile, ncname, DimIDs)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!-----------------------------------------------------------------------
!  Define variables and their attributes.
!-----------------------------------------------------------------------
!
!  Initialize local information variable arrays.
!
        DO i=1,Natt
          DO j=1,LEN(Vinfo(1))
            Vinfo(i)(j:j)=' '
          END DO
        END DO
        DO i=1,6
          Aval(i)=0.0_r8
        END DO
!
!  Define number of eigenvalues to compute.
!
        Vinfo( 1)='NEV'
        Vinfo( 2)='number of eigenvalues to compute'
        status=def_var(ng, model, GST(ng)%pioFile, varDesc, PIO_int,    &
     &                 1, (/0/), Aval, Vinfo, ncname,                   &
     &                 SetParAccess = .FALSE.)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Define number of Lanczos vectors to compute.
!
        Vinfo( 1)='NCV'
        Vinfo( 2)='number of Lanczos vectors to compute'
        status=def_var(ng, model, GST(ng)%pioFile, varDesc, PIO_int,    &
     &                 1, (/0/), Aval, Vinfo, ncname,                   &
     &                 SetParAccess = .FALSE.)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Define size of eigenvalue problem.
!
        Vinfo( 1)='Mstate'
        Vinfo( 2)='total size of eigenvalue problem'
        status=def_var(ng, model, GST(ng)%pioFile, varDesc, PIO_int,    &
     &                 1, (/0/), Aval, Vinfo, ncname,                   &
     &                 SetParAccess = .FALSE.)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Define iteration number.
!
        Vinfo( 1)='iter'
        Vinfo( 2)='iteration number'
        status=def_var(ng, model, GST(ng)%pioFile, varDesc, PIO_int,    &
     &                 1, (/0/), Aval, Vinfo, ncname,                   &
     &                 SetParAccess = .FALSE.)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Define reverse communications flag.
!
        Vinfo( 1)='ido'
        Vinfo( 2)='reverse communications flag'
        status=def_var(ng, model, GST(ng)%pioFile, varDesc, PIO_int,    &
     &                 1, (/0/), Aval, Vinfo, ncname,                   &
     &                 SetParAccess = .FALSE.)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Define information and error flag.
!
        Vinfo( 1)='info'
        Vinfo( 2)='information and error flag'
        status=def_var(ng, model, GST(ng)%pioFile, varDesc, PIO_int,    &
     &                 1, (/0/), Aval, Vinfo, ncname,                   &
     &                 SetParAccess = .FALSE.)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Define eigenvalue problem type.
!
        Vinfo( 1)='bmat'
        Vinfo( 2)='eigenvalue problem type'
        status=def_var(ng, model, GST(ng)%pioFile, varDesc, PIO_char,   &
     &                 1, (/0/), Aval, Vinfo, ncname,                   &
     &                 SetParAccess = .FALSE.)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Define Ritz eigenvalues to compute.
!
        Vinfo( 1)='which'
        Vinfo( 2)='Ritz eigenvalues to compute'
        status=def_var(ng, model, GST(ng)%pioFile, varDesc, PIO_char,   &
     &                 1, (/char2dim/), Aval, Vinfo, ncname,            &
     &                 SetParAccess = .FALSE.)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Define form of basis function.
!
        Vinfo( 1)='howmany'
        Vinfo( 2)='form of basis function'
        status=def_var(ng, model, GST(ng)%pioFile, varDesc, PIO_char,   &
     &                 1, (/0/), Aval, Vinfo, ncname,                   &
     &                 SetParAccess = .FALSE.)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Define relative accuracy of computed Ritz values.
!
        Vinfo( 1)='Ritz_tol'
        Vinfo( 2)='relative accuracy of computed Ritz values'
        status=def_var(ng, model, GST(ng)%pioFile, varDesc, PIO_FRST,   &
     &                 1, (/0/), Aval, Vinfo, ncname,                   &
     &                 SetParAccess = .FALSE.)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Define eigenproblem parameters.
!
        Vinfo( 1)='iparam'
        Vinfo( 2)='eigenproblem parameters'
        status=def_var(ng, model, GST(ng)%pioFile, varDesc, PIO_int,    &
     &                 1, (/iparamdim/), Aval, Vinfo, ncname,           &
     &                 SetParAccess = .FALSE.)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Define pointers to mark starting location in work arrays.
!
        Vinfo( 1)='ipntr'
        Vinfo( 2)='pointers to mark starting location in work arrays'
        status=def_var(ng, model, GST(ng)%pioFile,varDesc, PIO_int,     &
     &                 1, (/ipntrdim/), Aval, Vinfo, ncname,            &
     &                 SetParAccess = .FALSE.)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Define ARPACK internal integer parameters to _aupd routines.
!
        Vinfo( 1)='iaupd'
        Vinfo( 2)='ARPACK internal integer parameters to _aupd routines'
        status=def_var(ng, model, GST(ng)%pioFile, varDesc, PIO_int,    &
     &                 1, (/iaupddim/), Aval, Vinfo, ncname,            &
     &                 SetParAccess = .FALSE.)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Define ARPACK internal integer parameters to _aitr routines.
!
        Vinfo( 1)='iaitr'
        Vinfo( 2)='ARPACK internal integer parameters to _aitr routines'
        status=def_var(ng, model, GST(ng)%pioFile, varDesc, PIO_int,    &
     &                 1, (/iaitrdim/), Aval, Vinfo, ncname,            &
     &                 SetParAccess = .FALSE.)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Define ARPACK internal integer parameters to _aup2 routines.
!
        Vinfo( 1)='iaup2'
        Vinfo( 2)='ARPACK internal integer parameters to _aup2 routines'
        status=def_var(ng, model, GST(ng)%pioFile, varDesc, PIO_int,    &
     &                 1, (/iaup2dim/), Aval, Vinfo, ncname,            &
     &                 SetParAccess = .FALSE.)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Define ARPACK internal logical parameters to _aitr routines.
!
        Vinfo( 1)='laitr'
        Vinfo( 2)='ARPACK internal logical parameters to _aitr routines'
        Vinfo( 9)='.FALSE.'
        Vinfo(10)='.TRUE.'
        status=def_var(ng, model, GST(ng)%pioFile, varDesc, PIO_int,    &
     &                 1, (/laitrdim/), Aval, Vinfo, ncname,            &
     &                 SetParAccess = .FALSE.)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Define ARPACK internal logical parameters to _aup2 routines.
!
        Vinfo( 1)='laup2'
        Vinfo( 2)='ARPACK internal logical parameters to _aup2 routines'
        Vinfo( 9)='.FALSE.'
        Vinfo(10)='.TRUE.'
        status=def_var(ng, model, GST(ng)%pioFile, varDesc, PIO_int,    &
     &                 1, (/laup2dim/), Aval, Vinfo, ncname,            &
     &                 SetParAccess = .FALSE.)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Define ARPACK internal real parameters to _aitr routines.
!
        Vinfo( 1)='raitr'
        Vinfo( 2)='ARPACK internal real parameters to _aitr routines'
        status=def_var(ng, model, GST(ng)%pioFile, varDesc, PIO_FRST,   &
     &                 1, (/raitrdim/), Aval, Vinfo, ncname,            &
     &                 SetParAccess = .FALSE.)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Define ARPACK internal real parameters to _aup2 routines.
!
        Vinfo( 1)='raup2'
        Vinfo( 2)='ARPACK internal real parameters to _aup2 routines'
        status=def_var(ng, model, GST(ng)%pioFile, varDesc, PIO_FRST,   &
     &                 1, (/raup2dim/), Aval, Vinfo, ncname,            &
     &                 SetParAccess = .FALSE.)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Define Lanczos/Arnoldi basis vectors.
!
        Vinfo( 1)='Bvec'
        Vinfo( 2)='Lanczos/Arnoldi basis vectors'
        vardim(1)=Mstatedim
        vardim(2)=NCVdim
        status=def_var(ng, model, GST(ng)%pioFile, varDesc, PIO_FRST,   &
     &                 2, vardim, Aval, Vinfo, ncname)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Define eigenproblem residual vector.
!
        Vinfo( 1)='resid'
        Vinfo( 2)='eigenproblem residual vector'
        status=def_var(ng, model, GST(ng)%pioFile, varDesc, PIO_FRST,   &
     &                 1, (/Mstatedim/), Aval, Vinfo, ncname)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Define state reverse communications work array.
!
        Vinfo( 1)='SworkD'
        Vinfo( 2)='reverse communications state array'
        status=def_var(ng, model, GST(ng)%pioFile, varDesc, PIO_FRST,   &
     &                 1, (/SworkDdim/), Aval, Vinfo, ncname)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Define eigenproblem work array.
!
        Vinfo( 1)='SworkL'
        Vinfo( 2)='eigenproblem work array'
        status=def_var(ng, model, GST(ng)%pioFile, varDesc, PIO_FRST,   &
     &                 1, (/SworkLdim/), Aval, Vinfo, ncname)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!-----------------------------------------------------------------------
!  Leave definition mode.
!-----------------------------------------------------------------------
!
        CALL pio_netcdf_enddef (ng, model, ncname, GST(ng)%pioFile)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      END IF DEFINE
!
  10  FORMAT (/,2x,'DEF_GST_PIO      - creating checkpointing file,',   &
     &        t56,'Grid ',i2.2,': ',a)
  20  FORMAT (/,2x,'DEF_GST_PIO      - inquiring checkpointing file,',  &
     &        t56,'Grid ',i2.2,': ',a)
  30  FORMAT (/,' DEF_GST_PIO - unable to create checkpointing',        &
     &        ' NetCDF file: ',a)
!
      RETURN
      END SUBROUTINE def_gst_pio
# endif
#endif
      END MODULE def_gst_mod
