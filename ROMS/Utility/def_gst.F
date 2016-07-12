#include "cppdefs.h"
#if defined PROPAGATOR && defined CHECKPOINTING
      SUBROUTINE def_gst (ng, model)
!
!svn $Id: def_gst.F 795 2016-05-11 01:42:43Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine creates the checkpointing  restart  NetCDF file for    !
!  the Generalized Stability Theory (GST) analysis. It contains all    !
!  the necessary eigenproblem field for restating ARPACK.              !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars
      USE mod_storage
!
      USE def_var_mod, ONLY : def_var
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
!
!  Local variable declarations.
!
      integer, parameter :: Natt = 25

      integer :: Mstatedim, NCVdim, SworkDdim, SworkLdim
      integer :: char1dim, char2dim, iaitrdim, iaupddim, iaup2dim
      integer :: iparamdim, ipntrdim, laitrdim, laup2dim, raitrdim
      integer :: raup2dim
      integer :: i, j, status, varid
      integer :: DimIDs(32), vardim(2)

      integer :: def_dim

      real(r8) :: Aval(6)

      character (len=120) :: Vinfo(Natt)
      character (len=256) :: ncname
!
!-----------------------------------------------------------------------
!  Set and report file name.
!-----------------------------------------------------------------------
!
      IF (exit_flag.ne.NoError) RETURN
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
        IF (exit_flag.ne.NoError) THEN
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
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, model, GST(ng)%ncid, ncname, 'NCV',          &
     &                 NCV, NCVdim)
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, model, GST(ng)%ncid, ncname, 'LworkD',       &
     &                 3*Mstate(ng), SworkDdim)
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, model, GST(ng)%ncid, ncname, 'LworkL',       &
     &                 LworkL, SworkLdim)
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, model, GST(ng)%ncid, ncname, 'iparam',       &
     &                 SIZE(iparam), iparamdim)
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, model, GST(ng)%ncid, ncname, 'ipntr',        &
     &                 SIZE(ipntr), ipntrdim)
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, model, GST(ng)%ncid, ncname, 'iaupd',        &
     &                 SIZE(iaupd), iaupddim)
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, model, GST(ng)%ncid, ncname, 'laitr',        &
     &                 SIZE(laitr), laitrdim)
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, model, GST(ng)%ncid, ncname, 'iaitr',        &
     &                 SIZE(iaitr), iaitrdim)
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, model, GST(ng)%ncid, ncname, 'raitr',        &
     &                 SIZE(raitr), raitrdim)
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, model, GST(ng)%ncid, ncname, 'laup2',        &
     &                 SIZE(laup2), laup2dim)
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, model, GST(ng)%ncid, ncname, 'iaup2',        &
     &                 SIZE(iaup2), iaup2dim)
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, model, GST(ng)%ncid, ncname, 'raup2',        &
     &                 SIZE(raup2), raup2dim)
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, model, GST(ng)%ncid, ncname, 'char2',        &
     &                 2, char2dim)
        IF (exit_flag.ne.NoError) RETURN
!
!-----------------------------------------------------------------------
!  Define global attributes.
!-----------------------------------------------------------------------
!
        CALL def_info (ng, model, GST(ng)%ncid, ncname, DimIDs)
        IF (exit_flag.ne.NoError) RETURN
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

        IF (exit_flag.ne.NoError) RETURN
!
!  Define number of Lanczos vectors to compute.
!
        Vinfo( 1)='NCV'
        Vinfo( 2)='number of Lanczos vectors to compute'
        status=def_var(ng, model, GST(ng)%ncid, varid, nf90_int,        &
     &                 1, (/0/), Aval, Vinfo, ncname,                   &
     &                 SetParAccess = .FALSE.)

        IF (exit_flag.ne.NoError) RETURN
!
!  Define size of eigenvalue problem.
!
        Vinfo( 1)='Mstate'
        Vinfo( 2)='total size of eigenvalue problem'
        status=def_var(ng, model, GST(ng)%ncid, varid, nf90_int,        &
     &                 1, (/0/), Aval, Vinfo, ncname,                   &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define iteration number.
!
        Vinfo( 1)='iter'
        Vinfo( 2)='iteration number'
        status=def_var(ng, model, GST(ng)%ncid, varid, nf90_int,        &
     &                 1, (/0/), Aval, Vinfo, ncname,                   &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define reverse communications flag.
!
        Vinfo( 1)='ido'
        Vinfo( 2)='reverse communications flag'
        status=def_var(ng, model, GST(ng)%ncid, varid, nf90_int,        &
     &                 1, (/0/), Aval, Vinfo, ncname,                   &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define information and error flag.
!
        Vinfo( 1)='info'
        Vinfo( 2)='information and error flag'
        status=def_var(ng, model, GST(ng)%ncid, varid, nf90_int,        &
     &                 1, (/0/), Aval, Vinfo, ncname,                   &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define eigenvalue problem type.
!
        Vinfo( 1)='bmat'
        Vinfo( 2)='eigenvalue problem type'
        status=def_var(ng, model, GST(ng)%ncid, varid, nf90_char,       &
     &                 1, (/0/), Aval, Vinfo, ncname,                   &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define Ritz eigenvalues to compute.
!
        Vinfo( 1)='which'
        Vinfo( 2)='Ritz eigenvalues to compute'
        status=def_var(ng, model, GST(ng)%ncid, varid, nf90_char,       &
     &                 1, (/char2dim/), Aval, Vinfo, ncname,            &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define form of basis function.
!
        Vinfo( 1)='howmany'
        Vinfo( 2)='form of basis function'
        status=def_var(ng, model, GST(ng)%ncid, varid, nf90_char,       &
     &                 1, (/0/), Aval, Vinfo, ncname,                   &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define relative accuracy of computed Ritz values.
!
        Vinfo( 1)='Ritz_tol'
        Vinfo( 2)='relative accuracy of computed Ritz values'
        status=def_var(ng, model, GST(ng)%ncid, varid, NF_FRST,         &
     &                 1, (/0/), Aval, Vinfo, ncname,                   &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define eigenproblem parameters.
!
        Vinfo( 1)='iparam'
        Vinfo( 2)='eigenproblem parameters'
        status=def_var(ng, model, GST(ng)%ncid, varid, nf90_int,        &
     &                 1, (/iparamdim/), Aval, Vinfo, ncname,           &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define pointers to mark starting location in work arrays.
!
        Vinfo( 1)='ipntr'
        Vinfo( 2)='pointers to mark starting location in work arrays'
        status=def_var(ng, model, GST(ng)%ncid,varid, nf90_int,         &
     &                 1, (/ipntrdim/), Aval, Vinfo, ncname,            &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define ARPACK internal integer parameters to _aupd routines.
!
        Vinfo( 1)='iaupd'
        Vinfo( 2)='ARPACK internal integer parameters to _aupd routines'
        status=def_var(ng, model, GST(ng)%ncid, varid, nf90_int,        &
     &                 1, (/iaupddim/), Aval, Vinfo, ncname,            &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define ARPACK internal integer parameters to _aitr routines.
!
        Vinfo( 1)='iaitr'
        Vinfo( 2)='ARPACK internal integer parameters to _aitr routines'
        status=def_var(ng, model, GST(ng)%ncid, varid, nf90_int,        &
     &                 1, (/iaitrdim/), Aval, Vinfo, ncname,            &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define ARPACK internal integer parameters to _aup2 routines.
!
        Vinfo( 1)='iaup2'
        Vinfo( 2)='ARPACK internal integer parameters to _aup2 routines'
        status=def_var(ng, model, GST(ng)%ncid, varid, nf90_int,        &
     &                 1, (/iaup2dim/), Aval, Vinfo, ncname,            &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN
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
        IF (exit_flag.ne.NoError) RETURN
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
        IF (exit_flag.ne.NoError) RETURN
!
!  Define ARPACK internal real parameters to _aitr routines.
!
        Vinfo( 1)='raitr'
        Vinfo( 2)='ARPACK internal real parameters to _aitr routines'
        status=def_var(ng, model, GST(ng)%ncid, varid, NF_FRST,         &
     &                 1, (/raitrdim/), Aval, Vinfo, ncname,            &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define ARPACK internal real parameters to _aup2 routines.
!
        Vinfo( 1)='raup2'
        Vinfo( 2)='ARPACK internal real parameters to _aup2 routines'
        status=def_var(ng, model, GST(ng)%ncid, varid, NF_FRST,         &
     &                 1, (/raup2dim/), Aval, Vinfo, ncname,            &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define Lanczos/Arnoldi basis vectors.
!
        Vinfo( 1)='Bvec'
        Vinfo( 2)='Lanczos/Arnoldi basis vectors'
        vardim(1)=Mstatedim
        vardim(2)=NCVdim
        status=def_var(ng, model, GST(ng)%ncid, varid, NF_FRST,         &
     &                 2, vardim, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define eigenproblem residual vector.
!
        Vinfo( 1)='resid'
        Vinfo( 2)='eigenproblem residual vector'
        status=def_var(ng, model, GST(ng)%ncid, varid, NF_FRST,         &
     &                 1, (/Mstatedim/), Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define state reverse communications work array.
!
        Vinfo( 1)='SworkD'
        Vinfo( 2)='reverse communications state array'
        status=def_var(ng, model, GST(ng)%ncid, varid, NF_FRST,         &
     &                 1, (/SworkDdim/), Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define eigenproblem work array.
!
        Vinfo( 1)='SworkL'
        Vinfo( 2)='eigenproblem work array'
        status=def_var(ng, model, GST(ng)%ncid, varid, NF_FRST,         &
     &                 1, (/SworkLdim/), Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!-----------------------------------------------------------------------
!  Leave definition mode.
!-----------------------------------------------------------------------
!
        CALL netcdf_enddef (ng, model, ncname, GST(ng)%ncid)
        IF (exit_flag.ne.NoError) RETURN

      END IF DEFINE
!
  10  FORMAT (/,6x,'DEF_GST   - creating checkpointing file, Grid ',    &
     &        ': ', a)
  20  FORMAT (/,6x,'DEF_GST   - inquiring checkpointing file, Grid ',   &
     &        ': ', a)
  30  FORMAT (/,' DEF_GST - unable to create checkpointing NetCDF',     &
     &          ' file: ', a)

      RETURN
      END SUBROUTINE def_gst
#else
      SUBROUTINE def_gst
      RETURN
      END SUBROUTINE def_gst
#endif
