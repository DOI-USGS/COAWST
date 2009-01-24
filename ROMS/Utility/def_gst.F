#include "cppdefs.h"
#if defined PROPAGATOR && defined CHECKPOINTING
      SUBROUTINE def_gst (ng, model)
!
!svn $Id: def_gst.F 733 2008-09-07 01:56:45Z jcwarner $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
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
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
!
!  Local variable declarations.
!
      integer, parameter :: Natt = 24

      integer :: Mstatedim, NCVdim, Nnodesdim, SworkDdim, SworkLdim
      integer :: char2dim, iaitrdim, iaupddim, iaup2dim, iparamdim
      integer :: ipntrdim, laitrdim, laup2dim, raitrdim, raup2dim
      integer :: i, j, status, varid
      integer :: DimIDs(29), vardim(2)

      integer :: def_info, def_var

      real(r8) :: Aval(6)

      character (len=80) :: Vinfo(Natt)
      character (len=80) :: ncname
!
!-----------------------------------------------------------------------
!  Set and report file name.
!-----------------------------------------------------------------------
!
      IF (exit_flag.ne.NoError) RETURN
      ncname=GSTname(ng)
!
      IF (Master) THEN
        IF (.not.LrstGST) THEN
          WRITE (stdout,10) TRIM(ncname)
        ELSE
          WRITE (stdout,20) TRIM(ncname)
        END IF
      END IF
!
!=======================================================================
!  Create a new averages NetCDF file.
!=======================================================================
!
      IF (.not.LrstGST.and.OutThread) THEN
        status=nf90_create(TRIM(ncname),nf90_clobber,ncGSTid(ng))
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,30) TRIM(ncname)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Define the dimensions of staggered fields.
!-----------------------------------------------------------------------
!
      IF (.not.LrstGST.and.OutThread) THEN
        status=nf90_def_dim(ncGSTid(ng), 'Mstate', Mstate(ng),          &
     &                      Mstatedim)
        status=nf90_def_dim(ncGSTid(ng), 'NCV'   , NCV,                 &
     &                      NCVdim)
        status=nf90_def_dim(ncGSTid(ng), 'LworkD', 3*Nstate(ng),        &
     &                      SworkDdim)
        status=nf90_def_dim(ncGSTid(ng), 'LworkL', LworkL,              &
     &                      SworkLdim)
# ifdef DISTRIBUTE
        status=nf90_def_dim(ncGSTid(ng), 'Nnodes', numnodes,            &
     &                      Nnodesdim)
# endif
        status=nf90_def_dim(ncGSTid(ng), 'iparam', SIZE(iparam),        &
     &                      iparamdim)
        status=nf90_def_dim(ncGSTid(ng), 'ipntr' , SIZE(ipntr),         &
     &                      ipntrdim)
        status=nf90_def_dim(ncGSTid(ng), 'iaupd' , SIZE(iaupd),         &
     &                      iaupddim)
        status=nf90_def_dim(ncGSTid(ng), 'laitr' , SIZE(laitr),         &
     &                      laitrdim)
        status=nf90_def_dim(ncGSTid(ng), 'iaitr' , SIZE(iaitr),         &
     &                      iaitrdim)
        status=nf90_def_dim(ncGSTid(ng), 'raitr' , SIZE(raitr),         &
     &                      raitrdim)
        status=nf90_def_dim(ncGSTid(ng), 'laup2' , SIZE(laup2),         &
     &                      laup2dim)
        status=nf90_def_dim(ncGSTid(ng), 'iaup2' , SIZE(iaup2),         &
     &                      iaup2dim)
        status=nf90_def_dim(ncGSTid(ng), 'raup2' , SIZE(raup2),         &
     &                      raup2dim)
        status=nf90_def_dim(ncGSTid(ng), 'char2'   ,2,                  &
     &                      char2dim)
!
!-----------------------------------------------------------------------
!  Define global attributes.
!-----------------------------------------------------------------------
!
        status=def_info(ng, ncGSTid(ng), ncname, DimIDs)
        IF (exit_flag.ne.NoError) RETURN
!
!-----------------------------------------------------------------------
!  Define variables and their attributes.
!-----------------------------------------------------------------------
!
!  Initialize local information variable arrays.
!
        DO i=1,Natt
          DO j=1,80
            Vinfo(i)=' '
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
        status=def_var(ncGSTid(ng),varid,nf90_int,1,0,                  &
     &                 Aval,Vinfo,ncname)
!
!  Define number of Lanczos vectors to compute.
!
        Vinfo( 1)='NCV'
        Vinfo( 2)='number of Lanczos vectors to compute'
        status=def_var(ncGSTid(ng),varid,nf90_int,1,0,                  &
     &                 Aval,Vinfo,ncname)
!
!  Define size of eigenvalue problem.
!
        Vinfo( 1)='Mstate'
        Vinfo( 2)='total size of eigenvalue problem'
        status=def_var(ncGSTid(ng),varid,nf90_int,1,0,                  &
     &                 Aval,Vinfo,ncname)

# ifdef DISTRIBUTE
!
!  Define number of distribute memory nodes.
!
        Vinfo( 1)='Nnodes'
        Vinfo( 2)='number of distribute-memory nodes'
        status=def_var(ncGSTid(ng),varid,nf90_int,1,0,                  &
     &                 Aval,Vinfo,ncname)
# endif
!
!  Define iteration number.
!
        Vinfo( 1)='iter'
        Vinfo( 2)='iteration number'
        status=def_var(ncGSTid(ng),varid,nf90_int,1,0,                  &
     &                 Aval,Vinfo,ncname)
!
!  Define reverse communications flag.
!
        Vinfo( 1)='ido'
        Vinfo( 2)='reverse communications flag'
        status=def_var(ncGSTid(ng),varid,nf90_int,1,0,                  &
     &                 Aval,Vinfo,ncname)
!
!  Define information and error flag.
!
        Vinfo( 1)='info'
        Vinfo( 2)='information and error flag'
        status=def_var(ncGSTid(ng),varid,nf90_int,1,0,                  &
     &                 Aval,Vinfo,ncname)
!
!  Define eigenvalue problem type.
!
        Vinfo( 1)='bmat'
        Vinfo( 2)='eigenvalue problem type'
        status=def_var(ncGSTid(ng),varid,nf90_char,1,0,                 &
     &                 Aval,Vinfo,ncname)
!
!  Define Ritz eigenvalues to compute.
!
        Vinfo( 1)='which'
        Vinfo( 2)='Ritz eigenvalues to compute'
        status=def_var(ncGSTid(ng),varid,nf90_char,1,char2dim,          &
     &                 Aval,Vinfo,ncname)
!
!  Define form of basis function.
!
        Vinfo( 1)='howmany'
        Vinfo( 2)='form of basis function'
        status=def_var(ncGSTid(ng),varid,nf90_char,1,0,                 &
     &                 Aval,Vinfo,ncname)
!
!  Define relative accuracy of computed Ritz values.
!
        Vinfo( 1)='Ritz_tol'
        Vinfo( 2)='relative accuracy of computed Ritz values'
        status=def_var(ncGSTid(ng),varid,NF_FRST,1,0,                   &
     &                 Aval,Vinfo,ncname)
!
!  Define eigenproblem parameters.
!
        Vinfo( 1)='iparam'
        Vinfo( 2)='eigenproblem parameters'
        status=def_var(ncGSTid(ng),varid,nf90_int,1,iparamdim,          &
     &                 Aval,Vinfo,ncname)
!
!  Define pointers to mark starting location in work arrays.
!
        Vinfo( 1)='ipntr'
        Vinfo( 2)='pointers to mark starting location in work arrays'
        status=def_var(ncGSTid(ng),varid,nf90_int,1,ipntrdim,           &
     &                 Aval,Vinfo,ncname)
!
!  Define ARPACK internal integer parameters to _aupd routines.
!
        Vinfo( 1)='iaupd'
        Vinfo( 2)='ARPACK internal integer parameters to _aupd routines'
        status=def_var(ncGSTid(ng),varid,nf90_int,1,iaupddim,           &
     &                 Aval,Vinfo,ncname)
!
!  Define ARPACK internal integer parameters to _aitr routines.
!
        Vinfo( 1)='iaitr'
        Vinfo( 2)='ARPACK internal integer parameters to _aitr routines'
        status=def_var(ncGSTid(ng),varid,nf90_int,1,iaitrdim,           &
     &                 Aval,Vinfo,ncname)
!
!  Define ARPACK internal integer parameters to _aup2 routines.
!
        Vinfo( 1)='iaup2'
        Vinfo( 2)='ARPACK internal integer parameters to _aup2 routines'
        status=def_var(ncGSTid(ng),varid,nf90_int,1,iaup2dim,           &
     &                 Aval,Vinfo,ncname)
!
!  Define ARPACK internal logical parameters to _aitr routines.
!
        Vinfo( 1)='laitr'
        Vinfo( 2)='ARPACK internal logical parameters to _aitr routines'
        status=def_var(ncGSTid(ng),varid,nf90_char,1,laitrdim,          &
     &                 Aval,Vinfo,ncname)
!
!  Define ARPACK internal logical parameters to _aup2 routines.
!
        Vinfo( 1)='laup2'
        Vinfo( 2)='ARPACK internal logical parameters to _aup2 routines'
        status=def_var(ncGSTid(ng),varid,nf90_char,1,laup2dim,          &
     &                 Aval,Vinfo,ncname)
!
!  Define ARPACK internal real parameters to _aitr routines.
!
        Vinfo( 1)='raitr'
        Vinfo( 2)='ARPACK internal real parameters to _aitr routines'
        status=def_var(ncGSTid(ng),varid,NF_FRST,1,raitrdim,            &
     &                 Aval,Vinfo,ncname)
!
!  Define ARPACK internal real parameters to _aup2 routines.
!
        Vinfo( 1)='raup2'
        Vinfo( 2)='ARPACK internal real parameters to _aup2 routines'
        status=def_var(ncGSTid(ng),varid,NF_FRST,1,raup2dim,            &
     &                 Aval,Vinfo,ncname)
!
!  Define Lanczos/Arnoldi basis vectors.
!
        Vinfo( 1)='Bvec'
        Vinfo( 2)='Lanczos/Arnoldi basis vectors'
        vardim(1)=Mstatedim
        vardim(2)=NCVdim
        status=def_var(ncGSTid(ng),varid,NF_FRST,2,vardim,              &
     &                 Aval,Vinfo,ncname)
!
!  Define eigenproblem residual vector.
!
        Vinfo( 1)='resid'
        Vinfo( 2)='eigenproblem residual vector'
        status=def_var(ncGSTid(ng),varid,NF_FRST,1,Mstatedim,           &
     &                 Aval,Vinfo,ncname)
!
!  Define state reverse communications work array.
!
        Vinfo( 1)='SworkD'
        Vinfo( 2)='reverse communications state array'
# ifdef DISTRIBUTE
        vardim(1)=SworkDdim
        vardim(2)=Nnodesdim
        status=def_var(ncGSTid(ng),varid,NF_FRST,2,vardim,              &
     &                 Aval,Vinfo,ncname)
# else
        status=def_var(ncGSTid(ng),varid,NF_FRST,1,SworkDdim,           &
     &                 Aval,Vinfo,ncname)
# endif
!
!  Define eigenproblem work array.
!
        Vinfo( 1)='SworkL'
        Vinfo( 2)='eigenproblem work array'
# ifdef DISTRIBUTE
        vardim(1)=SworkLdim
        vardim(2)=Nnodesdim
        status=def_var(ncGSTid(ng),varid,NF_FRST,2,vardim,              &
     &                 Aval,Vinfo,ncname)
# else
        status=def_var(ncGSTid(ng),varid,NF_FRST,1,SworkLdim,           &
     &                 Aval,Vinfo,ncname)
# endif
!
!-----------------------------------------------------------------------
!  Leave definition mode.
!-----------------------------------------------------------------------
!
        status=nf90_enddef(ncGSTid(ng))
      END IF

  10  FORMAT (/,6x,'DEF_GST   - creating checkpointing file: ', a)
  20  FORMAT (/,6x,'DEF_GST   - inquiring checkpointing file: ', a)
  30  FORMAT (/,' DEF_GST - unable to create checkpointing NetCDF',     &
     &          ' file: ', a)

      RETURN
      END SUBROUTINE def_gst
#else
      SUBROUTINE def_gst
      RETURN
      END SUBROUTINE def_gst
#endif
