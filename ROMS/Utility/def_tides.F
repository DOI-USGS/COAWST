#include "cppdefs.h"
#if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
      SUBROUTINE def_tides (ng,ldef)
!
!svn $Id: def_tides.F 733 2008-09-07 01:56:45Z jcwarner $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine creates creates or updates tidal forcing NetCDF file.  !
!  If a file already exists and appropriate,  it will define several   !
!  additional fields needed to detide time-averaged fields.            !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars
      USE mod_stepping
      USE mod_tides
# ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_bcasti
# endif
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng

      logical, intent(in) :: ldef
!
!  Local variable declarations.
!
      integer, parameter :: Natt = 24

      logical :: got_var(NV)

      integer :: i, j, nrec, nvd, nvd3, nvd4
      integer :: status, varid

      integer :: DimIDs(29), tharm(2), t2dgrd(3), u2dgrd(3), v2dgrd(3)
# ifdef SOLVE3D
      integer :: u3dgrd(4), v3dgrd(4)
# endif
      integer :: Vsize(4)

      integer :: def_info, def_var, wrt_info

      real(r8) :: Aval(6)

      character (len=80) :: Vinfo(Natt)
      character (len=80) :: fname, ncname
!
!=======================================================================
!  Open existing tidal forcing file and define new variables.
!=======================================================================
!
      IF (exit_flag.ne.NoError) RETURN
!
      ncname=TIDEname(ng)
      IF (ldef.and.OutThread) THEN
        status=nf90_open(TRIM(ncname),nf90_write,ncTIDEid(ng))
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,10) TRIM(ncname)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# ifdef DISTRIBUTE
      IF (ldef) THEN
        CALL mp_bcasti (ng, iNLM, ncTIDEid(ng), 1)
      END IF
# endif
!
!  Put existing file into define mode so new variables can be added.
!
      IF (ldef.and.OutThread) THEN
        status=nf90_redef(ncTIDEid(ng))
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,20) TRIM(ncname)
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
      IF (ldef.and.OutThread) THEN
        status=nf90_inq_dimid(ncTIDEid(ng),'xi_rho',     DimIDs( 1))
        status=nf90_inq_dimid(ncTIDEid(ng),'xi_u'       ,DimIDs( 2))
        status=nf90_inq_dimid(ncTIDEid(ng),'xi_v'       ,DimIDs( 3))
        status=nf90_inq_dimid(ncTIDEid(ng),'eta_rho'    ,DimIDs( 5))
        status=nf90_inq_dimid(ncTIDEid(ng),'eta_u'      ,DimIDs( 6))
        status=nf90_inq_dimid(ncTIDEid(ng),'eta_v'      ,DimIDs( 7))
        status=nf90_inq_dimid(ncTIDEid(ng),'tide_period',DimIDs( 8))
# ifdef SOLVE3D
        status=nf90_def_dim(ncTIDEid(ng),'s_rho'    ,N(ng),             &
     &                      DimIDs( 9))
# endif
        status=nf90_def_dim(ncTIDEid(ng),'harmonics',2*NTC(ng)+1,       &
     &                      DimIDs(12))
!
!  Set number of dimensions for output variables.
!
        nvd3=3
        nvd4=4
!
!  Define dimension vectors for tide harmonics variables.
!
        tharm(1)=DimIDs(8)
        tharm(2)=DimIDs(8)
!
!  Define dimension vectors for staggered tracer type variables.
!
        t2dgrd(1)=DimIDs( 1)
        t2dgrd(2)=DimIDs( 5)
        t2dgrd(3)=DimIDs(12)
!
!  Define dimension vectors for staggered u-momemtum type variables.
!
        u2dgrd(1)=DimIDs( 2)
        u2dgrd(2)=DimIDs( 6)
        u2dgrd(3)=DimIDs(12)
# ifdef SOLVE3D
        u3dgrd(1)=DimIDs( 2)
        u3dgrd(2)=DimIDs( 6)
        u3dgrd(3)=DimIDs( 9)
        u3dgrd(4)=DimIDs(12)
# endif
!
!  Define dimension vectors for staggered v-momemtum type variables.
!
        v2dgrd(1)=DimIDs( 3)
        v2dgrd(2)=DimIDs( 7)
        v2dgrd(3)=DimIDs(12)
# ifdef SOLVE3D
        v3dgrd(1)=DimIDs( 3)
        v3dgrd(2)=DimIDs( 7)
        v3dgrd(3)=DimIDs( 9)
        v3dgrd(4)=DimIDs(12)
# endif
!
!  Initialize local information variable arrays.
!
        DO i=1,Natt
          DO j=1,80
            Vinfo(i)(j:j)=' '
          END DO
        END DO
        DO i=1,6
          Aval(i)=0.0_r8
        END DO
!
!-----------------------------------------------------------------------
!  Define variables and their attributes.
!-----------------------------------------------------------------------
!
!  Define number of time-accumulated harmonics.
!
        Vinfo( 1)='Hcount'
        Vinfo( 2)='number of time-accumulated tide harmonics'
        status=def_var(ncTIDEid(ng),varid,nf90_int,1,0,                 &
     &                 Aval,Vinfo,ncname)
!  
!  Define model time for accumulated tide harmonic fields.
!
        Vinfo( 1)=Vname(1,idtime)
        WRITE (Vinfo( 2),'(a,1x,a)') 'accumulated harmonics',           &
     &                               TRIM(Vname(2,idtime))
        IF (INT(time_ref).eq.-2) THEN
          Vinfo( 3)='seconds since 1968-05-23 00:00:00 GMT'
          Vinfo( 4)='gregorian'
        ELSE IF (INT(time_ref).eq.-1) THEN
          Vinfo( 3)='seconds since 0001-01-01 00:00:00'
          Vinfo( 4)='360_day'
        ELSE IF (INT(time_ref).eq.0) THEN
          Vinfo( 3)='seconds since 0001-01-01 00:00:00'
          Vinfo( 4)='julian'
        ELSE IF (time_ref.gt.0.0_r8) THEN
          WRITE (Vinfo( 3),'(a,1x,a)') 'seconds since', TRIM(r_text)
          Vinfo( 4)='gregorian'
        END IF
        status=def_var(ncTIDEid(ng),tideVid(idtime,ng),NF_FRST,         &
     &                 1,0,Aval,Vinfo,ncname)
!
!  Define time-accumulated COS(omega(k)*t) harmonics.
!
        Vinfo( 1)=Vname(1,idCosW)
        Vinfo( 2)=Vname(2,idCosW)
        Vinfo( 3)=Vname(3,idCosW)
        status=def_var(ncTIDEid(ng),tideVid(idCosW,ng),NF_FRST,         &
     &                 1,tharm(1),Aval,Vinfo,ncname)
!
!  Define time-accumulated SIN(omega(k)*t) harmonics.
!
        Vinfo( 1)=Vname(1,idSinW)
        Vinfo( 2)=Vname(2,idSinW)
        Vinfo( 3)=Vname(3,idSinW)
        status=def_var(ncTIDEid(ng),tideVid(idSinW,ng),NF_FRST,         &
     &                 1,tharm(1),Aval,Vinfo,ncname)
!
!  Define time-accumulated COS(omega(k)*t)*COS(omega(l)*t) harmonics.
!
        Vinfo( 1)=Vname(1,idCos2)
        Vinfo( 2)=Vname(2,idCos2)
        Vinfo( 3)=Vname(3,idCos2)
        status=def_var(ncTIDEid(ng),tideVid(idCos2,ng),NF_FRST,         &
     &                 2,tharm,Aval,Vinfo,ncname)
!
!  Define time-accumulated SIN(omega(k)*t)*SIN(omega(l)*t) harmonics.
!
        Vinfo( 1)=Vname(1,idSin2)
        Vinfo( 2)=Vname(2,idSin2)
        Vinfo( 3)=Vname(3,idSin2)
        status=def_var(ncTIDEid(ng),tideVid(idSin2,ng),NF_FRST,         &
     &                 2,tharm,Aval,Vinfo,ncname)
!
!  Define time-accumulated SIN(omega(k)*t)*COS(omega(l)*t) harmonics.
!
        Vinfo( 1)=Vname(1,idSWCW)
        Vinfo( 2)=Vname(2,idSWCW)
        Vinfo( 3)=Vname(3,idSWCW)
        status=def_var(ncTIDEid(ng),tideVid(idSWCW,ng),NF_FRST,         &
     &                 2,tharm,Aval,Vinfo,ncname)
!
!  Define free-surface time-accumulated tide harmonics.
!
        Vinfo( 1)=Vname(1,idFsuH)
        Vinfo( 2)=Vname(2,idFsuH)
        Vinfo( 3)=Vname(3,idFsuH)
        Vinfo(14)=Vname(4,idFsuH)
# if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_rho'
# endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idFsuH,ng),r8)
        status=def_var(ncTIDEid(ng),tideVid(idFsuH,ng),NF_FRST,         &
     &                 nvd3,t2dgrd,Aval,Vinfo,ncname)
!
!  Define 2D u-momentum time-accumulated tide harmonics.
!
        Vinfo( 1)=Vname(1,idu2dH)
        Vinfo( 2)=Vname(2,idu2dH)
        Vinfo( 3)=Vname(3,idu2dH)
        Vinfo(14)=Vname(4,idu2dH)
# if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_u'
# endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idu2dH,ng),r8)
        status=def_var(ncTIDEid(ng),tideVid(idu2dH,ng),NF_FRST,         &
     &                 nvd3,u2dgrd,Aval,Vinfo,ncname)
!
!  Define 2D v-momentum time-accumulated tide harmonics.
!
        Vinfo( 1)=Vname(1,idv2dH)
        Vinfo( 2)=Vname(2,idv2dH)
        Vinfo( 3)=Vname(3,idv2dH)
        Vinfo(14)=Vname(4,idv2dH)
# if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_v'
# endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idv2dH,ng),r8)
        status=def_var(ncTIDEid(ng),tideVid(idv2dH,ng),NF_FRST,         &
     &                 nvd3,v2dgrd,Aval,Vinfo,ncname)

# ifdef SOLVE3D
!
!  Define 3D u-momentum time-accumulated tide harmonics.
!
        Vinfo( 1)=Vname(1,idu3dH)
        Vinfo( 2)=Vname(2,idu3dH)
        Vinfo( 3)=Vname(3,idu3dH)
        Vinfo(14)=Vname(4,idu3dH)
#  if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_u'
#  endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idu3dH,ng),r8)
        status=def_var(ncTIDEid(ng),tideVid(idu3dH,ng),NF_FRST,         &
     &                 nvd4,u3dgrd,Aval,Vinfo,ncname)
!
!  Define 3D v-momentum time-accumulated tide harmonics.
!
        Vinfo( 1)=Vname(1,idv3dH)
        Vinfo( 2)=Vname(2,idv3dH)
        Vinfo( 3)=Vname(3,idv3dH)
        Vinfo(14)=Vname(4,idv3dH)
#  if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_v'
#  endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idv3dH,ng),r8)
        status=def_var(ncTIDEid(ng),tideVid(idv3dH,ng),NF_FRST,         &
     &                 nvd4,v3dgrd,Aval,Vinfo,ncname)
# endif
!
!-----------------------------------------------------------------------
!  Leave definition mode.
!-----------------------------------------------------------------------
!
        status=nf90_enddef(ncTIDEid(ng))
      END IF
!
!=======================================================================
!  Open an existing averages file, check its contents, and prepare
!  for appending data.
!=======================================================================
!
      IF (.not.ldef.and.InpThread) THEN
!
!  Inquire about the contents of averages NetCDF file:  Inquire about
!  the dimensions and variables.  Check for consistency.
!
        CALL opencdf (ng, 1, ncname, fname, N(ng), 0, nrec, nvd, Vsize)
        IF (exit_flag.ne.NoError) RETURN
!
!  Open  file for read/write.
!
        status=nf90_open(TRIM(ncname), nf90_write, ncTIDEid(ng))
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,30) TRIM(ncname)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
!
!  Initialize logical switches.
!
        DO i=1,NV
          got_var(i)=.FALSE.
        END DO
!
!  Scan variable list from input NetCDF and activate switches for
!  tide variables. Get variable IDs.
!
        DO i=1,nvars
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idtime))) THEN
            got_var(idtime)=.TRUE.
            status=nf90_inq_varid(ncTIDEid(ng),TRIM(Vname(1,idtime)),   &
     &                            tideVid(idtime,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idCosW))) THEN
            got_var(idCosW)=.TRUE.
            status=nf90_inq_varid(ncTIDEid(ng),TRIM(Vname(1,idCosW)),   &
     &                            tideVid(idCosW,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idSinW))) THEN
            got_var(idSinW)=.TRUE.
            status=nf90_inq_varid(ncTIDEid(ng),TRIM(Vname(1,idSinW)),   &
     &                            tideVid(idSinW,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idCos2))) THEN
            got_var(idCos2)=.TRUE.
            status=nf90_inq_varid(ncTIDEid(ng),TRIM(Vname(1,idCos2)),   &
     &                            tideVid(idCos2,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idSin2))) THEN
            got_var(idSin2)=.TRUE.
            status=nf90_inq_varid(ncTIDEid(ng),TRIM(Vname(1,idSin2)),   &
     &                            tideVid(idSin2,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idSWCW))) THEN
            got_var(idSWCW)=.TRUE.
            status=nf90_inq_varid(ncTIDEid(ng),TRIM(Vname(1,idSWCW)),   &
     &                            tideVid(idSWCW,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idFsuH))) THEN
            got_var(idFsuH)=.TRUE.
            status=nf90_inq_varid(ncTIDEid(ng),TRIM(Vname(1,idFsuH)),   &
     &                            tideVid(idFsuH,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idu2dH))) THEN
            got_var(idu2dH)=.TRUE.
            status=nf90_inq_varid(ncTIDEid(ng),TRIM(Vname(1,idu2dH)),   &
     &                            tideVid(idu2dH,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idv2dH))) THEN
            got_var(idv2dH)=.TRUE.
            status=nf90_inq_varid(ncTIDEid(ng),TRIM(Vname(1,idv2dH)),   &
     &                            tideVid(idv2dH,ng))
          END IF
# ifdef SOLVE3D
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idu3dH))) THEN
            got_var(idu3dH)=.TRUE.
            status=nf90_inq_varid(ncTIDEid(ng),TRIM(Vname(1,idu3dH)),   &
     &                            tideVid(idu3dH,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idv3dH))) THEN
            got_var(idv3dH)=.TRUE.
            status=nf90_inq_varid(ncTIDEid(ng),TRIM(Vname(1,idv3dH)),    &
     &                            tideVid(idv3dH,ng))
          END IF
# endif
        END DO
!
!  Check if tide variables are available in input NetCDF file.
!
        IF (.not.got_var(idtime)) THEN
          WRITE (stdout,30) TRIM(Vname(1,idtime)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idCosW)) THEN
          WRITE (stdout,30) TRIM(Vname(1,idCosW)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idSinW)) THEN
          WRITE (stdout,30) TRIM(Vname(1,idSinW)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idCos2)) THEN
          WRITE (stdout,30) TRIM(Vname(1,idCos2)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idSin2)) THEN
          WRITE (stdout,30) TRIM(Vname(1,idSin2)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idSWCW)) THEN
          WRITE (stdout,30) TRIM(Vname(1,idSWCW)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idFsuH)) THEN
          WRITE (stdout,30) TRIM(Vname(1,idFsuH)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idu2dH)) THEN
          WRITE (stdout,30) TRIM(Vname(1,idu2dH)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idv2dH)) THEN
          WRITE (stdout,30) TRIM(Vname(1,idv2dH)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# ifdef SOLVE3D
        IF (.not.got_var(idu3dH)) THEN
          WRITE (stdout,30) TRIM(Vname(1,idu3dH)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idv3dH)) THEN
          WRITE (stdout,30) TRIM(Vname(1,idv3dH)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
      END IF
!
  10  FORMAT (/,' DEF_TIDES - unable to open tide NetCDF file: ',a)
  20  FORMAT (/,' DEF_TIDES - unable to put in define mode tide NetCDF',  &
     &        ' file: ',a)
  30  FORMAT (/,' DEF_TIDES - unable to find variable: ',a,2x,            &
     &        ' in tide NetCDF file: ',a)
      RETURN
      END SUBROUTINE def_tides
#else
      SUBROUTINE def_tides
      RETURN
      END SUBROUTINE def_tides
#endif
