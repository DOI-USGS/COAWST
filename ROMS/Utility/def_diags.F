#include "cppdefs.h"
#ifdef DIAGNOSTICS
      SUBROUTINE def_diags (ng,ldef)
!
!svn $Id: def_diags.F 733 2008-09-07 01:56:45Z jcwarner $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine creates diagnostics NetCDF file, it defines its        !
!  dimensions, attributes, and variables.                              !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
# ifdef FOUR_DVAR
      USE mod_fourdvar
# endif
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars
# ifdef SEDIMENT
      USE mod_sediment
# endif
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

      integer :: i, ifield, itrc, ivar, j, nrec, nvd, nvd3, nvd4
      integer :: recdim, status, ofst

      integer :: DimIDs(29), t2dgrd(3), u2dgrd(3), v2dgrd(3)
      integer :: Vsize(4)

      integer :: def_info, def_var, wrt_info

# ifdef SOLVE3D
#  ifdef SEDIMENT
      integer :: b3dgrd(4)
#  endif
      integer :: t3dgrd(4), u3dgrd(4), v3dgrd(4), w3dgrd(4)
# endif

      real(r8) :: Aval(6)

      character (len=80) :: Vinfo(Natt)
      character (len=80) :: fname, ncname
!
!-----------------------------------------------------------------------
!  Set and report file name.
!-----------------------------------------------------------------------
!
      IF (exit_flag.ne.NoError) RETURN
      ncname=DIAname(ng)
!
      IF (Master) THEN
        IF (ldef) THEN
          WRITE (stdout,10) TRIM(ncname)
        ELSE
          WRITE (stdout,20) TRIM(ncname)
        END IF
      END IF
!
!=======================================================================
!  Create a new diagnostics NetCDF file.
!=======================================================================
!
      IF (ldef.and.OutThread) THEN
        status=nf90_create(TRIM(ncname),nf90_clobber,ncDIAid(ng))
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,30) TRIM(ncname)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#ifdef DISTRIBUTE
      IF (ldef) THEN
        CALL mp_bcasti (ng, iNLM, ncDIAid(ng), 1)
      END IF
#endif
!
!-----------------------------------------------------------------------
!  Define the dimensions of staggered fields.
!-----------------------------------------------------------------------
!
#ifdef COMPOSED_GRID
            ofst=6
#elif defined REFINED_GRID
      IF (ng.eq.1) THEN
        ofst=1
      ELSE
        ofst=6
      ENDIF
#else
            ofst=1
#endif
      IF (ldef.and.OutThread) THEN
        status=nf90_def_dim(ncDIAid(ng),'xi_rho'  ,Lm(ng)+ofst+1,       &
     &                      DimIDs( 1))
        status=nf90_def_dim(ncDIAid(ng),'xi_u'    ,Lm(ng)+ofst,         &
     &                      DimIDs( 2))
        status=nf90_def_dim(ncDIAid(ng),'xi_v'    ,Lm(ng)+ofst+1,       &
     &                      DimIDs( 3))
        status=nf90_def_dim(ncDIAid(ng),'xi_psi'  ,Lm(ng)+ofst,         &
     &                      DimIDs( 4))
        status=nf90_def_dim(ncDIAid(ng),'eta_rho' ,Mm(ng)+ofst+1,       &
     &                      DimIDs( 5))
        status=nf90_def_dim(ncDIAid(ng),'eta_u'   ,Mm(ng)+ofst+1,       &
     &                      DimIDs( 6))
        status=nf90_def_dim(ncDIAid(ng),'eta_v'   ,Mm(ng)+ofst,         &
     &                      DimIDs( 7))
        status=nf90_def_dim(ncDIAid(ng),'eta_psi' ,Mm(ng)+ofst,         &
     &                      DimIDs( 8))
# if defined WRITE_WATER && defined MASKING
        status=nf90_def_dim(ncDIAid(ng),'xy_rho'  ,Nxyr(ng),            &
     &                      DimIDs(17))
        status=nf90_def_dim(ncDIAid(ng),'xy_u'    ,Nxyu(ng),            &
     &                      DimIDs(18))
        status=nf90_def_dim(ncDIAid(ng),'xy_v'    ,Nxyv(ng),            &
     &                      DimIDs(19))
# endif
# ifdef SOLVE3D
#  if defined WRITE_WATER && defined MASKING
        status=nf90_def_dim(ncDIAid(ng),'xyz_rho' ,Nxyr(ng)*N(ng),      &
     &                      DimIDs(20))
        status=nf90_def_dim(ncDIAid(ng),'xyz_u'   ,Nxyu(ng)*N(ng),      &
     &                      DimIDs(21))
        status=nf90_def_dim(ncDIAid(ng),'xyz_v'   ,Nxyv(ng)*N(ng),      &
     &                      DimIDs(22))
        status=nf90_def_dim(ncDIAid(ng),'xyz_w'   ,Nxyr(ng)*(N(ng)+1),  &
     &                      DimIDs(23))
#  endif
        status=nf90_def_dim(ncDIAid(ng),'s_rho'   ,N(ng),               &
     &                      DimIDs( 9))
        status=nf90_def_dim(ncDIAid(ng),'s_w'     ,N(ng)+1,             &
     &                      DimIDs(10))
        status=nf90_def_dim(ncDIAid(ng),'tracer'  ,NT(ng),              &
     &                      DimIDs(11))
#  ifdef SEDIMENT
        status=nf90_def_dim(ncDIAid(ng),'Nbed'    ,Nbed,                &
     &                      DimIDs(16))
#   if defined WRITE_WATER && defined MASKING
        status=nf90_def_dim(ncDIAid(ng),'xybed'   ,Nxyr(ng)*Nbed,       &
     &                      DimIDs(24))
#   endif
#  endif
#  ifdef ECOSIM
        status=nf90_def_dim(ncDIAid(ng),'Nphy'    ,Nphy,                &
     &                      DimIDs(25))
        status=nf90_def_dim(ncDIAid(ng),'Nbac'    ,Nbac,                &
     &                      DimIDs(26))
        status=nf90_def_dim(ncDIAid(ng),'Ndom'    ,Ndom,                &
     &                      DimIDs(27))
        status=nf90_def_dim(ncDIAid(ng),'Nfec'    ,Nfec,                &
     &                      DimIDs(28))
#  endif
# endif
        status=nf90_def_dim(ncDIAid(ng),'boundary',4,                   &
     &                      DimIDs(14))
# ifdef FOUR_DVAR
        status=nf90_def_dim(ncDIAid(ng),'Nstate  ',NstateVar(ng),       &
     &                      DimIDs(29))
# endif
        status=nf90_def_dim(ncDIAid(ng),TRIM(ADJUSTL(Vname(5,idtime))), &
     &                      nf90_unlimited,DimIDs(12))
        recdim=DimIDs(12)
!
!  Set number of dimensions for output variables.
!
#if defined WRITE_WATER && defined MASKING
        nvd3=2
        nvd4=2
#else
        nvd3=3
        nvd4=4
#endif
!
!  Define dimension vectors for staggered tracer type variables.
!
# if defined WRITE_WATER && defined MASKING
        t2dgrd(1)=DimIDs(17)
        t2dgrd(2)=DimIDs(12)
#  ifdef SOLVE3D
        t3dgrd(1)=DimIDs(20)
        t3dgrd(2)=DimIDs(12)
#  endif
# else
        t2dgrd(1)=DimIDs( 1)
        t2dgrd(2)=DimIDs( 5)
        t2dgrd(3)=DimIDs(12)
#  ifdef SOLVE3D
        t3dgrd(1)=DimIDs( 1)
        t3dgrd(2)=DimIDs( 5)
        t3dgrd(3)=DimIDs( 9)
        t3dgrd(4)=DimIDs(12)
#  endif
# endif
!
!  Define dimension vectors for staggered u-momemtum type variables.
!
# if defined WRITE_WATER && defined MASKING
        u2dgrd(1)=DimIDs(18)
        u2dgrd(2)=DimIDs(12)
#  ifdef SOLVE3D
        u3dgrd(1)=DimIDs(21)
        u3dgrd(2)=DimIDs(12)
#  endif
# else
        u2dgrd(1)=DimIDs( 2)
        u2dgrd(2)=DimIDs( 6)
        u2dgrd(3)=DimIDs(12)
#  ifdef SOLVE3D
        u3dgrd(1)=DimIDs( 2)
        u3dgrd(2)=DimIDs( 6)
        u3dgrd(3)=DimIDs( 9)
        u3dgrd(4)=DimIDs(12)
#  endif
# endif
!
!  Define dimension vectors for staggered v-momemtum type variables.
!
# if defined WRITE_WATER && defined MASKING
        v2dgrd(1)=DimIDs(19)
        v2dgrd(2)=DimIDs(12)
#  ifdef SOLVE3D
        v3dgrd(1)=DimIDs(22)
        v3dgrd(2)=DimIDs(12)
#  endif
# else
        v2dgrd(1)=DimIDs( 3)
        v2dgrd(2)=DimIDs( 7)
        v2dgrd(3)=DimIDs(12)
#  ifdef SOLVE3D
        v3dgrd(1)=DimIDs( 3)
        v3dgrd(2)=DimIDs( 7)
        v3dgrd(3)=DimIDs( 9)
        v3dgrd(4)=DimIDs(12)
#  endif
# endif
# ifdef SOLVE3D
!
!  Define dimension vector for staggered w-momemtum type variables.
!
#  if defined WRITE_WATER && defined MASKING
        w3dgrd(1)=DimIDs(23)
        w3dgrd(2)=DimIDs(12)
#  else
        w3dgrd(1)=DimIDs( 1)
        w3dgrd(2)=DimIDs( 5)
        w3dgrd(3)=DimIDs(10)
        w3dgrd(4)=DimIDs(12)
#  endif
#  ifdef SEDIMENT
!
!  Define dimension vector for sediment bed layer type variables.
!
#   if defined WRITE_WATER && defined MASKING
        b3dgrd(1)=DimIDs(24)
        b3dgrd(2)=DimIDs(12)
#   else
        b3dgrd(1)=DimIDs( 1)
        b3dgrd(2)=DimIDs( 5)
        b3dgrd(3)=DimIDs(16)
        b3dgrd(4)=DimIDs(12)
#   endif
#  endif
# endif
!
!  Initialize unlimited time record dimension.
!
        tDIAindx(ng)=0
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
!  Define time-recordless information variables.
!-----------------------------------------------------------------------
!
        status=def_info(ng,ncDIAid(ng),ncname,DimIDs)
        IF (exit_flag.ne.NoError) RETURN
!
!-----------------------------------------------------------------------
!  Define variables and their attributes.
!-----------------------------------------------------------------------
!
!  Define model time.
!
        Vinfo( 1)=Vname(1,idtime)
        WRITE (Vinfo( 2),'(a,1x,a)') 'averaged', TRIM(Vname(2,idtime))
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
        Vinfo(14)=Vname(4,idtime)
        status=def_var(ncDIAid(ng),diaVid(idtime,ng),NF_TYPE,           &
     &                 1,recdim,Aval,Vinfo,ncname)

# ifdef DIAGNOSTICS_UV
!
!  Define 2D momentum diagnostic fields.
!
        DO ivar=1,NDM2d
          ifield=idDu2d(ivar)
          Vinfo( 1)=Vname(1,ifield)
          Vinfo( 2)=Vname(2,ifield)
          Vinfo( 3)=Vname(3,ifield)
          Vinfo(14)=Vname(4,ifield)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(u2dvar,r8)
          status=def_var(ncDIAid(ng),diaVid(ifield,ng),NF_FOUT,         &
     &                   nvd3,u2dgrd,Aval,Vinfo,ncname)
          ifield=idDv2d(ivar)
          Vinfo( 1)=Vname(1,ifield)
          Vinfo( 2)=Vname(2,ifield)
          Vinfo( 3)=Vname(3,ifield)
          Vinfo(14)=Vname(4,ifield)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(v2dvar,r8)
          status=def_var(ncDIAid(ng),diaVid(ifield,ng),NF_FOUT,         &
     &                   nvd3,v2dgrd,Aval,Vinfo,ncname)
        END DO
#  ifdef SOLVE3D
!
!  Define 3D momentum diagnostic fields.
!
        DO ivar=1,NDM3d
          ifield=idDu3d(ivar)
          Vinfo( 1)=Vname(1,ifield)
          Vinfo( 2)=Vname(2,ifield)
          Vinfo( 3)=Vname(3,ifield)
          Vinfo(14)=Vname(4,ifield)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
#   endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(u3dvar,r8)
          status=def_var(ncDIAid(ng),diaVid(ifield,ng),NF_FOUT,         &
     &                   nvd4,u3dgrd,Aval,Vinfo,ncname)
          ifield=idDv3d(ivar)
          Vinfo( 1)=Vname(1,ifield)
          Vinfo( 2)=Vname(2,ifield)
          Vinfo( 3)=Vname(3,ifield)
          Vinfo(14)=Vname(4,ifield)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
#   endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(v3dvar,r8)
          status=def_var(ncDIAid(ng),diaVid(ifield,ng),NF_FOUT,         &
     &                   nvd4,v3dgrd,Aval,Vinfo,ncname)
        END DO
#  endif
# endif

# ifdef DIAGNOSTICS_TS
!
!  Define tracer diagnostic fields.
!
        DO itrc=1,NT(ng)
          DO ivar=1,NDT
            ifield=idDtrc(itrc,ivar)
            Vinfo( 1)=Vname(1,ifield)
            Vinfo( 2)=Vname(2,ifield)
            Vinfo( 3)=Vname(3,ifield)
            Vinfo(14)=Vname(4,ifield)
            Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
            Vinfo(20)='mask_rho'
#  endif
            Vinfo(22)='coordinates'
            Aval(5)=REAL(r3dvar,r8)
            status=def_var(ncDIAid(ng),diaVid(ifield,ng),NF_FOUT,       &
     &                     nvd4,t3dgrd,Aval,Vinfo,ncname)
          END DO
        END DO
# endif

# ifdef BIO_FASHAM
!
!  Define 2D biological diagnostic fields.
!
        DO ivar=1,NDbio2d
          ifield=iDbio2(ivar)
          IF (Hout(ifield,ng)) THEN
            Vinfo( 1)=Vname(1,ifield)
            WRITE (Vinfo( 2),'(a,1x,a)') 'averaged',                    &
     &                                   TRIM(Vname(2,ifield))
            Vinfo( 3)=Vname(3,ifield)
            Vinfo(14)=Vname(4,ifield)
            Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
            Vinfo(20)='mask_rho'
#  endif
            Vinfo(22)='coordinates'
            Aval(5)=REAL(r3dvar,r8)
            status=def_var(ncDIAid(ng),diaVid(ifield,ng),NF_FOUT,       &
     &                     nvd3,t2dgrd,Aval,Vinfo,ncname)
          END IF
        END DO
!
!  Define 3D biological diagnostic fields.
!
        DO ivar=1,NDbio3d
          ifield=iDbio3(ivar)
          IF (Hout(ifield,ng)) THEN
            Vinfo( 1)=Vname(1,ifield)
            WRITE (Vinfo( 2),'(a,1x,a)') 'averaged',                    &
     &                                   TRIM(Vname(2,ifield))
            Vinfo( 3)=Vname(3,ifield)
            Vinfo(14)=Vname(4,ifield)
            Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
            Vinfo(20)='mask_rho'
#  endif
            Vinfo(22)='coordinates'
            Aval(5)=REAL(r3dvar,r8)
            status=def_var(ncDIAid(ng),diaVid(ifield,ng),NF_FOUT,       &
     &                     nvd4,t3dgrd,Aval,Vinfo,ncname)
          END IF
        END DO
# endif
!
!-----------------------------------------------------------------------
!  Leave definition mode.
!-----------------------------------------------------------------------
!
        status=nf90_enddef(ncDIAid(ng))
      END IF
!
!-----------------------------------------------------------------------
!  Write out time-recordless, information variables.
!-----------------------------------------------------------------------
!
      IF (ldef) THEN
        status=wrt_info(ng, iNLM, ncDIAid(ng), OutThread, ncname)
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!=======================================================================
!  Open an existing diagnostics file, check its contents, and prepare
!  for appending data.
!=======================================================================
!
      IF (.not.ldef.and.OutThread) THEN
!
!  Inquire about the contents of diagnostics NetCDF file: Inquire about
!  the dimensions and variables.  Check for consistency.
!
        CALL opencdf (ng, 1, ncname, fname, N(ng), 0, nrec, nvd, Vsize)
        IF (exit_flag.ne.NoError) RETURN
!
!  Open diagnostics file for read/write.
!
        status=nf90_open(TRIM(ncname), nf90_write, ncDIAid(ng))
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,50) TRIM(ncname)
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
!  diagnostics variables. Get variable IDs.
!
        DO i=1,nvars
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idtime))) THEN
            got_var(idtime)=.TRUE.
            status=nf90_inq_varid(ncDIAid(ng),TRIM(Vname(1,idtime)),    &
     &                            diaVid(idtime,ng))
          END IF
# ifdef DIAGNOSTICS_UV
          DO ivar=1,NDM2d
            ifield=idDu2d(ivar)
            IF (TRIM(varnam(i)).eq.TRIM(Vname(1,ifield))) THEN
             got_var(ifield)=.TRUE.
             status=nf90_inq_varid(ncDIAid(ng),TRIM(Vname(1,ifield)),   &
     &                             diaVid(ifield,ng))
            END IF
          END DO
          DO ivar=1,NDM2d
            ifield=idDv2d(ivar)
            IF (TRIM(varnam(i)).eq.TRIM(Vname(1,ifield))) THEN
             got_var(ifield)=.TRUE.
             status=nf90_inq_varid(ncDIAid(ng),TRIM(Vname(1,ifield)),   &
     &                             diaVid(ifield,ng))
            END IF
          END DO
#  ifdef SOLVE3D
          DO ivar=1,NDM3d
            ifield=idDu3d(ivar)
            IF (TRIM(varnam(i)).eq.TRIM(Vname(1,ifield))) THEN
             got_var(ifield)=.TRUE.
             status=nf90_inq_varid(ncDIAid(ng),TRIM(Vname(1,ifield)),   &
     &                             diaVid(ifield,ng))
            END IF
          END DO
          DO ivar=1,NDM3d
            ifield=idDv3d(ivar)
            IF (TRIM(varnam(i)).eq.TRIM(Vname(1,ifield))) THEN
             got_var(ifield)=.TRUE.
             status=nf90_inq_varid(ncDIAid(ng),TRIM(Vname(1,ifield)),   &
     &                             diaVid(ifield,ng))
            END IF
          END DO
#  endif
# endif
# ifdef DIAGNOSTICS_TS
          DO itrc=1,NT(ng)
            DO ivar=1,NDT
              ifield=idDtrc(itrc,ivar)
              IF (TRIM(varnam(i)).eq.TRIM(Vname(1,ifield))) THEN
               got_var(ifield)=.TRUE.
               status=nf90_inq_varid(ncDIAid(ng),TRIM(Vname(1,ifield)), &
     &                               diaVid(ifield,ng))
              END IF
            END DO
          END DO
# endif
# ifdef BIO_FASHAM
          DO ivar=1,NDbio2d
            ifield=iDbio2(ivar)
            IF (TRIM(varnam(i)).eq.TRIM(Vname(1,ifield))) THEN
             got_var(ifield)=.TRUE.
             status=nf90_inq_varid(ncDIAid(ng),TRIM(Vname(1,ifield)),   &
     &                             diaVid(ifield,ng))
            END IF
          END DO
          DO ivar=1,NDbio3d
            ifield=iDbio3(ivar)
            IF (TRIM(varnam(i)).eq.TRIM(Vname(1,ifield))) THEN
             got_var(ifield)=.TRUE.
             status=nf90_inq_varid(ncDIAid(ng),TRIM(Vname(1,ifield)),   &
     &                             diaVid(ifield,ng))
            END IF
          END DO
# endif
        END DO
!
!  Check if diagnostics variables are available in input NetCDF file.
!
        IF (.not.got_var(idtime)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idtime)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# ifdef DIAGNOSTICS_UV
        DO ivar=1,NDM2d
          ifield=idDu2d(ivar)
          IF (.not.got_var(ifield)) THEN
            WRITE (stdout,60) TRIM(Vname(1,ifield)), TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
        DO ivar=1,NDM2d
          ifield=idDv2d(ivar)
          IF (.not.got_var(ifield)) THEN
            WRITE (stdout,60) TRIM(Vname(1,ifield)), TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
#  ifdef SOLVE3D
        DO ivar=1,NDM3d
          ifield=idDu3d(ivar)
          IF (.not.got_var(ifield)) THEN
            WRITE (stdout,60) TRIM(Vname(1,ifield)), TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
        DO ivar=1,NDM3d
          ifield=idDv3d(ivar)
          IF (.not.got_var(ifield)) THEN
            WRITE (stdout,60) TRIM(Vname(1,ifield)), TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
#  endif
# endif
# ifdef DIAGNOSTICS_TS
        DO itrc=1,NT(ng)
          DO ivar=1,NDT
            ifield=idDtrc(itrc,ivar)
            IF (.not.got_var(ifield)) THEN
              WRITE (stdout,60) TRIM(Vname(1,ifield)), TRIM(ncname)
              exit_flag=3
              RETURN
            END IF
          END DO
        END DO
# endif
# ifdef BIO_FASHAM
        DO ivar=1,NDbio2d
          ifield=iDbio2(ivar)
          IF (.not.got_var(ifield).and.Hout(ifield,ng)) THEN
            WRITE (stdout,60) TRIM(Vname(1,ifield)), TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
        DO ivar=1,NDbio3d
          ifield=iDbio3(ivar)
          IF (.not.got_var(ifield).and.Hout(ifield,ng)) THEN
            WRITE (stdout,60) TRIM(Vname(1,ifield)), TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
# endif
!
!  Set unlimited time record dimension to the appropriate value.
!
        IF (nRST(ng).eq.nDIA(ng)) THEN
          IF (ndefDIA(ng).gt.0) THEN
            tDIAindx(ng)=((ntstart(ng)-1)-                              &
     &                    ndefDIA(ng)*((ntstart(ng)-1)/ndefDIA(ng)))/   &
     &                   nDIA(ng)
          ELSE
            tDIAindx(ng)=(ntstart(ng)-1)/nDIA(ng)
          END IF
        ELSE
          tDIAindx(ng)=tsize
        END IF
      END IF
!
!  Set initial averaged time.
!
      IF (ntsDIA(ng).eq.1) THEN
        DIAtime(ng)=time(ng)-0.5_r8*REAL(nDIA(ng),r8)*dt(ng)
      ELSE
        DIAtime(ng)=time(ng)+REAL(ntsDIA(ng),r8)*dt(ng)-                &
     &              0.5_r8*REAL(nDIA(ng),r8)*dt(ng)
      END IF
!
  10  FORMAT (6x,'DEF_DIAGS - creating diagnostics file: ',a)
  20  FORMAT (6x,'DEF_DIAGS - inquiring diagnostics file: ',a)
  30  FORMAT (/,' DEF_DIAGS - unable to create diagnostics NetCDF',     &
     &        ' file: ',a)
  40  FORMAT (1pe11.4,1x,'millimeter')
  50  FORMAT (/,' DEF_DIAGS - unable to open diagnostics NetCDF',       &
     &        ' file: ',a)
  60  FORMAT (/,' DEF_DIAGS - unable to find variable: ',a,2x,          &
     &        ' in diagnostics NetCDF file: ',a)
      RETURN
      END SUBROUTINE def_diags
#else
      SUBROUTINE def_diags
      RETURN
      END SUBROUTINE def_diags
#endif
