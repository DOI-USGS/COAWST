#include "cppdefs.h"
#if defined IS4DVAR
      SUBROUTINE def_hessian (ng)
!
!svn $Id: def_hessian.F 733 2008-09-07 01:56:45Z jcwarner $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine creates the Lanczos eigenvectors NetCDF file used for  !
!  the pre-conditioning of the 4DVar conjugate gradient algorithm. It  !
!  defines its dimensions, attributes, and variables.                  !
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
# if defined SEDIMENT || defined BBL_MODEL
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
!
!  Local variable declarations.
!
      logical :: got_var(NV)

      integer, parameter :: Natt = 24

      integer :: i, j, itrc, nrec, nvd, nvd3, nvd4
      integer :: recdim, status, varid
# if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
      integer :: frecdim
# endif
      integer :: DimIDs(29), t2dgrd(3), u2dgrd(3), v2dgrd(3)
      integer :: Vsize(4)

      integer :: def_info, def_var, wrt_info

# ifdef SOLVE3D
      integer :: t3dgrd(4), u3dgrd(4), v3dgrd(4), w3dgrd(4)
#  ifdef ADJUST_STFLUX
      integer :: t3dfrc(4)
#  endif
# endif
# ifdef ADJUST_WSTRESS
      integer :: u3dfrc(4), v3dfrc(4)
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
      ncname=HSSname(ng)
!
      IF (Master) THEN
        IF (LdefHSS(ng)) THEN
          WRITE (stdout,10) TRIM(ncname)
        ELSE
          WRITE (stdout,20) TRIM(ncname)
        END IF
      END IF
!
!=======================================================================
!  Create a new Lanczos eigenvectors file.
!=======================================================================
!
      IF (LdefHSS(ng).and.OutThread) THEN
        status=nf90_create(TRIM(ncname),nf90_clobber,ncHSSid(ng))
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,30) TRIM(ncname)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# ifdef DISTRIBUTE
      IF (LdefHSS(ng)) THEN
        CALL mp_bcasti (ng, iADM, ncHSSid(ng), 1)
      END IF
# endif
!
!-----------------------------------------------------------------------
!  Define the dimensions of staggered fields.
!-----------------------------------------------------------------------
!
      IF (LdefHSS(ng).and.OutThread) THEN
        status=nf90_def_dim(ncHSSid(ng),'xi_rho'  ,Lm(ng)+2,            &
     &                    DimIDs( 1))
        status=nf90_def_dim(ncHSSid(ng),'xi_u'    ,Lm(ng)+1,            &
     &                    DimIDs( 2))
        status=nf90_def_dim(ncHSSid(ng),'xi_v'    ,Lm(ng)+2,            &
     &                    DimIDs( 3))
        status=nf90_def_dim(ncHSSid(ng),'xi_psi'  ,Lm(ng)+1,            &
     &                    DimIDs( 4))
        status=nf90_def_dim(ncHSSid(ng),'eta_rho' ,Mm(ng)+2,            &
     &                    DimIDs( 5))
        status=nf90_def_dim(ncHSSid(ng),'eta_u'   ,Mm(ng)+2,            &
     &                    DimIDs( 6))
        status=nf90_def_dim(ncHSSid(ng),'eta_v'   ,Mm(ng)+1,            &
     &                    DimIDs( 7))
        status=nf90_def_dim(ncHSSid(ng),'eta_psi' ,Mm(ng)+1,            &
     &                    DimIDs( 8))
# if defined WRITE_WATER && defined MASKING
        status=nf90_def_dim(ncHSSid(ng),'xy_rho'  ,Nxyr(ng),            &
     &                    DimIDs(17))
        status=nf90_def_dim(ncHSSid(ng),'xy_u'    ,Nxyu(ng),            &
     &                    DimIDs(18))
        status=nf90_def_dim(ncHSSid(ng),'xy_v'    ,Nxyv(ng),            &
     &                    DimIDs(19))
# endif
# ifdef SOLVE3D
#  if defined WRITE_WATER && defined MASKING
        status=nf90_def_dim(ncHSSid(ng),'xyz_rho' ,Nxyr(ng)*N(ng),      &
     &                      DimIDs(20))
        status=nf90_def_dim(ncHSSid(ng),'xyz_u'   ,Nxyu(ng)*N(ng),      &
     &                      DimIDs(21))
        status=nf90_def_dim(ncHSSid(ng),'xyz_v'   ,Nxyv(ng)*N(ng),      &
     &                      DimIDs(22))
        status=nf90_def_dim(ncHSSid(ng),'xyz_w'   ,Nxyr(ng)*(N(ng)+1),  &
     &                      DimIDs(23))
#  endif
        status=nf90_def_dim(ncHSSid(ng),'N'       ,N(ng),               &
     &                      DimIDs( 9))
        status=nf90_def_dim(ncHSSid(ng),'s_rho'   ,N(ng),               &
     &                      DimIDs( 9))
        status=nf90_def_dim(ncHSSid(ng),'s_w'     ,N(ng)+1,             &
     &                      DimIDs(10))
        status=nf90_def_dim(ncHSSid(ng),'tracer'  ,NT(ng),              &
     &                      DimIDs(11))
#  ifdef SEDIMENT_NOT_YET
        status=nf90_def_dim(ncHSSid(ng),'Nbed'    ,Nbed,                &
     &                      DimIDs(16))
#   if defined WRITE_WATER && defined MASKING
        status=nf90_def_dim(ncHSSid(ng),'xybed'   ,Nxyr(ng)*Nbed,       &
     &                      DimIDs(24))
#   endif
#  endif
#  ifdef ECOSIM
        status=nf90_def_dim(ncHSSid(ng),'Nphy'    ,Nphy,                &
     &                      DimIDs(25))
        status=nf90_def_dim(ncHSSid(ng),'Nbac'    ,Nbac,                &
     &                      DimIDs(26))
        status=nf90_def_dim(ncHSSid(ng),'Ndom'    ,Ndom,                &
     &                      DimIDs(27))
        status=nf90_def_dim(ncHSSid(ng),'Nfec'    ,Nfec,                &
     &                      DimIDs(28))
#  endif
# endif
        status=nf90_def_dim(ncHSSid(ng),'boundary',4,                   &
     &                      DimIDs(14))
# ifdef FOUR_DVAR
        status=nf90_def_dim(ncHSSid(ng),'Nstate  ',NstateVar(ng),       &
     &                      DimIDs(29))
# endif
# if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
        status=nf90_def_dim(ncHSSid(ng),'Nadjust ',Nfrec(ng),           &
     &                      frecdim)
# endif
        status=nf90_def_dim(ncHSSid(ng),TRIM(ADJUSTL(Vname(5,idtime))), &
     &                      nf90_unlimited,DimIDs(12))
        recdim=DimIDs(12)
!
!  Set number of dimensions for output variables.
!
# if defined WRITE_WATER && defined MASKING
        nvd3=2
        nvd4=2
# else
        nvd3=3
        nvd4=4
# endif
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
#  ifdef ADJUST_STFLUX
        t3dfrc(1)=DimIDs( 1)
        t3dfrc(2)=DimIDs( 5)
        t3dfrc(3)=frecdim
        t3dfrc(4)=DimIDs(12)
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
#  ifdef ADJUST_WSTRESS
        u3dfrc(1)=DimIDs( 2)
        u3dfrc(2)=DimIDs( 6)
        u3dfrc(3)=frecdim
        u3dfrc(4)=DimIDs(12)
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
#  ifdef ADJUST_WSTRESS
        v3dfrc(1)=DimIDs( 3)
        v3dfrc(2)=DimIDs( 7)
        v3dfrc(3)=frecdim
        v3dfrc(4)=DimIDs(12)
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
# endif
!
!  Initialize unlimited time record dimension.
!
        tHSSindx(ng)=0
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
        status=def_info(ng,ncHSSid(ng),ncname,DimIDs)
        IF (exit_flag.ne.NoError) RETURN
!
!-----------------------------------------------------------------------
!  Define time-varying variables.
!-----------------------------------------------------------------------
!
!  Define number of converged Ritz eigenvalues.
!
        Vinfo( 1)='nConvRitz'
        Vinfo( 2)='number of converged Ritz eigenvalues'
        status=def_var(ncHSSid(ng),varid,nf90_int,1,0,                  &
     &                 Aval,Vinfo,ncname)
!
!  Define Ritz eigenvalues.
!
        Vinfo( 1)='Ritz'
        Vinfo( 2)='Ritz eigenvalues'
        status=def_var(ncHSSid(ng),varid,NF_TYPE,1,recdim,              &
     &                 Aval,Vinfo,ncname) 
!
!  Define accuracy of Ritz eigenvalues .
!
        Vinfo( 1)='Ritz_error'
        Vinfo( 2)='accuracy of Ritz eigenvalues'
        status=def_var(ncHSSid(ng),varid,NF_TYPE,1,recdim,              &
     &                 Aval,Vinfo,ncname)
!
!  Define model time.
!
        Vinfo( 1)=Vname(1,idtime)
        Vinfo( 2)=Vname(2,idtime)
        IF (INT(time_ref).eq.-2) THEN
          Vinfo( 3)='seconds since 1968-05-23 00:00:00 GMT'
          Vinfo( 4)='gregorian'
        ELSE IF (INT(time_ref).eq.-1) THEN
          Vinfo( 3)='seconds since 0001-01-01 00:00:00'
          Vinfo( 4)='360_day'
        ELSE IF (INT(time_ref).eq.0) THEN
          Vinfo( 3)='seconds since 0001-01-01 00:00:00'
          Vinfo( 4)='365.25_day'
        ELSE IF (time_ref.gt.0.0_r8) THEN
          WRITE (Vinfo( 3),'(a,1x,a)') 'seconds since', TRIM(r_text)
          Vinfo( 4)='standard'
        END IF
        Vinfo(14)=Vname(4,idtime)
        status=def_var(ncHSSid(ng),hssVid(idtime,ng),NF_TYPE,           &
     &                 1,recdim,Aval,Vinfo,ncname)
!
!  Define free-surface.
!
        Vinfo( 1)=Vname(1,idFsur)
        WRITE (Vinfo( 2),40) TRIM(Vname(2,idFsur))
        Vinfo( 3)='nondimensional'
        Vinfo(14)=Vname(4,idFsur)
        Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_rho'
# endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idFsur,ng),r8)
        status=def_var(ncHSSid(ng),hssVid(idFsur,ng),NF_FOUT,           &
     &                 nvd3,t2dgrd,Aval,Vinfo,ncname)
!
!  Define 2D momentum in the XI-direction.
!
        Vinfo( 1)=Vname(1,idUbar)
        WRITE (Vinfo( 2),40) TRIM(Vname(2,idUbar))
        Vinfo( 3)='nondimensional'
        Vinfo(14)=Vname(4,idUbar)
        Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_u'
# endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idUbar,ng),r8)
        status=def_var(ncHSSid(ng),hssVid(idUbar,ng),NF_FOUT,           &
     &                 nvd3,u2dgrd,Aval,Vinfo,ncname)
!
!  Define 2D momentum in the ETA-direction.
!
        Vinfo( 1)=Vname(1,idVbar)
        WRITE (Vinfo( 2),40) TRIM(Vname(2,idVbar))
        Vinfo( 3)='nondimensional'
        Vinfo(14)=Vname(4,idVbar)
        Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_v'
# endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idVbar,ng),r8)
        status=def_var(ncHSSid(ng),hssVid(idVbar,ng),NF_FOUT,           &
     &                 nvd3,v2dgrd,Aval,Vinfo,ncname)
# ifdef SOLVE3D
!
!  Define 3D momentum component in the XI-direction.
!
        Vinfo( 1)=Vname(1,idUvel)
        WRITE (Vinfo( 2),40) TRIM(Vname(2,idUvel))
        Vinfo( 3)='nondimensional'
        Vinfo(14)=Vname(4,idUvel)
        Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_u'
#  endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idUvel,ng),r8)
        status=def_var(ncHSSid(ng),hssVid(idUvel,ng),NF_FOUT,           &
     &                 nvd4,u3dgrd,Aval,Vinfo,ncname)
!
!  Define 3D momentum component in the ETA-direction.
!
        Vinfo( 1)=Vname(1,idVvel)
        WRITE (Vinfo( 2),40) TRIM(Vname(2,idVvel))
        Vinfo( 3)='nondimensional'
        Vinfo(14)=Vname(4,idVvel)
        Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_v'
#  endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idVvel,ng),r8)
        status=def_var(ncHSSid(ng),hssVid(idVvel,ng),NF_FOUT,           &
     &                 nvd4,v3dgrd,Aval,Vinfo,ncname)
!
!  Define tracer type variables.
!
        DO itrc=1,NT(ng)
          Vinfo( 1)=Vname(1,idTvar(itrc))
          WRITE (Vinfo( 2),40) TRIM(Vname(2,idTvar(itrc)))
          Vinfo( 3)='nondimensional'
          Vinfo(14)=Vname(4,idTvar(itrc))
          Vinfo(16)=Vname(1,idtime)
#  ifdef SEDIMENT
          DO i=1,NST
            IF (itrc.eq.idsed(i)) THEN
              WRITE (Vinfo(19),50) 1000.0_r8*Sd50(i,ng)
            END IF
          END DO
#  endif
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(r3dvar,r8)
          status=def_var(ncHSSid(ng),hssTid(itrc,ng),NF_FOUT,           &
     &                     nvd4,t3dgrd,Aval,Vinfo,ncname)
        END DO
#  ifdef ADJUST_STFLUX
!
!  Define surface tracer fluxes.
!
        DO itrc=1,NT(ng)
          Vinfo( 1)=Vname(1,idTsur(itrc))
          WRITE (Vinfo( 2),40) TRIM(Vname(2,idTsur(itrc)))
          Vinfo( 3)='nondimensional'
          IF (itrc.eq.itemp) THEN
            Vinfo(11)='upward flux, cooling'
            Vinfo(12)='downward flux, heating'
          ELSE IF (itrc.eq.isalt) THEN
            Vinfo(11)='upward flux, freshening (net precipitation)'
            Vinfo(12)='downward flux, salting (net evaporation)'
          END IF
          Vinfo(14)=Vname(4,idTsur(itrc))
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(r2dvar,r8)
          status=def_var(ncHSSid(ng),hssVid(idTsur(itrc),ng),           &
     &                   NF_FOUT,nvd4,t3dfrc,Aval,Vinfo,ncname)
        END DO
#  endif
# endif
# ifdef ADJUST_WSTRESS
!
!  Define surface U-momentum stress.
!
        Vinfo( 1)=Vname(1,idUsms)
        WRITE (Vinfo( 2),40) TRIM(Vname(2,idUsms))
        Vinfo( 3)='nondimensional'
        Vinfo(14)=Vname(4,idUsms)
        Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_u'
#  endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(u2dvar,r8)
        status=def_var(ncHSSid(ng),hssVid(idUsms,ng),NF_FOUT,           &
     &                 nvd4,u3dfrc,Aval,Vinfo,ncname)
!
!  Define surface V-momentum stress.
!
        Vinfo( 1)=Vname(1,idVsms)
        WRITE (Vinfo( 2),40) TRIM(Vname(2,idVsms))
        Vinfo( 2)=Vname(2,idVsms)
        Vinfo( 3)='nondimensional'
        Vinfo(14)=Vname(4,idVsms)
        Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_v'
#  endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(v2dvar,r8)
        status=def_var(ncHSSid(ng),hssVid(idVsms,ng),NF_FOUT,           &
     &                 nvd4,v3dfrc,Aval,Vinfo,ncname)
# endif
!
!-----------------------------------------------------------------------
!  Leave definition mode.
!-----------------------------------------------------------------------
!
        status=nf90_enddef(ncHSSid(ng))
      END IF
!
!-----------------------------------------------------------------------
!  Write out time-recordless, information variables.
!-----------------------------------------------------------------------
!
      IF (LdefHSS(ng)) THEN
        status=wrt_info(ng, iADM, ncHSSid(ng), OutThread, ncname)
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!=======================================================================
!  Open an existing Hessian eigenvectors file, check its contents, and
!  prepare for appending data.
!=======================================================================
!
      IF (.not.LdefHSS(ng).and.OutThread) THEN
!
!  Inquire about the contents of Hessian NetCDF file:  Inquire about
!  the dimensions and variables.  Check for consistency.
!
        CALL opencdf (ng, 1, ncname, fname, N(ng), 0, nrec, nvd, Vsize)
        IF (exit_flag.ne.NoError) RETURN
!
!  Open Lanczos eigenvectors file for read/write.
!
        status=nf90_open(TRIM(ncname),nf90_write,ncHSSid(ng))
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,60) TRIM(ncname)
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
!  Lanczos eigenvectors variables. Get variable IDs.
!
        DO i=1,nvars
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idtime))) THEN
            got_var(idtime)=.TRUE.
            status=nf90_inq_varid(ncHSSid(ng),TRIM(Vname(1,idtime)),    &
     &                            hssVid(idtime,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idFsur))) THEN
            got_var(idFsur)=.TRUE.
            status=nf90_inq_varid(ncHSSid(ng),TRIM(Vname(1,idFsur)),    &
     &                            hssVid(idFsur,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idUbar))) THEN
            got_var(idUbar)=.TRUE.
            status=nf90_inq_varid(ncHSSid(ng),TRIM(Vname(1,idUbar)),    &
     &                            hssVid(idUbar,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idVbar))) THEN
            got_var(idVbar)=.TRUE.
            status=nf90_inq_varid(ncHSSid(ng),TRIM(Vname(1,idVbar)),    &
     &                            hssVid(idVbar,ng))
          END IF
# ifdef SOLVE3D
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idUvel))) THEN
            got_var(idUvel)=.TRUE.
            status=nf90_inq_varid(ncHSSid(ng),TRIM(Vname(1,idUvel)),    &
     &                            hssVid(idUvel,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idVvel))) THEN
            got_var(idVvel)=.TRUE.
            status=nf90_inq_varid(ncHSSid(ng),TRIM(Vname(1,idVvel)),    &
     &                            hssVid(idVvel,ng))
          END IF
          DO itrc=1,NT(ng)
            IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idTvar(itrc)))) THEN
              got_var(idTvar(itrc))=.TRUE.
              status=nf90_inq_varid(ncHSSid(ng),                        &
     &                              TRIM(Vname(1,idTvar(itrc))),        &
     &                              hssTid(itrc,ng))
            END IF
          END DO
#  ifdef ADJUST_STFLUX
          DO itrc=1,NT(ng)
            IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idTsur(itrc)))) THEN
              got_var(idTsur(itrc))=.TRUE.
              status=nf90_inq_varid(ncHSSid(ng),                        &
     &                              TRIM(Vname(1,idTsur(itrc))),        &
     &                              hssVid(idTsur(itrc),ng))
            END IF
          END DO
#  endif
# endif
# ifdef ADJUST_WSTRESS
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idUsms))) THEN
            got_var(idUsms)=.TRUE.
            status=nf90_inq_varid(ncHSSid(ng),TRIM(Vname(1,idUsms)),    &
     &                            hssVid(idUsms,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idVsms))) THEN
            got_var(idVsms)=.TRUE.
            status=nf90_inq_varid(ncHSSid(ng),TRIM(Vname(1,idVsms)),    &
     &                            hssVid(idVsms,ng))
          END IF
# endif
        END DO
!
!  Check if Lanczos eigenvectors variables are available in input
!  NetCDF file.
!
        IF (.not.got_var(idtime)) THEN
          WRITE (stdout,70) TRIM(Vname(1,idtime)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idFsur)) THEN
          WRITE (stdout,70) TRIM(Vname(1,idFsur)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUbar)) THEN
          WRITE (stdout,70) TRIM(Vname(1,idUbar)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVbar)) THEN
          WRITE (stdout,70) TRIM(Vname(1,idVbar)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# ifdef SOLVE3D
        IF (.not.got_var(idUvel)) THEN
          WRITE (stdout,70) TRIM(Vname(1,idUvel)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVvel)) THEN
          WRITE (stdout,70) TRIM(Vname(1,idVvel)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        DO itrc=1,NT(ng)
          IF (.not.got_var(idTvar(itrc))) THEN
            WRITE (stdout,70) TRIM(Vname(1,idTvar(itrc))), TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
#  ifdef ADJUST_STFLUX
        DO itrc=1,NT(ng)
          IF (.not.got_var(idTsur(itrc))) THEN
            WRITE (stdout,70) TRIM(Vname(1,idTsur(itrc))), TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
#  endif
# endif
# ifdef ADJUST_WSTRESS
        IF (.not.got_var(idUsms)) THEN
          WRITE (stdout,70) TRIM(Vname(1,idUsms)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVsms)) THEN
          WRITE (stdout,70) TRIM(Vname(1,idVsms)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
!
!  Set unlimited time record dimension to the appropriate value.
!
        tHSSindx(ng)=tsize
      END IF
!
  10  FORMAT (3x,'DEF_HESSIAN  - creating Hessian  file: ',a)
  20  FORMAT (3x,'DEF_HESSIAN  - inquiring Hessian file: ',a)
  30  FORMAT (/,' DEF_HESSIAN - unable to create Hessian NetCDF file:', &
     &        1x,a)
  40  FORMAT (a,', Hessian eigenvectors')
  50  FORMAT (1pe11.4,1x,'millimeter')
  60  FORMAT (/,' DEF_HESSIAN - unable to open Hessian NetCDF file: ',a)
  70  FORMAT (/,' DEF_HESSIAN - unable to find variable: ',a,2x,        &
     &        ' in Hessian NetCDF file: ',a)
      RETURN
      END SUBROUTINE def_hessian
#else
      SUBROUTINE def_hessian
      RETURN
      END SUBROUTINE def_hessian
#endif
