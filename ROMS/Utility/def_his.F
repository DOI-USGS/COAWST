#include "cppdefs.h"
      SUBROUTINE def_his (ng,ldef)
!
!svn $Id: def_his.F 733 2008-09-07 01:56:45Z jcwarner $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine creates history NetCDF file, it defines its            !
!  dimensions, attributes, and variables.                              !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
#ifdef FOUR_DVAR
      USE mod_fourdvar
#endif
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars
#if defined SEDIMENT || defined BBL_MODEL
      USE mod_sediment
#endif
#ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_bcasti
#endif
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
      logical :: got_var(NV)

      integer, parameter :: Natt = 24

      integer :: i, j, itrc, nrec, nvd, nvd3, nvd4, varid
      integer :: recdim, status, ofst
#if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
      integer :: frecdim
#endif
      integer :: DimIDs(29), t2dgrd(3), u2dgrd(3), v2dgrd(3)
      integer :: Vsize(4)

      integer :: def_info, def_var, wrt_info

#ifdef SOLVE3D
# ifdef SEDIMENT
      integer :: b3dgrd(4)
# endif
      integer :: t3dgrd(4), u3dgrd(4), v3dgrd(4), w3dgrd(4)
# ifdef ADJUST_STFLUX
      integer :: t3dfrc(3)
# endif
#endif
#ifdef ADJUST_WSTRESS
      integer :: u3dfrc(3), v3dfrc(3)
#endif

      real(r8) :: Aval(6)

      character (len=80) :: Vinfo(Natt)
      character (len=80) :: fname, ncname
!
!-----------------------------------------------------------------------
!  Set and report file name.
!-----------------------------------------------------------------------
!
      IF (exit_flag.ne.NoError) RETURN
      ncname=HISname(ng)
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
!  Create a new history file.
!=======================================================================
!
      IF (ldef.and.OutThread) THEN
        status=nf90_create(TRIM(ncname),nf90_clobber,ncHISid(ng))
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,30) TRIM(ncname)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#ifdef DISTRIBUTE
      IF (ldef) THEN
        CALL mp_bcasti (ng, iNLM, ncHISid(ng), 1)
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
        status=nf90_def_dim(ncHISid(ng),'xi_rho'  ,Lm(ng)+ofst+1,       &
     &                      DimIDs( 1))
        status=nf90_def_dim(ncHISid(ng),'xi_u'    ,Lm(ng)+ofst,         &
     &                      DimIDs( 2))
        status=nf90_def_dim(ncHISid(ng),'xi_v'    ,Lm(ng)+ofst+1,       &
     &                      DimIDs( 3))
        status=nf90_def_dim(ncHISid(ng),'xi_psi'  ,Lm(ng)+ofst,         &
     &                      DimIDs( 4))
        status=nf90_def_dim(ncHISid(ng),'eta_rho' ,Mm(ng)+ofst+1,       &
     &                      DimIDs( 5))
        status=nf90_def_dim(ncHISid(ng),'eta_u'   ,Mm(ng)+ofst+1,       &
     &                      DimIDs( 6))
        status=nf90_def_dim(ncHISid(ng),'eta_v'   ,Mm(ng)+ofst,         &
     &                      DimIDs( 7))
        status=nf90_def_dim(ncHISid(ng),'eta_psi' ,Mm(ng)+ofst,         &
     &                      DimIDs( 8))
#if defined WRITE_WATER && defined MASKING
        status=nf90_def_dim(ncHISid(ng),'xy_rho'  ,Nxyr(ng),            &
     &                      DimIDs(17))
        status=nf90_def_dim(ncHISid(ng),'xy_u'    ,Nxyu(ng),            &
     &                      DimIDs(18))
        status=nf90_def_dim(ncHISid(ng),'xy_v'    ,Nxyv(ng),            &
     &                      DimIDs(19))
#endif
#ifdef SOLVE3D
# if defined WRITE_WATER && defined MASKING
        status=nf90_def_dim(ncHISid(ng),'xyz_rho' ,Nxyr(ng)*N(ng),      &
     &                      DimIDs(20))
        status=nf90_def_dim(ncHISid(ng),'xyz_u'   ,Nxyu(ng)*N(ng),      &
     &                      DimIDs(21))
        status=nf90_def_dim(ncHISid(ng),'xyz_v'   ,Nxyv(ng)*N(ng),      &
     &                      DimIDs(22))
        status=nf90_def_dim(ncHISid(ng),'xyz_w'   ,Nxyr(ng)*(N(ng)+1),  &
     &                      DimIDs(23))
# endif
        status=nf90_def_dim(ncHISid(ng),'N'       ,N(ng),               &
     &                      DimIDs( 9))
        status=nf90_def_dim(ncHISid(ng),'s_rho'   ,N(ng),               &
     &                      DimIDs( 9))
        status=nf90_def_dim(ncHISid(ng),'s_w'     ,N(ng)+1,             &
     &                      DimIDs(10))
        status=nf90_def_dim(ncHISid(ng),'tracer'  ,NT(ng),              &
     &                      DimIDs(11))
# ifdef SEDIMENT
        status=nf90_def_dim(ncHISid(ng),'Nbed'    ,Nbed,                &
     &                      DimIDs(16))
#  if defined WRITE_WATER && defined MASKING
        status=nf90_def_dim(ncHISid(ng),'xybed'   ,Nxyr(ng)*Nbed,       &
     &                      DimIDs(24))
#  endif
# endif
# ifdef ECOSIM
        status=nf90_def_dim(ncHISid(ng),'Nphy'    ,Nphy,                &
     &                      DimIDs(25))
        status=nf90_def_dim(ncHISid(ng),'Nbac'    ,Nbac,                &
     &                      DimIDs(26))
        status=nf90_def_dim(ncHISid(ng),'Ndom'    ,Ndom,                &
     &                      DimIDs(27))
        status=nf90_def_dim(ncHISid(ng),'Nfec'    ,Nfec,                &
     &                      DimIDs(28))
# endif
#endif
        status=nf90_def_dim(ncHISid(ng),'boundary',4,                   &
     &                      DimIDs(14))
#ifdef FOUR_DVAR
        status=nf90_def_dim(ncHISid(ng),'Nstate'  ,NstateVar(ng),       &
     &                      DimIDs(29))
#endif
#if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
        status=nf90_def_dim(ncHISid(ng),'Nadjust ',Nfrec(ng),           &
     &                      frecdim)
#endif
        status=nf90_def_dim(ncHISid(ng),TRIM(ADJUSTL(Vname(5,idtime))), &
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
#if defined WRITE_WATER && defined MASKING
        t2dgrd(1)=DimIDs(17)
        t2dgrd(2)=DimIDs(12)
# ifdef SOLVE3D
        t3dgrd(1)=DimIDs(20)
        t3dgrd(2)=DimIDs(12)
# endif
#else
        t2dgrd(1)=DimIDs( 1)
        t2dgrd(2)=DimIDs( 5)
        t2dgrd(3)=DimIDs(12)
# ifdef SOLVE3D
        t3dgrd(1)=DimIDs( 1)
        t3dgrd(2)=DimIDs( 5)
        t3dgrd(3)=DimIDs( 9)
        t3dgrd(4)=DimIDs(12)
# endif
# ifdef ADJUST_STFLUX
        t3dfrc(1)=DimIDs( 1)
        t3dfrc(2)=DimIDs( 5)
        t3dfrc(3)=frecdim
# endif
#endif
!
!  Define dimension vectors for staggered u-momemtum type variables.
!
#if defined WRITE_WATER && defined MASKING
        u2dgrd(1)=DimIDs(18)
        u2dgrd(2)=DimIDs(12)
# ifdef SOLVE3D
        u3dgrd(1)=DimIDs(21)
        u3dgrd(2)=DimIDs(12)
# endif
#else
        u2dgrd(1)=DimIDs( 2)
        u2dgrd(2)=DimIDs( 6)
        u2dgrd(3)=DimIDs(12)
# ifdef SOLVE3D
        u3dgrd(1)=DimIDs( 2)
        u3dgrd(2)=DimIDs( 6)
        u3dgrd(3)=DimIDs( 9)
        u3dgrd(4)=DimIDs(12)
# endif
# ifdef ADJUST_WSTRESS
        u3dfrc(1)=DimIDs( 2)
        u3dfrc(2)=DimIDs( 6)
        u3dfrc(3)=frecdim
# endif
#endif
!
!  Define dimension vectors for staggered v-momemtum type variables.
!
#if defined WRITE_WATER && defined MASKING
        v2dgrd(1)=DimIDs(19)
        v2dgrd(2)=DimIDs(12)
# ifdef SOLVE3D
        v3dgrd(1)=DimIDs(22)
        v3dgrd(2)=DimIDs(12)
# endif
#else
        v2dgrd(1)=DimIDs( 3)
        v2dgrd(2)=DimIDs( 7)
        v2dgrd(3)=DimIDs(12)
# ifdef SOLVE3D
        v3dgrd(1)=DimIDs( 3)
        v3dgrd(2)=DimIDs( 7)
        v3dgrd(3)=DimIDs( 9)
        v3dgrd(4)=DimIDs(12)
# endif
# ifdef ADJUST_WSTRESS
        v3dfrc(1)=DimIDs( 3)
        v3dfrc(2)=DimIDs( 7)
        v3dfrc(3)=frecdim
# endif
#endif
#ifdef SOLVE3D
!
!  Define dimension vector for staggered w-momemtum type variables.
!
# if defined WRITE_WATER && defined MASKING
        w3dgrd(1)=DimIDs(23)
        w3dgrd(2)=DimIDs(12)
# else
        w3dgrd(1)=DimIDs( 1)
        w3dgrd(2)=DimIDs( 5)
        w3dgrd(3)=DimIDs(10)
        w3dgrd(4)=DimIDs(12)
# endif
# ifdef SEDIMENT
!
!  Define dimension vector for sediment bed layer type variables.
!
#  if defined WRITE_WATER && defined MASKING
        b3dgrd(1)=DimIDs(24)
        b3dgrd(2)=DimIDs(12)
#  else
        b3dgrd(1)=DimIDs( 1)
        b3dgrd(2)=DimIDs( 5)
        b3dgrd(3)=DimIDs(16)
        b3dgrd(4)=DimIDs(12)
#  endif
# endif
#endif
!
!  Initialize unlimited time record dimension.
!
        tHISindx(ng)=0
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
        status=def_info(ng,ncHISid(ng),ncname,DimIDs)
        IF (exit_flag.ne.NoError) RETURN
!
!-----------------------------------------------------------------------
!  Define time-varying variables.
!-----------------------------------------------------------------------
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
          Vinfo( 4)='julian'
        ELSE IF (time_ref.gt.0.0_r8) THEN
          WRITE (Vinfo( 3),'(a,1x,a)') 'seconds since', TRIM(r_text)
          Vinfo( 4)='gregorian'
        END IF
        Vinfo(14)=Vname(4,idtime)
        status=def_var(ncHISid(ng),hisVid(idtime,ng),NF_TYPE,           &
     &                 1,recdim,Aval,Vinfo,ncname)

#ifdef SO_SEMI
!
!  Define Ritz eigenvalues and Ritz eigenvectors Euclidean norm.
!
        Vinfo( 1)='Ritz_rvalue'
        Vinfo( 2)='real Ritz eigenvalues'
        status=def_var(ncHISid(ng),varid,NF_TYPE,                       &
     &                 1,recdim,Aval,Vinfo,ncname)

        Vinfo( 1)='Ritz_norm'
        Vinfo( 2)='Ritz eigenvectors Euclidean norm'
        status=def_var(ncHISid(ng),varid,NF_TYPE,                       &
     &                 1,recdim,Aval,Vinfo,ncname)
#endif
#ifdef ADJUST_WSTRESS
!
!  Define surface U-momentum stress.  Notice that the stress has its
!  own fixed time-dimension (of size Nfrec) to allow 4DVAR adjustments
!  at other times in addition to initialization time.
!
        Vinfo( 1)=Vname(1,idUsms)
        Vinfo( 2)=Vname(2,idUsms)
        Vinfo( 3)=Vname(3,idUsms)
# if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_u'
# endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idUsms,ng),r8)
        status=def_var(ncHISid(ng),hisVid(idUsms,ng),NF_FOUT,           &
     &                 nvd3,u3dfrc,Aval,Vinfo,ncname)
!
!  Define surface V-momentum stress.
!
        Vinfo( 1)=Vname(1,idVsms)
        Vinfo( 2)=Vname(2,idVsms)
        Vinfo( 3)=Vname(3,idVsms)
# if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_v'
# endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idVsms,ng),r8)
        status=def_var(ncHISid(ng),hisVid(idVsms,ng),NF_FOUT,           &
     &                 nvd3,v3dfrc,Aval,Vinfo,ncname)
#endif
#if defined ADJUST_STFLUX && defined SOLVE3D
!
!  Define surface net heat flux. Notice that different tracer fluxes
!  are written at their own fixed time-dimension (of size Nfrec) to
!  allow 4DVAR adjustments at other times in addition to initial time.
!
        DO itrc=1,NT(ng)
          Vinfo( 1)=Vname(1,idTsur(itrc))
          Vinfo( 2)=Vname(2,idTsur(itrc))
          Vinfo( 3)=Vname(3,idTsur(itrc))
          IF (itrc.eq.itemp) THEN
            Vinfo(11)='upward flux, cooling'
            Vinfo(12)='downward flux, heating'
          ELSE IF (itrc.eq.isalt) THEN
            Vinfo(11)='upward flux, freshening (net precipitation)'
            Vinfo(12)='downward flux, salting (net evaporation)'
          END IF
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
# endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idTsur(itrc),ng),r8)
          status=def_var(ncHISid(ng),hisVid(idTsur(itrc),ng),           &
     &                   NF_FOUT,nvd3,t3dfrc,Aval,Vinfo,ncname)
        END DO
#endif
#if defined SEDIMENT && defined SED_MORPH
!
!  Define bathymetry.
!
        IF (Hout(idbath,ng)) THEN
          Vinfo( 1)=Vname(1,idbath)
          Vinfo( 2)=Vname(2,idbath)
          Vinfo( 3)=Vname(3,idbath)
          Vinfo(14)=Vname(4,idbath)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
# endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idbath,ng),r8)
          status=def_var(ncHISid(ng),hisVid(idbath,ng),NF_FOUT,         &
     &                   nvd3,t2dgrd,Aval,Vinfo,ncname)
        END IF
#endif
!
!  Define free-surface.
!
        IF (Hout(idFsur,ng)) THEN
          Vinfo( 1)=Vname(1,idFsur)
          Vinfo( 2)=Vname(2,idFsur)
          Vinfo( 3)=Vname(3,idFsur)
          Vinfo(14)=Vname(4,idFsur)
          Vinfo(16)=Vname(1,idtime)
#if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idFsur,ng),r8)
          status=def_var(ncHISid(ng),hisVid(idFsur,ng),NF_FOUT,         &
     &                   nvd3,t2dgrd,Aval,Vinfo,ncname)
#if defined FORWARD_WRITE && defined FORWARD_RHS
          Vinfo( 1)=Vname(1,idRzet)
          Vinfo( 2)=Vname(2,idRzet)
          Vinfo( 3)=Vname(3,idRzet)
          Vinfo(14)=Vname(4,idRzet)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
# endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(r2dvar,r8)
          status=def_var(ncHISid(ng),hisVid(idRzet,ng),NF_FOUT,         &
     &                   nvd3,t2dgrd,Aval,Vinfo,ncname)
#endif
        END IF
!
!  Define 2D momentum in the XI-direction.
!
        IF (Hout(idUbar,ng)) THEN
          Vinfo( 1)=Vname(1,idUbar)
          Vinfo( 2)=Vname(2,idUbar)
          Vinfo( 3)=Vname(3,idUbar)
          Vinfo(14)=Vname(4,idUbar)
          Vinfo(16)=Vname(1,idtime)
#if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
#endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUbar,ng),r8)
          status=def_var(ncHISid(ng),hisVid(idUbar,ng),NF_FOUT,         &
     &                   nvd3,u2dgrd,Aval,Vinfo,ncname)
#ifdef FORWARD_WRITE
# ifdef FORWARD_RHS
          Vinfo( 1)=Vname(1,idRu2d)
          Vinfo( 2)=Vname(2,idRu2d)
          Vinfo( 3)=Vname(3,idRu2d)
          Vinfo(14)=Vname(4,idRu2d)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(u2dvar,r8)

          status=def_var(ncHISid(ng),hisVid(idRu2d,ng),NF_FOUT,         &
     &                   nvd3,u2dgrd,Aval,Vinfo,ncname)
# endif
# ifdef SOLVE3D
#  ifdef FORWARD_RHS
          Vinfo( 1)=Vname(1,idRuct)
          Vinfo( 2)=Vname(2,idRuct)
          Vinfo( 3)=Vname(3,idRuct)
          Vinfo(14)=Vname(4,idRuct)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
#   endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(u2dvar,r8)
          status=def_var(ncHISid(ng),hisVid(idRuct,ng),NF_FOUT,         &
     &                   nvd3,u2dgrd,Aval,Vinfo,ncname)
#  endif
          Vinfo( 1)=Vname(1,idUfx1)
          Vinfo( 2)=Vname(2,idUfx1)
          Vinfo( 3)=Vname(3,idUfx1)
          Vinfo(14)=Vname(4,idUfx1)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(u2dvar,r8)
          status=def_var(ncHISid(ng),hisVid(idUfx1,ng),NF_FOUT,         &
     &                   nvd3,u2dgrd,Aval,Vinfo,ncname)
          Vinfo( 1)=Vname(1,idUfx2)
          Vinfo( 2)=Vname(2,idUfx2)
          Vinfo( 3)=Vname(3,idUfx2)
          Vinfo(14)=Vname(4,idUfx2)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(u2dvar,r8)
          status=def_var(ncHISid(ng),hisVid(idUfx2,ng),NF_FOUT,         &
     &                   nvd3,u2dgrd,Aval,Vinfo,ncname)
# endif
#endif
        END IF
!
!  Define 2D momentum in the ETA-direction.
!
        IF (Hout(idVbar,ng)) THEN
          Vinfo( 1)=Vname(1,idVbar)
          Vinfo( 2)=Vname(2,idVbar)
          Vinfo( 3)=Vname(3,idVbar)
          Vinfo(14)=Vname(4,idVbar)
          Vinfo(16)=Vname(1,idtime)
#if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
#endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVbar,ng),r8)
          status=def_var(ncHISid(ng),hisVid(idVbar,ng),NF_FOUT,         &
     &                   nvd3,v2dgrd,Aval,Vinfo,ncname)
#ifdef FORWARD_WRITE
# ifdef FORWARD_RHS
          Vinfo( 1)=Vname(1,idRv2d)
          Vinfo( 2)=Vname(2,idRv2d)
          Vinfo( 3)=Vname(3,idRv2d)
          Vinfo(14)=Vname(4,idRv2d)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(v2dvar,r8)
          status=def_var(ncHISid(ng),hisVid(idRv2d,ng),NF_FOUT,         &
     &                   nvd3,v2dgrd,Aval,Vinfo,ncname)
# endif
# ifdef SOLVE3D
#  ifdef FORWARD_RHS
          Vinfo( 1)=Vname(1,idRvct)
          Vinfo( 2)=Vname(2,idRvct)
          Vinfo( 3)=Vname(3,idRvct)
          Vinfo(14)=Vname(4,idRvct)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
#   endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(v2dvar,r8)
          status=def_var(ncHISid(ng),hisVid(idRvct,ng),NF_FOUT,         &
     &                   nvd3,v2dgrd,Aval,Vinfo,ncname)
#  endif
          Vinfo( 1)=Vname(1,idVfx1)
          Vinfo( 2)=Vname(2,idVfx1)
          Vinfo( 3)=Vname(3,idVfx1)
          Vinfo(14)=Vname(4,idVfx1)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(v2dvar,r8)
          status=def_var(ncHISid(ng),hisVid(idVfx1,ng),NF_FOUT,         &
     &                   nvd3,v2dgrd,Aval,Vinfo,ncname)
          Vinfo( 1)=Vname(1,idVfx2)
          Vinfo( 2)=Vname(2,idVfx2)
          Vinfo( 3)=Vname(3,idVfx2)
          Vinfo(14)=Vname(4,idVfx2)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(v2dvar,r8)
          status=def_var(ncHISid(ng),hisVid(idVfx2,ng),NF_FOUT,         &
     &                   nvd3,v2dgrd,Aval,Vinfo,ncname)
# endif
#endif
        END IF
#ifdef SOLVE3D
!
!  Define 3D momentum component in the XI-direction.
!
        IF (Hout(idUvel,ng)) THEN
          Vinfo( 1)=Vname(1,idUvel)
          Vinfo( 2)=Vname(2,idUvel)
          Vinfo( 3)=Vname(3,idUvel)
          Vinfo(14)=Vname(4,idUvel)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
# endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUvel,ng),r8)
          status=def_var(ncHISid(ng),hisVid(idUvel,ng),NF_FOUT,         &
     &                   nvd4,u3dgrd,Aval,Vinfo,ncname)
# if defined FORWARD_WRITE && defined FORWARD_RHS
          Vinfo( 1)=Vname(1,idRu3d)
          Vinfo( 2)=Vname(2,idRu3d)
          Vinfo( 3)=Vname(3,idRu3d)
          Vinfo(14)=Vname(4,idRu3d)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(u3dvar,r8)
          status=def_var(ncHISid(ng),hisVid(idRu3d,ng),NF_FOUT,         &
     &                   nvd4,u3dgrd,Aval,Vinfo,ncname)
# endif
        END IF
!
!  Define 3D momentum component in the ETA-direction.
!
        IF (Hout(idVvel,ng)) THEN
          Vinfo( 1)=Vname(1,idVvel)
          Vinfo( 2)=Vname(2,idVvel)
          Vinfo( 3)=Vname(3,idVvel)
          Vinfo(14)=Vname(4,idVvel)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
# endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVvel,ng),r8)
          status=def_var(ncHISid(ng),hisVid(idVvel,ng),NF_FOUT,         &
     &                   nvd4,v3dgrd,Aval,Vinfo,ncname)
# if defined FORWARD_WRITE && defined FORWARD_RHS
          Vinfo( 1)=Vname(1,idRv3d)
          Vinfo( 2)=Vname(2,idRv3d)
          Vinfo( 3)=Vname(3,idRv3d)
          Vinfo(14)=Vname(4,idRv3d)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(v3dvar,r8)
          status=def_var(ncHISid(ng),hisVid(idRv3d,ng),NF_FOUT,         &
     &                   nvd4,v3dgrd,Aval,Vinfo,ncname)
# endif
        END IF
!
!  Define 3D momentum component in the S-direction.
!
        IF (Hout(idWvel,ng)) THEN
          Vinfo( 1)=Vname(1,idWvel)
          Vinfo( 2)=Vname(2,idWvel)
          Vinfo( 3)=Vname(3,idWvel)
          Vinfo(14)=Vname(4,idWvel)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
# endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idWvel,ng),r8)
          status=def_var(ncHISid(ng),hisVid(idWvel,ng),NF_FOUT,         &
     &                   nvd4,w3dgrd,Aval,Vinfo,ncname)
        END IF
!
!  Define S-coordinate vertical "omega" momentum component.
!
        IF (Hout(idOvel,ng)) THEN
          Vinfo( 1)=Vname(1,idOvel)
          Vinfo( 2)=Vname(2,idOvel)
          Vinfo( 3)=Vname(3,idOvel)
          Vinfo(14)=Vname(4,idOvel)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
# endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idOvel,ng),r8)
          status=def_var(ncHISid(ng),hisVid(idOvel,ng),NF_FOUT,         &
     &                   nvd4,w3dgrd,Aval,Vinfo,ncname)
        END IF
!
!  Define tracer type variables.
!
        DO itrc=1,NT(ng)
          IF (Hout(idTvar(itrc),ng)) THEN
            Vinfo( 1)=Vname(1,idTvar(itrc))
            Vinfo( 2)=Vname(2,idTvar(itrc))
            Vinfo( 3)=Vname(3,idTvar(itrc))
            Vinfo(14)=Vname(4,idTvar(itrc))
            Vinfo(16)=Vname(1,idtime)
# ifdef SEDIMENT
            DO i=1,NST
              IF (itrc.eq.idsed(i)) THEN
                WRITE (Vinfo(19),40) 1000.0_r8*Sd50(i,ng)
              END IF
            END DO
# endif
# if defined WRITE_WATER && defined MASKING
            Vinfo(20)='mask_rho'
# endif
            Vinfo(22)='coordinates'
            Aval(5)=REAL(r3dvar,r8)
            status=def_var(ncHISid(ng),hisTid(itrc,ng),NF_FOUT,         &
     &                     nvd4,t3dgrd,Aval,Vinfo,ncname)
          END IF
        END DO
!
!  Define density anomaly.
!
        IF (Hout(idDano,ng)) THEN
          Vinfo( 1)=Vname(1,idDano)
          Vinfo( 2)=Vname(2,idDano)
          Vinfo( 3)=Vname(3,idDano)
          Vinfo(14)=Vname(4,idDano)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
# endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idDano,ng),r8)
          status=def_var(ncHISid(ng),hisVid(idDano,ng),NF_FOUT,         &
     &                   nvd4,t3dgrd,Aval,Vinfo,ncname)
        END IF
# ifdef LMD_SKPP
!
!  Define depth of surface boundary layer.
!
        IF (Hout(idHsbl,ng)) THEN
          Vinfo( 1)=Vname(1,idHsbl)
          Vinfo( 2)=Vname(2,idHsbl)
          Vinfo( 3)=Vname(3,idHsbl)
          Vinfo(14)=Vname(4,idHsbl)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idHsbl,ng),r8)
          status=def_var(ncHISid(ng),hisVid(idHsbl,ng),NF_FOUT,         &
     &                   nvd3,t2dgrd,Aval,Vinfo,ncname)
        END IF
# endif
# ifdef LMD_BKPP
!
!  Define depth of bottom boundary layer.
!
        IF (Hout(idHbbl,ng)) THEN
          Vinfo( 1)=Vname(1,idHbbl)
          Vinfo( 2)=Vname(2,idHbbl)
          Vinfo( 3)=Vname(3,idHbbl)
          Vinfo(14)=Vname(4,idHbbl)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idHbbl,ng),r8)
          status=def_var(ncHISid(ng),hisVid(idHbbl,ng),NF_FOUT,         &
     &                   nvd3,t2dgrd,Aval,Vinfo,ncname)
        END IF
# endif
# if defined FORWARD_WRITE && defined LMD_NONLOCAL
!
!  Define out KPP nonlocal transport.
!
        DO itrc=1,NAT
          IF (Hout(idGhat(itrc),ng)) THEN
            Vinfo( 1)=Vname(1,idGhat(itrc))
            Vinfo( 2)=Vname(2,idGhat(itrc))
            Vinfo( 3)=Vname(3,idGhat(itrc))
            Vinfo(14)=Vname(4,idGhat(itrc))
            Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
            Vinfo(20)='mask_rho'
#  endif
            Vinfo(22)='coordinates'
            Aval(5)=REAL(Iinfo(1,idGhat(itrc),ng),r8)
            status=def_var(ncHISid(ng),hisVid(idGhat(itrc),ng),NF_FOUT, &
     &                     nvd4,w3dgrd,Aval,Vinfo,ncname)
          END IF
        END DO
# endif
!
!  Define vertical viscosity coefficient.
!
        IF (Hout(idVvis,ng)) THEN
          Vinfo( 1)=Vname(1,idVvis)
          Vinfo( 2)=Vname(2,idVvis)
          Vinfo( 3)=Vname(3,idVvis)
          Vinfo(14)=Vname(4,idVvis)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
# endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVvis,ng),r8)
          status=def_var(ncHISid(ng),hisVid(idVvis,ng),NF_FOUT,         &
     &                   nvd4,w3dgrd,Aval,Vinfo,ncname)
        END IF
!
!  Define vertical diffusion coefficient for potential temperature.
!
        IF (Hout(idTdif,ng)) THEN
          Vinfo( 1)=Vname(1,idTdif)
          Vinfo( 2)=Vname(2,idTdif)
          Vinfo( 3)=Vname(3,idTdif)
          Vinfo(14)=Vname(4,idTdif)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
# endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idTdif,ng),r8)
          status=def_var(ncHISid(ng),hisVid(idTdif,ng),NF_FOUT,         &
     &                   nvd4,w3dgrd,Aval,Vinfo,ncname)
        END IF
# ifdef SALINITY
!
!  Define vertical diffusion coefficient for salinity.
!
        IF (Hout(idSdif,ng)) THEN
          Vinfo( 1)=Vname(1,idSdif)
          Vinfo( 2)=Vname(2,idSdif)
          Vinfo( 3)=Vname(3,idSdif)
          Vinfo(14)=Vname(4,idSdif)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idSdif,ng),r8)
          status=def_var(ncHISid(ng),hisVid(idSdif,ng),NF_FOUT,         &
     &                   nvd4,w3dgrd,Aval,Vinfo,ncname)
        END IF
# endif
# if defined GLS_MIXING || defined MY25_MIXING
!
!  Define turbulent kinetic energy.
!
        IF (Hout(idMtke,ng)) THEN
          Vinfo( 1)=Vname(1,idMtke)
          Vinfo( 2)=Vname(2,idMtke)
          Vinfo( 3)=Vname(3,idMtke)
          Vinfo(14)=Vname(4,idMtke)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idMtke,ng),r8)
          status=def_var(ncHISid(ng),hisVid(idMtke,ng),NF_FOUT,         &
     &                   nvd4,w3dgrd,Aval,Vinfo,ncname)
#  ifdef FORWARD_WRITE
          Vinfo( 1)=Vname(1,idVmKK)
          Vinfo( 2)=Vname(2,idVmKK)
          Vinfo( 3)=Vname(3,idVmKK)
          Vinfo(14)=Vname(4,idVmKK)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVmKK,ng),r8)
          status=def_var(ncHISid(ng),hisVid(idVmKK,ng),NF_FOUT,         &
     &                   nvd4,w3dgrd,Aval,Vinfo,ncname)
#  endif
        END IF
!
!  Define turbulent kinetic energy time length scale.
!
        IF (Hout(idMtls,ng)) THEN
          Vinfo( 1)=Vname(1,idMtls)
          Vinfo( 2)=Vname(2,idMtls)
          Vinfo( 3)=Vname(3,idMtls)
          Vinfo(14)=Vname(4,idMtls)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idMtls,ng),r8)
          status=def_var(ncHISid(ng),hisVid(idMtls,ng),NF_FOUT,         &
     &                   nvd4,w3dgrd,Aval,Vinfo,ncname)
#  ifdef FORWARD_WRITE
          Vinfo( 1)=Vname(1,idVmLS)
          Vinfo( 2)=Vname(2,idVmLS)
          Vinfo( 3)=Vname(3,idVmLS)
          Vinfo(14)=Vname(4,idVmLS)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVmLS,ng),r8)
          status=def_var(ncHISid(ng),hisVid(idVmLS,ng),NF_FOUT,         &
     &                   nvd4,w3dgrd,Aval,Vinfo,ncname)
#  endif
#  if defined FORWARD_WRITE && defined GLS_MIXING
          Vinfo( 1)=Vname(1,idVmKP)
          Vinfo( 2)=Vname(2,idVmKP)
          Vinfo( 3)=Vname(3,idVmKP)
          Vinfo(14)=Vname(4,idVmKP)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVmKP,ng),r8)
          status=def_var(ncHISid(ng),hisVid(idVmKP,ng),NF_FOUT,         &
     &                   nvd4,w3dgrd,Aval,Vinfo,ncname)
#  endif
        END IF
# endif
# ifndef ADJUST_STFLUX
!
!  Define surface active tracer fluxes.
!
        DO itrc=1,NAT
          IF (Hout(idTsur(itrc),ng)) THEN
            Vinfo( 1)=Vname(1,idTsur(itrc))
            Vinfo( 2)=Vname(2,idTsur(itrc))
            Vinfo( 3)=Vname(3,idTsur(itrc))
            IF (itrc.eq.itemp) THEN
              Vinfo(11)='upward flux, cooling'
              Vinfo(12)='downward flux, heating'
            ELSE IF (itrc.eq.isalt) THEN
              Vinfo(11)='upward flux, freshening (net precipitation)'
              Vinfo(12)='downward flux, salting (net evaporation)'
            END IF
            Vinfo(14)=Vname(4,idTsur(itrc))
            Vinfo(16)=Vname(1,idtime)
            Vinfo(22)='coordinates'
            Aval(5)=REAL(Iinfo(1,idTsur(itrc),ng),r8)
            status=def_var(ncHISid(ng),hisVid(idTsur(itrc),ng),         &
     &                     NF_FOUT,nvd3,t2dgrd,Aval,Vinfo,ncname)
          END IF
        END DO
# endif
# ifdef BULK_FLUXES
!
!  Define latent heat flux.
!
        IF (Hout(idLhea,ng)) THEN
          Vinfo( 1)=Vname(1,idLhea)
          Vinfo( 2)=Vname(2,idLhea)
          Vinfo( 3)=Vname(3,idLhea)
          Vinfo(11)='upward flux, cooling'
          Vinfo(12)='downward flux, heating'
          Vinfo(14)=Vname(4,idLhea)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idLhea,ng),r8)
          status=def_var(ncHISid(ng),hisVid(idLhea,ng),NF_FOUT,         &
     &                   nvd3,t2dgrd,Aval,Vinfo,ncname)
        END IF
!
!  Define sensible heat flux.
!
        IF (Hout(idShea,ng)) THEN
          Vinfo( 1)=Vname(1,idShea)
          Vinfo( 2)=Vname(2,idShea)
          Vinfo( 3)=Vname(3,idShea)
          Vinfo(11)='upward flux, cooling'
          Vinfo(12)='downward flux, heating'
          Vinfo(14)=Vname(4,idShea)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idShea,ng),r8)
          status=def_var(ncHISid(ng),hisVid(idShea,ng),NF_FOUT,         &
     &                   nvd3,t2dgrd,Aval,Vinfo,ncname)
        END IF
!
!  Define longwave radiation flux.
!
        IF (Hout(idLrad,ng)) THEN
          Vinfo( 1)=Vname(1,idLrad)
          Vinfo( 2)=Vname(2,idLrad)
          Vinfo( 3)=Vname(3,idLrad)
          Vinfo(11)='upward flux, cooling'
          Vinfo(12)='downward flux, heating'
          Vinfo(14)=Vname(4,idLrad)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idLrad,ng),r8)
          status=def_var(ncHISid(ng),hisVid(idLrad,ng),NF_FOUT,         &
     &                   nvd3,t2dgrd,Aval,Vinfo,ncname)
        END IF
!
!  Define atmospheric air temperature.
!
        IF (Hout(idTair,ng)) THEN
          Vinfo( 1)=Vname(1,idTair)
          Vinfo( 2)=Vname(2,idTair)
          Vinfo( 3)=Vname(3,idTair)
          Vinfo(14)=Vname(4,idTair)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idTair,ng),r8)
          status=def_var(ncHISid(ng),hisVid(idTair,ng),NF_FOUT,         &
     &                   nvd3,t2dgrd,Aval,Vinfo,ncname)
        END IF
#  ifdef EMINUSP
!
!  Define evaporation rate.
!
        IF (Hout(idevap,ng)) THEN
          Vinfo( 1)=Vname(1,idevap)
          Vinfo( 2)=Vname(2,idevap)
          Vinfo( 3)=Vname(3,idevap)
          Vinfo(11)='downward flux, freshening (condensation)'
          Vinfo(12)='upward flux, salting (evaporation)'
          Vinfo(14)=Vname(4,idevap)
          Vinfo(16)=Vname(1,idtime)
#     if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#     endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idevap,ng),r8)
          status=def_var(ncHISid(ng),hisVid(idevap,ng),NF_FOUT,         &
     &                   nvd3,t2dgrd,Aval,Vinfo,ncname)
        END IF
!
!  Define precipitation rate.
!
        IF (Hout(idrain,ng)) THEN
          Vinfo( 1)=Vname(1,idrain)
          Vinfo( 2)=Vname(2,idrain)
          Vinfo( 3)=Vname(3,idrain)
          Vinfo(11)='upward flux, salting (NOT POSSIBLE)'
          Vinfo(12)='downward flux, freshening (precipitation)'
          Vinfo(14)=Vname(4,idrain)
          Vinfo(16)=Vname(1,idtime)
#     if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#     endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idrain,ng),r8)
          status=def_var(ncHISid(ng),hisVid(idrain,ng),NF_FOUT,         &
     &                   nvd3,t2dgrd,Aval,Vinfo,ncname)
        END IF
#  endif
# endif
# ifdef SHORTWAVE
!
!  Define shortwave radiation flux.
!
        IF (Hout(idSrad,ng)) THEN
          Vinfo( 1)=Vname(1,idSrad)
          Vinfo( 2)=Vname(2,idSrad)
          Vinfo( 3)=Vname(3,idSrad)
          Vinfo(11)='upward flux, cooling'
          Vinfo(12)='downward flux, heating'
          Vinfo(14)=Vname(4,idSrad)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idSrad,ng),r8)
          status=def_var(ncHISid(ng),hisVid(idSrad,ng),NF_FOUT,         &
     &                   nvd3,t2dgrd,Aval,Vinfo,ncname)
        END IF
# endif
#endif
#ifndef ADJUST_WSTRESS
!
!  Define surface U-momentum stress.
!
        IF (Hout(idUsms,ng)) THEN
          Vinfo( 1)=Vname(1,idUsms)
          Vinfo( 2)=Vname(2,idUsms)
          Vinfo( 3)=Vname(3,idUsms)
          Vinfo(14)=Vname(4,idUsms)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
# endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUsms,ng),r8)
          status=def_var(ncHISid(ng),hisVid(idUsms,ng),NF_FOUT,         &
     &                   nvd3,u2dgrd,Aval,Vinfo,ncname)
        END IF
!
!  Define surface V-momentum stress.
!
        IF (Hout(idVsms,ng)) THEN
          Vinfo( 1)=Vname(1,idVsms)
          Vinfo( 2)=Vname(2,idVsms)
          Vinfo( 3)=Vname(3,idVsms)
          Vinfo(14)=Vname(4,idVsms)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
# endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVsms,ng),r8)
          status=def_var(ncHISid(ng),hisVid(idVsms,ng),NF_FOUT,         &
     &                   nvd3,v2dgrd,Aval,Vinfo,ncname)
        END IF
#endif
!
!  Define bottom U-momentum stress.
!
        IF (Hout(idUbms,ng)) THEN
          Vinfo( 1)=Vname(1,idUbms)
          Vinfo( 2)=Vname(2,idUbms)
          Vinfo( 3)=Vname(3,idUbms)
          Vinfo(14)=Vname(4,idUbms)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
# endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUbms,ng),r8)
          status=def_var(ncHISid(ng),hisVid(idUbms,ng),NF_FOUT,         &
     &                   nvd3,u2dgrd,Aval,Vinfo,ncname)
        END IF
!
!  Define bottom V-momentum stress.
!
        IF (Hout(idVbms,ng)) THEN
          Vinfo( 1)=Vname(1,idVbms)
          Vinfo( 2)=Vname(2,idVbms)
          Vinfo( 3)=Vname(3,idVbms)
          Vinfo(14)=Vname(4,idVbms)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
# endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVbms,ng),r8)
          status=def_var(ncHISid(ng),hisVid(idVbms,ng),NF_FOUT,         &
     &                   nvd3,v2dgrd,Aval,Vinfo,ncname)
        END IF
#ifdef SOLVE3D
# ifdef BBL_MODEL
!
!  Define bottom U-current stress.
!
        IF (Hout(idUbrs,ng)) THEN
          Vinfo( 1)=Vname(1,idUbrs)
          Vinfo( 2)=Vname(2,idUbrs)
          Vinfo( 3)=Vname(3,idUbrs)
          Vinfo(14)=Vname(4,idUbrs)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_r'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUbrs,ng),r8)
          status=def_var(ncHISid(ng),hisVid(idUbrs,ng),NF_FOUT,         &
     &                   nvd3,t2dgrd,Aval,Vinfo,ncname)
        END IF
!
!  Define bottom V-current stress.
!
        IF (Hout(idVbrs,ng)) THEN
          Vinfo( 1)=Vname(1,idVbrs)
          Vinfo( 2)=Vname(2,idVbrs)
          Vinfo( 3)=Vname(3,idVbrs)
          Vinfo(14)=Vname(4,idVbrs)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_r'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVbrs,ng),r8)
          status=def_var(ncHISid(ng),hisVid(idVbrs,ng),NF_FOUT,         &
     &                   nvd3,t2dgrd,Aval,Vinfo,ncname)
        END IF
!
!  Define wind-induced, bottom U-wave stress.
!
        IF (Hout(idUbws,ng)) THEN
          Vinfo( 1)=Vname(1,idUbws)
          Vinfo( 2)=Vname(2,idUbws)
          Vinfo( 3)=Vname(3,idUbws)
          Vinfo(14)=Vname(4,idUbws)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_r'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUbws,ng),r8)
          status=def_var(ncHISid(ng),hisVid(idUbws,ng),NF_FOUT,         &
     &                   nvd3,t2dgrd,Aval,Vinfo,ncname)
        END IF
!
!  Define bottom wind-induced, bottom V-wave stress.
!
        IF (Hout(idVbws,ng)) THEN
          Vinfo( 1)=Vname(1,idVbws)
          Vinfo( 2)=Vname(2,idVbws)
          Vinfo( 3)=Vname(3,idVbws)
          Vinfo(14)=Vname(4,idVbws)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_r'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVbws,ng),r8)
          status=def_var(ncHISid(ng),hisVid(idVbws,ng),NF_FOUT,         &
     &                   nvd3,t2dgrd,Aval,Vinfo,ncname)
        END IF
!
!  Define maximum wind and current, bottom U-wave stress.
!
        IF (Hout(idUbcs,ng)) THEN
          Vinfo( 1)=Vname(1,idUbcs)
          Vinfo( 2)=Vname(2,idUbcs)
          Vinfo( 3)=Vname(3,idUbcs)
          Vinfo(14)=Vname(4,idUbcs)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_r'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUbcs,ng),r8)
          status=def_var(ncHISid(ng),hisVid(idUbcs,ng),NF_FOUT,         &
     &                   nvd3,t2dgrd,Aval,Vinfo,ncname)
        END IF
!
!  Define maximum wind and current, bottom V-wave stress.
!
        IF (Hout(idVbcs,ng)) THEN
          Vinfo( 1)=Vname(1,idVbcs)
          Vinfo( 2)=Vname(2,idVbcs)
          Vinfo( 3)=Vname(3,idVbcs)
          Vinfo(14)=Vname(4,idVbcs)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_r'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVbcs,ng),r8)
          status=def_var(ncHISid(ng),hisVid(idVbcs,ng),NF_FOUT,         &
     &                   nvd3,t2dgrd,Aval,Vinfo,ncname)
        END IF
!
!  Define wind-induced, bed wave orbital U-velocity.
!
        IF (Hout(idUbot,ng)) THEN
          Vinfo( 1)=Vname(1,idUbot)
          Vinfo( 2)=Vname(2,idUbot)
          Vinfo( 3)=Vname(3,idUbot)
          Vinfo(14)=Vname(4,idUbot)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_r'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUbot,ng),r8)
          status=def_var(ncHISid(ng),hisVid(idUbot,ng),NF_FOUT,         &
     &                   nvd3,t2dgrd,Aval,Vinfo,ncname)
        END IF
!
!  Define wind-induced, bed wave orbital V-velocity.
!
        IF (Hout(idVbot,ng)) THEN
          Vinfo( 1)=Vname(1,idVbot)
          Vinfo( 2)=Vname(2,idVbot)
          Vinfo( 3)=Vname(3,idVbot)
          Vinfo(14)=Vname(4,idVbot)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_r'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVbot,ng),r8)
          status=def_var(ncHISid(ng),hisVid(idVbot,ng),NF_FOUT,         &
     &                   nvd3,t2dgrd,Aval,Vinfo,ncname)
        END IF
!
!  Define bottom U-momentum above bed.
!
        IF (Hout(idUbur,ng)) THEN
          Vinfo( 1)=Vname(1,idUbur)
          Vinfo( 2)=Vname(2,idUbur)
          Vinfo( 3)=Vname(3,idUbur)
          Vinfo(14)=Vname(4,idUbur)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_r'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUbur,ng),r8)
          status=def_var(ncHISid(ng),hisVid(idUbur,ng),NF_FOUT,         &
     &                   nvd3,t2dgrd,Aval,Vinfo,ncname)
        END IF
!
!  Define bottom V-momentum above bed.
!
        IF (Hout(idVbvr,ng)) THEN
          Vinfo( 1)=Vname(1,idVbvr)
          Vinfo( 2)=Vname(2,idVbvr)
          Vinfo( 3)=Vname(3,idVbvr)
          Vinfo(14)=Vname(4,idVbvr)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_r'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVbvr,ng),r8)
          status=def_var(ncHISid(ng),hisVid(idVbvr,ng),NF_FOUT,         &
     &                   nvd3,t2dgrd,Aval,Vinfo,ncname)
        END IF
# endif
# ifdef SEDIMENT
#  ifdef BEDLOAD
!
!  Define Bedload U-direction.
!
        DO i=1,NST
          IF (Hout(idUbld(i),ng)) THEN
            Vinfo( 1)=Vname(1,idUbld(i))
            Vinfo( 2)=Vname(2,idUbld(i))
            Vinfo( 3)=Vname(3,idUbld(i))
            Vinfo(14)=Vname(4,idUbld(i))
            Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
            Vinfo(20)='mask_u'
#   endif
            Vinfo(22)='coordinates'
            Aval(5)=REAL(Iinfo(1,idUbld(i),ng),r8)
            status=def_var(ncHISid(ng),hisVid(idUbld(i),ng),NF_FOUT,    &
     &                     nvd3,u2dgrd,Aval,Vinfo,ncname)
          END IF
!
!  Define Bedload V-direction.
!
          IF (Hout(idVbld(i),ng)) THEN
            Vinfo( 1)=Vname(1,idVbld(i))
            Vinfo( 2)=Vname(2,idVbld(i))
            Vinfo( 3)=Vname(3,idVbld(i))
            Vinfo(14)=Vname(4,idVbld(i))
            Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
            Vinfo(20)='mask_v'
#   endif
            Vinfo(22)='coordinates'
            Aval(5)=REAL(Iinfo(1,idVbld(i),ng),r8)
            status=def_var(ncHISid(ng),hisVid(idVbld(i),ng),NF_FOUT,    &
     &                     nvd3,v2dgrd,Aval,Vinfo,ncname)
          END IF
        END DO
#  endif
!
!  Define sediment fraction of each size class in each bed layer.
!
        DO i=1,NST
          IF (Hout(idfrac(i),ng)) THEN
            Vinfo( 1)=Vname(1,idfrac(i))
            Vinfo( 2)=Vname(2,idfrac(i))
            Vinfo( 3)=Vname(3,idfrac(i))
            Vinfo(14)=Vname(4,idfrac(i))
            Vinfo(16)=Vname(1,idtime)
            WRITE (Vinfo(19),40) 1000.0_r8*Sd50(i,ng)
#  if defined WRITE_WATER && defined MASKING
            Vinfo(20)='mask_rho'
#  endif
            Vinfo(22)='coordinates'
            Aval(5)=REAL(Iinfo(1,idfrac(i),ng))
            status=def_var(ncHISid(ng),hisVid(idfrac(i),ng),NF_FOUT,    &
     &                     nvd4,b3dgrd,Aval,Vinfo,ncname)
          END IF
        END DO
!
!  Define sediment mass of each size class in each bed layer.
!
        DO i=1,NST
          IF (Hout(idBmas(i),ng)) THEN
            Vinfo( 1)=Vname(1,idBmas(i))
            Vinfo( 2)=Vname(2,idBmas(i))
            Vinfo( 3)=Vname(3,idBmas(i))
            Vinfo(14)=Vname(4,idBmas(i))
            Vinfo(16)=Vname(1,idtime)
            WRITE (Vinfo(19),40) 1000.0_r8*Sd50(i,ng)
#  if defined WRITE_WATER && defined MASKING
            Vinfo(20)='mask_rho'
#  endif
            Vinfo(22)='coordinates'
            Aval(5)=REAL(Iinfo(1,idBmas(i),ng),r8)
            status=def_var(ncHISid(ng),hisVid(idBmas(i),ng),NF_FOUT,    &
     &                     nvd4,b3dgrd,Aval,Vinfo,ncname)
          END IF
        END DO
!
!  Define sediment properties in each bed layer.
!
        DO i=1,MBEDP
          IF (Hout(idSbed(i),ng)) THEN
            Vinfo( 1)=Vname(1,idSbed(i))
            Vinfo( 2)=Vname(2,idSbed(i))
            Vinfo( 3)=Vname(3,idSbed(i))
            Vinfo(14)=Vname(4,idSbed(i))
            Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
            Vinfo(20)='mask_rho'
#  endif
            Vinfo(22)='coordinates'
            Aval(5)=REAL(Iinfo(1,idSbed(i),ng),r8)
            status=def_var(ncHISid(ng),hisVid(idSbed(i),ng),NF_FOUT,    &
     &                     nvd4,b3dgrd,Aval,Vinfo,ncname)
          END IF
        END DO
# endif
# if defined SEDIMENT || defined BBL_MODEL
!
!  Define exposed sediment layer properties.
!
        DO i=1,MBOTP
          IF (Hout(idBott(i),ng)) THEN
            Vinfo( 1)=Vname(1,idBott(i))
            Vinfo( 2)=Vname(2,idBott(i))
            Vinfo( 3)=Vname(3,idBott(i))
            Vinfo(14)=Vname(4,idBott(i))
            Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
            Vinfo(22)='coordinates'
            Aval(5)=REAL(r2dvar,r8)
            status=def_var(ncHISid(ng),hisVid(idBott(i),ng),NF_FOUT,    &
     &                     nvd3,t2dgrd,Aval,Vinfo,ncname)
          END IF
        END DO
# endif
#endif
#ifdef NEARSHORE_MELLOR
!
!  Define 2D radiation stress, Sxx-component.
!
        IF (Hout(idW2xx,ng)) THEN
          Vinfo( 1)=Vname(1,idW2xx)
          Vinfo( 2)=Vname(2,idW2xx)
          Vinfo( 3)=Vname(3,idW2xx)
          Vinfo(14)=Vname(4,idW2xx)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
# endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(r2dvar,r8)
          status=def_var(ncHISid(ng),hisVid(idW2xx,ng),NF_FOUT,         &
     &                   nvd3,t2dgrd,Aval,Vinfo,ncname)
        END IF
!
!  Define 2D radiation stress, Sxy-component.
!
        IF (Hout(idW2xy,ng)) THEN
          Vinfo( 1)=Vname(1,idW2xy)
          Vinfo( 2)=Vname(2,idW2xy)
          Vinfo( 3)=Vname(3,idW2xy)
          Vinfo(14)=Vname(4,idW2xy)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
# endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(r2dvar,r8)
          status=def_var(ncHISid(ng),hisVid(idW2xy,ng),NF_FOUT,         &
     &                   nvd3,t2dgrd,Aval,Vinfo,ncname)
        END IF
!
!  Define 2D radiation stress, Syy-component.
!
        IF (Hout(idW2yy,ng)) THEN
          Vinfo( 1)=Vname(1,idW2yy)
          Vinfo( 2)=Vname(2,idW2yy)
          Vinfo( 3)=Vname(3,idW2yy)
          Vinfo(14)=Vname(4,idW2yy)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
# endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(r2dvar,r8)
          status=def_var(ncHISid(ng),hisVid(idW2yy,ng),NF_FOUT,         &
     &                   nvd3,t2dgrd,Aval,Vinfo,ncname)
        END IF
!
!  Define 2D total u-radiation stress.
!
        IF (Hout(idU2rs,ng)) THEN
          Vinfo( 1)=Vname(1,idU2rs)
          Vinfo( 2)=Vname(2,idU2rs)
          Vinfo( 3)=Vname(3,idU2rs)
          Vinfo(14)=Vname(4,idU2rs)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
# endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(u2dvar,r8)
          status=def_var(ncHISid(ng),hisVid(idU2rs,ng),NF_FOUT,         &
     &                   nvd3,u2dgrd,Aval,Vinfo,ncname)
        END IF
!
!  Define 2D total v-radiation stress.
!
        IF (Hout(idV2rs,ng)) THEN
          Vinfo( 1)=Vname(1,idV2rs)
          Vinfo( 2)=Vname(2,idV2rs)
          Vinfo( 3)=Vname(3,idV2rs)
          Vinfo(14)=Vname(4,idV2rs)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
# endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(v2dvar,r8)
          status=def_var(ncHISid(ng),hisVid(idV2rs,ng),NF_FOUT,         &
     &                   nvd3,v2dgrd,Aval,Vinfo,ncname)
        END IF
!
!  Define 2D u-Stokes drift velocity.
!
        IF (Hout(idU2Sd,ng)) THEN
          Vinfo( 1)=Vname(1,idU2Sd)
          Vinfo( 2)=Vname(2,idU2Sd)
          Vinfo( 3)=Vname(3,idU2Sd)
          Vinfo(14)=Vname(4,idU2Sd)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
# endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idU2Sd,ng),r8)
          status=def_var(ncHISid(ng),hisVid(idU2Sd,ng),NF_FOUT,         &
     &                   nvd3,u2dgrd,Aval,Vinfo,ncname)
        END IF
!
!  Define 2D v-Stokes drift velocity.
!
        IF (Hout(idV2Sd,ng)) THEN
          Vinfo( 1)=Vname(1,idV2Sd)
          Vinfo( 2)=Vname(2,idV2Sd)
          Vinfo( 3)=Vname(3,idV2Sd)
          Vinfo(14)=Vname(4,idV2Sd)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
# endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idV2Sd,ng),r8)
          status=def_var(ncHISid(ng),hisVid(idV2Sd,ng),NF_FOUT,         &
     &                   nvd3,v2dgrd,Aval,Vinfo,ncname)
        END IF

# ifdef SOLVE3D
!
!  Define 3D radiation stress, Sxx-component.
!
        IF (Hout(idW3xx,ng)) THEN
          Vinfo( 1)=Vname(1,idW3xx)
          Vinfo( 2)=Vname(2,idW3xx)
          Vinfo( 3)=Vname(3,idW3xx)
          Vinfo(14)=Vname(4,idW3xx)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
# endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(r3dvar,r8)
          status=def_var(ncHISid(ng),hisVid(idW3xx,ng),NF_FOUT,         &
     &                   nvd4,t3dgrd,Aval,Vinfo,ncname)
        END IF
!
!  Define 3D radiation stress, Sxy-component.
!
        IF (Hout(idW3xy,ng)) THEN
          Vinfo( 1)=Vname(1,idW3xy)
          Vinfo( 2)=Vname(2,idW3xy)
          Vinfo( 3)=Vname(3,idW3xy)
          Vinfo(14)=Vname(4,idW3xy)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
# endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(r3dvar,r8)
          status=def_var(ncHISid(ng),hisVid(idW3xy,ng),NF_FOUT,         &
     &                   nvd4,t3dgrd,Aval,Vinfo,ncname)
        END IF
!
!  Define 3D radiation stress, Syy-component.
!
        IF (Hout(idW3yy,ng)) THEN
          Vinfo( 1)=Vname(1,idW3yy)
          Vinfo( 2)=Vname(2,idW3yy)
          Vinfo( 3)=Vname(3,idW3yy)
          Vinfo(14)=Vname(4,idW3yy)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
# endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(r3dvar,r8)
          status=def_var(ncHISid(ng),hisVid(idW3yy,ng),NF_FOUT,         &
     &                   nvd4,t3dgrd,Aval,Vinfo,ncname)
        END IF
!
!  Define 3D radiation stress, Spx-component.
!
        IF (Hout(idW3zx,ng)) THEN
          Vinfo( 1)=Vname(1,idW3zx)
          Vinfo( 2)=Vname(2,idW3zx)
          Vinfo( 3)=Vname(3,idW3zx)
          Vinfo(14)=Vname(4,idW3zx)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
# endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(r3dvar,r8)
          status=def_var(ncHISid(ng),hisVid(idW3zx,ng),NF_FOUT,         &
     &                   nvd4,t3dgrd,Aval,Vinfo,ncname)
        END IF
!
!  Define 3D radiation stress, Spy-component.
!
        IF (Hout(idW3zy,ng)) THEN
          Vinfo( 1)=Vname(1,idW3zy)
          Vinfo( 2)=Vname(2,idW3zy)
          Vinfo( 3)=Vname(3,idW3zy)
          Vinfo(14)=Vname(4,idW3zy)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
# endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(r3dvar,r8)
          status=def_var(ncHISid(ng),hisVid(idW3zy,ng),NF_FOUT,         &
     &                   nvd4,t3dgrd,Aval,Vinfo,ncname)
        END IF
!
!  Define 3D total u-radiation stress.
!
        IF (Hout(idU3rs,ng)) THEN
          Vinfo( 1)=Vname(1,idU3rs)
          Vinfo( 2)=Vname(2,idU3rs)
          Vinfo( 3)=Vname(3,idU3rs)
          Vinfo(14)=Vname(4,idU3rs)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
# endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(u3dvar,r8)
          status=def_var(ncHISid(ng),hisVid(idU3rs,ng),NF_FOUT,         &
     &                   nvd4,u3dgrd,Aval,Vinfo,ncname)
        END IF
!
!  Define 3D total v-radiation stress.
!
        IF (Hout(idV3rs,ng)) THEN
          Vinfo( 1)=Vname(1,idV3rs)
          Vinfo( 2)=Vname(2,idV3rs)
          Vinfo( 3)=Vname(3,idV3rs)
          Vinfo(14)=Vname(4,idV3rs)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
# endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(v3dvar,r8)
          status=def_var(ncHISid(ng),hisVid(idV3rs,ng),NF_FOUT,         &
     &                   nvd4,v3dgrd,Aval,Vinfo,ncname)
        END IF
!
!  Define 3D u-Stokes velocity.
!
        IF (Hout(idU3Sd,ng)) THEN
          Vinfo( 1)=Vname(1,idU3Sd)
          Vinfo( 2)=Vname(2,idU3Sd)
          Vinfo( 3)=Vname(3,idU3Sd)
          Vinfo(14)=Vname(4,idU3Sd)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idU3Sd,ng),r8)
          status=def_var(ncHISid(ng),hisVid(idU3Sd,ng),NF_FOUT,         &
     &                   nvd4,u3dgrd,Aval,Vinfo,ncname)
        END IF
!
!  Define 3D v-Stokes velocity.
!
        IF (Hout(idV3Sd,ng)) THEN
          Vinfo( 1)=Vname(1,idV3Sd)
          Vinfo( 2)=Vname(2,idV3Sd)
          Vinfo( 3)=Vname(3,idV3Sd)
          Vinfo(14)=Vname(4,idV3Sd)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idV3Sd,ng),r8)
          status=def_var(ncHISid(ng),hisVid(idV3Sd,ng),NF_FOUT,         &
     &                   nvd4,v3dgrd,Aval,Vinfo,ncname)
        END IF
# endif
#endif
#ifdef WET_DRY
!
!  Define wet_dry masking.
!
!        IF (Hout(idWmsk,ng)) THEN
          Vinfo( 1)=Vname(1,idWmsk)
          Vinfo( 2)=Vname(2,idWmsk)
          Vinfo( 3)=Vname(3,idWmsk)
          Vinfo(14)=Vname(4,idWmsk)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
# endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idWmsk,ng),r8)
          status=def_var(ncHISid(ng),hisVid(idWmsk,ng),NF_FOUT,         &
     &                   nvd3,t2dgrd,Aval,Vinfo,ncname)
!        END IF
#endif
#ifdef WAVES_HEIGHT
!
!  Define Hwave.
!
        IF (Hout(idWamp,ng)) THEN
          Vinfo( 1)=Vname(1,idWamp)
          Vinfo( 2)=Vname(2,idWamp)
          Vinfo( 3)=Vname(3,idWamp)
          Vinfo(14)=Vname(4,idWamp)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
# endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idWamp,ng),r8)
          status=def_var(ncHISid(ng),hisVid(idWamp,ng),NF_FOUT,         &
     &                   nvd3,t2dgrd,Aval,Vinfo,ncname)
        END IF
#endif
#ifdef WAVES_LENGTH
!
!  Define Lwave.
!
        IF (Hout(idWlen,ng)) THEN
          Vinfo( 1)=Vname(1,idWlen)
          Vinfo( 2)=Vname(2,idWlen)
          Vinfo( 3)=Vname(3,idWlen)
          Vinfo(14)=Vname(4,idWlen)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
# endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idWlen,ng),r8)
          status=def_var(ncHISid(ng),hisVid(idWlen,ng),NF_FOUT,         &
     &                   nvd3,t2dgrd,Aval,Vinfo,ncname)
        END IF
#endif
#ifdef WAVES_DIR
!
!  Define Dwave.
!
        IF (Hout(idWdir,ng)) THEN
          Vinfo( 1)=Vname(1,idWdir)
          Vinfo( 2)=Vname(2,idWdir)
          Vinfo( 3)=Vname(3,idWdir)
          Vinfo(14)=Vname(4,idWdir)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
# endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idWdir,ng),r8)
          status=def_var(ncHISid(ng),hisVid(idWdir,ng),NF_FOUT,         &
     &                   nvd3,t2dgrd,Aval,Vinfo,ncname)
        END IF
#endif
!
!-----------------------------------------------------------------------
!  Leave definition mode.
!-----------------------------------------------------------------------
!
        status=nf90_enddef(ncHISid(ng))
      END IF
!
!-----------------------------------------------------------------------
!  Write out time-recordless, information variables.
!-----------------------------------------------------------------------
!
      IF (ldef) THEN
        status=wrt_info(ng, iNLM, ncHISid(ng), OutThread, ncname)
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!=======================================================================
!  Open an existing history file, check its contents, and prepare for
!  appending data.
!=======================================================================
!
      IF (.not.ldef.and.OutThread) THEN
!
!  Inquire about the contents of history NetCDF file:  Inquire about
!  the dimensions and variables.  Check for consistency.
!
        CALL opencdf (ng, 1, ncname, fname, N(ng), 0, nrec, nvd, Vsize)
        IF (exit_flag.ne.NoError) RETURN
!
!  Open history file for read/write.
!
        status=nf90_open(TRIM(ncname),nf90_write,ncHISid(ng))
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
!  history variables. Get variable IDs.
!
        DO i=1,nvars
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idtime))) THEN
            got_var(idtime)=.TRUE.
            status=nf90_inq_varid(ncHISid(ng),TRIM(Vname(1,idtime)),    &
     &                            hisVid(idtime,ng))
          END IF
#if defined SEDIMENT && defined SED_MORPH
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idbath))) THEN
            got_var(idbath)=.TRUE.
            status=nf90_inq_varid(ncHISid(ng),TRIM(Vname(1,idbath)),    &
     &                            hisVid(idbath,ng))
          END IF
#endif
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idFsur))) THEN
            got_var(idFsur)=.TRUE.
            status=nf90_inq_varid(ncHISid(ng),TRIM(Vname(1,idFsur)),    &
     &                            hisVid(idFsur,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idUbar))) THEN
            got_var(idUbar)=.TRUE.
            status=nf90_inq_varid(ncHISid(ng),TRIM(Vname(1,idUbar)),    &
     &                            hisVid(idUbar,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idVbar))) THEN
            got_var(idVbar)=.TRUE.
            status=nf90_inq_varid(ncHISid(ng),TRIM(Vname(1,idVbar)),    &
     &                            hisVid(idVbar,ng))
          END IF
#ifdef FORWARD_WRITE
# ifdef FORWARD_RHS
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idRzet))) THEN
            got_var(idRzet)=.TRUE.
            status=nf90_inq_varid(ncHISid(ng),TRIM(Vname(1,idRzet)),    &
     &                            hisVid(idRzet,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idRu2d))) THEN
            got_var(idRu2d)=.TRUE.
            status=nf90_inq_varid(ncHISid(ng),TRIM(Vname(1,idRu2d)),    &
     &                            hisVid(idRu2d,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idRv2d))) THEN
            got_var(idRv2d)=.TRUE.
            status=nf90_inq_varid(ncHISid(ng),TRIM(Vname(1,idRv2d)),    &
     &                            hisVid(idRv2d,ng))
          END IF
# endif
# ifdef SOLVE3D
#  ifdef FORWARD_RHS
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idRuct))) THEN
            got_var(idRuct)=.TRUE.
            status=nf90_inq_varid(ncHISid(ng),TRIM(Vname(1,idRuct)),    &
     &                            hisVid(idRuct,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idRvct))) THEN
            got_var(idRvct)=.TRUE.
            status=nf90_inq_varid(ncHISid(ng),TRIM(Vname(1,idRvct)),    &
     &                            hisVid(idRvct,ng))
          END IF
#  endif
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idUfx1))) THEN
            got_var(idUfx1)=.TRUE.
            status=nf90_inq_varid(ncHISid(ng),TRIM(Vname(1,idUfx1)),    &
     &                            hisVid(idUfx1,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idUfx1))) THEN
            got_var(idUfx1)=.TRUE.
            status=nf90_inq_varid(ncHISid(ng),TRIM(Vname(1,idUfx1)),    &
     &                            hisVid(idUfx1,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idUfx2))) THEN
            got_var(idUfx2)=.TRUE.
            status=nf90_inq_varid(ncHISid(ng),TRIM(Vname(1,idUfx2)),    &
     &                            hisVid(idUfx2,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idVfx1))) THEN
            got_var(idVfx1)=.TRUE.
            status=nf90_inq_varid(ncHISid(ng),TRIM(Vname(1,idVfx1)),    &
     &                            hisVid(idVfx1,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idVfx2))) THEN
            got_var(idVfx2)=.TRUE.
            status=nf90_inq_varid(ncHISid(ng),TRIM(Vname(1,idVfx2)),    &
     &                            hisVid(idVfx2,ng))
          END IF
#  ifdef FORWARD_RHS
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idRu3d))) THEN
            got_var(idRu3d)=.TRUE.
            status=nf90_inq_varid(ncHISid(ng),TRIM(Vname(1,idRu3d)),    &
     &                            hisVid(idRu3d,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idRv3d))) THEN
            got_var(idRv3d)=.TRUE.
            status=nf90_inq_varid(ncHISid(ng),TRIM(Vname(1,idRv3d)),    &
     &                            hisVid(idRv3d,ng))
          END IF
#  endif
# endif
#endif
#ifdef SOLVE3D
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idUvel))) THEN
            got_var(idUvel)=.TRUE.
            status=nf90_inq_varid(ncHISid(ng),TRIM(Vname(1,idUvel)),    &
     &                            hisVid(idUvel,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idVvel))) THEN
            got_var(idVvel)=.TRUE.
            status=nf90_inq_varid(ncHISid(ng),TRIM(Vname(1,idVvel)),    &
     &                            hisVid(idVvel,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idWvel))) THEN
            got_var(idWvel)=.TRUE.
            status=nf90_inq_varid(ncHISid(ng),TRIM(Vname(1,idWvel)),    &
     &                            hisVid(idWvel,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idOvel))) THEN
            got_var(idOvel)=.TRUE.
            status=nf90_inq_varid(ncHISid(ng),TRIM(Vname(1,idOvel)),    &
     &                            hisVid(idOvel,ng))
          END IF
          DO itrc=1,NT(ng)
            IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idTvar(itrc)))) THEN
              got_var(idTvar(itrc))=.TRUE.
              status=nf90_inq_varid(ncHISid(ng),                        &
     &                              TRIM(Vname(1,idTvar(itrc))),        &
     &                              hisTid(itrc,ng))
            END IF
          END DO
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idDano))) THEN
            got_var(idDano)=.TRUE.
            status=nf90_inq_varid(ncHISid(ng),TRIM(Vname(1,idDano)),    &
     &                            hisVid(idDano,ng))
          END IF
# ifdef LMD_SKPP
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idHsbl))) THEN
            got_var(idHsbl)=.TRUE.
            status=nf90_inq_varid(ncHISid(ng),TRIM(Vname(1,idHsbl)),    &
     &                            hisVid(idHsbl,ng))
          END IF
# endif
# ifdef LMD_BKPP
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idHbbl))) THEN
            got_var(idHbbl)=.TRUE.
            status=nf90_inq_varid(ncHISid(ng),TRIM(Vname(1,idHbbl)),    &
     &                            hisVid(idHbbl,ng))
          END IF
# endif
# if defined FORWARD_WRITE && defined LMD_NONLOCAL
          DO itrc=1,NAT
            IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idGhat(itrc)))) THEN
              got_var(idGhat(itrc))=.TRUE.
              status=nf90_inq_varid(ncHISid(ng),                        &
     &                              TRIM(Vname(1,idGhat(itrc))),        &
     &                              hisVid(idGhat(itrc),ng))
            END IF
          END DO
# endif
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idVvis))) THEN
            got_var(idVvis)=.TRUE.
            status=nf90_inq_varid(ncHISid(ng),TRIM(Vname(1,idVvis)),    &
     &                            hisVid(idVvis,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idTdif))) THEN
            got_var(idTdif)=.TRUE.
            status=nf90_inq_varid(ncHISid(ng),TRIM(Vname(1,idTdif)),    &
     &                            hisVid(idTdif,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idSdif))) THEN
            got_var(idSdif)=.TRUE.
            status=nf90_inq_varid(ncHISid(ng),TRIM(Vname(1,idSdif)),    &
     &                            hisVid(idSdif,ng))
          END IF
# if defined GLS_MIXING || defined MY25_MIXING
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idMtke))) THEN
            got_var(idMtke)=.TRUE.
            status=nf90_inq_varid(ncHISid(ng),TRIM(Vname(1,idMtke)),    &
     &                            hisVid(idMtke,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idMtls))) THEN
            got_var(idMtls)=.TRUE.
            status=nf90_inq_varid(ncHISid(ng),TRIM(Vname(1,idMtls)),    &
     &                            hisVid(idMtls,ng))
          END IF
# endif
          DO itrc=1,NAT
            IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idTsur(itrc)))) THEN
              got_var(idTsur(itemp))=.TRUE.
              status=nf90_inq_varid(ncHISid(ng),                        &
     &                              TRIM(Vname(1,idTsur(itrc))),        &
     &                              hisVid(idTsur(itrc),ng))
            END IF
          END DO
# ifdef BULK_FLUXES
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idLhea))) THEN
            got_var(idLhea)=.TRUE.
            status=nf90_inq_varid(ncHISid(ng),TRIM(Vname(1,idLhea)),    &
     &                            hisVid(idLhea,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idShea))) THEN
            got_var(idShea)=.TRUE.
            status=nf90_inq_varid(ncHISid(ng),TRIM(Vname(1,idShea)),    &
     &                            hisVid(idShea,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idLrad))) THEN
            got_var(idLrad)=.TRUE.
            status=nf90_inq_varid(ncHISid(ng),TRIM(Vname(1,idLrad)),    &
     &                            hisVid(idLrad,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idTair))) THEN
            got_var(idTair)=.TRUE.
            status=nf90_inq_varid(ncHISid(ng),TRIM(Vname(1,idTair)),    &
     &                            hisVid(idTair,ng))
          END IF
#  ifdef EMINUSP
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idevap))) THEN
            got_var(idevap)=.TRUE.
            status=nf90_inq_varid(ncHISid(ng),TRIM(Vname(1,idevap)),    &
     &                            hisVid(idevap,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idrain))) THEN
            got_var(idrain)=.TRUE.
            status=nf90_inq_varid(ncHISid(ng),TRIM(Vname(1,idrain)),    &
     &                            hisVid(idrain,ng))
          END IF
#  endif
# endif
# ifdef SHORTWAVE
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idSrad))) THEN
            got_var(idSrad)=.TRUE.
            status=nf90_inq_varid(ncHISid(ng),TRIM(Vname(1,idSrad)),    &
     &                            hisVid(idSrad,ng))
          END IF
# endif
#endif
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idUsms))) THEN
            got_var(idUsms)=.TRUE.
            status=nf90_inq_varid(ncHISid(ng),TRIM(Vname(1,idUsms)),    &
     &                            hisVid(idUsms,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idVsms))) THEN
            got_var(idVsms)=.TRUE.
            status=nf90_inq_varid(ncHISid(ng),TRIM(Vname(1,idVsms)),    &
     &                            hisVid(idVsms,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idUbms))) THEN
            got_var(idUbms)=.TRUE.
            status=nf90_inq_varid(ncHISid(ng),TRIM(Vname(1,idUbms)),    &
     &                            hisVid(idUbms,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idVbms))) THEN
            got_var(idVbms)=.TRUE.
            status=nf90_inq_varid(ncHISid(ng),TRIM(Vname(1,idVbms)),    &
     &                            hisVid(idVbms,ng))
          END IF
#ifdef SOLVE3D
# ifdef BBL_MODEL
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idUbrs))) THEN
            got_var(idUbrs)=.TRUE.
            status=nf90_inq_varid(ncHISid(ng),TRIM(Vname(1,idUbrs)),    &
     &                            hisVid(idUbrs,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idVbrs))) THEN
            got_var(idVbrs)=.TRUE.
            status=nf90_inq_varid(ncHISid(ng),TRIM(Vname(1,idVbrs)),    &
     &                            hisVid(idVbrs,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idUbws))) THEN
            got_var(idUbws)=.TRUE.
            status=nf90_inq_varid(ncHISid(ng),TRIM(Vname(1,idUbws)),    &
     &                            hisVid(idUbws,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idVbws))) THEN
            got_var(idVbws)=.TRUE.
            status=nf90_inq_varid(ncHISid(ng),TRIM(Vname(1,idVbws)),    &
     &                            hisVid(idVbws,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idUbcs))) THEN
            got_var(idUbcs)=.TRUE.
            status=nf90_inq_varid(ncHISid(ng),TRIM(Vname(1,idUbcs)),    &
     &                            hisVid(idUbcs,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idVbcs))) THEN
            got_var(idVbcs)=.TRUE.
            status=nf90_inq_varid(ncHISid(ng),TRIM(Vname(1,idVbcs)),    &
     &                            hisVid(idVbcs,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idUbot))) THEN
            got_var(idUbot)=.TRUE.
            status=nf90_inq_varid(ncHISid(ng),TRIM(Vname(1,idUbot)),    &
     &                            hisVid(idUbot,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idVbot))) THEN
            got_var(idVbot)=.TRUE.
            status=nf90_inq_varid(ncHISid(ng),TRIM(Vname(1,idVbot)),    &
     &                            hisVid(idVbot,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idUbur))) THEN
            got_var(idUbur)=.TRUE.
            status=nf90_inq_varid(ncHISid(ng),TRIM(Vname(1,idUbur)),    &
     &                            hisVid(idUbur,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idVbvr))) THEN
            got_var(idVbvr)=.TRUE.
            status=nf90_inq_varid(ncHISid(ng),TRIM(Vname(1,idVbvr)),    &
     &                            hisVid(idVbvr,ng))
          END IF
# endif
# ifdef SEDIMENT
          DO itrc=1,NST
#  ifdef BEDLOAD
            IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idUbld(itrc)))) THEN
              got_var(idUbld(itrc))=.true.
              status=nf90_inq_varid(ncHISid(ng),                        &
     &                              TRIM(Vname(1,idUbld(itrc))),        &
     &                              hisVid(idUbld(itrc),ng))
            END IF
            IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idVbld(itrc)))) THEN
              got_var(idVbld(itrc))=.true.
              status=nf90_inq_varid(ncHISid(ng),                        &
     &                              TRIM(Vname(1,idVbld(itrc))),        &
     &                              hisVid(idVbld(itrc),ng))
            END IF
#  endif
            IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idfrac(itrc)))) THEN
              got_var(idfrac(itrc))=.TRUE.
              status=nf90_inq_varid(ncHISid(ng),                        &
     &                              TRIM(Vname(1,idfrac(itrc))),        &
     &                              hisVid(idfrac(itrc),ng))
            END IF
            IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idBmas(itrc)))) THEN
              got_var(idBmas(itrc))=.TRUE.
              status=nf90_inq_varid(ncHISid(ng),                        &
     &                              TRIM(Vname(1,idBmas(itrc))),        &
     &                              hisVid(idBmas(itrc),ng))
            END IF
          END DO
          DO itrc=1,MBEDP
            IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idSbed(itrc)))) THEN
              got_var(idSbed(itrc))=.TRUE.
              status=nf90_inq_varid(ncHISid(ng),                        &
     &                              TRIM(Vname(1,idSbed(itrc))),        &
     &                              hisVid(idSbed(itrc),ng))
            END IF
          END DO
# endif
# if defined SEDIMENT || defined BBL_MODEL
          DO itrc=1,MBOTP
            IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idBott(itrc)))) THEN
              got_var(idBott(itrc))=.TRUE.
              status=nf90_inq_varid(ncHISid(ng),                        &
     &                              TRIM(Vname(1,idBott(itrc))),        &
     &                              hisVid(idBott(itrc),ng))
            END IF
          END DO
# endif
#endif
#if defined WET_DRY
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idWmsk))) THEN
            got_var(idWmsk)=.TRUE.
            status=nf90_inq_varid(ncHISid(ng),                          &
     &                            TRIM(Vname(1,idWmsk)),                &
     &                            hisVid(idWmsk,ng))
          END IF
#endif
#if defined NEARSHORE_MELLOR
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idW2xx))) THEN
            got_var(idW2xx)=.TRUE.
            status=nf90_inq_varid(ncHISid(ng),                          &
     &                            TRIM(Vname(1,idW2xx)),                &
     &                            hisVid(idW2xx,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idW2xy))) THEN
            got_var(idW2xy)=.TRUE.
            status=nf90_inq_varid(ncHISid(ng),                          &
     &                            TRIM(Vname(1,idW2xy)),                &
     &                            hisVid(idW2xy,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idW2yy))) THEN
            got_var(idW2yy)=.TRUE.
            status=nf90_inq_varid(ncHISid(ng),                          &
     &                            TRIM(Vname(1,idW2yy)),                &
     &                            hisVid(idW2yy,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idU2rs))) THEN
            got_var(idU2rs)=.TRUE.
            status=nf90_inq_varid(ncHISid(ng),                          &
     &                            TRIM(Vname(1,idU2rs)),                &
     &                            hisVid(idU2rs,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idV2rs))) THEN
            got_var(idV2rs)=.TRUE.
            status=nf90_inq_varid(ncHISid(ng),                          &
     &                            TRIM(Vname(1,idV2rs)),                &
     &                            hisVid(idV2rs,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idU2Sd))) THEN
            got_var(idU2Sd)=.TRUE.
            status=nf90_inq_varid(ncHISid(ng),                          &
     &                            TRIM(Vname(1,idU2Sd)),                &
     &                            hisVid(idU2Sd,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idV2Sd))) THEN
            got_var(idV2Sd)=.TRUE.
            status=nf90_inq_varid(ncHISid(ng),                          &
     &                            TRIM(Vname(1,idV2Sd)),                &
     &                            hisVid(idV2Sd,ng))
          END IF
# ifdef SOLVE3D
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idW3xx))) THEN
            got_var(idW3xx)=.TRUE.
            status=nf90_inq_varid(ncHISid(ng),                          &
     &                            TRIM(Vname(1,idW3xx)),                &
     &                            hisVid(idW3xx,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idW3xy))) THEN
            got_var(idW3xy)=.TRUE.
            status=nf90_inq_varid(ncHISid(ng),                          &
     &                            TRIM(Vname(1,idW3xy)),                &
     &                            hisVid(idW3xy,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idW3yy))) THEN
            got_var(idW3yy)=.TRUE.
            status=nf90_inq_varid(ncHISid(ng),                          &
     &                            TRIM(Vname(1,idW3yy)),                &
     &                            hisVid(idW3yy,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idW3zx))) THEN
            got_var(idW3zx)=.TRUE.
            status=nf90_inq_varid(ncHISid(ng),                          &
     &                            TRIM(Vname(1,idW3zx)),                &
     &                            hisVid(idW3zx,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idW3zy))) THEN
            got_var(idW3zy)=.TRUE.
            status=nf90_inq_varid(ncHISid(ng),                          &
     &                            TRIM(Vname(1,idW3zy)),                &
     &                            hisVid(idW3zy,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idU3rs))) THEN
            got_var(idU3rs)=.TRUE.
            status=nf90_inq_varid(ncHISid(ng),                          &
     &                            TRIM(Vname(1,idU3rs)),                &
     &                            hisVid(idU3rs,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idV3rs))) THEN
            got_var(idV3rs)=.TRUE.
            status=nf90_inq_varid(ncHISid(ng),                          &
     &                            TRIM(Vname(1,idV3rs)),                &
     &                            hisVid(idV3rs,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idU3Sd))) THEN
            got_var(idU3Sd)=.TRUE.
            status=nf90_inq_varid(ncHISid(ng),                          &
     &                            TRIM(Vname(1,idU3Sd)),                &
     &                            hisVid(idU3Sd,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idV3Sd))) THEN
            got_var(idV3Sd)=.TRUE.
            status=nf90_inq_varid(ncHISid(ng),                          &
     &                            TRIM(Vname(1,idV3Sd)),                &
     &                            hisVid(idV3Sd,ng))
          END IF
# endif
#endif
#ifdef WAVES_HEIGHT
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idWamp))) THEN
            got_var(idWamp)=.TRUE.
            status=nf90_inq_varid(ncHISid(ng),                          &
     &                            TRIM(Vname(1,idWamp)),                &
     &                            hisVid(idWamp,ng))
          END IF
#endif
#ifdef WAVES_LENGTH
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idWlen))) THEN
            got_var(idWlen)=.TRUE.
            status=nf90_inq_varid(ncHISid(ng),                          &
     &                            TRIM(Vname(1,idWlen)),                &
     &                            hisVid(idWlen,ng))
          END IF
#endif
#ifdef WAVES_DIR
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idWdir))) THEN
            got_var(idWdir)=.TRUE.
            status=nf90_inq_varid(ncHISid(ng),                          &
     &                            TRIM(Vname(1,idWdir)),                &
     &                            hisVid(idWdir,ng))
          END IF
#endif
        END DO
!
!  Check if history variables are available in input NetCDF file.
!
        IF (.not.got_var(idtime)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idtime)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#if defined SEDIMENT && defined SED_MORPH
        IF (.not.got_var(idbath).and.Hout(idbath,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idbath)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#endif
        IF (.not.got_var(idFsur).and.Hout(idFsur,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idFsur)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUbar).and.Hout(idUbar,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idUbar)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVbar).and.Hout(idVbar,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idVbar)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#ifdef SOLVE3D
        IF (.not.got_var(idUvel).and.Hout(idUvel,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idUvel)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVvel).and.Hout(idVvel,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idVvel)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idWvel).and.Hout(idWvel,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idWvel)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idOvel).and.Hout(idOvel,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idOvel)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        DO itrc=1,NT(ng)
          IF (.not.got_var(idTvar(itrc)).and.Hout(idTvar(itrc),ng)) THEN
            WRITE (stdout,60) TRIM(Vname(1,idTvar(itrc))), TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
        IF (.not.got_var(idDano).and.Hout(idDano,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idDano)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# ifdef LMD_SKPP
        IF (.not.got_var(idHsbl).and.Hout(idHsbl,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idHsbl)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
# ifdef LMD_BKPP
        IF (.not.got_var(idHbbl).and.Hout(idHbbl,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idHbbl)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
# if defined FORWARD_WRITE && defined LMD_NONLOCAL
        DO itrc=1,NAT
          IF (.not.got_var(idGhat(itrc)).and.Hout(idGhat(itrc),ng)) THEN
            WRITE (stdout,60) TRIM(Vname(1,idGhat(itrc))), TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
# endif
        IF (.not.got_var(idVvis).and.Hout(idVvis,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idVvis)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idTdif).and.Hout(idTdif,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idTdif)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# ifdef SALINITY
        IF (.not.got_var(idSdif).and.Hout(idSdif,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idSdif)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
# if defined GLS_MIXING || defined MY25_MIXING
        IF (.not.got_var(idMtke).and.Hout(idMtke,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idMtke)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idMtls).and.Hout(idMtls,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idMtls)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
        DO itrc=1,NAT
          IF (.not.got_var(idTsur(itrc)).and.Hout(idTsur(itrc),ng)) THEN
            WRITE (stdout,60) TRIM(Vname(1,idTsur(itrc))), TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
# ifdef BULK_FLUXES
        IF (.not.got_var(idLhea).and.Hout(idLhea,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idLhea)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idShea).and.Hout(idShea,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idShea)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idLrad).and.Hout(idLrad,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idLrad)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idTair).and.Hout(idTair,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idTair)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#    ifdef EMINUSP
        IF (.not.got_var(idevap).and.Hout(idevap,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idevap)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idrain).and.Hout(idrain,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idrain)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#    endif
# endif
# ifdef SHORTWAVE
        IF (.not.got_var(idSrad).and.Hout(idSrad,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idSrad)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
#endif
        IF (.not.got_var(idUsms).and.Hout(idUsms,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idUsms)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVsms).and.Hout(idVsms,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idVsms)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUbms).and.Hout(idUbms,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idUbms)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVbms).and.Hout(idVbms,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idVbms)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#ifdef SOLVE3D
# ifdef BBL_MODEL
        IF (.not.got_var(idUbrs).and.Hout(idUbrs,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idUbrs)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVbrs).and.Hout(idVbrs,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idVbrs)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUbws).and.Hout(idUbws,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idUbws)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVbws).and.Hout(idVbws,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idVbws)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUbcs).and.Hout(idUbcs,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idUbcs)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVbcs).and.Hout(idVbcs,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idVbcs)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUbot).and.Hout(idUbot,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idUbot)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVbot).and.Hout(idVbot,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idVbot)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUbur).and.Hout(idUbur,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idUbur)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVbvr).and.Hout(idVbvr,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idVbvr)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
# ifdef SEDIMENT
        DO i=1,NST
#  ifdef BEDLOAD
          IF (.not.got_var(idUbld(i)).and.Hout(idUbld(i),ng)) THEN
            WRITE (stdout,60) TRIM(Vname(1,idUbld(i))), TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
          IF (.not.got_var(idVbld(i)).and.Hout(idVbld(i),ng)) THEN
            WRITE (stdout,60) TRIM(Vname(1,idVbld(i))), TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
#  endif
          IF (.not.got_var(idfrac(i)).and.Hout(idfrac(i),ng)) THEN
            WRITE (stdout,60) TRIM(Vname(1,idfrac(i))), TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
          IF(.not.got_var(idBmas(i)).and.Hout(idBmas(i),ng)) THEN
            WRITE (stdout,60) TRIM(Vname(1,idBmas(i))), TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
        DO i=1,MBEDP
          IF (.not.got_var(idSbed(i)).and.Hout(idSbed(i),ng)) THEN
            WRITE (stdout,60) TRIM(Vname(1,idSbed(i))), TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
# endif
# if defined SEDIMENT || defined BBL_MODEL
        DO i=1,MBOTP
          IF (.not.got_var(idBott(i)).and.Hout(idBott(i),ng)) THEN
            WRITE (stdout,60) TRIM(Vname(1,idBott(i))), TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
# endif
#endif
#if defined WET_DRY
!        IF (.not.got_var(idWmsk).and.Hout(idWmsk,ng)) THEN
        IF (.not.got_var(idWmsk)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idWmsk)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#endif
#if defined NEARSHORE_MELLOR
        IF (.not.got_var(idW2xx).and.Hout(idW2xx,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idW2xx)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idW2xy).and.Hout(idW2xy,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idW2xy)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idW2yy).and.Hout(idW2yy,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idW2yy)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idU2rs).and.Hout(idU2rs,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idU2rs)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idV2rs).and.Hout(idV2rs,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idV2rs)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idU2Sd).and.Hout(idU2Sd,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idU2Sd)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idV2Sd).and.Hout(idV2Sd,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idV2Sd)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# ifdef SOLVE3D
        IF (.not.got_var(idW3xx).and.Hout(idW3xx,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idW3xx)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idW3xy).and.Hout(idW3xy,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idW3xy)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idW3yy).and.Hout(idW3yy,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idW3yy)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idW3zx).and.Hout(idW3zx,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idW3zx)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idW3zy).and.Hout(idW3zy,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idW3zy)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idU3Sd).and.Hout(idU3rs,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idU3rs)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idV3rs).and.Hout(idV3rs,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idV3rs)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idU3Sd).and.Hout(idU3Sd,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idU3Sd)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idV3Sd).and.Hout(idV3Sd,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idV3Sd)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
#endif
#ifdef WAVES_HEIGHT
        IF (.not.got_var(idWamp).and.Hout(idWamp,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idWamp)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#endif
#ifdef WAVES_LENGTH
        IF (.not.got_var(idWlen).and.Hout(idWlen,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idWlen)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#endif
#ifdef WAVES_DIR
        IF (.not.got_var(idWdir).and.Hout(idWdir,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idWdir)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#endif
!
!  Set unlimited time record dimension to the appropriate value.
!
!!      tHISindx(ng)=tsize
        IF (ndefHIS(ng).gt.0) THEN
          tHISindx(ng)=((ntstart(ng)-1)-                                &
     &                  ndefHIS(ng)*((ntstart(ng)-1)/ndefHIS(ng)))/     &
     &                 nHIS(ng)
        ELSE
          tHISindx(ng)=(ntstart(ng)-1)/nHIS(ng)
        END IF
        tHISindx(ng)=MIN(tHISindx(ng),tsize)
      END IF
!
  10  FORMAT (6x,'DEF_HIS   - creating history file: ',a)
  20  FORMAT (6x,'DEF_HIS   - inquiring history file: ',a)
  30  FORMAT (/,' DEF_HIS - unable to create history NetCDF file: ',a)
  40  FORMAT (1pe11.4,1x,'millimeter')
  50  FORMAT (/,' DEF_HIS - unable to open history NetCDF file: ',a)
  60  FORMAT (/,' DEF_HIS - unable to find variable: ',a,2x,            &
     &        ' in history NetCDF file: ',a)
      RETURN
      END SUBROUTINE def_his
