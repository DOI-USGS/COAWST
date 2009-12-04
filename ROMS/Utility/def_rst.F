#include "cppdefs.h"
      SUBROUTINE def_rst (ng)
!
!svn $Id: def_rst.F 733 2008-09-07 01:56:45Z jcwarner $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine creates restart NetCDF file, it defines its            !
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
!
!  Local variable declarations.
!
      integer, parameter :: Natt = 24

      logical :: Ldefine, got_var(NV)

      integer :: i, j, nrec, nvd, nvd3, nvd4, nvd5
      integer :: recdim, status, varid, ofst

      integer :: DimIDs(31)
      integer :: r2dgrd(4), ru2dgrd(4), rv2dgrd(4)
      integer :: sr2dgrd(3), su2dgrd(3), sv2dgrd(3)
      integer :: sr3dgrd(4), su3dgrd(4), sv3dgrd(4)
      integer :: t2dgrd(4), u2dgrd(4), v2dgrd(4)
      integer :: Vsize(4)

      integer :: def_info, def_var, wrt_info

#ifdef SOLVE3D
      integer :: itrc

      integer :: k3dgrd(5), t3dgrd(5)
      integer :: r3dgrd(4), ru3dgrd(5), rv3dgrd(5) 
      integer :: u3dgrd(5), v3dgrd(5), w3dgrd(4)
#endif

      real(r8) :: Aval(6)

      character (len=80) :: Vinfo(Natt)
      character (len=80) :: fname, ncname
!
!=======================================================================
!  Create a new restart NetCDF file.
!=======================================================================
!
!  Activate creation of restart NetCDF file.  Create a new restart
!  file if during a restart run, the restart filename "RSTname"
!  is different than the initial filename "INIname".
!
      IF (exit_flag.ne.NoError) RETURN
      ncname=RSTname(ng)
      Ldefine=.FALSE.
      IF (((nrrec(ng).eq.0).and.(iic(ng).eq.ntstart(ng))).or.           &
     &    ((nrrec(ng).ne.0).and.                                        &
     &     (TRIM(ncname).ne.TRIM(INIname(ng))))) THEN
        Ldefine=.TRUE.
      END IF
!
      IF (Ldefine.and.OutThread) THEN
        status=nf90_create(TRIM(ncname),nf90_clobber,ncRSTid(ng))
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,10) TRIM(ncname)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#ifdef DISTRIBUTE
      IF (Ldefine) THEN
        CALL mp_bcasti (ng, iNLM, ncRSTid(ng), 1)
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
      IF (Ldefine.and.OutThread) THEN
        status=nf90_def_dim(ncRSTid(ng),'xi_rho'  ,Lm(ng)+ofst+1,       &
     &                      DimIDs( 1))
        status=nf90_def_dim(ncRSTid(ng),'xi_u'    ,Lm(ng)+ofst,         &
     &                      DimIDs( 2))
        status=nf90_def_dim(ncRSTid(ng),'xi_v'    ,Lm(ng)+ofst+1,       &
     &                      DimIDs( 3))
        status=nf90_def_dim(ncRSTid(ng),'xi_psi'  ,Lm(ng)+ofst,         &
     &                      DimIDs( 4))
        status=nf90_def_dim(ncRSTid(ng),'eta_rho' ,Mm(ng)+ofst+1,       &
     &                      DimIDs( 5))
        status=nf90_def_dim(ncRSTid(ng),'eta_u'   ,Mm(ng)+ofst+1,       &
     &                      DimIDs( 6))
        status=nf90_def_dim(ncRSTid(ng),'eta_v'   ,Mm(ng)+ofst,         &
     &                      DimIDs( 7))
        status=nf90_def_dim(ncRSTid(ng),'eta_psi' ,Mm(ng)+ofst,         &
     &                      DimIDs( 8))
#if defined WRITE_WATER && defined MASKING
        status=nf90_def_dim(ncRSTid(ng),'xy_rho'  ,Nxyr(ng),            &
     &                      DimIDs(17))
        status=nf90_def_dim(ncRSTid(ng),'xy_u'    ,Nxyu(ng),            &
     &                      DimIDs(18))
        status=nf90_def_dim(ncRSTid(ng),'xy_v'    ,Nxyv(ng),            &
     &                      DimIDs(19))
#endif
#ifdef SOLVE3D
# if defined WRITE_WATER && defined MASKING
        status=nf90_def_dim(ncRSTid(ng),'xyz_rho' ,Nxyr(ng)*N(ng),      &
     &                      DimIDs(20))
        status=nf90_def_dim(ncRSTid(ng),'xyz_u'   ,Nxyu(ng)*N(ng),      &
     &                      DimIDs(21))
        status=nf90_def_dim(ncRSTid(ng),'xyz_v'   ,Nxyv(ng)*N(ng),      &
     &                      DimIDs(22))
        status=nf90_def_dim(ncRSTid(ng),'xyz_w'   ,Nxyr(ng)*(N(ng)+1),  &
     &                      DimIDs(23))
# endif
        status=nf90_def_dim(ncRSTid(ng),'s_rho'   ,N(ng),               &
     &                      DimIDs( 9))
        status=nf90_def_dim(ncRSTid(ng),'s_w'     ,N(ng)+1,             &
     &                      DimIDs(10))
        status=nf90_def_dim(ncRSTid(ng),'tracer'  ,NT(ng),              &
     &                      DimIDs(11))
# ifdef SEDIMENT
        status=nf90_def_dim(ncRSTid(ng),'Nbed'    ,Nbed,                &
     &                      DimIDs(16))
#  if defined WRITE_WATER && defined MASKING
        status=nf90_def_dim(ncRSTid(ng),'xybed'   ,Nxyr(ng)*Nbed,       &
     &                      DimIDs(24))
#  endif
# endif
# ifdef ECOSIM
        status=nf90_def_dim(ncRSTid(ng),'Nphy'    ,Nphy,                &
     &                      DimIDs(25))
        status=nf90_def_dim(ncRSTid(ng),'Nbac'    ,Nbac,                &
     &                      DimIDs(26))
        status=nf90_def_dim(ncRSTid(ng),'Ndom'    ,Ndom,                &
     &                      DimIDs(27))
        status=nf90_def_dim(ncRSTid(ng),'Nfec'    ,Nfec,                &
     &                      DimIDs(28))
# endif
#endif
        status=nf90_def_dim(ncRSTid(ng),'boundary',4,                   &
     &                      DimIDs(14))
#ifdef FOUR_DVAR
        status=nf90_def_dim(ncRSTid(ng),'Nstate  ',NstateVar(ng),       &
     &                      DimIDs(29))
#endif
#ifdef PERFECT_RESTART
        status=nf90_def_dim(ncRSTid(ng),'two     ',2            ,       &
     &                      DimIDs(30))
        status=nf90_def_dim(ncRSTid(ng),'three   ',3            ,       &
     &                      DimIDs(31))
#endif
        status=nf90_def_dim(ncRSTid(ng),TRIM(ADJUSTL(Vname(5,idtime))), &
     &                      nf90_unlimited,DimIDs(12))
        recdim=DimIDs(12)
!
!  Set number of dimensions for output variables.
!
#if defined WRITE_WATER && defined MASKING
        nvd3=2
        nvd4=2
        nvd5=2
#else
        nvd3=3
        nvd4=4
        nvd5=5
#endif
!
!  Define dimension vectors for staggered tracer type variables.
!
#if defined WRITE_WATER && defined MASKING
        t2dgrd(1)=DimIDs(17)
        t2dgrd(2)=DimIDs(12)
        sr2dgrd(1)=DimIDs(17)
        sr2dgrd(2)=DimIDs(12)
# ifdef SOLVE3D
        t3dgrd(1)=DimIDs(20)
        t3dgrd(2)=DimIDs(12)
# endif
#else
        t2dgrd(1)=DimIDs( 1)
        t2dgrd(2)=DimIDs( 5)
        sr2dgrd(1)=DimIDs( 1)
        sr2dgrd(2)=DimIDs( 5)
        sr2dgrd(3)=DimIDs(12)
# ifdef PERFECT_RESTART
        t2dgrd(3)=DimIDs(31)
        t2dgrd(4)=DimIDs(12)
# else
        t2dgrd(3)=DimIDs(12)
# endif
# ifdef SOLVE3D
        t3dgrd(1)=DimIDs( 1)
        t3dgrd(2)=DimIDs( 5)
        t3dgrd(3)=DimIDs( 9)
        r3dgrd(1)=DimIDs( 1)
        r3dgrd(2)=DimIDs( 5)
        r3dgrd(3)=DimIDs( 9)
#  ifdef PERFECT_RESTART
        t3dgrd(4)=DimIDs(30)
        t3dgrd(5)=DimIDs(12)
#  else
        t3dgrd(4)=DimIDs(12)
#  endif
        r3dgrd(4)=DimIDs(12)
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
# ifdef PERFECT_RESTART
        u2dgrd(3)=DimIDs(31)
        u2dgrd(4)=DimIDs(12)
# else
        u2dgrd(3)=DimIDs(12)
# endif
# ifdef SOLVE3D
        u3dgrd(1)=DimIDs( 2)
        u3dgrd(2)=DimIDs( 6)
        u3dgrd(3)=DimIDs( 9)
#  ifdef PERFECT_RESTART
        u3dgrd(4)=DimIDs(30)
        u3dgrd(5)=DimIDs(12)
#  else
        u3dgrd(4)=DimIDs(12)
#  endif
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
# ifdef PERFECT_RESTART
        v2dgrd(3)=DimIDs(31)
        v2dgrd(4)=DimIDs(12)
# else
        v2dgrd(3)=DimIDs(12)
# endif
# ifdef SOLVE3D
        v3dgrd(1)=DimIDs( 3)
        v3dgrd(2)=DimIDs( 7)
        v3dgrd(3)=DimIDs( 9)
#  ifdef PERFECT_RESTART
        v3dgrd(4)=DimIDs(30)
        v3dgrd(5)=DimIDs(12)
#  else
        v3dgrd(4)=DimIDs(12)
#  endif
# endif
#endif
#ifdef PERFECT_RESTART
!
!  Define dimension vectors for RHS free-surface equation.
!
# if defined WRITE_WATER && defined MASKING
        r2dgrd(1)=DimIDs(17)
        r2dgrd(2)=DimIDs(12)
# else
        r2dgrd(1)=DimIDs( 1)
        r2dgrd(2)=DimIDs( 5)
        r2dgrd(3)=DimIDs(30)
        r2dgrd(4)=DimIDs(12)
# endif
!
!  Define dimension vectors for RHS u-momemtum equation.
!
# if defined WRITE_WATER && defined MASKING
        ru2dgrd(1)=DimIDs(18)
        ru2dgrd(2)=DimIDs(12)
#  ifdef SOLVE3D
        ru2dgrd(1)=DimIDs(21)
        ru2dgrd(2)=DimIDs(12)
#  endif
# else
        ru2dgrd(1)=DimIDs( 2)
        ru2dgrd(2)=DimIDs( 6)
        ru2dgrd(3)=DimIDs(30)
        ru2dgrd(4)=DimIDs(12)
#  ifdef SOLVE3D
        ru3dgrd(1)=DimIDs( 2)
        ru3dgrd(2)=DimIDs( 6)
        ru3dgrd(3)=DimIDs(10)
        ru3dgrd(4)=DimIDs(30)
        ru3dgrd(5)=DimIDs(12)
#  endif
# endif
!
!  Define dimension vectors for RHS v-momemtum equation.
!
# if defined WRITE_WATER && defined MASKING
        rv2dgrd(1)=DimIDs(19)
        rv2dgrd(2)=DimIDs(12)
#  ifdef SOLVE3D
        rv3dgrd(1)=DimIDs(22)
        rv3dgrd(2)=DimIDs(12)
#  endif
# else
        rv2dgrd(1)=DimIDs( 3)
        rv2dgrd(2)=DimIDs( 7)
        rv2dgrd(3)=DimIDs(30)
        rv2dgrd(4)=DimIDs(12)
#  ifdef SOLVE3D
        rv3dgrd(1)=DimIDs( 3)
        rv3dgrd(2)=DimIDs( 7)
        rv3dgrd(3)=DimIDs(10)
        rv3dgrd(4)=DimIDs(30)
        rv3dgrd(5)=DimIDs(12)
#  endif
# endif
#endif
#ifdef SOLVE3D
!
!  Define dimension vector for staggered w-momemtum type variables.
!
# if defined WRITE_WATER && defined MASKING
        w3dgrd(1)=DimIDs(23)
        w3dgrd(2)=DimIDs(12)
#  ifdef PERFECT_RESTART
        k3dgrd(1)=DimIDs(23)
        k3dgrd(2)=DimIDs(12)
#  endif
# else
        w3dgrd(1)=DimIDs( 1)
        w3dgrd(2)=DimIDs( 5)
        w3dgrd(3)=DimIDs(10)
        w3dgrd(4)=DimIDs(12)
#  ifdef PERFECT_RESTART
        k3dgrd(1)=DimIDs( 1)
        k3dgrd(2)=DimIDs( 5)
        k3dgrd(3)=DimIDs(10)
        k3dgrd(4)=DimIDs(30)
        k3dgrd(5)=DimIDs(12)
#  endif
# endif
!
!  Define dimension vector for sediment, radiation stress variables.
!
# if defined WRITE_WATER && defined MASKING
        sr3dgrd(1)=DimIDs(24)
        sr3dgrd(2)=DimIDs(12)
        su3dgrd(1)=DimIDs(21)
        su3dgrd(2)=DimIDs(12)
        sv3dgrd(1)=DimIDs(22)
        sv3dgrd(2)=DimIDs(12)
# else
        sr3dgrd(1)=DimIDs( 1)
        sr3dgrd(2)=DimIDs( 5)
        sr3dgrd(3)=DimIDs(16)
        sr3dgrd(4)=DimIDs(12)
        su3dgrd(1)=DimIDs( 2)
        su3dgrd(2)=DimIDs( 6)
        su3dgrd(3)=DimIDs( 9)
        su3dgrd(4)=DimIDs(12)
        sv3dgrd(1)=DimIDs( 3)
        sv3dgrd(2)=DimIDs( 7)
        sv3dgrd(3)=DimIDs( 9)
        sv3dgrd(4)=DimIDs(12)
# endif
#endif
#if defined WRITE_WATER && defined MASKING
        su2dgrd(1)=DimIDs(18)
        su2dgrd(2)=DimIDs(12)
        sv2dgrd(1)=DimIDs(19)
        sv2dgrd(2)=DimIDs(12)
#else
        su2dgrd(1)=DimIDs( 2)
        su2dgrd(2)=DimIDs( 6)
        su2dgrd(3)=DimIDs(12)
        sv2dgrd(1)=DimIDs( 3)
        sv2dgrd(2)=DimIDs( 7)
        sv2dgrd(3)=DimIDs(12)
#endif
!
!  Initialize unlimited time record dimension.
!
        tRSTindx(ng)=0
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
!-----------------------------------------------------------------------
!  Define time-recordless information variables.
!-----------------------------------------------------------------------
!
        status=def_info(ng,ncRSTid(ng),ncname,DimIDs)
        IF (exit_flag.ne.NoError) RETURN
!
!-----------------------------------------------------------------------
!  Define variables and their attributes.
!-----------------------------------------------------------------------

#ifdef PERFECT_RESTART
!
!  Define time-stepping indices.
!
# ifdef SOLVE3D
        Vinfo( 1)='nstp'
        Vinfo( 2)='3D equations time level index, nstp'
        status=def_var(ncRSTid(ng),varid,nf90_int,1,recdim,             &
     &                 Aval,Vinfo,ncname)

        Vinfo( 1)='nrhs'
        Vinfo( 2)='3D equations time level index, nrhs'
        status=def_var(ncRSTid(ng),varid,nf90_int,1,recdim,             &
     &                 Aval,Vinfo,ncname)

        Vinfo( 1)='nnew'
        Vinfo( 2)='3D equations time level index, nnew'
        status=def_var(ncRSTid(ng),varid,nf90_int,1,recdim,             &
     &                 Aval,Vinfo,ncname)
# endif
        Vinfo( 1)='kstp'
        Vinfo( 2)='3D equations time level index, kstp'
        status=def_var(ncRSTid(ng),varid,nf90_int,1,recdim,             &
     &                 Aval,Vinfo,ncname)

        Vinfo( 1)='krhs'
        Vinfo( 2)='3D equations time level index, krhs'
        status=def_var(ncRSTid(ng),varid,nf90_int,1,recdim,             &
     &                 Aval,Vinfo,ncname)

        Vinfo( 1)='knew'
        Vinfo( 2)='3D equations time level index, knew'
        status=def_var(ncRSTid(ng),varid,nf90_int,1,recdim,             &
     &                 Aval,Vinfo,ncname)
#endif
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
        status=def_var(ncRSTid(ng),rstVid(idtime,ng),NF_TYPE,           &
     &                 1,recdim,Aval,Vinfo,ncname)
#if defined SEDIMENT && defined SED_MORPH
!
!  Define free-surface.
!
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
        status=def_var(ncRSTid(ng),rstVid(idbath,ng),NF_FRST,           &
     &                 nvd3,t2dgrd,Aval,Vinfo,ncname)
#endif
!
!  Define free-surface.
!
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
#ifdef PERFECT_RESTART
        status=def_var(ncRSTid(ng),rstVid(idFsur,ng),NF_FRST,           &
     &                 nvd4,t2dgrd,Aval,Vinfo,ncname)
#else
        status=def_var(ncRSTid(ng),rstVid(idFsur,ng),NF_FRST,           &
     &                 nvd3,t2dgrd,Aval,Vinfo,ncname)
#endif
#ifdef PERFECT_RESTART
!
!  Define RHS of free-surface equation.
!
        Vinfo( 1)=Vname(1,idRzet)
        Vinfo( 2)=Vname(2,idRzet)
        Vinfo( 3)=Vname(3,idRzet)
        Vinfo(14)=Vname(4,idRzet)
        Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_rho'
# endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idFsur,ng),r8)
        status=def_var(ncRSTid(ng),rstVid(idRzet,ng),NF_FRST,           &
     &                 nvd4,r2dgrd,Aval,Vinfo,ncname)
#endif
!
!  Define 2D momentum in the XI-direction.
!
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
#ifdef PERFECT_RESTART
        status=def_var(ncRSTid(ng),rstVid(idUbar,ng),NF_FRST,           &
     &                 nvd4,u2dgrd,Aval,Vinfo,ncname)
#else
        status=def_var(ncRSTid(ng),rstVid(idUbar,ng),NF_FRST,           &
     &                 nvd3,u2dgrd,Aval,Vinfo,ncname)
#endif
#ifdef PERFECT_RESTART
!
!  Define RHS of 2D momentum equation in the XI-direction.
!
        Vinfo( 1)=Vname(1,idRu2d)
        Vinfo( 2)=Vname(2,idRu2d)
        Vinfo( 3)=Vname(3,idRu2d)
        Vinfo(14)=Vname(4,idRu2d)
        Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_u'
# endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idRu2d,ng),r8)
        status=def_var(ncRSTid(ng),rstVid(idRu2d,ng),NF_FRST,           &
     &                 nvd4,ru2dgrd,Aval,Vinfo,ncname)
#endif
!
!  Define 2D momentum in the ETA-direction.
!
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
#ifdef PERFECT_RESTART
        status=def_var(ncRSTid(ng),rstVid(idVbar,ng),NF_FRST,           &
     &                 nvd4,v2dgrd,Aval,Vinfo,ncname)
#else
        status=def_var(ncRSTid(ng),rstVid(idVbar,ng),NF_FRST,           &
     &                 nvd3,v2dgrd,Aval,Vinfo,ncname)
#endif
#ifdef PERFECT_RESTART
!
!  Define RHS of 2D momentum equation in the ETA-direction.
!
        Vinfo( 1)=Vname(1,idRv2d)
        Vinfo( 2)=Vname(2,idRv2d)
        Vinfo( 3)=Vname(3,idRv2d)
        Vinfo(14)=Vname(4,idRv2d)
        Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_v'
# endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idRv2d,ng),r8)
        status=def_var(ncRSTid(ng),rstVid(idRv2d,ng),NF_FRST,           &
     &                 nvd4,rv2dgrd,Aval,Vinfo,ncname)
#endif
#ifdef SOLVE3D
!
!  Define 3D momentum component in the XI-direction.
!
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
# ifdef PERFECT_RESTART
        status=def_var(ncRSTid(ng),rstVid(idUvel,ng),NF_FRST,           &
     &                 nvd5,u3dgrd,Aval,Vinfo,ncname)
# else
        status=def_var(ncRSTid(ng),rstVid(idUvel,ng),NF_FRST,           &
     &                 nvd4,u3dgrd,Aval,Vinfo,ncname)
# endif
# ifdef PERFECT_RESTART
!
!  Define RHS of 3D momentum equation in the XI-direction.
!
        Vinfo( 1)=Vname(1,idRu3d)
        Vinfo( 2)=Vname(2,idRu3d)
        Vinfo( 3)=Vname(3,idRu3d)
        Vinfo(14)=Vname(4,idRu3d)
        Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_u'
#  endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idRu3d,ng),r8)
        status=def_var(ncRSTid(ng),rstVid(idRu3d,ng),NF_FRST,           &
     &                 nvd5,ru3dgrd,Aval,Vinfo,ncname)
# endif
!
!  Define 3D momentum component in the ETA-direction.
!
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
# ifdef PERFECT_RESTART
        status=def_var(ncRSTid(ng),rstVid(idVvel,ng),NF_FRST,           &
     &                 nvd5,v3dgrd,Aval,Vinfo,ncname)
# else
        status=def_var(ncRSTid(ng),rstVid(idVvel,ng),NF_FRST,           &
     &                 nvd4,v3dgrd,Aval,Vinfo,ncname)
# endif
# ifdef PERFECT_RESTART
!
!  Define RHS of 3D momentum equation in the ETA-direction.
!
        Vinfo( 1)=Vname(1,idRv3d)
        Vinfo( 2)=Vname(2,idRv3d)
        Vinfo( 3)=Vname(3,idRv3d)
        Vinfo(14)=Vname(4,idRv3d)
        Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_v'
#  endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idRv3d,ng),r8)
        status=def_var(ncRSTid(ng),rstVid(idRv3d,ng),NF_FRST,           &
     &                 nvd5,rv3dgrd,Aval,Vinfo,ncname)
# endif
!
!  Define tracer type variables.
!
        DO itrc=1,NT(ng)
          Vinfo( 1)=Vname(1,idTvar(itrc))
          Vinfo( 2)=Vname(2,idTvar(itrc))
          Vinfo( 3)=Vname(3,idTvar(itrc))
          Vinfo(14)=Vname(4,idTvar(itrc))
          Vinfo(16)=Vname(1,idtime)
# ifdef SEDIMENT
            DO i=1,NST
              IF (itrc.eq.idsed(i)) THEN
                WRITE (Vinfo(19),20) 1000.0_r8*Sd50(i,ng)
              END IF
            END DO
# endif
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
# endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(r3dvar,r8)
# ifdef PERFECT_RESTART
          status=def_var(ncRSTid(ng),rstTid(itrc,ng),NF_FRST,           &
     &                   nvd5,t3dgrd,Aval,Vinfo,ncname)
# else
          status=def_var(ncRSTid(ng),rstTid(itrc,ng),NF_FRST,           &
     &                   nvd4,t3dgrd,Aval,Vinfo,ncname)
# endif
        END DO
!
!  Define density anomaly.
!
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
        status=def_var(ncRSTid(ng),rstVid(idDano,ng),NF_FRST,           &
     &                 nvd4,r3dgrd,Aval,Vinfo,ncname)
# ifdef LMD_SKPP
!
!  Define depth of surface boundary layer.
!
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
        status=def_var(ncRSTid(ng),rstVid(idHsbl,ng),NF_FRST,           &
     &                 nvd3,t2dgrd,Aval,Vinfo,ncname)
# endif
# ifdef LMD_BKPP
!
!  Define depth of bottom boundary layer.
!
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
        status=def_var(ncRSTid(ng),rstVid(idHbbl,ng),NF_FRST,           &
     &                 nvd3,t2dgrd,Aval,Vinfo,ncname)
# endif
# if defined BVF_MIXING || defined LMD_MIXING || \
     defined GLS_MIXING || defined MY25_MIXING
!
!  Define vertical viscosity coefficient.
!
        Vinfo( 1)=Vname(1,idVvis)
        Vinfo( 2)=Vname(2,idVvis)
        Vinfo( 3)=Vname(3,idVvis)
        Vinfo(14)=Vname(4,idVvis)
        Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_rho'
#  endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idVvis,ng),r8)
        status=def_var(ncRSTid(ng),rstVid(idVvis,ng),NF_FRST,           &
     &                 nvd4,w3dgrd,Aval,Vinfo,ncname)
!
!  Define vertical diffusion coefficient for potential temperature.
!
        Vinfo( 1)=Vname(1,idTdif)
        Vinfo( 2)=Vname(2,idTdif)
        Vinfo( 3)=Vname(3,idTdif)
        Vinfo(14)=Vname(4,idTdif)
        Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_rho'
#  endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idTdif,ng),r8)
        status=def_var(ncRSTid(ng),rstVid(idTdif,ng),NF_FRST,           &
     &                 nvd4,w3dgrd,Aval,Vinfo,ncname)
#  ifdef SALINITY
!
!  Define vertical diffusion coefficient for salinity.
!
        Vinfo( 1)=Vname(1,idSdif)
        Vinfo( 2)=Vname(2,idSdif)
        Vinfo( 3)=Vname(3,idSdif)
        Vinfo(14)=Vname(4,idSdif)
        Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_rho'
#   endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idSdif,ng),r8)
        status=def_var(ncRSTid(ng),rstVid(idSdif,ng),NF_FRST,           &
     &                 nvd4,w3dgrd,Aval,Vinfo,ncname)
#  endif
# endif
# if defined PERFECT_RESTART && \
     (defined GLS_MIXING     || defined MY25_MIXING)
!
!  Define turbulent kinetic energy.
!
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
        status=def_var(ncRSTid(ng),rstVid(idMtke,ng),NF_FRST,           &
     &                 nvd5,k3dgrd,Aval,Vinfo,ncname)
!
!  Define turbulent kinetic energy time length scale.
!
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
        status=def_var(ncRSTid(ng),rstVid(idMtls,ng),NF_FRST,           &
     &                 nvd5,k3dgrd,Aval,Vinfo,ncname)
!
!  Define vertical mixing turbulent length scale.
!
        Vinfo( 1)=Vname(1,idVmLS)
        Vinfo( 2)=Vname(2,idVmLS)
        Vinfo( 3)=Vname(3,idVmLS)
        Vinfo(14)=Vname(4,idVmLS)
        Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_rho'
#  endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idVmLS,ng),r8)
        status=def_var(ncRSTid(ng),rstVid(idVmLS,ng),NF_FRST,           &
     &                 nvd4,w3dgrd,Aval,Vinfo,ncname)
!
!  Define turbulent kinetic energy vertical diffusion coefficient.
!
        Vinfo( 1)=Vname(1,idVmKK)
        Vinfo( 2)=Vname(2,idVmKK)
        Vinfo( 3)=Vname(3,idVmKK)
        Vinfo(14)=Vname(4,idVmKK)
        Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_rho'
#  endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idVmKK,ng),r8)
        status=def_var(ncRSTid(ng),rstVid(idVmKK,ng),NF_FRST,           &
     &                 nvd4,w3dgrd,Aval,Vinfo,ncname)
#  ifdef GLS_MIXING
!
!  Define turbulent length scale vertical diffusion coefficient.
!
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
        status=def_var(ncRSTid(ng),rstVid(idVmKP,ng),NF_FRST,           &
     &                 nvd4,w3dgrd,Aval,Vinfo,ncname)
#  endif
# endif
# ifdef SEDIMENT
#  ifdef BEDLOAD
!
!  Define Bedload U-direction.
!
        DO i=1,NST
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
          status=def_var(ncRSTid(ng),rstVid(idUbld(i),ng),NF_FRST,      &
     &                     nvd3,su2dgrd,Aval,Vinfo,ncname)
!
!  Define Bedload V-direction.
!
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
          status=def_var(ncRSTid(ng),rstVid(idVbld(i),ng),NF_FRST,      &
     &                   nvd3,sv2dgrd,Aval,Vinfo,ncname)
        END DO
#  endif
!
!  Define sediment fraction of each size class in each bed layer.
!
        DO i=1,NST
          Vinfo( 1)=Vname(1,idfrac(i))
          Vinfo( 2)=Vname(2,idfrac(i))
          Vinfo( 3)=Vname(3,idfrac(i))
          Vinfo(14)=Vname(4,idfrac(i))
          Vinfo(16)=Vname(1,idtime)
          WRITE (Vinfo(19),20) 1000.0_r8*Sd50(i,ng)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
# endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idfrac(i),ng),r8)
          status=def_var(ncRSTid(ng),rstVid(idfrac(i),ng),NF_FRST,      &
     &                   nvd4,sr3dgrd,Aval,Vinfo,ncname)
        END DO
!
!  Define sediment mass of each size class in each bed layer.
!
        DO i=1,NST
          Vinfo( 1)=Vname(1,idBmas(i))
          Vinfo( 2)=Vname(2,idBmas(i))
          Vinfo( 3)=Vname(3,idBmas(i))
          Vinfo(14)=Vname(4,idBmas(i))
          Vinfo(16)=Vname(1,idtime)
          WRITE (Vinfo(19),20) 1000.0_r8*Sd50(i,ng)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
# endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idBmas(i),ng),r8)
          status=def_var(ncRSTid(ng),rstVid(idBmas(i),ng),NF_FRST,      &
     &                   nvd4,sr3dgrd,Aval,Vinfo,ncname)
        END DO
!
!  Define sediment properties in each bed layer.
!
        DO i=1,MBEDP
          Vinfo( 1)=Vname(1,idSbed(i))
          Vinfo( 2)=Vname(2,idSbed(i))
          Vinfo( 3)=Vname(3,idSbed(i))
          Vinfo(14)=Vname(4,idSbed(i))
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
# endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idSbed(i),ng),r8)
          status=def_var(ncRSTid(ng),rstVid(idSbed(i),ng),NF_FRST,      &
     &                   nvd4,sr3dgrd,Aval,Vinfo,ncname)
        END DO
# endif
# if defined SEDIMENT || defined BBL_MODEL
!
!  define exposed sediment layer properties. Notice that only the
!  first four properties (mean grain diameter, mean grain density,
!  mean settling velocity, mean critical erosion stress, 
!  ripple length and ripple height) are written.
!
        DO i=1,6
          Vinfo( 1)=Vname(1,idBott(i))
          Vinfo( 2)=Vname(2,idBott(i))
          Vinfo( 3)=Vname(3,idBott(i))
          Vinfo(14)=Vname(4,idBott(i))
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
# endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idBott(i),ng),r8)
          status=def_var(ncRSTid(ng),rstVid(idBott(i),ng),NF_FRST,      &
     &                   nvd3,sr2dgrd,Aval,Vinfo,ncname)
        END DO
# endif
#endif
#ifdef NEARSHORE_MELLOR
!
!  Define 2D u-stokes velocity.
!
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
          status=def_var(ncRSTid(ng),rstVid(idU2Sd,ng),NF_FRST,         &
     &                   nvd3,su2dgrd,Aval,Vinfo,ncname)
!
!  Define 2D v-stokes velocity.
!
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
          status=def_var(ncRSTid(ng),rstVid(idV2Sd,ng),NF_FRST,         &
     &                   nvd3,sv2dgrd,Aval,Vinfo,ncname)
# ifdef SOLVE3D
!
!  Define 3D u-stokes velocity.
!
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
          status=def_var(ncRSTid(ng),rstVid(idU3Sd,ng),NF_FRST,         &
     &                   nvd4,su3dgrd,Aval,Vinfo,ncname)
!
!  Define 3D v-stokes velocity.
!
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
          status=def_var(ncRSTid(ng),rstVid(idV3Sd,ng),NF_FRST,         &
     &                   nvd4,sv3dgrd,Aval,Vinfo,ncname)
# endif
#endif
#ifdef WET_DRY
!
!  Define wet_dry masking.
!
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
          status=def_var(ncRSTid(ng),rstVid(idWmsk,ng),NF_FRST,         &
     &                   nvd3,sr2dgrd,Aval,Vinfo,ncname)
#endif
!
!-----------------------------------------------------------------------
!  Leave definition mode.
!-----------------------------------------------------------------------
!
        status=nf90_enddef(ncRSTid(ng))
      END IF
!
!-----------------------------------------------------------------------
!  Write out time-recordless, information variables.
!-----------------------------------------------------------------------
!
      IF (Ldefine) THEN
        status=wrt_info(ng, iNLM, ncRSTid(ng), OutThread, ncname)
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!=======================================================================
!  Open an existing restart file, check its contents, and prepare for
!  appending data.
!=======================================================================
!
      IF (.not.Ldefine.and.OutThread) THEN
!
!  Inquire about the contents of restart NetCDF file:  Inquire about
!  the dimensions and variables.  Check for consistency.
!
        ncname=RSTname(ng)
        CALL opencdf (ng, 1, ncname, fname, N(ng), 0, nrec, nvd, Vsize)
        IF (exit_flag.ne.NoError) RETURN
!
!  Open restart file for read/write.
!
        status=nf90_open(TRIM(ncname),nf90_write,ncRSTid(ng))
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,40) TRIM(ncname)
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
!  restart variables. Get variable IDs.
!
        DO i=1,nvars
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idtime))) THEN
            got_var(idtime)=.TRUE.
            status=nf90_inq_varid(ncRSTid(ng),TRIM(Vname(1,idtime)),    &
     &                            rstVid(idtime,ng))
          END IF
#if defined SEDIMENT && defined SED_MORPH
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idbath))) THEN
            got_var(idbath)=.TRUE.
            status=nf90_inq_varid(ncRSTid(ng),TRIM(Vname(1,idbath)),    &
     &                            rstVid(idbath,ng))
          END IF
#endif
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idFsur))) THEN
            got_var(idFsur)=.TRUE.
            status=nf90_inq_varid(ncRSTid(ng),TRIM(Vname(1,idFsur)),    &
     &                            rstVid(idFsur,ng))
          END IF
#ifdef PERFECT_RESTART
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idRzet))) THEN
            got_var(idRzet)=.TRUE.
            status=nf90_inq_varid(ncRSTid(ng),TRIM(Vname(1,idRzet)),    &
     &                            rstVid(idRzet,ng))
          END IF
#endif
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idUbar))) THEN
            got_var(idUbar)=.TRUE.
            status=nf90_inq_varid(ncRSTid(ng),TRIM(Vname(1,idUbar)),    &
     &                            rstVid(idUbar,ng))
          END IF
#ifdef PERFECT_RESTART
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idRu2d))) THEN
            got_var(idRu2d)=.TRUE.
            status=nf90_inq_varid(ncRSTid(ng),TRIM(Vname(1,idRu2d)),    &
     &                            rstVid(idRu2d,ng))
          END IF
#endif
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idVbar))) THEN
            got_var(idVbar)=.TRUE.
            status=nf90_inq_varid(ncRSTid(ng),TRIM(Vname(1,idVbar)),    &
     &                            rstVid(idVbar,ng))
          END IF
#ifdef PERFECT_RESTART
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idRv2d))) THEN
            got_var(idRv2d)=.TRUE.
            status=nf90_inq_varid(ncRSTid(ng),TRIM(Vname(1,idRv2d)),    &
     &                            rstVid(idRv2d,ng))
          END IF
#endif
#ifdef SOLVE3D
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idUvel))) THEN
            got_var(idUvel)=.TRUE.
            status=nf90_inq_varid(ncRSTid(ng),TRIM(Vname(1,idUvel)),    &
     &                            rstVid(idUvel,ng))
          END IF
# ifdef PERFECT_RESTART
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idRu3d))) THEN
            got_var(idRu3d)=.TRUE.
            status=nf90_inq_varid(ncRSTid(ng),TRIM(Vname(1,idRu3d)),    &
     &                            rstVid(idRu3d,ng))
          END IF
# endif
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idVvel))) THEN
            got_var(idVvel)=.TRUE.
            status=nf90_inq_varid(ncRSTid(ng),TRIM(Vname(1,idVvel)),    &
     &                            rstVid(idVvel,ng))
          END IF
# ifdef PERFECT_RESTART
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idRv3d))) THEN
            got_var(idRv3d)=.TRUE.
            status=nf90_inq_varid(ncRSTid(ng),TRIM(Vname(1,idRv3d)),    &
     &                            rstVid(idRv3d,ng))
          END IF
# endif
          DO itrc=1,NT(ng)
            IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idTvar(itrc)))) THEN
             got_var(idTvar(itrc))=.TRUE.
             status=nf90_inq_varid(ncRSTid(ng),                         &
     &                             TRIM(Vname(1,idTvar(itrc))),         &
     &                             rstTid(itrc,ng))
            END IF
          END DO
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idDano))) THEN
            got_var(idDano)=.TRUE.
            status=nf90_inq_varid(ncRSTid(ng),TRIM(Vname(1,idDano)),    &
     &                            rstVid(idDano,ng))
          END IF
# ifdef LMD_SKPP
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idHsbl))) THEN
            got_var(idHsbl)=.TRUE.
            status=nf90_inq_varid(ncRSTid(ng),TRIM(Vname(1,idHsbl)),    &
     &                            rstVid(idHsbl,ng))
          END IF
# endif
# ifdef LMD_BKPP
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idHbbl))) THEN
            got_var(idHbbl)=.TRUE.
            status=nf90_inq_varid(ncRSTid(ng),TRIM(Vname(1,idHbbl)),    &
     &                            rstVid(idHbbl,ng))
          END IF
# endif
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idVvis))) THEN
            got_var(idVvis)=.TRUE.
            status=nf90_inq_varid(ncRSTid(ng),TRIM(Vname(1,idVvis)),    &
     &                            rstVid(idVvis,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idTdif))) THEN
            got_var(idTdif)=.TRUE.
            status=nf90_inq_varid(ncRSTid(ng),TRIM(Vname(1,idTdif)),    &
     &                            rstVid(idTdif,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idSdif))) THEN
            got_var(idSdif)=.TRUE.
            status=nf90_inq_varid(ncRSTid(ng),TRIM(Vname(1,idSdif)),    &
     &                            rstVid(idSdif,ng))
          END IF
# if defined PERFECT_RESTART && \
     (defined GLS_MIXING     || defined MY25_MIXING)
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idMtke))) THEN
            got_var(idMtke)=.TRUE.
            status=nf90_inq_varid(ncRSTid(ng),TRIM(Vname(1,idMtke)),    &
     &                            rstVid(idMtke,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idMtls))) THEN
            got_var(idMtls)=.TRUE.
            status=nf90_inq_varid(ncRSTid(ng),TRIM(Vname(1,idMtls)),    &
     &                            rstVid(idMtls,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idVmLS))) THEN
            got_var(idVmLS)=.TRUE.
            status=nf90_inq_varid(ncRSTid(ng),TRIM(Vname(1,idVmLS)),    &
     &                            rstVid(idVmLS,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idVmKK))) THEN
            got_var(idVmKK)=.TRUE.
            status=nf90_inq_varid(ncRSTid(ng),TRIM(Vname(1,idVmKK)),    &
     &                            rstVid(idVmKK,ng))
          END IF
#  ifdef GLS_MIXING
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idVmKP))) THEN
            got_var(idVmKP)=.TRUE.
            status=nf90_inq_varid(ncRSTid(ng),TRIM(Vname(1,idVmKP)),    &
     &                            rstVid(idVmKP,ng))
          END IF
#  endif
# endif
# ifdef SEDIMENT
          DO itrc=1,NST
#  ifdef BEDLOAD
            IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idUbld(itrc)))) THEN
              got_var(idUbld(itrc))=.true.
              status=nf90_inq_varid(ncRSTid(ng),                        &
     &                              TRIM(Vname(1,idUbld(itrc))),        &
     &                              rstVid(idUbld(itrc),ng))
            END IF
            IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idVbld(itrc)))) THEN
              got_var(idVbld(itrc))=.true.
              status=nf90_inq_varid(ncRSTid(ng),                        &
     &                              TRIM(Vname(1,idVbld(itrc))),        &
     &                              rstVid(idVbld(itrc),ng))
            END IF
#  endif
            IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idfrac(itrc)))) THEN
              got_var(idfrac(itrc))=.TRUE.
              status=nf90_inq_varid(ncRSTid(ng),                        &
     &                              TRIM(Vname(1,idfrac(itrc))),        &
     &                              rstVid(idfrac(itrc),ng))
            END IF
            IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idBmas(itrc)))) THEN
              got_var(idBmas(itrc))=.TRUE.
              status=nf90_inq_varid(ncRSTid(ng),                        &
     &                              TRIM(Vname(1,idBmas(itrc))),        &
     &                              rstVid(idBmas(itrc),ng))
            END IF
          END DO
          DO itrc=1,MBEDP
            IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idSbed(itrc)))) THEN
              got_var(idSbed(itrc))=.TRUE.
              status=nf90_inq_varid(ncRSTid(ng),                        &
     &                              TRIM(Vname(1,idSbed(itrc))),        &
     &                              rstVid(idSbed(itrc),ng))
            END IF
          END DO
# endif
# if defined SEDIMENT || defined BBL_MODEL
          DO itrc=1,6
            IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idBott(itrc)))) THEN
              got_var(idBott(itrc))=.TRUE.
              status=nf90_inq_varid(ncRSTid(ng),                        &
     &                              TRIM(Vname(1,idBott(itrc))),        &
     &                              rstVid(idBott(itrc),ng))
            END IF
          END DO
# endif
#endif
#if defined WET_DRY
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idWmsk))) THEN
            got_var(idWmsk)=.TRUE.
            status=nf90_inq_varid(ncRSTid(ng),                          &
     &                            TRIM(Vname(1,idWmsk)),                &
     &                            rstVid(idWmsk,ng))
            END IF
#endif
#if defined NEARSHORE_MELLOR
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idU2Sd))) THEN
            got_var(idU2Sd)=.TRUE.
            status=nf90_inq_varid(ncRSTid(ng),                          &
     &                            TRIM(Vname(1,idU2Sd)),                &
     &                            rstVid(idU2Sd,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idV2Sd))) THEN
            got_var(idV2Sd)=.TRUE.
            status=nf90_inq_varid(ncRSTid(ng),                          &
     &                            TRIM(Vname(1,idV2Sd)),                &
     &                            rstVid(idV2Sd,ng))
          END IF
# ifdef SOLVE3D
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idU3Sd))) THEN
            got_var(idU3Sd)=.TRUE.
            status=nf90_inq_varid(ncRSTid(ng),                          &
     &                            TRIM(Vname(1,idU3Sd)),                &
     &                            rstVid(idU3Sd,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idV3Sd))) THEN
            got_var(idV3Sd)=.TRUE.
            status=nf90_inq_varid(ncRSTid(ng),                          &
     &                            TRIM(Vname(1,idV3Sd)),                &
     &                            rstVid(idV3Sd,ng))
          END IF
# endif
#endif
        END DO
!
!  Check if initialization variables are available in input NetCDF
!  file.
!
        IF (.not.got_var(idtime)) THEN
          WRITE (stdout,40) TRIM(Vname(1,idtime)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#if defined SEDIMENT && defined SED_MORPH
        IF (.not.got_var(idbath)) THEN
          WRITE (stdout,40) TRIM(Vname(1,idbath)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#endif
        IF (.not.got_var(idFsur)) THEN
          WRITE (stdout,40) TRIM(Vname(1,idFsur)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#ifdef PERFECT_RESTART
        IF (.not.got_var(idRzet)) THEN
          WRITE (stdout,40) TRIM(Vname(1,idRzet)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#endif
        IF (.not.got_var(idUbar)) THEN
          WRITE (stdout,40) TRIM(Vname(1,idUbar)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#ifdef PERFECT_RESTART
        IF (.not.got_var(idRu2d)) THEN
          WRITE (stdout,40) TRIM(Vname(1,idRu2d)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#endif
        IF (.not.got_var(idVbar)) THEN
          WRITE (stdout,40) TRIM(Vname(1,idVbar)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#ifdef PERFECT_RESTART
        IF (.not.got_var(idRv2d)) THEN
          WRITE (stdout,40) TRIM(Vname(1,idRv2d)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#endif
#ifdef SOLVE3D
        IF (.not.got_var(idUvel)) THEN
          WRITE (stdout,40) TRIM(Vname(1,idUvel)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# ifdef PERFECT_RESTART
        IF (.not.got_var(idRu3d)) THEN
          WRITE (stdout,40) TRIM(Vname(1,idRu3d)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
        IF (.not.got_var(idVvel)) THEN
          WRITE (stdout,40) TRIM(Vname(1,idVvel)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# ifdef PERFECT_RESTART
        IF (.not.got_var(idRv3d)) THEN
          WRITE (stdout,40) TRIM(Vname(1,idRv3d)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
        DO itrc=1,NT(ng)
          IF (.not.got_var(idTvar(itrc))) THEN
            WRITE (stdout,40) TRIM(Vname(1,idTvar(itrc))), TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
        IF (.not.got_var(idDano)) THEN
          WRITE (stdout,40) TRIM(Vname(1,idDano)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# if defined PERFECT_RESTART && \
     (defined GLS_MIXING     || defined MY25_MIXING)
        IF (.not.got_var(idMtke)) THEN
          WRITE (stdout,40) TRIM(Vname(1,idMtke)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idMtls)) THEN
          WRITE (stdout,40) TRIM(Vname(1,idMtls)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVmLS)) THEN
          WRITE (stdout,40) TRIM(Vname(1,idVmLS)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVmKK)) THEN
          WRITE (stdout,40) TRIM(Vname(1,idVmKK)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  ifdef GLS_MIXING
        IF (.not.got_var(idVmKP)) THEN
          WRITE (stdout,40) TRIM(Vname(1,idVmKP)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
# endif
# ifdef SEDIMENT
        DO i=1,NST
#  ifdef BEDLOAD
          IF (.not.got_var(idUbld(i))) THEN
            WRITE (stdout,40) TRIM(Vname(1,idUbld(i))), TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
          IF (.not.got_var(idVbld(i))) THEN
            WRITE (stdout,40) TRIM(Vname(1,idVbld(i))), TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
#  endif
          IF (.not.got_var(idfrac(i))) THEN
            WRITE (stdout,40) TRIM(Vname(1,idfrac(i))), TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
          IF (.not.got_var(idBmas(i))) THEN
            WRITE (stdout,40) TRIM(Vname(1,idBmas(i))), TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
        DO i=1,MBEDP
          IF (.not.got_var(idSbed(i))) THEN
            WRITE (stdout,40) TRIM(Vname(1,idSbed(i))), TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
# endif
# if defined SEDIMENT || defined BBL_MODEL
        DO i=1,6
          IF (.not.got_var(idBott(i))) THEN
            WRITE (stdout,40) TRIM(Vname(1,idBott(i))), TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
# endif
#endif
#if defined NEARSHORE_MELLOR
        IF (.not.got_var(idU2Sd)) THEN
          WRITE (stdout,40) TRIM(Vname(1,idU2Sd)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idV2Sd)) THEN
          WRITE (stdout,40) TRIM(Vname(1,idV2Sd)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# ifdef SOLVE3D
        IF (.not.got_var(idU3Sd)) THEN
          WRITE (stdout,40) TRIM(Vname(1,idU3Sd)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idV3Sd)) THEN
          WRITE (stdout,40) TRIM(Vname(1,idV3Sd)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
#endif
#if defined WET_DRY
        IF (.not.got_var(idWmsk)) THEN
          WRITE (stdout,40) TRIM(Vname(1,idWmsk)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#endif
!
!  Set unlimited time record dimension to current value.
!
        IF (LcycleRST(ng)) THEN
          tRSTindx(ng)=0
        ELSE
          tRSTindx(ng)=tsize
        END IF
      END IF
!
  10  FORMAT (/,' DEF_RST - unable to create restart NetCDF file: ',a)
  20  FORMAT (1pe11.4,1x,'millimeter')
  30  FORMAT (/,' DEF_RST - unable to open restart NetCDF file: ',a)
  40  FORMAT (/,' DEF_RST - unable to find variable: ',a,2x,            &
     &        ' in restart NetCDF file: ',a)
      RETURN
      END SUBROUTINE def_rst
