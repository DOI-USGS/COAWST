#include "cppdefs.h"
#if defined FOUR_DVAR || defined VERIFICATION
      SUBROUTINE read_AssPar (model, inp, out, Lwrite)
!
!svn $Id: read_asspar.F 807 2016-07-09 02:03:55Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine reads and reports input data assimilation parameters.  !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
# if defined FOUR_DVAR  || defined VERIFICATION || \
    (defined HESSIAN_SV && defined BNORM)
      USE mod_fourdvar
# endif
      USE mod_iounits
      USE mod_ncparam
      USE mod_scalars
!
      implicit none
!
!  Imported variable declarations
!
      logical, intent(in) :: Lwrite
      integer, intent(in) :: model, inp, out
!
!  Local variable declarations.
!
      logical :: find_file
      logical :: Lvalue(1)

      integer :: Mval, Npts, Nval
      integer :: i, ib, igrid, itrc, k, ng, status
      integer :: Cdim, Clen, Rdim
      integer :: Ivalue(1)

      integer :: decode_line, load_i, load_l, load_r, load_s1d

      integer :: Nfiles(Ngrids)
# if defined FOUR_DVAR  || defined VERIFICATION || \
    (defined HESSIAN_SV && defined BNORM)
      logical, dimension(MT) :: Ltracer
#  if defined ADJUST_STFLUX && defined SOLVE3D
      logical, dimension(MT,Ngrids) :: Ltsur
#  endif
#  ifdef ADJUST_BOUNDARY
      logical, dimension(4,Ngrids) :: Lbry
#   ifdef SOLVE3D
      logical, dimension(4,MT,Ngrids) :: Lbry_trc
      logical, dimension(MT,4) :: Lboundary
      real(r8), dimension(MT,4,Ngrids) :: Rboundary
#   endif
#  endif
      real(r8), dimension(MT,Ngrids) :: Rtracer
# endif
      real(r8), dimension(100) :: Rval
      real(r8) :: Rvalue(1)

      character (len=1 ), parameter :: blank = ' '
# ifdef ADJUST_BOUNDARY
      character (len=7) :: Text
# endif
      character (len=40 ) :: KeyWord
      character (len=50 ) :: label
      character (len=256) :: fname, line
      character (len=256), dimension(200) :: Cval
!
!-----------------------------------------------------------------------
!  Initialize.
!-----------------------------------------------------------------------
!
      igrid=1
      Nfiles(1:Ngrids)=0
      DO i=1,LEN(label)
        label(i:i)=blank
      END DO
      Cdim=SIZE(Cval,1)
      Clen=LEN(Cval(1))
      Rdim=SIZE(Rval,1)

# if defined FOUR_DVAR  || defined VERIFICATION || \
    (defined HESSIAN_SV && defined BNORM)
!
!-----------------------------------------------------------------------
!  Read in 4DVAR assimilation parameters. Then, load input data into
!  module.
!-----------------------------------------------------------------------
!
      DO WHILE (.TRUE.)
        READ (inp,'(a)',ERR=10,END=20) line
        status=decode_line(line, KeyWord, Nval, Cval, Rval)
        IF (status.gt.0) THEN
          SELECT CASE (TRIM(KeyWord))
            CASE ('dTdz_min')
              Npts=load_r(Nval, Rval, Ngrids, dTdz_min)
            CASE ('ml_depth')
              Npts=load_r(Nval, Rval, Ngrids, ml_depth)
              DO ng=1,Ngrids
                ml_depth(ng)=ABS(ml_depth(ng))
              END DO
#  if defined BALANCE_OPERATOR && defined ZETA_ELLIPTIC
            CASE ('Nbico')
              Npts=load_i(Nval, Rval, Ngrids, Nbico)
#  endif
            CASE ('LNM_depth')
              Npts=load_r(Nval, Rval, Ngrids, LNM_depth)
              DO ng=1,Ngrids
                LNM_depth(ng)=ABS(LNM_depth(ng))
              END DO
            CASE ('LNM_flag')
              Npts=load_i(Nval, Rval, 1, Ivalue)
              LNM_flag=Ivalue(1)
#  if defined WEAK_CONSTRAINT && \
     (defined ARRAY_MODES     || defined CLIPPING)
            CASE ('Nvct')
              Npts=load_i(Nval, Rval, 1, Ivalue)
              Nvct=Ivalue(1)
#  endif
            CASE ('GradErr')
              Npts=load_r(Nval, Rval, 1, Rvalue)
              GradErr=Rvalue(1)
            CASE ('HevecErr')
              Npts=load_r(Nval, Rval, 1, Rvalue)
              HevecErr=Rvalue(1)
            CASE ('LhessianEV')
              Npts=load_l(Nval, Cval, 1, Lvalue)
              LhessianEV=Lvalue(1)
            CASE ('LhotStart')
              Npts=load_l(Nval, Cval, 1, Lvalue)
              LhotStart=Lvalue(1)
            CASE ('Lprecond')
              Npts=load_l(Nval, Cval, 1, Lvalue)
              Lprecond=Lvalue(1)
#  if defined WEAK_CONSTRAINT
              IF ( LhessianEV.and.Lprecond ) THEN
                LhessianEV=.FALSE.
              END IF
#  endif
            CASE ('Lritz')
              Npts=load_l(Nval, Cval, 1, Lvalue)
              Lritz=Lvalue(1)
            CASE ('NritzEV')
              Npts=load_i(Nval, Rval, 1, Ivalue)
              NritzEV=Ivalue(1)
            CASE ('NpostI')
              Npts=load_i(Nval, Rval, 1, Ivalue)
              NpostI=Ivalue(1)
            CASE ('tl_M2diff')
              Npts=load_r(Nval, Rval, Ngrids, tl_M2diff)
            CASE ('tl_M3diff')
              Npts=load_r(Nval, Rval, Ngrids, tl_M3diff)
            CASE ('tl_Tdiff')
              Npts=load_r(Nval, Rval, MT*Ngrids, Rtracer)
              DO ng=1,Ngrids
                DO itrc=1,NT(ng)
                  tl_Tdiff(itrc,ng)=Rtracer(itrc,ng)
                END DO
              END DO
            CASE ('LdefNRM')
              Npts=load_l(Nval, Cval, 4*Ngrids, LdefNRM)
            CASE ('LwrtNRM')
              Npts=load_l(Nval, Cval, 4*Ngrids, LwrtNRM)
            CASE ('CnormM(isFsur)')
              IF (isFsur.eq.0) THEN
                IF (Master) WRITE (out,210) 'isFsur'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, 1, Cnorm(2,isFsur))
            CASE ('CnormM(isUbar)')
              IF (isUbar.eq.0) THEN
                IF (Master) WRITE (out,210) 'isUbar'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, 1, Cnorm(2,isUbar))
            CASE ('CnormM(isVbar)')
              IF (isVbar.eq.0) THEN
                IF (Master) WRITE (out,210) 'isVbar'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, 1, Cnorm(2,isVbar))
#  ifdef SOLVE3D
            CASE ('CnormM(isUvel)')
              IF (isUvel.eq.0) THEN
                IF (Master) WRITE (out,210) 'isUvel'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, 1, Cnorm(2,isUvel))
            CASE ('CnormM(isVvel)')
              IF (isVvel.eq.0) THEN
                IF (Master) WRITE (out,210) 'isVvel'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, 1, Cnorm(2,isVvel))
            CASE ('CnormM(isTvar)')
              IF (MAXVAL(isTvar).eq.0) THEN
                IF (Master) WRITE (out,210) 'isTvar'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, MT, Ltracer)
              DO itrc=1,MT
                i=isTvar(itrc)
                Cnorm(2,i)=Ltracer(itrc)
              END DO
#  endif
            CASE ('CnormI(isFsur)')
              Npts=load_l(Nval, Cval, 1, Cnorm(1,isFsur))
            CASE ('CnormI(isUbar)')
              Npts=load_l(Nval, Cval, 1, Cnorm(1,isUbar))
            CASE ('CnormI(isVbar)')
              Npts=load_l(Nval, Cval, 1, Cnorm(1,isVbar))
#  ifdef SOLVE3D
            CASE ('CnormI(isUvel)')
              Npts=load_l(Nval, Cval, 1, Cnorm(1,isUvel))
            CASE ('CnormI(isVvel)')
              Npts=load_l(Nval, Cval, 1, Cnorm(1,isVvel))
            CASE ('CnormI(isTvar)')
              Npts=load_l(Nval, Cval, MT, Ltracer)
              DO itrc=1,MT
                i=isTvar(itrc)
                Cnorm(1,i)=Ltracer(itrc)
              END DO
#  endif
#  ifdef ADJUST_BOUNDARY
            CASE ('CnormB(isFsur)')
              Npts=load_l(Nval, Cval, 4, CnormB(isFsur,:))
            CASE ('CnormB(isUbar)')
              Npts=load_l(Nval, Cval, 4, CnormB(isUbar,:))
            CASE ('CnormB(isVbar)')
              Npts=load_l(Nval, Cval, 4, CnormB(isVbar,:))
#   ifdef SOLVE3D
            CASE ('CnormB(isUvel)')
              Npts=load_l(Nval, Cval, 4, CnormB(isUvel,:))
            CASE ('CnormB(isVvel)')
              Npts=load_l(Nval, Cval, 4, CnormB(isVvel,:))
            CASE ('CnormB(isTvar)')
              Npts=load_l(Nval, Cval, MT*4, Lboundary)
              DO ib=1,4
                DO itrc=1,MT
                  i=isTvar(itrc)
                  CnormB(i,ib)=Lboundary(itrc,ib)
                END DO
              END DO
#   endif
#  endif
#  ifdef ADJUST_WSTRESS
            CASE ('CnormF(isUstr)')
              IF (isUstr.eq.0) THEN
                IF (Master) WRITE (out,210) 'isUstr'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, 1, Cnorm(1,isUstr))
            CASE ('CnormF(isVstr)')
              IF (isVstr.eq.0) THEN
                IF (Master) WRITE (out,210) 'isVstr'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, 1, Cnorm(1,isVstr))
#  endif
#  if defined ADJUST_STFLUX && defined SOLVE3D
            CASE ('CnormF(isTsur)')
              IF (MAXVAL(isTsur).eq.0) THEN
                IF (Master) WRITE (out,210) 'isTsur'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, MT, Ltracer)
              DO itrc=1,MT
                i=isTsur(itrc)
                Cnorm(1,i)=Ltracer(itrc)
              END DO
#  endif
            CASE ('balance(isSalt)')
              Npts=load_l(Nval, Cval, 1, balance(isTvar(isalt)))
            CASE ('balance(isFsur)')
              Npts=load_l(Nval, Cval, 1, balance(isFsur))
            CASE ('balance(isVbar)')
              Npts=load_l(Nval, Cval, 1, balance(isVbar))
            CASE ('balance(isVvel)')
              Npts=load_l(Nval, Cval, 1, balance(isVvel))
            CASE ('Nmethod')
              Npts=load_i(Nval, Rval, Ngrids, Nmethod)
            CASE ('Rscheme')
              Npts=load_i(Nval, Rval, Ngrids, Rscheme)
            CASE ('Nrandom')
              Npts=load_i(Nval, Rval, 1, Ivalue)
              Nrandom=Ivalue(1)
            CASE ('Hgamma')
              Npts=load_r(Nval, Rval, 4, Hgamma)
#  ifdef SOLVE3D
            CASE ('Vgamma')
              Npts=load_r(Nval, Rval, 4, Vgamma)
#  endif
            CASE ('HdecayM(isFsur)')
              Npts=load_r(Nval, Rval, Ngrids, Hdecay(2,isFsur,:))
            CASE ('HdecayM(isUbar)')
              Npts=load_r(Nval, Rval, Ngrids, Hdecay(2,isUbar,:))
            CASE ('HdecayM(isVbar)')
              Npts=load_r(Nval, Rval, Ngrids, Hdecay(2,isVbar,:))
#  ifdef SOLVE3D
            CASE ('HdecayM(isUvel)')
              Npts=load_r(Nval, Rval, Ngrids, Hdecay(2,isUvel,:))
            CASE ('HdecayM(isVvel)')
              Npts=load_r(Nval, Rval, Ngrids, Hdecay(2,isVvel,:))
            CASE ('HdecayM(isTvar)')
              Npts=load_r(Nval, Rval, MT*Ngrids, Rtracer)
              DO ng=1,Ngrids
                DO itrc=1,NT(ng)
                  Hdecay(2,isTvar(itrc),ng)=Rtracer(itrc,ng)
                END DO
              END DO
            CASE ('VdecayM(isUvel)')
              Npts=load_r(Nval, Rval, Ngrids, Vdecay(2,isUvel,:))
            CASE ('VdecayM(isVvel)')
              Npts=load_r(Nval, Rval, Ngrids, Vdecay(2,isVvel,:))
            CASE ('VdecayM(isTvar)')
              Npts=load_r(Nval, Rval, MT*Ngrids, Rtracer)
              DO ng=1,Ngrids
                DO itrc=1,NT(ng)
                  Vdecay(2,isTvar(itrc),ng)=Rtracer(itrc,ng)
                END DO
              END DO
#  endif
            CASE ('TdecayM(isFsur)')
              Npts=load_r(Nval, Rval, Ngrids, Tdecay(isFsur,:))
            CASE ('TdecayM(isUbar)')
              Npts=load_r(Nval, Rval, Ngrids, Tdecay(isUbar,:))
            CASE ('TdecayM(isVbar)')
              Npts=load_r(Nval, Rval, Ngrids, Tdecay(isVbar,:))
#  ifdef SOLVE3D
            CASE ('TdecayM(isUvel)')
              Npts=load_r(Nval, Rval, Ngrids, Tdecay(isUvel,:))
            CASE ('TdecayM(isVvel)')
              Npts=load_r(Nval, Rval, Ngrids, Tdecay(isVvel,:))
            CASE ('TdecayM(isTvar)')
              Npts=load_r(Nval, Rval, MT*Ngrids, Rtracer)
              DO ng=1,Ngrids
                DO itrc=1,NT(ng)
                  Tdecay(isTvar(itrc),ng)=Rtracer(itrc,ng)
                END DO
              END DO
#  endif
            CASE ('HdecayI(isFsur)')
              Npts=load_r(Nval, Rval, Ngrids, Hdecay(1,isFsur,:))
            CASE ('HdecayI(isUbar)')
              Npts=load_r(Nval, Rval, Ngrids, Hdecay(1,isUbar,:))
            CASE ('HdecayI(isVbar)')
              Npts=load_r(Nval, Rval, Ngrids, Hdecay(1,isVbar,:))
#  ifdef SOLVE3D
            CASE ('HdecayI(isUvel)')
              Npts=load_r(Nval, Rval, Ngrids, Hdecay(1,isUvel,:))
            CASE ('HdecayI(isVvel)')
              Npts=load_r(Nval, Rval, Ngrids, Hdecay(1,isVvel,:))
            CASE ('HdecayI(isTvar)')
              Npts=load_r(Nval, Rval, MT*Ngrids, Rtracer)
              DO ng=1,Ngrids
                DO itrc=1,NT(ng)
                  Hdecay(1,isTvar(itrc),ng)=Rtracer(itrc,ng)
                END DO
              END DO
            CASE ('VdecayI(isUvel)')
              Npts=load_r(Nval, Rval, Ngrids, Vdecay(1,isUvel,:))
            CASE ('VdecayI(isVvel)')
              Npts=load_r(Nval, Rval, Ngrids, Vdecay(1,isVvel,:))
            CASE ('VdecayI(isTvar)')
              Npts=load_r(Nval, Rval, MT*Ngrids, Rtracer)
              DO ng=1,Ngrids
                DO itrc=1,NT(ng)
                  Vdecay(1,isTvar(itrc),ng)=Rtracer(itrc,ng)
                END DO
              END DO
#  endif
#  ifdef ADJUST_BOUNDARY
            CASE ('HdecayB(isFsur)')
              Npts=load_r(Nval, Rval, 4*Ngrids, HdecayB(isFsur,:,:))
            CASE ('HdecayB(isUbar)')
              Npts=load_r(Nval, Rval, 4*Ngrids, HdecayB(isUbar,:,:))
            CASE ('HdecayB(isVbar)')
              Npts=load_r(Nval, Rval, 4*Ngrids, HdecayB(isVbar,:,:))
#   ifdef SOLVE3D
            CASE ('HdecayB(isUvel)')
              Npts=load_r(Nval, Rval, 4*Ngrids, HdecayB(isUvel,:,:))
            CASE ('HdecayB(isVvel)')
              Npts=load_r(Nval, Rval, 4*Ngrids, HdecayB(isVvel,:,:))
            CASE ('HdecayB(isTvar)')
              Npts=load_r(Nval, Rval, 4*MT*Ngrids, Rboundary)
              DO ng=1,Ngrids
                DO ib=1,4
                  DO itrc=1,NT(ng)
                    HdecayB(isTvar(itrc),ib,ng)=Rboundary(itrc,ib,ng)
                  END DO
                END DO
              END DO
            CASE ('VdecayB(isUvel)')
              Npts=load_r(Nval, Rval, 4*Ngrids, VdecayB(isUvel,:,:))
            CASE ('VdecayB(isVvel)')
              Npts=load_r(Nval, Rval, 4*Ngrids, VdecayB(isVvel,:,:))
            CASE ('VdecayB(isTvar)')
              Npts=load_r(Nval, Rval, 4*MT*Ngrids, Rboundary)
              DO ng=1,Ngrids
                DO ib=1,4
                  DO itrc=1,NT(ng)
                    VdecayB(isTvar(itrc),ib,ng)=Rboundary(itrc,ib,ng)
                  END DO
                END DO
              END DO
#   endif
#  endif
#  ifdef ADJUST_WSTRESS
            CASE ('HdecayF(isUstr)')
              Npts=load_r(Nval, Rval, Ngrids, Hdecay(1,isUstr,:))
            CASE ('HdecayF(isVstr)')
              Npts=load_r(Nval, Rval, Ngrids, Hdecay(1,isVstr,:))
#  endif
#  if defined ADJUST_STFLUX && defined SOLVE3D
            CASE ('HdecayF(isTsur)')
              Npts=load_r(Nval, Rval, MT*Ngrids, Rtracer)
              DO ng=1,Ngrids
                DO itrc=1,NT(ng)
                  Hdecay(1,isTsur(itrc),ng)=Rtracer(itrc,ng)
                END DO
              END DO
#  endif
#  ifdef BGQC
            CASE ('bgqc_type')
              Npts=load_i(Nval, Rval, Ngrids, bgqc_type)
            CASE ('S_bgqc(isFsur)')
              Npts=load_r(Nval, Rval, Ngrids, S_bgqc(isFsur,:))
            CASE ('S_bgqc(isUbar)')
              Npts=load_r(Nval, Rval, Ngrids, S_bgqc(isUbar,:))
            CASE ('S_bgqc(isVbar)')
              Npts=load_r(Nval, Rval, Ngrids, S_bgqc(isVbar,:))
            CASE ('S_bgqc(isUvel)')
              Npts=load_r(Nval, Rval, Ngrids, S_bgqc(isUvel,:))
            CASE ('S_bgqc(isVvel)')
              Npts=load_r(Nval, Rval, Ngrids, S_bgqc(isVvel,:))
            CASE ('S_bgqc(isTvar)')
              Npts=load_r(Nval, Rval, MT*Ngrids, Rtracer)
              DO ng=1,Ngrids
                DO itrc=1,NT(ng)
                  S_bgqc(isTvar(itrc),ng)=Rtracer(itrc,ng)
                END DO
              END DO
            CASE ('Nprovenance')
              Npts=load_i(Nval, Rval, Ngrids, Nprovenance)
              IF (.not.allocated(Iprovenance)) THEN
                allocate ( Iprovenance(MAXVAL(Nprovenance),Ngrids))
              END IF
              IF (.not.allocated(P_bgqc)) THEN
                allocate ( P_bgqc(MAXVAL(Nprovenance),Ngrids))
              END IF
            CASE ('Iprovenance')
              Mval=MAXVAL(Nprovenance)*Ngrids
              Npts=load_i(Nval, Rval, Mval, Iprovenance)
            CASE ('P_bgqc')
              Mval=MAXVAL(Nprovenance)*Ngrids
              Npts=load_r(Nval, Rval, Mval, P_bgqc)
#  endif
#  if defined ADJUST_STFLUX && defined SOLVE3D
            CASE ('Lstflux')
              Npts=load_l(Nval, Cval, MT*Ngrids, Ltsur)
              DO ng=1,Ngrids
                DO itrc=1,NT(ng)
                  Lstflux(itrc,ng)=Ltsur(itrc,ng)
                END DO
              END DO
#  endif
#  ifdef ADJUST_BOUNDARY
            CASE ('Lobc(isFsur)')
              Npts=load_l(Nval, Cval, 4*Ngrids, Lbry)
              DO ng=1,Ngrids
                DO ib=1,4
                  Lobc(ib,isFsur,ng)=Lbry(ib,ng)
                END DO
              END DO
            CASE ('Lobc(isUbar)')
              Npts=load_l(Nval, Cval, 4*Ngrids, Lbry)
              DO ng=1,Ngrids
                DO ib=1,4
                  Lobc(ib,isUbar,ng)=Lbry(ib,ng)
                END DO
              END DO
            CASE ('Lobc(isVbar)')
              Npts=load_l(Nval, Cval, 4*Ngrids, Lbry)
              DO ng=1,Ngrids
                DO ib=1,4
                  Lobc(ib,isVbar,ng)=Lbry(ib,ng)
                END DO
              END DO
#   ifdef SOLVE3D
            CASE ('Lobc(isUvel)')
              Npts=load_l(Nval, Cval, 4*Ngrids, Lbry)
              DO ng=1,Ngrids
                DO ib=1,4
                  Lobc(ib,isUvel,ng)=Lbry(ib,ng)
                END DO
              END DO
            CASE ('Lobc(isVvel)')
              Npts=load_l(Nval, Cval, 4*Ngrids, Lbry)
              DO ng=1,Ngrids
                DO ib=1,4
                  Lobc(ib,isVvel,ng)=Lbry(ib,ng)
                END DO
              END DO
            CASE ('Lobc(isTvar)')
              Npts=load_l(Nval, Cval, 4*MT*Ngrids, Lbry_trc)
              DO ng=1,Ngrids
                DO itrc=1,NT(ng)
                  i=isTvar(itrc)
                  DO ib=1,4
                    Lobc(ib,i,ng)=Lbry_trc(ib,itrc,ng)
                  END DO
                END DO
              END DO
#   endif
#  endif
            CASE ('STDnameI')
              label='STD - initial conditions standard deviation'
              Npts=load_s1d(Nval, Cval, Cdim, line, label, igrid,       &
     &                      Nfiles, STD(1,:))
            CASE ('STDnameM')
              label='STD - model error standard deviation'
              Npts=load_s1d(Nval, Cval, Cdim, line, label, igrid,       &
     &                      Nfiles, STD(2,:))
            CASE ('STDnameB')
              label='STD - boundary conditions standard deviation'
              Npts=load_s1d(Nval, Cval, Cdim, line, label, igrid,       &
     &                      Nfiles, STD(3,:))
            CASE ('STDnameF')
              label='STD - surface forcing standard deviation'
              Npts=load_s1d(Nval, Cval, Cdim, line, label, igrid,       &
     &                      Nfiles, STD(4,:))
            CASE ('NRMnameI')
              label='NRM - initial conditions normalization'
              Npts=load_s1d(Nval, Cval, Cdim, line, label, igrid,       &
     &                      Nfiles, NRM(1,:))
            CASE ('NRMnameM')
              label='NRM - model error normalization'
              Npts=load_s1d(Nval, Cval, Cdim, line, label, igrid,       &
     &                      Nfiles, NRM(2,:))
            CASE ('NRMnameB')
              label='NRM - boundary conditions normalization'
              Npts=load_s1d(Nval, Cval, Cdim, line, label, igrid,       &
     &                      Nfiles, NRM(3,:))
            CASE ('NRMnameF')
              label='NRM - surface forcing normalization'
              Npts=load_s1d(Nval, Cval, Cdim, line, label, igrid,       &
     &                      Nfiles, NRM(4,:))
            CASE ('OBSname')
              label='OBS - data assimilation observations'
              Npts=load_s1d(Nval, Cval, Cdim, line, label, igrid,       &
     &                      Nfiles, OBS)
            CASE ('HSSname')
              label='HSS - Hessian eigenvectors'
              Npts=load_s1d(Nval, Cval, Cdim, line, label, igrid,       &
     &                      Nfiles, HSS)
            CASE ('LCZname')
              label='LCZ - Lanczos vectors'
              Npts=load_s1d(Nval, Cval, Cdim, line, label, igrid,       &
     &                      Nfiles, LCZ)
            CASE ('LZEname')
              label='LZE - Time-evolved Lanczos vectors'
              Npts=load_s1d(Nval, Cval, Cdim, line, label, igrid,       &
     &                      Nfiles, LZE)
            CASE ('MODname')
              label='DAV - 4D-Var data assimilation variables'
              Npts=load_s1d(Nval, Cval, Cdim, line, label, igrid,       &
     &                      Nfiles, DAV)
            CASE ('ERRname')
              label='ERR - 4D-Var posterior error covariance'
              Npts=load_s1d(Nval, Cval, Cdim, line, label, igrid,       &
     &                      Nfiles, ERR)
          END SELECT
        END IF
      END DO
 10   IF (Master) WRITE (out,50) line
      exit_flag=4
      RETURN
 20   CONTINUE

#  ifdef ADJUST_BOUNDARY
!
!-----------------------------------------------------------------------
!  Check switches to adjust boundaries for consistency.
!-----------------------------------------------------------------------
!
!  Make sure that both momentum components are activated for processing.
!  If adjusting 2D momentum in 3D applications, make sure that the
!  free-surface and 3D momentum switches are activated. This is because
!  the 2D momentum adjustments are computed from the vertical integral
!  of the 3D momentum increments.
!
      DO ng=1,Ngrids
        DO ib=1,4
          IF (.not.Lobc(ib,isUbar,ng).and.Lobc(ib,isVbar,ng)) THEN
            Lobc(ib,isUbar,ng)=.TRUE.
          END IF
          IF (.not.Lobc(ib,isVbar,ng).and.Lobc(ib,isUbar,ng)) THEN
            Lobc(ib,isVbar,ng)=.TRUE.
          END IF
#   ifdef SOLVE3D
          IF (.not.Lobc(ib,isUvel,ng).and.Lobc(ib,isVvel,ng)) THEN
            Lobc(ib,isUvel,ng)=.TRUE.
          END IF
          IF (.not.Lobc(ib,isVvel,ng).and.Lobc(ib,isUvel,ng)) THEN
            Lobc(ib,isVvel,ng)=.TRUE.
          END IF
          IF (.not.Lobc(ib,isFsur,ng).and.Lobc(ib,isUbar,ng)) THEN
            Lobc(ib,isFsur,ng)=.TRUE.
          END IF
          IF (.not.Lobc(ib,isUvel,ng).and.Lobc(ib,isUbar,ng)) THEN
            Lobc(ib,isUvel,ng)=.TRUE.
            Lobc(ib,isVvel,ng)=.TRUE.
          END IF
#   endif
        END DO
      END DO
#  endif
#  if defined WEAK_CONSTRAINT && \
     (defined ARRAY_MODES     || defined CLIPPING)
!
!-----------------------------------------------------------------------
!  Check array modes parameter
!-----------------------------------------------------------------------
!
!  Array modes parameter must be greater than zero and less or equal
!  to Ninner.
!
      IF ((Nvct.lt.1).or.(Nvct.gt.Ninner)) THEN
        IF (Lwrite) THEN
          WRITE (out,55) 'Illegal parameter for array modes, Nvct = ',  &
     &                   Nvct, Ninner, '1 =< Nvct =< Ninner.'
        END IF
        exit_flag=6
        RETURN
      END IF
#  endif
#  ifdef BALANCE_OPERATOR
!
!-----------------------------------------------------------------------
!  Check balance operator switches for consitency.
!-----------------------------------------------------------------------
!
#   ifdef SOLVE3D
      IF (.not.balance(isTvar(isalt)).and.balance(isFsur)) THEN
        balance(isTvar(isalt))=.TRUE.
      END IF
      IF (.not.balance(isTvar(isalt)).and.balance(isVvel)) THEN
        balance(isTvar(isalt))=.TRUE.
      END IF
      IF (balance(isTvar(isalt))) THEN
        balance(isTvar(itemp))=.TRUE.
      END IF
      IF (balance(isVvel)) THEN
        balance(isUvel)=.TRUE.
      END IF
      IF (balance(isVbar)) THEN
        balance(isVbar)=.FALSE.
      END IF
#   else
      IF (balance(isVbar)) THEN
        balance(isUbar)=.TRUE.
      END IF
#   endif
#  endif
!
!-----------------------------------------------------------------------
!  Report input parameters.
!-----------------------------------------------------------------------
!
      IF (Lwrite) THEN
        DO ng=1,Ngrids
#  if defined FOUR_DVAR
          WRITE (out,60) ng
#   ifdef WEAK_CONSTRAINT
#    ifdef RPM_RELAXATION
          WRITE (out,120) tl_M2diff, 'tl_M2diff',                       &
     &            'RPM 2D momentum diffusive relaxation coefficient.'
#     ifdef SOLVE3D
          WRITE (out,130) tl_M3diff, 'tl_M3diff',                       &
     &            'RPM 3D momentum diffusive relaxation coefficient.'
          DO itrc=1,NT(ng)
            WRITE (out,130) tl_Tdiff(itrc,ng), 'tl_Tdiff',              &
     &            'RPM tracer diffusive relaxation coefficient, ',      &
     &            TRIM(Vname(1,idTvar(itrc)))
          END DO
#     endif
#    endif
#   endif
#   ifndef IS4DVAR_SENSITIVITY
#    ifdef IS4DVAR
          WRITE (out,100) GradErr, 'GradErr',                           &
     &            'Upper bound on relative error of the gradient.'
          WRITE (out,100) HevecErr, 'HevecErr',                         &
     &            'Accuracy required for Hessian eigenvectors.'
          WRITE (out,70) LhessianEV, 'LhessianEV',                      &
     &            'Switch to compute Hessian eigenvectors.'
#    endif
#    ifdef WEAK_CONSTRAINT
          WRITE (out,70) LhotStart, 'LhotStart',                        &
     &            'Switch for hot start of subsequent outer loops.'
#    endif
          WRITE (out,70) Lprecond, 'Lprecond',                          &
     &            'Switch for conjugate gradient preconditioning.'
          WRITE (out,70) Lritz, 'Lritz',                                &
     &            'Switch for Ritz limited-memory preconditioning.'
#    ifdef WEAK_CONSTRAINT
          IF (Lprecond.and.(NritzEV.gt.0)) THEN
            WRITE (out,80) NritzEV, 'NritzEV',                          &
     &            'Number of preconditioning eigenpairs to use.'
          END IF
#    endif
#   endif
#   ifdef BALANCE_OPERATOR
#    ifdef ZETA_ELLIPTIC
          WRITE (out,80) Nbico(ng), 'Nbico',                            &
     &            'Number of iterations in SSH elliptic equation.'
#    endif
          WRITE (out,100) dTdz_min(ng), 'dTdz_min',                     &
     &            'Minimum dTdz (C/m) used in balanced salinity.'
          WRITE (out,100) ml_depth(ng), 'ml_depth',                     &
     &            'Mixed-layer depth (m) used in balanced salinity.'
          IF (balance(isFsur)) THEN
            WRITE (out,100) LNM_depth(ng), 'LNM_depth',                 &
     &            'Level of no motion (m) in balanced free-sruface.'
            WRITE (out,80) LNM_flag, 'LNM_flag',                        &
     &            'Level of no motion integration flag.'
          END IF
          WRITE (out,70) balance(isFsur), 'balance(isFsur)',            &
                  'Switch to include free-surface in balance operator.'
#    ifdef SOLVE3D
          WRITE (out,70) balance(isTvar(isalt)), 'balance(isSalt)',     &
                  'Switch to include salinity in balance operator.'
          WRITE (out,70) balance(isVvel), 'balance(isVvel)',            &
                  'Switch to include 3D momentum in balance operator.'
#    else
          WRITE (out,70) balance(isVbar), 'balance(isVbar)',            &
                  'Switch to include 2D momentum in balance operator.'
#    endif
#   endif
#   ifdef WEAK_CONSTRAINT
#    if defined ARRAY_MODES
          WRITE (out,80) Nvct, 'Nvct',                                  &
     &            'Representer array mode eigenvector to process.'
#    elif defined CLIPPING
          WRITE (out,80) Nvct, 'Nvct',                                  &
     &            'Representer cut-off eigenvector to process.'
#    endif
#   endif
#   if defined POSTERIOR_EOFS && defined WEAK_CONSTRAINT
          WRITE (out,80) NpostI, 'NpostI',                              &
     &            'Number of Lanczos iterations in posterior analysis.'
#   endif
#   ifndef TLM_CHECK
#    ifndef IS4DVAR_SENSITIVITY
          WRITE (out,170) LdefNRM(1:4,ng), 'LdefNRM',                   &
     &            'Switch to create a normalization NetCDF file.'
          WRITE (out,170) LwrtNRM(1:4,ng), 'LwrtNRM',                   &
     &            'Switch to write out normalization factors.'
          IF (ANY(LwrtNRM(:,ng))) THEN
            IF (Nmethod(ng).eq.0) THEN
              WRITE (out,80) Nmethod(ng), 'Nmethod',                    &
     &            'Correlation normalization method: Exact.'
            ELSE IF (Nmethod(ng).eq.1) THEN
              WRITE (out,80) Nmethod(ng), 'Nmethod',                    &
     &            'Correlation normalization method: Randomization.'
              WRITE (out,80) Rscheme(ng), 'Rscheme',                    &
     &            'Random number generation scheme'
              WRITE (out,80) Nrandom, 'Nrandom',                        &
     &            'Number of iterations for randomization.'
            END IF
          END IF
#     if defined SENSITIVITY_4DVAR || \
         defined TL_W4DPSAS        || defined TL_W4DVAR || \
         defined W4DPSAS           || defined W4DVAR
          IF (ANY(LwrtNRM(:,ng))) THEN
            WRITE (out,70) Cnorm(2,isFsur), 'CnormM(isFsur)',           &
     &            'Compute model 2D RHO-normalization factors.'
            WRITE (out,70) Cnorm(2,isUbar), 'CnormM(isUbar)',           &
     &            'Compute model 2D U-normalization factors.'
            WRITE (out,70) Cnorm(2,isVbar), 'CnormM(isVbar)',           &
     &            'Compute model 2D V-normalization factors.'
#      ifdef SOLVE3D
            WRITE (out,70) Cnorm(2,isUvel), 'CnormM(isUvel)',           &
     &            'Compute model 3D U-normalization factors.'
            WRITE (out,70) Cnorm(2,isVvel), 'CnormM(isVvel)',           &
     &            'Compute model 3D V-normalization factors.'
            DO itrc=1,NT(ng)
              WRITE (out,110) Cnorm(2,isTvar(itrc)), 'CnormM(isTvar)',  &
     &            'Compute model normalization factors for tracer ',    &
     &            itrc, TRIM(Vname(1,idTvar(itrc)))
            END DO
#      endif
          END IF
#     endif
          IF (ANY(LwrtNRM(:,ng))) THEN
            WRITE (out,70) Cnorm(1,isFsur), 'CnormI(isFsur)',           &
     &            'Compute initial 2D RHO-normalization factors.'
            WRITE (out,70) Cnorm(1,isUbar), 'CnormI(isUbar)',           &
     &            'Compute initial 2D U-normalization factors.'
            WRITE (out,70) Cnorm(1,isVbar), 'CnormI(isVbar)',           &
     &            'Compute initial 2D V-normalization factors.'
#     ifdef SOLVE3D
            WRITE (out,70) Cnorm(1,isUvel), 'CnormI(isUvel)',           &
     &            'Compute initial 3D U-normalization factors.'
            WRITE (out,70) Cnorm(1,isVvel), 'CnormI(isVvel)',           &
     &            'Compute initial 3D V-normalization factors.'
            DO itrc=1,NT(ng)
              WRITE (out,110) Cnorm(1,isTvar(itrc)), 'CnormI(isTvar)',  &
     &            'Compute initial normalization factors for tracer ',  &
     &            itrc, TRIM(Vname(1,idTvar(itrc)))
            END DO
#     endif
          END IF
#     ifdef ADJUST_BOUNDARY
          IF (ANY(LwrtNRM(:,ng))) THEN
            WRITE (out,170) CnormB(isFsur,1:4), 'CnormB(isFsur)',       &
     &            'Compute boundary 2D RHO-normalization factors.'
            WRITE (out,170) CnormB(isUbar,1:4), 'CnormB(isUbar)',       &
     &            'Compute boundary 2D U-normalization factors.'
            WRITE (out,170) CnormB(isVbar,1:4), 'CnormB(isVbar)',       &
     &            'Compute initial 2D V-normalization factors.'
#      ifdef SOLVE3D
            WRITE (out,170) CnormB(isUvel,1:4), 'CnormB(isUvel)',       &
     &            'Compute initial 3D U-normalization factors.'
            WRITE (out,170) CnormB(isVvel,1:4), 'CnormI(isVvel)',       &
     &            'Compute initial 3D V-normalization factors.'
            DO itrc=1,NT(ng)
              WRITE (out,180) CnormB(isTvar(itrc),1:4),'CnormI(isTvar)',&
     &            'Compute initial normalization factors for tracer ',  &
     &            itrc, TRIM(Vname(1,idTvar(itrc)))
            END DO
#      endif
          END IF
#     endif
#     ifdef ADJUST_WSTRESS
          IF (ANY(LwrtNRM(:,ng))) THEN
            WRITE (out,70) Cnorm(1,isUstr), 'CnormF(isUstr)',           &
     &            'Compute normalization factors at surface U-stress.'
            WRITE (out,70) Cnorm(1,isVstr), 'CnormF(isVstr)',           &
     &            'Compute normalization factors at surface V-stress.'
          END IF
#     endif
#     if defined ADJUST_STFLUX && defined SOLVE3D
          IF (ANY(LwrtNRM(:,ng))) THEN
            DO itrc=1,NT(ng)
              WRITE (out,110) Cnorm(1,isTsur(itrc)), 'CnormF(isTsur)',  &
     &            'Compute normalization factors for flux of tracer ',  &
     &            itrc, TRIM(Vname(1,idTvar(itrc)))
            END DO
          END IF
#     endif
#    endif
          WRITE (out,100) Hgamma(1), 'Hgamma',                          &
     &            'Horizontal diffusion factor, initial conditions.'
#    ifdef WEAK_CONSTRAINT
          WRITE (out,100) Hgamma(2), 'HgammaM',                         &
     &            'Horizontal diffusion factor, model error.'
#    endif
#    ifdef ADJUST_BOUNDARY
          WRITE (out,100) Hgamma(3), 'HgammaB',                         &
     &            'Horizontal diffusion factor, boundary conditions.'
#    endif
#    ifdef ADJUST_STFLUX
          WRITE (out,100) Hgamma(4), 'HgammaF',                         &
     &            'Horizontal diffusion factor, surface forcing.'
#    endif
#    ifdef SOLVE3D
          WRITE (out,100) Vgamma(1), 'Vgamma',                          &
     &            'Vertical diffusion factor, initial conditions.'
#     ifdef WEAK_CONSTRAINT
          WRITE (out,100) Vgamma(2), 'VgammaM',                         &
     &            'Vertical diffusion factor, model error.'
#     endif
#     ifdef ADJUST_BOUNDARY
          WRITE (out,100) Vgamma(3), 'VgammaB',                         &
     &            'Vertical diffusion factor, boundary conditions.'
#     endif
#    endif
#    if defined SENSITIVITY_4DVAR || \
        defined TL_W4DPSAS        || defined TL_W4DVAR || \
        defined W4DPSAS           || defined W4DVAR
         IF (nADJ(ng).lt.ntimes(ng)) THEN
            WRITE (out,120) Hdecay(2,isFsur,ng), 'HdecayM(isFsur)',     &
     &            'Model decorrelation H-scale (m), free-surface.'
            WRITE (out,120) Hdecay(2,isUbar,ng), 'HdecayM(isUbar)',     &
     &            'Model decorrelation H-scale (m), 2D U-momentum.'
            WRITE (out,120) Hdecay(2,isVbar,ng), 'HdecayM(isVbar)',     &
     &            'Model decorrelation H-scale (m), 2D V-momentum.'
#     ifdef SOLVE3D
            WRITE (out,120) Hdecay(2,isUvel,ng), 'HdecayM(isUvel)',     &
     &            'Model decorrelation H-scale (m), 3D U-momentum.'
            WRITE (out,120) Hdecay(2,isVvel,ng), 'HdecayM(isVvel)',     &
     &            'Model decorrelation H-scale (m), 3D V-momentum.'
            DO itrc=1,NT(ng)
              WRITE (out,130) Hdecay(2,isTvar(itrc),ng),                &
     &            'HdecayM(idTvar)',                                    &
     &            'Model decorrelation H-scale (m), ',                  &
     &            TRIM(Vname(1,idTvar(itrc)))
            END DO
            WRITE (out,120) Vdecay(2,isUvel,ng), 'VdecayM(isUvel)',     &
     &            'Model decorrelation V-scale (m), 3D U-momentum.'
            WRITE (out,120) Vdecay(2,isVvel,ng), 'VdecayM(isVvel)',     &
     &            'Model decorrelation V-scale (m), 3D V-momentum.'
            DO itrc=1,NT(ng)
              WRITE (out,130) Vdecay(2,isTvar(itrc),ng),                &
     &            'VdecayM(idTvar)',                                    &
     &            'Model decorrelation V-scale (m), ',                  &
     &            TRIM(Vname(1,idTvar(itrc)))
            END DO
#     endif
          END IF
#    endif
#    if defined WEAK_CONSTRAINT && defined TIME_CONV
          WRITE (out,80) NrecTC(ng), 'NrecTC',                          &
     &            'Number of state records for time convolution.'
          WRITE (out,120) Tdecay(isFsur,ng), 'TdecayM(isFsur)',         &
     &          'Model decorrelation T-scale (day), free-surface.'
          WRITE (out,120) Tdecay(isUbar,ng), 'TdecayM(isUbar)',         &
     &          'Model decorrelation T-scale (day), 2D U-momentum.'
          WRITE (out,120) Tdecay(isVbar,ng), 'TdecayM(isVbar)',         &
     &          'Model decorrelation T-scale (day), 2D V-momentum.'
#     ifdef SOLVE3D
          WRITE (out,120) Tdecay(isUvel,ng), 'TdecayM(isUvel)',         &
     &          'Model decorrelation T-scale (day), 3D U-momentum.'
          WRITE (out,120) Tdecay(isVvel,ng), 'TdecayM(isVvel)',         &
     &          'Model decorrelation T-scale (day), 3D V-momentum.'
          DO itrc=1,NT(ng)
            WRITE (out,130) Tdecay(isTvar(itrc),ng),                    &
     &          'TdecayM(idTvar)',                                      &
     &          'Model decorrelation T-scale (day), ',                  &
     &          TRIM(Vname(1,idTvar(itrc)))
            END DO
#     endif
#    endif
          WRITE (out,120) Hdecay(1,isFsur,ng), 'HdecayI(isFsur)',       &
     &            'Initial decorrelation H-scale (m), free-surface.'
          WRITE (out,120) Hdecay(1,isUbar,ng), 'HdecayI(isUbar)',       &
     &            'Initial decorrelation H-scale (m), 2D U-momentum.'
          WRITE (out,120) Hdecay(1,isVbar,ng), 'HdecayI(isVbar)',       &
     &            'Initial decorrelation H-scale (m), 2D V-momentum.'
#    ifdef SOLVE3D
          WRITE (out,120) Hdecay(1,isUvel,ng), 'HdecayI(isUvel)',       &
     &            'Initial decorrelation H-scale (m), 3D U-momentum.'
          WRITE (out,120) Hdecay(1,isVvel,ng), 'HdecayI(isVvel)',       &
     &            'Initial decorrelation H-scale (m), 3D V-momentum.'
          DO itrc=1,NT(ng)
            WRITE (out,130) Hdecay(1,isTvar(itrc),ng),                  &
     &            'HdecayI(idTvar)',                                    &
     &            'Initial decorrelation H-scale (m), ',                &
     &            TRIM(Vname(1,idTvar(itrc)))
          END DO
          WRITE (out,120) Vdecay(1,isUvel,ng), 'VdecayI(isUvel)',       &
     &            'Initial decorrelation V-scale (m), 3D U-momentum.'
          WRITE (out,120) Vdecay(1,isVvel,ng), 'VdecayI(isVvel)',       &
     &            'Initial decorrelation V-scale (m), 3D V-momentum.'
          DO itrc=1,NT(ng)
            WRITE (out,130) Vdecay(1,isTvar(itrc),ng),                  &
     &            'VdecayI(idTvar)',                                    &
     &            'Initial decorrelation V-scale (m), ',                &
     &            TRIM(Vname(1,idTvar(itrc)))
          END DO
#    endif
#    ifdef ADJUST_BOUNDARY
          DO ib=1,4
            IF (ib.eq.iwest) THEN
              Text='W-bry  '
            ELSE IF (ib.eq.isouth) THEN
              Text='S-bry  '
            ELSE IF (ib.eq.ieast) THEN
              Text='E-bry  '
            ELSE IF (ib.eq.inorth) THEN
              Text='N-bry  '
            END IF
            IF (Lobc(ib,isFsur,ng)) THEN
              WRITE (out,120) HdecayB(isFsur,ib,ng), 'HdecayB(isFsur)', &
     &              Text//' decorrelation H-scale (m), free-surface.'
            END IF
            IF (Lobc(ib,isUbar,ng)) THEN
              WRITE (out,120) HdecayB(isUbar,ib,ng), 'HdecayB(isUbar)', &
     &              Text//' decorrelation H-scale (m), 2D U-momentum.'
            END IF
            IF (Lobc(ib,isVbar,ng)) THEN
              WRITE (out,120) HdecayB(isVbar,ib,ng), 'HdecayB(isVbar)', &
     &              Text//' decorrelation H-scale (m), 2D V-momentum.'
            END IF
#     ifdef SOLVE3D
            IF (Lobc(ib,isUvel,ng)) THEN
              WRITE (out,120) HdecayB(isUvel,ib,ng), 'HdecayB(isUvel)', &
     &              Text//' decorrelation H-scale (m), 3D U-momentum.'
            END IF
            IF (Lobc(ib,isVvel,ng)) THEN
              WRITE (out,120) HdecayB(isVvel,ib,ng), 'HdecayB(isVvel)', &
     &              Text//' decorrelation H-scale (m), 3D V-momentum.'
            END IF
            DO i=1,NT(ng)
              IF (Lobc(ib,isTvar(i),ng)) THEN
                WRITE(out,130) HdecayB(isTvar(i),ib,ng),                &
     &              'HdecayB(idTvar)',                                  &
     &              Text//' decorrelation H-scale (m), ',               &
     &              TRIM(Vname(1,idTvar(i)))
              END IF
            END DO
            IF (Lobc(ib,isUvel,ng)) THEN
              WRITE (out,120) VdecayB(isUvel,ib,ng), 'VdecayB(isUvel)', &
     &              Text//' decorrelation V-scale (m), 3D U-momentum.'
            END IF
            IF (Lobc(ib,isVvel,ng)) THEN
              WRITE (out,120) VdecayB(isVvel,ib,ng), 'VdecayB(isVvel)', &
     &              Text//' decorrelation V-scale (m), 3D V-momentum.'
            END IF
            DO i=1,NT(ng)
              IF (Lobc(ib,isTvar(i),ng)) THEN
                WRITE(out,130) VdecayB(isTvar(i),ib,ng),                &
     &              'VdecayB(idTvar)',                                  &
     &              Text//' decorrelation V-scale (m), ',               &
     &              TRIM(Vname(1,idTvar(i)))
              END IF
            END DO
#     endif
          END DO
#    endif
#    ifdef ADJUST_WSTRESS
          WRITE (out,120) Hdecay(1,isUstr,ng), 'HdecayF(isUstr)',       &
     &            'Forcing decorrelation H-scale (m), U-stress.'
          WRITE (out,120) Hdecay(1,isVstr,ng), 'HdecayF(isVstr)',       &
     &            'Forcing decorrelation H-scale (m), V-stress.'
#    endif
#    if defined ADJUST_STFLUX && defined SOLVE3D
          DO itrc=1,NT(ng)
            WRITE (out,130) Hdecay(1,isTsur(itrc),ng),                  &
     &            'HdecayF(idTsur)',                                    &
     &            'Forcing decorrelation H-scale (m), ',                &
     &            TRIM(Vname(1,idTvar(itrc)))
          END DO
#    endif
#    ifdef BGQC
          IF (bgqc_type(ng).eq.1) THEN
            WRITE (out,80) bgqc_type(ng), 'bgqc_type',                  &
     &            'Quality control in terms of state variable index'
          ELSE IF (bgqc_type(ng).eq.2) THEN
            WRITE (out,80) bgqc_type(ng), 'bgqc_type',                  &
     &            'Quality control in terms of observation provenance'
          END IF
          IF (bgqc_type(ng).eq.1) THEN
            WRITE (out,100) S_bgqc(isFsur,ng), 'S_bgqc(isFsur)',        &
     &            'Quality control reject threshold, free-surface.'
#     ifndef SOLVE3D
            WRITE (out,100) S_bgqc(isUbar,ng), 'S_bgqc(isUbar)',        &
     &            'Quality control reject threshold, 2D U-momentum.'
            WRITE (out,100) S_bgqc(isVbar,ng), 'S_bgqc(isVbar)',        &
     &            'Quality control reject threshold, 2D V-momentum.'
#     else
            WRITE (out,100) S_bgqc(isUvel,ng), 'S_bgqc(isUvel)',        &
     &            'Quality control reject threshold, 3D U-momentum.'
            WRITE (out,100) S_bgqc(isVbar,ng), 'S_bgqc(isVvel)',        &
     &            'Quality control reject threshold, 3D V-momentum.'
            DO itrc=1,NT(ng)
              WRITE (out,190) S_bgqc(isTvar(itrc),ng),                  &
     &            'S_bgqc(idTvar)',                                     &
     &            'Quality control reject threshold, ',                 &
     &            TRIM(Vname(1,idTvar(itrc)))
            END DO
#     endif
          ELSE IF (bgqc_type(ng).eq.2) THEN
            WRITE (out,80) Nprovenance(ng), 'Nprovenance',              &
     &            'Number of provenances to Quality control.'
            DO i=1,Nprovenance(ng)
              WRITE (out,200) P_bgqc(i,ng), 'P_bgqc', i,                &
     &            'Quality control reject threshold for provenance',    &
     &            Iprovenance(i,ng)
            END DO
          END IF
#    endif
#    if defined ADJUST_STFLUX && defined SOLVE3D
          DO itrc=1,NT(ng)
            WRITE (out,110) Lstflux(itrc,ng), 'Lstflux(itrc)',          &
     &            'Adjusting surface flux of tracer ', itrc,            &
     &            TRIM(Vname(1,idTvar(itrc)))
          END DO
#    endif
#    ifdef ADJUST_BOUNDARY
          WRITE (out,170) Lobc(1:4,isFsur,ng), 'Lobc(isFsur)',          &
     &            'Adjusting free-surface boundaries.'
          WRITE (out,170) Lobc(1:4,isUbar,ng), 'Lobc(isUbar)',          &
     &            'Adjusting 2D U-momentum boundaries.'
          WRITE (out,170) Lobc(1:4,isVbar,ng), 'Lobc(isVbar)',          &
     &            'Adjusting 2D V-momentum boundaries.'
#     ifdef SOLVE3D
          WRITE (out,170) Lobc(1:4,isUvel,ng), 'Lobc(isUvel)',          &
     &            'Adjusting 3D U-momentum boundaries.'
          WRITE (out,170) Lobc(1:4,isVvel,ng), 'Lobc(isVvel)',          &
     &            'Adjusting 3D V-momentum boundaries.'
          DO itrc=1,NT(ng)
            WRITE (out,180) Lobc(1:4,isTvar(itrc),ng),'Lobc(isTvar)',   &
     &            'Adjusting boundaries for tracer ', itrc,             &
     &            TRIM(Vname(1,idTvar(itrc)))
          END DO
#     endif
#    endif
#   endif
#  endif
!
!-----------------------------------------------------------------------
!  Report input files and check availability of input files.
!-----------------------------------------------------------------------
!
          WRITE (out,150)
#  ifdef VERIFICATION
          WRITE (out,160) ' Verification Parameters File:  ',           &
     &                    TRIM(aparnam)
#  else
          WRITE (out,160) ' Assimilation Parameters File:  ',           &
     &                    TRIM(aparnam)
#  endif
#  if defined FOUR_DVAR || (defined HESSIAN_SV && defined BNORM)
#   if defined IS4DVAR          || defined OBS_SENSITIBITY || \
       defined OPT_OBSERVATIONS || defined WEAK_CONSTRAINT
#    if defined SENSITIVITY_4DVAR || \
        defined TL_W4DPSAS        || defined TL_W4DVAR || \
        defined W4DPSAS           || defined W4DVAR
          fname=STD(2,ng)%name
          IF (.not.find_file(ng, fname, 'STDnameM')) GO TO 30
          WRITE (out,160) '               Model STD File:  ',           &
     &                    TRIM(fname)
#    endif
          fname=STD(1,ng)%name
          IF (.not.find_file(ng, fname, 'STDnameI')) GO TO 30
          WRITE (out,160) '  Initial Conditions STD File:  ',           &
     &                    TRIM(fname)
#    ifdef ADJUST_BOUNDARY
          fname=STD(3,ng)%name
          IF (.not.find_file(ng, fname, 'STDnameB')) GO TO 30
          WRITE (out,160) ' Boundary Conditions STD File:  ',           &
     &                    TRIM(fname)
#    endif
#    if defined ADJUST_WSTRESS || defined ADJUST_STFLUX
          fname=STD(4,ng)%name
          IF (.not.find_file(ng, fname, 'STDnameF')) GO TO 30
          WRITE (out,160) '     Surface Forcing STD File:  ',           &
     &                    TRIM(fname)
#    endif
#   endif
#   if defined SENSITIVITY_4DVAR || defined TL_W4DPSAS  || \
       defined TL_W4DVAR         || defined W4DPSAS     || \
       defined W4DVAR
          fname=NRM(2,ng)%name
          WRITE (out,160) '              Model Norm File:  ',           &
     &                    TRIM(fname)
          IF (.not.LdefNRM(2,ng)) THEN
            IF (.not.find_file(ng, fname, 'NRMnameM')) GO TO 30
          END IF
#   elif defined CORRELATION
          fname=NRM(2,ng)%name
          WRITE (out,160) '              Model Norm File:  ',           &
     &                    TRIM(fname)
          IF (.not.LdefNRM(2,ng).and.LwrtNRM(2,ng)) THEN
            IF (.not.find_file(ng, fname, 'NRMnameM')) GO TO 30
          END IF
#   endif
#   if defined CORRELATION
          fname=NRM(1,ng)%name
          WRITE (out,160) ' Initial Conditions Norm File:  ',           &
     &                    TRIM(fname)
          IF (.not.LdefNRM(1,ng).and.LwrtNRM(1,ng)) THEN
            IF (.not.find_file(ng, fname, 'NRMnameI')) GO TO 30
          END IF
#   else
          fname=NRM(1,ng)%name
          WRITE (out,160) ' Initial Conditions Norm File:  ',           &
     &                    TRIM(fname)
          IF (.not.LdefNRM(1,ng)) THEN
            IF (.not.find_file(ng, fname, 'NRMnameI')) GO TO 30
          END IF
#   endif
#   ifdef ADJUST_BOUNDARY
#    ifdef CORRELATION
          fname=NRM(3,ng)%name
          WRITE (out,160) 'Boundary Conditions Norm File:  ',           &
     &                    TRIM(fname)
          IF (.not.LdefNRM(3,ng).and.LwrtNRM(3,ng)) THEN
            IF (.not.find_file(ng, fname, 'NRMnameB')) GO TO 30
          END IF
#    else
          fname=NRM(3,ng)%name
          WRITE (out,160) 'Boundary Conditions Norm File:  ',           &
     &                    TRIM(fname)
          IF (.not.LdefNRM(3,ng)) THEN
            IF (.not.find_file(ng, fname, 'NRMnameB')) GO TO 30
          END IF
#    endif
#   endif
#   if defined ADJUST_WSTRESS || defined ADJUST_STFLUX
#    ifdef CORRELATION
          fname=NRM(4,ng)%name
          WRITE (out,160) '    Surface Forcing Norm File:  ',           &
     &                    TRIM(fname)
          IF (.not.LdefNRM(4,ng).and.LwrtNRM(4,ng)) THEN
            IF (.not.find_file(ng, fname, 'NRMnameF')) GO TO 30
          END IF
#    else
          fname=NRM(4,ng)%name
          WRITE (out,160) '    Surface Forcing Norm File:  ',           &
     &                    TRIM(fname)
          IF (.not.LdefNRM(4,ng)) THEN
            IF (.not.find_file(ng, fname, 'NRMnameF')) GO TO 30
          END IF
#    endif
#   endif
#   if !(defined CORRELATION || defined OPT_OBSERVATIONS)
          fname=OBS(ng)%name
          IF (.not.find_file(ng, fname, 'OBSname')) GO TO 30
          WRITE (out,160) '            Observations File:  ',           &
     &                    TRIM(fname)
#   endif
#  endif
          WRITE (out,160) '    Input/Output Lanczos File:  ',           &
     &                    TRIM(LCZ(ng)%name)
#  ifndef IS4DVAR_SENSITIVITY
          WRITE (out,160) '    Input/Output Hessian File:  ',           &
     &                    TRIM(HSS(ng)%name)
#  endif
#  ifdef EVOLVED_LCZ
          WRITE (out,160) '  Output evolved Lanczos File:  ',           &
     &                    TRIM(LZE(ng)%name)
#  endif
#  ifdef VERIFICATION
          WRITE (out,160) '     Output verification File:  ',           &
     &                    TRIM(DAV(ng)%name)
#  else
          WRITE (out,160) '            Output 4DVAR File:  ',           &
     &                    TRIM(DAV(ng)%name)
#  endif
#  if defined WEAK_CONSTRAINT   && \
     (defined POSTERIOR_ERROR_F || defined POSTERIOR_ERROR_I)
          WRITE (out,160) '  Output Posterior Error File:  ',           &
     &                    TRIM(ERR(ng)%name)
#  endif
          GO TO 40
  30      WRITE (out,220) TRIM(fname)
          exit_flag=4
          RETURN
  40      CONTINUE
        END DO
      END IF

#  if defined WEAK_CONSTRAINT && defined TIME_CONV
!
!-----------------------------------------------------------------------
!  Convert time decorrelation scales to seconds.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
        DO i=1,MstateVar
          Tdecay(i,ng)=Tdecay(i,ng)*86400.0_r8
        END DO
      END DO
#  endif
!
  50  FORMAT (/,' READ_AssPar - Error while processing line: ',/,a)
  55  FORMAT (/,' READ_AssPar - ',a,i4,2x,i4,/,15x,a)
  60  FORMAT (/,/,' Assimilation Parameters, Grid: ',i2.2,              &
     &        /,  ' =================================',/)
  70  FORMAT (10x,l1,2x,a,t30,a)
  80  FORMAT (1x,i10,2x,a,t30,a)
  90  FORMAT (1x,i10,2x,a,t30,a,/,t32,a)
 100  FORMAT (1p,e11.4,2x,a,t30,a)
 110  FORMAT (10x,l1,2x,a,t30,a,i2.2,':',1x,a)
 120  FORMAT (f11.3,2x,a,t30,a)
 130  FORMAT (f11.3,2x,a,t30,a,a,'.')
 150  FORMAT (/,' Input Assimilation Files:',/)
 160  FORMAT (2x,a,a)
 170  FORMAT (3x,4(1x,l1),2x,a,t30,a)
 180  FORMAT (3x,4(1x,l1),2x,a,t30,a,i2.2,':',1x,a)
 190  FORMAT (1p,e11.4,2x,a,t30,a,a,'.')
 200  FORMAT (1p,e11.4,2x,a,'(',i2.2,')',t30,a,':',1x,i6)
 210  FORMAT (/,' READ_ASSPAR - variable info not yet loaded, ', a)
 220  FORMAT (/,' READ_ASSPAR - could not find input file:  ',a)
# endif

      RETURN
      END SUBROUTINE read_AssPar
#else
      SUBROUTINE read_AssPar
      END SUBROUTINE read_AssPar
#endif
