#include "cppdefs.h"
#if defined FOUR_DVAR || defined VERIFICATION
      SUBROUTINE read_AssPar (model, inp, out, Lwrite)
!
!git $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2026 The ROMS Group                              !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.md                                               !
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
      USE inp_decode_mod
!
# ifdef DISTRIBUTE
      USE distribute_mod, ONLY : mp_bcasti
# endif
      USE strings_mod,    ONLY : FoundError
      USE strings_mod,    ONLY : uppercase
!
      implicit none
!
!  Imported variable declarations
!
      logical, intent(in) :: Lwrite
!
      integer, intent(in) :: model, inp, out
!
!  Local variable declarations.
!
      logical :: Lvalue(4)

      integer :: Mval, Npts, Nval
      integer :: i, ib, igrid, itrc, k, ng, ns, status
      integer :: Cdim, Clen, Mscale, Rdim
      integer :: obs_lib
      integer :: Ivalue(1)
!
      integer :: Nfiles(Ngrids)
      integer :: delta2(2,Ngrids), delta3(3,Ngrids)
# ifdef SOLVE3D
      integer :: delta2T(2,MT,Ngrids), delta3T(3,MT,Ngrids)
# endif
# ifdef MULTI_SCALE_B
      integer, allocatable :: Miter(:,:)
#  ifdef SOLVE3D
      integer, allocatable :: MiterT(:,:,:)
#  endif
# endif
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
      real(r8), dimension(4,MT,Ngrids) :: Rboundary
#   endif
#  endif
      real(r8), dimension(MT,Ngrids) :: Rtracer
# endif
      real(dp), dimension(nRval) :: Rval
      real(r8) :: Rvalue(1)
# ifdef MULTI_SCALE_B
      real(r8), allocatable :: Bscale(:,:)
      real(r8), allocatable :: BscaleB(:,:,:)
#  ifdef SOLVE3D
      real(r8), allocatable :: BscaleBT(:,:,:,:)
      real(r8), allocatable :: BscaleT(:,:,:)
#  endif
# endif
!
      character (len=1 ), parameter :: blank = ' '
# ifdef ADJUST_BOUNDARY
      character (len=7) :: Text
# endif
      character (len=40 ) :: KeyWord
      character (len=50 ) :: label
      character (len=256) :: fname, line
      character (len=256), dimension(nCval) :: Cval

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__
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
      Nval=0
      delta2=0
      delta3=0
# ifdef SOLVE3D
      delta2T=0
      delta3T=0
# endif 
!
!  Set NetCDF library to process reading and writing the observation
!  vectors. Since the input OBS(ng)%name file can be rewritten, we
!  need to ensure that the output DAV(ng)%name is processed with the
!  same library.
!
      IF (inp_lib.ne.out_lib) THEN
        obs_lib=io_nf90             ! use standard NetCDF library
      ELSE
        obs_lib=inp_lib             ! use specified input NetCDF library
      END IF

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
            CASE ('Nimpact')
              Npts=load_i(Nval, Rval, 1, Ivalue)
              Nimpact=Ivalue(1)
#  if defined SPLIT_4DVAR
            CASE ('OuterLoop')
              Npts=load_i(Nval, Rval, 1, Ivalue)
              OuterLoop=Ivalue(1)
            CASE ('Phase4DVAR')
              DO i=1,LEN(Phase4DVAR)
                Phase4DVAR(i:i)=blank
              END DO
              Phase4DVAR=TRIM(ADJUSTL(Cval(Nval)))
#  endif
            CASE ('NextraObs')
              Npts=load_i(Nval, Rval, 1, Ivalue)
              NextraObs=Ivalue(1)
              IF (NextraObs.gt.0) THEN
                IF (.not.allocated(ExtraIndex)) THEN
                  allocate ( ExtraIndex(NextraObs) )
                END IF
                IF (.not.allocated(ExtraName)) THEN
                  allocate ( ExtraName(NextraObs) )
                END IF
              END IF
            CASE ('ExtraIndex')
              IF (NextraObs.gt.0) THEN
                Npts=load_i(Nval, Rval, NextraObs, ExtraIndex)
              END IF
            CASE ('ExtraName')
              IF (NextraObs.gt.0) THEN
                ExtraName(Nval)=TRIM(Cval(Nval))
              END IF
            CASE ('tl_M2diff')
              Npts=load_r(Nval, Rval, Ngrids, tl_M2diff)
            CASE ('tl_M3diff')
              Npts=load_r(Nval, Rval, Ngrids, tl_M3diff)
            CASE ('tl_Tdiff')
              Npts=load_r(Nval, Rval, MT, Ngrids, Rtracer)
              DO ng=1,Ngrids
                DO itrc=1,NT(ng)
                  tl_Tdiff(itrc,ng)=Rtracer(itrc,ng)
                END DO
              END DO
            CASE ('timeIAU')
              Npts=load_r(Nval, Rval, Ngrids, timeIAU)
            CASE ('LdefNRM')
              Npts=load_l(Nval, Cval, 4, Ngrids, LdefNRM)
            CASE ('LwrtNRM')
              Npts=load_l(Nval, Cval, 4, Ngrids, LwrtNRM)
            CASE ('CnormM(isFsur)')
              IF (isFsur.eq.0) THEN
                IF (Master) WRITE (out,210) 'isFsur'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, 1, Lvalue)
              Cnorm(2,isFsur)=Lvalue(1)
            CASE ('CnormM(isUbar)')
              IF (isUbar.eq.0) THEN
                IF (Master) WRITE (out,210) 'isUbar'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, 1, Lvalue)
              Cnorm(2,isUbar)=Lvalue(1)
            CASE ('CnormM(isVbar)')
              IF (isVbar.eq.0) THEN
                IF (Master) WRITE (out,210) 'isVbar'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, 1, Lvalue)
              Cnorm(2,isVbar)=Lvalue(1)
#  ifdef SOLVE3D
            CASE ('CnormM(isUvel)')
              IF (isUvel.eq.0) THEN
                IF (Master) WRITE (out,210) 'isUvel'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, 1, Lvalue)
              Cnorm(2,isUvel)=Lvalue(1)
            CASE ('CnormM(isVvel)')
              IF (isVvel.eq.0) THEN
                IF (Master) WRITE (out,210) 'isVvel'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, 1, Lvalue)
              Cnorm(2,isVvel)=Lvalue(1)
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
              Npts=load_l(Nval, Cval, 1, Lvalue)
              Cnorm(1,isFsur)=Lvalue(1)
            CASE ('CnormI(isUbar)')
              Npts=load_l(Nval, Cval, 1, Lvalue)
              Cnorm(1,isUbar)=Lvalue(1)
            CASE ('CnormI(isVbar)')
              Npts=load_l(Nval, Cval, 1, Lvalue)
              Cnorm(1,isVbar)=Lvalue(1)
#  ifdef SOLVE3D
            CASE ('CnormI(isUvel)')
              Npts=load_l(Nval, Cval, 1, Lvalue)
              Cnorm(1,isUvel)=Lvalue(1)
            CASE ('CnormI(isVvel)')
              Npts=load_l(Nval, Cval, 1, Lvalue)
              Cnorm(1,isVvel)=Lvalue(1)
            CASE ('CnormI(isTvar)')
              Npts=load_l(Nval, Cval, MT, Ltracer)
              DO itrc=1,MT
                i=isTvar(itrc)
                Cnorm(1,i)=Ltracer(itrc)
              END DO
#  endif
#  ifdef ADJUST_BOUNDARY
            CASE ('CnormB(isFsur)')
              Npts=load_l(Nval, Cval, 4, Lvalue)
              CnormB(isFsur,1:4)=Lvalue(1:4)
            CASE ('CnormB(isUbar)')
              Npts=load_l(Nval, Cval, 4, Lvalue)
              CnormB(isUbar,1:4)=Lvalue(1:4)
            CASE ('CnormB(isVbar)')
              Npts=load_l(Nval, Cval, 4, Lvalue)
              CnormB(isVbar,1:4)=Lvalue(1:4)
#   ifdef SOLVE3D
            CASE ('CnormB(isUvel)')
              Npts=load_l(Nval, Cval, 4, Lvalue)
              CnormB(isUvel,1:4)=Lvalue(1:4)
            CASE ('CnormB(isVvel)')
              Npts=load_l(Nval, Cval, 4, Lvalue)
              CnormB(isVvel,1:4)=Lvalue(1:4)
            CASE ('CnormB(isTvar)')
              Npts=load_l(Nval, Cval, MT, 4, Lboundary)
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
              Npts=load_l(Nval, Cval, 1, Lvalue)
              Cnorm(1,isUstr)=Lvalue(1)
            CASE ('CnormF(isVstr)')
              IF (isVstr.eq.0) THEN
                IF (Master) WRITE (out,210) 'isVstr'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, 1, Lvalue)
              Cnorm(1,isVstr)=Lvalue(1)
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
#  ifdef SALINITY
            CASE ('balance(isSalt)')
              Npts=load_l(Nval, Cval, 1, Lvalue)
              balance(isTvar(isalt))=Lvalue(1)
#  endif
            CASE ('balance(isFsur)')
              Npts=load_l(Nval, Cval, 1, Lvalue)
              balance(isFsur)=Lvalue(1)
            CASE ('balance(isVbar)')
              Npts=load_l(Nval, Cval, 1, Lvalue)
              balance(isVbar)=Lvalue(1)
            CASE ('balance(isVvel)')
              Npts=load_l(Nval, Cval, 1, Lvalue)
              balance(isVvel)=Lvalue(1)
            CASE ('Nmethod')
              Npts=load_i(Nval, Rval, Ngrids, Nmethod)
            CASE ('Rscheme')
              Npts=load_i(Nval, Rval, Ngrids, Rscheme)
            CASE ('Nrandom')
              Npts=load_i(Nval, Rval, 1, Ivalue)
              Nrandom=Ivalue(1)
            CASE ('Dirac(isFsur)')
              Npts=load_i(Nval, Rval, 2, Ngrids, delta2)
              DO ng=1,Ngrids
                DO i=1,2
                  Dirac(i,isFsur,ng)=delta2(i,ng)
                END DO
              END DO            
            CASE ('Dirac(isUbar)')
              Npts=load_i(Nval, Rval, 2, Ngrids, delta2) 
              DO ng=1,Ngrids
                DO i=1,2
                  Dirac(i,isUbar,ng)=delta2(i,ng)
                END DO
              END DO            
            CASE ('Dirac(isVbar)')
              Npts=load_i(Nval, Rval, 2, Ngrids, delta2)
              DO ng=1,Ngrids
                DO i=1,2
                  Dirac(i,isVbar,ng)=delta2(i,ng)
                END DO
              END DO            
#  ifdef SOLVE3D
            CASE ('Dirac(isUvel)')
              Npts=load_i(Nval, Rval, 3, Ngrids, delta3)
              DO ng=1,Ngrids
                DO i=1,3
                  Dirac(i,isUvel,ng)=delta3(i,ng)
                END DO
              END DO            
            CASE ('Dirac(isVvel)')
              Npts=load_i(Nval, Rval, 3, Ngrids, delta3)
              DO ng=1,Ngrids
                DO i=1,3
                  Dirac(i,isVvel,ng)=delta3(i,ng)
                END DO
              END DO            
            CASE ('Dirac(isTvar)')
              Npts=load_i(Nval, Rval, 3, MT, Ngrids, delta3T)
              DO ng=1,Ngrids
                DO itrc=1,NT(ng)
                  k=isTvar(itrc)
                  DO i=1,3
                    Dirac(i,k,ng)=delta3T(i,itrc,ng)
                  END DO
                END DO
              END DO
#  endif
#  ifdef ADJUST_WSTRESS
            CASE ('Dirac(isUstr)')
              Npts=load_i(Nval, Rval, 2, Ngrids, delta2)
              DO ng=1,Ngrids
                DO i=1,2
                  Dirac(i,isUstr,ng)=delta2(i,ng)
                END DO
              END DO            
            CASE ('Dirac(isVstr)')
              Npts=load_i(Nval, Rval, 2, Ngrids, delta2)
              DO ng=1,Ngrids
                DO i=1,2
                  Dirac(i,isVstr,ng)=delta2(i,ng)
                END DO
              END DO            
#  endif
#  if defined ADJUST_STFLUX && defined SOLVE3D
            CASE ('Dirac(isTsur)')
              Npts=load_i(Nval, Rval, 2, MT, Ngrids, delta2T)
              DO ng=1,Ngrids
                DO itrc=1,NT(ng)
                  k=isTsur(itrc)
                  DO i=1,2
                    Dirac(i,k,ng)=delta2T(i,itrc,ng)
                  END DO
                END DO
              END DO
#  endif
            CASE ('Hgamma')
              Npts=load_r(Nval, Rval, 4, Hgamma)
#  ifdef SOLVE3D
            CASE ('Vgamma')
              Npts=load_r(Nval, Rval, 4, Vgamma)
#  endif
#  ifdef MULTI_SCALE_B
            CASE ('Nscale')
              Npts=load_i(Nval, Rval, Ngrids, Nscale)
              Mscale=MAXVAL(Nscale)
!                                                  local arrays
              IF (.not.allocated(Miter)) THEN
                allocate ( Miter(Mscale,Ngrids) )
                Miter=0
              END IF
              IF (.not.allocated(MiterT)) THEN
                allocate ( MiterT(MT,Mscale,Ngrids) )
                MiterT=0
              END IF
              IF (.not.allocated(Bscale)) THEN
                allocate ( Bscale(Mscale,Ngrids) )
                Bscale=0.0_r8
              END IF
              IF (.not.allocated(BscaleT)) THEN
                allocate ( BscaleT(MT,Mscale,Ngrids) )
                BscaleT=0.0_r8
              END IF
              IF (.not.allocated(BscaleB)) THEN
                allocate ( BscaleB(Mscale,4,Ngrids) )
                BscaleB=0.0_r8
              END IF
              IF (.not.allocated(BscaleBT)) THEN
                allocate ( BscaleBT(Mscale,4,MT,Ngrids) )
                BscaleBT=0.0_r8
              END IF
!                                                  global arrays
              IF (.not.allocated(Nlap)) THEN
                allocate ( Nlap(MstateVar,Mscale,Ngrids) )
                Nlap=0
              END IF
              IF (.not.allocated(Bwgt)) THEN
                allocate ( Bwgt(MstateVar,Mscale,Ngrids) )
                Bwgt=0.0_r8
              END IF
              IF (.not.allocated(NiterCG)) THEN
                allocate ( NiterCG(Mscale,Ngrids) )
                NiterCG=0
              END IF
              IF (.not.allocated(NiterCI)) THEN
                allocate ( NiterCI(Mscale,Ngrids) )
                NiterCI=0
              END IF
              IF (.not.allocated(HdecayX)) THEN
                allocate ( HdecayX(2,MstateVar,Mscale,Ngrids) )
                HdecayX=0.0_r8
              END IF
              IF (.not.allocated(HdecayY)) THEN
                allocate ( HdecayY(2,MstateVar,Mscale,Ngrids) )
                HdecayY=0.0_r8
              END IF
              IF (.not.allocated(HdecayB)) THEN
                allocate ( HdecayB(MstateVar,Mscale,4,Ngrids) )
                HdecayB=0.0_r8
              END IF
            CASE ('Mlap(isFsur)')
              Npts=load_i(Nval, Rval, Mscale, Ngrids, Miter)
              DO ng=1,Ngrids
                DO ns=1,Nscale(ng)
                  Nlap(isFsur,ns,ng)=Miter(ns,ng)
                END DO
              END DO
            CASE ('Mlap(isUbar)')
              Npts=load_i(Nval, Rval, Mscale, Ngrids, Miter)
              DO ng=1,Ngrids
                DO ns=1,Nscale(ng)
                  Nlap(isUbar,ns,ng)=Miter(ns,ng)
                END DO
              END DO
            CASE ('Mlap(isVbar)')
              Npts=load_i(Nval, Rval, Mscale, Ngrids, Miter)
              DO ng=1,Ngrids
                DO ns=1,Nscale(ng)
                  Nlap(isVbar,ns,ng)=Miter(ns,ng)
                END DO
              END DO
#   ifdef SOLVE3D
            CASE ('Mlap(isUvel)')
              Npts=load_i(Nval, Rval, Mscale, Ngrids, Miter)
              DO ng=1,Ngrids
                DO ns=1,Nscale(ng)
                  Nlap(isUvel,ns,ng)=Miter(ns,ng)
                END DO
              END DO
            CASE ('Mlap(isVvel)')
              Npts=load_i(Nval, Rval, Mscale, Ngrids, Miter)
              DO ng=1,Ngrids
                DO ns=1,Nscale(ng)
                  Nlap(isVvel,ns,ng)=Miter(ns,ng)
                END DO
              END DO
            CASE ('Mlap(isTvar)')
              Npts=load_i(Nval, Rval, MT, Mscale, Ngrids, MiterT)
              DO ng=1,Ngrids
                DO ns=1,Nscale(ng)
                  DO itrc=1,NT(ng)
                    i=isTvar(itrc)
                    Nlap(i,ns,ng)=MiterT(itrc,ns,ng)
                  END DO
                END DO
              END DO
#   endif
#   ifdef ADJUST_WSTRESS
            CASE ('Mlap(isUstr)')
              Npts=load_i(Nval, Rval, Mscale, Ngrids, Miter)
              DO ng=1,Ngrids
                DO ns=1,Nscale(ng)
                  Nlap(isUstr,ns,ng)=Miter(ns,ng)
                END DO
              END DO
            CASE ('Mlap(isVstr)')
              Npts=load_i(Nval, Rval, Mscale, Ngrids, Miter)
              DO ng=1,Ngrids
                DO ns=1,Nscale(ng)
                  Nlap(isVstr,ns,ng)=Miter(ns,ng)
                END DO
              END DO
#   endif
#   if defined ADJUST_STFLUX && defined SOLVE3D
            CASE ('Mlap(isTsur)')
              Npts=load_i(Nval, Rval, MT, Mscale, Ngrids, MiterT)
              DO ng=1,Ngrids
                DO ns=1,Nscale(ng)
                  DO itrc=1,NT(ng)
                    i=isTsur(itrc)
                    Nlap(i,ns,ng)=MiterT(itrc,ns,ng)
                  END DO
                END DO
              END DO
#   endif
            CASE ('Bwgt(isFsur)')
              Npts=load_r(Nval, Rval, Mscale, Ngrids, Bscale)
              DO ng=1,Ngrids
                DO ns=1,Nscale(ng)
                  Bwgt(isFsur,ns,ng)=Bscale(ns,ng)
                END DO
              END DO
            CASE ('Bwgt(isUbar)')
              Npts=load_r(Nval, Rval, Mscale, Ngrids, Bscale)
              DO ng=1,Ngrids
                DO ns=1,Nscale(ng)
                  Bwgt(isUbar,ns,ng)=Bscale(ns,ng)
                END DO
              END DO
            CASE ('Bwgt(isVbar)')
              Npts=load_r(Nval, Rval, Mscale, Ngrids, Bscale)
              DO ng=1,Ngrids
                DO ns=1,Nscale(ng)
                  Bwgt(isVbar,ns,ng)=Bscale(ns,ng)
                END DO
              END DO
#   ifdef SOLVE3D
            CASE ('Bwgt(isUvel)')
              Npts=load_r(Nval, Rval, Mscale, Ngrids, Bscale)
              DO ng=1,Ngrids
                DO ns=1,Nscale(ng)
                  Bwgt(isUvel,ns,ng)=Bscale(ns,ng)
                END DO
              END DO
            CASE ('Bwgt(isVvel)')
              Npts=load_r(Nval, Rval, Mscale, Ngrids, Bscale)
              DO ng=1,Ngrids
                DO ns=1,Nscale(ng)
                  Bwgt(isVvel,ns,ng)=Bscale(ns,ng)
                END DO
              END DO
            CASE ('Bwgt(isTvar)')
              Npts=load_r(Nval, Rval, MT, Mscale, Ngrids, BscaleT)
              DO ng=1,Ngrids
                DO ns=1,Nscale(ng)
                  DO itrc=1,NT(ng)
                    i=isTvar(itrc)
                    Bwgt(i,ns,ng)=BscaleT(itrc,ns,ng)
                  END DO
                END DO
              END DO
#   endif
#   ifdef ADJUST_WSTRESS
            CASE ('Bwgt(isUstr)')
              Npts=load_r(Nval, Rval, Mscale, Ngrids, Bscale)
              DO ng=1,Ngrids
                DO ns=1,Nscale(ng)
                  Bwgt(isUstr,ns,ng)=Bscale(ns,ng)
                END DO
              END DO
            CASE ('Bwgt(isVstr)')
              Npts=load_r(Nval, Rval, Mscale, Ngrids, Bscale)
              DO ng=1,Ngrids
                DO ns=1,Nscale(ng)
                  Bwgt(isVstr,ns,ng)=Bscale(ns,ng)
                END DO
              END DO
#   endif
#   if defined ADJUST_STFLUX && defined SOLVE3D
            CASE ('Bwgt(isTsur)')
              Npts=load_r(Nval, Rval, MT, Mscale, Ngrids, BscaleT)
              DO ng=1,Ngrids
                DO ns=1,Nscale(ng)
                  DO itrc=1,NT(ng)
                    i=isTsur(itrc)
                    Bwgt(i,ns,ng)=BscaleT(itrc,ns,ng)
                  END DO
                END DO
              END DO
#   endif
            CASE ('NiterCG')
              Npts=load_i(Nval, Rval, Mscale, Ngrids, Miter)
              DO ng=1,Ngrids
                DO ns=1,Nscale(ng)
                  NiterCG(ns,ng)=Miter(ns,ng)
                END DO
              END DO
            CASE ('NiterCI')
              Npts=load_i(Nval, Rval, Mscale, Ngrids, Miter)
              DO ng=1,Ngrids
                DO ns=1,Nscale(ng)
                  NiterCI(ns,ng)=Miter(ns,ng)
                END DO
              END DO
#  endif
#  ifdef MULTI_SCALE_B
            CASE ('HdecayMX(isFsur)')
              Npts=load_r(Nval, Rval, Mscale, Ngrids, Bscale)
              DO ng=1,Ngrids
                DO ns=1,Nscale(ng)
                  HdecayX(2,isFsur,ns,ng)=Bscale(ns,ng)
                END DO
              END DO
            CASE ('HdecayMY(isFsur)')
              Npts=load_r(Nval, Rval, Mscale, Ngrids, Bscale)
              DO ng=1,Ngrids
                DO ns=1,Nscale(ng)
                  HdecayY(2,isFsur,ns,ng)=Bscale(ns,ng)
                END DO
              END DO
            CASE ('HdecayMX(isUbar)')
              Npts=load_r(Nval, Rval, Mscale, Ngrids, Bscale)
              DO ng=1,Ngrids
                DO ns=1,Nscale(ng)
                  HdecayX(2,isUbar,ns,ng)=Bscale(ns,ng)
                END DO
              END DO
            CASE ('HdecayMY(isUbar)')
              Npts=load_r(Nval, Rval, Mscale, Ngrids, Bscale)
              DO ng=1,Ngrids
                DO ns=1,Nscale(ng)
                  HdecayY(2,isUbar,ns,ng)=Bscale(ns,ng)
                END DO
              END DO
            CASE ('HdecayMX(isVbar)')
              Npts=load_r(Nval, Rval, Mscale, Ngrids, Bscale)
              DO ng=1,Ngrids
                DO ns=1,Nscale(ng)
                  HdecayX(2,isVbar,ns,ng)=Bscale(ns,ng)
                END DO
              END DO
            CASE ('HdecayMY(isVbar)')
              Npts=load_r(Nval, Rval, Mscale, Ngrids, Bscale)
              DO ng=1,Ngrids
                DO ns=1,Nscale(ng)
                  HdecayY(2,isVbar,ns,ng)=Bscale(ns,ng)
                END DO
              END DO
#  else
            CASE ('HdecayM(isFsur)')
              Npts=load_r(Nval, Rval, Ngrids, Hdecay(2,isFsur,:))
            CASE ('HdecayM(isUbar)')
              Npts=load_r(Nval, Rval, Ngrids, Hdecay(2,isUbar,:))
            CASE ('HdecayM(isVbar)')
              Npts=load_r(Nval, Rval, Ngrids, Hdecay(2,isVbar,:))
#  endif
#  ifdef SOLVE3D
#   ifdef MULTI_SCALE_B
            CASE ('HdecayMX(isUvel)')
              Npts=load_r(Nval, Rval, Mscale, Ngrids, Bscale)
              DO ng=1,Ngrids
                DO ns=1,Nscale(ng)
                  HdecayX(2,isUvel,ns,ng)=Bscale(ns,ng)
                END DO
              END DO
            CASE ('HdecayMY(isUvel)')
              Npts=load_r(Nval, Rval, Mscale, Ngrids, Bscale)
              DO ng=1,Ngrids
                DO ns=1,Nscale(ng)
                  HdecayY(2,isUvel,ns,ng)=Bscale(ns,ng)
                END DO
              END DO
            CASE ('HdecayMX(isVvel)')
              Npts=load_r(Nval, Rval, Mscale, Ngrids, Bscale)
              DO ng=1,Ngrids
                DO ns=1,Nscale(ng)
                  HdecayX(2,isVvel,ns,ng)=Bscale(ns,ng)
                END DO
              END DO
            CASE ('HdecayMY(isVvel)')
              Npts=load_r(Nval, Rval, Mscale, Ngrids, Bscale)
              DO ng=1,Ngrids
                DO ns=1,Nscale(ng)
                  HdecayY(2,isVvel,ns,ng)=Bscale(ns,ng)
                END DO
              END DO
            CASE ('HdecayMX(isTvar)')
              Npts=load_r(Nval, Rval, MT, Mscale, Ngrids, BscaleT)
              DO ng=1,Ngrids
                DO ns=1,Nscale(ng)
                  DO itrc=1,NT(ng)
                    i=isTvar(itrc)
                    HdecayX(2,i,ns,ng)=BscaleT(itrc,ns,ng)
                  END DO
                END DO
              END DO
            CASE ('HdecayMY(isTvar)')
              Npts=load_r(Nval, Rval, MT, Mscale, Ngrids, BscaleT)
              DO ng=1,Ngrids
                DO ns=1,Nscale(ng)
                  DO itrc=1,NT(ng)
                    i=isTvar(itrc)
                    HdecayY(2,i,ns,ng)=BscaleT(itrc,ns,ng)
                  END DO
                END DO
              END DO
#   else
            CASE ('HdecayM(isUvel)')
              Npts=load_r(Nval, Rval, Ngrids, Hdecay(2,isUvel,:))
            CASE ('HdecayM(isVvel)')
              Npts=load_r(Nval, Rval, Ngrids, Hdecay(2,isVvel,:))
            CASE ('HdecayM(isTvar)')
              Npts=load_r(Nval, Rval, MT, Ngrids, Rtracer)
              DO ng=1,Ngrids
                DO itrc=1,NT(ng)
                  Hdecay(2,isTvar(itrc),ng)=Rtracer(itrc,ng)
                END DO
              END DO
#   endif
            CASE ('VdecayM(isUvel)')
              Npts=load_r(Nval, Rval, Ngrids, Vdecay(2,isUvel,:))
            CASE ('VdecayM(isVvel)')
              Npts=load_r(Nval, Rval, Ngrids, Vdecay(2,isVvel,:))
            CASE ('VdecayM(isTvar)')
              Npts=load_r(Nval, Rval, MT, Ngrids, Rtracer)
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
              Npts=load_r(Nval, Rval, MT, Ngrids, Rtracer)
              DO ng=1,Ngrids
                DO itrc=1,NT(ng)
                  Tdecay(isTvar(itrc),ng)=Rtracer(itrc,ng)
                END DO
              END DO
#  endif
#  ifdef MULTI_SCALE_B
            CASE ('HdecayIX(isFsur)')
              Npts=load_r(Nval, Rval, Mscale, Ngrids, Bscale)
              DO ng=1,Ngrids
                DO ns=1,Nscale(ng)
                  HdecayX(1,isFsur,ns,ng)=Bscale(ns,ng)
                END DO
              END DO
            CASE ('HdecayIY(isFsur)')
              Npts=load_r(Nval, Rval, Mscale, Ngrids, Bscale)
              DO ng=1,Ngrids
                DO ns=1,Nscale(ng)
                  HdecayY(1,isFsur,ns,ng)=Bscale(ns,ng)
                END DO
              END DO
            CASE ('HdecayIX(isUbar)')
              Npts=load_r(Nval, Rval, Mscale, Ngrids, Bscale)
              DO ng=1,Ngrids
                DO ns=1,Nscale(ng)
                  HdecayX(1,isUbar,ns,ng)=Bscale(ns,ng)
                END DO
              END DO
            CASE ('HdecayIY(isUbar)')
              Npts=load_r(Nval, Rval, Mscale, Ngrids, Bscale)
              DO ng=1,Ngrids
                DO ns=1,Nscale(ng)
                  HdecayY(1,isUbar,ns,ng)=Bscale(ns,ng)
                END DO
              END DO
            CASE ('HdecayIX(isVbar)')
              Npts=load_r(Nval, Rval, Mscale, Ngrids, Bscale)
              DO ng=1,Ngrids
                DO ns=1,Nscale(ng)
                  HdecayX(1,isVbar,ns,ng)=Bscale(ns,ng)
                END DO
              END DO
            CASE ('HdecayIY(isVbar)')
              Npts=load_r(Nval, Rval, Mscale, Ngrids, Bscale)
              DO ng=1,Ngrids
                DO ns=1,Nscale(ng)
                  HdecayY(1,isVbar,ns,ng)=Bscale(ns,ng)
                END DO
              END DO
#  else
            CASE ('HdecayI(isFsur)')
              Npts=load_r(Nval, Rval, Ngrids, Hdecay(1,isFsur,:))
            CASE ('HdecayI(isUbar)')
              Npts=load_r(Nval, Rval, Ngrids, Hdecay(1,isUbar,:))
            CASE ('HdecayI(isVbar)')
              Npts=load_r(Nval, Rval, Ngrids, Hdecay(1,isVbar,:))
#  endif
#  ifdef SOLVE3D
#   ifdef MULTI_SCALE_B
            CASE ('HdecayIX(isUvel)')
              Npts=load_r(Nval, Rval, Mscale, Ngrids, Bscale)
              DO ng=1,Ngrids
                DO ns=1,Nscale(ng)
                  HdecayX(1,isUvel,ns,ng)=Bscale(ns,ng)
                END DO
              END DO
            CASE ('HdecayIY(isUvel)')
              Npts=load_r(Nval, Rval, Mscale, Ngrids, Bscale)
              DO ng=1,Ngrids
                DO ns=1,Nscale(ng)
                  HdecayY(1,isUvel,ns,ng)=Bscale(ns,ng)
                END DO
              END DO
            CASE ('HdecayIX(isVvel)')
              Npts=load_r(Nval, Rval, Mscale, Ngrids, Bscale)
              DO ng=1,Ngrids
                DO ns=1,Nscale(ng)
                  HdecayX(1,isVvel,ns,ng)=Bscale(ns,ng)
                END DO
              END DO
            CASE ('HdecayIY(isVvel)')
              Npts=load_r(Nval, Rval, Mscale, Ngrids, Bscale)
              DO ng=1,Ngrids
                DO ns=1,Nscale(ng)
                  HdecayY(1,isVvel,ns,ng)=Bscale(ns,ng)
                END DO
              END DO
            CASE ('HdecayIX(isTvar)')
              Npts=load_r(Nval, Rval, MT, Mscale, Ngrids, BscaleT)
              DO ng=1,Ngrids
                DO ns=1,Nscale(ng)
                  DO itrc=1,NT(ng)
                    i=isTvar(itrc)
                    HdecayX(1,i,ns,ng)=BscaleT(itrc,ns,ng)
                  END DO
                END DO
              END DO
            CASE ('HdecayIY(isTvar)')
              Npts=load_r(Nval, Rval, MT, Mscale, Ngrids, BscaleT)
              DO ng=1,Ngrids
                DO ns=1,Nscale(ng)
                  DO itrc=1,NT(ng)
                    i=isTvar(itrc)
                    HdecayY(1,i,ns,ng)=BscaleT(itrc,ns,ng)
                  END DO
                END DO
              END DO
#   else
            CASE ('HdecayI(isUvel)')
              Npts=load_r(Nval, Rval, Ngrids, Hdecay(1,isUvel,:))
            CASE ('HdecayI(isVvel)')
              Npts=load_r(Nval, Rval, Ngrids, Hdecay(1,isVvel,:))
            CASE ('HdecayI(isTvar)')
              Npts=load_r(Nval, Rval, MT, Ngrids, Rtracer)
              DO ng=1,Ngrids
                DO itrc=1,NT(ng)
                  Hdecay(1,isTvar(itrc),ng)=Rtracer(itrc,ng)
                END DO
              END DO
#   endif
            CASE ('VdecayI(isUvel)')
              Npts=load_r(Nval, Rval, Ngrids, Vdecay(1,isUvel,:))
            CASE ('VdecayI(isVvel)')
              Npts=load_r(Nval, Rval, Ngrids, Vdecay(1,isVvel,:))
            CASE ('VdecayI(isTvar)')
              Npts=load_r(Nval, Rval, MT, Ngrids, Rtracer)
              DO ng=1,Ngrids
                DO itrc=1,NT(ng)
                  Vdecay(1,isTvar(itrc),ng)=Rtracer(itrc,ng)
                END DO
              END DO
#  endif
#  ifdef ADJUST_BOUNDARY
#   ifdef MULTI_SCALE_B
            CASE ('HdecayB(isFsur)')
              Npts=load_r(Nval, Rval, Mscale, 4, Ngrids, BscaleB)
              DO ng=1,Ngrids
                DO ib=1,4
                  DO ns=1,Nscale(ng)
                    HdecayB(isFsur,ns,ib,ng)=BscaleB(ns,ib,ng)
                  END DO
                END DO
              END DO
            CASE ('HdecayB(isUbar)')
              Npts=load_r(Nval, Rval, Mscale, 4, Ngrids, BscaleB)
              DO ng=1,Ngrids
                DO ib=1,4
                  DO ns=1,Nscale(ng)
                    HdecayB(isUbar,ns,ib,ng)=BscaleB(ns,ib,ng)
                  END DO
                END DO
              END DO
            CASE ('HdecayB(isVbar)')
              Npts=load_r(Nval, Rval, Mscale, 4, Ngrids, BscaleB)
              DO ng=1,Ngrids
                DO ib=1,4
                  DO ns=1,Nscale(ng)
                    HdecayB(isVbar,ns,ib,ng)=BscaleB(ns,ib,ng)
                  END DO
                END DO
              END DO
#   else
            CASE ('HdecayB(isFsur)')
              Npts=load_r(Nval, Rval, 4, Ngrids, HdecayB(isFsur,:,:))
            CASE ('HdecayB(isUbar)')
              Npts=load_r(Nval, Rval, 4, Ngrids, HdecayB(isUbar,:,:))
            CASE ('HdecayB(isVbar)')
              Npts=load_r(Nval, Rval, 4, Ngrids, HdecayB(isVbar,:,:))
#   endif
#   ifdef SOLVE3D
#    ifdef MULTI_SCALE_B
            CASE ('HdecayB(isUvel)')
              Npts=load_r(Nval, Rval, Mscale, 4, Ngrids, BscaleB)
              DO ng=1,Ngrids
                DO ib=1,4
                  DO ns=1,Nscale(ng)
                    HdecayB(isUvel,ns,ib,ng)=BscaleB(ns,ib,ng)
                  END DO
                END DO
              END DO
            CASE ('HdecayB(isVvel)')
              Npts=load_r(Nval, Rval, Mscale, 4, Ngrids, BscaleB)
              DO ng=1,Ngrids
                DO ib=1,4
                  DO ns=1,Nscale(ng)
                    HdecayB(isVvel,ns,ib,ng)=BscaleB(ns,ib,ng)
                  END DO
                END DO
              END DO
            CASE ('HdecayB(isTvar)')
              Npts=load_r(Nval, Rval, Mscale, 4, MT, Ngrids, BscaleBT)
              DO ng=1,Ngrids
                DO itrc=1,NT(ng)
                  DO ib=1,4
                    DO ns=1,Nscale(ng)
                      i=isTvar(itrc)
                      HdecayB(i,ns,ib,ng)=BscaleBT(ns,ib,itrc,ng)
                    END DO
                  END DO
                END DO
              END DO
#    else
            CASE ('HdecayB(isUvel)')
              Npts=load_r(Nval, Rval, 4, Ngrids, HdecayB(isUvel,:,:))
            CASE ('HdecayB(isVvel)')
              Npts=load_r(Nval, Rval, 4, Ngrids, HdecayB(isVvel,:,:))
            CASE ('HdecayB(isTvar)')
              Npts=load_r(Nval, Rval, 4, MT, Ngrids, Rboundary)
              DO ng=1,Ngrids
                DO itrc=1,NT(ng)
                  DO ib=1,4
                    HdecayB(isTvar(itrc),ib,ng)=Rboundary(ib,itrc,ng)
                  END DO
                END DO
              END DO
#    endif
            CASE ('VdecayB(isUvel)')
              Npts=load_r(Nval, Rval, 4, Ngrids, VdecayB(isUvel,:,:))
            CASE ('VdecayB(isVvel)')
              Npts=load_r(Nval, Rval, 4, Ngrids, VdecayB(isVvel,:,:))
            CASE ('VdecayB(isTvar)')
              Npts=load_r(Nval, Rval, 4, MT, Ngrids, Rboundary)
              DO ng=1,Ngrids
                DO itrc=1,NT(ng)
                  DO ib=1,4
                    VdecayB(isTvar(itrc),ib,ng)=Rboundary(ib,itrc,ng)
                  END DO
                END DO
              END DO
#   endif
#  endif
#  ifdef ADJUST_WSTRESS
#   ifdef MULTI_SCALE_B
            CASE ('HdecayFX(isUstr)')
              Npts=load_r(Nval, Rval, Mscale, Ngrids, Bscale)
              DO ng=1,Ngrids
                DO ns=1,Nscale(ng)
                  HdecayX(1,isUstr,ns,ng)=Bscale(ns,ng)
                END DO
              END DO
            CASE ('HdecayFY(isUstr)')
              Npts=load_r(Nval, Rval, Mscale, Ngrids, Bscale)
              DO ng=1,Ngrids
                DO ns=1,Nscale(ng)
                  HdecayY(1,isUstr,ns,ng)=Bscale(ns,ng)
                END DO
              END DO
            CASE ('HdecayFX(isVstr)')
              Npts=load_r(Nval, Rval, Mscale, Ngrids, Bscale)
              DO ng=1,Ngrids
                DO ns=1,Nscale(ng)
                  HdecayX(1,isVstr,ns,ng)=Bscale(ns,ng)
                END DO
              END DO
            CASE ('HdecayFY(isVstr)')
              Npts=load_r(Nval, Rval, Mscale, Ngrids, Bscale)
              DO ng=1,Ngrids
                DO ns=1,Nscale(ng)
                  HdecayY(1,isVstr,ns,ng)=Bscale(ns,ng)
                END DO
              END DO
#   else
            CASE ('HdecayF(isUstr)')
              Npts=load_r(Nval, Rval, Ngrids, Hdecay(1,isUstr,:))
            CASE ('HdecayF(isVstr)')
              Npts=load_r(Nval, Rval, Ngrids, Hdecay(1,isVstr,:))
#   endif
#  endif
#  if defined ADJUST_STFLUX && defined SOLVE3D
#   ifdef MULTI_SCALE_B
            CASE ('HdecayFX(isTsur)')
              Npts=load_r(Nval, Rval, MT, Mscale, Ngrids, BscaleT)
              DO ng=1,Ngrids
                DO ns=1,Nscale(ng)
                  DO itrc=1,NT(ng)
                    i=isTsur(itrc)
                    HdecayX(1,i,ns,ng)=BscaleT(itrc,ns,ng)
                  END DO
                END DO
              END DO
            CASE ('HdecayFY(isTsur)')
              Npts=load_r(Nval, Rval, MT, Mscale, Ngrids, BscaleT)
              DO ng=1,Ngrids
                DO ns=1,Nscale(ng)
                  DO itrc=1,NT(ng)
                    i=isTsur(itrc)
                    HdecayY(1,i,ns,ng)=BscaleT(itrc,ns,ng)
                  END DO
                END DO
              END DO
#   else
            CASE ('HdecayF(isTsur)')
              Npts=load_r(Nval, Rval, MT, Ngrids, Rtracer)
              DO ng=1,Ngrids
                DO itrc=1,NT(ng)
                  Hdecay(1,isTsur(itrc),ng)=Rtracer(itrc,ng)
                END DO
              END DO
#   endif
#  endif

#  ifdef STD_MODEL
#   ifndef COMPUTE_MLD
            CASE ('mld_uniform')
              Npts=load_r(Nval, Rval, Ngrids, mld_uniform(ng))
#   endif
            CASE ('Sigma_max(isFsur)')
              Npts=load_r(Nval, Rval, Ngrids, Sigma_max(isFsur,:))
#   ifdef SOLVE3D
            CASE ('Sigma_max(isUvel)')
              Npts=load_r(Nval, Rval, Ngrids, Sigma_max(isUvel,:))
            CASE ('Sigma_ml(isUvel)')
              Npts=load_r(Nval, Rval, Ngrids, Sigma_ml(isUvel,:))
            CASE ('Sigma_do(isUvel)')
              Npts=load_r(Nval, Rval, Ngrids, Sigma_do(isUvel,:))
            CASE ('Sigma_dz(isUvel)')
              Npts=load_r(Nval, Rval, Ngrids, Sigma_dz(isUvel,:))
            CASE ('Sigma_max(isVvel)')
              Npts=load_r(Nval, Rval, Ngrids, Sigma_max(isVvel,:))
            CASE ('Sigma_ml(isVvel)')
              Npts=load_r(Nval, Rval, Ngrids, Sigma_ml(isVvel,:))
            CASE ('Sigma_do(isVvel)')
              Npts=load_r(Nval, Rval, Ngrids, Sigma_do(isVvel,:))
            CASE ('Sigma_dz(isVvel)')
              Npts=load_r(Nval, Rval, Ngrids, Sigma_dz(isVvel,:))
            CASE ('Sigma_max(isTvar)')
              Npts=load_r(Nval, Rval, MT, Ngrids, Rtracer)
              DO ng=1,Ngrids
                DO itrc=1,NT(ng)
                  Sigma_max(isTvar(itrc),ng)=Rtracer(itrc,ng)
                END DO
              END DO
            CASE ('Sigma_ml(isTvar)')
              Npts=load_r(Nval, Rval, MT, Ngrids, Rtracer)
              DO ng=1,Ngrids
                DO itrc=1,NT(ng)
                  Sigma_ml(isTvar(itrc),ng)=Rtracer(itrc,ng)
                END DO
              END DO
            CASE ('Sigma_do(isTvar)')
              Npts=load_r(Nval, Rval, MT, Ngrids, Rtracer)
              DO ng=1,Ngrids
                DO itrc=1,NT(ng)
                  Sigma_do(isTvar(itrc),ng)=Rtracer(itrc,ng)
                END DO
              END DO
            CASE ('Sigma_dz(isTvar)')
              Npts=load_r(Nval, Rval, MT, Ngrids, Rtracer)
              DO ng=1,Ngrids
                DO itrc=1,NT(ng)
                  Sigma_dz(isTvar(itrc),ng)=Rtracer(itrc,ng)
                END DO
              END DO
#   endif
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
              Npts=load_r(Nval, Rval, MT, Ngrids, Rtracer)
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
              Mval=MAXVAL(Nprovenance)
              Npts=load_i(Nval, Rval, Mval, Ngrids, Iprovenance)
            CASE ('P_bgqc')
              Mval=MAXVAL(Nprovenance)
              Npts=load_r(Nval, Rval, Mval, Ngrids, P_bgqc)
#  endif
#  if defined ADJUST_STFLUX && defined SOLVE3D
            CASE ('Lstflux')
              Npts=load_l(Nval, Cval, MT, Ngrids, Ltsur)
              DO ng=1,Ngrids
                DO itrc=1,NT(ng)
                  Lstflux(itrc,ng)=Ltsur(itrc,ng)
                END DO
              END DO
#  endif
#  ifdef ADJUST_BOUNDARY
            CASE ('Lobc(isFsur)')
              Npts=load_l(Nval, Cval, 4, Ngrids, Lbry)
              DO ng=1,Ngrids
                DO ib=1,4
                  Lobc(ib,isFsur,ng)=Lbry(ib,ng)
                END DO
              END DO
            CASE ('Lobc(isUbar)')
              Npts=load_l(Nval, Cval, 4, Ngrids, Lbry)
              DO ng=1,Ngrids
                DO ib=1,4
                  Lobc(ib,isUbar,ng)=Lbry(ib,ng)
                END DO
              END DO
            CASE ('Lobc(isVbar)')
              Npts=load_l(Nval, Cval, 4, Ngrids, Lbry)
              DO ng=1,Ngrids
                DO ib=1,4
                  Lobc(ib,isVbar,ng)=Lbry(ib,ng)
                END DO
              END DO
#   ifdef SOLVE3D
            CASE ('Lobc(isUvel)')
              Npts=load_l(Nval, Cval, 4, Ngrids, Lbry)
              DO ng=1,Ngrids
                DO ib=1,4
                  Lobc(ib,isUvel,ng)=Lbry(ib,ng)
                END DO
              END DO
            CASE ('Lobc(isVvel)')
              Npts=load_l(Nval, Cval, 4, Ngrids, Lbry)
              DO ng=1,Ngrids
                DO ib=1,4
                  Lobc(ib,isVvel,ng)=Lbry(ib,ng)
                END DO
              END DO
            CASE ('Lobc(isTvar)')
              Npts=load_l(Nval, Cval, 4, MT, Ngrids, Lbry_trc)
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
     &                      Ngrids, Nfiles, 5, 1, inp_lib, STD)
            CASE ('STDnameM')
              label='STD - model error standard deviation'
              Npts=load_s1d(Nval, Cval, Cdim, line, label, igrid,       &
     &                      Ngrids, Nfiles, 5, 2, inp_lib, STD)
            CASE ('STDnameB')
              label='STD - boundary conditions standard deviation'
              Npts=load_s1d(Nval, Cval, Cdim, line, label, igrid,       &
     &                      Ngrids, Nfiles, 5, 3, inp_lib, STD)
            CASE ('STDnameF')
              label='STD - surface forcing standard deviation'
              Npts=load_s1d(Nval, Cval, Cdim, line, label, igrid,       &
     &                      Ngrids, Nfiles, 5, 4, inp_lib, STD)
#  ifdef STD_MODEL
            CASE ('STDnameC')
              label='STD - standard deviation computed from background'
              Npts=load_s1d(Nval, Cval, Cdim, line, label, igrid,       &
     &                      Ngrids, Nfiles, 5, 5, inp_lib, STD)
#  endif
            CASE ('NRMnameI')
              label='NRM - initial conditions normalization'
              Npts=load_s1d(Nval, Cval, Cdim, line, label, igrid,       &
     &                      Ngrids, Nfiles, 4, 1, inp_lib, NRM)
            CASE ('NRMnameM')
              label='NRM - model error normalization'
              Npts=load_s1d(Nval, Cval, Cdim, line, label, igrid,       &
     &                      Ngrids, Nfiles, 4, 2, inp_lib, NRM)
            CASE ('NRMnameB')
              label='NRM - boundary conditions normalization'
              Npts=load_s1d(Nval, Cval, Cdim, line, label, igrid,       &
     &                      Ngrids, Nfiles, 4, 3, inp_lib, NRM)
            CASE ('NRMnameF')
              label='NRM - surface forcing normalization'
              Npts=load_s1d(Nval, Cval, Cdim, line, label, igrid,       &
     &                      Ngrids, Nfiles, 4, 4, inp_lib, NRM)
            CASE ('OBSname')
              label='OBS - data assimilation observations'
              Npts=load_s1d(Nval, Cval, Cdim, line, label, igrid,       &
     &                      Ngrids, Nfiles, obs_lib, OBS)
            CASE ('HSSname')
              label='HSS - Hessian eigenvectors'
              Npts=load_s1d(Nval, Cval, Cdim, line, label, igrid,       &
     &                      Ngrids, Nfiles, out_lib, HSS)
            CASE ('LCZname')
              label='LCZ - Lanczos vectors'
              Npts=load_s1d(Nval, Cval, Cdim, line, label, igrid,       &
     &                      Ngrids, Nfiles, out_lib, LCZ)
            CASE ('LZEname')
              label='LZE - Time-evolved Lanczos vectors'
              Npts=load_s1d(Nval, Cval, Cdim, line, label, igrid,       &
     &                      Ngrids, Nfiles, out_lib, LZE)
            CASE ('MODname')
              label='DAV - 4D-Var data assimilation variables'
              Npts=load_s1d(Nval, Cval, Cdim, line, label, igrid,       &
     &                      Ngrids, Nfiles, obs_lib, DAV)
            CASE ('ERRname')
              label='ERR - 4D-Var posterior error covariance'
              Npts=load_s1d(Nval, Cval, Cdim, line, label, igrid,       &
     &                      Ngrids, Nfiles, out_lib, ERR)
#  ifdef SP4DVAR
            CASE ('SPTname')
              label='SPT - TLM Arnoldi vectors'
              Npts=load_s1d(Nval, Cval, Cdim, line, label, igrid,       &
     &                      Ngrids, Nfiles, out_lib, SPT)
            CASE ('SPAname')
              label='SPA - ADM Arnoldi vectors'
              Npts=load_s1d(Nval, Cval, Cdim, line, label, igrid,       &
     &                      Ngrids, Nfiles, out_lib, SPA)
            CASE ('SCTname')
              label='SCT - TLM Arnoldi vectors'
              Npts=load_s1d(Nval, Cval, Cdim, line, label, igrid,       &
     &                      Ngrids, Nfiles, out_lib, SCT)
            CASE ('SCAname')
              label='SCA - ADM Arnoldi vectors'
              Npts=load_s1d(Nval, Cval, Cdim, line, label, igrid,       &
     &                      Ngrids, Nfiles, out_lib, SCA)
#  endif
#  if defined MULTI_SCALE_B && defined NONUNIFORM_SCALES
            CASE ('SVCname')
              label='SVC - Spatially-varying correlation length scales'
              Npts=load_s1d(Nval, Cval, Cdim, line, label, igrid,       &
     &                      Ngrids, Nfiles, inp_lib, SVC)
#  endif
#  if defined RBL4DVAR_FCT_SENSITIVITY && defined OBS_SPACE
            CASE ('OIFnameA')
              label='OIFA - observation impacts forecast, analysis'
              Npts=load_s1d(Nval, Cval, Cdim, line, label, igrid,       &
     &                      Ngrids, Nfiles, inp_lib, OIFA)
            CASE ('OIFnameB')
              label='OIFB - observation impacts forecast, background'
              Npts=load_s1d(Nval, Cval, Cdim, line, label, igrid,       &
     &                      Ngrids, Nfiles, inp_lib, OIFB)
#  endif
            CASE ('IMPname')
              label='IMP - Coarser inners loop impulse forcing'
              Npts=load_s1d(Nval, Cval, Cdim, line, label, igrid,       &
     &                      Ngrids, Nfiles, out_lib, IMP)
            CASE ('INCname')
              label='INC - Coarser inners loops increments'
              Npts=load_s1d(Nval, Cval, Cdim, line, label, igrid,       &
     &                      Ngrids, Nfiles, out_lib, INC)
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
        IF (Master) THEN
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
#    ifdef SALINITY
      IF (.not.balance(isTvar(isalt)).and.balance(isFsur)) THEN
        balance(isTvar(isalt))=.TRUE.
      END IF
      IF (.not.balance(isTvar(isalt)).and.balance(isVvel)) THEN
        balance(isTvar(isalt))=.TRUE.
      END IF
      IF (balance(isTvar(isalt))) THEN
        balance(isTvar(itemp))=.TRUE.
      END IF
#    endif
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
      REPORT_LOOP : DO ng=1,Ngrids
        IF (Master.and.Lwrite) THEN
#  if defined FOUR_DVAR || defined VERIFICATION
          WRITE (out,60) ng
#  endif
#  if defined FOUR_DVAR
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
          WRITE (out,100) timeIAU(ng), 'timeIAU',                       &
     &            'Duration of the incremental analysis update (days).'
#   ifndef I4DVAR_ANA_SENSITIVITY
#    ifdef I4DVAR
          WRITE (out,100) GradErr, 'GradErr',                           &
     &            'Upper bound on relative error of the gradient.'
          WRITE (out,70) LhessianEV, 'LhessianEV',                      &
     &            'Switch to compute Hessian eigenvectors.'
#    endif
          WRITE (out,100) HevecErr, 'HevecErr',                         &
     &            'Accuracy required for eigenvectors.'
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
#     ifdef SALINITY
          WRITE (out,70) balance(isTvar(isalt)), 'balance(isSalt)',     &
                  'Switch to include salinity in balance operator.'
#     endif
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
#   if defined RBL4DVAR_ANA_SENSITIVITY || \
       defined RBL4DVAR_FCT_SENSITIVITY || \
       defined R4DVAR_ANA_SENSITIVITY
          WRITE (out,80) Nimpact, 'Nimpact',                            &
     &            'Observations impact/sensitivity outer loop to use.'
          IF (Nimpact.gt.Nouter) THEN
            IF (Master) THEN
              WRITE (out,240) 'Nimpact', Nimpact,                       &
     &                        'must be less or equal than Nouter',      &
     &                        TRIM(aparnam)
            END IF
            exit_flag=5
            RETURN
          END IF
#   endif
#   if defined ARRAY_MODES
          WRITE (out,80) Nimpact, 'Nimpact',                            &
     &            'Array mode outer loop.'
          IF (Nimpact.gt.Nouter) THEN
            IF (Master) THEN
              WRITE (out,240) 'Nimpact', Nimpact,                       &
     &                        'must be less or equal than Nouter',      &
     &                        TRIM(aparnam)
            END IF
            exit_flag=5
            RETURN
          END IF
#   endif
#   if defined SPLIT_4DVAR
          WRITE (out,80) OuterLoop, 'OuterLoop',                        &
     &            'Current outer loop counter.'
          IF (OuterLoop.gt.Nouter) THEN
            IF (Master) THEN
              WRITE (out,240) 'OuterLoop', OuterLoop,                   &
     &                        'must be less or equal than Nouter',      &
     &                        TRIM(aparnam)
            END IF
            exit_flag=5
            RETURN
          END IF
          WRITE (out,95) Phase4DVAR(1:10), 'Phase4DVAR',                &
     &            'Current 4D-Var phase (first 10 characters).'
#   endif
#   ifndef TLM_CHECK
#    ifndef I4DVAR_ANA_SENSITIVITY
          WRITE (out,170) (LdefNRM(i,ng),i=1,4), 'LdefNRM',             &
     &            'Switch to create a normalization NetCDF file.'
          WRITE (out,170) (LwrtNRM(i,ng),i=1,4), 'LwrtNRM',             &
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
#     if defined RBL4DVAR           || defined R4DVAR    || \
         defined SENSITIVITY_4DVAR  || defined SP4DVAR   || \
         defined TL_RBL4DVAR        || defined TL_R4DVAR
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
            WRITE (out,170) (CnormB(isFsur,i),i=1,4), 'CnormB(isFsur)', &
     &            'Compute boundary 2D RHO-normalization factors.'
            WRITE (out,170) (CnormB(isUbar,i),i=1,4), 'CnormB(isUbar)', &
     &            'Compute boundary 2D U-normalization factors.'
            WRITE (out,170) (CnormB(isVbar,i),i=1,4), 'CnormB(isVbar)', &
     &            'Compute boundary 2D V-normalization factors.'
#      ifdef SOLVE3D
            WRITE (out,170) (CnormB(isUvel,i),i=1,4), 'CnormB(isUvel)', &
     &            'Compute boundary 3D U-normalization factors.'
            WRITE (out,170) (CnormB(isVvel,i),i=1,4), 'CnormB(isVvel)', &
     &            'Compute boundary 3D V-normalization factors.'
            DO itrc=1,NT(ng)
              WRITE (out,180) (CnormB(isTvar(itrc),i),i=1,4),           &
     &                        'CnormB(isTvar)',                         &
     &            'Compute boundary normalization factors for tracer ', &
     &            itrc, TRIM(Vname(1,idTvar(itrc)))
            END DO
#      endif
          END IF
#     endif
#     ifdef ADJUST_WSTRESS
          IF (ANY(LwrtNRM(:,ng))) THEN
            WRITE (out,70) Cnorm(1,isUstr), 'CnormF(isUstr)',           &
     &            'Compute surface U-stress normalization factors.'
            WRITE (out,70) Cnorm(1,isVstr), 'CnormF(isVstr)',           &
     &            'Compute surface V-stress normalization factors.'
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
#    if defined NORMALIZATION || defined RBL4DVAR           || \
        defined R4DVAR        || defined SENSITIVITY_4DVAR  || \
        defined SP4DVAR       || defined TL_RBL4DVAR        || \
        defined TL_R4DVAR
#     ifdef MULTI_SCALE_B
          WRITE (out,80) Nscale(ng), 'Nscale',                          &
     &          'Number of multiscales for modeling B-Error Covariance.'
          DO ns=1,Nscale(ng)
            WRITE (out,90) NiterCG(ns,ng), 'NiterCG .........',         &
     &            'B-scale ', ns, ', Number of multiscale Conjugate '// &
     &            'Gradient solver iterations.'
            WRITE (out,90) NiterCI(ns,ng), 'NiterCI',                   &
     &            'B-scale ', ns, ', Number of multiscale implicit '//  &
     &            'Chebyshev solver iterations.'
            WRITE (out,90) Nlap(isFsur,ns,ng), 'Mlap(isFsur)',          &
     &            'B-scale ', ns, ', smoother applications: '//         &
     &            'free-surface.'
            IF ((Nlap(isFsur,ns,ng).lt.2).or.                           &
     &          (MOD(Nlap(isFsur,ns,ng),2).ne.0)) THEN
              WRITE (out,240) 'Mlap(isFsur)', Nlap(isFsur,ns,ng),       &
     &                        'Must be greater than 2 and even!',       &
     &                        TRIM(aparnam)
              exit_flag=5
            END IF
            WRITE (out,90) Nlap(isUbar,ns,ng), 'Mlap(isUbar)',          &
     &            'B-scale ', ns, ', smoother applications: '//         &
     &            '2D U-momentum.'
            IF ((Nlap(isUbar,ns,ng).lt.2).or.                           &
     &          (MOD(Nlap(isUbar,ns,ng),2).ne.0)) THEN
              WRITE (out,240) 'Mlap(isUbar)', Nlap(isUbar,ns,ng),       &
     &                        'Must be greater than 2 and even!',       &
     &                        TRIM(aparnam)
              exit_flag=5
            END IF
            WRITE (out,90) Nlap(isVbar,ns,ng), 'Mlap(isVbar)',          &
     &            'B-scale ', ns, ', smoother applications: '//         &
     &            '2D V-momentum.'
            IF ((Nlap(isVbar,ns,ng).lt.2).or.                           &
     &          (MOD(Nlap(isVbar,ns,ng),2).ne.0)) THEN
              WRITE (out,240) 'Mlap(isVbar)', Nlap(isVbar,ns,ng),       &
     &                        'Must be greater than 2 and even!',       &
     &                        TRIM(aparnam)
              exit_flag=5
            END IF
#      ifdef SOLVE3D
            WRITE (out,90) Nlap(isUvel,ns,ng), 'Mlap(isUvel)',          &
     &            'B-scale ', ns, ', smoother applications: '//         &
     &            '3D U-momentum.'
            IF ((Nlap(isUvel,ns,ng).lt.2).or.                           &
     &          (MOD(Nlap(isUvel,ns,ng),2).ne.0)) THEN
              WRITE (out,240) 'Mlap(isUvel)', Nlap(isUvel,ns,ng),       &
     &                        'Must be greater than 2 and even!',       &
     &                        TRIM(aparnam)
              exit_flag=5
            END IF
            WRITE (out,90) Nlap(isVvel,ns,ng), 'Mlap(isVvel)',          &
     &            'B-scale ', ns, ', smoother applications: '//         &
     &            '3D V-momentum.'
            IF ((Nlap(isVvel,ns,ng).lt.2).or.                           &
     &          (MOD(Nlap(isVvel,ns,ng),2).ne.0)) THEN
              WRITE (out,240) 'Mlap(isVvel)', Nlap(isVvel,ns,ng),       &
     &                        'Must be greater than 2 and even!',       &
     &                        TRIM(aparnam)
              exit_flag=5
            END IF
            DO itrc=1,NT(ng)
              WRITE (out,90) Nlap(isTvar(itrc),ns,ng), 'Mlap(idTvar)',  &
     &            'B-scale ', ns, ', smoother applications: '//         &
     &            TRIM(Vname(1,idTvar(itrc)))
              IF ((Nlap(isTvar(itrc),ns,ng).lt.2).or.                   &
     &            (MOD(Nlap(isTvar(itrc),ns,ng),2).ne.0)) THEN
                WRITE (out,240) 'Mlap(isTvar)',                         &
     &                          Nlap(isTvar(itrc),ns,ng),               &
     &                          'Must be greater than 2 and even!',     &
     &                          TRIM(aparnam)
                exit_flag=5
              END IF
            END DO
#      endif
#      ifdef ADJUST_WSTRESS
            WRITE (out,90) Nlap(isUstr,ns,ng), 'Mlap(isUstr)',          &
     &            'B-scale ', ns, ', smoother applications: '//         &
     &            'Surface U-stress.'
            IF ((Nlap(isUstr,ns,ng).lt.2).or.                           &
     &          (MOD(Nlap(isUstr,ns,ng),2).ne.0)) THEN
              WRITE (out,240) 'Mlap(isUstr)', Nlap(isUstr,ns,ng),       &
     &                        'Must be greater than 2 and even!',       &
     &                        TRIM(aparnam)
              exit_flag=5
            END IF
            WRITE (out,90) Nlap(isVstr,ns,ng), 'Mlap(isVstr)',          &
     &            'B-scale ', ns, ', smoother applications: '//         &
     &            'Surface V-stress.'
            IF ((Nlap(isVstr,ns,ng).lt.2).or.                           &
     &          (MOD(Nlap(isVstr,ns,ng),2).ne.0)) THEN
              WRITE (out,240) 'Mlap(isVstr)', Nlap(isVstr,ns,ng),       &
     &                        'Must be greater than 2 and even!',       &
     &                        TRIM(aparnam)
              exit_flag=5
            END IF
#      endif
#      if defined ADJUST_STFLUX && defined SOLVE3D
            DO itrc=1,NT(ng)
              WRITE (out,90) Nlap(isTsur(itrc),ns,ng), 'Mlap(idTsur)',  &
     &            'B-scale ', ns, ', smoother applications: '//         &
     &            'Surface flux, '//TRIM(Vname(1,idTvar(itrc)))//'.'
              IF ((Nlap(isTsur(itrc),ns,ng).lt.2).or.                   &
     &            (MOD(Nlap(isTsur(itrc),ns,ng),2).ne.0)) THEN
                WRITE (out,240) 'Mlap(isTsur)',                         &
     &                          Nlap(isTsur(itrc),ns,ng),               &
     &                          'Must be greater than 2 and even!',     &
     &                          TRIM(aparnam)
                exit_flag=5
              END IF
            END DO
#      endif
            WRITE (out,125) Bwgt(isFsur,ns,ng), 'Bwgt(isFsur)',         &
     &            'B-scale ', ns, ', multiscale weight: free-surface.'
            WRITE (out,125) Bwgt(isUbar,ns,ng), 'Bwgt(isUbar)',         &
     &            'B-scale ', ns, ', multiscale weight: 2D U-momentum.'
            WRITE (out,125) Bwgt(isVbar,ns,ng), 'Bwgt(isVbar)',         &
     &            'B-scale ', ns, ', multiscale weight: 2D V-momentum.'
#      ifdef SOLVE3D
            WRITE (out,125) Bwgt(isUvel,ns,ng), 'Bwgt(isUvel)',         &
     &            'B-scale ', ns, ', multiscale weight: 3D U-momentum.'
            WRITE (out,125) Bwgt(isVvel,ns,ng), 'Bwgt(isVvel)',         &
     &            'B-scale ', ns, ', multiscale weight: 3D V-momentum.'
            DO itrc=1,NT(ng)
              WRITE (out,125) Bwgt(isTvar(itrc),ns,ng), 'Bwgt(idTvar)', &
     &            'B-scale ', ns, ', multiscale weight: '//             &
     &            TRIM(Vname(1,idTvar(itrc)))//'.'
            END DO
#      endif
#      ifdef ADJUST_WSTRESS
            WRITE (out,125) Bwgt(isUstr,ns,ng), 'Bwgt(isUstr)',         &
     &            'B-scale ', ns, ', multiscale weight: '//             &
     &            'surface U-stress.'
            WRITE (out,125) Bwgt(isVstr,ns,ng), 'Bwgt(isVstr)',         &
     &            'B-scale ', ns, ', multiscale weight: '//             &
     &            'surface V-stress.'
#      endif
#      if defined ADJUST_STFLUX && defined SOLVE3D
            DO itrc=1,NT(ng)
              WRITE (out,125) Bwgt(isTsur(itrc),ns,ng), 'Bwgt(idTvar)', &
     &            'B-scale ', ns, ', multiscale weight, '//             &
     &            'surface flux: '//TRIM(Vname(1,idTvar(itrc)))//'.'
            END DO
#      endif
          END DO
#     endif
          IF (nADJ(ng).lt.ntimes(ng)) THEN
#     ifdef MULTI_SCALE_B
#      ifndef NONUNIFORM_SCALES
            DO ns=1,Nscale(ng)
              WRITE (out,125) HdecayX(2,isFsur,ns,ng),                  &
     &                        'HdecayMX(isFsur)', 'B-scale ', ns,       &
     &            ', Model Error X-decorrelation (m), free-surface.'
              WRITE (out,125) HdecayY(2,isFsur,ns,ng),                  &
     &                        'HdecayMY(isFsur)', 'B-scale ', ns,       &
     &            ', Model Error Y-decorrelation (m), free-surface.'
              WRITE (out,125) HdecayX(2,isUbar,ns,ng),                  &
     &                        'HdecayMX(isUbar)', 'B-scale ', ns,       &
     &            ', Model Error X-decorrelation (m), 2D U-momentum.'
              WRITE (out,125) HdecayY(2,isUbar,ns,ng),                  &
     &                        'HdecayMY(isUbar)', 'B-scale ', ns,       &
     &            ', Model Error Y-decorrelation (m), 2D U-momentum.'
              WRITE (out,125) HdecayX(2,isVbar,ns,ng),                  &
     &                        'HdecayMX(isVbar)', 'B-scale ', ns,       &
     &            ', Model Error X-decorrelation (m), 2D V-momentum.'
              WRITE (out,125) HdecayY(2,isVbar,ns,ng),                  &
     &                        'HdecayMY(isVbar)', 'B-scale ', ns,       &
     &            ', Model Error Y-decorrelation (m), 2D V-momentum.'
#       ifdef SOLVE3D
              WRITE (out,125) HdecayX(2,isUvel,ns,ng),                  &
     &                        'HdecayMX(isUvel)', 'B-scale ', ns,       &
     &            ', Model Error X-decorrelation (m), 3D U-momentum.'
              WRITE (out,125) HdecayY(2,isUvel,ns,ng),                  &
     &                        'HdecayMY(isUvel)', 'B-scale ', ns,       &
     &            ', Model Error Y-decorrelation (m), 3D U-momentum.'
              WRITE (out,125) HdecayX(2,isVvel,ns,ng),                  &
     &                        'HdecayMX(isVvel)', 'B-scale ', ns,       &
     &            'Model Error X-decorrelation (m), 3D V-momentum.'
              WRITE (out,125) HdecayY(2,isVvel,ns,ng),                  &
     &                        'HdecayMY(isVvel)', 'B-scale ', ns,       &
     &            'Model Error Y-decorrelation (m), 3D V-momentum.'
              DO itrc=1,NT(ng)
                WRITE (out,125) HdecayX(2,isTvar(itrc),ns,ng),          &
     &                          'HdecayMX(idTvar)', 'B-scale ', ns,     &
     &            'Model Error X-decorrelation (m), '//                 &
     &            TRIM(Vname(1,idTvar(itrc)))//'.'
                WRITE (out,125) HdecayY(2,isTvar(itrc),ns,ng),          &
     &                          'HdecayMY(idTvar)', 'Bscale ', ns,      &
     &            'Model Error Y-decorrelation (m), '//                 &
     &            TRIM(Vname(1,idTvar(itrc)))//'.'
              END DO
#       endif
            END DO
#      endif
#     else
            WRITE (out,120) Hdecay(2,isFsur,ng), 'HdecayM(isFsur)',     &
     &        'Model Error decorrelation H-scale (m), free-surface.'
            WRITE (out,120) Hdecay(2,isUbar,ng), 'HdecayM(isUbar)',     &
     &        'Model Error decorrelation H-scale (m), 2D U-momentum.'
            WRITE (out,120) Hdecay(2,isVbar,ng), 'HdecayM(isVbar)',     &
     &        'Model Error decorrelation H-scale (m), 2D V-momentum.'
#      ifdef SOLVE3D
            WRITE (out,120) Hdecay(2,isUvel,ng), 'HdecayM(isUvel)',     &
     &        'Model Error decorrelation H-scale (m), 3D U-momentum.'
            WRITE (out,120) Hdecay(2,isVvel,ng), 'HdecayM(isVvel)',     &
     &        'Model Error decorrelation H-scale (m), 3D V-momentum.'
            DO itrc=1,NT(ng)
              WRITE (out,130) Hdecay(2,isTvar(itrc),ng),                &
     &            'HdecayM(idTvar)',                                    &
     &            'Model Error decorrelation H-scale (m), ',            &
     &            TRIM(Vname(1,idTvar(itrc)))
            END DO
#      endif
#     endif
#     ifdef SOLVE3D
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
#    ifdef MULTI_SCALE_B
#     ifndef NONUNIFORM_SCALES
          DO ns=1,Nscale(ng)
            WRITE (out,125) HdecayX(1,isFsur,ns,ng),                    &
     &                      'HdecayIX(isFsur)', 'B-scale ', ns,         &
     &            ', IC Error X-decorrelation (m), free-surface.'
            WRITE (out,125) HdecayY(1,isFsur,ns,ng),                    &
     &                        'HdecayIY(isFsur)', 'B-scale ', ns,       &
     &            ', IC Error Y-decorrelation (m), free-surface.'
            WRITE (out,125) HdecayX(1,isUbar,ns,ng),                    &
     &                      'HdecayIX(isUbar)', 'B-scale ', ns,         &
     &            ', IC Error X-decorrelation (m), 2D U-momentum.'
            WRITE (out,125) HdecayY(1,isUbar,ns,ng),                    &
     &                      'HdecayIY(isUbar)', 'B-scale ', ns,         &
     &            ', IC Error Y-decorrelation (m), 2D U-momentum.'
            WRITE (out,125) HdecayX(1,isVbar,ns,ng),                    &
     &                      'HdecayIX(isVbar)', 'B-scale ', ns,         &
     &            ', IC Error X-decorrelation (m), 2D V-momentum.'
            WRITE (out,125) HdecayY(1,isVbar,ns,ng),                    &
     &                      'HdecayIY(isVbar)', 'B-scale ', ns,         &
     &            ', IC Error Y-decorrelation (m), 2D V-momentum.'
#      ifdef SOLVE3D
            WRITE (out,125) HdecayX(1,isUvel,ns,ng),                    &
     &                      'HdecayIX(isUvel)', 'B-scale ', ns,         &
     &            ', IC Error X-decorrelation (m), 3D U-momentum.'
            WRITE (out,125) HdecayY(1,isUvel,ns,ng),                    &
     &                      'HdecayIY(isUvel)', 'B-scale ', ns,         &
     &            ', IC Error Y-decorrelation (m), 3D U-momentum.'
            WRITE (out,125) HdecayX(1,isVvel,ns,ng),                    &
     &                      'HdecayIX(isVvel)', 'B-scale ', ns,         &
     &            ', IC Error X-decorrelation (m), 3D V-momentum.'
            WRITE (out,125) HdecayY(1,isVvel,ns,ng),                    &
     &                      'HdecayIY(isVvel)', 'B-scale ', ns,         &
     &            ', IC Error Y-decorrelation (m), 3D V-momentum.'
            DO itrc=1,NT(ng)
              WRITE (out,125) HdecayX(1,isTvar(itrc),ns,ng),            &
     &                        'HdecayIX(idTvar)', 'B-scale ', ns,       &
     &              ', IC Error X-decorrelation (m), '//                &
     &              TRIM(Vname(1,idTvar(itrc)))//'.'
              WRITE (out,125) HdecayY(1,isTvar(itrc),ns,ng),            &
     &                        'HdecayIY(idTvar)', 'B-scale ', ns,       &
     &              ', IC Error Y-decorrelation (m), '//                &
     &              TRIM(Vname(1,idTvar(itrc)))//'.'
            END DO
#      endif
          END DO
#     endif
#    else
          WRITE (out,120) Hdecay(1,isFsur,ng), 'HdecayI(isFsur)',       &
     &            'IC Error decorrelation H-scale (m), free-surface.'
          WRITE (out,120) Hdecay(1,isUbar,ng), 'HdecayI(isUbar)',       &
     &            'IC Error decorrelation H-scale (m), 2D U-momentum.'
          WRITE (out,120) Hdecay(1,isVbar,ng), 'HdecayI(isVbar)',       &
     &            'IC Error decorrelation H-scale (m), 2D V-momentum.'
#     ifdef SOLVE3D
          WRITE (out,120) Hdecay(1,isUvel,ng), 'HdecayI(isUvel)',       &
     &            'IC Error decorrelation H-scale (m), 3D U-momentum.'
          WRITE (out,120) Hdecay(1,isVvel,ng), 'HdecayI(isVvel)',       &
     &            'IC Error decorrelation H-scale (m), 3D V-momentum.'
          DO itrc=1,NT(ng)
            WRITE (out,130) Hdecay(1,isTvar(itrc),ng),                  &
     &            'HdecayI(idTvar)',                                    &
     &            'IC Error decorrelation H-scale (m), ',               &
     &            TRIM(Vname(1,idTvar(itrc)))
          END DO
#     endif
#    endif
#    ifdef SOLVE3D      
          WRITE (out,120) Vdecay(1,isUvel,ng), 'VdecayI(isUvel)',       &
     &            'IC Error decorrelation V-scale (m), 3D U-momentum.'
          WRITE (out,120) Vdecay(1,isVvel,ng), 'VdecayI(isVvel)',       &
     &            'IC Error decorrelation V-scale (m), 3D V-momentum.'
          DO itrc=1,NT(ng)
            WRITE (out,130) Vdecay(1,isTvar(itrc),ng),                  &
     &            'VdecayI(idTvar)',                                    &
     &            'IC Error decorrelation V-scale (m), ',               &
     &            TRIM(Vname(1,idTvar(itrc)))
          END DO
#    endif
#    ifdef ADJUST_BOUNDARY
#     ifdef MULTI_SCALE_B
          DO ib=1,4
            DO ns=1,Nscale(ng)
              IF (ib.eq.iwest) THEN
                Text=', W-bry'
              ELSE IF (ib.eq.isouth) THEN
                Text=', S-bry'
              ELSE IF (ib.eq.ieast) THEN
                Text=', E-bry'
              ELSE IF (ib.eq.inorth) THEN
                Text=', N-bry'
              END IF
              IF (Lobc(ib,isFsur,ng)) THEN
                WRITE (out,125) HdecayB(isFsur,ns,ib,ng),               &
     &                          'HdecayB(isFsur)', 'B-scale ', ns,      &
     &              Text//' decorrelation H-scale (m), free-surface.'
              END IF
              IF (Lobc(ib,isUbar,ng)) THEN
                WRITE (out,125) HdecayB(isUbar,ns,ib,ng),               &
     &                          'HdecayB(isUbar)', 'B-scale ', ns,      &
     &              Text//' decorrelation H-scale (m), 2D U-momentum.'
              END IF
              IF (Lobc(ib,isVbar,ng)) THEN
                WRITE (out,125) HdecayB(isVbar,ns,ib,ng),               &
     &                          'HdecayB(isVbar)', 'B-scale ', ns,      &
     &              Text//' decorrelation H-scale (m), 2D V-momentum.'
              END IF
#      ifdef SOLVE3D
              IF (Lobc(ib,isUvel,ng)) THEN
                WRITE (out,125) HdecayB(isUvel,ns,ib,ng),               &
     &                          'HdecayB(isUvel)', 'B-scale ', ns,      &
     &              Text//' decorrelation H-scale (m), 3D U-momentum.'
              END IF
              IF (Lobc(ib,isVvel,ng)) THEN
                WRITE (out,125) HdecayB(isVvel,ns,ib,ng),               &
     &                          'HdecayB(isVvel)', 'B-scale ', ns,      &
     &              Text//' decorrelation H-scale (m), 2D V-momentum.'
              END IF
              DO i=1,NT(ng)
                IF (Lobc(ib,isTvar(i),ng)) THEN
                  WRITE(out,125) HdecayB(isTvar(i),ns, ib,ng),          &
     &                           'HdecayB(idTvar)', 'B-scale ', ns,     &
     &              Text//' decorrelation H-scale (m), '//              &
     &              TRIM(Vname(1,idTvar(i)))//'.'
                END IF
              END DO
#      endif
            END DO
          END DO
#     else
          DO ib=1,4
            IF (ib.eq.iwest) THEN
              Text='W-bry'
            ELSE IF (ib.eq.isouth) THEN
              Text='S-bry'
            ELSE IF (ib.eq.ieast) THEN
              Text='E-bry'
            ELSE IF (ib.eq.inorth) THEN
              Text='N-bry'
            END IF
            IF (Lobc(ib,isFsur,ng)) THEN
              WRITE (out,120) HdecayB(isFsur,ib,ng),                    &
     &                        'HdecayB(isFsur)', TRIM(Text)//           &
     &              ', decorrelation H-scale (m), free-surface.'
            END IF
            IF (Lobc(ib,isUbar,ng)) THEN
              WRITE (out,120) HdecayB(isUbar,ib,ng),                    &
     &                        'HdecayB(isUbar)', TRIM(Text)//           &
     &              ', decorrelation H-scale (m), 2D U-momentum.'
            END IF
            IF (Lobc(ib,isVbar,ng)) THEN
              WRITE (out,120) HdecayB(isVbar,ib,ng),                    &
     &                        'HdecayB(isVbar)', TRIM(Text)//           &
     &              ', decorrelation H-scale (m), 2D V-momentum.'
            END IF
#      ifdef SOLVE3D
            IF (Lobc(ib,isUvel,ng)) THEN
              WRITE (out,120) HdecayB(isUvel,ib,ng),                    &
     &                        'HdecayB(isUvel)', TRIM(Text)//           &
     &              ', decorrelation H-scale (m), 3D U-momentum.'
            END IF
            IF (Lobc(ib,isVvel,ng)) THEN
              WRITE (out,120) HdecayB(isVvel,ib,ng),                    &
     &                        'HdecayB(isVvel)', TRIM(Text)//           &
     &              ', decorrelation H-scale (m), 3D V-momentum.'
            END IF
            DO i=1,NT(ng)
              IF (Lobc(ib,isTvar(i),ng)) THEN
                WRITE(out,130) HdecayB(isTvar(i),ib,ng),                &
     &                         'HdecayB(idTvar)', TRIM(Text)//          &
     &               ', decorrelation H-scale (m), ',                   &
     &              TRIM(Vname(1,idTvar(i)))//'.'
              END IF
            END DO
#      endif
          END DO
#     endif
#     ifdef SOLVE3D
          DO ib=1,4
            IF (ib.eq.iwest) THEN
              Text='W-bry'
            ELSE IF (ib.eq.isouth) THEN
              Text='S-bry'
            ELSE IF (ib.eq.ieast) THEN
              Text='E-bry'
            ELSE IF (ib.eq.inorth) THEN
              Text='N-bry'
            END IF
            IF (Lobc(ib,isUvel,ng)) THEN
              WRITE (out,120) VdecayB(isUvel,ib,ng),                    &
     &                        'VdecayB(isUvel)', TRIM(Text)//           &
     &              ', decorrelation V-scale (m), 3D U-momentum.'
            END IF
            IF (Lobc(ib,isVvel,ng)) THEN
              WRITE (out,120) VdecayB(isVvel,ib,ng),                    &
     &                        'VdecayB(isVvel)', TRIM(Text)//           &
     &              ', decorrelation V-scale (m), 3D V-momentum.'
            END IF
            DO i=1,NT(ng)
              IF (Lobc(ib,isTvar(i),ng)) THEN
                WRITE(out,130) VdecayB(isTvar(i),ib,ng),                &
     &                         'VdecayB(idTvar)', TRIM(Text)//          &
     &               ', decorrelation V-scale (m), ',                   &
     &               TRIM(Vname(1,idTvar(i)))
              END IF
            END DO
#     endif
          END DO
#    endif
#    ifdef ADJUST_WSTRESS
#     ifdef MULTI_SCALE_B
#      ifndef NONUNIFORM_SCALES
          DO ns=1,Nscale(ng)
            WRITE (out,125) HdecayX(1,isUstr,ns,ng),                    &
     &                      'HdecayFX(isUstr)', 'B-scale ', ns,         &
     &            ', Forcing Error X-decorrelation (m), U-stress.'
            WRITE (out,125) HdecayY(1,isVstr,ns,ng),                    &
     &                      'HdecayFY(isVstr)', 'B-scale ', ns,         &
     &            ', Forcing Error Y-decorrelation (m), V-stress.'
          END DO
#      endif
#     else
          WRITE (out,120) Hdecay(1,isUstr,ng), 'HdecayF(isUstr)',       &
     &            'Forcing Error decorrelation H-scale (m), U-stress.'
          WRITE (out,120) Hdecay(1,isVstr,ng), 'HdecayF(isVstr)',       &
     &            'Forcing Error decorrelation H-scale (m), V-stress.'
#     endif
#    endif
#    if defined ADJUST_STFLUX && defined SOLVE3D
#     ifdef MULTI_SCALE_B
#      ifndef NONUNIFORM_SCALES
          DO ns=1,Nscale(ng)
            DO itrc=1,NT(ng)
              WRITE (out,125) HdecayX(1,isTsur(itrc),ns,ng),            &
     &                        'HdecayFX(isTsur)', 'B-scale ', ns,       &
     &              ', Forcing Error X-decorrelation (m), '//           &
     &              TRIM(Vname(1,idTvar(itrc)))//'.'
              WRITE (out,125) HdecayY(1,isTsur(itrc),ns,ng),            &
     &                        'HdecayFY(isTsur)', 'B-scale ', ns,       &
     &              ', Forcing Error Y-decorrelation (m), '//           &
     &              TRIM(Vname(1,idTvar(itrc)))//'.'
            END DO
          END DO
#      endif
#     else
          DO itrc=1,NT(ng)
            WRITE (out,130) Hdecay(1,isTsur(itrc),ng),                  &
     &            'HdecayF(idTsur)',                                    &
     &            'Forcing Error decorrelation H-scale (m), ',          &
     &            TRIM(Vname(1,idTvar(itrc)))//'.'
          END DO
#     endif
#    endif
#    ifdef STD_MODEL
#     ifndef COMPUTE_MLD
          WRITE (out,120) mld_uniform(ng), 'mld_uniform',               &
     &            'Uniform mixed layer depth value (m).'
#     endif
          WRITE (out,190) Sigma_max(isFsur,ng), 'Sigma_max(isFsur)',    &
     &            'Maximum STD value (m)', ', free-surface.'
#     ifdef SOLVE3D
          WRITE (out,190) Sigma_max(isUvel,ng), 'Sigma_max(isUvel)',    &
     &            'Maximum STD value (m/s)', ', 3D U-momentum.'
          WRITE (out,190) Sigma_ml(isUvel,ng), 'Sigma_ml(isUvel)',      &
     &            'Minimum STD at mixed layer (m/s)', ', 3D U-momentum.'
          WRITE (out,190) Sigma_do(isUvel,ng), 'Sigma_do(isUvel)',      &
     &            'Minimum STD in deep ocean (m/s)', ', 3D U-momentum.'
          WRITE (out,130) Sigma_dz(isUvel,ng), 'Sigma_dz(isUvel)',      &
     &            'STD vertical displacement (m)', ', 3D U-momentum.'
          WRITE (out,190) Sigma_max(isUvel,ng), 'Sigma_max(isVvel)',    &
     &            'Maximum STD value (m/s)', ', 3D V-momentum.'
          WRITE (out,190) Sigma_ml(isVvel,ng), 'Sigma_ml(isVvel)',      &
     &            'Minimum STD at mixed layer (m/s)', ', 3D V-momentum.'
          WRITE (out,190) Sigma_do(isVvel,ng), 'Sigma_do(isVvel)',      &
     &            'Minimum STD in deep ocean (m/s)', ', 3D V-momentum.'
          WRITE (out,130) Sigma_dz(isVvel,ng), 'Sigma_dz(isVvel)',      &
     &            'STD vertical displacement (m)', ', 3D V-momentum.'
          DO itrc=1,NT(ng)
            WRITE (out,190) Sigma_max(isTvar(itrc),ng),                 &
     &            'Sigma_max(idTvar)',                                  &
     &            'Maximum STD Value, ',                                &
     &            TRIM(Vname(1,idTvar(itrc)))
            WRITE (out,190) Sigma_ml(isTvar(itrc),ng),                  &
     &            'Sigma_ml(idTvar)',                                   &
     &            'Minimum STD at mixed layer, ',                       &
     &            TRIM(Vname(1,idTvar(itrc)))
            WRITE (out,190) Sigma_do(isTvar(itrc),ng),                  &
     &            'Sigma_do(idTvar)',                                   &
     &            'Minimum STD in deep ocean, ',                        &
     &            TRIM(Vname(1,idTvar(itrc)))
            WRITE (out,130) Sigma_dz(isTvar(itrc),ng),                  &
     &            'Sigma_dz(idTvar)',                                   &
     &            'STD vertical displacement (m), ',                    &
     &            TRIM(Vname(1,idTvar(itrc)))
          END DO
#     endif
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
     &        'Quality control reject squared threshold, free-surface.'
#     ifndef SOLVE3D
            WRITE (out,100) S_bgqc(isUbar,ng), 'S_bgqc(isUbar)',        &
     &        'Quality control reject squared threshold, 2D U-momentum.'
            WRITE (out,100) S_bgqc(isVbar,ng), 'S_bgqc(isVbar)',        &
     &        'Quality control reject squared threshold, 2D V-momentum.'
#     else
            WRITE (out,100) S_bgqc(isUvel,ng), 'S_bgqc(isUvel)',        &
     &        'Quality control reject squared threshold, 3D U-momentum.'
            WRITE (out,100) S_bgqc(isVvel,ng), 'S_bgqc(isVvel)',        &
     &        'Quality control reject squared threshold, 3D V-momentum.'
            DO itrc=1,NT(ng)
              WRITE (out,190) S_bgqc(isTvar(itrc),ng),                  &
     &            'S_bgqc(idTvar)',                                     &
     &            'Quality control reject squared threshold, ',         &
     &            TRIM(Vname(1,idTvar(itrc)))
            END DO
#     endif
          ELSE IF (bgqc_type(ng).eq.2) THEN
            WRITE (out,80) Nprovenance(ng), 'Nprovenance',              &
     &            'Number of provenances to Quality control.'
            DO i=1,Nprovenance(ng)
              WRITE (out,200) P_bgqc(i,ng), 'P_bgqc', i,                &
     &       'Quality control reject squared threshold for provenance', &
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
          WRITE (out,170) (Lobc(i,isFsur,ng),i=1,4), 'Lobc(isFsur)',    &
     &            'Adjusting free-surface boundaries.'
          WRITE (out,170) (Lobc(i,isUbar,ng),i=1,4), 'Lobc(isUbar)',    &
     &            'Adjusting 2D U-momentum boundaries.'
          WRITE (out,170) (Lobc(i,isVbar,ng),i=1,4), 'Lobc(isVbar)',    &
     &            'Adjusting 2D V-momentum boundaries.'
#     ifdef SOLVE3D
          WRITE (out,170) (Lobc(i,isUvel,ng),i=1,4), 'Lobc(isUvel)',    &
     &            'Adjusting 3D U-momentum boundaries.'
          WRITE (out,170) (Lobc(i,isVvel,ng),i=1,4), 'Lobc(isVvel)',    &
     &            'Adjusting 3D V-momentum boundaries.'
          DO itrc=1,NT(ng)
            WRITE (out,180) (Lobc(i,isTvar(itrc),ng),i=1,4),            &
     &                      'Lobc(isTvar)',                             &
     &            'Adjusting boundaries for tracer ', itrc,             &
     &            TRIM(Vname(1,idTvar(itrc)))
          END DO
#     endif
#    endif
#   endif
#  endif
          IF (NextraObs.gt.0) THEN
            WRITE (out,80) NextraObs, 'NextraObs',                      &
     &            'Number of extra-observations classes to process.'
            DO i=1,NextraObs
              WRITE (out,85) ExtraIndex(i), 'ExtraIndex', i,            &
     &            'Extra-observation type index for: '//                &
     &            TRIM(ExtraName(i))
            END DO
          END IF

#  ifdef DIRAC
          WRITE (out,'(/,a,/)')                                         &
     &          ' Dirac Delta Functions (i,j,k) Locations:'
          WRITE (out,205) TRIM(StateVarName(isFsur)),                   &
     &                    Dirac(1:2,isFsur,ng)
          WRITE (out,205) TRIM(StateVarName(isUbar)),                   &
     &                    Dirac(1:2,isUbar,ng)
          WRITE (out,205) TRIM(StateVarName(isVbar)),                   &
     &                    Dirac(1:2,isVbar,ng)
#   ifdef ADJUST_WSTRESS
          WRITE (out,205) TRIM(StateVarName(isUstr)),                   &
     &                    Dirac(1:2,isUstr,ng)
          WRITE (out,205) TRIM(StateVarName(isVstr)),                   &
     &                    Dirac(1:2,isVstr,ng)
#   endif
#   ifdef SOLVE3D
          WRITE (out,205) TRIM(StateVarName(isUvel)),                   &
     &                    Dirac(1:3,isUvel,ng)
          WRITE (out,205) TRIM(StateVarName(isVvel)),                   &
     &                    Dirac(1:3,isVvel,ng)
          DO itrc=1,NT(ng)
            WRITE (out,205) TRIM(StateVarName(isTvar(itrc))),           &
     &                      Dirac(1:3,isTvar(itrc),ng)
          END DO
#    ifdef ADJUST_STFLUX
          DO itrc=1,NT(ng)
            WRITE (out,205) TRIM(StateVarName(isTsur(itrc))),           &
     &                      Dirac(1:2,isTsur(itrc),ng)
          END DO
#    endif
#   endif
#  endif
        END IF

#  if defined DISTRIBUTE && defined MULTI_SCALE_B
        CALL mp_bcasti (ng, model, exit_flag)
#  endif
      END DO REPORT_LOOP
!
!-----------------------------------------------------------------------
!  Report input files and check availability of input files.
!-----------------------------------------------------------------------
!
      IF (Master.and.Lwrite) THEN
        WRITE (out,150)
#  ifdef VERIFICATION
        WRITE (out,160) '     Verification Parameters File:  ',         &
     &                  TRIM(aparnam)
#  else
        WRITE (out,160) '     Assimilation Parameters File:  ',         &
     &                  TRIM(aparnam)
#  endif
      END IF
!
      DO ng=1,Ngrids

#  if defined FOUR_DVAR           || \
     (defined HESSIAN_SV          && defined BNORM)
#   if defined I4DVAR             || defined OBS_SENSITIVITY || \
       defined OPT_OBSERVATIONS   || defined WEAK_CONSTRAINT
#    if defined RBL4DVAR          || defined R4DVAR          || \
        defined SENSITIVITY_4DVAR || defined SP4DVAR         || \
        defined TL_RBL4DVAR       || defined TL_R4DVAR
        fname=STD(2,ng)%name
        IF (NSA.eq.2) THEN
          IF (.not.find_file(ng, out, fname, 'STDnameM')) THEN
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          ELSE
            IF (Master.and.Lwrite) WRITE (out,160)                      &
     &        '                   Model STD File:  ', TRIM(fname)
          END IF
        END IF
#    endif
#    ifdef STD_MODEL
        fname=STD(5,ng)%name
        IF (Master.and.Lwrite) WRITE (out,160)                          &
     &      '     Computed/Modeled IC STD File:  ', TRIM(fname)
#    else
        fname=STD(1,ng)%name
        IF (.not.find_file(ng, out, fname, 'STDnameI')) THEN
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        ELSE
          IF (Master.and.Lwrite) WRITE (out,160)                        &
     &      '      Initial Conditions STD File:  ', TRIM(fname)
        END IF
#    endif
#    ifdef ADJUST_BOUNDARY
        fname=STD(3,ng)%name
        IF (.not.find_file(ng, out, fname, 'STDnameB')) THEN
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        ELSE
          IF (Master.and.Lwrite) WRITE (out,160)                        &
     &      '     Boundary Conditions STD File:  ', TRIM(fname)
        END IF
#    endif
#    if defined ADJUST_WSTRESS || defined ADJUST_STFLUX
        fname=STD(4,ng)%name
        IF (.not.find_file(ng, out, fname, 'STDnameF')) THEN
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        ELSE
          IF (Master.and.Lwrite) WRITE (out,160)                        &
     &      '         Surface Forcing STD File:  ', TRIM(fname)
        END IF
#    endif
#   endif
#   if defined RBL4DVAR          || defined R4DVAR    || \
       defined SENSITIVITY_4DVAR || defined SP4DVAR   || \
       defined TL_RBL4DVAR       || defined TL_R4DVAR
        fname=NRM(2,ng)%name
        IF (Master.and.Lwrite) WRITE (out,160)                          &
     &    '                  Model Norm File:  ', TRIM(fname)
        IF (.not.LdefNRM(2,ng)) THEN
          IF (.not.find_file(ng, out, fname, 'NRMnameM')) THEN
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END IF
#   elif defined CORRELATION
        fname=NRM(2,ng)%name
        IF (Master.and.Lwrite) WRITE (out,160)                          &
     &    '                  Model Norm File:  ', TRIM(fname)
        IF (.not.LdefNRM(2,ng).and.LwrtNRM(2,ng)) THEN
          IF (.not.find_file(ng, out, fname, 'NRMnameM')) THEN
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END IF
#   endif
#   if defined CORRELATION
        fname=NRM(1,ng)%name
        IF (Master.and.Lwrite) WRITE (out,160)                          &
     &    '     Initial Conditions Norm File:  ', TRIM(fname)
        IF (.not.LdefNRM(1,ng).and.LwrtNRM(1,ng)) THEN
          IF (.not.find_file(ng, out, fname, 'NRMnameI')) THEN
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END IF
#   else
        fname=NRM(1,ng)%name
        IF (Master.and.Lwrite) WRITE (out,160)                          &
     &    '     Initial Conditions Norm File:  ', TRIM(fname)
        IF (.not.LdefNRM(1,ng)) THEN
          IF (.not.find_file(ng, out, fname, 'NRMnameI')) THEN
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END IF
#   endif
#   ifdef ADJUST_BOUNDARY
#    ifdef CORRELATION
        fname=NRM(3,ng)%name
        IF (Master.and.Lwrite) WRITE (out,160)                          &
     &    '    Boundary Conditions Norm File:  ', TRIM(fname)
        IF (.not.LdefNRM(3,ng).and.LwrtNRM(3,ng)) THEN
          IF (.not.find_file(ng, out, fname, 'NRMnameB')) THEN
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END IF
#    else
        fname=NRM(3,ng)%name
        IF (Master.and.Lwrite) WRITE (out,160)                          &
     &    '    Boundary Conditions Norm File:  ', TRIM(fname)
        IF (.not.LdefNRM(3,ng)) THEN
          IF (.not.find_file(ng, out, fname, 'NRMnameB')) THEN
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END IF
#    endif
#   endif
#   if defined ADJUST_WSTRESS || defined ADJUST_STFLUX
#    ifdef CORRELATION
        fname=NRM(4,ng)%name
        IF (Master.and.Lwrite) WRITE (out,160)                          &
     &    '        Surface Forcing Norm File:  ', TRIM(fname)
        IF (.not.LdefNRM(4,ng).and.LwrtNRM(4,ng)) THEN
          IF (.not.find_file(ng, out, fname, 'NRMnameF')) THEN
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END IF
#    else
        fname=NRM(4,ng)%name
        IF (Master.and.Lwrite) WRITE (out,160)                          &
     &    '        Surface Forcing Norm File:  ', TRIM(fname)
        IF (.not.LdefNRM(4,ng)) THEN
          IF (.not.find_file(ng, out, fname, 'NRMnameF')) THEN
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END IF
#    endif
#   endif
#   if !(defined CORRELATION || defined OPT_OBSERVATIONS)
        fname=OBS(ng)%name
        IF (.not.find_file(ng, out, fname, 'OBSname')) THEN
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        ELSE
          IF (Master.and.Lwrite) WRITE (out,160)                        &
     &      '          Input observations File:  ', TRIM(fname)
        END IF
#   endif
#   if !defined CORRELATION
        IF (Master.and.Lwrite) WRITE (out,160)                          &
     &    '        Input/Output Lanczos File:  ', TRIM(LCZ(ng)%name)
#    ifndef I4DVAR_ANA_SENSITIVITY
        IF (Master.and.Lwrite) WRITE (out,160)                          &
     &    '        Input/Output Hessian File:  ', TRIM(HSS(ng)%name)
#    endif
#    ifdef EVOLVED_LCZ
        IF (Master.and.Lwrite) WRITE (out,160)                          &
          '      Output evolved Lanczos File:  ', TRIM(LZE(ng)%name)
#    endif
#   endif
#  endif
#  if !defined CORRELATION
#   ifdef VERIFICATION
        fname=OBS(ng)%name
        IF (.not.find_file(ng, out, fname, 'OBSname')) THEN
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        ELSE
          IF (Master.and.Lwrite) WRITE (out,160)                        &
     &      '          Input observations File:  ', TRIM(fname)
        END IF
        IF (Master.and.Lwrite) WRITE (out,160)                          &
     &    '         Output verification File:  ', TRIM(DAV(ng)%name)
#   else
#    ifdef SP4DVAR
        IF (Master.and.Lwrite) THEN
          WRITE (out,160) '    Input/Output TLM Arnoldi File:  ',       &
     &                    TRIM(SPT(ng)%name)
          WRITE (out,160) '    Input/Output ADM Arnoldi File:  ',       &
     &                    TRIM(SPA(ng)%name)
          WRITE (out,160) '    Input/Output TLM Scratch File:  ',       &
     &                    TRIM(SCT(ng)%name)
          WRITE (out,160) '    Input/Output ADM Scratch File:  ',       &
     &                    TRIM(SCA(ng)%name)
        END IF
#    endif
        IF (Master.and.Lwrite) WRITE (out,160)                          &
     &    '               Output 4D-Var File:  ', TRIM(DAV(ng)%name)
#   endif
#  endif
#  if defined WEAK_CONSTRAINT   && \
     (defined POSTERIOR_ERROR_F || defined POSTERIOR_ERROR_I)
        IF (Master.and.Lwrite) WRITE (out,160)                          &
          '      Output Posterior Error File:  ', TRIM(ERR(ng)%name)
#  endif
#  if defined RBL4DVAR_FCT_SENSITIVITY && defined OBS_SPACE
        fname=OIFA(ng)%name
        IF (Master.and.Lwrite) WRITE (out,160)                          &
     &    '  Input Obs Impacts Analysis File:  ', TRIM(fname)
!
        fname=OIFB(ng)%name
        IF (Master.and.Lwrite) WRITE (out,160)                          &
     &    'Input Obs Impacts Background File:  ', TRIM(fname)
#  endif
#  if defined SPLIT_RBL4DVAR
        fname=IMP(ng)%name
        IF (Master.and.Lwrite) WRITE (out,160)                          &
     &    ' Input Coarser Grid Impulses File:  ', TRIM(fname)
!
        fname=INC(ng)%name
        IF (Master.and.Lwrite) WRITE (out,160)                          &
     &    'Input Coarser Grid Icrements File:  ', TRIM(fname)

#  endif
      END DO

#  if defined WEAK_CONSTRAINT && defined RPCG
!
!  Stop if activating pre-conditioning for the RBLanczos minimization
!  algorithm. It does not work yet.
!
      DO ng=1,Ngrids
        IF (Lprecond) THEN
          IF (Master) THEN
            WRITE (out,230) 'Lprecond', Lprecond,                       &
     &            'pre-conditioning does not work yet with ' //         &
     &            uppercase('rpcg') // ' algorithm.',                   &
     &            'Set Lprecond to F in ' // TRIM(aparnam)
          END IF
          exit_flag=5
          RETURN
        END IF
      END DO
#  endif
!
!-----------------------------------------------------------------------
!  Convert time scales to seconds.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
#  ifdef FOUR_DVAR
         timeIAU(ng)=timeIAU(ng)*86400.0_dp
#  endif
#  if defined WEAK_CONSTRAINT && defined TIME_CONV
        DO i=1,MstateVar
          Tdecay(i,ng)=Tdecay(i,ng)*86400.0_r8
        END DO
#  endif
      END DO
!
  50  FORMAT (/,' READ_AssPar - Error while processing line: ',/,a)
  55  FORMAT (/,' READ_AssPar - ',a,i4,2x,i4,/,15x,a)
#  ifdef VERIFICATION
  60  FORMAT (/,/,' Observation Parameters, Grid: ',i2.2,               &
     &        /,  ' ================================',/)
#  else
  60  FORMAT (/,/,' Assimilation Parameters, Grid: ',i2.2,              &
     &        /,  ' =================================',/)
#  endif
  70  FORMAT (10x,l1,2x,a,t32,a)
  80  FORMAT (1x,i10,2x,a,t32,a)
  85  FORMAT (1x,i10,2x,a,'(',i2.2,')',t32,a)
  90  FORMAT (1x,i10,2x,a,t32,a,i0,a)
  95  FORMAT (1x,a,2x,a,t32,a)
 100  FORMAT (1p,e11.4,2x,a,t32,a)
 110  FORMAT (10x,l1,2x,a,t32,a,i2.2,':',1x,a)
 120  FORMAT (f11.3,2x,a,t32,a)
 125  FORMAT (f11.3,2x,a,t32,a,i0,a)
 130  FORMAT (f11.3,2x,a,t32,a,a,'.')
#  ifdef VERIFICATION
 150  FORMAT (/,' Input/Output Verification Files:',/)
#  else
 150  FORMAT (/,' Input/Output Assimilation Files:',/)
#  endif
 160  FORMAT (2x,a,a)
 170  FORMAT (3x,4(1x,l1),2x,a,t32,a)
 180  FORMAT (3x,4(1x,l1),2x,a,t32,a,i2.2,':',1x,a)
 190  FORMAT (1p,e11.4,2x,a,t32,a,a,'.')
 200  FORMAT (1p,e11.4,2x,a,'(',i2.2,')',t32,a,':',1x,i6)
 205  FORMAT (t14,a,':',t32,3i6)
 210  FORMAT (/,' READ_ASSPAR - variable info not yet loaded, ', a)
 220  FORMAT (/,' READ_ASSPAR - Grid ', i2.2,                           &
     &        ', could not find input file:  ',a)
 230  FORMAT (/,' READ_ASSPAR - Illegal parameter, ', a, ' = ', 1x, l1, &
     &        /,15x,a,/,15x,a)
 240  FORMAT (/,' READ_ASSPAR - Illegal parameter, ', a, ' = ', 1x, i0, &
     &        /,15x,a,/,15x,'correct input script: ',1x,'''',a,'''',/)
# endif
!
      RETURN
      END SUBROUTINE read_AssPar
#else
      SUBROUTINE read_AssPar
      END SUBROUTINE read_AssPar
#endif
