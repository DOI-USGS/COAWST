#include "cppdefs.h"
#undef REF_DIFFUSION
#undef REF_FLUX_FILTER
      MODULE step3dref_t_mod
#if !defined TS_FIXED && (defined NONLINEAR && defined SOLVE3D) && \
     defined REFINED_GRID
!
!svn $Id: step3d_t.F 732 2008-09-07 01:55:51Z jcwarner $
!=======================================================================
!  Copyright (c) 2002-2010 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                           Hernan G. Arango   !
!========================================== Alexander F. Shchepetkin ===
!                                                                      !
!  This routine time-steps tracer equations. Notice that advective     !
!  and diffusive terms are time-stepped differently.                   !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: step3dref_t

      CONTAINS
!
!***********************************************************************
      SUBROUTINE step3dref_t (ng, tile)
!***********************************************************************
!
      USE mod_param
# ifdef CLIMATOLOGY
      USE mod_clima
# endif
# ifdef DIAGNOSTICS_TS
      USE mod_diags
# endif
      USE mod_grid
      USE mod_mixing
      USE mod_ocean
      USE mod_refined
      USE mod_stepping
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
      integer :: ngp, ngc, LBic, UBic, LBjc, UBjc
# include "tile.h"
      ngc=mychild(ng)
      LBic=BOUNDS(ngc)%LBi(tile)
      UBic=BOUNDS(ngc)%UBi(tile)
      LBjc=BOUNDS(ngc)%LBj(tile)
      UBjc=BOUNDS(ngc)%UBj(tile)
      ngp=myparent(ng)
!
# ifdef PROFILE
      CALL wclock_on (ng, iNLM, 35)
# endif
      CALL step3dref_t_tile (ng, tile,                                  &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    LBic, UBic, LBjc, UBjc,                       &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    nrhs(ng), nstp(ng), nnew(ng),                 &
     &                    REFINED(ngc) % ref_scale,                     &
# if defined TCLM_NUDGING && defined TCLIMATOLOGY
     &                    CLIMA(ng) % Tnudgcof,                         &
     &                    CLIMA(ng) % tclm,                             &
# endif
# ifdef MASKING
     &                    GRID(ng) % rmask,                             &
# endif
     &                    GRID(ng) % pm,                                &
     &                    GRID(ng) % pn,                                &
     &                    GRID(ng) % Hz,                                &
     &                    GRID(ng) % z_r,                               &
     &                    MIXING(ng) % Akt,                             &
     &                    OCEAN(ng) % hindcu,    OCEAN(ng) % hindcv,    &
     &                    OCEAN(ng) % hindci_u,  OCEAN(ng) % hindcj_u,  &
     &                    OCEAN(ng) % hindci_v,  OCEAN(ng) % hindcj_v,  &
     &                    OCEAN(ng) % FXp_west,  OCEAN(ng) % FXp_east,  &
     &                    OCEAN(ng) % FEp_south, OCEAN(ng) % FEp_north, &
     &                    OCEAN(ngc) % FXc_west, OCEAN(ngc) % FXc_east, &
     &                    OCEAN(ngc) % FEc_south,OCEAN(ngc) % FEc_north,&
# ifdef DIAGNOSTICS_TS
     &                    DIAGS(ng) % DiaTwrk,                          &
# endif
     &                    OCEAN(ng) % t)

# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 35)
# endif
      RETURN
      END SUBROUTINE step3dref_t
!
!***********************************************************************
      SUBROUTINE step3dref_t_tile (ng, tile,                            &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          LBic, UBic, LBjc, UBjc,                 &
     &                          IminS, ImaxS, JminS, JmaxS,             &
     &                          nrhs, nstp, nnew,                       &
     &                          ref_scale,                              &
# if defined TCLM_NUDGING && defined TCLIMATOLOGY
     &                          Tnudgcof, tclm,                         &
# endif
# ifdef MASKING
     &                          rmask,                                  &
# endif
     &                          pm, pn,                                 &
     &                          Hz, z_r,                                &
     &                          Akt,                                    &
     &                          hindcu,    hindcv,                      &
     &                          hindci_u,  hindcj_u,                    &
     &                          hindci_v,  hindcj_v,                    &
     &                          FXp_west,  FXp_east,                    &
     &                          FEp_south, FEp_north,                   &
     &                          FXc_west,  FXc_east,                    &
     &                          FEc_south, FEc_north,                   &
# ifdef DIAGNOSTICS_TS
     &                          DiaTwrk,                                &
# endif
     &                          t)
!***********************************************************************
!
      USE mod_ncparam
      USE mod_parallel
      USE mod_param
      USE mod_scalars
!
# if defined EW_PERIODIC || defined NS_PERIODIC
      USE exchange_3d_mod, ONLY : exchange_r3d_tile
# endif
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange4d
# endif
      USE t3dbc_mod, ONLY : t3dbc_tile
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: LBic, UBic, LBjc, UBjc
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nrhs, nstp, nnew
      integer, intent(in) :: ref_scale
!
# ifdef ASSUMED_SHAPE
#  if defined TCLM_NUDGING && defined TCLIMATOLOGY
      real(r8), intent(in) :: Tnudgcof(LBi:,LBj:,:)
      real(r8), intent(in) :: tclm(LBi:,LBj:,:,:)
#  endif
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
#  ifdef SUN
      real(r8), intent(in) :: Akt(LBi:UBi,LBj:UBj,0:N(ng),NAT)
#  else
      real(r8), intent(in) :: Akt(LBi:,LBj:,0:,:)
#  endif
      integer, intent(in) :: hindcu(:)
      integer, intent(in) :: hindcv(:)
      real(r8), intent(in) :: hindci_u(LBi:,LBj:)
      real(r8), intent(in) :: hindcj_u(LBi:,LBj:)
      real(r8), intent(in) :: hindci_v(LBi:,LBj:)
      real(r8), intent(in) :: hindcj_v(LBi:,LBj:)
      real(r8), intent(in) :: FXp_west(LBj:,:,:)
      real(r8), intent(in) :: FXp_east(LBj:,:,:)
      real(r8), intent(in) :: FEp_south(LBi:,:,:)
      real(r8), intent(in) :: FEp_north(LBi:,:,:)
      real(r8), intent(in) :: FXc_west(LBjc:,:)
      real(r8), intent(in) :: FXc_east(LBjc:,:)
      real(r8), intent(in) :: FEc_south(LBic:,:)
      real(r8), intent(in) :: FEc_north(LBic:,:)
#  ifdef DIAGNOSTICS_TS
      real(r8), intent(inout) :: DiaTwrk(LBi:,LBj:,:,:,:)
#  endif
#  ifdef SUN
      real(r8), intent(inout) :: t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
#  else
      real(r8), intent(inout) :: t(LBi:,LBj:,:,:,:)
#  endif

# else

#  if defined TCLM_NUDGING && defined TCLIMATOLOGY
      real(r8), intent(in) :: Tnudgcof(LBi:UBi,LBj:UBj,NT(ng))
      real(r8), intent(in) :: tclm(LBi:UBi,LBj:UBj,N(ng),NT(ng))
#  endif
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pn(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: Akt(LBi:UBi,LBj:UBj,0:N(ng),NAT)
      integer, intent(in) :: hindcu(4)
      integer, intent(in) :: hindcv(4)
      real(r8), intent(in) :: hindci_u(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: hindcj_u(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: hindci_v(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: hindcj_v(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: FXp_west(LBj:UBj,N(ng),NT(ng))
      real(r8), intent(in) :: FXp_east(LBj:UBj,N(ng),NT(ng))
      real(r8), intent(in) :: FEp_south(LBi:UBi,N(ng),NT(ng))
      real(r8), intent(in) :: FEp_north(LBi:UBi,N(ng),NT(ng))
      real(r8), intent(in) :: FXc_west(LBjc:UBjc,NT(ng))
      real(r8), intent(in) :: FXc_east(LBjc:UBjc,NT(ng))
      real(r8), intent(in) :: FEc_south(LBic:UBic,NT(ng))
      real(r8), intent(in) :: FEc_north(LBic:UBic,NT(ng))
#  ifdef DIAGNOSTICS_TS
      real(r8), intent(inout) :: DiaTwrk(LBi:UBi,LBj:UBj,N(ng),NT(ng),  &
     &                                   NDT)
#  endif
      real(r8), intent(inout) :: t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))

# endif
!
!  Local variable declarations.
!
# ifdef DISTRIBUTE
#  ifdef EW_PERIODIC
      logical :: EWperiodic=.TRUE.
#  else
      logical :: EWperiodic=.FALSE.
#  endif
#  ifdef NS_PERIODIC
      logical :: NSperiodic=.TRUE.
#  else
      logical :: NSperiodic=.FALSE.
#  endif
# endif
      integer :: i, is, itrc, j, k, ltrc
      integer :: ie, js, je, kk, Npts
      integer :: ic, jc, ref, ngc
# ifdef DIAGNOSTICS_TS
      integer :: idiag
# endif
      integer :: ISTAT(MPI_STATUS_SIZE)
      integer :: MySize, MyError, MyTag
      real(r8), parameter :: eps = 1.0E-16_r8
      real(r8) :: cff, cff1, cff2, cff3

      real(r8), dimension(0:N(ng)) :: CF
      real(r8), dimension(0:N(ng)) :: BC
      real(r8), dimension(0:N(ng)) :: DC
      real(r8), dimension(0:N(ng)) :: FC
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,N(ng)) :: oHz
      real(r8), allocatable :: FXc_west_all(:)
      real(r8), allocatable :: FXc_east_all(:)
      real(r8), allocatable :: FEc_south_all(:)
      real(r8), allocatable :: FEc_north_all(:)
      real(r8), allocatable :: fx_exch(:)

# include "set_bounds.h"
      ngc=mychild(ng)
      Istrc=BOUNDS(ngc)%Istr(tile)
      Iendc=BOUNDS(ngc)%Iend(tile)
      Jstrc=BOUNDS(ngc)%Jstr(tile)
      Jendc=BOUNDS(ngc)%Jend(tile)
!
!-----------------------------------------------------------------------
!  Time-step horizontal advection term.
!-----------------------------------------------------------------------
!
!  Compute inverse thickness.
!
      DO k=1,N(ng)
        DO j=Jstr-2,Jend+2
          DO i=Istr-2,Iend+2
            oHz(i,j,k)=1.0_r8/Hz(i,j,k)
          END DO
        END DO
      END DO
      T_LOOP : DO itrc=1,NT(ng)
!
!
!  NOTE : this mp gather send stuff should go to another routine.
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  CHILD WESTERN EDGE - Accumlate and scatter child perimeter flux
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
        Npts=Mm(ngc)
        allocate (FXc_west_all(Npts))
        allocate (fx_exch(Npts))
!
!  Pack FX child data for send
!
        FXc_west_all=0.0_r8
        fx_exch=0.0_r8
        DO j=Jstrc,Jendc
          fx_exch(j)=FXc_west(j,itrc)
        END DO

        IF (MyRank.eq.MyMaster) THEN
!
! Gather all the FX child data from all processors and add them up. 
! If a processor did not contribute, then its value is zero.
!
          DO is=1,NtileI(ngc)*NtileJ(ngc)-1
            MyTag=is
            CALL MPI_RECV(fx_exch,Npts,MP_FLOAT,is,                      &
     &                    MyTag,OCN_COMM_WORLD,ISTAT,MyError)
            DO j=1,Npts
              FXc_west_all(j)=FXc_west_all(j)+fx_exch(j)
            END DO
          END DO
!  add master node contribution
          DO j=Jstrc,Jendc
            FXc_west_all(j)=FXc_west_all(j)+FXc_west(j,itrc)
          END DO
        ELSE
          MyTag=tile
          CALL MPI_SEND(fx_exch,Npts,MP_FLOAT,0,MyTag,          &
     &                  OCN_COMM_WORLD,MyError)
        END IF
!
!  Scatter FX data to all the nodes.
!
        IF (MyRank.eq.MyMaster) THEN
          DO is=1,Npts
            fx_exch(is)=FXc_west_all(is)
          END DO
        ENDIF
        CALL MPI_BCAST(fx_exch,Npts,MP_FLOAT,0,                   &
     &                 OCN_COMM_WORLD,MyError)
        DO is=1,Npts
          FXc_west_all(is)=fx_exch(is)
        END DO
# ifdef REF_FLUX_FILTER
!
!  Introduce some spatial averaging.
!
        ref=(ref_scale-1)/2
        cff1=0.0_r8
        DO is=1,Npts
          cff1=cff1+FXc_west_all(is)
        END DO
        DO is=1+ref,Npts-ref
          cff=0.0_r8
          DO k=-ref,ref
            cff=cff+FXc_west_all(is+k)
          END DO
          fx_exch(is)=cff/REAL(ref_scale,r8)
        END DO
        cff=0.0_r8
        DO is=1,Npts
          cff=cff+fx_exch(is)
        END DO
        DO is=1,Npts
          FXc_west_all(is)=fx_exch(is)*cff1/cff
        END DO
# endif
        deallocate (fx_exch)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  CHILD EASTERN EDGE - Accumlate and scatter child perimeter flux
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
        Npts=Mm(ngc)
        allocate (FXc_east_all(Npts))
        allocate (fx_exch(Npts))
!
!  Pack FX child data for send
!
        FXc_east_all=0.0_r8
        fx_exch=0.0_r8
        DO j=Jstrc,Jendc
          fx_exch(j)=FXc_east(j,itrc)
        END DO

        IF (MyRank.eq.MyMaster) THEN
!
! Gather all the FX child data from all processors and add them up. 
! If a processor did not contribute, then its value is zero.
!
          DO is=1,NtileI(ngc)*NtileJ(ngc)-1
            MyTag=is
            CALL MPI_RECV(fx_exch,Npts,MP_FLOAT,is,                      &
     &                    MyTag,OCN_COMM_WORLD,ISTAT,MyError)
            DO j=1,Npts
              FXc_east_all(j)=FXc_east_all(j)+fx_exch(j)
            END DO
          END DO
!  add master node contribution
          DO j=Jstrc,Jendc
            FXc_east_all(j)=FXc_east_all(j)+FXc_east(j,itrc)
          END DO
        ELSE
          MyTag=tile
          CALL MPI_SEND(fx_exch,Npts,MP_FLOAT,0,MyTag,                  &
     &                  OCN_COMM_WORLD,MyError)
        END IF
!
!  Scatter FX data to all the nodes.
!
        IF (MyRank.eq.MyMaster) THEN
          DO is=1,Npts
            fx_exch(is)=FXc_east_all(is)
          END DO
        ENDIF
        CALL MPI_BCAST(fx_exch,Npts,MP_FLOAT,0,                         &
     &                 OCN_COMM_WORLD,MyError)
        DO is=1,Npts
          FXc_east_all(is)=fx_exch(is)
        END DO
# ifdef REF_FLUX_FILTER
!
!  Introduce some spatial averaging.
!
        ref=(ref_scale-1)/2
        cff1=0.0_r8
        DO is=1,Npts
          cff1=cff1+FXc_east_all(is)
        END DO
        DO is=1+ref,Npts-ref
          cff=0.0_r8
          DO k=-ref,ref
            cff=cff+FXc_east_all(is+k)
          END DO
          fx_exch(is)=cff/REAL(ref_scale,r8)
        END DO
        cff=0.0_r8
        DO is=1,Npts
          cff=cff+fx_exch(is)
        END DO
        DO is=1,Npts
          FXc_east_all(is)=fx_exch(is)*cff1/cff
        END DO
# endif
        deallocate (fx_exch)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  CHILD SOUTHERN EDGE - Accumlate and scatter child perimeter flux
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
        Npts=Lm(ngc)
        allocate (FEc_south_all(Npts))
        allocate (fx_exch(Npts))
!
!  Pack FX child data for send
!
        FEc_south_all=0.0_r8
        fx_exch=0.0_r8
        DO i=Istrc,Iendc
          fx_exch(i)=FEc_south(i,itrc)
        END DO

        IF (MyRank.eq.MyMaster) THEN
!
! Gather all the FX child data from all processors and add them up. 
! If a processor did not contribute, then its value is zero.
!
          DO is=1,NtileI(ngc)*NtileJ(ngc)-1
            MyTag=is
            CALL MPI_RECV(fx_exch,Npts,MP_FLOAT,is,                      &
     &                    MyTag,OCN_COMM_WORLD,ISTAT,MyError)
            DO i=1,Npts
              FEc_south_all(i)=FEc_south_all(i)+fx_exch(i)
            END DO
          END DO
!  add master node contribution
          DO i=Istrc,Iendc
            FEc_south_all(i)=FEc_south_all(i)+FEc_south(i,itrc)
          END DO
        ELSE
          MyTag=tile
          CALL MPI_SEND(fx_exch,Npts,MP_FLOAT,0,MyTag,          &
     &                  OCN_COMM_WORLD,MyError)
        END IF
!
!  Scatter FX data to all the nodes.
!
        IF (MyRank.eq.MyMaster) THEN
          DO is=1,Npts
            fx_exch(is)=FEc_south_all(is)
          END DO
        ENDIF
        CALL MPI_BCAST(fx_exch,Npts,MP_FLOAT,0,                   &
     &                 OCN_COMM_WORLD,MyError)
        DO is=1,Npts
          FEc_south_all(is)=fx_exch(is)
        END DO
# ifdef REF_FLUX_FILTER
!
!  Introduce some spatial averaging.
!
        ref=(ref_scale-1)/2
        cff1=0.0_r8
        DO is=1,Npts
          cff1=cff1+FEc_south_all(is)
        END DO
        DO is=1+ref,Npts-ref
          cff=0.0_r8
          DO k=-ref,ref
            cff=cff+FEc_south_all(is+k)
          END DO
          fx_exch(is)=cff/REAL(ref_scale,r8)
        END DO
        cff=0.0_r8
        DO is=1,Npts
          cff=cff+fx_exch(is)
        END DO
        DO is=1,Npts
          FEc_south_all(is)=fx_exch(is)*cff1/cff
        END DO
# endif
        deallocate (fx_exch)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  CHILD NORTHERN EDGE - Accumlate and scatter child perimeter flux
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
        Npts=Lm(ngc)
        allocate (FEc_north_all(Npts))
        allocate (fx_exch(Npts))
!
!  Pack FX child data for send
!
        FEc_north_all=0.0_r8
        fx_exch=0.0_r8
        DO i=Istrc,Iendc
          fx_exch(i)=FEc_north(i,itrc)
        END DO

        IF (MyRank.eq.MyMaster) THEN
!
! Gather all the FX child data from all processors and add them up. 
! If a processor did not contribute, then its value is zero.
!
          DO is=1,NtileI(ngc)*NtileJ(ngc)-1
            MyTag=is
            CALL MPI_RECV(fx_exch,Npts,MP_FLOAT,is,                      &
     &                    MyTag,OCN_COMM_WORLD,ISTAT,MyError)
            DO i=1,Npts
              FEc_north_all(i)=FEc_north_all(i)+fx_exch(i)
            END DO
          END DO
!  add master node contribution
          DO i=Istrc,Iendc
            FEc_north_all(i)=FEc_north_all(i)+FEc_north(i,itrc)
          END DO
        ELSE
          MyTag=tile
          CALL MPI_SEND(fx_exch,Npts,MP_FLOAT,0,MyTag,          &
     &                  OCN_COMM_WORLD,MyError)
        END IF
!
!  Scatter FX data to all the nodes.
!
        IF (MyRank.eq.MyMaster) THEN
          DO is=1,Npts
            fx_exch(is)=FEc_north_all(is)
          END DO
        ENDIF
        CALL MPI_BCAST(fx_exch,Npts,MP_FLOAT,0,                   &
     &                 OCN_COMM_WORLD,MyError)
        DO is=1,Npts
          FEc_north_all(is)=fx_exch(is)
        END DO
# ifdef REF_FLUX_FILTER
!
!  Introduce some spatial averaging.
!
        ref=(ref_scale-1)/2
        cff1=0.0_r8
        DO is=1,Npts
          cff1=cff1+FEc_north_all(is)
        END DO
        DO is=1+ref,Npts-ref
          cff=0.0_r8
          DO k=-ref,ref
            cff=cff+FEc_north_all(is+k)
          END DO
          fx_exch(is)=cff/REAL(ref_scale,r8)
        END DO
        cff=0.0_r8
        DO is=1,Npts
          cff=cff+fx_exch(is)
        END DO
        DO is=1,Npts
          FEc_north_all(is)=fx_exch(is)*cff1/cff
        END DO
# endif
        deallocate (fx_exch)
!
!   END OF MPI GATHER NOTE
!
!
!  Correct the flux of the parent grid at the location of the child boundary.
!  Compare vertical sum from last time step of parent FXp and FEp
!  to cumulated child fluxes FXb(ng) and FEb(ng) (these are the 
!  child grid fluxes FXc(ngc) and FEc(ngc)).
!
!
!  western and eastern child edges.
!
        js=hindcu(3)
        je=hindcu(4)
        DO j=js,je
!
!  Parent grid at child WESTERN EDGE
!
          is=hindcu(1)-1
          IF (((Istr.le.is-1).and.(is-1.le.Iend)).and.                  &
     &        ((Jstr.le.j).and.(j.le.Jend)).and.                        &
     &        (NINT(hindci_u(is,j)).eq.1)) THEN
!
!    compute vertical sum of parent flux from last time step.
!
            cff1=FXp_west(j,1,itrc)
            DO k=2,N(ng)
              cff1=cff1+FXp_west(j,k,itrc)
            END DO
!
!    compute vertical sum of child flux for nrefined steps.
!    jc is the index of the child grid.
!
            jc=NINT(hindcj_u(is,j))
            ref=(ref_scale-1)/2
            cff3=0.0_r8
            DO k=-ref,ref
              cff3=cff3+FXc_west_all(jc+k)
            END DO
!
!    correct parent flux.
!
            cff=pm(is-1,j)*pn(is-1,j)
            DO k=1,N(ng)
              cff2=(cff3-cff1)*FXp_west(j,k,itrc)/(cff1+eps)
# ifdef REF_DIFFUSION
              t(is-1,j,k,nstp,itrc)=t(is-1,j,k,nstp,itrc)*Hz(is-1,j,k)- &
     &                              cff2*cff
# else
              t(is-1,j,k,nstp,itrc)=t(is-1,j,k,nstp,itrc)-              &
     &                              cff2*cff*oHz(is-1,j,k)
#  if defined TCLM_NUDGING && defined TCLIMATOLOGY
!
!  Nudge towards tracer climatology.
!
              t(is-1,j,k,nstp,itrc)=t(is-1,j,k,nstp,itrc)+              &
     &                              dt(ng)*Tnudgcof(i,j,itrc)*          &
     &                              (tclm(i,j,k,itrc)-                  &
     &                               t(is-1,j,k,nstp,itrc))
#  endif
# endif
# ifdef MASKING
              t(is-1,j,k,nstp,itrc)=t(is-1,j,k,nstp,itrc)*              &
     &                              rmask(is-1,j)
# endif
            END DO
          END IF
!
!  Parent grid at child EASTERN EDGE
!
          is=hindcu(2)+1
          IF (((Istr.le.is).and.(is.le.Iend)).and.                      &
     &        ((Jstr.le.j).and.(j.le.Jend)).and.                        &
     &        (NINT(hindci_u(is,j)).eq.(Lm(ngc)+1))) THEN
!
!    compute vertical sum of parent flux from last time.
!
            cff1=FXp_east(j,1,itrc)
            DO k=2,N(ng)
              cff1=cff1+FXp_east(j,k,itrc)
            END DO
!
!    compute vertical sum of child flux for nrefined steps.
!    jc is the index of the child grid.
!
            jc=NINT(hindcj_u(is,j))
            ref=(ref_scale-1)/2
            cff3=0.0_r8
            DO k=-ref,ref
              cff3=cff3+FXc_east_all(jc+k)
            END DO
!
!    correct parent flux.
!
            cff=pm(is,j)*pn(is,j)
            DO k=1,N(ng)
              cff2=(cff3-cff1)*FXp_east(j,k,itrc)/(cff1+eps)
# ifdef REF_DIFFUSION
              t(is,j,k,nstp,itrc)=t(is,j,k,nstp,itrc)*Hz(is,j,k)+       &
     &                            cff2*cff
# else
              t(is,j,k,nstp,itrc)=t(is,j,k,nstp,itrc)+                  &
     &                            cff2*cff*oHz(is,j,k)
#  if defined TCLM_NUDGING && defined TCLIMATOLOGY
!
!  Nudge towards tracer climatology.
!
              t(is,j,k,nstp,itrc)=t(is,j,k,nstp,itrc)+                  &
     &                              dt(ng)*Tnudgcof(i,j,itrc)*          &
     &                              (tclm(i,j,k,itrc)-                  &
     &                               t(is,j,k,nstp,itrc))
#  endif
# endif
# ifdef MASKING
              t(is,j,k,nstp,itrc)=t(is,j,k,nstp,itrc)*                  &
     &                            rmask(is,j)
# endif
            END DO
          END IF
        END DO
!
!  southern and northern child edges.
!
        is=hindcv(1)
        ie=hindcv(2)
        DO i=is,ie
!
!  Parent grid at child SOUTHERN EDGE
!
          js=hindcv(3)-1
          IF (((Istr.le.i).and.(i.le.Iend)).and.                        &
     &        ((Jstr.le.js-1).and.(js-1.le.Jend)).and.                  &
     &        (NINT(hindcj_v(i,js)).eq.1)) THEN
!
!    compute vertical sum of parent flux from last time.
!
            cff1=FEp_south(i,1,itrc)
            DO k=2,N(ng)
              cff1=cff1+FEp_south(i,k,itrc)
            END DO
!
!    Use cumulated child flux FXc and parent flux FXp from last time step
!    to compute correction for parent.  ic is the index of the child grid.
!
            ic=NINT(hindci_v(i,js))
            ref=(ref_scale-1)/2
            cff3=0.0_r8
            DO k=-ref,ref
              cff3=cff3+FEc_south_all(ic+k)
            END DO
!
!    correct parent flux.
!
            cff=pm(i,js-1)*pn(i,js-1)
            DO k=1,N(ng)
              cff2=(cff3-cff1)*FEp_south(i,k,itrc)/(cff1+eps)
# ifdef REF_DIFFUSION
              t(i,js-1,k,nstp,itrc)=t(i,js-1,k,nstp,itrc)*Hz(i,js-1,k)- &
     &                              cff2*cff
# else
              t(i,js-1,k,nstp,itrc)=t(i,js-1,k,nstp,itrc)-              &
     &                              cff2*cff*oHz(i,js-1,k)
#  if defined TCLM_NUDGING && defined TCLIMATOLOGY
!
!  Nudge towards tracer climatology.
!
              t(i,js-1,k,nstp,itrc)=t(i,js-1,k,nstp,itrc)+              &
     &                              dt(ng)*Tnudgcof(i,j,itrc)*          &
     &                              (tclm(i,j,k,itrc)-                  &
     &                               t(i,js-1,k,nstp,itrc))
#  endif
# endif
# ifdef MASKING
              t(i,js-1,k,nstp,itrc)=t(i,js-1,k,nstp,itrc)*              &
     &                              rmask(i,js-1)
# endif
            END DO
          END IF
!
!  Parent grid at child NORTHERN EDGE
!
          js=hindcv(4)+1
          IF (((Istr.le.i).and.(i.le.Iend)).and.                        &
     &        ((Jstr.le.js).and.(js.le.Jend)).and.                      &
     &        (NINT(hindcj_v(i,js)).eq.(Mm(ngc)+1))) THEN
!
!    compute vertical sum of parent flux from last time.
!
            cff1=FEp_north(i,1,itrc)
            DO k=2,N(ng)
              cff1=cff1+FEp_north(i,k,itrc)
            END DO
!
!    Use cumulated child flux FXc and parent flux FXp from last time step
!    to compute correction for parent.  ic is the index of the child grid.
!
            ic=NINT(hindci_v(i,js))
            ref=(ref_scale-1)/2
            cff3=0.0_r8
            DO k=-ref,ref
              cff3=cff3+FEc_north_all(ic+k)
            END DO
!
!    correct parent flux.
!
            cff=pm(i,js)*pn(i,js)
            DO k=1,N(ng)
              cff2=(cff3-cff1)*FEp_north(i,k,itrc)/(cff1+eps)
# ifdef REF_DIFFUSION
              t(i,js,k,nstp,itrc)=t(i,js,k,nstp,itrc)*Hz(i,js,k)+       &
     &                            cff2*cff
# else
              t(i,js,k,nstp,itrc)=t(i,js,k,nstp,itrc)+                  &
     &                            cff2*cff*oHz(i,js,k)
#  if defined TCLM_NUDGING && defined TCLIMATOLOGY
!
!  Nudge towards tracer climatology.
!
              t(i,js,k,nstp,itrc)=t(i,js,k,nstp,itrc)+                  &
     &                              dt(ng)*Tnudgcof(i,j,itrc)*          &
     &                              (tclm(i,j,k,itrc)-                  &
     &                               t(i,js,k,nstp,itrc))
#  endif
# endif
# ifdef MASKING
              t(i,js,k,nstp,itrc)=t(i,js,k,nstp,itrc)*                  &
     &                            rmask(i,js)
# endif
            END DO
          END IF
        END DO
!
        deallocate (FXc_west_all)
        deallocate (FXc_east_all)
        deallocate (FEc_south_all)
        deallocate (FEc_north_all)
      END DO T_LOOP
#ifdef REF_DIFFUSION
!-----------------------------------------------------------------------
!  Time-step vertical diffusion term.
!-----------------------------------------------------------------------
!
!  Compute off-diagonal coefficients FC [lambda*dt*Akt/Hz] for the
!  implicit vertical diffusion terms at future time step, located
!  at horizontal RHO-points and vertical W-points.
!  Also set FC at the top and bottom levels.
!
      DO itrc=1,NT(ng)
        ltrc=MIN(NAT,itrc)
        js=hindcu(3)
        je=hindcu(4)
        DO j=js,je
!
!  Parent grid at child WESTERN EDGE
!
          is=hindcu(1)-1
          IF (((Istr.le.is-1).and.(is-1.le.Iend)).and.                  &
     &        ((Jstr.le.j).and.(j.le.Jend))) THEN
            cff=-dt(ng)*lambda
            DO k=1,N(ng)-1
              cff1=1.0_r8/(z_r(is-1,j,k+1)-z_r(is-1,j,k))
              FC(k)=cff*cff1*Akt(is-1,j,k,ltrc)
            END DO
            FC(0)=0.0_r8
            FC(N(ng))=0.0_r8
!
!  Compute diagonal matrix coefficients BC and load right-hand-side
!  terms for the tracer equation into DC.
!
            DO k=1,N(ng)
              BC(k)=Hz(is-1,j,k)-FC(k)-FC(k-1)
              DC(k)=t(is-1,j,k,nstp,itrc)
            END DO
!
!  Solve the tridiagonal system.
!
            cff=1.0_r8/BC(1)
            CF(1)=cff*FC(1)
            DC(1)=cff*DC(1)
            DO k=2,N(ng)-1
              cff=1.0_r8/(BC(k)-FC(k-1)*CF(k-1))
              CF(k)=cff*FC(k)
              DC(k)=cff*(DC(k)-FC(k-1)*DC(k-1))
            END DO
!
!  Compute new solution by back substitution.
!
            DC(N(ng))=(DC(N(ng))-FC(N(ng)-1)*DC(N(ng)-1))/     &
     &                (BC(N(ng))-FC(N(ng)-1)*CF(N(ng)-1))
            t(is-1,j,N(ng),nstp,itrc)=DC(N(ng))
            DO k=N(ng)-1,1,-1
              DC(k)=DC(k)-CF(k)*DC(k+1)
              t(is-1,j,k,nstp,itrc)=DC(k)
            END DO
          END IF
!
!  Parent grid at child EASTERN EDGE
!
          is=hindcu(2)+1
          IF (((Istr.le.is).and.(is.le.Iend)).and.                      &
     &        ((Jstr.le.j).and.(j.le.Jend))) THEN
            cff=-dt(ng)*lambda
            DO k=1,N(ng)-1
              cff1=1.0_r8/(z_r(is,j,k+1)-z_r(is,j,k))
              FC(k)=cff*cff1*Akt(is,j,k,ltrc)
            END DO
            FC(0)=0.0_r8
            FC(N(ng))=0.0_r8
!
!  Compute diagonal matrix coefficients BC and load right-hand-side
!  terms for the tracer equation into DC.
!
            DO k=1,N(ng)
              BC(k)=Hz(is,j,k)-FC(k)-FC(k-1)
              DC(k)=t(is,j,k,nstp,itrc)
            END DO
!
!  Solve the tridiagonal system.
!
            cff=1.0_r8/BC(1)
            CF(1)=cff*FC(1)
            DC(1)=cff*DC(1)
            DO k=2,N(ng)-1
              cff=1.0_r8/(BC(k)-FC(k-1)*CF(k-1))
              CF(k)=cff*FC(k)
              DC(k)=cff*(DC(k)-FC(k-1)*DC(k-1))
            END DO
!
!  Compute new solution by back substitution.
!
            DC(N(ng))=(DC(N(ng))-FC(N(ng)-1)*DC(N(ng)-1))/     &
     &                (BC(N(ng))-FC(N(ng)-1)*CF(N(ng)-1))
            t(is,j,N(ng),nstp,itrc)=DC(N(ng))
            DO k=N(ng)-1,1,-1
              DC(k)=DC(k)-CF(k)*DC(k+1)
              t(is,j,k,nstp,itrc)=DC(k)
            END DO
          END IF
        END DO
!
!  southern and northern child edges.
!
        is=hindcv(1)
        ie=hindcv(2)
        DO i=is,ie
!
!  Parent grid at child SOUTHERN EDGE
!
          js=hindcv(3)-1
          IF (((Istr.le.i).and.(i.le.Iend)).and.                        &
     &        ((Jstr.le.js-1).and.(js-1.le.Jend))) THEN
            cff=-dt(ng)*lambda
            DO k=1,N(ng)-1
              cff1=1.0_r8/(z_r(i,js-1,k+1)-z_r(i,js-1,k))
              FC(k)=cff*cff1*Akt(i,js-1,k,ltrc)
            END DO
            FC(0)=0.0_r8
            FC(N(ng))=0.0_r8
!
!  Compute diagonal matrix coefficients BC and load right-hand-side
!  terms for the tracer equation into DC.
!
            DO k=1,N(ng)
              BC(k)=Hz(i,js-1,k)-FC(k)-FC(k-1)
              DC(k)=t(i,js-1,k,nstp,itrc)
            END DO
!
!  Solve the tridiagonal system.
!
            cff=1.0_r8/BC(1)
            CF(1)=cff*FC(1)
            DC(1)=cff*DC(1)
            DO k=2,N(ng)-1
              cff=1.0_r8/(BC(k)-FC(k-1)*CF(k-1))
              CF(k)=cff*FC(k)
              DC(k)=cff*(DC(k)-FC(k-1)*DC(k-1))
            END DO
!
!  Compute new solution by back substitution.
!
            DC(N(ng))=(DC(N(ng))-FC(N(ng)-1)*DC(N(ng)-1))/     &
     &                (BC(N(ng))-FC(N(ng)-1)*CF(N(ng)-1))
            t(i,js-1,N(ng),nstp,itrc)=DC(N(ng))
            DO k=N(ng)-1,1,-1
              DC(k)=DC(k)-CF(k)*DC(k+1)
              t(i,js-1,k,nstp,itrc)=DC(k)
            END DO
          END IF
!
!  Parent grid at child NORTHERN EDGE
!
          js=hindcv(4)+1
          IF (((Istr.le.i).and.(i.le.Iend)).and.                        &
     &        ((Jstr.le.js).and.(js.le.Jend))) THEN
            cff=-dt(ng)*lambda
            DO k=1,N(ng)-1
              cff1=1.0_r8/(z_r(i,js,k+1)-z_r(i,js,k))
              FC(k)=cff*cff1*Akt(i,js,k,ltrc)
            END DO
            FC(0)=0.0_r8
            FC(N(ng))=0.0_r8
!
!  Compute diagonal matrix coefficients BC and load right-hand-side
!  terms for the tracer equation into DC.
!
            DO k=1,N(ng)
              BC(k)=Hz(i,js,k)-FC(k)-FC(k-1)
              DC(k)=t(i,js,k,nstp,itrc)
            END DO
!
!  Solve the tridiagonal system.
!
            cff=1.0_r8/BC(1)
            CF(1)=cff*FC(1)
            DC(1)=cff*DC(1)
            DO k=2,N(ng)-1
              cff=1.0_r8/(BC(k)-FC(k-1)*CF(k-1))
              CF(k)=cff*FC(k)
              DC(k)=cff*(DC(k)-FC(k-1)*DC(k-1))
            END DO
!
!  Compute new solution by back substitution.
!
            DC(N(ng))=(DC(N(ng))-FC(N(ng)-1)*DC(N(ng)-1))/     &
     &                (BC(N(ng))-FC(N(ng)-1)*CF(N(ng)-1))
            t(i,js,N(ng),nstp,itrc)=DC(N(ng))
            DO k=N(ng)-1,1,-1
              DC(k)=DC(k)-CF(k)*DC(k+1)
              t(i,js,k,nstp,itrc)=DC(k)
            END DO
          END IF
        END DO
      END DO
# endif

!-----------------------------------------------------------------------
!  Apply some diagnostics and exchanges.
!-----------------------------------------------------------------------
!
# ifdef DISTRIBUTE
      CALL mp_exchange4d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 1, N(ng), 1, NT(ng),      &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    t(:,:,:,nstp,:))
# endif

      RETURN
      END SUBROUTINE step3dref_t_tile
#endif
      END MODULE step3dref_t_mod
