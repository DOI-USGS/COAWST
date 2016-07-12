      SUBROUTINE ana_psource (ng, tile, model)
!
!! svn $Id: ana_psource.h 795 2016-05-11 01:42:43Z arango $
!!======================================================================
!! Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!!   Licensed under a MIT/X style license                              !
!!   See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This subroutine sets analytical tracer and mass point Sources       !
!  and/or Sinks.  River runoff can be consider as a point source.      !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_grid
      USE mod_ncparam
      USE mod_ocean
      USE mod_stepping
!
      integer, intent(in) :: ng, tile, model

#include "tile.h"
!
      CALL ana_psource_tile (ng, tile, model,                           &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       IminS, ImaxS, JminS, JmaxS,                &
     &                       nnew(ng), knew(ng),                        &
     &                       OCEAN(ng) % zeta,                          &
     &                       OCEAN(ng) % ubar,                          &
     &                       OCEAN(ng) % vbar,                          &
#ifdef SOLVE3D
     &                       OCEAN(ng) % u,                             &
     &                       OCEAN(ng) % v,                             &
     &                       GRID(ng) % z_w,                            &
#endif
     &                       GRID(ng) % h,                              &
     &                       GRID(ng) % on_u,                           &
     &                       GRID(ng) % om_v)
!
! Set analytical header file name used.
!
#ifdef DISTRIBUTE
      IF (Lanafile) THEN
#else
      IF (Lanafile.and.(tile.eq.0)) THEN
#endif
        ANANAME(20)=__FILE__
      END IF

      RETURN
      END SUBROUTINE ana_psource
!
!***********************************************************************
      SUBROUTINE ana_psource_tile (ng, tile, model,                     &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             IminS, ImaxS, JminS, JmaxS,          &
     &                             nnew, knew,                          &
     &                             zeta, ubar, vbar,                    &
#ifdef SOLVE3D
     &                             u, v, z_w,                           &
#endif
     &                             h, on_u, om_v)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_scalars
#ifdef SEDIMENT
      USE mod_sediment
#endif
      USE mod_sources

#ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_bcastf, mp_bcasti
      USE distribute_mod, ONLY : mp_collect, mp_reduce
#endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nnew, knew
!
#ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: zeta(LBi:,LBj:,:)
      real(r8), intent(in) :: ubar(LBi:,LBj:,:)
      real(r8), intent(in) :: vbar(LBi:,LBj:,:)
# ifdef SOLVE3D
      real(r8), intent(in) :: u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: v(LBi:,LBj:,:,:)
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
# endif
      real(r8), intent(in) :: h(LBi:,LBj:)
      real(r8), intent(in) :: on_u(LBi:,LBj:)
      real(r8), intent(in) :: om_v(LBi:,LBj:)
#else
      real(r8), intent(in) :: zeta(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: vbar(LBi:UBi,LBj:UBj,3)
# ifdef SOLVE3D
      real(r8), intent(in) :: u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: v(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: z_w(LBi:UBi,LBj:UBj,0:N(ng))
# endif
      real(r8), intent(in) :: h(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: on_u(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: om_v(LBi:UBi,LBj:UBj)
#endif
!
!  Local variable declarations.
!
      integer :: Npts, NSUB, is, i, j, k, ised

      real(r8) :: Pspv = 0.0_r8
      real(r8) :: cff, fac

#if defined DISTRIBUTE && defined SOLVE3D
      real(r8), dimension(Msrc(ng)*N(ng)) :: Pwrk
#endif
#if defined DISTRIBUTE
      real(r8), dimension(2) :: buffer

      character (len=3), dimension(2) :: io_handle
#endif

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  If initialization, set point Sources and/or Sinks locations.
!-----------------------------------------------------------------------
!
!! WARNING:  This routine is extremely difficult in parallel
!!           applications.  See examples for how to do this in
!!           source file:  "ROMS/Functionals/ana_psource.h".
!!
      IF (iic(ng).eq.ntstart(ng)) THEN
!
!  Set-up point Sources/Sink number (Nsrc), direction (Dsrc), I- and
!  J-grid locations (Isrc,Jsrc). Currently, the direction can be along
!  XI-direction (Dsrc = 0) or along ETA-direction (Dsrc > 0).  The
!  mass sources are located at U- or V-points so the grid locations
!  should range from 1 =< Isrc =< L  and  1 =< Jsrc =< M.
!
#if defined MY_APPLICATION
        IF (Master.and.DOMAIN(ng)%SouthWest_Test(tile)) THEN
          Nsrc(ng)=???
          SOURCES(ng)%Dsrc(:)=???
          SOURCES(ng)%Isrc(:)=???
          SOURCES(ng)%Jsrc(:)=???
        END IF
#else
        ana_psource.h: No values provided for Nsrc, Dsrc, Isrc, Jsrc.
#endif
#ifdef DISTRIBUTE
!
!  Broadcast point sources/sinks information to all nodes.
!
        CALL mp_bcasti (ng, iNLM, Nsrc(ng))
        CALL mp_bcasti (ng, iNLM, SOURCES(ng)%Isrc)
        CALL mp_bcasti (ng, iNLM, SOURCES(ng)%Jsrc)
        CALL mp_bcastf (ng, iNLM, SOURCES(ng)%Dsrc)
#endif
      END IF
!
!-----------------------------------------------------------------------
!  Set momentum point Sources and/or Sinks.
!-----------------------------------------------------------------------
!
      MOMENTUM : IF (LuvSrc(ng).or.LwSrc(ng)) THEN

#ifdef SOLVE3D
!
!  If appropriate, set-up nondimensional shape function to distribute
!  mass point sources/sinks vertically.  It must add to unity!!.
!
# ifdef DISTRIBUTE
        SOURCES(ng)%Qshape=Pspv
# endif
        Npts=Msrc(ng)*N(ng)

!$OMP BARRIER

# if defined MY_APPLICATION
!!
!!  Notice that in the simple statement below there is not need for
!!  distributed-memory communications since the computation below
!!  does not have a parallel tile dependency. All the nodes are
!!  computing this simple statement.  See original source routine
!!  "ROMS/Functionals/ana_psource.h" for examples that has parallel
!!  tile dependency.
!!
        IF (DOMAIN(ng)%NorthEast_Test(tile)) THEN
          DO k=1,N(ng)
            DO is=1,Nsrc(ng)
              SOURCES(ng)%Qshape(is,k)=1.0_r8/REAL(N(ng),r8)
            END DO
          END DO
        END IF
# else
        ana_psource.h: No values provided for Qshape.
# endif
#endif
!
!  Set-up vertically integrated mass transport (m3/s) of point
!  Sources/Sinks (positive in the positive U- or V-direction and
!  viceversa).
!
#ifdef DISTRIBUTE
        SOURCES(ng)%Qbar=Pspv
#endif

!$OMP BARRIER

#if defined MY_APPLICATION
!!
!!  Notice that in the simple statement below there is not need for
!!  distributed-memory communications since the computation below
!!  does not have a parallel tile dependency. All the nodes are
!!  computing this simple statement.  See original source routine
!!  "ROMS/Functional/ana_psource.h" for examples that has parallel
!!  tile dependency.
!!
        IF (DOMAIN(ng)%SouthWest_Test(tile)) THEN
          DO is=1,Nsrc(ng)
            SOURCES(ng)%Qbar(is)=???
          END DO
        END IF
#else
        ana_psource.h: No values provided for Qbar.
#endif

#ifdef SOLVE3D
!
!  Set-up mass transport profile (m3/s) of point Sources/Sinks.
!
!$OMP BARRIER

        IF (DOMAIN(ng)%NorthEast_Test(tile)) THEN
          DO k=1,N(ng)
            DO is=1,Nsrc(ng)
              SOURCES(ng)%Qsrc(is,k)=SOURCES(ng)%Qbar(is)*              &
     &                               SOURCES(ng)%Qshape(is,k)
            END DO
          END DO
        END IF
#endif
      END IF MOMENTUM

#ifdef SOLVE3D
!
!-----------------------------------------------------------------------
!  Set tracers point Sources and/or Sinks.
!-----------------------------------------------------------------------
!
      TRACERS : IF (ANY(LtracerSrc(:,ng))) THEN
        SOURCES(ng)%Tsrc=0.0_r8                      ! initialize
!
!  Set-up tracer (tracer units) point Sources/Sinks.
!
# if defined MY_APPLICATION
        IF (DOMAIN(ng)%NorthWest_Test(tile)) THEN
#  ifdef ONE_TRACER_SOURCE
          SOURCES(ng)%Tsrc(itemp)=???
          SOURCES(ng)%Tsrc(isalt)=???
#  elif defined TWO_D_TRACER_SOURCE
          DO is=1,Nsrc(ng)
            SOURCES(ng)%Tsrc(is,itemp)=???
            SOURCES(ng)%Tsrc(is,isalt)=???
          END DO
#  else
          DO k=1,N(ng)
            DO is=1,Nsrc(ng)
              SOURCES(ng)%Tsrc(is,k,itemp)=???
              SOURCES(ng)%Tsrc(is,k,isalt)=???
            END DO
          END DO
#  endif
        END IF
# else
        ana_psource.h: No values provided for Tsrc.
# endif
      END IF TRACERS
#endif

      RETURN
      END SUBROUTINE ana_psource_tile
