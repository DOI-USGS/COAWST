#include "cppdefs.h"
      MODULE curr_inw_mod
#if defined INWAVE_COUPLING
# if defined DOPPLER
!
!svn $Id: celer_inw.F 732 2008-09-07 01:55:51Z jcwarner $

!======================================================================!
!                                                                      !
!  This routine computes the currents affecting the wave field         !
!  Needed to solve the action density equations considering            !
!  wave-current interaction effects                                    !
!                                                                      !
!======================================================================!
!
      implicit none
      PRIVATE
      PUBLIC  :: curr_inw
      CONTAINS
!
!***********************************************************************
      SUBROUTINE curr_inw (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_ocean
      USE mod_stepping
      USE mod_inwave_vars
      USE mod_inwave_params
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
#  include "tile.h"
!
!#  ifdef PROFILE
!      CALL wclock_on (ng, iNLM, 35)
!#  endif

      CALL curr_inw_tile(ng, tile,                                      &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   nrhs(ng), nstp(ng), nnew(ng),                  &
#  ifdef SOLVE3D
     &                   OCEAN(ng) % u,                                 &
     &                   OCEAN(ng) % v,                                 &
#  else
     &                   OCEAN(ng) % ubar,                              &
     &                   OCEAN(ng) % vbar,                              &
#  endif
     &                   WAVEP(ng) % kwc,                               &
     &                   WAVEP(ng) % h_tot,                             &
     &                   WAVEP(ng) % u_rho,                             &
     &                   WAVEP(ng) % v_rho)

!#  ifdef PROFILE
!      CALL wclock_off (ng, iNLM, 35)
!#  endif

      RETURN
      END SUBROUTINE curr_inw
!
!***********************************************************************
      SUBROUTINE curr_inw_tile(ng, tile,                                &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         nrhs, nstp, nnew,                        &
#  ifdef SOLVE3D
     &                         u, v,                                    &
#  else
     &                         ubar, vbar,                              &
#  endif
     &                         kwc, h_tot, u_rho, v_rho)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
      USE mod_inwave_params

#  ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
#  endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: nrhs, nstp, nnew

#  ifdef ASSUMED_SHAPE
#   ifdef SOLVE3D
      real(r8), intent(in) :: u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: v(LBi:,LBj:,:,:)
#   else
      real(r8), intent(in) :: ubar(LBi:,LBj:,:)
      real(r8), intent(in) :: vbar(LBi:,LBj:,:)
#   endif
      real(r8), intent(in) :: kwc(LBi:,LBj:,:)
      real(r8), intent(in) :: h_tot(LBi:,LBj:)
      real(r8), intent(inout) :: u_rho(LBi:,LBj:)
      real(r8), intent(inout) :: v_rho(LBi:,LBj:)
#  else
#   ifdef SOLVE3D
      real(r8), intent(in) :: u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: v(LBi:UBi,LBj:UBj,N(ng),2)
#   else
      real(r8), intent(in) :: ubar(LBi:UBi,LBj:UBj,2)
      real(r8), intent(in) :: vbar(LBi:UBi,LBj:UBj,2)
#   endif
      real(r8), intent(in) :: kwc(LBi:UBi,LBj:UBj,ND(ng))
      real(r8), intent(in) :: h_tot(LBi:UBi,LBj:UBj,ND(ng))
      real(r8), intent(inout) :: u_rho(LBi:UBi,LBj:UBj,ND(ng))
      real(r8), intent(inout) :: v_rho(LBi:UBi,LBj:UBj,ND(ng))
#  endif
!
!  Local variable declarations.
!
#  ifdef DISTRIBUTE
#   ifdef EW_PERIODIC
      logical :: EWperiodic=.TRUE.
#   else
      logical :: EWperiodic=.FALSE.
#   endif
#   ifdef NS_PERIODIC
      logical :: NSperiodic=.TRUE.
#   else
      logical :: NSperiodic=.FALSE.
#   endif
#  endif
!
      integer :: i, j, k, d

      real(r8) :: twopi, otwopi
      real(r8) :: cff1, cff2
      real(r8) :: cff1_u, cff2_u, cff3_u
      real(r8) :: cff1_v, cff2_v, cff3_v
      real(r8) :: kwc_mean
!
#  include "set_bounds.h"
!
      twopi=2.0_r8*pi
      otwopi=1.0_r8/twopi

      DO j=Jstr,Jend
        DO i=Istr,Iend

! For the different directional components when there are currents the Tr
! is different and therefore also the kwc and the currents affecting the 
! waves. Here I am assuming that this difference is very small.
! We can change this is the future.

          kwc_mean=0.0_r8
          DO d=1,ND(ng)
            kwc_mean=kwc_mean+kwc(i,j,d)
          END DO
          kwc_mean=kwc_mean/ND(ng)
!
!======================================================================!
!     Compute u and v currents affecting the wave field at rho points  !
!     these can be computed using Kirby and Chen (1989) or assume that !
!     the currents are the superficial currents; in the 2D case those  !
!     computed with the vertically averaged currents                   !
!======================================================================!
!
#  ifdef SOLVE3D
#   ifdef UV_KIRBY
          cff1=2.0_r8*kwc_mean*h_tot(i,j)
          IF (cff1.lt.700.0_r8) THEN
            cff1=2.0_r8*kwc_mean
          ELSE
            cff1=700.0_r8/h_tot(i,j)
          ENDIF
          cff2=0.0_r8
          cff1_u=0.0_r8
          cff2_u=0.0_r8
          cff3_u=0.0_r8
          cff1_v=0.0_r8
          cff2_v=0.0_r8
          cff3_v=0.0_r8
          DO k=1,N(ng)
            cff2=COSH(cff1*(GRID(ng)%h_tot(i,j)+GRID(ng)%z_r(i,j,k)))*&
     &           GRID(ng)%Hz(i,j,k)
            cff1_u=0.5_r8*(u(i,j,k,nnew)+u(i+1,j,k,nnew))
            cff2_u=cff2_u+cff2*cff1_u
            cff3_u=cff3_u+cff2
            cff1_v=0.5_r8*(v(i,j,k,nnew)+v(i,j+1,k,nnew))
            cff2_v=cff2_v+cff2*cff1_v
            cff3_v=cff3_v+cff2
          END DO
          u_rho(i,j)=cff2_u/cff3_u
          v_rho(i,j)=cff2_v/cff3_v
#   else
          u_rho(i,j)=0.5_r8*(u(i,j,N(ng),nnew)+u(i+1,j,N(ng),nnew))
          v_rho(i,j)=0.5_r8*(v(i,j,N(ng),nnew)+v(i,j+1,N(ng),nnew))
#   endif
#  else
          u_rho(i,j)=0.5_r8*(ubar(i,j,nnew)+ubar(i+1,j,nnew))
          v_rho(i,j)=0.5_r8*(vbar(i,j,nnew)+vbar(i,j+1,nnew))
#  endif
        END DO
      END DO
      
      

#  ifdef DISTRIBUTE
!
!   Exchange boundary data.
!
      CALL mp_exchange2d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    u_rho(:,:))

      CALL mp_exchange2d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    v_rho(:,:))
!
#  endif

      RETURN
      END SUBROUTINE curr_inw_tile
# endif
#endif
      END MODULE curr_inw_mod
