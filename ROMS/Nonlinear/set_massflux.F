#include "cppdefs.h"
      MODULE set_massflux_mod

#ifdef SOLVE3D
# ifdef EW_PERIODIC
#  define IU_RANGE Istr,Iend
#  define IV_RANGE Istr,Iend
# elif defined COMPOSED_GRID || defined REFINED_GRID
#  define IU_RANGE IstrTU+1,IendT
#  define IV_RANGE IstrT,IendT
# else
#  define IU_RANGE Istr,IendR
#  define IV_RANGE IstrR,IendR
# endif

# ifdef NS_PERIODIC
#  define JU_RANGE Jstr,Jend
#  define JV_RANGE Jstr,Jend
# elif defined COMPOSED_GRID || defined REFINED_GRID
#  define JU_RANGE JstrT,JendT
#  define JV_RANGE JstrTV+1,JendT
# else
#  define JU_RANGE JstrR,JendR
#  define JV_RANGE Jstr,JendR
# endif
!
!svn $Id: set_massflux.F 732 2008-09-07 01:55:51Z jcwarner $
!=======================================================================
!  Copyright (c) 2002-2010 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                           Hernan G. Arango   !
!========================================== Alexander F. Shchepetkin ===
!                                                                      !
!  This routine computes horizontal mass fluxes, Hz*u/n and Hz*v/m.    !
!                                                                      !
!=======================================================================
!
!
      implicit none

      PRIVATE
      PUBLIC  :: set_massflux
# ifdef ADJOINT
      PUBLIC  :: reset_massflux
# endif

      CONTAINS
!
!***********************************************************************
      SUBROUTINE set_massflux (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
      USE mod_ocean
      USE mod_stepping
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
# include "tile.h"
!
# ifdef PROFILE
      CALL wclock_on (ng, iNLM, 12)
# endif
      CALL set_massflux_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        IminS, ImaxS, JminS, JmaxS,               &
     &                        nrhs(ng),                                 &
     &                        OCEAN(ng) % u,                            &
     &                        OCEAN(ng) % v,                            &
# ifdef WEC
     &                        OCEAN(ng) % u_stokes,                     &
     &                        OCEAN(ng) % v_stokes,                     &
# endif
     &                        GRID(ng) % Hz,                            &
     &                        GRID(ng) % om_v,                          &
     &                        GRID(ng) % on_u,                          &
     &                        GRID(ng) % Huon,                          &
     &                        GRID(ng) % Hvom)
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 12)
# endif
      RETURN
      END SUBROUTINE set_massflux
!
!***********************************************************************
      SUBROUTINE set_massflux_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              IminS, ImaxS, JminS, JmaxS,         &
     &                              nrhs,                               &
     &                              u, v,                               &
# ifdef WEC
     &                              u_stokes, v_stokes,                 &
# endif
     &                              Hz, om_v, on_u,                     &
     &                              Huon, Hvom)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars

# if defined EW_PERIODIC || defined NS_PERIODIC
!
      USE exchange_3d_mod
# endif
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange3d
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nrhs
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: v(LBi:,LBj:,:,:)
#  ifdef WEC
      real(r8), intent(in) :: u_stokes(LBi:,LBj:,:)
      real(r8), intent(in) :: v_stokes(LBi:,LBj:,:)
#  endif
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: om_v(LBi:,LBj:)
      real(r8), intent(in) :: on_u(LBi:,LBj:)

      real(r8), intent(out) :: Huon(LBi:,LBj:,:)
      real(r8), intent(out) :: Hvom(LBi:,LBj:,:)
# else
      real(r8), intent(in) :: u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: v(LBi:UBi,LBj:UBj,N(ng),2)
#  ifdef WEC
      real(r8), intent(in) :: u_stokes(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: v_stokes(LBi:UBi,LBj:UBj,N(ng))
#  endif
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: om_v(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: on_u(LBi:UBi,LBj:UBj)

      real(r8), intent(out) :: Huon(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(out) :: Hvom(LBi:UBi,LBj:UBj,N(ng))
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
      integer :: i, j, k

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Compute horizontal mass fluxes, Hz*u/n and Hz*v/m.
!-----------------------------------------------------------------------
!
      DO k=1,N(ng)
        DO j=JU_RANGE
          DO i=IU_RANGE
            Huon(i,j,k)=0.5_r8*(Hz(i,j,k)+Hz(i-1,j,k))*u(i,j,k,nrhs)*   &
     &                  on_u(i,j)
# ifdef WEC
            Huon(i,j,k)=Huon(i,j,k)+                                    &
     &                  0.5_r8*(Hz(i,j,k)+Hz(i-1,j,k))*                 &
     &                  u_stokes(i,j,k)*on_u(i,j)
# endif
          END DO
        END DO
        DO j=JV_RANGE
          DO i=IV_RANGE
            Hvom(i,j,k)=0.5_r8*(Hz(i,j,k)+Hz(i,j-1,k))*v(i,j,k,nrhs)*   &
     &                  om_v(i,j)
# ifdef WEC
            Hvom(i,j,k)=Hvom(i,j,k)+                                    &
     &                  0.5_r8*(Hz(i,j,k)+Hz(i,j-1,k))*                 &
     &                  v_stokes(i,j,k)*om_v(i,j)
# endif
          END DO
        END DO
      END DO
# undef IU_RANGE
# undef JU_RANGE
# undef IV_RANGE
# undef JV_RANGE
# if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
!
!  Exchange boundary information.
!
#  if defined EW_PERIODIC || defined NS_PERIODIC
      CALL exchange_u3d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        Huon)
      CALL exchange_v3d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        Hvom)
#  endif
#  ifdef DISTRIBUTE
      CALL mp_exchange3d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    Huon, Hvom)
#  endif
# endif
      RETURN
      END SUBROUTINE set_massflux_tile

# ifdef ADJOINT
!
!***********************************************************************
      SUBROUTINE reset_massflux (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_coupling
      USE mod_grid
      USE mod_ocean
      USE mod_stepping
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
#  include "tile.h"
!
#  ifdef PROFILE
      CALL wclock_on (ng, iADM, 12)
#  endif
      CALL reset_massflux_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          IminS, ImaxS, JminS, JmaxS,             &
     &                          nnew(ng),                               &
     &                          COUPLING(ng) % DU_avg2,                 &
     &                          COUPLING(ng) % DV_avg2,                 &
     &                          OCEAN(ng) % u,                          &
     &                          OCEAN(ng) % v,                          &
#  ifdef WEC
     &                          OCEAN(ng) % u_stokes,                   &
     &                          OCEAN(ng) % v_stokes,                   &
#  endif
     &                          GRID(ng) % Hz,                          &
     &                          GRID(ng) % om_v,                        &
     &                          GRID(ng) % on_u,                        &
     &                          GRID(ng) % Huon,                        &
     &                          GRID(ng) % Hvom)
#  ifdef PROFILE
      CALL wclock_off (ng, iADM, 12)
#  endif
      RETURN
      END SUBROUTINE reset_massflux
!
!***********************************************************************
      SUBROUTINE reset_massflux_tile (ng, tile,                         &
     &                                LBi, UBi, LBj, UBj,               &
     &                                IminS, ImaxS, JminS, JmaxS,       &
     &                                nnew,                             &
     &                                DU_avg2, DV_avg2,                 &
     &                                u, v,                             &
#  ifdef WEC
     &                                u_stokes, v_stokes,               &
#  endif
     &                                Hz, om_v, on_u,                   &
     &                                Huon, Hvom)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
#  if defined EW_PERIODIC || defined NS_PERIODIC
      USE exchange_3d_mod
#  endif
#  ifdef DISTRIBUTE
      USE mp_exchange_mod
#  endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS

      integer, intent(in) :: nnew
!
#  ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: DU_avg2(LBi:,LBj:)
      real(r8), intent(in) :: DV_avg2(LBi:,LBj:)
      real(r8), intent(in) :: u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: v(LBi:,LBj:,:,:)
#  ifdef WEC
      real(r8), intent(in) :: u_stokes(LBi:,LBj:,:)
      real(r8), intent(in) :: v_stokes(LBi:,LBj:,:)
#  endif
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: om_v(LBi:,LBj:)
      real(r8), intent(in) :: on_u(LBi:,LBj:)

      real(r8), intent(inout) :: Huon(LBi:,LBj:,:)
      real(r8), intent(inout) :: Hvom(LBi:,LBj:,:)
#  else
      real(r8), intent(in) :: DU_avg2(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: DV_avg2(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: v(LBi:UBi,LBj:UBj,N(ng),2)
#  ifdef WEC
      real(r8), intent(in) :: u_stokes(LBi:UBi,LBj:UBj,N(Ng))
      real(r8), intent(in) :: v_stokes(LBi:UBi,LBj:UBj,N(ng))
#  endif
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: om_v(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: on_u(LBi:UBi,LBj:UBj)

      real(r8), intent(inout) :: Huon(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: Hvom(LBi:UBi,LBj:UBj,N(ng))
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
      integer :: i, j, k

      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: DC
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: FC

#  include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Compute intermediate values of mass fluxes Huon and Hvom used by the
!  adjoint model.  The original values can be reinstated by calling
!  "set_massflux" after "ad_omega".
!-----------------------------------------------------------------------

#  ifdef EW_PERIODIC
#   define IU_RANGE Istr,Iend
#   define IV_RANGE Istr,Iend
# elif defined COMPOSED_GRID || defined REFINED_GRID
#   define IU_RANGE IstrT+1,IendT
#   define IV_RANGE IstrT,IendT
#  else
#   define IU_RANGE Istr,IendR
#   define IV_RANGE IstrR,IendR
#  endif
#  ifdef NS_PERIODIC
#   define JU_RANGE Jstr,Jend
#   define JV_RANGE Jstr,Jend
# elif defined COMPOSED_GRID || defined REFINED_GRID
#   define JU_RANGE JstrT,JendT
#   define JV_RANGE JstrT+1,JendT
#  else
#   define JU_RANGE JstrR,JendR
#   define JV_RANGE Jstr,JendR
#  endif
!
!  Compute mass flux, Hz*u/n.
!
      DO j=JU_RANGE
        DO i=IU_RANGE
          DC(i,0)=0.0_r8
          FC(i,0)=0.0_r8
        END DO
        DO k=1,N(ng)
          DO i=IU_RANGE
            DC(i,k)=0.5_r8*(Hz(i,j,k)+Hz(i-1,j,k))*on_u(i,j)
            DC(i,0)=DC(i,0)+DC(i,k)
          END DO
        END DO
        DO k=N(ng),1,-1
          DO i=IU_RANGE
            Huon(i,j,k)=0.5_r8*(Huon(i,j,k)+u(i,j,k,nnew)*DC(i,k))
# ifdef WEC
            Huon(i,j,k)=Huon(i,j,k)+0.5_r8*u_stokes(i,j,k)*DC(i,k)
# endif
            FC(i,0)=FC(i,0)+Huon(i,j,k)
          END DO
        END DO
!
!  Replace with correct vertical mean, DU_avg2.
!
        DO i=IU_RANGE
          DC(i,0)=1.0_r8/DC(i,0)
          FC(i,0)=DC(i,0)*(FC(i,0)-DU_avg2(i,j))
        END DO
        DO k=1,N(ng)
          DO i=IU_RANGE
            Huon(i,j,k)=Huon(i,j,k)-DC(i,k)*FC(i,0)
          END DO
        END DO
!
!  Compute mass flux, Hz*v/m.
!
        IF (j.ge.Jstr) THEN
          DO i=IV_RANGE
            DC(i,0)=0.0_r8
            FC(i,0)=0.0_r8
          END DO
          DO k=1,N(ng)
            DO i=IV_RANGE
              DC(i,k)=0.5_r8*(Hz(i,j,k)+Hz(i,j-1,k))*om_v(i,j)
              DC(i,0)=DC(i,0)+DC(i,k)
            END DO
          END DO
          DO k=N(ng),1,-1
            DO i=IV_RANGE
              Hvom(i,j,k)=0.5_r8*(Hvom(i,j,k)+v(i,j,k,nnew)*DC(i,k))
# ifdef WEC
              Hvom(i,j,k)=Hvom(i,j,k)+0.5_r8*v_stokes(i,j,k)*DC(i,k)
# endif
              FC(i,0)=FC(i,0)+Hvom(i,j,k)
            END DO
          END DO
!
!  Replace with correct vertical mean, DV_avg2.
!
          DO i=IV_RANGE
            DC(i,0)=1.0_r8/DC(i,0)
            FC(i,0)=DC(i,0)*(FC(i,0)-DV_avg2(i,j))
          END DO
          DO k=1,N(ng)
            DO i=IV_RANGE
              Hvom(i,j,k)=Hvom(i,j,k)-DC(i,k)*FC(i,0)
            END DO
          END DO
        ENDIF
      END DO
#  undef IU_RANGE
#  undef JU_RANGE
#  undef IV_RANGE
#  undef JV_RANGE
#  if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
!
!  Exchange boundary information.
!
#   if defined EW_PERIODIC || defined NS_PERIODIC
      CALL exchange_u3d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        Huon)
      CALL exchange_v3d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        Hvom)
#   endif
#   ifdef DISTRIBUTE
      CALL mp_exchange3d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    Huon, Hvom)
#   endif
#  endif
      RETURN
      END SUBROUTINE reset_massflux_tile
# endif
#endif
      END MODULE set_massflux_mod
