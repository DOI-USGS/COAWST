#include "cppdefs.h"
      MODULE frc_weak_mod

#ifdef WEAK_CONSTRAINT
!
!svn $Id: frc_weak.F 795 2016-05-11 01:42:43Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  These routines are used during the computation of the weak          !
!  constraint forcing.                                                 !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC :: frc_ADgather
      PUBLIC :: frc_clear

      CONTAINS

      SUBROUTINE frc_ADgather (ng, tile)
!
!=======================================================================
!                                                                      !
!  This subroutine is the adjoint of the  weak constraint forcing      !
!  interpolation between snapshots used in the tangent linear and      !
!  representer models.                                                 !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng        Nested grid number.                                    !
!     tile      Domain partition.                                      !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_forces
      USE mod_ocean
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
      CALL wclock_on (ng, iADM, 7)
# endif
      CALL frc_ADgather_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        IminS, ImaxS, JminS, JmaxS,               &
     &                        OCEAN(ng) % ad_zeta,                      &
     &                        OCEAN(ng) % ad_zeta_sol,                  &
# ifdef SOLVE3D
     &                        OCEAN(ng) % ad_u,                         &
     &                        OCEAN(ng) % ad_v,                         &
     &                        OCEAN(ng) % ad_t,                         &
# else
     &                        OCEAN(ng) % ad_ubar,                      &
     &                        OCEAN(ng) % ad_vbar,                      &
     &                        OCEAN(ng) % ad_ubar_sol,                  &
     &                        OCEAN(ng) % ad_vbar_sol,                  &
# endif
     &                        OCEAN(ng) % f_zetaG,                      &
# ifdef SOLVE3D
     &                        OCEAN(ng) % f_uG,                         &
     &                        OCEAN(ng) % f_vG,                         &
     &                        OCEAN(ng) % f_tG)
# else
     &                        OCEAN(ng) % f_ubarG,                      &
     &                        OCEAN(ng) % f_vbarG)
# endif
# ifdef PROFILE
      CALL wclock_off (ng, iADM, 7)
# endif
      RETURN
      END SUBROUTINE frc_ADgather
!
!***********************************************************************
      SUBROUTINE frc_ADgather_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              IminS, ImaxS, JminS, JmaxS,         &
     &                              ad_zeta, ad_zeta_sol,               &
# ifdef SOLVE3D
     &                              ad_u, ad_v, ad_t,                   &
# else
     &                              ad_ubar, ad_vbar,                   &
     &                              ad_ubar_sol, ad_vbar_sol,           &
# endif
     &                              f_zetaG,                            &
# ifdef SOLVE3D
     &                              f_uG, f_vG, f_tG)
# else
     &                              f_ubarG, f_vbarG)
# endif
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
      USE mod_stepping
      USE mod_fourdvar
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: ad_zeta(LBi:,LBj:,:)
      real(r8), intent(in) :: ad_zeta_sol(LBi:,LBj:)
      real(r8), intent(inout) :: f_zetaG(LBi:,LBj:,:)
#  ifdef SOLVE3D
      real(r8), intent(in) :: ad_u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: ad_v(LBi:,LBj:,:,:)
      real(r8), intent(in) :: ad_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: f_uG(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: f_vG(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: f_tG(LBi:,LBj:,:,:,:)
#  else
      real(r8), intent(in) :: ad_ubar(LBi:,LBj:,:)
      real(r8), intent(in) :: ad_vbar(LBi:,LBj:,:)
      real(r8), intent(in) :: ad_ubar_sol(LBi:,LBj:)
      real(r8), intent(inout) :: ad_vbar_sol(LBi:,LBj:)
      real(r8), intent(inout) :: f_ubarG(LBi:,LBj:,:)
      real(r8), intent(inout) :: f_vbarG(LBi:,LBj:,:)
#  endif
# else
      real(r8), intent(in) :: ad_zeta(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: ad_zeta_sol(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: f_zetaG(LBi:UBi,LBj:UBj,2)
#  ifdef SOLVE3D
      real(r8), intent(in) :: ad_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: ad_v(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: ad_t(LBi:UBi,LBj:UBj,N(ng),2,NT(ng))
      real(r8), intent(inout) :: f_uG(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: f_vG(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: f_tG(LBi:UBi,LBj:UBj,N(ng),2,NT(ng))
#  else
      real(r8), intent(in) :: ad_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: ad_vbar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: ad_ubar_sol(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: ad_vbar_sol(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: f_ubarG(LBi:UBi,LBj:UBj,2)
      real(r8), intent(inout) :: f_vbarG(LBi:UBi,LBj:UBj,2)
#  endif
# endif
!
!  Local variable declarations.
!
      integer :: i, it1, it2, j, k, kout
# ifdef SOLVE3D
      integer :: itrc, nout
# endif
      real(r8) :: fac, fac1, fac2, time1, time2

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Compute weak-constraint forcing terms used by the tangent linear
!  and representer models. Perform the adjoint of the interpolation
!  between snapshots.
!-----------------------------------------------------------------------
!
!  Set time records and interpolation factor, if any.
!
      it1=MAX(0,(iic(ng)-1)/nADJ(ng))+1
      it2=it1+1
      time1=dstart*day2sec+(it1-1)*nADJ(ng)*dt(ng)
      time2=dstart*day2sec+(it2-1)*nADJ(ng)*dt(ng)
      fac1=time2-time(ng)
      fac2=time(ng)-time1
      fac=1.0_r8/(fac1+fac2)
      fac1=fac*fac1
      fac2=fac*fac2
!
!  Set weak-constraint force time.
!
      ForceTime(ng)=time2
!
!  Determine time index of adjoint variables to process.
!
      kout=kstp(ng)

# ifdef SOLVE3D
      IF (iic(ng).ne.ntend(ng)) THEN
        nout=nnew(ng)
      ELSE
        nout=nstp(ng)
      END IF
# endif
!
!  Clear forcing arrays on first call to this routine.
!
      IF (iic(ng).eq.ntstart(ng)-1) THEN
        DO j=JstrR,JendR
          DO i=IstrR,IendR
           f_zetaG(i,j,1)=0.0_r8
           f_zetaG(i,j,2)=0.0_r8
          END DO
        END DO
# ifndef SOLVE3D
        DO j=JstrR,JendR
          DO i=Istr,IendR
             f_ubarG(i,j,1)=0.0_r8
             f_ubarG(i,j,2)=0.0_r8
          END DO
        END DO
        DO j=Jstr,JendR
          DO i=IstrR,IendR
             f_vbarG(i,j,1)=0.0_r8
             f_vbarG(i,j,2)=0.0_r8
          END DO
        END DO
# else
        DO k=1,N(ng)
          DO j=JstrR,JendR
            DO i=Istr,IendR
             f_uG(i,j,k,1)=0.0_r8
             f_uG(i,j,k,2)=0.0_r8
            END DO
          END DO
        END DO
        DO k=1,N(ng)
          DO j=Jstr,JendR
            DO i=IstrR,IendR
             f_vG(i,j,k,1)=0.0_r8
             f_vG(i,j,k,2)=0.0_r8
            END DO
          END DO
        END DO
        DO itrc=1,NT(ng)
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
               f_tG(i,j,k,1,itrc)=0.0_r8
               f_tG(i,j,k,2,itrc)=0.0_r8
              END DO
            END DO
          END DO
        END DO
# endif
      END IF
!
!  Gather free-surface weak-constraint forcing terms.
!
      IF (LwrtState2d(ng)) THEN
        DO j=JstrR,JendR
          DO i=IstrR,IendR
           f_zetaG(i,j,1)=f_zetaG(i,j,1)+fac1*ad_zeta(i,j,kout)
           f_zetaG(i,j,2)=f_zetaG(i,j,2)+fac2*ad_zeta(i,j,kout)
          END DO
        END DO
      ELSE
        DO j=JstrR,JendR
          DO i=IstrR,IendR
           f_zetaG(i,j,1)=f_zetaG(i,j,1)+fac1*ad_zeta_sol(i,j)
           f_zetaG(i,j,2)=f_zetaG(i,j,2)+fac2*ad_zeta_sol(i,j)
          END DO
        END DO
      END IF

# ifndef SOLVE3D
!
!  Gather 2D-momentum weak-constraint forcing terms.
!
      IF (LwrtState2d(ng)) THEN
        DO j=JstrR,JendR
          DO i=Istr,IendR
             f_ubarG(i,j,1)=f_ubarG(i,j,1)+fac1*ad_ubar(i,j,kout)
             f_ubarG(i,j,2)=f_ubarG(i,j,2)+fac2*ad_ubar(i,j,kout)
          END DO
        END DO
        DO j=Jstr,JendR
          DO i=IstrR,IendR
             f_vbarG(i,j,1)=f_vbarG(i,j,1)+fac1*ad_vbar(i,j,kout)
             f_vbarG(i,j,2)=f_vbarG(i,j,2)+fac2*ad_vbar(i,j,kout)
          END DO
        END DO
      ELSE
        DO j=JstrR,JendR
          DO i=Istr,IendR
             f_ubarG(i,j,1)=f_ubarG(i,j,1)+fac1*ad_ubar_sol(i,j)
             f_ubarG(i,j,2)=f_ubarG(i,j,2)+fac2*ad_ubar_sol(i,j)
          END DO
        END DO
        DO j=Jstr,JendR
          DO i=IstrR,IendR
             f_vbarG(i,j,1)=f_vbarG(i,j,1)+fac1*ad_vbar_sol(i,j)
             f_vbarG(i,j,2)=f_vbarG(i,j,2)+fac2*ad_vbar_sol(i,j)
          END DO
        END DO
      END IF
# endif
# ifdef SOLVE3D
!
!  Gather 3D-momentum weak-constraint forcing terms.
!
      DO k=1,N(ng)
        DO j=JstrR,JendR
          DO i=Istr,IendR
           f_uG(i,j,k,1)=f_uG(i,j,k,1)+fac1*ad_u(i,j,k,nout)
           f_uG(i,j,k,2)=f_uG(i,j,k,2)+fac2*ad_u(i,j,k,nout)
          END DO
        END DO
      END DO
      DO k=1,N(ng)
        DO j=Jstr,JendR
          DO i=IstrR,IendR
           f_vG(i,j,k,1)=f_vG(i,j,k,1)+fac1*ad_v(i,j,k,nout)
           f_vG(i,j,k,2)=f_vG(i,j,k,2)+fac2*ad_v(i,j,k,nout)
          END DO
        END DO
      END DO
!
!  Gather tracer weak-constraint forcing terms.
!
      DO itrc=1,NT(ng)
        DO k=1,N(ng)
          DO j=JstrR,JendR
            DO i=IstrR,IendR
             f_tG(i,j,k,1,itrc)=f_tG(i,j,k,1,itrc)+                     &
     &                          fac1*ad_t(i,j,k,nout,itrc)
             f_tG(i,j,k,2,itrc)=f_tG(i,j,k,2,itrc)+                     &
     &                          fac2*ad_t(i,j,k,nout,itrc)
            END DO
          END DO
        END DO
      END DO
# endif

      RETURN
      END SUBROUTINE frc_ADgather_tile

      SUBROUTINE frc_clear (ng, tile)
!
!=======================================================================
!                                                                      !
!  This routine copy weak-constraint arrays (f_***G storage arrays)    !
!  index 1 into index 2 and then clear index 1.                        !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng        Nested grid number.                                    !
!     tile      Domain partition.                                      !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_forces
      USE mod_ocean
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
      CALL wclock_on (ng, iADM, 7)
# endif
      CALL frc_clear_tile (ng, tile,                                    &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
     &                     OCEAN(ng) % f_zetaG,                         &
# ifdef SOLVE3D
     &                     OCEAN(ng) % f_uG,                            &
     &                     OCEAN(ng) % f_vG,                            &
     &                     OCEAN(ng) % f_tG)
# else
     &                     OCEAN(ng) % f_ubarG,                         &
     &                     OCEAN(ng) % f_vbarG)
# endif
# ifdef PROFILE
      CALL wclock_off (ng, iADM, 7)
# endif
      RETURN
      END SUBROUTINE frc_clear
!
!***********************************************************************
      SUBROUTINE frc_clear_tile (ng, tile,                              &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           IminS, ImaxS, JminS, JmaxS,            &
     &                           f_zetaG,                               &
# ifdef SOLVE3D
     &                           f_uG, f_vG, f_tG)
# else
     &                           f_ubarG, f_vbarG)
# endif
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
      USE mod_stepping
      USE mod_fourdvar
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: f_zetaG(LBi:,LBj:,:)
#  ifdef SOLVE3D
      real(r8), intent(inout) :: f_uG(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: f_vG(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: f_tG(LBi:,LBj:,:,:,:)
#  else
      real(r8), intent(inout) :: f_ubarG(LBi:,LBj:,:)
      real(r8), intent(inout) :: f_vbarG(LBi:,LBj:,:)
#  endif
# else
      real(r8), intent(inout) :: f_zetaG(LBi:UBi,LBj:UBj,2)
#  ifdef SOLVE3D
      real(r8), intent(inout) :: f_uG(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: f_vG(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: f_tG(LBi:UBi,LBj:UBj,N(ng),2,NT(ng))
#  else
      real(r8), intent(inout) :: f_ubarG(LBi:UBi,LBj:UBj,2)
      real(r8), intent(inout) :: f_vbarG(LBi:UBi,LBj:UBj,2)
#  endif
# endif
!
!  Local variable declarations.
!
      integer :: i, it1, it2, j, k, kout, nout
# ifdef SOLVE3D
      integer :: itrc
# endif
      real(r8) :: fac, fac1, fac2, time1, time2

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Copy weak-constraint forcing arrays index 1 into index 2, and
!  clear index 1.
!-----------------------------------------------------------------------
!
!  Reset weak-constraing forcing time on last timestep.
!
      IF (iic(ng).eq.ntend(ng)) THEN
        ForceTime(ng)=dstart*day2sec
      END IF
!
!  Update free-surface weak-constraint forcing terms.
!
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          f_zetaG(i,j,2)=f_zetaG(i,j,1)
          f_zetaG(i,j,1)=0.0_r8
        END DO
      END DO

# ifndef SOLVE3D
!
!  Update 2D-momentum weak-constraint forcing terms.
!
      DO j=JstrR,JendR
        DO i=Istr,IendR
          f_ubarG(i,j,2)=f_ubarG(i,j,1)
          f_ubarG(i,j,1)=0.0_r8
        END DO
      END DO
      DO j=Jstr,JendR
        DO i=IstrR,IendR
          f_vbarG(i,j,2)=f_vbarG(i,j,1)
          f_vbarG(i,j,1)=0.0_r8
        END DO
      END DO
# endif
# ifdef SOLVE3D
!
!  Update 3D-momentum weak-constraint forcing terms.
!
      DO k=1,N(ng)
        DO j=JstrR,JendR
          DO i=Istr,IendR
           f_uG(i,j,k,2)=f_uG(i,j,k,1)
           f_uG(i,j,k,1)=0.0_r8
          END DO
        END DO
        DO j=Jstr,JendR
          DO i=IstrR,IendR
           f_vG(i,j,k,2)=f_vG(i,j,k,1)
           f_vG(i,j,k,1)=0.0_r8
          END DO
        END DO
      END DO
!
!  Update tracer weak-constraint forcing terms.
!
      DO itrc=1,NT(ng)
        DO k=1,N(ng)
          DO j=JstrR,JendR
            DO i=IstrR,IendR
             f_tG(i,j,k,2,itrc)=f_tG(i,j,k,1,itrc)
             f_tG(i,j,k,1,itrc)=0.0_r8
            END DO
          END DO
        END DO
      END DO
# endif

      RETURN
      END SUBROUTINE frc_clear_tile
#endif
      END MODULE frc_weak_mod
