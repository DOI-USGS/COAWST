#include "cppdefs.h"
      MODULE back_cost_mod
#ifdef BACKGROUND
!
!svn $Id: back_cost.F 795 2016-05-11 01:42:43Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine computes the background cost function (Jb) as the      !
!  misfit (squared difference) between model and background state:     !
!                                                                      !
!    Jb = 1/2 * transpose(X - Xb) * B^(-1) * (X - Xb)                  !
!                                                                      !
!  where                                                               !
!                                                                      !
!    Xb : background or reference state (first guess)                  !
!    X  : nonlinear model initial state                                !
!    B  : background-error covariance                                  !
!                                                                      !
!  In incremental 4DVAR, the initial conditions estimate is done in    !
!  in minimization space by defining the following transformation:     !
!                                                                      !
!               deltaV = B^(-1/2) deltaX      (v-space)                !
!  or                                                                  !
!               deltaX = B^(1/2) deltaV       (x-space)                !
!                                                                      !
!  where                                                               !
!                                                                      !
!                    B = tanspose{B^(1/2)} B^(1/2)                     !
!                                                                      !
!  Then, the background cost function becomes:                         !
!                                                                      !
!    Jb = 1/2 * transpose{(X - Xb)} * B^(-1) * (X - Xb)     x-space    !
!                                                                      !
!  or                                                                  !
!                                                                      !
!    Jb = 1/2 * transpose(deltaX) * B^(-1) * deltaX         x-space    !
!                                                                      !
!  or                                                                  !
!                                                                      !
!    Jb = 1/2 * transpose(deltaV) * deltaV                  v-space    !
!                                                                      !
!  Therefore, in v-space the background cost function gradient is:     !
!                                                                      !
!    GRADv(Jb) = deltaV                                                !
!                                                                      !
!  Notice that initially, Jb is zero since the model is initialized    !
!  with the background state (X=Xb) and the minimization increment,    !
!  deltaX is zero.                                                     !
!                                                                      !
!                                                                      !
!  WARNING:                                                            !
!  -------                                                             !
!                                                                      !
!  The background cost function term is computed in v-space. Recall    !
!  that in the inner loop the increment vector deltaX (x-space) and    !
!  deltaV (v-space) are written into ITL(ng)%name NetCDF file at       !
!  records 1 and 2, respectively.                                      !
!                                                                      !
!=======================================================================
!
      implicit none

      CONTAINS
!
!***********************************************************************
      SUBROUTINE back_cost (ng, tile, Lsum)
!***********************************************************************
!
      USE mod_param
# ifdef ADJUST_BOUNDARY
      USE mod_boundary
# endif
# if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
      USE mod_forces
# endif
# ifdef MASKING
      USE mod_grid
# endif
      USE mod_ocean
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, Lsum
!
!  Local variable declarations.
!
# include "tile.h"
!
      CALL back_cost_tile (ng, tile,                                    &
     &                     LBi, UBi, LBj, UBj, LBij, UBij,              &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
     &                     Lsum,                                        &
# ifdef MASKING
     &                     GRID(ng) % rmask,                            &
     &                     GRID(ng) % umask,                            &
     &                     GRID(ng) % vmask,                            &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                     BOUNDARY(ng) % tl_t_obc,                     &
     &                     BOUNDARY(ng) % tl_u_obc,                     &
     &                     BOUNDARY(ng) % tl_v_obc,                     &
#  endif
     &                     BOUNDARY(ng) % tl_ubar_obc,                  &
     &                     BOUNDARY(ng) % tl_vbar_obc,                  &
     &                     BOUNDARY(ng) % tl_zeta_obc,                  &
# endif
# ifdef ADJUST_WSTRESS
     &                     FORCES(ng) % tl_ustr,                        &
     &                     FORCES(ng) % tl_vstr,                        &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                     FORCES(ng) % tl_tflux,                       &
#  endif
     &                     OCEAN(ng) % tl_t,                            &
     &                     OCEAN(ng) % tl_u,                            &
     &                     OCEAN(ng) % tl_v,                            &
# else
     &                     OCEAN(ng) % tl_ubar,                         &
     &                     OCEAN(ng) % tl_vbar,                         &
# endif
     &                     OCEAN(ng) % tl_zeta)

      RETURN
      END SUBROUTINE back_cost
!
!***********************************************************************
      SUBROUTINE back_cost_tile (ng, tile,                              &
     &                           LBi, UBi, LBj, UBj, LBij, UBij,        &
     &                           IminS, ImaxS, JminS, JmaxS,            &
     &                           Lsum,                                  &
# ifdef MASKING
     &                           rmask, umask, vmask,                   &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                           tl_t_obc, tl_u_obc, tl_v_obc,          &
#  endif
     &                           tl_ubar_obc, tl_vbar_obc,              &
     &                           tl_zeta_obc,                           &
# endif
# ifdef ADJUST_WSTRESS
     &                           tl_ustr,                               &
     &                           tl_vstr,                               &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                           tl_tflux,                              &
#  endif
     &                           tl_t,                                  &
     &                           tl_u,                                  &
     &                           tl_v,                                  &
# else
     &                           tl_ubar,                               &
     &                           tl_vbar,                               &
# endif
     &                           tl_zeta)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_fourdvar
      USE mod_ncparam
      USE mod_scalars

# ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_reduce
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBij, UBij
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: Lsum
!
# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#  endif
#  ifdef ADJUST_BOUNDARY
#   ifdef SOLVE3D
      real(r8), intent(in) :: tl_t_obc(LBij:,:,:,:,:,:)
      real(r8), intent(in) :: tl_u_obc(LBij:,:,:,:,:)
      real(r8), intent(in) :: tl_v_obc(LBij:,:,:,:,:)
#   endif
      real(r8), intent(in) :: tl_ubar_obc(LBij:,:,:,:)
      real(r8), intent(in) :: tl_vbar_obc(LBij:,:,:,:)
      real(r8), intent(in) :: tl_zeta_obc(LBij:,:,:,:)
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(in) :: tl_ustr(LBi:,LBj:,:,:)
      real(r8), intent(in) :: tl_vstr(LBi:,LBj:,:,:)
#  endif
#  ifdef SOLVE3D
#   ifdef ADJUST_STFLUX
      real(r8), intent(in) :: tl_tflux(LBi:,LBj:,:,:,:)
#   endif
      real(r8), intent(in) :: tl_t(LBi:,LBj:,:,:,:)
      real(r8), intent(in) :: tl_u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: tl_v(LBi:,LBj:,:,:)
#  else
      real(r8), intent(in) :: tl_ubar(LBi:,LBj:,:)
      real(r8), intent(in) :: tl_vbar(LBi:,LBj:,:)
#  endif
      real(r8), intent(in) :: tl_zeta(LBi:,LBj:,:)

# else

#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#  endif
#  ifdef ADJUST_BOUNDARY
#   ifdef SOLVE3D
      real(r8), intent(in) :: tl_t_obc(LBij:UBij,N(ng),4,               &
     &                                 Nbrec(ng),2,NT(ng))
      real(r8), intent(in) :: tl_u_obc(LBij:UBij,N(ng),4,Nbrec(ng),2)
      real(r8), intent(in) :: tl_v_obc(LBij:UBij,N(ng),4,Nbrec(ng),2)
#   endif
      real(r8), intent(in) :: tl_ubar_obc(LBij:UBij,4,Nbrec(ng),2)
      real(r8), intent(in) :: tl_vbar_obc(LBij:UBij,4,Nbrec(ng),2)
      real(r8), intent(in) :: tl_zeta_obc(LBij:UBij,4,Nbrec(ng),2)
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(in) :: tl_ustr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
      real(r8), intent(in) :: tl_vstr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
#  endif
#  ifdef SOLVE3D
#   ifdef ADJUST_STFLUX
      real(r8), intent(in) :: tl_tflux(LBi:UBi,LBj:UBj,                 &
     &                                 Nfrec(ng),2,NT(ng))
#   endif
      real(r8), intent(in) :: tl_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(in) :: tl_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: tl_v(LBi:UBi,LBj:UBj,N(ng),2)
#  else
      real(r8), intent(in) :: tl_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: tl_vbar(LBi:UBi,LBj:UBj,3)
#  endif
      real(r8), intent(in) :: tl_zeta(LBi:UBi,LBj:UBj,3)
# endif
!
!  Local variable declarations.
!
      integer :: NSUB, i, ib, ir, j, k
# ifdef SOLVE3D
      integer :: itrc
# endif

      real(r8), dimension(0:NstateVar(ng)) :: my_BackCost

      real(r8) :: cff1, cff2

# ifdef DISTRIBUTE
      character (len=3), dimension(0:NstateVar(ng)) :: op_handle
# endif

# include "set_bounds.h"
!
!----------------------------------------------------------------------
!  Compute v-space background cost function (Jb) as misfit between
!  model and background states at initial time of the assimilation
!  window using the sum of all the previous outer-loop increments
!  (Lsum).
!
!  Initially, the misfit innovation matrix (X-Xb) is zero. As the
!  assimilation algorithm iterates, Jb becomes greater than zero.
!----------------------------------------------------------------------
!
      DO i=0,NstateVar(ng)
        my_BackCost(i)=0.0_r8
      END DO
!
!  Free-surface contribution.
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          cff1=tl_zeta(i,j,Lsum)
# ifdef MASKING
          cff1=cff1*rmask(i,j)
# endif
          cff2=0.5_r8*cff1*cff1
          my_BackCost(0)=my_BackCost(0)+cff2
          my_BackCost(isFsur)=my_BackCost(isFsur)+cff2
        END DO
      END DO

# ifdef ADJUST_BOUNDARY
!
!  Free-surface open boundaries.
!
      IF (ANY(Lobc(:,isFsur,ng))) THEN
        DO ir=1,Nbrec(ng)
          IF ((Lobc(iwest,isFsur,ng)).and.                              &
     &        DOMAIN(ng)%Western_Edge(tile)) THEN
            ib=iwest
            DO j=Jstr,Jend
              cff1=tl_zeta_obc(j,ib,ir,Lsum)
#  ifdef MASKING
              cff1=cff1*rmask(Istr-1,j)
#  endif
              cff2=0.5_r8*cff1*cff1
              my_BackCost(0)=my_BackCost(0)+cff2
              my_BackCost(isFsur)=my_BackCost(isFsur)+cff2
            END DO
          END IF
          IF ((Lobc(ieast,isFsur,ng)).and.                              &
     &        DOMAIN(ng)%Eastern_Edge(tile)) THEN
            ib=ieast
            DO j=Jstr,Jend
              cff1=tl_zeta_obc(j,ib,ir,Lsum)
#  ifdef MASKING
              cff1=cff1*rmask(Iend+1,j)
#  endif
              cff2=0.5_r8*cff1*cff1
              my_BackCost(0)=my_BackCost(0)+cff2
              my_BackCost(isFsur)=my_BackCost(isFsur)+cff2
            END DO
          END IF
          IF ((Lobc(isouth,isFsur,ng)).and.                             &
     &        DOMAIN(ng)%Southern_Edge(tile)) THEN
            ib=isouth
            DO i=Istr,Iend
              cff1=tl_zeta_obc(i,ib,ir,Lsum)
#  ifdef MASKING
              cff1=cff1*rmask(i,Jstr-1)
#  endif
              cff2=0.5_r8*cff1*cff1
              my_BackCost(0)=my_BackCost(0)+cff2
              my_BackCost(isFsur)=my_BackCost(isFsur)+cff2
            END DO
          END IF
          IF ((Lobc(inorth,isFsur,ng)).and.                             &
     &        DOMAIN(ng)%Northern_Edge(tile)) THEN
            ib=inorth
            DO i=Istr,Iend
              cff1=tl_zeta_obc(i,ib,ir,Lsum)
#  ifdef MASKING
              cff1=cff1*rmask(i,Jend+1)
#  endif
              cff2=0.5_r8*cff1*cff1
              my_BackCost(0)=my_BackCost(0)+cff2
              my_BackCost(isFsur)=my_BackCost(isFsur)+cff2
            END DO
          END IF
        END DO
      END IF
# endif

# ifndef SOLVE3D
!
!  2D U-momentum contribution.
!
      DO j=Jstr,Jend
        DO i=IstrU,Iend
          cff1=tl_ubar(i,j,Lsum)
#  ifdef MASKING
          cff1=cff1*umask(i,j)
#  endif
          cff2=0.5_r8*cff1*cff1
          my_BackCost(0)=my_BackCost(0)+cff2
          my_BackCost(isUbar)=my_BackCost(isUbar)+cff2
        END DO
      END DO
# endif

# ifdef ADJUST_BOUNDARY
!
!  2D U-momentum open boundaries.
!
      IF (ANY(Lobc(:,isUbar,ng))) THEN
        DO ir=1,Nbrec(ng)
          IF ((Lobc(iwest,isUbar,ng)).and.                              &
     &        DOMAIN(ng)%Western_Edge(tile)) THEN
            ib=iwest
            DO j=Jstr,Jend
              cff1=tl_ubar_obc(j,ib,ir,Lsum)
#  ifdef MASKING
              cff1=cff1*umask(Istr,j)
#  endif
              cff2=0.5_r8*cff1*cff1
              my_BackCost(0)=my_BackCost(0)+cff2
              my_BackCost(isUbar)=my_BackCost(isUbar)+cff2
            END DO
          END IF
          IF ((Lobc(ieast,isUbar,ng)).and.                              &
     &        DOMAIN(ng)%Eastern_Edge(tile)) THEN
            ib=ieast
            DO j=Jstr,Jend
              cff1=tl_ubar_obc(j,ib,ir,Lsum)
#  ifdef MASKING
              cff1=cff1*umask(Iend+1,j)
#  endif
              cff2=0.5_r8*cff1*cff1
              my_BackCost(0)=my_BackCost(0)+cff2
              my_BackCost(isUbar)=my_BackCost(isUbar)+cff2
            END DO
          END IF
          IF ((Lobc(isouth,isUbar,ng)).and.                             &
     &        DOMAIN(ng)%Southern_Edge(tile)) THEN
            ib=isouth
            DO i=IstrU,Iend
              cff1=tl_ubar_obc(i,ib,ir,Lsum)
#  ifdef MASKING
              cff1=cff1*umask(i,Jstr-1)
#  endif
              cff2=0.5_r8*cff1*cff1
              my_BackCost(0)=my_BackCost(0)+cff2
              my_BackCost(isUbar)=my_BackCost(isUbar)+cff2
            END DO
          END IF
          IF ((Lobc(inorth,isUbar,ng)).and.                             &
     &        DOMAIN(ng)%Northern_Edge(tile)) THEN
            ib=inorth
            DO i=IstrU,Iend
              cff1=tl_ubar_obc(i,ib,ir,Lsum)
#  ifdef MASKING
              cff1=cff1*umask(i,Jend+1)
#  endif
              cff2=0.5_r8*cff1*cff1
              my_BackCost(0)=my_BackCost(0)+cff2
              my_BackCost(isUbar)=my_BackCost(isUbar)+cff2
            END DO
          END IF
        END DO
      END IF
# endif

# ifndef SOLVE3D
!
!  2D V-momentum contribution.
!
      DO j=JstrV,Jend
        DO i=Istr,Iend
          cff1=tl_vbar(i,j,Lsum)
#  ifdef MASKING
          cff1=cff1*vmask(i,j)
#  endif
          cff2=0.5_r8*cff1*cff1
          my_BackCost(0)=my_BackCost(0)+cff2
          my_BackCost(isVbar)=my_BackCost(isVbar)+cff2
        END DO
      END DO
# endif

# ifdef ADJUST_BOUNDARY
!
!  2D V-momentum open boundaries.
!
      IF (ANY(Lobc(:,isVbar,ng))) THEN
        DO ir=1,Nbrec(ng)
          IF ((Lobc(iwest,isVbar,ng)).and.                              &
     &        DOMAIN(ng)%Western_Edge(tile)) THEN
            ib=iwest
            DO j=JstrV,Jend
              cff1=tl_vbar_obc(j,ib,ir,Lsum)
#  ifdef MASKING
              cff1=cff1*vmask(Istr-1,j)
#  endif
              cff2=0.5_r8*cff1*cff1
              my_BackCost(0)=my_BackCost(0)+cff2
              my_BackCost(isVbar)=my_BackCost(isVbar)+cff2
            END DO
          END IF
          IF ((Lobc(ieast,isVbar,ng)).and.                              &
     &        DOMAIN(ng)%Eastern_Edge(tile)) THEN
            ib=ieast
            DO j=JstrV,Jend
              cff1=tl_vbar_obc(j,ib,ir,Lsum)
#  ifdef MASKING
              cff1=cff1*vmask(Iend+1,j)
#  endif
              cff2=0.5_r8*cff1*cff1
              my_BackCost(0)=my_BackCost(0)+cff2
              my_BackCost(isVbar)=my_BackCost(isVbar)+cff2
            END DO
          END IF
          IF ((Lobc(isouth,isVbar,ng)).and.                             &
     &        DOMAIN(ng)%Southern_Edge(tile)) THEN
            ib=isouth
            DO i=Istr,Iend
              cff1=tl_vbar_obc(i,ib,ir,Lsum)
#  ifdef MASKING
              cff1=cff1*vmask(i,Jstr)
#  endif
              cff2=0.5_r8*cff1*cff1
              my_BackCost(0)=my_BackCost(0)+cff2
              my_BackCost(isVbar)=my_BackCost(isVbar)+cff2
            END DO
          END IF
          IF ((Lobc(inorth,isVbar,ng)).and.                             &
     &        DOMAIN(ng)%Northern_Edge(tile)) THEN
            ib=inorth
            DO i=Istr,Iend
              cff1=tl_vbar_obc(i,ib,ir,Lsum)
#  ifdef MASKING
              cff1=cff1*vmask(i,Jend+1)
#  endif
              cff2=0.5_r8*cff1*cff1
              my_BackCost(0)=my_BackCost(0)+cff2
              my_BackCost(isVbar)=my_BackCost(isVbar)+cff2
            END DO
          END IF
        END DO
      END IF
# endif

# ifdef ADJUST_WSTRESS
!
!  Surface momentum stress contribution.
!
      DO ir=1,Nfrec(ng)
        DO j=JstrR,JendR
          DO i=Istr,IendR
            cff1=tl_ustr(i,j,ir,Lsum)
#  ifdef MASKING
            cff1=cff1*umask(i,j)
#  endif
            cff2=0.5_r8*cff1*cff1
            my_BackCost(0)=my_BackCost(0)+cff2
            my_BackCost(isUstr)=my_BackCost(isUstr)+cff2
          END DO
        END DO
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            cff1=tl_vstr(i,j,ir,Lsum)
#  ifdef MASKING
            cff1=cff1*vmask(i,j)
#  endif
            cff2=0.5_r8*cff1*cff1
            my_BackCost(0)=my_BackCost(0)+cff2
            my_BackCost(isVstr)=my_BackCost(isVstr)+cff2
          END DO
        END DO
      END DO
# endif

# ifdef SOLVE3D
!
!  3D U-momentum contribution.
!
      DO k=1,N(ng)
        DO j=Jstr,Jend
          DO i=IstrU,Iend
            cff1=tl_u(i,j,k,Lsum)
#  ifdef MASKING
            cff1=cff1*umask(i,j)
#  endif
            cff2=0.5_r8*cff1*cff1
            my_BackCost(0)=my_BackCost(0)+cff2
            my_BackCost(isUvel)=my_BackCost(isUvel)+cff2
          END DO
        END DO
      END DO

#  ifdef ADJUST_BOUNDARY
!
!  3D U-momentum open boundaries.
!
      IF (ANY(Lobc(:,isUvel,ng))) THEN
        DO ir=1,Nbrec(ng)
          IF ((Lobc(iwest,isUvel,ng)).and.                              &
     &        DOMAIN(ng)%Western_Edge(tile)) THEN
            ib=iwest
            DO k=1,N(ng)
              DO j=Jstr,Jend
                cff1=tl_u_obc(j,k,ib,ir,Lsum)
#   ifdef MASKING
                cff1=cff1*umask(Istr,j)
#   endif
                cff2=0.5_r8*cff1*cff1
                my_BackCost(0)=my_BackCost(0)+cff2
                my_BackCost(isUvel)=my_BackCost(isUvel)+cff2
              END DO
            END DO
          END IF
          IF ((Lobc(ieast,isUvel,ng)).and.                              &
     &        DOMAIN(ng)%Eastern_Edge(tile)) THEN
            ib=ieast
            DO k=1,N(ng)
              DO j=Jstr,Jend
                cff1=tl_u_obc(j,k,ib,ir,Lsum)
#   ifdef MASKING
                cff1=cff1*umask(Iend+1,j)
#   endif
                cff2=0.5_r8*cff1*cff1
                my_BackCost(0)=my_BackCost(0)+cff2
                my_BackCost(isUvel)=my_BackCost(isUvel)+cff2
              END DO
            END DO
          END IF
          IF ((Lobc(isouth,isUvel,ng)).and.                             &
     &        DOMAIN(ng)%Southern_Edge(tile)) THEN
            ib=isouth
            DO k=1,N(ng)
              DO i=IstrU,Iend
                cff1=tl_u_obc(i,k,ib,ir,Lsum)
#   ifdef MASKING
                cff1=cff1*umask(i,Jstr-1)
#   endif
                cff2=0.5_r8*cff1*cff1
                my_BackCost(0)=my_BackCost(0)+cff2
                my_BackCost(isUvel)=my_BackCost(isUvel)+cff2
              END DO
            END DO
          END IF
          IF ((Lobc(inorth,isUvel,ng)).and.                             &
     &        DOMAIN(ng)%Northern_Edge(tile)) THEN
            ib=inorth
            DO k=1,N(ng)
              DO i=IstrU,Iend
                cff1=tl_u_obc(i,k,ib,ir,Lsum)
#   ifdef MASKING
                cff1=cff1*umask(i,Jend+1)
#   endif
                cff2=0.5_r8*cff1*cff1
                my_BackCost(0)=my_BackCost(0)+cff2
                my_BackCost(isUvel)=my_BackCost(isUvel)+cff2
              END DO
            END DO
          END IF
        END DO
      END IF
#  endif
!
!  3D V-momentum contribution.
!
      DO k=1,N(ng)
        DO j=JstrV,Jend
          DO i=Istr,Iend
            cff1=tl_v(i,j,k,Lsum)
#  ifdef MASKING
            cff1=cff1*vmask(i,j)
#  endif
            cff2=0.5_r8*cff1*cff1
            my_BackCost(0)=my_BackCost(0)+cff2
            my_BackCost(isVvel)=my_BackCost(isVvel)+cff2
          END DO
        END DO
      END DO

#  ifdef ADJUST_BOUNDARY
!
!  3D V-momentum open boundaries.
!
      IF (ANY(Lobc(:,isVvel,ng))) THEN
        DO ir=1,Nbrec(ng)
          IF ((Lobc(iwest,isVvel,ng)).and.                              &
     &        DOMAIN(ng)%Western_Edge(tile)) THEN
            ib=iwest
            DO k=1,N(ng)
              DO j=JstrV,Jend
                cff1=tl_v_obc(j,k,ib,ir,Lsum)
#   ifdef MASKING
                cff1=cff1*vmask(Istr-1,j)
#   endif
                cff2=0.5_r8*cff1*cff1
                my_BackCost(0)=my_BackCost(0)+cff2
                my_BackCost(isVvel)=my_BackCost(isVvel)+cff2
              END DO
            END DO
          END IF
          IF ((Lobc(ieast,isVvel,ng)).and.                              &
     &        DOMAIN(ng)%Eastern_Edge(tile)) THEN
            ib=ieast
            DO k=1,N(ng)
              DO j=JstrV,Jend
                cff1=tl_v_obc(j,k,ib,ir,Lsum)
#   ifdef MASKING
                cff1=cff1*vmask(Iend+1,j)
#   endif
                cff2=0.5_r8*cff1*cff1
                my_BackCost(0)=my_BackCost(0)+cff2
                my_BackCost(isVvel)=my_BackCost(isVvel)+cff2
              END DO
            END DO
          END IF
          IF ((Lobc(isouth,isVvel,ng)).and.                             &
     &        DOMAIN(ng)%Southern_Edge(tile)) THEN
            ib=isouth
            DO k=1,N(ng)
              DO i=Istr,Iend
                cff1=tl_v_obc(i,k,ib,ir,Lsum)
#   ifdef MASKING
                cff1=cff1*vmask(i,Jstr)
#   endif
                cff2=0.5_r8*cff1*cff1
                my_BackCost(0)=my_BackCost(0)+cff2
                my_BackCost(isVvel)=my_BackCost(isVvel)+cff2
              END DO
            END DO
          END IF
          IF ((Lobc(inorth,isVvel,ng)).and.                             &
     &        DOMAIN(ng)%Northern_Edge(tile)) THEN
            ib=inorth
            DO k=1,N(ng)
              DO i=Istr,Iend
                cff1=tl_v_obc(i,k,ib,ir,Lsum)
#   ifdef MASKING
                cff1=cff1*vmask(i,Jend+1)
#   endif
                cff2=0.5_r8*cff1*cff1
                my_BackCost(0)=my_BackCost(0)+cff2
                my_BackCost(isVvel)=my_BackCost(isVvel)+cff2
              END DO
            END DO
          END IF
        END DO
      END IF
#  endif
!
!  Tracers contribution.
!
      DO itrc=1,NT(ng)
        DO k=1,N(ng)
          DO j=Jstr,Jend
            DO i=Istr,Iend
              cff1=tl_t(i,j,k,Lsum,itrc)
#  ifdef MASKING
              cff1=cff1*rmask(i,j)
#  endif
              cff2=0.5_r8*cff1*cff1
              my_BackCost(0)=my_BackCost(0)+cff2
              my_BackCost(isTvar(itrc))=my_BackCost(isTvar(itrc))+cff2
            END DO
          END DO
        END DO
      END DO

#  ifdef ADJUST_BOUNDARY
!
!  Tracers open boundaries.
!
      DO itrc=1,NT(ng)
        IF (ANY(Lobc(:,isTvar(itrc),ng))) THEN
          DO ir=1,Nbrec(ng)
            IF ((Lobc(iwest,isTvar(itrc),ng)).and.                      &
     &          DOMAIN(ng)%Western_Edge(tile)) THEN
              ib=iwest
              DO k=1,N(ng)
                DO j=Jstr,Jend
                  cff1=tl_t_obc(j,k,ib,ir,Lsum,itrc)
#   ifdef MASKING
                  cff1=cff1*rmask(Istr-1,j)
#   endif
                  cff2=0.5_r8*cff1*cff1
                  my_BackCost(0)=my_BackCost(0)+cff2
                  my_BackCost(isTvar(itrc))=my_BackCost(isTvar(itrc))+  &
     &                                      cff2
                END DO
              END DO
            END IF
            IF ((Lobc(ieast,isTvar(itrc),ng)).and.                      &
     &          DOMAIN(ng)%Eastern_Edge(tile)) THEN
              ib=ieast
              DO k=1,N(ng)
                DO j=Jstr,Jend
                  cff1=tl_t_obc(j,k,ib,ir,Lsum,itrc)
#   ifdef MASKING
                  cff1=cff1*rmask(Iend+1,j)
#   endif
                  cff2=0.5_r8*cff1*cff1
                  my_BackCost(0)=my_BackCost(0)+cff2
                  my_BackCost(isTvar(itrc))=my_BackCost(isTvar(itrc))+  &
     &                                      cff2
                END DO
              END DO
            END IF
            IF ((Lobc(isouth,isTvar(itrc),ng)).and.                     &
     &          DOMAIN(ng)%Southern_Edge(tile)) THEN
              ib=isouth
              DO k=1,N(ng)
                DO i=Istr,Iend
                  cff1=tl_t_obc(i,k,ib,ir,Lsum,itrc)
#   ifdef MASKING
                  cff1=cff1*rmask(i,Jstr-1)
#   endif
                  cff2=0.5_r8*cff1*cff1
                  my_BackCost(0)=my_BackCost(0)+cff2
                  my_BackCost(isTvar(itrc))=my_BackCost(isTvar(itrc))+  &
     &                                      cff2
                END DO
              END DO
            END IF
            IF ((Lobc(inorth,isTvar(itrc),ng)).and.                     &
     &          DOMAIN(ng)%Northern_Edge(tile)) THEN
              ib=inorth
              DO k=1,N(ng)
                DO i=Istr,Iend
                  cff1=tl_t_obc(i,k,ib,ir,Lsum,itrc)
#   ifdef MASKING
                  cff1=cff1*rmask(i,Jend+1)
#   endif
                  cff2=0.5_r8*cff1*cff1
                  my_BackCost(0)=my_BackCost(0)+cff2
                  my_BackCost(isTvar(itrc))=my_BackCost(isTvar(itrc))+  &
     &                                      cff2
                END DO
              END DO
            END IF
          END DO
        END IF
      END DO
#  endif

#  ifdef ADJUST_STFLUX
!
!  Surface tracers flux contribution.
!
      DO itrc=1,NT(ng)
        IF (Lstflux(itrc,ng)) THEN
          DO ir=1,Nfrec(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                cff1=tl_tflux(i,j,ir,Lsum,itrc)
#   ifdef MASKING
                cff1=cff1*rmask(i,j)
#   endif
                cff2=0.5_r8*cff1*cff1
                my_BackCost(0)=my_BackCost(0)+cff2
                my_BackCost(isTsur(itrc))=my_BackCost(isTsur(itrc))+cff2
              END DO
            END DO
          END DO
        END IF
      END DO
#  endif
# endif
!
!-----------------------------------------------------------------------
! Compute global background cost function.
!-----------------------------------------------------------------------
!
# ifdef DISTRIBUTE
      NSUB=1                             ! distributed-memory
# else
      IF (DOMAIN(ng)%SouthWest_Corner(tile).and.                        &
     &    DOMAIN(ng)%NorthEast_Corner(tile)) THEN
        NSUB=1                           ! non-tiled application
      ELSE
        NSUB=NtileX(ng)*NtileE(ng)       ! tiled application
      END IF
# endif
!$OMP CRITICAL (BACKCOST)
      IF (tile_count.eq.0) THEN
        DO i=0,NstateVar(ng)
          FOURDVAR(ng)%BackCost(i)=my_BackCost(i)
        END DO
      ELSE
        DO i=0,NstateVar(ng)
          FOURDVAR(ng)%BackCost(i)=FOURDVAR(ng)%BackCost(i)+            &
     &                             my_BackCost(i)
        END DO
      END IF
      tile_count=tile_count+1
      IF (tile_count.eq.NSUB) THEN
        tile_count=0
# ifdef DISTRIBUTE
        DO i=0,NstateVar(ng)
          op_handle(i)='SUM'
        END DO
        CALL mp_reduce (ng, iNLM, NstateVar(ng)+1,                      &
     &                  FOURDVAR(ng)%BackCost(0:),  op_handle(0:))
# endif
      END IF
!$OMP END CRITICAL (BACKCOST)

      RETURN
      END SUBROUTINE back_cost_tile
#endif
      END MODULE back_cost_mod
