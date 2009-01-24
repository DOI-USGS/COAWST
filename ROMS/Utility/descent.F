#include "cppdefs.h"
      MODULE descent_mod

#ifdef IS4DVAR_OLD
!
!svn $Id: descent.F 733 2008-09-07 01:56:45Z jcwarner $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine estimates the  "best" state initial conditions using   !
!  a descent algorithm. The scheme used to find the minimum function   !
!  is selected according to the values of parameter ICG:               !
!                                                                      !
!     ICG=0     Fletcher-Reeves scheme                                 !
!     ICG=1     Polak-Ribiere scheme                                   !
!                                                                      !
!  In 2D applications, the control vector at initialization contains   !
!  ZETA, UBAR, and VBAR.  In 3D application,  the control vector  at   !
!  initialization contains ZETA, U, V, and TRACERS.                    !
!                                                                      !
!  This routine assumes that the gradient solution, contained in the   !
!  adjoint state arrays, has been preconditioned with the background   !
!  error covariance matrix  via a space convolution of the diffusion   !
!  equation for each state variable.  That is, the gradient solution   !
!  is in V-space such that:                                            !
!                                                                      !
!        V = B^(1/2) X                                                 !
!                                                                      !
!  where B is the background error covariance matrix.  The  estimate   !
!  for the new state initial conditions are also in V-space.           !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     tile      Sub-domain partition.                                  !
!     lock      Lock key for parallel  reduction operations.           !
!     Iter      Current iteration.                                     !
!     step      Conjugate direction step size (nondimensional).        !
!                                                                      !
!  References:                                                         !
!                                                                      !
!     Fletcher, R. and. C.M. Reeves, 1964: Function minimization       !
!       by conjugate gradients, Comput. J., 7, 149-154.                !
!                                                                      !
!     Polak, E and G. Ribiere, 1969:  Note sur la convergence de       !
!       methodes de directions conjugees, Rev. Fr. Inform. Rech.       !
!       Oper., 16-R1, 35-43.                                           !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: descent

      CONTAINS
!
!***********************************************************************
      SUBROUTINE descent (ng, tile, model, Iter, step)
!***********************************************************************
!
      USE mod_param
# ifdef SOLVE3D
      USE mod_coupling
# endif
      USE mod_grid
      USE mod_ocean
      USE mod_stepping
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model, Iter

      real(r8), intent(in) :: step
!
!  Local variable declarations.
!
# include "tile.h"
!
# ifdef PROFILE
      CALL wclock_on (ng, model, 36)
# endif
      CALL descent_tile (ng, tile, model,                               &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   Lold(ng), Lnew(ng), Iter, step,                &
# ifdef MASKING
     &                   GRID(ng) % rmask,                              &
     &                   GRID(ng) % umask,                              &
     &                   GRID(ng) % vmask,                              &
# endif
# ifdef SOLVE3D
     &                   OCEAN(ng) % tl_t,                              &
     &                   OCEAN(ng) % tl_u,                              &
     &                   OCEAN(ng) % tl_v,                              &
# endif
     &                   OCEAN(ng) % tl_ubar,                           &
     &                   OCEAN(ng) % tl_vbar,                           &
     &                   OCEAN(ng) % tl_zeta,                           &
# ifdef SOLVE3D
     &                   OCEAN(ng) % d_t,                               &
     &                   OCEAN(ng) % d_u,                               &
     &                   OCEAN(ng) % d_v,                               &
# endif
     &                   OCEAN(ng) % d_ubar,                            &
     &                   OCEAN(ng) % d_vbar,                            &
     &                   OCEAN(ng) % d_zeta,                            &
# ifdef SOLVE3D
     &                   OCEAN(ng) % ad_t,                              &
     &                   OCEAN(ng) % ad_u,                              &
     &                   OCEAN(ng) % ad_v,                              &
# endif
     &                   OCEAN(ng) % ad_ubar,                           &
     &                   OCEAN(ng) % ad_vbar,                           &
     &                   OCEAN(ng) % ad_zeta)
# ifdef PROFILE
      CALL wclock_on (ng, model, 36)
# endif
      RETURN
      END SUBROUTINE descent
!
!***********************************************************************
      SUBROUTINE descent_tile (ng, tile, model,                         &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         Lold, Lnew, Iter, step,                  &
# ifdef MASKING
     &                         rmask, umask, vmask,                     &
# endif
# ifdef SOLVE3D
     &                         tl_t, tl_u, tl_v,                        &
# endif
     &                         tl_ubar, tl_vbar, tl_zeta,               &
# ifdef SOLVE3D
     &                         d_t, d_u, d_v,                           &
# endif
     &                         d_ubar, d_vbar, d_zeta,                  &
# ifdef SOLVE3D
     &                         ad_t, ad_u, ad_v,                        &
# endif
     &                         ad_ubar, ad_vbar, ad_zeta)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_fourdvar
      USE mod_scalars

# ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_reduce
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: Lold, Lnew, Iter

      real(r8), intent(in) :: step
!
# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#  endif
#  ifdef SOLVE3D
      real(r8), intent(in) :: ad_t(LBi:,LBj:,:,:,:)
      real(r8), intent(in) :: ad_u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: ad_v(LBi:,LBj:,:,:)
#  endif
      real(r8), intent(in) :: ad_ubar(LBi:,LBj:,:)
      real(r8), intent(in) :: ad_vbar(LBi:,LBj:,:)
      real(r8), intent(in) :: ad_zeta(LBi:,LBj:,:)
#  ifdef SOLVE3D
      real(r8), intent(inout) :: d_t(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: d_u(LBi:,LBj:,:)
      real(r8), intent(inout) :: d_v(LBi:,LBj:,:)
#  endif
      real(r8), intent(inout) :: d_ubar(LBi:,LBj:)
      real(r8), intent(inout) :: d_vbar(LBi:,LBj:)
      real(r8), intent(inout) :: d_zeta(LBi:,LBj:)
#  ifdef SOLVE3D
      real(r8), intent(inout) :: tl_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: tl_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: tl_v(LBi:,LBj:,:,:)
#  endif
      real(r8), intent(inout) :: tl_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: tl_vbar(LBi:,LBj:,:)
      real(r8), intent(inout) :: tl_zeta(LBi:,LBj:,:)

# else

#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#  endif
#  ifdef SOLVE3D
      real(r8), intent(in) :: ad_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(in) :: ad_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: ad_v(LBi:UBi,LBj:UBj,N(ng),2)
#  endif
      real(r8), intent(in) :: ad_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: ad_vbar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: ad_zeta(LBi:UBi,LBj:UBj,3)
#  ifdef SOLVE3D
      real(r8), intent(inout) :: d_t(LBi:UBi,LBj:UBj,N(ng),NT(ng))
      real(r8), intent(inout) :: d_u(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: d_v(LBi:UBi,LBj:UBj,N(ng))
#  endif
      real(r8), intent(inout) :: d_ubar(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: d_vbar(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: d_zeta(LBi:UBi,LBj:UBj)
#  ifdef SOLVE3D
      real(r8), intent(inout) :: tl_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: tl_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: tl_v(LBi:UBi,LBj:UBj,N(ng),2)
#  endif
      real(r8), intent(inout) :: tl_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: tl_vbar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: tl_zeta(LBi:UBi,LBj:UBj,3)
# endif
!
!  Local variable declarations.
!
      integer :: NSUB, i, j
# ifdef SOLVE3D
      integer :: itrc, k
# endif
      real(r8), save :: BetaK

      real(r8) :: CGscheme
      real(r8) :: beta1, beta2, my_beta1, my_beta2
      real(r8) :: dot1, dot2,  my_dot1, my_dot2
      real(r8) :: cff, cff1, cff2, val
# ifdef DISTRIBUTE
      real(r8), dimension(4) :: buffer

      character (len=3), dimension(4) :: op_handle
# endif
!
# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  On first pass, compute the dot product between previous and new
!  adjoint solutions. Compute the Polak-Ribiere scaling factor (BetaK)
!  for the conjugate vectors.
!-----------------------------------------------------------------------
!
      IF ((Ipass.eq.1).and.(Iter.gt.1)) THEN
        my_dot1=0.0_r8
        my_dot2=0.0_r8
        my_beta1=0.0_r8
        my_beta2=0.0_r8
        CGscheme=REAL(ICG,r8)
!
!  2D state variables.
!
# ifndef SOLVE3D
        DO j=JstrR,JendR
          DO i=Istr,IendR
            cff1=ad_ubar(i,j,Lnew)*ad_ubar(i,j,Lnew)
            cff2=ad_ubar(i,j,Lold)*ad_ubar(i,j,Lold)
            my_dot1=my_dot1+cff1
            my_dot2=my_dot2+                                            &
     &              ad_ubar(i,j,Lold)*ad_ubar(i,j,Lnew)
            my_beta1=my_beta1+cff2
            my_beta2=my_beta2+                                          &
     &               ad_ubar(i,j,Lnew)*                                 &
     &               (-CGscheme*ad_ubar(i,j,Lold)+                      &
     &                ad_ubar(i,j,Lnew))
          END DO
        END DO
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            cff1=ad_vbar(i,j,Lnew)*ad_vbar(i,j,Lnew)
            cff2=ad_vbar(i,j,Lold)*ad_vbar(i,j,Lold)
            my_dot1=my_dot1+cff1
            my_dot2=my_dot2+                                            &
     &              ad_vbar(i,j,Lold)*ad_vbar(i,j,Lnew)
            my_beta1=my_beta1+cff2
            my_beta2=my_beta2+                                          &
     &               ad_vbar(i,j,Lnew)*                                 &
     &               (-CGscheme*ad_vbar(i,j,Lold)+                      &
     &                ad_vbar(i,j,Lnew))
          END DO
        END DO
# endif
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            cff1=ad_zeta(i,j,Lnew)*ad_zeta(i,j,Lnew)
            cff2=ad_zeta(i,j,Lold)*ad_zeta(i,j,Lold)
            my_dot1=my_dot1+cff1
            my_dot2=my_dot2+                                            &
     &              ad_zeta(i,j,Lold)*ad_zeta(i,j,Lnew)
            my_beta1=my_beta1+cff2
            my_beta2=my_beta2+                                          &
     &               ad_zeta(i,j,Lnew)*                                 &
     &               (-CGscheme*ad_zeta(i,j,Lold)+                      &
     &                ad_zeta(i,j,Lnew))
          END DO
        END DO
# ifdef SOLVE3D
!
!  3D state variables.
!
        DO k=1,N(ng)
          DO j=JstrR,JendR
            DO i=Istr,IendR
              cff1=ad_u(i,j,k,Lnew)*ad_u(i,j,k,Lnew)
              cff2=ad_u(i,j,k,Lold)*ad_u(i,j,k,Lold)
              my_dot1=my_dot1+cff1
              my_dot2=my_dot2+                                          &
     &                ad_u(i,j,k,Lold)*ad_u(i,j,k,Lnew)
              my_beta1=my_beta1+cff2
              my_beta2=my_beta2+                                        &
     &                 ad_u(i,j,k,Lnew)*                                &
     &                 (-CGscheme*ad_u(i,j,k,Lold)+                     &
     &                  ad_u(i,j,k,Lnew))
            END DO
          END DO
        END DO
        DO k=1,N(ng)
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              cff1=ad_v(i,j,k,Lnew)*ad_v(i,j,k,Lnew)
              cff2=ad_v(i,j,k,Lold)*ad_v(i,j,k,Lold)
              my_dot1=my_dot1+cff1
              my_dot2=my_dot2+                                          &
     &                ad_v(i,j,k,Lold)*ad_v(i,j,k, Lnew)
              my_beta1=my_beta1+cff2
              my_beta2=my_beta2+                                        &
     &                 ad_v(i,j,k,Lnew)*                                &
     &                 (-CGscheme*ad_v(i,j,k,Lold)+                     &
     &                  ad_v(i,j,k,Lnew))
            END DO
          END DO
        END DO
        DO itrc=1,NT(ng)
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                cff1=ad_t(i,j,k,Lnew,itrc)*                             &
     &               ad_t(i,j,k,Lnew,itrc)
                cff2=ad_t(i,j,k,Lold,itrc)*                             &
     &               ad_t(i,j,k,Lold,itrc)
                my_dot1=my_dot1+cff1
                my_dot2=my_dot2+                                        &
     &                  ad_t(i,j,k,Lold,itrc)*                          &
     &                  ad_t(i,j,k,Lnew,itrc)
                my_beta1=my_beta1+cff2
                my_beta2=my_beta2+                                      &
     &                   ad_t(i,j,k,Lnew,itrc)*                         &
     &                   (-CGscheme*                                    &
     &                    ad_t(i,j,k,Lold,itrc)+                        &
     &                    ad_t(i,j,k,Lnew,itrc))
              END DO
            END DO
          END DO
        END DO
# endif
!
!  Perform parallel global reduction operations.
!
        IF (SOUTH_WEST_CORNER.and.                                      &
     &      NORTH_EAST_CORNER) THEN
          NSUB=1                         ! non-tiled application
        ELSE
          NSUB=NtileX(ng)*NtileE(ng)     ! tiled application
        END IF
!$OMP CRITICAL (TL_DOT)
        IF (tile_count.eq.0) THEN
          dot1=0.0_r8
          dot2=0.0_r8
          beta1=0.0_r8
          beta2=0.0_r8
        END IF
        dot1=dot1+my_dot1
        dot2=dot2+my_dot2
        beta1=beta1+my_beta1
        beta2=beta2+my_beta2
        tile_count=tile_count+1
        IF (tile_count.eq.NSUB) THEN
          tile_count=0
#ifdef DISTRIBUTE
          buffer(1)=dot1
          buffer(2)=dot2
          buffer(3)=beta1
          buffer(4)=beta2
          op_handle(1)='SUM'
          op_handle(2)='SUM'
          op_handle(3)='SUM'
          op_handle(4)='SUM'
          CALL mp_reduce (ng, model, 4, buffer, op_handle)
          dot1=buffer(1)
          dot2=buffer(2)
          beta1=buffer(3)
          beta2=buffer(4)
# endif
        END IF
!$OMP END CRITICAL (TL_DOT)
        dot1=CGtol*ABS(dot1)
        dot2=ABS(dot2)
        BetaK=beta2/beta1
!
!  Restart with steepest decent every "NiterSD" iterations.
!
        IF (MOD(Iter-IterSD,NiterSD).eq.0) THEN
          BetaK=0.0
          IterSD=Iter
        END IF
!
!  Perform conjugacy test and perform steepest descent if necessary.
!
        IF (dot2.gt.dot1) THEN
          BetaK=0.0_r8
          IterSD=Iter
        END IF
!
!  If first pass and first iteration, use steepest descent algorithm.
!
      ELSE IF ((Ipass.eq.1).and.(Iter.eq.1)) THEN
        dot1=0.0_r8
        dot2=0.0_r8
        BetaK=0.0_r8
      END IF
!
!-----------------------------------------------------------------------
!  Calculate conjugate vectors (descent directions) and new initial
!  conditions. Notice that the conjugate vectors are only saved on the
!  second pass.
!-----------------------------------------------------------------------
!
      IF (Master) THEN
        IF (Ipass.eq.1) THEN
          PRINT 10, dot1, dot2, IterSD
 10       FORMAT (/,' DESCENT - old state dot product, dot1 = ',        &
     &            1p,e15.8,                                             &
     &            /,11x,'new state dot product, dot2 = ',1p,e15.8,      &
     &            /,11x,'Last steepest descent iteration = ',i5.5,/)
        END IF
        PRINT 20, Nrun, Iter, Ipass, BetaK
 20     FORMAT (/,' DESCENT - conjugate vector scaling factor:',        &
     &          /,11x,'(Iter=',i4.4,', inner=',i3.3,', Ipass=',i1,')',  &
     &          ' BetaK = ',1p,e15.8)
      END IF      
!
!  First pass, 2D state variables.
!
      IF (Ipass.eq.1) THEN
# ifndef SOLVE3D
        DO j=JstrR,JendR
          DO i=Istr,IendR
            cff=-ad_ubar(i,j,Lnew)+BetaK*d_ubar(i,j)
#  ifdef MASKING
            cff=cff*umask(i,j)
#  endif
            tl_ubar(i,j,Lnew)=tl_ubar(i,j,Lold)+step*cff
          END DO
        END DO
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            cff=-ad_vbar(i,j,Lnew)+BetaK*d_vbar(i,j)
#  ifdef MASKING
            cff=cff*vmask(i,j)
#  endif
            tl_vbar(i,j,Lnew)=tl_vbar(i,j,Lold)+step*cff
          END DO
        END DO
# endif
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            cff=-ad_zeta(i,j,Lnew)+BetaK*d_zeta(i,j)
# ifdef MASKING
            cff=cff*rmask(i,j)
# endif          
            tl_zeta(i,j,Lnew)=tl_zeta(i,j,Lold)+step*cff
          END DO
        END DO
# ifdef SOLVE3D
!
!  First pass, 3D state variables.
!
        DO k=1,N(ng)
          DO j=JstrR,JendR
            DO i=Istr,IendR
              cff=-ad_u(i,j,k,Lnew)+BetaK*d_u(i,j,k)
#  ifdef MASKING
              cff=cff*umask(i,j)
#  endif
              tl_u(i,j,k,Lnew)=tl_u(i,j,k,Lold)+step*cff
            END DO
          END DO
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              cff=-ad_v(i,j,k,Lnew)+BetaK*d_v(i,j,k)
#  ifdef MASKING
              cff=cff*vmask(i,j)
#  endif
              tl_v(i,j,k,Lnew)=tl_v(i,j,k,Lold)+step*cff
            END DO
          END DO
        END DO
!
        DO itrc=1,NT(ng)
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                cff=-ad_t(i,j,k,Lnew,itrc)+BetaK*d_t(i,j,k,itrc)
#  ifdef MASKING
                cff=cff*rmask(i,j)
#  endif          
                tl_t(i,j,k,Lnew,itrc)=tl_t(i,j,k,Lold,itrc)+            &
     &                                step*cff
              END DO
            END DO          
          END DO
        END DO
# endif
      ELSE IF (Ipass.eq.2 ) THEN
!
!  Second pass, 2D state variables.
!
# ifndef SOLVE3D
        DO j=JstrR,JendR
          DO i=Istr,IendR
            cff=-ad_ubar(i,j,Lnew)+BetaK*d_ubar(i,j)
#  ifdef MASKING
            cff=cff*umask(i,j)
#  endif
            d_ubar(i,j)=cff
            tl_ubar(i,j,Lnew)=tl_ubar(i,j,Lold)+step*cff
          END DO
        END DO
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            cff=-ad_vbar(i,j,Lnew)+BetaK*d_vbar(i,j)
#  ifdef MASKING
            cff=cff*vmask(i,j)
#  endif
            d_vbar(i,j)=cff
            tl_vbar(i,j,Lnew)=tl_vbar(i,j,Lold)+step*cff
          END DO
        END DO
# endif
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            cff=-ad_zeta(i,j,Lnew)+BetaK*d_zeta(i,j)
# ifdef MASKING
            cff=cff*rmask(i,j)
# endif          
            d_zeta(i,j)=cff
            tl_zeta(i,j,Lnew)=tl_zeta(i,j,Lold)+step*cff
          END DO
        END DO
# ifdef SOLVE3D
!
!  Second pass, 3D state variables.
!
        DO k=1,N(ng)
          DO j=JstrR,JendR
            DO i=Istr,IendR
              cff=-ad_u(i,j,k,Lnew)+BetaK*d_u(i,j,k)
#  ifdef MASKING
              cff=cff*umask(i,j)
#  endif
              d_u(i,j,k)=cff
              tl_u(i,j,k,Lnew)=tl_u(i,j,k,Lold)+step*cff
            END DO
          END DO
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              cff=-ad_v(i,j,k,Lnew)+BetaK*d_v(i,j,k)
#  ifdef MASKING
              cff=cff*vmask(i,j)
#  endif
              d_v(i,j,k)=cff
              tl_v(i,j,k,Lnew)=tl_v(i,j,k,Lold)+step*cff
            END DO
          END DO
        END DO
!
        DO itrc=1,NT(ng)
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                cff=-ad_t(i,j,k,Lnew,itrc)+BetaK*d_t(i,j,k,itrc)
#  ifdef MASKING
                cff=cff*rmask(i,j)
#  endif          
                d_t(i,j,k,itrc)=cff
                tl_t(i,j,k,Lnew,itrc)=tl_t(i,j,k,Lold,itrc)+step*cff
              END DO
            END DO          
          END DO
        END DO
# endif
      END IF     

      RETURN
      END SUBROUTINE descent_tile
#endif
      END MODULE descent_mod
