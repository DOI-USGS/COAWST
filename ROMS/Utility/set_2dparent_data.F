#include "cppdefs.h"
      MODULE set_2dparent_data_mod
#ifdef REFINED_GRID
!
!svn $Id: set_data.F 1405 2008-02-20 01:44:07Z jcwarner $
!================================================== John C. Warner =====
!                                                                      !
!  This routine interpolates the parent data and places it             !
!  into the perimeter of the child grid.                               !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: set_2dparent_data, set_2dparent_data_tile

      CONTAINS
!
!***********************************************************************
      SUBROUTINE set_2dparent_data (ng, tile)
!***********************************************************************
!
      USE mod_param
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
      integer :: ngp, LBip, UBip, LBjp, UBjp
# include "tile.h"

      ngp=myparent(ng)
      LBip=BOUNDS(ngp)%LBi(tile)
      UBip=BOUNDS(ngp)%UBi(tile)
      LBjp=BOUNDS(ngp)%LBj(tile)
      UBjp=BOUNDS(ngp)%UBj(tile)
!
# ifdef PROFILE
      CALL wclock_on (ng, iNLM, 4)
# endif
      CALL set_2dparent_data_tile (ng, tile,                            &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           LBip, UBip, LBjp, UBjp)
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 4)
# endif
      RETURN
      END SUBROUTINE set_2dparent_data
!
!***********************************************************************
      SUBROUTINE set_2dparent_data_tile (ng, tile,                      &
     &                                   LBi, UBi, LBj, UBj,            &
     &                                   LBip, UBip, LBjp, UBjp)
!***********************************************************************
!
      USE mod_param
      USE mod_boundary
      USE mod_grid
      USE mod_iounits
      USE mod_ncparam
      USE mod_scalars
      USE mod_stepping
      USE mod_parallel
      USE mod_ocean
      USE mod_refined
      USE mod_scalars
      USE mod_coupling
      USE interpolate_mod

# if defined EW_PERIODIC || defined NS_PERIODIC
      USE exchange_2d_mod
# endif
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
# endif
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: LBip, UBip, LBjp, UBjp
!
!  Local variable declarations.
!
      logical, dimension(3) :: update =                                 &
     &         (/ .FALSE., .FALSE., .FALSE. /)

# ifdef OBC
      integer :: ILB, IUB, JLB, JUB
# endif
      integer :: i, j, ngp, ic, jc, t1, t2
      real(r8) :: my_val, fac, fac1, fac2, fac3, tscale
      real(r8) :: my_valu, my_valv, my_valh1, my_valh2, cff1
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
!
# include "set_bounds.h"
      ngp=myparent(ng)
!
# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn on input data time wall clock.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, iNLM, 3)
# endif
!
!-----------------------------------------------------------------------
!  Fill local arrays from the BC data.
!-----------------------------------------------------------------------
!
# ifdef SOLVE3D
      fac1=time(ngp)+dt(ngp)-time(ng)-dt(ng)
      fac2=time(ng)-time(ngp)+dt(ng)
      fac=1.0_r8/(fac1+fac2)
      fac1=fac*fac1
      fac2=fac*fac2
      IF (iic(1).eq.ntstart(1)) THEN
        t1=1
        t2=2
      ELSE
        t1=nstp(ngp)
        t2=nnew(ngp)
      END IF
# else
      IF (PREDICTOR_2D_STEP(ng)) THEN
        tscale=2.0_r8
      ELSE
        tscale=1.0_r8
      END IF
      fac1=time(ngp)+dt(ngp)-time(ng)-tscale*dt(ng)
      fac2=time(ng)-time(ngp)+tscale*dt(ng)
      fac=1.0_r8/(fac1+fac2)
      fac1=fac*fac1
      fac2=fac*fac2
      IF (iic(1).eq.ntstart(1)) THEN
        t1=1
        t2=2
      ELSE
        IF (MOD(iic(1)*3,2).eq.0) THEN
          t1=2
          t2=1
        ELSE
          t1=1
          t2=2
        ENDIF
      END IF
# endif

!
! zeta
!
      IF (WESTERN_EDGE) THEN
        DO j=JstrT,JendT
          DO i=IstrT,Istr-1
            ic=i+Nghostpoints+1
            my_val=fac1*REFINED(ng)%zetaR_west(ic,j,t1)+                &
     &             fac2*REFINED(ng)%zetaR_west(ic,j,t2)
# ifdef MASKING
            my_val=my_val*GRID(ng)%rmask(i,j)
# endif
# ifdef WET_DRY
            IF (my_val.le.(Dcrit(ng)-GRID(ng)%h(i,j))) THEN
              my_val=Dcrit(ng)-GRID(ng)%h(i,j)
            END IF
# endif
            OCEAN(ng)%zeta(i,j,knew(ng))=my_val
# ifdef SOLVE3D
            COUPLING(ng)%Zt_avg1(i,j)=my_val
# endif
          END DO
        END DO
      END IF
      IF (EASTERN_EDGE) THEN
        DO j=JstrT,JendT
          DO i=Iend+1,IendT
            ic=i-Iend+1
            my_val=fac1*REFINED(ng)%zetaR_east(ic,j,t1)+                &
     &             fac2*REFINED(ng)%zetaR_east(ic,j,t2)
# ifdef MASKING
            my_val=my_val*GRID(ng)%rmask(i,j)
# endif
# ifdef WET_DRY
            IF (my_val.le.(Dcrit(ng)-GRID(ng)%h(i,j))) THEN
              my_val=Dcrit(ng)-GRID(ng)%h(i,j)
            END IF
# endif
            OCEAN(ng)%zeta(i,j,knew(ng))=my_val
# ifdef SOLVE3D
            COUPLING(ng)%Zt_avg1(i,j)=my_val
# endif
          END DO
        END DO
      END IF
      IF (SOUTHERN_EDGE) THEN
        DO j=JstrT,Jstr-1
          DO i=IstrT,IendT
            jc=j+Nghostpoints+1
            my_val=fac1*REFINED(ng)%zetaR_south(i,jc,t1)+               &
     &             fac2*REFINED(ng)%zetaR_south(i,jc,t2)
# ifdef MASKING
            my_val=my_val*GRID(ng)%rmask(i,j)
# endif
# ifdef WET_DRY
            IF (my_val.le.(Dcrit(ng)-GRID(ng)%h(i,j))) THEN
              my_val=Dcrit(ng)-GRID(ng)%h(i,j)
            END IF
# endif
            OCEAN(ng)%zeta(i,j,knew(ng))=my_val
# ifdef SOLVE3D
            COUPLING(ng)%Zt_avg1(i,j)=my_val
# endif
          END DO
        END DO
      END IF
      IF (NORTHERN_EDGE) THEN
        DO j=Jend+1,JendT
          DO i=IstrT,IendT
            jc=j-Jend+1
            my_val=fac1*REFINED(ng)%zetaR_north(i,jc,t1)+               &
     &             fac2*REFINED(ng)%zetaR_north(i,jc,t2)
# ifdef MASKING
            my_val=my_val*GRID(ng)%rmask(i,j)
# endif
# ifdef WET_DRY
            IF (my_val.le.(Dcrit(ng)-GRID(ng)%h(i,j))) THEN
              my_val=Dcrit(ng)-GRID(ng)%h(i,j)
            END IF
# endif
            OCEAN(ng)%zeta(i,j,knew(ng))=my_val
# ifdef SOLVE3D
            COUPLING(ng)%Zt_avg1(i,j)=my_val
# endif
          END DO
        END DO
      END IF
!
! ubar
!
      IF (WESTERN_EDGE) THEN
        DO j=JstrT,JendT
          DO i=IstrT+1,Istr
            ic=i+Nghostpoints
            my_valu=fac1*REFINED(ng)%ubarR_west(ic,j,t1)+               &
     &              fac2*REFINED(ng)%ubarR_west(ic,j,t2)
# ifdef MASKING
            my_valu=my_valu*GRID(ng)%umask(i,j)
# endif
            OCEAN(ng)%ubar(i,j,knew(ng))=my_valu
          END DO
        END DO
      END IF
      IF (EASTERN_EDGE) THEN
        DO j=JstrT,JendT
          DO i=Iend+1,IendT
            ic=i-Iend
            my_valu=fac1*REFINED(ng)%ubarR_east(ic,j,t1)+               &
     &              fac2*REFINED(ng)%ubarR_east(ic,j,t2)
# ifdef MASKING
            my_valu=my_valu*GRID(ng)%umask(i,j)
# endif
            OCEAN(ng)%ubar(i,j,knew(ng))=my_valu
          END DO
        END DO
      END IF
      IF (SOUTHERN_EDGE) THEN
        DO j=JstrT,Jstr-1
          DO i=IstrTU+1,IendT
            jc=j+Nghostpoints+1
            my_valu=fac1*REFINED(ng)%ubarR_south(i,jc,t1)+              &
     &              fac2*REFINED(ng)%ubarR_south(i,jc,t2)
# ifdef MASKING
            my_valu=my_valu*GRID(ng)%umask(i,j)
# endif
            OCEAN(ng)%ubar(i,j,knew(ng))=my_valu
          END DO
        END DO
      END IF
      IF (NORTHERN_EDGE) THEN
        DO j=Jend+1,JendT
          DO i=IstrTU+1,IendT
            jc=j-Jend
            my_valu=fac1*REFINED(ng)%ubarR_north(i,jc,t1)+              &
     &              fac2*REFINED(ng)%ubarR_north(i,jc,t2)
# ifdef MASKING
            my_valu=my_valu*GRID(ng)%umask(i,j)
# endif
            OCEAN(ng)%ubar(i,j,knew(ng))=my_valu
          END DO
        END DO
      END IF
# ifdef WEC
!
! ubar_stokes
!
      IF (WESTERN_EDGE) THEN
        DO j=JstrT,JendT
          DO i=IstrT+1,Istr
            ic=i+Nghostpoints
            my_valu=fac1*REFINED(ng)%ubarsR_west(ic,j,t1)+              &
     &              fac2*REFINED(ng)%ubarsR_west(ic,j,t2)
#  ifdef MASKING
            my_valu=my_valu*GRID(ng)%umask(i,j)
#  endif
            OCEAN(ng)%ubar_stokes(i,j)=my_valu
          END DO
        END DO
      END IF
      IF (EASTERN_EDGE) THEN
        DO j=JstrT,JendT
          DO i=Iend+1,IendT
            ic=i-Iend
            my_valu=fac1*REFINED(ng)%ubarsR_east(ic,j,t1)+              &
     &              fac2*REFINED(ng)%ubarsR_east(ic,j,t2)
#  ifdef MASKING
            my_valu=my_valu*GRID(ng)%umask(i,j)
#  endif
            OCEAN(ng)%ubar_stokes(i,j)=my_valu
          END DO
        END DO
      END IF
      IF (SOUTHERN_EDGE) THEN
        DO j=JstrT,Jstr-1
          DO i=IstrTU+1,IendT
            jc=j+Nghostpoints+1
            my_valu=fac1*REFINED(ng)%ubarsR_south(i,jc,t1)+             &
     &              fac2*REFINED(ng)%ubarsR_south(i,jc,t2)
#  ifdef MASKING
            my_valu=my_valu*GRID(ng)%umask(i,j)
#  endif
            OCEAN(ng)%ubar_stokes(i,j)=my_valu
          END DO
        END DO
      END IF
      IF (NORTHERN_EDGE) THEN
        DO j=Jend+1,JendT
          DO i=IstrTU+1,IendT
            jc=j-Jend
            my_valu=fac1*REFINED(ng)%ubarsR_north(i,jc,t1)+             &
     &              fac2*REFINED(ng)%ubarsR_north(i,jc,t2)
#  ifdef MASKING
            my_valu=my_valu*GRID(ng)%umask(i,j)
#  endif
            OCEAN(ng)%ubar_stokes(i,j)=my_valu
          END DO
        END DO
      END IF
# endif
!
! vbar
!
      IF (WESTERN_EDGE) THEN
        DO j=JstrTV+1,JendT
          DO i=IstrT,Istr-1
            ic=i+Nghostpoints+1
            my_valv=fac1*REFINED(ng)%vbarR_west(ic,j,t1)+               &
     &              fac2*REFINED(ng)%vbarR_west(ic,j,t2)
# ifdef MASKING
            my_valv=my_valv*GRID(ng)%vmask(i,j)
# endif
            OCEAN(ng)%vbar(i,j,knew(ng))=my_valv
          END DO
        END DO
      END IF
      IF (EASTERN_EDGE) THEN
        DO j=JstrTV+1,JendT
          DO i=Iend+1,IendT
            ic=i-Iend
            my_valv=fac1*REFINED(ng)%vbarR_east(ic,j,t1)+               &
     &              fac2*REFINED(ng)%vbarR_east(ic,j,t2)
# ifdef MASKING
            my_valv=my_valv*GRID(ng)%vmask(i,j)
# endif
            OCEAN(ng)%vbar(i,j,knew(ng))=my_valv
          END DO
        END DO
      END IF
      IF (SOUTHERN_EDGE) THEN
        DO j=JstrT+1,Jstr
          DO i=IstrT,IendT
            jc=j+Nghostpoints
            my_valv=fac1*REFINED(ng)%vbarR_south(i,jc,t1)+              &
     &              fac2*REFINED(ng)%vbarR_south(i,jc,t2)
# ifdef MASKING
            my_valv=my_valv*GRID(ng)%vmask(i,j)
# endif
            OCEAN(ng)%vbar(i,j,knew(ng))=my_valv
          END DO
        END DO
      END IF
      IF (NORTHERN_EDGE) THEN
        DO j=Jend+1,JendT
          DO i=IstrT,IendT
            jc=j-Jend
            my_valv=fac1*REFINED(ng)%vbarR_north(i,jc,t1)+              &
     &              fac2*REFINED(ng)%vbarR_north(i,jc,t2)
# ifdef MASKING
            my_valv=my_valv*GRID(ng)%vmask(i,j)
# endif
            OCEAN(ng)%vbar(i,j,knew(ng))=my_valv
          END DO
        END DO
      END IF
# ifdef WEC
!
! vbar_stokes
!
      IF (WESTERN_EDGE) THEN
        DO j=JstrTV+1,JendT
          DO i=IstrT,Istr-1
            ic=i+Nghostpoints+1
            my_valv=fac1*REFINED(ng)%vbarsR_west(ic,j,t1)+              &
     &              fac2*REFINED(ng)%vbarsR_west(ic,j,t2)
#  ifdef MASKING
            my_valv=my_valv*GRID(ng)%vmask(i,j)
#  endif
            OCEAN(ng)%vbar_stokes(i,j)=my_valv
          END DO
        END DO
      END IF
      IF (EASTERN_EDGE) THEN
        DO j=JstrTV+1,JendT
          DO i=Iend+1,IendT
            ic=i-Iend
            my_valv=fac1*REFINED(ng)%vbarsR_east(ic,j,t1)+              &
     &              fac2*REFINED(ng)%vbarsR_east(ic,j,t2)
#  ifdef MASKING
            my_valv=my_valv*GRID(ng)%vmask(i,j)
#  endif
            OCEAN(ng)%vbar_stokes(i,j)=my_valv
          END DO
        END DO
      END IF
      IF (SOUTHERN_EDGE) THEN
        DO j=JstrT+1,Jstr
          DO i=IstrT,IendT
            jc=j+Nghostpoints
            my_valv=fac1*REFINED(ng)%vbarsR_south(i,jc,t1)+             &
     &              fac2*REFINED(ng)%vbarsR_south(i,jc,t2)
#  ifdef MASKING
            my_valv=my_valv*GRID(ng)%vmask(i,j)
#  endif
            OCEAN(ng)%vbar_stokes(i,j)=my_valv
          END DO
        END DO
      END IF
      IF (NORTHERN_EDGE) THEN
        DO j=Jend+1,JendT
          DO i=IstrT,IendT
            jc=j-Jend
            my_valv=fac1*REFINED(ng)%vbarsR_north(i,jc,t1)+             &
     &              fac2*REFINED(ng)%vbarsR_north(i,jc,t2)
#  ifdef MASKING
            my_valv=my_valv*GRID(ng)%vmask(i,j)
#  endif
            OCEAN(ng)%vbar_stokes(i,j)=my_valv
          END DO
        END DO
      END IF
# endif
!
# if defined EW_PERIODIC || defined NS_PERIODIC
!
!-----------------------------------------------------------------------
!  Exchange boundary information.
!-----------------------------------------------------------------------
!
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        OCEAN(ng)%zeta(:,:,knew(ng)))
      CALL exchange_u2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        OCEAN(ng)%ubar(:,:,knew(ng)))
      CALL exchange_v2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        OCEAN(ng)%vbar(:,:,knew(ng)))
#  ifdef SOLVE3D
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        COUPLING(ng)%Zt_avg1(:,:))
#  endif
#  ifdef WEC
      CALL exchange_u2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        OCEAN(ng)%ubar_stokes(:,:))
      CALL exchange_v2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        OCEAN(ng)%vbar_stokes(:,:))
#  endif
# endif
# ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, iNLM, 3,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    OCEAN(ng)%zeta(:,:,knew(ng)),                 &
     &                    OCEAN(ng)%ubar(:,:,knew(ng)),                 &
     &                    OCEAN(ng)%vbar(:,:,knew(ng)))
#  ifdef SOLVE3D
      CALL mp_exchange2d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    COUPLING(ng)%Zt_avg1(:,:))
#  endif
#  ifdef WEC
      CALL mp_exchange2d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    OCEAN(ng)%ubar_stokes(:,:),                   &
     &                    OCEAN(ng)%vbar_stokes(:,:))
#  endif
# endif

# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn off input data time wall clock.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, iNLM, 3)
# endif
      RETURN
      END SUBROUTINE set_2dparent_data_tile
#endif
      END MODULE set_2dparent_data_mod
