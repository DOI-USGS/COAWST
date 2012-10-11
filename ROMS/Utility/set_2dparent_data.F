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
     &                           LBip, UBip, LBjp, UBjp,                &
     &                           IminS, ImaxS, JminS, JmaxS)
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 4)
# endif
      RETURN
      END SUBROUTINE set_2dparent_data
!
!***********************************************************************
      SUBROUTINE set_2dparent_data_tile (ng, tile,                      &
     &                                   LBi, UBi, LBj, UBj,            &
     &                                   LBip, UBip, LBjp, UBjp,        &
     &                                   IminS, ImaxS, JminS, JmaxS)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
      USE mod_ncparam
      USE mod_scalars
      USE mod_stepping
      USE mod_parallel
      USE mod_ocean
      USE mod_refined
      USE mod_scalars
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
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
!  Local variable declarations.
!
      logical, dimension(3) :: update =                                 &
     &         (/ .FALSE., .FALSE., .FALSE. /)

# ifdef OBC
      integer :: ILB, IUB, JLB, JUB
# endif
      integer :: i, j, ngp, ic, jc, t1, t2
      real(r8) :: my_val, fac, fac1, fac2, fac3, ref_scale
      real(r8) :: my_valu, my_valv, cff1, cff2
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
!  Set refinement scale.
!
      ref_scale=REAL(REFINED(ng)%ref_scale,r8)
!
!  Compute time interpolation factors.
!
# ifdef SOLVE3D
      fac1=time(ngp)+dt(ngp)-time(ng)-dt(ng)
      fac2=time(ng)-time(ngp)+dt(ng)
      fac=1.0_r8/(fac1+fac2)
      fac1=fac*fac1
      fac2=fac*fac2
      IF (iic(ngp).eq.ntstart(ngp)) THEN
        t1=1
        t2=2
      ELSE
        t1=nstp(ngp)
        t2=nnew(ngp)
      END IF
# else
      fac1=time(ngp)+dt(ngp)-time(ng)-dt(ng)
      fac2=time(ng)-time(ngp)+dt(ng)
      fac=1.0_r8/(fac1+fac2)
      fac1=fac*fac1
      fac2=fac*fac2
      IF (iic(ngp).eq.ntstart(ngp)) THEN
        t1=1
        t2=2
      ELSE
        t1=kstp(ngp)
        t2=3-t1
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
            OCEAN(ng)%zeta(i,j,1)=my_val
            OCEAN(ng)%zeta(i,j,2)=my_val
            OCEAN(ng)%zeta(i,j,3)=my_val
# ifdef SOLVE3D
!            COUPLING(ng)%Zt_avg1(i,j)=my_val
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
            OCEAN(ng)%zeta(i,j,1)=my_val
            OCEAN(ng)%zeta(i,j,2)=my_val
            OCEAN(ng)%zeta(i,j,3)=my_val
# ifdef SOLVE3D
!            COUPLING(ng)%Zt_avg1(i,j)=my_val
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
            OCEAN(ng)%zeta(i,j,1)=my_val
            OCEAN(ng)%zeta(i,j,2)=my_val
            OCEAN(ng)%zeta(i,j,3)=my_val
# ifdef SOLVE3D
!            COUPLING(ng)%Zt_avg1(i,j)=my_val
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
            OCEAN(ng)%zeta(i,j,1)=my_val
            OCEAN(ng)%zeta(i,j,2)=my_val
            OCEAN(ng)%zeta(i,j,3)=my_val
# ifdef SOLVE3D
!           COUPLING(ng)%Zt_avg1(i,j)=my_val
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
# ifdef WET_DRY
            my_valu=my_valu*GRID(ng)%umask_wet(i,j)
# endif
            OCEAN(ng)%ubar(i,j,1)=my_valu
            OCEAN(ng)%ubar(i,j,2)=my_valu
            OCEAN(ng)%ubar(i,j,3)=my_valu
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
# ifdef WET_DRY
            my_valu=my_valu*GRID(ng)%umask_wet(i,j)
# endif
            OCEAN(ng)%ubar(i,j,1)=my_valu
            OCEAN(ng)%ubar(i,j,2)=my_valu
            OCEAN(ng)%ubar(i,j,3)=my_valu
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
# ifdef WET_DRY
            my_valu=my_valu*GRID(ng)%umask_wet(i,j)
# endif
            OCEAN(ng)%ubar(i,j,1)=my_valu
            OCEAN(ng)%ubar(i,j,2)=my_valu
            OCEAN(ng)%ubar(i,j,3)=my_valu
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
# ifdef WET_DRY
            my_valu=my_valu*GRID(ng)%umask_wet(i,j)
# endif
            OCEAN(ng)%ubar(i,j,1)=my_valu
            OCEAN(ng)%ubar(i,j,2)=my_valu
            OCEAN(ng)%ubar(i,j,3)=my_valu
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
#  ifdef WET_DRY
            my_valu=my_valu*GRID(ng)%umask_wet(i,j)
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
#  ifdef WET_DRY
            my_valu=my_valu*GRID(ng)%umask_wet(i,j)
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
#  ifdef WET_DRY
            my_valu=my_valu*GRID(ng)%umask_wet(i,j)
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
#  ifdef WET_DRY
            my_valu=my_valu*GRID(ng)%umask_wet(i,j)
#  endif
            OCEAN(ng)%ubar_stokes(i,j)=my_valu
          END DO
        END DO
      END IF
# endif
!
! DU_avg2
!
      IF (WESTERN_EDGE) THEN
        DO j=Jstr,Jend
!
! To ensure mass flux, set ubar at west edge.
!
          cff1=0.5_r8*GRID(ng)%on_u(Istr,j)
          cff2=cff1*((GRID(ng)%h(Istr-1,j)+                             &
     &                OCEAN(ng)%zeta(Istr-1,j,indx1(ng))+               &
     &                GRID(ng)%h(Istr  ,j)+                             &
     &                OCEAN(ng)%zeta(Istr  ,j,indx1(ng))))
!         my_val=REFINED(ng)%DUavg2R_west(j)/(ref_scale*cff2)
          my_val=REFINED(ng)%DUavg2R_west(j)/cff2
# ifdef WEC
          my_val=my_val-OCEAN(ng)%ubar_stokes(Istr,j)
# endif
# ifdef MASKING
          my_val=my_val*GRID(ng)%umask(Istr,j)
# endif
# ifdef WET_DRY
          my_val=my_val*GRID(ng)%umask_wet(Istr,j)
# endif
          OCEAN(ng)%ubar(Istr,j,indx1(ng))=my_val
        END DO
      END IF
      IF (EASTERN_EDGE) THEN
        DO j=Jstr,Jend
!
! To ensure mass flux, set ubar at east edge.
!
          cff2=0.5_r8*GRID(ng)%on_u(Iend+1,j)
          cff2=cff2*((GRID(ng)%h(Iend+1,j)+                             &
     &                OCEAN(ng)%zeta(Iend+1,j,indx1(ng))+               &
     &                GRID(ng)%h(Iend  ,j)+                             &
     &                OCEAN(ng)%zeta(Iend  ,j,indx1(ng))))
!         my_val=REFINED(ng)%DUavg2R_east(j)/(ref_scale*cff2)
          my_val=REFINED(ng)%DUavg2R_east(j)/cff2
# ifdef WEC
          my_val=my_val-OCEAN(ng)%ubar_stokes(Iend+1,j)
# endif
# ifdef MASKING
          my_val=my_val*GRID(ng)%umask(Iend+1,j)
# endif
# ifdef WET_DRY
          my_val=my_val*GRID(ng)%umask_wet(Iend+1,j)
# endif
          OCEAN(ng)%ubar(Iend+1,j,indx1(ng))=my_val
        END DO
      END IF
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
# ifdef WET_DRY
            my_valv=my_valv*GRID(ng)%vmask_wet(i,j)
# endif
            OCEAN(ng)%vbar(i,j,1)=my_valv
            OCEAN(ng)%vbar(i,j,2)=my_valv
            OCEAN(ng)%vbar(i,j,3)=my_valv
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
# ifdef WET_DRY
            my_valv=my_valv*GRID(ng)%vmask_wet(i,j)
# endif
            OCEAN(ng)%vbar(i,j,1)=my_valv
            OCEAN(ng)%vbar(i,j,2)=my_valv
            OCEAN(ng)%vbar(i,j,3)=my_valv
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
# ifdef WET_DRY
            my_valv=my_valv*GRID(ng)%vmask_wet(i,j)
# endif
            OCEAN(ng)%vbar(i,j,1)=my_valv
            OCEAN(ng)%vbar(i,j,2)=my_valv
            OCEAN(ng)%vbar(i,j,3)=my_valv
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
# ifdef WET_DRY
            my_valv=my_valv*GRID(ng)%vmask_wet(i,j)
# endif
            OCEAN(ng)%vbar(i,j,1)=my_valv
            OCEAN(ng)%vbar(i,j,2)=my_valv
            OCEAN(ng)%vbar(i,j,3)=my_valv
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
#  ifdef WET_DRY
            my_valv=my_valv*GRID(ng)%vmask_wet(i,j)
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
#  ifdef WET_DRY
            my_valv=my_valv*GRID(ng)%vmask_wet(i,j)
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
#  ifdef WET_DRY
            my_valv=my_valv*GRID(ng)%vmask_wet(i,j)
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
#  ifdef WET_DRY
            my_valv=my_valv*GRID(ng)%vmask_wet(i,j)
#  endif
            OCEAN(ng)%vbar_stokes(i,j)=my_valv
          END DO
        END DO
      END IF
# endif
!
! DV_avg2
!
      IF (SOUTHERN_EDGE) THEN
        DO i=Istr,Iend
!
! To ensure mass flux, set vbar at southern edge.
!
          cff1=0.5_r8*GRID(ng)%om_v(i,Jstr)
          cff2=cff1*((GRID(ng)%h(i,Jstr-1)+                             &
     &                OCEAN(ng)%zeta(i,Jstr-1,indx1(ng))+               &
     &                GRID(ng)%h(i,Jstr  )+                             &
     &                OCEAN(ng)%zeta(i,Jstr  ,indx1(ng))))
!         my_val=REFINED(ng)%DVavg2R_south(i)/(ref_scale*cff2)
          my_val=REFINED(ng)%DVavg2R_south(i)/cff2
# ifdef WEC
          my_val=my_val-OCEAN(ng)%vbar_stokes(i,Jstr)
# endif
# ifdef MASKING
          my_val=my_val*GRID(ng)%vmask(i,Jstr)
# endif
# ifdef WET_DRY
          my_val=my_val*GRID(ng)%vmask_wet(i,Jstr)
# endif
          OCEAN(ng)%vbar(i,Jstr,indx1(ng))=my_val
        END DO
      END IF
      IF (NORTHERN_EDGE) THEN
        DO i=Istr,Iend
!
! To ensure mass flux, set vbar at southern edge.
!
          cff1=0.5_r8*GRID(ng)%om_v(i,Jend+1)
          cff2=cff1*((GRID(ng)%h(i,Jend+1)+                             &
     &                OCEAN(ng)%zeta(i,Jend+1,indx1(ng))+               &
     &                GRID(ng)%h(i,Jend  )+                             &
     &                OCEAN(ng)%zeta(i,Jend  ,indx1(ng))))
!         my_val=REFINED(ng)%DVavg2R_north(i)/(ref_scale*cff2)
          my_val=REFINED(ng)%DVavg2R_north(i)/cff2
# ifdef WEC
          my_val=my_val-OCEAN(ng)%vbar_stokes(i,Jend+1)
# endif
# ifdef MASKING
          my_val=my_val*GRID(ng)%vmask(i,Jend+1)
# endif
# ifdef WET_DRY
          my_val=my_val*GRID(ng)%vmask_wet(i,Jend+1)
# endif
          OCEAN(ng)%vbar(i,Jend+1,indx1(ng))=my_val
        END DO
      END IF
!
# ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, iNLM, 3,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    OCEAN(ng)%zeta(:,:,1),                        &
     &                    OCEAN(ng)%zeta(:,:,2),                        &
     &                    OCEAN(ng)%zeta(:,:,3))
      CALL mp_exchange2d (ng, tile, iNLM, 3,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    OCEAN(ng)%ubar(:,:,1),                        &
     &                    OCEAN(ng)%ubar(:,:,2),                        &
     &                    OCEAN(ng)%ubar(:,:,3))
      CALL mp_exchange2d (ng, tile, iNLM, 3,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    OCEAN(ng)%vbar(:,:,1),                        &
     &                    OCEAN(ng)%vbar(:,:,2),                        &
     &                    OCEAN(ng)%vbar(:,:,3))
!#  ifdef SOLVE3D
!      CALL mp_exchange2d (ng, tile, iNLM, 1,                            &
!     &                    LBi, UBi, LBj, UBj,                           &
!     &                    NghostPoints, EWperiodic, NSperiodic,         &
!     &                    COUPLING(ng)%Zt_avg1(:,:))
!#  endif
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
