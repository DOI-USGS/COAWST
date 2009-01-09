#include "cppdefs.h"
#undef get_DUon
#undef get_DVom

      MODULE set_2dparent_data_mod
#ifdef REFINED_GRID
!
!svn $Id: set_data.F 1405 2008-02-20 01:44:07Z jcwarner $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine acquires date from a parent grid.                      !
!  If there is more than one time-record,                              !
!  data  is loaded  into global two-time record arrays.  The actual    !
!  interpolation is carried elsewhere.                                 !
!                                                                      !
!  Currently,  this routine is only executed in serial mode by the     !
!  main thread.                                                        !
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
      USE mod_scalars
      USE mod_coupling
      USE interpolate_mod

# ifdef DISTRIBUTE
      USE distribute_mod, ONLY : mp_bcasti, mp_gather
# endif
# ifdef DISTRIBUTE
      USE distribute_mod, ONLY : mp_bcasti, mp_gather
# endif
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
      real(r8) :: my_val, fac, fac1, fac2, fac3
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

!      fac1=0.0_r8
!      fac2=1.0_r8
      IF (iic(1).eq.ntstart(1)) THEN
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

!scale fac for DUon and DVom
      fac3=1.0_r8/nrefined(ng)

!       IF (myrank.eq.mymaster) then
!       write(*,*) 'set2dp ',t1,t2,fac1,fac2,time(ngp),time(ng),dt(ngp)
!       end if
!
! zeta
!
      IF (WESTERN_EDGE) THEN
        DO j=JstrT,JendT
          DO i=IstrT,Istr-1
!         DO i=IstrT,Istr-2
            ic=i+Nghostpoints+1
            my_val=fac1*BOUNDARY(ng)%zetaR_west(ic,j,t1)+               &
     &             fac2*BOUNDARY(ng)%zetaR_west(ic,j,t2)
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
            COUPLING(ng)%Zt_avg1(i,j)=my_val
# endif
          END DO
        END DO
      END IF
      IF (EASTERN_EDGE) THEN
        DO j=JstrT,JendT
          DO i=Iend+1,IendT
!         DO i=Iend+2,IendT
            ic=i-Iend                                                   !+1
            my_val=fac1*BOUNDARY(ng)%zetaR_east(ic,j,t1)+               &
     &             fac2*BOUNDARY(ng)%zetaR_east(ic,j,t2)
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
            COUPLING(ng)%Zt_avg1(i,j)=my_val
# endif
          END DO
        END DO
      END IF
      IF (SOUTHERN_EDGE) THEN
        DO j=JstrT,Jstr-1
!       DO j=JstrT,Jstr-2
          DO i=IstrT,IendT
            jc=j+Nghostpoints+1
            my_val=fac1*BOUNDARY(ng)%zetaR_south(i,jc,t1)+               &
     &             fac2*BOUNDARY(ng)%zetaR_south(i,jc,t2)
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
            COUPLING(ng)%Zt_avg1(i,j)=my_val
# endif
          END DO
        END DO
      END IF
      IF (NORTHERN_EDGE) THEN
        DO j=Jend+1,JendT
!       DO j=Jend+2,JendT
          DO i=IstrT,IendT
            jc=j-Jend                                                 !+1
            my_val=fac1*BOUNDARY(ng)%zetaR_north(i,jc,t1)+               &
     &             fac2*BOUNDARY(ng)%zetaR_north(i,jc,t2)
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
!         DO i=IstrT+1,Istr-1
            ic=i+Nghostpoints
            my_valu=fac1*BOUNDARY(ng)%ubarR_west(ic,j,t1)+              &
     &              fac2*BOUNDARY(ng)%ubarR_west(ic,j,t2)
# ifdef MASKING
            my_valu=my_valu*GRID(ng)%umask(i,j)
# endif
            my_valh1=fac1*BOUNDARY(ng)%hR_west(ic,j,t1)+                &
     &               fac2*BOUNDARY(ng)%hR_west(ic,j,t2)+                &
     &               fac1*BOUNDARY(ng)%zetaR_west(ic,j,t1)+             &
     &               fac2*BOUNDARY(ng)%zetaR_west(ic,j,t2)
            my_valh2=fac1*BOUNDARY(ng)%hR_west(ic+1,j,t1)+              &
     &               fac2*BOUNDARY(ng)%hR_west(ic+1,j,t2)+              &
     &               fac1*BOUNDARY(ng)%zetaR_west(ic+1,j,t1)+           &
     &               fac2*BOUNDARY(ng)%zetaR_west(ic+1,j,t2)
            cff1=0.5_r8*(GRID(ng)%h(i-1,j)+OCEAN(ng)%zeta(i-1,j,1)+     &
     &                   GRID(ng)%h(i  ,j)+OCEAN(ng)%zeta(i  ,j,1))
            my_val=my_valu*0.5_r8*(my_valh1+my_valh2)/cff1
            OCEAN(ng)%ubar(i,j,1)=my_val
            OCEAN(ng)%ubar(i,j,2)=my_val
            OCEAN(ng)%ubar(i,j,3)=my_val
# ifdef SOLVE3D
!            COUPLING(ng)%DU_avg1(i,j)=my_val
!            COUPLING(ng)%DU_avg2(i,j)=my_val
# endif
          END DO
        END DO
      END IF
      IF (EASTERN_EDGE) THEN
        DO j=JstrT,JendT
          DO i=Iend+1,IendT
!         DO i=Iend+2,IendT
            ic=i-Iend
            my_valu=fac1*BOUNDARY(ng)%ubarR_east(ic,j,t1)+              &
     &              fac2*BOUNDARY(ng)%ubarR_east(ic,j,t2)
# ifdef MASKING
            my_valu=my_valu*GRID(ng)%umask(i,j)
# endif
            my_valh1=fac1*BOUNDARY(ng)%hR_east(ic,j,t1)+                &
     &               fac2*BOUNDARY(ng)%hR_east(ic,j,t2)+                &
     &               fac1*BOUNDARY(ng)%zetaR_east(ic,j,t1)+             &
     &               fac2*BOUNDARY(ng)%zetaR_east(ic,j,t2)
            my_valh2=fac1*BOUNDARY(ng)%hR_east(ic+1,j,t1)+              &
     &               fac2*BOUNDARY(ng)%hR_east(ic+1,j,t2)+              &
     &               fac1*BOUNDARY(ng)%zetaR_east(ic+1,j,t1)+           &
     &               fac2*BOUNDARY(ng)%zetaR_east(ic+1,j,t2)
            cff1=0.5_r8*(GRID(ng)%h(i-1,j)+OCEAN(ng)%zeta(i-1,j,1)+     &
     &                   GRID(ng)%h(i  ,j)+OCEAN(ng)%zeta(i  ,j,1))
            my_val=my_valu*0.5_r8*(my_valh1+my_valh2)/cff1
            OCEAN(ng)%ubar(i,j,1)=my_val
            OCEAN(ng)%ubar(i,j,2)=my_val
            OCEAN(ng)%ubar(i,j,3)=my_val
# ifdef SOLVE3D
!            COUPLING(ng)%DU_avg1(i,j)=my_val
!            COUPLING(ng)%DU_avg2(i,j)=my_val
# endif
          END DO
        END DO
      END IF
      IF (SOUTHERN_EDGE) THEN
        DO j=JstrT,Jstr-1
!       DO j=JstrT,Jstr-2
          DO i=IstrTU+1,IendT
            jc=j+Nghostpoints+1
            my_valu=fac1*BOUNDARY(ng)%ubarR_south(i,jc,t1)+             &
     &              fac2*BOUNDARY(ng)%ubarR_south(i,jc,t2)
# ifdef MASKING
            my_valu=my_valu*GRID(ng)%umask(i,j)
# endif
            my_valh1=fac1*BOUNDARY(ng)%hR_south(i,jc,t1)+               &
     &               fac2*BOUNDARY(ng)%hR_south(i,jc,t2)+               &
     &               fac1*BOUNDARY(ng)%zetaR_south(i,jc,t1)+            &
     &               fac2*BOUNDARY(ng)%zetaR_south(i,jc,t2)
            my_valh2=fac1*BOUNDARY(ng)%hR_south(i-1,jc,t1)+             &
     &               fac2*BOUNDARY(ng)%hR_south(i-1,jc,t2)+             &
     &               fac1*BOUNDARY(ng)%zetaR_south(i-1,jc,t1)+          &
     &               fac2*BOUNDARY(ng)%zetaR_south(i-1,jc,t2)
            cff1=0.5_r8*(GRID(ng)%h(i-1,j)+OCEAN(ng)%zeta(i-1,j,1)+     &
     &                   GRID(ng)%h(i  ,j)+OCEAN(ng)%zeta(i  ,j,1))
            my_val=my_valu*0.5_r8*(my_valh1+my_valh2)/cff1
            OCEAN(ng)%ubar(i,j,1)=my_val
            OCEAN(ng)%ubar(i,j,2)=my_val
            OCEAN(ng)%ubar(i,j,3)=my_val
# ifdef SOLVE3D
!            COUPLING(ng)%DU_avg1(i,j)=my_val
!            COUPLING(ng)%DU_avg2(i,j)=my_val
# endif
          END DO
        END DO
      END IF
      IF (NORTHERN_EDGE) THEN
        DO j=Jend+1,JendT
!       DO j=Jend+2,JendT
          DO i=IstrTU+1,IendT
            jc=j-Jend
            my_valu=fac1*BOUNDARY(ng)%ubarR_north(i,jc,t1)+             &
     &             fac2*BOUNDARY(ng)%ubarR_north(i,jc,t2)
# ifdef MASKING
            my_valu=my_valu*GRID(ng)%umask(i,j)
# endif
            my_valh1=fac1*BOUNDARY(ng)%hR_north(i,jc+1,t1)+             &
     &               fac2*BOUNDARY(ng)%hR_north(i,jc+1,t2)+             &
     &               fac1*BOUNDARY(ng)%zetaR_north(i,jc+1,t1)+          &
     &               fac2*BOUNDARY(ng)%zetaR_north(i,jc+1,t2)
            my_valh2=fac1*BOUNDARY(ng)%hR_north(i-1,jc+1,t1)+           &
     &               fac2*BOUNDARY(ng)%hR_north(i-1,jc+1,t2)+           &
     &               fac1*BOUNDARY(ng)%zetaR_north(i-1,jc+1,t1)+        &
     &               fac2*BOUNDARY(ng)%zetaR_north(i-1,jc+1,t2)
            cff1=0.5_r8*(GRID(ng)%h(i-1,j)+OCEAN(ng)%zeta(i-1,j,1)+     &
     &                   GRID(ng)%h(i  ,j)+OCEAN(ng)%zeta(i  ,j,1))
            my_val=my_valu*0.5_r8*(my_valh1+my_valh2)/cff1
            OCEAN(ng)%ubar(i,j,1)=my_val
            OCEAN(ng)%ubar(i,j,2)=my_val
            OCEAN(ng)%ubar(i,j,3)=my_val
# ifdef SOLVE3D
!            COUPLING(ng)%DU_avg1(i,j)=my_val
!            COUPLING(ng)%DU_avg2(i,j)=my_val
# endif
          END DO
        END DO
      END IF

!
! DUon
!
# ifdef get_DUon
      IF (WESTERN_EDGE) THEN
        DO j=JstrT,JendT
          DO i=IstrT+1,Istr
            ic=i+Nghostpoints
            my_val=fac1*BOUNDARY(ng)%DUonR_west(ic,j,t1)+               &
     &             fac2*BOUNDARY(ng)%DUonR_west(ic,j,t2)
            my_val=my_val*fac3
#  ifdef MASKING
            my_val=my_val*GRID(ng)%umask(i,j)
#  endif
!            OCEAN(ng)%DUon(i,j)=my_val
          END DO
        END DO
      END IF
      IF (EASTERN_EDGE) THEN
        DO j=JstrT,JendT
          DO i=Iend+1,IendT
            ic=i-Iend
            my_val=fac1*BOUNDARY(ng)%DUonR_east(ic,j,t1)+               &
     &             fac2*BOUNDARY(ng)%DUonR_east(ic,j,t2)
            my_val=my_val*fac3
#  ifdef MASKING
            my_val=my_val*GRID(ng)%umask(i,j)
#  endif
!            OCEAN(ng)%DUon(i,j)=my_val
          END DO
        END DO
      END IF
      IF (SOUTHERN_EDGE) THEN
        DO j=JstrT,Jstr-1
          DO i=IstrT,IendT
            jc=j+Nghostpoints+1
            my_val=fac1*BOUNDARY(ng)%DUonR_south(i,jc,t1)+              &
     &             fac2*BOUNDARY(ng)%DUonR_south(i,jc,t2)
            my_val=my_val*fac3
#  ifdef MASKING
            my_val=my_val*GRID(ng)%umask(i,j)
#  endif
!            OCEAN(ng)%DUon(i,j)=my_val
          END DO
        END DO
      END IF
      IF (NORTHERN_EDGE) THEN
        DO j=Jend+1,JendT
          DO i=IstrT,IendT
            jc=j-Jend
            my_val=fac1*BOUNDARY(ng)%DUonR_north(i,jc,t1)+               &
     &             fac2*BOUNDARY(ng)%DUonR_north(i,jc,t2)
            my_val=my_val*fac3
#  ifdef MASKING
            my_val=my_val*GRID(ng)%umask(i,j)
#  endif
!            OCEAN(ng)%DUon(i,j)=my_val
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
!         DO i=IstrT,Istr-2
            ic=i+Nghostpoints+1
            my_valv=fac1*BOUNDARY(ng)%vbarR_west(ic,j,t1)+              &
     &              fac2*BOUNDARY(ng)%vbarR_west(ic,j,t2)
# ifdef MASKING
            my_valv=my_valv*GRID(ng)%vmask(i,j)
# endif
            my_valh1=fac1*BOUNDARY(ng)%hR_west(ic,j,t1)+                &
     &               fac2*BOUNDARY(ng)%hR_west(ic,j,t2)+                &
     &               fac1*BOUNDARY(ng)%zetaR_west(ic,j,t1)+             &
     &               fac2*BOUNDARY(ng)%zetaR_west(ic,j,t2)
            my_valh2=fac1*BOUNDARY(ng)%hR_west(ic,j-1,t1)+              &
     &               fac2*BOUNDARY(ng)%hR_west(ic,j-1,t2)+              &
     &               fac1*BOUNDARY(ng)%zetaR_west(ic,j-1,t1)+           &
     &               fac2*BOUNDARY(ng)%zetaR_west(ic,j-1,t2)
            cff1=0.5_r8*(GRID(ng)%h(i,j  )+OCEAN(ng)%zeta(i,j  ,1)+     &
     &                   GRID(ng)%h(i,j-1)+OCEAN(ng)%zeta(i,j  ,1))
            my_val=my_valv*0.5_r8*(my_valh1+my_valh2)/cff1
            OCEAN(ng)%vbar(i,j,1)=my_val
            OCEAN(ng)%vbar(i,j,2)=my_val
            OCEAN(ng)%vbar(i,j,3)=my_val
# ifdef SOLVE3D
!            COUPLING(ng)%DV_avg1(i,j)=my_val
!            COUPLING(ng)%DV_avg2(i,j)=my_val
# endif
          END DO
        END DO
      END IF
      IF (EASTERN_EDGE) THEN
        DO j=JstrTV+1,JendT
          DO i=Iend+1,IendT
!         DO i=Iend+2,IendT
            ic=i-Iend
            my_valv=fac1*BOUNDARY(ng)%vbarR_east(ic,j,t1)+              &
     &              fac2*BOUNDARY(ng)%vbarR_east(ic,j,t2)
# ifdef MASKING
            my_valv=my_valv*GRID(ng)%vmask(i,j)
# endif
            my_valh1=fac1*BOUNDARY(ng)%hR_east(ic+1,j,t1)+              &
     &               fac2*BOUNDARY(ng)%hR_east(ic+1,j,t2)+              &
     &               fac1*BOUNDARY(ng)%zetaR_east(ic+1,j,t1)+           &
     &               fac2*BOUNDARY(ng)%zetaR_east(ic+1,j,t2)
            my_valh2=fac1*BOUNDARY(ng)%hR_east(ic+1,j-1,t1)+            &
     &               fac2*BOUNDARY(ng)%hR_east(ic+1,j-1,t2)+            &
     &               fac1*BOUNDARY(ng)%zetaR_east(ic+1,j-1,t1)+         &
     &               fac2*BOUNDARY(ng)%zetaR_east(ic+1,j-1,t2)
            cff1=0.5_r8*(GRID(ng)%h(i,j  )+OCEAN(ng)%zeta(i,j  ,1)+     &
     &                   GRID(ng)%h(i,j-1)+OCEAN(ng)%zeta(i,j  ,1))
            my_val=my_valv*0.5_r8*(my_valh1+my_valh2)/cff1
            OCEAN(ng)%vbar(i,j,1)=my_val
            OCEAN(ng)%vbar(i,j,2)=my_val
            OCEAN(ng)%vbar(i,j,3)=my_val
# ifdef SOLVE3D
!            COUPLING(ng)%DV_avg1(i,j)=my_val
!            COUPLING(ng)%DV_avg2(i,j)=my_val
# endif
          END DO
        END DO
      END IF
      IF (SOUTHERN_EDGE) THEN
        DO j=JstrT+1,Jstr
!       DO j=JstrT+1,Jstr-1
          DO i=IstrT,IendT
            jc=j+Nghostpoints
            my_valv=fac1*BOUNDARY(ng)%vbarR_south(i,jc,t1)+             &
     &              fac2*BOUNDARY(ng)%vbarR_south(i,jc,t2)
# ifdef MASKING
            my_valv=my_valv*GRID(ng)%vmask(i,j)
# endif
            my_valh1=fac1*BOUNDARY(ng)%hR_south(i,jc,t1)+                &
     &               fac2*BOUNDARY(ng)%hR_south(i,jc,t2)+                &
     &               fac1*BOUNDARY(ng)%zetaR_south(i,jc,t1)+             &
     &               fac2*BOUNDARY(ng)%zetaR_south(i,jc,t2)
            my_valh2=fac1*BOUNDARY(ng)%hR_south(i,jc+1,t1)+              &
     &               fac2*BOUNDARY(ng)%hR_south(i,jc+1,t2)+              &
     &               fac1*BOUNDARY(ng)%zetaR_south(i,jc+1,t1)+           &
     &               fac2*BOUNDARY(ng)%zetaR_south(i,jc+1,t2)
            cff1=0.5_r8*(GRID(ng)%h(i,j-1)+OCEAN(ng)%zeta(i,j-1,1)+      &
     &                   GRID(ng)%h(i,j  )+OCEAN(ng)%zeta(i,j  ,1))
            my_val=my_valv*0.5_r8*(my_valh1+my_valh2)/cff1
            OCEAN(ng)%vbar(i,j,1)=my_val
            OCEAN(ng)%vbar(i,j,2)=my_val
            OCEAN(ng)%vbar(i,j,3)=my_val
# ifdef SOLVE3D
!            COUPLING(ng)%DV_avg1(i,j)=my_val
!            COUPLING(ng)%DV_avg2(i,j)=my_val
# endif
          END DO
        END DO
      END IF
      IF (NORTHERN_EDGE) THEN
        DO j=Jend+1,JendT
!       DO j=Jend+2,JendT
          DO i=IstrT,IendT
            jc=j-Jend
            my_valv=fac1*BOUNDARY(ng)%vbarR_north(i,jc,t1)+             &
     &              fac2*BOUNDARY(ng)%vbarR_north(i,jc,t2)
# ifdef MASKING
            my_valv=my_valv*GRID(ng)%vmask(i,j)
# endif
            my_valh1=fac1*BOUNDARY(ng)%hR_north(i,jc,t1)+               &
     &               fac2*BOUNDARY(ng)%hR_north(i,jc,t2)+               &
     &               fac1*BOUNDARY(ng)%zetaR_north(i,jc,t1)+            &
     &               fac2*BOUNDARY(ng)%zetaR_north(i,jc,t2)
            my_valh2=fac1*BOUNDARY(ng)%hR_north(i,jc+1,t1)+             &
     &               fac2*BOUNDARY(ng)%hR_north(i,jc+1,t2)+             &
     &               fac1*BOUNDARY(ng)%zetaR_north(i,jc+1,t1)+          &
     &               fac2*BOUNDARY(ng)%zetaR_north(i,jc+1,t2)
            cff1=0.5_r8*(GRID(ng)%h(i,j-1)+OCEAN(ng)%zeta(i,j-1,1)+     &
     &                   GRID(ng)%h(i,j  )+OCEAN(ng)%zeta(i,j  ,1))
            my_val=my_valv*0.5_r8*(my_valh1+my_valh2)/cff1
            OCEAN(ng)%vbar(i,j,1)=my_val
            OCEAN(ng)%vbar(i,j,2)=my_val
            OCEAN(ng)%vbar(i,j,3)=my_val
# ifdef SOLVE3D
!            COUPLING(ng)%DV_avg1(i,j)=my_val
!            COUPLING(ng)%DV_avg2(i,j)=my_val
# endif
          END DO
        END DO
      END IF

# ifdef get_DVom
!
! DVom
!
      IF (WESTERN_EDGE) THEN
        DO j=JstrT,JendT
          DO i=IstrT,Istr-1
            ic=i+Nghostpoints+1
            my_val=fac1*BOUNDARY(ng)%DVomR_west(ic,j,t1)+               &
     &             fac2*BOUNDARY(ng)%DVomR_west(ic,j,t2)
            my_val=my_val*fac3
#  ifdef MASKING
            my_val=my_val*GRID(ng)%vmask(i,j)
#  endif
!            OCEAN(ng)%DVom(i,j)=my_val
          END DO
        END DO
      END IF
      IF (EASTERN_EDGE) THEN
        DO j=JstrT,JendT
          DO i=Iend+1,IendT
            ic=i-Iend
            my_val=fac1*BOUNDARY(ng)%DVomR_east(ic,j,t1)+               &
     &             fac2*BOUNDARY(ng)%DVomR_east(ic,j,t2)
            my_val=my_val*fac3
#  ifdef MASKING
            my_val=my_val*GRID(ng)%vmask(i,j)
#  endif
!            OCEAN(ng)%DVom(i,j)=my_val
          END DO
        END DO
      END IF
      IF (SOUTHERN_EDGE) THEN
        DO j=JstrT+1,Jstr
          DO i=IstrT,IendT
            jc=j+Nghostpoints
            my_val=fac1*BOUNDARY(ng)%DVomR_south(i,jc,t1)+               &
     &             fac2*BOUNDARY(ng)%DVomR_south(i,jc,t2)
            my_val=my_val*fac3
#  ifdef MASKING
            my_val=my_val*GRID(ng)%vmask(i,j)
#  endif
!            OCEAN(ng)%DVom(i,j)=my_val
          END DO
        END DO
      END IF
      IF (NORTHERN_EDGE) THEN
        DO j=Jend+1,JendT
          DO i=IstrT,IendT
            jc=j-Jend
            my_val=fac1*BOUNDARY(ng)%DVomR_north(i,jc,t1)+               &
     &             fac2*BOUNDARY(ng)%DVomR_north(i,jc,t2)
            my_val=my_val*fac3
#  ifdef MASKING
            my_val=my_val*GRID(ng)%vmask(i,j)
#  endif
!            OCEAN(ng)%DVom(i,j)=my_val
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
      IF (ng.eq.1) THEN
        CALL exchange_u2d_tile (ng, tile,                               &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        OCEAN(ng)%ubar(:,:,1))
        CALL exchange_v2d_tile (ng, tile,                               &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        OCEAN(ng)%vbar(:,:,1))
        CALL exchange_u2d_tile (ng, tile,                               &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        OCEAN(ng)%ubar(:,:,2))
        CALL exchange_v2d_tile (ng, tile,                               &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        OCEAN(ng)%vbar(:,:,2))
        CALL exchange_u2d_tile (ng, tile,                               &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        OCEAN(ng)%ubar(:,:,3))
        CALL exchange_v2d_tile (ng, tile,                               &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        OCEAN(ng)%vbar(:,:,3))
        CALL exchange_v2d_tile (ng, tile,                               &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        OCEAN(ng)%zeta(:,:,1))
        CALL exchange_v2d_tile (ng, tile,                               &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        OCEAN(ng)%zeta(:,:,2))
        CALL exchange_v2d_tile (ng, tile,                               &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        OCEAN(ng)%zeta(:,:,3))
      END IF
#  endif
#  ifdef DISTRIBUTE
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
#  endif

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
