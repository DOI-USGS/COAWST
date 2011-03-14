#include "cppdefs.h"
      MODULE set_3dparent_data_mod
#ifdef SOLVE3D
# ifdef REFINED_GRID
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
      PUBLIC  :: set_3dparent_data, set_3dparent_data_tile

      CONTAINS
!
!***********************************************************************
      SUBROUTINE set_3dparent_data (ng, tile)
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
#  include "tile.h"

      ngp=myparent(ng)
      LBip=BOUNDS(ngp)%LBi(tile)
      UBip=BOUNDS(ngp)%UBi(tile)
      LBjp=BOUNDS(ngp)%LBj(tile)
      UBjp=BOUNDS(ngp)%UBj(tile)
!
#  ifdef PROFILE
      CALL wclock_on (ng, iNLM, 4)
#  endif
      CALL set_3dparent_data_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             LBip, UBip, LBjp, UBjp)
#  ifdef PROFILE
      CALL wclock_off (ng, iNLM, 4)
#  endif
      RETURN
      END SUBROUTINE set_3dparent_data
!
!***********************************************************************
      SUBROUTINE set_3dparent_data_tile (ng, tile,                      &
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
      USE interpolate_mod

#  ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange4d
      USE mp_exchange_mod, ONLY : mp_exchange2d, mp_exchange3d
#  endif
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

#  ifdef OBC
      integer :: ILB, IUB, JLB, JUB
#  endif
      integer :: i, j, k, ngp, ic, jc, t1, t2, itrc
      real(r8) :: my_val, fac, fac1, fac2

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
#  include "set_bounds.h"
      ngp=myparent(ng)
!

#  ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn on input data time wall clock.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, iNLM, 3)
#  endif
!set the parent grid to get data from.
!
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

!      if (myrank.eq.mymaster) then
!      write(*,*) 'set3dp ',t1,t2,fac1,fac2
!      end if
!
! tracer
!
      IF (WESTERN_EDGE) THEN
        DO itrc=1,NT(ng)
          DO k=1,N(ng)
            DO j=JstrT,JendT
              DO i=IstrT,Istr-1
                ic=i+Nghostpoints+1
                my_val=fac1*REFINED(ng)%tR_west(ic,j,k,t1,itrc)+        &
     &                 fac2*REFINED(ng)%tR_west(ic,j,k,t2,itrc)
#  ifdef MASKING
                my_val=my_val*GRID(ng)%rmask(i,j)
#  endif
                OCEAN(ng)%t(i,j,k,3,itrc)=my_val
                OCEAN(ng)%t(i,j,k,2,itrc)=my_val
                OCEAN(ng)%t(i,j,k,1,itrc)=my_val
               END DO
            END DO
          END DO
        END DO
      END IF
      IF (EASTERN_EDGE) THEN
        DO itrc=1,NT(ng)
          DO k=1,N(ng)
            DO j=JstrT,JendT
              DO i=Iend+1,IendT
                ic=i-Iend
                my_val=fac1*REFINED(ng)%tR_east(ic,j,k,t1,itrc)+        &
     &                 fac2*REFINED(ng)%tR_east(ic,j,k,t2,itrc)
#  ifdef MASKING
                my_val=my_val*GRID(ng)%rmask(i,j)
#  endif
                OCEAN(ng)%t(i,j,k,3,itrc)=my_val
                OCEAN(ng)%t(i,j,k,2,itrc)=my_val
                OCEAN(ng)%t(i,j,k,1,itrc)=my_val
              END DO
            END DO
          END DO
        END DO
      END IF
      IF (SOUTHERN_EDGE) THEN
        DO itrc=1,NT(ng)
          DO k=1,N(ng)
            DO j=JstrT,Jstr-1
              DO i=IstrT,IendT
                jc=j+Nghostpoints+1
                my_val=fac1*REFINED(ng)%tR_south(i,jc,k,t1,itrc)+       &
     &                 fac2*REFINED(ng)%tR_south(i,jc,k,t2,itrc)
#  ifdef MASKING
                my_val=my_val*GRID(ng)%rmask(i,j)
#  endif
                OCEAN(ng)%t(i,j,k,3,itrc)=my_val
                OCEAN(ng)%t(i,j,k,2,itrc)=my_val
                OCEAN(ng)%t(i,j,k,1,itrc)=my_val
              END DO
            END DO
          END DO
        END DO
      END IF
      IF (NORTHERN_EDGE) THEN
        DO itrc=1,NT(ng)
          DO k=1,N(ng)
            DO j=Jend+1,JendT
              DO i=IstrT,IendT
                jc=j-Jend
                my_val=fac1*REFINED(ng)%tR_north(i,jc,k,t1,itrc)+       &
     &                 fac2*REFINED(ng)%tR_north(i,jc,k,t2,itrc)
#  ifdef MASKING
                my_val=my_val*GRID(ng)%rmask(i,j)
#  endif
                OCEAN(ng)%t(i,j,k,3,itrc)=my_val
                OCEAN(ng)%t(i,j,k,2,itrc)=my_val
                OCEAN(ng)%t(i,j,k,1,itrc)=my_val
              END DO
            END DO
          END DO
        END DO
      END IF
!
! u
!
      IF (WESTERN_EDGE) THEN
        DO k=1,N(ng)
          DO j=JstrT,JendT
            DO i=IstrT+1,Istr
              ic=i+Nghostpoints
              my_val=fac1*REFINED(ng)%uR_west(ic,j,k,t1)+               &
     &               fac2*REFINED(ng)%uR_west(ic,j,k,t2)
#  ifdef MASKING
              my_val=my_val*GRID(ng)%umask(i,j)
#  endif
              OCEAN(ng)%u(i,j,k,2)=my_val
              OCEAN(ng)%u(i,j,k,1)=my_val
            END DO
          END DO
        END DO
      END IF
      IF (EASTERN_EDGE) THEN
        DO k=1,N(ng)
          DO j=JstrT,JendT
            DO i=Iend+1,IendT
              ic=i-Iend
              my_val=fac1*REFINED(ng)%uR_east(ic,j,k,t1)+               &
     &               fac2*REFINED(ng)%uR_east(ic,j,k,t2)
#  ifdef MASKING
              my_val=my_val*GRID(ng)%umask(i,j)
#  endif
              OCEAN(ng)%u(i,j,k,1)=my_val
              OCEAN(ng)%u(i,j,k,2)=my_val
            END DO
          END DO
        END DO
      END IF
      IF (SOUTHERN_EDGE) THEN
        DO k=1,N(ng)
          DO j=JstrT,Jstr-1
            DO i=IstrTU+1,IendT
              jc=j+Nghostpoints+1
              my_val=fac1*REFINED(ng)%uR_south(i,jc,k,t1)+              &
     &               fac2*REFINED(ng)%uR_south(i,jc,k,t2)
#  ifdef MASKING
              my_val=my_val*GRID(ng)%umask(i,j)
#  endif
              OCEAN(ng)%u(i,j,k,1)=my_val
              OCEAN(ng)%u(i,j,k,2)=my_val
            END DO
          END DO
        END DO
      END IF
      IF (NORTHERN_EDGE) THEN
        DO k=1,N(ng)
          DO j=Jend+1,JendT
            DO i=IstrTU+1,IendT
              jc=j-Jend
              my_val=fac1*REFINED(ng)%uR_north(i,jc,k,t1)+              &
     &               fac2*REFINED(ng)%uR_north(i,jc,k,t2)
#  ifdef MASKING
              my_val=my_val*GRID(ng)%umask(i,j)
#  endif
              OCEAN(ng)%u(i,j,k,1)=my_val
              OCEAN(ng)%u(i,j,k,2)=my_val
            END DO
          END DO
        END DO
      END IF
#  ifdef WEC
!
! u_stokes
!
      IF (WESTERN_EDGE) THEN
        DO k=1,N(ng)
          DO j=JstrT,JendT
            DO i=IstrT+1,Istr
              ic=i+Nghostpoints
              my_val=fac1*REFINED(ng)%usR_west(ic,j,k,t1)+              &
     &               fac2*REFINED(ng)%usR_west(ic,j,k,t2)
#   ifdef MASKING
              my_val=my_val*GRID(ng)%umask(i,j)
#   endif
              OCEAN(ng)%u_stokes(i,j,k)=my_val
            END DO
          END DO
        END DO
      END IF
      IF (EASTERN_EDGE) THEN
        DO k=1,N(ng)
          DO j=JstrT,JendT
            DO i=Iend+1,IendT
              ic=i-Iend
              my_val=fac1*REFINED(ng)%usR_east(ic,j,k,t1)+              &
     &               fac2*REFINED(ng)%usR_east(ic,j,k,t2)
#   ifdef MASKING
              my_val=my_val*GRID(ng)%umask(i,j)
#   endif
              OCEAN(ng)%u_stokes(i,j,k)=my_val
            END DO
          END DO
        END DO
      END IF
      IF (SOUTHERN_EDGE) THEN
        DO k=1,N(ng)
          DO j=JstrT,Jstr-1
            DO i=IstrTU+1,IendT
              jc=j+Nghostpoints+1
              my_val=fac1*REFINED(ng)%usR_south(i,jc,k,t1)+             &
     &               fac2*REFINED(ng)%usR_south(i,jc,k,t2)
#   ifdef MASKING
              my_val=my_val*GRID(ng)%umask(i,j)
#   endif
              OCEAN(ng)%u_stokes(i,j,k)=my_val
            END DO
          END DO
        END DO
      END IF
      IF (NORTHERN_EDGE) THEN
        DO k=1,N(ng)
          DO j=Jend+1,JendT
            DO i=IstrTU+1,IendT
              jc=j-Jend
              my_val=fac1*REFINED(ng)%usR_north(i,jc,k,t1)+             &
     &               fac2*REFINED(ng)%usR_north(i,jc,k,t2)
#   ifdef MASKING
              my_val=my_val*GRID(ng)%umask(i,j)
#   endif
              OCEAN(ng)%u_stokes(i,j,k)=my_val
            END DO
          END DO
        END DO
      END IF
#  endif
!
! v
!
      IF (WESTERN_EDGE) THEN
        DO k=1,N(ng)
          DO j=JstrTV+1,JendT
            DO i=IstrT,Istr-1
              ic=i+Nghostpoints+1
              my_val=fac1*REFINED(ng)%vR_west(ic,j,k,t1)+               &
     &               fac2*REFINED(ng)%vR_west(ic,j,k,t2)
#  ifdef MASKING
              my_val=my_val*GRID(ng)%vmask(i,j)
#  endif
              OCEAN(ng)%v(i,j,k,2)=my_val
              OCEAN(ng)%v(i,j,k,1)=my_val
            END DO
          END DO
        END DO
      END IF
      IF (EASTERN_EDGE) THEN
        DO k=1,N(ng)
          DO j=JstrTV+1,JendT
            DO i=Iend+1,IendT
              ic=i-Iend
              my_val=fac1*REFINED(ng)%vR_east(ic,j,k,t1)+               &
     &               fac2*REFINED(ng)%vR_east(ic,j,k,t2)
#  ifdef MASKING
              my_val=my_val*GRID(ng)%vmask(i,j)
#  endif
              OCEAN(ng)%v(i,j,k,1)=my_val
              OCEAN(ng)%v(i,j,k,2)=my_val
            END DO
          END DO
        END DO
      END IF
      IF (SOUTHERN_EDGE) THEN
        DO k=1,N(ng)
          DO j=JstrT+1,Jstr
            DO i=IstrT,IendT
              jc=j+Nghostpoints
              my_val=fac1*REFINED(ng)%vR_south(i,jc,k,t1)+              &
     &               fac2*REFINED(ng)%vR_south(i,jc,k,t2)
#  ifdef MASKING
              my_val=my_val*GRID(ng)%vmask(i,j)
#  endif
              OCEAN(ng)%v(i,j,k,1)=my_val
              OCEAN(ng)%v(i,j,k,2)=my_val
            END DO
          END DO
        END DO
      END IF
      IF (NORTHERN_EDGE) THEN
        DO k=1,N(ng)
          DO j=Jend+1,JendT
            DO i=IstrT,IendT
              jc=j-Jend
              my_val=fac1*REFINED(ng)%vR_north(i,jc,k,t1)+              &
     &               fac2*REFINED(ng)%vR_north(i,jc,k,t2)
#  ifdef MASKING
              my_val=my_val*GRID(ng)%vmask(i,j)
#  endif
              OCEAN(ng)%v(i,j,k,1)=my_val
              OCEAN(ng)%v(i,j,k,2)=my_val
            END DO
          END DO
        END DO
      END IF
#  ifdef WEC
!
! v_stokes
!
      IF (WESTERN_EDGE) THEN
        DO k=1,N(ng)
          DO j=JstrTV+1,JendT
            DO i=IstrT,Istr-1
              ic=i+Nghostpoints+1
              my_val=fac1*REFINED(ng)%vsR_west(ic,j,k,t1)+              &
     &               fac2*REFINED(ng)%vsR_west(ic,j,k,t2)
#   ifdef MASKING
              my_val=my_val*GRID(ng)%vmask(i,j)
#   endif
              OCEAN(ng)%v_stokes(i,j,k)=my_val
            END DO
          END DO
        END DO
      END IF
      IF (EASTERN_EDGE) THEN
        DO k=1,N(ng)
          DO j=JstrTV+1,JendT
            DO i=Iend+1,IendT
              ic=i-Iend
              my_val=fac1*REFINED(ng)%vsR_east(ic,j,k,t1)+              &
     &               fac2*REFINED(ng)%vsR_east(ic,j,k,t2)
#   ifdef MASKING
              my_val=my_val*GRID(ng)%vmask(i,j)
#   endif
              OCEAN(ng)%v_stokes(i,j,k)=my_val
            END DO
          END DO
        END DO
      END IF
      IF (SOUTHERN_EDGE) THEN
        DO k=1,N(ng)
          DO j=JstrT+1,Jstr
            DO i=IstrT,IendT
              jc=j+Nghostpoints
              my_val=fac1*REFINED(ng)%vsR_south(i,jc,k,t1)+             &
     &               fac2*REFINED(ng)%vsR_south(i,jc,k,t2)
#   ifdef MASKING
              my_val=my_val*GRID(ng)%vmask(i,j)
#   endif
              OCEAN(ng)%v_stokes(i,j,k)=my_val
            END DO
          END DO
        END DO
      END IF
      IF (NORTHERN_EDGE) THEN
        DO k=1,N(ng)
          DO j=Jend+1,JendT
            DO i=IstrT,IendT
              jc=j-Jend
              my_val=fac1*REFINED(ng)%vsR_north(i,jc,k,t1)+             &
     &               fac2*REFINED(ng)%vsR_north(i,jc,k,t2)
#   ifdef MASKING
              my_val=my_val*GRID(ng)%vmask(i,j)
#   endif
              OCEAN(ng)%v_stokes(i,j,k)=my_val
            END DO
          END DO
        END DO
      END IF
#  endif
!
#  ifdef DISTRIBUTE
#   if !defined TS_FIXED
      CALL mp_exchange4d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 1, N(ng), 1, NT(ng),      &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    OCEAN(ng)%t(:,:,:,1,:))
      CALL mp_exchange4d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 1, N(ng), 1, NT(ng),      &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    OCEAN(ng)%t(:,:,:,2,:))
      CALL mp_exchange4d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 1, N(ng), 1, NT(ng),      &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    OCEAN(ng)%t(:,:,:,3,:))
#   endif
      CALL mp_exchange3d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    OCEAN(ng)%u(:,:,:,1),                         &
     &                    OCEAN(ng)%u(:,:,:,2))
      CALL mp_exchange3d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    OCEAN(ng)%v(:,:,:,1),                         &
     &                    OCEAN(ng)%v(:,:,:,2))
#   ifdef WEC
      CALL mp_exchange3d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    OCEAN(ng)%u_stokes(:,:,:),                    &
     &                    OCEAN(ng)%v_stokes(:,:,:))
#   endif
#  endif
#  ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn off input data time wall clock.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, iNLM, 3)
#  endif
      RETURN
      END SUBROUTINE set_3dparent_data_tile
# endif
#endif
      END MODULE set_3dparent_data_mod
