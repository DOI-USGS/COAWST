      SUBROUTINE ana_nudgcoef (ng, tile, model)
!
!! svn $Id: ana_nudgcoef.h 34 2007-04-27 04:40:21Z arango $
!!================================================= Hernan G. Arango ===
!! Copyright (c) 2002-2007 The ROMS/TOMS Group                         !
!!   Licensed under a MIT/X style license                              !
!!   See License_ROMS.txt                                              !
!!======================================================================
!                                                                      !
!  This routine set nudging coefficients time-scales (1/s).            !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_ncparam
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
!
!  Local variable declarations.
!
#include "tile.h"
!
      CALL ana_nudgcoef_tile (ng, model, Istr, Iend, Jstr, Jend)
!
! Set analytical header file name used.
!
      IF (Lanafile) THEN
        WRITE (ANANAME(16),'(a,a)') TRIM(Adir), '/ana_nudgcoef.h'
      END IF

      RETURN
      END SUBROUTINE ana_nudgcoef
!
!***********************************************************************
      SUBROUTINE ana_nudgcoef_tile (ng, model, Istr, Iend, Jstr, Jend)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_boundary
#ifdef CLIMATOLOGY
      USE mod_clima
#endif
      USE mod_grid
      USE mod_ncparam
      USE mod_scalars
#ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_bcastf
#endif
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, Iend, Istr, Jend, Jstr
!
!  Local variable declarations.
!
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: Iwrk, i, itrc, j

      real(r8) :: cff1, cff2

      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY) :: wrk

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Set up nudging towards data time-scale coefficients (1/s).
!-----------------------------------------------------------------------
!
!  Initialize.
!
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          wrk(i,j)=0.0_r8
        END DO
      END DO

#if defined NATL
!
!  Set tracer nudging coefficients in the southern and northern edges
!  from a five days time scale at the boundary point to decrease
!  linearly to 60 days seven grids points away from the boundary.
!
      cff1=1.0_r8/(5.0_r8*86400.0_r8)
      cff2=1.0_r8/(60.0_r8*86400.0_r8)
      DO j=JstrR,MIN(8,JendR)
        DO i=IstrR,IendR
          wrk(i,j)=cff2+REAL(8-j,r8)*(cff1-cff2)/7.0_r8
        END DO
      END DO
      DO j=MAX(JstrR,Mm(ng)-7),JendR
        DO i=IstrR,IendR
          wrk(i,j)=cff1+REAL(Mm(ng)-j,r8)*(cff2-cff1)/7.0_r8
        END DO
      END DO
# ifdef ZCLM_NUDGING
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          CLIMA(ng)%Znudgcof(i,j)=wrk(i,j)
        END DO
      END DO
# endif
# ifdef TCLM_NUDGING
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          CLIMA(ng)%Tnudgcof(i,j,itemp)=wrk(i,j)
          CLIMA(ng)%Tnudgcof(i,j,isalt)=wrk(i,j)
        END DO
      END DO
# endif
#elif defined SCB
!
!  Set tracer nudging coefficients in the southern and northern edges
!  from a five days time scale at the boundary point to decrease
!  linearly to 60 days six grids points away from the boundary.
!
      cff1=1.0_r8/(5.0_r8*86400.0_r8)
      cff2=1.0_r8/(60.0_r8*86400.0_r8)
      DO j=JstrR,MIN(6,JendR)
        DO i=IstrR,IendR
          wrk(i,j)=cff2+REAL(6-j,r8)*(cff1-cff2)/6.0_r8
        END DO
      END DO
      DO j=MAX(JstrR,Mm(ng)+1-6),JendR
        DO i=IstrR,IendR
          wrk(i,j)=cff1+REAL(Mm(ng)+1-j,r8)*(cff2-cff1)/6.0_r8
        END DO
      END DO
      DO i=IstrR,MIN(6,IendR)
        DO j=MAX(JstrR,i),MIN(Mm(ng)+1-i,JendR)
          wrk(i,j)=cff2+REAL(6-i,r8)*(cff1-cff2)/6.0_r8
        END DO
      END DO
# ifdef ZCLM_NUDGING
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          CLIMA(ng)%Znudgcof(i,j)=wrk(i,j)
        END DO
      END DO
# endif
# ifdef M2CLM_NUDGING
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          CLIMA(ng)%M2nudgcof(i,j)=wrk(i,j)
        END DO
      END DO
# endif
# ifdef M3CLM_NUDGING
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          CLIMA(ng)%M3nudgcof(i,j)=wrk(i,j)
        END DO
      END DO
# endif
# ifdef TCLM_NUDGING
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          CLIMA(ng)%Tnudgcof(i,j,itemp)=wrk(i,j)
          CLIMA(ng)%Tnudgcof(i,j,isalt)=wrk(i,j)
        END DO
      END DO
# endif
!!
!!  Specify your spatially varying nudging coefficients here.
!!
!!#elif defined ???

#else
!
!  Default nudging coefficients.  Set nudging coefficients uniformly to
!  the values specified in the standard input file.
!
# ifdef ZCLM_NUDGING
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          CLIMA(ng)%Znudgcof(i,j)=Znudg(ng)
        END DO
      END DO
# endif
# ifdef M2CLM_NUDGING
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          CLIMA(ng)%M2nudgcof(i,j)=M2nudg(ng)
        END DO
      END DO
# endif
# ifdef SOLVE3D
#  ifdef TCLM_NUDGING
      DO itrc=1,NT(ng)
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            CLIMA(ng)%Tnudgcof(i,j,itrc)=Tnudg(itrc,ng)
          END DO
        END DO
      END DO
#  endif
#  ifdef M3CLM_NUDGING
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          CLIMA(ng)%M3nudgcof(i,j)=M3nudg(ng)
        END DO
      END DO
#  endif
# endif
#endif
!
!-----------------------------------------------------------------------
!  Set nudging coefficients (1/s) for passive/active (outflow/inflow)
!  open boundary conditions.  Weak nudging is expected in passive
!  outflow conditions and strong nudging is expected in active inflow
!  conditions.  Notice that interior nudging coefficient defined
!  above are zero out when boundary condition nudging.  The USER needs
!  to adapt this to his/her application!
!-----------------------------------------------------------------------
!
!! WARNING:  This section is generic for all applications. Please do not
!!           change the code below.
!!            
!  Free-surface nudging coefficients.
!
#ifdef WEST_FSNUDGING
# ifdef ZCLM_NUDGING
      IF (SOUTH_WEST_CORNER) THEN
        FSobc_out(ng,iwest)=CLIMA(ng)%Znudgcof(0,1)
        FSobc_in (ng,iwest)=obcfac(ng)*FSobc_out(ng,iwest)
      END IF
      IF (WESTERN_EDGE) THEN
        DO j=JstrR,JendR
          CLIMA(ng)%Znudgcof(0,j)=0.0_r8
        END DO
      END IF
#  ifdef DISTRIBUTE
      IF (ng.eq.Ngrids) THEN
        CALL mp_bcastf (ng, model, FSobc_out(:,iwest), Ngrids)
        CALL mp_bcastf (ng, model, FSobc_in (:,iwest), Ngrids)
      ENDIF
#  endif
# else
      IF (SOUTH_WEST_TEST) THEN
        FSobc_out(ng,iwest)=Znudg(ng)
        FSobc_in (ng,iwest)=obcfac(ng)*Znudg(ng)
      END IF
# endif
#endif
#ifdef EAST_FSNUDGING
# ifdef ZCLM_NUDGING
      IF (NORTH_EAST_CORNER) THEN
        FSobc_out(ng,ieast)=CLIMA(ng)%Znudgcof(Lm(ng)+1,Mm(ng))
        FSobc_in (ng,ieast)=obcfac(ng)*FSobc_out(ng,ieast)
      END IF
      IF (EASTERN_EDGE) THEN
        DO j=JstrR,JendR
          CLIMA(ng)%Znudgcof(Lm(ng)+1,j)=0.0_r8
        END DO
      END IF
#  ifdef DISTRIBUTE
      IF (ng.eq.Ngrids) THEN
        CALL mp_bcastf (ng, model, FSobc_out(:,ieast), Ngrids)
        CALL mp_bcastf (ng, model, FSobc_in (:,ieast), Ngrids)
      ENDIF
#  endif
# else
      IF (NORTH_EAST_TEST) THEN
        FSobc_out(ng,ieast)=Znudg(ng)
        FSobc_in (ng,ieast)=obcfac(ng)*Znudg(ng)
      END IF
# endif
#endif
#ifdef SOUTH_FSNUDGING
# ifdef ZCLM_NUDGING
      IF (SOUTH_WEST_CORNER) THEN
        FSobc_out(ng,isouth)=CLIMA(ng)%Znudgcof(1,0)
        FSobc_in (ng,isouth)=obcfac(ng)*FSobc_out(ng,isouth)
      END IF
      IF (SOUTHERN_EDGE) THEN
        DO i=IstrR,IendR
          CLIMA(ng)%Znudgcof(i,0)=0.0_r8
        END DO
      END IF
#  ifdef DISTRIBUTE
      IF (ng.eq.Ngrids) THEN
        CALL mp_bcastf (ng, model, FSobc_out(:,iwest), Ngrids)
        CALL mp_bcastf (ng, model, FSobc_in (:,iwest), Ngrids)
      ENDIF
#  endif
# else
      IF (SOUTH_WEST_TEST) THEN
        FSobc_out(ng,isouth)=Znudg(ng)
        FSobc_in (ng,isouth)=obcfac(ng)*Znudg(ng)
      END IF
# endif
#endif
#ifdef NORTH_FSNUDGING
# ifdef ZCLM_NUDGING
      IF (NORTH_EAST_CORNER) THEN
        FSobc_out(ng,inorth)=CLIMA(ng)%Znudgcof(Lm(ng),Mm(ng)+1)
        FSobc_in (ng,inorth)=obcfac(ng)*FSobc_out(ng,inorth)
      END IF
      IF (NORTHERN_EDGE) THEN
        DO i=IstrR,IendR
          CLIMA(ng)%Znudgcof(i,Mm(ng)+1)=0.0_r8
        END DO
      END IF
#  ifdef DISTRIBUTE
      IF (ng.eq.Ngrids) THEN
        CALL mp_bcastf (ng, model, FSobc_out(:,inorth), Ngrids)
        CALL mp_bcastf (ng, model, FSobc_in (:,inorth), Ngrids)
      ENDIF
#  endif
# else
      IF (NORTH_EAST_TEST) THEN
        FSobc_out(ng,inorth)=Znudg(ng)
        FSobc_in (ng,inorth)=obcfac(ng)*Znudg(ng)
      END IF
# endif
#endif
!
!  2D momentum nudging coefficients.
!
#ifdef WEST_M2NUDGING
# ifdef M2CLM_NUDGING
      IF (SOUTH_WEST_CORNER) THEN
        M2obc_out(ng,iwest)=0.5_r8*(CLIMA(ng)%M2nudgcof(0,1)+           &
     &                              CLIMA(ng)%M2nudgcof(1,1))
        M2obc_in (ng,iwest)=obcfac(ng)*M2obc_out(ng,iwest)
      END IF
      IF (WESTERN_EDGE) THEN
        DO j=JstrR,JendR
          CLIMA(ng)%M2nudgcof(0,j)=-CLIMA(ng)%M2nudgcof(1,j)
        END DO
      END IF
#  ifdef DISTRIBUTE
      IF (ng.eq.Ngrids) THEN
        CALL mp_bcastf (ng, model, M2obc_out(:,iwest), Ngrids)
        CALL mp_bcastf (ng, model, M2obc_in (:,iwest), Ngrids)
      ENDIF
#  endif
# else
      IF (SOUTH_WEST_TEST) THEN
        M2obc_out(ng,iwest)=M2nudg(ng)
        M2obc_in (ng,iwest)=obcfac(ng)*M2nudg(ng)
      END IF
# endif
#endif
#ifdef EAST_M2NUDGING
# ifdef M2CLM_NUDGING
      IF (NORTH_EAST_CORNER) THEN
        M2obc_out(ng,ieast)=0.5_r8*                                     &
     &                      (CLIMA(ng)%M2nudgcof(Lm(ng)  ,Mm(ng))+      &
     &                       CLIMA(ng)%M2nudgcof(Lm(ng)+1,Mm(ng)))
        M2obc_in (ng,ieast)=obcfac(ng)*M2obc_out(ng,ieast)
      END IF
      IF (EASTERN_EDGE) THEN
        DO j=JstrR,JendR
          CLIMA(ng)%M2nudgcof(Lm(ng)+1,j)=-CLIMA(ng)%M2nudgcof(Lm(ng),j)
        END DO
      END IF
#  ifdef DISTRIBUTE
      IF (ng.eq.Ngrids) THEN
        CALL mp_bcastf (ng, model, M2obc_out(:,ieast), Ngrids)
        CALL mp_bcastf (ng, model, M2obc_in (:,ieast), Ngrids)
      ENDIF
#  endif
# else
      IF (NORTH_EAST_TEST) THEN
        M2obc_out(ng,ieast)=M2nudg(ng)
        M2obc_in (ng,ieast)=obcfac(ng)*M2nudg(ng)
      END IF
# endif
#endif
#ifdef SOUTH_M2NUDGING
# ifdef M2CLM_NUDGING
      IF (SOUTH_WEST_CORNER) THEN
        M2obc_out(ng,isouth)=0.5_r8*(CLIMA(ng)%M2nudgcof(1,0)+          &
     &                               CLIMA(ng)%M2nudgcof(1,1))
        M2obc_in (ng,isouth)=obcfac(ng)*M2obc_out(ng,isouth)
      END IF
      IF (SOUTHERN_EDGE) THEN
        DO i=IstrR,IendR
          CLIMA(ng)%M2nudgcof(i,0)=-CLIMA(ng)%M2nudgcof(i,1)
        END DO
      END IF
#  ifdef DISTRIBUTE
      IF (ng.eq.Ngrids) THEN
        CALL mp_bcastf (ng, model, M2obc_out(:,isouth), Ngrids)
        CALL mp_bcastf (ng, model, M2obc_in (:,isouth), Ngrids)
      ENDIF
#  endif
# else
      IF (SOUTH_WEST_TEST) THEN
        M2obc_out(ng,isouth)=M2nudg(ng)
        M2obc_in (ng,isouth)=obcfac(ng)*M2nudg(ng)
      END IF
# endif
#endif
#ifdef NORTH_M2NUDGING
# ifdef M2CLM_NUDGING
      IF (NORTH_EAST_CORNER) THEN
        M2obc_out(ng,inorth)=0.5_r8*                                    &
     &                       (CLIMA(ng)%M2nudgcof(Lm(ng),Mm(ng)  )+     &
     &                        CLIMA(ng)%M2nudgcof(Lm(ng),Mm(ng)+1))
        M2obc_in (ng,inorth)=obcfac(ng)*M2obc_out(ng,inorth)
      END IF
      IF (NORTHERN_EDGE) THEN
        DO i=IstrR,IendR
          CLIMA(ng)%M2nudgcof(i,Mm(ng)+1)=-CLIMA(ng)%M2nudgcof(i,Mm(ng))
        END DO
      END IF
#  ifdef DISTRIBUTE
      IF (ng.eq.Ngrids) THEN
        CALL mp_bcastf (ng, model, M2obc_out(:,inorth), Ngrids)
        CALL mp_bcastf (ng, model, M2obc_in (:,inorth), Ngrids)
      ENDIF
#  endif
# else
      IF (NORTH_EAST_TEST) THEN
        M2obc_out(ng,inorth)=M2nudg(ng)
        M2obc_in (ng,inorth)=obcfac(ng)*M2nudg(ng)
      END IF
# endif
#endif
#ifdef SOLVE3D
!
!  Tracers nudging coefficients.
!
# ifdef WEST_TNUDGING
#  ifdef TCLM_NUDGING
      DO itrc=1,NT(ng)
        IF (SOUTH_WEST_CORNER) THEN
          Tobc_out(itrc,ng,iwest)=CLIMA(ng)%Tnudgcof(0,1,itrc)
          Tobc_in (itrc,ng,iwest)=obcfac(ng)*Tobc_out(itrc,ng,iwest)
        END IF
        IF (WESTERN_EDGE) THEN
          DO j=JstrR,JendR
            CLIMA(ng)%Tnudgcof(0,j,itrc)=0.0_r8
          END DO
        END IF
      END DO
#   ifdef DISTRIBUTE
      CALL mp_bcastf (ng, model, Tobc_out(:,ng,iwest), MT)
      CALL mp_bcastf (ng, model, Tobc_in (:,ng,iwest), MT)
#   endif
#  else
      DO itrc=1,NT(ng)
        IF (SOUTH_WEST_TEST) THEN
          Tobc_out(itrc,ng,iwest)=Tnudg(itrc,ng)
          Tobc_in (itrc,ng,iwest)=obcfac(ng)*Tnudg(itrc,ng)
        END IF
      END DO
#  endif
# endif
# ifdef EAST_TNUDGING
#  ifdef TCLM_NUDGING
      DO itrc=1,NT(ng)
        IF (NORTH_EAST_CORNER) THEN
          Tobc_out(itrc,ng,ieast)=                                      &
     &             CLIMA(ng)%Tnudgcof(Lm(ng)+1,Mm(ng),itrc)
          Tobc_in (itrc,ng,ieast)=obcfac(ng)*Tobc_out(itrc,ng,ieast)
        END IF
        IF (EASTERN_EDGE) THEN
          DO j=JstrR,JendR
            CLIMA(ng)%Tnudgcof(Lm(ng)+1,j,itrc)=0.0_r8
          END DO
        END IF
      END DO
#   ifdef DISTRIBUTE
      CALL mp_bcastf (ng, model, Tobc_out(:,ng,ieast), MT)
      CALL mp_bcastf (ng, model, Tobc_in (:,ng,ieast), MT)
#   endif
#  else
      DO itrc=1,NT(ng)
        IF (NORTH_EAST_TEST) THEN
          Tobc_out(itrc,ng,ieast)=Tnudg(itrc,ng)
          Tobc_in (itrc,ng,ieast)=obcfac(ng)*Tnudg(itrc,ng)
        END IF
      END DO
#  endif
# endif
# ifdef SOUTH_TNUDGING
#  ifdef TCLM_NUDGING
      DO itrc=1,NT(ng)
        IF (SOUTH_WEST_CORNER) THEN
          Tobc_out(itrc,ng,isouth)=CLIMA(ng)%Tnudgcof(1,0,itrc)
          Tobc_in (itrc,ng,isouth)=obcfac(ng)*Tobc_out(itrc,ng,isouth)
        END IF
        IF (SOUTHERN_EDGE) THEN
          DO i=IstrR,IendR
            CLIMA(ng)%Tnudgcof(i,0,itrc)=0.0_r8
          END DO
        END IF
      END DO
#   ifdef DISTRIBUTE
      CALL mp_bcastf (ng, model, Tobc_out(:,ng,isouth), MT)
      CALL mp_bcastf (ng, model, Tobc_in (:,ng,isouth), MT)
#   endif
#  else
      DO itrc=1,NT(ng)
        IF (SOUTH_WEST_TEST) THEN
          Tobc_out(itrc,ng,isouth)=Tnudg(itrc,ng)
          Tobc_in (itrc,ng,isouth)=obcfac(ng)*Tnudg(itrc,ng)
        END IF
      END DO
#  endif
# endif
# ifdef NORTH_TNUDGING
#  ifdef TCLM_NUDGING
      DO itrc=1,NT(ng)
        IF (NORTH_EAST_CORNER) THEN
          Tobc_out(itrc,ng,inorth)=                                     &
     &             CLIMA(ng)%Tnudgcof(Lm(ng),Mm(ng)+1,itrc)
          Tobc_in (itrc,ng,inorth)=obcfac(ng)*Tobc_out(itrc,ng,inorth)
        END IF
        IF (NORTHERN_EDGE) THEN
          DO i=IstrR,IendR
            CLIMA(ng)%Tnudgcof(i,Mm(ng)+1,itrc)=0.0_r8
          END DO
        END IF
      END DO
#   ifdef DISTRIBUTE
      CALL mp_bcastf (ng, model, Tobc_out(:,ng,inorth), MT)
      CALL mp_bcastf (ng, model, Tobc_in (:,ng,inorth), MT)
#   endif
#  else
      DO itrc=1,NT(ng)
        IF (NORTH_EAST_TEST) THEN
          Tobc_out(itrc,ng,inorth)=Tnudg(itrc,ng)
          Tobc_in (itrc,ng,inorth)=obcfac(ng)*Tnudg(itrc,ng)
        END IF
      END DO
#  endif
# endif
!
!  3D momentum nudging coefficients.
!
# ifdef WEST_M3NUDGING
#  ifdef M3CLM_NUDGING
      IF (SOUTH_WEST_CORNER) THEN
        M3obc_out(ng,iwest)=0.5_r8*                                     &
     &                      (CLIMA(ng)%M3nudgcof(0,1)+                  &
     &                       CLIMA(ng)%M3nudgcof(1,1))
        M3obc_in (ng,iwest)=obcfac(ng)*M3obc_out(ng,iwest)
      END IF
      IF (WESTERN_EDGE) THEN
        DO j=JstrR,JendR
          CLIMA(ng)%M3nudgcof(0,j)=-CLIMA(ng)%M3nudgcof(1,j)
        END DO
      END IF
#   ifdef DISTRIBUTE
      IF (ng.eq.Ngrids) THEN
        CALL mp_bcastf (ng, model, M3obc_out(:,iwest), Ngrids)
        CALL mp_bcastf (ng, model, M3obc_in (:,iwest), Ngrids)
      ENDIF
#   endif
#  else
      IF (SOUTH_WEST_TEST) THEN
        M3obc_out(ng,iwest)=M3nudg(ng)
        M3obc_in (ng,iwest)=obcfac(ng)*M3nudg(ng)
      END IF
#  endif
# endif
# ifdef EAST_M3NUDGING
#  ifdef M3CLM_NUDGING
      IF (NORTH_EAST_CORNER) THEN
        M3obc_out(ng,ieast)=0.5_r8*                                     &
     &                      (CLIMA(ng)%M3nudgcof(Lm(ng)  ,Mm(ng))+      &
     &                       CLIMA(ng)%M3nudgcof(Lm(ng)+1,Mm(ng)))
        M3obc_in (ng,ieast)=obcfac(ng)*M3obc_out(ng,ieast)
      END IF
      IF (EASTERN_EDGE) THEN
        DO j=JstrR,JendR
          CLIMA(ng)%M3nudgcof(Lm(ng)+1,j)=-CLIMA(ng)%M3nudgcof(Lm(ng),j)
        END DO
      END IF
#   ifdef DISTRIBUTE
      IF (ng.eq.Ngrids) THEN
        CALL mp_bcastf (ng, model, M3obc_out(:,ieast), Ngrids)
        CALL mp_bcastf (ng, model, M3obc_in (:,ieast), Ngrids)
      ENDIF
#   endif
#  else
      IF (NORTH_EAST_TEST) THEN
        M3obc_out(ng,ieast)=M3nudg(ng)
        M3obc_in (ng,ieast)=obcfac(ng)*M3nudg(ng)
      END IF
#  endif
# endif
# ifdef SOUTH_M3NUDGING
#  ifdef M3CLM_NUDGING
      IF (SOUTH_WEST_CORNER) THEN
        M3obc_out(ng,isouth)=0.5_r8*                                    &
     &                       (CLIMA(ng)%M3nudgcof(1,0)+                 &
     &                        CLIMA(ng)%M3nudgcof(1,1))
        M3obc_in (ng,isouth)=obcfac(ng)*M3obc_out(ng,isouth)
      END IF
      IF (SOUTHERN_EDGE) THEN
        DO i=IstrR,IendR
          CLIMA(ng)%M3nudgcof(i,0)=-CLIMA(ng)%M3nudgcof(i,1)
        END DO
      END IF
#   ifdef DISTRIBUTE
      IF (ng.eq.Ngrids) THEN
        CALL mp_bcastf (ng, model, M3obc_out(:,isouth), Ngrids)
        CALL mp_bcastf (ng, model, M3obc_in (:,isouth), Ngrids)
      ENDIF
#   endif
#  else
      IF (SOUTH_WEST_TEST) THEN
        M3obc_out(ng,isouth)=M3nudg(ng)
        M3obc_in (ng,isouth)=obcfac(ng)*M3nudg(ng)
      END IF
#  endif
# endif
# ifdef NORTH_M3NUDGING
#  ifdef M3CLM_NUDGING
      IF (NORTH_EAST_CORNER) THEN
        M3obc_out(ng,inorth)=0.5_r8*                                    &
     &                       (CLIMA(ng)%M3nudgcof(Lm(ng),Mm(ng)  )+     &
     &                        CLIMA(ng)%M3nudgcof(Lm(ng),Mm(ng)+1))
        M3obc_in (ng,inorth)=obcfac(ng)*M3obc_out(ng,inorth)
      END IF
      IF (NORTHERN_EDGE) THEN
        DO i=IstrR,IendR
          CLIMA(ng)%M3nudgcof(i,Mm(ng)+1)=-CLIMA(ng)%M3nudgcof(i,Mm(ng))
        END DO
      END IF
#   ifdef DISTRIBUTE
      IF (ng.eq.Ngrids) THEN
        CALL mp_bcastf (ng, model, M3obc_out(:,inorth), Ngrids)
        CALL mp_bcastf (ng, model, M3obc_in (:,inorth), Ngrids)
      ENDIF
#   endif
#  else
      IF (NORTH_EAST_TEST) THEN
        M3obc_out(ng,inorth)=M3nudg(ng)
        M3obc_in (ng,inorth)=obcfac(ng)*M3nudg(ng)
      END IF
#  endif
# endif
#endif

      RETURN
      END SUBROUTINE ana_nudgcoef_tile

