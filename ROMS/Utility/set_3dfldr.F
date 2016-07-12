#include "cppdefs.h"
      MODULE set_3dfldr_mod
#if defined ADJOINT && defined SOLVE3D
!
!svn $Id: set_3dfldr.F 795 2016-05-11 01:42:43Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine time-interpolates backwards in time requested 3D field !
!  from snapshots of input data.                                       !
!                                                                      !
!=======================================================================
!
      implicit none

      CONTAINS
!
!***********************************************************************
      SUBROUTINE set_3dfldr_tile (ng, tile, model, ifield,              &
     &                            LBi, UBi, LBj, UBj, LBk, UBk,         &
     &                            Finp, Fout, update,                   &
     &                            SetBC)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_ncparam
      USE mod_scalars
!
      USE exchange_3d_mod
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange3d
# endif
!
!  Imported variable declarations.
!
      logical, intent(in), optional :: SetBC

      logical, intent(out) :: update

      integer, intent(in) :: ng, tile, model, ifield
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBk, UBk
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: Finp(LBi:,LBj:,LBk:,:)
      real(r8), intent(out) :: Fout(LBi:,LBj:,LBk:)
# else
      real(r8), intent(in) :: Finp(LBi:UBi,LBj:UBj,LBk:UBk,2)
      real(r8), intent(out) :: Fout(LBi:UBi,LBj:UBj,LBk:UBk)
# endif
!
!  Local variable declarations.
!
      logical :: LapplyBC, Lgrided, Lonerec

      integer :: Tindex, gtype, i, it1, it2, j, k

      real(r8) :: Fval, fac, fac1, fac2

# include "set_bounds.h"
!
!--------------------------------------------------------------------
!  Set-up requested field for current tile.
!--------------------------------------------------------------------
!
!  Set switch to apply boundary conditions.
!
      IF (PRESENT(SetBC)) THEN
        LapplyBC=SetBC
      ELSE
        LapplyBC=.TRUE.
      END IF
!
!  Get requested field information from global storage.
!
      Lgrided=Linfo(1,ifield,ng)
      Lonerec=Linfo(3,ifield,ng)
      gtype  =Iinfo(1,ifield,ng)
      Tindex =Iinfo(8,ifield,ng)
      update=.TRUE.
!
!  Set linear-interpolation factors.
!
      it1=3-Tindex
      it2=Tindex
# if defined CRAY || defined SGI
      fac1=ANINT(time(ng)-Tintrp(it2,ifield,ng))
      fac2=ANINT(Tintrp(it1,ifield,ng)-time(ng))
# else
      fac1=ANINT(time(ng)-Tintrp(it2,ifield,ng),r8)
      fac2=ANINT(Tintrp(it1,ifield,ng)-time(ng),r8)
# endif
!
!  Load time-invariant data. Time interpolation is not necessary.
!
      IF (Lonerec) THEN
        IF (Lgrided) THEN
          DO k=LBk,UBk
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                Fout(i,j,k)=Finp(i,j,k,Tindex)
              END DO
            END DO
          END DO
        ELSE
          Fval=Fpoint(Tindex,ifield,ng)
          DO k=LBk,UBk
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                Fout(i,j,k)=Fval
              END DO
            END DO
          END DO
        END IF
!
!  Time-interpolate from gridded or point data.
!
      ELSE IF (((fac1*fac2).ge.0.0_r8).and.                             &
     &        ((fac1+fac2).gt.0.0_r8)) THEN
        fac=1.0_r8/(fac1+fac2)
        fac1=fac*fac1
        fac2=fac*fac2
        IF (Lgrided) THEN
          DO k=LBk,UBk
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                Fout(i,j,k)=fac1*Finp(i,j,k,it1)+fac2*Finp(i,j,k,it2)
              END DO
            END DO
          END DO
        ELSE
          Fval=fac1*Fpoint(it1,ifield,ng)+fac2*Fpoint(it2,ifield,ng)
          DO k=LBk,UBk
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                Fout(i,j,k)=Fval
              END DO
            END DO
          END DO
        END IF
!
!  Activate synchronization flag if a new time record needs to be
!  read in at the next time step.
!
        IF ((time(ng)-dt(ng)).lt.Tintrp(it2,ifield,ng)) THEN
          IF (DOMAIN(ng)%SouthWest_Test(tile)) synchro_flag(ng)=.TRUE.
        END IF
!
!  Unable to set-up requested field.  Activate error flag to quit.
!
      ELSE
        IF (DOMAIN(ng)%SouthWest_Test(tile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,ifield)), tdays(ng),         &
     &                        Finfo(1,ifield,ng), Finfo(2,ifield,ng),   &
     &                        Finfo(3,ifield,ng), Finfo(4,ifield,ng),   &
     &                        Tintrp(it1,ifield,ng)*sec2day,            &
     &                        Tintrp(it2,ifield,ng)*sec2day,            &
     &                        fac1*sec2day, fac2*sec2day
          END IF
  10      FORMAT (/,' SET_3DFLDR - current model time',                 &
     &            ' exceeds ending value for variable: ',a,             &
     &            /,14x,'TDAYS     = ',f15.4,                           &
     &            /,14x,'Data Tmin = ',f15.4,2x,'Data Tmax = ',f15.4,   &
     &            /,14x,'Data Tstr = ',f15.4,2x,'Data Tend = ',f15.4,   &
     &            /,14x,'TINTRP1   = ',f15.4,2x,'TINTRP2   = ',f15.4,   &
     &            /,14x,'FAC1      = ',f15.4,2x,'FAC2      = ',f15.4)
          exit_flag=2
          update=.FALSE.
        END IF
      END IF
!
!  Exchange boundary data.
!
      IF (update) THEN
        IF (LapplyBC.and.(EWperiodic(ng).or.NSperiodic(ng))) THEN
          IF (gtype.eq.r3dvar) THEN
            CALL exchange_r3d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, LBk, UBk,       &
     &                              Fout)
          ELSE IF (gtype.eq.u3dvar) THEN
            CALL exchange_u3d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, LBk, UBk,       &
     &                              Fout)
          ELSE IF (gtype.eq.v3dvar) THEN
            CALL exchange_v3d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, LBk, UBk,       &
     &                              Fout)
          ELSE IF (gtype.eq.w3dvar) THEN
            CALL exchange_w3d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, LBk, UBk,       &
     &                              Fout)
          END IF
        END IF

# ifdef DISTRIBUTE
        IF (.not.LapplyBC) THEN
          CALL mp_exchange3d (ng, tile, model, 1,                       &
     &                        LBi, UBi, LBj, UBj, LBk, UBk,             &
     &                        NghostPoints,                             &
     &                        .FALSE., .FALSE.,                         &
     &                        Fout)
        ELSE
          CALL mp_exchange3d (ng, tile, model, 1,                       &
     &                        LBi, UBi, LBj, UBj, LBk, UBk,             &
     &                        NghostPoints,                             &
     &                        EWperiodic(ng), NSperiodic(ng),           &
     &                        Fout)
        END IF
# endif
      END IF

      RETURN
      END SUBROUTINE set_3dfldr_tile
#endif
      END MODULE set_3dfldr_mod
