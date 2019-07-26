#include "cppdefs.h"
#undef VARY_ACBC
#if defined INWAVE_MODEL & defined INWAVE_SWAN_COUPLING
      SUBROUTINE set_inwave_swan_data (ng, tile)
!
!svn $Id: set_inwave_swan_data.F 799 2009-12-08 20:38:55Z jcwarner $
!=======================================================================
!                                                                      !
!  This routine computes the action density for the boundary condition !
!  The action density is derived from the directional                  !
!  wave spectra given by SWAN                                          !
!=======================================================================
!
      USE mod_param
      USE mod_inwave_bound
!
!  Imported variable declarations.

      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
# include "tile.h"
!
# ifdef PROFILE
!      CALL wclock_on (ng, iNLM, 4)
# endif
      CALL set_inwave_swan_data_tile (ng, tile,                         &
     &                                LBi, UBi, LBj, UBj)
# ifdef PROFILE
!      CALL wclock_off (ng, iNLM, 4)
# endif
      RETURN
      END SUBROUTINE set_inwave_swan_data
!
!***********************************************************************
      SUBROUTINE set_inwave_swan_data_tile (ng, tile,                   &
     &                                      LBi, UBi, LBj, UBj)
!***********************************************************************
!
      USE mod_param
      USE mod_ncparam
      USE mod_scalars
      USE mod_inwave_params
      USE mod_inwave_vars
      USE mod_inwave_bound
      USE mod_inwave_swan
      USE mod_grid
!
# if defined EW_PERIODIC || defined NS_PERIODIC
      USE exchange_2d_mod
# endif
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
#  ifdef SOLVE3D
      USE mp_exchange_mod, ONLY : mp_exchange3d
#  endif
# endif
      USE set_2dfld_mod
# ifdef SOLVE3D
      USE set_3dfld_mod
# endif
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
!  Local variable declarations.
!
      integer :: ILB, IUB, JLB, JUB
      integer :: Imin, Imax, Jmin, Jmax
      integer :: i, j, d, my_tile, tidx
      integer :: nt1, nt2

      real(r8):: twopi, otwopi, cff1, cff2, cff3
      real(r8):: phiw, cw, dist, toff, alpha
      real(r8), parameter :: eps = 0.0001_r8

      logical :: update = .FALSE.

# include "set_bounds.h"

      twopi=2.0_r8*pi
      otwopi=1.0_r8/(2.0_r8*pi)
!
!     Compute the time index for the BC data. The Ampzeta array
!     is Insteps long, but it repeats every 1/df.
!
      tidx=MOD(iic(ng),WAVES(ng)%Insteps)+1
!
      IF (DOMAIN(ng)%Western_Edge(tile)) THEN
        IF (LBC(iwest,isAC3d,ng)%acquire) THEN
          DO d=1,ND
            DO j=Jstr,Jend
#ifdef VARY_ACBC
            cff1=0.5_r8*(WAVEP(ng)%cy(Istr-1,j  ,d)+                    &
                         WAVEP(ng)%cy(Istr-1,j+1,d))
            cff2=WAVEP(ng)%cx(Istr-1,j,d)
            cw=sqrt(cff1**2+cff2**2)
            dist=REAL(j,r8)/GRID(ng)%pn(Istr-1,j)*sin(WAVEG(ng)%wd(d)-  &
     &           GRID(ng)%angler(Istr-1,j))
            toff=ABS(dist/(cw+eps))
            toff=MIN(toff,2.0_r8*WAVES(ng)%dur)
            cff1=toff/dt(ng)
            nt1=tidx+INT(cff1)
            IF (nt1.le.0) nt1=MOD(nt1,WAVES(ng)%Insteps)
            IF (nt1.le.0) nt1=nt1+WAVES(ng)%Insteps
            nt2=nt1+1
            IF (nt2>WAVES(ng)%Insteps) nt2=nt2-WAVES(ng)%Insteps
            phiw=toff/dt(ng)-INT(toff/dt(ng))
            cff2=(1.0_r8-phiw)*WAVES(ng)%Ampzeta(nt1,d)+                &
     &           phiw*WAVES(ng)%Ampzeta(nt2,d)
            WAVEB(ng)%AC_west(j,d)=cff2*WAVEP(ng)%Tr(Istr-1,j,d)*otwopi
#else
            WAVEB(ng)%AC_west(j,d)=WAVES(ng)%Ampzeta(tidx,d)*           &
     &                             WAVEP(ng)%Tr(Istr-1,j,d)*otwopi
#endif
            END DO
          END DO
        END IF
      END IF
!
      IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
        IF (LBC(ieast,isAC3d,ng)%acquire) THEN
          DO d=1,ND
            DO j=Jstr,Jend
#ifdef VARY_ACBC
            cff1=0.5_r8*(WAVEP(ng)%cy(Iend+1,j  ,d)+                    &
                         WAVEP(ng)%cy(Iend+1,j+1,d))
            cff2=WAVEP(ng)%cx(Iend+1,j,d)
            cw=sqrt(cff1**2+cff2**2)
            dist=REAL(j,r8)/GRID(ng)%pn(Iend+1,j)*sin(WAVEG(ng)%wd(d)-  &
     &           GRID(ng)%angler(Iend+1,j))
            toff=ABS(dist/(cw+eps))
            toff=MIN(toff,2.0_r8*WAVES(ng)%dur)
            cff1=toff/dt(ng)
            nt1=tidx+INT(cff1)
            IF (nt1.le.0) nt1=MOD(nt1,WAVES(ng)%Insteps)
            IF (nt1.le.0) nt1=nt1+WAVES(ng)%Insteps
            nt2=nt1+1
            IF (nt2>WAVES(ng)%Insteps) nt2=nt2-WAVES(ng)%Insteps
            phiw=toff/dt(ng)-INT(toff/dt(ng))
            cff2=(1.0_r8-phiw)*WAVES(ng)%Ampzeta(nt1,d)+                &
     &           phiw*WAVES(ng)%Ampzeta(nt2,d)
            WAVEB(ng)%AC_east(j,d)=cff2*WAVEP(ng)%Tr(Iend+1,j,d)*otwopi
#else
            WAVEB(ng)%AC_east(j,d)=WAVES(ng)%Ampzeta(tidx,d)*           &
     &                             WAVEP(ng)%Tr(Iend+1,j,d)*otwopi
#endif
            END DO
          END DO
        END IF
      END IF
!
      IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
        IF (LBC(inorth,isAC3d,ng)%acquire) THEN
          DO d=1,ND
            DO i=Istr,Iend
#ifdef VARY_ACBC
            cff1=0.5_r8*(WAVEP(ng)%cx(i  ,Jend+1,d)+                    &
                         WAVEP(ng)%cx(i+1,Jend+1,d))
            cff2=WAVEP(ng)%cy(i  ,Jend+1,d)
            cw=sqrt(cff1**2+cff2**2)
            dist=REAL(i,r8)/GRID(ng)%pm(i,Jend+1)*sin(WAVEG(ng)%wd(d)-  &
     &           GRID(ng)%angler(i,Jend+1))
            toff=ABS(dist/(cw+eps))
            toff=MIN(toff,2.0_r8*WAVES(ng)%dur)
            cff1=toff/dt(ng)
            nt1=tidx+INT(cff1)
            IF (nt1.le.0) nt1=MOD(nt1,WAVES(ng)%Insteps)
            IF (nt1.le.0) nt1=nt1+WAVES(ng)%Insteps
            nt2=nt1+1
            IF (nt2>WAVES(ng)%Insteps) nt2=nt2-WAVES(ng)%Insteps
            phiw=toff/dt(ng)-INT(toff/dt(ng))
            cff2=(1.0_r8-phiw)*WAVES(ng)%Ampzeta(nt1,d)+                &
     &           phiw*WAVES(ng)%Ampzeta(nt2,d)
            WAVEB(ng)%AC_north(i,d)=cff2*WAVEP(ng)%Tr(i,Jend+1,d)*otwopi
#else
            WAVEB(ng)%AC_north(i,d)=WAVES(ng)%Ampzeta(tidx,d)*          &
     &                              WAVEP(ng)%Tr(i,Jend+1,d)*otwopi
#endif
            END DO
          END DO
        END IF
      END IF
!
      IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
        IF (LBC(isouth,isAC3d,ng)%acquire) THEN
          DO d=1,ND
            DO i=Istr,Iend
#ifdef VARY_ACBC
            cff1=0.5_r8*(WAVEP(ng)%cx(i  ,Jstr-1,d)+                    &
                         WAVEP(ng)%cx(i+1,Jstr-1,d))
            cff2=WAVEP(ng)%cy(i  ,Jstr-1,d)
            alpha=GRID(ng)%angler(i,Jstr-1)
            cw=cff1*cos(alpha)-cff2*sin(alpha)
            dist=REAL(i,r8)/GRID(ng)%pm(i,Jstr-1)*                      &
     &           sin(pi-WAVEG(ng)%wd(d))
            toff=ABS(dist/(cw+eps))
            cff1=toff/dt(ng)
            nt1=tidx+INT(cff1)
            nt1=MAX(MOD(nt1,WAVES(ng)%Insteps),1)
            nt2=nt1+1
            IF (nt2>WAVES(ng)%Insteps) nt2=nt2-WAVES(ng)%Insteps
            phiw=toff/dt(ng)-INT(toff/dt(ng))
            cff2=(1.0_r8-phiw)*WAVES(ng)%Ampzeta(nt1,d)+                &
     &           phiw*WAVES(ng)%Ampzeta(nt2,d)
            WAVEB(ng)%AC_south(i,d)=cff2*WAVEP(ng)%Tr(i,Jstr-1,d)*otwopi
#else
            WAVEB(ng)%AC_south(i,d)=WAVES(ng)%Ampzeta(tidx,d)*          &
     &                              WAVEP(ng)%Tr(i,Jstr-1,d)*otwopi
#endif
            END DO
          END DO
        END IF
      END IF

      RETURN
      END SUBROUTINE set_inwave_swan_data_tile
#else
      SUBROUTINE set_inwave_swan_data
      RETURN
      END SUBROUTINE set_inwave_swan_data
#endif
