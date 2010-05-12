#include "cppdefs.h"
      MODULE oi_update_mod
#ifdef ASSIMILATION
# ifdef EW_PERIODIC
#  define IU_RANGE Istr,Iend
#  define IV_RANGE Istr,Iend
# else
#  define IU_RANGE Istr,IendR
#  define IV_RANGE IstrR,IendR
# endif
# ifdef NS_PERIODIC
#  define JU_RANGE Jstr,Jend
#  define JV_RANGE Jstr,Jend
# else
#  define JU_RANGE JstrR,JendR
#  define JV_RANGE Jstr,JendR
# endif
!
!svn $Id: oi_update.F 733 2008-09-07 01:56:45Z jcwarner $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2010 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine assimilates fields as a  pointwise  linear  optimal    !
!  combination between model and observations. During initilization    !
!  (at first assimilation cycle),  the initial model error variance    !
!  has the same shape distribution as observations.                    !
!                                                                      !
!  Reference:                                                          !
!                                                                      !
!    Dombrowsky, E. and P. De May, 1992:  Continuous assimilation      !
!      in an open domain of the Northeast Atlantic 1. Methodology      !
!      and application to AtheA-88, JGR, 97, 9719-9731.                !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: oi_update

      CONTAINS
!
!***********************************************************************
      SUBROUTINE oi_update (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
      USE mod_obs
      USE mod_ocean
      USE mod_stepping
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
      CALL wclock_on (ng, iNLM, 2)
# endif
# ifdef ASSIMILATION_SSH
      CALL oi_ssh_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  IminS, ImaxS, JminS, JmaxS,                     &
     &                  knew(ng),                                       &
#  ifdef MASKING
     &                  GRID(ng) % rmask,                               &
#  endif
     &                  OBS(ng) % SSHobs,                               &
     &                  OBS(ng) % EobsSSH,                              &
     &                  OBS(ng) % EmodSSH,                              &
     &                  OCEAN(ng) % zeta)
# endif

# ifdef ASSIMILATION_T
      CALL oi_trc_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  IminS, ImaxS, JminS, JmaxS,                     &
     &                  nnew(ng),                                       &
#  ifdef MASKING
     &                  GRID(ng) % rmask,                               &
#  endif
     &                  OBS(ng) % Tobs,                                 &
     &                  OBS(ng) % EobsT,                                &
     &                  OBS(ng) % EmodT,                                &
     &                  OCEAN(ng) % t)
# endif

# if defined ASSIMILATION_UV || defined ASSIMILATION_UVsur
      CALL oi_vel_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  IminS, ImaxS, JminS, JmaxS,                     &
     &                  knew(ng), nnew(ng),                             &
#  ifdef MASKING
     &                  GRID(ng) % umask,                               &
     &                  GRID(ng) % vmask,                               &
#  endif
     &                  GRID(ng) % Hz,                                  &
     &                  OBS(ng) % Uobs,                                 &
     &                  OBS(ng) % Vobs,                                 &
     &                  OBS(ng) % EobsUV,                               &
     &                  OBS(ng) % EmodU,                                &
     &                  OBS(ng) % EmodV,                                &
#  ifndef UV_BAROCLINIC
     &                  OCEAN(ng) % ubar,                               &
     &                  OCEAN(ng) % vbar,                               &
#  endif
     &                  OCEAN(ng) % u,                                  &
     &                  OCEAN(ng) % v)
# endif

# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 2)
# endif
      RETURN
      END SUBROUTINE oi_update
# ifdef ASSIMILATION_SSH
!
!***********************************************************************
      SUBROUTINE oi_ssh_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        IminS, ImaxS, JminS, JmaxS,               &
     &                        knew,                                     &
#  ifdef MASKING
     &                        rmask,                                    &
#  endif
     &                        SSHobs, EobsSSH, EmodSSH,                 &
     &                        zeta)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: knew
!
#  ifdef ASSUMED_SHAPE
#   ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
#   endif
      real(r8), intent(in) :: SSHobs(LBi:,LBj:)
      real(r8), intent(in) :: EobsSSH(LBi:,LBj:)

      real(r8), intent(inout) :: EmodSSH(LBi:,LBj:)
      real(r8), intent(inout) :: zeta(LBi:,LBj:,:)
#  else
#   ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
#   endif
      real(r8), intent(in) :: SSHobs(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: EobsSSH(LBi:UBi,LBj:UBj)

      real(r8), intent(inout) :: EmodSSH(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: zeta(LBi:UBi,LBj:UBj,3)
#  endif
!
!  Local variable declarations.
!
      integer :: i, j

      real(r8), parameter :: eps = 1.0E-4_r8

      real(r8) :: arg, cff, cff1, cff2, delta
      real(r8) :: decay, mu, ratio, Aweight, Tval

      character (len=14) :: t_code

#  include "set_bounds.h"
!
!---------------------------------------------------------------------
!  Assimilate sea surface height data.
!---------------------------------------------------------------------
!
      IF (assi_SSH(ng).and.update_SSH(ng)) THEN
!
!  On first pass, initialize model error variance.
!
        IF (first_SSH(ng)) THEN
          delta=MAX(EobsSSHmax(ng)-EobsSSHmin(ng),eps)
          DO j=JU_RANGE
            DO i=IV_RANGE
              ratio=(EobsSSHmax(ng)-EobsSSH(i,j))/delta
              mu=MIN(1.0_r8,ratio)*Emod0(ng)
              cff1=1.0_r8-2.0_r8*mu
              cff2=(cor(ng)*cff1+                                       &
     &              SQRT(1.0_r8+cff1*cff1*(cor(ng)**2-1.0_r8)))/        &
     &             MAX(2.0_r8-2.0_r8*mu,eps)
              EmodSSH(i,j)=cff2*cff2*EobsSSH(i,j)
            END DO
          END DO
          IF (SOUTH_WEST_TEST)                                          &
     &      tSSHobs(2,ng)=tSSHobs(1,ng)
          IF (NORTH_EAST_TEST)                                          &
     &      first_SSH(ng)=.FALSE.
        END IF
!
!  Determine assimilation weights and meld model and observations.
!
        IF ((time(ng).le.tsSSHobs(ng)).and.                             &
     &      (tsSSHobs(ng).lt.(time(ng)+dt(ng)))) THEN
          arg=ABS(tSSHobs(1,ng)-tSSHobs(2,ng))/Tgrowth(ng)
          decay=2.0_r8*(1.0_r8-EXP(-arg*arg))
          DO j=JU_RANGE
            DO i=IV_RANGE
              EmodSSH(i,j)=EmodSSH(i,j)+decay
              cff1=cor(ng)*SQRT(EobsSSH(i,j)*EmodSSH(i,j))
              cff2=EobsSSH(i,j)+EmodSSH(i,j)-2.0_r8*cff1
              Aweight=(EmodSSH(i,j)-cff1)/MAX(cff2,eps)
              Aweight=MAX(0.0_r8,MIN(1.0_r8,Aweight))
              zeta(i,j,1)=(Aweight*SSHobs(i,j)+                         &
     &                     (1.0_r8-Aweight)*zeta(i,j,knew))
#  ifdef MASKING
              zeta(i,j,1)=zeta(i,j,1)*rmask(i,j)
#  endif
              zeta(i,j,2)=zeta(i,j,1)
              EmodSSH(i,j)=(1.0_r8-cor(ng))*EobsSSH(i,j)*               &
     &                                      EmodSSH(i,j)/cff2
            END DO
          END DO
          IF (NORTH_EAST_TEST) THEN
            tSSHobs(2,ng)=tSSHobs(1,ng)
            synchro_flag(ng)=.TRUE.
            update_SSH(ng)=.FALSE.
            Tval=tSSHobs(1,ng)*day2sec
            CALL time_string (Tval, t_code)
            IF (Master) WRITE (stdout,10) 'SSH', t_code
 10         FORMAT (' OI_UPDATE   - Assimilating ',a,' data,',t62,       &
     &             't = ',a)
          END IF
        END IF
      END IF
      RETURN
      END SUBROUTINE oi_ssh_tile
# endif /* ASSIMILATION_SSH */
# ifdef ASSIMILATION_T
!
!***********************************************************************
      SUBROUTINE oi_trc_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        IminS, ImaxS, JminS, JmaxS,               &
     &                        nnew,                                     &
#  ifdef MASKING
     &                        rmask,                                    &
#  endif
     &                        Tobs, EobsT, EmodT,                       &
     &                        t)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nnew
!
#  ifdef ASSUMED_SHAPE
#   ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
#   endif
      real(r8), intent(in) :: Tobs(LBi:,LBj:,:,:)
      real(r8), intent(in) :: EobsT(LBi:,LBj:,:,:)

      real(r8), intent(inout) :: EmodT(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: t(LBi:,LBj:,:,:,:)
#  else
#   ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
#   endif
      real(r8), intent(in) :: Tobs(LBi:UBi,LBj:UBj,N(ng),NT(ng))
      real(r8), intent(in) :: EobsT(LBi:UBi,LBj:UBj,N(ng),NT(ng))

      real(r8), intent(inout) :: EmodT(LBi:UBi,LBj:UBj,N(ng),NT(ng))
      real(r8), intent(inout) :: t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
#  endif
!
!  Local variable declarations.
!
      integer :: i, itrc, j, k

      real(r8), parameter :: eps = 1.0E-4_r8

      real(r8) :: arg, cff, cff1, cff2, delta
      real(r8) :: decay, mu, ratio, Aweight, Tval

      character (len=14) :: t_code

# include "set_bounds.h"
!
!---------------------------------------------------------------------
!  Assimilate tracer data.
!---------------------------------------------------------------------
!
      DO itrc=1,NT(ng)
        IF (assi_T(itrc,ng).and.update_T(itrc,ng)) THEN
!
!  On first pass, initialize model error variance.
!
          IF (first_T(itrc,ng)) THEN
            delta=MAX(EobsTmax(itrc,ng)-EobsTmin(itrc,ng),eps)
            DO k=1,N(ng)
              DO j=JU_RANGE
                DO i=IV_RANGE
                  ratio=(EobsTmax(itrc,ng)-EobsT(i,j,k,itrc))/delta
                  mu=MIN(1.0_r8,ratio)*Emod0(ng)
                  cff1=1.0_r8-2.0_r8*mu
                  cff2=(cor(ng)*cff1+                                   &
     &                  SQRT(1.0_r8+cff1*cff1*(cor(ng)**2-1.0_r8)))/    &
     &                 MAX(2.0_r8-2.0_r8*mu,eps)
                  EmodT(i,j,k,itrc)=cff2*cff2*EobsT(i,j,k,itrc)
                END DO
              END DO
            END DO
            IF (SOUTH_WEST_TEST)                                        &
     &        tTobs(2,itrc,ng)=tTobs(1,itrc,ng)
            IF (NORTH_EAST_TEST)                                        &
     &        first_T(itrc,ng)=.FALSE.
          END IF
!
!  Determine assimilation weights and meld model and observations.
!
          IF ((time(ng).le.tsTobs(itrc,ng)).and.                        &
     &        (tsTobs(itrc,ng).lt.(time(ng)+dt(ng)))) THEN
            arg=ABS(tTobs(1,itrc,ng)-tTobs(2,itrc,ng))/Tgrowth(ng)
            decay=2.0_r8*(1.0_r8-EXP(-arg*arg))
            DO k=1,N(ng)
              DO j=JU_RANGE
                DO i=IV_RANGE
                  EmodT(i,j,k,itrc)=EmodT(i,j,k,itrc)+decay
                  cff1=cor(ng)*SQRT(EobsT(i,j,k,itrc)*EmodT(i,j,k,itrc))
                  cff2=EobsT(i,j,k,itrc)+EmodT(i,j,k,itrc)-2.0_r8*cff1
                  Aweight=(EmodT(i,j,k,itrc)-cff1)/MAX(cff2,eps)
                  Aweight=MAX(0.0_r8,MIN(1.0_r8,Aweight))
                  t(i,j,k,1,itrc)=(Aweight*Tobs(i,j,k,itrc)+            &
     &                             (1.0_r8-Aweight)*t(i,j,k,nnew,itrc))
#  ifdef MASKING
                  t(i,j,k,1,itrc)=t(i,j,k,1,itrc)*rmask(i,j)
#  endif
                  t(i,j,k,2,itrc)=t(i,j,k,1,itrc)
                  EmodT(i,j,k,itrc)=(1.0_r8-cor(ng))*                   &
     &                              EobsT(i,j,k,itrc)*EmodT(i,j,k,itrc)/&
     &                              cff2
                END DO
              END DO
            END DO
            IF (NORTH_EAST_TEST) THEN
              tTobs(2,itrc,ng)=tTobs(1,itrc,ng)
              synchro_flag(ng)=.TRUE.
              update_T(itrc,ng)=.FALSE.
              ntfirst(ng)=iic(ng)+1
              Tval=tTobs(1,itrc,ng)*day2sec
              CALL time_string (Tval, t_code)
              IF (update_SST(ng).and.(itrc.eq.itemp)) THEN
                IF (Master) WRITE (stdout,10) 'SST', t_code
 10             FORMAT (' OI_UPDATE   - Assimilating ',a,' data,',t62,  &
     &                  't = ',a)
              ELSE
                IF (Master) WRITE (stdout,20) 'TRACER', itrc, t_code
 20             FORMAT (' OI_UPDATE   - Assimilating ',a,1x,i2.2,       &
     &                  ' data,',t62,'t = ',a)
              END IF
            END IF
          END IF
        END IF
      END DO
      RETURN
      END SUBROUTINE oi_trc_tile
# endif /* ASSIMILATION_T */
# if defined ASSIMILATION_UV || defined ASSIMILATION_UVsur
!
!***********************************************************************
      SUBROUTINE oi_vel_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        IminS, ImaxS, JminS, JmaxS,               &
     &                        knew, nnew,                               &
#  ifdef MASKING
     &                        umask, vmask,                             &
#  endif
     &                        Hz,                                       &
     &                        Uobs, Vobs, EobsUV, EmodU, EmodV          &
#  ifndef UV_BAROCLINIC
     &                        ubar, vbar,                               &
#  endif
     &                        u, v)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: knew, nnew
!
#  ifdef ASSUMED_SHAPE
#   ifdef MASKING
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#   endif
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: Uobs(LBi:,LBj:,:)
      real(r8), intent(in) :: Vobs(LBi:,LBj:,:)
      real(r8), intent(in) :: EobsUV(LBi:,LBj:,:)

      real(r8), intent(inout) :: EmodU(LBi:,LBj:,:)
      real(r8), intent(inout) :: Emodv(LBi:,LBj:,:)
#   ifndef UV_BAROCLINIC
      real(r8), intent(inout) :: ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: vbar(LBi:,LBj:,:)
#   endif
      real(r8), intent(inout) :: u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: v(LBi:,LBj:,:,:)
#  else
#   ifdef MASKING
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#   endif
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: Uobs(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: Vobs(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: EobsUV(LBi:UBi,LBj:UBj,N(ng))

      real(r8), intent(inout) :: EmodU(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: Emodv(LBi:UBi,LBj:UBj,N(ng))
#   ifndef UV_BAROCLINIC
      real(r8), intent(inout) :: ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: vbar(LBi:UBi,LBj:UBj,3)
#   endif
      real(r8), intent(inout) :: u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: v(LBi:UBi,LBj:UBj,N(ng),2)
#  endif
!
!  Local variable declarations.
!
      integer :: i, j, k

      real(r8), parameter :: eps = 1.0E-4_r8

      real(r8) :: Eobs, arg, cff, cff1, cff2, delta
      real(r8) :: decay, mu, ratio, Aweight, Tval

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: Uwrk
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: Vwrk

      character (len=14) :: t_code

#  include "set_bounds.h"
!
!---------------------------------------------------------------------
!  Assimilate velocity data.
!---------------------------------------------------------------------
!
      IF ((assi_UVsur(ng).or.assi_UV(ng)).and.update_UV(ng)) THEN
!
!  On first pass, initialize model error variance.
!
        IF (first_UV(ng)) THEN
          delta=MAX(EobsUVmax(ng)-EobsUVmin(ng),eps)
          DO k=1,N(ng)
            DO j=JU_RANGE
              DO i=IU_RANGE
                Eobs=0.5_r8*(EobsUV(i-1,j,k)+EobsUV(i,j,k))
                ratio=(EobsUVmax(ng)-Eobs)/delta
                mu=MIN(1.0_r8,ratio)*Emod0(ng)
                cff1=1.0_r8-2.0_r8*mu
                cff2=(cor(ng)*cff1+                                     &
     &                SQRT(1.0_r8+cff1*cff1*(cor(ng)**2-1.0_r8)))/      &
     &               MAX(2.0_r8-2.0_r8*mu,eps)
                EmodU(i,j,k)=cff2*cff2*Eobs
              END DO
            END DO
            DO j=JV_RANGE
              DO i=IV_RANGE
                Eobs=0.5_r8*(EobsUV(i,j-1,k)+EobsUV(i,j,k))
                ratio=(EobsUVmax(ng)-Eobs)/delta
                mu=MIN(1.0_r8,ratio)*Emod0(ng)
                cff1=1.0_r8-2.0_r8*mu
                cff2=(cor(ng)*cff1+                                     &
     &                SQRT(1.0_r8+cff1*cff1*(cor(ng)**2-1.0_r8)))/      &
     &               MAX(2.0_r8-2.0_r8*mu,eps)
                EmodU(i,j,k)=cff2*cff2*Eobs
              END DO
            END DO
          END DO
          IF (SOUTH_WEST_TEST)                                          &
     &      tVobs(2,ng)=tVobs(1,ng)
          IF (NORTH_EAST_TEST)                                          &
     &      first_UV(ng)=.FALSE.
        END IF
!
!  Determine assimilation weights and meld model and observations.
!
        IF ((time(ng).le.tsVobs).and.                                   &
     &      (tsVobs(ng).lt.(time(ng)+dt(ng)))) THEN
          DO j=JU_RANGE
            DO i=IU_RANGE
              Uwrk(i,j,1)=0.0_r8
              Uwrk(i,j,2)=0.5_r8*(Hz(i-1,j,N(ng))+Hz(i,j,N(ng)))
            END DO
          END DO
          DO j=JV_RANGE
            DO i=IV_RANGE
              Vwrk(i,j,1)=0.0_r8
              Vwrk(i,j,2)=0.5_r8*(Hz(i,j-1,N(ng))+Hz(i,j,N(ng)))
            END DO
          END DO
#  ifdef UV_BAROCLINIC
          DO k=1,N(ng)
            DO j=JU_RANGE
              DO i=IU_RANGE
                cff=0.5_r8*(Hz(i-1,j,k)+Hz(i,j,k))
                Uwrk(i,j,1)=Uwrk(i,j,1)+cff*Uobs(i,j,k)
              END DO
            END DO
            DO j=JV_RANGE
              DO i=IV_RANGE
                cff=0.5_r8*(Hz(i,j-1,k)+Hz(i,j,k))
                Vwrk(i,j,1)=Vwrk(i,j,1)+cff*Vobs(i,j,k)
              END DO
            END DO
          END DO
          DO j=JU_RANGE
            DO i=IU_RANGE
              Uwrk(i,j,1)=Uwrk(i,j,1)/Uwrk(i,j,2)
            END DO
          END DO
          DO j=JV_RANGE
            DO i=IV_RANGE
              Vwrk(i,j,1)=Vwrk(i,j,1)/Vwrk(i,j,2)
            END DO
          END DO
#  endif
!
          arg=ABS(tVobs(1)-tVobs(2))/Tgrowth(ng)
          decay=2.0_r8*(1.0_r8-EXP(-arg*arg))
          DO k=1,N(ng)
            DO j=JU_RANGE
              DO i=IU_RANGE
                EmodU(i,j,k)=EmodU(i,j,k)+decay
                Eobs=0.5_r8*(EobsUV(i-1,j,k)+EobsUV(i,j,k))
                cff1=cor(ng)*SQRT(Eobs*EmodU(i,j,k))
                cff2=Eobs+EmodU(i,j,k)-2.0_r8*cff1
                Aweight=(EmodU(i,j,k)-cff1)/MAX(cff2,eps)
                Aweight=MAX(0.0_r8,MIN(1.0_r8,Aweight))
#  ifdef UV_BAROCLINIC
                u(i,j,k,1)=(ubar(i,j,knew)+                             &
     &                      Aweight*(Uobs(i,j,k)-Uwrk(i,j,1))+          &
     &                      (1.0_r8-Aweight)*(u(i,j,k,nnew)-            &
     &                                        ubar(i,j,knew)))
#   ifdef MASKING
                u(i,j,k,1)=u(i,j,k,1)*umask(i,j)
#   endif
                u(i,j,k,2)=u(i,j,k,1)
#  else
                u(i,j,k,1)=(Aweight*Uobs(i,j,k)+                        &
     &                      (1.0_r8-Aweight)*u(i,j,k,nnew))
#   ifdef MASKING
                u(i,j,k,1)=u(i,j,k,1)**umask(i,j)
#   endif
                u(i,j,k,2)=u(i,j,k,1)
                cff=0.5_r8*(Hz(i-1,j,k)+Hz(i,j,k))
                Uwrk(i,j,1)=Uwrk(i,j,1)+cff*u(i,j,k,1)
                Uwrk(i,j,2)=Uwrk(i,j,2)+cff
#  endif
                EmodU(i,j,k)=(1.0_r8-cor(ng))*Eobs*EmodU(i,j,k)/cff2
              END DO
            END DO
            DO j=JV_RANGE
              DO i=IV_RANGE
                EmodV(i,j,k)=EmodV(i,j,k)+decay
                Eobs=0.5_r8*(EobsUV(i,j-1,k)+EobsUV(i,j,k))
                cff1=cor(ng)*SQRT(Eobs*EmodV(i,j,k))
                cff2=Eobs+EmodV(i,j,k)-2.0_r8*cff1
                Aweight=(EmodV(i,j,k)-cff1)/MAX(cff2,eps)
                Aweight=MAX(0.0_r8,MIN(1.0_r8,Aweight))
#  ifdef UV_BAROCLINIC
                v(i,j,k,1)=(vbar(i,j,1)+                                &
     &                      Aweight*(Vobs(i,j,k)-Vwrk(i,j,1))+          &
     &                      (1.0_r8-Aweight)*(v(i,j,k,nnew)-            &
     &                                        vbar(i,j,1)))
#   ifdef MASKING
                v(i,j,k,1)=v(i,j,k,1)*vmask(i,j)
#   endif
                v(i,j,k,2)=v(i,j,k,1)
#  else
                v(i,j,k,1)=(Aweight*Vobs(i,j,k)+                        &
     &                      (1.0_r8-Aweight)*v(i,j,k,nnew))
#   ifdef MASKING
                v(i,j,k,1)=v(i,j,k,1)*vmask(i,j)
#   endif
                v(i,j,k,2)=v(i,j,k,1)
                cff=0.5_r8*(Hz(i,j-1,k)+Hz(i,j,k))
                Vwrk(i,j,1)=Vwrk(i,j,1)+cff*v(i,j,k,1)
                Vwrk(i,j,2)=Vwrk(i,j,2)+cff
#  endif
                EmodV(i,j,k)=(1.0_r8-cor(ng))*Eobs*EmodV(i,j,k)/cff2
              END DO
            END DO
          END DO
          IF (NORTH_EAST_TEST) THEN
            tVobs(2,ng)=tVobs(1,ng)
            synchro_flag(ng)=.TRUE.
            ntfirst(ng)=iic(ng)+1
            Tval=tVobs(1,ng)*day2sec
            CALL time_string (Tval, t_code)
            IF (update_UVsur(ng)) THEN
              update_UVsur(ng)=.FALSE.
              IF (Master) WRITE (stdout,10) 'UVsur', t_code
 10           FORMAT (' OI_UPDATE   - Assimilating ',a,' data,',t62,    &
     &                't = ',a)
            ELSE
              update_UV(ng)=.FALSE.
              IF (Master) WRITE (stdout,10) 'UV', tVobs(1,ng)
            END IF
          END IF
#  ifndef UV_BAROCLINIC
!
!  Set barotropic momentum to the vertically integrated values.
!
          DO j=JU_RANGE
            DO i=IU_RANGE
              ubar(i,j,1)=Uwrk(i,j,1)/Uwrk(i,j,2)
#   ifdef MASKING
              ubar(i,j,1)=ubar(i,j,1)*umask(i,j)
#   endif
              ubar(i,j,2)=ubar(i,j,1)
            END DO
          END DO
          DO j=JV_RANGE
            DO i=IV_RANGE
              vbar(i,j,1)=Vwrk(i,j,1)/Vwrk(i,j,2)
#   ifdef MASKING
              vbar(i,j,1)=vbar(i,j,1)*vmask(i,j)
#   endif
              vbar(i,j,2)=vbar(i,j,1)
            END DO
          END DO
#  endif
        END IF
      END IF
      RETURN
      END SUBROUTINE oi_vel_tile
# endif /* ASSIMILATION_UV || ASSIMILATION_UVsur */
#endif
      END MODULE oi_update_mod
