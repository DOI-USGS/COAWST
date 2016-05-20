      SUBROUTINE prsgrd (ng, tile)
!
!svn $Id: prsgrd31.h 732 2008-09-07 01:55:51Z jcwarner $
!***********************************************************************
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                           Hernan G. Arango   !
!****************************************** Alexander F. Shchepetkin ***
!                                                                      !
!  This subroutine evaluates the  baroclinic  hydrostatic  pressure    !
!  gradient term using  the STANDARD density Jacobian  or  WEIGHTED    !
!  density Jacobian scheme of Song (1998). Both of these approaches    !
!  compute horizontal differences of density before of the vertical    !
!  integration.                                                        !
!                                                                      !
!  The pressure gradient terms (m4/s2) are loaded into right-hand-     !
!  side arrays "ru" and "rv".                                          !
!                                                                      !
!  Reference:                                                          !
!                                                                      !
!    Song, Y.T., 1998:  A general pressure gradient formulation for    !
!      numerical ocean models. Part I: Scheme design and diagnostic    !
!      analysis, Monthly Weather Rev., 126, 3213-3230.                 !
!                                                                      !
!***********************************************************************
!
      USE mod_param
#ifdef DIAGNOSTICS
      USE mod_diags
#endif
#ifdef ATM_PRESS
      USE mod_forces
#endif
      USE mod_grid
      USE mod_ocean
      USE mod_stepping
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
#include "tile.h"
!
#ifdef PROFILE
      CALL wclock_on (ng, iNLM, 23)
#endif
      CALL prsgrd_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  IminS, ImaxS, JminS, JmaxS,                     &
     &                  nrhs(ng),                                       &
#ifdef WET_DRY
     &                  GRID(ng)%umask_wet,                             &
     &                  GRID(ng)%vmask_wet,                             &
#endif
     &                  GRID(ng) % Hz,                                  &
     &                  GRID(ng) % om_v,                                &
     &                  GRID(ng) % on_u,                                &
     &                  GRID(ng) % z_r,                                 &
     &                  GRID(ng) % z_w,                                 &
     &                  OCEAN(ng) % rho,                                &
#ifdef WEC_VF
     &                  OCEAN(ng) % zetat,                              &
#endif
#ifdef ATM_PRESS
     &                  FORCES(ng) % Pair,                              &
#endif
#ifdef DIAGNOSTICS_UV
     &                  DIAGS(ng) % DiaRU,                              &
     &                  DIAGS(ng) % DiaRV,                              &
#endif
     &                  OCEAN(ng) % ru,                                 &
     &                  OCEAN(ng) % rv)
#ifdef PROFILE
      CALL wclock_off (ng, iNLM, 23)
#endif
      RETURN
      END SUBROUTINE prsgrd
!
!***********************************************************************
      SUBROUTINE prsgrd_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        IminS, ImaxS, JminS, JmaxS,               &
     &                        nrhs,                                     &
#ifdef WET_DRY
     &                        umask_wet, vmask_wet,                     &
#endif
     &                        Hz, om_v, on_u, z_r, z_w,                 &
     &                        rho,                                      &
#ifdef WEC_VF
     &                        zetat,                                    &
#endif
#ifdef ATM_PRESS
     &                        Pair,                                     &
#endif
#ifdef DIAGNOSTICS_UV
     &                        DiaRU, DiaRV,                             &
#endif
     &                        ru, rv)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nrhs

#ifdef ASSUMED_SHAPE
# ifdef WET_DRY
      real(r8), intent(in) :: umask_wet(LBi:,LBj:)
      real(r8), intent(in) :: vmask_wet(LBi:,LBj:)
# endif
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: om_v(LBi:,LBj:)
      real(r8), intent(in) :: on_u(LBi:,LBj:)
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
      real(r8), intent(in) :: rho(LBi:,LBj:,:)
# ifdef WEC_VF
      real(r8), intent(in) :: zetat(LBi:,LBj:)
# endif
# ifdef ATM_PRESS
      real(r8), intent(in) :: Pair(LBi:,LBj:)
# endif
# ifdef DIAGNOSTICS_UV
      real(r8), intent(inout) :: DiaRU(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: DiaRV(LBi:,LBj:,:,:,:)
# endif
      real(r8), intent(inout) :: ru(LBi:,LBj:,0:,:)
      real(r8), intent(inout) :: rv(LBi:,LBj:,0:,:)
#else
# ifdef WET_DRY
      real(r8), intent(in) :: umask_wet(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask_wet(LBi:UBi,LBj:UBj)
# endif
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: om_v(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: on_u(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: z_w(LBi:UBi,LBj:UBj,0:N(ng))
      real(r8), intent(in) :: rho(LBi:UBi,LBj:UBj,N(ng))
# ifdef WEC_VF
      real(r8), intent(in) :: zetat(LBi:UBi,LBj:UBj)
# endif
# ifdef ATM_PRESS
      real(r8), intent(in) :: Pair(LBi:UBi,LBj:UBj)
# endif
# ifdef DIAGNOSTICS_UV
      real(r8), intent(inout) :: DiaRU(LBi:UBi,LBj:UBj,N(ng),2,NDrhs)
      real(r8), intent(inout) :: DiaRV(LBi:UBi,LBj:UBj,N(ng),2,NDrhs)
# endif
      real(r8), intent(inout) :: ru(LBi:UBi,LBj:UBj,0:N(ng),2)
      real(r8), intent(inout) :: rv(LBi:UBi,LBj:UBj,0:N(ng),2)
#endif
!
!  Local variable declarations.
!
      integer :: i, j, k
      real(r8) :: fac, fac1, fac2, fac3
      real(r8) :: cff1, cff2, cff3, cff4
#ifdef WJ_GRADP
      real(r8) :: gamma
#endif

      real(r8), dimension(IminS:ImaxS) :: phie
      real(r8), dimension(IminS:ImaxS) :: phix

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Calculate pressure gradient in the XI-direction (m4/s2).
!-----------------------------------------------------------------------
!
!  Compute surface baroclinic pressure gradient.
!
#ifdef ATM_PRESS
      fac=100.0_r8/rho0
#endif
      fac1=0.5_r8*g/rho0
      fac2=1000.0_r8*g/rho0
      fac3=0.25_r8*g/rho0

      DO j=Jstr,Jend
        DO i=IstrU,Iend
          cff1=z_w(i  ,j,N(ng))-z_r(i  ,j,N(ng))+                       &
     &         z_w(i-1,j,N(ng))-z_r(i-1,j,N(ng))
          phix(i)=fac1*(rho(i,j,N(ng))-rho(i-1,j,N(ng)))*cff1
#ifdef WEC_VF
          phix(i)=phix(i)+zetat(i,j)-zetat(i-1,j)
#endif
#ifdef ATM_PRESS
          phix(i)=phix(i)+fac*(Pair(i,j)-Pair(i-1,j))
#endif
#ifdef RHO_SURF
          phix(i)=phix(i)+                                              &
     &            (fac2+fac1*(rho(i,j,N(ng))+rho(i-1,j,N(ng))))*        &
     &            (z_w(i,j,N(ng))-z_w(i-1,j,N(ng)))
#endif
          ru(i,j,N(ng),nrhs)=-0.5_r8*(Hz(i,j,N(ng))+Hz(i-1,j,N(ng)))*   &
     &                       phix(i)*on_u(i,j)
#ifdef WET_DRY
          ru(i,j,N(ng),nrhs)=ru(i,j,N(ng),nrhs)*umask_wet(i,j)
#endif
#ifdef DIAGNOSTICS_UV
          DiaRU(i,j,N(ng),nrhs,M3pgrd)=ru(i,j,N(ng),nrhs)
#endif
        END DO
!
!  Compute interior baroclinic pressure gradient.  Differentiate and
!  then vertically integrate.
!
        DO k=N(ng)-1,1,-1
          DO i=IstrU,Iend
#ifdef WJ_GRADP
            cff1=1.0_r8/((z_r(i  ,j,k+1)-z_r(i  ,j,k))*                 &
     &                   (z_r(i-1,j,k+1)-z_r(i-1,j,k)))
            cff2=z_r(i  ,j,k  )-z_r(i-1,j,k  )+                         &
     &           z_r(i  ,j,k+1)-z_r(i-1,j,k+1)
            cff3=z_r(i  ,j,k+1)-z_r(i  ,j,k  )-                         &
     &           z_r(i-1,j,k+1)+z_r(i-1,j,k  )
            gamma=0.125_r8*cff1*cff2*cff3

            cff1=(1.0_r8+gamma)*(rho(i,j,k+1)-rho(i-1,j,k+1))+          &
     &           (1.0_r8-gamma)*(rho(i,j,k  )-rho(i-1,j,k  ))
            cff2=rho(i,j,k+1)+rho(i-1,j,k+1)-                           &
     &           rho(i,j,k  )-rho(i-1,j,k  )
            cff3=z_r(i,j,k+1)+z_r(i-1,j,k+1)-                           &
     &           z_r(i,j,k  )-z_r(i-1,j,k  )
            cff4=(1.0_r8+gamma)*(z_r(i,j,k+1)-z_r(i-1,j,k+1))+          &
     &           (1.0_r8-gamma)*(z_r(i,j,k  )-z_r(i-1,j,k  ))
            phix(i)=phix(i)+                                            &
     &              fac3*(cff1*cff3-cff2*cff4)
#else
            cff1=rho(i,j,k+1)-rho(i-1,j,k+1)+                           &
     &           rho(i,j,k  )-rho(i-1,j,k  )
            cff2=rho(i,j,k+1)+rho(i-1,j,k+1)-                           &
     &           rho(i,j,k  )-rho(i-1,j,k  )
            cff3=z_r(i,j,k+1)+z_r(i-1,j,k+1)-                           &
     &           z_r(i,j,k  )-z_r(i-1,j,k  )
            cff4=z_r(i,j,k+1)-z_r(i-1,j,k+1)+                           &
     &           z_r(i,j,k  )-z_r(i-1,j,k  )
            phix(i)=phix(i)+                                            &
     &              fac3*(cff1*cff3-cff2*cff4)
#endif
            ru(i,j,k,nrhs)=-0.5_r8*(Hz(i,j,k)+Hz(i-1,j,k))*             &
     &                     phix(i)*on_u(i,j)
#ifdef WET_DRY
            ru(i,j,k,nrhs)=ru(i,j,k,nrhs)*umask_wet(i,j)
#endif
#ifdef DIAGNOSTICS_UV
            DiaRU(i,j,k,nrhs,M3pgrd)=ru(i,j,k,nrhs)
#endif
          END DO
        END DO
!
!-----------------------------------------------------------------------
!  Calculate pressure gradient in the ETA-direction (m4/s2).
!-----------------------------------------------------------------------
!
!  Compute surface baroclinic pressure gradient.
!
        IF (j.ge.JstrV) THEN
          DO i=Istr,Iend
            cff1=z_w(i,j  ,N(ng))-z_r(i,j  ,N(ng))+                     &
     &           z_w(i,j-1,N(ng))-z_r(i,j-1,N(ng))
            phie(i)=fac1*(rho(i,j,N(ng))-rho(i,j-1,N(ng)))*cff1
#ifdef WEC_VF
            phie(i)=phie(i)+zetat(i,j)-zetat(i,j-1)
#endif
#ifdef ATM_PRESS
            phie(i)=phie(i)+fac*(Pair(i,j)-Pair(i,j-1))
#endif
#ifdef RHO_SURF
            phie(i)=phie(i)+                                            &
     &              (fac2+fac1*(rho(i,j,N(ng))+rho(i,j-1,N(ng))))*      &
     &              (z_w(i,j,N(ng))-z_w(i,j-1,N(ng)))
#endif
            rv(i,j,N(ng),nrhs)=-0.5_r8*(Hz(i,j,N(ng))+Hz(i,j-1,N(ng)))* &
     &                         phie(i)*om_v(i,j)
#ifdef WET_DRY
            rv(i,j,N(ng),nrhs)=rv(i,j,N(ng),nrhs)*vmask_wet(i,j)
#endif
#ifdef DIAGNOSTICS_UV
            DiaRV(i,j,N(ng),nrhs,M3pgrd)=rv(i,j,N(ng),nrhs)
#endif
          END DO
!
!  Compute interior baroclinic pressure gradient.  Differentiate and
!  then vertically integrate.
!
          DO k=N(ng)-1,1,-1
            DO i=Istr,Iend
#ifdef WJ_GRADP
              cff1=1.0_r8/((z_r(i,j  ,k+1)-z_r(i,j  ,k))*               &
     &                     (z_r(i,j-1,k+1)-z_r(i,j-1,k)))
              cff2=z_r(i,j  ,k  )-z_r(i,j-1,k  )+                       &
     &             z_r(i,j  ,k+1)-z_r(i,j-1,k+1)
              cff3=z_r(i,j  ,k+1)-z_r(i,j  ,k  )-                       &
     &             z_r(i,j-1,k+1)+z_r(i,j-1,k  )
              gamma=0.125_r8*cff1*cff2*cff3

              cff1=(1.0_r8+gamma)*(rho(i,j,k+1)-rho(i,j-1,k+1))+        &
     &             (1.0_r8-gamma)*(rho(i,j,k  )-rho(i,j-1,k  ))
              cff2=rho(i,j,k+1)+rho(i,j-1,k+1)-                         &
     &             rho(i,j,k  )-rho(i,j-1,k  )
              cff3=z_r(i,j,k+1)+z_r(i,j-1,k+1)-                         &
     &             z_r(i,j,k  )-z_r(i,j-1,k  )
              cff4=(1.0_r8+gamma)*(z_r(i,j,k+1)-z_r(i,j-1,k+1))+        &
     &             (1.0_r8-gamma)*(z_r(i,j,k  )-z_r(i,j-1,k  ))
              phie(i)=phie(i)+                                          &
     &                fac3*(cff1*cff3-cff2*cff4)
#else
              cff1=rho(i,j,k+1)-rho(i,j-1,k+1)+                         &
     &             rho(i,j,k  )-rho(i,j-1,k  )
              cff2=rho(i,j,k+1)+rho(i,j-1,k+1)-                         &
     &             rho(i,j,k  )-rho(i,j-1,k  )
              cff3=z_r(i,j,k+1)+z_r(i,j-1,k+1)-                         &
     &             z_r(i,j,k  )-z_r(i,j-1,k  )
              cff4=z_r(i,j,k+1)-z_r(i,j-1,k+1)+                         &
     &             z_r(i,j,k  )-z_r(i,j-1,k  )
              phie(i)=phie(i)+                                          &
     &                fac3*(cff1*cff3-cff2*cff4)
#endif
              rv(i,j,k,nrhs)=-0.5_r8*(Hz(i,j,k)+Hz(i,j-1,k))*           &
     &                       phie(i)*om_v(i,j)
#ifdef WET_DRY
              rv(i,j,k,nrhs)=rv(i,j,k,nrhs)*vmask_wet(i,j)
#endif
#ifdef DIAGNOSTICS_UV
              DiaRV(i,j,k,nrhs,M3pgrd)=rv(i,j,k,nrhs)
#endif
            END DO
          END DO
        END IF
      END DO
      RETURN
      END SUBROUTINE prsgrd_tile
