#include "cppdefs.h"
      MODULE ini_fields_mod
#ifdef NONLINEAR
!
!svn $Id: ini_fields.F 838 2008-11-17 04:22:18Z jcwarner $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine initializes other time levels for 2D fields. It also   !
!  couples 3D and 2D momentum equations:  it initializes 2D momentum   !
!  (ubar,vbar) to the vertical integral of initial 3D momentum (u,v).  !
!                                                                      !
!=======================================================================
!
      implicit none
!
      PRIVATE
      PUBLIC :: ini_fields
      PUBLIC :: ini_zeta
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE ini_fields (ng, tile, model)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
# ifdef SOLVE3D
      USE mod_coupling
# endif
      USE mod_ocean
# if defined SOLVE3D && (defined SEDIMENT || defined BBL_MODEL)
      USE mod_sedbed
# endif
# ifdef ICE_MODEL
      USE mod_ice
# endif
      USE mod_stepping
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
!
!  Local variable declarations.
!
# include "tile.h"
!
# ifdef PROFILE
      CALL wclock_on (ng, iNLM, 2)
# endif
      CALL ini_fields_tile (ng, tile, model,                            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      kstp(ng), krhs(ng), knew(ng),               &
# ifdef SOLVE3D
     &                      nstp(ng), nnew(ng),                         &
# endif
# ifdef MASKING
     &                      GRID(ng) % rmask,                           &
     &                      GRID(ng) % umask,                           &
     &                      GRID(ng) % vmask,                           &
# endif
# ifdef WET_DRY
     &                      GRID(ng) % rmask_wet,                       &
     &                      GRID(ng) % umask_wet,                       &
     &                      GRID(ng) % vmask_wet,                       &
# endif
# ifdef PERFECT_RESTART
#  ifdef SOLVE3D
     &                      OCEAN(ng) % ru,                             &
     &                      OCEAN(ng) % rv,                             &
#  endif
     &                      OCEAN(ng) % rubar,                          &
     &                      OCEAN(ng) % rvbar,                          &
     &                      OCEAN(ng) % rzeta,                          &
# endif
# ifdef SOLVE3D
#  if defined SEDIMENT || defined BBL_MODEL
     &                      SEDBED(ng) % bottom,                        &
#  endif
#  if defined SEDIMENT
     &                      SEDBED(ng) % bed,                           &
     &                      SEDBED(ng) % bed_frac,                      &
     &                      SEDBED(ng) % bed_mass,                      &
#  endif
     &                      GRID(ng) % Hz,                              &
     &                      OCEAN(ng) % t,                              &
     &                      OCEAN(ng) % u,                              &
     &                      OCEAN(ng) % v,                              &
#  ifdef ICE_MODEL
     &                      ICE(ng) % ui,                               &
     &                      ICE(ng) % vi,                               &
     &                      ICE(ng) % ai,                               &
     &                      ICE(ng) % hi,                               &
     &                      ICE(ng) % hsn,                              &
     &                      ICE(ng) % ti,                               &
#   ifdef MELT_PONDS
     &                      ICE(ng) % apond,                            &
     &                      ICE(ng) % hpond,                            &
#   endif
     &                      ICE(ng) % ageice,                           &
#  endif
# endif
     &                      OCEAN(ng) % ubar,                           &
     &                      OCEAN(ng) % vbar,                           &
     &                      OCEAN(ng) % zeta)
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 2)
# endif
      RETURN
      END SUBROUTINE ini_fields
!
!***********************************************************************
      SUBROUTINE ini_fields_tile (ng, tile, model,                      &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            IminS, ImaxS, JminS, JmaxS,           &
     &                            kstp, krhs, knew,                     &
# ifdef SOLVE3D
     &                            nstp, nnew,                           &
# endif
# ifdef MASKING
     &                            rmask, umask, vmask,                  &
# endif
# ifdef WET_DRY
     &                            rmask_wet, umask_wet, vmask_wet,      &
# endif
# ifdef PERFECT_RESTART
#  ifdef SOLVE3D
     &                            ru, rv,                               &
#  endif
     &                            rubar, rvbar, rzeta,                  &
# endif
# ifdef SOLVE3D
#  if defined SEDIMENT || defined BBL_MODEL
     &                            bottom,                               &
#  endif
#  if defined SEDIMENT
     &                            bed, bed_frac, bed_mass,              &
#  endif
     &                            Hz,                                   &
     &                            t, u, v,                              &
#  ifdef ICE_MODEL
     &                            ui, vi, ai, hi, hsn,                  &
     &                            ti, ageice,                           &
#   ifdef MELT_PONDS
     &                            apond, hpond,                         &
#   endif
#  endif
# endif
     &                            ubar, vbar, zeta)
!***********************************************************************
!
      USE mod_param
      USE mod_ncparam
      USE mod_scalars
# if defined SEDIMENT || defined BBL_MODEL
      USE mod_sediment
# endif
!
      USE exchange_2d_mod
# ifdef SOLVE3D
      USE exchange_3d_mod
# endif
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
#  ifdef SOLVE3D
      USE mp_exchange_mod, ONLY : mp_exchange3d, mp_exchange4d
#  else
#   ifdef PERFECT_RESTART
      USE mp_exchange_mod, ONLY : mp_exchange3d
#   endif
#  endif
# endif
# ifdef SOLVE3D
      USE t3dbc_mod, ONLY : t3dbc_tile
#  ifdef T_PASSIVE
      USE pt3dbc_mod, ONLY : pt3dbc_tile
#  endif
      USE u3dbc_mod, ONLY : u3dbc_tile
      USE v3dbc_mod, ONLY : v3dbc_tile
# endif
      USE u2dbc_mod, ONLY : u2dbc_tile
      USE v2dbc_mod, ONLY : v2dbc_tile
# ifdef ICE_MODEL
      USE ice_limit_mod, ONLY : ice_limit
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: kstp, krhs, knew
# ifdef SOLVE3D
      integer, intent(in) :: nstp, nnew
# endif
!
# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#  endif
#  ifdef WET_DRY
      real(r8), intent(in)    :: rmask_wet(LBi:,LBj:)
      real(r8), intent(inout) :: umask_wet(LBi:,LBj:)
      real(r8), intent(inout) :: vmask_wet(LBi:,LBj:)
#  endif
#  ifdef SOLVE3D
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
#  endif
      real(r8), intent(in) :: zeta(LBi:,LBj:,:)
#  ifdef SOLVE3D
#   if defined SEDIMENT || defined BBL_MODEL
      real(r8), intent(inout) :: bottom(LBi:,LBj:,:)
#   endif
#   if defined SEDIMENT
      real(r8), intent(inout) :: bed(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: bed_frac(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: bed_mass(LBi:,LBj:,:,:,:)
#   endif
#  endif
#  ifdef PERFECT_RESTART
#   ifdef SOLVE3D
      real(r8), intent(inout) :: ru(LBi:,LBj:,0:,:)
      real(r8), intent(inout) :: rv(LBi:,LBj:,0:,:)
#   endif
      real(r8), intent(inout) :: rubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: rvbar(LBi:,LBj:,:)
      real(r8), intent(inout) :: rzeta(LBi:,LBj:,:)
#  endif
#  ifdef SOLVE3D
      real(r8), intent(inout) :: t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: v(LBi:,LBj:,:,:)
#   ifdef ICE_MODEL
      real(r8), intent(inout) :: ui(LBi:,LBj:,:)
      real(r8), intent(inout) :: vi(LBi:,LBj:,:)
      real(r8), intent(inout) :: ai(LBi:,LBj:,:)
      real(r8), intent(inout) :: hi(LBi:,LBj:,:)
      real(r8), intent(inout) :: hsn(LBi:,LBj:,:)
      real(r8), intent(inout) :: ti(LBi:,LBj:,:)
      real(r8), intent(inout) :: ageice(LBi:,LBj:,:)
#    ifdef MELT_PONDS
      real(r8), intent(inout) :: apond(LBi:,LBj:,:)
      real(r8), intent(inout) :: hpond(LBi:,LBj:,:)
#    endif
#   endif
#  endif
      real(r8), intent(inout) :: ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: vbar(LBi:,LBj:,:)

# else

#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#  endif
#  ifdef WET_DRY
      real(r8), intent(in)    :: rmask_wet(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: umask_wet(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: vmask_wet(LBi:UBi,LBj:UBj)
#  endif
#  ifdef SOLVE3D
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))
#  endif
      real(r8), intent(in) :: zeta(LBi:UBi,LBj:UBj,3)
#  ifdef SOLVE3D
#   if defined SEDIMENT || defined BBL_MODEL
      real(r8), intent(inout) :: bottom(LBi:UBi,LBj:UBj,MBOTP)
#   endif
#   if defined SEDIMENT
      real(r8), intent(inout) :: bed(LBi:UBi,LBj:UBj,Nbed,MBEDP)
      real(r8), intent(inout) :: bed_frac(LBi:UBi,LBj:UBj,Nbed,NST)
      real(r8), intent(inout) :: bed_mass(LBi:UBi,LBj:UBj,Nbed,2,NST)
#   endif
#  endif
#  ifdef PERFECT_RESTART
#   ifdef SOLVE3D
      real(r8), intent(inout) :: ru(LBi:UBi,LBj:UBj,0:N(ng),2)
      real(r8), intent(inout) :: rv(LBi:UBi,LBj:UBj,0:N(ng),2)
#   endif
      real(r8), intent(inout) :: rubar(LBi:UBi,LBj:UBj,2)
      real(r8), intent(inout) :: rvbar(LBi:UBi,LBj:UBj,2)
      real(r8), intent(inout) :: rzeta(LBi:UBi,LBj:UBj,2)
#  endif
#  ifdef SOLVE3D
      real(r8), intent(inout) :: t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: v(LBi:UBi,LBj:UBj,N(ng),2)
#   ifdef ICE_MODEL
      real(r8), intent(inout) :: ui(LBi:UBi,LBj:UBj,2)
      real(r8), intent(inout) :: vi(LBi:UBi,LBj:UBj,2)
      real(r8), intent(inout) :: ai(LBi:UBi,LBj:UBj,2)
      real(r8), intent(inout) :: hi(LBi:UBi,LBj:UBj,2)
      real(r8), intent(inout) :: hsn(LBi:UBi,LBj:UBj,2)
      real(r8), intent(inout) :: ti(LBi:UBi,LBj:UBj,2)
      real(r8), intent(inout) :: ageice(LBi:UBi,LBj:UBj,2)
#    ifdef MELT_PONDS
      real(r8), intent(inout) :: apond(LBi:UBi,LBj:UBj,2)
      real(r8), intent(inout) :: hpond(LBi:UBi,LBj:UBj,2)
#    endif
#   endif
#  endif
      real(r8), intent(inout) :: ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: vbar(LBi:UBi,LBj:UBj,3)
# endif
!
!  Local variable declarations.
!
      integer :: i, ic, itrc, j, k
# if defined SEDIMENT || defined BBL_MODEL
      integer :: ised
# endif

      real(r8) :: cff1, cff2, cff3, cff4, cff5, cff6, cff7
# ifdef SOLVE3D
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: CF
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: DC
# endif

# include "set_bounds.h"

# ifdef SOLVE3D
!
!-----------------------------------------------------------------------
!  If not perfect restart, initialize other time levels for 3D momentum.
!-----------------------------------------------------------------------
!
      IF (.not.PerfectRST(ng)) THEN
        DO j=JstrB,JendB
          DO k=1,N(ng)
            DO i=IstrM,IendB
              cff1=u(i,j,k,nstp)
#  ifdef MASKING
              cff1=cff1*umask(i,j)
#  endif
#  ifdef WET_DRY
              cff1=cff1*umask_wet(i,j)
#  endif
              u(i,j,k,nstp)=cff1
              u(i,j,k,nnew)=cff1
            END DO
          END DO
!
          IF (j.ge.JstrM) THEN
            DO k=1,N(ng)
              DO i=IstrB,IendB
                cff2=v(i,j,k,nstp)
#  ifdef MASKING
                cff2=cff2*vmask(i,j)
#  endif
#  ifdef WET_DRY
                cff2=cff2*vmask_wet(i,j)
#  endif
                v(i,j,k,nstp)=cff2
                v(i,j,k,nnew)=cff2
              END DO
            END DO
          END IF
        END DO
!
!  Apply boundary conditions.
!
        CALL u3dbc_tile (ng, tile,                                      &
     &                   LBi, UBi, LBj, UBj, N(ng),                     &
     &                   IminS, ImaxS, JminS, JmaxS,                    &
     &                   nstp, nstp,                                    &
     &                   u)
        CALL v3dbc_tile (ng, tile,                                      &
     &                   LBi, UBi, LBj, UBj, N(ng),                     &
     &                   IminS, ImaxS, JminS, JmaxS,                    &
     &                   nstp, nstp,                                    &
     &                   v)

        CALL u3dbc_tile (ng, tile,                                      &
     &                   LBi, UBi, LBj, UBj, N(ng),                     &
     &                   IminS, ImaxS, JminS, JmaxS,                    &
     &                   nstp, nnew,                                    &
     &                   u)
        CALL v3dbc_tile (ng, tile,                                      &
     &                   LBi, UBi, LBj, UBj, N(ng),                     &
     &                   IminS, ImaxS, JminS, JmaxS,                    &
     &                   nstp, nnew,                                    &
     &                   v)
      END IF
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_u3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          u(:,:,:,nstp))
        CALL exchange_v3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          v(:,:,:,nstp))

        CALL exchange_u3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          u(:,:,:,nnew))
        CALL exchange_v3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          v(:,:,:,nnew))
      END IF

#  ifdef DISTRIBUTE
!
      CALL mp_exchange3d (ng, tile, model, 4,                           &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    u(:,:,:,nstp), v(:,:,:,nstp),                 &
     &                    u(:,:,:,nnew), v(:,:,:,nnew))
#  endif
# endif

# ifdef SOLVE3D
!
!-----------------------------------------------------------------------
!  If not perfect restart, compute vertically-integrated momentum
!  (ubar, vbar) from initial 3D momentum (u, v).
!-----------------------------------------------------------------------
!
!  Here DC(i,1:N) are the grid cell thicknesses, DC(i,0) is the total
!  depth of the water column, and CF(i,0) is the vertical integral.
!
      IF (.not.PerfectRST(ng)) THEN
        DO j=JstrB,JendB
          DO i=IstrM,IendB
            DC(i,0)=0.0_r8
            CF(i,0)=0.0_r8
          END DO
          DO k=1,N(ng)
            DO i=IstrM,IendB
              DC(i,k)=0.5_r8*(Hz(i,j,k)+Hz(i-1,j,k))
              DC(i,0)=DC(i,0)+DC(i,k)
              CF(i,0)=CF(i,0)+DC(i,k)*u(i,j,k,nstp)
            END DO
          END DO
          DO i=IstrM,IendB
            cff1=1.0_r8/DC(i,0)
            cff2=CF(i,0)*cff1
#  ifdef MASKING
            cff2=cff2*umask(i,j)
#  endif
#  ifdef WET_DRY
            cff2=cff2*umask_wet(i,j)
#  endif
            ubar(i,j,kstp)=cff2
            ubar(i,j,knew)=cff2
          END DO
!
          IF (j.ge.JstrM) THEN
            DO i=IstrB,IendB
              DC(i,0)=0.0_r8
              CF(i,0)=0.0_r8
            END DO
            DO k=1,N(ng)
              DO i=IstrB,IendB
                DC(i,k)=0.5_r8*(Hz(i,j,k)+Hz(i,j-1,k))
                DC(i,0)=DC(i,0)+DC(i,k)
                CF(i,0)=CF(i,0)+DC(i,k)*v(i,j,k,nstp)
              END DO
            END DO
            DO i=IstrB,IendB
              cff1=1.0_r8/DC(i,0)
              cff2=CF(i,0)*cff1
#  ifdef MASKING
              cff2=cff2*vmask(i,j)
#  endif
#  ifdef WET_DRY
              cff2=cff2*vmask_wet(i,j)
#  endif
              vbar(i,j,kstp)=cff2
              vbar(i,j,knew)=cff2
            END DO
          END IF
        END DO
!
!  Apply boundary conditions.
!
        IF (.not.(ANY(LBC(:,isUbar,ng)%radiation).or.                   &
     &            ANY(LBC(:,isVbar,ng)%radiation).or.                   &
     &            ANY(LBC(:,isUbar,ng)%Flather).or.                     &
     &            ANY(LBC(:,isVbar,ng)%Flather))) THEN
          CALL u2dbc_tile (ng, tile,                                    &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
     &                     krhs, kstp, kstp,                            &
     &                     ubar, vbar, zeta)
          CALL v2dbc_tile (ng, tile,                                    &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
     &                     krhs, kstp, kstp,                            &
     &                     ubar, vbar, zeta)

          CALL u2dbc_tile (ng, tile,                                    &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
     &                     krhs, kstp, knew,                            &
     &                     ubar, vbar, zeta)
          CALL v2dbc_tile (ng, tile,                                    &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
     &                     krhs, kstp, knew,                            &
     &                     ubar, vbar, zeta)
        END IF
      END IF
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_u2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          ubar(:,:,kstp))
        CALL exchange_v2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          vbar(:,:,kstp))

        CALL exchange_u2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          ubar(:,:,knew))
        CALL exchange_v2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          vbar(:,:,knew))
      END IF

#  ifdef DISTRIBUTE
!
      CALL mp_exchange2d (ng, tile, model, 4,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    ubar(:,:,kstp), vbar(:,:,kstp),               &
     &                    ubar(:,:,knew), vbar(:,:,knew))
#  endif

# else
!
!-----------------------------------------------------------------------
!  If not perfect restart, initialize other time levels for 2D momentum.
!-----------------------------------------------------------------------
!
      IF (.not.PerfectRST(ng)) THEN
        DO j=JstrB,JendB
          DO i=IstrM,IendB
            cff1=ubar(i,j,kstp)
#  ifdef MASKING
            cff1=cff1*umask(i,j)
#  endif
#  ifdef WET_DRY
            cff1=cff1*umask_wet(i,j)
#  endif
            ubar(i,j,kstp)=cff1
            ubar(i,j,krhs)=cff1
          END DO
          IF (j.ge.JstrM) THEN
            DO i=IstrB,IendB
              cff2=vbar(i,j,kstp)
#  ifdef MASKING
              cff2=cff2*vmask(i,j)
#  endif
#  ifdef WET_DRY
              cff2=cff2*vmask_wet(i,j)
#  endif
              vbar(i,j,kstp)=cff2
              vbar(i,j,krhs)=cff2
            END DO
          END IF
        END DO
!
!  Apply boundary conditions.
!
        IF (.not.(ANY(LBC(:,isUbar,ng)%radiation).or.                   &
     &            ANY(LBC(:,isVbar,ng)%radiation).or.                   &
     &            ANY(LBC(:,isUbar,ng)%Flather).or.                     &
     &            ANY(LBC(:,isVbar,ng)%Flather))) THEN
          CALL u2dbc_tile (ng, tile,                                    &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
     &                     krhs, kstp, kstp,                            &
     &                     ubar, vbar, zeta)
          CALL v2dbc_tile (ng, tile,                                    &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
     &                     krhs, kstp, kstp,                            &
     &                     ubar, vbar, zeta)

          CALL u2dbc_tile (ng, tile,                                    &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
     &                     krhs, kstp, krhs,                            &
     &                     ubar, vbar, zeta)
          CALL v2dbc_tile (ng, tile,                                    &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
     &                     krhs, kstp, krhs,                            &
     &                     ubar, vbar, zeta)
        END IF
      END IF
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_u2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          ubar(:,:,kstp))
        CALL exchange_v2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          vbar(:,:,kstp))

        CALL exchange_u2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          ubar(:,:,krhs))
        CALL exchange_v2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          vbar(:,:,krhs))

        IF (PerfectRST(ng)) THEN
          CALL exchange_u2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            ubar(:,:,knew))
          CALL exchange_v2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            vbar(:,:,knew))
        END IF
      END IF

#  ifdef DISTRIBUTE
!
      CALL mp_exchange2d (ng, tile, model, 4,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    ubar(:,:,kstp), vbar(:,:,kstp),               &
     &                    ubar(:,:,krhs), vbar(:,:,krhs))

      IF (PerfectRST(ng)) THEN
        CALL mp_exchange2d (ng, tile, model, 2,                         &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      NghostPoints,                               &
     &                      EWperiodic(ng), NSperiodic(ng),             &
     &                      ubar(:,:,knew), vbar(:,:,knew))
      END IF
#  endif
# endif

# ifdef SOLVE3D
!
!-----------------------------------------------------------------------
!  If not perfect restart, initialize other time levels for tracers.
!-----------------------------------------------------------------------
!
      ic=0
      IF (.not.PerfectRST(ng)) THEN
        DO itrc=1,NT(ng)
          IF (LtracerCLM(itrc,ng).and.LnudgeTCLM(itrc,ng)) THEN
            ic=ic+1
          END IF
          DO k=1,N(ng)
            DO j=JstrB,JendB
              DO i=IstrB,IendB
                cff1=t(i,j,k,nstp,itrc)
#  ifdef MASKING
                cff1=cff1*rmask(i,j)
#  endif
                t(i,j,k,nstp,itrc)=cff1
                t(i,j,k,nnew,itrc)=cff1
              END DO
            END DO
          END DO
!
!  Apply boundary conditions.
!
          CALL t3dbc_tile (ng, tile, itrc, ic,                          &
     &                     LBi, UBi, LBj, UBj, N(ng), NT(ng),           &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
     &                     nstp, nstp,                                  &
     &                     t)
          CALL t3dbc_tile (ng, tile, itrc, ic,                          &
     &                     LBi, UBi, LBj, UBj, N(ng), NT(ng),           &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
     &                     nstp, nnew,                                  &
     &                     t)
!!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!#  ifdef T_PASSIVE
!         ELSE
!!         else condition added by Weifeng Zhang for passive tracer
!!         boundary condition
!          CALL pt3dbc_tile (ng, tile, itrc,                             &
!     &                      LBi, UBi, LBj, UBj, N(ng), NT(ng),          &
!     &                      IminS, ImaxS, JminS, JmaxS,                 &
!     &                      nstp, nstp,                                 &
!     &                      t)
!          CALL pt3dbc_tile (ng, tile, itrc,                             &
!     &                      LBi, UBi, LBj, UBj, N(ng), NT(ng),          &
!     &                      IminS, ImaxS, JminS, JmaxS,                 &
!     &                      nstp, nnew,                                 &
!     &                      t)
!#  endif
!!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!         END IF
        END DO
      END IF
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        DO itrc=1,NT(ng)
          CALL exchange_r3d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj, 1, N(ng),         &
     &                            t(:,:,:,nstp,itrc))
          CALL exchange_r3d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj, 1, N(ng),         &
     &                            t(:,:,:,nnew,itrc))
        END DO
      END IF

#  ifdef DISTRIBUTE
!
      CALL mp_exchange4d (ng, tile, model, 2,                           &
     &                    LBi, UBi, LBj, UBj, 1, N(ng), 1, NT(ng),      &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    t(:,:,:,nstp,:),                              &
     &                    t(:,:,:,nnew,:))
#  endif

#  if defined BBL_MODEL || defined SEDIMENT
!
!-----------------------------------------------------------------------
!  Initialize sediment bed properties.
!-----------------------------------------------------------------------
!
#   ifdef SEDIMENT
      IF (.not.PerfectRST(ng)) THEN
        DO k=1,Nbed
          DO j=JstrT,JendT
            DO i=IstrT,IendT
              DO ised=1,NST
                cff1=bed_mass(i,j,k,1,ised)
                bed_mass(i,j,k,2,ised)=cff1
              END DO
            END DO
          END DO
        END DO
      END IF
#   endif
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
#   ifdef SEDIMENT
        DO ised=1,NST
          CALL exchange_r3d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj, 1, Nbed,          &
     &                            bed_frac(:,:,:,ised))
          CALL exchange_r3d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj, 1, Nbed,          &
     &                            bed_mass(:,:,:,1,ised))
          CALL exchange_r3d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj, 1, Nbed,          &
     &                            bed_mass(:,:,:,2,ised))
        END DO
        DO ised=1,MBEDP
          CALL exchange_r3d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj, 1, Nbed,          &
     &                            bed(:,:,:,ised))
        END DO
#   endif
        CALL exchange_r3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, MBOTP,           &
     &                          bottom)
      END IF

#   ifdef DISTRIBUTE
#    ifdef SEDIMENT
      CALL mp_exchange4d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj, 1, Nbed, 1, NST,          &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    bed_frac)
      CALL mp_exchange4d (ng, tile, model, 2,                           &
     &                    LBi, UBi, LBj, UBj, 1, Nbed, 1, NST,          &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    bed_mass(:,:,:,1,:),bed_mass(:,:,:,2,:))
      CALL mp_exchange4d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj, 1, Nbed, 1, MBEDP,        &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    bed)
#    endif
      CALL mp_exchange3d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj, 1, MBOTP,                 &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    bottom)
#   endif
#  endif
# endif

# ifdef ICE_MODEL
      IF (.not.PerfectRST(ng)) THEN
        CALL ice_limit(ng, TILE)
        DO j=Jstr,Jend
          DO i=IstrU,Iend
            cff1=ui(i,j,nstp)
#  ifdef MASKING
            cff1=cff1*umask(i,j)
#  endif
            ui(i,j,nstp)=cff1
            ui(i,j,nnew)=cff1
          END DO
        END DO
        DO j=JstrV,Jend
          DO i=Istr,Iend
            cff1=vi(i,j,nstp)
#  ifdef MASKING
            cff1=cff1*vmask(i,j)
#  endif
            vi(i,j,nstp)=cff1
            vi(i,j,nnew)=cff1
          END DO
        END DO
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            cff1=ai(i,j,nstp)
            cff2=hi(i,j,nstp)
            cff3=hsn(i,j,nstp)
            cff4=ti(i,j,nstp)
#  ifdef MELT_PONDS
            cff5=apond(i,j,nstp)
            cff6=hpond(i,j,nstp)
#  endif
            cff7=ageice(i,j,nstp)
#  ifdef MASKING
            cff1=cff1*rmask(i,j)
            cff2=cff2*rmask(i,j)
            cff3=cff3*rmask(i,j)
            cff4=cff4*rmask(i,j)
#   ifdef MELT_PONDS
            cff5=cff5*rmask(i,j)
            cff6=cff6*rmask(i,j)
#   endif
            cff7=cff7*rmask(i,j)
#  endif
            ai(i,j,nstp)=cff1
            ai(i,j,nnew)=cff1
            hi(i,j,nstp)=cff2
            hi(i,j,nnew)=cff2
            hsn(i,j,nstp)=cff3
            hsn(i,j,nnew)=cff3
            ti(i,j,nstp)=cff4
            ti(i,j,nnew)=cff4
#  ifdef MELT_PONDS
            apond(i,j,nstp)=cff5
            apond(i,j,nnew)=cff5
            hpond(i,j,nstp)=cff6
            hpond(i,j,nnew)=cff6
#  endif
            ageice(i,j,nstp)=cff7
            ageice(i,j,nnew)=cff7
          END DO
        END DO
      END IF
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        DO i=1,2
          CALL exchange_u2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            ui(:,:,i))
          CALL exchange_v2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            vi(:,:,i))
          CALL exchange_r2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            ai(:,:,i))
          CALL exchange_r2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            hi(:,:,i))
          CALL exchange_r2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            hsn(:,:,i))
          CALL exchange_r2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            ti(:,:,i))
#  ifdef MELT_PONDS
          CALL exchange_r2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            apond(:,:,i))
          CALL exchange_r2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            hpond(:,:,i))
#  endif
          CALL exchange_r2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            ageice(:,:,i))
        END DO
      END IF
#  ifdef DISTRIBUTE
        CALL mp_exchange3d (ng, tile, model, 4,                         &
     &                      LBi, UBi, LBj, UBj, 1, 2,                   &
     &                      NghostPoints,                               &
     &                      EWperiodic(ng), NSperiodic(ng),             &
     &                      ui, vi, ai, hi)
        CALL mp_exchange3d (ng, tile, model, 3,                         &
     &                      LBi, UBi, LBj, UBj, 1, 2,                   &
     &                      NghostPoints,                               &
     &                      EWperiodic(ng), NSperiodic(ng),             &
     &                      hsn, ti, ageice)
#   ifdef MELT_PONDS
        CALL mp_exchange3d (ng, tile, model, 2,                         &
     &                      LBi, UBi, LBj, UBj, 1, 2,                   &
     &                      NghostPoints,                               &
     &                      EWperiodic(ng), NSperiodic(ng),             &
     &                      apond, hpond)
#   endif
#  endif
# endif

# ifdef PERFECT_RESTART
!
!-----------------------------------------------------------------------
!  If perfect restart, apply periodic boundary condition to
!  right-hand-side terms.
!-----------------------------------------------------------------------
!
      IF (PerfectRST(ng)) THEN
        IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
          DO i=1,2
            CALL exchange_u2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              rubar(:,:,i))
            CALL exchange_v2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              rvbar(:,:,i))
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              rzeta(:,:,i))
#  ifdef SOLVE3D
            CALL exchange_u3d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, 0, N(ng),       &
     &                              ru(:,:,:,i))
            CALL exchange_v3d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, 0, N(ng),       &
     &                              rv(:,:,:,i))
#  endif
          END DO
        END IF

#  ifdef DISTRIBUTE
        CALL mp_exchange3d (ng, tile, model, 3,                         &
     &                      LBi, UBi, LBj, UBj, 1, 2,                   &
     &                      NghostPoints,                               &
     &                      EWperiodic(ng), NSperiodic(ng),             &
     &                      rubar, rvbar, rzeta)
#   ifdef SOLVE3D
        CALL mp_exchange4d (ng, tile, model, 2,                         &
     &                      LBi, UBi, LBj, UBj, 0, N(ng), 1, 2,         &
     &                      NghostPoints,                               &
     &                      EWperiodic(ng), NSperiodic(ng),             &
     &                      ru, rv)
#   endif
#  endif
      END IF
# endif

      RETURN
      END SUBROUTINE ini_fields_tile
!
!***********************************************************************
      SUBROUTINE ini_zeta (ng, tile, model)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
# ifdef SOLVE3D
      USE mod_coupling
# endif
      USE mod_ocean
# if defined SOLVE3D && defined SEDIMENT && defined SED_MORPH
      USE mod_sedbed
# endif
      USE mod_stepping
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
!
!  Local variable declarations.
!
# include "tile.h"
!
# ifdef PROFILE
      CALL wclock_on (ng, iNLM, 2)
# endif
      CALL ini_zeta_tile (ng, tile, model,                              &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    kstp(ng), krhs(ng), knew(ng),                 &
# ifdef MASKING
     &                    GRID(ng) % rmask,                             &
# endif
# ifdef WET_DRY
     &                    GRID(ng) % rmask_wet,                         &
     &                    GRID(ng) % h,                                 &
# endif
# ifdef SOLVE3D
#  if defined SEDIMENT && defined SED_MORPH
     &                    SEDBED(ng) % bed,                             &
     &                    SEDBED(ng) % bed_thick0,                      &
     &                    SEDBED(ng) % bed_thick,                       &
#  endif
     &                    COUPLING(ng) % Zt_avg1,                       &
# endif
     &                    OCEAN(ng) % zeta)
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 2)
# endif

      RETURN
      END SUBROUTINE ini_zeta
!
!***********************************************************************
      SUBROUTINE ini_zeta_tile (ng, tile, model,                        &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          IminS, ImaxS, JminS, JmaxS,             &
     &                          kstp, krhs, knew,                       &
# ifdef MASKING
     &                          rmask,                                  &
# endif
# ifdef WET_DRY
     &                          rmask_wet,                              &
     &                          h,                                      &
# endif
# ifdef SOLVE3D
#  if defined SEDIMENT && defined SED_MORPH
     &                          bed, bed_thick0, bed_thick,             &
#  endif
     &                          Zt_avg1,                                &
# endif
     &                          zeta)
!***********************************************************************
!
      USE mod_param
      USE mod_ncparam
      USE mod_scalars
# if defined SEDIMENT && defined SED_MORPH
      USE mod_sediment
# endif
!
      USE exchange_2d_mod, ONLY : exchange_r2d_tile
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
# endif
      USE zetabc_mod, ONLY : zetabc_tile
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: kstp, krhs, knew
!
# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
#  endif
#  ifdef WET_DRY
      real(r8), intent(in) :: h(LBi:,LBj:)
#  endif
#  if defined SOLVE3D && defined SEDIMENT && defined SED_MORPH
      real(r8), intent(in) :: bed(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: bed_thick0(LBi:,LBj:)
      real(r8), intent(inout) :: bed_thick(LBi:,LBj:,:)
#  endif
#  ifdef WET_DRY
      real(r8), intent(inout) :: rmask_wet(LBi:,LBj:)
#  endif
#  ifdef SOLVE3D
      real(r8), intent(inout) :: Zt_avg1(LBi:,LBj:)
#  endif
      real(r8), intent(inout) :: zeta(LBi:,LBj:,:)

# else

#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
#  endif
#  ifdef WET_DRY
      real(r8), intent(in) :: h(LBi:UBi,LBj:UBj)
#  endif
#  if defined SOLVE3D && defined SEDIMENT && defined SED_MORPH
      real(r8), intent(in) :: bed(LBi:UBi,LBj:UBj,Nbed,MBEDP)
      real(r8), intent(inout) :: bed_thick0(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: bed_thick(LBi:UBi,LBj:UBj,2)
#  endif
#  ifdef WET_DRY
      real(r8), intent(inout) :: rmask_wet(LBi:UBi,LBj:UBj)
#  endif
#  ifdef SOLVE3D
      real(r8), intent(inout) :: Zt_avg1(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(inout) :: zeta(LBi:UBi,LBj:UBj,3)
# endif
!
!  Local variable declarations.
!
      integer :: Imin, Imax, Jmin, Jmax
      integer :: i, j, kbed

      real(r8) :: cff1

# include "set_bounds.h"

#ifndef OFFLINE_FLOATS
!
!-----------------------------------------------------------------------
!  Initialize other time levels for free-surface.
!-----------------------------------------------------------------------
!
      IF (.not.PerfectRST(ng)) THEN
        IF (.not.(ANY(LBC(:,isFsur,ng)%radiation).or.                   &
     &            ANY(LBC(:,isFsur,ng)%Chapman_explicit).or.            &
     &            ANY(LBC(:,isFsur,ng)%Chapman_implicit))) THEN
          Imin=IstrB
          Imax=IendB
          Jmin=JstrB
          Jmax=JendB
        ELSE
          Imin=IstrT
          Imax=IendT
          Jmin=JstrT
          Jmax=JendT
        END IF
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            cff1=zeta(i,j,kstp)
# ifdef MASKING
            cff1=cff1*rmask(i,j)
# endif
# ifdef WET_DRY
            IF (cff1.le.(Dcrit(ng)-h(i,j))) THEN
              cff1=Dcrit(ng)-h(i,j)
!!            rmask_wet(i,j)=0.0_r8
!!          ELSE
!!            rmask_wet(i,j)=1.0_r8*rmask(i,j)
            END IF
# endif
            zeta(i,j,kstp)=cff1
# ifdef SOLVE3D
            zeta(i,j,knew)=cff1
# else
            zeta(i,j,krhs)=cff1
# endif
          END DO
        END DO
!
!  Apply boundary conditions.
!
        IF (.not.(ANY(LBC(:,isFsur,ng)%radiation).or.                   &
     &            ANY(LBC(:,isFsur,ng)%Chapman_explicit).or.            &
     &            ANY(LBC(:,isFsur,ng)%Chapman_implicit))) THEN
          CALL zetabc_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      krhs, kstp, kstp,                           &
     &                      zeta)
# ifdef SOLVE3D
          CALL zetabc_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      krhs, kstp, knew,                           &
     &                      zeta)
# else
          CALL zetabc_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      krhs, kstp, krhs,                           &
     &                      zeta)
# endif
        END IF
      END IF
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          zeta(:,:,kstp))
# ifdef SOLVE3D
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          zeta(:,:,knew))
# else
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          zeta(:,:,krhs))
# endif

        IF (PerfectRST(ng)) THEN
# ifdef SOLVE3D
          CALL exchange_r2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            zeta(:,:,krhs))
# else
          CALL exchange_r2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            zeta(:,:,knew))
# endif
        END IF

# ifdef WET_DRY
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          rmask_wet)
# endif
      END IF

# ifdef DISTRIBUTE
#  ifdef SOLVE3D
      CALL mp_exchange2d (ng, tile, model, 2,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    zeta(:,:,kstp),                               &
     &                    zeta(:,:,knew))
#  else
      CALL mp_exchange2d (ng, tile, model, 2,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    zeta(:,:,kstp),                               &
     &                    zeta(:,:,krhs))
#  endif

      IF (PerfectRST(ng)) THEN
#  ifdef SOLVE3D
        CALL mp_exchange2d (ng, tile, model, 1,                         &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      NghostPoints,                               &
     &                      EWperiodic(ng), NSperiodic(ng),             &
     &                      zeta(:,:,krhs))
#  else
        CALL mp_exchange2d (ng, tile, model, 1,                         &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      NghostPoints,                               &
     &                      EWperiodic(ng), NSperiodic(ng),             &
     &                      zeta(:,:,knew))
#  endif
      END IF

#  ifdef WET_DRY
      CALL mp_exchange2d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    rmask_wet)
#  endif
# endif

#endif  /* OFFLINE_FLOATS */

# ifdef SOLVE3D
!
!-----------------------------------------------------------------------
!  Initialize fast-time averaged free-surface (Zt_avg1) with the inital
!  free-surface
!-----------------------------------------------------------------------
!
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          Zt_avg1(i,j)=zeta(i,j,kstp)
        END DO
      END DO
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Zt_avg1)
      END IF

#  ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    Zt_avg1)
#  endif

#  if defined SEDIMENT && defined SED_MORPH
!
!-----------------------------------------------------------------------
!  Compute initial total thickness for all sediment bed layers.
!-----------------------------------------------------------------------
!
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          bed_thick0(i,j)=0.0_r8
          DO kbed=1,Nbed
            bed_thick0(i,j)=bed_thick0(i,j)+bed(i,j,kbed,ithck)
          END DO
          bed_thick(i,j,1)=bed_thick0(i,j)
          bed_thick(i,j,2)=bed_thick0(i,j)
        END DO
      END DO
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          bed_thick0)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          bed_thick(:,:,1))
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          bed_thick(:,:,2))
      END IF

#   ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, model, 3,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    bed_thick0,                                   &
     &                    bed_thick(:,:,1), bed_thick(:,:,2))
#   endif
#  endif
# endif

      RETURN
      END SUBROUTINE ini_zeta_tile
#endif
      END MODULE ini_fields_mod
