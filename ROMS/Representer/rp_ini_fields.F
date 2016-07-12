#include "cppdefs.h"
      MODULE rp_ini_fields_mod
#ifdef TL_IOMS
!
!svn $Id: rp_ini_fields.F 795 2016-05-11 01:42:43Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine initializes other time levels for tangent 2D fields.   !
!  It also couples 3D and 2D momentum equations:  it initializes 2D    !
!  momentum (ubar,vbar) to  the  vertical  integral  of  initial 3D    !
!  momentum (u,v).                                                     !
!                                                                      !
!=======================================================================
!
      implicit none
!
      PRIVATE
      PUBLIC :: rp_ini_fields
      PUBLIC :: rp_ini_zeta
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE rp_ini_fields (ng, tile, model)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
# ifdef SOLVE3D
      USE mod_coupling
# endif
      USE mod_ocean
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
      CALL wclock_on (ng, iRPM, 2)
# endif
      CALL rp_ini_fields_tile (ng, tile, model,                         &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         IminS, ImaxS, JminS, JmaxS,              &
     &                         kstp(ng), krhs(ng), knew(ng),            &
# ifdef SOLVE3D
     &                         nstp(ng), nnew(ng),                      &
# endif
# ifdef MASKING
     &                         GRID(ng) % rmask,                        &
     &                         GRID(ng) % umask,                        &
     &                         GRID(ng) % vmask,                        &
# endif
# ifdef SOLVE3D
     &                         GRID(ng) % Hz,                           &
     &                         GRID(ng) % tl_Hz,                        &
     &                         OCEAN(ng) % tl_t,                        &
     &                         OCEAN(ng) % u,                           &
     &                         OCEAN(ng) % tl_u,                        &
     &                         OCEAN(ng) % v,                           &
     &                         OCEAN(ng) % tl_v,                        &
# endif
     &                         OCEAN(ng) % ubar,                        &
     &                         OCEAN(ng) % tl_ubar,                     &
     &                         OCEAN(ng) % vbar,                        &
     &                         OCEAN(ng) % tl_vbar,                     &
     &                         OCEAN(ng) % zeta,                        &
     &                         OCEAN(ng) % tl_zeta)
# ifdef PROFILE
      CALL wclock_off (ng, iRPM, 2)
# endif

      RETURN
      END SUBROUTINE rp_ini_fields
!
!***********************************************************************
      SUBROUTINE rp_ini_fields_tile (ng, tile, model,                   &
     &                               LBi, UBi, LBj, UBj,                &
     &                               IminS, ImaxS, JminS, JmaxS,        &
     &                               kstp, krhs, knew,                  &
# ifdef SOLVE3D
     &                               nstp, nnew,                        &
# endif
# ifdef MASKING
     &                               rmask, umask, vmask,               &
# endif
# ifdef SOLVE3D
     &                               Hz, tl_Hz,                         &
     &                               tl_t,                              &
     &                               u, tl_u,                           &
     &                               v, tl_v,                           &
# endif
     &                               ubar, tl_ubar,                     &
     &                               vbar, tl_vbar,                     &
     &                               zeta, tl_zeta)
!***********************************************************************
!
      USE mod_param
      USE mod_ncparam
      USE mod_scalars
!
      USE exchange_2d_mod
# ifdef SOLVE3D
      USE exchange_3d_mod
# endif
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
#  ifdef SOLVE3D
      USE mp_exchange_mod, ONLY : mp_exchange3d, mp_exchange4d
#  endif
# endif
# ifdef SOLVE3D
      USE rp_t3dbc_mod, ONLY : rp_t3dbc_tile
      USE rp_u3dbc_mod, ONLY : rp_u3dbc_tile
      USE rp_v3dbc_mod, ONLY : rp_v3dbc_tile
# endif
      USE rp_u2dbc_mod, ONLY : rp_u2dbc_tile
      USE rp_v2dbc_mod, ONLY : rp_v2dbc_tile
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
#  ifdef SOLVE3D
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: tl_Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: v(LBi:,LBj:,:,:)
#  endif
      real(r8), intent(in) :: ubar(LBi:,LBj:,:)
      real(r8), intent(in) :: vbar(LBi:,LBj:,:)
      real(r8), intent(in) :: zeta(LBi:,LBj:,:)
      real(r8), intent(in) :: tl_zeta(LBi:,LBj:,:)
#  ifdef SOLVE3D
      real(r8), intent(inout) :: tl_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: tl_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: tl_v(LBi:,LBj:,:,:)
#  endif
      real(r8), intent(inout) :: tl_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: tl_vbar(LBi:,LBj:,:)

# else

#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#  endif
#  ifdef SOLVE3D
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: tl_Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: v(LBi:UBi,LBj:UBj,N(ng),2)
#  endif
      real(r8), intent(in) :: ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: vbar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: zeta(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: tl_zeta(LBi:UBi,LBj:UBj,3)
#  ifdef SOLVE3D
      real(r8), intent(inout) :: tl_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: tl_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: tl_v(LBi:UBi,LBj:UBj,N(ng),2)
#  endif
      real(r8), intent(inout) :: tl_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: tl_vbar(LBi:UBi,LBj:UBj,3)
# endif
!
!  Local variable declarations.
!
      integer :: i, ic, itrc, j, k

      real(r8) :: cff1
      real(r8) :: tl_cff1, tl_cff2

# ifdef SOLVE3D
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: CF
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: DC

      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: tl_CF
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: tl_DC
# endif

# include "set_bounds.h"

# ifdef SOLVE3D
!
!-----------------------------------------------------------------------
!  Initialize other time levels for 3D momentum.
!-----------------------------------------------------------------------
!
      DO j=JstrB,JendB
        DO k=1,N(ng)
          DO i=IstrM,IendB
!>          cff1=u(i,j,k,nstp)
!>
            tl_cff1=tl_u(i,j,k,nstp)
#  ifdef MASKING
!>          cff1=cff1*umask(i,j)
!>
            tl_cff1=tl_cff1*umask(i,j)
#  endif
!>          u(i,j,k,nstp)=cff1
!>          u(i,j,k,nnew)=cff1
!>
            tl_u(i,j,k,nstp)=tl_cff1
            tl_u(i,j,k,nnew)=tl_cff1
          END DO
        END DO
!
        IF (j.ge.JstrM) THEN
          DO k=1,N(ng)
            DO i=IstrB,IendB
!>            cff2=v(i,j,k,nstp)
!>
              tl_cff2=tl_v(i,j,k,nstp)
#  ifdef MASKING
!>            cff2=cff2*vmask(i,j)
!>
              tl_cff2=tl_cff2*vmask(i,j)
#  endif
!>            v(i,j,k,nstp)=cff2
!>            v(i,j,k,nnew)=cff2
!>
              tl_v(i,j,k,nstp)=tl_cff2
              tl_v(i,j,k,nnew)=tl_cff2
            END DO
          END DO
        END IF
      END DO
!
!  Apply boundary conditions.
!
!>    CALL u3dbc_tile (ng, tile,                                        &
!>   &                 LBi, UBi, LBj, UBj, N(ng),                       &
!>   &                 IminS, ImaxS, JminS, JmaxS,                      &
!>   &                 nstp, nstp,                                      &
!>   &                 u)
!>
      CALL rp_u3dbc_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj, N(ng),                    &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    nstp, nstp,                                   &
     &                    tl_u)
!>    CALL v3dbc_tile (ng, tile,                                        &
!>   &                 LBi, UBi, LBj, UBj, N(ng),                       &
!>   &                 IminS, ImaxS, JminS, JmaxS,                      &
!>   &                 nstp, nstp,                                      &
!>   &                 v)
!>
      CALL rp_v3dbc_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj, N(ng),                    &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    nstp, nstp,                                   &
     &                    tl_v)

!>    CALL u3dbc_tile (ng, tile,                                        &
!>   &                 LBi, UBi, LBj, UBj, N(ng),                       &
!>   &                 IminS, ImaxS, JminS, JmaxS,                      &
!>   &                 nstp, nnew,                                      &
!>   &                 u)
!>
      CALL rp_u3dbc_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj, N(ng),                    &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    nstp, nnew,                                   &
     &                    tl_u)
!>    CALL v3dbc_tile (ng, tile,                                        &
!>   &                 LBi, UBi, LBj, UBj, N(ng),                       &
!>   &                 IminS, ImaxS, JminS, JmaxS,                      &
!>   &                 nstp, nnew,                                      &
!>   &                 v)
!>
      CALL rp_v3dbc_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj, N(ng),                    &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    nstp, nnew,                                   &
     &                    tl_v)
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
!>      CALL exchange_u3d_tile (ng, tile,                               &
!>   &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
!>   &                          u(:,:,:,nstp))
!>
        CALL exchange_u3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          tl_u(:,:,:,nstp))
!>      CALL exchange_v3d_tile (ng, tile,                               &
!>   &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
!>   &                          v(:,:,:,nstp))
!>
        CALL exchange_v3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          tl_v(:,:,:,nstp))

!>      CALL exchange_u3d_tile (ng, tile,                               &
!>   &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
!>   &                          u(:,:,:,nnew))
!>
        CALL exchange_u3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          tl_u(:,:,:,nnew))
!>      CALL exchange_v3d_tile (ng, tile,                               &
!>   &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
!>   &                          v(:,:,:,nnew))
!>
        CALL exchange_v3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          tl_v(:,:,:,nnew))
      END IF

#  ifdef DISTRIBUTE
!
!>    CALL mp_exchange3d (ng, tile, model, 4,                           &
!>   &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
!>   &                    NghostPoints,                                 &
!>   &                    EWperiodic(ng), NSperiodic(ng),               &
!>   &                    u(:,:,:,nstp), v(:,:,:,nstp),                 &
!>   &                    u(:,:,:,nnew), v(:,:,:,nnew))
!>
      CALL mp_exchange3d (ng, tile, model, 4,                           &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    tl_u(:,:,:,nstp), tl_v(:,:,:,nstp),           &
     &                    tl_u(:,:,:,nnew), tl_v(:,:,:,nnew))
#  endif
# endif

# ifdef SOLVE3D
!
!-----------------------------------------------------------------------
!  Compute vertically-integrated momentum (tl_ubar, tl_vbar) from
!  initial 3D momentum (tl_u, tl_v).
!-----------------------------------------------------------------------
!
!  Here DC(i,1:N) are the grid cell thicknesses, DC(i,0) is the total
!  depth of the water column, and CF(i,0) is the vertical integral.
!
      DO j=JstrB,JendB
        DO i=IstrM,IendB
          DC(i,0)=0.0_r8
          tl_DC(i,0)=0.0_r8
          CF(i,0)=0.0_r8
          tl_CF(i,0)=0.0_r8
        END DO
        DO k=1,N(ng)
          DO i=IstrM,IendB
            DC(i,k)=0.5_r8*(Hz(i,j,k)+Hz(i-1,j,k))
            tl_DC(i,k)=0.5_r8*(tl_Hz(i,j,k)+tl_Hz(i-1,j,k))
            DC(i,0)=DC(i,0)+DC(i,k)
            tl_DC(i,0)=tl_DC(i,0)+tl_DC(i,k)
            CF(i,0)=CF(i,0)+DC(i,k)*u(i,j,k,nstp)
            tl_CF(i,0)=tl_CF(i,0)+tl_DC(i,k)*u(i,j,k,nstp)+             &
     &                            DC(i,k)*tl_u(i,j,k,nstp)-             &
#  ifdef TL_IOMS
     &                            DC(i,k)*u(i,j,k,nstp)
#  endif
          END DO
        END DO
        DO i=IstrM,IendB
          cff1=1.0_r8/DC(i,0)
          tl_cff1=-cff1*cff1*tl_DC(i,0)+                                &
#  ifdef TL_IOMS
     &            2.0_r8*cff1
#  endif
!>        cff2=CF(i,0)*cff1
!>
          tl_cff2=tl_CF(i,0)*cff1+CF(i,0)*tl_cff1-                      &
#  ifdef TL_IOMS
     &            CF(i,0)*cff1
#  endif
#  ifdef MASKING
!>        cff2=cff2*umask(i,j)
!>
          tl_cff2=tl_cff2*umask(i,j)
#  endif
!>        ubar(i,j,kstp)=cff2
!>        ubar(i,j,knew)=cff2
!>
          tl_ubar(i,j,kstp)=tl_cff2
          tl_ubar(i,j,knew)=tl_cff2
        END DO
!
        IF (j.ge.JstrM) THEN
          DO i=IstrB,IendB
            DC(i,0)=0.0_r8
            tl_DC(i,0)=0.0_r8
            CF(i,0)=0.0_r8
            tl_CF(i,0)=0.0_r8
          END DO
          DO k=1,N(ng)
            DO i=IstrB,IendB
              DC(i,k)=0.5_r8*(Hz(i,j,k)+Hz(i,j-1,k))
              tl_DC(i,k)=0.5_r8*(tl_Hz(i,j,k)+tl_Hz(i,j-1,k))
              DC(i,0)=DC(i,0)+DC(i,k)
              tl_DC(i,0)=tl_DC(i,0)+tl_DC(i,k)
              CF(i,0)=CF(i,0)+DC(i,k)*v(i,j,k,nstp)
              tl_CF(i,0)=tl_CF(i,0)+tl_DC(i,k)*v(i,j,k,nstp)+           &
     &                              DC(i,k)*tl_v(i,j,k,nstp)-           &
#  ifdef TL_IOMS
     &                              DC(i,k)*v(i,j,k,nstp)
#  endif
            END DO
          END DO
          DO i=IstrB,IendB
            cff1=1.0_r8/DC(i,0)
            tl_cff1=-cff1*cff1*tl_DC(i,0)+                              &
#  ifdef TL_IOMS
     &              2.0_r8*cff1
#  endif
!>          cff2=CF(i,0)*cff1
!>
            tl_cff2=tl_CF(i,0)*cff1+CF(i,0)*tl_cff1-                    &
#  ifdef TL_IOMS
     &              CF(i,0)*cff1
#  endif
#  ifdef MASKING
!>          cff2=cff2*vmask(i,j)
!>
            tl_cff2=tl_cff2*vmask(i,j)
#  endif
!>          vbar(i,j,kstp)=cff2
!>          vbar(i,j,knew)=cff2
!>
            tl_vbar(i,j,kstp)=tl_cff2
            tl_vbar(i,j,knew)=tl_cff2
          END DO
        END IF
      END DO
!
!  Apply boundary conditions.
!
      IF (.not.(ANY(tl_LBC(:,isUbar,ng)%radiation).or.                  &
     &          ANY(tl_LBC(:,isVbar,ng)%radiation).or.                  &
     &          ANY(tl_LBC(:,isUbar,ng)%Flather).or.                    &
     &          ANY(tl_LBC(:,isVbar,ng)%Flather))) THEN
!>      CALL u2dbc_tile (ng, tile,                                      &
!>   &                   LBi, UBi, LBj, UBj,                            &
!>   &                   IminS, ImaxS, JminS, JmaxS,                    &
!>   &                   krhs, kstp, kstp,                              &
!>   &                   ubar, vbar, zeta)
!>
        CALL rp_u2dbc_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      krhs, kstp, kstp,                           &
     &                      ubar, vbar, zeta,                           &
     &                      tl_ubar, tl_vbar, tl_zeta)
!>      CALL v2dbc_tile (ng, tile,                                      &
!>   &                   LBi, UBi, LBj, UBj,                            &
!>   &                   IminS, ImaxS, JminS, JmaxS,                    &
!>   &                   krhs, kstp, kstp,                              &
!>   &                   ubar, vbar, zeta)
!>
        CALL rp_v2dbc_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      krhs, kstp, kstp,                           &
     &                      ubar, vbar, zeta,                           &
     &                      tl_ubar, tl_vbar, tl_zeta)

!>      CALL u2dbc_tile (ng, tile,                                      &
!>   &                   LBi, UBi, LBj, UBj,                            &
!>   &                   IminS, ImaxS, JminS, JmaxS,                    &
!>   &                   krhs, kstp, knew,                              &
!>   &                   ubar, vbar, zeta)
!>
        CALL rp_u2dbc_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      krhs, kstp, knew,                           &
     &                      ubar, vbar, zeta,                           &
     &                      tl_ubar, tl_vbar, tl_zeta)
!>      CALL v2dbc_tile (ng, tile,                                      &
!>   &                   LBi, UBi, LBj, UBj,                            &
!>   &                   IminS, ImaxS, JminS, JmaxS,                    &
!>   &                   krhs, kstp, knew,                              &
!>   &                   ubar, vbar, zeta)
!>
        CALL rp_v2dbc_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      krhs, kstp, knew,                           &
     &                      ubar, vbar, zeta,                           &
     &                      tl_ubar, tl_vbar, tl_zeta)
      END IF
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
!>      CALL exchange_u2d_tile (ng, tile,                               &
!>   &                          LBi, UBi, LBj, UBj,                     &
!>   &                          ubar(:,:,kstp))
!>
        CALL exchange_u2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          tl_ubar(:,:,kstp))
!>      CALL exchange_v2d_tile (ng, tile,                               &
!>   &                          LBi, UBi, LBj, UBj,                     &
!>   &                          vbar(:,:,kstp))
!>
        CALL exchange_v2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          tl_vbar(:,:,kstp))

!>      CALL exchange_u2d_tile (ng, tile,                               &
!>   &                          LBi, UBi, LBj, UBj,                     &
!>   &                          ubar(:,:,knew))
!>
        CALL exchange_u2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          tl_ubar(:,:,knew))
!>      CALL exchange_v2d_tile (ng, tile,                               &
!>   &                          LBi, UBi, LBj, UBj,                     &
!>   &                          vbar(:,:,knew))
!>
        CALL exchange_v2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          tl_vbar(:,:,knew))
      END IF

#  ifdef DISTRIBUTE
!
!>    CALL mp_exchange2d (ng, tile, model, 4,                           &
!>   &                    LBi, UBi, LBj, UBj,                           &
!>   &                    NghostPoints,                                 &
!>   &                    EWperiodic(ng), NSperiodic(ng),               &
!>   &                    ubar(:,:,kstp), vbar(:,:,kstp),               &
!>   &                    ubar(:,:,knew), vbar(:,:,knew))
!>
      CALL mp_exchange2d (ng, tile, model, 4,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    tl_ubar(:,:,kstp), tl_vbar(:,:,kstp),         &
     &                    tl_ubar(:,:,knew), tl_vbar(:,:,knew))
#  endif

# else
!
!-----------------------------------------------------------------------
!  Initialize other time levels for 2D momentum.
!-----------------------------------------------------------------------
!
      DO j=JstrB,JendB
        DO i=IstrM,IendB
!>        cff1=ubar(i,j,kstp)
!>
          tl_cff1=tl_ubar(i,j,kstp)
#  ifdef MASKING
!>        cff1=cff1*umask(i,j)
!>
          tl_cff1=tl_cff1*umask(i,j)
#  endif
!>        ubar(i,j,kstp)=cff1
!>        ubar(i,j,krhs)=cff1
!>
          tl_ubar(i,j,kstp)=tl_cff1
          tl_ubar(i,j,krhs)=tl_cff1
        END DO
!
        IF (j.ge.JstrM) THEN
          DO i=IstrB,IendB
!>          cff2=vbar(i,j,kstp)
!>
            tl_cff2=tl_vbar(i,j,kstp)
#  ifdef MASKING
!>          cff2=cff2*vmask(i,j)
!>
            tl_cff2=tl_cff2*vmask(i,j)
#  endif
!>          vbar(i,j,kstp)=cff2
!>          vbar(i,j,krhs)=cff2
!>
            tl_vbar(i,j,kstp)=tl_cff2
            tl_vbar(i,j,krhs)=tl_cff2
          END DO
        END IF
      END DO
!
!  Apply boundary conditions.
!
      IF (.not.(ANY(tl_LBC(:,isUbar,ng)%radiation).or.                  &
     &          ANY(tl_LBC(:,isVbar,ng)%radiation).or.                  &
     &          ANY(tl_LBC(:,isUbar,ng)%Flather).or.                    &
     &          ANY(tl_LBC(:,isVbar,ng)%Flather))) THEN
!>      CALL u2dbc_tile (ng, tile,                                      &
!>   &                   LBi, UBi, LBj, UBj,                            &
!>   &                   IminS, ImaxS, JminS, JmaxS,                    &
!>   &                   krhs, kstp, kstp,                              &
!>   &                   ubar, vbar, zeta)
!>
        CALL rp_u2dbc_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      krhs, kstp, kstp,                           &
     &                      ubar, vbar, zeta,                           &
     &                      tl_ubar, tl_vbar, tl_zeta)
!>      CALL v2dbc_tile (ng, tile,                                      &
!>   &                   LBi, UBi, LBj, UBj,                            &
!>   &                   IminS, ImaxS, JminS, JmaxS,                    &
!>   &                   krhs, kstp, kstp,                              &
!>   &                   ubar, vbar, zeta)
!>
        CALL rp_v2dbc_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      krhs, kstp, kstp,                           &
     &                      ubar, vbar, zeta,                           &
     &                      tl_ubar, tl_vbar, tl_zeta)

!>      CALL u2dbc_tile (ng, tile,                                      &
!>   &                   LBi, UBi, LBj, UBj,                            &
!>   &                   IminS, ImaxS, JminS, JmaxS,                    &
!>   &                   krhs, kstp, krhs,                              &
!>   &                   ubar, vbar, zeta)
!>
        CALL rp_u2dbc_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      krhs, kstp, krhs,                           &
     &                      ubar, vbar, zeta,                           &
     &                      tl_ubar, tl_vbar, tl_zeta)
!>      CALL v2dbc_tile (ng, tile,                                      &
!>   &                   LBi, UBi, LBj, UBj,                            &
!>   &                   IminS, ImaxS, JminS, JmaxS,                    &
!>   &                   krhs, kstp, krhs,                              &
!>   &                   ubar, vbar, zeta)
!>
        CALL rp_v2dbc_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      krhs, kstp, krhs,                           &
     &                      ubar, vbar, zeta,                           &
     &                      tl_ubar, tl_vbar, tl_zeta)
      END IF
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
!>      CALL exchange_u2d_tile (ng, tile,                               &
!>   &                          LBi, UBi, LBj, UBj,                     &
!>   &                          ubar(:,:,kstp))
!>
        CALL exchange_u2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          tl_ubar(:,:,kstp))
!>      CALL exchange_v2d_tile (ng, tile,                               &
!>   &                          LBi, UBi, LBj, UBj,                     &
!>   &                          vbar(:,:,kstp))
!>
        CALL exchange_v2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          tl_vbar(:,:,kstp))

!>      CALL exchange_u2d_tile (ng, tile,                               &
!>   &                          LBi, UBi, LBj, UBj,                     &
!>   &                          ubar(:,:,krhs))
!>
        CALL exchange_u2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          tl_ubar(:,:,krhs))
!>      CALL exchange_v2d_tile (ng, tile,                               &
!>   &                          LBi, UBi, LBj, UBj,                     &
!>   &                          vbar(:,:,krhs))
!>
        CALL exchange_v2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          tl_vbar(:,:,krhs))
      END IF

#  ifdef DISTRIBUTE
!
!>    CALL mp_exchange2d (ng, tile, model, 4,                           &
!>   &                    LBi, UBi, LBj, UBj,                           &
!>   &                    NghostPoints,                                 &
!>   &                    EWperiodic(ng), NSperiodic(ng),               &
!>   &                    ubar(:,:,kstp), vbar(:,:,kstp),               &
!>   &                    ubar(:,:,krhs), vbar(:,:,krhs))
!>
      CALL mp_exchange2d (ng, tile, model, 4,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    tl_ubar(:,:,kstp), tl_vbar(:,:,kstp),         &
     &                    tl_ubar(:,:,krhs), tl_vbar(:,:,krhs))
#  endif
# endif

# ifdef SOLVE3D
!
!-----------------------------------------------------------------------
!  Initialize other time levels for tracers.
!-----------------------------------------------------------------------
!
      ic=0
      DO itrc=1,NT(ng)
        IF (LtracerCLM(itrc,ng).and.LnudgeTCLM(itrc,ng)) THEN
          ic=ic+1
        END IF
        DO k=1,N(ng)
          DO j=JstrB,JendB
            DO i=IstrB,IendB
!>            cff1=t(i,j,k,nstp,itrc)
!>
              tl_cff1=tl_t(i,j,k,nstp,itrc)
#  ifdef MASKING
              tl_cff1=tl_cff1*rmask(i,j)
#  endif
!>            t(i,j,k,nstp,itrc)=cff1
!>            t(i,j,k,nnew,itrc)=cff1
!>
              tl_t(i,j,k,nstp,itrc)=tl_cff1
              tl_t(i,j,k,nnew,itrc)=tl_cff1
            END DO
          END DO
        END DO
!
!  Apply boundary conditions.
!
!>      CALL t3dbc_tile (ng, tile, itrc, ic,                            &
!>   &                   LBi, UBi, LBj, UBj, N(ng), NT(ng),             &
!>   &                   IminS, ImaxS, JminS, JmaxS,                    &
!>   &                   nstp, nstp,                                    &
!>   &                   t)
!>
        CALL rp_t3dbc_tile (ng, tile, itrc, ic,                         &
     &                      LBi, UBi, LBj, UBj, N(ng), NT(ng),          &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      nstp, nstp,                                 &
     &                      tl_t)
!>      CALL t3dbc_tile (ng, tile, itrc, ic,                            &
!>   &                   LBi, UBi, LBj, UBj, N(ng), NT(ng),             &
!>   &                   IminS, ImaxS, JminS, JmaxS,                    &
!>   &                   nstp, nnew,                                    &
!>   &                   t)
!>
        CALL rp_t3dbc_tile (ng, tile, itrc, ic,                         &
     &                      LBi, UBi, LBj, UBj, N(ng), NT(ng),          &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      nstp, nnew,                                 &
     &                      tl_t)
!
        IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
!>        CALL exchange_r3d_tile (ng, tile,                             &
!>   &                            LBi, UBi, LBj, UBj, 1, N(ng),         &
!>   &                            t(:,:,:,nstp,itrc))
!>
          CALL exchange_r3d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj, 1, N(ng),         &
     &                            tl_t(:,:,:,nstp,itrc))
!>        CALL exchange_r3d_tile (ng, tile,                             &
!>   &                            LBi, UBi, LBj, UBj, 1, N(ng),         &
!>   &                            t(:,:,:,nnew,itrc))
!>
          CALL exchange_r3d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj, 1, N(ng),         &
     &                            tl_t(:,:,:,nnew,itrc))
        END IF
      END DO

#  ifdef DISTRIBUTE
!
!>    CALL mp_exchange4d (ng, tile, model, 2,                           &
!>   &                    LBi, UBi, LBj, UBj, 1, N(ng), 1, NT(ng),      &
!>   &                    NghostPoints,                                 &
!>   &                    EWperiodic(ng), NSperiodic(ng),               &
!>   &                    t(:,:,:,nstp,:),                              &
!>   &                    t(:,:,:,nnew,:))
!>
      CALL mp_exchange4d (ng, tile, model, 2,                           &
     &                    LBi, UBi, LBj, UBj, 1, N(ng), 1, NT(ng),      &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    tl_t(:,:,:,nstp,:),                           &
     &                    tl_t(:,:,:,nnew,:))
#  endif
# endif

      RETURN
      END SUBROUTINE rp_ini_fields_tile

!
!***********************************************************************
      SUBROUTINE rp_ini_zeta (ng, tile, model)
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
      CALL wclock_on (ng, iRPM, 2)
# endif
      CALL rp_ini_zeta_tile (ng, tile, model,                           &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       IminS, ImaxS, JminS, JmaxS,                &
     &                       kstp(ng), krhs(ng), knew(ng),              &
# ifdef MASKING
     &                       GRID(ng) % rmask,                          &
# endif
# ifdef WET_DRY_NOT_YET
     &                       GRID(ng) % h,                              &
# endif
# ifdef SOLVE3D
#  if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
     &                       SEDBED(ng) % tl_bed,                       &
     &                       SEDBED(ng) % tl_bed_thick0,                &
     &                       SEDBED(ng) % tl_bed_thick,                 &
#  endif
     &                       COUPLING(ng) % tl_Zt_avg1,                 &
# endif
     &                       OCEAN(ng) % zeta,                          &
     &                       OCEAN(ng) % tl_zeta)
# ifdef PROFILE
      CALL wclock_off (ng, iRPM, 2)
# endif

      RETURN
      END SUBROUTINE rp_ini_zeta
!
!***********************************************************************
      SUBROUTINE rp_ini_zeta_tile (ng, tile, model,                     &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             IminS, ImaxS, JminS, JmaxS,          &
     &                             kstp, krhs, knew,                    &
# ifdef MASKING
     &                             rmask,                               &
# endif
# ifdef WET_DRY_NOT_YET
     &                             h,                                   &
# endif
# ifdef SOLVE3D
#  if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
     &                             tl_bed,                              &
     &                             tl_bed_thick0, tl_bed_thick,         &
#  endif
     &                             tl_Zt_avg1,                          &
# endif
     &                             zeta, tl_zeta)
!***********************************************************************
!
      USE mod_param
      USE mod_ncparam
      USE mod_scalars
# if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
      USE mod_sediment
# endif
!
      USE exchange_2d_mod, ONLY : exchange_r2d_tile
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
# endif
      USE rp_zetabc_mod, ONLY : rp_zetabc_tile
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
#  ifdef WET_DRY_NOT_YET
      real(r8), intent(in) :: h(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: zeta(LBi:,LBj:,:)
#  ifdef SOLVE3D
#   if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
      real(r8), intent(in) :: tl_bed(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: tl_bed_thick0(LBi:,LBj:)
      real(r8), intent(inout) :: tl_bed_thick(LBi:,LBj:,:)
#   endif
      real(r8), intent(inout) :: tl_Zt_avg1(LBi:,LBj:)
#  endif
      real(r8), intent(inout) :: tl_zeta(LBi:,LBj:,:)

# else

#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
#  endif
#  ifdef WET_DRY_NOT_YET
      real(r8), intent(in) :: h(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: zeta(LBi:UBi,LBj:UBj,3)
#  ifdef SOLVE3D
#   if defined SOLVE3D && defined SEDIMENT && defined SED_MORPH
      real(r8), intent(in) :: tl_bed(LBi:UBi,LBj:UBj,Nbed,MBEDP)
      real(r8), intent(inout) :: tl_bed_thick0(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: tl_bed_thick(LBi:UBi,LBj:UBj,2)
#   endif
      real(r8), intent(inout) :: tl_Zt_avg1(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(inout) :: tl_zeta(LBi:UBi,LBj:UBj,3)
# endif
!
!  Local variable declarations.
!
      integer :: Imin, Imax, Jmin, Jmax
      integer :: i, j, kbed

      real(r8) :: cff1
      real(r8) :: tl_cff1

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Initialize other time levels for free-surface.
!-----------------------------------------------------------------------
!
      IF (.not.(ANY(tl_LBC(:,isFsur,ng)%radiation).or.                  &
     &          ANY(tl_LBC(:,isFsur,ng)%Chapman_explicit).or.           &
     &          ANY(tl_LBC(:,isFsur,ng)%Chapman_implicit))) THEN
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
!>        cff1=zeta(i,j,kstp)
!>
          tl_cff1=tl_zeta(i,j,kstp)
# ifdef MASKING
!>        cff1=cff1*rmask(i,j)
!>
          tl_cff1=tl_cff1*rmask(i,j)
# endif
!>        zeta(i,j,kstp)=cff1
!>
          tl_zeta(i,j,kstp)=tl_cff1
# ifdef SOLVE3D
!>        zeta(i,j,knew)=cff1
!>
          tl_zeta(i,j,knew)=tl_cff1
# else
!>        zeta(i,j,krhs)=cff1
!>
          tl_zeta(i,j,krhs)=tl_cff1
# endif
        END DO
      END DO
!
!  Apply boundary conditions.
!
      IF (.not.(ANY(tl_LBC(:,isFsur,ng)%radiation).or.                  &
     &          ANY(tl_LBC(:,isFsur,ng)%Chapman_explicit).or.           &
     &          ANY(tl_LBC(:,isFsur,ng)%Chapman_implicit))) THEN
!>      CALL zetabc_tile (ng, tile,                                     &
!>   &                    LBi, UBi, LBj, UBj,                           &
!>   &                    IminS, ImaxS, JminS, JmaxS,                   &
!>   &                    krhs, kstp, kstp,                             &
!>   &                    zeta)
!>
        CALL rp_zetabc_tile (ng, tile,                                  &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       IminS, ImaxS, JminS, JmaxS,                &
     &                       krhs, kstp, kstp,                          &
     &                       zeta,                                      &
     &                       tl_zeta)
# ifdef SOLVE3D
!>      CALL zetabc_tile (ng, tile,                                     &
!>   &                    LBi, UBi, LBj, UBj,                           &
!>   &                    IminS, ImaxS, JminS, JmaxS,                   &
!>   &                    krhs, kstp, knew,                             &
!>   &                    zeta)
!>
        CALL rp_zetabc_tile (ng, tile,                                  &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       IminS, ImaxS, JminS, JmaxS,                &
     &                       krhs, kstp, knew,                          &
     &                       zeta,                                      &
     &                       tl_zeta)
# else
!>      CALL zetabc_tile (ng, tile,                                     &
!>   &                    LBi, UBi, LBj, UBj,                           &
!>   &                    IminS, ImaxS, JminS, JmaxS,                   &
!>   &                    krhs, kstp, krhs,                             &
!>   &                    zeta)
!>
        CALL rp_zetabc_tile (ng, tile,                                  &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       IminS, ImaxS, JminS, JmaxS,                &
     &                       krhs, kstp, krhs,                          &
     &                       zeta,                                      &
     &                       tl_zeta)
# endif
      END IF
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
!>      CALL exchange_r2d_tile (ng, tile,                               &
!>   &                          LBi, UBi, LBj, UBj,                     &
!>   &                          zeta(:,:,kstp))
!>
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          tl_zeta(:,:,kstp))
# ifdef SOLVE3D
!>      CALL exchange_r2d_tile (ng, tile,                               &
!>   &                          LBi, UBi, LBj, UBj,                     &
!>   &                          zeta(:,:,knew))
!>
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          tl_zeta(:,:,knew))
# else
!>      CALL exchange_r2d_tile (ng, tile,                               &
!>   &                          LBi, UBi, LBj, UBj,                     &
!>   &                          zeta(:,:,krhs))
!>
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          tl_zeta(:,:,krhs))
# endif
      END IF

# ifdef DISTRIBUTE
#  ifdef SOLVE3D
!>    CALL mp_exchange2d (ng, tile, model, 2,                           &
!>   &                    LBi, UBi, LBj, UBj,                           &
!>   &                    NghostPoints,                                 &
!>   &                    EWperiodic(ng), NSperiodic(ng),               &
!>   &                    zeta(:,:,kstp),                               &
!>   &                    zeta(:,:,knew))
!>
      CALL mp_exchange2d (ng, tile, model, 2,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    tl_zeta(:,:,kstp),                            &
     &                    tl_zeta(:,:,knew))
#  else
!>    CALL mp_exchange2d (ng, tile, model, 2,                           &
!>   &                    LBi, UBi, LBj, UBj,                           &
!>   &                    NghostPoints,                                 &
!>   &                    EWperiodic(ng), NSperiodic(ng),               &
!>   &                    zeta(:,:,kstp),                               &
!>   &                    zeta(:,:,krhs))
!>
      CALL mp_exchange2d (ng, tile, model, 2,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    tl_zeta(:,:,kstp),                            &
     &                    tl_zeta(:,:,krhs))
#  endif
# endif

# ifdef SOLVE3D
!
!-----------------------------------------------------------------------
!  Initialize fast-time averaged free-surface (Zt_avg1) with the inital
!  free-surface
!-----------------------------------------------------------------------
!
      DO j=JstrT,JendT
        DO i=IstrT,IendT
!>        Zt_avg1(i,j)=zeta(i,j,kstp)
!>
          tl_Zt_avg1(i,j)=tl_zeta(i,j,kstp)
        END DO
      END DO
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
!>      CALL exchange_r2d_tile (ng, tile,                               &
!>   &                          LBi, UBi, LBj, UBj,                     &
!>   &                          Zt_avg1)
!>
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          tl_Zt_avg1)
      END IF

#  ifdef DISTRIBUTE
!>    CALL mp_exchange2d (ng, tile, model, 1,                           &
!>   &                    LBi, UBi, LBj, UBj,                           &
!>   &                    NghostPoints,                                 &
!>   &                    EWperiodic(ng), NSperiodic(ng),               &
!>   &                    Zt_avg1)
!>
      CALL mp_exchange2d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    tl_Zt_avg1)
#  endif

#  if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
!
!-----------------------------------------------------------------------
!  Compute initial total thickness for all sediment bed layers.
!-----------------------------------------------------------------------
!
      DO j=JstrT,JendT
        DO i=IstrT,IendT
!>        bed_thick0(i,j)=0.0_r8
!>
          tl_bed_thick0(i,j)=0.0_r8
          DO kbed=1,Nbed
!>          bed_thick0(i,j)=bed_thick0(i,j)+bed(i,j,kbed,ithck)
!>
            tl_bed_thick0(i,j)=tl_bed_thick0(i,j)+tl_bed(i,j,kbed,ithck)
          END DO
!>        bed_thick(i,j,1)=bed_thick0(i,j)
!>        bed_thick(i,j,2)=bed_thick0(i,j)
!>
          tl_bed_thick(i,j,1)=tl_bed_thick0(i,j)
          tl_bed_thick(i,j,2)=tl_bed_thick0(i,j)
        END DO
      END DO
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
!>      CALL exchange_r2d_tile (ng, tile,                               &
!>   &                          LBi, UBi, LBj, UBj,                     &
!>   &                          bed_thick0)
!>
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          tl_bed_thick0)
!>      CALL exchange_r2d_tile (ng, tile,                               &
!>   &                          LBi, UBi, LBj, UBj,                     &
!>   &                          bed_thick(:,:,1))
!>
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          tl_bed_thick(:,:,1))
!>      CALL exchange_r2d_tile (ng, tile,                               &
!>   &                          LBi, UBi, LBj, UBj,                     &
!>   &                          bed_thick(:,:,2))
!>
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          tl_bed_thick(:,:,2))
      END IF

#   ifdef DISTRIBUTE
!>    CALL mp_exchange2d (ng, tile, model, 3,                           &
!>   &                    LBi, UBi, LBj, UBj, 1, 1,                     &
!>   &                    NghostPoints,                                 &
!>   &                    EWperiodic(ng), NSperiodic(ng),               &
!>   &                    bed_thick0,                                   &
!>   &                    bed_thick(:,:,1), bed_thick(:,:,2))
!>
      CALL mp_exchange2d (ng, tile, model, 3,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    tl_bed_thick0,                                &
     &                    tl_bed_thick(:,:,1), tl_bed_thick(:,:,2))
#   endif
#  endif
# endif

      RETURN
      END SUBROUTINE rp_ini_zeta_tile
#endif
      END MODULE rp_ini_fields_mod
