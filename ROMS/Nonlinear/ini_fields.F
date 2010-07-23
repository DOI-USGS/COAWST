#include "cppdefs.h"
      MODULE ini_fields_mod
#ifdef NONLINEAR
!
!svn $Id: ini_fields.F 838 2008-11-17 04:22:18Z jcwarner $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2010 The ROMS/TOMS Group       Andrew M. Moore   !
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

      PRIVATE
      PUBLIC :: ini_fields
      PUBLIC :: ini_zeta

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
      CALL ini_fields_tile (ng, tile, model,                            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      kstp(ng), krhs(ng), knew(ng),               &
# ifdef SOLVE3D
     &                      nstp(ng), nnew(ng),                         &
# endif
     &                      GRID(ng) % h,                               &
# ifdef MASKING
     &                      GRID(ng) % rmask,                           &
     &                      GRID(ng) % umask,                           &
     &                      GRID(ng) % vmask,                           &
# endif
# if defined PERFECT_RESTART && \
    (defined EW_PERIODIC || defined NS_PERIODIC)
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
#   ifdef SED_MORPH
     &                      SEDBED(ng) % bed_thick0,                    &
     &                      SEDBED(ng) % bed_thick,                     &
#   endif
#  endif
     &                      GRID(ng) % Hz,                              &
#  ifdef ICESHELF
     &                      GRID(ng) % zice,                            &
#  endif
     &                      GRID(ng) % z_r,                             &
     &                      GRID(ng) % z_w,                             &
     &                      COUPLING(ng) % Zt_avg1,                     &
     &                      OCEAN(ng) % t,                              &
     &                      OCEAN(ng) % u,                              &
     &                      OCEAN(ng) % v,                              &
# endif
     &                      OCEAN(ng) % ubar,                           &
     &                      OCEAN(ng) % vbar,                           &
     &                      OCEAN(ng) % zeta)

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
     &                            h,                                    &
# ifdef MASKING
     &                            rmask, umask, vmask,                  &
# endif
# if defined PERFECT_RESTART && \
    (defined EW_PERIODIC || defined NS_PERIODIC)
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
#   ifdef SED_MORPH
     &                            bed_thick0, bed_thick,                &
#   endif
#  endif
     &                            Hz,                                   &
#  ifdef ICESHELF
     &                            zice,                                 &
#  endif
     &                            z_r, z_w,                             &
     &                            Zt_avg1,                              &
     &                            t, u, v,                              &
# endif
     &                            ubar, vbar, zeta)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
# ifdef COMPOSED_GRID
      USE mod_ncparam
      USE mod_nesting
# endif
# if defined SEDIMENT || defined BBL_MODEL
      USE mod_sediment
# endif
!
# if defined EW_PERIODIC || defined NS_PERIODIC
      USE exchange_2d_mod
#  ifdef SOLVE3D
      USE exchange_3d_mod
#  endif
# endif
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
#  ifdef SOLVE3D
      USE mp_exchange_mod, ONLY : mp_exchange3d, mp_exchange4d
#  else
#   if defined PERFECT_RESTART && \
      (defined EW_PERIODIC || defined NS_PERIODIC)
      USE mp_exchange_mod, ONLY : mp_exchange3d
#   endif
#  endif
# endif
# ifdef COMPOSED_GRID
      USE nesting_mod, ONLY : allocate_nesting2d, interp_nesting2d
#  ifdef SOLVE3D
      USE nesting_mod, ONLY : allocate_nesting3d, interp_nesting3d
      USE nesting_mod, ONLY : allocate_nesting4d, interp_nesting4d
#  endif
# endif
# ifdef SOLVE3D
      USE set_depth_mod, ONLY : set_depth_tile
      USE t3dbc_mod, ONLY : t3dbc_tile
      USE u3dbc_mod, ONLY : u3dbc_tile
      USE v3dbc_mod, ONLY : v3dbc_tile
# endif
# ifndef OBC_M2RADIATION
      USE u2dbc_mod, ONLY : u2dbc_tile
      USE v2dbc_mod, ONLY : v2dbc_tile
# endif
!
!  Imported variable declarations.
!
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
      real(r8), intent(in) :: zeta(LBi:,LBj:,:)
#  ifdef SOLVE3D
      real(r8), intent(in) :: Zt_avg1(LBi:,LBj:)
#   ifdef ICESHELF
      real(r8), intent(in) :: zice(LBi:,LBj:)
#   endif
#   if defined SEDIMENT || defined BBL_MODEL
      real(r8), intent(inout) :: bottom(LBi:,LBj:,:)
#   endif
#   if defined SEDIMENT
      real(r8), intent(inout) :: bed(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: bed_frac(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: bed_mass(LBi:,LBj:,:,:,:)
#   endif
#  endif
      real(r8), intent(inout) :: h(LBi:,LBj:)
#  if defined PERFECT_RESTART && \
     (defined EW_PERIODIC || defined NS_PERIODIC)
#   ifdef SOLVE3D
      real(r8), intent(inout) :: ru(LBi:,LBj:,0:,:)
      real(r8), intent(inout) :: rv(LBi:,LBj:,0:,:)
#   endif
      real(r8), intent(inout) :: rubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: rvbar(LBi:,LBj:,:)
      real(r8), intent(inout) :: rzeta(LBi:,LBj:,:)
#  endif
#  ifdef SOLVE3D
#   if defined SEDIMENT && defined SED_MORPH
      real(r8), intent(inout) :: bed_thick0(LBi:,LBj:)
      real(r8), intent(inout) :: bed_thick(LBi:,LBj:,:)
#   endif
      real(r8), intent(inout) :: t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: v(LBi:,LBj:,:,:)

      real(r8), intent(inout) :: Hz(LBi:,LBj:,:)
      real(r8), intent(inout) :: z_r(LBi:,LBj:,:)
      real(r8), intent(inout) :: z_w(LBi:,LBj:,0:)
#  endif
      real(r8), intent(inout) :: ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: vbar(LBi:,LBj:,:)

# else

#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: zeta(LBi:UBi,LBj:UBj,3)
#  ifdef SOLVE3D
      real(r8), intent(in) :: Zt_avg1(LBi:UBi,LBj:UBj)
#   ifdef ICESHELF
      real(r8), intent(in) :: zice(LBi:UBi,LBj:UBj)
#   endif
#   if defined SEDIMENT || defined BBL_MODEL
      real(r8), intent(inout) :: bottom(LBi:UBi,LBj:UBj,MBOTP)
#   endif
#   if defined SEDIMENT
      real(r8), intent(inout) :: bed(LBi:UBi,LBj:UBj,Nbed,MBEDP)
      real(r8), intent(inout) :: bed_frac(LBi:UBi,LBj:UBj,Nbed,NST)
      real(r8), intent(inout) :: bed_mass(LBi:UBi,LBj:UBj,Nbed,2,NST)
#   endif
#  endif
      real(r8), intent(inout) :: h(LBi:UBi,LBj:UBj)
#  if defined PERFECT_RESTART && \
     (defined EW_PERIODIC || defined NS_PERIODIC)
#   ifdef SOLVE3D
      real(r8), intent(inout) :: ru(LBi:UBi,LBj:UBj,0:N(ng),2)
      real(r8), intent(inout) :: rv(LBi:UBi,LBj:UBj,0:N(ng),2)
#   endif
      real(r8), intent(inout) :: rubar(LBi:UBi,LBj:UBj,2)
      real(r8), intent(inout) :: rvbar(LBi:UBi,LBj:UBj,2)
      real(r8), intent(inout) :: rzeta(LBi:UBi,LBj:UBj,2)
#  endif
#  ifdef SOLVE3D
#   if defined SEDIMENT && defined SED_MORPH
      real(r8), intent(inout) :: bed_thick0(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: bed_thick(LBi:UBi,LBj:UBj,2)
#   endif
      real(r8), intent(inout) :: t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: v(LBi:UBi,LBj:UBj,N(ng),2)

      real(r8), intent(inout) :: Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: z_w(LBi:UBi,LBj:UBj,0:N(ng))
#  endif
      real(r8), intent(inout) :: ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: vbar(LBi:UBi,LBj:UBj,3)
# endif
!
!  Local variable declarations.
!
      integer :: i, itrc, j, k, kbed
# if defined SEDIMENT || defined BBL_MODEL
      integer :: ised
# endif

      real(r8) :: cff1, cff2
# ifdef SOLVE3D
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: CF
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: DC
# endif
# ifdef COMPOSED_GRID
      TYPE (T_NEST2D), pointer :: field2d(:)
#  ifdef SOLVE3D
      TYPE (T_NEST3D), pointer :: field3d(:)
      TYPE (T_NEST4D), pointer :: field4d(:)
#  endif
# endif

# include "set_bounds.h"

# ifdef SOLVE3D
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

#   if defined EW_PERIODIC || defined NS_PERIODIC
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        bed_thick0)
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        bed_thick(:,:,1))
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        bed_thick(:,:,2))
#   endif
#   ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, model, 3,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    bed_thick0,                                   &
     &                    bed_thick(:,:,1), bed_thick(:,:,2))
#   endif
#  endif
# endif

# ifndef SOLVE3D
!
!-----------------------------------------------------------------------
!  Initialize other 2D momentum time-levels.
!-----------------------------------------------------------------------
!
      IF (.not.PerfectRST(ng)) THEN
#  ifdef REFINED_GRID
        DO j=JstrTV,JendT
          DO i=IstrTU+1,IendT
#  else
        DO j=Jstr,Jend
          DO i=IstrU,Iend
#  endif
            cff1=ubar(i,j,kstp)
#  ifdef MASKING
            cff1=cff1*umask(i,j)
#  endif
            ubar(i,j,kstp)=cff1
            ubar(i,j,krhs)=cff1
          END DO
#  ifdef REFINED_GRID
          IF (j.ge.JstrTV+1) THEN
            DO i=IstrTU,IendT
#  else
          IF (j.ge.JstrV) THEN
            DO i=Istr,Iend
#  endif
              cff2=vbar(i,j,kstp)
#  ifdef MASKING
              cff2=cff2*vmask(i,j)
#  endif
              vbar(i,j,kstp)=cff2
              vbar(i,j,krhs)=cff2
            END DO
          END IF
        END DO

#  ifndef OBC_M2RADIATION
!
!  Apply boundary conditions.
!
        CALL u2dbc_tile (ng, tile,                                      &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   IminS, ImaxS, JminS, JmaxS,                    &
     &                   krhs, kstp, kstp,                              &
     &                   ubar, vbar, zeta)
        CALL v2dbc_tile (ng, tile,                                      &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   IminS, ImaxS, JminS, JmaxS,                    &
     &                   krhs, kstp, kstp,                              &
     &                   ubar, vbar, zeta)
        CALL u2dbc_tile (ng, tile,                                      &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   IminS, ImaxS, JminS, JmaxS,                    &
     &                   krhs, kstp, krhs,                              &
     &                   ubar, vbar, zeta)
        CALL v2dbc_tile (ng, tile,                                      &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   IminS, ImaxS, JminS, JmaxS,                    &
     &                   krhs, kstp, krhs,                              &
     &                   ubar, vbar, zeta)
#  endif
      END IF

#  if defined EW_PERIODIC || defined NS_PERIODIC
!
      CALL exchange_u2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        ubar(:,:,kstp))
      CALL exchange_v2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        vbar(:,:,kstp))
      CALL exchange_u2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        ubar(:,:,krhs))
      CALL exchange_v2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        vbar(:,:,krhs))
      IF (PerfectRST(ng)) THEN
        CALL exchange_u2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          ubar(:,:,knew))
        CALL exchange_v2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          vbar(:,:,knew))
      END IF
#  endif
#  ifdef DISTRIBUTE
!
      CALL mp_exchange2d (ng, tile, model, 4,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    ubar(:,:,kstp), vbar(:,:,kstp),               &
     &                    ubar(:,:,krhs), vbar(:,:,krhs))
      IF (PerfectRST(ng)) THEN
        CALL mp_exchange2d (ng, tile, model, 2,                         &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      NghostPoints, EWperiodic, NSperiodic,       &
     &                      ubar(:,:,knew), vbar(:,:,knew))
      END IF
#  endif
#  ifdef COMPOSED_GRID
      IF ((ng.eq.Ngrids).and.(Ngrids.gt.1)) THEN
        DO i=1,3
          CALL allocate_nesting2d (field2d, idUbar, u2dvar, i)
          CALL interp_nesting2d (field2d)
          CALL allocate_nesting2d (field2d, idVbar, v2dvar, i)
          CALL interp_nesting2d (field2d)
        END DO
      END IF
#  endif

# else
!
!-----------------------------------------------------------------------
!  Initialize other 3D momentum time-levels.
!-----------------------------------------------------------------------
!
      IF (.not.PerfectRST(ng)) THEN
#  ifdef REFINED_GRID
        DO j=JstrT,JendT
          DO k=1,N(ng)
            DO i=IstrT,IendT
#  else
        DO j=Jstr,Jend
          DO k=1,N(ng)
            DO i=IstrU,Iend
#  endif
              cff1=u(i,j,k,nstp)
#  ifdef MASKING
              cff1=cff1*umask(i,j)
#  endif
              u(i,j,k,nstp)=cff1
              u(i,j,k,nnew)=cff1
            END DO
          END DO
#  ifdef REFINED_GRID
            DO k=1,N(ng)
              DO i=IstrT,IendT
#  else
          IF (j.ge.JstrV) THEN
            DO k=1,N(ng)
              DO i=Istr,Iend
#  endif
                cff2=v(i,j,k,nstp)
#  ifdef MASKING
                cff2=cff2*vmask(i,j)
#  endif
                v(i,j,k,nstp)=cff2
                v(i,j,k,nnew)=cff2
              END DO
            END DO
#  ifdef REFINED_GRID
#  else
          END IF
#  endif
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

#  if defined EW_PERIODIC || defined NS_PERIODIC
!
      CALL exchange_u3d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        u(:,:,:,nstp))
      CALL exchange_v3d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        v(:,:,:,nstp))

      CALL exchange_u3d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        u(:,:,:,nnew))
      CALL exchange_v3d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        v(:,:,:,nnew))
#  endif
#  ifdef DISTRIBUTE
!
      CALL mp_exchange3d (ng, tile, model, 4,                           &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    u(:,:,:,nstp), v(:,:,:,nstp),                 &
     &                    u(:,:,:,nnew), v(:,:,:,nnew))
#  endif
#  ifdef COMPOSED_GRID
      IF ((ng.eq.Ngrids).and.(Ngrids.gt.1)) THEN
        DO i=1,2
          CALL allocate_nesting3d (field3d, idUvel, u3dvar, i)
          CALL interp_nesting3d (field3d)
          CALL allocate_nesting3d (field3d, idVvel, v3dvar, i)
          CALL interp_nesting3d (field3d)
        END DO
      END IF
#  endif
!
!-----------------------------------------------------------------------
!  Initialize other tracers time-levels.
!-----------------------------------------------------------------------
!
      IF (.not.PerfectRST(ng)) THEN
        DO itrc=1,NT(ng)
          DO k=1,N(ng)
#  ifdef REFINED_GRID
            DO j=JstrT,JendT
              DO i=IstrT,IendT
#  else
            DO j=Jstr,Jend
              DO i=Istr,Iend
#  endif
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
          CALL t3dbc_tile (ng, tile, itrc,                              &
     &                     LBi, UBi, LBj, UBj, N(ng), NT(ng),           &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
     &                     nstp, nstp,                                  &
     &                     t)
          CALL t3dbc_tile (ng, tile, itrc,                              &
     &                     LBi, UBi, LBj, UBj, N(ng), NT(ng),           &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
     &                     nstp, nnew,                                  &
     &                     t)
        END DO
      END IF

#  if defined EW_PERIODIC || defined NS_PERIODIC
!
      DO itrc=1,NT(ng)
        CALL exchange_r3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          t(:,:,:,nstp,itrc))
        CALL exchange_r3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          t(:,:,:,nnew,itrc))
      END DO
#  endif
#  ifdef DISTRIBUTE
      CALL mp_exchange4d (ng, tile, model, 2,                           &
     &                    LBi, UBi, LBj, UBj, 1, N(ng), 1, NT(ng),      &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    t(:,:,:,nstp,:),                              &
     &                    t(:,:,:,nnew,:))
#  endif
#  ifdef COMPOSED_GRID
      IF ((ng.eq.Ngrids).and.(Ngrids.gt.1)) THEN
!       DO itrc=1,NT(ng)
          DO i=1,2
            CALL allocate_nesting4d (field4d, idTvar(itemp), r3dvar, i)
            CALL interp_nesting4d (field4d)
          END DO
!       END DO
      END IF
#  endif
!
!-----------------------------------------------------------------------
!  Compute initial depths and thicknesses.
!-----------------------------------------------------------------------
!
      CALL set_depth_tile (ng, tile,                                    &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
     &                     nstp, nnew,                                  &
     &                     h,                                           &
#  ifdef ICESHELF
     &                     zice,                                        &
#  endif
#  if defined SEDIMENT && defined SED_MORPH
     &                     bed_thick,                                   &
#  endif
     &                     Zt_avg1,                                     &
     &                     Hz, z_r, z_w)
!
!-----------------------------------------------------------------------
!  Compute vertically integrated momentum (ubar,vbar) from initial 3D
!  momentum (u,v).
!-----------------------------------------------------------------------
!
!  Compute adjoint 2D velocity component in the XI-direction.  Here
!  DC(i,1:N) are the thicknesses of U-boxes, DC(i,0) is total depth of
!  the water column, and CF(i,0) is the vertical integral.
!
      IF (.not.PerfectRST(ng)) THEN
#  ifdef REFINED_GRID
        DO j=JstrTV,JendT
          DO i=IstrTU+1,IendT
#  else
        DO j=Jstr,Jend
          DO i=IstrU,Iend
#  endif
            DC(i,0)=0.0_r8
            CF(i,0)=0.0_r8
          END DO
          DO k=1,N(ng)
#  ifdef REFINED_GRID
            DO i=IstrTU+1,IendT
#  else
            DO i=IstrU,Iend
#  endif
              DC(i,k)=0.5_r8*(Hz(i,j,k)+Hz(i-1,j,k))
              DC(i,0)=DC(i,0)+DC(i,k)
              CF(i,0)=CF(i,0)+DC(i,k)*u(i,j,k,nstp)
            END DO
          END DO
#  ifdef REFINED_GRID
          DO i=IstrTU+1,IendT
#  else
          DO i=IstrU,Iend
#  endif
            cff1=1.0_r8/DC(i,0)
            cff2=CF(i,0)*cff1
#  ifdef MASKING
            cff2=cff2*umask(i,j)
#  endif
            ubar(i,j,kstp)=cff2
            ubar(i,j,knew)=cff2
          END DO
!
!  Compute adjoint 2D velocity component in the ETA-direction.  Here
!  DC(i,1:N) are the thicknesses of V-boxes, DC(i,0) is total depth of
!  the water column, and CF(i,0) is the vertical integral.
!
#  ifdef REFINED_GRID
          IF (j.ge.JstrTV+1) THEN
            DO i=IstrTU,IendT
#  else
          IF (j.ge.Jstr) THEN
            DO i=Istr,Iend
#  endif
              DC(i,0)=0.0_r8
              CF(i,0)=0.0_r8
            END DO
            DO k=1,N(ng)
#  ifdef REFINED_GRID
              DO i=IstrTU,IendT
#  else
              DO i=Istr,Iend
#  endif
                DC(i,k)=0.5_r8*(Hz(i,j,k)+Hz(i,j-1,k))
                DC(i,0)=DC(i,0)+DC(i,k)
                CF(i,0)=CF(i,0)+DC(i,k)*v(i,j,k,nstp)
              END DO
            END DO
#  ifdef REFINED_GRID
            DO i=IstrTU,IendT
#  else
            DO i=Istr,Iend
#  endif
              cff1=1.0_r8/DC(i,0)
              cff2=CF(i,0)*cff1
#  ifdef MASKING
              cff2=cff2*vmask(i,j)
#  endif
              vbar(i,j,kstp)=cff2
              vbar(i,j,knew)=cff2
            END DO
          END IF
        END DO
#  ifndef OBC_M2RADIATION
!
!  Apply boundary conditions.
!
        CALL u2dbc_tile (ng, tile,                                      &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   IminS, ImaxS, JminS, JmaxS,                    &
     &                   krhs, kstp, kstp,                              &
     &                   ubar, vbar, zeta)
        CALL v2dbc_tile (ng, tile,                                      &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   IminS, ImaxS, JminS, JmaxS,                    &
     &                   krhs, kstp, kstp,                              &
     &                   ubar, vbar, zeta)
        CALL u2dbc_tile (ng, tile,                                      &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   IminS, ImaxS, JminS, JmaxS,                    &
     &                   krhs, kstp, knew,                              &
     &                   ubar, vbar, zeta)
        CALL v2dbc_tile (ng, tile,                                      &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   IminS, ImaxS, JminS, JmaxS,                    &
     &                   krhs, kstp, knew,                              &
     &                   ubar, vbar, zeta)
#  endif
      END IF

#  if defined EW_PERIODIC || defined NS_PERIODIC
!
      CALL exchange_u2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        ubar(:,:,kstp))
      CALL exchange_v2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        vbar(:,:,kstp))
      CALL exchange_u2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        ubar(:,:,knew))
      CALL exchange_v2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        vbar(:,:,knew))
#  endif
#  ifdef DISTRIBUTE
!
      CALL mp_exchange2d (ng, tile, model, 4,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    ubar(:,:,kstp), vbar(:,:,kstp),               &
     &                    ubar(:,:,knew), vbar(:,:,knew))
#  endif
#  ifdef COMPOSED_GRID
      IF ((ng.eq.Ngrids).and.(Ngrids.gt.1)) THEN
        DO i=1,3
          CALL allocate_nesting2d (field2d, idUbar, u2dvar, i)
          CALL interp_nesting2d (field2d)
          CALL allocate_nesting2d (field2d, idVbar, v2dvar, i)
          CALL interp_nesting2d (field2d)
        END DO
      END IF
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
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              DO ised=1,NST
                cff1=bed_mass(i,j,k,1,ised)
                bed_mass(i,j,k,2,ised)=cff1
              END DO
            END DO
          END DO
        END DO
      END IF
#   endif
#   if defined EW_PERIODIC || defined NS_PERIODIC
#    ifdef SEDIMENT
      DO ised=1,NST
        CALL exchange_r3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, Nbed,            &
     &                          bed_frac(:,:,:,ised))
        CALL exchange_r3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, Nbed,            &
     &                          bed_mass(:,:,:,1,ised))
        CALL exchange_r3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, Nbed,            &
     &                          bed_mass(:,:,:,2,ised))
      END DO
      DO ised=1,MBEDP
        CALL exchange_r3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, Nbed,            &
     &                          bed(:,:,:,ised))
      END DO
#    endif
      CALL exchange_r3d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj, 1, MBOTP,             &
     &                        bottom)
#   endif
#   ifdef DISTRIBUTE
#    ifdef SEDIMENT
      CALL mp_exchange4d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj, 1, Nbed, 1, NST,          &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    bed_frac)
      CALL mp_exchange4d (ng, tile, model, 2,                           &
     &                    LBi, UBi, LBj, UBj, 1, Nbed, 1, NST,          &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    bed_mass(:,:,:,1,:),bed_mass(:,:,:,2,:))
      CALL mp_exchange4d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj, 1, Nbed, 1, MBEDP,        &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    bed)
#    endif
      CALL mp_exchange3d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj, 1, MBOTP,                 &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    bottom)
#   endif
#  endif
# endif
# if defined PERFECT_RESTART && \
    (defined EW_PERIODIC || defined NS_PERIODIC)
!
!-----------------------------------------------------------------------
!  If perfect restart, apply periodic boundary condition to
!  right-hand-side terms.
!-----------------------------------------------------------------------
!
      IF (PerfectRST(ng)) THEN
        DO i=1,2
          CALL exchange_u2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            rubar(:,:,i))
          CALL exchange_v2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            rvbar(:,:,i))
          CALL exchange_r2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            rzeta(:,:,i))
#  ifdef SOLVE3D
          CALL exchange_u3d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj, 0, N(ng),         &
     &                            ru(:,:,:,i))
          CALL exchange_v3d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj, 0, N(ng),         &
     &                            rv(:,:,:,i))
#  endif
        END DO
#  ifdef DISTRIBUTE
        CALL mp_exchange3d (ng, tile, model, 3,                         &
     &                      LBi, UBi, LBj, UBj, 1, 2,                   &
     &                      NghostPoints, EWperiodic, NSperiodic,       &
     &                      rubar, rvbar, rzeta)
#   ifdef SOLVE3D
        CALL mp_exchange4d (ng, tile, model, 2,                         &
     &                      LBi, UBi, LBj, UBj, 0, N(ng), 1, 2,         &
     &                      NghostPoints, EWperiodic, NSperiodic,       &
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
     &                    COUPLING(ng) % Zt_avg1,                       &
# endif
     &                    OCEAN(ng) % zeta)

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
     &                          Zt_avg1,                                &
# endif
     &                          zeta)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
# ifdef COMPOSED_GRID
      USE mod_ncparam
      USE mod_nesting
# endif
!
# if defined EW_PERIODIC || defined NS_PERIODIC
      USE exchange_2d_mod, ONLY : exchange_r2d_tile
# endif
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
# endif
# ifndef OBC_FSRADIATION
      USE zetabc_mod, ONLY : zetabc_tile
# endif
# ifdef COMPOSED_GRID
      USE nesting_mod, ONLY : allocate_nesting2d, interp_nesting2d
#  ifdef SOLVE3D
      USE nesting_mod, ONLY : allocate_nesting3d, interp_nesting3d
      USE nesting_mod, ONLY : allocate_nesting4d, interp_nesting4d
#  endif
# endif
!
!  Imported variable declarations.
!
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
      integer :: i, j

      real(r8) :: cff1, cff2
# ifdef COMPOSED_GRID
      TYPE (T_NEST2D), pointer :: field2d(:)
# endif
# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Initialize other free-surface time-levels.
!-----------------------------------------------------------------------
!
      IF (.not.PerfectRST(ng)) THEN
# ifndef OBC_FSRADIATION
        DO j=Jstr,Jend
          DO i=Istr,Iend
# else
        DO j=JstrT,JendT
          DO i=IstrT,IendT
# endif
            cff1=zeta(i,j,kstp)
# ifdef MASKING
            cff1=cff1*rmask(i,j)
# endif
# ifdef WET_DRY
            IF (cff1.le.(Dcrit(ng)-h(i,j))) THEN
              cff1=Dcrit(ng)-h(i,j)
              rmask_wet(i,j)=0.0_r8
            ELSE
              rmask_wet(i,j)=1.0_r8*rmask(i,j)
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

# ifndef OBC_FSRADIATION
!
!  Apply boundary conditions.
!
        CALL zetabc_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    krhs, kstp, kstp,                             &
     &                    zeta)
#  ifdef SOLVE3D
        CALL zetabc_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    krhs, kstp, knew,                             &
     &                    zeta)
#  else
        CALL zetabc_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    krhs, kstp, krhs,                             &
     &                    zeta)
#  endif
# endif
      END IF

# if defined EW_PERIODIC || defined NS_PERIODIC
!
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        zeta(:,:,kstp))
#  ifdef SOLVE3D
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        zeta(:,:,knew))
#  else
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        zeta(:,:,krhs))
#  endif
      IF (PerfectRST(ng)) THEN
#  ifdef SOLVE3D
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          zeta(:,:,krhs))
#  else
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          zeta(:,:,knew))
#  endif
      END IF
#  ifdef WET_DRY
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          rmask_wet)
#  endif
# endif

# ifdef DISTRIBUTE
#  ifdef SOLVE3D
      CALL mp_exchange2d (ng, tile, model, 2,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    zeta(:,:,kstp),                               &
     &                    zeta(:,:,knew))
#  else
      CALL mp_exchange2d (ng, tile, model, 2,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    zeta(:,:,kstp),                               &
     &                    zeta(:,:,krhs))
#  endif
      IF (PerfectRST(ng)) THEN
#  ifdef SOLVE3D
        CALL mp_exchange2d (ng, tile, model, 1,                         &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      NghostPoints, EWperiodic, NSperiodic,       &
     &                      zeta(:,:,krhs))
#  else
        CALL mp_exchange2d (ng, tile, model, 1,                         &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      NghostPoints, EWperiodic, NSperiodic,       &
     &                      zeta(:,:,knew))
#  endif
      END IF
#  ifdef WET_DRY
        CALL mp_exchange2d (ng, tile, model, 1,                         &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      NghostPoints, EWperiodic, NSperiodic,       &
     &                      rmask_wet)
#  endif
# endif
# ifdef COMPOSED_GRID
      IF ((ng.eq.Ngrids).and.(Ngrids.gt.1)) THEN
        DO i=1,3
          CALL allocate_nesting2d (field2d, idFsur, r2dvar, i)
          CALL interp_nesting2d (field2d)
        END DO
      END IF
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
          Zt_avg1(i,j)=zeta(i,j,kstp)
        END DO
      END DO

#  if defined EW_PERIODIC || defined NS_PERIODIC
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Zt_avg1)
#  endif
#  ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    Zt_avg1)
#  endif
# endif
# ifdef COMPOSED_GRID
      IF ((ng.eq.Ngrids).and.(Ngrids.gt.1)) THEN
        CALL allocate_nesting2d (field2d, idZavg, r2dvar)
        CALL interp_nesting2d (field2d)
      END IF
# endif

      RETURN
      END SUBROUTINE ini_zeta_tile
#endif
      END MODULE ini_fields_mod
