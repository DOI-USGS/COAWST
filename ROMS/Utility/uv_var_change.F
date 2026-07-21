#include "cppdefs.h"
      MODULE uv_var_change_mod
#ifdef SOLVE3D
!
!git $Id$
!=======================================================================
!  Copyright (c) 2002-2026 The ROMS Group                              !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.md                            Hernan G. Arango   !
!=======================================================================
!                                                                      !
!  These routines are used for ocean current variable changes from     !
!  C-grid to A-grid and vice versa. It is done for output purposes     !
!  and data assimilation where the state vector is located at the      !
!  cell-center (Arakawa A-grid).                                       !
!                                                                      !
!  * If transforming vector components from C-grid to A-grid, rotate   !
!    to geographical Eastward and Northward directions.                !
!                                                                      !
!  Ur(i,j,k) = 0.5 * [u(i,j,k,ninp) + u(i+1,j,k,ninp)]   i=Istr:Iend   !
!  Vr(i,j,k) = 0.5 * [v(i,j,k,ninp) + v(i,j+1,k,ninp)]   j=Jstr:Jend   !
!                                                                      !
!  Apply lateral boundary conditions (gradient) via 'bc_r3d_tile'      !
!                                                                      !
!  ua(i,j,k) = Ur(i,j,k) * CosAngler(i,j) - Vr(i,j,k) * SinAngler(i,j) !
!  va(i,j,k) = Vr(i,j,k) * CosAngler(i,j) + Ur(i,j,k) * SinAngler(i,j) !
!                                                                      !
!  PUBLIC:  uv_C2A_grid, ad_uv_C2A_grid, tl_uv_C2A_grid                !
!  PRIVATE: uv_C2A_grid_tile, ad_uv_C2A_grid_tile, tl_uv_C2A_grid_tile !
!                                                                      !
!  * If transforming vector components from A-grid to C-grid, rotate   !
!    to computational XI and ETA directions.                           !
!                                                                      !
!  Ur(i,j,k) = ua(i,j,k) * CosAngler(i,j) + va(i,j,k) * SinAngler(i,j) !
!  Vr(i,j,k) = va(i,j,k) * CosAngler(i,j) - ua(i,j,k) * SinAngler(i,j) !
!                                                                      !
!  u(i,j,k,nout) = 0.5 * [Ur(i-1,j,k) + Ur(i,j,k)]      i=Istr:IendR   !
!  v(i,j,k<nout) = 0.5 * [Vr(i,j-1,k) + Vr(i,j,k)]      j=Jstr:JendR   !
!                                                                      !
!  PUBLIC:  uv_A2C_grid, ad_uv_A2C_grid, tl_uv_A2C_grid                !
!  PRIVATE: uv_A2C_grid_tile, ad_uv_A2C_grid_tile, tl_uv_A2C_grid_tile !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_grid
      USE mod_ocean
      USE mod_scalars
!
# ifdef DISTRIBUTE
#  ifdef ADJOINT
      USE mp_exchange_mod, ONLY : ad_mp_exchange3d
#  endif
      USE mp_exchange_mod, ONLY : mp_exchange3d
# endif
!
      implicit none
!
      PRIVATE
      PUBLIC  :: uv_A2C_grid
      PUBLIC  :: uv_C2A_grid
# ifdef ADJOINT
      PUBLIC  :: ad_uv_A2C_grid
      PUBLIC  :: ad_uv_C2A_grid
# endif
# if defined TANGENT || defined TL_IOMS
      PUBLIC  :: tl_uv_A2C_grid
      PUBLIC  :: tl_uv_C2A_grid
# endif
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE uv_C2A_grid (ng, tile, model, ninp)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model, ninp
!
!  Local variable declarations.
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__
!
# include "tile.h"
!
# ifdef PROFILE
      CALL wclock_on (ng, model, 34, __LINE__, MyFile)
# endif
      CALL uv_C2A_grid_tile (ng, tile, model, ninp,                     &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       IminS, ImaxS, JminS, JmaxS,                &
# ifdef MASKING
     &                       GRID(ng) % rmask_full,                     &
# endif
     &                       GRID(ng) % CosAngler,                      &
     &                       GRID(ng) % SinAngler,                      &
     &                       OCEAN(ng) % u,                             &
     &                       OCEAN(ng) % v,                             &
     &                       OCEAN(ng) % ua,                            &
     &                       OCEAN(ng) % va)
# ifdef PROFILE
      CALL wclock_off (ng, model, 34, __LINE__, MyFile)
# endif
!
      RETURN
      END SUBROUTINE uv_C2A_grid
!
!***********************************************************************
      SUBROUTINE uv_C2A_grid_tile (ng, tile, model, ninp,               &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             IminS, ImaxS, JminS, JmaxS,          &
# ifdef MASKING
     &                             rmask,                               &
# endif
     &                             CosAngler, SinAngler,                &
     &                             u, v, ua, va)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model, ninp
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: CosAngler(LBi:,LBj:)
      real(r8), intent(in) :: SinAngler(LBi:,LBj:)
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: v(LBi:,LBj:,:,:)

      real(r8), intent(inout) :: ua(LBi:,LBj:,:)
      real(r8), intent(inout) :: va(LBi:,LBj:,:)
# else
      real(r8), intent(in) :: CosAngler(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: SinAngler(LBi:UBi,LBj:UBj)
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: v(LBi:UBi,LBj:UBj,N(ng),2)

      real(r8), intent(inout) :: ua(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: va(LBi:UBi,LBj:UBj,N(ng))
# endif
!
!  Local variable declarations.
!
      integer :: i, j, k

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Urho, Vrho

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Tranform vector components from C-grid to A-grid.
!-----------------------------------------------------------------------
!
      K_LOOP : DO k=1,N(ng)
!
!  Compute A-grid (cell center) components. Apply gradient condition.
!
        DO j=JstrR,JendR
          DO i=Istr,Iend
            Urho(i,j)=0.5_r8*(u(i,j,k,ninp)+u(i+1,j,k,ninp))
            IF (.not.EWperiodic(ng)) THEN
              IF (DOMAIN(ng)%Western_Edge(tile)) THEN
                Urho(Istr-1,j)=Urho(Istr,j)
              END IF
              IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
                Urho(Iend+1,j)=Urho(Iend,j)
              END IF
            END IF
          END DO
        END DO
!
        DO j=Jstr,Jend
          DO i=IstrR,IendR
            Vrho(i,j)=0.5_r8*(v(i,j,k,ninp)+v(i,j+1,k,ninp))
            IF (.not.NSperiodic(ng)) THEN
              IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
                Vrho(i,Jstr-1) = Vrho(i,Jstr)
              END IF
              IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
                Vrho(i,Jend+1) = Vrho(i,Jend)
              END IF
            END IF
          END DO
        END DO
!
!  Rotate from computational to gegraphical Eastward and Northward
!  directions.
!
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            ua(i,j,k)=Urho(i,j)*CosAngler(i,j)-                         &
     &                Vrho(i,j)*SinAngler(i,j)
            va(i,j,k)=Vrho(i,j)*CosAngler(i,j)+                         &
     &                Urho(i,j)*SinAngler(i,j)
# ifdef MASKING
            ua(i,j,k)=ua(i,j,k)*rmask(i,j)
            va(i,j,k)=va(i,j,k)*rmask(i,j)
# endif
          END DO
        END DO
      END DO K_LOOP

# ifdef DISTRIBUTE
!
      CALL mp_exchange3d (ng, tile, model, 2,                           &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    ua, va)
# endif
!
      RETURN
      END SUBROUTINE uv_C2A_grid_tile    
!
!***********************************************************************
      SUBROUTINE uv_A2C_grid (ng, tile, model, nout)
!***********************************************************************
!
      USE mod_stepping
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model, nout
!
!  Local variable declarations.
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__
!
# include "tile.h"
!
# ifdef PROFILE
      CALL wclock_on (ng, model, 34, __LINE__, MyFile)
# endif
      CALL uv_A2C_grid_tile (ng, tile, model, nout,                     &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       IminS, ImaxS, JminS, JmaxS,                &
# ifdef MASKING
     &                       GRID(ng) % umask_full,                     &
     &                       GRID(ng) % vmask_full,                     &
# endif
     &                       GRID(ng) % CosAngler,                      &
     &                       GRID(ng) % SinAngler,                      &
     &                       OCEAN(ng) % ua,                            &
     &                       OCEAN(ng) % va,                            &
     &                       OCEAN(ng) % u,                             &
     &                       OCEAN(ng) % v)
# ifdef PROFILE
      CALL wclock_off (ng, model, 34, __LINE__, MyFile)
# endif
!
      RETURN
      END SUBROUTINE uv_A2C_grid
!
!***********************************************************************
      SUBROUTINE uv_A2C_grid_tile (ng, tile, model, nout,               &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             IminS, ImaxS, JminS, JmaxS,          &
# ifdef MASKING
     &                             umask, vmask,                        &
# endif
     &                             CosAngler, SinAngler,                &
     &                             ua, va, u, v)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model, nout
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: CosAngler(LBi:,LBj:)
      real(r8), intent(in) :: SinAngler(LBi:,LBj:)
#  ifdef MASKING
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: ua(LBi:,LBj:,:)
      real(r8), intent(in) :: va(LBi:,LBj:,:)

      real(r8), intent(inout) :: u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: v(LBi:,LBj:,:,:)
# else
      real(r8), intent(in) :: CosAngler(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: SinAngler(LBi:UBi,LBj:UBj)
#  ifdef MASKING
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: ua(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: va(LBi:UBi,LBj:UBj,N(ng))

      real(r8), intent(inout) :: u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: v(LBi:UBi,LBj:UBj,N(ng),2)
# endif
!
!  Local variable declarations.
!
      integer :: i, j, k
!
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Urho, Vrho

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Transform vector components from A-grid to C-grid.
!-----------------------------------------------------------------------
!
      K_LOOP : DO k=1,N(ng)
!
!  Rotate vector components to computational (XI,ETA) directions.
!
        DO j=Jstr-1,JendR
          DO i=Istr-1,IendR
            Urho(i,j)=ua(i,j,k)*CosAngler(i,j)+                         &
     &                va(i,j,k)*SinAngler(i,j)
            Vrho(i,j)=va(i,j,k)*CosAngler(i,j)-                         &
     &                ua(i,j,k)*SinAngler(i,j)
          END DO
        END DO
!
!  Compute staggered C-grid components.
!
        DO j=JstrR,JendR
          DO i=Istr,IendR
            u(i,j,k,nout)=0.5_r8*(Urho(i-1,j)+Urho(i,j))
# ifdef MASKING
            u(i,j,k,nout)=u(i,j,k,nout)*umask(i,j)
# endif
          END DO
        END DO
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            v(i,j,k,nout)=0.5_r8*(Vrho(i,j-1)+Vrho(i,j))
# ifdef MASKING
            v(i,j,k,nout)=v(i,j,k,nout)*vmask(i,j)
# endif
          END DO
        END DO
      END DO K_LOOP

# ifdef DISTRIBUTE
!
!  Exchange boundary data.
!
      CALL mp_exchange3d (ng, tile, model, 2,                           &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    u(:,:,:,nout), v(:,:,:,nout))
# endif
!
      RETURN
      END SUBROUTINE uv_A2C_grid_tile    

# ifdef TANGENT
!
!***********************************************************************
      SUBROUTINE tl_uv_C2A_grid (ng, tile, model, ninp)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model, ninp
!
!  Local variable declarations.
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__
!
#  include "tile.h"
!
#  ifdef PROFILE
      CALL wclock_on (ng, model, 34, __LINE__, MyFile)
#  endif
      CALL tl_uv_C2A_grid_tile (ng, tile, model, ninp,                  &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          IminS, ImaxS, JminS, JmaxS,             &
#  ifdef MASKING
     &                          GRID(ng) % rmask_full,                  &
#  endif
     &                          GRID(ng) % CosAngler,                   &
     &                          GRID(ng) % SinAngler,                   &
     &                          OCEAN(ng) % tl_u,                       &
     &                          OCEAN(ng) % tl_v,                       &
     &                          OCEAN(ng) % tl_ua,                      &
     &                          OCEAN(ng) % tl_va)
#  ifdef PROFILE
      CALL wclock_off (ng, model, 34, __LINE__, MyFile)
#  endif
!
      RETURN
      END SUBROUTINE tl_uv_C2A_grid
!
!***********************************************************************
      SUBROUTINE tl_uv_C2A_grid_tile (ng, tile, model, ninp,            &
     &                                LBi, UBi, LBj, UBj,               &
     &                                IminS, ImaxS, JminS, JmaxS,       &
#  ifdef MASKING
     &                                rmask,                            &
#  endif
     &                                CosAngler, SinAngler,             &
     &                                tl_u, tl_v, tl_ua, tl_va)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model, ninp
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
#  ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: CosAngler(LBi:,LBj:)
      real(r8), intent(in) :: SinAngler(LBi:,LBj:)
#   ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
#   endif
      real(r8), intent(in) :: tl_u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: tl_v(LBi:,LBj:,:,:)

      real(r8), intent(inout) :: tl_ua(LBi:,LBj:,:)
      real(r8), intent(inout) :: tl_va(LBi:,LBj:,:)
#  else
      real(r8), intent(in) :: CosAngler(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: SinAngler(LBi:UBi,LBj:UBj)
#   ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
#   endif
      real(r8), intent(in) :: tl_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: tl_v(LBi:UBi,LBj:UBj,N(ng),2)

      real(r8), intent(inout) :: tl_ua(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: tl_va(LBi:UBi,LBj:UBj,N(ng))
#  endif
!
!  Local variable declarations.
!
      integer :: i, j, k
!
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: tl_Urho, tl_Vrho

#  include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Transform vector components from C-grid to A-grid.
!-----------------------------------------------------------------------
!
      K_LOOP : DO k=1,N(ng)
!
!  Compute A-grid (cell center) components. Apply gradient condition.
!
        DO j=JstrR,JendR
          DO i=Istr,Iend
!>          Urho(i,j)=0.5_r8*(u(i,j,k,ninp)+u(i+1,j,k,ninp))
!>
            tl_Urho(i,j)=0.5_r8*(tl_u(i,j,k,ninp)+tl_u(i+1,j,k,ninp))
            IF (.not.EWperiodic(ng)) THEN
              IF (DOMAIN(ng)%Western_Edge(tile)) THEN
!>              Urho(Istr-1,j)=Urho(Istr,j)
!>
                tl_Urho(Istr-1,j)=tl_Urho(Istr,j)
              END IF
              IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
!>              Urho(Iend+1,j)=Urho(Iend,j)
!>
                tl_Urho(Iend+1,j)=tl_Urho(Iend,j)
              END IF
            END IF
          END DO
        END DO
!
        DO j=Jstr,Jend
          DO i=IstrR,IendR
!>          Vrho(i,j)=0.5_r8*(v(i,j,k,ninp)+v(i,j+1,k,ninp))
!>
            tl_Vrho(i,j)=0.5_r8*(tl_v(i,j,k,ninp)+tl_v(i,j+1,k,ninp))
            IF (.not.NSperiodic(ng)) THEN
              IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
!>              Vrho(i,Jstr-1) = Vrho(i,Jstr)
!>
                tl_Vrho(i,Jstr-1) = tl_Vrho(i,Jstr)
              END IF
              IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
!>              Vrho(i,Jend+1) = Vrho(i,Jend)
!>
                tl_Vrho(i,Jend+1) = tl_Vrho(i,Jend)
              END IF
            END IF
          END DO
        END DO
!
!  Rotate from computational to gegraphical Eastward and Northward
!  directions.
!
        DO j=JstrR,JendR
          DO i=IstrR,IendR
!>          ua(i,j,k)=Urho(i,j)*CosAngler(i,j)-                         &
!>   &                Vrho(i,j)*SinAngler(i,j)
!>
            tl_ua(i,j,k)=tl_Urho(i,j)*CosAngler(i,j)-                   &
     &                   tl_Vrho(i,j)*SinAngler(i,j)
!>          va(i,j,k)=Vrho(i,j)*CosAngler(i,j)+                         &
!>   &                Urho(i,j)*SinAngler(i,j)
!>
            tl_va(i,j,k)=tl_Vrho(i,j)*CosAngler(i,j)+                   &
     &                   tl_Urho(i,j)*SinAngler(i,j)
# ifdef MASKING
!>          ua(i,j,k)=ua(i,j,k)*rmask(i,j)
!>
            tl_ua(i,j,k)=tl_ua(i,j,k)*rmask(i,j)
!>          va(i,j,k)=va(i,j,k)*rmask(i,j)
!>
            tl_va(i,j,k)=tl_va(i,j,k)*rmask(i,j)
# endif
          END DO
        END DO
      END DO K_LOOP

#  ifdef DISTRIBUTE
!
!  Exchange boundary data.
!
!>    CALL mp_exchange3d (ng, tile, model, 2,                           &
!>   &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
!>   &                    NghostPoints,                                 &
!>   &                    EWperiodic(ng), NSperiodic(ng),               &
!>   &                    ua, va)
!>
      CALL mp_exchange3d (ng, tile, model, 2,                           &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    tl_ua, tl_va)
#  endif
!
      RETURN
      END SUBROUTINE tl_uv_C2A_grid_tile    
!
!***********************************************************************
      SUBROUTINE tl_uv_A2C_grid (ng, tile, model, nout)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model, nout
!
!  Local variable declarations.
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__
!
#  include "tile.h"
!
#  ifdef PROFILE
      CALL wclock_on (ng, model, 34, __LINE__, MyFile)
#  endif
      CALL tl_uv_A2C_grid_tile (ng, tile, model, nout,                  &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          IminS, ImaxS, JminS, JmaxS,             &
#  ifdef MASKING
     &                          GRID(ng) % umask_full,                  &
     &                          GRID(ng) % vmask_full,                  &
#  endif
     &                          GRID(ng) % CosAngler,                   &
     &                          GRID(ng) % SinAngler,                   &
     &                          OCEAN(ng) % tl_ua,                      &
     &                          OCEAN(ng) % tl_va,                      &
     &                          OCEAN(ng) % tl_u,                       &
     &                          OCEAN(ng) % tl_v)
#  ifdef PROFILE
      CALL wclock_off (ng, model, 34, __LINE__, MyFile)
#  endif
!
      RETURN
      END SUBROUTINE tl_uv_A2C_grid
!
!***********************************************************************
      SUBROUTINE tl_uv_A2C_grid_tile (ng, tile, model, nout,            &
     &                                LBi, UBi, LBj, UBj,               &
     &                                IminS, ImaxS, JminS, JmaxS,       &
#  ifdef MASKING
     &                                umask, vmask,                     &
#  endif
     &                                CosAngler, SinAngler,             &
     &                                tl_ua, tl_va, tl_u, tl_v)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model, nout
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
#  ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: CosAngler(LBi:,LBj:)
      real(r8), intent(in) :: SinAngler(LBi:,LBj:)
#   ifdef MASKING
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#   endif
      real(r8), intent(in) :: tl_ua(LBi:,LBj:,:)
      real(r8), intent(in) :: tl_va(LBi:,LBj:,:)

      real(r8), intent(inout) :: tl_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: tl_v(LBi:,LBj:,:,:)
#  else
      real(r8), intent(in) :: CosAngler(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: SinAngler(LBi:UBi,LBj:UBj)
#   ifdef MASKING
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#   endif
      real(r8), intent(in) :: tl_ua(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: tl_va(LBi:UBi,LBj:UBj,N(ng))

      real(r8), intent(inout) :: tl_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: tl_v(LBi:UBi,LBj:UBj,N(ng),2)
#  endif
!
!  Local variable declarations.
!
      integer :: i, j, k
!
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: tl_Urho, tl_Vrho

#  include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Transform vector components from A-grid to C-grid.
!-----------------------------------------------------------------------
!
      K_LOOP : DO k=1,N(ng)
!
!  Rotate vector components to computational (XI,ETA) directions.
!
        DO j=Jstr-1,JendR
          DO i=Istr-1,IendR
!>          Urho(i,j)=ua(i,j,k)*CosAngler(i,j)+                         &
!>   &                va(i,j,k)*SinAngler(i,j)
!>
            tl_Urho(i,j)=tl_ua(i,j,k)*CosAngler(i,j)+                   &
     &                   tl_va(i,j,k)*SinAngler(i,j)
!>          Vrho(i,j)=va(i,j,k)*CosAngler(i,j)-                         &
!>   &                ua(i,j,k)*SinAngler(i,j)
!>
            tl_Vrho(i,j)=tl_va(i,j,k)*CosAngler(i,j)-                   &
     &                   tl_ua(i,j,k)*SinAngler(i,j)
          END DO
        END DO
!
!  Compute staggered C-grid components.
!
        DO j=JstrR,JendR
          DO i=Istr,IendR
!>          u(i,j,k,nout)=0.5_r8*(Urho(i-1,j)+Urho(i,j))
!>
            tl_u(i,j,k,nout)=0.5_r8*(tl_Urho(i-1,j)+tl_Urho(i,j))
#  ifdef MASKING
!>          u(i,j,k,nout)=u(i,j,k,nout)*umask(i,j)
!>
            tl_u(i,j,k,nout)=tl_u(i,j,k,nout)*umask(i,j)
#  endif
          END DO
        END DO
!
        DO j=Jstr,JendR
          DO i=IstrR,IendR
!>          v(i,j,k,nout)=0.5_r8*(Vrho(i,j-1)+Vrho(i,j))
!>
            tl_v(i,j,k,nout)=0.5_r8*(tl_Vrho(i,j-1)+tl_Vrho(i,j))
#  ifdef MASKING
!>          v(i,j,k,nout)=v(i,j,k,nout)*vmask(i,j)
!>
            tl_v(i,j,k,nout)=tl_v(i,j,k,nout)*vmask(i,j)
#  endif
          END DO
        END DO
      END DO K_LOOP

#  ifdef DISTRIBUTE
!
!  Exchange boundary data.
!
!>    CALL mp_exchange3d (ng, tile, model, 2,                           &
!>   &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
!>   &                    NghostPoints,                                 &
!>   &                    EWperiodic(ng), NSperiodic(ng),               &
!>   &                    u(:,:,:,nout), v(:,:,:,nout))
!>
      CALL mp_exchange3d (ng, tile, model, 2,                           &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    tl_u(:,:,:,nout), tl_v(:,:,:,nout))
#  endif
!
      RETURN
      END SUBROUTINE tl_uv_A2C_grid_tile    
# endif

# ifdef ADJOINT
!
!***********************************************************************
      SUBROUTINE ad_uv_C2A_grid (ng, tile, model, ninp)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model, ninp
!
!  Local variable declarations.
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__
!
#  include "tile.h"
!
#  ifdef PROFILE
      CALL wclock_on (ng, model, 34, __LINE__, MyFile)
#  endif
      CALL ad_uv_C2A_grid_tile (ng, tile, model, ninp,                  &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          IminS, ImaxS, JminS, JmaxS,             &
#  ifdef MASKING
     &                          GRID(ng) % rmask_full,                  &
#  endif
     &                          GRID(ng) % CosAngler,                   &
     &                          GRID(ng) % SinAngler,                   &
     &                          OCEAN(ng) % ad_u,                       &
     &                          OCEAN(ng) % ad_v,                       &
     &                          OCEAN(ng) % ad_ua,                      &
     &                          OCEAN(ng) % ad_va)
#  ifdef PROFILE
      CALL wclock_off (ng, model, 34, __LINE__, MyFile)
#  endif
!
      RETURN
      END SUBROUTINE ad_uv_C2A_grid
!
!***********************************************************************
      SUBROUTINE ad_uv_C2A_grid_tile (ng, tile, model, ninp,            &
     &                                LBi, UBi, LBj, UBj,               &
     &                                IminS, ImaxS, JminS, JmaxS,       &
#  ifdef MASKING
     &                                rmask,                            &
#  endif
     &                                CosAngler, SinAngler,             &
     &                                ad_u, ad_v, ad_ua, ad_va)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model, ninp
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
#  ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: CosAngler(LBi:,LBj:)
      real(r8), intent(in) :: SinAngler(LBi:,LBj:)
#   ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
#   endif
      real(r8), intent(inout) :: ad_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: ad_v(LBi:,LBj:,:,:)

      real(r8), intent(inout) :: ad_ua(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_va(LBi:,LBj:,:)
#  else
      real(r8), intent(in) :: CosAngler(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: SinAngler(LBi:UBi,LBj:UBj)
#   ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
#   endif
      real(r8), intent(inout) :: ad_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: ad_v(LBi:UBi,LBj:UBj,N(ng),2)

      real(r8), intent(inout) :: ad_ua(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: tl_va(LBi:UBi,LBj:UBj,N(ng))
#  endif
!
!  Local variable declarations.
!
      integer :: i, j, k
!
      real(r8) :: adfac, adfac1, adfac2
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ad_Urho, ad_Vrho

#  include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Adjoint of transform vector components from C-grid to A-grid.
!-----------------------------------------------------------------------
!
!  Initialize.
!
      adfac=0.0_r8
      adfac1=0.0_r8
      adfac2=0.0_r8
      ad_Urho=0.0_r8
      ad_Vrho=0.0_r8

#  ifdef DISTRIBUTE
!
!  Adjoint of exchange boundary data.
!
!>    CALL mp_exchange3d (ng, tile, model, 2,                           &
!>   &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
!>   &                    NghostPoints,                                 &
!>   &                    EWperiodic(ng), NSperiodic(ng),               &
!>   &                    tl_ua, tl_va)
!>
      CALL ad_mp_exchange3d (ng, tile, model, 2,                        &
     &                       LBi, UBi, LBj, UBj, 1, N(ng),              &
     &                       NghostPoints,                              &
     &                       EWperiodic(ng), NSperiodic(ng),            &
     &                       ad_ua, ad_va)
#  endif
!
!  Adjoint of variable change from C-grid to A-grid.
!
      K_LOOP : DO k=1,N(ng)
!
!  Adjoint of rotate from computational to gegraphical Eastward and
!  Northward directions.
!
        DO j=JstrR,JendR
          DO i=IstrR,IendR
# ifdef MASKING
!>          tl_va(i,j,k)=tl_va(i,j,k)*rmask(i,j)
!>
            ad_va(i,j,k)=ad_va(i,j,k)*rmask(i,j)
!>          tl_ua(i,j,k)=tl_ua(i,j,k)*rmask(i,j)
!>
            ad_ua(i,j,k)=ad_ua(i,j,k)*rmask(i,j)
# endif
!>          tl_va(i,j,k)=tl_Vrho(i,j)*CosAngler(i,j)+                   &
!>   &                   tl_Urho(i,j)*SinAngler(i,j)
!>
            adfac1=CosAngler(i,j)*ad_va(i,j,k)
            adfac2=SinAngler(i,j)*ad_va(i,j,k)
            ad_Vrho(i,j)=ad_Vrho(i,j)+adfac1
            ad_Urho(i,j)=ad_Urho(i,j)+adfac2
            ad_va(i,j,k)=0.0_r8
!>          tl_ua(i,j,k)=tl_Urho(i,j)*CosAngler(i,j)                    &
!>   &                   tl_Vrho(i,j)*SinAngler(i,j)
!>
            adfac1=CosAngler(i,j)*ad_ua(i,j,k)
            adfac2=SinAngler(i,j)*ad_ua(i,j,k)
            ad_Urho=ad_Urho+adfac1
            ad_Vrho=ad_Vrho-adfac2
            ad_ua(i,j,k)=0.0_r8
          END DO
        END DO
!
!  Adjoint of compute A-grid (cell center) vector components.
!
        DO j=Jstr,Jend
          DO i=IstrR,IendR
            IF (.not.NSperiodic(ng)) THEN
              IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
!>              tl_Vrho(i,Jend+1)=tl_Vrho(i,Jend)
!>
                ad_Vrho(i,Jend  )=ad_Vrho(i,Jend)+ad_Vrho(i,Jend+1)
                ad_Vrho(i,Jend+1)=0.0_r8
              END IF
              IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
!>              tl_Vrho(i,Jstr-1)=tl_Vrho(i,Jstr)
!>
                ad_Vrho(i,Jstr  )=ad_Vrho(i,Jstr)+ad_Vrho(i,Jstr-1)
                ad_Vrho(i,Jstr-1)=0.0_r8
              END IF
            END IF
!>          tl_Vrho(i,j)=0.5_r8*(tl_v(i,j,k,ninp)+tl_v(i,j+1,k,ninp))
!>
            adfac=0.5_r8*ad_Vrho(i,j)
            ad_v(i,j  ,k,ninp)=ad_v(i,j  ,k,ninp)+adfac
            ad_v(i,j+1,k,ninp)=ad_v(i,j+1,k,ninp)+adfac
            ad_Vrho(i,j)=0.0_r8
          END DO
        END DO
!
        DO j=JstrR,JendR
          DO i=Istr,Iend
            IF (.not.EWperiodic(ng)) THEN
              IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
!>              tl_Urho(Iend+1,j)=tl_Urho(Iend,j)
!>
                ad_Urho(Iend  ,j)=ad_Urho(Iend,j)+ad_Urho(Iend+1,j)
                ad_Urho(Iend+1,j)=0.0_r8
              END IF
              IF (DOMAIN(ng)%Western_Edge(tile)) THEN
!>              tl_Urho(Istr-1,j)=tl_Urho(Istr,j)
!>
                ad_Urho(Istr  ,j)=ad_Urho(Istr,j)+ad_Urho(Istr-1,j)
                ad_Urho(Istr-1,j)=0.0_r8
              END IF
            END IF
!>          tl_Urho(i,j)=0.5_r8*(tl_u(i,j,k,ninp)+tl_u(i+1,j,k,ninp))
!>
            adfac=0.5_r8*ad_Urho(i,j)
            ad_u(i  ,j,k,ninp)=ad_u(i  ,j,k,ninp)+adfac
            ad_u(i+1,j,k,ninp)=ad_u(i+1,j,k,ninp)+adfac
            ad_Urho(i,j)=0.0_r8
          END DO
        END DO
      END DO K_LOOP
!
      RETURN
      END SUBROUTINE ad_uv_C2A_grid_tile    
!
!***********************************************************************
      SUBROUTINE ad_uv_A2C_grid (ng, tile, model, nout)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model, nout
!
!  Local variable declarations.
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__
!
#  include "tile.h"
!
#  ifdef PROFILE
      CALL wclock_on (ng, model, 34, __LINE__, MyFile)
#  endif
      CALL ad_uv_A2C_grid_tile (ng, tile, model, nout,                  &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          IminS, ImaxS, JminS, JmaxS,             &
#  ifdef MASKING
     &                          GRID(ng) % umask,                       &
     &                          GRID(ng) % vmask,                       &
#  endif
     &                          GRID(ng) % CosAngler,                   &
     &                          GRID(ng) % SinAngler,                   &
     &                          OCEAN(ng) % ad_ua,                      &
     &                          OCEAN(ng) % ad_va,                      &
     &                          OCEAN(ng) % ad_u,                       &
     &                          OCEAN(ng) % ad_v)
#  ifdef PROFILE
      CALL wclock_off (ng, model, 34, __LINE__, MyFile)
#  endif
!
      RETURN
      END SUBROUTINE ad_uv_A2C_grid
!
!***********************************************************************
      SUBROUTINE ad_uv_A2C_grid_tile (ng, tile, model, nout,            &
     &                                LBi, UBi, LBj, UBj,               &
     &                                IminS, ImaxS, JminS, JmaxS,       &
#  ifdef MASKING
     &                                umask, vmask,                     &
#  endif
     &                                CosAngler, SinAngler,             &
     &                                ad_ua, ad_va, ad_u, ad_v)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model, nout
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
#  ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: CosAngler(LBi:,LBj:)
      real(r8), intent(in) :: SinAngler(LBi:,LBj:)
#   ifdef MASKING
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#   endif
      real(r8), intent(inout) :: ad_ua(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_va(LBi:,LBj:,:)

      real(r8), intent(inout) :: ad_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: ad_v(LBi:,LBj:,:,:)
#  else
      real(r8), intent(in) :: CosAngler(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: SinAngler(LBi:UBi,LBj:UBj)
#   ifdef MASKING
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#   endif
      real(r8), intent(inout) :: ad_ua(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: ad_va(LBi:UBi,LBj:UBj,N(ng))

      real(r8), intent(inout) :: ad_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: ad_v(LBi:UBi,LBj:UBj,N(ng),2)
#  endif
!
!  Local variable declarations.
!
      integer :: i, j, k
!
      real(r8) :: adfac, adfac1, adfac2
!
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ad_Urho, ad_Vrho

#  include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Adjoint of transform vector components from A-grid to C-grid.
!-----------------------------------------------------------------------
!
!  Initialize.
!
      adfac=0.0_r8
      adfac1=0.0_r8
      adfac2=0.0_r8
      ad_Urho=0.0_r8
      ad_Vrho=0.0_r8

#  ifdef DISTRIBUTE
!
!  Adjoint of exchange boundary data.
!
!>    CALL mp_exchange3d (ng, tile, model, 2,                           &
!>   &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
!>   &                    NghostPoints,                                 &
!>   &                    EWperiodic(ng), NSperiodic(ng),               &
!>   &                    tl_u(:,:,:,nout), tl_v(:,:,:,nout))
!>
      CALL ad_mp_exchange3d (ng, tile, model, 2,                        &
     &                       LBi, UBi, LBj, UBj, 1, N(ng),              &
     &                       NghostPoints,                              &
     &                       EWperiodic(ng), NSperiodic(ng),            &
     &                       ad_u(:,:,:,nout), ad_v(:,:,:,nout))
#  endif
!
!  Adjoint of variable change from A-grid to C-grid.
!
      K_LOOP :DO k=1,N(ng)
!
!  Adjoint of compute staggered C-grid components.
!
        DO j=Jstr,JendR
          DO i=IstrR,IendR
#  ifdef MASKING
!>          tl_v(i,j,k,nout)=tl_v(i,j,k,nout)*vmask(i,j)
!>
            ad_v(i,j,k,nout)=ad_v(i,j,k,nout)*vmask(i,j)
#  endif
!>          tl_v(i,j,k,nout)=0.5_r8*(tl_Vrho(i,j-1)+tl_Vrho(i,j))
!>
            adfac=0.5_r8*ad_v(i,j,k,nout)
            ad_Vrho(i,j-1)=ad_Vrho(i,j-1)+adfac
            ad_Vrho(i,j  )=ad_Vrho(i,j  )+adfac
            ad_v(i,j,k,nout)=0.0_r8
          END DO
        END DO
!
        DO j=JstrR,JendR
          DO i=Istr,IendR
#  ifdef MASKING
!>          tl_u(i,j,k,nout)=tl_u(i,j,k,nout)*umask(i,j)
!>
            ad_u(i,j,k,nout)=ad_u(i,j,k,nout)*umask(i,j)
#  endif
!>          tl_u(i,j,k,nout)=0.5_r8*(tl_Urho(i-1,j)+tl_Urho(i,j))
!>
            adfac=0.5_r8*ad_u(i,j,k,nout)
            ad_Urho(i-1,j)=ad_Urho(i-1,j)+adfac
            ad_Urho(i  ,j)=ad_Urho(i  ,j)+adfac
            ad_u(i,j,k,nout)=0.0_r8
          END DO
        END DO
!
!  Adjoint of rotate vector components to computations (XI,ETA)
!  directions.
!
        DO j=Jstr-1,JendR
          DO i=Istr-1,IendR
!>          tl_Vrho(i,j)=tl_va(i,j,k)*CosAngler(i,j)-                   &
!>   &                   tl_ua(i,j,k)*SinAngler(i,j)
!>
            adfac1=CosAngler(i,j)*ad_Vrho(i,j)
            adfac2=SinAngler(i,j)*ad_Vrho(i,j)
            ad_va(i,j,k)=ad_va(i,j,k)+adfac1
            ad_ua(i,j,k)=ad_ua(i,j,k)-adfac2
            ad_Vrho(i,j)=0.0_r8
!>          tl_Urho(i,j)=tl_ua(i,j,k)*CosAngler(i,j)+                   &
!>   &                   tl_va(i,j,k)*SinAngler(i,j)
!>
            adfac1=CosAngler(i,j)*ad_Urho(i,j)
            adfac2=SinAngler(i,j)*ad_Urho(i,j)
            ad_ua(i,j,k)=ad_ua(i,j,k)+adfac1
            ad_va(i,j,k)=ad_va(i,j,k)+adfac2
            ad_Urho(i,j)=0.0_r8
          END DO
        END DO
      END DO K_LOOP
!
      RETURN
      END SUBROUTINE ad_uv_A2C_grid_tile    
# endif
#endif
      END MODULE uv_var_change_mod
