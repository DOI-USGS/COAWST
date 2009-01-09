#include "cppdefs.h"
      MODULE set_3dchild_data_mod
#if defined REFINED_GRID && defined SOLVE3D
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
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: set_3dchild_data, set_3dchild_data_tile

      CONTAINS
!
!***********************************************************************
      SUBROUTINE set_3dchild_data (ng, tile)
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
      integer :: ngc, LBic, UBic, LBjc, UBjc
# include "tile.h"

      ngc=mychild(ng)
      LBic=BOUNDS(ngc)%LBi(tile)
      UBic=BOUNDS(ngc)%UBi(tile)
      LBjc=BOUNDS(ngc)%LBj(tile)
      UBjc=BOUNDS(ngc)%UBj(tile)
!
# ifdef PROFILE
      CALL wclock_on (ng, iNLM, 4)
# endif
      CALL set_3dchild_data_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            LBic, UBic, LBjc, UBjc)
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 4)
# endif
      RETURN
      END SUBROUTINE set_3dchild_data
!
!***********************************************************************
      SUBROUTINE set_3dchild_data_tile (ng, tile,                       &
     &                                  LBi, UBi, LBj, UBj,             &
     &                                  LBic, UBic, LBjc, UBjc)
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
      USE mod_scalars
      USE mod_coupling
      USE interpolate_mod

# ifdef DISTRIBUTE
      USE distribute_mod, ONLY : mp_bcasti, mp_gather
# endif
# if defined EW_PERIODIC || defined NS_PERIODIC
      USE exchange_2d_mod
# endif
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
      USE mp_exchange_mod, ONLY : mp_exchange3d
      USE mp_exchange_mod, ONLY : mp_exchange4d
# endif
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: LBic, UBic, LBjc, UBjc
!
!  Local variable declarations.
!

# ifdef OBC
      integer :: ILB, IUB, JLB, JUB
# endif
      integer :: i, j, k, ip, jp, ngc, ic, jc, t1, t2, itrc
      real(r8) :: my_val, fac, fac1, fac2, fac3
      integer :: ii, gtype, Npts
      integer :: Imin, Imax, Jmin, Jmax, Kmin, Kmax, model
      integer :: Ilen, Jlen, IJlen, Klen, MyType, status
      integer :: MyError, MySize, ncvarid
      integer :: Xmin, Xmax, Ymin, Ymax
      integer :: Istrp, Iendp, Jstrp, Jendp
      real(r8) :: my_min, my_max
!
      real(r8), allocatable :: Aout(:)
      real(r8), allocatable :: A(:,:,:)
      real(r8), allocatable :: AA(:,:,:)
      real(r8), allocatable :: Amask(:,:)
      real(r8), allocatable :: A2(:,:,:)

      logical :: rectangular
      real(r8) :: IJspv
      character (len=20) :: Xname, Yname
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

# include "set_bounds.h"

      ngc=mychild(ng)
      model=iNLM

      Istrc=BOUNDS(ngc)%Istr(tile)
      Iendc=BOUNDS(ngc)%Iend(tile)
      Jstrc=BOUNDS(ngc)%Jstr(tile)
      Jendc=BOUNDS(ngc)%Jend(tile)

!
# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn on input data time wall clock.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, iNLM, 3)
# endif
!
!-----------------------------------------------------------------------
!  Fill local arrays from the BC data.
!-----------------------------------------------------------------------

!
! loop thru tracer, u3, v3 points
!
      DO ii=1,3   ! 3 gtypes
        IF (ii.eq.1) THEN
          gtype=r3dvar
        ELSE IF (ii.eq.2) THEN
          gtype=u3dvar
        ELSE IF (ii.eq.3) THEN
          gtype=v3dvar
        END IF

!       now get the bounds
        IF (ngc.eq.1) THEN
          Imin=0
          Imax=Lm(ngc)+1
          Jmin=0
          Jmax=Mm(ngc)+1
          MyType=gtype
          IF ((ABS(MyType).eq.p3dvar).or.(ABS(MyType).eq.u3dvar)) THEN
            Imin=1
          END IF
          IF ((ABS(MyType).eq.p3dvar).or.(ABS(MyType).eq.v3dvar)) THEN
            Jmin=1
          END IF
        ELSE
          Imin=LOWER_BOUND_I
          Imax=Lm(ngc)+GHOST_POINTS
          Jmin=LOWER_BOUND_J
          Jmax=Mm(ngc)+GHOST_POINTS
          MyType=gtype
          IF ((ABS(MyType).eq.p3dvar).or.(ABS(MyType).eq.u3dvar)) THEN
            Imin=LOWER_BOUND_I+1
          END IF
          IF ((ABS(MyType).eq.p3dvar).or.(ABS(MyType).eq.v3dvar)) THEN
            Jmin=LOWER_BOUND_J+1
          END IF
        END IF
        Kmin=1
        Kmax=N(ngc)
        Ilen=Imax-Imin+1
        Jlen=Jmax-Jmin+1
        Klen=Kmax-Kmin+1
        IJlen=Ilen*Jlen
        MySize=IJlen*Klen
        Npts=MySize
        rectangular=.FALSE.
        IJspv=-9999.0_r8
        ncvarid=1  !not used

        IF (.not.allocated(Aout)) THEN
          allocate (Aout(Npts))
        END IF
        IF (.not.allocated(A)) THEN
          allocate(A(Imin:Imax,Jmin:Jmax,Kmin:Kmax))
        END IF
        IF (.not.allocated(AA)) THEN
          allocate(AA(Imin:Imax,Jmin:Jmax,Kmin:Kmax))
        END IF
        IF (.not.allocated(Amask)) THEN
          allocate(Amask(Imin:Imax,Jmin:Jmax))
        END IF
        IF (.not.allocated(A2)) THEN
          allocate(A2(LBi:UBi,LBj:UBj,Kmin:Kmax))
        END IF
!
!-----------------------------------------------------------------------
!  Loop thru ii = 1 : rho points, 2: u points, 3: vpoints              !
!-----------------------------------------------------------------------
!
        IF (ii.eq.1) THEN
!!!!!!!!!!!!!!!!!!!
!  tracer  !!!!!!!!
!!!!!!!!!!!!!!!!!!!
!
! get mask
!
          ic=0
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              ic=ic+1
# ifdef MASKING
              Amask(i,j)=GRID(ng)%rmask_c(i,j)
# else
              Amask(i,j)=1.0_r8
# endif
            END DO
          END DO
          DO itrc=1,NT(ng)
            Npts=MySize
!
! get tracer
!
            CALL mp_gather (ngc, model, LBic, UBic, LBjc, UBjc,         &
     &                      1, N(ngc), gtype, 1.0_r8,                   &
# ifdef MASKING
     &                      GRID(ngc) % rmask(LBic,LBjc),               &
# endif
     &                      OCEAN(ngc) % t(LBic,LBjc,1,NOUTC,itrc),     &
     &                      Npts, Aout)
!
! distribute the Aout (1D packed) field to all the nodes for ng=2
            CALL mpi_bcast (Aout, MySize, MP_FLOAT, MyMaster,           &
     &                      OCN_COMM_WORLD, MyError)
! unpack and place in global ng=1 2D array
            ic=0
            DO k=Kmin,Kmax
              DO j=Jmin,Jmax
                DO i=Imin,Imax
                  ic=ic+1
                  AA(i,j,k)=Aout(ic)
                  A(i,j,k)=Aout(ic)
                END DO
              END DO
            END DO
            DO k=Kmin,Kmax
            DO j=Jmin+2,Jmax-2
              DO i=Imin+2,Imax-2
                A(i,j,k)=0.04_r8*                                       &
     &         (AA(i-2,j+2,k)+AA(i-1,j+2,k)+AA(i  ,j+2,k)+AA(i+1,j+2,k)+&
     &          AA(i+2,j+2,k)+AA(i-2,j+1,k)+AA(i-1,j+1,k)+AA(i  ,j+1,k)+&
     &          AA(i+1,j+1,k)+AA(i+2,j+1,k)+AA(i-2,j  ,k)+AA(i-1,j  ,k)+&
     &          AA(i  ,j  ,k)+AA(i+1,j  ,k)+AA(i+2,j  ,k)+AA(i-2,j-1,k)+&
     &          AA(i-1,j-1,k)+AA(i  ,j-1,k)+AA(i+1,j-1,k)+AA(i+2,j-1,k)+&
     &          AA(i-2,j-2,k)+AA(i-1,j-2,k)+AA(i  ,j-2,k)+AA(i+1,j-2,k)+&
     &          AA(i+2,j-2,k))
              END DO
            END DO
            END DO
!
! Now interpolate the child data to the parent grid.
!
            Xmin=OCEAN(ng)%hindcr(1)
            Xmax=OCEAN(ng)%hindcr(2)
            Ymin=OCEAN(ng)%hindcr(3)
            Ymax=OCEAN(ng)%hindcr(4)
            IF (spherical) THEN
              CALL linterp3d_ref (ng, Imin, Imax, Jmin, Jmax,           &
     &                            Kmin, Kmax,                           &
     &                            GRID(ng)%lonr_c(Imin:Imax,Jmin:Jmax), &
     &                            GRID(ng)%latr_c(Imin:Imax,Jmin:Jmax), &
     &                            A, Amask,                             &
     &                            LBi, UBi, LBj, UBj, 1, N(ngc),        &
     &                            Xmin, Xmax, Ymin, Ymax,               &
     &                            OCEAN(ng)%hindci_r,OCEAN(ng)%hindcj_r,&
     &                            GRID(ng)%lonr, GRID(ng)%latr, A2,     &
     &                            my_min, my_max)
            ELSE
              CALL linterp3d_ref (ng, Imin, Imax, Jmin, Jmax,           &
     &                            Kmin, Kmax,                           &
     &                            GRID(ng)%xr_c(Imin:Imax,Jmin:Jmax),   &
     &                            GRID(ng)%yr_c(Imin:Imax,Jmin:Jmax),   &
     &                            A, Amask,                             &
     &                            LBi, UBi, LBj, UBj, 1, N(ngc),        &
     &                            Xmin, Xmax, Ymin, Ymax,               &
     &                            OCEAN(ng)%hindci_r,OCEAN(ng)%hindcj_r,&
     &                            GRID(ng)%xr, GRID(ng)%yr, A2,         &
     &                            my_min, my_max)
            END IF
!
! Set the chld data into the parent locations.
!
            DO k=Kmin,Kmax
              DO j=Ymin,Ymax
                DO i=Xmin,Xmax
                  my_val=A2(i,j,k)
# ifdef MASKING
                  my_val=my_val*GRID(ng)%rmask(i,j)
# endif
!                  OCEAN(ng)%t(i,j,k,1,itrc)=my_val
                  OCEAN(ng)%t(i,j,k,nstp(ng),itrc)=my_val
                END DO
              END DO
            END DO
          END DO
        ELSE IF (ii.eq.2) THEN
!!!!!!!!!!!!!!!!!!!
!  u     !!!!!!!!!!
!!!!!!!!!!!!!!!!!!!
!
! get mask
!
          ic=0
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              ic=ic+1
# ifdef MASKING
              Amask(i,j)=GRID(ng)%umask_c(i,j)
# else
              Amask(i,j)=1.0_r8
# endif
            END DO
          END DO
!
! get u
!
          Npts=MySize
          CALL mp_gather (ngc, model, LBic, UBic, LBjc, UBjc, 1, N(ngc),&
     &                    gtype, 1.0_r8,                                &
# ifdef MASKING
     &                    GRID(ngc) % umask(LBic,LBjc),                 &
# endif
     &                    OCEAN(ngc) % u(LBic,LBjc,1,NOUTC),            &
     &                    Npts, Aout)
!
! distribute the Aout (1D packed) field to all the nodes for ng=2
          CALL mpi_bcast (Aout, MySize, MP_FLOAT, MyMaster,             &
     &                    OCN_COMM_WORLD, MyError)
! unpack and place in global ng=1 2D array
          ic=0
          DO k=Kmin,Kmax
            DO j=Jmin,Jmax
              DO i=Imin,Imax
                ic=ic+1
                AA(i,j,k)=Aout(ic)
                A(i,j,k)=Aout(ic)
              END DO
            END DO
          END DO
          DO k=Kmin,Kmax
            DO j=Jmin+2,Jmax-2
              DO i=Imin+2,Imax-2
                A(i,j,k)=0.04_r8*                                       &
     &         (AA(i-2,j+2,k)+AA(i-1,j+2,k)+AA(i  ,j+2,k)+AA(i+1,j+2,k)+&
     &          AA(i+2,j+2,k)+AA(i-2,j+1,k)+AA(i-1,j+1,k)+AA(i  ,j+1,k)+&
     &          AA(i+1,j+1,k)+AA(i+2,j+1,k)+AA(i-2,j  ,k)+AA(i-1,j  ,k)+&
     &          AA(i  ,j  ,k)+AA(i+1,j  ,k)+AA(i+2,j  ,k)+AA(i-2,j-1,k)+&
     &          AA(i-1,j-1,k)+AA(i  ,j-1,k)+AA(i+1,j-1,k)+AA(i+2,j-1,k)+&
     &          AA(i-2,j-2,k)+AA(i-1,j-2,k)+AA(i  ,j-2,k)+AA(i+1,j-2,k)+&
     &          AA(i+2,j-2,k))
              END DO
            END DO
          END DO
!
! Now interpolate the child data to the parent grid.
!
          Xmin=OCEAN(ng)%hindcu(1)
          Xmax=OCEAN(ng)%hindcu(2)
          Ymin=OCEAN(ng)%hindcu(3)
          Ymax=OCEAN(ng)%hindcu(4)
          IF (spherical) THEN
            CALL linterp3d_ref (ng, Imin, Imax, Jmin, Jmax,             &
     &                          Kmin, Kmax,                             &
     &                          GRID(ng)%lonu_c(Imin:Imax,Jmin:Jmax),   &
     &                          GRID(ng)%latu_c(Imin:Imax,Jmin:Jmax),   &
     &                          A, Amask,                               &
     &                          LBi, UBi, LBj, UBj, 1, N(ngc),          &
     &                          Xmin, Xmax, Ymin, Ymax,                 &
     &                          OCEAN(ng)%hindci_u,OCEAN(ng)%hindcj_u,  &
     &                          GRID(ng)%lonu, GRID(ng)%latu, A2,       &
     &                          my_min, my_max)
          ELSE
            CALL linterp3d_ref (ng, Imin, Imax, Jmin, Jmax,             &
     &                          Kmin, Kmax,                             &
     &                          GRID(ng)%xu_c(Imin:Imax,Jmin:Jmax),     &
     &                          GRID(ng)%yu_c(Imin:Imax,Jmin:Jmax),     &
     &                          A, Amask,                               &
     &                          LBi, UBi, LBj, UBj, 1, N(ngc),          &
     &                          Xmin, Xmax, Ymin, Ymax,                 &
     &                          OCEAN(ng)%hindci_u,OCEAN(ng)%hindcj_u,  &
     &                          GRID(ng)%xu, GRID(ng)%yu, A2,           &
     &                          my_min, my_max)
          END IF
          DO k=Kmin,Kmax
            DO j=Ymin,Ymax
              DO i=Xmin,Xmax
                my_val=A2(i,j,k)
# ifdef MASKING
                my_val=my_val*GRID(ng)%umask(i,j)
# endif
!                OCEAN(ng)%u(i,j,k,1)=my_val
                OCEAN(ng)%u(i,j,k,nstp(ng))=my_val
              END DO
            END DO
          END DO

        ELSE IF (ii.eq.3) THEN
!!!!!!!!!!!!!!!!!!!
!  v     !!!!!!!!!!
!!!!!!!!!!!!!!!!!!!
!
! get mask
!
          ic=0
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              ic=ic+1
# ifdef MASKING
              Amask(i,j)=GRID(ng)%vmask_c(i,j)
# else
              Amask(i,j)=1.0_r8
# endif
            END DO
          END DO
!
! get v
!
          Npts=MySize
          CALL mp_gather (ngc, model, LBic, UBic, LBjc, UBjc, 1, N(ngc),&
     &                    gtype, 1.0_r8,                                &
# ifdef MASKING
     &                    GRID(ngc) % vmask(LBic,LBjc),                 &
# endif
     &                    OCEAN(ngc) % v(LBic,LBjc,1,NOUTC),            &
     &                    Npts, Aout)
!
! distribute the Aout (1D packed) field to all the nodes for ng=2
          CALL mpi_bcast (Aout, MySize, MP_FLOAT, MyMaster,             &
     &                    OCN_COMM_WORLD, MyError)
! unpack and place in global ng=1 2D array
          ic=0
          DO k=Kmin,Kmax
            DO j=Jmin,Jmax
              DO i=Imin,Imax
                ic=ic+1
                AA(i,j,k)=Aout(ic)
                A(i,j,k)=Aout(ic)
              END DO
            END DO
          END DO
          DO k=Kmin,Kmax
            DO j=Jmin+2,Jmax-2
              DO i=Imin+2,Imax-2
                A(i,j,k)=0.04_r8*                                       &
     &         (AA(i-2,j+2,k)+AA(i-1,j+2,k)+AA(i  ,j+2,k)+AA(i+1,j+2,k)+&
     &          AA(i+2,j+2,k)+AA(i-2,j+1,k)+AA(i-1,j+1,k)+AA(i  ,j+1,k)+&
     &          AA(i+1,j+1,k)+AA(i+2,j+1,k)+AA(i-2,j  ,k)+AA(i-1,j  ,k)+&
     &          AA(i  ,j  ,k)+AA(i+1,j  ,k)+AA(i+2,j  ,k)+AA(i-2,j-1,k)+&
     &          AA(i-1,j-1,k)+AA(i  ,j-1,k)+AA(i+1,j-1,k)+AA(i+2,j-1,k)+&
     &          AA(i-2,j-2,k)+AA(i-1,j-2,k)+AA(i  ,j-2,k)+AA(i+1,j-2,k)+&
     &          AA(i+2,j-2,k))
              END DO
            END DO
          END DO
!
! Now interpolate the child data to the parent grid.
!
          Xmin=OCEAN(ng)%hindcv(1)
          Xmax=OCEAN(ng)%hindcv(2)
          Ymin=OCEAN(ng)%hindcv(3)
          Ymax=OCEAN(ng)%hindcv(4)
          IF (spherical) THEN
            CALL linterp3d_ref (ng, Imin, Imax, Jmin, Jmax,             &
     &                          Kmin, Kmax,                             &
     &                          GRID(ng)%lonv_c(Imin:Imax,Jmin:Jmax),   &
     &                          GRID(ng)%latv_c(Imin:Imax,Jmin:Jmax),   &
     &                          A, Amask,                               &
     &                          LBi, UBi, LBj, UBj, 1, N(ngc),          &
     &                          Xmin, Xmax, Ymin, Ymax,                 &
     &                          OCEAN(ng)%hindci_v,OCEAN(ng)%hindcj_v,  &
     &                          GRID(ng)%lonv, GRID(ng)%latv, A2,       &
     &                          my_min, my_max)
          ELSE
            CALL linterp3d_ref (ng, Imin, Imax, Jmin, Jmax,             &
     &                          Kmin, Kmax,                             &
     &                          GRID(ng)%xv_c(Imin:Imax,Jmin:Jmax),     &
     &                          GRID(ng)%yv_c(Imin:Imax,Jmin:Jmax),     &
     &                          A, Amask,                               &
     &                          LBi, UBi, LBj, UBj, 1, N(ngc),          &
     &                          Xmin, Xmax, Ymin, Ymax,                 &
     &                          OCEAN(ng)%hindci_v,OCEAN(ng)%hindcj_v,  &
     &                          GRID(ng)%xv, GRID(ng)%yv, A2,           &
     &                          my_min, my_max)
          END IF
          DO k=Kmin,Kmax
            DO j=Ymin,Ymax
              DO i=Xmin,Xmax
                my_val=A2(i,j,k)
# ifdef MASKING
                my_val=my_val*GRID(ng)%vmask(i,j)
# endif
                OCEAN(ng)%v(i,j,k,nstp(ng))=my_val
!                OCEAN(ng)%v(i,j,k,2)=my_val
              END DO
            END DO
          END DO
        END IF
!
        IF (allocated(Aout)) THEN
          deallocate (Aout)
        END IF
        IF (allocated(A)) THEN
          deallocate (A)
        END IF
        IF (allocated(AA)) THEN
          deallocate (AA)
        END IF
        IF (allocated(Amask)) THEN
          deallocate (Amask)
        END IF
        IF (allocated(A2)) THEN
          deallocate (A2)
        END IF
      END DO

# if defined EW_PERIODIC || defined NS_PERIODIC
!
!-----------------------------------------------------------------------
!  Exchange boundary information.
!-----------------------------------------------------------------------
!
      CALL exchange_r3d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        OCEAN(ng)%t(:,:,:,1,itrc))
      CALL exchange_r3d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        OCEAN(ng)%t(:,:,:,2,itrc))
      CALL exchange_u3d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        OCEAN(ng)%u(:,:,:,1))
      CALL exchange_u3d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        OCEAN(ng)%u(:,:,:,2))
      CALL exchange_u3d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        OCEAN(ng)%v(:,:,:,1))
      CALL exchange_u3d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        OCEAN(ng)%v(:,:,:,2))
#  endif
#  ifdef DISTRIBUTE
      CALL mp_exchange4d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 1, N(ng), 1, NT(ng),      &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    OCEAN(ng)%t(:,:,:,1,:))
      CALL mp_exchange4d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 1, N(ng), 1, NT(ng),      &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    OCEAN(ng)%t(:,:,:,2,:))
      CALL mp_exchange3d (ng, tile, iNLM, 4,                            &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    OCEAN(ng)%u(:,:,:,1), OCEAN(ng)%u(:,:,:,2),   &
     &                    OCEAN(ng)%v(:,:,:,1), OCEAN(ng)%v(:,:,:,2))
#  endif
!
# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn off input data time wall clock.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, iNLM, 3)
# endif
      RETURN
      END SUBROUTINE set_3dchild_data_tile
#endif
      END MODULE set_3dchild_data_mod
