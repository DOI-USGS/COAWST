#include "cppdefs.h"
      MODULE set_2dchild_data_mod
#ifdef REFINED_GRID
!
!svn $Id: set_data.F 1405 2008-02-20 01:44:07Z jcwarner $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine acquires date from a child grid.                       !
!  That data is interpolated to the parent grid.                       !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: set_2dchild_data, set_2dchild_data_tile

      CONTAINS
!
!***********************************************************************
      SUBROUTINE set_2dchild_data (ng, tile)
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
      CALL set_2dchild_data_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            LBic, UBic, LBjc, UBjc)
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 4)
# endif
      RETURN
      END SUBROUTINE set_2dchild_data
!
!***********************************************************************
      SUBROUTINE set_2dchild_data_tile (ng, tile,                       &
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
      USE mod_refined
      USE mod_scalars
      USE mod_coupling
      USE interpolate_mod

# ifdef DISTRIBUTE
      USE distribute_mod, ONLY : mp_bcasti, mp_gather2d
# endif
# if defined EW_PERIODIC || defined NS_PERIODIC
      USE exchange_2d_mod
# endif
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
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
      integer :: i, j, ip, jp, ngc, ic, jc, t1, t2
      real(r8) :: my_val, fac, fac1, fac2, fac3
      integer :: ii, gtype, Npts
      integer :: Imin, Imax, Jmin, Jmax, model
      integer :: Ilen, Jlen, IJlen, MyType, status
      integer :: MyError, MySizer, ncvarid
      integer :: Xmin, Xmax, Ymin, Ymax
      integer :: Istrp, Iendp, Jstrp, Jendp
      real(r8) :: my_min, my_max
!
      real(r8), allocatable :: Aout(:)
      real(r8), allocatable :: A(:,:)
      real(r8), allocatable :: AA(:,:)
      real(r8), allocatable :: Amask(:,:)
      real(r8), allocatable :: A2(:,:)

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

#define YESUSE
#ifdef YESUSE
!
! loop thru rho, u2, v2 points
!
      DO ii=1,3   ! 3 gtypes
        IF (ii.eq.1) THEN
          gtype=r2dvar
        ELSE IF (ii.eq.2) THEN
          gtype=u2dvar
        ELSE IF (ii.eq.3) THEN
          gtype=v2dvar
        END IF

!       now get the bounds
        IF (ngc.eq.1) THEN
          Imin=0
          Imax=Lm(ngc)+1
          Jmin=0
          Jmax=Mm(ngc)+1
          MyType=gtype
          IF ((ABS(MyType).eq.p2dvar).or.(ABS(MyType).eq.u2dvar)) THEN
            Imin=1
          END IF
          IF ((ABS(MyType).eq.p2dvar).or.(ABS(MyType).eq.v2dvar)) THEN
            Jmin=1
          END IF
        ELSE
          Imin=LOWER_BOUND_I
          Imax=Lm(ngc)+GHOST_POINTS
          Jmin=LOWER_BOUND_J
          Jmax=Mm(ngc)+GHOST_POINTS
          MyType=gtype
          IF ((ABS(MyType).eq.p2dvar).or.(ABS(MyType).eq.u2dvar)) THEN
            Imin=LOWER_BOUND_I+1
          END IF
          IF ((ABS(MyType).eq.p2dvar).or.(ABS(MyType).eq.v2dvar)) THEN
            Jmin=LOWER_BOUND_J+1
          END IF
        END IF
        Ilen=Imax-Imin+1
        Jlen=Jmax-Jmin+1
        IJlen=Ilen*Jlen
        MySizer=IJlen
        Npts=IJlen

        rectangular=.FALSE.
        IJspv=-9999.0_r8
        ncvarid=1  !not used

        IF (.not.allocated(Aout)) THEN
          allocate (Aout(IJlen))
        END IF
        IF (.not.allocated(A)) THEN
          allocate(A(Imin:Imax,Jmin:Jmax))
        END IF
        IF (.not.allocated(AA)) THEN
          allocate(AA(Imin:Imax,Jmin:Jmax))
        END IF
        IF (.not.allocated(Amask)) THEN
          allocate(Amask(Imin:Imax,Jmin:Jmax))
        END IF
        IF (.not.allocated(A2)) THEN
          allocate(A2(LBi:UBi,LBj:UBj))
        END IF
!
!-----------------------------------------------------------------------
!  Loop thru ii = 1 : rho points, 2: u points, 3: vpoints              !
!-----------------------------------------------------------------------
!
        IF (ii.eq.1) THEN
!!!!!!!!!!!!!!!!!!!
!  zeta  !!!!!!!!!!
!!!!!!!!!!!!!!!!!!!
!
! get mask
!
          ic=0
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              ic=ic+1
# ifdef MASKING
              Amask(i,j)=REFINED(ng)%rmask_c(i,j)
# else
              Amask(i,j)=1.0_r8
# endif
            END DO
          END DO
!
! get zeta
!
          CALL mp_gather2d (ngc, model, LBic, UBic, LBjc, UBjc, 0,      &
     &                    gtype, 1.0_r8,                                &
# ifdef MASKING
     &                    GRID(ngc) % rmask(LBic:UBic,LBjc:UBjc),       &
# endif
# ifdef SOLVE3D
     &                    COUPLING(ngc) % Zt_avg1(LBic:UBic,LBjc:UBjc), &
# else
     &               OCEAN(ngc) % zeta(LBic:UBic,LBjc:UBjc,KOUTC),      &
# endif
     &                    Npts, Aout)
!
! distribute the Aout (1D packed) field to all the nodes for ng=2
          CALL mpi_bcast (Aout, MySizer, MP_FLOAT, MyMaster,            &
     &                    OCN_COMM_WORLD, MyError)
! unpack and place in global ng=1 2D array
          ic=0
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              ic=ic+1
              AA(i,j)=Aout(ic)
              A(i,j)=Aout(ic)
            END DO
          END DO
          DO j=Jmin+2,Jmax-2
            DO i=Imin+2,Imax-2
              A(i,j)=0.04_r8*                                           &
     &     (AA(i-2,j+2)+AA(i-1,j+2)+AA(i  ,j+2)+AA(i+1,j+2)+AA(i+2,j+2)+&
     &      AA(i-2,j+1)+AA(i-1,j+1)+AA(i  ,j+1)+AA(i+1,j+1)+AA(i+2,j+1)+&
     &      AA(i-2,j  )+AA(i-1,j  )+AA(i  ,j  )+AA(i+1,j  )+AA(i+2,j  )+&
     &      AA(i-2,j-1)+AA(i-1,j-1)+AA(i  ,j-1)+AA(i+1,j-1)+AA(i+2,j-1)+&
     &      AA(i-2,j-2)+AA(i-1,j-2)+AA(i  ,j-2)+AA(i+1,j-2)+AA(i+2,j-2))
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
            CALL linterp2d_ref (ng, Imin, Imax, Jmin, Jmax,             &
     &                          REFINED(ng)%lonr_c(Imin:Imax,Jmin:Jmax),&
     &                          REFINED(ng)%latr_c(Imin:Imax,Jmin:Jmax),&
     &                          A, Amask,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Xmin, Xmax, Ymin, Ymax,                 &
     &                          OCEAN(ng)%hindci_r, OCEAN(ng)%hindcj_r, &
     &                          GRID(ng)%lonr, GRID(ng)%latr, A2,       &
     &                          my_min, my_max)
          ELSE
            CALL linterp2d_ref (ng, Imin, Imax, Jmin, Jmax,             &
     &                          REFINED(ng)%xr_c(Imin:Imax,Jmin:Jmax),  &
     &                          REFINED(ng)%yr_c(Imin:Imax,Jmin:Jmax),  &
     &                          A, Amask,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Xmin, Xmax, Ymin, Ymax,                 &
     &                          OCEAN(ng)%hindci_r, OCEAN(ng)%hindcj_r, &
     &                          GRID(ng)%xr, GRID(ng)%yr, A2,           &
     &                          my_min, my_max)
          END IF
!
! Set the chld data into the parent locations.
!
          DO j=Ymin,Ymax
            DO i=Xmin,Xmax
              my_val=A2(i,j)
# ifdef MASKING
              my_val=my_val*GRID(ng)%rmask(i,j)
# endif
# ifdef SOLVE3D
!              COUPLING(ng)%Zt_avg1(i,j)=my_val
# else
              OCEAN(ng)%zeta(i,j,1)=my_val
              OCEAN(ng)%zeta(i,j,2)=my_val
              OCEAN(ng)%zeta(i,j,3)=my_val
#endif
            END DO
          END DO
        ELSE IF (ii.eq.2) THEN
!!!!!!!!!!!!!!!!!!!
!  ubar  !!!!!!!!!!
!!!!!!!!!!!!!!!!!!!
!
! get mask
!
          ic=0
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              ic=ic+1
# ifdef MASKING
              Amask(i,j)=REFINED(ng)%umask_c(i,j)
# else
              Amask(i,j)=1.0_r8
# endif
            END DO
          END DO
!
! get ubar
!
          CALL mp_gather2d (ngc, model, LBic, UBic, LBjc, UBjc, 0,      &
     &                    gtype, 1.0_r8,                                &
# ifdef MASKING
     &                    GRID(ngc) % umask(LBic:UBic,LBjc:UBjc),       &
# endif
     &              OCEAN(ngc) % ubar(LBic:UBic,LBjc:UBjc,KOUTC),       &
     &                    Npts, Aout)
!
! distribute the Aout (1D packed) field to all the nodes for ng=2
          CALL mpi_bcast (Aout, MySizer, MP_FLOAT, MyMaster,            &
     &                    OCN_COMM_WORLD, MyError)
! unpack and place in global ng=1 2D array
          ic=0
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              ic=ic+1
              AA(i,j)=Aout(ic)
              A(i,j)=Aout(ic)
            END DO
          END DO
          DO j=Jmin+2,Jmax-2
            DO i=Imin+2,Imax-2
              A(i,j)=0.04_r8*                                           &
     &     (AA(i-2,j+2)+AA(i-1,j+2)+AA(i  ,j+2)+AA(i+1,j+2)+AA(i+2,j+2)+&
     &      AA(i-2,j+1)+AA(i-1,j+1)+AA(i  ,j+1)+AA(i+1,j+1)+AA(i+2,j+1)+&
     &      AA(i-2,j  )+AA(i-1,j  )+AA(i  ,j  )+AA(i+1,j  )+AA(i+2,j  )+&
     &      AA(i-2,j-1)+AA(i-1,j-1)+AA(i  ,j-1)+AA(i+1,j-1)+AA(i+2,j-1)+&
     &      AA(i-2,j-2)+AA(i-1,j-2)+AA(i  ,j-2)+AA(i+1,j-2)+AA(i+2,j-2))
            END DO
          END DO
!
! Now interpolate the child data to the local parent tile.
!
          Xmin=OCEAN(ng)%hindcu(1)
          Xmax=OCEAN(ng)%hindcu(2)
          Ymin=OCEAN(ng)%hindcu(3)
          Ymax=OCEAN(ng)%hindcu(4)
          IF (spherical) THEN
            CALL linterp2d_ref (ng, Imin, Imax, Jmin, Jmax,             &
     &                          REFINED(ng)%lonu_c(Imin:Imax,Jmin:Jmax),&
     &                          REFINED(ng)%latu_c(Imin:Imax,Jmin:Jmax),&
     &                          A, Amask,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Xmin, Xmax, Ymin, Ymax,                 &
     &                          OCEAN(ng)%hindci_u, OCEAN(ng)%hindcj_u, &
     &                          GRID(ng)%lonu, GRID(ng)%latu, A2,       &
     &                          my_min, my_max)
          ELSE
            CALL linterp2d_ref (ng, Imin, Imax, Jmin, Jmax,             &
     &                          REFINED(ng)%xu_c(Imin:Imax,Jmin:Jmax),  &
     &                          REFINED(ng)%yu_c(Imin:Imax,Jmin:Jmax),  &
     &                          A, Amask,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Xmin, Xmax, Ymin, Ymax,                 &
     &                          OCEAN(ng)%hindci_u, OCEAN(ng)%hindcj_u, &
     &                          GRID(ng)%xu, GRID(ng)%yu, A2,           &
     &                          my_min, my_max)
          END IF
          DO j=Ymin,Ymax
            DO i=Xmin,Xmax
              my_val=A2(i,j)
# ifdef MASKING
              my_val=my_val*GRID(ng)%umask(i,j)
# endif
              OCEAN(ng)%ubar(i,j,1)=my_val
              OCEAN(ng)%ubar(i,j,2)=my_val
              OCEAN(ng)%ubar(i,j,3)=my_val
            END DO
          END DO

        ELSE IF (ii.eq.3) THEN
!!!!!!!!!!!!!!!!!!!
!  vbar  !!!!!!!!!!
!!!!!!!!!!!!!!!!!!!
!
! get mask
!
          ic=0
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              ic=ic+1
# ifdef MASKING
              Amask(i,j)=REFINED(ng)%vmask_c(i,j)
# else
              Amask(i,j)=1.0_r8
# endif
            END DO
          END DO
!
! get vbar
!
          CALL mp_gather2d (ngc, model, LBic, UBic, LBjc, UBjc, 0,      &
     &                    gtype, 1.0_r8,                                &
# ifdef MASKING
     &                    GRID(ngc) % vmask(LBic:UBic,LBjc:UBjc),       &
# endif
     &              OCEAN(ngc) % vbar(LBic:UBic,LBjc:UBjc,KOUTC),       &
     &                    Npts, Aout)
!
! distribute the Aout (1D packed) field to all the nodes for ng=2
          CALL mpi_bcast (Aout, MySizer, MP_FLOAT, MyMaster,            &
     &                    OCN_COMM_WORLD, MyError)
! unpack and place in global ng=1 2D array
          ic=0
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              ic=ic+1
              AA(i,j)=Aout(ic)
              A(i,j)=Aout(ic)
            END DO
          END DO
          DO j=Jmin+2,Jmax-2
            DO i=Imin+2,Imax-2
              A(i,j)=0.04_r8*                                           &
     &     (AA(i-2,j+2)+AA(i-1,j+2)+AA(i  ,j+2)+AA(i+1,j+2)+AA(i+2,j+2)+&
     &      AA(i-2,j+1)+AA(i-1,j+1)+AA(i  ,j+1)+AA(i+1,j+1)+AA(i+2,j+1)+&
     &      AA(i-2,j  )+AA(i-1,j  )+AA(i  ,j  )+AA(i+1,j  )+AA(i+2,j  )+&
     &      AA(i-2,j-1)+AA(i-1,j-1)+AA(i  ,j-1)+AA(i+1,j-1)+AA(i+2,j-1)+&
     &      AA(i-2,j-2)+AA(i-1,j-2)+AA(i  ,j-2)+AA(i+1,j-2)+AA(i+2,j-2))
            END DO
          END DO
!
! Now interpolate the child data to the local parent tile.
!
          Xmin=OCEAN(ng)%hindcv(1)
          Xmax=OCEAN(ng)%hindcv(2)
          Ymin=OCEAN(ng)%hindcv(3)
          Ymax=OCEAN(ng)%hindcv(4)
          IF (spherical) THEN
            CALL linterp2d_ref (ng, Imin, Imax, Jmin, Jmax,             &
     &                          REFINED(ng)%lonv_c(Imin:Imax,Jmin:Jmax),&
     &                          REFINED(ng)%latv_c(Imin:Imax,Jmin:Jmax),&
     &                          A, Amask,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Xmin, Xmax, Ymin, Ymax,                 &
     &                          OCEAN(ng)%hindci_v, OCEAN(ng)%hindcj_v, &
     &                          GRID(ng)%lonv, GRID(ng)%latv, A2,       &
     &                          my_min, my_max)
          ELSE
            CALL linterp2d_ref (ng, Imin, Imax, Jmin, Jmax,             &
     &                          REFINED(ng)%xv_c(Imin:Imax,Jmin:Jmax),  &
     &                          REFINED(ng)%yv_c(Imin:Imax,Jmin:Jmax),  &
     &                          A, Amask,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Xmin, Xmax, Ymin, Ymax,                 &
     &                          OCEAN(ng)%hindci_v, OCEAN(ng)%hindcj_v, &
     &                          GRID(ng)%xv, GRID(ng)%yv, A2,           &
     &                          my_min, my_max)
          END IF
          DO j=Ymin,Ymax
            DO i=Xmin,Xmax
              my_val=A2(i,j)
# ifdef MASKING
              my_val=my_val*GRID(ng)%vmask(i,j)
# endif
              OCEAN(ng)%vbar(i,j,1)=my_val
              OCEAN(ng)%vbar(i,j,2)=my_val
              OCEAN(ng)%vbar(i,j,3)=my_val
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
#endif

#undef SOLITON_HACK
#ifdef SOLITON_HACK
!  Replace parent data with child data.
!  This is a hack for SOLITON test case
!
! zeta
        DO jp=Jstrc,Jendc
          DO ip=Istrc,Iendc
            DO j=Jstr,Jend
              DO i=Istr,Iend
                IF ((GRID(ngc)%xr(ip,jp).eq.GRID(ng)%xr(i,j)).and.      &
     &              (GRID(ngc)%yr(ip,jp).eq.GRID(ng)%yr(i,j))) THEN
!                    my_val=OCEAN(ng)%zeta(i,j,1)
!# ifdef MASKING
!                    my_val=my_val*GRID(ngp)%rmask(ip,jp)
!# endif
!                    OCEAN(ngp)%zeta(ip,jp,1)=my_val
!                    OCEAN(ngp)%zeta(ip,jp,2)=my_val
!                    OCEAN(ngp)%zeta(ip,jp,3)=my_val
!if call from the top
                    my_val=OCEAN(2)%zeta(ip,jp,1)
                    OCEAN(1)%zeta(i,j,1)=my_val
                    OCEAN(1)%zeta(i,j,2)=my_val
                    OCEAN(1)%zeta(i,j,3)=my_val
                END IF
              END DO
            END DO
          END DO
        END DO
! ubar
        DO jp=Jstrc,Jendc
          DO ip=Istrc+1,Iendc
            DO j=Jstr,Jend
              DO i=IstrU,Iend
                IF ((GRID(ngc)%xu(ip,jp).eq.GRID(ng)%xu(i,j)).and.      &
     &              (GRID(ngc)%yu(ip,jp).eq.GRID(ng)%yu(i,j))) THEN
!                    my_val=OCEAN(ng)%ubar(i,j,1)
!# ifdef MASKING
!                    my_val=my_val*GRID(ngp)%umask(ip,jp)
!# endif
!                    OCEAN(ngp)%ubar(ip,jp,1)=my_val
!                    OCEAN(ngp)%ubar(ip,jp,2)=my_val
!                    OCEAN(ngp)%ubar(ip,jp,3)=my_val
                    my_val=OCEAN(2)%ubar(ip,jp,1)
                    OCEAN(1)%ubar(i,j,1)=my_val
                    OCEAN(1)%ubar(i,j,2)=my_val
                    OCEAN(1)%ubar(i,j,3)=my_val
                END IF
              END DO
            END DO
          END DO
        END DO
! vbar
        DO jp=Jstrc+1,Jendc
          DO ip=Istrc,Iendc
            DO j=JstrV,Jend
              DO i=Istr,Iend
                IF ((GRID(ngc)%xv(ip,jp).eq.GRID(ng)%xv(i,j)).and.      &
     &              (GRID(ngc)%yv(ip,jp).eq.GRID(ng)%yv(i,j))) THEN
!                    my_val=OCEAN(ng)%vbar(i,j,1)
!# ifdef MASKING
!                    my_val=my_val*GRID(ngp)%vmask(ip,jp)
!# endif
!                    OCEAN(ngp)%vbar(ip,jp,1)=my_val
!                    OCEAN(ngp)%vbar(ip,jp,2)=my_val
!                    OCEAN(ngp)%vbar(ip,jp,3)=my_val
                     my_val=OCEAN(2)%vbar(ip,jp,1)
                    OCEAN(1)%vbar(i,j,1)=my_val
                    OCEAN(1)%vbar(i,j,2)=my_val
                    OCEAN(1)%vbar(i,j,3)=my_val
                END IF
              END DO
            END DO
          END DO
        END DO
#endif


# if defined EW_PERIODIC || defined NS_PERIODIC
!
!-----------------------------------------------------------------------
!  Exchange boundary information.
!-----------------------------------------------------------------------
!
      CALL exchange_u2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        OCEAN(ng)%zeta(:,:,1))
      CALL exchange_u2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        OCEAN(ng)%zeta(:,:,2))
      CALL exchange_u2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        OCEAN(ng)%zeta(:,:,3))
      CALL exchange_u2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        OCEAN(ng)%ubar(:,:,1))
      CALL exchange_v2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        OCEAN(ng)%vbar(:,:,1))
      CALL exchange_u2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        OCEAN(ng)%ubar(:,:,2))
      CALL exchange_v2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        OCEAN(ng)%vbar(:,:,2))
      CALL exchange_u2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        OCEAN(ng)%ubar(:,:,3))
      CALL exchange_v2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        OCEAN(ng)%vbar(:,:,3))
#  endif
#  ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, iNLM, 3,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    OCEAN(ng)%zeta(:,:,1),                        &
     &                    OCEAN(ng)%zeta(:,:,2),                        &
     &                    OCEAN(ng)%zeta(:,:,3))
      CALL mp_exchange2d (ng, tile, iNLM, 3,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    OCEAN(ng)%ubar(:,:,1),                        &
     &                    OCEAN(ng)%ubar(:,:,2),                        &
     &                    OCEAN(ng)%ubar(:,:,3))
      CALL mp_exchange2d (ng, tile, iNLM, 3,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    OCEAN(ng)%vbar(:,:,1),                        &
     &                    OCEAN(ng)%vbar(:,:,2),                        &
     &                    OCEAN(ng)%vbar(:,:,3))
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
      END SUBROUTINE set_2dchild_data_tile
#endif
      END MODULE set_2dchild_data_mod
