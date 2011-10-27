#include "cppdefs.h"
      MODULE init_child_hindices_mod
#ifdef REFINED_GRID
!
!svn $Id: set_data.F 1405 2008-02-20 01:44:07Z jcwarner $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  For the entire parent grid, this routine determines the             !
!  horizontal indices of the child grid                                !
!  at the locations of the parent grid.                                !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: init_child_hindices, init_child_hindices_tile

      CONTAINS
!
!***********************************************************************
      SUBROUTINE init_child_hindices (ng, tile)
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
# include "tile.h"
!
# ifdef PROFILE
      CALL wclock_on (ng, iNLM, 4)
# endif
      CALL init_child_hindices_tile (ng, tile,                                &
     &                    LBi, UBi, LBj, UBj)
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 4)
# endif
      RETURN
      END SUBROUTINE init_child_hindices
!
!***********************************************************************
      SUBROUTINE init_child_hindices_tile (ng, tile,                    &
     &                                     LBi, UBi, LBj, UBj)
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
      USE interpolate_mod
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

# ifdef OBC
      integer :: ILB, IUB, JLB, JUB
# endif
      integer :: LBip, UBip, LBjp, UBjp
      integer :: Imin, Imax, Jmin, Jmax, model
      integer :: i, ii, j, ngp, ngc, ic
      integer :: gtype, MyType, Ilen, Jlen, IJlen
      integer :: Npts

      real(r8) :: IJspv, my_min, my_max
      integer :: Xmin, Xmax, Ymin, Ymax, MyError

      logical :: rectangular
      real(r8), pointer :: Ipos(:,:)
      real(r8), pointer :: Jpos(:,:)
      real(r8), allocatable :: A(:)

      character (len=20) :: var_name

!
# include "set_bounds.h"
!
      ngc=mychild(ng)
      model=iNLM
!
!  Get the full grid of my parent, and compute the indices 
!  of the local tile in the parent grid.
!  Loop thru 3 times. For each loop we compute:
!  ii=1: rho var, ii=2:u var, ii=3:v var
!
      DO ii=1,3
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
        Npts=IJlen

        rectangular=.FALSE.
        IJspv=-9999.0_r8
!
! Get the grid information for the entire parent grid. 
! Use new routine "get_gridcoords.F".
!
        IF (ii.eq.1) THEN
!
! Now for rho points
!
          IF (.not.allocated(A)) THEN
            allocate (A(1:IJlen))
          END IF
          DO i=1,IJlen
            A(i)=0.0_r8
          END DO
!
! Case h
!
          var_name='h'
          CALL get_gridcoords (ngc, ncGRDid(ngc), Ilen, Jlen,           &
     &                         var_name, A)
          CALL mpi_bcast (A, IJlen, MP_FLOAT, MyMaster,                 &
     &                    OCN_COMM_WORLD, MyError)
          ic=0
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              ic=ic+1
              REFINED(ng)%h_c(i,j)=A(ic)
            END DO
          END DO
!
! Case angle
!
          var_name='angle'
          CALL get_gridcoords (ngc, ncGRDid(ngc), Ilen, Jlen,           &
     &                         var_name, A)
          CALL mpi_bcast (A, IJlen, MP_FLOAT, MyMaster,                 &
     &                    OCN_COMM_WORLD, MyError)
          ic=0
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              ic=ic+1
              REFINED(ng)%angler_c(i,j)=A(ic)
            END DO
          END DO
!
! Case lon_rho
!
          var_name='lon_rho'
          CALL get_gridcoords (ngc, ncGRDid(ngc), Ilen, Jlen,           &
     &                         var_name, A)
          CALL mpi_bcast (A, IJlen, MP_FLOAT, MyMaster,                 &
     &                    OCN_COMM_WORLD, MyError)
          ic=0
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              ic=ic+1
              REFINED(ng)%lonr_c(i,j)=A(ic)
            END DO
          END DO
!
! Case lat_rho
!
          var_name='lat_rho'
          CALL get_gridcoords (ngc, ncGRDid(ngc), Ilen, Jlen,           &
     &                         var_name, A)
          CALL mpi_bcast (A, IJlen, MP_FLOAT, MyMaster,                 &
     &                    OCN_COMM_WORLD, MyError)
          ic=0
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              ic=ic+1
              REFINED(ng)%latr_c(i,j)=A(ic)
            END DO
          END DO
!
! Case x_rho
!
          var_name='x_rho'
          CALL get_gridcoords (ngc, ncGRDid(ngc), Ilen, Jlen,           &
     &                         var_name, A)
          CALL mpi_bcast (A, IJlen, MP_FLOAT, MyMaster,                 &
     &                    OCN_COMM_WORLD, MyError)
          ic=0
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              ic=ic+1
              REFINED(ng)%xr_c(i,j)=A(ic)
            END DO
          END DO
!
! Case y_rho
!
          var_name='y_rho'
          CALL get_gridcoords (ngc, ncGRDid(ngc), Ilen, Jlen,           &
     &                         var_name, A)
          CALL mpi_bcast (A, IJlen, MP_FLOAT, MyMaster,                 &
     &                    OCN_COMM_WORLD, MyError)
          ic=0
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              ic=ic+1
              REFINED(ng)%yr_c(i,j)=A(ic)
            END DO
          END DO
# ifdef MASKING
!
! Case mask_rho
!
          var_name='mask_rho'
          CALL get_gridcoords (ngc, ncGRDid(ngc), Ilen, Jlen,           &
     &                         var_name, A)
          CALL mpi_bcast (A, IJlen, MP_FLOAT, MyMaster,                 &
     &                    OCN_COMM_WORLD, MyError)
          ic=0
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              ic=ic+1
              REFINED(ng)%rmask_c(i,j)=A(ic)
            END DO
          END DO
# endif
          IF (allocated(A)) THEN
            deallocate (A)
          END IF
        ELSE IF (ii.eq.2) THEN
!
! Now for u points
!
          IF (.not.allocated(A)) THEN
            allocate (A(1:IJlen))
          END IF
          DO i=1,IJlen
            A(i)=0.0_r8
          END DO
!
! Case lon_u
!
          var_name='lon_u'
          CALL get_gridcoords (ngc, ncGRDid(ngc), Ilen, Jlen,           &
     &                         var_name, A)
          CALL mpi_bcast (A, IJlen, MP_FLOAT, MyMaster,                 &
     &                    OCN_COMM_WORLD, MyError)
          ic=0
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              ic=ic+1
              REFINED(ng)%lonu_c(i,j)=A(ic)
            END DO
          END DO
!
! Case lat_u
!
          var_name='lat_u'
          CALL get_gridcoords (ngc, ncGRDid(ngc), Ilen, Jlen,           &
     &                         var_name, A)
          CALL mpi_bcast (A, IJlen, MP_FLOAT, MyMaster,                 &
     &                    OCN_COMM_WORLD, MyError)
          ic=0
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              ic=ic+1
              REFINED(ng)%latu_c(i,j)=A(ic)
            END DO
          END DO
!
! Case x_u
!
          var_name='x_u'
          CALL get_gridcoords (ngc, ncGRDid(ngc), Ilen, Jlen,           &
     &                         var_name, A)
          CALL mpi_bcast (A, IJlen, MP_FLOAT, MyMaster,                 &
     &                    OCN_COMM_WORLD, MyError)
          ic=0
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              ic=ic+1
              REFINED(ng)%xu_c(i,j)=A(ic)
            END DO
          END DO
!
! Case y_u
!
          var_name='y_u'
          CALL get_gridcoords (ngc, ncGRDid(ngc), Ilen, Jlen,           &
     &                         var_name, A)
          CALL mpi_bcast (A, IJlen, MP_FLOAT, MyMaster,                 &
     &                    OCN_COMM_WORLD, MyError)
          ic=0
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              ic=ic+1
              REFINED(ng)%yu_c(i,j)=A(ic)
            END DO
          END DO
# ifdef MASKING
!
! Case mask_u
!
          var_name='mask_u'
          CALL get_gridcoords (ngc, ncGRDid(ngc), Ilen, Jlen,           &
     &                         var_name, A)
          CALL mpi_bcast (A, IJlen, MP_FLOAT, MyMaster,                 &
     &                    OCN_COMM_WORLD, MyError)
          ic=0
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              ic=ic+1
              REFINED(ng)%umask_c(i,j)=A(ic)
            END DO
          END DO
# endif
          IF (allocated(A)) THEN
            deallocate (A)
          END IF
        ELSE IF (ii.eq.3) THEN
!
! Now for v points
!
          IF (.not.allocated(A)) THEN
            allocate (A(1:IJlen))
          END IF
          DO i=1,IJlen
            A(i)=0.0_r8
          END DO
!
! Case lon_v
!
          var_name='lon_v'
          CALL get_gridcoords (ngc, ncGRDid(ngc), Ilen, Jlen,           &
     &                         var_name, A)
          CALL mpi_bcast (A, IJlen, MP_FLOAT, MyMaster,                 &
     &                    OCN_COMM_WORLD, MyError)
          ic=0
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              ic=ic+1
              REFINED(ng)%lonv_c(i,j)=A(ic)
            END DO
          END DO
!
! Case lat_v
!
          var_name='lat_v'
          CALL get_gridcoords (ngc, ncGRDid(ngc), Ilen, Jlen,           &
     &                         var_name, A)
          CALL mpi_bcast (A, IJlen, MP_FLOAT, MyMaster,                 &
     &                    OCN_COMM_WORLD, MyError)
          ic=0
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              ic=ic+1
              REFINED(ng)%latv_c(i,j)=A(ic)
            END DO
          END DO
!
! Case x_v
!
          var_name='x_v'
          CALL get_gridcoords (ngc, ncGRDid(ngc), Ilen, Jlen,           &
     &                         var_name, A)
          CALL mpi_bcast (A, IJlen, MP_FLOAT, MyMaster,                 &
     &                    OCN_COMM_WORLD, MyError)
          ic=0
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              ic=ic+1
              REFINED(ng)%xv_c(i,j)=A(ic)
            END DO
          END DO
!
! Case y_v
!
          var_name='y_v'
          CALL get_gridcoords (ngc, ncGRDid(ngc), Ilen, Jlen,           &
     &                         var_name, A)
          CALL mpi_bcast (A, IJlen, MP_FLOAT, MyMaster,                 &
     &                    OCN_COMM_WORLD, MyError)
          ic=0
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              ic=ic+1
              REFINED(ng)%yv_c(i,j)=A(ic)
            END DO
          END DO
# ifdef MASKING
!
! Case mask_v
!
          var_name='mask_v'
          CALL get_gridcoords (ngc, ncGRDid(ngc), Ilen, Jlen,           &
     &                         var_name, A)
          CALL mpi_bcast (A, IJlen, MP_FLOAT, MyMaster,                 &
     &                    OCN_COMM_WORLD, MyError)
          ic=0
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              ic=ic+1
              REFINED(ng)%vmask_c(i,j)=A(ic)
            END DO
          END DO
# endif
          IF (allocated(A)) THEN
            deallocate (A)
          END IF
        END IF
!
! Get the horizontal indices of the grid locations in the parent
! grid for this local child tile.
!
! Set up a temp 2D array to fill the hindices for this tile.
!
        allocate ( Ipos(LBi:UBi,LBj:UBj) )
        allocate ( Jpos(LBi:UBi,LBj:UBj) )
!
        IF (ii.eq.1) THEN
          IF (spherical) THEN
            CALL hindices (ng, Imin, Imax, Jmin, Jmax,                  &
     &                   Imin, Imax, Jmin, Jmax,                        &
     &                   REFINED(ng)%angler_c(Imin:Imax,Jmin:Jmax),     &
     &                   REFINED(ng)%lonr_c(Imin:Imax,Jmin:Jmax),       &
     &                   REFINED(ng)%latr_c(Imin:Imax,Jmin:Jmax),       &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   GRID(ng)%lonr, GRID(ng)%latr,                  &
     &                   Ipos, Jpos,                                    &
     &                   IJspv, rectangular, Xmin, Xmax, Ymin, Ymax)
          ELSE
            CALL hindices (ng, Imin, Imax, Jmin, Jmax,                  &
     &                   Imin, Imax, Jmin, Jmax,                        &
     &                   REFINED(ng)%angler_c(Imin:Imax,Jmin:Jmax),     &
     &                   REFINED(ng)%xr_c(Imin:Imax,Jmin:Jmax),         &
     &                   REFINED(ng)%yr_c(Imin:Imax,Jmin:Jmax),         &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   GRID(ng)%xr, GRID(ng)%yr,                      &
     &                   Ipos, Jpos,                                    &
     &                   IJspv, rectangular, Xmin, Xmax, Ymin, Ymax)
          END IF
!
!  Recompute Xmin, Xmax, Ymin, and Ymax to determine the limits for 
!  child grid replacement into the parent.
!
          Xmin= 9999
          Xmax=-9999
          Ymin= 9999
          Ymax=-9999
          DO j=LBj,UBj
            DO i=LBi,UBi
              IF (Ipos(i,j).gt.1.0_r8) THEN
                Xmin=MIN(Xmin,i)
                IF (Ipos(i,j).lt.Lm(ngc)) THEN
                  Xmax=MAX(Xmax,i)
                END IF
              END IF
              IF (Jpos(i,j).gt.1.0_r8) THEN
                Ymin=MIN(Ymin,j)
                IF (Jpos(i,j).lt.Mm(ngc)) THEN
                  Ymax=MAX(Ymax,j)
                END IF
              END IF
            END DO
          END DO
          OCEAN(ng)%hindci_r=IPOS
          OCEAN(ng)%hindcj_r=JPOS
          OCEAN(ng)%hindcr(1)=Xmin
          OCEAN(ng)%hindcr(2)=Xmax
          OCEAN(ng)%hindcr(3)=Ymin
          OCEAN(ng)%hindcr(4)=Ymax
        ELSE IF (ii.eq.2) THEN
          IF (spherical) THEN
            CALL hindices (ng, Imin, Imax, Jmin, Jmax,                  &
     &                   Imin, Imax, Jmin, Jmax,                        &
     &                   REFINED(ng)%angler_c(Imin:Imax,Jmin:Jmax),     &
     &                   REFINED(ng)%lonu_c(Imin:Imax,Jmin:Jmax),       &
     &                   REFINED(ng)%latu_c(Imin:Imax,Jmin:Jmax),       &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   GRID(ng)%lonu, GRID(ng)%latu,                  &
     &                   Ipos, Jpos,                                    &
     &                   IJspv, rectangular, Xmin, Xmax, Ymin, Ymax)
          ELSE
            CALL hindices (ng, Imin, Imax, Jmin, Jmax,                  &
     &                   Imin, Imax, Jmin, Jmax,                        &
     &                   REFINED(ng)%angler_c(Imin:Imax,Jmin:Jmax),     &
     &                   REFINED(ng)%xu_c(Imin:Imax,Jmin:Jmax),         &
     &                   REFINED(ng)%yu_c(Imin:Imax,Jmin:Jmax),         &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   GRID(ng)%xu, GRID(ng)%yu,                      &
     &                   Ipos, Jpos,                                    &
     &                   IJspv, rectangular, Xmin, Xmax, Ymin, Ymax)
          END IF
!
!  Recompute Xmin, Xmax, Ymin, and Ymax to determine the limits for 
!  child grid replacement into the parent.
!
          Xmin= 9999
          Xmax=-9999
          Ymin= 9999
          Ymax=-9999
          DO j=LBj,UBj
            DO i=LBi,UBi
              IF (Ipos(i,j).gt.1.0_r8) THEN
                Xmin=MIN(Xmin,i)
                IF (Ipos(i,j).lt.Lm(ngc)) THEN
                  Xmax=MAX(Xmax,i)
                END IF
              END IF
              IF (Jpos(i,j).gt.1.0_r8) THEN
                Ymin=MIN(Ymin,j)
                IF (Jpos(i,j).lt.Mm(ngc)) THEN
                  Ymax=MAX(Ymax,j)
                END IF
              END IF
            END DO
          END DO
          OCEAN(ng)%hindci_u=IPOS
          OCEAN(ng)%hindcj_u=JPOS
          OCEAN(ng)%hindcu(1)=Xmin
          OCEAN(ng)%hindcu(2)=Xmax
          OCEAN(ng)%hindcu(3)=Ymin
          OCEAN(ng)%hindcu(4)=Ymax
        ELSE IF (ii.eq.3) THEN
          IF (spherical) THEN
            CALL hindices (ng, Imin, Imax, Jmin, Jmax,                  &
     &                   Imin, Imax, Jmin, Jmax,                        &
     &                   REFINED(ng)%angler_c(Imin:Imax,Jmin:Jmax),     &
     &                   REFINED(ng)%lonv_c(Imin:Imax,Jmin:Jmax),       &
     &                   REFINED(ng)%latv_c(Imin:Imax,Jmin:Jmax),       &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   GRID(ng)%lonv, GRID(ng)%latv,                  &
     &                   Ipos, Jpos,                                    &
     &                   IJspv, rectangular, Xmin, Xmax, Ymin, Ymax)
          ELSE
            CALL hindices (ng, Imin, Imax, Jmin, Jmax,                  &
     &                   Imin, Imax, Jmin, Jmax,                        &
     &                   REFINED(ng)%angler_c(Imin:Imax,Jmin:Jmax),     &
     &                   REFINED(ng)%xv_c(Imin:Imax,Jmin:Jmax),         &
     &                   REFINED(ng)%yv_c(Imin:Imax,Jmin:Jmax),         &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   GRID(ng)%xv, GRID(ng)%yv,                      &
     &                   Ipos, Jpos,                                    &
     &                   IJspv, rectangular, Xmin, Xmax, Ymin, Ymax)
          END IF
!
!  Recompute Xmin, Xmax, Ymin, and Ymax to determine the limits for 
!  child grid replacement into the parent.
!
          Xmin= 9999
          Xmax=-9999
          Ymin= 9999
          Ymax=-9999
          DO j=LBj,UBj
            DO i=LBi,UBi
              IF (Ipos(i,j).gt.1.0_r8) THEN
                Xmin=MIN(Xmin,i)
                IF (Ipos(i,j).lt.Lm(ngc)) THEN
                  Xmax=MAX(Xmax,i)
                END IF
              END IF
              IF (Jpos(i,j).gt.1.0_r8) THEN
                Ymin=MIN(Ymin,j)
                IF (Jpos(i,j).lt.Mm(ngc)) THEN
                  Ymax=MAX(Ymax,j)
                END IF
              END IF
            END DO
          END DO
          OCEAN(ng)%hindci_v=IPOS
          OCEAN(ng)%hindcj_v=JPOS
          OCEAN(ng)%hindcv(1)=Xmin
          OCEAN(ng)%hindcv(2)=Xmax
          OCEAN(ng)%hindcv(3)=Ymin
          OCEAN(ng)%hindcv(4)=Ymax
        END IF
!
        IF (ii.eq.1) THEN
!
!  Compute refined ratio for the child grids.
!
          CALL mpi_allreduce (OCEAN(ng)%hindcr(1), OCEAN(ng)%hindxmin,  &
     &                        1, MPI_INTEGER, MPI_MIN, OCN_COMM_WORLD,  &
     &                        MyError)
          CALL mpi_allreduce (OCEAN(ng)%hindcr(2), OCEAN(ng)%hindxmax,  &
     &                        1, MPI_INTEGER, MPI_MAX, OCN_COMM_WORLD,  &
     &                        MyError)
          REFINED(ngc)%ref_scale=(Ilen-7)/                              &
                              (OCEAN(ng)%hindxmax-OCEAN(ng)%hindxmin+1)
        END IF
!
        deallocate ( Ipos, Jpos)
      END DO

      RETURN
      END SUBROUTINE init_child_hindices_tile
#endif
      END MODULE init_child_hindices_mod
