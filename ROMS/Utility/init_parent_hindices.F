#include "cppdefs.h"
      MODULE init_parent_hindices_mod
#ifdef REFINED_GRID
!
!svn $Id: set_data.F 1405 2008-02-20 01:44:07Z jcwarner $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  For the entire child grid, this routine determines the              !
!  horizontal indices of the parent grid                               !
!  at the locations of the child grid.                                 !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: init_parent_hindices, init_parent_hindices_tile

      CONTAINS
!
!***********************************************************************
      SUBROUTINE init_parent_hindices (ng, tile)
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
      CALL init_parent_hindices_tile (ng, tile,                                &
     &                    LBi, UBi, LBj, UBj)
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 4)
# endif
      RETURN
      END SUBROUTINE init_parent_hindices
!
!***********************************************************************
      SUBROUTINE init_parent_hindices_tile (ng, tile,                          &
     &                                      LBi, UBi, LBj, UBj)
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
      ngp=myparent(ng)
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
        IF (ngp.eq.1) THEN
          Imin=0
          Imax=Lm(ngp)+1
          Jmin=0
          Jmax=Mm(ngp)+1
          MyType=gtype
          IF ((ABS(MyType).eq.p2dvar).or.(ABS(MyType).eq.u2dvar)) THEN
            Imin=1
          END IF
          IF ((ABS(MyType).eq.p2dvar).or.(ABS(MyType).eq.v2dvar)) THEN
            Jmin=1
          END IF
        ELSE
          Imin=LOWER_BOUND_I
          Imax=Lm(ngp)+GHOST_POINTS
          Jmin=LOWER_BOUND_J
          Jmax=Mm(ngp)+GHOST_POINTS
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
          CALL get_gridcoords (ngp, ncGRDid(ngp), Ilen, Jlen,           &
     &                         var_name, A)
          CALL mpi_bcast (A, IJlen, MP_FLOAT, MyMaster,                 &
     &                    OCN_COMM_WORLD, MyError)
          ic=0
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              ic=ic+1
              REFINED(ng)%h_p(i,j)=A(ic)
            END DO
          END DO
!
! Case angle
!
          var_name='angle'
          CALL get_gridcoords (ngp, ncGRDid(ngp), Ilen, Jlen,           &
     &                         var_name, A)
          CALL mpi_bcast (A, IJlen, MP_FLOAT, MyMaster,                 &
     &                    OCN_COMM_WORLD, MyError)
          ic=0
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              ic=ic+1
              REFINED(ng)%angler_p(i,j)=A(ic)
            END DO
          END DO
!
! Case lon_rho
!
          var_name='lon_rho'
          CALL get_gridcoords (ngp, ncGRDid(ngp), Ilen, Jlen,           &
     &                         var_name, A)
          CALL mpi_bcast (A, IJlen, MP_FLOAT, MyMaster,                 &
     &                    OCN_COMM_WORLD, MyError)
          ic=0
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              ic=ic+1
              REFINED(ng)%lonr_p(i,j)=A(ic)
            END DO
          END DO
!
! Case lat_rho
!
          var_name='lat_rho'
          CALL get_gridcoords (ngp, ncGRDid(ngp), Ilen, Jlen,           &
     &                         var_name, A)
          CALL mpi_bcast (A, IJlen, MP_FLOAT, MyMaster,                 &
     &                    OCN_COMM_WORLD, MyError)
          ic=0
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              ic=ic+1
              REFINED(ng)%latr_p(i,j)=A(ic)
            END DO
          END DO
!
! Case x_rho
!
          var_name='x_rho'
          CALL get_gridcoords (ngp, ncGRDid(ngp), Ilen, Jlen,           &
     &                         var_name, A)
          CALL mpi_bcast (A, IJlen, MP_FLOAT, MyMaster,                 &
     &                    OCN_COMM_WORLD, MyError)
          ic=0
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              ic=ic+1
              REFINED(ng)%xr_p(i,j)=A(ic)
            END DO
          END DO
!
! Case y_rho
!
          var_name='y_rho'
          CALL get_gridcoords (ngp, ncGRDid(ngp), Ilen, Jlen,           &
     &                         var_name, A)
          CALL mpi_bcast (A, IJlen, MP_FLOAT, MyMaster,                 &
     &                    OCN_COMM_WORLD, MyError)
          ic=0
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              ic=ic+1
              REFINED(ng)%yr_p(i,j)=A(ic)
            END DO
          END DO
# ifdef MASKING
!
! Case mask_rho
!
          var_name='mask_rho'
          CALL get_gridcoords (ngp, ncGRDid(ngp), Ilen, Jlen,           &
     &                         var_name, A)
          CALL mpi_bcast (A, IJlen, MP_FLOAT, MyMaster,                 &
     &                    OCN_COMM_WORLD, MyError)
          ic=0
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              ic=ic+1
              REFINED(ng)%rmask_p(i,j)=A(ic)
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
          CALL get_gridcoords (ngp, ncGRDid(ngp), Ilen, Jlen,           &
     &                         var_name, A)
          CALL mpi_bcast (A, IJlen, MP_FLOAT, MyMaster,                 &
     &                    OCN_COMM_WORLD, MyError)
          ic=0
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              ic=ic+1
              REFINED(ng)%lonu_p(i,j)=A(ic)
            END DO
          END DO
!
! Case lat_u
!
          var_name='lat_u'
          CALL get_gridcoords (ngp, ncGRDid(ngp), Ilen, Jlen,           &
     &                         var_name, A)
          CALL mpi_bcast (A, IJlen, MP_FLOAT, MyMaster,                 &
     &                    OCN_COMM_WORLD, MyError)
          ic=0
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              ic=ic+1
              REFINED(ng)%latu_p(i,j)=A(ic)
            END DO
          END DO
!
! Case x_u
!
          var_name='x_u'
          CALL get_gridcoords (ngp, ncGRDid(ngp), Ilen, Jlen,           &
     &                         var_name, A)
          CALL mpi_bcast (A, IJlen, MP_FLOAT, MyMaster,                 &
     &                    OCN_COMM_WORLD, MyError)
          ic=0
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              ic=ic+1
              REFINED(ng)%xu_p(i,j)=A(ic)
            END DO
          END DO
!
! Case y_u
!
          var_name='y_u'
          CALL get_gridcoords (ngp, ncGRDid(ngp), Ilen, Jlen,           &
     &                         var_name, A)
          CALL mpi_bcast (A, IJlen, MP_FLOAT, MyMaster,                 &
     &                    OCN_COMM_WORLD, MyError)
          ic=0
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              ic=ic+1
              REFINED(ng)%yu_p(i,j)=A(ic)
            END DO
          END DO
# ifdef MASKING
!
! Case mask_u
!
          var_name='mask_u'
          CALL get_gridcoords (ngp, ncGRDid(ngp), Ilen, Jlen,           &
     &                         var_name, A)
          CALL mpi_bcast (A, IJlen, MP_FLOAT, MyMaster,                 &
     &                    OCN_COMM_WORLD, MyError)
          ic=0
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              ic=ic+1
              REFINED(ng)%umask_p(i,j)=A(ic)
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
          CALL get_gridcoords (ngp, ncGRDid(ngp), Ilen, Jlen,           &
     &                         var_name, A)
          CALL mpi_bcast (A, IJlen, MP_FLOAT, MyMaster,                 &
     &                    OCN_COMM_WORLD, MyError)
          ic=0
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              ic=ic+1
              REFINED(ng)%lonv_p(i,j)=A(ic)
            END DO
          END DO
!
! Case lat_v
!
          var_name='lat_v'
          CALL get_gridcoords (ngp, ncGRDid(ngp), Ilen, Jlen,           &
     &                         var_name, A)
          CALL mpi_bcast (A, IJlen, MP_FLOAT, MyMaster,                 &
     &                    OCN_COMM_WORLD, MyError)
          ic=0
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              ic=ic+1
              REFINED(ng)%latv_p(i,j)=A(ic)
            END DO
          END DO
!
! Case x_v
!
          var_name='x_v'
          CALL get_gridcoords (ngp, ncGRDid(ngp), Ilen, Jlen,           &
     &                         var_name, A)
          CALL mpi_bcast (A, IJlen, MP_FLOAT, MyMaster,                 &
     &                    OCN_COMM_WORLD, MyError)
          ic=0
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              ic=ic+1
              REFINED(ng)%xv_p(i,j)=A(ic)
            END DO
          END DO
!
! Case y_v
!
          var_name='y_v'
          CALL get_gridcoords (ngp, ncGRDid(ngp), Ilen, Jlen,           &
     &                         var_name, A)
          CALL mpi_bcast (A, IJlen, MP_FLOAT, MyMaster,                 &
     &                    OCN_COMM_WORLD, MyError)
          ic=0
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              ic=ic+1
              REFINED(ng)%yv_p(i,j)=A(ic)
            END DO
          END DO
# ifdef MASKING
!
! Case mask_v
!
          var_name='mask_v'
          CALL get_gridcoords (ngp, ncGRDid(ngp), Ilen, Jlen,           &
     &                         var_name, A)
          CALL mpi_bcast (A, IJlen, MP_FLOAT, MyMaster,                 &
     &                    OCN_COMM_WORLD, MyError)
          ic=0
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              ic=ic+1
              REFINED(ng)%vmask_p(i,j)=A(ic)
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
     &                   REFINED(ng)%angler_p(Imin:Imax,Jmin:Jmax),     &
     &                   REFINED(ng)%lonr_p(Imin:Imax,Jmin:Jmax),       &
     &                   REFINED(ng)%latr_p(Imin:Imax,Jmin:Jmax),       &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   GRID(ng)%lonr, GRID(ng)%latr,                  &
     &                   Ipos, Jpos,                                    &
     &                   IJspv, rectangular, Xmin, Xmax, Ymin, Ymax)
          ELSE
            CALL hindices (ng, Imin, Imax, Jmin, Jmax,                  &
     &                   Imin, Imax, Jmin, Jmax,                        &
     &                   REFINED(ng)%angler_p(Imin:Imax,Jmin:Jmax),     &
     &                   REFINED(ng)%xr_p(Imin:Imax,Jmin:Jmax),         &
     &                   REFINED(ng)%yr_p(Imin:Imax,Jmin:Jmax),         &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   GRID(ng)%xr, GRID(ng)%yr,                      &
     &                   Ipos, Jpos,                                    &
     &                   IJspv, rectangular, Xmin, Xmax, Ymin, Ymax)
          END IF
          OCEAN(ng)%hindpi_r=IPOS
          OCEAN(ng)%hindpj_r=JPOS
          OCEAN(ng)%hindpr(1)=Xmin
          OCEAN(ng)%hindpr(2)=Xmax
          OCEAN(ng)%hindpr(3)=Ymin
          OCEAN(ng)%hindpr(4)=Ymax
        ELSE IF (ii.eq.2) THEN
          IF (spherical) THEN
            CALL hindices (ng, Imin, Imax, Jmin, Jmax,                  &
     &                   Imin, Imax, Jmin, Jmax,                        &
     &                   REFINED(ng)%angler_p(Imin:Imax,Jmin:Jmax),     &
     &                   REFINED(ng)%lonu_p(Imin:Imax,Jmin:Jmax),       &
     &                   REFINED(ng)%latu_p(Imin:Imax,Jmin:Jmax),       &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   GRID(ng)%lonu, GRID(ng)%latu,                  &
     &                   Ipos, Jpos,                                    &
     &                   IJspv, rectangular, Xmin, Xmax, Ymin, Ymax)
          ELSE
            CALL hindices (ng, Imin, Imax, Jmin, Jmax,                  &
     &                   Imin, Imax, Jmin, Jmax,                        &
     &                   REFINED(ng)%angler_p(Imin:Imax,Jmin:Jmax),     &
     &                   REFINED(ng)%xu_p(Imin:Imax,Jmin:Jmax),         &
     &                   REFINED(ng)%yu_p(Imin:Imax,Jmin:Jmax),         &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   GRID(ng)%xu, GRID(ng)%yu,                      &
     &                   Ipos, Jpos,                                    &
     &                   IJspv, rectangular, Xmin, Xmax, Ymin, Ymax)
          END IF
          OCEAN(ng)%hindpi_u=IPOS
          OCEAN(ng)%hindpj_u=JPOS
          OCEAN(ng)%hindpu(1)=Xmin
          OCEAN(ng)%hindpu(2)=Xmax
          OCEAN(ng)%hindpu(3)=Ymin
          OCEAN(ng)%hindpu(4)=Ymax
        ELSE IF (ii.eq.3) THEN
          IF (spherical) THEN
            CALL hindices (ng, Imin, Imax, Jmin, Jmax,                  &
     &                   Imin, Imax, Jmin, Jmax,                        &
     &                   REFINED(ng)%angler_p(Imin:Imax,Jmin:Jmax),     &
     &                   REFINED(ng)%lonv_p(Imin:Imax,Jmin:Jmax),       &
     &                   REFINED(ng)%latv_p(Imin:Imax,Jmin:Jmax),       &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   GRID(ng)%lonv, GRID(ng)%latv,                  &
     &                   Ipos, Jpos,                                    &
     &                   IJspv, rectangular, Xmin, Xmax, Ymin, Ymax)
          ELSE
            CALL hindices (ng, Imin, Imax, Jmin, Jmax,                  &
     &                   Imin, Imax, Jmin, Jmax,                        &
     &                   REFINED(ng)%angler_p(Imin:Imax,Jmin:Jmax),     &
     &                   REFINED(ng)%xv_p(Imin:Imax,Jmin:Jmax),         &
     &                   REFINED(ng)%yv_p(Imin:Imax,Jmin:Jmax),         &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   GRID(ng)%xv, GRID(ng)%yv,                      &
     &                   Ipos, Jpos,                                    &
     &                   IJspv, rectangular, Xmin, Xmax, Ymin, Ymax)
          END IF
          OCEAN(ng)%hindpi_v=IPOS
          OCEAN(ng)%hindpj_v=JPOS
          OCEAN(ng)%hindpv(1)=Xmin
          OCEAN(ng)%hindpv(2)=Xmax
          OCEAN(ng)%hindpv(3)=Ymin
          OCEAN(ng)%hindpv(4)=Ymax
        END IF
!
        deallocate ( Ipos, Jpos )
      END DO
      RETURN
      END SUBROUTINE init_parent_hindices_tile
#endif
      END MODULE init_parent_hindices_mod
