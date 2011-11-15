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
# ifdef DISTRIBUTE
      USE distribute_mod, ONLY : mp_bcasti, mp_gather2d
# endif
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
      integer :: Iminc, Imaxc, Jminc, Jmaxc

      real(r8) :: IJspv, my_min, my_max, offset
      integer :: Xmin, Xmax, Ymin, Ymax, MyError, MySizer

      logical :: rectangular
      real(r8), allocatable :: A(:)
      real(r8), allocatable :: Ipos(:,:)
      real(r8), allocatable :: Jpos(:,:)
      real(r8), allocatable :: lonpsi_par(:,:)
      real(r8), allocatable :: latpsi_par(:,:)

      character (len=20) :: var_name
!
# include "set_bounds.h"
!
      ngc=mychild(ng)
      model=iNLM
!
!  Loop thru 4 times. For each loop we compute:
!  ii=1: rho var, ii=2:u var, ii=3:v var, ii=4:psi var
!
      DO ii=1,4
        IF (ii.eq.1) THEN
          gtype=r2dvar
        ELSE IF (ii.eq.2) THEN
          gtype=u2dvar
        ELSE IF (ii.eq.3) THEN
          gtype=v2dvar
        ELSE IF (ii.eq.4) THEN
          gtype=p2dvar
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
        ELSE IF (ii.eq.4) THEN
!
! Now for psi points
!
          IF (.not.allocated(A)) THEN
            allocate (A(1:IJlen))
          END IF
          DO i=1,IJlen
            A(i)=0.0_r8
          END DO
!
! Case lon_psi
!
          var_name='lon_psi'
          CALL get_gridcoords (ngc, ncGRDid(ngc), Ilen, Jlen,           &
     &                         var_name, A)
          CALL mpi_bcast (A, IJlen, MP_FLOAT, MyMaster,                 &
     &                    OCN_COMM_WORLD, MyError)
          ic=0
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              ic=ic+1
              REFINED(ng)%lonp_c(i,j)=A(ic)
            END DO
          END DO
!
! Case lat_psi
!
          var_name='lat_psi'
          CALL get_gridcoords (ngc, ncGRDid(ngc), Ilen, Jlen,           &
     &                         var_name, A)
          CALL mpi_bcast (A, IJlen, MP_FLOAT, MyMaster,                 &
     &                    OCN_COMM_WORLD, MyError)
          ic=0
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              ic=ic+1
              REFINED(ng)%latp_c(i,j)=A(ic)
            END DO
          END DO
!
! Case x_psi
!
          var_name='x_psi'
          CALL get_gridcoords (ngc, ncGRDid(ngc), Ilen, Jlen,           &
     &                         var_name, A)
          CALL mpi_bcast (A, IJlen, MP_FLOAT, MyMaster,                 &
     &                    OCN_COMM_WORLD, MyError)
          ic=0
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              ic=ic+1
              REFINED(ng)%xp_c(i,j)=A(ic)
            END DO
          END DO
!
! Case y_psi
!
          var_name='y_psi'
          CALL get_gridcoords (ngc, ncGRDid(ngc), Ilen, Jlen,           &
     &                         var_name, A)
          CALL mpi_bcast (A, IJlen, MP_FLOAT, MyMaster,                 &
     &                    OCN_COMM_WORLD, MyError)
          ic=0
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              ic=ic+1
              REFINED(ng)%yp_c(i,j)=A(ic)
            END DO
          END DO
# ifdef MASKING
!
! Case mask_psi
!
          var_name='mask_psi'
          CALL get_gridcoords (ngc, ncGRDid(ngc), Ilen, Jlen,           &
     &                         var_name, A)
          CALL mpi_bcast (A, IJlen, MP_FLOAT, MyMaster,                 &
     &                    OCN_COMM_WORLD, MyError)
          ic=0
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              ic=ic+1
              REFINED(ng)%pmask_c(i,j)=A(ic)
            END DO
          END DO
# endif
          IF (allocated(A)) THEN
            deallocate (A)
          END IF
        END IF
      END DO
!
! Get the horizontal indices of the grid locations in the parent grid.
! Base this on the psi points to start.
! Gather all the lon and lat psi points of the parent and child.
! Use hindices to determine the location of the 4 corners of the child
! grid. Hindices is not accurate enough to give exact locations, so 
! we use integers increments to set all the hindices.
!
      gtype=p2dvar
      IF (ng.eq.1) THEN
        Imin=0
        Imax=Lm(ng)+1
        Jmin=0
        Jmax=Mm(ng)+1
        IF ((ABS(gtype).eq.p2dvar).or.(ABS(gtype).eq.u2dvar)) THEN
          Imin=1
        END IF
        IF ((ABS(gtype).eq.p2dvar).or.(ABS(gtype).eq.v2dvar)) THEN
          Jmin=1
        END IF
      ELSE
        Imin=LOWER_BOUND_I
        Imax=Lm(ng)+GHOST_POINTS
        Jmin=LOWER_BOUND_J
        Jmax=Mm(ng)+GHOST_POINTS
        IF ((ABS(gtype).eq.p2dvar).or.(ABS(gtype).eq.u2dvar)) THEN
          Imin=LOWER_BOUND_I+1
        END IF
        IF ((ABS(gtype).eq.p2dvar).or.(ABS(gtype).eq.v2dvar)) THEN
          Jmin=LOWER_BOUND_J+1
        END IF
      END IF
      Ilen=Imax-Imin+1
      Jlen=Jmax-Jmin+1
      IJlen=Ilen*Jlen
      MySizer=IJlen
      Npts=IJlen
      allocate (A(1:IJlen))
      allocate (lonpsi_par(Imin:Imax,Jmin:Jmax))
      allocate (latpsi_par(Imin:Imax,Jmin:Jmax))
      DO i=1,IJlen
        A(i)=0.0_r8
      END DO
      DO i=Imin,Imax
        DO j=Jmin,Jmax
          lonpsi_par(i,j)=0.0_r8
          latpsi_par(i,j)=0.0_r8
        END DO
      END DO
!
! Get all the parent grid lon lat psi data.
!
      IF (spherical) THEN
        CALL mp_gather2d (ng, model, LBi, UBi, LBj, UBj, 0,             &
     &                    gtype, 1.0_r8,                                &
# ifdef MASKING
     &                    GRID(ng) % pmask(LBi:UBi,LBj:UBj),            &
# endif
     &                    GRID(ng) % lonp(LBi:UBi,LBj:UBj),             &
     &                    Npts, A)
      ELSE
        CALL mp_gather2d (ng, model, LBi, UBi, LBj, UBj, 0,             &
     &                    gtype, 1.0_r8,                                &
# ifdef MASKING
     &                    GRID(ng) % pmask(LBi:UBi,LBj:UBj),            &
# endif
     &                    GRID(ng) % xp(LBi:UBi,LBj:UBj),               &
     &                    Npts, A)
      END IF
!
! distribute the A (1D packed) field to all the nodes.
!
      CALL mpi_bcast (A, MySizer, MP_FLOAT, MyMaster,                &
     &                OCN_COMM_WORLD, MyError)
!
! unpack and place in global ng=1 2D array
!
      ic=0
      DO j=Jmin,Jmax
        DO i=Imin,Imax
          ic=ic+1
          lonpsi_par(i,j)=A(ic)
        END DO
      END DO
      IF (spherical) THEN
        CALL mp_gather2d (ng, model, LBi, UBi, LBj, UBj, 0,             &
     &                    gtype, 1.0_r8,                                &
# ifdef MASKING
     &                    GRID(ng) % pmask(LBi:UBi,LBj:UBj),            &
# endif
     &                    GRID(ng) % latp(LBi:UBi,LBj:UBj),             &
     &                    Npts, A)
      ELSE
        CALL mp_gather2d (ng, model, LBi, UBi, LBj, UBj, 0,             &
     &                    gtype, 1.0_r8,                                &
# ifdef MASKING
     &                    GRID(ng) % pmask(LBi:UBi,LBj:UBj),            &
# endif
     &                    GRID(ng) % yp(LBi:UBi,LBj:UBj),               &
     &                    Npts, A)
      END IF
!
! distribute the A (1D packed) field to all the nodes.
!
      CALL mpi_bcast (A, MySizer, MP_FLOAT, MyMaster,                &
     &                OCN_COMM_WORLD, MyError)
!
! unpack and place in global ng=1 2D array
!
      ic=0
      DO j=Jmin,Jmax
        DO i=Imin,Imax
          ic=ic+1
          latpsi_par(i,j)=A(ic)
        END DO
      END DO
!
!  Now use hindices to get close enough for the child corners. then
!  int these values.
!
      Iminc=1                    !LOWER_BOUND_I+1
      Imaxc=Lm(ngc)+1            !GHOST_POINTS
      Jminc=1                    !LOWER_BOUND_J+1
      Jmaxc=Mm(ngc)+1            !GHOST_POINTS
!
      allocate ( Ipos(Imin:Imax,Jmin:Jmax) )
      allocate ( Jpos(Imin:Imax,Jmin:Jmax) )
!
      IF (spherical) THEN
        CALL hindices (ng, Iminc, Imaxc, Jminc, Jmaxc,                  &
     &                 Iminc, Imaxc, Jminc, Jmaxc,                      &
     &                 REFINED(ng)%angler_c(Iminc:Imaxc,Jminc:Jmaxc),   &
     &                 REFINED(ng)%lonp_c(Iminc:Imaxc,Jminc:Jmaxc),     &
     &                 REFINED(ng)%latp_c(Iminc:Imaxc,Jminc:Jmaxc),     &
     &                 Imin, Imax, Jmin, Jmax,                          &
     &                 Imin, Imax, Jmin, Jmax,                          &
     &                 lonpsi_par, latpsi_par,                          &
     &                 Ipos, Jpos,                                      &
     &                 IJspv, rectangular, Xmin, Xmax, Ymin, Ymax)
      ELSE
        CALL hindices (ng, Iminc, Imaxc, Jminc, Jmaxc,                  &
     &                 Iminc, Imaxc, Jminc, Jmaxc,                      &
     &                 REFINED(ng)%angler_c(Iminc:Imaxc,Jminc:Jmaxc),   &
     &                 REFINED(ng)%xp_c(Iminc:Imaxc,Jminc:Jmaxc),       &
     &                 REFINED(ng)%yp_c(Iminc:Imaxc,Jminc:Jmaxc),       &
     &                 Imin, Imax, Jmin, Jmax,                          &
     &                 Imin, Imax, Jmin, Jmax,                          &
     &                 lonpsi_par, latpsi_par,                          &
     &                 Ipos, Jpos,                                      &
     &                 IJspv, rectangular, Xmin, Xmax, Ymin, Ymax)
      END IF
!
!  Xmin Xmax Ymin Ymax are the bounding psi point limits 
!  of the parent grid that identify the boundary of the child grid.
!  Still working on the global parent grid (ie non-tiled).
!  The IPOS and JPOS are values of the child indices.
!  The XYmin maxs are the limits of the parent grid.
!
!  Compute refined ratio for the child grids.
!
      REFINED(ngc)%ref_scale=Lm(ngc)/(Xmax-Xmin)
!
      DO j=Jmin,Jmax
        DO i=Imin,Imax
          IPOS(i,j)=-9999.0_r8
          JPOS(i,j)=-9999.0_r8
        END DO
      END DO
      DO j=Ymin,Ymax
        DO i=Xmin,Xmax
          IPOS(i,j)=1.0_r8+REAL((i-Xmin)*REFINED(ngc)%ref_scale,r8)
          JPOS(i,j)=1.0_r8+REAL((j-Ymin)*REFINED(ngc)%ref_scale,r8)
        END DO
      END DO
      offset=(REAL(REFINED(ngc)%ref_scale,r8)-1.0_r8)/2.0_r8
!
!  Now do local rho points.
!
      OCEAN(ng)%hindcr(1)= 9999
      OCEAN(ng)%hindcr(2)=-9999
      OCEAN(ng)%hindcr(3)= 9999
      OCEAN(ng)%hindcr(4)=-9999
      DO j=Jstr,Jend
        DO i=Istr,Iend
          IF (Ipos(i,j).ge.1.0_r8) THEN
            OCEAN(ng)%hindcr(1)=MIN(OCEAN(ng)%hindcr(1),i)
            IF (Ipos(i,j).lt.REAL(Lm(ngc)+1,r8)) THEN
              OCEAN(ng)%hindcr(2)=MAX(OCEAN(ng)%hindcr(2),i)
            END IF
          END IF
          IF (Jpos(i,j).ge.1.0_r8) THEN
            OCEAN(ng)%hindcr(3)=MIN(OCEAN(ng)%hindcr(3),j)
            IF (Jpos(i,j).lt.REAL(Mm(ngc)+1,r8)) THEN
              OCEAN(ng)%hindcr(4)=MAX(OCEAN(ng)%hindcr(4),j)
            END IF
          END IF
        END DO
      END DO
      DO j=Jstr,Jend
        DO i=Istr,Iend
          OCEAN(ng)%hindci_r(i,j)=IJspv
          OCEAN(ng)%hindcj_r(i,j)=IJspv
          IF (((i.ge.Xmin).and.(i.le.Xmax-1)).and.                      &
              ((j.ge.Ymin).and.(j.le.Ymax-1)))  THEN
                OCEAN(ng)%hindci_r(i,j)=IPOS(i,j)+offset
                OCEAN(ng)%hindcj_r(i,j)=JPOS(i,j)+offset
          END IF
        END DO
      END DO
!
!  Now do local u points.
!
      OCEAN(ng)%hindcu(1)= 9999
      OCEAN(ng)%hindcu(2)=-9999
      OCEAN(ng)%hindcu(3)= 9999
      OCEAN(ng)%hindcu(4)=-9999
      DO j=Jstr,Jend
        DO i=Istr,Iend
          IF (Ipos(i,j).gt.1.0_r8) THEN
            OCEAN(ng)%hindcu(1)=MIN(OCEAN(ng)%hindcu(1),i)
            IF (Ipos(i,j).lt.REAL(Lm(ngc)+1,r8)) THEN
              OCEAN(ng)%hindcu(2)=MAX(OCEAN(ng)%hindcu(2),i)
            END IF
          END IF
          IF (Jpos(i,j).ge.1.0_r8) THEN
            OCEAN(ng)%hindcu(3)=MIN(OCEAN(ng)%hindcu(3),j)
            IF (Jpos(i,j).lt.REAL(Mm(ngc),r8)) THEN
              OCEAN(ng)%hindcu(4)=MAX(OCEAN(ng)%hindcu(4),j)
            END IF
          END IF
        END DO
      END DO
      DO j=Jstr,Jend
        DO i=Istr,Iend+1
          OCEAN(ng)%hindci_u(i,j)=IJspv
          OCEAN(ng)%hindcj_u(i,j)=IJspv
          IF (((i.ge.Xmin).and.(i.le.Xmax)).and.                      &
              ((j.ge.Ymin).and.(j.le.Ymax-1)))  THEN
                OCEAN(ng)%hindci_u(i,j)=IPOS(i,j)
                OCEAN(ng)%hindcj_u(i,j)=JPOS(i,j)+offset
          END IF
        END DO
      END DO
!
!  Now do local v points.
!
      OCEAN(ng)%hindcv(1)= 9999
      OCEAN(ng)%hindcv(2)=-9999
      OCEAN(ng)%hindcv(3)= 9999
      OCEAN(ng)%hindcv(4)=-9999
      DO j=Jstr,Jend
        DO i=Istr,Iend
          IF (Ipos(i,j).ge.1.0_r8) THEN
            OCEAN(ng)%hindcv(1)=MIN(OCEAN(ng)%hindcv(1),i)
            IF (Ipos(i,j).lt.REAL(Lm(ngc),r8)) THEN
              OCEAN(ng)%hindcv(2)=MAX(OCEAN(ng)%hindcv(2),i)
            END IF
          END IF
          IF (Jpos(i,j).gt.1.0_r8) THEN
            OCEAN(ng)%hindcv(3)=MIN(OCEAN(ng)%hindcv(3),j)
            IF (Jpos(i,j).lt.REAL(Mm(ngc)+1,r8)) THEN
              OCEAN(ng)%hindcv(4)=MAX(OCEAN(ng)%hindcv(4),j)
            END IF
          END IF
        END DO
      END DO
      DO j=Jstr,Jend+1
        DO i=Istr,Iend
          OCEAN(ng)%hindci_v(i,j)=IJspv
          OCEAN(ng)%hindcj_v(i,j)=IJspv
          IF (((i.ge.Xmin).and.(i.le.Xmax-1)).and.                      &
              ((j.ge.Ymin).and.(j.le.Ymax)))  THEN
                OCEAN(ng)%hindci_v(i,j)=IPOS(i,j)+offset
                OCEAN(ng)%hindcj_v(i,j)=JPOS(i,j)
          END IF
        END DO
      END DO
!
      deallocate ( Ipos, Jpos)
      deallocate ( A, lonpsi_par, latpsi_par)
!
      RETURN
      END SUBROUTINE init_child_hindices_tile
#endif
      END MODULE init_child_hindices_mod
