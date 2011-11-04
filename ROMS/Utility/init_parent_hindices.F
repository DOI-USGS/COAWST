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
      integer :: Iminp, Imaxp, Jminp, Jmaxp

      real(r8) :: IJspv, my_min, my_max, offset
      integer :: Xmin, Xmax, Ymin, Ymax, MyError, MySizer

      logical :: rectangular
      real(r8), allocatable :: A(:)
      real(r8), allocatable :: Ipos(:,:)
      real(r8), allocatable :: Jpos(:,:)
      real(r8), allocatable :: lonpsi_cld(:,:)
      real(r8), allocatable :: latpsi_cld(:,:)

      character (len=20) :: var_name
!
# include "set_bounds.h"
!
      ngp=myparent(ng)
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
          CALL get_gridcoords (ngp, ncGRDid(ngp), Ilen, Jlen,           &
     &                         var_name, A)
          CALL mpi_bcast (A, IJlen, MP_FLOAT, MyMaster,                 &
     &                    OCN_COMM_WORLD, MyError)
          ic=0
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              ic=ic+1
              REFINED(ng)%lonp_p(i,j)=A(ic)
            END DO
          END DO
!
! Case lat_psi
!
          var_name='lat_psi'
          CALL get_gridcoords (ngp, ncGRDid(ngp), Ilen, Jlen,           &
     &                         var_name, A)
          CALL mpi_bcast (A, IJlen, MP_FLOAT, MyMaster,                 &
     &                    OCN_COMM_WORLD, MyError)
          ic=0
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              ic=ic+1
              REFINED(ng)%latp_p(i,j)=A(ic)
            END DO
          END DO
!
! Case x_psi
!
          var_name='x_psi'
          CALL get_gridcoords (ngp, ncGRDid(ngp), Ilen, Jlen,           &
     &                         var_name, A)
          CALL mpi_bcast (A, IJlen, MP_FLOAT, MyMaster,                 &
     &                    OCN_COMM_WORLD, MyError)
          ic=0
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              ic=ic+1
              REFINED(ng)%xp_p(i,j)=A(ic)
            END DO
          END DO
!
! Case y_psi
!
          var_name='y_psi'
          CALL get_gridcoords (ngp, ncGRDid(ngp), Ilen, Jlen,           &
     &                         var_name, A)
          CALL mpi_bcast (A, IJlen, MP_FLOAT, MyMaster,                 &
     &                    OCN_COMM_WORLD, MyError)
          ic=0
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              ic=ic+1
              REFINED(ng)%yp_p(i,j)=A(ic)
            END DO
          END DO
# ifdef MASKING
!
! Case mask_psi
!
          var_name='mask_psi'
          CALL get_gridcoords (ngp, ncGRDid(ngp), Ilen, Jlen,           &
     &                         var_name, A)
          CALL mpi_bcast (A, IJlen, MP_FLOAT, MyMaster,                 &
     &                    OCN_COMM_WORLD, MyError)
          ic=0
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              ic=ic+1
              REFINED(ng)%pmask_p(i,j)=A(ic)
            END DO
          END DO
# endif
          IF (allocated(A)) THEN
            deallocate (A)
          END IF
        END IF
      END DO
!
! Get the horizontal indices of the grid locations in the child grid.
! Base this on the psi points to start.
! Gather all the lon and lat psi points of the parent and child.
! Use hindices to determine the location of the 4 corners of the parent
! grid. Hindices is not accurate enough to give exact locations, so 
! we use integers increments to set all the hindices.
!
      gtype=p2dvar
      IF (ng.eq.1) THEN
        Imin=0
        Imax=Lm(ng)+1
        Jmin=0
        Jmax=Mm(ng)+1
        IF ((ABS(MyType).eq.p2dvar).or.(ABS(MyType).eq.u2dvar)) THEN
          Imin=1
        END IF
        IF ((ABS(MyType).eq.p2dvar).or.(ABS(MyType).eq.v2dvar)) THEN
          Jmin=1
        END IF
      ELSE
        Imin=LOWER_BOUND_I
        Imax=Lm(ng)+GHOST_POINTS
        Jmin=LOWER_BOUND_J
        Jmax=Mm(ng)+GHOST_POINTS
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
      allocate (A(1:IJlen))
      allocate (lonpsi_cld(Imin:Imax,Jmin:Jmax))
      allocate (latpsi_cld(Imin:Imax,Jmin:Jmax))
      DO i=1,IJlen
        A(i)=0.0_r8
      END DO
      DO i=Imin,Imax
        DO j=Jmin,Jmax
          lonpsi_cld(i,j)=0.0_r8
          latpsi_cld(i,j)=0.0_r8
        END DO
      END DO
!
! Get all the child grid lon lat psi data.
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
          lonpsi_cld(i,j)=A(ic)
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
          latpsi_cld(i,j)=A(ic)
        END DO
      END DO
!
!  Now use hindices to get close enough for the parent corners. then
!  int these values.
!
      Iminp=1                    !LOWER_BOUND_I+1
      Imaxp=Lm(ngp)+1            !GHOST_POINTS
      Jminp=1                    !LOWER_BOUND_J+1
      Jmaxp=Mm(ngp)+1            !GHOST_POINTS
!
      allocate ( Ipos(Imin-1:Imax+1,Jmin-1:Jmax+1) )
      allocate ( Jpos(Imin-1:Imax+1,Jmin-1:Jmax+1) )
!
      IF (spherical) THEN
        CALL hindices (ng, Iminp, Imaxp, Jminp, Jmaxp,                  &
     &                 Iminp, Imaxp, Jminp, Jmaxp,                      &
     &                 REFINED(ng)%angler_p(Iminp:Imaxp,Jminp:Jmaxp),   &
     &                 REFINED(ng)%lonp_p(Iminp:Imaxp,Jminp:Jmaxp),     &
     &                 REFINED(ng)%latp_p(Iminp:Imaxp,Jminp:Jmaxp),     &
     &                 Imin, Imax, Jmin, Jmax,                          &
     &                 1, Lm(ng), 1, Mm(ng),                            &
     &                 lonpsi_cld, latpsi_cld,                          &
     &                 Ipos(Imin:Imax,Jmin:Jmax),                       &
     &                 Jpos(Imin:Imax,Jmin:Jmax),                       &
     &                 IJspv, rectangular, Xmin, Xmax, Ymin, Ymax)
      ELSE
        CALL hindices (ng, Iminp, Imaxp, Jminp, Jmaxp,                  &
     &                 Iminp, Imaxp, Jminp, Jmaxp,                      &
     &                 REFINED(ng)%angler_p(Iminp:Imaxp,Jminp:Jmaxp),   &
     &                 REFINED(ng)%xp_p(Iminp:Imaxp,Jminp:Jmaxp),       &
     &                 REFINED(ng)%yp_p(Iminp:Imaxp,Jminp:Jmaxp),       &
     &                 Imin, Imax, Jmin, Jmax,                          &
     &                 1, Lm(ng)+1, 1, Mm(ng)+1,                        &
     &                 lonpsi_cld, latpsi_cld,                          &
     &                 Ipos(Imin:Imax,Jmin:Jmax),                       &
     &                 Jpos(Imin:Imax,Jmin:Jmax),                       &
     &                 IJspv, rectangular, Xmin, Xmax, Ymin, Ymax)
      END IF

!  Xmin Xmax Ymin Ymax are the psi point limits of the child grid indices.
!  Now find the limits for the rho points on each tile.
!  Still global.
!
!  USe refined ratio for this child grid, computed in init_child_hind.
!
      offset=1.0_r8/REFINED(ng)%ref_scale
      Imin=IPOS(1,1)
      Imax=IPOS(Lm(ng)+1,1)
      Jmin=JPOS(1,1)
      Jmax=JPOS(1,Mm(ng)+1)
!
      DO j=Ymin-4,Ymax+3
        DO i=Xmin-4,Xmax+3
          IPOS(i,j)=Imin-offset*4.0_r8+REAL((i-(Xmin-4)),r8)*           &
     &              offset
          JPOS(i,j)=Jmin-offset*4.0_r8+REAL((j-(Ymin-4)),r8)*           &
     &              offset
        END DO
      END DO
!
!  Now do local rho points.
!
      DO j=LBj,UBj
        DO i=LBi,UBi
          OCEAN(ng)%hindpi_r(i,j)=IPOS(i,j)-offset
          OCEAN(ng)%hindpj_r(i,j)=JPOS(i,j)-offset
        END DO
      END DO
!
!  Now do local u points.
!
      DO j=LBj,UBj
        DO i=LBi,UBi
          OCEAN(ng)%hindpi_u(i,j)=IPOS(i,j)
          OCEAN(ng)%hindpj_u(i,j)=JPOS(i,j)-offset
        END DO
      END DO
!
!  Now do local v points.
!
      DO j=LBj,UBj
        DO i=LBi,UBi
          OCEAN(ng)%hindpi_v(i,j)=IPOS(i,j)-offset
          OCEAN(ng)%hindpj_v(i,j)=JPOS(i,j)
        END DO
      END DO
!
      deallocate ( Ipos, Jpos)
      deallocate ( A, lonpsi_cld, latpsi_cld)
!
      RETURN
      END SUBROUTINE init_parent_hindices_tile
#endif
      END MODULE init_parent_hindices_mod
