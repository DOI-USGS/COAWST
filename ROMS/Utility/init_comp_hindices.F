#include "cppdefs.h"
      MODULE init_comp_hindices_mod
#ifdef COMPOSED_GRID
!
!=======================================================================
!  For the entire parent grid, this routine determines the             !
!  horizontal indices of all the child composed grids                  !
!  at the locations of the parent grid.                                !
!=======================================================================
!                                                                      !
!
      implicit none

      PRIVATE
      PUBLIC  :: init_comp_hindices

      CONTAINS
!
!***********************************************************************
      SUBROUTINE init_comp_hindices (ng, tile)
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
      CALL init_comp_hindices_tile (ng, tile,                           &
     &                    LBi, UBi, LBj, UBj)
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 4)
# endif
      RETURN
      END SUBROUTINE init_comp_hindices
!
!***********************************************************************
      SUBROUTINE init_comp_hindices_tile (ng, tile,                     &
     &                                        LBi, UBi, LBj, UBj)
!***********************************************************************
!
      USE mod_param
      USE mod_composed
      USE mod_grid
      USE mod_ncparam
      USE mod_scalars
      USE mod_stepping
      USE mod_parallel
      USE mod_ocean
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
      integer :: i, j, ngp, ngc, ngg, ic, jc
      integer :: gtype, Ilen, Jlen, IJlen
      integer :: Npts
      integer :: Iminc, Imaxc, Jminc, Jmaxc
      integer :: Xmin, Xmax, Ymin, Ymax, MyError, MySizer

      real(r8) :: IJspv

      logical :: rectangular
      real(r8), parameter :: eps = 1.0E-12_r8
      real(r8), allocatable :: A(:)
      real(r8), allocatable :: x_p(:,:)
      real(r8), allocatable :: y_p(:,:)
      real(r8), allocatable :: Ipos(:,:)
      real(r8), allocatable :: Jpos(:,:)
      real(r8), allocatable :: angler_c(:,:)
      real(r8), allocatable :: x_c(:,:)
      real(r8), allocatable :: y_c(:,:)

      character (len=20) :: var_name
!
# include "set_bounds.h"
      model=iNLM
!      rectangular=.FALSE.
      rectangular=.TRUE.
!
!  This sets the child grid numbers. For example,
!  if we have 3 grids, then grid #2 has 
!  2 children that are mychild(1)=1 and mychld(2)=3.
!    
      i=0
      DO ngg=1,Ngrids
        IF (ngg.ne.ng) THEN
          i=i+1
          COMPOSED(ng)%mychild(i)=ngg
        END IF
      END DO
!
! Need to do this 3 times: rho, u, and v.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  For rho points.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! First get the horizontal indices of the parent (ie this ng grid) 
! locations.
!
      gtype=r2dvar
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
      Ilen=Imax-Imin+1
      Jlen=Jmax-Jmin+1
      IJlen=Ilen*Jlen
      MySizer=IJlen
      Npts=IJlen

      allocate (A(1:IJlen))
      allocate (x_p(Imin:Imax,Jmin:Jmax))
      allocate (y_p(Imin:Imax,Jmin:Jmax))
      DO i=1,IJlen
        A(i)=0.0_r8
      END DO
      DO i=Imin,Imax
        DO j=Jmin,Jmax
          x_p(i,j)=0.0_r8
          y_p(i,j)=0.0_r8
        END DO
      END DO
!
! Get all the parent grid lon lat psi data.
!
      IF (spherical) THEN
        CALL mp_gather2d (ng, model, LBi, UBi, LBj, UBj, 0,             &
     &                    gtype, 1.0_r8,                                &
# ifdef MASKING
     &                    GRID(ng) % rmask(LBi:UBi,LBj:UBj),            &
# endif
     &                    GRID(ng) % lonr(LBi:UBi,LBj:UBj),             &
     &                    Npts, A)
      ELSE
        CALL mp_gather2d (ng, model, LBi, UBi, LBj, UBj, 0,             &
     &                    gtype, 1.0_r8,                                &
# ifdef MASKING
     &                    GRID(ng) % rmask(LBi:UBi,LBj:UBj),            &
# endif
     &                    GRID(ng) % xr(LBi:UBi,LBj:UBj),               &
     &                    Npts, A)
      END IF
!
! distribute the A (1D packed) field to all the nodes.
!
      CALL mpi_bcast (A, MySizer, MP_FLOAT, MyMaster,                   &
     &                OCN_COMM_WORLD, MyError)
!
! unpack and place in global ng=1 2D array
!
      ic=0
      DO j=Jmin,Jmax
        DO i=Imin,Imax
          ic=ic+1
          x_p(i,j)=A(ic)
        END DO
      END DO
      IF (spherical) THEN
        CALL mp_gather2d (ng, model, LBi, UBi, LBj, UBj, 0,             &
     &                    gtype, 1.0_r8,                                &
# ifdef MASKING
     &                    GRID(ng) % rmask(LBi:UBi,LBj:UBj),            &
# endif
     &                    GRID(ng) % latr(LBi:UBi,LBj:UBj),             &
     &                    Npts, A)
      ELSE
        CALL mp_gather2d (ng, model, LBi, UBi, LBj, UBj, 0,             &
     &                    gtype, 1.0_r8,                                &
# ifdef MASKING
     &                    GRID(ng) % rmask(LBi:UBi,LBj:UBj),            &
# endif
     &                    GRID(ng) % yr(LBi:UBi,LBj:UBj),               &
     &                    Npts, A)
      END IF
!
! distribute the A (1D packed) field to all the nodes.
!
      CALL mpi_bcast (A, MySizer, MP_FLOAT, MyMaster,                   &
     &                OCN_COMM_WORLD, MyError)
!
! unpack and place in global ng=1 2D array
!
      ic=0
      DO j=Jmin,Jmax
        DO i=Imin,Imax
          ic=ic+1
          y_p(i,j)=A(ic)
        END DO
      END DO
      IF (allocated(A)) THEN
        deallocate (A)
      END IF
!
! Now for each child grid, get the hindices.
!
      DO ngg=1,Ngrids-1
        ngc=COMPOSED(ng)%mychild(ngg)
        model=iNLM
        IJspv=-9999.0_r8
!
! Get the bounds
!
        gtype=r2dvar
        Iminc=LOWER_BOUND_I
        Imaxc=Lm(ngc)+GHOST_POINTS
        Jminc=LOWER_BOUND_J
        Jmaxc=Mm(ngc)+GHOST_POINTS
        Ilen=Imaxc-Iminc+1
        Jlen=Jmaxc-Jminc+1
        IJlen=Ilen*Jlen
        Npts=IJlen
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
! Case angle
!
!
        IF (.not.allocated(angler_c)) THEN
          allocate (angler_c(Iminc:Imaxc,JminC:Jmaxc))
        END IF
        var_name='angle'
        CALL get_gridcoords (ngc, ncGRDid(ngc), Ilen, Jlen,           &
     &                       var_name, A)
        CALL mpi_bcast (A, IJlen, MP_FLOAT, MyMaster,                 &
     &                  OCN_COMM_WORLD, MyError)
        ic=0
        DO j=Jminc,Jmaxc
          DO i=Iminc,Imaxc
            ic=ic+1
            angler_c(i,j)=A(ic)
          END DO
        END DO
        IF (allocated(A)) THEN
          deallocate (A)
        END IF
!
! Now get rho/u/v point info.
!
        gtype=r2dvar
        Iminc=LOWER_BOUND_I
        Imaxc=Lm(ngc)+GHOST_POINTS
        Jminc=LOWER_BOUND_J
        Jmaxc=Mm(ngc)+GHOST_POINTS
        IF ((ABS(gtype).eq.p2dvar).or.(ABS(gtype).eq.u2dvar)) THEN
          Iminc=LOWER_BOUND_I+1
        END IF
        IF ((ABS(gtype).eq.p2dvar).or.(ABS(gtype).eq.v2dvar)) THEN
          Jminc=LOWER_BOUND_J+1
        END IF
        Ilen=Imaxc-Iminc+1
        Jlen=Jmaxc-Jminc+1
        IJlen=Ilen*Jlen
        Npts=IJlen
!
        IF (.not.allocated(A)) THEN
          allocate (A(1:IJlen))
        END IF
        DO i=1,IJlen
          A(i)=0.0_r8
        END DO
        IF (.not.allocated(x_c)) THEN
          allocate (x_c(Iminc:Imaxc,JminC:Jmaxc))
        END IF
        IF (.not.allocated(y_c)) THEN
          allocate (y_c(Iminc:Imaxc,JminC:Jmaxc))
        END IF
        IF (spherical) THEN
!
! Case lon_
!
          var_name='lon_rho'
          CALL get_gridcoords (ngc, ncGRDid(ngc), Ilen, Jlen,           &
     &                         var_name, A)
          CALL mpi_bcast (A, IJlen, MP_FLOAT, MyMaster,                 &
     &                    OCN_COMM_WORLD, MyError)
          ic=0
          DO j=Jminc,Jmaxc
            DO i=Iminc,Imaxc
              ic=ic+1
              x_c(i,j)=A(ic)
            END DO
          END DO
!
! Case lat_
!
          var_name='lat_rho'
          CALL get_gridcoords (ngc, ncGRDid(ngc), Ilen, Jlen,           &
     &                         var_name, A)
          CALL mpi_bcast (A, IJlen, MP_FLOAT, MyMaster,                 &
     &                    OCN_COMM_WORLD, MyError)
          ic=0
          DO j=Jminc,Jmaxc
            DO i=Iminc,Imaxc
              ic=ic+1
              y_c(i,j)=A(ic)
            END DO
          END DO
        ELSE
!
! Case x_
!
          var_name='x_rho'
          CALL get_gridcoords (ngc, ncGRDid(ngc), Ilen, Jlen,           &
     &                         var_name, A)
          CALL mpi_bcast (A, IJlen, MP_FLOAT, MyMaster,                 &
     &                    OCN_COMM_WORLD, MyError)
          ic=0
          DO j=Jminc,Jmaxc
            DO i=Iminc,Imaxc
              ic=ic+1
              x_c(i,j)=A(ic)
            END DO
          END DO
!
! Case y_
!
          var_name='y_rho'
          CALL get_gridcoords (ngc, ncGRDid(ngc), Ilen, Jlen,           &
     &                         var_name, A)
          CALL mpi_bcast (A, IJlen, MP_FLOAT, MyMaster,                 &
     &                    OCN_COMM_WORLD, MyError)
          ic=0
          DO j=Jminc,Jmaxc
            DO i=Iminc,Imaxc
              ic=ic+1
              y_c(i,j)=A(ic)
            END DO
          END DO
        END IF
        IF (allocated(A)) THEN
          deallocate (A)
        END IF
!
! Use hindices to determine the location of the child
! grid. Hindices is not very accurate, but for now the grids
! need to be co-incident.
!
        allocate ( Ipos(Imin:Imax,Jmin:Jmax) )
        allocate ( Jpos(Imin:Imax,Jmin:Jmax) )
!
        Xmin= 9999
        Xmax=-9999
        Ymin= 9999
        Ymax=-9999
        DO j=Jmin,Jmax         ! parent grid
          DO i=Imin,Imax       ! parent grid
            Ipos(i,j)=IJspv
            Jpos(i,j)=IJspv
            DO jc=Jminc,Jmaxc         ! child grid
              DO ic=Iminc,Imaxc       ! child grid
                IF (((x_p(i,j)-eps).lt.(x_c(ic,jc))).and.            &
     &              ((x_p(i,j)+eps).gt.(x_c(ic,jc)))) THEN
                  IF (((y_p(i,j)-eps).lt.(y_c(ic,jc))).and.          &
     &                ((y_p(i,j)+eps).gt.(y_c(ic,jc)))) THEN
                    Ipos(i,j)=ic
                    Jpos(i,j)=jc
                  END IF
                END IF
              END DO
            END DO
          END DO
        END DO
!
!  Xmin Xmax Ymin Ymax are the bounding indices
!  of the parent grid that identify the boundary of the child grid.
!  Still working on the global parent grid (ie non-tiled).
!
!  The IPOS and JPOS are values of the child indices. Only 
!  identify locations around the parent perimeter to be 
!  replaced with child data.
!
        DO j=JstrT,JendT
          DO i=IstrT,IendT
            OCEAN(ng)%hindi_r(i,j)=IJspv
            OCEAN(ng)%hindj_r(i,j)=IJspv
          END DO
        END DO
!
! Southern Edge.
!
        DO j=JstrT,Jstr-1
          DO i=IstrT,IendT
            OCEAN(ng)%hindi_r(i,j)=Ipos(i,j)
            OCEAN(ng)%hindj_r(i,j)=Jpos(i,j)
            IF (Ipos(i,j).ne.IJspv) THEN
              Xmin=MIN(Xmin,i)
              Xmax=MAX(Xmax,i)
            END IF
            IF (Jpos(i,j).ne.IJspv) THEN
              Ymin=MIN(Ymin,j)
              Ymax=MAX(Ymax,j)
            END IF
          END DO
        END DO
!
! Northern Edge.
!
        DO j=Jend+1,JendT
          DO i=IstrT,IendT
            OCEAN(ng)%hindi_r(i,j)=Ipos(i,j)
            OCEAN(ng)%hindj_r(i,j)=Jpos(i,j)
            IF (Ipos(i,j).ne.IJspv) THEN
              Xmin=MIN(Xmin,i)
              Xmax=MAX(Xmax,i)
            END IF
            IF (Jpos(i,j).ne.IJspv) THEN
              Ymin=MIN(Ymin,j)
              Ymax=MAX(Ymax,j)
            END IF
          END DO
        END DO
!
! Western Edge.
!
        DO j=JstrT,JendT
          DO i=IstrT,Istr-1
            OCEAN(ng)%hindi_r(i,j)=Ipos(i,j)
            OCEAN(ng)%hindj_r(i,j)=Jpos(i,j)
            IF (Ipos(i,j).ne.IJspv) THEN
              Xmin=MIN(Xmin,i)
              Xmax=MAX(Xmax,i)
            END IF
            IF (Jpos(i,j).ne.IJspv) THEN
              Ymin=MIN(Ymin,j)
              Ymax=MAX(Ymax,j)
            END IF
          END DO
        END DO
!
! Eastern Edge.
!
        DO j=JstrT,JendT
          DO i=Iend+1,IendT
            OCEAN(ng)%hindi_r(i,j)=Ipos(i,j)
            OCEAN(ng)%hindj_r(i,j)=Jpos(i,j)
            IF (Ipos(i,j).ne.IJspv) THEN
              Xmin=MIN(Xmin,i)
              Xmax=MAX(Xmax,i)
            END IF
            IF (Jpos(i,j).ne.IJspv) THEN
              Ymin=MIN(Ymin,j)
              Ymax=MAX(Ymax,j)
            END IF
          END DO
        END DO
!
        OCEAN(ng)%hindr(1,ngg)=Xmin
        OCEAN(ng)%hindr(2,ngg)=Xmax
        OCEAN(ng)%hindr(3,ngg)=Ymin
        OCEAN(ng)%hindr(4,ngg)=Ymax
!
        deallocate ( Ipos, Jpos)
        deallocate ( x_c, y_c)
!
      END DO
      deallocate ( x_p, y_p)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  For u points.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! First get the horizontal indices of the parent (ie this ng grid) 
! locations.
!
      gtype=u2dvar
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
      Ilen=Imax-Imin+1
      Jlen=Jmax-Jmin+1
      IJlen=Ilen*Jlen
      MySizer=IJlen
      Npts=IJlen

      allocate (A(1:IJlen))
      allocate (x_p(Imin:Imax,Jmin:Jmax))
      allocate (y_p(Imin:Imax,Jmin:Jmax))
      DO i=1,IJlen
        A(i)=0.0_r8
      END DO
      DO i=Imin,Imax
        DO j=Jmin,Jmax
          x_p(i,j)=0.0_r8
          y_p(i,j)=0.0_r8
        END DO
      END DO
!
! Get all the parent grid lon lat psi data.
!
      IF (spherical) THEN
        CALL mp_gather2d (ng, model, LBi, UBi, LBj, UBj, 0,             &
     &                    gtype, 1.0_r8,                                &
# ifdef MASKING
     &                    GRID(ng) % umask(LBi:UBi,LBj:UBj),            &
# endif
     &                    GRID(ng) % lonu(LBi:UBi,LBj:UBj),             &
     &                    Npts, A)
      ELSE
        CALL mp_gather2d (ng, model, LBi, UBi, LBj, UBj, 0,             &
     &                    gtype, 1.0_r8,                                &
# ifdef MASKING
     &                    GRID(ng) % umask(LBi:UBi,LBj:UBj),            &
# endif
     &                    GRID(ng) % xu(LBi:UBi,LBj:UBj),               &
     &                    Npts, A)
      END IF
!
! distribute the A (1D packed) field to all the nodes.
!
      CALL mpi_bcast (A, MySizer, MP_FLOAT, MyMaster,                   &
     &                OCN_COMM_WORLD, MyError)
!
! unpack and place in global ng=1 2D array
!
      ic=0
      DO j=Jmin,Jmax
        DO i=Imin,Imax
          ic=ic+1
          x_p(i,j)=A(ic)
        END DO
      END DO
      IF (spherical) THEN
        CALL mp_gather2d (ng, model, LBi, UBi, LBj, UBj, 0,             &
     &                    gtype, 1.0_r8,                                &
# ifdef MASKING
     &                    GRID(ng) % umask(LBi:UBi,LBj:UBj),            &
# endif
     &                    GRID(ng) % latu(LBi:UBi,LBj:UBj),             &
     &                    Npts, A)
      ELSE
        CALL mp_gather2d (ng, model, LBi, UBi, LBj, UBj, 0,             &
     &                    gtype, 1.0_r8,                                &
# ifdef MASKING
     &                    GRID(ng) % umask(LBi:UBi,LBj:UBj),            &
# endif
     &                    GRID(ng) % yu(LBi:UBi,LBj:UBj),               &
     &                    Npts, A)
      END IF
!
! distribute the A (1D packed) field to all the nodes.
!
      CALL mpi_bcast (A, MySizer, MP_FLOAT, MyMaster,                   &
     &                OCN_COMM_WORLD, MyError)
!
! unpack and place in global ng=1 2D array
!
      ic=0
      DO j=Jmin,Jmax
        DO i=Imin,Imax
          ic=ic+1
          y_p(i,j)=A(ic)
        END DO
      END DO
      IF (allocated(A)) THEN
        deallocate (A)
      END IF
!
! Now for each child grid, get the hindices.
!
      DO ngg=1,Ngrids-1
        ngc=COMPOSED(ng)%mychild(ngg)
        model=iNLM
        IJspv=-9999.0_r8
!
! Now get rho/u/v point info.
!
        gtype=u2dvar
        Iminc=LOWER_BOUND_I
        Imaxc=Lm(ngc)+GHOST_POINTS
        Jminc=LOWER_BOUND_J
        Jmaxc=Mm(ngc)+GHOST_POINTS
        IF ((ABS(gtype).eq.p2dvar).or.(ABS(gtype).eq.u2dvar)) THEN
          Iminc=LOWER_BOUND_I+1
        END IF
        IF ((ABS(gtype).eq.p2dvar).or.(ABS(gtype).eq.v2dvar)) THEN
          Jminc=LOWER_BOUND_J+1
        END IF
        Ilen=Imaxc-Iminc+1
        Jlen=Jmaxc-Jminc+1
        IJlen=Ilen*Jlen
        Npts=IJlen
!
        IF (.not.allocated(A)) THEN
          allocate (A(1:IJlen))
        END IF
        DO i=1,IJlen
          A(i)=0.0_r8
        END DO
        IF (.not.allocated(x_c)) THEN
          allocate (x_c(Iminc:Imaxc,JminC:Jmaxc))
        END IF
        IF (.not.allocated(y_c)) THEN
          allocate (y_c(Iminc:Imaxc,JminC:Jmaxc))
        END IF
        IF (spherical) THEN
!
! Case lon_
!
          var_name='lon_u'
          CALL get_gridcoords (ngc, ncGRDid(ngc), Ilen, Jlen,           &
     &                         var_name, A)
          CALL mpi_bcast (A, IJlen, MP_FLOAT, MyMaster,                 &
     &                    OCN_COMM_WORLD, MyError)
          ic=0
          DO j=Jminc,Jmaxc
            DO i=Iminc,Imaxc
              ic=ic+1
              x_c(i,j)=A(ic)
            END DO
          END DO
!
! Case lat_
!
          var_name='lat_u'
          CALL get_gridcoords (ngc, ncGRDid(ngc), Ilen, Jlen,           &
     &                         var_name, A)
          CALL mpi_bcast (A, IJlen, MP_FLOAT, MyMaster,                 &
     &                    OCN_COMM_WORLD, MyError)
          ic=0
          DO j=Jminc,Jmaxc
            DO i=Iminc,Imaxc
              ic=ic+1
              y_c(i,j)=A(ic)
            END DO
          END DO
        ELSE
!
! Case x_
!
          var_name='x_u'
          CALL get_gridcoords (ngc, ncGRDid(ngc), Ilen, Jlen,           &
     &                         var_name, A)
          CALL mpi_bcast (A, IJlen, MP_FLOAT, MyMaster,                 &
     &                    OCN_COMM_WORLD, MyError)
          ic=0
          DO j=Jminc,Jmaxc
            DO i=Iminc,Imaxc
              ic=ic+1
              x_c(i,j)=A(ic)
            END DO
          END DO
!
! Case y_
!
          var_name='y_u'
          CALL get_gridcoords (ngc, ncGRDid(ngc), Ilen, Jlen,           &
     &                         var_name, A)
          CALL mpi_bcast (A, IJlen, MP_FLOAT, MyMaster,                 &
     &                    OCN_COMM_WORLD, MyError)
          ic=0
          DO j=Jminc,Jmaxc
            DO i=Iminc,Imaxc
              ic=ic+1
              y_c(i,j)=A(ic)
            END DO
          END DO
        END IF
        IF (allocated(A)) THEN
          deallocate (A)
        END IF
!
! Use hindices to determine the location of the child
! grid. Hindices is not very accurate, but for now the grids
! need to be co-incident.
!
        allocate ( Ipos(Imin:Imax,Jmin:Jmax) )
        allocate ( Jpos(Imin:Imax,Jmin:Jmax) )
!
        Xmin= 9999
        Xmax=-9999
        Ymin= 9999
        Ymax=-9999
        DO j=Jmin,Jmax         ! parent grid
          DO i=Imin,Imax       ! parent grid
            Ipos(i,j)=IJspv
            Jpos(i,j)=IJspv
            DO jc=Jminc,Jmaxc         ! child grid
              DO ic=Iminc,Imaxc       ! child grid
                IF (((x_p(i,j)-eps).lt.(x_c(ic,jc))).and.            &
     &              ((x_p(i,j)+eps).gt.(x_c(ic,jc)))) THEN
                  IF (((y_p(i,j)-eps).lt.(y_c(ic,jc))).and.          &
     &                ((y_p(i,j)+eps).gt.(y_c(ic,jc)))) THEN
                    Ipos(i,j)=ic
                    Jpos(i,j)=jc
                  END IF
                END IF
              END DO
            END DO
          END DO
        END DO
!
!  Xmin Xmax Ymin Ymax are the bounding indices
!  of the parent grid that identify the boundary of the child grid.
!  Still working on the global parent grid (ie non-tiled).
!
!  The IPOS and JPOS are values of the child indices. Only 
!  identify locations around the parent perimeter to be 
!  replaced with child data.
!
        DO j=JstrT,JendT
          DO i=IstrT,IendT
            OCEAN(ng)%hindi_u(i,j)=IJspv
            OCEAN(ng)%hindj_u(i,j)=IJspv
          END DO
        END DO
!
! Southern Edge.
!
        DO j=JstrT,Jstr-1
          DO i=IstrTU+1,IendT
            OCEAN(ng)%hindi_u(i,j)=Ipos(i,j)
            OCEAN(ng)%hindj_u(i,j)=Jpos(i,j)
            IF (Ipos(i,j).ne.IJspv) THEN
              Xmin=MIN(Xmin,i)
              Xmax=MAX(Xmax,i)
            END IF
            IF (Jpos(i,j).ne.IJspv) THEN
              Ymin=MIN(Ymin,j)
              Ymax=MAX(Ymax,j)
            END IF
          END DO
        END DO
!
! Northern Edge.
!
        DO j=Jend+1,JendT
          DO i=IstrTU+1,IendT
            OCEAN(ng)%hindi_u(i,j)=Ipos(i,j)
            OCEAN(ng)%hindj_u(i,j)=Jpos(i,j)
            IF (Ipos(i,j).ne.IJspv) THEN
              Xmin=MIN(Xmin,i)
              Xmax=MAX(Xmax,i)
            END IF
            IF (Jpos(i,j).ne.IJspv) THEN
              Ymin=MIN(Ymin,j)
              Ymax=MAX(Ymax,j)
            END IF
          END DO
        END DO
!
! Western Edge.
!
        DO j=JstrT,JendT
          DO i=IstrTU+1,Istr-1
            OCEAN(ng)%hindi_u(i,j)=Ipos(i,j)
            OCEAN(ng)%hindj_u(i,j)=Jpos(i,j)
            IF (Ipos(i,j).ne.IJspv) THEN
              Xmin=MIN(Xmin,i)
              Xmax=MAX(Xmax,i)
            END IF
            IF (Jpos(i,j).ne.IJspv) THEN
              Ymin=MIN(Ymin,j)
              Ymax=MAX(Ymax,j)
            END IF
          END DO
        END DO
!
! Eastern Edge.
!
        DO j=JstrT,JendT
          DO i=Iend+1,IendT
            OCEAN(ng)%hindi_u(i,j)=Ipos(i,j)
            OCEAN(ng)%hindj_u(i,j)=Jpos(i,j)
            IF (Ipos(i,j).ne.IJspv) THEN
              Xmin=MIN(Xmin,i)
              Xmax=MAX(Xmax,i)
            END IF
            IF (Jpos(i,j).ne.IJspv) THEN
              Ymin=MIN(Ymin,j)
              Ymax=MAX(Ymax,j)
            END IF
          END DO
        END DO
!
        OCEAN(ng)%hindu(1,ngg)=Xmin
        OCEAN(ng)%hindu(2,ngg)=Xmax
        OCEAN(ng)%hindu(3,ngg)=Ymin
        OCEAN(ng)%hindu(4,ngg)=Ymax
!
        deallocate ( Ipos, Jpos)
        deallocate ( x_c, y_c)
!
      END DO
      deallocate ( x_p, y_p)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  For v points.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! First get the horizontal indices of the parent (ie this ng grid) 
! locations.
!
      gtype=v2dvar
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
      Ilen=Imax-Imin+1
      Jlen=Jmax-Jmin+1
      IJlen=Ilen*Jlen
      MySizer=IJlen
      Npts=IJlen

      allocate (A(1:IJlen))
      allocate (x_p(Imin:Imax,Jmin:Jmax))
      allocate (y_p(Imin:Imax,Jmin:Jmax))
      DO i=1,IJlen
        A(i)=0.0_r8
      END DO
      DO i=Imin,Imax
        DO j=Jmin,Jmax
          x_p(i,j)=0.0_r8
          y_p(i,j)=0.0_r8
        END DO
      END DO
!
! Get all the parent grid lon lat psi data.
!
      IF (spherical) THEN
        CALL mp_gather2d (ng, model, LBi, UBi, LBj, UBj, 0,             &
     &                    gtype, 1.0_r8,                                &
# ifdef MASKING
     &                    GRID(ng) % vmask(LBi:UBi,LBj:UBj),            &
# endif
     &                    GRID(ng) % lonv(LBi:UBi,LBj:UBj),             &
     &                    Npts, A)
      ELSE
        CALL mp_gather2d (ng, model, LBi, UBi, LBj, UBj, 0,             &
     &                    gtype, 1.0_r8,                                &
# ifdef MASKING
     &                    GRID(ng) % vmask(LBi:UBi,LBj:UBj),            &
# endif
     &                    GRID(ng) % xv(LBi:UBi,LBj:UBj),               &
     &                    Npts, A)
      END IF
!
! distribute the A (1D packed) field to all the nodes.
!
      CALL mpi_bcast (A, MySizer, MP_FLOAT, MyMaster,                   &
     &                OCN_COMM_WORLD, MyError)
!
! unpack and place in global ng=1 2D array
!
      ic=0
      DO j=Jmin,Jmax
        DO i=Imin,Imax
          ic=ic+1
          x_p(i,j)=A(ic)
        END DO
      END DO
      IF (spherical) THEN
        CALL mp_gather2d (ng, model, LBi, UBi, LBj, UBj, 0,             &
     &                    gtype, 1.0_r8,                                &
# ifdef MASKING
     &                    GRID(ng) % vmask(LBi:UBi,LBj:UBj),            &
# endif
     &                    GRID(ng) % latv(LBi:UBi,LBj:UBj),             &
     &                    Npts, A)
      ELSE
        CALL mp_gather2d (ng, model, LBi, UBi, LBj, UBj, 0,             &
     &                    gtype, 1.0_r8,                                &
# ifdef MASKING
     &                    GRID(ng) % vmask(LBi:UBi,LBj:UBj),            &
# endif
     &                    GRID(ng) % yv(LBi:UBi,LBj:UBj),               &
     &                    Npts, A)
      END IF
!
! distribute the A (1D packed) field to all the nodes.
!
      CALL mpi_bcast (A, MySizer, MP_FLOAT, MyMaster,                   &
     &                OCN_COMM_WORLD, MyError)
!
! unpack and place in global ng=1 2D array
!
      ic=0
      DO j=Jmin,Jmax
        DO i=Imin,Imax
          ic=ic+1
          y_p(i,j)=A(ic)
        END DO
      END DO
      IF (allocated(A)) THEN
        deallocate (A)
      END IF
!
! Now for each child grid, get the hindices.
!
      DO ngg=1,Ngrids-1
        ngc=COMPOSED(ng)%mychild(ngg)
        model=iNLM
        IJspv=-9999.0_r8
!
! Now get rho/u/v point info.
!
        gtype=v2dvar
        Iminc=LOWER_BOUND_I
        Imaxc=Lm(ngc)+GHOST_POINTS
        Jminc=LOWER_BOUND_J
        Jmaxc=Mm(ngc)+GHOST_POINTS
        IF ((ABS(gtype).eq.p2dvar).or.(ABS(gtype).eq.u2dvar)) THEN
          Iminc=LOWER_BOUND_I+1
        END IF
        IF ((ABS(gtype).eq.p2dvar).or.(ABS(gtype).eq.v2dvar)) THEN
          Jminc=LOWER_BOUND_J+1
        END IF
        Ilen=Imaxc-Iminc+1
        Jlen=Jmaxc-Jminc+1
        IJlen=Ilen*Jlen
        Npts=IJlen
!
        IF (.not.allocated(A)) THEN
          allocate (A(1:IJlen))
        END IF
        DO i=1,IJlen
          A(i)=0.0_r8
        END DO
        IF (.not.allocated(x_c)) THEN
          allocate (x_c(Iminc:Imaxc,JminC:Jmaxc))
        END IF
        IF (.not.allocated(y_c)) THEN
          allocate (y_c(Iminc:Imaxc,JminC:Jmaxc))
        END IF
        IF (spherical) THEN
!
! Case lon_
!
          var_name='lon_v'
          CALL get_gridcoords (ngc, ncGRDid(ngc), Ilen, Jlen,           &
     &                         var_name, A)
          CALL mpi_bcast (A, IJlen, MP_FLOAT, MyMaster,                 &
     &                    OCN_COMM_WORLD, MyError)
          ic=0
          DO j=Jminc,Jmaxc
            DO i=Iminc,Imaxc
              ic=ic+1
              x_c(i,j)=A(ic)
            END DO
          END DO
!
! Case lat_
!
          var_name='lat_v'
          CALL get_gridcoords (ngc, ncGRDid(ngc), Ilen, Jlen,           &
     &                         var_name, A)
          CALL mpi_bcast (A, IJlen, MP_FLOAT, MyMaster,                 &
     &                    OCN_COMM_WORLD, MyError)
          ic=0
          DO j=Jminc,Jmaxc
            DO i=Iminc,Imaxc
              ic=ic+1
              y_c(i,j)=A(ic)
            END DO
          END DO
        ELSE
!
! Case x_
!
          var_name='x_v'
          CALL get_gridcoords (ngc, ncGRDid(ngc), Ilen, Jlen,           &
     &                         var_name, A)
          CALL mpi_bcast (A, IJlen, MP_FLOAT, MyMaster,                 &
     &                    OCN_COMM_WORLD, MyError)
          ic=0
          DO j=Jminc,Jmaxc
            DO i=Iminc,Imaxc
              ic=ic+1
              x_c(i,j)=A(ic)
            END DO
          END DO
!
! Case y_
!
          var_name='y_v'
          CALL get_gridcoords (ngc, ncGRDid(ngc), Ilen, Jlen,           &
     &                         var_name, A)
          CALL mpi_bcast (A, IJlen, MP_FLOAT, MyMaster,                 &
     &                    OCN_COMM_WORLD, MyError)
          ic=0
          DO j=Jminc,Jmaxc
            DO i=Iminc,Imaxc
              ic=ic+1
              y_c(i,j)=A(ic)
            END DO
          END DO
        END IF
        IF (allocated(A)) THEN
          deallocate (A)
        END IF
!
! Use hindices to determine the location of the child
! grid. Hindices is not very accurate, but for now the grids
! need to be co-incident.
!
        allocate ( Ipos(Imin:Imax,Jmin:Jmax) )
        allocate ( Jpos(Imin:Imax,Jmin:Jmax) )
!
        Xmin= 9999
        Xmax=-9999
        Ymin= 9999
        Ymax=-9999
        DO j=Jmin,Jmax         ! parent grid
          DO i=Imin,Imax       ! parent grid
            Ipos(i,j)=IJspv
            Jpos(i,j)=IJspv
            DO jc=Jminc,Jmaxc         ! child grid
              DO ic=Iminc,Imaxc       ! child grid
                IF (((x_p(i,j)-eps).lt.(x_c(ic,jc))).and.            &
     &              ((x_p(i,j)+eps).gt.(x_c(ic,jc)))) THEN
                  IF (((y_p(i,j)-eps).lt.(y_c(ic,jc))).and.          &
     &                ((y_p(i,j)+eps).gt.(y_c(ic,jc)))) THEN
                    Ipos(i,j)=ic
                    Jpos(i,j)=jc
                  END IF
                END IF
              END DO
            END DO
          END DO
        END DO
!
!
!  Xmin Xmax Ymin Ymax are the bounding indices
!  of the parent grid that identify the boundary of the child grid.
!  Still working on the global parent grid (ie non-tiled).
!
!  The IPOS and JPOS are values of the child indices. Only 
!  identify locations around the parent perimeter to be 
!  replaced with child data.
!
        DO j=JstrT,JendT
          DO i=IstrT,IendT
            OCEAN(ng)%hindi_v(i,j)=IJspv
            OCEAN(ng)%hindj_v(i,j)=IJspv
          END DO
        END DO
!
! Southern Edge.
!
        DO j=JstrTV+1,Jstr-1
          DO i=IstrT,IendT
            OCEAN(ng)%hindi_v(i,j)=Ipos(i,j)
            OCEAN(ng)%hindj_v(i,j)=Jpos(i,j)
            IF (Ipos(i,j).ne.IJspv) THEN
              Xmin=MIN(Xmin,i)
              Xmax=MAX(Xmax,i)
            END IF
            IF (Jpos(i,j).ne.IJspv) THEN
              Ymin=MIN(Ymin,j)
              Ymax=MAX(Ymax,j)
            END IF
          END DO
        END DO
!
! Northern Edge.
!
        DO j=Jend+1,JendT
          DO i=IstrT,IendT
            OCEAN(ng)%hindi_v(i,j)=Ipos(i,j)
            OCEAN(ng)%hindj_v(i,j)=Jpos(i,j)
            IF (Ipos(i,j).ne.IJspv) THEN
              Xmin=MIN(Xmin,i)
              Xmax=MAX(Xmax,i)
            END IF
            IF (Jpos(i,j).ne.IJspv) THEN
              Ymin=MIN(Ymin,j)
              Ymax=MAX(Ymax,j)
            END IF
          END DO
        END DO
!
! Western Edge.
!
        DO j=JstrTV+1,JendT
          DO i=IstrT,Istr-1
            OCEAN(ng)%hindi_v(i,j)=Ipos(i,j)
            OCEAN(ng)%hindj_v(i,j)=Jpos(i,j)
            IF (Ipos(i,j).ne.IJspv) THEN
              Xmin=MIN(Xmin,i)
              Xmax=MAX(Xmax,i)
            END IF
            IF (Jpos(i,j).ne.IJspv) THEN
              Ymin=MIN(Ymin,j)
              Ymax=MAX(Ymax,j)
            END IF
          END DO
        END DO
!
! Eastern Edge.
!
        DO j=JstrTV+1,JendT
          DO i=Iend+1,IendT
            OCEAN(ng)%hindi_v(i,j)=Ipos(i,j)
            OCEAN(ng)%hindj_v(i,j)=Jpos(i,j)
            IF (Ipos(i,j).ne.IJspv) THEN
              Xmin=MIN(Xmin,i)
              Xmax=MAX(Xmax,i)
            END IF
            IF (Jpos(i,j).ne.IJspv) THEN
              Ymin=MIN(Ymin,j)
              Ymax=MAX(Ymax,j)
            END IF
          END DO
        END DO
!
        OCEAN(ng)%hindv(1,ngg)=Xmin
        OCEAN(ng)%hindv(2,ngg)=Xmax
        OCEAN(ng)%hindv(3,ngg)=Ymin
        OCEAN(ng)%hindv(4,ngg)=Ymax
!
        deallocate ( Ipos, Jpos)
        deallocate ( angler_c, x_c, y_c)
!
      END DO
      deallocate ( x_p, y_p)

      RETURN
      END SUBROUTINE init_comp_hindices_tile
#endif
      END MODULE init_comp_hindices_mod
