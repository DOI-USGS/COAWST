#include "cppdefs.h"
#define INTERP_BILINEAR
      MODULE get_2dparent_data_mod
#ifdef REFINED_GRID
!
!svn $Id: set_data.F 1405 2008-02-20 01:44:07Z jcwarner $
!================================================== John C. Warner =====
!                                                                      !
!  This routine acquires data from a parent grid and                   !
!  loaded into global two-time record arrays.                          !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: get_2dparent_data, get_2dparent_data_tile

      CONTAINS
!
!***********************************************************************
      SUBROUTINE get_2dparent_data (ng, tile)
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
      integer :: ngp, LBip, UBip, LBjp, UBjp
# include "tile.h"

      ngp=myparent(ng)
      LBip=BOUNDS(ngp)%LBi(tile)
      UBip=BOUNDS(ngp)%UBi(tile)
      LBjp=BOUNDS(ngp)%LBj(tile)
      UBjp=BOUNDS(ngp)%UBj(tile)
!
# ifdef PROFILE
      CALL wclock_on (ng, iNLM, 4)
# endif
      CALL get_2dparent_data_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             LBip, UBip, LBjp, UBjp)
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 4)
# endif
      RETURN
      END SUBROUTINE get_2dparent_data
!
!***********************************************************************
      SUBROUTINE get_2dparent_data_tile (ng, tile,                      &
     &                                   LBi, UBi, LBj, UBj,            &
     &                                   LBip, UBip, LBjp, UBjp)
!***********************************************************************
!
      USE mod_param
      USE mod_boundary
      USE mod_grid
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
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: LBip, UBip, LBjp, UBjp
!
!  Local variable declarations.
!
      logical, dimension(3) :: update =                                 &
     &         (/ .FALSE., .FALSE., .FALSE. /)

# ifdef OBC
      integer :: ILB, IUB, JLB, JUB
# endif
      integer :: i, j, ii, gtype, ic, jc, Npts, ngp, t1
      integer :: Imin, Imax, Jmin, Jmax, model
      integer :: Ilen, Jlen, IJlen, MyType, status
      integer :: MyError, MySize, ncvarid, refscale
      real(r8) :: my_min, my_max, cff, cff1, ic1, jc1

      real(r8), allocatable :: Aout(:)
      real(r8), allocatable :: A(:,:)
      real(r8), allocatable :: Amask(:,:)
      real(r8), allocatable :: A2(:,:)

      logical :: rectangular
      real(r8) :: IJspv
      character (len=20) :: Xname, Yname

!
# include "set_bounds.h"
!
      ngp=myparent(ng)
      model=iNLM

# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn on input data time wall clock.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, iNLM, 3)
# endif
!
!set time index to fill
      IF (iic(ngp).eq.ntstart(ngp)) THEN
        t1=2
      ELSE
# ifdef SOLVE3D
        t1=nnew(ngp)
# else
        t1=knew(ngp)
# endif
      END IF
!
!  Set refinement scale.
!
      refscale=REFINED(ng)%ref_scale

!  Loop thru 3 times. 
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
        MySize=IJlen
        Npts=IJlen

        rectangular=.FALSE.
        IJspv=-9999.0_r8
        ncvarid=1  !not used

        IF (.not.allocated(Aout))  allocate(Aout(IJlen))
        IF (.not.allocated(A))     allocate(A(Imin:Imax,Jmin:Jmax))
        IF (.not.allocated(Amask)) allocate(Amask(Imin:Imax,Jmin:Jmax))
        IF (.not.allocated(A2))    allocate(A2(LBi:UBi,LBj:UBj))
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
              Amask(i,j)=REFINED(ng)%rmask_p(i,j)
# else
              Amask(i,j)=1.0_r8
# endif
            END DO
          END DO
!
! get zeta
!
          IF (iic(1).eq.0) THEN
            CALL mp_gather2d (ngp, model, LBip, UBip, LBjp, UBjp, 0,    &
     &                    gtype, 1.0_r8,                                &
# ifdef MASKING
     &                    GRID(ngp) % rmask(LBip:UBip,LBjp:UBjp),       &
# endif
     &                    OCEAN(ngp) % zeta(LBip:UBip,LBjp:UBjp,KOUTP), &
     &                    Npts, Aout)
          ELSE
            CALL mp_gather2d (ngp, model, LBip, UBip, LBjp, UBjp, 0,    &
     &                    gtype, 1.0_r8,                                &
# ifdef MASKING
     &                    GRID(ngp) % rmask(LBip:UBip,LBjp:UBjp),       &
# endif
# ifdef SOLVE3D
     &                    COUPLING(ngp) % Zt_avg1(LBip:UBip,LBjp:UBjp), &
# else
     &                    OCEAN(ngp) % zeta(LBip:UBip,LBjp:UBjp,KOUTP), &
# endif
     &                    Npts, Aout)
          END IF
!
! distribute the Aout (1D packed) field to all the nodes.
          CALL mpi_bcast (Aout, MySize, MP_FLOAT, MyMaster,             &
     &                    OCN_COMM_WORLD, MyError)
! unpack and place in global ng=1 2D array
          ic=0
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              ic=ic+1
              A(i,j)=Aout(ic)
            END DO
          END DO
!
! Now interpolate to the local tile.
!
# ifdef INTERP_BILINEAR
          CALL linterp2d_ref (ng, Imin, Imax, Jmin, Jmax,               &
     &                        A, Amask,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        OCEAN(ng)%hindspi_r, OCEAN(ng)%hindspj_r, &
     &                        A2, my_min, my_max)
# else
          CALL consrvterp2d (ng, Imin, Imax, Jmin, Jmax,              &
     &                       A,                                       &
#  ifdef MASKING
     &                       Amask,                                   &
#  endif 
     &                       LBi, UBi, LBj, UBj,                      &
     &                       LBi, UBi, LBj, UBj, refscale,            &
     &                       OCEAN(ng)%hindspi_r(i,j),                &
     &                       OCEAN(ng)%hindspj_r(i,j),                &
     &                       A2, my_min, my_max)
# endif
!
          IF (WESTERN_EDGE) THEN
            DO j=JstrTV,JendT
              DO i=IstrT,Istr
                ic=i+Nghostpoints+1
                REFINED(ng)%zetaR_west(ic,j,t1)=A2(i,j)
              END DO
            END DO
          END IF
          IF (EASTERN_EDGE) THEN
            DO j=JstrTV,JendT
              DO i=Iend,IendT
                ic=i-Iend+1
                REFINED(ng)%zetaR_east(ic,j,t1)=A2(i,j)
              END DO
            END DO
          END IF
          IF (SOUTHERN_EDGE) THEN
            DO j=JstrT,Jstr
              DO i=IstrTU,IendT
                jc=j+Nghostpoints+1
                REFINED(ng)%zetaR_south(i,jc,t1)=A2(i,j)
              END DO
            END DO
          END IF
          IF (NORTHERN_EDGE) THEN
            DO j=Jend,JendT
              DO i=IstrTU,IendT
                jc=j-Jend+1
                REFINED(ng)%zetaR_north(i,jc,t1)=A2(i,j)
              END DO
            END DO
          END IF
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
              Amask(i,j)=REFINED(ng)%umask_p(i,j)
# else
              Amask(i,j)=1.0_r8
# endif
            END DO
          END DO
!
! get ubar
!
          CALL mp_gather2d (ngp, model, LBip, UBip, LBjp, UBjp, 0,      &
     &                    gtype, 1.0_r8,                                &
# ifdef MASKING
     &                    GRID(ngp) % umask(LBip:UBip,LBjp:UBjp),       &
# endif
     &                OCEAN(ngp) % ubar(LBip:UBip,LBjp:UBjp,KOUTP),     &
     &                    Npts, Aout)
!
! distribute the Aout (1D packed) field to all the nodes for ng=2
          CALL mpi_bcast (Aout, MySize, MP_FLOAT, MyMaster,             &
     &                    OCN_COMM_WORLD, MyError)
! unpack and place in global ng=1 2D array
          ic=0
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              ic=ic+1
              A(i,j)=Aout(ic)
            END DO
          END DO
!
! Now interpolate to the local tile for ng=2.
!
# ifdef INTERP_BILINEAR
          CALL linterp2d_ref (ng, Imin, Imax, Jmin, Jmax,               &
     &                        A, Amask,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        OCEAN(ng)%hindspi_u, OCEAN(ng)%hindspj_u, &
     &                        A2, my_min, my_max)
# else
          CALL consrvterp2d (ng, Imin, Imax, Jmin, Jmax,                &
     &                       A,                                         &
#  ifdef MASKING
     &                       Amask,                                     &
#  endif 
     &                       LBi, UBi, LBj, UBj,                        &
     &                       LBi, UBi, LBj, UBj, refscale,              &
     &                       OCEAN(ng)%hindspi_u, OCEAN(ng)%hindspj_u,  &
     &                       A2, my_min, my_max)
# endif
!
          IF (WESTERN_EDGE) THEN
            DO j=JstrT,JendT
              DO i=IstrT+1,Istr
                ic=i+Nghostpoints
                REFINED(ng)%ubarR_west(ic,j,t1)=A2(i,j)
              END DO
            END DO
          END IF
          IF (EASTERN_EDGE) THEN
            DO j=JstrT,JendT
              DO i=Iend+1,IendT
                ic=i-Iend
                REFINED(ng)%ubarR_east(ic,j,t1)=A2(i,j)
              END DO
            END DO
          END IF
          IF (SOUTHERN_EDGE) THEN
            DO j=JstrT,Jstr-1
              DO i=IstrT,IendT
                jc=j+Nghostpoints+1
                REFINED(ng)%ubarR_south(i,jc,t1)=A2(i,j)
              END DO
            END DO
          END IF
          IF (NORTHERN_EDGE) THEN
            DO j=Jend+1,JendT
              DO i=IstrT,IendT
                jc=j-Jend
                REFINED(ng)%ubarR_north(i,jc,t1)=A2(i,j)
              END DO
            END DO
          END IF
!!!!!!!!!!!!!!!!!!!
!  DU_avg2  !!!!!!!
!!!!!!!!!!!!!!!!!!!
!
! get mask
!
          ic=0
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              ic=ic+1
# ifdef MASKING
              Amask(i,j)=REFINED(ng)%umask_p(i,j)
# else
              Amask(i,j)=1.0_r8
# endif
            END DO
          END DO
!
! get DU_avg2
!
          CALL mp_gather2d (ngp, model, LBip, UBip, LBjp, UBjp, 0,      &
     &                    gtype, 1.0_r8,                                &
# ifdef MASKING
     &                    GRID(ngp) % umask(LBip:UBip,LBjp:UBjp),       &
# endif
# ifdef SOLVE3D
     &                COUPLING(ngp) % DU_avg2(LBip:UBip,LBjp:UBjp),     &
# else
     &                REFINED(ngp) % DU_avg2(LBip:UBip,LBjp:UBjp),      &
# endif
     &                    Npts, Aout)
!
! distribute the Aout (1D packed) field to all the nodes for ng=2
          CALL mpi_bcast (Aout, MySize, MP_FLOAT, MyMaster,             &
     &                    OCN_COMM_WORLD, MyError)
! unpack and place in global ng=1 2D array
          ic=0
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              ic=ic+1
              A(i,j)=Aout(ic)
            END DO
          END DO
!
! Vary the child depth avg flux within each parent block. Use a 
! linear variation to ensure the child sum equals the parent.
!
          IF (WESTERN_EDGE) THEN
            DO j=Jstr,Jend
              ic=NINT(OCEAN(ng)%hindpi_u(1,j))
              jc=NINT(OCEAN(ng)%hindpj_u(1,j))
              IF ((jc.ne.-9999).and.(ic.ne.-9999)) THEN
                REFINED(ng)%DUavg2R_west(j)=A(ic,jc)*                    &
     &                                      GRID(ng)%on_u(1,j)/          &
     &                                      REFINED(ng)%on_u_p(ic,jc)
              END IF
            END DO
          END IF
          IF (EASTERN_EDGE) THEN
            DO j=Jstr,Jend
              ic=NINT(OCEAN(ng)%hindpi_u(Iend+1,j))
              jc=NINT(OCEAN(ng)%hindpj_u(Iend+1,j))
              IF ((jc.ne.-9999).and.(ic.ne.-9999)) THEN
                REFINED(ng)%DUavg2R_east(j)=A(ic,jc)*                    &
     &                                      GRID(ng)%on_u(Iend+1,j)/     &
     &                                      REFINED(ng)%on_u_p(ic,jc)
              END IF
            END DO
          END IF
# ifdef WEC
!!!!!!!!!!!!!!!!!!!
!  ubar_stokes !!!!
!!!!!!!!!!!!!!!!!!!
!
! get ubar_stokes
!
          CALL mp_gather2d (ngp, model, LBip, UBip, LBjp, UBjp, 0,      &
     &                    gtype, 1.0_r8,                                &
#  ifdef MASKING
     &                    GRID(ngp) % umask(LBip:UBip,LBjp:UBjp),       &
#  endif
     &               OCEAN(ngp) % ubar_stokes(LBip:UBip,LBjp:UBjp),     &
     &                    Npts, Aout)
!
! distribute the Aout (1D packed) field to all the nodes for ng=2
          CALL mpi_bcast (Aout, MySize, MP_FLOAT, MyMaster,             &
     &                    OCN_COMM_WORLD, MyError)
! unpack and place in global ng=1 2D array
          ic=0
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              ic=ic+1
              A(i,j)=Aout(ic)
            END DO
          END DO
!
! Now interpolate to the local tile for ng=2.
!
# ifdef INTERP_BILINEAR
          CALL linterp2d_ref (ng, Imin, Imax, Jmin, Jmax,               &
     &                        A, Amask,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        OCEAN(ng)%hindspi_u, OCEAN(ng)%hindspj_u, &
     &                        A2, my_min, my_max)
# else
          CALL consrvterp2d (ng, Imin, Imax, Jmin, Jmax,                &
     &                       A,                                         &
#  ifdef MASKING
     &                       Amask,                                     &
#  endif 
     &                       LBi, UBi, LBj, UBj,                        &
     &                       LBi, UBi, LBj, UBj, refscale,              &
     &                       OCEAN(ng)%hindspi_u, OCEAN(ng)%hindspj_u,  &
     &                       A2, my_min, my_max)
# endif
!
          IF (WESTERN_EDGE) THEN
            DO j=JstrT,JendT
              DO i=IstrT+1,Istr
                ic=i+Nghostpoints
                REFINED(ng)%ubarsR_west(ic,j,t1)=A2(i,j)
              END DO
            END DO
          END IF
          IF (EASTERN_EDGE) THEN
            DO j=JstrT,JendT
              DO i=Iend+1,IendT
                ic=i-Iend
                REFINED(ng)%ubarsR_east(ic,j,t1)=A2(i,j)
              END DO
            END DO
          END IF
          IF (SOUTHERN_EDGE) THEN
            DO j=JstrT,Jstr-1
              DO i=IstrT,IendT
                jc=j+Nghostpoints+1
                REFINED(ng)%ubarsR_south(i,jc,t1)=A2(i,j)
              END DO
            END DO
          END IF
          IF (NORTHERN_EDGE) THEN
            DO j=Jend+1,JendT
              DO i=IstrT,IendT
                jc=j-Jend
                REFINED(ng)%ubarsR_north(i,jc,t1)=A2(i,j)
              END DO
            END DO
          END IF
# endif

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
              Amask(i,j)=REFINED(ng)%vmask_p(i,j)
# else
              Amask(i,j)=1.0_r8
# endif
            END DO
          END DO
!
! get vbar
!
          CALL mp_gather2d (ngp, model, LBip, UBip, LBjp, UBjp, 0,      &
     &                    gtype, 1.0_r8,                                &
# ifdef MASKING
     &                    GRID(ngp) % vmask(LBip:UBip,LBjp:UBjp),       &
# endif
     &                OCEAN(ngp) % vbar(LBip:UBip,LBjp:UBjp,KOUTP),     &
     &                    Npts, Aout)
!
! distribute the Aout (1D packed) field to all the nodes for ng=2
          CALL mpi_bcast (Aout, MySize, MP_FLOAT, MyMaster,             &
     &                    OCN_COMM_WORLD, MyError)
! unpack and place in global ng=1 2D array
          ic=0
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              ic=ic+1
              A(i,j)=Aout(ic)
            END DO
          END DO
!
! Now interpolate to the local tile for ng=2.
!
# ifdef INTERP_BILINEAR
          CALL linterp2d_ref (ng, Imin, Imax, Jmin, Jmax,               &
     &                        A, Amask,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        OCEAN(ng)%hindspi_v, OCEAN(ng)%hindspj_v, &
     &                        A2, my_min, my_max)
# else
          CALL consrvterp2d (ng, Imin, Imax, Jmin, Jmax,                &
     &                       A,                                         &
#  ifdef MASKING
     &                       Amask,                                     &
#  endif 
     &                       LBi, UBi, LBj, UBj,                        &
     &                       LBi, UBi, LBj, UBj, refscale,              &
     &                       OCEAN(ng)%hindspi_v, OCEAN(ng)%hindspj_v,  &
     &                       A2, my_min, my_max)
# endif
          IF (WESTERN_EDGE) THEN
            DO j=JstrT,JendT
              DO i=IstrT,Istr-1
                ic=i+Nghostpoints+1
                REFINED(ng)%vbarR_west(ic,j,t1)=A2(i,j)
              END DO
            END DO
          END IF
          IF (EASTERN_EDGE) THEN
            DO j=JstrT,JendT
              DO i=Iend+1,IendT
                ic=i-Iend
                REFINED(ng)%vbarR_east(ic,j,t1)=A2(i,j)
              END DO
            END DO
          END IF
          IF (SOUTHERN_EDGE) THEN
            DO j=JstrT+1,Jstr
              DO i=IstrT,IendT
                jc=j+Nghostpoints
                REFINED(ng)%vbarR_south(i,jc,t1)=A2(i,j)
              END DO
            END DO
          END IF
          IF (NORTHERN_EDGE) THEN
            DO j=Jend+1,JendT
              DO i=IstrT,IendT
                jc=j-Jend
                REFINED(ng)%vbarR_north(i,jc,t1)=A2(i,j)
              END DO
            END DO
          END IF
!!!!!!!!!!!!!!!!!!!
!  DV_avg2  !!!!!!!
!!!!!!!!!!!!!!!!!!!
!
! get mask
!
          ic=0
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              ic=ic+1
# ifdef MASKING
              Amask(i,j)=REFINED(ng)%vmask_p(i,j)
# else
              Amask(i,j)=1.0_r8
# endif
            END DO
          END DO
!
! get vbar
!
          CALL mp_gather2d (ngp, model, LBip, UBip, LBjp, UBjp, 0,      &
     &                    gtype, 1.0_r8,                                &
# ifdef MASKING
     &                    GRID(ngp) % vmask(LBip:UBip,LBjp:UBjp),       &
# endif
# ifdef SOLVE3D
     &                COUPLING(ngp) % DV_avg2(LBip:UBip,LBjp:UBjp),     &
# else
     &                REFINED(ngp) % DV_avg2(LBip:UBip,LBjp:UBjp),      &
# endif
     &                    Npts, Aout)
!
! distribute the Aout (1D packed) field to all the nodes for ng>2
          CALL mpi_bcast (Aout, MySize, MP_FLOAT, MyMaster,             &
     &                    OCN_COMM_WORLD, MyError)
! unpack and place in global ng=1 2D array
          ic=0
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              ic=ic+1
              A(i,j)=Aout(ic)
            END DO
          END DO
          IF (SOUTHERN_EDGE) THEN
            DO i=Istr,Iend
              ic=NINT(OCEAN(ng)%hindpi_v(i,1))
              jc=NINT(OCEAN(ng)%hindpj_v(i,1))
              IF ((jc.ne.-9999).and.(ic.ne.-9999)) THEN
                REFINED(ng)%DVavg2R_south(i)=A(ic,jc)*                   &
     &                                       GRID(ng)%om_v(i,1)/         &
     &                                       REFINED(ng)%om_v_p(ic,jc)
              END IF
            END DO
          END IF
          IF (NORTHERN_EDGE) THEN
            DO i=Istr,Iend
              ic=NINT(OCEAN(ng)%hindpi_v(i,Jend+1))
              jc=NINT(OCEAN(ng)%hindpj_v(i,Jend+1))
              IF ((jc.ne.-9999).and.(ic.ne.-9999)) THEN
                REFINED(ng)%DVavg2R_north(i)=A(ic,jc)*                   &
     &                                       GRID(ng)%om_v(i,Jend+1)/    &
     &                                       REFINED(ng)%om_v_p(ic,jc)
              END IF
            END DO
          END IF
# ifdef WEC
!!!!!!!!!!!!!!!!!!!
!  vbar_stokes !!!!
!!!!!!!!!!!!!!!!!!!
!
! get vbar_stokes
!
          CALL mp_gather2d (ngp, model, LBip, UBip, LBjp, UBjp, 0,      &
     &                    gtype, 1.0_r8,                                &
#  ifdef MASKING
     &                    GRID(ngp) % vmask(LBip:UBip,LBjp:UBjp),       &
#  endif
     &               OCEAN(ngp) % vbar_stokes(LBip:UBip,LBjp:UBjp),     &
     &                    Npts, Aout)
!
! distribute the Aout (1D packed) field to all the nodes for ng=2
          CALL mpi_bcast (Aout, MySize, MP_FLOAT, MyMaster,             &
     &                    OCN_COMM_WORLD, MyError)
! unpack and place in global ng=1 2D array
          ic=0
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              ic=ic+1
              A(i,j)=Aout(ic)
            END DO
          END DO
!
! Now interpolate to the local tile for ng=2.
!
# ifdef INTERP_BILINEAR
          CALL linterp2d_ref (ng, Imin, Imax, Jmin, Jmax,               &
     &                        A, Amask,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        OCEAN(ng)%hindspi_v, OCEAN(ng)%hindspj_v, &
     &                        A2, my_min, my_max)
# else
          CALL consrvterp2d (ng, Imin, Imax, Jmin, Jmax,                &
     &                       A,                                         &
#  ifdef MASKING
     &                       Amask,                                     &
#  endif 
     &                       LBi, UBi, LBj, UBj,                        &
     &                       LBi, UBi, LBj, UBj, refscale,              &
     &                       OCEAN(ng)%hindspi_v, OCEAN(ng)%hindspj_v,  &
     &                       A2, my_min, my_max)
# endif
          IF (WESTERN_EDGE) THEN
            DO j=JstrT,JendT
              DO i=IstrT,Istr-1
                ic=i+Nghostpoints+1
                REFINED(ng)%vbarsR_west(ic,j,t1)=A2(i,j)
              END DO
            END DO
          END IF
          IF (EASTERN_EDGE) THEN
            DO j=JstrT,JendT
              DO i=Iend+1,IendT
                ic=i-Iend
                REFINED(ng)%vbarsR_east(ic,j,t1)=A2(i,j)
              END DO
            END DO
          END IF
          IF (SOUTHERN_EDGE) THEN
            DO j=JstrT+1,Jstr
              DO i=IstrT,IendT
                jc=j+Nghostpoints
                REFINED(ng)%vbarsR_south(i,jc,t1)=A2(i,j)
              END DO
            END DO
          END IF
          IF (NORTHERN_EDGE) THEN
            DO j=Jend+1,JendT
              DO i=IstrT,IendT
                jc=j-Jend
                REFINED(ng)%vbarsR_north(i,jc,t1)=A2(i,j)
              END DO
            END DO
          END IF
# endif

        END IF
!
        IF (allocated(Aout))  deallocate (Aout)
        IF (allocated(A))     deallocate (A)
        IF (allocated(Amask)) deallocate (Amask)
        IF (allocated(A2))    deallocate (A2)
      END DO
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
      END SUBROUTINE get_2dparent_data_tile
#endif
      END MODULE get_2dparent_data_mod
