#include "cppdefs.h"
      MODULE set_2dchild_data_mod
#ifdef REFINED_GRID
!
!svn $Id: set_data.F 1405 2008-02-20 01:44:07Z jcwarner $
!================================================== John C. Warner =====
!                                                                      !
!  This routine acquires data from a child grid.                       !
!  That data is averaged and placed into the parent grid.              !
!  Need special treatment to use Amask for case of                     !
!  wet/dry to ensure dry cells are not used in the avg.                !
!  Also added ability for sum(child_cell_areas) may not be equal       !
!  to the parent area.                                                 !
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
      USE mp_exchange_mod, ONLY : mp_exchange2d
# endif
# if defined EW_PERIODIC || defined NS_PERIODIC
      USE exchange_2d_mod
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
      integer :: i, j, ngc, ic, jc, rs
      integer :: ii, ir, jr, gtype, Npts
      integer :: Imin, Imax, Jmin, Jmax, model
      integer :: Ilen, Jlen, IJlen, MyType, status
      integer :: MyError, MySizer, ncvarid
      integer :: Xmin, Xmax, Ymin, Ymax
      integer :: Istrp, Iendp, Jstrp, Jendp
      real(r8) :: my_val, fac
      real(r8) :: my_min, my_max
!
      real(r8), allocatable :: Aout(:)
      real(r8), allocatable :: A(:,:)
      real(r8), allocatable :: Amask(:,:)

      logical :: rectangular
      real(r8) :: cff, cff1, IJspv
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
        IF (.not.allocated(Amask)) THEN
          allocate(Amask(Imin:Imax,Jmin:Jmax))
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
# ifdef WET_DRY
          CALL mp_gather2d (ngc, model, LBic, UBic, LBjc, UBjc, 0,      &
     &               gtype, 1.0_r8,                                     &
     &               GRID(ngc) % rmask(LBic:UBic,LBjc:UBjc),            &
     &               GRID(ngc) % rmask_full(LBic:UBic,LBjc:UBjc),       &
     &               Npts, Aout)
!
! distribute the Aout (1D packed) field to all the nodes for ng=2
!
          CALL mpi_bcast (Aout, MySizer, MP_FLOAT, MyMaster,            &
     &                    OCN_COMM_WORLD, MyError)
!
! unpack and place in global ng=1 2D array
!
          ic=0
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              ic=ic+1
              Amask(i,j)=Aout(ic)
            END DO
          END DO
# endif
# if defined MASKING && !defined WET_DRY
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              Amask(i,j)=REFINED(ng)%rmask_c(i,j)
            END DO
          END DO
# endif
!
! get zeta
!
          CALL mp_gather2d (ngc, model, LBic, UBic, LBjc, UBjc, 0,      &
     &               gtype, 1.0_r8,                                     &
# ifdef MASKING
     &               GRID(ngc) % rmask(LBic:UBic,LBjc:UBjc),            &
# endif
# ifdef SOLVE3D
     &               COUPLING(ngc) % Zt_avg1(LBic:UBic,LBjc:UBjc),      &
# else
     &               OCEAN(ngc) % zeta(LBic:UBic,LBjc:UBjc,KOUTC),      &
# endif
     &               Npts, Aout)
!
! distribute the Aout (1D packed) field to all the nodes for ng=2
!
          CALL mpi_bcast (Aout, MySizer, MP_FLOAT, MyMaster,            &
     &                    OCN_COMM_WORLD, MyError)
!
! unpack and place in global ng=1 2D array
!
          ic=0
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              ic=ic+1
              A(i,j)=Aout(ic)
            END DO
          END DO
!
!  Set limits of the parent grid to fill with child data.
!
          Xmin=OCEAN(ng)%hindcr(1)
          Xmax=OCEAN(ng)%hindcr(2)
          Ymin=OCEAN(ng)%hindcr(3)
          Ymax=OCEAN(ng)%hindcr(4)
          IF (Xmax.gt.-9999) THEN
!
!  Average the child data back to parent resolution.
!
            DO j=Ymin,Ymax
              DO i=Xmin,Xmax
                ic=NINT(OCEAN(ng)%hindci_r(i,j))
                jc=NINT(OCEAN(ng)%hindcj_r(i,j))
                rs=(REFINED(ngc)%ref_scale-1)/2
                my_val=0.0_r8
                cff=0.0_r8
                DO ir=-rs,rs
                  DO jr=-rs,rs
                    cff1=A(ic+ir,jc+jr)*                                &
     &                   REFINED(ng)%om_r_c(ic+ir,jc+jr)*               &
     &                   REFINED(ng)%on_r_c(ic+ir,jc+jr)
# ifdef MASKING
                    my_val=my_val+cff1*Amask(ic+ir,jc+jr)
                    cff=cff+1.0_r8*Amask(ic+ir,jc+jr)
# else
                    my_val=my_val+cff1
# endif
                  END DO
                END DO
                my_val=my_val*GRID(ng)%pm(i,j)*GRID(ng)%pn(i,j)
# ifdef MASKING
                IF (cff.gt.0.0_r8) THEN
                  my_val=my_val*REAL(REFINED(ngc)%ref_scale,r8)*        &
     &                          REAL(REFINED(ngc)%ref_scale,r8)/cff
                END IF
                my_val=my_val*GRID(ng)%rmask(i,j)
# endif
# ifdef WET_DRY
                IF (my_val.le.(Dcrit(ng)-GRID(ng)%h(i,j))) THEN
                  my_val=Dcrit(ng)-GRID(ng)%h(i,j)
                END IF
# endif
# ifdef SOLVE3D
                COUPLING(ng)%Zt_avg1(i,j)=my_val
# else
                OCEAN(ng)%zeta(i,j,1)=my_val
                OCEAN(ng)%zeta(i,j,2)=my_val
# endif
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
# ifdef MASKING
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              Amask(i,j)=REFINED(ng)%umask_c(i,j)
            END DO
          END DO
# endif
!
! get ubar
!
          CALL mp_gather2d (ngc, model, LBic, UBic, LBjc, UBjc, 0,      &
     &              gtype, 1.0_r8,                                      &
# ifdef MASKING
     &              GRID(ngc) % umask(LBic:UBic,LBjc:UBjc),             &
# endif
     &              OCEAN(ngc) % ubar(LBic:UBic,LBjc:UBjc,KOUTC),       &
     &              Npts, Aout)
!
! distribute the Aout (1D packed) field to all the nodes
!
          CALL mpi_bcast (Aout, MySizer, MP_FLOAT, MyMaster,            &
     &                    OCN_COMM_WORLD, MyError)
!
! unpack and place in global ng=1 2D array
!
          ic=0
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              ic=ic+1
              A(i,j)=Aout(ic)
            END DO
          END DO
!
!  Set limits of the parent grid to fill with child data.
!
          Xmin=OCEAN(ng)%hindcu(1)
          Xmax=OCEAN(ng)%hindcu(2)
          Ymin=OCEAN(ng)%hindcu(3)
          Ymax=OCEAN(ng)%hindcu(4)
          IF (Xmax.gt.-9999) THEN
!
!  Average the child data back to parent resolution.
!
            DO j=Ymin,Ymax
              DO i=Xmin,Xmax
                ic=NINT(OCEAN(ng)%hindci_u(i,j))
                jc=NINT(OCEAN(ng)%hindcj_u(i,j))
                rs=(REFINED(ngc)%ref_scale-1)/2
                my_val=0.0_r8
                cff=0.0_r8
                DO ir=-rs,rs
                  DO jr=-rs,rs
                    cff1=A(ic+ir,jc+jr)
# ifdef MASKING
                    my_val=my_val+cff1*Amask(ic+ir,jc+jr)
                    cff=cff+1.0_r8*Amask(ic+ir,jc+jr)
# else
                    my_val=my_val+cff1
                    cff=cff+1.0_r8
# endif
                  END DO
                END DO
                IF (cff.gt.0.0_r8) my_val=my_val/cff
# ifdef MASKING
                my_val=my_val*GRID(ng)%umask(i,j)
# endif
                OCEAN(ng)%ubar(i,j,1)=my_val
                OCEAN(ng)%ubar(i,j,2)=my_val
              END DO
            END DO
          END IF
        ELSE IF (ii.eq.3) THEN
!!!!!!!!!!!!!!!!!!!
!  vbar  !!!!!!!!!!
!!!!!!!!!!!!!!!!!!!
!
! get mask
!
# ifdef MASKING
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              Amask(i,j)=REFINED(ng)%vmask_c(i,j)
            END DO
          END DO
# endif
!
! get vbar
!
          CALL mp_gather2d (ngc, model, LBic, UBic, LBjc, UBjc, 0,      &
     &              gtype, 1.0_r8,                                      &
# ifdef MASKING
     &              GRID(ngc) % vmask(LBic:UBic,LBjc:UBjc),             &
# endif
     &              OCEAN(ngc) % vbar(LBic:UBic,LBjc:UBjc,KOUTC),       &
     &              Npts, Aout)
!
! distribute the Aout (1D packed) field to all the nodes for ng=2
!
          CALL mpi_bcast (Aout, MySizer, MP_FLOAT, MyMaster,            &
     &                    OCN_COMM_WORLD, MyError)
!
! unpack and place in global ng=1 2D array
!
          ic=0
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              ic=ic+1
              A(i,j)=Aout(ic)
            END DO
          END DO
!
!  Set limits of the parent grid to fill with child data.
!
          Xmin=OCEAN(ng)%hindcv(1)
          Xmax=OCEAN(ng)%hindcv(2)
          Ymin=OCEAN(ng)%hindcv(3)
          Ymax=OCEAN(ng)%hindcv(4)
          IF (Xmax.gt.-9999) THEN
!
!  Average the child data back to parent resolution.
!
            DO j=Ymin,Ymax
              DO i=Xmin,Xmax
                ic=NINT(OCEAN(ng)%hindci_v(i,j))
                jc=NINT(OCEAN(ng)%hindcj_v(i,j))
                rs=(REFINED(ngc)%ref_scale-1)/2
                my_val=0.0_r8
                cff=0.0_r8
                DO ir=-rs,rs
                  DO jr=-rs,rs
                    cff1=A(ic+ir,jc+jr)
# ifdef MASKING
                    my_val=my_val+cff1*Amask(ic+ir,jc+jr)
                    cff=cff+1.0_r8*Amask(ic+ir,jc+jr)
# else
                    my_val=my_val+cff1
                    cff=cff+1.0_r8
# endif
                  END DO
                END DO
                IF (cff.gt.0.0_r8) my_val=my_val/cff
# ifdef MASKING
                my_val=my_val*GRID(ng)%vmask(i,j)
# endif
                OCEAN(ng)%vbar(i,j,1)=my_val
                OCEAN(ng)%vbar(i,j,2)=my_val
              END DO
            END DO
          END IF
        END IF
!
        IF (allocated(Aout)) THEN
          deallocate (Aout)
        END IF
        IF (allocated(A)) THEN
          deallocate (A)
        END IF
        IF (allocated(Amask)) THEN
          deallocate (Amask)
        END IF
      END DO
# if defined EW_PERIODIC || defined NS_PERIODIC
!
!-----------------------------------------------------------------------
!  Exchange boundary information.
!-----------------------------------------------------------------------
!
#   ifdef SOLVE3D
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        COUPLING(ng)%Zt_avg1(:,:))
#   else
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        OCEAN(ng)%zeta(:,:,1))
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        OCEAN(ng)%zeta(:,:,2))
#   endif
      CALL exchange_u2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        OCEAN(ng)%ubar(:,:,1))
      CALL exchange_u2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        OCEAN(ng)%ubar(:,:,2))
      CALL exchange_v2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        OCEAN(ng)%vbar(:,:,1))
      CALL exchange_v2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        OCEAN(ng)%vbar(:,:,2))
# endif
# ifdef DISTRIBUTE
#  ifdef SOLVE3D
      CALL mp_exchange2d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    COUPLING(ng)%Zt_avg1(:,:))
#  else
      CALL mp_exchange2d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    OCEAN(ng)%zeta(:,:,1),                        &
     &                    OCEAN(ng)%zeta(:,:,2))
#  endif
      CALL mp_exchange2d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    OCEAN(ng)%ubar(:,:,1),                        &
     &                    OCEAN(ng)%ubar(:,:,2))
      CALL mp_exchange2d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    OCEAN(ng)%vbar(:,:,1),                        &
     &                    OCEAN(ng)%vbar(:,:,2))
# endif
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
