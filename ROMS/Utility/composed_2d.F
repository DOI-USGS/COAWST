#include "cppdefs.h"
      MODULE composed_2d_mod
#ifdef COMPOSED_GRID
!=======================================================================
!                                                                      !
!  This routine allocates and exchanges data for composed grids.       !
!                                                                      !
!=======================================================================

      implicit none

      PRIVATE
      PUBLIC  :: composed_2d

      CONTAINS
!
!***********************************************************************
      SUBROUTINE composed_2d (tile, inflag)
!***********************************************************************
!
      USE mod_param
# ifdef SOLVE3D
      USE mod_coupling
# endif
      USE mod_forces
      USE mod_ncparam
      USE mod_scalars
      USE mod_stepping
      USE mod_composed
      USE mod_grid
      USE mod_ocean
      USE mod_parallel
# ifdef DISTRIBUTE
      USE distribute_mod, ONLY : mp_bcasti, mp_gather2d
      USE mp_exchange_mod, ONLY : mp_exchange2d
# endif

      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: tile, inflag
!
!  Local variable declarations.
!
      integer :: ng, LBi, UBi, LBj, UBj
      integer :: id, gtype, Tindex, m, itrc, idtim, ntim
      integer :: i, j, ic, jc, ngg, ngc
      integer :: LBic, UBic, LBjc, UBjc
      integer :: Imin, Imax, Jmin, Jmax, model
      integer :: Ilen, Jlen, IJlen, Npts
      integer :: MyError
      integer :: Xmin, Xmax, Ymin, Ymax

      real(r8) :: my_val, IJspv
      real(r8), allocatable :: Aout(:)
      real(r8), allocatable :: A(:,:)

      logical :: dointerp

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

      model=iNLM
      IJspv=-9999.0_r8



! inflag = 100 --> we are calling from main2d step2d.
      IF (inflag.eq.100) THEN
        DO m=1,4
          dointerp=.TRUE.
!          idtimm=knew(1)
          IF (m.eq.1) id=idFsur  !zeta
          IF (m.eq.2) id=idUBar  !ubar
          IF (m.eq.3) id=idVbar  !vbar
          IF (m.eq.4) THEN
            IF (PREDICTOR_2D_STEP(1)) THEN
              id=idRzet          !rzeta
!              idtimm=krhs(1)
            ELSE
              dointerp=.FALSE.
            END IF
          END IF
          gtype=Iinfo(1,id,1)
          IF (dointerp) THEN
!            CALL allocate_nesting2d (field2d, id, gtype, idtimm)
!            CALL interp_nesting2d (field2d)
          END IF
        END DO
      ENDIF


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!         inflag = 9 : bustr
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF (inflag.eq.9) THEN
!        ntim=1
!        DO m=1,ntim
          DO ng=1,Ngrids
            idtim=m
            LBi=BOUNDS(ng)%LBi(tile)
            UBi=BOUNDS(ng)%UBi(tile)
            LBj=BOUNDS(ng)%LBj(tile)
            UBj=BOUNDS(ng)%UBj(tile)
            DO ngg=1,Ngrids-1
              ngc=COMPOSED(ng)%mychild(ngg)
!
              LBic=BOUNDS(ngc)%LBi(tile)
              UBic=BOUNDS(ngc)%UBi(tile)
              LBjc=BOUNDS(ngc)%LBj(tile)
              UBjc=BOUNDS(ngc)%UBj(tile)
!
              gtype=u2dvar
              Imin=LOWER_BOUND_I
              Imax=Lm(ngc)+GHOST_POINTS
              Jmin=LOWER_BOUND_J
              Jmax=Mm(ngc)+GHOST_POINTS
              IF ((ABS(gtype).eq.p2dvar).or.(ABS(gtype).eq.u2dvar)) THEN
                Imin=LOWER_BOUND_I+1
              END IF
              IF ((ABS(gtype).eq.p2dvar).or.(ABS(gtype).eq.v2dvar)) THEN
                Jmin=LOWER_BOUND_J+1
              END IF
              Ilen=Imax-Imin+1
              Jlen=Jmax-Jmin+1
              IJlen=Ilen*Jlen
              Npts=IJlen

              IF (.not.allocated(Aout)) THEN
                allocate (Aout(IJlen))
              END IF
              IF (.not.allocated(A)) THEN
                allocate(A(Imin:Imax,Jmin:Jmax))
              END IF
!
! get ubar
!
              CALL mp_gather2d (ngc, model, LBic, UBic, LBjc, UBjc, 0,  &
     &                          gtype, 1.0_r8,                          &
# ifdef MASKING
     &               GRID(ngc) % umask(LBic:UBic,LBjc:UBjc),            &
# endif
     &               FORCES(ngc) % bustr(LBic:UBic,LBjc:UBjc),          &
     &               Npts, Aout)
!
! distribute the Aout (1D packed) field to all the nodes
!
              CALL mpi_bcast (Aout, IJlen, MP_FLOAT, MyMaster,          &
     &                        OCN_COMM_WORLD, MyError)
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
              Xmin=OCEAN(ng)%hindu(1,ngg)
              Xmax=OCEAN(ng)%hindu(2,ngg)
              Ymin=OCEAN(ng)%hindu(3,ngg)
              Ymax=OCEAN(ng)%hindu(4,ngg)
              IF (Xmax.gt.-9999) THEN
!
! Set the chld data into the parent locations. A is on the child grid.
! Use hindices to obtain the parent location on the child grid.
!
                DO j=Ymin,Ymax
                  DO i=Xmin,Xmax
                    ic=NINT(OCEAN(ng)%hindi_u(i,j))
                    jc=NINT(OCEAN(ng)%hindj_u(i,j))
                    IF ((ic.ne.IJspv).and.(jc.ne.IJspv)) THEN
                      my_val=A(ic,jc)
# ifdef MASKING
                      my_val=my_val*GRID(ng)%umask(i,j)
# endif
                      FORCES(ng)%bustr(i,j)=my_val
                    END IF
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
            END DO
# ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                            &
     &                          LBi, UBi, LBj, UBj,                           &
     &                          NghostPoints, EWperiodic, NSperiodic,         &
     &                          FORCES(ng)%bustr(:,:))
# endif
          END DO
!       END DO
      ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!         inflag = 9 : bvstr
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF (inflag.eq.9) THEN
!       ntim=1
!       DO m=1,ntim
          DO ng=1,Ngrids
            idtim=m
            LBi=BOUNDS(ng)%LBi(tile)
            UBi=BOUNDS(ng)%UBi(tile)
            LBj=BOUNDS(ng)%LBj(tile)
            UBj=BOUNDS(ng)%UBj(tile)
            DO ngg=1,Ngrids-1
              ngc=COMPOSED(ng)%mychild(ngg)
!
              LBic=BOUNDS(ngc)%LBi(tile)
              UBic=BOUNDS(ngc)%UBi(tile)
              LBjc=BOUNDS(ngc)%LBj(tile)
              UBjc=BOUNDS(ngc)%UBj(tile)
!
              gtype=v2dvar
              Imin=LOWER_BOUND_I
              Imax=Lm(ngc)+GHOST_POINTS
              Jmin=LOWER_BOUND_J
              Jmax=Mm(ngc)+GHOST_POINTS
              IF ((ABS(gtype).eq.p2dvar).or.(ABS(gtype).eq.u2dvar)) THEN
                Imin=LOWER_BOUND_I+1
              END IF
              IF ((ABS(gtype).eq.p2dvar).or.(ABS(gtype).eq.v2dvar)) THEN
                Jmin=LOWER_BOUND_J+1
              END IF
              Ilen=Imax-Imin+1
              Jlen=Jmax-Jmin+1
              IJlen=Ilen*Jlen
              Npts=IJlen

              IF (.not.allocated(Aout)) THEN
                allocate (Aout(IJlen))
              END IF
              IF (.not.allocated(A)) THEN
                allocate(A(Imin:Imax,Jmin:Jmax))
              END IF
!
! get ubar
!
              CALL mp_gather2d (ngc, model, LBic, UBic, LBjc, UBjc, 0,  &
     &                          gtype, 1.0_r8,                          &
# ifdef MASKING
     &               GRID(ngc) % vmask(LBic:UBic,LBjc:UBjc),            &
# endif
     &               FORCES(ngc) % bvstr(LBic:UBic,LBjc:UBjc),          &
     &               Npts, Aout)
!
! distribute the Aout (1D packed) field to all the nodes
!
              CALL mpi_bcast (Aout, IJlen, MP_FLOAT, MyMaster,          &
     &                        OCN_COMM_WORLD, MyError)
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
              Xmin=OCEAN(ng)%hindv(1,ngg)
              Xmax=OCEAN(ng)%hindv(2,ngg)
              Ymin=OCEAN(ng)%hindv(3,ngg)
              Ymax=OCEAN(ng)%hindv(4,ngg)
              IF (Xmax.gt.-9999) THEN
!
! Set the chld data into the parent locations. A is on the child grid.
! Use hindices to obtain the parent location on the child grid.
!
                DO j=Ymin,Ymax
                  DO i=Xmin,Xmax
                    ic=NINT(OCEAN(ng)%hindi_v(i,j))
                    jc=NINT(OCEAN(ng)%hindj_v(i,j))
                    IF ((ic.ne.IJspv).and.(jc.ne.IJspv)) THEN
                      my_val=A(ic,jc)
# ifdef MASKING
                      my_val=my_val*GRID(ng)%vmask(i,j)
# endif
                      FORCES(ng)%bvstr(i,j)=my_val
                    END IF
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
            END DO
# ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                            &
     &                          LBi, UBi, LBj, UBj,                           &
     &                          NghostPoints, EWperiodic, NSperiodic,         &
     &                          FORCES(ng)%bvstr(:,:))
# endif
          END DO
!       END DO
      ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!         inflag = 11, 14, or 15 : zeta
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF ((inflag.eq.11).or.(inflag.eq.14).or.(inflag.eq.15)) THEN
        IF (inflag.eq.11) THEN
          ntim=2
        ELSE
          ntim=1
        END IF
        DO m=1,ntim
          DO ng=1,Ngrids
            IF (inflag.eq.11) THEN
              idtim=m
            ELSE
              idtim=knew(ng)
            END IF
            LBi=BOUNDS(ng)%LBi(tile)
            UBi=BOUNDS(ng)%UBi(tile)
            LBj=BOUNDS(ng)%LBj(tile)
            UBj=BOUNDS(ng)%UBj(tile)
            DO ngg=1,Ngrids-1
              ngc=COMPOSED(ng)%mychild(ngg)
!
              LBic=BOUNDS(ngc)%LBi(tile)
              UBic=BOUNDS(ngc)%UBi(tile)
              LBjc=BOUNDS(ngc)%LBj(tile)
              UBjc=BOUNDS(ngc)%UBj(tile)

              gtype=r2dvar
              Imin=LOWER_BOUND_I
              Imax=Lm(ngc)+GHOST_POINTS
              Jmin=LOWER_BOUND_J
              Jmax=Mm(ngc)+GHOST_POINTS
              IF ((ABS(gtype).eq.p2dvar).or.(ABS(gtype).eq.u2dvar)) THEN
                Imin=LOWER_BOUND_I+1
              END IF
              IF ((ABS(gtype).eq.p2dvar).or.(ABS(gtype).eq.v2dvar)) THEN
                Jmin=LOWER_BOUND_J+1
              END IF
              Ilen=Imax-Imin+1
              Jlen=Jmax-Jmin+1
              IJlen=Ilen*Jlen
              Npts=IJlen

              IF (.not.allocated(Aout)) THEN
                allocate (Aout(IJlen))
              END IF
              IF (.not.allocated(A)) THEN
                allocate(A(Imin:Imax,Jmin:Jmax))
              END IF
!
! get zeta
!
              CALL mp_gather2d (ngc, model, LBic, UBic, LBjc, UBjc, 0,  &
     &                          gtype, 1.0_r8,                          &
# ifdef MASKING
     &               GRID(ngc) % rmask(LBic:UBic,LBjc:UBjc),            &
# endif
     &               OCEAN(ngc) % zeta(LBic:UBic,LBjc:UBjc,idtim),      &
     &               Npts, Aout)
!
! distribute the Aout (1D packed) field to all the nodes for ng=2
!
              CALL mpi_bcast (Aout, IJlen, MP_FLOAT, MyMaster,          &
     &                        OCN_COMM_WORLD, MyError)
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
              Xmin=OCEAN(ng)%hindr(1,ngg)
              Xmax=OCEAN(ng)%hindr(2,ngg)
              Ymin=OCEAN(ng)%hindr(3,ngg)
              Ymax=OCEAN(ng)%hindr(4,ngg)
              IF (Xmax.gt.-9999) THEN
!
! Set the chld data into the parent locations. A is on the child grid.
! Use hindices to obtain the parent location on the child grid.
!
                DO j=Ymin,Ymax
                  DO i=Xmin,Xmax
                    ic=NINT(OCEAN(ng)%hindi_r(i,j))
                    jc=NINT(OCEAN(ng)%hindj_r(i,j))
                    IF ((ic.ne.IJspv).and.(jc.ne.IJspv)) THEN
                      my_val=A(ic,jc)
# ifdef MASKING
                      my_val=my_val*GRID(ng)%rmask(i,j)
# endif
# ifdef WET_DRY
                      IF (my_val.le.(Dcrit(ng)-GRID(ng)%h(i,j))) THEN
                        my_val=Dcrit(ng)-GRID(ng)%h(i,j)
                      END IF
# endif
                      OCEAN(ng)%zeta(i,j,idtim)=my_val
                    END IF
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
            END DO
# ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints, EWperiodic, NSperiodic,   &
     &                          OCEAN(ng)%zeta(:,:,idtim))
# endif
          END DO
        END DO
      ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!         inflag = 14 : rzeta
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF (inflag.eq.14)  THEN
        ntim=1
        DO m=1,ntim
          DO ng=1,Ngrids
            idtim=krhs(ng)
            LBi=BOUNDS(ng)%LBi(tile)
            UBi=BOUNDS(ng)%UBi(tile)
            LBj=BOUNDS(ng)%LBj(tile)
            UBj=BOUNDS(ng)%UBj(tile)
            DO ngg=1,Ngrids-1
              ngc=COMPOSED(ng)%mychild(ngg)
!
              LBic=BOUNDS(ngc)%LBi(tile)
              UBic=BOUNDS(ngc)%UBi(tile)
              LBjc=BOUNDS(ngc)%LBj(tile)
              UBjc=BOUNDS(ngc)%UBj(tile)

              gtype=r2dvar
              Imin=LOWER_BOUND_I
              Imax=Lm(ngc)+GHOST_POINTS
              Jmin=LOWER_BOUND_J
              Jmax=Mm(ngc)+GHOST_POINTS
              IF ((ABS(gtype).eq.p2dvar).or.(ABS(gtype).eq.u2dvar)) THEN
                Imin=LOWER_BOUND_I+1
              END IF
              IF ((ABS(gtype).eq.p2dvar).or.(ABS(gtype).eq.v2dvar)) THEN
                Jmin=LOWER_BOUND_J+1
              END IF
              Ilen=Imax-Imin+1
              Jlen=Jmax-Jmin+1
              IJlen=Ilen*Jlen
              Npts=IJlen

              IF (.not.allocated(Aout)) THEN
                allocate (Aout(IJlen))
              END IF
              IF (.not.allocated(A)) THEN
                allocate(A(Imin:Imax,Jmin:Jmax))
              END IF
!
! get zeta
!
              CALL mp_gather2d (ngc, model, LBic, UBic, LBjc, UBjc, 0,  &
     &                          gtype, 1.0_r8,                          &
# ifdef MASKING
     &               GRID(ngc) % rmask(LBic:UBic,LBjc:UBjc),            &
# endif
     &               OCEAN(ngc) % rzeta(LBic:UBic,LBjc:UBjc,idtim),     &
     &               Npts, Aout)
!
! distribute the Aout (1D packed) field to all the nodes for ng=2
!
              CALL mpi_bcast (Aout, IJlen, MP_FLOAT, MyMaster,          &
     &                        OCN_COMM_WORLD, MyError)
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
              Xmin=OCEAN(ng)%hindr(1,ngg)
              Xmax=OCEAN(ng)%hindr(2,ngg)
              Ymin=OCEAN(ng)%hindr(3,ngg)
              Ymax=OCEAN(ng)%hindr(4,ngg)
              IF (Xmax.gt.-9999) THEN
!
! Set the chld data into the parent locations. A is on the child grid.
! Use hindices to obtain the parent location on the child grid.
!
                DO j=Ymin,Ymax
                  DO i=Xmin,Xmax
                    ic=NINT(OCEAN(ng)%hindi_r(i,j))
                    jc=NINT(OCEAN(ng)%hindj_r(i,j))
                    IF ((ic.ne.IJspv).and.(jc.ne.IJspv)) THEN
                      my_val=A(ic,jc)
# ifdef MASKING
                      my_val=my_val*GRID(ng)%rmask(i,j)
# endif
                      OCEAN(ng)%rzeta(i,j,idtim)=my_val
                    END IF
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
            END DO
# ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints, EWperiodic, NSperiodic,   &
     &                          OCEAN(ng)%rzeta(:,:,idtim))
# endif
          END DO
        END DO
      ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!         inflag = 21 : Zt_avg1
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF (inflag.eq.21)  THEN
!       ntim=1
!       DO m=1,ntim
          DO ng=1,Ngrids
            idtim=krhs(ng)
            LBi=BOUNDS(ng)%LBi(tile)
            UBi=BOUNDS(ng)%UBi(tile)
            LBj=BOUNDS(ng)%LBj(tile)
            UBj=BOUNDS(ng)%UBj(tile)
            DO ngg=1,Ngrids-1
              ngc=COMPOSED(ng)%mychild(ngg)
!
              LBic=BOUNDS(ngc)%LBi(tile)
              UBic=BOUNDS(ngc)%UBi(tile)
              LBjc=BOUNDS(ngc)%LBj(tile)
              UBjc=BOUNDS(ngc)%UBj(tile)

              gtype=r2dvar
              Imin=LOWER_BOUND_I
              Imax=Lm(ngc)+GHOST_POINTS
              Jmin=LOWER_BOUND_J
              Jmax=Mm(ngc)+GHOST_POINTS
              IF ((ABS(gtype).eq.p2dvar).or.(ABS(gtype).eq.u2dvar)) THEN
                Imin=LOWER_BOUND_I+1
              END IF
              IF ((ABS(gtype).eq.p2dvar).or.(ABS(gtype).eq.v2dvar)) THEN
                Jmin=LOWER_BOUND_J+1
              END IF
              Ilen=Imax-Imin+1
              Jlen=Jmax-Jmin+1
              IJlen=Ilen*Jlen
              Npts=IJlen

              IF (.not.allocated(Aout)) THEN
                allocate (Aout(IJlen))
              END IF
              IF (.not.allocated(A)) THEN
                allocate(A(Imin:Imax,Jmin:Jmax))
              END IF
!
! get zeta
!
              CALL mp_gather2d (ngc, model, LBic, UBic, LBjc, UBjc, 0,  &
     &                          gtype, 1.0_r8,                          &
# ifdef MASKING
     &               GRID(ngc) % rmask(LBic:UBic,LBjc:UBjc),            &
# endif
     &               COUPLING(ngc) % Zt_avg1(LBic:UBic,LBjc:UBjc),      &
     &               Npts, Aout)
!
! distribute the Aout (1D packed) field to all the nodes for ng=2
!
              CALL mpi_bcast (Aout, IJlen, MP_FLOAT, MyMaster,          &
     &                        OCN_COMM_WORLD, MyError)
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
              Xmin=OCEAN(ng)%hindr(1,ngg)
              Xmax=OCEAN(ng)%hindr(2,ngg)
              Ymin=OCEAN(ng)%hindr(3,ngg)
              Ymax=OCEAN(ng)%hindr(4,ngg)
              IF (Xmax.gt.-9999) THEN
!
! Set the chld data into the parent locations. A is on the child grid.
! Use hindices to obtain the parent location on the child grid.
!
                DO j=Ymin,Ymax
                  DO i=Xmin,Xmax
                    ic=NINT(OCEAN(ng)%hindi_r(i,j))
                    jc=NINT(OCEAN(ng)%hindj_r(i,j))
                    IF ((ic.ne.IJspv).and.(jc.ne.IJspv)) THEN
                      my_val=A(ic,jc)
# ifdef MASKING
                      my_val=my_val*GRID(ng)%rmask(i,j)
# endif
                      COUPLING(ng)%Zt_avg1(i,j)=my_val
                    END IF
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
            END DO
# ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints, EWperiodic, NSperiodic,   &
     &                          COUPLING(ng)%Zt_avg1(:,:))
# endif
          END DO
!       END DO
      ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!         inflag = 14, 15, or 16 : ubar
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF ((inflag.eq.14).or.(inflag.eq.15).or.(inflag.eq.16)) THEN
!
        IF (inflag.eq.16) THEN
          ntim=2
        ELSE
          ntim=1
        END IF
        DO m=1,ntim
          DO ng=1,Ngrids
            IF (inflag.eq.16) THEN
              idtim=m
            ELSE
              idtim=knew(ng)
            END IF
            LBi=BOUNDS(ng)%LBi(tile)
            UBi=BOUNDS(ng)%UBi(tile)
            LBj=BOUNDS(ng)%LBj(tile)
            UBj=BOUNDS(ng)%UBj(tile)
            DO ngg=1,Ngrids-1
              ngc=COMPOSED(ng)%mychild(ngg)
!
              LBic=BOUNDS(ngc)%LBi(tile)
              UBic=BOUNDS(ngc)%UBi(tile)
              LBjc=BOUNDS(ngc)%LBj(tile)
              UBjc=BOUNDS(ngc)%UBj(tile)
!
              gtype=u2dvar
              Imin=LOWER_BOUND_I
              Imax=Lm(ngc)+GHOST_POINTS
              Jmin=LOWER_BOUND_J
              Jmax=Mm(ngc)+GHOST_POINTS
              IF ((ABS(gtype).eq.p2dvar).or.(ABS(gtype).eq.u2dvar)) THEN
                Imin=LOWER_BOUND_I+1
              END IF
              IF ((ABS(gtype).eq.p2dvar).or.(ABS(gtype).eq.v2dvar)) THEN
                Jmin=LOWER_BOUND_J+1
              END IF
              Ilen=Imax-Imin+1
              Jlen=Jmax-Jmin+1
              IJlen=Ilen*Jlen
              Npts=IJlen

              IF (.not.allocated(Aout)) THEN
                allocate (Aout(IJlen))
              END IF
              IF (.not.allocated(A)) THEN
                allocate(A(Imin:Imax,Jmin:Jmax))
              END IF
!
! get ubar
!
              CALL mp_gather2d (ngc, model, LBic, UBic, LBjc, UBjc, 0,    &
     &                          gtype, 1.0_r8,                            &
# ifdef MASKING
     &               GRID(ngc) % umask(LBic:UBic,LBjc:UBjc),            &
# endif
     &               OCEAN(ngc) % ubar(LBic:UBic,LBjc:UBjc,idtim),      &
     &               Npts, Aout)
!
! distribute the Aout (1D packed) field to all the nodes
!
              CALL mpi_bcast (Aout, IJlen, MP_FLOAT, MyMaster,            &
     &                        OCN_COMM_WORLD, MyError)
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
              Xmin=OCEAN(ng)%hindu(1,ngg)
              Xmax=OCEAN(ng)%hindu(2,ngg)
              Ymin=OCEAN(ng)%hindu(3,ngg)
              Ymax=OCEAN(ng)%hindu(4,ngg)
              IF (Xmax.gt.-9999) THEN
!
! Set the chld data into the parent locations. A is on the child grid.
! Use hindices to obtain the parent location on the child grid.
!
                DO j=Ymin,Ymax
                  DO i=Xmin,Xmax
                    ic=NINT(OCEAN(ng)%hindi_u(i,j))
                    jc=NINT(OCEAN(ng)%hindj_u(i,j))
                    IF ((ic.ne.IJspv).and.(jc.ne.IJspv)) THEN
                      my_val=A(ic,jc)
# ifdef MASKING
                      my_val=my_val*GRID(ng)%umask(i,j)
# endif
                      OCEAN(ng)%ubar(i,j,idtim)=my_val
                    END IF
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
            END DO
# ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                            &
     &                          LBi, UBi, LBj, UBj,                           &
     &                          NghostPoints, EWperiodic, NSperiodic,         &
     &                          OCEAN(ng)%ubar(:,:,idtim))
# endif
          END DO
        END DO
      ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!         inflag = 21 : DU_avg1
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF (inflag.eq.21) THEN
!        ntim=1
!        DO m=1,ntim
          DO ng=1,Ngrids
            idtim=m
            LBi=BOUNDS(ng)%LBi(tile)
            UBi=BOUNDS(ng)%UBi(tile)
            LBj=BOUNDS(ng)%LBj(tile)
            UBj=BOUNDS(ng)%UBj(tile)
            DO ngg=1,Ngrids-1
              ngc=COMPOSED(ng)%mychild(ngg)
!
              LBic=BOUNDS(ngc)%LBi(tile)
              UBic=BOUNDS(ngc)%UBi(tile)
              LBjc=BOUNDS(ngc)%LBj(tile)
              UBjc=BOUNDS(ngc)%UBj(tile)
!
              gtype=u2dvar
              Imin=LOWER_BOUND_I
              Imax=Lm(ngc)+GHOST_POINTS
              Jmin=LOWER_BOUND_J
              Jmax=Mm(ngc)+GHOST_POINTS
              IF ((ABS(gtype).eq.p2dvar).or.(ABS(gtype).eq.u2dvar)) THEN
                Imin=LOWER_BOUND_I+1
              END IF
              IF ((ABS(gtype).eq.p2dvar).or.(ABS(gtype).eq.v2dvar)) THEN
                Jmin=LOWER_BOUND_J+1
              END IF
              Ilen=Imax-Imin+1
              Jlen=Jmax-Jmin+1
              IJlen=Ilen*Jlen
              Npts=IJlen

              IF (.not.allocated(Aout)) THEN
                allocate (Aout(IJlen))
              END IF
              IF (.not.allocated(A)) THEN
                allocate(A(Imin:Imax,Jmin:Jmax))
              END IF
!
! get ubar
!
              CALL mp_gather2d (ngc, model, LBic, UBic, LBjc, UBjc, 0,  &
     &                          gtype, 1.0_r8,                          &
# ifdef MASKING
     &               GRID(ngc) % umask(LBic:UBic,LBjc:UBjc),            &
# endif
     &               COUPLING(ngc) % DU_avg1(LBic:UBic,LBjc:UBjc),      &
     &               Npts, Aout)
!
! distribute the Aout (1D packed) field to all the nodes
!
              CALL mpi_bcast (Aout, IJlen, MP_FLOAT, MyMaster,          &
     &                        OCN_COMM_WORLD, MyError)
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
              Xmin=OCEAN(ng)%hindu(1,ngg)
              Xmax=OCEAN(ng)%hindu(2,ngg)
              Ymin=OCEAN(ng)%hindu(3,ngg)
              Ymax=OCEAN(ng)%hindu(4,ngg)
              IF (Xmax.gt.-9999) THEN
!
! Set the chld data into the parent locations. A is on the child grid.
! Use hindices to obtain the parent location on the child grid.
!
                DO j=Ymin,Ymax
                  DO i=Xmin,Xmax
                    ic=NINT(OCEAN(ng)%hindi_u(i,j))
                    jc=NINT(OCEAN(ng)%hindj_u(i,j))
                    IF ((ic.ne.IJspv).and.(jc.ne.IJspv)) THEN
                      my_val=A(ic,jc)
# ifdef MASKING
                      my_val=my_val*GRID(ng)%umask(i,j)
# endif
                      COUPLING(ng)%DU_avg1(i,j)=my_val
                    END IF
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
            END DO
# ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                            &
     &                          LBi, UBi, LBj, UBj,                           &
     &                          NghostPoints, EWperiodic, NSperiodic,         &
     &                          COUPLING(ng)%DU_avg1(:,:))
# endif
          END DO
!       END DO
      ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!         inflag = 14, 15, or 16 : vbar
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF ((inflag.eq.14).or.(inflag.eq.15).or.(inflag.eq.16)) THEN
!
        IF (inflag.eq.16) THEN
          ntim=2
        ELSE
          ntim=1
        END IF
        DO m=1,ntim
          DO ng=1,Ngrids
            IF (inflag.eq.16) THEN
              idtim=m
            ELSE
              idtim=knew(ng)
            END IF
            LBi=BOUNDS(ng)%LBi(tile)
            UBi=BOUNDS(ng)%UBi(tile)
            LBj=BOUNDS(ng)%LBj(tile)
            UBj=BOUNDS(ng)%UBj(tile)
            DO ngg=1,Ngrids-1
              ngc=COMPOSED(ng)%mychild(ngg)
!
              LBic=BOUNDS(ngc)%LBi(tile)
              UBic=BOUNDS(ngc)%UBi(tile)
              LBjc=BOUNDS(ngc)%LBj(tile)
              UBjc=BOUNDS(ngc)%UBj(tile)
!
              gtype=v2dvar
              Imin=LOWER_BOUND_I
              Imax=Lm(ngc)+GHOST_POINTS
              Jmin=LOWER_BOUND_J
              Jmax=Mm(ngc)+GHOST_POINTS
              IF ((ABS(gtype).eq.p2dvar).or.(ABS(gtype).eq.u2dvar)) THEN
                Imin=LOWER_BOUND_I+1
              END IF
              IF ((ABS(gtype).eq.p2dvar).or.(ABS(gtype).eq.v2dvar)) THEN
                Jmin=LOWER_BOUND_J+1
              END IF
              Ilen=Imax-Imin+1
              Jlen=Jmax-Jmin+1
              IJlen=Ilen*Jlen
              Npts=IJlen

              IF (.not.allocated(Aout)) THEN
                allocate (Aout(IJlen))
              END IF
              IF (.not.allocated(A)) THEN
                allocate(A(Imin:Imax,Jmin:Jmax))
              END IF
!
! get ubar
!
              CALL mp_gather2d (ngc, model, LBic, UBic, LBjc, UBjc, 0,    &
     &                          gtype, 1.0_r8,                            &
# ifdef MASKING
     &               GRID(ngc) % vmask(LBic:UBic,LBjc:UBjc),            &
# endif
     &               OCEAN(ngc) % vbar(LBic:UBic,LBjc:UBjc,idtim),      &
     &               Npts, Aout)
!
! distribute the Aout (1D packed) field to all the nodes
!
              CALL mpi_bcast (Aout, IJlen, MP_FLOAT, MyMaster,            &
     &                        OCN_COMM_WORLD, MyError)
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
              Xmin=OCEAN(ng)%hindv(1,ngg)
              Xmax=OCEAN(ng)%hindv(2,ngg)
              Ymin=OCEAN(ng)%hindv(3,ngg)
              Ymax=OCEAN(ng)%hindv(4,ngg)
              IF (Xmax.gt.-9999) THEN
!
! Set the chld data into the parent locations. A is on the child grid.
! Use hindices to obtain the parent location on the child grid.
!
                DO j=Ymin,Ymax
                  DO i=Xmin,Xmax
                    ic=NINT(OCEAN(ng)%hindi_v(i,j))
                    jc=NINT(OCEAN(ng)%hindj_v(i,j))
                    IF ((ic.ne.IJspv).and.(jc.ne.IJspv)) THEN
                      my_val=A(ic,jc)
# ifdef MASKING
                      my_val=my_val*GRID(ng)%vmask(i,j)
# endif
                      OCEAN(ng)%vbar(i,j,idtim)=my_val
                    END IF
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
            END DO
# ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                            &
     &                          LBi, UBi, LBj, UBj,                           &
     &                          NghostPoints, EWperiodic, NSperiodic,         &
     &                          OCEAN(ng)%vbar(:,:,idtim))
# endif
          END DO
        END DO
      ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!         inflag = 21 : DV_avg1
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF (inflag.eq.21) THEN
!       ntim=1
!       DO m=1,ntim
          DO ng=1,Ngrids
            idtim=m
            LBi=BOUNDS(ng)%LBi(tile)
            UBi=BOUNDS(ng)%UBi(tile)
            LBj=BOUNDS(ng)%LBj(tile)
            UBj=BOUNDS(ng)%UBj(tile)
            DO ngg=1,Ngrids-1
              ngc=COMPOSED(ng)%mychild(ngg)
!
              LBic=BOUNDS(ngc)%LBi(tile)
              UBic=BOUNDS(ngc)%UBi(tile)
              LBjc=BOUNDS(ngc)%LBj(tile)
              UBjc=BOUNDS(ngc)%UBj(tile)
!
              gtype=v2dvar
              Imin=LOWER_BOUND_I
              Imax=Lm(ngc)+GHOST_POINTS
              Jmin=LOWER_BOUND_J
              Jmax=Mm(ngc)+GHOST_POINTS
              IF ((ABS(gtype).eq.p2dvar).or.(ABS(gtype).eq.u2dvar)) THEN
                Imin=LOWER_BOUND_I+1
              END IF
              IF ((ABS(gtype).eq.p2dvar).or.(ABS(gtype).eq.v2dvar)) THEN
                Jmin=LOWER_BOUND_J+1
              END IF
              Ilen=Imax-Imin+1
              Jlen=Jmax-Jmin+1
              IJlen=Ilen*Jlen
              Npts=IJlen

              IF (.not.allocated(Aout)) THEN
                allocate (Aout(IJlen))
              END IF
              IF (.not.allocated(A)) THEN
                allocate(A(Imin:Imax,Jmin:Jmax))
              END IF
!
! get ubar
!
              CALL mp_gather2d (ngc, model, LBic, UBic, LBjc, UBjc, 0,  &
     &                          gtype, 1.0_r8,                          &
# ifdef MASKING
     &               GRID(ngc) % vmask(LBic:UBic,LBjc:UBjc),            &
# endif
     &               COUPLING(ngc) % DV_avg1(LBic:UBic,LBjc:UBjc),      &
     &               Npts, Aout)
!
! distribute the Aout (1D packed) field to all the nodes
!
              CALL mpi_bcast (Aout, IJlen, MP_FLOAT, MyMaster,          &
     &                        OCN_COMM_WORLD, MyError)
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
              Xmin=OCEAN(ng)%hindv(1,ngg)
              Xmax=OCEAN(ng)%hindv(2,ngg)
              Ymin=OCEAN(ng)%hindv(3,ngg)
              Ymax=OCEAN(ng)%hindv(4,ngg)
              IF (Xmax.gt.-9999) THEN
!
! Set the chld data into the parent locations. A is on the child grid.
! Use hindices to obtain the parent location on the child grid.
!
                DO j=Ymin,Ymax
                  DO i=Xmin,Xmax
                    ic=NINT(OCEAN(ng)%hindi_v(i,j))
                    jc=NINT(OCEAN(ng)%hindj_v(i,j))
                    IF ((ic.ne.IJspv).and.(jc.ne.IJspv)) THEN
                      my_val=A(ic,jc)
# ifdef MASKING
                      my_val=my_val*GRID(ng)%vmask(i,j)
# endif
                      COUPLING(ng)%DV_avg1(i,j)=my_val
                    END IF
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
            END DO
# ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                            &
     &                          LBi, UBi, LBj, UBj,                           &
     &                          NghostPoints, EWperiodic, NSperiodic,         &
     &                          COUPLING(ng)%DV_avg1(:,:))
# endif
          END DO
!       END DO
      ENDIF

      RETURN
      END SUBROUTINE composed_2d

#endif
      END MODULE composed_2d_mod
