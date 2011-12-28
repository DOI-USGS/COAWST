#include "cppdefs.h"
      MODULE composed_3d_mod
#if defined COMPOSED_GRID && defined SOLVE3D
!=======================================================================
!                                                                      !
!  This routine allocates and exchanges data for composed grids.       !
!                                                                      !
!=======================================================================

      implicit none

      PRIVATE
      PUBLIC  :: composed_3d

      CONTAINS
!
!***********************************************************************
      SUBROUTINE composed_3d (tile, inflag)
!***********************************************************************
!
      USE mod_param
      USE mod_coupling
      USE mod_ncparam
      USE mod_scalars
      USE mod_stepping
      USE mod_composed
      USE mod_grid
      USE mod_ocean
      USE mod_parallel
# ifdef DISTRIBUTE
      USE distribute_mod, ONLY : mp_bcasti, mp_gather3d
      USE mp_exchange_mod, ONLY : mp_exchange3d, mp_exchange4d
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
      integer :: i, j, ic, jc, k, ngg, ngc
      integer :: LBic, UBic, LBjc, UBjc
      integer :: Ilen, Jlen, IJlen, Klen
      integer :: Npts, MyError, MySize
      integer :: Xmin, Xmax, Ymin, Ymax
      integer :: Imin, Imax, Jmin, Jmax, Kmin, Kmax, model

      real :: my_val
      real(r8), allocatable :: Aout(:)
      real(r8), allocatable :: A(:,:,:)

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!         inflag = 12 or 18 : tracer
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# if !defined TS_FIXED
      IF ((inflag.eq.12).or.(inflag.eq.18)) THEN
        DO ng=1,Ngrids
          DO itrc=1,NT(ng)
            IF (inflag.eq.12) THEN
              idtim=3
            ELSE
              idtim=nnew(ng)
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

              gtype=r3dvar
              Imin=LOWER_BOUND_I
              Imax=Lm(ngc)+GHOST_POINTS
              Jmin=LOWER_BOUND_J
              Jmax=Mm(ngc)+GHOST_POINTS
              IF ((ABS(gtype).eq.p3dvar).or.(ABS(gtype).eq.u3dvar)) THEN
                Imin=LOWER_BOUND_I+1
              END IF
              IF ((ABS(gtype).eq.p3dvar).or.(ABS(gtype).eq.v3dvar)) THEN
                Jmin=LOWER_BOUND_J+1
              END IF
              Kmin=1
              Kmax=N(ngc)
              Ilen=Imax-Imin+1
              Jlen=Jmax-Jmin+1
              Klen=Kmax-Kmin+1
              IJlen=Ilen*Jlen
              Npts=IJlen*Klen
              MySize=Npts

              IF (.not.allocated(Aout)) THEN
                allocate (Aout(Npts))
              END IF
              IF (.not.allocated(A)) THEN
                allocate(A(Imin:Imax,Jmin:Jmax,Kmin:Kmax))
              END IF
!
! get tracer
!
              CALL mp_gather3d (ngc, model, LBic, UBic, LBjc, UBjc,     &
     &                          1, N(ngc), 0, gtype, 1.0_r8,            &
# ifdef MASKING
     &                        GRID(ngc) % rmask(LBic:UBic,LBjc:UBjc),   &
# endif
     &   OCEAN(ngc) % t(LBic:UBic,LBjc:UBjc,1:N(ngc),idtim,itrc),       &
     &                        Npts, Aout)
!
! distribute the Aout (1D packed) field to all the nodes for ng=2
!
              CALL mpi_bcast (Aout, MySize, MP_FLOAT, MyMaster,         &
     &                        OCN_COMM_WORLD, MyError)
!
! unpack and place in global 3D array
!
              ic=0
              DO k=Kmin,Kmax
                DO j=Jmin,Jmax
                  DO i=Imin,Imax
                    ic=ic+1
                    A(i,j,k)=Aout(ic)
                  END DO
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
                DO k=Kmin,Kmax
                  DO j=Ymin,Ymax
                    DO i=Xmin,Xmax
                      ic=NINT(OCEAN(ng)%hindi_r(i,j))
                      jc=NINT(OCEAN(ng)%hindj_r(i,j))
                      my_val=A(ic,jc,k)
# ifdef MASKING
                      my_val=my_val*GRID(ng)%rmask(i,j)
# endif
                      OCEAN(ng)%t(i,j,k,nstp(ng),itrc)=my_val
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
            END DO
# ifdef DISTRIBUTE
            CALL mp_exchange4d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj, 1, N(ng), 1, NT(ng),&
     &                          NghostPoints, EWperiodic, NSperiodic,   &
     &                          OCEAN(ng)%t(:,:,:,idtim,:))
# endif
          END DO
        END DO
      ENDIF
# endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!         inflag = 16 : u
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF (inflag.eq.16) THEN
        DO ng=1,Ngrids
          idtim=nnew(ng)
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
            gtype=u3dvar
            Imin=LOWER_BOUND_I
            Imax=Lm(ngc)+GHOST_POINTS
            Jmin=LOWER_BOUND_J
            Jmax=Mm(ngc)+GHOST_POINTS
            IF ((ABS(gtype).eq.p3dvar).or.(ABS(gtype).eq.u3dvar)) THEN
              Imin=LOWER_BOUND_I+1
            END IF
            IF ((ABS(gtype).eq.p3dvar).or.(ABS(gtype).eq.v3dvar)) THEN
              Jmin=LOWER_BOUND_J+1
            END IF
            Kmin=1
            Kmax=N(ngc)
            Ilen=Imax-Imin+1
            Jlen=Jmax-Jmin+1
            Klen=Kmax-Kmin+1
            IJlen=Ilen*Jlen
            Npts=IJlen*Klen
            MySize=Npts

            IF (.not.allocated(Aout)) THEN
              allocate (Aout(Npts))
            END IF
            IF (.not.allocated(A)) THEN
              allocate(A(Imin:Imax,Jmin:Jmax,Kmin:Kmax))
            END IF
!
! get u
!
            CALL mp_gather3d (ngc, model, LBic, UBic, LBjc, UBjc,       &
     &                        1, N(ngc), 0, gtype, 1.0_r8,              &
# ifdef MASKING
     &                        GRID(ngc) % umask(LBic:UBic,LBjc:UBjc),   &
# endif
     &            OCEAN(ngc) % u(LBic:UBic,LBjc:UBjc,1:N(ngc),idtim),   &
     &                        Npts, Aout)
!
! distribute the Aout (1D packed) field to all the nodes
!
            CALL mpi_bcast (Aout, MySize, MP_FLOAT, MyMaster,            &
     &                      OCN_COMM_WORLD, MyError)
!
! unpack and place in global 3D array
!
            ic=0
            DO k=Kmin,Kmax
              DO j=Jmin,Jmax
                DO i=Imin,Imax
                  ic=ic+1
                  A(i,j,k)=Aout(ic)
                END DO
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
              DO k=Kmin,Kmax
                DO j=Ymin,Ymax
                  DO i=Xmin,Xmax
                    ic=NINT(OCEAN(ng)%hindi_u(i,j))
                    jc=NINT(OCEAN(ng)%hindj_u(i,j))
                    my_val=A(ic,jc,k)
# ifdef MASKING
                    my_val=my_val*GRID(ng)%umask(i,j)
# endif
                    OCEAN(ng)%u(i,j,k,idtim)=my_val
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
          END DO
# ifdef DISTRIBUTE
          CALL mp_exchange3d (ng, tile, iNLM, 1,                        &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        NghostPoints, EWperiodic, NSperiodic,     &
     &                        OCEAN(ng)%u(:,:,:,idtim))
# endif
        END DO
      ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!         inflag = 16 : Huon
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF (inflag.eq.16) THEN
        DO ng=1,Ngrids
          idtim=nnew(ng)
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
            gtype=u3dvar
            Imin=LOWER_BOUND_I
            Imax=Lm(ngc)+GHOST_POINTS
            Jmin=LOWER_BOUND_J
            Jmax=Mm(ngc)+GHOST_POINTS
            IF ((ABS(gtype).eq.p3dvar).or.(ABS(gtype).eq.u3dvar)) THEN
              Imin=LOWER_BOUND_I+1
            END IF
            IF ((ABS(gtype).eq.p3dvar).or.(ABS(gtype).eq.v3dvar)) THEN
              Jmin=LOWER_BOUND_J+1
            END IF
            Kmin=1
            Kmax=N(ngc)
            Ilen=Imax-Imin+1
            Jlen=Jmax-Jmin+1
            Klen=Kmax-Kmin+1
            IJlen=Ilen*Jlen
            Npts=IJlen*Klen
            MySize=Npts

            IF (.not.allocated(Aout)) THEN
              allocate (Aout(Npts))
            END IF
            IF (.not.allocated(A)) THEN
              allocate(A(Imin:Imax,Jmin:Jmax,Kmin:Kmax))
            END IF
!
! get Huon
!
            CALL mp_gather3d (ngc, model, LBic, UBic, LBjc, UBjc,       &
     &                        1, N(ngc), 0, gtype, 1.0_r8,              &
# ifdef MASKING
     &                        GRID(ngc) % umask(LBic:UBic,LBjc:UBjc),   &
# endif
     &                GRID(ngc) % Huon(LBic:UBic,LBjc:UBjc,1:N(ngc)),   &
     &                        Npts, Aout)
!
! distribute the Aout (1D packed) field to all the nodes
!
            CALL mpi_bcast (Aout, MySize, MP_FLOAT, MyMaster,            &
     &                      OCN_COMM_WORLD, MyError)
!
! unpack and place in global 3D array
!
            ic=0
            DO k=Kmin,Kmax
              DO j=Jmin,Jmax
                DO i=Imin,Imax
                  ic=ic+1
                  A(i,j,k)=Aout(ic)
                END DO
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
              DO k=Kmin,Kmax
                DO j=Ymin,Ymax
                  DO i=Xmin,Xmax
                    ic=NINT(OCEAN(ng)%hindi_u(i,j))
                    jc=NINT(OCEAN(ng)%hindj_u(i,j))
                    my_val=A(ic,jc,k)
# ifdef MASKING
                    my_val=my_val*GRID(ng)%umask(i,j)
# endif
                    GRID(ng)%Huon(i,j,k)=my_val
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
          END DO
# ifdef DISTRIBUTE
          CALL mp_exchange3d (ng, tile, iNLM, 1,                        &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        NghostPoints, EWperiodic, NSperiodic,     &
     &                        GRID(ng)%Huon(:,:,:))
# endif
        END DO
      ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!         inflag = 16 : v
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF (inflag.eq.16) THEN
        DO ng=1,Ngrids
          idtim=nnew(ng)
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
            gtype=v3dvar
            Imin=LOWER_BOUND_I
            Imax=Lm(ngc)+GHOST_POINTS
            Jmin=LOWER_BOUND_J
            Jmax=Mm(ngc)+GHOST_POINTS
            IF ((ABS(gtype).eq.p3dvar).or.(ABS(gtype).eq.u3dvar)) THEN
              Imin=LOWER_BOUND_I+1
            END IF
            IF ((ABS(gtype).eq.p3dvar).or.(ABS(gtype).eq.v3dvar)) THEN
              Jmin=LOWER_BOUND_J+1
            END IF
            Kmin=1
            Kmax=N(ngc)
            Ilen=Imax-Imin+1
            Jlen=Jmax-Jmin+1
            Klen=Kmax-Kmin+1
            IJlen=Ilen*Jlen
            Npts=IJlen*Klen
            MySize=Npts

            IF (.not.allocated(Aout)) THEN
              allocate (Aout(Npts))
            END IF
            IF (.not.allocated(A)) THEN
              allocate(A(Imin:Imax,Jmin:Jmax,Kmin:Kmax))
            END IF
!
! get v
!
            CALL mp_gather3d (ngc, model, LBic, UBic, LBjc, UBjc,       &
     &                        1, N(ngc), 0, gtype, 1.0_r8,              &
# ifdef MASKING
     &                        GRID(ngc) % vmask(LBic:UBic,LBjc:UBjc),   &
# endif
     &            OCEAN(ngc) % v(LBic:UBic,LBjc:UBjc,1:N(ngc),idtim),   &
     &                        Npts, Aout)
!
! distribute the Aout (1D packed) field to all the nodes
!
            CALL mpi_bcast (Aout, MySize, MP_FLOAT, MyMaster,            &
     &                      OCN_COMM_WORLD, MyError)
!
! unpack and place in global 3D array
!
            ic=0
            DO k=Kmin,Kmax
              DO j=Jmin,Jmax
                DO i=Imin,Imax
                  ic=ic+1
                  A(i,j,k)=Aout(ic)
                END DO
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
              DO k=Kmin,Kmax
                DO j=Ymin,Ymax
                  DO i=Xmin,Xmax
                    ic=NINT(OCEAN(ng)%hindi_v(i,j))
                    jc=NINT(OCEAN(ng)%hindj_v(i,j))
                    my_val=A(ic,jc,k)
# ifdef MASKING
                    my_val=my_val*GRID(ng)%vmask(i,j)
# endif
                    OCEAN(ng)%v(i,j,k,idtim)=my_val
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
          END DO
# ifdef DISTRIBUTE
          CALL mp_exchange3d (ng, tile, iNLM, 1,                        &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        NghostPoints, EWperiodic, NSperiodic,     &
     &                        OCEAN(ng)%v(:,:,:,idtim))
# endif
        END DO
      ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!         inflag = 16 : Hvom
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF (inflag.eq.16) THEN
        DO ng=1,Ngrids
          idtim=nnew(ng)
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
            gtype=v3dvar
            Imin=LOWER_BOUND_I
            Imax=Lm(ngc)+GHOST_POINTS
            Jmin=LOWER_BOUND_J
            Jmax=Mm(ngc)+GHOST_POINTS
            IF ((ABS(gtype).eq.p3dvar).or.(ABS(gtype).eq.u3dvar)) THEN
              Imin=LOWER_BOUND_I+1
            END IF
            IF ((ABS(gtype).eq.p3dvar).or.(ABS(gtype).eq.v3dvar)) THEN
              Jmin=LOWER_BOUND_J+1
            END IF
            Kmin=1
            Kmax=N(ngc)
            Ilen=Imax-Imin+1
            Jlen=Jmax-Jmin+1
            Klen=Kmax-Kmin+1
            IJlen=Ilen*Jlen
            Npts=IJlen*Klen
            MySize=Npts

            IF (.not.allocated(Aout)) THEN
              allocate (Aout(Npts))
            END IF
            IF (.not.allocated(A)) THEN
              allocate(A(Imin:Imax,Jmin:Jmax,Kmin:Kmax))
            END IF
!
! get Hvom
!
            CALL mp_gather3d (ngc, model, LBic, UBic, LBjc, UBjc,       &
     &                        1, N(ngc), 0, gtype, 1.0_r8,              &
# ifdef MASKING
     &                        GRID(ngc) % vmask(LBic:UBic,LBjc:UBjc),   &
# endif
     &                GRID(ngc) % Hvom(LBic:UBic,LBjc:UBjc,1:N(ngc)),   &
     &                        Npts, Aout)
!
! distribute the Aout (1D packed) field to all the nodes
!
            CALL mpi_bcast (Aout, MySize, MP_FLOAT, MyMaster,            &
     &                      OCN_COMM_WORLD, MyError)
!
! unpack and place in global 3D array
!
            ic=0
            DO k=Kmin,Kmax
              DO j=Jmin,Jmax
                DO i=Imin,Imax
                  ic=ic+1
                  A(i,j,k)=Aout(ic)
                END DO
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
              DO k=Kmin,Kmax
                DO j=Ymin,Ymax
                  DO i=Xmin,Xmax
                    ic=NINT(OCEAN(ng)%hindi_v(i,j))
                    jc=NINT(OCEAN(ng)%hindj_v(i,j))
                    my_val=A(ic,jc,k)
# ifdef MASKING
                    my_val=my_val*GRID(ng)%vmask(i,j)
# endif
                    GRID(ng)%Hvom(i,j,k)=my_val
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
          END DO
# ifdef DISTRIBUTE
          CALL mp_exchange3d (ng, tile, iNLM, 1,                        &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        NghostPoints, EWperiodic, NSperiodic,     &
     &                        GRID(ng)%Hvom(:,:,:))
# endif
        END DO
      ENDIF







# ifdef SOLVE3D
! inflag = 12 rhs3d
      IF (inflag.eq.12) THEN
#  ifdef GLS_MIXING
        id=idMtke            !tke
        gtype=Iinfo(1,id,1)
!        CALL allocate_nesting3d (field3d, id, gtype, 3)
!        CALL interp_nesting3d (field3d)
!
        id=idMtls            !gls
        gtype=Iinfo(1,id,1)
!        CALL allocate_nesting3d (field3d, id, gtype, 3)
!        CALL interp_nesting3d (field3d)
#  endif
      ENDIF
# endif


# ifdef SOLVE3D
! inflag = 17 gls_corstep
      IF (inflag.eq.17) THEN
        id=idVvis            !AKv
        gtype=Iinfo(1,id,1)
!        CALL allocate_nesting3d (field3d, id, gtype)
!        CALL interp_nesting3d (field3d)
!
        id=idTdif            !AKt
        gtype=Iinfo(1,id,1)
!        CALL allocate_nesting4d (field4d, id, gtype)
!        CALL interp_nesting4d (field4d)
!
#  ifdef GLS_MIXING
        id=idMtke            !tke
        gtype=Iinfo(1,id,1)
!        CALL allocate_nesting3d (field3d, id, gtype, nnew(1))
!        CALL interp_nesting3d (field3d)
!
        id=idMtls            !gls
        gtype=Iinfo(1,id,1)
!        CALL allocate_nesting3d (field3d, id, gtype, nnew(1))
!        CALL interp_nesting3d (field3d)
#  endif
      ENDIF
# endif


      RETURN
      END SUBROUTINE composed_3d

#endif
      END MODULE composed_3d_mod
