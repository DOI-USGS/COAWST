#include "cppdefs.h"
      MODULE get_3dparent_data_mod
#ifdef SOLVE3D
# ifdef REFINED_GRID
!
!svn $Id: set_data.F 1405 2008-02-20 01:44:07Z jcwarner $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine acquires date from a parent grid.                      !
!  If there is more than one time-record,                              !
!  data  is loaded  into global two-time record arrays.  The actual    !
!  interpolation is carried elsewhere.                                 !
!                                                                      !
!  Currently,  this routine is only executed in serial mode by the     !
!  main thread.                                                        !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: get_3dparent_data, get_3dparent_data_tile

      CONTAINS
!
!***********************************************************************
      SUBROUTINE get_3dparent_data (ng, tile)
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
#  include "tile.h"

      ngp=myparent(ng)
      LBip=BOUNDS(ngp)%LBi(tile)
      UBip=BOUNDS(ngp)%UBi(tile)
      LBjp=BOUNDS(ngp)%LBj(tile)
      UBjp=BOUNDS(ngp)%UBj(tile)
!
#  ifdef PROFILE
      CALL wclock_on (ng, iNLM, 4)
#  endif
      CALL get_3dparent_data_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             LBip, UBip, LBjp, UBjp)
#  ifdef PROFILE
      CALL wclock_off (ng, iNLM, 4)
#  endif
      RETURN
      END SUBROUTINE get_3dparent_data
!
!***********************************************************************
      SUBROUTINE get_3dparent_data_tile (ng, tile,                      &
     &                                   LBi, UBi, LBj, UBj,            &
     &                                   LBip, UBip, LBjp, UBjp)
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

#  ifdef DISTRIBUTE
      USE distribute_mod, ONLY : mp_bcasti, mp_gather3d
#  endif
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

#  ifdef OBC
      integer :: ILB, IUB, JLB, JUB
#  endif
      integer :: i, ii, gtype, j, k, ic, jc, Npts, ngp, itrc, t1
      integer :: Imin, Imax, Jmin, Jmax, Kmin, Kmax, model
      integer :: Ilen, Jlen, Klen, IJlen, MyType, status
      integer :: MyError, MySize, ncvarid, MyMaskSize
      real(r8) :: Xmin, Xmax, Ymin, Ymax
      real(r8) :: my_min, my_max

      real(r8), allocatable :: Aout(:)
      real(r8), allocatable :: A(:,:,:)
      real(r8), allocatable :: Amask(:,:)
      real(r8), allocatable :: A2(:,:,:)

      logical :: rectangular
      real(r8) :: IJspv
      character (len=20) :: Xname, Yname

!
#  include "set_bounds.h"
!
      ngp=myparent(ng)
      model=iNLM

#  ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn on input data time wall clock.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, iNLM, 3)
#  endif
!
!  Lower and upper bounds of parent grid for tiled arrays.
      IF (iic(1).eq.ntstart(1)) THEN
        t1=2
      ELSE
        t1=nnew(ng)
      END IF

!      write(*,*) 'get3dp ',t1

      DO ii=1,3   ! 3 gtypes
        IF (ii.eq.1) THEN
          gtype=r3dvar
        ELSE IF (ii.eq.2) THEN
          gtype=u3dvar
        ELSE IF (ii.eq.3) THEN
          gtype=v3dvar
        END IF

!       now get the bounds
        IF (ngp.eq.1) THEN
          Imin=0
          Imax=Lm(ngp)+1
          Jmin=0
          Jmax=Mm(ngp)+1
          MyType=gtype
          IF ((ABS(MyType).eq.p3dvar).or.(ABS(MyType).eq.u3dvar)) THEN
            Imin=1
          END IF
          IF ((ABS(MyType).eq.p3dvar).or.(ABS(MyType).eq.v3dvar)) THEN
            Jmin=1
          END IF
        ELSE
          Imin=LOWER_BOUND_I
          Imax=Lm(ngp)+GHOST_POINTS
          Jmin=LOWER_BOUND_J
          Jmax=Mm(ngp)+GHOST_POINTS
          MyType=gtype
          IF ((ABS(MyType).eq.p3dvar).or.(ABS(MyType).eq.u3dvar)) THEN
            Imin=LOWER_BOUND_I+1
          END IF
          IF ((ABS(MyType).eq.p3dvar).or.(ABS(MyType).eq.v3dvar)) THEN
            Jmin=LOWER_BOUND_J+1
          END IF
        END IF
        Kmin=1
        Kmax=N(Ngp)
        Ilen=Imax-Imin+1
        Jlen=Jmax-Jmin+1
        Klen=Kmax-Kmin+1
        IJlen=Ilen*Jlen
        MyMaskSize=IJlen
        MySize=IJlen*Klen
        Npts=MySize

        rectangular=.FALSE.
        IJspv=-9999.0_r8
        ncvarid=1  !not used

        IF (.not.allocated(Aout)) THEN
          allocate (Aout(Npts))
        END IF
        IF (.not.allocated(A)) THEN
          allocate(A(Imin:Imax,Jmin:Jmax,Kmin:Kmax))
        END IF
        IF (.not.allocated(Amask)) THEN
          allocate(Amask(Imin:Imax,Jmin:Jmax))
        END IF
        IF (.not.allocated(A2)) THEN
          allocate(A2(LBi:UBi,LBj:UBj,Kmin:Kmax))
        END IF
!
        IF (ii.eq.1) THEN
!!!!!!!!!!!!!!!!!!!
!  tracer  !!!!!!!!
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
          DO itrc=1,NT(ng)
            Npts=MySize
            CALL mp_gather3d (ngp, model, LBip, UBip, LBjp, UBjp,         &
     &                      1, N(ngp), 0, gtype, 1.0_r8,                   &
#  ifdef MASKING
     &                      GRID(ngp) % rmask(LBip:UBip,LBjp:UBjp),     &
#  endif
     &    OCEAN(ngp) % t(LBip:UBip,LBjp:UBjp,1:N(ngp),NOUTP,itrc),      &
     &                      Npts, Aout)
!
! distribute the Aout (1D packed) field to all the nodes for ng=2
            CALL mpi_bcast (Aout, MySize, MP_FLOAT, MyMaster,           &
     &                      OCN_COMM_WORLD, MyError)
! unpack and place in global ng=1 3D array
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
! Now interpolate to the local tile for ng=2.
!
            IF (spherical) THEN 
              CALL linterp3d_ref (ng, Imin, Imax, Jmin, Jmax,           &
     &                          Kmin, Kmax,                             &
     &                          REFINED(ng)%lonr_p(Imin:Imax,Jmin:Jmax),&
     &                          REFINED(ng)%latr_p(Imin:Imax,Jmin:Jmax),&
     &                          A, Amask,                               &
     &                          LBi, UBi, LBj, UBj, 1, N(ngp),          &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          OCEAN(ng)%hindpi_r,OCEAN(ng)%hindpj_r,  &
     &                          GRID(ng)%lonr, GRID(ng)%latr, A2,       &
     &                          my_min, my_max)
            ELSE
              CALL linterp3d_ref (ng, Imin, Imax, Jmin, Jmax,           &
     &                          Kmin, Kmax,                             &
     &                          REFINED(ng)%xr_p(Imin:Imax,Jmin:Jmax),  &
     &                          REFINED(ng)%yr_p(Imin:Imax,Jmin:Jmax),  &
     &                          A, Amask,                               &
     &                          LBi, UBi, LBj, UBj, 1, N(ngp),          &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          OCEAN(ng)%hindpi_r,OCEAN(ng)%hindpj_r,  &
     &                          GRID(ng)%xr, GRID(ng)%yr, A2,           &
     &                          my_min, my_max)
            END IF
            IF (WESTERN_EDGE) THEN
              DO k=1,N(ng)
                DO j=JstrT,JendT
                  DO i=IstrT,Istr-1
                    ic=i+Nghostpoints+1
                    REFINED(ng)%tR_west(ic,j,k,t1,itrc)=A2(i,j,k)
                  END DO
                END DO
              END DO
            END IF
            IF (EASTERN_EDGE) THEN
              DO k=1,N(ng)
                DO j=JstrT,JendT
                  DO i=Iend+1,IendT
                    ic=i-Iend
                    REFINED(ng)%tR_east(ic,j,k,t1,itrc)=A2(i,j,k)
                  END DO
                END DO
              END DO
            END IF
            IF (SOUTHERN_EDGE) THEN
              DO k=1,N(ng)
                DO j=JstrT,Jstr-1
                  DO i=IstrT,IendT
                    jc=j+Nghostpoints+1
                    REFINED(ng)%tR_south(i,jc,k,t1,itrc)=A2(i,j,k)
                  END DO
                END DO
              END DO
            END IF
            IF (NORTHERN_EDGE) THEN
              DO k=1,N(ng)
                DO j=Jend+1,JendT
                  DO i=IstrT,IendT
                    jc=j-Jend
                    REFINED(ng)%tR_north(i,jc,k,t1,itrc)=A2(i,j,k)
                  END DO
                END DO
              END DO
            END IF
          END DO  !itrc
        ELSE IF (ii.eq.2) THEN
!!!!!!!!!!!!!!!!!!!
!  u     !!!!!!!!!!
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
          Npts=MySize
          CALL mp_gather3d (ngp, model, LBip, UBip, LBjp, UBjp,         &
     &                    1, N(ngp), 0, gtype, 1.0_r8,                  &
#  ifdef MASKING
     &                    GRID(ngp) % umask(LBip:UBip,LBjp:UBjp),       &
#  endif
     &        OCEAN(ngp) % u(LBip:UBip,LBjp:UBjp,1:N(ngp),NOUTP),       &
     &                    Npts, Aout)
!
! distribute the Aout (1D packed) field to all the nodes for ng=2
          CALL mpi_bcast (Aout, MySize, MP_FLOAT, MyMaster,             &
     &                    OCN_COMM_WORLD, MyError)
! unpack and place in global ng=1 3D array
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
! Now interpolate to the local tile for the child grid.
!
          IF (spherical) THEN
            CALL linterp3d_ref (ng, Imin, Imax, Jmin, Jmax, Kmin, Kmax, &
     &                        REFINED(ng)%lonu_p(Imin:Imax,Jmin:Jmax),  &
     &                        REFINED(ng)%latu_p(Imin:Imax,Jmin:Jmax),  &
     &                        A, Amask,                                 &
     &                        LBi, UBi, LBj, UBj, 1, N(ngp),            &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        OCEAN(ng)%hindpi_u, OCEAN(ng)%hindpj_u,   &
     &                        GRID(ng)%lonu, GRID(ng)%latu, A2,         &
     &                        my_min, my_max)
          ELSE
            CALL linterp3d_ref (ng, Imin, Imax, Jmin, Jmax, Kmin, Kmax, &
     &                        REFINED(ng)%xu_p(Imin:Imax,Jmin:Jmax),    &
     &                        REFINED(ng)%yu_p(Imin:Imax,Jmin:Jmax),    &
     &                        A, Amask,                                 &
     &                        LBi, UBi, LBj, UBj, 1, N(ngp),            &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        OCEAN(ng)%hindpi_u, OCEAN(ng)%hindpj_u,   &
     &                        GRID(ng)%xu, GRID(ng)%yu, A2,             &
     &                        my_min, my_max)
          END IF
          IF (WESTERN_EDGE) THEN
            DO k=1,N(ng)
              DO j=JstrT,JendT
                DO i=IstrT+1,Istr
                  ic=i+Nghostpoints
                  REFINED(ng)%uR_west(ic,j,k,t1)=A2(i,j,k)
                END DO
              END DO
            END DO
          END IF
          IF (EASTERN_EDGE) THEN
            DO k=1,N(ng)
              DO j=JstrT,JendT
                DO i=Iend+1,IendT
                  ic=i-Iend
                  REFINED(ng)%uR_east(ic,j,k,t1)=A2(i,j,k)
                END DO
              END DO
            END DO
          END IF
          IF (SOUTHERN_EDGE) THEN
            DO k=1,N(ng)
              DO j=JstrT,Jstr-1
                DO i=IstrT,IendT
                  jc=j+Nghostpoints+1
                  REFINED(ng)%uR_south(i,jc,k,t1)=A2(i,j,k)
                END DO
              END DO
            END DO
          END IF
          IF (NORTHERN_EDGE) THEN
            DO k=1,N(ng)
              DO j=Jend+1,JendT
                DO i=IstrT,IendT
                  jc=j-Jend
                  REFINED(ng)%uR_north(i,jc,k,t1)=A2(i,j,k)
                END DO
              END DO
            END DO
          END IF
        ELSE IF (ii.eq.3) THEN
!!!!!!!!!!!!!!!!!!!
!  v     !!!!!!!!!!
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
          Npts=MySize
          CALL mp_gather3d (ngp, model, LBip, UBip, LBjp, UBjp,         &
     &                    1, N(ngp), 0,  gtype, 1.0_r8,                 &
#  ifdef MASKING
     &                    GRID(ngp) % vmask(LBip:UBip,LBjp:UBjp),       &
#  endif
     &   OCEAN(ngp) % v(LBip:UBip,LBjp:UBjp,1:N(ngp),NOUTP),            &
     &                    Npts, Aout)
!
! distribute the Aout (1D packed) field to all the nodes for ng=2
          CALL mpi_bcast (Aout, MySize, MP_FLOAT, MyMaster,             &
     &                    OCN_COMM_WORLD, MyError)
! unpack and place in global ng=1 3D array
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
! Now interpolate to the local tile for ng=2.
!
          IF (spherical) THEN
            CALL linterp3d_ref (ng, Imin, Imax, Jmin, Jmax, Kmin, Kmax, &
     &                        REFINED(ng)%lonv_p(Imin:Imax,Jmin:Jmax),  &
     &                        REFINED(ng)%latv_p(Imin:Imax,Jmin:Jmax),  &
     &                        A, Amask,                                 &
     &                        LBi, UBi, LBj, UBj, 1, N(ngp),            &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        OCEAN(ng)%hindpi_v, OCEAN(ng)%hindpj_v,   &
     &                        GRID(ng)%lonv, GRID(ng)%latv, A2,         &
     &                        my_min, my_max)
          ELSE
            CALL linterp3d_ref (ng, Imin, Imax, Jmin, Jmax, Kmin, Kmax, &
     &                        REFINED(ng)%xv_p(Imin:Imax,Jmin:Jmax),    &
     &                        REFINED(ng)%yv_p(Imin:Imax,Jmin:Jmax),    &
     &                        A, Amask,                                 &
     &                        LBi, UBi, LBj, UBj, 1, N(ngp),            &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        OCEAN(ng)%hindpi_v, OCEAN(ng)%hindpj_v,   &
     &                        GRID(ng)%xv, GRID(ng)%yv, A2,             &
     &                        my_min, my_max)
          END IF
          IF (WESTERN_EDGE) THEN
            DO k=1,N(ng)
              DO j=JstrT,JendT
                DO i=IstrT,Istr-1
                  ic=i+Nghostpoints+1
                  REFINED(ng)%vR_west(ic,j,k,t1)=A2(i,j,k)
                END DO
              END DO
            END DO
          END IF
          IF (EASTERN_EDGE) THEN
            DO k=1,N(ng)
              DO j=JstrT,JendT
                DO i=Iend+1,IendT
                  ic=i-Iend
                  REFINED(ng)%vR_east(ic,j,k,t1)=A2(i,j,k)
                END DO
              END DO
            END DO
          END IF
          IF (SOUTHERN_EDGE) THEN
            DO k=1,N(ng)
              DO j=JstrT+1,Jstr
                DO i=IstrT,IendT
                  jc=j+Nghostpoints
                  REFINED(ng)%vR_south(i,jc,k,t1)=A2(i,j,k)
                END DO
              END DO
            END DO
          END IF
          IF (NORTHERN_EDGE) THEN
            DO k=1,N(ng)
              DO j=Jend+1,JendT
                DO i=IstrT,IendT
                  jc=j-Jend
                  REFINED(ng)%vR_north(i,jc,k,t1)=A2(i,j,k)
                END DO
              END DO
            END DO
          END IF
        END IF
!
        IF (allocated(Aout)) THEN
          deallocate (Aout)
        END IF
        IF (allocated(Amask)) THEN
          deallocate (Amask)
        END IF
        IF (allocated(A)) THEN
          deallocate (A)
        END IF
        IF (allocated(A2)) THEN
          deallocate (A2)
        END IF
      END DO
!
#  ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn off input data time wall clock.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, iNLM, 3)
#  endif
      RETURN
      END SUBROUTINE get_3dparent_data_tile
# endif
#endif
      END MODULE get_3dparent_data_mod
