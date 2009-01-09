#include "cppdefs.h"
      MODULE nesting_mod
#ifdef COMPOSED_GRID

      implicit none

      PRIVATE
      PUBLIC  :: nesting, interp_nesting2d
      PUBLIC  :: allocate_nesting2d
# ifdef SOLVE3D
      PUBLIC  :: allocate_nesting3d, allocate_nesting4d
      PUBLIC  :: interp_nesting3d, interp_nesting4d
# endif

      CONTAINS

      SUBROUTINE nesting (inflag)
!=======================================================================
!                                                                      !
!  This routine allocates and initializes nesting structure for 2D     !
!  state variables.                                                    !
!                                                                      !
!=======================================================================
!
      USE mod_param
# ifdef SOLVE3D
      USE mod_coupling
# endif
      USE mod_ncparam
      USE mod_scalars
      USE mod_stepping
      USE mod_nesting

      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: inflag
      integer :: id, gtype, Tindex, m, itrc, idtimm
      logical :: dointerp

      TYPE (T_NEST2D), pointer :: field2d(:)
#ifdef SOLVE3D
      TYPE (T_NEST3D), pointer :: field3d(:)
      TYPE (T_NEST4D), pointer :: field4d(:)
#endif
      IF (Ngrids.gt.1) THEN

! inflag = 100 --> we are calling from step2d.
      IF (inflag.eq.100) THEN
        DO m=1,4
          dointerp=.TRUE.
          idtimm=knew(1)
          IF (m.eq.1) id=idFsur  !zeta
          IF (m.eq.2) id=idUBar  !ubar
          IF (m.eq.3) id=idVbar  !vbar
          IF (m.eq.4) THEN
            IF (PREDICTOR_2D_STEP(1)) THEN
              id=idRzet          !rzeta
              idtimm=krhs(1)
            ELSE
              dointerp=.FALSE.
            END IF
          END IF
          gtype=Iinfo(1,id,1)
          IF (dointerp) THEN
            CALL allocate_nesting2d (field2d, id, gtype, idtimm)
            CALL interp_nesting2d (field2d)
          END IF
        END DO
      ENDIF

! inflag = 5 set_massflux, rho_eos
! should these be computed or interpolated !!!!!!!
!     Huon(:,:,:)
!     Hvom(:,:,:)
!     rho
!     rhoA
!     rhoS

! inflag = 9 set_vbc
      IF (inflag.eq.9) THEN
        DO m=1,2
          IF (m.eq.1) THEN
            id=idUbms          !bustr
          END IF
          IF (m.eq.2) THEN
            id=idVbms          !bvstr
          END IF
          gtype=Iinfo(1,id,1)
          CALL allocate_nesting2d (field2d, id, gtype)
          CALL interp_nesting2d (field2d)
        END DO
      END IF

! inflag = 10 ana_vmix, omega, wvelocity
!     # ifdef ana_vmix
!       Akv(:,:,:)
!       AKt(:,:,:,itrc)
!     # endif
      IF (inflag.eq.10) THEN
        DO m=1,2
          IF (m.eq.1) THEN
            id=idUfx1          !DU_avg1
            idtimm=1
          END IF
          IF (m.eq.2) THEN
            id=idVfx1          !DV_avg1
            idtimm=1
            dointerp=.TRUE.
          END IF
          gtype=Iinfo(1,id,1)
!          CALL allocate_nesting2d (field2d, id, gtype, idtimm)
!          CALL interp_nesting2d (field2d)
        END DO
      END IF

! inflag = 11 set_zeta
      IF (inflag.eq.11) THEN
        DO m=1,2
          IF (m.eq.1) THEN
            id=idFsur
            idtimm=1   !zeta
          ELSE IF (m.eq.2) THEN
            id=idFsur
            idtimm=2   !zeta
          END IF
          gtype=Iinfo(1,id,1)
          CALL allocate_nesting2d (field2d, id, gtype, idtimm)
          CALL interp_nesting2d (field2d)
        END DO
      ENDIF

# ifdef SOLVE3D
! inflag = 12 rhs3d
      IF (inflag.eq.12) THEN
        id=idTvar(itemp)     !t
        gtype=Iinfo(1,id,1)
        CALL allocate_nesting4d (field4d, id, gtype, 3)
        CALL interp_nesting4d (field4d)
!       CALL allocate_nesting4d (field4d, id, gtype, nnew(1))
!       CALL interp_nesting4d (field4d)
!
#  ifdef GLS_MIXING
        id=idMtke            !tke
        gtype=Iinfo(1,id,1)
        CALL allocate_nesting3d (field3d, id, gtype, 3)
        CALL interp_nesting3d (field3d)
!
        id=idMtls            !gls
        gtype=Iinfo(1,id,1)
        CALL allocate_nesting3d (field3d, id, gtype, 3)
        CALL interp_nesting3d (field3d)
#  endif
      ENDIF
# endif

! inflag = 14 or 15 step2d
      IF ((inflag.eq.14).or.(inflag.eq.15)) THEN
        DO m=1,4
          dointerp=.FALSE.
          idtimm=knew(1)
          IF (m.eq.1) THEN
            id=idFsur
            dointerp=.TRUE.  !zeta
          ELSE IF (m.eq.2) THEN
            id=idUBar
            dointerp=.TRUE.  !ubar
          ELSE IF (m.eq.3) THEN
            id=idVbar
            dointerp=.TRUE.  !vbar
          ELSE IF ((m.eq.4).and.(PREDICTOR_2D_STEP(1))) THEN
              id=idRzet          !rzeta
              idtimm=krhs(1)
              dointerp=.TRUE.
          END IF
          gtype=Iinfo(1,id,1)
          IF (dointerp) THEN
            CALL allocate_nesting2d (field2d, id, gtype, idtimm)
            CALL interp_nesting2d (field2d)
          END IF
        END DO
      ENDIF

# ifdef SOLVE3D
! inflag = 17 gls_corstep
      IF (inflag.eq.17) THEN
        id=idVvis            !AKv
        gtype=Iinfo(1,id,1)
        CALL allocate_nesting3d (field3d, id, gtype)
        CALL interp_nesting3d (field3d)
!
        id=idTdif            !AKt
        gtype=Iinfo(1,id,1)
        CALL allocate_nesting4d (field4d, id, gtype)
        CALL interp_nesting4d (field4d)
!
#  ifdef GLS_MIXING
        id=idMtke            !tke
        gtype=Iinfo(1,id,1)
        CALL allocate_nesting3d (field3d, id, gtype, nnew(1))
        CALL interp_nesting3d (field3d)
!
        id=idMtls            !gls
        gtype=Iinfo(1,id,1)
        CALL allocate_nesting3d (field3d, id, gtype, nnew(1))
        CALL interp_nesting3d (field3d)
#  endif
      ENDIF
# endif

! inflag = 21
      IF (inflag.eq.21) THEN
        DO m=5,7
          IF (m.eq.5) THEN
            id=idZavg          !Zt_avg1
            idtimm=1
          END IF
          IF (m.eq.6) THEN
            id=idUfx1          !DU_avg1
            idtimm=1
          END IF
          IF (m.eq.7) THEN
            id=idVfx1          !DV_avg1
            idtimm=1
          END IF
          gtype=Iinfo(1,id,1)
          CALL allocate_nesting2d (field2d, id, gtype, idtimm)
          CALL interp_nesting2d (field2d)
        END DO
      ENDIF

# ifdef SOLVE3D
! inflag = 16 step3d_uv
      IF (inflag.eq.16) THEN
        DO m=1,4
          IF (m.eq.1) THEN
            id=idUBar
            idtimm=1   !ubar
          ELSE IF (m.eq.2) THEN
            id=idVbar
            idtimm=1   !vbar
          ELSE IF (m.eq.3) THEN
            id=idUBar
            idtimm=2   !ubar
          ELSE IF (m.eq.4) THEN
            id=idVbar
            idtimm=2   !vbar
          END IF
          gtype=Iinfo(1,id,1)
          CALL allocate_nesting2d (field2d, id, gtype, idtimm)
          CALL interp_nesting2d (field2d)
        END DO
        DO m=1,2
          IF (m.eq.1) THEN
            id=idUvel
            idtimm=nnew(1)   !u
          ELSE IF (m.eq.2) THEN
            id=idVvel
            idtimm=nnew(1)   !v
          END IF
          gtype=Iinfo(1,id,1)
          CALL allocate_nesting3d (field3d, id, gtype, idtimm)
          CALL interp_nesting3d (field3d)
        END DO
        DO m=1,2
          IF (m.eq.1) THEN
            id=idHUav        !Huon
          ELSE IF (m.eq.2) THEN
            id=idHVav        !Hvom
          END IF
          gtype=Iinfo(1,id,1)
          CALL allocate_nesting3d (field3d, id, gtype)
          CALL interp_nesting3d (field3d)
        END DO
      ENDIF

! inflag = 17 omega

! inflag = 18 step3d_t
      IF (inflag.eq.18) THEN
        id=idTvar(itemp)     !t
        gtype=Iinfo(1,id,1)
        CALL allocate_nesting4d (field4d, id, gtype, nnew(1))
        CALL interp_nesting4d (field4d)
      ENDIF
# endif

      ENDIF  !ngrids

      RETURN
      END SUBROUTINE nesting

      SUBROUTINE interp_nesting2d (field2d)
!
!=======================================================================
!                                                                      !
!  This routine interpolates 2d fields for nesting.                    !
!                                                                      !
!=======================================================================
!
      USE mod_ncparam
      USE mod_nesting
      USE interpolate_mod
      USE mod_ocean
!
!  Imported variable declarations.
!
      TYPE (T_NEST2D), pointer :: field2d(:)
!
!  Local variable declarations.
!
      integer :: i, j, ng, ng1, ng2
      integer :: ILB1,IUB1,JLB1,JUB1
      integer :: ILB2,IUB2,JLB2,JUB2
      integer :: Ms, Me, Ns, Ne
      integer :: LBi, UBi, LBj, UBj
      integer :: Xmin, Xmax, Ymin, Ymax

      real(r8) :: IJspv, my_min, my_max

      logical :: EASTERN_EDGE_NESTING
      logical :: WESTERN_EDGE_NESTING
      logical :: rectangular, do_hindices
      real(r8), pointer :: Ipos(:,:)
      real(r8), pointer :: Jpos(:,:)

!     rectangular=.TRUE.
      rectangular=.FALSE.

      DO ng=1,Ngrids
        EASTERN_EDGE_NESTING=.FALSE.
        WESTERN_EDGE_NESTING=.FALSE.
# ifdef DOGBONE
        IF (ng.eq.1) THEN
          EASTERN_EDGE_NESTING=.TRUE.
          ng1=2    !other grid, ie gridded data.
          ng2=1    !local grid, grid that values are interpolated to. 
        END IF
        IF (ng.eq.2) THEN
          WESTERN_EDGE_NESTING=.TRUE.
          ng1=1    !other grid, ie gridded data.
          ng2=2    !local grid, grid that values are interpolated to. 
        END IF
# endif
        ILB1=LBOUND(field2d(ng1)%angler,DIM=1)
        IUB1=UBOUND(field2d(ng1)%angler,DIM=1)
        JLB1=LBOUND(field2d(ng1)%angler,DIM=2)
        JUB1=UBOUND(field2d(ng1)%angler,DIM=2)
        ILB2=LBOUND(field2d(ng2)%angler,DIM=1)
        IUB2=UBOUND(field2d(ng2)%angler,DIM=1)
        JLB2=LBOUND(field2d(ng2)%angler,DIM=2)
        JUB2=UBOUND(field2d(ng2)%angler,DIM=2)
        IJspv=0.0_r8
        IF (WESTERN_EDGE_NESTING) THEN
          Ms=ILB2
          Me=0
          Ns=JLB2
          Ne=JUB2
        END IF
        IF (EASTERN_EDGE_NESTING) THEN
          Ms=IUB2-GHOST_POINTS+1
          Me=IUB2
          Ns=JLB2
          Ne=JUB2
        END IF
        allocate ( Ipos(ILB2:IUB2,JLB2:JUB2) )
        allocate ( Jpos(ILB2:IUB2,JLB2:JUB2) )

        do_hindices=.FALSE.
        IF (field2d(ng2)%gtype.eq.u2dvar) THEN
          IF ( OCEAN(ng2)%hind_init(u2dvar) .eq. .FALSE.) THEN
            OCEAN(ng2)%hind_init(u2dvar)=.TRUE.
            do_hindices=.TRUE.
          ELSE
            IPOS=OCEAN(ng)%hind2di_u
            JPOS=OCEAN(ng)%hind2dj_u
          END IF
        ELSE IF (field2d(ng2)%gtype.eq.v2dvar) THEN
          IF ( OCEAN(ng2)%hind_init(v2dvar) .eq. .FALSE.) THEN
            OCEAN(ng2)%hind_init(v2dvar)=.TRUE.
            do_hindices=.TRUE.
          ELSE
            IPOS=OCEAN(ng)%hind2di_v
            JPOS=OCEAN(ng)%hind2dj_v
          END IF
        ELSE               !r2dvar
          IF ( OCEAN(ng2)%hind_init(r2dvar) .eq. .FALSE.) THEN
            OCEAN(ng2)%hind_init(r2dvar)=.TRUE.
            do_hindices=.TRUE.
          ELSE
            IPOS=OCEAN(ng)%hind2di_r
            JPOS=OCEAN(ng)%hind2dj_r
          END IF
        END IF

        IF (do_hindices) THEN
          CALL hindices (1, ILB1, IUB1, JLB1, JUB1,                     &
     &                   ILB1, IUB1, JLB1, JUB1,                        &
     &                   field2d(ng1)%angler,                           &
     &                   field2d(ng1)%x,   field2d(ng1)%y,              &
     &                   ILB2, IUB2, JLB2, JUB2,                        &
     &                   Ms, Me, Ns, Ne,                                &
     &                   field2d(ng2)%x,   field2d(ng2)%y,              &
     &                   Ipos, Jpos,                                    &
     &                   IJspv, rectangular, Xmin, Xmax, Ymin, Ymax)
          IF (field2d(ng2)%gtype.eq.u2dvar) THEN
            OCEAN(ng)%hind2di_u=IPOS
            OCEAN(ng)%hind2dj_u=JPOS
          ELSE IF (field2d(ng2)%gtype.eq.v2dvar) THEN
            OCEAN(ng)%hind2di_v=IPOS
            OCEAN(ng)%hind2dj_v=JPOS
          ELSE               !r2dvar
            OCEAN(ng)%hind2di_r=IPOS
            OCEAN(ng)%hind2dj_r=JPOS
          END IF
        END IF
        CALL linterp2d (1, ILB1, IUB1, JLB1, JUB1,                      &
     &                   field2d(ng1)%x,   field2d(ng1)%y,              &
     &                   field2d(ng1)%s,                                &
     &                   ILB2, IUB2, JLB2, JUB2,                        &
     &                   Ms, Me, Ns, Ne,                                &
     &                   Ipos,Jpos,                                     &
     &                   field2d(ng2)%x,   field2d(ng2)%y,              &
     &                   field2d(ng2)%s,                                &
     &                   my_min, my_max)
        deallocate ( Ipos, Jpos )
      END DO

      RETURN
      END SUBROUTINE interp_nesting2d

#ifdef SOLVE3D
      SUBROUTINE interp_nesting3d (field3d)
!
!=======================================================================
!                                                                      !
!  This routine interpolates 3d fields for nesting.                    !
!                                                                      !
!=======================================================================
!
      USE mod_ncparam
      USE mod_nesting
      USE interpolate_mod
      USE mod_ocean
!
!  Imported variable declarations.
!
      TYPE (T_NEST3D), pointer :: field3d(:)
!
!  Local variable declarations.
!
      integer :: i, j, ng, ng1, ng2
      integer :: ILB1,IUB1,JLB1,JUB1
      integer :: ILB2,IUB2,JLB2,JUB2
      integer :: KLB1,KUB1,KLB2,KUB2
      integer :: Ms, Me, Ns, Ne
      integer :: LBi, UBi, LBj, UBj
      integer :: Xmin, Xmax, Ymin, Ymax

      real(r8) :: IJspv, my_min, my_max

      logical :: EASTERN_EDGE_NESTING
      logical :: WESTERN_EDGE_NESTING
      logical :: rectangular, do_hindices
      real(r8), pointer :: Ipos(:,:)
      real(r8), pointer :: Jpos(:,:)

!     rectangular=.TRUE.
      rectangular=.FALSE.

      DO ng=1,Ngrids
        EASTERN_EDGE_NESTING=.FALSE.
        WESTERN_EDGE_NESTING=.FALSE.
# ifdef DOGBONE
        IF (ng.eq.1) THEN
          EASTERN_EDGE_NESTING=.TRUE.
          ng1=2    !other grid, ie gridded data.
          ng2=1    !local grid, grid that values are interpolated to. 
        END IF
        IF (ng.eq.2) THEN
          WESTERN_EDGE_NESTING=.TRUE.
          ng1=1    !other grid, ie gridded data.
          ng2=2    !local grid, grid that values are interpolated to. 
        END IF
# endif
        ILB1=LBOUND(field3d(ng1)%angler,DIM=1)
        IUB1=UBOUND(field3d(ng1)%angler,DIM=1)
        JLB1=LBOUND(field3d(ng1)%angler,DIM=2)
        JUB1=UBOUND(field3d(ng1)%angler,DIM=2)
        ILB2=LBOUND(field3d(ng2)%angler,DIM=1)
        IUB2=UBOUND(field3d(ng2)%angler,DIM=1)
        JLB2=LBOUND(field3d(ng2)%angler,DIM=2)
        JUB2=UBOUND(field3d(ng2)%angler,DIM=2)
!       KLB1=field3d(ng1)%LBk
!       KUB1=field3d(ng1)%UBk
!       KLB2=field3d(ng2)%LBk
!       KUB2=field3d(ng2)%UBk
        KUB1=N(ng)
        KUB2=N(ng)
        IF (field3d(ng2)%gtype.eq.w3dvar) THEN
          KLB1=0
          KLB2=0
        ELSE 
          KLB1=1
          KLB2=1
        END IF

        IJspv=0.0_r8
        IF (WESTERN_EDGE_NESTING) THEN
          Ms=ILB2
          Me=0
          Ns=JLB2
          Ne=JUB2
        END IF
        IF (EASTERN_EDGE_NESTING) THEN
          Ms=IUB2-GHOST_POINTS+1
          Me=IUB2
          Ns=JLB2
          Ne=JUB2
        END IF
        allocate ( Ipos(ILB2:IUB2,JLB2:JUB2) )
        allocate ( Jpos(ILB2:IUB2,JLB2:JUB2) )

        do_hindices=.FALSE.
        IF (field3d(ng2)%gtype.eq.u3dvar) THEN
          IF ( OCEAN(ng2)%hind_init(u3dvar) .eq. .FALSE.) THEN
            OCEAN(ng2)%hind_init(u3dvar)=.TRUE.
            do_hindices=.TRUE.
          ELSE
            IPOS=OCEAN(ng)%hind3di_u
            JPOS=OCEAN(ng)%hind3dj_u
          END IF
        ELSE IF (field3d(ng2)%gtype.eq.v3dvar) THEN
          IF ( OCEAN(ng2)%hind_init(v3dvar) .eq. .FALSE.) THEN
            OCEAN(ng2)%hind_init(v3dvar)=.TRUE.
            do_hindices=.TRUE.
          ELSE
            IPOS=OCEAN(ng)%hind3di_v
            JPOS=OCEAN(ng)%hind3dj_v
          END IF
        ELSE               !r3dvar
          IF ( OCEAN(ng2)%hind_init(r3dvar) .eq. .FALSE.) THEN
            OCEAN(ng2)%hind_init(r3dvar)=.TRUE.
            do_hindices=.TRUE.
          ELSE
            IPOS=OCEAN(ng)%hind3di_r
            JPOS=OCEAN(ng)%hind3dj_r
          END IF
        END IF

        IF (do_hindices) THEN
          CALL hindices (ng, ILB1, IUB1, JLB1, JUB1,                    &
     &                   ILB1, IUB1, JLB1, JUB1,                        &
     &                   field3d(ng1)%angler,                           &
     &                   field3d(ng1)%x,   field3d(ng1)%y,              &
     &                   ILB2, IUB2, JLB2, JUB2,                        &
     &                   Ms, Me, Ns, Ne,                                &
     &                   field3d(ng2)%x,   field3d(ng2)%y,              &
     &                   Ipos, Jpos,                                    &
     &                   IJspv, rectangular, Xmin, Xmax, Ymin, Ymax)
          IF (field3d(ng2)%gtype.eq.u3dvar) THEN
            OCEAN(ng)%hind3di_u=IPOS
            OCEAN(ng)%hind3dj_u=JPOS
          ELSE IF (field3d(ng2)%gtype.eq.v3dvar) THEN
            OCEAN(ng)%hind3di_v=IPOS
            OCEAN(ng)%hind3dj_v=JPOS
          ELSE               !r3dvar
            OCEAN(ng)%hind3di_r=IPOS
            OCEAN(ng)%hind3dj_r=JPOS
          END IF
        END IF
        CALL linterp3d (ng, ILB1, IUB1, JLB1, JUB1, KLB1, KUB1,         &
     &                   field3d(ng1)%x,   field3d(ng1)%y,              &
     &                   field3d(ng1)%s,                                &
     &                   ILB2, IUB2, JLB2, JUB2, KLB2, KUB2,            &
     &                   Ms, Me, Ns, Ne,                                &
     &                   Ipos,Jpos,                                     &
     &                   field3d(ng2)%x,   field3d(ng2)%y,              &
     &                   field3d(ng2)%s,                                &
     &                   my_min, my_max)
        deallocate ( Ipos, Jpos )
      END DO

      RETURN
      END SUBROUTINE interp_nesting3d

      SUBROUTINE interp_nesting4d (field4d)
!
!=======================================================================
!                                                                      !
!  This routine interpolates 4d fields for nesting.                    !
!                                                                      !
!=======================================================================
!
      USE mod_ncparam
      USE mod_nesting
      USE interpolate_mod
      USE mod_ocean
!
!  Imported variable declarations.
!
      TYPE (T_NEST4D), pointer :: field4d(:)
!
!  Local variable declarations.
!
      integer :: i, j, ng, ng1, ng2
      integer :: ILB1,IUB1,JLB1,JUB1
      integer :: ILB2,IUB2,JLB2,JUB2
      integer :: KLB1,KUB1,KLB2,KUB2
      integer :: LLB1,LUB1,LLB2,LUB2
      integer :: Ms, Me, Ns, Ne
      integer :: LBi, UBi, LBj, UBj
      integer :: Xmin, Xmax, Ymin, Ymax

      real(r8) :: IJspv, my_min, my_max

      logical :: EASTERN_EDGE_NESTING
      logical :: WESTERN_EDGE_NESTING
      logical :: rectangular, do_hindices
      real(r8), pointer :: Ipos(:,:)
      real(r8), pointer :: Jpos(:,:)

!     rectangular=.TRUE.
      rectangular=.FALSE.

      DO ng=1,Ngrids
        EASTERN_EDGE_NESTING=.FALSE.
        WESTERN_EDGE_NESTING=.FALSE.
# ifdef DOGBONE
        IF (ng.eq.1) THEN
          EASTERN_EDGE_NESTING=.TRUE.
          ng1=2    !other grid, ie gridded data.
          ng2=1    !local grid, grid that values are interpolated to. 
        END IF
        IF (ng.eq.2) THEN
          WESTERN_EDGE_NESTING=.TRUE.
          ng1=1    !other grid, ie gridded data.
          ng2=2    !local grid, grid that values are interpolated to. 
        END IF
# endif
        ILB1=LBOUND(field4d(ng1)%angler,DIM=1)
        IUB1=UBOUND(field4d(ng1)%angler,DIM=1)
        JLB1=LBOUND(field4d(ng1)%angler,DIM=2)
        JUB1=UBOUND(field4d(ng1)%angler,DIM=2)
        ILB2=LBOUND(field4d(ng2)%angler,DIM=1)
        IUB2=UBOUND(field4d(ng2)%angler,DIM=1)
        JLB2=LBOUND(field4d(ng2)%angler,DIM=2)
        JUB2=UBOUND(field4d(ng2)%angler,DIM=2)
!       KLB1=field4d(ng1)%LBk
!       KUB1=field4d(ng1)%UBk
!       KLB2=field4d(ng2)%LBk
!       KUB2=field4d(ng2)%UBk
        KUB1=N(ng)
        KUB2=N(ng)
        IF (field4d(ng2)%gtype.eq.w3dvar) THEN
          KLB1=0
          KLB2=0
        ELSE 
          KLB1=1
          KLB2=1
        END IF

!       LLB1=field4d(ng1)%LBl
!       LUB1=field4d(ng1)%UBl
!       LLB2=field4d(ng2)%LBl
!       LUB2=field4d(ng2)%UBl
        LLB1=1
        LUB1=NT(ng)
        LLB2=1
        LUB2=NT(ng)

        IJspv=0.0_r8
        IF (WESTERN_EDGE_NESTING) THEN
          Ms=ILB2
          Me=0
          Ns=JLB2
          Ne=JUB2
        END IF
        IF (EASTERN_EDGE_NESTING) THEN
          Ms=IUB2-GHOST_POINTS+1
          Me=IUB2
          Ns=JLB2
          Ne=JUB2
        END IF
        allocate ( Ipos(ILB2:IUB2,JLB2:JUB2) )
        allocate ( Jpos(ILB2:IUB2,JLB2:JUB2) )

        do_hindices=.FALSE.
        IF (field4d(ng2)%gtype.eq.u3dvar) THEN
          IF ( OCEAN(ng2)%hind_init(u3dvar) .eq. .FALSE.) THEN
            OCEAN(ng2)%hind_init(u3dvar)=.TRUE.
            do_hindices=.TRUE.
          ELSE
            IPOS=OCEAN(ng)%hind4di_u
            JPOS=OCEAN(ng)%hind4dj_u
          END IF
        ELSE IF (field4d(ng2)%gtype.eq.v3dvar) THEN
          IF ( OCEAN(ng2)%hind_init(v3dvar) .eq. .FALSE.) THEN
            OCEAN(ng2)%hind_init(v3dvar)=.TRUE.
            do_hindices=.TRUE.
          ELSE
            IPOS=OCEAN(ng)%hind4di_v
            JPOS=OCEAN(ng)%hind4dj_v
          END IF
        ELSE               !r4dvar
          IF ( OCEAN(ng2)%hind_init(r3dvar) .eq. .FALSE.) THEN
            OCEAN(ng2)%hind_init(r3dvar)=.TRUE.
            do_hindices=.TRUE.
          ELSE
            IPOS=OCEAN(ng)%hind3di_r
            JPOS=OCEAN(ng)%hind3dj_r
          END IF
        END IF
        IF (do_hindices) THEN
          CALL hindices (ng, ILB1, IUB1, JLB1, JUB1,                    &
     &                   ILB1, IUB1, JLB1, JUB1,                        &
     &                   field4d(ng1)%angler,                           &
     &                   field4d(ng1)%x,   field4d(ng1)%y,              &
     &                   ILB2, IUB2, JLB2, JUB2,                        &
     &                   Ms, Me, Ns, Ne,                                &
     &                   field4d(ng2)%x,   field4d(ng2)%y,              &
     &                   Ipos, Jpos,                                    &
     &                   IJspv, rectangular, Xmin, Xmax, Ymin, Ymax)
          IF (field4d(ng2)%gtype.eq.u3dvar) THEN
            OCEAN(ng)%hind4di_u=IPOS
            OCEAN(ng)%hind4dj_u=JPOS
          ELSE IF (field4d(ng2)%gtype.eq.v3dvar) THEN
            OCEAN(ng)%hind4di_v=IPOS
            OCEAN(ng)%hind4dj_v=JPOS
          ELSE               !r3dvar
            OCEAN(ng)%hind3di_r=IPOS
            OCEAN(ng)%hind3dj_r=JPOS
          END IF
        END IF
        CALL linterp4d (ng, ILB1, IUB1, JLB1, JUB1, KLB1, KUB1, LLB1,   &
     &                     LUB1,                                        &
     &                   field4d(ng1)%x,   field4d(ng1)%y,              &
     &                   field4d(ng1)%s,                                &
     &                   ILB2, IUB2, JLB2, JUB2, KLB2, KUB2, LLB2, LUB2,&
     &                   Ms, Me, Ns, Ne,                                &
     &                   Ipos,Jpos,                                     &
     &                   field4d(ng2)%x,   field4d(ng2)%y,              &
     &                   field4d(ng2)%s,                                &
     &                   my_min, my_max)
        deallocate ( Ipos, Jpos )
      END DO

      RETURN
      END SUBROUTINE interp_nesting4d
#endif
      SUBROUTINE allocate_nesting2d (field2d, id, gtype, Tindex)
!
!=======================================================================
!                                                                      !
!  This routine allocates and initializes nesting structure for 2D     !
!  state variables.                                                    !
!                                                                      !
!=======================================================================
!
      USE mod_param
# ifdef SOLVE3D
      USE mod_coupling
# endif
      USE mod_forces
      USE mod_grid
      USE mod_mixing
      USE mod_ncparam
      USE mod_ocean
      USE mod_scalars
      USE mod_nesting
!
!  Imported variable declarations.
!
      integer, target :: id, gtype
      integer, intent(in), optional :: Tindex

      TYPE (T_NEST2D), pointer :: field2d(:)
!
!  Local variable declarations.
!
      integer :: ng
!
!-----------------------------------------------------------------------
!  Generic 2D structure.
!-----------------------------------------------------------------------
!
!  Allocate.
!
      IF (ASSOCIATED(field2d)) THEN
        deallocate ( field2d )
      END IF
      allocate ( field2d(Ngrids) )
!
!  Load field descriptors.
!
      DO ng=1,Ngrids
!
!  Load id, grid type, tile, array bounds, and starting/ending
!  computational indices.
!
        field2d(ng) % id    => id
        field2d(ng) % gtype => gtype
        field2d(ng) % tile  => BOUNDS(ng) % tile
        field2d(ng) % LBi   => BOUNDS(ng) % LBi
        field2d(ng) % UBi   => BOUNDS(ng) % UBi
        field2d(ng) % LBj   => BOUNDS(ng) % LBj
        field2d(ng) % UBj   => BOUNDS(ng) % UBj
        field2d(ng) % Istr  => BOUNDS(ng) % Istr
        field2d(ng) % Iend  => BOUNDS(ng) % Iend
        field2d(ng) % Jstr  => BOUNDS(ng) % Jstr
        field2d(ng) % Jend  => BOUNDS(ng) % Jend
        field2d(ng) % IstrR => BOUNDS(ng) % IstrR
        field2d(ng) % IstrU => BOUNDS(ng) % IstrU
        field2d(ng) % IendR => BOUNDS(ng) % IendR
        field2d(ng) % JstrR => BOUNDS(ng) % JstrR
        field2d(ng) % JstrV => BOUNDS(ng) % JstrV
        field2d(ng) % JendR => BOUNDS(ng) % JendR
!
!  Associate the appropriate grid arrays.
!
        field2d(ng) % angler => GRID(ng) % angler
        IF (spherical) THEN
          IF (gtype.eq.u2dvar) THEN
            field2d(ng) % x => GRID(ng) % lonu
            field2d(ng) % y => GRID(ng) % latu
          ELSE IF (gtype.eq.v2dvar) THEN
            field2d(ng) % x => GRID(ng) % lonv
            field2d(ng) % y => GRID(ng) % latv
          ELSE
            field2d(ng) % x => GRID(ng) % lonr
            field2d(ng) % y => GRID(ng) % latr
          END IF
        ELSE
          IF (gtype.eq.u2dvar) THEN
            field2d(ng) % x => GRID(ng) % xu
            field2d(ng) % y => GRID(ng) % yu
          ELSE IF (gtype.eq.v2dvar) THEN
            field2d(ng) % x => GRID(ng) % xv
            field2d(ng) % y => GRID(ng) % yv
          ELSE
            field2d(ng) % x => GRID(ng) % xr
            field2d(ng) % y => GRID(ng) % yr
          END IF
        END IF

# ifdef MASKING
!
!  Associate the appropriate Land/Sea mask.
!
        IF (gtype.eq.u2dvar) THEN
          field2d(ng) % mask  => GRID(ng) % umask
        ELSE IF (gtype.eq.v2dvar) THEN
          field2d(ng) % mask  => GRID(ng) % vmask
        ELSE
          field2d(ng) % mask  => GRID(ng) % rmask
        END IF
# endif
!
!  Associate the appropriate state 2D array.
!
        IF (id.eq.idFsur) THEN
          field2d(ng) % s => OCEAN(ng) % zeta(:,:,Tindex)
        ELSE IF (id.eq.idUbar) THEN
          field2d(ng) % s => OCEAN(ng) % ubar(:,:,Tindex)
        ELSE IF (id.eq.idVbar) THEN
          field2d(ng) % s => OCEAN(ng) % vbar(:,:,Tindex)
        ELSE IF (id.eq.idRzet) THEN
          field2d(ng) % s => OCEAN(ng) % rzeta(:,:,Tindex)
        ELSE IF (id.eq.idUbms) THEN
          field2d(ng) % s => FORCES(ng) % bustr(:,:)
        ELSE IF (id.eq.idVbms) THEN
          field2d(ng) % s => FORCES(ng) % bvstr(:,:)
# ifdef VAR_RHO_2D
        ELSE IF (id.eq.510) THEN
          field2d(ng) % s => COUPLING(ng) % rhoA(:,:)
        ELSE IF (id.eq.511) THEN
          field2d(ng) % s => COUPLING(ng) % rhoS(:,:)
# endif
# ifdef SOLVE3D
        ELSE IF (id.eq.idZavg) THEN
          field2d(ng) % s => COUPLING(ng) % Zt_avg1
        ELSE IF (id.eq.idUfx1) THEN
          field2d(ng) % s => COUPLING(ng) % DU_avg1
        ELSE IF (id.eq.idVfx1) THEN
          field2d(ng) % s => COUPLING(ng) % DV_avg1
# endif
        END IF

      END DO

      RETURN
      END SUBROUTINE allocate_nesting2d

# ifdef SOLVE3D

      SUBROUTINE allocate_nesting3d (field3d, id, gtype, Tindex)
!
!=======================================================================
!                                                                      !
!  This routine allocates and initializes nesting structure for 3D     !
!  state variables.                                                    !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_grid
      USE mod_mixing
      USE mod_ncparam
      USE mod_ocean
      USE mod_scalars
      USE mod_nesting
!
!  Imported variable declarations.
!
      integer, target :: id, gtype

      integer, intent(in), optional :: Tindex

      TYPE (T_NEST3D), pointer :: field3d(:)
!
!  Local variable declarations.
!
      integer :: ng
      integer, target :: LBk, UBk, Kstr, Kend
!
!-----------------------------------------------------------------------
!  Generic 3D structure.
!-----------------------------------------------------------------------
!
!  Allocate.
!
      IF (ASSOCIATED(field3d)) THEN
        deallocate ( field3d )
      END IF
      allocate ( field3d(Ngrids) )
!
!  Load field descriptors.
!
      DO ng=1,Ngrids
!       gtype=Iinfo(1,id,ng)
!
!  Load id, grid type, tile, array bounds, and starting/ending
!  computational indices.
!
        field3d(ng) % id    => id
        field3d(ng) % gtype => gtype
        field3d(ng) % tile  => BOUNDS(ng) % tile
        field3d(ng) % LBi   => BOUNDS(ng) % LBi
        field3d(ng) % UBi   => BOUNDS(ng) % UBi
        field3d(ng) % LBj   => BOUNDS(ng) % LBj
        field3d(ng) % UBj   => BOUNDS(ng) % UBj
        field3d(ng) % Istr  => BOUNDS(ng) % Istr
        field3d(ng) % Iend  => BOUNDS(ng) % Iend
        field3d(ng) % Jstr  => BOUNDS(ng) % Jstr
        field3d(ng) % Jend  => BOUNDS(ng) % Jend
        field3d(ng) % IstrR => BOUNDS(ng) % IstrR
        field3d(ng) % IstrU => BOUNDS(ng) % IstrU
        field3d(ng) % IendR => BOUNDS(ng) % IendR
        field3d(ng) % JstrR => BOUNDS(ng) % JstrR
        field3d(ng) % JstrV => BOUNDS(ng) % JstrV
        field3d(ng) % JendR => BOUNDS(ng) % JendR
!
!  Load third-dimension bounds.
!
        IF (gtype.eq.w3dvar) THEN
          LBk=0
          UBk=N(ng)
          Kstr=0
          Kend=N(ng)
          field3d(ng) % LBk   => LBk
          field3d(ng) % UBk   => UBk
          field3d(ng) % Kstr  => Kstr
          field3d(ng) % Kend  => Kend
        ELSE
          LBk=1
          UBk=N(ng)
          Kstr=1
          Kend=N(ng)
          field3d(ng) % LBk   => LBk
          field3d(ng) % UBk   => UBk
          field3d(ng) % Kstr  => Kstr
          field3d(ng) % Kend  => Kend
        END IF
!
!  Associate the appropriate grid arrays.
!
        field3d(ng) % angler => GRID(ng) % angler
        IF (spherical) THEN
          IF (gtype.eq.u3dvar) THEN
            field3d(ng) % x => GRID(ng) % lonu
            field3d(ng) % y => GRID(ng) % latu
          ELSE IF (gtype.eq.v3dvar) THEN
            field3d(ng) % x => GRID(ng) % lonv
            field3d(ng) % y => GRID(ng) % latv
          ELSE
            field3d(ng) % x => GRID(ng) % lonr
            field3d(ng) % y => GRID(ng) % latr
          END IF
        ELSE
          IF (gtype.eq.u3dvar) THEN
            field3d(ng) % x => GRID(ng) % xu
            field3d(ng) % y => GRID(ng) % yu
          ELSE IF (gtype.eq.v3dvar) THEN
            field3d(ng) % x => GRID(ng) % xv
            field3d(ng) % y => GRID(ng) % yv
          ELSE
            field3d(ng) % x => GRID(ng) % xr
            field3d(ng) % y => GRID(ng) % yr
          END IF
        END IF

#  ifdef MASKING
!
!  Associate the appropriate Land/Sea mask.
!
        IF (gtype.eq.u3dvar) THEN
          field3d(ng) % mask  => GRID(ng) % umask
        ELSE IF (gtype.eq.v3dvar) THEN
          field3d(ng) % mask  => GRID(ng) % vmask
        ELSE
          field3d(ng) % mask  => GRID(ng) % rmask
        END IF
#  endif
!
!  Associate the appropriate state 3D array.
!
        IF (id.eq.idUvel) THEN
          field3d(ng) % s => OCEAN(ng) % u(:,:,:,Tindex)
        ELSE IF (id.eq.idVvel) THEN
          field3d(ng) % s => OCEAN(ng) % v(:,:,:,Tindex)
        ELSE IF (id.eq.idRu3d) THEN
          field3d(ng) % s => OCEAN(ng) % ru(:,:,:,Tindex)
        ELSE IF (id.eq.idRv3d) THEN
          field3d(ng) % s => OCEAN(ng) % rv(:,:,:,Tindex)
        ELSE IF (id.eq.idOvel) THEN
          field3d(ng) % s => OCEAN(ng) % W
        ELSE IF (id.eq.idWvel) THEN
          field3d(ng) % s => OCEAN(ng) % wvel
        ELSE IF (id.eq.idHUav) THEN
          field3d(ng) % s => GRID(ng) % Huon
        ELSE IF (id.eq.idHVav) THEN
          field3d(ng) % s => GRID(ng) % Hvom
        ELSE IF (id.eq.idDano) THEN
          field3d(ng) % s => OCEAN(ng) % rho
        ELSE IF (id.eq.idVvis) THEN
          field3d(ng) % s => MIXING(ng) % Akv
#  if defined GLS_MIXING || defined MY25_MIXING
        ELSE IF (id.eq.idMtke) THEN
          field3d(ng) % s => MIXING(ng) % tke(:,:,:,Tindex)
        ELSE IF (id.eq.idMtls) THEN
          field3d(ng) % s => MIXING(ng) % gls(:,:,:,Tindex)
        ELSE IF (id.eq.idVmKK) THEN
          field3d(ng) % s => MIXING(ng) % Akk
#   ifdef GLS_MIXING
        ELSE IF (id.eq.idVvel) THEN
          field3d(ng) % s => MIXING(ng) % Akp
#   endif
#  endif
        END IF

      END DO

      RETURN
      END SUBROUTINE allocate_nesting3d

      SUBROUTINE allocate_nesting4d (field4d, id, gtype, Tindex)
!
!=======================================================================
!                                                                      !
!  This routine allocates and initializes nesting structure for 4D     !
!  state variables.                                                    !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_grid
      USE mod_mixing
      USE mod_ncparam
      USE mod_ocean
      USE mod_scalars
      USE mod_nesting
!
!  Imported variable declarations.
!
      integer, target :: id, gtype

      integer, intent(in), optional :: Tindex

      TYPE (T_NEST4D), pointer :: field4d(:)
!
!  Local variable declarations.
!
      integer, target :: LBk, UBk, LBl, UBl
      integer, target :: Kstr, Kend, Lstr, Lend
      integer :: ng
!      integer, target :: gtype
!
!-----------------------------------------------------------------------
!  Generic 4D structure.
!-----------------------------------------------------------------------
!
!  Allocate.
!
      IF (associated(field4d)) THEN
        deallocate ( field4d )
      END IF
      allocate ( field4d(Ngrids) )
!
!  Load field descriptors.
!
      DO ng=1,Ngrids
!        gtype=Iinfo(1,id,ng)
!
!  Load id, grid type, tile, array bounds, and starting/ending
!  computational indices.
!
        field4d(ng) % id    => id
        field4d(ng) % gtype => gtype
        field4d(ng) % tile  => BOUNDS(ng) % tile
        field4d(ng) % LBi   => BOUNDS(ng) % LBi
        field4d(ng) % UBi   => BOUNDS(ng) % UBi
        field4d(ng) % LBj   => BOUNDS(ng) % LBj
        field4d(ng) % UBj   => BOUNDS(ng) % UBj
        field4d(ng) % Istr  => BOUNDS(ng) % Istr
        field4d(ng) % Iend  => BOUNDS(ng) % Iend
        field4d(ng) % Jstr  => BOUNDS(ng) % Jstr
        field4d(ng) % Jend  => BOUNDS(ng) % Jend
        field4d(ng) % IstrR => BOUNDS(ng) % IstrR
        field4d(ng) % IstrU => BOUNDS(ng) % IstrU
        field4d(ng) % IendR => BOUNDS(ng) % IendR
        field4d(ng) % JstrR => BOUNDS(ng) % JstrR
        field4d(ng) % JstrV => BOUNDS(ng) % JstrV
        field4d(ng) % JendR => BOUNDS(ng) % JendR
!
!  Load third-dimension bounds.
!
        IF (gtype.eq.r3dvar) THEN
          LBk=1
          UBk=N(ng)
          Kstr=1
          Kend=N(ng)
        ELSE IF (gtype.eq.w3dvar) THEN
          LBk=0
          UBk=N(ng)
          Kstr=0
          Kend=N(ng)
        END IF
        field4d(ng) % LBk   => LBk
        field4d(ng) % UBk   => UBk
        field4d(ng) % Kstr  => Kstr
        field4d(ng) % Kend  => Kend
!
!  Associate the appropriate grid arrays.
!
        field4d(ng) % angler => GRID(ng) % angler
        IF (spherical) THEN
          IF (gtype.eq.u3dvar) THEN
            field4d(ng) % x => GRID(ng) % lonu
            field4d(ng) % y => GRID(ng) % latu
          ELSE IF (gtype.eq.v3dvar) THEN
            field4d(ng) % x => GRID(ng) % lonv
            field4d(ng) % y => GRID(ng) % latv
          ELSE
            field4d(ng) % x => GRID(ng) % lonr
            field4d(ng) % y => GRID(ng) % latr
          END IF
        ELSE
          IF (gtype.eq.u3dvar) THEN
            field4d(ng) % x => GRID(ng) % xu
            field4d(ng) % y => GRID(ng) % yu
          ELSE IF (gtype.eq.v3dvar) THEN
            field4d(ng) % x => GRID(ng) % xv
            field4d(ng) % y => GRID(ng) % yv
          ELSE
            field4d(ng) % x => GRID(ng) % xr
            field4d(ng) % y => GRID(ng) % yr
          END IF
        END IF

#  ifdef MASKING
!
!  Associate the appropriate Land/Sea mask.
!
        IF (gtype.eq.u3dvar) THEN
          field4d(ng) % mask  => GRID(ng) % umask
        ELSE IF (gtype.eq.v3dvar) THEN
          field4d(ng) % mask  => GRID(ng) % vmask
        ELSE
          field4d(ng) % mask  => GRID(ng) % rmask
        END IF
#  endif
!
!  Associate the appropriate state 4D array.  Notice that temperature
!  triggers associating all tracers.
!
        IF (id.eq.idTvar(itemp)) THEN
          LBl=1
          UBl=NT(ng)
          Lstr=1
          Lend=NT(ng)
          field4d(ng) % LBl   => LBl
          field4d(ng) % UBl   => UBl
          field4d(ng) % Lstr  => Lstr
          field4d(ng) % Lend  => Lend
          field4d(ng) % s     => OCEAN(ng) % t(:,:,:,Tindex,:)
        ELSE IF (id.eq.idTdif) THEN       
          LBl=1
          UBl=NT(ng)
          Lstr=1
          Lend=NT(ng)
          field4d(ng) % LBl   => LBl
          field4d(ng) % UBl   => UBl
          field4d(ng) % Lstr  => Lstr
          field4d(ng) % Lend  => Lend
          field4d(ng) % s     => MIXING(ng) % AKt(:,:,:,:)
        END IF

      END DO

      RETURN
      END SUBROUTINE allocate_nesting4d
# endif

#endif
      END MODULE nesting_mod
