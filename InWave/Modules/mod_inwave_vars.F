#include "cppdefs.h"
      MODULE mod_inwave_vars
!
!svn $Id: driver_inwave.F 732 2008-09-07 01:55:51Z jcwarner $
! LAST CHANGE: mai 12/28/2010
!
!=======================================================================
!                                                                      !
!  AC        Action density                                            !
!  cx        wave celerity xi direction                                !
!  cy        wave celerity eta direction                               !
!  ct        wave celerity direction direction                         !
!  pd        wave direction bins, per radians                          !
!  Tr        relative short wave period                                !
!  Ta        absolute short wave period                                !
!  kwc       wave number in wave current flows                         !
!  cwc       wave celerity in wave current flows                       !
!  h_tot     total water depth at two different instants               !
!  u_rho     xi component of the currents affecting waves              !
!  v_rho     etai component of the currents affecting waves            !
!  wd        wave directions                                           !
!  pd        incremental wave directions for each directional bin      !
!  Trep      representative absolute boundary wave period.
!                                                                      !
!=======================================================================
!
        USE mod_kinds
        implicit none

        TYPE T_INWAVEP
          real(r8), pointer :: AC(:,:,:,:)
          real(r8), pointer :: cx(:,:,:)
          real(r8), pointer :: cy(:,:,:)
          real(r8), pointer :: ct(:,:,:)
          real(r8), pointer :: Tr(:,:,:)
          real(r8), pointer :: Ta(:,:,:)
          real(r8), pointer :: kwc(:,:,:)
          real(r8), pointer :: cwc(:,:,:)
          real(r8), pointer :: h_tot(:,:)
          real(r8), pointer :: u_rho(:,:)
          real(r8), pointer :: v_rho(:,:)
        END TYPE T_INWAVEP
        TYPE (T_INWAVEP), allocatable :: WAVEP(:)
!
        TYPE T_INWAVEG
          integer, pointer :: Inwavecircle
          real(r8), pointer :: wd(:)
          real(r8), pointer :: pd
          real(r8), pointer :: Trep

        END TYPE T_INWAVEG
        TYPE (T_INWAVEG), allocatable :: WAVEG(:)

      CONTAINS
!
!***********************************************************************
      SUBROUTINE allocate_inwave_vars (ng, LBi, UBi, LBj, UBj)
!***********************************************************************
!
!======================================================================!
!                                                                      !
!  This routine allocates all variables in the module for all nested   !
!  grids.                                                              !
!                                                                      !
!======================================================================!
!
      USE mod_inwave_params
!
!  Local variable declarations.
!
      integer, intent(in) :: ng, LBi, UBi, LBj, UBj
!
!-----------------------------------------------------------------------
!  Allocate and initialize module variables.
!-----------------------------------------------------------------------
!
      IF (ng.eq.1) allocate ( WAVEP(Ngrids) )
      IF (ng.eq.1) allocate ( WAVEG(Ngrids) )

      allocate ( WAVEP(ng) % AC(LBi:UBi,LBj:UBj,ND,3) )
      allocate ( WAVEP(ng) % cx(LBi:UBi,LBj:UBj,ND) )
      allocate ( WAVEP(ng) % cy(LBi:UBi,LBj:UBj,ND) )
      allocate ( WAVEP(ng) % ct(LBi:UBi,LBj:UBj,ND+1) )
      allocate ( WAVEP(ng) % Tr(LBi:UBi,LBj:UBj,ND) )
      allocate ( WAVEP(ng) % Ta(LBi:UBi,LBj:UBj,ND) )
      allocate ( WAVEP(ng) % kwc(LBi:UBi,LBj:UBj,ND) )
      allocate ( WAVEP(ng) % cwc(LBi:UBi,LBj:UBj,ND) )
      allocate ( WAVEP(ng) % h_tot(LBi:UBi,LBj:UBj) )
      allocate ( WAVEP(ng) % u_rho(LBi:UBi,LBj:UBj) )
      allocate ( WAVEP(ng) % v_rho(LBi:UBi,LBj:UBj) )

      allocate ( WAVEG(ng) % Inwavecircle )
      allocate ( WAVEG(ng) % wd(ND) )
      allocate ( WAVEG(ng) % pd )
      allocate ( WAVEG(ng) % Trep )

      RETURN
      END SUBROUTINE allocate_inwave_vars
!
!***********************************************************************
      SUBROUTINE initialize_inwave_vars (ng, LBi, UBi, LBj, UBj, tile,  &
     &                                   IniRec)
!***********************************************************************
!
      USE mod_kinds
      USE mod_grid
      USE mod_param
      USE mod_scalars
      USE mod_inwave_params
!     USE mod_inwave_bound
      USE inwave_iounits
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, IniRec
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
!  Local variable declarations.
!
      integer :: Imin, Imax, Jmin, Jmax
      integer :: i, j, d

      real(r8), parameter :: IniVal = 0.0_r8
      real(r8), allocatable :: ndirs(:)
      real(r8), allocatable :: pdirs(:)

      real(r8) :: cff
      real(r8) :: cnvrad
      real(r8) :: twopi
      real(r8), allocatable :: ac_ini(:,:,:)
      real(r8), allocatable :: cx_ini(:,:,:)
      real(r8), allocatable :: cy_ini(:,:,:)
      real(r8), allocatable :: ct_ini(:,:,:)
      real(r8), allocatable :: ta_ini(:,:,:)
!     real(r8) :: Ta_ini
      real(r8) :: end_ang

      character (len=80) :: ncini_name

#include "set_bounds.h"
!
!  Set array initialization range.
!
#ifdef DISTRIBUTE
      Imin=BOUNDS(ng)%LBi(tile)
      Imax=BOUNDS(ng)%UBi(tile)
      Jmin=BOUNDS(ng)%LBj(tile)
      Jmax=BOUNDS(ng)%UBj(tile)
#else
      IF (DOMAIN(ng)%Western_Edge(tile)) THEN
        Imin=BOUNDS(ng)%LBi(tile)
      ELSE
        Imin=Istr
      END IF
      IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
        Imax=BOUNDS(ng)%UBi(tile)
      ELSE
        Imax=Iend
      END IF
      IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
        Jmin=BOUNDS(ng)%LBj(tile)
      ELSE
        Jmin=Jstr
      END IF
      IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
        Jmax=BOUNDS(ng)%UBj(tile)
      ELSE
        Jmax=Jend
      END IF
#endif
!
!-----------------------------------------------------------------------
! Compute degree bins.
!-----------------------------------------------------------------------
!
      IF (.not.allocated(ndirs)) allocate (ndirs(ND))
      ncini_name=IWINIname(ng)
      call get_inwave_grid(ng, ncini_name, ndirs)
      cnvrad=pi/180.0_r8
      twopi=2.0_r8*pi
!
!-----------------------------------------------------------------------
! Pass the directions of the bins to radians
!-----------------------------------------------------------------------
!
      DO d=1,ND
        WAVEG(ng)%wd(d)=ndirs(d)*cnvrad
      END DO
!
!-----------------------------------------------------------------------
! Compute the angle increment
!-----------------------------------------------------------------------
!
      IF (ND.gt.1) THEN
        cff=WAVEG(ng)%wd(2)-WAVEG(ng)%wd(1)
        IF (cff.lt.0.0_r8) THEN
          cff=cff+2.0_r8*pi
        ENDIF
        WAVEG(ng)%pd=cff
      ELSE
        WAVEG(ng)%pd=2.0_r8*pi/16.0_r8
      END IF

!     If user sets computational grid for a full circle or just a sector.
      IF (INT(2.0_r8*pi/WAVEG(ng)%pd).eq.ND) THEN
        WAVEG(ng)%Inwavecircle=1
      ELSE
        WAVEG(ng)%Inwavecircle=0
      END IF

#if defined THETA_AC_PERIODIC_NOUSE
!
!-----------------------------------------------------------------------
! check if we are working on a segment or if we have to wrap around
!-----------------------------------------------------------------------
!
      end_ang=WAVEG(ng)%wd(ND)+WAVEG(ng)%pd
      IF(end_ang.ge.twopi) end_ang=end_ang-twopi
      IF(end_ang.ne.WAVEG(ng)%wd(1))THEN
        write (stdout,30)
        exit_flag=4
      ENDIF
#endif
!
! need to distribute pd and wd
!
      IF (allocated(ndirs)) deallocate(ndirs)
!
!-----------------------------------------------------------------------
!  Read initial condition.
!-----------------------------------------------------------------------

      allocate (ac_ini(LBi:UBi,LBj:UBj,ND))
      allocate (cx_ini(LBi:UBi,LBj:UBj,ND))
      allocate (cy_ini(LBi:UBi,LBj:UBj,ND))
      allocate (ct_ini(LBi:UBi,LBj:UBj,ND+1))
      allocate (ta_ini(LBi:UBi,LBj:UBj,ND))
!
      ncini_name=IWINIname(ng)
      call get_inwave_ini(ng, LBi, UBi, LBj, UBj, IniRec,               &
     &                    ncini_name, ac_ini, ta_ini)
!
!-----------------------------------------------------------------------
!  Initialize module variables.
!-----------------------------------------------------------------------
!
!  Nonlinear model state.
!
      DO j=Jmin,Jmax
        DO i=Imin,Imax
          WAVEP(ng) % u_rho(i,j)=IniVal
          WAVEP(ng) % v_rho(i,j)=IniVal
          WAVEP(ng) % h_tot(i,j)=GRID(ng) % h(i,j)
          DO d=1,ND
            WAVEP(ng) % AC(i,j,d,1)=ac_ini(i,j,d)
            WAVEP(ng) % AC(i,j,d,2)=ac_ini(i,j,d)
            WAVEP(ng) % AC(i,j,d,3)=ac_ini(i,j,d)
            WAVEP(ng) % cx(i,j,d)=IniVal
            WAVEP(ng) % cy(i,j,d)=IniVal
            WAVEP(ng) % ct(i,j,d)=IniVal
            WAVEP(ng) % Ta(i,j,d)=ta_ini(i,j,d)
            WAVEP(ng) % Tr(i,j,d)=ta_ini(i,j,d)
            WAVEP(ng) % kwc(i,j,d)=IniVal
            WAVEP(ng) % cwc(i,j,d)=IniVal
          END DO
          WAVEP(ng) % ct(i,j,ND+1)=IniVal
        END DO
      END DO
 
      IF (allocated(ac_ini)) deallocate(ac_ini)
      IF (allocated(cx_ini)) deallocate(cx_ini)
      IF (allocated(cy_ini)) deallocate(cy_ini)
      IF (allocated(ct_ini)) deallocate(ct_ini)
      IF (allocated(ta_ini)) deallocate(ta_ini)

 30   FORMAT (/,' THE DIRECTIONAL SPACE IS A SECTOR, NEED TO SELECT     &
     &        ANOTHER BOUNDARY CONDITION THAT IS NOT PERIODIC')

      RETURN
      END SUBROUTINE initialize_inwave_vars

      END MODULE mod_inwave_vars
