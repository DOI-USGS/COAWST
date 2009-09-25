#include "cppdefs.h"
      MODULE mod_inwave_vars
!
!svn $Id: mod_inwave.F 790 2008-10-14 19:32:07Z jcwarner $
!================================================== John C. Warner =====
!                                                                      !
!  AC        Action density                                            !
!  cx        wave celerity xi direction                                !
!  cy        wave celerity eta direction                               !
!  ct        wave celerity direction direction                         !
!  pd        wave direction bins, per degrees                          !
!  WD        wave directions                                           !
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
          real(r8), pointer :: AC_ini(:,:,:)
          real(r8), pointer :: cx_ini(:,:,:)
          real(r8), pointer :: cy_ini(:,:,:)
          real(r8), pointer :: ct_ini(:,:,:)
          real(r8) Ta
        END TYPE T_INWAVEP
        TYPE (T_INWAVEP), allocatable :: WAVEP(:)
!
        TYPE T_INWAVEG
          real(r8), pointer :: WD(:)
          real(r8), pointer :: pd(:)
        END TYPE T_INWAVEG
        TYPE (T_INWAVEG), allocatable :: WAVEG(:)

      CONTAINS

      SUBROUTINE allocate_inwave_vars (ng, LBi, UBi, LBj, UBj)
!
!=======================================================================
!                                                                      !
!  This routine allocates all variables in the module for all nested   !
!  grids.                                                              !
!                                                                      !
!=======================================================================
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
!
!  Nonlinear model state.
!
      allocate ( WAVEP(ng) % AC(LBi:UBi,LBj:UBj,ND(ng),3) )
      allocate ( WAVEP(ng) % cx(LBi:UBi,LBj:UBj,ND(ng)) )
      allocate ( WAVEP(ng) % cy(LBi:UBi,LBj:UBj,ND(ng)) )
      allocate ( WAVEP(ng) % ct(LBi:UBi,LBj:UBj,ND(ng)) )

      allocate (WAVEP(ng) % ac_ini(LBi:UBi,LBj:UBj,ND(ng)) )
      allocate (WAVEP(ng) % cx_ini(LBi:UBi,LBj:UBj,ND(ng)) )
      allocate (WAVEP(ng) % cy_ini(LBi:UBi,LBj:UBj,ND(ng)) )
      allocate (WAVEP(ng) % ct_ini(LBi:UBi,LBj:UBj,ND(ng)) )
!
      allocate ( WAVEG(ng) % WD(ND(ng)) )
      allocate ( WAVEG(ng) % pd(ND(ng)) )

      RETURN
      END SUBROUTINE allocate_inwave_vars

      SUBROUTINE initialize_inwave_vars (ng, LBi, UBi, LBj, UBj, tile)
!
!=======================================================================
!                                                                      !
!  This routine initialize all variables in the module using first     !
!  touch distribution policy. In shared-memory configuration, this     !
!  operation actually performs propagation of the  "shared arrays"     !
!  across the cluster, unless another policy is specified to           !
!  override the default.                                               !
!                                                                      !
!=======================================================================
!
      USE mod_kinds
      USE mod_param
      USE mod_scalars
      USE mod_inwave_params
      USE mod_inwave_bound
      USE mod_iounits
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
      integer :: Imin, Imax, Jmin, Jmax
      integer :: i, j, d

      real(r8), parameter :: IniVal = 0.0_r8
      real(r8), allocatable :: ndirs(:)
      real(r8), allocatable :: pdirs(:)
      real(r8), parameter :: eps = 1.0E-16_r8

      real(r8) :: cff, cff1, cff2, cff3, cff4
      real(r8), allocatable :: ac_ini(:,:,:)
      real(r8), allocatable :: cx_ini(:,:,:)
      real(r8), allocatable :: cy_ini(:,:,:)
      real(r8), allocatable :: ct_ini(:,:,:)
      real(r8) :: Ta

      character (len=80) :: ncini_name
!
#include "tile.h"
!
!  Set array initialization range.
!
#ifdef _OPENMP
      IF (WESTERN_EDGE) THEN
        Imin=LBi
      ELSE
        Imin=Istr
      END IF
      IF (EASTERN_EDGE) THEN
        Imax=UBi
      ELSE
        Imax=Iend
      END IF
      IF (SOUTHERN_EDGE) THEN
        Jmin=LBj
      ELSE
        Jmin=Jstr
      END IF
      IF (NORTHERN_EDGE) THEN
        Jmax=UBj
      ELSE
        Jmax=Jend
      END IF
#else
      Imin=LBi
      Imax=UBi
      Jmin=LBj
      Jmax=UBj
#endif

!
! Compute degree bins. This may need to be moved.
!
      IF (.not.allocated(ndirs)) allocate (ndirs(ND(ng)))
!      IF (Master) THEN
        ncini_name=IWINIname(ng)
        call get_inwave_grid(ng, ncini_name, ndirs)
!      END IF

      DO d=1,ND(ng)
        WAVEG(ng)%WD(d)=ndirs(d)
      END DO

      DO d=1,ND(ng)
        if(d.gt.1)then
          cff=WAVEG(ng)%WD(d)-WAVEG(ng)%WD(d-1)
        else
          cff=WAVEG(ng)%WD(ND(ng))-WAVEG(ng)%WD(1)
        endif
        IF (cff.lt.0.0_r8) THEN
          cff=cff+360.0_r8
        END IF
        WAVEG(ng)%WD(d)=WAVEG(ng)%WD(d)*pi/180.
        WAVEG(ng)%pd(d)=180.0_r8/(pi*cff)
      END DO

! need to distribute pd and WD

      if (allocated(ndirs)) deallocate(ndirs)

!
!-----------------------------------------------------------------------
!  Read initial condition.
!-----------------------------------------------------------------------

      allocate (ac_ini(LBi:UBi,LBj:UBj,ND(ng)))
      allocate (cx_ini(LBi:UBi,LBj:UBj,ND(ng)))
      allocate (cy_ini(LBi:UBi,LBj:UBj,ND(ng)))
      allocate (ct_ini(LBi:UBi,LBj:UBj,ND(ng)))
!
!      IF (Master) THEN
        ncini_name=IWINIname(ng)
        call get_inwave_ini(ng, LBi, UBi, LBj, UBj, ncini_name, ac_ini, &
     &                      cx_ini, cy_ini, ct_ini, Ta)
!      END IF
!
!-----------------------------------------------------------------------
!  Initialize module variables.
!-----------------------------------------------------------------------
!
!  Nonlinear model state.
!
      DO j=Jmin,Jmax
        DO i=Imin,Imax
          DO d=1,ND(ng)
            WAVEP(ng) % AC(i,j,d,1) =AC_ini(i,j,d)
            WAVEP(ng) % AC(i,j,d,2) =AC_ini(i,j,d)
            WAVEP(ng) % AC(i,j,d,3) =AC_ini(i,j,d)
            WAVEP(ng) % cx(i,j,d) =cx_ini(i,j,d)
            WAVEP(ng) % cy(i,j,d) =cy_ini(i,j,d)
            WAVEP(ng) % ct(i,j,d) =ct_ini(i,j,d)
          END DO
        END DO
      END DO
      WAVEP(ng) % Ta= Ta

      RETURN
      END SUBROUTINE initialize_inwave_vars

      END MODULE mod_inwave_vars
