#include "cppdefs.h"
      MODULE mod_tides
#if defined SSH_TIDES || defined UV_TIDES
!
!git $Id$
!svn $Id: mod_tides.F 1099 2022-01-06 21:01:01Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2022 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  Tidal Components:                                                   !
!                                                                      !
!  Each of the following arrays has a dimension in tidal components    !
!  classified by period:                                               !
!                                                                      !
!    semi-diurnal:  M2, S2, N2, K2  (12.42, 12.00, 12.66, 11.97h)      !
!         diurnal:  K1, O1, P1, Q1  (23.93, 25.82, 24.07, 26.87h)      !
!                                                                      !
!  and other longer periods. The order of these tidal components is    !
!  irrelevant here.  The number of components to USE is depends on     !
!  the regional application.                                           !
!                                                                      !
!  CosOmega     Cosine tidal harmonics for current omega(t).           !
!  SinOmega     Sine tidal harmonics for current omega(t).             !
!  SSH_Tamp     Tidal elevation amplitude (m) at RHO-points.           !
!  SSH_Tphase   Tidal elevation phase (degrees/360) at RHO-points.     !
!  Tperiod      Tidal period (s).                                      !
!  UV_Tangle    Tidal current angle (radians; counterclockwise         !
!                 from EAST and rotated to curvilinear grid) at        !
!                 RHO-points.                                          !
!  UV_Tmajor    Maximum tidal current: tidal ellipse major axis        !
!                 (m/s) at RHO-points.                                 !
!  UV_Tminor    Minimum tidal current: tidal ellipse minor axis        !
!                 (m/s) at RHO-points.                                 !
!  UV_Tphase    Tidal current phase (degrees/360) at RHO-points.       !
!                                                                      !
# if defined AVERAGES && defined AVERAGES_DETIDE
!                                                                      !
!  Detided time-averaged fields via least-squares fitting. Notice that !
!  the harmonics for the state variable have an extra dimension of     !
!  size (0:2*NTC) to store several terms:                              !
!                                                                      !
!               index 0               mean term (accumulated sum)      !
!                     1:NTC           accumulated sine terms           !
!                     NTC+1:2*NTC     accumulated cosine terms         !
!                                                                      !
!  CosW_avg     Current time-average window COS(omega(k)*t).           !
!  CosW_sum     Time-accumulated COS(omega(k)*t).                      !
!  SinW_avg     Current time-average window SIN(omega(k)*t).           !
!  SinW_sum     Time-accumulated SIN(omega(k)*t).                      !
!  CosWCosW     Time-accumulated COS(omega(k)*t)*COS(omega(l)*t).      !
!  SinWSinW     Time-accumulated SIN(omega(k)*t)*SIN(omega(l)*t).      !
!  SinWCosW     Time-accumulated SIN(omega(k)*t)*COS(omega(l)*t).      !
!                                                                      !
!  ubar_detided Time-averaged and detided 2D u-momentum.               !
!  ubar_tide    Time-accumulated 2D u-momentum tide harmonics.         !
!  vbar_detided Time-averaged and detided 2D v-momentum.               !
!  vbar_tide    Time-accumulated 2D v-momentum tide harmonics.         !
!  zeta_detided Time-averaged and detided free-surface.                !
!  zeta_tide    Time-accumulated free-surface tide harmonics.          !
#  ifdef SOLVE3D
!  t_detided    Time-averaged and detided tracers (T,S).               !
!  t_tide       Time-accumulated 3D tracers (T,S) tide harmonics.      !
!  u_detided    Time-averaged and detided 3D u-momentum.               !
!  u_tide       Time-accumulated 3D u-momentum tide harmonics.         !
!  v_detided    Time-averaged and detided 3D v-momentum.               !
!  v_tide       Time-accumulated 3D v-momentum tide harmonics.         !
#  endif
!                                                                      !
# endif
# ifdef UV_WAVEDRAG
!  We need to accumulate a detided estimate of the bottom velocities   !
!  which is the subtracted from the bottom velocity before the         !
!  wave drag due to conversion of tidal energy into internal tides     !
!  is applied.                                                         !
!  sum_ubot      Time-accumulated bottom u-momentum.                   !
!  sum_vbot      Time-accumulated bottom v-momentum.                   !
!  filt_ubot     Time-mean bottom u-momentum.                          !
!  filt_vbot     Time-mean bottom v-momentum.                          !
# endif
!=======================================================================
!
        USE mod_kinds
        USE mod_param
        USE mod_stepping

        implicit none
!
        PUBLIC :: allocate_tides
        PUBLIC :: deallocate_tides
        PUBLIC :: initialize_tides
!
!-----------------------------------------------------------------------
!  Define T_TIDES structure.
!-----------------------------------------------------------------------
!
        TYPE T_TIDES

          real(r8), pointer :: Tperiod(:)
# if defined AVERAGES && defined AVERAGES_DETIDE
          real(r8), pointer :: CosOmega(:)
          real(r8), pointer :: SinOmega(:)
          real(r8), pointer :: CosW_avg(:)
          real(r8), pointer :: CosW_sum(:)
          real(r8), pointer :: SinW_avg(:)
          real(r8), pointer :: SinW_sum(:)
          real(r8), pointer :: CosWCosW(:,:)
          real(r8), pointer :: SinWSinW(:,:)
          real(r8), pointer :: SinWCosW(:,:)
# endif
# if defined SSH_TIDES
          real(r8), pointer :: SSH_Tamp(:,:,:)
          real(r8), pointer :: SSH_Tphase(:,:,:)
# endif
# if defined UV_TIDES
          real(r8), pointer :: UV_Tangle(:,:,:)
          real(r8), pointer :: UV_Tmajor(:,:,:)
          real(r8), pointer :: UV_Tminor(:,:,:)
          real(r8), pointer :: UV_Tphase(:,:,:)
# endif
# if defined AVERAGES && defined AVERAGES_DETIDE
          real(r8), pointer :: ubar_detided(:,:)
          real(r8), pointer :: ubar_tide(:,:,:)

          real(r8), pointer :: vbar_detided(:,:)
          real(r8), pointer :: vbar_tide(:,:,:)

          real(r8), pointer :: zeta_detided(:,:)
          real(r8), pointer :: zeta_tide(:,:,:)
#  ifdef SOLVE3D
          real(r8), pointer :: t_detided(:,:,:,:)
          real(r8), pointer :: t_tide(:,:,:,:,:)

          real(r8), pointer :: u_detided(:,:,:)
          real(r8), pointer :: u_tide(:,:,:,:)

          real(r8), pointer :: v_detided(:,:,:)
          real(r8), pointer :: v_tide(:,:,:,:)
#  endif
# endif
# ifdef UV_WAVEDRAG
          real(r8), pointer :: sum_ubot(:,:)
          real(r8), pointer :: sum_vbot(:,:)
          real(r8), pointer :: filt_ubot(:,:)
          real(r8), pointer :: filt_vbot(:,:)
# endif

        END TYPE T_TIDES
!
        TYPE (T_TIDES), allocatable :: TIDES(:)

# ifdef UV_WAVEDRAG
        integer :: nstep_hour
        integer :: hour_steps = 0
        real(r8), parameter :: drag_scale = 0.375
# endif

      CONTAINS
!
      SUBROUTINE allocate_tides (ng, LBi, UBi, LBj, UBj)
!
!=======================================================================
!                                                                      !
!  This routine allocates all variables in the module for all nested   !
!  grids.                                                              !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
# if defined PIO_LIB && defined DISTRIBUTE
      USE mod_pio_netcdf
# endif
      USE mod_scalars
      USE mod_stepping
!
      USE strings_mod, ONLY : FoundError
!
! Inported variable declarations.
!
      integer, intent(in) :: ng, LBi, UBi, LBj, UBj
!
!  Local variable declarations.
!
      logical :: foundit

      integer :: Nfiles, Vid, i, ifile, mg, nvatt, nvdim

      real(r8) :: size2d
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", allocate_tides"

# if defined PIO_LIB && defined DISTRIBUTE
!
      TYPE (Var_desc_t)  :: my_pioVar
# endif
!
!-----------------------------------------------------------------------
!  Allocate module variables.
!-----------------------------------------------------------------------
!
!  Inquire about the maximum number of tidal components. Notice that
!  currently we only support nested applications where the tidal
!  forcing is applied to the main coarser grid (RefineScale(ng)=0)
!  and the other grids get the tidal forcing from the contact areas.
!
      IF (LprocessTides(ng)) THEN
        MTC=0
        foundit=.FALSE.
!        SELECT CASE (TIDE(ng)%IOtype)
!          CASE (io_nf90)
            CALL netcdf_inq_var (ng, iNLM, TIDE(ng)%name,               &
     &                           MyVarName = Vname(1,idTper),           &
     &                           SearchVar = foundit,                   &
     &                           VarID = Vid,                           &
     &                           nVardim = nvdim,                       &
     &                           nVarAtt = nvatt)

# if defined PIO_LIB && defined DISTRIBUTE
          CASE (io_pio)
            CALL pio_netcdf_inq_var (ng, iNLM, TIDE(ng)%name,           &
     &                               MyVarName = Vname(1,idTper),       &
     &                               SearchVar = foundit,               &
     &                               pioVar = my_pioVar,                &
     &                               nVardim = nvdim,                   &
     &                               nVarAtt = nvatt)
# endif
!        END SELECT
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Set maximum number of tidal components.  Allocate and initialize
!  TIDE I/O structure. Notice that in nested applications, all the
!  nested grids need to have the same number of tidal component.
!
        IF (foundit) THEN
          MTC=MAX(MTC,var_Dsize(1))            ! first dimension
          DO mg=1,Ngrids
            NTC(mg)=var_Dsize(1)
          END DO
        END IF
      END IF
!
!  Allocate structure.
!
      IF (ng.eq.1) allocate ( TIDES(Ngrids) )
!
!  Set horizontal array size.
!
      size2d=REAL((UBi-LBi+1)*(UBj-LBj+1),r8)
!
!  Allocate tidal forcing variables.
!
      allocate ( TIDES(ng) % Tperiod(MTC)  )
      Dmem(ng)=Dmem(ng)+REAL(MTC,r8)

# if defined SSH_TIDES
      allocate ( TIDES(ng) % SSH_Tamp(LBi:UBi,LBj:UBj,MTC) )
      Dmem(ng)=Dmem(ng)+REAL(MTC,r8)*size2d

      allocate ( TIDES(ng) % SSH_Tphase(LBi:UBi,LBj:UBj,MTC) )
      Dmem(ng)=Dmem(ng)+REAL(MTC,r8)*size2d
# endif

# if defined UV_TIDES
      allocate ( TIDES(ng) % UV_Tangle(LBi:UBi,LBj:UBj,MTC) )
      Dmem(ng)=Dmem(ng)+REAL(MTC,r8)*size2d

      allocate ( TIDES(ng) % UV_Tmajor(LBi:UBi,LBj:UBj,MTC) )
      Dmem(ng)=Dmem(ng)+REAL(MTC,r8)*size2d

      allocate ( TIDES(ng) % UV_Tminor(LBi:UBi,LBj:UBj,MTC) )
      Dmem(ng)=Dmem(ng)+REAL(MTC,r8)*size2d

      allocate ( TIDES(ng) % UV_Tphase(LBi:UBi,LBj:UBj,MTC) )
      Dmem(ng)=Dmem(ng)+REAL(MTC,r8)*size2d
# endif

# if defined AVERAGES && defined AVERAGES_DETIDE
!
!  Allocate variables used for the least-squares detiding of
!  time-averaged fields.
!
      allocate ( TIDES(ng) % CosOmega(MTC) )
      Dmem(ng)=Dmem(ng)+REAL(MTC,r8)

      allocate ( TIDES(ng) % SinOmega(MTC) )
      Dmem(ng)=Dmem(ng)+REAL(MTC,r8)

      allocate ( TIDES(ng) % CosW_avg(MTC) )
      Dmem(ng)=Dmem(ng)+REAL(MTC,r8)

      allocate ( TIDES(ng) % CosW_sum(MTC) )
      Dmem(ng)=Dmem(ng)+REAL(MTC,r8)

      allocate ( TIDES(ng) % SinW_avg(MTC) )
      Dmem(ng)=Dmem(ng)+REAL(MTC,r8)

      allocate ( TIDES(ng) % SinW_sum(MTC) )
      Dmem(ng)=Dmem(ng)+REAL(MTC,r8)

      allocate ( TIDES(ng) % CosWCosW(MTC,MTC) )
      Dmem(ng)=Dmem(ng)+REAL(MTC*MTC,r8)

      allocate ( TIDES(ng) % SinWSinW(MTC,MTC) )
      Dmem(ng)=Dmem(ng)+REAL(MTC*MTC,r8)

      allocate ( TIDES(ng) % SinWCosW(MTC,MTC) )
      Dmem(ng)=Dmem(ng)+REAL(MTC*MTC,r8)

      IF (Aout(idFsuD,ng)) THEN
        allocate ( TIDES(ng) % zeta_detided(LBi:UBi,LBj:UBj) )
        Dmem(ng)=Dmem(ng)+size2d

        allocate ( TIDES(ng) % zeta_tide(LBi:UBi,LBj:UBj,0:2*MTC) )
        Dmem(ng)=Dmem(ng)+REAL(2*MTC+1,r8)*size2d
      END IF

      IF (Aout(idu2dD,ng)) THEN
        allocate ( TIDES(ng) % ubar_detided(LBi:UBi,LBj:UBj) )

        allocate ( TIDES(ng) % ubar_tide(LBi:UBi,LBj:UBj,0:2*MTC) )
        Dmem(ng)=Dmem(ng)+REAL(2*MTC+1,r8)*size2d
      END IF

      IF (Aout(idv2dD,ng)) THEN
        allocate ( TIDES(ng) % vbar_detided(LBi:UBi,LBj:UBj) )

        allocate ( TIDES(ng) % vbar_tide(LBi:UBi,LBj:UBj,0:2*MTC) )
        Dmem(ng)=Dmem(ng)+REAL(2*MTC+1,r8)*size2d
      END IF

#  ifdef SOLVE3D
      IF (Aout(idu3dD,ng)) THEN
        allocate ( TIDES(ng) % u_detided(LBi:UBi,LBj:UBj,N(ng)) )
        Dmem(ng)=Dmem(ng)+REAL(N(ng),r8)*size2d

        allocate ( TIDES(ng) % u_tide(LBi:UBi,LBj:UBj,N(ng),0:2*MTC) )
        Dmem(ng)=Dmem(ng)+REAL(N(ng)*(2*MTC+1),r8)*size2d
      END IF

      IF (Aout(idv3dD,ng)) THEN
        allocate ( TIDES(ng) % v_detided(LBi:UBi,LBj:UBj,N(ng)) )
        Dmem(ng)=Dmem(ng)+REAL(N(ng),r8)*size2d

        allocate ( TIDES(ng) % v_tide(LBi:UBi,LBj:UBj,N(ng),0:2*MTC) )
        Dmem(ng)=Dmem(ng)+REAL(N(ng)*(2*MTC+1),r8)*size2d
      END IF

      IF (ANY(Aout(idTrcD(:),ng))) THEN
        allocate ( TIDES(ng) % t_detided(LBi:UBi,LBj:UBj,N(ng),NAT) )
        Dmem(ng)=Dmem(ng)+REAL(N(ng)*NAT,r8)*size2d

        allocate ( TIDES(ng) % t_tide(LBi:UBi,LBj:UBj,N(ng),            &
     &                                0:2*MTC,NAT) )
        Dmem(ng)=Dmem(ng)+REAL(N(ng)*(2*MTC+1)*NAT,r8)*size2d
      END IF
#  endif
# endif
# ifdef UV_WAVEDRAG
      allocate ( TIDES(ng) % sum_ubot(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( TIDES(ng) % sum_vbot(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( TIDES(ng) % filt_ubot(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( TIDES(ng) % filt_vbot(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
# endif

      RETURN
      END SUBROUTINE allocate_tides
!
      SUBROUTINE deallocate_tides (ng)
!
!=======================================================================
!                                                                      !
!  This routine allocates all variables in the module for all nested   !
!  grids.                                                              !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_ncparam

# ifdef SUBOBJECT_DEALLOCATION
!
      USE destroy_mod, ONLY : destroy
# endif
!
! Inported variable declarations.
!
      integer, intent(in) :: ng
!
!  Local variable declarations.
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", deallocate_tides"


# ifdef SUBOBJECT_DEALLOCATION
!
!-----------------------------------------------------------------------
!  Deallocate each variable in the derived-type T_TIDES structure
!  separately.
!-----------------------------------------------------------------------
!
      IF (.not.destroy(ng, TIDES(ng)%Tperiod, MyFile,                   &
     &                 __LINE__, 'TIDES(ng)%Tperiod')) RETURN

#  if defined SSH_TIDES
      IF (.not.destroy(ng, TIDES(ng)%SSH_Tamp, MyFile,                  &
     &                 __LINE__, 'TIDES(ng)%SSH_Tamp')) RETURN

      IF (.not.destroy(ng, TIDES(ng)%SSH_Tphase, MyFile,                &
     &                 __LINE__, 'TIDES(ng)%SSH_Tphase')) RETURN
#  endif

#  if defined UV_TIDES
      IF (.not.destroy(ng, TIDES(ng)%UV_Tangle, MyFile,                 &
     &                 __LINE__, 'TIDES(ng)%UV_Tangle')) RETURN

      IF (.not.destroy(ng, TIDES(ng)%UV_Tmajor, MyFile,                 &
     &                 __LINE__, 'TIDES(ng)%UV_Tmajor')) RETURN

      IF (.not.destroy(ng, TIDES(ng)%UV_Tminor, MyFile,                 &
     &                 __LINE__, 'TIDES(ng)%UV_Tminor')) RETURN

      IF (.not.destroy(ng, TIDES(ng)%UV_Tphase, MyFile,                 &
     &                 __LINE__, 'TIDES(ng)%UV_Tphase')) RETURN
#  endif

#  if defined AVERAGES && defined AVERAGES_DETIDE
!
!  Allocate variables used for the least-squares detiding of
!  time-averaged fields.
!
      IF (.not.destroy(ng, TIDES(ng)%CosOmega, MyFile,                  &
     &                 __LINE__, 'TIDES(ng)%CosOmega')) RETURN

      IF (.not.destroy(ng, TIDES(ng)%SinOmega, MyFile,                  &
     &                 __LINE__, 'TIDES(ng)%SinOmega')) RETURN

      IF (.not.destroy(ng, TIDES(ng)%CosW_avg, MyFile,                  &
     &                 __LINE__, 'TIDES(ng)%CosW_avg')) RETURN

      IF (.not.destroy(ng, TIDES(ng)%CosW_sum, MyFile,                  &
     &                 __LINE__, 'TIDES(ng)%CosW_sum')) RETURN

      IF (.not.destroy(ng, TIDES(ng)%SinW_avg, MyFile,                  &
     &                 __LINE__, 'TIDES(ng)%SinW_avg')) RETURN

      IF (.not.destroy(ng, TIDES(ng)%SinW_sum, MyFile,                  &
     &                 __LINE__, 'TIDES(ng)%SinW_sum')) RETURN

      IF (.not.destroy(ng, TIDES(ng)%CosWCosW, MyFile,                  &
     &                 __LINE__, 'TIDES(ng)%CosWCosW')) RETURN

      IF (.not.destroy(ng, TIDES(ng)%SinWSinW, MyFile,                  &
     &                 __LINE__, 'TIDES(ng)%SinWSinW')) RETURN

      IF (.not.destroy(ng, TIDES(ng)%SinWCosW, MyFile,                  &
     &                 __LINE__, 'TIDES(ng)%SinWCosW')) RETURN

      IF (Aout(idFsuD,ng)) THEN
        IF (.not.destroy(ng, TIDES(ng)%zeta_detided, MyFile,            &
     &                   __LINE__, 'TIDES(ng)%zeta_detided')) RETURN

        IF (.not.destroy(ng, TIDES(ng)%zeta_tide, MyFile,               &
     &                   __LINE__, 'TIDES(ng)%zeta_tide')) RETURN
      END IF

      IF (Aout(idu2dD,ng)) THEN
        IF (.not.destroy(ng, TIDES(ng)%ubar_detided, MyFile,            &
     &                   __LINE__, 'TIDES(ng)%ubar_detided')) RETURN

        IF (.not.destroy(ng, TIDES(ng)%ubar_tide, MyFile,               &
     &                   __LINE__, 'TIDES(ng)%ubar_tide')) RETURN
      END IF

      IF (Aout(idv2dD,ng)) THEN
        IF (.not.destroy(ng, TIDES(ng)%vbar_detided, MyFile,            &
     &                   __LINE__, 'TIDES(ng)%vbar_detided')) RETURN

        IF (.not.destroy(ng, TIDES(ng)%vbar_tide, MyFile,               &
     &                   __LINE__, 'TIDES(ng)%vbar_tide')) RETURN
      END IF

#   ifdef SOLVE3D
      IF (Aout(idu3dD,ng)) THEN
        IF (.not.destroy(ng, TIDES(ng)%u_detided, MyFile,               &
     &                   __LINE__, 'TIDES(ng)%u_detided')) RETURN

        IF (.not.destroy(ng, TIDES(ng)%u_tide, MyFile,                  &
     &                   __LINE__, 'TIDES(ng)%u_tide')) RETURN
      END IF

      IF (Aout(idv3dD,ng)) THEN
        IF (.not.destroy(ng, TIDES(ng)%v_detided, MyFile,               &
     &                   __LINE__, 'TIDES(ng)%v_detided')) RETURN

        IF (.not.destroy(ng, TIDES(ng)%v_tide, MyFile,                  &
     &                   __LINE__, 'TIDES(ng)%v_tide')) RETURN
      END IF

      IF (ANY(Aout(idTrcD(:),ng))) THEN
        IF (.not.destroy(ng, TIDES(ng)%t_detided, MyFile,               &
     &                   __LINE__, 'TIDES(ng)%t_detided')) RETURN

        IF (.not.destroy(ng, TIDES(ng)%t_tide, MyFile,                  &
     &                   __LINE__, 'TIDES(ng)%t_tide')) RETURN
      END IF
#   endif
#  endif
# endif
!
!-----------------------------------------------------------------------
!  Deallocate derived-type TIDES structure.
!-----------------------------------------------------------------------
!
      IF (ng.eq.Ngrids) THEN
        IF (allocated(TIDES)) deallocate ( TIDES )
      END IF
!
      RETURN
      END SUBROUTINE deallocate_tides
!
      SUBROUTINE initialize_tides (ng, tile)
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
      USE mod_param
      USE mod_ncparam
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
      integer :: Imin, Imax, Jmin, Jmax
      integer :: i, itide, itrc, j, jtide, k

      real(r8), parameter :: IniVal = 0.0_r8

# include "set_bounds.h"
!
!  Set array initialization range.
!
# ifdef DISTRIBUTE
      Imin=BOUNDS(ng)%LBi(tile)
      Imax=BOUNDS(ng)%UBi(tile)
      Jmin=BOUNDS(ng)%LBj(tile)
      Jmax=BOUNDS(ng)%UBj(tile)
# else
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
# endif
!
!-----------------------------------------------------------------------
!  Initialize module variables.
!-----------------------------------------------------------------------
!
!  Initialize tidal forcing variables.
!
      IF (DOMAIN(ng)%SouthWest_Test(tile)) THEN
        DO itide=1,MTC
          TIDES(ng) % Tperiod(itide) = IniVal
        END DO
      END IF

      DO itide=1,MTC
# if defined SSH_TIDES
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            TIDES(ng) % SSH_Tamp(i,j,itide) = IniVal
            TIDES(ng) % SSH_Tphase(i,j,itide) = IniVal
          END DO
        END DO
# endif
# if defined UV_TIDES
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            TIDES(ng) % UV_Tangle(i,j,itide) = IniVal
            TIDES(ng) % UV_Tmajor(i,j,itide) = IniVal
            TIDES(ng) % UV_Tminor(i,j,itide) = IniVal
            TIDES(ng) % UV_Tphase(i,j,itide) = IniVal
          END DO
        END DO
# endif
      END DO

# if defined AVERAGES && defined AVERAGES_DETIDE
!
!  Initialize cariables used for the least-squares detiding of
!  time-averaged fields.
!
      IF (DOMAIN(ng)%SouthWest_Test(tile)) THEN
        DO jtide=1,MTC
          TIDES(ng) % CosOmega(jtide) = IniVal
          TIDES(ng) % SinOmega(jtide) = IniVal
          TIDES(ng) % CosW_avg(jtide) = IniVal
          TIDES(ng) % CosW_sum(jtide) = IniVal
          TIDES(ng) % SinW_avg(jtide) = IniVal
          TIDES(ng) % SinW_sum(jtide) = IniVal
          DO itide=1,MTC
            TIDES(ng) % CosWCosW(itide,jtide) = IniVal
            TIDES(ng) % SinWSinW(itide,jtide) = IniVal
            TIDES(ng) % SinWCosW(itide,jtide) = IniVal
          END DO
        END DO
      END IF

      IF (Aout(idFsuD,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            TIDES(ng) % zeta_detided(i,j) = IniVal
          END DO
        END DO
        DO itide=0,2*MTC
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              TIDES(ng) % zeta_tide(i,j,itide) = IniVal
            END DO
          END DO
        END DO
      END IF

      IF (Aout(idu2dD,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            TIDES(ng) % ubar_detided(i,j) = IniVal
          END DO
        END DO
        DO itide=0,2*MTC
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              TIDES(ng) % ubar_tide(i,j,itide) = IniVal
            END DO
          END DO
        END DO
      END IF

      IF (Aout(idv2dD,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            TIDES(ng) % vbar_detided(i,j) = IniVal
          END DO
        END DO
        DO itide=0,2*MTC
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              TIDES(ng) % vbar_tide(i,j,itide) = IniVal
            END DO
          END DO
        END DO
      END IF

#  ifdef SOLVE3D
      IF (Aout(idu3dD,ng)) THEN
        DO k=1,N(ng)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              TIDES(ng) % u_detided(i,j,k) = IniVal
            END DO
          END DO
        END DO
        DO itide=0,2*MTC
          DO k=1,N(ng)
            DO j=Jmin,Jmax
              DO i=Imin,Imax
                TIDES(ng) % u_tide(i,j,k,itide) = IniVal
              END DO
            END DO
          END DO
        END DO
      END IF

      IF (Aout(idv3dD,ng)) THEN
        DO k=1,N(ng)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              TIDES(ng) % v_detided(i,j,k) = IniVal
            END DO
          END DO
        END DO
        DO itide=0,2*MTC
          DO k=1,N(ng)
            DO j=Jmin,Jmax
              DO i=Imin,Imax
                TIDES(ng) % v_tide(i,j,k,itide) = IniVal
              END DO
            END DO
          END DO
        END DO
      END IF

      IF (ANY(Aout(idTrcD(:),ng))) THEN
        DO itrc=1,NAT
          DO k=1,N(ng)
            DO j=Jmin,Jmax
              DO i=Imin,Imax
                TIDES(ng) % t_detided(i,j,k,itrc) = IniVal
              END DO
            END DO
          END DO
        END DO
        DO itrc=1,NAT
          DO itide=0,2*MTC
            DO k=1,N(ng)
              DO j=Jmin,Jmax
                DO i=Imin,Imax
                  TIDES(ng) % t_tide(i,j,k,itide,itrc) = IniVal
                END DO
              END DO
            END DO
          END DO
        END DO
      END IF
#  endif
# endif
# ifdef UV_WAVEDRAG
      DO j=Jmin,Jmax
        DO i=Imin,Imax
          TIDES(ng) % sum_ubot(i,j) = IniVal
          TIDES(ng) % sum_vbot(i,j) = IniVal
          TIDES(ng) % filt_ubot(i,j) = IniVal
          TIDES(ng) % filt_vbot(i,j) = IniVal
        END DO
      END DO
# endif

      RETURN
      END SUBROUTINE initialize_tides
#endif
      END MODULE mod_tides
