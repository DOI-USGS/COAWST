#include "cppdefs.h"
      MODULE mod_tides
#if defined SSH_TIDES || defined UV_TIDES || defined POT_TIDES
!
!svn $Id: mod_tides.F 1382 2011-02-19 01:44:37Z kate $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
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
# if defined AVERAGES_DETIDE && defined AVERAGES
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

        TYPE T_TIDES

          real(r8), pointer :: Tperiod(:)
# if defined AVERAGES_DETIDE && defined AVERAGES
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
# if defined POT_TIDES
          real(r8), pointer :: POT_Tamp(:,:,:)
          real(r8), pointer :: POT_Tphase(:,:,:)
          real(r8), pointer :: Ptide(:,:)
# endif
# if defined TIDES_ASTRO
          real(r8), pointer :: Vu_sat(:,:,:)
          real(r8), pointer :: f_sat(:,:,:)
# endif
# if defined AVERAGES_DETIDE && defined AVERAGES
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

        TYPE (T_TIDES), allocatable :: TIDES(:)

# if defined TIDES_ASTRO
        integer, parameter :: MAX_SAT = 10
        integer, allocatable  :: doodson(:,:)
        real(r8), allocatable :: phase_corr(:)
        integer, allocatable  :: num_sat(:)
        integer, allocatable  :: sat_doodson(:,:,:)
        real(r8), allocatable :: sat_phase_corr(:,:)
        real(r8), allocatable :: sat_amp(:,:)
        real(r8), allocatable :: sat_flag(:,:)
# endif
# ifdef UV_WAVEDRAG
        integer :: nstep_hour
        integer :: hour_steps = 0
        real(r8), parameter :: drag_scale = 0.375
# endif

      CONTAINS

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
      USE mod_scalars
      USE mod_stepping
!
! Inported variable declarations.
!
      integer, intent(in) :: ng, LBi, UBi, LBj, UBj
!
!  Local variable declarations.
!
      logical :: foundit

      integer :: Nfiles, Vid, i, ifile, mg, nvatt, nvdim
!
!-----------------------------------------------------------------------
!  Allocate module variables.
!-----------------------------------------------------------------------
!
!  Inquire about the maximum number of tidal components.
!
      IF (ng.eq.1) THEN
        MTC=0
        DO mg=1,Ngrids
          IF (.not.(RefinedGrid(ng).and.RefineScale(ng).gt.0)) THEN
            foundit=.FALSE.

            QUERY : DO ifile=1,nFfiles(mg)
              CALL netcdf_inq_var (ng, iNLM, FRC(ifile,mg)%name,        &
     &                             MyVarName = TRIM(Vname(1,idTper)),   &
     &                             SearchVar = foundit,                 &
     &                             VarID = Vid,                         &
     &                             nVardim = nvdim,                     &
     &                             nVarAtt = nvatt)
              IF (exit_flag.ne.NoError) RETURN

# if defined AVERAGES_DETIDE && defined AVERAGES
!
!  Check if detiding variables are available and set definition switch.
!
              IF (foundit) THEN
                DO i=1,n_var
                  IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idFsuH))) THEN
                    LdefTIDE(ng)=.FALSE.
                  END IF
                END DO
              END IF
# endif
!
!  Set maximum number of tidal components.  Allocate and initialize
!  TIDE I/O structure.
!
              IF (foundit) THEN
                MTC=MAX(MTC,var_Dsize(1))            ! first dimension
                NTC(mg)=var_Dsize(1)
!
                Nfiles=FRC(ifile,mg)%Nfiles
                allocate ( TIDE(ng)%Nrec(Nfiles) )
                allocate ( TIDE(ng)%time_min(Nfiles) )
                allocate ( TIDE(ng)%time_max(Nfiles) )
                allocate ( TIDE(ng)%Vid(NV) )
                allocate ( TIDE(ng)%Tid(MT) )
                allocate ( TIDE(ng)%files(Nfiles) )
                TIDE(ng)%Nfiles=Nfiles
                TIDE(ng)%Fcount=1
                TIDE(ng)%Rindex=0
                TIDE(ng)%ncid=-1
                TIDE(ng)%Vid=-1
                TIDE(ng)%Tid=-1
                TIDE(ng)%Nrec(1:Nfiles)=0
                TIDE(ng)%time_min(1:Nfiles)=0.0_r8
                TIDE(ng)%time_max(1:Nfiles)=0.0_r8
                TIDE(ng)%name=TRIM(FRC(ifile,mg)%name)
                TIDE(ng)%label='TIDE - tidal forcing'
                EXIT QUERY
              END IF
            END DO QUERY
          END IF
        END DO
      END IF
!
!  Allocate structure.
!
      IF (ng.eq.1) allocate ( TIDES(Ngrids) )
!
!  Allocate tidal forcing variables.
!
      allocate ( TIDES(ng) % Tperiod(MTC)  )

# if defined SSH_TIDES
      allocate ( TIDES(ng) % SSH_Tamp(LBi:UBi,LBj:UBj,MTC) )
      allocate ( TIDES(ng) % SSH_Tphase(LBi:UBi,LBj:UBj,MTC) )
# endif

# if defined UV_TIDES
      allocate ( TIDES(ng) % UV_Tangle(LBi:UBi,LBj:UBj,MTC) )
      allocate ( TIDES(ng) % UV_Tmajor(LBi:UBi,LBj:UBj,MTC) )
      allocate ( TIDES(ng) % UV_Tminor(LBi:UBi,LBj:UBj,MTC) )
      allocate ( TIDES(ng) % UV_Tphase(LBi:UBi,LBj:UBj,MTC) )
# endif

# if defined POT_TIDES
      allocate ( TIDES(ng) % POT_Tamp(LBi:UBi,LBj:UBj,MTC) )
      allocate ( TIDES(ng) % POT_Tphase(LBi:UBi,LBj:UBj,MTC) )
      allocate ( TIDES(ng) % Ptide(LBi:UBi,LBj:UBj) )
# endif

# ifdef TIDES_ASTRO
      allocate( TIDES(ng) % Vu_sat(LBi:UBi,LBj:UBj,MTC))
      allocate( TIDES(ng) % f_sat(LBi:UBi,LBj:UBj,MTC))
# endif

# if defined TIDES_ASTRO
       allocate ( doodson(6,MTC) )
       allocate ( phase_corr(MTC) )
       allocate ( num_sat(MTC) )
       allocate ( sat_doodson(3,MAX_SAT,MTC) )
       allocate ( sat_phase_corr(MAX_SAT,MTC) )
       allocate ( sat_amp(MAX_SAT,MTC) )
       allocate ( sat_flag(MAX_SAT,MTC) )
# endif

# if defined AVERAGES_DETIDE && defined AVERAGES
!
!  Allocate variables used for the least-squares detiding of
!  time-averaged fields.
!
      allocate ( TIDES(ng) % CosOmega(MTC) )
      allocate ( TIDES(ng) % SinOmega(MTC) )
      allocate ( TIDES(ng) % CosW_avg(MTC) )
      allocate ( TIDES(ng) % CosW_sum(MTC) )
      allocate ( TIDES(ng) % SinW_avg(MTC) )
      allocate ( TIDES(ng) % SinW_sum(MTC) )
      allocate ( TIDES(ng) % CosWCosW(MTC,MTC) )
      allocate ( TIDES(ng) % SinWSinW(MTC,MTC) )
      allocate ( TIDES(ng) % SinWCosW(MTC,MTC) )

      IF (Aout(idFsuD,ng)) THEN
        allocate ( TIDES(ng) % zeta_detided(LBi:UBi,LBj:UBj) )
        allocate ( TIDES(ng) % zeta_tide(LBi:UBi,LBj:UBj,0:2*MTC) )
      END IF

      IF (Aout(idu2dD,ng)) THEN
        allocate ( TIDES(ng) % ubar_detided(LBi:UBi,LBj:UBj) )
        allocate ( TIDES(ng) % ubar_tide(LBi:UBi,LBj:UBj,0:2*MTC) )
      END IF

      IF (Aout(idv2dD,ng)) THEN
        allocate ( TIDES(ng) % vbar_detided(LBi:UBi,LBj:UBj) )
        allocate ( TIDES(ng) % vbar_tide(LBi:UBi,LBj:UBj,0:2*MTC) )
      END IF

#  ifdef SOLVE3D
      IF (Aout(idu3dD,ng)) THEN
        allocate ( TIDES(ng) % u_detided(LBi:UBi,LBj:UBj,N(ng)) )
        allocate ( TIDES(ng) % u_tide(LBi:UBi,LBj:UBj,N(ng),0:2*MTC) )
      END IF

      IF (Aout(idv3dD,ng)) THEN
        allocate ( TIDES(ng) % v_detided(LBi:UBi,LBj:UBj,N(ng)) )
        allocate ( TIDES(ng) % v_tide(LBi:UBi,LBj:UBj,N(ng),0:2*MTC) )
      END IF

      IF (ANY(Aout(idTrcD(:),ng))) THEN
        allocate ( TIDES(ng) % t_detided(LBi:UBi,LBj:UBj,N(ng),NAT) )
        allocate ( TIDES(ng) % t_tide(LBi:UBi,LBj:UBj,N(ng),            &
     &                                0:2*MTC,NAT) )
      END IF
#  endif
# endif
# ifdef UV_WAVEDRAG
      allocate ( TIDES(ng) % sum_ubot(LBi:UBi,LBj:UBj) )
      allocate ( TIDES(ng) % sum_vbot(LBi:UBi,LBj:UBj) )
      allocate ( TIDES(ng) % filt_ubot(LBi:UBi,LBj:UBj) )
      allocate ( TIDES(ng) % filt_vbot(LBi:UBi,LBj:UBj) )
# endif

      RETURN
      END SUBROUTINE allocate_tides

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

# if defined AVERAGES_DETIDE && defined AVERAGES
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
# endif

# if defined POT_TIDES
      DO itide=1,MTC
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            TIDES(ng) % POT_Tamp(i,j,itide) = IniVal
            TIDES(ng) % POT_Tphase(i,j,itide) = IniVal
          END DO
        END DO
      END DO
      DO j=Jmin,Jmax
        DO i=Imin,Imax
          TIDES(ng) % Ptide(i,j) = IniVal
        END DO
      END DO
# endif

# ifdef TIDES_ASTRO
      DO itide=1,MTC
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            TIDES(ng) % Vu_sat(i,j,itide) = 0.0_r8
            TIDES(ng) % f_sat(i,j,itide) = 1.0_r8
          END DO
        END DO
      END DO
# endif

# if defined TIDES_ASTRO
! Assumes tides are in the order Q1, O1, P1, K1, N2, M2, S2, K2

      doodson = reshape( (/ 1, -2, 0, 1, 0, 0,                          &
     &                      1, -1, 0, 0, 0, 0,                          &
     &                      1, 1, -2, 0, 0, 0,                          &
     &                      1, 1,  0, 0, 0, 0,                          &
     &                      2, -1, 0, 1, 0, 0,                          &
     &                      2, 0,  0, 0, 0, 0,                          &
     &                      2, 2, -2, 0, 0, 0,                          &
     &                      2, 2,  0, 0, 0, 0/), (/ 6,MTC /) )

      phase_corr = (/ -0.25, -0.25, -0.25, -0.75,                       &
     &                             0.0, 0.0, 0.0, 0.0 /)

      num_sat = (/ 10, 8, 6, 10, 4, 9, 3, 5 /)

      sat_doodson(:,:,1) = reshape( (/ -2, -3, 0, -2, -2, 0,            &
     &                                 -1, -2, 0, -1, -1, 0,            &
     &                                 -1,  0, 0,  0, -2, 0,            &
     &                                 -1,  0, 1,  0, -1, 0,            &
     &                                  1,  0, 0,  2,  0, 0 /),         &
     &                               (/ 3,MAX_SAT /) )
      sat_doodson(:,:,2) = reshape( (/ -1, 0, 0, 0, -2, 0,              &
     &                                 0, -1, 0, 1, -1, 0,              &
     &                                 1,  0, 0, 1,  1, 0,              &
     &                                 2,  0, 0, 2,  1, 0,              &
     &                                 -99, -99, -99, -99, -99, -99 /), &
     &                                (/ 3,MAX_SAT /) )
      sat_doodson(:,:,3) = reshape( (/ 0, -2, 0, 0, -1, 0,              &
     &                                 0,  0, 2, 1,  0, 0,              &
     &                                 2,  0, 0, 2,  1, 0,              &
     &                                 -99, -99, -99, -99, -99, -99,    &
     &                                 -99, -99, -99, -99, -99, -99 /), &
     &                                (/ 3,MAX_SAT /) )
      sat_doodson(:,:,4) = reshape( (/ -2, -1, 0, -1, -1, 0,            &
     &                                 -1,  0, 0, -1,  1, 0,            &
     &                                 0, -2, 0,  0, -1, 0,             &
     &                                 0,  1, 0,  0,  2, 0,             &
     &                                 1,  0, 0,  1,  1, 0 /),          &
     &                                (/ 3,MAX_SAT /) )
      sat_doodson(:,:,5) = reshape( (/ -2, -2, 0, -1,  0, 1,            &
     &                                  0, -2, 0,  0, -1, 0,            &
     &                                 -99, -99, -99, -99, -99, -99,    &
     &                                 -99, -99, -99, -99, -99, -99,    &
     &                                 -99, -99, -99, -99, -99, -99 /), &
     &                                (/ 3,MAX_SAT /) )
      sat_doodson(:,:,6) = reshape( (/ -1, -1, 0, -1,  0, 0,            &
     &                                  0, -2, 0,  0, -1, 0,            &
     &                                  1, -1, 0,  1,  0, 0,            &
     &                                  1,  1, 0,  2,  0, 0,            &
     &                                  2,  1, 0, -99, -99, -99 /),     &
     &                                (/ 3,MAX_SAT /) )
      sat_doodson(:,:,7) = reshape( (/ 0, -1, 0, 1, 0, 0,               &
     &                                 2,  0, 0, -99, -99, -99,         &
     &                                 -99, -99, -99, -99, -99, -99,    &
     &                                 -99, -99, -99, -99, -99, -99,    &
     &                                 -99, -99, -99, -99, -99, -99 /), &
     &                                (/ 3,MAX_SAT /) )
      sat_doodson(:,:,8) = reshape( (/ -1, 0, 0, -1, 1, 0,              &
     &                                 0, -1, 0,  0, 1, 0,              &
     &                                 0,  2, 0, -99, -99, -99,         &
     &                                 -99, -99, -99, -99, -99, -99,    &
     &                                 -99, -99, -99, -99, -99, -99 /), &
     &                                (/ 3,MAX_SAT /) )

      sat_phase_corr = reshape(                                         &
     & (/ 0.5, 0.5, 0.75, 0.75, 0.75, 0.5, 0.0, 0.0, 0.75, 0.5,         &
     &    0.25, 0.5, 0.0, 0.25, 0.75, 0.25, 0.5, 0.5, -99., -99.,       &
     &    0.0, 0.5, 0.5, 0.75, 0.5, 0.5, -99., -99., -99., -99.,        &
     &    0.0, 0.75, 0.25, 0.75, 0.0, 0.5, 0.0, 0.5, 0.25, 0.25,        &
     &    0.5, 0.0, 0.0, 0.5, -99., -99., -99., -99., -99., -99.,       &
     &    0.75, 0.75, 0.0, 0.5, 0.25, 0.75, 0.75, 0.0, 0.0, -99.,       &
     &    0.0, 0.75, 0.0, -99., -99., -99., -99., -99., -99., -99.,     &
     &    0.75, 0.75, 0.5, 0.0, 0.0, -99., -99., -99., -99., -99. /),   &
     &      (/ MAX_SAT,MTC /) )

      sat_amp = reshape(                                                &
     &      (/ 0.0007, 0.0038, 0.0010, 0.0115, 0.0292,                  &
     &         0.0057, 0.0008, 0.1884, 0.0018, 0.0028,                  &
     &         0.0003, 0.0058, 0.1885, 0.0004, 0.0029,                  &
     &         0.0004, 0.0064, 0.0010, 9999.0, 9999.0,                  &
     &         0.0008, 0.0112, 0.0004, 0.0004, 0.0015,                  &
     &         0.0003, 9999.0, 9999.0, 9999.0, 9999.0,                  &
     &         0.0002, 0.0001, 0.0007, 0.0001, 0.0001,                  &
     &         0.0198, 0.1356, 0.0029, 0.0002, 0.0001,                  &
     &         0.0039, 0.0008, 0.0005, 0.0373, 9999.0,                  &
     &         9999.0, 9999.0, 9999.0, 9999.0, 9999.0,                  &
     &         0.0001, 0.0004, 0.0005, 0.0373, 0.0001,                  &
     &         0.0009, 0.0002, 0.0006, 0.0002, 9999.0,                  &
     &         0.0022, 0.0001, 0.0001, 9999.0, 9999.0,                  &
     &         9999.0, 9999.0, 9999.0, 9999.0, 9999.0,                  &
     &         0.0024, 0.0004, 0.0128, 0.2980, 0.0324,                  &
     &         9999.0, 9999.0, 9999.0, 9999.0, 9999.0 /),               &
     &      (/ MAX_SAT,MTC /) )

      sat_flag = reshape(                                               &
     &       (/ 0, 0, 1, 1, 1, 0, 0, 0, 1, 0,                           &
     &          1, 0, 0, 1, 1, 1, 0, 0, -99, -99,                       &
     &          0, 0, 0, 1, 0, 0, -99, -99, -99, -99,                   &
     &          0, 1, 1, 1, 0, 0, 0, 0, 1, 1,                           &
     &          0, 0, 0, 0, -99, -99, -99, -99, -99, -99,               &
     &          2, 2, 0, 0, 2, 2, 2, 0, 0, -99,                         &
     &          0, 2, 0, -99, -99, -99, -99, -99, -99, -99,             &
     &          2, 2, 0, 0, 0, -99, -99, -99, -99, -99 /),              &
     &      (/ MAX_SAT,MTC /) )

# endif
# if defined AVERAGES_DETIDE && defined AVERAGES
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

#ifdef TIDES_ASTRO
      SUBROUTINE tide_astro(ttime,vux,fx,xlat,ng,tile,                  &
     &                      LBi, UBi, LBj, UBj)
      USE mod_scalars
      implicit none

      integer, intent(in)  :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      real(r8), intent(in)  :: ttime
# ifdef ASSUMED_SHAPE
      real(r8), intent(in)  :: xlat(LBi:,LBj:)
      real(r8), intent(out) :: vux(LBi:,LBj:,:)
      real(r8), intent(out) :: fx(LBi:,LBj:,:)
# else
      real(r8), intent(in)  :: xlat(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: vux(LBi:UBi,LBj:UBj,MTC)
      real(r8), intent(out) :: fx(LBi:UBi,LBj:UBj,MTC)
# endif

! Local variables
      real(r8) :: slat, hh, d1, tau, dtau_2, freq(MTC)
      real(r8) :: h, pp, s, p, enp, dh, dpp, ds, dp, dnp
      real(r8) :: uu, twopi, sumc, sums, v, rr, vdbl, uudbl
      integer  :: i, iv, j, k, iuu, itide, intdays
      integer  :: Imin, Imax, Jmin, Jmax
      real(r8), parameter :: eps = 1.e-20_r8
!
!  Set array initialization range.
!
#include "set_bounds.h"

#ifdef _OPENMP
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
#else
      Imin=BOUNDS(ng)%LBi(tile)
      Imax=BOUNDS(ng)%UBi(tile)
      Jmin=BOUNDS(ng)%LBj(tile)
      Jmax=BOUNDS(ng)%UBj(tile)
#endif
!
!***********************************************************************
!*  THIS SUBROUTINE CALCULATES V (ASTRONOMICAL PHASE ARGUMENT), U AND F
!*  (NODAL MODULATION PHASE AND AMPLITUDE CORRECTIONS) FOR ALL CONSTITU-
!*  ENTS.
!***********************************************************************
!*  NTIDAL IS THE NUMBER OF MAIN CONSTITUENTS
!*  NTOTAL IS THE NUMBER OF CONSTITUENTS (MAIN + SHALLOW WATER)
!*  FOR  THE GIVEN TIME KH, THE TABLE OF F AND V+U VALUES IS
!*  CALCULATED FOR ALL THE CONSTITUENTS.
!*     F IS THE NODAL MODULATION ADJUSTMENT FACTOR FOR AMPLITUDE
!*     U IS THE NODAL MODULATION ADJUSTMENT FACTOR FOR PHASE
!*     V IS THE ASTRONOMICAL ARGUMENT ADJUSTMENT FOR PHASE.
!
!
!***********************************************************************
!  The astromical arguments are computed at ttime
!
!     d1 is days measured from 1200 UT December 31, 1899
!
      d1=ttime*sec2day
      twopi = 2.0_r8*pi
      call astr(d1,h,pp,s,p,enp,dh,dpp,ds,dp,dnp)
      intdays=int(ttime*sec2day)
      hh=real(ttime-intdays*day2sec,r8)/3600._r8
      tau = hh/24._r8 + h - s
      dtau_2 = 365.0_r8 + dh - ds
!
!***********************************************************************
!*  ONLY THE FRACTIONAL PART OF A SOLAR DAY NEED BE RETAINED FOR COMPU-
!*  TING THE LUNAR TIME TAU.
!
      DO itide=1,MTC
        freq(itide)=(doodson(1,itide)*dtau_2+doodson(2,itide)*ds+       &
     &               doodson(3,itide)*dh+doodson(4,itide)*dp+           &
     &               doodson(5,itide)*dnp+doodson(6,itide)*dpp)         &
     &               /(365._r8*24._r8)
        vdbl=doodson(1,itide)*tau+doodson(2,itide)*s+                   &
     &               doodson(3,itide)*h+doodson(4,itide)*p+             &
     &               doodson(5,itide)*enp+doodson(6,itide)*pp+          &
     &               phase_corr(itide)
        iv=vdbl
        iv=(iv/2)*2
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            v=vdbl-iv
            sumc=1.
            sums=0.
            slat=SIN(deg2rad*xlat(i,j))
            DO k=1,num_sat(itide)
!
!***********************************************************************
!*  HERE THE SATELLITE AMPLITUDE RATIO ADJUSTMENT FOR LATITUDE IS MADE
!
              rr=sat_amp(k,itide)
              IF (sat_flag(k,itide) == 1) THEN
                rr=sat_amp(k,itide)*0.36309*(1.-5.*slat*slat)/(eps+slat)
              ELSE IF (sat_flag(k,itide) == 2) THEN
                rr=sat_amp(k,itide)*2.59808*slat
              END IF
              uudbl=sat_doodson(1,k,itide)*p+sat_doodson(2,k,itide)*enp+&
     &              sat_doodson(3,k,itide)*pp+sat_phase_corr(k,itide)
              iuu=uudbl
              uu=uudbl-iuu
              sumc=sumc+rr*COS(uu*twopi)
              sums=sums+rr*SIN(uu*twopi)
            END DO
            fx(i,j,itide)=SQRT(sumc*sumc+sums*sums)
            vux(i,j,itide)=v+ATAN2(sums,sumc)/twopi
            v=v-int(v)
            if (v.lt.0.) v=v+1.
            vux(i,j,itide)=vux(i,j,itide)-int(vux(i,j,itide))
            if (vux(i,j,itide).lt.0.) vux(i,j,itide)=vux(i,j,itide)+1.
!       write(6,998) kon(k),v*360.,f(k),360.*(vux(k)-v)
!       write(7,998) kon(k),v*360.,f(k),360.*(vux(k)-v)
!998    format(' ',a5,5x,3f12.5)
          END DO
        END DO
      END DO
!
      RETURN
      END SUBROUTINE tide_astro
      SUBROUTINE astr(d1,h,pp,s,p,np,dh,dpp,ds,dp,dnp)
!        this subroutine calculates the following five ephermides
!        of the sun and moon
!        h = mean longitude of the sum
!        pp = mean longitude of the solar perigee
!        s = mean longitude of the moon
!        p = mean longitude of the lunar perigee
!        np = negative of the longitude of the mean ascending node
!        and their rates of change.
!        Units for the ephermides are cycles and for their derivatives
!        are cycles/365 days
!        The formulae for calculating this ephermides were taken from
!        pages 98 and 107 of the Explanatory Supplement to the
!        Astronomical Ephermeris and the American Ephermis and
!        Nautical Almanac (1961)
!
        implicit none
      real(r8), intent(in)  :: d1
      real(r8), intent(out) :: h, pp, s, p, np, dh, dpp, ds, dp, dnp
      real(r8)  :: d2, f, f2

        d2=d1*1.d-4
        f=360.d0
        f2=f/365.d0
        h=279.696678d0+.9856473354d0*d1+.00002267d0*d2*d2
        pp=281.220833d0+.0000470684d0*d1+.0000339d0*d2*d2+              &
     &  .00000007d0*d2**3
        s=270.434164d0+13.1763965268d0*d1-.000085d0*d2*d2+              &
     &  .000000039d0*d2**3
        p=334.329556d0+.1114040803d0*d1-.0007739d0*d2*d2-               &
     &  .00000026d0*d2**3
        np=-259.183275d0+.0529539222d0*d1-.0001557d0*d2*d2-             &
     &  .00000005d0*d2**3
        h=h/f
        pp=pp/f
        s=s/f
        p=p/f
        np=np/f
        h=h-int(h)
        pp=pp-int(pp)
        s=s-int(s)
        p=p-int(p)
        np=np-int(np)
        dh=.9856473354d0+2.d-8*.00002267d0*d1
        dpp=.0000470684d0+2.d-8*.0000339d0*d1                           &
     &  +3.d-12*.00000007d0*d1**2
        ds=13.1763965268d0-2.d-8*.000085d0*d1+                          &
     &  3.d-12*.000000039d0*d1**2
        dp=.1114040803d0-2.d-8*.0007739d0*d1-                           &
     &  3.d-12*.00000026d0*d1**2
        dnp=+.0529539222d0-2.d-8*.0001557d0*d1-                         &
     &  3.d-12*.00000005d0*d1**2
        dh=dh/f2
        dpp=dpp/f2
        ds=ds/f2
        dp=dp/f2
        dnp=dnp/f2
        RETURN
      END SUBROUTINE astr
#endif
      END MODULE mod_tides
