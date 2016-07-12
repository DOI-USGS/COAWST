#include "cppdefs.h"
      MODULE mod_scalars
!
!svn $Id: mod_scalars.F 807 2016-07-09 02:03:55Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!
        USE mod_param
#ifdef NO_4BYTE_REALS
        USE netcdf
#endif
!
        implicit none
!
!-----------------------------------------------------------------------
! Multiple grid structure.
!-----------------------------------------------------------------------
!
#if defined READ_WATER && defined MASKING && defined DISTRIBUTE
!    IJwater       IJ-indices of water points.
#endif
!    Fstate        Logical switches to control computations of the
!                    Forcing Singular Vectors or Stochastic Optimals.
!    Lstate        Logical switches to control computations of the
!                    model state.
#ifdef STATIONS
!    Sflag         Station extraction special flag:
!                    Sflag = 0  => locations in terms of (I,J) pairs.
!                    Sflag = 1  => locations in terms of (lon,lat) pairs.
!    SposX         Longitude or frational I-coordinate station location.
!    SposY         Latitude  or frational J-coordinate station location.
#endif
!    Cs_r          Set of S-curves used to stretch the vertical grid
!                    that follows the bathymetry at vertical RHO-points.
!    Cs_w          Set of S-curves used to stretch the vertical grid
!                    that follows the bathymetry at vertical W-points.
!    sc_r          S-coordinate independent variable, [-1 < sc < 0] at
!                    vertical RHO-points.
!    sc_w          S-coordinate independent variable, [-1 < sc < 0] at
!                    vertical W-points.
!
        TYPE T_SCALARS

          logical, pointer :: Fstate(:)
          logical, pointer :: Lstate(:)
#if defined READ_WATER && defined MASKING && defined DISTRIBUTE
          integer , pointer :: IJwater(:,:)
#endif
#ifdef STATIONS
          integer,  pointer :: Sflag(:)

          real(r8), pointer :: SposX(:)
          real(r8), pointer :: SposY(:)
#endif
          real(r8), pointer :: Cs_r(:)
          real(r8), pointer :: Cs_w(:)
          real(r8), pointer :: sc_r(:)
          real(r8), pointer :: sc_w(:)

        END TYPE T_SCALARS

        TYPE (T_SCALARS), allocatable :: SCALARS(:)
!
!-----------------------------------------------------------------------
!  Tracer identification indices.
!-----------------------------------------------------------------------
!
        integer :: itemp              ! Potential temperature
        integer :: isalt              ! Salinity
#ifdef T_PASSIVE
        integer, pointer :: inert(:)  ! inert tracers
#endif
#ifdef DIAGNOSTICS
!
!-----------------------------------------------------------------------
!  Diagnostic fields identification indices.
!-----------------------------------------------------------------------
!
# ifdef DIAGNOSTICS_TS
        integer :: iTrate          ! Tracer, time rate of change
        integer :: iTvadv          ! Tracer, vertical advection
        integer :: iThadv          ! Tracer, horizontal advection
        integer :: iTxadv          ! Tracer, horizontal X-advection
        integer :: iTyadv          ! Tracer, horizontal Y-advection
        integer :: iTvdif          ! Tracer, vertical diffusion
        integer :: iThdif          ! Tracer, horizontal diffusion
        integer :: iTxdif          ! Tracer, horizontal X-diffusion
        integer :: iTydif          ! Tracer, horizontal Y-diffusion
        integer :: iTsdif          ! Tracer, horizontal S-diffusion
# endif
# ifdef DIAGNOSTICS_UV
        integer :: M2fcor          ! 2D momentum, Coriolis
        integer :: M2hadv          ! 2D momentum, horizontal advection
        integer :: M2xadv          ! 2D momentum, horizontal X-advection
        integer :: M2yadv          ! 2D momentum, horizontal Y-advection
#  ifdef WEC_MELLOR
        integer :: M2hrad          ! 2D momentum, horizontal stresses
#  endif
#  ifdef WEC_VF
        integer :: M2hjvf          ! 2D momentum, horizontal J vortex force
        integer :: M2kvrf          ! 2D momentum, K vortex force
        integer :: M2fsco          ! 2D momentum, coriolis-stokes
        integer :: M2bstm          ! 2D momentum, bottom streaming
        integer :: M2sstm          ! 2D momentum, surface streaming
        integer :: M2wrol          ! 2D momentum, wave roller accel
        integer :: M2wbrk          ! 2D momentum, wave breaking
        integer :: M2zeta          ! 2D momentum, Eulerian sea level adjustment
        integer :: M2zetw          ! 2D momentum, quasi-static sea level adjustment
        integer :: M2zqsp          ! 2D momentum, quasi-static pressure
        integer :: M2zbeh          ! 2D momentum, Bernoulli head
#  endif
        integer :: M2pgrd          ! 2D momentum, pressure gradient
        integer :: M2hvis          ! 2D momentum, horizontal viscosity
        integer :: M2xvis          ! 2D momentum, horizontal X-viscosity
        integer :: M2yvis          ! 2D momentum, horizontal Y-viscosity
        integer :: M2sstr          ! 2D momentum, surface stress
        integer :: M2bstr          ! 2D momentum, bottom stress
        integer :: M2rate          ! 2D momentum, time rate of change
#  ifdef SOLVE3D
        integer :: M3fcor          ! 3D momentum, Coriolis
        integer :: M3vadv          ! 3D momentum, vertical advection
        integer :: M3hadv          ! 3D momentum, horizontal advection
        integer :: M3xadv          ! 3D momentum, horizontal X-advection
        integer :: M3yadv          ! 3D momentum, horizontal Y-advection
#   ifdef WEC_MELLOR
        integer :: M3hrad          ! 3D momentum, horizontal stresses
        integer :: M3vrad          ! 3D momentum, vertical stresses
#   endif
#   ifdef WEC_VF
        integer :: M3vjvf          ! 3D momentum, vertical J vortex force
        integer :: M3hjvf          ! 3D momentum, horizontal J vortex force
        integer :: M3kvrf          ! 3D momentum, K vortex force
        integer :: M3fsco          ! 3D momentum, coriolis-stokes
        integer :: M3bstm          ! 3D momentum, bottom streaming
        integer :: M3sstm          ! 3D momentum, surface streaming
        integer :: M3wrol          ! 3D momentum, wave roller accel
        integer :: M3wbrk          ! 3D momentum, wave breaking
#   endif
#   if defined VEGETATION && defined VEG_DRAG
        integer :: M3fveg          ! 3D momentum, vegetation drag force
        integer :: M2fveg          ! 2D momentum, vegetation drag force
#   endif
        integer :: M3pgrd          ! 3D momentum, pressure gradient
        integer :: M3vvis          ! 3D momentum, vertical viscosity
        integer :: M3hvis          ! 3D momentum, horizontal viscosity
        integer :: M3xvis          ! 3D momentum, horizontal X-viscosity
        integer :: M3yvis          ! 3D momentum, horizontal Y-viscosity
        integer :: M3rate          ! 3D momentum, time rate of change
#  endif
# endif
#endif

!
!-----------------------------------------------------------------------
!  Time stepping indices, variables, and clocks.
!-----------------------------------------------------------------------
!
!    indx1         2D timestep rolling counter.
!    iic           Timestep counter for 3D primitive equations.
!    iif           Timestep counter for 2D primitive equations.
!    ndtfast       Number of barotropic timesteps between each
!                    baroclinic timestep.
!    nfast         Number of barotropic timesteps needed to compute
!                    time-averaged barotropic variables centered at
!                    time level n+1.
!    dt            Size baroclinic timestep (s).
!    dtfast        Size barotropic timestep (s).
!    dtau          Size of age increment
!    run_time      Total run time for all nested grids (s).
!    io_time       Current I/O time (s) processed in "get_state".
!    tdays         Model time clock (days).
!    time          Model time clock (s).
!    time_code     Model time clock (string, Day HH:MM:SS)
!    AVGtime       Model time clock for averages output (s).
!    AVG2time      Model time clock for averages output (s).
!    DIAtime       Model time clock for diagnostics output (s).
!    IMPtime       Impulse forcing time (s) to process.
!    ObsTime       Observation time (s) to process.
!    FrcTime       Adjoint or tangent linear Impulse forcing time (s).
!    r_date        Model initialization reference date (vector):
!                    r_date(1) => reference date (YYYYMMDD.dd).
!                    r_date(2) => year.
!                    r_date(3) => year day.
!                    r_date(4) => month.
!                    r_date(5) => day.
!                    r_date(6) => hour.
!                    r_date(7) => minute.
!                    r_date(8) => second.
!    dstart        Time stamp assigned to model initialization (usually
!                    a Calendar day, like modified Julian Day).
!    tide_start    Reference time for tidal forcing (days).
!    time_ref      Reference time for "units" attribute (YYYYMMDD.dd).
!    r_text        Model initialization reference text (string).
#ifdef CICE_MODEL
!    tspy          timesteps per year
!    tspd          timesteps per day
!    dleftinSep    days in September before Oct1 restart
!    Jan,Feb,...,Nov,Dec  days after restart until first of each month
#endif
!
        logical, allocatable :: PerfectRST(:)
        logical, allocatable :: PREDICTOR_2D_STEP(:)
!$OMP THREADPRIVATE (PREDICTOR_2D_STEP)

        integer, allocatable :: indx1(:)
        integer, allocatable :: iic(:)
        integer, allocatable :: iif(:)
!$OMP THREADPRIVATE (indx1, iic, iif)

        integer, allocatable :: ndtfast(:)
        integer, allocatable :: nfast(:)

        real(r8), allocatable :: tdays(:)                ! days
        real(r8), allocatable :: time(:)                 ! seconds
!$OMP THREADPRIVATE (tdays, time)

        real(r8), allocatable :: dt(:)                   ! seconds
        real(r8), allocatable :: dtfast(:)               ! seconds

        real(r8), allocatable :: TimeEnd(:)              ! seconds
        real(r8), allocatable :: AVGtime(:)              ! seconds
        real(r8), allocatable :: AVG2time(:)             ! seconds
        real(r8), allocatable :: DIAtime(:)              ! seconds
        real(r8), allocatable :: IMPtime(:)              ! seconds
        real(r8), allocatable :: ObsTime(:)              ! seconds
        real(r8), allocatable :: FrcTime(:)              ! seconds
#ifdef TIDES_ASTRO
        logical, allocatable :: FIRST_TIDES_ASTRO(:)
#endif
#ifdef AGE_DISTRIBUTION
        real(r8), allocatable :: dtau(:)              ! seconds
#endif

        real(r8), dimension(8) :: r_date
        real(r8) :: dstart = 0.0_r8                      ! days
        real(r8) :: io_time = 0.0_r8                     ! seconds
        real(r8) :: run_time = 0.0_r8                    ! seconds
        real(r8) :: tide_start = 0.0_r8                  ! days
        real(r8) :: time_ref = 0.0_r8                    ! YYYYMMDD.dd
        real(r8) :: days_year = 365.25_r8                ! days per year
#ifdef CICE_MODEL
        real(r8), allocatable :: tspy(:)                 ! timesteps
        real(r8), allocatable :: tspd(:)                 ! timesteps
        real(r8) :: dleftinSep = 0.0_r8                  ! days
        real(r8) :: Jan = 0.0_r8                         ! days
        real(r8) :: Feb = 0.0_r8                         ! days
        real(r8) :: Mar = 0.0_r8                         ! days
        real(r8) :: Apr = 0.0_r8                         ! days
        real(r8) :: May = 0.0_r8                         ! days
        real(r8) :: Jun = 0.0_r8                         ! days
        real(r8) :: Jul = 0.0_r8                         ! days
        real(r8) :: Aug = 0.0_r8                         ! days
        real(r8) :: Sep = 0.0_r8                         ! days
        real(r8) :: Oct = 0.0_r8                         ! days
        real(r8) :: Nov = 0.0_r8                         ! days
        real(r8) :: Dec = 0.0_r8                         ! days
        real(r8), allocatable :: rhoice(:)
        real(r8), allocatable :: min_a(:)
#endif

        character (len=14), allocatable :: time_code(:)  ! DD HH:MM:SS
!$OMP THREADPRIVATE (time_code)

        character (len=19) :: r_text

#if defined POWER_LAW && defined SOLVE3D
!
!  Power-law shape filter parameters for time-averaging of barotropic
!  Fields.  The power-law shape filters are given by:
!
!     F(xi)=xi^Falpha*(1-xi^Fbeta)-Fgamma*xi
!
!  Possible settings of parameters to yield the second-order accuracy:
!
!     Falpha  Fbeta      Fgamma
!     ------------------------------
!      2.0     1.0    0.1181  0.169     The problem here is setting
!      2.0     2.0    0.1576  0.234     Fgamma. Its value here is
!      2.0     3.0    0.1772  0.266     understood as the MAXIMUM
!      2.0     4.0    0.1892  0.284     allowed. It is computed using
!      2.0     5.0    0.1976  0.296     a Newton iteration scheme.
!      2.0     6.0    0.2039  0.304
!      2.0     8.0    0.2129  0.314
!
!  NOTE: Theoretical values of Fgamma presented in the table above are
!  derived assuming "exact" barotropic mode stepping. Consequently, it
!  does not account for effects caused by Forward-Euler (FE) startup
!  of the barotropic mode at every 3D time step.  As the result, the
!  code may become unstable if the theoretical value of Fgamma is used
!  when mode splitting ratio "ndtfast" is small, thus yielding non-
!  negligible start up effects.  To compensate this, the accepted
!  value of Fgamma is reduced relatively to theoretical one, depending
!  on splitting ratio "ndtfast".  This measure is empirical. It is
!  shown to work with setting of "ndtfast" as low as 15, which is
!  more robust that the Hamming Window the squared cosine weights
!  options in "set_weights".
!
        real(r8) :: Falpha = 2.0_r8
        real(r8) :: Fbeta  = 4.0_r8
        real(r8) :: Fgamma = 0.284_r8
#endif
!
!  Total number timesteps in current run. In 3D configurations, "ntimes"
!  is the total of baroclinic timesteps. In 2D configuration, "ntimes"
!  is the total of barotropic timesteps.
!
        integer, allocatable :: ntimes(:)
!
!  Time-step counter for current execution time-window.
!
        integer, allocatable :: step_counter(:)

!$OMP THREADPRIVATE (step_counter)
!
!  Number of time interval divisions for Stochastic Optimals.  It must
!  a multiple of "ntimes".
!
        integer :: Nintervals = 1
!
!  Starting, current, and ending ensemble run parameters.
!
        integer :: ERstr = 1                    ! Starting value
        integer :: ERend = 1                    ! Ending value
        integer :: Ninner = 1                   ! number of inner loops
        integer :: Nouter = 1                   ! number of outer loops
        integer :: Nrun = 1                     ! Current counter
        integer :: inner = 0                    ! inner loop counter
        integer :: outer = 0                    ! outer loop counter
#ifdef SENSITIVITY_4DVAR
        integer :: NrunSAVE = 0                 ! Loop counter
#endif
!
!  First, starting, and ending timestepping parameters
!
        integer, allocatable :: ntfirst(:)      ! Forward-Euler step
        integer, allocatable :: ntstart(:)      ! Start step
        integer, allocatable :: ntend(:)        ! End step
!!$OMP THREADPRIVATE (ntfirst, ntstart, ntend)
!
!  Adjoint model or tangent linear model impulse forcing time record
!  counter and number of records available.
!
        integer, allocatable :: FrcRec(:)
!$OMP THREADPRIVATE (FrcRec)

        integer, allocatable :: NrecFrc(:)
!
!-----------------------------------------------------------------------
!  Control switches.
!-----------------------------------------------------------------------
!
!  Switch to proccess nudging coefficients for radiation open boundary
!  conditions.
!
        logical, allocatable :: NudgingCoeff(:)
!
!  Switch to proccess input boundary data.
!
        logical, allocatable :: ObcData(:)
!
!  These switches are designed to control computational options within
!  nested and/or multiple connected grids.  They are .TRUE. by default.
!  They can turned off for a particular grind in input scripts.
!
        logical, allocatable :: Lbiology(:)
        logical, allocatable :: Lfloats(:)
#ifdef ICE_MODEL
        logical, allocatable :: Lice(:)
#endif
        logical, allocatable :: Lsediment(:)
        logical, allocatable :: Lstations(:)
!
!-----------------------------------------------------------------------
!  Physical constants.
!-----------------------------------------------------------------------
!
!    Cp            Specific heat for seawater (Joules/Kg/degC).
!    Csolar        Solar irradiantion constant (W/m2).
!    Eradius       Earth equatorial radius (m).
!    Infinity      Value resulting when dividing by zero.
!    StefBo        Stefan-Boltzmann constant (W/m2/K4).
!    emmiss        Infrared emmissivity.
!    g             Acceleration due to gravity (m/s2).
!    gorho0        gravity divided by mean density anomaly.
!    rhow          fresh water density (kg/m3).
!    vonKar        von Karman constant.
!
#ifdef ISOMIP
        real(r8) :: Cp = 3974.0_r8              ! mod (4/24/06) for ICETEST
#else
        real(r8) :: Cp = 3985.0_r8              ! Joules/kg/degC
#endif
        real(r8) :: Csolar = 1353.0_r8          ! 1360-1380 W/m2
        real(r8) :: Infinity                    ! Infinity = 1.0/0.0
        real(r8) :: Eradius = 6371315.0_r8      ! m
#ifdef ICE_BOX
        real(r8) :: StefBo = 5.78E-8_r8         ! Watts/m2/K4 (Match MU)
#else
        real(r8) :: StefBo = 5.67E-8_r8         ! Watts/m2/K4
#endif
        real(r8) :: emmiss = 0.97_r8            ! non_dimensional
        real(r8) :: rhow = 1000.0_r8            ! kg/m3
#ifdef SOLITON
        real(r8) :: g = 1.0_r8                  ! non-dimensional
# elif defined WBC_1 || defined WBC_2 || defined WBC_3
        real(r8) :: g = 9.8_r8                  ! m/s2
# elif defined CIRCLE
        real(r8) :: g = 3.92e-2_r8                  ! m/s2
#else
        real(r8) :: g = 9.81_r8                 ! m/s2
#endif
        real(r8) :: gorho0                      ! m4/s2/kg
        real(r8) :: vonKar = 0.41_r8            ! non-dimensional
!
!-----------------------------------------------------------------------
!  Various model parameters.  Some of these parameters are overwritten
!  with the values provided from model standard input script.
!-----------------------------------------------------------------------
!
!  Switch for spherical grid (lon,lat) configurations.
!
        logical :: spherical = .FALSE.
!
!  Switch to compute the grid stiffness.
!
        logical :: Lstiffness = .TRUE.
!$OMP THREADPRIVATE (Lstiffness)
!
!  Composite grid a refined grids switches. They are .FALSE. by default.
!
        logical, allocatable :: CompositeGrid(:,:)
        logical, allocatable :: RefinedGrid(:)
!
!  Refinement grid scale factor from donor grid.
!
        integer, allocatable :: RefineScale(:)
!
!  Switch to extract donor grid (coarse) data at the refinement grid
!  contact point locations. The coarse data is extracted at the first
!  sub-refined time step.  Recall that the finer grid time-step is
!  smaller than the coarser grid by a factor of RefineScale(:). This
!  switch is only relevant during refinement nesting.
!
        logical, allocatable :: GetDonorData(:)
!
!  Periodic boundary swiches for distributed-memory exchanges.
!
        logical, allocatable :: EWperiodic(:)
        logical, allocatable :: NSperiodic(:)
!
!  Lateral open boundary edges volume conservation switches.
!
        logical, allocatable :: VolCons(:,:)
#if defined ADJOINT || defined TANGENT || defined TL_IOMS
        logical, allocatable :: ad_VolCons(:,:)
        logical, allocatable :: tl_VolCons(:,:)
#endif
!
!  Switches to read and process climatology fields.
!
        logical, allocatable :: CLM_FILE(:)          ! Process NetCDF
        logical, allocatable :: Lclimatology(:)      ! any field
        logical, allocatable :: LsshCLM(:)           ! free-surface
        logical, allocatable :: Lm2CLM(:)            ! 2D momentum
        logical, allocatable :: Lm3CLM(:)            ! 3D momentum
        logical, allocatable :: LtracerCLM(:,:)      ! tracers
!
!  Switched to nudge to climatology fields.
!
        logical, allocatable :: Lnudging(:)          ! any field
        logical, allocatable :: LnudgeM2CLM(:)       ! 2D momentum
        logical, allocatable :: LnudgeM3CLM(:)       ! 3D momentum
        logical, allocatable :: LnudgeTCLM(:,:)      ! tracers
!
!  Switches to activate point Source/Sinks in an application:
!    * Horizontal momentum transport (u or v)
!    * Vertical mass transport (w)
!    * Tracer transport
!
        logical, allocatable :: LuvSrc(:)            ! momentum
        logical, allocatable :: LwSrc(:)             ! mass
        logical, allocatable :: LtracerSrc(:,:)      ! tracers
!
!  Execution termination flag.
!
!    exit_flag = 0   No error
!    exit_flag = 1   Blows up
!    exit_flag = 2   Input error
!    exit_flag = 3   Output error
!    exit_flag = 4   IO error
!    exit_flag = 5   Configuration error
!    exit_flag = 6   Partition error
!    exit_flag = 7   Illegal input parameter
!    exit_flag = 8   Fatal algorithm result
!    exit_flag = 9   Frazil ice error
!
        integer :: exit_flag = 0
        integer :: blowup = 0
        integer :: NoError = 0
!
!  Set threshold maximum speed (m/s) and density anomaly (kg/m3) to
!  test if the model is blowing-up.
!
        real(r8), allocatable :: maxspeed(:)
        real(r8), allocatable :: maxrho(:)
!
        real(r8) :: max_speed = 20.0_r8         ! m/s
        real(r8) :: max_rho = 200.0_r8          ! kg/m3
#ifdef BIOLOGY
        real(r8), allocatable :: maxbio(:,:)
        real(r8), allocatable :: max_bio(:)
#endif
!
!  Interpolation scheme.
!
        integer, parameter :: linear = 0        ! linear interpolation
        integer, parameter :: cubic  = 1        ! cubic  interpolation
!
        integer :: InterpFlag = linear          ! interpolation flag
!
!  Shallowest and Deepest levels to apply bottom momentum stresses as
!  a bodyforce
!
        integer, allocatable :: levsfrc(:)
        integer, allocatable :: levbfrc(:)
!
!  Vertical coordinates transform.  Currently, there are two vertical
!  transformation equations (see set_scoord.F for details):
!
!    Original transform (Vtransform=1):
!
!         z_r(x,y,s,t) = Zo_r + zeta(x,y,t) * [1.0 + Zo_r / h(x,y)]
!
!                 Zo_r = hc * [s(k) - C(k)] + C(k) * h(x,y)
!
!    New transform (Vtransform=2):
!
!         z_r(x,y,s,t) = zeta(x,y,t) + [zeta(x,y,t)+ h(x,y)] * Zo_r
!
!                 Zo_r = [hc * s(k) + C(k) * h(x,y)] / [hc + h(x,y)]
!
        integer, allocatable :: Vtransform(:)
!
!  Vertical grid stretching function flag:
!
!    Vstretcing = 1   Original function (Song and Haidvogel, 1994)
!               = 2   A. Shchepetkin (ROMS-UCLA) function
!               = 3   R. Geyer BBL function
!
        integer, allocatable :: Vstretching(:)
!
!  Vertical grid stretching parameters.
!
!    Tcline        Width (m) of surface or bottom boundary layer in
!                    which higher vertical resolution is required
!                    during stretching.
!    hc            S-coordinate critical depth, hc=MIN(hmin,Tcline).
!    theta_s       S-coordinate surface control parameter.
!    theta_b       S-coordinate bottom control parameter.
!
        real(r8), allocatable :: Tcline(:)      ! m, positive
        real(r8), allocatable :: hc(:)          ! m, positive
        real(r8), allocatable :: theta_s(:)     ! 0 < theta_s < 20
        real(r8), allocatable :: theta_b(:)     ! 0 < theta_b < 1
!
!  Bathymetry range values.
!
        real(r8), allocatable :: hmin(:)        ! m, positive
        real(r8), allocatable :: hmax(:)        ! m, positive
!
!  Length (m) of domain box in the XI- and ETA-directions.
!
        real(r8), allocatable :: xl(:)          ! m
        real(r8), allocatable :: el(:)          ! m
!
!  Minimum and Maximum longitude and latitude at RHO-points
!
        real(r8), allocatable :: LonMin(:)      ! degrees east
        real(r8), allocatable :: LonMax(:)      ! degrees east
        real(r8), allocatable :: LatMin(:)      ! degrees north
        real(r8), allocatable :: LatMax(:)      ! degrees north
!
!  Constant used in the Shchepetkin boundary conditions for 2D momentum,
!  Co = 1.0_r8/(2.0_r8+SQRT(2.0_r8)).
!
        real(r8) :: Co
!
!  Number of digits in grid size for format statements.
!
        integer, allocatable :: Idigits(:)
        integer, allocatable :: Jdigits(:)
#ifdef SOLVE3D
        integer, allocatable :: Kdigits(:)
#endif
!
!  Diagnostic volume averaged variables.
!
        integer, allocatable :: first_time(:)

        real(r8) :: avgke = 0.0_r8              ! Kinetic energy
        real(r8) :: avgpe = 0.0_r8              ! Potential energy
        real(r8) :: avgkp = 0.0_r8              ! Total energy
        real(r8) :: volume = 0.0_r8             ! diagnostics volume
        real(r8) :: ad_volume                   ! adjoint volume
        real(r8) :: TotVolume                   ! Total volume
        real(r8) :: MaxVolume                   ! Minimum cell volume
        real(r8) :: MinVolume                   ! Maximum cell volume
!
!  Minimun and maximum grid spacing
!
        real(r8), allocatable :: DXmin(:)
        real(r8), allocatable :: DXmax(:)
        real(r8), allocatable :: DYmin(:)
        real(r8), allocatable :: DYmax(:)
#ifdef SOLVE3D
        real(r8), allocatable :: DZmin(:)
        real(r8), allocatable :: DZmax(:)
#endif
!$OMP THREADPRIVATE (DXmin, DXmax, DYmin, DYmax)
#ifdef SOLVE3D
!$OMP THREADPRIVATE (DZmin, DZmax)
#endif
!
!  Maximum size of a grid node (m) over the whole curvilinear grid
!  application. Used for scaling horizontal mixing by the grid size.
!
        real(r8), allocatable :: grdmax(:)
#ifdef DIFF_3DCOEF
        real(r8), allocatable :: DiffMin(:)     ! Minimun diffusion
        real(r8), allocatable :: DiffMax(:)     ! Maximum diffusion
#endif
#ifdef VISC_3DCOEF
        real(r8), allocatable :: ViscMin(:)     ! Minimum viscosity
        real(r8), allocatable :: ViscMax(:)     ! Maximum viscosity
#endif
!
!  Courant Numbers due to gravity wave speed limits.
!
        real(r8) :: Cg_min                      ! Minimun barotropic
        real(r8) :: Cg_max                      ! Maximun barotropic
        real(r8) :: Cg_Cor                      ! Maximun Coriolis
!
!  Time dependent Counrant Numbers due to velocity components and
!  indices location of maximum value.
!
        integer :: max_Ci, max_Cj, max_Ck       ! maximum location
        real(r8) :: max_C                       ! maximum total
        real(r8) :: max_Cu                      ! maximum I-component
        real(r8) :: max_Cv                      ! maximum J-component
#ifdef SOLVE3D
        real(r8) :: max_Cw                      ! maximum K-component
#endif
!
!  Linear equation of state parameters.
!
!    R0            Background constant density anomaly (kg/m3).
!    Tcoef         Thermal expansion coefficient (1/Celsius).
!    Scoef         Saline contraction coefficient (1/PSU).
!
        real(r8), allocatable :: R0(:)
        real(r8), allocatable :: Tcoef(:)
        real(r8), allocatable :: Scoef(:)
!
!  Background potential temperature (Celsius) and salinity (PSU) values
!  used in analytical initializations.
!
        real(r8), allocatable :: T0(:)
        real(r8), allocatable :: S0(:)
!
!  Slipperiness variable, either 1.0 (free slip) or -1.0 (no slip).
!
        real(r8), allocatable :: gamma2(:)
!
!  Weighting coefficient for the newest (implicit) time step derivatives
!  using either a Crack-Nicolson implicit scheme (lambda=0.5) or a
!  backward implicit scheme (lambda=1.0).
!
#if defined SPLINES_VDIFF || defined SPLINES_VVISC
        real(r8) :: lambda = 1.0_r8
#else
!!      real(r8) :: lambda = 0.5_r8
        real(r8) :: lambda = 1.0_r8
#endif
!
!  Jerlov water type to assign everywhere, range values: 1 - 5.
!
        integer, allocatable :: lmd_Jwt(:)
!
!  Grid r-factor (non-dimensional).
!
        real(r8) :: rx0 = 0.0_r8                ! Beckmann and Haidvogel
        real(r8) :: rx1 = 0.0_r8                ! Haney
!
!  Linear (m/s) and quadratic (nondimensional) bottom drag coefficients.
!
        real(r8), allocatable :: rdrg(:)
        real(r8), allocatable :: rdrg2(:)
!
!  Minimum and maximum threshold for transfer coefficient of momentum.
!
        real(r8) :: Cdb_min = 0.000001_r8
        real(r8) :: Cdb_max = 0.5_r8
!
!  Surface and bottom roughness (m)
!
        real(r8), allocatable :: Zos(:)
        real(r8), allocatable :: Zob(:)
!
!  Minimum depth for wetting and drying (m).
!
        real(r8), allocatable :: Dcrit(:)
!
!  Mean density (Kg/m3) used when the Boussinesq approximation is
!  inferred.
!
        real(r8) :: rho0 = 1025.0_r8
!
!  Background Brunt-Vaisala frequency (1/s2)
!
        real(r8) :: bvf_bak = 0.00001_r8

#ifdef PROPAGATOR
!
!  Number of converged Ritz values and relative accuracy of computed
!  Ritz values.
!
        integer, allocatable :: Nconv(:)
        real(r8) :: Ritz_tol = 1.0E-15_r8
#endif
!
!  Vector containing USER generic parameters.
!
        integer :: Nuser
        real(r8), dimension(25) :: user(25)
!
!  Weights for the time average of 2D fields.
!
        real(r8), allocatable :: weight(:,:,:)
!
!  Constants.
!
        real(r8), parameter :: pi = 3.14159265358979323846_r8
        real(r8), parameter :: deg2rad = pi / 180.0_r8
        real(r8), parameter :: rad2deg = 180.0_r8 / pi
        real(r8), parameter :: day2sec = 86400.0_r8
        real(r8), parameter :: sec2day = 1.0_r8 / 86400.0_r8
#ifdef NO_4BYTE_REALS
        real(r8), parameter :: spval = NF90_FILL_DOUBLE
#else
        real(r8), parameter :: spval = 1.0E+37_r8
#endif
        real(r8), parameter :: jul_off = 2440000.0_r8
!
!  Set special check value.  Notice that a smaller value is assigned
!  to account for both NetCDF fill value and roundoff. There are
!  many Matlab scripts out there that do not inquire correctly
!  the spval from the _FillValue attribute in single/double
!  precision.
!
        real(r8), parameter :: spval_check = 1.0E+35_r8
!
!-----------------------------------------------------------------------
!  Horizontal and vertical constant mixing coefficients.
!-----------------------------------------------------------------------
!
!    Akk_bak       Background vertical mixing coefficient (m2/s) for
!                    turbulent energy.
!    Akp_bak       Background vertical mixing coefficient (m2/s) for
!                    generic statistical field "psi".
!    Akt_bak       Background vertical mixing coefficient (m2/s) for
!                    tracers.
!    Akv_bak       Background vertical mixing coefficient (m2/s) for
!                    momentum.
!    Kdiff         Isopycnal mixing thickness diffusivity (m2/s) for
!                    tracers.
!    ad_visc2      ADM lateral harmonic constant mixing coefficient
!                    (m2/s) for momentum.
!    nl_visc2      NLM lateral harmonic constant mixing coefficient
!                    (m2/s) for momentum.
!    tl_visc2      TLM lateral harmonic constant mixing coefficient
!                    (m2/s) for momentum.
!    visc2         Current lateral harmonic constant mixing coefficient
!                    (m2/s) for momentum.
!    ad_visc4      ADM lateral biharmonic (squared root) constant
!                     mixing coefficient (m2 s^-1/2) for momentum.
!    nl_visc4      NLM lateral biharmonic (squared root) constant
!                     mixing coefficient (m2 s^-1/2) for momentum.
!    tl_visc4      TLM lateral biharmonic (squared root) constant
!                     mixing coefficient (m2 s^-1/2) for momentum.
!    visc4         Current lateral biharmonic (squared root) constant
!                     mixing coefficient (m2 s^-1/2) for momentum.
!    ad_tnu2       ADM lateral harmonic constant mixing coefficient
!                    (m2/s) for tracer type variables.
!    nl_tnu2       NLM lateral harmonic constant mixing coefficient
!                    (m2/s) for tracer type variables.
!    tl_tnu2       TLM lateral harmonic constant mixing coefficient
!                    (m2/s) for tracer type variables.
!    tnu2          Current lateral harmonic constant mixing coefficient
!                    (m2/s) for tracer type variables.
!    ad_tnu4       ADM lateral biharmonic (squared root) constant
!                     mixing coefficient (m2 s^-1/2) for tracers.
!    nl_tnu4       NLM lateral biharmonic (squared root) constant
!                     mixing coefficient (m2 s^-1/2) for tracers.
!    tl_tnu4       TLM lateral biharmonic (squared root) constant
!                     mixing coefficient (m2 s^-1/2) for tracers.
!    tnu4          Current lateral biharmonic (squared root) constant
!                     mixing coefficient (m2 s^-1/2) for tracers.
!    tkenu2        Lateral harmonic constant mixing coefficient
!                    (m2/s) for turbulent energy.
!    tkenu4        Lateral biharmonic (squared root) constant mixing
!                    coefficient (m2 s^-1/2) for turbulent energy.
!
        real(r8), allocatable :: Akk_bak(:)          ! m2/s
        real(r8), allocatable :: Akp_bak(:)          ! m2/s
        real(r8), allocatable :: Akv_bak(:)          ! m2/s

        real(r8), allocatable :: ad_visc2(:)         ! m2/s
        real(r8), allocatable :: nl_visc2(:)         ! m2/s
        real(r8), allocatable :: tl_visc2(:)         ! m2/s
        real(r8), allocatable :: visc2(:)            ! m2/s

        real(r8), allocatable :: ad_visc4(:)         ! m2 s-1/2
        real(r8), allocatable :: nl_visc4(:)         ! m2 s-1/2
        real(r8), allocatable :: tl_visc4(:)         ! m2 s-1/2
        real(r8), allocatable :: visc4(:)            ! m2 s-1/2

        real(r8), allocatable :: tkenu2(:)           ! m2/s
        real(r8), allocatable :: tkenu4(:)           ! m2 s-1/2

        real(r8), allocatable :: Akt_bak(:,:)        ! m2/s
        real(r8), allocatable :: Kdiff(:,:)          ! m2/s

        real(r8), allocatable :: ad_tnu2(:,:)        ! m2/s
        real(r8), allocatable :: nl_tnu2(:,:)        ! m2/s
        real(r8), allocatable :: tl_tnu2(:,:)        ! m2/s
        real(r8), allocatable :: tnu2(:,:)           ! m2/s

        real(r8), allocatable :: ad_tnu4(:,:)        ! m2 s-1/2
        real(r8), allocatable :: nl_tnu4(:,:)        ! m2 s-1/2
        real(r8), allocatable :: tl_tnu4(:,:)        ! m2 s-1/2
        real(r8), allocatable :: tnu4(:,:)           ! m2 s-1/2
!
!  Horizontal diffusive relaxation coefficients (m2/s) used to smooth
!  representer tangent linear solution during Picard iterations to
!  improve stability and convergence.
!
        real(r8), allocatable :: tl_M2diff(:)        ! 2D momentum
        real(r8), allocatable :: tl_M3diff(:)        ! 3D momentum

        real(r8), allocatable :: tl_Tdiff(:,:)       ! tracers
!
!  Basic state vertical mixing coefficient scale factors for adjoint
!  based algorithms. In some applications, a smaller/larger values of
!  vertical mixing are necessary for stability.
!
        real(r8), allocatable :: ad_Akv_fac(:)       ! ADM momentum
        real(r8), allocatable :: tl_Akv_fac(:)       ! TLM momentum

        real(r8), allocatable :: ad_Akt_fac(:,:)     ! ADM tracers
        real(r8), allocatable :: tl_Akt_fac(:,:)     ! TLM tracers

!
!  Switches to increase/decrease horizontal viscosity and/or diffusion
!  in specific areas of the application domain (like sponge areas).
!
        logical, allocatable :: Lsponge(:)
        logical, allocatable :: LuvSponge(:)         ! viscosity
        logical, allocatable :: LtracerSponge(:,:)   ! diffusion
!
!-----------------------------------------------------------------------
!  IO parameters.
!-----------------------------------------------------------------------
!
!  Switches to activate creation and writing of output NetCDF files.
!
        logical, allocatable :: LdefADJ(:)       ! Adjoint file
        logical, allocatable :: LdefAVG(:)       ! Average file
        logical, allocatable :: LdefDAI(:)       ! DA initial/restart
        logical, allocatable :: LdefAVG2(:)      ! Average2 file
        logical, allocatable :: LdefDIA(:)       ! Diagnostics file
        logical, allocatable :: LdefERR(:)       ! 4DVar error file
        logical, allocatable :: LdefFLT(:)       ! Floats file
        logical, allocatable :: LdefHIS(:)       ! History file
        logical, allocatable :: LdefHIS2(:)      ! History2 file
        logical, allocatable :: LdefHSS(:)       ! Hessian file
        logical, allocatable :: LdefINI(:)       ! Initial file
        logical, allocatable :: LdefIRP(:)       ! Initial RPM file
        logical, allocatable :: LdefITL(:)       ! Initial TLM file
        logical, allocatable :: LdefLCZ(:)       ! Lanczos Vectors file
        logical, allocatable :: LdefLZE(:)       ! Evolved Lanczos file
        logical, allocatable :: LdefMOD(:)       ! 4DVAR file
        logical, allocatable :: LdefRST(:)       ! Restart file
        logical, allocatable :: LdefSTA(:)       ! Stations file
        logical, allocatable :: LdefTIDE(:)      ! tide forcing file
        logical, allocatable :: LdefTLM(:)       ! Tangent linear file
        logical, allocatable :: LdefTLF(:)       ! TLM/RPM impulse file

        logical, allocatable :: LwrtADJ(:)       ! Write adjoint file
        logical, allocatable :: LwrtAVG(:)       ! Write average file
        logical, allocatable :: LwrtAVG2(:)      ! Write average2 file
        logical, allocatable :: LwrtDIA(:)       ! Write diagnostic file
        logical, allocatable :: LwrtHIS(:)       ! Write history file
        logical, allocatable :: LwrtHIS2(:)      ! Write history2 file
        logical, allocatable :: LwrtPER(:)       ! Write during ensemble
        logical, allocatable :: LwrtRST(:)       ! Write restart file
        logical, allocatable :: LwrtTLM(:)       ! Write tangent file
        logical, allocatable :: LwrtTLF(:)       ! Write impulse file

        logical, allocatable :: LdefNRM(:,:)     ! Norm file
        logical, allocatable :: LwrtNRM(:,:)     ! Write norm file
!
!  Switch to write out adjoint 2D state arrays instead of IO solution
!  arrays and adjoint ocean time. This is used in 4DVAR for IO
!  maniputations.
!
#if defined STOCHASTIC_OPT && !defined STOCH_OPT_WHITE
        logical, allocatable :: LwrtState3d(:)
        logical, allocatable :: SOinitial(:)
#endif
        logical, allocatable :: LwrtState2d(:)
        logical, allocatable :: LwrtTime(:)
!
!  Switch to write out adjoint surface forcing fields adjusted by the
!  4DVAR algorithms.
!
        logical, allocatable :: Ladjusted(:)
!
!  Switch to write application set-up information to standard output.
!
        logical, allocatable :: LwrtInfo(:)
!
!  Switch used to create new output NetCDF files. If TRUE, new output
!  files are created. If FALSE, data is appended to an existing output
!  files.  Used only for history, average and station files.
!
        logical, allocatable :: ldefout(:)       ! New output files
!
!  Number of timesteps between creation of new output files.
!
        integer, allocatable :: ndefADJ(:)       ! Adjoint file
        integer, allocatable :: ndefAVG(:)       ! Average file
        integer, allocatable :: ndefAVG2(:)      ! Average2 file
        integer, allocatable :: ndefDIA(:)       ! Diagnostics file
        integer, allocatable :: ndefHIS(:)       ! History file
        integer, allocatable :: ndefHIS2(:)      ! History2 file
        integer, allocatable :: ndefTLM(:)       ! Tangent linear file
!
!  Starting timestep for accumulation of output.
!
        integer, allocatable :: ntsAVG(:)        ! Average file
        integer, allocatable :: ntsAVG2(:)       ! Average2 file
        integer, allocatable :: ntsDIA(:)        ! Diagnostics file
!
!  Number of timesteps between writing of output data.
!
        integer, allocatable :: nADJ(:)          ! Adjoint file
        integer, allocatable :: nAVG(:)          ! Average file
        integer, allocatable :: nAVG2(:)         ! Average2 file
        integer, allocatable :: nDIA(:)          ! Diagnostics file
        integer, allocatable :: nFLT(:)          ! Floats file
        integer, allocatable :: nHIS(:)          ! History file
        integer, allocatable :: nHIS2(:)         ! History2 file
        integer, allocatable :: nRST(:)          ! Restart file
        integer, allocatable :: nSTA(:)          ! Stations file
        integer, allocatable :: nTLM(:)          ! Tangent linear file
!
!  Number of timesteps between print of single line information to
!  standard output.
!
        integer, allocatable :: ninfo(:)
!
!  Number of timesteps between 4DVAR adjustment of open boundaries.
!  In strong constraint 4DVAR, it is possible to open bounadies at
!  other intervals in addition to initial time. These parameters are
!  used to store the appropriate number of open boundary records in
!  output history NetCDF files.
!
!    Nbrec(:) = 1 + ntimes(:) / nOBC(:)
!
!  Here, it is assumed that nOBC is a multiple of NTIMES or greater
!  than NTIMES. If nOBC > NTIMES, only one record is stored in the
!  output history NetCDF files and the adjustment is for constant
!  open boundaries with constant correction.
!
        integer, allocatable :: nOBC(:)          ! number of timesteps
        integer, allocatable :: Nbrec(:)         ! number of records
        integer, allocatable :: OBCcount(:)      ! record counter

#ifdef ADJUST_BOUNDARY
!
!  Logical switches to process open boundary arrays during 4DVar
!  adjustments.
!
        logical, allocatable :: Lobc(:,:,:)
!
!  Time (s) of surface forcing adjustment.
!
        real(r8), allocatable :: OBC_time(:,:)
#endif
!
!  Number of timesteps between adjustment of 4DVAR surface forcing
!  fields. In strong constraint 4DVAR, it is possible to adjust surface
!  forcing fields at other intervals in addition to initial time.
!  These parameters are used to store the appropriate number of
!  surface forcing records in output history NetCDF files.
!
!    Nfrec(:) = 1 + ntimes(:) / nSFF(:)
!
!  Here, it is assumed that nSFF is a multiple of NTIMES or greater
!  than NTIMES. If nSFF > NTIMES, only one record is stored in the
!  output history NetCDF files and the adjustment is for constant
!  forcing with constant correction.
!
        integer, allocatable :: nSFF(:)          ! number of timesteps
        integer, allocatable :: Nfrec(:)         ! number of records
        integer, allocatable :: SFcount(:)       ! record counter

#ifdef ADJUST_STFLUX
!
!  Logical switches to process surface tracer fluxes during 4DVar
!  adjustments.
!
        logical, allocatable :: Lstflux(:,:)
#endif
#if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
!
!  Time (s) of surface forcing adjustment.
!
        real(r8), allocatable :: SF_time(:,:)
#endif
!
!  Restart time record to read from disk and use as the initial
!  conditions. Use nrrec=0 for new solutions. If nrrec is negative
!  (say, nrrec=-1), the model will restart from the most recent
!  time record. That is, the initialization record is assigned
!  internally.
!
        integer, allocatable :: nrrec(:)
!
!  Switch to activate processing of input data.  This switch becomes
!  very useful when reading input data serially in parallel
!  applications.
!
        logical, allocatable :: synchro_flag(:)
!$OMP THREADPRIVATE (synchro_flag)
!
!  Switch to inialize model with latest time record from initial
!  (restart/history) NetCDF file.
!
        logical, allocatable :: LastRec(:)
!
!  Generalized Statbility Theory (GST) parameters.
!
        logical :: LmultiGST          ! multiple eigenvector file switch
        logical :: LrstGST            ! restart switch
        integer :: MaxIterGST         ! Number of iterations
        integer :: nGST               ! check pointing interval
!
!  Switches used to recycle time records in some output file. If TRUE,
!  only the latest two time records are maintained.  If FALSE, all
!  field records are saved.
!
        logical, allocatable :: LcycleADJ(:)
        logical, allocatable :: LcycleRST(:)
        logical, allocatable :: LcycleTLM(:)

#if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
!
!  Counter storing the number of accumulated harmonic records used
!  for detiding.
!
        integer, allocatable :: Hcount(:)
#endif
!
!-----------------------------------------------------------------------
!  Adjoint sensitivity parameters.
!-----------------------------------------------------------------------
!
!  Starting and ending vertical levels of the 3D adjoint state whose
!  sensitivity is required.
!
        integer, allocatable :: KstrS(:)           ! starting level
        integer, allocatable :: KendS(:)           ! ending level
!
!  Starting and ending day for adjoint sensitivity forcing.
!
        real(r8), allocatable :: DstrS(:)          ! starting day
        real(r8), allocatable :: DendS(:)          ! ending day
!
!-----------------------------------------------------------------------
!  Stochastic optimals parameters.
!-----------------------------------------------------------------------
!
!  Stochastic optimals forcing records counter.
!
        integer, allocatable :: SOrec(:)
!$OMP THREADPRIVATE (SOrec)
!
!  Trace of stochastic optimals matrix.
!
        real(r8), allocatable :: TRnorm(:)
!
!  Stochastic optimals time decorrelation scale (days) assumed for
!  red noise processes.
!
        real(r8), allocatable :: SO_decay(:)
!
!  Stochastic optimals surface forcing standard deviation for
!  dimensionalization.
!
        real(r8), allocatable :: SO_sdev(:,:)

#if defined FOUR_DVAR || defined VERIFICATION
!
!------------------------------------------------------------------------
!  Background/model error covariance parameters.
!------------------------------------------------------------------------
!
!  Maximum number of model state variables to process.
!
        integer :: MstateVar
!
!  Logical switch to compute initial conditions, model and surface
!  forcing error covariance normalization factors.
!
        logical,  allocatable :: Cnorm(:,:)
!
!  Logical switch to compute boundary conditions error covariance
!  normalization factors.
!
        logical,  allocatable :: CnormB(:,:)
!
!  Logical switches to process weak constraint forcing as intermittent
!  or continuous impulses.
!
        logical, allocatable :: SporadicImpulse(:)     ! intermittent
        logical, allocatable :: FrequentImpulse(:)     ! continuous
!
!  Stability and accuracy factor used to scale the time-step of the
!  horizontal and vertical convolution operator below its theoretical
!  (CFL) limit. Notice that four values are needed for these factors to
!  facilitate the error covariance modeling for initial conditions (1),
!  model (2), boundary conditions (3), and surface forcing (4).
!
        real(r8), dimension(4) :: Hgamma
        real(r8), dimension(4) :: Vgamma
!
!  Parameters used to compute balanced salinity in terms of temperature
!  using empirical T-S relationships in error covariance balance
!  operator.
!
        real(r8), allocatable :: dTdz_min(:)     ! minimum dT/ds (C/m)
        real(r8), allocatable :: ml_depth(:)     ! mixed-layer depth (m)
!
!  Balance operator level of no motion depth (m; positive) used when
!  computing balanced free-surface contribution.
!
        real(r8), allocatable :: LNM_depth(:)
!
! Balance operator level of no motion flag used to compute balanced
! free-surface contribution:
!
!      [0] Integrate from local bottom to the surface
!      [1] Integrate from LNM_depth to surface or integrate from local bottom
!          if shallower than LNM_depth
!
        integer :: LNM_flag
!
! Balance operator logical switches for state variables to consider in the
! error covariance multivariate constraints.
!
        logical, allocatable :: balance(:)
!
!  Initial conditions, model and surface forcing error covariance
!  horizontal decorrelation scales (m).
!
        real(r8), allocatable :: Hdecay(:,:,:)
!
!  Initial conditions, model and surface forcing error covariance
!  temporal decorrelation scales (second).
!
        real(r8), allocatable :: Tdecay(:,:)
!
!  Initial conditions, model and surface forcing error covariance
!  vertical decorrelation scales (m).
!
        real(r8), allocatable :: Vdecay(:,:,:)
!
!  Boundary conditions error covariance horizontal decorrelation
!  scales (m).
!
        real(r8), allocatable :: HdecayB(:,:,:)
!
!  Boundary conditions error covariance vertical decorrelation
!  scales (m).
!
        real(r8), allocatable :: VdecayB(:,:,:)
!
!  Method for background quality control of observations, [Ngrids].
!
!      [1] Quality control in terms of state variable indices
!      [2] Quality control in terms of observation provenance
!
        integer, allocatable :: bgqc_type(:)
!
!  Background quality control standard deviation value for not
!  rejection of observations.
!
        real(r8) :: bgqc_large = 1.0E+5_r8
!
!  Number of observation provenances used in background quality control,
!  [Ngrids]. Only used when bgqc_type(ng)=2.
!
        integer, allocatable :: Nprovenance(:)
!
!  Observation provenance indices to process during background quality
!  control of observations, [MAXVAL(Nprovenance),Ngrids].
!
        integer, allocatable :: Iprovenance(:,:)
!
!  Background quality control threshold standard deviations in terms
!  of observation provenance indices, [MAXVAL(Nprovenance),Ngrids].
!
        real(r8), allocatable :: P_bgqc(:,:)
!
!  Background quality control threshold standard deviations in terms
!  of state variable indices, [MstateVar,Ngrids]
!
        real(r8), allocatable :: S_bgqc(:,:)
#endif
!
!-----------------------------------------------------------------------
!  Nudging variables for passive (outflow) and active (inflow) oepn
!  boundary conditions.
!-----------------------------------------------------------------------
!
!    iwest         West  identification index in boundary arrays.
!    isouth        South identification index in boundary arrays.
!    ieast         East  identification index in boundary arrays.
!    inorth        North identification index in boundary arrays.
!    obcfac        Factor between passive and active open boundary
!                    conditions (nondimensional and greater than one).
!                    The nudging time scales for the active conditions
!                    are obtained by multiplying the passive values by
!                    factor.
!    FSobc_in      Active and strong time-scale (1/sec) coefficients
!                    for nudging towards free-surface data at  inflow.
!    FSobc_out     Passive and weak  time-scale (1/sec) coefficients
!                    for nudging towards free-surface data at outflow.
!    M2obc_in      Active and strong time-scale (1/sec) coefficients
!                    for nudging towards 2D momentum data at  inflow.
!    M2obc_out     Passive and weak  time-scale (1/sec) coefficients
!                    for nudging towards 2D momentum data at outflow.
!    M3obc_in      Active and strong time-scale (1/sec) coefficients
!                    for nudging towards 3D momentum data at  inflow.
!    M3obc_out     Passive and weak  time-scale (1/sec) coefficients
!                    for nudging towards 3D momentum data at outflow.
!    Tobc_in       Active and strong time-scale (1/sec) coefficients
!                    for nudging towards tracer data at  inflow.
!    Tobc_out      Passive and weak  time-scale (1/sec) coefficients
!                    for nudging towards tracer data at outflow.
!
        integer, parameter :: iwest = 1
        integer, parameter :: isouth = 2
        integer, parameter :: ieast = 3
        integer, parameter :: inorth = 4

        real(r8), allocatable :: obcfac(:)
        real(r8), allocatable :: FSobc_in(:,:)
        real(r8), allocatable :: FSobc_out(:,:)
        real(r8), allocatable :: M2obc_in(:,:)
        real(r8), allocatable :: M2obc_out(:,:)
#ifdef SOLVE3D
        real(r8), allocatable :: M3obc_in(:,:)
        real(r8), allocatable :: M3obc_out(:,:)
        real(r8), allocatable :: Tobc_in(:,:,:)
        real(r8), allocatable :: Tobc_out(:,:,:)
#endif
!
!  Inverse time-scales (1/s) for nudging at open boundaries and sponge
!  areas.
!
        real(r8), allocatable :: Znudg(:)          ! Free-surface
        real(r8), allocatable :: M2nudg(:)         ! 2D momentum
        real(r8), allocatable :: M3nudg(:)         ! 3D momentum
        real(r8), allocatable :: Tnudg(:,:)        ! Tracers
        real(r8), allocatable :: Tnudg_SSS(:)      ! Sea surface salinity
        real(r8), allocatable :: sss_mismatch_threshold(:)! Sea surface salinity
!
!  Variables used to impose mass flux conservation in open boundary
!  configurations.
!
        real(r8) :: bc_area = 0.0_r8
        real(r8) :: bc_flux = 0.0_r8
        real(r8) :: ubar_xs = 0.0_r8
#if defined TANGENT || defined TL_IOMS
        real(r8) :: tl_bc_area = 0.0_r8
        real(r8) :: tl_bc_flux = 0.0_r8
        real(r8) :: tl_ubar_xs = 0.0_r8
#endif
#ifdef ADJOINT
        real(r8) :: ad_bc_area = 0.0_r8
        real(r8) :: ad_bc_flux = 0.0_r8
        real(r8) :: ad_ubar_xs = 0.0_r8
#endif

#if defined BULK_FLUXES || defined BULK_FLUXES2D
!
!-----------------------------------------------------------------------
! Constants used in surface fluxes bulk parameterization.
!-----------------------------------------------------------------------
!
!    blk_Cpa       Specific heat capacity for dry air (J/kg/K).
!    blk_Cpw       Specific heat capacity for seawater (J/kg/K).
!    blk_Rgas      Gas constant for dry air (J/kg/K).
!    blk_Zabl      Height (m) of atmospheric boundary layer.
!    blk_ZQ        Height (m) of surface air humidity measurement.
!    blk_ZT        Height (m) of surface air temperature measurement.
!    blk_ZW        Height (m) of surface winds measurement.
!    blk_beta      Beta parameter evaluated from Fairall low windspeed
!                    turbulence data.
!    blk_dter      Temperature change.
!    blk_tcw       Thermal conductivity of water (W/m/K).
!    blk_visw      Kinematic viscosity water (m2/s).
!
        real(r8) :: blk_Cpa = 1004.67_r8      ! (J/kg/K), Businger 1982
        real(r8) :: blk_Cpw = 4000.0_r8       ! (J/kg/K)
        real(r8) :: blk_Rgas = 287.1_r8       ! (J/kg/K)
        real(r8) :: blk_Zabl = 600.0_r8       ! (m)
        real(r8) :: blk_beta = 1.2_r8         ! non-dimensional
        real(r8) :: blk_dter = 0.3_r8         ! (K)
        real(r8) :: blk_tcw = 0.6_r8          ! (W/m/K)
        real(r8) :: blk_visw = 0.000001_r8    ! (m2/s)

        real(r8), allocatable :: blk_ZQ(:)        ! (m)
        real(r8), allocatable :: blk_ZT(:)        ! (m)
        real(r8), allocatable :: blk_ZW(:)        ! (m)
#endif

#ifdef ICESHELF
!
!---------------------------------------------------------------------
!  Define parameters associated with the transfer of heat and salt
!  beneath the ice shelf (values from Beckmann, et al. (1999) and
!  Hellmer and Olbers (1989)).  Note that gamma_s and gamma_t
!  can be changed if ICESHELF_3EQ defined.
!---------------------------------------------------------------------
!
!    gamma_s       Turbulent salt exchange coefficient underneath
!                    the ice (m/s)
!    gamma_t       Turbulent heat exchange coefficient underneath
!                    the ice (m/s)
!    Hlfreeze      Latent heat for freezing salt water (J/kg)
!    rho_ice       Density of sea-ice (kg/m^3)
!    blk_Cpi       Specific heat capicity of ice (J/kg/deg C)
!
        real(r8) :: gamma_s = 5.05E-7_r8       ! (m/s)
        real(r8) :: gamma_t = 1.00E-4_r8       ! (m/s)
        real(r8) :: Hlfreeze = 334000.0_r8     ! (J/kg)
        real(r8) :: rho_ice = 930.0_r8         ! (kg/m^3)
        real(r8) :: blk_Cpi = 2000.0_r8        ! (J/kg/C)
#endif

# if defined SG_BBL || defined SSW_BBL
!
!-----------------------------------------------------------------------
!  Closure parameters associated with Styles and Glenn (1999) bottom
!  currents and waves boundary layer.
!-----------------------------------------------------------------------
!
!    sg_Cdmax      Upper limit on bottom darg coefficient.
!    sg_alpha      Free parameter indicating the constant stress
!                    region of the wave boundary layer.
!    sg_g          Acceleration of gravity (m/s2).
!    sg_kappa      Von Karman constant.
!    sg_mp         Nondimensional closure constant.
!    sg_n          Maximum number of iterations for bisection method.
!    sg_nu         Kinematic viscosity of seawater (m2/s).
!    sg_pi         Ratio of circumference to diameter.
!    sg_tol        Convergence criterion.
!    sg_ustarcdef  Default bottom stress (m/s).
!    sg_z100       Depth (m), 100 cm above bottom.
!    sg_z1p        Nondimensional closure constant.
!    sg_zrmin      Minimum allowed height (m) of current above bed.
!                    Otherwise, logarithmic interpolation is used.
!    sg_znotcdef   Default apparent hydraulic roughness (m).
!    sg_znotdef    Default hydraulic roughness (m).
!
        integer, parameter :: sg_n = 20

        real(r8), parameter :: sg_pi = pi

        real(r8) :: sg_Cdmax = 0.01_r8         ! non-dimensional
        real(r8) :: sg_alpha = 1.0_r8          ! non-dimensional
        real(r8) :: sg_g = 9.81_r8             ! (m/s2)
        real(r8) :: sg_kappa = 0.41_r8         ! non-dimensional
        real(r8) :: sg_nu = 0.00000119_r8      ! (m2/s)
        real(r8) :: sg_tol = 0.0001_r8         ! non-dimensional
        real(r8) :: sg_ustarcdef = 0.01_r8     ! (m/s)
        real(r8) :: sg_z100 = 1.0_r8           ! (m)
        real(r8) :: sg_z1min = 0.20_r8         ! (m)
        real(r8) :: sg_z1p                     ! non-dimensional
        real(r8) :: sg_znotcdef = 0.01_r8      ! (m)
        real(r8) :: sg_znotdef                 ! (m)

        complex(c8) :: sg_mp
# endif
#if defined LMD_SKPP || defined SOLAR_SOURCE
!
!-----------------------------------------------------------------------
!  Water clarity parameters.
!-----------------------------------------------------------------------
!
!  The water type classification is based on Jerlov water type using
!  a double exponential function for light absorption:
!
!    Array
!    Index   WaterType   Examples
!    -----   ---------   --------
!
!      1         I       Open Pacific
!      2         IA      Eastern Mediterranean, Indian Ocean
!      3         IB      Western Mediterranean, Open Atlantic
!      4         II      Coastal waters, Azores
!      5         III     Coastal waters, North Sea
!      6         1       Skagerrak Strait
!      7         3       Baltic
!      8         5       Black Sea
!      9         7       Dark coastal water
!
!    lmd_mu1       Reciprocal of the absorption coefficient for solar
!                    wavelength band 1 as a function of the Jerlov
!                    water type.
!    lmd_mu2       Reciprocal of the absorption coefficient for solar
!                    wavelength band 2 as a function of the Jerlov
!                    water type.
!    lmd_r1        Fraction of total radiance for wavelength band 1 as
!                    a function of the Jerlov water type.
!
        real(r8), dimension(9) :: lmd_mu1 =                             &
     &            (/ 0.35_r8, 0.6_r8,  1.0_r8,  1.5_r8, 1.4_r8,         &
     &               0.42_r8, 0.37_r8, 0.33_r8, 0.0_r8 /)

        real(r8), dimension(9) :: lmd_mu2 =                             &
     &            (/ 23.0_r8,  20.0_r8,  17.0_r8,  14.0_r8, 7.9_r8,     &
     &                5.13_r8,  3.54_r8,  2.34_r8,  1.51_r8 /)

        real(r8), dimension(9) :: lmd_r1 =                              &
     &            (/ 0.58_r8, 0.62_r8, 0.67_r8, 0.77_r8, 0.78_r8,       &
     &               0.57_r8, 0.57_r8, 0.57_r8, 0.55_r8 /)
#endif

#ifdef LMD_MIXING
!
!-----------------------------------------------------------------------
!  Large et al. (1994) K-profile parameterization.
!-----------------------------------------------------------------------
!
!    lmd_Ri0       Critical gradient Richardson number below which
!                    turbulent mixing occurs.
!    lmd_Rrho0     Value of double-diffusive density ratio where
!                    mixing goes to zero in salt fingering.
!    lmd_bvfcon    Brunt-Vaisala frequency (1/s2) limit for convection.
!    lmd_fdd       Scaling factor for double diffusion of temperature
!                    in salt fingering case (lmd_fdd=0.7).
!    lmd_nu        Molecular viscosity (m2/s).
!    lmd_nu0c      Maximum interior convective viscosity and diffusivity
!                    due to shear instability.
!    lmd_nu0m      Maximum interior viscosity (m2/s) due shear
!                    instability.
!    lmd_nu0s      Maximum interior diffusivity (m2/s) due shear
!                    instability.
!    lmd_nuf       Scaling factor for double diffusion in salt
!                    fingering.
!    lmd_nuwm      Interior viscosity (m2/s) due to wave breaking.
!    lmd_nuws      Interior diffusivity (m2/s) due to wave breaking.
!    lmd_sdd1      Double diffusion constant for salinity in diffusive
!                    convection case (lmd_sdd1=0.15).
!    lmd_sdd2      Double diffusion constant for salinity in diffusive
!                    convection case (lmd_sdd2=1.85).
!    lmd_sdd3      Double diffusion constant for salinity in diffusive
!                    convection case (lmd_sdd3=0.85).
!    lmd_tdd1      Double diffusion constant for temperature
!                    in diffusive convection case (lmd_tdd1=0.909).
!    lmd_tdd2      Double diffusion constant for temperature in
!                    diffusive convection case (lmd_tdd2=4.6).
!    lmd_tdd3      Double diffusion constant for temperature in
!                    diffusive convection case (lmd_tdd3=0.54).
!
        real(r8) :: lmd_Ri0 = 0.7_r8          ! non-dimensional
        real(r8) :: lmd_Rrho0 = 1.9_r8        ! m2/s
        real(r8) :: lmd_bvfcon = -2.0E-5_r8   ! 1/s2
        real(r8) :: lmd_fdd = 0.7_r8          ! non-dimensional
        real(r8) :: lmd_nu = 1.5E-6_r8        ! m2/s
        real(r8) :: lmd_nu0c = 0.01_r8        ! m2/s
!!      real(r8) :: lmd_nu0c = 0.1_r8         ! m2/s
!!      real(r8) :: lmd_nu0c = 0.05_r8        ! m2/s
        real(r8) :: lmd_nu0m = 10.0E-4_r8     ! m2/s
        real(r8) :: lmd_nu0s = 10.0E-4_r8     ! m2/s
!!      real(r8) :: lmd_nu0m = 50.0E-4_r8     ! m2/s
!!      real(r8) :: lmd_nu0s = 50.0E-4_r8     ! m2/s
        real(r8) :: lmd_nuf = 10.0E-4_r8      ! m2/s
# ifdef DAMEE_4
        real(r8) :: lmd_nuwm = 1.0E-4_r8      ! m2/s
        real(r8) :: lmd_nuws = 1.0E-5_r8      ! m2/s
# else
        real(r8) :: lmd_nuwm = 1.0E-5_r8      ! m2/s
        real(r8) :: lmd_nuws = 1.0E-6_r8      ! m2/s
# endif
        real(r8) :: lmd_sdd1 = 0.15_r8        ! non-dimensional
        real(r8) :: lmd_sdd2 = 1.85_r8        ! non-dimensional
        real(r8) :: lmd_sdd3 = 0.85_r8        ! non-dimensional
        real(r8) :: lmd_tdd1 = 0.909_r8       ! non-dimensional
        real(r8) :: lmd_tdd2 = 4.6_r8         ! non-dimensional
        real(r8) :: lmd_tdd3 = 0.54_r8        ! non-dimensional

# if defined LMD_SKPP || defined LMD_BKPP
!
!-----------------------------------------------------------------------
!  Large et al. (1994) oceanic boundary layer parameters.
!-----------------------------------------------------------------------
!
!    lmd_Cg        Proportionality coefficient parameterizing nonlocal
!                    transport.
!    lmd_Cstar     Proportionality coefficient parameterizing nonlocal
!                    transport.
!    lmd_Cv        Ratio of interior Brunt-Vaisala frequency to that
!                    at entrainment depth "he".
!    lmd_Ric       Critical bulk Richardson number.
!    lmd_am        Coefficient of flux profile for momentum in their
!                    1/3 power law regimes.
!    lmd_as        Coefficient of flux profile for tracers in their
!                    1/3 power law regimes.
!    lmd_betaT     Ratio of entrainment flux to surface buoyancy flux.
!    lmd_cekman    Constant used in the computation of Ekman depth.
!    lmd_cmonob    Constant used in the computation of Monin-Obukhov
!                    depth.
!    lmd_cm        Coefficient of flux profile for momentum in their
!                    1/3 power law regimes.
!    lmd_cs        Coefficient of flux profile for tracers in their
!                    1/3 power law regimes.
!    lmd_epsilon   Non-dimensional extent of the surface layer.
!    lmd_zetam     Maximum stability parameter "zeta" value of the 1/3
!                    power law regime of flux profile for momentum.
!    lmd_zetas     Maximum stability parameter "zeta" value of the 1/3
!                    power law regime of flux profile for tracers.
!
        real(r8) :: lmd_Cg
        real(r8) :: lmd_Cstar = 10.0_r8
        real(r8) :: lmd_Cv = 1.25_r8
!!      real(r8) :: lmd_Cv = 1.4_r8
!!      real(r8) :: lmd_Cv = 1.6_r8
!!      real(r8) :: lmd_Cv = 1.8_r8
!!      real(r8) :: lmd_Ric = 0.25_r8
        real(r8) :: lmd_Ric = 0.3_r8
!!      real(r8) :: lmd_Ric = 0.5_r8
!!      real(r8) :: lmd_Ric = 0.75_r8
        real(r8) :: lmd_am = 1.257_r8
        real(r8) :: lmd_as = -28.86_r8
        real(r8) :: lmd_betaT = -0.2_r8
        real(r8) :: lmd_cekman = 0.7_r8
        real(r8) :: lmd_cmonob = 1.0_r8
        real(r8) :: lmd_cm = 8.36_r8
        real(r8) :: lmd_cs = 98.96_r8
        real(r8) :: lmd_epsilon = 0.1_r8
        real(r8) :: lmd_zetam = -0.2_r8
        real(r8) :: lmd_zetas = -1.0_r8
# endif
#endif
!
!-----------------------------------------------------------------------
!  Generic Length Scale parameters.
!-----------------------------------------------------------------------
!
!    gls_Gh0
!    gls_Ghcri
!    gls_Ghmin
!    gls_Kmin      Minimum value of specific turbulent kinetic energy.
!    gls_Pmin      Minimum Value of dissipation.
!    gls_cmu0      Stability coefficient (non-dimensional).
!    gls_c1        Shear production coefficient (non-dimensional).
!    gls_c2        Dissipation coefficient (non-dimensional).
!    gls_c3m       Buoyancy production coefficient (minus).
!    gls_c3p       Buoyancy production coefficient (plus).
!    gls_E2
!    gls_m         Turbulent kinetic energy exponent (non-dimensional).
!    gls_n         Turbulent length scale exponent (non-dimensional).
!    gls_p         Stability exponent (non-dimensional).
!    gls_sigk      Constant Schmidt number (non-dimensional) for
!                    turbulent kinetic energy diffusivity.
!    gls_sigp      Constant Schmidt number (non-dimensional) for
!                    turbulent generic statistical field, "psi".
!
        real(r8), allocatable :: gls_m(:)
        real(r8), allocatable :: gls_n(:)
        real(r8), allocatable :: gls_p(:)
        real(r8), allocatable :: gls_sigk(:)
        real(r8), allocatable :: gls_sigp(:)
        real(r8), allocatable :: gls_cmu0(:)
        real(r8), allocatable :: gls_cmupr(:)
        real(r8), allocatable :: gls_c1(:)
        real(r8), allocatable :: gls_c2(:)
        real(r8), allocatable :: gls_c3m(:)
        real(r8), allocatable :: gls_c3p(:)
        real(r8), allocatable :: gls_Kmin(:)
        real(r8), allocatable :: gls_Pmin(:)
#ifdef GLS_MIXING
# if defined CANUTO_A || defined CANUTO_B
        real(r8) :: gls_s0
        real(r8) :: gls_s1
        real(r8) :: gls_s2
        real(r8) :: gls_s3
        real(r8) :: gls_s4
        real(r8) :: gls_s5
        real(r8) :: gls_s6
        real(r8) :: gls_b0
        real(r8) :: gls_b1
        real(r8) :: gls_b2
        real(r8) :: gls_b3
        real(r8) :: gls_b4
        real(r8) :: gls_b5
# endif
# ifdef CANUTO_A
        real(r8), parameter :: gls_Gh0 = 0.0329_r8 ! 0.0329 GOTM, 0.0673 Burchard
        real(r8), parameter :: gls_Ghcri = 0.03_r8
        real(r8), parameter :: gls_L1 = 0.107_r8
        real(r8), parameter :: gls_L2 = 0.0032_r8
        real(r8), parameter :: gls_L3 = 0.0864_r8
        real(r8), parameter :: gls_L4 = 0.12_r8
        real(r8), parameter :: gls_L5 = 11.9_r8
        real(r8), parameter :: gls_L6 = 0.4_r8
        real(r8), parameter :: gls_L7 = 0.0_r8
        real(r8), parameter :: gls_L8 = 0.48_r8
# elif defined CANUTO_B
        real(r8), parameter :: gls_Gh0 = 0.0444_r8 ! 0.044 GOTM, 0.0673 Burchard
        real(r8), parameter :: gls_Ghcri = 0.0414_r8
        real(r8), parameter :: gls_L1 = 0.127_r8
        real(r8), parameter :: gls_L2 = 0.00336_r8
        real(r8), parameter :: gls_L3 = 0.0906_r8
        real(r8), parameter :: gls_L4 = 0.101_r8
        real(r8), parameter :: gls_L5 = 11.2_r8
        real(r8), parameter :: gls_L6 = 0.4_r8
        real(r8), parameter :: gls_L7 = 0.0_r8
        real(r8), parameter :: gls_L8 = 0.318_r8
# else
        real(r8), parameter :: gls_Gh0 = 0.028_r8
        real(r8), parameter :: gls_Ghcri = 0.02_r8
# endif
        real(r8), parameter :: gls_Ghmin = -0.28_r8
        real(r8), parameter :: gls_E2 = 1.33_r8
#endif
!
! Constants used in the various formulation of surface flux boundary
! conditions for the GLS vertical turbulence closure in terms of
! Charnok surface roughness (CHARNOK_ALPHA), roughness from wave
! amplitude (zos_hsig_alpha), wave dissipation (SZ_ALPHA), and
! Craig and Banner wave breaking (CRGBAN_CW).
! Wec_alpha partitions energy to roller or breaking.
        real(r8), allocatable :: charnok_alpha(:)
        real(r8), allocatable :: zos_hsig_alpha(:)
        real(r8), allocatable :: sz_alpha(:)
        real(r8), allocatable :: crgban_cw(:)
        real(r8), allocatable :: wec_alpha(:)

#if defined MY25_MIXING || defined GLS_MIXING
!
!-----------------------------------------------------------------------
!  Mellor-Yamada (1982) Level 2.5 vertical mixing variables.
!-----------------------------------------------------------------------
!
!    my_A1         Turbulent closure A1 constant.
!    my_A2         Turbulent closure A2 constant.
!    my_B1         Turbulent closure B1 constant.
!    my_B1p2o3     B1**(2/3).
!    my_B1pm1o3    B1**(-1/3).
!    my_B2         Turbulent closure B2 constant.
!    my_C1         Turbulent closure C1 constant.
!    my_C2         Turbulent closure C2 constant.
!    my_C3         Turbulent closure C3 constant.
!    my_E1         Turbulent closure E1 constant.
!    my_E1o2       0.5*E1
!    my_E2         Turbulent closure E2 constant.
!    my_Gh0        Lower bound on Galperin et al. stability function.
!    my_Sh1        Tracers stability function constant factor.
!    my_Sh2        Tracers stability function constant factor.
!    my_Sm1        Momentum stability function constant factor.
!    my_Sm2        Momentum stability function constant factor.
!    my_Sm3        Momentum stability function constant factor.
!    my_Sm4        Momentum stability function constant factor.
!    my_Sq         Scale for vertical mixing of turbulent energy.
!    my_dtfac      Asselin time filter coefficient.
!    my_lmax       Upper bound on the turbulent length scale.
!    my_qmin       Lower bound on turbulent energy "tke" and "gls".
!
        real(r8), parameter :: my_A1 = 0.92_r8
        real(r8), parameter :: my_A2 = 0.74_r8
        real(r8), parameter :: my_B1 = 16.6_r8
        real(r8), parameter :: my_B2 = 10.1_r8
        real(r8), parameter :: my_C1 = 0.08_r8
        real(r8), parameter :: my_C2 = 0.7_r8
        real(r8), parameter :: my_C3 = 0.2_r8
        real(r8), parameter :: my_E1 = 1.8_r8
        real(r8), parameter :: my_E2 = 1.33_r8
        real(r8), parameter :: my_Gh0 = 0.0233_r8
        real(r8), parameter :: my_Sq = 0.2_r8
        real(r8), parameter :: my_dtfac = 0.05_r8
        real(r8), parameter :: my_lmax = 0.53_r8
        real(r8), parameter :: my_qmin = 1.0E-8_r8

        real(r8) :: my_B1p2o3
        real(r8) :: my_B1pm1o3
        real(r8) :: my_E1o2
        real(r8) :: my_Sh1
        real(r8) :: my_Sh2
        real(r8) :: my_Sm1
        real(r8) :: my_Sm2
        real(r8) :: my_Sm3
        real(r8) :: my_Sm4
#endif
#ifdef BVF_MIXING
!
!-----------------------------------------------------------------------
!  Brunt-Vaisala frequency based vertical mixing variables.
!-----------------------------------------------------------------------
!
!    bvf_numax     Upper bound vertical diffusion (m2/s).
!    bvf_numin     Lower bound vertical diffusion (m2/s).
!    bvf_nu0       Proportionality constant (m2/s2).
!    bvf_nu0c      Convective diffusion (m2/s) in static unstable
!                    regime.
!
        real(r8) :: bvf_numax = 4.0E-4_r8      ! m2/s
        real(r8) :: bvf_numin = 3.0E-5_r8      ! m2/s
        real(r8) :: bvf_nu0 = 1.0E-7_r8        ! m2/s2
        real(r8) :: bvf_nu0c = 1.0_r8          ! m2/s
#endif
!
!-----------------------------------------------------------------------
!  Tangent linear and adjoint model parameters.
!-----------------------------------------------------------------------
!
!  Tangent linear and adjoint model control switches.
!
        logical :: TLmodel = .FALSE.
        logical :: ADmodel = .FALSE.

#ifdef FILTERED
        integer :: filtindx
#endif

#ifdef ICE_MODEL
        integer, allocatable :: nstrs(:)
        integer, allocatable :: nevp(:)
        integer, allocatable :: ievp(:)

        real(r8) :: tol     =  1.0E-16
        real(r8) :: ice_emiss =  0.97_r8
        real(r8) :: spec_heat_air =  1004._r8
        real(r8) :: trans_coeff = 1.75E-3
        real(r8) :: sublim_latent_heat = 2.834E+6
        real(r8) :: t0deg = 273.15_r8

        real(r8), allocatable :: dtice(:)
        real(r8), allocatable :: dte(:)
        real(r8), allocatable :: rhoice(:)
        real(r8), allocatable :: cdiw(:)
        real(r8), allocatable :: cdai(:)
        real(r8), allocatable :: rho_air(:)
        real(r8), allocatable :: rhosnow_dry(:)
        real(r8), allocatable :: rhosnow_wet(:)
        real(r8), allocatable :: pstar(:)
        real(r8), allocatable :: pstar_e(:)
        real(r8), allocatable :: zetamax(:)
        real(r8), allocatable :: zetamin(:)
        real(r8), allocatable :: ellip_sq(:)
        real(r8), allocatable :: astren(:)
        real(r8), allocatable :: alphai(:)
        real(r8), allocatable :: min_h(:)
        real(r8), allocatable :: min_a(:)
        real(r8), allocatable :: max_a(:)
        real(r8), allocatable :: stressang(:)
# ifdef ICE_LANDFAST
        real(r8), allocatable :: lf_k1(:)
        real(r8), allocatable :: lf_k2(:)
        real(r8), allocatable :: lf_u0(:)
# endif
#endif

      CONTAINS
!
      SUBROUTINE allocate_scalars
!
!=======================================================================
!                                                                      !
!  This routine allocates structure and  several variables in module   !
!  that depend on the number of nested grids.                          !
!                                                                      !
!=======================================================================
!
!  Local variable declarations.
!
      integer :: ng

      real(r8), parameter :: IniVal = 0.0_r8
!
!-----------------------------------------------------------------------
!  Allocate and initialize variables in module structure.
!-----------------------------------------------------------------------
!
      allocate ( SCALARS(Ngrids) )

      DO ng=1,Ngrids

#ifdef SOLVE3D
        allocate ( SCALARS(ng) % Fstate(7+2*MT) )
        SCALARS(ng) % Fstate(1:7+2*MT) = .FALSE.

        allocate ( SCALARS(ng) % Lstate(5+MT) )
        SCALARS(ng) % Lstate(1:5+MT) = .FALSE.
#else
        allocate ( SCALARS(ng) % Fstate(5) )
        SCALARS(ng) % Fstate(1:5) = .FALSE.

        allocate ( SCALARS(ng) % Lstate(3) )
        SCALARS(ng) % Lstate(1:3) = .FALSE.
#endif
#if defined READ_WATER && defined MASKING && defined DISTRIBUTE
        allocate ( SCALARS(ng) % IJwater((Lm(ng)+2)*(Mm(ng)+2),4) )
        SCALARS(ng) % IJwater = 0
#endif

        allocate ( SCALARS(ng) % Cs_r(N(ng)) )
        SCALARS(ng) % Cs_r(1:N(ng)) = IniVal

        allocate ( SCALARS(ng) % Cs_w(0:N(ng)) )
        SCALARS(ng) % Cs_w(0:N(ng)) = IniVal

        allocate ( SCALARS(ng) % sc_r(N(ng)) )
        SCALARS(ng) % sc_r(1:N(ng)) = IniVal

        allocate ( SCALARS(ng) % sc_w(0:N(ng)) )
        SCALARS(ng) % sc_w(0:N(ng)) = IniVal

      END DO
!
!  Allocate variables that require special treatment in shared-memory.
!  These variables are private for each thread to avoid collisions.
!
!$OMP PARALLEL
      allocate ( PREDICTOR_2D_STEP(Ngrids) )

      allocate ( indx1(Ngrids) )
      allocate ( iic(Ngrids) )
      allocate ( iif(Ngrids) )

      allocate ( DXmin(Ngrids) )
      allocate ( DXmax(Ngrids) )
      allocate ( DYmin(Ngrids) )
      allocate ( DYmax(Ngrids) )
#ifdef SOLVE3D
      allocate ( DZmin(Ngrids) )
      allocate ( DZmax(Ngrids) )
#endif

      allocate ( FrcRec(Ngrids) )
      allocate ( SOrec(Ngrids) )
!$OMP END PARALLEL
!
!-----------------------------------------------------------------------
!  Allocate variables.
!-----------------------------------------------------------------------
!
#if defined FOUR_DVAR || defined VERIFICATION
      MstateVar=5+MT
# ifdef ADJUST_WSTRESS
      MstateVar=MstateVar+2
# endif
# ifdef ADJUST_STFLUX
      MstateVar=MstateVar+MT
# endif
!
#endif

#ifdef T_PASSIVE
      allocate ( inert(NPT) )
#endif
      allocate ( PerfectRST(Ngrids) )

      allocate ( ndtfast(Ngrids) )
      allocate ( nfast(Ngrids) )

      allocate ( dt(Ngrids) )
      allocate ( dtfast(Ngrids) )
      allocate ( TimeEnd(Ngrids) )
      allocate ( AVGtime(Ngrids) )
      allocate ( AVG2time(Ngrids) )
      allocate ( DIAtime(Ngrids) )
      allocate ( IMPtime(Ngrids) )
      allocate ( ObsTime(Ngrids) )
      allocate ( FrcTime(Ngrids) )
#ifdef TIDES_ASTRO
      allocate ( FIRST_TIDES_ASTRO(Ngrids) )
#endif
#ifdef AGE_DISTRIBUTION
      allocate ( dtau(Ngrids) )
#endif
#ifdef CICE_MODEL
      allocate ( tspy(Ngrids) )
      allocate ( tspd(Ngrids) )
      allocate ( rhoice(Ngrids) )
      allocate ( min_a(Ngrids) )
#endif

      allocate ( ntimes(Ngrids) )
      allocate ( first_time(Ngrids) )

      allocate ( ntfirst(Ngrids) )
      allocate ( ntstart(Ngrids) )
      allocate ( ntend(Ngrids) )
!$OMP PARALLEL
      allocate ( synchro_flag(Ngrids) )
      allocate ( step_counter(Ngrids) )
      allocate ( tdays(Ngrids) )
      allocate ( time(Ngrids) )
      allocate ( time_code(Ngrids) )
!$OMP END PARALLEL

      allocate ( NrecFrc(Ngrids) )

      allocate ( NudgingCoeff(Ngrids) )

      allocate ( ObcData(Ngrids) )

      allocate ( Lbiology(Ngrids) )
      allocate ( Lfloats(Ngrids) )
#ifdef ICE_MODEL
      allocate ( Lice(Ngrids) )
#endif
      allocate ( Lsediment(Ngrids) )
      allocate ( Lstations(Ngrids) )

      allocate ( CompositeGrid(4,Ngrids) )
      allocate ( RefinedGrid(Ngrids) )
      allocate ( RefineScale(Ngrids) )
      allocate ( GetDonorData(Ngrids) )

      allocate ( EWperiodic(Ngrids) )
      allocate ( NSperiodic(Ngrids) )

      allocate ( VolCons(4,Ngrids) )
#if defined ADJOINT || defined TANGENT || defined TL_IOMS
      allocate ( ad_VolCons(4,Ngrids) )
      allocate ( tl_VolCons(4,Ngrids) )
#endif

      allocate ( Lsponge(Ngrids) )
      allocate ( LuvSponge(Ngrids) )
      allocate ( LtracerSponge(MT,Ngrids) )

      allocate ( CLM_FILE(Ngrids) )
      allocate ( Lclimatology(Ngrids) )
      allocate ( LsshCLM(Ngrids) )
      allocate ( Lm2CLM(Ngrids) )
      allocate ( Lm3CLM(Ngrids) )
      allocate ( LtracerCLM(MT,Ngrids) )

      allocate ( Lnudging(Ngrids) )
      allocate ( LnudgeM2CLM(Ngrids) )
      allocate ( LnudgeM3CLM(Ngrids) )
      allocate ( LnudgeTCLM(MT,Ngrids) )

      allocate ( LuvSrc(Ngrids) )
      allocate ( LwSrc(Ngrids) )
      allocate ( LtracerSrc(MT,Ngrids) )

      allocate ( maxspeed(Ngrids) )
      allocate ( maxrho(Ngrids) )

      allocate ( levsfrc(Ngrids) )
      allocate ( levbfrc(Ngrids) )

      allocate ( Vtransform(Ngrids) )
      allocate ( Vstretching(Ngrids) )

      allocate ( Tcline(Ngrids) )
      allocate ( hc(Ngrids) )
      allocate ( theta_s(Ngrids) )
      allocate ( theta_b(Ngrids) )

      allocate ( hmin(Ngrids) )
      allocate ( hmax(Ngrids) )

      allocate ( xl(Ngrids) )
      allocate ( el(Ngrids) )

      allocate ( LonMin(Ngrids) )
      allocate ( LonMax(Ngrids) )
      allocate ( LatMin(Ngrids) )
      allocate ( LatMax(Ngrids) )

      allocate ( Idigits(Ngrids) )
      allocate ( Jdigits(Ngrids) )
#ifdef SOLVE3D
      allocate ( Kdigits(Ngrids) )
#endif

      allocate ( grdmax(Ngrids) )

#ifdef DIFF_3DCOEF
      allocate ( DiffMin(Ngrids) )
      allocate ( DiffMax(Ngrids) )
#endif

#ifdef VISC_3DCOEF
      allocate ( ViscMin(Ngrids) )
      allocate ( ViscMax(Ngrids) )
#endif

      allocate ( R0(Ngrids) )
      allocate ( Tcoef(Ngrids) )
      allocate ( Scoef(Ngrids) )

      allocate ( T0(Ngrids) )
      allocate ( S0(Ngrids) )

      allocate ( gamma2(Ngrids) )

      allocate ( lmd_Jwt(Ngrids) )

      allocate ( rdrg(Ngrids) )
      allocate ( rdrg2(Ngrids) )

      allocate ( Zos(Ngrids) )
      allocate ( Zob(Ngrids) )

      allocate ( Dcrit(Ngrids) )

#ifdef PROPAGATOR
      allocate ( Nconv(Ngrids) )
#endif

      allocate ( weight(2,0:256,Ngrids) )

      allocate ( Akk_bak(Ngrids) )
      allocate ( Akp_bak(Ngrids) )
      allocate ( Akv_bak(Ngrids) )

      allocate ( ad_visc2(Ngrids) )
      allocate ( nl_visc2(Ngrids) )
      allocate ( tl_visc2(Ngrids) )
      allocate ( visc2(Ngrids) )

      allocate ( ad_visc4(Ngrids) )
      allocate ( nl_visc4(Ngrids) )
      allocate ( tl_visc4(Ngrids) )
      allocate ( visc4(Ngrids) )

      allocate ( tkenu2(Ngrids) )
      allocate ( tkenu4(Ngrids) )

      allocate ( Akt_bak(MT,Ngrids) )
      allocate ( Kdiff(MT,Ngrids) )

      allocate ( ad_tnu2(MT,Ngrids) )
      allocate ( nl_tnu2(MT,Ngrids) )
      allocate ( tl_tnu2(MT,Ngrids) )
      allocate ( tnu2(MT,Ngrids) )

      allocate ( ad_tnu4(MT,Ngrids) )
      allocate ( nl_tnu4(MT,Ngrids) )
      allocate ( tl_tnu4(MT,Ngrids) )
      allocate ( tnu4(MT,Ngrids) )

      allocate ( tl_M2diff(Ngrids) )
      allocate ( tl_M3diff(Ngrids) )
      allocate ( tl_Tdiff(MT,Ngrids) )

      allocate ( ad_Akv_fac(Ngrids) )
      allocate ( tl_Akv_fac(Ngrids) )

      allocate ( ad_Akt_fac(MT,Ngrids) )
      allocate ( tl_Akt_fac(MT,Ngrids) )

      allocate ( LdefADJ(Ngrids) )
      allocate ( LdefAVG(Ngrids) )
      allocate ( LdefDAI(Ngrids) )
      allocate ( LdefAVG2(Ngrids) )
      allocate ( LdefDIA(Ngrids) )
      allocate ( LdefERR(Ngrids) )
      allocate ( LdefFLT(Ngrids) )
      allocate ( LdefHIS(Ngrids) )
      allocate ( LdefHIS2(Ngrids) )
      allocate ( LdefHSS(Ngrids) )
      allocate ( LdefINI(Ngrids) )
      allocate ( LdefIRP(Ngrids) )
      allocate ( LdefITL(Ngrids) )
      allocate ( LdefLCZ(Ngrids) )
      allocate ( LdefLZE(Ngrids) )
      allocate ( LdefMOD(Ngrids) )
      allocate ( LdefRST(Ngrids) )
      allocate ( LdefSTA(Ngrids) )
      allocate ( LdefTIDE(Ngrids) )
      allocate ( LdefTLM(Ngrids) )
      allocate ( LdefTLF(Ngrids) )

      allocate ( LwrtADJ(Ngrids) )
      allocate ( LwrtAVG(Ngrids) )
      allocate ( LwrtAVG2(Ngrids) )
      allocate ( LwrtDIA(Ngrids) )
      allocate ( LwrtHIS(Ngrids) )
      allocate ( LwrtHIS2(Ngrids) )
      allocate ( LwrtPER(Ngrids) )
      allocate ( LwrtRST(Ngrids) )
      allocate ( LwrtTLM(Ngrids) )
      allocate ( LwrtTLF(Ngrids) )

      allocate ( LdefNRM(4,Ngrids) )
      allocate ( LwrtNRM(4,Ngrids) )

#if defined STOCHASTIC_OPT && !defined STOCH_OPT_WHITE
      allocate ( LwrtState3d(Ngrids) )
      allocate ( SOinitial(Ngrids) )
#endif

      allocate ( LwrtState2d(Ngrids) )
      allocate ( LwrtTime(Ngrids) )

      allocate ( Ladjusted(Ngrids) )

      allocate ( LwrtInfo(Ngrids) )

      allocate ( ldefout(Ngrids) )

      allocate ( ndefADJ(Ngrids) )
      allocate ( ndefAVG(Ngrids) )
      allocate ( ndefAVG2(Ngrids) )
      allocate ( ndefDIA(Ngrids) )
      allocate ( ndefHIS(Ngrids) )
      allocate ( ndefHIS2(Ngrids) )
      allocate ( ndefTLM(Ngrids) )

      allocate ( ntsAVG(Ngrids) )
      allocate ( ntsAVG2(Ngrids) )
      allocate ( ntsDIA(Ngrids) )

      allocate ( nADJ(Ngrids) )
      allocate ( nAVG(Ngrids) )
      allocate ( nAVG2(Ngrids) )
      allocate ( nDIA(Ngrids) )
      allocate ( nFLT(Ngrids) )
      allocate ( nHIS(Ngrids) )
      allocate ( nHIS2(Ngrids) )
      allocate ( nRST(Ngrids) )
      allocate ( nSTA(Ngrids) )
      allocate ( nTLM(Ngrids) )

      allocate ( ninfo(Ngrids) )

      allocate ( nOBC(Ngrids) )
      allocate ( Nbrec(Ngrids) )
      allocate ( OBCcount(Ngrids) )

#ifdef ADJUST_BOUNDARY
      allocate ( Lobc(4,MstateVar,Ngrids) )
#endif

      allocate ( nSFF(Ngrids) )
      allocate ( Nfrec(Ngrids) )
      allocate ( SFcount(Ngrids) )

#ifdef ADJUST_STFLUX
      allocate ( Lstflux(MT,Ngrids) )
#endif

      allocate ( nrrec(Ngrids) )

      allocate ( LastRec(Ngrids) )

      allocate ( LcycleADJ(Ngrids) )
      allocate ( LcycleRST(Ngrids) )
      allocate ( LcycleTLM(Ngrids) )

#if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
      allocate ( Hcount(Ngrids) )
#endif

      allocate ( KstrS(Ngrids) )
      allocate ( KendS(Ngrids) )

      allocate ( DstrS(Ngrids) )
      allocate ( DendS(Ngrids) )

      allocate ( TRnorm(Ngrids) )

      allocate ( SO_decay(Ngrids) )
#ifdef SOLVE3D
      allocate ( SO_sdev(7+2*MT,Ngrids) )
#else
      allocate ( SO_sdev(5,Ngrids) )
#endif

#if defined FOUR_DVAR || defined VERIFICATION
      allocate ( Cnorm(2,MstateVar) )
      allocate ( CnormB(MstateVar,4) )

      allocate ( SporadicImpulse(Ngrids) )
      allocate ( FrequentImpulse(Ngrids) )

      allocate ( dTdz_min(Ngrids) )
      allocate ( ml_depth(Ngrids) )

      allocate ( LNM_depth(Ngrids) )

      allocate ( balance(MstateVar) )

      allocate ( Hdecay(2,MstateVar,Ngrids) )
      allocate ( Vdecay(2,MstateVar,Ngrids) )

      allocate ( Tdecay(MstateVar,Ngrids) )

      allocate ( HdecayB(MstateVar,4,Ngrids) )
      allocate ( VdecayB(MstateVar,4,Ngrids) )

      allocate ( bgqc_type(Ngrids) )
      allocate ( Nprovenance(Ngrids) )
      allocate ( S_bgqc(MstateVar,Ngrids) )
#endif

      allocate ( obcfac(Ngrids) )
      allocate ( FSobc_in(Ngrids,4) )
      allocate ( FSobc_out(Ngrids,4) )
      allocate ( M2obc_in(Ngrids,4) )
      allocate ( M2obc_out(Ngrids,4) )

#ifdef SOLVE3D
      allocate ( M3obc_in(Ngrids,4) )
      allocate ( M3obc_out(Ngrids,4) )
      allocate ( Tobc_in(MT,Ngrids,4) )
      allocate ( Tobc_out(MT,Ngrids,4) )
#endif

      allocate ( Znudg(Ngrids) )
      allocate ( M2nudg(Ngrids) )
      allocate ( M3nudg(Ngrids) )
      allocate ( Tnudg(MT,Ngrids) )
      allocate ( Tnudg_SSS(Ngrids) )
      allocate ( sss_mismatch_threshold(Ngrids) )

#if defined BULK_FLUXES || defined BULK_FLUXES2D
      allocate ( blk_ZQ(Ngrids) )
      allocate ( blk_ZT(Ngrids) )
      allocate ( blk_ZW(Ngrids) )
#endif
#ifdef BIOLOGY
      allocate ( maxbio(NBT, Ngrids) )
      allocate ( max_bio(NBT) )
#endif

      allocate ( gls_m(Ngrids) )
      allocate ( gls_n(Ngrids) )
      allocate ( gls_p(Ngrids) )
      allocate ( gls_sigk(Ngrids) )
      allocate ( gls_sigp(Ngrids) )
      allocate ( gls_cmu0(Ngrids) )
      allocate ( gls_cmupr(Ngrids) )
      allocate ( gls_c1(Ngrids) )
      allocate ( gls_c2(Ngrids) )
      allocate ( gls_c3m(Ngrids) )
      allocate ( gls_c3p(Ngrids) )
      allocate ( gls_Kmin(Ngrids) )
      allocate ( gls_Pmin(Ngrids) )

      allocate ( charnok_alpha(Ngrids) )
      allocate ( zos_hsig_alpha(Ngrids) )
      allocate ( sz_alpha(Ngrids) )
      allocate ( crgban_cw(Ngrids) )
      allocate ( wec_alpha(Ngrids) )

#ifdef ICE_MODEL
     allocate ( nstrs(Ngrids) )
     allocate ( nevp(Ngrids) )
     allocate ( ievp(Ngrids) )

     allocate ( dtice(Ngrids) )
     allocate ( dte(Ngrids) )
     allocate ( rhoice(Ngrids) )
     allocate ( cdiw(Ngrids) )
     allocate ( cdai(Ngrids) )
     allocate ( rho_air(Ngrids) )
     allocate ( rhosnow_dry(Ngrids) )
     allocate ( rhosnow_wet(Ngrids) )
     allocate ( pstar(Ngrids) )
     allocate ( pstar_e(Ngrids) )
     allocate ( zetamax(Ngrids) )
     allocate ( zetamin(Ngrids) )
     allocate ( ellip_sq(Ngrids) )
     allocate ( astren(Ngrids) )
     allocate ( alphai(Ngrids) )
     allocate ( min_h(Ngrids) )
     allocate ( min_a(Ngrids) )
     allocate ( max_a(Ngrids) )
     allocate ( stressang(Ngrids) )
# ifdef ICE_LANDFAST
     allocate ( lf_k1(Ngrids) )
     allocate ( lf_k2(Ngrids) )
     allocate ( lf_u0(Ngrids) )
# endif
#endif

      RETURN
      END SUBROUTINE allocate_scalars

      SUBROUTINE initialize_scalars
!
!=======================================================================
!                                                                      !
!  This routine initializes several variables in module for all nested !
!  grids.                                               !              !
!                                                                      !
!=======================================================================
!
!  Local variable declarations.
!
      integer :: i, ic, j, ng, itrc

      real(r8) :: one, zero
      real(r8), parameter :: IniVal = 0.0_r8
!
!---------------------------------------------------------------------
!  Set tracer identification indices.
!---------------------------------------------------------------------
!
      itemp=1
      isalt=2
      ic=NAT

#ifdef T_PASSIVE
!
!  Indices for inert passive tracers to advect and diffuse.
!
      DO i=1,NPT
        ic=ic+1
        inert(i)=ic
      END DO
#endif

#ifdef DIAGNOSTICS
!
!---------------------------------------------------------------------
!  Set diagnostic fields identification indices.
!---------------------------------------------------------------------

# ifdef DIAGNOSTICS_TS
!
!  Indices for tracer diagnostic variables.
!
      iThadv=1
      iTxadv=2
      iTyadv=3
      iTvadv=4
      ic=4
#  if defined TS_DIF2 || defined TS_DIF4
      iThdif=ic+1
      iTxdif=ic+2
      iTydif=ic+3
      ic=ic+3
#   if defined MIX_GEO_TS || defined MIX_ISO_TS
      iTsdif=ic+1
      ic=ic+1
#   endif
#  endif
      iTvdif=ic+1
      iTrate=ic+2
# endif
# ifdef DIAGNOSTICS_UV
!
!  Indices for 2D momentum diagnostic variables.  In some places in
!  the code a compact DO-loop (idiag=1:M2pgrd) is used to improve
!  flexibility. Therefore, the order of indices is very important.
!  Only those fields that require special treatment are set below
!  the M2pgrd index.
!
      ic=0
#  if defined UV_COR
      M2fcor=ic+1
      ic=ic+1
#  endif
#  if defined UV_ADV
      M2hadv=ic+1
      M2xadv=ic+2
      M2yadv=ic+3
      ic=ic+3
#  endif
#  if defined WEC_MELLOR
      M2hrad=ic+1
      ic=ic+1
#  endif
#  if defined WEC_VF
#   if defined UV_COR
      M2fsco=ic+1
      ic=ic+1
#   endif
#   ifdef BOTTOM_STREAMING
      M2bstm=ic+1
      ic=ic+1
#   endif
#   ifdef SURFACE_STREAMING
      M2sstm=ic+1
      ic=ic+1
#   endif
      M2hjvf=ic+1
      M2kvrf=ic+2
      M2wrol=ic+3
      M2wbrk=ic+4
      ic=ic+4
#  endif
#  if defined VEGETATION && defined VEG_DRAG
      M2fveg=ic+1
      ic=ic+1
#  endif
#  if defined UV_VIS2 || defined UV_VIS4
      M2hvis=ic+1
      M2xvis=ic+2
      M2yvis=ic+3
      ic=ic+3
#  endif
      M2pgrd=ic+1
      M2sstr=ic+2             ! These indices need to be
      M2bstr=ic+3             ! specified last to allow a
#  if defined WEC_VF
      M2zeta=ic+4
      M2zetw=ic+5
      M2zqsp=ic+6
      M2zbeh=ic+7
#  endif
      M2rate=NDM2d            ! compact DO-loop structure
#  ifdef SOLVE3D
!
!  Indices for 3D momentum diagnostic variables.  In some places in
!  the code a compact DO-loop (idiag=1:M3pgrd) is used to improve
!  flexibility. Therefore, the order of indices is very important.
!  Only those fields that require special treatment are set below
!  the M3pgrd index.
!
      ic=0
#   if defined UV_COR
      M3fcor=ic+1
      ic=ic+1
#   endif
#   if defined UV_ADV
      M3vadv=ic+1
      M3hadv=ic+2
      M3xadv=ic+3
      M3yadv=ic+4
      ic=ic+4
#   endif
#   if defined WEC_MELLOR
      M3hrad=ic+1
      M3vrad=ic+2
      ic=ic+2
#   endif
#   if defined WEC_VF
#    if defined UV_COR
      M3fsco=ic+1
      ic=ic+1
#    endif
#    ifdef BOTTOM_STREAMING
      M3bstm=ic+1
      ic=ic+1
#    endif
#    ifdef SURFACE_STREAMING
      M3sstm=ic+1
      ic=ic+1
#    endif
      M3vjvf=ic+1
      M3hjvf=ic+2
      M3kvrf=ic+3
      M3wrol=ic+4
      M3wbrk=ic+5
      ic=ic+5
#   endif
#   if defined VEGETATION && defined VEG_DRAG
      M3fveg=ic+1
      ic=ic+1
#   endif
      M3pgrd=ic+1             ! needs to be here, indices below
      M3vvis=ic+2             ! require special treatment
#   if defined UV_VIS2 || defined UV_VIS4
      M3hvis=ic+3
      M3xvis=ic+4
      M3yvis=ic+5
#   endif
      M3rate=NDM3d
#  endif
# endif
#endif
!
!-----------------------------------------------------------------------
!  Activate all computation control switches.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
        LastRec(ng)=.FALSE.
        CompositeGrid(1:4,ng)=.FALSE.
        RefinedGrid(ng)=.FALSE.
        GetDonorData(ng)=.FALSE.
        Lbiology(ng)=.TRUE.
        LcycleADJ(ng)=.FALSE.
        LcycleRST(ng)=.FALSE.
        LcycleTLM(ng)=.FALSE.
        Lfloats(ng)=.TRUE.
#ifdef ICE_MODEL
        Lice(ng)=.TRUE.
#endif
        Lsediment(ng)=.TRUE.
        Lstations(ng)=.TRUE.
#ifdef TIDES_ASTRO
        FIRST_TIDES_ASTRO(ng) = .TRUE.
#endif
#if defined FOUR_DVAR || defined VERIFICATION
        FrequentImpulse(ng)=.FALSE.
        SporadicImpulse(ng)=.FALSE.
#endif
      END DO

#ifdef ADJUST_BOUNDARY
!
!-----------------------------------------------------------------------
!  Initilize switches to process open boundary arrays for 4DVar
!  adjustments.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
        DO j=1,MstateVar
          DO i=1,4
            Lobc(i,j,ng)=.FALSE.
          END DO
        END DO
      END DO
#endif
#ifdef ADJUST_STFLUX
!
!-----------------------------------------------------------------------
!  Initilize switches to process surface tracer fluexes for 4DVar
!  adjustments.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
        DO i=1,MT
          Lstflux(i,ng)=.FALSE.
        END DO
      END DO
#endif
!
!-----------------------------------------------------------------------
!  Initialize several scalar variables.
!-----------------------------------------------------------------------
!
      one=1.0_r8
      zero=0.0_r8
#ifdef SWAN_COUPLING
      Infinity=one/zero
#endif
      Co=1.0_r8/(2.0_r8+SQRT(2.0_r8))
      gorho0=g/rho0
      DO ng=1,Ngrids
        EWperiodic(ng)=.FALSE.
        NSperiodic(ng)=.FALSE.
        NudgingCoeff(ng)=.FALSE.
        ObcData(ng)=.FALSE.
        RefineScale(ng)=0
        gamma2(ng)=-1.0_r8
        Vtransform(ng)=1
        Vstretching(ng)=1
#if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
        Hcount(ng)=0
#endif
#ifdef ADJUST_BOUNDARY
        OBCcount(ng)=0
#endif
#if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
        SFcount(ng)=0
#endif
        first_time(ng)=0
        Idigits(ng)=INT(LOG10(REAL(Lm(ng),r8)))+1
        Jdigits(ng)=INT(LOG10(REAL(Mm(ng),r8)))+1
#ifdef SOLVE3D
        Kdigits(ng)=INT(LOG10(REAL(N (ng),r8)))+1
#endif
        CLM_FILE(ng)=.FALSE.
        Lnudging(ng)=.FALSE.
        LnudgeM2CLM(ng)=.FALSE.
        LnudgeM3CLM(ng)=.FALSE.
        Lclimatology(ng)=.FALSE.
        Lm2CLM(ng)=.FALSE.
        Lm3CLM(ng)=.FALSE.
        LsshCLM(ng)=.FALSE.
        Lsponge(ng)=.FALSE.
        LuvSponge(ng)=.FALSE.
        LuvSrc(ng)=.FALSE.
        LwSrc(ng)=.FALSE.
        DO itrc=1,MT
          LnudgeTCLM(itrc,ng)=.FALSE.
          LtracerCLM(itrc,ng)=.FALSE.
          LtracerSrc(itrc,ng)=.FALSE.
          LtracerSponge(itrc,ng)=.FALSE.
          ad_Akt_fac(itrc,ng)=1.0_r8
          tl_Akt_fac(itrc,ng)=1.0_r8
          ad_tnu2(itrc,ng)=IniVal
          nl_tnu2(itrc,ng)=IniVal
          tl_tnu2(itrc,ng)=IniVal
          tnu2(itrc,ng)=IniVal
          ad_tnu4(itrc,ng)=IniVal
          nl_tnu4(itrc,ng)=IniVal
          tl_tnu4(itrc,ng)=IniVal
          tnu4(itrc,ng)=IniVal
        END DO
        ad_Akv_fac(ng)=1.0_r8
        tl_Akv_fac(ng)=1.0_r8
        ad_visc2(ng)=IniVal
        nl_visc2(ng)=IniVal
        tl_visc2(ng)=IniVal
        visc2(ng)=IniVal
        ad_visc4(ng)=IniVal
        nl_visc4(ng)=IniVal
        tl_visc4(ng)=IniVal
        visc4(ng)=IniVal
#if defined BULK_FLUXES || defined BULK_FLUXES2D
        blk_ZQ(ng)=10.0_r8
        blk_ZT(ng)=10.0_r8
        blk_ZW(ng)=10.0_r8
#endif
        DO i=1,4
          VolCons(i,ng)=.FALSE.
#if defined ADJOINT || defined TANGENT || defined TL_IOMS
          ad_VolCons(i,ng)=.FALSE.
          tl_VolCons(i,ng)=.FALSE.
#endif
          FSobc_in (ng,i)=IniVal
          FSobc_out(ng,i)=IniVal
          M2obc_in (ng,i)=IniVal
          M2obc_out(ng,i)=IniVal
#ifdef SOLVE3D
          M3obc_in (ng,i)=IniVal
          M3obc_out(ng,i)=IniVal
#endif
        END DO
      END DO
#ifdef SOLVE3D
      Tobc_in = IniVal
      Tobc_out = IniVal
#endif
!
!  Initialize thread private variables.
!
!$OMP PARALLEL
      synchro_flag=.FALSE.

      ntfirst=1
      ntstart=1
      ntend=0
      step_counter=0
!$OMP END PARALLEL

#if defined LMD_SKPP || defined LMD_BKPP
!
!  Proportionality coefficient parameterizing boundary layer
!  nonlocal transport.
!
      lmd_Cg=lmd_Cstar*                                                 &
     &       vonKar*(lmd_cs*vonKar*lmd_epsilon)**(1.0_r8/3.0_r8)
#endif

#if defined FOUR_DVAR || defined VERIFICATION
!
!  Initialize error covariace variables.
!
      balance=.FALSE.

      Cnorm=.FALSE.
      Hdecay=IniVal
      Tdecay=IniVal
      Vdecay=IniVal

      CnormB=.FALSE.
      HdecayB=IniVal
      VdecayB=IniVal
#endif
!
!  Initialize several IO flags.
!
      LmultiGST=.FALSE.
      LrstGST=.FALSE.

      DO ng=1,Ngrids

        PerfectRST(ng)=.FALSE.

        Ladjusted(ng)=.FALSE.

        LdefADJ(ng)=.FALSE.
        LdefAVG(ng)=.TRUE.
        LdefDAI(ng)=.FALSE.
        LdefAVG2(ng)=.TRUE.
        LdefDIA(ng)=.TRUE.
        LdefERR(ng)=.FALSE.
#ifdef FLOATS
        LdefFLT(ng)=.TRUE.
#endif
        LdefHIS(ng)=.TRUE.
        LdefHIS2(ng)=.TRUE.
        LdefINI(ng)=.FALSE.
        LdefIRP(ng)=.FALSE.
        LdefITL(ng)=.FALSE.
        LdefMOD(ng)=.FALSE.
        LdefRST(ng)=.TRUE.
        LdefSTA(ng)=.TRUE.
        LdefTLM(ng)=.FALSE.
#if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
        LdefTIDE(ng)=.TRUE.
#endif
        LwrtADJ(ng)=.FALSE.
        LwrtAVG(ng)=.FALSE.
        LwrtAVG2(ng)=.FALSE.
        LwrtDIA(ng)=.FALSE.
        LwrtHIS(ng)=.FALSE.
        LwrtHIS2(ng)=.FALSE.
        LwrtPER(ng)=.FALSE.
        LwrtRST(ng)=.FALSE.
        LwrtTLM(ng)=.FALSE.
        LwrtInfo(ng)=.TRUE.
        LwrtState2d(ng)=.FALSE.
#if defined STOCHASTIC_OPT && !defined STOCH_OPT_WHITE
        LwrtState3d(ng)=.FALSE.
        SOinitial(ng)=.FALSE.
#endif
        LwrtTime(ng)=.TRUE.
        ldefout(ng)=.FALSE.

      END DO
#ifdef BIOLOGY
      DO ng=1,Ngrids
        DO itrc=1,NBT
          maxbio(itrc, ng) = IniVal
        END DO
      END DO
      DO itrc=1,NBT
        max_bio(itrc) = 10000.0_r8   ! bio units
      END DO
#endif

# if defined SG_BBL || defined SSW_BBL
!
!  Nondimensional closure parameters associated with Styles and Glenn
!  (1999) bottom currents and waves boundary layer.
!
      sg_z1p=sg_alpha
      sg_mp=CMPLX(SQRT(1.0_r8/(2.0_r8*sg_z1p)),                         &
     &            SQRT(1.0_r8/(2.0_r8*sg_z1p)))
#endif
#if defined GLS_MIXING || defined MY25_MIXING
# if defined CANUTO_A || defined CANUTO_B
!
! Compute parameters for Canuto et al. (2001) stability functions.
!   (Canuto, V.M., Cheng, H.Y., and Dubovikov, M.S., 2001: Ocean
!    turbulence. Part I: One-point closure model - momentum and
!    heat vertical diffusivities, JPO, 1413-1426).
!
      gls_s0=3.0_r8/2.0_r8*gls_L1*gls_L5**2
      gls_s1=-gls_L4*(gls_L6+gls_L7)+2.0_r8*gls_L4*gls_L5*              &
     &        (gls_L1-1.0_r8/3.0_r8*gls_L2-gls_L3)+3.0_r8/2.0_r8*       &
     &        gls_L1*gls_L5*gls_L8
      gls_s2=-3.0_r8/8.0_r8*gls_L1*(gls_L6**2-gls_L7**2)
      gls_s4=2.0_r8*gls_L5
      gls_s5=2.0_r8*gls_L4
      gls_s6=2.0_r8/3.0_r8*gls_L5*(3.0_r8*gls_L3**2-gls_L2**2)-         &
     &       1.0_r8/2.0_r8*gls_L5*gls_L1*(3.0_r8*gls_L3-gls_L2)+        &
     &       3.0_r8/4.0_r8*gls_L1*(gls_L6-gls_L7)
      gls_b0=3.0_r8*gls_L5**2
      gls_b1=gls_L5*(7.0_r8*gls_L4+3.0_r8*gls_L8)
      gls_b2=gls_L5**2*(3.0_r8*gls_L3**2-gls_L2**2)-                    &
     &       3.0_r8/4.0_r8*(gls_L6**2-gls_L7**2)
      gls_b3=gls_L4*(4.0_r8*gls_L4+3.0_r8*gls_L8)
      gls_b5=1.0_r8/4.0_r8*(gls_L2**2-3.0_r8*gls_L3**2)*                &
     &       (gls_L6**2-gls_L7**2)
      gls_b4=gls_L4*(gls_L2*gls_L6-3.0_r8*gls_L3*gls_L7-                &
     &       gls_L5*(gls_L2**2-gls_L3**2))+gls_L5*gls_L8*               &
     &       (3.0_r8*gls_L3**2-gls_L2**2)
# endif
!
!  Coefficients used to compute stability functions for tracer and
!  momentum.
!
      my_B1p2o3=my_B1**(2.0_r8/3.0_r8)
      my_B1pm1o3=1.0_r8/(my_B1**(1.0_r8/3.0_r8))
      my_E1o2=0.5_r8*my_E1
      my_Sm1=my_A1*my_A2*((my_B2-3.0_r8*my_A2)*                         &
     &                    (1.0_r8-6.0_r8*my_A1/my_B1)-                  &
     &                    3.0_r8*my_C1*(my_B2+6.0_r8*my_A1))
      my_Sm2=9.0_r8*my_A1*my_A2
      my_Sh1=my_A2*(1.0_r8-6.0_r8*my_A1/my_B1)
# ifdef KANTHA_CLAYSON
      my_Sh2=3.0_r8*my_A2*(6.0_r8*my_A1+my_B2*(1.0_r8-my_C3))
      my_Sm4=18.0_r8*my_A1*my_A1+9.0_r8*my_A1*my_A2*(1.0_r8-my_C2)
# else
      my_Sh2=3.0_r8*my_A2*(6.0_r8*my_A1+my_B2)
      my_Sm3=my_A1*(1.0_r8-3.0_r8*my_C1-6.0_r8*my_A1/my_B1)
      my_Sm4=18.0_r8*my_A1*my_A1+9.0_r8*my_A1*my_A2
# endif
#endif
#ifdef CICE_MODEL
      rhoice = 900.0_r8
      min_a = 1.0e-11_r8
#endif

      RETURN
      END SUBROUTINE initialize_scalars

      END MODULE mod_scalars
