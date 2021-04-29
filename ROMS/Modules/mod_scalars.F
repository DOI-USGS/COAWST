#include "cppdefs.h"
      MODULE mod_scalars
!
!svn $Id: mod_scalars.F 1054 2021-03-06 19:47:12Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2021 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!
        USE mod_param
#ifdef NO_4BYTE_REALS
        USE netcdf
#endif
!        implicit none
!
!-----------------------------------------------------------------------
!  Multiple grid structure.
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
          real(dp), pointer :: Cs_r(:)
          real(dp), pointer :: Cs_w(:)
          real(dp), pointer :: sc_r(:)
          real(dp), pointer :: sc_w(:)
        END TYPE T_SCALARS
!
        TYPE (T_SCALARS), allocatable :: SCALARS(:)
!
!-----------------------------------------------------------------------
!  Time clock structure.
!-----------------------------------------------------------------------
!
!  Reference time (yyyymmdd.f) used to compute relative time. The
!  application date clock is measured ad elapsed time interval since
!  reference-time. This parameter also provides information about the
!  calendar used:
!
!    If TIME_REF = -2, the model time and DSTART are in modified Julian
!                  days units.  The time "units" attribute is:
!
!                  'time-units since 1968-05-23 00:00:00 GMT'
!
!    If TIME_REF = -1, the model time and DSTART are in a calendar
!                  with 360 days in every year (30 days each month).
!                  The time "units" attribute is:
!
!                  'time-units since 0001-01-01 00:00:00'
!
!    If TIME_REF = 0, the model time and DSTART are in a common year
!                  calendar with 365.2524 days.  The "units" attribute
!                  is:
!
!                  'time-units since 0001-01-01 00:00:00'
!
!    If TIME_REF > 0, the model time and DSTART are the elapsed time
!                  units since specified reference time.  For example,
!                  TIME_REF=20020115.5 will yield the following
!                  time "units" attribute:
!
!                  'time-units since 2002-01-15 12:00:00'
!
        real(dp) :: time_ref = 0.0_dp                    ! YYYYMMDD.dd
!
      TYPE T_CLOCK
        integer :: yday                ! day of the year
        integer :: year                ! year including century (YYYY)
        integer :: month               ! month of the year (1,...,12)
        integer :: day                 ! day of the month
        integer :: hour                ! hour of the day (1,...,23)
        integer :: minutes             ! minutes of the hour

        real(dp) :: seconds            ! frational seconds of the minute
        real(dp) :: base               ! reference date (YYYYMMDD.dd)
        real(dp) :: DateNumber(2)      ! date number, [1]: days
                                       !              [2]: seconds

        character (len=22) :: string   ! YYYY-MM-DD hh:mm:ss.ss
        character (len=25) :: calendar ! date calendar
      END TYPE T_CLOCK
!
        TYPE (T_CLOCK) :: Rclock       ! reference/base date
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
!    run_time      Total run time for all nested grids (s), it is
!                    set in Masters/ocean.h
!    MyRunInterval Total run time for all nested grids (s), it is
!                    set in Drivers/nl_ocean.h (coupling window)
!    io_time       Current I/O time (s) processed in "get_state".
!    tdays         Model time clock (days).
!    time          Model time clock (s).
!    time_code     Model time clock (string, YYYY-MM-DD hh:mm:ss.ss)
!    AVGtime       Model time clock for averages output (s).
!    AVG2time      Model time clock for averages output (s).
!    DIAtime       Model time clock for diagnostics output (s).
!    F_code        Final time string for simulation
!    I_code        Initial time string for simulation
!    INItime       Nonlinear model initial conditions time (s).
!    IMPtime       Impulse forcing time (s) to process.
!    ObsTime       Observation time (s) to process.
!    FrcTime       Adjoint or tangent linear Impulse forcing time (s).
!    dstart        Time stamp assigned to model initialization (usually
!                    a Calendar day, like modified Julian Day).
!    tide_start    Reference time for tidal forcing (days).
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

        real(dp), allocatable :: tdays(:)                ! days
        real(dp), allocatable :: time(:)                 ! seconds
!$OMP THREADPRIVATE (tdays, time)

        real(dp), allocatable :: dt(:)                   ! seconds
        real(dp), allocatable :: dtfast(:)               ! seconds

        real(dp), allocatable :: TimeEnd(:)              ! seconds
        real(dp), allocatable :: AVGtime(:)              ! seconds
        real(dp), allocatable :: DIAtime(:)              ! seconds
        real(dp), allocatable :: IMPtime(:)              ! seconds
        real(dp), allocatable :: INItime(:)              ! seconds
        real(dp), allocatable :: INItimeS(:)             ! seconds
        real(dp), allocatable :: ObsTime(:)              ! seconds
        real(dp), allocatable :: FrcTime(:)              ! seconds
#ifdef TIDES_ASTRO
        logical, allocatable :: FIRST_TIDES_ASTRO(:)
#endif
#ifdef AGE_DISTRIBUTION
        real(r8), allocatable :: dtau(:)              ! seconds
#endif

        real(dp) :: dstart = 0.0_dp                      ! days
#ifdef SP4DVAR
        real(dp) :: dstartS = 0.0_dp                     ! days
#endif
        real(dp) :: io_time = 0.0_dp                     ! seconds
        real(dp) :: run_time = 0.0_dp                    ! seconds
        real(dp) :: MyRunInterval = 0.0_dp               ! seconds
        real(dp) :: tide_start = 0.0_dp                  ! days

        character (len=22) :: F_code, I_code

        character (len=22), allocatable :: time_code(:)  ! date string
!$OMP THREADPRIVATE (time_code)

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
        real(dp) :: Falpha = 2.0_dp
        real(dp) :: Fbeta  = 4.0_dp
        real(dp) :: Fgamma = 0.284_dp
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
#ifdef SENSITIVITY_4DVAR
        integer :: NrunSAVE = 0                 ! Loop counter
#endif
        integer :: OuterLoop = 0                ! split outer loop
        integer :: inner = 0                    ! inner loop counter
        integer :: outer = 0                    ! outer loop counter
!
! Set checksum algorithm: "adler32" or "crc32"
!
        character (len=*), parameter :: HashMethod = "crc32"
!
! Split 4D-Var phase.
!
        character (len=20) :: Phase4DVAR
!
! Number of sadde point 4D-Var intervals.
!
        integer :: Nsaddle = 1
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
!  HSIMT tracer advection coefficients for the TVD limiter (Wu and Zhu,
!  2010).
!
        real(r8) :: cc1 = 0.25_r8
        real(r8) :: cc2 = 0.5_r8
        real(r8) :: cc3 = 1.0_r8/12.0_r8
!
!-----------------------------------------------------------------------
!  Control switches.
!-----------------------------------------------------------------------
!
!  Switch to use three ghost-points in the halo region.
!
        logical :: ThreeGhostPoints = .FALSE.
!
!  Switch to set-up application grid, metrics, and associated variables
!  and parameters.
!
        logical, allocatable :: SetGridConfig(:)
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
!    emmiss        Infrared emissivity.
!    g             Acceleration due to gravity (m/s2).
!    gorho0        gravity divided by mean density anomaly.
!    rhow          fresh water density (kg/m3).
!    vonKar        von Karman constant.
!
#ifdef ISOMIP
        real(dp) :: Cp = 3974.0_dp              ! mod (4/24/06) for ICETEST
#else
        real(dp) :: Cp = 3985.0_dp              ! Joules/kg/degC
#endif
        real(dp) :: Csolar = 1353.0_dp          ! 1360-1380 W/m2
        real(dp) :: Infinity                    ! Infinity = 1.0/0.0
        real(dp) :: Eradius = 6371315.0_dp      ! m
        real(dp) :: StefBo = 5.67E-8_dp         ! Watts/m2/K4
        real(dp) :: emmiss = 0.97_dp            ! non_dimensional
        real(dp) :: rhow = 1000.0_dp            ! kg/m3
#ifdef SOLITON
        real(dp) :: g = 1.0_dp                  ! non-dimensional
#elif defined WBC_1 || defined WBC_2 || defined WBC_3
        real(dp) :: g = 9.8_dp                  ! m/s2
#elif defined CIRCLE
        real(dp) :: g = 3.92e-2_dp                  ! m/s2
#else
        real(dp) :: g = 9.81_dp                 ! m/s2
#endif
        real(dp) :: gorho0                      ! m4/s2/kg
        real(dp) :: vonKar = 0.41_dp            ! non-dimensional
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
!    exit_flag = 9   coupling error
!
        integer :: exit_flag = 0
        integer :: blowup = 0
        integer :: NoError = 0
!
!  Blow-up string.
!
        character (len=80) :: blowup_string
!
!  Set threshold maximum speed (m/s) and density anomaly (kg/m3) to
!  test if the model is blowing-up.
!
        real(dp), allocatable :: maxspeed(:)
        real(dp), allocatable :: maxrho(:)
!
        real(dp) :: max_speed = 20.0_dp         ! m/s
        real(dp) :: max_rho = 200.0_dp          ! kg/m3

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
        real(dp), allocatable :: Tcline(:)      ! m, positive
        real(dp), allocatable :: hc(:)          ! m, positive
        real(dp), allocatable :: theta_s(:)     ! 0 < theta_s < 20
        real(dp), allocatable :: theta_b(:)     ! 0 < theta_b < 1
!
!  Bathymetry range values.
!
        real(dp), allocatable :: hmin(:)        ! m, positive
        real(dp), allocatable :: hmax(:)        ! m, positive
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

        real(dp) :: avgke = 0.0_dp              ! Kinetic energy
        real(dp) :: avgpe = 0.0_dp              ! Potential energy
        real(dp) :: avgkp = 0.0_dp              ! Total energy
        real(dp) :: volume = 0.0_dp             ! diagnostics volume
        real(dp) :: ad_volume = 0.0_dp          ! adjoint volume

        real(dp), allocatable :: TotVolume(:)   ! Total volume
        real(dp), allocatable :: MinVolume(:)   ! Minimum cell volume
        real(dp), allocatable :: MaxVolume(:)   ! Maximum cell volume
!
!  Minimun and maximum grid spacing
!
        real(dp), allocatable :: DXmin(:)       ! all grid points
        real(dp), allocatable :: DXmax(:)
        real(dp), allocatable :: DYmin(:)
        real(dp), allocatable :: DYmax(:)
#ifdef MASKING
        real(dp), allocatable :: DXminW(:)      ! only water points
        real(dp), allocatable :: DXmaxW(:)
        real(dp), allocatable :: DYminW(:)
        real(dp), allocatable :: DYmaxW(:)
#endif
#ifdef SOLVE3D
        real(dp), allocatable :: DZmin(:)       ! all grid points
        real(dp), allocatable :: DZmax(:)
# ifdef MASKING
        real(dp), allocatable :: DZminW(:)      ! only water points
        real(dp), allocatable :: DZmaxW(:)
# endif
#endif
!
!  Maximum size of a grid node (m) over the whole curvilinear grid
!  application. Used for scaling horizontal mixing by the grid size.
!
        real(dp), allocatable :: grdmax(:)
#ifdef DIFF_3DCOEF
        real(dp), allocatable :: DiffMin(:)     ! Minimun diffusion
        real(dp), allocatable :: DiffMax(:)     ! Maximum diffusion
#endif
#ifdef VISC_3DCOEF
        real(dp), allocatable :: ViscMin(:)     ! Minimum viscosity
        real(dp), allocatable :: ViscMax(:)     ! Maximum viscosity
#endif
!
!  Courant Numbers due to gravity wave speed limits.
!
        real(dp), allocatable :: Cg_min(:)      ! Minimun barotropic
        real(dp), allocatable :: Cg_max(:)      ! Maximun barotropic
        real(dp), allocatable :: Cg_Cor(:)      ! Maximun Coriolis
!
!  Time dependent Counrant Numbers due to velocity components and
!  indices location of maximum value.
!
        integer :: max_Ci = 0                   ! maximum I-location
        integer :: max_Cj = 0                   ! maximum J-location
        integer :: max_Ck = 0                   ! maximum K-location
        real(r8) :: max_C = 0.0_r8              ! maximum total
        real(r8) :: max_Cu = 0.0_r8             ! maximum I-component
        real(r8) :: max_Cv = 0.0_r8             ! maximum J-component
#ifdef SOLVE3D
        real(r8) :: max_Cw = 0.0_r8             ! maximum K-component
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
        real(dp), allocatable :: rx0(:)         ! Beckmann and Haidvogel
        real(dp), allocatable :: rx1(:)         ! Haney
!
!  Linear (m/s) and quadratic (nondimensional) bottom drag coefficients.
!
        real(r8), allocatable :: rdrg(:)
        real(r8), allocatable :: rdrg2(:)
!
!  Minimum and maximum threshold for transfer coefficient of momentum.
!
        real(dp) :: Cdb_min = 0.000001_dp
        real(dp) :: Cdb_max = 0.5_dp
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
        real(dp) :: rho0 = 1025.0_dp
!
!  Background Brunt-Vaisala frequency (1/s2).
!
        real(dp) :: bvf_bak = 0.00001_dp

#ifdef PROPAGATOR
!
!  Number of converged Ritz values and relative accuracy of computed
!  Ritz values.
!
        integer, allocatable :: Nconv(:)
        real(dp) :: Ritz_tol = 1.0E-15_dp
#endif
!
!  Vector containing USER generic parameters.
!
        integer :: Nuser
        real(r8), dimension(25) :: user(25)
!
!  Weights for the time average of 2D fields.
!
        real(dp), allocatable :: weight(:,:,:)
!
!  Constants.
!
        real(dp), parameter :: pi = 3.14159265358979323846_dp
        real(dp), parameter :: deg2rad = pi / 180.0_dp
        real(dp), parameter :: rad2deg = 180.0_dp / pi
        real(dp), parameter :: day2sec = 86400.0_dp
        real(dp), parameter :: sec2day = 1.0_dp / 86400.0_dp
#ifdef NO_4BYTE_REALS
        real(dp), parameter :: spval = NF90_FILL_DOUBLE
#else
        real(dp), parameter :: spval = 1.0E+37_dp
#endif
        real(dp), parameter :: Large = 1.0E+20_dp
        real(dp), parameter :: jul_off = 2440000.0_dp
!
!  Set special check value.  Notice that a smaller value is assigned
!  to account for both NetCDF fill value and roundoff. There are
!  many Matlab scripts out there that do not inquire correctly
!  the spval from the _FillValue attribute in single/double
!  precision.
!
        real(dp), parameter :: spval_check = 1.0E+35_dp
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
!    Akt_limit     Upper threshold vertical mixing coefficient (m2/s)
!                    for tracers.
!    Akv_limit     Upper threshold vertical mixing coefficient (m2/s)
!                    for momentum.
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
        real(r8), allocatable :: Akv_limit(:)        ! m2/s

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
        real(r8), allocatable :: Akt_limit(:,:)      ! m2/s
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
        logical, allocatable :: LdefDIA(:)       ! Diagnostics file
        logical, allocatable :: LdefERR(:)       ! 4DVar error file
        logical, allocatable :: LdefFLT(:)       ! Floats file
        logical, allocatable :: LdefHIS(:)       ! History file
        logical, allocatable :: LdefHSS(:)       ! Hessian file
        logical, allocatable :: LdefINI(:)       ! Initial file
        logical, allocatable :: LdefIRP(:)       ! Initial RPM file
        logical, allocatable :: LdefITL(:)       ! Initial TLM file
        logical, allocatable :: LdefLCZ(:)       ! Lanczos Vectors file
        logical, allocatable :: LdefLZE(:)       ! Evolved Lanczos file
        logical, allocatable :: LdefMOD(:)       ! 4DVAR file
        logical, allocatable :: LdefQCK(:)       ! Quicksave file
        logical, allocatable :: LdefRST(:)       ! Restart file
#ifdef SP4DVAR
        logical, allocatable :: LdefSCA(:)       ! AD scratch state file
        logical, allocatable :: LdefSCT(:)       ! TL scratch state file
        logical, allocatable :: LdefSPA(:)       ! AD Arnoldi state file
        logical, allocatable :: LdefSPT(:)       ! TL Arnoldi state file
#endif
        logical, allocatable :: LdefSTA(:)       ! Stations file
        logical, allocatable :: LdefTIDE(:)      ! tide forcing file
        logical, allocatable :: LdefTLM(:)       ! Tangent linear file
        logical, allocatable :: LdefTLF(:)       ! TLM/RPM impulse file

        logical, allocatable :: LreadADM(:)      ! Read ADM multi-files
        logical, allocatable :: LreadBLK(:)      ! Read NLM bulk fluxes
        logical, allocatable :: LreadFRC(:)      ! Read FRC files
        logical, allocatable :: LreadFWD(:)      ! Read FWD trajectory
        logical, allocatable :: LreadQCK(:)      ! Read QCK trajectory
        logical, allocatable :: LreadTLM(:)      ! Read TLM multi-files

        logical, allocatable :: LwrtADJ(:)       ! Write adjoint file
        logical, allocatable :: LwrtAVG(:)       ! Write average file
        logical, allocatable :: LwrtDIA(:)       ! Write diagnostic file
        logical, allocatable :: LwrtHIS(:)       ! Write history file
        logical, allocatable :: LwrtPER(:)       ! Write during ensemble
        logical, allocatable :: LwrtQCK(:)       ! write quicksave file
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
!  Switch to write 4D-Var cost functions to DAV NetCDF.
!
        logical, allocatable :: LwrtCost(:)
!
!  Switch to write out adjoint surface forcing fields adjusted by the
!  4DVAR algorithms.
!
        logical, allocatable :: Ladjusted(:)
!
!  Switch to append information to an existing ROMS standard output
!  log file.
!
        logical :: Lappend = .FALSE.
!
!  Switch to read input open boundary conditions data.
!
        logical, allocatable :: LprocessOBC(:)
!
!  Switch to read input tidal forcing data.
!
        logical, allocatable :: LprocessTides(:)
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
        integer, allocatable :: ndefDIA(:)       ! Diagnostics file
        integer, allocatable :: ndefHIS(:)       ! History file
        integer, allocatable :: ndefQCK(:)       ! Quicksave file
        integer, allocatable :: ndefTLM(:)       ! Tangent linear file
!
!  Starting timestep for accumulation of output.
!
        integer, allocatable :: ntsAVG(:)        ! Average file
        integer, allocatable :: ntsDIA(:)        ! Diagnostics file
!
!  Number of timesteps between writing of output data.
!
        integer, allocatable :: nADJ(:)          ! Adjoint file
        integer, allocatable :: nAVG(:)          ! Average file
        integer, allocatable :: nDIA(:)          ! Diagnostics file
        integer, allocatable :: nFLT(:)          ! Floats file
        integer, allocatable :: nHIS(:)          ! History file
        integer, allocatable :: nQCK(:)          ! Quicksave file
        integer, allocatable :: nRST(:)          ! Restart file
        integer, allocatable :: nSTA(:)          ! Stations file
        integer, allocatable :: nTLM(:)          ! Tangent linear file
#ifdef SP4DVAR
        integer, allocatable :: nSCA(:)          ! AD scratch state file
        integer, allocatable :: nSCT(:)          ! TL scratch state file
        integer, allocatable :: nSPA(:)          ! AD Arnoldi state file
        integer, allocatable :: nSPT(:)          ! TL Arnoldi state file
        integer, allocatable :: Crec(:)          ! Convolution records
#endif
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
        real(dp), allocatable :: OBC_time(:,:)
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
        real(dp), allocatable :: SF_time(:,:)
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

#if defined AVERAGES  && defined AVERAGES_DETIDE && \
   (defined SSH_TIDES || defined UV_TIDES)
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

        real(dp), allocatable :: obcfac(:)
        real(dp), allocatable :: FSobc_in(:,:)
        real(dp), allocatable :: FSobc_out(:,:)
        real(dp), allocatable :: M2obc_in(:,:)
        real(dp), allocatable :: M2obc_out(:,:)
#ifdef SOLVE3D
        real(dp), allocatable :: M3obc_in(:,:)
        real(dp), allocatable :: M3obc_out(:,:)
        real(dp), allocatable :: Tobc_in(:,:,:)
        real(dp), allocatable :: Tobc_out(:,:,:)
#endif
!
!  Inverse time-scales (1/s) for nudging at open boundaries and sponge
!  areas.
!
        real(dp), allocatable :: Znudg(:)          ! Free-surface
        real(dp), allocatable :: M2nudg(:)         ! 2D momentum
        real(dp), allocatable :: M3nudg(:)         ! 3D momentum
        real(dp), allocatable :: Tnudg(:,:)        ! Tracers
        real(dp), allocatable :: Tnudg_SSS(:)      ! Sea surface salinity
        real(dp), allocatable :: sss_mismatch_threshold(:)! Sea surface salinity
!
!  Variables used to impose mass flux conservation in open boundary
!  configurations.
!
        real(dp) :: bc_area = 0.0_dp
        real(dp) :: bc_flux = 0.0_dp
        real(dp) :: ubar_xs = 0.0_dp
#if defined TANGENT || defined TL_IOMS
        real(dp) :: tl_bc_area = 0.0_dp
        real(dp) :: tl_bc_flux = 0.0_dp
        real(dp) :: tl_ubar_xs = 0.0_dp
#endif
#ifdef ADJOINT
        real(dp) :: ad_bc_area = 0.0_dp
        real(dp) :: ad_bc_flux = 0.0_dp
        real(dp) :: ad_ubar_xs = 0.0_dp
#endif
#ifdef SP4DVAR
        real(r8) :: sp_sumdot = 0.0_r8
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
        real(dp) :: blk_Cpa = 1004.67_dp      ! (J/kg/K), Businger 1982
        real(dp) :: blk_Cpw = 4000.0_dp       ! (J/kg/K)
        real(dp) :: blk_Rgas = 287.1_dp       ! (J/kg/K)
        real(dp) :: blk_Zabl = 600.0_dp       ! (m)
        real(dp) :: blk_beta = 1.2_dp         ! non-dimensional
        real(dp) :: blk_dter = 0.3_dp         ! (K)
        real(dp) :: blk_tcw = 0.6_dp          ! (W/m/K)
        real(dp) :: blk_visw = 0.000001_dp    ! (m2/s)

        real(r8), allocatable :: blk_ZQ(:)        ! (m)
        real(r8), allocatable :: blk_ZT(:)        ! (m)
        real(r8), allocatable :: blk_ZW(:)        ! (m)
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

        real(dp), parameter :: sg_pi = pi

        real(dp) :: sg_Cdmax = 0.01_dp         ! non-dimensional
        real(dp) :: sg_alpha = 1.0_dp          ! non-dimensional
        real(dp) :: sg_g = 9.81_dp             ! (m/s2)
        real(dp) :: sg_kappa = 0.41_dp         ! non-dimensional
        real(dp) :: sg_nu = 0.00000119_dp      ! (m2/s)
        real(dp) :: sg_tol = 0.0001_dp         ! non-dimensional
        real(dp) :: sg_ustarcdef = 0.01_dp     ! (m/s)
        real(dp) :: sg_z100 = 1.0_dp           ! (m)
        real(dp) :: sg_z1min = 0.20_dp         ! (m)
        real(dp) :: sg_z1p                     ! non-dimensional
        real(dp) :: sg_znotcdef = 0.01_dp      ! (m)
        real(dp) :: sg_znotdef                 ! (m)

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
     &               0.42_r8, 0.37_r8, 0.33_r8, 0.00468592_r8 /)

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
      IF (.not.allocated(SCALARS)) THEN
        allocate ( SCALARS(Ngrids) )

        DO ng=1,Ngrids

#ifdef SOLVE3D
          allocate ( SCALARS(ng) % Fstate(8+2*MT) )
          Dmem(ng)=Dmem(ng)+REAL(8+2*MT,r8)
          SCALARS(ng) % Fstate = .FALSE.

          allocate ( SCALARS(ng) % Lstate(8+2*MT) )
          Dmem(ng)=Dmem(ng)+REAL(8+2*MT,r8)
          SCALARS(ng) % Lstate = .FALSE.
#else
          allocate ( SCALARS(ng) % Fstate(5) )
          Dmem(ng)=Dmem(ng)+5.0_r8
          SCALARS(ng) % Fstate = .FALSE.

          allocate ( SCALARS(ng) % Lstate(3) )
          Dmem(ng)=Dmem(ng)+3.0_r8
          SCALARS(ng) % Lstate = .FALSE.
#endif

#if defined READ_WATER && defined MASKING && defined DISTRIBUTE
          allocate ( SCALARS(ng) % IJwater((Lm(ng)+2)*(Mm(ng)+2),4) )
          Dmem(ng)=Dmem(ng)+4.0_r8*REAL((Lm(ng)+2)*(Mm(ng)+2),r8)
          SCALARS(ng) % IJwater = 0
#endif

          allocate ( SCALARS(ng) % Cs_r(N(ng)) )
          Dmem(ng)=Dmem(ng)+REAL(N(ng),r8)
          SCALARS(ng) % Cs_r = IniVal

          allocate ( SCALARS(ng) % Cs_w(0:N(ng)) )
          Dmem(ng)=Dmem(ng)+REAL(N(ng)+1,r8)
          SCALARS(ng) % Cs_w = IniVal

          allocate ( SCALARS(ng) % sc_r(N(ng)) )
          Dmem(ng)=Dmem(ng)+REAL(N(ng),r8)
          SCALARS(ng) % sc_r = IniVal

          allocate ( SCALARS(ng) % sc_w(0:N(ng)) )
          Dmem(ng)=Dmem(ng)+REAL(N(ng)+1,r8)
          SCALARS(ng) % sc_w = IniVal

        END DO
      END IF
!
!  Allocate variables that require special treatment in shared-memory.
!  These variables are private for each thread to avoid collisions.
!
!$OMP PARALLEL
      IF (.not.allocated(PREDICTOR_2D_STEP)) THEN
        allocate ( PREDICTOR_2D_STEP(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

      IF (.not.allocated(indx1)) THEN
        allocate ( indx1(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(iic)) THEN
        allocate ( iic(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(iif)) THEN
        allocate ( iif(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

      IF (.not.allocated(FrcRec)) THEN
        allocate ( FrcRec(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(SOrec)) THEN
        allocate ( SOrec(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
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
      IF (.not.associated(inert)) THEN
        allocate ( inert(NPT) )
        Dmem(1)=Dmem(1)+REAL(NPT,r8)
      END IF
#endif
      IF (.not.allocated(PerfectRST)) THEN
        allocate ( PerfectRST(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

      IF (.not.allocated(ndtfast)) THEN
        allocate ( ndtfast(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(nfast)) THEN
        allocate ( nfast(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

      IF (.not.allocated(dt)) THEN
        allocate ( dt(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(dtfast)) THEN
        allocate ( dtfast(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(TimeEnd)) THEN
        allocate ( TimeEnd(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(AVGtime)) THEN
        allocate ( AVGtime(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(DIAtime)) THEN
        allocate ( DIAtime(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(IMPtime)) THEN
        allocate ( IMPtime(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(INItime)) THEN
        allocate ( INItime(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(INItimeS)) THEN
        allocate ( INItimeS(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(ObsTime)) THEN
        allocate ( ObsTime(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(FrcTime)) THEN
        allocate ( FrcTime(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

      IF (.not.allocated(ntimes)) THEN
        allocate ( ntimes(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(first_time)) THEN
        allocate ( first_time(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

      IF (.not.allocated(ntfirst)) THEN
        allocate ( ntfirst(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(ntstart)) THEN
        allocate ( ntstart(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(ntend)) THEN
        allocate ( ntend(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
#ifdef TIDES_ASTRO
      IF (.not.allocated(FIRST_TIDES_ASTRO)) THEN
        allocate ( FIRST_TIDES_ASTRO(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
#endif
#ifdef AGE_DISTRIBUTION
      IF (.not.allocated(dtau)) THEN
        allocate ( dtau(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
#endif

!$OMP PARALLEL
      IF (.not.allocated(synchro_flag)) THEN
        allocate ( synchro_flag(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(step_counter)) THEN
        allocate ( step_counter(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(tdays)) THEN
        allocate ( tdays(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(time)) THEN
        allocate ( time(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(time_code)) THEN
        allocate ( time_code(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
!$OMP END PARALLEL

      IF (.not.allocated(NrecFrc)) THEN
        allocate ( NrecFrc(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

      IF (.not.allocated(SetGridConfig)) THEN
        allocate ( SetGridConfig(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

      IF (.not.allocated(NudgingCoeff)) THEN
        allocate ( NudgingCoeff(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

      IF (.not.allocated(ObcData)) THEN
        allocate ( ObcData(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

      IF (.not.allocated(Lbiology)) THEN
        allocate ( Lbiology(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(Lfloats)) THEN
        allocate ( Lfloats(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(Lsediment)) THEN
        allocate ( Lsediment(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(Lstations)) THEN
        allocate ( Lstations(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

      IF (.not.allocated(CompositeGrid)) THEN
        allocate ( CompositeGrid(4,Ngrids) )
        Dmem(1)=Dmem(1)+4.0_r8*REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(RefinedGrid)) THEN
        allocate ( RefinedGrid(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(RefineScale)) THEN
        allocate ( RefineScale(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(GetDonorData)) THEN
        allocate ( GetDonorData(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

      IF (.not.allocated(EWperiodic)) THEN
        allocate ( EWperiodic(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(NSperiodic)) THEN
        allocate ( NSperiodic(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

      IF (.not.allocated(VolCons)) THEN
        allocate ( VolCons(4,Ngrids) )
        Dmem(1)=Dmem(1)+4.0_r8*REAL(Ngrids,r8)
      END IF
#if defined ADJOINT || defined TANGENT || defined TL_IOMS
      IF (.not.allocated(ad_VolCons)) THEN
        allocate ( ad_VolCons(4,Ngrids) )
        Dmem(1)=Dmem(1)+4.0_r8*REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(tl_VolCons)) THEN
        allocate ( tl_VolCons(4,Ngrids) )
        Dmem(1)=Dmem(1)+4.0_r8*REAL(Ngrids,r8)
      END IF
#endif

      IF (.not.allocated(Lsponge)) THEN
        allocate ( Lsponge(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(LuvSponge)) THEN
        allocate ( LuvSponge(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(LtracerSponge)) THEN
        allocate ( LtracerSponge(MT,Ngrids) )
        Dmem(1)=Dmem(1)+REAL(MT*Ngrids,r8)
      END IF

      IF (.not.allocated(CLM_FILE)) THEN
        allocate ( CLM_FILE(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(Lclimatology)) THEN
        allocate ( Lclimatology(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(LsshCLM)) THEN
        allocate ( LsshCLM(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(Lm2CLM)) THEN
        allocate ( Lm2CLM(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(Lm3CLM)) THEN
        allocate ( Lm3CLM(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(LtracerCLM)) THEN
        allocate ( LtracerCLM(MT,Ngrids) )
        Dmem(1)=Dmem(1)+REAL(MT*Ngrids,r8)
      END IF

      IF (.not.allocated(Lnudging)) THEN
        allocate ( Lnudging(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(LnudgeM2CLM)) THEN
        allocate ( LnudgeM2CLM(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(LnudgeM3CLM)) THEN
        allocate ( LnudgeM3CLM(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(LnudgeTCLM)) THEN
        allocate ( LnudgeTCLM(MT,Ngrids) )
        Dmem(1)=Dmem(1)+REAL(MT*Ngrids,r8)
      END IF

      IF (.not.allocated(LuvSrc)) THEN
        allocate ( LuvSrc(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(LwSrc)) THEN
        allocate ( LwSrc(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(LtracerSrc)) THEN
        allocate ( LtracerSrc(MT,Ngrids) )
        Dmem(1)=Dmem(1)+REAL(MT*Ngrids,r8)
      END IF

      IF (.not.allocated(maxspeed)) THEN
        allocate ( maxspeed(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(maxrho)) THEN
        allocate ( maxrho(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

      IF (.not.allocated(levsfrc)) THEN
        allocate ( levsfrc(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(levbfrc)) THEN
        allocate ( levbfrc(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

      IF (.not.allocated(Vtransform)) THEN
        allocate ( Vtransform(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(Vstretching)) THEN
        allocate ( Vstretching(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

      IF (.not.allocated(Tcline)) THEN
        allocate ( Tcline(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(hc)) THEN
        allocate ( hc(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(theta_s)) THEN
        allocate ( theta_s(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(theta_b)) THEN
        allocate ( theta_b(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

      IF (.not.allocated(hmin)) THEN
        allocate ( hmin(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(hmax)) THEN
        allocate ( hmax(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

      IF (.not.allocated(xl)) THEN
        allocate ( xl(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(el)) THEN
        allocate ( el(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

      IF (.not.allocated(LonMin)) THEN
        allocate ( LonMin(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(LonMax)) THEN
        allocate ( LonMax(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(LatMin)) THEN
        allocate ( LatMin(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(LatMax)) THEN
        allocate ( LatMax(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

      IF (.not.allocated(Idigits)) THEN
        allocate ( Idigits(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(Jdigits)) THEN
        allocate ( Jdigits(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
#ifdef SOLVE3D
      IF (.not.allocated(Kdigits)) THEN
        allocate ( Kdigits(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
#endif
      IF (.not.allocated(TotVolume)) THEN
        allocate ( TotVolume(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(MinVolume)) THEN
        allocate ( MinVolume(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(MaxVolume)) THEN
        allocate ( MaxVolume(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

      IF (.not.allocated(DXmin)) THEN
        allocate ( DXmin(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(DXmax)) THEN
        allocate ( DXmax(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(DYmin)) THEN
        allocate ( DYmin(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(DYmax)) THEN
        allocate ( DYmax(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

#ifdef MASKING
      IF (.not.allocated(DXminW)) THEN
        allocate ( DXminW(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(DXmaxW)) THEN
        allocate ( DXmaxW(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(DYminW)) THEN
        allocate ( DYminW(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(DYmaxW)) THEN
        allocate ( DYmaxW(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
#endif


#ifdef SOLVE3D
      IF (.not.allocated(DZmin)) THEN
        allocate ( DZmin(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(DZmax)) THEN
        allocate ( DZmax(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
# ifdef MASKING
      IF (.not.allocated(DZminW)) THEN
        allocate ( DZminW(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(DZmaxW)) THEN
        allocate ( DZmaxW(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
# endif
#endif

      IF (.not.allocated(grdmax)) THEN
        allocate ( grdmax(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

#ifdef DIFF_3DCOEF
      IF (.not.allocated(DiffMin)) THEN
        allocate ( DiffMin(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(DiffMax)) THEN
        allocate ( DiffMax(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
#endif

      IF (.not.allocated(Cg_min)) THEN
        allocate ( Cg_min(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(Cg_max)) THEN
        allocate ( Cg_max(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(Cg_Cor)) THEN
        allocate ( Cg_Cor(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

#ifdef VISC_3DCOEF
      IF (.not.allocated(ViscMin)) THEN
        allocate ( ViscMin(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(ViscMax)) THEN
        allocate ( ViscMax(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
#endif

      IF (.not.allocated(R0)) THEN
        allocate ( R0(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(Tcoef)) THEN
        allocate ( Tcoef(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(Scoef)) THEN
        allocate ( Scoef(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

      IF (.not.allocated(T0)) THEN
        allocate ( T0(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(S0)) THEN
        allocate ( S0(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

      IF (.not.allocated(gamma2)) THEN
        allocate ( gamma2(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

      IF (.not.allocated(lmd_Jwt)) THEN
        allocate ( lmd_Jwt(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

      IF (.not.allocated(rx0)) THEN
        allocate ( rx0(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(rx1)) THEN
        allocate ( rx1(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

      IF (.not.allocated(rdrg)) THEN
        allocate ( rdrg(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(rdrg2)) THEN
        allocate ( rdrg2(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

      IF (.not.allocated(Zos)) THEN
        allocate ( Zos(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(Zob)) THEN
        allocate ( Zob(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

      IF (.not.allocated(Dcrit)) THEN
        allocate ( Dcrit(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

#ifdef PROPAGATOR
      IF (.not.allocated(Nconv)) THEN
        allocate ( Nconv(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
#endif

      IF (.not.allocated(weight)) THEN
        allocate ( weight(2,0:256,Ngrids) )
        Dmem(1)=Dmem(1)+514.0_r8*REAL(Ngrids,r8)
      END IF

      IF (.not.allocated(Akk_bak)) THEN
        allocate ( Akk_bak(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(Akp_bak)) THEN
        allocate ( Akp_bak(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(Akv_bak)) THEN
        allocate ( Akv_bak(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(Akv_limit)) THEN
        allocate ( Akv_limit(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

      IF (.not.allocated(ad_visc2)) THEN
        allocate ( ad_visc2(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(nl_visc2)) THEN
        allocate ( nl_visc2(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(tl_visc2)) THEN
        allocate ( tl_visc2(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(visc2)) THEN
        allocate ( visc2(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

      IF (.not.allocated(ad_visc4)) THEN
        allocate ( ad_visc4(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(nl_visc4)) THEN
        allocate ( nl_visc4(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(tl_visc4)) THEN
        allocate ( tl_visc4(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(visc4)) THEN
        allocate ( visc4(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

      IF (.not.allocated(tkenu2)) THEN
        allocate ( tkenu2(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(tkenu4)) THEN
        allocate ( tkenu4(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

      IF (.not.allocated(Akt_bak)) THEN
        allocate ( Akt_bak(MT,Ngrids) )
        Dmem(1)=Dmem(1)+REAL(MT*Ngrids,r8)
      END IF
      IF (.not.allocated(Akt_limit)) THEN
        allocate ( Akt_limit(NAT,Ngrids) )
        Dmem(1)=Dmem(1)+REAL(NAT*Ngrids,r8)
      END IF
      IF (.not.allocated(Kdiff)) THEN
        allocate ( Kdiff(MT,Ngrids) )
        Dmem(1)=Dmem(1)+REAL(MT*Ngrids,r8)
      END IF

      IF (.not.allocated(ad_tnu2)) THEN
        allocate ( ad_tnu2(MT,Ngrids) )
        Dmem(1)=Dmem(1)+REAL(MT*Ngrids,r8)
      END IF
      IF (.not.allocated(nl_tnu2)) THEN
        allocate ( nl_tnu2(MT,Ngrids) )
        Dmem(1)=Dmem(1)+REAL(MT*Ngrids,r8)
      END IF
      IF (.not.allocated(tl_tnu2)) THEN
        allocate ( tl_tnu2(MT,Ngrids) )
        Dmem(1)=Dmem(1)+REAL(MT*Ngrids,r8)
      END IF
      IF (.not.allocated(tnu2)) THEN
        allocate ( tnu2(MT,Ngrids) )
        Dmem(1)=Dmem(1)+REAL(MT*Ngrids,r8)
      END IF

      IF (.not.allocated(ad_tnu4)) THEN
        allocate ( ad_tnu4(MT,Ngrids) )
        Dmem(1)=Dmem(1)+REAL(MT*Ngrids,r8)
      END IF
      IF (.not.allocated(nl_tnu4)) THEN
        allocate ( nl_tnu4(MT,Ngrids) )
        Dmem(1)=Dmem(1)+REAL(MT*Ngrids,r8)
      END IF
      IF (.not.allocated(tl_tnu4)) THEN
        allocate ( tl_tnu4(MT,Ngrids) )
        Dmem(1)=Dmem(1)+REAL(MT*Ngrids,r8)
      END IF
      IF (.not.allocated(tnu4)) THEN
        allocate ( tnu4(MT,Ngrids) )
        Dmem(1)=Dmem(1)+REAL(MT*Ngrids,r8)
      END IF

      IF (.not.allocated(tl_M2diff)) THEN
        allocate ( tl_M2diff(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(tl_M3diff)) THEN
        allocate ( tl_M3diff(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(tl_Tdiff)) THEN
        allocate ( tl_Tdiff(MT,Ngrids) )
        Dmem(1)=Dmem(1)+REAL(MT*Ngrids,r8)
      END IF

      IF (.not.allocated(ad_Akv_fac)) THEN
        allocate ( ad_Akv_fac(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(tl_Akv_fac)) THEN
        allocate ( tl_Akv_fac(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

      IF (.not.allocated(ad_Akt_fac)) THEN
        allocate ( ad_Akt_fac(MT,Ngrids) )
        Dmem(1)=Dmem(1)+REAL(MT*Ngrids,r8)
      END IF
      IF (.not.allocated(tl_Akt_fac)) THEN
        allocate ( tl_Akt_fac(MT,Ngrids) )
        Dmem(1)=Dmem(1)+REAL(MT*Ngrids,r8)
      END IF

      IF (.not.allocated(LdefADJ)) THEN
        allocate ( LdefADJ(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(LdefAVG)) THEN
        allocate ( LdefAVG(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(LdefDAI)) THEN
        allocate ( LdefDAI(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(LdefDIA)) THEN
        allocate ( LdefDIA(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(LdefERR)) THEN
        allocate ( LdefERR(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(LdefFLT)) THEN
        allocate ( LdefFLT(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(LdefHIS)) THEN
        allocate ( LdefHIS(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(LdefHSS)) THEN
        allocate ( LdefHSS(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(LdefINI)) THEN
        allocate ( LdefINI(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(LdefIRP)) THEN
        allocate ( LdefIRP(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(LdefITL)) THEN
        allocate ( LdefITL(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(LdefLCZ)) THEN
        allocate ( LdefLCZ(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(LdefLZE)) THEN
        allocate ( LdefLZE(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(LdefMOD)) THEN
        allocate ( LdefMOD(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(LdefQCK)) THEN
        allocate ( LdefQCK(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(LdefRST)) THEN
        allocate ( LdefRST(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
#ifdef SP4DVAR
      IF (.not.allocated(LdefSCA)) THEN
        allocate ( LdefSCA(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(LdefSCT)) THEN
        allocate ( LdefSCT(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(LdefSPA)) THEN
        allocate ( LdefSPA(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(LdefSPT)) THEN
        allocate ( LdefSPT(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
#endif
      IF (.not.allocated(LdefSTA)) THEN
        allocate ( LdefSTA(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(LdefTIDE)) THEN
        allocate ( LdefTIDE(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(LdefTLM)) THEN
        allocate ( LdefTLM(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(LdefTLF)) THEN
        allocate ( LdefTLF(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

      IF (.not.allocated(LreadADM)) THEN
        allocate ( LreadADM(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(LreadBLK)) THEN
        allocate ( LreadBLK(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(LreadFRC)) THEN
        allocate ( LreadFRC(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(LreadFWD)) THEN
        allocate ( LreadFWD(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(LreadQCK)) THEN
        allocate ( LreadQCK(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(LreadTLM)) THEN
        allocate ( LreadTLM(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

      IF (.not.allocated(LwrtADJ)) THEN
        allocate ( LwrtADJ(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(LwrtAVG)) THEN
        allocate ( LwrtAVG(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(LwrtDIA)) THEN
        allocate ( LwrtDIA(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(LwrtHIS)) THEN
        allocate ( LwrtHIS(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(LwrtPER)) THEN
        allocate ( LwrtPER(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(LwrtQCK)) THEN
        allocate ( LwrtQCK(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(LwrtRST)) THEN
        allocate ( LwrtRST(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(LwrtTLM)) THEN
        allocate ( LwrtTLM(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(LwrtTLF)) THEN
        allocate ( LwrtTLF(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

      IF (.not.allocated(LdefNRM)) THEN
        allocate ( LdefNRM(4,Ngrids) )
        Dmem(1)=Dmem(1)+4.0_r8*REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(LwrtNRM)) THEN
        allocate ( LwrtNRM(4,Ngrids) )
        Dmem(1)=Dmem(1)+4.0_r8*REAL(Ngrids,r8)
      END IF

#if defined STOCHASTIC_OPT && !defined STOCH_OPT_WHITE
      IF (.not.allocated(LwrtState3d)) THEN
        allocate ( LwrtState3d(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(SOinitial)) THEN
        allocate ( SOinitial(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
#endif

      IF (.not.allocated(LwrtState2d)) THEN
        allocate ( LwrtState2d(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(LwrtTime)) THEN
        allocate ( LwrtTime(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(LwrtCost)) THEN
        allocate ( LwrtCost(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

      IF (.not.allocated(Ladjusted)) THEN
        allocate ( Ladjusted(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

      IF (.not.allocated(LprocessOBC)) THEN
        allocate ( LprocessOBC(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(LprocessTides)) THEN
        allocate ( LprocessTides(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

      IF (.not.allocated(LwrtInfo)) THEN
        allocate ( LwrtInfo(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

      IF (.not.allocated(ldefout)) THEN
        allocate ( ldefout(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

      IF (.not.allocated(ndefADJ)) THEN
        allocate ( ndefADJ(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(ndefAVG)) THEN
        allocate ( ndefAVG(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(ndefDIA)) THEN
        allocate ( ndefDIA(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(ndefHIS)) THEN
        allocate ( ndefHIS(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(ndefQCK)) THEN
        allocate ( ndefQCK(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(ndefTLM)) THEN
        allocate ( ndefTLM(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

      IF (.not.allocated(ntsAVG)) THEN
        allocate ( ntsAVG(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(ntsDIA)) THEN
        allocate ( ntsDIA(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

      IF (.not.allocated(nADJ)) THEN
        allocate ( nADJ(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(nAVG)) THEN
        allocate ( nAVG(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(nDIA)) THEN
        allocate ( nDIA(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(nFLT)) THEN
        allocate ( nFLT(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(nHIS)) THEN
        allocate ( nHIS(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(nQCK)) THEN
        allocate ( nQCK(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(nRST)) THEN
        allocate ( nRST(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(nSTA)) THEN
        allocate ( nSTA(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(nTLM)) THEN
        allocate ( nTLM(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
#ifdef SP4DVAR
      IF (.not.allocated(nSCA)) THEN
        allocate ( nSCA(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(nSCT)) THEN
        allocate ( nSCT(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(nSPA)) THEN
        allocate ( nSPA(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(nSPT)) THEN
        allocate ( nSPT(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(Crec)) THEN
        allocate ( Crec(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
#endif

      IF (.not.allocated(ninfo)) THEN
        allocate ( ninfo(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

      IF (.not.allocated(nOBC)) THEN
        allocate ( nOBC(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(Nbrec)) THEN
        allocate ( Nbrec(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(OBCcount)) THEN
        allocate ( OBCcount(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

#ifdef ADJUST_BOUNDARY
      IF (.not.allocated(Lobc)) THEN
        allocate ( Lobc(4,MstateVar,Ngrids) )
        Dmem(1)=Dmem(1)+4.0_r8*REAL(MstateVar*Ngrids,r8)
      END IF
#endif

      IF (.not.allocated(nSFF)) THEN
        allocate ( nSFF(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(Nfrec)) THEN
        allocate ( Nfrec(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(SFcount)) THEN
        allocate ( SFcount(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

#ifdef ADJUST_STFLUX
      IF (.not.allocated(Lstflux)) THEN
        allocate ( Lstflux(MT,Ngrids) )
        Dmem(1)=Dmem(1)+REAL(MT*Ngrids,r8)
      END IF
#endif

      IF (.not.allocated(nrrec)) THEN
        allocate ( nrrec(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

      IF (.not.allocated(LastRec)) THEN
        allocate ( LastRec(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

      IF (.not.allocated(LcycleADJ)) THEN
        allocate ( LcycleADJ(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(LcycleRST)) THEN
        allocate ( LcycleRST(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(LcycleTLM)) THEN
        allocate ( LcycleTLM(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

#if defined AVERAGES  && defined AVERAGES_DETIDE && \
   (defined SSH_TIDES || defined UV_TIDES)
      IF (.not.allocated(Hcount)) THEN
        allocate ( Hcount(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
#endif

      IF (.not.allocated(KstrS)) THEN
        allocate ( KstrS(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(KendS)) THEN
        allocate ( KendS(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

      IF (.not.allocated(DstrS)) THEN
        allocate ( DstrS(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(DendS)) THEN
        allocate ( DendS(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

      IF (.not.allocated(TRnorm)) THEN
        allocate ( TRnorm(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

      IF (.not.allocated(SO_decay)) THEN
        allocate ( SO_decay(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
#ifdef SOLVE3D
      IF (.not.allocated(SO_sdev)) THEN
        allocate ( SO_sdev(7+2*MT,Ngrids) )
        Dmem(1)=Dmem(1)+REAL((7+2*MT)*Ngrids,r8)
      END IF
#else
      IF (.not.allocated(SO_sdev)) THEN
        allocate ( SO_sdev(5,Ngrids) )
        Dmem(1)=Dmem(1)+5.0_r8*REAL(Ngrids,r8)
      END IF
#endif

#if defined FOUR_DVAR || defined VERIFICATION
      IF (.not.allocated(Cnorm)) THEN
        allocate ( Cnorm(2,MstateVar) )
        Dmem(1)=Dmem(1)+2.0_r8*REAL(MstateVar,r8)
      END IF
      IF (.not.allocated(CnormB)) THEN
        allocate ( CnormB(MstateVar,4) )
        Dmem(1)=Dmem(1)+4.0_r8*REAL(MstateVar,r8)
      END IF

      IF (.not.allocated(SporadicImpulse)) THEN
        allocate ( SporadicImpulse(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(FrequentImpulse)) THEN
        allocate ( FrequentImpulse(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

      IF (.not.allocated(dTdz_min)) THEN
        allocate ( dTdz_min(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(ml_depth)) THEN
        allocate ( ml_depth(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

      IF (.not.allocated(LNM_depth)) THEN
        allocate ( LNM_depth(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

      IF (.not.allocated(balance)) THEN
        allocate ( balance(MstateVar) )
        Dmem(1)=Dmem(1)+REAL(MstateVar,r8)
      END IF

      IF (.not.allocated(Hdecay)) THEN
        allocate ( Hdecay(2,MstateVar,Ngrids) )
        Dmem(1)=Dmem(1)+2.0_r8*REAL(MstateVar*Ngrids,r8)
      END IF
      IF (.not.allocated(Vdecay)) THEN
        allocate ( Vdecay(2,MstateVar,Ngrids) )
        Dmem(1)=Dmem(1)+2.0_r8*REAL(MstateVar*Ngrids,r8)
      END IF

      IF (.not.allocated(Tdecay)) THEN
        allocate ( Tdecay(MstateVar,Ngrids) )
        Dmem(1)=Dmem(1)+REAL(MstateVar*Ngrids,r8)
      END IF

      IF (.not.allocated(HdecayB)) THEN
        allocate ( HdecayB(MstateVar,4,Ngrids) )
        Dmem(1)=Dmem(1)+4.0_r8*REAL(MstateVar*Ngrids,r8)
      END IF
      IF (.not.allocated(VdecayB)) THEN
        allocate ( VdecayB(MstateVar,4,Ngrids) )
        Dmem(1)=Dmem(1)+4.0_r8*REAL(MstateVar*Ngrids,r8)
      END IF

      IF (.not.allocated(bgqc_type)) THEN
        allocate ( bgqc_type(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(Nprovenance)) THEN
        allocate ( Nprovenance(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(S_bgqc)) THEN
        allocate ( S_bgqc(MstateVar,Ngrids) )
        Dmem(1)=Dmem(1)+REAL(MstateVar*Ngrids,r8)
      END IF
#endif

      IF (.not.allocated(obcfac)) THEN
        allocate ( obcfac(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(FSobc_in)) THEN
        allocate ( FSobc_in(Ngrids,4) )
        Dmem(1)=Dmem(1)+4.0_r8*REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(FSobc_out)) THEN
        allocate ( FSobc_out(Ngrids,4) )
        Dmem(1)=Dmem(1)+4.0_r8*REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(M2obc_in)) THEN
        allocate ( M2obc_in(Ngrids,4) )
        Dmem(1)=Dmem(1)+4.0_r8*REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(M2obc_out)) THEN
        allocate ( M2obc_out(Ngrids,4) )
        Dmem(1)=Dmem(1)+4.0_r8*REAL(Ngrids,r8)
      END IF

#ifdef SOLVE3D
      IF (.not.allocated(M3obc_in)) THEN
        allocate ( M3obc_in(Ngrids,4) )
        Dmem(1)=Dmem(1)+4.0_r8*REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(M3obc_out)) THEN
        allocate ( M3obc_out(Ngrids,4) )
        Dmem(1)=Dmem(1)+4.0_r8*REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(Tobc_in)) THEN
        allocate ( Tobc_in(MT,Ngrids,4) )
        Dmem(1)=Dmem(1)+4.0_r8*REAL(MT*Ngrids,r8)
      END IF
      IF (.not.allocated(Tobc_out)) THEN
        allocate ( Tobc_out(MT,Ngrids,4) )
        Dmem(1)=Dmem(1)+4.0_r8*REAL(MT*Ngrids,r8)
      END IF
#endif

      IF (.not.allocated(Znudg)) THEN
        allocate ( Znudg(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(M2nudg)) THEN
        allocate ( M2nudg(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(M3nudg)) THEN
        allocate ( M3nudg(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(Tnudg)) THEN
        allocate ( Tnudg(MT,Ngrids) )
        Dmem(1)=Dmem(1)+REAL(MT*Ngrids,r8)
      END IF
      IF (.not.allocated(Tnudg_SSS)) THEN
        allocate ( Tnudg_SSS(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not. allocated(sss_mismatch_threshold)) THEN
        allocate ( sss_mismatch_threshold(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

#if defined BULK_FLUXES || defined BULK_FLUXES2D
      IF (.not.allocated(blk_ZQ)) THEN
        allocate ( blk_ZQ(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(blk_ZT)) THEN
        allocate ( blk_ZT(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(blk_ZW)) THEN
        allocate ( blk_ZW(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
#endif
#ifdef BIOLOGY
      IF (.not.allocated(maxbio)) THEN
        allocate ( maxbio(NBT, Ngrids) )
        Dmem(1)=Dmem(1)+REAL(NBT*Ngrids,r8)
      END IF
      IF (.not.allocated(max_bio)) THEN
        allocate ( max_bio(NBT) )
        Dmem(1)=Dmem(1)+REAL(NBT,r8)
      END IF
#endif

      IF (.not.allocated(gls_m)) THEN
        allocate ( gls_m(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(gls_n)) THEN
        allocate ( gls_n(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(gls_p)) THEN
        allocate ( gls_p(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(gls_sigk)) THEN
        allocate ( gls_sigk(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(gls_sigp)) THEN
        allocate ( gls_sigp(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(gls_cmu0)) THEN
        allocate ( gls_cmu0(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(gls_cmupr)) THEN
        allocate ( gls_cmupr(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(gls_c1)) THEN
        allocate ( gls_c1(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(gls_c2)) THEN
        allocate ( gls_c2(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(gls_c3m)) THEN
        allocate ( gls_c3m(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(gls_c3p)) THEN
        allocate ( gls_c3p(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(gls_Kmin)) THEN
        allocate ( gls_Kmin(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(gls_Pmin)) THEN
        allocate ( gls_Pmin(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

      IF (.not.allocated(charnok_alpha)) THEN
        allocate ( charnok_alpha(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(zos_hsig_alpha)) THEN
        allocate ( zos_hsig_alpha(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(sz_alpha)) THEN
        allocate ( sz_alpha(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(crgban_cw)) THEN
        allocate ( crgban_cw(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(wec_alpha)) THEN
        allocate ( wec_alpha(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF


      RETURN
      END SUBROUTINE allocate_scalars

      SUBROUTINE initialize_scalars
!
!=======================================================================
!                                                                      !
!  This routine initializes several variables in module for all nested !
!  grids.                                                              !
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
!#ifdef SWAN_COUPLING
!      Infinity=one/zero
!#endif
      Co=1.0_r8/(2.0_r8+SQRT(2.0_r8))
      gorho0=g/rho0
      DO ng=1,Ngrids
        EWperiodic(ng)=.FALSE.
        NSperiodic(ng)=.FALSE.
        NudgingCoeff(ng)=.FALSE.
        ObcData(ng)=.FALSE.
        SetGridConfig(ng)=.TRUE.
        RefineScale(ng)=0
        gamma2(ng)=-1.0_r8
        Vtransform(ng)=1
        Vstretching(ng)=1
#if defined AVERAGES  && defined AVERAGES_DETIDE && \
   (defined SSH_TIDES || defined UV_TIDES)
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
        maxspeed(ng)=-Large
        maxrho(ng)=-Large
        TotVolume(ng)=0.0_dp
        MinVolume(ng)= Large
        MaxVolume(ng)=-Large
        DXmin(ng)= Large
        DXmax(ng)=-Large
        DYmin(ng)= Large
        DYmax(ng)=-Large
#ifdef MASKING
        DXminW(ng)= Large
        DXmaxW(ng)=-Large
        DYminW(ng)= Large
        DYmaxW(ng)=-Large
#endif
#ifdef SOLVE3D
        DZmin(ng)= Large
        DZmax(ng)=-Large
# ifdef MASKING
        DZminW(ng)= Large
        DZmaxW(ng)=-Large
# endif
#endif
        grdmax(ng)=-Large
#ifdef DIFF_3DCOEF
        DiffMin(ng)= Large
        DiffMax(ng)=-Large
#endif
        Cg_min(ng)= Large
        Cg_max(ng)=-Large
        Cg_Cor(ng)=-Large
#ifdef VISC_3DCOEF
        ViscMin(ng)= Large
        ViscMax(ng)=-Large
#endif
        rx0(ng)=-Large
        rx1(ng)=-Large
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
        DO itrc=1,NAT
          Akt_limit(itrc,ng)=1.0E-3_r8
        END DO
        Akv_limit(ng)=1.0E-3_r8
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
          FSobc_in (ng,i)=0.0_dp
          FSobc_out(ng,i)=0.0_dp
          M2obc_in (ng,i)=0.0_dp
          M2obc_out(ng,i)=0.0_dp
#ifdef SOLVE3D
          M3obc_in (ng,i)=0.0_dp
          M3obc_out(ng,i)=0.0_dp
#endif
        END DO
      END DO
#ifdef SOLVE3D
      Tobc_in = 0.0_dp
      Tobc_out = 0.0_dp
#endif
!
!  Initialize blowup string.
!
      DO i=1,LEN(blowup_string)
        blowup_string(i:i)=' '
      END DO
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

        LprocessOBC(ng)=.FALSE.
        LprocessTides(ng)=.FALSE.

        LdefADJ(ng)=.FALSE.
        LdefAVG(ng)=.TRUE.
        LdefDAI(ng)=.FALSE.
        LdefDIA(ng)=.TRUE.
        LdefERR(ng)=.FALSE.
#ifdef FLOATS
        LdefFLT(ng)=.TRUE.
#endif
        LdefHIS(ng)=.TRUE.
#ifdef FOUR_DVAR
        LdefINI(ng)=.TRUE.
#else
        LdefINI(ng)=.FALSE.
#endif
        LdefIRP(ng)=.FALSE.
        LdefITL(ng)=.FALSE.
        LdefMOD(ng)=.FALSE.
        LdefQCK(ng)=.FALSE.
        LdefRST(ng)=.TRUE.
        LdefSTA(ng)=.TRUE.
        LdefTLM(ng)=.FALSE.
#if defined AVERAGES  && defined AVERAGES_DETIDE && \
   (defined SSH_TIDES || defined UV_TIDES)
        LdefTIDE(ng)=.TRUE.
#else
        LdefTIDE(ng)=.FALSE.
#endif
        LreadADM(ng)=.FALSE.
        LreadBLK(ng)=.FALSE.
#ifdef FRC_FILE
        LreadFRC(ng)=.TRUE.
#else
        LreadFRC(ng)=.FALSE.
#endif
        LreadFWD(ng)=.FALSE.
        LreadQCK(ng)=.FALSE.
        LreadTLM(ng)=.FALSE.
        LwrtADJ(ng)=.FALSE.
        LwrtAVG(ng)=.FALSE.
        LwrtDIA(ng)=.FALSE.
        LwrtHIS(ng)=.FALSE.
        LwrtPER(ng)=.FALSE.
        LwrtQCK(ng)=.FALSE.
        LwrtRST(ng)=.FALSE.
        LwrtTLM(ng)=.FALSE.
        LwrtInfo(ng)=.TRUE.
        LwrtState2d(ng)=.FALSE.
#if defined STOCHASTIC_OPT && !defined STOCH_OPT_WHITE
        LwrtState3d(ng)=.FALSE.
        SOinitial(ng)=.FALSE.
#endif
        LwrtTime(ng)=.TRUE.
        LwrtCost(ng)=.FALSE.
        ldefout(ng)=.FALSE.

      END DO
!
!  Initialize the NLM initial conditions time to a negative number to
!  check if its value was assigned elsewhere.  It can be used during
!  the initialization of the adjoint model when DSTART is not the
!  same as the start of the simulation.
!
      DO ng=1,Ngrids
        INItime(ng)=-1.0_dp
        INItimeS(ng)=-1.0_dp
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

      RETURN
      END SUBROUTINE initialize_scalars

      END MODULE mod_scalars
