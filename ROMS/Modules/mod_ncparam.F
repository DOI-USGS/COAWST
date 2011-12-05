#include "cppdefs.h"
      MODULE mod_ncparam
!
!svn $Id: mod_ncparam.F 746 2008-09-10 04:09:44Z jcwarner $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2010 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This MODULE contains all the variables associated with input and    !
!  output  NetCDF  files.  The IO model is very generic and easy to    !
!  change or expand.  The NetCDF files can be in any language.  All    !
!  the IO information is managed using the following variables:        !
!                                                                      !
!  Vname      Input/output variables names and attributes:             !
!               Vname(1,*)  => field variable name.                    !
!               Vname(2,*)  => long-name attribute.                    !
!               Vname(3,*)  => units attribute.                        !
!               Vname(4,*)  => field type attribute.                   !
!               Vname(5,*)  => associated time variable name.          !
!  Tname      Input/output associated time variables names.            !
!                                                                      !
!  Cinfo      Input/output file names associated with each field       !
!                                                                      !
!  Linfo      Input/output fields logical information:                 !
!               Linfo(1,*)  => switch indicating grided data.          !
!               Linfo(2,*)  => switch indicating time cycling.         !
!               Linfo(3,*)  => switch indicating only one-time         !
!                              record available.                       !
!               Linfo(4,*)  => switch indication special record        !
!                              processing (like tides)                 !
!                                                                      !
!  Iinfo      Input/output fields integer information:                 !
!               Iinfo(1,*)  => variable grid type.                     !
!               Iinfo(2,*)  => field variable NetCDF ID.               !
!               Iinfo(3,*)  => associated time variable NetCDF ID.     !
!               Iinfo(4,*)  => number of time records.                 !
!               Iinfo(5,*)  => size of first  spatial dimension.       !
!               Iinfo(6,*)  => size of second spatial dimension.       !
!               Iinfo(7,*)  => size of third  spatial dimension.       !
!               Iinfo(8,*)  => rolling two-time levels index.          !
!               Iinfo(9,*)  => latest processed time record .          !
!                                                                      !
!  Finfo      Input/output field floating-point information:           !
!               Finfo( 1,*) => Starting time (days) of data.           !
!               Finfo( 2,*) => Ending time (days) of data.             !
!               Finfo( 3,*) => Data time lower bound (days) enclosing  !
!                                model starting time.                  !
!               Finfo( 4,*) => Data time upper bound (days) enclosing  !
!                                model starting time.                  !
!               Finfo( 5,*) => length (days) of time cycling.          !
!               Finfo( 6,*) => Scale to convert time to day units.     !
!               Finfo( 7,*) => latest monotonic time (sec).            !
!               Finfo( 8,*) => minimum value for current data.         !
!               Finfo( 9,*) => maximum value for current data.         !
!               Finfo(10,*) => value of scale_factor attribute if any. !
!  Fscale     Scale to convert input data to model units.              !
!  Fpoint     Latest two-time records of input point data.             !
!  Tintrp     Time (sec) of latest field snapshots used for            !
!               interpolation.                                         !
!  Vtime      Latest two-time values of processed input data.          !
!                                                                      !
!=======================================================================
!
        USE mod_param

        implicit none
!
!  Maximum number of variables in a generic NetCDF file (MV) and
!  maximum number of variables in information arrays (NV).
!
        integer, parameter :: MV = 800
        integer, parameter :: NV = 800
!
!  Input/output grid-type variables.
!
        integer, parameter :: p2dvar = 1         ! 2D PSI-variable
        integer, parameter :: r2dvar = 2         ! 2D RHO-variable
        integer, parameter :: u2dvar = 3         ! 2D U-variable
        integer, parameter :: v2dvar = 4         ! 2D V-variable
        integer, parameter :: p3dvar = 5         ! 3D PSI-variable
        integer, parameter :: r3dvar = 6         ! 3D RHO-variable
        integer, parameter :: u3dvar = 7         ! 3D U-variable
        integer, parameter :: v3dvar = 8         ! 3D V-variable
        integer, parameter :: w3dvar = 9         ! 3D W-variable
        integer, parameter :: b3dvar = 10        ! 3D BED-sediment
!
!  Number of horizontal interior and boundary water points.
!
        integer  :: Nxyp(Ngrids)        ! PSI water points
        integer  :: Nxyr(Ngrids)        ! RHO water points
        integer  :: Nxyu(Ngrids)        ! U water points
        integer  :: Nxyv(Ngrids)        ! V water points
!
!  Number of horizontal interior water points only.
!
        integer  :: NwaterR(Ngrids)     ! RHO water points
        integer  :: NwaterU(Ngrids)     ! U water points
        integer  :: NwaterV(Ngrids)     ! V water points
!
!  Lower and upper bound ranges for RHO-type variables for processing
!  state vector and observations.
!
        integer, dimension(Ngrids) :: rILB
        integer, dimension(Ngrids) :: rIUB
        integer, dimension(Ngrids) :: rJLB
        integer, dimension(Ngrids) :: rJUB
!
        real(r8), dimension(Ngrids) :: rXmin
        real(r8), dimension(Ngrids) :: rXmax
        real(r8), dimension(Ngrids) :: rYmin
        real(r8), dimension(Ngrids) :: rYmax
!
!  Lower and upper bound ranges for U-type variables for processing
!  state vector and observations.
!
        integer, dimension(Ngrids) :: uILB
        integer, dimension(Ngrids) :: uIUB
        integer, dimension(Ngrids) :: uJLB
        integer, dimension(Ngrids) :: uJUB
!
        real(r8), dimension(Ngrids) :: uXmin
        real(r8), dimension(Ngrids) :: uXmax
        real(r8), dimension(Ngrids) :: uYmin
        real(r8), dimension(Ngrids) :: uYmax
!
!  Lower and upper bound ranges for V-type variables for processing
!  state vector and observations.
!
        integer, dimension(Ngrids) :: vILB
        integer, dimension(Ngrids) :: vIUB
        integer, dimension(Ngrids) :: vJLB
        integer, dimension(Ngrids) :: vJUB
!
        real(r8), dimension(Ngrids) :: vXmin
        real(r8), dimension(Ngrids) :: vXmax
        real(r8), dimension(Ngrids) :: vYmin
        real(r8), dimension(Ngrids) :: vYmax
!
!  Switches indicating which variables are written to output files.
!
        logical  :: Hout(NV,Ngrids) ! history file switches
        logical  :: Sout(NV,Ngrids) ! station file switches
!
!  Grid identification indices.
!
        integer  :: idXgrd = -1   ! XI-grid position
        integer  :: idYgrd = -2   ! ETA-grid position
        integer  :: idZgrd = -3   ! S-grid position
        integer  :: iddpth = -4   ! depth
        integer  :: idglon = -5   ! longitude
        integer  :: idglat = -6   ! latitude
!
!  Input/output identification indices.
!
        integer  :: idAice        ! fraction of cell covered by ice
        integer  :: idbath        ! bathymetry
        integer  :: idCfra        ! cloud fraction
        integer  :: idCosW        ! COS(omega(k)*t)
        integer  :: idCos2        ! COS(omega(k)*t)*COS(omega(l)*t)
        integer  :: idDano        ! density anomaly
        integer  :: idDiff(2)     ! temperature and salinity diffusion
        integer  :: iddQdT        ! heat flux sensitivity to SST
        integer  :: idEmPf        ! E-P from bulk_flux.F
        integer  :: idevap        ! evaporation rate
        integer  :: idFsur        ! free-surface
        integer  :: idFsuD        ! detided free-surface
        integer  :: idFsuH        ! free-surface tide harmonics
        integer  :: idGhat(2)     ! KPP nonlocal transport
        integer  :: idHbbl        ! depth of bottom boundary layer
        integer  :: idHice        ! depth of ice cover
        integer  :: idHsbl        ! depth of surface boundary layer
        integer  :: idHsno        ! depth of snow cover
        integer  :: idKhor        ! convolution horizontal diffusion
        integer  :: idKver        ! convolution vertical diffusion
        integer  :: idLdwn        ! downwelling longwave radiation flux
        integer  :: idLhea        ! net latent heat flux
        integer  :: idLrad        ! net longwave radiation flux
        integer  :: idMadH        ! ADM interpolation weights
        integer  :: idMOMi        ! Initial model-observation misfit
        integer  :: idMOMf        ! final model-observation misfit
        integer  :: idMtke        ! turbulent kinetic energy
        integer  :: idMtls        ! turbulent length scale
        integer  :: idNLmi        ! initial NLM at observation locations
        integer  :: idNLmo        ! NLM at observations locations
        integer  :: idNobs        ! number of observations
        integer  :: idObsD        ! observations depth
        integer  :: idObsS        ! observations screening scale
        integer  :: idObsT        ! observations time
        integer  :: idObsX        ! observations X-grid location
        integer  :: idObsY        ! observations Y-grid location
        integer  :: idObsZ        ! observations Z-grid location
        integer  :: idOday        ! observations survey time
        integer  :: idOerr        ! observations error
        integer  :: idOtyp        ! observations type
        integer  :: idOval        ! observations value
        integer  :: idOvar        ! observations global variance
        integer  :: idOvel        ! omega vertical velocity
        integer  :: idQair        ! surface air humidity
        integer  :: idPair        ! surface air pressure
        integer  :: idPbar        ! streamfunction
        integer  :: idRdir        ! river runoff direction
        integer  :: idRepo        ! river runoff ETA-positions
        integer  :: idRflg        ! river runoff flag
        integer  :: idRtra        ! river runoff mass transport
        integer  :: idRuct        ! RHS of U-momentum coupling term
        integer  :: idRu2d        ! RHS of 2D U-momentum
        integer  :: idRu3d        ! RHS of total U-momentum
        integer  :: idRvct        ! RHS of V-momentum coupling term
        integer  :: idRv2d        ! RHS of 2D V-momentum
        integer  :: idRv3d        ! RHS of total V-momentum
        integer  :: idRxpo        ! river runoff XI-positions
        integer  :: idRvsh        ! river runoff transport profile
        integer  :: idRwet        ! wet/dry mask on RHO-points
        integer  :: idRzet        ! RHS of free-surface
        integer  :: idrain        ! rainfall rate
        integer  :: idSdif        ! vertical S-diffusion coefficient
        integer  :: idSinW        ! SIN(omega(k)*t)
        integer  :: idSin2        ! SIN(omega(k)*t)*SIN(omega(l)*t)
        integer  :: idSrad        ! net shortwave radiation flux
        integer  :: idSSHc        ! SSH climatology
        integer  :: idSSHe        ! SSH error variance
        integer  :: idSSHo        ! SSH observations
        integer  :: idSSSc        ! SSS climatology
        integer  :: idSSSf        ! SSS flux correction
        integer  :: idSSTc        ! SST climatology
        integer  :: idSSTe        ! SST error variance
        integer  :: idSSTo        ! SST observations
        integer  :: idShea        ! net sensible heat flux
        integer  :: idSWCW        ! SIN(omega(k)*t)*COS(omega(l)*t)
        integer  :: idsfwf        ! surface freswater flux
        integer  :: idTLmo        ! TLM at observation locations
        integer  :: idTair        ! surface air temperature
        integer  :: idTdif        ! vertical T-diffusion coefficient
        integer  :: idTice        ! temperature of ice surface
        integer  :: idtime        ! ocean time
        integer  :: idTpam        ! tidal potential amplitude
        integer  :: idTper        ! tidal period
        integer  :: idTpph        ! tidal potential phase
        integer  :: idTvan        ! tidal current angle
        integer  :: idTvma        ! maximum tidal current
        integer  :: idTvmi        ! minimum tidal current
        integer  :: idTvph        ! tidal current phase
        integer  :: idTzam        ! tidal elevation amplitude
        integer  :: idTzph        ! tidal elevation phase
        integer  :: idu2dA        ! accumulated 2D U-velocity
        integer  :: idU2rs        ! 2D total U-radiation stress
        integer  :: idU3rs        ! 3D total U-radiation stress
        integer  :: idU2Sd        ! 2D U-Stokes drift velocity
        integer  :: idU3Sd        ! 3D U-Stokes drift velocity
        integer  :: idUads        ! 3D U-velocity adjoint sensitivity
        integer  :: idUair        ! surface U-wind
        integer  :: idUbar        ! 2D U-velocity
        integer  :: idUwav        ! 2D U-velocity Kirby and Chen
        integer  :: idUbas        ! 2D U-velocity adjoint sensitivity
        integer  :: idUbcl        ! 2D U-velocity climatology
        integer  :: idUbcs        ! bottom max U-momentum-wave stress
        integer  :: idUVwc        ! bottom max wave-current stress magnitude
        integer  :: idUbed        ! bed load U-direction
        integer  :: idUbms        ! bottom U-momentum stress
        integer  :: idUbot        ! bed wave orbital U-velocity
        integer  :: idUbrs        ! bottom U-current stress
        integer  :: idUbtf        ! 2D U-velocity impulse forcing
        integer  :: idUbur        ! bottom U-velocity above bed
        integer  :: idUbws        ! bottom U-wave stress
        integer  :: idUclm        ! 3D U-velocity climatology
        integer  :: idUfx1        ! time averaged U-flux for 2D
        integer  :: idUfx2        ! time averaged U-flux for 3D
        integer  :: idUice        ! ice U-velocity
        integer  :: idUobs        ! 3D U-velocity observations
        integer  :: idUsms        ! surface U-momentum stress
        integer  :: idUsur        ! surface U-velocity observations
        integer  :: idUtlf        ! 3D U-velocity impulse forcing
        integer  :: idUVer        ! 3D velocity error variance
        integer  :: idUVse        ! surface velocity error variance
        integer  :: idUvel        ! 3D U-velocity
        integer  :: idUwet        ! wet/dry mask on U-points
        integer  :: idu2dD        ! detided 2D U-velocity
        integer  :: idu2dH        ! 2D U-velocity tide harmonics
        integer  :: idu3dD        ! detided 3D U-velocity
        integer  :: idu3dH        ! 3D U-velocity tide harmonics
        integer  :: idV2rs        ! 2D total V-radiation stress
        integer  :: idV3rs        ! 3D total V-radiation stress
        integer  :: idV2Sd        ! 2D U-Stokes drift velocity
        integer  :: idV3Sd        ! 3D U-Stokes drift velocity
        integer  :: idVads        ! 3D V-velocity adjoint sensitivity
        integer  :: idVair        ! surface V-wind
        integer  :: idVbar        ! 2D V-velocity
        integer  :: idVwav        ! 2D V-velocity Kirby and Chen
        integer  :: idVbas        ! 2D V-velocity adjoint sensitivity
        integer  :: idVbcl        ! 2D V-velocity climatology
        integer  :: idVbcs        ! bottom max V-current-wave stress
        integer  :: idVbed        ! bed load V-direction
        integer  :: idVbms        ! bottom V-momentum stress
        integer  :: idVbot        ! bed wave orbital V-velocity
        integer  :: idVbrs        ! bottom V-current stress
        integer  :: idVbtf        ! 2D V-velocity impulse forcing
        integer  :: idVbvr        ! bottom V-velocity above bed
        integer  :: idVbws        ! bottom V-wave stress
        integer  :: idVclm        ! 3D V-velocity climatology
        integer  :: idVfx1        ! 2D momentum time-averaged V-flux
        integer  :: idVfx2        ! 3D momentum time-averaged V-flux
        integer  :: idVice        ! ice V-velocity
        integer  :: idVmLS        ! vertical mixing length scale
        integer  :: idVmKK        ! Kinetic energy vertical mixing
        integer  :: idVmKP        ! Length scale vertical mixing
        integer  :: idVobs        ! 3D V-velocity observations
        integer  :: idVsms        ! surface V-momentum stress
        integer  :: idVsur        ! surface V-velocity observations
        integer  :: idVtlf        ! 3D V-velocity impulse forcing
        integer  :: idVvel        ! 3D V-velocity
        integer  :: idVvis        ! vertical viscosity coefficient
        integer  :: idVwet        ! wet/dry mask on V-points
        integer  :: idv2dD        ! detided 2D U-velocity
        integer  :: idv2dH        ! 2D U-velocity tide harmonics
        integer  :: idv3dD        ! detided 3D U-velocity
        integer  :: idv3dH        ! 3D U-velocity tide harmonics
        integer  :: idW2xx        ! 2D radiation stress, Sxx-component
        integer  :: idW2xy        ! 2D radiation stress, Sxy-component
        integer  :: idW2yy        ! 2D radiation stress, Syy-component
        integer  :: idW3xx        ! 3D radiation stress, Sxx-component
        integer  :: idW3xy        ! 3D radiation stress, Sxy-component
        integer  :: idW3yy        ! 3D radiation stress, Syy-component
        integer  :: idW3zx        ! 3D radiation stress, Szx-component
        integer  :: idW3zy        ! 3D radiation stress, Szy-component
        integer  :: idW3Sd        ! 3D W-Stokes drift velocity
        integer  :: idW3St        ! 3D W-Stokes drift velocity test
        integer  :: idWamp        ! wave amplitude
        integer  :: idWbeh        ! WEC Bernoulli head
        integer  :: idWbrk        ! percent wave breaking
        integer  :: idWdif        ! wave dissipation from bottom friction 
        integer  :: idWdib        ! wave dissipation from surface breaking
        integer  :: idWdiw        ! wave dissipation from white capping
        integer  :: idWdir        ! wave direction
        integer  :: idWdis        ! wave roller dissipation
        integer  :: idWlen        ! mean wave length
        integer  :: idWlep        ! peak wave length
        integer  :: idWorb        ! wave orbital velocity
        integer  :: idWptp        ! surface wave period
        integer  :: idWpbt        ! bottom wave period
        integer  :: idWqsp        ! WEC quasi-static pressure
        integer  :: idWrol        ! wave roller action density
        integer  :: idWvel        ! true vertical velocity
        integer  :: idWztw        ! WEC quasi-static sea level adjustment
        integer  :: idZads        ! Free-surface adjoint sensitivity
        integer  :: idZavg        ! avg of free-surface
        integer  :: idZtlf        ! Free-surface impulse forcing
        integer  :: id2dPV        ! 2D potential vorticity
        integer  :: id2dRV        ! 2D relative vorticity
        integer  :: id3dPV        ! 3D potential vorticity
        integer  :: id3dRV        ! 3D relative vorticity
!
!  Input/output identification tracer indices.
!
        integer  :: idRunoff       ! Surface freshwater runoff rate
        integer  :: idIcec         ! observed (NCEP) ice concentration
        integer  :: idSkt          ! observed (NCEP) skin temperature
        integer  :: idUnms         ! surface U-momentum stress (NCEP)
        integer  :: idVnms         ! surface V-momentum stress (NCEP)

        integer  :: idTimid        ! interior ice temperature
        integer  :: idTauiw        ! ice-water friction velocity
        integer  :: idChuiw        ! ice-water momentum transfer coefficient
        integer  :: idS0mk         ! salinity of molecular sub-layer under ice
        integer  :: idT0mk         ! temperature of molecular sub-layer under ice
        integer  :: idWfr          ! frazil ice growth rate
        integer  :: idWai          ! ice growth/melt rate
        integer  :: idWao          ! ice growth/melt rate
        integer  :: idWio          ! ice growth/melt rate
        integer  :: idWro          ! ice pond runoff
        integer  :: idSfwat        ! surface melt water thickness on ice
        integer  :: idAgeice       ! age of sea ice (time since formation)
        integer  :: idIomflx       ! ice-ocean mass flux
        integer  :: idWg2d         ! wind gustiness from NCEP
        integer  :: idCdd          ! momentum transfer coefficient from NCEP
        integer  :: idChd          ! sensible heat trans. coef. from NCEP
        integer  :: idCed          ! latent heat transfer coef. from NCEP
        integer  :: idWg2m         ! wind gustiness from model
        integer  :: idCdm          ! momentum transfer coefficient from model
        integer  :: idChm          ! sensible heat trans. coef. from model
        integer  :: idCem          ! latent heat transfer coef. from model
        integer  :: idRhoa         ! near-surface air density from NCEP
        integer  :: idSig11        ! internal ice stress component 11
        integer  :: idSig22        ! internal ice stress component 22
        integer  :: idSig12        ! internal ice stress component 12

        integer, allocatable :: idRtrc(:)    ! river runoff for tracers
        integer, allocatable :: idTads(:)    ! tracers adjoint sentivity
        integer, allocatable :: idTbot(:)    ! bottom flux for tracers
        integer, allocatable :: idTbry(:,:)  ! tracers boundary
        integer, allocatable :: idTclm(:)    ! tracers climatology
        integer, allocatable :: idTerr(:)    ! tracers error variance
        integer, allocatable :: idTobs(:)    ! tracers observations
        integer, allocatable :: idTsur(:)    ! surface flux for tracers
        integer, allocatable :: idTtlf(:)    ! tracers impulse forcing
!
!  Boundary conditions identification indices.
!
        integer  :: idU2bc(4)      ! 2D U-velocity boundary conditions
        integer  :: idU3bc(4)      ! 3D U-velocity boundary conditions
        integer  :: idV2bc(4)      ! 2D V-velocity boundary conditions
        integer  :: idV3bc(4)      ! 3D V-velocity boundary conditions
        integer  :: idZbry(4)      ! free-surface boundary conditions
#ifdef INWAVE_COUPLING
        integer  :: idACbc(4)      ! wave energy boundary condition
        integer  :: idTa           ! absolute peak period
        integer  :: idACen         ! wave energy
        integer  :: idACcx         ! wave energy celerity in xi dir
        integer  :: idACcy         ! wave energy celerity in etai dir
        integer  :: idACct         ! wave energy celerity in theta dir
#endif
!
!  Time-averaged quadratic terms IDs.
!
        integer  :: idU2av                    ! <ubar*ubar>
        integer  :: idV2av                    ! <vbar*vbar>
        integer  :: idZZav                    ! <zeta*zeta>

        integer  :: idHUav                    ! <Huon>
        integer  :: idHVav                    ! <Hvom>
        integer  :: idUUav                    ! <u*u>
        integer  :: idUVav                    ! <u*v>
        integer  :: idVVav                    ! <v*v>

        integer, allocatable :: iHUTav(:)     ! <Huon*t> for active tracers
        integer, allocatable :: iHVTav(:)     ! <Hvom*t> for active tracers
        integer, allocatable :: idTTav(:)     ! <t*t> for active tracers
        integer, allocatable :: idUTav(:)     ! <u*t> for active tracers
        integer, allocatable :: idVTav(:)     ! <v*t> for active tracers

#ifdef ICE_MODEL
 !
 !  Ice variable IDs.
 !
         integer  :: idUibc(4)      ! ice U-velocity boundary conditions
         integer  :: idVibc(4)      ! ice V-velocity boundary conditions
         integer  :: idAibc(4)      ! ice concentration boundary conditions
         integer  :: idHibc(4)      ! ice thickness boundary conditions
         integer  :: idHsnbc(4)     ! snow thickness boundary conditions
         integer  :: idTibc(4)      ! ice interior temp. boundary conditions
         integer  :: idSfwbc(4)     ! surface melt water boundary conditions
         integer  :: idAgeicebc(4)  ! ice age boundary conditions
         integer  :: idS11bc(4)     ! ice stress Sig11 boundary conditions
         integer  :: idS22bc(4)     ! ice stress Sig22 boundary conditions
         integer  :: idS12bc(4)     ! ice stress Sig12 boundary conditions
#endif

#ifdef DIAGNOSTICS
!
!  Tracer/Momentum Diagnostic variable IDs.
!
        integer, allocatable :: idDtrc(:,:)   ! tracers terms
        integer, allocatable :: idDu2d(:)     ! 2D u-momentum terms
        integer, allocatable :: idDv2d(:)     ! 2D v-momentum terms
        integer, allocatable :: idDu3d(:)     ! 3D u-momentum terms
        integer, allocatable :: idDv3d(:)     ! 3D v-momentum terms
#endif
!
!  Assimilation state variables indices (order is important).
!
        integer  :: isFsur = 1                ! free-surface
        integer  :: isUbar = 2                ! 2D U-velocity
        integer  :: isVbar = 3                ! 2D V-velocity
        integer  :: isUvel = 4                ! 3D U-velocity
        integer  :: isVvel = 5                ! 3D V-velocity
        integer  :: isUstr                    ! surface u-stress
        integer  :: isVstr                    ! surface v-stress
        integer, allocatable :: isTsur(:)     ! surface tracer flux
        integer, allocatable :: isTvar(:)     ! tracers
        integer, allocatable :: idSvar(:)     ! state vector indices
        integer, allocatable :: idSbry(:)     ! state boundaries indices
!
!  Input/Output NetCDF files IDs.
!
        integer  :: ncADJid(Ngrids)     ! input/output adjoint
        integer  :: ncADSid(Ngrids)     ! input adjoint sensitivity
        integer  :: ncAVGid(Ngrids)     ! output averages
        integer  :: ncBLKid(Ngrids)     ! Basic state surface fluxes
        integer  :: ncBRYid(Ngrids)     ! input boundary conditions
#ifdef INWAVE_COUPLING
        integer  :: ncIWBRYid(Ngrids)     ! input boundary conditions for Ac
        integer  :: ncIWINIid(Ngrids)     ! input initial conditions for Ac,cx,cy,ct
#endif
        integer  :: ncCLMid(Ngrids)     ! input climatology
        integer  :: ncDIAid(Ngrids)     ! output diagnostics
        integer  :: ncFLTid(Ngrids)     ! output floats
        integer  :: ncERRid(Ngrids)     ! output 4DVar posterior error
        integer  :: ncFRCid(NV,Ngrids)  ! input forcing
        integer  :: ncFWDid(Ngrids)     ! forward solution
        integer  :: ncGRDid(Ngrids)     ! input grid
        integer  :: ncGSTid(Ngrids)     ! input/output GST restart
        integer  :: ncHISid(Ngrids)     ! output history
        integer  :: ncHSSid(Ngrids)     ! input/output Hessian vectors
        integer  :: ncINIid(Ngrids)     ! input/output NLM initial
        integer  :: ncIRPid(Ngrids)     ! input/output RPM initial
        integer  :: ncITLid(Ngrids)     ! input/output TLM initial
        integer  :: ncLCZid(Ngrids)     ! input/output Lanczos vectors
        integer  :: ncMODid(Ngrids)     ! output 4DVAR fields
        integer  :: ncNRMid(4,Ngrids)   ! input/output covariance norm
        integer  :: ncOBSid(Ngrids)     ! input/output observations
        integer  :: ncRSTid(Ngrids)     ! input/output restart
        integer  :: ncSSHid(Ngrids)     ! SSH observations
        integer  :: ncSSTid(Ngrids)     ! SST observations
        integer  :: ncSTAid(Ngrids)     ! output stations
        integer  :: ncTIDEid(Ngrids)    ! input/output tide forcing
        integer  :: ncTLFid(Ngrids)     ! input/output TLM/RPM impulses
        integer  :: ncTLMid(Ngrids)     ! input/output tangent linear
        integer  :: ncTOBSid(Ngrids)    ! tracer observations
        integer  :: ncVOBSid(Ngrids)    ! currents observations
        integer  :: ncVSURid(Ngrids)    ! surface currents observations

        integer  :: idefADJ(Ngrids)     ! adjoint file creation flag
        integer  :: idefAVG(Ngrids)     ! averages file creation flag
        integer  :: idefDIA(Ngrids)     ! diagnostics file creation flag
        integer  :: idefHIS(Ngrids)     ! history file creation flag
        integer  :: idefTLM(Ngrids)     ! tangent file creation flag
!
!  Output NetCDF variables IDs.
!
        integer, allocatable :: idTvar(:)     ! tracers variables
        integer, allocatable :: adjTid(:,:)   ! adjoint tracers IDs
        integer, allocatable :: avgTid(:,:)   ! averages tracers IDs
        integer, allocatable :: errTid(:,:)   ! error tracer IDs
        integer, allocatable :: fltTid(:,:)   ! floats tracers IDs
        integer, allocatable :: hisTid(:,:)   ! history tracers IDs
        integer, allocatable :: hssTid(:,:)   ! Hessian tracers IDs
        integer, allocatable :: iniTid(:,:)   ! initial NLM tracers IDs
        integer, allocatable :: irpTid(:,:)   ! initial RPM tracers IDs
        integer, allocatable :: itlTid(:,:)   ! initial TLM tracers IDs
        integer, allocatable :: lczTid(:,:)   ! Lanczos tracers IDs
        integer, allocatable :: obsTid(:,:)   ! observations tracers IDs
        integer, allocatable :: rstTid(:,:)   ! restart tracers IDs
        integer, allocatable :: staTid(:,:)   ! stations tracers IDs
        integer, allocatable :: tlfTid(:,:)   ! TLM impulse tracers IDs
        integer, allocatable :: tlmTid(:,:)   ! tangent tracers IDs

        integer  :: adjVid(NV,Ngrids)    ! adjoint variables IDs
        integer  :: avgVid(NV,Ngrids)    ! averages variables IDs
        integer  :: diaVid(NV,Ngrids)    ! diagnostics variables IDs
        integer  :: errVid(NV,Ngrids)    ! error variables IDs
        integer  :: fltVid(-6:NV,Ngrids) ! floats variables IDs
        integer  :: hisVid(NV,Ngrids)    ! history variables IDs
        integer  :: hssVid(NV,Ngrids)    ! Hessian variables IDs
        integer  :: iniVid(NV,Ngrids)    ! initial NLM variables IDs
        integer  :: irpVid(NV,Ngrids)    ! initial RPM variables IDs
        integer  :: itlVid(NV,Ngrids)    ! initial TLM variables IDs
        integer  :: lczVid(NV,Ngrids)    ! Lanczos variables IDs
        integer  :: modVid(NV,Ngrids)    ! 4DVAR variables IDs
        integer  :: nrmVid(4,NV,Ngrids)  ! norm variables IDs
        integer  :: obsVid(NV,Ngrids)    ! observations variables IDs
        integer  :: rstVid(NV,Ngrids)    ! restart variables IDs
        integer  :: staVid(NV,Ngrids)    ! stations variables IDs
        integer  :: tideVid(NV,Ngrids)   ! tide variables IDs
        integer  :: tlfVid(NV,Ngrids)    ! TLM impulse variables IDs
        integer  :: tlmVid(NV,Ngrids)    ! tangent variables IDs

        integer  :: tADJindx(Ngrids)     ! adjoint time record index
        integer  :: tAVGindx(Ngrids)     ! averages time record index
        integer  :: tDIAindx(Ngrids)     ! diagnostics time record index
        integer  :: tERRindx(Ngrids)     ! error time record index
        integer  :: tFLTindx(Ngrids)     ! floats time record index
        integer  :: tHISindx(Ngrids)     ! history time record index
        integer  :: tHSSindx(Ngrids)     ! Hessian time record index
        integer  :: tINIindx(Ngrids)     ! initial NLM time record index
        integer  :: tIRPindx(Ngrids)     ! initial RPM time record index
        integer  :: tITLindx(Ngrids)     ! initial TLM time record index
        integer  :: tLCZindx(Ngrids)     ! Lanczos time record index
        integer  :: tNRMindx(4,Ngrids)   ! norm time record index
        integer  :: tRSTindx(Ngrids)     ! restart time record index
        integer  :: tSTAindx(Ngrids)     ! stations time record index
        integer  :: tTLFindx(Ngrids)     ! TLM impulse time record index
        integer  :: tTLMindx(Ngrids)     ! tangent time record index
!
!  Input/Output information variables.
!
        logical  :: Linfo(4,NV,Ngrids)

        integer  :: Iinfo(9,NV,Ngrids)

        real(r8) :: Finfo(10,NV,Ngrids)
        real(r8) :: Fpoint(2,NV,Ngrids)
        real(r8) :: Fscale(NV,Ngrids)
        real(r8) :: Tintrp(2,NV,Ngrids)

        real(r8), allocatable :: Vtime(:,:,:)

        character (len=5  ) :: version = '3.4  '
        character (len=40 ) :: varnam(MV)
        character (len=44 ) :: date_str
        character (len=46 ) :: Tname(0:NV)
        character (len=80 ) :: Cinfo(NV,Ngrids)
        character (len=100) :: Vname(5,0:NV)
        character (len=120) :: history
!
!  Source code root directory, cpp header file and directory, and
!  analytical expression directory.
!
        character (len=80) :: Rdir
        character (len=80) :: Hdir
        character (len=80) :: Hfile
        character (len=80) :: Adir
!
!  Analyical header file logical and names.
!
        logical :: Lanafile

        character (len=200), dimension(49) :: ANANAME

#ifdef BIOLOGY
!
!  Biology models file logical and names.
!
        logical, dimension(4) :: Lbiofile

        character (len=200), dimension(4) :: BIONAME
#endif
!
!  SVN revision and repository root URL.
!
        character (len=40 ) :: svn_rev
        character (len=120) :: svn_url

      CONTAINS

      SUBROUTINE initialize_ncparam
!
!=======================================================================
!                                                                      !
!  This routine allocates and initializes all variables in module      !
!  "mod_ncparam" for all nested grids.                                 !
!                                                                      !
!=======================================================================
!
        USE mod_parallel
#ifdef BIOLOGY
        USE mod_biology
#endif
        USE mod_iounits
        USE mod_scalars
#if defined SEDIMENT || defined BBL_MODEL
        USE mod_sediment
#endif
!
!  Local variable declarations.
!
        logical :: load

        integer, parameter :: inp = 10
#ifdef DISTRIBUTE
        integer :: Itile, Jtile
#endif
        integer :: Lvar, Ntiles, i, ic, ie, is, j, ng
        integer :: gtype, tile, varid

        real(r8), parameter :: epsilon = 1.0E-8_r8
        real(r8), parameter :: spv = 0.0_r8
        real(r8) :: offset, scale

        character (len=120), dimension(7) :: Vinfo
!
!-----------------------------------------------------------------------
!  Allocate several indices variables.
!-----------------------------------------------------------------------
!
        allocate ( idRtrc(MT) )
        allocate ( idTbot(MT) )
        allocate ( idTads(MT) )
        allocate ( idTbry(4,MT) )
        allocate ( idTclm(MT) )
        allocate ( idTerr(MT) )
        allocate ( idTobs(MT) )
        allocate ( idTsur(MT) )
        allocate ( idTtlf(MT) )
        allocate ( idTvar(MT) )
        allocate ( isTvar(MT) )
        allocate ( isTsur(MT) )

        allocate ( adjTid(MT,Ngrids) )
        allocate ( avgTid(MT,Ngrids) )
        allocate ( fltTid(MT,Ngrids) )
        allocate ( errTid(MT,Ngrids) )
        allocate ( hisTid(MT,Ngrids) )
        allocate ( hssTid(MT,Ngrids) )
        allocate ( iniTid(MT,Ngrids) )
        allocate ( irpTid(MT,Ngrids) )
        allocate ( itlTid(MT,Ngrids) )
        allocate ( lczTid(MT,Ngrids) )
        allocate ( obsTid(MT,Ngrids) )
        allocate ( rstTid(MT,Ngrids) )
        allocate ( staTid(MT,Ngrids) )
        allocate ( tlfTid(MT,Ngrids) )
        allocate ( tlmTid(MT,Ngrids) )

        allocate ( iHUTav(NAT) )
        allocate ( iHVTav(NAT) )
        allocate ( idTTav(NAT) )
        allocate ( idUTav(NAT) )
        allocate ( idVTav(NAT) )

#ifdef DIAGNOSTICS
        allocate ( idDtrc(MT,NDT) )
        allocate ( idDu2d(NDM2d) )
        allocate ( idDv2d(NDM2d) )
        allocate ( idDu3d(NDM3d) )
        allocate ( idDv3d(NDM3d) )
#endif
        allocate ( idSvar(MAXVAL(NSV)+1) )
        allocate ( idSbry(MAXVAL(NSV)+1) )

#if defined SSH_TIDES || defined UV_TIDES
        allocate ( Vtime(MAX(2,MTC),NV,Ngrids) )
#else
        allocate ( Vtime(2,NV,Ngrids) )
#endif
!
!-----------------------------------------------------------------------
!  Set minimum and maximum fractional coordinates for processing
!  observations. Either the full grid or only interior points will
!  be considered.  The strategy here is to add a small value (epsilon)
!  to the eastern and northern boundary values of Xmax and Ymax so
!  observations at such boundaries locations are processed. This
!  is needed because the .lt. operator in the following conditional:
!
!     IF (...
!    &    ((Xmin.le.Xobs(iobs)).and.(Xobs(iobs).lt.Xmax)).and.          &
!    &    ((Ymin.le.Yobs(iobs)).and.(Yobs(iobs).lt.Ymax))) THEN
!-----------------------------------------------------------------------
!
!  Allocate fractional grid lower and upper bounds structure.
!
      IF (.not.allocated(DOMAIN)) THEN
        allocate ( DOMAIN(Ngrids) )
        DO ng=1,Ngrids
          Ntiles=NtileI(ng)*NtileJ(ng)-1
          allocate ( DOMAIN(ng) % Xmin_psi (0:Ntiles) )
          allocate ( DOMAIN(ng) % Xmax_psi (0:Ntiles) )
          allocate ( DOMAIN(ng) % Ymin_psi (0:Ntiles) )
          allocate ( DOMAIN(ng) % Ymax_psi (0:Ntiles) )
          allocate ( DOMAIN(ng) % Xmin_rho (0:Ntiles) )
          allocate ( DOMAIN(ng) % Xmax_rho (0:Ntiles) )
          allocate ( DOMAIN(ng) % Ymin_rho (0:Ntiles) )
          allocate ( DOMAIN(ng) % Ymax_rho (0:Ntiles) )
          allocate ( DOMAIN(ng) % Xmin_u   (0:Ntiles) )
          allocate ( DOMAIN(ng) % Xmax_u   (0:Ntiles) )
          allocate ( DOMAIN(ng) % Ymin_u   (0:Ntiles) )
          allocate ( DOMAIN(ng) % Ymax_u   (0:Ntiles) )
          allocate ( DOMAIN(ng) % Xmin_v   (0:Ntiles) )
          allocate ( DOMAIN(ng) % Xmax_v   (0:Ntiles) )
          allocate ( DOMAIN(ng) % Ymin_v   (0:Ntiles) )
          allocate ( DOMAIN(ng) % Ymax_v   (0:Ntiles) )
        END DO
      END IF
!
!  Initialize.
!
      DO ng=1,Ngrids
        DOMAIN(ng) % Xmin_psi = spv
        DOMAIN(ng) % Xmax_psi = spv
        DOMAIN(ng) % Ymin_psi = spv
        DOMAIN(ng) % Ymax_psi = spv
        DOMAIN(ng) % Xmin_rho = spv
        DOMAIN(ng) % Xmax_rho = spv
        DOMAIN(ng) % Ymin_rho = spv
        DOMAIN(ng) % Ymax_rho = spv
        DOMAIN(ng) % Xmin_u   = spv
        DOMAIN(ng) % Xmax_u   = spv
        DOMAIN(ng) % Ymin_u   = spv
        DOMAIN(ng) % Ymax_u   = spv
        DOMAIN(ng) % Xmin_v   = spv
        DOMAIN(ng) % Xmax_v   = spv
        DOMAIN(ng) % Ymin_v   = spv
        DOMAIN(ng) % Ymax_v   = spv
      END DO
!
!  Set RHO-points domain lower and upper bounds (integer).
!
      DO ng=1,Ngrids
# ifdef DISTRIBUTE
        CALL get_bounds (ng, MyRank, r2dvar, 0, Itile, Jtile,           &
     &                   rILB(ng), rIUB(ng), rJLB(ng), rJUB(ng))
#  ifndef FULL_GRID
        IF (Itile.eq.0) THEN
          rILB(ng)=rILB(ng)+1
        END IF
        IF (Itile.eq.(NtileI(ng)-1)) THEN
          rIUB(ng)=rIUB(ng)-1
        END IF
        IF (Jtile.eq.0) THEN
          rJLB(ng)=rJLB(ng)+1
        END IF
        IF (Jtile.eq.(NtileJ(ng)-1)) THEN
          rJUB(ng)=rJUB(ng)-1
        END IF
#  endif
# else
#  ifdef FULL_GRID
        rILB(ng)=0
        rIUB(ng)=Lm(ng)+1
        rJLB(ng)=0
        rJUB(ng)=Mm(ng)+1
#  else
        rILB(ng)=1
        rIUB(ng)=Lm(ng)
        rJLB(ng)=1
        rJUB(ng)=Mm(ng)
#  endif
# endif
!
!  Minimum and maximum fractional coordinates for RHO-points.
!
        DO tile=0,NtileI(ng)*NtileJ(ng)-1
          CALL get_domain (ng, tile, r2dvar, 0, epsilon,                &
# ifdef FULL_GRID
     &                     .TRUE.,                                      &
# else
     &                     .FALSE.,                                     &
# endif
     &                     DOMAIN(ng) % Xmin_rho(tile),                 &
     &                     DOMAIN(ng) % Xmax_rho(tile),                 &
     &                     DOMAIN(ng) % Ymin_rho(tile),                 &
     &                     DOMAIN(ng) % Ymax_rho(tile))
        END DO
# ifdef DISTRIBUTE
        rXmin(ng)=DOMAIN(ng)%Xmin_rho(MyRank)
        rXmax(ng)=DOMAIN(ng)%Xmax_rho(MyRank)
        rYmin(ng)=DOMAIN(ng)%Ymin_rho(MyRank)
        rYmax(ng)=DOMAIN(ng)%Ymax_rho(MyRank)
# else
        rXmin(ng)=DOMAIN(ng)%Xmin_rho(0)
        rXmax(ng)=DOMAIN(ng)%Xmax_rho(0)
        rYmin(ng)=DOMAIN(ng)%Ymin_rho(0)
        rYmax(ng)=DOMAIN(ng)%Ymax_rho(0)
# endif
      END DO
!
!  Set U-points domain lower and upper bounds (integer).
!
      DO ng=1,Ngrids
# ifdef DISTRIBUTE
        CALL get_bounds (ng, MyRank, u2dvar, 0, Itile, Jtile,           &
     &                   uILB(ng), uIUB(ng), uJLB(ng), uJUB(ng))
#  ifndef FULL_GRID
        IF (Itile.eq.0) THEN
          uILB(ng)=uILB(ng)+1
        END IF
        IF (Itile.eq.(NtileI(ng)-1)) THEN
          uIUB(ng)=uIUB(ng)-1
        END IF
        IF (Jtile.eq.0) THEN
          uJLB(ng)=uJLB(ng)+1
        END IF
        IF (Jtile.eq.(NtileJ(ng)-1)) THEN
          uJUB(ng)=uJUB(ng)-1
        END IF
#  endif
# else
#  ifdef FULL_GRID
        uILB(ng)=1
        uIUB(ng)=Lm(ng)+1
        uJLB(ng)=0
        uJUB(ng)=Mm(ng)+1
#  else
        uILB(ng)=2
        uIUB(ng)=Lm(ng)
        uJLB(ng)=1
        uJUB(ng)=Mm(ng)
#  endif
# endif
!
!  Minimum and maximum fractional coordinates for U-points.
!
        DO tile=0,NtileI(ng)*NtileJ(ng)-1
          CALL get_domain (ng, tile, u2dvar, 0, epsilon,                &
# ifdef FULL_GRID
     &                     .TRUE.,                                      &
# else
     &                     .FALSE.,                                     &
# endif
     &                     DOMAIN(ng) % Xmin_u(tile),                   &
     &                     DOMAIN(ng) % Xmax_u(tile),                   &
     &                     DOMAIN(ng) % Ymin_u(tile),                   &
     &                     DOMAIN(ng) % Ymax_u(tile))
        END DO
# ifdef DISTRIBUTE
        uXmin(ng)=DOMAIN(ng)%Xmin_u(MyRank)
        uXmax(ng)=DOMAIN(ng)%Xmax_u(MyRank)
        uYmin(ng)=DOMAIN(ng)%Ymin_u(MyRank)
        uYmax(ng)=DOMAIN(ng)%Ymax_u(MyRank)
# else
        uXmin(ng)=DOMAIN(ng)%Xmin_u(0)
        uXmax(ng)=DOMAIN(ng)%Xmax_u(0)
        uYmin(ng)=DOMAIN(ng)%Ymin_u(0)
        uYmax(ng)=DOMAIN(ng)%Ymax_u(0)
# endif
      END DO
!
!  Set V-points domain lower and upper bounds (integer).
!
      DO ng=1,Ngrids
# ifdef DISTRIBUTE
        CALL get_bounds (ng, MyRank, v2dvar, 0, Itile, Jtile,           &
     &                   vILB(ng), vIUB(ng), vJLB(ng), vJUB(ng))
#  ifndef FULL_GRID
        IF (Itile.eq.0) THEN
          vILB(ng)=vILB(ng)+1
        END IF
        IF (Itile.eq.(NtileI(ng)-1)) THEN
          vIUB(ng)=vIUB(ng)-1
        END IF
        IF (Jtile.eq.0) THEN
          vJLB(ng)=vJLB(ng)+1
        END IF
        IF (Jtile.eq.(NtileJ(ng)-1)) THEN
          vJUB(ng)=vJUB(ng)-1
        END IF
#  endif
# else
#  ifdef FULL_GRID
        vILB(ng)=0
        vIUB(ng)=Lm(ng)+1
        vJLB(ng)=1
        vJUB(ng)=Mm(ng)+1
#  else
        vILB(ng)=1
        vIUB(ng)=Lm(ng)
        vJLB(ng)=2
        vJUB(ng)=Mm(ng)
#  endif
# endif
!
!  Minimum and maximum fractional coordinates for V-points.
!
        DO tile=0,NtileI(ng)*NtileJ(ng)-1
          CALL get_domain (ng, tile, v2dvar, 0, epsilon,                &
# ifdef FULL_GRID
     &                     .TRUE.,                                      &
# else
     &                     .FALSE.,                                     &
# endif
     &                     DOMAIN(ng) % Xmin_v(tile),                   &
     &                     DOMAIN(ng) % Xmax_v(tile),                   &
     &                     DOMAIN(ng) % Ymin_v(tile),                   &
     &                     DOMAIN(ng) % Ymax_v(tile))
        END DO
# ifdef DISTRIBUTE
        vXmin(ng)=DOMAIN(ng)%Xmin_v(MyRank)
        vXmax(ng)=DOMAIN(ng)%Xmax_v(MyRank)
        vYmin(ng)=DOMAIN(ng)%Ymin_v(MyRank)
        vYmax(ng)=DOMAIN(ng)%Ymax_v(MyRank)
# else
        vXmin(ng)=DOMAIN(ng)%Xmin_v(0)
        vXmax(ng)=DOMAIN(ng)%Xmax_v(0)
        vYmin(ng)=DOMAIN(ng)%Ymin_v(0)
        vYmax(ng)=DOMAIN(ng)%Ymax_v(0)
# endif
      END DO
!
!-----------------------------------------------------------------------
!  Initialize several variables.
!-----------------------------------------------------------------------
!
!  Initialize NetCDF files IDs to closed state.
!
        DO ng=1,Ngrids
          DO i=1,NV
            ncFRCid(i,ng)=-1
          END DO
          ncADJid(ng)=-1
          ncADSid(ng)=-1
          ncAVGid(ng)=-1
          ncBLKid(ng)=-1
          ncBRYid(ng)=-1
          ncCLMid(ng)=-1
          ncDIAid(ng)=-1
          ncERRid(ng)=-1
          ncFLTid(ng)=-1
          ncGSTid(ng)=-1
          ncFWDid(ng)=-1
          ncGRDid(ng)=-1
          ncHISid(ng)=-1
          ncHSSid(ng)=-1
          ncINIid(ng)=-1
          ncIRPid(ng)=-1
          ncITLid(ng)=-1
          ncLCZid(ng)=-1
          ncMODid(ng)=-1
          ncNRMid(1:4,ng)=-1
          ncOBSid(ng)=-1
          ncRSTid(ng)=-1
          ncSSHid(ng)=-1
          ncSSTid(ng)=-1
          ncSTAid(ng)=-1
          ncTIDEid(ng)=-1
          ncTLFid(ng)=-1
          ncTLMid(ng)=-1
          ncTOBSid(ng)=-1
          ncVOBSid(ng)=-1
          ncVSURid(ng)=-1
#ifdef INWAVE_COUPLING
          ncIWBRYid(ng)=-1
          ncIWINIid(ng)=-1
#endif
          tADJindx(ng)=0
          tAVGindx(ng)=0
          tDIAindx(ng)=0
          tFLTindx(ng)=0
          tHISindx(ng)=0
          tINIindx(ng)=0
          tIRPindx(ng)=1
          tITLindx(ng)=1
          tNRMindx(1:4,ng)=0
          tRSTindx(ng)=0
          tSTAindx(ng)=0
          tTLFindx(ng)=0
          tTLMindx(ng)=0

          idefADJ(ng)=-1
          idefAVG(ng)=-1
          idefDIA(ng)=-1
          idefHIS(ng)=-1
          idefTLM(ng)=-1
        END DO
!
!  Analytical files switch and names.
!
        Lanafile=.TRUE.
        DO i=1,49
          DO j=1,200
            ANANAME(i)(j:j)=' '
          END DO
        END DO

#ifdef BIOLOGY
!
!  Biology model header names.
!
        DO i=1,4
          Lbiofile(i)=.TRUE.
          DO j=1,200
            BIONAME(i)(j:j)=' '
          END DO
        END DO
#endif
!
!  Set IDs for state some state variables.
!
#if defined FORCING_SV || defined SO_SEMI
        isUstr=1
        isVstr=2
        ic=2
        DO i=1,MT
          ic=ic+1
          isTsur(i)=ic
        END DO
#else
        ic=5
        DO i=1,MT
          ic=ic+1
          isTvar(i)=ic
        END DO
# ifdef ADJUST_WSTRESS
        isUstr=ic+1
        isVstr=ic+2
        ic=ic+2
# endif
# ifdef ADJUST_STFLUX
        DO i=1,MT
          ic=ic+1
          isTsur(i)=ic
        END DO
# endif
#endif
!
!  Initialize IO information variables.
!
        DO ng=1,Ngrids
          DO i=1,NV
            Linfo(1,i,ng)=.FALSE.
            Linfo(2,i,ng)=.FALSE.
            Linfo(3,i,ng)=.FALSE.
            Linfo(4,i,ng)=.FALSE.
            Hout(i,ng)=.FALSE.
            Sout(i,ng)=.FALSE.
            Iinfo(1,i,ng)=0
            Iinfo(2,i,ng)=-1
            Iinfo(3,i,ng)=-1
            Iinfo(4,i,ng)=0
            Iinfo(5,i,ng)=0
            Iinfo(6,i,ng)=0
            Iinfo(7,i,ng)=0
            Iinfo(8,i,ng)=0
            Iinfo(9,i,ng)=0
            Finfo(1,i,ng)=0.0_r8
            Finfo(2,i,ng)=0.0_r8
            Finfo(3,i,ng)=0.0_r8
            Finfo(5,i,ng)=0.0_r8
            Finfo(6,i,ng)=0.0_r8
            Finfo(7,i,ng)=0.0_r8
            Finfo(10,i,ng)=1.0_r8
            Fscale(i,ng)=1.0_r8
            Fpoint(1,i,ng)=0.0_r8
            Fpoint(2,i,ng)=0.0_r8
            Tintrp(1,i,ng)=0.0_r8
            Tintrp(2,i,ng)=0.0_r8
            Vtime(1,i,ng)=0.0_r8
            Vtime(2,i,ng)=0.0_r8
          END DO
        END DO
!
!  Set source code root directory, cpp header file and directory, and
!  analytical expression directory.
!
#ifdef ROOT_DIR
        Rdir=ROOT_DIR
#endif
#ifdef HEADER_DIR
        Hdir=HEADER_DIR
#endif
#ifdef ROMS_HEADER
        Hfile=HEADER
#endif
#ifdef ANALYTICAL_DIR
        Adir=ANALYTICAL_DIR
#endif
!
!-----------------------------------------------------------------------
!  Define names of variables for Input/Output NetCDF files.
!-----------------------------------------------------------------------
!
!  Open input variable information file.
!
        OPEN (inp, FILE=TRIM(varname), FORM='formatted', STATUS='old',  &
     &        ERR=10)
        GOTO 20
  10    IF (Master) WRITE(stdout,50) TRIM(varname)
        STOP
  20    CONTINUE
!
!  Read in variable information.  Ignore blank and comment [char(33)=!]
!  input lines.
!
        varid=0
        DO WHILE (.TRUE.)
          READ (inp,*,ERR=30,END=40) Vinfo(1)
          Lvar=LEN_TRIM(Vinfo(1))
!
!  Extract SVN Repository Root URL.
!
          IF ((Lvar.gt.0).and.(Vinfo(1)(1:1).eq.CHAR(36))) THEN
            is=INDEX(Vinfo(1),'https')
            ie=INDEX(Vinfo(1),'/ROMS')-1
            svn_url=Vinfo(1)(is:ie)
#ifdef SVN_REV
            svn_rev=SVN_REV
#endif
!
!  Read in other variable information.
!
          ELSE IF ((Lvar.gt.0).and.(Vinfo(1)(1:1).ne.CHAR(33))) THEN
            READ (inp,*,ERR=30) Vinfo(2)
            READ (inp,*,ERR=30) Vinfo(3)
            READ (inp,*,ERR=30) Vinfo(4)
            READ (inp,*,ERR=30) Vinfo(5)
            READ (inp,*,ERR=30) Vinfo(6)
            READ (inp,*,ERR=30) Vinfo(7)
            READ (inp,*,ERR=30) scale
!
!  Determine staggered C-grid variable.
!
            SELECT CASE (TRIM(ADJUSTL(Vinfo(7))))
              CASE ('p2dvar')
                gtype=p2dvar
              CASE ('r2dvar')
                gtype=r2dvar
              CASE ('u2dvar')
                gtype=u2dvar
              CASE ('v2dvar')
                gtype=v2dvar
              CASE ('p3dvar')
                gtype=p3dvar
              CASE ('r3dvar')
                gtype=r3dvar
              CASE ('u3dvar')
                gtype=u3dvar
              CASE ('v3dvar')
                gtype=v3dvar
              CASE ('w3dvar')
                gtype=w3dvar
              CASE ('b3dvar')
                gtype=b3dvar
              CASE DEFAULT
                gtype=0
            END SELECT
!
!  Assign identification indices.
!
            load=.TRUE.
            varid=varid+1
            SELECT CASE (TRIM(ADJUSTL(Vinfo(6))))
              CASE ('idtime')
                idtime=varid
              CASE ('idbath')
                idbath=varid
              CASE ('idFsur')
                idFsur=varid
              CASE ('idRzet')
                idRzet=varid
              CASE ('idZavg')
                idZavg=varid
              CASE ('idUbar')
                idUbar=varid
              CASE ('idUwav')
                idUwav=varid
              CASE ('idRu2d')
                idRu2d=varid
              CASE ('idVbar')
                idVbar=varid
              CASE ('idVwav')
                idVwav=varid
              CASE ('idRv2d')
                idRv2d=varid
              CASE ('idUvel')
                idUvel=varid
              CASE ('idRu3d')
                idRu3d=varid
              CASE ('idVvel')
                idVvel=varid
              CASE ('idRv3d')
                idRv3d=varid
              CASE ('idWvel')
                idWvel=varid
              CASE ('idOvel')
                idOvel=varid
              CASE ('idDano')
                idDano=varid
              CASE ('idTvar(itemp)')
                idTvar(itemp)=varid
              CASE ('idTvar(isalt)')
                idTvar(isalt)=varid
              CASE ('idUsms')
                idUsms=varid
              CASE ('idVsms')
                idVsms=varid
              CASE ('idUbms')
                idUbms=varid
              CASE ('idVbms')
                idVbms=varid
              CASE ('idUbws')
                idUbws=varid
              CASE ('idUbcs')
                idUbcs=varid
              CASE ('idVbws')
                idVbws=varid
              CASE ('idVbcs')
                idVbcs=varid
              CASE ('idUVwc')
                idUVwc=varid
              CASE ('idTsur(itemp)')
                idTsur(itemp)=varid
              CASE ('iddQdT')
                iddQdT=varid
              CASE ('idsfwf')
                idsfwf=varid
              CASE ('idTsur(isalt)')
                idTsur(isalt)=varid
              CASE ('idTbot(itemp)')
                idTbot(itemp)=varid
              CASE ('idTbot(isalt)')
                idTbot(isalt)=varid
              CASE ('idGhat(itemp)')
                idGhat(itemp)=varid
              CASE ('idGhat(isalt)')
                idGhat(isalt)=varid
              CASE ('idMtke')
                idMtke=varid
              CASE ('idMtls')
                idMtls=varid
              CASE ('idVvis')
                idVvis=varid
              CASE ('idTdif')
                idTdif=varid
                idDiff(itemp)=idTdif
              CASE ('idSdif')
                idSdif=varid
                idDiff(isalt)=idSdif
              CASE ('idVmLS')
                idVmLS=varid
              CASE ('idVmKK')
                idVmKK=varid
              CASE ('idVmKP')
                idVmKP=varid
              CASE ('idZbry(iwest)')
                idZbry(iwest)=varid
              CASE ('idZbry(ieast)')
                idZbry(ieast)=varid
              CASE ('idZbry(isouth)')
                idZbry(isouth)=varid
              CASE ('idZbry(inorth)')
                idZbry(inorth)=varid
              CASE ('idU2bc(iwest)')
                idU2bc(iwest)=varid
              CASE ('idU2bc(ieast)')
                idU2bc(ieast)=varid
              CASE ('idU2bc(isouth)')
                idU2bc(isouth)=varid
              CASE ('idU2bc(inorth)')
                idU2bc(inorth)=varid
              CASE ('idV2bc(iwest)')
                idV2bc(iwest)=varid
              CASE ('idV2bc(ieast)')
                idV2bc(ieast)=varid
              CASE ('idV2bc(isouth)')
                idV2bc(isouth)=varid
              CASE ('idV2bc(inorth)')
                idV2bc(inorth)=varid
              CASE ('idU3bc(iwest)')
                idU3bc(iwest)=varid
              CASE ('idU3bc(ieast)')
                idU3bc(ieast)=varid
              CASE ('idU3bc(isouth)')
                idU3bc(isouth)=varid
              CASE ('idU3bc(inorth)')
                idU3bc(inorth)=varid
              CASE ('idV3bc(iwest)')
                idV3bc(iwest)=varid
              CASE ('idV3bc(ieast)')
                idV3bc(ieast)=varid
              CASE ('idV3bc(isouth)')
                idV3bc(isouth)=varid
              CASE ('idV3bc(inorth)')
                idV3bc(inorth)=varid
              CASE ('idTbry(iwest,itemp)')
                idTbry(iwest,itemp)=varid
              CASE ('idTbry(ieast,itemp)')
                idTbry(ieast,itemp)=varid
              CASE ('idTbry(isouth,itemp)')
                idTbry(isouth,itemp)=varid
              CASE ('idTbry(inorth,itemp)')
                idTbry(inorth,itemp)=varid
              CASE ('idTbry(iwest,isalt)')
                idTbry(iwest,isalt)=varid
              CASE ('idTbry(ieast,isalt)')
                idTbry(ieast,isalt)=varid
              CASE ('idTbry(isouth,isalt)')
                idTbry(isouth,isalt)=varid
              CASE ('idTbry(inorth,isalt)')
                idTbry(inorth,isalt)=varid
#ifdef ICE_MODEL
              CASE ('idAibc(iwest)')
                idAibc(iwest)=varid
              CASE ('idAibc(ieast)')
                idAibc(ieast)=varid
              CASE ('idAibc(isouth)')
                idAibc(isouth)=varid
              CASE ('idAibc(inorth)')
                idAibc(inorth)=varid
              CASE ('idHibc(iwest)')
                idHibc(iwest)=varid
              CASE ('idHibc(ieast)')
                idHibc(ieast)=varid
              CASE ('idHibc(isouth)')
                idHibc(isouth)=varid
              CASE ('idHibc(inorth)')
                idHibc(inorth)=varid
              CASE ('idHsnbc(iwest)')
                idHsnbc(iwest)=varid
              CASE ('idHsnbc(ieast)')
                idHsnbc(ieast)=varid
              CASE ('idHsnbc(isouth)')
                idHsnbc(isouth)=varid
              CASE ('idHsnbc(inorth)')
                idHsnbc(inorth)=varid
              CASE ('idTibc(iwest)')
                idTibc(iwest)=varid
              CASE ('idTibc(ieast)')
                idTibc(ieast)=varid
              CASE ('idTibc(isouth)')
                idTibc(isouth)=varid
              CASE ('idTibc(inorth)')
                idTibc(inorth)=varid
              CASE ('idSfwbc(iwest)')
                idSfwbc(iwest)=varid
              CASE ('idSfwbc(ieast)')
                idSfwbc(ieast)=varid
              CASE ('idSfwbc(isouth)')
                idSfwbc(isouth)=varid
              CASE ('idSfwbc(inorth)')
                idSfwbc(inorth)=varid
              CASE ('idAgeicebc(iwest)')
                idAgeicebc(iwest)=varid
              CASE ('idAgeicebc(ieast)')
                idAgeicebc(ieast)=varid
              CASE ('idAgeicebc(isouth)')
                idAgeicebc(isouth)=varid
              CASE ('idAgeicebc(inorth)')
                idAgeicebc(inorth)=varid
              CASE ('idUibc(iwest)')
                idUibc(iwest)=varid
              CASE ('idUibc(ieast)')
                idUibc(ieast)=varid
              CASE ('idUibc(isouth)')
                idUibc(isouth)=varid
              CASE ('idUibc(inorth)')
                idUibc(inorth)=varid
              CASE ('idVibc(iwest)')
                idVibc(iwest)=varid
              CASE ('idVibc(ieast)')
                idVibc(ieast)=varid
              CASE ('idVibc(isouth)')
                idVibc(isouth)=varid
              CASE ('idVibc(inorth)')
                idVibc(inorth)=varid
              CASE ('idS11bc(iwest)')
                idS11bc(iwest)=varid
              CASE ('idS11bc(ieast)')
                idS11bc(ieast)=varid
              CASE ('idS11bc(isouth)')
                idS11bc(isouth)=varid
              CASE ('idS11bc(inorth)')
                idS11bc(inorth)=varid
              CASE ('idS22bc(iwest)')
                idS22bc(iwest)=varid
              CASE ('idS22bc(ieast)')
                idS22bc(ieast)=varid
              CASE ('idS22bc(isouth)')
                idS22bc(isouth)=varid
              CASE ('idS22bc(inorth)')
                idS22bc(inorth)=varid
              CASE ('idS12bc(iwest)')
                idS12bc(iwest)=varid
              CASE ('idS12bc(ieast)')
                idS12bc(ieast)=varid
              CASE ('idS12bc(isouth)')
                idS12bc(isouth)=varid
              CASE ('idS12bc(inorth)')
                idS12bc(inorth)=varid
#endif
              CASE ('idRwet')
                idRwet=varid
              CASE ('idUwet')
                idUwet=varid
              CASE ('idVwet')
                idVwet=varid
              CASE ('idPair')
                idPair=varid
              CASE ('idTair')
                idTair=varid
              CASE ('idQair')
                idQair=varid
              CASE ('idCfra')
                idCfra=varid
              CASE ('idSrad')
                idSrad=varid
              CASE ('idLdwn')
                idLdwn=varid
              CASE ('idLrad')
                idLrad=varid
              CASE ('idLhea')
                idLhea=varid
              CASE ('idShea')
                idShea=varid
              CASE ('idrain')
                idrain=varid
              CASE ('idEmPf')
                idEmPf=varid
              CASE ('idevap')
                idevap=varid
              CASE ('idRunoff')
                idRunoff=varid
              CASE ('idUair')
                idUair=varid
              CASE ('idVair')
                idVair=varid
              CASE ('idWamp')
                idWamp=varid
              CASE ('idWbeh')
                idWbeh=varid
              CASE ('idWbrk')
                idWbrk=varid
              CASE ('idWdif')
                idWdif=varid
              CASE ('idWdib')
                idWdib=varid
              CASE ('idWdiw')
                idWdiw=varid
              CASE ('idWdir')
                idWdir=varid
              CASE ('idWdis')
                idWdis=varid
              CASE ('idWlen')
                idWlen=varid
              CASE('idWlep')
                idWlep=varid
              CASE ('idWptp')
                idWptp=varid
              CASE ('idWpbt')
                idWpbt=varid
              CASE ('idWqsp')
                idWqsp=varid
              CASE ('idWorb')
                idWorb=varid
              CASE ('idWrol')
                idWrol=varid
              CASE ('idW2xx')
                idW2xx=varid
              CASE ('idW2xy')
                idW2xy=varid
              CASE ('idW2yy')
                idW2yy=varid
              CASE ('idW3xx')
                idW3xx=varid
              CASE ('idW3xy')
                idW3xy=varid
              CASE ('idW3yy')
                idW3yy=varid
              CASE ('idW3zx')
                idW3zx=varid
              CASE ('idW3zy')
                idW3zy=varid
              CASE ('idWztw')
                idWztw=varid
              CASE ('idU2rs')
                idU2rs=varid
              CASE ('idV2rs')
                idV2rs=varid
              CASE ('idU2Sd')
                idU2Sd=varid
              CASE ('idV2Sd')
                idV2Sd=varid
              CASE ('idU3rs')
                idU3rs=varid
              CASE ('idV3rs')
                idV3rs=varid
              CASE ('idU3Sd')
                idU3Sd=varid
              CASE ('idV3Sd')
                idV3Sd=varid
              CASE ('idW3Sd')
                idW3Sd=varid
              CASE ('idW3St')
                idW3St=varid
              CASE ('idTper')
                idTper=varid
              CASE ('idTzam')
                idTzam=varid
              CASE ('idTzph')
                idTzph=varid
              CASE ('idTvph')
                idTvph=varid
              CASE ('idTvan')
                idTvan=varid
              CASE ('idTvma')
                idTvma=varid
              CASE ('idTvmi')
                idTvmi=varid
              CASE ('idTpam')
                idTpam=varid
              CASE ('idTpph')
                idTpph=varid
              CASE ('idRxpo')
                idRxpo=varid
              CASE ('idRepo')
                idRepo=varid
              CASE ('idRdir')
                idRdir=varid
              CASE ('idRvsh')
                idRvsh=varid
              CASE ('idRtra')
                idRtra=varid
              CASE ('idRflg')
                idRflg=varid
              CASE ('idRtrc(itemp)')
                idRtrc(itemp)=varid
              CASE ('idRtrc(isalt)')
                idRtrc(isalt)=varid
#ifdef INWAVE_COUPLING
              CASE ('idACen')
                idACen=varid
              CASE ('idACcx')
                idACcx=varid
              CASE ('idACcy')
                idACcy=varid
              CASE ('idACct')
                idACct=varid
              CASE ('idACbc(iwest)')
                idACbc(iwest)=varid
              CASE ('idACbc(ieast)')
                idACbc(ieast)=varid
              CASE ('idACbc(isouth)')
                idACbc(isouth)=varid
              CASE ('idACbc(inorth)')
                idACbc(inorth)=varid
              CASE ('idTa')
                idTa=varid
#endif
              CASE ('idHsbl')
                idHsbl=varid
              CASE ('idHbbl')
                idHbbl=varid
              CASE ('idUbot')
                idUbot=varid
              CASE ('idVbot')
                idVbot=varid
              CASE ('idUbur')
                idUbur=varid
              CASE ('idVbvr')
                idVbvr=varid
              CASE ('idUbrs')
                idUbrs=varid
              CASE ('idVbrs')
                idVbrs=varid
              CASE ('idSSHc')
                idSSHc=varid
              CASE ('idUbcl')
                idUbcl=varid
              CASE ('idVbcl')
                idVbcl=varid
              CASE ('idUclm')
                idUclm=varid
              CASE ('idVclm')
                idVclm=varid
              CASE ('idSSSc')
                idSSSc=varid
              CASE ('idSSSf')
                idSSSf=varid
              CASE ('idSSTc')
                idSSTc=varid
              CASE ('idTclm(itemp)')
                idTclm(itemp)=varid
              CASE ('idTclm(isalt)')
                idTclm(isalt)=varid
#if defined AD_SENSITIVITY   || defined IS4DVAR_SENSITIVITY || \
    defined OPT_OBSERVATIONS || defined SENSITIVITY_4DVAR   || \
    defined SO_SEMI
              CASE ('idZads')
                idZads=varid
              CASE ('idUbas')
                idUbas=varid
              CASE ('idVbas')
                idVbas=varid
              CASE ('idUads')
                idUads=varid
              CASE ('idVads')
                idVads=varid
              CASE ('idTads(itemp)')
                idTads(itemp)=varid
              CASE ('idTads(isalt)')
                idTads(isalt)=varid
#endif
#ifdef WEAK_CONSTRAINT
              CASE ('idZtlf')
                idZtlf=varid
              CASE ('idUbtf')
                idUbtf=varid
              CASE ('idVbtf')
                idVbtf=varid
              CASE ('idUtlf')
                idUtlf=varid
              CASE ('idVtlf')
                idVtlf=varid
              CASE ('idTtlf(itemp)')
                idTtlf(itemp)=varid
              CASE ('idTtlf(isalt)')
                idTtlf(isalt)=varid
#endif
              CASE ('idSSHo')
                idSSHo=varid
              CASE ('idSSHe')
                idSSHe=varid
              CASE ('idUobs')
                idUobs=varid
              CASE ('idVobs')
                idVobs=varid
              CASE ('idUVer')
                idUVer=varid
              CASE ('idUsur')
                idUsur=varid
              CASE ('idVsur')
                idVsur=varid
              CASE ('idUVse')
                idUVse=varid
              CASE ('idSSTo')
                idSSTo=varid
              CASE ('idSSTe')
                idSSTe=varid
              CASE ('idTobs(itemp)')
                idTobs(itemp)=varid
              CASE ('idTerr(itemp)')
                idTerr(itemp)=varid
              CASE ('idTobs(isalt)')
                idTobs(isalt)=varid
              CASE ('idTerr(isalt)')
                idTerr(isalt)=varid
#ifdef AVERAGES_DETIDE
              CASE ('idCosW')
                idCosW=varid
              CASE ('idCos2')
                idCos2=varid
              CASE ('idSinW')
                idSinW=varid
              CASE ('idSin2')
                idSin2=varid
              CASE ('idSWCW')
                idSWCW=varid
              CASE ('idFsuD')
                idFsuD=varid
              CASE ('idFsuH')
                idFsuH=varid
              CASE ('idu2dD')
                idu2dD=varid
              CASE ('idu2dH')
                idu2dH=varid
              CASE ('idv2dD')
                idv2dD=varid
              CASE ('idv2dH')
                idv2dH=varid
# ifdef SOLVE3D
              CASE ('idu3dD')
                idu3dD=varid
              CASE ('idu3dH')
                idu3dH=varid
              CASE ('idv3dD')
                idv3dD=varid
              CASE ('idv3dH')
                idv3dH=varid
# endif
#endif
#ifdef AVERAGES_QUADRATIC
              CASE ('idU2av')
                idU2av=varid
              CASE ('idV2av')
                idV2av=varid
              CASE ('idZZav')
                idZZav=varid
#endif
#ifdef SOLVE3D
# if defined AVERAGES_QUADRATIC || defined COMPOSED_GRID

              CASE ('idTTav(itemp)')
                idTTav(itemp)=varid
              CASE ('iHUTav(itemp)')
                iHUTav(itemp)=varid
              CASE ('iHVTav(itemp)')
                iHVTav(itemp)=varid
              CASE ('idUTav(itemp)')
                idUTav(itemp)=varid
              CASE ('idVTav(itemp)')
                idVTav(itemp)=varid
              CASE ('idHUav')
                idHUav=varid
              CASE ('idHVav')
                idHVav=varid
              CASE ('idUUav')
                idUUav=varid
              CASE ('idUVav')
                idUVav=varid
              CASE ('idVVav')
                idVVav=varid
# endif
# ifdef SALINITY
#  ifdef AVERAGES_QUADRATIC
              CASE ('idTTav(isalt)')
                idTTav(isalt)=varid
              CASE ('iHUTav(isalt)')
                iHUTav(isalt)=varid
              CASE ('iHVTav(isalt)')
                iHVTav(isalt)=varid
              CASE ('idUTav(isalt)')
                idUTav(isalt)=varid
              CASE ('idVTav(isalt)')
                idVTav(isalt)=varid
#  endif
# endif
#endif
#ifdef T_PASSIVE
              CASE ('idTvar(inert(i))')
                load=.TRUE.
              CASE ('idTbry(iwest,inert(i))')
                load=.TRUE.
              CASE ('idTbry(ieast,inert(i))')
                load=.TRUE.
              CASE ('idTbry(isouth,inert(i))')
                load=.TRUE.
              CASE ('idTbry(inorth,inert(i))')
                load=.TRUE.
# ifdef TS_PSOURCE
              CASE ('idRtrc(inert(i))')
                load=.TRUE.
# endif
#endif
#ifdef DIAGNOSTICS
              CASE ('id2dPV')
                id2dPV=varid
              CASE ('id2dRV')
                id2dRV=varid
# ifdef SOLVE3D
              CASE ('id3dPV')
                id3dPV=varid
              CASE ('id3dRV')
                id3dRV=varid
# endif
#endif
#ifdef DIAGNOSTICS_UV
              CASE ('idDu2d(M2pgrd)')
                idDu2d(M2pgrd)=varid
              CASE ('idDv2d(M2pgrd)')
                idDv2d(M2pgrd)=varid
              CASE ('idDu2d(M2sstr)')
                idDu2d(M2sstr)=varid
              CASE ('idDu2d(M2bstr)')
                idDu2d(M2bstr)=varid
              CASE ('idDv2d(M2sstr)')
                idDv2d(M2sstr)=varid
              CASE ('idDv2d(M2bstr)')
                idDv2d(M2bstr)=varid
              CASE ('idDu2d(M2rate)')
                idDu2d(M2rate)=varid
              CASE ('idDv2d(M2rate)')
                idDv2d(M2rate)=varid
# ifdef UV_ADV
              CASE ('idDu2d(M2xadv)')
                idDu2d(M2xadv)=varid
              CASE ('idDu2d(M2yadv)')
                idDu2d(M2yadv)=varid
              CASE ('idDu2d(M2hadv)')
                idDu2d(M2hadv)=varid
              CASE ('idDv2d(M2xadv)')
                idDv2d(M2xadv)=varid
              CASE ('idDv2d(M2yadv)')
                idDv2d(M2yadv)=varid
              CASE ('idDv2d(M2hadv)')
                idDv2d(M2hadv)=varid
# endif
# ifdef WEC_MELLOR
              CASE ('idDu2d(M2hrad)')
                idDu2d(M2hrad)=varid
              CASE ('idDv2d(M2hrad)')
                idDv2d(M2hrad)=varid
# endif
# ifdef WEC_VF
              CASE ('idDu2d(M2hjvf)')
                idDu2d(M2hjvf)=varid
              CASE ('idDv2d(M2hjvf)')
                idDv2d(M2hjvf)=varid
              CASE ('idDu2d(M2kvrf)')
                idDu2d(M2kvrf)=varid
              CASE ('idDv2d(M2kvrf)')
                idDv2d(M2kvrf)=varid
#  ifdef UV_COR
              CASE ('idDu2d(M2fsco)')
                idDu2d(M2fsco)=varid
              CASE ('idDv2d(M2fsco)')
                idDv2d(M2fsco)=varid
#  endif
#  ifdef BOTTOM_STREAMING
              CASE ('idDu2d(M2bstm)')
                idDu2d(M2bstm)=varid
              CASE ('idDv2d(M2bstm)')
                idDv2d(M2bstm)=varid
#  endif
#  ifdef SURFACE_STREAMING
              CASE ('idDu2d(M2sstm)')
                idDu2d(M2sstm)=varid
              CASE ('idDv2d(M2sstm)')
                idDv2d(M2sstm)=varid
#  endif
              CASE ('idDu2d(M2wrol)')
                idDu2d(M2wrol)=varid
              CASE ('idDv2d(M2wrol)')
                idDv2d(M2wrol)=varid
              CASE ('idDu2d(M2wbrk)')
                idDu2d(M2wbrk)=varid
              CASE ('idDv2d(M2wbrk)')
                idDv2d(M2wbrk)=varid
              CASE ('idDu2d(M2zeta)')
                idDu2d(M2zeta)=varid
              CASE ('idDv2d(M2zeta)')
                idDv2d(M2zeta)=varid
              CASE ('idDu2d(M2zetw)')
                idDu2d(M2zetw)=varid
              CASE ('idDv2d(M2zetw)')
                idDv2d(M2zetw)=varid
              CASE ('idDu2d(M2zqsp)')
                idDu2d(M2zqsp)=varid
              CASE ('idDv2d(M2zqsp)')
                idDv2d(M2zqsp)=varid
              CASE ('idDu2d(M2zbeh)')
                idDu2d(M2zbeh)=varid
              CASE ('idDv2d(M2zbeh)')
                idDv2d(M2zbeh)=varid
# endif
# ifdef UV_COR
              CASE ('idDu2d(M2fcor)')
                idDu2d(M2fcor)=varid
              CASE ('idDv2d(M2fcor)')
                idDv2d(M2fcor)=varid
# endif
# if defined UV_VIS2 || defined UV_VIS4
              CASE ('idDu2d(M2hvis)')
                idDu2d(M2hvis)=varid
              CASE ('idDu2d(M2xvis)')
                idDu2d(M2xvis)=varid
              CASE ('idDu2d(M2yvis)')
                idDu2d(M2yvis)=varid
              CASE ('idDv2d(M2hvis)')
                idDv2d(M2hvis)=varid
              CASE ('idDv2d(M2xvis)')
                idDv2d(M2xvis)=varid
              CASE ('idDv2d(M2yvis)')
                idDv2d(M2yvis)=varid
# endif
# ifdef SOLVE3D
              CASE ('idDu3d(M3pgrd)')
                idDu3d(M3pgrd)=varid
              CASE ('idDv3d(M3pgrd)')
                idDv3d(M3pgrd)=varid
              CASE ('idDu3d(M3vvis)')
                idDu3d(M3vvis)=varid
              CASE ('idDv3d(M3vvis)')
                idDv3d(M3vvis)=varid
              CASE ('idDu3d(M3rate)')
                idDu3d(M3rate)=varid
              CASE ('idDv3d(M3rate)')
                idDv3d(M3rate)=varid
#  ifdef UV_ADV
              CASE ('idDu3d(M3xadv)')
                idDu3d(M3xadv)=varid
              CASE ('idDu3d(M3yadv)')
                idDu3d(M3yadv)=varid
              CASE ('idDu3d(M3hadv)')
                idDu3d(M3hadv)=varid
              CASE ('idDv3d(M3xadv)')
                idDv3d(M3xadv)=varid
              CASE ('idDv3d(M3yadv)')
                idDv3d(M3yadv)=varid
              CASE ('idDv3d(M3hadv)')
                idDv3d(M3hadv)=varid
              CASE ('idDu3d(M3vadv)')
                idDu3d(M3vadv)=varid
              CASE ('idDv3d(M3vadv)')
                idDv3d(M3vadv)=varid
#  endif
#  ifdef WEC_MELLOR
              CASE ('idDu3d(M3hrad)')
                idDu3d(M3hrad)=varid
              CASE ('idDv3d(M3hrad)')
                idDv3d(M3hrad)=varid
              CASE ('idDu3d(M3vrad)')
                idDu3d(M3vrad)=varid
              CASE ('idDv3d(M3vrad)')
                idDv3d(M3vrad)=varid
#  endif
#  ifdef WEC_VF
              CASE ('idDu3d(M3vjvf)')
                idDu3d(M3vjvf)=varid
              CASE ('idDv3d(M3vjvf)')
                idDv3d(M3vjvf)=varid
              CASE ('idDu3d(M3hjvf)')
                idDu3d(M3hjvf)=varid
              CASE ('idDv3d(M3hjvf)')
                idDv3d(M3hjvf)=varid
              CASE ('idDu3d(M3kvrf)')
                idDu3d(M3kvrf)=varid
              CASE ('idDv3d(M3kvrf)')
                idDv3d(M3kvrf)=varid
#   ifdef UV_COR
              CASE ('idDu3d(M3fsco)')
                idDu3d(M3fsco)=varid
              CASE ('idDv3d(M3fsco)')
                idDv3d(M3fsco)=varid
#   endif
#   ifdef BOTTOM_STREAMING
              CASE ('idDu3d(M3bstm)')
                idDu3d(M3bstm)=varid
              CASE ('idDv3d(M3bstm)')
                idDv3d(M3bstm)=varid
#   endif
#   ifdef SURFACE_STREAMING
              CASE ('idDu3d(M3sstm)')
                idDu3d(M3sstm)=varid
              CASE ('idDv3d(M3sstm)')
                idDv3d(M3sstm)=varid
#  endif
              CASE ('idDu3d(M3wrol)')
                idDu3d(M3wrol)=varid
              CASE ('idDv3d(M3wrol)')
                idDv3d(M3wrol)=varid
              CASE ('idDu3d(M3wbrk)')
                idDu3d(M3wbrk)=varid
              CASE ('idDv3d(M3wbrk)')
                idDv3d(M3wbrk)=varid
#  endif
#  ifdef UV_COR
              CASE ('idDu3d(M3fcor)')
                idDu3d(M3fcor)=varid
              CASE ('idDv3d(M3fcor)')
                idDv3d(M3fcor)=varid
#  endif
#  if defined UV_VIS2 || defined UV_VIS4
              CASE ('idDu3d(M3hvis)')
                idDu3d(M3hvis)=varid
              CASE ('idDu3d(M3xvis)')
                idDu3d(M3xvis)=varid
              CASE ('idDu3d(M3yvis)')
                idDu3d(M3yvis)=varid
              CASE ('idDv3d(M3hvis)')
                idDv3d(M3hvis)=varid
              CASE ('idDv3d(M3xvis)')
                idDv3d(M3xvis)=varid
              CASE ('idDv3d(M3yvis)')
                idDv3d(M3yvis)=varid
#  endif
# endif
#endif
#ifdef DIAGNOSTICS_TS
              CASE ('idDtrc(iTrate)')
                load=.TRUE.
              CASE ('idDtrc(iThadv)')
                load=.TRUE.
              CASE ('idDtrc(iTxadv)')
                load=.TRUE.
              CASE ('idDtrc(iTyadv)')
                load=.TRUE.
              CASE ('idDtrc(iTvadv)')
                load=.TRUE.
# if defined TS_DIF2 || defined TS_DIF4
              CASE ('idDtrc(iThdif)')
                load=.TRUE.
              CASE ('idDtrc(iTxdif)')
                load=.TRUE.
              CASE ('idDtrc(iTydif)')
                load=.TRUE.
#  if defined MIX_GEO_TS || defined MIX_ISO_TS
              CASE ('idDtrc(iTsdif)')
                load=.TRUE.
#  endif
# endif
              CASE ('idDtrc(iTvdif)')
                load=.TRUE.
#endif
#if defined FORWARD_READ || defined FORWARD_WRITE
              CASE ('idRuct')
                idRuct=varid
              CASE ('idRvct')
                idRvct=varid
              CASE ('idUfx1')
                idUfx1=varid
              CASE ('idUfx2')
                idUfx2=varid
              CASE ('idVfx1')
                idVfx1=varid
              CASE ('idVfx2')
                idVfx2=varid
#endif
#if defined FOUR_DVAR || defined VERIFICATION
              CASE ('idNLmi')
                idNLmi=varid
              CASE ('idNLmo')
                idNLmo=varid
              CASE ('idNobs')
                idNobs=varid
              CASE ('idObsD')
                idObsD=varid
              CASE ('idObsS')
                idObsS=varid
              CASE ('idObsT')
                idObsT=varid
              CASE ('idObsX')
                idObsX=varid
              CASE ('idObsY')
                idObsY=varid
              CASE ('idObsZ')
                idObsZ=varid
              CASE ('idOday')
                idOday=varid
              CASE ('idOerr')
                idOerr=varid
              CASE ('idOtyp')
                idOtyp=varid
              CASE ('idOval')
                idOval=varid
              CASE ('idOvar')
                idOvar=varid
              CASE ('idKhor')
                idKhor=varid
              CASE ('idKver')
                idKver=varid
              CASE ('idTLmo')
                idTLmo=varid
              CASE ('idMOMi')
                idMOMi=varid
              CASE ('idMOMf')
                idMOMf=varid
#endif

#ifdef NCEP_FLUXES
              CASE ('idWg2d')
                idWg2d=varid
              CASE ('idCdd')
                idCdd=varid
              CASE ('idChd')
                idChd=varid
              CASE ('idCed')
                idCed=varid
              CASE ('idWg2m')
                idWg2m=varid
              CASE ('idCdm')
                idCdm=varid
              CASE ('idChm')
                idChm=varid
              CASE ('idCem')
                idCem=varid
              CASE ('idRhoa')
                idRhoa=varid
              CASE ('idIcec')
                idIcec=varid
              CASE ('idSkt')
                idSkt=varid
              CASE ('idUnms')
                idUnms=varid
              CASE ('idVnms')
                idVnms=varid
#endif
#ifdef ICE_MODEL
              CASE ('idUice')
                idUice=varid
              CASE ('idVice')
                idVice=varid
              CASE ('idAice')
                idAice=varid
              CASE ('idHice')
                idHice=varid
              CASE ('idHsno')
                idHsno=varid
              CASE ('idTice')
                idTice=varid
              CASE ('idTimid')
                idTimid=varid
              CASE ('idTauiw')
                idTauiw=varid
              CASE ('idChuiw')
                idChuiw=varid
              CASE ('idS0mk')
                idS0mk=varid
              CASE ('idT0mk')
                idT0mk=varid
              CASE ('idWfr')
                idWfr=varid
              CASE ('idWai')
                idWai=varid
              CASE ('idWao')
                idWao=varid
              CASE ('idWio')
                idWio=varid
              CASE ('idWro')
                idWro=varid
              CASE ('idSfwat')
                idSfwat=varid
              CASE ('idAgeice')
                idAgeice=varid
              CASE ('idIomflx')
                idIomflx=varid
              CASE ('idSig11')
                idSig11=varid
              CASE ('idSig22')
                idSig22=varid
              CASE ('idSig12')
                idSig12=varid
#endif

#ifdef BIOLOGY
# if defined BIO_FENNEL
#  include <fennel_var.h>
# elif defined ECOSIM
#  include <ecosim_var.h>
# elif defined NEMURO
#  include <nemuro_var.h>
# elif defined NPZD_FRANKS
#  include <npzd_Franks_var.h>
# elif defined NPZD_POWELL
#  include <npzd_Powell_var.h>
# elif defined NPZD_IRON
#  include <npzd_iron_var.h>
# endif
#endif

#if defined SEDIMENT || defined BBL_MODEL
#  include <sediment_var.h>
#endif

              CASE DEFAULT
                load=.FALSE.
            END SELECT
!
!  Load variable data into information arrays.
!
            IF (load) THEN
              load=.FALSE.
              IF (varid.gt.MV) THEN
                WRITE (stdout,60) MV, varid
                STOP
              END IF
              DO i=1,5
                Vname(i,varid)=TRIM(ADJUSTL(Vinfo(i)))
              END DO
              DO ng=1,Ngrids
                Iinfo(1,varid,ng)=gtype
                Fscale(varid,ng)=scale
              END DO

#ifdef T_PASSIVE
!
!  Adjust information for all inert passive tracers.
!
              SELECT CASE (TRIM(ADJUSTL(Vinfo(6))))
                CASE ('idTvar(inert(i))')
                  IF (NPT.gt.0) THEN
                    varid=varid-1
                    DO i=1,NPT
                      varid=varid+1
                      idTvar(inert(i))=varid
                      DO ng=1,Ngrids
                        Fscale(varid,ng)=scale
                        Iinfo(1,varid,ng)=gtype
                      END DO
                      WRITE (Vname(1,varid),'(a,i2.2)')                 &
     &                      TRIM(ADJUSTL(Vinfo(1))), i
                      WRITE (Vname(2,varid),'(a,a,i2.2)')               &
     &                      TRIM(ADJUSTL(Vinfo(2))), ', type ', i
                      WRITE (Vname(3,varid),'(a)')                      &
     &                      TRIM(ADJUSTL(Vinfo(3)))
                      WRITE (Vname(4,varid),'(a,a)')                    &
     &                      TRIM(Vname(1,varid)), ', scalar, series'
                      WRITE (Vname(5,varid),'(a)')                      &
     &                      TRIM(ADJUSTL(Vinfo(5)))
                    END DO
                  END IF
                CASE ('idTbry(iwest,inert(i))')
                  IF (NPT.gt.0) THEN
                    varid=varid-1
                    DO i=1,NPT
                      varid=varid+1
                      idTbry(iwest,inert(i))=varid
                      DO ng=1,Ngrids
                        Fscale(varid,ng)=scale
                        Iinfo(1,varid,ng)=gtype
                      END DO
                      WRITE (Vname(1,varid),'(a,i2.2)')                 &
     &                      TRIM(ADJUSTL(Vinfo(1))), i
                      WRITE (Vname(2,varid),'(a,a,i2.2)')               &
     &                      TRIM(ADJUSTL(Vinfo(2))), ', type ', i
                      WRITE (Vname(3,varid),'(a)')                      &
     &                      TRIM(ADJUSTL(Vinfo(3)))
                      WRITE (Vname(4,varid),'(a,a)')                    &
     &                      TRIM(Vname(1,varid)), ', scalar, series'
                      WRITE (Vname(5,varid),'(a)')                      &
     &                      TRIM(ADJUSTL(Vinfo(5)))
                    END DO
                  END IF
                CASE ('idTbry(ieast,inert(i))')
                  IF (NPT.gt.0) THEN
                    varid=varid-1
                    DO i=1,NPT
                      varid=varid+1
                      idTbry(ieast,inert(i))=varid
                      DO ng=1,Ngrids
                        Fscale(varid,ng)=scale
                        Iinfo(1,varid,ng)=gtype
                      END DO
                      WRITE (Vname(1,varid),'(a,i2.2)')                 &
     &                      TRIM(ADJUSTL(Vinfo(1))), i
                      WRITE (Vname(2,varid),'(a,a,i2.2)')               &
     &                      TRIM(ADJUSTL(Vinfo(2))), ', type ', i
                      WRITE (Vname(3,varid),'(a)')                      &
     &                      TRIM(ADJUSTL(Vinfo(3)))
                      WRITE (Vname(4,varid),'(a,a)')                    &
     &                      TRIM(Vname(1,varid)), ', scalar, series'
                      WRITE (Vname(5,varid),'(a)')                      &
     &                      TRIM(ADJUSTL(Vinfo(5)))
                    END DO
                  END IF
                CASE ('idTbry(isouth,inert(i))')
                  IF (NPT.gt.0) THEN
                    varid=varid-1
                    DO i=1,NPT
                      varid=varid+1
                      idTbry(isouth,inert(i))=varid
                      DO ng=1,Ngrids
                        Fscale(varid,ng)=scale
                        Iinfo(1,varid,ng)=gtype
                      END DO
                      WRITE (Vname(1,varid),'(a,i2.2)')                 &
     &                      TRIM(ADJUSTL(Vinfo(1))), i
                      WRITE (Vname(2,varid),'(a,a,i2.2)')               &
     &                      TRIM(ADJUSTL(Vinfo(2))), ', type ', i
                      WRITE (Vname(3,varid),'(a)')                      &
     &                      TRIM(ADJUSTL(Vinfo(3)))
                      WRITE (Vname(4,varid),'(a,a)')                    &
     &                      TRIM(Vname(1,varid)), ', scalar, series'
                      WRITE (Vname(5,varid),'(a)')                      &
     &                      TRIM(ADJUSTL(Vinfo(5)))
                    END DO
                  END IF
                CASE ('idTbry(inorth,inert(i))')
                  IF (NPT.gt.0) THEN
                    varid=varid-1
                    DO i=1,NPT
                      varid=varid+1
                      idTbry(inorth,inert(i))=varid
                      DO ng=1,Ngrids
                        Fscale(varid,ng)=scale
                        Iinfo(1,varid,ng)=gtype
                      END DO
                      WRITE (Vname(1,varid),'(a,i2.2)')                 &
     &                      TRIM(ADJUSTL(Vinfo(1))), i
                      WRITE (Vname(2,varid),'(a,a,i2.2)')               &
     &                      TRIM(ADJUSTL(Vinfo(2))), ', type ', i
                      WRITE (Vname(3,varid),'(a)')                      &
     &                      TRIM(ADJUSTL(Vinfo(3)))
                      WRITE (Vname(4,varid),'(a,a)')                    &
     &                      TRIM(Vname(1,varid)), ', scalar, series'
                      WRITE (Vname(5,varid),'(a)')                      &
     &                      TRIM(ADJUSTL(Vinfo(5)))
                    END DO
                  END IF
# ifdef TS_PSOURCE
                CASE ('idRtrc(inert(i))')
                  IF (NPT.gt.0) THEN
                    varid=varid-1
                    DO i=1,NPT
                      varid=varid+1
                      idRtrc(inert(i))=varid
                      DO ng=1,Ngrids
                        Fscale(varid,ng)=scale
                        Iinfo(1,varid,ng)=gtype
                      END DO
                      WRITE (Vname(1,varid),'(a,i2.2)')                 &
     &                      TRIM(ADJUSTL(Vinfo(1))), i
                      WRITE (Vname(2,varid),'(a,a,i2.2)')               &
     &                      TRIM(ADJUSTL(Vinfo(2))), ', type ', i
                      WRITE (Vname(3,varid),'(a)')                      &
     &                      TRIM(ADJUSTL(Vinfo(3)))
                      WRITE (Vname(4,varid),'(a,a)')                    &
     &                      TRIM(Vname(1,varid)), ', scalar, series'
                      WRITE (Vname(5,varid),'(a)')                      &
     &                      TRIM(ADJUSTL(Vinfo(5)))
                    END DO
                  END IF
# endif
              END SELECT
#endif
#ifdef DIAGNOSTICS_TS
!
!  Adjust information for tracer diagnostic variables.  This needs to be
!  done last because it needs all the tracers variable names.
!

              SELECT CASE (Vinfo(1))
                CASE ('_rate')
                  varid=varid-1
                  DO i=1,MT
                    varid=varid+1
                    idDtrc(i,iTrate)=varid
                    DO ng=1,Ngrids
                      Fscale(varid,ng)=scale
                      Iinfo(1,varid,ng)=gtype
                    END DO
                    WRITE (Vname(1,varid),'(a,a)')                      &
     &                    TRIM(ADJUSTL(Vname(1,idTvar(i)))),            &
     &                    TRIM(ADJUSTL(Vinfo(1)))
                    WRITE (Vname(2,varid),'(a,", ",a)')                 &
     &                    TRIM(ADJUSTL(Vname(2,idTvar(i)))),            &
     &                    TRIM(ADJUSTL(Vinfo(2)))
                    WRITE (Vname(3,varid),'(a,1x,a)')                   &
     &                    TRIM(ADJUSTL(Vname(3,idTvar(i)))),            &
     &                    TRIM(ADJUSTL(Vinfo(3)))
                    WRITE (Vname(4,varid),'(a,a)')                      &
     &                    TRIM(Vname(1,varid)), ', scalar, series'
                    WRITE (Vname(5,varid),'(a)')                        &
     &                    TRIM(ADJUSTL(Vinfo(5)))
                  END DO
                CASE ('_hadv')
                  varid=varid-1
                  DO i=1,MT
                    varid=varid+1
                    idDtrc(i,iThadv)=varid
                    DO ng=1,Ngrids
                      Fscale(varid,ng)=scale
                      Iinfo(1,varid,ng)=gtype
                    END DO
                    WRITE (Vname(1,varid),'(a,a)')                      &
     &                    TRIM(ADJUSTL(Vname(1,idTvar(i)))),            &
     &                    TRIM(ADJUSTL(Vinfo(1)))
                    WRITE (Vname(2,varid),'(a,", ",a)')                 &
     &                    TRIM(ADJUSTL(Vname(2,idTvar(i)))),            &
     &                    TRIM(ADJUSTL(Vinfo(2)))
                    WRITE (Vname(3,varid),'(a,1x,a)')                   &
     &                    TRIM(ADJUSTL(Vname(3,idTvar(i)))),            &
     &                    TRIM(ADJUSTL(Vinfo(3)))
                    WRITE (Vname(4,varid),'(a,a)')                      &
     &                    TRIM(Vname(1,varid)), ', scalar, series'
                    WRITE (Vname(5,varid),'(a)')                        &
     &                    TRIM(ADJUSTL(Vinfo(5)))
                  END DO
                CASE ('_xadv')
                  varid=varid-1
                  DO i=1,MT
                    varid=varid+1
                    idDtrc(i,iTxadv)=varid
                    DO ng=1,Ngrids
                      Fscale(varid,ng)=scale
                      Iinfo(1,varid,ng)=gtype
                    END DO
                    WRITE (Vname(1,varid),'(a,a)')                      &
     &                    TRIM(ADJUSTL(Vname(1,idTvar(i)))),            &
     &                    TRIM(ADJUSTL(Vinfo(1)))
                    WRITE (Vname(2,varid),'(a,", ",a)')                 &
     &                    TRIM(ADJUSTL(Vname(2,idTvar(i)))),            &
     &                    TRIM(ADJUSTL(Vinfo(2)))
                    WRITE (Vname(3,varid),'(a,1x,a)')                   &
     &                    TRIM(ADJUSTL(Vname(3,idTvar(i)))),            &
     &                    TRIM(ADJUSTL(Vinfo(3)))
                    WRITE (Vname(4,varid),'(a,a)')                      &
     &                    TRIM(Vname(1,varid)), ', scalar, series'
                    WRITE (Vname(5,varid),'(a)')                        &
     &                    TRIM(ADJUSTL(Vinfo(5)))
                  END DO
                CASE ('_yadv')
                  varid=varid-1
                  DO i=1,MT
                    varid=varid+1
                    idDtrc(i,iTyadv)=varid
                    DO ng=1,Ngrids
                      Fscale(varid,ng)=scale
                      Iinfo(1,varid,ng)=gtype
                    END DO
                    WRITE (Vname(1,varid),'(a,a)')                      &
     &                    TRIM(ADJUSTL(Vname(1,idTvar(i)))),            &
     &                    TRIM(ADJUSTL(Vinfo(1)))
                    WRITE (Vname(2,varid),'(a,", ",a)')                 &
     &                    TRIM(ADJUSTL(Vname(2,idTvar(i)))),            &
     &                    TRIM(ADJUSTL(Vinfo(2)))
                    WRITE (Vname(3,varid),'(a,1x,a)')                   &
     &                    TRIM(ADJUSTL(Vname(3,idTvar(i)))),            &
     &                    TRIM(ADJUSTL(Vinfo(3)))
                    WRITE (Vname(4,varid),'(a,a)')                      &
     &                    TRIM(Vname(1,varid)), ', scalar, series'
                    WRITE (Vname(5,varid),'(a)')                        &
     &                    TRIM(ADJUSTL(Vinfo(5)))
                  END DO
                CASE ('_vadv')
                  varid=varid-1
                  DO i=1,MT
                    varid=varid+1
                    idDtrc(i,iTvadv)=varid
                    DO ng=1,Ngrids
                      Fscale(varid,ng)=scale
                      Iinfo(1,varid,ng)=gtype
                    END DO
                    WRITE (Vname(1,varid),'(a,a)')                      &
     &                    TRIM(ADJUSTL(Vname(1,idTvar(i)))),            &
     &                    TRIM(ADJUSTL(Vinfo(1)))
                    WRITE (Vname(2,varid),'(a,", ",a)')                 &
     &                    TRIM(ADJUSTL(Vname(2,idTvar(i)))),            &
     &                    TRIM(ADJUSTL(Vinfo(2)))
                    WRITE (Vname(3,varid),'(a,1x,a)')                   &
     &                    TRIM(ADJUSTL(Vname(3,idTvar(i)))),            &
     &                    TRIM(ADJUSTL(Vinfo(3)))
                    WRITE (Vname(4,varid),'(a,a)')                      &
     &                    TRIM(Vname(1,varid)), ', scalar, series'
                    WRITE (Vname(5,varid),'(a)')                        &
     &                    TRIM(ADJUSTL(Vinfo(5)))
                  END DO
# if defined TS_DIF2 || defined TS_DIF4
                CASE ('_hdiff')
                  varid=varid-1
                  DO i=1,MT
                    varid=varid+1
                    idDtrc(i,iThdif)=varid
                    DO ng=1,Ngrids
                      Fscale(varid,ng)=scale
                      Iinfo(1,varid,ng)=gtype
                    END DO
                    WRITE (Vname(1,varid),'(a,a)')                      &
     &                    TRIM(ADJUSTL(Vname(1,idTvar(i)))),            &
     &                    TRIM(ADJUSTL(Vinfo(1)))
                    WRITE (Vname(2,varid),'(a,", ",a)')                 &
     &                    TRIM(ADJUSTL(Vname(2,idTvar(i)))),            &
     &                    TRIM(ADJUSTL(Vinfo(2)))
                    WRITE (Vname(3,varid),'(a,1x,a)')                   &
     &                    TRIM(ADJUSTL(Vname(3,idTvar(i)))),            &
     &                    TRIM(ADJUSTL(Vinfo(3)))
                    WRITE (Vname(4,varid),'(a,a)')                      &
     &                    TRIM(Vname(1,varid)), ', scalar, series'
                    WRITE (Vname(5,varid),'(a)')                        &
     &                    TRIM(ADJUSTL(Vinfo(5)))
                  END DO
                CASE ('_xdiff')
                  varid=varid-1
                  DO i=1,MT
                    varid=varid+1
                    idDtrc(i,iTxdif)=varid
                    DO ng=1,Ngrids
                      Fscale(varid,ng)=scale
                      Iinfo(1,varid,ng)=gtype
                    END DO
                    WRITE (Vname(1,varid),'(a,a)')                      &
     &                    TRIM(ADJUSTL(Vname(1,idTvar(i)))),            &
     &                    TRIM(ADJUSTL(Vinfo(1)))
                    WRITE (Vname(2,varid),'(a,", ",a)')                 &
     &                    TRIM(ADJUSTL(Vname(2,idTvar(i)))),            &
     &                    TRIM(ADJUSTL(Vinfo(2)))
                    WRITE (Vname(3,varid),'(a,1x,a)')                   &
     &                    TRIM(ADJUSTL(Vname(3,idTvar(i)))),            &
     &                    TRIM(ADJUSTL(Vinfo(3)))
                    WRITE (Vname(4,varid),'(a,a)')                      &
     &                    TRIM(Vname(1,varid)), ', scalar, series'
                    WRITE (Vname(5,varid),'(a)')                        &
     &                    TRIM(ADJUSTL(Vinfo(5)))
                  END DO
                CASE ('_ydiff')
                  varid=varid-1
                  DO i=1,MT
                    varid=varid+1
                    idDtrc(i,iTydif)=varid
                    DO ng=1,Ngrids
                      Fscale(varid,ng)=scale
                      Iinfo(1,varid,ng)=gtype
                    END DO
                    WRITE (Vname(1,varid),'(a,a)')                      &
     &                    TRIM(ADJUSTL(Vname(1,idTvar(i)))),            &
     &                    TRIM(ADJUSTL(Vinfo(1)))
                    WRITE (Vname(2,varid),'(a,", ",a)')                 &
     &                    TRIM(ADJUSTL(Vname(2,idTvar(i)))),            &
     &                    TRIM(ADJUSTL(Vinfo(2)))
                    WRITE (Vname(3,varid),'(a,1x,a)')                   &
     &                    TRIM(ADJUSTL(Vname(3,idTvar(i)))),            &
     &                    TRIM(ADJUSTL(Vinfo(3)))
                    WRITE (Vname(4,varid),'(a,a)')                      &
     &                    TRIM(Vname(1,varid)), ', scalar, series'
                    WRITE (Vname(5,varid),'(a)')                        &
     &                    TRIM(ADJUSTL(Vinfo(5)))
                  END DO
#  if defined MIX_GEO_TS || defined MIX_ISO_TS
                CASE ('_sdiff')
                  varid=varid-1
                  DO i=1,MT
                    varid=varid+1
                    idDtrc(i,iTsdif)=varid
                    DO ng=1,Ngrids
                      Fscale(varid,ng)=scale
                      Iinfo(1,varid,ng)=gtype
                    END DO
                    WRITE (Vname(1,varid),'(a,a)')                      &
     &                    TRIM(ADJUSTL(Vname(1,idTvar(i)))),            &
     &                    TRIM(ADJUSTL(Vinfo(1)))
                    WRITE (Vname(2,varid),'(a,", ",a)')                 &
     &                    TRIM(ADJUSTL(Vname(2,idTvar(i)))),            &
     &                    TRIM(ADJUSTL(Vinfo(2)))
                    WRITE (Vname(3,varid),'(a,1x,a)')                   &
     &                    TRIM(ADJUSTL(Vname(3,idTvar(i)))),            &
     &                    TRIM(ADJUSTL(Vinfo(3)))
                    WRITE (Vname(4,varid),'(a,a)')                      &
     &                    TRIM(Vname(1,varid)), ', scalar, series'
                    WRITE (Vname(5,varid),'(a)')                        &
     &                    TRIM(ADJUSTL(Vinfo(5)))
                  END DO
#  endif
# endif
                CASE ('_vdiff')
                  varid=varid-1
                  DO i=1,MT
                    varid=varid+1
                    idDtrc(i,iTvdif)=varid
                    DO ng=1,Ngrids
                      Fscale(varid,ng)=scale
                      Iinfo(1,varid,ng)=gtype
                    END DO
                    WRITE (Vname(1,varid),'(a,a)')                      &
     &                    TRIM(ADJUSTL(Vname(1,idTvar(i)))),            &
     &                    TRIM(ADJUSTL(Vinfo(1)))
                    WRITE (Vname(2,varid),'(a,", ",a)')                 &
     &                    TRIM(ADJUSTL(Vname(2,idTvar(i)))),            &
     &                    TRIM(ADJUSTL(Vinfo(2)))
                    WRITE (Vname(3,varid),'(a,1x,a)')                   &
     &                    TRIM(ADJUSTL(Vname(3,idTvar(i)))),            &
     &                    TRIM(ADJUSTL(Vinfo(3)))
                    WRITE (Vname(4,varid),'(a,a)')                      &
     &                    TRIM(Vname(1,varid)), ', scalar, series'
                    WRITE (Vname(5,varid),'(a)')                        &
     &                    TRIM(ADJUSTL(Vinfo(5)))
                  END DO
              END SELECT
#endif
            ELSE
              varid=varid-1
!!            WRITE (stdout,70) TRIM(ADJUSTL(Vinfo(1))),                &
!!   &                          TRIM(ADJUSTL(Vinfo(6)))
            END IF
          END IF
        END DO
        GOTO 40
  30    WRITE (stdout,80) TRIM(ADJUSTL(Vinfo(1)))
        STOP
  40    CLOSE (inp)

#ifdef SOLVE3D
!
!-----------------------------------------------------------------------
!  Set passive tracers surface flux metadata. The variable name is the
!  same as the basic tracer but with the _sflux suffix.
!-----------------------------------------------------------------------
!
        DO i=NAT+1,MT
          varid=varid+1
          IF (varid.gt.MV) THEN
            WRITE (stdout,60) MV, varid
            STOP
          END IF
          idTsur(i)=varid
          DO ng=1,Ngrids
            Fscale(varid,ng)=1.0_r8
            Iinfo(1,varid,ng)=r2dvar
          END DO
          WRITE (Vname(1,varid),'(a,a)')                                &
     &          TRIM(ADJUSTL(Vname(1,idTvar(i)))), '_sflux'
          WRITE (Vname(2,varid),'(a,a)')                                &
     &          TRIM(ADJUSTL(Vname(2,idTvar(i)))), ', surface flux'
          WRITE (Vname(3,varid),'(a,1x,a)')                             &
     &          TRIM(ADJUSTL(Vname(3,idTvar(i)))), 'meter second-1'
          WRITE (Vname(4,varid),'(a,a)')                                &
     &          TRIM(Vname(1,varid)), ', scalar, series'
          WRITE (Vname(5,varid),'(a)')                                  &
     &          TRIM(ADJUSTL(Vname(5,idTvar(i))))
        END DO
#endif
#ifdef WEAK_CONSTRAINT
!
!-----------------------------------------------------------------------
!  Set passive tracers impulse forcing. The metadata values are the
!  same as the basic tracer but with different indices sinc they are
!  used at the same time. This allows flexibility in the time
!  interpolation.
!-----------------------------------------------------------------------
!
        DO i=NAT+1,MT
          varid=varid+1
          IF (varid.gt.MV) THEN
            WRITE (stdout,60) MV, varid
            STOP
          END IF
          idTtlf(i)=varid
          DO ng=1,Ngrids
            Fscale(varid,ng)=1.0_r8
            Iinfo(1,varid,ng)=r3dvar
          END DO
          WRITE (Vname(1,varid),'(a)')                                  &
     &          TRIM(ADJUSTL(Vname(1,idTvar(i))))
          WRITE (Vname(2,varid),'(a,a)')                                &
     &          TRIM(ADJUSTL(Vname(2,idTvar(i)))), ' impulse forcing'
          WRITE (Vname(3,varid),'(a)')                                  &
     &          TRIM(ADJUSTL(Vname(3,idTvar(i))))
          WRITE (Vname(4,varid),'(a,a)')                                &
     &          TRIM(Vname(1,varid)), ', scalar, series'
          WRITE (Vname(5,varid),'(a)')                                  &
     &          TRIM(ADJUSTL(Vname(5,idTvar(i))))
        END DO
#endif
!
!-----------------------------------------------------------------------
!  Set model state variables indices.
!-----------------------------------------------------------------------
!
        idSvar(isFsur)=idFsur
        idSvar(isUbar)=idUbar
        idSvar(isVbar)=idVbar
        ic=3
#ifdef SOLVE3D
        idSvar(isUvel)=idUvel
        idSvar(isVvel)=idVvel
        ic=ic+2
        DO i=1,MT
          idSvar(isTvar(i))=idTvar(i)
          ic=ic+1
        END DO
#endif
#ifdef ADJUST_WSTRESS
        idSvar(isUstr)=idUsms
        idSvar(isVstr)=idVsms
#endif
#ifdef ADJUST_STFLUX
        DO i=1,MT
          idSvar(isTsur(i))=idTsur(i)
        END DO
#endif
#ifdef ADJUST_BOUNDARY
!
!-----------------------------------------------------------------------
!  If adjusting open boundaries, set metadata variables.  The variable
!  name is the same as the state variable but with the _obc suffix.
!-----------------------------------------------------------------------
!
        DO i=1,ic
          varid=varid+1
          IF (varid.gt.MV) THEN
            WRITE (stdout,60) MV, varid
            STOP
          END IF
          idSbry(i)=varid
          DO ng=1,Ngrids
            Fscale(varid,ng)=1.0_r8
            Iinfo(1,varid,ng)=0
          END DO
          WRITE (Vname(1,varid),'(a,a)')                                &
     &          TRIM(ADJUSTL(Vname(1,idSvar(i)))), '_obc'
          WRITE (Vname(2,varid),'(a,a)')                                &
     &          TRIM(ADJUSTL(Vname(2,idSvar(i)))), ', open boundaries'
          WRITE (Vname(3,varid),'(a)')                                  &
     &          TRIM(ADJUSTL(Vname(3,idSvar(i))))
          WRITE (Vname(4,varid),'(a,a)')                                &
     &          TRIM(Vname(1,varid)), ', scalar, series'
          WRITE (Vname(5,varid),'(a)')                                  &
     &          TRIM(ADJUSTL(Vname(5,idSvar(i))))
        END DO
#endif
!
  50    FORMAT (/,' MOD_NCPARAM - Unable to open variable information', &
     &          ' file: ',/,15x,a,/,15x,'Default file is located in',   &
     &          ' source directory.')
  60    FORMAT (/,' MOD_NCPARAM - too small dimension ',                &
     &          'parameter, MV = ',2i5,/,15x,                           &
     &          'change file  mod_ncparam.F  and recompile.')
  70    FORMAT (/,' MOD_NCPARM - Cannot load information for ',         &
     &          'variable: ',a,/,15x,'Need CASE construct for: ',a)
  80    FORMAT (/,' MOD_NCPARM - Error while reading information ',     &
     &          'for variable: ',a)
  90    FORMAT (a,i2.2)
 100    FORMAT (a,a,i2.2)
 110    FORMAT (a)
 120    FORMAT (a,a)

        RETURN
        END SUBROUTINE initialize_ncparam

      END MODULE mod_ncparam
