#include "cppdefs.h"
      MODULE mod_ncparam
!
!svn $Id: mod_ncparam.F 746 2008-09-10 04:09:44Z jcwarner $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
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
!               Finfo(1,*)  => Starting time (days) of data.           !
!               Finfo(2,*)  => Ending time (days) of data.             !
!               Finfo(3,*)  => Data time lower bound (days) enclosing  !
!                                model starting time.                  !
!               Finfo(4,*)  => Data time upper bound (days) enclosing  !
!                                model starting time.                  !
!               Finfo(5,*)  => length (days) of time cycling.          !
!               Finfo(6,*)  => Scale to convert time to day units.     !
!               Finfo(7,*)  => latest monotonic time (sec).            !
!               Finfo(8,*)  => minimum value for current data.         !
!               Finfo(9,*)  => maximum value for current data.         !
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
        integer, parameter :: MV = 500
        integer, parameter :: NV = 500
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
!  Lower and upper bound ranges for U-type variables for processing
!  state vector and observations.
!
        integer, dimension(Ngrids) :: uILB
        integer, dimension(Ngrids) :: uIUB
        integer, dimension(Ngrids) :: uJLB
        integer, dimension(Ngrids) :: uJUB
!
!  Lower and upper bound ranges for V-type variables for processing
!  state vector and observations.
!
        integer, dimension(Ngrids) :: vILB
        integer, dimension(Ngrids) :: vIUB
        integer, dimension(Ngrids) :: vJLB
        integer, dimension(Ngrids) :: vJUB
!
!  Switches indicating which variables are written to output files.
!
        logical  :: Hout(NV,Ngrids) ! history file switches
        logical  :: Sout(NV,Ngrids) ! station file switches
!
!  Input/output identification indices.
!
        integer  :: idXgrd = -1   ! XI-grid position
        integer  :: idYgrd = -2   ! ETA-grid position
        integer  :: idZgrd = -3   ! S-grid position
        integer  :: iddpth = -4   ! depth
        integer  :: idglon = -5   ! longitude
        integer  :: idglat = -6   ! latitude

        integer  :: idAice        ! fraction of cell covered by ice
        integer  :: idbath        ! bathymetry
        integer  :: idCfra        ! cloud fraction
        integer  :: idCG_P        ! conjugate gradient p-vector
        integer  :: idCG_R        ! conjugate gradient r-vector
        integer  :: idCG_S        ! conjugate gradient s-vector
        integer  :: idCG_V        ! conjugate gradient v-vector
        integer  :: idCG_X        ! conjugate gradient x-vector
        integer  :: idCG_Z        ! conjugate gradient z-vector
        integer  :: idCosW        ! COS(omega(k)*t)
        integer  :: idCos2        ! COS(omega(k)*t)*COS(omega(l)*t)
        integer  :: idDano        ! density anomaly
        integer  :: idDiff(2)     ! temperature and salinity diffusion
        integer  :: iddQdT        ! heat flux sensitivity to SST
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
        integer  :: idRepC        ! representer coefficients
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
        integer  :: idTper        ! tidal period
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
        integer  :: idUbur        ! bottom U-velocity above bed
        integer  :: idUbws        ! bottom U-wave stress
        integer  :: idUclm        ! 3D U-velocity climatology
        integer  :: idUfx1        ! time averaged U-flux for 2D
        integer  :: idUfx2        ! time averaged U-flux for 3D
        integer  :: idUice        ! ice U-velocity
        integer  :: idUobs        ! 3D U-velocity observations
        integer  :: idUsms        ! surface U-momentum stress
        integer  :: idUsur        ! surface U-velocity observations
        integer  :: idUVer        ! 3D velocity error variance
        integer  :: idUVse        ! surface velocity error variance
        integer  :: idUvel        ! 3D U-velocity
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
        integer  :: idVvel        ! 3D V-velocity
        integer  :: idVvis        ! vertical viscosity coefficient
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
        integer  :: idWamp        ! wind-induced wave amplitude
        integer  :: idWbrk        ! wind-induced wave breaking
        integer  :: idWdis        ! wind-induced wave dissipation
        integer  :: idWdir        ! wind-induced wave direction
        integer  :: idWlen        ! wind-induced wave length
        integer  :: idWmsk        ! wet-dry masking
        integer  :: idWptp        ! wind-induced surface wave period
        integer  :: idWpbt        ! wind-induced bottom wave period
        integer  :: idWorb        ! wind-induced wave orbital velocity
        integer  :: idWvel        ! true vertical velocity
        integer  :: idZads        ! Free-surface adjoint sensitivity
        integer  :: idZavg        ! avg of free-surface
        integer, allocatable :: idRtrc(:)    ! river runoff for tracers
        integer, allocatable :: idTads(:)    ! tracers adjoint sentivity
        integer, allocatable :: idTbot(:)    ! bottom flux for tracers
        integer, allocatable :: idTbry(:,:)  ! tracers boundary
        integer, allocatable :: idTclm(:)    ! tracers climatology
        integer, allocatable :: idTerr(:)    ! tracers error variance
        integer, allocatable :: idTobs(:)    ! tracers observations
        integer, allocatable :: idTsur(:)    ! surface flux for tracers

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

#ifdef DIAGNOSTICS
!
!  Biological 2D diagnostic variable IDs.
!
        integer, allocatable :: iDbio2(:)     ! 2D biological terms

        integer  :: iCOfx                     ! air-sea CO2 flux
        integer  :: iDNIT                     ! denitrification flux
        integer  :: ipCO2                     ! partial pressure of CO2
        integer  :: iO2fx                     ! air-sea O2 flux
!
!  Biological 3D diagnostic variable IDs.
!
        integer, allocatable :: iDbio3(:)     ! 3D biological terms

        integer  :: iPPro = 1                 ! primary productivity
        integer  :: iNO3u = 2                 ! NO3 uptake
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
!
!  Input/Output NetCDF files IDs.
!
        integer  :: ncADJid(Ngrids)     ! input/output adjoint
        integer  :: ncADSid(Ngrids)     ! input adjoint sensitivity
        integer  :: ncAVGid(Ngrids)     ! output averages
        integer  :: ncBRYid(Ngrids)     ! input boundary conditions
#ifdef INWAVE_COUPLING
        integer  :: ncIWBRYid(Ngrids)     ! input boundary conditions for Ac
        integer  :: ncIWINIid(Ngrids)     ! input initial conditions for Ac,cx,cy,ct
#endif
        integer  :: ncCLMid(Ngrids)     ! input climatology
        integer  :: ncDIAid(Ngrids)     ! output diagnostics
        integer  :: ncFLTid(Ngrids)     ! output floats
        integer  :: ncFRCid(NV,Ngrids)  ! input forcing
        integer  :: ncFWDid(Ngrids)     ! forward solution
        integer  :: ncGRDid(Ngrids)     ! input grid
        integer  :: ncGSTid(Ngrids)     ! input/output GST restart
        integer  :: ncHISid(Ngrids)     ! output history
        integer  :: ncINIid(Ngrids)     ! input/output NLM initial
        integer  :: ncIRPid(Ngrids)     ! input/output RPM initial
        integer  :: ncITLid(Ngrids)     ! input/output TLM initial
        integer  :: ncHSSid(Ngrids)     ! input/output Hessian vectors
        integer  :: ncMODid(Ngrids)     ! output 4DVAR fields
        integer  :: ncNRMid(Ngrids)     ! input/output covariance norm
        integer  :: ncOBSid(Ngrids)     ! input/output observations
        integer  :: ncRSTid(Ngrids)     ! input/output restart
        integer  :: ncSPSid(Ngrids)     ! Sparse matrix data
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
        integer, allocatable :: fltTid(:,:)   ! floats tracers IDs
        integer, allocatable :: hisTid(:,:)   ! history tracers IDs
        integer, allocatable :: iniTid(:,:)   ! initial NLM tracers IDs
        integer, allocatable :: irpTid(:,:)   ! initial RPM tracers IDs
        integer, allocatable :: itlTid(:,:)   ! initial TLM tracers IDs
        integer, allocatable :: hssTid(:,:)   ! Hessian vectors IDs
        integer, allocatable :: nrmTid(:,:)   ! norm tracers IDs
        integer, allocatable :: obsTid(:,:)   ! observations tracers IDs
        integer, allocatable :: rstTid(:,:)   ! restart tracers IDs
        integer, allocatable :: staTid(:,:)   ! stations tracers IDs
        integer, allocatable :: tlfTid(:,:)   ! TLM impulse tracers IDs
        integer, allocatable :: tlmTid(:,:)   ! tangent tracers IDs

        integer  :: adjVid(NV,Ngrids)    ! adjoint variables IDs
        integer  :: avgVid(NV,Ngrids)    ! averages variables IDs
        integer  :: diaVid(NV,Ngrids)    ! diagnostics variables IDs
        integer  :: fltVid(-6:NV,Ngrids) ! floats variables IDs
        integer  :: hisVid(NV,Ngrids)    ! history variables IDs
        integer  :: iniVid(NV,Ngrids)    ! initial NLM variables IDs
        integer  :: irpVid(NV,Ngrids)    ! initial RPM variables IDs
        integer  :: itlVid(NV,Ngrids)    ! initial TLM variables IDs
        integer  :: hssVid(NV,Ngrids)    ! Hessian variables IDs
        integer  :: modVid(NV,Ngrids)    ! 4DVAR variables IDs
        integer  :: nrmVid(NV,Ngrids)    ! norm variables IDs
        integer  :: obsVid(NV,Ngrids)    ! observations variables IDs
        integer  :: rstVid(NV,Ngrids)    ! restart variables IDs
        integer  :: staVid(NV,Ngrids)    ! stations variables IDs
        integer  :: tideVid(NV,Ngrids)   ! tide variables IDs
        integer  :: tlfVid(NV,Ngrids)    ! TLM impulse variables IDs
        integer  :: tlmVid(NV,Ngrids)    ! tangent variables IDs

        integer  :: tADJindx(Ngrids)     ! adjoint time record index
        integer  :: tAVGindx(Ngrids)     ! averages time record index
        integer  :: tDIAindx(Ngrids)     ! diagnostics time record index
        integer  :: tFLTindx(Ngrids)     ! floats time record index
        integer  :: tHISindx(Ngrids)     ! history time record index
        integer  :: tINIindx(Ngrids)     ! initial NLM time record index
        integer  :: tIRPindx(Ngrids)     ! initial RPM time record index
        integer  :: tITLindx(Ngrids)     ! initial TLM time record index
        integer  :: tHSSindx(Ngrids)     ! Hessian time record index
        integer  :: tNRMindx(Ngrids)     ! norm time record index
        integer  :: tRSTindx(Ngrids)     ! restart time record index
        integer  :: tSTAindx(Ngrids)     ! stations time record index
        integer  :: tTLFindx(Ngrids)     ! TLM impulse time record index
        integer  :: tTLMindx(Ngrids)     ! tangent time record index
!
!  Current NetCDF file information.
!
        integer  :: nvars          ! number of variables
        integer  :: nvdims(MV)     ! number of variable dimensions
        integer  :: tsize          ! size of unlimited dimension
        integer  :: vdims(5,MV)    ! variables dimensions IDs
        integer  :: vflag(MV)      ! variable flag: 1: full field
                                   !               -1: water points only
        character (len=80 ) :: type
!
!  Input/Output information variables.
!
        logical  :: Linfo(4,NV,Ngrids)

        integer  :: Iinfo(9,NV,Ngrids)

        real(r8) :: Finfo(9,NV,Ngrids)
        real(r8) :: Fpoint(2,NV,Ngrids)
        real(r8) :: Fscale(NV,Ngrids)
        real(r8) :: Tintrp(2,NV,Ngrids)

        real(r8), allocatable :: Vtime(:,:,:)

        character (len=5  ) :: version = '3.0  '
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
        character (len=80 ) :: Rdir
        character (len=80 ) :: Hdir
        character (len=80 ) :: Hfile
        character (len=80)  :: Adir
!
!  Analyical header file logical and names.
!
        logical :: Lanafile

        character (len=200), dimension(37) :: ANANAME
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
        integer :: Lvar, gtype, i, ic, ie, is, j, ng, varid

        real(r8) :: scale

        character (len=80), dimension(7) :: Vinfo
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
        allocate ( idTvar(MT) )
        allocate ( isTvar(MT) )
        allocate ( isTsur(MT) )

        allocate ( adjTid(MT,Ngrids) )
        allocate ( avgTid(MT,Ngrids) )
        allocate ( fltTid(MT,Ngrids) )
        allocate ( hisTid(MT,Ngrids) )
        allocate ( iniTid(MT,Ngrids) )
        allocate ( irpTid(MT,Ngrids) )
        allocate ( itlTid(MT,Ngrids) )
        allocate ( hssTid(MT,Ngrids) )
        allocate ( nrmTid(MT,Ngrids) )
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

#ifdef SEDIMENT
        allocate ( idBmas(NST) )
        allocate ( idfrac(NST) )
        allocate ( idUbld(NST) )
        allocate ( idVbld(NST) )
#endif
#ifdef DIAGNOSTICS
# if defined BIO_FASHAM && defined DIAGNOSTICS_BIO
        allocate ( iDbio2(NDbio2d) )
        allocate ( iDbio3(NDbio3d) )
# endif
        allocate ( idDtrc(MT,NDT) )
        allocate ( idDu2d(NDM2d) )
        allocate ( idDv2d(NDM2d) )
        allocate ( idDu3d(NDM3d) )
        allocate ( idDv3d(NDM3d) )
#endif
        allocate ( idSvar(MAXVAL(NSV)+1) )

#if defined SSH_TIDES || defined UV_TIDES
        allocate ( Vtime(MAX(2,MTC),NV,Ngrids) )
#else
        allocate ( Vtime(2,NV,Ngrids) )
#endif
!
!-----------------------------------------------------------------------
!  Set lower and upper observation processing (I,J) ranges: Either the
!  full grid or only interior points will be considered.
!-----------------------------------------------------------------------
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
          ncBRYid(ng)=-1
          ncCLMid(ng)=-1
          ncDIAid(ng)=-1
          ncFLTid(ng)=-1
          ncGSTid(ng)=-1
          ncFWDid(ng)=-1
          ncGRDid(ng)=-1
          ncHISid(ng)=-1
          ncHSSid(ng)=-1
          ncINIid(ng)=-1
          ncIRPid(ng)=-1
          ncITLid(ng)=-1
          ncMODid(ng)=-1
          ncNRMid(ng)=-1
          ncOBSid(ng)=-1
          ncRSTid(ng)=-1
          ncSPSid(ng)=-1
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
          tNRMindx(ng)=0
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
!  Analytical files switch and names
!
        Lanafile=.TRUE.
        DO i=1,37
          DO j=1,200
            ANANAME(i)(j:j)=' '
          END DO
        END DO
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
#if defined BIO_FASHAM && defined DIAGNOSTICS_BIO
!
!  Initialize biology diagnostic indices.
!
        i=0
# ifdef DENITRIFICATION
        iDNIT=i+1
        i=i+1
# endif
# ifdef CARBON
        iCOfx=i+1
        ipCO2=i+2
        i=i+2
# endif
# ifdef OXYGEN
        iO2fx=i+1
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
              CASE ('idevap')
                idevap=varid
              CASE ('idUair')
                idUair=varid
              CASE ('idVair')
                idVair=varid
              CASE ('idWamp')
                idWamp=varid
              CASE ('idWbrk')
                idWbrk=varid
              CASE ('idWdis')
                idWdis=varid
              CASE ('idWdir')
                idWdir=varid
              CASE ('idWlen')
                idWlen=varid
              CASE ('idWmsk')
                idWmsk=varid
              CASE ('idWptp')
                idWptp=varid
              CASE ('idWpbt')
                idWpbt=varid
              CASE ('idWorb')
                idWorb=varid
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
#if defined NPZD_FRANKS || defined NPZD_POWELL
              CASE ('idRtrc(iNO3_)')
                idRtrc(iNO3_)=varid
              CASE ('idRtrc(iPhyt)')
                idRtrc(iPhyt)=varid
              CASE ('idRtrc(iZoop)')
                idRtrc(iZoop)=varid
              CASE ('idRtrc(iSDet)')
                idRtrc(iSDet)=varid
#endif
#ifdef BIO_FASHAM
              CASE ('idRtrc(iNO3_)')
                idRtrc(iNO3_)=varid
              CASE ('idRtrc(iNH4_)')
                idRtrc(iNH4_)=varid
              CASE ('idRtrc(iPhyt)')
                idRtrc(iPhyt)=varid
              CASE ('idRtrc(iZoop)')
                idRtrc(iZoop)=varid
              CASE ('idRtrc(iLDeN)')
                idRtrc(iLDeN)=varid
              CASE ('idRtrc(iSDeN)')
                idRtrc(iSDeN)=varid
              CASE ('idRtrc(iChlo)')
                idRtrc(iChlo)=varid
# ifdef CARBON
              CASE ('idRtrc(iTIC_)')
                idRtrc(iTIC_)=varid
              CASE ('idRtrc(iTAlk)')
                idRtrc(iTAlk)=varid
              CASE ('idRtrc(iLDeC)')
                idRtrc(iLDeC)=varid
              CASE ('idRtrc(iSDeC)')
                idRtrc(iSDeC)=varid
# endif
# ifdef OXYGEN
              CASE ('idRtrc(iOxyg)')
                idRtrc(iOxyg)=varid
# endif
#endif
#ifdef NEMURO
              CASE ('idRtrc(iNO3_)')
                idRtrc(iNO3_)=varid
              CASE ('idRtrc(iNH4_)')
                idRtrc(iNH4_)=varid
              CASE ('idRtrc(iDON_)')
                idRtrc(iDON_)=varid
              CASE ('idRtrc(iPON_)')
                idRtrc(iPON_)=varid
#endif
#ifdef SEDIMENT
              CASE ('idRtrc(idmud)')
                load=.TRUE.
              CASE ('idRtrc(isand)')
                load=.TRUE.
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
              CASE ('idSSTc')
                idSSTc=varid
              CASE ('idTclm(itemp)')
                idTclm(itemp)=varid
              CASE ('idTclm(isalt)')
                idTclm(isalt)=varid
#if defined AD_SENSITIVITY   || defined OBS_SENSITIVITY || \
    defined OPT_OBSERVATIONS || defined SO_SEMI
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
#endif
#if defined SEDIMENT || defined BBL_MODEL
              CASE ('idTvar(idmud(i))')
                load=.TRUE.
              CASE ('idTvar(isand(i))')
                load=.TRUE.                
              CASE ('idfrac')
                load=.TRUE.
              CASE ('idBmas')
                load=.TRUE.
# ifdef BEDLOAD
              CASE ('idUbld')
                load=.TRUE.
              CASE ('idVbld')
                load=.TRUE.
# endif
              CASE ('idSbed(ithck)')
                idSbed(ithck)=varid
              CASE ('idSbed(iaged)')
                idSbed(iaged)=varid
              CASE ('idSbed(iporo)')
                idSbed(iporo)=varid
              CASE ('idSbed(idiff)')
                idSbed(idiff)=varid
# if defined COHESIVE_BED || defined SED_BIODIFF || defined MIXED_BED
              CASE ('idSbed(ibtcr)')
                idSbed(ibtcr)=varid
# endif
              CASE ('idBott(isd50)')
                idBott(isd50)=varid
              CASE ('idBott(idens)')
                idBott(idens)=varid
              CASE ('idBott(iwsed)')
                idBott(iwsed)=varid
              CASE ('idBott(itauc)')
                idBott(itauc)=varid
              CASE ('idBott(irlen)')
                idBott(irlen)=varid
              CASE ('idBott(irhgt)')
                idBott(irhgt)=varid
              CASE ('idBott(ibwav)')
                idBott(ibwav)=varid
              CASE ('idBott(izdef)')
                idBott(izdef)=varid
              CASE ('idBott(izapp)')
                idBott(izapp)=varid
              CASE ('idBott(izNik)')
                idBott(izNik)=varid
              CASE ('idBott(izbio)')
                idBott(izbio)=varid
              CASE ('idBott(izbfm)')
                idBott(izbfm)=varid
              CASE ('idBott(izbld)')
                idBott(izbld)=varid
              CASE ('idBott(izwbl)')
                idBott(izwbl)=varid
              CASE ('idBott(iactv)')
                idBott(iactv)=varid
              CASE ('idBott(ishgt)')
                idBott(ishgt)=varid
              CASE ('idBott(idflux)')
                idBott(idflux)=varid
              CASE ('idBott(idnet)')
                idBott(idnet)=varid
# if defined COHESIVE_BED || defined SED_BIODIFF || defined MIXED_BED
              CASE ('idBott(idoff)')
                idBott(idoff)=varid
              CASE ('idBott(idslp)')
                idBott(idslp)=varid
              CASE ('idBott(idtim)')
                idBott(idtim)=varid
              CASE ('idBott(idbmx)')
                idBott(idbmx)=varid
              CASE ('idBott(idbmm)')
                idBott(idbmm)=varid
              CASE ('idBott(idbzs)')
                idBott(idbzs)=varid
              CASE ('idBott(idbzm)')
                idBott(idbzm)=varid
              CASE ('idBott(idbzp)')
                idBott(idbzp)=varid
# endif
# if defined MIXED_BED
              CASE ('idBott(idprp)')
                idBott(idprp)=varid
# endif
#endif
#if defined NPZD_FRANKS || defined NPZD_POWELL
              CASE ('idTvar(iNO3_)')
                idTvar(iNO3_)=varid
              CASE ('idTvar(iPhyt)')
                idTvar(iPhyt)=varid
              CASE ('idTvar(iZoop)')
                idTvar(iZoop)=varid
              CASE ('idTvar(iSDet)')
                idTvar(iSDet)=varid
# if defined AD_SENSITIVITY   || defined OBS_SENSITIVITY || \
     defined OPT_OBSERVATIONS || defined SO_SEMI
              CASE ('idTads(iNO3_)')
                idTads(iNO3_)=varid
              CASE ('idTads(iPhyt)')
                idTads(iPhyt)=varid
              CASE ('idTads(iZoop)')
                idTads(iZoop)=varid
              CASE ('idTads(iSDet)')
                idTads(iSDet)=varid
# endif
              CASE ('idTbry(iwest,iNO3_)')
                idTbry(iwest,iNO3_)=varid
              CASE ('idTbry(ieast,iNO3_)')
                idTbry(ieast,iNO3_)=varid
              CASE ('idTbry(isouth,iNO3_)')
                idTbry(isouth,iNO3_)=varid
              CASE ('idTbry(inorth,iNO3_)')
                idTbry(inorth,iNO3_)=varid
              CASE ('idTbry(iwest,iPhyt)')
                idTbry(iwest,iPhyt)=varid
              CASE ('idTbry(ieast,iPhyt)')
                idTbry(ieast,iPhyt)=varid
              CASE ('idTbry(isouth,iPhyt)')
                idTbry(isouth,iPhyt)=varid
              CASE ('idTbry(inorth,iPhyt)')
                idTbry(inorth,iPhyt)=varid
              CASE ('idTbry(iwest,iZoop)')
                idTbry(iwest,iZoop)=varid
              CASE ('idTbry(ieast,iZoop)')
                idTbry(ieast,iZoop)=varid
              CASE ('idTbry(isouth,iZoop)')
                idTbry(isouth,iZoop)=varid
              CASE ('idTbry(inorth,iZoop)')
                idTbry(inorth,iZoop)=varid
              CASE ('idTbry(iwest,iSDet)')
                idTbry(iwest,iSDet)=varid
              CASE ('idTbry(ieast,iSDet)')
                idTbry(ieast,iSDet)=varid
              CASE ('idTbry(isouth,iSDet)')
                idTbry(isouth,iSDet)=varid
              CASE ('idTbry(inorth,iSDet)')
                idTbry(inorth,iSDet)=varid
#endif
#ifdef BIO_FASHAM
              CASE ('idTvar(iNO3_)')
                idTvar(iNO3_)=varid
              CASE ('idTvar(iNH4_)')
                idTvar(iNH4_)=varid
              CASE ('idTvar(iPhyt)')
                idTvar(iPhyt)=varid
              CASE ('idTvar(iZoop)')
                idTvar(iZoop)=varid
              CASE ('idTvar(iLDeN)')
                idTvar(iLDeN)=varid
              CASE ('idTvar(iSDeN)')
                idTvar(iSDeN)=varid
              CASE ('idTvar(iChlo)')
                idTvar(iChlo)=varid
# if defined AD_SENSITIVITY   || defined OBS_SENSITIVITY || \
     defined OPT_OBSERVATIONS || defined SO_SEMI
              CASE ('idTads(iNO3_)')
                idTads(iNO3_)=varid
              CASE ('idTads(iNH4_)')
                idTads(iNH4_)=varid
              CASE ('idTads(iPhyt)')
                idTads(iPhyt)=varid
              CASE ('idTads(iZoop)')
                idTads(iZoop)=varid
              CASE ('idTads(iLDeN)')
                idTads(iLDeN)=varid
              CASE ('idTads(iSDeN)')
                idTads(iSDeN)=varid
              CASE ('idTads(iChlo)')
                idTads(iChlo)=varid
# endif
# ifdef CARBON
              CASE ('idTvar(iTIC_)')
                idTvar(iTIC_)=varid
              CASE ('idTvar(iTAlk)')
                idTvar(iTAlk)=varid
              CASE ('idTvar(iLDeC)')
                idTvar(iLDeC)=varid
              CASE ('idTvar(iSDeC)')
                idTvar(iSDeC)=varid
#  if defined AD_SENSITIVITY   || defined OBS_SENSITIVITY || \
      defined OPT_OBSERVATIONS || defined SO_SEMI
              CASE ('idTads(iTIC_)')
                idTads(iTIC_)=varid
              CASE ('idTads(iTAlk)')
                idTads(iTAlk)=varid
              CASE ('idTads(iLDeC)')
                idTads(iLDeC)=varid
              CASE ('idTads(iSDeC)')
                idTads(iSDeC)=varid
#  endif
# endif
# ifdef OXYGEN
              CASE ('idTvar(iOxyg)')
                idTvar(iOxyg)=varid
#  if defined AD_SENSITIVITY   || defined OBS_SENSITIVITY || \
      defined OPT_OBSERVATIONS || defined SO_SEMI
              CASE ('idTads(iOxyg)')
                idTads(iOxyg)=varid
#  endif
# endif
              CASE ('idTbry(iwest,iNO3_)')
                idTbry(iwest,iNO3_)=varid
              CASE ('idTbry(ieast,iNO3_)')
                idTbry(ieast,iNO3_)=varid
              CASE ('idTbry(isouth,iNO3_)')
                idTbry(isouth,iNO3_)=varid
              CASE ('idTbry(inorth,iNO3_)')
                idTbry(inorth,iNO3_)=varid
              CASE ('idTbry(iwest,iNH4_)')
                idTbry(iwest,iNH4_)=varid
              CASE ('idTbry(ieast,iNH4_)')
                idTbry(ieast,iNH4_)=varid
              CASE ('idTbry(isouth,iNH4_)')
                idTbry(isouth,iNH4_)=varid
              CASE ('idTbry(inorth,iNH4_)')
                idTbry(inorth,iNH4_)=varid
              CASE ('idTbry(iwest,iPhyt)')
                idTbry(iwest,iPhyt)=varid
              CASE ('idTbry(ieast,iPhyt)')
                idTbry(ieast,iPhyt)=varid
              CASE ('idTbry(isouth,iPhyt)')
                idTbry(isouth,iPhyt)=varid
              CASE ('idTbry(inorth,iPhyt)')
                idTbry(inorth,iPhyt)=varid
              CASE ('idTbry(iwest,iZoop)')
                idTbry(iwest,iZoop)=varid
              CASE ('idTbry(ieast,iZoop)')
                idTbry(ieast,iZoop)=varid
              CASE ('idTbry(isouth,iZoop)')
                idTbry(isouth,iZoop)=varid
              CASE ('idTbry(inorth,iZoop)')
                idTbry(inorth,iZoop)=varid
              CASE ('idTbry(iwest,iSDeN)')
                idTbry(iwest,iSDeN)=varid
              CASE ('idTbry(ieast,iSDeN)')
                idTbry(ieast,iSDeN)=varid
              CASE ('idTbry(isouth,iSDeN)')
                idTbry(isouth,iSDeN)=varid
              CASE ('idTbry(inorth,iSDeN)')
                idTbry(inorth,iSDeN)=varid
              CASE ('idTbry(iwest,iLDeN)')
                idTbry(iwest,iLDeN)=varid
              CASE ('idTbry(ieast,iLDeN)')
                idTbry(ieast,iLDeN)=varid
              CASE ('idTbry(isouth,iLDeN)')
                idTbry(isouth,iLDeN)=varid
              CASE ('idTbry(inorth,iLDeN)')
                idTbry(inorth,iLDeN)=varid
              CASE ('idTbry(iwest,iChlo)')
                idTbry(iwest,iChlo)=varid
              CASE ('idTbry(ieast,iChlo)')
                idTbry(ieast,iChlo)=varid
              CASE ('idTbry(isouth,iChlo)')
                idTbry(isouth,iChlo)=varid
              CASE ('idTbry(inorth,iChlo)')
                idTbry(inorth,iChlo)=varid
# ifdef CARBON
              CASE ('idTbry(iwest,iSDeC)')
                idTbry(iwest,iSDeC)=varid
              CASE ('idTbry(ieast,iSDeC)')
                idTbry(ieast,iSDeC)=varid
              CASE ('idTbry(isouth,iSDeC)')
                idTbry(isouth,iSDeC)=varid
              CASE ('idTbry(inorth,iSDeC)')
                idTbry(inorth,iSDeC)=varid
              CASE ('idTbry(iwest,iLDeC)')
                idTbry(iwest,iLDeC)=varid
              CASE ('idTbry(ieast,iLDeC)')
                idTbry(ieast,iLDeC)=varid
              CASE ('idTbry(isouth,iLDeC)')
                idTbry(isouth,iLDeC)=varid
              CASE ('idTbry(inorth,iLDeC)')
                idTbry(inorth,iLDeC)=varid
              CASE ('idTbry(iwest,iTIC_)')
                idTbry(iwest,iTIC_)=varid
              CASE ('idTbry(ieast,iTIC_)')
                idTbry(ieast,iTIC_)=varid
              CASE ('idTbry(isouth,iTIC_)')
                idTbry(isouth,iTIC_)=varid
              CASE ('idTbry(inorth,iTIC_)')
                idTbry(inorth,iTIC_)=varid
              CASE ('idTbry(iwest,iTAlk)')
                idTbry(iwest,iTAlk)=varid
              CASE ('idTbry(ieast,iTAlk)')
                idTbry(ieast,iTAlk)=varid
              CASE ('idTbry(isouth,iTAlk)')
                idTbry(isouth,iTAlk)=varid
              CASE ('idTbry(inorth,iTAlk)')
                idTbry(inorth,iTAlk)=varid
# endif
# ifdef OXYGEN
              CASE ('idTbry(iwest,iOxyg)')
                idTbry(iwest,iOxyg)=varid
              CASE ('idTbry(ieast,iOxyg)')
                idTbry(ieast,iOxyg)=varid
              CASE ('idTbry(isouth,iOxyg)')
                idTbry(isouth,iOxyg)=varid
              CASE ('idTbry(inorth,iOxyg)')
                idTbry(inorth,iOxyg)=varid
# endif
# ifdef DIAGNOSTICS_BIO
              CASE ('iDbio2(iDNIT)')
                iDbio2(iDNIT)=varid
#  ifdef CARBON
              CASE ('iDbio2(iCOfx)')
                iDbio2(iCOfx)=varid
              CASE ('iDbio2(ipCO2)')
                iDbio2(ipCO2)=varid
#  endif
#  ifdef OXYGEN
              CASE ('iDbio2(iO2fx)')
                iDbio2(iO2fx)=varid
#  endif
              CASE ('iDbio3(iPPro)')
                iDbio3(iPPro)=varid
              CASE ('iDbio3(iNO3u)')
                iDbio3(iNO3u)=varid
# endif
#endif
#ifdef BIO_LIMADONEY
              CASE ('idTvar(iNO3_)')
                idTvar(iNO3_)=varid
              CASE ('idTvar(iNH4_)')
                idTvar(iNH4_)=varid
              CASE ('idTvar(iPhy1)')
                idTvar(iPhy1)=varid
              CASE ('idTvar(iPhy2)')
                idTvar(iPhy2)=varid
              CASE ('idTvar(iPhC1)')
                idTvar(iPhC1)=varid
              CASE ('idTvar(iPhC2)')
                idTvar(iPhC2)=varid
              CASE ('idTvar(iChl1)')
                idTvar(iChl1)=varid
              CASE ('idTvar(iChl2)')
                idTvar(iChl2)=varid
              CASE ('idTvar(iZoo_)')
                idTvar(iZoo_)=varid
              CASE ('idTvar(iZooC)')
                idTvar(iZooC)=varid
              CASE ('idTvar(iDet1)')
                idTvar(iDet1)=varid
              CASE ('idTvar(iDet2)')
                idTvar(iDet2)=varid
              CASE ('idTvar(iDeC1)')
                idTvar(iDeC1)=varid
              CASE ('idTvar(iDeC2)')
                idTvar(iDeC2)=varid
              CASE ('idTbry(iwest,iNO3_)')
                idTbry(iwest,iNO3_)=varid
              CASE ('idTbry(ieast,iNO3_)')
                idTbry(ieast,iNO3_)=varid
              CASE ('idTbry(isouth,iNO3_)')
                idTbry(isouth,iNO3_)=varid
              CASE ('idTbry(inorth,iNO3_)')
                idTbry(inorth,iNO3_)=varid
              CASE ('idTbry(iwest,iNH4_)')
                idTbry(iwest,iNH4_)=varid
              CASE ('idTbry(ieast,iNH4_)')
                idTbry(ieast,iNH4_)=varid
              CASE ('idTbry(isouth,iNH4_)')
                idTbry(isouth,iNH4_)=varid
              CASE ('idTbry(inorth,iNH4_)')
                idTbry(inorth,iNH4_)=varid
              CASE ('idTbry(iwest,iPhy1)')
                idTbry(iwest,iPhy1)=varid
              CASE ('idTbry(ieast,iPhy1)')
                idTbry(ieast,iPhy1)=varid
              CASE ('idTbry(isouth,iPhy1)')
                idTbry(isouth,iPhy1)=varid
              CASE ('idTbry(inorth,iPhy1)')
                idTbry(inorth,iPhy1)=varid
              CASE ('idTbry(iwest,iPhy2)')
                idTbry(iwest,iPhy2)=varid
              CASE ('idTbry(ieast,iPhy2)')
                idTbry(ieast,iPhy2)=varid
              CASE ('idTbry(isouth,iPhy2)')
                idTbry(isouth,iPhy2)=varid
              CASE ('idTbry(inorth,iPhy2)')
                idTbry(inorth,iPhy2)=varid
              CASE ('idTbry(iwest,iPhC1)')
                idTbry(iwest,iPhC1)=varid
              CASE ('idTbry(ieast,iPhC1)')
                idTbry(ieast,iPhC1)=varid
              CASE ('idTbry(isouth,iPhC1)')
                idTbry(isouth,iPhC1)=varid
              CASE ('idTbry(inorth,iPhC1)')
                idTbry(inorth,iPhC1)=varid
              CASE ('idTbry(iwest,iPhC2)')
                idTbry(iwest,iPhC2)=varid
              CASE ('idTbry(ieast,iPhC2)')
                idTbry(ieast,iPhC2)=varid
              CASE ('idTbry(isouth,iPhC2)')
                idTbry(isouth,iPhC2)=varid
              CASE ('idTbry(inorth,iPhC2)')
                idTbry(inorth,iPhC2)=varid
              CASE ('idTbry(iwest,iZoo_)')
                idTbry(iwest,iZoo_)=varid
              CASE ('idTbry(ieast,iZoo_)')
                idTbry(ieast,iZoo_)=varid
              CASE ('idTbry(isouth,iZoo_)')
                idTbry(isouth,iZoo_)=varid
              CASE ('idTbry(inorth,iZoo_)')
                idTbry(inorth,iZoo_)=varid
              CASE ('idTbry(iwest,iZooC)')
                idTbry(iwest,iZooC)=varid
              CASE ('idTbry(ieast,iZooC)')
                idTbry(ieast,iZooC)=varid
              CASE ('idTbry(isouth,iZooC)')
                idTbry(isouth,iZooC)=varid
              CASE ('idTbry(inorth,iZooC)')
                idTbry(inorth,iZooC)=varid
              CASE ('idTbry(iwest,iDet1)')
                idTbry(iwest,iDet1)=varid
              CASE ('idTbry(ieast,iDet1)')
                idTbry(ieast,iDet1)=varid
              CASE ('idTbry(isouth,iDet1)')
                idTbry(isouth,iDet1)=varid
              CASE ('idTbry(inorth,iDet1)')
                idTbry(inorth,iDet1)=varid
              CASE ('idTbry(iwest,iDet2)')
                idTbry(iwest,iDet2)=varid
              CASE ('idTbry(ieast,iDet2)')
                idTbry(ieast,iDet2)=varid
              CASE ('idTbry(isouth,iDet2)')
                idTbry(isouth,iDet2)=varid
              CASE ('idTbry(inorth,iDet2)')
                idTbry(inorth,iDet2)=varid
              CASE ('idTbry(iwest,iDeC1)')
                idTbry(iwest,iDeC1)=varid
              CASE ('idTbry(ieast,iDeC1)')
                idTbry(ieast,iDeC1)=varid
              CASE ('idTbry(isouth,iDeC1)')
                idTbry(isouth,iDeC1)=varid
              CASE ('idTbry(inorth,iDeC1)')
                idTbry(inorth,iDeC1)=varid
              CASE ('idTbry(iwest,iDeC2)')
                idTbry(iwest,iDeC2)=varid
              CASE ('idTbry(ieast,iDeC2)')
                idTbry(ieast,iDeC2)=varid
              CASE ('idTbry(isouth,iDeC2)')
                idTbry(isouth,iDeC2)=varid
              CASE ('idTbry(inorth,iDeC2)')
                idTbry(inorth,iDeC2)=varid
              CASE ('idTbry(iwest,iChl1)')
                idTbry(iwest,iChl1)=varid
              CASE ('idTbry(ieast,iChl1)')
                idTbry(ieast,iChl1)=varid
              CASE ('idTbry(isouth,iChl1)')
                idTbry(isouth,iChl1)=varid
              CASE ('idTbry(inorth,iChl1)')
                idTbry(inorth,iChl1)=varid
              CASE ('idTbry(iwest,iChl2)')
                idTbry(iwest,iChl2)=varid
              CASE ('idTbry(ieast,iChl2)')
                idTbry(ieast,iChl2)=varid
              CASE ('idTbry(isouth,iChl2)')
                idTbry(isouth,iChl2)=varid
              CASE ('idTbry(inorth,iChl2)')
                idTbry(inorth,iChl2)=varid
#endif
#ifdef NEMURO
              CASE ('idTvar(iLphy)')
                idTvar(iLphy)=varid
              CASE ('idTvar(iSphy)')
                idTvar(iSphy)=varid
              CASE ('idTvar(iLzoo)')
                idTvar(iLzoo)=varid
              CASE ('idTvar(iSzoo)')
                idTvar(iSzoo)=varid
              CASE ('idTvar(iPzoo)')
                idTvar(iPzoo)=varid
              CASE ('idTvar(iNO3_)')
                idTvar(iNO3_)=varid
              CASE ('idTvar(iNH4_)')
                idTvar(iNH4_)=varid
              CASE ('idTvar(iPON_)')
                idTvar(iPON_)=varid
              CASE ('idTvar(iDON_)')
                idTvar(iDON_)=varid
              CASE ('idTvar(iSiOH)')
                idTvar(iSiOH)=varid
              CASE ('idTvar(iopal)')
                idTvar(iopal)=varid
              CASE ('idTbry(iwest,iLphy)')
                idTbry(iwest,iLphy)=varid
              CASE ('idTbry(ieast,iLphy)')
                idTbry(ieast,iLphy)=varid
              CASE ('idTbry(isouth,iLphy)')
                idTbry(isouth,iLphy)=varid
              CASE ('idTbry(inorth,iLphy)')
                idTbry(inorth,iLphy)=varid
              CASE ('idTbry(iwest,iSphy)')
                idTbry(iwest,iSphy)=varid
              CASE ('idTbry(ieast,iSphy)')
                idTbry(ieast,iSphy)=varid
              CASE ('idTbry(isouth,iSphy)')
                idTbry(isouth,iSphy)=varid
              CASE ('idTbry(inorth,iSphy)')
                idTbry(inorth,iSphy)=varid
              CASE ('idTbry(iwest,iLzoo)')
                idTbry(iwest,iLzoo)=varid
              CASE ('idTbry(ieast,iLzoo)')
                idTbry(ieast,iLzoo)=varid
              CASE ('idTbry(isouth,iLzoo)')
                idTbry(isouth,iLzoo)=varid
              CASE ('idTbry(inorth,iLzoo)')
                idTbry(inorth,iLzoo)=varid
              CASE ('idTbry(iwest,iSzoo)')
                idTbry(iwest,iSzoo)=varid
              CASE ('idTbry(ieast,iSzoo)')
                idTbry(ieast,iSzoo)=varid
              CASE ('idTbry(isouth,iSzoo)')
                idTbry(isouth,iSzoo)=varid
              CASE ('idTbry(inorth,iSzoo)')
                idTbry(inorth,iSzoo)=varid
              CASE ('idTbry(iwest,iPzoo)')
                idTbry(iwest,iPzoo)=varid
              CASE ('idTbry(ieast,iPzoo)')
                idTbry(ieast,iPzoo)=varid
              CASE ('idTbry(isouth,iPzoo)')
                idTbry(isouth,iPzoo)=varid
              CASE ('idTbry(inorth,iPzoo)')
                idTbry(inorth,iPzoo)=varid
              CASE ('idTbry(iwest,iNO3_)')
                idTbry(iwest,iNO3_)=varid
              CASE ('idTbry(ieast,iNO3_)')
                idTbry(ieast,iNO3_)=varid
              CASE ('idTbry(isouth,iNO3_)')
                idTbry(isouth,iNO3_)=varid
              CASE ('idTbry(inorth,iNO3_)')
                idTbry(inorth,iNO3_)=varid
              CASE ('idTbry(iwest,iNH4_)')
                idTbry(iwest,iNH4_)=varid
              CASE ('idTbry(ieast,iNH4_)')
                idTbry(ieast,iNH4_)=varid
              CASE ('idTbry(isouth,iNH4_)')
                idTbry(isouth,iNH4_)=varid
              CASE ('idTbry(inorth,iNH4_)')
                idTbry(inorth,iNH4_)=varid
              CASE ('idTbry(iwest,iPON_)')
                idTbry(iwest,iPON_)=varid
              CASE ('idTbry(ieast,iPON_)')
                idTbry(ieast,iPON_)=varid
              CASE ('idTbry(isouth,iPON_)')
                idTbry(isouth,iPON_)=varid
              CASE ('idTbry(inorth,iPON_)')
                idTbry(inorth,iPON_)=varid
              CASE ('idTbry(iwest,iDON_)')
                idTbry(iwest,iDON_)=varid
              CASE ('idTbry(ieast,iDON_)')
                idTbry(ieast,iDON_)=varid
              CASE ('idTbry(isouth,iDON_)')
                idTbry(isouth,iDON_)=varid
              CASE ('idTbry(inorth,iDON_)')
                idTbry(inorth,iDON_)=varid
              CASE ('idTbry(iwest,iSiOH)')
                idTbry(iwest,iSiOH)=varid
              CASE ('idTbry(ieast,iSiOH)')
                idTbry(ieast,iSiOH)=varid
              CASE ('idTbry(isouth,iSiOH)')
                idTbry(isouth,iSiOH)=varid
              CASE ('idTbry(inorth,iSiOH)')
                idTbry(inorth,iSiOH)=varid
              CASE ('idTbry(iwest,iopal)')
                idTbry(iwest,iopal)=varid
              CASE ('idTbry(ieast,iopal)')
                idTbry(ieast,iopal)=varid
              CASE ('idTbry(isouth,iopal)')
                idTbry(isouth,iopal)=varid
              CASE ('idTbry(inorth,iopal)')
                idTbry(inorth,iopal)=varid
#endif
#ifdef ECOSIM
              CASE ('idTvar(iNO3_)')
                idTvar(iNO3_)=varid
              CASE ('idTvar(iNH4_)')
                idTvar(iNH4_)=varid
              CASE ('idTvar(iPO4_)')
                idTvar(iPO4_)=varid
              CASE ('idTvar(iFeO_)')
                idTvar(iFeO_)=varid
              CASE ('idTvar(iSiO_)')
                idTvar(iSiO_)=varid
              CASE ('idTvar(iDIC_)')
                idTvar(iDIC_)=varid
              CASE ('idTvar(iBacC)')
                load=.TRUE.
              CASE ('idTvar(iBacF)')
                load=.TRUE.
              CASE ('idTvar(iBacN)')
                load=.TRUE.
              CASE ('idTvar(iBacP)')
                load=.TRUE.
              CASE ('idTvar(iCDMC)')
                load=.TRUE.
              CASE ('idTvar(iDOMC)')
                load=.TRUE.
              CASE ('idTvar(iDOMN)')
                load=.TRUE.
              CASE ('idTvar(iDOMP)')
                load=.TRUE.
              CASE ('idTvar(iFecC)')
                load=.TRUE.
              CASE ('idTvar(iFecF)')
                load=.TRUE.
              CASE ('idTvar(iFecN)')
                load=.TRUE.
              CASE ('idTvar(iFecP)')
                load=.TRUE.
              CASE ('idTvar(iFecS)')
                load=.TRUE.
              CASE ('idTvar(iPhyC)')
                load=.TRUE.
              CASE ('idTvar(iPhyF)')
                load=.TRUE.
              CASE ('idTvar(iPhyN)')
                load=.TRUE.
              CASE ('idTvar(iPhyP)')
                load=.TRUE.
              CASE ('idTvar(iPhyS)')
                load=.TRUE.
              CASE ('idTvar(iPigs)')
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
              CASE ('idDu2d(M2hadv)')
                idDu2d(M2hadv)=varid
              CASE ('idDv2d(M2hadv)')
                idDv2d(M2hadv)=varid
# endif
# ifdef NEARSHORE_MELLOR
              CASE ('idDu2d(M2hrad)')
                idDu2d(M2hrad)=varid
              CASE ('idDv2d(M2hrad)')
                idDv2d(M2hrad)=varid
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
              CASE ('idDv2d(M2hvis)')
                idDv2d(M2hvis)=varid
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
              CASE ('idDu3d(M3hadv)')
                idDu3d(M3hadv)=varid
              CASE ('idDv3d(M3hadv)')
                idDv3d(M3hadv)=varid
              CASE ('idDu3d(M3vadv)')
                idDu3d(M3vadv)=varid
              CASE ('idDv3d(M3vadv)')
                idDv3d(M3vadv)=varid
#  endif
#  ifdef NEARSHORE_MELLOR
              CASE ('idDu3d(M3hrad)')
                idDu3d(M3hrad)=varid
              CASE ('idDv3d(M3hrad)')
                idDv3d(M3hrad)=varid
              CASE ('idDu3d(M3vrad)')
                idDu3d(M3vrad)=varid
              CASE ('idDv3d(M3vrad)')
                idDv3d(M3vrad)=varid
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
              CASE ('idDv3d(M3hvis)')
                idDv3d(M3hvis)=varid
#  endif
# endif
#endif
#ifdef DIAGNOSTICS_TS
              CASE ('idDtrc(iTrate)')
                load=.TRUE.
              CASE ('idDtrc(iThadv)')
                load=.TRUE.
              CASE ('idDtrc(iTvadv)')
                load=.TRUE.
# if defined TS_DIF2 || defined TS_DIF4
              CASE ('idDtrc(iThdif)')
                load=.TRUE.
# endif
              CASE ('idDtrc(iTvdif)')
                load=.TRUE.
#endif
#if defined FORWARD_READ || defined FORWARD_WRITE || \
    defined COMPOSED_GRID
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
# if defined WEAK_CONSTRAINT || defined IOM
              CASE ('idCG_P')
                idCG_P=varid
              CASE ('idCG_R')
                idCG_R=varid
              CASE ('idCG_S')
                idCG_S=varid
              CASE ('idCG_V')
                idCG_V=varid
              CASE ('idCG_X')
                idCG_X=varid
              CASE ('idCG_Z')
                idCG_Z=varid
# endif
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
# if defined WEAK_CONSTRAINT || defined IOM
              CASE ('idRepC')
                idRepC=varid
# endif
              CASE ('idTLmo')
                idTLmo=varid
              CASE ('idMOMi')
                idMOMi=varid
              CASE ('idMOMf')
                idMOMf=varid
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
              END SELECT
#endif
#ifdef SEDIMENT
!
!  Adjust information for all sediment tracers.
!
              SELECT CASE (Vinfo(1))
                CASE ('mud_')
                  IF (NCS.gt.0) THEN
                    varid=varid-1
                    DO i=1,NCS
                      varid=varid+1
                      idTvar(idmud(i))=varid
                      DO ng=1,Ngrids
                        Fscale(varid,ng)=scale
                        Iinfo(1,varid,ng)=gtype
                      END DO
                      WRITE (Vname(1,varid),'(a,i2.2)')                 &
     &                      TRIM(ADJUSTL(Vinfo(1))), i
                      WRITE (Vname(2,varid),'(a,a,i2.2)')               &
     &                      TRIM(ADJUSTL(Vinfo(2))), ', size class ', i
                      WRITE (Vname(3,varid),'(a)')                      &
     &                      TRIM(ADJUSTL(Vinfo(3)))
                      WRITE (Vname(4,varid),'(a,a)')                    &
     &                      TRIM(Vname(1,varid)), ', scalar, series'
                      WRITE (Vname(5,varid),'(a)')                      &
     &                      TRIM(ADJUSTL(Vinfo(5)))
                    END DO
                  END IF
                CASE ('sand_')
                  IF (NNS.gt.0) THEN
                    varid=varid-1
                    DO i=1,NNS
                      varid=varid+1
                      idTvar(isand(i))=varid
                      DO ng=1,Ngrids
                        Fscale(varid,ng)=scale
                        Iinfo(1,varid,ng)=gtype
                      END DO
                      WRITE (Vname(1,varid),'(a,i2.2)')                 &
     &                      TRIM(ADJUSTL(Vinfo(1))), i
                      WRITE (Vname(2,varid),'(a,a,i2.2)')               &
     &                      TRIM(ADJUSTL(Vinfo(2))), ', size class ', i
                      WRITE (Vname(3,varid),'(a)')                      &
     &                      TRIM(ADJUSTL(Vinfo(3)))
                      WRITE (Vname(4,varid),'(a,a)')                    &
     &                      TRIM(Vname(1,varid)), ', scalar, series'
                      WRITE (Vname(5,varid),'(a)')                      &
     &                      TRIM(ADJUSTL(Vinfo(5)))
                    END DO
                  END IF
                CASE ('mudfrac_')
                  IF (NCS.gt.0) THEN
                    varid=varid-1
                    DO i=1,NCS
                      varid=varid+1
                      idfrac(i)=varid
                      DO ng=1,Ngrids
                        Fscale(varid,ng)=scale
                        Iinfo(1,varid,ng)=gtype
                      END DO
                      WRITE (Vname(1,varid),'(a,i2.2)')                 &
     &                      TRIM(ADJUSTL(Vinfo(1))), i
                      WRITE (Vname(2,varid),'(a,a,i2.2)')               &
     &                      TRIM(ADJUSTL(Vinfo(2))), ', size class ', i
                      WRITE (Vname(3,varid),'(a)')                      &
     &                      TRIM(ADJUSTL(Vinfo(3)))
                      WRITE (Vname(4,varid),'(a,a)')                    &
     &                      TRIM(Vname(1,varid)), ', scalar, series'
                      WRITE (Vname(5,varid),'(a)')                      &
     &                      TRIM(ADJUSTL(Vinfo(5)))
                    END DO
                  END IF
                CASE ('sandfrac_')
                  IF (NNS.gt.0) THEN
                    varid=varid-1
                    DO i=1,NNS
                      varid=varid+1
                      idfrac(NCS+i)=varid
                      DO ng=1,Ngrids
                        Fscale(varid,ng)=scale
                        Iinfo(1,varid,ng)=gtype
                      END DO
                      WRITE (Vname(1,varid),'(a,i2.2)')                 &
     &                      TRIM(ADJUSTL(Vinfo(1))), i
                      WRITE (Vname(2,varid),'(a,a,i2.2)')               &
     &                      TRIM(ADJUSTL(Vinfo(2))), ', size class ', i
                      WRITE (Vname(3,varid),'(a)')                      &
     &                      TRIM(ADJUSTL(Vinfo(3)))
                      WRITE (Vname(4,varid),'(a,a)')                    &
     &                      TRIM(Vname(1,varid)), ', scalar, series'
                      WRITE (Vname(5,varid),'(a)')                      &
     &                      TRIM(ADJUSTL(Vinfo(5)))
                    END DO
                  END IF
                CASE ('mudmass_')
                  IF (NCS.gt.0) THEN
                    varid=varid-1
                    DO i=1,NCS
                      varid=varid+1
                      idBmas(i)=varid
                      DO ng=1,Ngrids
                        Fscale(varid,ng)=scale
                        Iinfo(1,varid,ng)=gtype
                      END DO
                      WRITE (Vname(1,varid),'(a,i2.2)')                 &
     &                      TRIM(ADJUSTL(Vinfo(1))), i
                      WRITE (Vname(2,varid),'(a,a,i2.2)')               &
     &                      TRIM(ADJUSTL(Vinfo(2))), ', size class ', i
                      WRITE (Vname(3,varid),'(a)')                      &
     &                      TRIM(ADJUSTL(Vinfo(3)))
                      WRITE (Vname(4,varid),'(a,a)')                    &
     &                      TRIM(Vname(1,varid)), ', scalar, series'
                      WRITE (Vname(5,varid),'(a)')                      &
     &                      TRIM(ADJUSTL(Vinfo(5)))
                    END DO
                  END IF
                CASE ('sandmass_')
                  IF (NNS.gt.0) THEN
                    varid=varid-1
                    DO i=1,NNS
                      varid=varid+1
                      idBmas(NCS+i)=varid
                      DO ng=1,Ngrids
                        Fscale(varid,ng)=scale
                        Iinfo(1,varid,ng)=gtype
                      END DO
                      WRITE (Vname(1,varid),'(a,i2.2)')                 &
     &                      TRIM(ADJUSTL(Vinfo(1))), i
                      WRITE (Vname(2,varid),'(a,a,i2.2)')               &
     &                      TRIM(ADJUSTL(Vinfo(2))), ', size class ', i
                      WRITE (Vname(3,varid),'(a)')                      &
     &                      TRIM(ADJUSTL(Vinfo(3)))
                      WRITE (Vname(4,varid),'(a,a)')                    &
     &                      TRIM(Vname(1,varid)), ', scalar, series'
                      WRITE (Vname(5,varid),'(a)')                      &
     &                      TRIM(ADJUSTL(Vinfo(5)))
                    END DO
                  END IF
                CASE ('bedload_Umud_')
                  IF (NCS.gt.0) THEN
                    varid=varid-1
                    DO i=1,NCS
                      varid=varid+1
                      idUbld(i)=varid
                      DO ng=1,Ngrids
                        Fscale(varid,ng)=scale
                        Iinfo(1,varid,ng)=gtype
                      END DO
                      WRITE (Vname(1,varid),'(a,i2.2)')                 &
     &                      TRIM(ADJUSTL(Vinfo(1))), i
                      WRITE (Vname(2,varid),'(a,a,i2.2)')               &
     &                      TRIM(ADJUSTL(Vinfo(2))), ', size class ', i
                      WRITE (Vname(3,varid),'(a)')                      &
     &                      TRIM(ADJUSTL(Vinfo(3)))
                      WRITE (Vname(4,varid),'(a,a)')                    &
     &                      TRIM(Vname(1,varid)), ', scalar, series'
                      WRITE (Vname(5,varid),'(a)')                      &
     &                      TRIM(ADJUSTL(Vinfo(5)))
                    END DO
                  END IF
                CASE ('bedload_Usand_')
                  IF (NNS.gt.0) THEN
                    varid=varid-1
                    DO i=1,NNS
                      varid=varid+1
                      idUbld(NCS+i)=varid
                      DO ng=1,Ngrids
                        Fscale(varid,ng)=scale
                        Iinfo(1,varid,ng)=gtype
                      END DO
                      WRITE (Vname(1,varid),'(a,i2.2)')                 &
     &                      TRIM(ADJUSTL(Vinfo(1))), i
                      WRITE (Vname(2,varid),'(a,a,i2.2)')               &
     &                      TRIM(ADJUSTL(Vinfo(2))), ', size class ', i
                      WRITE (Vname(3,varid),'(a)')                      &
     &                      TRIM(ADJUSTL(Vinfo(3)))
                      WRITE (Vname(4,varid),'(a,a)')                    &
     &                      TRIM(Vname(1,varid)), ', scalar, series'
                      WRITE (Vname(5,varid),'(a)')                      &
     &                      TRIM(ADJUSTL(Vinfo(5)))
                    END DO
                  END IF
                CASE ('bedload_Vmud_')
                  IF (NCS.gt.0) THEN
                    varid=varid-1
                    DO i=1,NCS
                      varid=varid+1
                      idVbld(i)=varid
                      DO ng=1,Ngrids
                        Fscale(varid,ng)=scale
                        Iinfo(1,varid,ng)=gtype
                      END DO
                      WRITE (Vname(1,varid),'(a,i2.2)')                 &
     &                      TRIM(ADJUSTL(Vinfo(1))), i
                      WRITE (Vname(2,varid),'(a,a,i2.2)')               &
     &                      TRIM(ADJUSTL(Vinfo(2))), ', size class ', i
                      WRITE (Vname(3,varid),'(a)')                      &
     &                      TRIM(ADJUSTL(Vinfo(3)))
                      WRITE (Vname(4,varid),'(a,a)')                    &
     &                      TRIM(Vname(1,varid)), ', scalar, series'
                      WRITE (Vname(5,varid),'(a)')                      &
     &                      TRIM(ADJUSTL(Vinfo(5)))
                    END DO
                  END IF
                CASE ('bedload_Vsand_')
                  IF (NNS.gt.0) THEN
                    varid=varid-1
                    DO i=1,NNS
                      varid=varid+1
                      idVbld(NCS+i)=varid
                      DO ng=1,Ngrids
                        Fscale(varid,ng)=scale
                        Iinfo(1,varid,ng)=gtype
                      END DO
                      WRITE (Vname(1,varid),'(a,i2.2)')                 &
     &                      TRIM(ADJUSTL(Vinfo(1))), i
                      WRITE (Vname(2,varid),'(a,a,i2.2)')               &
     &                      TRIM(ADJUSTL(Vinfo(2))), ', size class ', i
                      WRITE (Vname(3,varid),'(a)')                      &
     &                      TRIM(ADJUSTL(Vinfo(3)))
                      WRITE (Vname(4,varid),'(a,a)')                    &
     &                      TRIM(Vname(1,varid)), ', scalar, series'
                      WRITE (Vname(5,varid),'(a)')                      &
     &                      TRIM(ADJUSTL(Vinfo(5)))
                    END DO
                  END IF
                CASE ('river_mud_')
                  IF (NCS.gt.0) THEN
                    varid=varid-1
                    DO i=1,NCS
                      varid=varid+1
                      idRtrc(idmud(i))=varid
                      DO ng=1,Ngrids
                        Fscale(varid,ng)=scale
                        Iinfo(1,varid,ng)=gtype
                      END DO
                      WRITE (Vname(1,varid),'(a,i2.2)')                 &
     &                      TRIM(ADJUSTL(Vinfo(1))), i
                      WRITE (Vname(2,varid),'(a,a,i2.2)')               &
     &                      TRIM(ADJUSTL(Vinfo(2))), ', size class ', i
                      WRITE (Vname(3,varid),'(a)')                      &
     &                      TRIM(ADJUSTL(Vinfo(3)))
                      WRITE (Vname(4,varid),'(a,a)')                    &
     &                      TRIM(Vname(1,varid)), ', scalar, series'
                      WRITE (Vname(5,varid),'(a)')                      &
     &                      TRIM(ADJUSTL(Vinfo(5)))
                    END DO
                  END IF
                CASE ('river_sand_')
                  IF (NNS.gt.0) THEN
                    varid=varid-1
                    DO i=1,NNS
                      varid=varid+1
                      idRtrc(isand(i))=varid
                      DO ng=1,Ngrids
                        Fscale(varid,ng)=scale
                        Iinfo(1,varid,ng)=gtype
                      END DO
                      WRITE (Vname(1,varid),'(a,i2.2)')                 &
     &                      TRIM(ADJUSTL(Vinfo(1))), i
                      WRITE (Vname(2,varid),'(a,a,i2.2)')               &
     &                      TRIM(ADJUSTL(Vinfo(2))), ', size class ', i
                      WRITE (Vname(3,varid),'(a)')                      &
     &                      TRIM(ADJUSTL(Vinfo(3)))
                      WRITE (Vname(4,varid),'(a,a)')                    &
     &                      TRIM(Vname(1,varid)), ', scalar, series'
                      WRITE (Vname(5,varid),'(a)')                      &
     &                      TRIM(ADJUSTL(Vinfo(5)))
                    END DO
                  END IF
              END SELECT
#endif
#ifdef ECOSIM
!
!  Adjust information for all EcoSim tracers.
!
              SELECT CASE (Vinfo(1))
                CASE ('Bac_C')
                  varid=varid-1
                  DO i=1,Nbac
                    varid=varid+1
                    idTvar(iBacC(i))=varid
                    DO ng=1,Ngrids
                      Fscale(varid,ng)=scale
                      Iinfo(1,varid,ng)=gtype
                    END DO
                    WRITE (Vname(1,varid),'(a,i1)') 'Bac_C', i
                    WRITE (Vname(2,varid),'(a,i1)')                     &
     &                    'bacteria, carbon group ', i
                    WRITE (Vname(3,varid),'(a)')                        &
     &                    TRIM(ADJUSTL(Vinfo(3)))
                    WRITE (Vname(4,varid),'(a,a)')                      &
     &                    TRIM(Vname(1,varid)), ', scalar, series'
                    WRITE (Vname(5,varid),'(a)')                        &
     &                    TRIM(ADJUSTL(Vinfo(5)))
                  END DO
                CASE ('Bac_F')
                  varid=varid-1
                  DO i=1,Nbac
                    varid=varid+1
                    idTvar(iBacF(i))=varid
                    DO ng=1,Ngrids
                      Fscale(varid,ng)=scale
                      Iinfo(1,varid,ng)=gtype
                    END DO
                    WRITE (Vname(1,varid),'(a,i1)') 'Bac_F', i
                    WRITE (Vname(2,varid),'(a,i1)')                     &
     &                    'bacteria, iron group ', i
                    WRITE (Vname(3,varid),'(a)')                        &
     &                    TRIM(ADJUSTL(Vinfo(3)))
                    WRITE (Vname(4,varid),'(a,a)')                      &
     &                    TRIM(Vname(1,varid)), ', scalar, series'
                    WRITE (Vname(5,varid),'(a)')                        &
     &                    TRIM(ADJUSTL(Vinfo(5)))
                  END DO
                CASE ('Bac_N')
                  varid=varid-1
                  DO i=1,Nbac
                    varid=varid+1
                    idTvar(iBacN(i))=varid
                    DO ng=1,Ngrids
                      Fscale(varid,ng)=scale
                      Iinfo(1,varid,ng)=gtype
                    END DO
                    WRITE (Vname(1,varid),'(a,i1)') 'Bac_N', i
                    WRITE (Vname(2,varid),'(a,i1)')                     &
     &                    'bacteria, nitrogen group ', i
                    WRITE (Vname(3,varid),'(a)')                        &
     &                    TRIM(ADJUSTL(Vinfo(3)))
                    WRITE (Vname(4,varid),'(a,a)')                      &
     &                    TRIM(Vname(1,varid)), ', scalar, series'
                    WRITE (Vname(5,varid),'(a)')                        &
     &                    TRIM(ADJUSTL(Vinfo(5)))
                  END DO
                CASE ('Bac_P')
                  varid=varid-1
                  DO i=1,Nbac
                    varid=varid+1
                    idTvar(iBacP(i))=varid
                    DO ng=1,Ngrids
                      Fscale(varid,ng)=scale
                      Iinfo(1,varid,ng)=gtype
                    END DO
                    WRITE (Vname(1,varid),'(a,i1)') 'Bac_P', i
                    WRITE (Vname(2,varid),'(a,i1)')                     &
     &                    'bacteria, phosphorus group ', i
                    WRITE (Vname(3,varid),'(a)')                        &
     &                    TRIM(ADJUSTL(Vinfo(3)))
                    WRITE (Vname(4,varid),'(a,a)')                      &
     &                    TRIM(Vname(1,varid)), ', scalar, series'
                    WRITE (Vname(5,varid),'(a)')                        &
     &                    TRIM(ADJUSTL(Vinfo(5)))
                  END DO
                CASE ('CDM_C')
                  varid=varid-1
                  DO i=1,Ndom
                    varid=varid+1
                    idTvar(iCDMC(i))=varid
                    DO ng=1,Ngrids
                      Fscale(varid,ng)=scale
                      Iinfo(1,varid,ng)=gtype
                    END DO
                    WRITE (Vname(1,varid),'(a,i1)') 'CDM_C', i
                    WRITE (Vname(2,varid),'(a,i1)')                     &
     &                    'color degradational matter, carbon group ', i
                    WRITE (Vname(3,varid),'(a)')                        &
     &                    TRIM(ADJUSTL(Vinfo(3)))
                    WRITE (Vname(4,varid),'(a,a)')                      &
     &                    TRIM(Vname(1,varid)), ', scalar, series'
                    WRITE (Vname(5,varid),'(a)')                        &
     &                    TRIM(ADJUSTL(Vinfo(5)))
                  END DO
                CASE ('DOM_C')
                  varid=varid-1
                  DO i=1,Ndom
                    varid=varid+1
                    idTvar(iDOMC(i))=varid
                    DO ng=1,Ngrids
                      Fscale(varid,ng)=scale
                      Iinfo(1,varid,ng)=gtype
                    END DO
                    WRITE (Vname(1,varid),'(a,i1)') 'DOM_C', i
                    WRITE (Vname(2,varid),'(a,i1)')                     &
     &                    'dissolved organic matter, carbon group ', i
                    WRITE (Vname(3,varid),'(a)')                        &
     &                    TRIM(ADJUSTL(Vinfo(3)))
                    WRITE (Vname(4,varid),'(a,a)')                      &
     &                    TRIM(Vname(1,varid)), ', scalar, series'
                    WRITE (Vname(5,varid),'(a)')                        &
     &                    TRIM(ADJUSTL(Vinfo(5)))
                  END DO
                CASE ('DOM_N')
                  varid=varid-1
                  DO i=1,Ndom
                    varid=varid+1
                    idTvar(iDOMN(i))=varid
                    DO ng=1,Ngrids
                      Fscale(varid,ng)=scale
                      Iinfo(1,varid,ng)=gtype
                    END DO
                    WRITE (Vname(1,varid),'(a,i1)') 'DOM_N', i
                    WRITE (Vname(2,varid),'(a,i1)')                     &
     &                    'dissolved organic matter, nitrogen group ', i
                    WRITE (Vname(3,varid),'(a)')                        &
     &                    TRIM(ADJUSTL(Vinfo(3)))
                    WRITE (Vname(4,varid),'(a,a)')                      &
     &                    TRIM(Vname(1,varid)), ', scalar, series'
                    WRITE (Vname(5,varid),'(a)')                        &
     &                    TRIM(ADJUSTL(Vinfo(5)))
                  END DO
                CASE ('DOM_P')
                  varid=varid-1
                  DO i=1,Ndom
                    varid=varid+1
                    idTvar(iDOMP(i))=varid
                    DO ng=1,Ngrids
                      Fscale(varid,ng)=scale
                      Iinfo(1,varid,ng)=gtype
                    END DO
                    WRITE (Vname(1,varid),'(a,i1)') 'DOM_P', i
                    WRITE (Vname(2,varid),'(a,i1)')                     &
     &                  'dissolved organic matter, phosphorus group ', i
                    WRITE (Vname(3,varid),'(a)')                        &
     &                    TRIM(ADJUSTL(Vinfo(3)))
                    WRITE (Vname(4,varid),'(a,a)')                      &
     &                    TRIM(Vname(1,varid)), ', scalar, series'
                    WRITE (Vname(5,varid),'(a)')                        &
     &                    TRIM(ADJUSTL(Vinfo(5)))
                  END DO
                CASE ('Fec_C')
                  varid=varid-1
                  DO i=1,Nfec
                    varid=varid+1
                    idTvar(iFecC(i))=varid
                    DO ng=1,Ngrids
                      Fscale(varid,ng)=scale
                      Iinfo(1,varid,ng)=gtype
                    END DO
                    WRITE (Vname(1,varid),'(a,i1)') 'Fec_C', i
                    WRITE (Vname(2,varid),'(a,i1)')                     &
     &                    'fecal matter, carbon group ', i
                    WRITE (Vname(3,varid),'(a)')                        &
     &                    TRIM(ADJUSTL(Vinfo(3)))
                    WRITE (Vname(4,varid),'(a,a)')                      &
     &                    TRIM(Vname(1,varid)), ', scalar, series'
                    WRITE (Vname(5,varid),'(a)')                        &
     &                    TRIM(ADJUSTL(Vinfo(5)))
                  END DO
                CASE ('Fec_F')
                  varid=varid-1
                  DO i=1,Nfec
                    varid=varid+1
                    idTvar(iFecF(i))=varid
                    DO ng=1,Ngrids
                      Fscale(varid,ng)=scale
                      Iinfo(1,varid,ng)=gtype
                    END DO
                    WRITE (Vname(1,varid),'(a,i1)') 'Fec_F', i
                    WRITE (Vname(2,varid),'(a,i1)')                     &
     &                    'fecal matter, iron group ', i
                    WRITE (Vname(3,varid),'(a)')                        &
     &                    TRIM(ADJUSTL(Vinfo(3)))
                    WRITE (Vname(4,varid),'(a,a)')                      &
     &                    TRIM(Vname(1,varid)), ', scalar, series'
                    WRITE (Vname(5,varid),'(a)')                        &
     &                    TRIM(ADJUSTL(Vinfo(5)))
                  END DO
                CASE ('Fec_N')
                  varid=varid-1
                  DO i=1,Nfec
                    varid=varid+1
                    idTvar(iFecN(i))=varid
                    DO ng=1,Ngrids
                      Fscale(varid,ng)=scale
                      Iinfo(1,varid,ng)=gtype
                    END DO
                    WRITE (Vname(1,varid),'(a,i1)') 'Fec_N', i
                    WRITE (Vname(2,varid),'(a,i1)')                     &
     &                    'fecal matter, nitrogen group ', i
                    WRITE (Vname(3,varid),'(a)')                        &
     &                    TRIM(ADJUSTL(Vinfo(3)))
                    WRITE (Vname(4,varid),'(a,a)')                      &
     &                    TRIM(Vname(1,varid)), ', scalar, series'
                    WRITE (Vname(5,varid),'(a)')                        &
     &                    TRIM(ADJUSTL(Vinfo(5)))
                  END DO
                CASE ('Fec_P')
                  varid=varid-1
                  DO i=1,Nfec
                    varid=varid+1
                    idTvar(iFecP(i))=varid
                    DO ng=1,Ngrids
                      Fscale(varid,ng)=scale
                      Iinfo(1,varid,ng)=gtype
                    END DO
                    WRITE (Vname(1,varid),'(a,i1)') 'Fec_P', i
                    WRITE (Vname(2,varid),'(a,i1)')                     &
     &                    'fecal matter, phosphorus group ', i
                    WRITE (Vname(3,varid),'(a)')                        &
     &                    TRIM(ADJUSTL(Vinfo(3)))
                    WRITE (Vname(4,varid),'(a,a)')                      &
     &                    TRIM(Vname(1,varid)), ', scalar, series'
                    WRITE (Vname(5,varid),'(a)')                        &
     &                    TRIM(ADJUSTL(Vinfo(5)))
                  END DO
                CASE ('Fec_S')
                  varid=varid-1
                  DO i=1,Nfec
                    varid=varid+1
                    idTvar(iFecS(i))=varid
                    DO ng=1,Ngrids
                      Fscale(varid,ng)=scale
                      Iinfo(1,varid,ng)=gtype
                    END DO
                    WRITE (Vname(1,varid),'(a,i1)') 'Fec_S', i
                    WRITE (Vname(2,varid),'(a,i1)')                     &
     &                    'fecal matter, silica group ', i
                    WRITE (Vname(3,varid),'(a)')                        &
     &                    TRIM(ADJUSTL(Vinfo(3)))
                    WRITE (Vname(4,varid),'(a,a)')                      &
     &                    TRIM(Vname(1,varid)), ', scalar, series'
                    WRITE (Vname(5,varid),'(a)')                        &
     &                    TRIM(ADJUSTL(Vinfo(5)))
                  END DO
                CASE ('Phy_C')
                  varid=varid-1
                  DO i=1,Nphy
                    varid=varid+1
                    idTvar(iPhyC(i))=varid
                    DO ng=1,Ngrids
                      Fscale(varid,ng)=scale
                      Iinfo(1,varid,ng)=gtype
                    END DO
                    WRITE (Vname(1,varid),'(a,i1)') 'Phy_C', i
                    WRITE (Vname(2,varid),'(a,a)')                      &
     &                    TRIM(PhyName(i)), ', carbon group'
                    WRITE (Vname(3,varid),'(a)')                        &
     &                    TRIM(ADJUSTL(Vinfo(3)))
                    WRITE (Vname(4,varid),'(a,a)')                      &
     &                    TRIM(Vname(1,varid)), ', scalar, series'
                    WRITE (Vname(5,varid),'(a)')                        &
     &                    TRIM(ADJUSTL(Vinfo(5)))
                  END DO
                CASE ('Phy_F')
                  varid=varid-1
                  DO i=1,Nphy
                    varid=varid+1
                    idTvar(iPhyF(i))=varid
                    DO ng=1,Ngrids
                      Fscale(varid,ng)=scale
                      Iinfo(1,varid,ng)=gtype
                    END DO
                    WRITE (Vname(1,varid),'(a,i1)') 'Phy_F', i
                    WRITE (Vname(2,varid),'(a,a)')                      &
     &                    TRIM(PhyName(i)), ', iron group'
                    WRITE (Vname(3,varid),'(a)')                        &
     &                    TRIM(ADJUSTL(Vinfo(3)))
                    WRITE (Vname(4,varid),'(a,a)')                      &
     &                    TRIM(Vname(1,varid)), ', scalar, series'
                    WRITE (Vname(5,varid),'(a)')                        &
     &                    TRIM(ADJUSTL(Vinfo(5)))
                  END DO
                CASE ('Phy_N')
                  varid=varid-1
                  DO i=1,Nphy
                    varid=varid+1
                    idTvar(iPhyN(i))=varid
                    DO ng=1,Ngrids
                      Fscale(varid,ng)=scale
                      Iinfo(1,varid,ng)=gtype
                    END DO
                    WRITE (Vname(1,varid),'(a,i1)') 'Phy_N', i
                    WRITE (Vname(2,varid),'(a,a)')                      &
     &                    TRIM(PhyName(i)), ', nitrogen group'
                    WRITE (Vname(3,varid),'(a)')                        &
     &                    TRIM(ADJUSTL(Vinfo(3)))
                    WRITE (Vname(4,varid),'(a,a)')                      &
     &                    TRIM(Vname(1,varid)), ', scalar, series'
                    WRITE (Vname(5,varid),'(a)')                        &
     &                    TRIM(ADJUSTL(Vinfo(5)))
                  END DO
                CASE ('Phy_P')
                  varid=varid-1
                  DO i=1,Nphy
                    varid=varid+1
                    idTvar(iPhyP(i))=varid
                    DO ng=1,Ngrids
                      Fscale(varid,ng)=scale
                      Iinfo(1,varid,ng)=gtype
                    END DO
                    WRITE (Vname(1,varid),'(a,i1)') 'Phy_P', i
                    WRITE (Vname(2,varid),'(a,a)')                      &
     &                    TRIM(PhyName(i)), ', phosphorus group'
                    WRITE (Vname(3,varid),'(a)')                        &
     &                    TRIM(ADJUSTL(Vinfo(3)))
                    WRITE (Vname(4,varid),'(a,a)')                      &
     &                    TRIM(Vname(1,varid)), ', scalar, series'
                    WRITE (Vname(5,varid),'(a)')                        &
     &                    TRIM(ADJUSTL(Vinfo(5)))
                  END DO
                CASE ('Phy_S')
                  j=0
                  DO i=1,Nphy
                    IF (PHY(i).le.2) j=j+1
                  END DO                
                  IF (j.gt.0) varid=varid-1
                  DO i=1,Nphy
                    IF (PHY(i).le.2) THEN
                      varid=varid+1
                      idTvar(iPhyS(i))=varid
                      DO ng=1,Ngrids
                        Fscale(varid,ng)=scale
                        Iinfo(1,varid,ng)=gtype
                      END DO
                      WRITE (Vname(1,varid),'(a,i1)') 'Phy_S', i
                      WRITE (Vname(2,varid),'(a,a)')                    &
     &                      TRIM(PhyName(i)), ', silica group'
                      WRITE (Vname(3,varid),'(a)')                      &
     &                      TRIM(ADJUSTL(Vinfo(3)))
                      WRITE (Vname(4,varid),'(a,a)')                    &
     &                      TRIM(Vname(1,varid)), ', scalar, series'
                      WRITE (Vname(5,varid),'(a)')                      &
     &                      TRIM(ADJUSTL(Vinfo(5)))
                    END IF
                  END DO
                CASE ('Pigs_XY')
                  varid=varid-1
                  DO j=1,Npig 
                    DO i=1,Nphy
                      IF (PIG(PHY(i),j).gt.0) THEN
                        varid=varid+1
                        idTvar(iPigs(i,j))=varid
                        DO ng=1,Ngrids
                          Fscale(varid,ng)=scale
                          Iinfo(1,varid,ng)=gtype
                        END DO
                        WRITE (Vname(1,varid),'(a,i1,i1)') 'Pigs_', i, j
                        WRITE (Vname(2,varid),'(a,a,a)')                &
     &                        TRIM(PhyName(i)), ', ', TRIM(PigName(j))
                        WRITE (Vname(3,varid),'(a)')                    &
     &                        TRIM(ADJUSTL(Vinfo(3)))
                        WRITE (Vname(4,varid),'(a,a)')                  &
     &                        TRIM(Vname(1,varid)), ', scalar, series'
                        WRITE (Vname(5,varid),'(a)')                    &
     &                        TRIM(ADJUSTL(Vinfo(5)))
                      END IF
                    END DO
                  END DO
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
                    WRITE (Vname(5,varid),'(a)')                      &
     &                    TRIM(ADJUSTL(Vinfo(5)))
                  END DO
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
!
!-----------------------------------------------------------------------
!  Set model state variables indices.
!-----------------------------------------------------------------------
!
        idSvar(isFsur)=idFsur
        idSvar(isUbar)=idUbar
        idSvar(isVbar)=idVbar
#ifdef SOLVE3D
        idSvar(isUvel)=idUvel
        idSvar(isVvel)=idVvel
        DO i=1,MT
          idSvar(isTvar(i))=idTvar(i)
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
