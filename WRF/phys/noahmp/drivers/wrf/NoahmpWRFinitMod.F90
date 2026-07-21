module NoahmpWRFinitMod

! --------------------------------------------------------------------------
! this is for NoahmpIO variable mapping & initialization in WRF physics_init
! and calling the main NoahMP initialization module: NoahmpInitMain(NoahmpIO) 
! adapted from NOAHMP_INIT in original module_sf_noahmpdrv.F file
!
! Coder: Cenlin He (NCAR), December 2025
! ---------------------------------------------------------------------------

contains

  subroutine NoahmpWRFinit(NoahmpIO, MMINLU, SNOW, SNOWH, CANWAT, ISLTYP, IVGTYP, XLAT, &
                   TSLB,  SMOIS, SH2O,   DZS, FNDSOILW, FNDSNOWH,                       &
                   TSK, isnowxy, tvxy,  tgxy, canicexy,      TMN,   XICE,               &
                   canliqxy,    eahxy, tahxy,     cmxy,     chxy,                       &
                   fwetxy, sneqvoxy, alboldxy, qsnowxy, qrainxy, wslakexy, zwtxy, waxy, &
                   wtxy, tsnoxy, zsnsoxy, snicexy, snliqxy, lfmassxy, rtmassxy,         &
                   stmassxy, woodxy, stblcpxy, fastcpxy, xsaixy, lai,                   &
                   grainxy,   gddxy,                                                    &
                   croptype, cropcat,                                                   &
                   irnumsi, irnummi, irnumfi, irwatsi,                                  &
                   irwatmi, irwatfi, ireloss, irsivol,                                  &
                   irmivol, irfivol, irrsplh,                                           &
                   t2mvxy,   t2mbxy, chstarxy, fsatxy, wsurfxy,                         &
                   snrdsxy, snfrxy, bcphixy, bcphoxy, ocphixy, ocphoxy, dust1xy,        &
                   dust2xy, dust3xy, dust4xy, dust5xy, massconcbcphixy, massconcbcphoxy,&
                   massconcocphixy, massconcocphoxy, massconcdust1xy, massconcdust2xy,  &
                   massconcdust3xy, massconcdust4xy, massconcdust5xy,                   &
                   ALBSOILDIRXY, ALBSOILDIFXY,                                          &
                   NSOIL,   restart,                                                    &
                   allowed_to_read , IOPT_RUNSUB, IOPT_CROP, IOPT_IRR, IOPT_IRRM,       &
                   SF_URBAN_PHYSICS, IOPT_SOIL, IOPT_ALB, IOPT_WETLAND,                 &
                   SNICAR_SNOWOPTICS_OPT, SNICAR_DUSTOPTICS_OPT, SNICAR_SOLARSPEC_OPT,  & ! optional SNICAR option
                   SNICAR_BANDNUMBER_OPT,                                               & ! optional SNICAR option
                   ids,ide, jds,jde, kds,kde,                                           &
                   ims,ime, jms,jme, kms,kme,                                           &
                   its,ite, jts,jte, kts,kte,                                           &
                   smoiseq,smcwtdxy,rechxy,deeprechxy,qtdrain,areaxy,dx,dy,msftx,msfty, & ! Optional groundwater
                   wtddt,   stepwtd, dt, qrfsxy, qspringsxy, qslatxy,                   & ! Optional groundwater
                   fdepthxy, ht, riverbedxy, eqzwt, rivercondxy, pexpxy, rechclim       ) ! Optional groundwater

! ---------------------------------------------------------------------------------------

    use NoahmpIOVarType, only : NoahmpIO_type
    use NoahmpIOVarInitMod
    use NoahmpReadTableMod
    use NoahmpInitMainMod
    use SnowInputSnicarMod

    implicit none

    type(NoahmpIO_type), intent(inout)                         :: NoahmpIO

    ! input only
    INTEGER, INTENT(IN)                                        :: ids,ide, jds,jde, kds,kde,  &
                                                                  ims,ime, jms,jme, kms,kme,  &
                                                                  its,ite, jts,jte, kts,kte
    INTEGER, INTENT(IN)                                        :: NSOIL, IOPT_RUNSUB, IOPT_CROP, &
                                                                  IOPT_IRR, IOPT_IRRM, IOPT_SOIL,&
                                                                  IOPT_ALB, IOPT_WETLAND
    LOGICAL, INTENT(IN)                                        :: restart, allowed_to_read
    LOGICAL, INTENT(IN)                                        :: FNDSOILW, FNDSNOWH
    INTEGER, INTENT(IN)                                        :: SF_URBAN_PHYSICS                      
    CHARACTER(LEN=*),                    INTENT(IN)            :: MMINLU
    REAL,    DIMENSION(NSOIL), INTENT(IN)                      :: DZS                 ! Thickness of the soil layers [m]
    INTEGER, DIMENSION(ims:ime,jms:jme), INTENT(IN)            :: ISLTYP, IVGTYP
    REAL   , DIMENSION(ims:ime,5,jms:jme),INTENT(IN)           :: croptype            ! crop type fraction
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(IN)            :: XLAT                ! latitude
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(IN)            :: TSK                 ! skin temperature (k)
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(IN)            :: XICE                ! sea ice fraction
    REAL,                                INTENT(IN), OPTIONAL  :: DT, WTDDT
    REAL,                                INTENT(IN), OPTIONAL  :: DX, DY
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(IN), OPTIONAL  :: FDEPTHXY            ! efolding depth for transmissivity (m)
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(IN), OPTIONAL  :: HT                  ! terrain height (m)
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(IN), OPTIONAL  :: MSFTX, MSFTY
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(IN), OPTIONAL  :: rechclim
    INTEGER, INTENT(IN),                             OPTIONAL  :: SNICAR_SNOWOPTICS_OPT, &
                                                                  SNICAR_DUSTOPTICS_OPT, &
                                                                  SNICAR_SOLARSPEC_OPT,  &
                                                                  SNICAR_BANDNUMBER_OPT
    ! in/out
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: TMN                 ! deep soil temperature (k)
    INTEGER, DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: isnowxy             ! actual no. of snow layers
    REAL,    DIMENSION(ims:ime,1:NSOIL,jms:jme), INTENT(INOUT) :: SMOIS, SH2O, TSLB
    REAL,    DIMENSION(ims:ime, jms:jme), INTENT(INOUT)        :: SNOW, SNOWH, CANWAT
    REAL,    DIMENSION(ims:ime,-2:NSOIL,jms:jme),INTENT(INOUT) :: zsnsoxy             ! snow layer depth [m]
    REAL,    DIMENSION(ims:ime,-2:0,jms:jme), INTENT(INOUT)    :: tsnoxy              ! snow temperature [K]
    REAL,    DIMENSION(ims:ime,-2:0,jms:jme), INTENT(INOUT)    :: snicexy             ! snow layer ice [mm]
    REAL,    DIMENSION(ims:ime,-2:0,jms:jme), INTENT(INOUT)    :: snliqxy             ! snow layer liquid water [mm]
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: tvxy                ! vegetation canopy temperature
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: tgxy                ! ground surface temperature
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: canicexy            ! canopy-intercepted ice (mm)
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: canliqxy            ! canopy-intercepted liquid water (mm)
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: eahxy               ! canopy air vapor pressure (pa)
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: tahxy               ! canopy air temperature (k)
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: cmxy                ! momentum drag coefficient
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: chxy                ! sensible heat exchange coefficient
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: fwetxy              ! wetted or snowed fraction of the canopy (-)
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: sneqvoxy            ! snow mass at last time step(mm h2o)
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: alboldxy            ! snow albedo at last time step (-)
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: qsnowxy             ! snowfall on the ground [mm/s]
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: qrainxy             ! rainfall on the ground [mm/s]
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: wslakexy            ! lake water storage [mm]
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: zwtxy               ! water table depth [m]
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: waxy                ! water in the "aquifer" [mm]
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: wtxy                ! groundwater storage [mm]
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: lfmassxy            ! leaf mass [g/m2]
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: rtmassxy            ! mass of fine roots [g/m2]
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: stmassxy            ! stem mass [g/m2]
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: woodxy              ! mass of wood (incl. woody roots) [g/m2]
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: grainxy             ! mass of grain [g/m2] !XING
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: gddxy               ! growing degree days !XING
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: stblcpxy            ! stable carbon in deep soil [g/m2]
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: fastcpxy            ! short-lived carbon, shallow soil [g/m2]
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: xsaixy              ! stem area index
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: lai                 ! leaf area index
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: qtdrain             ! tile drainage (mm)
    INTEGER, DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: irnumsi             ! irrigation number
    INTEGER, DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: irnummi             ! irrigation number
    INTEGER, DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: irnumfi             ! irrigation number
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: irwatsi             ! irrigation water amount
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: irwatmi             ! irrigation water amount
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: irwatfi             ! irrigation water amount
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: ireloss             ! irrigation loss
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: irsivol             ! irrigation water volume
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: irmivol             ! irrigation water volume
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: irfivol             ! irrigation water volume
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: irrsplh             ! irrigation evaporation heat
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: t2mvxy              ! 2m temperature vegetation part (k)
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: t2mbxy              ! 2m temperature bare ground part (k)
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: chstarxy            ! dummy
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: fsatxy              ! saturation fraction of grid (-)
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: wsurfxy             ! wetland water storage (mm)
    REAL,    DIMENSION(ims:ime,-2:0,jms:jme), INTENT(INOUT)    :: snrdsxy             ! SNICAR snow radius
    REAL,    DIMENSION(ims:ime,-2:0,jms:jme), INTENT(INOUT)    :: snfrxy              ! SNICAR snow freezing rate
    REAL,    DIMENSION(ims:ime,-2:0,jms:jme), INTENT(INOUT)    :: bcphixy             ! SNICAR BCPHI mass in snow
    REAL,    DIMENSION(ims:ime,-2:0,jms:jme), INTENT(INOUT)    :: bcphoxy             ! SNICAR BCPHO mass in snow
    REAL,    DIMENSION(ims:ime,-2:0,jms:jme), INTENT(INOUT)    :: ocphixy             ! SNICAR OCPHI mass in snow
    REAL,    DIMENSION(ims:ime,-2:0,jms:jme), INTENT(INOUT)    :: ocphoxy             ! SNICAR OCPHO mass in snow
    REAL,    DIMENSION(ims:ime,-2:0,jms:jme), INTENT(INOUT)    :: dust1xy             ! SNICAR DUST1 mass in snow
    REAL,    DIMENSION(ims:ime,-2:0,jms:jme), INTENT(INOUT)    :: dust2xy             ! SNICAR DUST2 mass in snow
    REAL,    DIMENSION(ims:ime,-2:0,jms:jme), INTENT(INOUT)    :: dust3xy             ! SNICAR DUST3 mass in snow
    REAL,    DIMENSION(ims:ime,-2:0,jms:jme), INTENT(INOUT)    :: dust4xy             ! SNICAR DUST4 mass in snow
    REAL,    DIMENSION(ims:ime,-2:0,jms:jme), INTENT(INOUT)    :: dust5xy             ! SNICAR DUST5 mass in snow
    REAL,    DIMENSION(ims:ime,-2:0,jms:jme), INTENT(INOUT)    :: massconcbcphixy     ! SNICAR BCPHI mass conc in snow
    REAL,    DIMENSION(ims:ime,-2:0,jms:jme), INTENT(INOUT)    :: massconcbcphoxy     ! SNICAR BCPHO mass conc in snow
    REAL,    DIMENSION(ims:ime,-2:0,jms:jme), INTENT(INOUT)    :: massconcocphixy     ! SNICAR OCPHI mass conc in snow
    REAL,    DIMENSION(ims:ime,-2:0,jms:jme), INTENT(INOUT)    :: massconcocphoxy     ! SNICAR OCPHO mass conc in snow
    REAL,    DIMENSION(ims:ime,-2:0,jms:jme), INTENT(INOUT)    :: massconcdust1xy     ! SNICAR DUST1 mass conc in snow
    REAL,    DIMENSION(ims:ime,-2:0,jms:jme), INTENT(INOUT)    :: massconcdust2xy     ! SNICAR DUST2 mass conc in snow
    REAL,    DIMENSION(ims:ime,-2:0,jms:jme), INTENT(INOUT)    :: massconcdust3xy     ! SNICAR DUST3 mass conc in snow
    REAL,    DIMENSION(ims:ime,-2:0,jms:jme), INTENT(INOUT)    :: massconcdust4xy     ! SNICAR DUST4 mass conc in snow
    REAL,    DIMENSION(ims:ime,-2:0,jms:jme), INTENT(INOUT)    :: massconcdust5xy     ! SNICAR DUST5 mass conc in snow
    REAL,    DIMENSION(ims:ime,1:2,jms:jme),  INTENT(INOUT)    :: ALBSOILDIRXY        ! soil albedo direct
    REAL,    DIMENSION(ims:ime,1:2,jms:jme),  INTENT(INOUT)    :: ALBSOILDIFXY        ! soil albedo diffuse
    REAL,    DIMENSION(ims:ime,jms:jme),         INTENT(INOUT), OPTIONAL :: smcwtdxy     ! deep soil moisture content [m3m-3]
    REAL,    DIMENSION(ims:ime,jms:jme),         INTENT(INOUT), OPTIONAL :: deeprechxy   ! deep recharge [m]
    REAL,    DIMENSION(ims:ime,jms:jme),         INTENT(INOUT), OPTIONAL :: rechxy       ! accumulated recharge [mm]
    REAL,    DIMENSION(ims:ime,jms:jme),         INTENT(INOUT), OPTIONAL :: qrfsxy       ! accumulated flux from groundwater to rivers [mm]
    REAL,    DIMENSION(ims:ime,jms:jme),         INTENT(INOUT), OPTIONAL :: qspringsxy   ! accumulated seeping water [mm]
    REAL,    DIMENSION(ims:ime,jms:jme),         INTENT(INOUT), OPTIONAL :: qslatxy      ! accumulated lateral flow [mm]
    REAL,    DIMENSION(ims:ime,jms:jme),         INTENT(INOUT), OPTIONAL :: areaxy       ! grid cell area [m2]
    REAL,    DIMENSION(ims:ime,jms:jme),         INTENT(INOUT), OPTIONAL :: RIVERBEDXY   ! riverbed depth (m)
    REAL,    DIMENSION(ims:ime,jms:jme),         INTENT(INOUT), OPTIONAL :: EQZWT        ! equilibrium water table depth (m)
    REAL,    DIMENSION(ims:ime,jms:jme),         INTENT(INOUT), OPTIONAL :: RIVERCONDXY  ! river conductance
    REAL,    DIMENSION(ims:ime,jms:jme),         INTENT(INOUT), OPTIONAL :: PEXPXY       ! factor for river conductance
    REAL,    DIMENSION(ims:ime,1:NSOIL,jms:jme), INTENT(INOUT), OPTIONAL :: smoiseq      ! equilibrium soil moisture content [m3m-3]
    ! output only
    INTEGER, DIMENSION(ims:ime,jms:jme), INTENT(OUT)           :: cropcat             ! crop type
    INTEGER,                             INTENT(OUT), OPTIONAL :: STEPWTD
    ! local
    integer :: itf, jtf, I, J
! ----------------------------------------------------------------------------------

    ! initialize NoahmpIO dimension and key config variables
    NoahmpIO%xstart           = ims
    NoahmpIO%xend             = ime
    NoahmpIO%ystart           = jms
    NoahmpIO%yend             = jme
    NoahmpIO%ids              = ids
    NoahmpIO%ide              = ide
    NoahmpIO%jds              = jds
    NoahmpIO%jde              = jde
    NoahmpIO%kds              = kds
    NoahmpIO%kde              = kde
    NoahmpIO%ims              = ims
    NoahmpIO%ime              = ime
    NoahmpIO%jms              = jms
    NoahmpIO%jme              = jme
    NoahmpIO%kms              = kms
    NoahmpIO%kme              = kme
    NoahmpIO%its              = its
    NoahmpIO%ite              = ite
    NoahmpIO%jts              = jts
    NoahmpIO%jte              = jte
    NoahmpIO%kts              = kts
    NoahmpIO%kte              = kte
    NoahmpIO%NSOIL            = NSOIL
    NoahmpIO%LLANDUSE         = MMINLU
    NoahmpIO%IOPT_CROP        = IOPT_CROP
    NoahmpIO%IOPT_IRR         = IOPT_IRR
    NoahmpIO%IOPT_IRRM        = IOPT_IRRM
    NoahmpIO%SF_URBAN_PHYSICS = SF_URBAN_PHYSICS
    NoahmpIO%IOPT_SOIL        = IOPT_SOIL
    NoahmpIO%IOPT_ALB         = IOPT_ALB
    NoahmpIO%IOPT_WETLAND     = IOPT_WETLAND
    NoahmpIO%IOPT_RUNSUB      = IOPT_RUNSUB
    if ( NoahmpIO%IOPT_ALB == 3 ) then
       NoahmpIO%SNICAR_SNOWOPTICS_OPT = SNICAR_SNOWOPTICS_OPT
       NoahmpIO%SNICAR_DUSTOPTICS_OPT = SNICAR_DUSTOPTICS_OPT
       NoahmpIO%SNICAR_SOLARSPEC_OPT  = SNICAR_SOLARSPEC_OPT
       NoahmpIO%SNICAR_BANDNUMBER_OPT = SNICAR_BANDNUMBER_OPT
    endif

    ! initialze all NoahmpIO variables with default values
    call NoahmpIOVarInitDefault(NoahmpIO)

    ! read in Noahmp table parameters
    call NoahmpReadTable(NoahmpIO)

    ! read in SNICAR parameter netcdif file
    if ( NoahmpIO%IOPT_ALB == 3 ) then
       NoahmpIO%snicar_optic_flnm = "snicar_optics_5bnd_c013122.nc" 
       NoahmpIO%snicar_age_flnm   = "snicar_drdt_bst_fit_60_c070416.nc"
       call SnowInputSnicar(NoahmpIO)
    endif

    !--------- WRF variables mapped to NoahmpIO variables

    ! non-2D variables
    itf = min0(ite, ide-1)
    jtf = min0(jte, jde-1)
    NoahmpIO%DZS                = DZS
    NoahmpIO%FNDSNOWH           = FNDSNOWH
    NoahmpIO%restart_flag       = restart
    NoahmpIO%DTBL               = DT
    NoahmpIO%WTDDT              = WTDDT
    NoahmpIO%DX                 = DX
    NoahmpIO%DY                 = DY

    ! 2D/3D variables
    do J = jts, jtf
    do I = its, itf
    
    ! input only
    NoahmpIO%IVGTYP(I,J)               = IVGTYP(I,J)
    NoahmpIO%ISLTYP(I,J)               = ISLTYP(I,J)
    NoahmpIO%XLAT(I,J)                 = XLAT(I,J)
    NoahmpIO%TSK(I,J)                  = TSK(I,J)
    NoahmpIO%XICE(I,J)                 = XICE(I,J)
    NoahmpIO%CROPTYPE(I,:,J)           = CROPTYPE(I,:,J)
    NoahmpIO%FDEPTHXY(I,J)             = FDEPTHXY(I,J)
    NoahmpIO%MSFTX(I,J)                = MSFTX(I,J)
    NoahmpIO%MSFTY(I,J)                = MSFTY(I,J)
    NoahmpIO%TERRAIN(I,J)              = HT(I,J)
    NoahmpIO%RECHCLIM(I,J)             = RECHCLIM(I,J)
    ! in/out variables
    NoahmpIO%SMOIS(I,:,J)              = SMOIS(I,:,J)
    NoahmpIO%SH2O(I,:,J)               = SH2O(I,:,J)
    NoahmpIO%TSLB(I,:,J)               = TSLB(I,:,J)
    NoahmpIO%SNOW(I,J)                 = SNOW(I,J)    
    NoahmpIO%SNOWH(I,J)                = SNOWH(I,J)   
    NoahmpIO%CANWAT(I,J)               = CANWAT(I,J)  
    NoahmpIO%CANICEXY(I,J)             = CANICEXY(I,J)
    NoahmpIO%CANLIQXY(I,J)             = CANLIQXY(I,J)
    NoahmpIO%TMN(I,J)                  = TMN(I,J)
    NoahmpIO%ISNOWXY(I,J)              = ISNOWXY(I,J)
    NoahmpIO%ZSNSOXY(I,:,J)            = ZSNSOXY(I,:,J) 
    NoahmpIO%TSNOXY(I,:,J)             = TSNOXY(I,:,J)
    NoahmpIO%SNICEXY(I,:,J)            = SNICEXY(I,:,J)
    NoahmpIO%SNLIQXY(I,:,J)            = SNLIQXY(I,:,J)
    NoahmpIO%TVXY(I,J)                 = TVXY(I,J)
    NoahmpIO%TGXY(I,J)                 = TGXY(I,J)
    NoahmpIO%EAHXY(I,J)                = EAHXY(I,J)
    NoahmpIO%TAHXY(I,J)                = TAHXY(I,J)
    NoahmpIO%CMXY(I,J)                 = CMXY(I,J)
    NoahmpIO%CHXY(I,J)                 = CHXY(I,J)
    NoahmpIO%FWETXY(I,J)               = FWETXY(I,J)
    NoahmpIO%SNEQVOXY(I,J)             = SNEQVOXY(I,J)
    NoahmpIO%ALBOLDXY(I,J)             = ALBOLDXY(I,J)
    NoahmpIO%QSNOWXY(I,J)              = QSNOWXY(I,J)
    NoahmpIO%QRAINXY(I,J)              = QRAINXY(I,J)
    NoahmpIO%WSLAKEXY(I,J)             = WSLAKEXY(I,J)
    NoahmpIO%ZWTXY(I,J)                = ZWTXY(I,J)
    NoahmpIO%WAXY(I,J)                 = WAXY(I,J)
    NoahmpIO%WTXY(I,J)                 = WTXY(I,J)
    NoahmpIO%LFMASSXY(I,J)             = LFMASSXY(I,J)
    NoahmpIO%RTMASSXY(I,J)             = RTMASSXY(I,J)
    NoahmpIO%STMASSXY(I,J)             = STMASSXY(I,J)
    NoahmpIO%WOODXY(I,J)               = WOODXY(I,J)
    NoahmpIO%GRAINXY(I,J)              = GRAINXY(I,J)
    NoahmpIO%GDDXY(I,J)                = GDDXY(I,J)
    NoahmpIO%STBLCPXY(I,J)             = STBLCPXY(I,J)
    NoahmpIO%FASTCPXY(I,J)             = FASTCPXY(I,J)
    NoahmpIO%LAI(I,J)                  = LAI(I,J)
    NoahmpIO%XSAIXY(I,J)               = XSAIXY(I,J)
    NoahmpIO%QTDRAIN(I,J)              = QTDRAIN(I,J)
    NoahmpIO%IRNUMSI(I,J)              = IRNUMSI(I,J)
    NoahmpIO%IRNUMMI(I,J)              = IRNUMMI(I,J)
    NoahmpIO%IRNUMFI(I,J)              = IRNUMFI(I,J)
    NoahmpIO%IRWATSI(I,J)              = IRWATSI(I,J)
    NoahmpIO%IRWATMI(I,J)              = IRWATMI(I,J)
    NoahmpIO%IRWATFI(I,J)              = IRWATFI(I,J)
    NoahmpIO%IRELOSS(I,J)              = IRELOSS(I,J)
    NoahmpIO%IRSIVOL(I,J)              = IRSIVOL(I,J)
    NoahmpIO%IRMIVOL(I,J)              = IRMIVOL(I,J)
    NoahmpIO%IRFIVOL(I,J)              = IRFIVOL(I,J)
    NoahmpIO%IRRSPLH(I,J)              = IRRSPLH(I,J)
    NoahmpIO%T2MVXY(I,J)               = T2MVXY(I,J)
    NoahmpIO%T2MBXY(I,J)               = T2MBXY(I,J)
    NoahmpIO%SMCWTDXY(I,J)             = SMCWTDXY(I,J)
    NoahmpIO%DEEPRECHXY(I,J)           = DEEPRECHXY(I,J)
    NoahmpIO%RECHXY(I,J)               = RECHXY(I,J)
    NoahmpIO%QRFSXY(I,J)               = QRFSXY(I,J)
    NoahmpIO%QSPRINGSXY(I,J)           = QSPRINGSXY(I,J)
    NoahmpIO%QSLATXY(I,J)              = QSLATXY(I,J)
    NoahmpIO%AREAXY(I,J)               = AREAXY(I,J)
    NoahmpIO%RIVERBEDXY(I,J)           = RIVERBEDXY(I,J)
    NoahmpIO%EQZWT(I,J)                = EQZWT(I,J)
    NoahmpIO%RIVERCONDXY(I,J)          = RIVERCONDXY(I,J)
    NoahmpIO%PEXPXY(I,J)               = PEXPXY(I,J)
    NoahmpIO%SMOISEQ(I,:,J)            = SMOISEQ(I,:,J)
    NoahmpIO%ALBSOILDIRXY(I,:,J)       = ALBSOILDIRXY(I,:,J)
    NoahmpIO%ALBSOILDIFXY(I,:,J)       = ALBSOILDIFXY(I,:,J)
    if ( NoahmpIO%IOPT_WETLAND > 0 ) then
       NoahmpIO%FSATXY(I,J)            = FSATXY(I,J)
       NoahmpIO%WSURFXY(I,J)           = WSURFXY(I,J)
    endif
    if ( NoahmpIO%IOPT_ALB == 3 ) then
       NoahmpIO%SNRDSXY(I,:,J)         = SNRDSXY(I,:,J)
       NoahmpIO%SNFRXY(I,:,J)          = SNFRXY(I,:,J)
       NoahmpIO%BCPHIXY(I,:,J)         = BCPHIXY(I,:,J)
       NoahmpIO%BCPHOXY(I,:,J)         = BCPHOXY(I,:,J)
       NoahmpIO%OCPHIXY(I,:,J)         = OCPHIXY(I,:,J)
       NoahmpIO%OCPHOXY(I,:,J)         = OCPHOXY(I,:,J)
       NoahmpIO%DUST1XY(I,:,J)         = DUST1XY(I,:,J)
       NoahmpIO%DUST2XY(I,:,J)         = DUST2XY(I,:,J)
       NoahmpIO%DUST3XY(I,:,J)         = DUST3XY(I,:,J)
       NoahmpIO%DUST4XY(I,:,J)         = DUST4XY(I,:,J)
       NoahmpIO%DUST5XY(I,:,J)         = DUST5XY(I,:,J)
       NoahmpIO%MassConcBCPHIXY(I,:,J) = MassConcBCPHIXY(I,:,J)
       NoahmpIO%MassConcBCPHOXY(I,:,J) = MassConcBCPHOXY(I,:,J)
       NoahmpIO%MassConcOCPHIXY(I,:,J) = MassConcOCPHIXY(I,:,J)
       NoahmpIO%MassConcOCPHOXY(I,:,J) = MassConcOCPHOXY(I,:,J)
       NoahmpIO%MassConcDUST1XY(I,:,J) = MassConcDUST1XY(I,:,J)
       NoahmpIO%MassConcDUST2XY(I,:,J) = MassConcDUST2XY(I,:,J)
       NoahmpIO%MassConcDUST3XY(I,:,J) = MassConcDUST3XY(I,:,J)
       NoahmpIO%MassConcDUST4XY(I,:,J) = MassConcDUST4XY(I,:,J)
       NoahmpIO%MassConcDUST5XY(I,:,J) = MassConcDUST5XY(I,:,J)
    endif

    enddo ! I
    enddo ! J

    !--------- WRF -> NoahmpIO variables mapping ends

    !--------- main Noahmp initialization module
    call NoahmpInitMain(NoahmpIO)
    !---------

    !--------- initialized NoahmpIO variable mapped to WRF variables
    do J = jts, jtf
    do I = its, itf

    ! in/out variables
    SMOIS(I,:,J)        = NoahmpIO%SMOIS(I,:,J) 
    SH2O(I,:,J)         = NoahmpIO%SH2O(I,:,J)
    TSLB(I,:,J)         = NoahmpIO%TSLB(I,:,J)
    SNOW(I,J)           = NoahmpIO%SNOW(I,J)
    SNOWH(I,J)          = NoahmpIO%SNOWH(I,J) 
    CANWAT(I,J)         = NoahmpIO%CANWAT(I,J)
    CANICEXY(I,J)       = NoahmpIO%CANICEXY(I,J)
    CANLIQXY(I,J)       = NoahmpIO%CANLIQXY(I,J)
    TMN(I,J)            = NoahmpIO%TMN(I,J)
    ISNOWXY(I,J)        = NoahmpIO%ISNOWXY(I,J)
    ZSNSOXY(I,:,J)      = NoahmpIO%ZSNSOXY(I,:,J)
    TSNOXY(I,:,J)       = NoahmpIO%TSNOXY(I,:,J)
    SNICEXY(I,:,J)      = NoahmpIO%SNICEXY(I,:,J)
    SNLIQXY(I,:,J)      = NoahmpIO%SNLIQXY(I,:,J)
    TVXY(I,J)           = NoahmpIO%TVXY(I,J)
    TGXY(I,J)           = NoahmpIO%TGXY(I,J)
    EAHXY(I,J)          = NoahmpIO%EAHXY(I,J)
    TAHXY(I,J)          = NoahmpIO%TAHXY(I,J)
    CMXY(I,J)           = NoahmpIO%CMXY(I,J)
    CHXY(I,J)           = NoahmpIO%CHXY(I,J)
    FWETXY(I,J)         = NoahmpIO%FWETXY(I,J)
    SNEQVOXY(I,J)       = NoahmpIO%SNEQVOXY(I,J)
    ALBOLDXY(I,J)       = NoahmpIO%ALBOLDXY(I,J)
    QSNOWXY(I,J)        = NoahmpIO%QSNOWXY(I,J)
    QRAINXY(I,J)        = NoahmpIO%QRAINXY(I,J)
    WSLAKEXY(I,J)       = NoahmpIO%WSLAKEXY(I,J)
    ZWTXY(I,J)          = NoahmpIO%ZWTXY(I,J)
    WAXY(I,J)           = NoahmpIO%WAXY(I,J)
    WTXY(I,J)           = NoahmpIO%WTXY(I,J)
    LFMASSXY(I,J)       = NoahmpIO%LFMASSXY(I,J)
    RTMASSXY(I,J)       = NoahmpIO%RTMASSXY(I,J)
    STMASSXY(I,J)       = NoahmpIO%STMASSXY(I,J)
    WOODXY(I,J)         = NoahmpIO%WOODXY(I,J)
    GRAINXY(I,J)        = NoahmpIO%GRAINXY(I,J)
    GDDXY(I,J)          = NoahmpIO%GDDXY(I,J)
    STBLCPXY(I,J)       = NoahmpIO%STBLCPXY(I,J)
    FASTCPXY(I,J)       = NoahmpIO%FASTCPXY(I,J)
    LAI(I,J)            = NoahmpIO%LAI(I,J)
    XSAIXY(I,J)         = NoahmpIO%XSAIXY(I,J)
    QTDRAIN(I,J)        = NoahmpIO%QTDRAIN(I,J)
    IRNUMSI(I,J)        = NoahmpIO%IRNUMSI(I,J)
    IRNUMMI(I,J)        = NoahmpIO%IRNUMMI(I,J)
    IRNUMFI(I,J)        = NoahmpIO%IRNUMFI(I,J)
    IRWATSI(I,J)        = NoahmpIO%IRWATSI(I,J)
    IRWATMI(I,J)        = NoahmpIO%IRWATMI(I,J)
    IRWATFI(I,J)        = NoahmpIO%IRWATFI(I,J)
    IRELOSS(I,J)        = NoahmpIO%IRELOSS(I,J)
    IRSIVOL(I,J)        = NoahmpIO%IRSIVOL(I,J)
    IRMIVOL(I,J)        = NoahmpIO%IRMIVOL(I,J)
    IRFIVOL(I,J)        = NoahmpIO%IRFIVOL(I,J)
    IRRSPLH(I,J)        = NoahmpIO%IRRSPLH(I,J)
    T2MVXY(I,J)         = NoahmpIO%T2MVXY(I,J)
    T2MBXY(I,J)         = NoahmpIO%T2MBXY(I,J)
    SMCWTDXY(I,J)       = NoahmpIO%SMCWTDXY(I,J)
    DEEPRECHXY(I,J)     = NoahmpIO%DEEPRECHXY(I,J)
    RECHXY(I,J)         = NoahmpIO%RECHXY(I,J)
    QRFSXY(I,J)         = NoahmpIO%QRFSXY(I,J)
    QSPRINGSXY(I,J)     = NoahmpIO%QSPRINGSXY(I,J)
    QSLATXY(I,J)        = NoahmpIO%QSLATXY(I,J)
    AREAXY(I,J)         = NoahmpIO%AREAXY(I,J)
    RIVERBEDXY(I,J)     = NoahmpIO%RIVERBEDXY(I,J)
    EQZWT(I,J)          = NoahmpIO%EQZWT(I,J)
    RIVERCONDXY(I,J)    = NoahmpIO%RIVERCONDXY(I,J)
    PEXPXY(I,J)         = NoahmpIO%PEXPXY(I,J)
    SMOISEQ(I,:,J)      = NoahmpIO%SMOISEQ(I,:,J)
    CHSTARXY(I,J)       = 0.1 ! dummy
    ALBSOILDIRXY(I,:,J) = NoahmpIO%ALBSOILDIRXY(I,:,J)
    ALBSOILDIFXY(I,:,J) = NoahmpIO%ALBSOILDIFXY(I,:,J)
    if ( NoahmpIO%IOPT_WETLAND > 0 ) then
       FSATXY(I,J)      = NoahmpIO%FSATXY(I,J)
       WSURFXY(I,J)     = NoahmpIO%WSURFXY(I,J)
    endif
    if ( NoahmpIO%IOPT_ALB == 3 ) then
       SNRDSXY(I,:,J)         = NoahmpIO%SNRDSXY(I,:,J)
       SNFRXY(I,:,J)          = NoahmpIO%SNFRXY(I,:,J)
       BCPHIXY(I,:,J)         = NoahmpIO%BCPHIXY(I,:,J)
       BCPHOXY(I,:,J)         = NoahmpIO%BCPHOXY(I,:,J)
       OCPHIXY(I,:,J)         = NoahmpIO%OCPHIXY(I,:,J)
       OCPHOXY(I,:,J)         = NoahmpIO%OCPHOXY(I,:,J)
       DUST1XY(I,:,J)         = NoahmpIO%DUST1XY(I,:,J)
       DUST2XY(I,:,J)         = NoahmpIO%DUST2XY(I,:,J)
       DUST3XY(I,:,J)         = NoahmpIO%DUST3XY(I,:,J)
       DUST4XY(I,:,J)         = NoahmpIO%DUST4XY(I,:,J)
       DUST5XY(I,:,J)         = NoahmpIO%DUST5XY(I,:,J)
       MassConcBCPHIXY(I,:,J) = NoahmpIO%MassConcBCPHIXY(I,:,J)
       MassConcBCPHOXY(I,:,J) = NoahmpIO%MassConcBCPHOXY(I,:,J)
       MassConcOCPHIXY(I,:,J) = NoahmpIO%MassConcOCPHIXY(I,:,J)
       MassConcOCPHOXY(I,:,J) = NoahmpIO%MassConcOCPHOXY(I,:,J)
       MassConcDUST1XY(I,:,J) = NoahmpIO%MassConcDUST1XY(I,:,J)
       MassConcDUST2XY(I,:,J) = NoahmpIO%MassConcDUST2XY(I,:,J)
       MassConcDUST3XY(I,:,J) = NoahmpIO%MassConcDUST3XY(I,:,J)
       MassConcDUST4XY(I,:,J) = NoahmpIO%MassConcDUST4XY(I,:,J)
       MassConcDUST5XY(I,:,J) = NoahmpIO%MassConcDUST5XY(I,:,J)
    endif

    ! out variables only
    CROPCAT(I,J) = NoahmpIO%CROPCAT(I,J)

    enddo ! I
    enddo ! J

    STEPWTD      = NoahmpIO%STEPWTD

    !--------- NoahmpIO -> WRF variables mapping ends

  end subroutine NoahmpWRFinit

end module NoahmpWRFinitMod
