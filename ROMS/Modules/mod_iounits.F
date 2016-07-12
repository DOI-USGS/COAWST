#include "cppdefs.h"
      MODULE mod_iounits
!
!svn $Id: mod_iounits.F 807 2016-07-09 02:03:55Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  The input and output files information is stored in compact derived !
!  type structure, TYPE(T_IO):                                         !
!                                                                      !
!    ADM       Output adjoint model history data                       !
!    ADS       Input adjoint sensitivity functional                    !
!    AVG       Output time-averaged data                               !
!    AVG2      Secondary time-averaged data                            !
!    BLK       Input bulk fluxes data in adjoint-based applications    !
!    BRY       Input open boundary conditions data                     !
!    CLM       Input climatology data                                  !
!    DAI       Output data assimilation next cycle initial conditions  !
!                (4D-Var analysis) or restart (Ensemble Kalman Filter) !
!    DAV       Output 4D-Var data assimilation variables               !
!    DIA       Output diagnostics fields                               !
!    ERR       Output 4DVar posterior error estimate                   !
!    FLT       Output Lagrangian trajectories data                     !
!    FRC       Input forcing data                                      !
!    FWD       Input basic state forward solution                      !
!    GRD       Input grid data                                         !
!    GST       Input/output GST analysis check pointing data           !
!    HIS       Output nonlinear model history data                     !
!    HSS       Input/output Hessian eigenvectors                       !
!    IAD       Input adjoint model initial conditions                  !
!    INI       Input nonlinear model initial conditions                !
!    IPR       Input representer model initial conditions              !
!    ITL       Input tangent linear model initial conditions           !
!    LCZ       Input/output Lanczos vectors                            !
!    LZE       Input/output time-evolved Lanczos vectors               !
!    NRM       Input/output error covariance normalization data        !
!                NRM(1,ng)  initial conditions                         !
!                NRM(2,ng)  model error                                !
!                NRM(3,ng)  lateral open boundary conditions           !
!                NRM(4,ng)  surface forcing                            !
!    NUD       Input climatology nudging coefficients                  !
!    OBS       Input/output 4D-Var observations                        !
!    RST       Output restart data                                     !
!    TIDE      Input tide forcing                                      !
!    TLF       Input/output tangent linear model impulse forcing       !
!    TLM       Output tangent linear model history                     !
!    SSF       Input Sources/Sinks forcing (river runoff)              !
!    STA       Output station data                                     !
!    STD       Input error covariance standard deviations              !
!                STD(1,ng)  initial conditions                         !
!                STD(2,ng)  model error                                !
!                STD(3,ng)  lateral open boundary conditions           !
!                STD(4,ng)  surface forcing                            !
!                                                                      !
!  Standard input files:                                               !
!                                                                      !
!  Iname       Physical parameters standard input script file name.    !
!  NGCname     Nested grids contact point information file name.       !
!  USRname     USER input/output generic file name.                    !
!  Wname       Wave model stadard input file name.                     !
!  aparnam     Input assimilation parameters file name.                !
!  bparnam     Input biology parameters file name.                     !
!  fbionam     Input floats biological behavior parameters file name.  !
!  fposnam     Input initial floats positions file name.               !
!  iparnam     Input ice parameters file name.                         !
!  sparnam     Input sediment transport parameters file name.          !
!  vegnam      Input vegetation input parameters file name.             !
!  sposnam     Input station positions file name.                      !
!  varname     Input IO variables information file name.               !
!                                                                      !
!  stdinp      Unit number for standard input (often 5).               !
!  stdout      Unit number for standard output (often 6).              !
!  usrout      Unit number for generic USER output.                    !
!                                                                      !
!  Miscellaneous variables:                                            !
!                                                                      !
!  BRYids      NetCDF file ID associated with each boundary field.     !
!  CLMids      NetCDF file ID associated with each climatology field.  !
!  FRCids      NetCDF file ID associated with each forcing field.      !
!  Rerror      Running error messages.                                 !
!  SourceFile  Current executed file name. It is used for IO error     !
!                purposes.                                             !
!  ioerror     IO error flag.                                          !
!  ncfile      Current NetCDF file name being processed.               !
!  nFfiles     Number of forcing files.                                !
!                                                                      !
!=======================================================================
!
      USE mod_param
!
      implicit none
!
!  I/O units.
!
      integer, parameter :: stdinp = 5      ! standard input
      integer, parameter :: stdout = 6      ! standard output
      integer, parameter :: usrout = 10     ! generic user unit
!
!  I/O files management, derived type structures.
!
      TYPE(T_IO), allocatable :: ADM(:)     ! ADM history fields
      TYPE(T_IO), allocatable :: ADS(:)     ! sensitivity functional
      TYPE(T_IO), allocatable :: AVG(:)     ! time-averaged fields
      TYPE(T_IO), allocatable :: AVG2(:)    ! secondary averaged fields
      TYPE(T_IO), allocatable :: BLK(:)     ! bulk fluxes fields
      TYPE(T_IO), allocatable :: DAI(:)     ! Data assilation IC/restart
      TYPE(T_IO), allocatable :: DAV(:)     ! 4D-Var variables
      TYPE(T_IO), allocatable :: DIA(:)     ! diagnostics fields
      TYPE(T_IO), allocatable :: ERR(:)     ! 4D-Var posterior error
      TYPE(T_IO), allocatable :: FIL(:,:)   ! Filter restart
      TYPE(T_IO), allocatable :: FLT(:)     ! Lagrangian trajectories
      TYPE(T_IO), allocatable :: FWD(:)     ! forward solution
      TYPE(T_IO), allocatable :: GRD(:)     ! grid data
      TYPE(T_IO), allocatable :: GST(:)     ! generalized stability
      TYPE(T_IO), allocatable :: HIS(:)     ! NLM history fields
      TYPE(T_IO), allocatable :: HIS2(:)    ! surface history fields
      TYPE(T_IO), allocatable :: HSS(:)     ! Hessian eigenvectors
      TYPE(T_IO), allocatable :: IAD(:)     ! ADM initial conditions
      TYPE(T_IO), allocatable :: INI(:)     ! NLM initial conditions
      TYPE(T_IO), allocatable :: IRP(:)     ! RPM initial conditions
      TYPE(T_IO), allocatable :: ITL(:)     ! TLM initial conditions
      TYPE(T_IO), allocatable :: LCZ(:)     ! Lanczos vectors
      TYPE(T_IO), allocatable :: LZE(:)     ! evolved Lanczos vectors
      TYPE(T_IO), allocatable :: NRM(:,:)   ! normalization
      TYPE(T_IO), allocatable :: NUD(:)     ! nudging coefficients
      TYPE(T_IO), allocatable :: OBS(:)     ! observations
      TYPE(T_IO), allocatable :: RST(:)     ! restart fields
      TYPE(T_IO), allocatable :: SSF(:)     ! Sources/Sinks forcing
      TYPE(T_IO), allocatable :: STA(:)     ! stations data
      TYPE(T_IO), allocatable :: STD(:,:)   ! standard deviation
      TYPE(T_IO), allocatable :: TIDE(:)    ! tidal forcing
      TYPE(T_IO), allocatable :: TLF(:)     ! TLM impulse fields
      TYPE(T_IO), allocatable :: TLM(:)     ! TLM history fields
!
!  Input boundary condition data.
!
      integer, allocatable :: nBCfiles(:)
      integer, allocatable :: BRYids(:,:)

      TYPE(T_IO), allocatable :: BRY(:,:)
!
!  Input climatology data.
!
      integer, allocatable :: nCLMfiles(:)
      integer, allocatable :: CLMids(:,:)

      TYPE(T_IO), allocatable :: CLM(:,:)   ! climatology fields
!
!  Input forcing data.
!
      integer, allocatable :: nFfiles(:)
      integer, allocatable :: FRCids(:,:)

      TYPE(T_IO), allocatable :: FRC(:,:)
!
!  Error messages.
!
      character (len=50), dimension(9) :: Rerror =                      &
     &       (/ ' ROMS/TOMS - Blows up ................ exit_flag: ',   &
     &          ' ROMS/TOMS - Input error ............. exit_flag: ',   &
     &          ' ROMS/TOMS - Output error ............ exit_flag: ',   &
     &          ' ROMS/TOMS - I/O error ............... exit_flag: ',   &
     &          ' ROMS/TOMS - Configuration error ..... exit_flag: ',   &
     &          ' ROMS/TOMS - Partition error ......... exit_flag: ',   &
     &          ' ROMS/TOMS - Illegal input parameter . exit_flag: ',   &
     &          ' ROMS/TOMS - Fatal algorithm result .. exit_flag: ',   &
     &          ' ROMS/TOMS - Frazil ice problem ...... exit_flag: ' /)
!
!  Standard input scripts file names.
!
      character (len=256) :: Iname          ! ROMS physical parameters
!     character (len=256) :: Wname          ! wave model standard input
      character (len=256) :: NGCname        ! contact points file name
      character (len=256) :: USRname        ! use generic file name
      character (len=256) :: aparnam        ! assimilation parameters
      character (len=256) :: bparnam        ! biology model parameters
      character (len=256) :: fbionam        ! floats behavior parameters
      character (len=256) :: fposnam        ! floats positions
      character (len=256) :: iparnam        ! ice model parameters
      character (len=256) :: sparnam        ! sediment model parameters
      character (len=256) :: vegnam         ! vegetation model parameters
      character (len=256) :: sposnam        ! station positions
      character (len=256) :: varname        ! I/O metadata
!
!  Miscelaneous variables.
!
      integer :: ioerror = 0                ! I/O error flag

      character (len=256) :: MyAppCPP       ! application CPP flag
      character (len=256) :: SourceFile     ! current executed ROMS file
      character (len=256) :: ncfile         ! current NetCDF file
!
      CONTAINS
!
      SUBROUTINE allocate_iounits
!
!=======================================================================
!                                                                      !
!  This routine allocates several variables in the module that depend  !
!  on the number of nested grids.                                      !
!                                                                      !
!=======================================================================
!
!  Local variable declarations.
!
      integer :: i

      character (len=1), parameter :: blank = ' '
!
!-----------------------------------------------------------------------
!  Allocate I/O files management, derived type structures.
!-----------------------------------------------------------------------
!
      allocate ( ADM(Ngrids) )
      allocate ( ADS(Ngrids) )
      allocate ( AVG(Ngrids) )
      allocate ( AVG2(Ngrids) )
      allocate ( BLK(Ngrids) )
      allocate ( DAI(Ngrids) )
      allocate ( DAV(Ngrids) )
      allocate ( DIA(Ngrids) )
      allocate ( ERR(Ngrids) )
#ifdef FILTERED_RST
      allocate ( FIL(nfile,Ngrids) )
#endif
      allocate ( FLT(Ngrids) )
      allocate ( FWD(Ngrids) )
      allocate ( GRD(Ngrids) )
      allocate ( GST(Ngrids) )
      allocate ( HIS(Ngrids) )
      allocate ( HIS2(Ngrids) )
      allocate ( HSS(Ngrids) )
      allocate ( IAD(Ngrids) )
      allocate ( INI(Ngrids) )
      allocate ( IRP(Ngrids) )
      allocate ( ITL(Ngrids) )
      allocate ( LCZ(Ngrids) )
      allocate ( LZE(Ngrids) )
      allocate ( NUD(Ngrids) )
      allocate ( OBS(Ngrids) )
      allocate ( RST(Ngrids) )
      allocate ( SSF(Ngrids) )
      allocate ( STA(Ngrids) )
      allocate ( TIDE(Ngrids) )
      allocate ( TLF(Ngrids) )
      allocate ( TLM(Ngrids) )

      allocate ( NRM(4,Ngrids) )
      allocate ( STD(4,Ngrids) )

!
!-----------------------------------------------------------------------
!  Allocate variables.
!-----------------------------------------------------------------------
!
      allocate ( nBCfiles(Ngrids) )
      allocate ( nCLMfiles(Ngrids) )
      allocate ( nFfiles(Ngrids) )
!
!-----------------------------------------------------------------------
!  Initialize file names to blanks.
!-----------------------------------------------------------------------
!
      DO i=1,LEN(Iname)
        Iname(i:i)=blank
!       Wname(i:i)=blank
        NGCname(i:i)=blank
        USRname(i:i)=blank
        aparnam(i:i)=blank
        bparnam(i:i)=blank
        fbionam(i:i)=blank
        fposnam(i:i)=blank
        iparnam(i:i)=blank
        sparnam(i:i)=blank
        vegnam(i:i)=blank
        sposnam(i:i)=blank
      END DO

      RETURN
      END SUBROUTINE allocate_iounits

      END MODULE mod_iounits
