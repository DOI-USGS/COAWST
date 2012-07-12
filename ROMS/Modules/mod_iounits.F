#include "cppdefs.h"
      MODULE mod_iounits
!
!svn $Id: mod_iounits.F 755 2008-09-14 19:07:08Z jcwarner $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2010 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  ADJbase     Output adjoint base file name.                          !
!  ADJname     Output adjoint file name.                               !
!  ADSname     Input adjoint sensitivity functional file.              !
!  AVGbase     Output averages base file name.                         !
!  AVGname     Output averages file name.                              !
!  Aname       Atmosphere model stadard input file name.               !
!  BLKname     Input bulk fluxes file name.                            !
!  BRYname     Input boundary data file name.                          !
!  IWBNDname     Input boundary data file name for InWave model        !
!  IWSWNname     Input spectral SWAN data file name for InWave model   !
!  CLMname     Input climatology file name.                            !
!  DIAbase     Output diagnostics base file name.                      !
!  DIAname     Output diagnostics file name.                           !
!  ERRname     Output 4DVar posterior error file name.                 !
!  FLTname     Output floats file name.                                !
!  FRCids      NetCDF ID descriptors associated with FRCname.          !
!  FRCname     Input forcing fields file name(s).                      !
!  FWDbase     Input forward solution base file name.                  !
!  FWDname     Input forward solution file name.                       !
!  GRDname     Input grid file name.                                   !
!  GSTname     Input/output GST analysis check pointing NetCDF file.   !
!  HISbase     Output history base file name.                          !
!  HISname     Output history file name.                               !
!  HSSbase     Input/output Hessian eigenvectors base file name.       !
!  HSSname     Input/output Hessian eigenvectors file name.            !
!  IADname     Input adjoint initial conditions file name.             !
!  INIname     Input nonlinear initial conditions file name.           !
!  IPRname     Input representer initial conditions file name.         !
!  ITLname     Input tangent linear initial conditions file name.      !
!  Iname       Physical parameters standard input script file name.    !
!  LCZname     Input/output Lanczos vectors file name.                 !
!  MODname     Output 4DVAR processed fields file name.                !
!  MyAppCPP    Application C-preprocessing flag.                       !
!  NRMname     Input/output error covariance normalization file name:  !
!                NRMname(1)  initial conditions                        !
!                NRMname(2)  model                                     !
!                NRMname(3)  boundary conditions                       !
!                NRMname(4)  surface forcing                           !
!  OBSname     Input/output datum observations file name.              !
!  REPname     Input/output representer coefficients file name.        !
!  Rerror      Running error messages.                                 !
!  RSTname     Output restart file name.                               !
!  SSHname     Input SSH observations file name.                       !
!  SSTname     Input SST observations file name.                       !
!  TIDEname    Input tide forcing file name.                           !
!  TLFname     Input/output tangent linear impulse forcing file name.  !
!  TLMbase     Output tangent linear base file name.                   !
!  TLMname     Output tangent linear file name.                        !
!  TOBSname    Input tracer observations file name.                    !
!  USRname     USER input/output generic file name.                    !
!  VSURname    Input surface currents observations file name.          !
!  VOBSname    Input horizontal currents observations file name.       !
!  Wname       Wave model stadard input file name.                     !
!  aparnam     Input assimilation parameters file name.                !
!  bparnam     Input biology parameters file name.                     !
!  fposnam     Input initial floats positions file name.               !
!  ioerror     IO error flag.                                          !
!  iparnam     Input ice parameters file name.                         !
!  ncfile      Current NetCDF file name being processed.               !
!  nFfiles     Number of forcing files.                                !
!  sparnam     Input sediment transport parameters file name.          !
!  sposnam     Input station positions file name.                      !
!  SourceFile  Current executed file name. It is used for IO error     !
!                purposes.                                             !
!  STAname     Output station data file name.                          !
!  STDname     Input error covariance standard deviations file name:   !
!                STDname(1)  initial conditions                        !
!                STDname(2)  model                                     !
!                STDname(3)  boundary conditions                       !
!                STDname(4)  surface forcing                           !
!  stdinp      Unit number for standard input (often 5).               !
!  stdout      Unit number for standard output (often 6).              !
!  usrout      Unit number for generic USER output.                    !
!  varname     Input IO variables information file name.               !
!                                                                      !
!=======================================================================
!
        USE mod_param

        implicit none

        integer, parameter :: stdinp = 5
        integer, parameter :: stdout = 6
        integer, parameter :: usrout = 10
        integer :: ioerror = 0

        integer, dimension(Ngrids) :: nFfiles

        integer, allocatable :: FRCids(:,:)

        character (len=80) :: SourceFile

        character (len=50), dimension(8) :: Rerror =                    &
     &       (/ ' ROMS/TOMS - Blows up ................ exit_flag: ',   &
     &          ' ROMS/TOMS - Input error ............. exit_flag: ',   &
     &          ' ROMS/TOMS - Output error ............ exit_flag: ',   &
     &          ' ROMS/TOMS - I/O error ............... exit_flag: ',   &
     &          ' ROMS/TOMS - Configuration error ..... exit_flag: ',   &
     &          ' ROMS/TOMS - Partition error ......... exit_flag: ',   &
     &          ' ROMS/TOMS - Illegal input parameter . exit_flag: ',   &
     &          ' ROMS/TOMS - Fatal algorithm result .. exit_flag: ' /)

        character (len=80), allocatable :: FRCname(:,:)

        character (len=80), dimension(Ngrids) :: ADJbase
        character (len=80), dimension(Ngrids) :: ADJname
        character (len=80), dimension(Ngrids) :: ADSname
        character (len=80), dimension(Ngrids) :: AVGbase
        character (len=80), dimension(Ngrids) :: AVGname
        character (len=80), dimension(Ngrids) :: BLKname
        character (len=80), dimension(Ngrids) :: BRYname
        character (len=80), dimension(Ngrids) :: CLMname
        character (len=80), dimension(Ngrids) :: DIAbase
        character (len=80), dimension(Ngrids) :: DIAname
        character (len=80), dimension(Ngrids) :: ERRname
        character (len=80), dimension(Ngrids) :: FLTname
        character (len=80), dimension(Ngrids) :: FWDbase
        character (len=80), dimension(Ngrids) :: FWDname
        character (len=80), dimension(Ngrids) :: GRDname
        character (len=80), dimension(Ngrids) :: GSTname
        character (len=80), dimension(Ngrids) :: HISbase
        character (len=80), dimension(Ngrids) :: HISname
        character (len=80), dimension(Ngrids) :: HSSbase
        character (len=80), dimension(Ngrids) :: HSSname
        character (len=80), dimension(Ngrids) :: IADname
        character (len=80), dimension(Ngrids) :: INIname
        character (len=80), dimension(Ngrids) :: IRPname
        character (len=80), dimension(Ngrids) :: ITLname
        character (len=80), dimension(Ngrids) :: LCZname
        character (len=80), dimension(Ngrids) :: MODname
        character (len=80), dimension(Ngrids) :: OBSname
        character (len=80), dimension(Ngrids) :: REPname
#ifdef MODEL_COUPLING
        character (len=80), dimension(Ngrids) :: O2Aname
        character (len=80), dimension(Ngrids) :: O2Wname
        character (len=80), dimension(Ngrids) :: W2Aname
        character (len=80), dimension(Ngrids) :: W2Oname
        character (len=80), dimension(Ngrids) :: A2Wname
        character (len=80), dimension(Ngrids) :: A2Oname
        character (len=80) :: Aname
# ifdef INWAVE_MODEL
        character (len=80), dimension(Ngrids) :: IWINIname
        character (len=80), dimension(Ngrids) :: IWBNDname
# endif
#endif
        character (len=80), dimension(Ngrids) :: SSHname
        character (len=80), dimension(Ngrids) :: SSTname
        character (len=80), dimension(Ngrids) :: TIDEname
        character (len=80), dimension(Ngrids) :: TLFname
        character (len=80), dimension(Ngrids) :: TLMbase
        character (len=80), dimension(Ngrids) :: TLMname
        character (len=80), dimension(Ngrids) :: TOBSname
        character (len=80), dimension(Ngrids) :: VSURname
        character (len=80), dimension(Ngrids) :: VOBSname
        character (len=80), dimension(Ngrids) :: RSTname
        character (len=80), dimension(Ngrids) :: STAname

        character (len=80), dimension(4,Ngrids) :: NRMname
        character (len=80), dimension(4,Ngrids) :: STDname

#ifdef INWAVE_MODEL
        character (len=80), dimension(Ngrids) :: IWINIname
        character (len=80), dimension(Ngrids) :: IWBNDname
        character (len=80), dimension(Ngrids) :: IWSWNname
#endif
        character (len=80) :: Iname
        character (len=80) :: MyAppCPP
        character (len=80) :: USRname
        character (len=80) :: aparnam
        character (len=80) :: bparnam
        character (len=80) :: fposnam
        character (len=80) :: iparnam
        character (len=80) :: ncfile
        character (len=80) :: sparnam
        character (len=80) :: sposnam
        character (len=80) :: varname

      END MODULE mod_iounits
