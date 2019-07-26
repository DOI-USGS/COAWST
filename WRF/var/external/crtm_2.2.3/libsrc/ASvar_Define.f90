! 
! ASvar_Define
!
! Module defining the CRTM AerosolScatter module internal
! variable object.
! 
!
! CREATION HISTORY:
!       Written by:     Paul van Delst, 14-Feb-2012
!                       paul.vandelst@noaa.gov
!                       

MODULE ASvar_Define

  ! -----------------
  ! Environment setup
  ! -----------------
  ! Module use
  USE Type_Kinds           , ONLY: fp
  USE Message_Handler      , ONLY: SUCCESS, FAILURE, INFORMATION, Display_Message
  USE Compare_Float_Numbers, ONLY: OPERATOR(.EqualTo.)
  USE File_Utility         , ONLY: File_Open, File_Exists
  USE Binary_File_Utility  , ONLY: Open_Binary_File      , &
                                   WriteGAtts_Binary_File, &
                                   ReadGAtts_Binary_File
  USE CRTM_Interpolation   , ONLY: NPTS      , &
                                   LPoly_type
  ! Disable implicit typing
  IMPLICIT NONE


  ! ------------
  ! Visibilities
  ! ------------
  ! Everything private by default
  PRIVATE
  ! Datatypes
  PUBLIC :: ASvar_type
  PUBLIC :: ASinterp_type
  ! Operators
  PUBLIC :: OPERATOR(==)
  ! Procedures
  PUBLIC :: ASvar_Associated
  PUBLIC :: ASvar_Destroy
  PUBLIC :: ASvar_Create
  PUBLIC :: ASvar_Inspect
  PUBLIC :: ASvar_ValidRelease
  PUBLIC :: ASvar_Info
  PUBLIC :: ASvar_DefineVersion
  PUBLIC :: ASvar_InquireFile
  PUBLIC :: ASvar_ReadFile
  PUBLIC :: ASvar_WriteFile


  ! ---------------------
  ! Procedure overloading
  ! ---------------------
  INTERFACE OPERATOR(==)
    MODULE PROCEDURE ASvar_Equal
  END INTERFACE OPERATOR(==)


  ! -----------------
  ! Module parameters
  ! -----------------
  CHARACTER(*), PARAMETER :: MODULE_VERSION_ID = &
    '$Id: ASvar_Define.f90 60152 2015-08-13 19:19:13Z paul.vandelst@noaa.gov $'
  ! Release and version
  INTEGER, PARAMETER :: ASVAR_RELEASE = 1  ! This determines structure and file formats.
  INTEGER, PARAMETER :: ASVAR_VERSION = 1  ! This is just the default data version.
  ! Close status for write errors
  CHARACTER(*), PARAMETER :: WRITE_ERROR_STATUS = 'DELETE'
  ! Literal constants
  REAL(fp), PARAMETER :: ZERO = 0.0_fp
  REAL(fp), PARAMETER :: ONE  = 1.0_fp
  ! String lengths
  INTEGER,  PARAMETER :: ML = 256 ! Message length
  INTEGER,  PARAMETER :: SL =  80 ! String length

  
  ! ---------------------
  ! Structure definitions
  ! ---------------------
  ! The interpolation routine structure
  TYPE :: ASinterp_type
    ! The interpolating polynomials
    TYPE(LPoly_type) :: wlp  ! Frequency
    TYPE(LPoly_type) :: xlp  ! Effective radius
    ! The LUT interpolation indices
    INTEGER :: i1, i2        ! Frequency
    INTEGER :: j1, j2        ! Effective radius
    ! The LUT interpolation boundary check
    LOGICAL :: f_outbound    ! Frequency
    LOGICAL :: r_outbound    ! Effective radius
    ! The interpolation input
    REAL(fp) :: f_int        ! Frequency
    REAL(fp) :: r_int        ! Effective radius
    ! The data to be interpolated
    REAL(fp) :: f(NPTS)      ! Frequency
    REAL(fp) :: r(NPTS)      ! Effective radius
  END TYPE ASinterp_type
  
  
  ! The internal variable definition to hold information
  ! between FWD, TL, AD, and K-matrix calls
  TYPE :: ASvar_type
    ! Allocation indicator
    LOGICAL :: Is_Allocated = .FALSE.
    ! Release and version information
    INTEGER :: Release = ASVAR_RELEASE
    INTEGER :: Version = ASVAR_VERSION
    ! Dimensions
    INTEGER :: n_Legendre_Terms = 0  ! I1
    INTEGER :: n_Phase_Elements = 0  ! I2
    INTEGER :: n_Layers         = 0  ! I3
    INTEGER :: n_Aerosols       = 0  ! I4
    ! The interpolating data
    TYPE(ASinterp_type), ALLOCATABLE :: asi(:,:)  ! I3 x I4
    ! The interpolation results
    REAL(fp), ALLOCATABLE :: ke(:,:)          ! I3 x I4  Mass extinction coefficient
    REAL(fp), ALLOCATABLE :: w(:,:)           ! I3 x I4  Single Scatter Albedo
    REAL(fp), ALLOCATABLE :: g(:,:)           ! I3 x I4  Asymmetry factor
    REAL(fp), ALLOCATABLE :: pcoeff(:,:,:,:)  ! 0:I1 x I2 x I3 x I4  Phase coefficients
    ! The accumulated scattering coefficient
    REAL(fp), ALLOCATABLE :: total_bs(:)      ! I3  Volume scattering coefficient
  END TYPE ASvar_type


CONTAINS


!################################################################################
!################################################################################
!##                                                                            ##
!##                           ## PUBLIC PROCEDURES ##                          ##
!##                                                                            ##
!################################################################################
!################################################################################

  ELEMENTAL FUNCTION ASvar_Associated( self ) RESULT( Status )
    TYPE(ASvar_type), INTENT(IN) :: self
    LOGICAL :: Status
    Status = self%Is_Allocated
  END FUNCTION ASvar_Associated

 
  ELEMENTAL SUBROUTINE ASvar_Destroy( self )
    TYPE(ASvar_type), INTENT(OUT) :: self
    self%Is_Allocated = .FALSE.
    self%n_Legendre_Terms = 0
    self%n_Phase_Elements = 0
    self%n_Layers         = 0
    self%n_Aerosols       = 0
  END SUBROUTINE ASvar_Destroy
  
  
  ELEMENTAL SUBROUTINE ASvar_Create( &
    self            , &  ! Output
    n_Legendre_Terms, &  ! Input
    n_Phase_Elements, &  ! Input
    n_Layers        , &  ! Input
    n_Aerosols        )  ! Input
    ! Arguments
    TYPE(ASvar_type), INTENT(OUT) :: self
    INTEGER              , INTENT(IN)  :: n_Legendre_Terms        
    INTEGER              , INTENT(IN)  :: n_Phase_Elements             
    INTEGER              , INTENT(IN)  :: n_Layers                
    INTEGER              , INTENT(IN)  :: n_Aerosols                   
    ! Local variables
    INTEGER :: alloc_stat

    ! Check input
    IF ( n_Legendre_Terms < 1 .OR. &
         n_Phase_Elements < 1 .OR. &
         n_Layers         < 1 .OR. &
         n_Aerosols       < 1 ) RETURN

    ! Perform the allocation
    ALLOCATE( self%asi(n_Layers, n_Aerosols), &
              self%ke(n_Layers, n_Aerosols), &
              self%w(n_Layers, n_Aerosols), &
              self%g(n_Layers, n_Aerosols), &
              self%pcoeff(0:n_Legendre_Terms,n_Phase_Elements,n_Layers, n_Aerosols), &
              self%total_bs(n_Layers), &
              STAT = alloc_stat )
    IF ( alloc_stat /= 0 ) RETURN


    ! Initialise dimensions only!
    self%n_Legendre_Terms = n_Legendre_Terms
    self%n_Phase_Elements = n_Phase_Elements
    self%n_Layers         = n_Layers        
    self%n_Aerosols       = n_Aerosols      

    ! Set allocation indicator
    self%Is_Allocated = .TRUE.
  END SUBROUTINE ASvar_Create
  
  
  SUBROUTINE ASvar_Inspect( self)
    TYPE(ASvar_type), INTENT(IN) :: self
    INTEGER :: i2, i3, i4
    WRITE(*,'(1x,"ASvar OBJECT")')

    ! Release/version info
    WRITE(*,'(3x,"Release.Version     :",1x,i0,".",i0)') self%Release, self%Version

    ! Dimensions
    WRITE(*,'(3x,"n_Legendre_Terms    :",1x,i0)') self%n_Legendre_Terms
    WRITE(*,'(3x,"n_Phase_Elements    :",1x,i0)') self%n_Phase_Elements
    WRITE(*,'(3x,"n_Layers            :",1x,i0)') self%n_Layers        
    WRITE(*,'(3x,"n_Aerosols          :",1x,i0)') self%n_Aerosols      
    IF ( .NOT. ASvar_Associated(self) ) RETURN

    ! Data
    WRITE(*,'(3x,"Mass extinction coefficient (ke) :")')
    DO i4 = 1, self%n_Aerosols
      WRITE(*,'(5x,"ke Aerosol index #",i0)') i4
      WRITE(*,'(5(1x,es13.6,:))') self%ke(:,i4)
    END DO
    WRITE(*,'(3x,"Single scatter albedo (w) :")')
    DO i4 = 1, self%n_Aerosols
      WRITE(*,'(5x,"w Aerosol index #",i0)') i4
      WRITE(*,'(5(1x,es13.6,:))') self%w(:,i4)
    END DO
    WRITE(*,'(3x,"Asymmetry factor (g) :")')
    DO i4 = 1, self%n_Aerosols
      WRITE(*,'(5x,"g Aerosol index #",i0)') i4
      WRITE(*,'(5(1x,es13.6,:))') self%g(:,i4)
    END DO
    WRITE(*,'(3x,"Phase coefficients (pcoeff) :")')
    DO i4 = 1, self%n_Aerosols
      WRITE(*,'(5x,"pcoeff Aerosol index #",i0)') i4
      DO i3 = 1, self%n_Layers
        WRITE(*,'(7x,"pcoeff Layer index #",i0)') i3
        DO i2 = 1, self%n_Phase_Elements
          WRITE(*,'(9x,"pcoeff Phase element index #",i0)') i2
          WRITE(*,'(5(1x,es13.6,:))') self%pcoeff(0:,i2,i3,i4)
        END DO
      END DO
    END DO
    WRITE(*,'(3x,"Volume scattering coefficient (total_bs) :")')
    WRITE(*,'(5(1x,es13.6,:))') self%total_bs
  END SUBROUTINE ASvar_Inspect


  FUNCTION ASvar_ValidRelease( self ) RESULT( IsValid )
    ! Arguments
    TYPE(ASvar_type), INTENT(IN) :: self
    ! Function result
    LOGICAL :: IsValid
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'ASvar_ValidRelease'
    ! Local variables
    CHARACTER(ML) :: msg

    ! Set up
    IsValid = .TRUE.

    ! Check release is not too old
    IF ( self%Release < ASVAR_RELEASE ) THEN
      IsValid = .FALSE.
      WRITE( msg,'("An ASvar data update is needed. ", &
                  &"ASvar release is ",i0,". Valid release is ",i0,"." )' ) &
                  self%Release, ASVAR_RELEASE
      CALL Display_Message( ROUTINE_NAME, msg, INFORMATION ); RETURN
    END IF


    ! Check release is not too new
    IF ( self%Release > ASVAR_RELEASE ) THEN
      IsValid = .FALSE.
      WRITE( msg,'("An ASvar software update is needed. ", &
                  &"ASvar release is ",i0,". Valid release is ",i0,"." )' ) &
                  self%Release, ASVAR_RELEASE
      CALL Display_Message( ROUTINE_NAME, msg, INFORMATION ); RETURN
    END IF
  END FUNCTION ASvar_ValidRelease


  SUBROUTINE ASvar_Info( self, Info )
    ! Arguments
    TYPE(ASvar_type), INTENT(IN)  :: self
    CHARACTER(*),     INTENT(OUT) :: Info
    ! Parameters
    INTEGER, PARAMETER :: CARRIAGE_RETURN = 13
    INTEGER, PARAMETER :: LINEFEED = 10
    ! Local variables
    CHARACTER(2000) :: Long_String

    ! Write the required data to the local string
    WRITE( Long_String, &
           '(a,1x,"ASvar RELEASE.VERSION: ",i2,".",i2.2,a,3x, &
           &"N_LEGENDRE_TERMS=",i0,2x,&
           &"N_PHASE_ELEMENTS=",i0,2x,&
           &"N_LAYERS=",i0,2x,&
           &"N_AEROSOLS=",i0 )' ) &
           ACHAR(CARRIAGE_RETURN)//ACHAR(LINEFEED), &
           self%Release, self%Version, &
           ACHAR(CARRIAGE_RETURN)//ACHAR(LINEFEED), &
           self%n_Legendre_Terms, &
           self%n_Phase_Elements, &
           self%n_Layers        , &
           self%n_Aerosols      
                       
    ! Trim the output based on the
    ! dummy argument string length
    Info = Long_String(1:MIN(LEN(Info), LEN_TRIM(Long_String)))
  END SUBROUTINE ASvar_Info


  SUBROUTINE ASvar_DefineVersion( Id )
    CHARACTER(*), INTENT(OUT) :: Id
    Id = MODULE_VERSION_ID
  END SUBROUTINE ASvar_DefineVersion


  FUNCTION ASvar_InquireFile( &
    Filename        , &  ! Input
    n_Legendre_Terms, &  ! Optional output  
    n_Phase_Elements, &  ! Optional output  
    n_Layers        , &  ! Optional output  
    n_Aerosols      , &  ! Optional output  
    Release         , &  ! Optional output
    Version         , &  ! Optional output
    Title           , &  ! Optional output
    History         , &  ! Optional output
    Comment         ) &  ! Optional output
  RESULT( err_stat )
    ! Arguments
    CHARACTER(*),           INTENT(IN)  :: Filename
    INTEGER     , OPTIONAL, INTENT(OUT) :: n_Legendre_Terms
    INTEGER     , OPTIONAL, INTENT(OUT) :: n_Phase_Elements
    INTEGER     , OPTIONAL, INTENT(OUT) :: n_Layers        
    INTEGER     , OPTIONAL, INTENT(OUT) :: n_Aerosols      
    INTEGER     , OPTIONAL, INTENT(OUT) :: Release        
    INTEGER     , OPTIONAL, INTENT(OUT) :: Version        
    CHARACTER(*), OPTIONAL, INTENT(OUT) :: Title           
    CHARACTER(*), OPTIONAL, INTENT(OUT) :: History         
    CHARACTER(*), OPTIONAL, INTENT(OUT) :: Comment         
    ! Function result
    INTEGER :: err_stat
    ! Function parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'ASvar_InquireFile'
    ! Function variables
    CHARACTER(ML) :: msg
    CHARACTER(ML) :: io_msg
    INTEGER :: io_stat
    INTEGER :: fid
    TYPE(ASvar_type) :: ASvar

 
    ! Setup
    err_stat = SUCCESS
    ! ...Check that the file exists
    IF ( .NOT. File_Exists( Filename ) ) THEN
      msg = 'File '//TRIM(Filename)//' not found.'
      CALL Inquire_Cleanup(); RETURN
    END IF


    ! Open the file
    err_stat = Open_Binary_File( Filename, fid )
    IF ( err_stat /= SUCCESS ) THEN
      msg = 'Error opening '//TRIM(Filename)
      CALL Inquire_Cleanup(); RETURN
    END IF


    ! Read the release and version
    READ( fid, IOSTAT=io_stat, IOMSG=io_msg ) &
      ASvar%Release, &
      ASvar%Version
    IF ( io_stat /= 0 ) THEN
      msg = 'Error reading Release/Version - '//TRIM(io_msg)
      CALL Inquire_Cleanup(); RETURN
    END IF
    IF ( .NOT. ASvar_ValidRelease( ASvar ) ) THEN
      msg = 'ASvar Release check failed.'
      CALL Inquire_Cleanup(); RETURN
    END IF


    ! Read the dimensions
    READ( fid, IOSTAT=io_stat, IOMSG=io_msg ) &
      ASvar%n_Legendre_Terms, &
      ASvar%n_Phase_Elements, &
      ASvar%n_Layers        , &
      ASvar%n_Aerosols      
    IF ( io_stat /= 0 ) THEN
      msg = 'Error reading dimension values from '//TRIM(Filename)//' - '//TRIM(io_msg)
      CALL Inquire_Cleanup(); RETURN
    END IF


    ! Read the global attributes
    err_stat = ReadGAtts_Binary_File( &
                 fid, &
                 Title   = Title  , &
                 History = History, &
                 Comment = Comment  )
    IF ( err_stat /= SUCCESS ) THEN
      msg = 'Error reading global attributes'
      CALL Inquire_Cleanup(); RETURN
    END IF


    ! Close the file
    CLOSE( fid, IOSTAT=io_stat, IOMSG=io_msg )
    IF ( io_stat /= 0 ) THEN
      msg = 'Error closing '//TRIM(Filename)//' - '//TRIM(io_msg)
      CALL Inquire_Cleanup(); RETURN
    END IF


    ! Assign the return arguments
    IF ( PRESENT(n_Legendre_Terms) ) n_Legendre_Terms = ASvar%n_Legendre_Terms
    IF ( PRESENT(n_Phase_Elements) ) n_Phase_Elements = ASvar%n_Phase_Elements    
    IF ( PRESENT(n_Layers        ) ) n_Layers         = ASvar%n_Layers        
    IF ( PRESENT(n_Aerosols      ) ) n_Aerosols       = ASvar%n_Aerosols          
    IF ( PRESENT(Release         ) ) Release          = ASvar%Release        
    IF ( PRESENT(Version         ) ) Version          = ASvar%Version        
    
  CONTAINS
  
    SUBROUTINE Inquire_CleanUp()
      ! Close file if necessary
      IF ( File_Open(fid) ) THEN
        CLOSE( fid, IOSTAT=io_stat, IOMSG=io_msg )
        IF ( io_stat /= 0 ) &
          msg = TRIM(msg)//'; Error closing input file during error cleanup - '//TRIM(io_msg)
      END IF
      ! Set error status and print error message
      err_stat = FAILURE
      CALL Display_Message( ROUTINE_NAME, msg, err_stat )
    END SUBROUTINE Inquire_CleanUp
    
  END FUNCTION ASvar_InquireFile


  FUNCTION ASvar_ReadFile( &
    ASvar     , &  ! Output
    Filename  , &  ! Input
    No_Close  , &  ! Optional input
    Quiet     , &  ! Optional input
    Title     , &  ! Optional output
    History   , &  ! Optional output
    Comment   , &  ! Optional output
    Debug     ) &  ! Optional input (Debug output control)
  RESULT( err_stat )
    ! Arguments
    TYPE(ASvar_type),       INTENT(OUT) :: ASvar
    CHARACTER(*),           INTENT(IN)  :: Filename
    LOGICAL,      OPTIONAL, INTENT(IN)  :: No_Close
    LOGICAL,      OPTIONAL, INTENT(IN)  :: Quiet
    CHARACTER(*), OPTIONAL, INTENT(OUT) :: Title
    CHARACTER(*), OPTIONAL, INTENT(OUT) :: History
    CHARACTER(*), OPTIONAL, INTENT(OUT) :: Comment
    LOGICAL,      OPTIONAL, INTENT(IN)  :: Debug
    ! Function result
    INTEGER :: err_stat
    ! Function parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'ASvar_ReadFile'
    ! Function variables
    CHARACTER(ML) :: msg
    CHARACTER(ML) :: io_msg
    LOGICAL :: close_file
    LOGICAL :: noisy
    INTEGER :: io_stat
    INTEGER :: fid
    TYPE(ASvar_type) :: dummy
    
    ! Setup
    err_stat = SUCCESS
    ! ...Check No_Close argument
    close_file = .TRUE.
    IF ( PRESENT(No_Close) ) close_file = .NOT. No_Close
    ! ...Check Quiet argument
    noisy = .TRUE.
    IF ( PRESENT(Quiet) ) noisy = .NOT. Quiet
    ! ...Override Quiet settings if debug set.
    IF ( PRESENT(Debug) ) THEN
      IF ( Debug ) noisy = .TRUE.
    END IF

   
    ! Check if the file is open.
    IF ( File_Open( Filename ) ) THEN
      ! ...Inquire for the logical unit number
      INQUIRE( FILE=Filename, NUMBER=fid )
      ! ...Ensure it's valid
      IF ( fid < 0 ) THEN
        msg = 'Error inquiring '//TRIM(Filename)//' for its FileID'
        CALL Read_CleanUp(); RETURN
      END IF
    ELSE
      ! ...Open the file if it exists
      IF ( File_Exists( Filename ) ) THEN
        err_stat = Open_Binary_File( Filename, fid )
        IF ( err_Stat /= SUCCESS ) THEN
          msg = 'Error opening '//TRIM(Filename)
          CALL Read_CleanUp(); RETURN
        END IF
      ELSE
        msg = 'File '//TRIM(Filename)//' not found.'
        CALL Read_CleanUp(); RETURN
      END IF
    END IF


    ! Read and check the release and version
    READ( fid, IOSTAT=io_stat, IOMSG=io_msg ) &
      dummy%Release, &
      dummy%Version
    IF ( io_stat /= 0 ) THEN
      msg = 'Error reading Release/Version - '//TRIM(io_msg)
      CALL Read_Cleanup(); RETURN
    END IF
    IF ( .NOT. ASvar_ValidRelease( dummy ) ) THEN
      msg = 'ASvar Release check failed.'
      CALL Read_Cleanup(); RETURN
    END IF


    ! Read the dimensions
    READ( fid, IOSTAT=io_stat, IOMSG=io_msg ) &
      dummy%n_Legendre_Terms, &
      dummy%n_Phase_Elements, &
      dummy%n_Layers        , &
      dummy%n_Aerosols      
    IF ( io_stat /= 0 ) THEN
      msg = 'Error reading data dimensions - '//TRIM(io_msg)
      CALL Read_Cleanup(); RETURN
    END IF
    ! ...Allocate the object
    CALL ASvar_Create( &
           ASvar                 , &
           dummy%n_Legendre_Terms, &        
           dummy%n_Phase_Elements, &        
           dummy%n_Layers        , &        
           dummy%n_Aerosols        )                  
    IF ( .NOT. ASvar_Associated( ASvar ) ) THEN
      msg = 'ASvar object allocation failed.'
      CALL Read_Cleanup(); RETURN
    END IF
    ! ...Explicitly assign the version number
    ASvar%Version = dummy%Version
        

    ! Read the global attributes
    err_stat = ReadGAtts_Binary_File( &
                 fid, &
                 Title   = Title  , &
                 History = History, &
                 Comment = Comment  )
    IF ( err_stat /= SUCCESS ) THEN
      msg = 'Error reading global attributes'
      CALL Read_Cleanup(); RETURN
    END IF


    ! Read the data
    ! ...Mass extinction coefficient
    READ( fid, IOSTAT=io_stat, IOMSG=io_msg ) &
      ASvar%ke
    IF ( io_stat /= 0 ) THEN
      msg = 'Error reading mass extinction coefficient - '//TRIM(io_msg)
      CALL Read_Cleanup(); RETURN
    END IF
    ! ...Single scatter albedo
    READ( fid, IOSTAT=io_stat, IOMSG=io_msg ) &
      ASvar%w
    IF ( io_stat /= 0 ) THEN
      msg = 'Error reading single scatter albedo - '//TRIM(io_msg)
      CALL Read_Cleanup(); RETURN
    END IF
    ! ...Asymmetry factor
    READ( fid, IOSTAT=io_stat, IOMSG=io_msg ) &
      ASvar%g
    IF ( io_stat /= 0 ) THEN
      msg = 'Error reading asymmetry factor - '//TRIM(io_msg)
      CALL Read_Cleanup(); RETURN
    END IF
    ! ...Phase coefficients
    READ( fid, IOSTAT=io_stat, IOMSG=io_msg ) &
      ASvar%pcoeff
    IF ( io_stat /= 0 ) THEN
      msg = 'Error reading phase coefficients - '//TRIM(io_msg)
      CALL Read_Cleanup(); RETURN
    END IF
    ! ...Total volume scattering coefficient
    READ( fid, IOSTAT=io_stat, IOMSG=io_msg ) &
      ASvar%total_bs
    IF ( io_stat /= 0 ) THEN
      msg = 'Error reading total volume scattering coefficient - '//TRIM(io_msg)
      CALL Read_Cleanup(); RETURN
    END IF


    ! Close the file
    IF ( close_file ) THEN
      CLOSE( fid, IOSTAT=io_stat, IOMSG=io_msg )
      IF ( io_stat /= 0 ) THEN
        msg = 'Error closing '//TRIM(Filename)//' - '//TRIM(io_msg)
        CALL Read_Cleanup(); RETURN
      END IF
    END IF


    ! Output an info message
     IF ( noisy ) THEN
       CALL ASvar_Info( ASvar, msg )
       CALL Display_Message( ROUTINE_NAME, 'FILE: '//TRIM(Filename)//'; '//TRIM(msg), INFORMATION )
     END IF

   CONTAINS
   
     SUBROUTINE Read_CleanUp()
       IF ( File_Open(Filename) ) THEN
         CLOSE( fid, IOSTAT=io_stat, IOMSG=io_msg )
         IF ( io_stat /= 0 ) &
           msg = TRIM(msg)//'; Error closing input file during error cleanup - '//TRIM(io_msg)
       END IF
       CALL ASvar_Destroy( ASvar )
       err_stat = FAILURE
       CALL Display_Message( ROUTINE_NAME, msg, err_stat )
     END SUBROUTINE Read_CleanUp

  END FUNCTION ASvar_ReadFile


  FUNCTION ASvar_WriteFile( &
    ASvar     , &  ! Input
    Filename  , &  ! Input
    No_Close  , &  ! Optional input
    Quiet     , &  ! Optional input
    Title     , &  ! Optional input
    History   , &  ! Optional input
    Comment   , &  ! Optional input
    Debug     ) &  ! Optional input (Debug output control)
  RESULT( err_stat )
    ! Arguments
    TYPE(ASvar_type),  INTENT(IN) :: ASvar
    CHARACTER(*),           INTENT(IN) :: Filename
    LOGICAL,      OPTIONAL, INTENT(IN) :: No_Close
    LOGICAL,      OPTIONAL, INTENT(IN) :: Quiet
    CHARACTER(*), OPTIONAL, INTENT(IN) :: Title
    CHARACTER(*), OPTIONAL, INTENT(IN) :: History
    CHARACTER(*), OPTIONAL, INTENT(IN) :: Comment
    LOGICAL,      OPTIONAL, INTENT(IN) :: Debug
    ! Function result
    INTEGER :: err_stat
    ! Function parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'ASvar_WriteFile'
    ! Function variables
    CHARACTER(ML) :: msg
    CHARACTER(ML) :: io_msg
    LOGICAL :: close_file
    LOGICAL :: noisy
    INTEGER :: io_stat
    INTEGER :: fid
    

    ! Setup
    err_stat = SUCCESS
    ! ...Check No_Close argument
    close_file = .TRUE.
    IF ( PRESENT(No_Close) ) close_file = .NOT. No_Close
    ! ...Check Quiet argument
    noisy = .TRUE.
    IF ( PRESENT(Quiet) ) noisy = .NOT. Quiet
    ! ...Override Quiet settings if debug set.
    IF ( PRESENT(Debug) ) THEN
      IF ( Debug ) noisy = .TRUE.
    END IF
    ! ...Check there is data to write
    IF ( .NOT. ASvar_Associated( ASvar ) ) THEN
      msg = 'ASvar object is empty.'
      CALL Write_Cleanup(); RETURN
    END IF

   
    ! Check if the file is open.
    IF ( File_Open( FileName ) ) THEN
      ! ...Inquire for the logical unit number
      INQUIRE( FILE=Filename, NUMBER=fid )
      ! ...Ensure it's valid
      IF ( fid < 0 ) THEN
        msg = 'Error inquiring '//TRIM(Filename)//' for its FileID'
        CALL Write_CleanUp(); RETURN
      END IF
    ELSE
      ! ...Open the file for output
      err_stat = Open_Binary_File( Filename, fid, For_Output=.TRUE. )
      IF ( err_Stat /= SUCCESS ) THEN
        msg = 'Error opening '//TRIM(Filename)
        CALL Write_CleanUp(); RETURN
      END IF
    END IF


    ! Write the release and version
    WRITE( fid, IOSTAT=io_stat, IOMSG=io_msg ) &
      ASvar%Release, &
      ASvar%Version
    IF ( io_stat /= 0 ) THEN
      msg = 'Error writing Release/Version - '//TRIM(io_msg)
      CALL Write_Cleanup(); RETURN
    END IF


    ! Write the dimensions
    WRITE( fid, IOSTAT=io_stat, IOMSG=io_msg ) &
      ASvar%n_Legendre_Terms, &
      ASvar%n_Phase_Elements, &
      ASvar%n_Layers        , &
      ASvar%n_Aerosols      
    IF ( io_stat /= 0 ) THEN
      msg = 'Error writing data dimensions - '//TRIM(io_msg)
      CALL Write_Cleanup(); RETURN
    END IF


    ! Write the global attributes
    err_stat = WriteGAtts_Binary_File( &
                 fid, &
                 Write_Module = MODULE_VERSION_ID, &
                 Title        = Title  , &
                 History      = History, &
                 Comment      = Comment  )
    IF ( err_stat /= SUCCESS ) THEN
      msg = 'Error writing global attributes'
      CALL Write_Cleanup(); RETURN
    END IF


    ! Write the data
    ! ...Mass extinction coefficient
    WRITE( fid, IOSTAT=io_stat, IOMSG=io_msg ) &
      ASvar%ke
    IF ( io_stat /= 0 ) THEN
      msg = 'Error writing mass extinction coefficient - '//TRIM(io_msg)
      CALL Write_Cleanup(); RETURN
    END IF
    ! ...Single scatter albedo
    WRITE( fid, IOSTAT=io_stat, IOMSG=io_msg ) &
      ASvar%w
    IF ( io_stat /= 0 ) THEN
      msg = 'Error writing single scatter albedo - '//TRIM(io_msg)
      CALL Write_Cleanup(); RETURN
    END IF
    ! ...Asymmetry factor
    WRITE( fid, IOSTAT=io_stat, IOMSG=io_msg ) &
      ASvar%g
    IF ( io_stat /= 0 ) THEN
      msg = 'Error writing asymmetry factor - '//TRIM(io_msg)
      CALL Write_Cleanup(); RETURN
    END IF
    ! ...Phase coefficients
    WRITE( fid, IOSTAT=io_stat, IOMSG=io_msg ) &
      ASvar%pcoeff
    IF ( io_stat /= 0 ) THEN
      msg = 'Error writing phase coefficients - '//TRIM(io_msg)
      CALL Write_Cleanup(); RETURN
    END IF
    ! ...Total volume scattering coefficient
    WRITE( fid, IOSTAT=io_stat, IOMSG=io_msg ) &
      ASvar%total_bs
    IF ( io_stat /= 0 ) THEN
      msg = 'Error writing total volume scattering coefficient - '//TRIM(io_msg)
      CALL Write_Cleanup(); RETURN
    END IF


    ! Close the file
    IF ( close_file ) THEN
      CLOSE( fid, IOSTAT=io_stat, IOMSG=io_msg )
      IF ( io_stat /= 0 ) THEN
        msg = 'Error closing '//TRIM(Filename)//' - '//TRIM(io_msg)
        CALL Write_Cleanup(); RETURN
      END IF
    END IF


    ! Output an info message
     IF ( noisy ) THEN
       CALL ASvar_Info( ASvar, msg )
       CALL Display_Message( ROUTINE_NAME, 'FILE: '//TRIM(Filename)//'; '//TRIM(msg), INFORMATION )
     END IF

   CONTAINS
   
     SUBROUTINE Write_Cleanup()
       IF ( File_Open(Filename) ) THEN
         CLOSE( fid, IOSTAT=io_stat, IOMSG=io_msg )
         IF ( io_stat /= 0 ) &
           msg = TRIM(msg)//'; Error closing output file during error cleanup - '//TRIM(io_msg)
       END IF
       err_stat = FAILURE
       CALL Display_Message( ROUTINE_NAME, msg, err_stat )
     END SUBROUTINE Write_Cleanup

  END FUNCTION ASvar_WriteFile

  
!################################################################################
!################################################################################
!##                                                                            ##
!##                          ## PRIVATE PROCEDURES ##                          ##
!##                                                                            ##
!################################################################################
!################################################################################

  ELEMENTAL FUNCTION ASvar_Equal( x, y ) RESULT( is_equal )
    TYPE(ASvar_type), INTENT(IN) :: x, y
    LOGICAL :: is_equal

    ! Set up
    is_equal = .FALSE.
   
    ! Check the object association status
    IF ( (.NOT. ASvar_Associated(x)) .OR. &
         (.NOT. ASvar_Associated(y))      ) RETURN

   ! Check contents
    ! ...Release/version info
    IF ( (x%Release /= y%Release) .OR. &
         (x%Version /= y%Version) ) RETURN
    ! ...Dimensions
    IF ( (x%n_Legendre_Terms /= y%n_Legendre_Terms ) .OR. &
         (x%n_Phase_Elements /= y%n_Phase_Elements ) .OR. &
         (x%n_Layers         /= y%n_Layers         ) .OR. &
         (x%n_Aerosols       /= y%n_Aerosols       ) ) RETURN
    ! ...Arrays
    IF ( ALL(x%ke       .EqualTo. y%ke       ) .AND. &
         ALL(x%w        .EqualTo. y%w        ) .AND. &
         ALL(x%g        .EqualTo. y%g        ) .AND. &
         ALL(x%pcoeff   .EqualTo. y%pcoeff   ) .AND. &
         ALL(x%total_bs .EqualTo. y%total_bs ) ) &
      is_equal = .TRUE.
  END FUNCTION ASvar_Equal

END MODULE ASvar_Define
