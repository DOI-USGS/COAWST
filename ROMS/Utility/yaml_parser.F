      MODULE yaml_parser_mod
!
!git $Id$
!svn $Id: yaml_parser.F 1151 2023-02-09 03:08:53Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2023 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This module contains several routines to process input YAML files.  !
!                                                                      !
!  Notice that several Fortran parsers exist for complex and simple    !
!  YAML files coded with Object-Oriented Programming (OOP) principles. !
!  For example, check:                                                 !
!                                                                      !
!   * FCKit  (https://github.com/ecmwf/fckit)                          !
!   * Fortran-YAML (https://github.com/BoldingBruggeman/fortran-yaml)  !
!   * yaFyaml (https://github.com/Goddard-Fortran-Ecosystem/yaFyaml)   !
!                                                                      !
!  However, this YAML parser is more uncomplicated with substantial    !
!  capabilities. It is a hybrid between standard and OOP principles    !
!  but without the need for recurrency, inheritance, polymorphism,     !
!  and containers.                                                     !
!                                                                      !
!  The only constraint is that the YAML file is read twice for         !
!  simplicity. The first read determines the number indentation of     !
!  blanks policy and the length of the collection list(:) pairs object !
!  (CLASS yaml_pair).  It supports:                                    !
!                                                                      !
!   * Single or multiple line comments start with a hash '#'. Also,    !
!     comment after a key/value pair is allowed. All comments are      !
!     skipped during processing.                                       !
!   * Unlimited nested structure (lists, mappings, hierarchies).       !
!     Indentation of whitespace is used to denote structure.           !
!   * Unrestricted schema indentation. However, some schema validators !
!     recommend or impose two whitespace indentations.                 !
!   * A key is followed by a colon to denote a mapping value (like     !
!     ocean_model: ROMS).                                              !
!   * Aliases and Anchors.                                             !
!   * Blocking lists: members are denoted by a leading hyphen and      !
!     space, which is considered as part of the indentation.           !
!   * Flow sequence: a vector list with values enclosed in square      !
!     brackets and separated by a comma-and-space: [val1, ..., valN].  !
!   * Keyword values are processed and stored as strings but converted !
!     to a logical, integer, or floating-point type when appropriate   !
!     during extraction. If derived-type values are needed, the caller !
!     can process such structure outside this module, as shown below.  !
!   * Remove unwanted control characters like tabs and separators      !
!     (ASCII character code 0-31)                                      !
!   * English uppercase and lowercase alphabet, but it can be expanded !
!     to other characters (see yaml_ValueType routine)                 !
!   * Module is self contained, but it has very minimal association to !
!     four ROMS modules.                                               !
!   * Multiple or continuation lines are supported, for example we can !
!     have:                                                            !
!                                                                      !
!       state variables: [sea_surface_height_above_geoid,              !
!                         barotropic_sea_water_x_velocity,             !
!                         barotropic_sea_water_y_velocity,             !
!                         sea_water_x_velocity,                        !
!                         sea_water_y_velocity,                        !
!                         sea_water_potential_temperature,             !
!                         sea_water_practical_salinity]                !
!                                                                      !
!  Usage:                                                              !
!                                                                      !
!    USE yaml_parser_mod, ONLY : yaml_initialize, yaml_get, yaml_tree  !
!                                                                      !
!    TYPE (yaml_tree) :: YML                                           !
!                                                                      !
!    CALL yaml_initialize (YML, 'ocn_coupling.yaml', report)           !
!    status = yaml_get(YML, 'Nimport', Nimport)                        !
!    status = yaml_get(YML, 'import.standard_name', Sstandard)         !
!    status = yaml_get(YML, 'import.standard_name.short_name', S%short)!
!    status = yaml_get(YML, 'import.standard_name.unit', S%units)      !
!    ...                                                               !
!                                                                      !
!    and so on for logical, integer, floating-point, and string        !
!    key/value pairs.                                                  !
!                                                                      !
!    Here, 'Sstandard(1:Nimport) will contain all 'standard_name'      !
!    values for the YAML block 'import:'. Notice that nested keywords  !
!    are concatenated with a period: 'key1.key2.key3' for a three-     !
!    level nested block, similar to how Matlab build structures. The   !
!    key can have more than one word separated by one space. For       !
!    example, we can have 'bulk_flux import.standard_name'. Similarly, !
!    any otherkey/value pair can be extrated from the YML object.      !
!                                                                      !
!=======================================================================
!
!  ROMS kernel association.
!
      USE mod_kinds,    ONLY : kind_real    => dp     ! double-precision
      USE mod_parallel, ONLY : yaml_Master  => Master ! master PET
      USE mod_scalars,  ONLY : yaml_ErrFlag => exit_flag  ! error flag
      USE mod_iounits,  ONLY : yaml_stdout  => stdout ! standard ouput
!
      implicit none
!
!-----------------------------------------------------------------------
!  Structures/Objects to hold YAML dictionary lists with theirs keys and
!  values.
!-----------------------------------------------------------------------
!
!  YAML file key/value pair.
!
      TYPE, PUBLIC :: yaml_pair

        logical :: has_alias                      ! alias '*' token
        logical :: has_anchor                     ! anchor '&' token
        logical :: is_block                       ! block '-' list
        logical :: is_sequence                    ! sequence '[]' tokens
!
        logical :: is_logical                     ! logical value
        logical :: is_integer                     ! integer value
        logical :: is_real                        ! numerical value
        logical :: is_string                      ! string value
!
        integer :: id                             ! key/value ID
        integer :: parent_id                      ! parent ID
        integer :: left_padding                   ! indent level: 0,1,..
!
        character (len=:), allocatable :: line    ! YAML line
        character (len=:), allocatable :: key     ! YAML keyword:
        character (len=:), allocatable :: value   ! YAML value(s)
        character (len=:), allocatable :: anchor  ! anchor keyword

      END TYPE yaml_pair
!
!  YAML file dictionary tree.
!
      TYPE, PUBLIC :: yaml_tree

        logical :: IsCreated = .FALSE. ! Object creation switch

        integer :: Nbranches           ! total number of branches
        integer :: Npairs              ! total number of pairs
        integer :: indent              ! blank indentation policy
!
        character (len=:),  allocatable :: filename  ! YAML file name
!
!  YAML file collection pairs, [1:Npairs].
!
        TYPE (yaml_pair), pointer :: list(:)
!
        CONTAINS                       ! CLASS objects
!
        PROCEDURE :: create       => yaml_tree_create
        PROCEDURE :: destroy      => yaml_tree_destroy
        PROCEDURE :: dump         => yaml_tree_dump
        PROCEDURE :: extract      => yaml_tree_extract
        PROCEDURE :: fill         => yaml_tree_fill
        PROCEDURE :: fill_aliases => yaml_tree_fill_aliases
        PROCEDURE :: has          => yaml_tree_has
        PROCEDURE :: read_line    => yaml_tree_read_line

      END TYPE yaml_tree
!
!  Public structures that can be used in applications to extract block
!  list YAML constructs. The key may represents a sequence flow [...]
!  with a vector of values. The values can be integers, logicals, reals,
!  or strings.  For example,
!
!  import:
!    - standard_name:       surface_downward_heat_flux_in_sea_water
!      long_name:           surface net heat flux
!      short_name:          shflux
!      data_variables:      [shf, shf_time]
!    - standard_name:       surface_wind_x_stress
!      long_name:           surface zonal wind stress component
!      short_name:          sustr
!      data_variables:      [taux, atm_time]
!    - standard_name:       surface_wind_y_stress
!      long_name:           surface meridional wind stress component
!      short_name:          svstr
!      data_variables:      [tauy, atm_time]
!
!  The extraction is loaded into V, which is a TYPE "yaml_Svec"
!  structure:
!
!    status = yaml_get(YML, 'import.data_variables', V)
!
!  or altenatively
!
!    IF (yaml_Error(yaml_get(YML, 'import.data_variables', V),          &
!   &               NoErr, __LINE__, MyFile)) RETURN
!
!  yielding the following block-list, string structure in a single
!  invocation of the overloaded function "yaml_get":
!
!    V(1)%vector(1)%value = 'shf'
!    V(1)%vector(2)%value = 'shf_time'
!    V(2)%vector(1)%value = 'taux'
!    V(2)%vector(2)%value = 'atm_time'
!    V(3)%vector(1)%value = 'tauy'
!    V(3)%vector(2)%value = 'atm_time'
!
!  It is a compact way to extract similar data blocks.
!
      TYPE, PUBLIC :: yaml_Ivec                      ! integer structure
        integer, allocatable :: vector(:)            ! vector values
      END TYPE yaml_Ivec
!
      TYPE, PUBLIC :: yaml_Lvec                      ! logical structure
        logical, allocatable :: vector(:)            ! vector values
      END TYPE yaml_Lvec
!
      TYPE, PUBLIC :: yaml_Rvec                      ! real structure
        real (kind_real), allocatable :: vector(:)   ! vector values
      END TYPE yaml_Rvec
!
      TYPE, PUBLIC :: yaml_Svec                      ! string structure
        character (len=:), allocatable :: value      ! scalar value
        TYPE (yaml_Svec),  pointer :: vector(:)      ! recursive vector
      END TYPE yaml_Svec                             ! values
!
!  Derived-type structure, extended/inherited from parent "yaml_Svec",
!  is used to extract hierarchies of keys and associated values from
!  YAML dictionary object. The calling program specifies a key-string
!  that may be generated by aggregating nested keys with a period.
!  Also, it can extract flow sequence string element values that are
!  separated by commas.
!
      TYPE, PRIVATE, EXTENDS(yaml_Svec) :: yaml_extract
        logical :: has_vector             ! true if loaded vector values
      END TYPE yaml_extract
!
!-----------------------------------------------------------------------
!  Define public overloading API to extract key/value pairs from YAML
!  tree dictionary object accounting to variable type.
!-----------------------------------------------------------------------
!
      INTERFACE yaml_get
        MODULE PROCEDURE yaml_Get_i_struc       ! Gets integer structure
        MODULE PROCEDURE yaml_Get_l_struc       ! Gets logical structure
        MODULE PROCEDURE yaml_Get_r_struc       ! Gets real    structure
        MODULE PROCEDURE yaml_Get_s_struc       ! Gets string  structure
!
        MODULE PROCEDURE yaml_Get_ivar_0d       ! Gets integer value
        MODULE PROCEDURE yaml_Get_ivar_1d       ! Gest integer values
        MODULE PROCEDURE yaml_Get_lvar_0d       ! Gets logical value
        MODULE PROCEDURE yaml_Get_lvar_1d       ! Gets logical values
        MODULE PROCEDURE yaml_Get_rvar_0d       ! Gets real    value
        MODULE PROCEDURE yaml_Get_rvar_1d       ! Gets real    values
        MODULE PROCEDURE yaml_Get_svar_0d       ! Gets string  value
        MODULE PROCEDURE yaml_Get_svar_1d       ! Gets string  values
      END INTERFACE yaml_get
!
      PUBLIC  :: yaml_AssignString
      PUBLIC  :: yaml_Error
      PUBLIC  :: yaml_get
      PUBLIC  :: yaml_initialize
!
      PRIVATE :: yaml_CountKeys
      PRIVATE :: yaml_LowerCase
      PRIVATE :: yaml_UpperCase
      PRIVATE :: yaml_ValueType
!
      PRIVATE
!
!-----------------------------------------------------------------------
!  Local module parameters.
!-----------------------------------------------------------------------
!
      logical, parameter :: LdebugYAML = .FALSE.  ! debugging switch
!
      integer, parameter :: Ldim = 8    ! logical Lswitch dimension
      integer, parameter :: Lkey = 254  ! Maximum characters per keyword
      integer, parameter :: Lmax = 2048 ! Maximum characters per line
      integer, parameter :: NoErr = 0   ! no error flag
      integer, parameter :: iunit = 222 ! Fortan unit for reading
!
      character (len=55), dimension(7) :: yaml_ErrMsg =                 &
     &  (/ ' YAML_PARSER - Blows up ................ yaml_ErrFlag: ',   &
     &     ' YAML_PARSER - Input error ............. yaml_ErrFlag: ',   &
     &     ' YAML_PARSER - Output error ............ yaml_ErrFlag: ',   &
     &     ' YAML_PARSER - I/O error ............... yaml_ErrFlag: ',   &
     &     ' YAML_PARSER - Configuration error ..... yaml_ErrFlag: ',   &
     &     ' YAML_PARSER - Partition error ......... yaml_ErrFlag: ',   &
     &     ' YAML_PARSER - Illegal input parameter . yaml_ErrFlag: ' /)
!
!-----------------------------------------------------------------------
      CONTAINS
!-----------------------------------------------------------------------
!
      FUNCTION yaml_initialize (self, filename, report) RESULT (status)
!
!***********************************************************************
!                                                                      !
!  It creates a YAML tree dictionary object. First, it reads the YAML  !
!  file to determine the indentation policy and length of list oject   !
!  (TYPE yaml_tree).                                                   !
!                                                                      !
!  After the object is allocated, the Fortran unit is rewinded and the !
!  YAML file is read again to populate the keyword values.             !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     self       YAML tree dictionary object (CLASS yaml_tree)         !
!     filename   YAML filename (string)                                !
!     report     Switch to dump to standard output (logical, OPTIONAL) !
!                                                                      !
!  On Ouptut:                                                          !
!                                                                      !
!     self       Allocated and populated YAML object.                  !
!     status     Error flag (integer)                                  !
!                                                                      !
!***********************************************************************
!
!  Imported variable declarations.
!
      CLASS (yaml_tree), intent(inout) :: self
      character (len=*), intent(in   ) :: filename
      logical, optional, intent(in   ) :: report
!
!  Local variable declarations.
!
      logical :: Ldump
!
      integer :: LenStr, status
!
      character (len=*), parameter :: MyFile =                          &
     &   __FILE__//", yaml_initialize"
!
!-----------------------------------------------------------------------
!  Initialize YAML object.
!-----------------------------------------------------------------------
!
      status=Noerr
!
!  Set switch to print the processed YAML key/value pairs to standard
!  ouput.
!
      IF (PRESENT(report)) THEN
        Ldump = report
      ELSE
        Ldump = .FALSE.
      END IF
!
!  Set YAML file path and name.
!
      IF (yaml_Error(yaml_AssignString(self%filename,                   &
     &                                 filename, LenStr),               &
     &               NoErr, __LINE__, MyFile)) THEN
        status=yaml_ErrFlag
        RETURN
      END IF
!
!  Create and populate YAML object.
!
      IF (self%IsCreated) CALL self%destroy ()
!
      CALL self%create ()
      IF (yaml_Error(yaml_ErrFlag, NoErr, __LINE__, MyFile)) THEN
        status=yaml_ErrFlag
        RETURN
      END IF
!
!  Report YAML tree dictionary, for debugging.
!
      IF (Ldump) CALL self%dump ()
!
      RETURN
      END FUNCTION yaml_initialize
!
      SUBROUTINE yaml_tree_create (self)
!
!***********************************************************************
!                                                                      !
!  It creates a YAML tree dictionary object. First, it reads the YAML  !
!  file to determine the dimensions of parent and children structures. !
!  After the structures are allocate, the Fortran unit is rewinded and !
!  YAML file is read again to populate the keyword values.             !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     self       YAML tree dictionary object (CLASS yaml_tree)         !
!                                                                      !
!  On Ouptut:                                                          !
!                                                                      !
!     self       Allocated and populated YAML object.                  !
!                                                                      !
!***********************************************************************
!
!  Imported variable declarations.
!
      CLASS (yaml_tree), intent(inout) :: self
!
!  Local variable declarations.
!
      logical :: FirstPass, Lswitch(Ldim)
!
      integer :: Nblanks, Nbranches, Npairs
      integer :: i, io_err, status
!
      character (len=Lkey) :: anchor, io_errmsg, key
      character (len=Lmax) :: line, value
!
      character (len=*), parameter :: MyFile =                          &
     &   __FILE__//", yaml_tree_create"
!
!-----------------------------------------------------------------------
!  Open YAML file for reading.
!-----------------------------------------------------------------------
!
      OPEN (UNIT=iunit, FILE=self%filename, FORM='formatted',           &
            STATUS='old', IOSTAT=io_err, IOMSG=io_errmsg)
      IF (io_err.ne.0) THEN
        yaml_ErrFlag=5
        IF (yaml_Error(yaml_ErrFlag, NoErr, __LINE__, MyFile)) THEN
          IF (yaml_Master) WRITE (yaml_stdout,10) self%filename,        &
    &                                             TRIM(io_errmsg)
          RETURN
        END IF
      END IF
!
!  Determine the total number of YAML key/value pairs.
!
      Nblanks   = 0
      Nbranches = 0
      Npairs    = 0
!
      FirstPass = .TRUE.
!
      YAML_LINE : DO WHILE (.TRUE.)
!
        status=self%read_line(Nblanks, line, key, value, anchor,       &
     &                        Lswitch)
!
        SELECT CASE (status)
          CASE (-1)       ! error condition during reading
            yaml_ErrFlag=5
            IF (yaml_Error(yaml_ErrFlag, NoErr, __LINE__, MyFile)) THEN
              IF (yaml_Master) WRITE (yaml_stdout,20) self%filename,    &
     &                                                TRIM(line)
              RETURN
            END IF
          CASE (0)        ! end-of-file
            EXIT
          CASE (1)        ! blank or comment line, move to the next line
            CYCLE
          CASE (2)        ! processed viable line
            Npairs=Npairs+1
        END SELECT
!
!  If no leading blanks, advance YAML tree branch counter.
!
        IF (Nblanks.eq.0) THEN
          Nbranches = Nbranches+1          ! hierarchy start
        END IF
!
!  On first pass, establish indentation policy: number of blanks.
!
        IF (FirstPass.and.(Nblanks.gt.0)) THEN
          FirstPass=.FALSE.
          self%indent=Nblanks
        END IF
!
!  Check for consitent identation policy.  Some YAML validators impose
!  a two-blank spacing policy.
!
        IF (Nblanks.gt.0) THEN
          IF (MOD(Nblanks, self%indent).ne.0) THEN
            yaml_ErrFlag=2
            IF (yaml_Error(yaml_ErrFlag, NoErr, __LINE__, MyFile)) THEN
              IF (yaml_Master) WRITE (yaml_stdout,30) nblanks,          &
     &                                                self%indent
              RETURN
            END IF
          END IF
        END IF

      END DO YAML_LINE
!
!  Rewind unit since we need to reprocess the file again to load the
!  data into the allocated "self%list(1:Npairs)" container.
!
      REWIND (iunit)
!
!-----------------------------------------------------------------------
!  Allocate YAML dictionary container.
!-----------------------------------------------------------------------
!
!  Set number of main branches and number of key/values in YAML 'list'.
!
      self%Nbranches=Nbranches
      self%Npairs=Npairs
!
!  Allocate YAML key/value pair list (TYPE 'yaml_pair') object.
!
      IF (.not.self%IsCreated) THEN
        ALLOCATE ( self%list(Npairs) )
        self%IsCreated=.TRUE.
      END IF
!
!-----------------------------------------------------------------------
!  Re-read YAML file again and populate dictionary.
!-----------------------------------------------------------------------
!
      CALL self%fill ()
!
  10  FORMAT (/,' YAML_TREE_CREATE - Unable to open input YAML file: ', &
     &        a,/,20x,'ERROR: ',a)
  20  FORMAT (/,' YAML_TREE_CREATE - Error while reading YAML file: ',  &
     &        a,/,20x,'LINE: ',a)
  30  FORMAT (/,' YAML_TREE_CREATE - Inconsistent indentation, ',       &
     &        'self%indent = ',i0,',  nblanks = ',i0)
  40  FORMAT (/,' YAML_TREE_CREATE - Inconsistent indentation, ',       &
     &        'nblanks = ',i0,', indent blank policy = ',i0,/,20x,      &
     &        'Number of nested collections = ',i0,                     &
     &        ' is greater than 3',/,20x,'Line: ',a)
  50  FORMAT (/,' YAML_TREE_CREATE - Cannot find branches in YAML ',    &
     &        'file: ',a)
!
      RETURN
      END SUBROUTINE yaml_tree_create
!
      SUBROUTINE yaml_tree_destroy (self)
!
!***********************************************************************
!                                                                      !
!  It deallocates/destroys the YAML dictionary object.                 !
!                                                                      !
!  Arguments:                                                          !
!                                                                      !
!     self       YAML tree dictionary object (CLASS yaml_tree)         !
!                                                                      !
!***********************************************************************
!
!  Imported variable declarations.
!
      CLASS (yaml_tree), intent(inout) :: self
!
!  Local variable declarations.
!
      logical :: Lopened
!
!-----------------------------------------------------------------------
!  Deallocate YAML dictionary object.
!-----------------------------------------------------------------------
!
!  If applicable, close YAML Fortran unit.
!
      INQUIRE (UNIT=iunit, OPENED=Lopened)
      IF (Lopened) THEN
        CLOSE (iunit)
      END IF
!
!  Recall that Fortran 2003 standard allows deallocating just the
!  parent object to deallocate all array variables within its scope
!  automatically.
!
      IF (ASSOCIATED(self%list)) THEN
        self%IsCreated=.FALSE.
        DEALLOCATE (self%list)
      END IF
!
      RETURN
      END SUBROUTINE yaml_tree_destroy
!
      SUBROUTINE yaml_tree_dump (self)
!
!***********************************************************************
!                                                                      !
!  It prints the YAML dictionary to standard output for debugging      !
!  purposes.                                                           !
!                                                                      !
!  Arguments:                                                          !
!                                                                      !
!     self       YAML tree dictionary object (CLASS yaml_tree)         !
!                                                                      !
!***********************************************************************
!
!  Imported variable declarations.
!
      CLASS (yaml_tree), intent(in) :: self
!
!  Local variable declarations.
!
      integer :: Lstr, Nblanks, i, indent, padding
!
      character (len=Lkey) :: key
      character (len=Lmax) :: string, value
!
!-----------------------------------------------------------------------
!  Report the contents of the YAML tree directory.
!-----------------------------------------------------------------------
!
      IF (yaml_Master) THEN
        WRITE (yaml_stdout,10) self%filename
        indent=self%indent
!
        DO i=1,self%Npairs
          padding=self%list(i)%left_padding
          Nblanks=indent*padding
          key=self%list(i)%key
          IF (ALLOCATED(self%list(i)%value)) THEN
            value=self%list(i)%value
            Lstr=LEN_TRIM(value)
          ELSE
            Lstr=0
          END IF
          IF (self%list(i)%is_block) THEN
            IF (Lstr.gt.0) THEN
              IF (self%list(i)%is_sequence) THEN
                WRITE (string,20) i, REPEAT(CHAR(32),Nblanks), '- ',    &
     &                            TRIM(key), ': [', TRIM(value), ']'
              ELSE
                WRITE (string,30) i, REPEAT(CHAR(32),Nblanks), '- ',    &
     &                            TRIM(key), ': ', TRIM(value)
              END IF
            ELSE
              WRITE (string,40) i, REPEAT(CHAR(32),Nblanks), '- ',      &
     &                          TRIM(key), ':'
            END IF
          ELSE
            IF (Lstr.gt.0) THEN
              IF (self%list(i)%is_sequence) THEN
                WRITE (string,30) i, REPEAT(CHAR(32),Nblanks),          &
     &                            TRIM(key), ': [', TRIM(value), ']'
              ELSE
                WRITE (string,40) i, REPEAT(CHAR(32),Nblanks),          &
     &                            TRIM(key), ': ', TRIM(value)
              END IF
            ELSE
              WRITE (string,50) i, REPEAT(CHAR(32),Nblanks),            &
     &                          TRIM(key), ':'
            END IF
          END IF
          WRITE (yaml_stdout,60) TRIM(string)
        END DO
      END IF
!
  10  FORMAT (/,"YAML Tree Dictinary, file: '",a,"'",/,                 &
     &          '==========================',/)
  20  FORMAT ('L=',i4.4,1x,'% ',6a)
  30  FORMAT ('L=',i4.4,1x,'% ',5a)
  40  FORMAT ('L=',i4.4,1x,'% ',4a)
  50  FORMAT ('L=',i4.4,1x,'% ',3a)
  60  FORMAT (a)
!
      RETURN
      END SUBROUTINE yaml_tree_dump
!
      FUNCTION yaml_tree_extract (self, keystring, S) RESULT (status)
!
!***********************************************************************
!                                                                      !
!  It extracts YAML value(s) from processing key-string. The key       !
!  string may be a set of keys aggregated with a period, CHAR(46).     !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     self       YAML tree dictionary object (CLASS yaml_tree)         !
!     keystring  Aggregated YAML keys (string)                         !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     S          Separated YAML key/pair value (TYPE yaml_extract)     !
!     nkeys      Number of extracted keys (integer)                    !
!                                                                      !
!***********************************************************************
!
!  Imported variable declarations.
!
      CLASS (yaml_tree), intent(in) :: self
      character (len=*), intent(in) :: keystring
      TYPE (yaml_extract), allocatable, intent(inout) :: S(:)
!
!  Local variable declarations.
!
      TYPE (yaml_extract), allocatable :: K(:)
!
      logical :: BlockFlow
!
      integer :: i, ib, ic, ie, ik, ipair, is, j, li, pID
      integer :: Lstr, LenStr, nkeys, npairs, nvalues
      integer :: icomma, idot
      integer :: status
!
      integer, allocatable :: P(:)                   ! pair index
!
      character (len=:), allocatable :: Kstring      ! key string
      character (len=:), allocatable :: Vstring      ! value string
!
      character (len=*), parameter :: MyFile =                          &
     &   __FILE__//", yaml_tree_extract"
!
!-----------------------------------------------------------------------
!  Extract YAML key(s) from key-string.
!-----------------------------------------------------------------------
!
      status=Noerr
!
!  Count the number keys in key-string, separated by period.
!
      Lstr=LEN_TRIM(keystring)
      IF (yaml_Error(yaml_AssignString(Kstring,                         &
     &                                 keystring, LenStr),              &
     &               NoErr, __LINE__, MyFile)) RETURN
!
      nkeys=yaml_CountKeys(Kstring, CHAR(46))
!
!  Allocate key structure.
!
      ALLOCATE ( K(nkeys) )
!
!  Extract keys.
!
      is=1
      DO i=1,nkeys
        idot=INDEX(Kstring,CHAR(46),BACK=.FALSE.)
        ie=idot
        IF (idot.eq.0) THEN
          ie=LEN_TRIM(Kstring)
        ELSE
          ie=ie-1
        END IF
        IF (yaml_Error(yaml_AssignString(K(i)%value,                    &
     &                                   Kstring(is:ie), LenStr),       &
     &                 NoErr, __LINE__, MyFile)) RETURN
        IF (idot.gt.0) Kstring(is:ie+1) = REPEAT(CHAR(32), ie-is+2)
        Kstring=TRIM(ADJUSTL(Kstring))
      END DO
!
!-----------------------------------------------------------------------
!  Determine the number of YAML tree pairs to process and assoicated
!  list array element.
!-----------------------------------------------------------------------
!
!  Check if key-string is a blocking list where any of the members has
!  a leading hyphen as indentation. If blocking, all the pairs with the
!  same key-string are extracted.
!
      BlockFlow=.FALSE.          ! Switch to process all block membert
      ib=0                       ! block list(s) counter
      ic=0                       ! key/value pair counter
      ik=1                       ! key counter to avoid double DO-loops
!
      DO i=1,self%Npairs
        Lstr=LEN_TRIM(self%list(i)%key)
        IF ((self%list(i)%key).eq.(K(ik)%value)) THEN

          IF (yaml_Master.and.LdebugYAML) THEN
            PRINT '(2(a,i0,2a))', 'key ',ik,' = ', TRIM(K(ik)%value),   &
     &                            ', YAML list ',i,' = ',               &
     &                            TRIM(self%list(i)%key)
          END IF

          pID=self%list(i)%parent_id
          IF (self%list(i)%is_block.or.self%list(pID)%is_block) THEN
            ib=ib+1                             ! advance block counter
          END IF
          IF (ik.eq.nkeys) THEN
            ic=ic+1                             ! advance pair counter
            IF (ib.eq.0) THEN
              li=i                              ! list index to extract
              EXIT                              ! no blocking list found
            ELSE
              BlockFlow=.TRUE.
            END IF
          ELSE
            ik=ik+1                             ! advance key counter
          END IF
        END IF
        IF (BlockFlow.and.(self%list(i)%left_padding.eq.0)) THEN
          EXIT                                  ! processed all blocks
        END IF
      END DO
      npairs=ic                                 ! pairs to process
!
!  Allocate pair index array, P.
!
      IF (npairs.ne.0) THEN
        IF (.not.ALLOCATED(P)) ALLOCATE ( P(npairs) )
      ELSE
        yaml_ErrFlag=7
        status=yaml_ErrFlag
        IF (yaml_Error(yaml_ErrFlag, NoErr, __LINE__, MyFile)) THEN
          IF (yaml_Master) WRITE (yaml_stdout,10) keystring,            &
     &                                            self%filename
          RETURN
        END IF
      END IF
!
!  Set pair element(s) to extract.
!
      IF (BlockFlow) THEN
        ic=0
        ik=1
        DO i=1,self%Npairs
          IF ((self%list(i)%key).eq.(K(ik)%value)) THEN
            IF (ik.eq.nkeys) THEN
              ic=ic+1
              P(ic)=i                           ! multiple pair index
            ELSE
              ik=ik+1
            END IF
          END IF
          IF ((ic.gt.0).and.(self%list(i)%left_padding.eq.0)) THEN
            EXIT                                ! processed all blocks
          END IF
        END DO
      ELSE
        P(1)=li                                 ! single pair index
      END IF
!
!-----------------------------------------------------------------------
!  Extract pair(s) value(s).
!-----------------------------------------------------------------------
!
      DO i=1,npairs
        ipair=P(i)
!
!  Get key-string value.
!
        IF (yaml_Error(yaml_AssignString(Vstring,                       &
     &                                   self%list(ipair)%value,        &
     &                                   LenStr),                       &
     &                 NoErr, __LINE__, MyFile)) RETURN
!
!-----------------------------------------------------------------------
!  Extract/load keys-tring value(s). In a sequence, values are separated
!  by commas.
!-----------------------------------------------------------------------
!
!  Extract pair from a blocking list. Currently, scalar value are
!  allowed. Nested blocking is not supported since it is complicated
!  to process. For example, the following nested blocking lists are
!  not supported.
!
!  branch:
!    - blocklist1: value
!      blocklist1_key1: value
!        - blocklist1A: value                          ! not supported
!
        IF (BlockFlow) THEN
!
!  Process a block list pair with vector values (flow sequence):
!  S(i)%vector(1)%value  to  S(i)%vector(nvalues)%value
!
          IF (self%list(ipair)%is_sequence) THEN
            Lstr=LEN_TRIM(Vstring)
            nvalues=COUNT((/(Vstring(j:j), j=1,Lstr)/) == CHAR(44)) + 1
!
            IF (.not.ALLOCATED(S)) THEN
              ALLOCATE ( S(npairs) )                   ! main structure
            END IF
            ALLOCATE ( S(i)%vector(nvalues) )          ! sub-structure
            S(i)%has_vector=.TRUE.
!
            is=1
            DO j=1,nvalues
              icomma=INDEX(Vstring,CHAR(44),BACK=.FALSE.)
              ie=icomma
              IF (icomma.eq.0) THEN
                ie=LEN_TRIM(Vstring)
              ELSE
                ie=ie-1
              END IF
              IF (yaml_Error(yaml_AssignString(S(i)%vector(j)%value,    &
     &                                         Vstring(is:ie),          &
     &                                         LenStr),                 &
     &                       NoErr, __LINE__, MyFile)) RETURN

              IF (yaml_Master.and.LdebugYAML) THEN
                PRINT '(3a,2(i0,a),a)', 'keystring = ',TRIM(keystring), &
                                     ', S(', i, ')%vector(', j, ') = ', &
     &                               TRIM(S(i)%vector(j)%value)
              END IF

              IF (icomma.gt.0) Vstring(is:ie+1)=REPEAT(CHAR(32),ie-is+2)
              Vstring=TRIM(ADJUSTL(Vstring))
            END DO
!
!  Process a block list pair with a scalar value, S(i)%value
!
          ELSE
!
            IF (.not.ALLOCATED(S)) THEN
              ALLOCATE ( S(npairs) )
            END IF
            S(i)%has_vector=.FALSE.
!
            IF (yaml_Error(yaml_AssignString(S(i)%value,                &
     &                                       Vstring, LenStr),          &
     &                     NoErr, __LINE__, MyFile)) RETURN

            IF (yaml_Master.and.LdebugYAML) THEN
              PRINT '(a,i0,4a)', 'keystring ',i,' = ', TRIM(keystring), &
     &                           ', value = ', TRIM(S(i)%value)
            END IF
          END IF
!
!  Otherwise, process a non-block list.
!
        ELSE
!
!         Process a flow sequence:  S(1)%value to S(nvalues)%value.
!
          IF (self%list(ipair)%is_sequence) THEN
            Lstr=LEN_TRIM(Vstring)
            nvalues=COUNT((/(Vstring(j:j), j=1,Lstr)/) == CHAR(44)) + 1
!
            ALLOCATE ( S(nvalues) )
            S(i)%has_vector=.FALSE.
!
            is=1
            DO j=1,nvalues
              icomma=INDEX(Vstring,CHAR(44),BACK=.FALSE.)
              ie=icomma
              IF (icomma.eq.0) THEN
                ie=LEN_TRIM(Vstring)
              ELSE
                ie=ie-1
              END IF
              IF (yaml_Error(yaml_AssignString(S(j)%value,              &
     &                                         Vstring(is:ie),          &
     &                                         LenStr),                 &
     &                       NoErr, __LINE__, MyFile)) RETURN
              IF (icomma.gt.0) Vstring(is:ie+1)=REPEAT(CHAR(32),ie-is+2)
              Vstring=TRIM(ADJUSTL(Vstring))

              IF (yaml_Master.and.LdebugYAML) THEN
                PRINT '(a,i0,4a)', 'keystring ',j,' = ',                &
     &                             TRIM(keystring),                     &
     &                             ', value = ', TRIM(S(j)%value)
              END IF
            END DO
!
!         Process a single scalar value, S(1)%value.
!
          ELSE
!
            ALLOCATE ( S(1) )
            S(1)%has_vector=.FALSE.
!
            IF (yaml_Error(yaml_AssignString(S(1)%value,                &
     &                                       Vstring, LenStr),          &
     &                     NoErr, __LINE__, MyFile)) RETURN

            IF (yaml_Master.and.LdebugYAML) THEN
               PRINT '(4a)', 'keystring = ', TRIM(keystring),           &
     &                       ', value = ', TRIM(S(1)%value)
            END IF
          END IF
        END IF
      END DO
!
!-----------------------------------------------------------------------
!  Deallocate local variables.
!-----------------------------------------------------------------------
!
      IF (ALLOCATED(K))       DEALLOCATE (K)
      IF (ALLOCATED(P))       DEALLOCATE (P)
      IF (ALLOCATED(Kstring)) DEALLOCATE (Kstring)
      IF (ALLOCATED(Vstring)) DEALLOCATE (Vstring)
!
  10  FORMAT (/," YAML_TREE_EXTRACT - Cannot find key-string: '",a,     &
     &        "'",/,21x,'File: ',a)
  20  FORMAT (/," YAML_TREE_EXTRACT - Not supported key-string: '",a,   &
     &        "'",/,21x,'nested sub-blocking in a leading blocking ',   &
     &        'list',/,21x,'File: ',a)
!
      RETURN
      END FUNCTION yaml_tree_extract
!
      SUBROUTINE yaml_tree_fill (self)
!
!***********************************************************************
!                                                                      !
!  It reads YAML file and fills structure with the key/value pairs.    !
!  Both the key and value pairs are strings. The numercal convertions  !
!  are done elsewhere when needed.                                     !
!                                                                      !
!  Arguments:                                                          !
!                                                                      !
!     self       YAML tree dictionary object (CLASS yaml_tree)         !
!                                                                      !
!***********************************************************************
!
!  Imported variable declarations.
!
      CLASS (yaml_tree), intent(inout) :: self
!
!  Local variable declarations.
!
      logical :: Lswitch(Ldim)
!
      integer :: Nblanks, LenStr, left_padding, new_parent, old_parent
      integer :: i, ibranch, icount, ic_alias, ic_anchor, status
!
      character (len=Lkey) :: anchor, key
      character (len=Lmax) :: line, value
!
      character (len=*), parameter :: MyFile =                          &
     &   __FILE__//", yaml_tree_fill"
!
!-----------------------------------------------------------------------
!  Read and populate YAML structure.
!-----------------------------------------------------------------------
!
      ibranch=0
      icount=0
      ic_alias=0
      ic_anchor=0
      new_parent=0
      old_parent=0
!
      YAML_LINE : DO WHILE (.TRUE.)
!
        status=self%read_line(Nblanks, line, key, value, anchor,        &
     &                        Lswitch)
!
        SELECT CASE (status)
          CASE (-1)       ! error condition during reading
            yaml_ErrFlag=4
            IF (yaml_Error(yaml_ErrFlag, NoErr, __LINE__, MyFile)) THEN
              IF (yaml_Master) WRITE (yaml_stdout,10) self%filename,    &
     &                                                TRIM(line)
              RETURN
            END IF
          CASE (0)        ! end-of-file
            EXIT
          CASE (1)        ! blank or comment line, move to the next line
            CYCLE
          CASE (2)        ! processed viable line
            icount=icount+1
        END SELECT
!
!  Determine structure indices according to the nested levels counters.
!
        IF (Nblanks.eq.0) THEN
          ibranch=ibranch+1
          new_parent=icount
          old_parent=icount
        END IF
!
!  Load YAML pair switch identifiers.
!
        self%list(icount)%has_alias   = Lswitch(1)
        self%list(icount)%has_anchor  = Lswitch(2)
        self%list(icount)%is_block    = Lswitch(3)
        self%list(icount)%is_sequence = Lswitch(4)
        self%list(icount)%is_logical  = Lswitch(5)
        self%list(icount)%is_integer  = Lswitch(6)
        self%list(icount)%is_real     = Lswitch(7)
        self%list(icount)%is_string   = Lswitch(8)
!
        IF (Lswitch(1)) ic_alias=ic_alias+1
        IF (Lswitch(2)) ic_anchor=ic_anchor+1
!
!  Set left-padding indentation level: 0, 1, 2, ...
!
        left_padding=Nblanks/self%indent
        self%list(icount)%left_padding=left_padding
!
!  Load key/value ID and parent ID.
!
        IF (Nblanks.gt.0) THEN
          IF (left_padding.gt.self%list(icount-1)%left_padding) THEN
            new_parent=old_parent
            old_parent=icount
          END IF
        END IF
        self%list(icount)%id=icount
        self%list(icount)%parent_id=new_parent
!
!  Allocate and load line, key, and value strings. If applicable, loal
!  anchor value. Notice that it is possible to have keyword without
!  value in the main branches.
!
        IF (yaml_Error(yaml_AssignString(self%list(icount)%line,        &
     &                                   line, LenStr),                 &
     &                 NoErr, __LINE__, MyFile)) RETURN
!
        IF (yaml_Error(yaml_AssignString(self%list(icount)%key,         &
     &                                   key, LenStr),                  &
     &                 NoErr, __LINE__, MyFile)) RETURN
!
        IF (LEN_TRIM(value).gt.0) THEN
          IF (yaml_Error(yaml_AssignString(self%list(icount)%value,     &
     &                                     value, LenStr),              &
     &                   NoErr, __LINE__, MyFile)) RETURN
        END IF
!
        IF (Lswitch(2).and.LEN_TRIM(anchor).gt.0) THEN
          IF (yaml_Error(yaml_AssignString(self%list(icount)%anchor,    &
     &                                     anchor, LenStr),             &
     &                   NoErr, __LINE__, MyFile)) RETURN
        END IF
!
      END DO YAML_LINE
!
!  Substitute 'Alias' with 'Anchor' values. Notice that Anchors and
!  Aliases let you identify an item with an 'anchor' in a YAML pair,
!  and then refer to that item with an alias later in the same tile.
!  It is done to avoid repetitions and errors.
!
      IF ((ic_alias.gt.0).and.(ic_anchor.gt.0)) THEN
        CALL self%fill_aliases (ic_alias, ic_anchor)
      END IF
!
!  Close YAML file.
!
      CLOSE (iunit)
!
  10  FORMAT (/,' YAML_TREE_FILL - Error while reading YAML file: ',    &
     &        a,/,18x,'LINE: ',a)
!
      RETURN
      END SUBROUTINE yaml_tree_fill
!
      SUBROUTINE yaml_tree_fill_aliases (self, Nalias, Nanchor)
!
!***********************************************************************
!                                                                      !
!  It replaces the Aliases items with its respective Anchors values in !
!  YAML dictionary.                                                    !
!                                                                      !
!  Arguments:                                                          !
!                                                                      !
!     self       YAML tree dictionary object (CLASS yaml_tree)         !
!     Nalias     Number of Aliases items in YAML file (integer)        !
!     Nanchors   Number of Anchors in YAML file (integer)              !
!                                                                      !
!***********************************************************************
!
!  Imported variable declarations.
!
      CLASS (yaml_tree), intent(inout) :: self
      integer, intent(in) :: Nalias
      integer, intent(in) :: Nanchor
!
!  Local variable declarations.
!
      logical :: Lswitch(Ldim)
!
      integer :: Ialias, LenStr, i, j, ic
!
      character (len=Lkey), dimension(Nanchor) :: AnchorKey, AnchorVal
      character (len=Lmax) :: AliasVal
!
      character (len=*), parameter :: MyFile =                          &
     &   __FILE__//", yaml_tree_fill_aliases"
!
!-----------------------------------------------------------------------
!  Extract Anchors keyword and value.
!-----------------------------------------------------------------------
!
      ic=0
      DO i=1,self%Npairs
        IF (self%list(i)%has_anchor) THEN
          ic=ic+1
          AnchorKey(ic)=self%list(i)%anchor
          AnchorVal(ic)=self%list(i)%value
        END IF
      END DO
!
!-----------------------------------------------------------------------
!  Replace Aliases with associate Anchors values.
!-----------------------------------------------------------------------
!
      DO j=1,self%Npairs
!
!  Get Aliases keyword.
!
        IF (self%list(j)%has_alias) THEN
          AliasVal=self%list(j)%value
          Ialias=INDEX(AliasVal,CHAR(42),BACK=.FALSE.)    ! alias '*'
          IF (Ialias.gt.0) THEN
            AliasVal(Ialias:Ialias)=CHAR(32)              ! blank
            AliasVal=TRIM(ADJUSTL(AliasVal))
          END IF
!
!  Replace Aliases with anchor values and update value type.
!
          DO i=1,ic
            IF (TRIM(AliasVal).eq.TRIM(AnchorKey(i))) THEN
              DEALLOCATE (self%list(j)%value)
              IF (yaml_Error(yaml_AssignString(self%list(j)%value,      &
     &                                         TRIM(AnchorVal(i)),      &
     &                                         LenStr),                 &
     &                       NoErr, __LINE__, MyFile)) RETURN
!
              Lswitch=.FALSE.
              CALL yaml_ValueType (self%list(j)%value, Lswitch)
              self%list(j)%is_logical = Lswitch(5)
              self%list(j)%is_integer = Lswitch(6)
              self%list(j)%is_real    = Lswitch(7)
              self%list(j)%is_string  = Lswitch(8)
            END IF
          END DO
        END IF
      END DO
!
      RETURN
      END SUBROUTINE yaml_tree_fill_aliases
!
      FUNCTION yaml_tree_has (self, keystring) RESULT (FoundIt)
!
!***********************************************************************
!                                                                      !
!  It checks if YAML dictionary has requested key string.              !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     self       YAML tree dictionary object (CLASS yaml_tree)         !
!                                                                      !
!     keystring  Requested YAML key-string (string)                    !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     FoundIt    Switch indicating if the key string was found or not  !
!                  in the YAML dictionary.                             !
!                                                                      !
!***********************************************************************
!
!  Imported variable declarations.
!
      CLASS (yaml_tree), intent(inout) :: self
      character (len=*), intent(in   ) :: keystring
!
!  Local variable declarations.
!
      TYPE (yaml_extract), allocatable :: K(:)
!
      logical :: foundit
!
      integer :: Lstr, LenStr
      integer :: i, idot, ie, is, j, nkeys
!
      character (len=:), allocatable :: Kstring

      character (len=*), parameter :: MyFile =                          &
     &   __FILE__//", yaml_tree_has"
!
!-----------------------------------------------------------------------
!  Check if requested key-string is available in YAML dictionary
!-----------------------------------------------------------------------
!
      FoundIt=.FALSE.
!
!  Count number of keys in key-string, separated by period, CHAR(46).
!
      Lstr=LEN_TRIM(keystring)
      IF (yaml_Error(yaml_AssignString(Kstring,                         &
     &                                 keystring, LenStr),              &
     &               NoErr, __LINE__, MyFile)) RETURN
!
      nkeys=yaml_CountKeys(Kstring, CHAR(46))
!
!  Check single key.
!
      IF (nkeys.eq.1) THEN
        DO i=1,SIZE(self%list)
          IF (self%list(i)%key.eq.Kstring) THEN
            FoundIt=.TRUE.
            EXIT
          END IF
        END DO
        DEALLOCATE (Kstring)
        RETURN
!
!  Otherwise, check for compound key separated by period.
!
      ELSE
        ALLOCATE ( K(nkeys) )
!
!  Extract keys.
!
        is=1
        DO j=1,nkeys
          idot=INDEX(Kstring,CHAR(46),BACK=.FALSE.)
          ie=idot
          IF (idot.eq.0) THEN
            ie=LEN_TRIM(Kstring)
          ELSE
            ie=ie-1
          END IF
          IF (yaml_Error(yaml_AssignString(K(j)%value,                  &
     &                                     Kstring(is:ie), LenStr),     &
     &                   NoErr, __LINE__, MyFile)) RETURN
          IF (idot.gt.0) Kstring(is:ie+1) = REPEAT(CHAR(32), ie-is+2)
          Kstring=TRIM(ADJUSTL(Kstring))
        END DO
!
!  Check for compound key: val1.val2 ...
!  (It fails if any of the compound keys is not found)
!
        DO j=1,nkeys
          FoundIt=.FALSE.
          DO i=1,SIZE(self%list)
            IF (self%list(i)%key.eq.K(j)%value) THEN
              FoundIt=.TRUE.
              EXIT
            END IF
          END DO
          IF (.not.FoundIt) EXIT
        END DO
        DEALLOCATE (K, Kstring)
        RETURN
      END IF
!
      END FUNCTION yaml_tree_has
!
      FUNCTION yaml_tree_read_line (self, Nblanks, line, key, value,    &
     &                              anchor, Lswitch)                    &
     &                      RESULT (status)
!
!***********************************************************************
!                                                                      !
!  It reads a single text line from current YAML file. It uses the     !
!  ASCII character set extensively.                                    !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     self       YAML tree dictionary object (CLASS yaml_tree)         !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     Nblanks    YAML leading blanks identation (integer)              !
!                                                                      !
!     line       Current YAML file line (string)                       !
!                                                                      !
!     key        YAML line keyword (string)                            !
!                                                                      !
!     value      YAML line keyword value (string)                      !
!                                                                      !
!     anchor     YAML value ancher keyword (string)                    !
!                                                                      !
!     Lswitch    YAML key/value switches (logical, vector)             !
!                                                                      !
!                  Lswitch(1) = T/F, value has an alias '*' token      !
!                  Lswitch(2) = T/F, value has an anchor '&' token     !
!                  Lswitch(3) = T/F, key/value starts a block '-'      !
!                  Lswitch(4) = T/F, value has sequence '[...]' tokens !
!                  Lswitch(5) = T/F, logical value(s)                  !
!                  Lswitch(6) = T/F, integer value(s)                  !
!                  Lswitch(7) = T/F, floating-point value(s)           !
!                  Lswitch(8) = T/F, string value(s)                   !
!                                                                      !
!***********************************************************************
!
!  Imported variable declarations.
!
      CLASS (yaml_tree), intent(inout) :: self
!
      logical, intent(out) :: Lswitch(:)
!
      integer, intent(out) :: Nblanks
!
      character (len=*), intent(out) :: line
      character (len=*), intent(out) :: key
      character (len=*), intent(out) :: value
      character (len=*), intent(out) :: anchor
!
!  Local variable declarations.
!
      logical :: Lbracket, Rbracket
!
      integer :: Ialias, Ianchor, Iblank, Icolon, Idash, Ihash, Ispace
      integer :: IbracketL, IbracketR
      integer :: Lstr, LstrNext, LstrVal, i, j, k
      integer :: status
!
      character (len=Lmax) :: linecopy, next_line

      character (len=*), parameter :: MyFile =                          &
     &   __FILE__//", yaml_tree_read_line"
!
!-----------------------------------------------------------------------
!  Read a single YAML file line
!-----------------------------------------------------------------------
!
!  Initialize.
!
      status=-1                 ! error condition
!
      DO i=1,LEN(key)
        key(i:i)=CHAR(32)
      END DO

      DO i=1,LEN(value)
        value(i:i)=CHAR(32)
      END DO

      DO i=1,LEN(anchor)
        anchor(i:i)=CHAR(32)
      END DO
!
      Lswitch=.FALSE.
      Lbracket=.FALSE.
      Rbracket=.FALSE.
!
!  Read in next YAML file line.
!
      READ (iunit,'(a)',ERR=10,END=20) line
!
!  Replace illegal control characters with a blank space CHAR(32)
!
      Lstr=LEN_TRIM(line)

      DO i=1,LEN_TRIM(line)
        j=ICHAR(line(i:i))
        IF (j.lt.32) THEN
          line(i:i)=CHAR(32)
        END IF
      END DO
      linecopy=TRIM(ADJUSTL(line))
!
!   Get length of "line". Remove comment after the KEYWORD, if  any.
!   In YAML, a comment starts with a hash (#), CHAR(35).
!
      IF ((Lstr.gt.0).and.(linecopy(1:1).ne.CHAR(35))) THEN
        Ihash=INDEX(line,CHAR(35),BACK=.FALSE.)
        IF (Ihash.gt.0) THEN
          Lstr=Ihash-1
          line=TRIM(line(1:Lstr))         ! remove trailing comment
          Lstr=LEN_TRIM(line)
        END IF
        status=2                          ! viable line
      ELSE
        status=1                          ! blank or commented line,
        RETURN                            ! move to the next line
      END IF
!
!  Find colon sign CHAR(58) and determine the number of leading blank
!  spaces. YAML uses indentations with one or more spaces to describe
!  nested collections (lists, mappings). YAML also uses dash plus space,
!  '- ', as enumerator of block lists.
!
!  It checks if the KEYWORD is a multi-word separated by space.
!
      Icolon=INDEX(line,CHAR(58),BACK=.FALSE.)
      IF (Icolon.eq.0) THEN
        status=-1
        yaml_ErrFlag=2
        IF (yaml_Master) THEN
          WRITE (yaml_stdout,30) TRIM(line)
        END IF
        IF (yaml_Error(yaml_ErrFlag, NoErr, __LINE__, MyFile)) RETURN
      END IF
!
      Idash =INDEX(line,CHAR(45),BACK=.FALSE.)
      IF ((Idash.gt.0).and.(Idash.lt.Icolon)) THEN
        Iblank=INDEX(line(1:Idash),CHAR(32),BACK=.TRUE.)
      ELSE
        Iblank=INDEX(line(1:Icolon),CHAR(32),BACK=.TRUE.)
        IF (Iblank.gt.0) THEN
          k=Iblank-1
          IF  ((65.le.ICHAR(line(1:1))).and.                            &
     &         (ICHAR(line(1:1)).le.122)) THEN   ! multi-word branch key
            Iblank=0
          ELSE IF ((65.le.ICHAR(line(k:k))).and.                        &
     &             (ICHAR(line(k:k)).le.122)) THEN ! multi-word pair key
            Iblank=INDEX(line(1:k),CHAR(32),BACK=.TRUE.)
          END IF
        END IF
      END IF
!
!  Set number of YAML line leading blacks.
!
      Nblanks=Iblank
!
!  Extract key and value pair and return.
!
      IF ((Idash.gt.0).and.(Idash.lt.Icolon)) THEN
        key=TRIM(ADJUSTL(line(Idash+1:Icolon-1)))
      ELSE
        key=TRIM(ADJUSTL(line(1:Icolon-1)))
      END IF
      value=TRIM(ADJUSTL(line(Icolon+1:Lstr)))
!
!  Check for special YAML tokens in value string and replace with blank.
!
      Ialias=INDEX(value,CHAR(42),BACK=.FALSE.)      ! alias '*'
      IF (Ialias.gt.0) THEN
        Lswitch(1)=.TRUE.
      END IF
!
      Ianchor=INDEX(value,CHAR(38),BACK=.FALSE.)     ! anchor '&'
      IF (Ianchor.gt.0) THEN
        Ispace=INDEX(value(Ianchor+1:),CHAR(32),BACK=.FALSE.)
        anchor=value(Ianchor+1:Ispace)               ! anchor value
        Lswitch(2)=.TRUE.
        value(Ianchor:Ispace)=REPEAT(CHAR(32),Ispace-Ianchor+1)
        value=TRIM(ADJUSTL(value))
      END IF
!
      IF ((Idash.gt.0).and.(Idash.lt.Icolon)) THEN
        Lswitch(3)=.TRUE.                            ! block pair '- '
      END IF
!
      IbracketL=INDEX(value,CHAR(91),BACK=.FALSE.)   ! left  bracket '['
      IF (IbracketL.gt.0) THEN
        Lbracket=.TRUE.
        value(IbracketL:IbracketL)=CHAR(32)
        value=TRIM(ADJUSTL(value))
      END IF
!
      IbracketR=INDEX(value,CHAR(93),BACK=.FALSE.)   ! right bracket ']'
      IF (IbracketR.gt.0) THEN
        Rbracket=.TRUE.
        value(IbracketR:IbracketR)=CHAR(32)
        value=TRIM(ADJUSTL(value))
      END IF
!
!-----------------------------------------------------------------------
!  If right square bracket is not found, the key values are broken into
!  multiple lines. Process the necessary lines.
!-----------------------------------------------------------------------
!
      IF (.not.Rbracket.and.Lbracket) THEN
        DO WHILE (.not.Rbracket)
          READ (iunit,'(a)',ERR=10,END=20) next_line
!
!  Replace illegal control characters with a blank space CHAR(32)
!
          DO i=1,LEN_TRIM(next_line)
            j=ICHAR(next_line(i:i))
            IF (j.lt.32) THEN
              next_line(i:i)=CHAR(32)
            END IF
          END DO
          next_line=TRIM(ADJUSTL(next_line))
!
!  If applicable, remove trailing comments starting with a hash (#),
!  CHAR(35).
!
          Ihash=INDEX(next_line,CHAR(35),BACK=.FALSE.)
          LstrNext=LEN_TRIM(next_line)
          IF ((LstrNext.gt.0).and.(Ihash.gt.0)) THEN
            LstrNext=Ihash-1
            next_line=TRIM(next_line(1:LstrNext))
            LstrNext=LEN_TRIM(next_line)
          END IF
!
!  Aggregate new 'next_line' to previous 'line' and 'value'.
!
          Lstr=LEN_TRIM(line)
          LstrVal=LEN_TRIM(value)
          line=line(1:Lstr)//CHAR(32)//next_line(1:LstrNext)
          value=value(1:LstrVal)//CHAR(32)//next_line(1:LstrNext)
!
          IbracketR=INDEX(value,CHAR(93),BACK=.FALSE.)
          IF (IbracketR.gt.0) THEN
            Rbracket=.TRUE.
            value(IbracketR:IbracketR)=CHAR(32)
            value=TRIM(ADJUSTL(value))
          END IF
        END DO
      END IF
      IF (Lbracket.and.Rbracket) Lswitch(4)=.TRUE.
!
!-----------------------------------------------------------------------
!  Determine value representation: logical (boolean), numerical
! (integer/reals), or string(s).  Others can be added as needed.
!-----------------------------------------------------------------------
!
      CALL yaml_ValueType (value, Lswitch)
!
      RETURN
!
!  Read flow control determines the status flag.
!
  10  status=-1     ! error during reading
  20  status=0      ! end-of-file encoutered
!
  30  FORMAT (/,' YAML_TREE_READ_LINE - Unable to find colon token ',   &
     &        'after key word',/,23x,'LINE: ',a)

      END FUNCTION yaml_tree_read_line
!
      FUNCTION yaml_AssignString (OutString, InpString, Lstr)           &
     &                    RESULT (ErrFlag)
!
!=======================================================================
!                                                                      !
!  It assigns allocatable strings. It allocates/reallocates output     !
!  string variable.                                                    !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     InpString    String to be assigned (character)                   !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     OutString    Assigned allocatable string (character)             !
!     Lstr         Length allocated string (integer)                   !
!     ErrFlag      Error flag (integer)                                !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      character (len=:), allocatable, intent(inout) :: OutString
      character (len=*), intent(in) :: InpString
      integer, intent(out) :: Lstr
!
!  Local variable declarations.
!
      integer :: ErrFlag
!
!-----------------------------------------------------------------------
!  Allocate output string to the size of input string.
!-----------------------------------------------------------------------
!
      ErrFlag = -1
!
      Lstr=LEN_TRIM(InpString)
      IF (.not.allocated(OutString)) THEN
        allocate ( character(LEN=Lstr) :: OutString, STAT=ErrFlag)
      ELSE
        deallocate (OutString)
        allocate ( character(LEN=Lstr) :: OutString, STAT=ErrFlag)
      END IF
!
!  Assign requested value.
!
      OutString = InpString
!
      RETURN
      END FUNCTION yaml_AssignString
!
      FUNCTION yaml_CountKeys (string, token) RESULT (nkeys)
!
!=======================================================================
!                                                                      !
!  It counts the number of concatenated key separated by the specified !
!  token during processing extraction from YAML object.  The same task !
!  can be done elegantly as:                                           !
!                                                                      !
!  nkeys=COUNT((/ (string(i:i), i=1,Lstr) /) == token) + 1             !
!                                                                      !
!  But compilier like 'gfortran' cannot handle such abstraction.       !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     string     Aggregated YAML keys (string)                         !
!     token      Key separator (string)                                !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     nkeys      Number of aggregated keys (integer)                   !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      character (len=*), intent(in) :: string
      character (len=*), intent(in) :: token
!
!  Local variable declarations.
!
      integer :: nkeys
      integer :: Lstr, i
!
!-----------------------------------------------------------------------
!  Count number of concatenated keys in input string.
!-----------------------------------------------------------------------
!
      nkeys=1
      Lstr=LEN_TRIM(string)
      DO i=1,Lstr
        IF (string(i:i).eq.token) THEN
          nkeys=nkeys+1
        END IF
      END DO
!
      RETURN
      END FUNCTION yaml_CountKeys
!
      FUNCTION yaml_Error (flag, NoErr, line, routine) RESULT (foundit)
!
!=======================================================================
!                                                                      !
!  It checks YAML execution flag against no-error code and issue a     !
!  message if they are not equal.                                      !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     flag         YAML execution flag (integer)                       !
!     NoErr        No Error code (integer)                             !
!     line         Calling model routine line (integer)                !
!     routine      Calling model routine (string)                      !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     foundit      The value of the result is TRUE/FALSE if the        !
!                    execution flag is in error.                       !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer, intent(in) :: flag, NoErr, line

      character (len=*), intent(in) :: routine
!
!  Local variable declarations.
!
      logical :: foundit
!
!-----------------------------------------------------------------------
!  Scan array for requested string.
!-----------------------------------------------------------------------
!
      foundit=.FALSE.
      IF (flag.ne.NoErr) THEN
        foundit=.TRUE.
        IF (yaml_Master) THEN
          WRITE (yaml_stdout,10) flag, line, TRIM(routine)
  10      FORMAT (' Found Error: ', i2.2, t20, 'Line: ', i0,            &
     &            t35, 'Source: ', a)
        END IF
        FLUSH (yaml_stdout)
      END IF

      RETURN
      END FUNCTION yaml_Error
!
      FUNCTION yaml_Get_i_struc (self, keystring, V) RESULT (status)
!
!***********************************************************************
!                                                                      !
!  It loads a vector set of integers in YAML block-list structure,     !
!  V(1:Nitems)%vector(1:Nvalues).  If the dummy argument V structure   !
!  is allocated, it deallocates/allocates the required Nitems and      !
!  Nvalues dimensions.                                                 !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     self       YAML tree dictionary object (CLASS yaml_tree)         !
!     keystring  YAML tree aggregated keys (string)                    !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     V          Vector of integers in block list (TYPE yaml_Ivec)     !
!     status     Error code (integer)                                  !
!                                                                      !
!***********************************************************************
!
!  Imported variable declarations.
!
      CLASS (yaml_tree), intent(in) :: self
      character (len=*), intent(in) :: keystring
      TYPE (yaml_Ivec), allocatable, intent(out) :: V(:)
!
!  Local variable declarations.
!
      TYPE (yaml_extract), allocatable :: S(:)
!
      integer :: LenStr, Nitems, Nvalues, i, n
      integer :: status
!
      character (len=Lmax) :: msg

      character (len=*), parameter :: MyFile =                          &
     &   __FILE__//", yaml_Get_i_struc"
!
!-----------------------------------------------------------------------
!  Extract requested key-string values.
!-----------------------------------------------------------------------
!
      status=NoErr
!
      IF (yaml_Error(self%extract(keystring, S),                        &
     &               NoErr, __LINE__, MyFile)) RETURN
!
!  Allocate output structure.
!
      Nitems=SIZE(S, DIM=1)
      IF (ALLOCATED(V)) DEALLOCATE (V)
      ALLOCATE ( V(Nitems) )
!
!  Convert string vector values to integers.
!
      DO n=1,Nitems
        IF (S(n)%has_vector) THEN
          Nvalues=SIZE(S(n)%vector)
          ALLOCATE ( V(n)%vector(Nvalues) )
          DO i=1,Nvalues
            READ (S(n)%vector(i)%value, * ,IOSTAT=status, IOMSG=msg)    &
     &           V(n)%vector(i)
            IF (yaml_Error(status, NoErr, __LINE__, MyFile)) THEN
              yaml_ErrFlag=5
              IF (yaml_Master) WRITE (yaml_stdout,10) TRIM(keystring),  &
     &                                            S(n)%vector(i)%value, &
     &                                                TRIM(msg)
              RETURN
            END IF
          END DO
        END IF
      END DO
!
!  Deallocate local extraction structure.
!
      IF (ALLOCATED(S)) DEALLOCATE (S)
!
  10  FORMAT (/,' YAML_GET_I_STRUC - Error while converting string to', &
     &        ' integer, key: ',a,', value = ',a,/,20x,'ErrMsg: ',a)
!
      RETURN
      END FUNCTION yaml_Get_i_struc
!
      FUNCTION yaml_Get_l_struc (self, keystring, V) RESULT (status)
!
!***********************************************************************
!                                                                      !
!  It loads a vector set of logicals in YAML block-list structure,     !
!  V(1:Nitems)%vector(1:Nvalues).  If the dummy argument V structure   !
!  is allocated, it deallocates/allocates the required Nitems and      !
!  Nvalues dimensions.                                                 !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     self       YAML tree dictionary object (CLASS yaml_tree)         !
!     keystring  YAML tree aggregated keys (string)                    !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     V          Vector of logicals in block list (TYPE yaml_Lvec)     !
!     status     Error code (integer)                                  !
!                                                                      !
!***********************************************************************
!
!  Imported variable declarations.
!
      CLASS (yaml_tree), intent(in) :: self
      character (len=*), intent(in) :: keystring
      TYPE (yaml_Lvec), allocatable, intent(out) :: V(:)
!
!  Local variable declarations.
!
      TYPE (yaml_extract), allocatable :: S(:)
!
      integer :: Nitems, Nvalues, n, i
      integer :: status
!
      character (len=*), parameter :: MyFile =                          &
     &   __FILE__//", yaml_Get_l_struc"
!
!-----------------------------------------------------------------------
!  Extract requested key-string values.
!-----------------------------------------------------------------------
!
      status=NoErr
!
      IF (yaml_Error(self%extract(keystring, S),                        &
     &               NoErr, __LINE__, MyFile)) RETURN
!
!  Allocate output structure.
!
      Nitems=SIZE(S, DIM=1)
      IF (ALLOCATED(V)) DEALLOCATE (V)
      ALLOCATE ( V(Nitems) )
!
!  Convert string vector values to logicals.
!
      DO n=1,Nitems
        IF (S(n)%has_vector) THEN
          Nvalues=SIZE(S(n)%vector)
          ALLOCATE ( V(n)%vector(Nvalues) )
          DO i=1,Nvalues
            IF (yaml_UpperCase(S(n)%vector(i)%value(1:1)).eq.'T') THEN
              V(n)%vector(i)=.TRUE.
            ELSE
              V(n)%vector(i)=.FALSE.
            END IF
          END DO
        END IF
      END DO
!
!  Deallocate local extraction structure.
!
      IF (ALLOCATED(S)) DEALLOCATE (S)
!
      RETURN
      END FUNCTION yaml_Get_l_struc
!
      FUNCTION yaml_Get_r_struc (self, keystring, V) RESULT (status)
!
!***********************************************************************
!                                                                      !
!  It loads a vector set of real values in YAML block-list structure,  !
!  V(1:Nitems)%vector(1:Nvalues).  If the dummy argument V structure   !
!  is allocated, it deallocates/allocates the required Nitems and      !
!  Nvalues dimensions.                                                 !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     self       YAML tree dictionary object (CLASS yaml_tree)         !
!     keystring  YAML tree aggregated keys (string)                    !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     V          Vector of reals in block list (TYPE yaml_Rvec)        !
!     status     Error code (integer)                                  !
!                                                                      !
!***********************************************************************
!
!  Imported variable declarations.
!
      CLASS (yaml_tree), intent(in) :: self
      character (len=*), intent(in) :: keystring
      TYPE (yaml_Rvec), allocatable, intent(out) :: V(:)
!
!  Local variable declarations.
!
      TYPE (yaml_extract), allocatable :: S(:)
!
      integer :: Nitems, Nvalues, i, n
      integer :: status
!
      character (len=Lmax) :: msg

      character (len=*), parameter :: MyFile =                          &
     &   __FILE__//", yaml_Get_r_struc"
!
!-----------------------------------------------------------------------
!  Extract requested key-string values.
!-----------------------------------------------------------------------
!
      status=NoErr
!
      IF (yaml_Error(self%extract(keystring, S),                        &
     &               NoErr, __LINE__, MyFile)) RETURN
!
!  Allocate output structure.
!
      Nitems=SIZE(S, DIM=1)
      IF (ALLOCATED(V)) DEALLOCATE (V)
      ALLOCATE ( V(Nitems) )
!
!  Convert string vector values to floating-point.
!
      DO n=1,Nitems
        IF (S(n)%has_vector) THEN
          Nvalues=SIZE(S(n)%vector)
          ALLOCATE ( V(n)%vector(Nvalues) )
          DO i=1,Nvalues
            READ (S(n)%vector(i)%value, * ,IOSTAT=status, IOMSG=msg)    &
     &           V(n)%vector(i)
            IF (yaml_Error(status, NoErr, __LINE__, MyFile)) THEN
              yaml_ErrFlag=5
              IF (yaml_Master) WRITE (yaml_stdout,10) TRIM(keystring),  &
     &                                            S(n)%vector(i)%value, &
     &                                                TRIM(msg)
              RETURN
            END IF
          END DO
        ELSE
          READ (S(n)%value, * ,IOSTAT=status, IOMSG=msg)                &
     &         V(n)%vector(i)
          IF (yaml_Error(status, NoErr, __LINE__, MyFile)) THEN
            yaml_ErrFlag=5
            IF (yaml_Master) WRITE (yaml_stdout,10) TRIM(keystring),    &
     &                                              S(n)%value,         &
     &                                              TRIM(msg)
            RETURN
          END IF
        END IF
      END DO
!
!  Deallocate local extraction structure.
!
      IF (ALLOCATED(S)) DEALLOCATE (S)
!
  10  FORMAT (/,' YAML_GET_R_STRUC - Error while converting string to', &
     &        ' real, key: ',a,', value = ',a,/,20x,'ErrMsg: ',a)
!
      RETURN
      END FUNCTION yaml_Get_r_struc
!
      FUNCTION yaml_Get_s_struc (self, keystring, V) RESULT (status)
!
!***********************************************************************
!                                                                      !
!  It loads a vector set of strings in YAML block-list structure,      !
!  V(1:Nitems)%vector(1:Nvalues)%value.  If the dummy argument V       !
!  structure is allocated, it deallocates/allocates the required       !
!  Nitems and Nvalues dimensions.                                      !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     self       YAML tree dictionary object (CLASS yaml_tree)         !
!     keystring  YAML tree aggregated keys (string)                    !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     V          Vector of strings in block list (TYPE yaml_Svec)      !
!     status     Error code (integer)                                  !
!                                                                      !
!***********************************************************************
!
!  Imported variable declarations.
!
      CLASS (yaml_tree), intent(in) :: self
      character (len=*), intent(in) :: keystring
      TYPE (yaml_Svec), allocatable, intent(out) :: V(:)
!
!  Local variable declarations.
!
      TYPE (yaml_extract), allocatable :: S(:)
!
      integer :: Nitems, Nvalues, n, i
      integer :: status
!
      character (len=*), parameter :: MyFile =                          &
     &   __FILE__//", yaml_Get_s_struc"
!
!-----------------------------------------------------------------------
!  Extract requested key-string values.
!-----------------------------------------------------------------------
!
      status=NoErr
!
      IF (yaml_Error(self%extract(keystring, S),                        &
     &               NoErr, __LINE__, MyFile)) RETURN
!
!  Allocate output structure.
!
      Nitems=SIZE(S, DIM=1)
      IF (ALLOCATED(V)) DEALLOCATE (V)
      ALLOCATE ( V(Nitems) )
!
!  Load string vector values to output structure.
!
      DO n=1,Nitems
        IF (S(n)%has_vector) THEN
          Nvalues=SIZE(S(n)%vector)
          ALLOCATE ( V(n)%vector(Nvalues) )
          DO i=1,Nvalues
            V(n)%vector(i)%value=S(n)%vector(i)%value
          END DO
        ELSE
          V(n)%value=S(n)%value
        END IF
      END DO
!
!  Deallocate local extraction structure.
!
      IF (ALLOCATED(S)) DEALLOCATE (S)
!
      RETURN
      END FUNCTION yaml_Get_s_struc
!
      FUNCTION yaml_Get_ivar_0d (self, keystring, value) RESULT (status)
!
!***********************************************************************
!                                                                      !
!  It gets scalar integer data from YAML dictionary object.            !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     self       YAML tree dictionary object (CLASS yaml_tree)         !
!     keystring  YAML tree aggregated keys (string)                    !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     value      YAML value (integer; scalar)                          !
!     status     Error code (integer)                                  !
!                                                                      !
!***********************************************************************
!
!  Imported variable declarations.
!
      CLASS (yaml_tree), intent(in) :: self
      character (len=*), intent(in) :: keystring
      integer, intent(out) :: value
!
!  Local variable declarations.
!
      TYPE (yaml_extract), allocatable :: S(:)
!
      integer :: Nvalues, status
!
      character (len=Lmax) :: msg

      character (len=*), parameter :: MyFile =                          &
     &   __FILE__//", yaml_Get_ivar_0d"
!
!-----------------------------------------------------------------------
!  Extract requested key-string value.
!-----------------------------------------------------------------------
!
      status=NoErr
!
      IF (yaml_Error(self%extract(keystring, S),                        &
     &               NoErr, __LINE__, MyFile)) RETURN
!
!  Make sure that extracted value is a scalar.
!
      Nvalues=SIZE(S)
!
      IF (Nvalues.gt.1) THEN
        status=7
        yaml_ErrFlag=status
        IF (yaml_Error(status, NoErr, __LINE__, MyFile)) THEN
          IF (yaml_Master) WRITE (yaml_stdout,10) keystring, Nvalues,   &
     &                                            self%filename
          RETURN
        END IF
!
!  Convert string value to integer.
!
      ELSE
        READ (S(1)%value, *, IOSTAT=status, IOMSG=msg) value
        IF (yaml_Error(status, NoErr, __LINE__, MyFile)) THEN
          yaml_ErrFlag=5
          IF (yaml_Master) WRITE (yaml_stdout,20) TRIM(keystring),      &
     &                                            S(1)%value, TRIM(msg)
          RETURN
        END IF
      END IF
!
!  Deallocate.
!
      IF (ALLOCATED(S)) DEALLOCATE (S)
!
  10  FORMAT (/,' YAML_GET_IVAR_0D - Extracted value is not a scalar,', &
     &        ' key-string: ',a,/,20x,'Nvalues = ',i0,/,20x,'File: ',a)
  20  FORMAT (/,' YAML_GET_IVAR_0D - Error while converting string to', &
     &        ' integer, key: ',a,', value = ',a,/,20x,'ErrMsg: ',a)
!
      RETURN
      END FUNCTION yaml_Get_ivar_0d
!
      FUNCTION yaml_Get_ivar_1d (self, keystring, value) RESULT (status)
!
!***********************************************************************
!                                                                      !
!  It gets 1D integer data from YAML dictionary object.                !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     self       YAML tree dictionary object (CLASS yaml_tree)         !
!     keystring  YAML tree aggregated keys (string)                    !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     value      YAML value (integer; 1D-array)                        !
!     status     Error code (integer)                                  !
!                                                                      !
!***********************************************************************
!
!  Imported variable declarations.
!
      CLASS (yaml_tree), intent(in) :: self
      character (len=*), intent(in) :: keystring
      integer, intent(out) :: value(:)
!
!  Local variable declarations.
!
      TYPE (yaml_extract), allocatable :: S(:)
!
      integer :: Nvalues, i, status
!
      character (len=Lmax) :: msg

      character (len=*), parameter :: MyFile =                          &
     &   __FILE__//", yaml_Get_ivar_1d"
!
!-----------------------------------------------------------------------
!  Extract requested key-string values.
!-----------------------------------------------------------------------
!
      status=NoErr
!
      IF (yaml_Error(self%extract(keystring, S),                        &
     &               NoErr, __LINE__, MyFile)) RETURN
!
!  Make sure that extracted value is a scalar.
!
      Nvalues=SIZE(S, DIM=1)
!
      IF (SIZE(value, DIM=1).lt.Nvalues) THEN
        status=7
        yaml_ErrFlag=status
        IF (yaml_Error(status, NoErr, __LINE__, MyFile)) THEN
          IF (yaml_Master) WRITE (yaml_stdout,10) keystring, Nvalues,   &
     &                                            SIZE(value, DIM=1),   &
     &                                            self%filename
          RETURN
        END IF
      END IF
!
!  Convert string values to integers.
!
      DO i=1,Nvalues
        READ (S(i)%value, *, IOSTAT=status, IOMSG=msg) value(i)
        IF (yaml_Error(status, NoErr, __LINE__, MyFile)) THEN
          yaml_ErrFlag=5
          IF (yaml_Master) WRITE (yaml_stdout,20) TRIM(keystring),      &
     &                                            S(i)%value, TRIM(msg)
          RETURN
        END IF
      END DO
!
!  Deallocate.
!
      IF (ALLOCATED(S)) DEALLOCATE (S)
!
  10  FORMAT (/,' YAML_GET_IVAR_1D - Inconsistent number of values,',   &
     &        ' key-string: ',a,/,20x,'YAML size = ',i0,                &
              ', Variable size = ',i0,/,20x,'File: ',a)
  20  FORMAT (/,' YAML_GET_IVAR_1D - Error while converting string to', &
     &        ' integer, key: ',a,', value = ',a,/,20x,'ErrMsg: ',a)
!
      RETURN
      END FUNCTION yaml_Get_ivar_1d
!
      FUNCTION yaml_Get_lvar_0d (self, keystring, value) RESULT (status)
!
!***********************************************************************
!                                                                      !
!  It gets scalar logical data from YAML disctionary object.           !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     self       YAML tree dictionary object (CLASS yaml_tree)         !
!     keystring  YAML tree aggregated keys (string)                    !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     value      YAML value (logical; scalar)                          !
!                                                                      !
!***********************************************************************
!
!  Imported variable declarations.
!
      CLASS (yaml_tree), intent(in) :: self
      character (len=*), intent(in) :: keystring
      logical, intent(out) :: value
!
!  Local variable declarations.
!
      TYPE (yaml_extract), allocatable :: S(:)
!
      integer :: Nvalues, status
!
      character (len=*), parameter :: MyFile =                          &
     &   __FILE__//", yaml_Get_lvar_0d"
!
!-----------------------------------------------------------------------
!  Extract requested key-string value.
!-----------------------------------------------------------------------
!
      status=NoErr
!
      IF (yaml_Error(self%extract(keystring, S),                        &
     &               NoErr, __LINE__, MyFile)) RETURN
!
!  Make sure that extracted value is a scalar.
!
      Nvalues=SIZE(S)
!
      IF (Nvalues.gt.1) THEN
        status=7
        yaml_ErrFlag=status
        IF (yaml_Error(status, NoErr, __LINE__, MyFile)) THEN
          IF (yaml_Master) WRITE (yaml_stdout,10) keystring, Nvalues,   &
     &                                            self%filename
          RETURN
        END IF
!
!  Convert string value to logical.
!
      ELSE
        IF (yaml_UpperCase(S(1)%value(1:1)).eq.'T') THEN
          value=.TRUE.
        ELSE
          value=.FALSE.
        END IF
      END IF
!
!  Deallocate.
!
      IF (ALLOCATED(S)) DEALLOCATE (S)
!
  10  FORMAT (/,' YAML_GET_LVAR_0D - Extracted value is not a scalar,', &
     &        ' key-string: ',a,/,20x,'Nvalues = ',i0,/,20x,'File: ',a)
!
      RETURN
      END FUNCTION yaml_Get_lvar_0d
!
      FUNCTION yaml_Get_lvar_1d (self, keystring, value) RESULT (status)
!
!***********************************************************************
!                                                                      !
!  It gets 1D logical data from YAML dictionary object.                !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     self       YAML tree dictionary object (CLASS yaml_tree)         !
!     keystring  YAML tree aggregated keys (string)                    !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     value      YAML value (logical; 1D-array)                        !
!     status     Error code (integer)                                  !
!                                                                      !
!***********************************************************************
!
!  Imported variable declarations.
!
      CLASS (yaml_tree), intent(in) :: self
      character (len=*), intent(in) :: keystring
      logical, intent(out) :: value(:)
!
!  Local variable declarations.
!
      TYPE (yaml_extract), allocatable :: S(:)
!
      integer :: Nvalues, i, status
!
      character (len=*), parameter :: MyFile =                          &
     &   __FILE__//", yaml_Get_lvar_1d"
!
!-----------------------------------------------------------------------
!  Extract requested key-string values.
!-----------------------------------------------------------------------
!
      status=NoErr
!
      IF (yaml_Error(self%extract(keystring, S),                        &
     &               NoErr, __LINE__, MyFile)) RETURN
!
!  Make sure that extracted value is a scalar.
!
      Nvalues=SIZE(S, DIM=1)
!
      IF (SIZE(value, DIM=1).lt.Nvalues) THEN
        yaml_ErrFlag=7
        status=yaml_ErrFlag
        IF (yaml_Error(yaml_ErrFlag, NoErr, __LINE__, MyFile)) THEN
          IF (yaml_Master) WRITE (yaml_stdout,10) keystring, Nvalues,   &
     &                                            SIZE(value, DIM=1),   &
     &                                            self%filename
          RETURN
        END IF
      END IF
!
!  Convert string values to logicals.
!
      DO i=1,Nvalues
        IF (yaml_UpperCase(S(i)%value(1:1)).eq.'T') THEN
          value(i)=.TRUE.
        ELSE
          value(i)=.FALSE.
        END IF
      END DO
!
!  Deallocate.
!
      IF (ALLOCATED(S)) DEALLOCATE (S)
!
  10  FORMAT (/,' YAML_GET_LVAR_1D - Inconsistent number of values,',   &
     &        ' key-string: ',a,/,20x,'YAML size = ',i0,                &
              ', Variable size = ',i0,/,20x,'File: ',a)
!
      RETURN
      END FUNCTION yaml_Get_lvar_1d
!
      FUNCTION yaml_Get_rvar_0d (self, keystring, value) RESULT (status)
!
!***********************************************************************
!                                                                      !
!  It gets scalar floating-point data from YAML dictionary object.     !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     self       YAML tree dictionary object (CLASS yaml_tree)         !
!     keystring  YAML tree aggregated keys (string)                    !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     value      YAML value (real; scalar)                             !
!     status     Error code (integer)                                  !
!                                                                      !
!***********************************************************************
!
!  Imported variable declarations.
!
      CLASS (yaml_tree), intent(in) :: self
      character (len=*), intent(in) :: keystring
      real (kind_real), intent(out) :: value
!
!  Local variable declarations.
!
      TYPE (yaml_extract), allocatable :: S(:)
!
      integer :: Nvalues, ie, status
!
      character (len=Lmax) :: msg

      character (len=*), parameter :: MyFile =                          &
     &   __FILE__//", yaml_Get_rvar_0d"
!
!-----------------------------------------------------------------------
!  Extract requested key-string value.
!-----------------------------------------------------------------------
!
      status=NoErr
!
      IF (yaml_Error(self%extract(keystring, S),                        &
     &               NoErr, __LINE__, MyFile)) RETURN
!
!  Make sure that extracted value is a scalar.
!
      Nvalues=SIZE(S)
!
      IF (Nvalues.gt.1) THEN
        yaml_ErrFlag=7
        status=yaml_ErrFlag
        IF (yaml_Error(yaml_ErrFlag, NoErr, __LINE__, MyFile)) THEN
          IF (yaml_Master) WRITE (yaml_stdout,10) keystring, Nvalues,   &
     &                                            self%filename
          RETURN
        END IF
!
!  Convert string value to real.
!
      ELSE
        S(1)%value=ADJUSTL(S(1)%value)
        ie=LEN_TRIM(S(1)%value)
        READ (S(1)%value(1:ie), *, IOSTAT=status, IOMSG=msg) value
        IF (yaml_Error(status, NoErr, __LINE__, MyFile)) THEN
          yaml_ErrFlag=5
          IF (yaml_Master) WRITE (yaml_stdout,20) TRIM(keystring),      &
     &                                            S(1)%value, TRIM(msg)
          RETURN
        END IF
      END IF
!
!  Deallocate.
!
      IF (ALLOCATED(S)) DEALLOCATE (S)
!
  10  FORMAT (/,' YAML_GET_RVAR_0D - Extracted value is not a scalar,', &
     &        ' key-string: ',a,/,20x,'Nvalues = ',i0,/,20x,'File: ',a)
  20  FORMAT (/,' YAML_GET_RVAR_0D - Error while converting string to', &
     &        ' integer, key: ',a,', value = ',a,/,20x,'ErrMsg: ',a)
!
      RETURN
      END FUNCTION yaml_Get_rvar_0d
!
      FUNCTION yaml_Get_rvar_1d (self, keystring, value) RESULT (status)
!
!***********************************************************************
!                                                                      !
!  It gets 1D floating-point data from YAML dictionary object.         !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     self       YAML tree dictionary object (CLASS yaml_tree)         !
!     keystring  YAML tree aggregated keys (string)                    !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     value      YAML value (real; 1D-array)                           !
!     status     Error code (integer)                                  !
!                                                                      !
!***********************************************************************
!
!  Imported variable declarations.
!
      CLASS (yaml_tree), intent(in) :: self
      character (len=*), intent(in) :: keystring
      real (kind_real), intent(out) :: value(:)
!
!  Local variable declarations.
!
      TYPE (yaml_extract), allocatable :: S(:)
!
      integer :: Nvalues, i, ie, status
!
      character (len=Lmax) :: msg

      character (len=*), parameter :: MyFile =                          &
     &   __FILE__//", yaml_Get_rvar_1d"
!
!-----------------------------------------------------------------------
!  Extract requested key-string values.
!-----------------------------------------------------------------------
!
      status=NoErr
!
      IF (yaml_Error(self%extract(keystring, S),                        &
     &               NoErr, __LINE__, MyFile)) RETURN
!
!  Make sure that extracted value is a scalar.
!
      Nvalues=SIZE(S, DIM=1)
!
      IF (SIZE(value, DIM=1).lt.Nvalues) THEN
        yaml_ErrFlag=7
        status=yaml_ErrFlag
        IF (yaml_Error(yaml_ErrFlag, NoErr, __LINE__, MyFile)) THEN
          IF (yaml_Master) WRITE (yaml_stdout,10) keystring, Nvalues,   &
     &                                            SIZE(value, DIM=1),   &
     &                                            self%filename
          RETURN
        END IF
      END IF
!
!  Convert string values to reals.
!
      DO i=1,Nvalues
        S(i)%value=ADJUSTL(S(i)%value)
        ie=LEN_TRIM(S(i)%value)
        READ (S(i)%value(1:ie), *, IOSTAT=status, IOMSG=msg) value(i)
        IF (yaml_Error(status, NoErr, __LINE__, MyFile)) THEN
          yaml_ErrFlag=5
          IF (yaml_Master) WRITE (yaml_stdout,20) TRIM(keystring),      &
     &                                            S(i)%value, TRIM(msg)
          RETURN
        END IF
      END DO
!
!  Deallocate.
!
      IF (ALLOCATED(S)) DEALLOCATE (S)
!
  10  FORMAT (/,' YAML_GET_RVAR_1D - Inconsistent number of values,',   &
     &        ' key-string: ',a,/,20x,'YAML size = ',i0,                &
              ', Variable size = ',i0,/,20x,'File: ',a)
  20  FORMAT (/,' YAML_GET_RVAR_1D - Error while converting string to', &
     &        ' real, key: ',a,', value = ',a,/,20x,'ErrMsg: ',a)
!
      RETURN
      END FUNCTION yaml_Get_rvar_1d
!
      FUNCTION yaml_Get_svar_0d (self, keystring, value) RESULT (status)
!
!***********************************************************************
!                                                                      !
!  It gets scalar string data from YAML dictionary object.             !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     self       YAML tree dictionary object (CLASS yaml_tree)         !
!     keystring  YAML tree aggregated keys (string)                    !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     value      YAML value (string; scalar)                           !
!     status     Error code (integer)                                  !
!                                                                      !
!***********************************************************************
!
!  Imported variable declarations.
!
      CLASS (yaml_tree), intent(in ) :: self
      character (len=*), intent(in ) :: keystring
      character (len=*), intent(out) :: value
!
!  Local variable declarations.
!
      TYPE (yaml_extract), allocatable :: S(:)
!
      integer :: Nvalues, status
!
      character (len=*), parameter :: MyFile =                          &
     &   __FILE__//", yaml_Get_lvar_0d"
!
!-----------------------------------------------------------------------
!  Extract requested key-string value.
!-----------------------------------------------------------------------
!
      status=NoErr
!
      IF (yaml_Error(self%extract(keystring, S),                        &
     &               NoErr, __LINE__, MyFile)) RETURN
!
!  Make sure that extracted value is a scalar.
!
      Nvalues=SIZE(S)
!
      IF (Nvalues.gt.1) THEN
        yaml_ErrFlag=7
        status=yaml_ErrFlag
        IF (yaml_Error(yaml_ErrFlag, NoErr, __LINE__, MyFile)) THEN
          IF (yaml_Master) WRITE (yaml_stdout,10) keystring, Nvalues,   &
     &                                            self%filename
          RETURN
        END IF
!
!  Load string value.
!
      ELSE
        value=S(1)%value
      END IF
!
!  Deallocate.
!
      IF (ALLOCATED(S)) DEALLOCATE (S)
!
  10  FORMAT (/,' YAML_GET_SVAR_0D - Extracted value is not a scalar,', &
     &        ' key-string: ',a,/,20x,'Nvalues = ',i0,/,20x,'File: ',a)
!
      RETURN
      END FUNCTION yaml_Get_svar_0d
!
      FUNCTION yaml_Get_svar_1d (self, keystring, value) RESULT (status)
!
!***********************************************************************
!                                                                      !
!  It gets 1D string data from YAML dictionary object.                 !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     self       YAML tree dictionary object (CLASS yaml_tree)         !
!     keystring  YAML tree aggregated keys (string)                    !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     value      YAML value (string; 1D-array)                         !
!     status     Error code (integer)                                  !
!                                                                      !
!***********************************************************************
!
!  Imported variable declarations.
!
      CLASS (yaml_tree), intent(in ) :: self
      character (len=*), intent(in ) :: keystring
      character (len=*), intent(out) :: value(:)
!
!  Local variable declarations.
!
      TYPE (yaml_extract), allocatable :: S(:)
!
      integer :: Nvalues, i, status
!
      character (len=*), parameter :: MyFile =                          &
     &   __FILE__//", yaml_Get_lvar_1d"
!
!-----------------------------------------------------------------------
!  Extract requested key-string values.
!-----------------------------------------------------------------------
!
      status=NoErr
!
      IF (yaml_Error(self%extract(keystring, S),                        &
     &               NoErr, __LINE__, MyFile)) RETURN
!
!  Make sure that extracted value is a scalar.
!
      Nvalues=SIZE(S, DIM=1)
!
      IF (SIZE(value, DIM=1).lt.Nvalues) THEN
        yaml_ErrFlag=7
        status=yaml_ErrFlag
        IF (yaml_Error(yaml_ErrFlag, NoErr, __LINE__, MyFile)) THEN
          IF (yaml_Master) WRITE (yaml_stdout,10) keystring, Nvalues,   &
     &                                            SIZE(value, DIM=1),   &
     &                                            self%filename
          RETURN
        END IF
      END IF
!
!  Load string values.
!
      DO i=1,Nvalues
        value(i)=S(i)%value
      END DO
!
!  Deallocate.
!
      IF (ALLOCATED(S)) DEALLOCATE (S)
!
  10  FORMAT (/,' YAML_GET_SVAR_1D - Inconsistent number of values,',   &
     &        ' key-string: ',a,/,20x,'YAML size = ',i0,                &
              ', Variable size = ',i0,/,20x,'File: ',a)
!
      RETURN
      END FUNCTION yaml_Get_svar_1d
!
      FUNCTION yaml_LowerCase (Sinp) RESULT (Sout)
!
!=======================================================================
!                                                                      !
!  It converts all string elements to lowercase.                       !
!                                                                      !
!  Reference:                                                          !
!                                                                      !
!     Cooper Redwine, 1995: "Upgrading to Fortran 90", Springer-       !
!       Verlag, New York, pp 416.                                      !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      character(len=*), intent(in) :: Sinp
!
!  Local variable definitions.
!
      integer :: Lstr, i, j

      character (len=LEN(Sinp)) :: Sout

      character (len=26), parameter :: L = 'abcdefghijklmnopqrstuvwxyz'
      character (len=26), parameter :: U = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
!
!-----------------------------------------------------------------------
!  Convert input string to lowercase.
!-----------------------------------------------------------------------
!
      Lstr=LEN(Sinp)
      Sout=Sinp
      DO i=1,Lstr
        j=INDEX(U, Sout(i:i))
        IF (j.ne.0) THEN
          Sout(i:i)=L(j:j)
        END IF
      END DO
!
      RETURN
      END FUNCTION yaml_LowerCase
!
      FUNCTION yaml_UpperCase (Sinp) RESULT (Sout)
!
!=======================================================================
!                                                                      !
!  It converts all string elements to uppercase.                       !
!                                                                      !
!  Reference:                                                          !
!                                                                      !
!     Cooper Redwine, 1995: "Upgrading to Fortran 90", Springer-       !
!       Verlag, New York, pp 416.                                      !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      character (len=*), intent(in) :: Sinp
!
!  Local variable definitions.
!
      integer :: Lstr, i, j
!
      character (len=LEN(Sinp)) :: Sout
!
      character (len=26), parameter :: L = 'abcdefghijklmnopqrstuvwxyz'
      character (len=26), parameter :: U = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
!
!-----------------------------------------------------------------------
!  Convert input string to uppercase.
!-----------------------------------------------------------------------
!
      Lstr=LEN(Sinp)
      Sout=Sinp
      DO i=1,Lstr
        j=INDEX(L, Sout(i:i))
        IF (j.ne.0) THEN
          Sout(i:i)=U(j:j)
        END IF
      END DO
!
      RETURN
      END FUNCTION yaml_UpperCase
!
      SUBROUTINE yaml_ValueType (value, Lswitch)
!
!***********************************************************************
!                                                                      !
!  It determines the YAML value type as logical, integer, real, or     !
!  string.                                                             !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     self       YAML pair value (string)                              !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     Lswitch    YAML key/value switches (logical, vector)             !
!                                                                      !
!                  Lswitch(1) = T/F, value has an alias '*' token      !
!                  Lswitch(2) = T/F, value has an anchor '&' token     !
!                  Lswitch(3) = T/F, key/value starts a block '-'      !
!                  Lswitch(4) = T/F, value has sequence '[...]' tokens !
!                  Lswitch(5) = T/F, logical value(s)                  !
!                  Lswitch(6) = T/F, integer value(s)                  !
!                  Lswitch(7) = T/F, floating-point value(s)           !
!                  Lswitch(8) = T/F, string value(s)                   !
!                                                                      !
!***********************************************************************
!
!  Imported variable declarations.
!
      logical, intent(inout) :: Lswitch(:)
!
      character (len=*), intent(inout)  :: value
!
!  Local variable declarations.
!
      integer :: Lstr, Schar
      integer :: colon, dot, exponent, letter, numeric, precision
      integer :: i, multiple
!
!-----------------------------------------------------------------------
!  Set keyword value type.
!-----------------------------------------------------------------------
!
!  Initialize.
!
      colon=0
      dot=0
      exponent=0
      letter=0
      multiple=0
      numeric=0
      precision=0
!
!  Check input value string.
!
      Lstr=LEN_TRIM(value)
      IF (Lstr.gt.0) THEN
        DO i=1,lstr
          Schar=ICHAR(value(i:i))
!
!  Check for numbers, plus, and minus signs. For example, value=-1.0E+37
!  is a floating-point numerical value.
!
          IF (((48.le.Schar).and.(Schar.le.57)).or.                     &
     &        (Schar.eq.43).or.(Schar.eq.45)) numeric=numeric+1
!
!  Check for dot/period, CHAR(46). If period is not found in a numerical
!  expression, identify value as an integer.
!
          IF (Schar.eq.46) dot=dot+1
!
!  Check for precision character: D=ICHAR(68) and d=ICHAR(100). For
!  example, value=0.0d0 and others in the future.
!
          IF ((Schar.eq.69).or.(Schar.eq.101)) precision=precision+1
!
!  Check for exponent character: E=CHAR(69) and e=CHAR(101).
!
          IF ((Schar.eq.69).or.(Schar.eq.101)) exponent=exponent+1
!
!  Check for multiple values separate by comma, CHAR(44), in a sequence
!  of values: [val1, val2, ..., valN].
!
          IF (Lswitch(4)) multiple=multiple+1
!
!  Check for colon, CHAR(58). We can have value=2020-01-01T00:00:00Z,
!  which has numbers, hyphen, and letters. It is the colon that makes
!  it a string variable (https://www.myroms.org).
!
          IF (Schar.eq.58) colon=colon+1
!
!  English uppercase and lowercase alphabet.
!
          IF (((65.le.Schar).and.(Schar.le.90)).or.                     &
     &        (Schar.eq.97).or.(Schar.eq.122)) letter=letter+1
        END DO
!
!  Set integer or floating-point type.
!
        IF ((numeric.gt.0).and.(colon.eq.0)) THEN
          IF ((dot.gt.0).or.(exponent.gt.0).or.(precision.gt.0)) THEN
            Lswitch(7)=.TRUE.              ! floating-point
          ELSE
            Lswitch(6)=.TRUE.              ! integer
          END IF
        ELSE
          IF ((value.eq.'true').or.(value.eq.'false').or.              &
     &        (value.eq.'TRUE').or.(value.eq.'FALSE')) THEN
            Lswitch(5)=.TRUE.              ! logical
          ELSE IF (letter.gt.0) THEN
            Lswitch(8)=.TRUE.              ! string
          END IF
        END IF
      END IF
!
      RETURN
      END SUBROUTINE yaml_ValueType

!-----------------------------------------------------------------------

      END MODULE yaml_parser_mod
