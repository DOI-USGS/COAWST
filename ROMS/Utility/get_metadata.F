#include "cppdefs.h"
      MODULE get_metadata_mod
!
!git $Id$
!svn $Id: get_metadata.F 1151 2023-02-09 03:08:53Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2023 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This module has functions to process ROMS I/O metadata or coupling  !
!  metadata files.  Two formats are supported: native and YAML files.  !
!  The native format is deprecated since YAML files are expandable and !
!  portable.                                                           !
!                                                                      !
!  io_metadata               It processes entries in ROMS input/output !
!                              variables metadata:                     !
!                                                                      !
!                              'varinfo.dat'  or  'varinfo.yaml'       !
!                                                                      !
!  coupling_metadata         It processes dictionary entries for the   !
!                              ESMF/NUOPC coupling system:             !
!                                                                      !
!                              'coupling_*.dat'  or  'coupling_*.yaml' !
!                                                                      !
!=======================================================================
!
      USE mod_kinds,       ONLY : dp
      USE mod_parallel,    ONLY : Master
      USE mod_iounits,     ONLY : stdout, varname
      USE mod_scalars,     ONLY : exit_flag, NoError
      USE strings_mod,     ONLY : FoundError, assign_string, lowercase
      USE yaml_parser_mod, ONLY : yaml_initialize,                      &
     &                            yaml_get,                             &
     &                            yaml_Svec,                            &
     &                            yaml_tree
!
      implicit none
!
!-----------------------------------------------------------------------
!  Define generic coupling field to process import and export states.
!-----------------------------------------------------------------------
!
      TYPE CouplingField
        logical :: connected
        logical :: debug_write
!
        real(dp) :: add_offset
        real(dp) :: scale
!
        character (len=:), allocatable :: connected_to
        character (len=:), allocatable :: data_netcdf_vname
        character (len=:), allocatable :: data_netcdf_tname
        character (len=:), allocatable :: destination_grid
        character (len=:), allocatable :: destination_units
        character (len=:), allocatable :: extrapolate_method
        character (len=:), allocatable :: long_name
        character (len=:), allocatable :: map_norm
        character (len=:), allocatable :: map_type
        character (len=:), allocatable :: regrid_method
        character (len=:), allocatable :: source_units
        character (len=:), allocatable :: source_grid
        character (len=:), allocatable :: short_name
        character (len=:), allocatable :: standard_name
!
      END TYPE CouplingField
!
!-----------------------------------------------------------------------
!  Define generic YAML dictionary, containers, and counters used
!  during processing.
!-----------------------------------------------------------------------
!
!  YAML dictionary object. It is destroyed after processing, so it
!  can be reused to operate on other input YAML files.
!
      TYPE (yaml_tree) :: YML
!
!  Metadata debugging and reporting switches
!
#ifdef METADATA_REPORT
      logical :: LdebugMetadata = .TRUE.
      logical :: LreportYAML    = .TRUE.
#else
      logical :: LdebugMetadata = .FALSE.
      logical :: LreportYAML    = .FALSE.
#endif
!
!  Counters.
!
      integer :: Ientry                            ! entry counter
      integer :: Nentries                          ! number of entries
!
!  logical scalar dummy values.
!
      logical, allocatable :: Ylogical1(:)
!
!  Real scalar dummy values.
!
      real(dp), allocatable :: Yreal1(:)
      real(dp), allocatable :: Yreal2(:)
!
!  Derived-type dummy structures for processing string value or set
!  of values from a sequence flow, [val1, ..., valN].
!
      TYPE (yaml_Svec), allocatable :: Ystring1 (:)
      TYPE (yaml_Svec), allocatable :: Ystring2 (:)
      TYPE (yaml_Svec), allocatable :: Ystring3 (:)
      TYPE (yaml_Svec), allocatable :: Ystring4 (:)
      TYPE (yaml_Svec), allocatable :: Ystring5 (:)
      TYPE (yaml_Svec), allocatable :: Ystring6 (:)
      TYPE (yaml_Svec), allocatable :: Ystring7 (:)
      TYPE (yaml_Svec), allocatable :: Ystring8 (:)
      TYPE (yaml_Svec), allocatable :: Ystring9 (:)
      TYPE (yaml_Svec), allocatable :: Ystring10(:)
      TYPE (yaml_Svec), allocatable :: Ystring11(:)
      TYPE (yaml_Svec), allocatable :: Ystring12(:)
!
      PUBLIC :: cmeps_metadata
      PUBLIC :: coupling_metadata
      PUBLIC :: io_metadata
      PUBLIC :: metadata_has
!
!-----------------------------------------------------------------------
      CONTAINS
!-----------------------------------------------------------------------
!
      SUBROUTINE cmeps_metadata (self, filename, key, S)
!
!=======================================================================
!                                                                      !
!  It process either import or export fields which are stored as block !
!  lists (leading key/value is hyphenated) in the YAML file. The YAML  !
!  file is used to configure ROMS ESMF/NUOPC 'cap' module to be run by !
!  the Community Mediator for Earth Prediction Systems (CMEPS).        !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     self      YAML tree dictionary (TYPE yaml_tree)                  !
!                                                                      !
!     filename  ROMS YAML configuration filename for CMEPS (string)    !
!                                                                      !
!     key       Leading blocking key to process (string), for example: !
!                 'export', 'import', or 'bulk_flux import'            !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     S         Import or Export coupling fields (TYPE CouplingField)  !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      TYPE (yaml_tree),  intent(inout) :: self
      character (len=*), intent(in   ) :: filename
      character (len=*), intent(in   ) :: key
      TYPE (CouplingField), allocatable, intent(out) :: S(:)
!
!  Local variable declarations.
!
      integer :: i
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", cmeps_metadata"
!
!-----------------------------------------------------------------------
!  Process coupling import or export metadata for CMEPS.
!-----------------------------------------------------------------------
!
!  If applicable, create YAML tree dictionary.
!
      IF (.not.ASSOCIATED(self%list)) THEN
        IF (FoundError(yaml_initialize(self, TRIM(filename),            &
     &                                 LreportYAML),                    &
     &                 NoError, __LINE__, MyFile)) THEN
          IF (Master) WRITE (stdout,10) TRIM(filename)
          RETURN
        END IF
      END IF
!
!  Extract requested blocking list.
!
      IF (FoundError(yaml_get(self, TRIM(key)//'.standard_name',        &
     &                        Ystring1),                                &
     &                 NoError, __LINE__, MyFile)) RETURN
      Nentries=SIZE(Ystring1)
!
      IF (FoundError(yaml_get(self, TRIM(key)//'.long_name',            &
     &                        Ystring2),                                &
     &               NoError, __LINE__, MyFile)) RETURN
!
      IF (FoundError(yaml_get(self, TRIM(key)//'.short_name',           &
     &                        Ystring3),                                &
     &                 NoError, __LINE__, MyFile)) RETURN
!
      IF (FoundError(yaml_get(self, TRIM(key)//'.data_variables',       &
     &                        Ystring4),                                &
     &               NoError, __LINE__, MyFile)) RETURN
!
      IF (FoundError(yaml_get(self, TRIM(key)//'.source_units',         &
     &                        Ystring5),                                &
     &               NoError, __LINE__, MyFile)) RETURN
!
      IF (FoundError(yaml_get(self, TRIM(key)//'.destination_units',    &
     &                        Ystring6),                                &
     &               NoError, __LINE__, MyFile)) RETURN
!
      IF (FoundError(yaml_get(self, TRIM(key)//'.source_grid',          &
     &                        Ystring7),                                &
     &               NoError, __LINE__, MyFile)) RETURN
!
      IF (FoundError(yaml_get(self, TRIM(key)//'.destination_grid',     &
     &                        Ystring8),                                &
     &               NoError, __LINE__, MyFile)) RETURN
!
      IF (FoundError(yaml_get(self, TRIM(key)//'.connected_to',         &
     &                        Ystring9),                                &
     &               NoError, __LINE__, MyFile)) RETURN
!
      IF (FoundError(yaml_get(self, TRIM(key)//'.map_type',             &
     &                        Ystring10),                               &
     &               NoError, __LINE__, MyFile)) RETURN
!
      IF (FoundError(yaml_get(self, TRIM(key)//'.map_norm',             &
     &                        Ystring11),                               &
     &               NoError, __LINE__, MyFile)) RETURN
!
      IF (.not.allocated(Yreal1)) THEN
        allocate ( Yreal1(Nentries) )
      END IF
      IF (FoundError(yaml_get(self, TRIM(key)//'.add_offset',           &
     &                        Yreal1),                                  &
     &               NoError, __LINE__, MyFile)) RETURN
!
      IF (.not.allocated(Yreal2)) THEN
        allocate ( Yreal2(Nentries) )
      END IF
      IF (FoundError(yaml_get(self, TRIM(key)//'.scale',                &
                              Yreal2),                                  &
     &               NoError, __LINE__, MyFile)) RETURN
      IF (.not.allocated(Yreal1)) THEN
        allocate ( Yreal1(Nentries) )
      END IF
!
      IF (.not.allocated(Ylogical1)) THEN
        allocate ( Ylogical1(Nentries) )
      END IF
      IF (FoundError(yaml_get(self, TRIM(key)//'.debug_write',          &
     &                        Ylogical1),                               &
     &               NoError, __LINE__, MyFile)) RETURN
!
!  Load metadata into output structure.
!
      IF (.not.allocated(S)) allocate ( S(Nentries) )
!
      DO i=1,Nentries
        S(i)%debug_write = Ylogical1(i)
        S(i)%add_offset  = Yreal1(i)
        S(i)%scale       = Yreal2(i)
!
        IF (FoundError(assign_string(S(i)%standard_name,                &
     &                               Ystring1(i)%value),                &
     &                 NoError, __LINE__, MyFile)) RETURN
!
        IF (FoundError(assign_string(S(i)%long_name,                    &
     &                               Ystring2(i)%value),                &
     &                 NoError, __LINE__, MyFile)) RETURN
!
        IF (FoundError(assign_string(S(i)%short_name,                   &
     &                               Ystring3(i)%value),                &
     &                 NoError, __LINE__, MyFile)) RETURN
!
        IF (FoundError(assign_string(S(i)%data_netcdf_vname,            &
     &                               Ystring4(i)%vector(1)%value),      &
     &                 NoError, __LINE__, MyFile)) RETURN
!
        IF (FoundError(assign_string(S(i)%data_netcdf_tname,            &
     &                               Ystring4(i)%vector(2)%value),      &
     &                 NoError, __LINE__, MyFile)) RETURN
!
        IF (FoundError(assign_string(S(i)%source_units,                 &
     &                               Ystring5(i)%value),                &
     &                 NoError, __LINE__, MyFile)) RETURN
!
        IF (FoundError(assign_string(S(i)%destination_units,            &
     &                               Ystring6(i)%value),                &
     &                 NoError, __LINE__, MyFile)) RETURN
!
        IF (FoundError(assign_string(S(i)%source_grid,                  &
     &                               Ystring7(i)%value),                &
     &                 NoError, __LINE__, MyFile)) RETURN
!
        IF (FoundError(assign_string(S(i)%destination_grid,             &
     &                               Ystring8(i)%value),                &
     &                 NoError, __LINE__, MyFile)) RETURN
!
        IF (FoundError(assign_string(S(i)%connected_to,                 &
     &                               Ystring9(i)%value),                &
     &                 NoError, __LINE__, MyFile)) RETURN
        IF (lowercase(S(i)%connected_to).eq.'false') THEN
          S(i)%connected=.FALSE.
        ELSE
          S(i)%connected=.TRUE.
        END IF
!
        IF (FoundError(assign_string(S(i)%map_type,                     &
     &                               Ystring10(i)%value),               &
     &                 NoError, __LINE__, MyFile)) RETURN
!
        IF (FoundError(assign_string(S(i)%map_norm,                     &
     &                               Ystring11(i)%value),               &
     &                 NoError, __LINE__, MyFile)) RETURN
      END DO
!
!  Deallocate generic structures.
!
      IF (allocated(Ystring1 ))  deallocate (Ystring1 )
      IF (allocated(Ystring2 ))  deallocate (Ystring2 )
      IF (allocated(Ystring3 ))  deallocate (Ystring3 )
      IF (allocated(Ystring4 ))  deallocate (Ystring4 )
      IF (allocated(Ystring5 ))  deallocate (Ystring5 )
      IF (allocated(Ystring6 ))  deallocate (Ystring6 )
      IF (allocated(Ystring7 ))  deallocate (Ystring7 )
      IF (allocated(Ystring8 ))  deallocate (Ystring8 )
      IF (allocated(Ystring9 ))  deallocate (Ystring9 )
      IF (allocated(Ystring10))  deallocate (Ystring10)
      IF (allocated(Ystring11))  deallocate (Ystring11)
      IF (allocated(Ylogical1))  deallocate (Ylogical1)
      IF (allocated(Yreal1))     deallocate (Yreal1)
      IF (allocated(Yreal2))     deallocate (Yreal2)
!
!  Report.
!
      IF (Master.and.LdebugMetadata) THEN
        WRITE (stdout,'(/,3a,/,3a)')                                    &
     &        "Coupling Metadata Dictionary, key: '", TRIM(key), "',",  &
     &        REPEAT('=',28), '  File: ', TRIM(filename)
        DO i=1,SIZE(S)
          WRITE (stdout,'(/,a,a)')      '  - standard_name:        ',   &
     &                                  TRIM(S(i)%standard_name)
          WRITE (stdout,'(a,a)')        '    long_name:            ',   &
     &                                  TRIM(S(i)%long_name)
          WRITE (stdout,'(a,a)')        '    short_name:           ',   &
     &                                  TRIM(S(i)%short_name)
          WRITE (stdout,'(a,a)')        '    data_netcdf_variable: ',   &
     &                                  TRIM(S(i)%data_netcdf_vname)
          WRITE (stdout,'(a,a)')        '    data_netcdf_time:     ',   &
     &                                  TRIM(S(i)%data_netcdf_tname)
          WRITE (stdout,'(a,a)')        '    source_units:         ',   &
     &                                  TRIM(S(i)%source_units)
          WRITE (stdout,'(a,a)')        '    destination_units:    ',   &
     &                                  TRIM(S(i)%destination_units)
          WRITE (stdout,'(a,a)')        '    source_grid:          ',   &
     &                                  TRIM(S(i)%source_grid)
          WRITE (stdout,'(a,a)')        '    destination_grid:     ',   &
     &                                  TRIM(S(i)%destination_grid)
          WRITE (stdout,'(a,1p,e15.8)') '    add_offset:           ',   &
     &                                  S(i)%add_offset
          WRITE (stdout,'(a,1p,e15.8)') '    scale:                ',   &
     &                                  S(i)%scale
          WRITE (stdout,'(a,l1)')       '    debug_write:          ',   &
     &                                  S(i)%debug_write
          WRITE (stdout,'(a,l1)')       '    connected:            ',   &
     &                                  S(i)%connected
          WRITE (stdout,'(a,a)')        '    connected_to:         ',   &
     &                                  TRIM(S(i)%connected_to)
          WRITE (stdout,'(a,a)')        '    map_type:             ',   &
     &                                  TRIM(S(i)%map_type)
          WRITE (stdout,'(a,a)')        '    map_norm:             ',   &
     &                                  TRIM(S(i)%map_norm)
        END DO
        FLUSH (stdout)
      END IF
!
  10  FORMAT (/,' CMEPS_METADATA - Unable to create YAML object for',    &
     &        ' ROMS/CMEPS configuration metadata file: ',/,21x,a,/,     &
     &        21x,'Default file is located in source directory.')
!
      RETURN
      END SUBROUTINE cmeps_metadata
!
      SUBROUTINE coupling_metadata (filename, S)
!
!=======================================================================
!                                                                      !
!  It processes import and export field dictionary for ROMS coupling   !
!  system with the ESMF/NUOPC library. If processes field metadata     !
!  entry-by-entry from 'coupling_*.dat'  or  'coupling_*.yaml'.        !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     filename  Coupling metadata filename (string)                    !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     S         Import/Export coupling fields (TYPE CouplingField)     !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      character (len=*), intent(in) :: filename
!
      TYPE (CouplingField), allocatable, intent(out) :: S(:)
!
!  Local variable declarations.
!
      logical :: IsDat, IsYaml, connected, debug_write
!
      real(dp) :: add_offset, scale
!
      integer, parameter :: iunit = 10
      integer :: Idot, Lstr, Lvar, i, io_err
!
      character (len=40 ) :: Smodel, Tname
      character (len=100) :: Cinfo(12)
      character (len=256) :: io_errmsg

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", coupling_metadata"
!
!-----------------------------------------------------------------------
!  Process coupling import/export metadata.
!-----------------------------------------------------------------------
!
!  Determine metadata file extension: 'coupling_*.dat'  or
!                                     'coupling_*.yaml'
!
      IsDat =.FALSE.
      IsYAML=.FALSE.

      Lstr=LEN_TRIM(filename)
      Idot=INDEX(filename(1:Lstr), CHAR(46), BACK=.TRUE.)
!
      SELECT CASE (lowercase(filename(Idot+1:Lstr)))
        CASE ('dat')
          IsDat=.TRUE.
        CASE ('yaml', 'yml')
          IsYaml=.TRUE.
      END SELECT
!
!  If YAML metadata, create dictionary.
!
      IF (IsYaml) THEN
        Ientry=0
!
        IF (FoundError(yaml_initialize(YML, TRIM(filename),             &
     &                                 LreportYAML),                    &
     &                 NoError, __LINE__, MyFile)) THEN
          IF (Master) WRITE (stdout,30) TRIM(filename)
          RETURN
        END IF
!
!  If YAML metadata, extract key/value pair (blocking list).
!
        IF (FoundError(yaml_get(YML, 'metadata.standard_name',          &
     &                          Ystring1),                              &
     &                 NoError, __LINE__, MyFile)) RETURN
        Nentries=SIZE(Ystring1)
!
        IF (FoundError(yaml_get(YML, 'metadata.long_name',              &
     &                          Ystring2),                              &
     &                 NoError, __LINE__, MyFile)) RETURN
!
        IF (FoundError(yaml_get(YML, 'metadata.short_name',             &
     &                          Ystring3),                              &
     &                   NoError, __LINE__, MyFile)) RETURN
!
        IF (FoundError(yaml_get(YML, 'metadata.data_variables',         &
     &                          Ystring4),                              &
     &                 NoError, __LINE__, MyFile)) RETURN
!
        IF (FoundError(yaml_get(YML, 'metadata.source_units',           &
     &                          Ystring5),                              &
     &                 NoError, __LINE__, MyFile)) RETURN
!
        IF (FoundError(yaml_get(YML, 'metadata.destination_units',      &
     &                          Ystring6),                              &
     &                 NoError, __LINE__, MyFile)) RETURN
!
        IF (FoundError(yaml_get(YML, 'metadata.source_grid',            &
     &                          Ystring7),                              &
     &                 NoError, __LINE__, MyFile)) RETURN
!
        IF (FoundError(yaml_get(YML, 'metadata.destination_grid',       &
     &                          Ystring8),                              &
     &                 NoError, __LINE__, MyFile)) RETURN
!
        IF (FoundError(yaml_get(YML, 'metadata.connected_to',           &
     &                          Ystring9),                              &
     &                 NoError, __LINE__, MyFile)) RETURN
!
        IF (FoundError(yaml_get(YML, 'metadata.regrid_method',          &
     &                          Ystring10),                             &
     &                 NoError, __LINE__, MyFile)) RETURN
!
        IF (FoundError(yaml_get(YML, 'metadata.extrapolate_method',     &
     &                          Ystring11),                             &
     &                 NoError, __LINE__, MyFile)) RETURN
!
        IF (allocated(Yreal1)) deallocate (Yreal1)
        allocate ( Yreal1(Nentries) )
        IF (FoundError(yaml_get(YML, 'metadata.add_offset',             &
     &                          Yreal1),                                &
     &                 NoError, __LINE__, MyFile)) RETURN
!
        IF (allocated(Yreal2)) deallocate (Yreal2)
        allocate ( Yreal2(Nentries) )
        IF (FoundError(yaml_get(YML, 'metadata.scale',                  &
                                Yreal2),                                &
     &                 NoError, __LINE__, MyFile)) RETURN
!
        IF (allocated(Ylogical1)) deallocate (Ylogical1)
        allocate ( Ylogical1(Nentries) )
        IF (FoundError(yaml_get(YML, 'metadata.debug_write',            &
     &                          Ylogical1),                             &
     &                 NoError, __LINE__, MyFile)) RETURN
!
!  Otherwise, open deprecated 'coupling_*.dat' file.
!
      ELSE IF (IsDat) THEN
        OPEN (UNIT=iunit, FILE=TRIM(filename), FORM='formatted',        &
     &        STATUS='old', IOSTAT=io_err, IOMSG=io_errmsg)
        IF (FoundError(io_err, NoError, __LINE__, MyFile)) THEN
          exit_flag=5
          IF (Master) WRITE(stdout,40) TRIM(filename), TRIM(io_errmsg)
          RETURN
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Load metadata information from YAML structures.
!-----------------------------------------------------------------------
!
      IF (IsYaml) THEN
        IF (.not.allocated(S)) allocate ( S(Nentries) )
!
        DO i=1,Nentries
          S(i)%debug_write = Ylogical1(i)
          S(i)%add_offset  = Yreal1(i)
          S(i)%scale       = Yreal2(i)
!
          IF (FoundError(assign_string(S(i)%standard_name,              &
     &                                 Ystring1(i)%value),              &
     &                   NoError, __LINE__, MyFile)) RETURN
!
          IF (FoundError(assign_string(S(i)%long_name,                  &
     &                                 Ystring2(i)%value),              &
     &                   NoError, __LINE__, MyFile)) RETURN
!
          IF (FoundError(assign_string(S(i)%short_name,                 &
     &                                 Ystring3(i)%value),              &
     &                   NoError, __LINE__, MyFile)) RETURN
!
          IF (FoundError(assign_string(S(i)%data_netcdf_vname,          &
     &                                 Ystring4(i)%vector(1)%value),    &
     &                   NoError, __LINE__, MyFile)) RETURN
!
          IF (FoundError(assign_string(S(i)%data_netcdf_tname,          &
     &                                 Ystring4(i)%vector(2)%value),    &
     &                   NoError, __LINE__, MyFile)) RETURN
!
          IF (FoundError(assign_string(S(i)%source_units,               &
     &                                 Ystring5(i)%value),              &
     &                   NoError, __LINE__, MyFile)) RETURN
!
          IF (FoundError(assign_string(S(i)%destination_units,          &
     &                                 Ystring6(i)%value),              &
     &                   NoError, __LINE__, MyFile)) RETURN
!
          IF (FoundError(assign_string(S(i)%source_grid,                &
     &                                 Ystring7(i)%value),              &
     &                   NoError, __LINE__, MyFile)) RETURN
!
          IF (FoundError(assign_string(S(i)%destination_grid,           &
     &                                 Ystring8(i)%value),              &
     &                   NoError, __LINE__, MyFile)) RETURN
!
          IF (FoundError(assign_string(S(i)%connected_to,               &
     &                                 Ystring9(i)%value),             &
     &                   NoError, __LINE__, MyFile)) RETURN
          IF (lowercase(S(i)%connected_to).eq.'false') THEN
            S(i)%connected=.FALSE.
          ELSE
            S(i)%connected=.TRUE.
          END IF
!
          IF (FoundError(assign_string(S(i)%regrid_method,              &
     &                                 Ystring10(i)%value),             &
     &                   NoError, __LINE__, MyFile)) RETURN
!
          IF (FoundError(assign_string(S(i)%extrapolate_method,         &
     &                                 Ystring11(i)%value),             &
     &                   NoError, __LINE__, MyFile)) RETURN
        END DO
!
!  Deallocate generic structures.
!
        CALL YML%destroy ()
        IF (allocated(Ystring1 ))  deallocate (Ystring1 )
        IF (allocated(Ystring2 ))  deallocate (Ystring2 )
        IF (allocated(Ystring3 ))  deallocate (Ystring3 )
        IF (allocated(Ystring4 ))  deallocate (Ystring4 )
        IF (allocated(Ystring5 ))  deallocate (Ystring5 )
        IF (allocated(Ystring6 ))  deallocate (Ystring6 )
        IF (allocated(Ystring7 ))  deallocate (Ystring7 )
        IF (allocated(Ystring8 ))  deallocate (Ystring8 )
        IF (allocated(Ystring9 ))  deallocate (Ystring9 )
        IF (allocated(Ystring10))  deallocate (Ystring10)
        IF (allocated(Ystring11))  deallocate (Ystring11)
        IF (allocated(Ylogical1))  deallocate (Ylogical1)
        IF (allocated(Yreal1))     deallocate (Yreal1)
        IF (allocated(Yreal2))     deallocate (Yreal2)
!
!-----------------------------------------------------------------------
!  Read in '*.dat' file and load metadata entries into output structure.
!-----------------------------------------------------------------------
!
      ELSE
!
!  Inquire number of valid entries in metadata file.
!
        Ientry=0
        DO WHILE (.TRUE.)
          READ (iunit,*,ERR=20,END=10) Cinfo( 1)
          Lvar=LEN_TRIM(Cinfo(1))
          IF ((Lvar.gt.0).and.(Cinfo(1)(1:1).ne.CHAR(33))) THEN
            Ientry=Ientry+1
            READ (iunit,*,ERR=20,END=10) Cinfo( 2)
            READ (iunit,*,ERR=20,END=10) Cinfo( 3)
            READ (iunit,*,ERR=20,END=10) Cinfo( 4)
            READ (iunit,*,ERR=20,END=10) Cinfo( 5)
            READ (iunit,*,ERR=20,END=10) Cinfo( 6)
            READ (iunit,*,ERR=20,END=10) Cinfo( 7)
            READ (iunit,*,ERR=20,END=10) Cinfo( 8)
            READ (iunit,*,ERR=20,END=10) Cinfo( 9)
            READ (iunit,*,ERR=20,END=10) Cinfo(10)
            READ (iunit,*,ERR=20,END=10) Cinfo(11)
            READ (iunit,*,ERR=20,END=10) Cinfo(12)
            READ (iunit,*,ERR=20,END=10) connected
            READ (iunit,*,ERR=20,END=10) debug_write
            READ (iunit,*,ERR=20,END=10) add_offset
            READ (iunit,*,ERR=20,END=10) scale
          END IF
        END DO
   10   CONTINUE
!
!  Allocate ouput structure.
!
        Nentries=Ientry
        IF (.not.allocated(S)) allocate ( S(Nentries) )
!
!  Rewind input unit, reread metadata information.
!
        REWIND (iunit)
!
        Ientry=0
        DO WHILE (Ientry.lt.Nentries)
          READ (iunit,*,ERR=20) Cinfo( 1)           ! short_name
          Lvar=LEN_TRIM(Cinfo(1))
          IF ((Lvar.gt.0).and.                                          &
              (Cinfo(1)(1:1).ne.CHAR(33))) THEN
            READ (iunit,*,ERR=20) Cinfo( 2)         ! standard_name
            READ (iunit,*,ERR=20) Cinfo( 3)         ! long_name
            READ (iunit,*,ERR=20) Cinfo( 4), Smodel ! connected_to
            READ (iunit,*,ERR=20) Cinfo( 5)         ! source_units
            READ (iunit,*,ERR=20) Cinfo( 6)         ! source_grid
            READ (iunit,*,ERR=20) Cinfo( 7)         ! data_short_name
            READ (iunit,*,ERR=20) Cinfo( 8)         ! destination_units
            READ (iunit,*,ERR=20) Cinfo( 9)         ! destination_grid
            READ (iunit,*,ERR=20) Cinfo(10), Tname  ! data_variables
            READ (iunit,*,ERR=20) Cinfo(11)         ! regrid_method
            READ (iunit,*,ERR=20) Cinfo(12)         ! extrapolate_method
            READ (iunit,*,ERR=20) connected
            READ (iunit,*,ERR=20) debug_write
            READ (iunit,*,ERR=20) add_offset
            READ (iunit,*,ERR=20) scale
            Ientry=Ientry+1
!
!  Load metadata into output structure.
!
            S(Ientry)%connected   = connected
            S(Ientry)%debug_write = debug_write
            S(Ientry)%add_offset  = add_offset
            S(Ientry)%scale       = scale
!
            IF (FoundError(assign_string(S(Ientry)%short_name,          &
     &                                   TRIM(ADJUSTL(Cinfo(1)))),      &
     &                     NoError, __LINE__, MyFile)) RETURN
!
            IF (FoundError(assign_string(S(Ientry)%standard_name,       &
     &                                   TRIM(ADJUSTL(Cinfo(2)))),      &
     &                     NoError, __LINE__, MyFile)) RETURN
!
            IF (FoundError(assign_string(S(Ientry)%long_name,           &
     &                                   TRIM(ADJUSTL(Cinfo(3)))),      &
     &                     NoError, __LINE__, MyFile)) RETURN
!
            IF (FoundError(assign_string(S(Ientry)%connected_to,        &
     &                                   TRIM(ADJUSTL(Smodel))),        &
     &                     NoError, __LINE__, MyFile)) RETURN
!
            IF (FoundError(assign_string(S(Ientry)%source_units,        &
     &                                   TRIM(ADJUSTL(Cinfo(5)))),      &
     &                     NoError, __LINE__, MyFile)) RETURN
!
            IF (FoundError(assign_string(S(Ientry)%source_grid,         &
     &                                   TRIM(ADJUSTL(Cinfo(6)))),      &
     &                     NoError, __LINE__, MyFile)) RETURN
!
            IF (FoundError(assign_string(S(Ientry)%destination_units,   &
     &                                   TRIM(ADJUSTL(Cinfo(8)))),      &
     &                     NoError, __LINE__, MyFile)) RETURN
!
            IF (FoundError(assign_string(S(Ientry)%destination_grid,    &
     &                                   TRIM(ADJUSTL(Cinfo(9)))),      &
     &                     NoError, __LINE__, MyFile)) RETURN
!
            IF (FoundError(assign_string(S(Ientry)%data_netcdf_vname,   &
     &                                   TRIM(ADJUSTL(Cinfo(10)))),     &
     &                     NoError, __LINE__, MyFile)) RETURN
!
            IF (FoundError(assign_string(S(Ientry)%data_netcdf_tname,   &
     &                                   TRIM(ADJUSTL(Tname))),         &
     &                     NoError, __LINE__, MyFile)) RETURN
!
            IF (FoundError(assign_string(S(Ientry)%regrid_method,       &
     &                                   TRIM(ADJUSTL(Cinfo(11)))),     &
     &                     NoError, __LINE__, MyFile)) RETURN
!
            IF (FoundError(assign_string(S(Ientry)%extrapolate_method,  &
     &                                   TRIM(ADJUSTL(Cinfo(12)))),     &
     &                     NoError, __LINE__, MyFile)) RETURN
          END IF
        END DO
        CLOSE (iunit)
      END IF
!
!  Report.
!
      IF (Master.and.LdebugMetadata) THEN
        WRITE (stdout,'(/,2a,/,a)')                                     &
     &        'Coupling Metadata Dictionary, File: ',                   &
     &        TRIM(filename), REPEAT('=',28)
        DO i=1,SIZE(S)
          WRITE (stdout,'(/,a,a)')      '  - standard_name:        ',   &
     &                                  TRIM(S(i)%standard_name)
          WRITE (stdout,'(a,a)')        '    long_name:            ',   &
     &                                  TRIM(S(i)%long_name)
          WRITE (stdout,'(a,a)')        '    short_name:           ',   &
     &                                  TRIM(S(i)%short_name)
          WRITE (stdout,'(a,a)')        '    data_netcdf_variable: ',   &
     &                                  TRIM(S(i)%data_netcdf_vname)
          WRITE (stdout,'(a,a)')        '    data_netcdf_time:     ',   &
     &                                  TRIM(S(i)%data_netcdf_tname)
          WRITE (stdout,'(a,a)')        '    source_units:         ',   &
     &                                  TRIM(S(i)%source_units)
          WRITE (stdout,'(a,a)')        '    destination_units:    ',   &
     &                                  TRIM(S(i)%destination_units)
          WRITE (stdout,'(a,a)')        '    source_grid:          ',   &
     &                                  TRIM(S(i)%source_grid)
          WRITE (stdout,'(a,a)')        '    destination_grid:     ',   &
     &                                  TRIM(S(i)%destination_grid)
          WRITE (stdout,'(a,1p,e15.8)') '    add_offset:           ',   &
     &                                  S(i)%add_offset
          WRITE (stdout,'(a,1p,e15.8)') '    scale:                ',   &
     &                                  S(i)%scale
          WRITE (stdout,'(a,l1)')       '    debug_write:          ',   &
     &                                  S(i)%debug_write
          WRITE (stdout,'(a,l1)')       '    connected:            ',   &
     &                                  S(i)%connected
          WRITE (stdout,'(a,a)')        '    connected_to:         ',   &
     &                                  TRIM(S(i)%connected_to)
          WRITE (stdout,'(a,a)')        '    regrid_method:        ',   &
     &                                  TRIM(S(i)%regrid_method)
          WRITE (stdout,'(a,a)')        '    extrapolate_method:   ',   &
     &                                  TRIM(S(i)%extrapolate_method)
        END DO
        FLUSH (stdout)
      END IF
!
      RETURN
  20  IF (Master) WRITE (stdout,50) TRIM(ADJUSTL(Cinfo(1))),            &
     &                              TRIM(filename)
!
  30  FORMAT (/,' COUPLING_METADATA - Unable to create YAML object',    &
     &        ' for ROMS I/O metadata file: ',/,21x,a,/,                &
     &        21x,'Default file is located in source directory.')
  40  FORMAT (/,' COUPLING_METADATA - Unable to open ROMS coupling',    &
     &        ' coupling file:',/,21x,a,/,21x,'ERROR: ',a,/,            &
     &        21x,'Default file is located in source directory.')
  50  FORMAT (/,' COUPLING_METADATA - Error while reading information', &
     &        'for metadata variable: ',a,/,21x,'File: ',a)
!
      END SUBROUTINE coupling_metadata
!
      FUNCTION io_metadata (FirstPass, Vinfo, scale, offset)            &
                    RESULT (Ldone)
!
!=======================================================================
!                                                                      !
!  It processes ROMS input/output fields metadata entry-by-entry from  !
!  'varinfo.dat'  or 'varinfo.yaml' dictionary.                        !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     FirsPass   Switch to initialize metadata processing (logical)    !
!                                                                      !
!     Vinfo      I/O Variable information (string array)               !
!                  Vinfo(1):  Field variable short-name                !
!                  Vinfo(2):  Long-name attribute                      !
!                  Vinfo(3):  Units attribute                          !
!                  Vinfo(4):  Field attribute                          !
!                  Vinfo(5):  Associated time variable name            !
!                  Vinfo(6):  Standard-name attribute                  !
!                  Vinfo(7):  Staggered C-grid variable type:          !
!                               'nulvar' => non-grided variable        !
!                               'p2dvar' => 2D PHI-variable            !
!                               'r2dvar' => 2D RHO-variable            !
!                               'u2dvar' => 2D U-variable              !
!                               'v2dvar' => 2D V-variable              !
!                               'p3dvar' => 3D PSI-variable            !
!                               'r3dvar' => 3D RHO-variable            !
!                               'u3dvar' => 3D U-variable              !
!                               'v3dvar' => 3D V-variable              !
!                               'w3dvar' => 3D W-variable              !
!                               'b3dvar' => 3D BED-sediment            !
!                               'l3dvar' => 3D spectral light variable !
!                               'l4dvar' => 4D spectral light variable !
!                  Vinfo(8):  Index code for information arrays        !
!                                                                      !
!     scale      Scale to convert input data to model units (real)     !
!                                                                      !
!     offeset    Value to add to input data (real)                     !
!                                                                      !
!     Ldone      True if end-of-file or dictionary found               !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      logical, intent(inout) :: FirstPass
!
      real(dp), intent(out) :: offset, scale
!
      character (len=*), intent(out) :: Vinfo(:)
!
!  Local variable declarations.
!
      logical, save :: IsDat     = .FALSE.
      logical, save :: IsYaml    = .FALSE.
      logical :: Ldone
!
      integer, parameter :: iunit = 10
      integer :: Idot, Lstr, Lvar
      integer :: i, j, io_err
!
      character (len=256) :: io_errmsg

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", io_metadata"
!
!-----------------------------------------------------------------------
!  On first pass, initialize metadata processing.
!-----------------------------------------------------------------------
!
!  Initialize.
!
      Ldone=.FALSE.
!
!  Determine metadata file extension: 'varinfo.dat'  or
!                                     'varinfo.yaml'
!
      IF (FirstPass) THEN
        FirstPass=.FALSE.
!
        Lstr=LEN_TRIM(varname)
        Idot=INDEX(varname(1:Lstr), CHAR(46), BACK=.TRUE.)
        SELECT CASE (lowercase(varname(Idot+1:Lstr)))
          CASE ('dat')
            IsDat=.TRUE.
          CASE ('yaml', 'yml')
            IsYaml=.TRUE.
        END SELECT
!
!  If YAML metadata, create dictionary and extract values.
!
        IF (IsYaml) THEN
          Ientry=0
!
          IF (FoundError(yaml_initialize(YML, TRIM(varname),            &
     &                                   LreportYAML),                  &
     &                   NoError, __LINE__, MyFile)) THEN
            Ldone=.TRUE.
            IF (Master) WRITE (stdout,30) TRIM(varname)
            RETURN
          END IF
!
          IF (FoundError(yaml_get(YML, 'metadata.variable',             &
     &                            Ystring1),                            &
     &                   NoError, __LINE__, MyFile)) RETURN
          Nentries=SIZE(Ystring1, DIM=1)
!
          IF (FoundError(yaml_get(YML, 'metadata.long_name',            &
     &                            Ystring2),                            &
     &                   NoError, __LINE__, MyFile)) RETURN
!
          IF (FoundError(yaml_get(YML, 'metadata.units',                &
     &                            Ystring3),                            &
     &                   NoError, __LINE__, MyFile)) RETURN
!
          IF (FoundError(yaml_get(YML, 'metadata.field',                &
     &                            Ystring4),                            &
     &                   NoError, __LINE__, MyFile)) RETURN
!
          IF (FoundError(yaml_get(YML, 'metadata.time',                 &
     &                            Ystring5),                            &
     &                   NoError, __LINE__, MyFile)) RETURN
!
          IF (FoundError(yaml_get(YML, 'metadata.standard_name',        &
     &                            Ystring6),                            &
     &                   NoError, __LINE__, MyFile)) RETURN
!
          IF (FoundError(yaml_get(YML, 'metadata.type',                 &
     &                            Ystring7),                            &
     &                   NoError, __LINE__, MyFile)) RETURN
!
          IF (FoundError(yaml_get(YML, 'metadata.index_code',           &
     &                            Ystring8),                            &
     &                   NoError, __LINE__, MyFile)) RETURN
!
          IF (allocated(Yreal1)) deallocate (Yreal1)
          allocate ( Yreal1(Nentries) )
          IF (FoundError(yaml_get(YML, 'metadata.add_offset',           &
     &                            Yreal1),                              &
     &                   NoError, __LINE__, MyFile)) RETURN
!
          IF (allocated(Yreal2)) deallocate (Yreal2)
          allocate ( Yreal2(Nentries) )
          IF (FoundError(yaml_get(YML, 'metadata.scale',                &
                                  Yreal2),                              &
     &                   NoError, __LINE__, MyFile)) RETURN
!
!  Otherwise, open deprecated 'varinfo.dat' file.
!
        ELSE IF (IsDat) THEN
          OPEN (UNIT=iunit, FILE=TRIM(varname), FORM='formatted',       &
     &          STATUS='old', IOSTAT=io_err, IOMSG=io_errmsg)
          IF (FoundError(io_err, NoError, __LINE__, MyFile)) THEN
            exit_flag=5
            Ldone=.TRUE.
            IF (Master) WRITE(stdout,40) TRIM(varname), TRIM(io_errmsg)
            RETURN
          END IF
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Process metadata entries.
!-----------------------------------------------------------------------
!
      DO j=1,SIZE(Vinfo)
        DO i=1,LEN(Vinfo(1))
          Vinfo(j)(i:i)=CHAR(32)
        END DO
      END DO
!
!  Extract metadata information from YAML structures.
!
      IF (IsYaml) THEN
        Ientry=Ientry+1                    ! advance variable counter
        IF (Ientry.le.Nentries) THEN
          Vinfo(1)=Ystring1(Ientry)%value  ! 'variable'      key
          Vinfo(2)=Ystring2(Ientry)%value  ! 'long_name'     key
          Vinfo(3)=Ystring3(Ientry)%value  ! 'units'         key
          Vinfo(4)=Ystring4(Ientry)%value  ! 'field'         key
          Vinfo(5)=Ystring5(Ientry)%value  ! 'time'          key
          Vinfo(6)=Ystring6(Ientry)%value  ! 'standard_name' key
          Vinfo(7)=Ystring7(Ientry)%value  ! 'type'          key
          Vinfo(8)=Ystring8(Ientry)%value  ! 'index_code'    key
          offset  =Yreal1(Ientry)          ! 'add_offset'    key
          scale   =Yreal2(Ientry)          ! 'scale'         key
        ELSE
          Ldone=.TRUE.
          CALL YML%destroy ()
          IF (allocated(Ystring1))  deallocate (Ystring1)
          IF (allocated(Ystring2))  deallocate (Ystring2)
          IF (allocated(Ystring3))  deallocate (Ystring3)
          IF (allocated(Ystring4))  deallocate (Ystring4)
          IF (allocated(Ystring5))  deallocate (Ystring5)
          IF (allocated(Ystring6))  deallocate (Ystring6)
          IF (allocated(Ystring7))  deallocate (Ystring7)
          IF (allocated(Ystring8))  deallocate (Ystring8)
          IF (allocated(Yreal1))    deallocate (Yreal1)
          IF (allocated(Yreal2))    deallocate (Yreal2)
          RETURN
        END IF
!
!  Otherwise, read in next metadata entry.  The 'standard_name' and
!  'add_offset' attributes are unavailable in 'varinfo.dat'.
!
      ELSE IF (IsDat) THEN
        DO WHILE (.TRUE.)
          READ (iunit,*,ERR=10,END=20) Vinfo(1)  ! variable
          Lvar=LEN_TRIM(Vinfo(1))
          IF ((Lvar.gt.0).and.(Vinfo(1)(1:1).ne.CHAR(33)).and.          &
     &                        (Vinfo(1)(1:1).ne.CHAR(36))) THEN
            READ (iunit,*,ERR=10) Vinfo(2)       ! long_name
            READ (iunit,*,ERR=10) Vinfo(3)       ! units
            READ (iunit,*,ERR=10) Vinfo(4)       ! field
            READ (iunit,*,ERR=10) Vinfo(5)       ! associated time
            Vinfo(6)='nulvar'                    ! standard_name
            READ (iunit,*,ERR=10) Vinfo(8)       ! index code
            READ (iunit,*,ERR=10) Vinfo(7)       ! C-grid type
            READ (iunit,*,ERR=10) scale
            offset  =0.0_dp                      ! add_offset
            Ldone=.FALSE.
            RETURN
          END IF
        END DO
  10    WRITE (stdout,50) TRIM(ADJUSTL(Vinfo(1)))
        STOP
  20    CLOSE (iunit)
        Ldone=.TRUE.
      END IF
!
  30  FORMAT (/,' IO_METADATA - Unable to create YAML object for',      &
     &        ' ROMS I/O metadata file: ',/,15x,a,/,                    &
     &        15x,'Default file is located in source directory.')
  40  FORMAT (/,' IO_METADATA - Unable to open ROMS I/O metadata ',     &
     &        'file:',/,15x,a,/,15x,'ERROR: ',a,/,                      &
     &        15x,'Default file is located in source directory.')
  50  FORMAT (/,' IO_METADATA - Error while reading information for ',  &
     &        'variable: ',a)
!
      RETURN
      END FUNCTION io_metadata
!
      FUNCTION metadata_has (S, short_name) RESULT (Findex)
!
!=======================================================================
!                                                                      !
!  It scans the fields metadata object (TYPE CouplingField) and        !
!  returns the index location in the block list of the requested       !
!  short-name keyword.                                                 !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     S            Fields metadata object (TYPE CouplingField)         !
!                                                                      !
!     short_name   Field short_name to find (string)                   !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     Findex       Index location in fields metadata list (integer)    !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      TYPE (CouplingField), allocatable, intent(in) :: S(:)
      character (len=*),     intent(in) :: short_name
!
!  Local variable declarations.
!
      integer :: Findex
      integer :: i
!
!-----------------------------------------------------------------------
!  Find index of specified field from names list.
!-----------------------------------------------------------------------
!
      Findex=-1
!
      DO i=1,SIZE(S)
        IF (S(i)%short_name.eq.short_name) THEN
          Findex=i
          EXIT
        END IF
      END DO
!
      RETURN
      END FUNCTION metadata_has
!
      END MODULE get_metadata_mod
