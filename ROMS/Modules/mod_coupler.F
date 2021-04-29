#include "cppdefs.h"
      MODULE mod_coupler

#if defined MODEL_COUPLING && defined MCT_LIB
!
!svn $Id: mod_coupler.F 1054 2021-03-06 19:47:12Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2021 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!
      USE mod_param
!
      implicit none
!
!-----------------------------------------------------------------------
!  Set several model coupling structures.
!-----------------------------------------------------------------------
!
!  Integer vector structure.
!
      TYPE T_INTEGER
        integer, pointer :: val(:)
      END TYPE
!
!  String vector structure.
!
      TYPE T_STRING
        character (len=10), pointer :: code(:)
      END TYPE
!
!  Coupling field attributes structure.
!
      TYPE T_FIELD
        integer :: FieldID              ! field ID
        integer :: GridID               ! associated grid ID
        integer :: MaskID               ! associated land/sea mask ID
        integer :: GridType             ! grid type (RHO-, U-, V-points)
        real(r8) :: AddOffset           ! number added to data
        real(r8) :: scale               ! scaling factor
        real(r8) :: ExpMin              ! exported minimum value
        real(r8) :: ExpMax              ! exported maximum value
        real(r8) :: ImpMin              ! imported minimum value
        real(r8) :: ImpMax              ! imported maximum value
        character (len=40) :: code      ! field code
        character (len=40) :: variable  ! associated model variable
        character (len=80) :: name      ! long descriptive field name
        character (len=80) :: units     ! field units
      END TYPE T_FIELD
!
!  Coupling model exchange mesh attributes structure.
!
      TYPE T_MESH
        integer :: GridID               ! grid ID
        integer :: GridType             ! grid type (RHO-, U-, V-points)
        character (len=40) :: code      ! grid code
        character (len=40) :: variable  ! associated model variable
        character (len=80) :: name      ! long descriptive grid name
        character (len=80) :: units     ! grid units
      END TYPE T_MESH
!
!  Time clock.
!
      TYPE T_CLOCK
        integer :: year                 ! year
        integer :: month                ! month
        integer :: day                  ! day
        integer :: hour                 ! hour
        integer :: minute               ! minute
        integer :: second               ! second
        integer :: YearDay              ! day of the year
        integer :: TimeZone             ! time zone, hours offset
        character (len=30) :: string    ! time string
      END TYPE T_CLOCK
!
!-----------------------------------------------------------------------
!  Set various variables used to couple ROMS/TOMS to other modeling
!  systems.
!-----------------------------------------------------------------------
!
!  Number of models to couple.
!
      integer :: Nmodels
!
!  Coupled model components IDs.
!
      integer :: ATMid = 2
      integer :: WAVid = 3
      integer :: OCNid = 1
!
!  Logical switch to report verbose import/export field ranges.
!
      logical :: Lreport = .FALSE.
!
!  Input coupled model order labels used to determine the values of
!  each model index in information variable.
!
      character (len=20), allocatable :: OrderLabel(:)
!
!  Coupled model indices. Values are initilized here to zero and
!  assigned in "inp_par" using order labels codes.
!
      integer :: Iatmos = 0            ! atmospheric model
      integer :: Iocean = 0            ! ocean model
      integer :: Iwaves = 0            ! wave model
!
!  Standard input file name for each coupled model.
!
      character (len=256), allocatable :: INPname(:)
!
!  Export/Import fields information file name.
!
      character (len=256) :: CPLname
!
!  Number of parallel nodes assigned to each model in the coupled model.
!  Their sum must be equal to the total number of processors.
!
      integer, allocatable :: Nthreads(:)
!
!  Assigned Pertsistent Execution Threads (PETs) for each coupled
!  model.
!
      TYPE (T_INTEGER), allocatable :: pets(:)
!
!  Time interval (seconds) between coupling of models.  This is a symmetric
!  matrix.  For example, the time interval coupling between ocean and
!  atmosphere models is:
!
!     TimeInterval(Iocean,Iatmos) = TimeInterval(Iocean,Iatmos)
!
      real(r8), allocatable :: TimeInterval(:,:)
!
!  Number of time-steps for how often to couple ROMS to other models.
!
!     CoupleSteps(:,ng) = MAX(1,INT(TimeInterval(Iocean,:)/dt(ng)))
!
      integer, allocatable :: CoupleSteps(:,:)
!
!  Export/Import fields information structure.  This information is read
!  from input CPLname file.
!
      integer, parameter :: MaxNumberFields = 50

      TYPE (T_FIELD) :: Fields(MaxNumberFields)
!
!  Number export and import fields for each coupled model.
!
      integer, allocatable :: Nexport(:)
      integer, allocatable :: Nimport(:)
!
!  Export/import fields IDs for each coupled model.
!
      TYPE (T_INTEGER), allocatable :: ExportID(:)
      TYPE (T_INTEGER), allocatable :: ImportID(:)
!
!  Export/import fields codes for each coupled model.
!
      TYPE (T_STRING), allocatable :: Export(:)
      TYPE (T_STRING), allocatable :: Import(:)
!
!  Export fields attribute string.
!
      character (len=240), allocatable :: ExportList(:)

      CONTAINS

      SUBROUTINE allocate_coupler (Nnodes)
!
!=======================================================================
!                                                                      !
!  This routine allocates all variables in the module for all coupled  !
!  models. It also initialize variable when appropriate.               !
!                                                                      !
!=======================================================================
!
      USE mod_parallel
      USE mod_iounits
      USE mod_ncparam
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: Nnodes
!
!  Local variable declarations.
!
      logical :: load

      integer, parameter :: inp = 20

      integer :: Nfields
      integer :: gtype, i, id, ifield, inode, lvar, model, ng, varid

      real(r8) :: add_offset, scale

      character (len=40 ) :: code
      character (len=80 ) :: string
      character (len=120), dimension(5) :: Vinfo
!
!-----------------------------------------------------------------------
!  Read in coupling export/import variable information.
!-----------------------------------------------------------------------
!
!  Open input coupling variable information file.
!
      OPEN (inp, FILE=TRIM(CPLname), FORM='formatted', STATUS='old',    &
     &      ERR=10)
      GO TO 20
  10  IF (MyRank.eq.0) WRITE(stdout,50) TRIM(CPLname)
      STOP
  20  CONTINUE
!
!  Read in variable information.  Ignore blank and comment [char(33)=!]
!  input lines.
!
      varid=0
      DO WHILE (.TRUE.)
        READ (inp,*,ERR=30,END=40) Vinfo(1)
        Lvar=LEN_TRIM(Vinfo(1))
!
!  Read in other variable information.
!
        IF ((Lvar.gt.0).and.(Vinfo(1)(1:1).ne.CHAR(33))) THEN
          READ (inp,*,ERR=30) Vinfo(2)
          READ (inp,*,ERR=30) Vinfo(3)
          READ (inp,*,ERR=30) Vinfo(4)
          READ (inp,*,ERR=30) Vinfo(5)
          READ (inp,*,ERR=30) add_offset
          READ (inp,*,ERR=30) scale
!
!  Determine staggered C-grid variable.
!
          SELECT CASE (TRIM(ADJUSTL(Vinfo(5))))
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
#ifdef SPECTRAL_LIGHT
            CASE ('s3dvar')
              gtype=s3dvar
#endif
            CASE DEFAULT
              gtype=0
          END SELECT
!
!  Load variable data into information arrays.
!
          varid=varid+1
          IF (varid.gt.MaxNumberFields) THEN
            WRITE (stdout,60) MaxNumberFields, varid
            STOP
          END IF
          Fields(varid) % code      = TRIM(ADJUSTL(Vinfo(1)))
          Fields(varid) % variable  = TRIM(ADJUSTL(Vinfo(2)))
          Fields(varid) % name      = TRIM(ADJUSTL(Vinfo(3)))
          Fields(varid) % units     = TRIM(ADJUSTL(Vinfo(4)))
          Fields(varid) % FieldID   = varid
          Fields(varid) % GridType  = gtype
          Fields(varid) % AddOffset = add_offset
          Fields(varid) % scale     = scale
        END IF
      END DO
      GO TO 40
  30  WRITE (stdout,80) TRIM(ADJUSTL(Vinfo(1)))
      STOP
  40  CLOSE (inp)
      Nfields=varid
!
!-----------------------------------------------------------------------
!  Determine identification index for export and import fields.
!-----------------------------------------------------------------------
!
!  Allocate IDs structures.
!
      IF (.not.allocated(ExportID)) THEN
        allocate ( ExportID(Nmodels) )
        DO model=1,Nmodels
          allocate ( ExportID(model)%val(Nexport(model)) )
          ExportID(model)%val=0
        END DO
      END IF
      IF (.not.allocated(ImportID)) THEN
        allocate ( ImportID(Nmodels) )
        DO model=1,Nmodels
          allocate ( ImportID(model)%val(Nimport(model)) )
          ImportID(model)%val=0
        END DO
      END IF
      IF (.not.allocated(ExportList)) THEN
        allocate ( ExportList(Nmodels) )
      END IF
!
!  Look fields information and extract Export/Import fields IDs for
!  each coupled model.
!
      DO model=1,Nmodels
        DO ifield=1,Nexport(model)
          DO i=1,Nfields
            IF (TRIM(ADJUSTL(Fields(i)%code)).eq.                       &
     &          TRIM(ADJUSTL(Export(model)%code(ifield)))) THEN
              ExportID(model)%val(ifield)=Fields(i)%FieldID
            END IF
          END DO
        END DO
        DO ifield=1,Nimport(model)
          DO i=1,Nfields
            IF (TRIM(ADJUSTL(Fields(i)%code)).eq.                       &
     &          TRIM(ADJUSTL(Import(model)%code(ifield)))) THEN
              ImportID(model)%val(ifield)=Fields(i)%FieldID
            END IF
          END DO
        END DO
      END DO
      DO model=1,Nmodels
        ExportList(model)=' '
        DO ifield=1,Nexport(model)
          id=ExportID(model)%val(ifield)
          IF (id.gt.0) THEN
            code=ADJUSTL(Fields(id)%code)
            IF (ifield.eq.1) THEN
              ExportList(model)=TRIM(ExportList(model))//TRIM(code)
            ELSE
              ExportList(model)=TRIM(ExportList(model))//':'//TRIM(code)
            END IF
          ELSE
            WRITE (stdout,70) model, TRIM(ExportList(model)),           &
     &                        TRIM(CPLname)
            STOP
          END IF
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Assign processors to coupled models.
!-----------------------------------------------------------------------
!
!  Allocate structure.
!
      IF (.not.allocated(pets)) THEN
        allocate ( pets(Nmodels) )
        DO model=1,Nmodels
          allocate ( pets(model)%val(Nthreads(model)) )
        END DO
      END IF
!
!  Assign parallel threads for each coupled model.  Start counting from
!  zero. That is, they are [0:Nnodes-1] available.
!
      inode=-1
      DO model=1,Nmodels
        DO i=1,Nthreads(model)
          inode=inode+1
          pets(model)%val(i)=inode
        END DO
      END DO
!
!  Report.
!
      IF ((inode+1).ne.Nnodes) THEN
        IF (MyRank.eq.0) THEN
          WRITE (stdout,80) inode, Nnodes
        END IF
        STOP
      ELSE
        IF (MyRank.eq.0) THEN
          WRITE (stdout,90)
          DO model=1,Nmodels
            IF (model.eq.Iocean) THEN
              string='Ocean Model MPI nodes:'
            ELSE IF (model.eq.Iwaves) THEN
              string='Waves Model MPI nodes:'
            ELSE IF (model.eq.Iatmos) THEN
              string='Atmos Model MPI nodes:'
            END IF
            WRITE (stdout,100) TRIM(string),                            &
     &                         pets(model)%val(1),                      &
     &                         pets(model)%val(Nthreads(model))
          END DO
        END IF
        WRITE (stdout,'(/)')
      END IF
      IF (MyRank.eq.0) THEN
        DO model=1,Nmodels
          IF (model.eq.Iocean) THEN
            string='Ocean Export:'
          ELSE IF (model.eq.Iwaves) THEN
            string='Waves Export:'
          ELSE IF (model.eq.Iatmos) THEN
            string='Atmos Export:'
          END IF
          WRITE (stdout,110) TRIM(string), TRIM(ExportList(model))
        END DO
        WRITE (stdout,'(/)')
      END IF

!
 50   FORMAT (/,' MOD_COUPLER - Unable to open variable information',   &
     &        ' file: ',/,15x,a,/,15x,'Default file is located in',     &
     &        ' source directory.')
 60   FORMAT (/,' MOD_COUPLER - too small dimension ',                  &
     &        'parameter, MV = ',2i5,/,15x,                             &
     &        'change file  mod_ncparam.F  and recompile.')
 70   FORMAT (/,' MOD_COUPLER - Unregistered export field for ',        &
     &          ' model = ',i1,/,15x,'ExportList = ',a,/,15x,           &
     &          ' check file = ',a)
 80   FORMAT (/,' MOD_COUPLER - Number assigned processors: ',          &
     &        i3.3,/,15x,'not equal to spawned MPI nodes: ',i3.3)
 90   FORMAT (/,' Model Coupling Parallel Threads:',/)
100   FORMAT (3x,a,3x,i3.3,' - ',i3.3)
110   FORMAT (3x,a,1x,a)

      END SUBROUTINE allocate_coupler
!
#endif
      END MODULE mod_coupler
