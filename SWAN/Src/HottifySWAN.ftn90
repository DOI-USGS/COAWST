!     HottifySWAN
!     Author: Casey Dietrich
!     Date  : 09-28-2010
!
!     This program globalizes the original set of SWAN hotstart files
!     and then localizes them over a different number of cores.
!     This program is specifically meant for unstructured SWAN grids.
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 1993-2015  Delft University of Technology
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     A copy of the GNU General Public License is available at
!     http://www.gnu.org/copyleft/gpl.html#SEC3
!     or by writing to the Free Software Foundation, Inc.,
!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
MODULE DATA

      IMPLICIT NONE

      INTEGER,ALLOCATABLE :: G2L(:,:)
      INTEGER             :: NumCores
      INTEGER             :: NumDirs
      INTEGER             :: NumFreqs
      INTEGER             :: NumGlobalVerts
      INTEGER, PARAMETER  :: UnitNumber = 67

      REAL(8),ALLOCATABLE :: Lat(:)
      REAL(8),ALLOCATABLE :: Lon(:)

      CHARACTER(LEN=50)   :: HotFile

      LOGICAL             :: Free

      TYPE Connectivity
         INTEGER,ALLOCATABLE :: Conn(:)
      END TYPE
      TYPE(Connectivity),ALLOCATABLE :: L2G(:)

END MODULE



PROGRAM HottifySWAN

      IMPLICIT NONE

      CHARACTER(LEN=1) :: UserSel

      WRITE(*,'(A)',ADVANCE='YES') " "
      WRITE(*,'(A)',ADVANCE='YES') "This program will globalize or localize a set of SWAN hot-start files."
      WRITE(*,'(A)',ADVANCE='YES') "... Do you want to:"
      WRITE(*,'(A)',ADVANCE='YES') "...... 1. Create a global file from an existing set of local files."
      WRITE(*,'(A)',ADVANCE='YES') "...... 2. Create a set of local files from an existing global file."
      WRITE(*,'(A)',ADVANCE='NO') "... Please enter your selection (1/2): "
      READ(*,'(A)') UserSel

      IF(UserSel.EQ."1")THEN

         WRITE(*,'(A)',ADVANCE='YES') " "
         WRITE(*,'(A)',ADVANCE='YES') "You have chosen to GLOBALIZE."

         CALL GlobalToLocal
         CALL WriteHeader
         CALL Globalize

      ELSEIF(UserSel.EQ."2")THEN

         WRITE(*,'(A)',ADVANCE='YES') " "
         WRITE(*,'(A)',ADVANCE='YES') "You have chosen to LOCALIZE."

         CALL LocalToGlobal
         CALL Localize

      ELSE

         WRITE(*,'(A)',ADVANCE='YES') " "
         WRITE(*,'(A)',ADVANCE='YES') "Your selection is not valid."

      ENDIF

      WRITE(*,'(A)',ADVANCE='YES') " "

END PROGRAM



SUBROUTINE Globalize

      USE DATA

      IMPLICIT NONE

      CHARACTER(LEN=1)          :: JunkC
      CHARACTER(LEN=5*NumDirs)  :: Line
      CHARACTER(LEN=50)         :: SwanFile

      CHARACTER(LEN=16)         :: RPROJID
      CHARACTER(LEN=4)          :: RPROJNR
      CHARACTER(LEN=20)         :: RVERTXT
      CHARACTER(LEN=20)         :: RCHTIME

      INTEGER                   :: IVAL
      INTEGER                   :: IVAL2
      REAL                      :: RVAL
      REAL                      :: RVAL2

      INTEGER                   :: Core
      INTEGER                   :: IC
      INTEGER                   :: ID
      INTEGER                   :: IS
      INTEGER                   :: IV
      INTEGER                   :: LocalVert
      INTEGER                   :: NumLocalVerts

      LOGICAL                   :: Nonstat

      TYPE ActionDensity
         CHARACTER(LEN=6) :: Type
         CHARACTER(LEN=24) :: Factor
         INTEGER,ALLOCATABLE :: Values(:,:)
      END TYPE
      TYPE LocalInfo
         TYPE(ActionDensity),ALLOCATABLE :: Ac(:)
      END TYPE
      TYPE(LocalInfo),ALLOCATABLE :: Local(:)

      WRITE(*,'(A)',ADVANCE='YES') " "
      WRITE(*,'(A)',ADVANCE='YES') "Globalizing the local SWAN hot-start files ..."
      WRITE(*,'(A,A,A)',ADVANCE='YES') "... Building the array with the local action densities."
      WRITE(*,'(A)',ADVANCE='NO') "... Progress: +"
      FLUSH(6)

      Nonstat = .false.

      ALLOCATE(Local(1:NumCores))

      DO IC=1,NumCores

         CALL ShowBar(IC,NumCores)

         WRITE(SwanFile,'(A,I4.4,A1,A)') "PE",IC-1,"/",TRIM(HotFile)
         IF (Free) THEN
            OPEN(UNIT=99,FILE=TRIM(SwanFile),ACTION='READ')
         ELSE
            OPEN(UNIT=99,FILE=TRIM(SwanFile),FORM='UNFORMATTED',ACTION='READ')
         ENDIF

         IF (Free) THEN

            ! SWAN standard file, with version
            READ(99,'(A)') Line

            ! time and location
  10        READ(99,'(A)') Line
            IF (Line(1:4) == 'TIME') Nonstat=.true.
            IF (Line(1:6) /= 'LONLAT' .and. Line(1:4) /= 'LOCA') GOTO 10

            ! number of locations
            READ(99,'(A)') Line
            READ(Line,*) NumLocalVerts
            DO IV=1,NumLocalVerts
               READ(99,'(A)') JunkC
            ENDDO

            ALLOCATE(Local(IC)%Ac(1:NumLocalVerts))

            ! relative frequencies in Hz
            READ(99,'(A)') Line

            ! number of frequencies
            READ(99,'(A)') Line
            DO IS=1,NumFreqs
               READ(99,'(A)') Line
            ENDDO

            ! spectral Cartesian directions in degr
            READ(99,'(A)') Line

            ! number of directions
            READ(99,'(A)') Line
            READ(Line,*) NumDirs
            DO ID=1,NumDirs
               READ(99,'(A)') Line
            ENDDO

            ! QUANT
            READ(99,'(A)') Line

            ! number of quantities in table
            READ(99,'(A)') Line

            ! action densities
            READ(99,'(A)') Line

            ! unit
            READ(99,'(A)') Line

            ! exception value
            READ(99,'(A)') Line

            ! date and time
            IF (Nonstat) READ(99,'(A)') Line

            DO IV=1,NumLocalVerts

               READ(99,'(A)') Line
               WRITE(Local(IC)%Ac(IV)%Type,'(A)') TRIM(Line)

               IF(INDEX(TRIM(Line),"FACTOR").GT.0)THEN

                 ALLOCATE(Local(IC)%Ac(IV)%Values(1:NumFreqs,1:NumDirs))

                 READ(99,'(A)') Line
                 WRITE(Local(IC)%Ac(IV)%Factor,'(A)') TRIM(Line)

                 DO IS=1,NumFreqs
                    READ(99,*) Local(IC)%Ac(IV)%Values(IS,:)
                 ENDDO

               ENDIF

            ENDDO

         ELSE

            READ (99) RVERTXT
            READ (99) RPROJID, RPROJNR
            READ (99) IVAL
            IF (IVAL==1) THEN
              READ (99) IVAL2
              Nonstat=.true.
            ENDIF
            READ (99) IVAL
            !
            READ (99) NumLocalVerts
            DO IV=1,NumLocalVerts
               READ (99) RVAL, RVAL2
            ENDDO
            ALLOCATE(Local(IC)%Ac(1:NumLocalVerts))
            !
            READ (99) IVAL
            DO IS=1,NumFreqs
               READ (99) RVAL
            ENDDO
            READ (99) IVAL
            DO ID=1,NumDirs
               READ (99) RVAL
            ENDDO
            !
            IF (Nonstat) READ(99) RCHTIME
            !
            DO IV=1,NumLocalVerts
               ALLOCATE(Local(IC)%Ac(IV)%Values(1:NumDirs,1:NumFreqs))
               READ (99) Local(IC)%Ac(IV)%Values(:,:)
            ENDDO

         ENDIF

         CLOSE(UNIT=99,STATUS='KEEP')

      ENDDO

      WRITE(*,'(A)',ADVANCE='YES')
      WRITE(*,'(A)',ADVANCE='YES') "... Done."

      WRITE(*,'(A)',ADVANCE='YES')
      WRITE(*,'(A)',ADVANCE='YES') "Writing the global SWAN hot-start file ..."

      WRITE(SwanFile,'(A)') TRIM(HotFile)
      IF (Free) THEN
         OPEN(UNIT=UnitNumber,FILE=TRIM(SwanFile),ACTION='WRITE',POSITION='APPEND')
      ELSE
         OPEN(UNIT=UnitNumber,FILE=TRIM(SwanFile),FORM='UNFORMATTED',ACTION='WRITE',POSITION='APPEND')
      ENDIF

      WRITE(*,'(A,A,A)',ADVANCE='YES') "... The hot-start information is being written to the ",TRIM(SwanFile)," file."
      WRITE(*,'(A)',ADVANCE='NO') "... Progress: +"
      FLUSH(6)

      IF (Free) THEN

         DO IV=1,NumGlobalVerts

            CALL ShowBar(IV,NumGlobalVerts)

            Core      = G2L(IV,1)
            LocalVert = G2L(IV,2)

            WRITE(UnitNumber,'(A)') TRIM(Local(Core)%Ac(LocalVert)%Type)

            IF(INDEX(Local(Core)%Ac(LocalVert)%Type,"FACTOR").GT.0)THEN

              WRITE(UnitNumber,'(A)') TRIM(Local(Core)%Ac(LocalVert)%Factor)
              DO IS=1,NumFreqs
                 WRITE(UnitNumber,'(200(1X,I4))') Local(Core)%Ac(LocalVert)%Values(IS,:)
              ENDDO

            ENDIF

         ENDDO

      ELSE

         DO IV=1,NumGlobalVerts

            CALL ShowBar(IV,NumGlobalVerts)

            Core      = G2L(IV,1)
            LocalVert = G2L(IV,2)

            WRITE(UnitNumber) Local(Core)%Ac(LocalVert)%Values(:,:)

         ENDDO

      ENDIF

      WRITE(*,'(A)',ADVANCE='YES')
      WRITE(*,'(A)',ADVANCE='YES') "... Done."

      CLOSE(UNIT=UnitNumber,STATUS='KEEP')

      RETURN

END SUBROUTINE



SUBROUTINE GlobalToLocal

      USE DATA

      IMPLICIT NONE

      CHARACTER(LEN=50) :: Fort18File
      CHARACTER(LEN=1)  :: JunkC

      INTEGER           :: Core
      INTEGER           :: IC
      INTEGER           :: IE
      INTEGER           :: IV
      INTEGER           :: JunkI
      INTEGER           :: NumLocalElems
      INTEGER           :: NumLocalVerts
      INTEGER           :: Vert

      WRITE(*,'(A)',ADVANCE='YES') " "
      WRITE(*,'(A)',ADVANCE='YES') "Developing the global-to-local connectivity ..."

      OPEN(UNIT=14,FILE='fort.14',ACTION='READ')

      READ(14,'(A)') JunkC
      READ(14,*) JunkI,NumGlobalVerts

      ALLOCATE(Lon(1:NumGlobalVerts))
      ALLOCATE(Lat(1:NumGlobalVerts))

      DO IV=1,NumGlobalVerts
         READ(14,*) JunkI,Lon(IV),Lat(IV)
      ENDDO

      CLOSE(UNIT=14,STATUS='KEEP')

      WRITE(*,'(A,I8.8,A)',ADVANCE='YES') "... The global mesh has ",NumGlobalVerts," vertices."

      ALLOCATE(G2L(1:NumGlobalVerts,1:2))

      OPEN(UNIT=99,FILE='partmesh.txt',ACTION='READ')

      NumCores = 0
      DO IV=1,NumGlobalVerts
         READ(99,*) Core
         IF(Core.GT.NumCores) NumCores = Core
         G2L(IV,1) = Core
      ENDDO

      CLOSE(UNIT=99,STATUS='KEEP')

      WRITE(*,'(A,I8.8,A)',ADVANCE='YES') "... It has been decomposed on ",NumCores," computational cores."
      WRITE(*,'(A)',ADVANCE='NO') "... Progress: +"

      DO IC=1,NumCores

         CALL ShowBar(IC,NumCores)

         WRITE(Fort18File,'(A,I4.4,A)') "PE",IC-1,"/fort.18"

         OPEN(UNIT=18,FILE=TRIM(Fort18File),ACTION='READ')

         READ(18,'(A)') JunkC
         READ(18,'(A8,I12,I12,I12)') JunkC,JunkI,JunkI,NumLocalElems

         DO IE=1,NumLocalElems
            READ(18,'(A)') JunkC
         ENDDO

         READ(18,'(A8,I12,I12,I12)') JunkC,JunkI,JunkI,NumLocalVerts

         DO IV=1,NumLocalVerts
            READ(18,'(I12)') Vert
            IF(Vert.GT.0)THEN
               G2L(Vert,2) = IV
            ENDIF
         ENDDO

         CLOSE(UNIT=18,STATUS='KEEP')

      ENDDO

      WRITE(*,'(A)',ADVANCE='YES')
      WRITE(*,'(A)',ADVANCE='YES') "... Done."

      RETURN

END SUBROUTINE



SUBROUTINE Localize

      USE DATA

      IMPLICIT NONE

      CHARACTER(LEN=1000) :: Line
      CHARACTER(LEN=50)   :: SwanFile
      CHARACTER(LEN=1)    :: UserSel

      CHARACTER(LEN=16)   :: RPROJID
      CHARACTER(LEN=4)    :: RPROJNR
      CHARACTER(LEN=20)   :: RVERTXT
      CHARACTER(LEN=20)   :: RCHTIME

      INTEGER             :: IVAL
      INTEGER             :: IVAL2
      REAL                :: RVAL
      REAL                :: RVAL2

      INTEGER             :: GlobalVert
      INTEGER             :: IC
      INTEGER             :: ID
      INTEGER             :: IS
      INTEGER             :: IV
      INTEGER             :: NumLocalVerts

      LOGICAL             :: Nonstat, Kspher

      TYPE ActionDensity
         CHARACTER(LEN=6) :: Type
         CHARACTER(LEN=24) :: Factor
         INTEGER,ALLOCATABLE :: Values(:,:)
      END TYPE
      TYPE(ActionDensity),ALLOCATABLE :: Ac(:)

      WRITE(*,'(A)',ADVANCE='YES')
      WRITE(*,'(A)',ADVANCE='YES') "Reading the global SWAN hot-start file ..."

      WRITE(*,'(A)',ADVANCE='NO') "... Please enter the name of the global hotfile: "
      READ(*,*) HotFile

      WRITE(*,'(A)',ADVANCE='YES') " "
      WRITE(*,'(A)',ADVANCE='YES') "Indicate the format of the hotfile ..."
      WRITE(*,'(A)',ADVANCE='YES') "...... F. The format of the hotfile is free, i.e. human readable."
      WRITE(*,'(A)',ADVANCE='YES') "...... U. The format of the hotfile is unformatted or binary."
      WRITE(*,'(A)',ADVANCE='NO') "... Please enter your selection (F/U): "
      READ(*,*) UserSel

      IF(UserSel.EQ."F" .OR. UserSel.EQ."f")THEN
         WRITE(*,'(A)',ADVANCE='YES') "The format of the hotfile is FREE."
         Free = .TRUE.
      ELSEIF(UserSel.EQ."U" .OR. UserSel.EQ."u")THEN
         WRITE(*,'(A)',ADVANCE='YES') "The format of the hotfile is UNFORMATTED."
         Free = .FALSE.
      ELSE
         WRITE(*,'(A)',ADVANCE='YES') "Your selection is not valid."
         WRITE(*,'(A)',ADVANCE='YES') "The format of the hotfile is assumed to be FREE."
         Free = .TRUE.
      ENDIF
      WRITE(*,'(A)',ADVANCE='YES') " "

      WRITE(SwanFile,'(A)') TRIM(HotFile)

      WRITE(*,'(A,A,A)',ADVANCE='YES') "... Now processing the ",TRIM(SwanFile)," file."
      WRITE(*,'(A)',ADVANCE='NO') "... Progress: +"

      Nonstat = .false.

      IF (Free) THEN
         OPEN(UNIT=UnitNumber,FILE=TRIM(SwanFile),ACTION='READ')
      ELSE
         OPEN(UNIT=UnitNumber,FILE=TRIM(SwanFile),FORM='UNFORMATTED',ACTION='READ')
      ENDIF

      IF (Free) THEN

         ! SWAN standard file, with version
         READ(UnitNumber,'(A)') Line

         ! time and location
  10     READ(UnitNumber,'(A)') Line
         IF (Line(1:4) == 'TIME') Nonstat=.true.
         IF (Line(1:6) /= 'LONLAT' .and. Line(1:4) /= 'LOCA') GOTO 10

         ! number of locations
         READ(UnitNumber,'(A)') Line
         READ(Line,*) NumGlobalVerts
         DO IV=1,NumGlobalVerts
            READ(UnitNumber,'(A)') Line
         ENDDO

         ALLOCATE(Ac(1:NumGlobalVerts))

         ! relative frequencies in Hz
         READ(UnitNumber,'(A)') Line

         ! number of frequencies
         READ(UnitNumber,'(A)') Line
         READ(Line,*) NumFreqs
         DO IS=1,NumFreqs
            READ(UnitNumber,'(A)') Line
         ENDDO

         ! spectral Cartesian directions in degr
         READ(UnitNumber,'(A)') Line

         ! number of directions
         READ(UnitNumber,'(A)') Line
         READ(Line,*) NumDirs
         DO ID=1,NumDirs
            READ(UnitNumber,'(A)') Line
         ENDDO

         ! QUANT
         READ(UnitNumber,'(A)') Line

         ! number of quantities in table
         READ(UnitNumber,'(A)') Line

         ! action densities
         READ(UnitNumber,'(A)') Line

         ! unit
         READ(UnitNumber,'(A)') Line

         ! exception value
         READ(UnitNumber,'(A)') Line

         ! date and time
         IF (Nonstat) READ(UnitNumber,'(A)') Line

         DO IV=1,NumGlobalVerts

            CALL ShowBar(IV,NumGlobalVerts)

            READ(UnitNumber,'(A)') Line
            WRITE(Ac(IV)%Type,'(A)') TRIM(Line)

            IF(INDEX(TRIM(Line),"FACTOR").GT.0)THEN

              ALLOCATE(Ac(IV)%Values(1:NumFreqs,1:NumDirs))

              READ(UnitNumber,'(A)') Line
              WRITE(Ac(IV)%Factor,'(A)') TRIM(Line)

              DO IS=1,NumFreqs
                 READ(UnitNumber,*) Ac(IV)%Values(IS,:)
              ENDDO

           ENDIF

         ENDDO

      ELSE

         READ (UnitNumber) RVERTXT
         READ (UnitNumber) RPROJID, RPROJNR
         READ (UnitNumber) IVAL
         IF (IVAL==1) THEN
            READ (UnitNumber) IVAL2
            Nonstat=.true.
         ENDIF
         READ (UnitNumber) IVAL
         !
         READ (UnitNumber) NumGlobalVerts
         DO IV=1,NumGlobalVerts
            READ (UnitNumber) RVAL, RVAL2
         ENDDO
         ALLOCATE(Ac(1:NumGlobalVerts))
         !
         READ (UnitNumber) NumFreqs
         DO IS=1,NumFreqs
            READ (UnitNumber) RVAL
         ENDDO
         READ (UnitNumber) NumDirs
         DO ID=1,NumDirs
            READ (UnitNumber) RVAL
         ENDDO
         !
         IF (Nonstat) READ(UnitNumber) RCHTIME
         !
         DO IV=1,NumGlobalVerts
            CALL ShowBar(IV,NumGlobalVerts)
            ALLOCATE(Ac(IV)%Values(1:NumDirs,1:NumFreqs))
            READ(UnitNumber) Ac(IV)%Values(:,:)
         ENDDO

      ENDIF

      CLOSE(UNIT=UnitNumber,STATUS='KEEP')

      WRITE(*,'(A)',ADVANCE='YES')
      WRITE(*,'(A)',ADVANCE='YES') "... Done."

      WRITE(*,'(A)',ADVANCE='YES')
      WRITE(*,'(A)',ADVANCE='YES') "Writing the local SWAN hot-start files ..."
      WRITE(*,'(A)',ADVANCE='NO') "... Progress: +"

      DO IC=1,NumCores

         CALL ShowBar(IC,NumCores)

         WRITE(SwanFile,'(A)') TRIM(HotFile)
         IF (Free) THEN
            OPEN(UNIT=UnitNumber,FILE=TRIM(SwanFile),ACTION='READ')
         ELSE
            OPEN(UNIT=UnitNumber,FILE=TRIM(SwanFile),FORM='UNFORMATTED',ACTION='READ')
         ENDIF

         WRITE(SwanFile,'(A,I4.4,A1,A)') "PE",IC-1,"/",TRIM(HotFile)
         IF (Free) THEN
            OPEN(UNIT=99,FILE=TRIM(SwanFile),ACTION='WRITE')
         ELSE
            OPEN(UNIT=99,FILE=TRIM(SwanFile),FORM='UNFORMATTED',ACTION='WRITE')
         ENDIF

         IF (Free) THEN

            ! SWAN standard file, with version
            READ(UnitNumber,'(A)') Line
            WRITE(99,'(A)') TRIM(Line)

            ! time and location
  20        READ(UnitNumber,'(A)') Line
            WRITE(99,'(A)') TRIM(Line)
            IF (Line(1:4) == 'TIME') Nonstat=.true.
            IF (Line(1:6) /= 'LONLAT' .and. Line(1:4) /= 'LOCA') GOTO 20
            IF (Line(1:6) == 'LONLAT') THEN
               Kspher = .true.
            ELSEIF (Line(1:4) == 'LOCA') THEN
               Kspher = .false.
            ENDIF

            ! number of locations
            READ(UnitNumber,'(A)') Line
            DO IV=1,NumGlobalVerts
               READ(UnitNumber,'(A)') Line
            ENDDO
            NumLocalVerts = SIZE(L2G(IC)%Conn)
            WRITE(99,*) NumLocalVerts
            IF (Kspher) THEN
               DO IV=1,NumLocalVerts
                  GlobalVert = ABS(L2G(IC)%Conn(IV))
                  WRITE(99,'(F12.6,F12.6)') Lon(GlobalVert),Lat(GlobalVert)
               ENDDO
            ELSE
               DO IV=1,NumLocalVerts
                  GlobalVert = ABS(L2G(IC)%Conn(IV))
                  WRITE(99,'(F14.4,F14.4)') Lon(GlobalVert),Lat(GlobalVert)
               ENDDO
            ENDIF

            ! relative frequencies in Hz
            READ(UnitNumber,'(A)') Line
            WRITE(99,'(A)') TRIM(Line)

            ! number of frequencies
            READ(UnitNumber,'(A)') Line
            WRITE(99,'(A)') TRIM(Line)
            READ(Line,*) NumFreqs
            DO IS=1,NumFreqs
               READ(UnitNumber,'(A)') Line
               WRITE(99,'(A)') TRIM(Line)
            ENDDO

            ! spectral Cartesian directions in degr
            READ(UnitNumber,'(A)') Line
            WRITE(99,'(A)') TRIM(Line)

            ! number of directions
            READ(UnitNumber,'(A)') Line
            WRITE(99,'(A)') TRIM(Line)
            READ(Line,*) NumDirs
            DO ID=1,NumDirs
               READ(UnitNumber,'(A)') Line
               WRITE(99,'(A)') TRIM(Line)
            ENDDO

            ! QUANT
            READ(UnitNumber,'(A)') Line
            WRITE(99,'(A)') TRIM(Line)

            ! number of quantities in table
            READ(UnitNumber,'(A)') Line
            WRITE(99,'(A)') TRIM(Line)

            ! action densities
            READ(UnitNumber,'(A)') Line
            WRITE(99,'(A)') TRIM(Line)

            ! unit
            READ(UnitNumber,'(A)') Line
            WRITE(99,'(A)') TRIM(Line)

            ! exception value
            READ(UnitNumber,'(A)') Line
            WRITE(99,'(A)') TRIM(Line)

            ! date and time
            IF (Nonstat) THEN
               READ(UnitNumber,'(A)') Line
               WRITE(99,'(A)') TRIM(Line)
            ENDIF

         ELSE

            READ (UnitNumber) RVERTXT
            WRITE (99) RVERTXT
            READ (UnitNumber) RPROJID, RPROJNR
            WRITE (99) RPROJID, RPROJNR
            READ (UnitNumber) IVAL
            WRITE (99) IVAL
            IF (IVAL==1) THEN
              READ (UnitNumber) IVAL2
              WRITE (99) IVAL2
              Nonstat=.true.
            ENDIF
            READ (UnitNumber) IVAL
            WRITE (99) IVAL
            !
            READ (UnitNumber) IVAL
            DO IV=1,NumGlobalVerts
               READ (UnitNumber) RVAL, RVAL2
            ENDDO
            NumLocalVerts = SIZE(L2G(IC)%Conn)
            WRITE(99) NumLocalVerts
            DO IV=1,NumLocalVerts
               GlobalVert = ABS(L2G(IC)%Conn(IV))
               WRITE(99) Lon(GlobalVert),Lat(GlobalVert)
            ENDDO
            !
            READ (UnitNumber) NumFreqs
            WRITE (99) NumFreqs
            DO IS=1,NumFreqs
               READ (UnitNumber) RVAL
               WRITE (99) RVAL
            ENDDO
            READ (UnitNumber) NumDirs
            WRITE (99) NumDirs
            DO ID=1,NumDirs
               READ (UnitNumber) RVAL
               WRITE (99) RVAL
            ENDDO
            !
            IF (Nonstat) THEN
               READ(UnitNumber) RCHTIME
               WRITE(99) RCHTIME
            ENDIF

         ENDIF

         CLOSE(UNIT=UnitNumber,STATUS='KEEP')

         IF (Free) THEN

            DO IV=1,NumLocalVerts

               GlobalVert = ABS(L2G(IC)%Conn(IV))

               WRITE(99,'(A)') TRIM(Ac(GlobalVert)%Type)

               IF(INDEX(Ac(GlobalVert)%Type,"FACTOR").GT.0)THEN

                 WRITE(99,'(A)') TRIM(Ac(GlobalVert)%Factor)
                 DO IS=1,NumFreqs
                    WRITE(99,'(200(1X,I4))') Ac(GlobalVert)%Values(IS,:)
                 ENDDO

               ENDIF

            ENDDO

         ELSE

            DO IV=1,NumLocalVerts
               GlobalVert = ABS(L2G(IC)%Conn(IV))
               WRITE(99) Ac(GlobalVert)%Values(:,:)
            ENDDO

         ENDIF

         CLOSE(UNIT=99,STATUS='KEEP')

      ENDDO

      WRITE(*,'(A)',ADVANCE='YES')
      WRITE(*,'(A)',ADVANCE='YES') "... Done."

      RETURN

END SUBROUTINE



SUBROUTINE LocalToGlobal

      USE DATA

      IMPLICIT NONE

      CHARACTER(LEN=50) :: Fort18File
      CHARACTER(LEN=1)  :: JunkC

      INTEGER           :: Core
      INTEGER           :: IC
      INTEGER           :: IE
      INTEGER           :: IV
      INTEGER           :: JunkI
      INTEGER           :: NumLocalElems
      INTEGER           :: NumLocalVerts
      INTEGER           :: Vert

      WRITE(*,'(A)',ADVANCE='YES') " "
      WRITE(*,'(A)',ADVANCE='YES') "Developing the local-to-global connectivity ..."

      OPEN(UNIT=14,FILE='fort.14',ACTION='READ')

      READ(14,'(A)') JunkC
      READ(14,*) JunkI,NumGlobalVerts

      ALLOCATE(Lon(1:NumGlobalVerts))
      ALLOCATE(Lat(1:NumGlobalVerts))

      DO IV=1,NumGlobalVerts
         READ(14,*) JunkI,Lon(IV),Lat(IV)
      ENDDO

      CLOSE(UNIT=14,STATUS='KEEP')

      WRITE(*,'(A,I8.8,A)',ADVANCE='YES') "... The global mesh has ",NumGlobalVerts," vertices."

      OPEN(UNIT=99,FILE='partmesh.txt',ACTION='READ')

      NumCores = 0
      DO IV=1,NumGlobalVerts
         READ(99,*) Core
         IF(Core.GT.NumCores) NumCores = Core
      ENDDO

      CLOSE(UNIT=99,STATUS='KEEP')

      ALLOCATE(L2G(1:NumCores))

      WRITE(*,'(A,I8.8,A)',ADVANCE='YES') "... It has been decomposed on ",NumCores," computational cores."
      WRITE(*,'(A)',ADVANCE='NO') "... Progress: +"

      DO IC=1,NumCores

         CALL ShowBar(IC,NumCores)

         WRITE(Fort18File,'(A,I4.4,A)') "PE",IC-1,"/fort.18"

         OPEN(UNIT=18,FILE=TRIM(Fort18File),ACTION='READ')

         READ(18,'(A)') JunkC
         READ(18,'(A8,I12,I12,I12)') JunkC,JunkI,JunkI,NumLocalElems

         DO IE=1,NumLocalElems
            READ(18,'(A)') JunkC
         ENDDO

         READ(18,'(A8,I12,I12,I12)') JunkC,JunkI,JunkI,NumLocalVerts

         ALLOCATE(L2G(IC)%Conn(1:NumLocalVerts))

         DO IV=1,NumLocalVerts
            READ(18,*) Vert
            L2G(IC)%Conn(IV) = Vert
         ENDDO

         CLOSE(UNIT=18,STATUS='KEEP')

      ENDDO

      WRITE(*,'(A)',ADVANCE='YES')
      WRITE(*,'(A)',ADVANCE='YES') "... Done."

      RETURN

END SUBROUTINE



SUBROUTINE WriteHeader

      USE DATA

      IMPLICIT NONE

      CHARACTER(LEN=15)   :: CurrentTime
      CHARACTER(LEN=1)    :: JunkC
      CHARACTER(LEN=1000) :: Line
      CHARACTER(LEN=50)   :: SwanFile
      CHARACTER(LEN=1)    :: UserSel

      CHARACTER(LEN=16)   :: RPROJID
      CHARACTER(LEN=4)    :: RPROJNR
      CHARACTER(LEN=20)   :: RVERTXT
      CHARACTER(LEN=20)   :: RCHTIME

      INTEGER             :: IVAL
      INTEGER             :: IVAL2
      REAL                :: RVAL
      REAL                :: RVAL2

      INTEGER             :: IC
      INTEGER             :: ID
      INTEGER             :: IS
      INTEGER             :: IV
      INTEGER             :: NumLocalVerts

      LOGICAL             :: Nonstat, Kspher

      Nonstat = .false.

      WRITE(*,'(A)',ADVANCE='YES') " "
      WRITE(*,'(A)',ADVANCE='YES') "Writing the header information ..."
      WRITE(*,'(A)',ADVANCE='NO') "... Enter the name of the hotfile: "
      READ(*,*) HotFile

      WRITE(*,'(A)',ADVANCE='YES') " "
      WRITE(*,'(A)',ADVANCE='YES') "Indicate the format of the hotfile ..."
      WRITE(*,'(A)',ADVANCE='YES') "...... F. The format of the hotfile is free, i.e. human readable."
      WRITE(*,'(A)',ADVANCE='YES') "...... U. The format of the hotfile is unformatted or binary."
      WRITE(*,'(A)',ADVANCE='NO') "... Please enter your selection (F/U): "
      READ(*,*) UserSel

      IF(UserSel.EQ."F" .OR. UserSel.EQ."f")THEN
         WRITE(*,'(A)',ADVANCE='YES') "The format of the hotfile is FREE."
         Free = .TRUE.
      ELSEIF(UserSel.EQ."U" .OR. UserSel.EQ."u")THEN
         WRITE(*,'(A)',ADVANCE='YES') "The format of the hotfile is UNFORMATTED."
         Free = .FALSE.
      ELSE
         WRITE(*,'(A)',ADVANCE='YES') "Your selection is not valid."
         WRITE(*,'(A)',ADVANCE='YES') "The format of the hotfile is assumed to be FREE."
         Free = .TRUE.
      ENDIF
      WRITE(*,'(A)',ADVANCE='YES') " "

      WRITE(SwanFile,'(A)') "PE0000/"//TRIM(HotFile)
      IF (Free) THEN
         OPEN(UNIT=99,FILE=TRIM(SwanFile),ACTION='READ')
      ELSE
         OPEN(UNIT=99,FILE=TRIM(SwanFile),FORM='UNFORMATTED',ACTION='READ')
      ENDIF

      WRITE(SwanFile,'(A)') TRIM(HotFile)
      IF (Free) THEN
         OPEN(UNIT=UnitNumber,FILE=TRIM(SwanFile),ACTION='WRITE')
      ELSE
         OPEN(UNIT=UnitNumber,FILE=TRIM(SwanFile),FORM='UNFORMATTED',ACTION='WRITE')
      ENDIF

      IF (Free) THEN

         ! SWAN standard file, with version
         READ(99,'(A)') Line
         WRITE(UnitNumber,'(A)') TRIM(Line)

         ! time and location
  10     READ(99,'(A)') Line
         WRITE(UnitNumber,'(A)') TRIM(Line)
         IF (Line(1:4) == 'TIME') Nonstat=.true.
         IF (Line(1:6) /= 'LONLAT' .and. Line(1:4) /= 'LOCA') GOTO 10
         IF (Line(1:6) == 'LONLAT') THEN
            Kspher = .true.
         ELSEIF (Line(1:4) == 'LOCA') THEN
            Kspher = .false.
         ENDIF

         ! number of locations
         READ(99,'(A)') Line
         READ(Line,*) NumLocalVerts
         DO IV=1,NumLocalVerts
            READ(99,'(A)') JunkC
         ENDDO
         WRITE(UnitNumber,*) NumGlobalVerts
         IF (Kspher) THEN
            DO IV=1,NumGlobalVerts
               WRITE(UnitNumber,'(F12.6,F12.6)') Lon(IV),Lat(IV)
            ENDDO
         ELSE
            DO IV=1,NumGlobalVerts
               WRITE(UnitNumber,'(F14.4,F14.4)') Lon(IV),Lat(IV)
            ENDDO
         ENDIF

         ! relative frequencies in Hz
         READ(99,'(A)') Line
         WRITE(UnitNumber,'(A)') TRIM(Line)

         ! number of frequencies
         READ(99,'(A)') Line
         WRITE(UnitNumber,'(A)') TRIM(Line)
         READ(Line,*) NumFreqs
         WRITE(*,'(A,I2.2,A)',ADVANCE='YES') "... These files contain ",NumFreqs," frequency bins."
         DO IS=1,NumFreqs
            READ(99,'(A)') Line
            WRITE(UnitNumber,'(A)') TRIM(Line)
         ENDDO

         ! spectral Cartesian directions in degr
         READ(99,'(A)') Line
         WRITE(UnitNumber,'(A)') TRIM(Line)

         ! number of directions
         READ(99,'(A)') Line
         WRITE(UnitNumber,'(A)') TRIM(Line)
         READ(Line,*) NumDirs
         WRITE(*,'(A,I2.2,A)',ADVANCE='YES') "... These files contain ",NumDirs," directional bins."
         DO ID=1,NumDirs
            READ(99,'(A)') Line
            WRITE(UnitNumber,'(A)') TRIM(Line)
         ENDDO

         ! QUANT
         READ(99,'(A)') Line
         WRITE(UnitNumber,'(A)') TRIM(Line)

         ! number of quantities in table
         READ(99,'(A)') Line
         WRITE(UnitNumber,'(A)') TRIM(Line)

         ! action densities
         READ(99,'(A)') Line
         WRITE(UnitNumber,'(A)') TRIM(Line)

         ! unit
         READ(99,'(A)') Line
         WRITE(UnitNumber,'(A)') TRIM(Line)

         ! exception value
         READ(99,'(A)') Line
         WRITE(UnitNumber,'(A)') TRIM(Line)

         ! date and time, if nonstationary
         IF (Nonstat) THEN
            READ(99,'(A)') Line
            WRITE(UnitNumber,'(A)') TRIM(Line)
            READ(Line,'(A15)') CurrentTime
            WRITE(*,'(A,A,A)',ADVANCE='YES') "... These files are time-stamped to ",TRIM(CurrentTime),"."
         ENDIF

      ELSE

         READ (99) RVERTXT
         WRITE (UnitNumber) RVERTXT
         READ (99) RPROJID, RPROJNR
         WRITE (UnitNumber) RPROJID, RPROJNR
         READ (99) IVAL
         WRITE (UnitNumber) IVAL
         IF (IVAL==1) THEN
           READ (99) IVAL2
           WRITE (UnitNumber) IVAL2
           Nonstat=.true.
         ENDIF
         READ (99) IVAL
         WRITE (UnitNumber) IVAL
         !
         READ (99) NumLocalVerts
         DO IV=1,NumLocalVerts
            READ (99) RVAL, RVAL2
         ENDDO
         WRITE (UnitNumber) NumGlobalVerts
         DO IV=1,NumGlobalVerts
            WRITE(UnitNumber) Lon(IV),Lat(IV)
         ENDDO
         !
         READ (99) NumFreqs
         WRITE (UnitNumber) NumFreqs
         WRITE(*,'(A,I2.2,A)',ADVANCE='YES') "... These files contain ",NumFreqs," frequency bins."
         DO IS=1,NumFreqs
            READ (99) RVAL
            WRITE (UnitNumber) RVAL
         ENDDO
         READ (99) NumDirs
         WRITE (UnitNumber) NumDirs
         WRITE(*,'(A,I2.2,A)',ADVANCE='YES') "... These files contain ",NumDirs," directional bins."
         DO ID=1,NumDirs
            READ (99) RVAL
            WRITE (UnitNumber) RVAL
         ENDDO
         !
         IF (Nonstat) THEN
            READ(99) RCHTIME
            WRITE(UnitNumber) RCHTIME
            WRITE(*,'(A,A,A)',ADVANCE='YES') "... These files are time-stamped to ",TRIM(RCHTIME),"."
         ENDIF

      ENDIF

      CLOSE(UNIT=99,STATUS='KEEP')

      CLOSE(UNIT=UnitNumber,STATUS='KEEP')

      WRITE(*,'(A)',ADVANCE='YES') "... Done."

      RETURN

END SUBROUTINE



SUBROUTINE ShowBar(Now,Total)

       IMPLICIT NONE

       INTEGER,INTENT(IN) :: Now
       INTEGER,INTENT(IN) :: Total

       INTEGER            :: N

       outer: DO N=1,20
          IF(Now.EQ.CEILING(N*0.05*REAL(Total)))THEN
             IF(MOD(N,5).EQ.0)THEN
                WRITE(*,'(A)',ADVANCE='NO') "+"
             ELSE
                WRITE(*,'(A)',ADVANCE='NO') "-"
             ENDIF
             FLUSH(6)
             EXIT outer
          ENDIF
       ENDDO outer

END SUBROUTINE
