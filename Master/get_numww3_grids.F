#include "cppdefs.h"

#if defined WW3_MODEL && defined NESTING && defined COAWST_COUPLING
      SUBROUTINE get_numww3_grids(tempname)
!
!=======================================================================
!                                                                      !
!  This routine reads the ww3 model input and gets the number          !
!  of grids and dt's for a nested ww3 application.                     !
!                                                                      !
!=======================================================================
!
      USE mod_coupler_kinds
      USE mct_coupler_params
      USE ww3_iounits

      implicit none
!
      include 'mpif.h'
!
!  Imported variable declarations.
!
      character (len=160), intent(in) :: tempname
!
!  Local variable declarations.
!
      integer :: i, inp, out, status, num, test, indx, iw
      integer :: MyRank, MyError, MyMaster

      character (len=12)  :: moddef_name
      character (len=160) :: line
      character (len=160) :: aline
      character (len=1 ), parameter :: blank = ' '

!     Grab some code from w3iogrmd.F90
      integer :: NX, NY,NSEA, MTH,MK, NBI, NFBPO, GTYPE, IERR
      integer :: ICLOSE, JCWI, TRFLAG
      CHARACTER(LEN=10)       :: VERTST
      CHARACTER(LEN=13)       :: TEMPXT
      CHARACTER(LEN=30)       :: TNAME0, TNAME1, TNAME2, TNAME3,        &
     &                           TNAME4, TNAME5, TNAME6,                &
     &                           TNAMEP, TNAMEG, TNAMEF, TNAMEI
      CHARACTER(LEN=30)       :: FNAME0, FNAME1, FNAME2, FNAME3,        &
     &                           FNAME4, FNAME5, FNAME6,                &
     &                           FNAMEP, FNAMEG, FNAMEF, FNAMEI,        &
     &                           GNAME
      CHARACTER(LEN=35)       :: IDTST
      INTEGER, PARAMETER      :: RLGTYPE = 1
      INTEGER, PARAMETER      :: CLGTYPE = 2
      INTEGER, PARAMETER      :: UNGTYPE = 3
      INTEGER, PARAMETER      :: SMCTYPE = 4
      LOGICAL                 :: FLAGLL, JCWL
      REAL                    :: JCWR
      REAL, ALLOCATABLE       :: XGRD4(:,:), YGRD4(:,:)
      INTEGER, ALLOCATABLE       :: NBO(:), NBO2(:)

      REAL, ALLOCATABLE       :: TRNX(:,:), TRNY(:,:)
      REAL, ALLOCATABLE       ::    ZB(:)
      INTEGER, ALLOCATABLE    :: MAPTMP(:,:), MAPFS(:,:), MAPSF(:,:)

      REAL             :: SX, SY, X0, Y0, DTCFL, DTCFLI, DTMAX,      &
         DTMIN, DMIN, CTMAX, FICE0, FICEN, FICEL,   &
         PFMOVE, STEXU, STEYU, STEDU, IICEHMIN,     &
         IICEHINIT, ICESCALES(4), IICEHFAC, IICEHDISP, &
         IICEDDISP, IICEFDISP, BTBETA, AAIRCMIN, AAIRGB
      LOGICAL          :: GINIT, FLDRY, FLCX, FLCY, FLCTH, FLCK,    &
     &                    FLSOU, IICEDISP, IICESMOOTH
      REAL   , ALLOCATABLE :: CLATS(:)
      REAL   , ALLOCATABLE :: CLATIS(:)
      REAL   , ALLOCATABLE :: CTHG0S(:)
      LOGICAL          :: CMPRTRCK
      CHARACTER    :: CALTYPE*8
      LOGICAL               :: FLBPI, FLBPO
      
!$UNSTR
      LOGICAL :: FSN,FSPSI,FSFCT,FSNIMP,FSTOTALIMP,FSTOTALEXP
      LOGICAL :: FSREFRACTION, FSFREQSHIFT, FSSOURCE, FSBCCFL
      LOGICAL :: DO_CHANGE_WLV
      REAL(8) :: SOLVERTHR_STP
      REAL(8) :: CRIT_DEP_STP
      LOGICAL :: B_JGS_TERMINATE_MAXITER
      LOGICAL :: B_JGS_TERMINATE_DIFFERENCE
      LOGICAL :: B_JGS_TERMINATE_NORM
      LOGICAL :: B_JGS_LIMITER
      LOGICAL :: B_JGS_USE_JACOBI
      LOGICAL :: B_JGS_BLOCK_GAUSS_SEIDEL
      INTEGER :: B_JGS_MAXITER
      INTEGER :: B_JGS_LIMITER_FUNC
      REAL*8  :: B_JGS_PMIN
      REAL*8  :: B_JGS_DIFF_THR
      REAL*8  :: B_JGS_NORM_THR
      INTEGER :: B_JGS_NLEVEL
      LOGICAL :: B_JGS_SOURCE_NONLINEAR
      INTEGER  :: NTRI, COUNTOT, COUNTRI, NNZ
!
      REAL :: MAXX, MAXY, DXYMAX
      DOUBLE PRECISION, POINTER  :: XGRD(:,:), YGRD(:,:) 
      INTEGER, POINTER      :: TRIGP(:,:) 
      REAL(8), POINTER      :: TRIA(:)
      REAL(8), POINTER      :: LEN(:,:), SI(:), IEN(:,:)
      REAL, POINTER         :: ANGLE(:,:), ANGLE0(:,:)
      INTEGER, POINTER      :: CCON(:), COUNTCON(:), IE_CELL(:),        &
     &                        POS_CELL(:),  IAA(:), JAA(:),             &
     &                        POSI(:,:), INDEX_CELL(:)
      INTEGER*2, POINTER    :: IOBP(:)
      INTEGER*1, POINTER    :: IOBPD(:,:), IOBDP(:), IOBPA(:)

!
      inp=2
      out=6
      MyMaster=0
      CALL mpi_comm_rank (MPI_COMM_WORLD, MyRank, MyError)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   First read the ww3_multi.inp file to get the number of grids.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      iw=1
      OPEN (inp, FILE=TRIM(Wname(iw)), FORM='formatted', STATUS='old',  &
     &      ERR=210)
      GO TO 230
 210  WRITE (out,220) Wname(iw)
      IF (MyRank.eq.MyMaster) WRITE(out,*) 'MyRank = ', MyRank,         &
     &                        TRIM(Wname(iw))
      RETURN
 220  FORMAT (/,' READ MODEL INPUTS - Unable to open ww3 input file.',  &
     &        /,a80)
 230  CONTINUE
!
!  now read the ww3_multi.inp file to get number of grids
!
      REWIND (inp)
      test=1
      DO WHILE (test.eq.1)
! Read the file until the first non-comment line :
        READ (inp,'(a)',ERR=116,END=240) line
        IF (line(1:1).ne.'$') test=0
      END DO
!
      test=1
      aline=ADJUSTL(line)
      indx=INDEX(aline,blank)
      aline=aline(1:indx-1)
      read(aline,*) Nwav_grids
      allocate (dtwav(Nwav_grids))
!
 240  CLOSE (inp)
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     Now read WW3 mod_def files to get the dt's
!     Grab some code from w3iogrmd.F90
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      DO num=1,Nwav_grids
        write(moddef_name,'(A11,I1)') 'mod_def.grd',num
        OPEN (inp, FILE=TRIM(moddef_name), FORM='unformatted',          &
     &        STATUS='old', convert='big_endian',ERR=110)
      GO TO 130
 110  WRITE (out,120) tempname
      IF (MyRank.eq.MyMaster) WRITE(out,*) 'MyRank = ', MyRank,         &
     &                        TRIM(tempname)
!     exit_flag=4
      RETURN
 120  FORMAT (/,' get num ww3 grids - Unable to open ww3 input file.',  &
     &        /,a80)
 130  CONTINUE
!
!      DO WHILE (.TRUE.)
    REWIND ( inp )

      READ (inp,END=801,ERR=802,IOSTAT=IERR)                    &
           IDTST, VERTST, NX, NY, NSEA, MTH, MK,                 &
           NBI, NFBPO, GNAME, FNAME0, FNAME1, FNAME2, FNAME3,    &
           FNAME4, FNAME5, FNAME6, FNAMEP, FNAMEG,               &
           FNAMEF, FNAMEI
      !

      ALLOCATE (NBO(0:NFBPO), NBO2(0:NFBPO)) 
      READ (inp,END=801,ERR=802,IOSTAT=IERR)                   &
           (NBO(I),I=0,NFBPO), (NBO2(I),I=0,NFBPO)
      DEALLOCATE (NBO, NBO2) 

      READ (inp,END=801,ERR=802,IOSTAT=IERR)                      &
           GTYPE, FLAGLL, ICLOSE

      SELECT CASE ( GTYPE )
      CASE ( RLGTYPE )
       READ (inp,END=801,ERR=802,IOSTAT=IERR)                 &
             SX, SY, X0, Y0

      CASE ( CLGTYPE )
        ALLOCATE(XGRD4(NY,NX),YGRD4(NY,NX)); XGRD4 = 0.; YGRD4 = 0.
        READ (inp,END=801,ERR=802,IOSTAT=IERR)               &
             XGRD4, YGRD4
        DEALLOCATE(XGRD4, YGRD4)

      CASE (UNGTYPE)

        ALLOCATE ( XGRD(1,NX), YGRD(1,NX) )
        ALLOCATE ( CCON(NX), COUNTCON(NX) )
        ALLOCATE ( IAA(NX+1) )
        ALLOCATE ( INDEX_CELL(NX+1) )
        ALLOCATE ( IOBP(NX) )
        ALLOCATE ( IOBPD(MTH,NX), IOBDP(NX), IOBPA(NX) )

        READ (inp,END=801,ERR=802,IOSTAT=IERR)                          &
     &        FSN, FSPSI,FSFCT,FSNIMP,FSTOTALIMP,FSTOTALEXP,            &
     &        FSBCCFL, FSREFRACTION, FSFREQSHIFT, FSSOURCE,             &
     &        DO_CHANGE_WLV, SOLVERTHR_STP, CRIT_DEP_STP,      &
     &        NTRI, COUNTOT, COUNTRI, NNZ,                     &
     &        B_JGS_TERMINATE_MAXITER,                         &
     &        B_JGS_TERMINATE_DIFFERENCE,                      &
     &        B_JGS_TERMINATE_NORM,                            &
     &        B_JGS_LIMITER,                                   &
     &        B_JGS_BLOCK_GAUSS_SEIDEL,                        &
     &        B_JGS_USE_JACOBI,                                &
     &        B_JGS_MAXITER,                                   &
     &        B_JGS_PMIN,                                      &
     &        B_JGS_DIFF_THR,                                  &
     &        B_JGS_NORM_THR,                                  &
     &        B_JGS_NLEVEL,                                    &
     &        B_JGS_SOURCE_NONLINEAR

        ALLOCATE ( TRIGP(3,NTRI) )
        ALLOCATE ( TRIA(NTRI), LEN(NTRI,3),SI(NX), IEN(NTRI,6) )
        ALLOCATE ( ANGLE(NTRI,3), ANGLE0(NTRI,3) )
        ALLOCATE ( IE_CELL(COUNTOT) )
        ALLOCATE ( POS_CELL(COUNTOT) )
        ALLOCATE ( POSI(3,COUNTOT) )
        ALLOCATE ( JAA(NNZ) )

        READ (inp,END=801,ERR=802,IOSTAT=IERR)               &
             X0, Y0, SX, SY, DXYMAX, XGRD, YGRD, TRIGP, TRIA, &
             LEN, IEN, ANGLE0, ANGLE, SI, MAXX, MAXY,         &
             DXYMAX, INDEX_CELL, CCON, COUNTCON, IE_CELL,     &
             POS_CELL, IOBP, IOBPA, IOBDP, IOBPD, IAA, JAA, POSI


        DEALLOCATE ( XGRD, YGRD )
        DEALLOCATE ( TRIGP )
        DEALLOCATE ( TRIA, LEN, SI, IEN )
        DEALLOCATE ( ANGLE, ANGLE0 )
        DEALLOCATE ( CCON, COUNTCON, IE_CELL )
        DEALLOCATE ( POS_CELL, IAA, JAA )
        DEALLOCATE ( POSI, INDEX_CELL )
        DEALLOCATE ( IOBP )
        DEALLOCATE ( IOBPD, IOBDP, IOBPA )



      END SELECT !GTYPE
 

      ALLOCATE ( ZB(1:NSEA) )
      ALLOCATE ( MAPTMP(NY,NX), MAPFS(NY,NX), MAPSF(NSEA,3) )
      READ (inp,END=801,ERR=802,IOSTAT=IERR)   &
     &            ZB, MAPTMP, MAPFS, MAPSF, TRFLAG
      DEALLOCATE ( ZB, MAPTMP, MAPFS, MAPSF )

      IF ( TRFLAG .NE. 0 ) THEN
         ALLOCATE ( TRNX(NY,NX), TRNY(NY,NX) )
         READ (inp,END=801,ERR=802,IOSTAT=IERR) TRNX, TRNY

        if (myrank.eq.0) then
!        write(*,*) 'trnx is ', trnx
        endif
        if (myrank.eq.0) then
!        write(*,*) 'trny is ', trny
        endif
         DEALLOCATE ( TRNX, TRNY )
      END IF

      ALLOCATE ( CLATS(0:NSEA), CLATIS(0:NSEA), CTHG0S(0:NSEA) )

      READ (inp,END=801,ERR=802,IOSTAT=IERR)                       &
     &           DTCFL, DTCFLI, DTMAX, DTMIN, DMIN, CTMAX,         &
     &      FICE0, FICEN, FICEL, PFMOVE, FLDRY, FLCX, FLCY,        &
     &      FLCTH, FLCK, FLSOU, FLBPI, FLBPO, CLATS, CLATIS,       &
     &      CTHG0S, STEXU, STEYU, STEDU, IICEHMIN, IICEHINIT,      &
     &      IICEDISP, ICESCALES(1:4), CALTYPE, CMPRTRCK, IICEHFAC, &
     &      IICEDDISP, IICEHDISP, IICEFDISP, BTBETA,               &
     &      AAIRCMIN, AAIRGB

      DEALLOCATE ( CLATS, CLATIS, CTHG0S )

!  assign dtwav
        dtwav(num)=DTCFL

      END DO
 140  CLOSE (inp)
      RETURN
 801  CONTINUE
 802  WRITE (out,803) tempname
 803  FORMAT (/,' get num ww3 grids - error reading  ww3 input file.',  &
     &        /,a80)
      IF (MyRank.eq.MyMaster) WRITE(out,*) 'MyRank = ', MyRank,         &
     &                        TRIM(tempname)
!     exit_flag=4
      RETURN
!
 116  IF (MyRank.eq.MyMaster) WRITE (out,60) line
!     exit_flag=4
      RETURN

  60  FORMAT (/,'read model inputs - Error while processing line: ',/,a)

      RETURN
      END SUBROUTINE get_numww3_grids
#else
      SUBROUTINE get_numww3_grids(tempname)
      character (len=160), intent(in) :: tempname
      RETURN
      END SUBROUTINE get_numww3_grids
#endif

