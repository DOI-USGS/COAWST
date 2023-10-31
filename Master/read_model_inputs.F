#include "cppdefs.h"

#ifdef COAWST_COUPLING

      SUBROUTINE read_model_inputs
!
!=======================================================================
!                                                                      !
!  This routine reads in model input parameters of dt and              !
!  number of grids for each model.                                     !
!                                                                      !
!=======================================================================
!
      USE mct_coupler_params
# ifdef ROMS_MODEL
      USE mod_iounits
# endif
# ifdef SWAN_MODEL
      USE swan_iounits
# endif
# ifdef WW3_MODEL
      USE ww3_iounits
# endif
# if defined WRF_MODEL || defined WRFHYDRO_MODEL
      USE mod_coupler_iounits
# endif

      implicit none
!
      include 'mpif.h'
!
!  Imported variable declarations.
!
!
!  Local variable declarations.
!
      integer :: Npts, Nval, i, isval, iw, ia, inp, out, status
      integer :: MyRank, MyError, MyMaster, DT, num, den
      integer :: cdecode_line, cload_i, cload_r, indx, indx2, test
      integer :: Ivalue
# if defined WRF_COUPLING
      integer :: sstupdate
!     integer, allocatable :: parentid(:)
# endif
      real(m8), dimension(100) :: Rval
      real(m8) :: FAC

      character (len=1 ), parameter :: blank = ' '
      character (len=1 ) :: KEY
      character (len=40) :: KeyWord
      character (len=160) :: line
      character (len=160) :: aline
      character (len=160) :: saveline1, saveline2, saveline3
      character (len=160), dimension(100) :: Cval
!
      inp=1
      out=stdout
      MyMaster=0
      CALL mpi_comm_rank (MPI_COMM_WORLD, MyRank, MyError)
!
# if defined ROMS_COUPLING
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Read ROMS input file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      OPEN (inp, FILE=TRIM(Iname), FORM='formatted', STATUS='old',      &
     &      ERR=10)
      GO TO 30
 10   WRITE (out,20) Iname
      IF (MyRank.eq.MyMaster) WRITE(out,*) 'MyRank = ', MyRank,         &
     &                        TRIM(Iname)
!     exit_flag=4
      RETURN
 20   FORMAT (/,' READ MODEL INPUTS - Unable to open roms input file.', &
     &        /,a80)
 30   CONTINUE
!
      DO WHILE (.TRUE.)
        READ (inp,'(a)',ERR=15,END=40) line
        status=cdecode_line(line, KeyWord, Nval, Cval, Rval)
        IF (status.gt.0) THEN
          SELECT CASE (TRIM(KeyWord))
            CASE ('Ngrids')
              Npts=cload_i(Nval, Rval, 1, Ivalue)
              Nocn_grids=Ivalue
              IF (Nocn_grids.le.0) THEN
                IF (MyRank.eq.MyMaster) WRITE (out,290)'Ngrids',        &
     &            Nocn_grids, 'Ngrids must be greater than zero.'
!                exit_flag=5
                RETURN
              END IF
              allocate (dtocn(Nocn_grids))
            CASE ('DT')
              Npts=cload_r(Nval, Rval, Nocn_grids, dtocn)
          END SELECT
!         IF (exit_flag.ne.NoError) RETURN
        END IF
      END DO
 15   IF (MyRank.eq.MyMaster) WRITE (out,60) line
!     exit_flag=4
      RETURN
 40   CLOSE (inp)
 290  FORMAT (/,'read model inputs - Invalid dimension parameter,',a,i4,&
     &        /,15x,a)
# endif
# if defined SWAN_COUPLING
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Read SWAN input file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      iw=1
      OPEN (inp, FILE=TRIM(Wname(iw)), FORM='formatted', STATUS='old',  &
     &      ERR=110)
      GO TO 130
 110  WRITE (out,120) Wname(iw)
      IF (MyRank.eq.MyMaster) WRITE(out,*) 'MyRank = ', MyRank,         &
     &                        TRIM(Wname(iw))
!     exit_flag=4
      RETURN
 120  FORMAT (/,' READ MODEL INPUTS - Unable to open swan input file.', &
     &        /,a80)
 130  CONTINUE
!
      DO WHILE (.TRUE.)
        READ (inp,'(a)',ERR=116,END=140) line
        aline=ADJUSTL(line)
!
!  dont read comment lines
!
        isval=INDEX(aline,'&')
        IF(isval.eq.0) THEN
!
!  here we look for number of swan grids
!
          IF(aline(1:7).eq.'NSGRIDS') THEN
            read(aline(8:12),'(i5)') Nwav_grids
            allocate (dtwav(Nwav_grids))
          END IF
!
!  here we look for model time step
!
          IF(aline(1:4).eq.'COMP') THEN
            DO i=1,3
              indx=INDEX(aline,blank)
              aline=aline(indx+1:LEN(aline))
            END DO
            DO i=1,1
              indx=INDEX(aline,blank)
              read(aline(1:indx),'(i10)') DT
              dtwav(iw)=REAL(DT,m8)
              aline=aline(indx+1:LEN(aline))
            END DO
            READ(aline,'(a1)') KEY
            IF (KEY.eq.'D') THEN
              FAC = 24.0_m8*3600.0_m8
            ELSE IF (KEY.eq.'H') THEN
              FAC = 3600.0_m8
            ELSE IF (KEY.eq.'M') THEN
              FAC = 60.0_m8
            ELSE
              FAC = 1.0_m8
            ENDIF
            dtwav(iw)=dtwav(iw)*FAC
          END IF
#  ifdef SPECTRUM_STOKES
!
!  here we look for MSCs for SPECTUM_STOKES
!  such as CIRCLE 36 0.04 1.0 20
!
          indx2=INDEX(aline,'CIR')
          IF(indx2.ge.1) THEN
            aline=aline(indx2+3:LEN(aline))
            aline=ADJUSTL(aline)
!
            indx2=INDEX(aline,'C')
            IF(indx2.ge.1) THEN
              aline=aline(indx2+1:LEN(aline))
              aline=ADJUSTL(aline)
            END IF
!
            indx2=INDEX(aline,'L')
            IF(indx2.ge.1) THEN
              aline=aline(indx2+1:LEN(aline))
              aline=ADJUSTL(aline)
            END IF
!
            indx2=INDEX(aline,'E')
            IF(indx2.ge.1) THEN
              aline=aline(indx2+1:LEN(aline))
              aline=ADJUSTL(aline)
            END IF
!
            DO i=1,3
              indx=INDEX(aline,blank)
              aline=aline(indx+1:LEN(aline))
              aline=ADJUSTL(aline)
            END DO
            aline=ADJUSTL(aline)
            read(aline(1:5),'(i5)') MSCs
            MSCs=MSCs+1
          END IF
#  endif
        END IF
      END DO
 116  IF (MyRank.eq.MyMaster) WRITE (out,60) line
!     exit_flag=4
      RETURN
 140  CLOSE (inp)
!
      DO iw=2,Nwav_grids
        OPEN (inp, FILE=TRIM(Wname(iw)), FORM='formatted', STATUS='old',&
     &        ERR=110)
        GO TO 135
        IF (MyRank.eq.MyMaster) WRITE(out,*) 'MyRank = ', MyRank,       &
     &                          TRIM(Wname(iw))
!       exit_flag=4
        RETURN
 135    CONTINUE

        DO WHILE (.TRUE.)
          READ (inp,'(a)',ERR=115,END=160) line
          aline=ADJUSTL(line)
          IF(aline(1:4).eq.'COMP') THEN
            DO i=1,3
              indx=INDEX(aline,blank)
              aline=aline(indx+1:LEN(aline))
            END DO
            DO i=1,1
              indx=INDEX(aline,blank)
              read(aline(1:indx),'(i10)') DT
              dtwav(iw)=REAL(DT,m8)
              aline=aline(indx+1:LEN(aline))
            END DO
            READ(aline,'(a1)') KEY
            IF (KEY.eq.'D') THEN
              FAC = 24.0_m8*3600.0_m8
            ELSE IF (KEY.eq.'H') THEN
              FAC = 3600.0_m8
            ELSE IF (KEY.eq.'M') THEN
              FAC = 60.0_m8
            ELSE
              FAC = 1.0_m8
            ENDIF
            dtwav(iw)=dtwav(iw)*FAC
          END IF
        END DO
 115    IF (MyRank.eq.MyMaster) WRITE (out,60) line
!       exit_flag=4
        RETURN
 160    CLOSE (inp)
      END DO

# endif
# if defined WW3_COUPLING
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Read WW3 input file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      iw=1
      OPEN (inp, FILE=TRIM(Wname(iw)), FORM='formatted', STATUS='old',  &
     &      ERR=110)
      GO TO 130
 110  WRITE (out,120) Wname(iw)
      IF (MyRank.eq.MyMaster) WRITE(out,*) 'MyRank = ', MyRank,         &
     &                        TRIM(Wname(iw))
      RETURN
 120  FORMAT (/,' READ MODEL INPUTS - Unable to open ww3 input file.',  &
     &        /,a80)
 130  CONTINUE
!
!  set number of WW3 grids to = 1 for now.
!
      Nwav_grids=1
      allocate (dtwav(Nwav_grids))
!
!  now read the ww3_grid.inp file to get the dt
!
      REWIND (inp)
      test=1
      DO WHILE (test.eq.1)
! Name of grid :
        READ (inp,'(a)',ERR=116,END=140) line
        IF (line(1:1).ne.'$') test=0
      END DO
!
      test=1
      DO WHILE (test.eq.1)
! Define spectrum
        READ (inp,'(a)',ERR=116,END=140) line
        IF (line(1:1).ne.'$') test=0
      END DO
!
      test=1
      DO WHILE (test.eq.1)
! Define model run flags
        READ (inp,'(a)',ERR=116,END=140) line
        IF (line(1:1).ne.'$') test=0
      END DO
!
      test=1
      DO WHILE (test.eq.1)
! Define model time steps
        READ (inp,'(a)',ERR=116,END=140) line
        IF (line(1:1).ne.'$') test=0
      END DO
      aline=ADJUSTL(line)
      indx=INDEX(aline,blank)
      aline=aline(1:indx-1)
      read(aline,*) dtwav(iw)
! Set minimum 1 sec dt
      dtwav(iw)=MAX(dtwav(iw),1.0_m8)
!
!     write(*,*) 'ww3 dt is ', dtwav(iw)
!
      test=1
      DO WHILE (test.eq.1)
! Read rest of file
        READ (inp,'(a)',ERR=116,END=140) line
      END DO
 116  IF (MyRank.eq.MyMaster) WRITE (out,60) line
!     exit_flag=4
      RETURN
 140  CLOSE (inp)
!
# endif
# ifdef WRF_COUPLING
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Read WRF input file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      OPEN (inp, FILE=TRIM(Aname), FORM='formatted', STATUS='old',      &
     &      ERR=210)
      GO TO 230
 210  WRITE (out,220) Aname
      IF (MyRank.eq.MyMaster) WRITE(out,*) 'MyRank = ', MyRank,         &
     &                        TRIM(Aname)
!     exit_flag=4
      RETURN
 220  FORMAT (/,' READ MODEL INPUTS - Unable to open wrf input file.',  &
     &        /,a80)
 230  CONTINUE
!
      sstupdate=0
      DO WHILE (.TRUE.)
        READ (inp,'(a)',ERR=215,END=240) line
        aline=ADJUSTL(line)
        IF(aline(1:10).eq.'time_step ') THEN
          saveline1=aline
        END IF
        IF(aline(1:19).eq.'time_step_fract_num') THEN
          saveline2=aline
        END IF
        IF(aline(1:19).eq.'time_step_fract_den') THEN
          saveline3=aline
        END IF
        IF(aline(1:7).eq.'max_dom') THEN
          indx=INDEX(aline,'=')
          aline=aline(indx+1:LEN(aline))
          indx=INDEX(aline,',')
          read(aline(1:indx-1),'(i5)') Natm_grids
          allocate (dtatm(Natm_grids))
          allocate (parentid(Natm_grids))
          allocate (wrf_e_we(Natm_grids))
!
! Process the time steps
!
!  get DT
          aline=saveline1
          indx=INDEX(aline,'=')
          aline=aline(indx+1:LEN(aline))
          indx=INDEX(aline,',')
          read(aline(1:indx-1),'(i5)') DT
!  get num
          aline=saveline2
          indx=INDEX(aline,'=')
          aline=aline(indx+1:LEN(aline))
          indx=INDEX(aline,',')
          read(aline(1:indx-1),'(i5)') num
!  get den
          aline=saveline3
          indx=INDEX(aline,'=')
          aline=aline(indx+1:LEN(aline))
          indx=INDEX(aline,',')
          read(aline(1:indx-1),'(i5)') den
!  compute dt = DT + num/den
          IF (den.eq.0) THEN
            dtatm(1)=REAL(DT,m8)
          ELSE
            dtatm(1)=REAL(DT,m8)+REAL(num,m8)/REAL(den,m8)
          END IF
        END IF
        IF(aline(1:9).eq.'parent_id') THEN
          indx=INDEX(aline,'=')
          DO ia=1,Natm_grids
            aline=aline(indx+1:LEN(aline))
            indx=INDEX(aline,',')
            saveline1=TRIM(ADJUSTL(aline(1:indx-1)))
            read(saveline1,'(i5)') parentid(ia)
            IF (parentid(ia).EQ.0) parentid(ia) = 1
          END DO
        END IF
        IF(aline(1:22).eq.'parent_time_step_ratio') THEN
          indx=INDEX(aline,'=')
          DO ia=1,Natm_grids
            aline=aline(indx+1:LEN(aline))
            indx=INDEX(aline,',')
            saveline1=TRIM(ADJUSTL(aline(1:indx-1)))
            read(saveline1,'(i5)') den
!           dtatm(ia)=dtatm(1)/REAL(den,m8)
            dtatm(ia)=dtatm(parentid(ia))/REAL(den,m8)
          END DO
        END IF
        IF(aline(1:10).eq.'sst_update') THEN
          indx=INDEX(aline,'=')
          aline=ADJUSTL(aline(indx+1:LEN(aline)))
          indx=MAX(INDEX(aline,','),LEN(aline))
          read(aline(1:indx-1),'(i1)') sstupdate
        END IF
        IF(aline(1:4).eq.'e_we') THEN
          indx=INDEX(aline,'=')
          DO ia=1,Natm_grids
            aline=aline(indx+1:LEN(aline))
            indx=INDEX(aline,',')
            saveline1=TRIM(ADJUSTL(aline(1:indx-1)))
            read(saveline1,'(i5)') wrf_e_we(ia)
          END DO
        END IF
      END DO
!
 215  IF (MyRank.eq.MyMaster) WRITE (out,60) line
!     exit_flag=4
      RETURN
 240  CLOSE (inp)
#  if !defined SST_CONSTANT
      IF (sstupdate.eq.0) THEN
        WRITE (stdout,65) sstupdate
 65     FORMAT (/,' Recommend that sst_update be set to = 1 in the '    &
     &            'namelist.input for model coupling, not = ',i5)
!       STOP
      END IF
#  endif
# endif
# ifdef HYDRO_COUPLING
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Read WRF_HYDRO input file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      OPEN (inp, FILE=TRIM(Hname), FORM='formatted', STATUS='old',      &
     &      ERR=310)
      GO TO 330
 310  WRITE (out,320) Hname
      IF (MyRank.eq.MyMaster) WRITE(out,*) 'MyRank = ', MyRank,         &
     &                        TRIM(Hname)
!     exit_flag=4
      RETURN
 320  FORMAT (/,' READ MODEL INPUTS - Unable to open hyd input file.',  &
     &        /,a80)
 330  CONTINUE
!
!     For now hard code number of hyd grids at 1.
!
      Nhyd_grids=1
      allocate (dthyd(Nhyd_grids))
!
      DO WHILE (.TRUE.)
        READ (inp,'(a)',ERR=315,END=340) line
        aline=ADJUSTL(line)
!
! Process the time steps
!
        IF(aline(1:4).eq.'DTRT') THEN
          indx=INDEX(aline,'=')
          aline=aline(indx+1:LEN(aline))
          aline=ADJUSTL(aline)
          read(aline(1:5),'(i5)') DT
          dthyd(1)=REAL(DT,m8)
        END IF
      END DO
!
 315  IF (MyRank.eq.MyMaster) WRITE (out,60) line
!     exit_flag=4
      RETURN
 340  CLOSE (inp)
# endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  60  FORMAT (/,'read model inputs - Error while processing line: ',/,a)

      RETURN
      END SUBROUTINE read_model_inputs
      
#else
      SUBROUTINE read_model_inputs
      RETURN
      END SUBROUTINE read_model_inputs
#endif

