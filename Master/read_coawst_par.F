#include "cppdefs.h"

#ifdef COAWST_COUPLING

      SUBROUTINE read_coawst_par(flag)
!
!=======================================================================
!                                                                      !
!  This routine reads in physical model input parameters.              !
!                                                                      !
!=======================================================================
!
      USE mct_coupler_params
      USE mod_coupler_iounits
      USE mct_coupler_utils_mod
# ifdef ROMS_MODEL
      USE mod_iounits
# endif
# ifdef SWAN_MODEL
      USE swan_iounits
      USE M_COUPLING
# endif
      implicit none
!
      include 'mpif.h'
!
!  Imported variable declarations.
!
      integer, intent(in) :: flag
!
!  Local variable declarations.
!
      integer :: Npts, Nval, i, iw, io, ia, j, inp, out, status
      integer :: MyRank, MyError, MyMaster, Nchars
      integer :: cdecode_line, cload_i, cload_r, count

      real(m8), dimension(100) :: Rval

      character (len=1 ), parameter :: blank = ' '
      character (len=40) :: KeyWord
      character (len=80) :: Cname
      character (len=160) :: tempname
      character (len=160) :: line
      character (len=160), dimension(100) :: Cval
      character (len=80) :: varnam
!
!-----------------------------------------------------------------------
!  Determine coupling standard input file name.  In distributed-memory,
!  this name is assigned at the executtion command line and processed
!  with the Unix routine GETARG.  The ROMS/TOMS input parameter script
!  name is specified in this coupling script.
!-----------------------------------------------------------------------
!
      inp=1
      out=stdout
      MyMaster=0
!
      CALL mpi_comm_rank (MPI_COMM_WORLD, MyRank, MyError)
      IF (MyRank.eq.MyMaster) CALL mct_getarg (1,Cname)
      Nchars=LEN(Cname)
      CALL mpi_bcast (Cname, Nchars, MPI_BYTE, MyMaster,                &
     &                MPI_COMM_WORLD, MyError)
      IF (MyRank.eq.MyMaster) THEN
        WRITE(out,*) 'Coupled Input File name = ', TRIM(Cname)
      END IF
      OPEN (inp, FILE=TRIM(Cname), FORM='formatted', STATUS='old',      &
     &      ERR=10)
      GO TO 30
 10   WRITE (out,20)
      IF (MyRank.eq.MyMaster) WRITE(out,*) 'MyRank = ', MyRank,         &
     &                        TRIM(Cname)
!     exit_flag=4
      RETURN
 20   FORMAT (/,' READ_COAWST_PAR - Unable to open coupling script.',   &
     &        /,11x,'In distributed-memory applications, the input',    &
     &        /,11x,'script file is processed in parallel. The Unix',   &
     &        /,11x,'routine GETARG is used to get script file name.',  &
     &        /,11x,'For example, in MPI applications make sure that',  &
     &        /,11x,'command line is something like:',/,                &
     &        /,11x,'mpirun -np 4 masterM coupling.in',/,               &
     &        /,11x,'and not',/,                                        &
     &        /,11x,'mpirun -np 4 masterM < coupling.in',/)
 30   CONTINUE
!
!-----------------------------------------------------------------------
!  Read in multiple models coupling parameters. Then, load input
!  data into module. Take into account nested grid configurations.
!-----------------------------------------------------------------------
!
      IF (flag.eq.1) THEN
        scrip_opt=0
        DO WHILE (.TRUE.)
          READ (inp,'(a)',ERR=40,END=50) line
          status=cdecode_line(line, KeyWord, Nval, Cval, Rval)
          IF (status.gt.0) THEN
            IF (TRIM(KeyWord).eq.'VARNAME') THEN
              DO i=1,LEN(Cname)
                varnam(i:i)=blank
              END DO
              varnam=TRIM(ADJUSTL(Cval(Nval)))
            ELSE IF (TRIM(KeyWord).eq.'NnodesATM') THEN
              Npts=cload_i(Nval, Rval, 1, NnodesATM)
            ELSE IF (TRIM(KeyWord).eq.'NnodesWAV') THEN
              Npts=cload_i(Nval, Rval, 1, NnodesWAV)
            ELSE IF (TRIM(KeyWord).eq.'NnodesOCN') THEN
              Npts=cload_i(Nval, Rval, 1, NnodesOCN)
            ELSE IF (TRIM(KeyWord).eq.'TI_ATM2WAV') THEN
              Npts=cload_r(Nval, Rval, 1, TI_ATM2WAV)
            ELSE IF (TRIM(KeyWord).eq.'TI_ATM2OCN') THEN
              Npts=cload_r(Nval, Rval, 1, TI_ATM2OCN)
            ELSE IF (TRIM(KeyWord).eq.'TI_WAV2ATM') THEN
              Npts=cload_r(Nval, Rval, 1, TI_WAV2ATM)
            ELSE IF (TRIM(KeyWord).eq.'TI_WAV2OCN') THEN
              Npts=cload_r(Nval, Rval, 1, TI_WAV2OCN)
            ELSE IF (TRIM(KeyWord).eq.'TI_OCN2WAV') THEN
              Npts=cload_r(Nval, Rval, 1, TI_OCN2WAV)
            ELSE IF (TRIM(KeyWord).eq.'TI_OCN2ATM') THEN
              Npts=cload_r(Nval, Rval, 1, TI_OCN2ATM)
# ifdef ROMS_MODEL
            ELSE IF (TRIM(KeyWord).eq.'OCN_name') THEN
              DO i=1,LEN(Iname)
                Iname(i:i)=blank
              END DO
              Iname=TRIM(ADJUSTL(Cval(Nval)))
# endif
# ifdef SWAN_MODEL
            ELSE IF (TRIM(KeyWord).eq.'WAV_name') THEN
              tempname=TRIM(ADJUSTL(Cval(Nval)))
              IF (Nval.eq.1) THEN
                call get_numswan_grids(tempname)
                call allocate_swan_iounits
              END IF
              Wname(Nval)=TRIM(ADJUSTL(Cval(Nval)))
# endif
# ifdef WRF_MODEL
            ELSE IF (TRIM(KeyWord).eq.'ATM_name') THEN
              DO i=1,LEN(Aname)
                Aname(i:i)=blank
              END DO
              Aname=TRIM(ADJUSTL(Cval(Nval)))
# endif
            END IF
            IF (TRIM(KeyWord).eq.'SCRIP_WEIGHT_OPTION') THEN
              Npts=cload_i(Nval, Rval, 1, scrip_opt)
            END IF
          END IF
        END DO
      END IF
      IF (flag.eq.2) THEN
        DO WHILE (.TRUE.)
          READ (inp,'(a)',ERR=40,END=50) line
          status=cdecode_line(line, KeyWord, Nval, Cval, Rval)
          IF (status.gt.0) THEN
            IF (scrip_opt.eq.1) THEN
              IF (TRIM(KeyWord).eq.'SCRIP_COAWST_NAME') THEN
                SCRIPname=TRIM(ADJUSTL(Cval(Nval)))
              END IF
            ELSE IF (scrip_opt.eq.2) THEN
# ifdef ROMS_COUPLING
#  ifdef MCT_INTERP_OC2AT
              IF (TRIM(KeyWord).eq.'O2ANAME') THEN
                io=CEILING(REAL(Nval)/REAL(Natm_grids))
                ia=Nval-(Natm_grids*(io-1))
                O2Aname(io,ia)=TRIM(ADJUSTL(Cval(Nval)))
              END IF
              IF (TRIM(KeyWord).eq.'A2ONAME') THEN
                ia=CEILING(REAL(Nval)/REAL(Nocn_grids))
                io=Nval-(Nocn_grids*(ia-1))
                A2Oname(ia,io)=TRIM(ADJUSTL(Cval(Nval)))
              END IF
#  endif
#  ifdef MCT_INTERP_OC2WV
              IF (TRIM(KeyWord).eq.'O2WNAME') THEN
                io=CEILING(REAL(Nval)/REAL(Nwav_grids))
                iw=Nval-(Nwav_grids*(io-1))
                O2Wname(io,iw)=TRIM(ADJUSTL(Cval(Nval)))
              END IF
              IF (TRIM(KeyWord).eq.'W2ONAME') THEN
                iw=CEILING(REAL(Nval)/REAL(Nocn_grids))
                io=Nval-(Nocn_grids*(iw-1))
                W2Oname(iw,io)=TRIM(ADJUSTL(Cval(Nval)))
              END IF
#  endif
# endif
# ifdef SWAN_COUPLING
#  ifdef MCT_INTERP_WV2AT
              IF (TRIM(KeyWord).eq.'A2WNAME') THEN
                ia=CEILING(REAL(Nval)/REAL(Nwav_grids))
                iw=Nval-(Nwav_grids*(ia-1))
                A2Wname(ia,iw)=TRIM(ADJUSTL(Cval(Nval)))
              END IF
              IF (TRIM(KeyWord).eq.'W2ANAME') THEN
                iw=CEILING(REAL(Nval)/REAL(Natm_grids))
                ia=Nval-(Natm_grids*(iw-1))
                W2Aname(iw,ia)=TRIM(ADJUSTL(Cval(Nval)))
              END IF
#  endif
# endif
            ELSE
              IF (MyRank.eq.MyMaster) WRITE (out,70)
            END IF
          END IF
        END DO
      END IF
  40  IF (MyRank.eq.MyMaster) WRITE (out,60) line
!     exit_flag=4
      RETURN
  50  CLOSE (inp)
  60  FORMAT (/,' read_coawst_par - Error while processing line: ',/,a)
  70  FORMAT (/,' read_coawst_par - Invalid SCRIP_WEIGHT_OPTION ',/,a)

      RETURN
      END SUBROUTINE read_coawst_par
      FUNCTION cdecode_line (line_text, KeyWord, Nval, Cval, Rval)
!
!=======================================================================
!                                                                      !
!  This function decodes lines of text from input script files.        !
!                                                                      !
!=======================================================================
!
      USE mod_coupler_kinds
!
      implicit none
!
! Imported variable declarations.
!
      character (len=*), intent(in) :: line_text

      character (len=40), intent(inout) :: KeyWord

      integer, intent(inout) :: Nval

      character (len=160), dimension(100), intent(inout) :: Cval

      real(m8), dimension(100), intent(inout) :: Rval
!
! Local variable declarations
!
      logical :: IsString, Kextract, decode, nested
      integer :: Iblank, Icmm, Kstr, Kend, Linp
      integer :: Lend, LenS, Lstr, Lval, Nmul, Schar
      integer :: copies, i, ic, ie, is, j, status

      integer, dimension(20) :: Imul

      integer :: cdecode_line

      character (len=1 ), parameter :: blank = ' '
      character (len=160) :: Vstring, line, string
!
!------------------------------------------------------------------------
!  Decode input line.
!------------------------------------------------------------------------
!
!  Initialize.
!
      DO i=1,LEN(line)
        line(i:i)=blank
        Vstring(i:i)=blank
        string(i:i)=blank
      END DO
!
!  Get length of "line". Remove comment after the KEYWORD, if any.
!  Then, remove leading and trailing blanks.
!
      Linp=LEN(line_text)
      IF ((Linp.gt.0).and.(line_text(1:1).ne.CHAR(33))) THEN
        Icmm=INDEX(line_text,CHAR(33),BACK=.FALSE.)
        IF (Icmm.gt.0) Linp=Icmm-1
        line=TRIM(ADJUSTL(line_text(1:Linp)))
        Linp=LEN_TRIM(line)
      ELSE
        line=TRIM(ADJUSTL(line_text))
        Linp=LEN_TRIM(line)
      END IF
!
!  If not a blank or comment line [char(33)=!], decode and extract input
!  values.  Find equal sign [char(61)].
!
      status=-1
      nested=.FALSE.
      IF ((Linp.gt.0).and.(line(1:1).ne.CHAR(33))) THEN
        status=1
        Kstr=1
        Kend=INDEX(line,CHAR(61),BACK=.FALSE.)-1
        Lstr=INDEX(line,CHAR(61),BACK=.TRUE.)+1
!
! Determine if KEYWORD is followed by double equal sign (==) indicating
! nested parameter.
!
        IF ((Lstr-Kend).eq.3) nested=.TRUE.
!
! Extract KEYWORD, trim leading and trailing blanks.
!
        Kextract=.FALSE.
        IF (Kend.gt.0) THEN
          Lend=Linp
          KeyWord=line(Kstr:Kend)
          Nval=0
          Kextract=.TRUE.
        ELSE
          Lstr=1
          Lend=Linp
          Kextract=.TRUE.
        END IF
!
! Extract parameter values string.  Remove continuation symbol
! [char(92)=\], if any.  Trim leading trailing blanks.
!
        IF (Kextract) THEN
          Icmm=INDEX(line,CHAR(92),BACK=.FALSE.)
          IF (Icmm.gt.0) Lend=Icmm-1
          Vstring=ADJUSTL(line(Lstr:Lend))
          Lval=LEN_TRIM(Vstring)
!
! The TITLE KEYWORD is a special one since it can include strings,
! numbers, spaces, and continuation symbol.
!
          IsString=.FALSE.
          IF (TRIM(KeyWord).eq.'TITLE') THEN
            Nval=Nval+1
            Cval(Nval)=Vstring(1:Lval)
            IsString=.TRUE.
          ELSE
!
! Check if there is a multiplication symbol [char(42)=*] in the variable
! string indicating repetition of input values.
!
            Nmul=0
            DO i=1,Lval
              IF (Vstring(i:i).eq.CHAR(42)) THEN
                Nmul=Nmul+1
                Imul(Nmul)=i
              END IF
            END DO
            ic=1
!
! Check for blank spaces [char(32)=' '] between entries and decode.
!
            is=1
            ie=Lval
            Iblank=0
            decode=.FALSE.
            DO i=1,Lval
              IF (Vstring(i:i).eq.CHAR(32)) THEN
                IF (Vstring(i+1:i+1).ne.CHAR(32)) decode=.TRUE.
                Iblank=i
              ELSE
                ie=i
              ENDIF
              IF (decode.or.(i.eq.Lval)) THEN
                Nval=Nval+1
!
! Processing numeric values.  Check starting character to determine
! if numeric or character values. It is possible to have both when
! processing repetitions via the multiplication symbol.
!
                Schar=ICHAR(Vstring(is:is))
                IF (((48.le.Schar).and.(Schar.le.57)).or.               &
     &              (Schar.eq.43).or.(Schar.eq.45)) THEN
                  IF ((Nmul.gt.0).and.                                  &
     &                (is.lt.Imul(ic)).and.(Imul(ic).lt.ie)) THEN
                    READ (Vstring(is:Imul(ic)-1),*) copies
                    Schar=ICHAR(Vstring(Imul(ic)+1:Imul(ic)+1))
                    IF ((43.le.Schar).and.(Schar.le.57)) THEN
                      READ (Vstring(Imul(ic)+1:ie),*) Rval(Nval)
                      DO j=1,copies-1
                        Rval(Nval+j)=Rval(Nval)
                      END DO
                    ELSE
                      string=Vstring(Imul(ic)+1:ie)
                      LenS=LEN_TRIM(string)
                      Cval(Nval)=string(1:LenS)
                      DO j=1,copies-1
                        Cval(Nval+j)=Cval(Nval)
                      END DO
                    END IF
                    Nval=Nval+copies-1
                    ic=ic+1
                  ELSE
                    string=Vstring(is:ie)
                    LenS=LEN_TRIM(string)
                    READ (string(1:LenS),*) Rval(Nval)
                  END IF
                ELSE
!
! Processing character values (logicals and strings).
!
                  IF ((Nmul.gt.0).and.                                  &
     &                (is.lt.Imul(ic)).and.(Imul(ic).lt.ie)) THEN
                    READ (Vstring(is:Imul(ic)-1),*) copies
                    Cval(Nval)=Vstring(Imul(ic)+1:ie)
                    DO j=1,copies-1
                      Cval(Nval+j)=Cval(Nval)
                    END DO
                    Nval=Nval+copies-1
                    ic=ic+1
                  ELSE
                    string=Vstring(is:ie)
                    Cval(Nval)=TRIM(ADJUSTL(string))
                  END IF
                  IsString=.TRUE.
                END IF
                is=Iblank+1
                ie=Lval
                decode=.FALSE.
              END IF
            END DO
          END IF
        END IF
        status=Nval
      END IF
      cdecode_line=status
      RETURN
      END FUNCTION cdecode_line

      FUNCTION cload_i (Ninp, Vinp, Nout, Vout)
!
!=======================================================================
!                                                                      !
!  This function loads input values into a requested model integer     !
!  variable.                                                           !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     Ninp       Size of input variable.                               !
!     Vinp       Input values                                          !
!     Nout       Number of output values.                              !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     Vout       Output integer variable.                              !
!     cload_i    Number of output values processed.                   !
!                                                                      !
!=======================================================================
!
      USE mod_coupler_kinds
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: Ninp, Nout
      real(m8), intent(in) :: Vinp(Ninp)
      integer, intent(out) :: Vout(Nout)
!
!  Local variable declarations.
!
      integer :: i, ic
      integer :: cload_i
!
!-----------------------------------------------------------------------
!  Load integer variable with input values.
!-----------------------------------------------------------------------
!
!  If not all values are provided for variable, assume the last value
!  for the rest of the array.
!
      ic=0
      IF (Ninp.le.Nout) THEN
        DO i=1,Ninp
          ic=ic+1
          Vout(i)=INT(Vinp(i))
        END DO
        DO i=Ninp+1,Nout
          ic=ic+1
          Vout(i)=INT(Vinp(Ninp))
        END DO
      ELSE
        DO i=1,Nout
          ic=ic+1
          Vout(i)=INT(Vinp(i))
        END DO
      END IF
      cload_i=ic

      RETURN
      END FUNCTION cload_i

      FUNCTION cload_r (Ninp, Vinp, Nout, Vout)
!
!=======================================================================
!                                                                      !
!  This function loads input values into a requested model real        !
!  variable.                                                           !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     Ninp       Size of input variable.                               !
!     Vinp       Input values                                          !
!     Nout       Number of output values.                              !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     Vout       Output real variable.                                 !
!     cload_r     Number of output values processed.                   !
!                                                                      !
!=======================================================================
!
      USE mod_coupler_kinds
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: Ninp, Nout
      real(m8), intent(in) :: Vinp(Ninp)
      real(m8), intent(out) :: Vout(Nout)
!
!  Local variable declarations.
!
      integer :: i, ic
      integer :: cload_r
!
!-----------------------------------------------------------------------
!  Load integer variable with input values.
!-----------------------------------------------------------------------
!
!  If not all values are provided for variable, assume the last value
!  for the rest of the array.
!
      ic=0
      IF (Ninp.le.Nout) THEN
        DO i=1,Ninp
          ic=ic+1
          Vout(i)=Vinp(i)
        END DO
        DO i=Ninp+1,Nout
          ic=ic+1
          Vout(i)=Vinp(Ninp)
        END DO
      ELSE
        DO i=1,Nout
          ic=ic+1
          Vout(i)=Vinp(i)
        END DO
      END IF
      cload_r=ic

      RETURN
      END FUNCTION cload_r

#else
      SUBROUTINE read_coawst_par
      RETURN
      END SUBROUTINE read_coawst_par
#endif

