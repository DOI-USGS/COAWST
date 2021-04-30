#include "cppdefs.h"
!
!************************************************************************
      SUBROUTINE swan_reader (ng, first, tile)
!************************************************************************
!
#if defined INWAVE_MODEL & defined INWAVE_SWAN_COUPLING
!
!svn $Id: swan_reader.F 1336 2008-01-24 02:45:56Z jcwarner $
! LAST CHANGE: mai 12/28/2010
!
!======================================================================!
!                                                                      !
!  This routine reads the output spectra from swan.................... !
!                                                                      !
!======================================================================!
!
      USE inwave_iounits
      USE mod_inwave_swan
      USE mod_inwave_vars
      USE mod_parallel
      USE mod_scalars
      USE Tr3dbc_mod
      USE dateclock_mod, ONLY : caldate
# ifdef DISTRIBUTE
      USE mp_exchange_mod
# endif
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, first, tile
!
!  Local variable declarations.
!
      character(6)          :: rtext
      character(72)         :: rline
      real(r8)              :: factor, m0, frep, cff, DTTIME
      integer               :: nfreq, ndir, switch
      integer               :: i, j, d, inp, min1, min2
      integer               :: nfreqsw, ndirsw
      integer               :: minpos, readerr
      integer               :: INTTIM(6)
      integer               :: YYYY, NN, DD, HH, MN, SS
      integer               :: LBi, UBi, LBj, UBj
!     integer, allocatable  :: posang(:)

      real(r8)            :: OceanTime
      real(r8), pointer   :: temp(:), findline(:)
      real(r8), pointer   :: tempA(:,:)
      real(r8), pointer   :: sdsw1(:,:)
      real(r8), pointer   :: theta1(:)

      real(r8), parameter :: rad = pi/180.0_r8
      real(r8), parameter :: orad = 180.0_r8/pi

# include "set_bounds.h"
!
!  Determine array lower and upper bounds in the I- and J-directions.
!
      LBi=BOUNDS(ng)%LBi(tile)
      UBi=BOUNDS(ng)%UBi(tile)
      LBj=BOUNDS(ng)%LBj(tile)
      UBj=BOUNDS(ng)%UBj(tile)
!
!-----------------------------------------------------------------------
! Open the spectral file
!-----------------------------------------------------------------------
!
      readerr=0
      inp=44
      open(inp,file=IWSWNname(ng),form='formatted',status='old')
!
!  If first, then we have not been here before.  read file and allocate
!  some arrays. for other entry times, dont need to do this.
!
      IF (first.eq.1) THEN
        switch=0
!
!-----------------------------------------------------------------------
! Read file until RFREQ or AFREQ is found
!-----------------------------------------------------------------------
!
        DO while (switch==0)
          read(inp,'(a)')rtext
          IF (rtext == 'RFREQ ') THEN
            switch = 1
          ELSEIF (rtext == 'AFREQ ') THEN
            switch = 2
          END IF
        ENDDO
!
!-----------------------------------------------------------------------
! Read nfreq and f
!-----------------------------------------------------------------------
!
        read(inp,*)nfreqsw
        allocate(WAVES(ng)%fsw(nfreqsw))

        DO i=1,nfreqsw
          read(inp,*) WAVES(ng)%fsw(i)
        ENDDO
        WAVES(ng)%nfreqsw=nfreqsw
!
!-----------------------------------------------------------------------
! Convert to absolute frequencies
!-----------------------------------------------------------------------
!
        IF (switch == 1) THEN
          WAVES(ng)%fsw = WAVES(ng)%fsw
        ELSE
          WAVES(ng)%fsw = WAVES(ng)%fsw
        ENDIF
!
!-----------------------------------------------------------------------
! Read CDIR or NDIR
!-----------------------------------------------------------------------
!
        read(inp,'(a)')rtext
        IF (rtext == 'NDIR  ') THEN
          switch = 1
        ELSEIF (rtext == 'CDIR  ') THEN
          switch = 2
        ELSE
          stop
        ENDIF
!
!-----------------------------------------------------------------------
! Read ndir, theta
!-----------------------------------------------------------------------
!
        read(inp,*) ndirsw
        allocate(theta1(ndirsw))
        allocate(WAVES(ng)%theta(ndirsw))
        allocate(WAVES(ng)%posang(ndirsw))

        DO i=1,ndirsw
          read(inp,*)theta1(i)
          IF(theta1(i).lt.0.0_r8)THEN
            theta1(i)=(360.0_r8+theta1(i))
          ENDIF
          theta1(i)=theta1(i)*rad
        ENDDO

        WAVES(ng)%ndirsw=ndirsw
        WAVES(ng)%dang=abs(theta1(2)-theta1(1))
        minpos=minloc(theta1, DIM = 1)
!       allocate(posang(ndirsw))

!       If swan was run on a full circle or just a sector.
        IF (NINT(2.0_r8*pi/WAVES(ng)%dang).eq.ndirsw) THEN
          WAVES(ng)%Swancircle=1
        ELSE
          WAVES(ng)%Swancircle=0
        END IF
!
!        IF((theta1(minpos+1)-theta1(minpos)).lt.                        &
!     &     (theta1(minpos-1)-theta1(minpos)))THEN
        min1=minpos
        min2=minpos
        IF (minpos.eq.1)      min1=ndirsw+1
        IF (minpos.eq.ndirsw) min2=0
        IF((theta1(min2+1)-theta1(minpos)).lt.                          &
     &     (theta1(min1-1)-theta1(minpos)))THEN
          DO i=1,ndirsw
            WAVES(ng)%posang(i)=minpos+(i-1)
            IF(WAVES(ng)%posang(i).gt.ndirsw)THEN
              WAVES(ng)%posang(i)=WAVES(ng)%posang(i)-ndirsw
            ENDIF
          ENDDO
        ELSE
          DO i=1,ndirsw
            WAVES(ng)%posang(i)=minpos-(i-1)
            IF(WAVES(ng)%posang(i).lt.1)THEN
              WAVES(ng)%posang(i)=ndirsw+WAVES(ng)%posang(i)
            ENDIF
          ENDDO
        ENDIF
!
!-----------------------------------------------------------------------
! Find exception value
!-----------------------------------------------------------------------
!
        i=0
        DO WHILE (i==0)
          READ (inp,'(a)',ERR=215,END=240) rline
          IF(rline(41:49).eq.'exception') THEN
            read(rline(1:15),*) WAVES(ng)%exc
            i=1
          END IF
        END DO
!
!-----------------------------------------------------------------------
! Get time
!-----------------------------------------------------------------------
!
        i=0
        DO WHILE (i==0)
          READ (inp,'(a)',ERR=215,END=240) rline
          IF(rline(41:44).eq.'date') THEN
            read(rline(1:15),*) WAVES(ng)%SpecTimeIso(1)
            read(rline(1:4),'(I4)') INTTIM(1)
            read(rline(5:6),'(I2)') INTTIM(2)
            read(rline(7:8),'(I2)') INTTIM(3)
            read(rline(10:11),'(I2)') INTTIM(4)
            read(rline(12:13),'(I2)') INTTIM(5)
            read(rline(14:15),'(I2)') INTTIM(6)
            WAVES(ng)%SpecTime(1)=DTTIME(INTTIM)
!  Make sure this first time is the ocean time. if not, keep reading.
!  Be within 10 seconds.
            CALL caldate (tdays(ng),yy_i=YYYY,mm_i=NN,dd_i=DD,          &
     &                    h_i=HH,m_i=MN,s_i=SS)
            INTTIM(1)=YYYY
            INTTIM(2)=NN
            INTTIM(3)=DD
            INTTIM(4)=HH
            INTTIM(5)=MN
            INTTIM(6)=SS
            OceanTime=DTTIME(INTTIM)
            IF ((OceanTime-WAVES(ng)%SpecTime(1)).lt.10.0_r8) THEN
              i=1
            END IF
          END IF
        END DO
!  Set dt SpecTime to a large valule in case the spec2d file 
!  has only 1 time instance.
        WAVES(ng)%SpecTimedt=99999999.0_r8
        WAVES(ng)%SpecTime(2)=WAVES(ng)%SpecTime(1)+WAVES(ng)%SpecTimedt
!
!-----------------------------------------------------------------------
! Find FACTOR keyword
!-----------------------------------------------------------------------
!
        i=0
        DO WHILE (i==0)
          read(inp,'(a)')rtext
          IF (rtext == 'FACTOR') THEN
            i=1
          ELSEIF (rtext == 'ZERO  ') THEN
            IF (Master) THEN
              WRITE(stdout,*) 'Zero energy density input for this point'
            END IF
            stop
          ELSEIF (rtext == 'NODATA') THEN
            IF (Master) THEN
              WRITE(stdout,*) 'SWAN file has no data for this point'
            END IF
            stop
          ENDIF
        END DO
        read(inp,*)factor
!
!-----------------------------------------------------------------------
! Read S_array
!-----------------------------------------------------------------------
!
        allocate(sdsw1(nfreqsw,ndirsw))
        allocate(WAVES(ng)%SDSW(nfreqsw,ndirsw))

        DO i=1,nfreqsw
          read(inp,*)sdsw1(i,:)
        END DO
        where (sdsw1 == WAVES(ng)%exc) sdsw1=0.0_r8
      END IF
!
!-----------------------------------------------------------------------
!  If first, then read the next time and close file for now.
!  If not first, then read until the next times step and get data.
!-----------------------------------------------------------------------
!
      j=0
      DO WHILE (j==0)
        READ (inp,'(a)',ERR=215,END=245) rline 
        IF(rline(41:44).eq.'date') THEN
          read(rline(1:15),*) WAVES(ng)%SpecTimeIso(2)
          read(rline(1:4),'(I4)') INTTIM(1)
          read(rline(5:6),'(I2)') INTTIM(2)
          read(rline(7:8),'(I2)') INTTIM(3)
          read(rline(10:11),'(I2)') INTTIM(4)
          read(rline(12:13),'(I2)') INTTIM(5)
          read(rline(14:15),'(I2)') INTTIM(6)
          WAVES(ng)%SpecTime(2)=DTTIME(INTTIM)
          IF (first.eq.1) THEN
!  Compute the dt SpecTime if we do get a second time stamp.
            WAVES(ng)%SpecTimedt=WAVES(ng)%SpecTime(2)-                 &
     &                           WAVES(ng)%SpecTime(1)
            j=1
            IF (Master) THEN
              WRITE(stdout,*) 'Using 2D spec data at ',                 &
     &                         WAVES(ng)%SpecTimeIso(1)
            END IF
          ELSE
!  If not first, then keep reading until we get the next time stamp.
!  Be within 10 seconds.
            CALL caldate (tdays(ng),yy_i=YYYY,mm_i=NN,dd_i=DD,          &
     &                    h_i=HH,m_i=MN,s_i=SS)
            INTTIM(1)=YYYY
            INTTIM(2)=NN
            INTTIM(3)=DD
            INTTIM(4)=HH
            INTTIM(5)=MN
            INTTIM(6)=SS
            OceanTime=DTTIME(INTTIM)
            IF ((OceanTime-WAVES(ng)%SpecTime(2)).lt.10.0_r8) THEN
              read(inp,'(a)')rtext
              IF (rtext == 'ZERO  ') THEN
                IF (Master) THEN
                  WRITE(stdout,*) 'Zero energy density input'
                END IF
                stop
              ELSEIF (rtext == 'NODATA') THEN
                IF (Master) THEN
                  WRITE(stdout,*) 'SWAN file has no data'
                END IF
                stop
              ELSEIF (rtext == 'FACTOR') THEN
                read(inp,*) factor
                allocate(sdsw1(WAVES(ng)%nfreqsw,WAVES(ng)%ndirsw))
                DO i=1,WAVES(ng)%nfreqsw
                  read(inp,*) sdsw1(i,:)
                END DO
                where (sdsw1 == WAVES(ng)%exc) sdsw1=0.0_r8
                IF (Master) THEN
                  WRITE(stdout,*) 'Using 2D spec data at ',             &
     &                             WAVES(ng)%SpecTimeIso(2)
                END IF
!
!  Read next time step and update SpecTimeIso of 2.
!  Temporarily set speciso to be +dt in case there is no more data.
!
                WAVES(ng)%SpecTime(2)=WAVES(ng)%SpecTime(2)+            &
     &                                WAVES(ng)%SpecTimedt
                READ (inp,'(a)',ERR=215,END=245) rline 
                IF(rline(41:44).eq.'date') THEN
                  read(rline(1:15),*) WAVES(ng)%SpecTimeIso(2)
                  read(rline(1:4),'(I4)') INTTIM(1)
                  read(rline(5:6),'(I2)') INTTIM(2)
                  read(rline(7:8),'(I2)') INTTIM(3)
                  read(rline(10:11),'(I2)') INTTIM(4)
                  read(rline(12:13),'(I2)') INTTIM(5)
                  read(rline(14:15),'(I2)') INTTIM(6)
                  WAVES(ng)%SpecTime(2)=DTTIME(INTTIM)
                END IF
                j=1
              END IF
            END IF
          END IF
        END IF
      END DO
 245  CONTINUE
      close(inp)
!
!-----------------------------------------------------------------------
! Convert to m2/Hz/rad
!-----------------------------------------------------------------------
!
      DO j=1,WAVES(ng)%ndirsw
        DO i=1,WAVES(ng)%nfreqsw
          sdsw1(i,j)=sdsw1(i,j)*factor*orad
        END DO
      END DO
!
!-----------------------------------------------------------------------
! Convert from energy density to variance density and order
!-----------------------------------------------------------------------
!
      DO j=1,WAVES(ng)%ndirsw
        IF (first.eq.1) WAVES(ng)%theta(j)=theta1(WAVES(ng)%posang(j))
        DO i=1,WAVES(ng)%nfreqsw
         WAVES(ng)%SDSW(i,j)=sdsw1(i,WAVES(ng)%posang(j))/(rho0*g)
        ENDDO
      ENDDO

      deallocate(sdsw1)
      IF (first.eq.1) deallocate(theta1)
!     deallocate(posang)
!
!-----------------------------------------------------------------------
! Integrate the directional spectra over the directions
!-----------------------------------------------------------------------
!
      IF (first.eq.1) allocate(WAVES(ng)%Sf(WAVES(ng)%nfreqsw))
      
      WAVES(ng)%Sf = sum(WAVES(ng)%SDSW, DIM = 2)*WAVES(ng)%dang
!
!-----------------------------------------------------------------------
! Find main wave direction
!-----------------------------------------------------------------------
!
      allocate (temp(WAVES(ng)%ndirsw))
!
!-----------------------------------------------------------------------
! Integrate the directional spectra over the frequencies
!-----------------------------------------------------------------------
!
      temp=sum(WAVES(ng)%SDSW, DIM = 1)
      i=maxval(maxloc(temp))
      WAVES(ng)%mainang=WAVES(ng)%theta(i)
      WAVES(ng)%mainangi=i
      IF (Master) THEN
        WRITE(stdout,*) 'Peak dir from ', WAVES(ng)%mainang*180.0_r8/pi
      END IF
      deallocate(temp)
!
!-----------------------------------------------------------------------
! Calculate zero-order moment
!-----------------------------------------------------------------------
!
      allocate (temp(WAVES(ng)%nfreqsw+1))
      temp(1)=0.0_r8
      DO i=2,WAVES(ng)%nfreqsw
        temp(i)=0.5_r8*WAVES(ng)%fsw(i-1)+                              &
     &          0.5_r8*WAVES(ng)%fsw(i)
      END DO
      temp(WAVES(ng)%nfreqsw+1)=WAVES(ng)%fsw(WAVES(ng)%nfreqsw)

      m0=0.0_r8
      DO i=1,WAVES(ng)%nfreqsw
        m0=m0+WAVES(ng)%Sf(i)*(temp(i+1)-temp(i))
      END DO
      deallocate (temp)
      WAVES(ng)%hm0gew=4.004_r8*sqrt(m0)
      IF (Master) THEN
        WRITE(stdout,*) 'Waves hm0 is ', WAVES(ng)%hm0gew
      END IF
!
!-----------------------------------------------------------------------
! Compute the representative period
!-----------------------------------------------------------------------
!
      call tpdcalc(WAVES(ng)%Sf,WAVES(ng)%fsw,frep)
!
      WAVEG(ng)%Trep=1.0_r8/frep
      WAVES(ng)%fp=frep
      IF (Master) THEN
        WRITE(stdout,*) 'Representative period is ',  WAVEG(ng)%Trep
      END IF
!
      GO TO 250
 215  IF (Master) THEN
        WRITE (stdout,*) 'swan_reader.F Error reading swan 2d spec file'
      END IF
!     exit_flag=4
 240  IF (Master) THEN
        WRITE(stdout,*) '2d specfile END: NO EXC val or incorrect time'
      END IF
 250  CONTINUE
#endif
      RETURN
      END SUBROUTINE swan_reader

!*******************************************************************
      FUNCTION DTTIME (INTTIM)
!                                                                  *
!     Borrowed from SWAN ocpmix.F DTTIME                           *
!*******************************************************************
!
      USE mod_scalars

      IMPLICIT NONE

      real(r8) ::  DTTIME
!
      INTEGER INTTIM(6)
!
      INTEGER IDYMON(12), IYEAR, IYRM1, IDNOW, I, II
!
      LOGICAL LEAPYR, LOGREF
!
      DATA IDYMON /31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/
!
      IYEAR = INTTIM(1)
      IYRM1 = IYEAR-1
      LEAPYR=(MOD(IYEAR,4).EQ.0.AND.MOD(IYEAR,100).NE.0).OR.            &
     &        MOD(IYEAR,400).EQ.0
      IDNOW=0
      DO I = 1,INTTIM(2)-1
        IDNOW=IDNOW+IDYMON(I)
      END DO
      IDNOW=IDNOW+INTTIM(3)
      IF (LEAPYR.AND.INTTIM(2).GT.2) IDNOW=IDNOW+1
      IDNOW = IDNOW + IYEAR*365 + IYRM1/4 - IYRM1/100 + IYRM1/400 + 1
      IF (IYEAR.EQ.0) IDNOW=IDNOW-1
!      IF (.NOT.LOGREF) THEN
!        REFDAY = IDNOW
!        LOGREF = .TRUE.
!        DTTIME = 0.
!      ELSE
!        DTTIME = REAL(IDNOW-REFDAY) * 24.*3600.
        DTTIME = REAL(IDNOW) * 24.*3600.
!      ENDIF
      DTTIME = DTTIME + 3600.*REAL(INTTIM(4)) + 60.*REAL(INTTIM(5)) +   &
     &                  REAL(INTTIM(6))
      RETURN
      END
