#include "cppdefs.h"
!
!************************************************************************
      SUBROUTINE swan_reader (ng)
!************************************************************************
!
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
      USE mod_scalars
      USE mod_inwave_swan

!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng

!
!  Local variable declarations.
!
      character(6)       :: rtext
      real(r8)           :: factor,exc,m0
      integer            :: nfreq, ndir, switch
      integer            :: i, j
      integer            :: nfreqsw, ndirsw
      integer            :: minpos
      integer,allocatable       :: posang(:)

      real(r8), pointer  :: temp(:), findline(:)
      real(r8), pointer  :: tempA(:,:)
      real(r8), pointer  :: sdsw1(:,:)
      real(r8), pointer  :: theta1(:)

      real(r8), parameter :: rad = pi/180.0_r8
      real(r8), parameter :: irad = 180.0_r8/pi

      integer :: ND_BND

!
!-----------------------------------------------------------------------
! Open the spectral file
!-----------------------------------------------------------------------
!
      open(44,file=IWSWNname(ng),form='formatted',status='old')

      switch=0
!
!-----------------------------------------------------------------------
! Read file until RFREQ or AFREQ is found
!-----------------------------------------------------------------------
!
      DO while (switch==0)
        read(44,'(a)')rtext
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
      read(44,*)nfreqsw
      allocate(WAVES(ng)%fsw(nfreqsw))
      DO i=1,nfreqsw
        read(44,*)WAVES(ng)%fsw(i)
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
      read(44,'(a)')rtext
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
      read(44,*)ndirsw
      allocate(theta1(ndirsw))
      allocate(WAVES(ng)%theta(ndirsw))

      DO i=1,ndirsw
        read(44,*)theta1(i)
        IF(theta1(i).lt.0.0_r8)THEN
          theta1(i)=(360.0_r8+theta1(i))
        ENDIF
        theta1(i)=theta1(i)*rad
      ENDDO

      WAVES(ng)%ndirsw=ndirsw
      WAVES(ng)%dang=abs(theta1(2)-theta1(1))
      minpos=minloc(theta1, DIM = 1)
      allocate(posang(ndirsw))

      IF((theta1(minpos+1)-theta1(minpos)).lt.                          &
     &   (theta1(minpos-1)-theta1(minpos)))THEN
        DO i=1,ndirsw
          posang(i)=minpos+(i-1)
          IF(posang(i).gt.ndirsw)THEN
            posang(i)=posang(i)-ndirsw
          ENDIF
        ENDDO
      ELSE
        DO i=1,ndirsw
          posang(i)=minpos-(i-1)
          IF(posang(i).lt.1)THEN
            posang(i)=ndirsw+posang(i)
          ENDIF
        ENDDO
      ENDIF

      read(44,'(a)')rtext
      read(44,*)exc
      i=0
!
!-----------------------------------------------------------------------
! Find FACTOR keyword
!-----------------------------------------------------------------------
!
      DO WHILE (i==0)
        read(44,'(a)')rtext
        IF (rtext == 'FACTOR') THEN
          i=1
        ELSEIF (rtext == 'ZERO  ') THEN
          write(*,*) 'Zero energy density input for this point'
          stop
        ELSEIF (rtext == 'NODATA') THEN
          write(*,*) 'SWAN file has no data for this point'
          stop
        ENDIF
      END DO
      read(44,*)factor
!
!-----------------------------------------------------------------------
! Read S_array
!-----------------------------------------------------------------------
!
      allocate(sdsw1(nfreqsw,ndirsw))
      allocate(WAVES(ng)%SDSW(nfreqsw,ndirsw))
     
      DO i=1,nfreqsw
        read(44,*)sdsw1(i,:)
      END DO
      where (sdsw1 == exc) sdsw1 =0.0_r8
!
!-----------------------------------------------------------------------
! Convert to m2/Hz/rad
!-----------------------------------------------------------------------
!
      DO j=1,ndirsw
        DO i=1,nfreqsw
            sdsw1(i,j)=sdsw1(i,j)*factor*irad
        END DO
      END DO
      close(44)
!
!-----------------------------------------------------------------------
! Convert from energy density to variance density and order
!-----------------------------------------------------------------------
!
      DO j=1,ndirsw
        WAVES(ng)%theta(j)=theta1(posang(j))
        DO i=1,nfreqsw
          WAVES(ng)%SDSW(i,j)=sdsw1(i,posang(j))/(rho0*g)
        ENDDO
      ENDDO

      deallocate(sdsw1)
      deallocate(theta1)
!
!-----------------------------------------------------------------------
! Integrate the directional spectra over the directions
!-----------------------------------------------------------------------
!
      allocate(WAVES(ng)%Sf(nfreqsw))
      
      WAVES(ng)%Sf = sum(WAVES(ng)%SDSW, DIM = 2)*WAVES(ng)%dang
!
!-----------------------------------------------------------------------
! Find main wave direction
!-----------------------------------------------------------------------
!
      allocate (temp(ndirsw))
!
!-----------------------------------------------------------------------
! Integrate the directional spectra over the frequencies
!-----------------------------------------------------------------------
!
      temp=sum(WAVES(ng)%SDSW, DIM = 1)
      i=maxval(maxloc(temp))
      WAVES(ng)%mainang=WAVES(ng)%theta(i)
      deallocate(temp)
!
!-----------------------------------------------------------------------
! Calculate zero-order moment  
!-----------------------------------------------------------------------
!
      allocate (temp(nfreqsw+1))
      temp(1)=0.0_r8
      DO i=2,nfreqsw 
        temp(i)=0.5_r8*WAVES(ng)%fsw(i-1)+                              &
     &          0.5_r8*WAVES(ng)%fsw(i)
      END DO
      temp(nfreqsw+1)=WAVES(ng)%fsw(nfreqsw)

      m0=0.0_r8
      DO i=1,nfreqsw
        m0=m0+WAVES(ng)%Sf(i)*(temp(i+1)-temp(i))
      END DO
      deallocate (temp)
      WAVES(ng)%hm0gew=4.0_r8*sqrt(m0)
!
!-----------------------------------------------------------------------
! Compute the representative period
!-----------------------------------------------------------------------
!
      call tpdcalc(WAVES(ng)%Sf,WAVES(ng)%fsw,WAVES(ng)%Trep)
      WAVES(ng)%fp=WAVES(ng)%Trep
      WAVES(ng)%Trep=1.0_r8/WAVES(ng)%Trep

#endif
      RETURN
      END SUBROUTINE swan_reader
