      program convrt2d
C
C     Update: Corrected the count in one of the do-loops (3-11-99)
C
C     PURPOSE
C
C        This program converts a pre-SWAN-40.00
C        2D-spectrum file to a file suitable for
C        SWAN version 40.51.
C
C     USAGE OF THIS PROGRAM
C
C     the program will work interactively
C     it will ask the user to provide the following data:
C     1: name of the input file (file containing spectra in old format)
C     2: type of this file, i.e. nest/spec2d, stationary/nonstationary
C     3: name of the output file (file containing spectra in new format)
C
C     In the case of a nesting file the program will write computational
C     grid data to the screen.
C
      IMPLICIT NONE
C
      integer   maxloc, maxfrq, maxdir
      parameter (maxloc=5000, maxfrq=100, maxdir=360)
      integer   id, ifr, iver, intype, outtyp, itim, iloc, time,
     &          ndirs, nfreqs, numl, numq, mxn, myn
      real      dir(maxdir), fact, freqhz(maxfrq), pi, tt, efac,
     &          xnlen, ynlen, alpcn, slow, shig, sector, exc
      REAL      vadens(maxdir,maxfrq)
      DOUBLE PRECISION dvadens(maxdir,maxfrq), xpcn, ypcn,
     &          locx(maxloc), locy(maxloc)
      character infile *40, outfil *40, chtime *20
      logical   first
C
      iver = 1
      exc  = -99.
      pi   = 4.*atan(1.)
      numl = 1
      numq = 1
      itim = 0
      if (intype.eq.2 .or. intype.eq.4) itim = 1
C
C     Old SWAN data files had the Variance denisity calculated per radian
C     The new files store it per degree. The factor pi/180 is the conversion
C     factor.
C
      fact = pi/180
C
      WRITE(6,*) '      -----------------------------------------------'
      WRITE(6,*) '      |   This program converts a pre-SWAN-40.00    |'
      WRITE(6,*) '      |   2D-spectrum file to a file suitable for   |'
      WRITE(6,*) '      |   SWAN version 40.51                        |'
      WRITE(6,*) '      -----------------------------------------------'
      WRITE(6,*) ' '
C
      write (6,*) 'Give the name of the input file: '
      read  (5, '(a)') infile
      write (6,*) 'Give type: 1=old Swan stat, 2=old Swan nonstat,'
      write (6,*) '           3=old Nest stat, 4=old Nest nonstat: '
      read  (5, *) intype
      write (6,*) 'Give the name of the output file: '
      read  (5, '(a)') outfil
      outtyp = 5
C
      open (8, file=infile, status='old', form='formatted')
      open (9, file=outfil, status='unknown', form='formatted')
*
*     write heading into the file
*     write keyword SWAN and version number
*
      write (9, 101) iver
 101  format ('SWAN', I4, T41, 'Swan standard spectral file, version')
      WRITE (9,'(a35)') '$   Data produced by pre-SWAN-40.00'
      if (intype.eq.2 .or. intype.eq.4) then
        itim = 1
      endif
      if (itim.gt.0) then
        write (9, 102) 'TIME', 'time-dependent data'
        write (9, 103) itim, 'time-coding option'
      endif
 102  format (a, t41, a)
 103  format (i6, t41, a)
 104  format (f11.2,f11.2)
 105  format (e13.4, T41, A)
C
C     first read to obtain general data
C
      if (intype.eq.1 .or. intype.eq.2) then
        first = .true.
  31    if (.not.first) goto 200
*       first read to find all locations
        do iloc = 1, maxloc
          read (8, *, end=38) locx(iloc), locy(iloc)
          read (8, *) ndirs, nfreqs
          read (8, *) (dir(id), id=1,ndirs)
          read (8, *) (freqhz(ifr), ifr=1,nfreqs)
          do ifr = 1,nfreqs
            read  (8, *) (dvadens(id,ifr), id=1,ndirs)
          end do
          numl = iloc
        enddo
        read (8, *, end=38) locx(maxloc), locy(maxloc)
        write (6, *) 'too many locations on file'
        stop
  38    rewind(8)
        first = .false.
        goto 31
      else if (intype.eq.3 .or. intype.eq.4) then
        first = .true.
  41    if (intype.eq.4) then
          read (8,'(1X,A)') CHTIME
          read (8,'(1X,F8.0)') tt
        endif
        read (8, *) XNLEN, YNLEN, SLOW, SHIG
        read (8, *) MXN, MYN
        read (8, *) XPCN, YPCN, ALPCN
        write (6, *) 'CGRID Regular ', xpcn, ypcn, alpcn*180./pi,
     &                '    &', xnlen, ynlen, mxn, myn, '    &'
        read (8, *) ndirs, nfreqs
        read (8, *) (dir(id), id=1,ndirs)
        read (8, *) (freqhz(ifr), ifr=1,nfreqs)
        sector = real(ndirs) * (dir(2)-dir(1))
        if (abs(sector-360.).lt.10.) then
          write (6, *) 'circle ', ndirs,
     &                  freqhz(1), freqhz(nfreqs), nfreqs
        else
          write (6, *) 'sector ', dir(1), dir(ndirs), ndirs,
     &                  freqhz(1), freqhz(nfreqs), nfreqs
        endif
        if (.not.first) goto 200
*       first read to find all locations
        do iloc = 1, maxloc+1
          read (8, *, end=48) locx(iloc), locy(iloc)
          do ifr = 1,nfreqs
            read  (8, *) (dvadens(id,ifr), id=1,ndirs)
          end do
          numl = iloc
        enddo
        write (6, *) 'too many locations on file'
        stop
  48    rewind(8)
        first = .false.
        goto 41
      endif
*
 200  continue
      if (outtyp.eq.5) then
        write (9, 102) 'LOCATIONS', 'locations in x-y-space'
        write (9, 103) numl, 'number of locations'
        do iloc = 1, numl
          write (9, 104) locx(iloc), locy(iloc)
        enddo
C
        write (9, 102) 'AFREQ', 'absolute frequency in Hz'
        write (9, 103) nfreqs,  'number of frequencies'
        write (9, 131) (freqhz(ifr), ifr=1,nfreqs)
 131    format (f11.4)
C
        write (9, 102) 'CDIR', 'spectral Cartesian directions in degr'
        write (9, 103) ndirs,  'number of directions'
        write (9, 131) (dir(id), id=1,ndirs)
C
        WRITE (9, 102) 'QUANT'
        WRITE (9, 103) numq,         'number of quantities in table'
        WRITE (9, 102) 'VaDens',     'variance densities in m2/Hz/degr'
        write (9, 102) 'm2/Hz/degr', 'unit'
        write (9, 105) exc,   'exception value'
C
        do time= 1, 999
          if (intype.eq.4) then
            read (8,'(1X,A)', end=800) CHTIME
            read (8,'(1X,F8.0)') tt
            write (9, '(A)') chtime
          endif
          do iloc = 1, numl
            if (intype.gt.0 .and. intype.le.4) then
              read (8, *) locx(iloc), locy(iloc)
            endif
            if (intype.eq.1 .or. intype.eq.2) then
*             ignore frequencies and directions
              read (8, *) ndirs, nfreqs
              read (8, *) (dir(id), id=1,ndirs)
              read (8, *) (freqhz(ifr), ifr=1,nfreqs)
            endif
C
            efac=0.
            do ifr = 1,nfreqs
              read (8, *) (dvadens(id,ifr), id=1,ndirs)
              do id=1,ndirs
                if (intype.eq.1 .or. intype.eq.2) then
                  vadens(id,ifr) = REAL(dvadens(id,ifr))
                else if (intype.eq.3 .or. intype.eq.4) then
*                 convert action density into energy density
                  vadens(id,ifr) = REAL(dvadens(id,ifr)) *
     &                             (2.*pi*freqhz(ifr))
                endif
                if (vadens(id,ifr).ge.0.) then
                   efac = max (efac, vadens(id,ifr))
                else
                   efac = max (efac, 10.*abs(vadens(id,ifr)))
                endif
              end do
            end do
            if (efac .le. 1.e-10) then
              WRITE (9, '(a4)') 'ZERO'
            else
              efac = 1.01 * efac * 10.**(-4)
!             factor pi/180 introduced to account for change from rad to degr
!             factor 2*pi to account for transition from rad/s to hz
              write (9, 95) efac * 2. * pi**2 / 180.
  95          format ('FACTOR', /, e18.8)
              do ifr = 1, nfreqs
!                write spectral energy densities to file
                 write (9, 133) (nint(vadens(id,ifr)/efac), id=1,ndirs)
              enddo
 133          format (200(1x,i4))
            endif
          enddo
          if (itim.eq.0) goto 800
        enddo
      endif
C
 800  write (6, *) ' Conversion finished'
C
      end
C
