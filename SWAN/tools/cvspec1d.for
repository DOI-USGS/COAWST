      program cvspec1d
C
C     PURPOSE
C
C        This program converts a pre-SWAN-40.00
C        1D-spectrum file to a file suitable for
C        SWAN version 40.51.
C
C     USAGE OF THIS PROGRAM
C
C     the program will work interactively
C     it will ask the user to provide the following data:
C     1: name of the input file (file containing spectra in old format)
C     2: number of spectral frequencies
C     3: average direction, or -999. if average direction is in the file
C     4: directional spread (in degr), or -999. if spreading is in the file
C     5: number of columns in the input file (excluding the first column
C        which always contains the frequencies)
C
      parameter (mxfreq=500, mxcols=1000)
      real freqs(mxfreq), values(mxfreq,mxcols)
      character infile *40, outfil *40
C
      WRITE(6,*) '      -----------------------------------------------'
      WRITE(6,*) '      |   This program converts a pre-SWAN-40.00    |'
      WRITE(6,*) '      |   1D-spectrum file to a file suitable for   |'
      WRITE(6,*) '      |   SWAN version 40.51                        |'
      WRITE(6,*) '      -----------------------------------------------'
      WRITE(6,*) ' '
C
      inrhog = 0
      kcols = 3
C
      write (6,*) 'Give the name of the input file: '
      read  (5, '(a)') infile
      write (6,*) 'Number of frequencies: '
      read  (5,*) nfreq
      write (6,*) 'Average direction (-999. if in file): '
      read  (5,*) avdir
      write (6,*) 'Directional spreading (-999. if in file): '
      read  (5,*) dd
      write (6,*) 'Number of columns in the input file (excl freq): '
      read  (5, *) ncols
C
      if (avdir.gt.-990.) kcols = kcols-1
      if (dd.gt.0.) kcols = kcols-1
C
      write (6,*) 'Give the name of the output file: '
      read  (5, '(a)') outfil
C
      if (infile.ne.'    ') then
        iread = 8
        open (iread, file=infile, form='formatted', status='old')
      else
        iread = 5
      endif
      open (9, file=outfil, form='formatted', status='unknown')
C
      iverf = 1
      WRITE (9, 101) IVERF
 101  FORMAT ('SWAN', I4, T41, 'Swan standard spectral file, version')
      WRITE (9,'(a35)') '$   Data produced by pre-SWAN-40.00'
C
 102  FORMAT (A, T41, A)
 103  FORMAT (I6, T41, A)
 105  format (e13.4, T41, A)
      nlocs = (ncols+kcols-1) / kcols
      if (nlocs.gt.1) then
        WRITE (9, 102) 'LOCATIONS', 'locations in x-y-space'
        WRITE (9, 103) nlocs, 'number of locations'
        DO 110 IP = 1, nlocs
          WRITE (9, *) 0., 0.
 110    CONTINUE
      endif
C
      WRITE (9, 102) 'AFREQ', 'absolute frequencies in Hz'
      WRITE (9, 103) nfreq, 'number of frequencies'
      do ifreq=1, nfreq
         read (iread, *) freqs(ifreq), (values(ifreq,jj), jj=1,ncols)
         write (9, 114) FREQS(IFREQ)
 114     FORMAT (F10.4)
      enddo
C
      WRITE (9, 102) 'QUANT', ' '
      WRITE (9, 103) 3, 'number of quantities'
C
      IF (INRHOG.EQ.1) THEN
        WRITE (9, 102) 'EnDens', 'energy densities'
        write (9, 102) 'J/m2/Hz/rad', 'unit'
        write (9, 105) -99., 'exception value'
      ELSE
        WRITE (9, 102) 'VaDens', 'variance densities'
        write (9, 102) 'm2/Hz', 'unit'
        write (9, 105) -99., 'exception value'
      ENDIF
      WRITE (9, 102) 'CDIR','average Cartesian direction in degr'
      write (9, 102) 'degr', 'unit'
      write (9, 105) -999., 'exception value'
      WRITE (9, 102) 'DSPRDEGR', 'directional spread in degr'
      write (9, 102) 'degr', 'unit'
      write (9, 105) -9., 'exception value'
C
      do kk = 1, ncols, kcols
        iloc = (kk+kcols-1) / kcols
        write (9,555) iloc
 555    FORMAT('LOCATION     ',i3)
        do ifreq = 1, nfreq
          ee = values(ifreq,kk)
          if (avdir.gt.-990.) then
            aa = avdir
            jj = kk
          else
            aa = values(ifreq,kk+1)
            jj = kk+1
          endif
          if (dd.gt.0.) then
            bb = dd
          else
            bb = values(ifreq,jj+1)
          endif
          write (9,202) ee, aa, bb
 202      format (e12.5, 2f9.2)
        enddo
      enddo
      close (8)
      close (9)
      write (*,*) 'conversion finished'
      stop
      end
