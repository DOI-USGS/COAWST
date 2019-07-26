#include "cppdefs.h"
!      MODULE refdif_mod
#ifdef REFDIF_COUPLING
!
!=================================================== John C. Warner ====
!=================================================== Kevin Haas     ====
!                              ROMS/TOMS Groups                      ===
!================================================== Hernan G. Arango ===
!                                                                      !
!  This module is for coupling ROMS to REFDIF.                         !
!                                                                      !
!=======================================================================
!

!      implicit none

!      PRIVATE

!      PUBLIC :: Refdif
!      PUBLIC :: REFDIF_INITIALIZE
!      PUBLIC :: REFDIF_RUN
!      PUBLIC :: REFDIF_FINALIZE

!      integer :: ifreq, ir

!      CONTAINS

!      subroutine Refdif
      SUBROUTINE REFDIF_INITIALIZE (MyCOMM)
      USE WAVES_COUPLER_MOD, ONLY: INITIALIZE_OCEAN_COUPLING
      USE WAVES_COUPLER_MOD, ONLY: OCEAN_COUPLING
      include 'mpif.h'

!c  Define array dimension bounds.
    
!      include 'param.h'

!c  Internal common blocks.

!      include 'common.h'

!c  Common block data passing to Master program.

!      include 'pass.h'

#include "param.h"
#include "common.h"
#include "pass.h"


      integer i,j

!      IF (PRESENT(MyCOMM)) THEN
        WAV_COMM_WORLD=MyCOMM
!      ELSE
!        WAV_COMM_WORLD=MPI_COMM_WORLD
!      END IF
      ncomps=2
      OcnId=1
      WavId=2

!c --- master_start=0 or 1 for initialization

      if(Master_Start.eq.1) then
        write(*,*) 'wave module initialization ...'
        else
        write(*,*) 'call wave module ...'
      endif

!c  Constants which provide for conversion between MKS and English units
!c  on input and output.

      dconv(1)=1.
      dconv(2)=0.30488
      dconv2(1)=1.
      dconv2(2)=14.594    

!c  Read control parameters and reference grid data.
     
      call inref
     
!c  Read control parameters and initializing wave data.

      if(Master_Start.ge.0)then    
        call inwave
        close(1)
      endif

! jcw call to init MCT World
      call initialize_ocean_coupling
      call ocean_coupling
!!jcw - add 2nd call to now use depth data
!      call inref
!jcw - add computations for mct exchanged data
         do i=1,mr
           do j=1,nr
             dr(i,j)=Depth_Wave(i,j) !+intp_eta_wave(i,j)
             ur(i,j)=Intp_U_Wave(i,j)
             vr(i,j)=Intp_V_Wave(i,j)
           end do
         end do
!
      return
      end subroutine REFDIF_INITIALIZE

!
      SUBROUTINE REFDIF_RUN (TI_WAV_OCN,Wname)
      USE WAVES_COUPLER_MOD, ONLY: OCEAN_COUPLING
#include "param.h"
#include "common.h"
#include "pass.h"

!c  Pass program control to subroutine |model|.

!c  For each frequency component specified in |inwave|, |model|
!c  executes the model throughout the entire grid and then
!c  reinitializes the model for the next frequency.
    
!      if(Master_Start.le.0)then
       DO irun=1,45000.*0.04/TI_WAV_OCN
!jcw - add computations for mct exchanged data
         do i=1,mr
           do j=1,nr
             dr(i,j)=Depth_Wave(i,j)+intp_eta_wave(i,j)
             ur(i,j)=Intp_U_Wave(i,j)
             vr(i,j)=Intp_V_Wave(i,j)
           end do
         end do

         call model
         call ocean_coupling
       END DO
!      endif
      return
      end subroutine REFDIF_RUN

!
      SUBROUTINE REFDIF_FINALIZE
      USE WAVES_COUPLER_MOD, ONLY: FINALIZE_OCEAN_COUPLING

      call finalize_ocean_coupling('normal end')

!c  All done.  Close output data files if |open| and |close| statements are 
!c  being used.

!c  Log file.

      close(5)

!c  Wave height.

      if(fname6.ne.' ') close(26)

!c  Wave angle.

      if(fname7.ne.' ') close(7)

!c  Water depth.

      if(fname8.ne.' ') close(8)

!c  Water surface.

      if(fname9.ne.' ') close(9)

!c  Radiation stresses.

      if(fname10.ne.' ') then
      close(10)
      close(11)
      close(12)
      endif

!c Depth-integrated forcing.

      if(fname13.ne.' ') then
      close(13)
      close(14)
      endif

!c Total wave-induced mass flux.

      if(fname15.ne.' ') then
      close(15)
      close(16)
      endif

!c Bottom velocity moments.

      if(fname17.ne.' ') then
      close(17)
      close(18)
      endif

!c Breaking index.

      if(fname19.ne.' ') close(19)

!c Stored complex amplitude on last row.

      if(fname20.ne.' ') close(20)

!clocal radiation stresses

      if(fname21.ne.' ') close(21)
      if(fname22.ne.' ') close(22)
      if(fname23.ne.' ') close(23)
      if(fname24.ne.' ') close(24)
      if(fname25.ne.' ') close(25)
      if(fname26.ne.' ') close(27)

!c All done.


      return

      end subroutine REFDIF_FINALIZE

  
      subroutine inref
  
#include "param.h"
#include "common.h"
#include "pass.h"
!      include 'param.h'
        
!      include 'common.h'

!      include 'pass.h'

      integer i,j

      namelist/fnames/fname2,fname3,fname4,fname5,fname6,               &
     &                 fname7,fname8,fname9,fname10,fname11,fname12,    &
     &                 fname13,fname14,fname15,fname16,fname17,fname18, &
     &                 fname19,fname20,fname21,fname22,fname23,fname24, &
     &                 fname25,fname26

      namelist/ingrid/mr, nr, iu, ntype, icur, ibc, ismooth, dxr,       &
     &                dyr, dt, ispace, nd, iff, isp, iinput, ioutput

      namelist/inmd/  md
 
      if(Master_Start.ge.0)then

!c Constants.
  
      g=9.80621

!c Specify name of namelist data file.

      fname1='c:/work/models/roms/roms_sed_rutgers_cygwin/jcw_branch/Projects/Rip_current/indat.dat'
      open(unit=1,file=fname1)
      read(1,nml=fnames)

!c Open standard input and output files.

!c     open(unit=2,file=fname2)

      endif
!c------------ skip above after the first call wave module

      open(unit=5,file=fname5)

!c Open optional output files.

!c Store wave height?

      if(fname6.ne.' ') open(26, file=fname6)

!c Store wave angle?

      if(fname7.ne.' ') open(7, file=fname7)

!c Store water depth?

      if(fname8.ne.' ') open(8, file=fname8)

!c Store water surface?

      if(fname9.ne.' ') open(9, file=fname9)

!c Store depth-integrated radiation stresses?

      if(fname10.ne.' ') then
      open(10, file=fname10)
      open(11, file=fname11)
      open(12, file=fname12)
      endif

!c Store depth-integrated forcing?

      if(fname13.ne.' ') then
      open(13, file=fname13)
      open(14, file=fname14)
      endif

!c Store total wave-induced mass flux?

      if(fname15.ne.' ') then
      open(15, file=fname15)
      open(16, file=fname16)
      endif

!c Store bottom velocity moments?

      if(fname17.ne.' ') then
      open(17, file=fname17)
      open(18, file=fname18)
      endif

!c Store breaking index?

      if(fname19.ne.' ') open(19, file=fname19)

!c Store complex amplitude on last row?

      if(fname20.ne.' ') open(20, file=fname20)

!cStore radiation stresses for 3d circulation model

      if(fname21.ne.' ') open(21, file=fname21)
      if(fname22.ne.' ') open(22, file=fname22)
      if(fname23.ne.' ') open(23, file=fname23)
      if(fname24.ne.' ') open(24, file=fname24)
      if(fname25.ne.' ') open(25, file=fname25)
      if(fname26.ne.' ') open(27, file=fname26)

!c print headers on output
    
      write(5,120)
        
      write(5,106)
  
      if(Master_Start.ge.0)then

!c Read control data from unit 1.

      read(1,nml=ingrid)

  
      if(ispace.eq.1) read(1,nml=inmd)
!       mr=74
!       nr=92
!       mr=Nx_Wave
!       nr=Ny_Wave
        Nx_Wave=mr
        Ny_Wave=nr
!       dxr=X_Wave(2,1)-X_Wave(1,1)
!       dyr=Y_Wave(1,2)-Y_Wave(1,1)

      write(5,107) mr,nr,dxr,dyr

      if(iu.eq.1) write(5,114) iu
      if(iu.eq.2) write(5,115) iu

      if(icur.eq.0) write(5,200) 
      if(icur.eq.1) write(5,201)

      if(ibc.eq.0) write(5,202)
      if(ibc.eq.1) write(5,203)

      if(ispace.eq.0)write(5,108)
      if(ispace.eq.1)write(5,109)

      write(5,119) nd

      if(ntype.eq.0) write(5,110)
      if(ntype.eq.1) write(5,111)
      if(ntype.eq.2) write(5,112)
  
!c Check input from unit 1.
  
      if((mr.gt.ixr).or.(nr.gt.iyr)) then
        write(5,*) 'dimensions for reference grid too large, stopping'
!c       call exit(1) ! not good for intel compiler
        stop
      end if

      if((iu.ne.1).and.(iu.ne.2)) iu=1

      dt=dt*dconv(iu)
      dxr=dxr*dconv(iu)
      dyr=dyr*dconv(iu)

      if(dt.eq.0.) dt=2.

      if(nd.gt.ifix(float(iy-1)/float(nr-1))) then
        write(5,102)
!c       call exit(1) ! not good for intel compiler
        stop
      endif

      if(ispace.eq.1) then
        test=0.
        do 1 i=1,mr-1
        if(md(i).gt.(ix-1)) then
          write(5,103) i
          test=1.
        endif
 1      continue
      if(test.eq.1.) then
!c      call exit(1) ! not good for intel compiler
        stop
      endif
      endif

      endif
!c--------------- skip above after the first call wave module
!c     write(*,*)'nxwave = ',nx_wave

      do 7 i=1,mr
      do 7 j=1,nr
        dr(i,j)=Depth_Wave(i,j)+intp_eta_wave(i,j)
        ur(i,j)=Intp_U_Wave(i,j)
        vr(i,j)=Intp_V_Wave(i,j)
 7    continue
 
!c Check for large depth changes and large currents in reference grid data.
  
      do 8 i=2,mr-1
      do 8 j=2,nr-1
      dcheck=(dr(i+1,j)+dr(i-1,j)+dr(i,j-1)+dr(i,j+1))/4.
      if(abs(dcheck-dr(i,j)).gt.10) write(5,104) dr(i,j), i,j,dt
  8   continue

      if(icur.eq.1) then
      do 9 i=1,mr
      do 9 j=1,nr
      if(dr(i,j).le.0.0) go to 9
      fr=(ur(i,j)*ur(i,j)+vr(i,j)*vr(i,j))/(g*dr(i,j))
      if(fr.gt.1.) write(5,105) i,j,fr
 9    continue

      endif
    
!c Establish coordinates for reference grid.
    
      do 10 ir=1,mr
      xr(ir)=dfloat(ir-1)*dxr
 10   continue

      do 11 jr=1,nr
      yr(jr)=dfloat(jr-1)*dyr
 11   continue
    
!c Establish |y| coordinates for interpolated grid.
    
      n=nd*(nr-1)+1
      dy=dyr/dfloat(nd)

      do 12 j=1,n
      y(j)=dfloat(j-1)*dy
 12   continue
       
!c Check friction values.
 
!c        |iff(1)=1|, turbulent boundary layer damping everywhere

!c        |iff(2)=1|, porous bottom damping everywhere

!c        |iff(3)=1|, laminar boundary layer damping everywhere
      
      do 13 i=1,3
      if((iff(i).ne.0).and.(iff(i).ne.1)) iff(i)=0
 13   continue

      write(5,116)(iff(i),i=1,3)

!c Specify whether or not user specified subgrids are to be read in during 
!c model operation.

!c    |isp=0|, no subgrids specified

!c    |isp=1|, subgrids to be read in later from unit 3.
      
      if(isp.eq.0) write(5,117)
      
      if(isp.eq.1) then
      write(5,118)
      open(unit=3,file=fname3)
      endif

      if((isp.eq.1).and.(ispace.eq.0))write(5,113)

      if(isp.eq.0)then

      do 14 ir=1,mr
      do 14 jr=1,nr
      isd(ir,jr)=0
 14   continue

      else

!C CHECK UNIT NUMBER HERE>  SUBDAT NEEDS TO BE OPENED< TOO>

      do 15 ir=1,mr-1
      read(3,100)(isd(ir,jr),jr=1,nr-1)
 15   continue

      endif
    
!c Input done, return to main program.
    
      return

 100  format(15i4)
 101  format(501(f10.4))
 102  format(' y-direction subdivision too fine.'/                      &
     &' maximum number of y grid points will be exceeded.'/             &
     &' execution terminating.')
 103  format(' x-direction subdivision too fine on grid block'          &
     &,2x,i3/' execution terminating')
 104  format(' depth',2x,f7.2,'(m) at reference grid location',         &
     &2(2x,i3)/' differs from the average of its neighbors by',         &
     &' more than',2x,f7.2,'(m).'/' execution continuing')
 105  format(' ambient current at reference grid location'              &
     &,2(2x,i3),' is supercritical with froude number =',f7.4/          &
     &' execution continuing')
 106  format('0',/,/,/,20x,'input section, reference grid values',/,/,/)
 107  format(' reference grid dimensions  mr=',i3/                      &
     &       '                            nr=',i3,/,/,/,                &
     &       ' reference grid spacings   dxr=',f8.4,/,                  &
     &       '                           dyr=',f8.4)
 108  format(' '/' ispace =0 chosen, program will attempt its own ',    &
     &'reference grid subdivisions')
 109  format(' '/' ispace =1 chosen, subdivision spacings will be',     &
     &' input as data')
 110  format(' '/' ntype = 0, linear model')
 111  format(' '/' ntype = 1, stokes model matched to hedges model')
 112  format(' '/' ntype = 2, stokes model')
 113  format(' warning: input specifies that user will be supplying',   &
     &' specified subgrids (isp=1),'/                                   &
     &' while program has been told to generate its own subgrid',       &
     &' spacings (ispace=0).'/                                          &
     &' possible incompatibility in any or all subgrid blocks')
 114  format(' '/' physical unit switch iu=',i1,                        &
     &',  input in mks units')
 115  format(' '/' physical unit switch iu=',i1,                        &
     &',  input in english units')
 116  format(' ',/,/,'   switches for dissipation terms',/,/,           &
     &' ',i1,'   turbulent boundary layer'/                             &
     &' ',i1,'   porous bottom'/                                        &
     &' ',i1,'   laminar boundary layer')
 120  format(/,/,/,/,/,/,20x,'Refraction-Diffraction Model for',/,      &
     &20x,'Weakly Nonlinear Surface Water Waves',/,/,/,                 &
     &20x,'REF/DIF 1,  Version 2.6, March 2002',/,/,/,                  &
     &20x,'Center for Applied Coastal Research',/                       &
     &20x,'Department of Civil Engineering',/,                          &
     &20x,'University of Delaware',/,                                   &
     &20x,'Newark, Delaware 19716',/,/,/,                               &
     &10x,'James T. Kirby, Robert A. Dalrymple and Fengyan Shi',/,/,/,  &
     &20x,'Copyright (C) 2002 James T. Kirby',/,                        &
     &20x,'REF/DIF 1 comes with ABSOLUTELY NO WARRANTY,',/,             &
     &20x,'and is copyrighted under the provisions of the GNU',/,       &
     &20x,'General Public License.')

 117  format(' '/' isp=0, no user defined subgrids')
 118  format(' '/' isp=1, user defined subgrids to be read')
 119  format(' '/' y-direction subdivision according to nd=',i3)
 200  format(' '/' icur=0, no current values read from input files')
 201  format(' '/' icur=1, current values read from data files')
 202  format(' '/' ibc=0, closed (reflective) lateral boundaries')
 203  format(' '/' ibc=1, open lateral boundaries')

      end subroutine

    

      subroutine inwave

!      include 'param.h'

!      include 'common.h'

!      include 'pass.h'
#include "param.h"
#include "common.h"
#include "pass.h"

      namelist /waves1a/ iwave, nfreqs
      namelist /waves1b/ update_interval,num_data,                      &
     &          freqs, tide, nwavs, amp, dir
      namelist /waves1c/ thet0, freqs, tide, edens, nwavs, nseed     
      namelist /waves2/ freqin, tidein

       real*8 dis(numdata)

       integer i,j

      pi=3.1415927
    
!c Values of |iinput|, |ioutput| already entered in namelist statement in 
!c |inref|.
    
      if((iinput.ne.1).and.(iinput.ne.2))then
      write(5,*)' invalid value chosen for iinput, check indat.dat'
      stop
      endif

      if((ioutput.ne.1).and.(ioutput.ne.2))then
      write(5,*)' invalid value chosen for ioutput, check indat.dat'
      stop
      endif

      if(iinput.eq.1)then
      write(5,*)' iinput = 1, program specifies initial row of a'
    
!c Enter |iwave|, |nfreqs| for |iinput = 1|.

      read(1, nml=waves1a)
    
      write(5,102)

!c Enter data for case of |iinput=1, iwave=1|.

      if(iwave.eq.1) then 

      read(1, nml=waves1b)
      write(5,103)

      endif

!c Enter data for case of |iinput=1, iwave=2|.

      if(iwave.eq.2) then
    
      read(1, nml=waves1c)
      write(5,104)

      endif
      write(5,105) nfreqs

      if(iwave.eq.2) then

      thet0=thet0*pi/180.
      endif

!c For each frequency, enter the wave period and tidal offset.
    
      do 3 ifreq=1,nfreqs

      do itime=1,num_data

      write(5,107) ifreq,freqs(ifreq,itime),tide(ifreq)

!c Convert angles to radians.

      freqs(ifreq,itime)=2.*pi/freqs(ifreq,itime)
      tide(ifreq)=tide(ifreq)*dconv(iu)
    
!c If |iwave = 1|, read the number of discrete components.
    
      if(iwave.eq.1) then

      do 1 iwavs=1,nwavs(ifreq)

      write(5,106) iwavs,amp(ifreq,itime),dir(ifreq,itime)

      dir(ifreq,itime)=dir(ifreq,itime)*pi/180.
      amp(ifreq,itime)=amp(ifreq,itime)*dconv(iu)
 1    continue

      endif

      enddo
    
!c If |iwave = 2|, read the parameters for each frequency.
    
      if(iwave.eq.2)then

      seed=dfloat(nseed)/9999.
      write(5,108) edens(ifreq),nwavs(ifreq),nseed
      dir(ifreq,1)=thet0
      edens(ifreq)=edens(ifreq)*(dconv(iu)**2.)
      endif

 3    continue

      endif
    
!c If |iinput = 2|, read in wave period and tidal offset.
    
      if(iinput.eq.2)then

      read(1,nml=waves2)

      freqs(1,1)=freqin
      tide(1)=tidein

      write(5,*)' iinput = 2, user specifies a in wave.dat'
      nfreqs=1

      write(5,102)
      write(5,*)' wave period =',freqs(1,1),' sec.'
      write(5,*)' tidal offset=',tide(1)
      freqs(1,1)=2.*pi/freqs(1,1)
      tide(1)=tide(1)*dconv(iu)
      endif

!c--- re-arrange freqs,amp, and angle

        num_total=int(total_time/(Master_dt*N_Interval_CallWave))

        if(num_total.lt.1) num_total=1
    
        print*,'total number of wave updated =', num_total

      if(num_total.ge.numdata)then
       print*,'you should make a larger numdata in param.h'
      endif

      do ire=1,num_data
        dis(ire)=freqs(1,ire)
      enddo 

      do ire=1,num_total
        between=(Master_dt*N_Interval_CallWave)*ire                     &
     &               /UPDATE_INTERVAL
        partial=between-int(between)
        nstart=int(between)+1
        if(nstart.lt.num_data)then
        freqs(1,ire)=dis(nstart)*(1.-partial)+dis(nstart+1)*partial
        else
        freqs(1,ire)=dis(num_data)
        endif
      enddo


!c-- amp

      do ire=1,num_data
        dis(ire)=amp(1,ire)
      enddo 

      do ire=1,num_total
        between=(Master_dt*N_Interval_CallWave)*ire                     &
     &               /UPDATE_INTERVAL
        partial=between-int(between)
        nstart=int(between)+1
        if(nstart.lt.num_data)then
        amp(1,ire)=dis(nstart)*(1.-partial)+dis(nstart+1)*partial
        else
        amp(1,ire)=dis(num_data)
        endif
      enddo


!c-- angle

      do ire=1,num_data
        dis(ire)=dir(1,ire)
      enddo 

      do ire=1,num_total
        between=(Master_dt*N_Interval_CallWave)*ire                     &
     &               /UPDATE_INTERVAL
        partial=between-int(between)
        nstart=int(between)+1
        if(nstart.lt.num_data)then
        dir(1,ire)=dis(nstart)*(1.-partial)+dis(nstart+1)*partial
        else
        dir(1,ire)=dis(num_data)
        endif
      enddo

      return

 100  format(15i4)
 101  format(501(f10.4))
 102  format('1',/,/,/,20x,' input section, wave data values',/,/,/)
 103  format(' ',/,/,/,' iwave=1, discrete wave amps and directions')
 104  format(' ',/,/,/,' iwave=2, directional spreading model chosen')
 105  format(' ',/,/,/,' the model is to be run for',i3,' separate',    &
     &' frequency components')
 106  format(' '/' wave component ',i2,', amplitude =',f8.4,            &
     &', direction=',f8.4)
 107  format(' ',/,/,' frequency component ',i2,/,/,                    &
     &' wave period=',f8.4,'sec., tidal offset=',f8.4)
 108  format(' ',/,' total variance density =',f8.4,', spreading factor &
     &n=',i2,' seed number =',i5)

      end subroutine
  
 
      subroutine model

!      include 'param.h'

!      include 'common.h'

!      include 'pass.h'
#include "param.h"
#include "common.h"
#include "pass.h"

      integer i,j

      real*8 dthi(31),thi(31),thet(iy),akx,aky,akx1,aky1

!C DON't THINK THESE ARE NEEDED

      real*8 sxy(iy),sxx(iy),syy(iy)
  
!c Constants.

      g=9.80621
      rho=1000.
      pi=3.1415927
      eps=1.0e-05

!c find updated wave information
     
!      itime=max(1.,Time_Master/(N_Interval_CallWave*Master_dt))
!      itime=min(itime,num_total)
!jcw with kh set itime = 1
      itime=1
      print*,'wave itime =', itime, '  height=',2*amp(1,itime)
   
!c Execute model once for each frequency.

!c    |ifreq| is the controlling index value.
    
      do 200 ifreq=1,nfreqs
      psibar=0.
      write(5,203) ifreq

      period=2.*pi/freqs(ifreq,itime)

!c If |ismooth=1|, compute window width.

      if(ismooth.eq.1)then

      al0=g*period*period/(2.*pi)
      nrs=int((20.*al0/dyr-1.)/2.)
      write(5,*)' window half-width = ', nrs

      endif
    
!c Specify initial nonlinear parameters for each run.
    
      if(ntype.eq.0) an=0.
      if(ntype.ne.0) an=1.
      if(ntype.ne.2) anl=0.
      if(ntype.eq.2) anl=1.

!c Calculate the mean |kb| on the first row, for use in 
!c specifying initial conditions.

      npts=0
      sumk=0.

      do 10 jr=1,nr
      d(1,jr)=dr(1,jr)+tide(ifreq)!+Intp_eta_Wave(1,jr)
      call wvnum(d(1,jr),ur(1,jr),freqs(ifreq,itime),k(1,jr),eps,       &
     &icdw,1,1)

      if(d(1,jr).gt.0.05) then
        sumk=sumk+k(1,jr)
        npts=npts+1
      endif

!c--- pass wave number the first row

      Pass_WaveNum(1,jr)=k(1,jr)

!c--- pass c cg on the first row

      sig(1,jr)=freqs(ifreq,itime)-k(1,jr)*ur(1,jr)
      akd=k(1,jr)*dr(1,jr)
      q(1,jr)=(1.+akd/(sinh(akd)*cosh(akd)))/2.
      p(1,jr)=q(1,jr)*g*tanh(akd)/k(1,jr)

      Pass_C(1,jr)=sig(1,jr)/k(1,jr)
      Pass_Cg(1,jr)=sqrt(p(1,jr)*q(1,jr))

 10   continue

      kb(1)=sumk/dfloat(npts)
    
!c Establish initial wave conditions for the |ifreq| frequency
    
      if(iinput.eq.1)then
    
!c Compute wave from data given in |indat.dat|.
    
      if(iwave.eq.1) then
      
!c |iwave.eq.1|, discrete components specified.
   
      do 3 j=1,n
      a(1,j)=dcmplx(0.,0.)

      do 2 iwavs=1,nwavs(ifreq)
      thet(j)=dir(ifreq,itime)*180./pi
      a(1,j)=a(1,j)+amp(ifreq,itime)*cdexp(dcmplx(0.,kb(1)*sin(         &
     &dir(ifreq,itime))*y(j)))
 2    continue

 3    continue

      if(fname7.ne.' ') write(7,202)(thet(j),j=1,n,nd)

!c-- pass angle -- Fengyan 02/04/2002
        do j=1,n,nd
        jj=(j-1)/nd+1
        Pass_Theta(1,jj)=thet(j)
        enddo

      else
      
!c |iwave.eq.2|, directional spreading model (not recommended).
      
      sp=dfloat(nwavs(ifreq))
      nsp=nwavs(ifreq)
      thmax=pi/4.
      call acalc(thmax,nsp,a1)
      edens(ifreq)=sqrt(edens(ifreq)/a1)
      nn=31
      ii=(nn-1)/2+1
!      seed=rand1(seed)
      seed=rand1(seed)     !jcw
    
!c Compute randomly distributed $\Delta\theta$'s.
    
      sum0=0.
      do 12 i=1,nn
      seed=rand1(seed)
      dthi(i)=seed
      sum0=sum0+seed
 12   continue
      xnorm=2.*thmax/sum0
      do 101 i=1,nn
        dthi(i)=dthi(i)*xnorm
 101  continue
      thi0=-thmax
      do 4 i=1,nn
      thi0=thi0+dthi(i)
      thi(i)=thi0-dthi(i)/2.
      dth=dthi(i)
      amp(ifreq,i)=edens(ifreq)*sqrt(dth)*sqrt(cos(thi(i)+dth/2.)       &
     &**(2*nsp)+cos(thi(i)-dth/2.)**(2*nsp))
 4    continue
      do 5 i=1,nn
      ip1=i+1
      seed=rand1(seed)
      dir(ifreq,ip1)=2.*pi*seed/100.
 5    continue
      do 7 j=1,n
      a(1,j)=dcmplx(0.,0.)
      do 6 i=1,nn
      a(1,j)=a(1,j)+amp(ifreq,i)*cdexp(dcmplx(0.,kb(1)*sin(thi(i)-thet0)&
     &*y(j)+dir(ifreq,i+1)))*2.
 6    continue
 7    continue
      endif
      endif
    
!c If |iinput=2|, read |a| from data file |fname4|.
    
      if(iinput.eq.2)then
      open(4,file=fname4)
      read(4,*)(a(1,j),j=1,n)
      close(4)

!c Calculate wave angle on the first row.

!CCC NO - THIS DOESN'T WORK - IT HAS A BACKWARDS DERIVATIVE


!CCC NEW -- Tony add this part that calculates angle at the first row.

       do 15 j=1,n
       if(a(m,j).EQ.(0.,0.))then
       akx2=0.
       else
       akx2=dimag(cdlog(a(m,j)))
       endif
       if(a(m-1,j).EQ.(0.,0.))then
       akx1=0.
       else
       akx1=dimag(cdlog(a(m-1,j)))
       endif
       if(abs(akx2-akx1).GT.pi)then
       akx=sign((2.*pi-(abs(akx1)+abs(akx2)))/dx,akx1)
       else
       akx=(akx2-akx1)/dx
       endif
       if(j.NE.n)then
       if(a(m,j+1).EQ.(0.,0.))then
       aky2=0.
       else
       aky2=dimag(cdlog(a(m,j+1)))
       endif
       if(a(m,j).EQ.(0.,0.))then
       aky1=0.
       else
       aky1=dimag(cdlog(a(m,j)))
       endif
       else
       if(a(m,j).EQ.(0.,0.))then
       aky2=0.
       else
       aky2=dimag(cdlog(a(m,j)))
       endif
       if(a(m,j-1).EQ.(0.,0.))then
       aky1=0.
       else
       aky1=dimag(cdlog(a(m,j-1)))
       endif
       endif
       if(abs(aky2-aky1).GT.pi)then
       aky=sign((2.*pi-(abs(aky1)+abs(aky2)))/dy,aky1)
       else
       aky=(aky2-aky1)/dy
       endif
       thet(j)=datan2(aky,(akx+kb(m)))
       thet(j)=180.*thet(j)/pi
15    continue
       write(7,202)(thet(j),j=1,n,nd)
!C END NEW

!c-- pass angle on first row -- Fengyan 02/04/2002.

        do j=1,n,nd
        jj=(j-1)/nd+1
        Pass_Theta(1,jj)=thet(j)
        enddo


        endif
      
!c Store first row of wave heights on unit   CHECK UNIT 12.

      if(fname6.ne.' ')then
      write(26,202)(2*cdabs(a(1,j))/dconv(iu),j=1,n,nd)
      endif

!c--- Pass wave height (mks unit) -- Fengyan 02/04/2002.

        do j=1,n,nd
        jj=(j-1)/nd+1
        Pass_height(1,jj)=2d0*cdabs(a(1,j))
        Pass_ibrk(1,jj)=0
        enddo

!c Store surface on file |fname9|.

      if(fname9.ne.' ') then
      x(1)=0

      write(9,*) n
      write(9,*) (y(j),j=1,n)
      write(9,*) x(1)
      write(9,*) (dble(a(1,j)),j=1,n)

      endif
    
!c Now execute model for the |ifreq| frequency over each of |mr|
!c grid blocks.  |ir| is the controlling index value.
      
      do 100 ir=1,(mr-1)
      
!c Establish interpolated grid block for segment |ir|.
      
      call grid(ifreq,ir)
  
!c If |ir=1| write initial values on |iun(3)|.
  
      if(ir.eq.1) then

      write(5,201) x(1)/dconv(iu),psibar

      endif
      
!c Calculate constants for each grid block.
      
      call con(ifreq,ir)
     
!c Perform finite difference calculations.
      
      call fdcalc(ifreq,ir)
      
!c Grid block |ir| done, print output and go to next grid.


 100  continue 

      if(ioutput.eq.2)then
      write(20,*)(a(m,j),j=1,n)
      endif

!c Termination for the |surface.dat| file.

      if(fname9.ne.' ') then

        x(1)=-100.
        write(9,*) x(1)

      endif

!c  Calculate radiation stresses using new formula (roller and spline
!c  transition) -- fyshi 02/04/2002

        Pass_period=2.*pi/freqs(1,itime)

      call calculate_wave_forcing
      
!c Model complete for the |ifreq| frequency component, go to the next frequency 
!c component.
      
 200  continue
      
!c  Runs completed for all frequencies.  Return to end of main program.
      
      return

 201  format(' x=',f10.2,'  psibar=',f20.4)
 202  format(501(f10.4))
 203  format('1',20x,'model execution, frequency',                      &
     &' component',i4,/,/)

      end subroutine


      subroutine grid(ifreq,ir)
    
!      include 'param.h'

!      include 'common.h'

!      include 'pass.h'
#include "param.h"
#include "common.h"
#include "pass.h"

      real*8 dref
      integer i,j

!c Constants.

      pi=3.1415927
      eps=1.0e-05

!c Perform $y$-interpolation on reference grid.

!c Interpolate first row.

      do 10 j=1,n,nd
      d(1,j)=dr(ir,((j-1)/nd+1))!+Intp_eta_Wave(ir,((j-1)/nd+1))
      u(1,j)=ur(ir,((j-1)/nd+1))
      v(1,j)=vr(ir,((j-1)/nd+1))
 10   continue

      do 12 jj=2,nr
      do 11 j=1,(nd-1)
      jjj=nd*(jj-2)+(j+1)
      d(1,jjj)=(dr(ir,jj)-dr(ir,jj-1))*y(jjj)/dyr                       &
     &+(yr(jj)*dr(ir,jj-1)-yr(jj-1)*dr(ir,jj))/dyr
      u(1,jjj)=(ur(ir,jj)-ur(ir,jj-1))*y(jjj)/dyr                       &
     &+(yr(jj)*ur(ir,jj-1)-yr(jj-1)*ur(ir,jj))/dyr
      v(1,jjj)=(vr(ir,jj)-vr(ir,jj-1))*y(jjj)/dyr                       &
     &+(yr(jj)*vr(ir,jj-1)-yr(jj-1)*vr(ir,jj))/dyr
 11   continue
 12   continue
    
!c Set number of $x$ points and define $x$ values.
    
      if(ispace.eq.0) then
      
!c |ispace=0|, program sets subdivisions.
  
      do 13 j=1,n
      dref=d(1,j)+tide(ifreq)
      if(dref.lt.0.001) dref=0.001
      call wvnum(dref,u(1,j),freqs(ifreq,itime),k(1,j),eps,icdw,1,j)
 13   continue

      npts=0
      sumk=0.
      do 14 j=1,n
      if(d(1,j).gt.0.05) then
        sumk=sumk+k(1,j)
        npts=npts+1
      endif
 14   continue

      kb(1)=sumk/dfloat(npts)
      alw=2.*pi/kb(1)
      anw=dxr/alw
      np=ifix(5.*anw)
      if(np.lt.1) np=1
      md(ir)=min((ix-1),np)
      if(np.gt.(ix-1)) write(5,100) ir

      endif
    
!c    |ispace=1|, user specified subdivision.
    
      m=md(ir)+1
      dx=dxr/dfloat(md(ir))
      do 15 i=1,m
      x(i)=xr(ir)+dfloat(i-1)*dx
 15   continue

!c    interpolate values on |m| row.
    
      do 16 j=1,n,nd
      d(m,j)=dr(ir+1,((j-1)/nd+1))!+Intp_eta_Wave(ir+1,((j-1)/nd+1))
      u(m,j)=ur(ir+1,((j-1)/nd+1))
      v(m,j)=vr(ir+1,((j-1)/nd+1))
 16   continue



      do 18 jj=2,nr
      do 17 j=1,(nd-1)
      jjj=nd*(jj-2)+(j+1)
      d(m,jjj)=(dr(ir+1,jj)-dr(ir+1,jj-1))*y(jjj)/dyr                   &
     &+(yr(jj)*dr(ir+1,jj-1)-yr(jj-1)*dr(ir+1,jj))/dyr
      u(m,jjj)=(ur(ir+1,jj)-ur(ir+1,jj-1))*y(jjj)/dyr                   &
     &+(yr(jj)*ur(ir+1,jj-1)-yr(jj-1)*ur(ir+1,jj))/dyr
      v(m,jjj)=(vr(ir+1,jj)-vr(ir+1,jj-1))*y(jjj)/dyr                   &
     &+(yr(jj)*vr(ir+1,jj-1)-yr(jj-1)*vr(ir+1,jj))/dyr
 17   continue
 18   continue
    
!c Interpolate values in |x|-direction.
    
      do 19 i=2,m-1
      do 19 j=1,n
      d(i,j)=(d(m,j)-d(1,j))*x(i)/dxr+(x(m)*d(1,j)-x(1)*d(m,j))/dxr
      u(i,j)=(u(m,j)-u(1,j))*x(i)/dxr+(x(m)*u(1,j)-x(1)*u(m,j))/dxr
      v(i,j)=(v(m,j)-v(1,j))*x(i)/dxr+(x(m)*v(1,j)-x(1)*v(m,j))/dxr
 19   continue
    
!c Add in user specified grid subdivisions (read from unit 3).
    
      do 30 jr=1,nr-1
      if(isd(ir,jr).eq.1) then
      js=nd*jr+(1-nd)
      jf=js+nd
      read(3,101)((d(i,j),j=js,jf),i=1,m)

      if(icur.eq.1)then
      read(3,101)((u(i,j),j=js,jf),i=1,m)
      read(3,101)((v(i,j),j=js,jf),i=1,m)
      endif

      do 31 i=1,m
      do 31 j=js,jf
      d(i,j)=d(i,j)*dconv(iu)
      u(i,j)=u(i,j)*dconv(iu)
      v(i,j)=v(i,j)*dconv(iu)
 31   continue

      end if

 30   continue
    
!c Add tidal offset to all rows and establish thin film.  Set
!c current speed to zero in thin film area.
    
      do 20 i=1,m
      do 20 j=1,n

      d(i,j)=d(i,j)+tide(ifreq)

      if(d(i,j).lt.0.001) then
      d(i,j)=0.001
      u(i,j)=0.0
      v(i,j)=0.0
      endif

 20   continue
    
!c Interpolation complete, return to |model|.
    
      return

 100  format(' model tried to put more spaces than allowed in',         &
     &' grid block ',i3)
 101  format(501f10.4)

      end subroutine


      subroutine con(ifreq,ir)
    
!      include 'param.h'

!      include 'common.h'
#include "param.h"
#include "common.h"


!c Constants.

      eps=1.0e-05
      g=9.80621
    
!c Calculate constants.
    
      do 1 i=1,m
      do 1 j=1,n
      call wvnum(d(i,j),u(i,j),freqs(ifreq,itime),k(i,j),eps,icdw,i,j)
      sig(i,j)=freqs(ifreq,itime)-k(i,j)*u(i,j)
      akd=k(i,j)*d(i,j)
      q(i,j)=(1.+akd/(sinh(akd)*cosh(akd)))/2.
      p(i,j)=q(i,j)*g*tanh(akd)/k(i,j)
      dd(i,j)=(cosh(4.*akd)+8.-2.*(tanh(akd)**2))/(8.*(sinh(akd)**4.))
      bottomu(i,j)=g*k(i,j)/(2*freqs(ifreq,itime)*cosh(akd))
 1    continue

!c We seem to have occasional problems in adverse currents when the water
!c is shallow, probably due to using artificial current fields.  We now
!c check to see if the total advection velocity is negative, make it a small
!c positive number if not, and flag the change to the user.

      if(icur.eq.1) then

      do 2 i=1,m
      do 2 j=1,n

      cgroup=sqrt(p(i,j)*q(i,j))
      advect=cgroup+u(i,j)

      if(advect.le.0.) then
      advect=0.1
      utemp=advect-cgroup
      write(5,100)ir,i,j,u(i,j),utemp
      u(i,j)=utemp
      endif

 2    continue

      endif
  
!c Calculate the dissipation term |w|.
  
      call diss
  
!c Calculate the mean |kb| on each row.
  
      do 11 i=1,m
      npts=0
      sumk=0.
      do 10 j=1,n
      if(d(i,j).gt.0.05) then
        sumk=sumk+k(i,j)
        npts=npts+1
      endif
 10   continue
      if(npts.eq.0)then
      kb(i)=k(i,1)
      else
      kb(i)=sumk/dfloat(npts)
      endif
 11   continue

      return

 100  format('at grid block ',i3,', location i=',i4,', j=',i4,' model   &
     &found a negative, blocking velocity u=',f8.2,' and reduced it to  & 
     &u=', f8.2)

      end subroutine

  
      subroutine fdcalc(ifreq,ir)

!      include 'param.h'

!      include 'common.h'

!      include 'pass.h'
#include "param.h"
#include "common.h"
#include "pass.h"

      integer i,j

      real*8 kap,ksth1,ksth2,akx,aky,akx1,aky1
      complex*16 c1,c2,c3,cp1,cp2,cp3,ci,damp
      complex*16 ac(iy),bc(iy),cc(iy),rhs(iy),sol(iy)
      real*8 thet(iy), urs(iy),height(iy),heightb(iy)
      real*8 sxy(iy),sxx(iy),syy(iy),sbxx(iy),sbxy(iy),sbyy(iy)
      real*8 cg, pv, bet, dv, deltap, f1    !jcw

      cg(i,j)=sqrt(p(i,j)*q(i,j))

      pv(i,j)=p(i,j)-v(i,j)*v(i,j)

      bet(i,j)=-4.*(k(i+1,j) - k(i,j))/(dx*((k(i+1,j) + k(i,j))**2))    &
     &  -4.*(k(i+1,j)*(p(i+1,j) -u(i+1,j)**2) -                         &
     &  k(i,j)*(p(i,j)-u(i,j)**2)) /                                    &
     &  (dx*((k(i+1,j) +  k(i,j))**2.)*(p(i+1,j) + p(i,j) -             &
     &  (u(i+1,j)**2 +                                                  &
     &  u(i,j)**2) ) )
  
      dv(i,j)=( cg(i+1,j) + u(i+1,j) )/sig(i+1,j) -                     &
     &  (cg(i,j)+u(i,j))/sig(i,j)                                       &
     &  -delta1*dx*( (v(i+1,j+1)/sig(i+1,j+1)) + (v(i,j+1)/sig(i,j+1) ) &
     &  -(v(i+1,j-1)/sig(i+1,j-1))-(v(i,j-1)/sig(i,j-1)) )/(2.*dy)

      damp(i,j)=2.*ci*cdamp*( (cg(i+1,j) + u(i+1,j)) + (cg(i,j) +       &
     &  u(i,j)) )/                                                      &
     & ( dy*dy*( k(i+1,j)**2 + k(i,j)**2 ) )

      deltap(i,j)=a1-b1*kb(i)/k(i,j)

      cp1(i,j)=(cg(i+1,j)+u(i+1,j))*dcmplx(1.,dx*(kb(i+1)-a0*k(i+1,j))) &
     &+dcmplx(1.,0.)*(cg(i,j)+u(i,j)+dv(i,j)*(sig(i+1,j)+sig(i,j))/4.)  &
     &+2.*omeg*dcmplx(0.,1.)*(-b1)*bet(i,j)*(u(i+1,j)+u(i,j))/sig(i+1,j)& 
     &+4.*omeg*(-b1)*dcmplx(0.,1.)*(3.*(u(i+1,j)-u(i,j))/dx+(v(i+1,j+1) &
     &+v(i,j+1)-v(i+1,j-1)-v(i,j-1))/(4.*dy))/(sig(i+1,j)*(k(i+1,j)+    &
     & k(i,j)))                                                         &
     &+dcmplx(-2.*(-b1)/(dy*dy*(k(i+1,j)+k(i,j)))+b1*bet(i,j)*dx/(2.*   &
     & dy*dy),                                                          &
     &-deltap(i,j)*dx/(2.*dy*dy))*(pv(i+1,j+1)+2.*pv(i+1,j)+            &
     &pv(i+1,j-1))/sig(i+1,j)-dcmplx(1.,0.)*omeg*delta2*(3.*u(i+1,j)    &
     &+u(i,j))/(2.*sig(i+1,j))+ci*omeg*(a0-1.)*k(i+1,j)*u(i+1,j)*       &
     &  dx/sig(i+1,j)                                                   &
     &+2.*ifilt*damp(i,j) +dcmplx(1.,0.)*alphn*dx

      cp2(i,j)=dcmplx( (-delta1*dx)*(v(i+1,j) + v(i,j))/(2.*dy)         &
     &+b1*u2*bet(i,j)*(u(i+1,j)*v(i+1,j) + u(i,j)*v(i,j))/              &
     &  (dy*sig(i+1,j+1)),                                              &
     &(-delta1*u2)*(u(i+1,j+1)*v(i+1,j+1) + u(i,j+1)*v(i,j+1)+          &
     &2.*u(i+1,j)*v(i+1,j))/(2.*dy*sig(i+1,j+1))+dx*(-b1)*bet(i,j)      &
     &*(sig(i+1,j)*v(i+1,j)+sig(i,j)*v(i,j))/(2.*dy*sig(i+1,j+1)))      &
     &+dcmplx(2.*(-b1)/(dy*dy*(k(i+1,j)+k(i,j)))+                       &
     &  (-b1)*bet(i,j)*dx/(2.*dy*dy),                                   &
     &+deltap(i,j)*dx/(2.*dy*dy))*(pv(i+1,j+1) +                        &
     &  pv(i+1,j))/sig(i+1,j+1) +                                       &
     &4.*dcmplx(0.,1.)*(-b1)*sig(i+1,j)*v(i+1,j)/(dy*sig(i+1,j+1)*      &
     &  (k(i+1,j) +                                                     &
     &k(i,j))) -ifilt*damp(i,j)

      cp3(i,j)=dcmplx(-(-delta1*dx)*(v(i+1,j)+v(i,j))/(2.*dy)           &
     &+(-b1)*u2*bet(i,j)*(u(i+1,j)*v(i+1,j)+u(i,j)*v(i,j))/             &
     &  (dy*sig(i+1,j-1)),                                              &
     & -(-delta1*u2)*(u(i+1,j-1)*v(i+1,j-1)+u(i,j-1)*v(i,j-1)+          &
     &2.*u(i+1,j)*v(i+1,j))/(2.*dy*sig(i+1,j-1))-                       &
     &dx*(-b1)*bet(i,j)*(sig(i+1,j)*v(i+1,j) + sig(i,j) * v(i,j)) /     &
     &(2.*dy*sig(i+1,j-1))) +dcmplx(2.*(-b1)/(dy*dy*(k(i+1,j) + k(i,j)))&
     & +                                                                &
     &(-b1)*bet(i,j)*dx/(2.*dy*dy), -(-deltap(i,j)*dx) / (2.*dy*dy))*   &
     &  (pv(i+1,j)                                                      &
     &+ pv(i+1,j-1))/sig(i+1,j-1)- 4.*dcmplx(0.,1.)*(-b1)*sig(i+1,j)*   &
     &  v(i+1,j) /                                                      &
     &(dy*sig(i+1,j-1)*(k(i+1,j)+k(i,j))) -ifilt*damp(i,j)

      c1(i,j)=dcmplx(cg(i+1,j)+u(i+1,j)-dv(i,j)*(sig(i+1,j) +           &
     &  sig(i,j))/4.,0.) +                                              &
     &dcmplx(1.,-dx*(kb(i)-a0*k(i,j))) * (cg(i,j)+u(i,j)) +             &
     &2.*dcmplx(0.,1.)*omeg*(-b1)*bet(i,j)*(u(i+1,j) + u(i,j))/sig(i,j)+&
     &4.*dcmplx(0.,1.)*omeg*(-b1)*(3.*(u(i+1,j) - u(i,j))/dx+           &
     &  (v(i+1,j+1)+                                                    &
     &v(i,j+1) - v(i+1,j-1) - v(i,j-1))/(4.*dy))/(sig(i,j)*(k(i+1,j) +  &
     &  k(i,j)))                                                        &
     &+ dcmplx(2.*b1/(dy*dy*(k(i+1,j)+k(i,j))) -b1 * bet(i,j) *         &
     &dx/(2.*dy*dy), +deltap(i,j)*dx/(2.*dy*dy))*(pv(i,j+1) +           &
     &2.*pv(i,j)+pv(i,j-1))/sig(i,j) - dcmplx(1.,0.)*omeg*delta2*       &
     &(3.*u(i+1,j)+u(i,j))/(2.*sig(i,j)) - ci*omeg* (a0-1.) * k(i,j) *  &
     &u(i,j)*dx/sig(i,j) + 2.*ifilt*damp(i,j) - dcmplx(1.,0.)*alphn*dx

      c2(i,j)=dcmplx(delta1*dx*( v(i+1,j) + v(i,j) )/(2.*dy)            &
     & + b1*u2*bet(i,j)*(u(i+1,j)*v(i+1,j) + u(i,j)*v(i,j))/            &
     &  (dy*sig(i,j+1)),                                                &
     &(-delta1*u2)*(u(i+1,j+1)*v(i+1,j+1) + u(i,j+1)*v(i,j+1)           &
     &+2.*u(i,j)*v(i,j))/(2.*dy*sig(i,j+1)) + 4.*(-b1)*sig(i,j)*        &
     &  v(i,j) /                                                        &
     &(dy*(k(i+1,j) +k(i,j))*sig(i,j+1)) - dx*(-b1)*bet(i,j)*(          &
     &sig(i+1,j)*v(i+1,j) + sig(i,j)*v(i,j)) / ( 2.*dy*sig(i,j+1)  )) + &
     &dcmplx(2.*(-b1)/(dy*dy*( k(i+1,j) + k(i,j) )) + b1*bet(i,j)*dx/   &
     &  (2.*dy*dy),                                                     &
     & (-deltap(i,j)*dx)/(2.*dy*dy))*(pv(i,j+1) + pv(i,j))/sig(i,j+1)   &
     &  -                                                               &
     &ifilt*damp(i,j)

      c3(i,j)=dcmplx((-delta1*dx)*(v(i+1,j)+v(i,j))/(2.*dy) -           &
     &  b1*u2*bet(i,j) *                                                &
     &(u(i+1,j)*v(i+1,j) + u(i,j)*v(i,j))/(dy*sig(i,j-1)), (delta1*u2)* &
     &(u(i+1,j-1) * v(i+1,j-1) + u(i,j-1)*v(i,j-1) +2.*u(i,j)*v(i,j))/( &
     &2.*dy*sig(i,j-1) ) - 4.*(-b1)*sig(i,j)*v(i,j)/(dy*(k(i+1,j) +     &
     &k(i,j))*sig(i,j-1)) + dx*(-b1)*bet(i,j)*(sig(i+1,j)*v(i+1,j) +    &
     &sig(i,j)*v(i,j))/(2.*dy*sig(i,j-1))) + dcmplx(-2.*b1/(dy*dy*(     &
     &  k(i+1,j) +                                                      &
     &k(i,j) )) + b1*bet(i,j)*dx/(2.*dy*dy), (-deltap(i,j)*dx) /        &
     &  (2.*dy*dy)) *                                                   &
     & ( pv(i,j) + pv(i,j-1) )/sig(i,j-1)-ifilt*damp(i,j)

      f1(i,j)=tanh(k(i,j)*d(i,j))**5.

      f2(i,j)=(k(i,j)*d(i,j)/sinh(k(i,j)*d(i,j)))**4.



!c Constants defining the paraboli!cmodel angular aperture. 
     
!c 70 degree minimax coefficients.
!c   
!c     a0=0.994733 
!c     a1=-0.890065  
!c     b1=-0.451641

!c Pad\'{e} coefficients.

      a0=1.0 
      a1=-0.75 
      b1=-0.25 

!c Additional constants.

      u2=1.0
      kap=0.78
      gam=0.4
      omeg=freqs(ifreq,itime)
      pi=3.1415927
      ci=dcmplx(0.,1.)

!c cdamp is the constant for the noise suppression formulation.  It can be
!c changed.

      cdamp=0.00025

      alphn=0.
      delta1=a1-b1
      delta2=1+2.*a1-2.*b1

!c Initialize breaking index if |ir = 1|.

      if(ir.eq.1) then
      ifilt=0
        do 100 j=1,n
        ibr(j)=0
        wb(1,j)=dcmplx(0.,0.)
 100    continue
      endif

!c Solution for |m| grid blocks in reference block |ir|.

      do 200 i=1,(m-1)
      ip1=i+1
      it=1
      ih=1

!c r.h.s. of matrix equation. 

      rhs(1)=dcmplx(0.,0.)

      do 1 j=2,(n-1)
      rhs(j)=c1(i,j)*a(i,j)+c2(i,j)*a(i,j+1)+c3(i,j)*a(i,j-1)           &
     &-dx*(w(i,j)+wb(1,j))*a(i,j)/2.                                    &
     &-dx*dcmplx(0.,1.)*an*anl*sig(i,j)*k(i,j)*k(i,j)*dd(i,j)*          &
     &(1.-dfloat(ibr(j)))*(cdabs(a(i,j))**2.)*a(i,j)/2.                 &
     &-dx*dcmplx(0.,1.)*(1.-dfloat(ibr(j)))*an*(1.-anl)                 &
     &*sig(i,j)*((1.+f1(i,j)*k(i,j)*k(i,j)*(cdabs(a(i,j))**2.)          &
     &*dd(i,j))*tanh(k(i,j)*d(i,j)+f2(i,j)*k(i,j)*cdabs(a(i,j)))        &
     &/tanh(k(i,j)*d(i,j))-1.)*a(i,j)/2.
 1    continue

      rhs(n)=dcmplx(0.,0.)

!c Return here for iterations. 
  
 2    if(it.eq.1)ii=i
      if(it.eq.2)ii=ip1

!c Establish boundary conditions.

      if(ibc.eq.1)then
      ksth1=dble((2.*(a(i,2)-a(i,1))/((a(i,2)+a(i,1))*dy))*dcmplx       &
     &(0.,-1.))
      ksth2=dble((2.*(a(i,n)-a(i,n-1))/((a(i,n)+a(i,n-1))*dy))*         &
     &dcmplx(0.,-1.))
      bc(1)=dcmplx(1.,ksth1*dy/2.)
      cc(1)=-dcmplx(1.,-ksth1*dy/2.)
      bc(n)=-dcmplx(1.,-ksth2*dy/2.)
      ac(n)=dcmplx(1.,ksth2*dy/2.)
      else
      bc(1)=dcmplx(1.,0.)
      cc(1)=-bc(1)
      bc(n)=dcmplx(1.,0.)
      ac(n)=-bc(n)
      endif

!c Calculate dissipation in rows where breaking occurs. 

      do 3 j=1,n

      if(ibr(j).eq.1) wb(2,j)=dcmplx(1.,0.)*0.15*cg(ip1,j)*             &
     &(1.-(gam*d(ip1,j)/(2.*cdabs(a(ii,j))))**2.)/d(ip1,j)

      if(ibr(j).eq.0) wb(2,j)=dcmplx(0.,0.)

 3    continue

!c Coefficients for forward row. 

      do 4 j=2,(n-1)

      ac(j)=cp3(i,j)

      bc(j)=cp1(i,j)+(dx/2.)*(w(i+1,j)+wb(2,j))+dcmplx(0.,an*anl)*      &
     &  sig(i+1,j)                                                      &
     &*k(i+1,j)*k(i+1,j)*dd(i+1,j)*(cdabs(a(ii,j))**2.)*(dx/2.)         &
     &+dcmplx(0.,an*(1.-anl))*sig(i+1,j)*(dx/2.)*((1.+f1(i+1,j)*        &
     &k(i+1,j)*k(i+1,j)*(cdabs(a(ii,j))**2.)*dd(i+1,j))*tanh(k(i+1,j)*  &
     &d(i+1,j)+f2(i+1,j)*k(i+1,j)*cdabs(a(ii,j)))/tanh(k(i+1,j)*d(i+1,j &
     &))-1.)

      cc(j)=cp2(i,j)

 4    continue

!c Update solution one step. 

      call ctrida(1,n,ac,bc,cc,rhs,sol)
      do 5 j=1,n
      a(ip1,j)=sol(j)
      sol(j)=dcmplx(0.,0.)
 5    continue
      if((it.eq.2).or.(ih.eq.2)) go to 8

!c     write(*,*)sol(6)

!c Check for start or stop of breaking in each row. 

      do 6 j=1,n
      urs(j)=2.*cdabs(a(ip1,j))
 6    continue

      isave1=0
      isave2=0
      do 7 j=1,n
      iset=0
      ireset=0

      if(((urs(j)/d(ip1,j)).gt.kap).and.(ibr(j).eq.0)) iset=1

      if(iset.eq.1) then
        ibr(j)=1
        isave1=1
      end if

      if(((urs(j)/d(ip1,j)).lt.gam).and.(ibr(j).eq.1)) ireset=1

      if(ireset.eq.1) then
        ibr(j)=0
        isave2=1
      end if

 7    continue

      ih=2

!c Redo initial calculation if breaking status changes. 

      if((isave1.eq.1).or.(isave2.eq.1)) go to 2
 8    continue
      if(it.eq.2) go to 9
      it=2
      go to 2

 9    continue 
 
!c For Stokes model alone (|ntype.eq.2|), test to see whether Ursell 
!c parameter is too large. 

      if(ntype.eq.2) then
      do 11 j=1,n
      urs(j)=(cdabs(a(ip1,j))/d(ip1,j))/((k(ip1,j)*d(ip1,j))**2)
      if(urs(j).gt.0.5) write(5,204) urs(j),i,j
 11   continue
      end if

!c Roll back breaking dissipation coefficient at each row.

      do 12 j=1,n
      wb(1,j)=wb(2,j)
 12   continue

!c Calculate reference phase function for surface plotting. 

      psibar=psibar+(kb(ip1)+kb(i))*dx/2.

!c Store plotted surface if requested.

      if(fname9.ne.' ') then

      write(9,*) x(ip1)
      write(9,*)(dble(a(ip1,j)*cdexp(dcmplx(0.,psibar))),j=1,n)

      endif

!c Start filter if breaking is occuring. 

      do 13 j=1,n
      if(ibr(j).eq.1)ifilt=1
 13   continue


 200  continue

!cCalculate wave angles at reference grid rows.   Note:  angles are not well 
!cdefined in a directional, multicomponent sea, or where waves become short 
!ccrested.  This routine was heavily modified by Raul Medina, University of 
!cCantabria. It is further modified by Arun Chawla to take out a one-sided
!cderivative that introduced an asymmetry bias.

      do 15 j=1,n

      if(a(m,j).eq.(0.,0.)) then
          akx2=0.
      else
          akx2=dimag(cdlog(a(m,j)))
      endif

      if(a(m-1,j).eq.(0.,0.)) then
          akx1=0.
      else
          akx1=dimag(cdlog(a(m-1,j)))
      endif

      if(abs(akx2-akx1).gt.pi)then
      akx=sign((2.*pi-(abs(akx1)+abs(akx2)))/dx,akx1)
      else
      akx=(akx2-akx1)/dx
      endif

      if(j.eq.1) then
          if(a(m,j+1).eq.(0.,0.)) then
              aky2=0.
          else
              aky2=dimag(cdlog(a(m,j+1)))
          endif

          if(a(m,j).eq.(0.,0.)) then
              aky1=0.
          else
              aky1=dimag(cdlog(a(m,j)))
          endif

          if(abs(aky2-aky1).gt.pi)then
              aky=sign((2.*pi-(abs(aky1)+abs(aky2)))/dy,aky1)
          else
              aky=(aky2-aky1)/dy
          endif

      endif

      if(j.eq.n) then

          if(a(m,j).eq.(0.,0.)) then
              aky2=0.
          else
              aky2=dimag(cdlog(a(m,j)))
          endif

          if(a(m,j-1).eq.(0.,0.))then
              aky1=0.
          else
              aky1=dimag(cdlog(a(m,j-1)))
          endif

          if(abs(aky2-aky1).gt.pi)then
              aky=sign((2.*pi-(abs(aky1)+abs(aky2)))/dy,aky1)
          else
              aky=(aky2-aky1)/dy
          endif

      endif

      if((j.gt.1).and.(j.lt.n)) then

          if(a(m,j+1).eq.(0.,0.)) then
              aky2=0.
          else
              aky2=dimag(cdlog(a(m,j+1)))
          endif

          if(a(m,j-1).eq.(0.,0.))then
              aky1=0.
          else
              aky1=dimag(cdlog(a(m,j-1)))
          endif

          if(abs(aky2-aky1).gt.pi)then
              aky=sign((2.*pi-(abs(aky1)+abs(aky2)))/(2.*dy),aky1)
          else
              aky=(aky2-aky1)/(2.*dy)
          endif

      endif

      thet(j)=datan2(aky,(akx+kb(m)))

 15   continue

!c Estimation of Wave Height

      do 17 j=1,n

      height(j)=2.*cdabs(a(m,j))

   17 continue
!c     write(*,*)a(m,6)

!c Spatial smoothing of wave height if |(nd.gt.1)|.

      if (nd.NE.1) then

      jh = int(nd/2) 

      do j = 1,n,nd

         heightb(j) = 0.

         if (j.eq.1) then

            do jj = 1,1+jh
                heightb(1) = heightb(1) + height(jj)**2.
            end do

            heightb(1) = sqrt( heightb(1)/(jh+1) )
 
         endif

         if (j.eq.n) then

            do jj = n-jh,n              
                heightb(n) = heightb(n) + height(jj)**2.
            end do

            heightb(n) = sqrt( heightb(n)/(jh+1) )

         endif

         if((j.gt.1).and.(j.lt.n)) then

            do jj = j-jh,j+jh
                heightb(j) = heightb(j) + height(jj)**2.
            end do

            heightb(j) = sqrt( heightb(j)/(2*jh+1) )

         endif

      end do  

      endif

      if(nd.eq.1) then

         do j=1,n
         heightb(j)=height(j)
         end do

      endif

!c If |ismooth|=1, perform a much more drastic smoothing over a wider
!c footprint.

      if(ismooth.eq.1) then

      do j=1,n,nd
      height(j)=heightb(j)
      end do

!c First, employ a simple average over a wide base to get rid of most noise
!c and most medium scale variability.

      do j=1+nrs*nd,n-nrs*nd,nd
      heightb(j)=0.

      do is=-nrs,nrs,1
      heightb(j)=heightb(j)+height(j+is*nd)**2.
      end do

      heightb(j)=dsqrt(heightb(j)/dfloat(2*nrs+1) )
      end do

!c Then, emply a weighted three point average to get rid of small scale
!c noise.

      do j=1,n,nd
      height(j)=heightb(j)
      end do

      do j=1+nd,n-nd,nd
      heightb(j)=sqrt(0.25*height(j-nd)**2. + 0.5*height(j)**2.         &
     & + 0.25*height(j+nd)**2.)
      end do

      endif


!c Print out |abs(a)| at grid reference points. 

      mm1=m-1
      write(5,205) (ir+1),mm1
      write(5,202) x(m)/dconv(iu), psibar

!c Wave heights on |height.dat|.

      write(26,203) (heightb(j)/dconv(iu),j=1,n,nd)

! --- Pass wave height -- mks unit Fengyan 02/04/2002
        do j=1,n,nd
        jj=(j-1)/nd+1
        Pass_Height(ir+1,jj)=heightb(j)
        enddo     


! --- Pass wave number

        do j=1,n,nd
        jj=(j-1)/nd+1
        Pass_WaveNum(ir+1,jj)=k(m,j)
        enddo     

!c Wave angles on |angle.dat|.

      do 18 j=1,n
      thet(j)=180.*thet(j)/pi
 18   continue      

      write(7,203) (thet(j),j=1,n,nd)

! --- Pass angle -- Fengyan 02/04/2002
        do j=1,n,nd
        jj=(j-1)/nd+1
        Pass_Theta(ir+1,jj)=thet(j)
        enddo


!c Water depths on |depth.dat|.

      write(8,203) (d(m,j)/dconv(iu),j=1,n,nd)

! --- Pass wave cg, c and ibrk -- Fengyan 02/04/2002
        do j=1,n,nd
        jj=(j-1)/nd+1
        Pass_Cg(ir+1,jj)=cg(m,j)
        Pass_C(ir+1,jj)=sig(m,j)/k(m,j)
        Pass_ibrk(ir+1,jj)=ibr(j)
        enddo

!c Bottom velocities on |bottomu.dat|.

      if(fname7.ne.' ') then
      do 19 j=1,n,nd
      bottomu(m,j)=bottomu(m,j)*cdabs(a(m,j))
 19   continue
      write(17,203) (bottomu(m,j)/dconv(iu),j=1,n,nd)
      endif

!c Write out reference grid data on disk file |iun(3)|. 

!c     write(iun(3),*) x(m)/dconv(iu),psibar
!c     write(iun(3),*)(a(m,j)/dconv(iu),j=1,n,nd)

!c Roll back solution to first grid level.
 
      do 201 j=1,n
      a(1,j)=a(m,j)
 201  continue

      return

 202  format(' x=',f10.2,'    reference phase psibar=',f20.4)
 203  format(501(f10.4))
 204  format(' ',/,/,' warning: Ursell number =',f10.4,' encountered at'&
     &,'grid location',i6,',',i6,/,                                     &
     &' should be using Stokes-Hedges model (ntype=1) due to shallow',  &
     &'water')
 205  format(' grid row ir=',i3,',  ',i3,' x-direction subdivisions',   &
     &' used')

      end subroutine


  
      subroutine ctrida(ii,l,a,b,c,d,v)

!      include 'param.h'
#include "param.h"


      complex*16 a(iy),b(iy),c(iy),d(iy),v(iy),beta(iy),gamma(iy)
  
!c Compute intermediate vectors |beta| and |gamma|.
  
      beta(ii)=b(ii)
      gamma(ii)=d(ii)/beta(ii)
      iip1=ii+1

      do 1 i=iip1,l
      beta(i)=b(i)-a(i)*c(i-1)/beta(i-1)
      gamma(i)=(d(i)-a(i)*gamma(i-1))/beta(i)
 1    continue
  
!c Compute solution vector |v|.
  
      v(l)=gamma(l)
      last=l-ii

      do 2 k=1,last
      i=l-k
      v(i)=gamma(i)-c(i)*v(i+1)/beta(i)
 2    continue

      return

      end subroutine
    
    
      subroutine diss

!      include 'param.h'

!      include 'common.h'
#include "param.h"
#include "common.h"

      real*8 nu,cp,kd

!c Statement function.

      sq(i,j)=sqrt(nu/(2.*sig(i,j)))

!c Constants.

      nu=1.3e-06
      cp=4.5e-11
      g=9.80621
      pi=3.1415927
    
!c    Value of |f| here is value assuming $\tau=(f/8)u^{2}$.

!c    $f=4 f_w$ ; $f_w$ is the wave friction factor

      f=0.01*4.0
    
      do 1 j=1,n
      do 1 i=1,m
      w(i,j)=dcmplx(0.,0.)
      kd=k(i,j)*d(i,j)
  
!c    If |iff(1) = 1|, use turbulent boundary layer damping.
  
      if(iff(1).eq.1) w(i,j)=2.*f*cdabs(a(1,j))*sig(i,j)*k(i,j)/        &
     &(sinh(2.*kd)*sinh(kd)*3.*pi)
    
!c If |iff(2) = 1|, add porous bottom damping.
    
      if(iff(2).eq.1) w(i,j)=w(i,j)+(g*k(i,j)*cp/(nu*(cosh(kd)**2)))    &
     & *dcmplx(1.,0.)
    
!c  If |iff(3) = 1|, add boundary layer damping.
    
      if(iff(3).eq.1) w(i,j)=w(i,j)+2.*k(i,j)*sig(i,j)*sq(i,j)          &
     &*(1.+(cosh(kd)**2))*dcmplx(1.,-1.)/sinh(2.*kd)
 1    continue

      return

      end  subroutine
  
      function rand1(x)
      real*8 x                        !jcw
      ix=int(32767.*x)                !jcw
!      ix=ifix(32767.*x)
      irand=mod(4573*ix+6923,32767)
      rand1=dfloat(irand)/32767.
      return
      end function rand1
    
      function rdfact(xi)
      real*8 xi, prod                       !jcw
      prod=1.
      if (xi.gt.1.)then
         do 17 ii=2,int(xi)
           prod=prod*dfloat(ii)
  17     continue
      endif

      rdfact=prod

      return
      end function rdfact

      subroutine bnum(in,n,bn)    
      real*8 xt, xin                        !jcw
!      xin=in
      xin=dfloat(in)                   !jcw
      xt=rdfact(xin)
!      xt=rdfact(dfloat(xin))          !jcw
      xb=rdfact(dfloat(n))*rdfact(dfloat(in-n))
      bn=xt/xb

      return

      end subroutine
    
    
      subroutine acalc(thmax,nsp,a)

      itn=2*nsp
      call bnum(itn,nsp,bn)
      a=thmax*bn/(2.**(itn-1))
      sum=0.
      do 10 ik=1,nsp
        ki=ik-1
        call bnum(itn,ki,bn)
  10    sum=sum+bn*sin(dfloat(nsp-ki)*thmax)/dfloat(nsp-ki)
      a=a+sum/(2.**(itn-2))

      return

      end subroutine


      subroutine wvnum(d,u,s,k,eps,icdw,i,j)  

!      include 'param.h'
#include "param.h"

      common/ref2/dr(ixr,iyr),ur(ixr,iyr),vr(ixr,iyr),iun(8),iinput,    &
     &      ioutput

      real*8 k,kn,dr,ur,vr,d,u,s

!c Constants.

      g=9.806
      pi=3.1415927
      k=s*s/(g*sqrt(tanh(s*s*d/g)))

!c Newton-Raphson iteration.

      do 1 ii=1,20
      f=s*s-2.*s*k*u+k*k*u*u-g*k*tanh(k*d)
      fp=-2.*s*u+2.*k*u*u-g*tanh(k*d)-g*k*d/(cosh(k*d)**2.)
      kn=k-f/fp
      if((abs(kn-k)/kn).lt.eps)go to 2
      k=kn
 1    continue
      t=2.*pi/(sqrt(g*k*tanh(k*d))+k*u)
      write(5,100)i,j,k,u,d,f,t
      icdw=1

      return

 2    k=kn

      return

 100  format(' wavenumber iter. failed to converge on row',i10,         &
     &'  column',i10,/,                                                 &
     &'      k=',f15.8,'      u=',f15.8,/,                              &
     &'      d=',f15.8,'      f=',f15.8,/,                              &
     &'      t=',f15.8)

      end subroutine

        subroutine calculate_wave_forcing

!        include 'param.h'
!        include 'common.h'
!        include 'pass.h'
#include "param.h"
#include "common.h"
#include "pass.h"

        integer i,j

        real*8 Sm,Sp                                                    &
     &  ,Sxxfake(Nx_Max),Sxyfake(Nx_Max),Syyfake(Nx_Max)                &
     &  ,Qwfake(Nx_Max),xrr(Nx_Max),y2(Nx_Max),Area(nx_max,ny_max),     &
     &   Afake(nx_max),qw(nx_max,ny_max),qwrol(nx_max)                  &
     &  ,smrol(nx_max),Afake2(nx_max),xx(nx_max)                        &
     &  ,xfake(Nx_Max),sxx(Nx_Max,Ny_Max),                              &
     &  sinth(nx_max,ny_max),costh(nx_max,ny_max),                      &
     &  sxy(Nx_Max,Ny_Max),dwf(nx_max,ny_max),                          &
     &  syy(Nx_Max,Ny_Max),qwtrs(Nx_Max,Ny_Max)                         &
     &  ,Thetar(Nx_Max,Ny_Max),ht(nx_max,ny_max),                       &
     &   dhodyn(Nx_Max,Ny_Max),dhodxn(Nx_Max,Ny_Max),                   &
     &  fake1(Nx_Max,Ny_Max),fake2(Nx_Max,Ny_Max),fake(nx_max,ny_max),  &
     &  uw(0:1000),H,period

        pi=3.1415926
        grav=9.8
        B0=1./8.

        rho=1000.0

        period=Pass_period

!c---  calculate Sxx Sxy Syy
! ---   Outside the surfzone

!cuse depth dwf and check for small or negative depths
!cKAH 7/28/05
        
        do j=1,Ny_wave
        do i=1,nx_wave
        Thetar(i,j)=pi/180.*Pass_Theta(i,j)
          sinth(i,j)= dsin(thetar(i,j))
          costh(i,j)= dcos(thetar(i,j))
          dwf(i,j)=dr(i,j)
!         if (dwf(i,j).lt..005) dwf(i,j)=.005
          xx(i)=X_Wave(i,j)
!c       pass_ibrk(i,j)=0
!c       if (pass_diss(i,j).gt.0.00001) Pass_ibrk(i,j)=1
        if (dwf(i,j).lt..005) dwf(i,j)=.005
        enddo
        enddo

! --- look for 1st breakpoint along each each j (suppose wave from left side)
        dx = dxr

        do j=1,ny_wave
          index = nx_wave
          iflag = 0
          do i=1,nx_wave
          Afake2(i)=0
            if(Pass_ibrk(i,j)==1.and.iflag==0)then
              index=i
              iflag=1
            endif
          enddo
          index=index-1
!c         index=nx_wave
          slope=(Depth_Wave(index+1,j)-Depth_wave(index,j))/dx
!jcw set slope
          slope=-0.0333

          jump=int(0.4/(-1.*slope)**1.1*Depth_wave(index,j)/dx)
!c         jump=3

!jcw  force jump bounded on 1 min
!           jump=1

!jcw
!      write(*,*) 'here 10',j,jump, slope,Depth_wave(index+1,j),dx,index,&
!     &                                   Depth_wave(index,j)

!cOutside the surfzone 
	 do i=1,index
            wavenumber = 2.d0*pi/period /pass_c(i,j)
            aux = 2*wavenumber*(dwf(i,j))
            capg = aux/SINH(aux)
            sm = 1.d0/16.d0*(1.d0+capg)
            aux = rho*pass_height(i,j)*pass_height(i,j)*grav
            sm  = sm*aux
            sp = .5d0*rho*grav*pass_height(i,j)                         &
     &           *pass_height(i,j)*capg*.125d0
            pass_sxx(i,j) = costh(i,j)*costh(i,j)*(sm) + sp 
            pass_sxy(i,j) = sinth(i,j)*costh(i,j)*(sm)
            pass_syy(i,j) = sinth(i,j)*sinth(i,j)*(sm)  + sp
            pass_massflux(i,j) =grav*pass_height(i,j)*pass_height(i,j)  &
     &           /pass_c(i,j)*B0
            Pass_MassFluxU(i,j)=Pass_MassFlux(i,j)*cos(Thetar(i,j))
            Pass_MassFluxV(i,j)=Pass_MassFlux(i,j)*sin(Thetar(i,j))
            qwrol(i)=0.d0
            smrol(i)=0.d0
            Area(i,j)=0.d0
            Afake(i)=0.d0
            Afake2(i)=0.d0
            xfake(i)=xx(i)
         end do

!cRoller, past transition
          do i=index+jump,nx_wave
            Area(i,j)=0.06d0/pass_height(i,j)*pass_c(i,j)*period
            Afake(i-jump+1)=Area(i,j)
            Afake2(i)=Area(i,j)
            xfake(i-jump+1)=xx(i)
	 end do

!c        write(*,*)index,jump,nx_wave
!c       write(*,*)(xfake(i),i=1,nx_wave)

!c------ Do the cubic spline to for roller transition

         call spline1(nx_wave,ny_wave,xfake,Afake,nx_wave-jump,y2)
         call splint1(nx_wave,ny_wave,xfake,Afake,y2,nx_wave-jump,      &
     &        xx,Area,index,jump,j)

!c        call spline(xfake,Afake,nx_wave-jump,y2)
!c        call splint(xfake,Afake,y2,nx_wave-jump,
!c    .        xx,Afake2,index,jump,j)

!c        write(*,*)(Afake2(i),i=53,nx_wave)

!c------ Final value of Radiation stress and Volume flux in surfzone
          do i=index,nx_wave
!c           Area(i,j)=Afake2(i)
            wavenumber = 2.d0*pi/period /pass_c(i,j)
	    aux = 2*wavenumber*(dwf(i,j))
	    capg = aux/SINH(aux)
            sm = 1.d0/16.d0*(1.d0+capg)*rho*pass_height(i,j)            &
     &           *pass_height(i,j)*grav
            qwrol(i)=Area(i,j)/period*pass_height(i,j)*pass_height(i,j)
            smrol(i)=rho*pass_c(i,j)*Area(i,j)/period                   &
     &           *pass_height(i,j)*pass_height(i,j)
            sp = .5d0*rho*grav*pass_height(i,j)*                        &
     &           pass_height(i,j)*capg*.125d0
	    pass_sxx(i,j) = costh(i,j)*costh(i,j)*(sm+smrol(i)) + sp 
	    pass_sxy(i,j) = sinth(i,j)*costh(i,j)*(sm+smrol(i))
	    pass_syy(i,j) = sinth(i,j)*sinth(i,j)*(sm+smrol(i))  + sp
            pass_massflux(i,j) = B0*grav*pass_height(i,j)               &
     &           *pass_height(i,j)/pass_c(i,j)+qwrol(i)
            Pass_MassFluxU(i,j)=Pass_MassFlux(i,j)*cos(Thetar(i,j))
            Pass_MassFluxV(i,j)=Pass_MassFlux(i,j)*sin(Thetar(i,j))

	 end do

         end do


! --- cdiss, cubott

        gam = 0.3

        do j=1,ny_wave
        do i=1,nx_wave
        Pass_Diss(i,j)=abs(0.125*grav*Pass_Height(i,j)                  &
     &    *Pass_Height(i,j)*0.15                                        &
!c     *    *Pass_Cg(i,j)/Depth_wave(i,j)                               &
!c     *    *(1.-(gam*Depth_wave(i,j)/Pass_Height(i,j))**2)             &
     &    *Pass_Cg(i,j)/dwf(i,j)                                        &
     &    *(1.-(gam*dwf(i,j)/Pass_Height(i,j))**2)                      &
     &    *Pass_ibrk(i,j))

        Pass_ubott(i,j)=Pass_Height(i,j)*pi/period                      &
!c    *       /sinh(Depth_wave(i,j)*2.*pi                               &
     &       /sinh(dwf(i,j)*2.*pi                                       &
     &       /period/Pass_C(i,j))
        enddo
        enddo

        
        do j=1,ny_wave
           do i=1,nx_wave
           pass_area(i,j)=area(i,j)*pass_height(i,j)**2
        end do
        end do

203    format(501(f10.4))
           
!          open(unit=19,file='ubwave.dat', status='old')
!        do j=1,ny_wave
!        read(19,*)(Pass_ubott(i,j),i=nx_wave,1,-1)
!        end do
!        close(19)
     





      ncode = 3
      nxflag = 1

         do i = 1, Nx_wave
          do j = 1, Ny_wave
             ht(i,j)=Depth_Sedi(i,j)
          end do
         end do

      call derivsed(ht,fake,fake,dhodxn,dhodyn,fake1,fake2,             &
     &           ncode,nxflag,nx_wave,ny_wave,dxr,dyr)



104    format(201(f12.4))

       

       do j=1,ny_wave
       do i=1,nx_wave
!c    Also determine wave velocity based on Elfrink method

!c         H=pass_height(i,j)

!c       if(H.ge.0.5*ht(i,j).and.i.ge.nx_wave-1) 
!c    .            H=0.5*ht(i,j)
!c        call orbit7(100,H,ht(i,j),period
!c    .            ,-dhodxn(i,j),uw)

      
         do ii=1,100
            pass_uw(i,j,ii)=pass_ubott(i,j)*sin(pi/50*ii)
!c           pass_uw(i,j,ii)=uw(ii)
           
         end do
         end do
         end do




!c--- print sxx sxy and syy, ibrk

      if(fname10.ne.'  ') then
        do i = 1, nx_wave 
         write(10,100)(Pass_Sxx(i,j), j=1,ny_wave)
        end do     

        do i = 1, nx_wave 
         write(11,100)(Pass_Sxy(i,j), j=1,ny_wave ) 
        end do     

        do i = 1, nx_wave
         write(12,100)(Pass_Syy(i,j), j=1,ny_wave ) 
        end do     

        endif

!c--- write Fx Fy

      if(fname13.ne.'  ') then
        do i = 1, nx_wave 
         write(13,100)(Pass_Wave_Fx(i,j), j=1,ny_wave)
        end do     

        do i = 1, nx_wave 
         write(14,100)(Pass_Wave_Fy(i,j), j=1,ny_wave ) 
        end do     

        endif

!c--- write ibrk

      if(fname19.ne.'  ') then  
        do i = 1, nx_wave
         write(19,100)(Pass_ibrk(i,j), j=1,ny_wave ) 
        end do  
       endif

!c--- write Quv

      if(fname15.ne.'  ') then  
        do i = 1, nx_wave
         write(15,100)(Pass_MassFluxU(i,j), j=1,ny_wave ) 
        end do  
        do i = 1, nx_wave
         write(16,100)(Pass_MassFluxV(i,j), j=1,ny_wave ) 
        end do  
       endif

!c--- write tb

      if(fname17.ne.'  ') then  
        do i = 1, nx_wave
         write(17,100)(Pass_ubott(i,j), j=1,ny_wave ) 
        end do  
      endif

!c--- write local radiation stresses for 3d circ model

       if(fname21.ne.'  ') then 

        do i = 1, nx_wave
         write(21,100)(Pass_Sxx_body(i,j), j=1,ny_wave ) 
        end do  
        do i = 1, nx_wave
         write(22,100)(Pass_Sxy_body(i,j), j=1,ny_wave ) 
        end do  
        do i = 1, nx_wave
         write(23,100)(Pass_Syy_body(i,j), j=1,ny_wave ) 
        end do  
        do i = 1, nx_wave
         write(24,100)(Pass_Sxx_surf(i,j), j=1,ny_wave ) 
        end do  
        do i = 1, nx_wave
         write(25,100)(Pass_Sxy_surf(i,j), j=1,ny_wave ) 
        end do  
        do i = 1, nx_wave
         write(27,100)(Pass_Syy_surf(i,j), j=1,ny_wave ) 
        end do  

       endif

  101  format(501(I3))
  100  format(501(f10.4)) 
        
        
        return
        end subroutine



        subroutine spline1(nx,ny,xin,yin,n,y2)

! --- This subroutine will compute the cubic spline for an array.
!        include 'pass.h'    
#include "pass.h"

        integer i,j

        real*8 xin(nx_max),yin(nx_max),y2(nx_max),u(nx_max)    

      y2(1)=0.
      u(1)=0.
      do i=2,n-1
         sig=(xin(i)-xin(i-1))/(xin(i+i)-xin(i))
         p=sig*y2(i-1)+2.
         y2(i)=(sig-1.)/p
         u(i)=(6.*((yin(i+1)-yin(i))/(xin(i+1)-xin(i))                  &
     &        -(yin(i)-yin(i-1))                                        &
     &        /(xin(i)-xin(i-1)))/(xin(i+1)-xin(i-1))-sig*u(i-1))/p
      end do
      qn=0.
      un=0.
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
      
      do kk=n-1,1,-1
         y2(kk)=y2(kk)*y2(kk+1)+u(kk)
      end do

        return
      end  subroutine



      SUBROUTINE SPLINT1(nx,ny,xin,yin,y2,n,xout,yout,index,jump,j)

!        include 'pass.h'
#include "pass.h"

        integer i,j

        real*8 xin(nx_max),yin(nx_max),y2(nx_max),xout(nx_max)
        real*8 yout(nx_max,ny_max)
        integer index,jump

        do i=1,index-1
        yout(i,j)=yin(i)
        enddo

        do i=index+jump+1,nx
        yout(i,j)=yin(i-jump+1)
        enddo

      do i=index,index+jump
      KLO=1
      KHI=N
1     IF (KHI-KLO.GT.1) THEN
        kk=(KHI+KLO)/2
        IF(xin(kk).GT.xout(i))THEN
          KHI=kk
        ELSE
          KLO=kk
        ENDIF
      GOTO 1
      ENDIF
      hh=xin(KHI)-xin(KLO)
      IF (hh.EQ.0.) then
        print*,'Bad xin input'
        stop
      endif
      aa=(xin(KHI)-xout(i))/hh
      bb=(xout(i)-xin(KLO))/hh
      yout(i,j)=aa*yin(KLO)+bb*yin(KHI)+                                &
     &      ((aa**3-aa)*y2(KLO)+(bb**3-bb)*y2(KHI))*(hh**2)/6.

      end do
        return
      end  subroutine



!c#####################################################################
!cDERIVATIVES
      subroutine derivsed(zeta,qx,qy,dzetadx,dzetady,dqxdx,             &
     &                 dqydy,ncode,nxflag,nx,ny,dx,dy)
!c	LAST MODIFICATION:						
!c	- Sep 15, 1999							
!c#####################################################################
!CDIR$ FLOW


!c--- Local variables -------------------------------------------------
      integer i, j, ncode, nxflag,nxj(300),minnxa,nx,ny,nxm,nym,nmax
       parameter (nxm = 300, nym = 300, nmax=300, nsensor=20)
      real*8 zeta(nxm,nym), qx(nxm,nym), qy(nxm,nym),                   &
     & dzetadx(nxm,nym), dzetady(nxm,nym),                              &
     & dqxdx(nxm,nym), dqydy(nxm,nym),                                  &
     & spz(nmax),spqx(nmax),spqy(nmax),dzeta(nmax),                     &
     & dqx(nmax),dqy(nmax),dx,dy

!c=====================================================================

!c##  set number of points where to take derivative based on flag
!c---- nxflag=0 : perform derivatives until nx
!c---- nxflag=1 : perform derivatives until last wet point
      if (nxflag .eq. 0) then
	 minnxa=nx
         do j = 1,ny
	    nxj(j) = nx
         end do
      else
	 minnxa=nx
         do j = 1,ny
	    nxj(j) = nx
         end do
      end if	 

!c##  Calculating the x-derivatives ----------------------------------

      if (ncode.eq.1) then
!c------- 1st derivatives of higher order : dx^4 ---
         do j = 1,ny
            do i = 1, nxj(j)
               spz(i)  = zeta(i,j)
               spqx(i) = qx(i,j)
               dzeta(i) = 0.d0
               dqx(i) = 0.d0
            end do
!c           call splinex(nxj(j), spz, dzeta)
!c           call splinex(nxj(j), spqx, dqx)
            do i = 1, nxj(j)
               dzetadx(i,j) = dzeta(i)
               dqxdx(i,j)   = dqx(i)
            end do
            do i = nxj(j)+1,nx
               dzetadx(i,j) = 0.d0
               dqxdx(i,j)   = 0.d0
            end do
	 end do
      elseif (ncode.eq.3) then
!c------- 1st derivatives of lower order : dx^2 ---

         do j = 1,ny
            do i = 1, nxj(j)
               spz(i)  = zeta(i,j)
               spqx(i) = qx(i,j)
               dzeta(i) = 0.d0
            end do




            call splinex2sed(nxj(j), spz, dzeta,dx)
            call splinex2sed(nxj(j), spqx, dqx,dx)

            do i = 1, nxj(j)
               dzetadx(i,j) = dzeta(i)
               dqxdx(i,j)   = dqx(i)
            end do

            do i = nxj(j)+1,nx
               dzetadx(i,j) = 0.d0
               dqxdx(i,j)   = 0.d0
            end do

	 end do

      elseif (ncode.eq.2) then
!c------- 3rd derivatives ---
         do j = 1,ny
            do i = 1, nxj(j)
               spz(i)  = zeta(i,j)
               spqx(i) = qx(i,j)
               dzeta(i) = 0.d0
               dqx(i) = 0.d0
            end do
!c           call triplex(nxj(j), spz, dzeta)
!c           call triplex(nxj(j), spqx, dqx)
            do i = 1, nxj(j)
               dzetadx(i,j) = dzeta(i)
               dqxdx(i,j)   = dqx(i)
            end do
            do i = nxj(j)+1,nx
               dzetadx(i,j) = 0.d0
               dqxdx(i,j)   = 0.d0
            end do
	 end do
      elseif (ncode.eq.4) then
!c------- 2nd derivatives ---
         do j = 1,ny
            do i = 1, nxj(j)
               spz(i)  = zeta(i,j)
               spqx(i) = qx(i,j)
               dzeta(i) = 0.d0
               dqx(i) = 0.d0
            end do
!c           call doublex(nxj(j), spz, dzeta)
!c           call doublex(nxj(j), spqx, dqx)
            do i = 1, nxj(j)
               dzetadx(i,j) = dzeta(i)
               dqxdx(i,j)   = dqx(i)
            end do
            do i = nxj(j)+1,nx
               dzetadx(i,j) = 0.d0
               dqxdx(i,j)   = 0.d0
            end do
	 end do
      end if


!c##  Calculating the y-derivatives ----------------------------------

      if (ncode.eq.1) then
!c------- 1st derivatives  of higher order : dy^4 ---
         do i = 1,minnxa
            do j = 1,ny
               spz(j)  = zeta(i,j)
               spqy(j) = qy(i,j)
               dzeta(i) = 0.d0
               dqy(i) = 0.d0
            end do
!c           call spliney(ny, spz, dzeta)
!c           call spliney(ny, spqy, dqy)
            do j = 1,ny
               dzetady(i,j) = dzeta(j)
               dqydy(i,j)   = dqy(j)
            end do
         end do
         do j = 1,ny
               do i = nxj(j)+1,nx
                dzetady(i,j) = 0.d0
                dqydy(i,j)   = 0.d0
               end do
         end do
      elseif (ncode.eq.3) then
!c------- 1st derivatives  of lower order : dy^2 ---
         do i = 1,minnxa
            do j = 1,ny
               spz(j)  = zeta(i,j)
               spqy(j) = qy(i,j)
               dzeta(i) = 0.d0
               dqy(i) = 0.d0
            end do
            call spliney2sed(ny, spz, dzeta,dy)
            call spliney2sed(ny, spqy, dqy,dy)
            do j = 1,ny
               dzetady(i,j) = dzeta(j)
               dqydy(i,j)   = dqy(j)
            end do
         end do
            do j = 1,ny
               do i =nxj(j)+1,nx
                 dzetady(i,j) = 0.d0
                 dqydy(i,j)   = 0.d0
            end do
         end do
      elseif (ncode.eq.2) then
!c------- 3rd derivatives ---
         do i = 1,minnxa
            do j = 1,ny
               spz(j)  = zeta(i,j)
               spqy(j) = qy(i,j)
               dzeta(i) = 0.d0
               dqy(i) = 0.d0
            end do
!c           call tripley(ny, spz, dzeta)
!c           call tripley(ny, spqy, dqy)
            do j = 1,ny
               dzetady(i,j) = dzeta(j)
               dqydy(i,j)   = dqy(j)
            end do
         end do
            do j = 1,ny
                do i = nxj(j)+1,nx
                 dzetady(i,j) = 0.d0
                 dqydy(i,j)   = 0.d0
            end do
         end do
      elseif (ncode.eq.4) then
!c------- 2nd derivatives ---
         do i = 1,minnxa
            do j = 1,ny
               spz(j)  = zeta(i,j)
               spqy(j) = qy(i,j)
               dzeta(i) = 0.d0
               dqy(i) = 0.d0
            end do
!c           call doubley(ny, spz, dzeta)
!c           call doubley(ny, spqy, dqy)
            do j = 1,ny
               dzetady(i,j) = dzeta(j)
               dqydy(i,j)   = dqy(j)
            end do
         end do
            do j = 1,ny
              do i = nxj(j)+1,nx
               dzetady(i,j) = 0.d0
               dqydy(i,j)   = 0.d0
            end do
         end do
      end if


      return
      end subroutine

!c#####################################################################
!cSPLINE X, 2nd Order
      subroutine splinex2sed(ndata,fdata,fxdata,dx)
!c#####################################################################
!CDIR$ FLOW

!c--- Local variables -------------------------------------------------
      integer i, ndata
      real*8 fdata(ndata), fxdata(ndata), dernum,dx

!c=====================================================================

      do i = 2,ndata-1
         dernum = fdata(i+1) - fdata(i-1)
         fxdata(i) = dernum/(2.d0*dx)
      end do

      if (iperx.eq.1) then
         dernum = fdata(2) - fdata(ndata-1)
         fxdata(1) = dernum/(2.d0*dx)
         fxdata(ndata) = fxdata(1)
      else
         i = 1
         dernum = -3.d0*fdata(i) + 4.d0*fdata(i+1) - fdata(i+2)
         fxdata(i) = dernum/(2.d0*dx)
!c        dernum = -fdata(i) + fdata(i+1)
!c        fxdata(i) = dernum/dx
         i = ndata
         dernum = 3.d0*fdata(i) - 4.d0*fdata(i-1) + fdata(i-2)
         fxdata(i) = dernum/(2.d0*dx)
!c        dernum = fdata(i) - fdata(i-1)
!c        fxdata(i) = dernum/dx
      end if

      return
      end subroutine
!c#####################################################################
!cSPLINE Y, 2nd Order
      subroutine spliney2sed(ndata,fdata,fxdata,dy)
!c#####################################################################
!CDIR$ FLOW
!c     include 'winc_std_common.inc'
!c--- Local variables -------------------------------------------------
      integer j, ndata
      real*8 fdata(ndata), fxdata(ndata), dernum,dy

!c=====================================================================

      do j = 2,ndata-1
         dernum = fdata(j+1) - fdata(j-1)
         fxdata(j) = dernum/(2.d0*dy)
      end do

      if (ipery.eq.1) then
         dernum = fdata(2) - fdata(ndata-1)
         fxdata(1) = dernum/(2.d0*dy)
         fxdata(ndata) = fxdata(1)
      else
         j=1
         dernum = -3.d0*fdata(j) + 4.d0*fdata(j+1) - fdata(j+2)
         fxdata(1) = dernum/(2.d0*dy)
         j = ndata
         dernum = 3.d0*fdata(j) - 4.d0*fdata(j-1) + fdata(j-2)
         fxdata(j) = dernum/(2.d0*dy)
      end if

      return
      end subroutine

#else
      MODULE refdif_mod
      END MODULE refdif_mod
#endif
!      END MODULE refdif_mod
