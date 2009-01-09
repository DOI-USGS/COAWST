      common/ref1/mr,nr,ispace,nd,md(ixr),iu,dconv(2),iff(3),icur,ibc,  &
     &          dconv2(2),ismooth,isp
      common/ref2/dr(ixr,iyr),ur(ixr,iyr),vr(ixr,iyr),iun(8),iinput,    &
     &        ioutput
      common/ref3/dxr,dyr,xr(ixr),yr(iyr),x(ix),y(iy)
      common/ref4/isd(ixr,iyr)
      common/block1/d(ix,iy),u(ix,iy),v(ix,iy),m,n,dx,dy,ibr(iy)
      common/con1/q(ix,iy),p(ix,iy),sig(ix,iy),bottomu(ix,iy)
      common/con2/k(ix,iy),kb(ix),w(ix,iy),dd(ix,iy),wb(2,iy)
      common/nlin/an,anl,ntype
      common/wav1/iwave,nfreqs,freqs(ncomp,numdata),                    &
     &        edens(ncomp),nwavs(ncomp),                                &
     &        nrs
	common/wav2/seed,amp(ncomp,numdata),dir(ncomp,numdata),         &
     &        tide(ncomp),                                              &
     &        thet0,update_interval,num_data,itime,num_total
      common/comp/a(ix,iy),psibar,ifilt
      common/names/fname1,fname2,fname3,fname4,fname5,fname6,fname7,    &
     &        fname8,                                                   &
     &        fname9,fname10,fname11,fname12,fname13,fname14,           &
     &        fname15,fname16,fname17,fname18,fname19,fname20,          &
     &        fname21,fname22,fname23,fname24,fname25,fname26,          &
     &        fnamein
	common/mpis/wav_comm_world,ocnid,wavid,ncomps
       real*8 k,kb,dr,ur,vr,dxr,dyr,xr,yr,x,y,d,u,v,dx,dy,dd,          &
     &         freqs,edens,amp,dir,seed
       complex*16 w, a, wb
      character*255 fname1,fname2,fname3,fname4,fname5,fname6,fname7,   &
     &              fname8,                                             &
     &              fname9,fname10,fname11,fname12,fname13,fname14,     &
     &              fname15,fname16,fname17,fname18,fname19,fname20,    &
     &              fname21,fname22,fname23,fname24,fname25,fname26,    &
     &              fnamein

!jcw
!	integer :: ifreq, ir, ncomps, mr, nr
        integer WAV_COMM_WORLD, OcnId, WavId
