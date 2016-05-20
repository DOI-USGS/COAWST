/*
** svn $Id: vegetation_def.h 429 2009-12-20 17:30:26Z arango $
*************************************************** Hernan G. Arango ***
** Copyright (c) 2002-2016 The ROMS/TOMS Group                        **
**   Licensed under a MIT/X style license                             **
**   See License_ROMS.txt                                             **
*************************************************** John C. Warner    **
*************************************************** Neil K. Ganju     **
*************************************************** Alexis Beudin     **
*************************************************** Tarandeep S. Kalra**
**                                                                    **
**  Defines vegetation module input parameters in output restart      ** 
**  NetCDF files.                                                     **
**  It is included in routine "def_rst.F".                            **
**                                                                    **
************************************************************************
*/
!
!  Define vegetation module parameters.
!
#if defined VEG_DRAG || defined VEG_BIOMASS
      DO i=1,NVEGP
           Vinfo( 1)=Vname(1,idvprp(i))
           Vinfo( 2)=Vname(2,idvprp(i))
           Vinfo( 3)=Vname(3,idvprp(i))
           Vinfo(14)=Vname(4,idvprp(i))
           Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
#   if defined PERFECT_RESTART
          Vinfo(24)='_FillValue'
          Aval(6)=spval
#   else
          Vinfo(20)='mask_rho'
#   endif
#  endif
!          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idvprp(i),ng),r8)
          status=def_var(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idvprp(i)),&
     &                  NF_FRST, nvd4, v3pgrd, Aval, Vinfo, ncname)
         IF (exit_flag.ne.NoError) RETURN
      END DO
#endif
!
#if defined VEG_STREAMING 
!
!  Define wave dissipation due to vegetation 
!
          Vinfo( 1)=Vname(1,idWdvg)
          Vinfo( 2)=Vname(2,idWdvg)
          Vinfo( 3)=Vname(3,idWdvg)
          Vinfo(14)=Vname(4,idWdvg)
          Vinfo(16)=Vname(1,idWdvg)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
# endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idWdvg,ng),r8)
          status=def_var(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idWdvg),   &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
#endif 
#ifdef MARSH_WAVE_THRUST
!
!  Store initial masking marsh 
!
          Vinfo( 1)=Vname(1,idTims)
          Vinfo( 2)=Vname(2,idTims)
          Vinfo( 3)=Vname(3,idTims)
          Vinfo(14)=Vname(4,idTims)
          Vinfo(16)=Vname(1,idTims)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
# endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idTims,ng),r8)
          status=def_var(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idTims),   &
     &                   NF_FRST, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
!
          Vinfo( 1)=Vname(1,idTmsk)
          Vinfo( 2)=Vname(2,idTmsk)
          Vinfo( 3)=Vname(3,idTmsk)
          Vinfo(14)=Vname(4,idTmsk)
          Vinfo(16)=Vname(1,idTmsk)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
# endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idTmsk,ng),r8)
          status=def_var(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idTmsk),   &
     &                   NF_FRST, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
!
!  Define maximum thrust due to waves.
!
          Vinfo( 1)=Vname(1,idTmax)
          Vinfo( 2)=Vname(2,idTmax)
          Vinfo( 3)=Vname(3,idTmax)
          Vinfo(14)=Vname(4,idTmax)
          Vinfo(16)=Vname(1,idTmax)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
# endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idTmax,ng),r8)
          status=def_var(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idTmax),   &
     &                   NF_FRST, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
!     
!  Define Tonelli masking based thrust due to waves.
!
          Vinfo( 1)=Vname(1,idTton)
          Vinfo( 2)=Vname(2,idTton)
          Vinfo( 3)=Vname(3,idTton)
          Vinfo(14)=Vname(4,idTton)
          Vinfo(16)=Vname(1,idTton)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
# endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idTton,ng),r8)
          status=def_var(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idTton),   &
     &                   NF_FRST, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
#endif 
