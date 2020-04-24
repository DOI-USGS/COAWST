/*
** svn $Id: vegetation_def.h 429 2009-12-20 17:30:26Z arango $
*************************************************** Hernan G. Arango ***
** Copyright (c) 2002-2017 The ROMS/TOMS Group                        **
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
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
# endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idWdvg,ng),r8)
          status=def_var(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idWdvg),   &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
#endif 
!
#ifdef MARSH_DYNAMICS
!
!  Store marsh masking marsh from marsh cells. 
!
          Vinfo( 1)=Vname(1,idTims)
          Vinfo( 2)=Vname(2,idTims)
          Vinfo( 3)=Vname(3,idTims)
          Vinfo(14)=Vname(4,idTims)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
#    if defined PERFECT_RESTART
        Vinfo(24)='_FillValue'
        Aval(6)=spval
#    else
          Vinfo(20)='mask_rho'
#    endif
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idTims,ng),r8)
          status=def_var(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idTims),   &
     &                   NF_FRST, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
!
# ifdef MARSH_WAVE_THRUST
!
!  Total thrust from all directions due to waves.
!
          Vinfo( 1)=Vname(1,idTtot)
          Vinfo( 2)=Vname(2,idTtot)
          Vinfo( 3)=Vname(3,idTtot)
          Vinfo(14)=Vname(4,idTtot)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
#    if defined PERFECT_RESTART
        Vinfo(24)='_FillValue'
        Aval(6)=spval
#    else
          Vinfo(20)='mask_rho'
#    endif
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idTtot,ng),r8)
          status=def_var(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idTtot),   &
     &                   NF_FRST, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
!
#  ifdef MARSH_SED_EROSION
!
!  Marsh sediment flux out from marsh cells from each sedclass.
!
        DO i=1,NST
          Vinfo( 1)=Vname(1,idTmfo(i))
          Vinfo( 2)=Vname(2,idTmfo(i))
          Vinfo( 3)=Vname(3,idTmfo(i))
          Vinfo(14)=Vname(4,idTmfo(i))
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
#    if defined PERFECT_RESTART
          Vinfo(24)='_FillValue'
          Aval(6)=spval
#    else
          Vinfo(20)='mask_rho'
#    endif
#   endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idTmfo(i),ng),r8)
          status=def_var(ng, iNLM, RST(ng)%ncid,                        &
     &                   RST(ng)%Vid(idTmfo(i)), NF_FRST,               &
     &                   nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END DO
!
#   ifdef MARSH_RETREAT
!
!  Amount of marsh retreat from all four directions.
!
          Vinfo( 1)=Vname(1,idTmmr)
          Vinfo( 2)=Vname(2,idTmmr)
          Vinfo( 3)=Vname(3,idTmmr)
          Vinfo(14)=Vname(4,idTmmr)
          Vinfo(16)=Vname(1,idtime)
#    if defined WRITE_WATER && defined MASKING
#     if defined PERFECT_RESTART
        Vinfo(24)='_FillValue'
        Aval(6)=spval
#     else
          Vinfo(20)='mask_rho'
#     endif
#    endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idTmmr,ng),r8)
          status=def_var(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idTmmr),   &
     &                   NF_FRST, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
#   endif 
#  endif  
# endif   
#endif
