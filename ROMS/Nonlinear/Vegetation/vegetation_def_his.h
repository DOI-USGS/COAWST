/*
** svn $Id: vegetation_def.h 429 2019-01-01 17:30:26Z arango $
*************************************************** Hernan G. Arango ***
** Copyright (c) 2002-2019 The ROMS/TOMS Group                        **
**   Licensed under a MIT/X style license                             **
**   See License_ROMS.txt                                             **
*************************************************** John C. Warner    **
*************************************************** Neil K. Ganju     **
*************************************************** Alexis Beudin     **
*************************************************** Tarandeep S. Kalra**
**                                                                    **
**  Defines vegetation module input parameters in output NetCDF files.**
**  It is included in routine "def_his.F".                            **
**                                                                    **
************************************************************************
*/
#if defined VEG_DRAG || defined VEG_BIOMASS
!
!  Define vegetation module parameters.
!
      DO i=1,NVEGP
        IF (Hout(idvprp(i),ng)) THEN
           Vinfo( 1)=Vname(1,idvprp(i))
           Vinfo( 2)=Vname(2,idvprp(i))
           Vinfo( 3)=Vname(3,idvprp(i))
           Vinfo(14)=Vname(4,idvprp(i))
           Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
           Vinfo(20)='mask_rho'
#  endif
           Vinfo(21)=Vname(6,idvprp(i))
!          Vinfo(22)='coordinates'
           Aval(5)=REAL(Iinfo(1,idvprp(i),ng),r8)
           status=def_var(ng, iNLM, HIS(ng)%ncid,HIS(ng)%Vid(idvprp(i))  &
     &                   ,NF_FOUT, nvd4, v3pgrd, Aval, Vinfo, ncname,    &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
      END DO
#endif
#if defined VEG_STREAMING
!
!  Define wave dissipation due to vegetation.
!
        IF (Hout(idWdvg,ng)) THEN
          Vinfo( 1)=Vname(1,idWdvg)
          Vinfo( 2)=Vname(2,idWdvg)
          Vinfo( 3)=Vname(3,idWdvg)
          Vinfo(14)=Vname(4,idWdvg)
!         Vinfo(16)=Vname(1,idWdvg)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
# endif
          Vinfo(21)=Vname(6,idWdvg)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idWdvg,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idWdvg),   &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname,    &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define spectral Cd due to wave vegetation.
!
        IF (Hout(idCdvg,ng)) THEN
          Vinfo( 1)=Vname(1,idCdvg)
          Vinfo( 2)=Vname(2,idCdvg)
          Vinfo( 3)=Vname(3,idCdvg)
          Vinfo(14)=Vname(4,idCdvg)
!         Vinfo(16)=Vname(1,idCdvg)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
# endif
          Vinfo(21)=Vname(6,idCdvg)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idCdvg,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idCdvg),   &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname,    &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        END IF
#endif
#ifdef MARSH_DYNAMICS
!
!  Store masking marsh of marsh cells.
!
       IF (Hout(idTims,ng)) THEN
         Vinfo( 1)=Vname(1,idTims)
         Vinfo( 2)=Vname(2,idTims)
         Vinfo( 3)=Vname(3,idTims)
         Vinfo(14)=Vname(4,idTims)
         Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
         Vinfo(20)='mask_rho'
#  endif
         Vinfo(21)=Vname(6,idTims)
         Vinfo(22)='coordinates'
         Aval(5)=REAL(Iinfo(1,idTims,ng),r8)
         status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idTims),    &
     &                  NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname ,    &
     &                  SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
       END IF
# ifdef MARSH_WAVE_THRUST
!
!  Total thrust from all directions due to waves.
!
        IF (Hout(idTtot,ng)) THEN
          Vinfo( 1)=Vname(1,idTtot)
          Vinfo( 2)=Vname(2,idTtot)
          Vinfo( 3)=Vname(3,idTtot)
          Vinfo(14)=Vname(4,idTtot)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(21)=Vname(6,idTtot)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idTtot,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idTtot),   &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname,    &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
#  ifdef MARSH_SED_EROSION
!
!  Marsh sediment flux out from marsh cells from each sedclass type.
!
        DO i=1,NST
          IF (Hout(idTmfo(i),ng)) THEN
            Vinfo( 1)=Vname(1,idTmfo(i))
            Vinfo( 2)=Vname(2,idTmfo(i))
            Vinfo( 3)=Vname(3,idTmfo(i))
            Vinfo(14)=Vname(4,idTmfo(i))
            Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
            Vinfo(20)='mask_rho'
#   endif
            Vinfo(21)=Vname(6,idTmfo(i))
            Vinfo(22)='coordinates'
            Aval(5)=REAL(Iinfo(1,idTmfo(i),ng))
            status=def_var(ng, iNLM, HIS(ng)%ncid,                      &
     &                     HIS(ng)%Vid(idTmfo(i)), NF_FOUT,             &
     &                     nvd3, t2dgrd, Aval, Vinfo, ncname,           &
     &                     SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END DO
!
#   ifdef MARSH_RETREAT
!
!  Amount of marsh retreat from all directions.
!
        IF (Hout(idTmmr,ng)) THEN
          Vinfo( 1)=Vname(1,idTmmr)
          Vinfo( 2)=Vname(2,idTmmr)
          Vinfo( 3)=Vname(3,idTmmr)
          Vinfo(14)=Vname(4,idTmmr)
          Vinfo(16)=Vname(1,idtime)
#    if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#    endif
          Vinfo(21)=Vname(6,idTmmr)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idTmmr,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idTmmr),   &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname,    &
     &                     SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#   endif
#  endif
# endif
# ifdef MARSH_TIDAL_RANGE_CALC
!
!  Amount of marsh tidal range over a given frequency.
!
        IF (Hout(idTmtr,ng)) THEN
          Vinfo( 1)=Vname(1,idTmtr)
          Vinfo( 2)=Vname(2,idTmtr)
          Vinfo( 3)=Vname(3,idTmtr)
          Vinfo(14)=Vname(4,idTmtr)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(21)=Vname(6,idTmtr)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idTmtr,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idTmtr),   &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname,    &
     &                     SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
# endif
!
# if defined MARSH_VERT_GROWTH
!
        IF (Hout(idTmhw,ng)) THEN
          Vinfo( 1)=Vname(1,idTmhw)
          Vinfo( 2)=Vname(2,idTmhw)
          Vinfo( 3)=Vname(3,idTmhw)
          Vinfo(14)=Vname(4,idTmhw)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(21)=Vname(6,idTmhw)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idTmhw,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idTmhw),   &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname,    &
     &                     SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
	END IF
!
!  Amount of marsh mean low water over a given frequency.
!
        IF (Hout(idTmlw,ng)) THEN
          Vinfo( 1)=Vname(1,idTmlw)
          Vinfo( 2)=Vname(2,idTmlw)
          Vinfo( 3)=Vname(3,idTmlw)
          Vinfo(14)=Vname(4,idTmlw)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(21)=Vname(6,idTmlw)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idTmlw,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idTmlw),   &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname,    &
     &                     SetFillVal = .FALSE.)
         IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
	END IF
!
!
!  Amount of marsh biomass peak (kg/sq.m).
!
        IF (Hout(idTmbp,ng)) THEN
          Vinfo( 1)=Vname(1,idTmbp)
          Vinfo( 2)=Vname(2,idTmbp)
          Vinfo( 3)=Vname(3,idTmbp)
          Vinfo(14)=Vname(4,idTmbp)
          Vinfo(16)=Vname(1,idtime)
#    if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#    endif
          Vinfo(21)=Vname(6,idTmbp)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idTmbp,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idTmbp),   &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname,    &
     &                     SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Rate of marsh vertical growth (m/year).
!
        IF (Hout(idTmvg,ng)) THEN
          Vinfo( 1)=Vname(1,idTmvg)
          Vinfo( 2)=Vname(2,idTmvg)
          Vinfo( 3)=Vname(3,idTmvg)
          Vinfo(14)=Vname(4,idTmvg)
          Vinfo(16)=Vname(1,idtime)
#    if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#    endif
          Vinfo(21)=Vname(6,idTmvg)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idTmvg,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idTmvg),   &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname,    &
     &                     SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Amount of marsh vertical growth (m).
!
        IF (Hout(idTmvt,ng)) THEN
          Vinfo( 1)=Vname(1,idTmvt)
          Vinfo( 2)=Vname(2,idTmvt)
          Vinfo( 3)=Vname(3,idTmvt)
          Vinfo(14)=Vname(4,idTmvt)
          Vinfo(16)=Vname(1,idtime)
#    if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#    endif
          Vinfo(21)=Vname(6,idTmvt)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idTmvt,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idTmvt),   &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname,    &
     &                     SetFillVal = .FALSE.)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
# endif
#endif
