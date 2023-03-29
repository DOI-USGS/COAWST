/*
** svn $Id: vegetation_var.h 439 2015-06-10 12:00:00Z arango $
*************************************************** Hernan G. Arango ***
** Copyright (c) 2002-2019 The ROMS/TOMS Group                        **
**   Licensed under a MIT/X style license                             **
**   See License_ROMS.txt                                             **
*************************************************** John C. Warner    **
*************************************************** Neil K. Ganju     **
*************************************************** Alexis Beudin     **
*************************************************** Tarandeep S. Kalra**
**                                                                    **
**  Assigns metadata indices for the vegetation module variables that **
**  are used in input and output NetCDF files.  The metadata          **
**  information is read from "varinfo.dat".                           **
**                                                                    **
**  This file is included in file "mod_ncparam.F", routine            **
**  "initialize_ncparm".                                              **
**                                                                    **
************************************************************************
*/
#if defined VEG_DRAG || defined VEG_BIOMASS  
            CASE ('idvprp(pdens)')
                idvprp(pdens)=varid
            CASE ('idvprp(pdiam)')
                idvprp(pdiam)=varid
            CASE ('idvprp(pthck)')
                idvprp(pthck)=varid
            CASE ('idvprp(phght)')
                idvprp(phght)=varid
!#if defined VEG_BIOMASS  
!            CASE ('idvprp(pabbm)')
!                idvprp(pabbm)=varid
!            CASE ('idvprp(pbgbm)')
!                idvprp(pbgbm)=varid
!#endif
#endif 
#if defined VEG_STREAMING 
            CASE ('idWdvg')
              idWdvg=varid
            CASE ('idCdvg')
              idCdvg=varid
#endif 
! 
#if defined MARSH_DYNAMICS
            CASE ('idTims')
              idTims=varid
# if defined MARSH_WAVE_THRUST
            CASE ('idTtot')
              idTtot=varid
# endif 
# if defined MARSH_SED_EROSION 
            CASE ('idTmfo')
              load=.FALSE.
                IF ((NST.gt.0).and.                                     &
                  (Vinfo(1)(1:15).eq.'marsh_flux_out_')) THEN
                  varid=varid-1
                  DO i=1,NST
                    varid=varid+1
                    idTmfo(i)=varid
                    DO ng=1,Ngrids
                      Fscale(varid,ng)=scale
                      Iinfo(1,varid,ng)=gtype
                    END DO
                    WRITE (Vname(1,varid),'(a,i2.2)')                   &
     &                    TRIM(ADJUSTL(Vinfo(1))), i
                    WRITE (Vname(2,varid),'(a,a,i2.2)')                 &
     &                    TRIM(ADJUSTL(Vinfo(2))), ', size class ', i
                    WRITE (Vname(3,varid),'(a)')                        &
     &                    TRIM(ADJUSTL(Vinfo(3)))
                    WRITE (Vname(4,varid),'(a,a)')                      &
     &                    TRIM(Vname(1,varid)), ', scalar, series'
                    WRITE (Vname(5,varid),'(a)')                        &
     &                    TRIM(ADJUSTL(Vinfo(5)))
                  END DO
                  varid=varid+1
                END IF
# endif 
# if defined MARSH_RETREAT
            CASE ('idTmmr')
              idTmmr=varid
# endif 
# if defined MARSH_TIDAL_RANGE_CALC
	    CASE('idTmtr')
	      idTmtr=varid
# endif  
# if defined MARSH_VERT_GROWTH
	    CASE('idTmhw')
	      idTmhw=varid
	    CASE('idTmlw')
	      idTmlw=varid
	    CASE('idTmbp')
	      idTmbp=varid 
	    CASE('idTmvg')
	      idTmvg=varid 
	    CASE('idTmvt')
	      idTmvt=varid 
# endif  
#endif   
