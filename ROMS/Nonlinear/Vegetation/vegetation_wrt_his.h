/*
** svn $Id: vegetation_wrt.h 429 2015-06-10 10:40:26Z arango $
*************************************************** Hernan G. Arango ***
** Copyright (c) 2002-2017 The ROMS/TOMS Group                        **
**    See License_ROMS.txt                                            **
*************************************************** John C. Warner    **
*************************************************** Neil K. Ganju     **
*************************************************** Alexis Beudin     **
*************************************************** Tarandeep S. Kalra**
**                                                                    **
**  Writes vegetation input parameters into output NetCDF files.      **
**  It is included in routine "wrt_his.F".                           **
**                                                                    **
************************************************************************
*/
# if defined VEG_DRAG || defined VEG_BIOMASS
!
!  Write out vegetation properties 
! 
      DO i=1,NVEGP 
        IF (Hout(idvprp(i),ng)) THEN 
          scale=1.0_r8
          gtype=gfactor*r3dvar
          status=nf_fwrite3d(ng, iNLM, HIS(ng)%ncid,                    &
     &                       HIS(ng)%Vid(idvprp(i)),                    &
     &                       HIS(ng)%Rindex, gtype,                     &
     &                       LBi, UBi, LBj, UBj, 1, NVEG, scale,        &
#  ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#  endif
     &                       VEG(ng) % plant(:,:,:,i),                  &
     &                       SetFillVal= .FALSE.)
	  IF (FoundError(status, nf90_noerr, __LINE__,                    &
       &                   __FILE__)) THEN
            IF (Master) THEN 
              WRITE (stdout,10) TRIM(Vname(1,idvprp(i))), HIS(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO 
# endif 
!
# ifdef VEG_STREAMING 
!
!  Write out wave dissipation due to vegetation 
! 
      IF (Hout(idWdvg,ng)) THEN 
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idWdvg), &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % rmask,                            &
# endif
     &                     VEG(ng)%Dissip_veg)
        IF (FoundError(status, nf90_noerr, __LINE__,                    &
     &                   __FILE__)) THEN
          IF (Master) THEN 
            WRITE (stdout,10) TRIM(Vname(1,idWdvg)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif
! 
# ifdef MARSH_DYNAMICS
!
!  Write out masking for marsh cells. 
! 
      IF (Hout(idTims,ng)) THEN 
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idTims), &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     VEG(ng)%marsh_mask) 
        IF (FoundError(status, nf90_noerr, __LINE__,                    &
     &                   __FILE__)) THEN
          IF (Master) THEN 
            WRITE (stdout,10) TRIM(Vname(1,idTims)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  ifdef MARSH_WAVE_THRUST 
!
!  Write total thrust in all directions due to waves.
!
      IF (Hout(idTtot,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idTtot), &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     VEG(ng)%Thrust_total)                         
        IF (FoundError(status, nf90_noerr, __LINE__,                    &
     &                   __FILE__)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTtot)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  endif 
!
#  if defined MARSH_SED_EROSION
!
!  Write out marsh sediment flux out from marsh cells from each sedclass.
!
      DO i=1,NST
        IF (Hout(idTmfo(i),ng)) THEN
          scale=1.0_r8
          gtype=gfactor*r2dvar
          status=nf_fwrite2d(ng, iNLM, HIS(ng)%ncid,                    &
     &                       HIS(ng)%Vid(idTmfo(i)),                    &
     &                       HIS(ng)%Rindex, gtype,                     &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#   ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#   endif
     &                       VEG(ng) % marsh_flux_out(:,:,i))
          IF (FoundError(status, nf90_noerr, __LINE__,                  &
     &                   __FILE__)) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idTmfo(i))), HIS(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
#  endif 
!
#  if defined MARSH_RETREAT
!
!  Amount of marsh lateral retreat from all directions.
!
      IF (Hout(idTmmr,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idTmmr), &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     VEG(ng)%marsh_retreat)
        IF (FoundError(status, nf90_noerr, __LINE__,                    &
     &                   __FILE__)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTmmr)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF 
#  endif 
!# endif 
!
#  if defined MARSH_TIDAL_RANGE
!
!  Write tidal range over a given frequency.
!
      IF (Hout(idTmtr,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idTmtr), &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     VEG(ng)%marsh_tidal_range)
        IF (FoundError(status, nf90_noerr, __LINE__,                    &
     &                   __FILE__)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTmtr)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write mean high high water over a given frequency.
!
      IF (Hout(idTmhw,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idTmhw), &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     VEG(ng)%marsh_high_water) 
        IF (FoundError(status, nf90_noerr, __LINE__,                    &
     &                   __FILE__)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTmhw)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
#   if  defined MARSH_VERT_GROWTH 
!
!  Write amount of marsh biomass peak.
!
      IF (Hout(idTmbp,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idTmbp), &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#    ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#    endif
     &                     VEG(ng)%marsh_biomass_peak)
        IF (FoundError(status, nf90_noerr, __LINE__,                    &
     &                   __FILE__)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTmbp)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write amount of marsh vertical growth.
!
      IF (Hout(idTmvg,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idTmvg), &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#    ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#    endif
     &                     VEG(ng)%marsh_vert)
        IF (FoundError(status, nf90_noerr, __LINE__,                    &
     &                   __FILE__)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTmvg)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#   endif 
#  endif  
# endif  
