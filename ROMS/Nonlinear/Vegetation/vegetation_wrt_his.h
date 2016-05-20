/*
** svn $Id: vegetation_wrt.h 429 2015-06-10 10:40:26Z arango $
*************************************************** Hernan G. Arango ***
** Copyright (c) 2002-2016 The ROMS/TOMS Group                        **
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
     &                       VEG(ng) % plant(:,:,:,i))
          IF (status.ne.nf90_noerr) THEN 
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
     &                     VEG(ng)%Dissip_veg )
        IF (status.ne.nf90_noerr) THEN 
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
# ifdef MARSH_WAVE_THRUST
!
!  Write out initial masking for marshes 
! 
      IF (Hout(idTims,ng)) THEN 
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idTims), &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % rmask,                            &
# endif
     &                     VEG(ng)%marsh_mask)
        IF (status.ne.nf90_noerr) THEN 
          IF (Master) THEN 
            WRITE (stdout,10) TRIM(Vname(1,idTims)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out wave thrust on marsh output 
! 
      IF (Hout(idTmsk,ng)) THEN 
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idTmsk), &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % rmask,                            &
# endif
     &                     VEG(ng)%mask_thrust)
        IF (status.ne.nf90_noerr) THEN 
          IF (Master) THEN 
            WRITE (stdout,10) TRIM(Vname(1,idTmsk)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Define maximum thrust due to waves. 
!
      IF (Hout(idTmax,ng)) THEN 
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idTmax), &
     &                     HIS(ng)%Rindex, gtype,                       &    
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % rmask,                            &
# endif
     &                     VEG(ng)%Thrust_max)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTmax)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Define maximum thrust due to waves. 
!
      IF (Hout(idTmax,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idTmax), &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % rmask,                            &
# endif
     &                     VEG(ng)%Thrust_max)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTmax)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Define Tonelli masking based thrust due to waves. 
!
      IF (Hout(idTton,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idTton), &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % rmask,                            &
# endif
     &                     VEG(ng)%Thrust_tonelli)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTton)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif 

