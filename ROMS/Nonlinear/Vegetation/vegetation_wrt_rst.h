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
**  Writes vegetation input parameters into output restart            **
**  NetCDF files.                                                     **
**  It is included in routine "wrt_rst.F".                            **
**                                                                    **
************************************************************************
*/
# if defined VEG_DRAG || defined VEG_BIOMASS
!
!  Write out vegetation properties 
! 
      DO i=1,NVEGP 
          scale=1.0_r8
          gtype=gfactor*r3dvar
          status=nf_fwrite3d(ng, iNLM, RST(ng)%ncid,                    &
     &                       RST(ng)%Vid(idvprp(i)),                    &
     &                       RST(ng)%Rindex, gtype,                     &
     &                       LBi, UBi, LBj, UBj, 1, NVEG, scale,        &
#  ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#  endif
     &                       VEG(ng) % plant(:,:,:,i))
          IF (status.ne.nf90_noerr) THEN 
            IF (Master) THEN 
              WRITE (stdout,10) TRIM(Vname(1,idvprp(i))), RST(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
      END DO 
# endif 
!
# ifdef VEG_STREAMING 
!
!  Write out wave dissipation due to vegetation 
! 
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idWdvg), &
     &                     RST(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % rmask,                            &
# endif
     &                     VEG(ng)%Dissip_veg )
        IF (status.ne.nf90_noerr) THEN 
          IF (Master) THEN 
            WRITE (stdout,10) TRIM(Vname(1,idWdvg)), RST(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
# endif
! 
# ifdef MARSH_WAVE_THRUST
!
!  Write out initial masking for marshes 
! 
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idTims), &
     &                     RST(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % rmask,                            &
# endif
     &                     VEG(ng)%marsh_mask)
        IF (status.ne.nf90_noerr) THEN 
          IF (Master) THEN 
            WRITE (stdout,10) TRIM(Vname(1,idTims)), RST(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
!
!  Write out wave thrust on marsh output 
! 
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idTmsk), &
     &                     RST(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % rmask,                            &
# endif
     &                     VEG(ng)%mask_thrust)
        IF (status.ne.nf90_noerr) THEN 
          IF (Master) THEN 
            WRITE (stdout,10) TRIM(Vname(1,idTmsk)), RST(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
!
!  Define maximum thrust due to waves. 
!
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idTmax), &
     &                     RST(ng)%Rindex, gtype,                       &    
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % rmask,                            &
# endif
     &                     VEG(ng)%Thrust_max)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTmax)), RST(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
!
!  Define maximum thrust due to waves. 
!
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idTmax), &
     &                     RST(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % rmask,                            &
# endif
     &                     VEG(ng)%Thrust_max)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTmax)), RST(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
!
!  Define Tonelli masking based thrust due to waves. 
!
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idTton), &
     &                     RST(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % rmask,                            &
# endif
     &                     VEG(ng)%Thrust_tonelli)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTton)), RST(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
# endif 

