/*
** svn $Id: vegetation_wrt.h 429 2015-06-10 10:40:26Z arango $
*************************************************** Hernan G. Arango ***
** Copyright (c) 2002-2019 The ROMS/TOMS Group                        **
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
        IF (Hout(idvprp(i),ng)) THEN 
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
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idvprp(i))), RST(ng)%Rindex
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
        status=nf_fwrite2d(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idWdvg), &
     &                     RST(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % rmask,                            &
# endif
     &                     VEG(ng)%Dissip_veg )
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idWdvg)), RST(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out spectral Cd due to wave vegetation. 
! 
      IF (Hout(idCdvg,ng)) THEN 
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idCdvg), &
     &                     RST(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % rmask,                            &
# endif
     &                     VEG(ng)%Cdwave_veg )
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idCdvg)), RST(ng)%Rindex
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
!  Store marsh masking from marsh cells.
!
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idTims), &
     &                     RST(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     VEG(ng)%marsh_mask)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTims)), RST(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
!      END IF
#  ifdef MARSH_WAVE_THRUST 
!
!  Write total thrust in all directions due to waves.
!
!      IF (Hout(idTtot,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idTtot), &
     &                     RST(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     VEG(ng)%Thrust_total)
         IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
           IF (Master) THEN
             WRITE (stdout,10) TRIM(Vname(1,idTtot)), RST(ng)%Rindex
           END IF
           exit_flag=3
           ioerror=status
           RETURN
         END IF
!       END IF
!     END DO 
#  endif 
!
#  if defined MARSH_SED_EROSION
!
!  Write out marsh sediment flux out from marsh cells from each sedclass.
!
      DO i=1,NST
!        IF (Hout(idTmfo(i),ng)) THEN
          scale=1.0_r8
          gtype=gfactor*r2dvar
          status=nf_fwrite2d(ng, iNLM, RST(ng)%ncid,                    &
     &                     RST(ng)%Vid(idTmfo(i)),                      &
     &                     RST(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#    ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#    endif
     &                     VEG(ng) % marsh_flux_out(:,:,i))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTmfo(i))), RST(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END DO
#   endif
!
#   ifdef MARSH_RETREAT
!
!  Amount of marsh retreat from all directions.
!
!      IF (Hout(idTmmr,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idTmmr), &
     &                     RST(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     VEG(ng)%marsh_retreat)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTmmr)), RST(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
!      END IF 
!
#   endif
#  endif
!# endif

