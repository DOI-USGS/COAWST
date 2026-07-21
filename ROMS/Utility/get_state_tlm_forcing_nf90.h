/*
** git $Id$
*************************************************** Hernan G. Arango ***
** Copyright (c) 2002-2026 The ROMS Group                             **
**    Licensed under a MIT/X style license                            **
**    See License_ROMS.md                                             **
************************************************************************
**                                                                    **
** Standard NetCDF Library Processing:                                **
**                                                                    **
** It process the reading of ROMS TLM_FORCING state from input NetCDF **
** files.                                                             **
**                                                                    **
** It is used to read in 4D-Var adjustment corrections to the tangent **
** linear forcing terms. Data is loaded into "tl_" arrays.            **
**                                                                    **
************************************************************************
*/

!
!  Set switch to process surface forcing and/or open boundaries during
!  4D-Var minimization.
!
        get_adjust=.TRUE.

#ifdef ADJUST_BOUNDARY
!
!  Read in free-surface open boundaries adjustments.
!
        IF (get_var(idSbry(isFsur)).and.get_adjust.and.                 &
     &      ANY(Lobc(:,isFsur,ng))) THEN
          ifield=idSbry(isFsur)
          foundit=find_string(var_name, n_var, TRIM(Vname(1,ifield)),   &
     &                        varid)
          IF (foundit) THEN
            status=nf_fread2d_bry(ng, IDmod, ncname, ncINPid,           &
     &                            Vname(1,ifield), varid,               &
     &                            InpRec, r2dvar,                       &
     &                            LBij, UBij, Nbrec(ng),                &
     &                            Fscl, Fmin, Fmax,                     &
# ifdef CHECKSUM
     &                            BOUNDARY(ng) % tl_zeta_obc(:,:,:,     &
     &                                                       Tindex),   &
     &                            checksum = Fhash)
# else
     &                            BOUNDARY(ng) % tl_zeta_obc(:,:,:,     &
     &                                                       Tindex))
# endif
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              IF (Master) THEN
                WRITE (stdout,60) string, TRIM(Vname(1,ifield)),        &
     &                            InpRec, TRIM(ncname)
              END IF
              exit_flag=2
              ioerror=status
              RETURN
            ELSE
              IF (Master) THEN
# ifdef CHECKSUM
                WRITE (stdout,70) TRIM(Vname(2,ifield)), Fmin, Fmax,    &
     &                            Fhash
# else
                WRITE (stdout,70) TRIM(Vname(2,ifield)), Fmin, Fmax
# endif
              END IF
            END IF
          ELSE
            IF (Master) THEN
              WRITE (stdout,80) string, TRIM(Vname(1,ifield)),          &
     &                          TRIM(ncname)
            END IF
            exit_flag=4
            IF (FoundError(exit_flag, nf90_noerr,                       &
     &                     __LINE__, MyFile)) THEN
              RETURN
            END IF
          END IF
        END IF
!
!  Read in 2D U-momentum component open boundaries adjustments.
!
        IF (get_var(idSbry(isUbar)).and.get_adjust.and.                 &
     &        ANY(Lobc(:,isUbar,ng))) THEN
          ifield=idSbry(isUbar)
          foundit=find_string(var_name, n_var, TRIM(Vname(1,ifield)),   &
     &                        varid)
          IF (foundit) THEN
            status=nf_fread2d_bry(ng, IDmod, ncname, ncINPid,           &
     &                            Vname(1,ifield), varid,               &
     &                            InpRec, u2dvar,                       &
     &                            LBij, UBij, Nbrec(ng),                &
     &                            Fscl, Fmin, Fmax,                     &
# ifdef CHECKSUM
     &                            BOUNDARY(ng) % tl_ubar_obc(:,:,:,     &
     &                                                       Tindex),   &
     &                            checksum = Fhash)
# else
     &                            BOUNDARY(ng) % tl_ubar_obc(:,:,:,     &
     &                                                       Tindex))
# endif
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              IF (Master) THEN
                WRITE (stdout,60) string, TRIM(Vname(1,ifield)),        &
     &                            InpRec, TRIM(ncname)
              END IF
              exit_flag=2
              ioerror=status
              RETURN
            ELSE
              IF (Master) THEN
# ifdef CHECKSUM
                WRITE (stdout,70) TRIM(Vname(2,ifield)), Fmin, Fmax,    &
     &                            Fhash
# else
                WRITE (stdout,70) TRIM(Vname(2,ifield)), Fmin, Fmax
# endif
              END IF
            END IF
          ELSE
            IF (Master) THEN
              WRITE (stdout,80) string, TRIM(Vname(1,ifield)),          &
     &                          TRIM(ncname)
            END IF
            exit_flag=4
            IF (FoundError(exit_flag, nf90_noerr,                       &
     &                     __LINE__, MyFile)) THEN
              RETURN
            END IF
          END IF
        END IF
!
!  Read in 2D V-momentum component open boundaries adjustments.
!
        IF (get_var(idSbry(isVbar)).and.get_adjust.and.                 &
     &      ANY(Lobc(:,isVbar,ng))) THEN
          ifield=idSbry(isVbar)
          foundit=find_string(var_name, n_var, TRIM(Vname(1,ifield)),   &
     &                        varid)
          IF (foundit) THEN
            status=nf_fread2d_bry(ng, IDmod, ncname, ncINPid,           &
     &                            Vname(1,ifield), varid,               &
     &                            InpRec, v2dvar,                       &
     &                            LBij, UBij, Nbrec(ng),                &
     &                            Fscl, Fmin, Fmax,                     &
# ifdef CHECKSUM
     &                            BOUNDARY(ng) % tl_vbar_obc(:,:,:,     &
     &                                                       Tindex),   &
     &                            checksum = Fhash)
# else
     &                            BOUNDARY(ng) % tl_vbar_obc(:,:,:,     &
     &                                                       Tindex))
# endif
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              IF (Master) THEN
                WRITE (stdout,60) string, TRIM(Vname(1,ifield)),        &
     &                            InpRec, TRIM(ncname)
              END IF
              exit_flag=2
              ioerror=status
              RETURN
            ELSE
              IF (Master) THEN
# ifdef CHECKSUM
                WRITE (stdout,70) TRIM(Vname(2,ifield)), Fmin, Fmax,    &
     &                            Fhash
# else
                WRITE (stdout,70) TRIM(Vname(2,ifield)), Fmin, Fmax
# endif
              END IF
            END IF
          ELSE
            IF (Master) THEN
              WRITE (stdout,80) string, TRIM(Vname(1,ifield)),          &
     &                          TRIM(ncname)
            END IF
            exit_flag=4
            IF (FoundError(exit_flag, nf90_noerr,                       &
     &                     __LINE__, MyFile)) THEN
              RETURN
            END IF
          END IF
        END IF

# ifdef SOLVE3D
!
!  Read in 3D U-momentum component open boundaries adjustments.
!
        IF (get_var(idSbry(isUvel)).and.get_adjust.and.                 &
     &      ANY(Lobc(:,isUvel,ng))) THEN
          ifield=idSbry(isUvel)
          foundit=find_string(var_name, n_var, TRIM(Vname(1,ifield)),   &
     &                        varid)
          IF (foundit) THEN
            status=nf_fread3d_bry(ng, IDmod, ncname, ncINPid,           &
     &                            Vname(1,ifield), varid,               &
     &                            InpRec, u3dvar,                       &
     &                            LBij, UBij, 1, N(ng), Nbrec(ng),      &
     &                            Fscl, Fmin, Fmax,                     &
#  ifdef CHECKSUM
     &                            BOUNDARY(ng) % tl_u_obc(:,:,:,:,      &
     &                                                    Tindex),      &
     &                            checksum = Fhash)
#  else
     &                            BOUNDARY(ng) % tl_u_obc(:,:,:,:,      &
     &                                                    Tindex))
#  endif
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              IF (Master) THEN
                WRITE (stdout,60) string, TRIM(Vname(1,ifield)),        &
     &                            InpRec, TRIM(ncname)
              END IF
              exit_flag=2
              ioerror=status
              RETURN
            ELSE
              IF (Master) THEN
#  ifdef CHECKSUM
                WRITE (stdout,70) TRIM(Vname(2,ifield)), Fmin, Fmax,    &
     &                            Fhash
#  else
                WRITE (stdout,70) TRIM(Vname(2,ifield)), Fmin, Fmax
#  endif
              END IF
            END IF
          ELSE
            IF (Master) THEN
              WRITE (stdout,80) string, TRIM(Vname(1,ifield)),          &
     &                          TRIM(ncname)
            END IF
            exit_flag=4
            IF (FoundError(exit_flag, nf90_noerr,                       &
     &                     __LINE__, MyFile)) THEN
              RETURN
            END IF
          END IF
        END IF
!
!  Read in 3D V-momentum component open boundaries adjustments.
!
        IF (get_var(idSbry(isVvel)).and.get_adjust.and.                 &
     &      ANY(Lobc(:,isVvel,ng))) THEN
          ifield=idSbry(isVvel)
          foundit=find_string(var_name, n_var, TRIM(Vname(1,ifield)),   &
     &                        varid)
          IF (foundit) THEN
            status=nf_fread3d_bry(ng, IDmod, ncname, ncINPid,           &
     &                            Vname(1,ifield), varid,               &
     &                            InpRec, v3dvar,                       &
     &                            LBij, UBij, 1, N(ng), Nbrec(ng),      &
     &                            Fscl, Fmin, Fmax,                     &
#  ifdef CHECKSUM
     &                            BOUNDARY(ng) % tl_v_obc(:,:,:,:,      &
     &                                                    Tindex),      &
     &                            checksum = Fhash)
#  else
     &                            BOUNDARY(ng) % tl_v_obc(:,:,:,:,      &
     &                                                    Tindex))
#  endif
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              IF (Master) THEN
                WRITE (stdout,60) string, TRIM(Vname(1,ifield)),        &
     &                            InpRec, TRIM(ncname)
              END IF
              exit_flag=2
              ioerror=status
              RETURN
            ELSE
              IF (Master) THEN
#  ifdef CHECKSUM
                WRITE (stdout,70) TRIM(Vname(2,ifield)), Fmin, Fmax,    &
     &                            Fhash
#  else
                WRITE (stdout,70) TRIM(Vname(2,ifield)), Fmin, Fmax
#  endif

              END IF
            END IF
          ELSE
            IF (Master) THEN
              WRITE (stdout,80) string, TRIM(Vname(1,ifield)),          &
     &                          TRIM(ncname)
            END IF
            exit_flag=4
            IF (FoundError(exit_flag, nf90_noerr,                       &
     &                     __LINE__, MyFile)) THEN
              RETURN
            END IF
          END IF
        END IF
!
!  Read in 3D tracers open boundaries adjustments.
!
        DO itrc=1,NT(ng)
          IF (get_var(idSbry(isTvar(itrc))).and.get_adjust.and.         &
     &        ANY(Lobc(:,isTvar(itrc),ng))) THEN
            ifield=idSbry(isTvar(itrc))
            foundit=find_string(var_name, n_var, TRIM(Vname(1,ifield)), &
     &                          varid)
            IF (foundit) THEN
              status=nf_fread3d_bry(ng, IDmod, ncname, ncINPid,         &
     &                              Vname(1,ifield), varid,             &
     &                              InpRec, r3dvar,                     &
     &                              LBij, UBij, 1, N(ng), Nbrec(ng),    &
     &                              Fscl, Fmin, Fmax,                   &
#  ifdef CHECKSUM
     &                              BOUNDARY(ng) % tl_t_obc(:,:,:,:,    &
     &                                                   Tindex,itrc),  &
     &                              checksum = Fhash)
#  else
     &                              BOUNDARY(ng) % tl_t_obc(:,:,:,:,    &
     &                                                   Tindex,itrc))
#  endif
              IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
                IF (Master) THEN
                  WRITE (stdout,60) string, TRIM(Vname(1,ifield)),      &
     &                              InpRec, TRIM(ncname)
                END IF
                exit_flag=2
                ioerror=status
                RETURN
              ELSE
                IF (Master) THEN
#  ifdef CHECKSUM
                  WRITE (stdout,70) TRIM(Vname(2,ifield)), Fmin, Fmax,  &
     &                            Fhash
#  else
                  WRITE (stdout,70) TRIM(Vname(2,ifield)), Fmin, Fmax
#  endif
                END IF
              END IF
            ELSE
              IF (Master) THEN
                WRITE (stdout,80) string, TRIM(Vname(1,ifield)),        &
     &                            TRIM(ncname)
              END IF
              exit_flag=4
              IF (FoundError(exit_flag, nf90_noerr,                     &
     &                       __LINE__, MyFile)) THEN
                RETURN
              END IF
            END IF
          END IF
        END DO
# endif
#endif
#ifdef ADJUST_WSTRESS
!
!  Read in tangent linear surface U-momentum stress.
!
        IF (get_var(idUsms).and.get_adjust) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idUsms)),   &
     &                        varid)
          IF (foundit) THEN
            gtype=var_flag(varid)*u3dvar
            scale=1.0_dp
            status=nf_fread3d(ng, IDmod, ncname, ncINPid,               &
     &                        Vname(1,idUsms), varid,                   &
     &                        InpRec, gtype, Vsize,                     &
     &                        LBi, UBi, LBj, UBj, 1, Nfrec(ng),         &
     &                        scale, Fmin, Fmax,                        &
# ifdef MASKING
     &                        GRID(ng) % umask,                         &
# endif
# ifdef CHECKSUM
     &                        FORCES(ng) % tl_ustr(:,:,:,Tindex),       &
     &                        checksum = Fhash)
# else
     &                        FORCES(ng) % tl_ustr(:,:,:,Tindex))
# endif
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              IF (Master) THEN
                WRITE (stdout,60) string, TRIM(Vname(1,idUsms)),        &
     &                            InpRec, TRIM(ncname)
              END IF
              exit_flag=2
              ioerror=status
              RETURN
            ELSE
              IF (Master) THEN
# ifdef CHECKSUM
                WRITE (stdout,70) TRIM(Vname(2,idUsms))//               &
     &                            ', adjusted tl_ustr', Fmin, Fmax,     &
     &                            Fhash
# else
                WRITE (stdout,70) TRIM(Vname(2,idUsms))//               &
     &                            ', adjusted tl_ustr', Fmin, Fmax
# endif
              END IF
            END IF
          ELSE
            IF (Master) THEN
              WRITE (stdout,80) string, TRIM(Vname(1,idUsms)),          &
     &                          TRIM(ncname)
            END IF
            exit_flag=4
            IF (FoundError(exit_flag, nf90_noerr,                       &
     &                     __LINE__, MyFile)) THEN
              RETURN
            END IF
          END IF
        END IF
!
!  Read in tangent linear surface V-momentum stress.
!
        IF (get_var(idVsms).and.get_adjust) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idVsms)),   &
     &                        varid)
          IF (foundit) THEN
            gtype=var_flag(varid)*v3dvar
            scale=1.0_dp
            status=nf_fread3d(ng, IDmod, ncname, ncINPid,               &
     &                        Vname(1,idVsms), varid,                   &
     &                        InpRec, gtype, Vsize,                     &
     &                        LBi, UBi, LBj, UBj, 1, Nfrec(ng),         &
     &                        scale, Fmin, Fmax,                        &
# ifdef MASKING
     &                        GRID(ng) % vmask,                         &
# endif
# ifdef CHECKSUM
     &                        FORCES(ng) % tl_vstr(:,:,:,Tindex),       &
     &                        checksum = Fhash)
# else
     &                        FORCES(ng) % tl_vstr(:,:,:,Tindex))
# endif
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              IF (Master) THEN
                WRITE (stdout,60) string, TRIM(Vname(1,idVsms)),        &
     &                            InpRec, TRIM(ncname)
              END IF
              exit_flag=2
              ioerror=status
              RETURN
            ELSE
              IF (Master) THEN
# ifdef CHECKSUM
                WRITE (stdout,70) TRIM(Vname(2,idVsms))//               &
     &                            ', adjusted tl_vstr', Fmin, Fmax,     &
     &                            Fhash
# else
                WRITE (stdout,70) TRIM(Vname(2,idVsms))//               &
     &                            ', adjusted tl_vstr', Fmin, Fmax
# endif
              END IF
            END IF
          ELSE
            IF (Master) THEN
              WRITE (stdout,80) string, TRIM(Vname(1,idVsms)),          &
     &                          TRIM(ncname)
            END IF
            exit_flag=4
            IF (FoundError(exit_flag, nf90_noerr,                       &
     &                     __LINE__, MyFile)) THEN
              RETURN
            END IF
          END IF
        END IF
#endif
#if defined ADJUST_STFLUX && defined SOLVE3D
!
!  Read in tangent linear surface tracers flux.
!
        DO itrc=1,NT(ng)
          IF (get_var(idTsur(itrc)).and.get_adjust.and.                 &
     &        Lstflux(itrc,ng)) THEN
            foundit=find_string(var_name, n_var,                        &
     &                          TRIM(Vname(1,idTsur(itrc))), varid)
            IF (foundit) THEN
              gtype=var_flag(varid)*r3dvar
              scale=1.0_dp
              status=nf_fread3d(ng, IDmod, ncname, ncINPid,             &
     &                          Vname(1,idTsur(itrc)), varid,           &
     &                          InpRec, gtype, Vsize,                   &
     &                          LBi, UBi, LBj, UBj, 1, Nfrec(ng),       &
     &                          scale, Fmin, Fmax,                      &
# ifdef MASKING
     &                          GRID(ng) % rmask,                       &
# endif
# ifdef CHECKSUM
     &                          FORCES(ng)% tl_tflux(:,:,:,             &
     &                                               Tindex,itrc),      &
     &                          checksum = Fhash)
# else
     &                          FORCES(ng)% tl_tflux(:,:,:,             &
     &                                               Tindex,itrc))
# endif
              IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
                IF (Master) THEN
                  WRITE (stdout,60) string, TRIM(Vname(1,idTsur(itrc))),&
     &                              InpRec, TRIM(ncname)
                END IF
                exit_flag=2
                ioerror=status
                RETURN
              ELSE
                IF (Master) THEN
# ifdef CHECKSUM
                  WRITE (stdout,70) TRIM(Vname(2,idTsur(itrc)))//       &
     &                              ', adjusted tl_tflux', Fmin, Fmax,  &
     &                              Fhash
# else
                  WRITE (stdout,70) TRIM(Vname(2,idTsur(itrc)))//       &
     &                              ', adjusted tl_tflux', Fmin, Fmax
# endif
                END IF
              END IF
            ELSE
              IF (Master) THEN
                WRITE (stdout,80) string, TRIM(Vname(1,idTsur(itrc))),  &
     &                            TRIM(ncname)
              END IF
              exit_flag=4
              IF (FoundError(exit_flag, nf90_noerr,                     &
     &                       __LINE__, MyFile)) THEN
                RETURN
              END IF
            END IF
          END IF
        END DO
#endif
