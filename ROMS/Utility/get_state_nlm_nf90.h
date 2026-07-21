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
** It process the reading of ROMS NLM state from input NetCDF files.  **
**                                                                    **
** It is used to read nonlinear kernel initial conditions, restart,   **
** and trajectory data.                                               **
**                                                                    **
************************************************************************
*/

#ifdef PERFECT_RESTART
!
!  Read in time-stepping indices.
!
        IF ((model.eq.0).and.(nrrec(ng).ne.0)) THEN
# ifdef SOLVE3D
          CALL netcdf_get_ivar (ng, IDmod, ncname,                      &
     &                          'nstp', nstp(ng:),                      &
     &                          ncid = ncINPid,                         &
     &                          start = (/InpRec/),                     &
     &                          total = (/1/))
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

          CALL netcdf_get_ivar (ng, IDmod, ncname,                      &
     &                          'nrhs', nrhs(ng:),                      &
     &                          ncid = ncINPid,                         &
     &                          start = (/InpRec/),                     &
     &                          total = (/1/))
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

          CALL netcdf_get_ivar (ng, IDmod, ncname,                      &
     &                          'nnew', nnew(ng:),                      &
     &                          ncid = ncINPid,                         &
     &                          start = (/InpRec/),                     &
     &                          total = (/1/))
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif
          CALL netcdf_get_ivar (ng, IDmod, ncname,                      &
     &                          'kstp', kstp(ng:),                      &
     &                          ncid = ncINPid,                         &
     &                          start = (/InpRec/),                     &
     &                          total = (/1/))
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

          CALL netcdf_get_ivar (ng, IDmod, ncname,                      &
     &                          'krhs', krhs(ng:),                      &
     &                          ncid = ncINPid,                         &
     &                          start = (/InpRec/),                     &
     &                          total = (/1/))
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

          CALL netcdf_get_ivar (ng, IDmod, ncname,                      &
     &                          'knew', knew(ng:),                      &
     &                          ncid = ncINPid,                         &
     &                          start = (/InpRec/),                     &
     &                          total = (/1/))
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#endif
#if defined SEDIMENT && defined SED_MORPH
!
!  Read in time-evolving bathymetry (m).
!
        IF (get_var(idbath)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idbath)),   &
     &                        varid)
          IF (foundit) THEN
            gtype=var_flag(varid)*r2dvar
            status=nf_fread2d(ng, IDmod, ncname, ncINPid,               &
     &                        Vname(1,idbath), varid,                   &
     &                        InpRec, gtype, Vsize,                     &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
# ifdef MASKING
     &                        GRID(ng) % rmask,                         &
# endif
# ifdef CHECKSUM
     &                        GRID(ng) % h,                             &
     &                        checksum = Fhash)
# else
     &                        GRID(ng) % h)
# endif
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              IF (Master) THEN
                WRITE (stdout,60) string, TRIM(Vname(1,idbath)),        &
     &                            InpRec, TRIM(ncname)
              END IF
              exit_flag=2
              ioerror=status
              RETURN
            ELSE
              IF (Master) THEN
# ifdef CHECKSUM
                WRITE (stdout,70) TRIM(Vname(2,idbath)), Fmin, Fmax,    &
     &                            Fhash
# else
                WRITE (stdout,70) TRIM(Vname(2,idbath)), Fmin, Fmax
# endif

              END IF
            END IF
          ELSE
            IF (Master) THEN
              WRITE (stdout,80) string, TRIM(Vname(1,idbath)),          &
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
!
!  Read in nonlinear free-surface (m).
!
        IF (get_var(idFsur)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idFsur)),   &
     &                        varid)
          IF (foundit) THEN
            IF (Perfect2D) THEN
              gtype=var_flag(varid)*r3dvar
            ELSE
              gtype=var_flag(varid)*r2dvar
            END IF
            IF (Perfect2D) THEN
              status=nf_fread3d(ng, IDmod, ncname, ncINPid,             &
     &                          Vname(1,idFsur), varid,                 &
     &                          InpRec, gtype, Vsize,                   &
     &                          LBi, UBi, LBj, UBj, 1, 3,               &
     &                          Fscl, Fmin, Fmax,                       &
#ifdef MASKING
     &                          GRID(ng) % rmask,                       &
#endif
#ifdef CHECKSUM
     &                          OCEAN(ng) % zeta,                       &
     &                          checksum = Fhash)
#else
     &                          OCEAN(ng) % zeta)
#endif
            ELSE
              status=nf_fread2d(ng, IDmod, ncname, ncINPid,             &
     &                          Vname(1,idFsur), varid,                 &
     &                          InpRec, gtype, Vsize,                   &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Fscl, Fmin, Fmax,                       &
#ifdef MASKING
     &                          GRID(ng) % rmask,                       &
#endif
#ifdef CHECKSUM
     &                          OCEAN(ng) % zeta(:,:,Tindex),           &
     &                          checksum = Fhash)
#else
     &                          OCEAN(ng) % zeta(:,:,Tindex))
#endif
            END IF
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              IF (Master) THEN
                WRITE (stdout,60) string, TRIM(Vname(1,idFsur)),        &
     &                            InpRec, TRIM(ncname)
              END IF
              exit_flag=2
              ioerror=status
              RETURN
            ELSE
              IF (Master) THEN
#ifdef CHECKSUM
                WRITE (stdout,70) TRIM(Vname(2,idFsur)), Fmin, Fmax,    &
     &                            Fhash
#else
                WRITE (stdout,70) TRIM(Vname(2,idFsur)), Fmin, Fmax
#endif

              END IF
            END IF
          ELSE
            IF (Master) THEN
              WRITE (stdout,80) string, TRIM(Vname(1,idFsur)),          &
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
!  Read in nonlinear RHS of free-surface.
!
        IF (get_var(idRzet).and.Perfect2D) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idRzet)),   &
     &                        varid)
          IF (foundit) THEN
            gtype=var_flag(varid)*r3dvar
            status=nf_fread3d(ng, IDmod, ncname, ncINPid,               &
     &                        Vname(1,idRzet), varid,                   &
     &                        InpRec, gtype, Vsize,                     &
     &                        LBi, UBi, LBj, UBj, 1, 2,                 &
     &                        Fscl, Fmin, Fmax,                         &
#ifdef MASKING
     &                        GRID(ng) % rmask,                         &
#endif
#ifdef CHECKSUM
     &                        OCEAN(ng) % rzeta,                        &
     &                        checksum = Fhash)
#else
     &                        OCEAN(ng) % rzeta)
#endif
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              IF (Master) THEN
                WRITE (stdout,60) string, TRIM(Vname(1,idRzet)),        &
     &                            InpRec, TRIM(ncname)
              END IF
              exit_flag=2
              ioerror=status
              RETURN
            ELSE
              IF (Master) THEN
#ifdef CHECKSUM
                WRITE (stdout,70) TRIM(Vname(2,idRzet)), Fmin, Fmax,    &
     &                            Fhash
#else
                WRITE (stdout,70) TRIM(Vname(2,idRzet)), Fmin, Fmax
#endif

              END IF
            END IF
          ELSE
            IF (Master) THEN
              WRITE (stdout,80) string, TRIM(Vname(1,idRzet)),          &
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
!  Read in nonlinear 2D U-momentum component (m/s).
!
        IF (get_var(idUbar)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idUbar)),   &
     &                        varid)
          IF (foundit) THEN
            IF (Perfect2D) THEN
              gtype=var_flag(varid)*u3dvar
            ELSE
              gtype=var_flag(varid)*u2dvar
            END IF
            IF (Perfect2D) THEN
              status=nf_fread3d(ng, IDmod, ncname, ncINPid,             &
     &                          Vname(1,idUbar), varid,                 &
     &                          InpRec, gtype, Vsize,                   &
     &                          LBi, UBi, LBj, UBj, 1, 3,               &
     &                          Fscl, Fmin, Fmax,                       &
#ifdef MASKING
     &                          GRID(ng) % umask,                       &
#endif
#ifdef CHECKSUM
     &                          OCEAN(ng) % ubar,                       &
     &                          checksum = Fhash)
#else
     &                          OCEAN(ng) % ubar)
#endif
            ELSE
              status=nf_fread2d(ng, IDmod, ncname, ncINPid,             &
     &                          Vname(1,idUbar), varid,                 &
     &                          InpRec, gtype, Vsize,                   &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Fscl, Fmin, Fmax,                       &
#ifdef MASKING
     &                          GRID(ng) % umask,                       &
#endif
#ifdef CHECKSUM
     &                          OCEAN(ng) % ubar(:,:,Tindex),           &
     &                          checksum = Fhash)
#else
     &                          OCEAN(ng) % ubar(:,:,Tindex))
#endif
            END IF
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              IF (Master) THEN
                WRITE (stdout,60) string, TRIM(Vname(1,idUbar)),        &
     &                            InpRec, TRIM(ncname)
              END IF
              exit_flag=2
              ioerror=status
              RETURN
            ELSE
              IF (Master) THEN
#ifdef CHECKSUM
                WRITE (stdout,70) TRIM(Vname(2,idUbar)), Fmin, Fmax,    &
     &                            Fhash
#else
                WRITE (stdout,70) TRIM(Vname(2,idUbar)), Fmin, Fmax
#endif

              END IF
            END IF
          ELSE
            IF (Master) THEN
              WRITE (stdout,80) string, TRIM(Vname(1,idUbar)),          &
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
!  Read in nonlinear RHS of 2D U-momentum component.
!
        IF (get_var(idRu2d).and.Perfect2D) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idRu2d)),   &
     &                        varid)
          IF (foundit) THEN
            gtype=var_flag(varid)*u3dvar
            status=nf_fread3d(ng, IDmod, ncname, ncINPid,               &
     &                        Vname(1,idRu2d), varid,                   &
     &                        InpRec, gtype, Vsize,                     &
     &                        LBi, UBi, LBj, UBj, 1, 2,                 &
     &                        Fscl, Fmin, Fmax,                         &
#ifdef MASKING
     &                        GRID(ng) % umask,                         &
#endif
#ifdef CHECKSUM
     &                        OCEAN(ng) % rubar,                        &
     &                        checksum = Fhash)
#else
     &                        OCEAN(ng) % rubar)
#endif
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              IF (Master) THEN
                WRITE (stdout,60) string, TRIM(Vname(1,idRu2d)),        &
     &                            InpRec, TRIM(ncname)
              END IF
              exit_flag=2
              ioerror=status
              RETURN
            ELSE
              IF (Master) THEN
#ifdef CHECKSUM
                WRITE (stdout,70) TRIM(Vname(2,idRu2d)), Fmin, Fmax,    &
     &                            Fhash
#else
                WRITE (stdout,70) TRIM(Vname(2,idRu2d)), Fmin, Fmax
#endif

              END IF
            END IF
          ELSE
            IF (Master) THEN
              WRITE (stdout,80) string, TRIM(Vname(1,idRu2d)),          &
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
!  Read in nonlinear 2D V-momentum component (m/s).
!
        IF (get_var(idVbar)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idVbar)),   &
     &                        varid)
          IF (foundit) THEN
            IF (Perfect2D) THEN
              gtype=var_flag(varid)*v3dvar
            ELSE
              gtype=var_flag(varid)*v2dvar
            END IF
            IF (Perfect2D) THEN
              status=nf_fread3d(ng, IDmod, ncname, ncINPid,             &
     &                          Vname(1,idVbar), varid,                 &
     &                          InpRec, gtype, Vsize,                   &
     &                          LBi, UBi, LBj, UBj, 1, 3,               &
     &                          Fscl, Fmin, Fmax,                       &
#ifdef MASKING
     &                          GRID(ng) % vmask,                       &
#endif
#ifdef CHECKSUM
     &                          OCEAN(ng) % vbar,                       &
     &                          checksum = Fhash)
#else
     &                          OCEAN(ng) % vbar)
#endif
            ELSE
              status=nf_fread2d(ng, IDmod, ncname, ncINPid,             &
     &                          Vname(1,idVbar), varid,                 &
     &                          InpRec, gtype, Vsize,                   &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Fscl, Fmin, Fmax,                       &
#ifdef MASKING
     &                          GRID(ng) % vmask,                       &
#endif
#ifdef CHECKSUM
     &                          OCEAN(ng) % vbar(:,:,Tindex),           &
     &                          checksum = Fhash)
#else
     &                          OCEAN(ng) % vbar(:,:,Tindex))
#endif
            END IF
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              IF (Master) THEN
                WRITE (stdout,60) string, TRIM(Vname(1,idVbar)),        &
     &                            InpRec, TRIM(ncname)
              END IF
              exit_flag=2
              ioerror=status
              RETURN
            ELSE
              IF (Master) THEN
#ifdef CHECKSUM
                WRITE (stdout,70) TRIM(Vname(2,idVbar)), Fmin, Fmax,    &
     &                            Fhash
#else
                WRITE (stdout,70) TRIM(Vname(2,idVbar)), Fmin, Fmax
#endif

              END IF
            END IF
          ELSE
            IF (Master) THEN
              WRITE (stdout,80) string, TRIM(Vname(1,idVbar)),          &
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
!  Read in nonlinear RHS 2D V-momentum component.
!
        IF (get_var(idRv2d).and.Perfect2D) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idRv2d)),   &
     &                        varid)
          IF (foundit) THEN
            gtype=var_flag(varid)*v3dvar
            status=nf_fread3d(ng, IDmod, ncname, ncINPid,               &
     &                        Vname(1,idRv2d), varid,                   &
     &                        InpRec, gtype, Vsize,                     &
     &                        LBi, UBi, LBj, UBj, 1, 2,                 &
     &                        Fscl, Fmin, Fmax,                         &
#ifdef MASKING
     &                        GRID(ng) % vmask,                         &
#endif
#ifdef CHECKSUM
     &                        OCEAN(ng) % rvbar,                        &
     &                        checksum = Fhash)
#else
     &                        OCEAN(ng) % rvbar)
#endif
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              IF (Master) THEN
                WRITE (stdout,60) string, TRIM(Vname(1,idRv2d)),        &
     &                            InpRec, TRIM(ncname)
              END IF
              exit_flag=2
              ioerror=status
              RETURN
            ELSE
              IF (Master) THEN
#ifdef CHECKSUM
                WRITE (stdout,70) TRIM(Vname(2,idRv2d)), Fmin, Fmax,    &
     &                            Fhash
#else
                WRITE (stdout,70) TRIM(Vname(2,idRv2d)), Fmin, Fmax
#endif

              END IF
            END IF
          ELSE
            IF (Master) THEN
              WRITE (stdout,80) string, TRIM(Vname(1,idRv2d)),          &
     &                          TRIM(ncname)
            END IF
            exit_flag=4
            IF (FoundError(exit_flag, nf90_noerr,                       &
     &                     __LINE__, MyFile)) THEN
              RETURN
            END IF
          END IF
        END IF

#ifdef SOLVE3D
!
!  Read in nonlinear 3D U-momentum component (m/s).
!
        IF (get_var(idUvel)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idUvel)),   &
     &                        varid)
          IF (foundit) THEN
            gtype=var_flag(varid)*u3dvar
            IF (Perfect3D) THEN
              status=nf_fread4d(ng, IDmod, ncname, ncINPid,             &
     &                          Vname(1,idUvel), varid,                 &
     &                          InpRec, gtype, Vsize,                   &
     &                          LBi, UBi, LBj, UBj, 1, N(ng), 1, 2,     &
     &                          Fscl, Fmin, Fmax,                       &
# ifdef MASKING
     &                          GRID(ng) % umask,                       &
# endif
# ifdef CHECKSUM
     &                          OCEAN(ng) % u,                          &
     &                          checksum = Fhash)
# else
     &                          OCEAN(ng) % u)
# endif
            ELSE
              status=nf_fread3d(ng, IDmod, ncname, ncINPid,             &
     &                          Vname(1,idUvel), varid,                 &
     &                          InpRec, gtype, Vsize,                   &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          Fscl, Fmin, Fmax,                       &
# ifdef MASKING
     &                          GRID(ng) % umask,                       &
# endif
# ifdef CHECKSUM
     &                          OCEAN(ng) % u(:,:,:,Tindex),            &
     &                          checksum = Fhash)
# else
     &                          OCEAN(ng) % u(:,:,:,Tindex))
# endif
            END IF
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              IF (Master) THEN
                WRITE (stdout,60) string, TRIM(Vname(1,idUvel)),        &
     &                            InpRec, TRIM(ncname)
              END IF
              exit_flag=2
              ioerror=status
              RETURN
            ELSE
              IF (Master) THEN
# ifdef CHECKSUM
                WRITE (stdout,70) TRIM(Vname(2,idUvel)), Fmin, Fmax,    &
     &                            Fhash
# else
                WRITE (stdout,70) TRIM(Vname(2,idUvel)), Fmin, Fmax
# endif
              END IF
            END IF
          ELSE
            IF (Master) THEN
              WRITE (stdout,80) string, TRIM(Vname(1,idUvel)),          &
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
!  Read in nonlinear RHS of 3D U-momentum component.
!
        IF (get_var(idRu3d).and.Perfect3D) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idRu3d)),   &
     &                        varid)
          IF (foundit) THEN
            gtype=var_flag(varid)*u3dvar
            status=nf_fread4d(ng, IDmod, ncname, ncINPid,               &
     &                        Vname(1,idRu3d), varid,                   &
     &                        InpRec, gtype, Vsize,                     &
     &                        LBi, UBi, LBj, UBj, 0, N(ng), 1, 2,       &
     &                        Fscl, Fmin, Fmax,                         &
# ifdef MASKING
     &                        GRID(ng) % umask,                         &
# endif
# ifdef CHECKSUM
     &                        OCEAN(ng) % ru,                           &
     &                        checksum = Fhash)
# else
     &                        OCEAN(ng) % ru)
# endif
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              IF (Master) THEN
                WRITE (stdout,60) string, TRIM(Vname(1,idRu3d)),        &
     &                            InpRec, TRIM(ncname)
              END IF
              exit_flag=2
              ioerror=status
              RETURN
            ELSE
              IF (Master) THEN
# ifdef CHECKSUM
                WRITE (stdout,70) TRIM(Vname(2,idRu3d)), Fmin, Fmax,    &
     &                            Fhash
# else
                WRITE (stdout,70) TRIM(Vname(2,idRu3d)), Fmin, Fmax
# endif
              END IF
            END IF
          ELSE
            IF (Master) THEN
              WRITE (stdout,80) string, TRIM(Vname(1,idRu3d)),          &
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
!  Read in nonlinear 3D V-momentum component (m/s).
!
        IF (get_var(idVvel)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idVvel)),   &
     &                        varid)
          IF (foundit) THEN
            gtype=var_flag(varid)*v3dvar
            IF (Perfect3D) THEN
              status=nf_fread4d(ng, IDmod, ncname, ncINPid,             &
     &                          Vname(1,idVvel), varid,                 &
     &                          InpRec, gtype, Vsize,                   &
     &                          LBi, UBi, LBj, UBj, 1, N(ng), 1, 2,     &
     &                          Fscl, Fmin, Fmax,                       &
# ifdef MASKING
     &                          GRID(ng) % vmask,                       &
# endif
# ifdef CHECKSUM
     &                          OCEAN(ng) % v,                          &
     &                          checksum = Fhash)
# else
     &                          OCEAN(ng) % v)
# endif
            ELSE
              status=nf_fread3d(ng, IDmod, ncname, ncINPid,             &
     &                          Vname(1,idVvel), varid,                 &
     &                          InpRec, gtype, Vsize,                   &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          Fscl, Fmin, Fmax,                       &
# ifdef MASKING
     &                          GRID(ng) % vmask,                       &
# endif
# ifdef CHECKSUM
     &                          OCEAN(ng) % v(:,:,:,Tindex),            &
     &                          checksum = Fhash)
# else
     &                          OCEAN(ng) % v(:,:,:,Tindex))
# endif
            END IF
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              IF (Master) THEN
                WRITE (stdout,60) string, TRIM(Vname(1,idVvel)),        &
     &                            InpRec, TRIM(ncname)
              END IF
              exit_flag=2
              ioerror=status
              RETURN
            ELSE
              IF (Master) THEN
# ifdef CHECKSUM
                WRITE (stdout,70) TRIM(Vname(2,idVvel)), Fmin, Fmax,    &
     &                            Fhash
# else
                WRITE (stdout,70) TRIM(Vname(2,idVvel)), Fmin, Fmax
# endif

              END IF
            END IF
          ELSE
            IF (Master) THEN
              WRITE (stdout,80) string, TRIM(Vname(1,idVvel)),          &
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
!  Read in nonlinear RHS of 3D V-momentum component.
!
        IF (get_var(idRv3d).and.Perfect3D) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idRv3d)),   &
     &                        varid)
          IF (foundit) THEN
            gtype=var_flag(varid)*v3dvar
            status=nf_fread4d(ng, IDmod, ncname, ncINPid,               &
     &                        Vname(1,idRv3d), varid,                   &
     &                        InpRec, gtype, Vsize,                     &
     &                        LBi, UBi, LBj, UBj, 0, N(ng), 1, 2,       &
     &                        Fscl, Fmin, Fmax,                         &
# ifdef MASKING
     &                        GRID(ng) % vmask,                         &
# endif
# ifdef CHECKSUM
     &                        OCEAN(ng) % rv,                           &
     &                        checksum = Fhash)
# else
     &                        OCEAN(ng) % rv)
# endif
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              IF (Master) THEN
                WRITE (stdout,60) string, TRIM(Vname(1,idRv3d)),        &
     &                            InpRec, TRIM(ncname)
              END IF
              exit_flag=2
              ioerror=status
              RETURN
            ELSE
              IF (Master) THEN
# ifdef CHECKSUM
                WRITE (stdout,70) TRIM(Vname(2,idRv3d)), Fmin, Fmax,    &
     &                            Fhash
# else
                WRITE (stdout,70) TRIM(Vname(2,idRv3d)), Fmin, Fmax
# endif

              END IF
            END IF
          ELSE
            IF (Master) THEN
              WRITE (stdout,80) string, TRIM(Vname(1,idRv3d)),          &
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
!  Read in nonlinear tracer type variables.
!
        DO itrc=1,NT(ng)
          IF (get_var(idTvar(itrc))) THEN
            foundit=find_string(var_name, n_var,                        &
     &                          TRIM(Vname(1,idTvar(itrc))), varid)
            IF (foundit) THEN
              gtype=var_flag(varid)*r3dvar
              IF (Perfect3D) THEN
                status=nf_fread4d(ng, IDmod, ncname, ncINPid,           &
     &                            Vname(1,idTvar(itrc)), varid,         &
     &                            InpRec, gtype, Vsize,                 &
     &                            LBi, UBi, LBj, UBj, 1, N(ng), 1, 2,   &
     &                            Fscl, Fmin, Fmax,                     &
# ifdef MASKING
     &                            GRID(ng) % rmask,                     &
# endif
# ifdef CHECKSUM
     &                            OCEAN(ng) % t(:,:,:,:,itrc),          &
     &                            checksum = Fhash)
# else
     &                            OCEAN(ng) % t(:,:,:,:,itrc))
# endif
              ELSE
                status=nf_fread3d(ng, IDmod, ncname, ncINPid,           &
     &                            Vname(1,idTvar(itrc)), varid,         &
     &                            InpRec, gtype, Vsize,                 &
     &                            LBi, UBi, LBj, UBj, 1, N(ng),         &
     &                            Fscl, Fmin, Fmax,                     &
# ifdef MASKING
     &                            GRID(ng) % rmask,                     &
# endif
# ifdef CHECKSUM
     &                            OCEAN(ng) % t(:,:,:,Tindex,itrc),     &
     &                            checksum = Fhash)
# else
     &                            OCEAN(ng) % t(:,:,:,Tindex,itrc))
# endif
              END IF
              IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
                IF (Master) THEN
                  WRITE (stdout,60) string, TRIM(Vname(1,idTvar(itrc))),&
     &                              InpRec, TRIM(ncname)
                END IF
                exit_flag=2
                ioerror=status
                RETURN
              ELSE
                IF (Master) THEN
# ifdef CHECKSUM
                  WRITE (stdout,70) TRIM(Vname(2,idTvar(itrc))),        &
     &                              Fmin, Fmax, Fhash
# else
                  WRITE (stdout,70) TRIM(Vname(2,idTvar(itrc))),        &
     &                              Fmin, Fmax
# endif
                END IF
              END IF
            ELSE
              IF (Master) THEN
                WRITE (stdout,80) string, TRIM(Vname(1,idTvar(itrc))),  &
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
#  ifdef NEMURO_SED1
!
!  Read in bio sediment variables.
!
        IF (get_var(idPONsed)) THEN
          foundit=find_string(var_name, n_var,                          &
     &                          TRIM(Vname(1,idPONsed)), varid)
          gtype=var_flag(varid)*r2dvar
          status=nf_fread2d(ng, IDmod, ncname, ncINPid,                 &
     &                      Vname(1,idPONsed), varid,                   &
     &                      InpRec, gtype, Vsize,                       &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      Fscl, Fmin, Fmax,                           &
#   ifdef MASKING
     &                      GRID(ng) % rmask,                           &
#   endif
     &                      OCEAN(ng) % PONsed)
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,60) string, TRIM(Vname(1,idPONsed)),        &
     &                          InpRec, TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,70) TRIM(Vname(2,idPONsed)),               &
     &                           Fmin, Fmax
            END IF
          END IF
        END IF
!
        IF (get_var(idOPALsed)) THEN
          foundit=find_string(var_name, n_var,                          &
     &                          TRIM(Vname(1,idOPALsed)), varid)
          gtype=var_flag(varid)*r2dvar
          status=nf_fread2d(ng, IDmod, ncname, ncINPid,                 &
     &                      Vname(1,idOPALsed), varid,                  &
     &                      InpRec, gtype, Vsize,                       &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      Fscl, Fmin, Fmax,                           &
#   ifdef MASKING
     &                      GRID(ng) % rmask,                           &
#   endif
     &                      OCEAN(ng) % OPALsed)
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,60) string, TRIM(Vname(1,idOPALsed)),       &
     &                          InpRec, TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,70) TRIM(Vname(2,idOPALsed)),              &
     &                           Fmin, Fmax
            END IF
          END IF
        END IF
!
        IF (get_var(idDENITsed)) THEN
          foundit=find_string(var_name, n_var,                          &
     &                          TRIM(Vname(1,idDENITsed)), varid)
          gtype=var_flag(varid)*r2dvar
          status=nf_fread2d(ng, IDmod, ncname, ncINPid,                 &
     &                      Vname(1,idDENITsed), varid,                 &
     &                      InpRec, gtype, Vsize,                       &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      Fscl, Fmin, Fmax,                           &
#   ifdef MASKING
     &                      GRID(ng) % rmask,                           &
#   endif
     &                      OCEAN(ng) % DENITsed)
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,60) string, TRIM(Vname(1,idDENITsed)),      &
     &                          InpRec, TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,70) TRIM(Vname(2,idDENITsed)),             &
     &                           Fmin, Fmax
            END IF
          END IF
        END IF
!
        IF (get_var(idPONbur)) THEN
          foundit=find_string(var_name, n_var,                          &
     &                          TRIM(Vname(1,idPONbur)), varid)
          gtype=var_flag(varid)*r2dvar
          status=nf_fread2d(ng, IDmod, ncname, ncINPid,                 &
     &                      Vname(1,idPONbur), varid,                   &
     &                      InpRec, gtype, Vsize,                       &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      Fscl, Fmin, Fmax,                           &
#   ifdef MASKING
     &                      GRID(ng) % rmask,                           &
#   endif
     &                      OCEAN(ng) % PON_burial)
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,60) string, TRIM(Vname(1,idPONbur)),        &
     &                          InpRec, TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,70) TRIM(Vname(2,idPONbur)),               &
     &                           Fmin, Fmax
            END IF
          END IF
        END IF
!
        IF (get_var(idOPALbur)) THEN
          foundit=find_string(var_name, n_var,                          &
     &                          TRIM(Vname(1,idOPALbur)), varid)
          gtype=var_flag(varid)*r2dvar
          status=nf_fread2d(ng, IDmod, ncname, ncINPid,                 &
     &                      Vname(1,idOPALbur), varid,                  &
     &                      InpRec, gtype, Vsize,                       &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      Fscl, Fmin, Fmax,                           &
#   ifdef MASKING
     &                      GRID(ng) % rmask,                           &
#   endif
     &                      OCEAN(ng) % OPAL_burial)
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,60) string, TRIM(Vname(1,idOPALbur)),       &
     &                          InpRec, TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,70) TRIM(Vname(2,idOPALbur)),              &
     &                           Fmin, Fmax
            END IF
          END IF
        END IF
#  endif

# if defined GLS_MIXING || defined MY25_MIXING || defined LMD_MIXING
!
!  Read in vertical viscosity.
!
        IF (have_var(idVvis)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idVvis)),   &
     &                        varid)
          IF (foundit) THEN
            gtype=var_flag(varid)*w3dvar
            status=nf_fread3d(ng, IDmod, ncname, ncINPid,               &
     &                        Vname(1,idVvis), varid,                   &
     &                        InpRec, gtype, Vsize,                     &
     &                        LBi, UBi, LBj, UBj, 0, N(ng),             &
     &                        Fscl, Fmin,Fmax,                          &
#  ifdef MASKING
     &                        GRID(ng) % rmask,                         &
#  endif
#  ifdef CHECKSUM
     &                        MIXING(ng) % AKv,                         &
     &                        checksum = Fhash)
#  else
     &                        MIXING(ng) % AKv)
#  endif
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              IF (Master) THEN
                WRITE (stdout,60) string, TRIM(Vname(1,idVvis)),        &
     &                            InpRec, TRIM(ncname)
              END IF
              exit_flag=2
              ioerror=status
              RETURN
            ELSE
              IF (Master) THEN
#  ifdef CHECKSUM
                WRITE (stdout,70) TRIM(Vname(2,idVvis)), Fmin, Fmax,    &
     &                            Fhash
#  else
                WRITE (stdout,70) TRIM(Vname(2,idVvis)), Fmin, Fmax
#  endif

              END IF
            END IF
#  ifdef DISTRIBUTE
            CALL mp_exchange3d (ng, MyRank, IDmod, 1,                   &
     &                          LBi, UBi, LBj, UBj, 0, N(ng),           &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          MIXING(ng) % AKv)
#  endif
          ELSE
            IF (Master) THEN
              WRITE (stdout,80) string, TRIM(Vname(1,idVvis)),          &
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
!  Read in temperature vertical diffusion.
!
        IF (have_var(idTdif)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idTdif)),   &
     &                        varid)
          IF (foundit) THEN
            gtype=var_flag(varid)*w3dvar
            status=nf_fread3d(ng, IDmod, ncname, ncINPid,               &
     &                        Vname(1,idTdif), varid,                   &
     &                        InpRec, gtype, Vsize,                     &
     &                        LBi, UBi, LBj, UBj, 0, N(ng),             &
     &                        Fscl, Fmin,Fmax,                          &
#  ifdef MASKING
     &                        GRID(ng) % rmask,                         &
#  endif
#  ifdef CHECKSUM
     &                        MIXING(ng) % AKt(:,:,:,itemp),            &
     &                        checksum = Fhash)
#  else
     &                        MIXING(ng) % AKt(:,:,:,itemp))
#  endif
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              IF (Master) THEN
                WRITE (stdout,60) string, TRIM(Vname(1,idTdif)),        &
     &                            InpRec, TRIM(ncname)
              END IF
              exit_flag=2
              ioerror=status
              RETURN
            ELSE
              IF (Master) THEN
#  ifdef CHECKSUM
                WRITE (stdout,70) TRIM(Vname(2,idTdif)), Fmin, Fmax,    &
     &                            Fhash
#  else
                WRITE (stdout,70) TRIM(Vname(2,idTdif)), Fmin, Fmax
#  endif

              END IF
            END IF
#  ifdef DISTRIBUTE
            CALL mp_exchange3d (ng, MyRank, IDmod, 1,                   &
     &                          LBi, UBi, LBj, UBj, 0, N(ng),           &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          MIXING(ng) % AKt(:,:,:,itemp))
#  endif
          ELSE
            IF (Master) THEN
              WRITE (stdout,80) string, TRIM(Vname(1,idTdif)),          &
     &                          TRIM(ncname)
            END IF
            exit_flag=4
            IF (FoundError(exit_flag, nf90_noerr,                       &
     &                     __LINE__, MyFile)) THEN
              RETURN
            END IF
          END IF
        END IF
#  ifdef SALINITY
!
!  Read in salinity vertical diffusion.
!
        IF (have_var(idSdif)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idSdif)),   &
     &                        varid)
          IF (foundit) THEN
            gtype=var_flag(varid)*w3dvar
            status=nf_fread3d(ng, IDmod, ncname, ncINPid,               &
     &                        Vname(1,idSdif), varid,                   &
     &                        InpRec, gtype, Vsize,                     &
     &                        LBi, UBi, LBj, UBj, 0, N(ng),             &
     &                        Fscl, Fmin,Fmax,                          &
#   ifdef MASKING
     &                        GRID(ng) % rmask,                         &
#   endif
#   ifdef CHECKSUM
     &                        MIXING(ng) % AKt(:,:,:,isalt),            &
     &                        checksum = Fhash)
#   else
     &                        MIXING(ng) % AKt(:,:,:,isalt))
#   endif
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              IF (Master) THEN
                WRITE (stdout,60) string, TRIM(Vname(1,idSdif)),        &
     &                            InpRec, TRIM(ncname)
              END IF
              exit_flag=2
              ioerror=status
              RETURN
            ELSE
              IF (Master) THEN
#   ifdef CHECKSUM
                WRITE (stdout,70) TRIM(Vname(2,idSdif)), Fmin, Fmax,    &
     &                            Fhash
#   else
                WRITE (stdout,70) TRIM(Vname(2,idSdif)), Fmin, Fmax
#   endif

              END IF
            END IF
#   ifdef DISTRIBUTE
            CALL mp_exchange3d (ng, MyRank, IDmod, 1,                   &
     &                          LBi, UBi, LBj, UBj, 0, N(ng),           &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          MIXING(ng) % AKt(:,:,:,isalt))
#   endif
          ELSE
            IF (Master) THEN
              WRITE (stdout,80) string, TRIM(Vname(1,idSdif)),          &
     &                          TRIM(ncname)
            END IF
            exit_flag=4
            IF (FoundError(exit_flag, nf90_noerr,                       &
     &                     __LINE__, MyFile)) THEN
              RETURN
            END IF
          END IF
        END IF
#  endif
# endif
# if defined LMD_SKPP
!
!  Read in Hsbl
!
        IF (have_var(idHsbl).and.Perfect3D) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idHsbl)),   &
     &                        varid)
          IF (foundit) THEN
            gtype=var_flag(varid)*r2dvar
            status=nf_fread2d(ng, IDmod, ncname, ncINPid,               &
     &                        Vname(1,idHsbl), varid,                   &
     &                        InpRec, gtype, Vsize,                     &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
#  ifdef MASKING
     &                        GRID(ng) % rmask,                         &
#  endif
#  ifdef CHECKSUM
     &                        MIXING(ng) % Hsbl,                        &
     &                        checksum = Fhash)
#  else
     &                        MIXING(ng) % Hsbl)
#  endif
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              IF (Master) THEN
                WRITE (stdout,60) string, TRIM(Vname(1,idHsbl)),        &
     &                            InpRec, TRIM(ncname)
              END IF
              exit_flag=2
              ioerror=status
              RETURN
            ELSE
              IF (Master) THEN
#  ifdef CHECKSUM
                WRITE (stdout,70) TRIM(Vname(2,idHsbl)), Fmin, Fmax,    &
     &                            Fhash
#  else
                WRITE (stdout,70) TRIM(Vname(2,idHsbl)), Fmin, Fmax
#  endif

              END IF
            END IF
          ELSE
            IF (Master) THEN
              WRITE (stdout,80) string, TRIM(Vname(1,idHsbl)),          &
     &                          TRIM(ncname)
            END IF
            exit_flag=4
            IF (FoundError(exit_flag, nf90_noerr,                       &
     &                     __LINE__, MyFile)) THEN
              RETURN
            END IF
          END IF
        END IF
# endif
# if defined LMD_BKPP
!
!  Read in Hbbl
!
        IF (have_var(idHbbl).and.Perfect3D) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idHbbl)),   &
     &                        varid)
          IF (foundit) THEN
            gtype=var_flag(varid)*r2dvar
            status=nf_fread2d(ng, IDmod, ncname, ncINPid,               &
     &                        Vname(1,idHbbl), varid,                   &
     &                        InpRec, gtype, Vsize,                     &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
#  ifdef MASKING
     &                        GRID(ng) % rmask,                         &
#  endif
#  ifdef CHECKSUM
     &                        MIXING(ng) % Hbbl,                        &
     &                        checksum = Fhash)
#  else
     &                        MIXING(ng) % Hbbl)
#  endif
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              IF (Master) THEN
                WRITE (stdout,60) string, TRIM(Vname(1,idHbbl)),        &
     &                            InpRec, TRIM(ncname)
              END IF
              exit_flag=2
              ioerror=status
              RETURN
            ELSE
              IF (Master) THEN
# ifdef CHECKSUM
                WRITE (stdout,70) TRIM(Vname(2,idHbbl)), Fmin, Fmax,    &
     &                            Fhash
# else
                WRITE (stdout,70) TRIM(Vname(2,idHbbl)), Fmin, Fmax
# endif

              END IF
            END IF
          ELSE
            IF (Master) THEN
              WRITE (stdout,80) string, TRIM(Vname(1,idHbbl)),          &
     &                          TRIM(ncname)
            END IF
            exit_flag=4
            IF (FoundError(exit_flag, nf90_noerr,                       &
     &                     __LINE__, MyFile)) THEN
              RETURN
            END IF
          END IF
        END IF
# endif
# if defined LMD_NONLOCAL && defined PERFECT_RESTART
!
!  Read in Ghats
!
      DO itrc=1,NAT
        IF (have_var(idGhat(itrc))) THEN
          foundit=find_string(var_name, n_var,                          &
     &                        TRIM(Vname(1,idGhat(itrc))), varid)
          IF (foundit) THEN
            gtype=var_flag(varid)*w3dvar
            status=nf_fread3d(ng, IDmod, ncname, ncINPid,               &
     &                        Vname(1,idGhat(itrc)), varid,             &
     &                        InpRec, gtype, Vsize,                     &
     &                        LBi, UBi, LBj, UBj, 0, N(ng),             &
     &                        Fscl, Fmin,Fmax,                          &
#  ifdef MASKING
     &                        GRID(ng) % rmask,                         &
#  endif
#  ifdef CHECKSUM
     &                        MIXING(ng) % Ghats(:,:,:,itrc),           &
     &                        checksum = Fhash)
#  else
     &                        MIXING(ng) % Ghats(:,:,:,itrc))
#  endif
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              IF (Master) THEN
                WRITE (stdout,60) string, TRIM(Vname(1,idGhat(itrc))),  &
     &                            InpRec, TRIM(ncname)
              END IF
              exit_flag=2
              ioerror=status
              RETURN
            ELSE
              IF (Master) THEN
#  ifdef CHECKSUM
                WRITE (stdout,70) TRIM(Vname(2,idGhat(itrc))),          &
     &                            Fmin, Fmax, Fhash
#  else
                WRITE (stdout,70) TRIM(Vname(2,idGhat(itrc))),          &
     &                            Fmin, Fmax
#  endif
              END IF
            END IF
          ELSE
            IF (Master) THEN
              WRITE (stdout,80) string, TRIM(Vname(1,idGhat(itrc))),    &
     &                          TRIM(ncname)
            END IF
            exit_flag=4
            IF (FoundError(exit_flag, nf90_noerr,                       &
     &                     __LINE__, MyFile)) THEN
              RETURN
            END IF
          END IF
        END IF
      END DO
# endif
# if defined GLS_MIXING || defined MY25_MIXING
!
!  Read in turbulent kinetic energy.
!
        IF (get_var(idMtke).and.Perfect3D) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idMtke)),   &
     &                        varid)
          IF (foundit) THEN
            gtype=var_flag(varid)*w3dvar
            status=nf_fread4d(ng, IDmod, ncname, ncINPid,               &
     &                        Vname(1,idMtke), varid,                   &
     &                        InpRec, gtype, Vsize,                     &
     &                        LBi, UBi, LBj, UBj, 0, N(ng), 1, 2,       &
     &                        Fscl, Fmin, Fmax,                         &
#  ifdef MASKING
     &                        GRID(ng) % rmask,                         &
#  endif
#  ifdef CHECKSUM
     &                        MIXING(ng) % tke,                         &
     &                        checksum = Fhash)
#  else
     &                        MIXING(ng) % tke)
#  endif
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              IF (Master) THEN
                WRITE (stdout,60) string, TRIM(Vname(1,idMtke)),        &
     &                            InpRec, TRIM(ncname)
              END IF
              exit_flag=2
              ioerror=status
              RETURN
            ELSE
              IF (Master) THEN
#  ifdef CHECKSUM
                WRITE (stdout,70) TRIM(Vname(2,idMtke)), Fmin, Fmax,    &
     &                            Fhash
#  else
                WRITE (stdout,70) TRIM(Vname(2,idMtke)), Fmin, Fmax
#  endif

              END IF
            END IF
          ELSE
            IF (Master) THEN
              WRITE (stdout,80) string, TRIM(Vname(1,idMtke)),          &
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
!  Read in turbulent kinetic energy time length scale.
!
        IF (get_var(idMtls).and.Perfect3D) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idMtls)),   &
     &                        varid)
          IF (foundit) THEN
            gtype=var_flag(varid)*w3dvar
            status=nf_fread4d(ng, IDmod, ncname, ncINPid,               &
     &                        Vname(1,idMtls), varid,                   &
     &                        InpRec, gtype, Vsize,                     &
     &                        LBi, UBi, LBj, UBj, 0, N(ng), 1, 2,       &
     &                        Fscl, Fmin, Fmax,                         &
#  ifdef MASKING
     &                        GRID(ng) % rmask,                         &
#  endif
#  ifdef CHECKSUM
     &                        MIXING(ng) % gls,                         &
     &                        checksum = Fhash)
#  else
     &                        MIXING(ng) % gls)
#  endif
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              IF (Master) THEN
                WRITE (stdout,60) string, TRIM(Vname(1,idMtls)),        &
     &                            InpRec, TRIM(ncname)
              END IF
              exit_flag=2
              ioerror=status
             RETURN
            ELSE
              IF (Master) THEN
#  ifdef CHECKSUM
                WRITE (stdout,70) TRIM(Vname(2,idMtls)), Fmin, Fmax,    &
     &                            Fhash
#  else
                WRITE (stdout,70) TRIM(Vname(2,idMtls)), Fmin, Fmax
#  endif

              END IF
            END IF
          ELSE
            IF (Master) THEN
              WRITE (stdout,80) string, TRIM(Vname(1,idMtls)),          &
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
!  Read in vertical mixing turbulent length scale.
!
        IF (get_var(idVmLS).and.Perfect3D) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idVmLS)),   &
     &                        varid)
          IF (foundit) THEN
            gtype=var_flag(varid)*w3dvar
            status=nf_fread3d(ng, IDmod, ncname, ncINPid,               &
     &                        Vname(1,idVmLS), varid,                   &
     &                        InpRec, gtype, Vsize,                     &
     &                        LBi, UBi, LBj, UBj, 0, N(ng),             &
     &                        Fscl, Fmin, Fmax,                         &
#  ifdef MASKING
     &                        GRID(ng) % rmask,                         &
#  endif
#  ifdef CHECKSUM
     &                        MIXING(ng) % Lscale,                      &
     &                        checksum = Fhash)
#  else
     &                        MIXING(ng) % Lscale)
#  endif
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              IF (Master) THEN
                WRITE (stdout,60) string, TRIM(Vname(1,idVmLS)),        &
     &                            InpRec, TRIM(ncname)
              END IF
              exit_flag=2
              ioerror=status
              RETURN
            ELSE
              IF (Master) THEN
#  ifdef CHECKSUM
                WRITE (stdout,70) TRIM(Vname(2,idVmLS)), Fmin, Fmax,    &
     &                            Fhash
#  else
                WRITE (stdout,70) TRIM(Vname(2,idVmLS)), Fmin, Fmax
#  endif

              END IF
            END IF
          ELSE
            IF (Master) THEN
              WRITE (stdout,80) string, TRIM(Vname(1,idVmLS)),          &
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
!  Read in turbulent kinetic energy vertical diffusion coefficient.
!
        IF (get_var(idVmKK).and.Perfect3D) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idVmKK)),   &
     &                        varid)
          IF (foundit) THEN
            gtype=var_flag(varid)*w3dvar
            status=nf_fread3d(ng, IDmod, ncname, ncINPid,               &
     &                        Vname(1,idVmKK), varid,                   &
     &                        InpRec, gtype, Vsize,                     &
     &                        LBi, UBi, LBj, UBj, 0, N(ng),             &
     &                        Fscl, Fmin, Fmax,                         &
#  ifdef MASKING
     &                        GRID(ng) % rmask,                         &
#  endif
#  ifdef CHECKSUM
     &                        MIXING(ng) % Akk,                         &
     &                        checksum = Fhash)
#  else
     &                        MIXING(ng) % Akk)
#  endif
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              IF (Master) THEN
                WRITE (stdout,60) string, TRIM(Vname(1,idVmKK)),        &
     &                            InpRec, TRIM(ncname)
              END IF
              exit_flag=2
              ioerror=status
              RETURN
            ELSE
              IF (Master) THEN
#  ifdef CHECKSUM
                WRITE (stdout,70) TRIM(Vname(2,idVmKK)), Fmin, Fmax,    &
     &                            Fhash
#  else
                WRITE (stdout,70) TRIM(Vname(2,idVmKK)), Fmin, Fmax
#  endif

              END IF
            END IF
          ELSE
            IF (Master) THEN
              WRITE (stdout,80) string, TRIM(Vname(1,idVmKK)),          &
     &                          TRIM(ncname)
            END IF
            exit_flag=4
            IF (FoundError(exit_flag, nf90_noerr,                       &
     &                     __LINE__, MyFile)) THEN
              RETURN
            END IF
          END IF
        END IF

#  ifdef GLS_MIXING
!
!  Read in turbulent length scale vertical diffusion coefficient.
!
        IF (get_var(idVmKP).and.Perfect3D) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idVmKP)),   &
     &                        varid)
          IF (foundit) THEN
            gtype=var_flag(varid)*w3dvar
            status=nf_fread3d(ng, IDmod, ncname, ncINPid,               &
     &                        Vname(1,idVmKP), varid,                   &
     &                        InpRec, gtype, Vsize,                     &
     &                        LBi, UBi, LBj, UBj, 0, N(ng),             &
     &                        Fscl, Fmin, Fmax,                         &
#   ifdef MASKING
     &                        GRID(ng) % rmask,                         &
#   endif
#   ifdef CHECKSUM
     &                        MIXING(ng) % Akp,                         &
     &                        checksum = Fhash)
#   else
     &                        MIXING(ng) % Akp)
#   endif
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              IF (Master) THEN
                WRITE (stdout,60) string, TRIM(Vname(1,idVmKP)),        &
     &                            InpRec, TRIM(ncname)
              END IF
              exit_flag=2
              ioerror=status
              RETURN
            ELSE
              IF (Master) THEN
#   ifdef CHECKSUM
                WRITE (stdout,70) TRIM(Vname(2,idVmKP)), Fmin, Fmax,    &
     &                            Fhash
#   else
                WRITE (stdout,70) TRIM(Vname(2,idVmKP)), Fmin, Fmax
#   endif

              END IF
            END IF
          ELSE
            IF (Master) THEN
              WRITE (stdout,80) string, TRIM(Vname(1,idVmKP)),          &
     &                          TRIM(ncname)
            END IF
            exit_flag=4
            IF (FoundError(exit_flag, nf90_noerr,                       &
     &                     __LINE__, MyFile)) THEN
              RETURN
            END IF
          END IF
        END IF
#  endif
# endif
# ifdef SEDIMENT
!
!  Read in nonlinear sediment fraction of each size class in each bed
!  layer.
!
        DO i=1,NST
          IF (get_var(idfrac(i))) THEN
            foundit=find_string(var_name, n_var,                        &
     &                          TRIM(Vname(1,idfrac(i))), varid)
            IF (foundit) THEN
              gtype=var_flag(varid)*b3dvar
              status=nf_fread3d(ng, IDmod, ncname, ncINPid,             &
     &                          Vname(1,idfrac(i)), varid,              &
     &                          InpRec, gtype, Vsize,                   &
     &                          LBi, UBi, LBj, UBj, 1, Nbed,            &
     &                          Fscl, Fmin, Fmax,                       &
#  ifdef MASKING
     &                          GRID(ng) % rmask,                       &
#  endif
#  ifdef CHECKSUM
     &                          SEDBED(ng) % bed_frac(:,:,:,i),         &
     &                          checksum = Fhash)
#  else
     &                          SEDBED(ng) % bed_frac(:,:,:,i))
#  endif
              IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
                IF (Master) THEN
                  WRITE (stdout,60) string, TRIM(Vname(1,idfrac(i))),   &
     &                              InpRec, TRIM(ncname)
                END IF
                exit_flag=2
                ioerror=status
                RETURN
              ELSE
                IF (Master) THEN
#  ifdef CHECKSUM
                  WRITE (stdout,70) TRIM(Vname(2,idfrac(i))),           &
     &                              Fmin, Fmax, Fhash
#  else
                  WRITE (stdout,70) TRIM(Vname(2,idfrac(i))),           &
     &                              Fmin, Fmax
#  endif

                END IF
              END IF
            ELSE
              IF (Master) THEN
                WRITE (stdout,80) string, TRIM(Vname(1,idfrac(i))),     &
     &                            TRIM(ncname)
              END IF
              exit_flag=4
              IF (FoundError(exit_flag, nf90_noerr,                     &
     &                       __LINE__, MyFile)) THEN
                RETURN
              END IF
            END IF
          END IF
!
!  Read in nonlinear sediment mass of each size class in each bed layer.
!
          IF (get_var(idBmas(i))) THEN
            foundit=find_string(var_name, n_var,                        &
     &                          TRIM(Vname(1,idBmas(i))), varid)
            IF (foundit) THEN
              gtype=var_flag(varid)*b3dvar
              status=nf_fread3d(ng, IDmod, ncname, ncINPid,             &
     &                          Vname(1,idBmas(i)), varid,              &
     &                          InpRec, gtype, Vsize,                   &
     &                          LBi, UBi, LBj, UBj, 1, Nbed,            &
     &                          Fscl, Fmin, Fmax,                       &
#  ifdef MASKING
     &                          GRID(ng) % rmask,                       &
#  endif
#  ifdef CHECKSUM
     &                          SEDBED(ng) % bed_mass(:,:,:,Tindex,i),  &
     &                          checksum = Fhash)
#  else
     &                          SEDBED(ng) % bed_mass(:,:,:,Tindex,i))
#  endif
              IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
                IF (Master) THEN
                  WRITE (stdout,60) string, TRIM(Vname(1,idBmas(i))),   &
     &                              InpRec, TRIM(ncname)
                END IF
                exit_flag=2
                ioerror=status
                RETURN
              ELSE
                IF (Master) THEN
# ifdef CHECKSUM
                  WRITE (stdout,70) TRIM(Vname(2,idBmas(i))),           &
     &                              Fmin, Fmax, Fhash
# else
                  WRITE (stdout,70) TRIM(Vname(2,idBmas(i))),           &
     &                              Fmin, Fmax
# endif

                END IF
              END IF
            ELSE
              IF (Master) THEN
                WRITE (stdout,80) string, TRIM(Vname(1,idBmas(i))),     &
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
!
!  Read in nonlinear sediment properties in each bed layer.
!
        DO i=1,MBEDP
          IF (get_var(idSbed(i))) THEN
            foundit=find_string(var_name, n_var,                        &
     &                          TRIM(Vname(1,idSbed(i))), varid)
            IF (foundit) THEN
              gtype=var_flag(varid)*b3dvar
              status=nf_fread3d(ng, IDmod, ncname, ncINPid,             &
     &                          Vname(1,idSbed(i)), varid,              &
     &                          InpRec, gtype, Vsize,                   &
     &                          LBi, UBi, LBj, UBj, 1, Nbed,            &
     &                          Fscl, Fmin, Fmax,                       &
#  ifdef MASKING
     &                          GRID(ng) % rmask,                       &
#  endif
#  ifdef CHECKSUM
     &                          SEDBED(ng) % bed(:,:,:,i),              &
     &                          checksum = Fhash)
#  else
     &                          SEDBED(ng) % bed(:,:,:,i))
#  endif
              IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
                IF (Master) THEN
                  WRITE (stdout,60) string, TRIM(Vname(1,idSbed(i))),   &
     &                              InpRec, TRIM(ncname)
                END IF
                exit_flag=2
                ioerror=status
                RETURN
              ELSE
                IF (Master) THEN
#  ifdef CHECKSUM
                  WRITE (stdout,70) TRIM(Vname(2,idSbed(i))),           &
     &                              Fmin, Fmax, Fhash
#  else
                  WRITE (stdout,70) TRIM(Vname(2,idSbed(i))),           &
     &                              Fmin, Fmax
#  endif

                END IF
              END IF
            ELSE
              IF (Master) THEN
                WRITE (stdout,80) string, TRIM(Vname(1,idSbed(i))),     &
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

#  ifdef BEDLOAD
!
!  Read in nonlinear sediment fraction of bed load.
!
        DO i=1,NST
          IF (get_var(idUbld(i))) THEN
            foundit=find_string(var_name, n_var,                        &
     &                          TRIM(Vname(1,idUbld(i))), varid)
            IF (foundit) THEN
              gtype=var_flag(varid)*u2dvar
              status=nf_fread2d(ng, IDmod, ncname, ncINPid,             &
     &                          Vname(1,idUbld(i)), varid,              &
     &                          InpRec, gtype, Vsize,                   &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Fscl, Fmin, Fmax,                       &
#   ifdef MASKING
     &                          GRID(ng) % umask,                       &
#   endif
#   ifdef CHECKSUM
     &                          SEDBED(ng) % bedldu(:,:,i),             &
     &                          checksum = Fhash)
#   else
     &                          SEDBED(ng) % bedldu(:,:,i))
#   endif
              IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
                IF (Master) THEN
                  WRITE (stdout,60) string, TRIM(Vname(1,idUbld(i))),   &
     &                              InpRec, TRIM(ncname)
                END IF
                exit_flag=2
                ioerror=status
                RETURN
              ELSE
                IF (Master) THEN
#   ifdef CHECKSUM
                  WRITE (stdout,70) TRIM(Vname(2,idUbld(i))),           &
     &                              Fmin, Fmax, Fhash
#   else
                  WRITE (stdout,70) TRIM(Vname(2,idUbld(i))),           &
     &                              Fmin, Fmax
#   endif

                END IF
              END IF
            ELSE
              IF (Master) THEN
                WRITE (stdout,80) string, TRIM(Vname(1,idUbld(i))),     &
     &                            TRIM(ncname)
              END IF
              exit_flag=4
              IF (FoundError(exit_flag, nf90_noerr,                     &
     &                       __LINE__, MyFile)) THEN
                RETURN
              END IF
            END IF
          END IF
!
          IF (get_var(idVbld(i))) THEN
            foundit=find_string(var_name, n_var,                        &
     &                          TRIM(Vname(1,idVbld(i))), varid)
            IF (foundit) THEN
              gtype=var_flag(varid)*v2dvar
              status=nf_fread2d(ng, IDmod, ncname, ncINPid,             &
     &                          Vname(1,idVbld(i)), varid,              &
     &                          InpRec, gtype, Vsize,                   &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Fscl, Fmin, Fmax,                       &
#   ifdef MASKING
     &                          GRID(ng) % vmask,                       &
#   endif
#   ifdef CHECKSUM
     &                          SEDBED(ng) % bedldv(:,:,i),             &
     &                          checksum = Fhash)
#   else
     &                          SEDBED(ng) % bedldv(:,:,i))
#   endif
              IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
                IF (Master) THEN
                  WRITE (stdout,60) string, TRIM(Vname(1,idVbld(i))),   &
     &                              InpRec, TRIM(ncname)
                END IF
                exit_flag=2
                ioerror=status
                RETURN
              ELSE
                IF (Master) THEN
#   ifdef CHECKSUM
                  WRITE (stdout,70) TRIM(Vname(2,idVbld(i))),           &
     &                              Fmin, Fmax, Fhash
#   else
                  WRITE (stdout,70) TRIM(Vname(2,idVbld(i))),           &
     &                              Fmin, Fmax
#   endif

                END IF
              END IF
            ELSE
              IF (Master) THEN
                WRITE (stdout,80) string, TRIM(Vname(1,idVbld(i))),     &
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
#  endif
# endif

# if defined SEDIMENT || defined BBL_MODEL
!
!  Read in nonlinear sediment properties in exposed bed layer.
!
        DO i=1,MBOTP
          IF (get_var(idBott(i)).and.have_var(idBott(i))) THEN
            foundit=find_string(var_name, n_var,                        &
     &                          TRIM(Vname(1,idBott(i))), varid)
            IF (foundit) THEN
              gtype=var_flag(varid)*r2dvar
              status=nf_fread2d(ng, IDmod, ncname, ncINPid,             &
     &                          Vname(1,idBott(i)), varid,              &
     &                          InpRec, gtype, Vsize,                   &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Fscl, Fmin, Fmax,                       &
#  ifdef MASKING
     &                          GRID(ng) % rmask,                       &
#  endif
#  ifdef CHECKSUM
     &                          SEDBED(ng) % bottom(:,:,i),             &
     &                          checksum = Fhash)
#  else
     &                          SEDBED(ng) % bottom(:,:,i))
#  endif
              IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
                IF (Master) THEN
                  WRITE (stdout,60) string, TRIM(Vname(1,idBott(i))),   &
     &                              InpRec, TRIM(ncname)
                END IF
                exit_flag=2
                ioerror=status
                RETURN
              ELSE
                IF (Master) THEN
#  ifdef CHECKSUM
                  WRITE (stdout,70) TRIM(Vname(2,idBott(i))),           &
     &                              Fmin, Fmax, Fhash
#  else
                  WRITE (stdout,70) TRIM(Vname(2,idBott(i))),           &
     &                              Fmin, Fmax
#  endif

                END IF
              END IF
            ELSE
              IF (Master) THEN
                WRITE (stdout,80) string, TRIM(Vname(1,idBott(i))),     &
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

#  if defined ESTUARYBGC && defined SAV_BIOMASS
!
!  Read in AGB and BGB values
!
        foundit=find_string(var_name, n_var, TRIM(Vname(1,idsagb)),     & 
     &                        varid)
        gtype=var_flag(varid)*r2dvar
        status=nf_fread2d(ng, IDmod, ncname, ncINPid,                   &
     &                Vname(1,idsagb), varid, InpRec, gtype,            &
     &                Vsize, LBi, UBi, LBj, UBj,                        &
     &                Fscl, Fmin, Fmax,                                 &
#   ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#   endif
     &                 OCEAN(ng)%AGB)
        IF (FoundError(exit_flag, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN 
            WRITE (stdout,60) string, TRIM(Vname(1,idsagb )),           &
     &                        InpRec, TRIM(ncname)
          END IF
          exit_flag=2
          ioerror=status
          RETURN
        ELSE
          IF (Master) THEN 
            WRITE (stdout,70) TRIM(Vname(2,idsagb)), Fmin, Fmax 
          END IF
        END IF
!
! BGB 
!  
        foundit=find_string(var_name, n_var, TRIM(Vname(1,idsbgb)),     & 
     &                        varid)
        gtype=var_flag(varid)*r2dvar
        status=nf_fread2d(ng, IDmod, ncname, ncINPid,                   &
     &                Vname(1,idsbgb), varid, InpRec, gtype,            &
     &                Vsize, LBi, UBi, LBj, UBj,                        &
     &                Fscl, Fmin, Fmax,                                 &
#   ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#   endif
     &               OCEAN(ng)%BGB)
        IF (FoundError(exit_flag, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN 
            WRITE (stdout,60) string, TRIM(Vname(1,idsbgb )),           &
     &                        InpRec, TRIM(ncname)
          END IF
          exit_flag=2
          ioerror=status
          RETURN
        ELSE
          IF (Master) THEN 
            WRITE (stdout,70) TRIM(Vname(2,idsbgb)), Fmin, Fmax 
          END IF
        END IF
!
! Read in EPB (Epiphyte Biomass) 
!  
        foundit=find_string(var_name, n_var, TRIM(Vname(1,idsepb)),     & 
     &                        varid)
        gtype=var_flag(varid)*r2dvar
        status=nf_fread2d(ng, IDmod, ncname, ncINPid,                   &
     &                Vname(1,idsepb), varid, InpRec, gtype,            &
     &                Vsize, LBi, UBi, LBj, UBj,                        &
     &                Fscl, Fmin, Fmax,                                 &
#   ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#   endif
     &               OCEAN(ng)%EPB)
        IF (FoundError(exit_flag, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN 
            WRITE (stdout,60) string, TRIM(Vname(1,idsepb )),           &
     &                        InpRec, TRIM(ncname)
          END IF
          exit_flag=2
          ioerror=status
          RETURN
        ELSE
          IF (Master) THEN 
            WRITE (stdout,70) TRIM(Vname(2,idsepb)), Fmin, Fmax 
          END IF
        END IF
!
! DINsed
!  
        foundit=find_string(var_name, n_var, TRIM(Vname(1,iddins)),     & 
     &                        varid)
        gtype=var_flag(varid)*r3dvar
        status=nf_fread3d(ng, IDmod, ncname, ncINPid,                   &
     &                Vname(1,iddins), varid, InpRec, gtype,            &
     &                Vsize, LBi, UBi, LBj, UBj,  1, N(ng),             &
     &                Fscl, Fmin, Fmax,                                 &
#   ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#   endif
     &               OCEAN(ng)%DINsed)
        IF (FoundError(exit_flag, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN 
            WRITE (stdout,60) string, TRIM(Vname(1,iddins )),           &
     &                        InpRec, TRIM(ncname)
          END IF
          exit_flag=2
          ioerror=status
          RETURN
        ELSE
          IF (Master) THEN 
            WRITE (stdout,70) TRIM(Vname(2,iddins)), Fmin, Fmax 
          END IF
        END IF
!
#  endif 


# ifdef VEGETATION 
#  if defined VEG_DRAG || defined VEG_BIOMASS
!
!  Read in submerged aquatic plant properties for each vegetation type.
!
        DO i=1,NVEGP
          IF (get_var(idvprp(i))) THEN
            foundit=find_string(var_name, n_var,                        &
     &                          TRIM(Vname(1,idvprp(i))), varid)
            IF (foundit) THEN
              gtype=var_flag(varid)*a3dvar
              status=nf_fread3d(ng, IDmod, ncname, ncINPid,             &
     &                          Vname(1,idvprp(i)), varid,              &
     &                          InpRec, gtype, Vsize,                   &
     &                          LBi, UBi, LBj, UBj, 1, NVEG,            &
     &                          Fscl, Fmin, Fmax,                       &
#   ifdef MASKING
     &                          GRID(ng) % rmask,                       &
#   endif
#   ifdef CHECKSUM
     &                          VEG(ng) % plant(:,:,:,i),               &
     &                          checksum = Fhash)
#   else
     &                          VEG(ng) % plant(:,:,:,i))
#   endif 
              IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
                IF (Master) THEN
                  WRITE (stdout,60) string, TRIM(Vname(1,idvprp(i))),   &
     &                              InpRec, TRIM(ncname)
                END IF
                exit_flag=2
                ioerror=status
                RETURN
              ELSE
                IF (Master) THEN
#   ifdef CHECKSUM
                  WRITE (stdout,70) TRIM(Vname(2,idvprp(i))),           &
     &                              Fmin, Fmax, Fhash
#   else
                  WRITE (stdout,70) TRIM(Vname(2,idvprp(i))),           &
     &                              Fmin, Fmax
#   endif
                END IF
              END IF
            ELSE
              IF (Master) THEN
                WRITE (stdout,80) string, TRIM(Vname(1,idvprp(i))),     &
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
#  endif 
#  ifdef MARSH_DYNAMICS
!
!  Read in masking for determining marsh interface. 
!
        IF (get_var(idTims)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idTims)),   & 
     &                        varid)
          IF (foundit) THEN
            gtype=var_flag(varid)*r2dvar
            status=nf_fread2d(ng, IDmod, ncname, ncINPid,               &
     &                        Vname(1,idTims), varid,                   &
     &                        InpRec, gtype,  Vsize,                    &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
#   ifdef MASKING
     &                        GRID(ng) % rmask,                         &
#   endif
#   ifdef CHECKSUM
     &                        VEG(ng) % marsh_mask,                     &
     &                        checksum = Fhash)
#   else
     &                        VEG(ng) % marsh_mask)
#   endif
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              IF (Master) THEN 
                WRITE (stdout,60) string, TRIM(Vname(1,idTims )),       &
     &                            InpRec, TRIM(ncname)
              END IF
              exit_flag=2
              ioerror=status
              RETURN
            ELSE
              IF (Master) THEN 
#   ifdef CHECKSUM
                WRITE (stdout,70) TRIM(Vname(2,idTims)),                &
     &                            Fmin, Fmax, Fhash
#   else
                WRITE (stdout,70) TRIM(Vname(2,idTims)), Fmin, Fmax 
#   endif
              END IF
            END IF
          ELSE
            IF (Master) THEN
              WRITE (stdout,80) string, TRIM(Vname(1,idTims)),          &
     &                          TRIM(ncname)
            END IF
            exit_flag=4
            IF (FoundError(exit_flag, nf90_noerr,                       &
     &                     __LINE__, MyFile)) THEN
              RETURN
            END IF
          END IF
        END IF
#  endif 
# endif 

# ifdef ICE_MODEL
!
!  Read sea ice model state variables.
!
        DO i=1,nIceS
          IF (iSice(i).gt.0) THEN
            ifield=iSice(i)
            IF (get_var(ifield)) THEN
              foundit=find_string(var_name, n_var,                      &
     &                            TRIM(Vname(1,ifield)), varid)
              IF (foundit) THEN
                SELECT CASE (i)
                  CASE (isUice)
                    gtype=var_flag(varid)*u2dvar
                    status=nf_fread2d(ng, IDmod, ncname, ncINPid,       &
     &                                Vname(1,ifield), varid,           &
     &                                InpRec, gtype, Vsize,             &
     &                                LBi, UBi, LBj, UBj,               &
     &                                Fscl, Fmin, Fmax,                 &
#  ifdef MASKING
     &                                GRID(ng) % umask,                 &
#  endif
#  ifdef CHECKSUM
     &                                ICE(ng) % Si(:,:,Tindex,i),       &
     &                                checksum = Fhash)
#  else
     &                                ICE(ng) % Si(:,:,Tindex,i))
#  endif
                  CASE (isVice)
                    gtype=var_flag(varid)*v2dvar
                    status=nf_fread2d(ng, IDmod, ncname, ncINPid,       &
     &                                Vname(1,ifield), varid,           &
     &                                InpRec, gtype, Vsize,             &
     &                                LBi, UBi, LBj, UBj,               &
     &                                Fscl, Fmin, Fmax,                 &
#  ifdef MASKING
     &                                GRID(ng) % vmask,                 &
#  endif
#  ifdef CHECKSUM
     &                                ICE(ng) % Si(:,:,Tindex,i),       &
     &                                checksum = Fhash)
#  else
     &                                ICE(ng) % Si(:,:,Tindex,i))
#  endif
                  CASE DEFAULT
                    gtype=var_flag(varid)*r2dvar
                    status=nf_fread2d(ng, IDmod, ncname, ncINPid,       &
     &                                Vname(1,ifield), varid,           &
     &                                InpRec, gtype, Vsize,             &
     &                                LBi, UBi, LBj, UBj,               &
     &                                Fscl, Fmin, Fmax,                 &
#  ifdef MASKING
     &                                GRID(ng) % rmask,                 &
#  endif
#  ifdef CHECKSUM
     &                                ICE(ng) % Si(:,:,Tindex,i),       &
     &                                checksum = Fhash)
#  else
     &                                ICE(ng) % Si(:,:,Tindex,i))
#  endif
                END SELECT
!
                IF (FoundError(status, nf90_noerr,                      &
     &                         __LINE__, MyFile)) THEN
                  IF (Master) THEN
                    WRITE (stdout,60) string, TRIM(Vname(1,ifield)),    &
     &                                InpRec, TRIM(ncname)
                  END IF
                  exit_flag=2
                  ioerror=status
                  RETURN
                ELSE
                  IF (Master) THEN
#  ifdef CHECKSUM
                    WRITE (stdout,70) TRIM(Vname(2,ifield)),            &
     &                                Fmin, Fmax, Fhash
#  else
                    WRITE (stdout,70) TRIM(Vname(2,ifield)), Fmin, Fmax
#  endif
                  END IF
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
