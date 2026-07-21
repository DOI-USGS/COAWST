/*
** git $Id$
*************************************************** Hernan G. Arango ***
** Copyright (c) 2002-2026 The ROMS Group                             **
**    Licensed under a MIT/X style license                            **
**    See License_ROMS.md                                             **
************************************************************************
**                                                                    **
** NetCDF PIO Library Processing:                                     **
**                                                                    **
** It process the reading of ROMS FRC state from input NetCDF files.  **
**                                                                    **
** It is used to read in 4D-Var impulse forcing terms for the adjoint **
** or tangent linear kernels. Data is loaded into "f_" arrays.        **
**                                                                    **
************************************************************************
*/

!
!  Set number of records available.
!
        NrecFrc(ng)=Nrec
!
!  Read in next impulse forcing time to process.
!
        CALL pio_netcdf_get_time (ng, IDmod, ncname, TRIM(tvarnam),     &
     &                            Rclock%DateNumber, FrcTime(ng:),      &
     &                            pioFile = pioFile,                    &
     &                            start = (/InpRec/),                   &
     &                            total = (/1/))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Read in free-surface impulse forcing.
!
        IF (get_var(idFsur)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idFsur)),   &
     &                        vindex)
          IF (foundit) THEN
            my_pioVar%vd=var_desc(vindex)
            my_pioVar%gtype=r2dvar
            IF (KIND(OCEAN(ng)%f_zeta).eq.8) THEN
              my_pioVar%dkind=PIO_double
              ioDesc => ioDesc_dp_r2dvar(ng)
            ELSE
              my_pioVar%dkind=PIO_real
              ioDesc => ioDesc_sp_r2dvar(ng)
            END IF
!
            status=nf_fread2d(ng, IDmod, ncname, pioFile,               &
     &                        Vname(1,idFsur), my_pioVar,               &
     &                        InpRec, ioDesc, Vsize,                    &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
#ifdef MASKING
     &                        GRID(ng) % rmask,                         &
#endif
#ifdef CHECKSUM
     &                        OCEAN(ng) % f_zeta,                       &
     &                        checksum = Fhash)
#else
     &                        OCEAN(ng) % f_zeta)
#endif
            IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
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
            IF (FoundError(exit_flag, PIO_noerr,                        &
     &                     __LINE__, MyFile)) THEN
              RETURN
            END IF
          END IF
        END IF

#ifndef SOLVE3D
!
!  Read in 2D U-momentum impulse forcing.
!
        IF (get_var(idUbar)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idUbar)),   &
     &                        vindex)
          IF (foundit) THEN
            my_pioVar%vd=var_desc(vindex)
            my_pioVar%gtype=u2dvar
            IF (KIND(OCEAN(ng)%f_ubar).eq.8) THEN
              my_pioVar%dkind=PIO_double
              ioDesc => ioDesc_dp_u2dvar(ng)
            ELSE
              my_pioVar%dkind=PIO_real
              ioDesc => ioDesc_sp_u2dvar(ng)
            END IF
!
            status=nf_fread2d(ng, IDmod, ncname, pioFile,               &
     &                        Vname(1,idUbar), my_pioVar,               &
     &                        InpRec, ioDesc, Vsize,                    &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
# ifdef MASKING
     &                        GRID(ng) % umask,                         &
# endif
# ifdef CHECKSUM
     &                        OCEAN(ng) % f_ubar,                       &
     &                        checksum = Fhash)
# else
     &                        OCEAN(ng) % f_ubar)
# endif
            IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
              IF (Master) THEN
                WRITE (stdout,60) string, TRIM(Vname(1,idUbar)),        &
     &                            InpRec, TRIM(ncname)
              END IF
              exit_flag=2
              ioerror=status
              RETURN
            ELSE
              IF (Master) THEN
# ifdef CHECKSUM
                WRITE (stdout,70) TRIM(Vname(2,idUbar)), Fmin, Fmax,    &
     &                            Fhash
# else
                WRITE (stdout,70) TRIM(Vname(2,idUbar)), Fmin, Fmax
# endif
              END IF
            END IF
          ELSE
            IF (Master) THEN
              WRITE (stdout,80) string, TRIM(Vname(1,idUbar)),          &
     &                          TRIM(ncname)
            END IF
            exit_flag=4
            IF (FoundError(exit_flag, PIO_noerr,                        &
     &                     __LINE__, MyFile)) THEN
              RETURN
            END IF
          END IF
        END IF
!
!  Read in 2D V-momentum impulse forcing.
!
        IF (get_var(idVbar)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idVbar)),   &
     &                        vindex)
          IF (foundit) THEN
            my_pioVar%vd=var_desc(vindex)
            my_pioVar%gtype=v2dvar
            IF (KIND(OCEAN(ng)%f_vbar).eq.8) THEN
              my_pioVar%dkind=PIO_double
              ioDesc => ioDesc_dp_v2dvar(ng)
            ELSE
              my_pioVar%dkind=PIO_real
              ioDesc => ioDesc_sp_v2dvar(ng)
            END IF
!
            status=nf_fread2d(ng, IDmod, ncname, pioFile,               &
     &                        Vname(1,idVbar), my_pioVar,               &
     &                        InpRec, ioDesc, Vsize,                    &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
# ifdef MASKING
     &                        GRID(ng) % vmask,                         &
# endif
# ifdef CHECKSUM
     &                        OCEAN(ng) % f_vbar,                       &
     &                        checksum = Fhash)
# else
     &                        OCEAN(ng) % f_vbar)
# endif
            IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
              IF (Master) THEN
                WRITE (stdout,60) string, TRIM(Vname(1,idVbar)),        &
     &                            InpRec, TRIM(ncname)
              END IF
              exit_flag=2
              ioerror=status
              RETURN
            ELSE
              IF (Master) THEN
# ifdef CHECKSUM
                WRITE (stdout,70) TRIM(Vname(2,idVbar)), Fmin, Fmax,    &
     &                            Fhash
# else
                WRITE (stdout,70) TRIM(Vname(2,idVbar)), Fmin, Fmax
# endif
              END IF
            END IF
          ELSE
            IF (Master) THEN
              WRITE (stdout,80) string, TRIM(Vname(1,idVbar)),          &
     &                          TRIM(ncname)
            END IF
            exit_flag=4
            IF (FoundError(exit_flag, PIO_noerr,                        &
     &                     __LINE__, MyFile)) THEN
              RETURN
            END IF
          END IF
        END IF
#endif
#ifdef SOLVE3D
!
!  Read in 3D U-momentum impulse forcing.
!
        IF (get_var(idUvel)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idUvel)),   &
     &                        vindex)
          IF (foundit) THEN
            my_pioVar%vd=var_desc(vindex)
            my_pioVar%gtype=u3dvar
            IF (KIND(OCEAN(ng)%f_u).eq.8) THEN
              my_pioVar%dkind=PIO_double
              ioDesc => ioDesc_dp_u3dvar(ng)
            ELSE
              my_pioVar%dkind=PIO_real
              ioDesc => ioDesc_sp_u3dvar(ng)
            END IF
!
            status=nf_fread3d(ng, IDmod, ncname, pioFile,               &
     &                        Vname(1,idUvel), my_pioVar,               &
     &                        InpRec, ioDesc, Vsize,                    &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        Fscl, Fmin, Fmax,                         &
# ifdef MASKING
     &                        GRID(ng) % umask,                         &
# endif
# ifdef CHECKSUM
     &                        OCEAN(ng) % f_u,                          &
     &                        checksum = Fhash)
# else
     &                        OCEAN(ng) % f_u)
# endif
            IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
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
            IF (FoundError(exit_flag, PIO_noerr,                        &
     &                     __LINE__, MyFile)) THEN
              RETURN
            END IF
          END IF
        END IF
!
!  Read in 3D V-momentum impulse forcing.
!
        IF (get_var(idVvel)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idVvel)),   &
     &                        vindex)
          IF (foundit) THEN
            my_pioVar%vd=var_desc(vindex)
            my_pioVar%gtype=v3dvar
            IF (KIND(OCEAN(ng)%f_v).eq.8) THEN
              my_pioVar%dkind=PIO_double
              ioDesc => ioDesc_dp_v3dvar(ng)
            ELSE
              my_pioVar%dkind=PIO_real
              ioDesc => ioDesc_sp_v3dvar(ng)
            END IF
!
            status=nf_fread3d(ng, IDmod, ncname, pioFile,               &
     &                        Vname(1,idVvel), my_pioVar,               &
     &                        InpRec, ioDesc, Vsize,                    &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        Fscl, Fmin, Fmax,                         &
# ifdef MASKING
     &                        GRID(ng) % vmask,                         &
# endif
# ifdef CHECKSUM
     &                        OCEAN(ng) % f_v,                          &
     &                        checksum = Fhash)
# else
     &                        OCEAN(ng) % f_v)
# endif
            IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
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
            IF (FoundError(exit_flag, PIO_noerr,                        &
     &                     __LINE__, MyFile)) THEN
              RETURN
            END IF
          END IF
        END IF
!
!  Read in tracer variables impulse forcing.
!
        DO itrc=1,NT(ng)
          IF (get_var(idTvar(itrc))) THEN
            foundit=find_string(var_name, n_var,                        &
     &                          TRIM(Vname(1,idTvar(itrc))), vindex)
            IF (foundit) THEN
              my_pioVar%vd=var_desc(vindex)
              my_pioVar%gtype=r3dvar
              IF (KIND(OCEAN(ng)%f_t).eq.8) THEN
                my_pioVar%dkind=PIO_double
                ioDesc => ioDesc_dp_r3dvar(ng)
              ELSE
                my_pioVar%dkind=PIO_real
                ioDesc => ioDesc_sp_r3dvar(ng)
              END IF
!
              status=nf_fread3d(ng, IDmod, ncname, pioFile,             &
     &                          Vname(1,idTvar(itrc)), my_pioVar,       &
     &                          InpRec, ioDesc, Vsize,                  &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          Fscl, Fmin, Fmax,                       &
# ifdef MASKING
     &                          GRID(ng) % rmask,                       &
# endif
# ifdef CHECKSUM
     &                          OCEAN(ng) % f_t(:,:,:,itrc),            &
     &                          checksum = Fhash)
# else
     &                          OCEAN(ng) % f_t(:,:,:,itrc))
# endif
              IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
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
              IF (FoundError(exit_flag, PIO_noerr,                      &
     &                       __LINE__, MyFile)) THEN
                RETURN
              END IF
            END IF
          END IF
        END DO
#endif
