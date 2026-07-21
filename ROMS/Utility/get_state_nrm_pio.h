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
** It process the reading of ROMS NRM state from input NetCDF files.  **
**                                                                    **
** It is used to read in 4D-Var background-error covariance matrix    **
** normalization coefficients. Data is loaded into "b_" arrays.       **
**                                                                    **
************************************************************************
*/

!
!  Read in free-surface normalization factor.
!
        IF (get_var(idFsur).and.((model.eq.14).or.(model.eq.15))) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idFsur)),   &
     &                        vindex)
          IF (foundit) THEN
            my_pioVar%vd=var_desc(vindex)
            my_pioVar%gtype=r2dvar
            IF (KIND(OCEAN(ng)%b_zeta).eq.8) THEN
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
     &                        OCEAN(ng) % b_zeta(:,:,Tindex),           &
     &                        checksum = Fhash)
#else
     &                        OCEAN(ng) % b_zeta(:,:,Tindex))
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
#ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, MyRank, IDmod, 1,                   &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          OCEAN(ng) % b_zeta(:,:,Tindex))
#endif
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
!
!  Read in 2D U-momentum component normalization factor.
!
        IF (get_var(idUbar).and.((model.eq.14).or.(model.eq.15))) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idUbar)),   &
     &                        vindex)
          IF (foundit) THEN
            my_pioVar%vd=var_desc(vindex)
            my_pioVar%gtype=u2dvar
            IF (KIND(OCEAN(ng)%b_ubar).eq.8) THEN
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
#ifdef MASKING
     &                        GRID(ng) % umask,                         &
#endif
#ifdef CHECKSUM
     &                        OCEAN(ng) % b_ubar(:,:,Tindex),           &
     &                        checksum = Fhash)
#else
     &                        OCEAN(ng) % b_ubar(:,:,Tindex))
#endif
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
#ifdef CHECKSUM
                WRITE (stdout,70) TRIM(Vname(2,idUbar)), Fmin, Fmax,    &
     &                            Fhash
#else
                WRITE (stdout,70) TRIM(Vname(2,idUbar)), Fmin, Fmax
#endif
              END IF
            END IF
#ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, MyRank, IDmod, 1,                   &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          OCEAN(ng) % b_ubar(:,:,Tindex))
#endif
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
!  Read in 2D V-momentum component normalization factor.
!
        IF (get_var(idVbar).and.((model.eq.14).or.(model.eq.15))) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idVbar)),   &
     &                        vindex)
          IF (foundit) THEN
            my_pioVar%vd=var_desc(vindex)
            my_pioVar%gtype=v2dvar
            IF (KIND(OCEAN(ng)%b_vbar).eq.8) THEN
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
#ifdef MASKING
     &                        GRID(ng) % vmask,                         &
#endif
#ifdef CHECKSUM
     &                        OCEAN(ng) % b_vbar(:,:,Tindex),           &
     &                        checksum = Fhash)
#else
     &                        OCEAN(ng) % b_vbar(:,:,Tindex))
#endif
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
#ifdef CHECKSUM
                WRITE (stdout,70) TRIM(Vname(2,idVbar)), Fmin, Fmax,    &
     &                            Fhash
#else
                WRITE (stdout,70) TRIM(Vname(2,idVbar)), Fmin, Fmax
#endif
              END IF
            END IF
#ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, MyRank, IDmod, 1,                   &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          OCEAN(ng) % b_vbar(:,:,Tindex))
#endif
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

#ifdef SOLVE3D
!
!  Read in 3D U-momentum component normalization factor.
!
        IF (get_var(idUvel).and.((model.eq.14).or.(model.eq.15))) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idUvel)),   &
     &                        vindex)
          IF (foundit) THEN
            my_pioVar%vd=var_desc(vindex)
            my_pioVar%gtype=u3dvar
            IF (KIND(OCEAN(ng)%b_u).eq.8) THEN
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
     &                        OCEAN(ng) % b_u(:,:,:,Tindex),            &
     &                        checksum = Fhash)
# else
     &                        OCEAN(ng) % b_u(:,:,:,Tindex))
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
# ifdef DISTRIBUTE
            CALL mp_exchange3d (ng, MyRank, IDmod, 1,                   &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          OCEAN(ng) % b_u(:,:,:,Tindex))
# endif
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
!  Read in 3D V-momentum component normalization factor.
!
        IF (get_var(idVvel).and.((model.eq.14).or.(model.eq.15))) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idVvel)),   &
     &                        vindex)
          IF (foundit) THEN
            my_pioVar%vd=var_desc(vindex)
            my_pioVar%gtype=v3dvar
            IF (KIND(OCEAN(ng)%b_v).eq.8) THEN
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
     &                        OCEAN(ng) % b_v(:,:,:,Tindex),            &
     &                        checksum = Fhash)
# else
     &                        OCEAN(ng) % b_v(:,:,:,Tindex))
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
# ifdef DISTRIBUTE
            CALL mp_exchange3d (ng, MyRank, IDmod, 1,                   &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          OCEAN(ng) % b_v(:,:,:,Tindex))
# endif
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
!  Read in tracer type variables normalization factor.
!
        DO itrc=1,NT(ng)
          IF (get_var(idTvar(itrc)).and.                                &
     &        ((model.eq.14).or.(model.eq.15))) THEN
            foundit=find_string(var_name, n_var,                        &
     &                          TRIM(Vname(1,idTvar(itrc))), vindex)
            IF (foundit) THEN
              my_pioVar%vd=var_desc(vindex)
              my_pioVar%gtype=r3dvar
              IF (KIND(OCEAN(ng)%b_t).eq.8) THEN
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
     &                          OCEAN(ng) % b_t(:,:,:,Tindex,itrc),     &
     &                          checksum = Fhash)
# else
     &                          OCEAN(ng) % b_t(:,:,:,Tindex,itrc))
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
# ifdef DISTRIBUTE
              CALL mp_exchange3d (ng, MyRank, IDmod, 1,                 &
     &                            LBi, UBi, LBj, UBj, 1, N(ng),         &
     &                            NghostPoints,                         &
     &                            EWperiodic(ng), NSperiodic(ng),       &
     &                            OCEAN(ng) % b_t(:,:,:,Tindex,itrc))
# endif
            ELSE
              IF (Master) THEN
                WRITE (stdout,80) string, TRIM(Vname(1,idTvar(itrc))),  &
     &                          TRIM(ncname)
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
#ifdef ADJUST_BOUNDARY
!
!  Read in free-surface open boundaries normalization factor.
!
        IF (get_var(idSbry(isFsur)).and.(model.eq.16).and.              &
     &      ANY(Lobc(:,isFsur,ng))) THEN
          CALL pio_netcdf_get_fvar (ng, IDmod, ncname,                  &
     &                              Vname(1,idSbry(isFsur)),            &
     &                              BOUNDARY(ng) % b_zeta_obc(LBij:,:), &
     &                              pioFile = pioFile,                  &
     &                              start = (/1,1,InpRec/),             &
     &                              total = (/IorJ,4,1/),               &
     &                              min_val = Fmin,                     &
     &                              max_val = Fmax)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          IF (Master) THEN
            WRITE (stdout,75) TRIM(Vname(1,idSbry(isFsur))),            &
     &                        Fmin, Fmax
          END IF
        END IF
!
!  Read in 2D U-momentum component open boundaries normalization factor.
!
        IF (get_var(idSbry(isUbar)).and.(model.eq.16).and.              &
     &        ANY(Lobc(:,isUbar,ng))) THEN
          CALL pio_netcdf_get_fvar (ng, IDmod, ncname,                  &
     &                              Vname(1,idSbry(isUbar)),            &
     &                              BOUNDARY(ng) % b_ubar_obc(LBij:,:), &
     &                              pioFile = pioFile,                  &
     &                              start = (/1,1,InpRec/),             &
     &                              total = (/IorJ,4,1/),               &
     &                              min_val = Fmin,                     &
     &                              max_val = Fmax)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          IF (Master) THEN
            WRITE (stdout,75) TRIM(Vname(1,idSbry(isUbar))),            &
     &                        Fmin, Fmax
          END IF
        END IF
!
!  Read in 2D V-momentum component open boundaries normalization factor.
!
        IF (get_var(idSbry(isVbar)).and.(model.eq.16).and.              &
     &      ANY(Lobc(:,isVbar,ng))) THEN
          CALL pio_netcdf_get_fvar (ng, IDmod, ncname,                  &
     &                              Vname(1,idSbry(isVbar)),            &
     &                              BOUNDARY(ng) % b_vbar_obc(LBij:,:), &
     &                              pioFile = pioFile,                  &
     &                              start = (/1,1,InpRec/),             &
     &                              total = (/IorJ,4,1/),               &
     &                              min_val = Fmin,                     &
     &                              max_val = Fmax)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          IF (Master) THEN
            WRITE (stdout,75) TRIM(Vname(1,idSbry(isVbar))),            &
     &                        Fmin, Fmax
          END IF
        END IF

# ifdef SOLVE3D
!
!  Read in 3D U-momentum component open boundaries normalization factor.
!
        IF (get_var(idSbry(isUvel)).and.(model.eq.16).and.              &
     &      ANY(Lobc(:,isUvel,ng))) THEN
          CALL pio_netcdf_get_fvar (ng, IDmod, ncname,                  &
     &                              Vname(1,idSbry(isUvel)),            &
     &                              BOUNDARY(ng) % b_u_obc(LBij:,:,:),  &
     &                              pioFile = pioFile,                  &
     &                              start = (/1,1,1,InpRec/),           &
     &                              total = (/IorJ,N(ng),4,1/),         &
     &                              min_val = Fmin,                     &
     &                              max_val = Fmax)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          IF (Master) THEN
            WRITE (stdout,75) TRIM(Vname(1,idSbry(isUvel))),            &
     &                        Fmin, Fmax
          END IF
        END IF
!
!  Read in 3D V-momentum component open boundaries normalization factor.
!
        IF (get_var(idSbry(isVvel)).and.(model.eq.16).and.              &
     &      ANY(Lobc(:,isVvel,ng))) THEN
          CALL pio_netcdf_get_fvar (ng, IDmod, ncname,                  &
     &                              Vname(1,idSbry(isVvel)),            &
     &                              BOUNDARY(ng) % b_v_obc(LBij:,:,:),  &
     &                              pioFile = pioFile,                  &
     &                              start = (/1,1,1,InpRec/),           &
     &                              total = (/IorJ,N(ng),4,1/),         &
     &                              min_val = Fmin,                     &
     &                              max_val = Fmax)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          IF (Master) THEN
            WRITE (stdout,75) TRIM(Vname(1,idSbry(isVvel))),            &
     &                        Fmin, Fmax
          END IF
        END IF
!
!  Read in 3D tracers open boundaries normalization factor.
!
        DO itrc=1,NT(ng)
          IF (get_var(idSbry(isTvar(itrc))).and.(model.eq.16).and.      &
     &        ANY(Lobc(:,isTvar(itrc),ng))) THEN
            CALL pio_netcdf_get_fvar (ng, IDmod, ncname,                &
     &                                Vname(1,idSbry(isTvar(itrc))),    &
     &                                BOUNDARY(ng) % b_t_obc(LBij:,:,:, &
     &                                                       itrc),     &
     &                                pioFile = pioFile,                &
     &                                start =(/1,1,1,InpRec/),          &
     &                                total =(/IorJ,N(ng),4,1/),        &
     &                                min_val = Fmin,                   &
     &                                max_val = Fmax)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
            IF (Master) THEN
              WRITE (stdout,75) TRIM(Vname(1,idSbry(isTvar(itrc)))),    &
     &                          Fmin, Fmax
            END IF
          END IF
        END DO
# endif
#endif
#ifdef ADJUST_WSTRESS
!
!  Read in surface U-momentum stress normalization factors.
!
        IF (get_var(idUsms).and.(model.eq.17)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idUsms)),   &
     &                        vindex)
          IF (foundit) THEN
            my_pioVar%vd=var_desc(vindex)
            my_pioVar%gtype=u2dvar
            IF (KIND(FORCES(ng)%b_sustr).eq.8) THEN
              my_pioVar%dkind=PIO_double
              ioDesc => ioDesc_dp_u2dvar(ng)
            ELSE
              my_pioVar%dkind=PIO_real
              ioDesc => ioDesc_sp_u2dvar(ng)
            END IF
!
            status=nf_fread2d(ng, IDmod, ncname, pioFile,               &
     &                        Vname(1,idUsms), my_pioVar,               &
     &                        InpRec, ioDesc, Vsize,                    &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
# ifdef MASKING
     &                        GRID(ng) % umask,                         &
# endif
# ifdef CHECKSUM
     &                        FORCES(ng) % b_sustr,                     &
     &                        checksum = Fhash)
# else
     &                        FORCES(ng) % b_sustr)
# endif
            IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
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
                WRITE (stdout,70) TRIM(Vname(2,idUsms)), Fmin, Fmax,    &
     &                            Fhash
# else
                WRITE (stdout,70) TRIM(Vname(2,idUsms)), Fmin, Fmax
# endif
              END IF
            END IF
# ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, MyRank, IDmod, 1,                   &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          FORCES(ng) % b_sustr)
# endif
          ELSE
            IF (Master) THEN
              WRITE (stdout,80) string, TRIM(Vname(1,idUsms)),          &
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
!  Read in surface V-momentum stress normalization factors.
!
        IF (get_var(idVsms).and.(model.eq.17)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idVsms)),   &
     &                        vindex)
          IF (foundit) THEN
            my_pioVar%vd=var_desc(vindex)
            my_pioVar%gtype=v2dvar
            IF (KIND(FORCES(ng)%b_svstr).eq.8) THEN
              my_pioVar%dkind=PIO_double
              ioDesc => ioDesc_dp_v2dvar(ng)
            ELSE
              my_pioVar%dkind=PIO_real
              ioDesc => ioDesc_sp_v2dvar(ng)
            END IF
!
            status=nf_fread2d(ng, IDmod, ncname, pioFile,               &
     &                        Vname(1,idVsms), my_pioVar,               &
     &                        InpRec, ioDesc, Vsize,                    &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
# ifdef MASKING
     &                        GRID(ng) % vmask,                         &
# endif
# ifdef CHECKSUM
     &                        FORCES(ng) % b_svstr,                     &
     &                        checksum = Fhash)
# else
     &                        FORCES(ng) % b_svstr)
# endif
            IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
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
                WRITE (stdout,70) TRIM(Vname(2,idVsms)), Fmin, Fmax,    &
     &                            Fhash
# else
                WRITE (stdout,70) TRIM(Vname(2,idVsms)), Fmin, Fmax
# endif
              END IF
            END IF
# ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, MyRank, IDmod, 1,                   &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          FORCES(ng) % b_svstr)
# endif
          ELSE
            IF (Master) THEN
              WRITE (stdout,80) string, TRIM(Vname(1,idVsms)),          &
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
#if defined ADJUST_STFLUX && defined SOLVE3D
!
!  Read in surface tracer flux normalization factors.
!
        DO itrc=1,NT(ng)
          IF (get_var(idTsur(itrc)).and.(model.eq.17).and.              &
     &        Lstflux(itrc,ng)) THEN
            foundit=find_string(var_name, n_var,                        &
     &                          TRIM(Vname(1,idTsur(itrc))), vindex)
            IF (foundit) THEN
              my_pioVar%vd=var_desc(vindex)
              my_pioVar%gtype=r2dvar
              IF (KIND(FORCES(ng)%b_stflx).eq.8) THEN
                my_pioVar%dkind=PIO_double
                ioDesc => ioDesc_dp_r2dvar(ng)
              ELSE
                my_pioVar%dkind=PIO_real
                ioDesc => ioDesc_sp_r2dvar(ng)
              END IF
!
              status=nf_fread2d(ng, IDmod, ncname, pioFile,             &
     &                          Vname(1,idTsur(itrc)), my_pioVar,       &
     &                          InpRec, ioDesc, Vsize,                  &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Fscl, Fmin, Fmax,                       &
# ifdef MASKING
     &                          GRID(ng) % rmask,                       &
# endif
# ifdef CHECKSUM
     &                          FORCES(ng) % b_stflx(:,:,itrc),         &
     &                          checksum = Fhash)
# else
     &                          FORCES(ng) % b_stflx(:,:,itrc))
# endif
              IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
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
                  WRITE (stdout,70) TRIM(Vname(2,idTsur(itrc))),        &
     &                              Fmin, Fmax, Fhash
# else
                  WRITE (stdout,70) TRIM(Vname(2,idTsur(itrc))),        &
     &                              Fmin, Fmax
# endif
                END IF
              END IF
# ifdef DISTRIBUTE
              CALL mp_exchange2d (ng, MyRank, IDmod, 1,                 &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            NghostPoints,                         &
     &                            EWperiodic(ng), NSperiodic(ng),       &
     &                            FORCES(ng) % b_stflx(:,:,itrc))
# endif
            ELSE
              IF (Master) THEN
                WRITE (stdout,80) string, TRIM(Vname(1,idTsur(itrc))),  &
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
