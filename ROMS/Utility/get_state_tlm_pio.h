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
** It process the reading of ROMS TLM or RPM state from input NetCDF  **
** files.                                                             **
**                                                                    **
** It is used to read the perturbation and finite-amplitude tangent   **
** linear kernels state variables. Data is loaded into "tl_" arrays.  **
**                                                                    **
************************************************************************
*/

#if defined ADJUST_BOUNDARY || \
    defined ADJUST_WSTRESS  || defined ADJUST_STFLUX
        IF (inner.eq.0.and.model.eq.iRPM) THEN
          get_adjust=.FALSE.
        ELSE
          get_adjust=.TRUE.
        END IF
#endif
!
!  Read in tangent linear free-surface (m).
!
        IF (get_var(idFsur)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idFsur)),   &
     &                        vindex)
          IF (foundit) THEN
            my_pioVar%vd=var_desc(vindex)
            my_pioVar%gtype=r2dvar
            IF (KIND(OCEAN(ng)%tl_zeta).eq.8) THEN
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
     &                        OCEAN(ng) % tl_zeta(:,:,Tindex),          &
     &                        checksum = Fhash)
#else
     &                        OCEAN(ng) % tl_zeta(:,:,Tindex))
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

#ifdef ADJUST_BOUNDARY
!
!  Read in free-surface open boundaries adjustments.
!
        IF (get_var(idSbry(isFsur)).and.get_adjust.and.                 &
     &      ANY(Lobc(:,isFsur,ng))) THEN
          ifield=idSbry(isFsur)
          foundit=find_string(var_name, n_var, TRIM(Vname(1,ifield)),   &
     &                        vindex)
          IF (foundit) THEN
            my_pioVar%vd=var_desc(vindex)
            my_pioVar%gtype=r2dobc
            IF (KIND(BOUNDARY(ng)%tl_zeta_obc).eq.8) THEN
              my_pioVar%dkind=PIO_double
            ELSE
              my_pioVar%dkind=PIO_real
            END IF
!
            status=nf_fread2d_bry(ng, IDmod, ncname, pioFile,           &
     &                            Vname(1,ifield), my_pioVar,           &
     &                            InpRec,                               &
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
            IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
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
            IF (FoundError(exit_flag, PIO_noerr,                        &
     &                     __LINE__, MyFile)) THEN
              RETURN
            END IF
          END IF
        END IF
#endif
!
!  Read in tangent linear 2D U-momentum component (m/s).
!
        IF (get_var(idUbar)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idUbar)),   &
     &                        vindex)
          IF (foundit) THEN
            my_pioVar%vd=var_desc(vindex)
            my_pioVar%gtype=u2dvar
            IF (KIND(OCEAN(ng)%tl_ubar).eq.8) THEN
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
     &                        OCEAN(ng) % tl_ubar(:,:,Tindex),          &
     &                        checksum = Fhash)
#else
     &                        OCEAN(ng) % tl_ubar(:,:,Tindex))
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

#ifdef ADJUST_BOUNDARY
!
!  Read in 2D U-momentum component open boundaries adjustments.
!
        IF (get_var(idSbry(isUbar)).and.get_adjust.and.                 &
     &        ANY(Lobc(:,isUbar,ng))) THEN
          ifield=idSbry(isUbar)
          foundit=find_string(var_name, n_var, TRIM(Vname(1,ifield)),   &
     &                        vindex)
          IF (foundit) THEN
            my_pioVar%vd=var_desc(vindex)
            my_pioVar%gtype=u2dobc
            IF (KIND(BOUNDARY(ng)%tl_ubar_obc).eq.8) THEN
              my_pioVar%dkind=PIO_double
            ELSE
              my_pioVar%dkind=PIO_real
            END IF
!
            status=nf_fread2d_bry(ng, IDmod, ncname, pioFile,           &
     &                            Vname(1,ifield), my_pioVar,           &
     &                            InpRec,                               &
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
            IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
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
            IF (FoundError(exit_flag, PIO_noerr,                        &
     &                     __LINE__, MyFile)) THEN
              RETURN
            END IF
          END IF
        END IF
#endif
!
!  Read in tangent linear 2D V-momentum component.
!
        IF (get_var(idVbar)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idVbar)),   &
     &                        vindex)
          IF (foundit) THEN
            my_pioVar%vd=var_desc(vindex)
            my_pioVar%gtype=v2dvar
            IF (KIND(OCEAN(ng)%tl_vbar).eq.8) THEN
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
     &                        OCEAN(ng) % tl_vbar(:,:,Tindex),          &
     &                        checksum = Fhash)
#else
     &                        OCEAN(ng) % tl_vbar(:,:,Tindex))
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

#ifdef ADJUST_BOUNDARY
!
!  Read in 2D V-momentum component open boundaries adjustments.
!
        IF (get_var(idSbry(isVbar)).and.get_adjust.and.                 &
     &      ANY(Lobc(:,isVbar,ng))) THEN
          ifield=idSbry(isVbar)
          foundit=find_string(var_name, n_var, TRIM(Vname(1,ifield)),   &
     &                        vindex)
          IF (foundit) THEN
            my_pioVar%vd=var_desc(vindex)
            my_pioVar%gtype=v2dobc
            IF (KIND(BOUNDARY(ng)%tl_vbar_obc).eq.8) THEN
              my_pioVar%dkind=PIO_double
            ELSE
              my_pioVar%dkind=PIO_real
            END IF
!
            status=nf_fread2d_bry(ng, IDmod, ncname, pioFile,           &
     &                            Vname(1,ifield), my_pioVar,           &
     &                            InpRec,                               &
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
            IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
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
            IF (FoundError(exit_flag, PIO_noerr,                        &
     &                     __LINE__, MyFile)) THEN
              RETURN
            END IF
          END IF
        END IF
#endif
#ifdef ADJUST_WSTRESS
!
!  Read in tangent linear surface U-momentum stress.
!
        IF (get_var(idUsms).and.get_adjust) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idUsms)),   &
     &                        vindex)
          IF (foundit) THEN
            scale=1.0_dp
            my_pioVar%vd=var_desc(vindex)
            my_pioVar%gtype=u2dvar
            IF (KIND(FORCES(ng)%tl_ustr).eq.8) THEN
              my_pioVar%dkind=PIO_double
              ioDesc => ioDesc_dp_u2dfrc(ng)
            ELSE
              my_pioVar%dkind=PIO_real
              ioDesc => ioDesc_sp_u2dfrc(ng)
            END IF
!
            status=nf_fread3d(ng, IDmod, ncname, pioFile,               &
     &                        Vname(1,idUsms), my_pioVar,               &
     &                        InpRec, ioDesc, Vsize,                    &
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
            IF (FoundError(exit_flag, PIO_noerr,                        &
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
     &                        vindex)
          IF (foundit) THEN
            scale=1.0_dp
            my_pioVar%vd=var_desc(vindex)
            my_pioVar%gtype=v2dvar
            IF (KIND(FORCES(ng)%tl_vstr).eq.8) THEN
              my_pioVar%dkind=PIO_double
              ioDesc => ioDesc_dp_v2dfrc(ng)
            ELSE
              my_pioVar%dkind=PIO_real
              ioDesc => ioDesc_sp_v2dfrc(ng)
            END IF
!
            status=nf_fread3d(ng, IDmod, ncname, pioFile,               &
     &                        Vname(1,idVsms), my_pioVar,               &
     &                        InpRec, ioDesc, Vsize,                    &
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
            IF (FoundError(exit_flag, PIO_noerr,                        &
     &                     __LINE__, MyFile)) THEN
              RETURN
            END IF
          END IF
        END IF
#endif
#ifdef SOLVE3D
!
!  Read in tangent linear 3D U-momentum component.

!
        IF (get_var(idUvel)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idUvel)),   &
     &                        vindex)
          IF (foundit) THEN
            my_pioVar%vd=var_desc(vindex)
            my_pioVar%gtype=u3dvar
            IF (KIND(OCEAN(ng)%tl_u).eq.8) THEN
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
     &                        OCEAN(ng) % tl_u(:,:,:,Tindex),           &
     &                        checksum = Fhash)
# else
     &                        OCEAN(ng) % tl_u(:,:,:,Tindex))
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

# ifdef ADJUST_BOUNDARY
!
!  Read in 3D U-momentum component open boundaries adjustments.
!
        IF (get_var(idSbry(isUvel)).and.get_adjust.and.                 &
     &      ANY(Lobc(:,isUvel,ng))) THEN
          ifield=idSbry(isUvel)
          foundit=find_string(var_name, n_var, TRIM(Vname(1,ifield)),   &
     &                        vindex)
          IF (foundit) THEN
            my_pioVar%vd=var_desc(vindex)
            my_pioVar%gtype=u3dobc
            IF (KIND(BOUNDARY(ng)%tl_u_obc).eq.8) THEN
              my_pioVar%dkind=PIO_double
            ELSE
              my_pioVar%dkind=PIO_real
            END IF
!
            status=nf_fread3d_bry(ng, IDmod, ncname, pioFile,           &
     &                            Vname(1,ifield), my_pioVar,           &
     &                            InpRec,                               &
     &                            LBij, UBij, 1, N(ng), Nbrec(ng),      &
     &                            Fscl, Fmin, Fmax,                     &
# ifdef CHECKSUM
     &                            BOUNDARY(ng) % tl_u_obc(:,:,:,:,      &
     &                                                    Tindex),      &
     &                            checksum = Fhash)
# else
     &                            BOUNDARY(ng) % tl_u_obc(:,:,:,:,      &
     &                                                    Tindex))
# endif
            IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
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
            IF (FoundError(exit_flag, PIO_noerr,                        &
     &                     __LINE__, MyFile)) THEN
              RETURN
            END IF
          END IF
        END IF
# endif
!
!  Read in tangent linear 3D V-momentum component.
!
        IF (get_var(idVvel)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idVvel)),   &
     &                        vindex)
          IF (foundit) THEN
            my_pioVar%vd=var_desc(vindex)
            my_pioVar%gtype=v3dvar
            IF (KIND(OCEAN(ng)%tl_v).eq.8) THEN
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
     &                        OCEAN(ng) % tl_v(:,:,:,Tindex),           &
     &                        checksum = Fhash)
# else
     &                        OCEAN(ng) % tl_v(:,:,:,Tindex))
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

# ifdef ADJUST_BOUNDARY
!
!  Read in 3D V-momentum component open boundaries adjustments.
!
        IF (get_var(idSbry(isVvel)).and.get_adjust.and.                 &
     &      ANY(Lobc(:,isVvel,ng))) THEN
          ifield=idSbry(isVvel)
          foundit=find_string(var_name, n_var, TRIM(Vname(1,ifield)),   &
     &                        vindex)
          IF (foundit) THEN
            my_pioVar%vd=var_desc(vindex)
            my_pioVar%gtype=v3dobc
            IF (KIND(BOUNDARY(ng)%tl_v_obc).eq.8) THEN
              my_pioVar%dkind=PIO_double
            ELSE
              my_pioVar%dkind=PIO_real
            END IF
!
            status=nf_fread3d_bry(ng, IDmod, ncname, pioFile,           &
     &                            Vname(1,ifield), my_pioVar,           &
     &                            InpRec,                               &
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
            IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
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
            IF (FoundError(exit_flag, PIO_noerr,                        &
     &                     __LINE__, MyFile)) THEN
              RETURN
            END IF
          END IF
        END IF
# endif
!
!  Read in tangent linear tracer type variables.
!
        DO itrc=1,NT(ng)
          IF (get_var(idTvar(itrc))) THEN
            foundit=find_string(var_name, n_var,                        &
     &                          TRIM(Vname(1,idTvar(itrc))), vindex)
            IF (foundit) THEN
              my_pioVar%vd=var_desc(vindex)
              my_pioVar%gtype=r3dvar
              IF (KIND(OCEAN(ng)%tl_t).eq.8) THEN
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
     &                          OCEAN(ng) % tl_t(:,:,:,Tindex,itrc),    &
     &                          checksum = Fhash)
# else
     &                          OCEAN(ng) % tl_t(:,:,:,Tindex,itrc))
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
#  ifdef CHECKSUM
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

# ifdef ADJUST_BOUNDARY
!
!  Read in 3D tracers open boundaries adjustments.
!
        DO itrc=1,NT(ng)
          IF (get_var(idSbry(isTvar(itrc))).and.get_adjust.and.         &
     &        ANY(Lobc(:,isTvar(itrc),ng))) THEN
            ifield=idSbry(isTvar(itrc))
            foundit=find_string(var_name, n_var, TRIM(Vname(1,ifield)), &
     &                          vindex)
            IF (foundit) THEN
              my_pioVar%vd=var_desc(vindex)
              my_pioVar%gtype=r3dobc
              IF (KIND(BOUNDARY(ng)%tl_t_obc).eq.8) THEN
                my_pioVar%dkind=PIO_double
              ELSE
                my_pioVar%dkind=PIO_real
              END IF
!
              status=nf_fread3d_bry(ng, IDmod, ncname, pioFile,         &
     &                              Vname(1,ifield), my_pioVar,         &
     &                              InpRec,                             &
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
              IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
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
     &                              Fhash
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
              IF (FoundError(exit_flag, PIO_noerr,                      &
     &                       __LINE__, MyFile)) THEN
                RETURN
              END IF
            END IF
          END IF
        END DO
# endif
# ifdef ADJUST_STFLUX
!
!  Read in tangent linear surface tracers flux.
!
        DO itrc=1,NT(ng)
          IF (get_var(idTsur(itrc)).and.get_adjust.and.                 &
     &        Lstflux(itrc,ng)) THEN
            foundit=find_string(var_name, n_var,                        &
     &                          TRIM(Vname(1,idTsur(itrc))), vindex)
            IF (foundit) THEN
              scale=1.0_dp
              my_pioVar%vd=var_desc(vindex)
              my_pioVar%gtype=r2dvar
              IF (KIND(FORCES(ng)%tl_tflux).eq.8) THEN
                my_pioVar%dkind=PIO_double
                ioDesc => ioDesc_dp_r2dfrc(ng)
              ELSE
                my_pioVar%dkind=PIO_real
                ioDesc => ioDesc_sp_r2dfrc(ng)
              END IF
!
              status=nf_fread3d(ng, IDmod, ncname, pioFile,             &
     &                          Vname(1,idTsur(itrc)), my_pioVar,       &
     &                          InpRec, ioDesc, Vsize,                  &
     &                          LBi, UBi, LBj, UBj, 1, Nfrec(ng),       &
     &                          scale, Fmin, Fmax,                      &
#  ifdef MASKING
     &                          GRID(ng) % rmask,                       &
#  endif
# ifdef CHECKSUM
     &                          FORCES(ng)% tl_tflux(:,:,:,             &
     &                                               Tindex,itrc),      &
     &                          checksum = Fhash)
# else
     &                          FORCES(ng)% tl_tflux(:,:,:,             &
     &                                               Tindex,itrc))
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
              IF (FoundError(exit_flag, PIO_noerr,                      &
     &                       __LINE__, MyFile)) THEN
                RETURN
              END IF
            END IF
          END IF
        END DO
# endif
# ifdef SEDIMENT
!
!  Read in tangent linear sediment fraction of each size class in each
!  bed layer.
!
        DO i=1,NST
          IF (get_var(idfrac(i))) THEN
            foundit=find_string(var_name, n_var,                        &
     &                          TRIM(Vname(1,idfrac(i))), vindex)
            IF (foundit) THEN
              my_pioVar%vd=var_desc(vindex)
              my_pioVar%gtype=b3dvar
              IF (KIND(OCEAN(ng)%tl_bed_frac).eq.8) THEN
                my_pioVar%dkind=PIO_double
                ioDesc => ioDesc_dp_b3dvar(ng)
              ELSE
                my_pioVar%dkind=PIO_real
                ioDesc => ioDesc_sp_b3dvar(ng)
              END IF
!
              status=nf_fread3d(ng, IDmod, ncname, pioFile,             &
     &                          Vname(1,idfrac(i)), my_pioVar,          &
     &                          InpRec, ioDesc, Vsize,                  &
     &                          LBi, UBi, LBj, UBj, 1, Nbed,            &
     &                          Fscl, Fmin, Fmax,                       &
#  ifdef MASKING
     &                          GRID(ng) % rmask,                       &
#  endif
# ifdef CHECKSUM
     &                          SEDBED(ng) % tl_bed_frac(:,:,:,i),      &
     &                          checksum = Fhash)
# else
     &                          SEDBED(ng) % tl_bed_frac(:,:,:,i))
# endif
              IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
                IF (Master) THEN
                  WRITE (stdout,60) string, TRIM(Vname(1,idfrac(i))),   &
     &                              InpRec, TRIM(ncname)
                END IF
                exit_flag=2
                ioerror=status
                RETURN
              ELSE
                IF (Master) THEN
# ifdef CHECKSUM
                  WRITE (stdout,70) TRIM(Vname(2,idfrac(i))),           &
     &                              Fmin, Fmax, Fhash
# else
                  WRITE (stdout,70) TRIM(Vname(2,idfrac(i))),           &
     &                              Fmin, Fmax
# endif
                END IF
              END IF
            ELSE
              IF (Master) THEN
                WRITE (stdout,80) string, TRIM(Vname(1,idfrac(i))),     &
     &                            TRIM(ncname)
              END IF
              exit_flag=4
              IF (FoundError(exit_flag, PIO_noerr,                      &
     &                       __LINE__, MyFile)) THEN
                RETURN
              END IF
            END IF
          END IF
!
!  Read in tangent linear sediment mass of each size class in each
!  bed layer.
!
          IF (get_var(idBmas(i))) THEN
            foundit=find_string(var_name, n_var,                        &
     &                          TRIM(Vname(1,idBmas(i))), vindex)
            IF (foundit) THEN
              my_pioVar%vd=var_desc(vindex)
              my_pioVar%gtype=b3dvar
              IF (KIND(OCEAN(ng)%tl_bed_mass).eq.8) THEN
                my_pioVar%dkind=PIO_double
                ioDesc => ioDesc_dp_b3dvar(ng)
              ELSE
                my_pioVar%dkind=PIO_real
                ioDesc => ioDesc_sp_b3dvar(ng)
              END IF
!
              status=nf_fread3d(ng, IDmod, ncname, pioFile,             &
     &                          Vname(1,idBmas(i)), my_pioVar,          &
     &                          InpRec, ioDesc, Vsize,                  &
     &                          LBi, UBi, LBj, UBj, 1, Nbed,            &
     &                          Fscl, Fmin, Fmax,                       &
#  ifdef MASKING
     &                          GRID(ng) % rmask,                       &
#  endif
#  ifdef CHECKSUM
     &                          SEDBED(ng) % tl_bed_mass(:,:,:,         &
     &                                                   Tindex,i),     &
     &                          checksum = Fhash)
#  else
     &                          SEDBED(ng) % tl_bed_mass(:,:,:,         &
     &                                                   Tindex,i))
#  endif
              IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
                IF (Master) THEN
                  WRITE (stdout,60) string, TRIM(Vname(1,idBmas(i))),   &
     &                              InpRec, TRIM(ncname)
                END IF
                exit_flag=2
                ioerror=status
                RETURN
              ELSE
                IF (Master) THEN
#  ifdef CHECKSUM
                  WRITE (stdout,70) TRIM(Vname(2,idBmas(i))),           &
     &                              Fmin, Fmax, Fhash
#  else
                  WRITE (stdout,70) TRIM(Vname(2,idBmas(i))),           &
     &                              Fmin, Fmax
#  endif
                END IF
              END IF
            ELSE
              IF (Master) THEN
                WRITE (stdout,80) string, TRIM(Vname(1,idBmas(i))),     &
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
!
!  Read in tangent linear sediment properties in each bed layer.
!
        DO i=1,MBEDP
          IF (get_var(idSbed(i))) THEN
            foundit=find_string(var_name, n_var,                        &
     &                          TRIM(Vname(1,idSbed(i))), vindex)
            IF (foundit) THEN
              my_pioVar%vd=var_desc(vindex)
              my_pioVar%gtype=b3dvar
              IF (KIND(OCEAN(ng)%tl_bed).eq.8) THEN
                my_pioVar%dkind=PIO_double
                ioDesc => ioDesc_dp_b3dvar(ng)
              ELSE
                my_pioVar%dkind=PIO_real
                ioDesc => ioDesc_sp_b3dvar(ng)
              END IF
!
              status=nf_fread3d(ng, IDmod, ncname, pioFile,             &
     &                          Vname(1,idSbed(i)), my_pioVar,          &
     &                          InpRec, ioDesc, Vsize,                  &
     &                          LBi, UBi, LBj, UBj, 1, Nbed,            &
     &                          Fscl, Fmin, Fmax,                       &
#  ifdef MASKING
     &                          GRID(ng) % rmask,                       &
#  endif
#  ifdef CHECKSUM
     &                          SEDBED(ng) % tl_bed(:,:,:,i),           &
     &                          checksum = Fhash)
#  else
     &                          SEDBED(ng) % tl_bed(:,:,:,i))
#  endif
              IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
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
              IF (FoundError(exit_flag, PIO_noerr,                      &
     &                       __LINE__, MyFile)) THEN
                RETURN
              END IF
            END IF
          END IF
        END DO

#  ifdef BEDLOAD
!
!  Read in tangent linear sediment fraction of bed load.
!
        DO i=1,NST
          IF (get_var(idUbld(i))) THEN
            foundit=find_string(var_name, n_var,                        &
     &                          TRIM(Vname(1,idUbld(i))), vindex)
            IF (foundit) THEN
              my_pioVar%vd=var_desc(vindex)
              my_pioVar%gtype=u2dvar
              IF (KIND(OCEAN(ng)%tl_bedldu).eq.8) THEN
                my_pioVar%dkind=PIO_double
                ioDesc => ioDesc_dp_u2dvar(ng)
              ELSE
                my_pioVar%dkind=PIO_real
                ioDesc => ioDesc_sp_u2dvar(ng)
              END IF
!
              status=nf_fread2d(ng, IDmod, ncname, pioFile,             &
     &                          Vname(1,idUbld(i)), my_pioVar,          &
     &                          InpRec, ioDesc, Vsize,                  &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Fscl, Fmin, Fmax,                       &
#   ifdef MASKING
     &                          GRID(ng) % umask,                       &
#   endif
#   ifdef CHECKSUM
     &                          SEDBED(ng) % tl_bedldu(:,:,i),          &
     &                          checksum = Fhash)
#   else
     &                          SEDBED(ng) % tl_bedldu(:,:,i))
#   endif
              IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
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
              IF (FoundError(exit_flag, PIO_noerr,                      &
     &                       __LINE__, MyFile)) THEN
                RETURN
              END IF
            END IF
          END IF
!
          IF (get_var(idVbld(i))) THEN
            foundit=find_string(var_name, n_var,                        &
     &                          TRIM(Vname(1,idVbld(i))), vindex)
            IF (foundit) THEN
              my_pioVar%vd=var_desc(vindex)
              my_pioVar%gtype=v2dvar
              IF (KIND(OCEAN(ng)%tl_bedldv).eq.8) THEN
                my_pioVar%dkind=PIO_double
                ioDesc => ioDesc_dp_v2dvar(ng)
              ELSE
                my_pioVar%dkind=PIO_real
                ioDesc => ioDesc_sp_v2dvar(ng)
              END IF
!
              status=nf_fread2d(ng, IDmod, ncname, pioFile,             &
     &                          Vname(1,idVbld(i)), my_pioVar,          &
     &                          InpRec, ioDesc, Vsize,                  &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Fscl, Fmin, Fmax,                       &
#   ifdef MASKING
     &                          GRID(ng) % vmask,                       &
#   endif
#   ifdef CHECKSUM
     &                          SEDBED(ng) % tl_bedldv(:,:,i),          &
     &                          checksum = Fhash)
#   else
     &                          SEDBED(ng) % tl_bedldv(:,:,i))
#   endif
              IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
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
              IF (FoundError(exit_flag, PIO_noerr,                      &
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
!  Read in tangent linear sediment properties in exposed bed layer.
!
        DO i=1,MBOTP
          IF (get_var(idBott(i)).and.have_var(idBott(i))) THEN
            foundit=find_string(var_name, n_var,                        &
     &                          TRIM(Vname(1,idBott(i))), vindex)
            IF (foundit) THEN
              my_pioVar%vd=var_desc(vindex)
              my_pioVar%gtype=r2dvar
              IF (KIND(OCEAN(ng)%tl_bottom).eq.8) THEN
                my_pioVar%dkind=PIO_double
                ioDesc => ioDesc_dp_r2dvar(ng)
              ELSE
                my_pioVar%dkind=PIO_real
                ioDesc => ioDesc_sp_r2dvar(ng)
              END IF
!
              status=nf_fread2d(ng, IDmod, ncname, pioFile,             &
     &                          Vname(1,idBott(i)), my_pioVar,          &
     &                          InpRec, ioDesc, Vsize,                  &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Fscl, Fmin, Fmax,                       &
#  ifdef MASKING
     &                          GRID(ng) % rmask,                       &
#  endif
#  ifdef CHECKSUM
     &                          SEDBED(ng) % tl_bottom(:,:,i),          &
     &                          checksum = Fhash)
#  else
     &                          SEDBED(ng) % tl_bottom(:,:,i))
#  endif
              IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
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
              IF (FoundError(exit_flag, PIO_noerr,                      &
     &                       __LINE__, MyFile)) THEN
                RETURN
              END IF
            END IF
          END IF
        END DO
# endif
#endif
