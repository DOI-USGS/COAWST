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
     &                        varid)
          IF (foundit) THEN
            gtype=var_flag(varid)*r2dvar
            status=nf_fread2d(ng, IDmod, ncname, ncINPid,               &
     &                        Vname(1,idFsur), varid,                   &
     &                        InpRec, gtype, Vsize,                     &
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
#endif
!
!  Read in tangent linear 2D U-momentum component (m/s).
!
        IF (get_var(idUbar)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idUbar)),   &
     &                        varid)
          IF (foundit) THEN
            gtype=var_flag(varid)*u2dvar
            status=nf_fread2d(ng, IDmod, ncname, ncINPid,               &
     &                        Vname(1,idUbar), varid,                   &
     &                        InpRec, gtype, Vsize,                     &
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

#ifdef ADJUST_BOUNDARY
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
#endif
!
!  Read in tangent linear 2D V-momentum component.
!
        IF (get_var(idVbar)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idVbar)),   &
     &                        varid)
          IF (foundit) THEN
            gtype=var_flag(varid)*v2dvar
            status=nf_fread2d(ng, IDmod, ncname, ncINPid,               &
     &                        Vname(1,idVbar), varid,                   &
     &                        InpRec, gtype, Vsize,                     &
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

#ifdef ADJUST_BOUNDARY
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
#ifdef SOLVE3D
!
!  Read in tangent linear 3D U-momentum component.

!
        IF (get_var(idUvel)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idUvel)),   &
     &                        varid)
          IF (foundit) THEN
            gtype=var_flag(varid)*u3dvar
            status=nf_fread3d(ng, IDmod, ncname, ncINPid,               &
     &                        Vname(1,idUvel), varid,                   &
     &                        InpRec, gtype, Vsize,                     &
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

# ifdef ADJUST_BOUNDARY
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
# ifdef CHECKSUM
     &                            BOUNDARY(ng) % tl_u_obc(:,:,:,:,      &
     &                                                    Tindex),      &
     &                            checksum = Fhash)
# else
     &                            BOUNDARY(ng) % tl_u_obc(:,:,:,:,      &
     &                                                    Tindex))
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
# endif
!
!  Read in tangent linear 3D V-momentum component.
!
        IF (get_var(idVvel)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idVvel)),   &
     &                        varid)
          IF (foundit) THEN
            gtype=var_flag(varid)*v3dvar
            status=nf_fread3d(ng, IDmod, ncname, ncINPid,               &
     &                        Vname(1,idVvel), varid,                   &
     &                        InpRec, gtype, Vsize,                     &
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

# ifdef ADJUST_BOUNDARY
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
# endif
!
!  Read in tangent linear tracer type variables.
!
        DO itrc=1,NT(ng)
          IF (get_var(idTvar(itrc))) THEN
            foundit=find_string(var_name, n_var,                        &
     &                          TRIM(Vname(1,idTvar(itrc))), varid)
            IF (foundit) THEN
              gtype=var_flag(varid)*r3dvar
              status=nf_fread3d(ng, IDmod, ncname, ncINPid,             &
     &                          Vname(1,idTvar(itrc)), varid,           &
     &                          InpRec, gtype, Vsize,                   &
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
              IF (FoundError(exit_flag, nf90_noerr,                     &
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
              IF (FoundError(exit_flag, nf90_noerr,                     &
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
     &                          TRIM(Vname(1,idTsur(itrc))), varid)
            IF (foundit) THEN
              gtype=var_flag(varid)*r3dvar
              scale=1.0_dp
              status=nf_fread3d(ng, IDmod, ncname, ncINPid,             &
     &                          Vname(1,idTsur(itrc)), varid,           &
     &                          InpRec, gtype, Vsize,                   &
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
# endif
# ifdef SEDIMENT
!
!  Read in tangent linear sediment fraction of each size class in each
!  bed layer.
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
# ifdef CHECKSUM
     &                          SEDBED(ng) % tl_bed_frac(:,:,:,i),      &
     &                          checksum = Fhash)
# else
     &                          SEDBED(ng) % tl_bed_frac(:,:,:,i))
# endif
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
              IF (FoundError(exit_flag, nf90_noerr,                     &
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
     &                          SEDBED(ng) % tl_bed_mass(:,:,:,         &
     &                                                   Tindex,i),     &
     &                          checksum = Fhash)
#  else
     &                          SEDBED(ng) % tl_bed_mass(:,:,:,         &
     &                                                   Tindex,i))
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
              IF (FoundError(exit_flag, nf90_noerr,                     &
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
     &                          SEDBED(ng) % tl_bed(:,:,:,i),           &
     &                          checksum = Fhash)
#  else
     &                          SEDBED(ng) % tl_bed(:,:,:,i))
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
!  Read in tangent linear sediment fraction of bed load.
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
     &                          SEDBED(ng) % tl_bedldu(:,:,i),          &
     &                          checksum = Fhash)
#   else
     &                          SEDBED(ng) % tl_bedldu(:,:,i))
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
     &                          SEDBED(ng) % tl_bedldv(:,:,i),          &
     &                          checksum = Fhash)
#   else
     &                          SEDBED(ng) % tl_bedldv(:,:,i))
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
!  Read in tangent linear sediment properties in exposed bed layer.
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
     &                          SEDBED(ng) % tl_bottom(:,:,i),          &
     &                          checksum = Fhash)
#  else
     &                          SEDBED(ng) % tl_bottom(:,:,i))
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
#endif
