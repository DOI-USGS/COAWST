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
** It process the reading of ROMS ADM state from input NetCDF files.  **
**                                                                    **
** It is used to read the adjoint kernel state variables. Data is     **       
** loaded into "ad_" arrays.                                          **
**                                                                    **
************************************************************************
*/

!
!  Read in adjoint free-surface.
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
     &                        OCEAN(ng) % ad_zeta(:,:,Tindex),          &
     &                        checksum = Fhash)
#else
     &                        OCEAN(ng) % ad_zeta(:,:,Tindex))
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
!  Read in adjoint free-surface open boundaries adjustments.
!
        IF (get_var(idSbry(isFsur)).and.                                &
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
     &                            BOUNDARY(ng) % ad_zeta_obc(:,:,:,     &
     &                                                       Tindex),   &
     &                            checksum = Fhash)
# else
     &                            BOUNDARY(ng) % ad_zeta_obc(:,:,:,     &
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
!  Read in adjoint 2D U-momentum component.
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
     &                        OCEAN(ng) % ad_ubar(:,:,Tindex),          &
     &                        checksum = Fhash)
#else
     &                        OCEAN(ng) % ad_ubar(:,:,Tindex))
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
!  Read in 2D adjoint U-momentum component open boundaries adjustments.
!
        IF (get_var(idSbry(isUbar)).and.                                &
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
     &                            BOUNDARY(ng) % ad_ubar_obc(:,:,:,     &
     &                                                       Tindex),   &
     &                            checksum = Fhash)
# else
     &                            BOUNDARY(ng) % ad_ubar_obc(:,:,:,     &
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
!  Read in adjoint 2D V-momentum component.
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
     &                        OCEAN(ng) % ad_vbar(:,:,Tindex),          &
     &                        checksum = Fhash)
#else
     &                        OCEAN(ng) % ad_vbar(:,:,Tindex))
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
        IF (get_var(idSbry(isVbar)).and.                                &
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
     &                            BOUNDARY(ng) % ad_vbar_obc(:,:,:,     &
     &                                                       Tindex),   &
     &                            checksum = Fhash)
# else
     &                            BOUNDARY(ng) % ad_vbar_obc(:,:,:,     &
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
!  Read in adjoint linear surface U-momentum stress.
!
        IF (get_var(idUsms)) THEN
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
     &                        FORCES(ng) % ad_ustr(:,:,:,Tindex),       &
     &                        checksum = Fhash)
# else
     &                        FORCES(ng) % ad_ustr(:,:,:,Tindex))
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
     &                            ', adjusted ad_ustr', Fmin, Fmax,     &
     &                            Fhash
# else
                WRITE (stdout,70) TRIM(Vname(2,idUsms))//               &
     &                            ', adjusted ad_ustr', Fmin, Fmax
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
!  Read in adjoint linear surface V-momentum stress.
!
        IF (get_var(idVsms)) THEN
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
     &                        FORCES(ng) % ad_vstr(:,:,:,Tindex),       &
     &                        checksum = Fhash)
# else
     &                        FORCES(ng) % ad_vstr(:,:,:,Tindex))
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
     &                            ', adjusted ad_vstr', Fmin, Fmax,     &
     &                            Fhash
# else
                WRITE (stdout,70) TRIM(Vname(2,idVsms))//               &
     &                            ', adjusted ad_vstr', Fmin, Fmax
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
!  Read in adjoint 3D U-momentum component.
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
     &                        OCEAN(ng) % ad_u(:,:,:,Tindex),           &
     &                        checksum = Fhash)
# else
     &                        OCEAN(ng) % ad_u(:,:,:,Tindex))
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
!  Read in adjoint 3D U-momentum component open boundaries adjustments.
!
        IF (get_var(idSbry(isUvel)).and.                                &
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
     &                            BOUNDARY(ng) % ad_u_obc(:,:,:,:,      &
     &                                                    Tindex),      &
     &                            checksum = Fhash)
#  else
     &                            BOUNDARY(ng) % ad_u_obc(:,:,:,:,      &
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
!  Read in adjoint 3D V-momentum component.
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
     &                        OCEAN(ng) % ad_v(:,:,:,Tindex),           &
     &                        checksum = Fhash)
# else
     &                        OCEAN(ng) % ad_v(:,:,:,Tindex))
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
        IF (get_var(idSbry(isVvel)).and.                                &
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
     &                            BOUNDARY(ng) % ad_v_obc(:,:,:,:,      &
     &                                                    Tindex),      &
     &                            checksum = Fhash)
#  else
     &                            BOUNDARY(ng) % ad_v_obc(:,:,:,:,      &
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
!  Read in adjoint tracer type variables.
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
     &                          OCEAN(ng) % ad_t(:,:,:,Tindex,itrc),    &
     &                          checksum = Fhash)
# else
     &                          OCEAN(ng) % ad_t(:,:,:,Tindex,itrc))
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
# ifdef CHECKSUM
                  WRITE (stdout,70) TRIM(Vname(2,idTvar(itrc))),        &
     &                              Fmin, Fmin, Fhash
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
!  Read in adjoint 3D tracers open boundaries adjustments.
!
        DO itrc=1,NT(ng)
          IF (get_var(idSbry(isTvar(itrc))).and.                        &
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
     &                              BOUNDARY(ng) % ad_t_obc(:,:,:,:,    &
     &                                                   Tindex,itrc),  &
     &                              checksum = Fhash)
#  else
     &                              BOUNDARY(ng) % ad_t_obc(:,:,:,:,    &
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
!  Read in adjoint surface tracers flux.
!
        DO itrc=1,NT(ng)
          IF (get_var(idTsur(itrc)).and.Lstflux(itrc,ng)) THEN
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
#  ifdef CHECKSUM
     &                          FORCES(ng) % ad_tflux(:,:,:,            &
     &                                                Tindex,itrc),     &
     &                          checksum = Fhash)
#  else
     &                          FORCES(ng) % ad_tflux(:,:,:,            &
     &                                                Tindex,itrc))
#  endif
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
#  ifdef CHECKSUM
                  WRITE (stdout,70) TRIM(Vname(2,idTsur(itrc)))//       &
     &                              ', adjusted ad_tflux', Fmin, Fmax,  &
     &                              Fhash
#  else
                  WRITE (stdout,70) TRIM(Vname(2,idTsur(itrc)))//       &
     &                              ', adjusted ad_tflux', Fmin, Fmax
#  endif
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
!  Read in adjoint sediment fraction of each size class in each bed
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
     &                          SEDBED(ng) % ad_bed_frac(:,:,:,i),      &
     &                          checksum = Fhash)
#  else
     &                          SEDBED(ng) % ad_bed_frac(:,:,:,i))
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
!  Read in adjoint sediment mass of each size class in each bed layer.
!
          IF (get_var(idBmas(i))) THEN
            foundit=find_string(var_name, n_var,
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
     &                          SEDBED(ng) % ad_bed_mass(:,:,:,         &
                                                         Tindex,i),     &
     &                          checksum = Fhash)
#  else
     &                          SEDBED(ng) % ad_bed_mass(:,:,:,         &
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
!  Read in adjoint sediment properties in each bed layer.
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
     &                          SEDBED(ng) % ad_bed(:,:,:,i),           &
     &                          checksum = Fhash)
#  else
     &                          SEDBED(ng) % ad_bed(:,:,:,i))
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
!  Read in adjoint sediment fraction of bed load.
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
     &                          SEDBED(ng) % ad_bedldu(:,:,i),          &
     &                          checksum = Fhash)
#   else
     &                          SEDBED(ng) % ad_bedldu(:,:,i))
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
     &                          SEDBED(ng) % ad_bedldv(:,:,i),          &
     &                          checksum = Fhash)
#   else
     &                          SEDBED(ng) % ad_bedldv(:,:,i))
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
!  Read in adjoint sediment properties in exposed bed layer.
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
     &                          SEDBED(ng) % ad_bottom(:,:,i),          &
     &                          checksum = Fhash)
#  else
     &                          SEDBED(ng) % ad_bottom(:,:,i))
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
