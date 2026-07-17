#include "cppdefs.h"
      MODULE esmf_coupler_mod

#if defined MODEL_COUPLING && defined ESMF_LIB
!
!git $Id$
!svn $Id: esmf_coupler.h 1151 2023-02-09 03:08:53Z arango $
!=======================================================================
!  Copyright (c) 2002-2023 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license         Hernan G. Arango     !
!    See License_ROMS.txt                         Ufuk Utku Turuncoglu !
!=======================================================================
!                                                                      !
!  This module includes the coupler routines for the computation,      !
!  execution, and release of the connectors between source and         !
!  destination fields. It uses ESMF/NUOPC regrid operators for         !
!  interpolation with or without extrapolation support.                !
!                                                                      !
!    Coupler_SetServices    Sets coupler shared-object entry points    !
!                           using NUOPC generic methods for the        !
!                           computation, execution, and release of     !
!                           RouteHandle connectors.                    !
!                                                                      !
!    Coupler_ComputeRH      Computes coupler RouteHandle connectors    !
!                           for interpolation between source and       !
!                           destination fields.                        !
!                                                                      !
!    Coupler_ExecuteRH      Uses coupler RouteHandle connectors to     !
!                           perform the regridding interpolation       !
!                           between source and destination fields.     !
!                                                                      !
!    Coupler_ReleaseRH      Releases coupler RouteHandle connectors.   !
!                                                                      !
!    Coupler_AdjustField    Adjust regrid field to ensure global       !
!                           area integral conservation over the        !
!                           matched regions.                           !
!                                                                      !
!    Coupler_AreaIntegral   Computes source or destination global      !
!                           area integral over the matched regions.    !
!                                                                      !
!    Coupler_FieldCreate    Creates and initializes and new field      !
!                           using specified input field attributes     !
!                           and mask.                                  !
!                                                                      !
!    Coupler_FindUnmapped   Modifies mask associated with regrid       !
!                           field to split masked and unmasked grid    !
!                           cells during regridding with extrapolation !
!                           support.                                   !
!                                                                      !
!  ESMF:   Earth System Modeling Framework (Version 7 or higher)       !
!            https://www.earthsystemcog.org/projects/esmf              !
!                                                                      !
!  NUOPC:  National Unified Operational Prediction Capability          !
!            https://www.earthsystemcog.org/projects/nuopc             !
!                                                                      !
!=======================================================================
!
      USE ESMF
      USE NUOPC
      USE NUOPC_Connector, ONLY:                                        &
     &    NUOPC_SetServices     => SetServices,                         &
     &    NUOPC_Label_ComputeRH => label_ComputeRouteHandle,            &
     &    NUOPC_Label_ExecuteRH => label_ExecuteRouteHandle,            &
     &    NUOPC_Label_ReleaseRH => label_ReleaseRouteHandle,            &
     &    NUOPC_ConnectorGet,                                           &
     &    NUOPC_ConnectorSet
!
      USE mod_esmf_esm          ! ESM coupling structures and variables
!
      implicit none
!
      PUBLIC  :: Coupler_SetServices

      PRIVATE :: Coupler_ComputeRH
      PRIVATE :: Coupler_ExecuteRH
      PRIVATE :: Coupler_ReleaseRH
!
      PRIVATE :: Coupler_AdjustField
      PRIVATE :: Coupler_AreaIntegral
      PRIVATE :: Coupler_FieldCreate
      PRIVATE :: Coupler_FindUnmapped
!
      CONTAINS
!
      SUBROUTINE Coupler_SetServices (coupler, rc)
!
!=======================================================================
!                                                                      !
!  Sets the coupler shared-object entry points using NUOPC generic     !
!  methods for the computation, execution, and release of RouteHandle  !
!  connectors between source (srcFields) and destination (dstFields)   !
!  fields.                                                             !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer, intent(out) :: rc
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", Coupler_SetServices"

      character (ESMF_MAXSTR) :: Cname
!
      TYPE (ESMF_CplComp) :: coupler
!
!-----------------------------------------------------------------------
!  Initialize return code flag to success state (no error).
!-----------------------------------------------------------------------
!
      rc=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!  Querry coupler component.
!-----------------------------------------------------------------------
!
      CALL ESMF_CplCompGet (coupler,                                    &
     &                      name=Cname,                                 &
     &                      rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '==> Entering Coupler_SetServices for ' &
     &                          // TRIM(Cname), ', PET', PETrank
        CALL my_flush (trac)
      END IF
!
!-----------------------------------------------------------------------
!  Register generic methods.
!-----------------------------------------------------------------------
!
      CALL NUOPC_CompDerive (coupler,                                   &
     &                       NUOPC_SetServices,                         &
     &                       rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!-----------------------------------------------------------------------
!  Attach specializing methods.
!-----------------------------------------------------------------------
!
!  Set shared-object entry point for the communication RouteHandle used
!  in the data transfer between connectors.
!
      CALL NUOPC_CompSpecialize (coupler,                               &
     &                           specLabel=NUOPC_Label_ComputeRH,       &
     &                           specRoutine=Coupler_ComputeRH,         &
     &                           rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!  Set shared-object entry point to execute connector operations between
!  source (srcFields) and destination (dstFields) data.
!
      CALL NUOPC_CompSpecialize (coupler,                               &
     &                           specLabel=NUOPC_Label_ExecuteRH,       &
     &                           specRoutine=Coupler_ExecuteRH,         &
     &                           rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!  Set share-object entry point for the release of connector operations.
!
      CALL NUOPC_CompSpecialize (coupler,                               &
     &                           specLabel=NUOPC_Label_ReleaseRH,       &
     &                           specRoutine=Coupler_ReleaseRH,         &
     &                           rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '<== Exiting  Coupler_SetServices for ' &
     &                          // TRIM(Cname), ', PET', PETrank
        CALL my_flush (trac)
      END IF
!
      RETURN
      END SUBROUTINE Coupler_SetServices
!
      SUBROUTINE Coupler_ComputeRH (coupler, rc)
!
!=======================================================================
!                                                                      !
!  Sets coupler RouteHandle connectors between source (srcFields) and  !
!  destination (dstFields) fields for ESMF/NUOPC regridding operators  !
!  with or without extrapolation support.                              !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer, intent(out) :: rc
!
      TYPE (ESMF_CplComp) :: coupler
!
!  Local variable declarations.
!
      logical :: rh1Exist, rh2Exist
!
      integer :: i, ic, j, localPET, PETcount, MyComm, NcplSets
      integer :: iSrc, iDst, idSrc, idDst, grSrc, grDst
      integer :: etSrc, etDst, itSrc, itDst
      integer :: srcCount, dstCount, itemCount, srcTerm
!
      integer (i4b) :: srcMaskVal, dstMaskVal
      integer (i4b) :: LandValue(1), SeaValue(1)
!
      integer (i4b), allocatable, dimension(:,:) :: tlw, tuw
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", Coupler_ComputeRH"

      character (ESMF_MAXSTR) :: msgString
      character (ESMF_MAXSTR) :: Cname, Dname, Fname, Rname, Sname
!
      character (ESMF_MAXSTR), pointer :: CplSetList(:) => NULL()
      character (ESMF_MAXSTR), pointer :: dstList(:) => NULL()
      character (ESMF_MAXSTR), pointer :: srcList(:) => NULL()
!
      TYPE (ESMF_ExtrapMethod_Flag)   :: extrapMethod
      TYPE (ESMF_Field)               :: dstField, srcField, tmpField
      TYPE (ESMF_FieldBundle)         :: dstFields, srcFields
      TYPE (ESMF_FieldStatus_Flag)    :: FieldStatus
      TYPE (ESMF_RegridMethod_Flag)   :: regridMethod
      TYPE (ESMF_RouteHandle)         :: routeHandle
      TYPE (ESMF_State)               :: state
      TYPE (ESMF_UnmappedAction_Flag) :: unmap
      TYPE (ESMF_VM)                  :: vm
!
!-----------------------------------------------------------------------
!  Initialize return code flag to success state (no error).
!-----------------------------------------------------------------------
!
      rc=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!  Querry coupler component.
!-----------------------------------------------------------------------
!
!  Querry the coupler for the Virtual Machine (VM) parallel environmemt.
!
      CALL ESMF_CplCompGet (coupler,                                    &
     &                      name=Cname,                                 &
     &                      vm=vm,                                      &
     &                      rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '==> Entering Coupler_ComputeRH for '   &
     &                          // TRIM(Cname), ', PET', PETrank
        CALL my_flush (trac)
      END IF
!
!  Get current parallel node rank and number of nodes.
!
      CALL ESMF_VMGet (vm,                                              &
     &                 localPet=localPET,                               &
     &                 petCount=PETcount,                               &
     &                 mpiCommunicator=MyComm,                          &
     &                 rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!  Set source and destination couple model indices.
!
      DO i=1,Nmodels
        DO j=1,Nmodels
          IF ((CONNECTORS(i,j)%IsActive).and.                           &
     &        (TRIM(CONNECTORS(i,j)%name).eq.TRIM(Cname))) THEN
            iSrc=i
            iDst=j
          END if
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Exchange land-sea mask flag.
!-----------------------------------------------------------------------
!
      LandValue(1)=MODELS(iSrc)%LandValue
      CALL ESMF_VMBroadcast (vm,                                        &
     &                       bcstData=LandValue,                        &
     &                       count=1,                                   &
     &                       rootPet=0,                                 &
     &                       rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
      MODELS(iSrc)%LandValue=LandValue(1)
!
      SeaValue(1)=MODELS(iSrc)%SeaValue
      CALL ESMF_VMBroadcast (vm,                                        &
     &                       bcstData=SeaValue,                         &
     &                       count=1,                                   &
     &                       rootPet=0,                                 &
     &                       rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
      MODELS(iSrc)%SeaValue=SeaValue(1)
!
      LandValue(1)=MODELS(iDst)%LandValue
      CALL ESMF_VMBroadcast (vm,                                        &
     &                       bcstData=LandValue,                        &
     &                       count=1,                                   &
     &                       rootPet=0,                                 &
     &                       rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
      MODELS(iDst)%LandValue=LandValue(1)
!
      SeaValue(1)=MODELS(iDst)%SeaValue
      CALL ESMF_VMBroadcast (vm,                                        &
     &                       bcstData=SeaValue,                         &
     &                       count=1,                                   &
     &                       rootPet=0,                                 &
     &                       rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
      MODELS(iDst)%SeaValue=SeaValue(1)
!
!-----------------------------------------------------------------------
!  Set source and destination masks for connector according to the
!  connector interaction flag.
!-----------------------------------------------------------------------
!
      SELECT CASE (CONNECTORS(iSrc,iDst)%MaskInteraction)
        CASE (OverOcean)
          srcMaskVal=MODELS(iSrc)%LandValue
          dstMaskVal=MODELS(iDst)%LandValue
        CASE (OverLand)
          srcMaskVal=MODELS(iSrc)%SeaValue
          dstMaskVal=MODELS(iDst)%SeaValue
      END SELECT
!
!-----------------------------------------------------------------------
!  Get coupled set list for connector (Cname).
!-----------------------------------------------------------------------
!
      IF ( associated(CplSetList) ) nullify (CplSetList)
      CALL NUOPC_ConnectorGet (coupler,                                 &
     &                         cplSetList=CplSetList,                   &
     &                         rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
      NcplSets=SIZE(CplSetList)
!
!-----------------------------------------------------------------------
!  Inquire size of source and destination field bundles.
!-----------------------------------------------------------------------
!
      CPLSET_LOOP : DO ic=1,NcplSets
!
!  Get source and destination fields for each coupled set.
!
        CALL NUOPC_ConnectorGet (coupler,                               &
     &                           srcFields=srcFields,                   &
     &                           dstFields=dstFields,                   &
     &                           state=state,                           &
     &                           cplSet=CplSetList(ic),                 &
     &                           rc=rc)
        IF (ESMF_LogFoundError(rcToCheck=rc,                            &
     &                         msg=ESMF_LOGERR_PASSTHRU,                &
     &                         line=__LINE__,                           &
     &                         file=MyFile)) THEN
          RETURN
        END IF
!
        CALL ESMF_FieldBundleGet (srcFields,                            &
     &                            fieldCount=srcCount,                  &
     &                            rc=rc)
        IF (ESMF_LogFoundError(rcToCheck=rc,                            &
     &                         msg=ESMF_LOGERR_PASSTHRU,                &
     &                         line=__LINE__,                           &
     &                         file=MyFile)) THEN
          RETURN
        END IF
!
        CALL ESMF_FieldBundleGet (dstFields,                            &
     &                            fieldCount=dstCount,                  &
     &                            rc=rc)
        IF (ESMF_LogFoundError(rcToCheck=rc,                            &
     &                         msg=ESMF_LOGERR_PASSTHRU,                &
     &                         line=__LINE__,                           &
     &                         file=MyFile)) THEN
          RETURN
        END IF
!
        IF ((DebugLevel.gt.0).and.(localPET.eq.0)) THEN
          WRITE (cplout,10) localPET, iSrc, iDst,                       &
     &                      srcMaskVal, dstMaskVal,                     &
     &                      TRIM(CplSetList(ic)),                       &
     &                      TRIM(CONNECTORS(iSrc,iDst)%name)
        END IF
!
!-----------------------------------------------------------------------
!  Get source and destination fields.
!-----------------------------------------------------------------------
!
        DEFINE : IF ((srcCount.eq.dstCount).and.(dstCount.gt. 0)) THEN
!
!  Allocate.
!
          allocate ( srcList(srcCount) )
          allocate ( dstList(dstCount) )
!
!  Get source and destination fields.
!
          CALL ESMF_FieldBundleGet (srcFields,                          &
     &                              fieldNameList=srcList,              &
     &                              rc=rc)
          IF (ESMF_LogFoundError(rcToCheck=rc,                          &
     &                           msg=ESMF_LOGERR_PASSTHRU,              &
     &                           line=__LINE__,                         &
     &                           file=MyFile)) THEN
            RETURN
          END IF
!
          CALL ESMF_FieldBundleGet (dstFields,                          &
     &                              fieldNameList=dstList,              &
     &                              rc=rc)
          IF (ESMF_LogFoundError(rcToCheck=rc,                          &
     &                           msg=ESMF_LOGERR_PASSTHRU,              &
     &                           line=__LINE__,                         &
     &                           file=MyFile)) THEN
            RETURN
          END IF
!
!=======================================================================
!  Create connector RootHandlers for two step interpolation.
!=======================================================================
!
          CREATE : DO i=1,srcCount
!
!  Get source and destination field index.
!
            idSrc=field_index(MODELS(iSrc)%ExportField, srcList(i))
            idDst=field_index(MODELS(iDst)%ImportField, dstList(i))
!
!  Get field name. Both source and destination should have the same
!  short name.
!
            Fname=TRIM(MODELS(iSrc)%ExportField(idSrc)%short_name)
!
!  Get regrid method for interpolation.
!
            itSrc=MODELS(iSrc)%ExportField(idSrc)%itype
            itDst=MODELS(iDst)%ImportField(idDst)%itype
!
            IF (itSrc.NE.itDst) THEN
              WRITE (msgString,'(a)') TRIM(Cname)//                     &
     &          ': SRC and DST field interpolation type does not match!'
              CALL ESMF_LogWrite (TRIM(msgString), ESMF_LOGMSG_ERROR)
              RETURN
            END IF
!
!  Get extrapolation method for unmapped destination points.
!
            etSrc=MODELS(iSrc)%ExportField(idSrc)%etype
            etDst=MODELS(iDst)%ImportField(idDst)%etype
!
            IF (etSrc.NE.etDst) THEN
              WRITE (msgString,'(a)') TRIM(Cname)//                     &
     &          ': SRC and DST field extrapolation type does not match!'
              CALL ESMF_LogWrite (TRIM(msgString), ESMF_LOGMSG_ERROR)
              RETURN
            END IF
!
!  Get grid type.
!
            grSrc=MODELS(iSrc)%ExportField(idSrc)%gtype
            grDst=MODELS(iDst)%ImportField(idDst)%gtype
!
!  Get source field object from bundle.
!
            CALL ESMF_FieldBundleGet (srcFields,                        &
     &                                TRIM(Fname),                      &
     &                                field=srcField,                   &
     &                                rc=rc)
            IF (ESMF_LogFoundError(rcToCheck=rc,                        &
     &                             msg=ESMF_LOGERR_PASSTHRU,            &
     &                             line=__LINE__,                       &
     &                             file=MyFile)) THEN
              RETURN
            END IF
!
            CALL ESMF_FieldGet (srcField,                               &
     &                          status=FieldStatus,                     &
     &                          rc=rc)
            IF (ESMF_LogFoundError(rcToCheck=rc,                        &
     &                             msg=ESMF_LOGERR_PASSTHRU,            &
     &                             line=__LINE__,                       &
     &                             file=MyFile)) THEN
              RETURN
            END IF
!
            IF (FieldStatus.ne.ESMF_FIELDSTATUS_COMPLETE) THEN
              rc=ESMF_RC_OBJ_BAD
              IF (localPET.eq.0) THEN
                IF (FieldStatus.eq.ESMF_FIELDSTATUS_EMPTY) THEN
                  msgString='ESMF_FIELDSTATUS_EMPTY'
                ELSE IF (FieldStatus.eq.ESMF_FIELDSTATUS_GRIDSET) THEN
                  msgString='ESMF_FIELDSTATUS_GRIDSET'
                END IF
                WRITE (cplout,20) 'Source Field: ', TRIM(Fname),        &
     &                            TRIM(msgString)
              END IF
              IF (ESMF_LogFoundError(rcToCheck=rc,                      &
     &                               msg=ESMF_LOGERR_PASSTHRU,          &
     &                               line=__LINE__,                     &
     &                               file=MyFile)) THEN
                RETURN
              END IF
            END IF
!
            IF (DebugLevel.gt.1) THEN
              CALL ESMF_FieldPrint (srcField, rc=rc)
              IF (ESMF_LogFoundError(rcToCheck=rc,                      &
     &                               msg=ESMF_LOGERR_PASSTHRU,          &
     &                               line=__LINE__,                     &
     &                               file=MyFile)) THEN
                RETURN
              END IF
            END IF
!
!  Get destination field object from bundle.
!
            CALL ESMF_FieldBundleGet (dstFields,                        &
     &                                TRIM(Fname),                      &
     &                                field=dstField,                   &
     &                                rc=rc)
            IF (ESMF_LogFoundError(rcToCheck=rc,                        &
     &                             msg=ESMF_LOGERR_PASSTHRU,            &
     &                             line=__LINE__,                       &
     &                             file=MyFile)) THEN
              RETURN
            END IF
!
            CALL ESMF_FieldGet (dstField,                               &
     &                          status=FieldStatus,                     &
     &                          rc=rc)
            IF (ESMF_LogFoundError(rcToCheck=rc,                        &
     &                             msg=ESMF_LOGERR_PASSTHRU,            &
     &                             line=__LINE__,                       &
     &                             file=MyFile)) THEN
              RETURN
            END IF
!
            IF (FieldStatus.ne.ESMF_FIELDSTATUS_COMPLETE) THEN
              rc=ESMF_RC_OBJ_BAD
              IF (localPET.eq.0) THEN
                IF (FieldStatus.eq.ESMF_FIELDSTATUS_EMPTY) THEN
                  msgString='ESMF_FIELDSTATUS_EMPTY'
                ELSE IF (FieldStatus.eq.ESMF_FIELDSTATUS_GRIDSET) THEN
                  msgString='ESMF_FIELDSTATUS_GRIDSET'
                END IF
                WRITE (cplout,20) 'Destination Field: ', TRIM(Fname),   &
     &                            TRIM(msgString)
              END IF
              IF (ESMF_LogFoundError(rcToCheck=rc,                      &
     &                               msg=ESMF_LOGERR_PASSTHRU,          &
     &                               line=__LINE__,                     &
     &                               file=MyFile)) THEN
                RETURN
              END IF
            END IF
!
            IF (DebugLevel.gt.1) THEN
              CALL ESMF_FieldPrint (dstField, rc=rc)
              IF (ESMF_LogFoundError(rcToCheck=rc,                      &
     &                               msg=ESMF_LOGERR_PASSTHRU,          &
     &                               line=__LINE__,                     &
     &                               file=MyFile)) THEN
                RETURN
              END IF
            END IF
!
!-----------------------------------------------------------------------
!  REGRIDDING: Set RouteHandle for two-steps extrapolation.
!-----------------------------------------------------------------------
!
            QUERRY : IF (etSrc.eq.E2steps) THEN
!
!  Check 1st RouteHandle. If the source field is the DATA component,
!  there is a RouteHandle for each export field. Examples:
!
!    rh_Center_Corner_BLIN_WRF_01-TO-ROMS_01               source WRF_01
!    rh_SST_Center_Corner_BLIN_DATA-TO-ROMS_01             source DATA
!
              IF (iSrc.eq.Idata) THEN
                Rname='rh_'//TRIM(srcList(i))//'_'//                    &
     &                TRIM(GridType (grSrc ))//'_'//                    &
     &                TRIM(GridType (grDst ))//'_'//                    &
     &                TRIM(IntrpType(Ibilin))//'_'//                    &
     &                TRIM(ExtrpType(etSrc ))//'_'//                    &
     &                TRIM(CplSetList(ic))//'_'//                       &
     &                TRIM(Cname)
              ELSE
                Rname='rh_'//                                           &
     &                TRIM(GridType (grSrc ))//'_'//                    &
     &                TRIM(GridType (grDst ))//'_'//                    &
     &                TRIM(IntrpType(Ibilin))//'_'//                    &
     &                TRIM(ExtrpType(etSrc ))//'_'//                    &
     &                TRIM(CplSetList(ic))//'_'//                       &
     &                TRIM(Cname)
              END IF
!
              CALL ESMF_StateGet (state,                                &
     &                            itemSearch=TRIM(Rname),               &
     &                            itemCount=itemCount,                  &
     &                            rc=rc)
              IF (ESMF_LogFoundError(rcToCheck=rc,                      &
     &                               msg=ESMF_LOGERR_PASSTHRU,          &
     &                               line=__LINE__,                     &
     &                               file=MyFile)) THEN
                RETURN
              END IF
!
              IF (itemCount.le.0) THEN
                rh1Exist=.FALSE.
              ELSE
                rh1Exist=.TRUE.
              END IF
              rh2Exist=.FALSE.
!
!  Debug: report exchanged fields before regridding.
!
              IF ((DebugLevel.gt.0).and.(localPET.eq.0)) THEN
                WRITE (cplout,30) TRIM(CplSetList(ic)), TRIM(Cname),    &
     &          TRIM(MODELS(iSrc)%ExportField(idSrc)%short_name),       &
     &          TRIM(GridType (MODELS(iSrc)%ExportField(idSrc)%gtype)), &
     &          TRIM(MODELS(iDst)%ImportField(idDst)%short_name),       &
     &          TRIM(GridType (MODELS(iDst)%ImportField(idDst)%gtype)), &
     &          TRIM(IntrpType(MODELS(iSrc)%ExportField(idSrc)%itype)), &
     &          rh1Exist, rh2Exist
                CALL my_flush (cplout)
              END IF
!
!  Create 1st RouteHandle.
!
              IF (.not.rh1Exist) THEN
                unmap=ESMF_UNMAPPEDACTION_IGNORE
                regridMethod=ESMF_REGRIDMETHOD_BILINEAR
!
                srcTerm=0
!
                CALL ESMF_FieldRegridStore (srcField=srcField,          &
     &                                      dstField=dstField,          &
     &                             srcMaskValues=(/srcMaskVal/),        &
     &                             dstMaskValues=(/dstMaskVal/),        &
     &                             unmappedaction=unmap,                &
     &                             routeHandle=routeHandle,             &
     &                             regridmethod=regridMethod,           &
     &                             ignoreDegenerate=.TRUE.,             &
     &                             srcTermProcessing=srcTerm,           &
     &                                      rc=rc)
                IF (ESMF_LogFoundError(rcToCheck=rc,                    &
     &                                 msg=ESMF_LOGERR_PASSTHRU,        &
     &                                 line=__LINE__,                   &
     &                                 file=MyFile)) THEN
                  RETURN
                END IF
!
!  Set name to 1st RouteHandle.
!
                CALL ESMF_RouteHandleSet (routeHandle,                  &
     &                                    name=TRIM(Rname),             &
     &                                    rc=rc)
                IF (ESMF_LogFoundError(rcToCheck=rc,                    &
     &                                 msg=ESMF_LOGERR_PASSTHRU,        &
     &                                 line=__LINE__,                   &
     &                                 file=MyFile)) THEN
                  RETURN
                END IF
!
!  Add 1st RouteHandle to the state.
!
                CALL ESMF_StateAdd (state,                              &
     &                              (/ routeHandle /),                  &
     &                              rc=rc)
                IF (ESMF_LogFoundError(rcToCheck=rc,                    &
     &                                 msg=ESMF_LOGERR_PASSTHRU,        &
     &                                 line=__LINE__,                   &
     &                                 file=MyFile)) THEN
                  RETURN
                END IF
!
!  Debug: report successful computation of regridding 1st RouteHandle.
!
                IF ((DebugLevel.gt.0).and.(localPET.eq.0)) THEN
                  WRITE (cplout,40) TRIM(Rname)
                  CALL my_flush (cplout)
                END IF
              END IF
!
!  Check 2nd routehandle (i.e. rh_Center_Corner_NS2D_ATM-OCN).
!
!    rh_Center_Corner_NS2D_WRF_01-TO-ROMS_01     source WRF_01 component
!    rh_SST_Center_Corner_NS2D_DATA-TO-ROMS_01   source DATA component
!
              IF (iSrc.eq.Idata) THEN
                Rname='rh_'//TRIM(srcList(i))//'_'//                    &
     &                TRIM(GridType (grSrc ))//'_'//                    &
     &                TRIM(GridType (grDst ))//'_'//                    &
     &                TRIM(IntrpType(InStoD))//'_'//                    &
     &                TRIM(ExtrpType(etSrc ))//'_'//                    &
     &                TRIM(CplSetList(ic))//'_'//                       &
     &                TRIM(Cname)
              ELSE
                Rname='rh_'//                                           &
     &                TRIM(GridType (grSrc ))//'_'//                    &
     &                TRIM(GridType (grDst ))//'_'//                    &
     &                TRIM(IntrpType(InStoD))//'_'//                    &
     &                TRIM(ExtrpType(etSrc ))//'_'//                    &
     &                TRIM(CplSetList(ic))//'_'//                       &
     &                TRIM(Cname)
              END IF
!
              CALL ESMF_StateGet (state,                                &
     &                            itemSearch=TRIM(Rname),               &
     &                            itemCount=itemCount,                  &
     &                            rc=rc)
              IF (ESMF_LogFoundError(rcToCheck=rc,                      &
     &                               msg=ESMF_LOGERR_PASSTHRU,          &
     &                               line=__LINE__,                     &
     &                               file=MyFile)) THEN
                RETURN
              END IF
!
              IF (itemCount.le.0) THEN
                rh2Exist=.FALSE.
              ELSE
                rh2Exist=.TRUE.
              END IF
!
!  Create temporary field in destination grid.
!
              IF (.not.rh2Exist) THEN
                tmpField=Coupler_FieldCreate(dstField, 'temp_field',    &
     &                                       1.0_dp, -1_i4b, rc)
!
!  Modify grid mask to split masked and unmasked grid cells.
!
                CALL Coupler_FindUnmapped (srcField, dstField,          &
     &                                     srcMaskVal, dstMaskVal,      &
     &                                     iSrc, iDst, rc)
                IF (ESMF_LogFoundError(rcToCheck=rc,                    &
     &                                 msg=ESMF_LOGERR_PASSTHRU,        &
     &                                 line=__LINE__,                   &
     &                                 file=MyFile)) THEN
                  RETURN
                END IF
!
!  Create 2nd RouteHandle.
!
                unmap=ESMF_UNMAPPEDACTION_IGNORE
                regridMethod=ESMF_REGRIDMETHOD_NEAREST_STOD
!
                srcTerm=0
!
                CALL ESMF_FieldRegridStore (srcField=tmpField,          &
     &                                      dstField=dstField,          &
     &                             srcMaskValues=(/dstMaskVal,          &
     &                                             UNMAPPED_MASK/),     &
     &                             dstMaskValues=(/dstMaskVal,          &
     &                                             MAPPED_MASK/),       &
     &                             unmappedaction=unmap,                &
     &                             routeHandle=routeHandle,             &
     &                             regridmethod=regridMethod,           &
     &                             srcTermProcessing=srcTerm,           &
     &                             ignoreDegenerate=.TRUE.,             &
     &                                      rc=rc)
                IF (ESMF_LogFoundError(rcToCheck=rc,                    &
     &                                 msg=ESMF_LOGERR_PASSTHRU,        &
     &                                 line=__LINE__,                   &
     &                                 file=MyFile)) THEN
                  RETURN
                END IF
!
!  Add name to 2nd RouteHandle.
!
                CALL ESMF_RouteHandleSet (routeHandle,                  &
     &                                    name=TRIM(Rname),             &
     &                                    rc=rc)
                IF (ESMF_LogFoundError(rcToCheck=rc,                    &
     &                                 msg=ESMF_LOGERR_PASSTHRU,        &
     &                                 line=__LINE__,                   &
     &                                 file=MyFile)) THEN
                  RETURN
                END IF
!
!  Add 2nd RouteHandle to the state.
!
                CALL ESMF_StateAdd (state,                              &
     &                              (/ routeHandle /),                  &
     &                              rc=rc)
                IF (ESMF_LogFoundError(rcToCheck=rc,                    &
     &                                 msg=ESMF_LOGERR_PASSTHRU,        &
     &                                 line=__LINE__,                   &
     &                                 file=MyFile)) THEN
                  RETURN
                END IF
!
!  Delete temporary field.
!
                CALL ESMF_FieldDestroy (tmpField, rc=rc)
                IF (ESMF_LogFoundError(rcToCheck=rc,                    &
     &                                 msg=ESMF_LOGERR_PASSTHRU,        &
     &                                 line=__LINE__,                   &
     &                                 file=MyFile)) THEN
                  RETURN
                END IF
!
!  Debug: report successful computation of regridding 2nd RouteHandle.
!
                IF ((DebugLevel.gt.0).and.(localPET.eq.0)) THEN
                  WRITE (cplout,40) TRIM(Rname)
                  CALL my_flush (cplout)
                END IF
              END IF
!
!-----------------------------------------------------------------------
!  REGRIDDING: create RouteHandle for one-step interpolation
!              (native extrapolation of unmapped points is possible)
!-----------------------------------------------------------------------
!
            ELSE
!
              IF (iSrc.eq.Idata) THEN
                Rname='rh_'//TRIM(srcList(i))//'_'//                    &
     &                TRIM(GridType  (grSrc))//'_'//                    &
     &                TRIM(GridType  (grDst))//'_'//                    &
     &                TRIM(IntrpType (itSrc))//'_'//                    &
     &                TRIM(ExtrpType (etSrc))//'_'//                    &
     &                TRIM(CplSetList(ic))//'_'//                       &
     &                TRIM(Cname)
              ELSE
                Rname='rh_'//                                           &
     &                TRIM(GridType  (grSrc))//'_'//                    &
     &                TRIM(GridType  (grDst))//'_'//                    &
     &                TRIM(IntrpType (itSrc))//'_'//                    &
     &                TRIM(ExtrpType (etSrc))//'_'//                    &
     &                TRIM(CplSetList(ic))//'_'//                       &
     &                TRIM(Cname)
              END IF
!
              CALL ESMF_StateGet (state,                                &
     &                            itemSearch=TRIM(Rname),               &
     &                            itemCount=itemCount,                  &
     &                            rc=rc)
              IF (ESMF_LogFoundError(rcToCheck=rc,                      &
     &                               msg=ESMF_LOGERR_PASSTHRU,          &
     &                               line=__LINE__,                     &
     &                               file=MyFile)) THEN
                RETURN
              END IF
!
              IF (itemCount.le.0) THEN
                rh1Exist=.FALSE.
              ELSE
                rh1Exist=.TRUE.
              END IF
              rh2Exist=.FALSE.
!
!  Debug: report exchanged fields before regridding.
!
              IF ((DebugLevel.gt.0).and.(localPET.eq.0)) THEN
                WRITE (cplout,30) TRIM(CplSetList(ic)), TRIM(Cname),    &
     &          TRIM(MODELS(iSrc)%ExportField(idSrc)%short_name),       &
     &          TRIM(GridType (MODELS(iSrc)%ExportField(idSrc)%gtype)), &
     &          TRIM(MODELS(iDst)%ImportField(idDst)%short_name),       &
     &          TRIM(GridType (MODELS(iDst)%ImportField(idDst)%gtype)), &
     &          TRIM(IntrpType(MODELS(iSrc)%ExportField(idSrc)%itype)), &
     &          rh1Exist, rh2Exist
                CALL my_flush (cplout)
              END IF
!
!  Regrid method from source to destination.
!
              IF (.not.rh1Exist) THEN
                unmap=ESMF_UNMAPPEDACTION_IGNORE
                IF (itSrc.eq.Ibilin) THEN
                  regridMethod=ESMF_REGRIDMETHOD_BILINEAR
                ELSE IF (itSrc.eq.Ipatch) THEN
                  regridMethod=ESMF_REGRIDMETHOD_PATCH
                ELSE IF (itSrc.eq.Iconsv1) THEN
                  regridMethod=ESMF_REGRIDMETHOD_CONSERVE
                ELSE IF (itSrc.eq.InStoD) THEN
                  regridMethod=ESMF_REGRIDMETHOD_NEAREST_STOD
                ELSE IF (itSrc.eq.InDtoS) THEN
                  regridMethod=ESMF_REGRIDMETHOD_NEAREST_DTOS
                ELSE
                  WRITE (msgString,'(a)') TRIM(Cname)//': selected '//  &
     &                       'interpolation type is not supported! '//  &
     &                                    IntrpType(itSrc)
                  CALL ESMF_LogWrite (TRIM(msgString),ESMF_LOGMSG_ERROR)
                  CALL ESMF_Finalize (endflag=ESMF_END_ABORT)
                END IF
!
                IF (etSrc.eq.Enone) THEN
                  extrapMethod=ESMF_EXTRAPMETHOD_NONE
                ELSE IF (etSrc.eq.ExStoD) THEN
                  extrapMethod=ESMF_EXTRAPMETHOD_NEAREST_STOD
                ELSE IF (etSrc.eq.Eidavg) THEN
                  extrapMethod=ESMF_EXTRAPMETHOD_NEAREST_IDAVG
                ELSE IF (etSrc.eq.Ecreep) THEN
                  extrapMethod=ESMF_EXTRAPMETHOD_CREEP
                ELSE
                  WRITE (msgString,'(a)') TRIM(Cname)//': selected '//  &
     &                       'extrapolation type is not supported! '//  &
     &                                    ExtrpType(etSrc)
                  CALL ESMF_LogWrite (TRIM(msgString),ESMF_LOGMSG_ERROR)
                  CALL ESMF_Finalize (endflag=ESMF_END_ABORT)
                END IF
!
                srcTerm=0
!
                SELECT CASE (Cmodel(iSrc))
                  CASE ('ROMS')
                    LandValue(1)=MODELS(iSrc)%LandValue
                    CALL ESMF_FieldRegridStore (srcField=srcField,      &
     &                                          dstField=dstField,      &
     &                                 srcMaskValues=LandValue,         &
     &                                 unmappedaction=unmap,            &
     &                                 routeHandle=routeHandle,         &
     &                                 regridmethod=regridMethod,       &
     &                                 extrapMethod=extrapMethod,       &
     &                                 extrapNumLevels=extrapNumLevels, &
     &                                 srcTermProcessing=srcTerm,       &
     &                                 ignoreDegenerate=.TRUE.,         &
     &                                          rc=rc)
                    IF (ESMF_LogFoundError(rcToCheck=rc,                &
     &                                     msg=ESMF_LOGERR_PASSTHRU,    &
     &                                     line=__LINE__,               &
     &                                     file=MyFile)) THEN
                      RETURN
                    END IF
                  CASE ('DATA')
                    LandValue(1)=MODELS(iSrc)%LandValue
                    CALL ESMF_FieldRegridStore (srcField=srcField,      &
     &                                          dstField=dstField,      &
     &                                 srcMaskValues=LandValue,         &
     &                                 unmappedaction=unmap,            &
     &                                 routeHandle=routeHandle,         &
     &                                 regridmethod=regridMethod,       &
     &                                 extrapMethod=extrapMethod,       &
     &                                 extrapNumLevels=extrapNumLevels, &
     &                                 srcTermProcessing=srcTerm,       &
     &                                 ignoreDegenerate=.TRUE.,         &
     &                                          rc=rc)
                    IF (ESMF_LogFoundError(rcToCheck=rc,                &
     &                                     msg=ESMF_LOGERR_PASSTHRU,    &
     &                                     line=__LINE__,               &
     &                                     file=MyFile)) THEN
                      RETURN
                    END IF
                  CASE DEFAULT
                    CALL ESMF_FieldRegridStore (srcField=srcField,      &
     &                                          dstField=dstField,      &
     &                                 unmappedaction=unmap,            &
     &                                 routeHandle=routeHandle,         &
     &                                 regridmethod=regridMethod,       &
     &                                 extrapMethod=extrapMethod,       &
     &                                 extrapNumLevels=extrapNumLevels, &
     &                                 srcTermProcessing=srcTerm,       &
     &                                 ignoreDegenerate=.TRUE.,         &
     &                                          rc=rc)
                    IF (ESMF_LogFoundError(rcToCheck=rc,                &
     &                                     msg=ESMF_LOGERR_PASSTHRU,    &
     &                                     line=__LINE__,               &
     &                                     file=MyFile)) THEN
                      RETURN
                    END IF
                END SELECT
!
!  Add name to RouteHandle.
!
                CALL ESMF_RouteHandleSet (routeHandle,                  &
     &                                    name=TRIM(Rname),             &
     &                                    rc=rc)
                IF (ESMF_LogFoundError(rcToCheck=rc,                    &
     &                                 msg=ESMF_LOGERR_PASSTHRU,        &
     &                                 line=__LINE__,                   &
     &                                 file=MyFile)) THEN
                  RETURN
                END IF
!
!  Add RouteHandle to the state.
!
                CALL ESMF_StateAdd (state,                              &
     &                              (/ routeHandle /),                  &
     &                              rc=rc)
                IF (ESMF_LogFoundError(rcToCheck=rc,                    &
     &                                 msg=ESMF_LOGERR_PASSTHRU,        &
     &                                 line=__LINE__,                   &
     &                                 file=MyFile)) THEN
                  RETURN
                END IF
                rh1Exist=.TRUE.
!
!  Debug: report successful computation of regridding 1st RouteHandle.
!
                IF ((DebugLevel.gt.0).and.(localPET.eq.0)) THEN
                  WRITE (cplout,40) TRIM(Rname)
                END IF
              END IF

            END IF QUERRY
          END DO CREATE
        END IF DEFINE
      END DO CPLSET_LOOP
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '<== Exiting  Coupler_ComputeRH for '   &
     &                          // TRIM(Cname), ', PET', PETrank
        CALL my_flush (trac)
      END IF
      IF (DebugLevel.gt.0) CALL my_flush (cplout)
!
 10   FORMAT (4x,'RouteHandle - PET = ',i0,' iSrc = ',i0,' iDst = ',i0, &
     &        ' srcMask = ',i0,' dstMask = ',i0,', cplSet = ',a,', ',a)
 20   FORMAT (' Coupler_ComputerRH - ',a,                               &
     &        '''',a,''' has an incorrect status',/,22x,a)
 30   FORMAT (4x,'RouteHandle - ESMF: ',a,', ',a,', ',a,                &
     &        ' [',a,'] to ',a,' [',a,']',' >> ',a, ' - ',l1,' - ',l1)
 40   FORMAT (18x,'Computed interpolant ''',a,''',  sucessfully')
!
      RETURN
      END SUBROUTINE Coupler_ComputeRH
!
      SUBROUTINE Coupler_ExecuteRH (coupler, rc)
!
!=======================================================================
!                                                                      !
!  Performs the interpolation between source and destination fields.   !
!  It uses REGRID with or without extrapolation support.               !
!                                                                      !
!=======================================================================
!
      USE strings_mod, ONLY : lowercase
!
!  Imported variable declarations.
!
      integer, intent(out) :: rc
!
      TYPE (ESMF_CplComp) :: coupler
!
!  Local variable declarations.
!
      logical :: IsValid
!
      integer :: localPET, PETcount, MyComm
      integer :: i, ic, is, j, srcCount, dstCount, NcplSets
      integer :: iSrc, iDst, idSrc, idDst, grSrc, grDst
      integer :: etSrc, etDst, itSrc, itDst
!
      real (dp) :: src_total, dst_total, rel_error
!
      real (dp), dimension(:,:), pointer :: ptr2d => NULL()
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", Coupler_ExecuteRH"

      character (len=19     ) :: dstFileString, srcFileString
      character (ESMF_MAXSTR) :: msgString, Cname, Fname, Rname
      character (ESMF_MAXSTR) :: dstTimeString, srcTimeString
      character (ESMF_MAXSTR) :: dstFile, srcFile
!
      character (ESMF_MAXSTR), pointer :: CplSetList(:) => NULL()
      character (ESMF_MAXSTR), pointer :: dstList(:) => NULL()
      character (ESMF_MAXSTR), pointer :: srcList(:) => NULL()
!
      TYPE (ESMF_Field)       :: srcField, dstField, tmpField
      TYPE (ESMF_FieldBundle) :: dstFields, srcFields
      TYPE (ESMF_RouteHandle) :: routeHandle
      TYPE (ESMF_State)       :: state
      TYPE (ESMF_Time)        :: dstTime, srcTime
      TYPE (ESMF_VM)          :: vm
!
!-----------------------------------------------------------------------
!  Initialize return code flag to success state (no error).
!-----------------------------------------------------------------------
!
      rc=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!  Query coupler component.
!-----------------------------------------------------------------------
!
!  Querry the coupler for the Virtual Machine (VM) parallel environmemt.
!
      CALL ESMF_CplCompGet (coupler,                                    &
     &                      name=Cname,                                 &
     &                      vm=vm,                                      &
     &                      rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '==> Entering Coupler_ExecuteRH for '   &
     &                          // TRIM(Cname), ', PET', PETrank
        CALL my_flush (trac)
      END IF
!
!  Get current parallel node rank and number of nodes.
!
      CALL ESMF_VMGet (vm,                                              &
     &                 localPet=localPET,                               &
     &                 petCount=PETcount,                               &
     &                 mpiCommunicator=MyComm,                          &
     &                 rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!  Set source and destination couple model indices.
!
      DO i=1,Nmodels
        DO j=1,Nmodels
          IF ((CONNECTORS(i,j)%IsActive).and.                           &
     &        (TRIM(CONNECTORS(i,j)%name).eq.TRIM(Cname))) THEN
            iSrc=i
            iDst=j
          END IF
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Get coupled set list for connector (Cname).
!-----------------------------------------------------------------------
!
      IF ( associated(CplSetList) ) nullify (CplSetList)
      CALL NUOPC_ConnectorGet (coupler,                                 &
     &                         cplSetList=CplSetList,                   &
     &                         rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
      NcplSets=SIZE(CplSetList)
!
!-----------------------------------------------------------------------
!  Inquire about source and destination fields.
!-----------------------------------------------------------------------
!
      CPLSET_LOOP : DO ic=1,NcplSets
!
!  Get source and destination fields for each coupled set.
!
        CALL NUOPC_ConnectorGet (coupler,                               &
     &                           srcFields=srcFields,                   &
     &                           dstFields=dstFields,                   &
     &                           state=state,                           &
     &                           cplSet=CplSetList(ic),                 &
     &                           rc=rc)
        IF (ESMF_LogFoundError(rcToCheck=rc,                            &
     &                         msg=ESMF_LOGERR_PASSTHRU,                &
     &                         line=__LINE__,                           &
     &                         file=MyFile)) THEN
          RETURN
        END IF
!
!  Number of source fields.
!
        CALL ESMF_FieldBundleGet (srcFields,                            &
     &                            fieldCount=srcCount,                  &
     &                            rc=rc)
        IF (ESMF_LogFoundError(rcToCheck=rc,                            &
     &                         msg=ESMF_LOGERR_PASSTHRU,                &
     &                         line=__LINE__,                           &
     &                         file=MyFile)) THEN
          RETURN
        END IF
!
!  Number of destination fields.
!
        CALL ESMF_FieldBundleGet (dstFields,                            &
     &                            fieldCount=dstCount,                  &
     &                            rc=rc)
        IF (ESMF_LogFoundError(rcToCheck=rc,                            &
     &                         msg=ESMF_LOGERR_PASSTHRU,                &
     &                         line=__LINE__,                           &
     &                         file=MyFile)) THEN
          RETURN
        END IF
!
!  Source fields names.
!
        allocate ( srcList(srcCount) )

        CALL ESMF_FieldBundleGet (srcFields,                            &
     &                            fieldNameList=srcList,                &
     &                            rc=rc)
        IF (ESMF_LogFoundError(rcToCheck=rc,                            &
     &                         msg=ESMF_LOGERR_PASSTHRU,                &
     &                         line=__LINE__,                           &
     &                         file=MyFile)) THEN
          RETURN
        END IF
!
!  Destination fields names.
!
        allocate ( dstList(dstCount) )

        CALL ESMF_FieldBundleGet (dstFields,                            &
     &                            fieldNameList=dstList,                &
     &                            rc=rc)
        IF (ESMF_LogFoundError(rcToCheck=rc,                            &
     &                         msg=ESMF_LOGERR_PASSTHRU,                &
     &                         line=__LINE__,                           &
     &                         file=MyFile)) THEN
          RETURN
        END IF
!
!=======================================================================
!  Interpolate or extrapolate between source and destination fields.
!=======================================================================
!
        EXCHANGE : DO i=1,srcCount
!
!  Set source and destination field index.
!
          idSrc=field_index(MODELS(iSrc)%ExportField, srcList(i))
          idDst=field_index(MODELS(iDst)%ImportField, dstList(i))
!
!  Get field name. Both source and destination should have the same
!  short name.
!
          Fname=TRIM(MODELS(iSrc)%ExportField(idSrc)%short_name)
!
!  Set interpolation type.
!
          itSrc=MODELS(iSrc)%ExportField(idSrc)%itype
          itDst=MODELS(iDst)%ImportField(idDst)%itype
!
          IF (itSrc.ne.itDst) THEN
            WRITE (msgString,'(a)') TRIM(Cname)//                       &
     &          ': SRC and DST field interpolation type does not match!'
            CALL ESMF_LogWrite (TRIM(msgString), ESMF_LOGMSG_ERROR)
            RETURN
          END IF
!
!  Get extrapolation method for unmapped destination points.
!
          etSrc=MODELS(iSrc)%ExportField(idSrc)%etype
          etDst=MODELS(iDst)%ImportField(idDst)%etype
!
          IF (etSrc.NE.etDst) THEN
            WRITE (msgString,'(a)') TRIM(Cname)//                       &
     &          ': SRC and DST field extrapolation type does not match!'
            CALL ESMF_LogWrite (TRIM(msgString), ESMF_LOGMSG_ERROR)
            RETURN
          END IF
!
!  Set grid type.
!
          grSrc=MODELS(iSrc)%ExportField(idSrc)%gtype
          grDst=MODELS(iDst)%ImportField(idDst)%gtype
!
!  Get source field object from bundle.
!
          CALL ESMF_FieldBundleGet (srcFields,                          &
     &                              TRIM(Fname),                        &
     &                              field=srcField,                     &
     &                              rc=rc)
          IF (ESMF_LogFoundError(rcToCheck=rc,                          &
     &                           msg=ESMF_LOGERR_PASSTHRU,              &
     &                           line=__LINE__,                         &
     &                           file=MyFile)) THEN
            RETURN
          END IF
!
!  Get destination field object from bundle.
!
          CALL ESMF_FieldBundleGet (dstFields,                          &
     &                              TRIM(Fname),                        &
     &                              field=dstField,                     &
     &                              rc=rc)
          IF (ESMF_LogFoundError(rcToCheck=rc,                          &
     &                           msg=ESMF_LOGERR_PASSTHRU,              &
     &                           line=__LINE__,                         &
     &                           file=MyFile)) THEN
            RETURN
          END IF
!
!-----------------------------------------------------------------------
!  Perform REGRID with two-steps extrapolation.
!-----------------------------------------------------------------------
!
          QUERRY : IF (etSrc.eq.E2steps) THEN
!
!  Check 1st RouteHandle. If the source field is the DATA component,
!  there is a RouteHandle for each export field. Examples:
!
!    rh_Center_Corner_BLIN_ATM-TO-ROMS_01        (source ATM  component)
!    rh_SST_Center_Corner_BLIN_DATA-TO-ROMS_01   (source DATA component)
!
            IF (iSrc.eq.Idata) THEN
              Rname='rh_'//TRIM(srcList(i))//'_'//                      &
     &              TRIM(GridType (grSrc ))//'_'//                      &
     &              TRIM(GridType (grDst ))//'_'//                      &
     &              TRIM(IntrpType(Ibilin))//'_'//                      &
     &              TRIM(ExtrpType(etSrc ))//'_'//                      &
     &              TRIM(CplSetList(ic))//'_'//                         &
     &              TRIM(Cname)
            ELSE
              Rname='rh_'//                                             &
     &              TRIM(GridType (grSrc ))//'_'//                      &
     &              TRIM(GridType (grDst ))//'_'//                      &
     &              TRIM(IntrpType(Ibilin))//'_'//                      &
     &              TRIM(ExtrpType(etSrc ))//'_'//                      &
     &              TRIM(CplSetList(ic))//'_'//                         &
     &              TRIM(Cname)
            END IF
!
            CALL ESMF_StateGet (state,                                  &
     &                          TRIM(Rname),                            &
     &                          routeHandle,                            &
     &                          rc=rc)
            IF (ESMF_LogFoundError(rcToCheck=rc,                        &
     &                             msg=ESMF_LOGERR_PASSTHRU,            &
     &                             line=__LINE__,                       &
     &                             file=MyFile)) THEN
              RETURN
            END IF
!
!  Create temporary field in destination grid.
!
            tmpField=Coupler_FieldCreate (dstField,                     &
     &                                    Fname,                        &
     &                                    MISSING_dp, -1, rc)
            IF (ESMF_LogFoundError(rcToCheck=rc,                        &
     &                             msg=ESMF_LOGERR_PASSTHRU,            &
     &                             line=__LINE__,                       &
     &                             file=MyFile)) THEN
              RETURN
            END IF
!
!  Perform 1st REGRID operation.
!
            CALL ESMF_FieldRegrid (srcField,                            &
     &                             tmpField,                            &
     &                             routeHandle,                         &
     &                             zeroregion=ESMF_REGION_SELECT,       &
     &                             termorderflag=ESMF_TERMORDER_SRCSEQ, &
     &                             rc=rc)
            IF (ESMF_LogFoundError(rcToCheck=rc,                        &
     &                             msg=ESMF_LOGERR_PASSTHRU,            &
     &                             line=__LINE__,                       &
     &                             file=MyFile)) THEN
              RETURN
            END IF
!
!  Copy content from temporary field to destination field.
!
            CALL ESMF_FieldCopy (dstField,                              &
     &                           tmpField,                              &
     &                           rc=rc)
            IF (ESMF_LogFoundError(rcToCheck=rc,                        &
     &                             msg=ESMF_LOGERR_PASSTHRU,            &
     &                             line=__LINE__,                       &
     &                             file=MyFile)) THEN
              RETURN
            END IF
!
!  Get 2nd RouteHandle from state.
!
            IF (iSrc.eq.Idata) THEN
              Rname='rh_'//TRIM(srcList(i)) //'_'//                     &
     &               TRIM(GridType (grSrc ))//'_'//                     &
     &               TRIM(GridType (grDst ))//'_'//                     &
     &               TRIM(IntrpType(InStoD))//'_'//                     &
     &               TRIM(ExtrpType(etSrc ))//'_'//                     &
     &               TRIM(CplSetList(ic))//'_'//                        &
     &               TRIM(Cname)
            ELSE
              Rname='rh_'//                                             &
     &               TRIM(GridType(grSrc  ))//'_'//                     &
     &               TRIM(GridType(grDst  ))//'_'//                     &
     &               TRIM(IntrpType(InStoD))//'_'//                     &
     &               TRIM(ExtrpType(etSrc ))//'_'//                     &
     &               TRIM(CplSetList(ic))//'_'//                        &
     &               TRIM(Cname)
            END IF
!
            CALL ESMF_StateGet (state,                                  &
     &                          TRIM(Rname),                            &
     &                          routeHandle,                            &
     &                          rc=rc)
            IF (ESMF_LogFoundError(rcToCheck=rc,                        &
     &                             msg=ESMF_LOGERR_PASSTHRU,            &
     &                             line=__LINE__,                       &
     &                             file=MyFile)) THEN
              RETURN
            END IF
!
!  Perform 2nd REGRID operation to fill unmapped grid points.
!
            CALL ESMF_FieldRegrid (tmpField,                            &
     &                             dstField,                            &
     &                             routeHandle,                         &
     &                             zeroregion=ESMF_REGION_SELECT,       &
     &                             termorderflag=ESMF_TERMORDER_SRCSEQ, &
     &                             rc=rc)
            IF (ESMF_LogFoundError(rcToCheck=rc,                        &
     &                             msg=ESMF_LOGERR_PASSTHRU,            &
     &                             line=__LINE__,                       &
     &                             file=MyFile)) THEN
              RETURN
            END IF
!
            CALL ESMF_FieldDestroy (tmpField, rc=rc)
            IF (ESMF_LogFoundError(rcToCheck=rc,                        &
     &                             msg=ESMF_LOGERR_PASSTHRU,            &
     &                             line=__LINE__,                       &
     &                             file=MyFile)) THEN
              RETURN
            END IF
!
!  If integral adjustment is activated, calculate integral.
!
            IF (MODELS(iSrc)%ExportField(idSrc)%                        &
     &                                  enable_integral_adj) THEN
              src_total=0.0_dp
              src_total=Coupler_AreaIntegral(vm,                        &
     &                                       srcField,                  &
     &                                  (/UNMAPPED_MASK, MAPPED_MASK/), &
     &                                       rc)
              IF (ESMF_LogFoundError(rcToCheck=rc,                      &
     &                               msg=ESMF_LOGERR_PASSTHRU,          &
     &                               line=__LINE__,                     &
     &                               file=MyFile)) THEN
                RETURN
              END IF
!
              IF (localPET.eq.0) THEN
                WRITE (cplout,10)                                       &
     &                localPET, 'SRC. INTEGRAL', src_total,             &
     &                TRIM(MODELS(iSrc)%ExportField(idSrc)%short_name)
              END IF
!
              dst_total=0.0_dp
              dst_total=Coupler_AreaIntegral(vm,                        &
     &                                       dstField,                  &
     &                                  (/UNMAPPED_MASK, MAPPED_MASK/), &
     &                                       rc)
              IF (ESMF_LogFoundError(rcToCheck=rc,                      &
     &                               msg=ESMF_LOGERR_PASSTHRU,          &
     &                               line=__LINE__,                     &
     &                               file=MyFile)) THEN
                RETURN
              END IF
!
              IF (localPET.eq.0) THEN
                WRITE (cplout,10)                                       &
     &                localPET, 'DST. INTEGRAL', dst_total,             &
     &                TRIM(MODELS(iSrc)%ExportField(idSrc)%short_name)
                rel_error=0.0_dp
                IF (src_total.ne.0.0_dp) THEN
                  rel_error=(dst_total-src_total)/src_total
                END IF
                WRITE (cplout,10)                                       &
     &                localPET, 'RELATIVE ERROR 1', rel_error,          &
     &                TRIM(MODELS(iSrc)%ExportField(idSrc)%short_name)
              END IF
!
!  Adjust destination field based on calculated integral.
!
              CALL Coupler_AdjustField (vm,                             &
     &                                  dstField,                       &
     &                                  (/UNMAPPED_MASK, MAPPED_MASK/), &
     &                                  dst_total-src_total,            &
     &                                  rc)
              IF (ESMF_LogFoundError(rcToCheck=rc,                      &
     &                               msg=ESMF_LOGERR_PASSTHRU,          &
     &                               line=__LINE__,                     &
     &                               file=MyFile)) THEN
                RETURN
              END IF
!
              dst_total=0.0_dp
              dst_total=Coupler_AreaIntegral(vm,                        &
     &                                       dstField,                  &
     &                                  (/UNMAPPED_MASK, MAPPED_MASK/), &
     &                                       rc)
              IF (ESMF_LogFoundError(rcToCheck=rc,                      &
     &                               msg=ESMF_LOGERR_PASSTHRU,          &
     &                               line=__LINE__,                     &
     &                               file=MyFile)) THEN
                RETURN
              END IF
!
              IF (localPET.eq.0) THEN
                WRITE (cplout,10)                                       &
     &                localPET, 'DST. INTEGRAL (CORR)', dst_total,      &
     &                TRIM(MODELS(iSrc)%ExportField(idSrc)%short_name)
                rel_error=0.0_dp
                IF (src_total.ne.0.0_dp) THEN
                  rel_error=(dst_total-src_total)/src_total
                END IF
                WRITE (cplout,10)                                       &
     &                localPET, 'RELATIVE ERROR 2', rel_error,          &
     &                TRIM(MODELS(iSrc)%ExportField(idSrc)%short_name)
              END IF
            END IF
!
!-----------------------------------------------------------------------
!  Perform REGRID without extrapolation support.
!-----------------------------------------------------------------------
!
          ELSE
!
!  Get RouteHandle from state.
!
            IF (iSrc.eq.Idata) THEN
              Rname='rh_'//TRIM(srcList(i))//'_'//                      &
     &              TRIM(GridType  (grSrc))//'_'//                      &
     &              TRIM(GridType  (grDst))//'_'//                      &
     &              TRIM(IntrpType (itSrc))//'_'//                      &
     &              TRIM(ExtrpType (etSrc))//'_'//                      &
     &              TRIM(CplSetList(ic))//'_'//                         &
     &              TRIM(Cname)
            ELSE
              Rname='rh_'//                                             &
     &              TRIM(GridType (grSrc))//'_'//                       &
     &              TRIM(GridType (grDst))//'_'//                       &
     &              TRIM(IntrpType(itSrc))//'_'//                       &
     &              TRIM(ExtrpType(etSrc))//'_'//                       &
     &              TRIM(CplSetList(ic))//'_'//                         &
     &              TRIM(Cname)
            END IF
!
            CALL ESMF_StateGet (state,                                  &
     &                          TRIM(Rname),                            &
     &                          routeHandle,                            &
     &                          rc=rc)
            IF (ESMF_LogFoundError(rcToCheck=rc,                        &
     &                             msg=ESMF_LOGERR_PASSTHRU,            &
     &                             line=__LINE__,                       &
     &                             file=MyFile)) THEN
              RETURN
            END IF
!
!  Perform REGRID operation.
!
            CALL ESMF_FieldRegrid (srcField,                            &
     &                             dstField,                            &
     &                             routeHandle,                         &
     &                             zeroregion=ESMF_REGION_SELECT,       &
     &                             termorderflag=ESMF_TERMORDER_SRCSEQ, &
     &                             rc=rc)
            IF (ESMF_LogFoundError(rcToCheck=rc,                        &
     &                             msg=ESMF_LOGERR_PASSTHRU,            &
     &                             line=__LINE__,                       &
     &                             file=MyFile)) THEN
              RETURN
            END IF
!
          END IF QUERRY
!
!-----------------------------------------------------------------------
!  Debugging: report exchanged fields.
!-----------------------------------------------------------------------
!
          IF ((DebugLevel.gt.0).and.(localPET.eq.0)) THEN
            WRITE (cplout,20)                                           &
                TRIM(CplSetList(ic)), TRIM(Cname),                      &
     &          TRIM(MODELS(iSrc)%ExportField(idSrc)%short_name),       &
     &          TRIM(GridType (MODELS(iSrc)%ExportField(idSrc)%gtype)), &
     &          TRIM(MODELS(iDst)%ImportField(idDst)%short_name),       &
     &          TRIM(GridType (MODELS(iDst)%ImportField(idDst)%gtype)), &
     &          TRIM(IntrpType(MODELS(iSrc)%ExportField(idSrc)%itype))
          END IF
!
!-----------------------------------------------------------------------
!  Debugging: print out import/export fields time stamp information.
!-----------------------------------------------------------------------
!
          IF (DebugLevel.gt.2) THEN
            CALL NUOPC_GetTimeStamp (srcField,                          &
     &                               isValid = IsValid,                 &
     &                               time = srcTime,                    &
     &                               rc = rc)
            IF (ESMF_LogFoundError(rcToCheck=rc,                        &
     &                             msg=ESMF_LOGERR_PASSTHRU,            &
     &                             line=__LINE__,                       &
     &                             file=MyFile)) THEN
              RETURN
            END IF
!
            IF (IsValid) THEN
              CALL ESMF_TimeGet (srcTime,                               &
     &                           timeStringISOFrac = srcTimeString,     &
     &                           rc=rc)
              IF (ESMF_LogFoundError(rcToCheck=rc,                      &
     &                               msg=ESMF_LOGERR_PASSTHRU,          &
     &                               line=__LINE__,                     &
     &                               file=MyFile)) THEN
                RETURN
              END IF
              is=INDEX(srcTimeString, 'T')                  ! remove 'T'
              IF (is.gt.0) srcTimeString(is:is)=' '
              srcFileString=srcTimeString(1:19)
            ELSE
              srcFileString='0000-00-00 00:00:00'
            END IF
            srcFileString(11:11)='_'
            srcFileString(14:14)='.'
            srcFileString(17:17)='.'
!
            CALL NUOPC_GetTimeStamp (dstField,                          &
     &                               isValid = IsValid,                 &
     &                               time = dstTime,                    &
     &                               rc = rc)
            IF (ESMF_LogFoundError(rcToCheck=rc,                        &
     &                             msg=ESMF_LOGERR_PASSTHRU,            &
     &                             line=__LINE__,                       &
     &                             file=MyFile)) THEN
              RETURN
            END IF
!
            IF (IsValid) THEN
              CALL ESMF_TimeGet (dstTime,                               &
     &                           timeStringISOFrac = dstTimeString,     &
     &                           rc=rc)
              IF (ESMF_LogFoundError(rcToCheck=rc,                      &
     &                               msg=ESMF_LOGERR_PASSTHRU,          &
     &                               line=__LINE__,                     &
     &                               file=MyFile)) THEN
                RETURN
              END IF
              is=INDEX(dstTimeString, 'T')                  ! remove 'T'
              IF (is.gt.0) dstTimeString(is:is)=' '
              dstFileString=dstTimeString(1:19)
            ELSE
              dstFileString='0000-00-00 00:00:00'
            END IF
            dstFileString(11:11)='_'
            dstFileString(14:14)='.'
            dstFileString(17:17)='.'
!
            IF (localPET.eq.0) THEN
              WRITE (cplout,30) TRIM(srcTimeString),                    &
     &                          TRIM(dstTimeString),                    &
     &              TRIM(MODELS(iSrc)%ExportField(idSrc)%short_name),   &
     &              TRIM(MODELS(iDst)%ImportField(idDst)%short_name)
            END IF
          END IF
!
!-----------------------------------------------------------------------
!  Debugging: write data into NetCDF file. It uses the source field
!  time stamp in both source and destination NetCDF files for easy
!  matching of files. Usually, the time stamps between source and
!  destination fields is different.
!-----------------------------------------------------------------------
!
          IF ((DebugLevel.ge.3).and.                                    &
     &        MODELS(iSrc)%ExportField(idSrc)%debug_write) THEN
            WRITE (srcFile,40) 'src_'//TRIM(srcList(i))//'_'//          &
     &                         TRIM(CplSetList(ic))//'_'//              &
     &                         TRIM(lowercase(Cname)),                  &
     &                         TRIM(srcFileString)
            CALL ESMF_FieldWrite (srcField,                             &
     &                            TRIM(srcFile),                        &
     &                            variableName=TRIM(srcList(i)),        &
     &                            overwrite=.TRUE.,                     &
     &                            rc=rc)
            IF (ESMF_LogFoundError(rcToCheck=rc,                        &
     &                             msg=ESMF_LOGERR_PASSTHRU,            &
     &                             line=__LINE__,                       &
     &                             file=MyFile)) THEN
              RETURN
            END IF
          END IF
!
          IF ((DebugLevel.ge.3).and.                                    &
     &        MODELS(iDst)%ImportField(idDst)%debug_write) THEN
            WRITE (dstFile,40) 'dst_'//TRIM(dstList(i))//'_'//          &
     &                         TRIM(CplSetList(ic))//'_'//              &
     &                         TRIM(lowercase(Cname)),                  &
     &                         TRIM(srcFileString)
            CALL ESMF_FieldWrite (dstField,                             &
     &                            TRIM(dstFile),                        &
     &                            variableName=TRIM(dstList(i)),        &
     &                            overwrite=.TRUE.,                     &
     &                            rc=rc)
            IF (ESMF_LogFoundError(rcToCheck=rc,                        &
     &                             msg=ESMF_LOGERR_PASSTHRU,            &
     &                             line=__LINE__,                       &
     &                             file=MyFile)) THEN
              RETURN
            END IF
          END IF
!
        END DO EXCHANGE
!
!-----------------------------------------------------------------------
!  Deallocate temporary arrays
!-----------------------------------------------------------------------
!
        deallocate (srcList)
        deallocate (dstList)
      END DO CPLSET_LOOP
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '<== Exiting  Coupler_ExecuteRH for '   &
     &                          // TRIM(Cname), ', PET', PETrank
        CALL my_flush (trac)
      END IF
      IF (DebugLevel.gt.0) CALL my_flush (cplout)
!
 10   FORMAT (3x,'ESMF Coupler - PET(',i3.3,') - ',a,' = ',e14.5,       &
     &        ' (',a,')')
 20   FORMAT (3x,'ESMF Coupler - ',a,', ',a,': Regridded ',a,           &
     &        ' [',a,'] to ',a,' [',a,']',' >> ',a)
 30   FORMAT (18x,'(SRC TimeStamp = ',a,', DST TimeStamp = ',a,')',     &
     &              2x,a,' to ',a)
 40   FORMAT (a,'_',a,'.nc')
!
      RETURN
      END SUBROUTINE Coupler_ExecuteRH
!
      SUBROUTINE Coupler_ReleaseRH (coupler, rc)
!
!=======================================================================
!                                                                      !
!  Releases coupler RouteHandle connectors between source (srcFields)  !
!  and destination (dstFields) fields.                                 !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer, intent(out) :: rc
!
      TYPE (ESMF_CplComp) :: coupler
!
!  Local variable declarations.
!
      logical :: rhExist, rh1Exist, rh2Exist
!
      integer :: i, ic, j, localPET, PETcount, MyComm
      integer :: itemCount, srcCount, dstCount, NcplSets
      integer :: iSrc, iDst, idSrc, idDst, grSrc, grDst
      integer :: etSrc, etDst, itSrc, itDst
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", Coupler_ReleaseRH"

      character(ESMF_MAXSTR) :: Cname, Rname

      character (ESMF_MAXSTR), pointer :: CplSetList(:) => NULL()
      character (ESMF_MAXSTR), pointer :: dstList(:) => NULL()
      character (ESMF_MAXSTR), pointer :: srcList(:) => NULL()
!
      TYPE (ESMF_VM)          :: vm
      TYPE (ESMF_State)       :: state
      TYPE (ESMF_FieldBundle) :: srcFields, dstFields
      TYPE (ESMF_RouteHandle) :: routeHandle
!
!-----------------------------------------------------------------------
!  Initialize return code flag to success state (no error).
!-----------------------------------------------------------------------
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '==> Entering Coupler_ReleaseRH',       &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
      rc=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!  Query coupler component.
!-----------------------------------------------------------------------
!
!  Querry the coupler for the Virtual Machine (VM) parallel environmemt.
!
      CALL ESMF_CplCompGet (coupler,                                    &
     &                      name=Cname,                                 &
     &                      vm=vm,                                      &
     &                      rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!  Get current parallel node rank and number of nodes.
!
      CALL ESMF_VMGet (vm,                                              &
     &                 localPet=localPET,                               &
     &                 petCount=PETcount,                               &
     &                 mpiCommunicator=MyComm,                          &
     &                 rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!  Set source and destination couple model indices.
!
      DO i=1,Nmodels
        DO j=1,Nmodels
          IF ((CONNECTORS(i,j)%IsActive).and.                           &
     &        (TRIM(CONNECTORS(i,j)%name).eq.TRIM(Cname))) THEN
            iSrc=i
            iDst=j
          END if
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Get coupled set list for connector (Cname).
!-----------------------------------------------------------------------
!
      IF ( associated(CplSetList) ) nullify (CplSetList)
      CALL NUOPC_ConnectorGet (coupler,                                 &
     &                         cplSetList=CplSetList,                   &
     &                         rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
      NcplSets=SIZE(CplSetList)
!
!-----------------------------------------------------------------------
!  Inquire about source and destination fields.
!-----------------------------------------------------------------------
!
      CPLSET_LOOP : DO ic=1,NcplSets
!
!  Get source and destination fields for each coupled set.
!
        CALL NUOPC_ConnectorGet (coupler,                               &
     &                           srcFields=srcFields,                   &
     &                           dstFields=dstFields,                   &
     &                           state=state,                           &
     &                           cplSet=CplSetList(ic),                 &
     &                           rc=rc)
        IF (ESMF_LogFoundError(rcToCheck=rc,                            &
     &                         msg=ESMF_LOGERR_PASSTHRU,                &
     &                         line=__LINE__,                           &
     &                         file=MyFile)) THEN
          RETURN
        END IF
!
!  Number of source fields.
!
        CALL ESMF_FieldBundleGet (srcFields,                            &
     &                            fieldCount=srcCount,                  &
     &                            rc=rc)
        IF (ESMF_LogFoundError(rcToCheck=rc,                            &
     &                         msg=ESMF_LOGERR_PASSTHRU,                &
     &                         line=__LINE__,                           &
     &                         file=MyFile)) THEN
          RETURN
        END IF
!
!  Number of destination fields.
!
        CALL ESMF_FieldBundleGet (dstFields,                            &
     &                            fieldCount=dstCount,                  &
     &                            rc=rc)
        IF (ESMF_LogFoundError(rcToCheck=rc,                            &
     &                         msg=ESMF_LOGERR_PASSTHRU,                &
     &                         line=__LINE__,                           &
     &                         file=MyFile)) THEN
          RETURN
        END IF
!
!  Source fields names.
!
        allocate ( srcList(srcCount) )

        CALL ESMF_FieldBundleGet (srcFields,                            &
     &                            fieldNameList=srcList,                &
     &                            rc=rc)
        IF (ESMF_LogFoundError(rcToCheck=rc,                            &
     &                         msg=ESMF_LOGERR_PASSTHRU,                &
     &                         line=__LINE__,                           &
     &                         file=MyFile)) THEN
          RETURN
        END IF
!
!  Destination fields names.
!
        allocate ( dstList(dstCount) )

        CALL ESMF_FieldBundleGet (dstFields,                            &
     &                            fieldNameList=dstList,                &
     &                            rc=rc)
        IF (ESMF_LogFoundError(rcToCheck=rc,                            &
     &                         msg=ESMF_LOGERR_PASSTHRU,                &
     &                         line=__LINE__,                           &
     &                         file=MyFile)) THEN
          RETURN
        END IF
!
!=======================================================================
!  Release coupling connector RouteHandle between source and destination
!  fields.
!=======================================================================
!
        EXCHANGE : DO i=1,srcCount
!
!  Set source and destination field index.
!
          idSrc=field_index(MODELS(iSrc)%ExportField, srcList(i))
          idDst=field_index(MODELS(iDst)%ImportField, dstList(i))
!
!  Set interpolation type.
!
          itSrc=MODELS(iSrc)%ExportField(idSrc)%itype
          itDst=MODELS(iDst)%ImportField(idDst)%itype
!
!  Get extrapolation method for unmapped destination points.
!
          etSrc=MODELS(iSrc)%ExportField(idSrc)%etype
          etDst=MODELS(iDst)%ImportField(idDst)%etype
!
!  Set grid type.
!
          grSrc=MODELS(iSrc)%ExportField(idSrc)%gtype
          grDst=MODELS(iDst)%ImportField(idDst)%gtype
!
!-----------------------------------------------------------------------
!  Release RouteHandle for REGRID with two-steps extrapolation.
!-----------------------------------------------------------------------
!
          QUERRY : IF (etSrc.eq.E2steps) THEN
!
!  Check 1st RouteHandle (i.e. rh_Center_Corner_BLIN_ATM-OCN).
!
            IF (iSrc.eq.Idata) THEN
              Rname='rh_'//TRIM(srcList(i))//'_'//                      &
     &              TRIM(GridType (grSrc ))//'_'//                      &
     &              TRIM(GridType (grDst ))//'_'//                      &
     &              TRIM(IntrpType(Ibilin))//'_'//                      &
     &              TRIM(ExtrpType(etSrc ))//'_'//                      &
     &              TRIM(CplSetList(ic))//'_'//                         &
     &              TRIM(Cname)
            ELSE
              Rname='rh_'//                                             &
     &              TRIM(GridType (grSrc ))//'_'//                      &
     &              TRIM(GridType (grDst ))//'_'//                      &
     &              TRIM(IntrpType(Ibilin))//'_'//                      &
     &              TRIM(ExtrpType(etSrc ))//'_'//                      &
     &              TRIM(CplSetList(ic))//'_'//                         &
     &              TRIM(Cname)
            END IF
!
            CALL ESMF_StateGet (state,                                  &
     &                          itemSearch=TRIM(Rname),                 &
     &                          itemCount=itemCount,                    &
     &                          rc=rc)
            IF (ESMF_LogFoundError(rcToCheck=rc,                        &
     &                             msg=ESMF_LOGERR_PASSTHRU,            &
     &                             line=__LINE__,                       &
     &                             file=MyFile)) THEN
              RETURN
            END IF
!
            IF (itemCount.le.0) THEN
              rh1Exist=.FALSE.
            ELSE
              rh1Exist=.TRUE.
            END IF
!
!  Release 1st RouteHandle.
!
            IF (rh1Exist) THEN
              CALL ESMF_StateGet (state,                                &
     &                            TRIM(Rname),                          &
     &                            routeHandle,                          &
     &                            rc=rc)
              IF (ESMF_LogFoundError(rcToCheck=rc,                      &
     &                               msg=ESMF_LOGERR_PASSTHRU,          &
     &                               line=__LINE__,                     &
     &                               file=MyFile)) THEN
                RETURN
              END IF
!
              CALL ESMF_FieldBundleRegridRelease (routeHandle,          &
     &                                            rc=rc)
              IF (ESMF_LogFoundError(rcToCheck=rc,                      &
     &                               msg=ESMF_LOGERR_PASSTHRU,          &
     &                               line=__LINE__,                     &
     &                               file=MyFile)) THEN
                RETURN
              END IF
            END IF
!
!  Check 2nd RouteHandle (i.e. rh_Center_Corner_NS2D_ATM-OCN).
!
            IF (iSrc.eq.Idata) THEN
              Rname='rh_'//TRIM(srcList(i))//'_'//                      &
     &              TRIM(GridType (grSrc ))//'_'//                      &
     &              TRIM(GridType (grDst ))//'_'//                      &
     &              TRIM(IntrpType(InStoD))//'_'//                      &
     &              TRIM(ExtrpType(etSrc ))//'_'//                      &
     &              TRIM(CplSetList(ic))//'_'//                         &
     &              TRIM(Cname)
            ELSE
              Rname='rh_'//                                             &
     &              TRIM(GridType (grSrc ))//'_'//                      &
     &              TRIM(GridType (grDst ))//'_'//                      &
     &              TRIM(IntrpType(InStoD))//'_'//                      &
     &              TRIM(ExtrpType(etSrc ))//'_'//                      &
     &              TRIM(CplSetList(ic))//'_'//                         &
     &              TRIM(Cname)
            END IF
!
            CALL ESMF_StateGet (state,                                  &
     &                          itemSearch=TRIM(Rname),                 &
     &                          itemCount=itemCount,                    &
     &                          rc=rc)
            IF (ESMF_LogFoundError(rcToCheck=rc,                        &
     &                             msg=ESMF_LOGERR_PASSTHRU,            &
     &                             line=__LINE__,                       &
     &                             file=MyFile)) THEN
              RETURN
            END IF
!
            IF (itemCount.le.0) THEN
              rh2Exist=.FALSE.
            ELSE
              rh2Exist=.TRUE.
            END IF
!
!  Release 2nd RouteHandle.
!
            IF (rh2Exist) THEN
              CALL ESMF_StateGet (state,                                &
     &                            TRIM(Rname),                          &
     &                            routeHandle,                          &
     &                            rc=rc)
              IF (ESMF_LogFoundError(rcToCheck=rc,                      &
     &                               msg=ESMF_LOGERR_PASSTHRU,          &
     &                               line=__LINE__,                     &
     &                               file=MyFile)) THEN
                RETURN
             END IF
!
              CALL ESMF_FieldBundleRegridRelease (routeHandle,          &
     &                                            rc=rc)
              IF (ESMF_LogFoundError(rcToCheck=rc,                      &
     &                               msg=ESMF_LOGERR_PASSTHRU,          &
     &                               line=__LINE__,                     &
     &                               file=MyFile)) THEN
                RETURN
              END IF
            END IF
!
!-----------------------------------------------------------------------
!  Release RouteHandle for REGRID without extrapolation support.
!-----------------------------------------------------------------------
!
          ELSE
!
!  Check RouteHandle for one step interpolation.
!
            IF (iSrc.eq.Idata) THEN
              Rname='rh_'//TRIM(srcList(i))//'_'//                      &
     &              TRIM(GridType  (grSrc))//'_'//                      &
     &              TRIM(GridType  (grDst))//'_'//                      &
     &              TRIM(IntrpType (itSrc))//'_'//                      &
     &              TRIM(ExtrpType (etSrc))//'_'//                      &
     &              TRIM(CplSetList(ic))//'_'//                         &
     &              TRIM(Cname)
            ELSE
              Rname='rh_'//                                             &
     &              TRIM(GridType  (grSrc))//'_'//                      &
     &              TRIM(GridType  (grDst))//'_'//                      &
     &              TRIM(IntrpType (itSrc))//'_'//                      &
     &              TRIM(ExtrpType (etSrc))//'_'//                      &
     &              TRIM(CplSetList(ic))//'_'//                         &
     &              TRIM(Cname)
            END IF
!
            CALL ESMF_StateGet (state,                                  &
     &                          itemSearch=TRIM(Rname),                 &
     &                          itemCount=itemCount,                    &
     &                          rc=rc)
            IF (ESMF_LogFoundError(rcToCheck=rc,                        &
     &                             msg=ESMF_LOGERR_PASSTHRU,            &
     &                             line=__LINE__,                       &
     &                             file=MyFile)) THEN
              RETURN
            END IF
!
            IF (itemCount.le.0) THEN
              rhExist=.FALSE.
            ELSE
              rhExist=.TRUE.
            END IF
!
!  Release RouteHandle.
!
            IF (rhExist) THEN
              CALL ESMF_StateGet (state,                                &
     &                            TRIM(Rname),                          &
     &                            routeHandle,                          &
     &                            rc=rc)
              IF (ESMF_LogFoundError(rcToCheck=rc,                      &
     &                               msg=ESMF_LOGERR_PASSTHRU,          &
     &                               line=__LINE__,                     &
     &                               file=MyFile)) THEN
                RETURN
              END IF
!
              CALL ESMF_FieldBundleRegridRelease (routeHandle,          &
     &                                            rc=rc)
              IF (ESMF_LogFoundError(rcToCheck=rc,                      &
     &                               msg=ESMF_LOGERR_PASSTHRU,          &
     &                               line=__LINE__,                     &
     &                               file=MyFile)) THEN
                RETURN
              END IF
            END IF
!
          END IF QUERRY
        END DO EXCHANGE
      END DO CPLSET_LOOP
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '<== Exiting  Coupler_ReleaseRH',       &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
!
      RETURN
      END SUBROUTINE Coupler_ReleaseRH
!
      SUBROUTINE Coupler_AdjustField (vm, field, maskval, error, rc)
!
!=======================================================================
!                                                                      !
!  Adjusts regrid field to ensure global area integral conservation    !
!  over the matched regions. The destination field is adjusted using   !
!  the global error of the area integrated difference between          !
!  destination and source field over the matched regions:              !
!                                                                      !
!    field = field - error/Area                                        !
!                                                                      !
!  where                                                               !
!                                                                      !
!    error = SumDstArea - SumSrcArea                                   !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer (i4b), intent(in) :: maskval(:)
      integer, intent(out) :: rc
!
      real (dp), intent(in) :: error
!
      TYPE (ESMF_VM), intent(in)       :: vm
      TYPE (ESMF_Field), intent(inout) :: field
!
!  Local variable declarations.
!
      integer :: i, j
      integer :: localDE, localDEcount, localPET, PETcount, MyComm
      integer :: cLbnd(2), cUbnd(2)
!
      integer (i4b), pointer :: ptrMask(:,:) => NULL()
!
      real (dp) :: MyAreaSum(1), AreaSum(1)
      real (dp) :: error_unit
!
      real (dp), pointer :: ptrField(:,:) => NULL()
      real (dp), pointer :: ptrArea(:,:) => NULL()
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", Coupler_AdjustedField"

      character(ESMF_MAXSTR) :: Fname
!
      TYPE (ESMF_Grid)       :: grid
      TYPE (ESMF_StaggerLoc) :: sLoc
!
!-----------------------------------------------------------------------
!  Initialize.
!-----------------------------------------------------------------------
!
!  Set return code flag to success state (no error).
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '==> Entering Coupler_AdjustField',     &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
      rc=ESMF_SUCCESS
!
!  Integral.
!
      MyAreaSum(1)=0.0_dp
!
!-----------------------------------------------------------------------
!  Querry the Virtual Machine (VM) parallel environmemt for the mpi
!  communicator handle and current node rank.
!-----------------------------------------------------------------------
!
      CALL ESMF_VMGet (vm,                                              &
     &                 localPet=localPET,                               &
     &                 petCount=PETcount,                               &
     &                 mpiCommunicator=MyComm,                          &
     &                 rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!-----------------------------------------------------------------------
!  Query input field.
!-----------------------------------------------------------------------
!
      CALL ESMF_FieldGet (field,                                        &
     &                    grid=grid,                                    &
     &                    name=Fname,                                   &
     &                    staggerloc=sLoc,                              &
     &                    rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!  Get number of local decomposition elements (DEs) in the grid.
!
      CALL ESMF_GridGet (grid,                                          &
     &                   localDECount=localDEcount,                     &
     &                   rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!  Get field pointer.
!
      DE_LOOP1 : DO localDE=0,localDEcount-1
        CALL ESMF_FieldGet (field,                                      &
     &                      localDe=localDE,                            &
     &                      farrayPtr=ptrField,                         &
     &                      computationalLBound=cLbnd,                  &
     &                      computationalUBound=cUbnd,                  &
     &                      rc=rc)
        IF (ESMF_LogFoundError(rcToCheck=rc,                            &
     &                         msg=ESMF_LOGERR_PASSTHRU,                &
     &                         line=__LINE__,                           &
     &                         file=MyFile)) THEN
          RETURN
        END IF
!
!  Get area pointer from grid.
!
        CALL ESMF_GridGetItem (grid,                                    &
     &                         ESMF_GRIDITEM_AREA,                      &
     &                         staggerloc=sLoc,                         &
     &                         localDe=localDE,                         &
     &                         farrayPtr=ptrArea,                       &
     &                         rc=rc)
        IF (ESMF_LogFoundError(rcToCheck=rc,                            &
     &                         msg=ESMF_LOGERR_PASSTHRU,                &
     &                         line=__LINE__,                           &
     &                         file=MyFile)) THEN
          RETURN
        END IF
!
!  Get mask pointer from grid.
!
        CALL ESMF_GridGetItem (grid,                                    &
     &                         ESMF_GRIDITEM_MASK,                      &
     &                         staggerloc=sLoc,                         &
     &                         localDe=localDE,                         &
     &                         farrayPtr=ptrMask,                       &
     &                         rc=rc)
        IF (ESMF_LogFoundError(rcToCheck=rc,                            &
     &                         msg=ESMF_LOGERR_PASSTHRU,                &
     &                         line=__LINE__,                           &
     &                         file=MyFile)) THEN
          RETURN
        END IF
!
!-----------------------------------------------------------------------
!  Calculate total area of matched region
!-----------------------------------------------------------------------
!
        cLbnd(1)=LBOUND(ptrMask, dim=1)
        cUbnd(1)=UBOUND(ptrMask, dim=1)
        cLbnd(2)=LBOUND(ptrMask, dim=2)
        cUbnd(2)=UBOUND(ptrMask, dim=2)
!
        DO j=cLbnd(2),cUbnd(2)
          DO i=cLbnd(1),cUbnd(1)
            IF (ANY(ptrMask(i,j).eq.maskval)) THEN
              MyAreaSum(1)=MyAreaSum(1)+ptrArea(i,j)
            END IF
          END DO
        END DO
!
!-----------------------------------------------------------------------
!  Nullify pointer to make sure that it does not point on a random
!  part in the memory.
!-----------------------------------------------------------------------
!
        IF (associated(ptrArea)) THEN
          nullify (ptrArea)
        END IF
        IF (associated(ptrMask)) THEN
          nullify (ptrMask)
        END IF
      END DO DE_LOOP1
!
!-----------------------------------------------------------------------
!  Collect calculated total area from PETs
!-----------------------------------------------------------------------
!
      AreaSum(1)=0.0_dp
      CALL ESMF_VMAllReduce (vm,                                        &
     &                       MyAreaSum,                                 &
     &                       AreaSum, 1,                                &
     &                       ESMF_REDUCE_SUM,                           &
     &                       rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!-----------------------------------------------------------------------
!  Calculate error per unit area.
!-----------------------------------------------------------------------
!
      error_unit=error/AreaSum(1)
      IF (localPET.eq.0) THEN
        WRITE (cplout,10) localPET, AreaSum(1), error_unit, TRIM(Fname)
 10     FORMAT (' PET(',i3.3,') - AVGERAGE DIFF = ',2e14.5,' (',a,')')
      END IF
!
!-----------------------------------------------------------------------
!  Adjust field using the global error of the area integrated
!  difference between source and destination field.
!-----------------------------------------------------------------------
!
      DE_LOOP2 : DO localDE=0,localDEcount-1
!
!  Get field pointers
!
        CALL ESMF_FieldGet (field,                                      &
     &                      localDe=localDE,                            &
     &                      farrayPtr=ptrField,                         &
     &                      computationalLBound=cLbnd,                  &
     &                      computationalUBound=cUbnd,                  &
     &                      rc=rc)
        IF (ESMF_LogFoundError(rcToCheck=rc,                            &
     &                         msg=ESMF_LOGERR_PASSTHRU,                &
     &                         line=__LINE__,                           &
     &                         file=MyFile)) THEN
          RETURN
        END IF
!
!  Get mask pointer from grid.
!
        CALL ESMF_GridGetItem (grid,                                    &
     &                         ESMF_GRIDITEM_MASK,                      &
     &                         staggerloc=sLoc,                         &
     &                         localDe=localDE,                         &
     &                         farrayPtr=ptrMask,                       &
     &                         rc=rc)
        IF (ESMF_LogFoundError(rcToCheck=rc,                            &
     &                         msg=ESMF_LOGERR_PASSTHRU,                &
     &                         line=__LINE__,                           &
     &                         file=MyFile)) THEN
          RETURN
        END IF
!
!  Adjust input destination field.
!
        DO j=cLbnd(2),cUbnd(2)
          DO i=cLbnd(1),cUbnd(1)
            IF (ANY(ptrMask(i,j).eq.maskval)) THEN
              ptrField(i,j)=ptrField(i,j)-error_unit
            END IF
          END DO
        END DO
!
!-----------------------------------------------------------------------
!  Nullify pointer to make sure that it does not point on a random
!  part in the memory.
!-----------------------------------------------------------------------
!
        IF (associated(ptrField)) THEN
          nullify (ptrField)
        END IF
        IF (associated(ptrMask)) THEN
          nullify (ptrMask)
        END IF
      END DO DE_LOOP2
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '<== Exiting  Coupler_AdjustField',     &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
!
      RETURN
      END SUBROUTINE Coupler_AdjustField
!
      FUNCTION Coupler_AreaIntegral (vm, field, maskval, rc)
!
!=======================================================================
!                                                                      !
!  Computes source or destination field global area integral over the  !
!  matched regions (cells having the specified mask values, maskval).  !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer (i4b), intent(in) :: maskval(:)
      integer, intent(out) :: rc
!
      real (dp) :: Coupler_AreaIntegral
!
      TYPE (ESMF_Field), intent(in) :: field
      TYPE (ESMF_VM), intent(in)    :: vm
!
!  Local variable declarations.
!
      integer :: i, j
      integer :: localDE, localDEcount, localPET, PETcount, MyComm
      integer :: cLbnd(2), cUbnd(2)
!
      integer (i4b), pointer :: ptrMask(:,:) => NULL()
!
      real (dp) :: MyAreaSum(1), AreaSum(1)
!
      real (dp), pointer :: ptrField(:,:) => NULL()
      real (dp), pointer :: ptrArea(:,:) => NULL()
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", Coupler_AreaIntegral"

      character (ESMF_MAXSTR) :: Fname
!
      TYPE (ESMF_Grid)       :: grid
      TYPE (ESMF_StaggerLoc) :: sLoc
!
!-----------------------------------------------------------------------
!  Initialize.
!-----------------------------------------------------------------------
!
!  Set return code flag to success state (no error).
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '==> Entering Coupler_AreaIntegral',    &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
      rc=ESMF_SUCCESS
!
!  Integral.
!
      Coupler_AreaIntegral=0.0_dp
      MyAreaSum(1)=0.0_dp
!
!-----------------------------------------------------------------------
!  Querry the Virtual Machine (VM) parallel environmemt for the mpi
!  communicator handle and current node rank.
!-----------------------------------------------------------------------
!
      CALL ESMF_VMGet (vm,                                              &
     &                 localPet=localPET,                               &
     &                 petCount=PETcount,                               &
     &                 mpiCommunicator=MyComm,                          &
     &                 rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!-----------------------------------------------------------------------
!  Query input field.
!-----------------------------------------------------------------------
!
      CALL ESMF_FieldGet (field,                                        &
     &                    grid=grid,                                    &
     &                    name=Fname,                                   &
     &                    staggerloc=sLoc,                              &
     &                    rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!  Get number of local decomposition elements (DEs) in the grid.
!
      CALL ESMF_GridGet (grid,                                          &
     &                   localDECount=localDEcount,                     &
     &                   rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!  Get field pointer.
!
      DE_LOOP : DO localDE=0,localDEcount-1
        CALL ESMF_FieldGet (field,                                      &
     &                      localDe=localDE,                            &
     &                      farrayPtr=ptrField,                         &
     &                      computationalLBound=cLbnd,                  &
     &                      computationalUBound=cUbnd,                  &
     &                      rc=rc)
        IF (ESMF_LogFoundError(rcToCheck=rc,                            &
     &                         msg=ESMF_LOGERR_PASSTHRU,                &
     &                         line=__LINE__,                           &
     &                         file=MyFile)) THEN
          RETURN
        END IF
!
!  Get area pointer from grid.
!
        CALL ESMF_GridGetItem (grid,                                    &
     &                         ESMF_GRIDITEM_AREA,                      &
     &                         staggerloc=sLoc,                         &
     &                         localDe=localDE,                         &
     &                         farrayPtr=ptrArea,                       &
     &                         rc=rc)
        IF (ESMF_LogFoundError(rcToCheck=rc,                            &
     &                         msg=ESMF_LOGERR_PASSTHRU,                &
     &                         line=__LINE__,                           &
     &                         file=MyFile)) THEN
          RETURN
        END IF
!
!  Get mask pointer from grid.
!
        CALL ESMF_GridGetItem (grid,                                    &
     &                         ESMF_GRIDITEM_MASK,                      &
     &                         staggerloc=sLoc,                         &
     &                         localDe=localDE,                         &
     &                         farrayPtr=ptrMask,                       &
     &                         rc=rc)
        IF (ESMF_LogFoundError(rcToCheck=rc,                            &
     &                         msg=ESMF_LOGERR_PASSTHRU,                &
     &                         line=__LINE__,                           &
     &                         file=MyFile)) THEN
          RETURN
        END IF
!
!-----------------------------------------------------------------------
!  Calculate field area integral for each local DE and PET.
!-----------------------------------------------------------------------
!
        DO j=cLbnd(2),cUbnd(2)
          DO i=cLbnd(1),cUbnd(1)
            IF (ANY(ptrMask(i,j).eq.maskval)) THEN
              MyAreaSum(1)=MyAreaSum(1)+ptrField(i,j)*ptrArea(i,j)
            END IF
          END DO
        END DO
!
!-----------------------------------------------------------------------
!  Nullify pointer to make sure that it does not point on a random
!  part in the memory.
!-----------------------------------------------------------------------
!
        IF (associated(ptrField)) THEN
          nullify (ptrField)
        END IF
        IF (associated(ptrArea)) THEN
          nullify (ptrArea)
        END IF
        IF (associated(ptrMask)) THEN
          nullify  (ptrMask)
        END IF
      END DO DE_LOOP
!
!-----------------------------------------------------------------------
!  Debugging: write sum of each PETs
!-----------------------------------------------------------------------
!
      IF (DebugLevel.gt.2) THEN
        WRITE (cplout,10) localPET, localDE, MyAreaSum(1), TRIM(Fname)
 10     FORMAT (' PET(',i3.3,') - DE(',i2.2,') - Area Integral = ',     &
     &          e14.5,' (',a,')')
        CALL ESMF_VMBarrier (vm, rc=rc)
        IF (ESMF_LogFoundError(rcToCheck=rc,                            &
     &                         msg=ESMF_LOGERR_PASSTHRU,                &
     &                         line=__LINE__,                           &
     &                         file=MyFile)) THEN
          RETURN
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Collect fiels area integral from each PET and calculate global value.
!-----------------------------------------------------------------------
!
      CALL ESMF_VMAllReduce (vm,                                        &
     &                       MyAreaSum, AreaSum, 1,                     &
     &                       ESMF_REDUCE_SUM,                           &
     &                       rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
      Coupler_AreaIntegral=AreaSum(1)
!
!-----------------------------------------------------------------------
!  Debugging: report global area integral.
!-----------------------------------------------------------------------
!
      IF (DebugLevel.eq.1) THEN
        IF (localPET.eq.0) THEN
          WRITE (cplout,20) localPET, AreaSum(1), TRIM(Fname)
 20       FORMAT (' PET(',I3.3,') - Global Area Integral = ',e14.5,     &
     &            ' (',a,')')
        END IF
        CALL ESMF_VMBarrier (vm, rc=rc)
        IF (ESMF_LogFoundError(rcToCheck=rc,                            &
     &                         msg=ESMF_LOGERR_PASSTHRU,                &
     &                         line=__LINE__,                           &
     &                         file=MyFile)) THEN
          RETURN
        END IF
      END IF
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '<== Exiting  Coupler_AreaIntegral',    &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
      IF (DebugLevel.gt.0) CALL my_flush (cplout)
!
      RETURN
      END FUNCTION Coupler_AreaIntegral
!
      FUNCTION Coupler_FieldCreate (field, Fname, IniVal, dstLandMask,  &
     &                              rc)
!
!=======================================================================
!                                                                      !
!  Creates and initializes a new field using input field attributes.   !
!  Masked grid cells are set to missing value (MISSING_dp).            !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer (i4b), intent(in) :: dstLandMask
      integer, intent(out) :: rc
!
      real (dp), intent(in) :: IniVal
!
      character (len=*), intent(in) :: Fname
!
      TYPE (ESMF_Field), intent(in) :: field
      TYPE (ESMF_Field) :: Coupler_FieldCreate
!
!  Local variable declarations.
!
      integer :: i, j, localDE, localDEcount
      integer :: cLbnd(2), cUbnd(2)
!
      integer (i4b), pointer :: msk2d(:,:) => NULL()
      integer (i4b), allocatable :: tlw(:,:), tuw(:,:)
!
      real (dp), pointer :: ptr2d(:,:) => NULL()
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", Coupler_FieldCreate"
!
      TYPE (ESMF_Grid)       :: grid
      TYPE (ESMF_DistGrid)   :: distGrid
      TYPE (ESMF_ArraySpec)  :: arraySpec
      TYPE (ESMF_StaggerLoc) :: staggerLoc
!
!-----------------------------------------------------------------------
!  Initialize return code flag to success state (no error).
!-----------------------------------------------------------------------
!
      rc=ESMF_SUCCESS
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '==> Entering Coupler_FieldCreate',     &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
!
!-----------------------------------------------------------------------
!  Query input field.
!-----------------------------------------------------------------------
!
      CALL ESMF_FieldGet (field,                                        &
     &                    arrayspec=arraySpec,                          &
     &                    grid=grid,                                    &
     &                    staggerloc=staggerLoc,                        &
     &                    rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!  Get information about field associated grid.
!
      CALL ESMF_GridGet (grid,                                          &
     &                   distgrid=distGrid,                             &
     &                   localDECount=localDEcount,                     &
     &                   rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!  Get lower and upper bound of halo region.
!
      IF (.not.allocated(tlw)) THEN
        allocate ( tlw(2,localDEcount) )
      END IF
      IF (.not.allocated(tuw)) THEN
        allocate ( tuw(2,localDEcount) )
      END IF
!
      CALL ESMF_FieldGet (field,                                        &
     &                    totalLWidth=tlw,                              &
     &                    totalUWidth=tuw,                              &
     &                    rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!-----------------------------------------------------------------------
!  Create new field using input field attributes.
!-----------------------------------------------------------------------
!
      IF (localDEcount.eq.1) THEN
        Coupler_FieldCreate=ESMF_FieldCreate(grid,                      &
     &                                       arraySpec,                 &
     &                                       staggerloc=staggerLoc,     &
     &                                       totalLWidth=tlw(:,1),      &
     &                                       totalUWidth=tuw(:,1),      &
     &                                       name=TRIM(Fname),          &
     &                                       rc=rc)
        IF (ESMF_LogFoundError(rcToCheck=rc,                            &
     &                         msg=ESMF_LOGERR_PASSTHRU,                &
     &                         line=__LINE__,                           &
     &                         file=MyFile)) THEN
          RETURN
        END IF
      ELSE
        Coupler_FieldCreate=ESMF_FieldCreate(grid,                      &
     &                                       arraySpec,                 &
     &                                       staggerloc=staggerLoc,     &
     &                                       name=TRIM(Fname),          &
     &                                       rc=rc)
        IF (ESMF_LogFoundError(rcToCheck=rc,                            &
     &                         msg=ESMF_LOGERR_PASSTHRU,                &
     &                         line=__LINE__,                           &
     &                         file=MyFile)) THEN
          RETURN
        END IF
      END IF
!
!  Get pointer from new field.

      DE_LOOP : DO localDE=0,localDEcount-1
        CALL ESMF_FieldGet (Coupler_FieldCreate,                        &
     &                      localDe=localDE,                            &
     &                      farrayPtr=ptr2d,                            &
     &                      computationalLBound=cLbnd,                  &
     &                      computationalUBound=cUbnd,                  &
     &                      rc=rc)
        IF (ESMF_LogFoundError(rcToCheck=rc,                            &
     &                         msg=ESMF_LOGERR_PASSTHRU,                &
     &                         line=__LINE__,                           &
     &                         file=MyFile)) THEN
          RETURN
        END IF
!
!  Get mask pointer from grid.
!
        CALL ESMF_GridGetItem (grid,                                    &
     &                         ESMF_GRIDITEM_MASK,                      &
     &                         staggerloc=staggerLoc,                   &
     &                         localDe=localDE,                         &
     &                         farrayPtr=msk2d,                         &
     &                         rc=rc)
        IF (ESMF_LogFoundError(rcToCheck=rc,                            &
     &                         msg=ESMF_LOGERR_PASSTHRU,                &
     &                         line=__LINE__,                           &
     &                         file=MyFile)) THEN
          RETURN
        END IF
!
!  Initialize pointer to new field. Masked grid cells are set to
!  missing value.
!
        DO j=cLbnd(2),cUbnd(2)
          DO i=cLbnd(1),cUbnd(1)
            IF (msk2d(i,j).ne.dstLandMask) THEN
              ptr2d(i,j)=IniVal
            ELSE
              ptr2d(i,j)=MISSING_dp
            END IF
          END DO
        END DO
!
!  Nullify pointers to make sure that it does not point to a random
!  part in the memory.
!
        IF (associated(ptr2d)) THEN
          nullify (ptr2d)
        END IF
        IF (associated(msk2d)) THEN
          nullify (msk2d)
        END IF
      END DO DE_LOOP
!
!-----------------------------------------------------------------------
!  Deallocate local arrays.
!-----------------------------------------------------------------------
!
      IF (allocated(tlw)) THEN
        deallocate (tlw)
      END IF
      IF (allocated(tuw)) THEN
        deallocate (tuw)
      END IF
!
!-----------------------------------------------------------------------
!  Check consistency of the created field.
!-----------------------------------------------------------------------
!
      CALL ESMF_FieldValidate (Coupler_FieldCreate, rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '<== Exiting  Coupler_FieldCreate',     &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
!
      RETURN
      END FUNCTION Coupler_FieldCreate
!
      SUBROUTINE Coupler_FindUnmapped (srcField, dstField,              &
     &                                 srcLandMask, dstLandMask,        &
     &                                 srcMId, dstMId, rc)
!
!=======================================================================
!                                                                      !
!  Modifies the grid mask associated with the destination field to     !
!  split the masked and unmasked grid cells. It is used during         !
!  regridding with extrapolation support.                              !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer (i4b), intent(in) :: srcLandMask, dstLandMask
      integer, intent(in) :: srcMId, dstMId
!
      integer, intent(out) :: rc
!
      TYPE (ESMF_Field), intent(in) :: srcField, dstField
!
!  Local variable declarations.
!
      integer :: i, j, k, srcTermProcessing
      integer :: localDE, localDEcount
      integer :: cLbnd(2), cUbnd(2)
!
      integer (i4b), pointer :: msk2d(:,:) => NULL()
!
      real (dp) :: IniVal
!
      real (dp), pointer :: bdy2d(:,:) => NULL()
      real (dp), pointer :: ptr2d(:,:) => NULL()
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", Coupler_FindUnmapped"

      character (ESMF_MAXSTR) :: Fname
!
      TYPE (ESMF_Grid)                :: grid
      TYPE (ESMF_Field)               :: aField, bField, cField
      TYPE (ESMF_UnmappedAction_Flag) :: unmap
      TYPE (ESMF_RegridMethod_Flag)   :: regridMethod
      TYPE (ESMF_RouteHandle)         :: routeHandle
      TYPE (ESMF_StaggerLoc)          :: sLoc
!
!-----------------------------------------------------------------------
!  Initialize return code flag to success state (no error).
!-----------------------------------------------------------------------
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '==> Entering Coupler_FindUnmapped',    &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
      rc=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!  Create dummy fields.
!-----------------------------------------------------------------------
!
      Fname='const_1'
      IniVal=1.0_dp
      aField=Coupler_FieldCreate(srcField, Fname, IniVal,               &
     &                           srcLandMask, rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
      Fname='const_2'
      IniVal=MISSING_dp
      bField=Coupler_FieldCreate(dstField, Fname, IniVal,               &
     &                           dstLandMask, rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
      Fname='const_3'
      IniVal=0.0_dp
      cField=Coupler_FieldCreate(dstField, Fname, IniVal,               &
     &                           -1_i4b, rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!-----------------------------------------------------------------------
!  Create 1st RouteHandle, which is used to find the boundary of the
!  destination grid.
!-----------------------------------------------------------------------
!
      unmap=ESMF_UNMAPPEDACTION_IGNORE
      IF (Iatmos.eq.srcMId) THEN                      ! HGA why Iatmos?
        regridMethod=ESMF_REGRIDMETHOD_NEAREST_STOD
      ELSE
        regridMethod=ESMF_REGRIDMETHOD_NEAREST_DTOS
      END IF
!
      srcTermProcessing=0
!
      CALL ESMF_FieldRegridStore (srcField=aField,                      &
     &                            dstField=bField,                      &
     &                            srcMaskValues=(/srcLandMask/),        &
     &                            dstMaskValues=(/dstLandMask/),        &
     &                            unmappedaction=unmap,                 &
     &                            routeHandle=routeHandle,              &
     &                            regridmethod=regridMethod,            &
     &                            srcTermProcessing=srcTermProcessing,  &
     &                            ignoreDegenerate=.TRUE.,              &
     &                            rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!-----------------------------------------------------------------------
!  Perform REGRID using 1st RouteHandle.
!-----------------------------------------------------------------------
!
      CALL ESMF_FieldRegrid (aField,                                    &
     &                       bField,                                    &
     &                       routeHandle,                               &
     &                       zeroregion=ESMF_REGION_EMPTY,              &
     &                       termorderflag=ESMF_TERMORDER_SRCSEQ,       &
     &                       rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!-----------------------------------------------------------------------
!  Create 2nd RouteHandle, which is used to find the unmapped grid
!  cells.
!-----------------------------------------------------------------------
!
      unmap=ESMF_UNMAPPEDACTION_IGNORE
      regridMethod=ESMF_REGRIDMETHOD_BILINEAR
!
      CALL ESMF_FieldRegridStore (srcField=aField,                      &
     &                            dstField=cField,                      &
     &                            srcMaskValues=(/srcLandMask/),        &
     &                            dstMaskValues=(/dstLandMask/),        &
     &                            unmappedaction=unmap,                 &
     &                            routeHandle=routeHandle,              &
     &                            regridmethod=regridMethod,            &
     &                            srcTermProcessing=srcTermProcessing,  &
     &                            ignoreDegenerate=.TRUE.,              &
     &                            rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!-----------------------------------------------------------------------
!  Perform REGRID using 2nd RouteHandle.
!-----------------------------------------------------------------------
!
      CALL ESMF_FieldRegrid (aField,                                    &
     &                       cField,                                    &
     &                       routeHandle,                               &
     &                       zeroregion=ESMF_REGION_TOTAL,              &
     &                       termorderflag=ESMF_TERMORDER_SRCSEQ,       &
     &                       rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!-----------------------------------------------------------------------
!  Query regridded field.
!-----------------------------------------------------------------------
!
      CALL ESMF_FieldGet (cField,                                       &
     &                    grid=grid,                                    &
     &                    staggerloc=sLoc,                              &
     &                    rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!  Get number of local decomposition elements (DEs) in the grid.
!  Usually, a single DE is associated with each Persistent Execution
!  Thread (PETs). Thus, localDEcount=1.
!
      CALL ESMF_GridGet (grid,                                          &
     &                   localDECount=localDEcount,                     &
     &                   rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!  Get pointer from fields.
!
      DE_LOOP : DO localDE=0,localDEcount-1
        CALL ESMF_FieldGet (bField,                                     &
     &                      localDe=localDE,                            &
     &                      farrayPtr=bdy2d,                            &
     &                      rc=rc)
        IF (ESMF_LogFoundError(rcToCheck=rc,                            &
     &                         msg=ESMF_LOGERR_PASSTHRU,                &
     &                         line=__LINE__,                           &
     &                         file=MyFile)) THEN
          RETURN
        END IF
!
        CALL ESMF_FieldGet (cField,                                     &
     &                      localDe=LocalDE,                            &
     &                      farrayPtr=ptr2d,                            &
     &                      computationalLBound=cLbnd,                  &
     &                      computationalUBound=cUbnd,                  &
     &                      rc=rc)
        IF (ESMF_LogFoundError(rcToCheck=rc,                            &
     &                         msg=ESMF_LOGERR_PASSTHRU,                &
     &                         line=__LINE__,                           &
     &                         file=MyFile)) THEN
          RETURN
        END IF
!
!  Get mask pointer from grid.
!
        CALL ESMF_GridGetItem (grid,                                    &
     &                         ESMF_GRIDITEM_MASK,                      &
     &                         staggerloc=sLoc,                         &
     &                         localDe=LocalDE,                         &
     &                         farrayPtr=msk2d,                         &
     &                         rc=rc)
        IF (ESMF_LogFoundError(rcToCheck=rc,                            &
     &                         msg=ESMF_LOGERR_PASSTHRU,                &
     &                         line=__LINE__,                           &
     &                         file=MyFile)) THEN
          RETURN
        END IF
!
!  Modify masking to split mapped and unmapped grid cells.
!
        DO j=cLbnd(2),cUbnd(2)
          DO i=cLbnd(1),cUbnd(1)
            IF ((bdy2d(i,j).lt.TOL_dp).and.                             &
     &          (msk2d(i,j).ne.dstLandMask)) THEN
              IF (ptr2d(i,j).lt.0.5_dp) THEN
                msk2d(i,j)=UNMAPPED_MASK
              ELSE
                msk2d(i,j)=MAPPED_MASK
              END IF
            END IF
          END DO
        END DO
!
!  Nullify pointer to make sure that it does not point on a random
!  part in the memory.
!
        IF (associated(ptr2d)) THEN
          nullify (ptr2d)
        END IF
        IF (associated(bdy2d)) THEN
          nullify (bdy2d)
        END IF
        IF (associated(msk2d)) THEN
          nullify (msk2d)
        END IF
      END DO DE_LOOP
!
!-----------------------------------------------------------------------
!  Remove temporary fields.
!-----------------------------------------------------------------------
!
      CALL ESMF_FieldDestroy (aField, rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
      CALL ESMF_FieldDestroy (bField, rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
      CALL ESMF_FieldDestroy (cField, rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '<== Exiting  Coupler_FindUnmapped',    &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
!
      RETURN
      END SUBROUTINE Coupler_FindUnmapped
!
#endif
      END MODULE esmf_coupler_mod
