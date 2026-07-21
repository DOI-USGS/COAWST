module NoahmpIO_fi

  use NoahmpIOVarType, ONLY: NoahmpIO_type
  use NoahmpIOVarInitMod, ONLY: NoahmpIOVarInitDefault
  use NoahmpInitMainMod, ONLY: NoahmpInitMain
  use NoahmpReadNamelistMod, ONLY: NoahmpReadNamelist
  use NoahmpReadTableMod, ONLY: NoahmpReadTable
  use NoahmpReadLandMod, ONLY: NoahmpReadLandHeader, NoahmpReadLandMain
  use NoahmpDriverMainMod, ONLY: NoahmpDriverMain
  use NoahmpWriteLandMod, ONLY: NoahmpWriteLand

  use  iso_c_binding

  implicit none

  ! ---------------------------------------------------------------------------
  ! Public type for NoahmpIO_level
  ! --------------------------------------------------------------------------- 
  type, public :: NoahmpIO_level
    type(NoahmpIO_type), allocatable :: NoahmpIO(:)
  end type NoahmpIO_level

  ! ---------------------------------------------------------------------------
  ! Public declaration for NoahmpIO array per level, hardcoding levels 0-2 
  ! ---------------------------------------------------------------------------
  type(NoahmpIO_level), save, target, public :: NoahmpIO_vect(0:2)

  ! ---------------------------------------------------------------------------
  ! Mirror of extern C struct
  !
  ! DEVNOTE :: The order variables between C struct and the binded Fortran type
  !            should be consistent for memory managment. 
  ! ---------------------------------------------------------------------------
  type, bind(c), public :: NoahmpIO_type_fi
    type(C_PTR)                                         ::  ids,ide, &          ! d -> domain
                                                            jds,jde, &          ! d -> domain
                                                            kds,kde, &          ! d -> domain
                                                            ims,ime, &          ! m -> memory
                                                            jms,jme, &          ! m -> memory
                                                            kms,kme, &          ! m -> memory
                                                            its,ite, &          ! t -> tile
                                                            jts,jte, &          ! t -> tile
                                                            kts,kte             ! t -> tile

    type(C_PTR) :: xstart, xend, ystart, yend
    type(C_PTR) :: nsoil, nsnow
    type(C_PTR) :: itimestep, ntime
    type(C_PTR) :: rank, blkid, level
    type(C_PTR) :: comm
    type(C_PTR) :: XLAT, WSLAKEXY
    type(C_PTR) :: U_PHY, T_PHY, V_PHY, QV_CURR
    type(C_PTR) :: HFX, LH
    type(C_PTR) :: SWDOWN, GLW, TSK, EMISS
    type(C_PTR) :: ALBSFCDIRXY, ALBSFCDIFXY
    type(C_PTR) :: COSZEN, P8W
    type(C_PTR) :: TAU_EW, TAU_NS
  end type NoahmpIO_type_fi

contains

  subroutine NoahmpIOScalarInitDefault_fi(NoahmpIO_cptr) bind(C, name="NoahmpIOScalarInitDefault_fi")
    use  iso_c_binding, only : C_INT
    implicit none
    type(NoahmpIO_type_fi), intent(inout) :: NoahmpIO_cptr
    integer(C_INT), pointer :: level, bid
   
    call C_F_POINTER(NoahmpIO_cptr%BLKID, bid) 
    call C_F_POINTER(NoahmpIO_cptr%LEVEL, level)

    call C_F_POINTER(NoahmpIO_cptr%BLKID, NoahmpIO_vect(level)%NoahmpIO(bid)%blkid) 
    call C_F_POINTER(NoahmpIO_cptr%LEVEL, NoahmpIO_vect(level)%NoahmpIO(bid)%level)

    call C_F_POINTER(NoahmpIO_cptr%XSTART, NoahmpIO_vect(level)%NoahmpIO(bid)%XSTART)
    call C_F_POINTER(NoahmpIO_cptr%XEND, NoahmpIO_vect(level)%NoahmpIO(bid)%XEND)
    call C_F_POINTER(NoahmpIO_cptr%YSTART, NoahmpIO_vect(level)%NoahmpIO(bid)%YSTART)
    call C_F_POINTER(NoahmpIO_cptr%YEND, NoahmpIO_vect(level)%NoahmpIO(bid)%YEND)

    call C_F_POINTER(NoahmpIO_cptr%NSOIL, NoahmpIO_vect(level)%NoahmpIO(bid)%NSOIL)
    call C_F_POINTER(NoahmpIO_cptr%NSNOW, NoahmpIO_vect(level)%NoahmpIO(bid)%NSNOW)

    call C_F_POINTER(NoahmpIO_cptr%IDS, NoahmpIO_vect(level)%NoahmpIO(bid)%IDS)
    call C_F_POINTER(NoahmpIO_cptr%IDE, NoahmpIO_vect(level)%NoahmpIO(bid)%IDE)
    call C_F_POINTER(NoahmpIO_cptr%JDS, NoahmpIO_vect(level)%NoahmpIO(bid)%JDS)
    call C_F_POINTER(NoahmpIO_cptr%JDE, NoahmpIO_vect(level)%NoahmpIO(bid)%JDE)
    call C_F_POINTER(NoahmpIO_cptr%KDS, NoahmpIO_vect(level)%NoahmpIO(bid)%KDS)
    call C_F_POINTER(NoahmpIO_cptr%KDE, NoahmpIO_vect(level)%NoahmpIO(bid)%KDE)

    call C_F_POINTER(NoahmpIO_cptr%IMS, NoahmpIO_vect(level)%NoahmpIO(bid)%IMS)
    call C_F_POINTER(NoahmpIO_cptr%IME, NoahmpIO_vect(level)%NoahmpIO(bid)%IME)
    call C_F_POINTER(NoahmpIO_cptr%JMS, NoahmpIO_vect(level)%NoahmpIO(bid)%JMS)
    call C_F_POINTER(NoahmpIO_cptr%JME, NoahmpIO_vect(level)%NoahmpIO(bid)%JME)
    call C_F_POINTER(NoahmpIO_cptr%KMS, NoahmpIO_vect(level)%NoahmpIO(bid)%KMS)
    call C_F_POINTER(NoahmpIO_cptr%KME, NoahmpIO_vect(level)%NoahmpIO(bid)%KME)

    call C_F_POINTER(NoahmpIO_cptr%ITS, NoahmpIO_vect(level)%NoahmpIO(bid)%ITS)
    call C_F_POINTER(NoahmpIO_cptr%ITE, NoahmpIO_vect(level)%NoahmpIO(bid)%ITE)
    call C_F_POINTER(NoahmpIO_cptr%JTS, NoahmpIO_vect(level)%NoahmpIO(bid)%JTS)
    call C_F_POINTER(NoahmpIO_cptr%JTE, NoahmpIO_vect(level)%NoahmpIO(bid)%JTE)
    call C_F_POINTER(NoahmpIO_cptr%KTS, NoahmpIO_vect(level)%NoahmpIO(bid)%KTS)
    call C_F_POINTER(NoahmpIO_cptr%KTE, NoahmpIO_vect(level)%NoahmpIO(bid)%KTE)

    call C_F_POINTER(NoahmpIO_cptr%ITIMESTEP, NoahmpIO_vect(level)%NoahmpIO(bid)%ITIMESTEP)
    call C_F_POINTER(NoahmpIO_cptr%NTIME, NoahmpIO_vect(level)%NoahmpIO(bid)%NTIME)

    call C_F_POINTER(NoahmpIO_cptr%RANK, NoahmpIO_vect(level)%NoahmpIO(bid)%RANK)
    call C_F_POINTER(NoahmpIO_cptr%COMM, NoahmpIO_vect(level)%NoahmpIO(bid)%COMM)

  end subroutine NoahmpIOScalarInitDefault_fi

  subroutine NoahmpIOVarInitDefault_fi(NoahmpIO_cptr) bind(C, name="NoahmpIOVarInitDefault_fi")
    use  iso_c_binding, only : C_INT
    implicit none
    type(NoahmpIO_type_fi), intent(inout) :: NoahmpIO_cptr
    integer(C_INT), pointer :: level, bid    

    call C_F_POINTER(NoahmpIO_cptr%BLKID, bid)
    call C_F_POINTER(NoahmpIO_cptr%LEVEL, level)

    call NoahmpIOVarInitDefault(NoahmpIO_vect(level)%NoahmpIO(bid))

    NoahmpIO_cptr%XLAT = C_LOC(NoahmpIO_vect(level)%NoahmpIO(bid)%XLAT)
    NoahmpIO_cptr%WSLAKEXY = C_LOC(NoahmpIO_vect(level)%NoahmpIO(bid)%WSLAKEXY)
    NoahmpIO_cptr%T_PHY = C_LOC(NoahmpIO_vect(level)%NoahmpIO(bid)%T_PHY)
    NoahmpIO_cptr%U_PHY = C_LOC(NoahmpIO_vect(level)%NoahmpIO(bid)%U_PHY)
    NoahmpIO_cptr%V_PHY = C_LOC(NoahmpIO_vect(level)%NoahmpIO(bid)%V_PHY)
    NoahmpIO_cptr%QV_CURR = C_LOC(NoahmpIO_vect(level)%NoahmpIO(bid)%QV_CURR)
    NoahmpIO_cptr%HFX = C_LOC(NoahmpIO_vect(level)%NoahmpIO(bid)%HFX)
    NoahmpIO_cptr%LH = C_LOC(NoahmpIO_vect(level)%NoahmpIO(bid)%LH)
    NoahmpIO_cptr%SWDOWN = C_LOC(NoahmpIO_vect(level)%NoahmpIO(bid)%SWDOWN)
    NoahmpIO_cptr%GLW = C_LOC(NoahmpIO_vect(level)%NoahmpIO(bid)%GLW)
    NoahmpIO_cptr%TSK = C_LOC(NoahmpIO_vect(level)%NoahmpIO(bid)%TSK)
    NoahmpIO_cptr%EMISS = C_LOC(NoahmpIO_vect(level)%NoahmpIO(bid)%EMISS)
    NoahmpIO_cptr%ALBSFCDIRXY = C_LOC(NoahmpIO_vect(level)%NoahmpIO(bid)%ALBSFCDIRXY)
    NoahmpIO_cptr%ALBSFCDIFXY = C_LOC(NoahmpIO_vect(level)%NoahmpIO(bid)%ALBSFCDIFXY)
    NoahmpIO_cptr%COSZEN = C_LOC(NoahmpIO_vect(level)%NoahmpIO(bid)%COSZEN)
    NoahmpIO_cptr%P8W = C_LOC(NoahmpIO_vect(level)%NoahmpIO(bid)%P8W)
    NoahmpIO_cptr%TAU_EW = C_LOC(NoahmpIO_vect(level)%NoahmpIO(bid)%TAU_EW)
    NoahmpIO_cptr%TAU_NS = C_LOC(NoahmpIO_vect(level)%NoahmpIO(bid)%TAU_NS)
  end subroutine NoahmpIOVarInitDefault_fi

  subroutine NoahmpInitMain_fi(NoahmpIO_cptr) bind(C, name="NoahmpInitMain_fi")
    use  iso_c_binding, only : C_INT
    implicit none 
    type(NoahmpIO_type_fi), intent(inout) :: NoahmpIO_cptr
    integer(C_INT), pointer :: level,bid
    call C_F_POINTER(NoahmpIO_cptr%BLKID, bid)
    call C_F_POINTER(NoahmpIO_cptr%LEVEL, level)
    call NoahmpInitMain(NoahmpIO_vect(level)%NoahmpIO(bid))
  end subroutine NoahmpInitMain_fi

  subroutine NoahmpReadTable_fi(NoahmpIO_cptr) bind(C, name="NoahmpReadTable_fi")
    use  iso_c_binding, only : C_INT
    implicit none 
    type(NoahmpIO_type_fi), intent(inout) :: NoahmpIO_cptr
    integer(C_INT), pointer :: level, bid
    call C_F_POINTER(NoahmpIO_cptr%BLKID, bid)
    call C_F_POINTER(NoahmpIO_cptr%LEVEL, level)
    call NoahmpReadTable(NoahmpIO_vect(level)%NoahmpIO(bid))
  end subroutine NoahmpReadTable_fi

  subroutine NoahmpReadNamelist_fi(NoahmpIO_cptr) bind(C, name="NoahmpReadNamelist_fi")
    use  iso_c_binding, only : C_INT 
    implicit none 
    type(NoahmpIO_type_fi), intent(inout) :: NoahmpIO_cptr
    integer(C_INT), pointer :: level, bid
    call C_F_POINTER(NoahmpIO_cptr%BLKID, bid)
    call C_F_POINTER(NoahmpIO_cptr%LEVEL, level)
    call NoahmpReadNamelist(NoahmpIO_vect(level)%NoahmpIO(bid))
  end subroutine NoahmpReadNamelist_fi

  subroutine NoahmpReadLandHeader_fi(NoahmpIO_cptr) bind(C, name="NoahmpReadLandHeader_fi")
    use  iso_c_binding, only : C_INT 
    implicit none 
    type(NoahmpIO_type_fi), intent(inout) :: NoahmpIO_cptr
    integer(C_INT), pointer :: level, bid
    call C_F_POINTER(NoahmpIO_cptr%BLKID, bid)
    call C_F_POINTER(NoahmpIO_cptr%LEVEL, level)
    call NoahmpReadLandHeader(NoahmpIO_vect(level)%NoahmpIO(bid))
  end subroutine NoahmpReadLandHeader_fi

  subroutine NoahmpReadLandMain_fi(NoahmpIO_cptr) bind(C, name="NoahmpReadLandMain_fi")
    use  iso_c_binding, only : C_INT 
    implicit none 
    type(NoahmpIO_type_fi), intent(inout) :: NoahmpIO_cptr
    integer(C_INT), pointer :: level, bid
    call C_F_POINTER(NoahmpIO_cptr%BLKID, bid)
    call C_F_POINTER(NoahmpIO_cptr%LEVEL, level)
    call NoahmpReadLandMain(NoahmpIO_vect(level)%NoahmpIO(bid))
  end subroutine NoahmpReadLandMain_fi

  subroutine NoahmpDriverMain_fi(NoahmpIO_cptr) bind(C, name="NoahmpDriverMain_fi")
    use iso_c_binding, only : C_INT 
    implicit none 
    type(NoahmpIO_type_fi), intent(inout) :: NoahmpIO_cptr
    integer(C_INT), pointer :: level, bid
    call C_F_POINTER(NoahmpIO_cptr%BLKID, bid)
    call C_F_POINTER(NoahmpIO_cptr%LEVEL, level)
    call NoahmpDriverMain(NoahmpIO_vect(level)%NoahmpIO(bid))
  end subroutine NoahmpDriverMain_fi

  subroutine NoahmpWriteLand_fi(NoahmpIO_cptr, filenum) bind(C, name="NoahmpWriteLand_fi")
    use iso_c_binding, only : C_INT 
    implicit none 
    type(NoahmpIO_type_fi), intent(inout) :: NoahmpIO_cptr
    integer(C_INT), intent(in) :: filenum
    integer(C_INT), pointer :: level, bid
    call C_F_POINTER(NoahmpIO_cptr%BLKID, bid)
    call C_F_POINTER(NoahmpIO_cptr%LEVEL, level)
    call NoahmpWriteLand(NoahmpIO_vect(level)%NoahmpIO(bid), filenum, SIZE(NoahmpIO_vect(level)%NoahmpIO))
  end subroutine NoahmpWriteLand_fi

  subroutine NoahmpIOTypeVectInit_fi(level, NBlocks) bind(C, name="NoahmpIOTypeVectInit_fi")
    use iso_c_binding, only : C_INT
    implicit none
    integer(C_INT), intent(in) :: level
    integer(C_INT), intent(in) :: NBlocks
    allocate(NoahmpIO_vect(level)%NoahmpIO(0:NBlocks-1))
  end subroutine NoahmpIOTypeVectInit_fi

end module NoahmpIO_fi
