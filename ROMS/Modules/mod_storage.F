#include "cppdefs.h"
      MODULE mod_storage

#if defined PROPAGATOR
!
!git $Id$
!svn $Id: mod_storage.F 1151 2023-02-09 03:08:53Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2023 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  ROMS/TOMS Generalized Stability Theory (GST) Analysis: ARPACK       !
!                                                                      !
!  LworkL       Size of Arnoldi iterations work array SworkL.          !
!  Lrvec        ARPACK logical to compute converged Ritz values.       !
!  NCV          Number of Lanczos vectors to compute.                  !
!  NEV          Number of eigenvalues to compute.                      !
!  Bvec         Lanczos/Arnoldi basis vectors.                         !
!  RvalueR      Real Ritz eigenvalues.                                 !
!  RvalueI      Imaginary Ritz eigenvalues.                            !
!  Rvector      Real Ritz eigenvectors.                                !
!  Swork        FULL state work array used in distributed memory       !
!                 communications.                                      !
!  SworkR       FULL state work array used to compute the ROMS state   !
!                 eigenvectors.                                        !
!  SworkD       State work array for ARPACK reverse communications.    !
!  SworkEV      ARPACK work array.                                     !
!  SworkL       ARPACK work array.                                     !
!  bmat         ARPACK eigenvalue value problem identifier:            !
!                 bmat='I'     standard eigenvalue problem             !
!                 bmat='G'     generalized eigenvalue problem          !
!  howmany      ARPACK form of basis functions identifier:             !
!                 howmany='A'  compute NEV Ritz vectors                !
!                 howmany='P'  compute NEV Schur vectors               !
!                 howmany='S'  compute some Ritz vectors using select  !
!  ido          ARPACK reverse communications flag (input/output).     !
!  iparam       ARPACK input/output integer parameters.                !
!  ipntr        ARPACK pointer to mark the starting location in SworkD !
!                 and SworkK arrays used in the Arnoldi iterations.    !
!  info         ARPACK Information (input) and error flag (output).    !
!  norm         Euclidean norm.                                        !
!  pick         ARPACK logical switch of Ritz vectors to compute.      !
!  resid        Initial/final residual vector.                         !
!  sigmaR       ARPACK real part of the shifts (not referenced).       !
!  sigmaI       ARPACK imaginary part of the shifts (not referenced).  !
# ifdef SO_SEMI
!  so_state     Stochastic optimals adjoint state surface forcing      !
!                 vector sample in time.                               !
# endif
!  which        ARPACK Ritz eigenvalues to compute identifier:         !
!                 which='LA'   compute NEV largest (algebraic)         !
!                 which='SA'   compute NEV smallest (algebraic)        !
!                 which='LM'   compute NEV largest in magnitude        !
!                 which='SM'   compute NEV smallest in magnitude       !
!                 which='BE'   compute NEV from each end of spectrum   !
!                                                                      !
!=======================================================================
!
      USE mod_param
!
      implicit none
!
      PUBLIC :: allocate_storage
      PUBLIC :: deallocate_storage
!
!-----------------------------------------------------------------------
!  Define T_STORAGE structure and other module variables.
!-----------------------------------------------------------------------
!
!  State work arrays structure for nested grids.
!
      TYPE T_STORAGE

        real(r8), pointer :: Bvec(:,:)              ! [Nstr:Nend,NCV]
        real(r8), pointer :: Rvector(:,:)           ! [Nstr:Nend,NEV+1]
        real(r8), pointer :: SworkD(:)              ! [3*Nstate]
        real(r8), pointer :: resid(:)               ! [Nstr:Nend]
# ifdef STOCHASTIC_OPT
        real(r8), pointer :: ad_Work(:)             ! [Mstate]
# endif
# if defined STOCHASTIC_OPT && defined HESSIAN_SO
        real(r8), pointer :: my_state(:)            ! [Nstr:Nend]
# endif
# ifdef SO_SEMI
        real(r8), pointer :: so_state(:,:)          ! [Nstr:Nend,Nsemi]
# endif

      END TYPE T_STORAGE
!
      TYPE (T_STORAGE), allocatable :: STORAGE(:)
!
!  ARPACK logical switches.
!
      logical :: Lrvec                    ! Ritz values switch
      logical, allocatable :: pick(:,:)   ! Ritz vectors, [NCV,Ngrids]
!
!  Eigenproblem parameters.
!
      integer :: NCV                      ! number of Lanczos vectors
      integer :: NEV                      ! number of eigenvalues
      integer :: LworkL                   ! size of array SworkL

      integer, allocatable :: LworkD(:)   ! size of array SworkD
!
!  ARPACK integer parameters and flags.
!
      integer, allocatable :: iparam(:,:) ! eigenproblem parameters
      integer, allocatable :: ipntr(:,:)  ! index for work arrays
      integer, allocatable :: ido(:)      ! reverse communication flag
      integer, allocatable :: info(:)     ! information and error flag
!
!  ARPACK floating-point parameters.
!
      real(r8) :: sigmaI                  ! real part shifts
      real(r8) :: sigmaR                  ! imaginary part shifts
!
!  ARPACK character parameters.
!
      character (len=1) :: bmat           ! eigenvalue problem ID
      character (len=1) :: howmany        ! basis functions ID
      character (len=2) :: which          ! Ritz values ID
!
!  ARPACK work arrays.
!
      real(r8), allocatable :: RvalueR(:,:)       ! [NEV+1,Ngrids]
      real(r8), allocatable :: RvalueI(:,:)       ! [NEV+1,Ngrids]
      real(r8), allocatable :: norm(:,:)          ! [NEV+1,Ngrids]
      real(r8), allocatable :: SworkL(:,:)        ! [LworkL,Ngrids]
      real(r8), allocatable :: SworkEV(:,:)       ! [3*NCV,Ngrids]
# ifdef DISTRIBUTE
      real(r8), allocatable :: Swork(:)           ! [Mstate]
# endif
      real(r8), pointer :: SworkR(:)              ! [Mstate]
!
!  ARPACK private common blocks containing parameters needed for
!  checkpoiniting. The original include files "i_aupd.h" and
!  "idaup2.h" have several parameters in their commom blocks. All
!  these values are compacted here in vector arrays to allow IO
!  manipulations during checkpointing.
!
      integer  :: iaitr(8), iaup2(8), iaupd(20)
      logical  :: laitr(5), laup2(5)
      real(r8) :: raitr(8), raup2(2)
!
      common /i_aupd/ iaupd
# ifdef DOUBLE_PRECISION
      common /idaitr/ iaitr
      common /ldaitr/ laitr
      common /rdaitr/ raitr
      common /idaup2/ iaup2
      common /ldaup2/ laup2
      common /rdaup2/ raup2
# else
      common /isaitr/ iaitr
      common /lsaitr/ laitr
      common /rsaitr/ raitr
      common /isaup2/ iaup2
      common /lsaup2/ laup2
      common /rsaup2/ raup2
# endif
!
!  ARPACK debugging common block.
!
      integer :: logfil, ndigit, mgetv0
      integer :: msaupd, msaup2, msaitr, mseigt, msapps, msgets, mseupd
      integer :: mnaupd, mnaup2, mnaitr, mneigh, mnapps, mngets, mneupd
      integer :: mcaupd, mcaup2, mcaitr, mceigh, mcapps, mcgets, mceupd
!
      common /debug/                                                    &
     &       logfil, ndigit, mgetv0,                                    &
     &       msaupd, msaup2, msaitr, mseigt, msapps, msgets, mseupd,    &
     &       mnaupd, mnaup2, mnaitr, mneigh, mnapps, mngets, mneupd,    &
     &       mcaupd, mcaup2, mcaitr, mceigh, mcapps, mcgets, mceupd
!
      CONTAINS
!
      SUBROUTINE allocate_storage
!
!=======================================================================
!                                                                      !
!  This routine allocates and initialize module variables. For now,    !
!  only non-nested applications are considered.                        !
!                                                                      !
!=======================================================================
!
      USE mod_scalars
!
!  Local variable declarations
!
      integer :: i, j, ng

      integer, allocatable :: Mstr(:), Mend(:)
!
      real(r8), parameter :: IniVal = 0.0_r8
!
!-----------------------------------------------------------------------
!  Allocate module variables.
!-----------------------------------------------------------------------
!
!  Allocate local dimension parameters.
!
      allocate ( Mstr(Ngrids) )
      allocate ( Mend(Ngrids) )
!
!  Allocate parameters associated with nested grids.
!
      allocate ( ido(Ngrids) )
      allocate ( info(Ngrids) )
      allocate ( LworkD(Ngrids) )

      allocate ( iparam(11,Ngrids) )
# if defined AFT_EIGENMODES || defined FT_EIGENMODES
      allocate ( ipntr(14,Ngrids) )
# else
      allocate ( ipntr(11,Ngrids) )
# endif
!
!  Determine size of work array SworkL:
!
# if defined OPT_PERTURBATION
      LworkL=NCV*(NCV+8)
# elif defined HESSIAN_SV
      LworkL=NCV*(NCV+8)
# elif defined HESSIAN_FSV
      LworkL=NCV*(NCV+8)
# elif defined FORCING_SV
      LworkL=NCV*(NCV+8)
# elif defined STOCHASTIC_OPT
      LworkL=NCV*(NCV+8)
# elif defined SO_SEMI
      LworkL=NCV*(NCV+8)
# elif defined FT_EIGENMODES || defined AFT_EIGENMODES
      LworkL=3*NCV*NCV+6*NCV
# endif
# ifdef SO_SEMI
      DO ng=1,Ngrids
        Nsemi(ng)=1+ntimes(ng)/nADJ(ng)
      END DO
# endif
!
!  Set local dimensions of storage arrays according to the propagator
!  driver.
!
      DO ng=1,Ngrids
# if defined HESSIAN_FSV || defined HESSIAN_SO || defined HESSIAN_SV
        Mstr(ng)=1
        Mend(ng)=Ninner
        LworkD(ng)=3*Ninner
# else
        Mstr(ng)=Nstr(ng)
        Mend(ng)=Nend(ng)
        LworkD(ng)=3*Nstate(ng)
# endif
      END DO
!
!  Allocate structure.
!
      allocate ( STORAGE(Ngrids) )

      DO ng=1,Ngrids
        allocate ( STORAGE(ng) % Bvec(Mstr(ng):Mend(ng),NCV) )
        Dmem(ng)=Dmem(ng)+REAL((Mend(ng)-Mstr(ng))*NCV,r8)

        allocate ( STORAGE(ng) % Rvector(Mstr(ng):Mend(ng),NEV+1) )
        Dmem(ng)=Dmem(ng)+REAL((Mend(ng)-Mstr(ng))*(NEV+1),r8)

        allocate ( STORAGE(ng) % SworkD(LworkD(ng)) )
        Dmem(ng)=Dmem(ng)+REAL(LworkD(ng),r8)

        allocate ( STORAGE(ng) % resid(Mstr(ng):Mend(ng)) )
        Dmem(ng)=Dmem(ng)+REAL(Mend(ng)-Mstr(ng),r8)

# ifdef STOCHASTIC_OPT
        allocate ( STORAGE(ng) % ad_Work(MAXVAL(Mstate)) )
        Dmem(ng)=Dmem(ng)+REAL(MAXVAL(Mstate),r8)
# endif

# if defined STOCHASTIC_OPT && defined HESSIAN_SO
        allocate ( STORAGE(ng) % my_state(Nstr(ng):Nend(ng)) )
        Dmem(ng)=Dmem(ng)+REAL(Nend(ng)-Nstr(ng),r8)
# endif

# ifdef SO_SEMI
        allocate ( STORAGE(ng) % so_state(Mstr(ng):Mend(ng),Nsemi(ng)) )
        Dmem(ng)=Dmem(ng)+REAL((Mend(ng)-Mstr(ng))*Nsemi(ng),r8)
# endif
      END DO
!
!  Allocate other arrays.
!
      ng=1
!
      allocate ( pick(NCV,Ngrids) )
      Dmem(ng)=Dmem(ng)+REAL(NCV*Ngrids,r8)

      allocate ( norm(NEV+1,Ngrids) )
      Dmem(ng)=Dmem(ng)+REAL((NEV+1)*Ngrids,r8)

      allocate ( RvalueR(NEV+1,Ngrids) )
      Dmem(ng)=Dmem(ng)+REAL((NEV+1)*Ngrids,r8)

      allocate ( RvalueI(NEV+1,Ngrids) )
      Dmem(ng)=Dmem(ng)+REAL((NEV+1)*Ngrids,r8)

      allocate ( SworkEV(3*NCV,Ngrids) )
      Dmem(ng)=Dmem(ng)+3.0_r8*REAL(NCV,r8)

      allocate ( SworkL(LworkL,Ngrids) )
      Dmem(ng)=Dmem(ng)+REAL(LworkL*Ngrids,r8)

      allocate ( SworkR(MAXVAL(Mstate)) )
      Dmem(ng)=Dmem(ng)+REAL(MAXVAL(Mstate),r8)

# ifdef DISTRIBUTE
      allocate ( Swork(MAXVAL(Mstate)) )
      Dmem(ng)=Dmem(ng)+REAL(MAXVAL(Mstate),r8)
# endif
!
!-----------------------------------------------------------------------
!  Initialize module variables.
!-----------------------------------------------------------------------
!
      iparam=0
      ipntr=0

      DO ng=1,Ngrids
        DO j=1,NCV
          pick(j,ng) = .TRUE.
          DO i=Mstr(ng),Mend(ng)
            STORAGE(ng) % Bvec(i,j) = IniVal
          END DO
        END DO
        DO j=1,NEV+1
          norm(j,ng) = IniVal
          RvalueR(j,ng) = IniVal
          RvalueI(j,ng) = IniVal
          DO i=Mstr(ng),Mend(ng)
            STORAGE(ng) % Rvector(i,j) = IniVal
          END DO
        END DO
        DO i=Mstr(ng),Mend(ng)
          STORAGE(ng) % resid(i) = IniVal
        END DO
# ifdef STOCHASTIC_OPT
        DO i=1,MAXVAL(Mstate)
          STORAGE(ng) % ad_Work(i) = IniVal
        END DO
# endif
# if defined STOCHASTIC_OPT && defined HESSIAN_SO
        DO i=Nstr(ng),Nend(ng)
          STORAGE(ng) % my_state(i) = IniVal
        END DO
# endif
# ifdef SO_SEMI
        DO j=1,Nsemi(ng)
          DO i=Mstr(ng),Mend(ng)
            STORAGE(ng) % so_state(i,j) = IniVal
          END DO
        END DO
# endif
      END DO
      DO i=1,MAXVAL(Mstate)
        SworkR(i) = IniVal
      END DO
# ifdef DISTRIBUTE
      DO i=1,MAXVAL(Mstate)
        Swork(i) = IniVal
      END DO
# endif
!
      DO ng=1,Ngrids
        DO i=1,LworkD(ng)
          STORAGE(ng) % SworkD(i) = IniVal
        END DO
        DO i=1,3*NCV
          SworkEV(i,ng) = IniVal
        END DO
        DO i=1,LworkL
          SworkL(i,ng) = IniVal
        END DO
      END DO
!
      RETURN
      END SUBROUTINE allocate_storage
!
      SUBROUTINE deallocate_storage (ng)
!
!=======================================================================
!                                                                      !
!  This routine deallocates and initialize module variables.           !
!                                                                      !
!=======================================================================
!
      USE mod_param, ONLY : Ngrids
!
!  Imported variable declarations
!
      integer, intent(in) :: ng

# ifdef SUBOBJECT_DEALLOCATION
!
!-----------------------------------------------------------------------
!  Deallocate each variable in the derived-type T_STORAGE structure
!  separately.
!-----------------------------------------------------------------------
!
      IF (associated(STORAGE(ng)%Bvec)) THEN
        deallocate ( STORAGE(ng)%Bvec )
      END IF

      IF (associated(STORAGE(ng)%Rvector)) THEN
        deallocate ( STORAGE(ng)%Rvector )
      END IF

      IF (associated(STORAGE(ng)%SworkD)) THEN
        deallocate ( STORAGE(ng)%SworkD )
      END IF

      IF (associated(STORAGE(ng)%resid)) THEN
        deallocate ( STORAGE(ng)%resid )
      END IF

#  ifdef STOCHASTIC_OPT
      IF (associated(STORAGE(ng)%ad_Work)) THEN
        deallocate ( STORAGE(ng)%ad_Work )
      END IF
#  endif

#  if defined STOCHASTIC_OPT && defined HESSIAN_SO
      IF (associated(STORAGE(ng)%my_state)) THEN
        deallocate ( STORAGE(ng)%my_state )
      END IF
#  endif

#  ifdef SO_SEMI
      IF (associated(STORAGE(ng)%so_state)) THEN
        deallocate ( STORAGE(ng)%so_state )
      END IF
#  endif
# endif
!
!-----------------------------------------------------------------------
!  Deallocate derived-type STORAGE structure.
!-----------------------------------------------------------------------
!
      IF (ng.eq.Ngrids) THEN
        IF (allocated(STORAGE)) deallocate ( STORAGE )
      END IF
!
!-----------------------------------------------------------------------
!  Deallocate other variables in module.
!-----------------------------------------------------------------------
!
      IF (allocated(Mstr))    deallocate ( Mstr )
      IF (allocated(Mend))    deallocate ( Mend )
      IF (allocated(ido))     deallocate ( ido )
      IF (allocated(info))    deallocate ( info )
      IF (allocated(LworkD))  deallocate ( LworkD )
      IF (allocated(iparam))  deallocate ( iparam )
      IF (allocated(ipntr))   deallocate ( ipntr )
      IF (allocated(pick))    deallocate ( pick )
      IF (allocated(norm))    deallocate ( norm )
      IF (allocated(RvalueR)) deallocate ( RvalueR )
      IF (allocated(RvalueI)) deallocate ( RvalueI )
      IF (allocated(SworkEV)) deallocate ( SworkEV )
      IF (allocated(SworkL))  deallocate ( SworkL )
      IF (allocated(SworkR))  deallocate ( SworkR )
# ifdef DISTRIBUTE
      IF (allocated(Swork))   deallocate ( Swork )
# endif
!
      RETURN
      END SUBROUTINE deallocate_storage
#endif
      END MODULE mod_storage
