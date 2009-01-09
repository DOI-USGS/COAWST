#include "cppdefs.h"
      MODULE mod_storage

#if defined PROPAGATOR
!
!svn $Id: mod_storage.F 1338 2008-01-24 02:47:35Z jcwarner $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  ROMS/TOMS Generalized Stability Theory (GST) Analysis: ARPACK       !
!                                                                      !
!  LworkL       Size of Arnoldi iterations work array SworkL.          !
!  Lrvec        ARPACK logical to compute converged Ritz values.       !
!  Mstate       Size of FULL state vector (water points only).         !
# ifdef SO_SEMI
!  Nsemi        Number of time record saved in seminorm strochastic    !
!                 optimals adjoint state vector.                       !
# endif
!  Nsize        Size of the eigenvalue problem: Nend-Nstr+1.           !
!  Nstate       Size of NODE partition state vector: Nstate=Mstate in  ! 
!                 serial applications.                                 !
!  Nstr         State vector node partition starting index.            !
!  Nend         State vector node partition ending   index.            !
!  NCV          Number of Lanczos vectors to compute.                  !
!  NEV          Number of eigenvalues to compute.                      !
!  Bvec         Lanczos/Arnoldi basis vectors.                         !
!  RvalueR      Real Ritz eigenvalues.                                 !
!  RvalueI      Imaginary Ritz eigenvalues.                            !
!  Rvector      Real Ritz eigenvectors.                                !
!  Swork        FULL state work array used in distributed memory       !
!                 communications.                                      !
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
!  ido          ARPACK reverse communications flag (input/outpu).      !
!  iparam       ARPACK input/output integer parameters.                !
!  ipntr        ARPACK pointer to mark the starting location in SworkD !
!                 and SworkK arrays used in the Arnoldi iterations.    !
!  info         ARPACK Information (input) and error flag (output).    !
!  norm         Euclidean norm.                                        !
!  resid        Initial/final residual vector.                         !
!  select       ARPACK logical switch of Ritz vectors to compute.      !
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
        integer, dimension(Ngrids) :: Mstate
        integer, dimension(Ngrids) :: Nstate
        integer, dimension(Ngrids) :: Nstr
        integer, dimension(Ngrids) :: Nend
!
        logical :: Lrvec
!
        integer :: NCV
        integer :: NEV
# ifdef SO_SEMI
        integer :: Nsemi
# endif
        integer :: Nsize
        integer :: LworkL
        integer :: ido
        integer :: info
!
        real(r8) :: sigmaI
        real(r8) :: sigmaR
!
        character (len=1) :: bmat, howmany
        character (len=2) :: which
!
        integer, dimension(11) :: iparam
# if defined AFT_EIGENMODES || defined FT_EIGENMODES
        integer, dimension(14) :: ipntr 
# else
        integer, dimension(11) :: ipntr 
# endif
!
        logical, allocatable :: select(:)           ! [NCV]
!
        real(r8), allocatable :: Bvec(:,:)          ! [Nstr:Nend,NCV]
        real(r8), allocatable :: RvalueR(:)         ! [NEV+1]
        real(r8), allocatable :: RvalueI(:)         ! [NEV+1]
        real(r8), allocatable :: Rvector(:,:)       ! [Nstr:Nend,NEV+1]
# ifdef DISTRIBUTE
        real(r8), allocatable :: Swork(:)           ! [Mstate]
# endif
        real(r8), allocatable :: SworkD(:)          ! [3*Nstate]
        real(r8), allocatable :: SworkEV(:)         ! [3*NCV]
        real(r8), allocatable :: SworkL(:)          ! [LworkL]
        real(r8), allocatable :: norm(:)            ! [NEV+1]
        real(r8), allocatable :: resid(:)           ! [Nstr:Nend]
# ifdef SO_SEMI
        real(r8), allocatable :: so_state(:,:)      ! [Nstr:Nend,Nsemi]
# endif
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

      common /debug/                                                    &
     &       logfil, ndigit, mgetv0,                                    &
     &       msaupd, msaup2, msaitr, mseigt, msapps, msgets, mseupd,    &
     &       mnaupd, mnaup2, mnaitr, mneigh, mnapps, mngets, mneupd,    &
     &       mcaupd, mcaup2, mcaitr, mceigh, mcapps, mcgets, mceupd

      CONTAINS

      SUBROUTINE allocate_storage (ng)
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
!  Imported variable declarations.
!
      integer, intent(in) :: ng
!
!  Local variable declarations
!
      integer :: i, j

      real(r8), parameter :: IniVal = 0.0_r8
!
!-----------------------------------------------------------------------
!  Allocate module variables.
!-----------------------------------------------------------------------
!
!  Determine size of work array SworkL:
!
# if defined OPT_PERTURBATION
      LworkL=NCV*(NCV+8)
# elif defined FORCING_SV
      LworkL=NCV*(NCV+8)
# elif defined SO_SEMI
      LworkL=NCV*(NCV+8)
# elif defined FT_EIGENMODES || defined AFT_EIGENMODES
      LworkL=3*NCV*NCV+6*NCV
# endif
# ifdef SO_SEMI
      Nsemi=1+ntimes(ng)/nADJ(ng)
# endif

      IF (ng.eq.1) THEN

        allocate ( select(NCV) )

        allocate ( Bvec(Nstr(ng):Nend(ng),NCV) )

        allocate ( RvalueR(NEV+1) )

        allocate ( RvalueI(NEV+1) )

        allocate ( Rvector(Nstr(ng):Nend(ng),NEV+1) )

# ifdef DISTRIBUTE
        allocate ( Swork(Mstate(ng)) )
# endif

        allocate ( SworkD(3*Nstate(ng)) )

        allocate ( SworkEV(3*NCV) )

        allocate ( SworkL(LworkL) )

        allocate ( norm(NEV+1) )

        allocate ( resid(Nstr(ng):Nend(ng)) )

# ifdef SO_SEMI
        allocate ( so_state(Nstr(ng):Nend(ng), Nsemi) )
# endif

      END IF

!
!-----------------------------------------------------------------------
!  Initialize module variables.
!-----------------------------------------------------------------------
!
      DO j=1,NCV
        select(j) = .TRUE.
        DO i=Nstr(ng),Nend(ng)
          Bvec(i,j) = IniVal
        END DO
      END DO
      DO j=1,NEV+1
        norm(j) = IniVal
        RvalueR(j) = IniVal
        RvalueI(j) = IniVal
        DO i=Nstr(ng),Nend(ng)
          Rvector(i,j) = IniVal
        END DO
      END DO
# ifdef DISTRIBUTE
      DO i=1,Mstate(ng)
        Swork(i) = IniVal
      END DO
# endif
      DO i=1,3*Nstate(ng)
        SworkD(i) = IniVal
      END DO
      DO i=1,3*NCV
        SworkEV(i) = IniVal
      END DO
      DO i=1,LworkL
        SworkL(i) = IniVal
      END DO
      DO i=Nstr(ng),Nend(ng)
        resid(i) = IniVal
      END DO
# ifdef SO_SEMI
      DO j=1,Nsemi
        DO i=Nstr(ng),Nend(ng)
          so_state(i,j) = IniVal
        END DO
      END DO
# endif

      RETURN
      END SUBROUTINE allocate_storage
#endif
      END MODULE mod_storage
