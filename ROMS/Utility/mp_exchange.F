#include "cppdefs.h"
      MODULE mp_exchange_mod
#ifdef DISTRIBUTE
!
!svn $Id: mp_exchange.F 795 2016-05-11 01:42:43Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  Distributed-memory tile exchange:                                   !
!                                                                      !
!  This routine updates the I,J tile overlap halo of NV variables.     !
!  It exchanges the specified number of "ghost-points".  In  order     !
!  to minimize the number send and receive calls, the ghost-points     !
!  are included in the buffers.  Therefore, the order of the pack,     !
!  send, receive, and unpack is crucial.                               !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng          Nested grid number.                                  !
!     model       Calling model identifier.                            !
!     Nvar        Number of variables for aggregated exchanges.        !
!     Istr        Starting tile index in the I-direction.              !
!     Iend        Ending   tile index in the I-direction.              !
!     Jstr        Starting tile index in the J-direction.              !
!     Jend        Ending   tile index in the J-direction.              !
!     LBi         I-dimension Lower bound.                             !
!     UBi         I-dimension Upper bound.                             !
!     LBj         J-dimension Lower bound.                             !
!     UBj         J-dimension Upper bound.                             !
!     LBk         K-dimension Lower bound.                             !
!     UBk         K-dimension Upper bound.                             !
!     LBt         T-dimension Lower bound.                             !
!     UBt         T-dimension Upper bound.                             !
!     Nghost      Number of ghost-points in the halo region.           !
!     EW_periodic Switch indicating EW periodicity exchanges.          !
!     NS_periodic Switch indicating NS periodicity exchanges.          !
!     A           2D tiled array to process.                           !
!     B           2D tiled array (optional) to process.                !
!     C           2D tiled array (optional) to process.                !
!     D           2D tiled array (optional) to process.                !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     A           Updated tiled array.                                 !
!     B           Updated tiled array (optional).                      !
!     C           Updated tiled array (optional).                      !
!     D           Updated tiled array (optional).                      !
!                                                                      !
!  Routines:                                                           !
!                                                                      !
!  mp_exchange2d         2D variables tile exchanges                   !
!  mp_exchange2d_bry     2D boundary variables tile exchanges          !
!  mp_exchange3d         3D variables tile exchanges                   !
!  mp_exchange3d_bry     3D boundary variables tile exchanges          !
!  mp_exchange4d         4D variables tile exchanges                   !
!                                                                      !
!  ad_mp_exchange2d      2D variables tile adjoint exchanges           !
!  ad_mp_exchange2d_bry  2D boundary variables tile adjoint exchanges  !
!  ad_mp_exchange3d      3D variables tile adjoint exchanges           !
!  ad_mp_exchange3d_bry  3D boundary variables tile adjoint exchanges  !
!  ad_mp_exchange4d      4D variables tile adjoint exchanges           !
!                                                                      !
!=======================================================================
!
      implicit none

      CONTAINS
!
!***********************************************************************
      SUBROUTINE tile_neighbors (ng, Nghost, EW_periodic, NS_periodic,  &
     &                           GrecvW, GsendW, Wtile, Wexchange,      &
     &                           GrecvE, GsendE, Etile, Eexchange,      &
     &                           GrecvS, GsendS, Stile, Sexchange,      &
     &                           GrecvN, GsendN, Ntile, Nexchange)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
!
      implicit none
!
!  Imported variable declarations.
!
      logical, intent(in) :: EW_periodic, NS_periodic

      integer, intent(in) :: ng, Nghost

      logical, intent(out) :: Wexchange, Eexchange
      logical, intent(out) :: Sexchange, Nexchange

      integer, intent(out) :: GrecvW, GsendW, Wtile
      integer, intent(out) :: GrecvE, GsendE, Etile
      integer, intent(out) :: GrecvS, GsendS, Stile
      integer, intent(out) :: GrecvN, GsendN, Ntile
!
!  Local variable declarations.
!
      integer :: i, j
      integer :: MyRankI, MyRankJ, Null_Value, rank

      integer, dimension(-1:NtileI(ng),-1:NtileJ(ng)) :: table
!
!-----------------------------------------------------------------------
!  Set tile partition table for looking up adjacent processes.
!-----------------------------------------------------------------------
!
!  Notice that a null value is used in places that data transmition is
!  not required.
!
# if defined MPI
      Null_Value=MPI_PROC_NULL
# else
      Null_Value=-1
# endif
      DO j=-1,NtileJ(ng)
        DO i=-1,NtileI(ng)
          table(i,j)=Null_Value
        END DO
      END DO
      rank=0
      DO j=0,NtileJ(ng)-1
        DO i=0,NtileI(ng)-1
          table(i,j)=rank
          IF (MyRank.eq.rank) THEN
            MyRankI=i
            MyRankJ=j
          END IF
          rank=rank+1
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Determine the rank of Western and Eastern tiles.  Then, determine
!  the number of ghost-points to send and receive in the West- and
!  East-directions.
!-----------------------------------------------------------------------
!
!  This logic only works for two and three ghost points. The number of
!  ghost-points changes when periodic boundary condition are activated.
!  The periodicity is as follows:
!
!  If two ghost-points:
!
!                      Lm-2  Lm-1  Lm   Lm+1  Lm+2
!                      -2    -1     0    1     2
!
!  If three ghost-points:
!
!                      Lm-2  Lm-1  Lm   Lm+1  Lm+2   Lm+3
!                      -2    -1     0    1     2      3
!
      IF (EW_periodic) THEN
        IF ((table(MyRankI-1,MyRankJ).eq.Null_Value).and.               &
     &      (NtileI(ng).gt.1)) THEN
          Wtile=table(NtileI(ng)-1,MyRankJ)
          Etile=table(MyRankI+1,MyRankJ)
          GsendW=Nghost
          GsendE=Nghost
          IF (NghostPoints.eq.3) THEN
            GrecvW=Nghost
          ELSE
            GrecvW=Nghost+1
          END IF
          GrecvE=Nghost
        ELSE IF ((table(MyRankI+1,MyRankJ).eq.Null_Value).and.          &
     &           (NtileI(ng).gt.1)) THEN
          Wtile=table(MyRankI-1,MyRankJ)
          Etile=table(0,MyRankJ)
          GsendW=Nghost
          IF (NghostPoints.eq.3) THEN
            GsendE=Nghost
          ELSE
            GsendE=Nghost+1
          END IF
          GrecvW=Nghost
          GrecvE=Nghost
        ELSE
          Wtile=table(MyRankI-1,MyRankJ)
          Etile=table(MyRankI+1,MyRankJ)
          GsendW=Nghost
          GsendE=Nghost
          GrecvW=Nghost
          GrecvE=Nghost
        END IF
      ELSE
        Wtile=table(MyRankI-1,MyRankJ)
        Etile=table(MyRankI+1,MyRankJ)
        GsendW=Nghost
        GsendE=Nghost
        GrecvW=Nghost
        GrecvE=Nghost
      END IF
!
!  Determine exchange switches.
!
      IF (Wtile.eq.Null_Value) THEN
        Wexchange=.FALSE.
      ELSE
        Wexchange=.TRUE.
      END IF
      IF (Etile.eq.Null_Value) THEN
        Eexchange=.FALSE.
      ELSE
        Eexchange=.TRUE.
      END IF
!
!-----------------------------------------------------------------------
!  Determine the rank of Southern and Northern tiles.  Then, determine
!  the number of ghost-points to send and receive in the South- and
!  North-directions.
!-----------------------------------------------------------------------
!
!  This logic only works for two and three ghost-points. The number of
!  ghost-points changes when periodic boundary condition are activated.
!  The periodicity is as follows:
!
!  If two ghost-points:
!
!                      Mm-2  Mm-1  Mm   Mm+1  Mm+2
!                      -2    -1     0    1     2
!
!  If three ghost-points:
!
!                      Mm-2  Mm-1  Mm   Mm+1  Mm+2  Mm+3
!                      -2    -1     0    1     2     3
!
      IF (NS_periodic) THEN
        IF ((table(MyRankI,MyRankJ-1).eq.Null_Value).and.               &
     &      (NtileJ(ng).gt.1)) THEN
          Stile=table(MyRankI,NtileJ(ng)-1)
          Ntile=table(MyRankI,MyRankJ+1)
          GsendS=Nghost
          GsendN=Nghost
          IF (NghostPoints.eq.3) THEN
            GrecvS=Nghost
          ELSE
            GrecvS=Nghost+1
          END IF
          GrecvN=Nghost
        ELSE IF ((table(MyRankI,MyRankJ+1).eq.Null_Value).and.          &
     &           (NtileJ(ng).gt.1)) then
          Stile=table(MyRankI,MyRankJ-1)
          Ntile=table(MyRankI,0)
          GsendS=Nghost
          IF (NghostPoints.eq.3) THEN
            GsendN=Nghost
          ELSE
            GsendN=Nghost+1
          END IF
          GrecvS=Nghost
          GrecvN=Nghost
        ELSE
          Stile=table(MyRankI,MyRankJ-1)
          Ntile=table(MyRankI,MyRankJ+1)
          GsendS=Nghost
          GsendN=Nghost
          GrecvS=Nghost
          GrecvN=Nghost
        END IF
      ELSE
        Stile=table(MyRankI,MyRankJ-1)
        Ntile=table(MyRankI,MyRankJ+1)
        GsendS=Nghost
        GsendN=Nghost
        GrecvS=Nghost
        GrecvN=Nghost
      END IF
!
!  Determine exchange switches.
!
      IF (Stile.eq.Null_Value) THEN
        Sexchange=.FALSE.
      ELSE
        Sexchange=.TRUE.
      END IF
      IF (Ntile.eq.Null_Value) THEN
        Nexchange=.FALSE.
      ELSE
        Nexchange=.TRUE.
      END IF

      RETURN
      END SUBROUTINE tile_neighbors

!
!***********************************************************************
      SUBROUTINE mp_exchange2d (ng, tile, model, Nvar,                  &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Nghost, EW_periodic, NS_periodic,       &
     &                          A, B, C, D)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
      implicit none
!
!  Imported variable declarations.
!
      logical, intent(in) :: EW_periodic, NS_periodic

      integer, intent(in) :: ng, tile, model, Nvar
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: Nghost
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: A(LBi:,LBj:)

      real(r8), intent(inout), optional :: B(LBi:,LBj:)
      real(r8), intent(inout), optional :: C(LBi:,LBj:)
      real(r8), intent(inout), optional :: D(LBi:,LBj:)
# else
      real(r8), intent(inout) :: A(LBi:UBi,LBj:UBj)

      real(r8), intent(inout), optional :: B(LBi:UBi,LBj:UBj)
      real(r8), intent(inout), optional :: C(LBi:UBi,LBj:UBj)
      real(r8), intent(inout), optional :: D(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
      logical :: Wexchange, Sexchange, Eexchange, Nexchange

      integer :: i, icS, icN, ioff, Imin, Imax, Ilen
      integer :: j, jcW, jcE, joff, Jmin, Jmax, Jlen
      integer :: m, mc, Ierror, Lstr, pp
      integer :: Wtile, GsendW, GrecvW, Wtag, Werror, Wrequest
      integer :: Stile, GsendS, GrecvS, Stag, Serror, Srequest
      integer :: Etile, GsendE, GrecvE, Etag, Eerror, Erequest
      integer :: Ntile, GsendN, GrecvN, Ntag, Nerror, Nrequest
      integer :: EWsize, sizeW, sizeE
      integer :: NSsize, sizeS, sizeN

# ifdef MPI
      integer, dimension(MPI_STATUS_SIZE,4) :: status
# endif
!
      real(r8), dimension(Nvar*HaloSizeJ(ng)) :: sendW, sendE
      real(r8), dimension(Nvar*HaloSizeJ(ng)) :: recvW, recvE

      real(r8), dimension(Nvar*HaloSizeI(ng)) :: sendS, sendN
      real(r8), dimension(Nvar*HaloSizeI(ng)) :: recvS, recvN

      character (len=MPI_MAX_ERROR_STRING) :: string

# include "set_bounds.h"

# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn on time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, model, 40)
# endif
!
!-----------------------------------------------------------------------
!  Determine rank of tile neighbors and number of ghost-points to
!  exchange.
!-----------------------------------------------------------------------
!
      CALL tile_neighbors (ng, Nghost, EW_periodic, NS_periodic,        &
     &                     GrecvW, GsendW, Wtile, Wexchange,            &
     &                     GrecvE, GsendE, Etile, Eexchange,            &
     &                     GrecvS, GsendS, Stile, Sexchange,            &
     &                     GrecvN, GsendN, Ntile, Nexchange)
!
!  Set communication tags.
!
      Wtag=1
      Stag=2
      Etag=3
      Ntag=4
!
!  Determine range and length of the distributed tile boundary segments.
!
      Imin=LBi
      Imax=UBi
      Jmin=LBj
      Jmax=UBj
      Ilen=Imax-Imin+1
      Jlen=Jmax-Jmin+1
      IF (EW_periodic.or.NS_periodic) THEN
        pp=1
      ELSE
        pp=0
      END IF
      EWsize=Nvar*(Nghost+pp)*Jlen
      NSsize=Nvar*(Nghost+pp)*Ilen
      IF (SIZE(sendE).lt.EWsize) THEN
        WRITE (stdout,10) 'EWsize = ', EWsize, SIZE(sendE)
 10     FORMAT (/,' MP_EXCHANGE2D - communication buffer too small, ',  &
     &          a, 2i8)
      END IF
      IF (SIZE(sendN).lt.NSsize) THEN
        WRITE (stdout,10) 'NSsize = ', NSsize, SIZE(sendN)
      END IF
!
!-----------------------------------------------------------------------
!  Pack Western and Eastern tile boundary data including ghost-points.
!-----------------------------------------------------------------------
!
      IF (Wexchange) THEN
        sizeW=0
        DO m=1,GsendW
          mc=(m-1)*Jlen
          i=Istr+m-1
          DO j=Jmin,Jmax
            sizeW=sizeW+1
            jcW=1+(j-Jmin)+mc
            sendW(jcW)=A(i,j)
          END DO
        END DO
        IF (PRESENT(B)) THEN
          joff=jcW
          DO m=1,GsendW
            mc=(m-1)*Jlen
            i=Istr+m-1
            DO j=Jmin,Jmax
              sizeW=sizeW+1
              jcW=joff+1+(j-Jmin)+mc
              sendW(jcW)=B(i,j)
            END DO
          END DO
        END IF
        IF (PRESENT(C)) THEN
          joff=jcW
          DO m=1,GsendW
            mc=(m-1)*Jlen
            i=Istr+m-1
            DO j=Jmin,Jmax
              sizeW=sizeW+1
              jcW=joff+1+(j-Jmin)+mc
              sendW(jcW)=C(i,j)
            END DO
          END DO
        END IF
        IF (PRESENT(D)) THEN
          joff=jcW
          DO m=1,GsendW
            mc=(m-1)*Jlen
            i=Istr+m-1
            DO j=Jmin,Jmax
              sizeW=sizeW+1
              jcW=joff+1+(j-Jmin)+mc
              sendW(jcW)=D(i,j)
            END DO
          END DO
        END IF
      END IF
!
      IF (Eexchange) THEN
        sizeE=0
        DO m=1,GsendE
          mc=(m-1)*Jlen
          i=Iend-GsendE+m
          DO j=Jmin,Jmax
            sizeE=sizeE+1
            jcE=1+(j-Jmin)+mc
            sendE(jcE)=A(i,j)
          END DO
        END DO
        IF (PRESENT(B)) THEN
          joff=jcE
          DO m=1,GsendE
            mc=(m-1)*Jlen
            i=Iend-GsendE+m
            DO j=Jmin,Jmax
              sizeE=sizeE+1
              jcE=joff+1+(j-Jmin)+mc
              sendE(jcE)=B(i,j)
            END DO
          END DO
        END IF
        IF (PRESENT(C)) THEN
          joff=jcE
          DO m=1,GsendE
            mc=(m-1)*Jlen
            i=Iend-GsendE+m
            DO j=Jmin,Jmax
              sizeE=sizeE+1
              jcE=joff+1+(j-Jmin)+mc
              sendE(jcE)=C(i,j)
            END DO
          END DO
        END IF
        IF (PRESENT(D)) THEN
          joff=jcE
          DO m=1,GsendE
            mc=(m-1)*Jlen
            i=Iend-GsendE+m
            DO j=Jmin,Jmax
              sizeE=sizeE+1
              jcE=joff+1+(j-Jmin)+mc
              sendE(jcE)=D(i,j)
            END DO
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Send and receive Western and Eastern segments.
!-----------------------------------------------------------------------
!
# if defined MPI
      IF (Wexchange) THEN
        CALL mpi_irecv (recvW, EWsize, MP_FLOAT, Wtile, Etag,           &
     &                  OCN_COMM_WORLD, Wrequest, Werror)
      END IF
      IF (Eexchange) THEN
        CALL mpi_irecv (recvE, EWsize, MP_FLOAT, Etile, Wtag,           &
     &                  OCN_COMM_WORLD, Erequest, Eerror)
      END IF
      IF (Wexchange) THEN
        CALL mpi_send  (sendW, sizeW, MP_FLOAT, Wtile, Wtag,            &
     &                  OCN_COMM_WORLD, Werror)
      END IF
      IF (Eexchange) THEN
        CALL mpi_send  (sendE, sizeE, MP_FLOAT, Etile, Etag,            &
     &                  OCN_COMM_WORLD, Eerror)
      END IF
# endif
!
!-----------------------------------------------------------------------
!  Unpack Western and Eastern segments.
!-----------------------------------------------------------------------
!
      IF (Wexchange) THEN
# ifdef MPI
        CALL mpi_wait (Wrequest, status(1,1), Werror)
        IF (Werror.ne.MPI_SUCCESS) THEN
          CALL mpi_error_string (Werror, string, Lstr, Ierror)
          Lstr=LEN_TRIM(string)
          WRITE (stdout,20) 'MPI_SEND/MPI_IRECV (Western Edge)',        &
     &                      MyRank, Werror, string(1:Lstr)
 20       FORMAT (/,' MP_EXCHANGE2D - error during ',a,                 &
     &            ' call, Node = ',i3.3,' Error = ',i3,/,15x,a)
          exit_flag=2
          RETURN
        END IF
# endif
        DO m=GrecvW,1,-1
          mc=(GrecvW-m)*Jlen
          i=Istr-m
          DO j=Jmin,Jmax
            jcW=1+(j-Jmin)+mc
            A(i,j)=recvW(jcW)
          END DO
        END DO
        IF (PRESENT(B)) THEN
          joff=jcW
          DO m=GrecvW,1,-1
            mc=(GrecvW-m)*Jlen
            i=Istr-m
            DO j=Jmin,Jmax
              jcW=joff+1+(j-Jmin)+mc
              B(i,j)=recvW(jcW)
            END DO
          END DO
        END IF
        IF (PRESENT(C)) THEN
          joff=jcW
          DO m=GrecvW,1,-1
            mc=(GrecvW-m)*Jlen
            i=Istr-m
            DO j=Jmin,Jmax
              jcW=joff+1+(j-Jmin)+mc
              C(i,j)=recvW(jcW)
            END DO
          END DO
        END IF
        IF (PRESENT(D)) THEN
          joff=jcW
          DO m=GrecvW,1,-1
            mc=(GrecvW-m)*Jlen
            i=Istr-m
            DO j=Jmin,Jmax
              jcW=joff+1+(j-Jmin)+mc
              D(i,j)=recvW(jcW)
            END DO
          END DO
        END IF
      END IF
!
      IF (Eexchange) THEN
# ifdef MPI
        CALL mpi_wait (Erequest, status(1,3), Eerror)
        IF (Eerror.ne.MPI_SUCCESS) THEN
          CALL mpi_error_string (Eerror, string, Lstr, Ierror)
          Lstr=LEN_TRIM(string)
          WRITE (stdout,20) 'MPI_SEND/MPI_IRECV (Eastern Edge)',        &
     &                      MyRank, Eerror, string(1:Lstr)
          exit_flag=2
          RETURN
        END IF
# endif
        DO m=1,GrecvE
          mc=(m-1)*Jlen
          i=Iend+m
          DO j=Jmin,Jmax
            jcE=1+(j-Jmin)+mc
            A(i,j)=recvE(jcE)
          ENDDO
        END DO
        IF (PRESENT(B)) THEN
          joff=jcE
          DO m=1,GrecvE
            mc=(m-1)*Jlen
            i=Iend+m
            DO j=Jmin,Jmax
              jcE=joff+1+(j-Jmin)+mc
              B(i,j)=recvE(jcE)
            ENDDO
          END DO
        END IF
        IF (PRESENT(C)) THEN
          joff=jcE
          DO m=1,GrecvE
            mc=(m-1)*Jlen
            i=Iend+m
            DO j=Jmin,Jmax
              jcE=joff+1+(j-Jmin)+mc
              C(i,j)=recvE(jcE)
            ENDDO
          END DO
        END IF
        IF (PRESENT(D)) THEN
          joff=jcE
          DO m=1,GrecvE
            mc=(m-1)*Jlen
            i=Iend+m
            DO j=Jmin,Jmax
              jcE=joff+1+(j-Jmin)+mc
              D(i,j)=recvE(jcE)
            ENDDO
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Pack Southern and Northern tile boundary data including ghost-points.
!-----------------------------------------------------------------------
!
      IF (Sexchange) THEN
        sizeS=0
        DO m=1,GsendS
          mc=(m-1)*Ilen
          j=Jstr+m-1
          DO i=Imin,Imax
            sizeS=sizeS+1
            icS=1+(i-Imin)+mc
            sendS(icS)=A(i,j)
          END DO
        END DO
        IF (PRESENT(B)) THEN
          ioff=icS
          DO m=1,GsendS
            mc=(m-1)*Ilen
            j=Jstr+m-1
            DO i=Imin,Imax
              sizeS=sizeS+1
              icS=ioff+1+(i-Imin)+mc
              sendS(icS)=B(i,j)
            END DO
          END DO
        END IF
        IF (PRESENT(C)) THEN
          ioff=icS
          DO m=1,GsendS
            mc=(m-1)*Ilen
            j=Jstr+m-1
            DO i=Imin,Imax
              sizeS=sizeS+1
              icS=ioff+1+(i-Imin)+mc
              sendS(icS)=C(i,j)
            END DO
          END DO
        END IF
        IF (PRESENT(D)) THEN
          ioff=icS
          DO m=1,GsendS
            mc=(m-1)*Ilen
            j=Jstr+m-1
            DO i=Imin,Imax
              sizeS=sizeS+1
              icS=ioff+1+(i-Imin)+mc
              sendS(icS)=D(i,j)
            END DO
          END DO
        END IF
      END IF
!
      IF (Nexchange) THEN
        sizeN=0
        DO m=1,GsendN
          mc=(m-1)*Ilen
          j=Jend-GsendN+m
          DO i=Imin,Imax
            sizeN=sizeN+1
            icN=1+(i-Imin)+mc
            sendN(icN)=A(i,j)
          END DO
        END DO
        IF (PRESENT(B)) THEN
          ioff=icN
          DO m=1,GsendN
            mc=(m-1)*Ilen
            j=Jend-GsendN+m
            DO i=Imin,Imax
              sizeN=sizeN+1
              icN=ioff+1+(i-Imin)+mc
              sendN(icN)=B(i,j)
            END DO
          END DO
        END IF
        IF (PRESENT(C)) THEN
          ioff=icN
          DO m=1,GsendN
            mc=(m-1)*Ilen
            j=Jend-GsendN+m
            DO i=Imin,Imax
              sizeN=sizeN+1
              icN=ioff+1+(i-Imin)+mc
              sendN(icN)=C(i,j)
            END DO
          END DO
        END IF
        IF (PRESENT(D)) THEN
          ioff=icN
          DO m=1,GsendN
            mc=(m-1)*Ilen
            j=Jend-GsendN+m
            DO i=Imin,Imax
              sizeN=sizeN+1
              icN=ioff+1+(i-Imin)+mc
              sendN(icN)=D(i,j)
            END DO
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Send and receive Southern and Northern segments.
!-----------------------------------------------------------------------
!
# if defined MPI
      IF (Sexchange) THEN
        CALL mpi_irecv (recvS, NSsize, MP_FLOAT, Stile, Ntag,           &
     &                  OCN_COMM_WORLD, Srequest, Serror)
      END IF
      IF (Nexchange) THEN
        CALL mpi_irecv (recvN, NSsize, MP_FLOAT, Ntile, Stag,           &
     &                  OCN_COMM_WORLD, Nrequest, Nerror)
      END IF
      IF (Sexchange) THEN
        CALL mpi_send  (sendS, sizeS, MP_FLOAT, Stile, Stag,            &
     &                  OCN_COMM_WORLD, Serror)
      END IF
      IF (Nexchange) THEN
        CALL mpi_send  (sendN, sizeN, MP_FLOAT, Ntile, Ntag,            &
     &                  OCN_COMM_WORLD, Nerror)
      END IF
# endif
!
!-----------------------------------------------------------------------
!  Unpack Northern and Southern segments.
!-----------------------------------------------------------------------
!
      IF (Sexchange) THEN
# ifdef MPI
        CALL mpi_wait (Srequest, status(1,2), Serror)
        IF (Serror.ne.MPI_SUCCESS) THEN
          CALL mpi_error_string (Serror, string, Lstr, Ierror)
          Lstr=LEN_TRIM(string)
          WRITE (stdout,20) 'MPI_SEND/MPI_IRECV (Southern Edge)',       &
     &                      MyRank, Serror, string(1:Lstr)
          exit_flag=2
          RETURN
        END IF
# endif
        DO m=GrecvS,1,-1
          mc=(GrecvS-m)*Ilen
          j=Jstr-m
          DO i=Imin,Imax
            icS=1+(i-Imin)+mc
            A(i,j)=recvS(icS)
          END DO
        END DO
        IF (PRESENT(B)) THEN
          ioff=icS
          DO m=GrecvS,1,-1
            mc=(GrecvS-m)*Ilen
            j=Jstr-m
            DO i=Imin,Imax
              icS=ioff+1+(i-Imin)+mc
              B(i,j)=recvS(icS)
            END DO
          END DO
        END IF
        IF (PRESENT(C)) THEN
          ioff=icS
          DO m=GrecvS,1,-1
            mc=(GrecvS-m)*Ilen
            j=Jstr-m
            DO i=Imin,Imax
              icS=ioff+1+(i-Imin)+mc
              C(i,j)=recvS(icS)
            END DO
          END DO
        END IF
        IF (PRESENT(D)) THEN
          ioff=icS
          DO m=GrecvS,1,-1
            mc=(GrecvS-m)*Ilen
            j=Jstr-m
            DO i=Imin,Imax
              icS=ioff+1+(i-Imin)+mc
              D(i,j)=recvS(icS)
            END DO
          END DO
        END IF
      END IF
!
      IF (Nexchange) THEN
# ifdef MPI
        CALL mpi_wait (Nrequest, status(1,4), Nerror)
        IF (Nerror.ne.MPI_SUCCESS) THEN
          CALL mpi_error_string (Nerror, string, Lstr, Ierror)
          Lstr=LEN_TRIM(string)
          WRITE (stdout,20) 'MPI_SEND/MPI_IRECV (Northern Edge)',       &
     &                      MyRank, Nerror, string(1:Lstr)
          exit_flag=2
          RETURN
        END IF
# endif
        DO m=1,GrecvN
          mc=(m-1)*Ilen
          j=Jend+m
          DO i=Imin,Imax
            icN=1+(i-Imin)+mc
            A(i,j)=recvN(icN)
          END DO
        END DO
        IF (PRESENT(B)) THEN
          ioff=icN
          DO m=1,GrecvN
            mc=(m-1)*Ilen
            j=Jend+m
            DO i=Imin,Imax
              icN=ioff+1+(i-Imin)+mc
              B(i,j)=recvN(icN)
            END DO
          END DO
        END IF
        IF (PRESENT(C)) THEN
          ioff=icN
          DO m=1,GrecvN
            mc=(m-1)*Ilen
            j=Jend+m
            DO i=Imin,Imax
              icN=ioff+1+(i-Imin)+mc
              C(i,j)=recvN(icN)
            END DO
          END DO
        END IF
        IF (PRESENT(D)) THEN
          ioff=icN
          DO m=1,GrecvN
            mc=(m-1)*Ilen
            j=Jend+m
            DO i=Imin,Imax
              icN=ioff+1+(i-Imin)+mc
              D(i,j)=recvN(icN)
            END DO
          END DO
        END IF
      END IF
# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn off time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, model, 40)
# endif

      RETURN
      END SUBROUTINE mp_exchange2d
!
!***********************************************************************
      SUBROUTINE mp_exchange2d_bry (ng, tile, model, Nvar, boundary,    &
     &                              LBij, UBij,                         &
     &                              Nghost, EW_periodic, NS_periodic,   &
     &                              A, B, C, D)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
      implicit none
!
!  Imported variable declarations.
!
      logical, intent(in) :: EW_periodic, NS_periodic

      integer, intent(in) :: ng, tile, model, Nvar, boundary
      integer, intent(in) :: LBij, UBij
      integer, intent(in) :: Nghost
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: A(LBij:)

      real(r8), intent(inout), optional :: B(LBij:)
      real(r8), intent(inout), optional :: C(LBij:)
      real(r8), intent(inout), optional :: D(LBij:)
# else
      real(r8), intent(inout) :: A(LBij:UBij)

      real(r8), intent(inout), optional :: B(LBij:UBij)
      real(r8), intent(inout), optional :: C(LBij:UBij)
      real(r8), intent(inout), optional :: D(LBij:UBij)
# endif
!
!  Local variable declarations.
!
      logical :: Wexchange, Sexchange, Eexchange, Nexchange

      integer :: i, icS, icN
      integer :: j, jcW, jcE
      integer :: m, Ierror, Lstr, pp
      integer :: Wtile, GsendW, GrecvW, Wtag, Werror, Wrequest
      integer :: Stile, GsendS, GrecvS, Stag, Serror, Srequest
      integer :: Etile, GsendE, GrecvE, Etag, Eerror, Erequest
      integer :: Ntile, GsendN, GrecvN, Ntag, Nerror, Nrequest
      integer :: EWsize, sizeW, sizeE
      integer :: NSsize, sizeS, sizeN

# ifdef MPI
      integer, dimension(MPI_STATUS_SIZE,4) :: status
# endif
!
      real(r8), dimension(Nvar*HaloBry(ng)) :: sendW, sendE
      real(r8), dimension(Nvar*HaloBry(ng)) :: recvW, recvE

      real(r8), dimension(Nvar*HaloBry(ng)) :: sendS, sendN
      real(r8), dimension(Nvar*HaloBry(ng)) :: recvS, recvN

      character (len=MPI_MAX_ERROR_STRING) :: string

# include "set_bounds.h"

# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn on time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, model, 43)
# endif
!
!-----------------------------------------------------------------------
!  Determine rank of tile neighbors and number of ghost-points to
!  exchange.
!-----------------------------------------------------------------------
!
      CALL tile_neighbors (ng, Nghost, EW_periodic, NS_periodic,        &
     &                     GrecvW, GsendW, Wtile, Wexchange,            &
     &                     GrecvE, GsendE, Etile, Eexchange,            &
     &                     GrecvS, GsendS, Stile, Sexchange,            &
     &                     GrecvN, GsendN, Ntile, Nexchange)
!
!  Adjust exchange swiches according to boundary edge to process.
!
      Wexchange=Wexchange.and.((boundary.eq.isouth).or.                 &
     &                         (boundary.eq.inorth))
      Eexchange=Eexchange.and.((boundary.eq.isouth).or.                 &
     &                         (boundary.eq.inorth))
      Sexchange=Sexchange.and.((boundary.eq.iwest).or.                  &
     &                         (boundary.eq.ieast))
      Nexchange=Nexchange.and.((boundary.eq.iwest).or.                  &
     &                         (boundary.eq.ieast))
!
!  Set communication tags.
!
      Wtag=1
      Stag=2
      Etag=3
      Ntag=4
!
!  Determine range and length of the distributed tile boundary segments.
!
      IF (EW_periodic.or.NS_periodic) THEN
        pp=1
      ELSE
        pp=0
      END IF
      EWsize=Nvar*(Nghost+pp)
      NSsize=Nvar*(Nghost+pp)
      IF (SIZE(sendE).lt.EWsize) THEN
        WRITE (stdout,10) 'EWsize = ', EWsize, SIZE(sendE)
 10     FORMAT (/,' MP_EXCHANGE2D_BRY - communication buffer too ',     &
     &          'small, ',a, 2i8)
      END IF
      IF (SIZE(sendN).lt.NSsize) THEN
        WRITE (stdout,10) 'NSsize = ', NSsize, SIZE(sendN)
      END IF
!
!-----------------------------------------------------------------------
!  Pack Western and Eastern tile boundary data including ghost-points.
!-----------------------------------------------------------------------
!
      IF (Wexchange) THEN
        jcW=0
        sizeW=0
        DO m=1,GsendW
          i=Istr+m-1
          sizeW=sizeW+1
          jcW=jcW+1
          sendW(jcW)=A(i)
        END DO
        IF (PRESENT(B)) THEN
          DO m=1,GsendW
            i=Istr+m-1
            sizeW=sizeW+1
            jcW=jcW+1
            sendW(jcW)=B(i)
          END DO
        END IF
        IF (PRESENT(C)) THEN
          DO m=1,GsendW
            i=Istr+m-1
            sizeW=sizeW+1
            jcW=jcW+1
            sendW(jcW)=C(i)
          END DO
        END IF
        IF (PRESENT(D)) THEN
          DO m=1,GsendW
            i=Istr+m-1
            sizeW=sizeW+1
            jcW=jcW+1
            sendW(jcW)=D(i)
          END DO
        END IF
      END IF
!
      IF (Eexchange) THEN
        jcE=0
        sizeE=0
        DO m=1,GsendE
          i=Iend-GsendE+m
          sizeE=sizeE+1
          jcE=jcE+1
          sendE(jcE)=A(i)
        END DO
        IF (PRESENT(B)) THEN
          DO m=1,GsendE
            i=Iend-GsendE+m
            sizeE=sizeE+1
            jcE=jcE+1
            sendE(jcE)=B(i)
          END DO
        END IF
        IF (PRESENT(C)) THEN
          DO m=1,GsendE
            i=Iend-GsendE+m
            sizeE=sizeE+1
            jcE=jcE+1
            sendE(jcE)=C(i)
          END DO
        END IF
        IF (PRESENT(D)) THEN
          DO m=1,GsendE
            i=Iend-GsendE+m
            sizeE=sizeE+1
            jcE=jcE+1
            sendE(jcE)=D(i)
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Send and receive Western and Eastern segments.
!-----------------------------------------------------------------------
!
# if defined MPI
      IF (Wexchange) THEN
        CALL mpi_irecv (recvW, EWsize, MP_FLOAT, Wtile, Etag,           &
     &                  OCN_COMM_WORLD, Wrequest, Werror)
      END IF
      IF (Eexchange) THEN
        CALL mpi_irecv (recvE, EWsize, MP_FLOAT, Etile, Wtag,           &
     &                  OCN_COMM_WORLD, Erequest, Eerror)
      END IF
      IF (Wexchange) THEN
        CALL mpi_send  (sendW, sizeW, MP_FLOAT, Wtile, Wtag,            &
     &                  OCN_COMM_WORLD, Werror)
      END IF
      IF (Eexchange) THEN
        CALL mpi_send  (sendE, sizeE, MP_FLOAT, Etile, Etag,            &
     &                  OCN_COMM_WORLD, Eerror)
      END IF
# endif
!
!-----------------------------------------------------------------------
!  Unpack Western and Eastern segments.
!-----------------------------------------------------------------------
!
      IF (Wexchange) THEN
# ifdef MPI
        CALL mpi_wait (Wrequest, status(1,1), Werror)
        IF (Werror.ne.MPI_SUCCESS) THEN
          CALL mpi_error_string (Werror, string, Lstr, Ierror)
          Lstr=LEN_TRIM(string)
          WRITE (stdout,20) 'MPI_SEND/MPI_IRECV (Western Edge)',        &
     &                      MyRank, Werror, string(1:Lstr)
 20       FORMAT (/,' MP_EXCHANGE2D_BRY - error during ',a,             &
     &            ' call, Node = ',i3.3,' Error = ',i3,/,15x,a)
          exit_flag=2
          RETURN
        END IF
# endif
        jcW=0
        DO m=GrecvW,1,-1
          i=Istr-m
          jcW=jcW+1
          A(i)=recvW(jcW)
        END DO
        IF (PRESENT(B)) THEN
          DO m=GrecvW,1,-1
            i=Istr-m
            jcW=jcW+1
            B(i)=recvW(jcW)
          END DO
        END IF
        IF (PRESENT(C)) THEN
          DO m=GrecvW,1,-1
            i=Istr-m
            jcW=jcW+1
            C(i)=recvW(jcW)
          END DO
        END IF
        IF (PRESENT(D)) THEN
          DO m=GrecvW,1,-1
            i=Istr-m
            jcW=jcW+1
            D(i)=recvW(jcW)
          END DO
        END IF
      END IF
!
      IF (Eexchange) THEN
# ifdef MPI
        CALL mpi_wait (Erequest, status(1,3), Eerror)
        IF (Eerror.ne.MPI_SUCCESS) THEN
          CALL mpi_error_string (Eerror, string, Lstr, Ierror)
          Lstr=LEN_TRIM(string)
          WRITE (stdout,20) 'MPI_SEND/MPI_IRECV (Eastern Edge)',        &
     &                      MyRank, Eerror, string(1:Lstr)
          exit_flag=2
          RETURN
        END IF
# endif
        jcE=0
        DO m=1,GrecvE
          i=Iend+m
          jcE=jcE+1
          A(i)=recvE(jcE)
        END DO
        IF (PRESENT(B)) THEN
          DO m=1,GrecvE
            i=Iend+m
            jcE=jcE+1
            B(i)=recvE(jcE)
          END DO
        END IF
        IF (PRESENT(C)) THEN
          DO m=1,GrecvE
            i=Iend+m
            jcE=jcE+1
            C(i)=recvE(jcE)
          END DO
        END IF
        IF (PRESENT(D)) THEN
          DO m=1,GrecvE
            i=Iend+m
            jcE=jcE+1
            D(i)=recvE(jcE)
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Pack Southern and Northern tile boundary data including ghost-points.
!-----------------------------------------------------------------------
!
      IF (Sexchange) THEN
        icS=0
        sizeS=0
        DO m=1,GsendS
          j=Jstr+m-1
          sizeS=sizeS+1
          icS=icS+1
          sendS(icS)=A(j)
        END DO
        IF (PRESENT(B)) THEN
          DO m=1,GsendS
            j=Jstr+m-1
            sizeS=sizeS+1
            icS=icS+1
            sendS(icS)=B(j)
          END DO
        END IF
        IF (PRESENT(C)) THEN
          DO m=1,GsendS
            j=Jstr+m-1
            sizeS=sizeS+1
            icS=icS+1
            sendS(icS)=C(j)
          END DO
        END IF
        IF (PRESENT(D)) THEN
          DO m=1,GsendS
            j=Jstr+m-1
            sizeS=sizeS+1
            icS=icS+1
            sendS(icS)=D(j)
          END DO
        END IF
      END IF
!
      IF (Nexchange) THEN
        icN=0
        sizeN=0
        DO m=1,GsendN
          j=Jend-GsendN+m
          sizeN=sizeN+1
          icN=icN+1
          sendN(icN)=A(j)
        END DO
        IF (PRESENT(B)) THEN
          DO m=1,GsendN
            j=Jend-GsendN+m
            sizeN=sizeN+1
            icN=icN+1
            sendN(icN)=B(j)
          END DO
        END IF
        IF (PRESENT(C)) THEN
          DO m=1,GsendN
            j=Jend-GsendN+m
            sizeN=sizeN+1
            icN=icN+1
            sendN(icN)=C(j)
          END DO
        END IF
        IF (PRESENT(D)) THEN
          DO m=1,GsendN
            j=Jend-GsendN+m
            sizeN=sizeN+1
            icN=icN+1
            sendN(icN)=D(j)
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Send and receive Southern and Northern segments.
!-----------------------------------------------------------------------
!
# if defined MPI
      IF (Sexchange) THEN
        CALL mpi_irecv (recvS, NSsize, MP_FLOAT, Stile, Ntag,           &
     &                  OCN_COMM_WORLD, Srequest, Serror)
      END IF
      IF (Nexchange) THEN
        CALL mpi_irecv (recvN, NSsize, MP_FLOAT, Ntile, Stag,           &
     &                  OCN_COMM_WORLD, Nrequest, Nerror)
      END IF
      IF (Sexchange) THEN
        CALL mpi_send  (sendS, sizeS, MP_FLOAT, Stile, Stag,            &
     &                  OCN_COMM_WORLD, Serror)
      END IF
      IF (Nexchange) THEN
        CALL mpi_send  (sendN, sizeN, MP_FLOAT, Ntile, Ntag,            &
     &                  OCN_COMM_WORLD, Nerror)
      END IF
# endif
!
!-----------------------------------------------------------------------
!  Unpack Northern and Southern segments.
!-----------------------------------------------------------------------
!
      IF (Sexchange) THEN
# ifdef MPI
        CALL mpi_wait (Srequest, status(1,2), Serror)
        IF (Serror.ne.MPI_SUCCESS) THEN
          CALL mpi_error_string (Serror, string, Lstr, Ierror)
          Lstr=LEN_TRIM(string)
          WRITE (stdout,20) 'MPI_SEND/MPI_IRECV (Southern Edge)',       &
     &                      MyRank, Serror, string(1:Lstr)
          exit_flag=2
          RETURN
        END IF
# endif
        icS=0
        DO m=GrecvS,1,-1
          j=Jstr-m
          icS=icS+1
          A(j)=recvS(icS)
        END DO
        IF (PRESENT(B)) THEN
          DO m=GrecvS,1,-1
            j=Jstr-m
            icS=icS+1
            B(j)=recvS(icS)
          END DO
        END IF
        IF (PRESENT(C)) THEN
          DO m=GrecvS,1,-1
            j=Jstr-m
            icS=icS+1
            C(j)=recvS(icS)
          END DO
        END IF
        IF (PRESENT(D)) THEN
          DO m=GrecvS,1,-1
            j=Jstr-m
            icS=icS+1
            D(j)=recvS(icS)
          END DO
        END IF
      END IF
!
      IF (Nexchange) THEN
# ifdef MPI
        CALL mpi_wait (Nrequest, status(1,4), Nerror)
        IF (Nerror.ne.MPI_SUCCESS) THEN
          CALL mpi_error_string (Nerror, string, Lstr, Ierror)
          Lstr=LEN_TRIM(string)
          WRITE (stdout,20) 'MPI_SEND/MPI_IRECV (Northern Edge)',       &
     &                      MyRank, Nerror, string(1:Lstr)
          exit_flag=2
          RETURN
        END IF
# endif
        icN=0
        DO m=1,GrecvN
          j=Jend+m
          icN=icN+1
          A(j)=recvN(icN)
        END DO
        IF (PRESENT(B)) THEN
          DO m=1,GrecvN
            j=Jend+m
            icN=icN+1
            B(j)=recvN(icN)
          END DO
        END IF
        IF (PRESENT(C)) THEN
          DO m=1,GrecvN
            j=Jend+m
            icN=icN+1
            C(j)=recvN(icN)
          END DO
        END IF
        IF (PRESENT(D)) THEN
          DO m=1,GrecvN
            j=Jend+m
            icN=icN+1
            D(j)=recvN(icN)
          END DO
        END IF
      END IF

# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn off time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, model, 43)
# endif

      RETURN
      END SUBROUTINE mp_exchange2d_bry
!
!***********************************************************************
      SUBROUTINE mp_exchange3d (ng, tile, model, Nvar,                  &
     &                          LBi, UBi, LBj, UBj, LBk, UBk,           &
     &                          Nghost, EW_periodic, NS_periodic,       &
     &                          A, B, C, D)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
      implicit none
!
!  Imported variable declarations.
!
      logical, intent(in) :: EW_periodic, NS_periodic

      integer, intent(in) :: ng, tile, model, Nvar
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBk, UBk
      integer, intent(in) :: Nghost

# ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: A(LBi:,LBj:,LBk:)

      real(r8), intent(inout), optional :: B(LBi:,LBj:,LBk:)
      real(r8), intent(inout), optional :: C(LBi:,LBj:,LBk:)
      real(r8), intent(inout), optional :: D(LBi:,LBj:,LBk:)
# else
      real(r8), intent(inout) :: A(LBi:UBi,LBj:UBj,LBk:UBk)

      real(r8), intent(inout), optional :: B(LBi:UBi,LBj:UBj,LBk:UBk)
      real(r8), intent(inout), optional :: C(LBi:UBi,LBj:UBj,LBk:UBk)
      real(r8), intent(inout), optional :: D(LBi:UBi,LBj:UBj,LBk:UBk)
# endif
!
!  Local variable declarations.
!
      logical :: Wexchange, Sexchange, Eexchange, Nexchange

      integer :: i, ikS, ikN, ioff, Imin, Imax, Ilen, IKlen
      integer :: j, jkW, jkE, joff, Jmin, Jmax, Jlen, JKlen
      integer :: k, kc,  m, mc, Ierror, Klen, Lstr, pp
      integer :: Wtile, GsendW, GrecvW, Wtag, Werror, Wrequest
      integer :: Stile, GsendS, GrecvS, Stag, Serror, Srequest
      integer :: Etile, GsendE, GrecvE, Etag, Eerror, Erequest
      integer :: Ntile, GsendN, GrecvN, Ntag, Nerror, Nrequest
      integer :: EWsize, sizeW, sizeE
      integer :: NSsize, sizeS, sizeN

# ifdef MPI
      integer, dimension(MPI_STATUS_SIZE,4) :: status
# endif
!
      real(r8), dimension(Nvar*HaloSizeJ(ng)*                           &
     &                          (UBk-LBk+1)) :: sendW, sendE
      real(r8), dimension(Nvar*HaloSizeJ(ng)*                           &
     &                          (UBk-LBk+1)) :: recvW, recvE

      real(r8), dimension(Nvar*HaloSizeI(ng)*                           &
     &                          (UBk-LBk+1)) :: sendS, sendN
      real(r8), dimension(Nvar*HaloSizeI(ng)*                           &
     &                          (UBk-LBk+1)) :: recvS, recvN

      character (len=MPI_MAX_ERROR_STRING) :: string

# include "set_bounds.h"

# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn on time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, model, 41)
# endif
!
!-----------------------------------------------------------------------
!  Determine rank of tile neighbors and number of ghost-points to
!  exchange.
!-----------------------------------------------------------------------
!
      CALL tile_neighbors (ng, Nghost, EW_periodic, NS_periodic,        &
     &                     GrecvW, GsendW, Wtile, Wexchange,            &
     &                     GrecvE, GsendE, Etile, Eexchange,            &
     &                     GrecvS, GsendS, Stile, Sexchange,            &
     &                     GrecvN, GsendN, Ntile, Nexchange)
!
!  Set communication tags.
!
      Wtag=1
      Stag=2
      Etag=3
      Ntag=4
!
!  Determine range and length of the distributed tile boundary segments.
!
      Imin=LBi
      Imax=UBi
      Jmin=LBj
      Jmax=UBj
      Ilen=Imax-Imin+1
      Jlen=Jmax-Jmin+1
      Klen=UBk-LBk+1
      IKlen=Ilen*Klen
      JKlen=Jlen*Klen
      IF (EW_periodic.or.NS_periodic) THEN
        pp=1
      ELSE
        pp=0
      END IF
      EWsize=Nvar*(Nghost+pp)*JKlen
      NSsize=Nvar*(Nghost+pp)*IKlen
      IF (SIZE(sendE).lt.EWsize) THEN
        WRITE (stdout,10) 'EWsize = ', EWsize, SIZE(sendE)
 10     FORMAT (/,' MP_EXCHANGE3D - communication buffer too small, ',  &
     &          a, 2i8)
      END IF
      IF (SIZE(sendN).lt.NSsize) THEN
        WRITE (stdout,10) 'NSsize = ', NSsize, SIZE(sendN)
      END IF
!
!-----------------------------------------------------------------------
!  Pack Western and Eastern tile boundary data including ghost-points.
!-----------------------------------------------------------------------
!
      IF (Wexchange) THEN
        sizeW=0
        DO m=1,GsendW
          mc=(m-1)*JKlen
          i=Istr+m-1
          DO k=LBk,UBk
            kc=(k-LBk)*Jlen+mc
            DO j=Jmin,Jmax
              sizeW=sizeW+1
              jkW=1+(j-Jmin)+kc
              sendW(jkW)=A(i,j,k)
            END DO
          END DO
        END DO
        IF (PRESENT(B)) THEN
          joff=jkW
          DO m=1,GsendW
            mc=(m-1)*JKlen
            i=Istr+m-1
            DO k=LBk,UBk
              kc=(k-LBk)*Jlen+mc
              DO j=Jmin,Jmax
                sizeW=sizeW+1
                jkW=joff+1+(j-Jmin)+kc
                sendW(jkW)=B(i,j,k)
              END DO
            END DO
          END DO
        END IF
        IF (PRESENT(C)) THEN
          joff=jkW
          DO m=1,GsendW
            mc=(m-1)*JKlen
            i=Istr+m-1
            DO k=LBk,UBk
              kc=(k-LBk)*Jlen+mc
              DO j=Jmin,Jmax
                sizeW=sizeW+1
                jkW=joff+1+(j-Jmin)+kc
                sendW(jkW)=C(i,j,k)
              END DO
            END DO
          END DO
        END IF
        IF (PRESENT(D)) THEN
          joff=jkW
          DO m=1,GsendW
            mc=(m-1)*JKlen
            i=Istr+m-1
            DO k=LBk,UBk
              kc=(k-LBk)*Jlen+mc
              DO j=Jmin,Jmax
                sizeW=sizeW+1
                jkW=joff+1+(j-Jmin)+kc
                sendW(jkW)=D(i,j,k)
              END DO
            END DO
          END DO
        END IF
      END IF
!
      IF (Eexchange) THEN
        sizeE=0
        DO m=1,GsendE
          mc=(m-1)*JKlen
          i=Iend-GsendE+m
          DO k=LBk,UBk
            kc=(k-LBk)*Jlen+mc
            DO j=Jmin,Jmax
              sizeE=sizeE+1
              jkE=1+(j-Jmin)+kc
              sendE(jkE)=A(i,j,k)
            END DO
          END DO
        END DO
        IF (PRESENT(B)) THEN
          joff=jkE
          DO m=1,GsendE
            mc=(m-1)*JKlen
            i=Iend-GsendE+m
            DO k=LBk,UBk
              kc=(k-LBk)*Jlen+mc
              DO j=Jmin,Jmax
                sizeE=sizeE+1
                jkE=joff+1+(j-Jmin)+kc
                sendE(jkE)=B(i,j,k)
              END DO
            END DO
          END DO
        END IF
        IF (PRESENT(C)) THEN
          joff=jkE
          DO m=1,GsendE
            mc=(m-1)*JKlen
            i=Iend-GsendE+m
            DO k=LBk,UBk
              kc=(k-LBk)*Jlen+mc
              DO j=Jmin,Jmax
                sizeE=sizeE+1
                jkE=joff+1+(j-Jmin)+kc
                sendE(jkE)=C(i,j,k)
              END DO
            END DO
          END DO
        END IF
        IF (PRESENT(D)) THEN
          joff=jkE
          DO m=1,GsendE
            mc=(m-1)*JKlen
            i=Iend-GsendE+m
            DO k=LBk,UBk
              kc=(k-LBk)*Jlen+mc
              DO j=Jmin,Jmax
                sizeE=sizeE+1
                jkE=joff+1+(j-Jmin)+kc
                sendE(jkE)=D(i,j,k)
              END DO
            END DO
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Send and receive Western and Eastern segments.
!-----------------------------------------------------------------------
!
# if defined MPI
      IF (Wexchange) THEN
        CALL mpi_irecv (recvW, EWsize, MP_FLOAT, Wtile, Etag,           &
     &                  OCN_COMM_WORLD, Wrequest, Werror)
      END IF
      IF (Eexchange) THEN
        CALL mpi_irecv (recvE, EWsize, MP_FLOAT, Etile, Wtag,           &
     &                  OCN_COMM_WORLD, Erequest, Eerror)
      END IF
      IF (Wexchange) THEN
        CALL mpi_send  (sendW, sizeW, MP_FLOAT, Wtile, Wtag,            &
     &                  OCN_COMM_WORLD, Werror)
      END IF
      IF (Eexchange) THEN
        CALL mpi_send  (sendE, sizeE, MP_FLOAT, Etile, Etag,            &
     &                  OCN_COMM_WORLD, Eerror)
      END IF
# endif
!
!-----------------------------------------------------------------------
!  Unpack Eastern and Western segments.
!-----------------------------------------------------------------------
!
      IF (Wexchange) THEN
# ifdef MPI
        CALL mpi_wait (Wrequest, status(1,1), Werror)
        IF (Werror.ne.MPI_SUCCESS) THEN
          CALL mpi_error_string (Werror, string, Lstr, Ierror)
          Lstr=LEN_TRIM(string)
          WRITE (stdout,20) 'MPI_SEND/MPI_IRECV (Western Edge)',        &
     &                      MyRank, Werror, string(1:Lstr)
          exit_flag=2
          RETURN
        END IF
# endif
        DO m=GrecvW,1,-1
          mc=(GrecvW-m)*JKlen
          i=Istr-m
          DO k=LBk,UBk
            kc=(k-LBk)*Jlen+mc
            DO j=Jmin,Jmax
              jkW=1+(j-Jmin)+kc
              A(i,j,k)=recvW(jkW)
            END DO
          END DO
        END DO
        IF (PRESENT(B)) THEN
          joff=jkW
          DO m=GrecvW,1,-1
            mc=(GrecvW-m)*JKlen
            i=Istr-m
            DO k=LBk,UBk
              kc=(k-LBk)*Jlen+mc
              DO j=Jmin,Jmax
                jkW=joff+1+(j-Jmin)+kc
                B(i,j,k)=recvW(jkW)
              END DO
            END DO
          END DO
        END IF
        IF (PRESENT(C)) THEN
          joff=jkW
          DO m=GrecvW,1,-1
            mc=(GrecvW-m)*JKlen
            i=Istr-m
            DO k=LBk,UBk
              kc=(k-LBk)*Jlen+mc
              DO j=Jmin,Jmax
                jkW=joff+1+(j-Jmin)+kc
                C(i,j,k)=recvW(jkW)
              END DO
            END DO
          END DO
        END IF
        IF (PRESENT(D)) THEN
          joff=jkW
          DO m=GrecvW,1,-1
            mc=(GrecvW-m)*JKlen
            i=Istr-m
            DO k=LBk,UBk
              kc=(k-LBk)*Jlen+mc
              DO j=Jmin,Jmax
                jkW=joff+1+(j-Jmin)+kc
                D(i,j,k)=recvW(jkW)
              END DO
            END DO
          END DO
        END IF
      END IF
!
      IF (Eexchange) THEN
# ifdef MPI
        CALL mpi_wait (Erequest, status(1,3), Eerror)
        IF (Eerror.ne.MPI_SUCCESS) THEN
          CALL mpi_error_string (Eerror, string, Lstr, Ierror)
          Lstr=LEN_TRIM(string)
          WRITE (stdout,20) 'MPI_SEND/MPI_IRECV (Eastern Edge)',        &
     &                      MyRank, Eerror, string(1:Lstr)
 20       FORMAT (/,' MP_EXCHANGE3D - error during ',a,                 &
     &            ' call, Node = ',i3.3,' Error = ',i3,/,15x,a)
          exit_flag=2
          RETURN
        END IF
# endif
        DO m=1,GrecvE
          mc=(m-1)*JKlen
          i=Iend+m
          DO k=LBk,UBk
            kc=(k-LBk)*Jlen+mc
            DO j=Jmin,Jmax
              jkE=1+(j-Jmin)+kc
              A(i,j,k)=recvE(jkE)
            END DO
          ENDDO
        END DO
        IF (PRESENT(B)) THEN
          joff=jkE
          DO m=1,GrecvE
            mc=(m-1)*JKlen
            i=Iend+m
            DO k=LBk,UBk
              kc=(k-LBk)*Jlen+mc
              DO j=Jmin,Jmax
                jkE=joff+1+(j-Jmin)+kc
                B(i,j,k)=recvE(jkE)
              END DO
            ENDDO
          END DO
        END IF
        IF (PRESENT(C)) THEN
          joff=jkE
          DO m=1,GrecvE
            mc=(m-1)*JKlen
            i=Iend+m
            DO k=LBk,UBk
              kc=(k-LBk)*Jlen+mc
              DO j=Jmin,Jmax
                jkE=joff+1+(j-Jmin)+kc
                C(i,j,k)=recvE(jkE)
              END DO
            ENDDO
          END DO
        END IF
        IF (PRESENT(D)) THEN
          joff=jkE
          DO m=1,GrecvE
            mc=(m-1)*JKlen
            i=Iend+m
            DO k=LBk,UBk
              kc=(k-LBk)*Jlen+mc
              DO j=Jmin,Jmax
                jkE=joff+1+(j-Jmin)+kc
                D(i,j,k)=recvE(jkE)
              END DO
            ENDDO
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Pack Southern and Northern tile boundary data including ghost-points.
!-----------------------------------------------------------------------
!
      IF (Sexchange) THEN
        sizeS=0
        DO m=1,GsendS
          mc=(m-1)*IKlen
          j=Jstr+m-1
          DO k=LBk,UBk
            kc=(k-LBk)*Ilen+mc
            DO i=Imin,Imax
              sizeS=sizeS+1
              ikS=1+(i-Imin)+kc
              sendS(ikS)=A(i,j,k)
            END DO
          END DO
        END DO
        IF (PRESENT(B)) THEN
          ioff=ikS
          DO m=1,GsendS
            mc=(m-1)*IKlen
            j=Jstr+m-1
            DO k=LBk,UBk
              kc=(k-LBk)*Ilen+mc
              DO i=Imin,Imax
                sizeS=sizeS+1
                ikS=ioff+1+(i-Imin)+kc
                sendS(ikS)=B(i,j,k)
              END DO
            END DO
          END DO
        END IF
        IF (PRESENT(C)) THEN
          ioff=ikS
          DO m=1,GsendS
            mc=(m-1)*IKlen
            j=Jstr+m-1
            DO k=LBk,UBk
              kc=(k-LBk)*Ilen+mc
              DO i=Imin,Imax
                sizeS=sizeS+1
                ikS=ioff+1+(i-Imin)+kc
                sendS(ikS)=C(i,j,k)
              END DO
            END DO
          END DO
        END IF
        IF (PRESENT(D)) THEN
          ioff=ikS
          DO m=1,GsendS
            mc=(m-1)*IKlen
            j=Jstr+m-1
            DO k=LBk,UBk
              kc=(k-LBk)*Ilen+mc
              DO i=Imin,Imax
                sizeS=sizeS+1
                ikS=ioff+1+(i-Imin)+kc
                sendS(ikS)=D(i,j,k)
              END DO
            END DO
          END DO
        END IF
      END IF
!
      IF (Nexchange) THEN
        sizeN=0
        DO m=1,GsendN
          mc=(m-1)*IKlen
          j=Jend-GsendN+m
          DO k=LBk,UBk
            kc=(k-LBk)*Ilen+mc
            DO i=Imin,Imax
              sizeN=sizeN+1
              ikN=1+(i-Imin)+kc
              sendN(ikN)=A(i,j,k)
            END DO
          END DO
        END DO
        IF (PRESENT(B)) THEN
          ioff=ikN
          DO m=1,GsendN
            mc=(m-1)*IKlen
            j=Jend-GsendN+m
            DO k=LBk,UBk
              kc=(k-LBk)*Ilen+mc
              DO i=Imin,Imax
                sizeN=sizeN+1
                ikN=ioff+1+(i-Imin)+kc
                sendN(ikN)=B(i,j,k)
              END DO
            END DO
          END DO
        END IF
        IF (PRESENT(C)) THEN
          ioff=ikN
          DO m=1,GsendN
            mc=(m-1)*IKlen
            j=Jend-GsendN+m
            DO k=LBk,UBk
              kc=(k-LBk)*Ilen+mc
              DO i=Imin,Imax
                sizeN=sizeN+1
                ikN=ioff+1+(i-Imin)+kc
                sendN(ikN)=C(i,j,k)
              END DO
            END DO
          END DO
        END IF
        IF (PRESENT(D)) THEN
          ioff=ikN
          DO m=1,GsendN
            mc=(m-1)*IKlen
            j=Jend-GsendN+m
            DO k=LBk,UBk
              kc=(k-LBk)*Ilen+mc
              DO i=Imin,Imax
                sizeN=sizeN+1
                ikN=ioff+1+(i-Imin)+kc
                sendN(ikN)=D(i,j,k)
              END DO
            END DO
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Send and receive Southern and Northern segments.
!-----------------------------------------------------------------------
!
# if defined MPI
      IF (Sexchange) THEN
        CALL mpi_irecv (recvS, NSsize, MP_FLOAT, Stile, Ntag,           &
     &                  OCN_COMM_WORLD, Srequest, Serror)
      END IF
      IF (Nexchange) THEN
        CALL mpi_irecv (recvN, NSsize, MP_FLOAT, Ntile, Stag,           &
     &                  OCN_COMM_WORLD, Nrequest, Nerror)
      END IF
      IF (Sexchange) THEN
        CALL mpi_send  (sendS, sizeS, MP_FLOAT, Stile, Stag,            &
     &                  OCN_COMM_WORLD, Serror)
      END IF
      IF (Nexchange) THEN
        CALL mpi_send  (sendN, sizeN, MP_FLOAT, Ntile, Ntag,            &
     &                  OCN_COMM_WORLD, Nerror)
      END IF
# endif
!
!-----------------------------------------------------------------------
!  Unpack Northern and Southern segments.
!-----------------------------------------------------------------------
!
      IF (Sexchange) THEN
# ifdef MPI
        CALL mpi_wait (Srequest, status(1,2), Serror)
        IF (Serror.ne.MPI_SUCCESS) THEN
          CALL mpi_error_string (Serror, string, Lstr, Ierror)
          Lstr=LEN_TRIM(string)
          WRITE (stdout,20) 'MPI_SEND/MPI_IRECV (Southern Edge)',       &
     &                      MyRank, Serror, string(1:Lstr)
          exit_flag=2
          RETURN
        END IF
# endif
        DO m=GrecvS,1,-1
          mc=(GrecvS-m)*IKlen
          j=Jstr-m
          DO k=LBk,UBk
            kc=(k-LBk)*Ilen+mc
            DO i=Imin,Imax
              ikS=1+(i-Imin)+kc
              A(i,j,k)=recvS(ikS)
            END DO
          END DO
        END DO
        IF (PRESENT(B)) THEN
          ioff=ikS
          DO m=GrecvS,1,-1
            mc=(GrecvS-m)*IKlen
            j=Jstr-m
            DO k=LBk,UBk
              kc=(k-LBk)*Ilen+mc
              DO i=Imin,Imax
                ikS=ioff+1+(i-Imin)+kc
                B(i,j,k)=recvS(ikS)
              END DO
            END DO
          END DO
        END IF
        IF (PRESENT(C)) THEN
          ioff=ikS
          DO m=GrecvS,1,-1
            mc=(GrecvS-m)*IKlen
            j=Jstr-m
            DO k=LBk,UBk
              kc=(k-LBk)*Ilen+mc
              DO i=Imin,Imax
                ikS=ioff+1+(i-Imin)+kc
                C(i,j,k)=recvS(ikS)
              END DO
            END DO
          END DO
        END IF
        IF (PRESENT(D)) THEN
          ioff=ikS
          DO m=GrecvS,1,-1
            mc=(GrecvS-m)*IKlen
            j=Jstr-m
            DO k=LBk,UBk
              kc=(k-LBk)*Ilen+mc
              DO i=Imin,Imax
                ikS=ioff+1+(i-Imin)+kc
                D(i,j,k)=recvS(ikS)
              END DO
            END DO
          END DO
        END IF
      END IF
!
      IF (Nexchange) THEN
# ifdef MPI
        CALL mpi_wait (Nrequest, status(1,4), Nerror)
        IF (Nerror.ne.MPI_SUCCESS) THEN
          CALL mpi_error_string (Nerror, string, Lstr, Ierror)
          Lstr=LEN_TRIM(string)
          WRITE (stdout,20) 'MPI_SEND/MPI_IRECV (Northern Edge)',       &
     &                      MyRank, Nerror, string(1:Lstr)
          exit_flag=2
          RETURN
        END IF
# endif
        DO m=1,GrecvN
          mc=(m-1)*IKlen
          j=Jend+m
          DO k=LBk,UBk
            kc=(k-LBk)*Ilen+mc
            DO i=Imin,Imax
              ikN=1+(i-Imin)+kc
              A(i,j,k)=recvN(ikN)
            END DO
          END DO
        END DO
        IF (PRESENT(B)) THEN
          ioff=ikN
          DO m=1,GrecvN
            mc=(m-1)*IKlen
            j=Jend+m
            DO k=LBk,UBk
              kc=(k-LBk)*Ilen+mc
              DO i=Imin,Imax
                ikN=ioff+1+(i-Imin)+kc
                B(i,j,k)=recvN(ikN)
              END DO
            END DO
          END DO
        END IF
        IF (PRESENT(C)) THEN
          ioff=ikN
          DO m=1,GrecvN
            mc=(m-1)*IKlen
            j=Jend+m
            DO k=LBk,UBk
              kc=(k-LBk)*Ilen+mc
              DO i=Imin,Imax
                ikN=ioff+1+(i-Imin)+kc
                C(i,j,k)=recvN(ikN)
              END DO
            END DO
          END DO
        END IF
        IF (PRESENT(D)) THEN
          ioff=ikN
          DO m=1,GrecvN
            mc=(m-1)*IKlen
            j=Jend+m
            DO k=LBk,UBk
              kc=(k-LBk)*Ilen+mc
              DO i=Imin,Imax
                ikN=ioff+1+(i-Imin)+kc
                D(i,j,k)=recvN(ikN)
              END DO
            END DO
          END DO
        END IF
      END IF
# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn off time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, model, 41)
# endif

      RETURN
      END SUBROUTINE mp_exchange3d

!
!***********************************************************************
      SUBROUTINE mp_exchange3d_bry (ng, tile, model, Nvar, boundary,    &
     &                              LBij, UBij, LBk, UBk,               &
     &                              Nghost, EW_periodic, NS_periodic,   &
     &                              A, B, C, D)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
      implicit none
!
!  Imported variable declarations.
!
      logical, intent(in) :: EW_periodic, NS_periodic

      integer, intent(in) :: ng, tile, model, Nvar, boundary
      integer, intent(in) :: LBij, UBij, LBk, UBk
      integer, intent(in) :: Nghost
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: A(LBij:,LBk:)

      real(r8), intent(inout), optional :: B(LBij:,LBk:)
      real(r8), intent(inout), optional :: C(LBij:,LBk:)
      real(r8), intent(inout), optional :: D(LBij:,LBk:)
# else
      real(r8), intent(inout) :: A(LBij:UBij,LBk:UBk)

      real(r8), intent(inout), optional :: B(LBij:UBij,LBk:UBk)
      real(r8), intent(inout), optional :: C(LBij:UBij,LBk:UBk)
      real(r8), intent(inout), optional :: D(LBij:UBij,LBk:UBk)
# endif
!
!  Local variable declarations.
!
      logical :: Wexchange, Sexchange, Eexchange, Nexchange

      integer :: i, ikS, ikN, ioff
      integer :: j, jkW, jkE, joff
      integer :: k, m, mc, Ierror, Klen, Lstr, pp
      integer :: Wtile, GsendW, GrecvW, Wtag, Werror, Wrequest
      integer :: Stile, GsendS, GrecvS, Stag, Serror, Srequest
      integer :: Etile, GsendE, GrecvE, Etag, Eerror, Erequest
      integer :: Ntile, GsendN, GrecvN, Ntag, Nerror, Nrequest
      integer :: EWsize, sizeW, sizeE
      integer :: NSsize, sizeS, sizeN

# ifdef MPI
      integer, dimension(MPI_STATUS_SIZE,4) :: status
# endif
!
      real(r8), dimension(Nvar*HaloBry(ng)*(UBk-LBk+1)) :: sendW, sendE
      real(r8), dimension(Nvar*HaloBry(ng)*(UBk-LBk+1)) :: recvW, recvE
      real(r8), dimension(Nvar*HaloBry(ng)*(UBk-LBk+1)) :: sendS, sendN
      real(r8), dimension(Nvar*HaloBry(ng)*(UBk-LBk+1)) :: recvS, recvN

      character (len=MPI_MAX_ERROR_STRING) :: string

# include "set_bounds.h"

# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn on time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, model, 43)
# endif
!
!-----------------------------------------------------------------------
!  Determine rank of tile neighbors and number of ghost-points to
!  exchange.
!-----------------------------------------------------------------------
!
      CALL tile_neighbors (ng, Nghost, EW_periodic, NS_periodic,        &
     &                     GrecvW, GsendW, Wtile, Wexchange,            &
     &                     GrecvE, GsendE, Etile, Eexchange,            &
     &                     GrecvS, GsendS, Stile, Sexchange,            &
     &                     GrecvN, GsendN, Ntile, Nexchange)
!
!  Adjust exchange swiches according to boundary edge to process.
!
      Wexchange=Wexchange.and.((boundary.eq.isouth).or.                 &
     &                         (boundary.eq.inorth))
      Eexchange=Eexchange.and.((boundary.eq.isouth).or.                 &
     &                         (boundary.eq.inorth))
      Sexchange=Sexchange.and.((boundary.eq.iwest).or.                  &
     &                         (boundary.eq.ieast))
      Nexchange=Nexchange.and.((boundary.eq.iwest).or.                  &
     &                         (boundary.eq.ieast))
!
!  Set communication tags.
!
      Wtag=1
      Stag=2
      Etag=3
      Ntag=4
!
!  Determine range and length of the distributed tile boundary segments.
!
      Klen=UBk-LBk+1
      IF (EW_periodic.or.NS_periodic) THEN
        pp=1
      ELSE
        pp=0
      END IF
      EWsize=Nvar*(Nghost+pp)*Klen
      NSsize=Nvar*(Nghost+pp)*Klen
      IF (SIZE(sendE).lt.EWsize) THEN
        WRITE (stdout,10) 'EWsize = ', EWsize, SIZE(sendE)
 10     FORMAT (/,' MP_EXCHANGE3D_BRY - communication buffer too ',     &
     &          'small, ', a, 2i8)
      END IF
      IF (SIZE(sendN).lt.NSsize) THEN
        WRITE (stdout,10) 'NSsize = ', NSsize, SIZE(sendN)
      END IF
!
!-----------------------------------------------------------------------
!  Pack Western and Eastern tile boundary data including ghost-points.
!-----------------------------------------------------------------------
!
      IF (Wexchange) THEN
        sizeW=0
        DO m=1,GsendW
          mc=(m-1)*Klen
          i=Istr+m-1
          DO k=LBk,UBk
            sizeW=sizeW+1
            jkW=1+(k-LBk)+mc
            sendW(jkW)=A(i,k)
          END DO
        END DO
        IF (PRESENT(B)) THEN
          joff=jkW
          DO m=1,GsendW
            mc=(m-1)*Klen
            i=Istr+m-1
            DO k=LBk,UBk
              sizeW=sizeW+1
              jkW=joff+1+(k-LBk)+mc
              sendW(jkW)=B(i,k)
            END DO
          END DO
        END IF
        IF (PRESENT(C)) THEN
          joff=jkW
          DO m=1,GsendW
            mc=(m-1)*Klen
            i=Istr+m-1
            DO k=LBk,UBk
              sizeW=sizeW+1
              jkW=joff+1+(k-LBk)+mc
              sendW(jkW)=C(i,k)
            END DO
          END DO
        END IF
        IF (PRESENT(D)) THEN
          joff=jkW
          DO m=1,GsendW
            mc=(m-1)*Klen
            i=Istr+m-1
            DO k=LBk,UBk
              sizeW=sizeW+1
              jkW=joff+1+(k-LBk)+mc
              sendW(jkW)=D(i,k)
            END DO
          END DO
        END IF
      END IF
!
      IF (Eexchange) THEN
        sizeE=0
        DO m=1,GsendE
          mc=(m-1)*Klen
          i=Iend-GsendE+m
          DO k=LBk,UBk
            sizeE=sizeE+1
            jkE=1+(k-LBk)+mc
            sendE(jkE)=A(i,k)
          END DO
        END DO
        IF (PRESENT(B)) THEN
          joff=jkE
          DO m=1,GsendE
            mc=(m-1)*Klen
            i=Iend-GsendE+m
            DO k=LBk,UBk
              sizeE=sizeE+1
              jkE=joff+1+(k-LBk)+mc
              sendE(jkE)=B(i,k)
            END DO
          END DO
        END IF
        IF (PRESENT(C)) THEN
          joff=jkE
          DO m=1,GsendE
            mc=(m-1)*Klen
            i=Iend-GsendE+m
            DO k=LBk,UBk
              sizeE=sizeE+1
              jkE=joff+1+(k-LBk)+mc
              sendE(jkE)=C(i,k)
            END DO
          END DO
        END IF
        IF (PRESENT(D)) THEN
          joff=jkE
          DO m=1,GsendE
            mc=(m-1)*Klen
            i=Iend-GsendE+m
            DO k=LBk,UBk
              sizeE=sizeE+1
              jkE=joff+1+(k-LBk)+mc
              sendE(jkE)=D(i,k)
            END DO
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Send and receive Western and Eastern segments.
!-----------------------------------------------------------------------
!
# if defined MPI
      IF (Wexchange) THEN
        CALL mpi_irecv (recvW, EWsize, MP_FLOAT, Wtile, Etag,           &
     &                  OCN_COMM_WORLD, Wrequest, Werror)
      END IF
      IF (Eexchange) THEN
        CALL mpi_irecv (recvE, EWsize, MP_FLOAT, Etile, Wtag,           &
     &                  OCN_COMM_WORLD, Erequest, Eerror)
      END IF
      IF (Wexchange) THEN
        CALL mpi_send  (sendW, sizeW, MP_FLOAT, Wtile, Wtag,            &
     &                  OCN_COMM_WORLD, Werror)
      END IF
      IF (Eexchange) THEN
        CALL mpi_send  (sendE, sizeE, MP_FLOAT, Etile, Etag,            &
     &                  OCN_COMM_WORLD, Eerror)
      END IF
# endif
!
!-----------------------------------------------------------------------
!  Unpack Eastern and Western segments.
!-----------------------------------------------------------------------
!
      IF (Wexchange) THEN
# ifdef MPI
        CALL mpi_wait (Wrequest, status(1,1), Werror)
        IF (Werror.ne.MPI_SUCCESS) THEN
          CALL mpi_error_string (Werror, string, Lstr, Ierror)
          Lstr=LEN_TRIM(string)
          WRITE (stdout,20) 'MPI_SEND/MPI_IRECV (Western Edge)',        &
     &                      MyRank, Werror, string(1:Lstr)
          exit_flag=2
          RETURN
        END IF
# endif
        DO m=GrecvW,1,-1
          mc=(GrecvW-m)*Klen
          i=Istr-m
          DO k=LBk,UBk
            jkW=1+(k-LBk)+mc
            A(i,k)=recvW(jkW)
          END DO
        END DO
        IF (PRESENT(B)) THEN
          joff=jkW
          DO m=GrecvW,1,-1
            mc=(GrecvW-m)*Klen
            i=Istr-m
            DO k=LBk,UBk
              jkW=joff+1+(k-LBk)+mc
              B(i,k)=recvW(jkW)
            END DO
          END DO
        END IF
        IF (PRESENT(C)) THEN
          joff=jkW
          DO m=GrecvW,1,-1
            mc=(GrecvW-m)*Klen
            i=Istr-m
            DO k=LBk,UBk
              jkW=joff+1+(k-LBk)+mc
              C(i,k)=recvW(jkW)
            END DO
          END DO
        END IF
        IF (PRESENT(D)) THEN
          joff=jkW
          DO m=GrecvW,1,-1
            mc=(GrecvW-m)*Klen
            i=Istr-m
            DO k=LBk,UBk
              jkW=joff+1+(k-LBk)+mc
              D(i,k)=recvW(jkW)
            END DO
          END DO
        END IF
      END IF
!
      IF (Eexchange) THEN
# ifdef MPI
        CALL mpi_wait (Erequest, status(1,3), Eerror)
        IF (Eerror.ne.MPI_SUCCESS) THEN
          CALL mpi_error_string (Eerror, string, Lstr, Ierror)
          Lstr=LEN_TRIM(string)
          WRITE (stdout,20) 'MPI_SEND/MPI_IRECV (Eastern Edge)',        &
     &                      MyRank, Eerror, string(1:Lstr)
 20       FORMAT (/,' MP_EXCHANGE3D_BRY - error during ',a,             &
     &            ' call, Node = ',i3.3,' Error = ',i3,/,15x,a)
          exit_flag=2
          RETURN
        END IF
# endif
        DO m=1,GrecvE
          mc=(m-1)*Klen
          i=Iend+m
          DO k=LBk,UBk
            jkE=1+(k-LBk)+mc
            A(i,k)=recvE(jkE)
          END DO
        END DO
        IF (PRESENT(B)) THEN
          joff=jkE
          DO m=1,GrecvE
            mc=(m-1)*Klen
            i=Iend+m
            DO k=LBk,UBk
              jkE=joff+1+(k-LBk)+mc
              B(i,k)=recvE(jkE)
            END DO
          END DO
        END IF
        IF (PRESENT(C)) THEN
          joff=jkE
          DO m=1,GrecvE
            mc=(m-1)*Klen
            i=Iend+m
            DO k=LBk,UBk
              jkE=joff+1+(k-LBk)+mc
              C(i,k)=recvE(jkE)
            END DO
          END DO
        END IF
        IF (PRESENT(D)) THEN
          joff=jkE
          DO m=1,GrecvE
            mc=(m-1)*Klen
            i=Iend+m
            DO k=LBk,UBk
              jkE=joff+1+(k-LBk)+mc
              D(i,k)=recvE(jkE)
            END DO
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Pack Southern and Northern tile boundary data including ghost-points.
!-----------------------------------------------------------------------
!
      IF (Sexchange) THEN
        sizeS=0
        DO m=1,GsendS
          mc=(m-1)*Klen
          j=Jstr+m-1
          DO k=LBk,UBk
            sizeS=sizeS+1
            ikS=1+(k-LBk)+mc
            sendS(ikS)=A(j,k)
          END DO
        END DO
        IF (PRESENT(B)) THEN
          ioff=ikS
          DO m=1,GsendS
            mc=(m-1)*Klen
            j=Jstr+m-1
            DO k=LBk,UBk
              sizeS=sizeS+1
              ikS=ioff+1+(k-LBk)+mc
              sendS(ikS)=B(j,k)
            END DO
          END DO
        END IF
        IF (PRESENT(C)) THEN
          ioff=ikS
          DO m=1,GsendS
            mc=(m-1)*Klen
            j=Jstr+m-1
            DO k=LBk,UBk
              sizeS=sizeS+1
              ikS=ioff+1+(k-LBk)+mc
              sendS(ikS)=C(j,k)
            END DO
          END DO
        END IF
        IF (PRESENT(D)) THEN
          ioff=ikS
          DO m=1,GsendS
            mc=(m-1)*Klen
            j=Jstr+m-1
            DO k=LBk,UBk
              sizeS=sizeS+1
              ikS=ioff+1+(k-LBk)+mc
              sendS(ikS)=D(j,k)
            END DO
          END DO
        END IF
      END IF
!
      IF (Nexchange) THEN
        sizeN=0
        DO m=1,GsendN
          mc=(m-1)*Klen
          j=Jend-GsendN+m
          DO k=LBk,UBk
            sizeN=sizeN+1
            ikN=1+(k-LBk)+mc
            sendN(ikN)=A(j,k)
          END DO
        END DO
        IF (PRESENT(B)) THEN
          ioff=ikN
          DO m=1,GsendN
            mc=(m-1)*Klen
            j=Jend-GsendN+m
            DO k=LBk,UBk
              sizeN=sizeN+1
              ikN=ioff+1+(k-LBk)+mc
              sendN(ikN)=B(j,k)
            END DO
          END DO
        END IF
        IF (PRESENT(C)) THEN
          ioff=ikN
          DO m=1,GsendN
            mc=(m-1)*Klen
            j=Jend-GsendN+m
            DO k=LBk,UBk
              sizeN=sizeN+1
              ikN=ioff+1+(k-LBk)+mc
              sendN(ikN)=C(j,k)
            END DO
          END DO
        END IF
        IF (PRESENT(D)) THEN
          ioff=ikN
          DO m=1,GsendN
            mc=(m-1)*Klen
            j=Jend-GsendN+m
            DO k=LBk,UBk
              sizeN=sizeN+1
              ikN=ioff+1+(k-LBk)+mc
              sendN(ikN)=D(j,k)
            END DO
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Send and receive Southern and Northern segments.
!-----------------------------------------------------------------------
!
# if defined MPI
      IF (Sexchange) THEN
        CALL mpi_irecv (recvS, NSsize, MP_FLOAT, Stile, Ntag,           &
     &                  OCN_COMM_WORLD, Srequest, Serror)
      END IF
      IF (Nexchange) THEN
        CALL mpi_irecv (recvN, NSsize, MP_FLOAT, Ntile, Stag,           &
     &                  OCN_COMM_WORLD, Nrequest, Nerror)
      END IF
      IF (Sexchange) THEN
        CALL mpi_send  (sendS, sizeS, MP_FLOAT, Stile, Stag,            &
     &                  OCN_COMM_WORLD, Serror)
      END IF
      IF (Nexchange) THEN
        CALL mpi_send  (sendN, sizeN, MP_FLOAT, Ntile, Ntag,            &
     &                  OCN_COMM_WORLD, Nerror)
      END IF
# endif
!
!-----------------------------------------------------------------------
!  Unpack Northern and Southern segments.
!-----------------------------------------------------------------------
!
      IF (Sexchange) THEN
# ifdef MPI
        CALL mpi_wait (Srequest, status(1,2), Serror)
        IF (Serror.ne.MPI_SUCCESS) THEN
          CALL mpi_error_string (Serror, string, Lstr, Ierror)
          Lstr=LEN_TRIM(string)
          WRITE (stdout,20) 'MPI_SEND/MPI_IRECV (Southern Edge)',       &
     &                      MyRank, Serror, string(1:Lstr)
          exit_flag=2
          RETURN
        END IF
# endif
        DO m=GrecvS,1,-1
          mc=(GrecvS-m)*Klen
          j=Jstr-m
          DO k=LBk,UBk
            ikS=1+(k-LBk)+mc
            A(j,k)=recvS(ikS)
          END DO
        END DO
        IF (PRESENT(B)) THEN
          ioff=ikS
          DO m=GrecvS,1,-1
            mc=(GrecvS-m)*Klen
            j=Jstr-m
            DO k=LBk,UBk
              ikS=ioff+1+(k-LBk)+mc
              B(j,k)=recvS(ikS)
            END DO
          END DO
        END IF
        IF (PRESENT(C)) THEN
          ioff=ikS
          DO m=GrecvS,1,-1
            mc=(GrecvS-m)*Klen
            j=Jstr-m
            DO k=LBk,UBk
              ikS=ioff+1+(k-LBk)+mc
              C(j,k)=recvS(ikS)
            END DO
          END DO
        END IF
        IF (PRESENT(D)) THEN
          ioff=ikS
          DO m=GrecvS,1,-1
            mc=(GrecvS-m)*Klen
            j=Jstr-m
            DO k=LBk,UBk
              ikS=ioff+1+(k-LBk)+mc
              D(j,k)=recvS(ikS)
            END DO
          END DO
        END IF
      END IF
!
      IF (Nexchange) THEN
# ifdef MPI
        CALL mpi_wait (Nrequest, status(1,4), Nerror)
        IF (Nerror.ne.MPI_SUCCESS) THEN
          CALL mpi_error_string (Nerror, string, Lstr, Ierror)
          Lstr=LEN_TRIM(string)
          WRITE (stdout,20) 'MPI_SEND/MPI_IRECV (Northern Edge)',       &
     &                      MyRank, Nerror, string(1:Lstr)
          exit_flag=2
          RETURN
        END IF
# endif
        DO m=1,GrecvN
          mc=(m-1)*Klen
          j=Jend+m
          DO k=LBk,UBk
            ikN=1+(k-LBk)+mc
            A(j,k)=recvN(ikN)
          END DO
        END DO
        IF (PRESENT(B)) THEN
          ioff=ikN
          DO m=1,GrecvN
            mc=(m-1)*Klen
            j=Jend+m
            DO k=LBk,UBk
              ikN=ioff+1+(k-LBk)+mc
              B(j,k)=recvN(ikN)
            END DO
          END DO
        END IF
        IF (PRESENT(C)) THEN
          ioff=ikN
          DO m=1,GrecvN
            mc=(m-1)*Klen
            j=Jend+m
            DO k=LBk,UBk
              ikN=ioff+1+(k-LBk)+mc
              C(j,k)=recvN(ikN)
            END DO
          END DO
        END IF
        IF (PRESENT(D)) THEN
          ioff=ikN
          DO m=1,GrecvN
            mc=(m-1)*Klen
            j=Jend+m
            DO k=LBk,UBk
              ikN=ioff+1+(k-LBk)+mc
              D(j,k)=recvN(ikN)
            END DO
          END DO
        END IF
      END IF

# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn off time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, model, 43)
# endif

      RETURN
      END SUBROUTINE mp_exchange3d_bry

!
!***********************************************************************
      SUBROUTINE mp_exchange4d (ng, tile, model, Nvar,                  &
     &                          LBi, UBi, LBj, UBj, LBk, UBk, LBt, UBt, &
     &                          Nghost, EW_periodic, NS_periodic,       &
     &                          A, B, C)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
      implicit none
!
!  Imported variable declarations.
!
      logical, intent(in) :: EW_periodic, NS_periodic

      integer, intent(in) :: ng, tile, model, Nvar
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBk, UBk, LBt, UBt
      integer, intent(in) :: Nghost
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: A(LBi:,LBj:,LBk:,LBt:)

      real(r8), intent(inout), optional :: B(LBi:,LBj:,LBk:,LBt:)
      real(r8), intent(inout), optional :: C(LBi:,LBj:,LBk:,LBt:)

# else
      real(r8), intent(inout) :: A(LBi:UBi,LBj:UBj,LBk:UBk,LBt:UBt)

      real(r8), intent(inout), optional ::                              &
     &                           B(LBi:UBi,LBj:UBj,LBk:UBk,LBt:UBt)
      real(r8), intent(inout), optional ::                              &
     &                           C(LBi:UBi,LBj:UBj,LBk:UBk,LBt:UBt)
# endif
!
!  Local variable declarations.
!
      logical :: Wexchange, Sexchange, Eexchange, Nexchange

      integer :: i, ikS, ikN, ioff, Imin, Imax, Ilen, IKlen, IKTlen
      integer :: j, jkW, jkE, joff, Jmin, Jmax, Jlen, JKlen, JKTlen
      integer :: k, kc, m, mc, Ierror, Klen, Lstr, Tlen, pp
      integer :: l, lc
      integer :: Wtile, GsendW, GrecvW, Wtag, Werror, Wrequest
      integer :: Stile, GsendS, GrecvS, Stag, Serror, Srequest
      integer :: Etile, GsendE, GrecvE, Etag, Eerror, Erequest
      integer :: Ntile, GsendN, GrecvN, Ntag, Nerror, Nrequest
      integer :: EWsize, sizeW, sizeE
      integer :: NSsize, sizeS, sizeN

# ifdef MPI
      integer, dimension(MPI_STATUS_SIZE,4) :: status
# endif
!
      real(r8), dimension(Nvar*HaloSizeJ(ng)*                           &
     &                    (UBk-LBk+1)*(UBt-LBt+1)) :: sendW, sendE
      real(r8), dimension(Nvar*HaloSizeJ(ng)*                           &
     &                    (UBk-LBk+1)*(UBt-LBt+1)) :: recvW, recvE

      real(r8), dimension(Nvar*HaloSizeI(ng)*                           &
     &                    (UBk-LBk+1)*(UBt-LBt+1)) :: sendS, sendN
      real(r8), dimension(Nvar*HaloSizeI(ng)*                           &
     &                    (UBk-LBk+1)*(UBt-LBt+1)) :: recvS, recvN

      character (len=MPI_MAX_ERROR_STRING) :: string

# include "set_bounds.h"

# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn on time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, model, 42)
# endif
!
!-----------------------------------------------------------------------
!  Determine rank of tile neighbors and number of ghost-points to
!  exchange.
!-----------------------------------------------------------------------
!
      CALL tile_neighbors (ng, Nghost, EW_periodic, NS_periodic,        &
     &                     GrecvW, GsendW, Wtile, Wexchange,            &
     &                     GrecvE, GsendE, Etile, Eexchange,            &
     &                     GrecvS, GsendS, Stile, Sexchange,            &
     &                     GrecvN, GsendN, Ntile, Nexchange)
!
!  Set communication tags.
!
      Wtag=1
      Stag=2
      Etag=3
      Ntag=4
!
!  Determine range and length of the distributed tile boundary segments.
!
      Imin=LBi
      Imax=UBi
      Jmin=LBj
      Jmax=UBj
      Ilen=Imax-Imin+1
      Jlen=Jmax-Jmin+1
      Klen=UBk-LBk+1
      Tlen=UBt-LBt+1
      IKlen=Ilen*Klen
      JKlen=Jlen*Klen
      IKTlen=IKlen*Tlen
      JKTlen=JKlen*Tlen
      IF (EW_periodic.or.NS_periodic) THEN
        pp=1
      ELSE
        pp=0
      END IF
      EWsize=Nvar*(Nghost+pp)*JKTlen
      NSsize=Nvar*(Nghost+pp)*IKTlen
      IF (SIZE(sendE).lt.EWsize) THEN
        WRITE (stdout,10) 'EWsize = ', EWsize, SIZE(sendE)
 10     FORMAT (/,' MP_EXCHANGE4D - communication buffer too small, ',  &
     &          a, 2i8)
      END IF
      IF (SIZE(sendN).lt.NSsize) THEN
        WRITE (stdout,10) 'NSsize = ', NSsize, SIZE(sendN)
      END IF
!
!-----------------------------------------------------------------------
!  Pack Western and Eastern tile boundary data including ghost-points.
!-----------------------------------------------------------------------
!
      IF (Wexchange) THEN
        sizeW=0
        DO m=1,GsendW
          mc=(m-1)*JKTlen
          i=Istr+m-1
          DO l=LBt,UBt
           lc=(l-LBt)*JKlen+mc
            DO k=LBk,UBk
              kc=(k-LBk)*Jlen+lc
              DO j=Jmin,Jmax
                sizeW=sizeW+1
                jkW=1+(j-Jmin)+kc
                sendW(jkW)=A(i,j,k,l)
              END DO
            END DO
          END DO
        END DO
        IF (PRESENT(B)) THEN
          joff=jkW
          DO m=1,GsendW
            mc=(m-1)*JKTlen
            i=Istr+m-1
            DO l=LBt,UBt
              lc=(l-LBt)*JKlen+mc
              DO k=LBk,UBk
                kc=(k-LBk)*Jlen+lc
                DO j=Jmin,Jmax
                  sizeW=sizeW+1
                  jkW=joff+1+(j-Jmin)+kc
                  sendW(jkW)=B(i,j,k,l)
                END DO
              END DO
            END DO
          END DO
        END IF
        IF (PRESENT(C)) THEN
          joff=jkW
          DO m=1,GsendW
            mc=(m-1)*JKTlen
            i=Istr+m-1
            DO l=LBt,UBt
              lc=(l-LBt)*JKlen+mc
              DO k=LBk,UBk
                kc=(k-LBk)*Jlen+lc
                DO j=Jmin,Jmax
                  sizeW=sizeW+1
                  jkW=joff+1+(j-Jmin)+kc
                  sendW(jkW)=C(i,j,k,l)
                END DO
              END DO
            END DO
          END DO
        END IF
      END IF
!
      IF (Eexchange) THEN
        sizeE=0
        DO m=1,GsendE
          mc=(m-1)*JKTlen
          i=Iend-GsendE+m
          DO l=LBt,UBt
            lc=(l-LBt)*JKlen+mc
            DO k=LBk,UBk
              kc=(k-LBk)*Jlen+lc
              DO j=Jmin,Jmax
                sizeE=sizeE+1
                jkE=1+(j-Jmin)+kc
                sendE(jkE)=A(i,j,k,l)
              END DO
            END DO
          END DO
        END DO
        IF (PRESENT(B)) THEN
          joff=jkE
          DO m=1,GsendE
            mc=(m-1)*JKTlen
            i=Iend-GsendE+m
            DO l=LBt,UBt
              lc=(l-LBt)*JKlen+mc
              DO k=LBk,UBk
                kc=(k-LBk)*Jlen+lc
                DO j=Jmin,Jmax
                  sizeE=sizeE+1
                  jkE=joff+1+(j-Jmin)+kc
                  sendE(jkE)=B(i,j,k,l)
                END DO
              END DO
            END DO
          END DO
        END IF
        IF (PRESENT(C)) THEN
          joff=jkE
          DO m=1,GsendE
            mc=(m-1)*JKTlen
            i=Iend-GsendE+m
            DO l=LBt,UBt
              lc=(l-LBt)*JKlen+mc
              DO k=LBk,UBk
                kc=(k-LBk)*Jlen+lc
                DO j=Jmin,Jmax
                  sizeE=sizeE+1
                  jkE=joff+1+(j-Jmin)+kc
                  sendE(jkE)=C(i,j,k,l)
                END DO
              END DO
            END DO
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Send and receive Western and Eastern segments.
!-----------------------------------------------------------------------
!
# if defined MPI
      IF (Wexchange) THEN
        CALL mpi_irecv (recvW, EWsize, MP_FLOAT, Wtile, Etag,           &
     &                  OCN_COMM_WORLD, Wrequest, Werror)
      END IF
      IF (Eexchange) THEN
        CALL mpi_irecv (recvE, EWsize, MP_FLOAT, Etile, Wtag,           &
     &                  OCN_COMM_WORLD, Erequest, Eerror)
      END IF
      IF (Wexchange) THEN
        CALL mpi_send  (sendW, sizeW, MP_FLOAT, Wtile, Wtag,            &
     &                  OCN_COMM_WORLD, Werror)
      END IF
      IF (Eexchange) THEN
        CALL mpi_send  (sendE, sizeE, MP_FLOAT, Etile, Etag,            &
     &                  OCN_COMM_WORLD, Eerror)
      END IF
# endif
!
!-----------------------------------------------------------------------
!  Unpack Eastern and Western segments.
!-----------------------------------------------------------------------
!
      IF (Wexchange) THEN
# ifdef MPI
        CALL mpi_wait (Wrequest, status(1,1), Werror)
        IF (Werror.ne.MPI_SUCCESS) THEN
          CALL mpi_error_string (Werror, string, Lstr, Ierror)
          Lstr=LEN_TRIM(string)
          WRITE (stdout,20) 'MPI_SEND/MPI_IRECV (Western Edge)',        &
     &                      MyRank, Werror, string(1:Lstr)
 20       FORMAT (/,' MP_EXCHANGE4D - error during ',a,                 &
     &            ' call, Node = ',i3.3,' Error = ',i3,/,15x,a)
          exit_flag=2
          RETURN
        END IF
# endif
        DO m=GrecvW,1,-1
          mc=(GrecvW-m)*JKTlen
          i=Istr-m
          DO l=LBt,UBt
            lc=(l-LBt)*JKlen+mc
            DO k=LBk,UBk
              kc=(k-LBk)*Jlen+lc
              DO j=Jmin,Jmax
                jkW=1+(j-Jmin)+kc
                A(i,j,k,l)=recvW(jkW)
              END DO
            END DO
          END DO
        END DO
        IF (PRESENT(B)) THEN
          joff=jkW
          DO m=GrecvW,1,-1
            mc=(GrecvW-m)*JKTlen
            i=Istr-m
            DO l=LBt,UBt
              lc=(l-LBt)*JKlen+mc
              DO k=LBk,UBk
                kc=(k-LBk)*Jlen+lc
                DO j=Jmin,Jmax
                  jkW=joff+1+(j-Jmin)+kc
                  B(i,j,k,l)=recvW(jkW)
                END DO
              END DO
            END DO
          END DO
        END IF
        IF (PRESENT(C)) THEN
          joff=jkW
          DO m=GrecvW,1,-1
            mc=(GrecvW-m)*JKTlen
            i=Istr-m
            DO l=LBt,UBt
              lc=(l-LBt)*JKlen+mc
              DO k=LBk,UBk
                kc=(k-LBk)*Jlen+lc
                DO j=Jmin,Jmax
                  jkW=joff+1+(j-Jmin)+kc
                  C(i,j,k,l)=recvW(jkW)
                END DO
              END DO
            END DO
          END DO
        END IF
      END IF
!
      IF (Eexchange) THEN
# ifdef MPI
        CALL mpi_wait (Erequest, status(1,3), Eerror)
        IF (Eerror.ne.MPI_SUCCESS) THEN
          CALL mpi_error_string (Eerror, string, Lstr, Ierror)
          Lstr=LEN_TRIM(string)
          WRITE (stdout,20) 'MPI_SEND/MPI_IRECV (Eastern Edge)',        &
     &                      MyRank, Eerror, string(1:Lstr)
          exit_flag=2
          RETURN
        END IF
# endif
        DO m=1,GrecvE
          mc=(m-1)*JKTlen
          i=Iend+m
          DO l=LBt,UBt
            lc=(l-LBt)*JKlen+mc
            DO k=LBk,UBk
              kc=(k-LBk)*Jlen+lc
              DO j=Jmin,Jmax
                jkE=1+(j-Jmin)+kc
                A(i,j,k,l)=recvE(jkE)
              END DO
            END DO
          ENDDO
        END DO
        IF (PRESENT(B)) THEN
          joff=jkE
          DO m=1,GrecvE
            mc=(m-1)*JKTlen
            i=Iend+m
            DO l=LBt,UBt
              lc=(l-LBt)*JKlen+mc
              DO k=LBk,UBk
                kc=(k-LBk)*Jlen+lc
                DO j=Jmin,Jmax
                  jkE=joff+1+(j-Jmin)+kc
                  B(i,j,k,l)=recvE(jkE)
                END DO
              END DO
            ENDDO
          END DO
        END IF
        IF (PRESENT(C)) THEN
          joff=jkE
          DO m=1,GrecvE
            mc=(m-1)*JKTlen
            i=Iend+m
            DO l=LBt,UBt
              lc=(l-LBt)*JKlen+mc
              DO k=LBk,UBk
                kc=(k-LBk)*Jlen+lc
                DO j=Jmin,Jmax
                  jkE=joff+1+(j-Jmin)+kc
                  C(i,j,k,l)=recvE(jkE)
                END DO
              END DO
            ENDDO
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Pack Southern and Northern tile boundary data including ghost-points.
!-----------------------------------------------------------------------
!
      IF (Sexchange) THEN
        sizeS=0
        DO m=1,GsendS
          mc=(m-1)*IKTlen
          j=Jstr+m-1
          DO l=LBt,UBt
            lc=(l-LBt)*IKlen+mc
            DO k=LBk,UBk
              kc=(k-LBk)*Ilen+lc
              DO i=Imin,Imax
                sizeS=sizeS+1
                ikS=1+(i-Imin)+kc
                sendS(ikS)=A(i,j,k,l)
              END DO
            END DO
          END DO
        END DO
        IF (PRESENT(B)) THEN
          ioff=ikS
          DO m=1,GsendS
            mc=(m-1)*IKTlen
            j=Jstr+m-1
            DO l=LBt,UBt
              lc=(l-LBt)*IKlen+mc
              DO k=LBk,UBk
                kc=(k-LBk)*Ilen+lc
                DO i=Imin,Imax
                  sizeS=sizeS+1
                  ikS=ioff+1+(i-Imin)+kc
                  sendS(ikS)=B(i,j,k,l)
                END DO
              END DO
            END DO
          END DO
        END IF
        IF (PRESENT(C)) THEN
          ioff=ikS
          DO m=1,GsendS
            mc=(m-1)*IKTlen
            j=Jstr+m-1
            DO l=LBt,UBt
              lc=(l-LBt)*IKlen+mc
              DO k=LBk,UBk
                kc=(k-LBk)*Ilen+lc
                DO i=Imin,Imax
                  sizeS=sizeS+1
                  ikS=ioff+1+(i-Imin)+kc
                  sendS(ikS)=C(i,j,k,l)
                END DO
              END DO
            END DO
          END DO
        END IF
      END IF
!
      IF (Nexchange) THEN
        sizeN=0
        DO m=1,GsendN
          mc=(m-1)*IKTlen
          j=Jend-GsendN+m
          DO l=LBt,UBt
            lc=(l-LBt)*IKlen+mc
            DO k=LBk,UBk
              kc=(k-LBk)*Ilen+lc
              DO i=Imin,Imax
                sizeN=sizeN+1
                ikN=1+(i-Imin)+kc
                sendN(ikN)=A(i,j,k,l)
              END DO
            END DO
          END DO
        END DO
        IF (PRESENT(B)) THEN
          ioff=ikN
          DO m=1,GsendN
            mc=(m-1)*IKTlen
            j=Jend-GsendN+m
            DO l=LBt,UBt
              lc=(l-LBt)*IKlen+mc
              DO k=LBk,UBk
                kc=(k-LBk)*Ilen+lc
                DO i=Imin,Imax
                  sizeN=sizeN+1
                  ikN=ioff+1+(i-Imin)+kc
                  sendN(ikN)=B(i,j,k,l)
                END DO
              END DO
            END DO
          END DO
        END IF
        IF (PRESENT(C)) THEN
          ioff=ikN
          DO m=1,GsendN
            mc=(m-1)*IKTlen
            j=Jend-GsendN+m
            DO l=LBt,UBt
              lc=(l-LBt)*IKlen+mc
              DO k=LBk,UBk
                kc=(k-LBk)*Ilen+lc
                DO i=Imin,Imax
                  sizeN=sizeN+1
                  ikN=ioff+1+(i-Imin)+kc
                  sendN(ikN)=C(i,j,k,l)
                END DO
              END DO
            END DO
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Send and receive Southern and Northern segments.
!-----------------------------------------------------------------------
!
# if defined MPI
      IF (Sexchange) THEN
        CALL mpi_irecv (recvS, NSsize, MP_FLOAT, Stile, Ntag,           &
     &                  OCN_COMM_WORLD, Srequest, Serror)
      END IF
      IF (Nexchange) THEN
        CALL mpi_irecv (recvN, NSsize, MP_FLOAT, Ntile, Stag,           &
     &                  OCN_COMM_WORLD, Nrequest, Nerror)
      END IF
      IF (Sexchange) THEN
        CALL mpi_send  (sendS, sizeS, MP_FLOAT, Stile, Stag,            &
     &                  OCN_COMM_WORLD, Serror)
      END IF
      IF (Nexchange) THEN
        CALL mpi_send  (sendN, sizeN, MP_FLOAT, Ntile, Ntag,            &
     &                  OCN_COMM_WORLD, Nerror)
      END IF
# endif
!
!-----------------------------------------------------------------------
!  Unpack Northern and Southern segments.
!-----------------------------------------------------------------------
!
      IF (Sexchange) THEN
# ifdef MPI
        CALL mpi_wait (Srequest, status(1,2), Serror)
        IF (Serror.ne.MPI_SUCCESS) THEN
          CALL mpi_error_string (Serror, string, Lstr, Ierror)
          Lstr=LEN_TRIM(string)
          WRITE (stdout,20) 'MPI_SEND/MPI_IRECV (Southern Edge)',       &
     &                      MyRank, Serror, string(1:Lstr)
          exit_flag=2
          RETURN
        END IF
# endif
        DO m=GrecvS,1,-1
          mc=(GrecvS-m)*IKTlen
          j=Jstr-m
          DO l=LBt,UBt
            lc=(l-LBt)*IKlen+mc
            DO k=LBk,UBk
              kc=(k-LBk)*Ilen+lc
              DO i=Imin,Imax
                ikS=1+(i-Imin)+kc
                A(i,j,k,l)=recvS(ikS)
              END DO
            END DO
          END DO
        END DO
        IF (PRESENT(B)) THEN
          ioff=ikS
          DO m=GrecvS,1,-1
            mc=(GrecvS-m)*IKTlen
            j=Jstr-m
            DO l=LBt,UBt
              lc=(l-LBt)*IKlen+mc
              DO k=LBk,UBk
                kc=(k-LBk)*Ilen+lc
                DO i=Imin,Imax
                  ikS=ioff+1+(i-Imin)+kc
                  B(i,j,k,l)=recvS(ikS)
                END DO
              END DO
            END DO
          END DO
        END IF
        IF (PRESENT(C)) THEN
          ioff=ikS
          DO m=GrecvS,1,-1
            mc=(GrecvS-m)*IKTlen
            j=Jstr-m
            DO l=LBt,UBt
              lc=(l-LBt)*IKlen+mc
              DO k=LBk,UBk
                kc=(k-LBk)*Ilen+lc
                DO i=Imin,Imax
                  ikS=ioff+1+(i-Imin)+kc
                  C(i,j,k,l)=recvS(ikS)
                END DO
              END DO
            END DO
          END DO
        END IF
      END IF
!
      IF (Nexchange) THEN
# ifdef MPI
        CALL mpi_wait (Nrequest, status(1,4), Nerror)
        IF (Nerror.ne.MPI_SUCCESS) THEN
          CALL mpi_error_string (Nerror, string, Lstr, Ierror)
          Lstr=LEN_TRIM(string)
          WRITE (stdout,20) 'MPI_SEND/MPI_IRECV (Northern Edge)',       &
     &                      MyRank, Nerror, string(1:Lstr)
          exit_flag=2
          RETURN
        END IF
# endif
        DO m=1,GrecvN
          mc=(m-1)*IKTlen
          j=Jend+m
          DO l=LBt,UBt
            lc=(l-LBt)*IKlen+mc
            DO k=LBk,UBk
              kc=(k-LBk)*Ilen+lc
              DO i=Imin,Imax
                ikN=1+(i-Imin)+kc
                A(i,j,k,l)=recvN(ikN)
              END DO
            END DO
          END DO
        END DO
        IF (PRESENT(B)) THEN
          ioff=ikN
          DO m=1,GrecvN
            mc=(m-1)*IKTlen
            j=Jend+m
            DO l=LBt,UBt
              lc=(l-LBt)*IKlen+mc
              DO k=LBk,UBk
                kc=(k-LBk)*Ilen+lc
                DO i=Imin,Imax
                  ikN=ioff+1+(i-Imin)+kc
                  B(i,j,k,l)=recvN(ikN)
                END DO
              END DO
            END DO
          END DO
        END IF
        IF (PRESENT(C)) THEN
          ioff=ikN
          DO m=1,GrecvN
            mc=(m-1)*IKTlen
            j=Jend+m
            DO l=LBt,UBt
              lc=(l-LBt)*IKlen+mc
              DO k=LBk,UBk
                kc=(k-LBk)*Ilen+lc
                DO i=Imin,Imax
                  ikN=ioff+1+(i-Imin)+kc
                  C(i,j,k,l)=recvN(ikN)
                END DO
              END DO
            END DO
          END DO
        END IF
      END IF
# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn off time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, model, 42)
# endif

      RETURN
      END SUBROUTINE mp_exchange4d

# ifdef ADJOINT
!
!***********************************************************************
      SUBROUTINE ad_mp_exchange2d (ng, tile, model, Nvar,               &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             Nghost, EW_periodic, NS_periodic,    &
     &                             ad_A, ad_B, ad_C, ad_D)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
      implicit none
!
!  Imported variable declarations.
!
      logical, intent(in) :: EW_periodic, NS_periodic

      integer, intent(in) :: ng, tile, model, Nvar
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: Nghost

#  ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: ad_A(LBi:,LBj:)

      real(r8), intent(inout), optional :: ad_B(LBi:,LBj:)
      real(r8), intent(inout), optional :: ad_C(LBi:,LBj:)
      real(r8), intent(inout), optional :: ad_D(LBi:,LBj:)
#  else
      real(r8), intent(inout) :: ad_A(LBi:UBi,LBj:UBj)

      real(r8), intent(inout), optional :: ad_B(LBi:UBi,LBj:UBj)
      real(r8), intent(inout), optional :: ad_C(LBi:UBi,LBj:UBj)
      real(r8), intent(inout), optional :: ad_D(LBi:UBi,LBj:UBj)
#  endif
!
!  Local variable declarations.
!
      logical :: Wexchange, Sexchange, Eexchange, Nexchange

      integer :: i, icS, icN, ioff, Imin, Imax, Ilen
      integer :: j, jcW, jcE, joff, Jmin, Jmax, Jlen
      integer :: m, mc, Ierror, Lstr, pp
      integer :: Wtile, GsendW, GrecvW, Wtag, Werror, Wrequest
      integer :: Stile, GsendS, GrecvS, Stag, Serror, Srequest
      integer :: Etile, GsendE, GrecvE, Etag, Eerror, Erequest
      integer :: Ntile, GsendN, GrecvN, Ntag, Nerror, Nrequest
      integer :: BufferSizeEW, EWsize, sizeW, sizeE
      integer :: BufferSizeNS, NSsize, sizeS, sizeN

#  ifdef MPI
      integer, dimension(MPI_STATUS_SIZE,4) :: status
#  endif
!
      real(r8), dimension(Nvar*HaloSizeJ(ng)) :: sendW, sendE
      real(r8), dimension(Nvar*HaloSizeJ(ng)) :: recvW, recvE

      real(r8), dimension(Nvar*HaloSizeI(ng)) :: sendS, sendN
      real(r8), dimension(Nvar*HaloSizeI(ng)) :: recvS, recvN

      character (len=MPI_MAX_ERROR_STRING) :: string

#  include "set_bounds.h"
#  ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn on time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, model, 40)
#  endif
!
!-----------------------------------------------------------------------
!  Determine rank of tile neighbors and number of ghost-points to
!  exchange.
!-----------------------------------------------------------------------
!
      CALL tile_neighbors (ng, Nghost, EW_periodic, NS_periodic,        &
     &                     GrecvW, GsendW, Wtile, Wexchange,            &
     &                     GrecvE, GsendE, Etile, Eexchange,            &
     &                     GrecvS, GsendS, Stile, Sexchange,            &
     &                     GrecvN, GsendN, Ntile, Nexchange)
!
!  Set communication tags.
!
      Wtag=1
      Stag=2
      Etag=3
      Ntag=4
!
!  Determine range and length of the distributed tile boundary segments.
!
      Imin=LBi
      Imax=UBi
      Jmin=LBj
      Jmax=UBj
      Ilen=Imax-Imin+1
      Jlen=Jmax-Jmin+1
      IF (EW_periodic.or.NS_periodic) THEN
        pp=1
      ELSE
        pp=0
      END IF
      NSsize=Nvar*(Nghost+pp)*Ilen
      EWsize=Nvar*(Nghost+pp)*Jlen
      BufferSizeNS=Nvar*HaloSizeI(ng)
      BufferSizeEW=Nvar*HaloSizeJ(ng)
      IF (SIZE(sendE).lt.EWsize) THEN
        WRITE (stdout,10) 'EWsize = ', EWsize, SIZE(sendE)
 10     FORMAT (/,' AD_MP_EXCHANGE2D - communication buffer too',       &
     &          ' small, ',a, 2i8)
      END IF
      IF (SIZE(sendN).lt.NSsize) THEN
        WRITE (stdout,10) 'NSsize = ', NSsize, SIZE(sendN)
      END IF
!
!-----------------------------------------------------------------------
!  Adjoint of unpacking Northern and Southern segments.
!-----------------------------------------------------------------------
!
      IF (Nexchange) THEN
        DO i=1,BufferSizeNS
          recvN(i)=0.0_r8
          sendN(i)=0.0_r8
        END DO
        sizeN=0
        DO m=1,GrecvN
          mc=(m-1)*Ilen
          j=Jend+m
          DO i=Imin,Imax
            sizeN=sizeN+1
            icN=1+(i-Imin)+mc
!>          A(i,j)=recvN(icN)
!>
            recvN(icN)=ad_A(i,j)
            ad_A(i,j)=0.0_r8
          END DO
        END DO
        IF (PRESENT(ad_B)) THEN
          ioff=icN
          DO m=1,GrecvN
            mc=(m-1)*Ilen
            j=Jend+m
            DO i=Imin,Imax
              sizeN=sizeN+1
              icN=ioff+1+(i-Imin)+mc
!>            B(i,j)=recvN(icN)
!>
              recvN(icN)=ad_B(i,j)
              ad_B(i,j)=0.0_r8
            END DO
          END DO
        END IF
        IF (PRESENT(ad_C)) THEN
          ioff=icN
          DO m=1,GrecvN
            mc=(m-1)*Ilen
            j=Jend+m
            DO i=Imin,Imax
              sizeN=sizeN+1
              icN=ioff+1+(i-Imin)+mc
!>            C(i,j)=recvN(icN)
!>
              recvN(icN)=ad_C(i,j)
              ad_C(i,j)=0.0_r8
            END DO
          END DO
        END IF
        IF (PRESENT(ad_D)) THEN
          ioff=icN
          DO m=1,GrecvN
            mc=(m-1)*Ilen
            j=Jend+m
            DO i=Imin,Imax
              sizeN=sizeN+1
              icN=ioff+1+(i-Imin)+mc
!>            D(i,j)=recvN(icN)
!>
              recvN(icN)=ad_D(i,j)
              ad_D(i,j)=0.0_r8
            END DO
          END DO
        END IF
      END IF
!
      IF (Sexchange) THEN
        DO i=1,BufferSizeNS
          recvS(i)=0.0_r8
          sendS(i)=0.0_r8
        END DO
        sizeS=0
        DO m=GrecvS,1,-1
          mc=(GrecvS-m)*Ilen
          j=Jstr-m
          DO i=Imin,Imax
            sizeS=sizeS+1
            icS=1+(i-Imin)+mc
!>          A(i,j)=recvS(icS)
!>
            recvS(icS)=ad_A(i,j)
            ad_A(i,j)=0.0_r8
          END DO
        END DO
        IF (PRESENT(ad_B)) THEN
          ioff=icS
          DO m=GrecvS,1,-1
            mc=(GrecvS-m)*Ilen
            j=Jstr-m
            DO i=Imin,Imax
              sizeS=sizeS+1
              icS=ioff+1+(i-Imin)+mc
!>            B(i,j)=recvS(icS)
!>
              recvS(icS)=ad_B(i,j)
              ad_B(i,j)=0.0_r8
            END DO
          END DO
        END IF
        IF (PRESENT(ad_C)) THEN
          ioff=icS
          DO m=GrecvS,1,-1
            mc=(GrecvS-m)*Ilen
            j=Jstr-m
            DO i=Imin,Imax
              sizeS=sizeS+1
              icS=ioff+1+(i-Imin)+mc
!>            C(i,j)=recvS(icS)
!>
              recvS(icS)=ad_C(i,j)
              ad_C(i,j)=0.0_r8
            END DO
          END DO
        END IF
        IF (PRESENT(ad_D)) THEN
          ioff=icS
          DO m=GrecvS,1,-1
            mc=(GrecvS-m)*Ilen
            j=Jstr-m
            DO i=Imin,Imax
              sizeS=sizeS+1
              icS=ioff+1+(i-Imin)+mc
!>            D(i,j)=recvS(icS)
!>
              recvS(icS)=ad_D(i,j)
              ad_D(i,j)=0.0_r8
            END DO
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Adjoint of send and receive Southern and Northern segments.
!-----------------------------------------------------------------------
!
#  if defined MPI
      IF (Sexchange) THEN
!>      CALL mpi_irecv (recvS, NSsize, MP_FLOAT, Stile, Ntag,           &
!>   &                  OCN_COMM_WORLD, Srequest, Serror)
!>
        CALL mpi_irecv (sendS, NSsize, MP_FLOAT, Stile, Ntag,           &
     &                  OCN_COMM_WORLD, Srequest, Serror)
      END IF
      IF (Nexchange) THEN
!>      CALL mpi_irecv (recvN, NSsize, MP_FLOAT, Ntile, Stag,           &
!>   &                  OCN_COMM_WORLD, Nrequest, Nerror)
!>
        CALL mpi_irecv (sendN, NSsize, MP_FLOAT, Ntile, Stag,           &
     &                  OCN_COMM_WORLD, Nrequest, Nerror)
      END IF
      IF (Sexchange) THEN
!>      CALL mpi_send  (sendS, sizeS, MP_FLOAT, Stile, Stag,            &
!>   &                  OCN_COMM_WORLD, Serror)
!>
        CALL mpi_send  (recvS, sizeS, MP_FLOAT, Stile, Stag,            &
     &                  OCN_COMM_WORLD, Serror)
      END IF
      IF (Nexchange) THEN
!>      CALL mpi_send  (sendN, sizeN, MP_FLOAT, Ntile, Ntag,            &
!>   &                  OCN_COMM_WORLD, Nerror)
!>
        CALL mpi_send  (recvN, sizeN, MP_FLOAT, Ntile, Ntag,            &
     &                  OCN_COMM_WORLD, Nerror)
      END IF
#  endif
!
!  Adjoint of packing tile boundary data including ghost-points.
!
      IF (Sexchange) THEN
#  ifdef MPI
        CALL mpi_wait (Srequest, status(1,2), Serror)
        IF (Serror.ne.MPI_SUCCESS) THEN
          CALL mpi_error_string (Serror, string, Lstr, Ierror)
          Lstr=LEN_TRIM(string)
          WRITE (stdout,20) 'MPI_SEND/MPI_IRECV (Southern Edge)',       &
     &                      MyRank, Serror, string(1:Lstr)
 20       FORMAT (/,' AD_MP_EXCHANGE2D - error during ',a,' call,',     &
     &            ' Node = ', i3.3,' Error = ',i3,/,18x,a)
          exit_flag=2
          RETURN
        END IF
#  endif
        DO m=1,GsendS
          mc=(m-1)*Ilen
          j=Jstr+m-1
          DO i=Imin,Imax
            icS=1+(i-Imin)+mc
!>          sendS(icS)=A(i,j)
!>
            ad_A(i,j)=ad_A(i,j)+sendS(icS)
            sendS(icS)=0.0_r8
          END DO
        END DO
        IF (PRESENT(ad_B)) THEN
          ioff=icS
          DO m=1,GsendS
            mc=(m-1)*Ilen
            j=Jstr+m-1
            DO i=Imin,Imax
              icS=ioff+1+(i-Imin)+mc
!>            sendS(icS)=B(i,j)
!>
              ad_B(i,j)=ad_B(i,j)+sendS(icS)
              sendS(icS)=0.0_r8
            END DO
          END DO
        END IF
        IF (PRESENT(ad_C)) THEN
          ioff=icS
          DO m=1,GsendS
            mc=(m-1)*Ilen
            j=Jstr+m-1
            DO i=Imin,Imax
              icS=ioff+1+(i-Imin)+mc
!>            sendS(icS)=C(i,j)
!>
              ad_C(i,j)=ad_C(i,j)+sendS(icS)
              sendS(icS)=0.0_r8
            END DO
          END DO
        END IF
        IF (PRESENT(ad_D)) THEN
          ioff=icS
          DO m=1,GsendS
            mc=(m-1)*Ilen
            j=Jstr+m-1
            DO i=Imin,Imax
              icS=ioff+1+(i-Imin)+mc
!>            sendS(icS)=D(i,j)
!>
              ad_D(i,j)=ad_D(i,j)+sendS(icS)
              sendS(icS)=0.0_r8
            END DO
          END DO
        END IF
      END IF
!
      IF (Nexchange) THEN
#  ifdef MPI
        CALL mpi_wait (Nrequest, status(1,4), Nerror)
        IF (Nerror.ne.MPI_SUCCESS) THEN
          CALL mpi_error_string (Nerror, string, Lstr, Ierror)
          Lstr=LEN_TRIM(string)
          WRITE (stdout,20) 'MPI_SEND/MPI_IRECV (Northern Edge)',       &
     &                      MyRank, Nerror, string(1:Lstr)
          exit_flag=2
          RETURN
        END IF
#  endif
        DO m=1,GsendN
          mc=(m-1)*Ilen
          j=Jend-GsendN+m
          DO i=Imin,Imax
            icN=1+(i-Imin)+mc
!>          sendN(icN)=A(i,j)
!>
            ad_A(i,j)=ad_A(i,j)+sendN(icN)
            sendN(icN)=0.0_r8
          END DO
        END DO
        IF (PRESENT(ad_B)) THEN
          ioff=icN
          DO m=1,GsendN
            mc=(m-1)*Ilen
            j=Jend-GsendN+m
            DO i=Imin,Imax
              icN=ioff+1+(i-Imin)+mc
!>            sendN(icN)=B(i,j)
!>
              ad_B(i,j)=ad_B(i,j)+sendN(icN)
              sendN(icN)=0.0_r8
            END DO
          END DO
        END IF
        IF (PRESENT(ad_C)) THEN
          ioff=icN
          DO m=1,GsendN
            mc=(m-1)*Ilen
            j=Jend-GsendN+m
            DO i=Imin,Imax
              icN=ioff+1+(i-Imin)+mc
!>            sendN(icN)=C(i,Jend-GsendN+m)
!>
              ad_C(i,j)=ad_C(i,j)+sendN(icN)
              sendN(icN)=0.0_r8
            END DO
          END DO
        END IF
        IF (PRESENT(ad_D)) THEN
          ioff=icN
          DO m=1,GsendN
            mc=(m-1)*Ilen
            j=Jend-GsendN+m
            DO i=Imin,Imax
              icN=ioff+1+(i-Imin)+mc
!>            sendN(icN)=D(i,j)
!>
              ad_D(i,j)=ad_D(i,j)+sendN(icN)
              sendN(icN)=0.0_r8
            END DO
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Adjoint of unpack Eastern and Western segments.
!-----------------------------------------------------------------------
!
      IF (Eexchange) THEN
        DO i=1,BufferSizeEW
          recvE(i)=0.0_r8
          sendE(i)=0.0_r8
        END DO
        sizeE=0
        DO m=1,GrecvE
          mc=(m-1)*Jlen
          i=Iend+m
          DO j=Jmin,Jmax
            sizeE=sizeE+1
            jcE=1+(j-Jmin)+mc
!>          A(i,j)=recvE(jcE)
!>
            recvE(jcE)=ad_A(i,j)
            ad_A(i,j)=0.0_r8
          ENDDO
        END DO
        IF (PRESENT(ad_B)) THEN
          joff=jcE
          DO m=1,GrecvE
            mc=(m-1)*Jlen
            i=Iend+m
            DO j=Jmin,Jmax
              sizeE=sizeE+1
              jcE=joff+1+(j-Jmin)+mc
!>            B(i,j)=recvE(jcE)
!>
              recvE(jcE)=ad_B(i,j)
              ad_B(i,j)=0.0_r8
            END DO
          END DO
        END IF
        IF (PRESENT(ad_C)) THEN
          joff=jcE
          DO m=1,GrecvE
            mc=(m-1)*Jlen
            i=Iend+m
            DO j=Jmin,Jmax
              sizeE=sizeE+1
              jcE=joff+1+(j-Jmin)+mc
!>            C(i,j)=recvE(jcE)
!>
              recvE(jcE)=ad_C(i,j)
              ad_C(i,j)=0.0_r8
            END DO
          END DO
        END IF
        IF (PRESENT(ad_D)) THEN
          joff=jcE
          DO m=1,GrecvE
            mc=(m-1)*Jlen
            i=Iend+m
            DO j=Jmin,Jmax
              sizeE=sizeE+1
              jcE=joff+1+(j-Jmin)+mc
!>            D(i,j)=recvE(jcE)
!>
              recvE(jcE)=ad_D(i,j)
              ad_D(i,j)=0.0_r8
            END DO
          END DO
        END IF
      END IF
!
      IF (Wexchange) THEN
        DO i=1,BufferSizeEW
          recvW(i)=0.0_r8
          sendW(i)=0.0_r8
        END DO
        sizeW=0
        DO m=GrecvW,1,-1
          mc=(GrecvW-m)*Jlen
          i=Istr-m
          DO j=Jmin,Jmax
            sizeW=sizeW+1
            jcW=1+(j-Jmin)+mc
!>          A(i,j)=recvW(jcW)
!>
            recvW(jcW)=ad_A(i,j)
            ad_A(i,j)=0.0_r8
          END DO
        END DO
        IF (PRESENT(ad_B)) THEN
          joff=jcW
          DO m=GrecvW,1,-1
            mc=(GrecvW-m)*Jlen
            i=Istr-m
            DO j=Jmin,Jmax
              sizeW=sizeW+1
              jcW=joff+1+(j-Jmin)+mc
!>            B(i,j)=recvW(jcW)
!>
              recvW(jcW)=ad_B(i,j)
              ad_B(i,j)=0.0_r8
            END DO
          END DO
        END IF
        IF (PRESENT(ad_C)) THEN
          joff=jcW
          DO m=GrecvW,1,-1
            mc=(GrecvW-m)*Jlen
            i=Istr-m
            DO j=Jmin,Jmax
              sizeW=sizeW+1
              jcW=joff+1+(j-Jmin)+mc
!>            C(i,j)=recvW(jcW)
!>
              recvW(jcW)=ad_C(i,j)
              ad_C(i,j)=0.0_r8
            END DO
          END DO
        END IF
        IF (PRESENT(ad_D)) THEN
          joff=jcW
          DO m=GrecvW,1,-1
            mc=(GrecvW-m)*Jlen
            i=Istr-m
            DO j=Jmin,Jmax
              sizeW=sizeW+1
              jcW=joff+1+(j-Jmin)+mc
!>            D(i,j)=recvW(jcW)
!>
              recvW(jcW)=ad_D(i,j)
              ad_D(i,j)=0.0_r8
            END DO
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Send and receive Western and Eastern segments.
!-----------------------------------------------------------------------
!
#  if defined MPI
      IF (Wexchange) THEN
!>      CALL mpi_irecv (recvW, EWsize, MP_FLOAT, Wtile, Etag,           &
!>   &                  OCN_COMM_WORLD, Wrequest, Werror)
!>
        CALL mpi_irecv (sendW, EWsize, MP_FLOAT, Wtile, Etag,           &
     &                  OCN_COMM_WORLD, Wrequest, Werror)
      END IF
      IF (Eexchange) THEN
!>      CALL mpi_irecv (recvE, EWsize, MP_FLOAT, Etile, Wtag,           &
!>   &                  OCN_COMM_WORLD, Erequest, Eerror)
!>
        CALL mpi_irecv (sendE, EWsize, MP_FLOAT, Etile, Wtag,           &
     &                  OCN_COMM_WORLD, Erequest, Eerror)
      END IF
      IF (Wexchange) THEN
!>      CALL mpi_send  (sendW, sizeW, MP_FLOAT, Wtile, Wtag,            &
!>   &                  OCN_COMM_WORLD, Werror)
!>
        CALL mpi_send  (recvW, sizeW, MP_FLOAT, Wtile, Wtag,            &
     &                  OCN_COMM_WORLD, Werror)
      END IF
      IF (Eexchange) THEN
!>      CALL mpi_send  (sendE, sizeE, MP_FLOAT, Etile, Etag,            &
!>   &                  OCN_COMM_WORLD, Eerror)
!>
        CALL mpi_send  (recvE, sizeE, MP_FLOAT, Etile, Etag,            &
     &                  OCN_COMM_WORLD, Eerror)
      END IF
#  endif
!
!  Adjoint of packing tile boundary data including ghost-points.
!
      IF (Wexchange) THEN
#  ifdef MPI
        CALL mpi_wait (Wrequest, status(1,1), Werror)
        IF (Werror.ne.MPI_SUCCESS) THEN
          CALL mpi_error_string (Werror, string, Lstr, Ierror)
          Lstr=LEN_TRIM(string)
          WRITE (stdout,20) 'MPI_SEND/MPI_IRECV (Western Edge)',        &
     &                      MyRank, Werror, string(1:Lstr)
          exit_flag=2
          RETURN
        END IF
#  endif
        DO m=1,GsendW
          mc=(m-1)*Jlen
          i=Istr+m-1
          DO j=Jmin,Jmax
            jcW=1+(j-Jmin)+mc
!>          sendW(jcW)=A(i,j)
!>
            ad_A(i,j)=ad_A(i,j)+sendW(jcW)
            sendW(jcW)=0.0_r8
          END DO
        END DO
        IF (PRESENT(ad_B)) THEN
          joff=jcW
          DO m=1,GsendW
            mc=(m-1)*Jlen
            i=Istr+m-1
            DO j=Jmin,Jmax
              jcW=joff+1+(j-Jmin)+mc
!>            sendW(jcW)=B(i,j)
!>
              ad_B(i,j)=ad_B(i,j)+sendW(jcW)
              sendW(jcW)=0.0_r8
            END DO
          END DO
        END IF
        IF (PRESENT(ad_C)) THEN
          joff=jcW
          DO m=1,GsendW
            mc=(m-1)*Jlen
            i=Istr+m-1
            DO j=Jmin,Jmax
              jcW=joff+1+(j-Jmin)+mc
!>            sendW(jcW)=C(i,j)
!>
              ad_C(i,j)=ad_C(i,j)+sendW(jcW)
              sendW(jcW)=0.0_r8
            END DO
          END DO
        END IF
        IF (PRESENT(ad_D)) THEN
          joff=jcW
          DO m=1,GsendW
            mc=(m-1)*Jlen
            i=Istr+m-1
            DO j=Jmin,Jmax
              jcW=joff+1+(j-Jmin)+mc
!>            sendW(jcW)=D(i,j)
!>
              ad_D(i,j)=ad_D(i,j)+sendW(jcW)
              sendW(jcW)=0.0_r8
            END DO
          END DO
        END IF
      END IF
!
      IF (Eexchange) THEN
#  ifdef MPI
        CALL mpi_wait (Erequest, status(1,3), Eerror)
        IF (Eerror.ne.MPI_SUCCESS) THEN
          CALL mpi_error_string (Eerror, string, Lstr, Ierror)
          Lstr=LEN_TRIM(string)
          WRITE (stdout,20) 'MPI_SEND/MPI_IRECV (Eastern Edge)',        &
     &                      MyRank, Eerror, string(1:Lstr)
          exit_flag=2
          RETURN
        END IF
#  endif
        DO m=1,GsendE
          mc=(m-1)*Jlen
          i=Iend-GsendE+m
          DO j=Jmin,Jmax
            jcE=1+(j-Jmin)+mc
!>          sendE(jcE)=A(i,j)
!>
            ad_A(i,j)=ad_A(i,j)+sendE(jcE)
            sendE(jcE)=0.0_r8
          END DO
        END DO
        IF (PRESENT(ad_B)) THEN
          joff=jcE
          DO m=1,GsendE
            mc=(m-1)*Jlen
            i=Iend-GsendE+m
            DO j=Jmin,Jmax
              jcE=joff+1+(j-Jmin)+mc
!>            sendE(jcE)=B(i,j)
!>
              ad_B(i,j)=ad_B(i,j)+sendE(jcE)
              sendE(jcE)=0.0_r8
            END DO
          END DO
        END IF
        IF (PRESENT(ad_C)) THEN
          joff=jcE
          DO m=1,GsendE
            mc=(m-1)*Jlen
            i=Iend-GsendE+m
            DO j=Jmin,Jmax
              jcE=joff+1+(j-Jmin)+mc
!>            sendE(jcE)=C(i,j)
!>
              ad_C(i,j)=ad_C(i,j)+sendE(jcE)
              sendE(jcE)=0.0_r8
            END DO
          END DO
        END IF
        IF (PRESENT(ad_D)) THEN
          joff=jcE
          DO m=1,GsendE
            mc=(m-1)*Jlen
            i=Iend-GsendE+m
            DO j=Jmin,Jmax
              jcE=joff+1+(j-Jmin)+mc
!>            sendE(jcE)=D(i,j)
!>
              ad_D(i,j)=ad_D(i,j)+sendE(jcE)
              sendE(jcE)=0.0_r8
            END DO
          END DO
        END IF
      END IF
#  ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn off time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, model, 40)
#  endif

      RETURN
      END SUBROUTINE ad_mp_exchange2d

!
!***********************************************************************
      SUBROUTINE ad_mp_exchange2d_bry (ng, tile, model, Nvar, boundary, &
     &                                 LBij, UBij,                      &
     &                                 Nghost, EW_periodic, NS_periodic,&
     &                                 ad_A, ad_B, ad_C, ad_D)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
      implicit none
!
!  Imported variable declarations.
!
      logical, intent(in) :: EW_periodic, NS_periodic

      integer, intent(in) :: ng, tile, model, Nvar, boundary
      integer, intent(in) :: LBij, UBij
      integer, intent(in) :: Nghost
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: ad_A(LBij:)

      real(r8), intent(inout), optional :: ad_B(LBij:)
      real(r8), intent(inout), optional :: ad_C(LBij:)
      real(r8), intent(inout), optional :: ad_D(LBij:)
# else
      real(r8), intent(inout) :: ad_A(LBij:UBij)

      real(r8), intent(inout), optional :: ad_B(LBij:UBij)
      real(r8), intent(inout), optional :: ad_C(LBij:UBij)
      real(r8), intent(inout), optional :: ad_D(LBij:UBij)
# endif
!
!  Local variable declarations.
!
      logical :: Wexchange, Sexchange, Eexchange, Nexchange

      integer :: i, icS, icN
      integer :: j, jcW, jcE
      integer :: m, Ierror, Lstr, pp
      integer :: Wtile, GsendW, GrecvW, Wtag, Werror, Wrequest
      integer :: Stile, GsendS, GrecvS, Stag, Serror, Srequest
      integer :: Etile, GsendE, GrecvE, Etag, Eerror, Erequest
      integer :: Ntile, GsendN, GrecvN, Ntag, Nerror, Nrequest
      integer :: BufferSizeEW, EWsize, sizeW, sizeE
      integer :: BufferSizeNS, NSsize, sizeS, sizeN

#  ifdef MPI
      integer, dimension(MPI_STATUS_SIZE,4) :: status
#  endif
!
      real(r8), dimension(Nvar*HaloBry(ng)) :: sendW, sendE
      real(r8), dimension(Nvar*HaloBry(ng)) :: recvW, recvE

      real(r8), dimension(Nvar*HaloBry(ng)) :: sendS, sendN
      real(r8), dimension(Nvar*HaloBry(ng)) :: recvS, recvN

      character (len=MPI_MAX_ERROR_STRING) :: string

#  include "set_bounds.h"
#  ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn on time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, model, 43)
#  endif
!
!-----------------------------------------------------------------------
!  Determine rank of tile neighbors and number of ghost-points to
!  exchange.
!-----------------------------------------------------------------------
!
      CALL tile_neighbors (ng, Nghost, EW_periodic, NS_periodic,        &
     &                     GrecvW, GsendW, Wtile, Wexchange,            &
     &                     GrecvE, GsendE, Etile, Eexchange,            &
     &                     GrecvS, GsendS, Stile, Sexchange,            &
     &                     GrecvN, GsendN, Ntile, Nexchange)
!
!  Adjust exchange swiches according to boundary edge to process.
!
      Wexchange=Wexchange.and.((boundary.eq.isouth).or.                 &
     &                         (boundary.eq.inorth))
      Eexchange=Eexchange.and.((boundary.eq.isouth).or.                 &
     &                         (boundary.eq.inorth))
      Sexchange=Sexchange.and.((boundary.eq.iwest).or.                  &
     &                         (boundary.eq.ieast))
      Nexchange=Nexchange.and.((boundary.eq.iwest).or.                  &
     &                         (boundary.eq.ieast))
!
!  Set communication tags.
!
      Wtag=1
      Stag=2
      Etag=3
      Ntag=4
!
!  Determine range and length of the distributed tile boundary segments.
!
      IF (EW_periodic.or.NS_periodic) THEN
        pp=1
      ELSE
        pp=0
      END IF
      NSsize=Nvar*(Nghost+pp)
      EWsize=Nvar*(Nghost+pp)
      BufferSizeNS=Nvar*(Nghost+pp)
      BufferSizeEW=Nvar*(Nghost+pp)
      IF (SIZE(sendE).lt.EWsize) THEN
        WRITE (stdout,10) 'EWsize = ', EWsize, SIZE(sendE)
 10     FORMAT (/,' AD_MP_EXCHANGE2D_BRY - communication buffer too',   &
     &          ' small, ',a, 2i8)
      END IF
      IF (SIZE(sendN).lt.NSsize) THEN
        WRITE (stdout,10) 'NSsize = ', NSsize, SIZE(sendN)
      END IF
!
!-----------------------------------------------------------------------
!  Adjoint of unpacking Northern and Southern segments.
!-----------------------------------------------------------------------
!
      IF (Nexchange) THEN
        DO i=1,BufferSizeNS
          recvN(i)=0.0_r8
          sendN(i)=0.0_r8
        END DO
        icN=0
        sizeN=0
        DO m=1,GrecvN
          j=Jend+m
          sizeN=sizeN+1
          icN=icN+1
!>        A(j)=recvN(icN)
!>
          recvN(icN)=ad_A(j)
          ad_A(j)=0.0_r8
        END DO
        IF (PRESENT(ad_B)) THEN
          DO m=1,GrecvN
            j=Jend+m
            sizeN=sizeN+1
            icN=icN+1
!>          B(j)=recvN(icN)
!>
            recvN(icN)=ad_B(j)
            ad_B(j)=0.0_r8
          END DO
        END IF
        IF (PRESENT(ad_C)) THEN
          DO m=1,GrecvN
            j=Jend+m
            sizeN=sizeN+1
            icN=icN+1
!>          C(j)=recvN(icN)
!>
            recvN(icN)=ad_C(j)
            ad_C(j)=0.0_r8
          END DO
        END IF
        IF (PRESENT(ad_D)) THEN
          DO m=1,GrecvN
            j=Jend+m
            sizeN=sizeN+1
            icN=icN+1
!>          D(j)=recvN(icN)
!>
            recvN(icN)=ad_D(j)
            ad_D(j)=0.0_r8
          END DO
        END IF
      END IF
!
      IF (Sexchange) THEN
        DO i=1,BufferSizeNS
          recvS(i)=0.0_r8
          sendS(i)=0.0_r8
        END DO
        icS=0
        sizeS=0
        DO m=GrecvS,1,-1
          j=Jstr-m
          sizeS=sizeS+1
          icS=icS+1
!>        A(j)=recvS(icS)
!>
          recvS(icS)=ad_A(j)
          ad_A(j)=0.0_r8
        END DO
        IF (PRESENT(ad_B)) THEN
          DO m=GrecvS,1,-1
            j=Jstr-m
            sizeS=sizeS+1
            icS=icS+1
!>          B(j)=recvS(icS)
!>
            recvS(icS)=ad_B(j)
            ad_B(j)=0.0_r8
          END DO
        END IF
        IF (PRESENT(ad_C)) THEN
          DO m=GrecvS,1,-1
            j=Jstr-m
            sizeS=sizeS+1
            icS=icS+1
!>          C(j)=recvS(icS)
!>
            recvS(icS)=ad_C(j)
            ad_C(j)=0.0_r8
          END DO
        END IF
        IF (PRESENT(ad_D)) THEN
          DO m=GrecvS,1,-1
            j=Jstr-m
            sizeS=sizeS+1
            icS=icS+1
!>          D(j)=recvS(icS)
!>
            recvS(icS)=ad_D(j)
            ad_D(j)=0.0_r8
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Adjoint of send and receive Southern and Northern segments.
!-----------------------------------------------------------------------
!
#  if defined MPI
      IF (Sexchange) THEN
!>      CALL mpi_irecv (recvS, NSsize, MP_FLOAT, Stile, Ntag,           &
!>   &                  OCN_COMM_WORLD, Srequest, Serror)
!>
        CALL mpi_irecv (sendS, NSsize, MP_FLOAT, Stile, Ntag,           &
     &                  OCN_COMM_WORLD, Srequest, Serror)
      END IF
      IF (Nexchange) THEN
!>      CALL mpi_irecv (recvN, NSsize, MP_FLOAT, Ntile, Stag,           &
!>   &                  OCN_COMM_WORLD, Nrequest, Nerror)
!>
        CALL mpi_irecv (sendN, NSsize, MP_FLOAT, Ntile, Stag,           &
     &                  OCN_COMM_WORLD, Nrequest, Nerror)
      END IF
      IF (Sexchange) THEN
!>      CALL mpi_send  (sendS, sizeS, MP_FLOAT, Stile, Stag,            &
!>   &                  OCN_COMM_WORLD, Serror)
!>
        CALL mpi_send  (recvS, sizeS, MP_FLOAT, Stile, Stag,            &
     &                  OCN_COMM_WORLD, Serror)
      END IF
      IF (Nexchange) THEN
!>      CALL mpi_send  (sendN, sizeN, MP_FLOAT, Ntile, Ntag,            &
!>   &                  OCN_COMM_WORLD, Nerror)
!>
        CALL mpi_send  (recvN, sizeN, MP_FLOAT, Ntile, Ntag,            &
     &                  OCN_COMM_WORLD, Nerror)
      END IF
#  endif
!
!  Adjoint of packing tile boundary data including ghost-points.
!
      IF (Sexchange) THEN
#  ifdef MPI
        CALL mpi_wait (Srequest, status(1,2), Serror)
        IF (Serror.ne.MPI_SUCCESS) THEN
          CALL mpi_error_string (Serror, string, Lstr, Ierror)
          Lstr=LEN_TRIM(string)
          WRITE (stdout,20) 'MPI_SEND/MPI_IRECV (Southern Edge)',       &
     &                      MyRank, Serror, string(1:Lstr)
 20       FORMAT (/,' AD_MP_EXCHANGE2D_BRY - error during ',a,' call,', &
     &            ' Node = ', i3.3,' Error = ',i3,/,18x,a)
          exit_flag=2
          RETURN
        END IF
#  endif
        icS=0
        DO m=1,GsendS
          j=Jstr+m-1
          icS=icS+1
!>        sendS(icS)=A(j)
!>
          ad_A(j)=ad_A(j)+sendS(icS)
          sendS(icS)=0.0_r8
        END DO
        IF (PRESENT(ad_B)) THEN
          DO m=1,GsendS
            j=Jstr+m-1
            icS=icS+1
!>          sendS(icS)=B(j)
!>
            ad_B(j)=ad_B(j)+sendS(icS)
            sendS(icS)=0.0_r8
          END DO
        END IF
        IF (PRESENT(ad_C)) THEN
          DO m=1,GsendS
            j=Jstr+m-1
            icS=ics+1
!>          sendS(icS)=C(j)
!>
            ad_C(j)=ad_C(j)+sendS(icS)
            sendS(icS)=0.0_r8
          END DO
        END IF
        IF (PRESENT(ad_D)) THEN
          DO m=1,GsendS
            j=Jstr+m-1
            icS=icS+1
!>          sendS(icS)=D(j)
!>
            ad_D(j)=ad_D(j)+sendS(icS)
            sendS(icS)=0.0_r8
          END DO
        END IF
      END IF
!
      IF (Nexchange) THEN
#  ifdef MPI
        CALL mpi_wait (Nrequest, status(1,4), Nerror)
        IF (Nerror.ne.MPI_SUCCESS) THEN
          CALL mpi_error_string (Nerror, string, Lstr, Ierror)
          Lstr=LEN_TRIM(string)
          WRITE (stdout,20) 'MPI_SEND/MPI_IRECV (Northern Edge)',       &
     &                      MyRank, Nerror, string(1:Lstr)
          exit_flag=2
          RETURN
        END IF
#  endif
        icN=0
        DO m=1,GsendN
          j=Jend-GsendN+m
          icN=icN+1
!>        sendN(icN)=A(j)
!>
          ad_A(j)=ad_A(j)+sendN(icN)
          sendN(icN)=0.0_r8
        END DO
        IF (PRESENT(ad_B)) THEN
          DO m=1,GsendN
            j=Jend-GsendN+m
            icN=icN+1
!>          sendN(icN)=B(j)
!>
            ad_B(j)=ad_B(j)+sendN(icN)
            sendN(icN)=0.0_r8
          END DO
        END IF
        IF (PRESENT(ad_C)) THEN
          DO m=1,GsendN
            j=Jend-GsendN+m
            icN=icN+1
!>          sendN(icN)=C(j)
!>
            ad_C(j)=ad_C(j)+sendN(icN)
            sendN(icN)=0.0_r8
          END DO
        END IF
        IF (PRESENT(ad_D)) THEN
          DO m=1,GsendN
            j=Jend-GsendN+m
            icN=icN+1
!>          sendN(icN)=D(j)
!>
            ad_D(j)=ad_D(j)+sendN(icN)
            sendN(icN)=0.0_r8
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Adjoint of unpack Eastern and Western segments.
!-----------------------------------------------------------------------
!
      IF (Eexchange) THEN
        DO i=1,BufferSizeEW
          recvE(i)=0.0_r8
          sendE(i)=0.0_r8
        END DO
        jcE=0
        sizeE=0
        DO m=1,GrecvE
          i=Iend+m
          sizeE=sizeE+1
          jcE=jcE+1
!>        A(i)=recvE(jcE)
!>
          recvE(jcE)=ad_A(i)
          ad_A(i)=0.0_r8
        END DO
        IF (PRESENT(ad_B)) THEN
          DO m=1,GrecvE
            i=Iend+m
            sizeE=sizeE+1
            jcE=jcE+1
!>          B(i)=recvE(jcE)
!>
            recvE(jcE)=ad_B(i)
            ad_B(i)=0.0_r8
          END DO
        END IF
        IF (PRESENT(ad_C)) THEN
          DO m=1,GrecvE
            i=Iend+m
            sizeE=sizeE+1
            jcE=jcE+1
!>          C(i)=recvE(jcE)
!>
            recvE(jcE)=ad_C(i)
            ad_C(i)=0.0_r8
          END DO
        END IF
        IF (PRESENT(ad_D)) THEN
          DO m=1,GrecvE
            i=Iend+m
            sizeE=sizeE+1
            jcE=jcE+1
!>          D(i)=recvE(jcE)
!>
            recvE(jcE)=ad_D(i)
            ad_D(i)=0.0_r8
          END DO
        END IF
      END IF
!
      IF (Wexchange) THEN
        DO i=1,BufferSizeEW
          recvW(i)=0.0_r8
          sendW(i)=0.0_r8
        END DO
        jcW=0
        sizeW=0
        DO m=GrecvW,1,-1
          i=Istr-m
          sizeW=sizeW+1
          jcW=jcW+1
!>        A(i)=recvW(jcW)
!>
          recvW(jcW)=ad_A(i)
          ad_A(i)=0.0_r8
        END DO
        IF (PRESENT(ad_B)) THEN
          DO m=GrecvW,1,-1
            i=Istr-m
            sizeW=sizeW+1
            jcW=jcW+1
!>          B(i)=recvW(jcW)
!>
            recvW(jcW)=ad_B(i)
            ad_B(i)=0.0_r8
          END DO
        END IF
        IF (PRESENT(ad_C)) THEN
          DO m=GrecvW,1,-1
            i=Istr-m
            sizeW=sizeW+1
            jcW=jcW+1
!>          C(i)=recvW(jcW)
!>
            recvW(jcW)=ad_C(i)
            ad_C(i)=0.0_r8
          END DO
        END IF
        IF (PRESENT(ad_D)) THEN
          DO m=GrecvW,1,-1
            i=Istr-m
            sizeW=sizeW+1
            jcW=jcW+1
!>          D(i)=recvW(jcW)
!>
            recvW(jcW)=ad_D(i)
            ad_D(i)=0.0_r8
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Send and receive Western and Eastern segments.
!-----------------------------------------------------------------------
!
#  if defined MPI
      IF (Wexchange) THEN
!>      CALL mpi_irecv (recvW, EWsize, MP_FLOAT, Wtile, Etag,           &
!>   &                  OCN_COMM_WORLD, Wrequest, Werror)
!>
        CALL mpi_irecv (sendW, EWsize, MP_FLOAT, Wtile, Etag,           &
     &                  OCN_COMM_WORLD, Wrequest, Werror)
      END IF
      IF (Eexchange) THEN
!>      CALL mpi_irecv (recvE, EWsize, MP_FLOAT, Etile, Wtag,           &
!>   &                  OCN_COMM_WORLD, Erequest, Eerror)
!>
        CALL mpi_irecv (sendE, EWsize, MP_FLOAT, Etile, Wtag,           &
     &                  OCN_COMM_WORLD, Erequest, Eerror)
      END IF
      IF (Wexchange) THEN
!>      CALL mpi_send  (sendW, sizeW, MP_FLOAT, Wtile, Wtag,            &
!>   &                  OCN_COMM_WORLD, Werror)
!>
        CALL mpi_send  (recvW, sizeW, MP_FLOAT, Wtile, Wtag,            &
     &                  OCN_COMM_WORLD, Werror)
      END IF
      IF (Eexchange) THEN
!>      CALL mpi_send  (sendE, sizeE, MP_FLOAT, Etile, Etag,            &
!>   &                  OCN_COMM_WORLD, Eerror)
!>
        CALL mpi_send  (recvE, sizeE, MP_FLOAT, Etile, Etag,            &
     &                  OCN_COMM_WORLD, Eerror)
      END IF
#  endif
!
!  Adjoint of packing tile boundary data including ghost-points.
!
      IF (Wexchange) THEN
#  ifdef MPI
        CALL mpi_wait (Wrequest, status(1,1), Werror)
        IF (Werror.ne.MPI_SUCCESS) THEN
          CALL mpi_error_string (Werror, string, Lstr, Ierror)
          Lstr=LEN_TRIM(string)
          WRITE (stdout,20) 'MPI_SEND/MPI_IRECV (Western Edge)',        &
     &                      MyRank, Werror, string(1:Lstr)
          exit_flag=2
          RETURN
        END IF
#  endif
        jcW=0
        DO m=1,GsendW
          i=Istr+m-1
          jcW=jcW+1
!>        sendW(jcW)=A(i)
!>
          ad_A(i)=ad_A(i)+sendW(jcW)
          sendW(jcW)=0.0_r8
        END DO
        IF (PRESENT(ad_B)) THEN
          DO m=1,GsendW
            i=Istr+m-1
            jcW=jcW+1
!>          sendW(jcW)=B(i)
!>
            ad_B(i)=ad_B(i)+sendW(jcW)
            sendW(jcW)=0.0_r8
          END DO
        END IF
        IF (PRESENT(ad_C)) THEN
          DO m=1,GsendW
            i=Istr+m-1
            jcW=jcW+1
!>          sendW(jcW)=C(i)
!>
            ad_C(i)=ad_C(i)+sendW(jcW)
            sendW(jcW)=0.0_r8
          END DO
        END IF
        IF (PRESENT(ad_D)) THEN
          DO m=1,GsendW
            i=Istr+m-1
            jcW=jcW+1
!>          sendW(jcW)=D(i)
!>
            ad_D(i)=ad_D(i)+sendW(jcW)
            sendW(jcW)=0.0_r8
          END DO
        END IF
      END IF
!
      IF (Eexchange) THEN
#  ifdef MPI
        CALL mpi_wait (Erequest, status(1,3), Eerror)
        IF (Eerror.ne.MPI_SUCCESS) THEN
          CALL mpi_error_string (Eerror, string, Lstr, Ierror)
          Lstr=LEN_TRIM(string)
          WRITE (stdout,20) 'MPI_SEND/MPI_IRECV (Eastern Edge)',        &
     &                      MyRank, Eerror, string(1:Lstr)
          exit_flag=2
          RETURN
        END IF
#  endif
        jcE=0
        DO m=1,GsendE
          i=Iend-GsendE+m
          jcE=jcE+1
!>        sendE(jcE)=A(i)
!>
          ad_A(i)=ad_A(i)+sendE(jcE)
          sendE(jcE)=0.0_r8
        END DO
        IF (PRESENT(ad_B)) THEN
          DO m=1,GsendE
            i=Iend-GsendE+m
            jcE=jcE+1
!>          sendE(jcE)=B(i)
!>
            ad_B(i)=ad_B(i)+sendE(jcE)
            sendE(jcE)=0.0_r8
          END DO
        END IF
        IF (PRESENT(ad_C)) THEN
          DO m=1,GsendE
            i=Iend-GsendE+m
            jcE=jcE+1
!>          sendE(jcE)=C(i)
!>
            ad_C(i)=ad_C(i)+sendE(jcE)
            sendE(jcE)=0.0_r8
          END DO
        END IF
        IF (PRESENT(ad_D)) THEN
          DO m=1,GsendE
            i=Iend-GsendE+m
            jcE=jcE+1
!>          sendE(jcE)=D(i)
!>
            ad_D(i)=ad_D(i)+sendE(jcE)
            sendE(jcE)=0.0_r8
          END DO
        END IF
      END IF

#  ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn off time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, model, 43)
#  endif

      RETURN
      END SUBROUTINE ad_mp_exchange2d_bry

!
!***********************************************************************
      SUBROUTINE ad_mp_exchange3d (ng, tile, model, Nvar,               &
     &                             LBi, UBi, LBj, UBj, LBk, UBk,        &
     &                             Nghost, EW_periodic, NS_periodic,    &
     &                             ad_A, ad_B, ad_C, ad_D)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
      implicit none
!
!  Imported variable declarations.
!
      logical, intent(in) :: EW_periodic, NS_periodic

      integer, intent(in) :: ng, tile, model, Nvar
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBk, UBk
      integer, intent(in) :: Nghost
!
#  ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: ad_A(LBi:,LBj:,LBk:)

      real(r8), intent(inout), optional :: ad_B(LBi:,LBj:,LBk:)
      real(r8), intent(inout), optional :: ad_C(LBi:,LBj:,LBk:)
      real(r8), intent(inout), optional :: ad_D(LBi:,LBj:,LBk:)
#  else
      real(r8), intent(inout) :: ad_A(LBi:UBi,LBj:UBj,LBk:UBk)

      real(r8), intent(inout), optional :: ad_B(LBi:UBi,LBj:UBj,LBk:UBk)
      real(r8), intent(inout), optional :: ad_C(LBi:UBi,LBj:UBj,LBk:UBk)
      real(r8), intent(inout), optional :: ad_D(LBi:UBi,LBj:UBj,LBk:UBk)
#  endif
!
!  Local variable declarations.
!
      logical :: Wexchange, Sexchange, Eexchange, Nexchange

      integer :: i, ikS, ikN, ioff, Imin, Imax, Ilen, IKlen
      integer :: j, jkW, jkE, joff, Jmin, Jmax, Jlen, JKlen
      integer :: k, kc, m, mc, Ierror, Klen, Lstr, pp
      integer :: Wtile, GsendW, GrecvW, Wtag, Werror, Wrequest
      integer :: Stile, GsendS, GrecvS, Stag, Serror, Srequest
      integer :: Etile, GsendE, GrecvE, Etag, Eerror, Erequest
      integer :: Ntile, GsendN, GrecvN, Ntag, Nerror, Nrequest
      integer :: BufferSizeEW, EWsize, sizeW, sizeE
      integer :: BufferSizeNS, NSsize, sizeS, sizeN

#  ifdef MPI
      integer, dimension(MPI_STATUS_SIZE,4) :: status
#  endif
!
      real(r8), dimension(Nvar*HaloSizeJ(ng)*                           &
     &                           (UBk-LBk+1)) :: sendW, sendE
      real(r8), dimension(Nvar*HaloSizeI(ng)*                           &
     &                           (UBk-LBk+1)) :: sendS, sendN

      real(r8), dimension(Nvar*HaloSizeJ(ng)*                           &
     &                           (UBk-LBk+1)) :: recvW, recvE
      real(r8), dimension(Nvar*HaloSizeI(ng)*                           &
     &                           (UBk-LBk+1)) :: recvS, recvN

      character (len=MPI_MAX_ERROR_STRING) :: string

#  include "set_bounds.h"
#  ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn on time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, model, 41)
#  endif
!
!-----------------------------------------------------------------------
!  Determine rank of tile neighbors and number of ghost-points to
!  exchange.
!-----------------------------------------------------------------------
!
      CALL tile_neighbors (ng, Nghost, EW_periodic, NS_periodic,        &
     &                     GrecvW, GsendW, Wtile, Wexchange,            &
     &                     GrecvE, GsendE, Etile, Eexchange,            &
     &                     GrecvS, GsendS, Stile, Sexchange,            &
     &                     GrecvN, GsendN, Ntile, Nexchange)
!
!  Set communication tags.
!
      Wtag=1
      Stag=2
      Etag=3
      Ntag=4
!
!  Determine range and length of the distributed tile boundary segments.
!
      Imin=LBi
      Imax=UBi
      Jmin=LBj
      Jmax=UBj
      Ilen=Imax-Imin+1
      Jlen=Jmax-Jmin+1
      Klen=UBk-LBk+1
      IKlen=Ilen*Klen
      JKlen=Jlen*Klen
      IF (EW_periodic.or.NS_periodic) THEN
        pp=1
      ELSE
        pp=0
      END IF
      NSsize=Nvar*(Nghost+pp)*IKlen
      EWsize=Nvar*(Nghost+pp)*JKlen
      BufferSizeNS=Nvar*HaloSizeI(ng)*Klen
      BufferSizeEW=Nvar*HaloSizeJ(ng)*Klen
      IF (SIZE(sendE).lt.EWsize) THEN
        WRITE (stdout,10) 'EWsize = ', EWsize, SIZE(sendE)
 10     FORMAT (/,' AD_MP_EXCHANGE3D - communication buffer too',       &
     &          ' small, ',a, 2i8)
      END IF
      IF (SIZE(sendN).lt.NSsize) THEN
        WRITE (stdout,10) 'NSsize = ', NSsize, SIZE(sendN)
      END IF
!
!-----------------------------------------------------------------------
!  Adjoint of unpacking Northern and Southern segments.
!-----------------------------------------------------------------------
!
      IF (Nexchange) THEN
        DO i=1,BufferSizeNS
          recvN(i)=0.0_r8
          sendN(i)=0.0_r8
        END DO
        sizeN=0
        DO m=1,GrecvN
          mc=(m-1)*IKlen
          j=Jend+m
          DO k=LBk,UBk
            kc=(k-LBk)*Ilen+mc
            DO i=Imin,Imax
              sizeN=sizeN+1
              ikN=1+(i-Imin)+kc
!>            A(i,j,k)=recvN(ikN)
!>
              recvN(ikN)=ad_A(i,j,k)
              ad_A(i,j,k)=0.0_r8
            END DO
          END DO
        END DO
        IF (PRESENT(ad_B)) THEN
          ioff=ikN
          DO m=1,GrecvN
            mc=(m-1)*IKlen
            j=Jend+m
            DO k=LBk,UBk
              kc=(k-LBk)*Ilen+mc
              DO i=Imin,Imax
                sizeN=sizeN+1
                ikN=ioff+1+(i-Imin)+kc
!>              B(i,j,k)=recvN(ikN)
!>
                recvN(ikN)=ad_B(i,j,k)
                ad_B(i,j,k)=0.0_r8
              END DO
            END DO
          END DO
        END IF
        IF (PRESENT(ad_C)) THEN
          ioff=ikN
          DO m=1,GrecvN
            mc=(m-1)*IKlen
            j=Jend+m
            DO k=LBk,UBk
              kc=(k-LBk)*Ilen+mc
              DO i=Imin,Imax
                sizeN=sizeN+1
                ikN=ioff+1+(i-Imin)+kc
!>              C(i,j,k)=recvN(ikN)
!>
                recvN(ikN)=ad_C(i,j,k)
                ad_C(i,j,k)=0.0_r8
              END DO
            END DO
          END DO
        END IF
        IF (PRESENT(ad_D)) THEN
          ioff=ikN
          DO m=1,GrecvN
            mc=(m-1)*IKlen
            j=Jend+m
            DO k=LBk,UBk
              kc=(k-LBk)*Ilen+mc
              DO i=Imin,Imax
                sizeN=sizeN+1
                ikN=ioff+1+(i-Imin)+kc
!>              D(i,j,k)=recvN(ikN)
!>
                recvN(ikN)=ad_D(i,j,k)
                ad_D(i,j,k)=0.0_r8
              END DO
            END DO
          END DO
        END IF
      END IF
!
      IF (Sexchange) THEN
        DO i=1,BufferSizeNS
          recvS(i)=0.0_r8
          sendS(i)=0.0_r8
        END DO
        sizeS=0
        DO m=GrecvS,1,-1
          mc=(GrecvS-m)*IKlen
          j=Jstr-m
          DO k=LBk,UBk
            kc=(k-LBk)*Ilen+mc
            DO i=Imin,Imax
              sizeS=sizeS+1
              ikS=1+(i-Imin)+kc
!>            A(i,j,k)=recvS(ikS)
!>
              recvS(ikS)=ad_A(i,j,k)
              ad_A(i,j,k)=0.0_r8
            END DO
          END DO
        END DO
        IF (PRESENT(ad_B)) THEN
          ioff=ikS
          DO m=GrecvS,1,-1
            mc=(GrecvS-m)*IKlen
            j=Jstr-m
            DO k=LBk,UBk
              kc=(k-LBk)*Ilen+mc
              DO i=Imin,Imax
                sizeS=sizeS+1
                ikS=ioff+1+(i-Imin)+kc
!>              B(i,j,k)=recvS(ikS)
!>
                recvS(ikS)=ad_B(i,j,k)
                ad_B(i,j,k)=0.0_r8
              END DO
            END DO
          END DO
        END IF
        IF (PRESENT(ad_C)) THEN
          ioff=ikS
          DO m=GrecvS,1,-1
            mc=(GrecvS-m)*IKlen
            j=Jstr-m
            DO k=LBk,UBk
              kc=(k-LBk)*Ilen+mc
              DO i=Imin,Imax
                sizeS=sizeS+1
                ikS=ioff+1+(i-Imin)+kc
!>              C(i,j,k)=recvS(ikS)
!>
                recvS(ikS)=ad_C(i,j,k)
                ad_C(i,j,k)=0.0_r8
              END DO
            END DO
          END DO
        END IF
        IF (PRESENT(ad_D)) THEN
          ioff=ikS
          DO m=GrecvS,1,-1
            mc=(GrecvS-m)*IKlen
            j=Jstr-m
            DO k=LBk,UBk
              kc=(k-LBk)*Ilen+mc
              DO i=Imin,Imax
                sizeS=sizeS+1
                ikS=ioff+1+(i-Imin)+kc
!>              D(i,j,k)=recvS(ikS)
!>
                recvS(ikS)=ad_D(i,j,k)
                ad_D(i,j,k)=0.0_r8
              END DO
            END DO
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Adjoint of send and receive Southern and Northern segments.
!-----------------------------------------------------------------------
!
#  if defined MPI
      IF (Sexchange) THEN
!>      CALL mpi_irecv (recvS, NSsize, MP_FLOAT, Stile, Ntag,           &
!>   &                  OCN_COMM_WORLD, Srequest, Serror)
!>
        CALL mpi_irecv (sendS, NSsize, MP_FLOAT, Stile, Ntag,           &
     &                  OCN_COMM_WORLD, Srequest, Serror)
      END IF
      IF (Nexchange) THEN
!>      CALL mpi_irecv (recvN, NSsize, MP_FLOAT, Ntile, Stag,           &
!>   &                  OCN_COMM_WORLD, Nrequest, Nerror)
!>
        CALL mpi_irecv (sendN, NSsize, MP_FLOAT, Ntile, Stag,           &
     &                  OCN_COMM_WORLD, Nrequest, Nerror)
      END IF
      IF (Sexchange) THEN
!>      CALL mpi_send  (sendS, sizeS, MP_FLOAT, Stile, Stag,            &
!>   &                  OCN_COMM_WORLD, Serror)
!>
        CALL mpi_send  (recvS, sizeS, MP_FLOAT, Stile, Stag,            &
     &                  OCN_COMM_WORLD, Serror)
      END IF
      IF (Nexchange) THEN
!>      CALL mpi_send  (sendN, sizeN, MP_FLOAT, Ntile, Ntag,            &
!>   &                  OCN_COMM_WORLD, Nerror)
!>
        CALL mpi_send  (recvN, sizeN, MP_FLOAT, Ntile, Ntag,            &
     &                  OCN_COMM_WORLD, Nerror)
      END IF
#  endif
!
!  Adjoint of packing tile boundary data including ghost-points.
!
      IF (Sexchange) THEN
#  ifdef MPI
        CALL mpi_wait (Srequest, status(1,2), Serror)
        IF (Serror.ne.MPI_SUCCESS) THEN
          CALL mpi_error_string (Serror, string, Lstr, Ierror)
          Lstr=LEN_TRIM(string)
          WRITE (stdout,20) 'MPI_SEND/MPI_IRECV (Southern Edge)',       &
     &                      MyRank, Serror, string(1:Lstr)
 20       FORMAT (/,' AD_MP_EXCHANGE3D - error during ',a,' call,',     &
     &            ' Node = ', i3.3,' Error = ',i3,/,18x,a)
          exit_flag=2
          RETURN
        END IF
#  endif
        DO m=1,GsendS
          mc=(m-1)*IKlen
          j=Jstr+m-1
          DO k=LBk,UBk
            kc=(k-LBk)*Ilen+mc
            DO i=Imin,Imax
              ikS=1+(i-Imin)+kc
!>            sendS(ikS)=A(i,j,k)
!>
              ad_A(i,j,k)=ad_A(i,j,k)+sendS(ikS)
              sendS(ikS)=0.0_r8
            END DO
          END DO
        END DO
        IF (PRESENT(ad_B)) THEN
          ioff=ikS
          DO m=1,GsendS
            mc=(m-1)*IKlen
            j=Jstr+m-1
            DO k=LBk,UBk
              kc=(k-LBk)*Ilen+mc
              DO i=Imin,Imax
                ikS=ioff+1+(i-Imin)+kc
!>              sendS(ikS)=B(i,j,k)
!>
                ad_B(i,j,k)=ad_B(i,j,k)+sendS(ikS)
                sendS(ikS)=0.0_r8
              END DO
            END DO
          END DO
        END IF
        IF (PRESENT(ad_C)) THEN
          ioff=ikS
          DO m=1,GsendS
            mc=(m-1)*IKlen
            j=Jstr+m-1
            DO k=LBk,UBk
              kc=(k-LBk)*Ilen+mc
              DO i=Imin,Imax
                ikS=ioff+1+(i-Imin)+kc
!>              sendS(ikS)=C(i,j,k)
!>
                ad_C(i,j,k)=ad_C(i,j,k)+sendS(ikS)
                sendS(ikS)=0.0_r8
              END DO
            END DO
          END DO
        END IF
        IF (PRESENT(ad_D)) THEN
          ioff=ikS
          DO m=1,GsendS
            mc=(m-1)*IKlen
            j=Jstr+m-1
            DO k=LBk,UBk
              kc=(k-LBk)*Ilen+mc
              DO i=Imin,Imax
                ikS=ioff+1+(i-Imin)+kc
!>              sendS(ikS)=D(i,j,k)
!>
                ad_D(i,j,k)=ad_D(i,j,k)+sendS(ikS)
                sendS(ikS)=0.0_r8
              END DO
            END DO
          END DO
        END IF
      END IF
!
      IF (Nexchange) THEN
#  ifdef MPI
        CALL mpi_wait (Nrequest, status(1,4), Nerror)
        IF (Nerror.ne.MPI_SUCCESS) THEN
          CALL mpi_error_string (Nerror, string, Lstr, Ierror)
          Lstr=LEN_TRIM(string)
          WRITE (stdout,20) 'MPI_SEND/MPI_IRECV (Northern Edge)',       &
     &                      MyRank, Nerror, string(1:Lstr)
          exit_flag=2
          RETURN
        END IF
#  endif
        DO m=1,GsendN
          mc=(m-1)*IKlen
          j=Jend-GsendN+m
          DO k=LBk,UBk
            kc=(k-LBk)*Ilen+mc
            DO i=Imin,Imax
              ikN=1+(i-Imin)+kc
!>            sendN(ikN)=A(i,j,k)
!>
              ad_A(i,j,k)=ad_A(i,j,k)+sendN(ikN)
              sendN(ikN)=0.0_r8
            END DO
          END DO
        END DO
        IF (PRESENT(ad_B)) THEN
          ioff=ikN
          DO m=1,GsendN
            mc=(m-1)*IKlen
            j=Jend-GsendN+m
            DO k=LBk,UBk
              kc=(k-LBk)*Ilen+mc
              DO i=Imin,Imax
                ikN=ioff+1+(i-Imin)+kc
!>              sendN(ikN)=B(i,j,k)
!>
                ad_B(i,j,k)=ad_B(i,j,k)+sendN(ikN)
                sendN(ikN)=0.0_r8
              END DO
            END DO
          END DO
        END IF
        IF (PRESENT(ad_C)) THEN
          ioff=ikN
          DO m=1,GsendN
            mc=(m-1)*IKlen
            j=Jend-GsendN+m
            DO k=LBk,UBk
              kc=(k-LBk)*Ilen+mc
              DO i=Imin,Imax
                ikN=ioff+1+(i-Imin)+kc
!>              sendN(ikN)=C(i,j,k)
!>
                ad_C(i,j,k)=ad_C(i,j,k)+sendN(ikN)
                sendN(ikN)=0.0_r8
              END DO
            END DO
          END DO
        END IF
        IF (PRESENT(ad_D)) THEN
          ioff=ikN
          DO m=1,GsendN
            mc=(m-1)*IKlen
            j=Jend-GsendN+m
            DO k=LBk,UBk
              kc=(k-LBk)*Ilen+mc
              DO i=Imin,Imax
                ikN=ioff+1+(i-Imin)+kc
!>              sendN(ikN)=D(i,j,k)
!>
                ad_D(i,j,k)=ad_D(i,j,k)+sendN(ikN)
                sendN(ikN)=0.0_r8
              END DO
            END DO
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Adjoint of unpack Eastern and Western segments.
!-----------------------------------------------------------------------
!
      IF (Eexchange) THEN
        DO i=1,BufferSizeEW
          recvE(i)=0.0_r8
          sendE(i)=0.0_r8
        END DO
        sizeE=0
        DO m=1,GrecvE
          mc=(m-1)*JKlen
          i=Iend+m
          DO k=LBk,UBk
            kc=(k-LBk)*Jlen+mc
            DO j=Jmin,Jmax
              sizeE=sizeE+1
              jkE=1+(j-Jmin)+kc
!>            A(i,j,k)=recvE(jkE)
!>
              recvE(jkE)=ad_A(i,j,k)
              ad_A(i,j,k)=0.0_r8
            END DO
          ENDDO
        END DO
        IF (PRESENT(ad_B)) THEN
          joff=jkE
          DO m=1,GrecvE
            mc=(m-1)*JKlen
            i=Iend+m
            DO k=LBk,UBk
              kc=(k-LBk)*Jlen+mc
              DO j=Jmin,Jmax
                sizeE=sizeE+1
                jkE=joff+1+(j-Jmin)+kc
!>              B(i,j,k)=recvE(jkE)
!>
                recvE(jkE)=ad_B(i,j,k)
                ad_B(i,j,k)=0.0_r8
              END DO
            END DO
          END DO
        END IF
        IF (PRESENT(ad_C)) THEN
          joff=jkE
          DO m=1,GrecvE
            mc=(m-1)*JKlen
            i=Iend+m
            DO k=LBk,UBk
              kc=(k-LBk)*Jlen+mc
              DO j=Jmin,Jmax
                sizeE=sizeE+1
                jkE=joff+1+(j-Jmin)+kc
!>              C(i,j,k)=recvE(jkE)
!>
                recvE(jkE)=ad_C(i,j,k)
                ad_C(i,j,k)=0.0_r8
              END DO
            END DO
          END DO
        END IF
        IF (PRESENT(ad_D)) THEN
          joff=jkE
          DO m=1,GrecvE
            mc=(m-1)*JKlen
            i=Iend+m
            DO k=LBk,UBk
              kc=(k-LBk)*Jlen+mc
              DO j=Jmin,Jmax
                sizeE=sizeE+1
                jkE=joff+1+(j-Jmin)+kc
!>              D(i,j,k)=recvE(jkE)
!>
                recvE(jkE)=ad_D(i,j,k)
                ad_D(i,j,k)=0.0_r8
              END DO
            END DO
          END DO
        END IF
      END IF
!
      IF (Wexchange) THEN
        DO i=1,BufferSizeEW
          recvW(i)=0.0_r8
          sendW(i)=0.0_r8
        END DO
        sizeW=0
        DO m=GrecvW,1,-1
          mc=(GrecvW-m)*JKlen
          i=Istr-m
          DO k=LBk,UBk
            kc=(k-LBk)*Jlen+mc
            DO j=Jmin,Jmax
              sizeW=sizeW+1
              jkW=1+(j-Jmin)+kc
!>            A(i,j,k)=recvW(jkW)
!>
              recvW(jkW)=ad_A(i,j,k)
              ad_A(i,j,k)=0.0_r8
            END DO
          END DO
        END DO
        IF (PRESENT(ad_B)) THEN
          joff=jkW
          DO m=GrecvW,1,-1
            mc=(GrecvW-m)*JKlen
            i=Istr-m
            DO k=LBk,UBk
              kc=(k-LBk)*Jlen+mc
              DO j=Jmin,Jmax
                sizeW=sizeW+1
                jkW=joff+1+(j-Jmin)+kc
!>              B(i,j,k)=recvW(jkW)
!>
                recvW(jkW)=ad_B(i,j,k)
                ad_B(i,j,k)=0.0_r8
              END DO
            END DO
          END DO
        END IF
        IF (PRESENT(ad_C)) THEN
          joff=jkW
          DO m=GrecvW,1,-1
            mc=(GrecvW-m)*JKlen
            i=Istr-m
            DO k=LBk,UBk
              kc=(k-LBk)*Jlen+mc
              DO j=Jmin,Jmax
                sizeW=sizeW+1
                jkW=joff+1+(j-Jmin)+kc
!>              C(i,j,k)=recvW(jkW)
!>
                recvW(jkW)=ad_C(i,j,k)
                ad_C(i,j,k)=0.0_r8
              END DO
            END DO
          END DO
        END IF
        IF (PRESENT(ad_D)) THEN
          joff=jkW
          DO m=GrecvW,1,-1
            mc=(GrecvW-m)*JKlen
            i=Istr-m
            DO k=LBk,UBk
              kc=(k-LBk)*Jlen+mc
              DO j=Jmin,Jmax
                sizeW=sizeW+1
                jkW=joff+1+(j-Jmin)+kc
!>              D(i,j,k)=recvW(jkW)
!>
                recvW(jkW)=ad_D(i,j,k)
                ad_D(i,j,k)=0.0_r8
              END DO
            END DO
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Send and receive Western and Eastern segments.
!-----------------------------------------------------------------------
!
#  if defined MPI
      IF (Wexchange) THEN
!>      CALL mpi_irecv (recvW, EWsize, MP_FLOAT, Wtile, Etag,           &
!>   &                  OCN_COMM_WORLD, Wrequest, Werror)
!>
        CALL mpi_irecv (sendW, EWsize, MP_FLOAT, Wtile, Etag,           &
     &                  OCN_COMM_WORLD, Wrequest, Werror)
      END IF
      IF (Eexchange) THEN
!>      CALL mpi_irecv (recvE, EWsize, MP_FLOAT, Etile, Wtag,           &
!>   &                  OCN_COMM_WORLD, Erequest, Eerror)
!>
        CALL mpi_irecv (sendE, EWsize, MP_FLOAT, Etile, Wtag,           &
     &                  OCN_COMM_WORLD, Erequest, Eerror)
      END IF
      IF (Wexchange) THEN
!>      CALL mpi_send  (sendW, sizeW, MP_FLOAT, Wtile, Wtag,            &
!>   &                  OCN_COMM_WORLD, Werror)
!>
        CALL mpi_send  (recvW, sizeW, MP_FLOAT, Wtile, Wtag,            &
     &                  OCN_COMM_WORLD, Werror)
      END IF
      IF (Eexchange) THEN
!>      CALL mpi_send  (sendE, sizeE, MP_FLOAT, Etile, Etag,            &
!>   &                  OCN_COMM_WORLD, Eerror)
!>
        CALL mpi_send  (recvE, sizeE, MP_FLOAT, Etile, Etag,            &
     &                  OCN_COMM_WORLD, Eerror)
      END IF
#  endif
!
!  Adjoint of packing tile boundary data including ghost-points.
!
      IF (Wexchange) THEN
#  ifdef MPI
        CALL mpi_wait (Wrequest, status(1,1), Werror)
        IF (Werror.ne.MPI_SUCCESS) THEN
          CALL mpi_error_string (Werror, string, Lstr, Ierror)
          Lstr=LEN_TRIM(string)
          WRITE (stdout,20) 'MPI_SEND/MPI_IRECV (Western Edge)',        &
     &                      MyRank, Werror, string(1:Lstr)
          exit_flag=2
          RETURN
        END IF
#  endif
        DO m=1,GsendW
          mc=(m-1)*JKlen
          i=Istr+m-1
          DO k=LBk,UBk
            kc=(k-LBk)*Jlen+mc
            DO j=Jmin,Jmax
              jkW=1+(j-Jmin)+kc
!>            sendW(jkW)=A(i,j,k)
!>
              ad_A(i,j,k)=ad_A(i,j,k)+sendW(jkW)
              sendW(jkW)=0.0_r8
            END DO
          END DO
        END DO
        IF (PRESENT(ad_B)) THEN
          joff=jkW
          DO m=1,GsendW
            mc=(m-1)*JKlen
            i=Istr+m-1
            DO k=LBk,UBk
              kc=(k-LBk)*Jlen+mc
              DO j=Jmin,Jmax
                jkW=joff+1+(j-Jmin)+kc
!>              sendW(jkW)=B(i,j,k)
!>
                ad_B(i,j,k)=ad_B(i,j,k)+sendW(jkW)
                sendW(jkW)=0.0_r8
              END DO
            END DO
          END DO
        END IF
        IF (PRESENT(ad_C)) THEN
          joff=jkW
          DO m=1,GsendW
            mc=(m-1)*JKlen
            i=Istr+m-1
            DO k=LBk,UBk
              kc=(k-LBk)*Jlen+mc
              DO j=Jmin,Jmax
                jkW=joff+1+(j-Jmin)+kc
!>              sendW(jkW)=C(i,j,k)
!>
                ad_C(i,j,k)=ad_C(i,j,k)+sendW(jkW)
                sendW(jkW)=0.0_r8
              END DO
            END DO
          END DO
        END IF
        IF (PRESENT(ad_D)) THEN
          joff=jkW
          DO m=1,GsendW
            mc=(m-1)*JKlen
            i=Istr+m-1
            DO k=LBk,UBk
              kc=(k-LBk)*Jlen+mc
              DO j=Jmin,Jmax
                jkW=joff+1+(j-Jmin)+kc
!>              sendW(jkW)=D(i,j,k)
!>
                ad_D(i,j,k)=ad_D(i,j,k)+sendW(jkW)
                sendW(jkW)=0.0_r8
              END DO
            END DO
          END DO
        END IF
      END IF
!
      IF (Eexchange) THEN
#  ifdef MPI
        CALL mpi_wait (Erequest, status(1,3), Eerror)
        IF (Eerror.ne.MPI_SUCCESS) THEN
          CALL mpi_error_string (Eerror, string, Lstr, Ierror)
          Lstr=LEN_TRIM(string)
          WRITE (stdout,20) 'MPI_SEND/MPI_IRECV (Eastern Edge)',        &
     &                      MyRank, Eerror, string(1:Lstr)
          exit_flag=2
          RETURN
        END IF
#  endif
        DO m=1,GsendE
          mc=(m-1)*JKlen
          i=Iend-GsendE+m
          DO k=LBk,UBk
            kc=(k-LBk)*Jlen+mc
            DO j=Jmin,Jmax
              jkE=1+(j-Jmin)+kc
!>            sendE(jkE)=A(i,j,k)
!>
              ad_A(i,j,k)=ad_A(i,j,k)+sendE(jkE)
              sendE(jkE)=0.0_r8
            END DO
          END DO
        END DO
        IF (PRESENT(ad_B)) THEN
          joff=jkE
          DO m=1,GsendE
            mc=(m-1)*JKlen
            i=Iend-GsendE+m
            DO k=LBk,UBk
              kc=(k-LBk)*Jlen+mc
              DO j=Jmin,Jmax
                jkE=joff+1+(j-Jmin)+kc
!>              sendE(jkE)=B(i,j,k)
!>
                ad_B(i,j,k)=ad_B(i,j,k)+sendE(jkE)
                sendE(jkE)=0.0_r8
              END DO
            END DO
          END DO
        END IF
        IF (PRESENT(ad_C)) THEN
          joff=jkE
          DO m=1,GsendE
            mc=(m-1)*JKlen
            i=Iend-GsendE+m
            DO k=LBk,UBk
              kc=(k-LBk)*Jlen+mc
              DO j=Jmin,Jmax
                jkE=joff+1+(j-Jmin)+kc
!>              sendE(jkE)=C(i,j,k)
!>
                ad_C(i,j,k)=ad_C(i,j,k)+sendE(jkE)
                sendE(jkE)=0.0_r8
              END DO
            END DO
          END DO
        END IF
        IF (PRESENT(ad_D)) THEN
          joff=jkE
          DO m=1,GsendE
            mc=(m-1)*JKlen
            i=Iend-GsendE+m
            DO k=LBk,UBk
              kc=(k-LBk)*Jlen+mc
              DO j=Jmin,Jmax
                jkE=joff+1+(j-Jmin)+kc
!>              sendE(jkE)=D(i,j,k)
!>
                ad_D(i,j,k)=ad_D(i,j,k)+sendE(jkE)
                sendE(jkE)=0.0_r8
              END DO
            END DO
          END DO
        END IF
      END IF
#  ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn off time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, model, 41)
#  endif

      RETURN
      END SUBROUTINE ad_mp_exchange3d
!
!***********************************************************************
      SUBROUTINE ad_mp_exchange3d_bry (ng, tile, model, Nvar, boundary, &
     &                                 LBij, UBij, LBk, UBk,            &
     &                                 Nghost, EW_periodic, NS_periodic,&
     &                                 ad_A, ad_B, ad_C, ad_D)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
      implicit none
!
!  Imported variable declarations.
!
      logical, intent(in) :: EW_periodic, NS_periodic

      integer, intent(in) :: ng, tile, model, Nvar, boundary
      integer, intent(in) :: LBij, UBij, LBk, UBk
      integer, intent(in) :: Nghost
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: ad_A(LBij:,LBk:)

      real(r8), intent(inout), optional :: ad_B(LBij:,LBk:)
      real(r8), intent(inout), optional :: ad_C(LBij:,LBk:)
      real(r8), intent(inout), optional :: ad_D(LBij:,LBk:)
# else
      real(r8), intent(inout) :: ad_A(LBij:UBij,LBk:UBk)

      real(r8), intent(inout), optional :: ad_B(LBij:UBij,LBk:UBk)
      real(r8), intent(inout), optional :: ad_C(LBij:UBij,LBk:UBk)
      real(r8), intent(inout), optional :: ad_D(LBij:UBij,LBk:UBk)
# endif
!
!  Local variable declarations.
!
      logical :: Wexchange, Sexchange, Eexchange, Nexchange

      integer :: i, ikS, ikN, ioff
      integer :: j, jkW, jkE, joff
      integer :: k, m, mc, Ierror, Klen, Lstr, pp
      integer :: Wtile, GsendW, GrecvW, Wtag, Werror, Wrequest
      integer :: Stile, GsendS, GrecvS, Stag, Serror, Srequest
      integer :: Etile, GsendE, GrecvE, Etag, Eerror, Erequest
      integer :: Ntile, GsendN, GrecvN, Ntag, Nerror, Nrequest
      integer :: BufferSizeEW, EWsize, sizeW, sizeE
      integer :: BufferSizeNS, NSsize, sizeS, sizeN

#  ifdef MPI
      integer, dimension(MPI_STATUS_SIZE,4) :: status
#  endif
!
      real(r8), dimension(Nvar*HaloBry(ng)*(UBk-LBk+1)) :: sendW, sendE
      real(r8), dimension(Nvar*HaloBry(ng)*(UBk-LBk+1)) :: recvW, recvE
      real(r8), dimension(Nvar*HaloBry(ng)*(UBk-LBk+1)) :: sendS, sendN
      real(r8), dimension(Nvar*HaloBry(ng)*(UBk-LBk+1)) :: recvS, recvN

      character (len=MPI_MAX_ERROR_STRING) :: string

#  include "set_bounds.h"

#  ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn on time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, model, 43)
#  endif
!
!-----------------------------------------------------------------------
!  Determine rank of tile neighbors and number of ghost-points to
!  exchange.
!-----------------------------------------------------------------------
!
      CALL tile_neighbors (ng, Nghost, EW_periodic, NS_periodic,        &
     &                     GrecvW, GsendW, Wtile, Wexchange,            &
     &                     GrecvE, GsendE, Etile, Eexchange,            &
     &                     GrecvS, GsendS, Stile, Sexchange,            &
     &                     GrecvN, GsendN, Ntile, Nexchange)
!
!  Adjust exchange swiches according to boundary edge to process.
!
      Wexchange=Wexchange.and.((boundary.eq.isouth).or.                 &
     &                         (boundary.eq.inorth))
      Eexchange=Eexchange.and.((boundary.eq.isouth).or.                 &
     &                         (boundary.eq.inorth))
      Sexchange=Sexchange.and.((boundary.eq.iwest).or.                  &
     &                         (boundary.eq.ieast))
      Nexchange=Nexchange.and.((boundary.eq.iwest).or.                  &
     &                         (boundary.eq.ieast))
!
!  Set communication tags.
!
      Wtag=1
      Stag=2
      Etag=3
      Ntag=4
!
!  Determine range and length of the distributed tile boundary segments.
!
      Klen=UBk-LBk+1
      IF (EW_periodic.or.NS_periodic) THEN
        pp=1
      ELSE
        pp=0
      END IF
      NSsize=Nvar*(Nghost+pp)*Klen
      EWsize=Nvar*(Nghost+pp)*Klen
      BufferSizeNS=Nvar*(Nghost+pp)*Klen
      BufferSizeEW=Nvar*(Nghost+pp)*Klen
      IF (SIZE(sendE).lt.EWsize) THEN
        WRITE (stdout,10) 'EWsize = ', EWsize, SIZE(sendE)
 10     FORMAT (/,' AD_MP_EXCHANGE3D_BRY - communication buffer too',   &
     &          ' small, ',a, 2i8)
      END IF
      IF (SIZE(sendN).lt.NSsize) THEN
        WRITE (stdout,10) 'NSsize = ', NSsize, SIZE(sendN)
      END IF
!
!-----------------------------------------------------------------------
!  Adjoint of unpacking Northern and Southern segments.
!-----------------------------------------------------------------------
!
      IF (Nexchange) THEN
        DO i=1,BufferSizeNS
          recvN(i)=0.0_r8
          sendN(i)=0.0_r8
        END DO
        sizeN=0
        DO m=1,GrecvN
          mc=(m-1)*Klen
          j=Jend+m
          DO k=LBk,UBk
            sizeN=sizeN+1
            ikN=1+(k-LBk)+mc
!>          A(j,k)=recvN(ikN)
!>
            recvN(ikN)=ad_A(j,k)
            ad_A(j,k)=0.0_r8
          END DO
        END DO
        IF (PRESENT(ad_B)) THEN
          ioff=ikN
          DO m=1,GrecvN
            mc=(m-1)*Klen
            j=Jend+m
            DO k=LBk,UBk
              sizeN=sizeN+1
              ikN=ioff+1+(k-LBk)+mc
!>            B(j,k)=recvN(ikN)
!>
              recvN(ikN)=ad_B(j,k)
              ad_B(j,k)=0.0_r8
            END DO
          END DO
        END IF
        IF (PRESENT(ad_C)) THEN
          ioff=ikN
          DO m=1,GrecvN
            mc=(m-1)*Klen
            j=Jend+m
            DO k=LBk,UBk
              sizeN=sizeN+1
              ikN=ioff+1+(k-LBk)+mc
!>            C(j,k)=recvN(ikN)
!>
              recvN(ikN)=ad_C(j,k)
              ad_C(j,k)=0.0_r8
            END DO
          END DO
        END IF
        IF (PRESENT(ad_D)) THEN
          ioff=ikN
          DO m=1,GrecvN
            mc=(m-1)*Klen
            j=Jend+m
            DO k=LBk,UBk
              sizeN=sizeN+1
              ikN=ioff+1+(k-LBk)+mc
!>            D(j,k)=recvN(ikN)
!>
              recvN(ikN)=ad_D(j,k)
              ad_D(j,k)=0.0_r8
            END DO
          END DO
        END IF
      END IF
!
      IF (Sexchange) THEN
        DO i=1,BufferSizeNS
          recvS(i)=0.0_r8
          sendS(i)=0.0_r8
        END DO
        sizeS=0
        DO m=GrecvS,1,-1
          mc=(GrecvS-m)*Klen
          j=Jstr-m
          DO k=LBk,UBk
            sizeS=sizeS+1
            ikS=1+(k-LBk)+mc
!>          A(j,k)=recvS(ikS)
!>
            recvS(ikS)=ad_A(j,k)
            ad_A(j,k)=0.0_r8
          END DO
        END DO
        IF (PRESENT(ad_B)) THEN
          ioff=ikS
          DO m=GrecvS,1,-1
            mc=(GrecvS-m)*Klen
            j=Jstr-m
            DO k=LBk,UBk
              sizeS=sizeS+1
              ikS=ioff+1+(k-LBk)+mc
!>            B(j,k)=recvS(ikS)
!>
              recvS(ikS)=ad_B(j,k)
              ad_B(j,k)=0.0_r8
            END DO
          END DO
        END IF
        IF (PRESENT(ad_C)) THEN
          ioff=ikS
          DO m=GrecvS,1,-1
            mc=(GrecvS-m)*Klen
            j=Jstr-m
            DO k=LBk,UBk
              sizeS=sizeS+1
              ikS=ioff+1+(k-LBk)+mc
!>            C(j,k)=recvS(ikS)
!>
              recvS(ikS)=ad_C(j,k)
              ad_C(j,k)=0.0_r8
            END DO
          END DO
        END IF
        IF (PRESENT(ad_D)) THEN
          ioff=ikS
          DO m=GrecvS,1,-1
            mc=(GrecvS-m)*Klen
            j=Jstr-m
            DO k=LBk,UBk
              sizeS=sizeS+1
              ikS=ioff+1+(k-LBk)+mc
!>            D(j,k)=recvS(ikS)
!>
              recvS(ikS)=ad_D(j,k)
              ad_D(j,k)=0.0_r8
            END DO
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Adjoint of send and receive Southern and Northern segments.
!-----------------------------------------------------------------------
!
#  if defined MPI
      IF (Sexchange) THEN
!>      CALL mpi_irecv (recvS, NSsize, MP_FLOAT, Stile, Ntag,           &
!>   &                  OCN_COMM_WORLD, Srequest, Serror)
!>
        CALL mpi_irecv (sendS, NSsize, MP_FLOAT, Stile, Ntag,           &
     &                  OCN_COMM_WORLD, Srequest, Serror)
      END IF
      IF (Nexchange) THEN
!>      CALL mpi_irecv (recvN, NSsize, MP_FLOAT, Ntile, Stag,           &
!>   &                  OCN_COMM_WORLD, Nrequest, Nerror)
!>
        CALL mpi_irecv (sendN, NSsize, MP_FLOAT, Ntile, Stag,           &
     &                  OCN_COMM_WORLD, Nrequest, Nerror)
      END IF
      IF (Sexchange) THEN
!>      CALL mpi_send  (sendS, sizeS, MP_FLOAT, Stile, Stag,            &
!>   &                  OCN_COMM_WORLD, Serror)
!>
        CALL mpi_send  (recvS, sizeS, MP_FLOAT, Stile, Stag,            &
     &                  OCN_COMM_WORLD, Serror)
      END IF
      IF (Nexchange) THEN
!>      CALL mpi_send  (sendN, sizeN, MP_FLOAT, Ntile, Ntag,            &
!>   &                  OCN_COMM_WORLD, Nerror)
!>
        CALL mpi_send  (recvN, sizeN, MP_FLOAT, Ntile, Ntag,            &
     &                  OCN_COMM_WORLD, Nerror)
      END IF
#  endif
!
!  Adjoint of packing tile boundary data including ghost-points.
!
      IF (Sexchange) THEN
#  ifdef MPI
        CALL mpi_wait (Srequest, status(1,2), Serror)
        IF (Serror.ne.MPI_SUCCESS) THEN
          CALL mpi_error_string (Serror, string, Lstr, Ierror)
          Lstr=LEN_TRIM(string)
          WRITE (stdout,20) 'MPI_SEND/MPI_IRECV (Southern Edge)',       &
     &                      MyRank, Serror, string(1:Lstr)
 20       FORMAT (/,' AD_MP_EXCHANGE3D_BRY - error during ',a,' call,', &
     &            ' Node = ', i3.3,' Error = ',i3,/,18x,a)
          exit_flag=2
          RETURN
        END IF
#  endif
        DO m=1,GsendS
          mc=(m-1)*Klen
          j=Jstr+m-1
          DO k=LBk,UBk
            ikS=1+(k-LBk)+mc
!>          sendS(ikS)=A(j,k)
!>
            ad_A(j,k)=ad_A(j,k)+sendS(ikS)
            sendS(ikS)=0.0_r8
          END DO
        END DO
        IF (PRESENT(ad_B)) THEN
          ioff=ikS
          DO m=1,GsendS
            mc=(m-1)*Klen
            j=Jstr+m-1
            DO k=LBk,UBk
              ikS=ioff+1+(k-LBk)+mc
!>            sendS(ikS)=B(j,k)
!>
              ad_B(j,k)=ad_B(j,k)+sendS(ikS)
              sendS(ikS)=0.0_r8
            END DO
          END DO
        END IF
        IF (PRESENT(ad_C)) THEN
          ioff=ikS
          DO m=1,GsendS
            mc=(m-1)*Klen
            j=Jstr+m-1
            DO k=LBk,UBk
              ikS=ioff+1+(k-LBk)+mc
!>            sendS(ikS)=C(j,k)
!>
              ad_C(j,k)=ad_C(j,k)+sendS(ikS)
              sendS(ikS)=0.0_r8
            END DO
          END DO
        END IF
        IF (PRESENT(ad_D)) THEN
          ioff=ikS
          DO m=1,GsendS
            mc=(m-1)*Klen
            j=Jstr+m-1
            DO k=LBk,UBk
              ikS=ioff+1+(k-LBk)+mc
!>            sendS(ikS)=D(j,k)
!>
              ad_D(j,k)=ad_D(j,k)+sendS(ikS)
              sendS(ikS)=0.0_r8
            END DO
          END DO
        END IF
      END IF
!
      IF (Nexchange) THEN
#  ifdef MPI
        CALL mpi_wait (Nrequest, status(1,4), Nerror)
        IF (Nerror.ne.MPI_SUCCESS) THEN
          CALL mpi_error_string (Nerror, string, Lstr, Ierror)
          Lstr=LEN_TRIM(string)
          WRITE (stdout,20) 'MPI_SEND/MPI_IRECV (Northern Edge)',       &
     &                      MyRank, Nerror, string(1:Lstr)
          exit_flag=2
          RETURN
        END IF
#  endif
        DO m=1,GsendN
          mc=(m-1)*Klen
          j=Jend-GsendN+m
          DO k=LBk,UBk
            ikN=1+(k-LBk)+mc
!>          sendN(ikN)=A(j,k)
!>
            ad_A(j,k)=ad_A(j,k)+sendN(ikN)
            sendN(ikN)=0.0_r8
          END DO
        END DO
        IF (PRESENT(ad_B)) THEN
          ioff=ikN
          DO m=1,GsendN
            mc=(m-1)*Klen
            j=Jend-GsendN+m
            DO k=LBk,UBk
              ikN=ioff+1+(k-LBk)+mc
!>            sendN(ikN)=B(j,k)
!>
              ad_B(j,k)=ad_B(j,k)+sendN(ikN)
              sendN(ikN)=0.0_r8
            END DO
          END DO
        END IF
        IF (PRESENT(ad_C)) THEN
          ioff=ikN
          DO m=1,GsendN
            mc=(m-1)*Klen
            j=Jend-GsendN+m
            DO k=LBk,UBk
              ikN=ioff+1+(k-LBk)+mc
!>            sendN(ikN)=C(j,k)
!>
              ad_C(j,k)=ad_C(j,k)+sendN(ikN)
              sendN(ikN)=0.0_r8
            END DO
          END DO
        END IF
        IF (PRESENT(ad_D)) THEN
          ioff=ikN
          DO m=1,GsendN
            mc=(m-1)*Klen
            j=Jend-GsendN+m
            DO k=LBk,UBk
              ikN=ioff+1+(k-LBk)+mc
!>            sendN(ikN)=D(j,k)
!>
              ad_D(j,k)=ad_D(j,k)+sendN(ikN)
              sendN(ikN)=0.0_r8
            END DO
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Adjoint of unpack Eastern and Western segments.
!-----------------------------------------------------------------------
!
      IF (Eexchange) THEN
        DO i=1,BufferSizeEW
          recvE(i)=0.0_r8
          sendE(i)=0.0_r8
        END DO
        sizeE=0
        DO m=1,GrecvE
          mc=(m-1)*Klen
          i=Iend+m
          DO k=LBk,UBk
            sizeE=sizeE+1
            jkE=1+(k-LBk)+mc
!>          A(i,k)=recvE(jkE)
!>
            recvE(jkE)=ad_A(i,k)
            ad_A(i,k)=0.0_r8
          ENDDO
        END DO
        IF (PRESENT(ad_B)) THEN
          joff=jkE
          DO m=1,GrecvE
            mc=(m-1)*Klen
            i=Iend+m
            DO k=LBk,UBk
              sizeE=sizeE+1
              jkE=joff+1+(k-LBk)+mc
!>            B(i,k)=recvE(jkE)
!>
              recvE(jkE)=ad_B(i,k)
              ad_B(i,k)=0.0_r8
            END DO
          END DO
        END IF
        IF (PRESENT(ad_C)) THEN
          joff=jkE
          DO m=1,GrecvE
            mc=(m-1)*Klen
            i=Iend+m
            DO k=LBk,UBk
              sizeE=sizeE+1
              jkE=joff+1+(k-LBk)+mc
!>            C(i,k)=recvE(jkE)
!>
              recvE(jkE)=ad_C(i,k)
              ad_C(i,k)=0.0_r8
            END DO
          END DO
        END IF
        IF (PRESENT(ad_D)) THEN
          joff=jkE
          DO m=1,GrecvE
            mc=(m-1)*Klen
            i=Iend+m
            DO k=LBk,UBk
              sizeE=sizeE+1
              jkE=joff+1+(k-LBk)+mc
!>            D(i,k)=recvE(jkE)
!>
              recvE(jkE)=ad_D(i,k)
              ad_D(i,k)=0.0_r8
            END DO
          END DO
        END IF
      END IF
!
      IF (Wexchange) THEN
        DO i=1,BufferSizeEW
          recvW(i)=0.0_r8
          sendW(i)=0.0_r8
        END DO
        sizeW=0
        DO m=GrecvW,1,-1
          mc=(GrecvW-m)*Klen
          i=Istr-m
          DO k=LBk,UBk
            sizeW=sizeW+1
            jkW=1+(k-LBk)+mc
!>          A(i,k)=recvW(jkW)
!>
            recvW(jkW)=ad_A(i,k)
            ad_A(i,k)=0.0_r8
          END DO
        END DO
        IF (PRESENT(ad_B)) THEN
          joff=jkW
          DO m=GrecvW,1,-1
            mc=(GrecvW-m)*Klen
            i=Istr-m
            DO k=LBk,UBk
              sizeW=sizeW+1
              jkW=joff+1+(k-LBk)+mc
!>            B(i,k)=recvW(jkW)
!>
              recvW(jkW)=ad_B(i,k)
              ad_B(i,k)=0.0_r8
            END DO
          END DO
        END IF
        IF (PRESENT(ad_C)) THEN
          joff=jkW
          DO m=GrecvW,1,-1
            mc=(GrecvW-m)*Klen
            i=Istr-m
            DO k=LBk,UBk
              sizeW=sizeW+1
              jkW=joff+1+(k-LBk)+mc
!>            C(i,k)=recvW(jkW)
!>
              recvW(jkW)=ad_C(i,k)
              ad_C(i,k)=0.0_r8
            END DO
          END DO
        END IF
        IF (PRESENT(ad_D)) THEN
          joff=jkW
          DO m=GrecvW,1,-1
            mc=(GrecvW-m)*Klen
            i=Istr-m
            DO k=LBk,UBk
              sizeW=sizeW+1
              jkW=joff+1+(k-LBk)+mc
!>            D(i,k)=recvW(jkW)
!>
              recvW(jkW)=ad_D(i,k)
              ad_D(i,k)=0.0_r8
            END DO
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Send and receive Western and Eastern segments.
!-----------------------------------------------------------------------
!
#  if defined MPI
      IF (Wexchange) THEN
!>      CALL mpi_irecv (recvW, EWsize, MP_FLOAT, Wtile, Etag,           &
!>   &                  OCN_COMM_WORLD, Wrequest, Werror)
!>
        CALL mpi_irecv (sendW, EWsize, MP_FLOAT, Wtile, Etag,           &
     &                  OCN_COMM_WORLD, Wrequest, Werror)
      END IF
      IF (Eexchange) THEN
!>      CALL mpi_irecv (recvE, EWsize, MP_FLOAT, Etile, Wtag,           &
!>   &                  OCN_COMM_WORLD, Erequest, Eerror)
!>
        CALL mpi_irecv (sendE, EWsize, MP_FLOAT, Etile, Wtag,           &
     &                  OCN_COMM_WORLD, Erequest, Eerror)
      END IF
      IF (Wexchange) THEN
!>      CALL mpi_send  (sendW, sizeW, MP_FLOAT, Wtile, Wtag,            &
!>   &                  OCN_COMM_WORLD, Werror)
!>
        CALL mpi_send  (recvW, sizeW, MP_FLOAT, Wtile, Wtag,            &
     &                  OCN_COMM_WORLD, Werror)
      END IF
      IF (Eexchange) THEN
!>      CALL mpi_send  (sendE, sizeE, MP_FLOAT, Etile, Etag,            &
!>   &                  OCN_COMM_WORLD, Eerror)
!>
        CALL mpi_send  (recvE, sizeE, MP_FLOAT, Etile, Etag,            &
     &                  OCN_COMM_WORLD, Eerror)
      END IF
#  endif
!
!  Adjoint of packing tile boundary data including ghost-points.
!
      IF (Wexchange) THEN
#  ifdef MPI
        CALL mpi_wait (Wrequest, status(1,1), Werror)
        IF (Werror.ne.MPI_SUCCESS) THEN
          CALL mpi_error_string (Werror, string, Lstr, Ierror)
          Lstr=LEN_TRIM(string)
          WRITE (stdout,20) 'MPI_SEND/MPI_IRECV (Western Edge)',        &
     &                      MyRank, Werror, string(1:Lstr)
          exit_flag=2
          RETURN
        END IF
#  endif
        DO m=1,GsendW
          mc=(m-1)*Klen
          i=Istr+m-1
          DO k=LBk,UBk
            jkW=1+(k-LBk)+mc
!>          sendW(jkW)=A(i,k)
!>
            ad_A(i,k)=ad_A(i,k)+sendW(jkW)
            sendW(jkW)=0.0_r8
          END DO
        END DO
        IF (PRESENT(ad_B)) THEN
          joff=jkW
          DO m=1,GsendW
            mc=(m-1)*Klen
            i=Istr+m-1
            DO k=LBk,UBk
              jkW=joff+1+(k-LBk)+mc
!>            sendW(jkW)=B(i,k)
!>
              ad_B(i,k)=ad_B(i,k)+sendW(jkW)
              sendW(jkW)=0.0_r8
            END DO
          END DO
        END IF
        IF (PRESENT(ad_C)) THEN
          joff=jkW
          DO m=1,GsendW
            mc=(m-1)*Klen
            i=Istr+m-1
            DO k=LBk,UBk
              jkW=joff+1+(k-LBk)+mc
!>            sendW(jkW)=C(i,k)
!>
              ad_C(i,k)=ad_C(i,k)+sendW(jkW)
              sendW(jkW)=0.0_r8
            END DO
          END DO
        END IF
        IF (PRESENT(ad_D)) THEN
          joff=jkW
          DO m=1,GsendW
            mc=(m-1)*Klen
            i=Istr+m-1
            DO k=LBk,UBk
              jkW=joff+1+(k-LBk)+mc
!>            sendW(jkW)=D(i,k)
!>
              ad_D(i,k)=ad_D(i,k)+sendW(jkW)
              sendW(jkW)=0.0_r8
            END DO
          END DO
        END IF
      END IF
!
      IF (Eexchange) THEN
#  ifdef MPI
        CALL mpi_wait (Erequest, status(1,3), Eerror)
        IF (Eerror.ne.MPI_SUCCESS) THEN
          CALL mpi_error_string (Eerror, string, Lstr, Ierror)
          Lstr=LEN_TRIM(string)
          WRITE (stdout,20) 'MPI_SEND/MPI_IRECV (Eastern Edge)',        &
     &                      MyRank, Eerror, string(1:Lstr)
          exit_flag=2
          RETURN
        END IF
#  endif
        DO m=1,GsendE
          mc=(m-1)*Klen
          i=Iend-GsendE+m
          DO k=LBk,UBk
            jkE=1+(k-LBk)+mc
!>          sendE(jkE)=A(i,k)
!>
            ad_A(i,k)=ad_A(i,k)+sendE(jkE)
            sendE(jkE)=0.0_r8
          END DO
        END DO
        IF (PRESENT(ad_B)) THEN
          joff=jkE
          DO m=1,GsendE
            mc=(m-1)*Klen
            i=Iend-GsendE+m
            DO k=LBk,UBk
              jkE=joff+1+(k-LBk)+mc
!>            sendE(jkE)=B(i,k)
!>
              ad_B(i,k)=ad_B(i,k)+sendE(jkE)
              sendE(jkE)=0.0_r8
            END DO
          END DO
        END IF
        IF (PRESENT(ad_C)) THEN
          joff=jkE
          DO m=1,GsendE
            mc=(m-1)*Klen
            i=Iend-GsendE+m
            DO k=LBk,UBk
              jkE=joff+1+(k-LBk)+mc
!>            sendE(jkE)=C(i,k)
!>
              ad_C(i,k)=ad_C(i,k)+sendE(jkE)
              sendE(jkE)=0.0_r8
            END DO
          END DO
        END IF
        IF (PRESENT(ad_D)) THEN
          joff=jkE
          DO m=1,GsendE
            mc=(m-1)*Klen
            i=Iend-GsendE+m
            DO k=LBk,UBk
              jkE=joff+1+(k-LBk)+mc
!>            sendE(jkE)=D(i,k)
!>
              ad_D(i,k)=ad_D(i,k)+sendE(jkE)
              sendE(jkE)=0.0_r8
            END DO
          END DO
        END IF
      END IF

#  ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn off time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, model, 43)
#  endif

      RETURN
      END SUBROUTINE ad_mp_exchange3d_bry
!
!***********************************************************************
      SUBROUTINE ad_mp_exchange4d (ng, tile, model, Nvar,               &
     &                             LBi, UBi, LBj, UBj, LBk, UBk,        &
     &                             LBt, UBt,                            &
     &                             Nghost, EW_periodic, NS_periodic,    &
     &                             ad_A, ad_B, ad_C)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
      implicit none
!
!  Imported variable declarations.
!
      logical, intent(in) :: EW_periodic, NS_periodic

      integer, intent(in) :: ng, tile, model, Nvar
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBk, UBk, LBt, UBt
      integer, intent(in) :: Nghost
!
#  ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: ad_A(LBi:,LBj:,LBk:,LBt:)

      real(r8), intent(inout), optional :: ad_B(LBi:,LBj:,LBk:,LBt:)
      real(r8), intent(inout), optional :: ad_C(LBi:,LBj:,LBk:,LBt:)
#  else
      real(r8), intent(inout) :: ad_A(LBi:UBi,LBj:UBj,LBk:UBk,LBt:UBt)

      real(r8), intent(inout), optional ::                              &
     &                           ad_B(LBi:UBi,LBj:UBj,LBk:UBk,LBt:UBt)
      real(r8), intent(inout), optional ::                              &
     &                           ad_C(LBi:UBi,LBj:UBj,LBk:UBk,LBt:UBt)
#  endif
!
!  Local variable declarations.
!
      logical :: Wexchange, Sexchange, Eexchange, Nexchange

      integer :: i, ikS, ikN, ioff, Imin, Imax, Ilen, IKlen, IKTlen
      integer :: j, jkW, jkE, joff, Jmin, Jmax, Jlen, JKlen, JKTlen
      integer :: k, kc, m, mc, Ierror, Klen, Lstr, Tlen, pp
      integer :: l, lc
      integer :: Wtile, GsendW, GrecvW, Wtag, Werror, Wrequest
      integer :: Stile, GsendS, GrecvS, Stag, Serror, Srequest
      integer :: Etile, GsendE, GrecvE, Etag, Eerror, Erequest
      integer :: Ntile, GsendN, GrecvN, Ntag, Nerror, Nrequest
      integer :: BufferSizeEW, EWsize, sizeW, sizeE
      integer :: BufferSizeNS, NSsize, sizeS, sizeN

#  ifdef MPI
      integer, dimension(MPI_STATUS_SIZE,4) :: status
#  endif
!
      real(r8), dimension(Nvar*HaloSizeJ(ng)*(UBk-LBk+1)*               &
     &                                      (UBt-LBt+1)) :: sendW, sendE
      real(r8), dimension(Nvar*HaloSizeI(ng)*(UBk-LBk+1)*               &
     &                                      (UBt-LBt+1)) :: sendS, sendN

      real(r8), dimension(Nvar*HaloSizeJ(ng)*(UBk-LBk+1)*               &
     &                                      (UBt-LBt+1)) :: recvW, recvE
      real(r8), dimension(Nvar*HaloSizeI(ng)*(UBk-LBk+1)*               &
     &                                      (UBt-LBt+1)) :: recvS, recvN

      character (len=MPI_MAX_ERROR_STRING) :: string

#  include "set_bounds.h"
#  ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn on time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, model, 42)
#  endif
!
!-----------------------------------------------------------------------
!  Determine rank of tile neighbors and number of ghost-points to
!  exchange.
!-----------------------------------------------------------------------
!
      CALL tile_neighbors (ng, Nghost, EW_periodic, NS_periodic,        &
     &                     GrecvW, GsendW, Wtile, Wexchange,            &
     &                     GrecvE, GsendE, Etile, Eexchange,            &
     &                     GrecvS, GsendS, Stile, Sexchange,            &
     &                     GrecvN, GsendN, Ntile, Nexchange)
!
!  Set communication tags.
!
      Wtag=1
      Stag=2
      Etag=3
      Ntag=4
!
!  Determine range and length of the distributed tile boundary segments.
!
      Imin=LBi
      Imax=UBi
      Jmin=LBj
      Jmax=UBj
      Ilen=Imax-Imin+1
      Jlen=Jmax-Jmin+1
      Klen=UBk-LBk+1
      Tlen=UBt-LBt+1
      IKlen=Ilen*Klen
      JKlen=Jlen*Klen
      IKTlen=IKlen*Tlen
      JKTlen=JKlen*Tlen
      IF (EW_periodic.or.NS_periodic) THEN
        pp=1
      ELSE
        pp=0
      END IF
      NSsize=Nvar*(Nghost+pp)*IKTlen
      EWsize=Nvar*(Nghost+pp)*JKTlen
      BufferSizeNS=Nvar*HaloSizeI(ng)*Klen*Tlen
      BufferSizeEW=Nvar*HaloSizeJ(ng)*Klen*Tlen
      IF (SIZE(sendE).lt.EWsize) THEN
        WRITE (stdout,10) 'EWsize = ', EWsize, SIZE(sendE)
 10     FORMAT (/,' AD_MP_EXCHANGE4D - communication buffer too',       &
     &          ' small, ',a, 2i8)
      END IF
      IF (SIZE(sendN).lt.NSsize) THEN
        WRITE (stdout,10) 'NSsize = ', NSsize, SIZE(sendN)
      END IF
!
!-----------------------------------------------------------------------
!  Adjoint of unpacking Northern and Southern segments.
!-----------------------------------------------------------------------
!
      IF (Nexchange) THEN
        DO i=1,BufferSizeNS
          recvN(i)=0.0_r8
          sendN(i)=0.0_r8
        END DO
        sizeN=0
        DO m=1,GrecvN
          mc=(m-1)*IKTlen
          j=Jend+m
          DO l=LBt,UBt
            lc=(l-LBt)*IKlen+mc
            DO k=LBk,UBk
              kc=(k-LBk)*Ilen+lc
              DO i=Imin,Imax
                sizeN=sizeN+1
                ikN=1+(i-Imin)+kc
!>              A(i,j,k,l)=recvN(ikN)
!>
                recvN(ikN)=ad_A(i,j,k,l)
                ad_A(i,j,k,l)=0.0_r8
              END DO
            END DO
          END DO
        END DO
        IF (PRESENT(ad_B)) THEN
          ioff=ikN
          DO m=1,GrecvN
            mc=(m-1)*IKTlen
            j=Jend+m
            DO l=LBt,UBt
              lc=(l-LBt)*IKlen+mc
              DO k=LBk,UBk
                kc=(k-LBk)*Ilen+lc
                DO i=Imin,Imax
                  sizeN=sizeN+1
                  ikN=ioff+1+(i-Imin)+kc
!>                B(i,j,k,l)=recvN(ikN)
!>
                  recvN(ikN)=ad_B(i,j,k,l)
                  ad_B(i,j,k,l)=0.0_r8
                END DO
              END DO
            END DO
          END DO
        END IF
        IF (PRESENT(ad_C)) THEN
          ioff=ikN
          DO m=1,GrecvN
            mc=(m-1)*IKTlen
            j=Jend+m
            DO l=LBt,UBt
              lc=(l-LBt)*IKlen+mc
              DO k=LBk,UBk
                kc=(k-LBk)*Ilen+lc
                DO i=Imin,Imax
                  sizeN=sizeN+1
                  ikN=ioff+1+(i-Imin)+kc
!>                C(i,j,k,l)=recvN(ikN)
!>
                  recvN(ikN)=ad_C(i,j,k,l)
                  ad_C(i,j,k,l)=0.0_r8
                END DO
              END DO
            END DO
          END DO
        END IF
      END IF
!
      IF (Sexchange) THEN
        DO i=1,BufferSizeNS
          recvS(i)=0.0_r8
          sendS(i)=0.0_r8
        END DO
        sizeS=0
        DO m=GrecvS,1,-1
          mc=(GrecvS-m)*IKTlen
          j=Jstr-m
          DO l=LBt,UBt
            lc=(l-LBt)*IKlen+mc
            DO k=LBk,UBk
              kc=(k-LBk)*Ilen+lc
              DO i=Imin,Imax
                sizeS=sizeS+1
                ikS=1+(i-Imin)+kc
!>              A(i,j,k,l)=recvS(ikS)
!>
                recvS(ikS)=ad_A(i,j,k,l)
                ad_A(i,j,k,l)=0.0_r8
              END DO
            END DO
          END DO
        END DO
        IF (PRESENT(ad_B)) THEN
          ioff=ikS
          DO m=GrecvS,1,-1
            mc=(GrecvS-m)*IKTlen
            j=Jstr-m
            DO l=LBt,UBt
              lc=(l-LBt)*IKlen+mc
              DO k=LBk,UBk
                kc=(k-LBk)*Ilen+lc
                DO i=Imin,Imax
                  sizeS=sizeS+1
                  ikS=ioff+1+(i-Imin)+kc
!>                B(i,Jstr-m,k,l)=recvS(ikS)
!>
                  recvS(ikS)=ad_B(i,j,k,l)
                  ad_B(i,j,k,l)=0.0_r8
                END DO
              END DO
            END DO
          END DO
        END IF
        IF (PRESENT(ad_C)) THEN
          ioff=ikS
          DO m=GrecvS,1,-1
            mc=(GrecvS-m)*IKTlen
            j=Jstr-m
            DO l=LBt,UBt
              lc=(l-LBt)*IKlen+mc
              DO k=LBk,UBk
                kc=(k-LBk)*Ilen+lc
                DO i=Imin,Imax
                  sizeS=sizeS+1
                  ikS=ioff+1+(i-Imin)+kc
!>                C(i,Jstr-m,k,l)=recvS(ikS)
!>
                  recvS(ikS)=ad_C(i,j,k,l)
                  ad_C(i,j,k,l)=0.0_r8
                END DO
              END DO
            END DO
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Adjoint of send and receive Southern and Northern segments.
!-----------------------------------------------------------------------
!
#  if defined MPI
      IF (Sexchange) THEN
!>      CALL mpi_irecv (recvS, NSsize, MP_FLOAT, Stile, Ntag,           &
!>   &                  OCN_COMM_WORLD, Srequest, Serror)
!>
        CALL mpi_irecv (sendS, NSsize, MP_FLOAT, Stile, Ntag,           &
     &                  OCN_COMM_WORLD, Srequest, Serror)
      END IF
      IF (Nexchange) THEN
!>      CALL mpi_irecv (recvN, NSsize, MP_FLOAT, Ntile, Stag,           &
!>   &                  OCN_COMM_WORLD, Nrequest, Nerror)
!>
        CALL mpi_irecv (sendN, NSsize, MP_FLOAT, Ntile, Stag,           &
     &                  OCN_COMM_WORLD, Nrequest, Nerror)
      END IF
      IF (Sexchange) THEN
!>      CALL mpi_send  (sendS, sizeS, MP_FLOAT, Stile, Stag,            &
!>   &                  OCN_COMM_WORLD, Serror)
!>
        CALL mpi_send  (recvS, sizeS, MP_FLOAT, Stile, Stag,            &
     &                  OCN_COMM_WORLD, Serror)
      END IF
      IF (Nexchange) THEN
!>      CALL mpi_send  (sendN, sizeN, MP_FLOAT, Ntile, Ntag,            &
!>   &                  OCN_COMM_WORLD, Nerror)
!>
        CALL mpi_send  (recvN, sizeN, MP_FLOAT, Ntile, Ntag,            &
     &                  OCN_COMM_WORLD, Nerror)
      END IF
#  endif
!
!  Adjoint of packing tile boundary data including ghost-points.
!
      IF (Sexchange) THEN
#  ifdef MPI
        CALL mpi_wait (Srequest, status(1,2), Serror)
        IF (Serror.ne.MPI_SUCCESS) THEN
          CALL mpi_error_string (Serror, string, Lstr, Ierror)
          Lstr=LEN_TRIM(string)
          WRITE (stdout,20) 'MPI_SEND/MPI_IRECV (Southern Edge)',       &
     &                      MyRank, Serror, string(1:Lstr)
 20       FORMAT (/,' AD_MP_EXCHANGE4D - error during ',a,' call,',     &
     &            ' Node = ', i3.3,' Error = ',i3,/,18x,a)
          exit_flag=2
          RETURN
        END IF
#  endif
        DO m=1,GsendS
          mc=(m-1)*IKTlen
          j=Jstr+m-1
          DO l=LBt,UBt
            lc=(l-LBt)*IKlen+mc
            DO k=LBk,UBk
              kc=(k-LBk)*Ilen+lc
              DO i=Imin,Imax
                ikS=1+(i-Imin)+kc
!>              sendS(ikS)=A(i,j,k,l)
!>
                ad_A(i,j,k,l)=ad_A(i,j,k,l)+sendS(ikS)
                sendS(ikS)=0.0_r8
              END DO
            END DO
          END DO
        END DO
        IF (PRESENT(ad_B)) THEN
          ioff=ikS
          DO m=1,GsendS
            mc=(m-1)*IKTlen
            j=Jstr+m-1
            DO l=LBt,UBt
              lc=(l-LBt)*IKlen+mc
              DO k=LBk,UBk
                kc=(k-LBk)*Ilen+lc
                DO i=Imin,Imax
                  ikS=ioff+1+(i-Imin)+kc
!>                sendS(ikS)=B(i,j,k,l)
!>
                  ad_B(i,j,k,l)=ad_B(i,j,k,l)+sendS(ikS)
                  sendS(ikS)=0.0_r8
                END DO
              END DO
            END DO
          END DO
        END IF
        IF (PRESENT(ad_C)) THEN
          ioff=ikS
          DO m=1,GsendS
            mc=(m-1)*IKTlen
            j=Jstr+m-1
            DO l=LBt,UBt
              lc=(l-LBt)*IKlen+mc
              DO k=LBk,UBk
                kc=(k-LBk)*Ilen+lc
                DO i=Imin,Imax
                  ikS=ioff+1+(i-Imin)+kc
!>                sendS(ikS)=C(i,j,k,l)
!>
                  ad_C(i,j,k,l)=ad_C(i,j,k,l)+sendS(ikS)
                  sendS(ikS)=0.0_r8
                END DO
              END DO
            END DO
          END DO
        END IF
      END IF
!
      IF (Nexchange) THEN
#  ifdef MPI
        CALL mpi_wait (Nrequest, status(1,4), Nerror)
        IF (Nerror.ne.MPI_SUCCESS) THEN
          CALL mpi_error_string (Nerror, string, Lstr, Ierror)
          Lstr=LEN_TRIM(string)
          WRITE (stdout,20) 'MPI_SEND/MPI_IRECV (Northern Edge)',       &
     &                      MyRank, Nerror, string(1:Lstr)
          exit_flag=2
          RETURN
        END IF
#  endif
        DO m=1,GsendN
          mc=(m-1)*IKTlen
          j=Jend-GsendN+m
          DO l=LBt,UBt
            lc=(l-LBt)*IKlen+mc
            DO k=LBk,UBk
              kc=(k-LBk)*Ilen+lc
              DO i=Imin,Imax
                ikN=1+(i-Imin)+kc
!>              sendN(ikN)=A(i,j,k,l)
!>
                ad_A(i,j,k,l)=ad_A(i,j,k,l)+sendN(ikN)
                sendN(ikN)=0.0_r8
              END DO
            END DO
          END DO
        END DO
        IF (PRESENT(ad_B)) THEN
          ioff=ikN
          DO m=1,GsendN
            mc=(m-1)*IKlen
            j=Jend-GsendN+m
            DO l=LBt,UBt
              lc=(l-LBt)*IKlen+mc
              DO k=LBk,UBk
                kc=(k-LBk)*Ilen+lc
                DO i=Imin,Imax
                  ikN=ioff+1+(i-Imin)+kc
!>                sendN(ikN)=B(i,j,k,l)
!>
                  ad_B(i,j,k,l)=ad_B(i,j,k,l)+sendN(ikN)
                  sendN(ikN)=0.0_r8
                END DO
              END DO
            END DO
          END DO
        END IF
        IF (PRESENT(ad_C)) THEN
          ioff=ikN
          DO m=1,GsendN
            mc=(m-1)*IKlen
            j=Jend-GsendN+m
            DO l=LBt,UBt
              lc=(l-LBt)*IKlen+mc
              DO k=LBk,UBk
                kc=(k-LBk)*Ilen+lc
                DO i=Imin,Imax
                  ikN=ioff+1+(i-Imin)+kc
!>                sendN(ikN)=C(i,j,k,l)
!>
                  ad_C(i,j,k,l)=ad_C(i,j,k,l)+sendN(ikN)
                  sendN(ikN)=0.0_r8
                END DO
              END DO
            END DO
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Adjoint of unpack Eastern and Western segments.
!-----------------------------------------------------------------------
!
      IF (Eexchange) THEN
        DO i=1,BufferSizeEW
          recvE(i)=0.0_r8
          sendE(i)=0.0_r8
        END DO
        sizeE=0
        DO m=1,GrecvE
          mc=(m-1)*JKTlen
          i=Iend+m
          DO l=LBt,UBt
            lc=(l-LBt)*JKlen+mc
            DO k=LBk,UBk
              kc=(k-LBk)*Jlen+lc
              DO j=Jmin,Jmax
                sizeE=sizeE+1
                jkE=1+(j-Jmin)+kc
!>              A(i,j,k,l)=recvE(jkE)
!>
                recvE(jkE)=ad_A(i,j,k,l)
                ad_A(i,j,k,l)=0.0_r8
              END DO
            END DO
          ENDDO
        END DO
        IF (PRESENT(ad_B)) THEN
          joff=jkE
          DO m=1,GrecvE
            mc=(m-1)*JKTlen
            i=Iend+m
            DO l=LBt,UBt
              lc=(l-LBt)*JKlen+mc
              DO k=LBk,UBk
                kc=(k-LBk)*Jlen+lc
                DO j=Jmin,Jmax
                  sizeE=sizeE+1
                  jkE=joff+1+(j-Jmin)+kc
!>                B(i,j,k,l)=recvE(jkE)
!>
                  recvE(jkE)=ad_B(i,j,k,l)
                  ad_B(i,j,k,l)=0.0_r8
                END DO
              END DO
            END DO
          END DO
        END IF
        IF (PRESENT(ad_C)) THEN
          joff=jkE
          DO m=1,GrecvE
            mc=(m-1)*JKTlen
            i=Iend+m
            DO l=LBt,UBt
              lc=(l-LBt)*JKlen+mc
              DO k=LBk,UBk
                kc=(k-LBk)*Jlen+lc
                DO j=Jmin,Jmax
                  sizeE=sizeE+1
                  jkE=joff+1+(j-Jmin)+kc
!>                C(i,j,k,l)=recvE(jkE)
!>
                  recvE(jkE)=ad_C(i,j,k,l)
                  ad_C(i,j,k,l)=0.0_r8
                END DO
              END DO
            END DO
          END DO
        END IF
      END IF
!
      IF (Wexchange) THEN
        DO i=1,BufferSizeEW
          recvW(i)=0.0_r8
          sendW(i)=0.0_r8
        END DO
        sizeW=0
        DO m=GrecvW,1,-1
          mc=(GrecvW-m)*JKTlen
          i=Istr-m
          DO l=LBt,UBt
            lc=(l-LBt)*JKlen+mc
            DO k=LBk,UBk
              kc=(k-LBk)*Jlen+lc
              DO j=Jmin,Jmax
                sizeW=sizeW+1
                jkW=1+(j-Jmin)+kc
!>              A(i,j,k,l)=recvW(jkW)
!>
                recvW(jkW)=ad_A(i,j,k,l)
                ad_A(i,j,k,l)=0.0_r8
              END DO
            END DO
          END DO
        END DO
        IF (PRESENT(ad_B)) THEN
          joff=jkW
          DO m=GrecvW,1,-1
            mc=(GrecvW-m)*JKTlen
            i=Istr-m
            DO l=LBt,UBt
              lc=(l-LBt)*JKlen+mc
              DO k=LBk,UBk
                kc=(k-LBk)*Jlen+lc
                DO j=Jmin,Jmax
                  sizeW=sizeW+1
                  jkW=joff+1+(j-Jmin)+kc
!>                B(i,j,k,l)=recvW(jkW)
!>
                  recvW(jkW)=ad_B(i,j,k,l)
                  ad_B(i,j,k,l)=0.0_r8
                END DO
              END DO
            END DO
          END DO
        END IF
        IF (PRESENT(ad_C)) THEN
          joff=jkW
          DO m=GrecvW,1,-1
            mc=(GrecvW-m)*JKTlen
            i=Istr-m
            DO l=LBt,UBt
              lc=(l-LBt)*JKlen+mc
              DO k=LBk,UBk
                kc=(k-LBk)*Jlen+lc
                DO j=Jmin,Jmax
                  sizeW=sizeW+1
                  jkW=joff+1+(j-Jmin)+kc
!>                C(i,j,k,l)=recvW(jkW)
!>
                  recvW(jkW)=ad_C(i,j,k,l)
                  ad_C(i,j,k,l)=0.0_r8
                END DO
              END DO
            END DO
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Send and receive Western and Eastern segments.
!-----------------------------------------------------------------------
!
#  if defined MPI
      IF (Wexchange) THEN
!>      CALL mpi_irecv (recvW, EWsize, MP_FLOAT, Wtile, Etag,           &
!>   &                  OCN_COMM_WORLD, Wrequest, Werror)
!>
        CALL mpi_irecv (sendW, EWsize, MP_FLOAT, Wtile, Etag,           &
     &                  OCN_COMM_WORLD, Wrequest, Werror)
      END IF
      IF (Eexchange) THEN
!>      CALL mpi_irecv (recvE, EWsize, MP_FLOAT, Etile, Wtag,           &
!>   &                  OCN_COMM_WORLD, Erequest, Eerror)
!>
        CALL mpi_irecv (sendE, EWsize, MP_FLOAT, Etile, Wtag,           &
     &                  OCN_COMM_WORLD, Erequest, Eerror)
      END IF
      IF (Wexchange) THEN
!>      CALL mpi_send  (sendW, sizeW, MP_FLOAT, Wtile, Wtag,            &
!>   &                  OCN_COMM_WORLD, Werror)
!>
        CALL mpi_send  (recvW, sizeW, MP_FLOAT, Wtile, Wtag,            &
     &                  OCN_COMM_WORLD, Werror)
      END IF
      IF (Eexchange) THEN
!>      CALL mpi_send  (sendE, sizeE, MP_FLOAT, Etile, Etag,            &
!>   &                  OCN_COMM_WORLD, Eerror)
!>
        CALL mpi_send  (recvE, sizeE, MP_FLOAT, Etile, Etag,            &
     &                  OCN_COMM_WORLD, Eerror)
      END IF
#  endif
!
!  Adjoint of packing tile boundary data including ghost-points.
!
      IF (Wexchange) THEN
#  ifdef MPI
        CALL mpi_wait (Wrequest, status(1,1), Werror)
        IF (Werror.ne.MPI_SUCCESS) THEN
          CALL mpi_error_string (Werror, string, Lstr, Ierror)
          Lstr=LEN_TRIM(string)
          WRITE (stdout,20) 'MPI_SEND/MPI_IRECV (Western Edge)',        &
     &                      MyRank, Werror, string(1:Lstr)
          exit_flag=2
          RETURN
        END IF
#  endif
        DO m=1,GsendW
          mc=(m-1)*JKTlen
          i=Istr+m-1
          DO l=LBt,UBt
            lc=(l-LBt)*JKlen+mc
            DO k=LBk,UBk
              kc=(k-LBk)*Jlen+lc
              DO j=Jmin,Jmax
                jkW=1+(j-Jmin)+kc
!>              sendW(jkW)=A(i,j,k,l)
!>
                ad_A(i,j,k,l)=ad_A(i,j,k,l)+sendW(jkW)
                sendW(jkW)=0.0_r8
              END DO
            END DO
          END DO
        END DO
        IF (PRESENT(ad_B)) THEN
          joff=jkW
          DO m=1,GsendW
            mc=(m-1)*JKTlen
            i=Istr+m-1
            DO l=LBt,UBt
              lc=(l-LBt)*JKlen+mc
              DO k=LBk,UBk
                kc=(k-LBk)*Jlen+lc
                DO j=Jmin,Jmax
                  jkW=joff+1+(j-Jmin)+kc
!>                sendW(jkW)=B(i,j,k,l)
!>
                  ad_B(i,j,k,l)=ad_B(i,j,k,l)+sendW(jkW)
                  sendW(jkW)=0.0_r8
                END DO
              END DO
            END DO
          END DO
        END IF
        IF (PRESENT(ad_C)) THEN
          joff=jkW
          DO m=1,GsendW
            mc=(m-1)*JKTlen
            i=Istr+m-1
            DO l=LBt,UBt
              lc=(l-LBt)*JKlen+mc
              DO k=LBk,UBk
                kc=(k-LBk)*Jlen+lc
                DO j=Jmin,Jmax
                  jkW=joff+1+(j-Jmin)+kc
!>                sendW(jkW)=C(i,j,k,l)
!>
                  ad_C(i,j,k,l)=ad_C(i,j,k,l)+sendW(jkW)
                  sendW(jkW)=0.0_r8
                END DO
              END DO
            END DO
          END DO
        END IF
      END IF
!
      IF (Eexchange) THEN
#  ifdef MPI
        CALL mpi_wait (Erequest, status(1,3), Eerror)
        IF (Eerror.ne.MPI_SUCCESS) THEN
          CALL mpi_error_string (Eerror, string, Lstr, Ierror)
          Lstr=LEN_TRIM(string)
          WRITE (stdout,20) 'MPI_SEND/MPI_IRECV (Eastern Edge)',        &
     &                      MyRank, Eerror, string(1:Lstr)
          exit_flag=2
          RETURN
        END IF
#  endif
        DO m=1,GsendE
          mc=(m-1)*JKTlen
          i=Iend-GsendE+m
          DO l=LBt,UBt
            lc=(l-LBt)*JKlen+mc
            DO k=LBk,UBk
              kc=(k-LBk)*Jlen+lc
              DO j=Jmin,Jmax
                jkE=1+(j-Jmin)+kc
!>              sendE(jkE)=A(i,j,k,l)
!>
                ad_A(i,j,k,l)=ad_A(i,j,k,l)+sendE(jkE)
                sendE(jkE)=0.0_r8
              END DO
            END DO
          END DO
        END DO
        IF (PRESENT(ad_B)) THEN
          joff=jkE
          DO m=1,GsendE
            mc=(m-1)*JKTlen
            i=Iend-GsendE+m
            DO l=LBt,UBt
              lc=(l-LBt)*JKlen+mc
              DO k=LBk,UBk
                kc=(k-LBk)*Jlen+lc
                DO j=Jmin,Jmax
                  jkE=joff+1+(j-Jmin)+kc
!>                sendE(jkE)=B(i,j,k,l)
!>
                  ad_B(i,j,k,l)=ad_B(i,j,k,l)+sendE(jkE)
                  sendE(jkE)=0.0_r8
                END DO
              END DO
            END DO
          END DO
        END IF
        IF (PRESENT(ad_C)) THEN
          joff=jkE
          DO m=1,GsendE
            mc=(m-1)*JKTlen
            i=Iend-GsendE+m
            DO l=LBt,UBt
              lc=(l-LBt)*JKlen+mc
              DO k=LBk,UBk
                kc=(k-LBk)*Jlen+lc
                DO j=Jmin,Jmax
                  jkE=joff+1+(j-Jmin)+kc
!>                sendE(jkE)=C(i,j,k,l)
!>
                  ad_C(i,j,k,l)=ad_C(i,j,k,l)+sendE(jkE)
                  sendE(jkE)=0.0_r8
                END DO
              END DO
            END DO
          END DO
        END IF
      END IF
#  ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn off time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, model, 42)
#  endif

      RETURN
      END SUBROUTINE ad_mp_exchange4d
# endif
#endif
      END MODULE mp_exchange_mod
