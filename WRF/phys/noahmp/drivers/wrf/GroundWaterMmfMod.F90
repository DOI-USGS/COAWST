module GroundWaterMmfMod

!!! Module to calculate lateral groundwater flow and the flux between groundwater and rivers
!!! plus the routine to update soil moisture and water table due to those two fluxes
!!! according to the Miguez-Macho & Fan groundwater scheme (Miguez-Macho et al., JGR 2007).
!!! Module written by Gonzalo Miguez-Macho , U. de Santiago de Compostela, Galicia, Spain
!!! November 2012 

! ------------------------ Code history -----------------------------------
! Original Noah-MP subroutine: module_sf_groundwater.F
! Original code: Miguez-Macho&Fan (Miguez-Macho et al 2007, Fan et al 2007)
! Note: this MMF scheme is not refactored
! It is directly adapted from the original WRF module_sf_noahmp_groundwater code
! Coder: Cenlin He (NCAR), December 2025
! -------------------------------------------------------------------------

  use NoahmpIOVarType, only : NoahmpIO_type

  implicit none

contains

  subroutine WTABLE_mmf_noahmp (NoahmpIO, NSOIL     ,XLAND    ,XICE    ,XICE_THRESHOLD  ,ISICE ,& !in
                                ISLTYP    ,SMOISEQ  ,DZS     ,WTDDT                  ,& !in
                                FDEPTH    ,AREA     ,TOPO    ,ISURBAN ,IVGTYP        ,& !in
                                RIVERCOND ,RIVERBED ,EQWTD   ,PEXP                   ,& !in
                                SMOIS     ,SH2OXY   ,SMCWTD  ,WTD  , QLAT, QRF       ,& !inout
                                DEEPRECH  ,QSPRING  ,QSLAT   ,QRFS ,QSPRINGS  ,RECH  ,& !inout
                                ids,ide, jds,jde, kds,kde,                            &
                                ims,ime, jms,jme, kms,kme,                            &
                                its,ite, jts,jte, kts,kte                             )
! ----------------------------------------------------------------------

    implicit none

    type(NoahmpIO_type), intent(in) :: NoahmpIO

    ! IN only
    INTEGER, INTENT(IN)                             :: ids,ide, jds,jde, kds,kde,  &
                                                       ims,ime, jms,jme, kms,kme,  &
                                                       its,ite, jts,jte, kts,kte
    REAL,    INTENT(IN)                             :: WTDDT
    REAL,    INTENT(IN)                             :: XICE_THRESHOLD
    INTEGER, INTENT(IN)                             :: ISICE
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(IN) :: XLAND, XICE
    INTEGER, DIMENSION(ims:ime,jms:jme), INTENT(IN) :: ISLTYP, IVGTYP
    INTEGER, INTENT(IN)                             :: NSOIL
    INTEGER, INTENT(IN)                             :: ISURBAN
    REAL,    DIMENSION(1:NSOIL),         INTENT(IN) :: DZS
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(IN) :: FDEPTH, AREA, TOPO, EQWTD,  &
                                                       PEXP, RIVERBED, RIVERCOND
    REAL,    DIMENSION(ims:ime,1:NSOIL,jms:jme), INTENT(IN) :: SMOISEQ
    ! IN and OUT 
    REAL,    DIMENSION(ims:ime,1:NSOIL,jms:jme), INTENT(INOUT) :: SMOIS, SH2OXY 
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT) :: WTD, SMCWTD, DEEPRECH,   &
                                                          QSLAT, QRFS, QSPRINGS, RECH
    ! OUT
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(OUT)   :: QRF, &  ! groundwater - river water flux
                                                          QSPRING ! water springing at the surface from groundwater convergence in the column
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(OUT)   :: QLAT
    ! LOCAL   
    INTEGER                             :: I,J,K  
    REAL                                :: BEXP,DKSAT,PSISAT,SMCMAX,SMCWLT
    REAL                                :: DELTAT,RCOND,TOTWATER,PSI,WPLUS,WMINUS, &
                                           WFLUXDEEP,WCNDDEEP,DDZ,SMCWTDMID
    REAL,    DIMENSION(0:NSOIL)         :: ZSOIL    ! depth of soil layer-bottom [m]
    REAL,    DIMENSION(1:NSOIL)         :: SMCEQ    ! equilibrium soil water  content [m3/m3]
    REAL,    DIMENSION(1:NSOIL)         :: SMC,SH2O
    INTEGER, DIMENSION(ims:ime,jms:jme) :: LANDMASK  
! ----------------------------------------------------------------------

    DELTAT   = WTDDT * 60.0 !timestep in seconds for this calculation
    ZSOIL(0) = 0.0
    ZSOIL(1) = -DZS(1)
    do K = 2, NSOIL
       ZSOIL(K) = -DZS(K) + ZSOIL(K-1)
    enddo

    where ( XLAND-1.5 < 0.0 .and. XICE < XICE_THRESHOLD .and. IVGTYP .NE. ISICE )
         LANDMASK = 1
    elsewhere
         LANDMASK = -1
    endwhere

    ! Calculate lateral flow
    QLAT = 0.0
    call LATERALFLOW(NoahmpIO,ISLTYP,WTD,QLAT,FDEPTH,TOPO,LANDMASK,DELTAT,AREA &
                     ,ids,ide,jds,jde,kds,kde                         &
                     ,ims,ime,jms,jme,kms,kme                         &
                     ,its,ite,jts,jte,kts,kte                         )


    ! compute flux from grounwater to rivers in the cell
    do J = jts, jte
       do I = its, ite
          if ( LANDMASK(I,J) > 0 ) then
             if (WTD(I,J) > RIVERBED(I,J) .and. EQWTD(I,J) > RIVERBED(I,J)) then
               RCOND = RIVERCOND(I,J) * exp(PEXP(I,J) * (WTD(I,J) - EQWTD(I,J)))
             else    
               RCOND = RIVERCOND(I,J)       
             endif
             QRF(I,J) = RCOND * (WTD(I,J) - RIVERBED(I,J)) * DELTAT / AREA(I,J)
             ! for now, dont allow it to go from river to groundwater
             QRF(I,J) = max(QRF(I,J), 0.0)
          else
             QRF(I,J) = 0.0
          endif
       enddo
    enddo

    do J = jts, jte
       do I = its, ite
          if ( LANDMASK(I,J) > 0 ) then
             BEXP   = NoahmpIO%BEXP_TABLE   (ISLTYP(I,J))
             DKSAT  = NoahmpIO%DKSAT_TABLE  (ISLTYP(I,J))
             PSISAT = -1.0 * NoahmpIO%PSISAT_TABLE (ISLTYP(I,J))
             SMCMAX = NoahmpIO%SMCMAX_TABLE (ISLTYP(I,J))
             SMCWLT = NoahmpIO%SMCWLT_TABLE (ISLTYP(I,J))

             if ( IVGTYP(I,J) == ISURBAN ) then
                 SMCMAX = 0.45
                 SMCWLT = 0.40
             endif

             ! for deep water table calculate recharge
             if ( WTD(I,J) < ZSOIL(NSOIL)-DZS(NSOIL) ) then
                ! assume all liquid if the wtd is deep
                DDZ       = ZSOIL(NSOIL) - WTD(I,J)
                SMCWTDMID = 0.5 * (SMCWTD(I,J) + SMCMAX)
                PSI       = PSISAT * (SMCMAX / SMCWTD(I,J)) ** BEXP
                WCNDDEEP  = DKSAT * (SMCWTDMID / SMCMAX) ** (2.0*BEXP + 3.0)
                WFLUXDEEP =  - DELTAT * WCNDDEEP * ((PSISAT-PSI) / DDZ - 1.0)
                ! update deep soil moisture
                SMCWTD(I,J)   = SMCWTD(I,J) + (DEEPRECH(I,J) - WFLUXDEEP) / DDZ
                WPLUS         = max((SMCWTD(I,J)-SMCMAX), 0.0) * DDZ
                WMINUS        = max((1.0e-4 - SMCWTD(I,J)), 0.0) * DDZ
                SMCWTD(I,J)   = max(min(SMCWTD(I,J),SMCMAX), 1.0e-4)
                WFLUXDEEP     = WFLUXDEEP + WPLUS - WMINUS
                DEEPRECH(I,J) = WFLUXDEEP
             endif

             ! Total water flux to or from groundwater in the cell
             TOTWATER       = QLAT(I,J) - QRF(I,J) + DEEPRECH(I,J)
             SMC(1:NSOIL)   = SMOIS(I,1:NSOIL,J)
             SH2O(1:NSOIL)  = SH2OXY(I,1:NSOIL,J)
             SMCEQ(1:NSOIL) = SMOISEQ(I,1:NSOIL,J)

             ! Update the water table depth and soil moisture
             call UPDATEWTD (NSOIL, DZS , ZSOIL, SMCEQ, SMCMAX, SMCWLT, PSISAT, BEXP, I, J, & ! in
                             TOTWATER, WTD(I,J), SMC, SH2O, SMCWTD(I,J),                    & ! inout
                             QSPRING(I,J) )                                                   ! out

             ! now update soil moisture
             SMOIS(I,1:NSOIL,J)  = SMC(1:NSOIL)
             SH2OXY(I,1:NSOIL,J) = SH2O(1:NSOIL)
          endif
       enddo
    enddo

    ! accumulate fluxes for output
    do J = jts, jte
       do I = its, ite
          QSLAT(I,J)    = QSLAT(I,J) + QLAT(I,J) * 1.0e3
          QRFS(I,J)     = QRFS(I,J) + QRF(I,J) * 1.0e3
          QSPRINGS(I,J) = QSPRINGS(I,J) + QSPRING(I,J) * 1.0e3
          RECH(I,J)     = RECH(I,J) + DEEPRECH(I,J) * 1.0e3
          ! zero out DEEPRECH
          DEEPRECH(I,J) =0.0
       enddo
    enddo

  end subroutine WTABLE_mmf_noahmp

! ==================================================================================================
  subroutine LATERALFLOW (NoahmpIO, ISLTYP,WTD,QLAT,FDEPTH,TOPO,LANDMASK,DELTAT,AREA &
                          ,ids,ide,jds,jde,kds,kde                         &
                          ,ims,ime,jms,jme,kms,kme                         &
                          ,its,ite,jts,jte,kts,kte                         )
! ----------------------------------------------------------------------

    implicit none

    ! input
    type(NoahmpIO_type), intent(in)    :: NoahmpIO
    INTEGER,                            INTENT(IN)   :: ids,ide, jds,jde, kds,kde,  &
                                                        ims,ime, jms,jme, kms,kme,  &
                                                        its,ite, jts,jte, kts,kte
    REAL,                                INTENT(IN)  :: DELTAT                                 
    INTEGER, DIMENSION(ims:ime,jms:jme), INTENT(IN)  :: ISLTYP, LANDMASK
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(IN)  :: FDEPTH,WTD,TOPO,AREA
    ! output
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(OUT) :: QLAT
    ! local
    INTEGER                                          :: I, J, itsh,iteh,jtsh,jteh
    REAL                                             :: Q,KLAT
    REAL,    DIMENSION(ims:ime,jms:jme)              :: KCELL, HEAD
    REAL,    DIMENSION(19)                           :: KLATFACTOR
    REAL,    PARAMETER                               :: PI = 3.14159265
    REAL,    PARAMETER                               :: FANGLE = 0.22754493   ! = 0.5*sqrt(0.5*tan(pi/8))
    DATA     KLATFACTOR /2.,3.,4.,10.,10.,12.,14.,20.,24.,28.,40.,48.,2.,0.,10.,0.,20.,2.,2./
! ----------------------------------------------------------------------

    itsh = max(its-1, ids)
    iteh = min(ite+1, ide-1)
    jtsh = max(jts-1, jds)
    jteh = min(jte+1, jde-1)

    do J = jtsh, jteh
       do I = itsh, iteh
          if ( FDEPTH(I,J) > 0.0 ) then
             KLAT = NoahmpIO%DKSAT_TABLE(ISLTYP(I,J)) * KLATFACTOR(ISLTYP(I,J))
             if ( WTD(I,J) < -1.5 ) then
                KCELL(I,J) = FDEPTH(I,J) * KLAT * exp((WTD(I,J) + 1.5) / FDEPTH(I,J))
             else
                KCELL(I,J) = KLAT * (WTD(I,J) + 1.5 + FDEPTH(I,J))  
             endif
           else
             KCELL(i,J) = 0.0
           endif
           HEAD(I,J) = TOPO(I,J) + WTD(I,J)
       enddo
    enddo

    itsh = max(its, ids+1)
    iteh = min(ite, ide-2)
    jtsh = max(jts, jds+1)
    jteh = min(jte, jde-2)

    do J = jtsh, jteh
       do I = itsh, iteh
          if ( LANDMASK(I,J) > 0 ) then
             Q = 0.0       
             Q = Q + (KCELL(I-1,J+1)+KCELL(I,J)) &
                   * (HEAD(I-1,J+1)-HEAD(I,J))/SQRT(2.0)                  
             Q = Q + (KCELL(I-1,J)+KCELL(I,J)) &
                   * (HEAD(I-1,J)-HEAD(I,J))
             Q = Q + (KCELL(I-1,J-1)+KCELL(I,J)) &
                   * (HEAD(I-1,J-1)-HEAD(I,J))/SQRT(2.0)
             Q = Q + (KCELL(I,J+1)+KCELL(I,J)) &
                   * (HEAD(I,J+1)-HEAD(I,J))
             Q = Q + (KCELL(I,J-1)+KCELL(I,J)) &
                   * (HEAD(I,J-1)-HEAD(I,J))
             Q = Q + (KCELL(I+1,J+1)+KCELL(I,J)) &
                   * (HEAD(I+1,J+1)-HEAD(I,J))/SQRT(2.0)
             Q = Q + (KCELL(I+1,J)+KCELL(I,J)) &
                   * (HEAD(I+1,J)-HEAD(I,J))
             Q = Q + (KCELL(I+1,J-1)+KCELL(I,J)) &
                   * (HEAD(I+1,J-1)-HEAD(I,J))/SQRT(2.0)
             QLAT(I,J) = FANGLE* Q * DELTAT / AREA(I,J)
          endif
       enddo
    enddo

  end subroutine LATERALFLOW

! ==================================================================================================
  subroutine UPDATEWTD (NSOIL,  DZS,  ZSOIL, SMCEQ,               & !in
                        SMCMAX, SMCWLT, PSISAT, BEXP, ILOC, JLOC, & !in
                        TOTWATER, WTD, SMC, SH2O, SMCWTD,         & !inout
                        QSPRING                                   )  !out
! ----------------------------------------------------------------------

    implicit none

    ! input
    INTEGER,                  INTENT(IN)    :: NSOIL        ! no. of soil layers
    INTEGER,                  INTENT(IN)    :: ILOC, JLOC
    REAL,                     INTENT(IN)    :: SMCMAX
    REAL,                     INTENT(IN)    :: SMCWLT
    REAL,                     INTENT(IN)    :: PSISAT
    REAL,                     INTENT(IN)    :: BEXP
    REAL, DIMENSION(0:NSOIL), INTENT(IN)    :: ZSOIL        ! depth of soil layer-bottom [m]
    REAL, DIMENSION(1:NSOIL), INTENT(IN)    :: SMCEQ        ! equilibrium soil water  content [m3/m3]
    REAL, DIMENSION(1:NSOIL), INTENT(IN)    :: DZS          ! soil layer thickness [m]
    ! input-output
    REAL,                     INTENT(INOUT) :: TOTWATER
    REAL,                     INTENT(INOUT) :: WTD
    REAL,                     INTENT(INOUT) :: SMCWTD
    REAL, DIMENSION(1:NSOIL), INTENT(INOUT) :: SMC
    REAL, DIMENSION(1:NSOIL), INTENT(INOUT) :: SH2O
    ! output
    REAL,                     INTENT(OUT)   :: QSPRING
    ! local
    INTEGER                                 :: K
    INTEGER                                 :: K1
    INTEGER                                 :: IWTD
    INTEGER                                 :: KWTD
    REAL                                    :: MAXWATUP, MAXWATDW ,WTDOLD
    REAL                                    :: WGPMID
    REAL                                    :: SYIELDDW
    REAL                                    :: DZUP
    REAL                                    :: SMCEQDEEP
    REAL, DIMENSION(1:NSOIL)                :: SICE
! -------------------------------------------------------------

    QSPRING = 0.0
    SICE = SMC - SH2O
    iwtd = 1

    ! case 1: totwater > 0 (water table going up):
    if ( totwater > 0.0 ) then
       if ( wtd >= zsoil(nsoil) ) then
          do k = nsoil-1, 1, -1
             if (wtd < zsoil(k)) exit
          enddo
          iwtd = k
          kwtd = iwtd + 1

          ! max water that fits in the layer
          maxwatup = dzs(kwtd) * (smcmax-smc(kwtd))
          if (totwater <= maxwatup) then
             smc(kwtd) = smc(kwtd) + totwater / dzs(kwtd)
             smc(kwtd) = min(smc(kwtd),smcmax)
             if (smc(kwtd) > smceq(kwtd)) wtd = min( (smc(kwtd)*dzs(kwtd) &
                       - smceq(kwtd)*zsoil(iwtd) + smcmax*zsoil(kwtd) ) / &
                       (smcmax-smceq(kwtd)), zsoil(iwtd) )
             totwater = 0.0
          else   ! water enough to saturate the layer
             smc(kwtd) = smcmax
             totwater = totwater - maxwatup
             k1 = iwtd
             do k = k1, 0, -1
                wtd = zsoil(k)
                iwtd = k-1
                if (k == 0) exit
                maxwatup = dzs(k) * (smcmax-smc(k))
                if (totwater <= maxwatup) then
                   smc(k) = smc(k) + totwater / dzs(k)
                   smc(k) = min(smc(k),smcmax)
                   if (smc(k) > smceq(k)) wtd = min( (smc(k)*dzs(k) &
                       - smceq(k)*zsoil(iwtd) + smcmax*zsoil(k) ) / &
                       (smcmax-smceq(k)), zsoil(iwtd) )
                   totwater = 0.0
                   exit
                else
                   smc(k) = smcmax
                   totwater = totwater - maxwatup
                endif
             enddo
          endif

       elseif ( wtd >= zsoil(nsoil)-dzs(nsoil) ) then ! wtd below bottom of soil model
          ! gmm equilibrium soil moisture content
          smceqdeep = smcmax * ( psisat / &
                      (psisat - dzs(nsoil)) ) ** (1.0/bexp)
          ! smceqdeep = max(smceqdeep, smcwlt)
          smceqdeep = max(smceqdeep, 1.0e-4)
          maxwatup = (smcmax-smcwtd) * dzs(nsoil)
          if ( totwater <= maxwatup ) then
             smcwtd = smcwtd + totwater / dzs(nsoil)
             smcwtd = min(smcwtd,smcmax)
             if (smcwtd > smceqdeep) wtd = min( (smcwtd*dzs(nsoil) &
                     - smceqdeep*zsoil(nsoil) + smcmax*(zsoil(nsoil)-dzs(nsoil))) / &
                     (smcmax-smceqdeep), zsoil(nsoil) )
             totwater = 0.0
          else
             smcwtd = smcmax
             totwater = totwater - maxwatup
             do k = nsoil, 0, -1
                wtd = zsoil(k)
                iwtd = k-1
                if (k == 0) exit
                maxwatup = dzs(k) * (smcmax-smc(k))
                if (totwater <= maxwatup) then
                   smc(k) = min(smc(k) + totwater / dzs(k),smcmax)
                   if (smc(k) > smceq(k)) wtd = min ( (smc(k)*dzs(k) &
                        - smceq(k)*zsoil(iwtd) + smcmax*zsoil(k) ) / &
                        (smcmax-smceq(k)), zsoil(iwtd) )
                   totwater = 0.0
                   exit
                else
                   smc(k) = smcmax
                   totwater = totwater - maxwatup
                endif
             enddo
          endif

       ! deep water table
       else
          maxwatup = (smcmax - smcwtd) * (zsoil(nsoil)-dzs(nsoil)-wtd)
          if (totwater <= maxwatup) then
             wtd = wtd + totwater / (smcmax - smcwtd)
             totwater = 0.0
          else
             totwater = totwater - maxwatup
             wtd = zsoil(nsoil) - dzs(nsoil)
             maxwatup = (smcmax - smcwtd) * dzs(nsoil)
             if (totwater <= maxwatup) then
                ! gmm equilibrium soil moisture content
                smceqdeep = smcmax * (psisat / &
                            (psisat - dzs(nsoil))) ** (1.0/bexp)
                ! smceqdeep = max(smceqdeep, smcwlt)
                smceqdeep = max(smceqdeep, 1.0e-4)
                smcwtd = smcwtd + totwater / dzs(nsoil)
                smcwtd = min(smcwtd, smcmax)
                wtd = ( smcwtd*dzs(nsoil) &
                      - smceqdeep*zsoil(nsoil) + smcmax*(zsoil(nsoil)-dzs(nsoil)) ) / &
                     (smcmax - smceqdeep)
                totwater = 0.0
             else
                smcwtd = smcmax
                totwater = totwater - maxwatup
                do k = nsoil, 0, -1
                   wtd = zsoil(k)
                   iwtd = k - 1
                   if (k == 0) exit
                   maxwatup = dzs(k) * (smcmax-smc(k))
                   if (totwater <= maxwatup) then
                      smc(k) = smc(k) + totwater / dzs(k)
                      smc(k) = min(smc(k),smcmax)
                      if (smc(k) > smceq(k)) wtd = (smc(k)*dzs(k) &
                           - smceq(k)*zsoil(iwtd) + smcmax*zsoil(k) ) / &
                           (smcmax-smceq(k))
                      totwater = 0.0
                      exit
                   else
                      smc(k) = smcmax
                      totwater = totwater - maxwatup
                   endif
                enddo
             endif
          endif
       endif

       ! water springing at the surface
       qspring = totwater

    ! case 2: totwater < 0 (water table going down):
    elseif ( totwater < 0.0 ) then
       if (wtd >= zsoil(nsoil)) then ! wtd in the resolved layers
          do k = nsoil-1, 1, -1
             if (wtd < zsoil(k)) exit
          enddo
          iwtd = k
          k1 = iwtd + 1
          do kwtd = k1, nsoil
             ! max water that the layer can yield
             maxwatdw = dzs(kwtd) * (smc(kwtd) - max(smceq(kwtd),sice(kwtd)))
             if (-totwater <= maxwatdw) then
                smc(kwtd) = smc(kwtd) + totwater / dzs(kwtd)
                if (smc(kwtd) > smceq(kwtd)) then
                   wtd = ( smc(kwtd) * dzs(kwtd) &
                         - smceq(kwtd)*zsoil(iwtd) + smcmax*zsoil(kwtd) ) / &
                         (smcmax-smceq(kwtd))
                else
                   wtd = zsoil(kwtd)
                   iwtd = iwtd + 1
                endif
                totwater = 0.0
                exit
             else
                wtd = zsoil(kwtd)
                iwtd = iwtd + 1
                if (maxwatdw >= 0.0) then
                   smc(kwtd) = smc(kwtd) + maxwatdw / dzs(kwtd)
                   totwater = totwater + maxwatdw
                endif
             endif
          enddo

          if (iwtd == nsoil .and. totwater < 0.0) then
             ! gmm equilibrium soil moisture content
             smceqdeep = smcmax * (psisat / &
                         (psisat - dzs(nsoil))) ** (1.0/bexp)
             ! smceqdeep = max(smceqdeep, smcwlt)
             smceqdeep = max(smceqdeep, 1.0e-4)
             maxwatdw = dzs(nsoil) * (smcwtd - smceqdeep)
             if (-totwater <= maxwatdw) then
                smcwtd = smcwtd + totwater / dzs(nsoil)
                wtd = max( (smcwtd*dzs(nsoil) &
                      - smceqdeep*zsoil(nsoil) + smcmax*(zsoil(nsoil)-dzs(nsoil)) ) / &
                      (smcmax-smceqdeep), zsoil(nsoil)-dzs(nsoil) )
             else
                wtd = zsoil(nsoil) - dzs(nsoil)
                smcwtd = smcwtd + totwater / dzs(nsoil)
                ! and now even further down
                dzup = (smceqdeep-smcwtd) * dzs(nsoil) / (smcmax-smceqdeep)
                wtd = wtd - dzup
                smcwtd = smceqdeep
             endif
          endif

       elseif ( wtd >= zsoil(nsoil)-dzs(nsoil) ) then
          ! if wtd was already below the bottom of the resolved soil crust
          ! gmm equilibrium soil moisture content
          smceqdeep = smcmax * ( psisat / &
                      (psisat - dzs(nsoil)) ) ** (1./bexp)
          ! smceqdeep = max(smceqdeep,smcwlt)
          smceqdeep = max(smceqdeep, 1.0e-4)
          maxwatdw = dzs(nsoil) * (smcwtd - smceqdeep)
          if (-totwater <= maxwatdw) then
             smcwtd = smcwtd + totwater / dzs(nsoil)
             wtd = max( (smcwtd*dzs(nsoil) &
                   - smceqdeep*zsoil(nsoil) + smcmax*(zsoil(nsoil)-dzs(nsoil)) ) / &
                   (smcmax-smceqdeep), zsoil(nsoil)-dzs(nsoil) )
          else
             wtd = zsoil(nsoil) - dzs(nsoil)
             smcwtd = smcwtd + totwater / dzs(nsoil)
             ! and now even further down
             dzup = (smceqdeep-smcwtd) * dzs(nsoil) / (smcmax-smceqdeep)
             wtd = wtd - dzup
             smcwtd = smceqdeep
          endif
       else
          ! gmm equilibrium soil moisture content
          wgpmid = smcmax * (psisat / &
                   (psisat - (zsoil(nsoil)-wtd))) ** (1.0/bexp)
          ! wgpmid = max(wgpmid, smcwlt)
          wgpmid = max(wgpmid, 1.0e-4)
          syielddw = smcmax - wgpmid
          wtdold = wtd
          wtd = wtdold + totwater / syielddw
          ! update wtdwgp
          smcwtd = (smcwtd * (zsoil(nsoil)-wtdold) + wgpmid * (wtdold-wtd)) / (zsoil(nsoil)-wtd)
       endif

       qspring=0.0

    endif
         
    SH2O = SMC - SICE

  end subroutine UPDATEWTD
! ----------------------------------------------------------------------

end module GroundWaterMmfMod
