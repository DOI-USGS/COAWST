#include "cppdefs.h"
#ifdef STATIONS
      SUBROUTINE read_StaPar (model, inp, out, Lwrite)
!
!svn $Id: read_stapar.F 1054 2021-03-06 19:47:12Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2021 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine reads and reports stations input parameters.           !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_ncparam
# if defined SEDIMENT || defined BBL_MODEL
      USE mod_sediment
# endif
      USE mod_scalars
!
      USE inp_decode_mod
!
      implicit none
!
!  Imported variable declarations
!
      logical, intent(in) :: Lwrite
      integer, intent(in) :: model, inp, out
!
!  Local variable declarations.
!
      integer :: Mstation, Npts, Nval
      integer :: flag, i, igrid, ista, itrc, ng, status

      real(r8) :: Xpos, Ypos

# if defined SOLVE3D && (defined BBL_MODEL || defined SEDIMENT)
      logical, dimension(MBOTP,Ngrids) :: Lbottom
# endif
      logical, dimension(Ngrids) :: Lswitch
      logical, dimension(MT,Ngrids) :: Ltracer

      integer, dimension(Ngrids) :: is

      real(dp), dimension(nRval) :: Rval

      character (len=40 ) :: KeyWord
      character (len=256) :: line
      character (len=256), dimension(nCval) :: Cval
!
!-----------------------------------------------------------------------
!  Read in stations parameters.
!-----------------------------------------------------------------------
!
      DO WHILE (.TRUE.)
        READ (inp,'(a)',ERR=20,END=30) line
        status=decode_line(line, KeyWord, Nval, Cval, Rval)
        IF (status.gt.0) THEN
          SELECT CASE (TRIM(KeyWord))
            CASE ('Lstations')
              Npts=load_l(Nval, Cval, Ngrids, Lstations)
            CASE ('Sout(idUvel)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idUvel,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Sout(idVvel)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idVvel,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Sout(idu3dE)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idu3dE,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Sout(idv3dN)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idv3dN,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Sout(idWvel)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idWvel,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Sout(idOvel)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idOvel,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Sout(idUbar)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idUbar,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Sout(idVbar)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idVbar,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Sout(idu2dE)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idu2dE,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Sout(idv2dN)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idv2dN,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Sout(idFsur)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idFsur,1:Ngrids)=Lswitch(1:Ngrids)
# if defined SEDIMENT && defined SED_MORPH
            CASE ('Sout(idBath)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idBath,1:Ngrids)=Lswitch(1:Ngrids)
# endif
# ifdef PRIMARY_PROD
            CASE ('Sout(idNPP)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idNPP,1:Ngrids)=Lswitch(1:Ngrids)
# endif
            CASE ('Sout(idTvar)')
              Npts=load_l(Nval, Cval, MT, Ngrids, Ltracer)
              DO ng=1,Ngrids
                DO itrc=1,NT(ng)
                  Sout(idTvar(itrc),ng)=Ltracer(itrc,ng)
                END DO
              END DO
            CASE ('Sout(idUsms)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idUsms,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Sout(idVsms)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idVsms,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Sout(idUbms)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idUbms,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Sout(idVbms)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idVbms,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Sout(idUbrs)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idUbrs,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Sout(idVbrs)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idVbrs,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Sout(idUbws)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idUbws,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Sout(idVbws)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idVbws,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Sout(idUbcs)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idUbcs,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Sout(idVbcs)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idVbcs,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Sout(idUbot)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idUbot,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Sout(idVbot)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idVbot,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Sout(idUbur)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idUbur,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Sout(idVbvr)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idVbvr,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Sout(idW2xx)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idW2xx,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Sout(idW2xy)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idW2xy,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Sout(idW2yy)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idW2yy,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Sout(idU2rs)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idU2rs,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Sout(idV2rs)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idV2rs,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Sout(idU2Sd)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idU2Sd,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Sout(idV2Sd)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idV2Sd,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Sout(idW3xx)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idW3xx,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Sout(idW3xy)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idW3xy,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Sout(idW3yy)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idW3yy,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Sout(idW3zx)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idW3zx,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Sout(idW3zy)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idW3zy,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Sout(idU3rs)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idU3rs,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Sout(idV3rs)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idV3rs,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Sout(idU3Sd)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idU3Sd,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Sout(idV3Sd)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idV3Sd,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Sout(idWamp)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idWamp,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Sout(idWlen)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idWlen,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Sout(idWdir)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idWdir,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Sout(idWdip)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idWdip,1:Ngrids)=Lswitch(1:Ngrids)
# ifdef WAVES_TOP_PERIOD
            CASE ('Sout(idWptp)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idWptp,1:Ngrids)=Lswitch(1:Ngrids)
# endif
# ifdef WAVES_BOT_PERIOD
            CASE ('Sout(idWpbt)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idWpbt,1:Ngrids)=Lswitch(1:Ngrids)
# endif
# ifdef WAVES_UB
            CASE ('Sout(idWorb)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idWorb,1:Ngrids)=Lswitch(1:Ngrids)
# endif
# if defined TKE_WAVEDISS || defined WAVES_OCEAN
            CASE ('Sout(idWdis)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idWdis,1:Ngrids)=Lswitch(1:Ngrids)
# endif
# ifdef SOLVE3D
#  if defined BULK_FLUXES || defined ECOSIM || defined ATM_PRESS
            CASE ('Sout(idPair)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idPair,1:Ngrids)=Lswitch(1:Ngrids)
#  endif
#  if defined BULK_FLUXES || defined ECOSIM
            CASE ('Sout(idUair)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idUair,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Sout(idVair)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idVair,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Sout(idUairE)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idUairE,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Sout(idVairN)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idVairN,1:Ngrids)=Lswitch(1:Ngrids)
#  endif
            CASE ('Sout(idTsur)')
              Npts=load_l(Nval, Cval, MT, Ngrids, Ltracer)
              DO ng=1,Ngrids
                DO itrc=1,NAT
                  Sout(idTsur(itrc),ng)=Ltracer(itrc,ng)
                END DO
              END DO
            CASE ('Sout(idLhea)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idLhea,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Sout(idShea)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idShea,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Sout(idLrad)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idLrad,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Sout(idSrad)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idSrad,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Sout(idEmPf)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idEmPf,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Sout(idevap)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idevap,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Sout(idrain)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idrain,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Sout(idDano)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idDano,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Sout(idVvis)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idVvis,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Sout(idTdif)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idTdif,1:Ngrids)=Lswitch(1:Ngrids)
#  ifdef SALINITY
            CASE ('Sout(idSdif)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idSdif,1:Ngrids)=Lswitch(1:Ngrids)
#  endif
            CASE ('Sout(idHsbl)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idHsbl,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Sout(idHbbl)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idHbbl,1:Ngrids)=Lswitch(1:Ngrids)
#  ifdef NEMURO_SED1
            CASE ('Sout(idPONsed)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idPONsed,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Sout(idOPALsed)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idOPALsed,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Sout(idDENITsed)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idDENITsed,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Sout(idPONbur)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idPONbur,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Sout(idOPALbur)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idOPALbur,1:Ngrids)=Lswitch(1:Ngrids)
#  endif
            CASE ('Sout(idMtke)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idMtke,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Sout(idMtls)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idMtls,1:Ngrids)=Lswitch(1:Ngrids)
#  ifdef WET_DRY
            CASE ('Sout(idRwet)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idRwet,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Sout(idUwet)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idUwet,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Sout(idVwet)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idVwet,1:Ngrids)=Lswitch(1:Ngrids)
#  endif
#  if defined BBL_MODEL || defined SEDIMENT
            CASE ('Sout(isd50)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              i=idBott(isd50)
              DO ng=1,Ngrids
                Sout(i,ng)=Lswitch(ng)
              END DO
            CASE ('Sout(idens)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              i=idBott(idens)
              DO ng=1,Ngrids
                Sout(i,ng)=Lswitch(ng)
             END DO
            CASE ('Sout(iwsed)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              i=idBott(iwsed)
              DO ng=1,Ngrids
                Sout(i,ng)=Lswitch(ng)
              END DO
            CASE ('Sout(itauc)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              i=idBott(itauc)
              DO ng=1,Ngrids
                Sout(i,ng)=Lswitch(ng)
              END DO
            CASE ('Sout(irlen)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              i=idBott(irlen)
              DO ng=1,Ngrids
                Sout(i,ng)=Lswitch(ng)
              END DO
            CASE ('Sout(irhgt)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              i=idBott(irhgt)
              DO ng=1,Ngrids
                Sout(i,ng)=Lswitch(ng)
              END DO
            CASE ('Sout(ibwav)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              i=idBott(ibwav)
              DO ng=1,Ngrids
                Sout(i,ng)=Lswitch(ng)
              END DO
            CASE ('Sout(izdef)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              i=idBott(izdef)
              DO ng=1,Ngrids
                Sout(i,ng)=Lswitch(ng)
              END DO
            CASE ('Sout(izapp)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              i=idBott(izapp)
              DO ng=1,Ngrids
                Sout(i,ng)=Lswitch(ng)
              END DO
            CASE ('Sout(izNik)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              i=idBott(izNik)
             DO ng=1,Ngrids
                Sout(i,ng)=Lswitch(ng)
              END DO
            CASE ('Sout(izbio)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              i=idBott(izbio)
              DO ng=1,Ngrids
                Sout(i,ng)=Lswitch(ng)
              END DO
            CASE ('Sout(izbfm)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              i=idBott(izbfm)
              DO ng=1,Ngrids
               Sout(i,ng)=Lswitch(ng)
              END DO
            CASE ('Sout(izbld)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              i=idBott(izbld)
              DO ng=1,Ngrids
                Sout(i,ng)=Lswitch(ng)
              END DO
            CASE ('Sout(izwbl)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              i=idBott(izwbl)
              DO ng=1,Ngrids
                Sout(i,ng)=Lswitch(ng)
              END DO
            CASE ('Sout(iactv)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              i=idBott(iactv)
              DO ng=1,Ngrids
                Sout(i,ng)=Lswitch(ng)
              END DO
            CASE ('Sout(ishgt)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              i=idBott(ishgt)
              DO ng=1,Ngrids
                Sout(i,ng)=Lswitch(ng)
              END DO
#  endif
# endif
            CASE ('NSTATION')
              Npts=load_i(Nval, Rval, Ngrids, Nstation)
              DO ng=1,Ngrids
                IF (.not.Lstations(ng)) THEN
                  Nstation(ng)=0
                ELSE
                  IF (Nstation(ng).le.0) THEN
                    IF (Master) WRITE (out,40) ng, 'Nstation',          &
     &                                         Nstation(ng),            &
     &                 'Must be positive and greater than zero.'
                    exit_flag=4
                    RETURN
                  END IF
                END IF
              END DO
            CASE ('POS')
              DO ng=1,Ngrids
                IF (Lstations(ng)) THEN
                  allocate ( SCALARS(ng) % Sflag(Nstation(ng)) )
                  allocate ( SCALARS(ng) % SposX(Nstation(ng)) )
                  allocate ( SCALARS(ng) % SposY(Nstation(ng)) )
                  Dmem(ng)=Dmem(ng)+3.0_r8*REAL(Nstation(ng),r8)
                END IF
              END DO
              is(1:Ngrids)=0
              DO WHILE (.TRUE.)
                READ (inp,*,ERR=10,END=10) igrid, flag, Xpos, Ypos
                ng=MAX(1,MIN(ABS(igrid),Ngrids))
                IF (Lstations(ng)) THEN
                  is(ng)=is(ng)+1
                  SCALARS(ng)%Sflag(is(ng))=flag
                  SCALARS(ng)%SposX(is(ng))=Xpos
                  SCALARS(ng)%SposY(is(ng))=Ypos
                END IF
              END DO
 10           DO ng=1,Ngrids
                IF (Lstations(ng).and.(Nstation(ng).ne.is(ng))) THEN
                  IF (Master) WRITE (out,50) Nstation(ng), is(ng)
                  exit_flag=4
                  RETURN
                END IF
              END DO
          END SELECT
        END IF
      END DO
 20   IF (Master) WRITE (out,60) line
      exit_flag=4
      RETURN
 30   CONTINUE
!
!-----------------------------------------------------------------------
!  Process input parameters.
!-----------------------------------------------------------------------
!
!  Turn off the processing of stations if not running long enough to
!  create a stations file (LdefSTA=.FALSE. because nSTA < ntimes or
!  nSTA = 0 when nrrec = 0).
!
      DO ng=1,Ngrids
        IF (.not.LdefSTA(ng).and.Lstations(ng)) THEN
          Lstations(ng)=.FALSE.
        END IF
      END DO
!
!  Make sure that both component switches are activated when processing
!  (Eastward,Northward) momentum components at RHO-points.
!
      DO ng=1,Ngrids
        IF (.not.Sout(idu2dE,ng).and.Sout(idv2dN,ng)) THEN
          Sout(idu2dE,ng)=.TRUE.
        END IF
        IF (.not.Sout(idv2dN,ng).and.Sout(idu2dE,ng)) THEN
          Sout(idv2dN,ng)=.TRUE.
        END IF
# ifdef SOLVE3D
        IF (.not.Sout(idu3dE,ng).and.Sout(idv3dN,ng)) THEN
          Sout(idu3dE,ng)=.TRUE.
        END IF
        IF (.not.Sout(idv3dN,ng).and.Sout(idu3dE,ng)) THEN
          Sout(idv3dN,ng)=.TRUE.
        END IF
# endif
      END DO
!
!-----------------------------------------------------------------------
!  Report input parameters.
!-----------------------------------------------------------------------
!
      IF (Master.and.Lwrite) THEN
        DO ng=1,Ngrids
          IF (Lstations(ng)) THEN
            WRITE (out,70) ng
            WRITE (out,80) Nstation(ng), 'Nstation',                    &
     &            'Number of stations to write out into stations file.'
#if defined SEDIMENT && defined SED_MORPH
            IF (Sout(idbath,ng)) WRITE (out,90) Sout(idbath,ng),        &
     &          'Sout(idbath)',                                         &
     &          'Write out free-surface.'
#endif
# ifdef PRIMARY_PROD
            IF (Sout(idNPP,ng)) WRITE (out,90) Sout(idNPP,ng),          &
     &          'Sout(idNPP)',                                          &
     &          'Write out net primary productivity.'
# endif
            IF (Sout(idFsur,ng)) WRITE (out,90) Sout(idFsur,ng),        &
     &          'Sout(idFsur)',                                         &
     &          'Write out free-surface.'
            IF (Sout(idUbar,ng)) WRITE (out,90) Sout(idUbar,ng),        &
     &          'Sout(idUbar)',                                         &
     &          'Write out 2D U-momentum component.'
            IF (Sout(idVbar,ng)) WRITE (out,90) Sout(idVbar,ng),        &
     &          'Sout(idVbar)',                                         &
     &          'Write out 2D V-momentum component.'
            IF (Sout(idu2dE,ng)) WRITE (out,90) Sout(idu2dE,ng),        &
     &          'Sout(idu2dE)',                                         &
     &          'Write out 2D U-eastward  at RHO-points.'
            IF (Sout(idv2dN,ng)) WRITE (out,90) Sout(idv2dN,ng),        &
     &          'Sout(idv2dN)',                                         &
     &          'Write out 2D V-northward at RHO-points.'
# ifdef SOLVE3D
            IF (Sout(idUvel,ng)) WRITE (out,90) Sout(idUvel,ng),        &
     &          'Sout(idUvel)',                                         &
     &          'Write out 3D U-momentum component.'
            IF (Sout(idVvel,ng)) WRITE (out,90) Sout(idVvel,ng),        &
     &          'Sout(idVvel)',                                         &
     &          'Write out 3D V-momentum component.'
            IF (Sout(idu3dE,ng)) WRITE (out,90) Sout(idu3dE,ng),        &
     &          'Sout(idu3dE)',                                         &
     &          'Write out 3D U-eastward  at RHO-points.'
            IF (Sout(idv3dN,ng)) WRITE (out,90) Sout(idv3dN,ng),        &
     &          'Sout(idv3dN)',                                         &
     &          'Write out 3D V-northward at RHO-points.'
            IF (Sout(idWvel,ng)) WRITE (out,90) Sout(idWvel,ng),        &
     &          'Sout(idWvel)',                                         &
     &          'Write out W-momentum component.'
            IF (Sout(idOvel,ng)) WRITE (out,90) Sout(idOvel,ng),        &
     &          'Sout(idOvel)',                                         &
     &          'Write out omega vertical velocity.'
            DO itrc=1,NT(ng)
              IF (Sout(idTvar(itrc),ng)) WRITE (out,100)                &
     &            Sout(idTvar(itrc),ng), 'Sout(idTvar)',                &
     &            'Write out tracer ', itrc, TRIM(Vname(1,idTvar(itrc)))
            END DO
# endif
            IF (Sout(idUsms,ng)) WRITE (out,90) Sout(idUsms,ng),        &
     &          'Sout(idUsms)',                                         &
     &          'Write out surface U-momentum stress.'
            IF (Sout(idVsms,ng)) WRITE (out,90) Sout(idVsms,ng),        &
     &          'Sout(idVsms)',                                         &
     &          'Write out surface V-momentum stress.'
            IF (Sout(idUbms,ng)) WRITE (out,90) Sout(idUbms,ng),        &
     &          'Sout(idUbms)',                                         &
     &          'Write out bottom U-momentum stress.'
            IF (Sout(idVbms,ng)) WRITE (out,90) Sout(idVbms,ng),        &
     &          'Sout(idVbms)',                                         &
     &          'Write out bottom V-momentum stress.'
# ifdef BBL_MODEL
            IF (Sout(idUbrs,ng)) WRITE (out,90) Sout(idUbrs,ng),        &
     &          'Sout(idUbrs)',                                         &
     &          'Write out bottom U-current stress.'
            IF (Sout(idVbrs,ng)) WRITE (out,90) Sout(idVbrs,ng),        &
     &          'Sout(idVbrs)',                                         &
     &          'Write out bottom V-current stress.'
            IF (Sout(idUbws,ng)) WRITE (out,90) Sout(idUbws,ng),        &
     &          'Sout(idUbws)',                                         &
     &          'Write out wind-induced, bottom U-wave stress.'
            IF (Sout(idVbws,ng)) WRITE (out,90) Sout(idVbws,ng),        &
     &          'Sout(idVbws)',                                         &
     &          'Write out wind-induced, bottom V-wave stress.'
            IF (Sout(idUbcs,ng)) WRITE (out,90) Sout(idUbcs,ng),        &
     &          'Sout(idUbcs)',                                         &
     &          'Write out max wind + current, bottom U-wave stress.'
            IF (Sout(idVbcs,ng)) WRITE (out,90) Sout(idVbcs,ng),        &
     &          'Sout(idVbcs)',                                         &
     &          'Write out max wind + current, bottom V-wave stress.'
            IF (Sout(idUbot,ng)) WRITE (out,90) Sout(idUbot,ng),        &
     &          'Sout(idUbot)',                                         &
     &          'Write out bed wave orbital U-velocity.'
            IF (Sout(idVbot,ng)) WRITE (out,90) Sout(idVbot,ng),        &
     &          'Sout(idVbot)',                                         &
     &          'Write out bed wave orbital V-velocity.'
            IF (Sout(idUbur,ng)) WRITE (out,90) Sout(idUbur,ng),        &
     &          'Sout(idUbur)',                                         &
     &          'Write out bottom U-velocity above bed.'
            IF (Sout(idVbvr,ng)) WRITE (out,90) Sout(idVbvr,ng),        &
     &          'Sout(idVbvr)',                                         &
     &          'Write out bottom V-velocity above bed.'
# endif
# if defined NEARSHORE_MELLOR
            IF (Sout(idW2xx,ng)) WRITE (out,90) Sout(idW2xx,ng),        &
     &         'Sout(idW2xx)',                                          &
     &         'Write out 2D radiation stress, Sxx.'
            IF (Sout(idW2xy,ng)) WRITE (out,90) Sout(idW2xy,ng),        &
     &         'Sout(idW2xy)',                                          &
     &         'Write out 2D radiation stress, Sxy.'
            IF (Sout(idW2yy,ng)) WRITE (out,90) Sout(idW2yy,ng),        &
     &         'Sout(idW2yy)',                                          &
     &         'Write out 2D radiation stress, Syy.'
            IF (Sout(idU2rs,ng)) WRITE (out,90) Sout(idU2rs,ng),        &
     &         'Sout(idU2rs)',                                          &
     &         'Write out total 2D u-radiation stress.'
            IF (Sout(idV2rs,ng)) WRITE (out,90) Sout(idV2rs,ng),        &
     &         'Sout(idV2rs)',                                          &
     &         'Write out total 2D v-radiation stress.'
            IF (Sout(idU2Sd,ng)) WRITE (out,90) Sout(idU2Sd,ng),        &
     &         'Sout(idU2Sd)',                                          &
     &         'Write out 2D u-momentum stokes velocity.'
            IF (Sout(idV2Sd,ng)) WRITE (out,90) Sout(idV2Sd,ng),        &
     &         'Sout(idV2Sd)',                                          &
     &         'Write out 2D v-momentum stokes velocity.'
#  ifdef SOLVE3D
            IF (Sout(idW3xx,ng)) WRITE (out,90) Sout(idW3xx,ng),        &
     &         'Sout(idW3xx)',                                          &
     &         'Write out 3D horizonrtal radiation stress, Sxx.'
            IF (Sout(idW3xy,ng)) WRITE (out,90) Sout(idW3xy,ng),        &
     &         'Sout(idW3xy)',                                          &
     &         'Write out 3D horizonrtal radiation stress, Sxy.'
            IF (Sout(idW3yy,ng)) WRITE (out,90) Sout(idW3yy,ng),        &
     &         'Sout(idW3yy)',                                          &
     &         'Write out 3D horizonrtal radiation stress, Syy.'
            IF (Sout(idW3zx,ng)) WRITE (out,90) Sout(idW3zx,ng),        &
     &         'Sout(idW3zx)',                                          &
     &         'Write out 3D vertical radiation stress, Spx.'
            IF (Sout(idW3zy,ng)) WRITE (out,90) Sout(idW3zy,ng),        &
     &         'Sout(idW3zy)',                                          &
     &         'Write out 3D vertical radiation stress, Spy.'
            IF (Sout(idU3rs,ng)) WRITE (out,90) Sout(idU3rs,ng),        &
     &         'Sout(idU3rs)',                                          &
     &         'Write out total 3D u-radiation stress.'
            IF (Sout(idV3rs,ng)) WRITE (out,90) Sout(idV3rs,ng),        &
     &         'Sout(idV3rs)',                                          &
     &         'Write out total 3D v-radiation stress.'
            IF (Sout(idU3Sd,ng)) WRITE (out,90) Sout(idU3Sd,ng),        &
     &         'Sout(idU3Sd)',                                          &
     &         'Write out 3D u-momentum stokes velocity.'
            IF (Sout(idV3Sd,ng)) WRITE (out,90) Sout(idV3Sd,ng),        &
     &         'Sout(idV3Sd)',                                          &
     &         'Write out 3D v-momentum stokes velocity.'
#  endif
# endif
# ifdef WAVES_HEIGHT
            IF (Sout(idWamp,ng)) WRITE (out,90) Sout(idWamp,ng),        &
     &         'Sout(idWamp)',                                          &
     &         'Write out wave height.'
# endif
# ifdef WAVES_LENGTH
            IF (Sout(idWlen,ng)) WRITE (out,90) Sout(idWlen,ng),        &
     &         'Sout(idWlen)',                                          &
     &         'Write out wave length.'
# endif
# ifdef WAVES_DIR
            IF (Sout(idWdir,ng)) WRITE (out,90) Sout(idWdir,ng),        &
     &         'Sout(idWdir)',                                          &
     &         'Write out mean wave direction.'
# endif
# ifdef WAVES_DIRP
            IF (Sout(idWdip,ng)) WRITE (out,90) Sout(idWdip,ng),        &
     &         'Sout(idWdip)',                                          &
     &         'Write out peak wave direction.'
# endif
# ifdef WAVES_TOP_PERIOD
            IF (Sout(idWptp,ng)) WRITE (out,90) Sout(idWptp,ng),        &
     &         'Sout(idWptp)',                                          &
     &         'Write out wave surface period.'
# endif
# ifdef WAVES_BOT_PERIOD
            IF (Sout(idWpbt,ng)) WRITE (out,90) Sout(idWpbt,ng),        &
     &         'Sout(idWpbt)',                                          &
     &         'Write out wave bottom period.'
# endif
# ifdef WAVES_UB
            IF (Sout(idWorb,ng)) WRITE (out,90) Sout(idWorb,ng),        &
     &         'Sout(idWorb)',                                          &
     &         'Write out wave bottom orbital velocity.'
# endif
# if defined TKE_WAVEDISS || defined WAVES_OCEAN
            IF (Sout(idWdis,ng)) WRITE (out,90) Sout(idWdis,ng),        &
     &         'Sout(idWdis)',                                          &
     &         'Write out wave dissipation.'
# endif
# if defined SOLVE3D && (defined BBL_MODEL || defined SEDIMENT)
            DO itrc=1,MBOTP
              IF (Sout(idBott(itrc),ng)) WRITE (out,100)                &
     &            Sout(idBott(itrc),ng), 'Sout(idBott)',                &
     &            'Write out bottom property ', itrc,                   &
     &            TRIM(Vname(1,idBott(itrc)))
            END DO
# endif
# ifdef SOLVE3D
#  if defined BULK_FLUXES || defined ECOSIM || defined ATM_PRESS
            IF (Sout(idPair,ng)) WRITE (out,90) Sout(idPair,ng),        &
     &          'Sout(idPair)',                                         &
     &          'Write out surface air pressure.'
#  endif
#  if defined BULK_FLUXES || defined ECOSIM
            IF (Sout(idUair,ng)) WRITE (out,90) Sout(idUair,ng),        &
     &          'Sout(idUair)',                                         &
     &          'Write out surface U-wind component.'
            IF (Sout(idVair,ng)) WRITE (out,90) Sout(idVair,ng),        &
     &          'Sout(idVair)',                                         &
     &          'Write out surface V-wind component.'
            IF (Sout(idUairE,ng)) WRITE (out,90) Sout(idUairE,ng),       &
     &          'Sout(idUairE)',                                         &
     &          'Write out surface Eastward U-wind component.'
            IF (Sout(idVairN,ng)) WRITE (out,90) Sout(idVairN,ng),       &
     &          'Sout(idVairN)',                                         &
     &          'Write out surface Northward V-wind component.'
#  endif
            IF (Sout(idTsur(itemp),ng)) WRITE (out,90)                  &
     &          Sout(idTsur(itemp),ng), 'Sout(idTsur)',                 &
     &          'Write out surface net heat flux.'
#  ifdef SALINITY
            IF (Sout(idTsur(isalt),ng)) WRITE (out,90)                  &
     &          Sout(idTsur(isalt),ng), 'Sout(idTsur)',                 &
     &          'Write out surface net salt flux.'
#  endif
#  ifdef SHORTWAVE
            IF (Sout(idSrad,ng)) WRITE (out,90) Sout(idSrad,ng),        &
     &          'Sout(idSrad)',                                         &
     &          'Write out shortwave radiation flux.'
#  endif
#  ifdef BULK_FLUXES
            IF (Sout(idLrad,ng)) WRITE (out,90) Sout(idLrad,ng),        &
     &          'Sout(idLrad)',                                         &
     &          'Write out longwave radiation flux.'
            IF (Sout(idLhea,ng)) WRITE (out,90) Sout(idLhea,ng),        &
     &          'Sout(idLhea)',                                         &
     &          'Write out latent heat flux.'
            IF (Sout(idShea,ng)) WRITE (out,90) Sout(idShea,ng),        &
     &          'Sout(idShea)',                                         &
     &          'Write out sensible heat flux.'
#   ifdef EMINUSP
            IF (Sout(idEmPf,ng)) WRITE (out,90) Sout(idEmPf,ng),        &
     &         'Sout(idEmPf)',                                          &
     &         'Write out E-P flux.'
            IF (Sout(idevap,ng)) WRITE (out,90) Sout(idevap,ng),        &
     &         'Sout(idevap)',                                          &
     &         'Write out evaporation rate.'
            IF (Sout(idrain,ng)) WRITE (out,90) Sout(idrain,ng),        &
     &         'Sout(idrain)',                                          &
     &         'Write out rain rate.'
#   endif
#  endif
            IF (Sout(idDano,ng)) WRITE (out,90) Sout(idDano,ng),        &
     &          'Sout(idDano)',                                         &
     &          'Write out density anomaly.'
            IF (Sout(idVvis,ng)) WRITE (out,90) Sout(idVvis,ng),        &
     &          'Sout(idVvis)',                                         &
     &          'Write out vertical viscosity coefficient.'
            IF (Sout(idTdif,ng)) WRITE (out,90) Sout(idTdif,ng),        &
     &          'Sout(idTdif)',                                         &
     &          'Write out vertical T-diffusion coefficient.'
#  ifdef SALINITY
            IF (Sout(idSdif,ng)) WRITE (out,90) Sout(idSdif,ng),        &
     &          'Sout(idSdif)',                                         &
     &          'Write out vertical S-diffusion coefficient.'
#  endif
#  ifdef LMD_SKPP
            IF (Sout(idHsbl,ng)) WRITE (out,90) Sout(idHsbl,ng),        &
     &          'Sout(idHsbl)',                                         &
     &          'Write out depth of surface boundary layer.'
#  endif
#  ifdef LMD_BKPP
            IF (Sout(idHbbl,ng)) WRITE (out,90) Sout(idHbbl,ng),        &
     &          'Sout(idHbbl)',                                         &
     &          'Write out depth of bottom boundary layer.'
#  endif
#  if defined GLS_MIXING || defined MY25_MIXING
            IF (Sout(idMtke,ng)) WRITE (out,90) Sout(idMtke,ng),        &
     &          'Sout(idMtke)',                                         &
     &          'Write out turbulent kinetic energy.'
            IF (Sout(idMtls,ng)) WRITE (out,90) Sout(idMtls,ng),        &
     &          'Sout(idMtls)',                                         &
     &          'Write out turbulent generic length-scale.'
#  endif
#  ifdef WET_DRY
            IF (Sout(idRwet,ng)) WRITE (out,90) Sout(idRwet,ng),        &
     &          'Sout(idRwet)',                                         &
     &          'Write out rho point wet-dry mask.'
            IF (Sout(idUwet,ng)) WRITE (out,90) Sout(idUwet,ng),        &
     &          'Sout(idUwet)',                                         &
     &          'Write out u point wet-dry mask.'
            IF (Sout(idVwet,ng)) WRITE (out,90) Sout(idVwet,ng),        &
     &          'Sout(idVwet)',                                         &
     &          'Write out v point wet-dry mask.'
#  endif
# endif
            WRITE (out,*)
            DO i=1,Nstation(ng)
              WRITE (out,110) i, SCALARS(ng)%Sflag(i),                  &
     &                           SCALARS(ng)%SposX(i),                  &
     &                           SCALARS(ng)%SposY(i)
            END DO
          END IF
        END DO
      END IF

  40  FORMAT (/,' READ_StaPar - Grid = ',i2.2,',',3x,                   &
     &        'Illegal value for ',a,' = ', i8,/,15x,a)
  50  FORMAT (/,' READ_StaPar - Inconsistent number of stations, ',     &
     &        'Nstation = ',2i8,/,15x,                                  &
     &        'change stations input script values.')
  60  FORMAT (/,' READ_StaPar - Error while processing line: ',/,a)

  70  FORMAT (/,/,' Stations Parameters, Grid: ',i2.2,                  &
     &        /,  ' =============================',/)
  80  FORMAT (1x,i10,2x,a,t30,a)
  90  FORMAT (10x,l1,2x,a,t30,a)
 100  FORMAT (10x,l1,2x,a,t30,a,i2.2,':',1x,a)
 110  FORMAT (13x,'Flag and positions for station ',i4.4,':',           &
     &        i3,1x,2f10.4)

      RETURN
      END SUBROUTINE read_StaPar
#else
      SUBROUTINE read_StaPar
      END SUBROUTINE read_StaPar
#endif
