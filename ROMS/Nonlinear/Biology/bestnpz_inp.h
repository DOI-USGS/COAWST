      SUBROUTINE read_BioPar (model, inp, out, Lwrite)
!
!=======================================================================
!                                                                      !
!  This routine reads in biological model input parameters.            !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_biology
      USE mod_ncparam
      USE mod_scalars
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
      integer :: Npts, Nval
      integer :: iTrcStr, iTrcEnd
      integer :: i, ifield, igrid, is, itracer, itrc, ng, nline, status

      integer :: decode_line, load_i, load_l, load_r

      logical, dimension(NBT,Ngrids) :: Ltrc

#ifdef ICE_BIO      
      logical, dimension(NIB,Ngrids) :: LtrcI
#endif
      real(r8), dimension(NBT,Ngrids) :: Rbio

      real(r8), dimension(200) :: Rval

      character (len=40 ) :: KeyWord
      character (len=256) :: line
      character (len=256), dimension(200) :: Cval
!
!-----------------------------------------------------------------------
!  Initialize.
!-----------------------------------------------------------------------
!
      igrid=1                            ! nested grid counter
      itracer=0                          ! LBC tracer counter
      iTrcStr=1                          ! first LBC tracer to process
      iTrcEnd=NBT                        ! last  LBC tracer to process
      nline=0                            ! LBC multi-line counter

! ==================================================================== !
! READ Bering Sea BEST_NPZ PARAMS
! ==================================================================== !

      DO WHILE (.TRUE.)
        READ (inp,'(a)',ERR=10,END=20) line
        status=decode_line(line, KeyWord, Nval, Cval, Rval)
        IF (status.gt.0) THEN
          IF (TRIM(KeyWord).eq.'Lbiology') THEN
            Npts=load_l(Nval, Cval, Ngrids, Lbiology)
          ELSE IF (TRIM(KeyWord).eq.'BioIter') THEN
            Npts=load_i(Nval, Rval, Ngrids, BioIter)
          ELSE IF (TRIM(KeyWord).eq.'PARfrac') THEN
            Npts=load_r(Nval, Rval, Ngrids, Parfrac)
!----------------------------------
!  Vertical mixing tuning parameter
!----------------------------------
          ELSE IF (TRIM(KeyWord).eq.'VertMixIncr') THEN 
            Npts=load_r(Nval, Rval, 1, VertMixIncr)         
!------------------
!  Bio- conversions
!------------------
          ELSE IF (TRIM(KeyWord).eq.'xi') THEN
            Npts=load_r(Nval, Rval, 1, xi)
          ELSE IF (TRIM(KeyWord).eq.'ccr') THEN
            Npts=load_r(Nval, Rval, 1, ccr)
          ELSE IF (TRIM(KeyWord).eq.'ccrPhL') THEN 
            Npts=load_r(Nval, Rval, 1, ccrPhL)
!-------------------------
!  extinction coefficients
!-------------------------
          ELSE IF (TRIM(KeyWord).eq.'k_ext') THEN
            Npts=load_r(Nval, Rval, 1, k_ext)
          ELSE IF (TRIM(KeyWord).eq.'k_chl') THEN
            Npts=load_r(Nval, Rval, 1, k_chl)
!-------------------
!  PhS growth params
!-------------------
          ELSE IF (TRIM(KeyWord).eq.'DiS') THEN
            Npts=load_r(Nval, Rval, 1, DiS)
          ELSE IF (TRIM(KeyWord).eq.'DpS') THEN
            Npts=load_r(Nval, Rval, 1, DpS)
          ELSE IF (TRIM(KeyWord).eq.'alphaPhS') THEN
            Npts=load_r(Nval, Rval, 1, alphaPhS)
          ELSE IF (TRIM(KeyWord).eq.'psiPhS') THEN
            Npts=load_r(Nval, Rval, 1, psiPhS)
          ELSE IF (TRIM(KeyWord).eq.'k1PhS') THEN
            Npts=load_r(Nval, Rval, 1, k1PhS)
          ELSE IF (TRIM(KeyWord).eq.'k2PhS') THEN
            Npts=load_r(Nval, Rval, 1, k2PhS)
          ELSE IF (TRIM(KeyWord).eq.'aPS') THEN
            Npts=load_r(Nval, Rval, 1, aPS)
!-------------------
!  PhL growth params
!-------------------
          ELSE IF (TRIM(KeyWord).eq.'DiL') THEN
            Npts=load_r(Nval, Rval, 1, DiL)
          ELSE IF (TRIM(KeyWord).eq.'DpL') THEN
            Npts=load_r(Nval, Rval, 1, DpL)
          ELSE IF (TRIM(KeyWord).eq.'alphaPhL') THEN
            Npts=load_r(Nval, Rval, 1, alphaPhL)
          ELSE IF (TRIM(KeyWord).eq.'psiPhL') THEN
            Npts=load_r(Nval, Rval, 1, psiPhL)
          ELSE IF (TRIM(KeyWord).eq.'k1PhL') THEN
            Npts=load_r(Nval, Rval, 1, k1PhL)
          ELSE IF (TRIM(KeyWord).eq.'k2PhL') THEN
            Npts=load_r(Nval, Rval, 1, k2PhL)
          ELSE IF (TRIM(KeyWord).eq.'aPL') THEN
            Npts=load_r(Nval, Rval, 1, aPL)
!-----------------------
!  MZS preference params
!-----------------------
          ELSE IF (TRIM(KeyWord).eq.'fpPhSMZS') THEN
            Npts=load_r(Nval, Rval, 1, fpPhSMZS)
          ELSE IF (TRIM(KeyWord).eq.'fpPhLMZS') THEN
            Npts=load_r(Nval, Rval, 1, fpPhLMZS)
!-------------------------------
!  MZS growth and feeding params
!-------------------------------
          ELSE IF (TRIM(KeyWord).eq.'eMZS') THEN
            Npts=load_r(Nval, Rval, 1, eMZS)
          ELSE IF (TRIM(KeyWord).eq.'Q10MZS') THEN
            Npts=load_r(Nval, Rval, 1, Q10MZS)
          ELSE IF (TRIM(KeyWord).eq.'Q10MZST') THEN
            Npts=load_r(Nval, Rval, 1, Q10MZST)
          ELSE IF (TRIM(KeyWord).eq.'fMZS') THEN
            Npts=load_r(Nval, Rval, 1, fMZS)
          ELSE IF (TRIM(KeyWord).eq.'kMZS') THEN
            Npts=load_r(Nval, Rval, 1, kMZS)
          ELSE IF (TRIM(KeyWord).eq.'gammaMZS') THEN
            Npts=load_r(Nval, Rval, 1, gammaMZS)
!-----------------------
!  MZL preference params
!-----------------------
          ELSE IF (TRIM(KeyWord).eq.'fpPhSMZL') THEN
            Npts=load_r(Nval, Rval, 1, fpPhSMZL)
          ELSE IF (TRIM(KeyWord).eq.'fpPhLMZL') THEN
            Npts=load_r(Nval, Rval, 1, fpPhLMZL)
!-------------------------------
!  MZL growth and feeding params
!-------------------------------
          ELSE IF (TRIM(KeyWord).eq.'eMZL') THEN
            Npts=load_r(Nval, Rval, 1, eMZL)
          ELSE IF (TRIM(KeyWord).eq.'Q10MZL') THEN
            Npts=load_r(Nval, Rval, 1, Q10MZL)
          ELSE IF (TRIM(KeyWord).eq.'Q10MZLT') THEN
            Npts=load_r(Nval, Rval, 1, Q10MZLT)
          ELSE IF (TRIM(KeyWord).eq.'fMZL') THEN
            Npts=load_r(Nval, Rval, 1, fMZL)
          ELSE IF (TRIM(KeyWord).eq.'kMZL') THEN
            Npts=load_r(Nval, Rval, 1, kMZL)
          ELSE IF (TRIM(KeyWord).eq.'gammaMZL') THEN
            Npts=load_r(Nval, Rval, 1, gammaMZL)
!-----------------------
!  Cop preference params
!-----------------------
          ELSE IF (TRIM(KeyWord).eq.'fpPhSCop') THEN
            Npts=load_r(Nval, Rval, 1, fpPhSCop)
          ELSE IF (TRIM(KeyWord).eq.'fpPhLCop') THEN
            Npts=load_r(Nval, Rval, 1, fpPhLCop)
          ELSE IF (TRIM(KeyWord).eq.'fpMZSCop') THEN
            Npts=load_r(Nval, Rval, 1, fpMZSCop)
          ELSE IF (TRIM(KeyWord).eq.'fpMZLCop') THEN
            Npts=load_r(Nval, Rval, 1, fpMZLCop)
!-------------------------------
!  Cop growth and feeding params
!-------------------------------
          ELSE IF (TRIM(KeyWord).eq.'eCop') THEN
            Npts=load_r(Nval, Rval, 1, eCop)
          ELSE IF (TRIM(KeyWord).eq.'Q10Cop') THEN
            Npts=load_r(Nval, Rval, 1, Q10Cop)
          ELSE IF (TRIM(KeyWord).eq.'Q10CopT') THEN
            Npts=load_r(Nval, Rval, 1, Q10CopT)
          ELSE IF (TRIM(KeyWord).eq.'fCop') THEN
            Npts=load_r(Nval, Rval, 1, fCop)
          ELSE IF (TRIM(KeyWord).eq.'gammaCop') THEN
            Npts=load_r(Nval, Rval, 1, gammaCop)
          ELSE IF (TRIM(KeyWord).eq.'kCop') THEN
            Npts=load_r(Nval, Rval, 1, kCop)
!-----------------------
!  NCa preference params
!-----------------------
          ELSE IF (TRIM(KeyWord).eq.'fpPhSNCa') THEN
            Npts=load_r(Nval, Rval, 1, fpPhSNCa)
          ELSE IF (TRIM(KeyWord).eq.'fpPhLNCa') THEN
            Npts=load_r(Nval, Rval, 1, fpPhLNCa)
          ELSE IF (TRIM(KeyWord).eq.'fpMZSNCa') THEN
            Npts=load_r(Nval, Rval, 1, fpMZSNCa)
          ELSE IF (TRIM(KeyWord).eq.'fpMZLNCa') THEN
            Npts=load_r(Nval, Rval, 1, fpMZLNCa)
!-------------------------------
!  NCa growth and feeding params
!-------------------------------
          ELSE IF (TRIM(KeyWord).eq.'eNCa') THEN
            Npts=load_r(Nval, Rval, 1, eNCa)
          ELSE IF (TRIM(KeyWord).eq.'Q10NCa') THEN
            Npts=load_r(Nval, Rval, 1, Q10NCa)
          ELSE IF (TRIM(KeyWord).eq.'Q10NCaT') THEN
            Npts=load_r(Nval, Rval, 1, Q10NCaT)
          ELSE IF (TRIM(KeyWord).eq.'fNCa') THEN
            Npts=load_r(Nval, Rval, 1, fNCa)
          ELSE IF (TRIM(KeyWord).eq.'gammaNCa') THEN
            Npts=load_r(Nval, Rval, 1, gammaNCa)
          ELSE IF (TRIM(KeyWord).eq.'kNCa') THEN
            Npts=load_r(Nval, Rval, 1, kNCa)
!-----------------------
!  Eup preference params
!-----------------------
          ELSE IF (TRIM(KeyWord).eq.'fpPhSEup') THEN
            Npts=load_r(Nval, Rval, 1, fpPhSEup)
          ELSE IF (TRIM(KeyWord).eq.'fpPhLEup') THEN
            Npts=load_r(Nval, Rval, 1, fpPhLEup)
          ELSE IF (TRIM(KeyWord).eq.'fpMZSEup') THEN
            Npts=load_r(Nval, Rval, 1, fpMZSEup)
          ELSE IF (TRIM(KeyWord).eq.'fpMZLEup') THEN
            Npts=load_r(Nval, Rval, 1, fpMZLEup)
          ELSE IF (TRIM(KeyWord).eq.'fpCopEup') THEN
            Npts=load_r(Nval, Rval, 1, fpCopEup)
!-------------------------------
!  Eup growth and feeding params
!-------------------------------
          ELSE IF (TRIM(KeyWord).eq.'eEup') THEN
            Npts=load_r(Nval, Rval, 1, eEup)
          ELSE IF (TRIM(KeyWord).eq.'Q10Eup') THEN
            Npts=load_r(Nval, Rval, 1, Q10Eup)
          ELSE IF (TRIM(KeyWord).eq.'Q10EupT') THEN
            Npts=load_r(Nval, Rval, 1, Q10EupT)
          ELSE IF (TRIM(KeyWord).eq.'fEup') THEN
            Npts=load_r(Nval, Rval, 1, fEup)
          ELSE IF (TRIM(KeyWord).eq.'gammaEup') THEN
            Npts=load_r(Nval, Rval, 1, gammaEup)
          ELSE IF (TRIM(KeyWord).eq.'kEup') THEN
            Npts=load_r(Nval, Rval, 1, kEup)
#if defined JELLY
!--------------------------
!  Jellyfish  param
!--------------------------
          ELSE IF (TRIM(KeyWord).eq.'fpCopJel') THEN
            Npts=load_r(Nval, Rval, 1, fpCopJel)
          ELSE IF (TRIM(KeyWord).eq.'fpEupJel') THEN
            Npts=load_r(Nval, Rval, 1, fpEupJel)
          ELSE IF (TRIM(KeyWord).eq.'fpNCaJel') THEN
            Npts=load_r(Nval, Rval, 1, fpNCaJel)
          ELSE IF (TRIM(KeyWord).eq.'eJel') THEN
            Npts=load_r(Nval, Rval, 1, eJel)
          ELSE IF (TRIM(KeyWord).eq.'Q10Jelr') THEN
            Npts=load_r(Nval, Rval, 1, Q10Jelr) 
          ELSE IF (TRIM(KeyWord).eq.'Q10JelTr') THEN
            Npts=load_r(Nval, Rval, 1, Q10JelTr)
          ELSE IF (TRIM(KeyWord).eq.'Q10Jele') THEN
            Npts=load_r(Nval, Rval, 1, Q10Jele) 
          ELSE IF (TRIM(KeyWord).eq.'Q10JelTe') THEN
            Npts=load_r(Nval, Rval, 1, Q10JelTe)
          ELSE IF (TRIM(KeyWord).eq.'gammaJel') THEN
            Npts=load_r(Nval, Rval, 1, gammaJel) 
          ELSE IF (TRIM(KeyWord).eq.'mpredJel') THEN
            Npts=load_r(Nval, Rval, 1, mpredJel) 
          ELSE IF (TRIM(KeyWord).eq.'respJel') THEN
            Npts=load_r(Nval, Rval, 1,respJel)
          ELSE IF (TRIM(KeyWord).eq.'bmJ') THEN
            Npts=load_r(Nval, Rval, 1,bmJ)
          ELSE IF (TRIM(KeyWord).eq.'ktbmJ') THEN
            Npts=load_r(Nval, Rval, 1,ktbmJ)
          ELSE IF (TRIM(KeyWord).eq.'TrefJ') THEN
            Npts=load_r(Nval, Rval, 1,TrefJ)
          ELSE IF (TRIM(KeyWord).eq.'fJel') THEN
            Npts=load_r(Nval, Rval, 1,fJel)
#endif

!--------------------------
!  Phytoplankton senescence
!--------------------------
          ELSE IF (TRIM(KeyWord).eq.'mPhS') THEN
            Npts=load_r(Nval, Rval, 1, mPhS)
          ELSE IF (TRIM(KeyWord).eq.'mPhL') THEN
            Npts=load_r(Nval, Rval, 1, mPhL)
          ELSE IF (TRIM(KeyWord).eq.'NcritPhS') THEN
            Npts=load_r(Nval, Rval, 1, NcritPhS)
          ELSE IF (TRIM(KeyWord).eq.'minmPhL') THEN
            Npts=load_r(Nval, Rval, 1, minmPhL)
          ELSE IF (TRIM(KeyWord).eq.'maxmPhL') THEN
            Npts=load_r(Nval, Rval, 1, maxmPhL)
          ELSE IF (TRIM(KeyWord).eq.'NcritPhL') THEN
            Npts=load_r(Nval, Rval, 1, NcritPhL)
!-----------------------
!  Zoopkankton mortality
!-----------------------
          ELSE IF (TRIM(KeyWord).eq.'mMZS') THEN
            Npts=load_r(Nval, Rval, 1, mMZS)
          ELSE IF (TRIM(KeyWord).eq.'mMZL') THEN
            Npts=load_r(Nval, Rval, 1, mMZL)
          ELSE IF (TRIM(KeyWord).eq.'mCop') THEN
            Npts=load_r(Nval, Rval, 1, mCop)
          ELSE IF (TRIM(KeyWord).eq.'mNCa') THEN
            Npts=load_r(Nval, Rval, 1, mNCa)
          ELSE IF (TRIM(KeyWord).eq.'mEup') THEN
            Npts=load_r(Nval, Rval, 1, mEup)
!-------------------
!  predation closure
!-------------------
          ELSE IF (TRIM(KeyWord).eq.'mpredMZL') THEN
            Npts=load_r(Nval, Rval, 1, mpredMZL)
          ELSE IF (TRIM(KeyWord).eq.'mpredCop') THEN
            Npts=load_r(Nval, Rval, 1, mpredCop)
          ELSE IF (TRIM(KeyWord).eq.'mpredNCa') THEN
            Npts=load_r(Nval, Rval, 1, mpredNCa)
          ELSE IF (TRIM(KeyWord).eq.'mpredEup') THEN
            Npts=load_r(Nval, Rval, 1, mpredEup)
!--------------------------------
!  sinking 
!--------------------------------
        
          ELSE IF (TRIM(KeyWord).eq.'wPhS') THEN
            Npts=load_r(Nval, Rval, 1, wPhS)
          ELSE IF (TRIM(KeyWord).eq.'wPhL') THEN
            Npts=load_r(Nval, Rval, 1, wPhL)
          ELSE IF (TRIM(KeyWord).eq.'wDet') THEN
            Npts=load_r(Nval, Rval, 1, wDet)
          ELSE IF (TRIM(KeyWord).eq.'wDetF') THEN
            Npts=load_r(Nval, Rval, 1, wDetF)
!------------------------
!  Respiration parameters
!------------------------
          ELSE IF (TRIM(KeyWord).eq.'respMZS') THEN
            Npts=load_r(Nval, Rval, 1, respMZS)
          ELSE IF (TRIM(KeyWord).eq.'respMZL') THEN
            Npts=load_r(Nval, Rval, 1, respMZL)
          ELSE IF (TRIM(KeyWord).eq.'respPhS') THEN
            Npts=load_r(Nval, Rval, 1, respPhS)
          ELSE IF (TRIM(KeyWord).eq.'respPhL') THEN
            Npts=load_r(Nval, Rval, 1, respPhL)
          ELSE IF (TRIM(KeyWord).eq.'respCop') THEN
            Npts=load_r(Nval, Rval, 1, respCop)
          ELSE IF (TRIM(KeyWord).eq.'respNCa') THEN
            Npts=load_r(Nval, Rval, 1, respNCa)
          ELSE IF (TRIM(KeyWord).eq.'respEup') THEN
            Npts=load_r(Nval, Rval, 1, respEup)

          ELSE IF (TRIM(KeyWord).eq.'ktbmC') THEN
            Npts=load_r(Nval, Rval, 1,ktbmC)
          ELSE IF (TRIM(KeyWord).eq.'TrefC') THEN
            Npts=load_r(Nval, Rval, 1,TrefC)
          ELSE IF (TRIM(KeyWord).eq.'ktbmN') THEN
            Npts=load_r(Nval, Rval, 1,ktbmN)
          ELSE IF (TRIM(KeyWord).eq.'TrefN') THEN
            Npts=load_r(Nval, Rval, 1,TrefN)

          ELSE IF (TRIM(KeyWord).eq.'ktbmE') THEN
            Npts=load_r(Nval, Rval, 1,ktbmE)
          ELSE IF (TRIM(KeyWord).eq.'TrefE') THEN
            Npts=load_r(Nval, Rval, 1,TrefE)

          ELSE IF (TRIM(KeyWord).eq.'TmaxPhS') THEN
            Npts=load_r(Nval, Rval, 1, TmaxPhS)
          ELSE IF (TRIM(KeyWord).eq.'TminPhS') THEN
            Npts=load_r(Nval, Rval, 1,TminPhS)
          ELSE IF (TRIM(KeyWord).eq.'Topt_PhS') THEN
            Npts=load_r(Nval, Rval, 1, Topt_PhS)
          ELSE IF (TRIM(KeyWord).eq.'KtBm_PhS') THEN
            Npts=load_r(Nval, Rval, 1,KtBm_PhS)
          ELSE IF (TRIM(KeyWord).eq.'TmaxPhL') THEN
            Npts=load_r(Nval, Rval, 1, TmaxPhL)
          ELSE IF (TRIM(KeyWord).eq.'TminPhL') THEN
            Npts=load_r(Nval, Rval, 1,TminPhL)
          ELSE IF (TRIM(KeyWord).eq.'Topt_PhL') THEN
            Npts=load_r(Nval, Rval, 1, Topt_PhL)
          ELSE IF (TRIM(KeyWord).eq.'KtBm_PhL') THEN
            Npts=load_r(Nval, Rval, 1,KtBm_PhL)

          ELSE IF (TRIM(KeyWord).eq.'TmaxMZS') THEN
            Npts=load_r(Nval, Rval, 1, TmaxMZS)
          ELSE IF (TRIM(KeyWord).eq.'KtBm_MZS') THEN
            Npts=load_r(Nval, Rval, 1,KtBm_MZS)
          ELSE IF (TRIM(KeyWord).eq.'TmaxMZL') THEN
            Npts=load_r(Nval, Rval, 1, TmaxMZL)
          ELSE IF (TRIM(KeyWord).eq.'KtBm_MZL') THEN
            Npts=load_r(Nval, Rval, 1,KtBm_MZL)
!------------------------
!  Iron climatology terms
!------------------------
          ELSE IF (TRIM(KeyWord).eq.'Feinlo') THEN
            Npts=load_r(Nval, Rval, 1, Feinlo)
          ELSE IF (TRIM(KeyWord).eq.'Feinhi') THEN
            Npts=load_r(Nval, Rval, 1, Feinhi)
          ELSE IF (TRIM(KeyWord).eq.'Feinh') THEN
            Npts=load_r(Nval, Rval, 1, Feinh)
          ELSE IF (TRIM(KeyWord).eq.'Feofflo') THEN
            Npts=load_r(Nval, Rval, 1, Feofflo)
          ELSE IF (TRIM(KeyWord).eq.'Feoffhi') THEN
            Npts=load_r(Nval, Rval, 1, Feoffhi)
          ELSE IF (TRIM(KeyWord).eq.'Feoffh') THEN
            Npts=load_r(Nval, Rval, 1, Feoffh)
!-----------------------
!  Iron limitation terms
!-----------------------
          ELSE IF (TRIM(KeyWord).eq.'kfePhS') THEN
            Npts=load_r(Nval, Rval, 1, kfePhS)
          ELSE IF (TRIM(KeyWord).eq.'kfePhL') THEN
            Npts=load_r(Nval, Rval, 1, kfePhL)
          ELSE IF (TRIM(KeyWord).eq.'FeC') THEN
            Npts=load_r(Nval, Rval, 1, FeC)
!----------
!  Diapause
!----------
          ELSE IF (TRIM(KeyWord).eq.'NCmaxz') THEN
            Npts=load_r(Nval, Rval, 1, NCmaxz)
          ELSE IF (TRIM(KeyWord).eq.'wNCrise') THEN
            Npts=load_r(Nval, Rval, 1, wNCrise)
          ELSE IF (TRIM(KeyWord).eq.'wNCsink') THEN
            Npts=load_r(Nval, Rval, 1, wNCsink)
          ELSE IF (TRIM(KeyWord).eq.'RiseStart') THEN
            Npts=load_r(Nval, Rval, 1, RiseStart)
          ELSE IF (TRIM(KeyWord).eq.'RiseEnd') THEN
            Npts=load_r(Nval, Rval, 1, RiseEnd)
          ELSE IF (TRIM(KeyWord).eq.'SinkStart') THEN
            Npts=load_r(Nval, Rval, 1, SinkStart)
          ELSE IF (TRIM(KeyWord).eq.'SinkEnd') THEN
            Npts=load_r(Nval, Rval, 1, SinkEnd)


!-----------------------------------
!Remineralization and Nitrification
!----------------------------------

          ELSE IF (TRIM(KeyWord).eq.'regen') THEN
            Npts=load_r(Nval, Rval, 1, regen)
          ELSE IF (TRIM(KeyWord).eq.'dgrad') THEN
            Npts=load_r(Nval, Rval, 1, dgrad)
          ELSE IF (TRIM(KeyWord).eq.'Pv0') THEN
            Npts=load_r(Nval, Rval, 1, Pv0)
          ELSE IF (TRIM(KeyWord).eq.'PvT') THEN
            Npts=load_r(Nval, Rval, 1, PvT)
          ELSE IF (TRIM(KeyWord).eq.'Nitr0') THEN
            Npts=load_r(Nval, Rval, 1,Nitr0)
          ELSE IF (TRIM(KeyWord).eq.'KnT') THEN
            Npts=load_r(Nval, Rval, 1,KnT)
          ELSE IF (TRIM(KeyWord).eq.'ToptNtr') THEN
            Npts=load_r(Nval, Rval, 1,ToptNtr)
          ELSE IF (TRIM(KeyWord).eq.'ktntr') THEN
            Npts=load_r(Nval, Rval, 1,ktntr)
          ELSE IF (TRIM(KeyWord).eq.'KNH4Nit') THEN
            Npts=load_r(Nval, Rval, 1,KNH4Nit)
          ELSE IF (TRIM(KeyWord).eq.'tI0') THEN
            Npts=load_r(Nval, Rval, 1,tI0)
          ELSE IF (TRIM(KeyWord).eq.'KI') THEN
            Npts=load_r(Nval, Rval, 1,KI)  

!-------------------
!Benthic Parameters
!-------------------
#ifdef BENTHIC
          ELSE IF (TRIM(KeyWord).eq.'Hout(idBeTvar)') THEN
            Npts=load_l(Nval, Cval,NBEN*Ngrids, Ltrc)
            DO ng=1,Ngrids
              DO itrc=1,NBEN
                i=idBeTvar(idben(itrc))
                Hout(i,ng)=Ltrc(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'iremin') THEN
            Npts=load_r(Nval, Rval, 1, iremin)
          ELSE IF (TRIM(KeyWord).eq.'q10') THEN
            Npts=load_r(Nval, Rval, 1, q10)
          ELSE IF (TRIM(KeyWord).eq.'q10r') THEN
            Npts=load_r(Nval, Rval, 1, q10r)
          ELSE IF (TRIM(KeyWord).eq.'Rup') THEN
            Npts=load_r(Nval, Rval, 1, Rup)
          ELSE IF (TRIM(KeyWord).eq.'KupD') THEN
            Npts=load_r(Nval, Rval, 1, KupD)
          ELSE IF (TRIM(KeyWord).eq.'KupP') THEN
            Npts=load_r(Nval, Rval, 1, KupP)
          ELSE IF (TRIM(KeyWord).eq.'LupD') THEN
            Npts=load_r(Nval, Rval, 1, LupD)
          ELSE IF (TRIM(KeyWord).eq.'LupP') THEN
            Npts=load_r(Nval, Rval, 1, LupP)
          ELSE IF (TRIM(KeyWord).eq.'Qres') THEN
            Npts=load_r(Nval, Rval, 1, Qres)
          ELSE IF (TRIM(KeyWord).eq.'Rres') THEN
            Npts=load_r(Nval, Rval, 1, Rres)
          ELSE IF (TRIM(KeyWord).eq.'rmort') THEN
            Npts=load_r(Nval, Rval, 1, rmort)
          ELSE IF (TRIM(KeyWord).eq.'eex') THEN
            Npts=load_r(Nval, Rval, 1, eex)
          ELSE IF (TRIM(KeyWord).eq.'eexD') THEN
            Npts=load_r(Nval, Rval, 1, eexD)
          ELSE IF (TRIM(KeyWord).eq.'prefD') THEN
            Npts=load_r(Nval, Rval, 1, prefD)
          ELSE IF (TRIM(KeyWord).eq.'prefPL') THEN
            Npts=load_r(Nval, Rval, 1, prefPL)
          ELSE IF (TRIM(KeyWord).eq.'prefPS') THEN
            Npts=load_r(Nval, Rval, 1, prefPS)
          ELSE IF (TRIM(KeyWord).eq.'T0ben') THEN
            Npts=load_r(Nval, Rval, 1, T0ben)
          ELSE IF (TRIM(KeyWord).eq.'T0benr') THEN
            Npts=load_r(Nval, Rval, 1, T0benr)
          ELSE IF (TRIM(KeyWord).eq.'BenPred') THEN
            Npts=load_r(Nval, Rval, 1, BenPred)
          ELSE IF (TRIM(KeyWord).eq.'bmB') THEN
            Npts=load_r(Nval, Rval, 1,bmB)
          ELSE IF (TRIM(KeyWord).eq.'ktbmB') THEN
            Npts=load_r(Nval, Rval, 1,ktbmB)
          ELSE IF (TRIM(KeyWord).eq.'TrefB') THEN
            Npts=load_r(Nval, Rval, 1,TrefB)
#endif

!--------------
!ice bio params
!--------------
#ifdef ICE_BIO
# ifdef CLIM_ICE_1D
          ELSE IF (TRIM(KeyWord).eq.'Hout(idIceBvar)') THEN
            Npts=load_l(Nval, Cval,NIB*Ngrids, Ltrc)
            DO ng=1,Ngrids
              DO itrc=1,NIB
                i=idIceBvar(idice(itrc))
                Hout(i,ng)=Ltrc(itrc,ng)
              END DO
            END DO
# elif defined BERING_10K  
          ELSE IF (TRIM(KeyWord).eq.'Hout(idIcePhL)') THEN
            Npts=load_l(Nval, Cval,NIB*Ngrids, LtrcI)
            DO ng=1,Ngrids
!g            Hout(idIcePhL,ng)=LtrcI(idIcePhL,ng)
              Hout(idIcePhL,ng)=LtrcI(1,ng)
            END DO
          ELSE IF (TRIM(KeyWord).eq.'Hout(idIceNO3)') THEN
!            Npts=load_l(Nval, Cval,NIB*Ngrids, LtrcI)
             DO ng=1,Ngrids
               Hout(idIceNO3,ng)=LtrcI(2,ng)
!g             Hout(idIceNO3,ng)=LtrcI(idIceNO3,ng)
             END DO
           ELSE IF (TRIM(KeyWord).eq.'Hout(idIceNH4)') THEN
!            Npts=load_l(Nval, Cval,NIB*Ngrids, LtrcI)
             DO ng=1,Ngrids
!g             Hout(idIceNH4,ng)=LtrcI(idIceNH4,ng)
               Hout(idIceNH4,ng)=LtrcI(3,ng)
             END DO
           ELSE IF (TRIM(KeyWord).eq.'Hout(idIceLog)') THEN
!            Npts=load_l(Nval, Cval,NIB*Ngrids, LtrcI)
             DO ng=1,Ngrids
!              Hout(idIceLog,ng)=LtrcI(idIceLog,ng)
               Hout(idIceLog,ng)=LtrcI(4,ng)
             END DO
# endif
           ELSE IF (TRIM(KeyWord).eq.'alphaIb') THEN
             Npts=load_r(Nval, Rval, 1,alphaIb)
           ELSE IF (TRIM(KeyWord).eq.'betaI') THEN
             Npts=load_r(Nval, Rval, 1,betaI)
           ELSE IF (TRIM(KeyWord).eq.'inhib') THEN
             Npts=load_r(Nval, Rval, 1,inhib)
           ELSE IF (TRIM(KeyWord).eq.'ksnut1') THEN
             Npts=load_r(Nval, Rval, 1,ksnut1)
           ELSE IF (TRIM(KeyWord).eq.'ksnut2') THEN
             Npts=load_r(Nval, Rval, 1,ksnut2)
           ELSE IF (TRIM(KeyWord).eq.'mu0') THEN
             Npts=load_r(Nval, Rval, 1,mu0)
           ELSE IF (TRIM(KeyWord).eq.'R0i') THEN
             Npts=load_r(Nval, Rval, 1,R0i)
           ELSE IF (TRIM(KeyWord).eq.'rg0') THEN
             Npts=load_r(Nval, Rval, 1,rg0)
           ELSE IF (TRIM(KeyWord).eq.'rg') THEN
             Npts=load_r(Nval, Rval, 1,rg)
           ELSE IF (TRIM(KeyWord).eq.'annit') THEN
             Npts=load_r(Nval, Rval, 1,annit)
           ELSE IF (TRIM(KeyWord).eq.'aidz') THEN
             Npts=load_r(Nval, Rval, 1,aidz)
#endif

          ELSE IF (TRIM(KeyWord).eq.'TNU2') THEN
            Npts=load_r(Nval, Rval, NBT*Ngrids, Rbio)
            DO ng=1,Ngrids
              DO itrc=1,NBT
                i=idbio(itrc)
                tnu2(i,ng)=Rbio(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'TNU4') THEN
            Npts=load_r(Nval, Rval, NBT*Ngrids, Rbio)
            DO ng=1,Ngrids
              DO itrc=1,NBT
                i=idbio(itrc)
                tnu4(i,ng)=Rbio(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'AKT_BAK') THEN
            Npts=load_r(Nval, Rval, NBT*Ngrids, Rbio)
            DO ng=1,Ngrids
              DO itrc=1,NBT
                i=idbio(itrc)
                Akt_bak(i,ng)=Rbio(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'TNUDG') THEN
            Npts=load_r(Nval, Rval, NBT*Ngrids, Rbio)

! Hard wired this in for now as wasnt reading correctly from input file    
            Rbio(15,1)=360_r8

            DO ng=1,Ngrids
              DO itrc=1,NBT
                i=idbio(itrc)
                Tnudg(i,ng)=Rbio(itrc,ng)
              END DO
            END DO
            CASE ('LBC(isTvar)')
              IF (itracer.lt.NBT) THEN
                itracer=itracer+1
              ELSE
                itracer=1                      ! next nested grid
              END IF
              ifield=isTvar(idbio(itracer))
              Npts=load_lbc(Nval, Cval, line, nline, ifield, igrid,     &
     &                      idbio(iTrcStr), idbio(iTrcEnd),             &
     &                      Vname(1,idTvar(idbio(itracer))), LBC)
#ifdef TCLIMATOLOGY
            CASE ('LtracerCLM')
              Npts=load_l(Nval, Cval, NBT*Ngrids, Ltrc)
              DO ng=1,Ngrids
                DO itrc=1,NBT
                  i=idbio(itrc)
                  LtracerCLM(i,ng)=Ltrc(itrc,ng)
                END DO
              END DO
#endif
#ifdef TS_PSOURCE
            CASE ('LtracerSrc')
              Npts=load_l(Nval, Cval, NBT*Ngrids, Ltrc)
              DO ng=1,Ngrids
                DO itrc=1,NBT
                  i=idbio(itrc)
                  LtracerSrc(i,ng)=Ltrc(itrc,ng)
                END DO
              END DO
#endif
          ELSE IF (TRIM(KeyWord).eq.'Hout(idTvar)') THEN
            Npts=load_l(Nval, Cval, NBT*Ngrids, Ltrc)
            DO ng=1,Ngrids
              DO itrc=1,NBT
                i=idTvar(idbio(itrc))
                Hout(i,ng)=Ltrc(itrc,ng)
              END DO
            END DO
#ifdef AVERAGES
          ELSE IF (TRIM(KeyWord).eq.'Aout(idTvar)') THEN
            Npts=load_l(Nval, Cval, NBT*Ngrids, Ltrc)
            DO ng=1,Ngrids
              DO itrc=1,NBT
                i=idTvar(idbio(itrc))
                Aout(i,ng)=Ltrc(itrc,ng)
              END DO
            END DO
# ifdef BENTHIC
          ELSE IF (TRIM(KeyWord).eq.'Aout(idBeTvar)') THEN
            Npts=load_l(Nval, Cval,NBEN*Ngrids, Ltrc)
            DO ng=1,Ngrids
              DO itrc=1,NBEN
                i=idBeTvar(idben(itrc))
                Aout(i,ng)=Ltrc(itrc,ng)
              END DO
            END DO
# endif
# ifdef ICE_BIO
#  ifdef CLIM_ICE_1D
          ELSE IF (TRIM(KeyWord).eq.'Aout(idIceBvar)') THEN
            Npts=load_l(Nval, Cval,NIB*Ngrids, Ltrc)
            DO ng=1,Ngrids
              DO itrc=1,NIB
                i=idIceBvar(idice(itrc))
                Aout(i,ng)=Ltrc(itrc,ng)
              END DO
            END DO
#  elif defined BERING_10K  
          ELSE IF (TRIM(KeyWord).eq.'Aout(idIcePhL)') THEN
            Npts=load_l(Nval, Cval,NIB*Ngrids, LtrcI)
            DO ng=1,Ngrids
!g            Aout(idIcePhL,ng)=LtrcI(idIcePhL,ng)
              Aout(idIcePhL,ng)=LtrcI(1,ng)
            END DO
          ELSE IF (TRIM(KeyWord).eq.'Aout(idIceNO3)') THEN
!           Npts=load_l(Nval, Cval,NIB*Ngrids, LtrcI)
            DO ng=1,Ngrids
              Aout(idIceNO3,ng)=LtrcI(2,ng)
!g            Aout(idIceNO3,ng)=LtrcI(idIceNO3,ng)
            END DO
          ELSE IF (TRIM(KeyWord).eq.'Aout(idIceNH4)') THEN
!           Npts=load_l(Nval, Cval,NIB*Ngrids, LtrcI)
            DO ng=1,Ngrids
!g            Aout(idIceNH4,ng)=LtrcI(idIceNH4,ng)
              Aout(idIceNH4,ng)=LtrcI(3,ng)
            END DO
          ELSE IF (TRIM(KeyWord).eq.'Aout(idIceLog)') THEN
!           Npts=load_l(Nval, Cval,NIB*Ngrids, LtrcI)
            DO ng=1,Ngrids
!             Aout(idIceLog,ng)=LtrcI(idIceLog,ng)
              Aout(idIceLog,ng)=LtrcI(4,ng)
            END DO
#  endif
# endif
#endif
#ifdef STATIONARY
          ELSE IF (TRIM(KeyWord).eq.'Hout(idTSvar)') THEN
            Npts=load_l(Nval, Cval, NTS*Ngrids, Ltrc)
            DO ng=1,Ngrids
              DO itrc=1, NBTS   
                i=idTSvar(itrc)
                Hout(i,ng)=Ltrc(itrc,ng)
              END DO
            END DO
         
#endif
#ifdef BIOFLUX
          ELSE IF (TRIM(KeyWord).eq.'Hout(idTBFvar)') THEN
            Npts=load_l(Nval, Cval, NT*Ngrids, Ltrc)
            DO ng=1,Ngrids
!              DO itrc=1, NT(ng)   
                i=idTBFvar(iBF)
                Hout(i,ng)=Ltrc(iBF,ng)
!              END DO
            END DO
         
#endif
#ifdef STATIONARY2
          ELSE IF (TRIM(KeyWord).eq.'Hout(idTS2var)') THEN
            Npts=load_l(Nval, Cval, NTS2*Ngrids, Ltrc)
            DO ng=1,Ngrids
              DO itrc=1, NBTS2   
                i=idTS2var(itrc)
                Hout(i,ng)=Ltrc(itrc,ng)
              END DO
            END DO
         
#endif

#ifdef PROD3
          ELSE IF (TRIM(KeyWord).eq.'Hout(idPT3var)') THEN
            Npts=load_l(Nval, Cval, NPT3*Ngrids, Ltrc)
            DO ng=1,Ngrids
              DO itrc=1, NBPT3   
                i=idPT3var(itrc)
                Hout(i,ng)=Ltrc(itrc,ng)
              END DO
            END DO
          
#endif
#ifdef PROD2
          ELSE IF (TRIM(KeyWord).eq.'Hout(idPT2var)') THEN
            Npts=load_l(Nval, Cval, NPT2*Ngrids, Ltrc)
            DO ng=1,Ngrids
              DO itrc=1, NBPT2   
                i=idPT2var(itrc)
                Hout(i,ng)=Ltrc(itrc,ng)
              END DO
            END DO
#endif

          END IF 

          IF ( Lwrite .and.                                             &
     &         TRIM(KeyWord).ne.'TNU2' .and.                            &
     &         TRIM(KeyWord).ne.'TNU4' .and.                            &
     &         TRIM(KeyWord).ne.'AKT_BAK' .and.                         &
     &         TRIM(KeyWord).ne.'TNUDG' .and.                           &
     &         TRIM(KeyWord).ne.'Hout(idTvar)' ) THEN
            write(6,'(a15,i3,20e12.5)') TRIM(KeyWord),Nval,Rval(1:Nval)
          END IF
        END IF
      END DO
  10  IF (Master) WRITE (out,30) line
      exit_flag=4
      RETURN
  20  CLOSE (inp)
!

!-----------------------------------------------------------------------
!  Report input parameters.
!-----------------------------------------------------------------------
!
      IF (Lwrite) THEN
        DO ng=1,Ngrids
          IF (Lbiology(ng)) THEN
            WRITE (out,40) ng
            WRITE (out,50) BioIter(ng), 'BioIter',                      &
     &            'Number of iterations for nonlinear convergence.'
#ifdef TS_DIF2
            DO itrc=1,NBT
              i=idbio(itrc)
              WRITE (out,90) tnu2(i,ng), 'tnu2', i,                     &
     &              'Horizontal, harmonic mixing coefficient (m2/s)',   &
     &              'for tracer ', i, TRIM(Vname(1,idTvar(i)))
            END DO
#endif
#ifdef TS_DIF4
            DO itrc=1,NBT
              i=idbio(itrc)
              WRITE (out,90) tnu4(i,ng), 'tnu4', i,                     &
     &              'Horizontal, biharmonic mixing coefficient (m4/s)', &
     &              'for tracer ', i, TRIM(Vname(1,idTvar(i)))
            END DO
#endif
            DO itrc=1,NBT
              i=idbio(itrc)
              WRITE(out,90) Akt_bak(i,ng), 'Akt_bak', i,                &
     &             'Background vertical mixing coefficient (m2/s)',     &
     &             'for tracer ', i, TRIM(Vname(1,idTvar(i)))
            END DO
            DO itrc=1,NBT
              i=idbio(itrc)
              WRITE (out,90) Tnudg(i,ng), 'Tnudg', i,                   &
     &              'Nudging/relaxation time scale (days)',             &
     &              'for tracer ', i, TRIM(Vname(1,idTvar(i)))
            END DO
#ifdef TCLIMATOLOGY
            DO itrc=1,NBT
              i=idbio(itrc)
	      IF (LtracerCLM(i,ng)) THEN
                WRITE (out,110) LtracerCLM(i,ng), 'LtracerCLM', i,      &
     &              'Turning ON processing of climatology tracer ', i,  &
     &              TRIM(Vname(1,idTvar(i)))
              ELSE
                WRITE (out,110) LtracerCLM(i,ng), 'LtracerCLM', i,      &
     &              'Turning OFF processing of climatology tracer ', i, &
     &              TRIM(Vname(1,idTvar(i)))
              END IF
            END DO
            DO itrc=1,NBT
              i=idbio(itrc)
              IF (LnudgeTCLM(i,ng)) THEN
                WRITE (out,110) LnudgeTCLM(i,ng), 'LnudgeTCLM', i,      &
     &              'Turning ON  nudging of climatology tracer ', i,    &
     &              TRIM(Vname(1,idTvar(i)))
              ELSE
                WRITE (out,110) LnudgeTCLM(i,ng), 'LnudgeTCLM', i,      &
     &              'Turning OFF nudging of climatology tracer ', i,    &
     &              TRIM(Vname(1,idTvar(i)))
              END IF
            END DO
#endif
#ifdef TS_PSOURCE
            DO itrc=1,NBT
              i=idbio(itrc)
              WRITE (out,100) LtracerSrc(i,ng), 'LtracerSrc',           &
     &              i, 'Processing point sources/Sink on tracer ', i,   &
     &              TRIM(Vname(1,idTvar(i)))
            END DO
#endif
            DO itrc=1,NBT
              i=idbio(itrc)
              IF (Hout(idTvar(i),ng)) WRITE (out,60)                    &
     &            Hout(idTvar(i),ng), 'Hout(idTvar)',                   &
     &            'Write out tracer ', i, TRIM(Vname(1,idTvar(i)))
            END DO
#ifdef AVERAGES
            WRITE (out,'(1x)')
            DO itrc=1,NBT
              i=idbio(itrc)
              IF (Aout(idTvar(i),ng)) WRITE (out,110)                   &
     &            Aout(idTvar(i),ng), 'Aout(idTvar)',                   &
     &            'Write out averaged tracer ', i,                      &
     &            TRIM(Vname(1,idTvar(i)))
            END DO
#endif
          END IF
        END DO
      END IF
      DO ng=1,Ngrids
        DO itrc=1,NBT
          i=idbio(itrc)
          tnu4(i,ng)=SQRT(ABS(tnu4(i,ng)))
!
!  Compute inverse nudging coefficients (1/s) used in various tasks.
!
          IF (Tnudg(i,ng).gt.0.0_r8) THEN
            Tnudg(i,ng)=1.0_r8/(Tnudg(i,ng)*86400.0_r8)
          ELSE
            Tnudg(i,ng)=0.0_r8
          END IF
        END DO
      END DO

  30  FORMAT (/,' READ_BioPar - Error while processing line: ',/,a)
  40  FORMAT (/,/,' Biology Parameters, Grid: ',i2.2,                   &
     &        /,  ' ============================',/)
  50  FORMAT (1x,i10,2x,a,t28,a)
  60  FORMAT (10x,l1,2x,a,t28,a,i2.2,':',1x,a)
  90  FORMAT (1p,e11.4,2x,a,'(',i2.2,')',t28,a,/,t30,a,i2.2,':',1x,a)
 100  FORMAT (10x,l1,2x,a,'(',i2.2,')',t30,a,i2.2,':',1x,a)
 110  FORMAT (10x,l1,2x,a,t30,a,i2.2,':',1x,a)

      RETURN
      END SUBROUTINE read_BioPar
