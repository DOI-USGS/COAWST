!-------------------------------------------------------------------------------
! NUOPC CPP Macros
!-------------------------------------------------------------------------------
#define CONTEXT  line=__LINE__,file=__FILE__
#define PASSTHRU msg=ESMF_LOGERR_PASSTHRU,CONTEXT
#define ESMF_STDERRORCHECK(rc) ESMF_LogFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)

!-------------------------------------------------------------------------------
! Define ESMF real kind to match Appplications single/double precision
!-------------------------------------------------------------------------------
#if defined(REAL4) || defined(SINGLE)
#define ESMF_KIND_COORD ESMF_KIND_R4
#elif defined(REAL8)
#define ESMF_KIND_COORD ESMF_KIND_R8
#else
#define ESMF_KIND_COORD ESMF_KIND_R8
#endif

