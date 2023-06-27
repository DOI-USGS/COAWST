#include "cppdefs.h"
#if defined ADJOINT && defined BIOLOGY

/*
** svn $Is$
************************************************************************
**                                                                    **
** This routine computes the adjoint source and sink terms for the    **
** selected biology model.                                            **
**                                                                    **
************************************************************************
*/

# if defined BIO_FENNEL
!!#  include <ad_fennel.h>
# elif defined NEMURO
!!#  include <ad_nemuro.h>
# elif defined NPZD_FRANKS
#  include <ad_npzd_Franks.h>
# elif defined NPZD_IRON
#  include <ad_npzd_iron.h>
# elif defined NPZD_POWELL
#  include <ad_npzd_Powell.h>
# endif
#else
      MODULE ad_biology_mod
      END MODULE ad_biology_mod
#endif
