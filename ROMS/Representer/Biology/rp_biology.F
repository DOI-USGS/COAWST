#include "cppdefs.h"
#if defined TL_IOMS && defined BIOLOGY

/*
** svn $Is$
************************************************************************
**                                                                    **
** This routine computes the representers tangent linear source and   **
** sink terms for the selected biology model.                         **
**                                                                    **
************************************************************************
*/

# if defined BIO_FENNEL
!!#  include <rp_fennel.h>
# elif defined NEMURO
!!#  include <rp_nemuro.h>
# elif defined NPZD_FRANKS
#  include <rp_npzd_Franks.h>
# elif defined NPZD_IRON
#  include <rp_npzd_iron.h>
# elif defined NPZD_POWELL
#  include <rp_npzd_Powell.h>
# endif
#else
      MODULE rp_biology_mod
      END MODULE rp_biology_mod
#endif
