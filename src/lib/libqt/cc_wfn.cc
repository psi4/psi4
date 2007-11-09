/*!
  \file cc_wfn.c
  \ingroup (QT)
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <psifiles.h>

extern "C" {
	
/*!
** cc_wfn()
** This function takes a WFN string and returns
** 1 if the WFN is a CC method and
** 0 if the WFN is not a CC method.
**
**  \param char *wfn: wavefunction string
** \ingroup (QT)
*/
int cc_wfn(char *wfn)
{

  if ( !strcmp(wfn, "CCSD")     || !strcmp(wfn, "CCSD_T") || !strcmp(wfn, "BCCD") ||
       !strcmp(wfn, "BCCD_T")   || !strcmp(wfn, "CC2")    || !strcmp(wfn, "CC3")  ||
       !strcmp(wfn, "EOM_CCSD") || !strcmp(wfn, "LEOM_CCSD") ||
       !strcmp(wfn, "EOM_CC2")  || !strcmp(wfn, "EOM_CC3") ) {
			 return 1;
	}
	else {
	  return 0;
	}
}

} /* extern "C" */