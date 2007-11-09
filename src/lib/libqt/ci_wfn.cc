/*!
  \file ci_wfn.c
  \ingroup (QT)
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <psifiles.h>

extern "C" {
	
/*!
** ci_wfn()
** This function takes a WFN string and returns
** 1 if the WFN is a CI/CAS/RAS method and
** 0 if the WFN is not a CI/CAS/RAS method.
**
** \param char *wfn: wavefunction string
** \ingroup (QT)
*/
int ci_wfn(char *wfn)
{

  if (strcmp(wfn, "CI")==0     || strcmp(wfn, "DETCAS")==0 || 
      strcmp(wfn, "CASSCF")==0 || strcmp(wfn, "RASSCF")==0 ||
      strcmp(wfn, "DETCI")==0 ) 
  {
    return(1);
  }
  else  {
    return 0;
  }
}

} /* extern "C" */