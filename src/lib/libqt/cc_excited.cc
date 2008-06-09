/*!
  \file
  \brief Determine if a wavefunction is a CC-excited wavefunction type
  \ingroup QT
*/

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <psifiles.h>

namespace psi {
	
/*!
** cc_excited(): This function takes a WFN string and returns 1 if the WFN 
** is an excited-state method and 0 if the WFN is a ground-state method.
**
** \param *wfn = wavefunction string
** 
** Returns: 1 if an excited state method, else 0
** \ingroup QT
*/
int cc_excited(char *wfn)
{

  if ( !strcmp(wfn, "CCSD")   || !strcmp(wfn, "CCSD_T") || !strcmp(wfn, "BCCD") 
    || !strcmp(wfn, "BCCD_T") || !strcmp(wfn, "CC2")    || !strcmp(wfn, "CC3") 
    || !strcmp(wfn, "CCSD_MVD") ) {
    return 0;
  }
  else if ( !strcmp(wfn, "EOM_CCSD") || !strcmp(wfn, "LEOM_CCSD") ||
            !strcmp(wfn, "EOM_CC2")  || !strcmp(wfn, "EOM_CC3") ) {
    return 1;
  }
  else {
    printf("Invalid value of input keyword WFN: %s\n", wfn);
    exit(PSI_RETURN_FAILURE);
  }
}

}

