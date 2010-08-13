/*!
  \file
  \brief Check if wavefunction is coupled-cluster type
  \ingroup QT
*/

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <psifiles.h>

namespace psi {
	
/*!
** cc_wfn(): Checks if the given wavefunction string is a coupled-cluster
** type and returns 1 if yes and 0 if no.  
**
** Note: "coupled-cluster type" means it is handled by PSI like the 
** coupled-cluster codes, not necessarily that it is literally a
** coupled-cluster wavefunction
**
** \param *wfn = wavefunction string
**
** Returns: 1 if the WFN is a CC method, 0 otherwise
**
** \ingroup QT
*/
int cc_wfn(const char *wfn)
{
  if ( !strcmp(wfn, "CCSD")     || !strcmp(wfn, "CCSD_T") || 
       !strcmp(wfn, "BCCD")     || !strcmp(wfn, "BCCD_T") || 
       !strcmp(wfn, "CC2")      || !strcmp(wfn, "CC3")    ||
       !strcmp(wfn, "EOM_CCSD") || !strcmp(wfn, "LEOM_CCSD") ||
       !strcmp(wfn, "EOM_CC2")  || !strcmp(wfn, "EOM_CC3") ||
       !strcmp(wfn, "CIS") ) {
    return 1;
  }
  else {
    return 0;
  }
}

/*!
** cc_wfn(): Checks if the given wavefunction string is a coupled-cluster
** type and returns 1 if yes and 0 if no.
**
** Note: "coupled-cluster type" means it is handled by PSI like the
** coupled-cluster codes, not necessarily that it is literally a
** coupled-cluster wavefunction
**
** \param wfn = wavefunction string
**
** Returns: 1 if the WFN is a CC method, 0 otherwise
**
** \ingroup QT
*/
int cc_wfn(std::string wfn)
{
  return cc_wfn(wfn.c_str());
}

}
