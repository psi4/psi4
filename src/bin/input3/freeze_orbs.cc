/*! \file
    \ingroup INPUT
    \brief Enter brief description of file here
*/
#define EXTERN
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <libipv1/ip_lib.h>
#include <libciomr/libciomr.h>
#include <libchkpt/chkpt.h>
#include "input.h"
#include "global.h"
#include "defines.h"

namespace psi { namespace input {

void freeze_core()
{
  int large = 0;
  int atom;

  if (!strcmp(frozen_core,"FALSE") ||
      !strcmp(frozen_core,"NO")) {
    nfzc = 0;
    return;
  }
  else if (!strcmp(frozen_core,"TRUE") ||
	   !strcmp(frozen_core,"YES") ||
	   !strcmp(frozen_core,"SMALL") ||
	   !strcmp(frozen_core,"LARGE")) {

      if (!strcmp(frozen_core,"LARGE"))
	large = 1;

      nfzc = 0;
      for(atom=0; atom<num_atoms; atom++) {
	/* H - He */
	if (nuclear_charges[atom] < 2.1)
	  continue;
	/* Li - Ne */
	else if (nuclear_charges[atom] > 2.9 && nuclear_charges[atom] < 10.1)
	  nfzc++;
	/* Na - Ar */
	else if (nuclear_charges[atom] > 10.9 && nuclear_charges[atom] < 18.1)
	  if (large)
	    nfzc += 5;
	  else
	    nfzc++;
	else
	  punt("Cannot freeze core automatically for fourth and higher row elements yet");
      }
  }
  else if (frozen_core[0] >= '0' && frozen_core[0] <= '9') {
    nfzc = atoi(frozen_core);
  }
  else
    punt("Invalid value for FREEZE_CORE");

}


void freeze_virt()
{
  if (!strcmp(frozen_virt,"FALSE") ||
      !strcmp(frozen_virt,"NO")) {
    nfzv = 0;
    return;
  }
  else if (frozen_virt[0] >= '0' && frozen_virt[0] <= '9') {
    nfzv = atoi(frozen_virt);
  }
  else
    punt("Invalid value for FREEZE_VIRT");

}

}} // namespace psi::input
