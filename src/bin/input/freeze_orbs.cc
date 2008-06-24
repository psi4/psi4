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
#include "input.h"
#include "global.h"
#include "defines.h"

namespace psi { namespace input {

void freeze_core()
{
  int large = 0, atom;

  if (frozen_core == "FALSE" || frozen_core == "NO") {
    nfzc = 0;
  }
  else if (frozen_core == "TRUE" || frozen_core == "YES" ||
	frozen_core == "SMALL" || frozen_core == "LARGE") {

    if (frozen_core == "LARGE") large = 1;

    nfzc = 0;
    for(atom=0; atom<num_atoms; atom++) {
      /* H - Be */
      if (nuclear_charges[atom] < 4.1)
        continue;
      /* B - Ne */
      else if (nuclear_charges[atom] > 4.9 && nuclear_charges[atom] < 10.1)
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

  return;
}

}} // namespace psi::input
