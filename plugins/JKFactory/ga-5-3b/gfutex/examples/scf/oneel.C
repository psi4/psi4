#include <ga.h>

#include "cscc.h"
#include "scf.h"
#include "oneel.h"

void oneel(double schwmax, double *eone, int nproc)
{
  int lo[2], hi[2], i, j, iloc, jloc, ld;
  double work[ichunk][ichunk],tfock[ichunk][ichunk];
  int dotask;

  //     fill in the one-electron part of the fock matrix and;
  //     compute the one-electron energy contribution;

  GA_Zero(g_counter);
  dotask = next_chunk(lo, hi);
  ld = ichunk;
  while (dotask) {
    NGA_Get(g_schwarz, lo, hi, work, &ld);

    for (i = lo[0];i <= hi[0];i++) {
      iloc = i - lo[0];
      for (j = lo[1]; j <= hi[1]; j++) {
	jloc = j - lo[1];

	tfock[iloc][jloc] = 0.00;
	if ((work[iloc][jloc] * schwmax) > tol2e) {
	  tfock[iloc][jloc] = h(i, j);
	}
      }
    }
    NGA_Put(g_fock, lo, hi, tfock, &ld);

    dotask = next_chunk(lo, hi);
  }
  *eone = 0.50 * contract_matrices(g_fock, g_dens);
} // oneel
