/*! \file
    \ingroup CCTRIPLES
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cctriples {

void count_ijk(void)
{
  int nirreps;
  int Gi, Gj, Gk;
  int i, j, k;
  int I, J, K;
  int *occpi, *aoccpi, *boccpi;
  int *occ_off, *aocc_off, *bocc_off;
  int nijk;

  nirreps = moinfo.nirreps;

  if(params.ref == 0) { /** RHF **/
    occpi = moinfo.occpi;
    occ_off = moinfo.occ_off;

    nijk = 0;
    for(Gi=0; Gi < nirreps; Gi++) {
      for(Gj=0; Gj < nirreps; Gj++) {
	for(Gk=0; Gk < nirreps; Gk++) {
	  for(i=0; i < occpi[Gi]; i++) {
	    I = occ_off[Gi] + i;
	    for(j=0; j < occpi[Gj]; j++) {
	      J = occ_off[Gj] + j;
	      for(k=0; k < occpi[Gk]; k++) {
		K = occ_off[Gk] + k;

		if(I >= J && J >= K) nijk++;
	      }
	    }
	  }
	}
      }
    }

    fprintf(outfile, "\n\tNumber of ijk index combinations: %d\n", nijk);
  }
  else if(params.ref == 2) { /** UHF **/

    aoccpi = moinfo.aoccpi;
    aocc_off = moinfo.aocc_off;
    boccpi = moinfo.boccpi;
    bocc_off = moinfo.bocc_off;

    fprintf(outfile, "\n\tNumber of ijk index combinations:\n");

    /** AAA **/
    nijk = 0;
    for(Gi=0; Gi < nirreps; Gi++) {
      for(Gj=0; Gj < nirreps; Gj++) {
	for(Gk=0; Gk < nirreps; Gk++) {
	  for(i=0; i < aoccpi[Gi]; i++) {
	    I = aocc_off[Gi] + i;
	    for(j=0; j < aoccpi[Gj]; j++) {
	      J = aocc_off[Gj] + j;
	      for(k=0; k < aoccpi[Gk]; k++) {
		K = aocc_off[Gk] + k;

		if(I > J && J > K) nijk++;
	      }
	    }
	  }
	}
      }
    }

    fprintf(outfile, "\tSpin Case AAA: %d\n", nijk);

    /** BBB **/
    nijk = 0;
    for(Gi=0; Gi < nirreps; Gi++) {
      for(Gj=0; Gj < nirreps; Gj++) {
	for(Gk=0; Gk < nirreps; Gk++) {
	  for(i=0; i < boccpi[Gi]; i++) {
	    I = bocc_off[Gi] + i;
	    for(j=0; j < boccpi[Gj]; j++) {
	      J = bocc_off[Gj] + j;
	      for(k=0; k < boccpi[Gk]; k++) {
		K = bocc_off[Gk] + k;

		if(I > J && J > K) nijk++;
	      }
	    }
	  }
	}
      }
    }

    fprintf(outfile, "\tSpin Case BBB: %d\n", nijk);

    /** AAB **/
    nijk = 0;
    for(Gi=0; Gi < nirreps; Gi++) {
      for(Gj=0; Gj < nirreps; Gj++) {
	for(Gk=0; Gk < nirreps; Gk++) {
	  for(i=0; i < aoccpi[Gi]; i++) {
	    I = aocc_off[Gi] + i;
	    for(j=0; j < aoccpi[Gj]; j++) {
	      J = aocc_off[Gj] + j;
	      for(k=0; k < boccpi[Gk]; k++) {
		K = bocc_off[Gk] + k;

		if(I > J) nijk++;
	      }
	    }
	  }
	}
      }
    }

    fprintf(outfile, "\tSpin Case AAB: %d\n", nijk);

    /** ABB **/
    nijk = 0;
    for(Gi=0; Gi < nirreps; Gi++) {
      for(Gj=0; Gj < nirreps; Gj++) {
	for(Gk=0; Gk < nirreps; Gk++) {
	  for(i=0; i < aoccpi[Gi]; i++) {
	    I = aocc_off[Gi] + i;
	    for(j=0; j < boccpi[Gj]; j++) {
	      J = bocc_off[Gj] + j;
	      for(k=0; k < boccpi[Gk]; k++) {
		K = bocc_off[Gk] + k;

		if(J > K) nijk++;
	      }
	    }
	  }
	}
      }
    }

    fprintf(outfile, "\tSpin Case ABB: %d\n", nijk);

  }

  fprintf(outfile, "\n");

}

}} // namespace psi::cctriples
