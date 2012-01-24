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
  int Ga, Gb, Gc;
  int a, b, c;
  int A, B, C;
  int *occpi, *aoccpi, *boccpi;
  int *virtpi, *avirtpi, *bvirtpi;
  int *occ_off, *aocc_off, *bocc_off;
  int *vir_off, *avir_off, *bvir_off;
  int nijk;
  int nabc;

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

                if(params.dertype == 1) nijk++;
		else if(I >= J && J >= K) nijk++;
	      }
	    }
	  }
	}
      }
    }

    fprintf(outfile, "\n\tNumber of ijk index combinations: %d\n", nijk);

    if(params.dertype == 1) {
      virtpi = moinfo.virtpi;
      vir_off = moinfo.vir_off;
      nabc = 0;
      for (Ga=0; Ga < nirreps; ++Ga) {
        for (a=0; a< virtpi[Ga]; ++a) {
          for (Gb=0; Gb < nirreps; ++Gb) {
            for (b=0; b< virtpi[Gb]; ++b) {
              for (Gc=0; Gc < nirreps; ++Gc) {
                for (c=0; c < virtpi[Gc]; ++c) {
                  nabc++;
                }
              }
            }
          }
        }
      }
      fprintf(outfile, "\n\tNumber of abc index combinations: %d\n", nabc);
    } // dertype==1

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

                if(params.dertype == 1) nijk++;
		else if(I > J && J > K) nijk++;
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

                if(params.dertype == 1) nijk++;
		else if(I > J && J > K) nijk++;
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

                if(params.dertype == 1) nijk++;
		else if(I > J) nijk++;
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

                if(params.dertype == 1) nijk++;
		else if(J > K) nijk++;
	      }
	    }
	  }
	}
      }
    }

    fprintf(outfile, "\tSpin Case ABB: %d\n", nijk);

    if(params.dertype == 1) {
      avirtpi = moinfo.avirtpi;
      avir_off = moinfo.avir_off;
      bvirtpi = moinfo.bvirtpi;
      bvir_off = moinfo.bvir_off;

      fprintf(outfile, "\n\tNumber of abc index combinations:\n");
      nabc = 0;
      for (Ga=0; Ga < nirreps; ++Ga)
        for (a=0; a < avirtpi[Ga]; ++a)
          for (Gb=0; Gb < nirreps; ++Gb)
            for (b=0; b < avirtpi[Gb]; ++b)
              for (Gc=0; Gc < nirreps; ++Gc)
                for (c=0; c < avirtpi[Gc]; ++c) nabc++;
      fprintf(outfile, "\tSpin Case AAA: %d\n", nabc);

      nabc = 0;
      for (Ga=0; Ga < nirreps; ++Ga)
        for (a=0; a < bvirtpi[Ga]; ++a)
          for (Gb=0; Gb < nirreps; ++Gb)
            for (b=0; b < bvirtpi[Gb]; ++b)
              for (Gc=0; Gc < nirreps; ++Gc)
                for (c=0; c < bvirtpi[Gc]; ++c) nabc++;
      fprintf(outfile, "\tSpin Case BBB: %d\n", nabc);

      nabc = 0;
      for (Ga=0; Ga < nirreps; ++Ga)
        for (a=0; a < avirtpi[Ga]; ++a)
          for (Gb=0; Gb < nirreps; ++Gb)
            for (b=0; b < avirtpi[Gb]; ++b)
              for (Gc=0; Gc < nirreps; ++Gc)
                for (c=0; c < bvirtpi[Gc]; ++c) nabc++;
      fprintf(outfile, "\tSpin Case AAB: %d\n", nabc);

      nabc = 0;
      for (Ga=0; Ga < nirreps; ++Ga)
        for (a=0; a < avirtpi[Ga]; ++a)
          for (Gb=0; Gb < nirreps; ++Gb)
            for (b=0; b < bvirtpi[Gb]; ++b)
              for (Gc=0; Gc < nirreps; ++Gc)
                for (c=0; c < bvirtpi[Gc]; ++c) nabc++;
      fprintf(outfile, "\tSpin Case ABB: %d\n", nabc);

    } // dertype == 1

  } // UHF

}

}} // namespace psi::CCTRIPLES
