/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

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

    outfile->Printf( "\n    Number of ijk index combinations:   %14d\n", nijk);

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
      outfile->Printf( "\n    Number of abc index combinations:   %14d\n", nabc);
    } // dertype==1

  }
  else if(params.ref == 2) { /** UHF **/

    aoccpi = moinfo.aoccpi;
    aocc_off = moinfo.aocc_off;
    boccpi = moinfo.boccpi;
    bocc_off = moinfo.bocc_off;

    outfile->Printf( "\n    Number of ijk index combinations:\n");

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

    outfile->Printf( "    Spin Case AAA:                      %14d\n", nijk);

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

    outfile->Printf( "    Spin Case BBB:                      %14d\n", nijk);

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

    outfile->Printf( "    Spin Case AAB:                      %14d\n", nijk);

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

    outfile->Printf( "    Spin Case ABB:                      %14d\n", nijk);

    if(params.dertype == 1) {
      avirtpi = moinfo.avirtpi;
      avir_off = moinfo.avir_off;
      bvirtpi = moinfo.bvirtpi;
      bvir_off = moinfo.bvir_off;

      outfile->Printf( "\n    Number of abc index combinations:\n");
      nabc = 0;
      for (Ga=0; Ga < nirreps; ++Ga)
        for (a=0; a < avirtpi[Ga]; ++a)
          for (Gb=0; Gb < nirreps; ++Gb)
            for (b=0; b < avirtpi[Gb]; ++b)
              for (Gc=0; Gc < nirreps; ++Gc)
                for (c=0; c < avirtpi[Gc]; ++c) nabc++;
      outfile->Printf( "    Spin Case AAA:                      %14d\n", nabc);

      nabc = 0;
      for (Ga=0; Ga < nirreps; ++Ga)
        for (a=0; a < bvirtpi[Ga]; ++a)
          for (Gb=0; Gb < nirreps; ++Gb)
            for (b=0; b < bvirtpi[Gb]; ++b)
              for (Gc=0; Gc < nirreps; ++Gc)
                for (c=0; c < bvirtpi[Gc]; ++c) nabc++;
      outfile->Printf( "    Spin Case BBB:                      %14d\n", nabc);

      nabc = 0;
      for (Ga=0; Ga < nirreps; ++Ga)
        for (a=0; a < avirtpi[Ga]; ++a)
          for (Gb=0; Gb < nirreps; ++Gb)
            for (b=0; b < avirtpi[Gb]; ++b)
              for (Gc=0; Gc < nirreps; ++Gc)
                for (c=0; c < bvirtpi[Gc]; ++c) nabc++;
      outfile->Printf( "    Spin Case AAB:                      %14d\n", nabc);

      nabc = 0;
      for (Ga=0; Ga < nirreps; ++Ga)
        for (a=0; a < avirtpi[Ga]; ++a)
          for (Gb=0; Gb < nirreps; ++Gb)
            for (b=0; b < bvirtpi[Gb]; ++b)
              for (Gc=0; Gc < nirreps; ++Gc)
                for (c=0; c < bvirtpi[Gc]; ++c) nabc++;
      outfile->Printf( "    Spin Case ABB:                      %14d\n", nabc);

    } // dertype == 1

  } // UHF

}

}} // namespace psi::CCTRIPLES
