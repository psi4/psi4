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
    \ingroup CCDENSITY
    \brief Enter brief description of file here
*/
#include <cstdio>
#include "psi4/libdpd/dpd.h"
#include "MOInfo.h"
#include "Params.h"
#include "Frozen.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccdensity {


/*
** NB dpd_list[0] contains only active orbitals and dpd_list[1] contains
** all orbitals.
*/

void resort_gamma(void)
{
  int h, nirreps, row, col, nfzc, nfzv;
  int i,j,k,l,I,J,K,L, ij, kl, ib, ja, ci;
  int a,b,c,d,A,B,C,D, ab, cd;
  int ka;
  int *qt_occ, *qt_vir;
  int *cc_occ, *cc_vir;
  dpdfile2 g, gnew;
  dpdbuf4 G, Gnew;

  nirreps = moinfo.nirreps;
  qt_occ = frozen.qt_occ;  qt_vir = frozen.qt_vir;
  cc_occ = frozen.cc_occ;  cc_vir = frozen.cc_vir;
  nfzc = moinfo.nfzc;  nfzv = moinfo.nfzv;

  /* G(I,J) --> G'(I,J) */
  dpd_set_default(0);
  global_dpd_->file2_init(&g, PSIF_CC_OEI, 0, 0, 0, "DIJ");

  dpd_set_default(1);
  global_dpd_->file2_init(&gnew, PSIF_CC_OEI_NEW, 0, 0, 0, "DIJ");

  global_dpd_->file2_mat_init(&g);
  global_dpd_->file2_mat_rd(&g);
  global_dpd_->file2_mat_init(&gnew);

  for(h=0; h < nirreps; h++) {
      for(row=0; row < gnew.params->rowtot[h]; row++) {
	  i = gnew.params->roworb[h][row];

	  I = qt_occ[i]; i = cc_occ[I];

	  for(col=0; col < gnew.params->coltot[h]; col++) {
	      j = gnew.params->colorb[h][col];

	      J = qt_occ[j]; j = cc_occ[J];

	      if(i<0 || j<0)
		  gnew.matrix[h][row][col] = 0.0;
	      else {
		  i = g.params->rowidx[i];
		  j = g.params->colidx[j];

		  gnew.matrix[h][row][col] = g.matrix[h][i][j];
		}
	    }
	}
    }

  global_dpd_->file2_mat_wrt(&gnew);
  global_dpd_->file2_mat_close(&gnew);
  global_dpd_->file2_mat_close(&g);

  global_dpd_->file2_close(&g);
  global_dpd_->file2_close(&gnew);

  /* G(i,j) --> G'(i,j) */
  dpd_set_default(0);
  global_dpd_->file2_init(&g, PSIF_CC_OEI, 0, 0, 0, "Dij");

  dpd_set_default(1);
  global_dpd_->file2_init(&gnew, PSIF_CC_OEI_NEW, 0, 0, 0, "Dij");

  global_dpd_->file2_mat_init(&g);
  global_dpd_->file2_mat_rd(&g);
  global_dpd_->file2_mat_init(&gnew);

  for(h=0; h < nirreps; h++) {
      for(row=0; row < gnew.params->rowtot[h]; row++) {
	  i = gnew.params->roworb[h][row];

	  I = qt_occ[i]; i = cc_occ[I];

	  for(col=0; col < gnew.params->coltot[h]; col++) {
	      j = gnew.params->colorb[h][col];

	      J = qt_occ[j]; j = cc_occ[J];

	      if(i<0 || j<0)
		  gnew.matrix[h][row][col] = 0.0;
	      else {
		  i = g.params->rowidx[i];
		  j = g.params->colidx[j];

		  gnew.matrix[h][row][col] = g.matrix[h][i][j];
		}
	    }
	}
    }

  global_dpd_->file2_mat_wrt(&gnew);
  global_dpd_->file2_mat_close(&gnew);
  global_dpd_->file2_mat_close(&g);

  global_dpd_->file2_close(&g);
  global_dpd_->file2_close(&gnew);

  /* G(A,B) --> G'(A,B) */
  dpd_set_default(0);
  global_dpd_->file2_init(&g, PSIF_CC_OEI, 0, 1, 1, "DAB");

  dpd_set_default(1);
  global_dpd_->file2_init(&gnew, PSIF_CC_OEI_NEW, 0, 1, 1, "DAB");

  global_dpd_->file2_mat_init(&g);
  global_dpd_->file2_mat_rd(&g);
  global_dpd_->file2_mat_init(&gnew);

  for(h=0; h < nirreps; h++) {
      for(row=0; row < gnew.params->rowtot[h]; row++) {
	  a = gnew.params->roworb[h][row];

	  A = qt_vir[a]; a = cc_vir[A];

	  for(col=0; col < gnew.params->coltot[h]; col++) {
	      b = gnew.params->colorb[h][col];

	      B = qt_vir[b]; b = cc_vir[B];

	      if(a<0 || b<0)
		  gnew.matrix[h][row][col] = 0.0;
	      else {
		  a = g.params->rowidx[a];
		  b = g.params->colidx[b];

		  gnew.matrix[h][row][col] = g.matrix[h][a][b];
		}
	    }
	}
    }

  global_dpd_->file2_mat_wrt(&gnew);
  global_dpd_->file2_mat_close(&gnew);
  global_dpd_->file2_mat_close(&g);

  global_dpd_->file2_close(&g);
  global_dpd_->file2_close(&gnew);

  /* G(a,b) --> G'(a,b) */
  dpd_set_default(0);
  global_dpd_->file2_init(&g, PSIF_CC_OEI, 0, 1, 1, "Dab");

  dpd_set_default(1);
  global_dpd_->file2_init(&gnew, PSIF_CC_OEI_NEW, 0, 1, 1, "Dab");

  global_dpd_->file2_mat_init(&g);
  global_dpd_->file2_mat_rd(&g);
  global_dpd_->file2_mat_init(&gnew);

  for(h=0; h < nirreps; h++) {
      for(row=0; row < gnew.params->rowtot[h]; row++) {
	  a = gnew.params->roworb[h][row];

	  A = qt_vir[a]; a = cc_vir[A];

	  for(col=0; col < gnew.params->coltot[h]; col++) {
	      b = gnew.params->colorb[h][col];

	      B = qt_vir[b]; b = cc_vir[B];

	      if(a<0 || b<0)
		  gnew.matrix[h][row][col] = 0.0;
	      else {
		  a = g.params->rowidx[a];
		  b = g.params->colidx[b];

		  gnew.matrix[h][row][col] = g.matrix[h][a][b];
		}
	    }
	}
    }

  global_dpd_->file2_mat_wrt(&gnew);
  global_dpd_->file2_mat_close(&gnew);
  global_dpd_->file2_mat_close(&g);

  global_dpd_->file2_close(&g);
  global_dpd_->file2_close(&gnew);

  /* G(A,I) --> G'(A,I) */
  dpd_set_default(0);
  global_dpd_->file2_init(&g, PSIF_CC_OEI, 0, 0, 1, "DAI");

  dpd_set_default(1);
  global_dpd_->file2_init(&gnew, PSIF_CC_OEI_NEW, 0, 0, 1, "DAI");

  global_dpd_->file2_mat_init(&g);
  global_dpd_->file2_mat_rd(&g);
  global_dpd_->file2_mat_init(&gnew);

  for(h=0; h < nirreps; h++) {
      for(row=0; row < gnew.params->rowtot[h]; row++) {
	  i = gnew.params->roworb[h][row];

	  I = qt_occ[i]; i = cc_occ[I];

	  for(col=0; col < gnew.params->coltot[h]; col++) {
	      a = gnew.params->colorb[h][col];

	      A = qt_vir[a]; a = cc_vir[A];


	      if(i<0 || a<0)
		  gnew.matrix[h][row][col] = 0.0;
	      else {
		  i = g.params->rowidx[i];
		  a = g.params->colidx[a];

		  gnew.matrix[h][row][col] = g.matrix[h][i][a];
		}
	    }
	}
    }

  global_dpd_->file2_mat_wrt(&gnew);
  global_dpd_->file2_mat_close(&gnew);
  global_dpd_->file2_mat_close(&g);

  global_dpd_->file2_close(&g);
  global_dpd_->file2_close(&gnew);

  /* G(a,i) --> G'(a,i) */
  dpd_set_default(0);
  global_dpd_->file2_init(&g, PSIF_CC_OEI, 0, 0, 1, "Dai");

  dpd_set_default(1);
  global_dpd_->file2_init(&gnew, PSIF_CC_OEI_NEW, 0, 0, 1, "Dai");

  global_dpd_->file2_mat_init(&g);
  global_dpd_->file2_mat_rd(&g);
  global_dpd_->file2_mat_init(&gnew);

  for(h=0; h < nirreps; h++) {
      for(row=0; row < gnew.params->rowtot[h]; row++) {
	  i = gnew.params->roworb[h][row];

	  I = qt_occ[i]; i = cc_occ[I];

	  for(col=0; col < gnew.params->coltot[h]; col++) {
	      a = gnew.params->colorb[h][col];

	      A = qt_vir[a]; a = cc_vir[A];


	      if(i<0 || a<0)
		  gnew.matrix[h][row][col] = 0.0;
	      else {
		  i = g.params->rowidx[i];
		  a = g.params->colidx[a];

		  gnew.matrix[h][row][col] = g.matrix[h][i][a];
		}
	    }
	}
    }

  global_dpd_->file2_mat_wrt(&gnew);
  global_dpd_->file2_mat_close(&gnew);
  global_dpd_->file2_mat_close(&g);

  global_dpd_->file2_close(&g);
  global_dpd_->file2_close(&gnew);

  /* G(I,A) --> G'(I,A) */
  dpd_set_default(0);
  global_dpd_->file2_init(&g, PSIF_CC_OEI, 0, 0, 1, "DIA");

  dpd_set_default(1);
  global_dpd_->file2_init(&gnew, PSIF_CC_OEI_NEW, 0, 0, 1, "DIA");

  global_dpd_->file2_mat_init(&g);
  global_dpd_->file2_mat_rd(&g);
  global_dpd_->file2_mat_init(&gnew);

  for(h=0; h < nirreps; h++) {
      for(row=0; row < gnew.params->rowtot[h]; row++) {
	  i = gnew.params->roworb[h][row];

	  I = qt_occ[i]; i = cc_occ[I];

	  for(col=0; col < gnew.params->coltot[h]; col++) {
	      a = gnew.params->colorb[h][col];

	      A = qt_vir[a]; a = cc_vir[A];

	      if(i<0 || a<0)
		  gnew.matrix[h][row][col] = 0.0;
	      else {
		  i = g.params->rowidx[i];
		  a = g.params->colidx[a];

		  gnew.matrix[h][row][col] = g.matrix[h][i][a];
		}
	    }
	}
    }

  global_dpd_->file2_mat_wrt(&gnew);
  global_dpd_->file2_mat_close(&gnew);
  global_dpd_->file2_mat_close(&g);

  global_dpd_->file2_close(&g);
  global_dpd_->file2_close(&gnew);

  /* G(i,a) --> G'(i,a) */
  dpd_set_default(0);
  global_dpd_->file2_init(&g, PSIF_CC_OEI, 0, 0, 1, "Dia");

  dpd_set_default(1);
  global_dpd_->file2_init(&gnew, PSIF_CC_OEI_NEW, 0, 0, 1, "Dia");

  global_dpd_->file2_mat_init(&g);
  global_dpd_->file2_mat_rd(&g);
  global_dpd_->file2_mat_init(&gnew);

  for(h=0; h < nirreps; h++) {
      for(row=0; row < gnew.params->rowtot[h]; row++) {
	  i = gnew.params->roworb[h][row];

	  I = qt_occ[i]; i = cc_occ[I];

	  for(col=0; col < gnew.params->coltot[h]; col++) {
	      a = gnew.params->colorb[h][col];

	      A = qt_vir[a]; a = cc_vir[A];

	      if(i<0 || a<0)
		  gnew.matrix[h][row][col] = 0.0;
	      else {
		  i = g.params->rowidx[i];
		  a = g.params->colidx[a];

		  gnew.matrix[h][row][col] = g.matrix[h][i][a];
		}
	    }
	}
    }

  global_dpd_->file2_mat_wrt(&gnew);
  global_dpd_->file2_mat_close(&gnew);
  global_dpd_->file2_mat_close(&g);

  global_dpd_->file2_close(&g);
  global_dpd_->file2_close(&gnew);


  /* G(IJ,KL) --> G'(IJ,KL) */
  dpd_set_default(0);
  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 2, 2, 2, 2, 0, "GIJKL");

  dpd_set_default(1);
  global_dpd_->buf4_init(&Gnew, PSIF_CC_GAMMA_NEW, 0, 2, 2, 2, 2, 0, "GIJKL");

  for(h=0; h < nirreps; h++) {
      global_dpd_->buf4_mat_irrep_init(&G, h);
      global_dpd_->buf4_mat_irrep_rd(&G, h);

      global_dpd_->buf4_mat_irrep_init(&Gnew, h);

      for(row=0; row < Gnew.params->rowtot[h]; row++) {
	  i = Gnew.params->roworb[h][row][0];
	  j = Gnew.params->roworb[h][row][1];

	  /* Compute the active CC row index */
	  I = qt_occ[i];  J = qt_occ[j];
	  i = cc_occ[I];  j = cc_occ[J];

	  for(col=0; col < Gnew.params->coltot[h]; col++) {
	      k = Gnew.params->colorb[h][col][0];
	      l = Gnew.params->colorb[h][col][1];

	      /* Compute the active CC column index */
	      K = qt_occ[k];  L = qt_occ[l];
	      k = cc_occ[K];  l = cc_occ[L];


	      if(i<0 || j<0 || k<0 || l<0)
		  Gnew.matrix[h][row][col] = 0.0;
	      else {
		  ij = G.params->rowidx[i][j];
		  kl = G.params->colidx[k][l];

		  Gnew.matrix[h][row][col] = G.matrix[h][ij][kl];
		}
	    }
	}

      global_dpd_->buf4_mat_irrep_wrt(&Gnew, h);

      global_dpd_->buf4_mat_irrep_close(&G, h);
    }

  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_close(&Gnew);

  /* G(ij,kl) --> G'(ij,kl) */
  dpd_set_default(0);
  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 2, 2, 2, 2, 0, "Gijkl");

  dpd_set_default(1);
  global_dpd_->buf4_init(&Gnew, PSIF_CC_GAMMA_NEW, 0, 2, 2, 2, 2, 0, "Gijkl");

  for(h=0; h < nirreps; h++) {
      global_dpd_->buf4_mat_irrep_init(&G, h);
      global_dpd_->buf4_mat_irrep_rd(&G, h);

      global_dpd_->buf4_mat_irrep_init(&Gnew, h);

      for(row=0; row < Gnew.params->rowtot[h]; row++) {
	  i = Gnew.params->roworb[h][row][0];
	  j = Gnew.params->roworb[h][row][1];

	  /* Compute the active CC row index */
	  I = qt_occ[i];  J = qt_occ[j];
	  i = cc_occ[I];  j = cc_occ[J];

	  for(col=0; col < Gnew.params->coltot[h]; col++) {
	      k = Gnew.params->colorb[h][col][0];
	      l = Gnew.params->colorb[h][col][1];

	      /* Compute the active CC column index */
	      K = qt_occ[k];  L = qt_occ[l];
	      k = cc_occ[K];  l = cc_occ[L];

	      if(i<0 || j<0 || k<0 || l<0)
		  Gnew.matrix[h][row][col] = 0.0;
	      else {
		  ij = G.params->rowidx[i][j];
		  kl = G.params->colidx[k][l];

		  Gnew.matrix[h][row][col] = G.matrix[h][ij][kl];
		}
	    }
	}

      global_dpd_->buf4_mat_irrep_wrt(&Gnew, h);

      global_dpd_->buf4_mat_irrep_close(&G, h);
    }

  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_close(&Gnew);

  /* G(Ij,Kl) --> G'(Ij,Kl) */
  dpd_set_default(0);
  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 0, 0, 0, 0, 0, "GIjKl");

  dpd_set_default(1);
  global_dpd_->buf4_init(&Gnew, PSIF_CC_GAMMA_NEW, 0, 0, 0, 0, 0, 0, "GIjKl");

  for(h=0; h < nirreps; h++) {
      global_dpd_->buf4_mat_irrep_init(&G, h);
      global_dpd_->buf4_mat_irrep_rd(&G, h);

      global_dpd_->buf4_mat_irrep_init(&Gnew, h);

      for(row=0; row < Gnew.params->rowtot[h]; row++) {
	  i = Gnew.params->roworb[h][row][0];
	  j = Gnew.params->roworb[h][row][1];

	  /* Compute the active CC row index */
	  I = qt_occ[i];  J = qt_occ[j];
	  i = cc_occ[I];  j = cc_occ[J];


	  for(col=0; col < Gnew.params->coltot[h]; col++) {
	      k = Gnew.params->colorb[h][col][0];
	      l = Gnew.params->colorb[h][col][1];

	      /* Compute the active CC column index */
	      K = qt_occ[k];  L = qt_occ[l];
	      k = cc_occ[K];  l = cc_occ[L];

	      if(i<0 || j<0 || k<0 || l<0)
		  Gnew.matrix[h][row][col] = 0.0;
	      else {
		  ij = G.params->rowidx[i][j];
		  kl = G.params->colidx[k][l];

		  Gnew.matrix[h][row][col] = G.matrix[h][ij][kl];
		}
	    }
	}

      global_dpd_->buf4_mat_irrep_wrt(&Gnew, h);

      global_dpd_->buf4_mat_irrep_close(&G, h);
    }

  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_close(&Gnew);

  /* G(IJ,KA) --> G'(IJ,KA) */
  dpd_set_default(0);
  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 2, 10, 2, 10, 0, "GIJKA");

  dpd_set_default(1);
  global_dpd_->buf4_init(&Gnew, PSIF_CC_GAMMA_NEW, 0, 2, 10, 2, 10, 0, "GIJKA");

  for(h=0; h < nirreps; h++) {
      global_dpd_->buf4_mat_irrep_init(&G, h);
      global_dpd_->buf4_mat_irrep_rd(&G, h);

      global_dpd_->buf4_mat_irrep_init(&Gnew, h);

      for(row=0; row < Gnew.params->rowtot[h]; row++) {
	  i = Gnew.params->roworb[h][row][0];
	  j = Gnew.params->roworb[h][row][1];

	  /* Compute the active CC row index */
	  I = qt_occ[i];  J = qt_occ[j];
	  i = cc_occ[I];  j = cc_occ[J];


	  for(col=0; col < Gnew.params->coltot[h]; col++) {
	      k = Gnew.params->colorb[h][col][0];
	      a = Gnew.params->colorb[h][col][1];

	      /* Compute the active CC column index */
	      K = qt_occ[k];  A = qt_vir[a];
	      k = cc_occ[K];  a = cc_vir[A];

	      if(i<0 || j<0 || k<0 || a<0)
		  Gnew.matrix[h][row][col] = 0.0;
	      else {
		  ij = G.params->rowidx[i][j];
		  ka = G.params->colidx[k][a];

		  Gnew.matrix[h][row][col] = G.matrix[h][ij][ka];
		}
	    }
	}

      global_dpd_->buf4_mat_irrep_wrt(&Gnew, h);

      global_dpd_->buf4_mat_irrep_close(&G, h);
    }

  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_close(&Gnew);

  /* G(ij,ka) --> G'(ij,ka) */
  dpd_set_default(0);
  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 2, 10, 2, 10, 0, "Gijka");

  dpd_set_default(1);
  global_dpd_->buf4_init(&Gnew, PSIF_CC_GAMMA_NEW, 0, 2, 10, 2, 10, 0, "Gijka");

  for(h=0; h < nirreps; h++) {
      global_dpd_->buf4_mat_irrep_init(&G, h);
      global_dpd_->buf4_mat_irrep_rd(&G, h);

      global_dpd_->buf4_mat_irrep_init(&Gnew, h);

      for(row=0; row < Gnew.params->rowtot[h]; row++) {
	  i = Gnew.params->roworb[h][row][0];
	  j = Gnew.params->roworb[h][row][1];

	  /* Compute the active CC row index */
	  I = qt_occ[i];  J = qt_occ[j];
	  i = cc_occ[I];  j = cc_occ[J];


	  for(col=0; col < Gnew.params->coltot[h]; col++) {
	      k = Gnew.params->colorb[h][col][0];
	      a = Gnew.params->colorb[h][col][1];

	      /* Compute the active CC column index */
	      K = qt_occ[k];  A = qt_vir[a];
	      k = cc_occ[K];  a = cc_vir[A];

	      if(i<0 || j<0 || k<0 || a<0)
		  Gnew.matrix[h][row][col] = 0.0;
	      else {
		  ij = G.params->rowidx[i][j];
		  ka = G.params->colidx[k][a];

		  Gnew.matrix[h][row][col] = G.matrix[h][ij][ka];
		}
	    }
	}

      global_dpd_->buf4_mat_irrep_wrt(&Gnew, h);

      global_dpd_->buf4_mat_irrep_close(&G, h);
    }

  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_close(&Gnew);

  /* G(Ij,Ka) --> G'(Ij,Ka) */
  dpd_set_default(0);
  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 0, 10, 0, 10, 0, "GIjKa");

  dpd_set_default(1);
  global_dpd_->buf4_init(&Gnew, PSIF_CC_GAMMA_NEW, 0, 0, 10, 0, 10, 0, "GIjKa");

  for(h=0; h < nirreps; h++) {
      global_dpd_->buf4_mat_irrep_init(&G, h);
      global_dpd_->buf4_mat_irrep_rd(&G, h);

      global_dpd_->buf4_mat_irrep_init(&Gnew, h);

      for(row=0; row < Gnew.params->rowtot[h]; row++) {
	  i = Gnew.params->roworb[h][row][0];
	  j = Gnew.params->roworb[h][row][1];

	  /* Compute the active CC row index */
	  I = qt_occ[i];  J = qt_occ[j];
	  i = cc_occ[I];  j = cc_occ[J];


	  for(col=0; col < Gnew.params->coltot[h]; col++) {
	      k = Gnew.params->colorb[h][col][0];
	      a = Gnew.params->colorb[h][col][1];

	      /* Compute the active CC column index */
	      K = qt_occ[k];  A = qt_vir[a];
	      k = cc_occ[K];  a = cc_vir[A];

	      if(i<0 || j<0 || k<0 || a<0)
		  Gnew.matrix[h][row][col] = 0.0;
	      else {
		  ij = G.params->rowidx[i][j];
		  ka = G.params->colidx[k][a];

		  Gnew.matrix[h][row][col] = G.matrix[h][ij][ka];
		}
	    }
	}

      global_dpd_->buf4_mat_irrep_wrt(&Gnew, h);

      global_dpd_->buf4_mat_irrep_close(&G, h);
    }

  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_close(&Gnew);

  /* G(IJ,AB) --> G'(IJ,AB) */
  dpd_set_default(0);
  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 2, 7, 2, 7, 0, "GIJAB");

  dpd_set_default(1);
  global_dpd_->buf4_init(&Gnew, PSIF_CC_GAMMA_NEW, 0, 2, 7, 2, 7, 0, "GIJAB");

  for(h=0; h < nirreps; h++) {
      global_dpd_->buf4_mat_irrep_init(&G, h);
      global_dpd_->buf4_mat_irrep_rd(&G, h);

      global_dpd_->buf4_mat_irrep_init(&Gnew, h);

      for(row=0; row < Gnew.params->rowtot[h]; row++) {
	  i = Gnew.params->roworb[h][row][0];
	  j = Gnew.params->roworb[h][row][1];

	  /* Compute the active CC row index */
	  I = qt_occ[i];  J = qt_occ[j];
	  i = cc_occ[I];  j = cc_occ[J];

	  for(col=0; col < Gnew.params->coltot[h]; col++) {
	      a = Gnew.params->colorb[h][col][0];
	      b = Gnew.params->colorb[h][col][1];

	      /* Compute the active CC column index */
	      A = qt_vir[a];  B = qt_vir[b];
	      a = cc_vir[A];  b = cc_vir[B];

	      if(i<0 || j<0 || a<0 || b<0)
		  Gnew.matrix[h][row][col] = 0.0;
	      else {
		  ij = G.params->rowidx[i][j];
		  ab = G.params->colidx[a][b];

		  Gnew.matrix[h][row][col] = G.matrix[h][ij][ab];
		}
	    }
	}

      global_dpd_->buf4_mat_irrep_wrt(&Gnew, h);

      global_dpd_->buf4_mat_irrep_close(&G, h);
    }

  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_close(&Gnew);

  /* G(ij,ab) --> G'(ij,ab) */
  dpd_set_default(0);
  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 2, 7, 2, 7, 0, "Gijab");

  dpd_set_default(1);
  global_dpd_->buf4_init(&Gnew, PSIF_CC_GAMMA_NEW, 0, 2, 7, 2, 7, 0, "Gijab");

  for(h=0; h < nirreps; h++) {
      global_dpd_->buf4_mat_irrep_init(&G, h);
      global_dpd_->buf4_mat_irrep_rd(&G, h);

      global_dpd_->buf4_mat_irrep_init(&Gnew, h);

      for(row=0; row < Gnew.params->rowtot[h]; row++) {
	  i = Gnew.params->roworb[h][row][0];
	  j = Gnew.params->roworb[h][row][1];

	  /* Compute the active CC row index */
	  I = qt_occ[i];  J = qt_occ[j];
	  i = cc_occ[I];  j = cc_occ[J];

	  for(col=0; col < Gnew.params->coltot[h]; col++) {
	      a = Gnew.params->colorb[h][col][0];
	      b = Gnew.params->colorb[h][col][1];

	      /* Compute the active CC column index */
	      A = qt_vir[a];  B = qt_vir[b];
	      a = cc_vir[A];  b = cc_vir[B];

	      if(i<0 || j<0 || a<0 || b<0)
		  Gnew.matrix[h][row][col] = 0.0;
	      else {
		  ij = G.params->rowidx[i][j];
		  ab = G.params->colidx[a][b];

		  Gnew.matrix[h][row][col] = G.matrix[h][ij][ab];
		}
	    }
	}

      global_dpd_->buf4_mat_irrep_wrt(&Gnew, h);

      global_dpd_->buf4_mat_irrep_close(&G, h);
    }

  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_close(&Gnew);

  /* G(Ij,Ab) --> G'(Ij,Ab) */
  dpd_set_default(0);
  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 0, 5, 0, 5, 0, "GIjAb");

  dpd_set_default(1);
  global_dpd_->buf4_init(&Gnew, PSIF_CC_GAMMA_NEW, 0, 0, 5, 0, 5, 0, "GIjAb");

  for(h=0; h < nirreps; h++) {
      global_dpd_->buf4_mat_irrep_init(&G, h);
      global_dpd_->buf4_mat_irrep_rd(&G, h);

      global_dpd_->buf4_mat_irrep_init(&Gnew, h);

      for(row=0; row < Gnew.params->rowtot[h]; row++) {
	  i = Gnew.params->roworb[h][row][0];
	  j = Gnew.params->roworb[h][row][1];

	  /* Compute the active CC row index */
	  I = qt_occ[i];  J = qt_occ[j];
	  i = cc_occ[I];  j = cc_occ[J];

	  for(col=0; col < Gnew.params->coltot[h]; col++) {
	      a = Gnew.params->colorb[h][col][0];
	      b = Gnew.params->colorb[h][col][1];

	      /* Compute the active CC column index */
	      A = qt_vir[a];  B = qt_vir[b];
	      a = cc_vir[A];  b = cc_vir[B];

	      if(i<0 || j<0 || a<0 || b<0)
		  Gnew.matrix[h][row][col] = 0.0;
	      else {
		  ij = G.params->rowidx[i][j];
		  ab = G.params->colidx[a][b];

		  Gnew.matrix[h][row][col] = G.matrix[h][ij][ab];
		}
	    }
	}

      global_dpd_->buf4_mat_irrep_wrt(&Gnew, h);

      global_dpd_->buf4_mat_irrep_close(&G, h);
    }

  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_close(&Gnew);

  /* G(IB,JA) --> G'(IB,JA) */
  dpd_set_default(0);
  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 10, 10, 10, 10, 0, "GIBJA");

  dpd_set_default(1);
  global_dpd_->buf4_init(&Gnew, PSIF_CC_GAMMA_NEW, 0, 10, 10, 10, 10, 0, "GIBJA");

  for(h=0; h < nirreps; h++) {
      global_dpd_->buf4_mat_irrep_init(&G, h);
      global_dpd_->buf4_mat_irrep_rd(&G, h);

      global_dpd_->buf4_mat_irrep_init(&Gnew, h);

      for(row=0; row < Gnew.params->rowtot[h]; row++) {
	  i = Gnew.params->roworb[h][row][0];
	  b = Gnew.params->roworb[h][row][1];

	  /* Compute the active CC row index */
	  I = qt_occ[i];  B = qt_vir[b];
	  i = cc_occ[I];  b = cc_vir[B];

	  for(col=0; col < Gnew.params->coltot[h]; col++) {
	      j = Gnew.params->colorb[h][col][0];
	      a = Gnew.params->colorb[h][col][1];

	      /* Compute the active CC column index */
	      J = qt_occ[j];  A = qt_vir[a];
	      j = cc_occ[J];  a = cc_vir[A];

	      if(i<0 || j<0 || a<0 || b<0)
		  Gnew.matrix[h][row][col] = 0.0;
	      else {
		  ib = G.params->rowidx[i][b];
		  ja = G.params->colidx[j][a];

		  Gnew.matrix[h][row][col] = G.matrix[h][ib][ja];
		}
	    }
	}

      global_dpd_->buf4_mat_irrep_wrt(&Gnew, h);

      global_dpd_->buf4_mat_irrep_close(&G, h);
    }

  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_close(&Gnew);

  /* G(ib,ja) --> G'(ib,ja) */
  dpd_set_default(0);
  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 10, 10, 10, 10, 0, "Gibja");

  dpd_set_default(1);
  global_dpd_->buf4_init(&Gnew, PSIF_CC_GAMMA_NEW, 0, 10, 10, 10, 10, 0, "Gibja");

  for(h=0; h < nirreps; h++) {
      global_dpd_->buf4_mat_irrep_init(&G, h);
      global_dpd_->buf4_mat_irrep_rd(&G, h);

      global_dpd_->buf4_mat_irrep_init(&Gnew, h);

      for(row=0; row < Gnew.params->rowtot[h]; row++) {
	  i = Gnew.params->roworb[h][row][0];
	  b = Gnew.params->roworb[h][row][1];

	  /* Compute the active CC row index */
	  I = qt_occ[i];  B = qt_vir[b];
	  i = cc_occ[I];  b = cc_vir[B];

	  for(col=0; col < Gnew.params->coltot[h]; col++) {
	      j = Gnew.params->colorb[h][col][0];
	      a = Gnew.params->colorb[h][col][1];

	      /* Compute the active CC column index */
	      J = qt_occ[j];  A = qt_vir[a];
	      j = cc_occ[J];  a = cc_vir[A];

	      if(i<0 || j<0 || a<0 || b<0)
		  Gnew.matrix[h][row][col] = 0.0;
	      else {
		  ib = G.params->rowidx[i][b];
		  ja = G.params->colidx[j][a];

		  Gnew.matrix[h][row][col] = G.matrix[h][ib][ja];
		}
	    }
	}

      global_dpd_->buf4_mat_irrep_wrt(&Gnew, h);

      global_dpd_->buf4_mat_irrep_close(&G, h);
    }

  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_close(&Gnew);

  /* G(Ib,Ja) --> G'(Ib,Ja) */
  dpd_set_default(0);
  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 10, 10, 10, 10, 0, "GIbJa");

  dpd_set_default(1);
  global_dpd_->buf4_init(&Gnew, PSIF_CC_GAMMA_NEW, 0, 10, 10, 10, 10, 0, "GIbJa");

  for(h=0; h < nirreps; h++) {
      global_dpd_->buf4_mat_irrep_init(&G, h);
      global_dpd_->buf4_mat_irrep_rd(&G, h);

      global_dpd_->buf4_mat_irrep_init(&Gnew, h);

      for(row=0; row < Gnew.params->rowtot[h]; row++) {
	  i = Gnew.params->roworb[h][row][0];
	  b = Gnew.params->roworb[h][row][1];

	  /* Compute the active CC row index */
	  I = qt_occ[i];  B = qt_vir[b];
	  i = cc_occ[I];  b = cc_vir[B];

	  for(col=0; col < Gnew.params->coltot[h]; col++) {
	      j = Gnew.params->colorb[h][col][0];
	      a = Gnew.params->colorb[h][col][1];

	      /* Compute the active CC column index */
	      J = qt_occ[j];  A = qt_vir[a];
	      j = cc_occ[J];  a = cc_vir[A];

	      if(i<0 || j<0 || a<0 || b<0)
		  Gnew.matrix[h][row][col] = 0.0;
	      else {
		  ib = G.params->rowidx[i][b];
		  ja = G.params->colidx[j][a];

		  Gnew.matrix[h][row][col] = G.matrix[h][ib][ja];
		}
	    }
	}

      global_dpd_->buf4_mat_irrep_wrt(&Gnew, h);

      global_dpd_->buf4_mat_irrep_close(&G, h);
    }

  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_close(&Gnew);

  /* G(iB,jA) --> G'(iB,jA) */
  dpd_set_default(0);
  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 10, 10, 10, 10, 0, "GiBjA");

  dpd_set_default(1);
  global_dpd_->buf4_init(&Gnew, PSIF_CC_GAMMA_NEW, 0, 10, 10, 10, 10, 0, "GiBjA");

  for(h=0; h < nirreps; h++) {
      global_dpd_->buf4_mat_irrep_init(&G, h);
      global_dpd_->buf4_mat_irrep_rd(&G, h);

      global_dpd_->buf4_mat_irrep_init(&Gnew, h);

      for(row=0; row < Gnew.params->rowtot[h]; row++) {
	  i = Gnew.params->roworb[h][row][0];
	  b = Gnew.params->roworb[h][row][1];

	  /* Compute the active CC row index */
	  I = qt_occ[i];  B = qt_vir[b];
	  i = cc_occ[I];  b = cc_vir[B];

	  for(col=0; col < Gnew.params->coltot[h]; col++) {
	      j = Gnew.params->colorb[h][col][0];
	      a = Gnew.params->colorb[h][col][1];

	      /* Compute the active CC column index */
	      J = qt_occ[j];  A = qt_vir[a];
	      j = cc_occ[J];  a = cc_vir[A];

	      if(i<0 || j<0 || a<0 || b<0)
		  Gnew.matrix[h][row][col] = 0.0;
	      else {
		  ib = G.params->rowidx[i][b];
		  ja = G.params->colidx[j][a];

		  Gnew.matrix[h][row][col] = G.matrix[h][ib][ja];
		}
	    }
	}

      global_dpd_->buf4_mat_irrep_wrt(&Gnew, h);

      global_dpd_->buf4_mat_irrep_close(&G, h);
    }

  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_close(&Gnew);

  /* G(iB,Ja) --> G'(iB,Ja) */
  dpd_set_default(0);
  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 10, 10, 10, 10, 0, "GiBJa");

  dpd_set_default(1);
  global_dpd_->buf4_init(&Gnew, PSIF_CC_GAMMA_NEW, 0, 10, 10, 10, 10, 0, "GiBJa");

  for(h=0; h < nirreps; h++) {
      global_dpd_->buf4_mat_irrep_init(&G, h);
      global_dpd_->buf4_mat_irrep_rd(&G, h);

      global_dpd_->buf4_mat_irrep_init(&Gnew, h);

      for(row=0; row < Gnew.params->rowtot[h]; row++) {
	  i = Gnew.params->roworb[h][row][0];
	  b = Gnew.params->roworb[h][row][1];

	  /* Compute the active CC row index */
	  I = qt_occ[i];  B = qt_vir[b];
	  i = cc_occ[I];  b = cc_vir[B];

	  for(col=0; col < Gnew.params->coltot[h]; col++) {
	      j = Gnew.params->colorb[h][col][0];
	      a = Gnew.params->colorb[h][col][1];

	      /* Compute the active CC column index */
	      J = qt_occ[j];  A = qt_vir[a];
	      j = cc_occ[J];  a = cc_vir[A];

	      if(i<0 || j<0 || a<0 || b<0)
		  Gnew.matrix[h][row][col] = 0.0;
	      else {
		  ib = G.params->rowidx[i][b];
		  ja = G.params->colidx[j][a];

		  Gnew.matrix[h][row][col] = G.matrix[h][ib][ja];
		}
	    }
	}

      global_dpd_->buf4_mat_irrep_wrt(&Gnew, h);

      global_dpd_->buf4_mat_irrep_close(&G, h);
    }

  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_close(&Gnew);

  /* G(Ib,jA) --> G'(Ib,jA) */
  dpd_set_default(0);
  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 10, 10, 10, 10, 0, "GIbjA");

  dpd_set_default(1);
  global_dpd_->buf4_init(&Gnew, PSIF_CC_GAMMA_NEW, 0, 10, 10, 10, 10, 0, "GIbjA");

  for(h=0; h < nirreps; h++) {
      global_dpd_->buf4_mat_irrep_init(&G, h);
      global_dpd_->buf4_mat_irrep_rd(&G, h);

      global_dpd_->buf4_mat_irrep_init(&Gnew, h);

      for(row=0; row < Gnew.params->rowtot[h]; row++) {
	  i = Gnew.params->roworb[h][row][0];
	  b = Gnew.params->roworb[h][row][1];

	  /* Compute the active CC row index */
	  I = qt_occ[i];  B = qt_vir[b];
	  i = cc_occ[I];  b = cc_vir[B];

	  for(col=0; col < Gnew.params->coltot[h]; col++) {
	      j = Gnew.params->colorb[h][col][0];
	      a = Gnew.params->colorb[h][col][1];

	      /* Compute the active CC column index */
	      J = qt_occ[j];  A = qt_vir[a];
	      j = cc_occ[J];  a = cc_vir[A];

	      if(i<0 || j<0 || a<0 || b<0)
		  Gnew.matrix[h][row][col] = 0.0;
	      else {
		  ib = G.params->rowidx[i][b];
		  ja = G.params->colidx[j][a];

		  Gnew.matrix[h][row][col] = G.matrix[h][ib][ja];
		}
	    }
	}

      global_dpd_->buf4_mat_irrep_wrt(&Gnew, h);

      global_dpd_->buf4_mat_irrep_close(&G, h);
    }

  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_close(&Gnew);

  /* G(CI,AB) --> G'(CI,AB) */
  dpd_set_default(0);
  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 11, 7, 11, 7, 0, "GCIAB");

  dpd_set_default(1);
  global_dpd_->buf4_init(&Gnew, PSIF_CC_GAMMA_NEW, 0, 11, 7, 11, 7, 0, "GCIAB");

  for(h=0; h < nirreps; h++) {
      global_dpd_->buf4_mat_irrep_init(&G, h);
      global_dpd_->buf4_mat_irrep_rd(&G, h);

      global_dpd_->buf4_mat_irrep_init(&Gnew, h);

      for(row=0; row < Gnew.params->rowtot[h]; row++) {
	  c = Gnew.params->roworb[h][row][0];
	  i = Gnew.params->roworb[h][row][1];

	  /* Compute the active CC row index */

	  C = qt_vir[c];  I = qt_occ[i];
	  c = cc_vir[C];  i = cc_occ[I];

	  for(col=0; col < Gnew.params->coltot[h]; col++) {
	      a = Gnew.params->colorb[h][col][0];
	      b = Gnew.params->colorb[h][col][1];

	      /* Compute the active CC column index */
	      A = qt_vir[a];  B = qt_vir[b];
	      a = cc_vir[A];  b = cc_vir[B];

	      if(c<0 || i<0 || a<0 || b<0)
		  Gnew.matrix[h][row][col] = 0.0;
	      else {
		  ci = G.params->rowidx[c][i];
		  ab = G.params->colidx[a][b];

		  Gnew.matrix[h][row][col] = G.matrix[h][ci][ab];
		}
	    }
	}

      global_dpd_->buf4_mat_irrep_wrt(&Gnew, h);

      global_dpd_->buf4_mat_irrep_close(&G, h);
    }

  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_close(&Gnew);

  /* G(ci,ab) --> G'(ci,ab) */
  dpd_set_default(0);
  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 11, 7, 11, 7, 0, "Gciab");

  dpd_set_default(1);
  global_dpd_->buf4_init(&Gnew, PSIF_CC_GAMMA_NEW, 0, 11, 7, 11, 7, 0, "Gciab");

  for(h=0; h < nirreps; h++) {
      global_dpd_->buf4_mat_irrep_init(&G, h);
      global_dpd_->buf4_mat_irrep_rd(&G, h);

      global_dpd_->buf4_mat_irrep_init(&Gnew, h);

      for(row=0; row < Gnew.params->rowtot[h]; row++) {
	  c = Gnew.params->roworb[h][row][0];
	  i = Gnew.params->roworb[h][row][1];

	  /* Compute the active CC row index */

	  C = qt_vir[c];  I = qt_occ[i];
	  c = cc_vir[C];  i = cc_occ[I];

	  for(col=0; col < Gnew.params->coltot[h]; col++) {
	      a = Gnew.params->colorb[h][col][0];
	      b = Gnew.params->colorb[h][col][1];

	      /* Compute the active CC column index */
	      A = qt_vir[a];  B = qt_vir[b];
	      a = cc_vir[A];  b = cc_vir[B];

	      if(c<0 || i<0 || a<0 || b<0)
		  Gnew.matrix[h][row][col] = 0.0;
	      else {
		  ci = G.params->rowidx[c][i];
		  ab = G.params->colidx[a][b];

		  Gnew.matrix[h][row][col] = G.matrix[h][ci][ab];
		}
	    }
	}

      global_dpd_->buf4_mat_irrep_wrt(&Gnew, h);

      global_dpd_->buf4_mat_irrep_close(&G, h);
    }

  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_close(&Gnew);

  /* G(cI,aB) --> G'(cI,aB) */
  dpd_set_default(0);
  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 11, 5, 11, 5, 0, "GcIaB");

  dpd_set_default(1);
  global_dpd_->buf4_init(&Gnew, PSIF_CC_GAMMA_NEW, 0, 11, 5, 11, 5, 0, "GcIaB");

  for(h=0; h < nirreps; h++) {
      global_dpd_->buf4_mat_irrep_init(&G, h);
      global_dpd_->buf4_mat_irrep_rd(&G, h);

      global_dpd_->buf4_mat_irrep_init(&Gnew, h);

      for(row=0; row < Gnew.params->rowtot[h]; row++) {
	  c = Gnew.params->roworb[h][row][0];
	  i = Gnew.params->roworb[h][row][1];

	  /* Compute the active CC row index */

	  C = qt_vir[c];  I = qt_occ[i];
	  c = cc_vir[C];  i = cc_occ[I];

	  for(col=0; col < Gnew.params->coltot[h]; col++) {
	      a = Gnew.params->colorb[h][col][0];
	      b = Gnew.params->colorb[h][col][1];

	      /* Compute the active CC column index */
	      A = qt_vir[a];  B = qt_vir[b];
	      a = cc_vir[A];  b = cc_vir[B];

	      if(c<0 || i<0 || a<0 || b<0)
		  Gnew.matrix[h][row][col] = 0.0;
	      else {
		  ci = G.params->rowidx[c][i];
		  ab = G.params->colidx[a][b];

		  Gnew.matrix[h][row][col] = G.matrix[h][ci][ab];
		}
	    }
	}

      global_dpd_->buf4_mat_irrep_wrt(&Gnew, h);

      global_dpd_->buf4_mat_irrep_close(&G, h);
    }

  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_close(&Gnew);

  /* G(Ci,Ab) --> G'(Ci,Ab) */
  dpd_set_default(0);
  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 11, 5, 11, 5, 0, "GCiAb");

  dpd_set_default(1);
  global_dpd_->buf4_init(&Gnew, PSIF_CC_GAMMA_NEW, 0, 11, 5, 11, 5, 0, "GCiAb");

  for(h=0; h < nirreps; h++) {
      global_dpd_->buf4_mat_irrep_init(&G, h);
      global_dpd_->buf4_mat_irrep_rd(&G, h);

      global_dpd_->buf4_mat_irrep_init(&Gnew, h);

      for(row=0; row < Gnew.params->rowtot[h]; row++) {
	  c = Gnew.params->roworb[h][row][0];
	  i = Gnew.params->roworb[h][row][1];

	  /* Compute the active CC row index */

	  C = qt_vir[c];  I = qt_occ[i];
	  c = cc_vir[C];  i = cc_occ[I];

	  for(col=0; col < Gnew.params->coltot[h]; col++) {
	      a = Gnew.params->colorb[h][col][0];
	      b = Gnew.params->colorb[h][col][1];

	      /* Compute the active CC column index */
	      A = qt_vir[a];  B = qt_vir[b];
	      a = cc_vir[A];  b = cc_vir[B];

	      if(c<0 || i<0 || a<0 || b<0)
		  Gnew.matrix[h][row][col] = 0.0;
	      else {
		  ci = G.params->rowidx[c][i];
		  ab = G.params->colidx[a][b];

		  Gnew.matrix[h][row][col] = G.matrix[h][ci][ab];
		}
	    }
	}

      global_dpd_->buf4_mat_irrep_wrt(&Gnew, h);

      global_dpd_->buf4_mat_irrep_close(&G, h);
    }

  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_close(&Gnew);

  /* G(AB,CD) --> G'(AB,CD) */
  dpd_set_default(0);
  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 7, 7, 7, 7, 0, "GABCD");

  dpd_set_default(1);
  global_dpd_->buf4_init(&Gnew, PSIF_CC_GAMMA_NEW, 0, 7, 7, 7, 7, 0, "GABCD");

  for(h=0; h < nirreps; h++) {
      global_dpd_->buf4_mat_irrep_init(&G, h);
      global_dpd_->buf4_mat_irrep_rd(&G, h);

      global_dpd_->buf4_mat_irrep_init(&Gnew, h);

      for(row=0; row < Gnew.params->rowtot[h]; row++) {
	  a = Gnew.params->roworb[h][row][0];
	  b = Gnew.params->roworb[h][row][1];

	  /* Compute the active CC row index */

	  A = qt_vir[a];  B = qt_vir[b];
	  a = cc_vir[A];  b = cc_vir[B];

	  for(col=0; col < Gnew.params->coltot[h]; col++) {
	      c = Gnew.params->colorb[h][col][0];
	      d = Gnew.params->colorb[h][col][1];

	      /* Compute the active CC column index */
	      C = qt_vir[c];  D = qt_vir[d];
	      c = cc_vir[C];  d = cc_vir[D];

	      if(a<0 || b<0 || c<0 || d<0)
		  Gnew.matrix[h][row][col] = 0.0;
	      else {
		  ab = G.params->rowidx[a][b];
		  cd = G.params->colidx[c][d];

		  Gnew.matrix[h][row][col] = G.matrix[h][ab][cd];
		}
	    }
	}

      global_dpd_->buf4_mat_irrep_wrt(&Gnew, h);

      global_dpd_->buf4_mat_irrep_close(&G, h);
    }

  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_close(&Gnew);

  /* G(ab,cd) --> G'(ab,cd) */
  dpd_set_default(0);
  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 7, 7, 7, 7, 0, "Gabcd");

  dpd_set_default(1);
  global_dpd_->buf4_init(&Gnew, PSIF_CC_GAMMA_NEW, 0, 7, 7, 7, 7, 0, "Gabcd");

  for(h=0; h < nirreps; h++) {
      global_dpd_->buf4_mat_irrep_init(&G, h);
      global_dpd_->buf4_mat_irrep_rd(&G, h);

      global_dpd_->buf4_mat_irrep_init(&Gnew, h);

      for(row=0; row < Gnew.params->rowtot[h]; row++) {
	  a = Gnew.params->roworb[h][row][0];
	  b = Gnew.params->roworb[h][row][1];

	  /* Compute the active CC row index */

	  A = qt_vir[a];  B = qt_vir[b];
	  a = cc_vir[A];  b = cc_vir[B];

	  for(col=0; col < Gnew.params->coltot[h]; col++) {
	      c = Gnew.params->colorb[h][col][0];
	      d = Gnew.params->colorb[h][col][1];

	      /* Compute the active CC column index */
	      C = qt_vir[c];  D = qt_vir[d];
	      c = cc_vir[C];  d = cc_vir[D];

	      if(a<0 || b<0 || c<0 || d<0)
		  Gnew.matrix[h][row][col] = 0.0;
	      else {
		  ab = G.params->rowidx[a][b];
		  cd = G.params->colidx[c][d];

		  Gnew.matrix[h][row][col] = G.matrix[h][ab][cd];
		}
	    }
	}

      global_dpd_->buf4_mat_irrep_wrt(&Gnew, h);

      global_dpd_->buf4_mat_irrep_close(&G, h);
    }

  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_close(&Gnew);

  /* G(Ab,Cd) --> G'(Ab,Cd) */
  dpd_set_default(0);
  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 5, 5, 5, 5, 0, "GAbCd");

  dpd_set_default(1);
  global_dpd_->buf4_init(&Gnew, PSIF_CC_GAMMA_NEW, 0, 5, 5, 5, 5, 0, "GAbCd");

  for(h=0; h < nirreps; h++) {
      global_dpd_->buf4_mat_irrep_init(&G, h);
      global_dpd_->buf4_mat_irrep_rd(&G, h);

      global_dpd_->buf4_mat_irrep_init(&Gnew, h);

      for(row=0; row < Gnew.params->rowtot[h]; row++) {
	  a = Gnew.params->roworb[h][row][0];
	  b = Gnew.params->roworb[h][row][1];

	  /* Compute the active CC row index */

	  A = qt_vir[a];  B = qt_vir[b];
	  a = cc_vir[A];  b = cc_vir[B];

	  for(col=0; col < Gnew.params->coltot[h]; col++) {
	      c = Gnew.params->colorb[h][col][0];
	      d = Gnew.params->colorb[h][col][1];

	      /* Compute the active CC column index */
	      C = qt_vir[c];  D = qt_vir[d];
	      c = cc_vir[C];  d = cc_vir[D];

	      if(a<0 || b<0 || c<0 || d<0)
		  Gnew.matrix[h][row][col] = 0.0;
	      else {
		  ab = G.params->rowidx[a][b];
		  cd = G.params->colidx[c][d];

		  Gnew.matrix[h][row][col] = G.matrix[h][ab][cd];
		}
	    }
	}

      global_dpd_->buf4_mat_irrep_wrt(&Gnew, h);

      global_dpd_->buf4_mat_irrep_close(&G, h);
    }

  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_close(&Gnew);
}

}} // namespace psi::ccdensity
