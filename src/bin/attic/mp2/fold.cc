/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
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
 *@END LICENSE
 */

/*! \file
    \ingroup MP2
    \brief Enter brief description of file here 
*/
#include <libdpd/dpd.h>
#define EXTERN
#include "globals.h"

namespace psi{ namespace mp2{

/* FOLD(): Fold the Fock matrix contributions to the energy (or energy
** derivative) into the two-particle density matrix.  Here we are
** trying to convert from an energy expression of the form:
**
** E = sum_pq Dpq fpq + 1/4 sum_pqrs Gpqrs <pq||rs>
**
** to the form:
**
** E = sum_pq Dpq hpq + 1/4 sum_pqrs Gpqrs <pq||rs>
**
** We do this by shifting some one-particle density matrix components
** into appropriate two-particle density matrix components:
**
** G'pmrm = Dpr + 4 * Gpmrm
**
** One problem is that we need to make sure the resulting density,
** G'pqrs, is still antisymmetric to permutation of p and q or r and
** s.  So, for example, for the Gimkm component we compute:
**
** G'pmrm = Dpr + Gpmrm
** G'mprm = Dpr - Gmprm
** G'pmmr = Dpr - Gpmmr
** G'mpmr = Dpr + Gmpmr
** */

void rhf_sf_fold(void);
void uhf_sf_fold(void);

void fold(void)
{
  if(params.ref == 0) rhf_sf_fold();
  else if(params.ref == 2) uhf_sf_fold();
}

void rhf_sf_fold(void)
{
  int h, nirreps;
  int i, j, k, l, m, a, b;
  int I, J, K, L, M, A, B;
  int IM, JM, MI, MJ, MK, ML, MA, MB;
  int Gi, Gj, Gk, Gl, Gm, Ga, Gb;
  int *occpi, *virtpi;
  int *occ_off, *vir_off;
  int *occ_sym, *vir_sym;
  dpdfile2 D, D1, D2, F;
  dpdbuf4 G, Aints, E, C, DInts, FInts, BInts, G1, G2;
  double one_energy=0.0, two_energy=0.0, total_two_energy=0.0;
  double test_energy = 0.0, tmp;
  double this_energy;

  nirreps = mo.nirreps;
  occpi = mo.occpi; 
  virtpi = mo.virtpi;
  occ_off = mo.occ_off; 
  vir_off = mo.vir_off;
  occ_sym = mo.occ_sym; 
  vir_sym = mo.vir_sym;

  outfile->Printf( "\n\tEnergies re-computed from Fock-adjusted CC density:\n");
  outfile->Printf(   "\t---------------------------------------------------\n");

  global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 0, "DIJ");
  global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, 0, 0, "h(i,j)");
  this_energy = global_dpd_->file2_dot(&D, &F);
  global_dpd_->file2_close(&F);
  global_dpd_->file2_close(&D);

  /*    outfile->Printf( "\tDIJ = %20.15f\n", this_energy); */
  one_energy += this_energy;

  global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 0, "Dij");
  global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, 0, 0, "h(i,j)");
  this_energy = global_dpd_->file2_dot(&D, &F);
  global_dpd_->file2_close(&F);
  global_dpd_->file2_close(&D);

  /*    outfile->Printf( "\tDij = %20.15f\n", this_energy); */
  one_energy += this_energy;

  global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 1, 1, "DAB");
  global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, 1, 1, "h(a,b)");
  this_energy = global_dpd_->file2_dot(&D, &F);
  global_dpd_->file2_close(&F);
  global_dpd_->file2_close(&D);

  /*    outfile->Printf( "\tDAB = %20.15f\n", this_energy); */
  one_energy += this_energy;

  global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 1, 1, "Dab");
  global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, 1, 1, "h(a,b)");
  this_energy = global_dpd_->file2_dot(&D, &F);
  global_dpd_->file2_close(&F);
  global_dpd_->file2_close(&D);

  /*    outfile->Printf( "\tDab = %20.15f\n", this_energy); */
  one_energy += this_energy;

  global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 1, "DIA");
  global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, 0, 1, "h(i,a)");
  this_energy = global_dpd_->file2_dot(&D, &F);
  global_dpd_->file2_close(&F);
  global_dpd_->file2_close(&D);

  /*    outfile->Printf( "\tDIA = %20.15f\n", this_energy); */
  one_energy += this_energy;

  global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 1, "Dia");
  global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, 0, 1, "h(i,a)");
  this_energy = global_dpd_->file2_dot(&D, &F);
  global_dpd_->file2_close(&F);
  global_dpd_->file2_close(&D);

  /*    outfile->Printf( "\tDia = %20.15f\n", this_energy); */
  one_energy += this_energy;

  global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 1, "DAI");
  global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, 0, 1, "h(i,a)");
  this_energy = global_dpd_->file2_dot(&D, &F);
  global_dpd_->file2_close(&F);
    global_dpd_->file2_close(&D);

  /*    outfile->Printf( "\tDAI = %20.15f\n", this_energy); */
  one_energy += this_energy;

  global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 1, "Dai");
  global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, 0, 1, "h(i,a)");
  this_energy = global_dpd_->file2_dot(&D, &F);
  global_dpd_->file2_close(&F);
  global_dpd_->file2_close(&D);

  /*    outfile->Printf( "\tDai = %20.15f\n", this_energy); */
  one_energy += this_energy;

  outfile->Printf( "\tOne-electron energy        = %20.15f\n", one_energy);
  

  global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 0, "DIJ");
  global_dpd_->file2_mat_init(&D);
  global_dpd_->file2_mat_rd(&D);

  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 0, 0, 2, 2, 0, "GIJKL");
  global_dpd_->buf4_scm(&G,0.0);
  global_dpd_->buf4_close(&G);

  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 0, 0, 2, 2, 0, "GIJKL");
  for(h=0; h < nirreps; h++) {
    global_dpd_->buf4_mat_irrep_init(&G, h);
    global_dpd_->buf4_mat_irrep_rd(&G, h);

    for(Gm=0; Gm < nirreps; Gm++) {
      Gi = Gj = h^Gm;

      for(i=0; i < occpi[Gi]; i++) {
	I = occ_off[Gi] + i;
	for(j=0; j < occpi[Gj]; j++) {
	  J = occ_off[Gj] + j;
	  for(m=0; m < occpi[Gm]; m++) {
	    M = occ_off[Gm] + m;

	    IM = G.params->rowidx[I][M];
	    JM = G.params->colidx[J][M];
	    MI = G.params->rowidx[M][I];
	    MJ = G.params->colidx[M][J];

	    G.matrix[h][IM][JM] += D.matrix[Gi][i][j];
	    G.matrix[h][IM][MJ] -= D.matrix[Gi][i][j];
	    G.matrix[h][MI][MJ] += D.matrix[Gi][i][j];
	    G.matrix[h][MI][JM] -= D.matrix[Gi][i][j];
	  }
	}
      }
    }

    global_dpd_->buf4_mat_irrep_wrt(&G, h);
    global_dpd_->buf4_mat_irrep_close(&G, h);
  }
      
  two_energy = 0.0;
  global_dpd_->buf4_init(&Aints, PSIF_CC_AINTS, 0, 0, 0, 0, 0, 1, "A <ij|kl>");
  two_energy += 0.25 * global_dpd_->buf4_dot(&Aints, &G);
  global_dpd_->buf4_close(&Aints);

  global_dpd_->buf4_close(&G);

  global_dpd_->file2_mat_close(&D);
  global_dpd_->file2_close(&D);

  global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 0, "Dij");
  global_dpd_->file2_mat_init(&D);
  global_dpd_->file2_mat_rd(&D);

  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 0, 0, 2, 2, 0, "Gijkl");
  global_dpd_->buf4_scm(&G,0.0);
  global_dpd_->buf4_close(&G);

  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 0, 0, 2, 2, 0, "Gijkl");
  for(h=0; h < nirreps; h++) {
    global_dpd_->buf4_mat_irrep_init(&G, h);
    global_dpd_->buf4_mat_irrep_rd(&G, h);

    for(Gm=0; Gm < nirreps; Gm++) {
      Gi = Gj = h^Gm;

      for(i=0; i < occpi[Gi]; i++) {
	I = occ_off[Gi] + i;
	for(j=0; j < occpi[Gj]; j++) {
	  J = occ_off[Gj] + j;
	  for(m=0; m < occpi[Gm]; m++) {
	    M = occ_off[Gm] + m;

	    IM = G.params->rowidx[I][M];
	    JM = G.params->colidx[J][M];
	    MI = G.params->rowidx[M][I];
	    MJ = G.params->colidx[M][J];

	    G.matrix[h][IM][JM] += D.matrix[Gi][i][j];
	    G.matrix[h][MI][JM] -= D.matrix[Gi][i][j];
	    G.matrix[h][MI][MJ] += D.matrix[Gi][i][j];
	    G.matrix[h][IM][MJ] -= D.matrix[Gi][i][j];
	  }
	}
      }
    }

    global_dpd_->buf4_mat_irrep_wrt(&G, h);
    global_dpd_->buf4_mat_irrep_close(&G, h);
  }

  global_dpd_->buf4_init(&Aints, PSIF_CC_AINTS, 0, 0, 0, 0, 0, 1, "A <ij|kl>");
  two_energy += 0.25 * global_dpd_->buf4_dot(&Aints, &G);
  global_dpd_->buf4_close(&Aints);

  global_dpd_->buf4_close(&G);

  global_dpd_->file2_mat_close(&D);
  global_dpd_->file2_close(&D);

  global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 0, "DIJ");
  global_dpd_->file2_mat_init(&D);
  global_dpd_->file2_mat_rd(&D);

  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 0, 0, 0, 0, 0, "GIjKl");
  global_dpd_->buf4_scm(&G,0.0);
  global_dpd_->buf4_close(&G);

  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 0, 0, 0, 0, 0, "GIjKl");
  for(h=0; h < nirreps; h++) {
    global_dpd_->buf4_mat_irrep_init(&G, h);
    global_dpd_->buf4_mat_irrep_rd(&G, h);

    for(Gm=0; Gm < nirreps; Gm++) {
      Gi = Gj = h^Gm;

      for(i=0; i < occpi[Gi]; i++) {
	I = occ_off[Gi] + i;
	for(j=0; j < occpi[Gj]; j++) {
	  J = occ_off[Gj] + j;
	  for(m=0; m < occpi[Gm]; m++) {
	    M = occ_off[Gm] + m;

	    IM = G.params->rowidx[I][M];
	    JM = G.params->colidx[J][M];

	    G.matrix[h][IM][JM] += D.matrix[Gi][i][j];
	  }
	}
      }
    }

    global_dpd_->buf4_mat_irrep_wrt(&G, h); 
    global_dpd_->buf4_mat_irrep_close(&G, h);
  }
  global_dpd_->buf4_close(&G);

  global_dpd_->file2_mat_close(&D);
  global_dpd_->file2_close(&D);


  global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 0, "Dij");
  global_dpd_->file2_mat_init(&D);
  global_dpd_->file2_mat_rd(&D);

  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 0, 0, 0, 0, 0, "GIjKl");
  for(h=0; h < nirreps; h++) {
    global_dpd_->buf4_mat_irrep_init(&G, h);
    global_dpd_->buf4_mat_irrep_rd(&G, h);

    for(Gm=0; Gm < nirreps; Gm++) {
      Gk = Gl = h^Gm;

      for(k=0; k < occpi[Gk]; k++) {
	K = occ_off[Gk] + k;
	for(l=0; l < occpi[Gl]; l++) {
	  L = occ_off[Gl] + l;
	  for(m=0; m < occpi[Gm]; m++) {
	    M = occ_off[Gm] + m;

	    MK = G.params->rowidx[M][K];
	    ML = G.params->colidx[M][L];

	    G.matrix[h][MK][ML] += D.matrix[Gk][k][l];
	  }
	}
      }
    }

    global_dpd_->buf4_mat_irrep_wrt(&G, h);
    global_dpd_->buf4_mat_irrep_close(&G, h);
  }

  global_dpd_->buf4_init(&Aints, PSIF_CC_AINTS, 0, 0, 0, 0, 0, 0, "A <ij|kl>");
  two_energy += global_dpd_->buf4_dot(&Aints, &G);
  global_dpd_->buf4_close(&Aints);

  global_dpd_->buf4_close(&G);

  total_two_energy += two_energy;
  outfile->Printf( "\tIJKL energy                = %20.15f\n", two_energy);
  

  global_dpd_->file2_mat_close(&D);
  global_dpd_->file2_close(&D);

  global_dpd_->file2_init(&D1, PSIF_CC_OEI, 0, 0, 1, "DIA");
  global_dpd_->file2_mat_init(&D1);
  global_dpd_->file2_mat_rd(&D1);
  global_dpd_->file2_init(&D2, PSIF_CC_OEI, 0, 0, 1, "DAI");
  global_dpd_->file2_mat_init(&D2);
  global_dpd_->file2_mat_rd(&D2);

  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 0, 10, 2, 10, 0, "GIJKA");
  global_dpd_->buf4_scm(&G,0.0);
  global_dpd_->buf4_close(&G);

  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 0, 10, 2, 10, 0, "GIJKA");
  for(h=0; h < nirreps; h++) {
    global_dpd_->buf4_mat_irrep_init(&G, h);
    global_dpd_->buf4_mat_irrep_rd(&G, h);

    for(Gm=0; Gm < nirreps; Gm++) {
      Gi = Ga = h^Gm;

      for(i=0; i < occpi[Gi]; i++) {
	I = occ_off[Gi] + i;
	for(a=0; a < virtpi[Ga]; a++) {
	  A = vir_off[Ga] + a;
	  for(m=0; m < occpi[Gm]; m++) {
	    M = occ_off[Gm] + m;

	    MI = G.params->rowidx[M][I];
	    IM = G.params->rowidx[I][M];
	    MA = G.params->colidx[M][A];

	    G.matrix[h][MI][MA] += 0.5 * (D1.matrix[Gi][i][a] +
					  D2.matrix[Gi][i][a]);
	    G.matrix[h][IM][MA] -= 0.5 * (D1.matrix[Gi][i][a] +
					  D2.matrix[Gi][i][a]);
	  }
	}
      }
    }

    global_dpd_->buf4_mat_irrep_wrt(&G, h);
    global_dpd_->buf4_mat_irrep_close(&G, h);
  }

  two_energy = 0.0;
  global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 0, 10, 2, 10, 0, "E <ij||ka> (i>j,ka)");
  two_energy += global_dpd_->buf4_dot(&E, &G);
  global_dpd_->buf4_close(&E);

  global_dpd_->buf4_close(&G);

  global_dpd_->file2_mat_close(&D1);
  global_dpd_->file2_close(&D1);
  global_dpd_->file2_mat_close(&D2);
  global_dpd_->file2_close(&D2);

  global_dpd_->file2_init(&D1, PSIF_CC_OEI, 0, 0, 1, "Dia");
  global_dpd_->file2_mat_init(&D1);
  global_dpd_->file2_mat_rd(&D1);
  global_dpd_->file2_init(&D2, PSIF_CC_OEI, 0, 0, 1, "Dai");
  global_dpd_->file2_mat_init(&D2);
  global_dpd_->file2_mat_rd(&D2);

  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 0, 10, 2, 10, 0, "Gijka");
  global_dpd_->buf4_scm(&G,0.0);
  global_dpd_->buf4_close(&G);

  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 0, 10, 2, 10, 0, "Gijka");
  for(h=0; h < nirreps; h++) {
    global_dpd_->buf4_mat_irrep_init(&G, h);
    global_dpd_->buf4_mat_irrep_rd(&G, h);

    for(Gm=0; Gm < nirreps; Gm++) {
      Gi = Ga = h^Gm;

      for(i=0; i < occpi[Gi]; i++) {
	I = occ_off[Gi] + i;
	for(a=0; a < virtpi[Ga]; a++) {
	  A = vir_off[Ga] + a;
	  for(m=0; m < occpi[Gm]; m++) {
	    M = occ_off[Gm] + m;

	    MI = G.params->rowidx[M][I];
	    IM = G.params->rowidx[I][M];
	    MA = G.params->colidx[M][A];

	    G.matrix[h][MI][MA] += 0.5 * (D1.matrix[Gi][i][a] +
					  D2.matrix[Gi][i][a]);
	    G.matrix[h][IM][MA] -= 0.5 * (D1.matrix[Gi][i][a] +
					  D2.matrix[Gi][i][a]);
	  }
	}
      }
    }

    global_dpd_->buf4_mat_irrep_wrt(&G, h);
    global_dpd_->buf4_mat_irrep_close(&G, h);
  }

  global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 0, 10, 2, 10, 0, "E <ij||ka> (i>j,ka)");
  two_energy += global_dpd_->buf4_dot(&E, &G);
  global_dpd_->buf4_close(&E);

  global_dpd_->buf4_close(&G);

  global_dpd_->file2_mat_close(&D1);
  global_dpd_->file2_close(&D1);
  global_dpd_->file2_mat_close(&D2);
  global_dpd_->file2_close(&D2);

  global_dpd_->file2_init(&D1, PSIF_CC_OEI, 0, 0, 1, "DIA");
  global_dpd_->file2_mat_init(&D1);
  global_dpd_->file2_mat_rd(&D1);
  global_dpd_->file2_init(&D2, PSIF_CC_OEI, 0, 0, 1, "DAI");
  global_dpd_->file2_mat_init(&D2);
  global_dpd_->file2_mat_rd(&D2);

  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 0, 10, 0, 10, 0, "GiJkA");
  global_dpd_->buf4_scm(&G,0.0);
  global_dpd_->buf4_close(&G);

  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 0, 10, 0, 10, 0, "GiJkA");
  for(h=0; h < nirreps; h++) {
    global_dpd_->buf4_mat_irrep_init(&G, h);
    global_dpd_->buf4_mat_irrep_rd(&G, h);

    for(Gm=0; Gm < nirreps; Gm++) {
      Gi = Ga = h^Gm;

      for(i=0; i < occpi[Gi]; i++) {
	I = occ_off[Gi] + i;
	for(a=0; a < virtpi[Ga]; a++) {
	  A = vir_off[Ga] + a;
	  for(m=0; m < occpi[Gm]; m++) {
	    M = occ_off[Gm] + m;

	    MI = G.params->rowidx[M][I];
	    MA = G.params->colidx[M][A];

	    G.matrix[h][MI][MA] += 0.5 * (D1.matrix[Gi][i][a] +
					  D2.matrix[Gi][i][a]);
	  }
	}
      }
    }

    global_dpd_->buf4_mat_irrep_wrt(&G, h);
    global_dpd_->buf4_mat_irrep_close(&G, h);
  }

  global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
  two_energy += 2 * global_dpd_->buf4_dot(&E, &G);
  global_dpd_->buf4_close(&E);

  global_dpd_->buf4_close(&G);

  global_dpd_->file2_mat_close(&D1);
  global_dpd_->file2_close(&D1);
  global_dpd_->file2_mat_close(&D2);
  global_dpd_->file2_close(&D2);

  global_dpd_->file2_init(&D1, PSIF_CC_OEI, 0, 0, 1, "Dia");
  global_dpd_->file2_mat_init(&D1);
  global_dpd_->file2_mat_rd(&D1);
  global_dpd_->file2_init(&D2, PSIF_CC_OEI, 0, 0, 1, "Dai");
  global_dpd_->file2_mat_init(&D2);
  global_dpd_->file2_mat_rd(&D2);

  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 0, 10, 0, 10, 0, "GIjKa");
  global_dpd_->buf4_scm(&G,0.0);
  global_dpd_->buf4_close(&G);

  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 0, 10, 0, 10, 0, "GIjKa");
  for(h=0; h < nirreps; h++) {
    global_dpd_->buf4_mat_irrep_init(&G, h);
    global_dpd_->buf4_mat_irrep_rd(&G, h);

    for(Gm=0; Gm < nirreps; Gm++) {
      Gi = Ga = h^Gm;

      for(i=0; i < occpi[Gi]; i++) {
	I = occ_off[Gi] + i;
	for(a=0; a < virtpi[Ga]; a++) {
	  A = vir_off[Ga] + a;
	  for(m=0; m < occpi[Gm]; m++) {
	    M = occ_off[Gm] + m;

	    MI = G.params->rowidx[M][I];
	    MA = G.params->colidx[M][A];

	    G.matrix[h][MI][MA] += 0.5 * (D1.matrix[Gi][i][a] +
					  D2.matrix[Gi][i][a]);
	  }
	}
      }
    }

    global_dpd_->buf4_mat_irrep_wrt(&G, h);
    global_dpd_->buf4_mat_irrep_close(&G, h);
  }

  global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
  two_energy += 2 * global_dpd_->buf4_dot(&E, &G);
  global_dpd_->buf4_close(&E);

  global_dpd_->buf4_close(&G);

  total_two_energy += two_energy;
  outfile->Printf( "\tIJKA energy                = %20.15f\n", two_energy);
  

  global_dpd_->file2_mat_close(&D1);
  global_dpd_->file2_close(&D1);
  global_dpd_->file2_mat_close(&D2);
  global_dpd_->file2_close(&D2);

  two_energy = 0.0;
  global_dpd_->buf4_init(&DInts, PSIF_CC_DINTS, 0, 2, 7, 2, 7, 0, "D <ij||ab> (i>j,a>b)");
  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 2, 7, 2, 7, 0, "GIJAB");
  two_energy += global_dpd_->buf4_dot(&G, &DInts);
  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 2, 7, 2, 7, 0, "Gijab");
  two_energy += global_dpd_->buf4_dot(&G, &DInts); 
  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_close(&DInts);
  global_dpd_->buf4_init(&DInts, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 0, 5, 0, 5, 0, "GIjAb");
  two_energy += global_dpd_->buf4_dot(&G, &DInts);
  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_close(&DInts);

  two_energy *= 2;
  total_two_energy += two_energy;
  outfile->Printf( "\tIJAB energy                = %20.15f\n", two_energy);
  

  global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 1, 1, "DAB");
  global_dpd_->file2_mat_init(&D);
  global_dpd_->file2_mat_rd(&D);

  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 10, 10, 10, 10, 0, "GIBJA");
  global_dpd_->buf4_scm(&G, 0.0);
  global_dpd_->buf4_close(&G);

  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 10, 10, 10, 10, 0, "GIBJA");
  for(h=0; h < nirreps; h++) {
    global_dpd_->buf4_mat_irrep_init(&G, h);
    global_dpd_->buf4_mat_irrep_rd(&G, h);

    for(Gm=0; Gm < nirreps; Gm++) {
      Ga = Gb = h^Gm;

      for(b=0; b < virtpi[Gb]; b++) {
	B = vir_off[Gb] + b;
	for(a=0; a < virtpi[Ga]; a++) {
	  A = vir_off[Ga] + a;
	  for(m=0; m < occpi[Gm]; m++) {
	    M = occ_off[Gm] + m;

	    MB = G.params->rowidx[M][B];
	    MA = G.params->colidx[M][A];

	    G.matrix[h][MB][MA] += D.matrix[Ga][a][b];
	  }
	}
      }
    }

    global_dpd_->buf4_mat_irrep_wrt(&G, h);
    global_dpd_->buf4_mat_irrep_close(&G, h);
  }

  two_energy = 0.0;
  global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia||jb>");
  two_energy += global_dpd_->buf4_dot(&C, &G);
  global_dpd_->buf4_close(&C);

  global_dpd_->buf4_close(&G);

  global_dpd_->file2_mat_close(&D);
  global_dpd_->file2_close(&D);


  global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 1, 1, "Dab");
  global_dpd_->file2_mat_init(&D);
  global_dpd_->file2_mat_rd(&D);

  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 10, 10, 10, 10, 0, "Gibja");
  global_dpd_->buf4_scm(&G, 0.0);
  global_dpd_->buf4_close(&G);

  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 10, 10, 10, 10, 0, "Gibja");
  for(h=0; h < nirreps; h++) {
    global_dpd_->buf4_mat_irrep_init(&G, h);
    global_dpd_->buf4_mat_irrep_rd(&G, h);

    for(Gm=0; Gm < nirreps; Gm++) {
      Ga = Gb = h^Gm;

      for(b=0; b < virtpi[Gb]; b++) {
	B = vir_off[Gb] + b;
	for(a=0; a < virtpi[Ga]; a++) {
	  A = vir_off[Ga] + a;
	  for(m=0; m < occpi[Gm]; m++) {
	    M = occ_off[Gm] + m;

	    MB = G.params->rowidx[M][B];
	    MA = G.params->colidx[M][A];

	    G.matrix[h][MB][MA] += D.matrix[Ga][a][b];
	  }
	}
      }
    }

    global_dpd_->buf4_mat_irrep_wrt(&G, h);
    global_dpd_->buf4_mat_irrep_close(&G, h);
  }

  global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia||jb>");
  two_energy += global_dpd_->buf4_dot(&C, &G);
  global_dpd_->buf4_close(&C);

  global_dpd_->buf4_close(&G);

  global_dpd_->file2_mat_close(&D);
  global_dpd_->file2_close(&D);

  global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 1, 1, "Dab");
  global_dpd_->file2_mat_init(&D);
  global_dpd_->file2_mat_rd(&D);

  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 10, 10, 10, 10, 0, "GIbJa");
  global_dpd_->buf4_scm(&G, 0.0);
  global_dpd_->buf4_close(&G);

  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 10, 10, 10, 10, 0, "GIbJa");
  for(h=0; h < nirreps; h++) {
    global_dpd_->buf4_mat_irrep_init(&G, h);
    global_dpd_->buf4_mat_irrep_rd(&G, h);

    for(Gm=0; Gm < nirreps; Gm++) {
      Ga = Gb = h^Gm;

      for(b=0; b < virtpi[Gb]; b++) {
	B = vir_off[Gb] + b;
	for(a=0; a < virtpi[Ga]; a++) {
	  A = vir_off[Ga] + a;
	  for(m=0; m < occpi[Gm]; m++) {
	    M = occ_off[Gm] + m;

	    MB = G.params->rowidx[M][B];
	    MA = G.params->colidx[M][A];

	    G.matrix[h][MB][MA] += D.matrix[Ga][a][b];
	  }
	}
      }
    }

    global_dpd_->buf4_mat_irrep_wrt(&G, h);
    global_dpd_->buf4_mat_irrep_close(&G, h);
  }

  global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb>");
  two_energy += global_dpd_->buf4_dot(&C, &G);
  global_dpd_->buf4_close(&C);

  global_dpd_->buf4_close(&G);

  global_dpd_->file2_mat_close(&D);
  global_dpd_->file2_close(&D);

  global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 1, 1, "DAB");
  global_dpd_->file2_mat_init(&D);
  global_dpd_->file2_mat_rd(&D);

  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 10, 10, 10, 10, 0, "GiBjA");
  global_dpd_->buf4_scm(&G, 0.0);
  global_dpd_->buf4_close(&G);

  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 10, 10, 10, 10, 0, "GiBjA");
  for(h=0; h < nirreps; h++) {
    global_dpd_->buf4_mat_irrep_init(&G, h);
    global_dpd_->buf4_mat_irrep_rd(&G, h);

    for(Gm=0; Gm < nirreps; Gm++) {
      Ga = Gb = h^Gm;

      for(b=0; b < virtpi[Gb]; b++) {
	B = vir_off[Gb] + b;
	for(a=0; a < virtpi[Ga]; a++) {
	  A = vir_off[Ga] + a;
	  for(m=0; m < occpi[Gm]; m++) {
	    M = occ_off[Gm] + m;

	    MB = G.params->rowidx[M][B];
	    MA = G.params->colidx[M][A];

	    G.matrix[h][MB][MA] += D.matrix[Ga][a][b];
	  }
	}
      }
    }
    global_dpd_->buf4_mat_irrep_wrt(&G, h);
    global_dpd_->buf4_mat_irrep_close(&G, h);
  }

  global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb>");
  two_energy += global_dpd_->buf4_dot(&C, &G); 
  global_dpd_->buf4_close(&C);

  global_dpd_->buf4_close(&G);

  global_dpd_->file2_mat_close(&D);
  global_dpd_->file2_close(&D);

  global_dpd_->buf4_init(&DInts, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ib,ja)");
  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 10, 10, 10, 10, 0, "GIbjA");
  two_energy -= global_dpd_->buf4_dot(&G, &DInts);
  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 10, 10, 10, 10, 0, "GiBJa");
  two_energy -= global_dpd_->buf4_dot(&G, &DInts);
  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_close(&DInts);

  total_two_energy += two_energy;
  outfile->Printf( "\tIBJA energy                = %20.15f\n", two_energy);
  

  outfile->Printf( "\tMP2 correlation energy    = %20.15f\n",
	  one_energy + total_two_energy);
  outfile->Printf( "\tTotal MP2 energy          = %20.15f\n",
	  one_energy + total_two_energy + mo.Escf);
}

void uhf_sf_fold(void)
{

}

}} /* End namespaces */
