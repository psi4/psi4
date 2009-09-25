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

  fprintf(outfile, "\n\tEnergies re-computed from Fock-adjusted CC density:\n");
  fprintf(outfile,   "\t---------------------------------------------------\n");

  dpd_file2_init(&D, CC_OEI, 0, 0, 0, "DIJ");
  dpd_file2_init(&F, CC_OEI, 0, 0, 0, "h(i,j)");
  this_energy = dpd_file2_dot(&D, &F);
  dpd_file2_close(&F);
  dpd_file2_close(&D);

  /*    fprintf(outfile, "\tDIJ = %20.15f\n", this_energy); */
  one_energy += this_energy;

  dpd_file2_init(&D, CC_OEI, 0, 0, 0, "Dij");
  dpd_file2_init(&F, CC_OEI, 0, 0, 0, "h(i,j)");
  this_energy = dpd_file2_dot(&D, &F);
  dpd_file2_close(&F);
  dpd_file2_close(&D);

  /*    fprintf(outfile, "\tDij = %20.15f\n", this_energy); */
  one_energy += this_energy;

  dpd_file2_init(&D, CC_OEI, 0, 1, 1, "DAB");
  dpd_file2_init(&F, CC_OEI, 0, 1, 1, "h(a,b)");
  this_energy = dpd_file2_dot(&D, &F);
  dpd_file2_close(&F);
  dpd_file2_close(&D);

  /*    fprintf(outfile, "\tDAB = %20.15f\n", this_energy); */
  one_energy += this_energy;

  dpd_file2_init(&D, CC_OEI, 0, 1, 1, "Dab");
  dpd_file2_init(&F, CC_OEI, 0, 1, 1, "h(a,b)");
  this_energy = dpd_file2_dot(&D, &F);
  dpd_file2_close(&F);
  dpd_file2_close(&D);

  /*    fprintf(outfile, "\tDab = %20.15f\n", this_energy); */
  one_energy += this_energy;

  dpd_file2_init(&D, CC_OEI, 0, 0, 1, "DIA");
  dpd_file2_init(&F, CC_OEI, 0, 0, 1, "h(i,a)");
  this_energy = dpd_file2_dot(&D, &F);
  dpd_file2_close(&F);
  dpd_file2_close(&D);

  /*    fprintf(outfile, "\tDIA = %20.15f\n", this_energy); */
  one_energy += this_energy;

  dpd_file2_init(&D, CC_OEI, 0, 0, 1, "Dia");
  dpd_file2_init(&F, CC_OEI, 0, 0, 1, "h(i,a)");
  this_energy = dpd_file2_dot(&D, &F);
  dpd_file2_close(&F);
  dpd_file2_close(&D);

  /*    fprintf(outfile, "\tDia = %20.15f\n", this_energy); */
  one_energy += this_energy;

  dpd_file2_init(&D, CC_OEI, 0, 0, 1, "DAI");
  dpd_file2_init(&F, CC_OEI, 0, 0, 1, "h(i,a)");
  this_energy = dpd_file2_dot(&D, &F);
  dpd_file2_close(&F);
    dpd_file2_close(&D);

  /*    fprintf(outfile, "\tDAI = %20.15f\n", this_energy); */
  one_energy += this_energy;

  dpd_file2_init(&D, CC_OEI, 0, 0, 1, "Dai");
  dpd_file2_init(&F, CC_OEI, 0, 0, 1, "h(i,a)");
  this_energy = dpd_file2_dot(&D, &F);
  dpd_file2_close(&F);
  dpd_file2_close(&D);

  /*    fprintf(outfile, "\tDai = %20.15f\n", this_energy); */
  one_energy += this_energy;

  fprintf(outfile, "\tOne-electron energy        = %20.15f\n", one_energy);
  fflush(outfile);

  dpd_file2_init(&D, CC_OEI, 0, 0, 0, "DIJ");
  dpd_file2_mat_init(&D);
  dpd_file2_mat_rd(&D);

  dpd_buf4_init(&G, CC_GAMMA, 0, 0, 0, 2, 2, 0, "GIJKL");
  dpd_buf4_scm(&G,0.0);
  dpd_buf4_close(&G);

  dpd_buf4_init(&G, CC_GAMMA, 0, 0, 0, 2, 2, 0, "GIJKL");
  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&G, h);
    dpd_buf4_mat_irrep_rd(&G, h);

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

    dpd_buf4_mat_irrep_wrt(&G, h);
    dpd_buf4_mat_irrep_close(&G, h);
  }
      
  two_energy = 0.0;
  dpd_buf4_init(&Aints, CC_AINTS, 0, 0, 0, 0, 0, 1, "A <ij|kl>");
  two_energy += 0.25 * dpd_buf4_dot(&Aints, &G);
  dpd_buf4_close(&Aints);

  dpd_buf4_close(&G);

  dpd_file2_mat_close(&D);
  dpd_file2_close(&D);

  dpd_file2_init(&D, CC_OEI, 0, 0, 0, "Dij");
  dpd_file2_mat_init(&D);
  dpd_file2_mat_rd(&D);

  dpd_buf4_init(&G, CC_GAMMA, 0, 0, 0, 2, 2, 0, "Gijkl");
  dpd_buf4_scm(&G,0.0);
  dpd_buf4_close(&G);

  dpd_buf4_init(&G, CC_GAMMA, 0, 0, 0, 2, 2, 0, "Gijkl");
  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&G, h);
    dpd_buf4_mat_irrep_rd(&G, h);

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

    dpd_buf4_mat_irrep_wrt(&G, h);
    dpd_buf4_mat_irrep_close(&G, h);
  }

  dpd_buf4_init(&Aints, CC_AINTS, 0, 0, 0, 0, 0, 1, "A <ij|kl>");
  two_energy += 0.25 * dpd_buf4_dot(&Aints, &G);
  dpd_buf4_close(&Aints);

  dpd_buf4_close(&G);

  dpd_file2_mat_close(&D);
  dpd_file2_close(&D);

  dpd_file2_init(&D, CC_OEI, 0, 0, 0, "DIJ");
  dpd_file2_mat_init(&D);
  dpd_file2_mat_rd(&D);

  dpd_buf4_init(&G, CC_GAMMA, 0, 0, 0, 0, 0, 0, "GIjKl");
  dpd_buf4_scm(&G,0.0);
  dpd_buf4_close(&G);

  dpd_buf4_init(&G, CC_GAMMA, 0, 0, 0, 0, 0, 0, "GIjKl");
  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&G, h);
    dpd_buf4_mat_irrep_rd(&G, h);

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

    dpd_buf4_mat_irrep_wrt(&G, h); 
    dpd_buf4_mat_irrep_close(&G, h);
  }
  dpd_buf4_close(&G);

  dpd_file2_mat_close(&D);
  dpd_file2_close(&D);


  dpd_file2_init(&D, CC_OEI, 0, 0, 0, "Dij");
  dpd_file2_mat_init(&D);
  dpd_file2_mat_rd(&D);

  dpd_buf4_init(&G, CC_GAMMA, 0, 0, 0, 0, 0, 0, "GIjKl");
  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&G, h);
    dpd_buf4_mat_irrep_rd(&G, h);

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

    dpd_buf4_mat_irrep_wrt(&G, h);
    dpd_buf4_mat_irrep_close(&G, h);
  }

  dpd_buf4_init(&Aints, CC_AINTS, 0, 0, 0, 0, 0, 0, "A <ij|kl>");
  two_energy += dpd_buf4_dot(&Aints, &G);
  dpd_buf4_close(&Aints);

  dpd_buf4_close(&G);

  total_two_energy += two_energy;
  fprintf(outfile, "\tIJKL energy                = %20.15f\n", two_energy);
  fflush(outfile);

  dpd_file2_mat_close(&D);
  dpd_file2_close(&D);

  dpd_file2_init(&D1, CC_OEI, 0, 0, 1, "DIA");
  dpd_file2_mat_init(&D1);
  dpd_file2_mat_rd(&D1);
  dpd_file2_init(&D2, CC_OEI, 0, 0, 1, "DAI");
  dpd_file2_mat_init(&D2);
  dpd_file2_mat_rd(&D2);

  dpd_buf4_init(&G, CC_GAMMA, 0, 0, 10, 2, 10, 0, "GIJKA");
  dpd_buf4_scm(&G,0.0);
  dpd_buf4_close(&G);

  dpd_buf4_init(&G, CC_GAMMA, 0, 0, 10, 2, 10, 0, "GIJKA");
  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&G, h);
    dpd_buf4_mat_irrep_rd(&G, h);

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

    dpd_buf4_mat_irrep_wrt(&G, h);
    dpd_buf4_mat_irrep_close(&G, h);
  }

  two_energy = 0.0;
  dpd_buf4_init(&E, CC_EINTS, 0, 0, 10, 2, 10, 0, "E <ij||ka> (i>j,ka)");
  two_energy += dpd_buf4_dot(&E, &G);
  dpd_buf4_close(&E);

  dpd_buf4_close(&G);

  dpd_file2_mat_close(&D1);
  dpd_file2_close(&D1);
  dpd_file2_mat_close(&D2);
  dpd_file2_close(&D2);

  dpd_file2_init(&D1, CC_OEI, 0, 0, 1, "Dia");
  dpd_file2_mat_init(&D1);
  dpd_file2_mat_rd(&D1);
  dpd_file2_init(&D2, CC_OEI, 0, 0, 1, "Dai");
  dpd_file2_mat_init(&D2);
  dpd_file2_mat_rd(&D2);

  dpd_buf4_init(&G, CC_GAMMA, 0, 0, 10, 2, 10, 0, "Gijka");
  dpd_buf4_scm(&G,0.0);
  dpd_buf4_close(&G);

  dpd_buf4_init(&G, CC_GAMMA, 0, 0, 10, 2, 10, 0, "Gijka");
  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&G, h);
    dpd_buf4_mat_irrep_rd(&G, h);

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

    dpd_buf4_mat_irrep_wrt(&G, h);
    dpd_buf4_mat_irrep_close(&G, h);
  }

  dpd_buf4_init(&E, CC_EINTS, 0, 0, 10, 2, 10, 0, "E <ij||ka> (i>j,ka)");
  two_energy += dpd_buf4_dot(&E, &G);
  dpd_buf4_close(&E);

  dpd_buf4_close(&G);

  dpd_file2_mat_close(&D1);
  dpd_file2_close(&D1);
  dpd_file2_mat_close(&D2);
  dpd_file2_close(&D2);

  dpd_file2_init(&D1, CC_OEI, 0, 0, 1, "DIA");
  dpd_file2_mat_init(&D1);
  dpd_file2_mat_rd(&D1);
  dpd_file2_init(&D2, CC_OEI, 0, 0, 1, "DAI");
  dpd_file2_mat_init(&D2);
  dpd_file2_mat_rd(&D2);

  dpd_buf4_init(&G, CC_GAMMA, 0, 0, 10, 0, 10, 0, "GiJkA");
  dpd_buf4_scm(&G,0.0);
  dpd_buf4_close(&G);

  dpd_buf4_init(&G, CC_GAMMA, 0, 0, 10, 0, 10, 0, "GiJkA");
  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&G, h);
    dpd_buf4_mat_irrep_rd(&G, h);

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

    dpd_buf4_mat_irrep_wrt(&G, h);
    dpd_buf4_mat_irrep_close(&G, h);
  }

  dpd_buf4_init(&E, CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
  two_energy += 2 * dpd_buf4_dot(&E, &G);
  dpd_buf4_close(&E);

  dpd_buf4_close(&G);

  dpd_file2_mat_close(&D1);
  dpd_file2_close(&D1);
  dpd_file2_mat_close(&D2);
  dpd_file2_close(&D2);

  dpd_file2_init(&D1, CC_OEI, 0, 0, 1, "Dia");
  dpd_file2_mat_init(&D1);
  dpd_file2_mat_rd(&D1);
  dpd_file2_init(&D2, CC_OEI, 0, 0, 1, "Dai");
  dpd_file2_mat_init(&D2);
  dpd_file2_mat_rd(&D2);

  dpd_buf4_init(&G, CC_GAMMA, 0, 0, 10, 0, 10, 0, "GIjKa");
  dpd_buf4_scm(&G,0.0);
  dpd_buf4_close(&G);

  dpd_buf4_init(&G, CC_GAMMA, 0, 0, 10, 0, 10, 0, "GIjKa");
  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&G, h);
    dpd_buf4_mat_irrep_rd(&G, h);

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

    dpd_buf4_mat_irrep_wrt(&G, h);
    dpd_buf4_mat_irrep_close(&G, h);
  }

  dpd_buf4_init(&E, CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
  two_energy += 2 * dpd_buf4_dot(&E, &G);
  dpd_buf4_close(&E);

  dpd_buf4_close(&G);

  total_two_energy += two_energy;
  fprintf(outfile, "\tIJKA energy                = %20.15f\n", two_energy);
  fflush(outfile);

  dpd_file2_mat_close(&D1);
  dpd_file2_close(&D1);
  dpd_file2_mat_close(&D2);
  dpd_file2_close(&D2);

  two_energy = 0.0;
  dpd_buf4_init(&DInts, CC_DINTS, 0, 2, 7, 2, 7, 0, "D <ij||ab> (i>j,a>b)");
  dpd_buf4_init(&G, CC_GAMMA, 0, 2, 7, 2, 7, 0, "GIJAB");
  two_energy += dpd_buf4_dot(&G, &DInts);
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, CC_GAMMA, 0, 2, 7, 2, 7, 0, "Gijab");
  two_energy += dpd_buf4_dot(&G, &DInts); 
  dpd_buf4_close(&G);
  dpd_buf4_close(&DInts);
  dpd_buf4_init(&DInts, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
  dpd_buf4_init(&G, CC_GAMMA, 0, 0, 5, 0, 5, 0, "GIjAb");
  two_energy += dpd_buf4_dot(&G, &DInts);
  dpd_buf4_close(&G);
  dpd_buf4_close(&DInts);

  two_energy *= 2;
  total_two_energy += two_energy;
  fprintf(outfile, "\tIJAB energy                = %20.15f\n", two_energy);
  fflush(outfile);

  dpd_file2_init(&D, CC_OEI, 0, 1, 1, "DAB");
  dpd_file2_mat_init(&D);
  dpd_file2_mat_rd(&D);

  dpd_buf4_init(&G, CC_GAMMA, 0, 10, 10, 10, 10, 0, "GIBJA");
  dpd_buf4_scm(&G, 0.0);
  dpd_buf4_close(&G);

  dpd_buf4_init(&G, CC_GAMMA, 0, 10, 10, 10, 10, 0, "GIBJA");
  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&G, h);
    dpd_buf4_mat_irrep_rd(&G, h);

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

    dpd_buf4_mat_irrep_wrt(&G, h);
    dpd_buf4_mat_irrep_close(&G, h);
  }

  two_energy = 0.0;
  dpd_buf4_init(&C, CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia||jb>");
  two_energy += dpd_buf4_dot(&C, &G);
  dpd_buf4_close(&C);

  dpd_buf4_close(&G);

  dpd_file2_mat_close(&D);
  dpd_file2_close(&D);


  dpd_file2_init(&D, CC_OEI, 0, 1, 1, "Dab");
  dpd_file2_mat_init(&D);
  dpd_file2_mat_rd(&D);

  dpd_buf4_init(&G, CC_GAMMA, 0, 10, 10, 10, 10, 0, "Gibja");
  dpd_buf4_scm(&G, 0.0);
  dpd_buf4_close(&G);

  dpd_buf4_init(&G, CC_GAMMA, 0, 10, 10, 10, 10, 0, "Gibja");
  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&G, h);
    dpd_buf4_mat_irrep_rd(&G, h);

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

    dpd_buf4_mat_irrep_wrt(&G, h);
    dpd_buf4_mat_irrep_close(&G, h);
  }

  dpd_buf4_init(&C, CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia||jb>");
  two_energy += dpd_buf4_dot(&C, &G);
  dpd_buf4_close(&C);

  dpd_buf4_close(&G);

  dpd_file2_mat_close(&D);
  dpd_file2_close(&D);

  dpd_file2_init(&D, CC_OEI, 0, 1, 1, "Dab");
  dpd_file2_mat_init(&D);
  dpd_file2_mat_rd(&D);

  dpd_buf4_init(&G, CC_GAMMA, 0, 10, 10, 10, 10, 0, "GIbJa");
  dpd_buf4_scm(&G, 0.0);
  dpd_buf4_close(&G);

  dpd_buf4_init(&G, CC_GAMMA, 0, 10, 10, 10, 10, 0, "GIbJa");
  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&G, h);
    dpd_buf4_mat_irrep_rd(&G, h);

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

    dpd_buf4_mat_irrep_wrt(&G, h);
    dpd_buf4_mat_irrep_close(&G, h);
  }

  dpd_buf4_init(&C, CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb>");
  two_energy += dpd_buf4_dot(&C, &G);
  dpd_buf4_close(&C);

  dpd_buf4_close(&G);

  dpd_file2_mat_close(&D);
  dpd_file2_close(&D);

  dpd_file2_init(&D, CC_OEI, 0, 1, 1, "DAB");
  dpd_file2_mat_init(&D);
  dpd_file2_mat_rd(&D);

  dpd_buf4_init(&G, CC_GAMMA, 0, 10, 10, 10, 10, 0, "GiBjA");
  dpd_buf4_scm(&G, 0.0);
  dpd_buf4_close(&G);

  dpd_buf4_init(&G, CC_GAMMA, 0, 10, 10, 10, 10, 0, "GiBjA");
  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&G, h);
    dpd_buf4_mat_irrep_rd(&G, h);

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
    dpd_buf4_mat_irrep_wrt(&G, h);
    dpd_buf4_mat_irrep_close(&G, h);
  }

  dpd_buf4_init(&C, CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb>");
  two_energy += dpd_buf4_dot(&C, &G); 
  dpd_buf4_close(&C);

  dpd_buf4_close(&G);

  dpd_file2_mat_close(&D);
  dpd_file2_close(&D);

  dpd_buf4_init(&DInts, CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ib,ja)");
  dpd_buf4_init(&G, CC_GAMMA, 0, 10, 10, 10, 10, 0, "GIbjA");
  two_energy -= dpd_buf4_dot(&G, &DInts);
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, CC_GAMMA, 0, 10, 10, 10, 10, 0, "GiBJa");
  two_energy -= dpd_buf4_dot(&G, &DInts);
  dpd_buf4_close(&G);
  dpd_buf4_close(&DInts);

  total_two_energy += two_energy;
  fprintf(outfile, "\tIBJA energy                = %20.15f\n", two_energy);
  fflush(outfile);

  fprintf(outfile, "\tMP2 correlation energy    = %20.15f\n",
	  one_energy + total_two_energy);
  fprintf(outfile, "\tTotal MP2 energy          = %20.15f\n",
	  one_energy + total_two_energy + mo.Escf);
}

void uhf_sf_fold(void)
{

}

}} /* End namespaces */
