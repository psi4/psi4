/*! \file
    \ingroup CCDENSITY
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <libdpd/dpd.h>
#include "MOInfo.h"
#include "Params.h"
#include "Frozen.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccdensity {

/* FOLD_ROHF(): Fold the ROHF Fock matrix contributions to the energy
** (or energy derivative) into the two-particle density matrix.  Here
** we are trying to convert from an energy expression of the form:
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

void fold_ROHF(struct RHO_Params rho_params)
{
  int h, nirreps;
  int i, j, k, l, m, a, b;
  int I, J, K, L, M, A, B;
  int IM, JM, MI, MJ, MK, ML, MA, MB;
  int Gi, Gj, Gk, Gl, Gm, Ga, Gb;
  int *occpi, *virtpi;
  int *occ_off, *vir_off;
  int *occ_sym, *vir_sym;
  int *openpi;
  dpdfile2 D, D1, D2, F;
  dpdbuf4 G, Aints, E, C, DInts, FInts, BInts, G1, G2;
  double one_energy=0.0, two_energy=0.0, total_two_energy=0.0;
  double test_energy = 0.0, tmp;
  double this_energy;

  nirreps = moinfo.nirreps;
  occpi = moinfo.occpi; virtpi = moinfo.virtpi;
  occ_off = moinfo.occ_off; vir_off = moinfo.vir_off;
  occ_sym = moinfo.occ_sym; vir_sym = moinfo.vir_sym;
  openpi = moinfo.openpi;

  if(!params.aobasis) {
    fprintf(outfile, "\n\tEnergies re-computed from Fock-adjusted CC density:\n");
    fprintf(outfile,   "\t---------------------------------------------------\n");

    dpd_file2_init(&D, CC_OEI, 0, 0, 0, rho_params.DIJ_lbl);
    dpd_file2_init(&F, CC_OEI, 0, 0, 0, "h(i,j)");
    this_energy = dpd_file2_dot(&D, &F);
    dpd_file2_close(&F);
    dpd_file2_close(&D);

    /*    fprintf(outfile, "\tDIJ = %20.15f\n", this_energy); */
    one_energy += this_energy;

    dpd_file2_init(&D, CC_OEI, 0, 0, 0, rho_params.Dij_lbl);
    dpd_file2_init(&F, CC_OEI, 0, 0, 0, "h(i,j)");
    this_energy = dpd_file2_dot(&D, &F);
    dpd_file2_close(&F);
    dpd_file2_close(&D);

    /*    fprintf(outfile, "\tDij = %20.15f\n", this_energy); */
    one_energy += this_energy;

    dpd_file2_init(&D, CC_OEI, 0, 1, 1, rho_params.DAB_lbl);
    dpd_file2_init(&F, CC_OEI, 0, 1, 1, "h(a,b)");
    this_energy = dpd_file2_dot(&D, &F);
    dpd_file2_close(&F);
    dpd_file2_close(&D);

    /*    fprintf(outfile, "\tDAB = %20.15f\n", this_energy); */
    one_energy += this_energy;

    dpd_file2_init(&D, CC_OEI, 0, 1, 1, rho_params.Dab_lbl);
    dpd_file2_init(&F, CC_OEI, 0, 1, 1, "h(a,b)");
    this_energy = dpd_file2_dot(&D, &F);
    dpd_file2_close(&F);
    dpd_file2_close(&D);

    /*    fprintf(outfile, "\tDab = %20.15f\n", this_energy); */
    one_energy += this_energy;

    dpd_file2_init(&D, CC_OEI, 0, 0, 1, rho_params.DIA_lbl);
    dpd_file2_init(&F, CC_OEI, 0, 0, 1, "h(i,a)");
    this_energy = dpd_file2_dot(&D, &F);
    dpd_file2_close(&F);
    dpd_file2_close(&D);

    /*    fprintf(outfile, "\tDIA = %20.15f\n", this_energy); */
    one_energy += this_energy;

    dpd_file2_init(&D, CC_OEI, 0, 0, 1, rho_params.Dia_lbl);
    dpd_file2_init(&F, CC_OEI, 0, 0, 1, "h(i,a)");
    this_energy = dpd_file2_dot(&D, &F);
    dpd_file2_close(&F);
    dpd_file2_close(&D);

    /*    fprintf(outfile, "\tDia = %20.15f\n", this_energy); */
    one_energy += this_energy;

    dpd_file2_init(&D, CC_OEI, 0, 0, 1, rho_params.DAI_lbl);
    dpd_file2_init(&F, CC_OEI, 0, 0, 1, "h(i,a)");
    this_energy = dpd_file2_dot(&D, &F);
    dpd_file2_close(&F);
    dpd_file2_close(&D);

    /*    fprintf(outfile, "\tDAI = %20.15f\n", this_energy); */
    one_energy += this_energy;

    dpd_file2_init(&D, CC_OEI, 0, 0, 1, rho_params.Dai_lbl);
    dpd_file2_init(&F, CC_OEI, 0, 0, 1, "h(i,a)");
    this_energy = dpd_file2_dot(&D, &F);
    dpd_file2_close(&F);
    dpd_file2_close(&D);

    /*    fprintf(outfile, "\tDai = %20.15f\n", this_energy); */
    one_energy += this_energy;

    fprintf(outfile, "\tOne-electron energy        = %20.15f\n", one_energy);
    fflush(outfile);
  }

  dpd_file2_init(&D, CC_OEI, 0, 0, 0, rho_params.DIJ_lbl);
  dpd_file2_mat_init(&D);
  dpd_file2_mat_rd(&D);

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
      
  if(!params.aobasis) {
    two_energy = 0.0;
    dpd_buf4_init(&Aints, CC_AINTS, 0, 0, 0, 0, 0, 1, "A <ij|kl>");
    two_energy += 0.25 * dpd_buf4_dot(&Aints, &G);
    dpd_buf4_close(&Aints);
  }

  dpd_buf4_close(&G);

  dpd_file2_mat_close(&D);
  dpd_file2_close(&D);

  dpd_file2_init(&D, CC_OEI, 0, 0, 0, rho_params.Dij_lbl);
  dpd_file2_mat_init(&D);
  dpd_file2_mat_rd(&D);

  dpd_buf4_init(&G, CC_GAMMA, 0, 0, 0, 2, 2, 0, "Gijkl");
  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&G, h);
    dpd_buf4_mat_irrep_rd(&G, h);

    for(Gm=0; Gm < nirreps; Gm++) {
      Gi = Gj = h^Gm;

      for(i=0; i < (occpi[Gi] - openpi[Gi]); i++) {
	I = occ_off[Gi] + i;
	for(j=0; j < (occpi[Gj] - openpi[Gj]); j++) {
	  J = occ_off[Gj] + j;
	  for(m=0; m < (occpi[Gm] - openpi[Gm]); m++) {
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

  if(!params.aobasis) {
    dpd_buf4_init(&Aints, CC_AINTS, 0, 0, 0, 0, 0, 1, "A <ij|kl>");
    two_energy += 0.25 * dpd_buf4_dot(&Aints, &G);
    dpd_buf4_close(&Aints);
  }

  dpd_buf4_close(&G);

  dpd_file2_mat_close(&D);
  dpd_file2_close(&D);

  dpd_file2_init(&D, CC_OEI, 0, 0, 0, rho_params.DIJ_lbl);
  dpd_file2_mat_init(&D);
  dpd_file2_mat_rd(&D);

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
	  for(m=0; m < (occpi[Gm] - openpi[Gm]); m++) {
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


  dpd_file2_init(&D, CC_OEI, 0, 0, 0, rho_params.Dij_lbl);
  dpd_file2_mat_init(&D);
  dpd_file2_mat_rd(&D);

  dpd_buf4_init(&G, CC_GAMMA, 0, 0, 0, 0, 0, 0, "GIjKl");
  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&G, h);
    dpd_buf4_mat_irrep_rd(&G, h);

    for(Gm=0; Gm < nirreps; Gm++) {
      Gk = Gl = h^Gm;

      for(k=0; k < (occpi[Gk] - openpi[Gk]); k++) {
	K = occ_off[Gk] + k;
	for(l=0; l < (occpi[Gl] - openpi[Gl]); l++) {
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

  if(!params.aobasis) {
    dpd_buf4_init(&Aints, CC_AINTS, 0, 0, 0, 0, 0, 0, "A <ij|kl>");
    two_energy += dpd_buf4_dot(&Aints, &G);
    dpd_buf4_close(&Aints);
  }

  dpd_buf4_close(&G);

  if(!params.aobasis) {
    total_two_energy += two_energy;
    fprintf(outfile, "\tIJKL energy                = %20.15f\n", two_energy);
    fflush(outfile);
  }

  dpd_file2_mat_close(&D);
  dpd_file2_close(&D);

  dpd_file2_init(&D1, CC_OEI, 0, 0, 1, rho_params.DIA_lbl);
  dpd_file2_mat_init(&D1);
  dpd_file2_mat_rd(&D1);
  dpd_file2_init(&D2, CC_OEI, 0, 0, 1, rho_params.DAI_lbl);
  dpd_file2_mat_init(&D2);
  dpd_file2_mat_rd(&D2);

  dpd_buf4_init(&G, CC_GAMMA, 0, 0, 10, 2, 10, 0, "GIJKA");
  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&G, h);
    dpd_buf4_mat_irrep_rd(&G, h);

    for(Gm=0; Gm < nirreps; Gm++) {
      Gi = Ga = h^Gm;

      for(i=0; i < occpi[Gi]; i++) {
	I = occ_off[Gi] + i;
	for(a=0; a < (virtpi[Ga] - openpi[Ga]); a++) {
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

  if(!params.aobasis) {
    two_energy = 0.0;
    dpd_buf4_init(&E, CC_EINTS, 0, 0, 10, 2, 10, 0, "E <ij||ka> (i>j,ka)");
    two_energy += dpd_buf4_dot(&E, &G);
    dpd_buf4_close(&E);
  }

  dpd_buf4_close(&G);

  dpd_file2_mat_close(&D1);
  dpd_file2_close(&D1);
  dpd_file2_mat_close(&D2);
  dpd_file2_close(&D2);

  dpd_file2_init(&D1, CC_OEI, 0, 0, 1, rho_params.Dia_lbl);
  dpd_file2_mat_init(&D1);
  dpd_file2_mat_rd(&D1);
  dpd_file2_init(&D2, CC_OEI, 0, 0, 1, rho_params.Dai_lbl);
  dpd_file2_mat_init(&D2);
  dpd_file2_mat_rd(&D2);

  dpd_buf4_init(&G, CC_GAMMA, 0, 0, 10, 2, 10, 0, "Gijka");
  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&G, h);
    dpd_buf4_mat_irrep_rd(&G, h);

    for(Gm=0; Gm < nirreps; Gm++) {
      Gi = Ga = h^Gm;

      for(i=0; i < (occpi[Gi] - openpi[Gi]); i++) {
	I = occ_off[Gi] + i;
	for(a=0; a < virtpi[Ga]; a++) {
	  A = vir_off[Ga] + a;
	  for(m=0; m < (occpi[Gm] - openpi[Gm]); m++) {
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
  if(!params.aobasis) {
    dpd_buf4_init(&E, CC_EINTS, 0, 0, 10, 2, 10, 0, "E <ij||ka> (i>j,ka)");
    two_energy += dpd_buf4_dot(&E, &G);
    dpd_buf4_close(&E);
  }

  dpd_buf4_close(&G);

  dpd_file2_mat_close(&D1);
  dpd_file2_close(&D1);
  dpd_file2_mat_close(&D2);
  dpd_file2_close(&D2);

  dpd_file2_init(&D1, CC_OEI, 0, 0, 1, rho_params.DIA_lbl);
  dpd_file2_mat_init(&D1);
  dpd_file2_mat_rd(&D1);
  dpd_file2_init(&D2, CC_OEI, 0, 0, 1, rho_params.DAI_lbl);
  dpd_file2_mat_init(&D2);
  dpd_file2_mat_rd(&D2);

  dpd_buf4_init(&G, CC_GAMMA, 0, 0, 10, 0, 10, 0, "GiJkA");
  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&G, h);
    dpd_buf4_mat_irrep_rd(&G, h);

    for(Gm=0; Gm < nirreps; Gm++) {
      Gi = Ga = h^Gm;

      for(i=0; i < occpi[Gi]; i++) {
	I = occ_off[Gi] + i;
	for(a=0; a < (virtpi[Ga] - openpi[Ga]); a++) {
	  A = vir_off[Ga] + a;
	  for(m=0; m < (occpi[Gm] - openpi[Gm]); m++) {
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

  if(!params.aobasis) {
    dpd_buf4_init(&E, CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
    two_energy += 2 * dpd_buf4_dot(&E, &G);
    dpd_buf4_close(&E);
  }

  dpd_buf4_close(&G);

  dpd_file2_mat_close(&D1);
  dpd_file2_close(&D1);
  dpd_file2_mat_close(&D2);
  dpd_file2_close(&D2);


  dpd_file2_init(&D1, CC_OEI, 0, 0, 1, rho_params.Dia_lbl);
  dpd_file2_mat_init(&D1);
  dpd_file2_mat_rd(&D1);
  dpd_file2_init(&D2, CC_OEI, 0, 0, 1, rho_params.Dai_lbl);
  dpd_file2_mat_init(&D2);
  dpd_file2_mat_rd(&D2);

  dpd_buf4_init(&G, CC_GAMMA, 0, 0, 10, 0, 10, 0, "GIjKa");
  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&G, h);
    dpd_buf4_mat_irrep_rd(&G, h);

    for(Gm=0; Gm < nirreps; Gm++) {
      Gi = Ga = h^Gm;

      for(i=0; i < (occpi[Gi] - openpi[Gi]); i++) {
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

  if(!params.aobasis) {
    dpd_buf4_init(&E, CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
    two_energy += 2 * dpd_buf4_dot(&E, &G);
    dpd_buf4_close(&E);
  }

  dpd_buf4_close(&G);

  if(!params.aobasis) {
    total_two_energy += two_energy;
    fprintf(outfile, "\tIJKA energy                = %20.15f\n", two_energy);
    fflush(outfile);
  }

  dpd_file2_mat_close(&D1);
  dpd_file2_close(&D1);
  dpd_file2_mat_close(&D2);
  dpd_file2_close(&D2);

  if(!params.aobasis) {
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
  }

  dpd_file2_init(&D, CC_OEI, 0, 1, 1, rho_params.DAB_lbl);
  dpd_file2_mat_init(&D);
  dpd_file2_mat_rd(&D);

  dpd_buf4_init(&G, CC_GAMMA, 0, 10, 10, 10, 10, 0, "GIBJA");
  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&G, h);
    dpd_buf4_mat_irrep_rd(&G, h);

    for(Gm=0; Gm < nirreps; Gm++) {
      Ga = Gb = h^Gm;

      for(b=0; b < (virtpi[Gb] - openpi[Gb]); b++) {
	B = vir_off[Gb] + b;
	for(a=0; a < (virtpi[Ga] - openpi[Ga]); a++) {
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

  if(!params.aobasis) {
    two_energy = 0.0;
    dpd_buf4_init(&C, CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia||jb>");
    two_energy += dpd_buf4_dot(&C, &G);
    dpd_buf4_close(&C);
  }

  dpd_buf4_close(&G);

  dpd_file2_mat_close(&D);
  dpd_file2_close(&D);


  dpd_file2_init(&D, CC_OEI, 0, 1, 1, rho_params.Dab_lbl);
  dpd_file2_mat_init(&D);
  dpd_file2_mat_rd(&D);

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
	  for(m=0; m < (occpi[Gm] - openpi[Gm]); m++) {
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

  if(!params.aobasis) {
    dpd_buf4_init(&C, CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia||jb>");
    two_energy += dpd_buf4_dot(&C, &G);
    dpd_buf4_close(&C);
  }

  dpd_buf4_close(&G);

  dpd_file2_mat_close(&D);
  dpd_file2_close(&D);


  dpd_file2_init(&D, CC_OEI, 0, 1, 1, rho_params.Dab_lbl);
  dpd_file2_mat_init(&D);
  dpd_file2_mat_rd(&D);

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

  if(!params.aobasis) {
    dpd_buf4_init(&C, CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb>");
    two_energy += dpd_buf4_dot(&C, &G);
    dpd_buf4_close(&C);
  }

  dpd_buf4_close(&G);

  dpd_file2_mat_close(&D);
  dpd_file2_close(&D);


  dpd_file2_init(&D, CC_OEI, 0, 1, 1, rho_params.DAB_lbl);
  dpd_file2_mat_init(&D);
  dpd_file2_mat_rd(&D);

  dpd_buf4_init(&G, CC_GAMMA, 0, 10, 10, 10, 10, 0, "GiBjA");
  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&G, h);
    dpd_buf4_mat_irrep_rd(&G, h);

    for(Gm=0; Gm < nirreps; Gm++) {
      Ga = Gb = h^Gm;

      for(b=0; b < (virtpi[Gb] - openpi[Gb]); b++) {
	B = vir_off[Gb] + b;
	for(a=0; a < (virtpi[Ga] - openpi[Ga]); a++) {
	  A = vir_off[Ga] + a;
	  for(m=0; m < (occpi[Gm] - openpi[Gm]); m++) {
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

  if(!params.aobasis) {
    dpd_buf4_init(&C, CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb>");
    two_energy += dpd_buf4_dot(&C, &G); 
    dpd_buf4_close(&C);
  }

  dpd_buf4_close(&G);

  dpd_file2_mat_close(&D);
  dpd_file2_close(&D);


  if(!params.aobasis) {
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
  }

  if(!params.aobasis) {
    two_energy = 0.0;
    dpd_buf4_init(&FInts, CC_FINTS, 0, 10, 7, 10, 5, 1, "F <ia|bc>");
    dpd_buf4_sort(&FInts, CC_TMP0, qprs, 11, 7, "F(CI,AB)");
    dpd_buf4_close(&FInts);
    dpd_buf4_init(&FInts, CC_TMP0, 0, 11, 7, 11, 7, 0, "F(CI,AB)");
    dpd_buf4_init(&G, CC_GAMMA, 0, 11, 7, 11, 7, 0, "GCIAB");
    two_energy = -1.0 * dpd_buf4_dot(&G, &FInts);
    dpd_buf4_close(&G);
    dpd_buf4_init(&G, CC_GAMMA, 0, 11, 7, 11, 7, 0, "Gciab");
    two_energy -= dpd_buf4_dot(&G, &FInts); 
    dpd_buf4_close(&G);
    dpd_buf4_close(&FInts);
    dpd_buf4_init(&FInts, CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
    dpd_buf4_sort(&FInts, CC_TMP0, qprs, 11, 5, "F(cI,Ba)");
    dpd_buf4_close(&FInts);
    dpd_buf4_init(&FInts, CC_TMP0, 0, 11, 5, 11, 5, 0, "F(cI,Ba)");
    dpd_buf4_sort(&FInts, CC_TMP1, pqsr, 11, 5, "F(cI,aB)");
    dpd_buf4_close(&FInts);
    dpd_buf4_init(&FInts, CC_TMP1, 0, 11, 5, 11, 5, 0, "F(cI,aB)");
    dpd_buf4_init(&G, CC_GAMMA, 0, 11, 5, 11, 5, 0, "GcIaB");
    two_energy += dpd_buf4_dot(&G, &FInts);
    dpd_buf4_close(&G);
    dpd_buf4_init(&G, CC_GAMMA, 0, 11, 5, 11, 5, 0, "GCiAb");
    two_energy += dpd_buf4_dot(&G, &FInts); 
    dpd_buf4_close(&G);
    dpd_buf4_close(&FInts);

    two_energy *= 2;
    total_two_energy += two_energy;
    fprintf(outfile, "\tCIAB energy                = %20.15f\n", two_energy);
    fflush(outfile);
  }

  if(!params.aobasis) {
    two_energy = 0.0;
    dpd_buf4_init(&BInts, CC_BINTS, 0, 7, 7, 5, 5, 1, "B <ab|cd>");
    dpd_buf4_init(&G, CC_GAMMA, 0, 7, 7, 7, 7, 0, "GABCD");
    two_energy += dpd_buf4_dot(&G, &BInts);
    dpd_buf4_close(&G);
    dpd_buf4_init(&G, CC_GAMMA, 0, 7, 7, 7, 7, 0, "Gabcd");
    two_energy += dpd_buf4_dot(&G, &BInts);
    dpd_buf4_close(&G);
    dpd_buf4_close(&BInts);
    dpd_buf4_init(&BInts, CC_BINTS, 0, 5, 5, 5, 5, 0, "B <ab|cd>");
    dpd_buf4_init(&G, CC_GAMMA, 0, 5, 5, 5, 5, 0, "GAbCd");
    two_energy += dpd_buf4_dot(&G, &BInts);
    dpd_buf4_close(&G);
    dpd_buf4_close(&BInts);
  }

  if(!params.aobasis) {
    total_two_energy += two_energy;
    fprintf(outfile, "\tABCD energy                = %20.15f\n", two_energy);
    fprintf(outfile, "\tTotal two-electron energy  = %20.15f\n", total_two_energy);
    if (params.ground) {
      fprintf(outfile, "\tCCSD correlation energy    = %20.15f\n",
	      one_energy + total_two_energy);
      fprintf(outfile, "\tTotal CCSD energy          = %20.15f\n",
	      one_energy + total_two_energy + moinfo.eref);
    }
    else {
      fprintf(outfile, "\tTotal EOM CCSD correlation energy        = %20.15f\n",
          one_energy + total_two_energy);
      fprintf(outfile, "\tCCSD correlation + EOM excitation energy = %20.15f\n",
          moinfo.ecc + params.cceom_energy);
      fprintf(outfile, "\tTotal EOM CCSD energy                    = %20.15f\n",
          one_energy + total_two_energy + moinfo.eref);
    }
  }
}

}} // namespace psi::ccdensity
