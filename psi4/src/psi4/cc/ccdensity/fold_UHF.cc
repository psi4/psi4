/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2018 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This file is part of Psi4.
 *
 * Psi4 is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * Psi4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with Psi4; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

/*! \file
    \ingroup CCDENSITY
    \brief Enter brief description of file here
*/
#include <cstdio>
#include <cstring>
#include "psi4/libdpd/dpd.h"
#include "MOInfo.h"
#include "Params.h"
#include "Frozen.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccdensity {

/* FOLD_UHF(): Fold the UHF Fock matrix contributions to the energy
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

void fold_UHF(struct RHO_Params rho_params)
{
  int h, nirreps;
  int i, j, k, l, m, a, b;
  int I, J, K, L, M, A, B;
  int IM, JM, MI, MJ, MK, ML, MA, MB;
  int Gi, Gj, Gk, Gl, Gm, Ga, Gb;
  int *aoccpi, *avirtpi;
  int *boccpi, *bvirtpi;
  int *aocc_off, *avir_off;
  int *bocc_off, *bvir_off;
  int *aocc_sym, *avir_sym;
  int *bocc_sym, *bvir_sym;
  dpdfile2 D, D1, D2, F;
  dpdbuf4 G, Aints, E, C, DInts, FInts, BInts, G1, G2;
  double one_energy=0.0, two_energy=0.0, total_two_energy=0.0;
  double this_energy;

  nirreps = moinfo.nirreps;
  aoccpi = moinfo.aoccpi; avirtpi = moinfo.avirtpi;
  boccpi = moinfo.boccpi; bvirtpi = moinfo.bvirtpi;
  aocc_off = moinfo.aocc_off; avir_off = moinfo.avir_off;
  bocc_off = moinfo.bocc_off; bvir_off = moinfo.bvir_off;
  aocc_sym = moinfo.aocc_sym; avir_sym = moinfo.avir_sym;
  bocc_sym = moinfo.bocc_sym; bvir_sym = moinfo.bvir_sym;


  if(params.debug_) {
    outfile->Printf( "\n\tEnergies re-computed from Fock-adjusted CC density:\n");
    outfile->Printf(   "\t---------------------------------------------------\n");

    global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 0, rho_params.DIJ_lbl);
    global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, 0, 0, "h(I,J)");
    this_energy = global_dpd_->file2_dot(&D, &F);
    global_dpd_->file2_close(&F);
    global_dpd_->file2_close(&D);
    /*  outfile->Printf( "\tDIJ = %20.15f\n", this_energy); */
    one_energy += this_energy;

    global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 2, 2, rho_params.Dij_lbl);
    global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, 2, 2, "h(i,j)");
    this_energy = global_dpd_->file2_dot(&D, &F);
    global_dpd_->file2_close(&F);
    global_dpd_->file2_close(&D);

    /*  outfile->Printf( "\tDij = %20.15f\n", this_energy); */
    one_energy += this_energy;

    global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 1, 1, rho_params.DAB_lbl);
    global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, 1, 1, "h(A,B)");
    this_energy = global_dpd_->file2_dot(&D, &F);
    global_dpd_->file2_close(&F);
    global_dpd_->file2_close(&D);

    /*  outfile->Printf( "\tDAB = %20.15f\n", this_energy); */
    one_energy += this_energy;

    global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 3, 3, rho_params.Dab_lbl);
    global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, 3, 3, "h(a,b)");
    this_energy = global_dpd_->file2_dot(&D, &F);
    global_dpd_->file2_close(&F);
    global_dpd_->file2_close(&D);

    /*  outfile->Printf( "\tDab = %20.15f\n", this_energy); */
    one_energy += this_energy;

    global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 1, rho_params.DIA_lbl);
    global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, 0, 1, "h(I,A)");
    this_energy = global_dpd_->file2_dot(&D, &F);
    global_dpd_->file2_close(&F);
    global_dpd_->file2_close(&D);

    /*  outfile->Printf( "\tDIA = %20.15f\n", this_energy); */
    one_energy += this_energy;

    global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 2, 3, rho_params.Dia_lbl);
    global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, 2, 3, "h(i,a)");
    this_energy = global_dpd_->file2_dot(&D, &F);
    global_dpd_->file2_close(&F);
    global_dpd_->file2_close(&D);

    /*  outfile->Printf( "\tDia = %20.15f\n", this_energy); */
    one_energy += this_energy;

    global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 1, rho_params.DAI_lbl);
    global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, 0, 1, "h(I,A)");
    this_energy = global_dpd_->file2_dot(&D, &F);
    global_dpd_->file2_close(&F);
    global_dpd_->file2_close(&D);

    /*  outfile->Printf( "\tDAI = %20.15f\n", this_energy); */
    one_energy += this_energy;

    global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 2, 3, rho_params.Dai_lbl);
    global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, 2, 3, "h(i,a)");
    this_energy = global_dpd_->file2_dot(&D, &F);
    global_dpd_->file2_close(&F);
    global_dpd_->file2_close(&D);

    /*  outfile->Printf( "\tDai = %20.15f\n", this_energy); */
    one_energy += this_energy;

    outfile->Printf( "\tOne-electron energy        = %20.15f\n", one_energy);
  }

  global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 0, rho_params.DIJ_lbl);
  global_dpd_->file2_mat_init(&D);
  global_dpd_->file2_mat_rd(&D);

  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 0, 0, 2, 2, 0, "GIJKL");
  for(h=0; h < nirreps; h++) {
    global_dpd_->buf4_mat_irrep_init(&G, h);
    global_dpd_->buf4_mat_irrep_rd(&G, h);

    for(Gm=0; Gm < nirreps; Gm++) {
      Gi = Gj = h^Gm;

      for(i=0; i < aoccpi[Gi]; i++) {
	    I = aocc_off[Gi] + i;
	    for(j=0; j < aoccpi[Gj]; j++) {
	      J = aocc_off[Gj] + j;
	      for(m=0; m < aoccpi[Gm]; m++) {
	        M = aocc_off[Gm] + m;

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

  if(params.debug_) {
    two_energy = 0.0;
    global_dpd_->buf4_init(&Aints, PSIF_CC_AINTS, 0, 0, 0, 0, 0, 1, "A <IJ|KL>");
    two_energy += 0.25 * global_dpd_->buf4_dot(&Aints, &G);
    global_dpd_->buf4_close(&Aints);
  }

  global_dpd_->buf4_close(&G);

  global_dpd_->file2_mat_close(&D);
  global_dpd_->file2_close(&D);

  global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 2, 2, rho_params.Dij_lbl);
  global_dpd_->file2_mat_init(&D);
  global_dpd_->file2_mat_rd(&D);

  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 10, 10, 12, 12, 0, "Gijkl");
  for(h=0; h < nirreps; h++) {
    global_dpd_->buf4_mat_irrep_init(&G, h);
    global_dpd_->buf4_mat_irrep_rd(&G, h);

    for(Gm=0; Gm < nirreps; Gm++) {
      Gi = Gj = h^Gm;
      for(i=0; i < boccpi[Gi]; i++) {
	    I = bocc_off[Gi] + i;
	    for(j=0; j < boccpi[Gj]; j++) {
	      J = bocc_off[Gj] + j;
	      for(m=0; m < boccpi[Gm]; m++) {
	        M = bocc_off[Gm] + m;

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

  if(params.debug_){
    global_dpd_->buf4_init(&Aints, PSIF_CC_AINTS, 0, 10, 10, 10, 10, 1, "A <ij|kl>");
    two_energy += 0.25 * global_dpd_->buf4_dot(&Aints, &G);
    global_dpd_->buf4_close(&Aints);
  }
  global_dpd_->buf4_close(&G);

  global_dpd_->file2_mat_close(&D);
  global_dpd_->file2_close(&D);

  global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 0, rho_params.DIJ_lbl);
  global_dpd_->file2_mat_init(&D);
  global_dpd_->file2_mat_rd(&D);

  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 22, 22, 22, 22, 0, "GIjKl");
  for(h=0; h < nirreps; h++) {
    global_dpd_->buf4_mat_irrep_init(&G, h);
    global_dpd_->buf4_mat_irrep_rd(&G, h);

    for(Gm=0; Gm < nirreps; Gm++) {
      Gi = Gj = h^Gm;

      for(i=0; i < aoccpi[Gi]; i++) {
	    I = aocc_off[Gi] + i;
	    for(j=0; j < aoccpi[Gj]; j++) {
	      J = aocc_off[Gj] + j;
	      for(m=0; m < boccpi[Gm]; m++) {
	        M = bocc_off[Gm] + m;

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


  global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 2, 2, rho_params.Dij_lbl);
  global_dpd_->file2_mat_init(&D);
  global_dpd_->file2_mat_rd(&D);

  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 22, 22, 22, 22, 0, "GIjKl");
  for(h=0; h < nirreps; h++) {
    global_dpd_->buf4_mat_irrep_init(&G, h);
    global_dpd_->buf4_mat_irrep_rd(&G, h);

    for(Gm=0; Gm < nirreps; Gm++) {
      Gk = Gl = h^Gm;

      for(k=0; k < boccpi[Gk]; k++) {
	    K = bocc_off[Gk] + k;
	    for(l=0; l < boccpi[Gl]; l++) {
	      L = bocc_off[Gl] + l;
	      for(m=0; m < aoccpi[Gm]; m++) {
	        M = aocc_off[Gm] + m;

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

  if(params.debug_){
    global_dpd_->buf4_init(&Aints, PSIF_CC_AINTS, 0, 22, 22, 22, 22, 0, "A <Ij|Kl>");
    two_energy += global_dpd_->buf4_dot(&Aints, &G);
    global_dpd_->buf4_close(&Aints);

    total_two_energy += two_energy;
    outfile->Printf( "\tIJKL energy                = %20.15f\n", two_energy);
  }

  global_dpd_->buf4_close(&G);
  global_dpd_->file2_mat_close(&D);
  global_dpd_->file2_close(&D);

  global_dpd_->file2_init(&D1, PSIF_CC_OEI, 0, 0, 1, rho_params.DIA_lbl);
  global_dpd_->file2_mat_init(&D1);
  global_dpd_->file2_mat_rd(&D1);
  global_dpd_->file2_init(&D2, PSIF_CC_OEI, 0, 0, 1, rho_params.DAI_lbl);
  global_dpd_->file2_mat_init(&D2);
  global_dpd_->file2_mat_rd(&D2);

  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 0, 20, 2, 20, 0, "GIJKA");
  for(h=0; h < nirreps; h++) {
    global_dpd_->buf4_mat_irrep_init(&G, h);
    global_dpd_->buf4_mat_irrep_rd(&G, h);

    for(Gm=0; Gm < nirreps; Gm++) {
      Gi = Ga = h^Gm;

      for(i=0; i < aoccpi[Gi]; i++) {
	    I = aocc_off[Gi] + i;
	    for(a=0; a < avirtpi[Ga]; a++) {
	      A = avir_off[Ga] + a;
	      for(m=0; m < aoccpi[Gm]; m++) {
	        M = aocc_off[Gm] + m;

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

  if(params.debug_){
    two_energy = 0.0;
    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 0, 20, 2, 20, 0, "E <IJ||KA> (I>J,KA)");
    two_energy += global_dpd_->buf4_dot(&E, &G);
    global_dpd_->buf4_close(&E);
  }

  global_dpd_->buf4_close(&G);

  global_dpd_->file2_mat_close(&D1);
  global_dpd_->file2_close(&D1);
  global_dpd_->file2_mat_close(&D2);
  global_dpd_->file2_close(&D2);

  global_dpd_->file2_init(&D1, PSIF_CC_OEI, 0, 2, 3, rho_params.Dia_lbl);
  global_dpd_->file2_mat_init(&D1);
  global_dpd_->file2_mat_rd(&D1);
  global_dpd_->file2_init(&D2, PSIF_CC_OEI, 0, 2, 3, rho_params.Dai_lbl);
  global_dpd_->file2_mat_init(&D2);
  global_dpd_->file2_mat_rd(&D2);

  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 10, 30, 12, 30, 0, "Gijka");
  for(h=0; h < nirreps; h++) {
    global_dpd_->buf4_mat_irrep_init(&G, h);
    global_dpd_->buf4_mat_irrep_rd(&G, h);

    for(Gm=0; Gm < nirreps; Gm++) {
      Gi = Ga = h^Gm;

      for(i=0; i < boccpi[Gi]; i++) {
	    I = bocc_off[Gi] + i;
	    for(a=0; a < bvirtpi[Ga]; a++) {
	      A = bvir_off[Ga] + a;
	      for(m=0; m < boccpi[Gm]; m++) {
	        M = bocc_off[Gm] + m;

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
  if(params.debug_){
    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 10, 30, 12, 30, 0, "E <ij||ka> (i>j,ka)");
    two_energy += global_dpd_->buf4_dot(&E, &G);
    global_dpd_->buf4_close(&E);
  }
  global_dpd_->buf4_close(&G);

  global_dpd_->file2_mat_close(&D1);
  global_dpd_->file2_close(&D1);
  global_dpd_->file2_mat_close(&D2);
  global_dpd_->file2_close(&D2);

  global_dpd_->file2_init(&D1, PSIF_CC_OEI, 0, 0, 1, rho_params.DIA_lbl);
  global_dpd_->file2_mat_init(&D1);
  global_dpd_->file2_mat_rd(&D1);
  global_dpd_->file2_init(&D2, PSIF_CC_OEI, 0, 0, 1, rho_params.DAI_lbl);
  global_dpd_->file2_mat_init(&D2);
  global_dpd_->file2_mat_rd(&D2);

  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 23, 27, 23, 27, 0, "GiJkA");
  for(h=0; h < nirreps; h++) {
    global_dpd_->buf4_mat_irrep_init(&G, h);
    global_dpd_->buf4_mat_irrep_rd(&G, h);

    for(Gm=0; Gm < nirreps; Gm++) {
      Gi = Ga = h^Gm;

      for(i=0; i < aoccpi[Gi]; i++) {
	    I = aocc_off[Gi] + i;
	    for(a=0; a < avirtpi[Ga]; a++) {
	      A = avir_off[Ga] + a;
	      for(m=0; m < boccpi[Gm]; m++) {
	        M = bocc_off[Gm] + m;

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

  if(params.debug_){
    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 23, 27, 23, 27, 0, "E <iJ|kA>");
    two_energy += 2 * global_dpd_->buf4_dot(&E, &G);
    global_dpd_->buf4_close(&E);
  }

  global_dpd_->buf4_close(&G);

  global_dpd_->file2_mat_close(&D1);
  global_dpd_->file2_close(&D1);
  global_dpd_->file2_mat_close(&D2);
  global_dpd_->file2_close(&D2);


  global_dpd_->file2_init(&D1, PSIF_CC_OEI, 0, 2, 3, rho_params.Dia_lbl);
  global_dpd_->file2_mat_init(&D1);
  global_dpd_->file2_mat_rd(&D1);
  global_dpd_->file2_init(&D2, PSIF_CC_OEI, 0, 2, 3, rho_params.Dai_lbl);
  global_dpd_->file2_mat_init(&D2);
  global_dpd_->file2_mat_rd(&D2);

  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 22, 24, 22, 24, 0, "GIjKa");
  for(h=0; h < nirreps; h++) {
    global_dpd_->buf4_mat_irrep_init(&G, h);
    global_dpd_->buf4_mat_irrep_rd(&G, h);

    for(Gm=0; Gm < nirreps; Gm++) {
      Gi = Ga = h^Gm;

      for(i=0; i < boccpi[Gi]; i++) {
	    I = bocc_off[Gi] + i;
	    for(a=0; a < bvirtpi[Ga]; a++) {
	      A = bvir_off[Ga] + a;
	      for(m=0; m < aoccpi[Gm]; m++) {
	        M = aocc_off[Gm] + m;

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

  if(params.debug_){
    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 22, 24, 22, 24, 0, "E <Ij|Ka>");
    two_energy += 2 * global_dpd_->buf4_dot(&E, &G);
    global_dpd_->buf4_close(&E);
  }

  global_dpd_->buf4_close(&G);

  if(params.debug_){
    total_two_energy += two_energy;
    outfile->Printf( "\tIJKA energy                = %20.15f\n", two_energy);
  }

  global_dpd_->file2_mat_close(&D1);
  global_dpd_->file2_close(&D1);
  global_dpd_->file2_mat_close(&D2);
  global_dpd_->file2_close(&D2);

  if(params.debug_){
    two_energy = 0.0;
    global_dpd_->buf4_init(&DInts, PSIF_CC_DINTS, 0, 2, 7, 2, 7, 0, "D <IJ||AB> (I>J,A>B)");
    global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 2, 7, 2, 7, 0, "GIJAB");
    two_energy += global_dpd_->buf4_dot(&G, &DInts);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_close(&DInts);


    global_dpd_->buf4_init(&DInts, PSIF_CC_DINTS, 0, 12, 17, 12, 17, 0, "D <ij||ab> (i>j,a>b)");
    global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 12, 17, 12, 17, 0, "Gijab");
    two_energy += global_dpd_->buf4_dot(&G, &DInts);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_close(&DInts);

    global_dpd_->buf4_init(&DInts, PSIF_CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");
    global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 22, 28, 22, 28, 0, "GIjAb");
    two_energy += global_dpd_->buf4_dot(&G, &DInts);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_close(&DInts);

    two_energy *= 2;
    total_two_energy += two_energy;
    outfile->Printf( "\tIJAB energy                = %20.15f\n", two_energy);
  }


  global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 1, 1, rho_params.DAB_lbl);
  global_dpd_->file2_mat_init(&D);
  global_dpd_->file2_mat_rd(&D);

  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 20, 20, 20, 20, 0, "GIBJA");
  for(h=0; h < nirreps; h++) {
    global_dpd_->buf4_mat_irrep_init(&G, h);
    global_dpd_->buf4_mat_irrep_rd(&G, h);

    for(Gm=0; Gm < nirreps; Gm++) {
      Ga = Gb = h^Gm;

      for(b=0; b < avirtpi[Gb]; b++) {
	    B = avir_off[Gb] + b;
	    for(a=0; a < avirtpi[Ga]; a++) {
	      A = avir_off[Ga] + a;
	      for(m=0; m < aoccpi[Gm]; m++) {
	        M = aocc_off[Gm] + m;

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

  if(params.debug_){
    two_energy = 0.0;
    global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 20, 20, 20, 20, 0, "C <IA||JB>");
    two_energy += global_dpd_->buf4_dot(&C, &G);
    global_dpd_->buf4_close(&C);
  }

  global_dpd_->buf4_close(&G);

  global_dpd_->file2_mat_close(&D);
  global_dpd_->file2_close(&D);


  global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 3, 3, rho_params.Dab_lbl);
  global_dpd_->file2_mat_init(&D);
  global_dpd_->file2_mat_rd(&D);

  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 30, 30, 30, 30, 0, "Gibja");
  for(h=0; h < nirreps; h++) {
    global_dpd_->buf4_mat_irrep_init(&G, h);
    global_dpd_->buf4_mat_irrep_rd(&G, h);

    for(Gm=0; Gm < nirreps; Gm++) {
      Ga = Gb = h^Gm;

      for(b=0; b < bvirtpi[Gb]; b++) {
	    B = bvir_off[Gb] + b;
	    for(a=0; a < bvirtpi[Ga]; a++) {
	      A = bvir_off[Ga] + a;
	      for(m=0; m < boccpi[Gm]; m++) {
	        M = bocc_off[Gm] + m;

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

  if(params.debug_){
    global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 30, 30, 30, 30, 0, "C <ia||jb>");
    two_energy += global_dpd_->buf4_dot(&C, &G);
    global_dpd_->buf4_close(&C);
  }

  global_dpd_->buf4_close(&G);

  global_dpd_->file2_mat_close(&D);
  global_dpd_->file2_close(&D);


  global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 3, 3, rho_params.Dab_lbl);
  global_dpd_->file2_mat_init(&D);
  global_dpd_->file2_mat_rd(&D);

  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 24, 24, 24, 24, 0, "GIbJa");
  for(h=0; h < nirreps; h++) {
    global_dpd_->buf4_mat_irrep_init(&G, h);
    global_dpd_->buf4_mat_irrep_rd(&G, h);

    for(Gm=0; Gm < nirreps; Gm++) {
      Ga = Gb = h^Gm;

      for(b=0; b < bvirtpi[Gb]; b++) {
	    B = bvir_off[Gb] + b;
	    for(a=0; a < bvirtpi[Ga]; a++) {
	      A = bvir_off[Ga] + a;
	      for(m=0; m < aoccpi[Gm]; m++) {
	        M = aocc_off[Gm] + m;

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

  if(params.debug_){
    global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 24, 24, 24, 24, 0, "C <Ia|Jb>");
    two_energy += global_dpd_->buf4_dot(&C, &G);
    global_dpd_->buf4_close(&C);
  }

  global_dpd_->buf4_close(&G);

  global_dpd_->file2_mat_close(&D);
  global_dpd_->file2_close(&D);


  global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 1, 1, rho_params.DAB_lbl);
  global_dpd_->file2_mat_init(&D);
  global_dpd_->file2_mat_rd(&D);

  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 27, 27, 27, 27, 0, "GiBjA");
  for(h=0; h < nirreps; h++) {
    global_dpd_->buf4_mat_irrep_init(&G, h);
    global_dpd_->buf4_mat_irrep_rd(&G, h);

    for(Gm=0; Gm < nirreps; Gm++) {
      Ga = Gb = h^Gm;

      for(b=0; b < avirtpi[Gb]; b++) {
	    B = avir_off[Gb] + b;
	    for(a=0; a < avirtpi[Ga]; a++) {
	      A = avir_off[Ga] + a;
	      for(m=0; m < boccpi[Gm]; m++) {
	        M = bocc_off[Gm] + m;

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

  if(params.debug_){
    global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 27, 27, 27, 27, 0, "C <iA|jB>");
    two_energy += global_dpd_->buf4_dot(&C, &G);
    global_dpd_->buf4_close(&C);
  }

  global_dpd_->buf4_close(&G);

  global_dpd_->file2_mat_close(&D);
  global_dpd_->file2_close(&D);


  if(params.debug_){
    global_dpd_->buf4_init(&DInts, PSIF_CC_DINTS, 0, 24, 27, 24, 27, 0, "D <Ij|Ab> (Ib,jA)");
    global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 24, 27, 24, 27, 0, "GIbjA");
    two_energy -= global_dpd_->buf4_dot(&G, &DInts);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_close(&DInts);
  
    global_dpd_->buf4_init(&DInts, PSIF_CC_DINTS, 0, 27, 24, 27, 24, 0, "D <iJ|aB> (iB,Ja)");
    global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 27, 24, 27, 24, 0, "GiBJa");
    two_energy -= global_dpd_->buf4_dot(&G, &DInts);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_close(&DInts);
  
    total_two_energy += two_energy;
    outfile->Printf( "\tIBJA energy                = %20.15f\n", two_energy);
  
  
    two_energy = 0.0;
  
    global_dpd_->buf4_init(&FInts, PSIF_CC_FINTS, 0, 21, 7, 21, 5, 1, "F <AI|BC>");
    global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 21, 7, 21, 7, 0, "GCIAB");
    two_energy += global_dpd_->buf4_dot(&G, &FInts);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_close(&FInts);
  
    global_dpd_->buf4_init(&FInts, PSIF_CC_FINTS, 0, 31, 17, 31, 15, 1, "F <ai|bc>");
    global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 31, 17, 31, 17, 0, "Gciab");
    two_energy += global_dpd_->buf4_dot(&G, &FInts);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_close(&FInts);
  
    global_dpd_->buf4_init(&FInts, PSIF_CC_FINTS, 0, 25, 29, 25, 29, 0, "F <aI|bC>");
    global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 25, 29, 25, 29, 0, "GcIaB");
    two_energy += global_dpd_->buf4_dot(&G, &FInts);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_close(&FInts);
  
    global_dpd_->buf4_init(&FInts, PSIF_CC_FINTS, 0, 26, 28, 26, 28, 0, "F <Ai|Bc>");
    global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 26, 28, 26, 28, 0, "GCiAb");
    two_energy += global_dpd_->buf4_dot(&G, &FInts);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_close(&FInts);
  
    two_energy *= 2;
    total_two_energy += two_energy;
    outfile->Printf( "\tCIAB energy                = %20.15f\n", two_energy);
  
  
    two_energy = 0.0;
  
    global_dpd_->buf4_init(&BInts, PSIF_CC_BINTS, 0, 7, 7, 5, 5, 1, "B <AB|CD>");
    global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 7, 7, 7, 7, 0, "GABCD");
    two_energy += global_dpd_->buf4_dot(&G, &BInts);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_close(&BInts);
  
    global_dpd_->buf4_init(&BInts, PSIF_CC_BINTS, 0, 17, 17, 15, 15, 1, "B <ab|cd>");
    global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 17, 17, 17, 17, 0, "Gabcd");
    two_energy += global_dpd_->buf4_dot(&G, &BInts);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_close(&BInts);
  
    global_dpd_->buf4_init(&BInts, PSIF_CC_BINTS, 0, 28, 28, 28, 28, 0, "B <Ab|Cd>");
    global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 28, 28, 28, 28, 0, "GAbCd");
    two_energy += global_dpd_->buf4_dot(&G, &BInts);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_close(&BInts);
  
  
    total_two_energy += two_energy;
    outfile->Printf( "\tABCD energy                = %20.15f\n", two_energy);
  
    outfile->Printf( "\tTotal two-electron energy  = %20.15f\n", total_two_energy);
    outfile->Printf( "\t%7s correlation energy = %20.15f\n", params.wfn == "CCSD_T" ? "CCSD(T)" : params.wfn.c_str(),
  	  one_energy + total_two_energy);
    outfile->Printf( "\tTotal %7s energy       = %20.15f\n", params.wfn == "CCSD_T" ? "CCSD(T)" : params.wfn.c_str(),
  	  one_energy + total_two_energy + moinfo.eref);
    }
}

}} // namespace psi::ccdensity
