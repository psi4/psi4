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
#include <stdio.h>
#include "psi4/libdpd/dpd.h"
#include "MOInfo.h"
#include "Params.h"
#include "Frozen.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccdensity {

    /* FOLD_RHF(): Fold the RHF Fock matrix contributions to the energy
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

    void fold_RHF(struct RHO_Params rho_params)
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
	outfile->Printf( "\n\tEnergies re-computed from Fock-adjusted CC density:\n");
	outfile->Printf(   "\t---------------------------------------------------\n");

	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 0, rho_params.DIJ_lbl);
	global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, 0, 0, "h(i,j)");
	this_energy = 2.0 * global_dpd_->file2_dot(&D, &F);
	global_dpd_->file2_close(&F);
	global_dpd_->file2_close(&D);
	one_energy += this_energy;

	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 1, 1, rho_params.DAB_lbl);
	global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, 1, 1, "h(a,b)");
	this_energy = 2.0 * global_dpd_->file2_dot(&D, &F);
	global_dpd_->file2_close(&F);
	global_dpd_->file2_close(&D);
	one_energy += this_energy;

	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 1, rho_params.DIA_lbl);
	global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, 0, 1, "h(i,a)");
	this_energy = 2.0 * global_dpd_->file2_dot(&D, &F);
	global_dpd_->file2_close(&F);
	global_dpd_->file2_close(&D);
	one_energy += this_energy;

	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 1, rho_params.DAI_lbl);
	global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, 0, 1, "h(i,a)");
	this_energy = 2.0 * global_dpd_->file2_dot(&D, &F);
	global_dpd_->file2_close(&F);
	global_dpd_->file2_close(&D);
	one_energy += this_energy;

	outfile->Printf( "\tOne-electron energy        = %20.15f\n", one_energy);

      }

      global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 0, rho_params.DIJ_lbl);
      global_dpd_->file2_mat_init(&D);
      global_dpd_->file2_mat_rd(&D);

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
		MI = G.params->rowidx[M][I];
		MJ = G.params->colidx[M][J];

		G.matrix[h][IM][JM] += D.matrix[Gi][i][j];
		G.matrix[h][MI][MJ] += D.matrix[Gi][i][j];
	      }
	    }
	  }
	}

	global_dpd_->buf4_mat_irrep_wrt(&G, h);
	global_dpd_->buf4_mat_irrep_close(&G, h);
      }

      /* Generate spin-adapted Gijkl just for the energy calculation */
      global_dpd_->buf4_scmcopy(&G, PSIF_CC_GAMMA, "2 Gijkl - Gijlk", 2);
      global_dpd_->buf4_sort_axpy(&G, PSIF_CC_GAMMA, pqsr, 0, 0, "2 Gijkl - Gijlk", -1);
      global_dpd_->buf4_close(&G);

      global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 0, 0, 0, 0, 0, "2 Gijkl - Gijlk");
      if(!params.aobasis) {
	global_dpd_->buf4_init(&Aints, PSIF_CC_AINTS, 0, 0, 0, 0, 0, 0, "A <ij|kl>");
	two_energy += global_dpd_->buf4_dot(&Aints, &G);
	global_dpd_->buf4_close(&Aints);
      }
      global_dpd_->buf4_close(&G);

      if(!params.aobasis) {
	total_two_energy += two_energy;
	outfile->Printf( "\tIJKL energy                = %20.15f\n", two_energy);

      }

      global_dpd_->file2_mat_close(&D);
      global_dpd_->file2_close(&D);

      global_dpd_->file2_init(&D1, PSIF_CC_OEI, 0, 0, 1, rho_params.DIA_lbl);
      global_dpd_->file2_mat_init(&D1);
      global_dpd_->file2_mat_rd(&D1);
      global_dpd_->file2_init(&D2, PSIF_CC_OEI, 0, 0, 1, rho_params.DAI_lbl);
      global_dpd_->file2_mat_init(&D2);
      global_dpd_->file2_mat_rd(&D2);

      global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 0, 10, 0, 10, 0, "GIjKa");
      for(h=0; h < nirreps; h++) {
	global_dpd_->buf4_mat_irrep_init(&G, h);
	global_dpd_->buf4_mat_irrep_rd(&G, h);

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

	global_dpd_->buf4_mat_irrep_wrt(&G, h);
	global_dpd_->buf4_mat_irrep_close(&G, h);
      }

      /* Generate spin-adapted Gijka just for the energy calculation */
      global_dpd_->buf4_scmcopy(&G, PSIF_CC_GAMMA, "2 Gijka - Gjika", 2);
      global_dpd_->buf4_sort_axpy(&G, PSIF_CC_GAMMA, qprs, 0, 10, "2 Gijka - Gjika", -1);
      global_dpd_->buf4_close(&G);

      global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 0, 10, 0, 10, 0, "2 Gijka - Gjika");
      if(!params.aobasis) {
	global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
	/* The factor of 4 here is necessary because Gijka is multiplied by 1/2 in Gijka.cc */
	two_energy = 4 * global_dpd_->buf4_dot(&E, &G);
	global_dpd_->buf4_close(&E);
      }
      global_dpd_->buf4_close(&G);

      if(!params.aobasis) {
	total_two_energy += two_energy;
	outfile->Printf( "\tIJKA energy                = %20.15f\n", two_energy);

      }

      global_dpd_->file2_mat_close(&D1);
      global_dpd_->file2_close(&D1);
      global_dpd_->file2_mat_close(&D2);
      global_dpd_->file2_close(&D2);

      if(!params.aobasis) {
	/* Generate spin-adapted Gijab jut for energy calculation */
	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 0, 5, 0, 5, 0, "GIjAb");
	global_dpd_->buf4_scmcopy(&G, PSIF_CC_GAMMA, "2 Gijab - Gijba", 2);
	global_dpd_->buf4_sort_axpy(&G, PSIF_CC_GAMMA, pqsr, 0, 5, "2 Gijab - Gijba", -1);
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 0, 5, 0, 5, 0, "2 Gijab - Gijba");
	global_dpd_->buf4_init(&DInts, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
	two_energy = 2 * global_dpd_->buf4_dot(&G, &DInts);
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_close(&DInts);

	total_two_energy += two_energy;
	outfile->Printf( "\tIJAB energy                = %20.15f\n", two_energy);

      }

      global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 1, 1, rho_params.DAB_lbl);
      global_dpd_->file2_mat_init(&D);
      global_dpd_->file2_mat_rd(&D);

      global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 10, 10, 10, 10, 0, "GIBJA");
      for(h=0; h < nirreps; h++) {
	global_dpd_->buf4_mat_irrep_init(&G, h);
	global_dpd_->buf4_mat_irrep_rd(&G, h);

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

	global_dpd_->buf4_mat_irrep_wrt(&G, h);
	global_dpd_->buf4_mat_irrep_close(&G, h);
      }

      global_dpd_->buf4_close(&G);
      global_dpd_->file2_mat_close(&D);
      global_dpd_->file2_close(&D);


      global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 1, 1, rho_params.DAB_lbl);
      global_dpd_->file2_mat_init(&D);
      global_dpd_->file2_mat_rd(&D);

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

      global_dpd_->buf4_close(&G);

      global_dpd_->file2_mat_close(&D);
      global_dpd_->file2_close(&D);

      if(!params.aobasis) {
	two_energy = 0.0;
	global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia||jb>");
	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 10, 10, 10, 10, 0, "GIBJA");
	two_energy += 2.0 * global_dpd_->buf4_dot(&G, &C);
	global_dpd_->buf4_close(&G);

	global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb>");
	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 10, 10, 10, 10, 0, "GIbJa");
	two_energy += 2.0 * global_dpd_->buf4_dot(&G, &C);
	global_dpd_->buf4_close(&G);

	global_dpd_->buf4_init(&DInts, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ib,ja)");
	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 10, 10, 10, 10, 0, "GIbjA");
	two_energy -= 2.0 * global_dpd_->buf4_dot(&G, &DInts);
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_close(&DInts);

	total_two_energy += two_energy;
	outfile->Printf( "\tIBJA energy                = %20.15f\n", two_energy);

      }

      if(!params.aobasis) {
	/* Generate spin-adapted Gciab just for energy calculation */
	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 11, 5, 11, 5, 0, "GCiAb");
	global_dpd_->buf4_scmcopy(&G, PSIF_CC_GAMMA, "2 Gciab - Gciba", 2);
	global_dpd_->buf4_sort_axpy(&G, PSIF_CC_GAMMA, pqsr, 11, 5, "2 Gciab - Gciba", -1);
	global_dpd_->buf4_close(&G);

	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 11, 5, 11, 5, 0, "2 Gciab - Gciba");
	global_dpd_->buf4_init(&FInts, PSIF_CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
	global_dpd_->buf4_sort(&FInts, PSIF_CC_FINTS, qpsr, 11, 5, "F <ai|bc>");
	global_dpd_->buf4_close(&FInts);
	global_dpd_->buf4_init(&FInts, PSIF_CC_FINTS, 0, 11, 5, 11, 5, 0, "F <ai|bc>");
	/* The factor of 4 here is necessary because Gciab is multiplied by 1/2 in Gciab.cc */
	two_energy = 4*global_dpd_->buf4_dot(&FInts, &G);
	global_dpd_->buf4_close(&FInts);
	global_dpd_->buf4_close(&G);
	total_two_energy += two_energy;
	outfile->Printf( "\tCIAB energy                = %20.15f\n", two_energy);

      }

      if(!params.aobasis) {
	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 5, 5, 5, 5, 0, "GAbCd");

	global_dpd_->buf4_scmcopy(&G, PSIF_CC_GAMMA, "2 Gabcd - Gabdc", 2);
	global_dpd_->buf4_sort_axpy(&G, PSIF_CC_GAMMA, pqsr, 5, 5, "2 Gabcd - Gabdc", -1);
	global_dpd_->buf4_close(&G);

	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 5, 5, 5, 5, 0, "2 Gabcd - Gabdc");
	global_dpd_->buf4_init(&BInts, PSIF_CC_BINTS, 0, 5, 5, 5, 5, 0, "B <ab|cd>");
	two_energy = global_dpd_->buf4_dot(&BInts, &G);
	global_dpd_->buf4_close(&BInts);
	global_dpd_->buf4_close(&G);
	total_two_energy += two_energy;
	outfile->Printf( "\tABCD energy                = %20.15f\n", two_energy);
      }

      if(!params.aobasis) {
	outfile->Printf( "\tTotal two-electron energy  = %20.15f\n", total_two_energy);
	if (params.ground) {
	  //outfile->Printf( "\tCCSD correlation energy    = %20.15f\n",
	  //	  one_energy + total_two_energy);
	  //outfile->Printf( "\tTotal CCSD energy          = %20.15f\n",
	  //	  one_energy + total_two_energy + moinfo.eref);
	  outfile->Printf( "\t%7s correlation energy = %20.15f\n", params.wfn == "CCSD_T" ? "CCSD(T)" : params.wfn.c_str(),
          	one_energy + total_two_energy);
          outfile->Printf( "\tTotal %7s energy       = %20.15f\n", params.wfn == "CCSD_T" ? "CCSD(T)" : params.wfn.c_str(),
          	one_energy + total_two_energy + moinfo.eref);
	}
	else {
	  outfile->Printf( "\tTotal EOM CCSD correlation energy        = %20.15f\n",
		  one_energy + total_two_energy);
	  outfile->Printf( "\tCCSD correlation + EOM excitation energy = %20.15f\n",
		  moinfo.ecc + params.cceom_energy);
	  outfile->Printf( "\tTotal EOM CCSD energy                    = %20.15f\n",
		  one_energy + total_two_energy + moinfo.eref);
	}
      }
    }

  }} // namespace psi::ccdensity
