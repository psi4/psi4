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
    \brief Calculates the one- and two-electron CC energies using the
    coresponding one- and two-particle density matrices.
*/
#include <stdio.h>
#include "psi4/libdpd/dpd.h"
#include "MOInfo.h"
#include "Params.h"
#include "Frozen.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccdensity {

    /* ENERGY_RHF(): Compute the RHF CC energy using the one- and two-particle
    ** density matrices.
    **
    */

    void energy_RHF(struct RHO_Params rho_params)
    {
      dpdfile2 D, F;
      dpdbuf4 G, A, B, C, DInts, E, FInts;
      double one_energy=0.0, two_energy=0.0, total_two_energy = 0.0;
      double this_energy;

      outfile->Printf( "\n\tEnergies re-computed from CC density:\n");
      outfile->Printf(   "\t-------------------------------------\n");

      global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 0, rho_params.DIJ_lbl);
      global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, 0, 0, "fIJ");
      this_energy = 2.0 * global_dpd_->file2_dot(&D, &F);
      global_dpd_->file2_close(&F);
      global_dpd_->file2_close(&D);
      one_energy += this_energy;

      global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 1, 1, rho_params.DAB_lbl);
      global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, 1, 1, "fAB");
      this_energy = 2.0 * global_dpd_->file2_dot(&D, &F);
      global_dpd_->file2_close(&F);
      global_dpd_->file2_close(&D);
      one_energy += this_energy;

      global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 1, rho_params.DIA_lbl);
      global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, 0, 1, "fIA");
      this_energy = 2.0 * global_dpd_->file2_dot(&D, &F);
      global_dpd_->file2_close(&F);
      global_dpd_->file2_close(&D);
      one_energy += this_energy;

      global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 1, rho_params.DAI_lbl);
      global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, 0, 1, "fIA");
      this_energy = 2.0 * global_dpd_->file2_dot(&D, &F);
      global_dpd_->file2_close(&F);
      global_dpd_->file2_close(&D);
      one_energy += this_energy;

      outfile->Printf( "\tOne-electron energy        = %20.15f\n", one_energy);


      total_two_energy = 0.0;


      /* E_ijkl = (2 Gijkl = Gijlk) <ij|kl> */
      global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 0, 0, 0, 0, 0, "GIjKl");
      global_dpd_->buf4_scmcopy(&G, PSIF_CC_GAMMA, "2 Gijkl - Gijlk", 2);
      global_dpd_->buf4_sort_axpy(&G, PSIF_CC_GAMMA, pqsr, 0, 0, "2 Gijkl - Gijlk", -1);
      global_dpd_->buf4_close(&G);
      global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 0, 0, 0, 0, 0, "2 Gijkl - Gijlk");
      global_dpd_->buf4_init(&A, PSIF_CC_AINTS, 0, 0, 0, 0, 0, 0, "A <ij|kl>");
      two_energy = global_dpd_->buf4_dot(&A, &G);
      global_dpd_->buf4_close(&A);
      global_dpd_->buf4_close(&G);
      total_two_energy += two_energy;
      outfile->Printf( "\tIJKL energy                = %20.15f\n", two_energy);


      two_energy = 0.0;

      /* E_ijka = 2 [ (2 Gijka - Gjika) <ij|ka> ] */
      /* NB: GIjKa is scaled by 1/2 in Gijka.cc. */
      global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 0, 10, 0, 10, 0, "GIjKa");
      global_dpd_->buf4_scmcopy(&G, PSIF_CC_GAMMA, "2 Gijka - Gjika", 2);
      global_dpd_->buf4_sort_axpy(&G, PSIF_CC_GAMMA, qprs, 0, 10, "2 Gijka - Gjika", -1);
      global_dpd_->buf4_close(&G);

      global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 0, 10, 0, 10, 0, "2 Gijka - Gjika");
      global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
      /* The factor of 4 here is necessary because Gijka is multiplied by 1/2 in Gijka.cc */
      two_energy = 4.0 * global_dpd_->buf4_dot(&E, &G);
      global_dpd_->buf4_close(&E);
      global_dpd_->buf4_close(&G);
      total_two_energy += two_energy;
      outfile->Printf( "\tIJKA energy                = %20.15f\n", two_energy);


      two_energy = 0.0;

      /* Generate spin-adapted Gijab jut for energy calculation */
      global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 0, 5, 0, 5, 0, "GIjAb");
      global_dpd_->buf4_scmcopy(&G, PSIF_CC_GAMMA, "2 Gijab - Gijba", 2);
      global_dpd_->buf4_sort_axpy(&G, PSIF_CC_GAMMA, pqsr, 0, 5, "2 Gijab - Gijba", -1);
      global_dpd_->buf4_close(&G);
      global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 0, 5, 0, 5, 0, "2 Gijab - Gijba");
      global_dpd_->buf4_init(&DInts, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
      two_energy = 2.0 * global_dpd_->buf4_dot(&G, &DInts);
      global_dpd_->buf4_close(&G);
      global_dpd_->buf4_close(&DInts);

      total_two_energy += two_energy;
      outfile->Printf( "\tIJAB energy                = %20.15f\n", two_energy);


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


      two_energy = 0.0;

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


      two_energy = 0.0;
      global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 5, 5, 5, 5, 0, "GAbCd");

      global_dpd_->buf4_scmcopy(&G, PSIF_CC_GAMMA, "2 Gabcd - Gabdc", 2);
      global_dpd_->buf4_sort_axpy(&G, PSIF_CC_GAMMA, pqsr, 5, 5, "2 Gabcd - Gabdc", -1);
      global_dpd_->buf4_close(&G);

      global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 5, 5, 5, 5, 0, "2 Gabcd - Gabdc");
      global_dpd_->buf4_init(&B, PSIF_CC_BINTS, 0, 5, 5, 5, 5, 0, "B <ab|cd>");
      two_energy = global_dpd_->buf4_dot(&B, &G);
      global_dpd_->buf4_close(&B);
      global_dpd_->buf4_close(&G);

      total_two_energy += two_energy;
      outfile->Printf( "\tABCD energy                = %20.15f\n", two_energy);

      outfile->Printf( "\tTotal two-electron energy  = %20.15f\n", total_two_energy);
      if (params.ground) {
	//outfile->Printf( "\tCCSD correlation energy    = %20.15f\n",
	//	one_energy + total_two_energy);
	//outfile->Printf( "\tTotal CCSD energy          = %20.15f\n",
	//	one_energy + total_two_energy + moinfo.eref);
	outfile->Printf( "\t%-7s correlation energy = %20.15f\n", params.wfn == "CCSD_T" ? "CCSD(T)" : params.wfn.c_str(),
                one_energy + total_two_energy);
        outfile->Printf( "\tTotal %-7s energy       = %20.15f\n", params.wfn == "CCSD_T" ? "CCSD(T)" : params.wfn.c_str(),
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

  }} // namespace psi::ccdensity
