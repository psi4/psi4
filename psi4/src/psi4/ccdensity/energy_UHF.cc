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
#include <strings.h>
#include <string.h>
#include "psi4/libdpd/dpd.h"
#include "MOInfo.h"
#include "Params.h"
#include "Frozen.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccdensity {

    /* ENERGY_UHF(): Compute the UHF CC energy using the one- and two-particle
    ** density matrices.
    */

    void energy_UHF(struct RHO_Params rho_params)
    {
      dpdfile2 D, F;
      dpdbuf4 G, A, B, C, DInts, E, FInts;
      double one_energy=0.0, two_energy=0.0, total_two_energy = 0.0;
      double this_energy, test_energy;

      outfile->Printf( "\n\tEnergies re-computed from CC density:\n");
      outfile->Printf(   "\t-------------------------------------\n");

      global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 0, rho_params.DIJ_lbl);
      global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, 0, 0, "fIJ");
      this_energy = global_dpd_->file2_dot(&D, &F);
      global_dpd_->file2_close(&F);
      global_dpd_->file2_close(&D);

      /*    outfile->Printf( "\tDIJ = %20.15f\n", this_energy); */
      one_energy += this_energy;

      global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 2, 2, rho_params.Dij_lbl);
      global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, 2, 2, "fij");
      this_energy = global_dpd_->file2_dot(&D, &F);
      global_dpd_->file2_close(&F);
      global_dpd_->file2_close(&D);

      /*    outfile->Printf( "\tDij = %20.15f\n", this_energy); */
      one_energy += this_energy;

      global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 1, 1, rho_params.DAB_lbl);
      global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, 1, 1, "fAB");
      this_energy = global_dpd_->file2_dot(&D, &F);
      global_dpd_->file2_close(&F);
      global_dpd_->file2_close(&D);

      /*    outfile->Printf( "\tDAB = %20.15f\n", this_energy); */
      one_energy += this_energy;

      global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 3, 3, rho_params.Dab_lbl);
      global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, 3, 3, "fab");
      this_energy = global_dpd_->file2_dot(&D, &F);
      global_dpd_->file2_close(&F);
      global_dpd_->file2_close(&D);

      /*    outfile->Printf( "\tDab = %20.15f\n", this_energy); */
      one_energy += this_energy;

      global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 1, rho_params.DIA_lbl);
      global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, 0, 1, "fIA");
      this_energy = global_dpd_->file2_dot(&D, &F);
      global_dpd_->file2_close(&F);
      global_dpd_->file2_close(&D);

      /*    outfile->Printf( "\tDIA = %20.15f\n", this_energy); */
      one_energy += this_energy;

      global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 2, 3, rho_params.Dia_lbl);
      global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, 2, 3, "fia");
      this_energy = global_dpd_->file2_dot(&D, &F);
      global_dpd_->file2_close(&F);
      global_dpd_->file2_close(&D);

      /*    outfile->Printf( "\tDia = %20.15f\n", this_energy); */
      one_energy += this_energy;

      global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 1, rho_params.DAI_lbl);
      global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, 0, 1, "fIA");
      this_energy = global_dpd_->file2_dot(&D, &F);
      global_dpd_->file2_close(&F);
      global_dpd_->file2_close(&D);

      /*    outfile->Printf( "\tDAI = %20.15f\n", this_energy); */
      one_energy += this_energy;

      global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 2, 3, rho_params.Dai_lbl);
      global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, 2, 3, "fia");
      this_energy = global_dpd_->file2_dot(&D, &F);
      global_dpd_->file2_close(&F);
      global_dpd_->file2_close(&D);
      /*    outfile->Printf( "\tDai = %20.15f\n", this_energy); */
      one_energy += this_energy;

      outfile->Printf( "\tOne-electron energy        = %20.15f\n", one_energy);

      if (params.onepdm) return;

      total_two_energy = 0.0;


      two_energy = 0.0;

      global_dpd_->buf4_init(&A, PSIF_CC_AINTS, 0, 2, 2, 0, 0, 1, "A <IJ|KL>");
      global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 2, 2, 2, 2, 0, "GIJKL");
      two_energy += global_dpd_->buf4_dot(&G, &A);
      global_dpd_->buf4_close(&G);
      global_dpd_->buf4_close(&A);

      global_dpd_->buf4_init(&A, PSIF_CC_AINTS, 0, 12, 12, 10, 10, 1, "A <ij|kl>");
      global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 12, 12, 12, 12, 0, "Gijkl");
      two_energy += global_dpd_->buf4_dot(&G, &A);
      global_dpd_->buf4_close(&G);
      global_dpd_->buf4_close(&A);

      global_dpd_->buf4_init(&A, PSIF_CC_AINTS, 0, 22, 22, 22, 22, 0, "A <Ij|Kl>");
      global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 22, 22, 22, 22, 0, "GIjKl");
      two_energy += global_dpd_->buf4_dot(&G, &A);
      global_dpd_->buf4_close(&G);
      global_dpd_->buf4_close(&A);

      total_two_energy += two_energy;
      outfile->Printf( "\tIJKL energy                = %20.15f\n", two_energy);


      two_energy = 0.0;

      global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 2, 20, 2, 20, 0, "E <IJ||KA> (I>J,KA)");
      global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 2, 20, 2, 20, 0, "GIJKA");
      two_energy += global_dpd_->buf4_dot(&G, &E);
      global_dpd_->buf4_close(&G);
      global_dpd_->buf4_close(&E);

      global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 12, 30, 12, 30, 0, "E <ij||ka> (i>j,ka)");
      global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 12, 30, 12, 30, 0, "Gijka");
      two_energy += global_dpd_->buf4_dot(&G, &E);
      global_dpd_->buf4_close(&G);
      global_dpd_->buf4_close(&E);

      global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 22, 24, 22, 24, 0, "E <Ij|Ka>");
      global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 22, 24, 22, 24, 0, "GIjKa");
      two_energy += global_dpd_->buf4_dot(&G, &E);
      global_dpd_->buf4_close(&G);
      global_dpd_->buf4_close(&E);

      global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 23, 27, 23, 27, 0, "E <iJ|kA>");
      global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 23, 27, 23, 27, 0, "GiJkA");
      two_energy += global_dpd_->buf4_dot(&G, &E);
      global_dpd_->buf4_close(&G);
      global_dpd_->buf4_close(&E);

      two_energy *= 2;
      total_two_energy += two_energy;
      outfile->Printf( "\tIJKA energy                = %20.15f\n", two_energy);



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


      /*
      ** Compute the Gibja contribution to the two-electron energy.  By
      ** spin-case this contribution looks like:
      **
      **  E(AA) <-- sum_IBJA G(IB,JA) <JA||IB>
      **  E(BB) <-- sum_ibja G(ib,ja) <ja||ib>
      **  E(AB) <-- sum_IbJa ( G(Ib,Ja) <Ja|Ib> + G(iB,jA) <jA|iB> -
      **                         G(Ib,jA) <jA|bI> - G(iB,Ja) <Ja|Bi> )
      **
      **  See Gibja.c for the definition of G here.
      */

      two_energy = 0.0;

      global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 20, 20, 20, 20, 0, "C <IA||JB>");
      global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 20, 20, 20, 20, 0, "GIBJA");
      two_energy += global_dpd_->buf4_dot(&G, &C);
      global_dpd_->buf4_close(&G);
      global_dpd_->buf4_close(&C);

      global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 30, 30, 30, 30, 0, "C <ia||jb>");
      global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 30, 30, 30, 30, 0, "Gibja");
      two_energy += global_dpd_->buf4_dot(&G, &C);
      global_dpd_->buf4_close(&G);
      global_dpd_->buf4_close(&C);

      global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 24, 24, 24, 24, 0, "C <Ia|Jb>");
      global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 24, 24, 24, 24, 0, "GIbJa");
      two_energy += global_dpd_->buf4_dot(&G, &C);
      global_dpd_->buf4_close(&G);
      global_dpd_->buf4_close(&C);

      global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 27, 27, 27, 27, 0, "C <iA|jB>");
      global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 27, 27, 27, 27, 0, "GiBjA");
      two_energy += global_dpd_->buf4_dot(&G, &C);
      global_dpd_->buf4_close(&G);
      global_dpd_->buf4_close(&C);

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

      global_dpd_->buf4_init(&B, PSIF_CC_BINTS, 0, 7, 7, 5, 5, 1, "B <AB|CD>");
      global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 7, 7, 7, 7, 0, "GABCD");
      two_energy += global_dpd_->buf4_dot(&G, &B);
      global_dpd_->buf4_close(&G);
      global_dpd_->buf4_close(&B);

      global_dpd_->buf4_init(&B, PSIF_CC_BINTS, 0, 17, 17, 15, 15, 1, "B <ab|cd>");
      global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 17, 17, 17, 17, 0, "Gabcd");
      two_energy += global_dpd_->buf4_dot(&G, &B);
      global_dpd_->buf4_close(&G);
      global_dpd_->buf4_close(&B);

      global_dpd_->buf4_init(&B, PSIF_CC_BINTS, 0, 28, 28, 28, 28, 0, "B <Ab|Cd>");
      global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 28, 28, 28, 28, 0, "GAbCd");
      two_energy += global_dpd_->buf4_dot(&G, &B);
      global_dpd_->buf4_close(&G);
      global_dpd_->buf4_close(&B);

      total_two_energy += two_energy;
      outfile->Printf( "\tABCD energy                = %20.15f\n", two_energy);

      outfile->Printf( "\tTotal two-electron energy  = %20.15f\n", total_two_energy);
      if (params.ground) {
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
