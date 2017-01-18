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
    \ingroup CCENERGY
    \brief Enter brief description of file here
*/
#include <cstdio>
#include <cstdlib>
#include "psi4/libdpd/dpd.h"
#include "psi4/psifiles.h"
#include "Params.h"
#include "MOInfo.h"
#include "ccwave.h"
#include "psi4/libciomr/libciomr.h"

namespace psi { namespace ccenergy {

/* pair_energies(): For RHF references, compute pair energies. Spin-adapt
** pair energies if SPINADAPT_ENERGIES is set to true.
**
** E(IJ) = T2(IJ,AB) * (<ij|ab> - <ij|ba>)
** E(Ij) = T2(Ij,Ab) * <ij|ab>
**
*/

void CCEnergyWavefunction::pair_energies(double** epair_aa, double** epair_ab)
{
  dpdbuf4 tau, D, E;

  if(params_.ref == 0) { /** RHF **/

    int i, j, ij;
    int irrep;
    int nocc_act = 0;
    int naa, nab;
    for(irrep=0; irrep<moinfo_.nirreps; irrep++)
      nocc_act += moinfo_.clsdpi[irrep];
    naa = nocc_act * (nocc_act-1)/2;
    nab = nocc_act * nocc_act;

    /* Compute alpha-alpha pair energies */
    if (naa) {
      double* eaa = init_array(naa);

      global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 2, 5, 0, 5, 1, "D <ij|ab>");
      global_dpd_->buf4_init(&tau, PSIF_CC_TAMPS, 0, 2, 5, 0, 5, 1, "tauIjAb");
      global_dpd_->buf4_init(&E, PSIF_CC_TMP0, 0, 2, 2, 2, 2, 0, "E <ij|kl>");
      global_dpd_->contract444(&D, &tau, &E, 0, 0, 1.0, 0.0);

      //dpd_buf4_print(&E, outfile, 1);

      /* Extract diagonal elements (i.e. pair energies) and print them out nicely */
      for(irrep=0; irrep<moinfo_.nirreps; irrep++) {
        double **block;
        dpdparams4 *Params = E.params;
        int p;
        int np = Params->rowtot[irrep];

        global_dpd_->buf4_mat_irrep_init(&E, irrep);
        global_dpd_->buf4_mat_irrep_rd(&E, irrep);
        block = E.matrix[irrep];

        for(p=0; p<np; p++) {
          int i, j, ij;

          i = Params->roworb[irrep][p][0];
          j = Params->roworb[irrep][p][1];

          ij = (i > j) ? i*(i-1)/2 + j : j*(j-1)/2 + i;
          eaa[ij] = block[p][p];
        }
	global_dpd_->buf4_mat_irrep_close(&E, irrep);
      }

      *epair_aa = eaa;

      global_dpd_->buf4_close(&tau);
      global_dpd_->buf4_close(&D);
      global_dpd_->buf4_close(&E);
    }

    /* Compute alpha-beta pair energies */
    if (nab) {
      double* eab = init_array(nab);

      global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
      global_dpd_->buf4_init(&tau, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauIjAb");
      global_dpd_->buf4_init(&E, PSIF_CC_TMP0, 0, 0, 0, 0, 0, 0, "E <ij|kl>");
      global_dpd_->contract444(&D, &tau, &E, 0, 0, 1.0, 0.0);

      //dpd_buf4_print(&E, outfile, 1);

      /* Extract diagonal elements (i.e. pair energies) and print them out nicely */
      for(irrep=0; irrep<moinfo_.nirreps; irrep++) {
        double **block;
        dpdparams4 *Params = E.params;
        int p;
        int np = Params->rowtot[irrep];

        global_dpd_->buf4_mat_irrep_init(&E, irrep);
        global_dpd_->buf4_mat_irrep_rd(&E, irrep);
        block = E.matrix[irrep];

        for(p=0; p<np; p++) {
          int i, j, ij;

          i = Params->roworb[irrep][p][0];
          j = Params->roworb[irrep][p][1];

          ij = i*nocc_act + j;
          eab[ij] = block[p][p];
        }
	global_dpd_->buf4_mat_irrep_close(&E, irrep);
      }

      *epair_ab = eab;

      global_dpd_->buf4_close(&tau);
      global_dpd_->buf4_close(&D);
      global_dpd_->buf4_close(&E);
    }
  }

}

void CCEnergyWavefunction::print_pair_energies(double* emp2_aa, double* emp2_ab, double* ecc_aa, double* ecc_ab)
{
  if(params_.ref == 0) { /** RHF **/

    int i, j, ij;
    int irrep;
    int nocc_act = 0;
    int naa, nab;
    for(irrep=0; irrep<moinfo_.nirreps; irrep++)
      nocc_act += moinfo_.clsdpi[irrep];
    naa = nocc_act * (nocc_act-1)/2;
    nab = nocc_act * nocc_act;

    if (!params_.spinadapt_energies) {
      double emp2_aa_tot = 0.0;
      double emp2_ab_tot = 0.0;
      double ecc_aa_tot = 0.0;
      double ecc_ab_tot = 0.0;

      outfile->Printf( "    Alpha-alpha pair energies\n");
      outfile->Printf( "        i       j         MP2             %s\n",params_.wfn.c_str());
      outfile->Printf( "      -----   -----   ------------   ------------\n");
      if (naa) {
        ij = 0;
        for(i=0; i<nocc_act; i++)
          for(j=0; j<i; j++,ij++) {
            outfile->Printf( "      %3d     %3d     %12.9lf   %12.9lf\n", i+1, j+1, emp2_aa[ij], ecc_aa[ij]);
            emp2_aa_tot += emp2_aa[ij];
            ecc_aa_tot += ecc_aa[ij];
          }
      }
      outfile->Printf( "      -------------   ------------   ------------\n");
      outfile->Printf( "          Total       %12.9lf   %12.9lf\n\n", emp2_aa_tot, ecc_aa_tot);

      outfile->Printf( "    Alpha-beta pair energies\n");
      outfile->Printf( "        i       j         MP2             %s\n",params_.wfn.c_str());
      outfile->Printf( "      -----   -----   ------------   ------------\n");
      if (nab) {
        ij = 0;
        for(i=0; i<nocc_act; i++)
          for(j=0; j<nocc_act; j++,ij++) {
            outfile->Printf( "      %3d     %3d     %12.9lf   %12.9lf\n", i+1, j+1, emp2_ab[ij], ecc_ab[ij]);
            emp2_ab_tot += emp2_ab[ij];
            ecc_ab_tot += ecc_ab[ij];
          }
      }
      outfile->Printf( "      -------------   ------------   ------------\n");
      outfile->Printf( "          Total       %12.9lf   %12.9lf\n", emp2_ab_tot, ecc_ab_tot);

    }
    else {
      /* Spin-adapt pair energies */
      double emp2_s_tot = 0.0;
      double emp2_t_tot = 0.0;
      double ecc_s_tot = 0.0;
      double ecc_t_tot = 0.0;

      outfile->Printf( "    Singlet pair energies\n");
      outfile->Printf( "        i       j         MP2             %s\n",params_.wfn.c_str());
      outfile->Printf( "      -----   -----   ------------   ------------\n");
      ij = 0;
      for(i=0; i<nocc_act; i++)
        for(j=0; j<=i; j++,ij++) {
          double emp2_s, ecc_s;
          int ij_ab = i*nocc_act + j;
          int ij_aa = i*(i-1)/2 + j;
          double eab, eaa;

          /* MP2 */
          eab = emp2_ab[ij_ab];
          if (i != j)
            eaa = emp2_aa[ij_aa];
          else
            eaa = 0.0;
          emp2_s = (i != j ? 2.0 : 1.0) * eab - 0.5 * eaa;

          /* CC */
          eab = ecc_ab[ij_ab];
          if (i != j)
            eaa = ecc_aa[ij_aa];
          else
            eaa = 0.0;
          ecc_s = (i != j ? 2.0 : 1.0) * eab - 0.5 * eaa;

          outfile->Printf( "      %3d     %3d     %12.9lf   %12.9lf\n", i+1, j+1, emp2_s, ecc_s);
          emp2_s_tot += emp2_s;
          ecc_s_tot += ecc_s;
        }
      outfile->Printf( "      -------------   ------------   ------------\n");
      outfile->Printf( "          Total       %12.9lf   %12.9lf\n\n", emp2_s_tot, ecc_s_tot);

      outfile->Printf( "    Triplet pair energies\n");
      outfile->Printf( "        i       j         MP2             %s\n",params_.wfn.c_str());
      outfile->Printf( "      -----   -----   ------------   ------------\n");
      if (naa) {
        ij = 0;
        for(i=0; i<nocc_act; i++)
          for(j=0; j<i; j++,ij++) {
            outfile->Printf( "      %3d     %3d     %12.9lf   %12.9lf\n", i+1, j+1, 1.5*emp2_aa[ij], 1.5*ecc_aa[ij]);
            emp2_t_tot += 1.5*emp2_aa[ij];
            ecc_t_tot += 1.5*ecc_aa[ij];
          }
      }
      outfile->Printf( "      -------------   ------------   ------------\n");
      outfile->Printf( "          Total       %12.9lf   %12.9lf\n", emp2_t_tot, ecc_t_tot);

    }

    outfile->Printf( "\n");
  }
}
}} // namespace psi::ccenergy
