/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2023 The Psi4 Developers.
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
    \ingroup CCENERGY
    \brief Print and compute pair energies.
*/
#include <cstdio>
#include <cstdlib>
#include "psi4/libdpd/dpd.h"
#include "psi4/psifiles.h"
#include "Params.h"
#include "MOInfo.h"
#include "psi4/cc/ccwave.h"
#include "psi4/libciomr/libciomr.h"
#include "psi4/libmints/matrix.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libpsi4util/PsiOutStream.h"

namespace psi {
namespace ccenergy {

/* pair_energies(): For RHF references, compute pair energies.
**
** E(IJ) = T2(IJ,AB) * (<ij|ab> - <ij|ba>)
** E(Ij) = T2(Ij,Ab) * <ij|ab>
**
*/

void CCEnergyWavefunction::pair_energies(std::vector<double>& epair_aa, std::vector<double>& epair_ab) const {
    dpdbuf4 tau, D, E;

    if (params_.ref == 0) { /** RHF **/

        auto nocc_act = moinfo_.clsdpi.sum();
        auto naa = nocc_act * (nocc_act - 1) / 2;
        auto nab = nocc_act * nocc_act;

        /* Compute alpha-alpha pair energies */
        if (naa) {
            global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 2, 5, 0, 5, 1, "D <ij|ab>");
            global_dpd_->buf4_init(&tau, PSIF_CC_TAMPS, 0, 2, 5, 0, 5, 1, "tauIjAb");
            global_dpd_->buf4_init(&E, PSIF_CC_TMP0, 0, 2, 2, 2, 2, 0, "E <ij|kl>");
            global_dpd_->contract444(&D, &tau, &E, 0, 0, 1.0, 0.0);

            // dpd_buf4_print(&E, outfile, 1);

            /* Extract diagonal elements (i.e. pair energies) and print them out nicely */
            for (int irrep = 0; irrep < moinfo_.nirreps; irrep++) {
                dpdparams4* Params = E.params;

                global_dpd_->buf4_mat_irrep_init(&E, irrep);
                global_dpd_->buf4_mat_irrep_rd(&E, irrep);
                auto block = E.matrix[irrep];

                for (int p = 0; p < Params->rowtot[irrep]; p++) {
                    auto i = Params->roworb[irrep][p][0];
                    auto j = Params->roworb[irrep][p][1];

                    auto ij = (i > j) ? i * (i - 1) / 2 + j : j * (j - 1) / 2 + i;
                    epair_aa[ij] = block[p][p];
                }
                global_dpd_->buf4_mat_irrep_close(&E, irrep);
            }

            global_dpd_->buf4_close(&tau);
            global_dpd_->buf4_close(&D);
            global_dpd_->buf4_close(&E);
        }

        /* Compute alpha-beta pair energies */
        if (nab) {
            auto eab = init_array(nab);

            global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
            global_dpd_->buf4_init(&tau, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauIjAb");
            global_dpd_->buf4_init(&E, PSIF_CC_TMP0, 0, 0, 0, 0, 0, 0, "E <ij|kl>");
            global_dpd_->contract444(&D, &tau, &E, 0, 0, 1.0, 0.0);

            // dpd_buf4_print(&E, outfile, 1);

            /* Extract diagonal elements (i.e. pair energies) and print them out nicely */
            for (int irrep = 0; irrep < moinfo_.nirreps; irrep++) {
                dpdparams4* Params = E.params;

                global_dpd_->buf4_mat_irrep_init(&E, irrep);
                global_dpd_->buf4_mat_irrep_rd(&E, irrep);
                auto block = E.matrix[irrep];

                for (int p = 0; p < Params->rowtot[irrep]; p++) {
                    auto i = Params->roworb[irrep][p][0];
                    auto j = Params->roworb[irrep][p][1];

                    auto ij = i * nocc_act + j;
                    epair_ab[ij] = block[p][p];
                }
                global_dpd_->buf4_mat_irrep_close(&E, irrep);
            }

            global_dpd_->buf4_close(&tau);
            global_dpd_->buf4_close(&D);
            global_dpd_->buf4_close(&E);
        }
    }
}

void CCEnergyWavefunction::print_pair_energies(const std::vector<double>& emp2_aa, const std::vector<double>& emp2_ab, const std::vector<double>& ecc_aa, const std::vector<double>& ecc_ab) {
    if (params_.ref != 0) { return; }

    // Be warned that this code is heavily boilerplate, but classes of pair energies differ in fussy indexing details.
    // Also, while it's possible to combine the storage and pairing loops, separating the tasks reads cleaner.

    // => Prep <=
    auto nocc_act = moinfo_.clsdpi.sum();
    double emp_tot, ecc_tot;
    std::string full_name;
    if (options_.get_str("WFN").find("CC2") != std::string::npos) {
        full_name = "CC2";
    } else if (options_.get_str("WFN").find("CC3") != std::string::npos) {
        full_name = "CC3";
    } else if (options_.get_str("WFN").find("CCSD") != std::string::npos) {
        full_name = "CCSD";
    } else {
        throw PSIEXCEPTION("Unknown wfn type in print_pair_energies.");
    }

    // => Alpha-alpha <=
    // ==> Store pair energies <==
    auto mp2_mat = std::make_shared<Matrix>("MP2 Alpha-Alpha Pair Energies", nocc_act, nocc_act);
    auto cc_mat = std::make_shared<Matrix>("CC Alpha-Alpha Pair Energies", nocc_act, nocc_act);
    for (int i = 0, ij = 0; i < nocc_act; i++) {
        for (int j = 0; j < i; j++, ij++) {
            mp2_mat->set(i, j, emp2_aa[ij]);
            mp2_mat->set(j, i, emp2_aa[ij]);
            cc_mat->set(i, j, ecc_aa[ij]);
            cc_mat->set(j, i, ecc_aa[ij]);
        }
    }
    set_array_variable("MP2 ALPHA-ALPHA PAIR ENERGIES", mp2_mat);
    set_array_variable("CC ALPHA-ALPHA PAIR ENERGIES", cc_mat);
    auto cc_str = full_name + " ALPHA-ALPHA PAIR ENERGIES";
    set_array_variable(cc_str, cc_mat);
    // Process::environment.globals["CCSD ALPHA-ALPHA PAIR ENERGIES"]
    // Process::environment.globals["CC2 ALPHA-ALPHA PAIR ENERGIES"]
    // Process::environment.globals["CC3 ALPHA-ALPHA PAIR ENERGIES"]
    // ==> Print pair energies <==
    emp_tot = 0;
    ecc_tot = 0;
    outfile->Printf("    Alpha-alpha pair energies\n");
    outfile->Printf("        i       j         MP2             %s\n", params_.wfn.c_str());
    outfile->Printf("      -----   -----   ------------   ------------\n");
    for (int i = 0, ij = 0; i < nocc_act; i++) {
        for (int j = 0; j < i; j++, ij++) {
            outfile->Printf("      %3d     %3d     %12.9lf   %12.9lf\n", i + 1, j + 1, mp2_mat->get(i, j), cc_mat->get(i, j));
            emp_tot += mp2_mat->get(i, j);
            ecc_tot += cc_mat->get(i, j);
        }
    }
    outfile->Printf("      -------------   ------------   ------------\n");
    outfile->Printf("          Total       %12.9lf   %12.9lf\n\n", emp_tot, ecc_tot);

    // => Alpha-beta <=
    // ==> Store pair energies <==
    mp2_mat = std::make_shared<Matrix>("MP2 Alpha-Beta Pair Energies", nocc_act, nocc_act);
    cc_mat = std::make_shared<Matrix>("CC Alpha-Beta Pair Energies", nocc_act, nocc_act);
    for (int i = 0, ij = 0; i < nocc_act; i++) {
        for (int j = 0; j < nocc_act; j++, ij++) {
            mp2_mat->set(i, j, emp2_ab[ij]);
            cc_mat->set(i, j, ecc_ab[ij]);
        }
    }
    set_array_variable("MP2 ALPHA-BETA PAIR ENERGIES", mp2_mat);
    set_array_variable("CC ALPHA-BETA PAIR ENERGIES", cc_mat);
    cc_str = full_name +  " ALPHA-BETA PAIR ENERGIES";
    set_array_variable(cc_str, cc_mat);
    // Process::environment.globals["CCSD ALPHA-BETA PAIR ENERGIES"]
    // Process::environment.globals["CC2 ALPHA-BETA PAIR ENERGIES"]
    // Process::environment.globals["CC3 ALPHA-BETA PAIR ENERGIES"]
    // ==> Print pair energies <==
    emp_tot = 0;
    ecc_tot = 0;
    outfile->Printf("    Alpha-beta pair energies\n");
    outfile->Printf("        i       j         MP2             %s\n", params_.wfn.c_str());
    outfile->Printf("      -----   -----   ------------   ------------\n");
    for (int i = 0, ij = 0; i < nocc_act; i++) {
        for (int j = 0; j < nocc_act; j++, ij++) {
            outfile->Printf("      %3d     %3d     %12.9lf   %12.9lf\n", i + 1, j + 1, mp2_mat->get(i, j), cc_mat->get(i, j));
            emp_tot += mp2_mat->get(i, j);
            ecc_tot += cc_mat->get(i, j);
        }
    }
    outfile->Printf("      -------------   ------------   ------------\n");
    outfile->Printf("          Total       %12.9lf   %12.9lf\n\n", emp_tot, ecc_tot);

    // => Singlet <=
    // ==> Store pair energies <==
    mp2_mat = std::make_shared<Matrix>("MP2 Singlet Pair Energies", nocc_act, nocc_act);
    cc_mat = std::make_shared<Matrix>("CC Singlet Pair Energies", nocc_act, nocc_act);
    for (int i = 0; i < nocc_act; i++) {
        for (int j = 0; j <= i; j++) {
            int ij_ab = i * nocc_act + j;
            int ij_aa = i * (i - 1) / 2 + j;;
            if (i == j) {
                mp2_mat->set(i, j, emp2_ab[ij_ab]);
                cc_mat->set(j, i, ecc_ab[ij_ab]);
            } else {
                auto emp2 = 2 * emp2_ab[ij_ab] - 0.5 * emp2_aa[ij_aa];
                auto ecc = 2 * ecc_ab[ij_ab] - 0.5 * ecc_aa[ij_aa];
                mp2_mat->set(i, j, emp2);
                mp2_mat->set(j, i, emp2);
                cc_mat->set(i, j, ecc);
                cc_mat->set(j, i, ecc);
            }
        }
    }
    set_array_variable("MP2 SINGLET PAIR ENERGIES", mp2_mat);
    set_array_variable("CC SINGLET PAIR ENERGIES", cc_mat);
    cc_str = full_name +  " SINGLET PAIR ENERGIES";
    set_array_variable(cc_str, cc_mat);
    // Process::environment.globals["CCSD SINGLET PAIR ENERGIES"]
    // Process::environment.globals["CC2 SINGLET PAIR ENERGIES"]
    // Process::environment.globals["CC3 SINGLET PAIR ENERGIES"]
    // ==> Print pair energies <==
    emp_tot = 0;
    ecc_tot = 0;
    outfile->Printf("    Singlet pair energies\n");
    outfile->Printf("        i       j         MP2             %s\n", params_.wfn.c_str());
    outfile->Printf("      -----   -----   ------------   ------------\n");
    for (int i = 0; i < nocc_act; i++) {
        for (int j = 0; j <= i; j++) {
            outfile->Printf("      %3d     %3d     %12.9lf   %12.9lf\n", i + 1, j + 1, mp2_mat->get(i, j), cc_mat->get(i, j));
            emp_tot += mp2_mat->get(i, j);
            ecc_tot += cc_mat->get(i, j);
        }
    }
    outfile->Printf("      -------------   ------------   ------------\n");
    outfile->Printf("          Total       %12.9lf   %12.9lf\n\n", emp_tot, ecc_tot);

    // => Triplet <=
    // ==> Store pair energies <==
    mp2_mat = std::make_shared<Matrix>("MP2 Triplet Pair Energies", nocc_act, nocc_act);
    cc_mat = std::make_shared<Matrix>("CC Triplet Pair Energies", nocc_act, nocc_act);
    for (int i = 0, ij = 0; i < nocc_act; i++) {
        for (int j = 0; j < i; j++, ij++) {
            mp2_mat->set(i, j, 1.5 * emp2_aa[ij]);
            mp2_mat->set(j, i, 1.5 * emp2_aa[ij]);
            cc_mat->set(i, j, 1.5 * ecc_aa[ij]);
            cc_mat->set(j, i, 1.5 * ecc_aa[ij]);
        }
    }
    set_array_variable("MP2 TRIPLET PAIR ENERGIES", mp2_mat);
    set_array_variable("CC TRIPLET PAIR ENERGIES", cc_mat);
    cc_str = full_name +  " TRIPLET PAIR ENERGIES";
    set_array_variable(cc_str, cc_mat);
    // Process::environment.globals["CCSD TRIPLET PAIR ENERGIES"]
    // Process::environment.globals["CC2 TRIPLET PAIR ENERGIES"]
    // Process::environment.globals["CC3 TRIPLET PAIR ENERGIES"]
    // ==> Print pair energies <==
    emp_tot = 0;
    ecc_tot = 0;
    outfile->Printf("    Triplet pair energies\n");
    outfile->Printf("        i       j         MP2             %s\n", params_.wfn.c_str());
    outfile->Printf("      -----   -----   ------------   ------------\n");
    for (int i = 0, ij = 0; i < nocc_act; i++) {
        for (int j = 0; j < i; j++, ij++) {
            outfile->Printf("      %3d     %3d     %12.9lf   %12.9lf\n", i + 1, j + 1, mp2_mat->get(i, j), cc_mat->get(i, j));
            emp_tot += mp2_mat->get(i, j);
            ecc_tot += cc_mat->get(i, j);
        }
    }
    outfile->Printf("      -------------   ------------   ------------\n");
    outfile->Printf("          Total       %12.9lf   %12.9lf\n\n", emp_tot, ecc_tot);
}
}  // namespace ccenergy
}  // namespace psi
