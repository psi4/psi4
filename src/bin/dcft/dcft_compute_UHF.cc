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

#include "dcft.h"
#include <cmath>
#include <libdpd/dpd.h>
#include <libtrans/integraltransform.h>
#include <libdiis/diismanager.h>
#include <libpsio/psio.hpp>
#include <libpsio/psio.h>

#include "defines.h"

using namespace boost;

namespace psi{ namespace dcft{

/**
  * Compute DCFT energy using unrestricted HF reference
  */
double DCFTSolver::compute_energy_UHF()
{
    orbitalsDone_    = false;
    cumulantDone_ = false;
    densityConverged_ = false;
    energyConverged_ = false;
    // Perform SCF guess for the orbitals
    scf_guess();

    // If DCFT computation type is density fitting, build b(Q|mn) in AO basis
    if (options_.get_str("DCFT_TYPE") == "DF") df_build_b_ao();

    // Perform MP2 guess for the cumulant
    mp2_guess();

    // Print out information about the job
    outfile->Printf( "\n\tDCFT Functional:    \t\t %s", options_.get_str("DCFT_FUNCTIONAL").c_str());
    outfile->Printf( "\n\tAlgorithm:          \t\t %s", options_.get_str("ALGORITHM").c_str());
    outfile->Printf( "\n\tAO-Basis Integrals: \t\t %s", options_.get_str("AO_BASIS").c_str());
    if (options_.get_str("ALGORITHM") == "QC") {
        outfile->Printf( "\n\tQC type:            \t\t %s", options_.get_str("QC_TYPE").c_str());
        outfile->Printf( "\n\tQC coupling:        \t\t %s", options_.get_bool("QC_COUPLING") ? "TRUE" : "FALSE");
    }
    if (energy_level_shift_ > 1E-6) {
        outfile->Printf( "\n\tUsing level shift of %5.3f a.u.            ", energy_level_shift_);
    }

    // Things that are not implemented yet...
    if (options_.get_str("DERTYPE") == "FIRST" && (options_.get_str("DCFT_FUNCTIONAL") == "DC-12"))
        throw FeatureNotImplemented("DC-12 functional", "Analytic gradients", __FILE__, __LINE__);
    if (options_.get_str("AO_BASIS") == "DISK" && options_.get_str("DCFT_FUNCTIONAL") == "CEPA0")
        throw FeatureNotImplemented("CEPA0", "AO_BASIS = DISK", __FILE__, __LINE__);
    if (options_.get_str("AO_BASIS") == "DISK" && options_.get_str("ALGORITHM") == "QC" && options_.get_str("QC_TYPE") == "SIMULTANEOUS")
        throw FeatureNotImplemented("Simultaneous QC", "AO_BASIS = DISK", __FILE__, __LINE__);
    if (!(options_.get_str("ALGORITHM") == "TWOSTEP") && options_.get_str("DCFT_FUNCTIONAL") == "CEPA0")
        throw FeatureNotImplemented("CEPA0", "Requested DCFT algorithm", __FILE__, __LINE__);
    if (!(options_.get_str("DCFT_FUNCTIONAL") == "ODC-06" || options_.get_str("DCFT_FUNCTIONAL") == "ODC-12") && options_.get_str("DCFT_TYPE") == "DF")
        throw FeatureNotImplemented("DC-06/DC-12/ODC-13/CEPA0", "Density Fitting", __FILE__, __LINE__);
    if (options_.get_str("THREE_PARTICLE") == "PERTURBATIVE" && options_.get_str("DCFT_TYPE") == "DF")
        throw FeatureNotImplemented("Three-particle energy correction", "Density Fitting", __FILE__, __LINE__);

    // Orbital-optimized stuff
    if (options_.get_str("ALGORITHM") == "TWOSTEP" && orbital_optimized_)
        throw PSIEXCEPTION("Two-step algorithm cannot be run for the orbital-optimized DCFT methods");

    // Choose a paricular algorithm and solve the equations
    if(options_.get_str("ALGORITHM") == "TWOSTEP") {
        run_twostep_dcft();
    }
    else if (options_.get_str("ALGORITHM") == "SIMULTANEOUS") {
        if (!orbital_optimized_) {
            run_simult_dcft();
        }
        else {
            run_simult_dcft_oo();

            // Check spin-free equations
            if (0){
                // ************* Start *************
                // Build spinfree Lambda_IJAB = Lambda_IJAB + Lambda_ijab + 2 Lambda_IjAb
                dpdbuf4 Laa, Lab, Lbb, Lsf;
                global_dpd_->buf4_init(&Laa, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                              ID("[O>O]-"), ID("[V>V]-"), 0, "Lambda <OO|VV>");
                global_dpd_->buf4_copy(&Laa, PSIF_DCFT_DPD, "Lambda SF <OO|VV>");
                global_dpd_->buf4_close(&Laa);
                global_dpd_->buf4_init(&Lsf, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                                       ID("[O,O]"), ID("[V,V]"), 0, "Lambda SF <OO|VV>");
                global_dpd_->buf4_init(&Lab, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                              ID("[O,o]"), ID("[V,v]"), 0, "Lambda <Oo|Vv>");
                global_dpd_->buf4_init(&Lbb, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                              ID("[o>o]-"), ID("[v>v]-"), 0, "Lambda <oo|vv>");
                dpd_buf4_add(&Lsf, &Lab, 2.0);
                dpd_buf4_add(&Lsf, &Lbb, 1.0);
                global_dpd_->buf4_close(&Lab);
                global_dpd_->buf4_close(&Lbb);
                global_dpd_->buf4_close(&Lsf); // Passed!

                // According to Kutzelnigg's notes,
                // spinfree F_PQ = h_PQ + (g<PR|QS> - 1/2 g<RP|QS>) Gamma_RS
                // Following is to code F_IJ:
                if (0){
                    // Build spinfree Gamma_IJ = Gamma_IJ + Gamma_ij
                    //       spinfree Gamma_AB = Gamma_AB + Gamma_ab
                    dpdfile2 GIJ, Gij, GAB, Gab, GIJ_sf, GAB_sf;
                    global_dpd_->file2_init(&GIJ, PSIF_DCFT_DPD, 0, ID('O'), ID('O'), "Tau <O|O>");
                    global_dpd_->file2_copy(&GIJ, PSIF_DCFT_DPD, "Gamma SF <O|O>");
                    global_dpd_->file2_close(&GIJ);
                    global_dpd_->file2_init(&Gij, PSIF_DCFT_DPD, 0, ID('o'), ID('o'), "Tau <o|o>");
                    global_dpd_->file2_init(&GIJ_sf, PSIF_DCFT_DPD, 0, ID('O'), ID('O'), "Gamma SF <O|O>");
                    global_dpd_->file2_axpy(&Gij, &GIJ_sf, 1.0, 0);
                    global_dpd_->file2_close(&Gij);
                    global_dpd_->file2_close(&GIJ_sf);

                    global_dpd_->file2_init(&GIJ_sf, PSIF_DCFT_DPD, 0, ID('O'), ID('O'), "Gamma SF <O|O>");
                    global_dpd_->file2_mat_init(&GIJ_sf);
                    global_dpd_->file2_mat_rd(&GIJ_sf);
                    for(int h = 0; h < nirrep_; ++h){
                        for(int i = 0; i < naoccpi_[h]; ++i){
                                GIJ_sf.matrix[h][i][i] += 2.0;
                        }
                    }
                    global_dpd_->file2_mat_wrt(&GIJ_sf);
                    global_dpd_->file2_mat_close(&GIJ_sf);
                    global_dpd_->file2_close(&GIJ_sf);

                    global_dpd_->file2_init(&GIJ_sf, PSIF_DCFT_DPD, 0, ID('O'), ID('O'), "Gamma SF <O|O>");
                    global_dpd_->file2_print(&GIJ_sf, "outfile");
                    global_dpd_->file2_close(&GIJ_sf);

                    global_dpd_->file2_init(&GAB, PSIF_DCFT_DPD, 0, ID('V'), ID('V'), "Tau <V|V>");
                    global_dpd_->file2_copy(&GAB, PSIF_DCFT_DPD, "Gamma SF <V|V>");
                    global_dpd_->file2_close(&GAB);
                    global_dpd_->file2_init(&Gab, PSIF_DCFT_DPD, 0, ID('v'), ID('v'), "Tau <v|v>");
                    global_dpd_->file2_init(&GAB_sf, PSIF_DCFT_DPD, 0, ID('V'), ID('V'), "Gamma SF <V|V>");
                    global_dpd_->file2_axpy(&Gab, &GAB_sf, 1.0, 0);
                    global_dpd_->file2_close(&Gab);
                    global_dpd_->file2_close(&GAB_sf); // Passed!

                    psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);
                    // Build spinfree generalized fock F_IJ = h_IJ + (g<IP|JQ> - 1/2 g<PI|JQ>) Gamma_PQ
                    dpdfile2 hIJ, FIJ;
                    dpdbuf4 g1, g2, gtemp;
                    // Spinfree F_IJ = h_IJ
                    global_dpd_->file2_init(&hIJ, PSIF_LIBTRANS_DPD, 0, ID('O'), ID('O'), "H <O|O>");
                    global_dpd_->file2_copy(&hIJ, PSIF_LIBTRANS_DPD, "F SF <O|O>");
                    global_dpd_->file2_close(&hIJ);
                    // gtemp<IK|JL> = 1/2 (g<IK|JL> - g<KI|JL>) + 1/2 g<IK|JL>
                    global_dpd_->buf4_init(&g1, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[O,O]"),
                                           ID("[O,O]"), ID("[O,O]"), 1, "MO Ints <OO|OO>");
                    global_dpd_->buf4_copy(&g1, PSIF_LIBTRANS_DPD, "MO Ints Temp <OO|OO>");
                    global_dpd_->buf4_close(&g1);
                    global_dpd_->buf4_init(&gtemp, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[O,O]"),
                                           ID("[O,O]"), ID("[O,O]"), 0, "MO Ints Temp <OO|OO>");
                    global_dpd_->buf4_scm(&gtemp, 0.5);
                    global_dpd_->buf4_init(&g1, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[O,O]"),
                                           ID("[O,O]"), ID("[O,O]"), 0, "MO Ints <OO|OO>");
                    dpd_buf4_add(&gtemp, &g1, 0.5);
                    global_dpd_->buf4_close(&g1);
                    global_dpd_->buf4_close(&gtemp);
                    // gtemp<IK|JL> -> gtemp(IJ|KL)
                    global_dpd_->buf4_init(&gtemp, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[O,O]"),
                                           ID("[O,O]"), ID("[O,O]"), 0, "MO Ints Temp <OO|OO>");
                    global_dpd_->buf4_sort(&gtemp, PSIF_LIBTRANS_DPD, prqs, ID("[O,O]"), ID("[O,O]"), "MO Ints Temp (OO|OO)");
                    global_dpd_->buf4_close(&gtemp);

                    global_dpd_->buf4_init(&gtemp, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[O,O]"),
                                           ID("[O,O]"), ID("[O,O]"), 0, "MO Ints Temp (OO|OO)");
                    global_dpd_->file2_init(&GIJ_sf, PSIF_DCFT_DPD, 0, ID('O'), ID('O'), "Gamma SF <O|O>");
                    global_dpd_->file2_init(&FIJ, PSIF_LIBTRANS_DPD, 0, ID('O'), ID('O'), "F SF <O|O>");
                    // Spinfree F_IJ += gtemp(IJ|KL) Gamma_KL
                    global_dpd_->contract422(&gtemp, &GIJ_sf, &FIJ, 0, 0, 1.0, 1.0);
                    global_dpd_->buf4_close(&gtemp);
                    global_dpd_->file2_close(&GIJ_sf);
                    global_dpd_->file2_close(&FIJ);

                    // Spinfree F_IJ += (g <IA|JB> - 1/2 g<AI|JB>) Gamma_AB
                    // g <IA|JB> Gamma_AB
                    global_dpd_->file2_init(&FIJ, PSIF_LIBTRANS_DPD, 0, ID('O'), ID('O'), "F SF <O|O>");
                    global_dpd_->buf4_init(&g1, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                                           ID("[O>=O]+"), ID("[V>=V]+"), 0, "MO Ints (OO|VV)");
                    global_dpd_->file2_init(&GAB_sf, PSIF_DCFT_DPD, 0, ID('V'), ID('V'), "Gamma SF <V|V>");
                    global_dpd_->contract422(&g1, &GAB_sf, &FIJ, 0, 0, 1.0, 1.0);
                    global_dpd_->buf4_close(&g1);
                    // - 1/2 g<AI|JB> Gamma_AB
                    global_dpd_->buf4_init(&g2, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                                           ID("[O,O]"), ID("[V,V]"),0, "MO Ints <OO|VV>");
                    global_dpd_->contract422(&g2, &GAB_sf, &FIJ, 0, 1, -0.5, 1.0);
                    global_dpd_->buf4_close(&g2);
                    global_dpd_->file2_close(&GAB_sf);
                    global_dpd_->file2_close(&FIJ); // Passed! Actually, spinfree F_IJ = moFa_ = moFb_

                    psio_->close(PSIF_LIBTRANS_DPD, 1);
                }
                // The truth is, for closed-shell system,
                // F_PQ = moFa_ = moFb_,
                // so "F <O|O>" and "F <V|V>" work.

                // Build tilde_Lambda <IJ|AB> = 2 Lambda SF <IJ|AB> + Lambda SF <JI|AB>
                dpdbuf4 tL;
                global_dpd_->buf4_init(&Lsf, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                                       ID("[O,O]"), ID("[V,V]"), 0, "Lambda SF <OO|VV>");
                global_dpd_->buf4_copy(&Lsf, PSIF_DCFT_DPD, "Tilde Lambda <OO|VV>");
                global_dpd_->buf4_sort(&Lsf, PSIF_DCFT_DPD, qprs, ID("[O,O]"), ID("[V,V]"), "Lambda SF <JI|AB>");
                global_dpd_->buf4_close(&Lsf);
                global_dpd_->buf4_init(&tL, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                                       ID("[O,O]"), ID("[V,V]"), 0, "Tilde Lambda <OO|VV>");
                global_dpd_->buf4_scm(&tL, 2.0);
                global_dpd_->buf4_init(&Lsf, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                                       ID("[O,O]"), ID("[V,V]"), 0, "Lambda SF <JI|AB>");
                dpd_buf4_add(&tL, &Lsf, 1.0);
                global_dpd_->buf4_close(&Lsf);
                global_dpd_->buf4_close(&tL); // Passed!

                psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);
                // Build tilde_g <IJ|KL> = 2 g<IJ|KL> + g<JI|KL>
                dpdbuf4 tG, G;
                global_dpd_->buf4_init(&G, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[O,O]"),
                                       ID("[O,O]"), ID("[O,O]"), 0, "MO Ints <OO|OO>");
                global_dpd_->buf4_copy(&G, PSIF_LIBTRANS_DPD, "Tilde MO Ints <OO|OO>");
                global_dpd_->buf4_sort(&G, PSIF_LIBTRANS_DPD, qprs, ID("[O,O]"), ID("[O,O]"), "MO Ints <JI|KL>");
                global_dpd_->buf4_close(&G);
                global_dpd_->buf4_init(&tG, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[O,O]"),
                                       ID("[O,O]"), ID("[O,O]"), 0, "Tilde MO Ints <OO|OO>");
                global_dpd_->buf4_scm(&tG, 2.0);
                global_dpd_->buf4_init(&G, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[O,O]"),
                                       ID("[O,O]"), ID("[O,O]"), 0, "MO Ints <JI|KL>");
                dpd_buf4_add(&tG, &G, 1.0);
                global_dpd_->buf4_close(&G);
                global_dpd_->buf4_close(&tG); // Passed!

                // Build tilde_g <AB|CD> = 2 g<AB|CD> + g<BA|CD>
                global_dpd_->buf4_init(&G, PSIF_LIBTRANS_DPD, 0, ID("[V,V]"), ID("[V,V]"),
                                       ID("[V,V]"), ID("[V,V]"), 0, "MO Ints <VV|VV>");
                global_dpd_->buf4_copy(&G, PSIF_LIBTRANS_DPD, "Tilde MO Ints <VV|VV>");
                global_dpd_->buf4_sort(&G, PSIF_LIBTRANS_DPD, qprs, ID("[V,V]"), ID("[V,V]"), "MO Ints <BA|CD>");
                global_dpd_->buf4_close(&G);
                global_dpd_->buf4_init(&tG, PSIF_LIBTRANS_DPD, 0, ID("[V,V]"), ID("[V,V]"),
                                       ID("[V,V]"), ID("[V,V]"), 0, "Tilde MO Ints <VV|VV>");
                global_dpd_->buf4_scm(&tG, 2.0);
                global_dpd_->buf4_init(&G, PSIF_LIBTRANS_DPD, 0, ID("[V,V]"), ID("[V,V]"),
                                       ID("[V,V]"), ID("[V,V]"), 0, "MO Ints <BA|CD>");
                dpd_buf4_add(&tG, &G, 1.0);
                global_dpd_->buf4_close(&G);
                global_dpd_->buf4_close(&tG); // Passed!

                // Now let's do stationarity condition for cumulant
                // Build F intermediates
                dpdfile2 FOO, FVV;
                dpdbuf4 T, F;

                /*
                 * F_IJAB = 1/3 P+(IJ,AB) (f_CA tilde_Lambda_IJCB - f_IK tilde_Lambda_KJAB)
                 */

                global_dpd_->buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                                       ID("[O,O]"), ID("[V,V]"), 0, "Temp <OO|VV>");
                global_dpd_->buf4_init(&Lsf, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                                       ID("[O,O]"), ID("[V,V]"), 0, "Tilde Lambda <OO|VV>");
                // Temp_IJAB = f_AC tilde_Lambda<IJ|CB>
                global_dpd_->file2_init(&FVV, PSIF_LIBTRANS_DPD, 0, ID('V'), ID('V'), "F <V|V>");
                global_dpd_->contract244(&FVV, &Lsf, &T, 1, 2, 1, 1.0, 0.0);
                global_dpd_->file2_close(&FVV);
                global_dpd_->buf4_close(&Lsf);
                global_dpd_->buf4_close(&T);
                // F_IJAB = Temp_IJAB
                global_dpd_->buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                                       ID("[O,O]"), ID("[V,V]"), 0, "Temp <OO|VV>");
                global_dpd_->buf4_copy(&T, PSIF_DCFT_DPD, "F SF <OO|VV>");
                global_dpd_->buf4_sort(&T, PSIF_DCFT_DPD, qpsr, ID("[O,O]"), ID("[V,V]"), "P+(Temp) <OO|VV>");
                global_dpd_->buf4_close(&T);
                // F_IJAB += Temp_JIBA
                global_dpd_->buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                                       ID("[O,O]"), ID("[V,V]"), 0, "P+(Temp) <OO|VV>");
                global_dpd_->buf4_init(&F, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                                       ID("[O,O]"), ID("[V,V]"), 0, "F SF <OO|VV>");
                dpd_buf4_add(&F, &T, 1.0);
                global_dpd_->buf4_close(&T);
                global_dpd_->buf4_close(&F); // Passed!
                // Temp_IJAB = -f_IK tilde_Lambda<KJ|AB>
                global_dpd_->buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                                       ID("[O,O]"), ID("[V,V]"), 0, "Temp <OO|VV>");
                global_dpd_->file2_init(&FOO, PSIF_LIBTRANS_DPD, 0, ID('O'), ID('O'), "F <O|O>");
                global_dpd_->buf4_init(&Lsf, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                                       ID("[O,O]"), ID("[V,V]"), 0, "Tilde Lambda <OO|VV>");
                global_dpd_->contract244(&FOO, &Lsf, &T, 1, 0, 0, -1.0, 0.0);
                global_dpd_->file2_close(&FOO);
                global_dpd_->buf4_close(&T);
                // F_IJAB += Temp_IJAB
                global_dpd_->buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                                       ID("[O,O]"), ID("[V,V]"), 0, "Temp <OO|VV>");
                global_dpd_->buf4_init(&F, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                                       ID("[O,O]"), ID("[V,V]"), 0, "F SF <OO|VV>");
                dpd_buf4_add(&F, &T, 1.0);
                global_dpd_->buf4_close(&T);
                // F_IJAB += Temp_JIBA
                global_dpd_->buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                                       ID("[O,O]"), ID("[V,V]"), 0, "Temp <OO|VV>");
                global_dpd_->buf4_sort(&T, PSIF_DCFT_DPD, qpsr, ID("[O,O]"), ID("[V,V]"), "P+(Temp) <OO|VV>");
                global_dpd_->buf4_close(&T);
                global_dpd_->buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                                       ID("[O,O]"), ID("[V,V]"), 0, "P+(Temp) <OO|VV>");
                dpd_buf4_add(&F, &T, 1.0);
                // F_IJAB *= 1/3
                global_dpd_->buf4_scm(&F, 1.0/3.0);
                global_dpd_->buf4_close(&Lsf);
                global_dpd_->buf4_close(&F); // Passed!

                // Build G intermediates
                dpdbuf4 I;
                /*
                 * G_IJAB = 2 g_IJAB
                 */
                global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                                       ID("[O,O]"), ID("[V,V]"), 0, "MO Ints <OO|VV>");
                global_dpd_->buf4_copy(&I, PSIF_DCFT_DPD, "G SF <OO|VV>");
                global_dpd_->buf4_close(&I);
                global_dpd_->buf4_init(&G, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                                       ID("[O,O]"), ID("[V,V]"), 0, "G SF <OO|VV>");
                global_dpd_->buf4_scm(&G, 2.0);
                global_dpd_->buf4_close(&G); // Passed!

                /*
                 * G_IJAB += 1/3 tilde_g_IJKL spinfree_Lambda_KLAB
                 */
                global_dpd_->buf4_init(&G, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                                       ID("[O,O]"), ID("[V,V]"), 0, "G SF <OO|VV>");
                global_dpd_->buf4_init(&tG, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[O,O]"),
                                       ID("[O,O]"), ID("[O,O]"), 0, "Tilde MO Ints <OO|OO>");
                global_dpd_->buf4_init(&Lsf, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                                       ID("[O,O]"), ID("[V,V]"), 0, "Lambda SF <OO|VV>");
                global_dpd_->contract444(&tG, &Lsf, &G, 0, 1, 1.0/3.0, 1.0);
                global_dpd_->buf4_close(&Lsf);
                global_dpd_->buf4_close(&tG);
                global_dpd_->buf4_close(&G); // Passed!

                /*
                 * G_IJAB += 1/3 tilde_g_CDAB spinfree_Lambda_IJCD
                 */
                global_dpd_->buf4_init(&G, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                                       ID("[O,O]"), ID("[V,V]"), 0, "G SF <OO|VV>");
                global_dpd_->buf4_init(&tG, PSIF_LIBTRANS_DPD, 0, ID("[V,V]"), ID("[V,V]"),
                                       ID("[V,V]"), ID("[V,V]"), 0, "Tilde MO Ints <VV|VV>");
                global_dpd_->buf4_init(&Lsf, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                                       ID("[O,O]"), ID("[V,V]"), 0, "Lambda SF <OO|VV>");
                global_dpd_->contract444(&Lsf, &tG, &G, 0, 1, 1.0/3.0, 1.0);
                global_dpd_->buf4_close(&Lsf);
                global_dpd_->buf4_close(&tG);
                global_dpd_->buf4_close(&G);

                /*
                 * G_IJAB += 1/3 P+(IJ,AB) (3 g_ICAK spinfree_Lambda_KJCB
                 *                          - g_CJAK tilde_Lambda_KIBC
                 *                          - g_CIAK tilde_Lambda_KJCB)
                 */

                // g_ICAK -> g_IACK
                global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                                       ID("[O,O]"), ID("[V,V]"), 0, "MO Ints <OO|VV>");
                global_dpd_->buf4_sort(&I, PSIF_LIBTRANS_DPD, psrq, ID("[O,V]"), ID("[V,O]"), "MO Ints <OV|VO>");
                global_dpd_->buf4_close(&I);
                global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[V,O]"),
                                       ID("[O,V]"), ID("[V,O]"), 0, "MO Ints <OV|VO>");
                global_dpd_->buf4_sort(&I, PSIF_LIBTRANS_DPD, prqs, ID("[O,V]"), ID("[V,O]"), "MO Ints <IA|CK>");
                global_dpd_->buf4_close(&I);
                global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[V,O]"),
                                       ID("[O,V]"), ID("[V,O]"), 0, "MO Ints <IA|CK>");
                // spinfree_Lambda_KJCB -> spinfree_Lambda_JBCK
                global_dpd_->buf4_init(&Lsf, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                                       ID("[O,O]"), ID("[V,V]"), 0, "Lambda SF <OO|VV>");
                global_dpd_->buf4_sort(&Lsf, PSIF_DCFT_DPD, qsrp, ID("[O,V]"), ID("[V,O]"), "Lambda SF <JB|CK>");
                global_dpd_->buf4_close(&Lsf);
                global_dpd_->buf4_init(&Lsf, PSIF_DCFT_DPD, 0, ID("[O,V]"), ID("[V,O]"),
                                       ID("[O,V]"), ID("[V,O]"), 0, "Lambda SF <JB|CK>");
                // Temp_IAJB = 3 g_IACK spinfree_Lambda_JBCK
                global_dpd_->buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                                       ID("[O,V]"), ID("[O,V]"), 0, "Temp (OV|OV)");
                global_dpd_->contract444(&I, &Lsf, &T, 0, 0, 3.0, 0.0);
                global_dpd_->buf4_close(&T);
                global_dpd_->buf4_close(&Lsf);
                global_dpd_->buf4_close(&I);
                // Temp_IAJB -> Temp_IJAB
                global_dpd_->buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                                       ID("[O,V]"), ID("[O,V]"), 0, "Temp (OV|OV)");
                global_dpd_->buf4_sort(&T, PSIF_DCFT_DPD, prqs, ID("[O,O]"), ID("[V,V]"), "Temp <OO|VV>");
                global_dpd_->buf4_close(&T);
                // threeterms_IJAB = Temp_IJAB
                global_dpd_->buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                                       ID("[O,O]"), ID("[V,V]"), 0, "Temp <OO|VV>");
                global_dpd_->buf4_copy(&T, PSIF_DCFT_DPD, "threeterms <IJ|AB>");
                global_dpd_->buf4_close(&T);

                // g_CJAK -> g_CKJA
                global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                                       ID("[O>=O]+"), ID("[V>=V]+"), 0, "MO Ints (OO|VV)");
                global_dpd_->buf4_sort(&I, PSIF_LIBTRANS_DPD, rpsq, ID("[V,O]"), ID("[V,O]"), "MO Ints <VO|VO>");
                global_dpd_->buf4_close(&I);
                global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,O]"), ID("[V,O]"),
                                       ID("[V,O]"), ID("[V,O]"), 0, "MO Ints <VO|VO>");
                global_dpd_->buf4_sort(&I, PSIF_LIBTRANS_DPD, psqr, ID("[V,O]"), ID("[O,V]"), "MO Ints <CK|JA>");
                global_dpd_->buf4_close(&I);
                global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,O]"), ID("[O,V]"),
                                       ID("[V,O]"), ID("[O,V]"), 0, "MO Ints <CK|JA>");
                // tilde_Lambda_KIBC -> tilde_Lambda_CKIB
                global_dpd_->buf4_init(&Lsf, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                                       ID("[O,O]"), ID("[V,V]"), 0, "Tilde Lambda <OO|VV>");
                global_dpd_->buf4_sort(&Lsf, PSIF_DCFT_DPD, spqr, ID("[V,O]"), ID("[O,V]"), "Tilde Lambda <CK|IB>");
                global_dpd_->buf4_close(&Lsf);
                global_dpd_->buf4_init(&Lsf,PSIF_DCFT_DPD, 0, ID("[V,O]"), ID("[O,V]"),
                                       ID("[V,O]"), ID("[O,V]"), 0, "Tilde Lambda <CK|IB>");
                // Temp_JAIB = - g_CKJA tilde_Lambda_CKIB
                global_dpd_->buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                                       ID("[O,V]"), ID("[O,V]"), 0, "Temp (OV|OV)");
                global_dpd_->contract444(&I, &Lsf, &T, 1, 1, -1.0, 0.0);
                global_dpd_->buf4_close(&T);
                global_dpd_->buf4_close(&I);
                global_dpd_->buf4_close(&Lsf);
                // Temp_JAIB -> Temp_IJAB
                global_dpd_->buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                                       ID("[O,V]"), ID("[O,V]"), 0, "Temp (OV|OV)");
                global_dpd_->buf4_sort(&T, PSIF_DCFT_DPD, rpqs, ID("[O,O]"), ID("[V,V]"), "Temp <OO|VV>");
                global_dpd_->buf4_close(&T);
                // threeterms_IJAB += Temp_IJAB
                dpdbuf4 three;
                global_dpd_->buf4_init(&three, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                                       ID("[O,O]"), ID("[V,V]"), 0, "threeterms <IJ|AB>");
                global_dpd_->buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                                       ID("[O,O]"), ID("[V,V]"), 0, "Temp <OO|VV>");
                dpd_buf4_add(&three, &T, 1.0);
                global_dpd_->buf4_close(&T);
                global_dpd_->buf4_close(&three);

                // g_CIAK -> g_IACK
                global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,O]"), ID("[V,O]"),
                                       ID("[V,O]"), ID("[V,O]"), 0, "MO Ints <VO|VO>");
                global_dpd_->buf4_sort(&I, PSIF_LIBTRANS_DPD, qrps, ID("[O,V]"), ID("[V,O]"), "MO Ints <IA|CK>");
                global_dpd_->buf4_close(&I);
                global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[V,O]"),
                                       ID("[O,V]"), ID("[V,O]"), 0, "MO Ints <IA|CK>");
                // tilde_Lambda_KJCB -> tilde_Lambda_CKJB
                global_dpd_->buf4_init(&Lsf, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                                       ID("[O,O]"), ID("[V,V]"), 0, "Tilde Lambda <OO|VV>");
                global_dpd_->buf4_sort(&Lsf, PSIF_DCFT_DPD, rpqs, ID("[V,O]"), ID("[O,V]"), "Tilde Lambda <CK|JB>");
                global_dpd_->buf4_close(&Lsf);
                global_dpd_->buf4_init(&Lsf,PSIF_DCFT_DPD, 0, ID("[V,O]"), ID("[O,V]"),
                                       ID("[V,O]"), ID("[O,V]"), 0, "Tilde Lambda <CK|JB>");
                // Temp_IAJB = - g_IACK tilde_Lambda_CKJB
                global_dpd_->buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                                       ID("[O,V]"), ID("[O,V]"), 0, "Temp (OV|OV)");
                global_dpd_->contract444(&I, &Lsf, &T, 0, 1, -1.0, 0.0);
                global_dpd_->buf4_close(&T);
                global_dpd_->buf4_close(&I);
                global_dpd_->buf4_close(&Lsf);
                // Temp_IAJB -> Temp_IJAB
                global_dpd_->buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                                       ID("[O,V]"), ID("[O,V]"), 0, "Temp (OV|OV)");
                global_dpd_->buf4_sort(&T, PSIF_DCFT_DPD, prqs, ID("[O,O]"), ID("[V,V]"), "Temp <OO|VV>");
                global_dpd_->buf4_close(&T);
                // threeterms_IJAB += Temp_IJAB
                global_dpd_->buf4_init(&three, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                                       ID("[O,O]"), ID("[V,V]"), 0, "threeterms <IJ|AB>");
                global_dpd_->buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                                       ID("[O,O]"), ID("[V,V]"), 0, "Temp <OO|VV>");
                dpd_buf4_add(&three, &T, 1.0);
                global_dpd_->buf4_close(&T);
                global_dpd_->buf4_close(&three);

                // G_IJAB += 1/3 threeterms_IJAB
                global_dpd_->buf4_init(&G, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                                       ID("[O,O]"), ID("[V,V]"), 0, "G SF <OO|VV>");
                global_dpd_->buf4_init(&three, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                                       ID("[O,O]"), ID("[V,V]"), 0, "threeterms <IJ|AB>");
                dpd_buf4_add(&G, &three, 1.0/3.0);
                global_dpd_->buf4_close(&G);
                global_dpd_->buf4_close(&three);

                // G_IJAB += 1/3 threeterms_JIBA
                global_dpd_->buf4_init(&three, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                                       ID("[O,O]"), ID("[V,V]"), 0, "threeterms <IJ|AB>");
                global_dpd_->buf4_sort(&three, PSIF_DCFT_DPD, qpsr, ID("[O,O]"), ID("[V,V]"), "threeterms <JI|BA>");
                global_dpd_->buf4_close(&three);
                global_dpd_->buf4_init(&three, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                                       ID("[O,O]"), ID("[V,V]"), 0, "threeterms <JI|BA>");
                global_dpd_->buf4_init(&G, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                                       ID("[O,O]"), ID("[V,V]"), 0, "G SF <OO|VV>");
                dpd_buf4_add(&G, &three, 1.0/3.0);
                global_dpd_->buf4_close(&G);
                global_dpd_->buf4_close(&three);


                // Finally, let's add F_IJAB and G_IJAB and see if we get a zero
                dpdbuf4 sum;
                global_dpd_->buf4_init(&G, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                                       ID("[O,O]"), ID("[V,V]"), 0, "G SF <OO|VV>");
//                global_dpd_->buf4_print(&G, "outfile", 1);

                global_dpd_->buf4_copy(&G, PSIF_DCFT_DPD, "G+F <OO|VV>");
                global_dpd_->buf4_close(&G);
                global_dpd_->buf4_init(&sum, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                                       ID("[O,O]"), ID("[V,V]"), 0, "G+F <OO|VV>");
                global_dpd_->buf4_init(&F, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                                       ID("[O,O]"), ID("[V,V]"), 0, "F SF <OO|VV>");
//                global_dpd_->buf4_print(&F, "outfile", 1);

                dpd_buf4_add(&sum, &F, 1.0);
                global_dpd_->buf4_close(&F);

//                global_dpd_->buf4_print(&sum, "outfile", 1);
                global_dpd_->buf4_close(&sum);

                psio_->close(PSIF_LIBTRANS_DPD, 1);


                // ************* End *************

            }
        }
    }
    else if (options_.get_str("ALGORITHM") == "QC") {
        run_qc_dcft();
    }
    else {
        throw PSIEXCEPTION("Unknown DCFT algoritm");
    }

    // If not converged -> Break
    if(!orbitalsDone_ || !cumulantDone_ || !densityConverged_)
        throw ConvergenceError<int>("DCFT", maxiter_, cumulant_threshold_, cumulant_convergence_, __FILE__, __LINE__);

    outfile->Printf("\n\t*%6s SCF Energy                                 = %20.15f\n", options_.get_str("DCFT_FUNCTIONAL").c_str(), scf_energy_);
    outfile->Printf("\t*%6s Lambda Energy                              = %20.15f\n", options_.get_str("DCFT_FUNCTIONAL").c_str(), lambda_energy_);
    outfile->Printf("\t*%6s Total Energy                               = %20.15f\n", options_.get_str("DCFT_FUNCTIONAL").c_str(), new_total_energy_);


    Process::environment.globals["DCFT SCF ENERGY"]    = scf_energy_;
    Process::environment.globals["DCFT LAMBDA ENERGY"] = lambda_energy_;
    Process::environment.globals["DCFT TOTAL ENERGY"]  = new_total_energy_;

    // Compute three-particle contribution to the DCFT energy
    if (options_.get_str("THREE_PARTICLE") == "PERTURBATIVE") {
        // Check options
        if (options_.get_str("DERTYPE") == "FIRST")
            throw FeatureNotImplemented("DCFT three-particle energy correction", "Analytic gradients", __FILE__, __LINE__);
        // Compute the three-particle energy
        double three_particle_energy = compute_three_particle_energy();
        outfile->Printf("\t*DCFT Three-particle Energy                        = %20.15f\n", three_particle_energy);
        outfile->Printf("\t*DCFT Total Energy                                 = %20.15f\n", new_total_energy_ + three_particle_energy);
        // Set global variables
        Process::environment.globals["DCFT THREE-PARTICLE ENERGY"] = three_particle_energy;
        Process::environment.globals["CURRENT ENERGY"]             = new_total_energy_ + three_particle_energy;
    }
    else {
        Process::environment.globals["CURRENT ENERGY"]             = new_total_energy_;
    }

    if(!options_.get_bool("MO_RELAX")){
        outfile->Printf( "Warning!  The orbitals were not relaxed\n");
    }

    // Print natural occupations
    print_opdm();

    if (orbital_optimized_) {
        // Compute one-electron properties
        compute_oe_properties();
        // Write to MOLDEN file if requested
        if (options_.get_bool("MOLDEN_WRITE")) write_molden_file();
    }

    if(options_.get_bool("TPDM")) dump_density();
//    check_n_representability();

    if (options_.get_str("DCFT_FUNCTIONAL") == "CEPA0") {
        compute_unrelaxed_density_OOOO();
        compute_unrelaxed_density_OVOV();
        compute_unrelaxed_density_VVVV();
        compute_TPDM_trace();
    }

    return(new_total_energy_);
}

void
DCFTSolver::run_twostep_dcft()
{

    // This is the two-step update - in each macro iteration, update the orbitals first, then update lambda
    // to self-consistency, until converged.  When lambda is converged and only one scf cycle is needed to reach
    // the desired cutoff, we're done

    int cycle = 0;
    outfile->Printf( "\n\n\t*=================================================================================*\n"
                         "\t* Cycle  RMS [F, Kappa]   RMS Lambda Error   delta E        Total Energy     DIIS *\n"
                         "\t*---------------------------------------------------------------------------------*\n");

    // Set up the DIIS manager for the density cumulant and SCF iterations
    old_ca_->copy(Ca_);
    old_cb_->copy(Cb_);
    // Save F0 = H + G * Kappa for the Fock intermediate update in lambda iterations
    moF0a_->copy(Fa_);
    moF0b_->copy(Fb_);
    moF0a_->transform(Ca_);
    moF0b_->transform(Cb_);
    // Just so the correct value is printed in the first macro iteration
    orbitals_convergence_ = compute_scf_error_vector();
    // Start macro-iterations
    while((!orbitalsDone_ || !cumulantDone_) && cycle++ < maxiter_){
        outfile->Printf( "\t                          *** Macro Iteration %d ***\n"
                         "\tCumulant Iterations\n",cycle);
        // If it's the first iteration and the user requested to relax guess orbitals, then skip the density cumulant update
        if ((cycle != 1) || !options_.get_bool("RELAX_GUESS_ORBITALS")) {
            run_twostep_dcft_cumulant_updates();
        }
        else outfile->Printf( "\tSkipping the cumulant update to relax guess orbitals\n");
        // Break if it's a CEPA0 computation
        if (options_.get_str("DCFT_FUNCTIONAL") == "CEPA0") {
            orbitalsDone_ = true;
            cumulantDone_ = true;
            densityConverged_ = true;
            break;
        }
        // Build new Tau from the density cumulant in the MO basis and transform it the SO basis
        build_tau();
        // Compute tau exactly if requested
        if (exact_tau_) {
            refine_tau();
        }
        transform_tau();
        run_twostep_dcft_orbital_updates();
    }

    outfile->Printf( "\t*=================================================================================*\n");

}

int
DCFTSolver::run_twostep_dcft_cumulant_updates() {

    // Set up DIIS
    dpdbuf4 Laa, Lab, Lbb;
    global_dpd_->buf4_init(&Laa, PSIF_LIBTRANS_DPD, 0, ID("[O>O]-"), ID("[V>V]-"),
                  ID("[O>O]-"), ID("[V>V]-"), 0, "Lambda <OO|VV>");
    global_dpd_->buf4_init(&Lab, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "Lambda <Oo|Vv>");
    global_dpd_->buf4_init(&Lbb, PSIF_LIBTRANS_DPD, 0, ID("[o>o]-"), ID("[v>v]-"),
                  ID("[o>o]-"), ID("[v>v]-"), 0, "Lambda <oo|vv>");
    DIISManager lambdaDiisManager(maxdiis_, "DCFT DIIS Lambdas",DIISManager::LargestError,DIISManager::InCore);
    if ((nalpha_ + nbeta_) > 1) {
        lambdaDiisManager.set_error_vector_size(3, DIISEntry::DPDBuf4, &Laa,
                                                DIISEntry::DPDBuf4, &Lab,
                                                DIISEntry::DPDBuf4, &Lbb);
        lambdaDiisManager.set_vector_size(3, DIISEntry::DPDBuf4, &Laa,
                                          DIISEntry::DPDBuf4, &Lab,
                                          DIISEntry::DPDBuf4, &Lbb);
    }
    global_dpd_->buf4_close(&Laa);
    global_dpd_->buf4_close(&Lab);
    global_dpd_->buf4_close(&Lbb);
    cumulantDone_ = false;
    int nLambdaIterations = 0;
    // Start density cumulant (lambda) iterations
    while((!cumulantDone_ || !energyConverged_) && nLambdaIterations++ < maxiter_){
        std::string diisString;
        // Build new Tau from current Lambda
        if (options_.get_str("DCFT_FUNCTIONAL") != "CEPA0") {
            // If not CEPA0
            if (options_.get_bool("RELAX_TAU")) {
                build_tau();
                // Compute tau exactly if requested
                if (exact_tau_) {
                    refine_tau();
                }
                if (options_.get_str("AO_BASIS") == "DISK") {
                    // Transform new Tau to the SO basis
                    transform_tau();
                    // Build SO basis tensors for the <VV||VV>, <vv||vv>, and <Vv|Vv> terms in the G intermediate
                    build_AO_tensors();
                }
                else {
                    // Compute GTau contribution for the Fock operator
                    build_gtau();
                }
                // Update Fock operator for the F intermediate
                update_fock();
            }
            else {
                if (options_.get_str("AO_BASIS") == "DISK") {
                    // Build SO basis tensors for the <VV||VV>, <vv||vv>, and <Vv|Vv> terms in the G intermediate
                    build_AO_tensors();
                }
            }
        }
        // Build G and F intermediates needed for the density cumulant residual equations and DCFT energy computation
        build_cumulant_intermediates();
        // Compute the residuals for density cumulant equations
        cumulant_convergence_ = compute_cumulant_residual();
        // Update density cumulant tensor
        update_cumulant_jacobi();
        if(cumulant_convergence_ < diis_start_thresh_ && (nalpha_ + nbeta_) > 1){
            //Store the DIIS vectors
            dpdbuf4 Laa, Lab, Lbb, Raa, Rab, Rbb;
            global_dpd_->buf4_init(&Raa, PSIF_DCFT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"),
                          ID("[O>O]-"), ID("[V>V]-"), 0, "R <OO|VV>");
            global_dpd_->buf4_init(&Rab, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                          ID("[O,o]"), ID("[V,v]"), 0, "R <Oo|Vv>");
            global_dpd_->buf4_init(&Rbb, PSIF_DCFT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"),
                          ID("[o>o]-"), ID("[v>v]-"), 0, "R <oo|vv>");
            global_dpd_->buf4_init(&Laa, PSIF_DCFT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"),
                          ID("[O>O]-"), ID("[V>V]-"), 0, "Lambda <OO|VV>");
            global_dpd_->buf4_init(&Lab, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                          ID("[O,o]"), ID("[V,v]"), 0, "Lambda <Oo|Vv>");
            global_dpd_->buf4_init(&Lbb, PSIF_DCFT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"),
                          ID("[o>o]-"), ID("[v>v]-"), 0, "Lambda <oo|vv>");

            if(lambdaDiisManager.add_entry(6, &Raa, &Rab, &Rbb, &Laa, &Lab, &Lbb)){
                diisString += "S";
            }
            if(lambdaDiisManager.subspace_size() >= mindiisvecs_ && maxdiis_ > 0){
                diisString += "/E";
                lambdaDiisManager.extrapolate(3, &Laa, &Lab, &Lbb);
            }
            global_dpd_->buf4_close(&Raa);
            global_dpd_->buf4_close(&Rab);
            global_dpd_->buf4_close(&Rbb);
            global_dpd_->buf4_close(&Laa);
            global_dpd_->buf4_close(&Lab);
            global_dpd_->buf4_close(&Lbb);
        }
        // Save old DCFT energy
        old_total_energy_ = new_total_energy_;
        // Compute new DCFT energy (lambda contribution)
        if (options_.get_str("DCFT_FUNCTIONAL") == "CEPA0") {
            compute_cepa0_energy();
        } else {
            compute_dcft_energy();
        }
        new_total_energy_ = scf_energy_ + lambda_energy_;
        // Check convergence for density cumulant iterations
        cumulantDone_ = cumulant_convergence_ < cumulant_threshold_;
        energyConverged_ = fabs(new_total_energy_ - old_total_energy_) < cumulant_threshold_;
        if (options_.get_str("ALGORITHM") == "TWOSTEP") {
            outfile->Printf( "\t* %-3d   %12.3e      %12.3e   %12.3e  %21.15f  %-3s *\n",
                    nLambdaIterations, orbitals_convergence_, cumulant_convergence_, new_total_energy_ - old_total_energy_,
                    new_total_energy_, diisString.c_str());
        }
        if (fabs(cumulant_convergence_) > 100.0) throw PSIEXCEPTION("DCFT density cumulant equations diverged");

    }

    return nLambdaIterations;

}


void
DCFTSolver::run_twostep_dcft_orbital_updates() {

    SharedMatrix tmp = SharedMatrix(new Matrix("temp", nirrep_, nsopi_, nsopi_));

    // Set up DIIS
    DIISManager scfDiisManager(maxdiis_, "DCFT DIIS Orbitals",DIISManager::LargestError,DIISManager::InCore);
    if ((nalpha_ + nbeta_) > 1) {
        scfDiisManager.set_error_vector_size(2, DIISEntry::Matrix, scf_error_a_.get(),
                                             DIISEntry::Matrix, scf_error_b_.get());
        scfDiisManager.set_vector_size(2, DIISEntry::Matrix, Fa_.get(),
                                       DIISEntry::Matrix, Fb_.get());
    }
    // Update the orbitals
    int nSCFCycles = 0;
    // Reset the booleans that control the convergence
    densityConverged_ = false;
    energyConverged_ = false;
    outfile->Printf( "\tOrbital Updates\n");
    while((!densityConverged_ || !orbitalsDone_ || !energyConverged_) && (nSCFCycles++ < maxiter_)){
        std::string diisString;
        // Copy core hamiltonian into the Fock matrix array: F = H
        Fa_->copy(so_h_);
        Fb_->copy(so_h_);
        // Build the new Fock matrix from the SO integrals: F += Gbar * Kappa
        process_so_ints();
        // Save F0 = H + G * Kappa for the Fock intermediate update in lambda iterations
        moF0a_->copy(Fa_);
        moF0b_->copy(Fb_);
        moF0a_->transform(Ca_);
        moF0b_->transform(Cb_);
        // Save old SCF energy
        old_total_energy_ = new_total_energy_;
        // Add non-idempotent density contribution (Tau) to the Fock matrix: F += Gbar * Tau
        Fa_->add(g_tau_a_);
        Fb_->add(g_tau_b_);
        // Back up the SO basis Fock before it is symmetrically orthogonalized to transform it to the MO basis
        moFa_->copy(Fa_);
        moFb_->copy(Fb_);
        // Compute new SCF energy
        compute_scf_energy();
        // Check SCF convergence
        orbitals_convergence_ = compute_scf_error_vector();
        orbitalsDone_ = orbitals_convergence_ < orbitals_threshold_;
        if(orbitals_convergence_ < diis_start_thresh_ && (nalpha_ + nbeta_) > 1){
            if(scfDiisManager.add_entry(4, scf_error_a_.get(), scf_error_b_.get(), Fa_.get(), Fb_.get()))
                diisString += "S";
            if(scfDiisManager.subspace_size() > mindiisvecs_ && (nalpha_ + nbeta_) > 1){
                diisString += "/E";
                scfDiisManager.extrapolate(2, Fa_.get(), Fb_.get());
            }
        }
        // Transform the Fock matrix to the symmetrically orhogonalized basis set and digonalize it
        // Obtain new orbitals
        Fa_->transform(s_half_inv_);
        Fa_->diagonalize(tmp, epsilon_a_);
        old_ca_->copy(Ca_);
        Ca_->gemm(false, false, 1.0, s_half_inv_, tmp, 0.0);
        Fb_->transform(s_half_inv_);
        Fb_->diagonalize(tmp, epsilon_b_);
        old_cb_->copy(Cb_);
        Cb_->gemm(false, false, 1.0, s_half_inv_, tmp, 0.0);
        // Make sure that the orbital phase is retained
        correct_mo_phases(false);
        // Update SCF density (Kappa) and check its RMS
        densityConverged_ = update_scf_density() < orbitals_threshold_;
        // Compute the DCFT energy
        new_total_energy_ = scf_energy_ + lambda_energy_;
        // Check convergence of the total DCFT energy
        energyConverged_ = fabs(new_total_energy_ - old_total_energy_) < cumulant_threshold_;
        outfile->Printf( "\t* %-3d   %12.3e      %12.3e   %12.3e  %21.15f  %-3s *\n",
                nSCFCycles, orbitals_convergence_, cumulant_convergence_, new_total_energy_ - old_total_energy_,
                new_total_energy_, diisString.c_str());
        if (fabs(orbitals_convergence_) > 100.0) throw PSIEXCEPTION("DCFT orbital updates diverged");

    }
    // Write orbitals to the checkpoint file
    write_orbitals_to_checkpoint();
    orbitalsDone_ = nSCFCycles == 1;
    energyConverged_ = false;
    // Transform the Fock matrix to the MO basis
    moFa_->transform(Ca_);
    moFb_->transform(Cb_);
    // Transform two-electron integrals to the MO basis using new orbitals, build denominators
    transform_integrals();
}

void
DCFTSolver::run_simult_dcft()
{
    // This is the simultaneous orbital/lambda update algorithm
    int cycle = 0;
    outfile->Printf( "\n\n\t*=================================================================================*\n"
                         "\t* Cycle  RMS [F, Kappa]   RMS Lambda Error   delta E        Total Energy     DIIS *\n"
                         "\t*---------------------------------------------------------------------------------*\n");

    SharedMatrix tmp = SharedMatrix(new Matrix("temp", nirrep_, nsopi_, nsopi_));
    // Set up the DIIS manager
    DIISManager diisManager(maxdiis_, "DCFT DIIS vectors");
    dpdbuf4 Laa, Lab, Lbb;
    global_dpd_->buf4_init(&Laa, PSIF_LIBTRANS_DPD, 0, ID("[O>O]-"), ID("[V>V]-"),
                  ID("[O>O]-"), ID("[V>V]-"), 0, "Lambda <OO|VV>");
    global_dpd_->buf4_init(&Lab, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "Lambda <Oo|Vv>");
    global_dpd_->buf4_init(&Lbb, PSIF_LIBTRANS_DPD, 0, ID("[o>o]-"), ID("[v>v]-"),
                  ID("[o>o]-"), ID("[v>v]-"), 0, "Lambda <oo|vv>");
    diisManager.set_error_vector_size(5, DIISEntry::Matrix, scf_error_a_.get(),
                                         DIISEntry::Matrix, scf_error_b_.get(),
                                         DIISEntry::DPDBuf4, &Laa,
                                         DIISEntry::DPDBuf4, &Lab,
                                         DIISEntry::DPDBuf4, &Lbb);
    diisManager.set_vector_size(5, DIISEntry::Matrix, Fa_.get(),
                                   DIISEntry::Matrix, Fb_.get(),
                                   DIISEntry::DPDBuf4, &Laa,
                                   DIISEntry::DPDBuf4, &Lab,
                                   DIISEntry::DPDBuf4, &Lbb);
    global_dpd_->buf4_close(&Laa);
    global_dpd_->buf4_close(&Lab);
    global_dpd_->buf4_close(&Lbb);
    while((!orbitalsDone_ || !cumulantDone_ || !densityConverged_ || !energyConverged_)
            && cycle++ < maxiter_){
        std::string diisString;
        // Save the old energy
        old_total_energy_ = new_total_energy_;
        // Build new Tau from the density cumulant in the MO basis and transform it the SO basis
        build_tau();
        if (exact_tau_) {
            refine_tau();
        }
        transform_tau();
        // Copy core hamiltonian into the Fock matrix array: F = H
        Fa_->copy(so_h_);
        Fb_->copy(so_h_);
        // Build the new Fock matrix from the SO integrals: F += Gbar * Kappa
        process_so_ints();
        // Add non-idempotent density contribution (Tau) to the Fock matrix: F += Gbar * Tau
        Fa_->add(g_tau_a_);
        Fb_->add(g_tau_b_);
        // Back up the SO basis Fock before it is symmetrically orthogonalized to transform it to the MO basis
        moFa_->copy(Fa_);
        moFb_->copy(Fb_);
        // Transform the Fock matrix to the MO basis
        moFa_->transform(Ca_);
        moFb_->transform(Cb_);
        // Compute new SCF energy
        compute_scf_energy();
        // Add SCF energy contribution to the total DCFT energy
        new_total_energy_ = scf_energy_;
        // Check SCF convergence
        orbitals_convergence_ = compute_scf_error_vector();
        orbitalsDone_ = orbitals_convergence_ < orbitals_threshold_;
        // Build G and F intermediates needed for the density cumulant residual equations and DCFT energy computation
        build_cumulant_intermediates();
        // Compute the residuals for density cumulant equations
        cumulant_convergence_ = compute_cumulant_residual();
        if (fabs(cumulant_convergence_) > 100.0) throw PSIEXCEPTION("DCFT density cumulant equations diverged");
        // Check convergence for density cumulant iterations
        cumulantDone_ = cumulant_convergence_ < cumulant_threshold_;
        // Update density cumulant tensor
        update_cumulant_jacobi();
        // Compute new DCFT energy (lambda contribution)
        compute_dcft_energy();
        // Add lambda energy to the DCFT total energy
        new_total_energy_ += lambda_energy_;
        // Check convergence of the total DCFT energy
        energyConverged_ = fabs(old_total_energy_ - new_total_energy_) < cumulant_threshold_;
        if(orbitals_convergence_ < diis_start_thresh_ && cumulant_convergence_ < diis_start_thresh_){
            //Store the DIIS vectors
            dpdbuf4 Laa, Lab, Lbb, Raa, Rab, Rbb;
            global_dpd_->buf4_init(&Raa, PSIF_DCFT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"),
                          ID("[O>O]-"), ID("[V>V]-"), 0, "R <OO|VV>");
            global_dpd_->buf4_init(&Rab, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                          ID("[O,o]"), ID("[V,v]"), 0, "R <Oo|Vv>");
            global_dpd_->buf4_init(&Rbb, PSIF_DCFT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"),
                          ID("[o>o]-"), ID("[v>v]-"), 0, "R <oo|vv>");
            global_dpd_->buf4_init(&Laa, PSIF_DCFT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"),
                          ID("[O>O]-"), ID("[V>V]-"), 0, "Lambda <OO|VV>");
            global_dpd_->buf4_init(&Lab, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                          ID("[O,o]"), ID("[V,v]"), 0, "Lambda <Oo|Vv>");
            global_dpd_->buf4_init(&Lbb, PSIF_DCFT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"),
                          ID("[o>o]-"), ID("[v>v]-"), 0, "Lambda <oo|vv>");
            if(diisManager.add_entry(10, scf_error_a_.get(), scf_error_b_.get(), &Raa, &Rab, &Rbb,
                                       Fa_.get(), Fb_.get(), &Laa, &Lab, &Lbb)){
                diisString += "S";
            }
            if(diisManager.subspace_size() > mindiisvecs_){
                diisString += "/E";
                diisManager.extrapolate(5, Fa_.get(), Fb_.get(), &Laa, &Lab, &Lbb);
            }
            global_dpd_->buf4_close(&Raa);
            global_dpd_->buf4_close(&Rab);
            global_dpd_->buf4_close(&Rbb);
            global_dpd_->buf4_close(&Laa);
            global_dpd_->buf4_close(&Lab);
            global_dpd_->buf4_close(&Lbb);
        }

        // Transform the Fock matrix to the symmetrically orhogonalized basis set and digonalize it
        // Obtain new orbitals
        Fa_->transform(s_half_inv_);
        Fa_->diagonalize(tmp, epsilon_a_);
        old_ca_->copy(Ca_);
        Ca_->gemm(false, false, 1.0, s_half_inv_, tmp, 0.0);
        Fb_->transform(s_half_inv_);
        Fb_->diagonalize(tmp, epsilon_b_);
        old_cb_->copy(Cb_);
        Cb_->gemm(false, false, 1.0, s_half_inv_, tmp, 0.0);
        // Make sure that the orbital phase is retained
        if(!correct_mo_phases(false)){
            outfile->Printf("\t\tThere was a problem correcting the MO phases.\n"
                            "\t\tIf this does not converge, try ALGORITHM=TWOSTEP\n");
        }
        // Write orbitals to the checkpoint file
        write_orbitals_to_checkpoint();
        // Transform two-electron integrals to the MO basis using new orbitals, build denominators
        transform_integrals();
        // Update SCF density (Kappa) and check its RMS
        densityConverged_ = update_scf_density() < orbitals_threshold_;
        // If we've performed enough lambda updates since the last orbitals
        // update, reset the counter so another SCF update is performed
        outfile->Printf( "\t* %-3d   %12.3e      %12.3e   %12.3e  %21.15f  %-3s *\n",
                cycle, orbitals_convergence_, cumulant_convergence_, new_total_energy_ - old_total_energy_,
                new_total_energy_, diisString.c_str());

    }

    outfile->Printf( "\t*=================================================================================*\n");
}

}} // Namespace

