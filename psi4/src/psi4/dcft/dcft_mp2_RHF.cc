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

#include "dcft.h"
#include "defines.h"
#include "psi4/psifiles.h"
#include <vector>
#include "psi4/liboptions/liboptions.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libtrans/integraltransform.h"
#include "psi4/libtrans/mospace.h"
#include "psi4/libdpd/dpd.h"
#include "psi4/libdiis/diismanager.h"



namespace psi{ namespace dcft{
/**
 * Computes the Hartree-Fock energy and then the MP2 energy as an initial guess.
 * This code is responible for initializing the integral transformation too.
 */
void
DCFTSolver::mp2_guess_RHF()
{
    dcft_timer_on("DCFTSolver::mp2_guess()");

    // Initialize the integral transformation object
    std::vector<std::shared_ptr<MOSpace> > spaces;
    spaces.push_back(MOSpace::occ);
    spaces.push_back(MOSpace::vir);
    spaces.push_back(MOSpace::all);

    // This wavefunction is really the global reference wavefunction
    _ints = new IntegralTransform(shared_from_this(), spaces, IntegralTransform::Restricted);
    _ints->set_keep_iwl_so_ints(true);
    _ints->set_keep_dpd_so_ints(true);
    dpd_set_default(_ints->get_dpd_id());

    outfile->Printf( "\n\n\tTransforming two-electron integrals (transformation type: restricted)...\n");
    transform_integrals_RHF();

    std::string guess = options_.get_str("DCFT_GUESS");

    if (guess == "MP2") {
        outfile->Printf( "\tComputing MP2 amplitude guess...\n\n");

        psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

        dpdbuf4 I, D;

        /*
         * In spin-adapted closed-shell system, only alpha-beta case is needed for computing energy
         */

        // L_IjAb = <Ij|Ab> / D_IjAb
        dcft_timer_on("DCFTSolver::g_IJAB / D_IJAB");
        global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                      ID("[O,O]"), ID("[V,V]"), 0, "MO Ints <OO|VV>"); // MO Ints <Oo|Vv>
        global_dpd_->buf4_copy(&I, PSIF_DCFT_DPD, "Lambda SF <OO|VV>"); // Lambda <Oo|Vv>
        global_dpd_->buf4_close(&I);

        global_dpd_->buf4_init(&D, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                      ID("[O>=O]+"), ID("[V>=V]+"), 0, "D <OO|VV>"); // D <Oo|Vv>
        global_dpd_->buf4_init(&I, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                      ID("[O,O]"), ID("[V,V]"), 0, "Lambda SF <OO|VV>"); // Lambda <Oo|Vv>
        global_dpd_->buf4_dirprd(&D, &I);
        global_dpd_->buf4_close(&I);
        global_dpd_->buf4_close(&D);
        dcft_timer_off("DCFTSolver::g_IJAB / D_IJAB");

        /* build lambda <OO|VV> for tau and G intermediates */
        dpdbuf4 T;
        // Lambda_IJAB = Lambda_IjAb - Lambda_JiAb
        global_dpd_->buf4_init(&I, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                               ID("[O,O]"), ID("[V,V]"), 1, "Lambda SF <OO|VV>");
        global_dpd_->buf4_copy(&I, PSIF_DCFT_DPD, "Lambda <OO|VV>");
        // The purpose of having Lambda <oo|vv> is for better performance of DIIS
        global_dpd_->buf4_copy(&I, PSIF_DCFT_DPD, "Lambda <oo|vv>");
        global_dpd_->buf4_close(&I);

        /*
        * E = lambda_IjAb * M_IjAb
        * where M_IjAb = 2 * gbar_IjAb - gbar_JiAb
        */
        dpdbuf4 L, M, temp;

        dcft_timer_on("DCFTSolver::2 * g_IJAB - g_JIAB");
        // M_IjAb = g_IjAb - g_JiAb
        global_dpd_->buf4_init(&M, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                               ID("[O,O]"), ID("[V,V]"), 1, "MO Ints <OO|VV>");
        global_dpd_->buf4_copy(&M, PSIF_LIBTRANS_DPD, "MO Ints Temp <OO|VV>");
        global_dpd_->buf4_close(&M);
        global_dpd_->buf4_init(&M, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                               ID("[O,O]"), ID("[V,V]"), 0, "MO Ints Temp <OO|VV>");
        global_dpd_->buf4_init(&temp, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                               ID("[O,O]"), ID("[V,V]"), 0, "MO Ints <OO|VV>");
        // M_IjAb += g_IjAb
        dpd_buf4_add(&M, &temp, 1.0);
        global_dpd_->buf4_close(&temp);
        dcft_timer_off("DCFTSolver::2 * g_IJAB - g_JIAB");

        dcft_timer_on("DCFTSolver::lambda_IjAb M_IjAb");
        global_dpd_->buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                               ID("[O,O]"), ID("[V,V]"), 0, "Lambda SF <OO|VV>");

        // E_MP2 = lambda_IjAb * M_IjAb
        double eMP2 = global_dpd_->buf4_dot(&L, &M);
        global_dpd_->buf4_close(&M);
        global_dpd_->buf4_close(&L);
        dcft_timer_off("DCFTSolver::lambda_IjAb M_IjAb");

        new_total_energy_ = scf_energy_ + eMP2;
        outfile->Printf( "\t*Total Hartree-Fock energy        = %20.15f\n", scf_energy_);
        outfile->Printf( "\t Total MP2 correlation energy     = %20.15f\n", eMP2);
        outfile->Printf( "\t*Total MP2 energy                 = %20.15f\n", new_total_energy_);

        Process::environment.globals["MP2 TOTAL ENERGY"] = new_total_energy_;
        Process::environment.globals["MP2 CORRELATION ENERGY"] = eMP2;

        psio_->close(PSIF_LIBTRANS_DPD, 1);

    }
    // Not implemented
    else if(guess == "CC" || guess == "BCC"){
        throw FeatureNotImplemented("Spin-adapted RHF-reference ODC-12", "DCFT_GUESS=CC/BCC", __FILE__, __LINE__);

    }

    dcft_timer_off("DCFTSolver::mp2_guess()");


}
}} // Namespace
