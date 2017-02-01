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
DCFTSolver::mp2_guess()
{
    dcft_timer_on("DCFTSolver::mp2_guess()");

    dpdbuf4 I, D;

    // Initialize the integral transformation object
    std::vector<std::shared_ptr<MOSpace> > spaces;
    spaces.push_back(MOSpace::occ);
    spaces.push_back(MOSpace::vir);
    _ints = new IntegralTransform(shared_from_this(), spaces, IntegralTransform::Unrestricted);
    _ints->set_keep_iwl_so_ints(true);
    _ints->set_keep_dpd_so_ints(true);
    dpd_set_default(_ints->get_dpd_id());

    outfile->Printf( "\n\n\tTransforming two-electron integrals (transformation type: unrestricted)...\n");
    transform_integrals();

    std::string guess = options_.get_str("DCFT_GUESS");

    if (guess == "MP2") {
        outfile->Printf( "\tComputing MP2 amplitude guess...\n\n");

        psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

        /*
        * L_ijab = <ij||ab> / D_ijab
        */

        // L_IJAB = <IJ||AB> / D_IJAB
        global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O>O]-"), ID("[V>V]-"),
                      ID("[O,O]"), ID("[V,V]"), 1, "MO Ints <OO|VV>");
        global_dpd_->buf4_copy(&I, PSIF_DCFT_DPD, "Lambda <OO|VV>");
        global_dpd_->buf4_close(&I);
        global_dpd_->buf4_init(&I, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                      ID("[O>O]-"), ID("[V>V]-"), 0, "Lambda <OO|VV>");
        global_dpd_->buf4_init(&D, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                      ID("[O>=O]+"), ID("[V>=V]+"), 0, "D <OO|VV>");
        global_dpd_->buf4_dirprd(&D, &I);
        global_dpd_->buf4_close(&I);
        global_dpd_->buf4_close(&D);

        // L_IjAb = <Ij|Ab> / D_IjAb
        global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                      ID("[O,o]"), ID("[V,v]"), 0, "MO Ints <Oo|Vv>");
        global_dpd_->buf4_copy(&I, PSIF_DCFT_DPD, "Lambda <Oo|Vv>");
        global_dpd_->buf4_close(&I);
        global_dpd_->buf4_init(&I, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                      ID("[O,o]"), ID("[V,v]"), 0, "Lambda <Oo|Vv>");
        global_dpd_->buf4_init(&D, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                      ID("[O,o]"), ID("[V,v]"), 0, "D <Oo|Vv>");
        global_dpd_->buf4_dirprd(&D, &I);
        global_dpd_->buf4_close(&I);
        global_dpd_->buf4_close(&D);

        // L_ijab = <ij||ab> / D_ijab
        global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o>o]-"), ID("[v>v]-"),
                      ID("[o,o]"), ID("[v,v]"), 1, "MO Ints <oo|vv>");
        global_dpd_->buf4_copy(&I, PSIF_DCFT_DPD, "Lambda <oo|vv>");
        global_dpd_->buf4_close(&I);
        global_dpd_->buf4_init(&I, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                      ID("[o>o]-"), ID("[v>v]-"), 0, "Lambda <oo|vv>");
        global_dpd_->buf4_init(&D, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                      ID("[o>=o]+"), ID("[v>=v]+"), 0, "D <oo|vv>");
        global_dpd_->buf4_dirprd(&D, &I);
        global_dpd_->buf4_close(&I);
        global_dpd_->buf4_close(&D);


        /*
     * E = 1/4 L_IJAB <IJ||AB>
     *        +L_IjAb <Ij|Ab>
     *    +1/4 L_ijab <ij||ab>
     */
        dpdbuf4 L;
        // Alpha - Alpha
        global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                      ID("[O,O]"), ID("[V,V]"), 1, "MO Ints <OO|VV>");
        global_dpd_->buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                      ID("[O>O]-"), ID("[V>V]-"), 0, "Lambda <OO|VV>");
        double eAA = 0.25 * global_dpd_->buf4_dot(&L, &I);
        global_dpd_->buf4_close(&I);
        global_dpd_->buf4_close(&L);

        // Alpha - Beta
        global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                      ID("[O,o]"), ID("[V,v]"), 0, "MO Ints <Oo|Vv>");
        global_dpd_->buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                      ID("[O,o]"), ID("[V,v]"), 0, "Lambda <Oo|Vv>");
        double eAB = global_dpd_->buf4_dot(&L, &I);
        global_dpd_->buf4_close(&I);
        global_dpd_->buf4_close(&L);

        // Beta - Beta
        global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                      ID("[o,o]"), ID("[v,v]"), 1, "MO Ints <oo|vv>");
        global_dpd_->buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                      ID("[o>o]-"), ID("[v>v]-"), 0, "Lambda <oo|vv>");
        double eBB = 0.25 * global_dpd_->buf4_dot(&L, &I);
        global_dpd_->buf4_close(&I);
        global_dpd_->buf4_close(&L);

        new_total_energy_ = scf_energy_ + eAA + eAB + eBB;
        outfile->Printf( "\t*Total Hartree-Fock energy        = %20.15f\n", scf_energy_);
        outfile->Printf( "\t Alpha - Alpha MP2 energy         = %20.15f\n", eAA);
        outfile->Printf( "\t Alpha - Beta  MP2 energy         = %20.15f\n", eAB);
        outfile->Printf( "\t Beta  - Beta  MP2 energy         = %20.15f\n", eBB);
        outfile->Printf( "\t Total MP2 correlation energy     = %20.15f\n", eAA + eAB + eBB);
        outfile->Printf( "\t*Total MP2 energy                 = %20.15f\n", new_total_energy_);

        Process::environment.globals["MP2 TOTAL ENERGY"] = new_total_energy_;
        Process::environment.globals["MP2 CORRELATION ENERGY"] = eAA + eAB + eBB;

        psio_->close(PSIF_LIBTRANS_DPD, 1);
    }
    else if(guess == "CC" || guess == "BCC"){
        outfile->Printf( "\tReading existing coupled cluster amplitudes\n\n");
        psio_->open(PSIF_CC_TAMPS, PSIO_OPEN_OLD);
        dpdbuf4 T2;
        // Copy the AA amplitudes from CCEnergy
        global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, ID("[O>O]-"), ID("[V>V]-"),
                      ID("[O>O]-"), ID("[V>V]-"), 0, "tIJAB");
        global_dpd_->buf4_copy(&T2, PSIF_DCFT_DPD, "Lambda <OO|VV>");
        global_dpd_->buf4_close(&T2);
        // Copy the AB amplitudes from CCEnergy
        global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, ID("[O,o]"), ID("[V,v]"),
                      ID("[O,o]"), ID("[V,v]"), 0, "tIjAb");
        global_dpd_->buf4_copy(&T2, PSIF_DCFT_DPD, "Lambda <Oo|Vv>");
        global_dpd_->buf4_close(&T2);
        // Copy the BB amplitudes from CCEnergy
        global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, ID("[o>o]-"), ID("[v>v]-"),
                      ID("[o>o]-"), ID("[v>v]-"), 0, "tijab");
        global_dpd_->buf4_copy(&T2, PSIF_DCFT_DPD, "Lambda <oo|vv>");
        global_dpd_->buf4_close(&T2);
        psio_->close(PSIF_CC_TAMPS, 1);
    }

    dcft_timer_off("DCFTSolver::mp2_guess()");
}


}} // Namespaces
