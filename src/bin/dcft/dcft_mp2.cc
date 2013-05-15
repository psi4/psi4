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
#include "defines.h"
#include <psifiles.h>
#include <vector>
#include <liboptions/liboptions.h>
#include <libpsio/psio.hpp>
#include <libtrans/integraltransform.h>
#include <libtrans/mospace.h>
#include <libdpd/dpd.h>
#include <libdiis/diismanager.h>
#include <libchkpt/chkpt.hpp>

using namespace boost;

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
    std::vector<boost::shared_ptr<MOSpace> > spaces;
    spaces.push_back(MOSpace::occ);
    spaces.push_back(MOSpace::vir);
    // This wavefunction is really the global reference wavefunction
    boost::shared_ptr<Wavefunction> wfn = Process::environment.wavefunction();
    _ints = new IntegralTransform(wfn, spaces, IntegralTransform::Unrestricted);
    _ints->set_keep_iwl_so_ints(true);
    _ints->set_keep_dpd_so_ints(true);
    dpd_set_default(_ints->get_dpd_id());

    fprintf(outfile, "\n\n\tTransforming two-electron integrals...\n");
    transform_integrals();

    std::string guess = options_.get_str("DCFT_GUESS");

    if (guess == "MP2") {
        fprintf(outfile, "\tComputing MP2 amplitude guess...\n\n"); fflush(outfile);

        psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

        /*
        * L_ijab = <ij||ab> / D_ijab
        */

        // L_IJAB = <IJ||AB> / D_IJAB
        dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O>O]-"), ID("[V>V]-"),
                      ID("[O,O]"), ID("[V,V]"), 1, "MO Ints <OO|VV>");
        dpd_buf4_copy(&I, PSIF_DCFT_DPD, "Lambda <OO|VV>");
        dpd_buf4_close(&I);
        dpd_buf4_init(&I, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                      ID("[O>O]-"), ID("[V>V]-"), 0, "Lambda <OO|VV>");
        dpd_buf4_init(&D, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                      ID("[O>=O]+"), ID("[V>=V]+"), 0, "D <OO|VV>");
        dpd_buf4_dirprd(&D, &I);
        dpd_buf4_close(&I);
        dpd_buf4_close(&D);

        // L_IjAb = <Ij|Ab> / D_IjAb
        dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                      ID("[O,o]"), ID("[V,v]"), 0, "MO Ints <Oo|Vv>");
        dpd_buf4_copy(&I, PSIF_DCFT_DPD, "Lambda <Oo|Vv>");
        dpd_buf4_close(&I);
        dpd_buf4_init(&I, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                      ID("[O,o]"), ID("[V,v]"), 0, "Lambda <Oo|Vv>");
        dpd_buf4_init(&D, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                      ID("[O,o]"), ID("[V,v]"), 0, "D <Oo|Vv>");
        dpd_buf4_dirprd(&D, &I);
        dpd_buf4_close(&I);
        dpd_buf4_close(&D);

        // L_ijab = <ij||ab> / D_ijab
        dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o>o]-"), ID("[v>v]-"),
                      ID("[o,o]"), ID("[v,v]"), 1, "MO Ints <oo|vv>");
        dpd_buf4_copy(&I, PSIF_DCFT_DPD, "Lambda <oo|vv>");
        dpd_buf4_close(&I);
        dpd_buf4_init(&I, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                      ID("[o>o]-"), ID("[v>v]-"), 0, "Lambda <oo|vv>");
        dpd_buf4_init(&D, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                      ID("[o>=o]+"), ID("[v>=v]+"), 0, "D <oo|vv>");
        dpd_buf4_dirprd(&D, &I);
        dpd_buf4_close(&I);
        dpd_buf4_close(&D);


        /*
     * E = 1/4 L_IJAB <IJ||AB>
     *        +L_IjAb <Ij|Ab>
     *    +1/4 L_ijab <ij||ab>
     */
        dpdbuf4 L;
        // Alpha - Alpha
        dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                      ID("[O,O]"), ID("[V,V]"), 1, "MO Ints <OO|VV>");
        dpd_buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                      ID("[O>O]-"), ID("[V>V]-"), 0, "Lambda <OO|VV>");
        double eAA = 0.25 * dpd_buf4_dot(&L, &I);
        dpd_buf4_close(&I);
        dpd_buf4_close(&L);

        // Alpha - Beta
        dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                      ID("[O,o]"), ID("[V,v]"), 0, "MO Ints <Oo|Vv>");
        dpd_buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                      ID("[O,o]"), ID("[V,v]"), 0, "Lambda <Oo|Vv>");
        double eAB = dpd_buf4_dot(&L, &I);
        dpd_buf4_close(&I);
        dpd_buf4_close(&L);

        // Beta - Beta
        dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                      ID("[o,o]"), ID("[v,v]"), 1, "MO Ints <oo|vv>");
        dpd_buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                      ID("[o>o]-"), ID("[v>v]-"), 0, "Lambda <oo|vv>");
        double eBB = 0.25 * dpd_buf4_dot(&L, &I);
        dpd_buf4_close(&I);
        dpd_buf4_close(&L);

        new_total_energy_ = scf_energy_ + eAA + eAB + eBB;
        fprintf(outfile, "\t*Total Hartree-Fock energy        = %20.15f\n", scf_energy_);
        fprintf(outfile, "\t Alpha - Alpha MP2 energy         = %20.15f\n", eAA);
        fprintf(outfile, "\t Alpha - Beta  MP2 energy         = %20.15f\n", eAB);
        fprintf(outfile, "\t Beta  - Beta  MP2 energy         = %20.15f\n", eBB);
        fprintf(outfile, "\t Total MP2 correlation energy     = %20.15f\n", eAA + eAB + eBB);
        fprintf(outfile, "\t*Total MP2 energy                 = %20.15f\n", new_total_energy_);

        Process::environment.globals["MP2 TOTAL ENERGY"] = new_total_energy_;
        Process::environment.globals["MP2 CORRELATION ENERGY"] = eAA + eAB + eBB;

        psio_->close(PSIF_LIBTRANS_DPD, 1);
    }
    else if(guess == "CC" || guess == "BCC"){
        fprintf(outfile, "\tReading existing coupled cluster amplitudes\n\n");
        psio_->open(PSIF_CC_TAMPS, PSIO_OPEN_OLD);
        dpdbuf4 T2;
        // Copy the AA amplitudes from CCEnergy
        dpd_buf4_init(&T2, PSIF_CC_TAMPS, 0, ID("[O>O]-"), ID("[V>V]-"),
                      ID("[O>O]-"), ID("[V>V]-"), 0, "tIJAB");
        dpd_buf4_copy(&T2, PSIF_DCFT_DPD, "Lambda <OO|VV>");
        dpd_buf4_close(&T2);
        // Copy the AB amplitudes from CCEnergy
        dpd_buf4_init(&T2, PSIF_CC_TAMPS, 0, ID("[O,o]"), ID("[V,v]"),
                      ID("[O,o]"), ID("[V,v]"), 0, "tIjAb");
        dpd_buf4_copy(&T2, PSIF_DCFT_DPD, "Lambda <Oo|Vv>");
        dpd_buf4_close(&T2);
        // Copy the BB amplitudes from CCEnergy
        dpd_buf4_init(&T2, PSIF_CC_TAMPS, 0, ID("[o>o]-"), ID("[v>v]-"),
                      ID("[o>o]-"), ID("[v>v]-"), 0, "tijab");
        dpd_buf4_copy(&T2, PSIF_DCFT_DPD, "Lambda <oo|vv>");
        dpd_buf4_close(&T2);
        psio_->close(PSIF_CC_TAMPS, 1);
    }

    dcft_timer_off("DCFTSolver::mp2_guess()");
}


}} // Namespaces
