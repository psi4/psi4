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
    fprintf(outfile, "\n\n\tComputing MP2 amplitude guess...\n\n"); fflush(outfile);

    // If NSO basis is requested - do the OOVV integral transform only
    if (options_.get_str("DCFT_BASIS") == "NSO") {
        _ints->update_orbitals();
        if(print_ > 1){
            fprintf(outfile, "\tTransforming integrals...\n");
            fflush(outfile);
        }
        _ints->set_print(print_ - 2 >= 0 ? print_ - 2 : 0);

        // Generate the integrals in various spaces in chemists' notation for OOVV
        _ints->transform_tei(MOSpace::occ, MOSpace::vir, MOSpace::occ, MOSpace::vir);

        /*
         * Re-sort the chemists' notation integrals to physisists' notation
         * (pq|rs) = <pr|qs>
         */

        // The integral object closes this file - we need to re-open it.
        psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

        /*
         * Re-sort the chemists' notation integrals to physisicts' notation
         * (pq|rs) = <pr|qs>
         */
        dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                      ID("[O,V]"), ID("[O,V]"),0, "MO Ints (OV|OV)");
        dpd_buf4_sort(&I, PSIF_LIBTRANS_DPD, prqs, ID("[O,O]"), ID("[V,V]"), "MO Ints <OO|VV>");
        dpd_buf4_close(&I);

        dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[o,v]"),
                      ID("[O,V]"), ID("[o,v]"), 0, "MO Ints (OV|ov)");
        dpd_buf4_sort(&I, PSIF_LIBTRANS_DPD, prqs, ID("[O,o]"), ID("[V,v]"), "MO Ints <Oo|Vv>");
        dpd_buf4_close(&I);

        dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[o,v]"),
                      ID("[o,v]"), ID("[o,v]"), 0, "MO Ints (ov|ov)");
        dpd_buf4_sort(&I, PSIF_LIBTRANS_DPD, prqs, ID("[o,o]"), ID("[v,v]"), "MO Ints <oo|vv>");
        dpd_buf4_close(&I);

        build_denominators();

        psio_->close(PSIF_LIBTRANS_DPD, 1);
    }
    else {
        transform_integrals();
    }
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

    std::string guess = options_.get_str("DCFT_GUESS");
    if(guess == "CC" || guess == "BCC"){
        fprintf(outfile, "\tReading existing coupled cluster amplitudes\n");
        psio_->open(CC_TAMPS, PSIO_OPEN_OLD);
        dpdbuf4 T2;
        // Copy the AA amplitudes from CCEnergy
        dpd_buf4_init(&T2, CC_TAMPS, 0, ID("[O>O]-"), ID("[V>V]-"),
                      ID("[O>O]-"), ID("[V>V]-"), 0, "tIJAB");
        dpd_buf4_copy(&T2, PSIF_DCFT_DPD, "Lambda <OO|VV>");
        dpd_buf4_close(&T2);
        // Copy the AB amplitudes from CCEnergy
        dpd_buf4_init(&T2, CC_TAMPS, 0, ID("[O,o]"), ID("[V,v]"),
                      ID("[O,o]"), ID("[V,v]"), 0, "tIjAb");
        dpd_buf4_copy(&T2, PSIF_DCFT_DPD, "Lambda <Oo|Vv>");
        dpd_buf4_close(&T2);
        // Copy the BB amplitudes from CCEnergy
        dpd_buf4_init(&T2, CC_TAMPS, 0, ID("[o>o]-"), ID("[v>v]-"),
                      ID("[o>o]-"), ID("[v>v]-"), 0, "tijab");
        dpd_buf4_copy(&T2, PSIF_DCFT_DPD, "Lambda <oo|vv>");
        dpd_buf4_close(&T2);
        psio_->close(CC_TAMPS, 1);
    }

    // If user specified the DCFT basis to be NSO, then form tau and transform everything to NSO basis
    if (options_.get_str("DCFT_BASIS") == "NSO") {
        build_tau();
        form_nso_basis();
    }

//    exit(1);

    psio_->close(PSIF_LIBTRANS_DPD, 1);

    dcft_timer_off("DCFTSolver::mp2_guess()");
}


}} // Namespaces
