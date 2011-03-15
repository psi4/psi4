#include "dcft.h"
#include "defines.h"
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
    // Start with good old SCF
    fprintf(outfile, "\n\tComputing Hartree-Fock orbital guess...\n\n"); fflush(outfile);
    perform_scf();
    // Now that SCF has run, we can get the resulting orbital info, and initialize the memory
    init_moinfo();

    // Print out orbital energies.
    print_orbital_energies();

    // Initialize the integral transformation object
    std::vector<shared_ptr<MOSpace> > spaces;
    spaces.push_back(MOSpace::occ);
    spaces.push_back(MOSpace::vir);
//    _ints = new IntegralTransform(spaces, IntegralTransform::Unrestricted,
//            IntegralTransform::DPDOnly, IntegralTransform::QTOrder,
//            IntegralTransform::None,false);
    _ints = new IntegralTransform(_chkpt, spaces, IntegralTransform::Unrestricted);
    _ints->set_keep_iwl_so_ints(true);
    _ints->set_keep_dpd_so_ints(true);
//    _ints->set_chkpt(_chkpt);
        // Make sure that we are using the correct DPD structure
    dpd_set_default(_ints->get_dpd_id());
    fprintf(outfile, "\n\n\tComputing MP2 amplitude guess...\n\n"); fflush(outfile);
    transform_integrals();
    _psio->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);


    dpdbuf4 I, D;
    /*
     * L_ijab = <ij||ab> / D_ijab
     */

    // L_IJAB = <IJ||AB> / D_IJAB
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 1, "MO Ints <OO|VV>");
    dpd_buf4_copy(&I, PSIF_DCFT_DPD, "Lambda <OO|VV>");
    dpd_buf4_close(&I);
    dpd_buf4_init(&I, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "Lambda <OO|VV>");
    dpd_buf4_init(&D, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "D <OO|VV>");
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
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 1, "MO Ints <oo|vv>");
    dpd_buf4_copy(&I, PSIF_DCFT_DPD, "Lambda <oo|vv>");
    dpd_buf4_close(&I);
    dpd_buf4_init(&I, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "Lambda <oo|vv>");
    dpd_buf4_init(&D, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "D <oo|vv>");
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
                  ID("[O,O]"), ID("[V,V]"), 0, "Lambda <OO|VV>");
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
                  ID("[o,o]"), ID("[v,v]"), 0, "Lambda <oo|vv>");
    double eBB = 0.25 * dpd_buf4_dot(&L, &I);
    dpd_buf4_close(&I);
    dpd_buf4_close(&L);

    _newTotalEnergy = _scfEnergy + eAA + eAB + eBB;
    fprintf(outfile, "\t*Total Hartree-Fock energy        = %20.15f\n", _scfEnergy);
    fprintf(outfile, "\t Alpha - Alpha MP2 energy         = %20.15f\n", eAA);
    fprintf(outfile, "\t Alpha - Beta  MP2 energy         = %20.15f\n", eAB);
    fprintf(outfile, "\t Beta  - Beta  MP2 energy         = %20.15f\n", eBB);
    fprintf(outfile, "\t Total MP2 correlation energy     = %20.15f\n", eAA + eAB + eBB);
    fprintf(outfile, "\t*Total MP2 energy                 = %20.15f\n", _newTotalEnergy);

    _psio->close(PSIF_LIBTRANS_DPD, 1);
}


}} // Namespaces
