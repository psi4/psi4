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

#include "psi4/libtrans/integraltransform.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libqt/qt.h"
#include "psi4/libiwl/iwl.h"
#include "dcft.h"
#include "defines.h"

namespace psi{ namespace dcft{

double
DCFTSolver::compute_three_particle_energy()
{

    outfile->Printf("\n\tEvaluating three-particle DCFT energy correction...\n\n");

    // Tansform OOOV and VVVO two-electron integrals and two-particle cumulants to semicanonical basis
    dcft_semicanonicalize();

    // Compute the three-particle energy contributions
    double triples_energy_aaa = compute_triples_aaa();
    outfile->Printf("\t*Lambda_3 Energy (AAA)                             = %20.15f\n", triples_energy_aaa);

    double triples_energy_aab = compute_triples_aab();
    outfile->Printf("\t*Lambda_3 Energy (AAB)                             = %20.15f\n", triples_energy_aab);

    double triples_energy_abb = compute_triples_abb();
    outfile->Printf("\t*Lambda_3 Energy (ABB)                             = %20.15f\n", triples_energy_abb);

    double triples_energy_bbb = compute_triples_bbb();
    outfile->Printf("\t*Lambda_3 Energy (BBB)                             = %20.15f\n\n", triples_energy_bbb);

    return triples_energy_aaa + triples_energy_aab + triples_energy_abb + triples_energy_bbb;

}

void
DCFTSolver::dcft_semicanonicalize()
{

    // If OOOV or VVVO integrals have not yet been transformed - do it
    if ((options_.get_str("ALGORITHM") != "QC" || !options_.get_bool("QC_COUPLING")
         || options_.get_str("QC_TYPE") != "SIMULTANEOUS") && !orbital_optimized_) {
        outfile->Printf("\tTransforming OVVV and OOOV integrals ... \t\t\t");
        transform_integrals_triples();
        outfile->Printf("DONE\n");
    }

    // Compute transformation to the semicanonical basis and dump it to disk
    dump_semicanonical();

    // Transform <OV||VV> integrals to the semicanonical basis
    outfile->Printf( "\tSemicanonicalizing OVVV integrals ... \t\t\t");
    semicanonicalize_gbar_ovvv();
    outfile->Printf( "DONE\n");

    outfile->Printf( "\tSemicanonicalizing OOOV integrals ... \t\t\t");
    // Transform <OO||OV> integrals to the semicanonical basis
    semicanonicalize_gbar_ooov();
    outfile->Printf( "DONE\n");

    outfile->Printf( "\tSemicanonicalizing density cumulant ...\t\t\t");
    // Transform OOVV cumulants to the semicanonical basis
    semicanonicalize_dc();
    outfile->Printf( "DONE\n\n");

}

void
DCFTSolver::transform_integrals_triples()
{
    _ints->update_orbitals();
    _ints->set_print(print_ - 2 >= 0 ? print_ - 2 : 0);

    // Compute the integrals needed for the MO Hessian
    _ints->transform_tei(MOSpace::vir, MOSpace::occ, MOSpace::occ, MOSpace::occ);
    _ints->transform_tei(MOSpace::occ, MOSpace::occ, MOSpace::vir, MOSpace::occ);
    _ints->transform_tei(MOSpace::occ, MOSpace::vir, MOSpace::vir, MOSpace::vir);
    _ints->transform_tei(MOSpace::vir, MOSpace::vir, MOSpace::occ, MOSpace::vir);

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

    // VVVO and OOOV integrals are needed for the QC algorithm
    sort_OOOV_integrals();

    sort_OVVV_integrals();

    psio_->close(PSIF_LIBTRANS_DPD, 1);

}


void
DCFTSolver::dump_semicanonical()
{

    // Zero out occupied-virtual part
    for(int h = 0; h < nirrep_; ++h){
        for(int i = 0; i < naoccpi_[h]; ++i){
            for(int a = naoccpi_[h]; a < nmopi_[h]; ++a){
                Ftilde_a_->set(h, i, a, 0.0);
                Ftilde_a_->set(h, a, i, 0.0);
            }
        }
        for(int i = 0; i < nboccpi_[h]; ++i){
            for(int a = nboccpi_[h]; a < nmopi_[h]; ++a){
                Ftilde_b_->set(h, i, a, 0.0);
                Ftilde_b_->set(h, a, i, 0.0);
            }
        }
    }

    // Diagonalize F0 to get transformation matrix to semicanonical basis
    SharedMatrix a_evecs (new Matrix ("F0 Eigenvectors (Alpha)", nirrep_, nmopi_, nmopi_));
    SharedMatrix b_evecs (new Matrix ("F0 Eigenvectors (Beta)", nirrep_, nmopi_, nmopi_));
    SharedVector a_evals (new Vector ("F0 Eigenvalues (Alpha)", nirrep_, nmopi_));
    SharedVector b_evals (new Vector ("F0 Eigenvalues (Beta)", nirrep_, nmopi_));

    Ftilde_a_->diagonalize(a_evecs, a_evals);
    Ftilde_a_->zero();
    Ftilde_a_->set_diagonal(a_evals);
    Ftilde_b_->diagonalize(b_evecs, b_evals);
    Ftilde_b_->zero();
    Ftilde_b_->set_diagonal(b_evals);

    // Write the transformation matrix to disk
    dpdfile2 U_OO, U_VV, U_oo, U_vv;
    global_dpd_->file2_init(&U_OO, PSIF_DCFT_DPD, 0, ID('O'), ID('O'), "U <O|O>");
    global_dpd_->file2_init(&U_oo, PSIF_DCFT_DPD, 0, ID('o'), ID('o'), "U <o|o>");
    global_dpd_->file2_init(&U_VV, PSIF_DCFT_DPD, 0, ID('V'), ID('V'), "U <V|V>");
    global_dpd_->file2_init(&U_vv, PSIF_DCFT_DPD, 0, ID('v'), ID('v'), "U <v|v>");

    global_dpd_->file2_mat_init(&U_OO);
    global_dpd_->file2_mat_init(&U_oo);
    global_dpd_->file2_mat_init(&U_VV);
    global_dpd_->file2_mat_init(&U_vv);

    for(int h = 0; h < nirrep_; ++h){
        if(nsopi_[h] == 0) continue;

        // Alpha occupied
        for(int ip = 0 ; ip < naoccpi_[h]; ++ip){
            for(int kp = 0 ; kp < naoccpi_[h]; ++kp){
                U_OO.matrix[h][ip][kp] = a_evecs->get(h, ip, kp);
            }
        }

        // Beta occupied
        for(int ip = 0 ; ip < nboccpi_[h]; ++ip){
            for(int kp = 0 ; kp < nboccpi_[h]; ++kp){
                U_oo.matrix[h][ip][kp] = b_evecs->get(h, ip, kp);
            }
        }

        // Alpha virtual
        for(int ap = 0 ; ap < navirpi_[h]; ++ap){
            for(int dp = 0 ; dp < navirpi_[h]; ++dp){
                U_VV.matrix[h][ap][dp] = a_evecs->get(h, ap + naoccpi_[h], dp + naoccpi_[h]);
            }
        }

        // Beta virtual
        for(int ap = 0 ; ap < nbvirpi_[h]; ++ap){
            for(int dp = 0 ; dp < nbvirpi_[h]; ++dp){
                U_vv.matrix[h][ap][dp] = b_evecs->get(h, ap + nboccpi_[h], dp + nboccpi_[h]);
            }
        }
    }

    global_dpd_->file2_mat_wrt(&U_OO);
    global_dpd_->file2_mat_wrt(&U_oo);
    global_dpd_->file2_mat_wrt(&U_VV);
    global_dpd_->file2_mat_wrt(&U_vv);

    global_dpd_->file2_close(&U_OO);
    global_dpd_->file2_close(&U_oo);
    global_dpd_->file2_close(&U_VV);
    global_dpd_->file2_close(&U_vv);
}

void
DCFTSolver::semicanonicalize_gbar_ovvv(){

    psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

    dpdbuf4 I, It;
    dpdfile2 U_VV, U_OO, U_vv, U_oo;

    // TODO: The transformed tensors should be packed on disk wherever possible

    global_dpd_->file2_init(&U_VV, PSIF_DCFT_DPD, 0, ID('V'), ID('V'), "U <V|V>");
    global_dpd_->file2_init(&U_OO, PSIF_DCFT_DPD, 0, ID('O'), ID('O'), "U <O|O>");
    global_dpd_->file2_init(&U_vv, PSIF_DCFT_DPD, 0, ID('v'), ID('v'), "U <v|v>");
    global_dpd_->file2_init(&U_oo, PSIF_DCFT_DPD, 0, ID('o'), ID('o'), "U <o|o>");

    //
    // OVVV
    //

    // <IA||BC> * U_CC' -> <IA||BC'>
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[V,V]"),
                  ID("[O,V]"), ID("[V,V]"), 1, "MO Ints <OV|VV>");
    global_dpd_->buf4_init(&It, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[V,V]"),
                  ID("[O,V]"), ID("[V,V]"), 0, "MO Ints <OV||VV'>");
    global_dpd_->contract424(&I, &U_VV, &It, 3, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&It);
    global_dpd_->buf4_close(&I);

    // U_BB' * <IA||BC'> -> <IA||B'C'>
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[V,V]"),
                  ID("[O,V]"), ID("[V,V]"), 0, "MO Ints <OV||VV'>");
    global_dpd_->buf4_init(&It, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[V,V]"),
                  ID("[O,V]"), ID("[V,V]"), 0, "MO Ints <OV||V'V'>");
    global_dpd_->contract244(&U_VV, &I, &It, 0, 2, 1, 1.0, 0.0);
    global_dpd_->buf4_close(&It);
    global_dpd_->buf4_close(&I);

    // <IA||B'C'> * U_AA' -> <IA'||B'C'>
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[V,V]"),
                  ID("[O,V]"), ID("[V,V]"), 0, "MO Ints <OV||V'V'>");
    global_dpd_->buf4_init(&It, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[V,V]"),
                  ID("[O,V]"), ID("[V,V]"), 0, "MO Ints <OV'||V'V'>");
    global_dpd_->contract424(&I, &U_VV, &It, 1, 0, 1, 1.0, 0.0);
    global_dpd_->buf4_close(&It);
    global_dpd_->buf4_close(&I);

    // U_II' * <IA'||B'C'> -> <I'A'||B'C'>
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[V,V]"),
                  ID("[O,V]"), ID("[V,V]"), 0, "MO Ints <OV'||V'V'>");
    global_dpd_->buf4_init(&It, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[V,V]"),
                  ID("[O,V]"), ID("[V,V]"), 0, "MO Ints <O'V'||V'V'>");
    global_dpd_->contract244(&U_OO, &I, &It, 0, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&It);
    global_dpd_->buf4_close(&I);

    //
    // OvVv
    //

    // <Ia|Bc> * U_cc' -> <Ia|Bc'>
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,v]"), ID("[V,v]"),
                  ID("[O,v]"), ID("[V,v]"), 0, "MO Ints <Ov|Vv>");
    global_dpd_->buf4_init(&It, PSIF_LIBTRANS_DPD, 0, ID("[O,v]"), ID("[V,v]"),
                  ID("[O,v]"), ID("[V,v]"), 0, "MO Ints <Ov|Vv'>");
    global_dpd_->contract424(&I, &U_vv, &It, 3, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&It);
    global_dpd_->buf4_close(&I);

    // U_BB' * <Ia|Bc'> -> <Ia|B'c'>
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,v]"), ID("[V,v]"),
                  ID("[O,v]"), ID("[V,v]"), 0, "MO Ints <Ov|Vv'>");
    global_dpd_->buf4_init(&It, PSIF_LIBTRANS_DPD, 0, ID("[O,v]"), ID("[V,v]"),
                  ID("[O,v]"), ID("[V,v]"), 0, "MO Ints <Ov|V'v'>");
    global_dpd_->contract244(&U_VV, &I, &It, 0, 2, 1, 1.0, 0.0);
    global_dpd_->buf4_close(&It);
    global_dpd_->buf4_close(&I);

    // <Ia|B'c'> * U_aa' -> <Ia'|B'c'>
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,v]"), ID("[V,v]"),
                  ID("[O,v]"), ID("[V,v]"), 0, "MO Ints <Ov|V'v'>");
    global_dpd_->buf4_init(&It, PSIF_LIBTRANS_DPD, 0, ID("[O,v]"), ID("[V,v]"),
                  ID("[O,v]"), ID("[V,v]"), 0, "MO Ints <Ov'|V'v'>");
    global_dpd_->contract424(&I, &U_vv, &It, 1, 0, 1, 1.0, 0.0);
    global_dpd_->buf4_close(&It);
    global_dpd_->buf4_close(&I);

    // U_II' * <Ia'|B'c'> -> <I'a'|B'c'>
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,v]"), ID("[V,v]"),
                  ID("[O,v]"), ID("[V,v]"), 0, "MO Ints <Ov'|V'v'>");
    global_dpd_->buf4_init(&It, PSIF_LIBTRANS_DPD, 0, ID("[O,v]"), ID("[V,v]"),
                  ID("[O,v]"), ID("[V,v]"), 0, "MO Ints <O'v'|V'v'>");
    global_dpd_->contract244(&U_OO, &I, &It, 0, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&It);
    global_dpd_->buf4_close(&I);

    //
    // oVvV
    //

    // <iA|bC> * U_CC' -> <iA|bC'>
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,V]"), ID("[v,V]"),
                  ID("[o,V]"), ID("[v,V]"), 0, "MO Ints <oV|vV>");
    global_dpd_->buf4_init(&It, PSIF_LIBTRANS_DPD, 0, ID("[o,V]"), ID("[v,V]"),
                  ID("[o,V]"), ID("[v,V]"), 0, "MO Ints <oV|vV'>");
    global_dpd_->contract424(&I, &U_VV, &It, 3, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&It);
    global_dpd_->buf4_close(&I);

    // U_bb' * <iA|bC'> -> <iA|b'C'>
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,V]"), ID("[v,V]"),
                  ID("[o,V]"), ID("[v,V]"), 0, "MO Ints <oV|vV'>");
    global_dpd_->buf4_init(&It, PSIF_LIBTRANS_DPD, 0, ID("[o,V]"), ID("[v,V]"),
                  ID("[o,V]"), ID("[v,V]"), 0, "MO Ints <oV|v'V'>");
    global_dpd_->contract244(&U_vv, &I, &It, 0, 2, 1, 1.0, 0.0);
    global_dpd_->buf4_close(&It);
    global_dpd_->buf4_close(&I);

    // <iA|b'C'> * U_AA' -> <iA'|b'C'>
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,V]"), ID("[v,V]"),
                  ID("[o,V]"), ID("[v,V]"), 0, "MO Ints <oV|v'V'>");
    global_dpd_->buf4_init(&It, PSIF_LIBTRANS_DPD, 0, ID("[o,V]"), ID("[v,V]"),
                  ID("[o,V]"), ID("[v,V]"), 0, "MO Ints <oV'|v'V'>");
    global_dpd_->contract424(&I, &U_VV, &It, 1, 0, 1, 1.0, 0.0);
    global_dpd_->buf4_close(&It);
    global_dpd_->buf4_close(&I);

    // U_ii' * <iA'|b'C'> -> <i'A'|b'C'>
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,V]"), ID("[v,V]"),
                  ID("[o,V]"), ID("[v,V]"), 0, "MO Ints <oV'|v'V'>");
    global_dpd_->buf4_init(&It, PSIF_LIBTRANS_DPD, 0, ID("[o,V]"), ID("[v,V]"),
                  ID("[o,V]"), ID("[v,V]"), 0, "MO Ints <o'V'|v'V'>");
    global_dpd_->contract244(&U_oo, &I, &It, 0, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&It);
    global_dpd_->buf4_close(&I);

    //
    // ovvv
    //

    // <ia||bc> * U_cc' -> <ia||bc'>
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[v,v]"),
                  ID("[o,v]"), ID("[v,v]"), 1, "MO Ints <ov|vv>");
    global_dpd_->buf4_init(&It, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[v,v]"),
                  ID("[o,v]"), ID("[v,v]"), 0, "MO Ints <ov||vv'>");
    global_dpd_->contract424(&I, &U_vv, &It, 3, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&It);
    global_dpd_->buf4_close(&I);

    // U_bb' * <ia||bc'> -> <ia||b'c'>
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[v,v]"),
                  ID("[o,v]"), ID("[v,v]"), 0, "MO Ints <ov||vv'>");
    global_dpd_->buf4_init(&It, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[v,v]"),
                  ID("[o,v]"), ID("[v,v]"), 0, "MO Ints <ov||v'v'>");
    global_dpd_->contract244(&U_vv, &I, &It, 0, 2, 1, 1.0, 0.0);
    global_dpd_->buf4_close(&It);
    global_dpd_->buf4_close(&I);

    // <ia||b'c'> * U_aa' -> <ia'||b'c'>
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[v,v]"),
                  ID("[o,v]"), ID("[v,v]"), 0, "MO Ints <ov||v'v'>");
    global_dpd_->buf4_init(&It, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[v,v]"),
                  ID("[o,v]"), ID("[v,v]"), 0, "MO Ints <ov'||v'v'>");
    global_dpd_->contract424(&I, &U_vv, &It, 1, 0, 1, 1.0, 0.0);
    global_dpd_->buf4_close(&It);
    global_dpd_->buf4_close(&I);

    // U_ii' * <ia'||b'c'> -> <i'a'||b'c'>
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[v,v]"),
                  ID("[o,v]"), ID("[v,v]"), 0, "MO Ints <ov'||v'v'>");
    global_dpd_->buf4_init(&It, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[v,v]"),
                  ID("[o,v]"), ID("[v,v]"), 0, "MO Ints <o'v'||v'v'>");
    global_dpd_->contract244(&U_oo, &I, &It, 0, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&It);
    global_dpd_->buf4_close(&I);

    global_dpd_->file2_close(&U_VV);
    global_dpd_->file2_close(&U_OO);
    global_dpd_->file2_close(&U_vv);
    global_dpd_->file2_close(&U_oo);

    psio_->close(PSIF_LIBTRANS_DPD, 1);

}

void
DCFTSolver::semicanonicalize_gbar_ooov(){

    psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

    dpdbuf4 I, It;
    dpdfile2 U_VV, U_OO, U_vv, U_oo;

    // TODO: The transformed tensors should be packed on disk wherever possible

    global_dpd_->file2_init(&U_VV, PSIF_DCFT_DPD, 0, ID('V'), ID('V'), "U <V|V>");
    global_dpd_->file2_init(&U_OO, PSIF_DCFT_DPD, 0, ID('O'), ID('O'), "U <O|O>");
    global_dpd_->file2_init(&U_vv, PSIF_DCFT_DPD, 0, ID('v'), ID('v'), "U <v|v>");
    global_dpd_->file2_init(&U_oo, PSIF_DCFT_DPD, 0, ID('o'), ID('o'), "U <o|o>");

    //
    // OOOV
    //

    // <AI||JK> * U_KK' -> <AI||JK'>
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,O]"), ID("[O,O]"),
                  ID("[V,O]"), ID("[O,O]"), 1, "MO Ints <VO|OO>");
    global_dpd_->buf4_init(&It, PSIF_LIBTRANS_DPD, 0, ID("[V,O]"), ID("[O,O]"),
                  ID("[V,O]"), ID("[O,O]"), 0, "MO Ints <VO||OO'>");
    global_dpd_->contract424(&I, &U_OO, &It, 3, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&It);
    global_dpd_->buf4_close(&I);

    // U_JJ' * <AI||JK'> -> <AI||J'K'>
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,O]"), ID("[O,O]"),
                  ID("[V,O]"), ID("[O,O]"), 0, "MO Ints <VO||OO'>");
    global_dpd_->buf4_init(&It, PSIF_LIBTRANS_DPD, 0, ID("[V,O]"), ID("[O,O]"),
                  ID("[V,O]"), ID("[O,O]"), 0, "MO Ints <VO||O'O'>");
    global_dpd_->contract244(&U_OO, &I, &It, 0, 2, 1, 1.0, 0.0);
    global_dpd_->buf4_close(&It);
    global_dpd_->buf4_close(&I);

    // <AI||J'K'> * U_II' -> <AI'||J'K'>
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,O]"), ID("[O,O]"),
                  ID("[V,O]"), ID("[O,O]"), 0, "MO Ints <VO||O'O'>");
    global_dpd_->buf4_init(&It, PSIF_LIBTRANS_DPD, 0, ID("[V,O]"), ID("[O,O]"),
                  ID("[V,O]"), ID("[O,O]"), 0, "MO Ints <VO'||O'O'>");
    global_dpd_->contract424(&I, &U_OO, &It, 1, 0, 1, 1.0, 0.0);
    global_dpd_->buf4_close(&It);
    global_dpd_->buf4_close(&I);

    // U_AA' * <AI'||J'K'> -> <A'I'||J'K'>
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,O]"), ID("[O,O]"),
                  ID("[V,O]"), ID("[O,O]"), 0, "MO Ints <VO'||O'O'>");
    global_dpd_->buf4_init(&It, PSIF_LIBTRANS_DPD, 0, ID("[V,O]"), ID("[O,O]"),
                  ID("[V,O]"), ID("[O,O]"), 0, "MO Ints <V'O'||O'O'>");
    global_dpd_->contract244(&U_VV, &I, &It, 0, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&It);
    global_dpd_->buf4_close(&I);

    // <A'I'||J'K'> -> <K'J'||I'A'>
    global_dpd_->buf4_init(&It, PSIF_LIBTRANS_DPD, 0, ID("[V,O]"), ID("[O,O]"),
                  ID("[V,O]"), ID("[O,O]"), 0, "MO Ints <V'O'||O'O'>");
    global_dpd_->buf4_sort(&It, PSIF_LIBTRANS_DPD, srqp, ID("[O,O]"), ID("[O,V]"), "MO Ints <O'O'||O'V'>");
    global_dpd_->buf4_close(&It);


    //
    // OoOv
    //

    // <Ia|Jk> * U_kk' -> <Ia|Jk'>
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,v]"), ID("[O,o]"),
                  ID("[O,v]"), ID("[O,o]"), 0, "MO Ints <Ov|Oo>");
    global_dpd_->buf4_init(&It, PSIF_LIBTRANS_DPD, 0, ID("[O,v]"), ID("[O,o]"),
                  ID("[O,v]"), ID("[O,o]"), 0, "MO Ints <Ov|Oo'>");
    global_dpd_->contract424(&I, &U_oo, &It, 3, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&It);
    global_dpd_->buf4_close(&I);

    // U_JJ' * <Ia|Jk'> -> <Ia|J'k'>
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,v]"), ID("[O,o]"),
                  ID("[O,v]"), ID("[O,o]"), 0, "MO Ints <Ov|Oo'>");
    global_dpd_->buf4_init(&It, PSIF_LIBTRANS_DPD, 0, ID("[O,v]"), ID("[O,o]"),
                  ID("[O,v]"), ID("[O,o]"), 0, "MO Ints <Ov|O'o'>");
    global_dpd_->contract244(&U_OO, &I, &It, 0, 2, 1, 1.0, 0.0);
    global_dpd_->buf4_close(&It);
    global_dpd_->buf4_close(&I);

    // <Ia|J'k'> * U_aa' -> <Ia'|J'k'>
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,v]"), ID("[O,o]"),
                  ID("[O,v]"), ID("[O,o]"), 0, "MO Ints <Ov|O'o'>");
    global_dpd_->buf4_init(&It, PSIF_LIBTRANS_DPD, 0, ID("[O,v]"), ID("[O,o]"),
                  ID("[O,v]"), ID("[O,o]"), 0, "MO Ints <Ov'|O'o'>");
    global_dpd_->contract424(&I, &U_vv, &It, 1, 0, 1, 1.0, 0.0);
    global_dpd_->buf4_close(&It);
    global_dpd_->buf4_close(&I);

    // U_II' * <Ia'|J'k'> -> <I'a'|J'k'>
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,v]"), ID("[O,o]"),
                  ID("[O,v]"), ID("[O,o]"), 0, "MO Ints <Ov'|O'o'>");
    global_dpd_->buf4_init(&It, PSIF_LIBTRANS_DPD, 0, ID("[O,v]"), ID("[O,o]"),
                  ID("[O,v]"), ID("[O,o]"), 0, "MO Ints <O'v'|O'o'>");
    global_dpd_->contract244(&U_OO, &I, &It, 0, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&It);
    global_dpd_->buf4_close(&I);

    // <I'a'|J'k'> -> <J'k'|I'a'>
    global_dpd_->buf4_init(&It, PSIF_LIBTRANS_DPD, 0, ID("[O,v]"), ID("[O,o]"),
                  ID("[O,v]"), ID("[O,o]"), 0, "MO Ints <O'v'|O'o'>");
    global_dpd_->buf4_sort(&It, PSIF_LIBTRANS_DPD, rspq, ID("[O,o]"), ID("[O,v]"), "MO Ints <O'o'|O'v'>");
    global_dpd_->buf4_close(&It);


    //
    // oOoV
    //

    // <iA|jK> * U_KK' -> <iA|jK'>
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,V]"), ID("[o,O]"),
                  ID("[o,V]"), ID("[o,O]"), 0, "MO Ints <oV|oO>");
    global_dpd_->buf4_init(&It, PSIF_LIBTRANS_DPD, 0, ID("[o,V]"), ID("[o,O]"),
                  ID("[o,V]"), ID("[o,O]"), 0, "MO Ints <oV|oO'>");
    global_dpd_->contract424(&I, &U_OO, &It, 3, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&It);
    global_dpd_->buf4_close(&I);

    // U_jj' * <iA|jK'> -> <iA|j'K'>
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,V]"), ID("[o,O]"),
                  ID("[o,V]"), ID("[o,O]"), 0, "MO Ints <oV|oO'>");
    global_dpd_->buf4_init(&It, PSIF_LIBTRANS_DPD, 0, ID("[o,V]"), ID("[o,O]"),
                  ID("[o,V]"), ID("[o,O]"), 0, "MO Ints <oV|o'O'>");
    global_dpd_->contract244(&U_oo, &I, &It, 0, 2, 1, 1.0, 0.0);
    global_dpd_->buf4_close(&It);
    global_dpd_->buf4_close(&I);

    // <iA|j'K'> * U_AA' -> <iA'|j'K'>
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,V]"), ID("[o,O]"),
                  ID("[o,V]"), ID("[o,O]"), 0, "MO Ints <oV|o'O'>");
    global_dpd_->buf4_init(&It, PSIF_LIBTRANS_DPD, 0, ID("[o,V]"), ID("[o,O]"),
                  ID("[o,V]"), ID("[o,O]"), 0, "MO Ints <oV'|o'O'>");
    global_dpd_->contract424(&I, &U_VV, &It, 1, 0, 1, 1.0, 0.0);
    global_dpd_->buf4_close(&It);
    global_dpd_->buf4_close(&I);

    // U_ii' * <iA'|j'K'> -> <i'A'|j'K'>
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,V]"), ID("[o,O]"),
                  ID("[o,V]"), ID("[o,O]"), 0, "MO Ints <oV'|o'O'>");
    global_dpd_->buf4_init(&It, PSIF_LIBTRANS_DPD, 0, ID("[o,V]"), ID("[o,O]"),
                  ID("[o,V]"), ID("[o,O]"), 0, "MO Ints <o'V'|o'O'>");
    global_dpd_->contract244(&U_oo, &I, &It, 0, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&It);
    global_dpd_->buf4_close(&I);

    // <i'A'|j'K'> -> <j'K'|i'A'>
    global_dpd_->buf4_init(&It, PSIF_LIBTRANS_DPD, 0, ID("[o,V]"), ID("[o,O]"),
                  ID("[o,V]"), ID("[o,O]"), 0, "MO Ints <o'V'|o'O'>");
    global_dpd_->buf4_sort(&It, PSIF_LIBTRANS_DPD, rspq, ID("[o,O]"), ID("[o,V]"), "MO Ints <o'O'|o'V'>");
    global_dpd_->buf4_close(&It);


    //
    // ooov
    //

    // <ai||jk> * U_kk' -> <ai||jk'>
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[v,o]"), ID("[o,o]"),
                  ID("[v,o]"), ID("[o,o]"), 1, "MO Ints <vo|oo>");
    global_dpd_->buf4_init(&It, PSIF_LIBTRANS_DPD, 0, ID("[v,o]"), ID("[o,o]"),
                  ID("[v,o]"), ID("[o,o]"), 0, "MO Ints <vo||oo'>");
    global_dpd_->contract424(&I, &U_oo, &It, 3, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&It);
    global_dpd_->buf4_close(&I);

    // U_jj' * <ai||jk> -> <ai||j'k'>
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[v,o]"), ID("[o,o]"),
                  ID("[v,o]"), ID("[o,o]"), 0, "MO Ints <vo||oo'>");
    global_dpd_->buf4_init(&It, PSIF_LIBTRANS_DPD, 0, ID("[v,o]"), ID("[o,o]"),
                  ID("[v,o]"), ID("[o,o]"), 0, "MO Ints <vo||o'o'>");
    global_dpd_->contract244(&U_oo, &I, &It, 0, 2, 1, 1.0, 0.0);
    global_dpd_->buf4_close(&It);
    global_dpd_->buf4_close(&I);

    // <ai||j'k'> * U_ii' -> <ai'||j'k'>
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[v,o]"), ID("[o,o]"),
                  ID("[v,o]"), ID("[o,o]"), 0, "MO Ints <vo||o'o'>");
    global_dpd_->buf4_init(&It, PSIF_LIBTRANS_DPD, 0, ID("[v,o]"), ID("[o,o]"),
                  ID("[v,o]"), ID("[o,o]"), 0, "MO Ints <vo'||o'o'>");
    global_dpd_->contract424(&I, &U_oo, &It, 1, 0, 1, 1.0, 0.0);
    global_dpd_->buf4_close(&It);
    global_dpd_->buf4_close(&I);

    // U_aa' * <ai'||j'k'> -> <a'i'||j'k'>
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[v,o]"), ID("[o,o]"),
                  ID("[v,o]"), ID("[o,o]"), 0, "MO Ints <vo'||o'o'>");
    global_dpd_->buf4_init(&It, PSIF_LIBTRANS_DPD, 0, ID("[v,o]"), ID("[o,o]"),
                  ID("[v,o]"), ID("[o,o]"), 0, "MO Ints <v'o'||o'o'>");
    global_dpd_->contract244(&U_vv, &I, &It, 0, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&It);
    global_dpd_->buf4_close(&I);

    // <a'i'||j'k'> -> <k'j'||i'a'>
    global_dpd_->buf4_init(&It, PSIF_LIBTRANS_DPD, 0, ID("[v,o]"), ID("[o,o]"),
                  ID("[v,o]"), ID("[o,o]"), 0, "MO Ints <v'o'||o'o'>");
    global_dpd_->buf4_sort(&It, PSIF_LIBTRANS_DPD, srqp, ID("[o,o]"), ID("[o,v]"), "MO Ints <o'o'||o'v'>");
    global_dpd_->buf4_close(&It);

    global_dpd_->file2_close(&U_VV);
    global_dpd_->file2_close(&U_OO);
    global_dpd_->file2_close(&U_vv);
    global_dpd_->file2_close(&U_oo);

    psio_->close(PSIF_LIBTRANS_DPD, 1);

}

void
DCFTSolver::semicanonicalize_dc(){

    dpdbuf4 L, Lt;
    dpdfile2 U_VV, U_OO, U_vv, U_oo;

    // TODO: The transformed tensors should be packed on disk wherever possible

    global_dpd_->file2_init(&U_VV, PSIF_DCFT_DPD, 0, ID('V'), ID('V'), "U <V|V>");
    global_dpd_->file2_init(&U_OO, PSIF_DCFT_DPD, 0, ID('O'), ID('O'), "U <O|O>");
    global_dpd_->file2_init(&U_vv, PSIF_DCFT_DPD, 0, ID('v'), ID('v'), "U <v|v>");
    global_dpd_->file2_init(&U_oo, PSIF_DCFT_DPD, 0, ID('o'), ID('o'), "U <o|o>");

    //
    // Lambda_OOVV
    //

    // Lambda <IJ|AB> * U_BB' -> Lambda <IJ|AB'>
    global_dpd_->buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O>O]-"), ID("[V>V]-"), 0, "Lambda <OO|VV>");
    global_dpd_->buf4_init(&Lt, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "Lambda <OO|VV'>");
    global_dpd_->contract424(&L, &U_VV, &Lt, 3, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&Lt);
    global_dpd_->buf4_close(&L);

    // Lambda <IJ|AB'> * U_AA' -> Lambda <IJ|A'B'>
    global_dpd_->buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "Lambda <OO|VV'>");
    global_dpd_->buf4_init(&Lt, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "Lambda <OO|V'V'>");
    global_dpd_->contract244(&U_VV, &L, &Lt, 0, 2, 1, 1.0, 0.0);
    global_dpd_->buf4_close(&Lt);
    global_dpd_->buf4_close(&L);

    // Lambda <IJ|A'B'> * U_JJ' -> Lambda <IJ'|A'B'>
    global_dpd_->buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "Lambda <OO|V'V'>");
    global_dpd_->buf4_init(&Lt, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "Lambda <OO'|V'V'>");
    global_dpd_->contract424(&L, &U_OO, &Lt, 1, 0, 1, 1.0, 0.0);
    global_dpd_->buf4_close(&Lt);
    global_dpd_->buf4_close(&L);

    // Lambda <IJ'|A'B'> * U_II' -> Lambda <I'J'|A'B'>
    global_dpd_->buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "Lambda <OO'|V'V'>");
    global_dpd_->buf4_init(&Lt, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "Lambda <O'O'|V'V'>");
    global_dpd_->contract244(&U_OO, &L, &Lt, 0, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&Lt);
    global_dpd_->buf4_close(&L);

    //
    // Lambda_OoVv & Lambda_oOvV
    //

    // Lambda <Ij|Ab> * U_bb' -> Lambda <Ij|Ab'>
    global_dpd_->buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "Lambda <Oo|Vv>");
    global_dpd_->buf4_init(&Lt, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "Lambda <Oo|Vv'>");
    global_dpd_->contract424(&L, &U_vv, &Lt, 3, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&Lt);
    global_dpd_->buf4_close(&L);

    // Lambda <Ij|Ab'> * U_AA' -> Lambda <Ij|A'b'>
    global_dpd_->buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "Lambda <Oo|Vv'>");
    global_dpd_->buf4_init(&Lt, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "Lambda <Oo|V'v'>");
    global_dpd_->contract244(&U_VV, &L, &Lt, 0, 2, 1, 1.0, 0.0);
    global_dpd_->buf4_close(&Lt);
    global_dpd_->buf4_close(&L);

    // Lambda <Ij|A'b'> * U_jj' -> Lambda <Ij'|A'b'>
    global_dpd_->buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "Lambda <Oo|V'v'>");
    global_dpd_->buf4_init(&Lt, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "Lambda <Oo'|V'v'>");
    global_dpd_->contract424(&L, &U_oo, &Lt, 1, 0, 1, 1.0, 0.0);
    global_dpd_->buf4_close(&Lt);
    global_dpd_->buf4_close(&L);

    // Lambda <Ij'|A'b'> * U_II' -> Lambda <I'j'|A'b'>
    global_dpd_->buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "Lambda <Oo'|V'v'>");
    global_dpd_->buf4_init(&Lt, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "Lambda <O'o'|V'v'>");
    global_dpd_->contract244(&U_OO, &L, &Lt, 0, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&Lt);
    global_dpd_->buf4_close(&L);

    // Lambda <I'j'|A'b'> -> Lambda <j'I'|b'A'>
    global_dpd_->buf4_init(&Lt, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "Lambda <O'o'|V'v'>");
    global_dpd_->buf4_sort(&Lt, PSIF_DCFT_DPD, qpsr, ID("[o,O]"), ID("[v,V]"), "Lambda <o'O'|v'V'>");
    global_dpd_->buf4_close(&Lt);

    global_dpd_->buf4_init(&Lt, PSIF_DCFT_DPD, 0, ID("[o,O]"), ID("[v,V]"),
                  ID("[o,O]"), ID("[v,V]"), 0, "Lambda <o'O'|v'V'>");
    global_dpd_->buf4_close(&Lt);

    //
    // Lambda_oovv
    //

    // Lambda <ij|ab> * U_bb' -> Lambda <ij|ab'>
    global_dpd_->buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o>o]-"), ID("[v>v]-"), 0, "Lambda <oo|vv>");
    global_dpd_->buf4_init(&Lt, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "Lambda <oo|vv'>");
    global_dpd_->contract424(&L, &U_vv, &Lt, 3, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&Lt);
    global_dpd_->buf4_close(&L);

    // Lambda <ij|ab'> * U_aa' -> Lambda <ij|a'b'>
    global_dpd_->buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "Lambda <oo|vv'>");
    global_dpd_->buf4_init(&Lt, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "Lambda <oo|v'v'>");
    global_dpd_->contract244(&U_vv, &L, &Lt, 0, 2, 1, 1.0, 0.0);
    global_dpd_->buf4_close(&Lt);
    global_dpd_->buf4_close(&L);

    // Lambda <ij|a'b'> * U_jj' -> Lambda <ij'|a'b'>
    global_dpd_->buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "Lambda <oo|v'v'>");
    global_dpd_->buf4_init(&Lt, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "Lambda <oo'|v'v'>");
    global_dpd_->contract424(&L, &U_oo, &Lt, 1, 0, 1, 1.0, 0.0);
    global_dpd_->buf4_close(&Lt);
    global_dpd_->buf4_close(&L);

    // Lambda <ij'|a'b'> * U_ii' -> Lambda <i'j'|a'b'>
    global_dpd_->buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "Lambda <oo'|v'v'>");
    global_dpd_->buf4_init(&Lt, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "Lambda <o'o'|v'v'>");
    global_dpd_->contract244(&U_oo, &L, &Lt, 0, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&Lt);
    global_dpd_->buf4_close(&L);

    ////////////////////////////////

    global_dpd_->file2_close(&U_VV);
    global_dpd_->file2_close(&U_OO);
    global_dpd_->file2_close(&U_vv);
    global_dpd_->file2_close(&U_oo);

}

double
DCFTSolver::compute_triples_aaa()
{

    int h;
    int Gi, Gj, Gk, Ga, Gb, Gc, Gd, Gl;
    int Gji, Gij, Gjk, Gkj, Gik, Gki, Gijk;
    int Gab, Gbc, Gac;
    int Gid, Gjd, Gkd;
    int Gil, Gjl, Gkl;
    int I, J, K, A, B, C;
    int i, j, k, a, b, c;
    int ij, ji, ik, ki, jk, kj;
    int ab;
    int cd, ad, bd;
    int id, jd, kd;
    int il, jl, kl;
    int lc, la, lb;
    double dijk, denom, ET_aaa;
    int nrows, ncols, nlinks;
    dpdbuf4 L, I_OVVV, I_OOOV;
    double ***WABC, ***WBCA, ***WACB, ***LABC;

    psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

    global_dpd_->buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                           ID("[O,O]"), ID("[V,V]"), 0, "Lambda <O'O'|V'V'>");
    global_dpd_->buf4_init(&I_OVVV, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[V,V]"),
                           ID("[O,V]"), ID("[V,V]"), 0, "MO Ints <O'V'||V'V'>");
    global_dpd_->buf4_init(&I_OOOV, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[O,V]"),
                           ID("[O,O]"), ID("[O,V]"), 0, "MO Ints <O'O'||O'V'>");

    for(h=0; h < nirrep_; h++) {
      global_dpd_->buf4_mat_irrep_init(&L, h);
      global_dpd_->buf4_mat_irrep_rd(&L, h);

      global_dpd_->buf4_mat_irrep_init(&I_OOOV, h);
      global_dpd_->buf4_mat_irrep_rd(&I_OOOV, h);
    }

    WABC = (double ***) malloc(nirrep_ * sizeof(double **));
    LABC = (double ***) malloc(nirrep_ * sizeof(double **));
    WBCA = (double ***) malloc(nirrep_ * sizeof(double **));
    WACB = (double ***) malloc(nirrep_ * sizeof(double **));

    ET_aaa = 0.0;

    for(Gi=0; Gi < nirrep_; Gi++) {
      for(Gj=0; Gj < nirrep_; Gj++) {
        for(Gk=0; Gk < nirrep_; Gk++) {

      Gij = Gji = Gi ^ Gj;
      Gjk = Gkj = Gj ^ Gk;
      Gik = Gki = Gi ^ Gk;

      Gijk = Gi ^ Gj ^ Gk;

      for(i=0; i < naoccpi_[Gi]; i++) {
        I = aocc_off_[Gi] + i;
        for(j=0; j < naoccpi_[Gj]; j++) {
          J = aocc_off_[Gj] + j;
          for(k=0; k < naoccpi_[Gk]; k++) {
            K = aocc_off_[Gk] + k;

            if(I > J && J > K) {

          ij = I_OOOV.params->rowidx[I][J];
          ji = I_OOOV.params->rowidx[J][I];
          jk = I_OOOV.params->rowidx[J][K];
          kj = I_OOOV.params->rowidx[K][J];
          ik = I_OOOV.params->rowidx[I][K];
          ki = I_OOOV.params->rowidx[K][I];

          dijk = 0.0;
          if(Ftilde_a_->rowspi()[Gi])
            dijk += Ftilde_a_->get(Gi, i, i);
          if(Ftilde_a_->rowspi()[Gj])
            dijk += Ftilde_a_->get(Gj, j, j);
          if(Ftilde_a_->rowspi()[Gk])
            dijk += Ftilde_a_->get(Gk, k, k);

          for(Gab=0; Gab < nirrep_; Gab++) {
            Gc = Gab ^ Gijk;

            WABC[Gab] = global_dpd_->dpd_block_matrix(I_OVVV.params->coltot[Gab], navirpi_[Gc]);
          }

          for(Gd=0; Gd < nirrep_; Gd++) {
            /* -lambda_jkcd * <id||ab> */
            Gab = Gid = Gi ^ Gd;
            Gc = Gjk ^ Gd;

            cd = L.col_offset[Gjk][Gc];
            id = I_OVVV.row_offset[Gid][I];

            I_OVVV.matrix[Gid] = global_dpd_->dpd_block_matrix(navirpi_[Gd], I_OVVV.params->coltot[Gid]);
            global_dpd_->buf4_mat_irrep_rd_block(&I_OVVV, Gid, id, navirpi_[Gd]);

            nrows = I_OVVV.params->coltot[Gid];
            ncols = navirpi_[Gc];
            nlinks = navirpi_[Gd];

            if(nrows && ncols && nlinks)
              C_DGEMM('t', 't', nrows, ncols, nlinks, -1.0,
                  &(I_OVVV.matrix[Gid][0][0]), nrows,
                  &(L.matrix[Gjk][jk][cd]), nlinks, 1.0,
                  &(WABC[Gab][0][0]), ncols);

            global_dpd_->free_dpd_block(I_OVVV.matrix[Gid], navirpi_[Gd], I_OVVV.params->coltot[Gid]);

            /* +lambda_ikcd * <jd||ab> */
            Gab = Gjd = Gj ^ Gd;
            Gc = Gik ^ Gd;

            cd = L.col_offset[Gik][Gc];
            jd = I_OVVV.row_offset[Gjd][J];

            I_OVVV.matrix[Gjd] = global_dpd_->dpd_block_matrix(navirpi_[Gd], I_OVVV.params->coltot[Gjd]);
            global_dpd_->buf4_mat_irrep_rd_block(&I_OVVV, Gjd, jd, navirpi_[Gd]);

            nrows = I_OVVV.params->coltot[Gjd];
            ncols = navirpi_[Gc];
            nlinks = navirpi_[Gd];

            if(nrows && ncols && nlinks)
              C_DGEMM('t', 't', nrows, ncols, nlinks, 1.0,
                  &(I_OVVV.matrix[Gjd][0][0]), nrows,
                  &(L.matrix[Gik][ik][cd]), nlinks, 1.0,
                  &(WABC[Gab][0][0]), ncols);

            global_dpd_->free_dpd_block(I_OVVV.matrix[Gjd], navirpi_[Gd], I_OVVV.params->coltot[Gjd]);

            /* +lambda_jicd * <kd||ab> */
            Gab = Gkd = Gk ^ Gd;
            Gc = Gji ^ Gd;

            cd = L.col_offset[Gji][Gc];
            kd = I_OVVV.row_offset[Gkd][K];

            I_OVVV.matrix[Gkd] = global_dpd_->dpd_block_matrix(navirpi_[Gd], I_OVVV.params->coltot[Gkd]);
            global_dpd_->buf4_mat_irrep_rd_block(&I_OVVV, Gkd, kd, navirpi_[Gd]);

            nrows = I_OVVV.params->coltot[Gkd];
            ncols = navirpi_[Gc];
            nlinks = navirpi_[Gd];

            if(nrows && ncols && nlinks)
              C_DGEMM('t', 't', nrows, ncols, nlinks, 1.0,
                  &(I_OVVV.matrix[Gkd][0][0]), nrows,
                  &(L.matrix[Gji][ji][cd]), nlinks, 1.0,
                  &(WABC[Gab][0][0]), ncols);

            global_dpd_->free_dpd_block(I_OVVV.matrix[Gkd], navirpi_[Gd], I_OVVV.params->coltot[Gkd]);

          }

          for(Gl=0; Gl < nirrep_; Gl++) {
            /* -lambda_ilab <jk||lc> */
            Gab = Gil = Gi ^ Gl;
            Gc = Gjk ^ Gl;

            lc = I_OOOV.col_offset[Gjk][Gl];
            il = L.row_offset[Gil][I];

            nrows = L.params->coltot[Gil];
            ncols = navirpi_[Gc];
            nlinks = naoccpi_[Gl];

            if(nrows && ncols && nlinks)
              C_DGEMM('t', 'n', nrows, ncols, nlinks, -1.0,
                  &(L.matrix[Gil][il][0]), nrows,
                  &(I_OOOV.matrix[Gjk][jk][lc]), ncols, 1.0,
                  &(WABC[Gab][0][0]), ncols);

            /* +lambda_jlab <ik||lc> */
            Gab = Gjl = Gj ^ Gl;
            Gc = Gik ^ Gl;

            lc = I_OOOV.col_offset[Gik][Gl];
            jl = L.row_offset[Gjl][J];

            nrows = L.params->coltot[Gjl];
            ncols = navirpi_[Gc];
            nlinks = naoccpi_[Gl];

            if(nrows && ncols && nlinks)
              C_DGEMM('t', 'n', nrows, ncols, nlinks, 1.0,
                  &(L.matrix[Gjl][jl][0]), nrows,
                  &(I_OOOV.matrix[Gik][ik][lc]), ncols, 1.0,
                  &(WABC[Gab][0][0]), ncols);

            /* +lambda_klab <ji||lc> */
            Gab = Gkl = Gk ^ Gl;
            Gc = Gji ^ Gl;

            lc = I_OOOV.col_offset[Gji][Gl];
            kl = L.row_offset[Gkl][K];

            nrows = L.params->coltot[Gkl];
            ncols = navirpi_[Gc];
            nlinks = naoccpi_[Gl];

            if(nrows && ncols && nlinks)
              C_DGEMM('t', 'n', nrows, ncols, nlinks, 1.0,
                  &(L.matrix[Gkl][kl][0]), nrows,
                  &(I_OOOV.matrix[Gji][ji][lc]), ncols, 1.0,
                  &(WABC[Gab][0][0]), ncols);
          }

          for(Gab=0; Gab < nirrep_; Gab++) {
            Gc = Gab ^ Gijk;

            WBCA[Gab] = global_dpd_->dpd_block_matrix(I_OVVV.params->coltot[Gab], navirpi_[Gc]);
          }

          for(Gd=0; Gd < nirrep_; Gd++) {
            /* -lambda_jkad * <id||bc> */
            Gbc = Gid = Gi ^ Gd;
            Ga = Gjk ^ Gd;

            ad = L.col_offset[Gjk][Ga];
            id = I_OVVV.row_offset[Gid][I];

            I_OVVV.matrix[Gid] = global_dpd_->dpd_block_matrix(navirpi_[Gd], I_OVVV.params->coltot[Gid]);
            global_dpd_->buf4_mat_irrep_rd_block(&I_OVVV, Gid, id, navirpi_[Gd]);

            nrows = I_OVVV.params->coltot[Gid];
            ncols = navirpi_[Ga];
            nlinks = navirpi_[Gd];

            if(nrows && ncols && nlinks)
              C_DGEMM('t', 't', nrows, ncols, nlinks, -1.0,
                  &(I_OVVV.matrix[Gid][0][0]), nrows,
                  &(L.matrix[Gjk][jk][ad]), nlinks, 1.0,
                  &(WBCA[Gbc][0][0]), ncols);

            global_dpd_->free_dpd_block(I_OVVV.matrix[Gid], navirpi_[Gd], I_OVVV.params->coltot[Gid]);

            /* +lambda_ikad * <jd||bc> */
            Gbc = Gjd = Gj ^ Gd;
            Ga = Gik ^ Gd;

            ad = L.col_offset[Gik][Ga];
            jd = I_OVVV.row_offset[Gjd][J];

            I_OVVV.matrix[Gjd] = global_dpd_->dpd_block_matrix(navirpi_[Gd], I_OVVV.params->coltot[Gjd]);
            global_dpd_->buf4_mat_irrep_rd_block(&I_OVVV, Gjd, jd, navirpi_[Gd]);

            nrows = I_OVVV.params->coltot[Gjd];
            ncols = navirpi_[Ga];
            nlinks = navirpi_[Gd];

            if(nrows && ncols && nlinks)
              C_DGEMM('t', 't', nrows, ncols, nlinks, 1.0,
                  &(I_OVVV.matrix[Gjd][0][0]), nrows,
                  &(L.matrix[Gik][ik][ad]), nlinks, 1.0,
                  &(WBCA[Gbc][0][0]), ncols);

            global_dpd_->free_dpd_block(I_OVVV.matrix[Gjd], navirpi_[Gd], I_OVVV.params->coltot[Gjd]);

            /* +lambda_jiad * <kd||bc> */
            Gbc = Gkd = Gk ^ Gd;
            Ga = Gji ^ Gd;

            ad = L.col_offset[Gji][Ga];
            kd = I_OVVV.row_offset[Gkd][K];

            I_OVVV.matrix[Gkd] = global_dpd_->dpd_block_matrix(navirpi_[Gd], I_OVVV.params->coltot[Gkd]);
            global_dpd_->buf4_mat_irrep_rd_block(&I_OVVV, Gkd, kd, navirpi_[Gd]);

            nrows = I_OVVV.params->coltot[Gkd];
            ncols = navirpi_[Ga];
            nlinks = navirpi_[Gd];

            if(nrows && ncols && nlinks)
              C_DGEMM('t', 't', nrows, ncols, nlinks, 1.0,
                  &(I_OVVV.matrix[Gkd][0][0]), nrows,
                  &(L.matrix[Gji][ji][ad]), nlinks, 1.0,
                  &(WBCA[Gbc][0][0]), ncols);

            global_dpd_->free_dpd_block(I_OVVV.matrix[Gkd], navirpi_[Gd], I_OVVV.params->coltot[Gkd]);

          }

          for(Gl=0; Gl < nirrep_; Gl++) {
            /* -lambda_ilbc * <jk||la> */
            Gbc = Gil = Gi ^ Gl;
            Ga = Gjk ^ Gl;

            la = I_OOOV.col_offset[Gjk][Gl];
            il = L.row_offset[Gil][I];

            nrows = L.params->coltot[Gil];
            ncols = navirpi_[Ga];
            nlinks = naoccpi_[Gl];

            if(nrows && ncols && nlinks)
              C_DGEMM('t', 'n', nrows, ncols, nlinks, -1.0,
                  &(L.matrix[Gil][il][0]), nrows,
                  &(I_OOOV.matrix[Gjk][jk][la]), ncols, 1.0,
                  &(WBCA[Gbc][0][0]), ncols);

            /* +lambda_jlbc <ik||la> */
            Gbc = Gjl = Gj ^ Gl;
            Ga = Gik ^ Gl;

            la = I_OOOV.col_offset[Gik][Gl];
            jl = L.row_offset[Gjl][J];

            nrows = L.params->coltot[Gjl];
            ncols = navirpi_[Ga];
            nlinks = naoccpi_[Gl];

            if(nrows && ncols && nlinks)
              C_DGEMM('t', 'n', nrows, ncols, nlinks, 1.0,
                  &(L.matrix[Gjl][jl][0]), nrows,
                  &(I_OOOV.matrix[Gik][ik][la]), ncols, 1.0,
                  &(WBCA[Gbc][0][0]), ncols);

            /* +lambda_klbc <ji||la> */
            Gbc = Gkl = Gk ^ Gl;
            Ga = Gji ^ Gl;

            la = I_OOOV.col_offset[Gji][Gl];
            kl = L.row_offset[Gkl][K];

            nrows = L.params->coltot[Gkl];
            ncols = navirpi_[Ga];
            nlinks = naoccpi_[Gl];

            if(nrows && ncols && nlinks)
              C_DGEMM('t', 'n', nrows, ncols, nlinks, 1.0,
                  &(L.matrix[Gkl][kl][0]), nrows,
                  &(I_OOOV.matrix[Gji][ji][la]), ncols, 1.0,
                  &(WBCA[Gbc][0][0]), ncols);
          }

          global_dpd_->sort_3d(WBCA, WABC, nirrep_, Gijk, I_OVVV.params->coltot, I_OVVV.params->colidx,
                 I_OVVV.params->colorb, I_OVVV.params->rsym, I_OVVV.params->ssym,
                 avir_off_, avir_off_, navirpi_, avir_off_, I_OVVV.params->colidx, cab, 1);

          for(Gab=0; Gab < nirrep_; Gab++) {
            Gc = Gab ^ Gijk;

            global_dpd_->free_dpd_block(WBCA[Gab], I_OVVV.params->coltot[Gab], navirpi_[Gc]);

            WACB[Gab] = global_dpd_->dpd_block_matrix(I_OVVV.params->coltot[Gab], navirpi_[Gc]);
          }

          for(Gd=0; Gd < nirrep_; Gd++) {
            /* +lambda_jkbd * <id||ac> */
            Gac = Gid = Gi ^ Gd;
            Gb = Gjk ^ Gd;

            bd = L.col_offset[Gjk][Gb];
            id = I_OVVV.row_offset[Gid][I];

            I_OVVV.matrix[Gid] = global_dpd_->dpd_block_matrix(navirpi_[Gd], I_OVVV.params->coltot[Gid]);
            global_dpd_->buf4_mat_irrep_rd_block(&I_OVVV, Gid, id, navirpi_[Gd]);

            nrows = I_OVVV.params->coltot[Gid];
            ncols = navirpi_[Gb];
            nlinks = navirpi_[Gd];

            if(nrows && ncols && nlinks)
              C_DGEMM('t', 't', nrows, ncols, nlinks, 1.0,
                  &(I_OVVV.matrix[Gid][0][0]), nrows,
                  &(L.matrix[Gjk][jk][bd]), nlinks, 1.0,
                  &(WACB[Gac][0][0]), ncols);

            global_dpd_->free_dpd_block(I_OVVV.matrix[Gid], navirpi_[Gd], I_OVVV.params->coltot[Gid]);

            /* -lambda_ikbd * <jd||ac> */
            Gac = Gjd = Gj ^ Gd;
            Gb = Gik ^ Gd;

            bd = L.col_offset[Gik][Gb];
            jd = I_OVVV.row_offset[Gjd][J];

            I_OVVV.matrix[Gjd] = global_dpd_->dpd_block_matrix(navirpi_[Gd], I_OVVV.params->coltot[Gjd]);
            global_dpd_->buf4_mat_irrep_rd_block(&I_OVVV, Gjd, jd, navirpi_[Gd]);

            nrows = I_OVVV.params->coltot[Gjd];
            ncols = navirpi_[Gb];
            nlinks = navirpi_[Gd];

            if(nrows && ncols && nlinks)
              C_DGEMM('t', 't', nrows, ncols, nlinks, -1.0,
                  &(I_OVVV.matrix[Gjd][0][0]), nrows,
                  &(L.matrix[Gik][ik][bd]), nlinks, 1.0,
                  &(WACB[Gac][0][0]), ncols);

            global_dpd_->free_dpd_block(I_OVVV.matrix[Gjd], navirpi_[Gd], I_OVVV.params->coltot[Gjd]);

            /* -lambda_jibd * <kd||ac> */
            Gac = Gkd = Gk ^ Gd;
            Gb = Gji ^ Gd;

            bd = L.col_offset[Gji][Gb];
            kd = I_OVVV.row_offset[Gkd][K];

            I_OVVV.matrix[Gkd] = global_dpd_->dpd_block_matrix(navirpi_[Gd], I_OVVV.params->coltot[Gkd]);
            global_dpd_->buf4_mat_irrep_rd_block(&I_OVVV, Gkd, kd, navirpi_[Gd]);

            nrows = I_OVVV.params->coltot[Gkd];
            ncols = navirpi_[Gb];
            nlinks = navirpi_[Gd];

            if(nrows && ncols && nlinks)
              C_DGEMM('t', 't', nrows, ncols, nlinks, -1.0,
                  &(I_OVVV.matrix[Gkd][0][0]), nrows,
                  &(L.matrix[Gji][ji][bd]), nlinks, 1.0,
                  &(WACB[Gac][0][0]), ncols);

            global_dpd_->free_dpd_block(I_OVVV.matrix[Gkd], navirpi_[Gd], I_OVVV.params->coltot[Gkd]);

          }

          for(Gl=0; Gl < nirrep_; Gl++) {
            /* +lambda_ilac * <jk||lb> */
            Gac = Gil = Gi ^ Gl;
            Gb = Gjk ^ Gl;

            lb = I_OOOV.col_offset[Gjk][Gl];
            il = L.row_offset[Gil][I];

            nrows = L.params->coltot[Gil];
            ncols = navirpi_[Gb];
            nlinks = naoccpi_[Gl];

            if(nrows && ncols && nlinks)
              C_DGEMM('t', 'n', nrows, ncols, nlinks, 1.0,
                  &(L.matrix[Gil][il][0]), nrows,
                  &(I_OOOV.matrix[Gjk][jk][lb]), ncols, 1.0,
                  &(WACB[Gac][0][0]), ncols);

            /* -lambda_jlac * <ik||lb> */
            Gac = Gjl = Gj ^ Gl;
            Gb = Gik ^ Gl;

            lb = I_OOOV.col_offset[Gik][Gl];
            jl = L.row_offset[Gjl][J];

            nrows = L.params->coltot[Gjl];
            ncols = navirpi_[Gb];
            nlinks = naoccpi_[Gl];

            if(nrows && ncols && nlinks)
              C_DGEMM('t', 'n', nrows, ncols, nlinks, -1.0,
                  &(L.matrix[Gjl][jl][0]), nrows,
                  &(I_OOOV.matrix[Gik][ik][lb]), ncols, 1.0,
                  &(WACB[Gac][0][0]), ncols);

            /* -lambda_klac * <ji||lb> */
            Gac = Gkl = Gk ^ Gl;
            Gb = Gji ^ Gl;

            lb = I_OOOV.col_offset[Gji][Gl];
            kl = L.row_offset[Gkl][K];

            nrows = L.params->coltot[Gkl];
            ncols = navirpi_[Gb];
            nlinks = naoccpi_[Gl];

            if(nrows && ncols && nlinks)
              C_DGEMM('t', 'n', nrows, ncols, nlinks, -1.0,
                  &(L.matrix[Gkl][kl][0]), nrows,
                  &(I_OOOV.matrix[Gji][ji][lb]), ncols, 1.0,
                  &(WACB[Gac][0][0]), ncols);
          }

          global_dpd_->sort_3d(WACB, WABC, nirrep_, Gijk, I_OVVV.params->coltot, I_OVVV.params->colidx,
                 I_OVVV.params->colorb, I_OVVV.params->rsym, I_OVVV.params->ssym,
                 avir_off_, avir_off_, navirpi_, avir_off_, I_OVVV.params->colidx, acb, 1);

          for(Gab=0; Gab < nirrep_; Gab++) {
            Gc = Gab ^ Gijk;

            global_dpd_->free_dpd_block(WACB[Gab], I_OVVV.params->coltot[Gab], navirpi_[Gc]);
          }

          /* Compute denominator and evaluate labda_ijkabc */
          for(Gab=0; Gab < nirrep_; Gab++) {
            Gc = Gab ^ Gijk;

            LABC[Gab] = global_dpd_->dpd_block_matrix(I_OVVV.params->coltot[Gab], navirpi_[Gc]);

            for(ab=0; ab < I_OVVV.params->coltot[Gab]; ab++) {
              A = I_OVVV.params->colorb[Gab][ab][0];
              Ga = I_OVVV.params->rsym[A];
              a = A - avir_off_[Ga];
              B = I_OVVV.params->colorb[Gab][ab][1];
              Gb = I_OVVV.params->ssym[B];
              b = B - avir_off_[Gb];

              Gbc = Gb ^ Gc;
              Gac = Ga ^ Gc;

              for(c=0; c < navirpi_[Gc]; c++) {
                C = avir_off_[Gc] + c;

                /* Copy W into labda_ijkabc */
                LABC[Gab][ab][c] = WABC[Gab][ab][c];

                /* Build the rest of the denominator and compute labda_ijkabc */
                denom = dijk;

                if(Ftilde_a_->rowspi()[Ga])
                    denom -= Ftilde_a_->get(Ga, a + naoccpi_[Ga], a + naoccpi_[Ga]);
                if(Ftilde_a_->rowspi()[Gb])
                    denom -= Ftilde_a_->get(Gb, b + naoccpi_[Gb], b + naoccpi_[Gb]);
                if(Ftilde_a_->rowspi()[Gc])
                    denom -= Ftilde_a_->get(Gc, c + naoccpi_[Gc], c + naoccpi_[Gc]);

                LABC[Gab][ab][c] /= denom;

              } /* c */
            } /* ab */
          } /* Gab */

          /* Compute the AAA energy contribution  */
          for(Gab=0; Gab < nirrep_; Gab++) {
            Gc = Gab ^ Gijk;
            ET_aaa += dot_block(LABC[Gab], WABC[Gab], I_OVVV.params->coltot[Gab], navirpi_[Gc], 1.0/6.0);
            global_dpd_->free_dpd_block(WABC[Gab], I_OVVV.params->coltot[Gab], navirpi_[Gc]);
            global_dpd_->free_dpd_block(LABC[Gab], I_OVVV.params->coltot[Gab], navirpi_[Gc]);
          }

            } /* I >= J >= K */

          } /* k */
        } /* j */
      } /* i */

        } /* Gk */
      } /* Gj */
    } /* Gi */

    free(WABC);
    free(LABC);
    free(WBCA);
    free(WACB);

    for(h=0; h < nirrep_; h++) {
      global_dpd_->buf4_mat_irrep_close(&L, h);
      global_dpd_->buf4_mat_irrep_close(&I_OOOV, h);
    }

    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&I_OVVV);
    global_dpd_->buf4_close(&I_OOOV);

    psio_->close(PSIF_LIBTRANS_DPD, 1);

    return ET_aaa;
}

double
DCFTSolver::compute_triples_aab()
{

    int h;
    int Gi, Gj, Gk, Ga, Gb, Gc, Gd, Gl;
    int Gji, Gij, Gjk, Gkj, Gik, Gki, Gijk;
    int Gab, Gbc, Gac, Gcb, Gca;
    int Gid, Gjd, Gkd;
    int Gil, Gjl, Gkl;
    int I, J, K, A, B, C;
    int i, j, k, a, b, c;
    int ij, ji, ik, ki, jk, kj;
    int ab;
    int dc, ad, bd;
    int lc, la, lb;
    int id, jd, kd;
    int il, jl, kl;
    double dijk, denom, ET_aab;
    int nrows, ncols, nlinks;
    dpdbuf4 L_AB, L_AA, L_BA;
    dpdbuf4 I_OVVV, I_OvVv, I_oVvV;
    dpdbuf4 I_OOOV, I_OoOv, I_oOoV;
    double ***WABc, ***WBcA, ***WAcB, ***WcAB, ***WcBA, ***LABc;

    psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

    global_dpd_->buf4_init(&L_AA, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                           ID("[O,O]"), ID("[V,V]"), 0, "Lambda <O'O'|V'V'>");
    global_dpd_->buf4_init(&L_AB, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                           ID("[O,o]"), ID("[V,v]"), 0, "Lambda <O'o'|V'v'>");
    global_dpd_->buf4_init(&L_BA, PSIF_DCFT_DPD, 0, ID("[o,O]"), ID("[v,V]"),
                           ID("[o,O]"), ID("[v,V]"), 0, "Lambda <o'O'|v'V'>");

    global_dpd_->buf4_init(&I_OVVV, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[V,V]"),
                           ID("[O,V]"), ID("[V,V]"), 0, "MO Ints <O'V'||V'V'>");
    global_dpd_->buf4_init(&I_OvVv, PSIF_LIBTRANS_DPD, 0, ID("[O,v]"), ID("[V,v]"),
                           ID("[O,v]"), ID("[V,v]"), 0, "MO Ints <O'v'|V'v'>");
    global_dpd_->buf4_init(&I_oVvV, PSIF_LIBTRANS_DPD, 0, ID("[o,V]"), ID("[v,V]"),
                           ID("[o,V]"), ID("[v,V]"), 0, "MO Ints <o'V'|v'V'>");

    global_dpd_->buf4_init(&I_OOOV, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[O,V]"),
                           ID("[O,O]"), ID("[O,V]"), 0, "MO Ints <O'O'||O'V'>");
    global_dpd_->buf4_init(&I_OoOv, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[O,v]"),
                           ID("[O,o]"), ID("[O,v]"), 0, "MO Ints <O'o'|O'v'>");
    global_dpd_->buf4_init(&I_oOoV, PSIF_LIBTRANS_DPD, 0, ID("[o,O]"), ID("[o,V]"),
                           ID("[o,O]"), ID("[o,V]"), 0, "MO Ints <o'O'|o'V'>");

    for(h=0; h < nirrep_; h++) {
      global_dpd_->buf4_mat_irrep_init(&L_AA, h);
      global_dpd_->buf4_mat_irrep_rd(&L_AA, h);

      global_dpd_->buf4_mat_irrep_init(&L_AB, h);
      global_dpd_->buf4_mat_irrep_rd(&L_AB, h);

      global_dpd_->buf4_mat_irrep_init(&L_BA, h);
      global_dpd_->buf4_mat_irrep_rd(&L_BA, h);

      global_dpd_->buf4_mat_irrep_init(&I_OOOV, h);
      global_dpd_->buf4_mat_irrep_rd(&I_OOOV, h);

      global_dpd_->buf4_mat_irrep_init(&I_OoOv, h);
      global_dpd_->buf4_mat_irrep_rd(&I_OoOv, h);

      global_dpd_->buf4_mat_irrep_init(&I_oOoV, h);
      global_dpd_->buf4_mat_irrep_rd(&I_oOoV, h);
    }

    WABc = (double ***) malloc(nirrep_ * sizeof(double **));
    WBcA = (double ***) malloc(nirrep_ * sizeof(double **));
    WAcB = (double ***) malloc(nirrep_ * sizeof(double **));
    WcAB = (double ***) malloc(nirrep_ * sizeof(double **));
    WcBA = (double ***) malloc(nirrep_ * sizeof(double **));
    LABc = (double ***) malloc(nirrep_ * sizeof(double **));

    ET_aab = 0.0;

    for(Gi=0; Gi < nirrep_; Gi++) {
      for(Gj=0; Gj < nirrep_; Gj++) {
        for(Gk=0; Gk < nirrep_; Gk++) {

      Gij = Gji = Gi ^ Gj;
      Gjk = Gkj = Gj ^ Gk;
      Gik = Gki = Gi ^ Gk;

      Gijk = Gi ^ Gj ^ Gk;

      for(i=0; i < naoccpi_[Gi]; i++) {
        I = aocc_off_[Gi] + i;
        for(j=0; j < naoccpi_[Gj]; j++) {
          J = aocc_off_[Gj] + j;
          for(k=0; k < nboccpi_[Gk]; k++) {
            K = bocc_off_[Gk] + k;

            if(I > J) {

          ij = I_OOOV.params->rowidx[I][J];
          ji = I_OOOV.params->rowidx[J][I];
          jk = I_OoOv.params->rowidx[J][K];
          kj = I_oOoV.params->rowidx[K][J];
          ik = I_OoOv.params->rowidx[I][K];
          ki = I_oOoV.params->rowidx[K][I];

          dijk = 0.0;
          if(Ftilde_a_->rowspi()[Gi])
            dijk += Ftilde_a_->get(Gi, i, i);
          if(Ftilde_a_->rowspi()[Gj])
            dijk += Ftilde_a_->get(Gj, j, j);
          if(Ftilde_b_->rowspi()[Gk])
            dijk += Ftilde_b_->get(Gk, k, k);

          /* Begin the W intermediate */

          for(Gab=0; Gab < nirrep_; Gab++) {
            Gc = Gab ^ Gijk;

            WABc[Gab] = global_dpd_->dpd_block_matrix(I_OVVV.params->coltot[Gab], nbvirpi_[Gc]);
          }

          for(Gd=0; Gd < nirrep_; Gd++) {
            /* +lambda_JkDc * <ID||AB> */
            Gab = Gid = Gi ^ Gd;
            Gc = Gjk ^ Gd;

            dc = L_AB.col_offset[Gjk][Gd];
            id = I_OVVV.row_offset[Gid][I];

            I_OVVV.matrix[Gid] = global_dpd_->dpd_block_matrix(navirpi_[Gd], I_OVVV.params->coltot[Gid]);
            global_dpd_->buf4_mat_irrep_rd_block(&I_OVVV, Gid, id, navirpi_[Gd]);

            nrows = I_OVVV.params->coltot[Gid];
            ncols = nbvirpi_[Gc];
            nlinks = navirpi_[Gd];

            if(nrows && ncols && nlinks)
              C_DGEMM('t', 'n', nrows, ncols, nlinks, 1.0,
                  &(I_OVVV.matrix[Gid][0][0]), nrows,
                  &(L_AB.matrix[Gjk][jk][dc]), ncols, 1.0,
                  &(WABc[Gab][0][0]), ncols);

            global_dpd_->free_dpd_block(I_OVVV.matrix[Gid], navirpi_[Gd], I_OVVV.params->coltot[Gid]);

            /* -lambda_IkDc * <JD||AB> */
            Gab = Gjd = Gj ^ Gd;
            Gc = Gik ^ Gd;

            dc = L_AB.col_offset[Gik][Gd];
            jd = I_OVVV.row_offset[Gjd][J];

            I_OVVV.matrix[Gjd] = global_dpd_->dpd_block_matrix(navirpi_[Gd], I_OVVV.params->coltot[Gjd]);
            global_dpd_->buf4_mat_irrep_rd_block(&I_OVVV, Gjd, jd, navirpi_[Gd]);

            nrows = I_OVVV.params->coltot[Gjd];
            ncols = nbvirpi_[Gc];
            nlinks = navirpi_[Gd];

            if(nrows && ncols && nlinks)
              C_DGEMM('t', 'n', nrows, ncols, nlinks, -1.0,
                  &(I_OVVV.matrix[Gjd][0][0]), nrows,
                  &(L_AB.matrix[Gik][ik][dc]), ncols, 1.0,
                  &(WABc[Gab][0][0]), ncols);

            global_dpd_->free_dpd_block(I_OVVV.matrix[Gjd], navirpi_[Gd], I_OVVV.params->coltot[Gjd]);

          }

          for(Gl=0; Gl < nirrep_; Gl++) {
            /* -lambda_ILAB * <Jk||Lc> */
            Gab = Gil = Gi ^ Gl;
            Gc = Gjk ^ Gl;

            lc = I_OoOv.col_offset[Gjk][Gl];
            il = L_AA.row_offset[Gil][I];

            nrows = L_AA.params->coltot[Gil];
            ncols = nbvirpi_[Gc];
            nlinks = naoccpi_[Gl];

            if(nrows && ncols && nlinks)
              C_DGEMM('t', 'n', nrows, ncols, nlinks, -1.0,
                  &(L_AA.matrix[Gil][il][0]), nrows,
                  &(I_OoOv.matrix[Gjk][jk][lc]), ncols, 1.0,
                  &(WABc[Gab][0][0]), ncols);


            /* +lambda_JLAB * <Ik||Lc> */
            Gab = Gjl = Gj ^ Gl;
            Gc = Gik ^ Gl;

            lc = I_OoOv.col_offset[Gik][Gl];
            jl = L_AA.row_offset[Gjl][J];

            nrows = L_AA.params->coltot[Gjl];
            ncols = nbvirpi_[Gc];
            nlinks = naoccpi_[Gl];

            if(nrows && ncols && nlinks)
              C_DGEMM('t', 'n', nrows, ncols, nlinks, 1.0,
                  &(L_AA.matrix[Gjl][jl][0]), nrows,
                  &(I_OoOv.matrix[Gik][ik][lc]), ncols, 1.0,
                  &(WABc[Gab][0][0]), ncols);
          }

          for(Gab=0; Gab < nirrep_; Gab++) {
            Gc = Gab ^ Gijk;

            WBcA[Gab] = global_dpd_->dpd_block_matrix(I_OvVv.params->coltot[Gab], navirpi_[Gc]);
          }

          for(Gd=0; Gd < nirrep_; Gd++) {

            /* -lambda_JkAd * <Id||Bc> */
            Gbc = Gid = Gi ^ Gd;
            Ga = Gjk ^ Gd;

            ad = L_AB.col_offset[Gjk][Ga];
            id = I_OvVv.row_offset[Gid][I];

            I_OvVv.matrix[Gid] = global_dpd_->dpd_block_matrix(nbvirpi_[Gd], I_OvVv.params->coltot[Gid]);
            global_dpd_->buf4_mat_irrep_rd_block(&I_OvVv, Gid, id, nbvirpi_[Gd]);

            nrows = I_OvVv.params->coltot[Gid];
            ncols = navirpi_[Ga];
            nlinks = nbvirpi_[Gd];

            if(nrows && ncols && nlinks)
              C_DGEMM('t', 't', nrows, ncols, nlinks, -1.0,
                  &(I_OvVv.matrix[Gid][0][0]), nrows,
                  &(L_AB.matrix[Gjk][jk][ad]), nlinks, 1.0,
                  &(WBcA[Gbc][0][0]), ncols);

            global_dpd_->free_dpd_block(I_OvVv.matrix[Gid], nbvirpi_[Gd], I_OvVv.params->coltot[Gid]);


            /* +lambda_IkAd * <Jd||Bc> */
            Gbc = Gjd = Gj ^ Gd;
            Ga = Gik ^ Gd;

            ad = L_AB.col_offset[Gik][Ga];
            jd = I_OvVv.row_offset[Gjd][J];

            I_OvVv.matrix[Gjd] = global_dpd_->dpd_block_matrix(nbvirpi_[Gd], I_OvVv.params->coltot[Gjd]);
            global_dpd_->buf4_mat_irrep_rd_block(&I_OvVv, Gjd, jd, nbvirpi_[Gd]);

            nrows = I_OvVv.params->coltot[Gjd];
            ncols = navirpi_[Ga];
            nlinks = nbvirpi_[Gd];

            if(nrows && ncols && nlinks)
              C_DGEMM('t', 't', nrows, ncols, nlinks, 1.0,
                  &(I_OvVv.matrix[Gjd][0][0]), nrows,
                  &(L_AB.matrix[Gik][ik][ad]), nlinks, 1.0,
                  &(WBcA[Gbc][0][0]), ncols);

            global_dpd_->free_dpd_block(I_OvVv.matrix[Gjd], nbvirpi_[Gd], I_OvVv.params->coltot[Gjd]);
          }

          for(Gl=0; Gl < nirrep_; Gl++) {

            /* +lambda_IlBc * <kJ||lA> */
            Gbc = Gil = Gi ^ Gl;
            Ga = Gkj ^ Gl;

            la = I_oOoV.col_offset[Gkj][Gl];
            il = L_AB.row_offset[Gil][I];

            nrows = L_AB.params->coltot[Gil];
            ncols = navirpi_[Ga];
            nlinks = nboccpi_[Gl];

            if(nrows && ncols && nlinks)
              C_DGEMM('t', 'n', nrows, ncols, nlinks, 1.0,
                  &(L_AB.matrix[Gil][il][0]), nrows,
                  &(I_oOoV.matrix[Gkj][kj][la]), ncols, 1.0,
                  &(WBcA[Gbc][0][0]), ncols);


            /* -lambda_JlBc * <kI||lA> */
            Gbc = Gjl = Gj ^ Gl;
            Ga = Gki ^ Gl;

            la = I_oOoV.col_offset[Gki][Gl];
            jl = L_AB.row_offset[Gjl][J];

            nrows = L_AB.params->coltot[Gjl];
            ncols = navirpi_[Ga];
            nlinks = nboccpi_[Gl];

            if(nrows && ncols && nlinks)
              C_DGEMM('t', 'n', nrows, ncols, nlinks, -1.0,
                  &(L_AB.matrix[Gjl][jl][0]), nrows,
                  &(I_oOoV.matrix[Gki][ki][la]), ncols, 1.0,
                  &(WBcA[Gbc][0][0]), ncols);

          }

          global_dpd_->sort_3d(WBcA, WABc, nirrep_, Gijk, I_OvVv.params->coltot, I_OvVv.params->colidx,
                 I_OvVv.params->colorb, I_OvVv.params->rsym, I_OvVv.params->ssym,
                 avir_off_, bvir_off_, navirpi_, avir_off_, I_OVVV.params->colidx, cab, 1);

          for(Gab=0; Gab < nirrep_; Gab++) {
            Gc = Gab ^ Gijk;

            global_dpd_->free_dpd_block(WBcA[Gab], I_OvVv.params->coltot[Gab], navirpi_[Gc]);

            WAcB[Gab] = global_dpd_->dpd_block_matrix(I_OvVv.params->coltot[Gab], navirpi_[Gc]);
          }

          for(Gd=0; Gd < nirrep_; Gd++) {

            /* +lambda_JkBd * <Id||Ac> */
            Gac = Gid = Gi ^ Gd;
            Gb = Gjk ^ Gd;

            bd = L_AB.col_offset[Gjk][Gb];
            id = I_OvVv.row_offset[Gid][I];

            I_OvVv.matrix[Gid] = global_dpd_->dpd_block_matrix(nbvirpi_[Gd], I_OvVv.params->coltot[Gid]);
            global_dpd_->buf4_mat_irrep_rd_block(&I_OvVv, Gid, id, nbvirpi_[Gd]);

            nrows = I_OvVv.params->coltot[Gid];
            ncols = navirpi_[Gb];
            nlinks = nbvirpi_[Gd];

            if(nrows && ncols && nlinks)
              C_DGEMM('t', 't', nrows, ncols, nlinks, 1.0,
                  &(I_OvVv.matrix[Gid][0][0]), nrows,
                  &(L_AB.matrix[Gjk][jk][bd]), nlinks, 1.0,
                  &(WAcB[Gac][0][0]), ncols);

            global_dpd_->free_dpd_block(I_OvVv.matrix[Gid], nbvirpi_[Gd], I_OvVv.params->coltot[Gid]);


            /* -lambda_IkBd * <Jd||Ac> */
            Gac = Gjd = Gj ^ Gd;
            Gb = Gik ^ Gd;

            bd = L_AB.col_offset[Gik][Gb];
            jd = I_OvVv.row_offset[Gjd][J];

            I_OvVv.matrix[Gjd] = global_dpd_->dpd_block_matrix(nbvirpi_[Gd], I_OvVv.params->coltot[Gjd]);
            global_dpd_->buf4_mat_irrep_rd_block(&I_OvVv, Gjd, jd, nbvirpi_[Gd]);

            nrows = I_OvVv.params->coltot[Gjd];
            ncols = navirpi_[Gb];
            nlinks = nbvirpi_[Gd];

            if(nrows && ncols && nlinks)
              C_DGEMM('t', 't', nrows, ncols, nlinks, -1.0,
                  &(I_OvVv.matrix[Gjd][0][0]), nrows,
                  &(L_AB.matrix[Gik][ik][bd]), nlinks, 1.0,
                  &(WAcB[Gac][0][0]), ncols);

            global_dpd_->free_dpd_block(I_OvVv.matrix[Gjd], nbvirpi_[Gd], I_OvVv.params->coltot[Gjd]);

          }

          for(Gl=0; Gl < nirrep_; Gl++) {

            /* -lambda_IlAc * <kJ||lB> */
            Gac = Gil = Gi ^ Gl;
            Gb = Gkj ^ Gl;

            lb = I_oOoV.col_offset[Gkj][Gl];
            il = L_AB.row_offset[Gil][I];

            nrows = L_AB.params->coltot[Gil];
            ncols = navirpi_[Gb];
            nlinks = nboccpi_[Gl];

            if(nrows && ncols && nlinks)
              C_DGEMM('t', 'n', nrows, ncols, nlinks, -1.0,
                  &(L_AB.matrix[Gil][il][0]), nrows,
                  &(I_oOoV.matrix[Gkj][kj][lb]), ncols, 1.0,
                  &(WAcB[Gac][0][0]), ncols);


            /* +lambda_JlAc * <kI||lB> */
            Gac = Gjl = Gj ^ Gl;
            Gb = Gki ^ Gl;

            lb = I_oOoV.col_offset[Gki][Gl];
            jl = L_AB.row_offset[Gjl][J];

            nrows = L_AB.params->coltot[Gjl];
            ncols = navirpi_[Gb];
            nlinks = nboccpi_[Gl];

            if(nrows && ncols && nlinks)
              C_DGEMM('t', 'n', nrows, ncols, nlinks, 1.0,
                  &(L_AB.matrix[Gjl][jl][0]), nrows,
                  &(I_oOoV.matrix[Gki][ki][lb]), ncols, 1.0,
                  &(WAcB[Gac][0][0]), ncols);
          }

          global_dpd_->sort_3d(WAcB, WABc, nirrep_, Gijk, I_OvVv.params->coltot, I_OvVv.params->colidx,
                 I_OvVv.params->colorb, I_OvVv.params->rsym, I_OvVv.params->ssym,
                 avir_off_, bvir_off_, navirpi_, avir_off_, I_OVVV.params->colidx, acb, 1);

          for(Gab=0; Gab < nirrep_; Gab++) {
            Gc = Gab ^ Gijk;

            global_dpd_->free_dpd_block(WAcB[Gab], I_OvVv.params->coltot[Gab], navirpi_[Gc]);

            WcBA[Gab] = global_dpd_->dpd_block_matrix(I_oVvV.params->coltot[Gab], navirpi_[Gc]);
          }

          for(Gd=0; Gd < nirrep_; Gd++) {

            /* -lambda_JIAD * <kD||cB> */
            Gcb = Gkd = Gk ^ Gd;
            Ga = Gji ^ Gd;

            ad = L_AA.col_offset[Gji][Ga];
            kd = I_oVvV.row_offset[Gkd][K];

            I_oVvV.matrix[Gkd] = global_dpd_->dpd_block_matrix(navirpi_[Gd], I_oVvV.params->coltot[Gkd]);
            global_dpd_->buf4_mat_irrep_rd_block(&I_oVvV, Gkd, kd, navirpi_[Gd]);

            nrows = I_oVvV.params->coltot[Gkd];
            ncols = navirpi_[Ga];
            nlinks = navirpi_[Gd];

            if(nrows && ncols && nlinks)
              C_DGEMM('t', 't', nrows, ncols, nlinks, -1.0,
                  &(I_oVvV.matrix[Gkd][0][0]), nrows,
                  &(L_AA.matrix[Gji][ji][ad]), nlinks, 1.0,
                  &(WcBA[Gcb][0][0]), ncols);

            global_dpd_->free_dpd_block(I_oVvV.matrix[Gkd], navirpi_[Gd], I_oVvV.params->coltot[Gkd]);
          }

          for(Gl=0; Gl < nirrep_; Gl++) {

            /* -lambda_kLcB * <JI||LA> */
            Gcb = Gkl = Gk ^ Gl;
            Ga = Gji ^ Gl;

            la = I_OOOV.col_offset[Gji][Gl];
            kl = L_BA.row_offset[Gkl][K];

            nrows = L_BA.params->coltot[Gkl];
            ncols = navirpi_[Ga];
            nlinks = naoccpi_[Gl];

            if(nrows && ncols && nlinks)
              C_DGEMM('t', 'n', nrows, ncols, nlinks, -1.0,
                  &(L_BA.matrix[Gkl][kl][0]), nrows,
                  &(I_OOOV.matrix[Gji][ji][la]), ncols, 1.0,
                  &(WcBA[Gcb][0][0]), ncols);
          }

          global_dpd_->sort_3d(WcBA, WABc, nirrep_, Gijk, I_oVvV.params->coltot, I_oVvV.params->colidx,
                 I_oVvV.params->colorb, I_oVvV.params->rsym, I_oVvV.params->ssym,
                 bvir_off_, avir_off_, navirpi_, avir_off_, I_OVVV.params->colidx, cba, 1);

          for(Gab=0; Gab < nirrep_; Gab++) {
            Gc = Gab ^ Gijk;

            global_dpd_->free_dpd_block(WcBA[Gab], I_oVvV.params->coltot[Gab], navirpi_[Gc]);

            WcAB[Gab] = global_dpd_->dpd_block_matrix(I_oVvV.params->coltot[Gab], navirpi_[Gc]);
          }

          for(Gd=0; Gd < nirrep_; Gd++) {

            /* +lambda_JIBD * <kD||cA> */
            Gca = Gkd = Gk ^ Gd;
            Gb = Gji ^ Gd;

            bd = L_AA.col_offset[Gji][Gb];
            kd = I_oVvV.row_offset[Gkd][K];

            I_oVvV.matrix[Gkd] = global_dpd_->dpd_block_matrix(navirpi_[Gd], I_oVvV.params->coltot[Gkd]);
            global_dpd_->buf4_mat_irrep_rd_block(&I_oVvV, Gkd, kd, navirpi_[Gd]);

            nrows = I_oVvV.params->coltot[Gkd];
            ncols = navirpi_[Gb];
            nlinks = navirpi_[Gd];

            if(nrows && ncols && nlinks)
              C_DGEMM('t', 't', nrows, ncols, nlinks, 1.0,
                  &(I_oVvV.matrix[Gkd][0][0]), nrows,
                  &(L_AA.matrix[Gji][ji][bd]), nlinks, 1.0,
                  &(WcAB[Gca][0][0]), ncols);

            global_dpd_->free_dpd_block(I_oVvV.matrix[Gkd], navirpi_[Gd], I_oVvV.params->coltot[Gkd]);
          }

          for(Gl=0; Gl < nirrep_; Gl++) {

            /* lambda_kLcA * <JI||LB> */
            Gca = Gkl = Gk ^ Gl;
            Gb = Gji ^ Gl;

            lb = I_OOOV.col_offset[Gji][Gl];
            kl = L_BA.row_offset[Gkl][K];

            nrows = L_BA.params->coltot[Gkl];
            ncols = navirpi_[Gb];
            nlinks = naoccpi_[Gl];

            if(nrows && ncols && nlinks)
              C_DGEMM('t', 'n', nrows, ncols, nlinks, 1.0,
                  &(L_BA.matrix[Gkl][kl][0]), nrows,
                  &(I_OOOV.matrix[Gji][ji][lb]), ncols, 1.0,
                  &(WcAB[Gca][0][0]), ncols);
          }

          global_dpd_->sort_3d(WcAB, WABc, nirrep_, Gijk, I_oVvV.params->coltot, I_oVvV.params->colidx,
                 I_oVvV.params->colorb, I_oVvV.params->rsym, I_oVvV.params->ssym,
                 bvir_off_, avir_off_, navirpi_, avir_off_, I_OVVV.params->colidx, bca, 1);


          for(Gab=0; Gab < nirrep_; Gab++) {
            Gc = Gab ^ Gijk;

            global_dpd_->free_dpd_block(WcAB[Gab], I_oVvV.params->coltot[Gab], navirpi_[Gc]);

            LABc[Gab] = global_dpd_->dpd_block_matrix(I_OVVV.params->coltot[Gab], nbvirpi_[Gc]);
          }

          /* Compute denominator and evaluate labda_ijkabc */
          for(Gab=0; Gab < nirrep_; Gab++) {
            Gc = Gab ^ Gijk;

            for(ab=0; ab < I_OVVV.params->coltot[Gab]; ab++) {
              A = I_OVVV.params->colorb[Gab][ab][0];
              Ga = I_OVVV.params->rsym[A];
              a = A - avir_off_[Ga];
              B = I_OVVV.params->colorb[Gab][ab][1];
              Gb = I_OVVV.params->ssym[B];
              b = B - avir_off_[Gb];

              Gbc = Gb ^ Gc;
              Gac = Ga ^ Gc;

              for(c=0; c < nbvirpi_[Gc]; c++) {
                C = bvir_off_[Gc] + c;

                /* Copy W into labda_ijkabc */
                LABc[Gab][ab][c] = WABc[Gab][ab][c];

                /* Build the rest of the denominator and compute labda_ijkabc */
                denom = dijk;
                if(Ftilde_a_->rowspi()[Ga])
                    denom -= Ftilde_a_->get(Ga, a + naoccpi_[Ga], a + naoccpi_[Ga]);
                if(Ftilde_a_->rowspi()[Gb])
                    denom -= Ftilde_a_->get(Gb, b + naoccpi_[Gb], b + naoccpi_[Gb]);
                if(Ftilde_b_->rowspi()[Gc])
                    denom -= Ftilde_b_->get(Gc, c + nboccpi_[Gc], c + nboccpi_[Gc]);

                LABc[Gab][ab][c] /= denom;

              } /* c */
            } /* ab */
          } /* Gab */

          /* Compute the AAB energy contribution  */
          for(Gab=0; Gab < nirrep_; Gab++) {
            Gc = Gab ^ Gijk;
            ET_aab += dot_block(LABc[Gab], WABc[Gab], I_OVVV.params->coltot[Gab], nbvirpi_[Gc], 0.5);
          }

          for(Gab=0; Gab < nirrep_; Gab++) {
            Gc = Gab ^ Gijk;
            global_dpd_->free_dpd_block(WABc[Gab], I_OVVV.params->coltot[Gab], nbvirpi_[Gc]);
            global_dpd_->free_dpd_block(LABc[Gab], I_OVVV.params->coltot[Gab], nbvirpi_[Gc]);
          }

            } /* I >= J */

          } /* k */
        } /* j */
      } /* i */

        } /* Gk */
      } /* Gj */
    } /* Gi */

    free(WABc);
    free(WBcA);
    free(WAcB);
    free(WcAB);
    free(WcBA);
    free(LABc);

    for(h=0; h < nirrep_; h++) {
      global_dpd_->buf4_mat_irrep_close(&L_AA, h);
      global_dpd_->buf4_mat_irrep_close(&L_AB, h);
      global_dpd_->buf4_mat_irrep_close(&L_BA, h);
      global_dpd_->buf4_mat_irrep_close(&I_OOOV, h);
      global_dpd_->buf4_mat_irrep_close(&I_OoOv, h);
      global_dpd_->buf4_mat_irrep_close(&I_oOoV, h);
    }

    global_dpd_->buf4_close(&L_AA);
    global_dpd_->buf4_close(&L_AB);
    global_dpd_->buf4_close(&L_BA);
    global_dpd_->buf4_close(&I_OVVV);
    global_dpd_->buf4_close(&I_OvVv);
    global_dpd_->buf4_close(&I_oVvV);
    global_dpd_->buf4_close(&I_OOOV);
    global_dpd_->buf4_close(&I_OoOv);
    global_dpd_->buf4_close(&I_oOoV);

    psio_->close(PSIF_LIBTRANS_DPD, 1);

    return ET_aab;

}

double
DCFTSolver::compute_triples_abb()
{

    int h;
    int Gi, Gj, Gk, Ga, Gb, Gc, Gd, Gl;
    int Gji, Gij, Gjk, Gkj, Gik, Gki, Gijk;
    int Gab, Gbc, Gac, Gca, Gba;
    int Gid, Gjd, Gkd;
    int Gil, Gjl, Gkl;
    int I, J, K, A, B, C;
    int i, j, k, a, b, c;
    int ij, ji, ik, ki, jk, kj;
    int ab;
    int cd, bd, ad, db, dc;
    int lc, lb, la;
    int id, jd, kd;
    int il, jl, kl;
    double dijk, denom, ET_abb;
    int nrows, ncols, nlinks;
    dpdbuf4 L_AB, L_BB, L_BA;
    dpdbuf4 I_ovvv, I_OvVv, I_oVvV;
    dpdbuf4 I_ooov, I_OoOv, I_oOoV;
    double ***WAbc, ***WAcb, ***WbAc, ***WcAb, ***WbcA, ***LAbc;

    psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

    global_dpd_->buf4_init(&L_BB, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                           ID("[o,o]"), ID("[v,v]"), 0, "Lambda <o'o'|v'v'>");
    global_dpd_->buf4_init(&L_AB, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                           ID("[O,o]"), ID("[V,v]"), 0, "Lambda <O'o'|V'v'>");
    global_dpd_->buf4_init(&L_BA, PSIF_DCFT_DPD, 0, ID("[o,O]"), ID("[v,V]"),
                           ID("[o,O]"), ID("[v,V]"), 0, "Lambda <o'O'|v'V'>");

    global_dpd_->buf4_init(&I_ovvv, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[v,v]"),
                           ID("[o,v]"), ID("[v,v]"), 0, "MO Ints <o'v'||v'v'>");
    global_dpd_->buf4_init(&I_OvVv, PSIF_LIBTRANS_DPD, 0, ID("[O,v]"), ID("[V,v]"),
                           ID("[O,v]"), ID("[V,v]"), 0, "MO Ints <O'v'|V'v'>");
    global_dpd_->buf4_init(&I_oVvV, PSIF_LIBTRANS_DPD, 0, ID("[o,V]"), ID("[v,V]"),
                           ID("[o,V]"), ID("[v,V]"), 0, "MO Ints <o'V'|v'V'>");

    global_dpd_->buf4_init(&I_ooov, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[o,v]"),
                           ID("[o,o]"), ID("[o,v]"), 0, "MO Ints <o'o'||o'v'>");
    global_dpd_->buf4_init(&I_OoOv, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[O,v]"),
                           ID("[O,o]"), ID("[O,v]"), 0, "MO Ints <O'o'|O'v'>");
    global_dpd_->buf4_init(&I_oOoV, PSIF_LIBTRANS_DPD, 0, ID("[o,O]"), ID("[o,V]"),
                           ID("[o,O]"), ID("[o,V]"), 0, "MO Ints <o'O'|o'V'>");

    for(h=0; h < nirrep_; h++) {
      global_dpd_->buf4_mat_irrep_init(&L_BB, h);
      global_dpd_->buf4_mat_irrep_rd(&L_BB, h);

      global_dpd_->buf4_mat_irrep_init(&L_AB, h);
      global_dpd_->buf4_mat_irrep_rd(&L_AB, h);

      global_dpd_->buf4_mat_irrep_init(&L_BA, h);
      global_dpd_->buf4_mat_irrep_rd(&L_BA, h);

      global_dpd_->buf4_mat_irrep_init(&I_ooov, h);
      global_dpd_->buf4_mat_irrep_rd(&I_ooov, h);

      global_dpd_->buf4_mat_irrep_init(&I_OoOv, h);
      global_dpd_->buf4_mat_irrep_rd(&I_OoOv, h);

      global_dpd_->buf4_mat_irrep_init(&I_oOoV, h);
      global_dpd_->buf4_mat_irrep_rd(&I_oOoV, h);
    }

    ET_abb = 0.0;

    WAbc = (double ***) malloc(nirrep_ * sizeof(double **));
    LAbc = (double ***) malloc(nirrep_ * sizeof(double **));
    WAcb = (double ***) malloc(nirrep_ * sizeof(double **));
    WbcA = (double ***) malloc(nirrep_ * sizeof(double **));
    WcAb = (double ***) malloc(nirrep_ * sizeof(double **));
    WbAc = (double ***) malloc(nirrep_ * sizeof(double **));

    for(Gi=0; Gi < nirrep_; Gi++) {
      for(Gj=0; Gj < nirrep_; Gj++) {
        for(Gk=0; Gk < nirrep_; Gk++) {

      Gij = Gji = Gi ^ Gj;
      Gjk = Gkj = Gj ^ Gk;
      Gik = Gki = Gi ^ Gk;

      Gijk = Gi ^ Gj ^ Gk;

      for(i=0; i < naoccpi_[Gi]; i++) {
        I = aocc_off_[Gi] + i;
        for(j=0; j < nboccpi_[Gj]; j++) {
          J = bocc_off_[Gj] + j;
          for(k=0; k < nboccpi_[Gk]; k++) {
            K = bocc_off_[Gk] + k;

            if(J > K) {

          ij = I_OoOv.params->rowidx[I][J];
          ji = I_oOoV.params->rowidx[J][I];
          jk = I_ooov.params->rowidx[J][K];
          kj = I_ooov.params->rowidx[K][J];
          ik = I_OoOv.params->rowidx[I][K];
          ki = I_oOoV.params->rowidx[K][I];

          dijk = 0.0;
          if(Ftilde_a_->rowspi()[Gi])
            dijk += Ftilde_a_->get(Gi, i, i);
          if(Ftilde_b_->rowspi()[Gj])
            dijk += Ftilde_b_->get(Gj, j, j);
          if(Ftilde_b_->rowspi()[Gk])
            dijk += Ftilde_b_->get(Gk, k, k);

          /* Begin the W intermediate */

          for(Gab=0; Gab < nirrep_; Gab++) {
            Gc = Gab ^ Gijk;

            WAbc[Gab] = global_dpd_->dpd_block_matrix(I_OvVv.params->coltot[Gab], nbvirpi_[Gc]);
          }

          for(Gd=0; Gd < nirrep_; Gd++) {
            /* -lambda_jkcd * <Id||Ab> */
            Gab = Gid = Gi ^ Gd;
            Gc = Gjk ^ Gd;

            cd = L_BB.col_offset[Gjk][Gc];
            id = I_OvVv.row_offset[Gid][I];

            I_OvVv.matrix[Gid] = global_dpd_->dpd_block_matrix(nbvirpi_[Gd], I_OvVv.params->coltot[Gid]);
            global_dpd_->buf4_mat_irrep_rd_block(&I_OvVv, Gid, id, nbvirpi_[Gd]);

            nrows = I_OvVv.params->coltot[Gid];
            ncols = nbvirpi_[Gc];
            nlinks = nbvirpi_[Gd];

            if(nrows && ncols && nlinks)
              C_DGEMM('t', 't', nrows, ncols, nlinks, -1.0,
                  &(I_OvVv.matrix[Gid][0][0]), nrows,
                  &(L_BB.matrix[Gjk][jk][cd]), nlinks, 1.0,
                  &(WAbc[Gab][0][0]), ncols);

            global_dpd_->free_dpd_block(I_OvVv.matrix[Gid], nbvirpi_[Gd], I_OvVv.params->coltot[Gid]);

          }

          for(Gl=0; Gl < nirrep_; Gl++) {
            /* -lambda_IlAb <jk||lc> */
            Gab = Gil = Gi ^ Gl;
            Gc = Gjk ^ Gl;

            lc = I_ooov.col_offset[Gjk][Gl];
            il = L_AB.row_offset[Gil][I];

            nrows = L_AB.params->coltot[Gil];
            ncols = nbvirpi_[Gc];
            nlinks = nboccpi_[Gl];

            if(nrows && ncols && nlinks)
              C_DGEMM('t', 'n', nrows, ncols, nlinks, -1.0,
                  &(L_AB.matrix[Gil][il][0]), nrows,
                  &(I_ooov.matrix[Gjk][jk][lc]), ncols, 1.0,
                  &(WAbc[Gab][0][0]), ncols);
          }

          for(Gab=0; Gab < nirrep_; Gab++) {
            Gc = Gab ^ Gijk;

            WAcb[Gab] = global_dpd_->dpd_block_matrix(I_OvVv.params->coltot[Gab], nbvirpi_[Gc]);
          }

          for(Gd=0; Gd < nirrep_; Gd++) {
            /* +lambda_jkbd * <Id||Ac> */
            Gac = Gid = Gi ^ Gd;
            Gb = Gjk ^ Gd;

            bd = L_BB.col_offset[Gjk][Gb];
            id = I_OvVv.row_offset[Gid][I];

            I_OvVv.matrix[Gid] = global_dpd_->dpd_block_matrix(nbvirpi_[Gd], I_OvVv.params->coltot[Gid]);
            global_dpd_->buf4_mat_irrep_rd_block(&I_OvVv, Gid, id, nbvirpi_[Gd]);

            nrows = I_OvVv.params->coltot[Gid];
            ncols = nbvirpi_[Gb];
            nlinks = nbvirpi_[Gd];

            if(nrows && ncols && nlinks)
              C_DGEMM('t', 't', nrows, ncols, nlinks, 1.0,
                  &(I_OvVv.matrix[Gid][0][0]), nrows,
                  &(L_BB.matrix[Gjk][jk][bd]), nlinks, 1.0,
                  &(WAcb[Gac][0][0]), ncols);

            global_dpd_->free_dpd_block(I_OvVv.matrix[Gid], nbvirpi_[Gd], I_OvVv.params->coltot[Gid]);

          }

          for(Gl=0; Gl < nirrep_; Gl++) {
            /* +lambda_IlAc <jk||lb> */
            Gac = Gil = Gi ^ Gl;
            Gb = Gjk ^ Gl;

            lb = I_ooov.col_offset[Gjk][Gl];
            il = L_AB.row_offset[Gil][I];

            nrows = L_AB.params->coltot[Gil];
            ncols = nbvirpi_[Gb];
            nlinks = nboccpi_[Gl];

            if(nrows && ncols && nlinks)
              C_DGEMM('t', 'n', nrows, ncols, nlinks, 1.0,
                  &(L_AB.matrix[Gil][il][0]), nrows,
                  &(I_ooov.matrix[Gjk][jk][lb]), ncols, 1.0,
                  &(WAcb[Gac][0][0]), ncols);
          }

          global_dpd_->sort_3d(WAcb, WAbc, nirrep_, Gijk, I_OvVv.params->coltot, I_OvVv.params->colidx,
                 I_OvVv.params->colorb, I_OvVv.params->rsym, I_OvVv.params->ssym,
                 avir_off_, bvir_off_, nbvirpi_, bvir_off_, I_OvVv.params->colidx, acb, 1);

          for(Gab=0; Gab < nirrep_; Gab++) {
            Gc = Gab ^ Gijk;

            global_dpd_->free_dpd_block(WAcb[Gab], I_OvVv.params->coltot[Gab], nbvirpi_[Gc]);

            WbcA[Gab] = global_dpd_->dpd_block_matrix(I_ovvv.params->coltot[Gab], navirpi_[Gc]);
          }

          for(Gd=0; Gd < nirrep_; Gd++) {
            /* +lambda_IkAd * <jd||bc> */
            Gbc = Gjd = Gj ^ Gd;
            Ga = Gik ^ Gd;

            ad = L_AB.col_offset[Gik][Ga];
            jd = I_ovvv.row_offset[Gjd][J];

            I_ovvv.matrix[Gjd] = global_dpd_->dpd_block_matrix(nbvirpi_[Gd], I_ovvv.params->coltot[Gjd]);
            global_dpd_->buf4_mat_irrep_rd_block(&I_ovvv, Gjd, jd, nbvirpi_[Gd]);

            nrows = I_ovvv.params->coltot[Gjd];
            ncols = navirpi_[Ga];
            nlinks = nbvirpi_[Gd];

            if(nrows && ncols && nlinks)
              C_DGEMM('t', 't', nrows, ncols, nlinks, 1.0,
                  &(I_ovvv.matrix[Gjd][0][0]), nrows,
                  &(L_AB.matrix[Gik][ik][ad]), nlinks, 1.0,
                  &(WbcA[Gbc][0][0]), ncols);

            global_dpd_->free_dpd_block(I_ovvv.matrix[Gjd], nbvirpi_[Gd], I_ovvv.params->coltot[Gjd]);

            /* -lambda_IjAd * <kd||bc> */
            Gbc = Gkd = Gk ^ Gd;
            Ga = Gij ^ Gd;

            ad = L_AB.col_offset[Gij][Ga];
            kd = I_ovvv.row_offset[Gkd][K];

            I_ovvv.matrix[Gkd] = global_dpd_->dpd_block_matrix(nbvirpi_[Gd], I_ovvv.params->coltot[Gkd]);
            global_dpd_->buf4_mat_irrep_rd_block(&I_ovvv, Gkd, kd, nbvirpi_[Gd]);

            nrows = I_ovvv.params->coltot[Gkd];
            ncols = navirpi_[Ga];
            nlinks = nbvirpi_[Gd];

            if(nrows && ncols && nlinks)
              C_DGEMM('t', 't', nrows, ncols, nlinks, -1.0,
                  &(I_ovvv.matrix[Gkd][0][0]), nrows,
                  &(L_AB.matrix[Gij][ij][ad]), nlinks, 1.0,
                  &(WbcA[Gbc][0][0]), ncols);

            global_dpd_->free_dpd_block(I_ovvv.matrix[Gkd], nbvirpi_[Gd], I_ovvv.params->coltot[Gkd]);

          }

          for(Gl=0; Gl < nirrep_; Gl++) {
            /* -lambda_jlbc <kI||lA> */
            Gbc = Gjl = Gj ^ Gl;
            Ga = Gki ^ Gl;

            la = I_oOoV.col_offset[Gki][Gl];
            jl = L_BB.row_offset[Gjl][J];

            nrows = L_BB.params->coltot[Gjl];
            ncols = navirpi_[Ga];
            nlinks = nboccpi_[Gl];

            if(nrows && ncols && nlinks)
              C_DGEMM('t', 'n', nrows, ncols, nlinks, -1.0,
                  &(L_BB.matrix[Gjl][jl][0]), nrows,
                  &(I_oOoV.matrix[Gki][ki][la]), ncols, 1.0,
                  &(WbcA[Gbc][0][0]), ncols);

            /* +lambda_klbc <jI||lA> */
            Gbc = Gkl = Gk ^ Gl;
            Ga = Gji ^ Gl;

            la = I_oOoV.col_offset[Gji][Gl];
            kl = L_BB.row_offset[Gkl][K];

            nrows = L_BB.params->coltot[Gkl];
            ncols = navirpi_[Ga];
            nlinks = nboccpi_[Gl];

            if(nrows && ncols && nlinks)
              C_DGEMM('t', 'n', nrows, ncols, nlinks, 1.0,
                  &(L_BB.matrix[Gkl][kl][0]), nrows,
                  &(I_oOoV.matrix[Gji][ji][la]), ncols, 1.0,
                  &(WbcA[Gbc][0][0]), ncols);
          }

          global_dpd_->sort_3d(WbcA, WAbc, nirrep_, Gijk, I_ovvv.params->coltot, I_ovvv.params->colidx,
                 I_ovvv.params->colorb, I_ovvv.params->rsym, I_ovvv.params->ssym,
                 bvir_off_, bvir_off_, navirpi_, avir_off_, I_OvVv.params->colidx, cab, 1);

          for(Gab=0; Gab < nirrep_; Gab++) {
            Gc = Gab ^ Gijk;

            global_dpd_->free_dpd_block(WbcA[Gab], I_ovvv.params->coltot[Gab], navirpi_[Gc]);

            WcAb[Gab] = global_dpd_->dpd_block_matrix(I_oVvV.params->coltot[Gab], nbvirpi_[Gc]);
          }

          for(Gd=0; Gd < nirrep_; Gd++) {
            /* -lambda_IkDb * <jD||cA> */
            Gca = Gjd = Gj ^ Gd;
            Gb = Gik ^ Gd;

            db = L_AB.col_offset[Gik][Gd];
            jd = I_oVvV.row_offset[Gjd][J];

            I_oVvV.matrix[Gjd] = global_dpd_->dpd_block_matrix(navirpi_[Gd], I_oVvV.params->coltot[Gjd]);
            global_dpd_->buf4_mat_irrep_rd_block(&I_oVvV, Gjd, jd, navirpi_[Gd]);

            nrows = I_oVvV.params->coltot[Gjd];
            ncols = nbvirpi_[Gb];
            nlinks = navirpi_[Gd];

            if(nrows && ncols && nlinks)
              C_DGEMM('t', 'n', nrows, ncols, nlinks, -1.0,
                  &(I_oVvV.matrix[Gjd][0][0]), nrows,
                  &(L_AB.matrix[Gik][ik][db]), ncols, 1.0,
                  &(WcAb[Gca][0][0]), ncols);

            global_dpd_->free_dpd_block(I_oVvV.matrix[Gjd], navirpi_[Gd], I_oVvV.params->coltot[Gjd]);

            /* +lambda_IjDb * <kD||cA> */
            Gca = Gkd = Gk ^ Gd;
            Gb = Gij ^ Gd;

            db = L_AB.col_offset[Gij][Gd];
            kd = I_oVvV.row_offset[Gkd][K];

            I_oVvV.matrix[Gkd] = global_dpd_->dpd_block_matrix(navirpi_[Gd], I_oVvV.params->coltot[Gkd]);
            global_dpd_->buf4_mat_irrep_rd_block(&I_oVvV, Gkd, kd, navirpi_[Gd]);

            nrows = I_oVvV.params->coltot[Gkd];
            ncols = nbvirpi_[Gb];
            nlinks = navirpi_[Gd];

            if(nrows && ncols && nlinks)
              C_DGEMM('t', 'n', nrows, ncols, nlinks, 1.0,
                  &(I_oVvV.matrix[Gkd][0][0]), nrows,
                  &(L_AB.matrix[Gij][ij][db]), ncols, 1.0,
                  &(WcAb[Gca][0][0]), ncols);

            global_dpd_->free_dpd_block(I_oVvV.matrix[Gkd], navirpi_[Gd], I_oVvV.params->coltot[Gkd]);

          }

          for(Gl=0; Gl < nirrep_; Gl++) {
            /* +lambda_jLcA <Ik||Lb> */
            Gca = Gjl = Gj ^ Gl;
            Gb = Gik ^ Gl;

            lb = I_OoOv.col_offset[Gik][Gl];
            jl = L_BA.row_offset[Gjl][J];

            nrows = L_BA.params->coltot[Gjl];
            ncols = nbvirpi_[Gb];
            nlinks = naoccpi_[Gl];

            if(nrows && ncols && nlinks)
              C_DGEMM('t', 'n', nrows, ncols, nlinks, 1.0,
                  &(L_BA.matrix[Gjl][jl][0]), nrows,
                  &(I_OoOv.matrix[Gik][ik][lb]), ncols, 1.0,
                  &(WcAb[Gca][0][0]), ncols);

            /* -lambda_kLcA <Ij||Lb> */
            Gca = Gkl = Gk ^ Gl;
            Gb = Gij ^ Gl;

            lb = I_OoOv.col_offset[Gij][Gl];
            kl = L_BA.row_offset[Gkl][K];

            nrows = L_BA.params->coltot[Gkl];
            ncols = nbvirpi_[Gb];
            nlinks = naoccpi_[Gl];

            if(nrows && ncols && nlinks)
              C_DGEMM('t', 'n', nrows, ncols, nlinks, -1.0,
                  &(L_BA.matrix[Gkl][kl][0]), nrows,
                  &(I_OoOv.matrix[Gij][ij][lb]), ncols, 1.0,
                  &(WcAb[Gca][0][0]), ncols);
          }

          global_dpd_->sort_3d(WcAb, WAbc, nirrep_, Gijk, I_oVvV.params->coltot, I_oVvV.params->colidx,
                 I_oVvV.params->colorb, I_oVvV.params->rsym, I_oVvV.params->ssym,
                 bvir_off_, avir_off_, nbvirpi_, bvir_off_, I_OvVv.params->colidx, bca, 1);

          for(Gab=0; Gab < nirrep_; Gab++) {
            Gc = Gab ^ Gijk;

            global_dpd_->free_dpd_block(WcAb[Gab], I_oVvV.params->coltot[Gab], nbvirpi_[Gc]);

            WbAc[Gab] = global_dpd_->dpd_block_matrix(I_oVvV.params->coltot[Gab], nbvirpi_[Gc]);
          }

          for(Gd=0; Gd < nirrep_; Gd++) {
            /* +lambda_IkDc * <jD||bA> */
            Gba = Gjd = Gj ^ Gd;
            Gc = Gik ^ Gd;

            dc = L_AB.col_offset[Gik][Gd];
            jd = I_oVvV.row_offset[Gjd][J];

            I_oVvV.matrix[Gjd] = global_dpd_->dpd_block_matrix(navirpi_[Gd], I_oVvV.params->coltot[Gjd]);
            global_dpd_->buf4_mat_irrep_rd_block(&I_oVvV, Gjd, jd, navirpi_[Gd]);

            nrows = I_oVvV.params->coltot[Gjd];
            ncols = nbvirpi_[Gc];
            nlinks = navirpi_[Gd];

            if(nrows && ncols && nlinks)
              C_DGEMM('t', 'n', nrows, ncols, nlinks, 1.0,
                  &(I_oVvV.matrix[Gjd][0][0]), nrows,
                  &(L_AB.matrix[Gik][ik][dc]), ncols, 1.0,
                  &(WbAc[Gba][0][0]), ncols);

            global_dpd_->free_dpd_block(I_oVvV.matrix[Gjd], navirpi_[Gd], I_oVvV.params->coltot[Gjd]);

            /* -lambda_IjDc * <kD||bA> */
            Gba = Gkd = Gk ^ Gd;
            Gc = Gij ^ Gd;

            dc = L_AB.col_offset[Gij][Gd];
            kd = I_oVvV.row_offset[Gkd][K];

            I_oVvV.matrix[Gkd] = global_dpd_->dpd_block_matrix(navirpi_[Gd], I_oVvV.params->coltot[Gkd]);
            global_dpd_->buf4_mat_irrep_rd_block(&I_oVvV, Gkd, kd, navirpi_[Gd]);

            nrows = I_oVvV.params->coltot[Gkd];
            ncols = nbvirpi_[Gc];
            nlinks = navirpi_[Gd];

            if(nrows && ncols && nlinks)
              C_DGEMM('t', 'n', nrows, ncols, nlinks, -1.0,
                  &(I_oVvV.matrix[Gkd][0][0]), nrows,
                  &(L_AB.matrix[Gij][ij][dc]), ncols, 1.0,
                  &(WbAc[Gba][0][0]), ncols);

            global_dpd_->free_dpd_block(I_oVvV.matrix[Gkd], navirpi_[Gd], I_oVvV.params->coltot[Gkd]);

          }

          for(Gl=0; Gl < nirrep_; Gl++) {
            /* -lambda_jLbA * <Ik||Lc> */
            Gba = Gjl = Gj ^ Gl;
            Gc = Gik ^ Gl;

            lc = I_OoOv.col_offset[Gik][Gl];
            jl = L_BA.row_offset[Gjl][J];

            nrows = L_BA.params->coltot[Gjl];
            ncols = nbvirpi_[Gc];
            nlinks = naoccpi_[Gl];

            if(nrows && ncols && nlinks)
              C_DGEMM('t', 'n', nrows, ncols, nlinks, -1.0,
                  &(L_BA.matrix[Gjl][jl][0]), nrows,
                  &(I_OoOv.matrix[Gik][ik][lc]), ncols, 1.0,
                  &(WbAc[Gba][0][0]), ncols);

            /* +lambda_kLbA * <Ij||Lc> */
            Gba = Gkl = Gk ^ Gl;
            Gc = Gij ^ Gl;

            lc = I_OoOv.col_offset[Gij][Gl];
            kl = L_BA.row_offset[Gkl][K];

            nrows = L_BA.params->coltot[Gkl];
            ncols = nbvirpi_[Gc];
            nlinks = naoccpi_[Gl];

            if(nrows && ncols && nlinks)
              C_DGEMM('t', 'n', nrows, ncols, nlinks, 1.0,
                  &(L_BA.matrix[Gkl][kl][0]), nrows,
                  &(I_OoOv.matrix[Gij][ij][lc]), ncols, 1.0,
                  &(WbAc[Gba][0][0]), ncols);
          }

          global_dpd_->sort_3d(WbAc, WAbc, nirrep_, Gijk, I_oVvV.params->coltot, I_oVvV.params->colidx,
                 I_oVvV.params->colorb, I_oVvV.params->rsym, I_oVvV.params->ssym,
                 bvir_off_, avir_off_, nbvirpi_, bvir_off_, I_OvVv.params->colidx, bac, 1);

          for(Gab=0; Gab < nirrep_; Gab++) {
            Gc = Gab ^ Gijk;

            global_dpd_->free_dpd_block(WbAc[Gab], I_oVvV.params->coltot[Gab], nbvirpi_[Gc]);
          }

          /* Compute denominator and evaluate labda_ijkabc */
          for(Gab=0; Gab < nirrep_; Gab++) {
            Gc = Gab ^ Gijk;

            LAbc[Gab] = global_dpd_->dpd_block_matrix(I_OvVv.params->coltot[Gab], nbvirpi_[Gc]);

            for(ab=0; ab < I_OvVv.params->coltot[Gab]; ab++) {
              A = I_OvVv.params->colorb[Gab][ab][0];
              Ga = I_OvVv.params->rsym[A];
              a = A - avir_off_[Ga];
              B = I_OvVv.params->colorb[Gab][ab][1];
              Gb = I_OvVv.params->ssym[B];
              b = B - bvir_off_[Gb];

              Gbc = Gb ^ Gc;
              Gac = Ga ^ Gc;

              for(c=0; c < nbvirpi_[Gc]; c++) {
                C = bvir_off_[Gc] + c;

                /* Copy W into labda_ijkabc */
                LAbc[Gab][ab][c] = WAbc[Gab][ab][c];

                /* Build the rest of the denominator and compute labda_ijkabc */
                denom = dijk;
                if(Ftilde_a_->rowspi()[Ga])
                    denom -= Ftilde_a_->get(Ga, a + naoccpi_[Ga], a + naoccpi_[Ga]);
                if(Ftilde_b_->rowspi()[Gb])
                    denom -= Ftilde_b_->get(Gb, b + nboccpi_[Gb], b + nboccpi_[Gb]);
                if(Ftilde_b_->rowspi()[Gc])
                    denom -= Ftilde_b_->get(Gc, c + nboccpi_[Gc], c + nboccpi_[Gc]);


                LAbc[Gab][ab][c] /= denom;

              } /* c */
            } /* ab */
          } /* Gab */

          /* Compute the ABB energy contribution  */
          for(Gab=0; Gab < nirrep_; Gab++) {
            Gc = Gab ^ Gijk;
            ET_abb += dot_block(LAbc[Gab], WAbc[Gab], I_OvVv.params->coltot[Gab], nbvirpi_[Gc], 0.5);
            global_dpd_->free_dpd_block(WAbc[Gab], I_OvVv.params->coltot[Gab], nbvirpi_[Gc]);
            global_dpd_->free_dpd_block(LAbc[Gab], I_OvVv.params->coltot[Gab], nbvirpi_[Gc]);
          }

            } /* J >= K */

          } /* k */
        } /* j */
      } /* i */

        } /* Gk */
      } /* Gj */
    } /* Gi */

    free(WAbc);
    free(LAbc);
    free(WAcb);
    free(WbcA);
    free(WcAb);
    free(WbAc);

    for(h=0; h < nirrep_; h++) {
      global_dpd_->buf4_mat_irrep_close(&L_BB, h);
      global_dpd_->buf4_mat_irrep_close(&L_AB, h);
      global_dpd_->buf4_mat_irrep_close(&L_BA, h);
      global_dpd_->buf4_mat_irrep_close(&I_ooov, h);
      global_dpd_->buf4_mat_irrep_close(&I_OoOv, h);
      global_dpd_->buf4_mat_irrep_close(&I_oOoV, h);
    }

    global_dpd_->buf4_close(&L_BB);
    global_dpd_->buf4_close(&L_AB);
    global_dpd_->buf4_close(&L_BA);
    global_dpd_->buf4_close(&I_ovvv);
    global_dpd_->buf4_close(&I_OvVv);
    global_dpd_->buf4_close(&I_oVvV);
    global_dpd_->buf4_close(&I_ooov);
    global_dpd_->buf4_close(&I_OoOv);
    global_dpd_->buf4_close(&I_oOoV);

    psio_->close(PSIF_LIBTRANS_DPD, 1);

    return ET_abb;

}

double
DCFTSolver::compute_triples_bbb()
{

    int h;
    int Gi, Gj, Gk, Ga, Gb, Gc, Gd, Gl;
    int Gji, Gij, Gjk, Gkj, Gik, Gki, Gijk;
    int Gab, Gbc, Gac;
    int Gid, Gjd, Gkd;
    int Gil, Gjl, Gkl;
    int I, J, K, A, B, C;
    int i, j, k, a, b, c;
    int ij, ji, ik, ki, jk, kj;
    int ab;
    int cd, ad, bd;
    int id, jd, kd;
    int il, jl, kl;
    int lc, la, lb;
    double dijk, denom, ET_bbb;
    int nrows, ncols, nlinks;
    dpdbuf4 L, I_ovvv, I_ooov;
    double ***WABC, ***WBCA, ***WACB, ***LABC;

    psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

    global_dpd_->buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                           ID("[o,o]"), ID("[v,v]"), 0, "Lambda <o'o'|v'v'>");
    global_dpd_->buf4_init(&I_ovvv, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[v,v]"),
                           ID("[o,v]"), ID("[v,v]"), 0, "MO Ints <o'v'||v'v'>");
    global_dpd_->buf4_init(&I_ooov, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[o,v]"),
                           ID("[o,o]"), ID("[o,v]"), 0, "MO Ints <o'o'||o'v'>");

    for(h=0; h < nirrep_; h++) {
      global_dpd_->buf4_mat_irrep_init(&L, h);
      global_dpd_->buf4_mat_irrep_rd(&L, h);

      global_dpd_->buf4_mat_irrep_init(&I_ooov, h);
      global_dpd_->buf4_mat_irrep_rd(&I_ooov, h);
    }

    WABC = (double ***) malloc(nirrep_ * sizeof(double **));
    LABC = (double ***) malloc(nirrep_ * sizeof(double **));
    WBCA = (double ***) malloc(nirrep_ * sizeof(double **));
    WACB = (double ***) malloc(nirrep_ * sizeof(double **));

    ET_bbb = 0.0;

    for(Gi=0; Gi < nirrep_; Gi++) {
      for(Gj=0; Gj < nirrep_; Gj++) {
        for(Gk=0; Gk < nirrep_; Gk++) {

      Gij = Gji = Gi ^ Gj;
      Gjk = Gkj = Gj ^ Gk;
      Gik = Gki = Gi ^ Gk;

      Gijk = Gi ^ Gj ^ Gk;

      for(i=0; i < nboccpi_[Gi]; i++) {
        I = bocc_off_[Gi] + i;
        for(j=0; j < nboccpi_[Gj]; j++) {
          J = bocc_off_[Gj] + j;
          for(k=0; k < nboccpi_[Gk]; k++) {
            K = bocc_off_[Gk] + k;

            if(I > J && J > K) {

          ij = I_ooov.params->rowidx[I][J];
          ji = I_ooov.params->rowidx[J][I];
          jk = I_ooov.params->rowidx[J][K];
          kj = I_ooov.params->rowidx[K][J];
          ik = I_ooov.params->rowidx[I][K];
          ki = I_ooov.params->rowidx[K][I];

          dijk = 0.0;
          if(Ftilde_b_->rowspi()[Gi])
            dijk += Ftilde_b_->get(Gi, i, i);
          if(Ftilde_b_->rowspi()[Gj])
            dijk += Ftilde_b_->get(Gj, j, j);
          if(Ftilde_b_->rowspi()[Gk])
            dijk += Ftilde_b_->get(Gk, k, k);

          for(Gab=0; Gab < nirrep_; Gab++) {
            Gc = Gab ^ Gijk;

            WABC[Gab] = global_dpd_->dpd_block_matrix(I_ovvv.params->coltot[Gab], nbvirpi_[Gc]);
          }

          for(Gd=0; Gd < nirrep_; Gd++) {
            /* -lambda_jkcd * <id||ab> */
            Gab = Gid = Gi ^ Gd;
            Gc = Gjk ^ Gd;

            cd = L.col_offset[Gjk][Gc];
            id = I_ovvv.row_offset[Gid][I];

            I_ovvv.matrix[Gid] = global_dpd_->dpd_block_matrix(nbvirpi_[Gd], I_ovvv.params->coltot[Gid]);
            global_dpd_->buf4_mat_irrep_rd_block(&I_ovvv, Gid, id, nbvirpi_[Gd]);

            nrows = I_ovvv.params->coltot[Gid];
            ncols = nbvirpi_[Gc];
            nlinks = nbvirpi_[Gd];

            if(nrows && ncols && nlinks)
              C_DGEMM('t', 't', nrows, ncols, nlinks, -1.0,
                  &(I_ovvv.matrix[Gid][0][0]), nrows,
                  &(L.matrix[Gjk][jk][cd]), nlinks, 1.0,
                  &(WABC[Gab][0][0]), ncols);

            global_dpd_->free_dpd_block(I_ovvv.matrix[Gid], nbvirpi_[Gd], I_ovvv.params->coltot[Gid]);

            /* +lambda_ikcd * <jd||ab> */
            Gab = Gjd = Gj ^ Gd;
            Gc = Gik ^ Gd;

            cd = L.col_offset[Gik][Gc];
            jd = I_ovvv.row_offset[Gjd][J];

            I_ovvv.matrix[Gjd] = global_dpd_->dpd_block_matrix(nbvirpi_[Gd], I_ovvv.params->coltot[Gjd]);
            global_dpd_->buf4_mat_irrep_rd_block(&I_ovvv, Gjd, jd, nbvirpi_[Gd]);

            nrows = I_ovvv.params->coltot[Gjd];
            ncols = nbvirpi_[Gc];
            nlinks = nbvirpi_[Gd];

            if(nrows && ncols && nlinks)
              C_DGEMM('t', 't', nrows, ncols, nlinks, 1.0,
                  &(I_ovvv.matrix[Gjd][0][0]), nrows,
                  &(L.matrix[Gik][ik][cd]), nlinks, 1.0,
                  &(WABC[Gab][0][0]), ncols);

            global_dpd_->free_dpd_block(I_ovvv.matrix[Gjd], nbvirpi_[Gd], I_ovvv.params->coltot[Gjd]);

            /* +lambda_jicd * <kd||ab> */
            Gab = Gkd = Gk ^ Gd;
            Gc = Gji ^ Gd;

            cd = L.col_offset[Gji][Gc];
            kd = I_ovvv.row_offset[Gkd][K];

            I_ovvv.matrix[Gkd] = global_dpd_->dpd_block_matrix(nbvirpi_[Gd], I_ovvv.params->coltot[Gkd]);
            global_dpd_->buf4_mat_irrep_rd_block(&I_ovvv, Gkd, kd, nbvirpi_[Gd]);

            nrows = I_ovvv.params->coltot[Gkd];
            ncols = nbvirpi_[Gc];
            nlinks = nbvirpi_[Gd];

            if(nrows && ncols && nlinks)
              C_DGEMM('t', 't', nrows, ncols, nlinks, 1.0,
                  &(I_ovvv.matrix[Gkd][0][0]), nrows,
                  &(L.matrix[Gji][ji][cd]), nlinks, 1.0,
                  &(WABC[Gab][0][0]), ncols);

            global_dpd_->free_dpd_block(I_ovvv.matrix[Gkd], nbvirpi_[Gd], I_ovvv.params->coltot[Gkd]);

          }

          for(Gl=0; Gl < nirrep_; Gl++) {
            /* -lambda_ilab <jk||lc> */
            Gab = Gil = Gi ^ Gl;
            Gc = Gjk ^ Gl;

            lc = I_ooov.col_offset[Gjk][Gl];
            il = L.row_offset[Gil][I];

            nrows = L.params->coltot[Gil];
            ncols = nbvirpi_[Gc];
            nlinks = nboccpi_[Gl];

            if(nrows && ncols && nlinks)
              C_DGEMM('t', 'n', nrows, ncols, nlinks, -1.0,
                  &(L.matrix[Gil][il][0]), nrows,
                  &(I_ooov.matrix[Gjk][jk][lc]), ncols, 1.0,
                  &(WABC[Gab][0][0]), ncols);

            /* +lambda_jlab <ik||lc> */
            Gab = Gjl = Gj ^ Gl;
            Gc = Gik ^ Gl;

            lc = I_ooov.col_offset[Gik][Gl];
            jl = L.row_offset[Gjl][J];

            nrows = L.params->coltot[Gjl];
            ncols = nbvirpi_[Gc];
            nlinks = nboccpi_[Gl];

            if(nrows && ncols && nlinks)
              C_DGEMM('t', 'n', nrows, ncols, nlinks, 1.0,
                  &(L.matrix[Gjl][jl][0]), nrows,
                  &(I_ooov.matrix[Gik][ik][lc]), ncols, 1.0,
                  &(WABC[Gab][0][0]), ncols);

            /* +lambda_klab <ji||lc> */
            Gab = Gkl = Gk ^ Gl;
            Gc = Gji ^ Gl;

            lc = I_ooov.col_offset[Gji][Gl];
            kl = L.row_offset[Gkl][K];

            nrows = L.params->coltot[Gkl];
            ncols = nbvirpi_[Gc];
            nlinks = nboccpi_[Gl];

            if(nrows && ncols && nlinks)
              C_DGEMM('t', 'n', nrows, ncols, nlinks, 1.0,
                  &(L.matrix[Gkl][kl][0]), nrows,
                  &(I_ooov.matrix[Gji][ji][lc]), ncols, 1.0,
                  &(WABC[Gab][0][0]), ncols);
          }

          for(Gab=0; Gab < nirrep_; Gab++) {
            Gc = Gab ^ Gijk;

            WBCA[Gab] = global_dpd_->dpd_block_matrix(I_ovvv.params->coltot[Gab], nbvirpi_[Gc]);
          }

          for(Gd=0; Gd < nirrep_; Gd++) {
            /* -lambda_jkad * <id||bc> */
            Gbc = Gid = Gi ^ Gd;
            Ga = Gjk ^ Gd;

            ad = L.col_offset[Gjk][Ga];
            id = I_ovvv.row_offset[Gid][I];

            I_ovvv.matrix[Gid] = global_dpd_->dpd_block_matrix(nbvirpi_[Gd], I_ovvv.params->coltot[Gid]);
            global_dpd_->buf4_mat_irrep_rd_block(&I_ovvv, Gid, id, nbvirpi_[Gd]);

            nrows = I_ovvv.params->coltot[Gid];
            ncols = nbvirpi_[Ga];
            nlinks = nbvirpi_[Gd];

            if(nrows && ncols && nlinks)
              C_DGEMM('t', 't', nrows, ncols, nlinks, -1.0,
                  &(I_ovvv.matrix[Gid][0][0]), nrows,
                  &(L.matrix[Gjk][jk][ad]), nlinks, 1.0,
                  &(WBCA[Gbc][0][0]), ncols);

            global_dpd_->free_dpd_block(I_ovvv.matrix[Gid], nbvirpi_[Gd], I_ovvv.params->coltot[Gid]);

            /* +lambda_ikad * <jd||bc> */
            Gbc = Gjd = Gj ^ Gd;
            Ga = Gik ^ Gd;

            ad = L.col_offset[Gik][Ga];
            jd = I_ovvv.row_offset[Gjd][J];

            I_ovvv.matrix[Gjd] = global_dpd_->dpd_block_matrix(nbvirpi_[Gd], I_ovvv.params->coltot[Gjd]);
            global_dpd_->buf4_mat_irrep_rd_block(&I_ovvv, Gjd, jd, nbvirpi_[Gd]);

            nrows = I_ovvv.params->coltot[Gjd];
            ncols = nbvirpi_[Ga];
            nlinks = nbvirpi_[Gd];

            if(nrows && ncols && nlinks)
              C_DGEMM('t', 't', nrows, ncols, nlinks, 1.0,
                  &(I_ovvv.matrix[Gjd][0][0]), nrows,
                  &(L.matrix[Gik][ik][ad]), nlinks, 1.0,
                  &(WBCA[Gbc][0][0]), ncols);

            global_dpd_->free_dpd_block(I_ovvv.matrix[Gjd], nbvirpi_[Gd], I_ovvv.params->coltot[Gjd]);

            /* +lambda_jiad * <kd||bc> */
            Gbc = Gkd = Gk ^ Gd;
            Ga = Gji ^ Gd;

            ad = L.col_offset[Gji][Ga];
            kd = I_ovvv.row_offset[Gkd][K];

            I_ovvv.matrix[Gkd] = global_dpd_->dpd_block_matrix(nbvirpi_[Gd], I_ovvv.params->coltot[Gkd]);
            global_dpd_->buf4_mat_irrep_rd_block(&I_ovvv, Gkd, kd, nbvirpi_[Gd]);

            nrows = I_ovvv.params->coltot[Gkd];
            ncols = nbvirpi_[Ga];
            nlinks = nbvirpi_[Gd];

            if(nrows && ncols && nlinks)
              C_DGEMM('t', 't', nrows, ncols, nlinks, 1.0,
                  &(I_ovvv.matrix[Gkd][0][0]), nrows,
                  &(L.matrix[Gji][ji][ad]), nlinks, 1.0,
                  &(WBCA[Gbc][0][0]), ncols);

            global_dpd_->free_dpd_block(I_ovvv.matrix[Gkd], nbvirpi_[Gd], I_ovvv.params->coltot[Gkd]);

          }

          for(Gl=0; Gl < nirrep_; Gl++) {
            /* -lambda_ilbc * <jk||la> */
            Gbc = Gil = Gi ^ Gl;
            Ga = Gjk ^ Gl;

            la = I_ooov.col_offset[Gjk][Gl];
            il = L.row_offset[Gil][I];

            nrows = L.params->coltot[Gil];
            ncols = nbvirpi_[Ga];
            nlinks = nboccpi_[Gl];

            if(nrows && ncols && nlinks)
              C_DGEMM('t', 'n', nrows, ncols, nlinks, -1.0,
                  &(L.matrix[Gil][il][0]), nrows,
                  &(I_ooov.matrix[Gjk][jk][la]), ncols, 1.0,
                  &(WBCA[Gbc][0][0]), ncols);

            /* +lambda_jlbc <ik||la> */
            Gbc = Gjl = Gj ^ Gl;
            Ga = Gik ^ Gl;

            la = I_ooov.col_offset[Gik][Gl];
            jl = L.row_offset[Gjl][J];

            nrows = L.params->coltot[Gjl];
            ncols = nbvirpi_[Ga];
            nlinks = nboccpi_[Gl];

            if(nrows && ncols && nlinks)
              C_DGEMM('t', 'n', nrows, ncols, nlinks, 1.0,
                  &(L.matrix[Gjl][jl][0]), nrows,
                  &(I_ooov.matrix[Gik][ik][la]), ncols, 1.0,
                  &(WBCA[Gbc][0][0]), ncols);

            /* +lambda_klbc <ji||la> */
            Gbc = Gkl = Gk ^ Gl;
            Ga = Gji ^ Gl;

            la = I_ooov.col_offset[Gji][Gl];
            kl = L.row_offset[Gkl][K];

            nrows = L.params->coltot[Gkl];
            ncols = nbvirpi_[Ga];
            nlinks = nboccpi_[Gl];

            if(nrows && ncols && nlinks)
              C_DGEMM('t', 'n', nrows, ncols, nlinks, 1.0,
                  &(L.matrix[Gkl][kl][0]), nrows,
                  &(I_ooov.matrix[Gji][ji][la]), ncols, 1.0,
                  &(WBCA[Gbc][0][0]), ncols);
          }

          global_dpd_->sort_3d(WBCA, WABC, nirrep_, Gijk, I_ovvv.params->coltot, I_ovvv.params->colidx,
                 I_ovvv.params->colorb, I_ovvv.params->rsym, I_ovvv.params->ssym,
                 bvir_off_, bvir_off_, nbvirpi_, bvir_off_, I_ovvv.params->colidx, cab, 1);

          for(Gab=0; Gab < nirrep_; Gab++) {
            Gc = Gab ^ Gijk;

            global_dpd_->free_dpd_block(WBCA[Gab], I_ovvv.params->coltot[Gab], nbvirpi_[Gc]);

            WACB[Gab] = global_dpd_->dpd_block_matrix(I_ovvv.params->coltot[Gab], nbvirpi_[Gc]);
          }

          for(Gd=0; Gd < nirrep_; Gd++) {
            /* +lambda_jkbd * <id||ac> */
            Gac = Gid = Gi ^ Gd;
            Gb = Gjk ^ Gd;

            bd = L.col_offset[Gjk][Gb];
            id = I_ovvv.row_offset[Gid][I];

            I_ovvv.matrix[Gid] = global_dpd_->dpd_block_matrix(nbvirpi_[Gd], I_ovvv.params->coltot[Gid]);
            global_dpd_->buf4_mat_irrep_rd_block(&I_ovvv, Gid, id, nbvirpi_[Gd]);

            nrows = I_ovvv.params->coltot[Gid];
            ncols = nbvirpi_[Gb];
            nlinks = nbvirpi_[Gd];

            if(nrows && ncols && nlinks)
              C_DGEMM('t', 't', nrows, ncols, nlinks, 1.0,
                  &(I_ovvv.matrix[Gid][0][0]), nrows,
                  &(L.matrix[Gjk][jk][bd]), nlinks, 1.0,
                  &(WACB[Gac][0][0]), ncols);

            global_dpd_->free_dpd_block(I_ovvv.matrix[Gid], nbvirpi_[Gd], I_ovvv.params->coltot[Gid]);

            /* -lambda_ikbd * <jd||ac> */
            Gac = Gjd = Gj ^ Gd;
            Gb = Gik ^ Gd;

            bd = L.col_offset[Gik][Gb];
            jd = I_ovvv.row_offset[Gjd][J];

            I_ovvv.matrix[Gjd] = global_dpd_->dpd_block_matrix(nbvirpi_[Gd], I_ovvv.params->coltot[Gjd]);
            global_dpd_->buf4_mat_irrep_rd_block(&I_ovvv, Gjd, jd, nbvirpi_[Gd]);

            nrows = I_ovvv.params->coltot[Gjd];
            ncols = nbvirpi_[Gb];
            nlinks = nbvirpi_[Gd];

            if(nrows && ncols && nlinks)
              C_DGEMM('t', 't', nrows, ncols, nlinks, -1.0,
                  &(I_ovvv.matrix[Gjd][0][0]), nrows,
                  &(L.matrix[Gik][ik][bd]), nlinks, 1.0,
                  &(WACB[Gac][0][0]), ncols);

            global_dpd_->free_dpd_block(I_ovvv.matrix[Gjd], nbvirpi_[Gd], I_ovvv.params->coltot[Gjd]);

            /* -lambda_jibd * <kd||ac> */
            Gac = Gkd = Gk ^ Gd;
            Gb = Gji ^ Gd;

            bd = L.col_offset[Gji][Gb];
            kd = I_ovvv.row_offset[Gkd][K];

            I_ovvv.matrix[Gkd] = global_dpd_->dpd_block_matrix(nbvirpi_[Gd], I_ovvv.params->coltot[Gkd]);
            global_dpd_->buf4_mat_irrep_rd_block(&I_ovvv, Gkd, kd, nbvirpi_[Gd]);

            nrows = I_ovvv.params->coltot[Gkd];
            ncols = nbvirpi_[Gb];
            nlinks = nbvirpi_[Gd];

            if(nrows && ncols && nlinks)
              C_DGEMM('t', 't', nrows, ncols, nlinks, -1.0,
                  &(I_ovvv.matrix[Gkd][0][0]), nrows,
                  &(L.matrix[Gji][ji][bd]), nlinks, 1.0,
                  &(WACB[Gac][0][0]), ncols);

            global_dpd_->free_dpd_block(I_ovvv.matrix[Gkd], nbvirpi_[Gd], I_ovvv.params->coltot[Gkd]);

          }

          for(Gl=0; Gl < nirrep_; Gl++) {
            /* +lambda_ilac * <jk||lb> */
            Gac = Gil = Gi ^ Gl;
            Gb = Gjk ^ Gl;

            lb = I_ooov.col_offset[Gjk][Gl];
            il = L.row_offset[Gil][I];

            nrows = L.params->coltot[Gil];
            ncols = nbvirpi_[Gb];
            nlinks = nboccpi_[Gl];

            if(nrows && ncols && nlinks)
              C_DGEMM('t', 'n', nrows, ncols, nlinks, 1.0,
                  &(L.matrix[Gil][il][0]), nrows,
                  &(I_ooov.matrix[Gjk][jk][lb]), ncols, 1.0,
                  &(WACB[Gac][0][0]), ncols);

            /* -lambda_jlac * <ik||lb> */
            Gac = Gjl = Gj ^ Gl;
            Gb = Gik ^ Gl;

            lb = I_ooov.col_offset[Gik][Gl];
            jl = L.row_offset[Gjl][J];

            nrows = L.params->coltot[Gjl];
            ncols = nbvirpi_[Gb];
            nlinks = nboccpi_[Gl];

            if(nrows && ncols && nlinks)
              C_DGEMM('t', 'n', nrows, ncols, nlinks, -1.0,
                  &(L.matrix[Gjl][jl][0]), nrows,
                  &(I_ooov.matrix[Gik][ik][lb]), ncols, 1.0,
                  &(WACB[Gac][0][0]), ncols);

            /* -lambda_klac * <ji||lb> */
            Gac = Gkl = Gk ^ Gl;
            Gb = Gji ^ Gl;

            lb = I_ooov.col_offset[Gji][Gl];
            kl = L.row_offset[Gkl][K];

            nrows = L.params->coltot[Gkl];
            ncols = nbvirpi_[Gb];
            nlinks = nboccpi_[Gl];

            if(nrows && ncols && nlinks)
              C_DGEMM('t', 'n', nrows, ncols, nlinks, -1.0,
                  &(L.matrix[Gkl][kl][0]), nrows,
                  &(I_ooov.matrix[Gji][ji][lb]), ncols, 1.0,
                  &(WACB[Gac][0][0]), ncols);
          }

          global_dpd_->sort_3d(WACB, WABC, nirrep_, Gijk, I_ovvv.params->coltot, I_ovvv.params->colidx,
                 I_ovvv.params->colorb, I_ovvv.params->rsym, I_ovvv.params->ssym,
                 bvir_off_, bvir_off_, nbvirpi_, bvir_off_, I_ovvv.params->colidx, acb, 1);

          for(Gab=0; Gab < nirrep_; Gab++) {
            Gc = Gab ^ Gijk;

            global_dpd_->free_dpd_block(WACB[Gab], I_ovvv.params->coltot[Gab], nbvirpi_[Gc]);
          }

          /* Compute denominator and evaluate labda_ijkabc */
          for(Gab=0; Gab < nirrep_; Gab++) {
            Gc = Gab ^ Gijk;

            LABC[Gab] = global_dpd_->dpd_block_matrix(I_ovvv.params->coltot[Gab], nbvirpi_[Gc]);

            for(ab=0; ab < I_ovvv.params->coltot[Gab]; ab++) {
              A = I_ovvv.params->colorb[Gab][ab][0];
              Ga = I_ovvv.params->rsym[A];
              a = A - bvir_off_[Ga];
              B = I_ovvv.params->colorb[Gab][ab][1];
              Gb = I_ovvv.params->ssym[B];
              b = B - bvir_off_[Gb];

              Gbc = Gb ^ Gc;
              Gac = Ga ^ Gc;

              for(c=0; c < nbvirpi_[Gc]; c++) {
                C = bvir_off_[Gc] + c;

                /* Copy W into labda_ijkabc */
                LABC[Gab][ab][c] = WABC[Gab][ab][c];

                /* Build the rest of the denominator and compute labda_ijkabc */
                denom = dijk;

                if(Ftilde_b_->rowspi()[Ga])
                    denom -= Ftilde_b_->get(Ga, a + nboccpi_[Ga], a + nboccpi_[Ga]);
                if(Ftilde_b_->rowspi()[Gb])
                    denom -= Ftilde_b_->get(Gb, b + nboccpi_[Gb], b + nboccpi_[Gb]);
                if(Ftilde_b_->rowspi()[Gc])
                    denom -= Ftilde_b_->get(Gc, c + nboccpi_[Gc], c + nboccpi_[Gc]);

                LABC[Gab][ab][c] /= denom;

              } /* c */
            } /* ab */
          } /* Gab */

          /* Compute the BBB energy contribution  */
          for(Gab=0; Gab < nirrep_; Gab++) {
            Gc = Gab ^ Gijk;
            ET_bbb += dot_block(LABC[Gab], WABC[Gab], I_ovvv.params->coltot[Gab], nbvirpi_[Gc], 1.0/6.0);
            global_dpd_->free_dpd_block(WABC[Gab], I_ovvv.params->coltot[Gab], nbvirpi_[Gc]);
            global_dpd_->free_dpd_block(LABC[Gab], I_ovvv.params->coltot[Gab], nbvirpi_[Gc]);
          }

            } /* I >= J >= K */

          } /* k */
        } /* j */
      } /* i */

        } /* Gk */
      } /* Gj */
    } /* Gi */

    free(WABC);
    free(LABC);
    free(WBCA);
    free(WACB);

    for(h=0; h < nirrep_; h++) {
      global_dpd_->buf4_mat_irrep_close(&L, h);
      global_dpd_->buf4_mat_irrep_close(&I_ooov, h);
    }

    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&I_ovvv);
    global_dpd_->buf4_close(&I_ooov);

    psio_->close(PSIF_LIBTRANS_DPD, 1);

    return ET_bbb;

}

}} //End namespaces
