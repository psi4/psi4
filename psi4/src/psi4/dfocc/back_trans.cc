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

/** Standard library includes */
#include "psi4/libqt/qt.h"
#include "defines.h"
#include "dfocc.h"
#include "psi4/libmints/oeprop.h"
#include "psi4/libmints/matrix.h"

using namespace std;


namespace psi{ namespace dfoccwave{

void DFOCC::back_trans()
{
    outfile->Printf("\tBacktransforming OPDM, TPDM, and GFM to the AO basis...\n");


    SharedTensor2d Gmo, Gao, Gao2, G, Gref, Gsep, Gcorr, cQso2;
    timer_on("back_trans");
if (reference_ == "RESTRICTED") {

    //=========================
    // OPDM back trans
    //=========================
    G1ao->back_transform(G1, CmoA);

    //=========================
    // GFM back trans
    //=========================
    GFao->back_transform(GF, CmoA);

    //=========================
    // Seprable TPDM
    //=========================
    // OO Block Sep + OO block Ref
    Gref = SharedTensor2d(new Tensor2d("Reference 3-Index TPDM <Q|OO>", nQ_ref, noccA, noccA));
    Gmo = SharedTensor2d(new Tensor2d("3-Index Separable TPDM <Q|OO>", nQ_ref, noccA, noccA));
    Gref->read(psio_, PSIF_DFOCC_DENS);
    Gmo->read(psio_, PSIF_DFOCC_DENS);
    Gmo->add(Gref);
    Gref.reset();
    G = SharedTensor2d(new Tensor2d("3-Index Separable TPDM (Q|On)", nQ_ref, noccA, nso_));
    G->contract(false, true, nQ_ref * noccA, nso_, noccA, Gmo, CoccA, 1.0, 0.0);
    Gmo.reset();
    Gao = SharedTensor2d(new Tensor2d("RefSep 3-Index TPDM (Q|nn)", nQ_ref, nso_, nso_));
    Gao->contract233(false, false, nso_, nso_, CoccA, G, 1.0, 0.0);
    G.reset();

    // OV Block
    Gmo = SharedTensor2d(new Tensor2d("3-Index Separable TPDM <Q|OV>", nQ_ref, noccA, nvirA));
    Gmo->read(psio_, PSIF_DFOCC_DENS);
    G = SharedTensor2d(new Tensor2d("3-Index Separable TPDM (Q|On)", nQ_ref, noccA, nso_));
    G->contract(false, true, nQ_ref * noccA, nso_, nvirA, Gmo, CvirA, 1.0, 0.0);
    Gmo.reset();
    Gao->contract233(false, false, nso_, nso_, CoccA, G, 2.0, 1.0);
    G.reset();

    // VV Block
    Gmo = SharedTensor2d(new Tensor2d("3-Index Separable TPDM <Q|VV>", nQ_ref, nvirA, nvirA));
    Gmo->read(psio_, PSIF_DFOCC_DENS, true, true);
    G = SharedTensor2d(new Tensor2d("3-Index Separable TPDM (Q|Vn)", nQ_ref, nvirA, nso_));
    G->contract(false, true, nQ_ref * nvirA, nso_, nvirA, Gmo, CvirA, 1.0, 0.0);
    Gmo.reset();
    Gao->contract233(false, false, nso_, nso_, CvirA, G, 1.0, 1.0);
    G.reset();

    // symmetrize : This is necessary since we only consider OV block
    Gao->symmetrize3(Gao);
    Gao->write(psio_, PSIF_DFOCC_DENS, true, true);

    /*
    // 2-Index TPDM
    bQso = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|mn)", nQ_ref, nso_, nso_));
    bQso->read(psio_, PSIF_DFOCC_INTS, true, true);
    Jmhalf = SharedTensor2d(new Tensor2d("DF_BASIS_SCF Jmhalf <P|Q>", nQ_ref, nQ_ref));
    Jmhalf->read(psio_, PSIF_DFOCC_INTS);
    cQso = SharedTensor2d(new Tensor2d("DF_BASIS_SCF C (Q|mn)", nQ_ref, nso_, nso_));
    cQso->gemm(true, false, Jmhalf, bQso, 1.0, 0.0);
    bQso.reset();
    Jmhalf.reset();
    G = SharedTensor2d(new Tensor2d("2-Index RefSep TPDM (P|Q)", nQ_ref, nQ_ref));
    G->gemm(false, true, cQso, Gao, 0.5, 0.0);
    //G->gemm(false, true, cQso, Gao, 0.25, 0.0);
    //G->gemm(false, true, Gao, cQso, 0.25, 1.0);
    cQso.reset();
    Gao.reset();
    //G->write(psio_, PSIF_DFOCC_DENS);
    G->write_symm(psio_, PSIF_DFOCC_DENS);
    G.reset();
    */

    // 2-Index TPDM
    bQso = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|mn)", nQ_ref, nso_, nso_));
    bQso->read(psio_, PSIF_DFOCC_INTS, true, true);
    Jmhalf = SharedTensor2d(new Tensor2d("DF_BASIS_SCF Jmhalf <P|Q>", nQ_ref, nQ_ref));
    Jmhalf->read(psio_, PSIF_DFOCC_INTS);
    cQso = SharedTensor2d(new Tensor2d("DF_BASIS_SCF C (Q|mn)", nQ_ref, nso_, nso_));
    cQso->gemm(true, false, Jmhalf, bQso, 1.0, 0.0);
    bQso.reset();
    Jmhalf.reset();

    // Packed c_mn^Q
    cQso2 = SharedTensor2d(new Tensor2d("DF_BASIS_SCF C (Q|m>=n)", nQ_ref, ntri_so));
    cQso2->symm_packed(cQso);
    cQso.reset();

    // LTM of G_mn^Q
    Gao2 = SharedTensor2d(new Tensor2d("RefSep 3-Index TPDM (Q|n>=n)", nQ_ref, ntri_so));
    Gao2->ltm(Gao);
    Gao.reset();

    // G_PQ = 1/2 \sum_{mn} c_mn^P G_mn^Q = 1/2 \sum_{m>=n} c_mn^P G_mn^Q (2 - \delta_{mn})
    G = SharedTensor2d(new Tensor2d("2-Index RefSep TPDM (P|Q)", nQ_ref, nQ_ref));
    G->gemm(false, true, cQso2, Gao2, 0.5, 0.0);
    Gao2.reset();
    cQso2.reset();
    G->write_symm(psio_, PSIF_DFOCC_DENS);
    G.reset();

    //=========================
    // Correlation TPDM
    //=========================
    // OV Block
    Gmo = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM <Q|OV>", nQ, noccA, nvirA));
    Gmo->read(psio_, PSIF_DFOCC_DENS);
    G = SharedTensor2d(new Tensor2d("3-Index Correlation TPDM (Q|On)", nQ, noccA, nso_));
    G->contract(false, true, nQ * noccA, nso_, nvirA, Gmo, CvirA, 1.0, 0.0);
    Gmo.reset();
    Gao = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|nn)", nQ, nso_, nso_));
    Gao->contract233(false, false, nso_, nso_, CoccA, G, 2.0, 0.0);
    G.reset();

    // symmetrize : This is necessary since we only consider OV block
    Gao->symmetrize3(Gao);
    Gao->write(psio_, PSIF_DFOCC_DENS, true, true);

    /*
    // 2-Index TPDM
    bQso = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|mn)", nQ, nso_, nso_));
    bQso->read(psio_, PSIF_DFOCC_INTS, true, true);
    Jmhalf = SharedTensor2d(new Tensor2d("DF_BASIS_CC Jmhalf <P|Q>", nQ, nQ));
    Jmhalf->read(psio_, PSIF_DFOCC_INTS);
    cQso = SharedTensor2d(new Tensor2d("DF_BASIS_CC C (Q|mn)", nQ, nso_, nso_));
    cQso->gemm(true, false, Jmhalf, bQso, 1.0, 0.0);
    bQso.reset();
    Jmhalf.reset();
    G = SharedTensor2d(new Tensor2d("2-Index Correlation TPDM (P|Q)", nQ, nQ));
    G->gemm(false, true, cQso, Gao, 0.5, 0.0);
    cQso.reset();
    Gao.reset();
    //G->write(psio_, PSIF_DFOCC_DENS);
    G->write_symm(psio_, PSIF_DFOCC_DENS);
    G.reset();
    */

    // 2-Index TPDM
    bQso = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|mn)", nQ, nso_, nso_));
    bQso->read(psio_, PSIF_DFOCC_INTS, true, true);
    Jmhalf = SharedTensor2d(new Tensor2d("DF_BASIS_CC Jmhalf <P|Q>", nQ, nQ));
    Jmhalf->read(psio_, PSIF_DFOCC_INTS);
    cQso = SharedTensor2d(new Tensor2d("DF_BASIS_CC C (Q|mn)", nQ, nso_, nso_));
    cQso->gemm(true, false, Jmhalf, bQso, 1.0, 0.0);
    bQso.reset();
    Jmhalf.reset();

    // Packed c_mn^Q
    cQso2 = SharedTensor2d(new Tensor2d("DF_BASIS_CC C (Q|m>=n)", nQ, ntri_so));
    cQso2->symm_packed(cQso);
    cQso.reset();

    // LTM of G_mn^Q
    Gao2 = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|n>=n)", nQ, ntri_so));
    Gao2->ltm(Gao);
    Gao.reset();

    // G_PQ = 1/2 \sum_{mn} c_mn^P G_mn^Q = 1/2 \sum_{m>=n} c_mn^P G_mn^Q (2 - \delta_{mn})
    G = SharedTensor2d(new Tensor2d("2-Index Correlation TPDM (P|Q)", nQ, nQ));
    G->gemm(false, true, cQso2, Gao2, 0.5, 0.0);
    Gao2.reset();
    cQso2.reset();
    G->write_symm(psio_, PSIF_DFOCC_DENS);
    G.reset();

}// end if (reference_ == "RESTRICTED")

else if (reference_ == "UNRESTRICTED") {

    //=========================
    // OPDM back trans
    //=========================
    G1ao->back_transform(G1A, CmoA);
    G1ao->back_transform(G1B, CmoB, 1.0, 1.0);

    //=========================
    // GFM back trans
    //=========================
    GFao->back_transform(GFA, CmoA);
    GFao->back_transform(GFB, CmoB, 1.0, 1.0);

    //=========================
    // Seprable TPDM
    //=========================
    // OO Block sep + OO Block ref
    Gref = SharedTensor2d(new Tensor2d("Reference 3-Index TPDM <Q|OO>", nQ_ref, noccA, noccA));
    Gmo = SharedTensor2d(new Tensor2d("3-Index Separable TPDM <Q|OO>", nQ_ref, noccA, noccA));
    Gref->read(psio_, PSIF_DFOCC_DENS);
    Gmo->read(psio_, PSIF_DFOCC_DENS);
    Gmo->add(Gref);
    Gref.reset();
    G = SharedTensor2d(new Tensor2d("3-Index Separable TPDM (Q|On)", nQ_ref, noccA, nso_));
    G->contract(false, true, nQ_ref * noccA, nso_, noccA, Gmo, CoccA, 1.0, 0.0);
    Gmo.reset();
    Gao = SharedTensor2d(new Tensor2d("RefSep 3-Index TPDM (Q|nn)", nQ_ref, nso_, nso_));
    Gao->contract233(false, false, nso_, nso_, CoccA, G, 1.0, 0.0);
    G.reset();

    // oo Block sep + oo block ref
    Gref = SharedTensor2d(new Tensor2d("Reference 3-Index TPDM <Q|oo>", nQ_ref, noccB, noccB));
    Gmo = SharedTensor2d(new Tensor2d("3-Index Separable TPDM <Q|oo>", nQ_ref, noccB, noccB));
    Gref->read(psio_, PSIF_DFOCC_DENS);
    Gmo->read(psio_, PSIF_DFOCC_DENS);
    Gmo->add(Gref);
    Gref.reset();
    G = SharedTensor2d(new Tensor2d("3-Index Separable TPDM (Q|on)", nQ_ref, noccB, nso_));
    G->contract(false, true, nQ_ref * noccB, nso_, noccB, Gmo, CoccB, 1.0, 0.0);
    Gmo.reset();
    Gao->contract233(false, false, nso_, nso_, CoccB, G, 1.0, 1.0);
    G.reset();

    // OV Block
    Gmo = SharedTensor2d(new Tensor2d("3-Index Separable TPDM <Q|OV>", nQ_ref, noccA, nvirA));
    Gmo->read(psio_, PSIF_DFOCC_DENS);
    G = SharedTensor2d(new Tensor2d("3-Index Separable TPDM (Q|On)", nQ_ref, noccA, nso_));
    G->contract(false, true, nQ_ref * noccA, nso_, nvirA, Gmo, CvirA, 1.0, 0.0);
    Gmo.reset();
    Gao->contract233(false, false, nso_, nso_, CoccA, G, 2.0, 1.0);
    G.reset();

    // ov Block
    Gmo = SharedTensor2d(new Tensor2d("3-Index Separable TPDM <Q|ov>", nQ_ref, noccB, nvirB));
    Gmo->read(psio_, PSIF_DFOCC_DENS);
    G = SharedTensor2d(new Tensor2d("3-Index Separable TPDM (Q|on)", nQ_ref, noccB, nso_));
    G->contract(false, true, nQ_ref * noccB, nso_, nvirB, Gmo, CvirB, 1.0, 0.0);
    Gmo.reset();
    Gao->contract233(false, false, nso_, nso_, CoccB, G, 2.0, 1.0);
    G.reset();

    // VV Block
    Gmo = SharedTensor2d(new Tensor2d("3-Index Separable TPDM <Q|VV>", nQ_ref, nvirA, nvirA));
    Gmo->read(psio_, PSIF_DFOCC_DENS, true, true);
    G = SharedTensor2d(new Tensor2d("3-Index Separable TPDM (Q|Vn)", nQ_ref, nvirA, nso_));
    G->contract(false, true, nQ_ref * nvirA, nso_, nvirA, Gmo, CvirA, 1.0, 0.0);
    Gmo.reset();
    Gao->contract233(false, false, nso_, nso_, CvirA, G, 1.0, 1.0);
    G.reset();

    // vv Block
    Gmo = SharedTensor2d(new Tensor2d("3-Index Separable TPDM <Q|vv>", nQ_ref, nvirB, nvirB));
    Gmo->read(psio_, PSIF_DFOCC_DENS, true, true);
    G = SharedTensor2d(new Tensor2d("3-Index Separable TPDM (Q|vn)", nQ_ref, nvirB, nso_));
    G->contract(false, true, nQ_ref * nvirB, nso_, nvirB, Gmo, CvirB, 1.0, 0.0);
    Gmo.reset();
    Gao->contract233(false, false, nso_, nso_, CvirB, G, 1.0, 1.0);
    G.reset();

    // symmetrize : This is necessary since we only consider OV block
    Gao->symmetrize3(Gao);
    Gao->write(psio_, PSIF_DFOCC_DENS, true, true);

    /*
    // 2-Index TPDM
    bQso = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|mn)", nQ_ref, nso_, nso_));
    bQso->read(psio_, PSIF_DFOCC_INTS, true, true);
    Jmhalf = SharedTensor2d(new Tensor2d("DF_BASIS_SCF Jmhalf <P|Q>", nQ_ref, nQ_ref));
    Jmhalf->read(psio_, PSIF_DFOCC_INTS);
    cQso = SharedTensor2d(new Tensor2d("DF_BASIS_SCF C (Q|mn)", nQ_ref, nso_, nso_));
    cQso->gemm(true, false, Jmhalf, bQso, 1.0, 0.0);
    bQso.reset();
    Jmhalf.reset();
    G = SharedTensor2d(new Tensor2d("2-Index RefSep TPDM (P|Q)", nQ_ref, nQ_ref));
    G->gemm(false, true, cQso, Gao, 0.5, 0.0);
    cQso.reset();
    Gao.reset();
    //G->write(psio_, PSIF_DFOCC_DENS);
    G->write_symm(psio_, PSIF_DFOCC_DENS);
    G.reset();
    */

    // 2-Index TPDM
    bQso = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|mn)", nQ_ref, nso_, nso_));
    bQso->read(psio_, PSIF_DFOCC_INTS, true, true);
    Jmhalf = SharedTensor2d(new Tensor2d("DF_BASIS_SCF Jmhalf <P|Q>", nQ_ref, nQ_ref));
    Jmhalf->read(psio_, PSIF_DFOCC_INTS);
    cQso = SharedTensor2d(new Tensor2d("DF_BASIS_SCF C (Q|mn)", nQ_ref, nso_, nso_));
    cQso->gemm(true, false, Jmhalf, bQso, 1.0, 0.0);
    bQso.reset();
    Jmhalf.reset();

    // Packed c_mn^Q
    cQso2 = SharedTensor2d(new Tensor2d("DF_BASIS_SCF C (Q|m>=n)", nQ_ref, ntri_so));
    cQso2->symm_packed(cQso);
    cQso.reset();

    // LTM of G_mn^Q
    Gao2 = SharedTensor2d(new Tensor2d("RefSep 3-Index TPDM (Q|n>=n)", nQ_ref, ntri_so));
    Gao2->ltm(Gao);
    Gao.reset();

    // G_PQ = 1/2 \sum_{mn} c_mn^P G_mn^Q = 1/2 \sum_{m>=n} c_mn^P G_mn^Q (2 - \delta_{mn})
    G = SharedTensor2d(new Tensor2d("2-Index RefSep TPDM (P|Q)", nQ_ref, nQ_ref));
    G->gemm(false, true, cQso2, Gao2, 0.5, 0.0);
    Gao2.reset();
    cQso2.reset();
    G->write_symm(psio_, PSIF_DFOCC_DENS);
    G.reset();

    //=========================
    // Correlation TPDM
    //=========================
    // OV Block
    Gmo = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM <Q|OV>", nQ, noccA, nvirA));
    Gmo->read(psio_, PSIF_DFOCC_DENS);
    G = SharedTensor2d(new Tensor2d("3-Index Correlation TPDM (Q|On)", nQ, noccA, nso_));
    G->contract(false, true, nQ * noccA, nso_, nvirA, Gmo, CvirA, 1.0, 0.0);
    Gmo.reset();
    Gao = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|nn)", nQ, nso_, nso_));
    Gao->contract233(false, false, nso_, nso_, CoccA, G, 2.0, 0.0);
    G.reset();

    // ov Block
    Gmo = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM <Q|ov>", nQ, noccB, nvirB));
    Gmo->read(psio_, PSIF_DFOCC_DENS);
    G = SharedTensor2d(new Tensor2d("3-Index Correlation TPDM (Q|on)", nQ, noccB, nso_));
    G->contract(false, true, nQ * noccB, nso_, nvirB, Gmo, CvirB, 1.0, 0.0);
    Gmo.reset();
    Gao->contract233(false, false, nso_, nso_, CoccB, G, 2.0, 1.0);
    G.reset();

    // symmetrize : This is necessary since we only consider OV block
    Gao->symmetrize3(Gao);
    Gao->write(psio_, PSIF_DFOCC_DENS, true, true);

    /*
    // 2-Index TPDM
    bQso = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|mn)", nQ, nso_, nso_));
    bQso->read(psio_, PSIF_DFOCC_INTS, true, true);
    Jmhalf = SharedTensor2d(new Tensor2d("DF_BASIS_CC Jmhalf <P|Q>", nQ, nQ));
    Jmhalf->read(psio_, PSIF_DFOCC_INTS);
    cQso = SharedTensor2d(new Tensor2d("DF_BASIS_CC C (Q|mn)", nQ, nso_, nso_));
    cQso->gemm(true, false, Jmhalf, bQso, 1.0, 0.0);
    bQso.reset();
    Jmhalf.reset();
    G = SharedTensor2d(new Tensor2d("2-Index Correlation TPDM (P|Q)", nQ, nQ));
    G->gemm(false, true, cQso, Gao, 0.5, 0.0);
    cQso.reset();
    Gao.reset();
    //G->write(psio_, PSIF_DFOCC_DENS);
    G->write_symm(psio_, PSIF_DFOCC_DENS);
    G.reset();
    */

    // 2-Index TPDM
    bQso = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|mn)", nQ, nso_, nso_));
    bQso->read(psio_, PSIF_DFOCC_INTS, true, true);
    Jmhalf = SharedTensor2d(new Tensor2d("DF_BASIS_CC Jmhalf <P|Q>", nQ, nQ));
    Jmhalf->read(psio_, PSIF_DFOCC_INTS);
    cQso = SharedTensor2d(new Tensor2d("DF_BASIS_CC C (Q|mn)", nQ, nso_, nso_));
    cQso->gemm(true, false, Jmhalf, bQso, 1.0, 0.0);
    bQso.reset();
    Jmhalf.reset();

    // Packed c_mn^Q
    cQso2 = SharedTensor2d(new Tensor2d("DF_BASIS_CC C (Q|m>=n)", nQ, ntri_so));
    cQso2->symm_packed(cQso);
    cQso.reset();

    // LTM of G_mn^Q
    Gao2 = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|n>=n)", nQ, ntri_so));
    Gao2->ltm(Gao);
    Gao.reset();

    // G_PQ = 1/2 \sum_{mn} c_mn^P G_mn^Q = 1/2 \sum_{m>=n} c_mn^P G_mn^Q (2 - \delta_{mn})
    G = SharedTensor2d(new Tensor2d("2-Index Correlation TPDM (P|Q)", nQ, nQ));
    G->gemm(false, true, cQso2, Gao2, 0.5, 0.0);
    Gao2.reset();
    cQso2.reset();
    G->write_symm(psio_, PSIF_DFOCC_DENS);
    G.reset();

}// else if (reference_ == "UNRESTRICTED")
    timer_off("back_trans");
    //outfile->Printf("\tBacktransformation is done.\n");
    //
} // end back_trans

//=========================
// OEPROP
//=========================
void DFOCC::oeprop()
{
    outfile->Printf("\tComputing one-electron properties...\n");

    timer_on("oeprop");
    SharedMatrix Da_ = SharedMatrix(new Matrix("MO-basis alpha OPDM", nmo_, nmo_));
    SharedMatrix Db_ = SharedMatrix(new Matrix("MO-basis beta OPDM", nmo_, nmo_));
    if (reference_ == "RESTRICTED") {
        G1->to_shared_matrix(Da_);
        Da_->scale(0.5);
        Db_->copy(Da_);
    }

    else if (reference_ == "UNRESTRICTED") {
        G1A->to_shared_matrix(Da_);
        G1B->to_shared_matrix(Db_);
    }

    // Compute oeprop
    std::shared_ptr<OEProp> oe(new OEProp(shared_from_this()));
    oe->set_Da_mo(Da_);
    if (reference_ == "UNRESTRICTED") oe->set_Db_mo(Db_);
    oe->add("DIPOLE");
    oe->add("QUADRUPOLE");
    oe->add("MULLIKEN_CHARGES");
    oe->add("NO_OCCUPATIONS");
    oe->set_title(wfn_type_);
    oe->compute();
    Da_.reset();
    Db_.reset();

    timer_off("oeprop");
} // end oeprop

//======================================================================
//    CCSD: Back Trans
//======================================================================
void DFOCC::back_trans_cc()
{
    outfile->Printf("\tBacktransforming OPDM, TPDM, and GFM to the AO basis...\n");


    SharedTensor2d Gmo, Gao, Gao2, G, Gref, Gsep, Gcorr, cQso2;
    timer_on("back_trans");
if (reference_ == "RESTRICTED") {

    //=========================
    // OPDM back trans
    //=========================
    G1ao->back_transform(G1, CmoA);

    //=========================
    // GFM back trans
    //=========================
    GFao->back_transform(GF, CmoA);

    //=========================
    // Seprable TPDM
    //=========================
    // OO Block Sep + OO block Ref
    Gref = SharedTensor2d(new Tensor2d("Reference 3-Index TPDM <Q|OO>", nQ_ref, noccA, noccA));
    Gmo = SharedTensor2d(new Tensor2d("3-Index Separable TPDM <Q|OO>", nQ_ref, noccA, noccA));
    Gref->read(psio_, PSIF_DFOCC_DENS);
    Gmo->read(psio_, PSIF_DFOCC_DENS);
    Gmo->add(Gref);
    Gref.reset();
    G = SharedTensor2d(new Tensor2d("3-Index Separable TPDM (Q|On)", nQ_ref, noccA, nso_));
    G->contract(false, true, nQ_ref * noccA, nso_, noccA, Gmo, CoccA, 1.0, 0.0);
    Gmo.reset();
    Gao = SharedTensor2d(new Tensor2d("RefSep 3-Index TPDM (Q|nn)", nQ_ref, nso_, nso_));
    Gao->contract233(false, false, nso_, nso_, CoccA, G, 1.0, 0.0);
    G.reset();

    // OV Block
    Gmo = SharedTensor2d(new Tensor2d("3-Index Separable TPDM <Q|OV>", nQ_ref, noccA, nvirA));
    Gmo->read(psio_, PSIF_DFOCC_DENS);
    G = SharedTensor2d(new Tensor2d("3-Index Separable TPDM (Q|On)", nQ_ref, noccA, nso_));
    G->contract(false, true, nQ_ref * noccA, nso_, nvirA, Gmo, CvirA, 1.0, 0.0);
    Gmo.reset();
    Gao->contract233(false, false, nso_, nso_, CoccA, G, 2.0, 1.0);
    G.reset();

    // VV Block
    Gmo = SharedTensor2d(new Tensor2d("3-Index Separable TPDM <Q|VV>", nQ_ref, nvirA, nvirA));
    Gmo->read(psio_, PSIF_DFOCC_DENS, true, true);
    G = SharedTensor2d(new Tensor2d("3-Index Separable TPDM (Q|Vn)", nQ_ref, nvirA, nso_));
    G->contract(false, true, nQ_ref * nvirA, nso_, nvirA, Gmo, CvirA, 1.0, 0.0);
    Gmo.reset();
    Gao->contract233(false, false, nso_, nso_, CvirA, G, 1.0, 1.0);
    G.reset();

    // symmetrize : This is necessary since we only consider OV block
    Gao->symmetrize3(Gao);
    Gao->write(psio_, PSIF_DFOCC_DENS, true, true);

    // 2-Index TPDM
    bQso = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|mn)", nQ_ref, nso_, nso_));
    bQso->read(psio_, PSIF_DFOCC_INTS, true, true);
    Jmhalf = SharedTensor2d(new Tensor2d("DF_BASIS_SCF Jmhalf <P|Q>", nQ_ref, nQ_ref));
    Jmhalf->read(psio_, PSIF_DFOCC_INTS);
    cQso = SharedTensor2d(new Tensor2d("DF_BASIS_SCF C (Q|mn)", nQ_ref, nso_, nso_));
    cQso->gemm(true, false, Jmhalf, bQso, 1.0, 0.0);
    bQso.reset();
    Jmhalf.reset();

    // Packed c_mn^Q
    cQso2 = SharedTensor2d(new Tensor2d("DF_BASIS_SCF C (Q|m>=n)", nQ_ref, ntri_so));
    cQso2->symm_packed(cQso);
    cQso.reset();

    // LTM of G_mn^Q
    Gao2 = SharedTensor2d(new Tensor2d("RefSep 3-Index TPDM (Q|n>=n)", nQ_ref, ntri_so));
    Gao2->ltm(Gao);
    Gao.reset();

    // G_PQ = 1/2 \sum_{mn} c_mn^P G_mn^Q = 1/2 \sum_{m>=n} c_mn^P G_mn^Q (2 - \delta_{mn})
    G = SharedTensor2d(new Tensor2d("2-Index RefSep TPDM (P|Q)", nQ_ref, nQ_ref));
    G->gemm(false, true, cQso2, Gao2, 0.5, 0.0);
    Gao2.reset();
    cQso2.reset();
    G->write_symm(psio_, PSIF_DFOCC_DENS);
    G.reset();

    //=========================
    // Correlation TPDM
    //=========================
    // OO Block
    Gmo = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM <Q|OO>", nQ, noccA, noccA));
    Gmo->read(psio_, PSIF_DFOCC_DENS);
    G = SharedTensor2d(new Tensor2d("3-Index Correlation TPDM (Q|On)", nQ, noccA, nso_));
    G->contract(false, true, nQ * noccA, nso_, noccA, Gmo, CoccA, 1.0, 0.0);
    Gmo.reset();
    Gao = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|nn)", nQ, nso_, nso_));
    Gao->contract233(false, false, nso_, nso_, CoccA, G, 1.0, 0.0);
    G.reset();

    // OV Block
    Gmo = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM <Q|OV>", nQ, noccA, nvirA));
    Gmo->read(psio_, PSIF_DFOCC_DENS);
    G = SharedTensor2d(new Tensor2d("3-Index Correlation TPDM (Q|On)", nQ, noccA, nso_));
    G->contract(false, true, nQ * noccA, nso_, nvirA, Gmo, CvirA, 1.0, 0.0);
    Gmo.reset();
    Gao->contract233(false, false, nso_, nso_, CoccA, G, 2.0, 1.0);
    G.reset();

    // VV Block
    Gmo = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM <Q|VV>", nQ, nvirA, nvirA));
    Gmo->read(psio_, PSIF_DFOCC_DENS, true, true);
    G = SharedTensor2d(new Tensor2d("3-Index Correlation TPDM (Q|Vn)", nQ, nvirA, nso_));
    G->contract(false, true, nQ * nvirA, nso_, nvirA, Gmo, CvirA, 1.0, 0.0);
    Gmo.reset();
    Gao->contract233(false, false, nso_, nso_, CvirA, G, 1.0, 1.0);
    G.reset();

    // symmetrize : This is necessary since we only consider OV block
    Gao->symmetrize3(Gao);
    Gao->write(psio_, PSIF_DFOCC_DENS, true, true);

    // 2-Index TPDM
    bQso = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|mn)", nQ, nso_, nso_));
    bQso->read(psio_, PSIF_DFOCC_INTS, true, true);
    Jmhalf = SharedTensor2d(new Tensor2d("DF_BASIS_CC Jmhalf <P|Q>", nQ, nQ));
    Jmhalf->read(psio_, PSIF_DFOCC_INTS);
    cQso = SharedTensor2d(new Tensor2d("DF_BASIS_CC C (Q|mn)", nQ, nso_, nso_));
    cQso->gemm(true, false, Jmhalf, bQso, 1.0, 0.0);
    bQso.reset();
    Jmhalf.reset();

    // Packed c_mn^Q
    cQso2 = SharedTensor2d(new Tensor2d("DF_BASIS_CC C (Q|m>=n)", nQ, ntri_so));
    cQso2->symm_packed(cQso);
    cQso.reset();

    // LTM of G_mn^Q
    Gao2 = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|n>=n)", nQ, ntri_so));
    Gao2->ltm(Gao);
    Gao.reset();

    // G_PQ = 1/2 \sum_{mn} c_mn^P G_mn^Q = 1/2 \sum_{m>=n} c_mn^P G_mn^Q (2 - \delta_{mn})
    G = SharedTensor2d(new Tensor2d("2-Index Correlation TPDM (P|Q)", nQ, nQ));
    G->gemm(false, true, cQso2, Gao2, 0.5, 0.0);
    Gao2.reset();
    cQso2.reset();
    G->write_symm(psio_, PSIF_DFOCC_DENS);
    G.reset();

}// end if (reference_ == "RESTRICTED")

else if (reference_ == "UNRESTRICTED") {

    //=========================
    // OPDM back trans
    //=========================
    G1ao->back_transform(G1A, CmoA);
    G1ao->back_transform(G1B, CmoB, 1.0, 1.0);

    //=========================
    // GFM back trans
    //=========================
    GFao->back_transform(GFA, CmoA);
    GFao->back_transform(GFB, CmoB, 1.0, 1.0);

    //=========================
    // Seprable TPDM
    //=========================
    // OO Block sep + OO Block ref
    Gref = SharedTensor2d(new Tensor2d("Reference 3-Index TPDM <Q|OO>", nQ_ref, noccA, noccA));
    Gmo = SharedTensor2d(new Tensor2d("3-Index Separable TPDM <Q|OO>", nQ_ref, noccA, noccA));
    Gref->read(psio_, PSIF_DFOCC_DENS);
    Gmo->read(psio_, PSIF_DFOCC_DENS);
    Gmo->add(Gref);
    Gref.reset();
    G = SharedTensor2d(new Tensor2d("3-Index Separable TPDM (Q|On)", nQ_ref, noccA, nso_));
    G->contract(false, true, nQ_ref * noccA, nso_, noccA, Gmo, CoccA, 1.0, 0.0);
    Gmo.reset();
    Gao = SharedTensor2d(new Tensor2d("RefSep 3-Index TPDM (Q|nn)", nQ_ref, nso_, nso_));
    Gao->contract233(false, false, nso_, nso_, CoccA, G, 1.0, 0.0);
    G.reset();

    // oo Block sep + oo block ref
    Gref = SharedTensor2d(new Tensor2d("Reference 3-Index TPDM <Q|oo>", nQ_ref, noccB, noccB));
    Gmo = SharedTensor2d(new Tensor2d("3-Index Separable TPDM <Q|oo>", nQ_ref, noccB, noccB));
    Gref->read(psio_, PSIF_DFOCC_DENS);
    Gmo->read(psio_, PSIF_DFOCC_DENS);
    Gmo->add(Gref);
    Gref.reset();
    G = SharedTensor2d(new Tensor2d("3-Index Separable TPDM (Q|on)", nQ_ref, noccB, nso_));
    G->contract(false, true, nQ_ref * noccB, nso_, noccB, Gmo, CoccB, 1.0, 0.0);
    Gmo.reset();
    Gao->contract233(false, false, nso_, nso_, CoccB, G, 1.0, 1.0);
    G.reset();

    // OV Block
    Gmo = SharedTensor2d(new Tensor2d("3-Index Separable TPDM <Q|OV>", nQ_ref, noccA, nvirA));
    Gmo->read(psio_, PSIF_DFOCC_DENS);
    G = SharedTensor2d(new Tensor2d("3-Index Separable TPDM (Q|On)", nQ_ref, noccA, nso_));
    G->contract(false, true, nQ_ref * noccA, nso_, nvirA, Gmo, CvirA, 1.0, 0.0);
    Gmo.reset();
    Gao->contract233(false, false, nso_, nso_, CoccA, G, 2.0, 1.0);
    G.reset();

    // ov Block
    Gmo = SharedTensor2d(new Tensor2d("3-Index Separable TPDM <Q|ov>", nQ_ref, noccB, nvirB));
    Gmo->read(psio_, PSIF_DFOCC_DENS);
    G = SharedTensor2d(new Tensor2d("3-Index Separable TPDM (Q|on)", nQ_ref, noccB, nso_));
    G->contract(false, true, nQ_ref * noccB, nso_, nvirB, Gmo, CvirB, 1.0, 0.0);
    Gmo.reset();
    Gao->contract233(false, false, nso_, nso_, CoccB, G, 2.0, 1.0);
    G.reset();

    // VV Block
    Gmo = SharedTensor2d(new Tensor2d("3-Index Separable TPDM <Q|VV>", nQ_ref, nvirA, nvirA));
    Gmo->read(psio_, PSIF_DFOCC_DENS, true, true);
    G = SharedTensor2d(new Tensor2d("3-Index Separable TPDM (Q|Vn)", nQ_ref, nvirA, nso_));
    G->contract(false, true, nQ_ref * nvirA, nso_, nvirA, Gmo, CvirA, 1.0, 0.0);
    Gmo.reset();
    Gao->contract233(false, false, nso_, nso_, CvirA, G, 1.0, 1.0);
    G.reset();

    // vv Block
    Gmo = SharedTensor2d(new Tensor2d("3-Index Separable TPDM <Q|vv>", nQ_ref, nvirB, nvirB));
    Gmo->read(psio_, PSIF_DFOCC_DENS, true, true);
    G = SharedTensor2d(new Tensor2d("3-Index Separable TPDM (Q|vn)", nQ_ref, nvirB, nso_));
    G->contract(false, true, nQ_ref * nvirB, nso_, nvirB, Gmo, CvirB, 1.0, 0.0);
    Gmo.reset();
    Gao->contract233(false, false, nso_, nso_, CvirB, G, 1.0, 1.0);
    G.reset();

    // symmetrize : This is necessary since we only consider OV block
    Gao->symmetrize3(Gao);
    Gao->write(psio_, PSIF_DFOCC_DENS, true, true);

    // 2-Index TPDM
    bQso = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|mn)", nQ_ref, nso_, nso_));
    bQso->read(psio_, PSIF_DFOCC_INTS, true, true);
    Jmhalf = SharedTensor2d(new Tensor2d("DF_BASIS_SCF Jmhalf <P|Q>", nQ_ref, nQ_ref));
    Jmhalf->read(psio_, PSIF_DFOCC_INTS);
    cQso = SharedTensor2d(new Tensor2d("DF_BASIS_SCF C (Q|mn)", nQ_ref, nso_, nso_));
    cQso->gemm(true, false, Jmhalf, bQso, 1.0, 0.0);
    bQso.reset();
    Jmhalf.reset();

    // Packed c_mn^Q
    cQso2 = SharedTensor2d(new Tensor2d("DF_BASIS_SCF C (Q|m>=n)", nQ_ref, ntri_so));
    cQso2->symm_packed(cQso);
    cQso.reset();

    // LTM of G_mn^Q
    Gao2 = SharedTensor2d(new Tensor2d("RefSep 3-Index TPDM (Q|n>=n)", nQ_ref, ntri_so));
    Gao2->ltm(Gao);
    Gao.reset();

    // G_PQ = 1/2 \sum_{mn} c_mn^P G_mn^Q = 1/2 \sum_{m>=n} c_mn^P G_mn^Q (2 - \delta_{mn})
    G = SharedTensor2d(new Tensor2d("2-Index RefSep TPDM (P|Q)", nQ_ref, nQ_ref));
    G->gemm(false, true, cQso2, Gao2, 0.5, 0.0);
    Gao2.reset();
    cQso2.reset();
    G->write_symm(psio_, PSIF_DFOCC_DENS);
    G.reset();

    //=========================
    // Correlation TPDM
    //=========================
    // OO Block
    Gmo = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM <Q|OO>", nQ, noccA, noccA));
    Gmo->read(psio_, PSIF_DFOCC_DENS);
    G = SharedTensor2d(new Tensor2d("3-Index Correlation TPDM (Q|On)", nQ, noccA, nso_));
    G->contract(false, true, nQ * noccA, nso_, noccA, Gmo, CoccA, 1.0, 0.0);
    Gmo.reset();
    Gao = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|nn)", nQ, nso_, nso_));
    Gao->contract233(false, false, nso_, nso_, CoccA, G, 1.0, 0.0);
    G.reset();

    // oo Block
    Gmo = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM <Q|oo>", nQ, noccB, noccB));
    Gmo->read(psio_, PSIF_DFOCC_DENS);
    G = SharedTensor2d(new Tensor2d("3-Index Correlation TPDM (Q|on)", nQ, noccB, nso_));
    G->contract(false, true, nQ * noccB, nso_, noccB, Gmo, CoccB, 1.0, 0.0);
    Gmo.reset();
    Gao->contract233(false, false, nso_, nso_, CoccB, G, 1.0, 1.0);
    G.reset();

    // OV Block
    Gmo = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM <Q|OV>", nQ, noccA, nvirA));
    Gmo->read(psio_, PSIF_DFOCC_DENS);
    G = SharedTensor2d(new Tensor2d("3-Index Correlation TPDM (Q|On)", nQ, noccA, nso_));
    G->contract(false, true, nQ * noccA, nso_, nvirA, Gmo, CvirA, 1.0, 0.0);
    Gmo.reset();
    Gao->contract233(false, false, nso_, nso_, CoccA, G, 2.0, 1.0);
    G.reset();

    // ov Block
    Gmo = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM <Q|ov>", nQ, noccB, nvirB));
    Gmo->read(psio_, PSIF_DFOCC_DENS);
    G = SharedTensor2d(new Tensor2d("3-Index Correlation TPDM (Q|on)", nQ, noccB, nso_));
    G->contract(false, true, nQ * noccB, nso_, nvirB, Gmo, CvirB, 1.0, 0.0);
    Gmo.reset();
    Gao->contract233(false, false, nso_, nso_, CoccB, G, 2.0, 1.0);
    G.reset();

    // VV Block
    Gmo = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM <Q|VV>", nQ, nvirA, nvirA));
    Gmo->read(psio_, PSIF_DFOCC_DENS, true, true);
    G = SharedTensor2d(new Tensor2d("3-Index Correlation TPDM (Q|Vn)", nQ, nvirA, nso_));
    G->contract(false, true, nQ * nvirA, nso_, nvirA, Gmo, CvirA, 1.0, 0.0);
    Gmo.reset();
    Gao->contract233(false, false, nso_, nso_, CvirA, G, 1.0, 1.0);
    G.reset();

    // vv Block
    Gmo = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM <Q|vv>", nQ, nvirB, nvirB));
    Gmo->read(psio_, PSIF_DFOCC_DENS, true, true);
    G = SharedTensor2d(new Tensor2d("3-Index Correlation TPDM (Q|vn)", nQ, nvirB, nso_));
    G->contract(false, true, nQ * nvirB, nso_, nvirB, Gmo, CvirB, 1.0, 0.0);
    Gmo.reset();
    Gao->contract233(false, false, nso_, nso_, CvirB, G, 1.0, 1.0);
    G.reset();

    // symmetrize : This is necessary since we only consider OV block
    Gao->symmetrize3(Gao);
    Gao->write(psio_, PSIF_DFOCC_DENS, true, true);

    // 2-Index TPDM
    bQso = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|mn)", nQ, nso_, nso_));
    bQso->read(psio_, PSIF_DFOCC_INTS, true, true);
    Jmhalf = SharedTensor2d(new Tensor2d("DF_BASIS_CC Jmhalf <P|Q>", nQ, nQ));
    Jmhalf->read(psio_, PSIF_DFOCC_INTS);
    cQso = SharedTensor2d(new Tensor2d("DF_BASIS_CC C (Q|mn)", nQ, nso_, nso_));
    cQso->gemm(true, false, Jmhalf, bQso, 1.0, 0.0);
    bQso.reset();
    Jmhalf.reset();

    // Packed c_mn^Q
    cQso2 = SharedTensor2d(new Tensor2d("DF_BASIS_CC C (Q|m>=n)", nQ, ntri_so));
    cQso2->symm_packed(cQso);
    cQso.reset();

    // LTM of G_mn^Q
    Gao2 = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|n>=n)", nQ, ntri_so));
    Gao2->ltm(Gao);
    Gao.reset();

    // G_PQ = 1/2 \sum_{mn} c_mn^P G_mn^Q = 1/2 \sum_{m>=n} c_mn^P G_mn^Q (2 - \delta_{mn})
    G = SharedTensor2d(new Tensor2d("2-Index Correlation TPDM (P|Q)", nQ, nQ));
    G->gemm(false, true, cQso2, Gao2, 0.5, 0.0);
    Gao2.reset();
    cQso2.reset();
    G->write_symm(psio_, PSIF_DFOCC_DENS);
    G.reset();

}// else if (reference_ == "UNRESTRICTED")
    timer_off("back_trans");
    //outfile->Printf("\tBacktransformation is done.\n");
    //
} // end back_trans_cc


}} // End Namespaces
