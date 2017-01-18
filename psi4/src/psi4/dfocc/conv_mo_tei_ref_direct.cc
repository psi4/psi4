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

using namespace psi;
using namespace std;

namespace psi{ namespace dfoccwave{

void DFOCC::tei_oooo_chem_ref_directAA(SharedTensor2d &K)
{
    timer_on("Build (OO|OO)");
    bQooA = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|OO)", nQ_ref, noccA * noccA));
    bQooA->read(psio_, PSIF_DFOCC_INTS);
    K->gemm(true, false, bQooA, bQooA, 1.0, 0.0);
    bQooA.reset();
    timer_off("Build (OO|OO)");
}

//=======================================================
//          (oo|oo)
//=======================================================
void DFOCC::tei_oooo_chem_ref_directBB(SharedTensor2d &K)
{
    timer_on("Build (oo|oo)");
    bQooB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|oo)", nQ_ref, noccB * noccB));
    bQooB->read(psio_, PSIF_DFOCC_INTS);
    K->gemm(true, false, bQooB, bQooB, 1.0, 0.0);
    bQooB.reset();
    timer_off("Build (oo|oo)");
}

//=======================================================
//          (OO|oo)
//=======================================================
void DFOCC::tei_oooo_chem_ref_directAB(SharedTensor2d &K)
{
    timer_on("Build (OO|oo)");
    bQooA = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|OO)", nQ_ref, noccA * noccA));
    bQooB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|oo)", nQ_ref, noccB * noccB));
    bQooA->read(psio_, PSIF_DFOCC_INTS);
    bQooB->read(psio_, PSIF_DFOCC_INTS);
    K->gemm(true, false, bQooA, bQooB, 1.0, 0.0);
    bQooA.reset();
    bQooB.reset();
    timer_off("Build (OO|oo)");
}

//=======================================================
//          (OO|OV)
//=======================================================
void DFOCC::tei_ooov_chem_ref_directAA(SharedTensor2d &K)
{
    timer_on("Build (OO|OV)");
    bQooA = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|OO)", nQ_ref, noccA * noccA));
    bQovA = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|OV)", nQ_ref, noccA * nvirA));
    bQooA->read(psio_, PSIF_DFOCC_INTS);
    bQovA->read(psio_, PSIF_DFOCC_INTS);
    K->gemm(true, false, bQooA, bQovA, 1.0, 0.0);
    bQooA.reset();
    bQovA.reset();
    timer_off("Build (OO|OV)");
}

//=======================================================
//          (oo|ov)
//=======================================================
void DFOCC::tei_ooov_chem_ref_directBB(SharedTensor2d &K)
{
    timer_on("Build (oo|ov)");
    bQooB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|oo)", nQ_ref, noccB * noccB));
    bQovB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|ov)", nQ_ref, noccB * nvirB));
    bQooB->read(psio_, PSIF_DFOCC_INTS);
    bQovB->read(psio_, PSIF_DFOCC_INTS);
    K->gemm(true, false, bQooB, bQovB, 1.0, 0.0);
    bQooB.reset();
    bQovB.reset();
    timer_off("Build (oo|ov)");
}

//=======================================================
//          (OO|ov)
//=======================================================
void DFOCC::tei_ooov_chem_ref_directAB(SharedTensor2d &K)
{
    timer_on("Build (OO|ov)");
    bQooA = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|OO)", nQ_ref, noccA * noccA));
    bQovB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|ov)", nQ_ref, noccB * nvirB));
    bQooA->read(psio_, PSIF_DFOCC_INTS);
    bQovB->read(psio_, PSIF_DFOCC_INTS);
    K->gemm(true, false, bQooA, bQovB, 1.0, 0.0);
    bQooA.reset();
    bQovB.reset();
    timer_off("Build (OO|ov)");
}

//=======================================================
//          (OV|oo)
//=======================================================
void DFOCC::tei_ovoo_chem_ref_directAB(SharedTensor2d &K)
{
    timer_on("Build (OV|oo)");
    bQooB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|oo)", nQ_ref, noccB * noccB));
    bQovA = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|OV)", nQ_ref, noccA * nvirA));
    bQooB->read(psio_, PSIF_DFOCC_INTS);
    bQovA->read(psio_, PSIF_DFOCC_INTS);
    K->gemm(true, false, bQovA, bQooB, 1.0, 0.0);
    bQooB.reset();
    bQovA.reset();
    timer_off("Build (OV|oo)");
}

//=======================================================
//          (OO|VV)
//=======================================================
void DFOCC::tei_oovv_chem_ref_directAA(SharedTensor2d &K)
{
    timer_on("Build (OO|VV)");
    bQooA = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|OO)", nQ_ref, noccA * noccA));
    bQvvA = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|VV)", nQ_ref, nvirA, nvirA));
    bQooA->read(psio_, PSIF_DFOCC_INTS);
    bQvvA->read(psio_, PSIF_DFOCC_INTS, true, true);
    K->gemm(true, false, bQooA, bQvvA, 1.0, 0.0);
    bQooA.reset();
    bQvvA.reset();
    timer_off("Build (OO|VV)");
}

//=======================================================
//          (oo|vv)
//=======================================================
void DFOCC::tei_oovv_chem_ref_directBB(SharedTensor2d &K)
{
    timer_on("Build (oo|vv)");
    bQooB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|oo)", nQ_ref, noccB * noccB));
    bQvvB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|vv)", nQ_ref, nvirB, nvirB));
    bQooB->read(psio_, PSIF_DFOCC_INTS);
    bQvvB->read(psio_, PSIF_DFOCC_INTS, true, true);
    K->gemm(true, false, bQooB, bQvvB, 1.0, 0.0);
    timer_off("Build (oo|vv)");
}

//=======================================================
//          (OO|vv)
//=======================================================
void DFOCC::tei_oovv_chem_ref_directAB(SharedTensor2d &K)
{
    timer_on("Build (OO|vv)");
    bQooA = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|OO)", nQ_ref, noccA * noccA));
    bQvvB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|vv)", nQ_ref, nvirB, nvirB));
    bQooA->read(psio_, PSIF_DFOCC_INTS);
    bQvvB->read(psio_, PSIF_DFOCC_INTS, true, true);
    K->gemm(true, false, bQooA, bQvvB, 1.0, 0.0);
    bQooA.reset();
    bQvvB.reset();
    timer_off("Build (OO|vv)");
}

//=======================================================
//          (VV|oo)
//=======================================================
void DFOCC::tei_vvoo_chem_ref_directAB(SharedTensor2d &K)
{
    timer_on("Build (VV|oo)");
    bQooB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|oo)", nQ_ref, noccB * noccB));
    bQvvA = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|VV)", nQ_ref, nvirA, nvirA));
    bQooB->read(psio_, PSIF_DFOCC_INTS);
    bQvvA->read(psio_, PSIF_DFOCC_INTS, true, true);
    K->gemm(true, false, bQvvA, bQooB, 1.0, 0.0);
    bQooB.reset();
    bQvvA.reset();
    timer_off("Build (VV|oo)");
}

//=======================================================
//          (OV|OV)
//=======================================================
void DFOCC::tei_ovov_chem_ref_directAA(SharedTensor2d &K)
{
    timer_on("Build (OV|OV)");
    bQovA = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|OV)", nQ_ref, noccA * nvirA));
    bQovA->read(psio_, PSIF_DFOCC_INTS);
    K->gemm(true, false, bQovA, bQovA, 1.0, 0.0);
    bQovA.reset();
    timer_off("Build (OV|OV)");
}

//=======================================================
//          (ov|ov)
//=======================================================
void DFOCC::tei_ovov_chem_ref_directBB(SharedTensor2d &K)
{
    timer_on("Build (ov|ov)");
    bQovB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|ov)", nQ_ref, noccB * nvirB));
    bQovB->read(psio_, PSIF_DFOCC_INTS);
    K->gemm(true, false, bQovB, bQovB, 1.0, 0.0);
    bQovB.reset();
    timer_off("Build (ov|ov)");
}

//=======================================================
//          (OV|ov)
//=======================================================
void DFOCC::tei_ovov_chem_ref_directAB(SharedTensor2d &K)
{
    timer_on("Build (OV|ov)");
    bQovA = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|OV)", nQ_ref, noccA * nvirA));
    bQovB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|ov)", nQ_ref, noccB * nvirB));
    bQovA->read(psio_, PSIF_DFOCC_INTS);
    bQovB->read(psio_, PSIF_DFOCC_INTS);
    K->gemm(true, false, bQovA, bQovB, 1.0, 0.0);
    bQovA.reset();
    bQovB.reset();
    timer_off("Build (OV|ov)");
}

//=======================================================
//          (VO|VO)
//=======================================================
void DFOCC::tei_vovo_chem_ref_directAA(SharedTensor2d &K)
{
    timer_on("Build (VO|VO)");
    bQovA = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|OV)", nQ_ref, noccA, nvirA));
    bQovA->read(psio_, PSIF_DFOCC_INTS);
    SharedTensor2d bQvo = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|VO)", nQ_ref, nvirA, noccA));
    bQvo->swap_3index_col(bQovA);
    bQovA.reset();
    K->gemm(true, false, bQvo, bQvo, 1.0, 0.0);
    bQvo.reset();
    timer_off("Build (VO|VO)");
}

//=======================================================
//          (vo|vo)
//=======================================================
void DFOCC::tei_vovo_chem_ref_directBB(SharedTensor2d &K)
{
    timer_on("Build (vo|vo)");
    bQovB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|ov)", nQ_ref, noccB, nvirB));
    bQovB->read(psio_, PSIF_DFOCC_INTS);
    SharedTensor2d bQvo = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|vo)", nQ_ref, nvirB, noccB));
    bQvo->swap_3index_col(bQovB);
    bQovB.reset();
    K->gemm(true, false, bQvo, bQvo, 1.0, 0.0);
    bQvo.reset();
    timer_off("Build (vo|vo)");
}

//=======================================================
//          (VO|vo)
//=======================================================
void DFOCC::tei_vovo_chem_ref_directAB(SharedTensor2d &K)
{
    timer_on("Build (VO|vo)");
    bQovA = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|OV)", nQ_ref, noccA, nvirA));
    bQovA->read(psio_, PSIF_DFOCC_INTS);
    SharedTensor2d bQvoA = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|VO)", nQ_ref, nvirA, noccA));
    bQvoA->swap_3index_col(bQovA);
    bQovA.reset();

    bQovB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|ov)", nQ_ref, noccB, nvirB));
    bQovB->read(psio_, PSIF_DFOCC_INTS);
    SharedTensor2d bQvoB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|vo)", nQ_ref, nvirB, noccB));
    bQvoB->swap_3index_col(bQovB);
    bQovB.reset();

    K->gemm(true, false, bQvoA, bQvoB, 1.0, 0.0);
    bQvoA.reset();
    bQvoB.reset();
    timer_off("Build (VO|vo)");
}

//=======================================================
//          <OO|OO>
//=======================================================
void DFOCC::tei_oooo_phys_ref_directAA(SharedTensor2d &K)
{
    timer_on("Build <OO|OO>");
    SharedTensor2d L = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (OO|OO)", noccA, noccA, noccA, noccA));
    tei_oooo_chem_ref_directAA(L);
    K->sort(1324, L, 1.0, 0.0);
    L.reset();
    timer_off("Build <OO|OO>");
}

//=======================================================
//          <oo|oo>
//=======================================================
void DFOCC::tei_oooo_phys_ref_directBB(SharedTensor2d &K)
{
    timer_on("Build <oo|oo>");
    SharedTensor2d L = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (oo|oo)", noccB, noccB, noccB, noccB));
    tei_oooo_chem_ref_directBB(L);
    K->sort(1324, L, 1.0, 0.0);
    L.reset();
    timer_off("Build <oo|oo>");
}

//=======================================================
//          <Oo|Oo>
//=======================================================
void DFOCC::tei_oooo_phys_ref_directAB(SharedTensor2d &K)
{
    timer_on("Build <Oo|Oo>");
    SharedTensor2d L = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (OO|oo)", noccA, noccA, noccB, noccB));
    tei_oooo_chem_ref_directAB(L);
    K->sort(1324, L, 1.0, 0.0);
    L.reset();
    timer_off("Build <Oo|Oo>");
}

//=======================================================
//          <OO|OV>
//=======================================================
void DFOCC::tei_ooov_phys_ref_directAA(SharedTensor2d &K)
{
    timer_on("Build <OO|OV>");
    SharedTensor2d L = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (OO|OV)", noccA, noccA, noccA, nvirA));
    tei_ooov_chem_ref_directAA(L);
    K->sort(1324, L, 1.0, 0.0);
    L.reset();
    timer_off("Build <OO|OV>");
}

//=======================================================
//          <oo|ov>
//=======================================================
void DFOCC::tei_ooov_phys_ref_directBB(SharedTensor2d &K)
{
    timer_on("Build <oo|ov>");
    SharedTensor2d L = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (oo|ov)", noccB, noccB, noccB, nvirB));
    tei_ooov_chem_ref_directBB(L);
    K->sort(1324, L, 1.0, 0.0);
    L.reset();
    timer_off("Build <oo|ov>");
}

//=======================================================
//          <Oo|Ov>
//=======================================================
void DFOCC::tei_ooov_phys_ref_directAB(SharedTensor2d &K)
{
    timer_on("Build <Oo|Ov>");
    SharedTensor2d L = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (OO|ov)", noccA, noccA, noccB, nvirB));
    tei_ooov_chem_ref_directAB(L);
    K->sort(1324, L, 1.0, 0.0);
    L.reset();
    timer_off("Build <Oo|Ov>");
}

//=======================================================
//          <Oo|Vo>
//=======================================================
void DFOCC::tei_oovo_phys_ref_directAB(SharedTensor2d &K)
{
    timer_on("Build <Oo|Vo>");
    SharedTensor2d L = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (OV|oo)", noccA, nvirA, noccB, noccB));
    tei_ovoo_chem_ref_directAB(L);
    K->sort(1324, L, 1.0, 0.0);
    L.reset();
    timer_off("Build <Oo|Vo>");
}

//=======================================================
//          <OO|VV>
//=======================================================
void DFOCC::tei_oovv_phys_ref_directAA(SharedTensor2d &K)
{
    timer_on("Build <OO|VV>");
    SharedTensor2d L = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (OV|OV)", noccA, nvirA, noccA, nvirA));
    tei_ovov_chem_ref_directAA(L);
    K->sort(1324, L, 1.0, 0.0);
    L.reset();
    timer_off("Build <OO|VV>");
}

//=======================================================
//          <oo|vv>
//=======================================================
void DFOCC::tei_oovv_phys_ref_directBB(SharedTensor2d &K)
{
    timer_on("Build <oo|vv>");
    SharedTensor2d L = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (ov|ov)", noccB, nvirB, noccB, nvirB));
    tei_ovov_chem_ref_directBB(L);
    K->sort(1324, L, 1.0, 0.0);
    L.reset();
    timer_off("Build <oo|vv>");
}

//=======================================================
//          <Oo|Vv>
//=======================================================
void DFOCC::tei_oovv_phys_ref_directAB(SharedTensor2d &K)
{
    timer_on("Build <Oo|Vv>");
    SharedTensor2d L = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (OV|ov)", noccA, nvirA, noccB, nvirB));
    tei_ovov_chem_ref_directAB(L);
    K->sort(1324, L, 1.0, 0.0);
    L.reset();
    timer_off("Build <Oo|Vv>");
}

//=======================================================
//          <OV|OV>
//=======================================================
void DFOCC::tei_ovov_phys_ref_directAA(SharedTensor2d &K)
{
    timer_on("Build <OV|OV>");
    SharedTensor2d L = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (OO|VV)", noccA, noccA, nvirA, nvirA));
    tei_oovv_chem_ref_directAA(L);
    K->sort(1324, L, 1.0, 0.0);
    L.reset();
    timer_off("Build <OV|OV>");
}

//=======================================================
//          <ov|ov>
//=======================================================
void DFOCC::tei_ovov_phys_ref_directBB(SharedTensor2d &K)
{
    timer_on("Build <ov|ov>");
    SharedTensor2d L = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (oo|vv)", noccB, noccB, nvirB, nvirB));
    tei_oovv_chem_ref_directBB(L);
    K->sort(1324, L, 1.0, 0.0);
    L.reset();
    timer_off("Build <ov|ov>");
}

//=======================================================
//          <Ov|Ov>
//=======================================================
void DFOCC::tei_ovov_phys_ref_directAB(SharedTensor2d &K)
{
    timer_on("Build <Ov|Ov>");
    SharedTensor2d L = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (OO|vv)", noccA, noccA, nvirB, nvirB));
    tei_oovv_chem_ref_directAB(L);
    K->sort(1324, L, 1.0, 0.0);
    L.reset();
    timer_off("Build <Ov|Ov>");
}

//=======================================================
//          <Vo|Vo>
//=======================================================
void DFOCC::tei_vovo_phys_ref_directAB(SharedTensor2d &K)
{
    timer_on("Build <Vo|Vo>");
    SharedTensor2d L = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (VV|oo)", nvirA, nvirA, noccB, noccB));
    tei_vvoo_chem_ref_directAB(L);
    K->sort(1324, L, 1.0, 0.0);
    L.reset();
    timer_off("Build <Vo|Vo>");
}

}} // End Namespaces
