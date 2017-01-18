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

void DFOCC::tei_ijkl_chem_directAA(SharedTensor2d &K)
{
    timer_on("Build (IJ|KL)");
    bQijA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|IJ)", nQ, naoccA * naoccA));
    bQijA->read(psio_, PSIF_DFOCC_INTS);
    K->gemm(true, false, bQijA, bQijA, 1.0, 0.0);
    bQijA.reset();
    timer_off("Build (IJ|KL)");
}

//=======================================================
//          (ij|kl)
//=======================================================
void DFOCC::tei_ijkl_chem_directBB(SharedTensor2d &K)
{
    timer_on("Build (ij|kl)");
    bQijB = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|ij)", nQ, naoccB * naoccB));
    bQijB->read(psio_, PSIF_DFOCC_INTS);
    K->gemm(true, false, bQijB, bQijB, 1.0, 0.0);
    bQijB.reset();
    timer_off("Build (ij|kl)");
}

//=======================================================
//          (IJ|kl)
//=======================================================
void DFOCC::tei_ijkl_chem_directAB(SharedTensor2d &K)
{
    timer_on("Build (IJ|kl)");
    bQijA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|IJ)", nQ, naoccA * naoccA));
    bQijB = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|ij)", nQ, naoccB * naoccB));
    bQijA->read(psio_, PSIF_DFOCC_INTS);
    bQijB->read(psio_, PSIF_DFOCC_INTS);
    K->gemm(true, false, bQijA, bQijB, 1.0, 0.0);
    bQijA.reset();
    bQijB.reset();
    timer_off("Build (IJ|kl)");
}

//=======================================================
//          (OO|OO)
//=======================================================
void DFOCC::tei_oooo_chem_directAA(SharedTensor2d &K)
{
    timer_on("Build (OO|OO)");
    bQooA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|OO)", nQ, noccA * noccA));
    bQooA->read(psio_, PSIF_DFOCC_INTS);
    K->gemm(true, false, bQooA, bQooA, 1.0, 0.0);
    bQooA.reset();
    timer_off("Build (OO|OO)");
}

//=======================================================
//          (oo|oo)
//=======================================================
void DFOCC::tei_oooo_chem_directBB(SharedTensor2d &K)
{
    timer_on("Build (oo|oo)");
    bQooB = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|oo)", nQ, noccB * noccB));
    bQooB->read(psio_, PSIF_DFOCC_INTS);
    K->gemm(true, false, bQooB, bQooB, 1.0, 0.0);
    bQooB.reset();
    timer_off("Build (oo|oo)");
}

//=======================================================
//          (OO|oo)
//=======================================================
void DFOCC::tei_oooo_chem_directAB(SharedTensor2d &K)
{
    timer_on("Build (OO|oo)");
    bQooA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|OO)", nQ, noccA * noccA));
    bQooB = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|oo)", nQ, noccB * noccB));
    bQooA->read(psio_, PSIF_DFOCC_INTS);
    bQooB->read(psio_, PSIF_DFOCC_INTS);
    K->gemm(true, false, bQooA, bQooB, 1.0, 0.0);
    bQooA.reset();
    bQooB.reset();
    timer_off("Build (OO|oo)");
}

//=======================================================
//          (IJ|KA)
//=======================================================
void DFOCC::tei_ijka_chem_directAA(SharedTensor2d &K)
{
    timer_on("Build (IJ|KA)");
    bQijA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|IJ)", nQ, naoccA * naoccA));
    bQiaA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|IA)", nQ, naoccA * navirA));
    bQijA->read(psio_, PSIF_DFOCC_INTS);
    bQiaA->read(psio_, PSIF_DFOCC_INTS);
    K->gemm(true, false, bQijA, bQiaA, 1.0, 0.0);
    bQijA.reset();
    bQiaA.reset();
    timer_off("Build (IJ|KA)");
}

//=======================================================
//          (ij|ka)
//=======================================================
void DFOCC::tei_ijka_chem_directBB(SharedTensor2d &K)
{
    timer_on("Build (ij|ka)");
    bQijB = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|ij)", nQ, naoccB * naoccB));
    bQiaB = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|ia)", nQ, naoccB * navirB));
    bQijB->read(psio_, PSIF_DFOCC_INTS);
    bQiaB->read(psio_, PSIF_DFOCC_INTS);
    K->gemm(true, false, bQijB, bQiaB, 1.0, 0.0);
    bQijB.reset();
    bQiaB.reset();
    timer_off("Build (ij|ka)");
}

//=======================================================
//          (IJ|ka)
//=======================================================
void DFOCC::tei_ijka_chem_directAB(SharedTensor2d &K)
{
    timer_on("Build (IJ|ka)");
    bQijA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|IJ)", nQ, naoccA * naoccA));
    bQiaB = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|ia)", nQ, naoccB * navirB));
    bQijA->read(psio_, PSIF_DFOCC_INTS);
    bQiaB->read(psio_, PSIF_DFOCC_INTS);
    K->gemm(true, false, bQijA, bQiaB, 1.0, 0.0);
    bQijA.reset();
    bQiaB.reset();
    timer_off("Build (IJ|ka)");
}

//=======================================================
//          (IA|jk)
//=======================================================
void DFOCC::tei_iajk_chem_directAB(SharedTensor2d &K)
{
    timer_on("Build (IA|jk)");
    bQijB = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|ij)", nQ, naoccB * naoccB));
    bQiaA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|IA)", nQ, naoccA * navirA));
    bQijB->read(psio_, PSIF_DFOCC_INTS);
    bQiaA->read(psio_, PSIF_DFOCC_INTS);
    K->gemm(true, false, bQiaA, bQijB, 1.0, 0.0);
    bQijB.reset();
    bQiaA.reset();
    timer_off("Build (IA|jk)");
}

//=======================================================
//          (OO|OV)
//=======================================================
void DFOCC::tei_ooov_chem_directAA(SharedTensor2d &K)
{
    timer_on("Build (OO|OV)");
    bQooA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|OO)", nQ, noccA * noccA));
    bQovA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|OV)", nQ, noccA * nvirA));
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
void DFOCC::tei_ooov_chem_directBB(SharedTensor2d &K)
{
    timer_on("Build (oo|ov)");
    bQooB = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|oo)", nQ, noccB * noccB));
    bQovB = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|ov)", nQ, noccB * nvirB));
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
void DFOCC::tei_ooov_chem_directAB(SharedTensor2d &K)
{
    timer_on("Build (OO|ov)");
    bQooA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|OO)", nQ, noccA * noccA));
    bQovB = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|ov)", nQ, noccB * nvirB));
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
void DFOCC::tei_ovoo_chem_directAB(SharedTensor2d &K)
{
    timer_on("Build (OV|oo)");
    bQooB = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|oo)", nQ, noccB * noccB));
    bQovA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|OV)", nQ, noccA * nvirA));
    bQooB->read(psio_, PSIF_DFOCC_INTS);
    bQovA->read(psio_, PSIF_DFOCC_INTS);
    K->gemm(true, false, bQovA, bQooB, 1.0, 0.0);
    bQooB.reset();
    bQovA.reset();
    timer_off("Build (OV|oo)");
}

//=======================================================
//          (IJ|AB)
//=======================================================
void DFOCC::tei_ijab_chem_directAA(SharedTensor2d &K)
{
    timer_on("Build (IJ|AB)");
    bQijA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|IJ)", nQ, naoccA, naoccA));
    bQabA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|AB)", nQ, navirA, navirA));
    bQijA->read(psio_, PSIF_DFOCC_INTS);
    bQabA->read(psio_, PSIF_DFOCC_INTS, true, true);
    K->gemm(true, false, bQijA, bQabA, 1.0, 0.0);
    bQijA.reset();
    bQabA.reset();
    timer_off("Build (IJ|AB)");
}

//=======================================================
//          (ij|ab)
//=======================================================
void DFOCC::tei_ijab_chem_directBB(SharedTensor2d &K)
{
    timer_on("Build (ij|ab)");
    bQijB = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|ij)", nQ, naoccB, naoccB));
    bQabB = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|ab)", nQ, navirB, navirB));
    bQijB->read(psio_, PSIF_DFOCC_INTS);
    bQabB->read(psio_, PSIF_DFOCC_INTS, true, true);
    K->gemm(true, false, bQijB, bQabB, 1.0, 0.0);
    bQijB.reset();
    bQabB.reset();
    timer_off("Build (ij|ab)");
}

//=======================================================
//          (IJ|ab)
//=======================================================
void DFOCC::tei_ijab_chem_directAB(SharedTensor2d &K)
{
    timer_on("Build (IJ|ab)");
    bQijA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|IJ)", nQ, naoccA, naoccA));
    bQabB = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|ab)", nQ, navirB, navirB));
    bQijA->read(psio_, PSIF_DFOCC_INTS);
    bQabB->read(psio_, PSIF_DFOCC_INTS, true, true);
    K->gemm(true, false, bQijA, bQabB, 1.0, 0.0);
    bQijA.reset();
    bQabB.reset();
    timer_off("Build (IJ|ab)");
}

//=======================================================
//          (AB|ij)
//=======================================================
void DFOCC::tei_abij_chem_directAB(SharedTensor2d &K)
{
    timer_on("Build (AB|ij)");
    bQijB = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|ij)", nQ, naoccB * naoccB));
    bQabA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|AB)", nQ, navirA * navirA));
    bQijB->read(psio_, PSIF_DFOCC_INTS);
    bQabA->read(psio_, PSIF_DFOCC_INTS, true, true);
    K->gemm(true, false, bQabA, bQijB, 1.0, 0.0);
    bQijB.reset();
    bQabA.reset();
    timer_off("Build (AB|ij)");
}

//=======================================================
//          (OO|VV)
//=======================================================
void DFOCC::tei_oovv_chem_directAA(SharedTensor2d &K)
{
    timer_on("Build (OO|VV)");
    bQooA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|OO)", nQ, noccA * noccA));
    bQvvA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|VV)", nQ, nvirA, nvirA));
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
void DFOCC::tei_oovv_chem_directBB(SharedTensor2d &K)
{
    timer_on("Build (oo|vv)");
    bQooB = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|oo)", nQ, noccB * noccB));
    bQvvB = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|vv)", nQ, nvirB, nvirB));
    bQooB->read(psio_, PSIF_DFOCC_INTS);
    bQvvB->read(psio_, PSIF_DFOCC_INTS, true, true);
    K->gemm(true, false, bQooB, bQvvB, 1.0, 0.0);
    timer_off("Build (oo|vv)");
}

//=======================================================
//          (OO|vv)
//=======================================================
void DFOCC::tei_oovv_chem_directAB(SharedTensor2d &K)
{
    timer_on("Build (OO|vv)");
    bQooA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|OO)", nQ, noccA * noccA));
    bQvvB = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|vv)", nQ, nvirB, nvirB));
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
void DFOCC::tei_vvoo_chem_directAB(SharedTensor2d &K)
{
    timer_on("Build (VV|oo)");
    bQooB = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|oo)", nQ, noccB * noccB));
    bQvvA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|VV)", nQ, nvirA, nvirA));
    bQooB->read(psio_, PSIF_DFOCC_INTS);
    bQvvA->read(psio_, PSIF_DFOCC_INTS, true, true);
    K->gemm(true, false, bQvvA, bQooB, 1.0, 0.0);
    bQooB.reset();
    bQvvA.reset();
    timer_off("Build (VV|oo)");
}

//=======================================================
//          (IA|JB)
//=======================================================
void DFOCC::tei_iajb_chem_directAA(SharedTensor2d &K)
{
    timer_on("Build (IA|JB)");
    bQiaA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|IA)", nQ, naoccA, navirA));
    bQiaA->read(psio_, PSIF_DFOCC_INTS);
    K->gemm(true, false, bQiaA, bQiaA, 1.0, 0.0);
    bQiaA.reset();
    timer_off("Build (IA|JB)");
}

//=======================================================
//          (ia|jb)
//=======================================================
void DFOCC::tei_iajb_chem_directBB(SharedTensor2d &K)
{
    timer_on("Build (ia|jb)");
    bQiaB = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|ia)", nQ, naoccB * navirB));
    bQiaB->read(psio_, PSIF_DFOCC_INTS);
    K->gemm(true, false, bQiaB, bQiaB, 1.0, 0.0);
    bQiaB.reset();
    timer_off("Build (ia|jb)");
}

//=======================================================
//          (IA|jb)
//=======================================================
void DFOCC::tei_iajb_chem_directAB(SharedTensor2d &K)
{
    timer_on("Build (IA|jb)");
    bQiaA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|IA)", nQ, naoccA * navirA));
    bQiaB = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|ia)", nQ, naoccB * navirB));
    bQiaA->read(psio_, PSIF_DFOCC_INTS);
    bQiaB->read(psio_, PSIF_DFOCC_INTS);
    K->gemm(true, false, bQiaA, bQiaB, 1.0, 0.0);
    bQiaA.reset();
    bQiaB.reset();
    timer_off("Build (IA|jb)");
}

//=======================================================
//          (OV|OV)
//=======================================================
void DFOCC::tei_ovov_chem_directAA(SharedTensor2d &K)
{
    timer_on("Build (OV|OV)");
    bQovA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|OV)", nQ, noccA * nvirA));
    bQovA->read(psio_, PSIF_DFOCC_INTS);
    K->gemm(true, false, bQovA, bQovA, 1.0, 0.0);
    bQovA.reset();
    timer_off("Build (OV|OV)");
}

//=======================================================
//          (ov|ov)
//=======================================================
void DFOCC::tei_ovov_chem_directBB(SharedTensor2d &K)
{
    timer_on("Build (ov|ov)");
    bQovB = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|ov)", nQ, noccB * nvirB));
    bQovB->read(psio_, PSIF_DFOCC_INTS);
    K->gemm(true, false, bQovB, bQovB, 1.0, 0.0);
    bQovB.reset();
    timer_off("Build (ov|ov)");
}

//=======================================================
//          (OV|ov)
//=======================================================
void DFOCC::tei_ovov_chem_directAB(SharedTensor2d &K)
{
    timer_on("Build (OV|ov)");
    bQovA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|OV)", nQ, noccA * nvirA));
    bQovB = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|ov)", nQ, noccB * nvirB));
    bQovA->read(psio_, PSIF_DFOCC_INTS);
    bQovB->read(psio_, PSIF_DFOCC_INTS);
    K->gemm(true, false, bQovA, bQovB, 1.0, 0.0);
    bQovA.reset();
    bQovB.reset();
    timer_off("Build (OV|ov)");
}

//=======================================================
//          <IJ|KL>
//=======================================================
void DFOCC::tei_ijkl_phys_directAA(SharedTensor2d &K)
{
    timer_on("Build <IJ|KL>");
    SharedTensor2d L = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (IJ|KL)", naoccA, naoccA, naoccA, naoccA));
    tei_ijkl_chem_directAA(L);
    K->sort(1324, L, 1.0, 0.0);
    L.reset();
    timer_off("Build <IJ|KL>");
}

//=======================================================
//          <ij|kl>
//=======================================================
void DFOCC::tei_ijkl_phys_directBB(SharedTensor2d &K)
{
    timer_on("Build <ij|kl>");
    SharedTensor2d L = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (ij|kl)", naoccB, naoccB, naoccB, naoccB));
    tei_ijkl_chem_directBB(L);
    K->sort(1324, L, 1.0, 0.0);
    L.reset();
    timer_off("Build <ij|kl>");
}

//=======================================================
//          <Ij|Kl>
//=======================================================
void DFOCC::tei_ijkl_phys_directAB(SharedTensor2d &K)
{
    timer_on("Build <Ij|Kl>");
    SharedTensor2d L = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (IJ|kl)", naoccA, naoccA, naoccB, naoccB));
    tei_ijkl_chem_directAB(L);
    K->sort(1324, L, 1.0, 0.0);
    L.reset();
    timer_off("Build <Ij|Kl>");
}

//=======================================================
//          <OO|OO>
//=======================================================
void DFOCC::tei_oooo_phys_directAA(SharedTensor2d &K)
{
    timer_on("Build <OO|OO>");
    SharedTensor2d L = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (OO|OO)", noccA, noccA, noccA, noccA));
    tei_oooo_chem_directAA(L);
    K->sort(1324, L, 1.0, 0.0);
    L.reset();
    timer_off("Build <OO|OO>");
}

//=======================================================
//          <oo|oo>
//=======================================================
void DFOCC::tei_oooo_phys_directBB(SharedTensor2d &K)
{
    timer_on("Build <oo|oo>");
    SharedTensor2d L = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (oo|oo)", noccB, noccB, noccB, noccB));
    tei_oooo_chem_directBB(L);
    K->sort(1324, L, 1.0, 0.0);
    L.reset();
    timer_off("Build <oo|oo>");
}

//=======================================================
//          <Oo|Oo>
//=======================================================
void DFOCC::tei_oooo_phys_directAB(SharedTensor2d &K)
{
    timer_on("Build <Oo|Oo>");
    SharedTensor2d L = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (OO|oo)", noccA, noccA, noccB, noccB));
    tei_oooo_chem_directAB(L);
    K->sort(1324, L, 1.0, 0.0);
    L.reset();
    timer_off("Build <Oo|Oo>");
}

//=======================================================
//          <IJ|KA>
//=======================================================
void DFOCC::tei_ijka_phys_directAA(SharedTensor2d &K)
{
    timer_on("Build <IJ|KA>");
    SharedTensor2d L = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (IJ|KA)", naoccA, naoccA, naoccA, navirA));
    tei_ijka_chem_directAA(L);
    K->sort(1324, L, 1.0, 0.0);
    L.reset();
    timer_off("Build <IJ|KA>");
}

//=======================================================
//          <ij|ka>
//=======================================================
void DFOCC::tei_ijka_phys_directBB(SharedTensor2d &K)
{
    timer_on("Build <ij|ka>");
    SharedTensor2d L = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (ij|ka)", naoccB, naoccB, naoccB, navirB));
    tei_ijka_chem_directBB(L);
    K->sort(1324, L, 1.0, 0.0);
    L.reset();
    timer_off("Build <ij|ka>");
}

//=======================================================
//          <Ij|Ka>
//=======================================================
void DFOCC::tei_ijka_phys_directAB(SharedTensor2d &K)
{
    timer_on("Build <Ij|Ka>");
    SharedTensor2d L = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (IJ|ka)", naoccA, naoccA, naoccB, navirB));
    tei_ijka_chem_directAB(L);
    K->sort(1324, L, 1.0, 0.0);
    L.reset();
    timer_off("Build <Ij|Ka>");
}

//=======================================================
//          <Ij|Ak>
//=======================================================
void DFOCC::tei_ijak_phys_directAB(SharedTensor2d &K)
{
    timer_on("Build <Ij|Ak>");
    SharedTensor2d L = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (IA|jk)", naoccA, navirA, naoccB, naoccB));
    tei_iajk_chem_directAB(L);
    K->sort(1324, L, 1.0, 0.0);
    L.reset();
    timer_off("Build <Ij|Ak>");
}

//=======================================================
//          <OO|OV>
//=======================================================
void DFOCC::tei_ooov_phys_directAA(SharedTensor2d &K)
{
    timer_on("Build <OO|OV>");
    SharedTensor2d L = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (OO|OV)", noccA, noccA, noccA, nvirA));
    tei_ooov_chem_directAA(L);
    K->sort(1324, L, 1.0, 0.0);
    L.reset();
    timer_off("Build <OO|OV>");
}

//=======================================================
//          <oo|ov>
//=======================================================
void DFOCC::tei_ooov_phys_directBB(SharedTensor2d &K)
{
    timer_on("Build <oo|ov>");
    SharedTensor2d L = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (oo|ov)", noccB, noccB, noccB, nvirB));
    tei_ooov_chem_directBB(L);
    K->sort(1324, L, 1.0, 0.0);
    L.reset();
    timer_off("Build <oo|ov>");
}

//=======================================================
//          <Oo|Ov>
//=======================================================
void DFOCC::tei_ooov_phys_directAB(SharedTensor2d &K)
{
    timer_on("Build <Oo|Ov>");
    SharedTensor2d L = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (OO|ov)", noccA, noccA, noccB, nvirB));
    tei_ooov_chem_directAB(L);
    K->sort(1324, L, 1.0, 0.0);
    L.reset();
    timer_off("Build <Oo|Ov>");
}

//=======================================================
//          <Oo|Vo>
//=======================================================
void DFOCC::tei_oovo_phys_directAB(SharedTensor2d &K)
{
    timer_on("Build <Oo|Vo>");
    SharedTensor2d L = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (OV|oo)", noccA, nvirA, noccB, noccB));
    tei_ovoo_chem_directAB(L);
    K->sort(1324, L, 1.0, 0.0);
    L.reset();
    timer_off("Build <Oo|Vo>");
}

//=======================================================
//          <IJ|AB>
//=======================================================
void DFOCC::tei_ijab_phys_directAA(SharedTensor2d &K)
{
    timer_on("Build <IJ|AB>");
    SharedTensor2d L = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (IA|JB)", naoccA, navirA, naoccA, navirA));
    tei_iajb_chem_directAA(L);
    K->sort(1324, L, 1.0, 0.0);
    L.reset();
    timer_off("Build <IJ|AB>");
}

//=======================================================
//          <ij|ab>
//=======================================================
void DFOCC::tei_ijab_phys_directBB(SharedTensor2d &K)
{
    timer_on("Build <ij|ab>");
    SharedTensor2d L = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (ia|jb)", naoccB, navirB, naoccB, navirB));
    tei_iajb_chem_directBB(L);
    K->sort(1324, L, 1.0, 0.0);
    L.reset();
    timer_off("Build <ij|ab>");
}

//=======================================================
//          <Ij|Ab>
//=======================================================
void DFOCC::tei_ijab_phys_directAB(SharedTensor2d &K)
{
    timer_on("Build <Ij|Ab>");
    SharedTensor2d L = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (IA|jb)", naoccA, navirA, naoccB, navirB));
    tei_iajb_chem_directAB(L);
    K->sort(1324, L, 1.0, 0.0);
    L.reset();
    timer_off("Build <Ij|Ab>");
}

//=======================================================
//          <OO|VV>
//=======================================================
void DFOCC::tei_oovv_phys_directAA(SharedTensor2d &K)
{
    timer_on("Build <OO|VV>");
    SharedTensor2d L = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (OV|OV)", noccA, nvirA, noccA, nvirA));
    tei_ovov_chem_directAA(L);
    K->sort(1324, L, 1.0, 0.0);
    L.reset();
    timer_off("Build <OO|VV>");
}

//=======================================================
//          <oo|vv>
//=======================================================
void DFOCC::tei_oovv_phys_directBB(SharedTensor2d &K)
{
    timer_on("Build <oo|vv>");
    SharedTensor2d L = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (ov|ov)", noccB, nvirB, noccB, nvirB));
    tei_ovov_chem_directBB(L);
    K->sort(1324, L, 1.0, 0.0);
    L.reset();
    timer_off("Build <oo|vv>");
}

//=======================================================
//          <Oo|Vv>
//=======================================================
void DFOCC::tei_oovv_phys_directAB(SharedTensor2d &K)
{
    timer_on("Build <Oo|Vv>");
    SharedTensor2d L = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (OV|ov)", noccA, nvirA, noccB, nvirB));
    tei_ovov_chem_directAB(L);
    K->sort(1324, L, 1.0, 0.0);
    L.reset();
    timer_off("Build <Oo|Vv>");
}

//=======================================================
//          <IA|JB>
//=======================================================
void DFOCC::tei_iajb_phys_directAA(SharedTensor2d &K)
{
    timer_on("Build <IA|JB>");
    SharedTensor2d L = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (IJ|AB)", naoccA, naoccA, navirA, navirA));
    tei_ijab_chem_directAA(L);
    K->sort(1324, L, 1.0, 0.0);
    L.reset();
    timer_off("Build <IA|JB>");
}

//=======================================================
//          <ia|jb>
//=======================================================
void DFOCC::tei_iajb_phys_directBB(SharedTensor2d &K)
{
    timer_on("Build <ia|jb>");
    SharedTensor2d L = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (ij|ab)", naoccB, naoccB, navirB, navirB));
    tei_ijab_chem_directBB(L);
    K->sort(1324, L, 1.0, 0.0);
    L.reset();
    timer_off("Build <ia|jb>");
}

//=======================================================
//          <Ia|Jb>
//=======================================================
void DFOCC::tei_iajb_phys_directAB(SharedTensor2d &K)
{
    timer_on("Build <Ia|Jb>");
    SharedTensor2d L = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (IJ|ab)", naoccA, naoccA, navirB, navirB));
    tei_ijab_chem_directAB(L);
    K->sort(1324, L, 1.0, 0.0);
    L.reset();
    timer_off("Build <Ia|Jb>");
}

//=======================================================
//          <Ai|Bj>
//=======================================================
void DFOCC::tei_aibj_phys_directAB(SharedTensor2d &K)
{
    timer_on("Build <Ai|Bj>");
    SharedTensor2d L = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (AB|ij)", navirA, navirA, naoccB, naoccB));
    tei_abij_chem_directAB(L);
    K->sort(1324, L, 1.0, 0.0);
    L.reset();
    timer_off("Build <Ai|Bj>");
}

//=======================================================
//          <OV|OV>
//=======================================================
void DFOCC::tei_ovov_phys_directAA(SharedTensor2d &K)
{
    timer_on("Build <OV|OV>");
    SharedTensor2d L = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (OO|VV)", noccA, noccA, nvirA, nvirA));
    tei_oovv_chem_directAA(L);
    K->sort(1324, L, 1.0, 0.0);
    L.reset();
    timer_off("Build <OV|OV>");
}

//=======================================================
//          <ov|ov>
//=======================================================
void DFOCC::tei_ovov_phys_directBB(SharedTensor2d &K)
{
    timer_on("Build <ov|ov>");
    SharedTensor2d L = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (oo|vv)", noccB, noccB, nvirB, nvirB));
    tei_oovv_chem_directBB(L);
    K->sort(1324, L, 1.0, 0.0);
    L.reset();
    timer_off("Build <ov|ov>");
}

//=======================================================
//          <Ov|Ov>
//=======================================================
void DFOCC::tei_ovov_phys_directAB(SharedTensor2d &K)
{
    timer_on("Build <Ov|Ov>");
    SharedTensor2d L = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (OO|vv)", noccA, noccA, nvirB, nvirB));
    tei_oovv_chem_directAB(L);
    K->sort(1324, L, 1.0, 0.0);
    L.reset();
    timer_off("Build <Ov|Ov>");
}

//=======================================================
//          <Vo|Vo>
//=======================================================
void DFOCC::tei_vovo_phys_directAB(SharedTensor2d &K)
{
    timer_on("Build <Vo|Vo>");
    SharedTensor2d L = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (VV|oo)", nvirA, nvirA, noccB, noccB));
    tei_vvoo_chem_directAB(L);
    K->sort(1324, L, 1.0, 0.0);
    L.reset();
    timer_off("Build <Vo|Vo>");
}

//=======================================================
//          <PQ||RS> : <IJ||KL>, <IJ||AB>
//=======================================================
void DFOCC::tei_pqrs_anti_symm_direct(SharedTensor2d &K, SharedTensor2d &L)
{
    timer_on("Build <PQ||RS>");
    // K = <PQ||RS>; L = <PQ|RS>
    // <PQ||RS> = <PQ|RS> - <PQ|SR>
    K->sort(1243, L, 1.0, 0.0);
    K->scale(-1.0);
    K->add(L);
    L.reset();
    timer_off("Build <PQ||RS>");
}

//=======================================================
//          <PQ||RS> : <IJ||KA>
//=======================================================
void DFOCC::tei_pqrs2_anti_symm_direct(SharedTensor2d &K, SharedTensor2d &L)
{
    timer_on("Build <PQ||RS>");
    // K = <PQ||RS>; L = <PQ|RS>
    // <PQ||RS> = <PQ|RS> - <QP|RS>
    K->sort(2134, L, 1.0, 0.0);
    K->scale(-1.0);
    K->add(L);
    L.reset();
    timer_off("Build <PQ||RS>");
}

//=======================================================
//          <PQ||RS>
//=======================================================
void DFOCC::tei_pqrs3_anti_symm_direct(SharedTensor2d &K, SharedTensor2d &L, SharedTensor2d &M)
{
    timer_on("Build <PQ||RS>");
    // K = <PQ||RS>; L = <PQ|RS>, M = (PS|RQ)
    // <PQ||RS> = <PQ|RS> - (PS|RQ)
    K->sort(1432, M, 1.0, 0.0);
    M.reset();
    K->scale(-1.0);
    K->add(L);
    L.reset();
    timer_off("Build <PQ||RS>");
}

//=======================================================
//          (PQ|RS) : General
//=======================================================
void DFOCC::tei_chem_direct(SharedTensor2d &K, SharedTensor2d &L, SharedTensor2d &M)
{
    timer_on("Build (PQ|RS)");
    // K = (pq|rs); L = B_pq^Q, M = B_rs^Q
    K->gemm(true, false, L, M, 1.0, 0.0);
    timer_off("Build (PQ|RS)");
}

//=======================================================
//          <PQ|RS> : General
//=======================================================
void DFOCC::tei_phys_direct(SharedTensor2d &I, SharedTensor2d &K, SharedTensor2d &L, SharedTensor2d &M)
{
    timer_on("Build <PQ|RS>");
    // K = (pr|qs); L = B_pr^Q, M = B_qs^Q
    K->gemm(true, false, L, M, 1.0, 0.0);
    I->sort(1324, K, 1.0, 0.0);
    timer_off("Build <PQ|RS>");
}

//=======================================================
//      Closed-Shell Anti-Symmetrized Integrals
//      L(pq,rs) = 2 <pq|rs> - <pq|sr>
//=======================================================
void DFOCC::tei_cs1_anti_symm_direct(SharedTensor2d &I, SharedTensor2d &J, SharedTensor2d &K)
{
    I->sort(1243, K, -1.0, 0.0);
    I->axpy(J, 2.0);
}

//=======================================================
//      Closed-Shell Anti-Symmetrized Integrals
//      L(pq,rs) = 2 <pq|rs> - <qp|rs>
//=======================================================
void DFOCC::tei_cs2_anti_symm_direct(SharedTensor2d &I, SharedTensor2d &J, SharedTensor2d &K)
{
    I->sort(2134, K, -1.0, 0.0);
    I->axpy(J, 2.0);
}

//=======================================================
//      Closed-Shell Anti-Symmetrized Integrals
//      L(pq,rs) = 2 (pq|rs) - (ps|rq)
//=======================================================
void DFOCC::tei_cs3_anti_symm_direct(SharedTensor2d &I, SharedTensor2d &J, SharedTensor2d &K)
{
    I->sort(1432, K, -1.0, 0.0);
    I->axpy(J, 2.0);
}

//=======================================================
//      Closed-Shell Anti-Symmetrized Integrals
//      L(pq,rs) = 2 (pq|rs) - (rq|ps)
//=======================================================
void DFOCC::tei_cs4_anti_symm_direct(SharedTensor2d &I, SharedTensor2d &J, SharedTensor2d &K)
{
    I->sort(3214, K, -1.0, 0.0);
    I->axpy(J, 2.0);
}

}} // End Namespaces
