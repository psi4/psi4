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

void DFOCC::t2_rmp2_direct(SharedTensor2d& T)
{
    SharedTensor2d K;
    timer_on("T2_MP2");
    // Build amplitudes in Mulliken order
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (IA|JB)", naoccA, navirA, naoccA, navirA));
    tei_iajb_chem_directAA(K);
    T->copy(K);
    T->apply_denom_chem(nfrzc, noccA, FockA);
    timer_off("T2_MP2");
} // end t2_rmp2_direct

//=======================================================
// U(ia,jb) = 2*T(ia,jb) - T(ib,ja): T(ia,jb)= T_ij^ab
//=======================================================
void DFOCC::u2_rmp2_direct(SharedTensor2d& T, SharedTensor2d& U)
{
    SharedTensor2d K;
    timer_on("T2_MP2");
    // Build amplitudes in Mulliken order
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (IA|JB)", naoccA, navirA, naoccA, navirA));
    tei_iajb_chem_directAA(K);
    T->copy(K);
    T->apply_denom_chem(nfrzc, noccA, FockA);

    // form U(ia,jb)
    U->sort(1432, T, 1.0, 0.0);
    U->scale(-1.0);
    U->axpy(T, 2.0);
    timer_off("T2_MP2");
} // end u2_rmp2_direct

//=======================================================
// U(ia,jb) = 2*T(ia,jb) - T(ib,ja): T(ia,jb)= T_ij^ab
//=======================================================
void DFOCC::u2_rmp2_direct(SharedTensor2d& U)
{
    SharedTensor2d K, T;
    timer_on("T2_MP2");
    // Build amplitudes in Mulliken order
    T = SharedTensor2d(new Tensor2d("T2_1 (ia|jb)", naoccA, navirA, naoccA, navirA));
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (IA|JB)", naoccA, navirA, naoccA, navirA));
    tei_iajb_chem_directAA(K);
    T->copy(K);
    T->apply_denom_chem(nfrzc, noccA, FockA);

    // form U(ia,jb)
    U->sort(1432, T, 1.0, 0.0);
    U->scale(-1.0);
    U->axpy(T, 2.0);
    T.reset();
    timer_off("T2_MP2");
} // end u2_rmp2_direct

//=======================================================
//          T2AA
//=======================================================
void DFOCC::t2AA_ump2_direct(SharedTensor2d& T)
{
    SharedTensor2d K, L, M;
    timer_on("T2AA_MP2");
    L = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (IA|JB)", naoccA, navirA, naoccA, navirA));
    tei_iajb_chem_directAA(L);
    M = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <IJ|AB>", naoccA, naoccA, navirA, navirA));
    M->sort(1324, L, 1.0, 0.0);
    L.reset();
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <IJ||AB>", naoccA, naoccA, navirA, navirA));
    tei_pqrs_anti_symm_direct(K, M);
    M.reset();
    T->copy(K);
    T->apply_denom(nfrzc, noccA, FockA);
    timer_off("T2AA_MP2");
} // end t2AA_ump2_direct

//=======================================================
//          T2BB
//=======================================================
void DFOCC::t2BB_ump2_direct(SharedTensor2d& T)
{
    SharedTensor2d K, L, M;
    timer_on("T2BB_MP2");
    L = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (ia|jb)", naoccB, navirB, naoccB, navirB));
    tei_iajb_chem_directBB(L);
    M = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <ij|ab>", naoccB, naoccB, navirB, navirB));
    M->sort(1324, L, 1.0, 0.0);
    L.reset();
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <ij||ab>", naoccB, naoccB, navirB, navirB));
    tei_pqrs_anti_symm_direct(K, M);
    M.reset();
    T->copy(K);
    T->apply_denom(nfrzc, noccB, FockB);
    timer_off("T2BB_MP2");
} // end t2BB_ump2_direct

//=======================================================
//          T2AB
//=======================================================
void DFOCC::t2AB_ump2_direct(SharedTensor2d& T)
{
    SharedTensor2d K, L;
    timer_on("T2AB_MP2");
    L = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (IA|jb)", naoccA, navirA, naoccB, navirB));
    tei_iajb_chem_directAB(L);
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <Ij|Ab>", naoccA, naoccB, navirA, navirB));
    K->sort(1324, L, 1.0, 0.0);
    L.reset();
    T->copy(K);
    T->apply_denom_os(nfrzc, noccA, noccB, FockA, FockB);
    timer_off("T2AB_MP2");
} // end t2AB_ump2_direct

}} // End Namespaces
