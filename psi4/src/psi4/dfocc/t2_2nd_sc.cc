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

#include "psi4/libqt/qt.h"
#include "defines.h"
#include "dfocc.h"

using namespace psi;
using namespace std;


namespace psi{ namespace dfoccwave{

void DFOCC::t2_2nd_sc()
{

    // defs
    SharedTensor2d K, L, M, I, T, Tnew, T1, T2, U, Tau, W, X, Y;

if (reference_ == "RESTRICTED") {
    // Read DF integrals
    bQijA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|IJ)", nQ, naoccA, naoccA));
    bQiaA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|IA)", nQ, naoccA, navirA));
    bQabA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|AB)", nQ, navirA, navirA));
    bQijA->read(psio_, PSIF_DFOCC_INTS);
    bQiaA->read(psio_, PSIF_DFOCC_INTS);
    bQabA->read(psio_, PSIF_DFOCC_INTS, true, true);

    // Write for the general code
    Tnew = SharedTensor2d(new Tensor2d("New T2_2 (IA|JB)", naoccA, navirA, naoccA, navirA));
    Tnew->write_symm(psio_, PSIF_DFOCC_AMPS);
    Tnew.reset();

    // Read T2_1
    t2 = SharedTensor2d(new Tensor2d("T2_1 (IA|JB)", naoccA, navirA, naoccA, navirA));
    t2->read_symm(psio_, PSIF_DFOCC_AMPS);

    // WmnijT2
    mp3_WmnijT2();

    // WmbejT2
    mp3_WmbejT2();

    // WabefT2
    mp3_WabefT2();

    // Denom
    Tnew = SharedTensor2d(new Tensor2d("New T2_2 (IA|JB)", naoccA, navirA, naoccA, navirA));
    Tnew->read_symm(psio_, PSIF_DFOCC_AMPS);
    Tnew->apply_denom_chem(nfrzc, noccA, FockA);
    //Tnew->print();

    // Form T2 = T2(1) + T2(2) : MP3
    if (wfn_type_ == "DF-OMP3" || wfn_type_ == "CD-OMP3") t2->axpy(Tnew, 1.0);
    // Form T2 = T2(1) + 1/2 T2(2) : MP2.5
    else if (wfn_type_ == "DF-OMP2.5" || wfn_type_ == "CD-OMP2.5") t2->axpy(Tnew, 0.5);

    // Write
    if (dertype == "FIRST" || orb_opt_ == "TRUE") {
	// Write T2_2
        T = SharedTensor2d(new Tensor2d("T2_2 (IA|JB)", naoccA, navirA, naoccA, navirA));
	T->copy(Tnew);
        T->write_symm(psio_, PSIF_DFOCC_AMPS);
        T.reset();

	// Write T2
        T = SharedTensor2d(new Tensor2d("T2 (IA|JB)", naoccA, navirA, naoccA, navirA));
	T->copy(t2);
        T->write_symm(psio_, PSIF_DFOCC_AMPS);
        T.reset();
    }
    Tnew.reset();

    // Form U(ia,jb) = 2*T(ia,jb) - T (ib,ja)
    U = SharedTensor2d(new Tensor2d("U2 (IA|JB)", naoccA, navirA, naoccA, navirA));
    ccsd_u2_amps(U,t2);
    if (dertype == "FIRST" || orb_opt_ == "TRUE") U->write_symm(psio_, PSIF_DFOCC_AMPS);
    t2.reset();

    // Energy
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (IA|JB)", naoccA, navirA, naoccA, navirA));
    K->gemm(true, false, bQiaA, bQiaA, 1.0, 0.0);
    Ecorr = U->vector_dot(K);
    U.reset();
    K.reset();
    Emp3 = Eref + Ecorr;

    // Free ints
    bQijA.reset();
    bQiaA.reset();
    bQabA.reset();
}// end if (reference_ == "RESTRICTED")

else if (reference_ == "UNRESTRICTED") {
    // Read DF integrals
    bQijA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|IJ)", nQ, naoccA, naoccA));
    bQijB = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|ij)", nQ, naoccB, naoccB));
    bQiaA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|IA)", nQ, naoccA, navirA));
    bQiaB = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|ia)", nQ, naoccB, navirB));
    bQabA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|AB)", nQ, navirA, navirA));
    bQabB = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|ab)", nQ, navirB, navirB));
    bQijA->read(psio_, PSIF_DFOCC_INTS);
    bQijB->read(psio_, PSIF_DFOCC_INTS);
    bQiaA->read(psio_, PSIF_DFOCC_INTS);
    bQiaB->read(psio_, PSIF_DFOCC_INTS);
    bQabA->read(psio_, PSIF_DFOCC_INTS, true, true);
    bQabB->read(psio_, PSIF_DFOCC_INTS, true, true);

    //=========================
    // T2AA
    //=========================

    // WmnijT2
    mp3_WmnijT2AA();

    // WmbejT2
    mp3_WmbejT2AA();

    // WabefT2
    mp3_WabefT2AA();

    // Denom
    Tnew = SharedTensor2d(new Tensor2d("New T2_2 <IJ|AB>", naoccA, naoccA, navirA, navirA));
    Tnew->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
    Tnew->apply_denom(nfrzc, noccA, FockA);

    // Form T2 = T2(1) + T2(2) : MP3
    U = SharedTensor2d(new Tensor2d("T2_1 <IJ|AB>", naoccA, naoccA, navirA, navirA));
    U->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
    if (wfn_type_ == "DF-OMP3" || wfn_type_ == "CD-OMP3") U->axpy(Tnew, 1.0);
    // Form T2 = T2(1) + 1/2 T2(2) : MP2.5
    else if (wfn_type_ == "DF-OMP2.5" || wfn_type_ == "CD-OMP2.5") U->axpy(Tnew, 0.5);

    // Write
    if (dertype == "FIRST" || orb_opt_ == "TRUE") {
	// Write T2_2
        T = SharedTensor2d(new Tensor2d("T2_2 <IJ|AB>", naoccA, naoccA, navirA, navirA));
	T->copy(Tnew);
        T->write_anti_symm(psio_, PSIF_DFOCC_AMPS);
        T.reset();

	// Write T2
        T = SharedTensor2d(new Tensor2d("T2 <IJ|AB>", naoccA, naoccA, navirA, navirA));
	T->copy(U);
        T->write_anti_symm(psio_, PSIF_DFOCC_AMPS);
        T.reset();
    }
    Tnew.reset();

    // Energy
    L = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (IA|JB)", naoccA, navirA, naoccA, navirA));
    L->gemm(true, false, bQiaA, bQiaA, 1.0, 0.0);
    M = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <IJ|AB>", naoccA, naoccA, navirA, navirA));
    M->sort(1324, L, 1.0, 0.0);
    L.reset();
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <IJ||AB>", naoccA, naoccA, navirA, navirA));
    tei_pqrs_anti_symm_direct(K, M);
    M.reset();
    Emp3AA = 0.25 * U->vector_dot(K);
    U.reset();
    K.reset();

    //=========================
    // T2BB
    //=========================

    // WmnijT2
    mp3_WmnijT2BB();

    // WmbejT2
    mp3_WmbejT2BB();

    // WabefT2
    mp3_WabefT2BB();

    // Denom
    Tnew = SharedTensor2d(new Tensor2d("New T2_2 <ij|ab>", naoccB, naoccB, navirB, navirB));
    Tnew->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
    Tnew->apply_denom(nfrzc, noccB, FockB);

    // Form T2 = T2(1) + T2(2): MP3
    U = SharedTensor2d(new Tensor2d("T2_1 <ij|ab>", naoccB, naoccB, navirB, navirB));
    U->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
    if (wfn_type_ == "DF-OMP3" || wfn_type_ == "CD-OMP3") U->axpy(Tnew, 1.0);
    // Form T2 = T2(1) + 1/2 T2(2) : MP2.5
    else if (wfn_type_ == "DF-OMP2.5" || wfn_type_ == "CD-OMP2.5") U->axpy(Tnew, 0.5);

    // Write
    if (dertype == "FIRST" || orb_opt_ == "TRUE") {
	// Write T2_2
        T = SharedTensor2d(new Tensor2d("T2_2 <ij|ab>", naoccB, naoccB, navirB, navirB));
	T->copy(Tnew);
        T->write_anti_symm(psio_, PSIF_DFOCC_AMPS);
        T.reset();

	// Write T2
        T = SharedTensor2d(new Tensor2d("T2 <ij|ab>", naoccB, naoccB, navirB, navirB));
	T->copy(U);
        T->write_anti_symm(psio_, PSIF_DFOCC_AMPS);
        T.reset();
    }
    Tnew.reset();

    // Energy BB
    L = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (ia|jb)", naoccB, navirB, naoccB, navirB));
    L->gemm(true, false, bQiaB, bQiaB, 1.0, 0.0);
    M = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <ij|ab>", naoccB, naoccB, navirB, navirB));
    M->sort(1324, L, 1.0, 0.0);
    L.reset();
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <ij||ab>", naoccB, naoccB, navirB, navirB));
    tei_pqrs_anti_symm_direct(K, M);
    M.reset();
    Emp3BB = 0.25 * U->vector_dot(K);
    U.reset();
    K.reset();

    //=========================
    // T2AB
    //=========================

    // WmnijT2
    mp3_WmnijT2AB();

    // WmbejT2
    mp3_WmbejT2AB();

    // WabefT2
    mp3_WabefT2AB();

    // Denom
    Tnew = SharedTensor2d(new Tensor2d("New T2_2 <Ij|Ab>", naoccA, naoccB, navirA, navirB));
    Tnew->read(psio_, PSIF_DFOCC_AMPS);
    Tnew->apply_denom_os(nfrzc, noccA, noccB, FockA, FockB);

    // Form T2 = T2(1) + T2(2) : MP3
    U = SharedTensor2d(new Tensor2d("T2_1 <Ij|Ab>", naoccA, naoccB, navirA, navirB));
    U->read(psio_, PSIF_DFOCC_AMPS);
    if (wfn_type_ == "DF-OMP3" || wfn_type_ == "CD-OMP3") U->axpy(Tnew, 1.0);
    // Form T2 = T2(1) + 1/2 T2(2) : MP2.5
    else if (wfn_type_ == "DF-OMP2.5" || wfn_type_ == "CD-OMP2.5") U->axpy(Tnew, 0.5);

    // Write
    if (dertype == "FIRST" || orb_opt_ == "TRUE") {
	// Write T2_2
        T = SharedTensor2d(new Tensor2d("T2_2 <Ij|Ab>", naoccA, naoccB, navirA, navirB));
	T->copy(Tnew);
        T->write(psio_, PSIF_DFOCC_AMPS);
        T.reset();

	// Write T2
        T = SharedTensor2d(new Tensor2d("T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB));
	T->copy(U);
        T->write(psio_, PSIF_DFOCC_AMPS);
        T.reset();
    }
    Tnew.reset();

    // Energy
    L = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (IA|jb)", naoccA, navirA, naoccB, navirB));
    L->gemm(true, false, bQiaA, bQiaB, 1.0, 0.0);
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <Ij|Ab>", naoccA, naoccB, navirA, navirB));
    K->sort(1324, L, 1.0, 0.0);
    L.reset();
    Emp3AB = U->vector_dot(K);
    U.reset();
    K.reset();

    // Overall energy
    Ecorr = Emp3AA + Emp3BB + Emp3AB;
    Emp3 = Escf + Ecorr;

    // Free ints
    bQijA.reset();
    bQijB.reset();
    bQiaA.reset();
    bQiaB.reset();
    bQabA.reset();
    bQabB.reset();

}// else if (reference_ == "UNRESTRICTED")

}// end t2_2nd_sc

}} // End Namespaces
