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

void DFOCC::mp2l_energy()
{

    timer_on("CCL Energy");

    EcorrL = 0.0;
    double Eoei = 0.0;
    double Eoo = 0.0;
    double Eov = 0.0;
    double Evv = 0.0;
    Emp2L_old = Emp2L;

if (reference_ == "RESTRICTED") {
    // DE = \sum_{p,q} G_pq f_pq
    EcorrL += G1c->vector_dot(FockA);
    Eoei = EcorrL;

    // DE += \sum_{Q} \sum_{i,a} G_ia^Q b_ia^Q
    G2c_ov = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|OV)", nQ, noccA * nvirA));
    bQovA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|OV)", nQ, noccA * nvirA));
    G2c_ov->read(psio_, PSIF_DFOCC_DENS);
    bQovA->read(psio_, PSIF_DFOCC_INTS);
    EcorrL += G2c_ov->vector_dot(bQovA);
    G2c_ov.reset();
    bQovA.reset();
    Eov = EcorrL - Eoei;

}// end if (reference_ == "RESTRICTED")

else if (reference_ == "UNRESTRICTED") {
    // DE = \sum_{p,q} G_pq f_pq
    EcorrL += G1cA->vector_dot(FockA);
    EcorrL += G1cB->vector_dot(FockB);
    Eoei = EcorrL;

    // DE += \sum_{Q} \sum_{I,A} G_IA^Q b_IA^Q
    G2c_ovA = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|OV)", nQ, noccA * nvirA));
    bQovA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|OV)", nQ, noccA * nvirA));
    G2c_ovA->read(psio_, PSIF_DFOCC_DENS);
    bQovA->read(psio_, PSIF_DFOCC_INTS);
    EcorrL += G2c_ovA->vector_dot(bQovA);
    G2c_ovA.reset();
    bQovA.reset();

    // DE += \sum_{Q} \sum_{i,a} G_ia^Q b_ia^Q
    G2c_ovB = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|ov)", nQ, noccB * nvirB));
    bQovB = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|ov)", nQ, noccB * nvirB));
    G2c_ovB->read(psio_, PSIF_DFOCC_DENS);
    bQovB->read(psio_, PSIF_DFOCC_INTS);
    EcorrL += G2c_ovB->vector_dot(bQovB);
    G2c_ovB.reset();
    bQovB.reset();
    Eov = EcorrL - Eoei;

}// else if (reference_ == "UNRESTRICTED")

    Emp2L = Eref + EcorrL;
    DE = Emp2L - Emp2L_old;

    /*
    outfile->Printf("\tDF-MP2L One-Electron Energy (a.u.) : %20.14f\n", Eoei);
    outfile->Printf("\tDF-MP2L OV Energy (a.u.)           : %20.14f\n", Eov);
    outfile->Printf("\tDF-MP2L Correlation Energy (a.u.)  : %20.14f\n", EcorrL);
    outfile->Printf("\tDF-MP2L Total Energy (a.u.)        : %20.14f\n", Emp2L);
    */

    timer_off("CCL Energy");
} // end mp2l_energy

//=======================================================
//          MP3-L Energy
//=======================================================
void DFOCC::mp3l_energy()
{

    timer_on("CCL Energy");

    SharedTensor2d G, K;

    EcorrL = 0.0;
    double Eoei = 0.0;
    double Eoo = 0.0;
    double Eov = 0.0;
    double Evv = 0.0;
    Emp3L_old = Emp3L;

if (reference_ == "RESTRICTED") {
    // DE = \sum_{p,q} G_pq f_pq
    EcorrL += G1c->vector_dot(FockA);
    Eoei = EcorrL;

    // DE += 1/2 \sum_{Q} \sum_{i,j} G_ij^Q b_ij^Q
    G = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|OO)", nQ, noccA, noccA));
    G->read(psio_, PSIF_DFOCC_DENS);
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|OO)", nQ, noccA, noccA));
    K->read(psio_, PSIF_DFOCC_INTS);
    EcorrL += 0.5 * G->vector_dot(K);
    G.reset();
    K.reset();
    Eoo = EcorrL - Eoei;

    // DE += \sum_{Q} \sum_{i,a} G_ia^Q b_ia^Q
    G = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|OV)", nQ, noccA, nvirA));
    G->read(psio_, PSIF_DFOCC_DENS);
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|OV)", nQ, noccA, nvirA));
    K->read(psio_, PSIF_DFOCC_INTS);
    EcorrL += G->vector_dot(K);
    G.reset();
    K.reset();
    Eov = EcorrL - Eoei - Eoo;

    // DE += 1/2 \sum_{Q} \sum_{a,b} G_ab^Q b_ab^Q
    G = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|VV)", nQ, nvirA, nvirA));
    G->read(psio_, PSIF_DFOCC_DENS, true, true);
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|VV)", nQ, nvirA, nvirA));
    K->read(psio_, PSIF_DFOCC_INTS, true, true);
    EcorrL += 0.5 * G->vector_dot(K);
    G.reset();
    K.reset();
    Evv = EcorrL - Eoei - Eoo - Eov;

}// end if (reference_ == "RESTRICTED")

else if (reference_ == "UNRESTRICTED") {
    // DE = \sum_{p,q} G_pq f_pq
    EcorrL += G1cA->vector_dot(FockA);
    EcorrL += G1cB->vector_dot(FockB);
    Eoei = EcorrL;

    // DE += 1/2 \sum_{Q} \sum_{I,J} G_IJ^Q b_IJ^Q
    G = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|OO)", nQ, noccA, noccA));
    G->read(psio_, PSIF_DFOCC_DENS);
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|OO)", nQ, noccA, noccA));
    K->read(psio_, PSIF_DFOCC_INTS);
    EcorrL += 0.5 * G->vector_dot(K);
    G.reset();
    K.reset();

    // DE += 1/2 \sum_{Q} \sum_{i,j} G_ij^Q b_ij^Q
    G = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|oo)", nQ, noccB, noccB));
    G->read(psio_, PSIF_DFOCC_DENS);
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|oo)", nQ, noccB, noccB));
    K->read(psio_, PSIF_DFOCC_INTS);
    EcorrL += 0.5 * G->vector_dot(K);
    G.reset();
    K.reset();
    Eoo = EcorrL - Eoei;

    // DE += \sum_{Q} \sum_{I,A} G_IA^Q b_IA^Q
    G2c_ovA = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|OV)", nQ, noccA * nvirA));
    bQovA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|OV)", nQ, noccA * nvirA));
    G2c_ovA->read(psio_, PSIF_DFOCC_DENS);
    bQovA->read(psio_, PSIF_DFOCC_INTS);
    EcorrL += G2c_ovA->vector_dot(bQovA);
    G2c_ovA.reset();
    bQovA.reset();

    // DE += \sum_{Q} \sum_{i,a} G_ia^Q b_ia^Q
    G2c_ovB = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|ov)", nQ, noccB * nvirB));
    bQovB = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|ov)", nQ, noccB * nvirB));
    G2c_ovB->read(psio_, PSIF_DFOCC_DENS);
    bQovB->read(psio_, PSIF_DFOCC_INTS);
    EcorrL += G2c_ovB->vector_dot(bQovB);
    G2c_ovB.reset();
    bQovB.reset();
    Eov = EcorrL - Eoei - Eoo;

    // DE += 1/2 \sum_{Q} \sum_{A,B} G_AB^Q b_AB^Q
    G = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|VV)", nQ, nvirA, nvirA));
    G->read(psio_, PSIF_DFOCC_DENS, true, true);
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|VV)", nQ, nvirA, nvirA));
    K->read(psio_, PSIF_DFOCC_INTS, true, true);
    EcorrL += 0.5 * G->vector_dot(K);
    G.reset();
    K.reset();

    // DE += 1/2 \sum_{Q} \sum_{a,b} G_ab^Q b_ab^Q
    G = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|vv)", nQ, nvirB, nvirB));
    G->read(psio_, PSIF_DFOCC_DENS, true, true);
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|vv)", nQ, nvirB, nvirB));
    K->read(psio_, PSIF_DFOCC_INTS, true, true);
    EcorrL += 0.5 * G->vector_dot(K);
    G.reset();
    K.reset();
    Evv = EcorrL - Eoei - Eoo - Eov;

}// else if (reference_ == "UNRESTRICTED")

    Emp3L = Eref + EcorrL;
    DE = Emp3L - Emp3L_old;

    /*
    outfile->Printf("\n\tEnergies re-computed from CC density: \n");
    outfile->Printf("\t------------------------------------- \n");
    outfile->Printf("\tReference Energy (a.u.)            : %20.14f\n", Eref);
    outfile->Printf("\tOne-Electron Energy (a.u.)         : %20.14f\n", Eoei);
    outfile->Printf("\tOO Energy (a.u.)                   : %20.14f\n", Eoo);
    outfile->Printf("\tOV Energy (a.u.)                   : %20.14f\n", Eov);
    outfile->Printf("\tVV Energy (a.u.)                   : %20.14f\n", Evv);
    outfile->Printf("\tLagrangian Corr. Energy (a.u.)     : %20.14f\n", EcorrL);
    outfile->Printf("\tLagrangian Total Energy (a.u.)     : %20.14f\n", EccsdL);
    */

    timer_off("CCL Energy");
} // end mp3l_energy

//=======================================================
//          LCCD-L Energy
//=======================================================
void DFOCC::lccdl_energy()
{

    timer_on("CCL Energy");

    SharedTensor2d G, K;

    EcorrL = 0.0;
    double Eoei = 0.0;
    double Eoo = 0.0;
    double Eov = 0.0;
    double Evv = 0.0;
    ElccdL_old = ElccdL;

if (reference_ == "RESTRICTED") {
    // DE = \sum_{p,q} G_pq f_pq
    EcorrL += G1c->vector_dot(FockA);
    Eoei = EcorrL;

    // DE += 1/2 \sum_{Q} \sum_{i,j} G_ij^Q b_ij^Q
    G = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|OO)", nQ, noccA, noccA));
    G->read(psio_, PSIF_DFOCC_DENS);
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|OO)", nQ, noccA, noccA));
    K->read(psio_, PSIF_DFOCC_INTS);
    EcorrL += 0.5 * G->vector_dot(K);
    G.reset();
    K.reset();
    Eoo = EcorrL - Eoei;

    // DE += \sum_{Q} \sum_{i,a} G_ia^Q b_ia^Q
    G = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|OV)", nQ, noccA, nvirA));
    G->read(psio_, PSIF_DFOCC_DENS);
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|OV)", nQ, noccA, nvirA));
    K->read(psio_, PSIF_DFOCC_INTS);
    EcorrL += G->vector_dot(K);
    G.reset();
    K.reset();
    Eov = EcorrL - Eoei - Eoo;

    // DE += 1/2 \sum_{Q} \sum_{a,b} G_ab^Q b_ab^Q
    G = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|VV)", nQ, nvirA, nvirA));
    G->read(psio_, PSIF_DFOCC_DENS, true, true);
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|VV)", nQ, nvirA, nvirA));
    K->read(psio_, PSIF_DFOCC_INTS, true, true);
    EcorrL += 0.5 * G->vector_dot(K);
    G.reset();
    K.reset();
    Evv = EcorrL - Eoei - Eoo - Eov;

}// end if (reference_ == "RESTRICTED")

else if (reference_ == "UNRESTRICTED") {
    // DE = \sum_{p,q} G_pq f_pq
    EcorrL += G1cA->vector_dot(FockA);
    EcorrL += G1cB->vector_dot(FockB);
    Eoei = EcorrL;

    // DE += 1/2 \sum_{Q} \sum_{I,J} G_IJ^Q b_IJ^Q
    G = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|OO)", nQ, noccA, noccA));
    G->read(psio_, PSIF_DFOCC_DENS);
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|OO)", nQ, noccA, noccA));
    K->read(psio_, PSIF_DFOCC_INTS);
    EcorrL += 0.5 * G->vector_dot(K);
    G.reset();
    K.reset();

    // DE += 1/2 \sum_{Q} \sum_{i,j} G_ij^Q b_ij^Q
    G = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|oo)", nQ, noccB, noccB));
    G->read(psio_, PSIF_DFOCC_DENS);
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|oo)", nQ, noccB, noccB));
    K->read(psio_, PSIF_DFOCC_INTS);
    EcorrL += 0.5 * G->vector_dot(K);
    G.reset();
    K.reset();
    Eoo = EcorrL - Eoei;

    // DE += \sum_{Q} \sum_{I,A} G_IA^Q b_IA^Q
    G2c_ovA = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|OV)", nQ, noccA * nvirA));
    bQovA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|OV)", nQ, noccA * nvirA));
    G2c_ovA->read(psio_, PSIF_DFOCC_DENS);
    bQovA->read(psio_, PSIF_DFOCC_INTS);
    EcorrL += G2c_ovA->vector_dot(bQovA);
    G2c_ovA.reset();
    bQovA.reset();

    // DE += \sum_{Q} \sum_{i,a} G_ia^Q b_ia^Q
    G2c_ovB = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|ov)", nQ, noccB * nvirB));
    bQovB = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|ov)", nQ, noccB * nvirB));
    G2c_ovB->read(psio_, PSIF_DFOCC_DENS);
    bQovB->read(psio_, PSIF_DFOCC_INTS);
    EcorrL += G2c_ovB->vector_dot(bQovB);
    G2c_ovB.reset();
    bQovB.reset();
    Eov = EcorrL - Eoei - Eoo;

    // DE += 1/2 \sum_{Q} \sum_{A,B} G_AB^Q b_AB^Q
    G = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|VV)", nQ, nvirA, nvirA));
    G->read(psio_, PSIF_DFOCC_DENS, true, true);
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|VV)", nQ, nvirA, nvirA));
    K->read(psio_, PSIF_DFOCC_INTS, true, true);
    EcorrL += 0.5 * G->vector_dot(K);
    G.reset();
    K.reset();

    // DE += 1/2 \sum_{Q} \sum_{a,b} G_ab^Q b_ab^Q
    G = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|vv)", nQ, nvirB, nvirB));
    G->read(psio_, PSIF_DFOCC_DENS, true, true);
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|vv)", nQ, nvirB, nvirB));
    K->read(psio_, PSIF_DFOCC_INTS, true, true);
    EcorrL += 0.5 * G->vector_dot(K);
    G.reset();
    K.reset();
    Evv = EcorrL - Eoei - Eoo - Eov;

}// else if (reference_ == "UNRESTRICTED")

    ElccdL = Eref + EcorrL;
    DE = ElccdL - ElccdL_old;

    /*
    outfile->Printf("\n\tEnergies re-computed from CC density: \n");
    outfile->Printf("\t------------------------------------- \n");
    outfile->Printf("\tReference Energy (a.u.)            : %20.14f\n", Eref);
    outfile->Printf("\tOne-Electron Energy (a.u.)         : %20.14f\n", Eoei);
    outfile->Printf("\tOO Energy (a.u.)                   : %20.14f\n", Eoo);
    outfile->Printf("\tOV Energy (a.u.)                   : %20.14f\n", Eov);
    outfile->Printf("\tVV Energy (a.u.)                   : %20.14f\n", Evv);
    outfile->Printf("\tLagrangian Corr. Energy (a.u.)     : %20.14f\n", EcorrL);
    outfile->Printf("\tLagrangian Total Energy (a.u.)     : %20.14f\n", EccsdL);
    */

    timer_off("CCL Energy");
} // end lccdl_energy

//=======================================================
//          CCL Energy
//=======================================================
void DFOCC::ccl_energy()
{

    timer_on("CCL Energy");

    SharedTensor2d G, K;

    EcorrL = 0.0;
    double Eoei = 0.0;
    double Eoo = 0.0;
    double Eov = 0.0;
    double Evv = 0.0;

if (reference_ == "RESTRICTED") {
    // DE = \sum_{p,q} G_pq f_pq
    EcorrL += G1c->vector_dot(FockA);
    Eoei = EcorrL;

    // DE += 1/2 \sum_{Q} \sum_{i,j} G_ij^Q b_ij^Q
    G = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|OO)", nQ, noccA, noccA));
    G->read(psio_, PSIF_DFOCC_DENS);
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|OO)", nQ, noccA, noccA));
    K->read(psio_, PSIF_DFOCC_INTS);
    EcorrL += 0.5 * G->vector_dot(K);
    G.reset();
    K.reset();
    Eoo = EcorrL - Eoei;

    // DE += \sum_{Q} \sum_{i,a} G_ia^Q b_ia^Q
    G = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|OV)", nQ, noccA, nvirA));
    G->read(psio_, PSIF_DFOCC_DENS);
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|OV)", nQ, noccA, nvirA));
    K->read(psio_, PSIF_DFOCC_INTS);
    EcorrL += G->vector_dot(K);
    G.reset();
    K.reset();
    Eov = EcorrL - Eoei - Eoo;

    // DE += 1/2 \sum_{Q} \sum_{a,b} G_ab^Q b_ab^Q
    G = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|VV)", nQ, nvirA, nvirA));
    G->read(psio_, PSIF_DFOCC_DENS, true, true);
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|VV)", nQ, nvirA, nvirA));
    K->read(psio_, PSIF_DFOCC_INTS, true, true);
    EcorrL += 0.5 * G->vector_dot(K);
    G.reset();
    K.reset();
    Evv = EcorrL - Eoei - Eoo - Eov;

}// end if (reference_ == "RESTRICTED")

else if (reference_ == "UNRESTRICTED") {
    // DE = \sum_{p,q} G_pq f_pq
    EcorrL += G1cA->vector_dot(FockA);
    EcorrL += G1cB->vector_dot(FockB);
    Eoei = EcorrL;

    // DE += 1/2 \sum_{Q} \sum_{I,J} G_IJ^Q b_IJ^Q
    G = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|OO)", nQ, noccA, noccA));
    G->read(psio_, PSIF_DFOCC_DENS);
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|OO)", nQ, noccA, noccA));
    K->read(psio_, PSIF_DFOCC_INTS);
    EcorrL += 0.5 * G->vector_dot(K);
    G.reset();
    K.reset();

    // DE += 1/2 \sum_{Q} \sum_{i,j} G_ij^Q b_ij^Q
    G = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|oo)", nQ, noccB, noccB));
    G->read(psio_, PSIF_DFOCC_DENS);
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|oo)", nQ, noccB, noccB));
    K->read(psio_, PSIF_DFOCC_INTS);
    EcorrL += 0.5 * G->vector_dot(K);
    G.reset();
    K.reset();
    Eoo = EcorrL - Eoei;

    // DE += \sum_{Q} \sum_{I,A} G_IA^Q b_IA^Q
    G2c_ovA = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|OV)", nQ, noccA * nvirA));
    bQovA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|OV)", nQ, noccA * nvirA));
    G2c_ovA->read(psio_, PSIF_DFOCC_DENS);
    bQovA->read(psio_, PSIF_DFOCC_INTS);
    EcorrL += G2c_ovA->vector_dot(bQovA);
    G2c_ovA.reset();
    bQovA.reset();

    // DE += \sum_{Q} \sum_{i,a} G_ia^Q b_ia^Q
    G2c_ovB = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|ov)", nQ, noccB * nvirB));
    bQovB = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|ov)", nQ, noccB * nvirB));
    G2c_ovB->read(psio_, PSIF_DFOCC_DENS);
    bQovB->read(psio_, PSIF_DFOCC_INTS);
    EcorrL += G2c_ovB->vector_dot(bQovB);
    G2c_ovB.reset();
    bQovB.reset();
    Eov = EcorrL - Eoei - Eoo;

    // DE += 1/2 \sum_{Q} \sum_{A,B} G_AB^Q b_AB^Q
    G = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|VV)", nQ, nvirA, nvirA));
    G->read(psio_, PSIF_DFOCC_DENS, true, true);
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|VV)", nQ, nvirA, nvirA));
    K->read(psio_, PSIF_DFOCC_INTS, true, true);
    EcorrL += 0.5 * G->vector_dot(K);
    G.reset();
    K.reset();

    // DE += 1/2 \sum_{Q} \sum_{a,b} G_ab^Q b_ab^Q
    G = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|vv)", nQ, nvirB, nvirB));
    G->read(psio_, PSIF_DFOCC_DENS, true, true);
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|vv)", nQ, nvirB, nvirB));
    K->read(psio_, PSIF_DFOCC_INTS, true, true);
    EcorrL += 0.5 * G->vector_dot(K);
    G.reset();
    K.reset();
    Evv = EcorrL - Eoei - Eoo - Eov;

}// else if (reference_ == "UNRESTRICTED")

    EccsdL = Eref + EcorrL;
    outfile->Printf("\n\tEnergies re-computed from CC density: \n");
    outfile->Printf("\t------------------------------------- \n");
    outfile->Printf("\tReference Energy (a.u.)            : %20.14f\n", Eref);
    outfile->Printf("\tOne-Electron Energy (a.u.)         : %20.14f\n", Eoei);
    outfile->Printf("\tOO Energy (a.u.)                   : %20.14f\n", Eoo);
    outfile->Printf("\tOV Energy (a.u.)                   : %20.14f\n", Eov);
    outfile->Printf("\tVV Energy (a.u.)                   : %20.14f\n", Evv);
    outfile->Printf("\tLagrangian Corr. Energy (a.u.)     : %20.14f\n", EcorrL);
    outfile->Printf("\tLagrangian Total Energy (a.u.)     : %20.14f\n", EccsdL);

    timer_off("CCL Energy");
} // end ccl_energy

//=======================================================
//          CCL Energy-2
//=======================================================
void DFOCC::ccl_energy2()
{

    timer_on("CCL Energy");

    SharedTensor2d G, K;

    EcorrL = 0.0;
    double Eoei = 0.0;
    double Eoo = 0.0;
    double Eov = 0.0;
    double Evv = 0.0;

if (reference_ == "RESTRICTED") {
    // DE = \sum_{p,q} G_pq h_pq
    EcorrL += G1c->vector_dot(HmoA);
    Eoei = EcorrL;

    // DE += 1/2 \sum_{Q} \sum_{i,j} G_ij^Q b_ij^Q
    // Sep
    G = SharedTensor2d(new Tensor2d("3-Index Separable TPDM (Q|OO)", nQ_ref, noccA, noccA));
    G->read(psio_, PSIF_DFOCC_DENS);
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|OO)", nQ_ref, noccA, noccA));
    K->read(psio_, PSIF_DFOCC_INTS);
    EcorrL += 0.5 * G->vector_dot(K);
    G.reset();
    K.reset();
    // Corr
    G = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|OO)", nQ, noccA, noccA));
    G->read(psio_, PSIF_DFOCC_DENS);
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|OO)", nQ, noccA, noccA));
    K->read(psio_, PSIF_DFOCC_INTS);
    EcorrL += 0.5 * G->vector_dot(K);
    G.reset();
    K.reset();
    Eoo = EcorrL - Eoei;

    // DE += \sum_{Q} \sum_{i,a} G_ia^Q b_ia^Q
    // Sep
    G = SharedTensor2d(new Tensor2d("3-Index Separable TPDM (Q|OV)", nQ_ref, noccA, nvirA));
    G->read(psio_, PSIF_DFOCC_DENS);
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|OV)", nQ_ref, noccA, nvirA));
    K->read(psio_, PSIF_DFOCC_INTS);
    EcorrL += G->vector_dot(K);
    G.reset();
    K.reset();
    // Corr
    G = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|OV)", nQ, noccA, nvirA));
    G->read(psio_, PSIF_DFOCC_DENS);
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|OV)", nQ, noccA, nvirA));
    K->read(psio_, PSIF_DFOCC_INTS);
    EcorrL += G->vector_dot(K);
    G.reset();
    K.reset();
    Eov = EcorrL - Eoei - Eoo;

    // DE += 1/2 \sum_{Q} \sum_{a,b} G_ab^Q b_ab^Q
    // Sep
    G = SharedTensor2d(new Tensor2d("3-Index Separable TPDM (Q|VV)", nQ_ref, nvirA, nvirA));
    G->read(psio_, PSIF_DFOCC_DENS, true, true);
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|VV)", nQ_ref, nvirA, nvirA));
    K->read(psio_, PSIF_DFOCC_INTS, true, true);
    EcorrL += 0.5 * G->vector_dot(K);
    G.reset();
    K.reset();
    // Corr
    G = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|VV)", nQ, nvirA, nvirA));
    G->read(psio_, PSIF_DFOCC_DENS, true, true);
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|VV)", nQ, nvirA, nvirA));
    K->read(psio_, PSIF_DFOCC_INTS, true, true);
    EcorrL += 0.5 * G->vector_dot(K);
    G.reset();
    K.reset();
    Evv = EcorrL - Eoei - Eoo - Eov;

}// end if (reference_ == "RESTRICTED")

else if (reference_ == "UNRESTRICTED") {

}// else if (reference_ == "UNRESTRICTED")

    EccsdL = Eref + EcorrL;
    outfile->Printf("\n\tEnergies re-computed from Fock-adjusted CC density: \n");
    outfile->Printf("\t--------------------------------------------------- \n");
    outfile->Printf("\tOne-Electron Energy (a.u.)         : %20.14f\n", Eoei);
    outfile->Printf("\tOO Energy (a.u.)                   : %20.14f\n", Eoo);
    outfile->Printf("\tOV Energy (a.u.)                   : %20.14f\n", Eov);
    outfile->Printf("\tVV Energy (a.u.)                   : %20.14f\n", Evv);
    outfile->Printf("\tLagrangian Corr. Energy (a.u.)     : %20.14f\n", EcorrL);
    outfile->Printf("\tLagrangian Total Energy (a.u.)     : %20.14f\n", EccsdL);

    timer_off("CCL Energy");
} // end ccl_energy2


}} // End Namespaces
