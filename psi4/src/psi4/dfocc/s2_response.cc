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

void DFOCC::s2_response()
{

    outfile->Printf("\tComputing <S**2>...\n");

    timer_on("s2_response");
    SharedTensor2d T, T2;

    //=========================
    // Read AO basis SO
    //=========================
    SharedTensor2d SmoAB = SharedTensor2d(new Tensor2d("MO-basis Alpha-Beta Overlap Ints", nmo_, nmo_));
    SharedTensor2d SmoBA = SharedTensor2d(new Tensor2d("MO-basis Beta-Alpha Overlap Ints", nmo_, nmo_));
    SharedTensor2d temp = SharedTensor2d(new Tensor2d("Temp", nso_, nmo_));

    // AB
    temp->gemm(false, false, Sso, CmoB, 1.0, 0.0);
    SmoAB->gemm(true, false, CmoA, temp, 1.0, 0.0);

    // BA
    temp->zero();
    temp->gemm(false, false, Sso, CmoA, 1.0, 0.0);
    SmoBA->gemm(true, false, CmoB, temp, 1.0, 0.0);
    temp.reset();

    // Form overlap blocks
    // S_Ij
    SharedTensor2d SooAB = SharedTensor2d(new Tensor2d("S <O|o>", noccA, noccB));
    SooAB->form_ooAB(SmoAB);

    // S_Ia
    SharedTensor2d SiaAB = SharedTensor2d(new Tensor2d("S <I|a>", naoccA, navirB));
    SiaAB->form_act_ov(nfrzc, noccB, SmoAB);

    // S_iA
    SharedTensor2d SiaBA = SharedTensor2d(new Tensor2d("S <i|A>", naoccB, navirA));
    SiaBA->form_act_ov(nfrzc, noccA, SmoBA);

    //=========================
    // Comput <S2>_ref
    //=========================
    // <S2>_ref = 1/4 (Na - Nb)^2 + 1/2 (Na + Nb) - \sum_{Ij} S_Ij * S_jI
    s2_ref = 0.0;
    s2_ref = 0.25 * (noccA - noccB) * (noccA - noccB);
    s2_ref += 0.5 * (noccA + noccB);
    s2_ref -= SooAB->vector_dot(SooAB);
    SooAB.reset();

    //=========================
    // Comput <S2>_resp
    //=========================
    // <S2> = <S2>_ref - 1/2 \sum_{Ij} \sum_{Ab} S_Ib S_jA t_Ij^Ab
    s2_resp = 0.0;
    s2_proj = 0.0;
    s2_resp += s2_ref;
    s2_proj += s2_ref;

    // Compute amplitude contribution
    T = SharedTensor2d(new Tensor2d("T2_1 <Ij|Ab>", naoccA, naoccB, navirA, navirB));
    T2 = SharedTensor2d(new Tensor2d("T2_1 (Ib,jA)", naoccA, navirB, naoccB, navirA));
    if (orb_opt_ == "FALSE" && mp2_amp_type_ == "DIRECT") t2AB_ump2_direct(T);
    else T->read(psio_, PSIF_DFOCC_AMPS);
    T2->sort(1423, T, 1.0, 0.0);
    T.reset();

    // X_Ib = -1/2 \sum_{jA} T2(Ib,jA) S_jA. NoTe: For projected value there is no 1/2 factor
    SharedTensor2d X = SharedTensor2d(new Tensor2d("X <O|v>", noccA, nvirB));
    X->gemv(false, T2, SiaBA, -1.0, 0.0);
    T2.reset();
    SiaBA.reset();

    // <S2> += \sum_{Ib} X_Ib S_Ib
    double value = 0.0;
    value = X->vector_dot(SiaAB);
    X.reset();
    SiaAB.reset();
    s2_proj += value;
    s2_resp += 0.5 * value;

    double value2 = 0.0;
    double value3 = 0.0;
    double s_qn_resp = 0.0;
    double s_qn_proj = 0.0;
    double mult_resp = 0.0;
    double mult_proj = 0.0;

    value2 = 1.0 + (4.0*s2_resp);
    value3 = 1.0 + (4.0*s2_proj);
    s_qn_resp = 0.5 * (sqrt(value2) - 1.0);
    s_qn_proj = 0.5 * (sqrt(value3) - 1.0);
    mult_resp = sqrt(value2);
    mult_proj = sqrt(value3);

    // Print
    outfile->Printf("\t<S**2>_reference (a.u.): %12.8f\n", s2_ref);
    outfile->Printf("\t<S**2>_response (a.u.) : %12.8f\n", s2_resp);
    outfile->Printf("\t<S**2>_projected (a.u.): %12.8f\n", s2_proj);
    outfile->Printf("\tSpin q.n. (resp)       : %12.8f\n", s_qn_resp);
    outfile->Printf("\tSpin q.n. (proj)       : %12.8f\n", s_qn_proj);
    outfile->Printf("\tSpin mult. (resp)      : %12.8f\n", mult_resp);
    outfile->Printf("\tSpin mult. (proj)      : %12.8f\n", mult_proj);


    timer_off("s2_response");
} // s2_response

//=========================
// <S**2> Lagrangian
//=========================
void DFOCC::s2_lagrangian()
{

    outfile->Printf("\tComputing <S**2>_lagrangian...\n");

    timer_on("s2_lagrangian");
    SharedTensor2d T, T2;

    //=========================
    // Read AO basis SO
    //=========================
    SharedTensor2d SmoAB = SharedTensor2d(new Tensor2d("MO-basis Alpha-Beta Overlap Ints", nmo_, nmo_));
    SharedTensor2d SmoBA = SharedTensor2d(new Tensor2d("MO-basis Beta-Alpha Overlap Ints", nmo_, nmo_));
    SharedTensor2d temp = SharedTensor2d(new Tensor2d("Temp", nso_, nmo_));

    // AB
    temp->gemm(false, false, Sso, CmoB, 1.0, 0.0);
    SmoAB->gemm(true, false, CmoA, temp, 1.0, 0.0);

    // BA
    temp->zero();
    temp->gemm(false, false, Sso, CmoA, 1.0, 0.0);
    SmoBA->gemm(true, false, CmoB, temp, 1.0, 0.0);
    temp.reset();

    // Form overlap blocks
    // S_Ij
    SharedTensor2d SooAB = SharedTensor2d(new Tensor2d("S <O|o>", noccA, noccB));
    SooAB->form_ooAB(SmoAB);

    // S_iJ
    SharedTensor2d SooBA = SharedTensor2d(new Tensor2d("S <o|O>", noccB, noccA));
    SooBA->form_ooAB(SmoBA);

    // S_Ia
    SharedTensor2d SiaAB = SharedTensor2d(new Tensor2d("S <I|a>", naoccA, navirB));
    SharedTensor2d SovAB = SharedTensor2d(new Tensor2d("S <O|v>", noccA, nvirB));
    SiaAB->form_act_ov(nfrzc, noccB, SmoAB);
    SovAB->form_ov(noccB, SmoAB);

    // S_iA
    SharedTensor2d SiaBA = SharedTensor2d(new Tensor2d("S <i|A>", naoccB, navirA));
    SharedTensor2d SovBA = SharedTensor2d(new Tensor2d("S <o|V>", noccB, nvirA));
    SiaBA->form_act_ov(nfrzc, noccA, SmoBA);
    SovBA->form_ov(noccA, SmoBA);

    // S_Vo
    SharedTensor2d SvoAB = SharedTensor2d(new Tensor2d("S <V|o>", nvirA, noccB));
    for (int a = 0; a < nvirA; a++) {
         for (int i = 0; i < noccB; i++) {
              SvoAB->set(a, i, SmoAB->get(a + noccA, i));
         }
    }

    // S_vO
    SharedTensor2d SvoBA = SharedTensor2d(new Tensor2d("S <v|O>", nvirB, noccA));
    for (int a = 0; a < nvirB; a++) {
         for (int i = 0; i < noccA; i++) {
              SvoBA->set(a, i, SmoBA->get(a + noccB, i));
         }
    }

    //=========================
    // Comput <S2>_ref
    //=========================
    // <S2>_ref = 1/4 (Na - Nb)^2 + 1/2 (Na + Nb) - \sum_{Ij} S_Ij * S_jI
    s2_ref = 0.0;
    s2_ref = 0.25 * (noccA - noccB) * (noccA - noccB);
    s2_ref += 0.5 * (noccA + noccB);
    s2_ref -= SooAB->vector_dot(SooAB);

    //=========================
    // Comput <S2>_resp
    //=========================
    // <S2> = <S2>_ref - 1/2 \sum_{Ij} \sum_{Ab} S_Ib S_jA t_Ij^Ab
    s2_resp = 0.0;
    s2_proj = 0.0;
    s2_resp += s2_ref;
    s2_proj += s2_ref;
    s2_lag += s2_ref;

    // Compute amplitude contribution
    T = SharedTensor2d(new Tensor2d("T2_1 <Ij|Ab>", naoccA, naoccB, navirA, navirB));
    T2 = SharedTensor2d(new Tensor2d("T2_1 (Ib,jA)", naoccA, navirB, naoccB, navirA));
    if (orb_opt_ == "FALSE" && mp2_amp_type_ == "DIRECT") t2AB_ump2_direct(T);
    else T->read(psio_, PSIF_DFOCC_AMPS);
    T->read(psio_, PSIF_DFOCC_AMPS);
    T2->sort(1423, T, 1.0, 0.0);
    T.reset();

    // X_Ib = -1/2 \sum_{jA} T2(Ib,jA) S_jA. For response
    // X_Ib = -\sum_{jA} T2(Ib,jA) S_jA. For projected
    // X_Ib = -2\sum_{jA} T2(Ib,jA) S_jA. For lagrangian
    SharedTensor2d X = SharedTensor2d(new Tensor2d("X <O|v>", noccA, nvirB));
    X->gemv(false, T2, SiaBA, -1.0, 0.0);
    T2.reset();
    SiaBA.reset();

    // <S2> += \sum_{Ib} X_Ib S_Ib
    double value = 0.0;
    value = X->vector_dot(SiaAB);
    X.reset();
    SiaAB.reset();
    s2_proj += value;
    s2_resp += 0.5 * value;
    //s2_lag += 2.0 * value;
    s2_lag += value;

    //=========================
    // Comput <S2>_lag terms
    //=========================
    // <S**2> -= \sum_{I,m,n} S_In Gc_nm S_mI
    SharedTensor2d SGS;
    SGS = SharedTensor2d(new Tensor2d("SGS <O|O>", noccA, noccA));
    SGS->triple_gemm(SooAB, G1c_ooB, SooBA);
    //s2_lag -= SGS->trace();
    s2_lag -= 0.5 * SGS->trace();

    // <S**2> -= \sum_{I,e,f} S_If Gc_fe S_eI
    SGS->zero();
    SGS->triple_gemm(SovAB, G1c_vvB, SvoBA);
    //s2_lag -= SGS->trace();
    s2_lag -= 0.5 * SGS->trace();
    SGS.reset();

    // <S**2> -= \sum_{i,M,N} S_iN Gc_NM S_Mi
    SGS = SharedTensor2d(new Tensor2d("SGS <o|o>", noccB, noccB));
    SGS->triple_gemm(SooBA, G1c_ooA, SooAB);
    //s2_lag -= SGS->trace();
    s2_lag -= 0.5 * SGS->trace();

    // <S**2> -= \sum_{i,E,F} S_iF Gc_FE S_Ei
    SGS->zero();
    SGS->triple_gemm(SovBA, G1c_vvA, SvoAB);
    //s2_lag -= SGS->trace();
    s2_lag -= 0.5 * SGS->trace();
    SGS.reset();

    // free
    SooAB.reset();
    SooBA.reset();
    SovAB.reset();
    SovBA.reset();
    SvoAB.reset();
    SvoBA.reset();
    //SaiAB.reset();
    //SaiBA.reset();

    //=========================
    // Comput Multiplicities
    //=========================
    double value2 = 0.0;
    double value3 = 0.0;
    double value4 = 0.0;
    double s_qn_resp = 0.0;
    double s_qn_proj = 0.0;
    double s_qn_lag = 0.0;
    double mult_resp = 0.0;
    double mult_proj = 0.0;
    double mult_lag = 0.0;

    value2 = 1.0 + (4.0*s2_resp);
    value3 = 1.0 + (4.0*s2_proj);
    value4 = 1.0 + (4.0*s2_lag);
    s_qn_resp = 0.5 * (sqrt(value2) - 1.0);
    s_qn_proj = 0.5 * (sqrt(value3) - 1.0);
    s_qn_lag = 0.5 * (sqrt(value4) - 1.0);
    mult_resp = sqrt(value2);
    mult_proj = sqrt(value3);
    mult_lag = sqrt(value4);

    // Print
    outfile->Printf("\t<S**2>_reference (a.u.)  : %12.8f\n", s2_ref);
    outfile->Printf("\t<S**2>_response (a.u.)   : %12.8f\n", s2_resp);
    outfile->Printf("\t<S**2>_projected (a.u.)  : %12.8f\n", s2_proj);
    outfile->Printf("\t<S**2>_lagrangian (a.u.) : %12.8f\n", s2_lag);
    outfile->Printf("\tSpin q.n. (resp)         : %12.8f\n", s_qn_resp);
    outfile->Printf("\tSpin q.n. (proj)         : %12.8f\n", s_qn_proj);
    outfile->Printf("\tSpin q.n. (lagr)         : %12.8f\n", s_qn_lag);
    outfile->Printf("\tSpin mult. (resp)        : %12.8f\n", mult_resp);
    outfile->Printf("\tSpin mult. (proj)        : %12.8f\n", mult_proj);
    outfile->Printf("\tSpin mult. (lagr)        : %12.8f\n", mult_lag);


    timer_off("s2_lagrangian");
} // s2_lagrangian

}} // End Namespaces
