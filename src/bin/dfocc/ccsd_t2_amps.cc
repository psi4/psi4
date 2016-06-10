/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2016 The Psi4 Developers.
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

#include <libqt/qt.h>
#include "defines.h"
#include "dfocc.h"

using namespace psi;
using namespace std;


namespace psi{ namespace dfoccwave{
  
void DFOCC::ccsd_t2_amps()
{

    // defs
    SharedTensor2d K, I, T, Tnew, U, Tau, W, X, Y;

    // t_ij^ab <= X(ia,jb) + X(jb,a) = 2Xt(ia,jb)
    // X(ia,jb) = \sum_{e} t_ij^ae Ft_be = \sum_{e} T(ia,je) Ft_be
    X = SharedTensor2d(new Tensor2d("X (IA|JB)", naoccA, navirA, naoccA, navirA));
    X->contract(false, true, naoccA * navirA * naoccA, navirA, navirA, t2, FtabA, 1.0, 0.0);

    // t_ij^ab <= X(ia,jb) + X(jb,a) = 2Xt(ia,jb)
    // X(ia,jb) = -\sum_{m} t_mj^ab Ft_mi = -\sum_{m} Ft(m,i) T(ma,jb)
    X->contract(true, false, naoccA, naoccA * navirA * navirA, naoccA, FtijA, t2, -1.0, 1.0);

    // t_ij^ab <= X(ia,jb) + X(jb,a) = 2Xt(ia,jb)
    // X(ia,jb) = \sum_{Q} t'_ia^Q b_jb^Q
    T = SharedTensor2d(new Tensor2d("T1p (Q|IA)", nQ, naoccA, navirA));
    T->read(psio_, PSIF_DFOCC_AMPS);
    X->gemm(true, false, T, bQiaA, 1.0, 1.0);
    T.reset();

    // t_ij^ab <= X(ia,jb) + X(jb,a) = 2Xt(ia,jb)
    // X(ia,jb) = -\sum_{Q} t_ai^Q t_jb^Q
    U = SharedTensor2d(new Tensor2d("T1 (Q|AI)", nQ, navirA, naoccA));
    U->read(psio_, PSIF_DFOCC_AMPS);
    K = SharedTensor2d(new Tensor2d("Temp (Q|IA)", nQ, naoccA, navirA));
    K->swap_3index_col(U);
    U.reset();
    T = SharedTensor2d(new Tensor2d("T1 (Q|IA)", nQ, naoccA, navirA));
    T->read(psio_, PSIF_DFOCC_AMPS);
    X->gemm(true, false, K, T, -1.0, 1.0);
    T.reset();
    X->symmetrize();

    // t_ij^ab <= <ij|ab>
    Tnew = SharedTensor2d(new Tensor2d("New T2 (IA|JB)", naoccA, navirA, naoccA, navirA));
    Tnew->gemm(true, false, bQiaA, bQiaA, 1.0, 0.0); 

    // Contributions of X
    Tnew->axpy(X, 2.0);
    X.reset();

    // Write and close
    Tnew->write_symm(psio_, PSIF_DFOCC_AMPS);
    Tnew.reset();

    // WmnijT2
    ccsd_WmnijT2();

    // WmbejT2
    ccsd_WmbejT2();

    // WijamT2
    //if (itr_occ > 1) ccsd_WijamT2();

    // WabefT2
    if (Wabef_type_ == "AUTO") {
	if (!do_ppl_hm) ccsd_Wabef2T2();
	else {
	    ccsd_WijamT2_high_mem();
	    ccsd_WabefT2_high_mem();
	}
    }
    else if (Wabef_type_ == "LOW_MEM") ccsd_Wabef2T2();
    else if (Wabef_type_ == "HIGH_MEM") {
	ccsd_WijamT2_high_mem();
	ccsd_WabefT2_high_mem();
    }
    else if (Wabef_type_ == "CD") {
	ccsd_WijamT2();
	ccsd_WabefT2_cd();
    }

    // Denom
    Tnew = SharedTensor2d(new Tensor2d("New T2 (IA|JB)", naoccA, navirA, naoccA, navirA));
    Tnew->read_symm(psio_, PSIF_DFOCC_AMPS);
    Tnew->apply_denom_chem(nfrzc, noccA, FockA);

    // Reset T1
    rms_t1 = t1newA->rms(t1A); 
    SharedTensor2d Rt1A = SharedTensor2d(new Tensor2d("RT1 <I|A>", naoccA, navirA));
    Rt1A->copy(t1newA);
    Rt1A->subtract(t1A);
    t1A->copy(t1newA);

    // Reset T2
    rms_t2 = Tnew->rms(t2); 
    // Error vector
    Tau = SharedTensor2d(new Tensor2d("RT2 (IA|JB)", naoccA, navirA, naoccA, navirA));
    Tau->copy(Tnew);
    Tau->subtract(t2);
    t2->copy(Tnew);
    Tnew.reset();

    // DIIS
    boost::shared_ptr<Matrix> RT2(new Matrix("RT2", naoccA*navirA, naoccA*navirA));
    Tau->to_matrix(RT2);
    Tau.reset();
    boost::shared_ptr<Matrix> T2(new Matrix("T2", naoccA*navirA, naoccA*navirA));
    t2->to_matrix(T2);
    boost::shared_ptr<Matrix> RT1(new Matrix("RT1", naoccA, navirA));
    Rt1A->to_matrix(RT1);
    Rt1A.reset();
    boost::shared_ptr<Matrix> T1(new Matrix("T1", naoccA, navirA));
    t1A->to_matrix(T1);

    // add entry
    if (do_diis_ == 1) ccsdDiisManager->add_entry(4, RT2.get(), RT1.get(), T2.get(), T1.get());
    RT2.reset();
    RT1.reset();

    // extrapolate
    if (do_diis_ == 1) {
        if (ccsdDiisManager->subspace_size() >= cc_mindiis_) ccsdDiisManager->extrapolate(2, T2.get(), T1.get());
        t2->set2(T2);
        t1A->set2(T1);
    }
    T1.reset();
    T2.reset();

    // Form Tau(ia,jb) = T(ia,jb) + t(ia) * t(jb)
    Tau = SharedTensor2d(new Tensor2d("Tau (IA|JB)", naoccA, navirA, naoccA, navirA));
    ccsd_tau_amps(Tau,t2);
 
    // Energy
    U = SharedTensor2d(new Tensor2d("2*Tau(ia,jb) - Tau(ib,ja)", naoccA, navirA, naoccA, navirA));
    U->sort(1432, Tau, 1.0, 0.0);
    U->scale(-1.0);
    U->axpy(Tau, 2.0);
    Tau.reset();
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (IA|JB)", naoccA, navirA, naoccA, navirA));
    K->gemm(true, false, bQiaA, bQiaA, 1.0, 0.0);
    Ecorr = U->vector_dot(K);
    U.reset();
    K.reset();
    Eccsd = Escf + Ecorr;

}// end ccsd_t2_amps

//======================================================================
//    CCSD: u_ij^ab = 2*t_ij^ab - t_ji^ab; T(ia,jb)
//======================================================================             
void DFOCC::ccsd_u2_amps(SharedTensor2d &U, SharedTensor2d &T)
{
    U->sort(1432, T, 1.0, 0.0);
    U->scale(-1.0);
    U->axpy(T, 2.0);
}// end ccsd_u2_amps

//======================================================================
//    CCSD: u_ij^ab = 2*t_ij^ab - t_ji^ab; T(ij,ab)
//======================================================================             
void DFOCC::ccsd_u2_amps2(SharedTensor2d &U, SharedTensor2d &T)
{
    U->sort(2134, T, 1.0, 0.0);
    U->scale(-1.0);
    U->axpy(T, 2.0);
}// end ccsd_u2_amps2

//======================================================================
//    CCSD: T'(ib,ja) = t_ij^ab = T(ia,jb)
//======================================================================             
void DFOCC::ccsd_t2_prime_amps(SharedTensor2d &U, SharedTensor2d &T)
{
    U->sort(1432, T, 1.0, 0.0);
}// end ccsd_t2_prime_amps

//======================================================================
//    CCSD: T'(ib,ja) = t_ij^ab = T(ij,ab) 
//======================================================================             
void DFOCC::ccsd_t2_prime_amps2(SharedTensor2d &U, SharedTensor2d &T)
{
    U->sort(1423, T, 1.0, 0.0);
}// end ccsd_t2_prime_amps2

//======================================================================
//    CCSD: Tau_ij^ab = t_ij^ab + t_i^a t_j^b
//======================================================================             
void DFOCC::ccsd_tau_amps(SharedTensor2d &U, SharedTensor2d &T)
{
    U->dirprd224(t1A, t1A);
    U->add(T);
}// end ccsd_tau_amps

//======================================================================
//    CCSD: \tilde{Tau}_ij^ab = t_ij^ab + 1/2 t_i^a t_j^b
//======================================================================             
void DFOCC::ccsd_tau_tilde_amps(SharedTensor2d &U, SharedTensor2d &T)
{
    U->dirprd224(t1A, t1A, 0.5, 0.0);
    U->add(T);
}// end ccsd_tau_tilde_amps

}} // End Namespaces

