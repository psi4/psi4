/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2021 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This file is part of Psi4.
 *
 * Psi4 is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * Psi4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with Psi4; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

/** Standard library includes */

#include "dfocc.h"
#include "defines.h"

#include "psi4/psifiles.h"
#include "psi4/libiwl/iwl.hpp"
#include "psi4/libpsio/psio.h"
#include "psi4/libciomr/libciomr.h"
#include "psi4/libqt/qt.h"
#include "psi4/libpsi4util/process.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/vector.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/wavefunction.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/mintshelper.h"
//#include "psi4/lib3index/df_helper.h"
//#include "psi4/lib3index/dftensor.h"
#include "psi4/libmints/twobody.h"
#include "psi4/libmints/integral.h"
#include "psi4/libmints/sieve.h"
#include "psi4/libfock/jk.h"
#include "psi4/psi4-dec.h"

#include <limits>
#include <fstream>

#ifdef _OPENMP
#include <omp.h>
#include "psi4/libpsi4util/process.h"
#endif

namespace psi {
namespace dfoccwave {

//======================================================================
//     CIS: MO-BASIS
//======================================================================
void DFOCC::cis() {

    timer_on("CIS");
    SharedTensor2d J, K, L, M, T, U, Tau;
    SharedTensor1d EciA;

    // malloc Effective hamiltonian
    //outfile->Printf("\tI am here\n");
    nconfigA = naoccA * navirA;
    nconfigB = naoccB * navirB;

  //====================================//
  // FULL DIAG
  //====================================//
  if (diag_method == "FULL_DIAG") {
      if (reference_ == "RESTRICTED") {

          // Malloc
          auto HciA = std::make_shared<Tensor2d>("H <IA|IA>", naoccA, navirA, naoccA, navirA);
          CciA = std::make_shared<Tensor2d>("C CIS", nconfigA, nconfigA);

          // Form Integrals
          // (IJ|AB)
          K = std::make_shared<Tensor2d>("DF_BASIS_CC MO Ints (IJ|AB)", naoccA, naoccA, navirA, navirA);
          L = std::make_shared<Tensor2d>("DF_BASIS_CC B (Q|IJ)", nQ, naoccA, naoccA);
          M = std::make_shared<Tensor2d>("DF_BASIS_CC B (Q|AB)", nQ, navirA, navirA);
          L->read(psio_, PSIF_DFOCC_INTS);
          M->read(psio_, PSIF_DFOCC_INTS, true, true);
          K->gemm(true, false, L, M, 1.0, 0.0);
          M.reset();
          L.reset();
          HciA->sort(1324, K, -1.0, 0.0);
          K.reset();

          // (IA|JB)
          J = std::make_shared<Tensor2d>("DF_BASIS_CC MO Ints (IA|JB)", naoccA, navirA, naoccA, navirA);
          K = std::make_shared<Tensor2d>("DF_BASIS_CC B (Q|IA)", nQ, naoccA, navirA);
          K->read(psio_, PSIF_DFOCC_INTS);
          J->gemm(true, false, K, K, 1.0, 0.0);
          K.reset();
          HciA->axpy(J,2.0);
          J.reset();

          // Form <Phi_mu | H | Phi_nu>, 
          outfile->Printf("\n\tComputing DF-CIS energy...\n");
          // Diagonal
          for(int i = 0; i < naoccA; ++i) {
              for(int a = 0; a < navirA; ++a) {
                  int ia = (i*navirA) + a;
                  //double value = Escf + FockA->get(a+noccA,a+noccA) - FockA->get(i+nfrzc,i+nfrzc);
                  double value = FockA->get(a+noccA,a+noccA) - FockA->get(i+nfrzc,i+nfrzc);
                  HciA->add(ia,ia,value);
              }
          }

          // diagonalize the effective Hamiltonian
          timer_on("HciA diagonalize");
          EciA = std::make_shared<Tensor1d>("CIS Alpha E vector", nconfigA);
          //auto EciB = std::make_shared<Tensor1d>("CIS Alpha E vector", nconfigA);

          // Diagonalize
          //HciA->print();
          HciA->diagonalize(CciA, EciA, cutoff);
          //EciA->print();
          Ecis = EciA->get(0);
          Ecorr = Ecis - Escf;

          // print excitation energies
          outfile->Printf("\n\tExcitation energies (in au)...\n");
          for(int i = 0; i < nroot; ++i) {
              //outfile->Printf("\tExcited-state: %2d, Omega(%2d): %12.10f\n", i+1, i+1, EciA->get(i)-Escf);
              outfile->Printf("\tExcited-state: %2d, Omega(%2d): %12.10f\n", i+1, i+1, EciA->get(i));
          }

          // print excitation energies
          outfile->Printf("\n\tExcitation energies (in eV)...\n");
          for(int i = 0; i < nroot; ++i) {
              //outfile->Printf("\tExcited-state: %2d, Omega(%2d): %12.6f\n", i+1, i+1, (EciA->get(i)-Escf)*hartree2ev);
              outfile->Printf("\tExcited-state: %2d, Omega(%2d): %12.6f\n", i+1, i+1, EciA->get(i)*hartree2ev);
          }

          timer_off("HciA diagonalize");

          // Reset
          HciA.reset();
      }//RHF
      else if (reference_ == "UNRESTRICTED") {
           outfile->Printf("\n\tFull diagonalization is not valid for UHF!\n");
           exit(1);
      }//UHF
  }//if (diag_method == "FULL_DIAG")

  //====================================//
  // DAVIDSON
  //====================================//
  else if (diag_method == "DAVIDSON") {
    outfile->Printf("\n\tComputing DF-CIS energy...\n");
    if (reference_ == "RESTRICTED") {
        CciA = std::make_shared<Tensor2d>("C <I|A>", nconfigA, nroot);
        EciA = std::make_shared<Tensor1d>("CIS E vector", nroot);

        // cis-sigma
        cis_davidson(nroot, EciA, CciA, tol_davidson);
        Ecis = EciA->get(0);
        Ecorr = Ecis - Escf;

        // print excitation energies
        outfile->Printf("\n\tExcitation energies (in au)...\n");
        for(int i = 0; i < nroot; ++i) {
            outfile->Printf("\tExcited-state: %2d, Omega(%2d): %12.10f\n", i+1, i+1, EciA->get(i));
        }

        // print excitation energies
        outfile->Printf("\n\tExcitation energies (in eV)...\n");
        for(int i = 0; i < nroot; ++i) {
            outfile->Printf("\tExcited-state: %2d, Omega(%2d): %12.6f\n", i+1, i+1, EciA->get(i)*hartree2ev);
        }
    }//RHF
    else if (reference_ == "UNRESTRICTED") {
    }//UHF

  }//else if (diag_method == "DAVIDSON")

  //====================================//
  // Oscillator Strength
  //====================================//
    if (reference_ == "RESTRICTED") {
        // Oscillator Strength
        auto fosA = std::make_shared<Tensor1d>("f <0|R>", nroot);
        cis_oscillator_strength(CciA, EciA, fosA);

        // Print all CI coefficients
        if (print_ci_vecs == "TRUE") CciA->print();

        // Print Leading CI coefficients
        outfile->Printf("\n\tLeading CI coefficients...\n");
        for(int j = 0; j < nroot; ++j) {
            outfile->Printf("\n\tRoot(%2d):  \n", j+1);
            outfile->Printf("\tExcitation Energy: %12.10f\n", EciA->get(j));
            outfile->Printf("\tOscillator Strength: %12.8f\n", fosA->get(j));
            for(int i = 0; i < naoccA; ++i) {
                for(int a = 0; a < navirA; ++a) {
                    int ia = (i*navirA) + a;
                    if (std::fabs(CciA->get(ia,j)) >= 1e-2) 
                        outfile->Printf("\tExcitation: %2d -> %2d  %12.6f\n", i+nfrzc+1, a+noccA+1, CciA->get(ia,j));
                }
            }
        }

        // Reset
        CciA.reset();
        EciA.reset();
    }//RHF
    else if (reference_ == "UNRESTRICTED") {
    }//UHF

    timer_off("CIS");

}  // end of cis

//======================================================================
//     Diagonal CIS H
//======================================================================
void DFOCC::cis_diagonal(SharedTensor1d H) {
    SharedTensor2d J, K, L, M, T, U, Tau;

    // H_ia = e_a - e_i
    for(int i = 0; i < naoccA; ++i) {
        for(int a = 0; a < navirA; ++a) {
	    int ia = (i*navirA) + a;
            double value = FockA->get(a+noccA,a+noccA) - FockA->get(i+nfrzc,i+nfrzc);
            H->set(ia,value);
	}
    }

    // Form Integrals
    // -(II|AA)
    K = std::make_shared<Tensor2d>("DF_BASIS_CC B (Q|IJ)", nQ, naoccA, naoccA);
    L = std::make_shared<Tensor2d>("DF_BASIS_CC B (Q|AB)", nQ, navirA, navirA);
    K->read(psio_, PSIF_DFOCC_INTS);
    L->read(psio_, PSIF_DFOCC_INTS, true, true);
    for(int i = 0; i < naoccA; ++i) {
	int ii = (i*naoccA) + i;
        for(int a = 0; a < navirA; ++a) {
	    int aa = (a*navirA) + a;
	    int ia = (i*navirA) + a;
	    double value = 0.0;
            for(int Q = 0; Q < nQ; ++Q) {
		value -= K->get(Q,ii) * L->get(Q,aa);
	    }
            H->add(ia,value);
	}
    }
    L.reset();
    K.reset();

    // 2*(IA|IA)
    K = std::make_shared<Tensor2d>("DF_BASIS_CC B (Q|IA)", nQ, naoccA, navirA);
    K->read(psio_, PSIF_DFOCC_INTS);
    for(int i = 0; i < naoccA; ++i) {
        for(int a = 0; a < navirA; ++a) {
	    int ia = (i*navirA) + a;
	    double value = 0.0;
            for(int Q = 0; Q < nQ; ++Q) {
		value += K->get(Q,ia) * K->get(Q,ia);
	    }
            H->add(ia,2.0*value);
	}
    }
    K.reset();
    //H->print();

}  // end of cis_diagonal

//======================================================================
//     APProximately diagonal CIS H
//======================================================================
void DFOCC::cis_diagonal_approx(SharedTensor1d H) {
    SharedTensor2d J, K, L, M, T, U, Tau;

    // H_ia = e_a - e_i
    for(int i = 0; i < naoccA; ++i) {
        for(int a = 0; a < navirA; ++a) {
	    int ia = (i*navirA) + a;
            double value = FockA->get(a+noccA,a+noccA) - FockA->get(i+nfrzc,i+nfrzc);
            H->set(ia,value);
	}
    }
    //H->print();

}  // end of cis_diagonal_approx

//======================================================================
//     CIS SIGMA: MO Basis
//======================================================================
void DFOCC::cis_sigma(SharedTensor1d C, SharedTensor1d sigma) {

    timer_on("CIS-Sigma");
    SharedTensor2d J, K, L, M, T, U, Tau;

    // Malloc
    auto CQ = std::make_shared<Tensor1d>("C <Q>", nQ);
    auto CiaQ = std::make_shared<Tensor2d>("C <Q|IA>", nQ, naoccA, navirA);
    auto sigmaA = std::make_shared<Tensor2d>("Sigma <I|A>", naoccA, navirA);
    auto CiaA = std::make_shared<Tensor2d>("C <I|A>", naoccA, navirA);

    // Expand C 
    for(int i = 0; i < naoccA; ++i) {
        for(int a = 0; a < navirA; ++a) {
	    int ia = (i*navirA) + a;
            CiaA->set(i, a, C->get(ia));
	}
    }

    // C(Q) = 2\sum_{jb} b_jb^Q c_j^b
    K = std::make_shared<Tensor2d>("DF_BASIS_CC B (Q|IA)", nQ, naoccA, navirA);
    K->read(psio_, PSIF_DFOCC_INTS);
    CQ->gemv(false, nQ, naoccA * navirA, K, CiaA, 2.0, 0.0);

    // C_ia^Q = \sum_{b} b_ab^Q c_i^b
    M = std::make_shared<Tensor2d>("DF_BASIS_CC B (Q|AB)", nQ, navirA, navirA);
    M->read(psio_, PSIF_DFOCC_INTS, true, true);
    auto CaiQ = std::make_shared<Tensor2d>("C <Q|AI>", nQ, navirA, naoccA);
    CaiQ->contract(false, true, nQ * navirA, naoccA, navirA, M, CiaA, 1.0, 0.0);
    M.reset();
    CiaQ->swap_3index_col(CaiQ);
    CaiQ.reset();

    // S(ia) = \sum(Q) C_Q b_ia^Q 
    sigmaA->gemv(true, K, CQ, 1.0, 0.0);
    K.reset();
    CQ.reset();

    // S(ia) -= \sum(Qj) b_ij^Q C_ja^Q
    K = std::make_shared<Tensor2d>("DF_BASIS_CC B (Q|IJ)", nQ, naoccA, naoccA);
    K->read(psio_, PSIF_DFOCC_INTS);
    sigmaA->contract(true, false, naoccA, navirA, nQ*naoccA, K, CiaQ, -1.0, 1.0);
    K.reset();

    // Diagonal
    for(int i = 0; i < naoccA; ++i) {
        for(int a = 0; a < navirA; ++a) {
            double value = FockA->get(a+noccA,a+noccA) - FockA->get(i+nfrzc,i+nfrzc);
            sigmaA->add(i, a, value*CiaA->get(i,a));
	}
    }
    //CiaA->print();
    //sigmaA->print();

    // Compress sigma
    for(int i = 0; i < naoccA; ++i) {
        for(int a = 0; a < navirA; ++a) {
	    int ia = (i*navirA) + a;
            sigma->set(ia, sigmaA->get(i,a));
	}
    }
    sigmaA.reset();
    CiaA.reset();

    timer_off("CIS-Sigma");

}  // end of cis_sigma

//======================================================================
//     CIS Oscillator strength
//======================================================================
void DFOCC::cis_oscillator_strength(SharedTensor2d Cci, SharedTensor1d Eci, SharedTensor1d fos) {
    // Build Hso
    std::vector<SharedMatrix> Xtemp;
    SharedTensor2d Xso_, Yso_, Zso_, Rso_, Ria_;
    Xso_ = std::make_shared<Tensor2d>("Dipole X <mu|nu>", nso_, nso_);
    Yso_ = std::make_shared<Tensor2d>("Dipole Y <mu|nu>", nso_, nso_);
    Zso_ = std::make_shared<Tensor2d>("Dipole Z <mu|nu>", nso_, nso_);
    Rso_ = std::make_shared<Tensor2d>("Dipole R <mu|nu>", nso_, nso_);
    Ria_ = std::make_shared<Tensor2d>("Dipole R <i|a>", naoccA, navirA);

    // Read SO-basis one-electron integrals
    MintsHelper mints(MintsHelper(reference_wavefunction_->basisset(), options_, 0));
    Xtemp = mints.ao_dipole();
    for(int mu = 0; mu < nso_; ++mu) {
        for(int nu = 0; nu < nso_; ++nu) {
            Xso_->set(mu,nu, Xtemp[0]->get(mu,nu));
            Yso_->set(mu,nu, Xtemp[1]->get(mu,nu));
            Zso_->set(mu,nu, Xtemp[2]->get(mu,nu));
        }
    }
    Rso_->copy(Xso_);
    Rso_->add(Yso_);
    Rso_->add(Zso_);
    //Rso_->print();

    // Form R_ia 
    auto Rtemp = std::make_shared<Tensor2d>("Dipole R <mu|a>", nso_, navirA);
    Rtemp->gemm(false, false, Rso_, CavirA, 1.0, 0.0);
    Ria_->gemm(true, false, CaoccA, Rtemp, 1.0, 0.0);
    //Ria_->print();

    // Compute transition dipole moments
    outfile->Printf("\n\tComputing transition dipole moments...\n");
    auto MuA = std::make_shared<Tensor1d>("Mu <0|R>", nroot);
    for(int j = 0; j < nroot; ++j) {
        double sum = 0.0;
        for(int i = 0; i < naoccA; ++i) {
            for(int a = 0; a < navirA; ++a) {
	        int ia = (i*navirA) + a;
                sum -= std::sqrt(2.0) * Cci->get(ia,j) * Ria_->get(i,a);
	    }
        }
        MuA->set(j,sum);
    }
    //MuA->print();

    // Compute OS
    outfile->Printf("\n\tComputing oscillator strengths...\n");
    for(int j = 0; j < nroot; ++j) {
        double value = (2.0/3.0) * Eci->get(j) * MuA->get(j) * MuA->get(j); 
        fos->set(j,value);
    }
    //fos->print();

}  // end of cis_os

//======================================================================
//     CIS: AO-BASIS
//======================================================================
void DFOCC::cis_ao() {

    timer_on("CIS");
    SharedTensor2d J, K, L, M, T, U, Tau;
    SharedTensor1d EciA;

    // malloc Effective hamiltonian
    //outfile->Printf("\tI am here\n");
    nconfigA = naoccA * navirA;
    nconfigB = naoccB * navirB;

  if (diag_method == "FULL_DIAG") {
     outfile->Printf("\tWarning: AO-basis CIS algorithm uses only davidson method!\n");
  }

  //====================================//
  // DAVIDSON
  //====================================//
    outfile->Printf("\n\tComputing DF-CIS energy...\n");
    CciA = std::make_shared<Tensor2d>("C <IA|root>", nconfigA, nroot);
    EciA = std::make_shared<Tensor1d>("CIS E vector", nroot);

    // cis-sigma
    cis_davidson(nroot, EciA, CciA, tol_davidson);
    Ecis = EciA->get(0);
    Ecorr = Ecis - Escf;

    // print excitation energies
    outfile->Printf("\n\tExcitation energies (in au)...\n");
    for(int i = 0; i < nroot; ++i) {
	outfile->Printf("\tExcited-state: %2d, Omega(%2d): %12.10f\n", i+1, i+1, EciA->get(i));
    }

    // print excitation energies
    outfile->Printf("\n\tExcitation energies (in eV)...\n");
    for(int i = 0; i < nroot; ++i) {
	outfile->Printf("\tExcited-state: %2d, Omega(%2d): %12.6f\n", i+1, i+1, EciA->get(i)*hartree2ev);
    }

    // Oscillator Strength
    auto fosA = std::make_shared<Tensor1d>("f <0|R>", nroot);
    cis_oscillator_strength(CciA, EciA, fosA);

    // Print all CI coefficients
    if (print_ci_vecs == "TRUE") CciA->print();
    
    // Print Leading CI coefficients
    outfile->Printf("\n\tLeading CI coefficients...\n");
    for(int j = 0; j < nroot; ++j) {
	outfile->Printf("\n\tRoot(%2d):  \n", j+1);
	outfile->Printf("\tExcitation Energy: %12.10f\n", EciA->get(j));
	outfile->Printf("\tOscillator Strength: %12.8f\n", fosA->get(j));
        for(int i = 0; i < naoccA; ++i) {
            for(int a = 0; a < navirA; ++a) {
	        int ia = (i*navirA) + a;
            if (std::fabs(CciA->get(ia,j)) >= 1e-2) 
	        outfile->Printf("\tExcitation: %2d -> %2d  %12.6f\n", i+nfrzc+1, a+noccA+1, CciA->get(ia,j));
	    }
        }
    }

    // Reset
    CciA.reset();
    EciA.reset();

    timer_off("CIS");

}  // end of cis_ao

//======================================================================
//     CIS SIGMA: AO Basis
//======================================================================
void DFOCC::cis_sigma_ao(SharedTensor1d C, SharedTensor1d sigma) {

    timer_on("CIS-Sigma");
    SharedTensor2d J, K, L, M, T, U, Tau;

    // Malloc
    auto sigmaA = std::make_shared<Tensor2d>("Sigma <I|A>", naoccA, navirA);
    auto CiaA = std::make_shared<Tensor2d>("C <I|A>", naoccA, navirA);

    // Expand C 
    for(int i = 0; i < naoccA; ++i) {
        for(int a = 0; a < navirA; ++a) {
	    int ia = (i*navirA) + a;
            CiaA->set(i, a, C->get(ia));
	}
    }

    // Form Ctilde
    auto CtA = std::make_shared<Tensor2d>("C <n|I>", nso_, naoccA);
    CtA->gemm(false, true, CavirA, CiaA, 1.0, 0.0);

  //====================================//
  // JK-terms
  //====================================//
    // init JK
    timer_on("JKobject");
    //static std::shared_ptr<JK> jk_ = JK::build_JK(reference_wavefunction_->basisset(), reference_wavefunction_->get_basisset("DF_BASIS_CC"), options_);
    static std::shared_ptr<JK> jk_= JK::build_JK(reference_wavefunction_->basisset(), reference_wavefunction_->get_basisset("DF_BASIS_CC"), options_, jk_type);

    // 8 GB Memory, 1 G doubles
    //jk_->set_memory(1000000000L);
    jk_->set_memory(jk_memory);

    // Cutoff of 1.0E-10
    jk_->set_cutoff(1.0E-10);

    // Do J/K, Not wK (superfluous)
    jk_->set_do_J(true);
    jk_->set_do_K(true);
    jk_->set_do_wK(false);

    // init
    jk_->initialize();

    // left C = Ctilde
    auto Cmat = std::make_shared<Matrix>("Cleft", nso_, naoccA);
    CtA->to_shared_matrix(Cmat);
    CtA.reset();
    std::vector<SharedMatrix>& C_left = jk_->C_left();
    C_left.clear();
    C_left.push_back(Cmat);
    Cmat.reset();

    // right C = CaoccA
    auto Cmat2 = std::make_shared<Matrix>("Cright", nso_, naoccA);
    CaoccA->to_shared_matrix(Cmat2);
    std::vector<SharedMatrix>& C_right = jk_->C_right();
    C_right.clear();
    C_right.push_back(Cmat2);
    Cmat2.reset();

    // Run the JK object
    jk_->compute();
    SharedMatrix Jmat = jk_->J()[0];
    SharedMatrix Kmat = jk_->K()[0];

    // finalize
    jk_->finalize();
    timer_off("JKobject");

    // set J
    auto JmnA = std::make_shared<Tensor2d>("J <m|n>", nso_, nso_);
    JmnA->set(Jmat);
    Jmat.reset();
    //JmnA->print();

    // set K
    auto KmnA = std::make_shared<Tensor2d>("K <m|n>", nso_, nso_);
    KmnA->set(Kmat);
    Kmat.reset();
    //KmnA->print();

    // Compute Fmn
    // J
    auto FmnA = std::make_shared<Tensor2d>("F <m|n>", nso_, nso_);
    FmnA->axpy(JmnA, 2.0);
    JmnA.reset();
    // K 
    FmnA->axpy(KmnA, -1.0);
    KmnA.reset();
    //FmnA->print();

    // Compute Fai part
    auto X = std::make_shared<Tensor2d>("X <m|I>", nso_, naoccA);
    X->gemm(false, false, FmnA, CaoccA, 1.0, 0.0);
    FmnA.reset();
    auto FaiA = std::make_shared<Tensor2d>("F <A|I>", navirA, naoccA);
    FaiA->gemm(true, false, CavirA, X, 1.0, 0.0);
    X.reset();
    sigmaA->trans(FaiA);
    FaiA.reset();

    // Diagonal
    for(int i = 0; i < naoccA; ++i) {
        for(int a = 0; a < navirA; ++a) {
            double value = FockA->get(a+noccA,a+noccA) - FockA->get(i+nfrzc,i+nfrzc);
            sigmaA->add(i, a, value*CiaA->get(i,a));
	}
    }
    //CiaA->print();
    //sigmaA->print();

    // Compress sigma
    for(int i = 0; i < naoccA; ++i) {
        for(int a = 0; a < navirA; ++a) {
	    int ia = (i*navirA) + a;
            sigma->set(ia, sigmaA->get(i,a));
	}
    }
    sigmaA.reset();
    CiaA.reset();

    timer_off("CIS-Sigma");

}  // end of cis_sigma_ao

//======================================================================
//     CIS SIGMA: AO Basis: MY J-version
//======================================================================
void DFOCC::cis_sigma_ao_v0(SharedTensor1d C, SharedTensor1d sigma) {

    timer_on("CIS-Sigma");
    SharedTensor2d J, K, L, M, T, U, Tau;

    // Malloc
    auto JpQ = std::make_shared<Tensor1d>("Jp <Q>", nQ);
    auto JdpQ = std::make_shared<Tensor1d>("Jdp <Q>", nQ);
    auto CiaQ = std::make_shared<Tensor2d>("C <Q|IA>", nQ, naoccA, navirA);
    auto sigmaA = std::make_shared<Tensor2d>("Sigma <I|A>", naoccA, navirA);
    auto CiaA = std::make_shared<Tensor2d>("C <I|A>", naoccA, navirA);

    // Expand C 
    for(int i = 0; i < naoccA; ++i) {
        for(int a = 0; a < navirA; ++a) {
	    int ia = (i*navirA) + a;
            CiaA->set(i, a, C->get(ia));
	}
    }

    // Form Ctilde
    auto CtA = std::make_shared<Tensor2d>("C <n|I>", nso_, naoccA);
    CtA->gemm(false, true, CavirA, CiaA, 1.0, 0.0);

    // Form Ctilde
    auto PaoA = std::make_shared<Tensor2d>("P <m|n>", nso_, nso_);
    PaoA->gemm(false, true, CaoccA, CtA, 1.0, 0.0);

  //====================================//
  // My J-term: WRONG needs debugging
  //====================================//
    // first half trans
    cis_ao_J_fht(PaoA, JpQ);

    // Form J''(P)
    auto Jm = std::make_shared<Tensor2d>("DF_BASIS_CC Jm <P|Q>", nQ, nQ);
    Jm->read(psio_, PSIF_DFOCC_INTS);
    JdpQ->gemv(false, Jm, JpQ, 1.0, 0.0);
    Jm.reset();
    JpQ.reset();

    // Compute J: second-half trans
    auto JmnA = std::make_shared<Tensor2d>("J <m|n>", nso_, nso_);
    cis_ao_J_sht(JdpQ, JmnA);
    JdpQ.reset();
    PaoA.reset();

  //====================================//
  // JK-terms
  //====================================//
    // init JK
    timer_on("JKobject");
    //static std::shared_ptr<JK> jk_ = JK::build_JK(reference_wavefunction_->basisset(), reference_wavefunction_->get_basisset("DF_BASIS_CC"), options_);
    static std::shared_ptr<JK> jk_= JK::build_JK(reference_wavefunction_->basisset(), reference_wavefunction_->get_basisset("DF_BASIS_CC"), options_, jk_type);

    // 8 GB Memory, 1 G doubles
    //jk_->set_memory(1000000000L);
    jk_->set_memory(jk_memory);

    // Cutoff of 1.0E-10
    jk_->set_cutoff(1.0E-10);

    // Do J/K, Not wK (superfluous)
    jk_->set_do_J(false);
    jk_->set_do_K(true);
    jk_->set_do_wK(false);

    // init
    jk_->initialize();

    // left C = Ctilde
    auto Cmat = std::make_shared<Matrix>("Cleft", nso_, naoccA);
    CtA->to_shared_matrix(Cmat);
    CtA.reset();
    std::vector<SharedMatrix>& C_left = jk_->C_left();
    C_left.clear();
    C_left.push_back(Cmat);
    Cmat.reset();

    // right C = CaoccA
    auto Cmat2 = std::make_shared<Matrix>("Cright", nso_, naoccA);
    CaoccA->to_shared_matrix(Cmat2);
    std::vector<SharedMatrix>& C_right = jk_->C_right();
    C_right.clear();
    C_right.push_back(Cmat2);
    Cmat2.reset();

    // Run the JK object
    jk_->compute();
    SharedMatrix Kmat = jk_->K()[0];

    // finalize
    jk_->finalize();
    timer_off("JKobject");

    // set K
    auto KmnA = std::make_shared<Tensor2d>("K <m|n>", nso_, nso_);
    KmnA->set(Kmat);
    Kmat.reset();
    //KmnA->print();

    // Compute Fmn
    // J
    auto FmnA = std::make_shared<Tensor2d>("F <m|n>", nso_, nso_);
    FmnA->axpy(JmnA, 2.0);
    JmnA.reset();
    // K 
    FmnA->axpy(KmnA, -1.0);
    KmnA.reset();
    //FmnA->print();

    // Compute Fai part
    auto X = std::make_shared<Tensor2d>("X <m|I>", nso_, naoccA);
    X->gemm(false, false, FmnA, CaoccA, 1.0, 0.0);
    FmnA.reset();
    auto FaiA = std::make_shared<Tensor2d>("F <A|I>", navirA, naoccA);
    FaiA->gemm(true, false, CavirA, X, 1.0, 0.0);
    X.reset();
    sigmaA->trans(FaiA);
    FaiA.reset();

    // Diagonal
    for(int i = 0; i < naoccA; ++i) {
        for(int a = 0; a < navirA; ++a) {
            double value = FockA->get(a+noccA,a+noccA) - FockA->get(i+nfrzc,i+nfrzc);
            sigmaA->add(i, a, value*CiaA->get(i,a));
	}
    }
    //CiaA->print();
    //sigmaA->print();

    // Compress sigma
    for(int i = 0; i < naoccA; ++i) {
        for(int a = 0; a < navirA; ++a) {
	    int ia = (i*navirA) + a;
            sigma->set(ia, sigmaA->get(i,a));
	}
    }
    sigmaA.reset();
    CiaA.reset();

    timer_off("CIS-Sigma");

}  // end of cis_sigma_ao_v0

//======================================================================
//     CIS AO Basis J: First-half trans
//======================================================================
void DFOCC::cis_ao_J_fht(SharedTensor2d Pao, SharedTensor1d JpQ) {

    // Read in the basis set informations
    std::shared_ptr<BasisSet> auxiliary_ = get_basisset("DF_BASIS_CC");
    std::shared_ptr<BasisSet> primary_ = get_basisset("ORBITAL");

    std::shared_ptr<BasisSet> zero(BasisSet::zero_ao_basis_set());

    int nthreads = 1;
#ifdef _OPENMP
    nthreads = Process::environment.get_n_threads();
#endif

    std::shared_ptr<ERISieve> sieve_ = std::shared_ptr<ERISieve>(new ERISieve(primary_, cutoff));
    const std::vector<std::pair<int, int> >& shell_pairs = sieve_->shell_pairs();
    int npairs = shell_pairs.size();

    // => Memory Constraints <= //
    int max_rows;
    max_rows = auxiliary_->nshell();

    // => Block Sizing <= //
    std::vector<int> Pstarts;
    int counter = 0;
    Pstarts.push_back(0);
    for (int P = 0; P < auxiliary_->nshell(); P++) {
        int nP = auxiliary_->shell(P).nfunction();
        if (counter + nP > max_rows) {
            counter = 0;
            Pstarts.push_back(P);
        }
        counter += nP;
    }
    Pstarts.push_back(auxiliary_->nshell());

    // => Integrals <= //
    std::shared_ptr<IntegralFactory> rifactory2(new IntegralFactory(auxiliary_, zero, primary_, primary_));
    std::vector<std::shared_ptr<TwoBodyAOInt> > eri;
    std::vector<const double*> buffer;
    for (int t = 0; t < nthreads; t++) {
        eri.push_back(std::shared_ptr<TwoBodyAOInt>(rifactory2->eri()));
        buffer.push_back(eri[t]->buffer());
    }

    // => Master Loop <= //

    for (int block = 0; block < Pstarts.size() - 1; block++) {
        // > Sizing < //

        int Pstart = Pstarts[block];
        int Pstop = Pstarts[block + 1];
        int NP = Pstop - Pstart;

        int pstart = auxiliary_->shell(Pstart).function_index();
        int pstop = (Pstop == auxiliary_->nshell() ? nQ : auxiliary_->shell(Pstop).function_index());
        int np = pstop - pstart;

// > Integrals < //
#pragma omp parallel for schedule(dynamic) num_threads(nthreads)
        for (long int PMN = 0L; PMN < NP * npairs; PMN++) {
            int thread = 0;
#ifdef _OPENMP
            thread = omp_get_thread_num();
#endif

            int P = PMN / npairs + Pstart;
            int MN = PMN % npairs;
            int M = shell_pairs[MN].first;
            int N = shell_pairs[MN].second;

            eri[thread]->compute_shell(P, 0, M, N);

            int nP = auxiliary_->shell(P).nfunction();
            int oP = auxiliary_->shell(P).function_index();

            int nM = primary_->shell(M).nfunction();
            int oM = primary_->shell(M).function_index();

            int nN = primary_->shell(N).nfunction();
            int oN = primary_->shell(N).function_index();

            int index = 0;
            for (int p = 0; p < nP; p++) {
                double sum = 0.0;
                for (int m = 0; m < nM; m++) {
                    for (int n = 0; n < nN; n++, index++) {
                        //Bp[p + oP][(m + oM) * nso_ + (n + oN)] = buffer[thread][index];
                        //Bp[p + oP][(n + oN) * nso_ + (m + oM)] = buffer[thread][index];
                        sum += Pao->get(m + oM, n + oN) * buffer[thread][index];
                    }
                }
                JpQ->set(p + oP, sum);
            }
        }
    }

}  // end of cis_ao_J_fht

//======================================================================
//     CIS AO Basis J: Second-half trans
//======================================================================
void DFOCC::cis_ao_J_sht(SharedTensor1d JdpQ, SharedTensor2d Jmn) {

    // Read in the basis set informations
    std::shared_ptr<BasisSet> auxiliary_ = get_basisset("DF_BASIS_CC");
    std::shared_ptr<BasisSet> primary_ = get_basisset("ORBITAL");

    std::shared_ptr<BasisSet> zero(BasisSet::zero_ao_basis_set());

    int nthreads = 1;
#ifdef _OPENMP
    nthreads = Process::environment.get_n_threads();
#endif

    std::shared_ptr<ERISieve> sieve_ = std::shared_ptr<ERISieve>(new ERISieve(primary_, cutoff));
    const std::vector<std::pair<int, int> >& shell_pairs = sieve_->shell_pairs();
    int npairs = shell_pairs.size();

    // => Memory Constraints <= //
    int max_rows;
    max_rows = auxiliary_->nshell();

    // => Block Sizing <= //
    std::vector<int> Pstarts;
    int counter = 0;
    Pstarts.push_back(0);
    for (int P = 0; P < auxiliary_->nshell(); P++) {
        int nP = auxiliary_->shell(P).nfunction();
        if (counter + nP > max_rows) {
            counter = 0;
            Pstarts.push_back(P);
        }
        counter += nP;
    }
    Pstarts.push_back(auxiliary_->nshell());

    // => Integrals <= //
    std::shared_ptr<IntegralFactory> rifactory2(new IntegralFactory(auxiliary_, zero, primary_, primary_));
    std::vector<std::shared_ptr<TwoBodyAOInt> > eri;
    std::vector<const double*> buffer;
    for (int t = 0; t < nthreads; t++) {
        eri.push_back(std::shared_ptr<TwoBodyAOInt>(rifactory2->eri()));
        buffer.push_back(eri[t]->buffer());
    }

    // => Master Loop <= //

    for (int block = 0; block < Pstarts.size() - 1; block++) {
        // > Sizing < //

        int Pstart = Pstarts[block];
        int Pstop = Pstarts[block + 1];
        int NP = Pstop - Pstart;

        int pstart = auxiliary_->shell(Pstart).function_index();
        int pstop = (Pstop == auxiliary_->nshell() ? nQ : auxiliary_->shell(Pstop).function_index());
        int np = pstop - pstart;

// > Integrals < //
#pragma omp parallel for schedule(dynamic) num_threads(nthreads)
        for (long int PMN = 0L; PMN < NP * npairs; PMN++) {
            int thread = 0;
#ifdef _OPENMP
            thread = omp_get_thread_num();
#endif

            int P = PMN / npairs + Pstart;
            int MN = PMN % npairs;
            int M = shell_pairs[MN].first;
            int N = shell_pairs[MN].second;

            eri[thread]->compute_shell(P, 0, M, N);

            int nP = auxiliary_->shell(P).nfunction();
            int oP = auxiliary_->shell(P).function_index();

            int nM = primary_->shell(M).nfunction();
            int oM = primary_->shell(M).function_index();

            int nN = primary_->shell(N).nfunction();
            int oN = primary_->shell(N).function_index();

            int index = 0;
            for (int m = 0; m < nM; m++) {
                for (int n = 0; n < nN; n++) {
                     double sum = 0.0;
                    for (int p = 0; p < nP; p++, index++) {
                        //Bp[p + oP][(m + oM) * nso_ + (n + oN)] = buffer[thread][index];
                        //Bp[p + oP][(n + oN) * nso_ + (m + oM)] = buffer[thread][index];
                        sum += JdpQ->get(p + oP) * buffer[thread][index];
                    }
                    Jmn->set(m + oM, n + oN, sum);
                }
            }
        }
    }

}  // end of cis_ao_J_sht



}  // namespace plugin_qdpt
}  // namespace psi
