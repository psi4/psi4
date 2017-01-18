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


#include "psi4/libmints/sieve.h"
#include "psi4/libmints/mintshelper.h"
#include "psi4/libmints/wavefunction.h"
#include "psi4/libciomr/libciomr.h"
#include "psi4/libqt/qt.h"
#include "psi4/lib3index/cholesky.h"
#include "psi4/psifiles.h"
#include "psi4/psi4-dec.h"
#include "psi4/libmints/integral.h"
#include "psi4/libmints/matrix.h"
#include "defines.h"
#include "dfocc.h"
#include "tensors.h"

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace psi;
using namespace std;

namespace psi{ namespace dfoccwave{

void DFOCC::trans_cd()
{
    // Read SO integrals
    bQso = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|mn)", nQ, nso_, nso_));
    bQso->read(psio_, PSIF_DFOCC_INTS, true, true);

    trans_ab = 1;
    if (orb_opt_ == "TRUE" || dertype == "FIRST" || ekt_ip_ == "TRUE") {
        // Form B(Q,ij)
        timer_on("Form B(Q,ij)");
        b_oo_cd();
        timer_off("Form B(Q,ij)");

        // Form B(Q,ia)
        timer_on("Form B(Q,ia)");
        b_ov_cd();
        timer_off("Form B(Q,ia)");

        // Form B(Q,ab)
        timer_on("Form B(Q,ab)");
        b_vv_cd();
        timer_off("Form B(Q,ab)");
    }

    else {
        // Form B(Q,ij)
        timer_on("Form B(Q,ij)");
        b_ij_cd();
        timer_off("Form B(Q,ij)");

        // Form B(Q,ia)
        timer_on("Form B(Q,ia)");
        b_ia_cd();
        timer_off("Form B(Q,ia)");

        // Form B(Q,ab)
        timer_on("Form B(Q,ab)");
        b_ab_cd();
        timer_off("Form B(Q,ab)");
    }
    bQso.reset();

    // Trans OEI
    timer_on("Trans OEI");
    trans_oei();
    timer_off("Trans OEI");
}

//=======================================================
//          trans for mp2 energy
//=======================================================
void DFOCC::trans_cd_mp2()
{
    // Read SO integrals
    bQso = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|mn)", nQ, nso_, nso_));
    bQso->read(psio_, PSIF_DFOCC_INTS, true, true);

    // Form B(Q,ia)
    trans_ab = 0;
    timer_on("Form B(Q,ia)");
    b_ia_cd();
    timer_off("Form B(Q,ia)");
    bQso.reset();
}

//=======================================================
//          DF CC
//=======================================================
void DFOCC::cd_ints()
{
    //outfile->Printf("\tComputing DF-BASIS-CC integrals... \n");

  // 1.  read scf 3-index integrals from disk

  // get ntri from sieve
  std::shared_ptr<ERISieve> sieve (new ERISieve(basisset_, options_.get_double("INTS_TOLERANCE")));
  const std::vector<std::pair<int, int> >& function_pairs = sieve->function_pairs();
  long int ntri_cd = function_pairs.size();

      // Cholesky

      // read integrals from disk if they were generated in the SCF
      if ( options_.get_str("SCF_TYPE") == "CD" ) {
          outfile->Printf("\tReading Cholesky vectors from disk ...\n");

          // ntri comes from sieve above
          psio_->open(PSIF_DFSCF_BJ,PSIO_OPEN_OLD);
          // Read the NAUX from the file
          psio_->read_entry(PSIF_DFSCF_BJ, "length", (char*)&nQ, sizeof(long int));
          std::shared_ptr<Matrix> Qmn = SharedMatrix(new Matrix("Qmn Integrals",nQ,ntri_cd));
          double** Qmnp = Qmn->pointer();
          psio_->read_entry(PSIF_DFSCF_BJ, "(Q|mn) Integrals", (char*) Qmnp[0], sizeof(double) * ntri_cd * nQ);
          psio_->close(PSIF_DFSCF_BJ,1);

          outfile->Printf("\tCholesky decomposition threshold: %8.2le\n", options_.get_double("CHOLESKY_TOLERANCE"));
          outfile->Printf("\tNumber of Cholesky vectors:   %5li\n",nQ);

          nQ_ref = nQ;
          bQso = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|mn)", nQ, nso_, nso_));
          for (long int mn = 0; mn < ntri_cd; mn++) {
              long int m = function_pairs[mn].first;
              long int n = function_pairs[mn].second;
              for (long int P = 0; P < nQ; P++) {
                  //Lp[P][m*nso_+n] = Qmnp[P][mn];
                  //Lp[P][n*nso_+m] = Qmnp[P][mn];
                  bQso->set(P, (m*nso_) + n, Qmnp[P][mn]);
                  bQso->set(P, (n*nso_) + m, Qmnp[P][mn]);
              }
          }
          bQso->write(psio_, PSIF_DFOCC_INTS, true, true);
      }// end if ( options_.get_str("SCF_TYPE") == "CD" )

      else {
          // generate Cholesky 3-index integrals
          outfile->Printf("\tGenerating Cholesky vectors ...\n");
          std::shared_ptr<BasisSet> primary = basisset();
          std::shared_ptr<IntegralFactory> integral (new IntegralFactory(primary,primary,primary,primary));
          double tol_cd = options_.get_double("CHOLESKY_TOLERANCE");
          std::shared_ptr<CholeskyERI> Ch (new CholeskyERI(std::shared_ptr<TwoBodyAOInt>(integral->eri()),cutoff,tol_cd,Process::environment.get_memory()));
          Ch->choleskify();
          nQ  = Ch->Q();
          nQ_ref = nQ;
          std::shared_ptr<Matrix> L = Ch->L();
          bQso = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|mn)", nQ, nso_, nso_));
          bQso->set(L);
          L.reset();
          bQso->write(psio_, PSIF_DFOCC_INTS, true, true);
          outfile->Printf("\tCholesky decomposition threshold: %8.2le\n", options_.get_double("CHOLESKY_TOLERANCE"));
          outfile->Printf("\tNumber of Cholesky vectors:   %5li\n",nQ);

      }

} // end df_corr

//=======================================================
//          form b(Q,ij) : active
//=======================================================
void DFOCC::b_ij_cd()
{
    bQnoA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|mI)", nQ, nso_ * naoccA));
    bQijA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|IJ)", nQ, naoccA * naoccA));
    bQnoA->contract(false, false, nQ * nso_, naoccA, nso_, bQso, CaoccA, 1.0, 0.0);
    bQijA->contract233(true, false, naoccA, naoccA, CaoccA, bQnoA, 1.0, 0.0);
    bQnoA.reset();
    bQijA->write(psio_, PSIF_DFOCC_INTS);
    bQijA.reset();

 if (reference_ == "UNRESTRICTED") {
    bQnoB = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|mi)", nQ, nso_ * naoccB));
    bQijB = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|ij)", nQ, naoccB * naoccB));
    bQnoB->contract(false, false, nQ * nso_, naoccB, nso_, bQso, CaoccB, 1.0, 0.0);
    bQijB->contract233(true, false, naoccB, naoccB, CaoccB, bQnoB, 1.0, 0.0);
    bQnoB.reset();
    bQijB->write(psio_, PSIF_DFOCC_INTS);
    bQijB.reset();
 }

} // end b_ij

//=======================================================
//          form b(Q,oo)
//=======================================================
void DFOCC::b_oo_cd()
{
    bQnoA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|mO)", nQ, nso_ * noccA));
    bQooA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|OO)", nQ, noccA, noccA));
    bQnoA->contract(false, false, nQ * nso_, noccA, nso_, bQso, CoccA, 1.0, 0.0);
    bQooA->contract233(true, false, noccA, noccA, CoccA, bQnoA, 1.0, 0.0);
    bQnoA.reset();
    bQooA->write(psio_, PSIF_DFOCC_INTS);
    bQooA->write(psio_, "DF_BASIS_SCF B (Q|OO)", PSIF_DFOCC_INTS);

    // Form active b(Q,ij)
    bQijA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|IJ)", nQ, naoccA, naoccA));
    bQijA->form_b_ij(nfrzc, bQooA);
    bQooA.reset();
    bQijA->write(psio_, PSIF_DFOCC_INTS);
    bQijA.reset();

 if (reference_ == "UNRESTRICTED") {
    bQnoB = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|mo)", nQ, nso_ * noccB));
    bQooB = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|oo)", nQ, noccB, noccB));
    bQnoB->contract(false, false, nQ * nso_, noccB, nso_, bQso, CoccB, 1.0, 0.0);
    bQooB->contract233(true, false, noccB, noccB, CoccB, bQnoB, 1.0, 0.0);
    bQnoB.reset();
    bQooB->write(psio_, PSIF_DFOCC_INTS);
    bQooB->write(psio_, "DF_BASIS_SCF B (Q|oo)", PSIF_DFOCC_INTS);

    // Form active b(Q,ij)
    bQijB = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|ij)", nQ, naoccB, naoccB));
    bQijB->form_b_ij(nfrzc, bQooB);
    bQooB.reset();
    bQijB->write(psio_, PSIF_DFOCC_INTS);
    bQijB.reset();
 }

} // end b_oo_cd

//=======================================================
//          form b(Q,ia) : active
//=======================================================
void DFOCC::b_ia_cd()
{
    bQnvA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|mA)", nQ, nso_ * navirA));
    bQiaA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|IA)", nQ, naoccA * navirA));
    bQnvA->contract(false, false, nQ * nso_, navirA, nso_, bQso, CavirA, 1.0, 0.0);
    bQiaA->contract233(true, false, naoccA, navirA, CaoccA, bQnvA, 1.0, 0.0);
    bQiaA->write(psio_, PSIF_DFOCC_INTS);
    bQnvA->write(psio_, PSIF_DFOCC_INTS);
    bQiaA.reset();
    bQnvA.reset();

 if (reference_ == "UNRESTRICTED") {
    bQnvB = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|ma)", nQ, nso_ * navirB));
    bQiaB = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|ia)", nQ, naoccB * navirB));
    bQnvB->contract(false, false, nQ * nso_, navirB, nso_, bQso, CavirB, 1.0, 0.0);
    bQiaB->contract233(true, false, naoccB, navirB, CaoccB, bQnvB, 1.0, 0.0);
    bQiaB->write(psio_, PSIF_DFOCC_INTS);
    bQnvB->write(psio_, PSIF_DFOCC_INTS);
    bQiaB.reset();
    bQnvB.reset();
 }

} // end b_ia

//=======================================================
//          form b(Q,ov) : all
//=======================================================
void DFOCC::b_ov_cd()
{
    bQnvA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|mV)", nQ, nso_ * nvirA));
    bQovA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|OV)", nQ, noccA, nvirA));
    bQnvA->contract(false, false, nQ * nso_, nvirA, nso_, bQso, CvirA, 1.0, 0.0);
    bQovA->contract233(true, false, noccA, nvirA, CoccA, bQnvA, 1.0, 0.0);
    bQovA->write(psio_, PSIF_DFOCC_INTS);
    bQovA->write(psio_, "DF_BASIS_SCF B (Q|OV)", PSIF_DFOCC_INTS);
    bQnvA->write(psio_, PSIF_DFOCC_INTS);
    bQnvA.reset();

    // Form active b(Q,ia)
    bQiaA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|IA)", nQ, naoccA, navirA));
    bQiaA->form_b_ia(nfrzc, bQovA);
    bQovA.reset();
    bQiaA->write(psio_, PSIF_DFOCC_INTS);
    bQiaA.reset();

 if (reference_ == "UNRESTRICTED") {
    bQnvB = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|mv)", nQ, nso_ * nvirB));
    bQovB = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|ov)", nQ, noccB, nvirB));
    bQnvB->contract(false, false, nQ * nso_, nvirB, nso_, bQso, CvirB, 1.0, 0.0);
    bQovB->contract233(true, false, noccB, nvirB, CoccB, bQnvB, 1.0, 0.0);
    bQovB->write(psio_, PSIF_DFOCC_INTS);
    bQovB->write(psio_, "DF_BASIS_SCF B (Q|ov)", PSIF_DFOCC_INTS);
    bQnvB->write(psio_, PSIF_DFOCC_INTS);
    bQnvB.reset();

    // Form active b(Q,ia)
    bQiaB = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|ia)", nQ, naoccB, navirB));
    bQiaB->form_b_ia(nfrzc, bQovB);
    bQovB.reset();
    bQiaB->write(psio_, PSIF_DFOCC_INTS);
    bQiaB.reset();
 }

} // end b_ov_cd

//=======================================================
//          form b(Q,ab) : active
//=======================================================
void DFOCC::b_ab_cd()
{
    bQabA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|AB)", nQ, navirA, navirA));
    bQnvA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|mA)", nQ, nso_ * navirA));
    bQnvA->read(psio_, PSIF_DFOCC_INTS);
    bQabA->contract233(true, false, navirA, navirA, CavirA, bQnvA, 1.0, 0.0);
    bQnvA.reset();
    bQabA->write(psio_, PSIF_DFOCC_INTS, true, true);
    bQabA.reset();

 if (reference_ == "UNRESTRICTED") {
    bQabB = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|ab)", nQ, navirB, navirB));
    bQnvB = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|ma)", nQ, nso_ * navirB));
    bQnvB->read(psio_, PSIF_DFOCC_INTS);
    bQabB->contract233(true, false, navirB, navirB, CavirB, bQnvB, 1.0, 0.0);
    bQnvB.reset();
    bQabB->write(psio_, PSIF_DFOCC_INTS, true, true);
    bQabB.reset();
 }

} // end b_ab

//=======================================================
//          form b(Q,vv) : all
//=======================================================
void DFOCC::b_vv_cd()
{
    bQvvA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|VV)", nQ, nvirA, nvirA));
    bQnvA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|mV)", nQ, nso_ * nvirA));
    bQnvA->read(psio_, PSIF_DFOCC_INTS);
    bQvvA->contract233(true, false, nvirA, nvirA, CvirA, bQnvA, 1.0, 0.0);
    bQnvA.reset();
    bQvvA->write(psio_, PSIF_DFOCC_INTS, true, true);
    bQvvA->write(psio_, "DF_BASIS_SCF B (Q|VV)", PSIF_DFOCC_INTS, true, true);

    // Form active b(Q,ab)
    bQabA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|AB)", nQ, navirA, navirA));
    bQabA->form_b_ab(bQvvA);
    bQvvA.reset();
    bQabA->write(psio_, PSIF_DFOCC_INTS, true, true);
    bQabA.reset();

 if (reference_ == "UNRESTRICTED") {
    bQvvB = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|vv)", nQ, nvirB, nvirB));
    bQnvB = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|mv)", nQ, nso_ * nvirB));
    bQnvB->read(psio_, PSIF_DFOCC_INTS);
    bQvvB->contract233(true, false, nvirB, nvirB, CvirB, bQnvB, 1.0, 0.0);
    bQnvB.reset();
    bQvvB->write(psio_, PSIF_DFOCC_INTS, true, true);
    bQvvB->write(psio_, "DF_BASIS_SCF B (Q|vv)", PSIF_DFOCC_INTS, true, true);

    // Form active b(Q,ab)
    bQabB = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|ab)", nQ, navirB, navirB));
    bQabB->form_b_ab(bQvvB);
    bQvvB.reset();
    bQabB->write(psio_, PSIF_DFOCC_INTS, true, true);
    bQabB.reset();
 }

} // end b_vv_cd

}} // Namespaces

