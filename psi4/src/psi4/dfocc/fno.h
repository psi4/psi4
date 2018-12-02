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

#ifndef _dfocc_fno_h_
#define _dfocc_fno_h_

#include "psi4/libpsio/psio.h"
#include "psi4/libmints/typedefs.h"
#include "tensors.h"

using namespace psi;

namespace psi {
namespace dfoccwave {

class FnoCC;
typedef std::shared_ptr<FnoCC> SharedFnoCC;

class FnoCC {
   private:
    std::string name_;

    int nocc_;
    int naocc_;
    int nvir_;
    int navir_;
    int norb_;
    int nfrzc_;
    int nfrzv_;
    int naux_;
    int naob_;
    int navir_fno_;

    double cutoff_;
    double fno_cutoff_;
    double fno_percentage_;
    double Ecorr_;

    bool fno_by_perct_;
    bool fno_by_user_;

    // 1D
    SharedTensor1d diag_n_;

    // 2D
    SharedTensor2d GF_;
    SharedTensor2d Tvv_;
    SharedTensor2d Tmat_;  // eigenvectors of opdm
    SharedTensor2d Vmat_;  // V = C*T
    SharedTensor2d Gamma_;
    SharedTensor2d Fock_;
    SharedTensor2d bQ_;
    SharedTensor2d Corb_;

    // fno
    void compute_fno();

   public:
    // FnoCC(SharedWavefunction ref_wfn, Options &options);
    FnoCC(std::string name, int naux, int nao, int nfrzc, int naocc, int nvir, const SharedTensor2d& Corb,
          const SharedTensor2d& Fock, const SharedTensor2d& bQ, double fno_cutoff, double fno_perct, int navir_fno,
          bool fno_by_perct, bool fno_by_user);
    ~FnoCC();

    // Return occupation numbers
    SharedTensor1d occupation_numbers() const { return diag_n_; }

    // Return natural orbitals
    SharedTensor2d nat_orbs() const { return Tmat_; }

    // Return VV block NOs
    SharedTensor2d nat_orbs_vv() const { return Tvv_; }

    // Return NO coeff matrix
    SharedTensor2d no_coeff() const { return Vmat_; }

    // Return NO Fock matrix
    SharedTensor2d no_fock() const { return Fock_; }

    // Return number of active virtuals
    int nactvir() const { return navir_; }

    // Return number of frozen virtuals
    int nfrzvir() const { return nfrzv_; }

    // Return MP2 correction in the MO space
    double mp2_corr() const { return Ecorr_; }
};

}  // namespace dfoccwave
}  // namespace psi
#endif  // _dfocc_fno_h_
