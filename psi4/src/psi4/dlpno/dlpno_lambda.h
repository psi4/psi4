/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2025 The Psi4 Developers.
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

#ifndef PSI4_SRC_DLPNO_H_
#define PSI4_SRC_DLPNO_H_

#include "sparse.h"
#include "dlpno.h"

#include "psi4/libmints/wavefunction.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/vector.h"
#include "psi4/libqt/qt.h"
#include "psi4/libpsio/psio.h"
#include "psi4/psifiles.h"

#include <map>
#include <tuple>
#include <string>
#include <unordered_map>

namespace psi {
namespace dlpno {

enum class DLPNOLambdaMethod { CCSD_L };

class DLPNOCCSD_Lambda : public DLPNOCCSD {
   protected:

    /// Lambda Amplitudes

    std::vector<SharedMatrix> lambda_ia_; // Dimensions: [nocc * (npno_[ii], 1)]
    std::vector<SharedMatrix> lambda_iajb_; // Dimensions: [n_lmo_pairs * (n_pno_[ij], n_pno_[ij])]

    // => Lambda CCSD Integrals (Fock matrices and ERIs) <= //
    std::vector<SharedMatrix> delta_imae_tilde_; // Toth eq. 26a
    SharedMatrix F_im_double_tilde_; // Toth Eq. 26b
    std::vector<SharedMatrix> F_vv_double_tilde_;

    std::vector<SharedMatrix> M_imae_; // M(i m | a_{mm} e_{mm})

    // => Lambda CCSD Intermediates <= //

    std::vector<SharedMatrix> zizi_is_delinquent_;
    std::vector<SharedMatrix> john_bigback_;

    SharedMatrix rho_oo_; // Toth Section IIIA
    std::vector<SharedMatrix> rho_vv_; // Toth Section IIIA
    std::vector<SharedMatrix> M_imae_tilde_; // Toth Eq. 24
    std::vector<std::vector<SharedMatrix>> F_fcia_hat_; // Toth Eq. 33
    std::vector<std::vector<SharedMatrix>> F_knia_hat_; // Toth Eq. 34

    // => Computing integrals <= //

    /// A function to estimate memory costs for lambda CCSD
    void estimate_memory();
    /// Compute some of the intermediates required in lambda DLPNO CCSD
    void compute_lambda_intermediates();
    /// Computes the specific integral types needed for lambda CCSD (PNO basis)
    void compute_lambda_pno_integrals();

    // => Lambda CCSD intermediates <= //

    /// New defined intermediate in Toth Eq. 50b
    std::vector<SharedMatrix> compute_alpha_ijkl();
    /// Toth Eq. XX
    std::vector<SharedMatrix> weird_stuff();
    /// Toth Eq. XY
    std::vector<SharedMatrix> strawberry_raspberry_shortcake();
    /// Toth Eq. ZeeZee
    std::vector<SharedMatrix> lasagne();

    /// compute singles residual in lambda CCSD equations, Toth Eq. 26-36
    void compute_L_ia(std::vector<SharedMatrix>& L_ia, std::vector<std::vector<SharedMatrix>>& L_ia_buffer);
    /// compute doubles residual in lambda CCSD equations, Toth Eq. 37-55
    void compute_L_iajb(std::vector<SharedMatrix> &L_iajb, std::vector<SharedMatrix>& Ln_iajb);

    /// Solves lambda_ia and lambda_iajb in the PNO basis
    void lambda_ccsd_iterations();

    void print_header();
    void print_results();

   public:
    DLPNOCCSD_Lambda(SharedWavefunction ref_wfn, Options& options);
    ~DLPNOCCSD_Lambda() override;

    double compute_energy();
};

}
}