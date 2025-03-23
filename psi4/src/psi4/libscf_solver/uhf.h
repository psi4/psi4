/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2024 The Psi4 Developers.
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

#ifndef __math_test_uhf_h__
#define __math_test_uhf_h__

#include "psi4/libpsio/psio.hpp"
#include "psi4/libfock/v.h"
#include "hf.h"

namespace psi {
namespace scf {

class UHF : public HF {
   protected:
    SharedMatrix Da_old_, Db_old_;
    SharedMatrix Ga_, Gb_, J_, Ka_, Kb_, wKa_, wKb_;

    std::shared_ptr<UV> potential_;

    double compute_initial_E() override;

    void common_init();
    void setup_potential() override;

    // Guess mix performed?
    bool mix_performed_;

    // Scaling factor for orbital rotation
    double step_scale_;
    // Increment to explore different scaling factors
    double step_increment_;
    // Stability eigenvalue, for doing smart eigenvector following
    double stab_val;

    // Compute UHF NOs
    void compute_nos();

    // Second-order convergence code
    void Hx(SharedMatrix x_a, SharedMatrix IFock_a, SharedMatrix Cocc_a, SharedMatrix Cvir_a, SharedMatrix ret_a,
            SharedMatrix x_b, SharedMatrix IFock_b, SharedMatrix Cocc_b, SharedMatrix Cvir_b, SharedMatrix ret_b);

   public:
    UHF(SharedWavefunction ref_wfn, std::shared_ptr<SuperFunctional> functional);
    UHF(SharedWavefunction ref_wfn, std::shared_ptr<SuperFunctional> functional, Options& options,
        std::shared_ptr<PSIO> psio);
    ~UHF() override;

    virtual bool same_a_b_orbs() const { return false; }
    virtual bool same_a_b_dens() const { return false; }

    void save_density_and_energy() override;

    void form_C(double shift = 0.0) override;
    void form_D() override;
    void form_F() override;
    void form_G() override;
    void form_V() override;
    double compute_E() override;
    void finalize() override;

    void openorbital_scf() override;

    void damping_update(double) override;
    int soscf_update(double soscf_conv, int soscf_min_iter, int soscf_max_iter, int soscf_print) override;

    std::shared_ptr<VBase> V_potential() const override { return potential_; };

    /// Hessian-vector computers and solvers
    std::vector<SharedMatrix> onel_Hx(std::vector<SharedMatrix> x) override;
    std::vector<SharedMatrix> twoel_Hx(std::vector<SharedMatrix> x, bool combine = true,
                                       std::string return_basis = "MO") override;
    std::vector<SharedMatrix> cphf_Hx(std::vector<SharedMatrix> x) override;
    std::vector<SharedMatrix> cphf_solve(std::vector<SharedMatrix> x_vec, double conv_tol = 1.e-4, int max_iter = 10,
                                         int print_lvl = 1) override;

    std::shared_ptr<UHF> c1_deep_copy(std::shared_ptr<BasisSet> basis);
};
}  // namespace scf
}  // namespace psi

#endif
