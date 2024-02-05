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

#ifndef __rohf_psi_h__
#define __rohf_psi_h__

#include <vector>
#include "psi4/libpsio/psio.hpp"
#include "hf.h"

namespace psi {
class VBase;
namespace scf {

class ROHF : public HF {
   protected:
    SharedMatrix moFeff_;
    SharedMatrix soFeff_;
    SharedMatrix Dt_;
    SharedMatrix Da_old_;
    SharedMatrix Db_old_;
    SharedMatrix Dt_old_;
    SharedMatrix Ct_;
    SharedMatrix Ga_;
    SharedMatrix Gb_;
    SharedMatrix Ka_;
    SharedMatrix Kb_;
    SharedMatrix moFa_;
    SharedMatrix moFb_;

    void form_initial_F() override;
    void form_initial_C() override;
    double compute_initial_E() override;
    void prepare_canonical_orthogonalization() override;
    void semicanonicalize() override;

    // Second-order convergence code
    void Hx(SharedMatrix x, SharedMatrix ret);

    void format_guess() override;

    void common_init();
    void setup_potential() override;

   public:
    ROHF(SharedWavefunction ref_wfn, std::shared_ptr<SuperFunctional> functional);
    ROHF(SharedWavefunction ref_wfn, std::shared_ptr<SuperFunctional> functional, Options& options,
         std::shared_ptr<PSIO> psio);
    ~ROHF() override;

    SharedMatrix soFeff() const { return soFeff_; }
    SharedMatrix moFeff() const { return moFeff_; }
    SharedMatrix moFa() const { return moFa_; }
    SharedMatrix moFb() const { return moFb_; }
    SharedMatrix Ct() const {return Ct_; }

    void save_density_and_energy() override;

    void form_C(double shift = 0.0) override;
    void form_D() override;
    void form_F() override;
    void form_G() override;
    double compute_E() override;
    void finalize() override;

    void compute_SAD_guess(bool natorb) override;

    void damping_update(double) override;
    int soscf_update(double soscf_conv, int soscf_min_iter, int soscf_max_iter, int soscf_print) override;
    bool stability_analysis() override;

    std::shared_ptr<VBase> V_potential() const override { return nullptr; };

    std::shared_ptr<ROHF> c1_deep_copy(std::shared_ptr<BasisSet> basis);
};
}  // namespace scf
}  // namespace psi

#endif
