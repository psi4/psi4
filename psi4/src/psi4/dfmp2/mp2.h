/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2023 The Psi4 Developers.
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

#ifndef TOTAL_DFMP2_H
#define TOTAL_DFMP2_H

#include "psi4/libmints/wavefunction.h"
#include <map>

namespace psi {

class PSIO;

namespace dfmp2 {

// References to DiStasio are to J Comput Chem 28: 839–856, 2007; DOI: 10.1002/jcc.20604

class DFMP2 : public Wavefunction {
   protected:
    // Auxiliary basis
    std::shared_ptr<BasisSet> ribasis_;
    // Gradients map
    std::map<std::string, SharedMatrix> gradients_;

    // Same-spin scale
    double sss_;
    // Opposite-spin scale
    double oss_;

    void common_init();
    // Common printing of energies/SCS
    virtual void print_energies();
    virtual void print_gradients();

    // Print header/reference information
    virtual void print_header() = 0;
    // Form the (A|ia) = (A|mn) C_mi C_na tensor(s)
    virtual void form_Aia() = 0;
    // Apply the fitting (Q|ia) = J_QA^-1/2 (A|ia)
    virtual void form_Bia() = 0;
    // Apply the fitting (Q|ia) = J_QA^-1/2 (A|ia) and J_QA^-1 (A|ia)
    virtual void form_Bia_Cia() = 0;
    // Transpose the B integrals to (ai|Q)
    virtual void form_Bia_transpose() = 0;
    // Form the energy contributions
    virtual void form_energy() = 0;
    // Form the VV block of the correlation OPDM (DiStasio 8) and Gamma_ia^Q (DiStasio 2)
    virtual void form_Pab() = 0;
    // Form the OO block of the correlation OPDM (DiStasio 7)
    virtual void form_Pij() = 0;
    // Form the small gamma; Eq. 3 of DiStasio
    virtual void form_gamma() = 0;
    // Transpose the G
    virtual void form_G_transpose() = 0;
    // Form the (A|B)^x contribution to the gradient; Term 2 of DiStasio 1
    virtual void form_AB_x_terms() = 0;
    // Form the (A|mn)^x contribution to the gradient; Term 1 of DiStasio 1
    virtual void form_Amn_x_terms() = 0;
    // Form the L_μa and L_μi matrices; DiStasio 19 and 20
    virtual void form_L() = 0;
    // Form the unrelaxed correlation OPDM; Compute DiStasio 6 and 9; Assemble DiStasio 6-9 into one matrix
    virtual void form_P() = 0;
    // Form part of the unrelaxed correlation EWDM; DiStasio 11-13... plus fudge factors
    virtual void form_W() = 0;
    // Form the full Lagrangian, solve the Z-vector equations, and apply final corrections to W and P (DiStasio 10, 16, 17)
    virtual void form_Z() = 0;
    // Manage the formation of W and P contributions to the gradient
    virtual void form_gradient() = 0;

    // Compute singles correction [nonzero for ROHF-MBPT(2) or dual-basis]
    virtual void form_singles();
    // Apply the fitting and transposition to a given disk entry Aia tensor
    virtual void apply_fitting(SharedMatrix Jm12, size_t file, size_t naux, size_t nia);
    // Apply the fitting again to a given disk entry Qia tensor
    virtual void apply_fitting_grad(SharedMatrix Jm12, size_t file, size_t naux, size_t nia);
    // Form the inverse square root of the fitting metric, or read it off disk
    virtual SharedMatrix form_inverse_metric();
    // Form an abstract gamma
    virtual void apply_gamma(size_t file, size_t naux, size_t nia);
    // Form a transposed copy of G_ia^P
    virtual void apply_G_transpose(size_t file, size_t naux, size_t nia);
    // Form a transposed copy of iaQ
    virtual void apply_B_transpose(size_t file, size_t naux, size_t naocc, size_t navir);

    // Debugging-routine: prints block sizing
    void block_status(std::vector<int> inds, const char* file, int line);
    void block_status(std::vector<size_t> inds, const char* file, int line);

    void compute_opdm_and_nos(const SharedMatrix Dnosym, SharedMatrix Dso, SharedMatrix Cno, SharedVector occ);

   public:
    DFMP2(SharedWavefunction ref_wfn, Options& options, std::shared_ptr<PSIO> psio);
    ~DFMP2() override;

    double compute_energy() override;
    SharedMatrix compute_gradient() override;
};

class RDFMP2 : public DFMP2 {
   protected:
    SharedMatrix Cfocc_;
    SharedMatrix Caocc_;
    SharedMatrix Cavir_;
    SharedMatrix Cfvir_;

    SharedVector eps_focc_;
    SharedVector eps_aocc_;
    SharedVector eps_avir_;
    SharedVector eps_fvir_;

    void common_init();

    // Print additional header
    void print_header() override;
    // Form the (A|ia) = (A|mn) C_mi C_ma tensor(s)
    void form_Aia() override;
    // Apply the fitting (Q|ia) = J_QA^-1/2 (A|ia)
    void form_Bia() override;
    // Apply the fitting (Q|ia) = J_QA^-1/2 (A|ia) and J_QA^-1 (A|ia)
    void form_Bia_Cia() override;
    // Transpose the integrals to (ai|Q)
    void form_Bia_transpose() override;
    // Form the energy contributions
    void form_energy() override;
    // Form the energy contributions and gradients
    void form_Pab() override;
    // Form the energy contributions and gradients
    void form_Pij() override;
    // Form the small gamma
    void form_gamma() override;
    // Transpose the G
    void form_G_transpose() override;
    // Form the (A|B)^x contribution to the gradient
    void form_AB_x_terms() override;
    // Form the (A|mn)^x contribution to the gradient
    void form_Amn_x_terms() override;
    // Form the Lma and Lmi matrices
    void form_L() override;
    // Form the unrelaxed OPDM
    void form_P() override;
    // Form the unrelaxed energy-weighted OPDM
    void form_W() override;
    // Form the full Lagrangian, solve the Z-vector equations, and apply final corrections to W and P
    void form_Z() override;
    // Manage the formation of W and P contributions to the gradient
    void form_gradient() override;

   public:
    RDFMP2(SharedWavefunction ref_wfn, Options& options, std::shared_ptr<PSIO> psio);
    ~RDFMP2() override;
};

class UDFMP2 : public DFMP2 {
   protected:
    SharedMatrix Caocc_a_;
    SharedMatrix Cavir_a_;
    SharedMatrix Caocc_b_;
    SharedMatrix Cavir_b_;

    SharedVector eps_aocc_a_;
    SharedVector eps_avir_a_;
    SharedVector eps_aocc_b_;
    SharedVector eps_avir_b_;

    void common_init();

    // Print additional header
    void print_header() override;
    // Form the (A|ia) = (A|mn) C_mi C_na tensor(s)
    void form_Aia() override;
    // Apply the fitting (Q|ia) = J_QA^-1/2 (A|ia)
    void form_Bia() override;
    // Apply the fitting (Q|ia) = J_QA^-1/2 (A|ia) and J_QA^-1 (A|ia)
    void form_Bia_Cia() override;
    // Transpose the integrals to (ai|Q)
    void form_Bia_transpose() override;
    // Form the energy contributions
    void form_energy() override;
    // Form the energy contributions and gradients
    void form_Pab() override;
    // Form the energy contributions and gradients
    void form_Pij() override;
    // Form the small gamma
    void form_gamma() override;
    // Transpose the G
    void form_G_transpose() override;
    // Form the (A|B)^x contribution to the gradient
    void form_AB_x_terms() override;
    // Form the (A|mn)^x contribution to the gradient
    void form_Amn_x_terms() override;
    // Form the Lma and Lmi matrices
    void form_L() override;
    // Form the unrelaxed OPDM
    void form_P() override;
    // Form the unrelaxed energy-weighted OPDM
    void form_W() override;
    // Form the full Lagrangian, solve the Z-vector equations, and apply final corrections to W and P
    void form_Z() override;
    // Manage the formation of W and P contributions to the gradient
    void form_gradient() override;

   public:
    UDFMP2(SharedWavefunction ref_wfn, Options& options, std::shared_ptr<PSIO> psio);
    ~UDFMP2() override;
};

class RODFMP2 : public UDFMP2 {
   protected:
    void common_init();

    // Print additional header
    void print_header() override;

   public:
    RODFMP2(SharedWavefunction ref_wfn, Options& options, std::shared_ptr<PSIO> psio);
    ~RODFMP2() override;
};

}  // namespace dfmp2
}  // namespace psi

#endif
