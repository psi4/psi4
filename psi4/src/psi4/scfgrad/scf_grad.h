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

#ifndef SCF_GRAD_H
#define SCF_GRAD_H

#include "psi4/libfock/jk.h"
#include "psi4/libmints/wavefunction.h"
#include "psi4/libmints/typedefs.h"
#include "psi4/libscf_solver/hf.h"

namespace psi {
class SuperFunctional;
class VBase;
namespace scf {
    class RHF;
    class UHF;
}

namespace scfgrad {

class SCFDeriv : public Wavefunction {

protected:

    /// Common initialization
    void common_init();
    std::shared_ptr<SuperFunctional> functional_;
    std::shared_ptr<VBase> potential_;
    std::map<std::string, SharedMatrix> gradients_;
    std::map<std::string, SharedMatrix> hessians_;

public:
    SCFDeriv(std::shared_ptr<scf::HF> ref_wfn, Options& options);
    ~SCFDeriv() override;

    double compute_energy() override { throw PSIEXCEPTION("SCFDeriv not implemented for the requested reference type."); }
    virtual SharedMatrix compute_gradient() override;
    virtual SharedMatrix compute_hessian() override;
    virtual SharedMatrix hessian_response() { throw PSIEXCEPTION("SCFDeriv not implemented for the requested reference type."); }
};

class RSCFDeriv : public SCFDeriv {
protected:
    std::shared_ptr<scf::RHF> rhf_wfn_;
public:
    RSCFDeriv(std::shared_ptr<scf::RHF> rhf_wfn, Options& options): SCFDeriv(std::static_pointer_cast<scf::HF>(rhf_wfn), options), rhf_wfn_(rhf_wfn) {}
    ~RSCFDeriv() override {}
    virtual SharedMatrix hessian_response() override;
};

class USCFDeriv : public SCFDeriv {
protected:
    std::shared_ptr<scf::UHF> uhf_wfn_;
    SharedMatrix hessian_response_spin(int spin);
    void overlap_deriv(std::shared_ptr<Matrix> C, 
                              std::shared_ptr<Matrix> Cocc,
                              std::shared_ptr<Matrix> Cvir,
                              int nso, int nocc, int nvir, bool alpha);

    void kinetic_deriv(std::shared_ptr<Matrix> C, 
                       std::shared_ptr<Matrix> Cocc,
                       int nso, int nocc, int nvir, bool alpha);

    void potential_deriv(std::shared_ptr<Matrix> C, 
                         std::shared_ptr<Matrix> Cocc,
                         int nso, int nocc, int nvir, bool alpha);

#ifdef USING_ecpint
    void ecp_deriv(std::shared_ptr<Matrix> C, 
                   std::shared_ptr<Matrix> Cocc,
                   int nso, int nocc, int nvir, bool alpha);
#endif

    // Compute the JK contribution to the the Fock derivative onthe
    //   right-side of the CP-SCF equations.
    void JK_deriv1(std::shared_ptr<Matrix> D1,
                   std::shared_ptr<Matrix> C1, 
                   std::shared_ptr<Matrix> C1occ,
                   std::shared_ptr<Matrix> D2, 
                   int nso, int nocc, int nvir, bool alpha);

    // Compute the JK contribution to the the overlap derivative *
    //   TEI term  on the right-side of the CP-SCF equations.
    void JK_deriv2(std::shared_ptr<JK> jk, int mem, 
                   std::shared_ptr<Matrix> Ca,
                   std::shared_ptr<Matrix> Caocc,
                   std::shared_ptr<Matrix> Cb,
                   std::shared_ptr<Matrix> Cbocc,
                   int nso, int naocc, int nbocc, int navir);

    void VXC_deriv(std::shared_ptr<Matrix> Ca,
                   std::shared_ptr<Matrix> Caocc,
                   std::shared_ptr<Matrix> Cb,
                   std::shared_ptr<Matrix> Cbocc,
                   int nso, int naocc, int nbocc, int navir);

    void assemble_Fock(int nocc, int nvir, bool alpha);

    void assemble_B(std::shared_ptr<Vector> eocc, int nocc, int nvir, bool alpha);
    void assemble_U(int nocc, int nvir, bool alpha);
    void assemble_Q(std::shared_ptr<JK> jk,
                    std::shared_ptr<Matrix> C1, 
                    std::shared_ptr<Matrix> C1occ,
                    std::shared_ptr<Matrix> C2, 
                    std::shared_ptr<Matrix> C2occ, 
                    int nso, int n1occ, int n2occ, int nvir);

    void dipole_derivatives(std::shared_ptr<Matrix> C1, 
                            std::shared_ptr<Matrix> C1occ,
                            std::shared_ptr<Matrix> C2, 
                            std::shared_ptr<Matrix> C2occ,
                            std::shared_ptr<Matrix> Dt,
                            int nso, int n1occ, int n2occ, int n1vir);

public:
    USCFDeriv(std::shared_ptr<scf::UHF> uhf_wfn, Options& options): SCFDeriv(std::static_pointer_cast<scf::HF>(uhf_wfn), options), uhf_wfn_(uhf_wfn) {}
    ~USCFDeriv() override {}
    virtual SharedMatrix hessian_response() override;
};

}} // Namespaces

#endif
