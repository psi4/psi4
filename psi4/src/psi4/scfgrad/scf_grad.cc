/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2022 The Psi4 Developers.
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

#include "scf_grad.h"

#include <algorithm>
#include <numeric>
#include <sstream>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "psi4/libqt/qt.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/integral.h"
#include "psi4/libmints/vector.h"
#include "psi4/libmints/extern.h"
#include "psi4/libmints/mintshelper.h"
#include "psi4/psi4-dec.h"
#include "psi4/libfock/v.h"
#include "psi4/libfunctional/superfunctional.h"
#include "psi4/libdisp/dispersion.h"
#include "psi4/libscf_solver/hf.h"
#include "psi4/libscf_solver/uhf.h"
#include "psi4/libscf_solver/rhf.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libpsi4util/process.h"

#include "jk_grad.h"

#ifdef USING_BrianQC

#include <brian_types.h>

extern bool brianCPHFFlag;
extern BrianCookie brianCookie;
extern bool brianEnable;
extern bool brianEnableDFT;

#endif

namespace psi {
namespace scfgrad {

SCFDeriv::SCFDeriv(SharedWavefunction ref_wfn, Options& options) :
    Wavefunction(options)
{
    shallow_copy(ref_wfn);
    common_init();
    scf::HF* scfwfn = (scf::HF*)ref_wfn.get();
    functional_ = scfwfn->functional();
    potential_ = scfwfn->V_potential();
    if (ref_wfn->has_array_variable("-D Gradient")) {
        gradients_["-D Gradient"] = ref_wfn->array_variable("-D Gradient");
    }
    if (ref_wfn->has_array_variable("-D Hessian")) {
        hessians_["-D Hessian"] = ref_wfn->array_variable("-D Hessian");
    }

}
SCFDeriv::~SCFDeriv()
{
}
void SCFDeriv::common_init()
{

    module_ = "scf";
    print_ = options_.get_int("PRINT");
    debug_ = options_.get_int("DEBUG");
}
SharedMatrix SCFDeriv::compute_gradient()
{
    // => Echo <= //

    outfile->Printf( "\n");
    outfile->Printf( "         ------------------------------------------------------------\n");
    outfile->Printf( "                                   SCF GRAD                          \n");
    outfile->Printf( "                          Rob Parrish, Justin Turney,                \n");
    outfile->Printf( "                       Andy Simmonett, and Alex Sokolov              \n");
    outfile->Printf( "         ------------------------------------------------------------\n\n");

    outfile->Printf( "  ==> Geometry <==\n\n");
    molecule_->print();

    outfile->Printf( "  Nuclear repulsion = %20.15f\n", basisset_->molecule()->nuclear_repulsion_energy(dipole_field_strength_));
    outfile->Printf( "\n");

    outfile->Printf( "  ==> Basis Set <==\n\n");
    basisset_->print_by_level("outfile", print_);

    // => Registers <= //

    // std::map<std::string, SharedMatrix> gradients;

    std::vector<std::string> gradient_terms;
    gradient_terms.push_back("Nuclear");
    gradient_terms.push_back("Core");
    gradient_terms.push_back("Overlap");
    gradient_terms.push_back("Coulomb");
    gradient_terms.push_back("Exchange");
    gradient_terms.push_back("Exchange,LR");
    gradient_terms.push_back("XC");
    gradient_terms.push_back("-D Gradient");
    gradient_terms.push_back("Total");

    // => Densities <= //
    SharedMatrix Da;
    SharedMatrix Db;
    SharedMatrix Dt;

    Da = Da_subset("AO");
    if (options_.get_str("REFERENCE") == "RHF" || options_.get_str("REFERENCE") == "RKS") {
        Db = Da;
    } else {
        Db = Db_subset("AO");
    }
    Dt = SharedMatrix(Da->clone());
    Dt->add(Db);
    Dt->set_name("Dt");

    // => Occupations (AO) <= //
    SharedMatrix Ca_occ;
    SharedMatrix Cb_occ;
    SharedVector eps_a_occ;
    SharedVector eps_b_occ;

    Ca_occ = Ca_subset("AO", "OCC");
    eps_a_occ = epsilon_a_subset("AO", "OCC");
    if (options_.get_str("REFERENCE") == "RHF" || options_.get_str("REFERENCE") == "RKS") {
        Cb_occ = Ca_occ;
        eps_b_occ = eps_a_occ;
    } else {
        Cb_occ = Cb_subset("AO", "OCC");
        eps_b_occ = epsilon_b_subset("AO", "OCC");
    }

    // => Potential/Functional <= //
    if (functional_->needs_xc()) {
        if (options_.get_str("REFERENCE") == "RKS") {
            potential_->set_D({Da_});
        } else {
            potential_->set_D({Da_, Db_});
        }
    }

    // => Sizings <= //
    int natom = molecule_->natom();
    int nso = basisset_->nbf();
    int nalpha = Ca_occ->colspi()[0];
    int nbeta = Cb_occ->colspi()[0];

    // => Nuclear Gradient <= //
    gradients_["Nuclear"] = SharedMatrix(molecule_->nuclear_repulsion_energy_deriv1(dipole_field_strength_).clone());
    gradients_["Nuclear"]->set_name("Nuclear Gradient");

    auto mints = std::make_shared<MintsHelper>(basisset_, options_);

    // => V T Perturbation Gradients <= //
    timer_on("Grad: V T Perturb");
    gradients_["Core"] = mints->core_hamiltonian_grad(Dt);
    timer_off("Grad: V T Perturb");

    // If an external field exists, add it to the one-electron Hamiltonian
    if (external_pot_) {
        gradient_terms.push_back("External Potential");
        timer_on("Grad: External");
        gradients_["External Potential"] = external_pot_->computePotentialGradients(basisset_, Dt);
        timer_off("Grad: External");
    }  // end external

    // => Overlap Gradient <= //
    timer_on("Grad: S");
    {
        // Energy weighted density matrix
        SharedMatrix W(Da->clone());
        W->set_name("W");

        // Alpha
        auto tmp = Ca_occ->clone();
        for (size_t i = 0; i < nalpha; i++){
            tmp->scale_column(0, i, eps_a_occ->get(i));
        }
        W->gemm(false, true, 1.0, tmp, Ca_occ, 0.0);

        // Beta
        tmp->copy(Cb_occ);
        for (size_t i = 0; i < nbeta; i++){
            tmp->scale_column(0, i, eps_b_occ->get(i));
        }
        W->gemm(false, true, 1.0, tmp, Cb_occ, 1.0);

        gradients_["Overlap"] = mints->overlap_grad(W);
        gradients_["Overlap"]->scale(-1.0);
    }
    timer_off("Grad: S");

    // => Two-Electron Gradient <= //
    timer_on("Grad: JK");

    auto jk = JKGrad::build_JKGrad(1, mintshelper_);
    jk->set_memory((size_t) (options_.get_double("SCF_MEM_SAFETY_FACTOR") * memory_ / 8L));

    jk->set_Ca(Ca_occ);
    jk->set_Cb(Cb_occ);
    jk->set_Da(Da);
    jk->set_Db(Db);
    jk->set_Dt(Dt);

    jk->set_do_J(true);
    if (functional_->is_x_hybrid()) {
        jk->set_do_K(true);
    } else {
        jk->set_do_K(false);
    }
    if (functional_->is_x_lrc()) {
        jk->set_do_wK(true);
        jk->set_omega(functional_->x_omega());
    } else {
        jk->set_do_wK(false);
    }

    jk->print_header();
    jk->compute_gradient();
    
    double alpha = functional_->x_alpha();
    double beta = functional_->x_beta();
    
#ifdef USING_BrianQC
    if (brianEnable and brianEnableDFT) {
        // BrianQC multiplies with the exact exchange factors inside the Fock building, so we must not do it here
        alpha = 1.0;
        beta = 1.0;
    }
#endif

    std::map<std::string, SharedMatrix>& jk_gradients = jk->gradients();
    gradients_["Coulomb"] = jk_gradients["Coulomb"];
    if (functional_->is_x_hybrid()) {
        gradients_["Exchange"] = jk_gradients["Exchange"];
        gradients_["Exchange"]->scale(-alpha);
    }
    if (functional_->is_x_lrc()) {
        gradients_["Exchange,LR"] = jk_gradients["Exchange,LR"];
        gradients_["Exchange,LR"]->scale(-beta);
    }
    timer_off("Grad: JK");

    // => XC Gradient <= //
    timer_on("Grad: XC");
    if (functional_->needs_xc()) {
        potential_->print_header();
        gradients_["XC"] = potential_->compute_gradient();
    }
    timer_off("Grad: XC");

    // => Total Gradient <= //
    SharedMatrix total = SharedMatrix(gradients_["Nuclear"]->clone());
    total->zero();

    for (int i = 0; i < gradient_terms.size(); i++) {
        if (gradients_.count(gradient_terms[i])) {
            total->add(gradients_[gradient_terms[i]]);
        }
    }

    // Symmetrize
    total->symmetrize_gradient(molecule_);

    gradients_["Total"] = total;
    gradients_["Total"]->set_name("Total Gradient");

    // => Final Printing <= //
    if (print_ > 1) {
        for (int i = 0; i < gradient_terms.size(); i++) {
            if (gradients_.count(gradient_terms[i])) {
                gradients_[gradient_terms[i]]->print_atom_vector();
            }
        }
    } else {
        gradients_["Total"]->print_atom_vector();
    }


    return gradients_["Total"];
}

void process_buffers(double **Hess, const std::vector<double> &Hvalues, int atom1, int atom2, int natoms,
                     bool shells_equivalent, bool is_potential) {
    // the number of buffers expected for each dimension; 2 for (P| and |Q), plus natoms more if nuclear
    size_t num_buffers = 2 + (is_potential ? natoms : 0);
    size_t address = 0;
    for (int c1 = 0; c1 < num_buffers; ++c1) {
        auto a1 = (c1 == 0 ? atom1 : (c1 == 1 ? atom2 : c1 - 2)); 
        for (int xyz1 = 0; xyz1 < 3; ++xyz1) {
            auto coord1 = 3 * a1 + xyz1;
            for (int c2 = c1; c2 < num_buffers; ++c2) {
                auto a2 = (c2 == 0 ? atom1 : (c2 == 1 ? atom2 : c2 - 2)); 
                auto xyz2_start = (c1 == c2 ? xyz1 : 0);
                for (int xyz2 = xyz2_start; xyz2 < 3; ++xyz2) {
                    auto coord2 = 3 * a2 + xyz2;
                    double scale = (c1 != c2 && coord1 == coord2 ? 2.0 : 1.0);
                    double val = scale * Hvalues[address];
                    Hess[coord1][coord2] += val;
                    if (!shells_equivalent) {
                        Hess[coord2][coord1] += val;
                    }
                    ++address;
                }
            }
        }
    }
}

SharedMatrix SCFDeriv::compute_hessian()
{
    // => Echo <= //

    outfile->Printf( "\n");
    outfile->Printf( "         ------------------------------------------------------------\n");
    outfile->Printf( "                                   SCF HESS                          \n");
    outfile->Printf( "                          Rob Parrish, Justin Turney,                \n");
    outfile->Printf( "                       Andy Simmonett, and Alex Sokolov              \n");
    outfile->Printf( "         ------------------------------------------------------------\n\n");

    outfile->Printf( "  ==> Geometry <==\n\n");
    molecule_->print();

    outfile->Printf( "  Nuclear repulsion = %20.15f\n", basisset_->molecule()->nuclear_repulsion_energy(dipole_field_strength_));
    outfile->Printf( "\n");

    outfile->Printf( "  ==> Basis Set <==\n\n");
    basisset_->print_by_level("outfile", print_);

    // => Registers <= //

    std::vector<std::string> hessian_terms;
    hessian_terms.push_back("Nuclear");
    hessian_terms.push_back("Kinetic");
    hessian_terms.push_back("Potential");
    hessian_terms.push_back("Overlap");
    hessian_terms.push_back("Coulomb");
    hessian_terms.push_back("Exchange");
    hessian_terms.push_back("Exchange,LR");
    hessian_terms.push_back("XC");
    hessian_terms.push_back("-D Hessian");
    hessian_terms.push_back("Response");

    // => Densities <= //
    SharedMatrix Da;
    SharedMatrix Db;
    SharedMatrix Dt;

    Da = Da_subset("AO");
    if (options_.get_str("REFERENCE") == "RHF" || options_.get_str("REFERENCE") == "RKS") {
        Db = Da;
    } else {
        Db = Db_subset("AO");
    }
    Dt = SharedMatrix(Da->clone());
    Dt->add(Db);
    Dt->set_name("Dt");

    // => Occupations (AO) <= //
    SharedMatrix Ca;
    SharedMatrix Cb;
    SharedVector eps_a;
    SharedVector eps_b;

    Ca = Ca_subset("AO", "OCC");
    eps_a = epsilon_a_subset("AO", "OCC");
    if (options_.get_str("REFERENCE") == "RHF" || options_.get_str("REFERENCE") == "RKS") {
        Cb = Ca;
        eps_b = eps_a;
    } else {
        Cb = Cb_subset("AO", "OCC");
        eps_b = epsilon_b_subset("AO", "OCC");
    }

    // => Potential/Functional <= //
    std::shared_ptr<VBase> potential;

    if (functional_->needs_xc()) {
        if (options_.get_str("REFERENCE") == "RKS") {
            potential_->set_D({Da_});
        } else {
            potential_->set_D({Da_, Db_});
        }
    }

    // => Sizings <= //
    int natom = molecule_->natom();
    int nso = basisset_->nbf();
    int nalpha = Ca->colspi()[0];
    int nbeta = Cb->colspi()[0];

    // => Nuclear Hessian <= //
    hessians_["Nuclear"] = SharedMatrix(molecule_->nuclear_repulsion_energy_deriv2().clone());
    hessians_["Nuclear"]->set_name("Nuclear Hessian");

    // => XC Hessian <= //
    timer_on("Hess: XC");
    if (functional_->needs_xc()) {
        potential_->print_header();
        hessians_["XC"] = potential_->compute_hessian();
    }
    timer_off("Hess: XC");

    // => Potential Hessian <= //
    timer_on("Hess: V");
    {
        std::vector<std::pair<double, std::array<double, 3>>> Zxyz;
        double** Dp = Dt->pointer();

        hessians_["Potential"] = SharedMatrix(hessians_["Nuclear"]->clone());
        hessians_["Potential"]->set_name("Potential Hessian");
        hessians_["Potential"]->zero();
        double** Vp = hessians_["Potential"]->pointer();

        // Potential energy derivatives
        std::shared_ptr<OneBodyAOInt> Vint(integral_->ao_potential(2));
        const double* buffer = Vint->buffer();

        for (int P = 0; P < basisset_->nshell(); P++) {
            const GaussianShell& s1 = basisset_->shell(P);
            const libint2::Shell &l2_s1 = basisset_->l2_shell(P);
            int nP = s1.nfunction();
            int oP = s1.function_index();
            int aP = s1.ncenter();
            int Px = 3 * aP + 0;
            int Py = 3 * aP + 1;
            int Pz = 3 * aP + 2;
            for (int Q = 0; Q <= P; Q++) {

                const GaussianShell& s2 = basisset_->shell(Q);
                const libint2::Shell &l2_s2 = basisset_->l2_shell(Q);
                int nQ = s2.nfunction();
                int oQ = s2.function_index();
                int aQ = s2.ncenter();

                int Qx = 3 * aQ + 0;
                int Qy = 3 * aQ + 1;
                int Qz = 3 * aQ + 2;

#define DEBUGINTS 0

#if DEBUGINTS
                outfile->Printf("AM1 %d AM2 %d a1 %f a2 %f center1 %d center2 %d\n", s1.am(), s2.am(), s1.exp(0), s2.exp(0), s1.ncenter(), s2.ncenter());
#endif
                Vint->compute_pair_deriv2(l2_s1, l2_s2);
                const auto &buffers = Vint->buffers();

                std::vector<double> Dvals;
                // find the D values against which this batch will be contracted
                for (int p = 0; p < nP; p++) {
                    for (int q = 0; q < nQ; q++) {
                        Dvals.push_back(Dp[p + oP][q + oQ]);
                    }
                }
                // build the Hessian contributions for each buffer entry
                std::vector<double> Hvals;
                for (int i = 0; i < buffers.size(); ++i) {
                    const double *buffer = buffers[i];
                    Hvals.push_back(std::inner_product(Dvals.begin(), Dvals.end(), buffer, 0.0));
                }
                process_buffers(Vp, Hvals, aP, aQ, natom, P==Q, true);
            }
        }
        // Symmetrize the result
        int dim = hessians_["Potential"]->rowdim();
        for (int row = 0; row < dim; ++row){
            for (int col = 0; col < row; ++col){
                Vp[row][col] = Vp[col][row] = (Vp[row][col] + Vp[col][row]);
            }
        }
        timer_off("Hess: V");
    }


    // => Kinetic Hessian <= //
    timer_on("Hess: T");
    {
        double** Dp = Dt->pointer();

        hessians_["Kinetic"] = SharedMatrix(hessians_["Nuclear"]->clone());
        hessians_["Kinetic"]->set_name("Kinetic Hessian");
        hessians_["Kinetic"]->zero();
        double** Tp = hessians_["Kinetic"]->pointer();

        // Kinetic energy derivatives
        std::shared_ptr<OneBodyAOInt> Tint(integral_->ao_kinetic(2));
        const double* buffer = Tint->buffer();

        for (int P = 0; P < basisset_->nshell(); P++) {
            const GaussianShell& s1 = basisset_->shell(P);
            const libint2::Shell& l2_s1 = basisset_->l2_shell(P);
            int nP = s1.nfunction();
            int oP = s1.function_index();
            int aP = s1.ncenter();
            int Px = 3 * aP + 0;
            int Py = 3 * aP + 1;
            int Pz = 3 * aP + 2;
            for (int Q = 0; Q <= P; Q++) {
                const GaussianShell& s2 = basisset_->shell(Q);
                const libint2::Shell& l2_s2 = basisset_->l2_shell(Q);
                int nQ = s2.nfunction();
                int oQ = s2.function_index();
                int aQ = s2.ncenter();

                Tint->compute_pair_deriv2(l2_s1, l2_s2);
                const auto &buffers = Tint->buffers();

                std::vector<double> Dvals;
                // find the D values against which this batch will be conctracted
                for (int p = 0; p < nP; p++) {
                    for (int q = 0; q < nQ; q++) {
                        Dvals.push_back(Dp[p + oP][q + oQ]);
                    }
                }
                // build the Hessian contributions for each buffer entry
                std::vector<double> Hvals;
                for (const double *buffer : buffers) {
                    Hvals.push_back(std::inner_product(Dvals.begin(), Dvals.end(), buffer, 0.0));
                }
                process_buffers(Tp, Hvals, aP, aQ, natom, P==Q, false);
            }
        }
        // Symmetrize the result
        int dim = hessians_["Kinetic"]->rowdim();
        for (int row = 0; row < dim; ++row){
            for (int col = 0; col < row; ++col){
                Tp[row][col] = Tp[col][row] = (Tp[row][col] + Tp[col][row]);
            }
        }
        timer_off("Hess: T");
    }

    // => Overlap Hessian <= //
    timer_on("Hess: S");
    {
        // Energy weighted density matrix
        SharedMatrix W(Da->clone());
        W->set_name("W");

        double** Wp = W->pointer();
        double** Cap = Ca->pointer();
        double** Cbp = Cb->pointer();
        double* eps_ap = eps_a->pointer();
        double* eps_bp = eps_b->pointer();

        auto* temp = new double[nso * (size_t) nalpha];

        ::memset((void*) temp, '\0', sizeof(double) * nso * nalpha);
        for (int i = 0; i < nalpha; i++) {
            C_DAXPY(nso,eps_ap[i], &Cap[0][i], nalpha, &temp[i], nalpha);
        }

        C_DGEMM('N','T',nso,nso,nalpha,-1.0,Cap[0],nalpha,temp,nalpha,0.0,Wp[0],nso);

        ::memset((void*) temp, '\0', sizeof(double) * nso * nbeta);
        for (int i = 0; i < nbeta; i++) {
            C_DAXPY(nso,eps_bp[i], &Cbp[0][i], nbeta, &temp[i], nbeta);
        }

        C_DGEMM('N','T',nso,nso,nbeta,-1.0,Cbp[0],nbeta,temp,nbeta,1.0,Wp[0],nso);

        delete[] temp;

        hessians_["Overlap"] = SharedMatrix(hessians_["Nuclear"]->clone());
        hessians_["Overlap"]->set_name("Overlap Hessian");
        hessians_["Overlap"]->zero();
        double** Sp = hessians_["Overlap"]->pointer();

        // Overlap derivatives
        std::shared_ptr<OneBodyAOInt> Sint(integral_->ao_overlap(2));
        const double* buffer = Sint->buffer();

        for (int P = 0; P < basisset_->nshell(); P++) {
            const libint2::Shell& l2_s1 = basisset_->l2_shell(P);
            const GaussianShell& s1 = basisset_->shell(P);
            int nP = s1.nfunction();
            int oP = s1.function_index();
            int aP = s1.ncenter();
            int Px = 3 * aP + 0;
            int Py = 3 * aP + 1;
            int Pz = 3 * aP + 2;
            for (int Q = 0; Q <= P; Q++) {
                const libint2::Shell& l2_s2 = basisset_->l2_shell(Q);
                const GaussianShell& s2 = basisset_->shell(Q);
                int nQ = s2.nfunction();
                int oQ = s2.function_index();
                int aQ = s2.ncenter();

                Sint->compute_pair_deriv2(l2_s1, l2_s2);
                const auto &buffers = Sint->buffers();

                std::vector<double> Wvals;
                // find the W values against which this batch will be contracted
                for (int p = 0; p < nP; p++) {
                    for (int q = 0; q < nQ; q++) {
                        Wvals.push_back(Wp[p + oP][q + oQ]);
                    }
                }
                // build the Hessian contributions for each buffer entry
                std::vector<double> Hvals;
                for (int i = 0; i < buffers.size(); ++i) {
                    const double *buffer = buffers[i];
                    Hvals.push_back(std::inner_product(Wvals.begin(), Wvals.end(), buffer, 0.0));
                }
                process_buffers(Sp, Hvals, aP, aQ, natom, P==Q, false);
            }
        }
        // Symmetrize the result
        int dim = hessians_["Overlap"]->rowdim();
        for (int row = 0; row < dim; ++row){
            for (int col = 0; col < row; ++col){
                Sp[row][col] = Sp[col][row] = (Sp[row][col] + Sp[col][row]);
            }
        }
        timer_off("Hess: S");
    }
    // => Two-Electron Hessian <= //

    timer_on("Hess: JK");

    auto jk = JKGrad::build_JKGrad(2, mintshelper_);
    jk->set_memory((size_t) (options_.get_double("SCF_MEM_SAFETY_FACTOR") * memory_ / 8L));

    jk->set_Ca(Ca);
    jk->set_Cb(Cb);
    jk->set_Da(Da);
    jk->set_Db(Db);
    jk->set_Dt(Dt);
    if (functional_) {
        jk->set_do_J(true);
        if (functional_->is_x_hybrid()) {
            jk->set_do_K(true);
        } else {
            jk->set_do_K(false);
        }
        if (functional_->is_x_lrc()) {
            jk->set_do_wK(true);
            jk->set_omega(functional_->x_omega());
        } else {
            jk->set_do_wK(false);
        }
    } else {
        jk->set_do_J(true);
        jk->set_do_K(true);
        jk->set_do_wK(false);
    }

    jk->print_header();
    jk->compute_hessian();

    std::map<std::string, SharedMatrix>& jk_hessians = jk->hessians();
    if (functional_) {
        hessians_["Coulomb"] = jk_hessians["Coulomb"];
        if (functional_->is_x_hybrid()) {
            hessians_["Exchange"] = jk_hessians["Exchange"];
            hessians_["Exchange"]->scale(-functional_->x_alpha());
        }
        if (functional_->is_x_lrc()) {
            hessians_["Exchange,LR"] = jk_hessians["Exchange,LR"];
            hessians_["Exchange,LR"]->scale(-functional_->x_beta());
        }
    } else {
        hessians_["Coulomb"] = jk_hessians["Coulomb"];
        hessians_["Exchange"] = jk_hessians["Exchange"];
        hessians_["Exchange"]->scale(-1.0);
    }
    timer_off("Hess: JK");

    // => Response Terms (Brace Yourself) <= //
#ifdef USING_BrianQC
    brianCPHFFlag = true;
#endif
    if (options_.get_str("REFERENCE") == "RHF" || 
        options_.get_str("REFERENCE") == "RKS" || 
        options_.get_str("REFERENCE") == "UHF") {
        hessians_["Response"] = hessian_response();
    } else {
        throw PSIEXCEPTION("SCFHessian: Response not implemented for this reference");
    }
#ifdef USING_BrianQC
    brianCPHFFlag = false;
#endif

    // => Total Hessian <= //
    SharedMatrix total = SharedMatrix(hessians_["Nuclear"]->clone());
    total->zero();
    for (int i = 0; i < hessian_terms.size(); i++) {
        if (hessians_.count(hessian_terms[i])) {
            total->add(hessians_[hessian_terms[i]]);
        }
    }
    // total->symmetrize_hessian(molecule_);

    hessians_["Total"] = total;
    hessians_["Total"]->set_name("Total Hessian");

    // => Final Printing <= //
    if (print_ > 1) {
        for (int i = 0; i < hessian_terms.size(); i++) {
            if (hessians_.count(hessian_terms[i])) {
                outfile->Printf("%s\n", hessian_terms[i].c_str());
                hessians_[hessian_terms[i]]->print();
            }
        }
    }
    hessians_["Total"]->print();

    return hessians_["Total"];
}

}} // Namespaces
