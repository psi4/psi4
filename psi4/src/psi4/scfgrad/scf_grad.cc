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

#include "scf_grad.h"

#include <algorithm>
#include <numeric>
#include <sstream>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "psi4/libqt/qt.h"
#include "psi4/libpsio/psio.hpp"
#ifdef USING_ecpint
#include "psi4/libmints/ecpint.h"
#endif
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

SCFDeriv::SCFDeriv(std::shared_ptr<scf::HF> ref_wfn, Options& options) :
    Wavefunction(options)
{
    shallow_copy(ref_wfn);
    common_init();
    functional_ = ref_wfn->functional();
    potential_ = ref_wfn->V_potential();
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
    auto Da = Da_subset("AO");
    SharedMatrix Db;
    if (options_.get_str("REFERENCE") == "RHF" || options_.get_str("REFERENCE") == "RKS") {
        Db = Da;
    } else {
        Db = Db_subset("AO");
    }
    auto Dt = Da->clone();
    Dt->add(Db);
    Dt->set_name("Dt");

    // => Occupations (AO) <= //
    auto Ca_occ = Ca_subset("AO", "OCC");
    auto eps_a_occ = epsilon_a_subset("AO", "OCC");
    SharedVector eps_b_occ;
    SharedMatrix Cb_occ;
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
        const auto& shell_pairs = Vint->shellpairs();
        size_t n_pairs = shell_pairs.size();

        for (size_t p = 0; p < n_pairs; ++p) {
            auto P = shell_pairs[p].first;
            auto Q = shell_pairs[p].second;
            const GaussianShell& s1 = basisset_->shell(P);
            const GaussianShell& s2 = basisset_->shell(Q);
            int nP = s1.nfunction();
            int oP = s1.function_index();
            int aP = s1.ncenter();
            int nQ = s2.nfunction();
            int oQ = s2.function_index();
            int aQ = s2.ncenter();
            int Px = 3 * aP + 0;
            int Py = 3 * aP + 1;
            int Pz = 3 * aP + 2;
            int Qx = 3 * aQ + 0;
            int Qy = 3 * aQ + 1;
            int Qz = 3 * aQ + 2;
#define DEBUGINTS 0

#if DEBUGINTS
            outfile->Printf("AM1 %d AM2 %d a1 %f a2 %f center1 %d center2 %d\n", s1.am(), s2.am(), s1.exp(0), s2.exp(0), s1.ncenter(), s2.ncenter());
#endif
            Vint->compute_shell_deriv2(P, Q);
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
        // Symmetrize the result
        int dim = hessians_["Potential"]->rowdim();
        for (int row = 0; row < dim; ++row){
            for (int col = 0; col < row; ++col){
                Vp[row][col] = Vp[col][row] = (Vp[row][col] + Vp[col][row]);
            }
        }
        timer_off("Hess: V");
    }

    if (basisset_->has_ECP() ) {
    // => Potential Hessian <= //
#ifdef USING_ecpint
    timer_on("Hess: ECP");
    {
        double** Dp = Dt->pointer();

        hessians_["Effective Core Potential"] = SharedMatrix(hessians_["Nuclear"]->clone());
        hessians_["Effective Core Potential"]->set_name("Effective Core Potential Hessian");
        hessians_["Effective Core Potential"]->zero();
        double** ECPp = hessians_["Effective Core Potential"]->pointer();
        hessian_terms.push_back("Effective Core Potential");

        // Potential energy derivatives
        std::shared_ptr<ECPInt> ecpint(dynamic_cast<ECPInt*>(integral_->ao_ecp(2).release()));
        const auto& buffers = ecpint->buffers();

        for (int P = 0; P < basisset_->nshell(); P++) {
            const GaussianShell& s1 = basisset_->shell(P);
            int nP = s1.nfunction();
            int oP = s1.function_index();
            int aP = s1.ncenter();
            int Ax = 3 * aP + 0;
            int Ay = 3 * aP + 1;
            int Az = 3 * aP + 2;
            for (int Q = 0; Q <= P; Q++) {

                const GaussianShell& s2 = basisset_->shell(Q);
                int nQ = s2.nfunction();
                int oQ = s2.function_index();
                int aQ = s2.ncenter();

                int Bx = 3 * aQ + 0;
                int By = 3 * aQ + 1;
                int Bz = 3 * aQ + 2;

                double perm = (P == Q ? 1.0 : 2.0);

                size_t offset = static_cast<size_t> (nP)*nQ;
#define DEBUGINTS 0

#if DEBUGINTS
                outfile->Printf("AM1 %d AM2 %d a1 %f a2 %f center1 %d center2 %d\n", s1.am(), s2.am(), s1.exp(0), s2.exp(0), s1.ncenter(), s2.ncenter());
#endif
                ecpint->setup_hessian_iterations();
                while(ecpint->next_hessian_ecp()) {
                    int ecp_center = ecpint->current_ecp_center();
                    int Cx = 3 * ecp_center + 0;
                    int Cy = 3 * ecp_center + 1;
                    int Cz = 3 * ecp_center + 2;

                    ecpint->compute_shell_deriv2(P,Q);

                    const double *AxAx = buffers[ 0];
                    const double *AxAy = buffers[ 1];
                    const double *AxAz = buffers[ 2];
                    const double *AyAy = buffers[ 3];
                    const double *AyAz = buffers[ 4];
                    const double *AzAz = buffers[ 5];
                    const double *AxBx = buffers[ 6];
                    const double *AxBy = buffers[ 7];
                    const double *AxBz = buffers[ 8];
                    const double *AyBx = buffers[ 9];
                    const double *AyBy = buffers[10];
                    const double *AyBz = buffers[11];
                    const double *AzBx = buffers[12];
                    const double *AzBy = buffers[13];
                    const double *AzBz = buffers[14];
                    const double *AxCx = buffers[15];
                    const double *AxCy = buffers[16];
                    const double *AxCz = buffers[17];
                    const double *AyCx = buffers[18];
                    const double *AyCy = buffers[19];
                    const double *AyCz = buffers[20];
                    const double *AzCx = buffers[21];
                    const double *AzCy = buffers[22];
                    const double *AzCz = buffers[23];
                    const double *BxBx = buffers[24];
                    const double *BxBy = buffers[25];
                    const double *BxBz = buffers[26];
                    const double *ByBy = buffers[27];
                    const double *ByBz = buffers[28];
                    const double *BzBz = buffers[29];
                    const double *BxCx = buffers[30];
                    const double *BxCy = buffers[31];
                    const double *BxCz = buffers[32];
                    const double *ByCx = buffers[33];
                    const double *ByCy = buffers[34];
                    const double *ByCz = buffers[35];
                    const double *BzCx = buffers[36];
                    const double *BzCy = buffers[37];
                    const double *BzCz = buffers[38];
                    const double *CxCx = buffers[39];
                    const double *CxCy = buffers[40];
                    const double *CxCz = buffers[41];
                    const double *CyCy = buffers[42];
                    const double *CyCz = buffers[43];
                    const double *CzCz = buffers[44];

                    double ABscale = (aP == aQ ? 2.0 : 1.0);
                    double ACscale = (aP == ecp_center ? 2.0 : 1.0);
                    double BCscale = (aQ == ecp_center ? 2.0 : 1.0);

                    for (int p = 0; p < nP; p++) {
                        for (int q = 0; q < nQ; q++) {
                            double Delem = perm * Dp[p + oP][q + oQ];
                            double tmpAxAx =  Delem * (*AxAx);
                            double tmpAxAy =  Delem * (*AxAy);
                            double tmpAxAz =  Delem * (*AxAz);
                            double tmpAyAy =  Delem * (*AyAy);
                            double tmpAyAz =  Delem * (*AyAz);
                            double tmpAzAz =  Delem * (*AzAz);
                            double tmpAxBx =  Delem * (*AxBx);
                            double tmpAxBy =  Delem * (*AxBy);
                            double tmpAxBz =  Delem * (*AxBz);
                            double tmpAyBx =  Delem * (*AyBx);
                            double tmpAyBy =  Delem * (*AyBy);
                            double tmpAyBz =  Delem * (*AyBz);
                            double tmpAzBx =  Delem * (*AzBx);
                            double tmpAzBy =  Delem * (*AzBy);
                            double tmpAzBz =  Delem * (*AzBz);
                            double tmpAxCx =  Delem * (*AxCx);
                            double tmpAxCy =  Delem * (*AxCy);
                            double tmpAxCz =  Delem * (*AxCz);
                            double tmpAyCx =  Delem * (*AyCx);
                            double tmpAyCy =  Delem * (*AyCy);
                            double tmpAyCz =  Delem * (*AyCz);
                            double tmpAzCx =  Delem * (*AzCx);
                            double tmpAzCy =  Delem * (*AzCy);
                            double tmpAzCz =  Delem * (*AzCz);
                            double tmpBxBx =  Delem * (*BxBx);
                            double tmpBxBy =  Delem * (*BxBy);
                            double tmpBxBz =  Delem * (*BxBz);
                            double tmpByBy =  Delem * (*ByBy);
                            double tmpByBz =  Delem * (*ByBz);
                            double tmpBzBz =  Delem * (*BzBz);
                            double tmpBxCx =  Delem * (*BxCx);
                            double tmpBxCy =  Delem * (*BxCy);
                            double tmpBxCz =  Delem * (*BxCz);
                            double tmpByCx =  Delem * (*ByCx);
                            double tmpByCy =  Delem * (*ByCy);
                            double tmpByCz =  Delem * (*ByCz);
                            double tmpBzCx =  Delem * (*BzCx);
                            double tmpBzCy =  Delem * (*BzCy);
                            double tmpBzCz =  Delem * (*BzCz);
                            double tmpCxCx =  Delem * (*CxCx);
                            double tmpCxCy =  Delem * (*CxCy);
                            double tmpCxCz =  Delem * (*CxCz);
                            double tmpCyCy =  Delem * (*CyCy);
                            double tmpCyCz =  Delem * (*CyCz);
                            double tmpCzCz =  Delem * (*CzCz);

                            // AxAx
                            ECPp[Ax][Ax] += tmpAxAx;
                            // AyAy
                            ECPp[Ay][Ay] += tmpAyAy;
                            // AzAz
                            ECPp[Az][Az] += tmpAzAz;
                            // AxAy
                            ECPp[Ax][Ay] += tmpAxAy;
                            // AxAz
                            ECPp[Ax][Az] += tmpAxAz;
                            // AyAz
                            ECPp[Ay][Az] += tmpAyAz;
                            // BxBx
                            ECPp[Bx][Bx] += tmpBxBx;
                            // ByBy
                            ECPp[By][By] += tmpByBy;
                            // BzBz
                            ECPp[Bz][Bz] += tmpBzBz;
                            // BxBy
                            ECPp[Bx][By] += tmpBxBy;
                            // BxBz
                            ECPp[Bx][Bz] += tmpBxBz;
                            // ByBz
                            ECPp[By][Bz] += tmpByBz;
                            // AxBx
                            ECPp[Ax][Bx] += ABscale* tmpAxBx;
                            // AxBy
                            ECPp[Ax][By] += tmpAxBy;
                            // AxBz
                            ECPp[Ax][Bz] += tmpAxBz;
                            // AyBx
                            ECPp[Ay][Bx] += tmpAyBx;
                            // AyBy
                            ECPp[Ay][By] += ABscale*tmpAyBy;
                            // AyBz
                            ECPp[Ay][Bz] += tmpAyBz;
                            // AzBx
                            ECPp[Az][Bx] += tmpAzBx;
                            // AzBy
                            ECPp[Az][By] += tmpAzBy;
                            // AzBz
                            ECPp[Az][Bz] += ABscale*tmpAzBz;
                            // CxAx
                            ECPp[Cx][Ax] += ACscale*tmpAxCx;
                            // CxAy
                            ECPp[Cx][Ay] += tmpAyCx;
                            // CxAz
                            ECPp[Cx][Az] += tmpAzCx;
                            // CyAx
                            ECPp[Cy][Ax] += tmpAxCy;
                            // CyAy
                            ECPp[Cy][Ay] += ACscale*tmpAyCy;
                            // CyAz
                            ECPp[Cy][Az] += tmpAzCy;
                            // CzAx
                            ECPp[Cz][Ax] += tmpAxCz;
                            // CzAy
                            ECPp[Cz][Ay] += tmpAyCz;
                            // CzAz
                            ECPp[Cz][Az] += ACscale*tmpAzCz;
                            // CxBx
                            ECPp[Cx][Bx] += BCscale*tmpBxCx;
                            // CxBy
                            ECPp[Cx][By] += tmpByCx;
                            // CxBz
                            ECPp[Cx][Bz] += tmpBzCx;
                            // CyBx
                            ECPp[Cy][Bx] += tmpBxCy;
                            // CyBy
                            ECPp[Cy][By] += BCscale*tmpByCy;
                            // CyBz
                            ECPp[Cy][Bz] += tmpBzCy;
                            // CzBx
                            ECPp[Cz][Bx] += tmpBxCz;
                            // CzBy
                            ECPp[Cz][By] += tmpByCz;
                            // CzBz
                            ECPp[Cz][Bz] += BCscale*tmpBzCz;
                            // CxCx
                            ECPp[Cx][Cx] += tmpCxCx;
                            // CyCy
                            ECPp[Cy][Cy] += tmpCyCy;
                            // CzCz
                            ECPp[Cz][Cz] += tmpCzCz;
                            // CxCy
                            ECPp[Cx][Cy] += tmpCxCy;
                            // CxCz
                            ECPp[Cx][Cz] += tmpCxCz;
                            // CyCz
                            ECPp[Cy][Cz] += tmpCyCz;

                            ++AxAx;
                            ++AxAy;
                            ++AxAz;
                            ++AyAy;
                            ++AyAz;
                            ++AzAz;
                            ++AxBx;
                            ++AxBy;
                            ++AxBz;
                            ++AyBx;
                            ++AyBy;
                            ++AyBz;
                            ++AzBx;
                            ++AzBy;
                            ++AzBz;
                            ++AxCx;
                            ++AxCy;
                            ++AxCz;
                            ++AyCx;
                            ++AyCy;
                            ++AyCz;
                            ++AzCx;
                            ++AzCy;
                            ++AzCz;
                            ++BxBx;
                            ++BxBy;
                            ++BxBz;
                            ++ByBy;
                            ++ByBz;
                            ++BzBz;
                            ++BxCx;
                            ++BxCy;
                            ++BxCz;
                            ++ByCx;
                            ++ByCy;
                            ++ByCz;
                            ++BzCx;
                            ++BzCy;
                            ++BzCz;
                            ++CxCx;
                            ++CxCy;
                            ++CxCz;
                            ++CyCy;
                            ++CyCz;
                            ++CzCz;
                        }
                    }
                }
            }
        }
        // Symmetrize the result
        int dim = hessians_["Effective Core Potential"]->rowdim();
        for (int row = 0; row < dim; ++row){
            for (int col = 0; col < row; ++col){
                ECPp[row][col] = ECPp[col][row] = (ECPp[row][col] + ECPp[col][row]);
            }
        }
        timer_off("Hess: ECP");
    }
#endif
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

        const auto& shell_pairs = Tint->shellpairs();
        size_t n_pairs = shell_pairs.size();

        for (size_t p = 0; p < n_pairs; ++p) {
            auto P = shell_pairs[p].first;
            auto Q = shell_pairs[p].second;
            const GaussianShell& s1 = basisset_->shell(P);
            const GaussianShell& s2 = basisset_->shell(Q);
            int nP = s1.nfunction();
            int oP = s1.function_index();
            int aP = s1.ncenter();
            int Px = 3 * aP + 0;
            int Py = 3 * aP + 1;
            int Pz = 3 * aP + 2;
            int nQ = s2.nfunction();
            int oQ = s2.function_index();
            int aQ = s2.ncenter();

            Tint->compute_shell_deriv2(P, Q);
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

        const auto& shell_pairs = Sint->shellpairs();
        size_t n_pairs = shell_pairs.size();

        for (size_t p = 0; p < n_pairs; ++p) {
            auto P = shell_pairs[p].first;
            auto Q = shell_pairs[p].second;
            const GaussianShell& s1 = basisset_->shell(P);
            const GaussianShell& s2 = basisset_->shell(Q);
            int nP = s1.nfunction();
            int oP = s1.function_index();
            int aP = s1.ncenter();
            int nQ = s2.nfunction();
            int oQ = s2.function_index();
            int aQ = s2.ncenter();
            int Px = 3 * aP + 0;
            int Py = 3 * aP + 1;
            int Pz = 3 * aP + 2;

            Sint->compute_shell_deriv2(P, Q);
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
