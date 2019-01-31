/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2019 The Psi4 Developers.
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
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libpsi4util/process.h"

#include "jk_grad.h"

namespace psi {
namespace scfgrad {

SCFGrad::SCFGrad(SharedWavefunction ref_wfn, Options& options) :
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
SCFGrad::~SCFGrad()
{
}
void SCFGrad::common_init()
{


    print_ = options_.get_int("PRINT");
    debug_ = options_.get_int("DEBUG");
}
SharedMatrix SCFGrad::compute_gradient()
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

    std::shared_ptr<JKGrad> jk = JKGrad::build_JKGrad(1, basisset_, basissets_["DF_BASIS_SCF"]);
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

    std::map<std::string, SharedMatrix>& jk_gradients = jk->gradients();
    gradients_["Coulomb"] = jk_gradients["Coulomb"];
    if (functional_->is_x_hybrid()) {
        gradients_["Exchange"] = jk_gradients["Exchange"];
        gradients_["Exchange"]->scale(-functional_->x_alpha());
    }
    if (functional_->is_x_lrc()) {
        gradients_["Exchange,LR"] = jk_gradients["Exchange,LR"];
        gradients_["Exchange,LR"]->scale(-functional_->x_beta());
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
SharedMatrix SCFGrad::compute_hessian()
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
    hessian_terms.push_back("Total");

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
    std::shared_ptr<SuperFunctional> functional;
    std::shared_ptr<VBase> potential;

    if (functional_->needs_xc()) {
        throw PSIEXCEPTION("Missing XC derivatives for Hessians");
        // if (options_.get_str("REFERENCE") == "RKS") {
        //     potential_->set_D({Da_});
        // } else {
        //     potential_->set_D({Da_, Db_});
        // }
    }

    // => Sizings <= //
    int natom = molecule_->natom();
    int nso = basisset_->nbf();
    int nalpha = Ca->colspi()[0];
    int nbeta = Cb->colspi()[0];

    // => Nuclear Hessian <= //
    hessians_["Nuclear"] = SharedMatrix(molecule_->nuclear_repulsion_energy_deriv2().clone());
    hessians_["Nuclear"]->set_name("Nuclear Hessian");

    auto Zxyz = std::make_shared<Matrix>("Zxyz", 1, 4);

    // => Potential Hessian <= //
    timer_on("Hess: V");
    {
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
            int nP = s1.nfunction();
            int oP = s1.function_index();
            int aP = s1.ncenter();
            int Px = 3 * aP + 0;
            int Py = 3 * aP + 1;
            int Pz = 3 * aP + 2;
            for (int Q = 0; Q <= P; Q++) {

                const GaussianShell& s2 = basisset_->shell(Q);
                int nQ = s2.nfunction();
                int oQ = s2.function_index();
                int aQ = s2.ncenter();

                int Qx = 3 * aQ + 0;
                int Qy = 3 * aQ + 1;
                int Qz = 3 * aQ + 2;

                double perm = (P == Q ? 1.0 : 2.0);

                size_t offset = nP*nQ;
#define DEBUGINTS 0

#if DEBUGINTS
                outfile->Printf("AM1 %d AM2 %d a1 %f a2 %f center1 %d center2 %d\n", s1.am(), s2.am(), s1.exp(0), s2.exp(0), s1.ncenter(), s2.ncenter());
#endif
                for(int atom = 0; atom < natom; ++atom){
                    int Cx = 3 * atom + 0;
                    int Cy = 3 * atom + 1;
                    int Cz = 3 * atom + 2;

                    double Z = molecule_->Z(atom);
                    Vint->set_origin(molecule_->xyz(atom));

                    Vint->compute_shell_deriv2(P,Q);

                    const double *CxAx = buffer +  0*offset;
                    const double *CxAy = buffer +  1*offset;
                    const double *CxAz = buffer +  2*offset;
                    const double *CyAx = buffer +  3*offset;
                    const double *CyAy = buffer +  4*offset;
                    const double *CyAz = buffer +  5*offset;
                    const double *CzAx = buffer +  6*offset;
                    const double *CzAy = buffer +  7*offset;
                    const double *CzAz = buffer +  8*offset;
                    const double *AxAx = buffer +  9*offset;
                    const double *AxAy = buffer + 10*offset;
                    const double *AxAz = buffer + 11*offset;
                    const double *AyAy = buffer + 12*offset;
                    const double *AyAz = buffer + 13*offset;
                    const double *AzAz = buffer + 14*offset;
                    const double *BxBx = buffer + 15*offset;
                    const double *BxBy = buffer + 16*offset;
                    const double *BxBz = buffer + 17*offset;
                    const double *ByBy = buffer + 18*offset;
                    const double *ByBz = buffer + 19*offset;
                    const double *BzBz = buffer + 20*offset;
                    const double *CxCx = buffer + 21*offset;
                    const double *CxCy = buffer + 22*offset;
                    const double *CxCz = buffer + 23*offset;
                    const double *CyCy = buffer + 24*offset;
                    const double *CyCz = buffer + 25*offset;
                    const double *CzCz = buffer + 26*offset;

                    double ABscale = (aP == aQ ? 2.0 : 1.0);
                    double ACscale = (aP == atom ? 2.0 : 1.0);
                    double BCscale = (aQ == atom ? 2.0 : 1.0);

                    for (int p = 0; p < nP; p++) {
                        for (int q = 0; q < nQ; q++) {
                            double Delem = perm * Z * Dp[p + oP][q + oQ];
                            double tmpCxAx = Delem * (*CxAx);
                            double tmpCxAy = Delem * (*CxAy);
                            double tmpCxAz = Delem * (*CxAz);
                            double tmpCyAx = Delem * (*CyAx);
                            double tmpCyAy = Delem * (*CyAy);
                            double tmpCyAz = Delem * (*CyAz);
                            double tmpCzAx = Delem * (*CzAx);
                            double tmpCzAy = Delem * (*CzAy);
                            double tmpCzAz = Delem * (*CzAz);
                            double tmpAxAx = Delem * (*AxAx);
                            double tmpAxAy = Delem * (*AxAy);
                            double tmpAxAz = Delem * (*AxAz);
                            double tmpAyAy = Delem * (*AyAy);
                            double tmpAyAz = Delem * (*AyAz);
                            double tmpAzAz = Delem * (*AzAz);
                            double tmpBxBx = Delem * (*BxBx);
                            double tmpBxBy = Delem * (*BxBy);
                            double tmpBxBz = Delem * (*BxBz);
                            double tmpByBy = Delem * (*ByBy);
                            double tmpByBz = Delem * (*ByBz);
                            double tmpBzBz = Delem * (*BzBz);
                            double tmpCxCx = Delem * (*CxCx);
                            double tmpCxCy = Delem * (*CxCy);
                            double tmpCxCz = Delem * (*CxCz);
                            double tmpCyCy = Delem * (*CyCy);
                            double tmpCyCz = Delem * (*CyCz);
                            double tmpCzCz = Delem * (*CzCz);

                            /*
                             * Translational invariance relationship for derivatives w.r.t. centers A, B and C:
                             *
                             *     ∂ S   ∂ S   ∂ S
                             *     --- + --- + ---  =  0
                             *     ∂ A   ∂ B   ∂ C
                             *
                             * Take the derivative again, w.r.t. A, B and C to get relationships like
                             *
                             *     ∂^2 S   ∂^2 S   ∂^2 S
                             *     ----- + ----- + -----  =  0
                             *     ∂A ∂A   ∂B ∂A   ∂C ∂A
                             *
                             * which leads to the identities
                             *
                             *     ∂^2 S     ∂^2 S     ∂^2 S     ∂^2 S
                             *     -----  =  -----  +  -----  -  -----
                             *     ∂A ∂B     ∂C ∂C     ∂C ∂A     ∂B ∂B
                             *
                             * and
                             *
                             *     ∂^2 S     ∂^2 S     ∂^2 S     ∂^2 S
                             *     -----  =  -----  +  -----  -  -----
                             *     ∂B ∂C     ∂A ∂C     ∂A ∂A     ∂B ∂B
                             *
                             * Currently we compute all of the following
                             *
                             *     ∂^2 S     ∂^2 S     ∂^2 S     ∂^2 S
                             *     -----     -----     -----     -----
                             *     ∂A ∂C     ∂A ∂A     ∂B ∂B     ∂C ∂C
                             *
                             * and use the identities above to fill in the gaps
                             *
                             */
                            // AxAx
                            Vp[Px][Px] += tmpAxAx;
                            // AyAy
                            Vp[Py][Py] += tmpAyAy;
                            // AzAz
                            Vp[Pz][Pz] += tmpAzAz;
                            // AxAy
                            Vp[Px][Py] += tmpAxAy;
                            // AxAz
                            Vp[Px][Pz] += tmpAxAz;
                            // AyAz
                            Vp[Py][Pz] += tmpAyAz;
                            // BxBx
                            Vp[Qx][Qx] += tmpBxBx;
                            // ByBy
                            Vp[Qy][Qy] += tmpByBy;
                            // BzBz
                            Vp[Qz][Qz] += tmpBzBz;
                            // BxBy
                            Vp[Qx][Qy] += tmpBxBy;
                            // BxBz
                            Vp[Qx][Qz] += tmpBxBz;
                            // ByBz
                            Vp[Qy][Qz] += tmpByBz;
                            // AxBx
                            Vp[Px][Qx] += ABscale*(tmpCxCx + tmpCxAx - tmpBxBx);
                            // AxBy
                            Vp[Px][Qy] += tmpCxCy + tmpCxAy - tmpBxBy;
                            // AxBz
                            Vp[Px][Qz] += tmpCxCz + tmpCxAz - tmpBxBz;
                            // AyBx
                            Vp[Py][Qx] += tmpCxCy + tmpCyAx - tmpBxBy;
                            // AyBy
                            Vp[Py][Qy] += ABscale*(tmpCyCy + tmpCyAy - tmpByBy);
                            // AyBz
                            Vp[Py][Qz] += tmpCyCz + tmpCyAz - tmpByBz;
                            // AzBx
                            Vp[Pz][Qx] += tmpCxCz + tmpCzAx - tmpBxBz;
                            // AzBy
                            Vp[Pz][Qy] += tmpCyCz + tmpCzAy - tmpByBz;
                            // AzBz
                            Vp[Pz][Qz] += ABscale*(tmpCzCz + tmpCzAz - tmpBzBz);
                            // CxAx
                            Vp[Cx][Px] += ACscale*tmpCxAx;
                            // CxAy
                            Vp[Cx][Py] += tmpCxAy;
                            // CxAz
                            Vp[Cx][Pz] += tmpCxAz;
                            // CyAx
                            Vp[Cy][Px] += tmpCyAx;
                            // CyAy
                            Vp[Cy][Py] += ACscale*tmpCyAy;
                            // CyAz
                            Vp[Cy][Pz] += tmpCyAz;
                            // CzAx
                            Vp[Cz][Px] += tmpCzAx;
                            // CzAy
                            Vp[Cz][Py] += tmpCzAy;
                            // CzAz
                            Vp[Cz][Pz] += ACscale*tmpCzAz;
                            // CxBx
                            Vp[Cx][Qx] += BCscale*(tmpCxAx + tmpAxAx - tmpBxBx);
                            // CxBy
                            Vp[Cx][Qy] += tmpCyAx + tmpAxAy - tmpBxBy;
                            // CxBz
                            Vp[Cx][Qz] += tmpCzAx + tmpAxAz - tmpBxBz;
                            // CyBx
                            Vp[Cy][Qx] += tmpCxAy + tmpAxAy - tmpBxBy;
                            // CyBy
                            Vp[Cy][Qy] += BCscale*(tmpCyAy + tmpAyAy - tmpByBy);
                            // CyBz
                            Vp[Cy][Qz] += tmpCzAy + tmpAyAz - tmpByBz;
                            // CzBx
                            Vp[Cz][Qx] += tmpCxAz + tmpAxAz - tmpBxBz;
                            // CzBy
                            Vp[Cz][Qy] += tmpCyAz + tmpAyAz - tmpByBz;
                            // CzBz
                            Vp[Cz][Qz] += BCscale*(tmpCzAz + tmpAzAz - tmpBzBz);
                            // CxCx
                            Vp[Cx][Cx] += tmpCxCx;
                            // CyCy
                            Vp[Cy][Cy] += tmpCyCy;
                            // CzCz
                            Vp[Cz][Cz] += tmpCzCz;
                            // CxCy
                            Vp[Cx][Cy] += tmpCxCy;
                            // CxCz
                            Vp[Cx][Cz] += tmpCxCz;
                            // CyCz
                            Vp[Cy][Cz] += tmpCyCz;

                            ++CxAx;
                            ++CxAy;
                            ++CxAz;
                            ++CyAx;
                            ++CyAy;
                            ++CyAz;
                            ++CzAx;
                            ++CzAy;
                            ++CzAz;
                            ++AxAx;
                            ++AxAy;
                            ++AxAz;
                            ++AyAy;
                            ++AyAz;
                            ++AzAz;
                            ++BxBx;
                            ++BxBy;
                            ++BxBz;
                            ++ByBy;
                            ++ByBz;
                            ++BzBz;
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
            int nP = s1.nfunction();
            int oP = s1.function_index();
            int aP = s1.ncenter();
            int Px = 3 * aP + 0;
            int Py = 3 * aP + 1;
            int Pz = 3 * aP + 2;
            for (int Q = 0; Q <= P; Q++) {

                Tint->compute_shell_deriv2(P,Q);

                const GaussianShell& s2 = basisset_->shell(Q);
                int nQ = s2.nfunction();
                int oQ = s2.function_index();
                int aQ = s2.ncenter();

                int Qx = 3 * aQ + 0;
                int Qy = 3 * aQ + 1;
                int Qz = 3 * aQ + 2;

                size_t offset = nP*nQ;

                double perm = (P == Q ? 1.0 : 2.0);

                const double *pxx = buffer + 0*offset;
                const double *pxy = buffer + 1*offset;
                const double *pxz = buffer + 2*offset;
                const double *pyy = buffer + 3*offset;
                const double *pyz = buffer + 4*offset;
                const double *pzz = buffer + 5*offset;

                double diagscale = (aP == aQ ? 2.0 : 1.0);

                for (int p = 0; p < nP; p++) {
                    for (int q = 0; q < nQ; q++) {
                        double Delem = perm * Dp[p + oP][q + oQ];
                        double tmpxx = Delem * (*pxx);
                        double tmpxy = Delem * (*pxy);
                        double tmpxz = Delem * (*pxz);
                        double tmpyy = Delem * (*pyy);
                        double tmpyz = Delem * (*pyz);
                        double tmpzz = Delem * (*pzz);

                        /*
                         * Translational invariance relationship for derivatives w.r.t. centers A and B:
                         *
                         *     ∂ S   ∂ S
                         *     --- + ---  =  0
                         *     ∂ A   ∂ B
                         *
                         * Take the derivative again, w.r.t. A and B to get
                         *
                         *     ∂^2 S   ∂^2 S
                         *     ----- + -----  =  0
                         *     ∂A ∂A   ∂B ∂A
                         *
                         *     ∂^2 S   ∂^2 S
                         *     ----- + -----  =  0
                         *     ∂A ∂A   ∂B ∂A
                         *
                         *  Therefore we have
                         *
                         *     ∂^2 S   ∂^2 S       ∂^2 S
                         *     ----- = -----  =  - -----
                         *     ∂A ∂A   ∂B ∂B       ∂A ∂B
                         *
                         *  Only the double derivative w.r.t. A is provided, so we need to fill in the blanks below.
                         */

                        // AxAx
                        Tp[Px][Px] += tmpxx;
                        // AxAy
                        Tp[Px][Py] += tmpxy;
                        // AxAz
                        Tp[Px][Pz] += tmpxz;
                        // AyAy
                        Tp[Py][Py] += tmpyy;
                        // AyAz
                        Tp[Py][Pz] += tmpyz;
                        // AzAz
                        Tp[Pz][Pz] += tmpzz;
                        // BxBx
                        Tp[Qx][Qx] += tmpxx;
                        // BxBy
                        Tp[Qx][Qy] += tmpxy;
                        // BxBz
                        Tp[Qx][Qz] += tmpxz;
                        // ByBy
                        Tp[Qy][Qy] += tmpyy;
                        // ByBz
                        Tp[Qy][Qz] += tmpyz;
                        // BzBz
                        Tp[Qz][Qz] += tmpzz;
                        // AxBx
                        Tp[Px][Qx] += -diagscale*tmpxx;
                        // AxBy
                        Tp[Px][Qy] += -tmpxy;
                        // AxBz
                        Tp[Px][Qz] += -tmpxz;
                        // AyBx
                        Tp[Py][Qx] += -tmpxy;
                        // AyBy
                        Tp[Py][Qy] += -diagscale*tmpyy;
                        // AyBz
                        Tp[Py][Qz] += -tmpyz;
                        // AzBx
                        Tp[Pz][Qx] += -tmpxz;
                        // AzBy
                        Tp[Pz][Qy] += -tmpyz;
                        // AzBz
                        Tp[Pz][Qz] += -diagscale*tmpzz;

                        ++pxx;
                        ++pxy;
                        ++pxz;
                        ++pyy;
                        ++pyz;
                        ++pzz;
                    }
                }
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
            const GaussianShell& s1 = basisset_->shell(P);
            int nP = s1.nfunction();
            int oP = s1.function_index();
            int aP = s1.ncenter();
            int Px = 3 * aP + 0;
            int Py = 3 * aP + 1;
            int Pz = 3 * aP + 2;
            for (int Q = 0; Q <= P; Q++) {

                Sint->compute_shell_deriv2(P,Q);

                const GaussianShell& s2 = basisset_->shell(Q);
                int nQ = s2.nfunction();
                int oQ = s2.function_index();
                int aQ = s2.ncenter();

                int Qx = 3 * aQ + 0;
                int Qy = 3 * aQ + 1;
                int Qz = 3 * aQ + 2;

                size_t offset = nP*nQ;

                double perm = (P == Q ? 1.0 : 2.0);

                const double *pxx = buffer + 0*offset;
                const double *pxy = buffer + 1*offset;
                const double *pxz = buffer + 2*offset;
                const double *pyy = buffer + 3*offset;
                const double *pyz = buffer + 4*offset;
                const double *pzz = buffer + 5*offset;

                double diagscale = (aP == aQ ? 2.0 : 1.0);

                for (int p = 0; p < nP; p++) {
                    for (int q = 0; q < nQ; q++) {
                        double Welem = perm * Wp[p + oP][q + oQ];
                        double tmpxx = Welem * (*pxx);
                        double tmpxy = Welem * (*pxy);
                        double tmpxz = Welem * (*pxz);
                        double tmpyy = Welem * (*pyy);
                        double tmpyz = Welem * (*pyz);
                        double tmpzz = Welem * (*pzz);

                        /*
                         * Translational invariance relationship for derivatives w.r.t. centers A and B:
                         *
                         *     ∂ S   ∂ S
                         *     --- + ---  =  0
                         *     ∂ A   ∂ B
                         *
                         * Take the derivative again, w.r.t. A and B to get
                         *
                         *     ∂^2 S   ∂^2 S
                         *     ----- + -----  =  0
                         *     ∂A ∂A   ∂B ∂A
                         *
                         *     ∂^2 S   ∂^2 S
                         *     ----- + -----  =  0
                         *     ∂A ∂A   ∂B ∂A
                         *
                         *  Therefore we have
                         *
                         *     ∂^2 S   ∂^2 S       ∂^2 S
                         *     ----- = -----  =  - -----
                         *     ∂A ∂A   ∂B ∂B       ∂A ∂B
                         *
                         *  Only the double derivative w.r.t. A is provided, so we need to fill in the blanks below.
                         */

                        // AxAx
                        Sp[Px][Px] += tmpxx;
                        // AxAy
                        Sp[Px][Py] += tmpxy;
                        // AxAz
                        Sp[Px][Pz] += tmpxz;
                        // AyAy
                        Sp[Py][Py] += tmpyy;
                        // AyAz
                        Sp[Py][Pz] += tmpyz;
                        // AzAz
                        Sp[Pz][Pz] += tmpzz;
                        // BxBx
                        Sp[Qx][Qx] += tmpxx;
                        // BxBy
                        Sp[Qx][Qy] += tmpxy;
                        // BxBz
                        Sp[Qx][Qz] += tmpxz;
                        // ByBy
                        Sp[Qy][Qy] += tmpyy;
                        // ByBz
                        Sp[Qy][Qz] += tmpyz;
                        // BzBz
                        Sp[Qz][Qz] += tmpzz;
                        // AxBx
                        Sp[Px][Qx] += -diagscale*tmpxx;
                        // AxBy
                        Sp[Px][Qy] += -tmpxy;
                        // AxBz
                        Sp[Px][Qz] += -tmpxz;
                        // AyBx
                        Sp[Py][Qx] += -tmpxy;
                        // AyBy
                        Sp[Py][Qy] += -diagscale*tmpyy;
                        // AyBz
                        Sp[Py][Qz] += -tmpyz;
                        // AzBx
                        Sp[Pz][Qx] += -tmpxz;
                        // AzBy
                        Sp[Pz][Qy] += -tmpyz;
                        // AzBz
                        Sp[Pz][Qz] += -diagscale*tmpzz;

                        ++pxx;
                        ++pxy;
                        ++pxz;
                        ++pyy;
                        ++pyz;
                        ++pzz;
                    }
                }
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

    std::shared_ptr<JKGrad> jk = JKGrad::build_JKGrad(2, basisset_, basissets_["DF_BASIS_SCF"]);
    jk->set_memory((size_t) (options_.get_double("SCF_MEM_SAFETY_FACTOR") * memory_ / 8L));

    jk->set_Ca(Ca);
    jk->set_Cb(Cb);
    jk->set_Da(Da);
    jk->set_Db(Db);
    jk->set_Dt(Dt);
    if (functional) {
        jk->set_do_J(true);
        if (functional->is_x_hybrid()) {
            jk->set_do_K(true);
        } else {
            jk->set_do_K(false);
        }
        if (functional->is_x_lrc()) {
            jk->set_do_wK(true);
            jk->set_omega(functional->x_omega());
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
    if (functional) {
        hessians_["Coulomb"] = jk_hessians["Coulomb"];
        if (functional->is_x_hybrid()) {
            hessians_["Exchange"] = jk_hessians["Exchange"];
            hessians_["Exchange"]->scale(-functional->x_alpha());
        }
        if (functional->is_x_lrc()) {
            hessians_["Exchange,LR"] = jk_hessians["Exchange,LR"];
            hessians_["Exchange,LR"]->scale(-functional->x_beta());
        }
    } else {
        hessians_["Coulomb"] = jk_hessians["Coulomb"];
        hessians_["Exchange"] = jk_hessians["Exchange"];
        hessians_["Exchange"]->scale(-1.0);
    }
    timer_off("Hess: JK");

    // => XC Hessian <= //
    timer_on("Hess: XC");
    if (functional) {
        potential->print_header();
        throw PSIEXCEPTION("KS Hessians not implemented");
        //hessians_["XC"] = potential->compute_hessian();
    }
    timer_off("Hess: XC");

    // => Response Terms (Brace Yourself) <= //
    if (options_.get_str("REFERENCE") == "RHF") {
        hessians_["Response"] = rhf_hessian_response();
    } else {
        throw PSIEXCEPTION("SCFHessian: Response not implemented for this reference");
    }

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
                printf("%s\n", hessian_terms[i].c_str());
                hessians_[hessian_terms[i]]->print();
            }
        }
    }
    // hessians_["Total"]->print();

    return hessians_["Total"];
}

}} // Namespaces
