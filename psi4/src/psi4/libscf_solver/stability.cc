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

#include "stability.h"
#include "psi4/psi4-dec.h"
#include "psi4/libparallel/parallel.h"
#include "psi4/liboptions/liboptions.h"

#include "psi4/libpsio/psio.hpp"
#include "psi4/libciomr/libciomr.h"
#include "psi4/libscf_solver/hf.h"
#include "psi4/libfock/jk.h"
#include "psi4/libfock/solver.h"
#include "psi4/libfock/hamiltonian.h"
#include "psi4/libmints/matrix.h"
#include "psi4/physconst.h"


using namespace std;

namespace psi{

namespace scf {

PsiReturnType stability(SharedWavefunction ref_wfn, Options& options)
{

    tstart();
    std::shared_ptr<UStab> stab = std::shared_ptr<UStab>(new UStab(ref_wfn, options));
    stab->compute_energy();
    tstop();

    return Success;
}

UStab::UStab(SharedWavefunction ref_wfn, Options& options) :
       options_(options)
{
    common_init();
    set_reference(ref_wfn);
}

UStab::~UStab() {
}

void UStab::common_init()
{

    print_ = options_.get_int("PRINT");
    debug_ = options_.get_int("DEBUG");
    bench_ = options_.get_int("BENCH");
    convergence_ = options_.get_double("SOLVER_CONVERGENCE");
    memory_ = Process::environment.get_memory();

}

void UStab::set_reference(std::shared_ptr<Wavefunction> wfn)
{
    reference_wavefunction_ = wfn;

    if (!reference_wavefunction_) {
        throw PSIEXCEPTION("UStab: Run SCF first");
    }

    if (reference_wavefunction_->same_a_b_dens()) {
        throw PSIEXCEPTION("UStab: Reference is restricted!");
    }

    Ca_ = wfn->Ca();
    Cb_ = wfn->Cb();
    Cocca_  = wfn->Ca_subset("SO","OCC");
    Coccb_  = wfn->Cb_subset("SO","OCC");
    Cvira_  = wfn->Ca_subset("SO","VIR");
    Cvirb_  = wfn->Cb_subset("SO","VIR");
    eps_occa_ = wfn->epsilon_a_subset("SO","OCC");
    eps_occb_ = wfn->epsilon_b_subset("SO","OCC");
    eps_vira_ = wfn->epsilon_a_subset("SO","VIR");
    eps_virb_ = wfn->epsilon_b_subset("SO","VIR");
    molecule_ = wfn->molecule();
    basis_ = wfn->basisset();
    Eref_ = wfn->reference_energy();

}

void UStab::print_header()
{
    std::shared_ptr<Wavefunction> wfn = reference_wavefunction_;
    outfile->Printf( "\n");
    outfile->Printf( "         ------------------------------------------------------------\n");
    outfile->Printf( "                              UHF Stability code                     \n");
    outfile->Printf( "                                Jérôme Gonthier                     \n");
    outfile->Printf( "               Strong inspiration from R. Parrish's CIS              \n");
    outfile->Printf( "         ------------------------------------------------------------\n\n");

    outfile->Printf( "  ==> Geometry <==\n\n");
    molecule_->print();
    outfile->Printf( "  Nuclear repulsion = %20.15f\n", molecule_->nuclear_repulsion_energy());
    //outfile->Printf( "  Reference energy  = %20.15f\n\n", Eref_);

    outfile->Printf( "  ==> Basis Set <==\n\n");
    basis_->print_by_level("outfile", print_);

    if (debug_ > 1) {
        outfile->Printf( "  ==> Fock Matrix (MO Basis) <==\n\n");
        eps_occa_->print();
        eps_occb_->print();
        eps_vira_->print();
        eps_virb_->print();
    }
}

double UStab::compute_energy()
{
    // Main UStability Header
    print_header();

    if (!jk_)
        preiterations();

    // Construct components
    std::shared_ptr<USTABHamiltonian> H(new USTABHamiltonian(jk_, Cocca_,Cvira_,Coccb_,Cvirb_,eps_occa_,eps_vira_,
                            eps_occb_,eps_virb_));
    std::shared_ptr<DLUSolver> solver;
    if (options_.get_str("SOLVER_TYPE") == "DL")
        solver = DLUSolver::build_solver(options_,H);
    else if (options_.get_str("SOLVER_TYPE") == "RAYLEIGH")
        throw PSIEXCEPTION("Rayleigh solver not implemented for UStab.");

    // Extra Knobs
    H->set_print(print_);
    H->set_debug(debug_);
    H->set_bench(bench_);
    H->set_exact_diagonal(options_.get_bool("SOLVER_EXACT_DIAGONAL"));
    solver->set_convergence(convergence_);

    // Initialization/Memory
    solver->initialize();

    // Component Headers
    solver->print_header();
    H->print_header();
    jk_->print_header();

    solver->solve();

    // Did we converge?
    if ( !solver->converged()) {
        throw PSIEXCEPTION("Error: Roots not converged.");
    }

    // Unpack
    const std::vector<std::shared_ptr<Vector> > stabvecs = solver->eigenvectors();
    const std::vector<std::vector<double> > stabvals = solver->eigenvalues();

    std::vector<std::pair<SharedMatrix, SharedMatrix > > evec_temp;
    std::vector<std::pair<double, int> > eval_temp;

    for (size_t N = 0, index = 0; N < stabvecs.size(); ++N) {
        std::vector< std::pair<SharedMatrix, SharedMatrix> > tpair = H->unpack_paired(stabvecs[N]);
        for (int h = 0; h < Cocca_->nirrep(); h++) {
            // Spurious zero eigenvalue due to not enough states
            if (N >= (size_t)stabvecs[N]->dimpi()[h]) continue;
            evec_temp.push_back(tpair[h]);
            eval_temp.push_back(make_pair(stabvals[N][h], index));
            index++;
        }
    }

    std::sort(eval_temp.begin(), eval_temp.end());

    vecs_.clear();
    vals_.clear();

    for (size_t i = 0; i < eval_temp.size(); i++) {
        vals_.push_back(eval_temp[i].first);
        vecs_.push_back(evec_temp[eval_temp[i].second]);
    }

    if (debug_ > 1 ) {
        for (size_t i = 0; i < eval_temp.size(); ++i) {
            outfile->Printf("Eigenvalue %4i: %.12f\n", i, vals_[i]);
        }
    }


    // Finalize solver
    solver->finalize();

    return 0.0;
}

SharedMatrix UStab::analyze()
{

    // We use the convergence criterion to eliminate zero eigenvalues
    // suffering from numerical noise.
    int nirrep = vecs_[0].first->nirrep();
    int eig_dims[nirrep];
    int col_dim[nirrep];

    for (int i = 0; i < nirrep; ++i) {
        eig_dims[i] = 0;
        col_dim[i] = 1;
    }

    for (int i = 0; i < vals_.size(); ++i) {
        ++(eig_dims[vecs_[i].first->symmetry()]);

    }

    SharedMatrix eval_sym(new Matrix("SCF STABILITY EIGENVALUES", nirrep, eig_dims, col_dim));
    for (int h = 0; h < nirrep; ++h)
    {
        eig_dims[h] = 0;
    }

    for (int i = 0; i < vals_.size(); ++i) {
        int h = vecs_[i].first->symmetry();
        eval_sym->set(h,eig_dims[h],0,vals_[i]);
        ++eig_dims[h];
        if ((vals_[i] < unstable_val) && (abs(vals_[i]) > convergence_) ) {
            if ( vecs_[i].first->symmetry() == 0) {
                unstable = true;
                unstable_val = vals_[i];
                unstable_vec = vecs_[i];
            }
        }
    }


    if (unstable) {
        outfile->Printf("    Negative totally symmetric eigenvalue detected: %f \n", unstable_val);
        outfile->Printf("    Wavefunction unstable!\n");
    } else {
        outfile->Printf("    Wavefunction stable under totally symmetric rotations.\n");
        outfile->Printf("    Lowest totally symmetric eigenvalue: %f \n", vals_[0]);
    }

    return eval_sym;
}

void UStab::rotate_orbs(double step_scale)
{
    double scale = pc_pi*step_scale/2.0;
    outfile->Printf("    Rotating orbitals by %f * pi / 2 radians along unstable eigenvector.\n", step_scale);

    int nirrep = unstable_vec.first->nirrep();

    SharedMatrix unveca = unstable_vec.first;
    SharedMatrix unvecb = unstable_vec.second;
    for (int h = 0; h < nirrep; ++h) {
        int nocca = unveca->rowdim(h);
        int nvira = unveca->coldim(h);
    // Rotate the alpha orbitals
        for (int i = 0; i < nocca; ++i) {
            for (int a = nocca; a < nvira + nocca; ++a) {
                Ca_->rotate_columns(h, i, a, scale*unveca->get(h,i,a - nocca));
            }
        }
        int noccb = unvecb->rowdim(h);
        int nvirb = unvecb->coldim(h);
    // Rotate the beta orbitals
        for (int i = 0; i < noccb; ++i) {
            for (int a = noccb; a < nvirb + noccb; ++a) {
                Cb_->rotate_columns(h, i, a, scale*unvecb->get(h,i,a - noccb));
            }
        }
    }
}

void UStab::preiterations()
{
    if (!jk_) {
        if (options_.get_bool("SAVE_JK")) {
            jk_ = (static_cast<psi::scf::HF*>(reference_wavefunction_.get()))->jk();
            outfile->Printf("    Reusing JK object from SCF.\n\n");
        } else {
            if (options_.get_str("SCF_TYPE") == "DF"){
                jk_ = JK::build_JK(basis_, reference_wavefunction_->get_basisset("DF_BASIS_SCF"), options_);
            } else {
                jk_ = JK::build_JK(basis_, BasisSet::zero_ao_basis_set(), options_);
            }
            unsigned long int effective_memory = (unsigned long int)(0.125 * options_.get_double("CPHF_MEM_SAFETY_FACTOR") * memory_);
            jk_->set_memory(effective_memory);
            jk_->initialize();
        }
    }

  }

    }} // End namespaces
