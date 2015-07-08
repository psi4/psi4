/*
 *@BEGIN LICENSE
 *
 * stability by Psi4 Developer, a plugin to:
 *
 * PSI4: an ab initio quantum chemistry software package
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
 *@END LICENSE
 */

#include "stability.h"
#include <libplugin/plugin.h>
#include <psi4-dec.h>
#include <libparallel/parallel.h>
#include <liboptions/liboptions.h>
#include <libmints/mints.h>
#include <libpsio/psio.hpp>
#include <libciomr/libciomr.h>
#include <libscf_solver/hf.h>
#include <libfock/jk.h>
#include <libfock/solver.h>
#include <libfock/hamiltonian.h>
#include <physconst.h>

using namespace boost;
using namespace std;

namespace psi{ 

namespace scf {

//extern "C"
//int read_options(std::string name, Options& options)
//{
//    if (name == "STABILITY"|| options.read_globals()) {
//        /*- The amount of information printed to the output file -*/
//        options.add_int("PRINT", 1);
//    }
//
//    return true;
//}

PsiReturnType stability(Options& options)
{

    tstart();
    boost::shared_ptr<UStab> stab = boost::shared_ptr<UStab>(new UStab());
    stab->compute_energy();
    tstop();

    return Success;
}

UStab::UStab() : options_(Process::environment.options) {
    common_init();
}

UStab::~UStab() {
}

void UStab::common_init()
{
    boost::shared_ptr<Wavefunction> ref = Process::environment.wavefunction();
    if (!ref) {
        throw PSIEXCEPTION("Need an SCF wavefunction in Process::environment !!");
    }

    set_reference(ref);

    print_ = options_.get_int("PRINT");
    debug_ = options_.get_int("DEBUG");
    bench_ = options_.get_int("BENCH");
    convergence_ = options_.get_double("SOLVER_CONVERGENCE");
    memory_ = options_.get_double("MEMORY");

}

void UStab::set_reference(boost::shared_ptr<Wavefunction> wfn)
{
    reference_wavefunction_ = wfn;

    if (!reference_wavefunction_) {
        throw PSIEXCEPTION("UStab: Run SCF first");
    }

    if (reference_wavefunction_->same_a_b_dens()) {
        throw PSIEXCEPTION("UStab: Reference is restricted!");
    }

    Ca_ = wfn->Ca_subset("SO","ALL");
    Cb_ = wfn->Cb_subset("SO","ALL");
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

    //std::vector<SharedMatrix> Cs;
    //Cs.push_back(Cocca_);
    //Cs.push_back(Coccb_);
    //Cs.push_back(Cvira_);
    //Cs.push_back(Cvirb_);
    //C_ = Matrix::horzcat(Cs);
}

void UStab::print_header()
{
    boost::shared_ptr<Wavefunction> wfn = reference_wavefunction_;
    outfile->Printf( "\n");
    outfile->Printf( "         ------------------------------------------------------------\n");
    outfile->Printf( "                              UHF Stability code                     \n");
    outfile->Printf( "                                Jérôme Gonthier                     \n");
    outfile->Printf( "               Strong inspiration from R. Parrish's CIS              \n");
    outfile->Printf( "         ------------------------------------------------------------\n\n");

    outfile->Printf( "  ==> Geometry <==\n\n");
    molecule_->print();
    outfile->Printf( "  Nuclear repulsion = %20.15f\n", molecule_->nuclear_repulsion_energy());
    outfile->Printf( "  Reference energy  = %20.15f\n\n", Eref_);

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
    boost::shared_ptr<USTABHamiltonian> H(new USTABHamiltonian(jk_, Cocca_,Cvira_,Coccb_,Cvirb_,eps_occa_,eps_vira_,
                            eps_occb_,eps_virb_));
    boost::shared_ptr<DLUSolver> solver;
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

    // Unpack
    const std::vector<boost::shared_ptr<Vector> > stabvecs = solver->eigenvectors();
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

    // Print wavefunctions and properties
/*    sort_states();
    print_wavefunctions();
    print_amplitudes();
    print_transitions();
    print_densities();
*/
    return 0.0;
}

void UStab::analyze()
{

    // We use the convergence criterion to eliminate zero eigenvalues
    // suffering from numerical noise.
    for (int i = 0; i < vals_.size(); ++i) {
        if ( vals_[i] > convergence_ ) {
            break;
        } else if ( abs(vals_[i]) > convergence_ ) {
            if ( vecs_[i].first->symmetry() == 0) {
                unstable = true;
                unstable_val = vals_[i];
                unstable_vec = vecs_[i];
                break;
            }
        }
    }

    if (unstable) {
        outfile->Printf("    Negative totally symmetric eigenvalue detected: %f \n", unstable_val);
        outfile->Printf("    Wavefunction unstable!\n");
    } else {
        outfile->Printf("    Wavefunction stable under totally symmetric rotations.\n");
        outfile->Printf("    Lowest eigenvalue: %f \n", vals_[0]);
    }
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
            for (int a = 0; a < nvira; ++a) {
                Ca_->rotate_columns(h, i, a, scale*unveca->get(i,a));
            }
        }
        int noccb = unvecb->rowdim(h);
        int nvirb = unvecb->coldim(h);
    // Rotate the beta orbitals
        for (int i = 0; i < noccb; ++i) {
            for (int a = 0; a < nvirb; ++a) {
                Cb_->rotate_columns(h, i, a, scale*unvecb->get(i,a));
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
            jk_ = JK::build_JK();
            unsigned long int effective_memory = (unsigned long int)(0.125 * options_.get_double("CPHF_MEM_SAFETY_FACTOR") * memory_);
            jk_->set_memory(effective_memory);
            jk_->initialize();
        }
    }

  }

    }} // End namespaces

