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


#include "psi4/libqt/qt.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/liboptions/liboptions_python.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/integral.h"
#include "psi4/libmints/vector.h"
#include "psi4/libmints/extern.h"
#include "psi4/psi4-dec.h"
#include "psi4/libfock/v.h"
#include "psi4/libfunctional/superfunctional.h"
#include "psi4/libdisp/dispersion.h"
#include "psi4/libscf_solver/hf.h"
#include "scf_grad.h"
#include "jk_grad.h"

#include <algorithm>

#include <sstream>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;
using namespace psi;


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
    if (ref_wfn->arrays().count("-D Gradient")){
        gradients_["-D"] = ref_wfn->get_array("-D Gradient");
    }

}
SCFGrad::~SCFGrad()
{
}
void SCFGrad::common_init()
{

    // if (!reference_wavefunction_) {
    //     throw PSIEXCEPTION("SCFGrad: Run SCF first");
    // }

    //if (options_.get_str("REFERENCE") == "ROHF" || options_.get_str("REFERENCE") == "CUHF")
    //    reference_wavefunction_->semicanonicalize();

    //copy(reference_wavefunction_);

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

    outfile->Printf( "  Nuclear repulsion = %20.15f\n", basisset_->molecule()->nuclear_repulsion_energy());
    outfile->Printf( "\n");

    outfile->Printf( "  ==> Basis Set <==\n\n");
    basisset_->print_by_level("outfile", print_);

    // => Registers <= //

    // std::map<std::string, SharedMatrix> gradients;

    std::vector<std::string> gradient_terms;
    gradient_terms.push_back("Nuclear");
    gradient_terms.push_back("Kinetic");
    gradient_terms.push_back("Potential");
    gradient_terms.push_back("Overlap");
    gradient_terms.push_back("Coulomb");
    if(options_.get_bool("PERTURB_H"))
        gradient_terms.push_back("Perturbation");
    gradient_terms.push_back("Exchange");
    gradient_terms.push_back("Exchange,LR");
    gradient_terms.push_back("XC");
    gradient_terms.push_back("-D");
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
    if (options_.get_str("REFERENCE") == "RKS") {
        std::vector<SharedMatrix>& C = potential_->C();
        C.clear();
        C.push_back(Ca_subset("SO", "OCC"));
    } else if (options_.get_str("REFERENCE") == "UKS") {
        std::vector<SharedMatrix>& C = potential_->C();
        C.clear();
        C.push_back(Ca_subset("SO", "OCC"));
        C.push_back(Cb_subset("SO", "OCC"));
    }

    // => Sizings <= //
    int natom = molecule_->natom();
    int nso = basisset_->nbf();
    int nalpha = Ca->colspi()[0];
    int nbeta = Cb->colspi()[0];

    // => Nuclear Gradient <= //
    gradients_["Nuclear"] = SharedMatrix(molecule_->nuclear_repulsion_energy_deriv1().clone());
    gradients_["Nuclear"]->set_name("Nuclear Gradient");

    // => Kinetic Gradient <= //
    timer_on("Grad: T");
    {
        double** Dp = Dt->pointer();

        gradients_["Kinetic"] = SharedMatrix(gradients_["Nuclear"]->clone());
        gradients_["Kinetic"]->set_name("Kinetic Gradient");
        gradients_["Kinetic"]->zero();
        double** Tp = gradients_["Kinetic"]->pointer();

        // Kinetic derivatives
        std::shared_ptr<OneBodyAOInt> Tint(integral_->ao_kinetic(1));
        const double* buffer = Tint->buffer();

        for (int P = 0; P < basisset_->nshell(); P++) {
            for (int Q = 0; Q <= P; Q++) {

                Tint->compute_shell_deriv1(P,Q);

                int nP = basisset_->shell(P).nfunction();
                int oP = basisset_->shell(P).function_index();
                int aP = basisset_->shell(P).ncenter();

                int nQ = basisset_->shell(Q).nfunction();
                int oQ = basisset_->shell(Q).function_index();
                int aQ = basisset_->shell(Q).ncenter();

                int offset = nP * nQ;
                const double* ref = buffer;
                double perm = (P == Q ? 1.0 : 2.0);

                // Px
                for (int p = 0; p < nP; p++) {
                    for (int q = 0; q < nQ; q++) {
                        Tp[aP][0] += perm * Dp[p + oP][q + oQ] * (*ref++);
                    }
                }

                // Py
                for (int p = 0; p < nP; p++) {
                    for (int q = 0; q < nQ; q++) {
                        Tp[aP][1] += perm * Dp[p + oP][q + oQ] * (*ref++);
                    }
                }

                // Pz
                for (int p = 0; p < nP; p++) {
                    for (int q = 0; q < nQ; q++) {
                        Tp[aP][2] += perm * Dp[p + oP][q + oQ] * (*ref++);
                    }
                }

                // Qx
                for (int p = 0; p < nP; p++) {
                    for (int q = 0; q < nQ; q++) {
                        Tp[aQ][0] += perm * Dp[p + oP][q + oQ] * (*ref++);
                    }
                }

                // Qy
                for (int p = 0; p < nP; p++) {
                    for (int q = 0; q < nQ; q++) {
                        Tp[aQ][1] += perm * Dp[p + oP][q + oQ] * (*ref++);
                    }
                }

                // Qz
                for (int p = 0; p < nP; p++) {
                    for (int q = 0; q < nQ; q++) {
                        Tp[aQ][2] += perm * Dp[p + oP][q + oQ] * (*ref++);
                    }
                }
            }
        }
    }
    timer_off("Grad: T");

    // => Potential Gradient <= //
    timer_on("Grad: V");
    {
        double** Dp = Dt->pointer();

        gradients_["Potential"] = SharedMatrix(gradients_["Nuclear"]->clone());
        gradients_["Potential"]->set_name("Potential Gradient");
        gradients_["Potential"]->zero();

            // Thread count
            int threads = 1;
            #ifdef _OPENMP
                threads = Process::environment.get_n_threads();
            #endif

            // Potential derivatives
            std::vector<std::shared_ptr<OneBodyAOInt> > Vint;
            std::vector<SharedMatrix> Vtemps;
            for (int t = 0; t < threads; t++) {
                Vint.push_back(std::shared_ptr<OneBodyAOInt>(integral_->ao_potential(1)));
                Vtemps.push_back(SharedMatrix(gradients_["Potential"]->clone()));
            }

            // Lower Triangle
            std::vector<std::pair<int,int> > PQ_pairs;
            for (int P = 0; P < basisset_->nshell(); P++) {
                for (int Q = 0; Q <= P; Q++) {
                    PQ_pairs.push_back(std::pair<int,int>(P,Q));
                }
            }

        #pragma omp parallel for schedule(dynamic) num_threads(threads)
        for (long int PQ = 0L; PQ < PQ_pairs.size(); PQ++) {

            int P = PQ_pairs[PQ].first;
            int Q = PQ_pairs[PQ].second;

            int thread = 0;
            #ifdef _OPENMP
                thread = omp_get_thread_num();
            #endif

            Vint[thread]->compute_shell_deriv1(P,Q);
            const double* buffer = Vint[thread]->buffer();

            int nP = basisset_->shell(P).nfunction();
            int oP = basisset_->shell(P).function_index();
            int aP = basisset_->shell(P).ncenter();

            int nQ = basisset_->shell(Q).nfunction();
            int oQ = basisset_->shell(Q).function_index();
            int aQ = basisset_->shell(Q).ncenter();

            double perm = (P == Q ? 1.0 : 2.0);

            double** Vp = Vtemps[thread]->pointer();

            for (int A = 0; A < natom; A++) {
                const double* ref0 = &buffer[3 * A * nP * nQ + 0 * nP * nQ];
                const double* ref1 = &buffer[3 * A * nP * nQ + 1 * nP * nQ];
                const double* ref2 = &buffer[3 * A * nP * nQ + 2 * nP * nQ];
                for (int p = 0; p < nP; p++) {
                    for (int q = 0; q < nQ; q++) {
                        double Vval = perm * Dp[p + oP][q + oQ];
                        Vp[A][0] += Vval * (*ref0++);
                        Vp[A][1] += Vval * (*ref1++);
                        Vp[A][2] += Vval * (*ref2++);
                    }
                }
            }
        }

        for (int t = 0; t < threads; t++) {
            gradients_["Potential"]->add(Vtemps[t]);
        }
    }
    timer_off("Grad: V");

    // If an external field exists, add it to the one-electron Hamiltonian
    py::object pyExtern = dynamic_cast<PythonDataType*>(options_["EXTERN"].get())->to_python();
    if (pyExtern) {
        std::shared_ptr<ExternalPotential> external = pyExtern.cast<std::shared_ptr<ExternalPotential>>();
        if (external) {
            gradient_terms.push_back("External Potential");
            timer_on("Grad: External");
            gradients_["External Potential"] = external->computePotentialGradients(basisset_, Dt);
            timer_off("Grad: External");
        }  // end external
    }

    // => Perturbation Gradient <= //
    if(options_.get_bool("PERTURB_H")) {
        timer_on("Grad: Perturbation");

        double xlambda = 0.0;
        double ylambda = 0.0;
        double zlambda = 0.0;

        std::string perturb_with = options_.get_str("PERTURB_WITH");
        if (perturb_with == "DIPOLE_X")
            xlambda = options_.get_double("PERTURB_MAGNITUDE");
        else if (perturb_with == "DIPOLE_Y")
            ylambda = options_.get_double("PERTURB_MAGNITUDE");
        else if (perturb_with == "DIPOLE_Z")
            zlambda = options_.get_double("PERTURB_MAGNITUDE");
        else if (perturb_with == "DIPOLE") {
            if(options_["PERTURB_DIPOLE"].size() !=3)
                throw PSIEXCEPTION("The PERTURB dipole should have exactly three floating point numbers.");
            xlambda = options_["PERTURB_DIPOLE"][0].to_double();
            ylambda = options_["PERTURB_DIPOLE"][1].to_double();
            zlambda = options_["PERTURB_DIPOLE"][2].to_double();
        } else {
            std::string msg("Gradients for a ");
            msg += perturb_with;
            msg += " perturbation are not available yet.\n";
            throw PSIEXCEPTION(msg);
        }

        double** Dp = Dt->pointer();

        gradients_["Perturbation"] = SharedMatrix(gradients_["Nuclear"]->clone());
        gradients_["Perturbation"]->set_name("Perturbation Gradient");
        gradients_["Perturbation"]->zero();
        double** Pp = gradients_["Perturbation"]->pointer();

        // Nuclear dipole perturbation derivatives
        for(int n = 0; n < natom; ++n){
            double charge = molecule_->Z(n);
            Pp[n][0] += xlambda*charge;
            Pp[n][1] += ylambda*charge;
            Pp[n][2] += zlambda*charge;
        }

        // Electronic dipole perturbation derivatives
        std::shared_ptr<OneBodyAOInt> Dint(integral_->ao_dipole(1));
        const double* buffer = Dint->buffer();

        for (int P = 0; P < basisset_->nshell(); P++) {
            for (int Q = 0; Q <= P; Q++) {

                Dint->compute_shell_deriv1(P,Q);

                int nP = basisset_->shell(P).nfunction();
                int oP = basisset_->shell(P).function_index();
                int aP = basisset_->shell(P).ncenter();

                int nQ = basisset_->shell(Q).nfunction();
                int oQ = basisset_->shell(Q).function_index();
                int aQ = basisset_->shell(Q).ncenter();

                const double* ref = buffer;
                double perm = (P == Q ? 1.0 : 2.0);
                double prefac;

                /*
                 * Mu X derivatives
                 */
                if (xlambda != 0.0) {
                    prefac = perm*xlambda;
                    // Px
                    for (int p = 0; p < nP; p++) {
                        for (int q = 0; q < nQ; q++) {
                            Pp[aP][0] += prefac * Dp[p + oP][q + oQ] * (*ref++);
                        }
                    }

                    // Py
                    for (int p = 0; p < nP; p++) {
                        for (int q = 0; q < nQ; q++) {
                            Pp[aP][1] += prefac * Dp[p + oP][q + oQ] * (*ref++);
                        }
                    }

                    // Pz
                    for (int p = 0; p < nP; p++) {
                        for (int q = 0; q < nQ; q++) {
                            Pp[aP][2] += prefac * Dp[p + oP][q + oQ] * (*ref++);
                        }
                    }

                    // Qx
                    for (int p = 0; p < nP; p++) {
                        for (int q = 0; q < nQ; q++) {
                            Pp[aQ][0] += prefac * Dp[p + oP][q + oQ] * (*ref++);
                        }
                    }

                    // Qy
                    for (int p = 0; p < nP; p++) {
                        for (int q = 0; q < nQ; q++) {
                            Pp[aQ][1] += prefac * Dp[p + oP][q + oQ] * (*ref++);
                        }
                    }

                    // Qz
                    for (int p = 0; p < nP; p++) {
                        for (int q = 0; q < nQ; q++) {
                            Pp[aQ][2] += prefac * Dp[p + oP][q + oQ] * (*ref++);
                        }
                    }
                } else {
                    // Xlambda is zero, so we just advance the pointer to the buffer
                    ref += 6*nP*nQ;
                }

                /*
                 * Mu Y derivatives
                 */
                if (ylambda != 0.0) {
                    prefac = perm*ylambda;
                    // Px
                    for (int p = 0; p < nP; p++) {
                        for (int q = 0; q < nQ; q++) {
                            Pp[aP][0] += prefac * Dp[p + oP][q + oQ] * (*ref++);
                        }
                    }

                    // Py
                    for (int p = 0; p < nP; p++) {
                        for (int q = 0; q < nQ; q++) {
                            Pp[aP][1] += prefac * Dp[p + oP][q + oQ] * (*ref++);
                        }
                    }

                    // Pz
                    for (int p = 0; p < nP; p++) {
                        for (int q = 0; q < nQ; q++) {
                            Pp[aP][2] += prefac * Dp[p + oP][q + oQ] * (*ref++);
                        }
                    }

                    // Qx
                    for (int p = 0; p < nP; p++) {
                        for (int q = 0; q < nQ; q++) {
                            Pp[aQ][0] += prefac * Dp[p + oP][q + oQ] * (*ref++);
                        }
                    }

                    // Qy
                    for (int p = 0; p < nP; p++) {
                        for (int q = 0; q < nQ; q++) {
                            Pp[aQ][1] += prefac * Dp[p + oP][q + oQ] * (*ref++);
                        }
                    }

                    // Qz
                    for (int p = 0; p < nP; p++) {
                        for (int q = 0; q < nQ; q++) {
                            Pp[aQ][2] += prefac * Dp[p + oP][q + oQ] * (*ref++);
                        }
                    }
                } else {
                    // Ylambda is zero, so we just advance the pointer to the buffer
                    ref += 6*nP*nQ;
                }

                /*
                 * Mu Z derivatives
                 */
                if (zlambda != 0.0) {
                    prefac = perm*zlambda;
                    // Px
                    for (int p = 0; p < nP; p++) {
                        for (int q = 0; q < nQ; q++) {
                            Pp[aP][0] += prefac * Dp[p + oP][q + oQ] * (*ref++);
                        }
                    }

                    // Py
                    for (int p = 0; p < nP; p++) {
                        for (int q = 0; q < nQ; q++) {
                            Pp[aP][1] += prefac * Dp[p + oP][q + oQ] * (*ref++);
                        }
                    }

                    // Pz
                    for (int p = 0; p < nP; p++) {
                        for (int q = 0; q < nQ; q++) {
                            Pp[aP][2] += prefac * Dp[p + oP][q + oQ] * (*ref++);
                        }
                    }

                    // Qx
                    for (int p = 0; p < nP; p++) {
                        for (int q = 0; q < nQ; q++) {
                            Pp[aQ][0] += prefac * Dp[p + oP][q + oQ] * (*ref++);
                        }
                    }

                    // Qy
                    for (int p = 0; p < nP; p++) {
                        for (int q = 0; q < nQ; q++) {
                            Pp[aQ][1] += prefac * Dp[p + oP][q + oQ] * (*ref++);
                        }
                    }

                    // Qz
                    for (int p = 0; p < nP; p++) {
                        for (int q = 0; q < nQ; q++) {
                            Pp[aQ][2] += prefac * Dp[p + oP][q + oQ] * (*ref++);
                        }
                    }
                }

            }
        }
        timer_off("Grad: Perturbation");
    }

    // => Overlap Gradient <= //
    timer_on("Grad: S");
    {
        // Energy weighted density matrix
        SharedMatrix W(Da->clone());
        W->set_name("W");

        double** Wp = W->pointer();
        double** Cap = Ca->pointer();
        double** Cbp = Cb->pointer();
        double* eps_ap = eps_a->pointer();
        double* eps_bp = eps_b->pointer();

        double* temp = new double[nso * (ULI) nalpha];

        ::memset((void*) temp, '\0', sizeof(double) * nso * nalpha);
        for (int i = 0; i < nalpha; i++) {
            C_DAXPY(nso,eps_ap[i], &Cap[0][i], nalpha, &temp[i], nalpha);
        }

        C_DGEMM('N','T',nso,nso,nalpha,1.0,Cap[0],nalpha,temp,nalpha,0.0,Wp[0],nso);

        ::memset((void*) temp, '\0', sizeof(double) * nso * nbeta);
        for (int i = 0; i < nbeta; i++) {
            C_DAXPY(nso,eps_bp[i], &Cbp[0][i], nbeta, &temp[i], nbeta);
        }

        C_DGEMM('N','T',nso,nso,nbeta,1.0,Cbp[0],nbeta,temp,nbeta,1.0,Wp[0],nso);

        delete[] temp;

        gradients_["Overlap"] = SharedMatrix(gradients_["Nuclear"]->clone());
        gradients_["Overlap"]->set_name("Overlap Gradient");
        gradients_["Overlap"]->zero();
        double** Sp = gradients_["Overlap"]->pointer();

        // Overlap derivatives
        std::shared_ptr<OneBodyAOInt> Sint(integral_->ao_overlap(1));
        const double* buffer = Sint->buffer();

        for (int P = 0; P < basisset_->nshell(); P++) {
            for (int Q = 0; Q <= P; Q++) {

                Sint->compute_shell_deriv1(P,Q);

                int nP = basisset_->shell(P).nfunction();
                int oP = basisset_->shell(P).function_index();
                int aP = basisset_->shell(P).ncenter();

                int nQ = basisset_->shell(Q).nfunction();
                int oQ = basisset_->shell(Q).function_index();
                int aQ = basisset_->shell(Q).ncenter();

                int offset = nP * nQ;
                const double* ref = buffer;
                double perm = (P == Q ? 1.0 : 2.0);

                // Px
                for (int p = 0; p < nP; p++) {
                    for (int q = 0; q < nQ; q++) {
                        Sp[aP][0] -= perm * Wp[p + oP][q + oQ] * (*ref++);
                    }
                }

                // Py
                for (int p = 0; p < nP; p++) {
                    for (int q = 0; q < nQ; q++) {
                        Sp[aP][1] -= perm * Wp[p + oP][q + oQ] * (*ref++);
                    }
                }

                // Pz
                for (int p = 0; p < nP; p++) {
                    for (int q = 0; q < nQ; q++) {
                        Sp[aP][2] -= perm * Wp[p + oP][q + oQ] * (*ref++);
                    }
                }

                // Qx
                for (int p = 0; p < nP; p++) {
                    for (int q = 0; q < nQ; q++) {
                        Sp[aQ][0] -= perm * Wp[p + oP][q + oQ] * (*ref++);
                    }
                }

                // Qy
                for (int p = 0; p < nP; p++) {
                    for (int q = 0; q < nQ; q++) {
                        Sp[aQ][1] -= perm * Wp[p + oP][q + oQ] * (*ref++);
                    }
                }

                // Qz
                for (int p = 0; p < nP; p++) {
                    for (int q = 0; q < nQ; q++) {
                        Sp[aQ][2] -= perm * Wp[p + oP][q + oQ] * (*ref++);
                    }
                }
            }
        }
    }
    timer_off("Grad: S");

    // => Two-Electron Gradient <= //
    timer_on("Grad: JK");

    std::shared_ptr<JKGrad> jk = JKGrad::build_JKGrad(1, basisset_, basissets_["DF_BASIS_SCF"]);
    jk->set_memory((ULI) (options_.get_double("SCF_MEM_SAFETY_FACTOR") * memory_ / 8L));

    jk->set_Ca(Ca);
    jk->set_Cb(Cb);
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
        gradients_["Exchange"]->scale(-(functional_->x_alpha()));
    }
    if (functional_->is_x_lrc()) {
        gradients_["Exchange,LR"] = jk_gradients["Exchange,LR"];
        gradients_["Exchange,LR"]->scale(-(1.0 - functional_->x_alpha()));
    }
    timer_off("Grad: JK");

    // => XC Gradient <= //
    timer_on("Grad: XC");
    if (potential_) {
        potential_->print_header();
        gradients_["XC"] = potential_->compute_gradient();
    }
    timer_off("Grad: XC");

    // => -D Gradient <= //
    // py::object
    // for (auto item : )
    // gradients_["-D"] = py_grad.cast<SharedMatrix>();
    // if (functional_ && functional_->dispersion()) {
    //     gradients_["-D"] = functional_->dispersion()->compute_gradient(basisset_->molecule());
    // }

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

    outfile->Printf( "  Nuclear repulsion = %20.15f\n", basisset_->molecule()->nuclear_repulsion_energy());
    outfile->Printf( "\n");

    outfile->Printf( "  ==> Basis Set <==\n\n");
    basisset_->print_by_level("outfile", print_);

    // => Registers <= //

    std::map<std::string, SharedMatrix> hessians;

    std::vector<std::string> hessian_terms;
    hessian_terms.push_back("Nuclear");
    hessian_terms.push_back("Kinetic");
    hessian_terms.push_back("Potential");
    hessian_terms.push_back("Overlap");
    hessian_terms.push_back("Coulomb");
    hessian_terms.push_back("Exchange");
    hessian_terms.push_back("Exchange,LR");
    hessian_terms.push_back("XC");
    hessian_terms.push_back("-D");
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

    if (options_.get_str("REFERENCE") == "RKS") {
        throw PSIEXCEPTION("Missing XC derivates for Hessians");
        //potential = VBase::build_V(basisset_, options_, "RV");
        //potential->initialize();
        //std::vector<SharedMatrix>& C = potential->C();
        //C.push_back(Ca_subset("SO", "OCC"));
        //functional = potential->functional();
    } else if (options_.get_str("REFERENCE") == "UKS") {
        throw PSIEXCEPTION("Missing XC derivates for Hessians");
        //potential = VBase::build_V(basisset_, options_, "UV");
        //potential->initialize();
        //std::vector<SharedMatrix>& C = potential->C();
        //C.push_back(Ca_subset("SO", "OCC"));
        //C.push_back(Cb_subset("SO", "OCC"));
        //functional = potential->functional();
    }

    // => Sizings <= //
    int natom = molecule_->natom();
    int nso = basisset_->nbf();
    int nalpha = Ca->colspi()[0];
    int nbeta = Cb->colspi()[0];

    // => Nuclear Hessian <= //
    hessians["Nuclear"] = SharedMatrix(molecule_->nuclear_repulsion_energy_deriv2().clone());
    hessians["Nuclear"]->set_name("Nuclear Hessian");

    SharedMatrix Zxyz(new Matrix("Zxyz", 1, 4));

    // => Potential Hessian <= //
    timer_on("Hess: V");
    {
        double** Dp = Dt->pointer();

        hessians["Potential"] = SharedMatrix(hessians["Nuclear"]->clone());
        hessians["Potential"]->set_name("Potential Hessian");
        hessians["Potential"]->zero();
        double** Vp = hessians["Potential"]->pointer();

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
        int dim = hessians["Potential"]->rowdim();
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

        hessians["Kinetic"] = SharedMatrix(hessians["Nuclear"]->clone());
        hessians["Kinetic"]->set_name("Kinetic Hessian");
        hessians["Kinetic"]->zero();
        double** Tp = hessians["Kinetic"]->pointer();

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
        int dim = hessians["Kinetic"]->rowdim();
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

        double* temp = new double[nso * (ULI) nalpha];

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

        hessians["Overlap"] = SharedMatrix(hessians["Nuclear"]->clone());
        hessians["Overlap"]->set_name("Overlap Hessian");
        hessians["Overlap"]->zero();
        double** Sp = hessians["Overlap"]->pointer();

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
        int dim = hessians["Overlap"]->rowdim();
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
    jk->set_memory((ULI) (options_.get_double("SCF_MEM_SAFETY_FACTOR") * memory_ / 8L));

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
        hessians["Coulomb"] = jk_hessians["Coulomb"];
        if (functional->is_x_hybrid()) {
            hessians["Exchange"] = jk_hessians["Exchange"];
            hessians["Exchange"]->scale(-(functional->x_alpha()));
        }
        if (functional->is_x_lrc()) {
            hessians["Exchange,LR"] = jk_hessians["Exchange,LR"];
            hessians["Exchange,LR"]->scale(-(1.0 - functional->x_alpha()));
        }
    } else {
        hessians["Coulomb"] = jk_hessians["Coulomb"];
        hessians["Exchange"] = jk_hessians["Exchange"];
        hessians["Exchange"]->scale(-1.0);
    }
    timer_off("Hess: JK");

    // => XC Hessian <= //
    timer_on("Hess: XC");
    if (functional) {
        potential->print_header();
        throw PSIEXCEPTION("KS Hessians not implemented");
        //hessians["XC"] = potential->compute_hessian();
    }
    timer_off("Hess: XC");

    // => Response Terms (Brace Yourself) <= //
    if (options_.get_str("REFERENCE") == "RHF") {
        hessians["Response"] = rhf_hessian_response();
    } else {
        throw PSIEXCEPTION("SCFHessian: Response not implemented for this reference");
    }

    // => Total Hessian <= //
    SharedMatrix total = SharedMatrix(hessians["Nuclear"]->clone());
    total->zero();

    for (int i = 0; i < hessian_terms.size(); i++) {
        if (hessians.count(hessian_terms[i])) {
            total->add(hessians[hessian_terms[i]]);
        }
    }

    hessians["Total"] = total;
    hessians["Total"]->set_name("Total Hessian");

    // => Final Printing <= //
    if (print_ > 1) {
        for (int i = 0; i < hessian_terms.size(); i++) {
            if (hessians.count(hessian_terms[i])) {
                hessians[hessian_terms[i]]->print();
            }
        }
    } else {
        hessians["Total"]->print();
    }

    return hessians["Total"];
}

}} // Namespaces
