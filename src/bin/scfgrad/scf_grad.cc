/*
 *@BEGIN LICENSE
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

#include <libmints/mints.h>
#include <libqt/qt.h>
#include <libpsio/psio.hpp>
#include <psi4-dec.h>
#include <libfock/v.h>
#include <libfunctional/superfunctional.h>
#include <libdisp/dispersion.h>
#include "scf_grad.h"
#include "jk_grad.h"

#include <algorithm>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>

#include <sstream>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;
using namespace psi;
using namespace boost;

namespace psi {
namespace scfgrad {

SCFGrad::SCFGrad() :
    Wavefunction(Process::environment.options,_default_psio_lib_) 
{
    common_init();
}
SCFGrad::~SCFGrad()
{
}
void SCFGrad::common_init()
{
    reference_wavefunction_ = Process::environment.wavefunction();
    
    if (!reference_wavefunction_) {
        throw PSIEXCEPTION("SCFGrad: Run SCF first");
    }

    if (options_.get_str("REFERENCE") == "ROHF" || options_.get_str("REFERENCE") == "CUHF")
        reference_wavefunction_->semicanonicalize();

    copy(reference_wavefunction_);
    
    print_ = options_.get_int("PRINT");
    debug_ = options_.get_int("DEBUG");
}
SharedMatrix SCFGrad::compute_gradient()
{
    // => Echo <= //

    fprintf(outfile, "\n");
    fprintf(outfile, "         ------------------------------------------------------------\n");
    fprintf(outfile, "                                   SCF GRAD                          \n"); 
    fprintf(outfile, "                          Rob Parrish, Justin Turney,                \n"); 
    fprintf(outfile, "                       Andy Simmonett, and Alex Sokolov              \n"); 
    fprintf(outfile, "         ------------------------------------------------------------\n\n");

    fprintf(outfile, "  ==> Geometry <==\n\n");
    molecule_->print();

    fprintf(outfile, "  Nuclear repulsion = %20.15f\n", basisset_->molecule()->nuclear_repulsion_energy());
    fprintf(outfile, "\n");

    fprintf(outfile, "  ==> Basis Set <==\n\n");
    basisset_->print_by_level(outfile, print_);

    // => Registers <= //

    std::map<std::string, SharedMatrix> gradients;

    std::vector<std::string> gradient_terms;
    gradient_terms.push_back("Nuclear");
    gradient_terms.push_back("Kinetic");
    gradient_terms.push_back("Potential");
    gradient_terms.push_back("Overlap");
    gradient_terms.push_back("Coulomb");
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
    boost::shared_ptr<SuperFunctional> functional;
    boost::shared_ptr<VBase> potential;

    if (options_.get_str("REFERENCE") == "RKS") {
        potential = VBase::build_V(options_, "RV"); 
        potential->initialize();
        std::vector<SharedMatrix>& C = potential->C();
        C.push_back(Ca_subset("SO", "OCC"));
        functional = potential->functional();
    } else if (options_.get_str("REFERENCE") == "UKS") { 
        potential = VBase::build_V(options_, "UV"); 
        potential->initialize();
        std::vector<SharedMatrix>& C = potential->C();
        C.push_back(Ca_subset("SO", "OCC"));
        C.push_back(Cb_subset("SO", "OCC"));
        functional = potential->functional();
    }

    // => Sizings <= //
    int natom = molecule_->natom();
    int nso = basisset_->nbf();
    int nalpha = Ca->colspi()[0];
    int nbeta = Cb->colspi()[0];
     
    // => Nuclear Gradient <= //
    gradients["Nuclear"] = SharedMatrix(molecule_->nuclear_repulsion_energy_deriv1().clone());
    gradients["Nuclear"]->set_name("Nuclear Gradient");

    // => Kinetic Gradient <= //
    timer_on("Grad: T");
    {
        double** Dp = Dt->pointer();

        gradients["Kinetic"] = SharedMatrix(gradients["Nuclear"]->clone());
        gradients["Kinetic"]->set_name("Kinetic Gradient");
        gradients["Kinetic"]->zero();
        double** Tp = gradients["Kinetic"]->pointer();

        // Kinetic derivatives
        boost::shared_ptr<OneBodyAOInt> Tint(integral_->ao_kinetic(1));
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

        gradients["Potential"] = SharedMatrix(gradients["Nuclear"]->clone());
        gradients["Potential"]->set_name("Potential Gradient");
        gradients["Potential"]->zero();

        // Thread count
        int threads = 1;
        #ifdef _OPENMP
            threads = omp_get_max_threads();
        #endif

        // Potential derivatives
        std::vector<boost::shared_ptr<OneBodyAOInt> > Vint;
        std::vector<SharedMatrix> Vtemps;
        for (int t = 0; t < threads; t++) { 
            Vint.push_back(boost::shared_ptr<OneBodyAOInt>(integral_->ao_potential(1)));
            Vtemps.push_back(SharedMatrix(gradients["Potential"]->clone()));
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
            gradients["Potential"]->add(Vtemps[t]);
        }
    }
    timer_off("Grad: V");

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

        gradients["Overlap"] = SharedMatrix(gradients["Nuclear"]->clone());
        gradients["Overlap"]->set_name("Overlap Gradient");
        gradients["Overlap"]->zero();
        double** Sp = gradients["Overlap"]->pointer();

        // Overlap derivatives
        boost::shared_ptr<OneBodyAOInt> Sint(integral_->ao_overlap(1));
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

    boost::shared_ptr<JKGrad> jk = JKGrad::build_JKGrad(1);
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
    jk->compute_gradient();    

    std::map<std::string, SharedMatrix>& jk_gradients = jk->gradients();
    if (functional) {
        gradients["Coulomb"] = jk_gradients["Coulomb"];
        if (functional->is_x_hybrid()) {
            gradients["Exchange"] = jk_gradients["Exchange"];
            gradients["Exchange"]->scale(-(functional->x_alpha()));
        } 
        if (functional->is_x_lrc()) {
            gradients["Exchange,LR"] = jk_gradients["Exchange,LR"];
            gradients["Exchange,LR"]->scale(-(1.0 - functional->x_alpha()));
        }
    } else {
        gradients["Coulomb"] = jk_gradients["Coulomb"];
        gradients["Exchange"] = jk_gradients["Exchange"];
        gradients["Exchange"]->scale(-1.0);
    }
    timer_off("Grad: JK");

    // => XC Gradient <= //
    timer_on("Grad: XC");
    if (functional) {
        potential->print_header();
        gradients["XC"] = potential->compute_gradient();
    }
    timer_off("Grad: XC");

    // => -D Gradient <= //
    if (functional && functional->dispersion()) {
        gradients["-D"] = functional->dispersion()->compute_gradient(basisset_->molecule());
    }

    // => Total Gradient <= //
    SharedMatrix total = SharedMatrix(gradients["Nuclear"]->clone());
    total->zero();

    for (int i = 0; i < gradient_terms.size(); i++) {
        if (gradients.count(gradient_terms[i])) {
            total->add(gradients[gradient_terms[i]]); 
        }
    }

    gradients["Total"] = total; 
    gradients["Total"]->set_name("Total Gradient");

    // => Final Printing <= //
    if (print_ > 1) {
        for (int i = 0; i < gradient_terms.size(); i++) {
            if (gradients.count(gradient_terms[i])) {
                gradients[gradient_terms[i]]->print_atom_vector(); 
            }
        }
    } else {
        gradients["Total"]->print_atom_vector();
    }

    return gradients["Total"]; 
}
SharedMatrix SCFGrad::compute_hessian()
{
    // => Echo <= //

    fprintf(outfile, "\n");
    fprintf(outfile, "         ------------------------------------------------------------\n");
    fprintf(outfile, "                                   SCF HESS                          \n"); 
    fprintf(outfile, "                          Rob Parrish, Justin Turney,                \n"); 
    fprintf(outfile, "                       Andy Simmonett, and Alex Sokolov              \n"); 
    fprintf(outfile, "         ------------------------------------------------------------\n\n");

    fprintf(outfile, "  ==> Geometry <==\n\n");
    molecule_->print();

    fprintf(outfile, "  Nuclear repulsion = %20.15f\n", basisset_->molecule()->nuclear_repulsion_energy());
    fprintf(outfile, "\n");

    fprintf(outfile, "  ==> Basis Set <==\n\n");
    basisset_->print_by_level(outfile, print_);

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
    boost::shared_ptr<SuperFunctional> functional;
    boost::shared_ptr<VBase> potential;

    if (options_.get_str("REFERENCE") == "RKS") {
        potential = VBase::build_V(options_, "RV"); 
        potential->initialize();
        std::vector<SharedMatrix>& C = potential->C();
        C.push_back(Ca_subset("SO", "OCC"));
        functional = potential->functional();
    } else if (options_.get_str("REFERENCE") == "UKS") { 
        potential = VBase::build_V(options_, "UV"); 
        potential->initialize();
        std::vector<SharedMatrix>& C = potential->C();
        C.push_back(Ca_subset("SO", "OCC"));
        C.push_back(Cb_subset("SO", "OCC"));
        functional = potential->functional();
    }

    // => Sizings <= //
    int natom = molecule_->natom();
    int nso = basisset_->nbf();
    int nalpha = Ca->colspi()[0];
    int nbeta = Cb->colspi()[0];
     
    // => Nuclear Hessian <= //
    hessians["Nuclear"] = SharedMatrix(molecule_->nuclear_repulsion_energy_deriv2().clone());
    hessians["Nuclear"]->set_name("Nuclear Hessian");

    // => Kinetic Hessian <= //
    timer_on("Hess: T");
    {
        double** Dp = Dt->pointer();

        hessians["Kinetic"] = SharedMatrix(hessians["Nuclear"]->clone());
        hessians["Kinetic"]->set_name("Kinetic Hessian");
        hessians["Kinetic"]->zero();
        double** Tp = hessians["Kinetic"]->pointer();

        // Kinetic derivatives
        boost::shared_ptr<OneBodyAOInt> Tint(integral_->ao_kinetic(2));
        const double* buffer = Tint->buffer();   

        for (int P = 0; P < basisset_->nshell(); P++) {
            for (int Q = 0; Q <= P; Q++) {

                Tint->compute_shell_deriv2(P,Q);
                                
                int nP = basisset_->shell(P).nfunction();
                int oP = basisset_->shell(P).function_index();
                int aP = basisset_->shell(P).ncenter();
 
                int nQ = basisset_->shell(Q).nfunction();
                int oQ = basisset_->shell(Q).function_index();
                int aQ = basisset_->shell(Q).ncenter();

                int offset = nP * nQ;
                const double* ref = buffer;
                double perm = (P == Q ? 1.0 : 2.0);

                int Px = 3 * aP + 0;
                int Py = 3 * aP + 1;
                int Pz = 3 * aP + 2;

                int Qx = 3 * aQ + 0;
                int Qy = 3 * aQ + 1;
                int Qz = 3 * aQ + 2;
               
                // Px Qx 
                for (int p = 0; p < nP; p++) {
                    for (int q = 0; q < nQ; q++) {
                        Tp[Px][Qx] += perm * Dp[p + oP][q + oQ] * (*ref);
                        if (aP != aQ)
                            Tp[Qx][Px] += perm * Dp[p + oP][q + oQ] * (*ref);
                        ref++;
                    }
                }
               
                // Px Qy 
                for (int p = 0; p < nP; p++) {
                    for (int q = 0; q < nQ; q++) {
                        Tp[Px][Qy] += perm * Dp[p + oP][q + oQ] * (*ref);
                        Tp[Qy][Px] += perm * Dp[p + oP][q + oQ] * (*ref);
                        ref++;
                    }
                }
               
                // Px Pz 
                for (int p = 0; p < nP; p++) {
                    for (int q = 0; q < nQ; q++) {
                        Tp[Px][Qz] += perm * Dp[p + oP][q + oQ] * (*ref);
                        Tp[Qz][Px] += perm * Dp[p + oP][q + oQ] * (*ref);
                        ref++;
                    }
                }
               
                // Py Qy 
                for (int p = 0; p < nP; p++) {
                    for (int q = 0; q < nQ; q++) {
                        Tp[Py][Qy] += perm * Dp[p + oP][q + oQ] * (*ref);
                        if (aP != aQ)
                            Tp[Qy][Py] += perm * Dp[p + oP][q + oQ] * (*ref);
                        ref++;
                    }
                }
               
                // Py Qz 
                for (int p = 0; p < nP; p++) {
                    for (int q = 0; q < nQ; q++) {
                        Tp[Py][Qz] += perm * Dp[p + oP][q + oQ] * (*ref);
                        Tp[Qz][Py] += perm * Dp[p + oP][q + oQ] * (*ref);
                        ref++;
                    }
                }
               
                // Pz Qz 
                for (int p = 0; p < nP; p++) {
                    for (int q = 0; q < nQ; q++) {
                        Tp[Pz][Qz] += perm * Dp[p + oP][q + oQ] * (*ref);
                        if (aP != aQ)
                            Tp[Qz][Pz] += perm * Dp[p + oP][q + oQ] * (*ref);
                        ref++;
                    }
                }
            }
        } 
    }
    timer_off("Hess: T");

    // => Potential Hessian <= //
    // TODO
    timer_on("Hess: V");
    {
        double** Dp = Dt->pointer();

        hessians["Potential"] = SharedMatrix(hessians["Nuclear"]->clone());
        hessians["Potential"]->set_name("Potential Hessian");
        hessians["Potential"]->zero();

        // Thread count
        int threads = 1;
        #ifdef _OPENMP
            threads = omp_get_max_threads();
        #endif

        // Potential derivatives
        std::vector<boost::shared_ptr<OneBodyAOInt> > Vint;
        std::vector<SharedMatrix> Vtemps;
        for (int t = 0; t < threads; t++) { 
            Vint.push_back(boost::shared_ptr<OneBodyAOInt>(integral_->ao_potential(2)));
            Vtemps.push_back(SharedMatrix(hessians["Potential"]->clone()));
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

            Vint[thread]->compute_shell_deriv2(P,Q);
            const double* buffer = Vint[thread]->buffer();
                            
            int nP = basisset_->shell(P).nfunction();
            int oP = basisset_->shell(P).function_index();
            int aP = basisset_->shell(P).ncenter();
 
            int nQ = basisset_->shell(Q).nfunction();
            int oQ = basisset_->shell(Q).function_index();
            int aQ = basisset_->shell(Q).ncenter();

            double perm = (P == Q ? 1.0 : 2.0);
                
            double** Vp = Vtemps[thread]->pointer();

            for (int alpha = 0; alpha < 3 * natom; alpha++) {
                for (int beta = 0; beta < 3 * natom; beta++) {
                    for (int p = 0; p < nP; p++) {
                        for (int q = 0; q < nQ; q++) {
                            double Vval = perm * Dp[p + oP][q + oQ];
                            Vp[alpha][beta] += Vval * (*buffer++);
                        }
                    }
                }
            }
        } 
    
        for (int t = 0; t < threads; t++) { 
            hessians["Potential"]->add(Vtemps[t]);
        }
    }
    timer_off("Hess: V");

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

        C_DGEMM('N','T',nso,nso,nalpha,1.0,Cap[0],nalpha,temp,nalpha,0.0,Wp[0],nso);

        ::memset((void*) temp, '\0', sizeof(double) * nso * nbeta);
        for (int i = 0; i < nbeta; i++) {
            C_DAXPY(nso,eps_bp[i], &Cbp[0][i], nbeta, &temp[i], nbeta);
        }

        C_DGEMM('N','T',nso,nso,nbeta,1.0,Cbp[0],nbeta,temp,nbeta,1.0,Wp[0],nso);
        
        delete[] temp;

        hessians["Overlap"] = SharedMatrix(hessians["Nuclear"]->clone());
        hessians["Overlap"]->set_name("Overlap Hessian");
        hessians["Overlap"]->zero();
        double** Sp = hessians["Overlap"]->pointer();

        // Overlap derivatives
        boost::shared_ptr<OneBodyAOInt> Sint(integral_->ao_overlap(2));
        const double* buffer = Sint->buffer();   

        for (int P = 0; P < basisset_->nshell(); P++) {
            for (int Q = 0; Q <= P; Q++) {

                Sint->compute_shell_deriv2(P,Q);
                                
                int nP = basisset_->shell(P).nfunction();
                int oP = basisset_->shell(P).function_index();
                int aP = basisset_->shell(P).ncenter();
 
                int nQ = basisset_->shell(Q).nfunction();
                int oQ = basisset_->shell(Q).function_index();
                int aQ = basisset_->shell(Q).ncenter();

                int offset = nP * nQ;
                const double* ref = buffer;
                double perm = (P == Q ? 1.0 : 2.0);
               
                int Px = 3 * aP + 0;
                int Py = 3 * aP + 1;
                int Pz = 3 * aP + 2;

                int Qx = 3 * aQ + 0;
                int Qy = 3 * aQ + 1;
                int Qz = 3 * aQ + 2;
               
                // Px Qx 
                for (int p = 0; p < nP; p++) {
                    for (int q = 0; q < nQ; q++) {
                        Sp[Px][Qx] += perm * Wp[p + oP][q + oQ] * (*ref);
                        if (aP != aQ)
                            Sp[Qx][Px] += perm * Wp[p + oP][q + oQ] * (*ref);
                        ref++;
                    }
                }
               
                // Px Qy 
                for (int p = 0; p < nP; p++) {
                    for (int q = 0; q < nQ; q++) {
                        Sp[Px][Qy] += perm * Wp[p + oP][q + oQ] * (*ref);
                        Sp[Qy][Px] += perm * Wp[p + oP][q + oQ] * (*ref);
                        ref++;
                    }
                }
               
                // Px Pz 
                for (int p = 0; p < nP; p++) {
                    for (int q = 0; q < nQ; q++) {
                        Sp[Px][Qz] += perm * Wp[p + oP][q + oQ] * (*ref);
                        Sp[Qz][Px] += perm * Wp[p + oP][q + oQ] * (*ref);
                        ref++;
                    }
                }
               
                // Py Qy 
                for (int p = 0; p < nP; p++) {
                    for (int q = 0; q < nQ; q++) {
                        Sp[Py][Qy] += perm * Wp[p + oP][q + oQ] * (*ref);
                        if (aP != aQ)
                            Sp[Qy][Py] += perm * Wp[p + oP][q + oQ] * (*ref);
                        ref++;
                    }
                }
               
                // Py Qz 
                for (int p = 0; p < nP; p++) {
                    for (int q = 0; q < nQ; q++) {
                        Sp[Py][Qz] += perm * Wp[p + oP][q + oQ] * (*ref);
                        Sp[Qz][Py] += perm * Wp[p + oP][q + oQ] * (*ref);
                        ref++;
                    }
                }
               
                // Pz Qz 
                for (int p = 0; p < nP; p++) {
                    for (int q = 0; q < nQ; q++) {
                        Sp[Pz][Qz] += perm * Wp[p + oP][q + oQ] * (*ref);
                        if (aP != aQ)
                            Sp[Qz][Pz] += perm * Wp[p + oP][q + oQ] * (*ref);
                        ref++;
                    }
                }
            }
        } 
    }
    timer_off("Hess: S");

    // => Two-Electron Hessian <= //
    /**

    timer_on("Hess: JK");

    boost::shared_ptr<JKGrad> jk = JKGrad::build_JKGrad(2);
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

    std::map<std::string, SharedMatrix>& jk_hessians = jk->gradients();
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
    
    **/

    // => XC Hessian <= //
    timer_on("Hess: XC");
    if (functional) {
        potential->print_header();
        throw PSIEXCEPTION("KS Hessians not implemented");
        //hessians["XC"] = potential->compute_hessian();
    }
    timer_off("Hess: XC");

    // => -D Hessian <= //
    if (functional && functional->dispersion()) {
        throw PSIEXCEPTION("-D Hessians not implemented");
        //hessians["-D"] = functional->dispersion()->compute_hessian(basisset_->molecule());
    }

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
