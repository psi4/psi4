#include "3index.h"

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <vector>
#include <utility>

#include <psifiles.h>
#include <libpsio/psio.h>
#include <libqt/qt.h>
#include <libciomr/libciomr.h>
#include <libmints/mints.h>
#include <libmints/psimath.h>

//MKL Header
#ifdef HAVE_MKL
#include <mkl.h>
#endif

//OpenMP Header
//_OPENMP is defined by the compiler if it exists
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;
using namespace psi;

namespace psi { 

FittingMetric::FittingMetric()
{
    Options& opt = Process::environment.options;
    shared_ptr<Molecule> molecule = Process::environment.molecule();

    if (molecule.get() == 0) {
        fprintf(outfile, "  Active molecule not set!");
        throw PSIEXCEPTION("Active molecule not set!");
    }

    // Make sure molecule is valid.
    molecule->update_geometry();

    // Grab the auxiliary bases (use the SCF tags for now)
    shared_ptr<BasisSetParser> parser(new Gaussian94BasisSetParser(opt.get_str("BASIS_PATH")));
    aux_ = BasisSet::construct(parser, molecule, opt.get_str("RI_BASIS_SCF"));
    if (opt.get_str("POISSON_BASIS_SCF") != "") {
        pois_ =  BasisSet::construct(parser, molecule, opt.get_str("POISSON_BASIS_SCF"));
        is_poisson_ = true;
    } else {
        is_poisson_ = false;
    }
    is_inverted_ = false;
}
FittingMetric::FittingMetric(shared_ptr<BasisSet> aux) : 
    aux_(aux), is_poisson_(false), is_inverted_(false)
{
}
FittingMetric::FittingMetric(shared_ptr<BasisSet> aux, shared_ptr<BasisSet> pois) : 
    aux_(aux), pois_(pois), is_poisson_(true), is_inverted_(false)
{
}

FittingMetric::~FittingMetric()
{
}

void FittingMetric::form_fitting_metric()
{
    is_inverted_ = false;
    algorithm_ = "NONE";

    // Sizing/symmetry indexing
    shared_ptr<IntegralFactory> auxfact(new IntegralFactory(aux_, aux_, aux_, aux_));   
    shared_ptr<PetiteList> auxpet(new PetiteList(aux_, auxfact));
    shared_ptr<IntegralFactory> poisfact;
    shared_ptr<PetiteList> poispet;
    if (is_poisson_) {
        poisfact = shared_ptr<IntegralFactory>(new IntegralFactory(pois_, pois_, pois_, pois_));
        poispet = shared_ptr<PetiteList>(new PetiteList(pois_, poisfact)); 
    }

    int naux = 0;
    int ngaussian = 0;
    Dimension nauxpi(auxpet->nirrep(), "Fitting Metric Dimensions");
    for (int h = 0; h < auxpet->nirrep(); h++) {
        naux += auxpet->SO_basisdim()[h];
        ngaussian += auxpet->SO_basisdim()[h];
        nauxpi[h] = auxpet->SO_basisdim()[h];
        if (is_poisson_) {
            naux += poispet->SO_basisdim()[h];
            nauxpi[h] = poispet->SO_basisdim()[h];
        } 
    }   
    Dimension ngauspi = auxpet->SO_basisdim();

    // Build the full DF/Poisson matrix in the AO basis first
    shared_ptr<Matrix> AOmetric(new Matrix("AO Basis DF Metric", naux, naux));  
    double** W = AOmetric->pointer(0);
    shared_ptr<BasisSet> zero = BasisSet::zero_ao_basis_set();

    // Only thread if not already in parallel (handy for local fitting)
    int nthread = 1;
    #ifdef _OPENMP
        if (!omp_in_parallel()) {
            nthread = omp_get_max_threads();
        }
    #endif

    // == (A|B) Block == //
    IntegralFactory rifactory_J(aux_, zero, aux_, zero);
    const double **Jbuffer = new const double*[nthread];
    shared_ptr<TwoBodyAOInt> *Jint = new shared_ptr<TwoBodyAOInt>[nthread];
    for (int Q = 0; Q<nthread; Q++) {
        Jint[Q] = shared_ptr<TwoBodyAOInt>(rifactory_J.eri());
        Jbuffer[Q] = Jint[Q]->buffer();
    }

    #pragma omp parallel for schedule (dynamic) num_threads(nthread) 
    for (int MU=0; MU < aux_->nshell(); ++MU) {
        int nummu = aux_->shell(MU)->nfunction();
    
        int thread = 0;
        #ifdef _OPENMP
            thread = omp_get_thread_num();
        #endif

        for (int NU=0; NU <= MU; ++NU) {
            int numnu = aux_->shell(NU)->nfunction();

            Jint[thread]->compute_shell(MU, 0, NU, 0);

            int index = 0;
            for (int mu=0; mu < nummu; ++mu) {
                int omu = aux_->shell(MU)->function_index() + mu;

                for (int nu=0; nu < numnu; ++nu, ++index) {
                    int onu = aux_->shell(NU)->function_index() + nu;

                    W[omu][onu] = Jbuffer[thread][index];
                    W[onu][omu] = Jbuffer[thread][index];
                }
            }
        }
    }
    delete[] Jbuffer; 
    delete[] Jint; 
    
    if (is_poisson_) {
        // == (AB) Block == //
        IntegralFactory rifactory_RP(pois_, aux_,  zero, zero);
        const double **Obuffer = new const double*[nthread];
        shared_ptr<OneBodyAOInt> *Oint = new shared_ptr<OneBodyAOInt>[nthread];
        for (int Q = 0; Q<nthread; Q++) {
            Oint[Q] = shared_ptr<OneBodyAOInt>(rifactory_RP.ao_overlap());
            Obuffer[Q] = Oint[Q]->buffer();
        }

        #pragma omp parallel for schedule (dynamic) num_threads(nthread) 
        for (int NU=0; NU < pois_->nshell(); ++NU) {
            int numnu = pois_->shell(NU)->nfunction();

            int thread = 0;
            #ifdef _OPENMP
                thread = omp_get_thread_num();
            #endif

            for (int MU=0; MU < aux_->nshell(); ++MU) {
                int nummu = aux_->shell(MU)->nfunction();

                Oint[thread]->compute_shell(NU, MU);

                int index = 0;
                for (int nu=0; nu < numnu; ++nu) {
                    int onu = pois_->shell(NU)->function_index() + nu;

                    for (int mu=0; mu < nummu; ++mu, ++index) {
                        int omu = aux_->shell(MU)->function_index() + mu;

                        W[omu][onu + ngaussian] = Obuffer[thread][index];
                        W[onu + ngaussian][omu] = Obuffer[thread][index];
                    }
                }
            }
        }

        // == (A|T|B) Block == //
        IntegralFactory rifactory_P(pois_, pois_,  zero, zero);
        const double **Tbuffer = new const double*[nthread];
        shared_ptr<OneBodyAOInt> *Tint = new shared_ptr<OneBodyAOInt>[nthread];
        for (int Q = 0; Q<nthread; Q++) {
            Tint[Q] = shared_ptr<OneBodyAOInt>(rifactory_P.ao_kinetic());
            Tbuffer[Q] = Tint[Q]->buffer();
        }

        #pragma omp parallel for schedule (dynamic) num_threads(nthread) 
        for (int MU=0; MU < pois_->nshell(); ++MU) {
            int nummu = pois_->shell(MU)->nfunction();

            int thread = 0;
            #ifdef _OPENMP
                thread = omp_get_thread_num();
            #endif

            for (int NU=0; NU <= MU ; ++NU) {
                int numnu = pois_->shell(NU)->nfunction();

                Tint[thread]->compute_shell(MU, NU);

                int index = 0;
                for (int mu=0; mu < nummu; ++mu) {
                    int omu = pois_->shell(MU)->function_index() + mu;

                    for (int nu=0; nu < numnu; ++nu, ++index) {
                        int onu = pois_->shell(NU)->function_index() + nu;

                        // These integrals are (A | -1/2 \nabla^2 | B), and should be (A | - 1 / (4 * pi) \nabla^2 |B)
                        // So a factor of 1 / (2 * PI) is the difference 
                        W[omu + ngaussian][onu + ngaussian] = 1.0 / (2.0 * M_PI) * Tbuffer[thread][index];
                        W[onu + ngaussian][omu + ngaussian] = 1.0 / (2.0 * M_PI) * Tbuffer[thread][index];
                    }
                }
            }
        }
        delete[] Tbuffer; 
        delete[] Tint; 
    }

    // If C1, form indexing and exit immediately
    if (auxpet->nirrep() == 1) {
        metric_ = AOmetric;
        pivots_ = shared_ptr<IntVector>(new IntVector(naux));
        int* piv = pivots_->pointer();
        for (int Q = 0; Q < naux; Q++)
            piv[Q] = Q;
        return;
    }    
}
void FittingMetric::form_cholesky_inverse()
{
    is_inverted_ = true;
    algorithm_ = "CHOLESKY";
}
void FittingMetric::form_QR_inverse(double tol)
{
    is_inverted_ = true;
    algorithm_ = "QR";
}
void FittingMetric::form_eig_inverse(double tol)
{
    is_inverted_ = true;
    algorithm_ = "EIG";
}
void FittingMetric::form_SVD_inverse(double tol)
{
    is_inverted_ = true;
    algorithm_ = "SVD";
}
void FittingMetric::form_full_inverse()
{
    is_inverted_ = true;
    algorithm_ = "FULL";
}

}
