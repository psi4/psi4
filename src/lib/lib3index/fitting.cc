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

shared_ptr<Matrix> DFTensor::form_fitting_metric()
{
    shared_ptr<Matrix> fitting_metric = shared_ptr<Matrix>(new Matrix("DF Metric", naux_, naux_)); 
    double** W = fitting_metric->get_pointer(0);

    // == (A|B) Block == //
    IntegralFactory rifactory_J(auxiliary_basis_, zero_, auxiliary_basis_, zero_);
    const double **Jbuffer = new const double*[nthread_];
    shared_ptr<TwoBodyAOInt> *Jint = new shared_ptr<TwoBodyAOInt>[nthread_];
    for (int Q = 0; Q<nthread_; Q++) {
        Jint[Q] = shared_ptr<TwoBodyAOInt>(rifactory_J.eri());
        Jbuffer[Q] = Jint[Q]->buffer();
    }

    #pragma omp parallel for schedule (dynamic) num_threads(nthread_) 
    for (int MU=0; MU < auxiliary_basis_->nshell(); ++MU) {
        int nummu = auxiliary_basis_->shell(MU)->nfunction();
    
        int thread = 0;
        #ifdef _OPENMP
            thread = omp_get_thread_num();
        #endif

        for (int NU=0; NU <= MU; ++NU) {
            int numnu = auxiliary_basis_->shell(NU)->nfunction();

            Jint[thread]->compute_shell(MU, 0, NU, 0);

            int index = 0;
            for (int mu=0; mu < nummu; ++mu) {
                int omu = auxiliary_basis_->shell(MU)->function_index() + mu;

                for (int nu=0; nu < numnu; ++nu, ++index) {
                    int onu = auxiliary_basis_->shell(NU)->function_index() + nu;

                    W[omu][onu] = Jbuffer[thread][index];
                    W[onu][omu] = Jbuffer[thread][index];
                }
            }
        }
    }
    delete[] Jbuffer; 
    delete[] Jint; 
    
    if (poisson_) {
        // == (AB) Block == //
        IntegralFactory rifactory_RP(poisson_basis_, auxiliary_basis_,  zero_, zero_);
        const double **Obuffer = new const double*[nthread_];
        shared_ptr<OneBodyAOInt> *Oint = new shared_ptr<OneBodyAOInt>[nthread_];
        for (int Q = 0; Q<nthread_; Q++) {
            Oint[Q] = shared_ptr<OneBodyAOInt>(rifactory_RP.ao_overlap());
            Obuffer[Q] = Oint[Q]->buffer();
        }

        #pragma omp parallel for schedule (dynamic) num_threads(nthread_) 
        for (int NU=0; NU < poisson_basis_->nshell(); ++NU) {
            int numnu = poisson_basis_->shell(NU)->nfunction();

            int thread = 0;
            #ifdef _OPENMP
                thread = omp_get_thread_num();
            #endif

            for (int MU=0; MU < auxiliary_basis_->nshell(); ++MU) {
                int nummu = auxiliary_basis_->shell(MU)->nfunction();

                Oint[thread]->compute_shell(NU, MU);

                int index = 0;
                for (int nu=0; nu < numnu; ++nu) {
                    int onu = poisson_basis_->shell(NU)->function_index() + nu;

                    for (int mu=0; mu < nummu; ++mu, ++index) {
                        int omu = auxiliary_basis_->shell(MU)->function_index() + mu;

                        W[omu][onu + ngaussian_] = Obuffer[thread][index];
                        W[onu + ngaussian_][omu] = Obuffer[thread][index];
                    }
                }
            }
        }

        // == (A|T|B) Block == //
        IntegralFactory rifactory_P(poisson_basis_, poisson_basis_,  zero_, zero_);
        const double **Tbuffer = new const double*[nthread_];
        shared_ptr<OneBodyAOInt> *Tint = new shared_ptr<OneBodyAOInt>[nthread_];
        for (int Q = 0; Q<nthread_; Q++) {
            Tint[Q] = shared_ptr<OneBodyAOInt>(rifactory_P.ao_kinetic());
            Tbuffer[Q] = Tint[Q]->buffer();
        }

        #pragma omp parallel for schedule (dynamic) num_threads(nthread_) 
        for (int MU=0; MU < poisson_basis_->nshell(); ++MU) {
            int nummu = poisson_basis_->shell(MU)->nfunction();

            int thread = 0;
            #ifdef _OPENMP
                thread = omp_get_thread_num();
            #endif

            for (int NU=0; NU <= MU ; ++NU) {
                int numnu = poisson_basis_->shell(NU)->nfunction();

                Tint[thread]->compute_shell(MU, NU);

                int index = 0;
                for (int mu=0; mu < nummu; ++mu) {
                    int omu = poisson_basis_->shell(MU)->function_index() + mu;

                    for (int nu=0; nu < numnu; ++nu, ++index) {
                        int onu = poisson_basis_->shell(NU)->function_index() + nu;

                        // These integrals are (A | -1/2 \nabla^2 | B), and should be (A | - 1 / (4 * pi) \nabla^2 |B)
                        // So a factor of 1 / (2 * PI) is the difference 
                        W[omu + ngaussian_][onu + ngaussian_] = 1.0 / (2.0 * M_PI) * Tbuffer[thread][index];
                        W[onu + ngaussian_][omu + ngaussian_] = 1.0 / (2.0 * M_PI) * Tbuffer[thread][index];
                    }
                }
            }
        }
        delete[] Tbuffer; 
        delete[] Tint; 
    }
    fitting_metric->save(psio_, PSIF_3INDEX);
    return fitting_metric;
}
shared_ptr<Matrix> DFTensor::form_cholesky_metric()
{
    // Grab the fitting metric
    shared_ptr<Matrix> M = form_fitting_metric();
    for (int h = 0; h < M->nirreps(); h++) {
        // Choleskify each irrep O(N^3/3), essentially free
        PSI_DPOTRF(h,'L', M->colspi()[h], M, M->colspi()[h]);
        // The result goes in the upper triangle 
        // Zero the lower triangle for stype 
        double** Mp = M->get_pointer(h);
        for (int A = 0; A < M->colspi()[h]; A++)
            for (int B = 0; B < A; B++)
                Mp[A][B] = 0.0; 
    }
    M->set_name("DF Metric Cholesky");
    M->save(psio_, PSIF_3INDEX);
    return M;
}
shared_ptr<Matrix> DFTensor::form_qr_metric(double cutoff)
{
    // Grab the fitting metric
    shared_ptr<Matrix> M = form_fitting_metric();
    M->set_name("DF Metric QR");
    return M;
}


}
