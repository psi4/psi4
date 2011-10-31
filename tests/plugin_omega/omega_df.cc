#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <algorithm>
#include <vector>
#include <utility>
#include <string>
#include <cstring>

#include <psifiles.h>
#include <physconst.h>
#include <libciomr/libciomr.h>
#include <libpsio/psio.hpp>
#include <libchkpt/chkpt.hpp>
#include <libiwl/iwl.hpp>
#include <libqt/qt.h>

#include <libmints/mints.h>
#include <libfunctional/superfunctional.h>
#include <libscf_solver/ks.h>
#include <libscf_solver/integralfunctors.h>
#include <libscf_solver/omegafunctors.h>
#include <lib3index/3index.h>

#include "omega.h"

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;
using namespace psi;
using namespace psi::functional;
using namespace boost;

namespace psi{ namespace scf {

OmegaDF::OmegaDF(boost::shared_ptr<PSIO> psio, boost::shared_ptr<BasisSet> primary, boost::shared_ptr<BasisSet> auxiliary) :
    psio_(psio), primary_(primary), auxiliary_(auxiliary)
{
    common_init();
}
OmegaDF::~OmegaDF()
{
}
void OmegaDF::common_init()
{
    omega_ = 0.0;

    nthread_ = 1;
    #ifdef _OPENMP
        nthread_ = omp_get_max_threads();
    #endif

    zero_ = BasisSet::zero_ao_basis_set();

    build_static();
}
void OmegaDF::set_omega(double omega)
{
    omega_ = omega;

    for (int i = 0; i < nthread_; i++) {
        erf_eri_[i]->setOmega(omega_);
    }
    build_dynamic();
}
void OmegaDF::build_static()
{
    int nso = primary_->nbf();
    int naux = auxiliary_->nbf();

    Cmn_ = SharedMatrix(new Matrix("Cmn", naux, nso * (ULI) nso)); 
    Amn_ = SharedMatrix(new Matrix("Amn", naux, nso * (ULI) nso)); 
    Wmn_ = SharedMatrix(new Matrix("Wmn", naux, nso * (ULI) nso)); 
    double** Cmnp = Cmn_->pointer();
    double** Wmnp = Wmn_->pointer();
    double** Amnp = Amn_->pointer();

    factory_ = boost::shared_ptr<IntegralFactory>(new IntegralFactory(auxiliary_, zero_, primary_, primary_));
    std::vector<boost::shared_ptr<TwoBodyAOInt> > eri;
    std::vector<const double*> buffer;

    for (int i = 0; i < nthread_; i++) {
        eri.push_back(boost::shared_ptr<TwoBodyAOInt>(factory_->eri()));
        buffer.push_back(eri[i]->buffer());
        erf_eri_.push_back(boost::shared_ptr<ErfERI>(static_cast<ErfERI*>(factory_->erf_eri(0.0))));    
        erf_buffer_.push_back(erf_eri_[i]->buffer());
    }

    #pragma omp parallel for num_threads(nthread_) schedule(dynamic)
    for (int P = 0; P < auxiliary_->nshell(); P++) {
        int thread = 0;
        #ifdef _OPENMP
            thread = omp_get_thread_num();
        #endif
        int nP = auxiliary_->shell(P)->nfunction();
        int p0 = auxiliary_->shell(P)->function_index();
        for (int M = 0; M < primary_->nshell(); M++) {
            int nM = primary_->shell(M)->nfunction();
            int m0 = primary_->shell(M)->function_index();
            for (int N = 0; N < primary_->nshell(); N++) {
                int nN = primary_->shell(N)->nfunction();
                int n0 = primary_->shell(N)->function_index();

                eri[thread]->compute_shell(P,0,M,N);

                for (int p = 0, index = 0; p < nP; p++) {
                    int op = p0 + p;
                    for (int m = 0; m < nM; m++) {
                        int om = m0 + m;
                        for (int n = 0; n < nN; n++, index++) {
                            int on = n0 + n;
                            Amnp[op][om * nso + on] = buffer[thread][index];
                        }    
                    }    
                }    
            }
        }
    }

    boost::shared_ptr<FittingMetric> fit(new FittingMetric(auxiliary_));
    fit->form_full_eig_inverse();
    SharedMatrix J = fit->get_metric();
    double** Jp = J->pointer();

    C_DGEMM('N','N',naux, nso * (ULI) nso, naux, 1.0, Jp[0], naux, Amnp[0], nso * (ULI) nso, 0.0,
        Cmnp[0], nso * (ULI) nso);

}
void OmegaDF::build_dynamic()
{
    int nso = primary_->nbf();
    double** Wmnp = Wmn_->pointer();
    #pragma omp parallel for num_threads(nthread_) schedule(dynamic)
    for (int P = 0; P < auxiliary_->nshell(); P++) {
        int thread = 0;
        #ifdef _OPENMP
            thread = omp_get_thread_num();
        #endif
        int nP = auxiliary_->shell(P)->nfunction();
        int p0 = auxiliary_->shell(P)->function_index();
        for (int M = 0; M < primary_->nshell(); M++) {
            int nM = primary_->shell(M)->nfunction();
            int m0 = primary_->shell(M)->function_index();
            for (int N = 0; N < primary_->nshell(); N++) {
                int nN = primary_->shell(N)->nfunction();
                int n0 = primary_->shell(N)->function_index();

                erf_eri_[thread]->compute_shell(P,0,M,N);

                for (int p = 0, index = 0; p < nP; p++) {
                    int op = p0 + p;
                    for (int m = 0; m < nM; m++) {
                        int om = m0 + m;
                        for (int n = 0; n < nN; n++, index++) {
                            int on = n0 + n;
                            Wmnp[op][om * nso + on] = erf_buffer_[thread][index];
                        }    
                    }    
                }    
            }
        }
    }
}
SharedMatrix OmegaDF::J(SharedMatrix D)
{
    int nso = primary_->nbf();
    int naux = auxiliary_->nbf();
    
    SharedMatrix J(new Matrix("J",nso,nso));
    double* d = new double[naux];

    double** Jp = J->pointer();
    double** Dp = D->pointer();

    double** Cmnp = Cmn_->pointer();
    double** Amnp = Amn_->pointer();

    C_DGEMV('N', naux, nso *(ULI) nso, 1.0, Cmnp[0], nso * (ULI) nso, Dp[0], 1, 0.0, d, 1);
    C_DGEMV('T', naux, nso *(ULI) nso, 1.0, Amnp[0], nso * (ULI) nso, d, 1, 1.0, Jp[0], 1);

    return J;
} 
SharedMatrix OmegaDF::K(SharedMatrix C, int nocc)
{
    int nso = primary_->nbf();
    int naux = auxiliary_->nbf();
    int nmo = C->colspi()[0];   
 
    SharedMatrix K(new Matrix("K",nso,nso));

    if (nocc == 0) return K;

    SharedMatrix W(new Matrix("W",naux, nocc * (ULI) nso));
    SharedMatrix Q(new Matrix("Q",naux, nocc * (ULI) nso));
    SharedMatrix T(new Matrix("T",nso, nocc));

    double** Wmnp = Amn_->pointer(); 
    double** Cmnp = Cmn_->pointer(); 
    double** Wp = W->pointer(); 
    double** Qp = Q->pointer(); 
    double** Tp = T->pointer(); 
    double** Kp = K->pointer(); 
    double** Cp = C->pointer(); 

    C_DGEMM('N','N',naux * (ULI) nso, nocc, nso, 1.0, Wmnp[0], nso, Cp[0], nmo, 0.0, Wp[0], nocc);
    C_DGEMM('N','N',naux * (ULI) nso, nocc, nso, 1.0, Cmnp[0], nso, Cp[0], nmo, 0.0, Qp[0], nocc);

    for (int P = 0; P < naux; P++) {
        C_DCOPY(nocc *(ULI) nso, Wp[P], 1, Tp[0], 1);
        for (int i = 0; i < nocc; i++) {
            C_DCOPY(nso, &Tp[0][i], nocc, &Wp[P][i * nso], 1);
        }
        C_DCOPY(nocc *(ULI) nso, Qp[P], 1, Tp[0], 1);
        for (int i = 0; i < nocc; i++) {
            C_DCOPY(nso, &Tp[0][i], nocc, &Qp[P][i * nso], 1);
        }
    }

    C_DGEMM('T','N',nso, nso, nocc * (ULI)naux, 1.0, Wp[0], nso, Qp[0], nso, 0.0, Kp[0], nso);

    return K; 
}
SharedMatrix OmegaDF::wK(SharedMatrix C, int nocc)
{
    int nso = primary_->nbf();
    int naux = auxiliary_->nbf();
    int nmo = C->colspi()[0];   
 
    SharedMatrix K(new Matrix("wK",nso,nso));

    if (nocc == 0) return K;

    SharedMatrix W(new Matrix("W",naux, nocc * (ULI) nso));
    SharedMatrix Q(new Matrix("Q",naux, nocc * (ULI) nso));
    SharedMatrix T(new Matrix("T",nso, nocc));

    double** Wmnp = Wmn_->pointer(); 
    double** Cmnp = Cmn_->pointer(); 
    double** Wp = W->pointer(); 
    double** Qp = Q->pointer(); 
    double** Tp = T->pointer(); 
    double** Kp = K->pointer(); 
    double** Cp = C->pointer(); 

    C_DGEMM('N','N',naux * (ULI) nso, nocc, nso, 1.0, Wmnp[0], nso, Cp[0], nmo, 0.0, Wp[0], nocc);
    C_DGEMM('N','N',naux * (ULI) nso, nocc, nso, 1.0, Cmnp[0], nso, Cp[0], nmo, 0.0, Qp[0], nocc);

    for (int P = 0; P < naux; P++) {
        C_DCOPY(nocc *(ULI) nso, Wp[P], 1, Tp[0], 1);
        for (int i = 0; i < nocc; i++) {
            C_DCOPY(nso, &Tp[0][i], nocc, &Wp[P][i * nso], 1);
        }
        C_DCOPY(nocc *(ULI) nso, Qp[P], 1, Tp[0], 1);
        for (int i = 0; i < nocc; i++) {
            C_DCOPY(nso, &Tp[0][i], nocc, &Qp[P][i * nso], 1);
        }
    }

    C_DGEMM('T','N',nso, nso, nocc * (ULI)naux, 1.0, Wp[0], nso, Qp[0], nso, 0.0, Kp[0], nso);

    return K; 
}

}} // End Namespaces
