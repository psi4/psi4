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

//MKL Header
#ifdef HAVE_MKL
#include <mkl.h>
#endif

//OpenMP Header
//_OPENMP is defined by the compiler if it exists
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace boost;
using namespace std;
using namespace psi;

namespace psi {

SchwarzSieve::SchwarzSieve(boost::shared_ptr<BasisSet> bas, double cut) :
    basis_(bas), schwarz_(cut), initialized_(false)
{
    form_schwarz_sieve(cut);
}
SchwarzSieve::~SchwarzSieve()
{
    if (initialized_) {
        free(schwarz_shells_);
        free(schwarz_funs_);
        free(schwarz_shells_reverse_);
        free(schwarz_funs_reverse_);
        free(schwarz_shell_vals_);
        free(schwarz_fun_vals_);
    }
}
void SchwarzSieve::form_schwarz_ints()
{
    int nshell = basis_->nshell();
    int nbf = basis_->nbf();

    schwarz_shell_vals_ = (double*) malloc(nshell*(nshell+1)/2L * sizeof(double));
    schwarz_fun_vals_ = (double*) malloc(nbf*(nbf+1)/2L * sizeof(double));

    max_global_val_ = 0.0;
    for (ULI Q = 0L; Q < nshell*(nshell+1)/2L; Q++)
        schwarz_shell_vals_[Q] = 0.0;
    for (ULI Q = 0L; Q < nbf*(nbf+1)/2L; Q++)
        schwarz_fun_vals_[Q] = 0.0;

    IntegralFactory schwarzfactory(basis_,basis_,basis_,basis_);
    boost::shared_ptr<TwoBodyAOInt> eri = boost::shared_ptr<TwoBodyAOInt>(schwarzfactory.eri());
    const double *buffer = eri->buffer();

    int MU, NU, mu, nu,omu,onu, nummu, numnu, index;
    ULI MUNU = 0L;
    ULI munu = 0L;
    for (MU=0; MU < nshell; ++MU) {
        nummu = basis_->shell(MU).nfunction();
        for (NU=0; NU <= MU; ++NU, ++MUNU) {
            numnu = basis_->shell(NU).nfunction();
            eri->compute_shell(MU,NU,MU,NU);
            for (mu=0; mu < nummu; ++mu) {
                omu = basis_->shell(MU).function_index() + mu;
                for (nu=0; nu < numnu; ++nu) {
                    onu = basis_->shell(NU).function_index() + nu;

                    if (omu>=onu) {
                        index = mu*(numnu*nummu*numnu+numnu)+nu*(nummu*numnu+1);
                        if (max_global_val_ < fabs(buffer[index]))
                            max_global_val_ = fabs(buffer[index]);
                        if (schwarz_shell_vals_[MUNU] < fabs(buffer[index]))
                            schwarz_shell_vals_[MUNU] = fabs(buffer[index]);
                        if (schwarz_fun_vals_[omu*(omu+1)/2+onu] < fabs(buffer[index]))
                            schwarz_fun_vals_[omu*(omu+1)/2+onu] = fabs(buffer[index]);
                    }
                }
            }
        }
    }
}
void SchwarzSieve::form_schwarz_sieve(double cut)
{
    schwarz_ = cut;
    int nshell = basis_->nshell();
    int nbf = basis_->nbf();

    if (!initialized_) {
        initialized_ = true;
        form_schwarz_ints();
        schwarz_shells_ = (int*) malloc(nshell * (nshell + 1L) * sizeof(int));
        schwarz_funs_ = (int*) malloc(nbf * (nbf + 1L) * sizeof(int));
        schwarz_shells_reverse_ = (long int*) malloc(nshell * (nshell + 1L) / 2L * sizeof(long int));
        schwarz_funs_reverse_ = (long int*) malloc(nbf * (nbf + 1L) / 2L * sizeof(long int));
    }

    nshell_pairs_ = 0L;
    nfun_pairs_ = 0L;

    double tol = (max_global_val_ > 0.0 ? schwarz_ * schwarz_ / max_global_val_ : 0.0);

    for (ULI Q = 0L; Q < nshell*(nshell+1)/2L; Q++) {
        if (schwarz_shell_vals_[Q] >= tol)
            nshell_pairs_++;
    }
    for (ULI Q = 0L; Q < nbf*(nbf+1)/2L; Q++) {
        if (schwarz_fun_vals_[Q] >= tol)
            nfun_pairs_++;
    }


    for (ULI Q = 0L; Q < nshell*(nshell+1)/2L; Q++)
        schwarz_shells_reverse_[Q] = -1L;

    for (ULI Q = 0L; Q < nbf*(nbf+1)/2L; Q++)
        schwarz_funs_reverse_[Q] = -1L;

    ULI counter = 0L; int MU;
    for (ULI Q = 0L, MU = 0; MU < nshell; MU++) {
        for (int NU = 0; NU <= MU;  NU++, Q++) {
            if (schwarz_shell_vals_[Q] >= tol) {
                schwarz_shells_[2*counter] = MU;
                schwarz_shells_[2*counter + 1] = NU;
                schwarz_shells_reverse_[Q] = counter;
                counter++;
            }
        }
    }
    counter = 0L;
    for (ULI Q = 0L, MU = 0; MU < nbf; MU++) {
        for (int NU = 0; NU <= MU;  NU++, Q++) {
            if (schwarz_fun_vals_[Q] >= tol) {
                schwarz_funs_[2*counter] = MU;
                schwarz_funs_[2*counter + 1] = NU;
                schwarz_funs_reverse_[Q] = counter;
                counter++;
            }
        }
    }
}

}

