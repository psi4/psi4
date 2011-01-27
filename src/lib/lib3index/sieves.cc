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


using namespace std;
using namespace psi;

namespace psi { 

void ThreeIndexTensor::form_schwarz_ints()
{
    int nshell = primary_basis_->nshell();

    schwarz_shell_vals_ = (double*) malloc(nshell*(nshell+1)/2L * sizeof(double)); 
    schwarz_fun_vals_ = (double*) malloc(nbf_*(nbf_+1)/2L * sizeof(double)); 

    max_global_val_ = 0.0;
    for (ULI Q = 0L; Q < nshell*(nshell+1)/2L; Q++)
        schwarz_shell_vals_[Q] = 0.0;   
    for (ULI Q = 0L; Q < nbf_*(nbf_+1)/2L; Q++)
        schwarz_fun_vals_[Q] = 0.0;   

    if (schwarz_cutoff_ > 0.0) {
 
        IntegralFactory schwarzfactory(primary_basis_,primary_basis_,primary_basis_,primary_basis_);
        shared_ptr<TwoBodyAOInt> eri = shared_ptr<TwoBodyAOInt>(schwarzfactory.eri());
        const double *buffer = eri->buffer();
        
        int MU, NU, mu, nu,omu,onu, nummu, numnu, index;
        ULI MUNU = 0L;
        ULI munu = 0L;
        for (MU=0; MU < nshell; ++MU) {
            nummu = primary_basis_->shell(MU)->nfunction();
            for (NU=0; NU <= MU; ++NU, ++MUNU) {
                numnu = primary_basis_->shell(NU)->nfunction();
                eri->compute_shell(MU,NU,MU,NU);
                for (mu=0; mu < nummu; ++mu) {
                    omu = primary_basis_->shell(MU)->function_index() + mu;
                    for (nu=0; nu < numnu; ++nu) {
                        onu = primary_basis_->shell(NU)->function_index() + nu;
        
                        if (omu>=onu) {
                            index = mu*(numnu*nummu*numnu+numnu)+nu*(nummu*numnu+1);
                            if (max_global_val_<abs(buffer[index]))
                                max_global_val_ = abs(buffer[index]);
                            if (schwarz_shell_vals_[MUNU]<abs(buffer[index]))
                                schwarz_shell_vals_[MUNU] = abs(buffer[index]);
                            if (schwarz_fun_vals_[omu*(omu+1)/2+onu]<abs(buffer[index]))
                                schwarz_fun_vals_[omu*(omu+1)/2+onu] = abs(buffer[index]);
                        }
                    }
                }
            }
        }
    }

}
void ThreeIndexTensor::form_schwarz_sieve(double cut)
{
    schwarz_cutoff_ = cut;
    int nshell = primary_basis_->nshell();

    if (schwarz_fun_vals_ == NULL) {
        form_schwarz_ints();
        schwarz_shells_ = (int*) malloc(nshell * (nshell + 1L) * sizeof(int)); 
        schwarz_funs_ = (int*) malloc(nbf_ * (nbf_ + 1L) * sizeof(int)); 
        schwarz_shells_reverse_ = (long int*) malloc(nshell * (nshell + 1L) / 2L * sizeof(long int));
        schwarz_funs_reverse_ = (long int*) malloc(nbf_ * (nbf_ + 1L) / 2L * sizeof(long int));
    }

    nshell_pairs_ = 0L;
    nfun_pairs_ = 0L;

    double tol = (max_global_val_ > 0.0 ? schwarz_cutoff_ * schwarz_cutoff_ / max_global_val_ : 0.0);

    for (ULI Q = 0L; Q < nshell*(nshell+1)/2L; Q++) {
        if (schwarz_shell_vals_[Q] >= tol)
            nshell_pairs_++; 
    }
    for (ULI Q = 0L; Q < nbf_*(nbf_+1)/2L; Q++) {
        if (schwarz_fun_vals_[Q] >= tol)
            nfun_pairs_++; 
    }

    
    for (ULI Q = 0L; Q < nshell*(nshell+1)/2L; Q++)
        schwarz_shells_reverse_[Q] = -1L; 

    for (ULI Q = 0L; Q < nbf_*(nbf_+1)/2L; Q++)
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
    for (ULI Q = 0L, MU = 0; MU < nbf_; MU++) {
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

