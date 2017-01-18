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

 #include "psi4/pragma.h"
 PRAGMA_WARNING_PUSH
 PRAGMA_WARNING_IGNORE_DEPRECATED_DECLARATIONS
 #include <memory>
 PRAGMA_WARNING_POP

#include "psi4/libqt/qt.h"
#include <math.h>
#include "fitter.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/integral.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/vector.h"

namespace psi {

DFChargeFitter::DFChargeFitter() :
    print_(0), debug_(0)
{
}
DFChargeFitter::~DFChargeFitter()
{
}
SharedVector DFChargeFitter::fit()
{
    int naux = auxiliary_->nbf();
    int nso  = primary_->nbf();

    SharedVector d (new Vector("d", naux));

    double* dp = d->pointer();
    double** Dp = D_->pointer();

    /* 3-index */ {

    std::shared_ptr<IntegralFactory> factory(new IntegralFactory(auxiliary_,BasisSet::zero_ao_basis_set(),
        primary_,primary_));
    std::shared_ptr<TwoBodyAOInt> eri(factory->eri());
    const double* buffer = eri->buffer();

    for (int Q = 0; Q < auxiliary_->nshell(); Q++) {
        for (int M = 0; M < primary_->nshell(); M++) {
            for (int N = 0; N < primary_->nshell(); N++) {
                eri->compute_shell(Q,0,M,N);
                int nq = auxiliary_->shell(Q).nfunction();
                int nm = primary_->shell(M).nfunction();
                int nn = primary_->shell(N).nfunction();
                int sq = auxiliary_->shell(Q).function_index();
                int sm = primary_->shell(M).function_index();
                int sn = primary_->shell(N).function_index();
                for (int oq = 0; oq < nq; oq++) {
                    for (int om = 0; om < nm; om++) {
                        for (int on = 0; on < nn; on++) {
                            dp[sq + oq] += Dp[sm + om][sn + on] * buffer[oq * nm * nn + om * nn + on];
                        }
                    }
                }
            }
        }
    }

    /* End 3-index */ }
    /* 2-index */ {

    SharedMatrix J(new Matrix("J", naux, naux));
    double** Jp = J->pointer();

    std::shared_ptr<IntegralFactory> factory(new IntegralFactory(auxiliary_,BasisSet::zero_ao_basis_set(),
                                                                   auxiliary_,BasisSet::zero_ao_basis_set()));
    std::shared_ptr<TwoBodyAOInt> eri(factory->eri());
    const double* buffer = eri->buffer();

    for (int Q = 0; Q < auxiliary_->nshell(); Q++) {
        for (int P = 0; P < auxiliary_->nshell(); P++) {
            eri->compute_shell(Q,0,P,0);
            int nq = auxiliary_->shell(Q).nfunction();
            int np = auxiliary_->shell(P).nfunction();
            int sq = auxiliary_->shell(Q).function_index();
            int sp = auxiliary_->shell(P).function_index();
            for (int oq = 0; oq < nq; oq++) {
                for (int op = 0; op < np; op++) {
                    Jp[sq + oq][sp + op] = buffer[oq * np + op];
                }
            }
        }
    }

    int info;
    info = C_DPOTRF('L',naux,Jp[0],naux);
    if (info) throw PSIEXCEPTION("DFChargeFitter: C_DPOTRF Failed");
    info = C_DPOTRS('L',naux,1,Jp[0],naux,dp,naux);
    if (info) throw PSIEXCEPTION("DFChargeFitter: C_DPOTRS Failed");

    /* End 2-index */ }

    d_ = d;
    return d;
}


} // Namespace psi
