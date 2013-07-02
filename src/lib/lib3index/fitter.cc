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

#include <boost/shared_ptr.hpp>
#include <libmints/mints.h>
#include <libqt/qt.h>
#include <math.h>
#include "fitter.h"

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

    boost::shared_ptr<IntegralFactory> factory(new IntegralFactory(auxiliary_,BasisSet::zero_ao_basis_set(),
        primary_,primary_));
    boost::shared_ptr<TwoBodyAOInt> eri(factory->eri());
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

    boost::shared_ptr<IntegralFactory> factory(new IntegralFactory(auxiliary_,BasisSet::zero_ao_basis_set(),
                                                                   auxiliary_,BasisSet::zero_ao_basis_set()));
    boost::shared_ptr<TwoBodyAOInt> eri(factory->eri());
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
