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

#include "3index.h"

#include "psi4/libqt/qt.h"

#include <string>
#include <sstream>
#include <iostream>
#include <cstdio>
#include <fstream>
#include <algorithm>
#include <utility>
#include <ctype.h>
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/mintshelper.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/twobody.h"
#include "psi4/libmints/integral.h"
#include "psi4/libmints/basisset_parser.h"

using ULI=unsigned long int;

using namespace std;

namespace psi {

DFTensor::DFTensor(std::shared_ptr<BasisSet> primary,
                   std::shared_ptr<BasisSet> auxiliary,
                   SharedMatrix C,
                   int nocc,
                   int nvir,
                   int naocc,
                   int navir,
                   Options& options) :
    primary_(primary), auxiliary_(auxiliary), C_(C), nocc_(nocc), nvir_(nvir),
    naocc_(naocc), navir_(navir), options_(options)
{
    common_init();
}
DFTensor::DFTensor(std::shared_ptr<BasisSet> primary,
                   std::shared_ptr<BasisSet> auxiliary,
                   SharedMatrix C,
                   int nocc,
                   int nvir) :
    primary_(primary), auxiliary_(auxiliary), C_(C), nocc_(nocc), nvir_(nvir),
    naocc_(nocc), navir_(nvir), options_(Process::environment.options)
{
    common_init();
}
DFTensor::~DFTensor()
{
}
void DFTensor::common_init()
{
    print_ = options_.get_int("PRINT");
    debug_ = options_.get_int("DEBUG");

    print_header();

    molecule_ = primary_->molecule();

    nfocc_ = nocc_ - naocc_;
    nfvir_ = nvir_ - navir_;

    nso_ = C_->rowspi()[0];
    nmo_ = C_->colspi()[0];

    Caocc_ = SharedMatrix(new Matrix("C active occupied", nso_, naocc_));
    Cavir_ = SharedMatrix(new Matrix("C active virtual", nso_, navir_));

    double** Cp = C_->pointer();
    double** Cop  = Caocc_->pointer();
    double** Cvp  = Cavir_->pointer();

    for (int m = 0; m < nso_; m++) {
        C_DCOPY(naocc_, &Cp[m][nfocc_],1, Cop[m], 1);
        C_DCOPY(navir_, &Cp[m][nocc_],1, Cvp[m], 1);
    }

    if (debug_) {
        C_->print();
        Caocc_->print();
        Cavir_->print();
    }

    naux_ = auxiliary_->nbf();

    build_metric();
}
void DFTensor::print_header()
{
    outfile->Printf("  ==> DF Tensor (by Rob Parrish) <==\n\n");

    outfile->Printf(" => Primary Basis Set <= \n\n");
    primary_->print_by_level("outfile",print_);


    outfile->Printf(" => Auxiliary Basis Set <= \n\n");
    auxiliary_->print_by_level("outfile",print_);

}
void DFTensor::build_metric()
{
    std::shared_ptr<FittingMetric> met(new FittingMetric(auxiliary_, true));
    met->form_eig_inverse();
    metric_ = met->get_metric();

    if (debug_) {
        metric_->print();
    }
}
SharedMatrix DFTensor::Qso()
{
    SharedMatrix B(new Matrix("Bso", naux_, nso_ * nso_));
    SharedMatrix A(new Matrix("Aso", naux_, nso_ * nso_));
    double** Ap = A->pointer();
    double** Bp = B->pointer();
    double** Jp = metric_->pointer();

    std::shared_ptr<BasisSet> zero = BasisSet::zero_ao_basis_set();

    std::shared_ptr<IntegralFactory> fact(new IntegralFactory(auxiliary_,zero,primary_,primary_));
    std::shared_ptr<TwoBodyAOInt> eri(fact->eri());
    const double* buffer = eri->buffer();

    for (int P = 0; P < auxiliary_->nshell(); P++) {
        int np = auxiliary_->shell(P).nfunction();
        int pstart = auxiliary_->shell(P).function_index();
        for (int M = 0; M < primary_->nshell(); M++) {
            int nm = primary_->shell(M).nfunction();
            int mstart = primary_->shell(M).function_index();
            for (int N = 0; N < primary_->nshell(); N++) {
                int nn = primary_->shell(N).nfunction();
                int nstart = primary_->shell(N).function_index();

                eri->compute_shell(P,0,M,N);

                for (int p = 0, index = 0; p < np; p++) {
                    for (int m = 0; m < nm; m++) {
                        for (int n = 0; n < nn; n++, index++) {
                            Bp[p + pstart][(m + mstart) * nso_ + (n + nstart)] = buffer[index];
                        }
                    }
                }
            }
        }
    }

    C_DGEMM('N','N',naux_, nso_ * nso_, naux_, 1.0, Jp[0], naux_, Bp[0], nso_ * nso_, 0.0,
        Ap[0], nso_ * nso_);

    if (debug_) {
        metric_->print();
        B->print();
        A->print();
    }
    // Build numpy and final matrix shape
    std::vector<int> nshape{naux_, nso_, nso_};
    A->set_numpy_shape(nshape);

    return A;
}
SharedMatrix DFTensor::Qoo()
{
    SharedMatrix Amn = Qso();
    SharedMatrix Ami(new Matrix("Ami", naux_, naocc_ * (ULI) nso_));

    double** Amnp = Amn->pointer();
    double** Amip = Ami->pointer();
    double** Cop = Caocc_->pointer();

    C_DGEMM('N','N', naux_ * (ULI) nso_, naocc_, nso_, 1.0, Amnp[0], nso_, Cop[0], naocc_,
        0.0, Amip[0], naocc_);

    Amn.reset();

    SharedMatrix Aia(new Matrix("Aij", naux_, naocc_ * (ULI) naocc_));
    double** Aiap = Aia->pointer();

    for (int Q = 0; Q < naux_; Q++) {
        C_DGEMM('T','N',naocc_,naocc_,nso_,1.0,Amip[Q],naocc_,Cop[0],naocc_, 0.0, Aiap[Q], naocc_);
    }

    if (debug_) {
        Caocc_->print();
        Ami->print();
        Aia->print();
    }
    // Build numpy and final matrix shape
    std::vector<int> nshape{naux_, naocc_, naocc_};
    Aia->set_numpy_shape(nshape);

    return Aia;
}
SharedMatrix DFTensor::Qov()
{
    SharedMatrix Amn = Qso();
    SharedMatrix Ami(new Matrix("Qmi", naux_, naocc_ * (ULI) nso_));

    double** Amnp = Amn->pointer();
    double** Amip = Ami->pointer();
    double** Cop = Caocc_->pointer();
    double** Cvp = Cavir_->pointer();

    C_DGEMM('N','N', naux_ * (ULI) nso_, naocc_, nso_, 1.0, Amnp[0], nso_, Cop[0], naocc_,
        0.0, Amip[0], naocc_);

    Amn.reset();

    outfile->Printf( "DFTensor::Qov: naux %d, naocc %d, navir %d\n", naux_, naocc_, navir_);
    SharedMatrix Aia(new Matrix("Qia", naux_, naocc_ * (ULI) navir_));
    double** Aiap = Aia->pointer();

    for (int Q = 0; Q < naux_; Q++) {
        C_DGEMM('T','N',naocc_,navir_,nso_,1.0,Amip[Q],naocc_,Cvp[0],navir_, 0.0, Aiap[Q], navir_);
    }

    if (debug_) {
        Caocc_->print();
        Cavir_->print();
        Ami->print();
        Aia->print();
    }
    // Build numpy and final matrix shape
    std::vector<int> nshape{naux_, naocc_, navir_};
    Aia->set_numpy_shape(nshape);

    return Aia;
}
SharedMatrix DFTensor::Qvv()
{
    SharedMatrix Amn = Qso();
    SharedMatrix Ami(new Matrix("Qmi", naux_, navir_ * (ULI) nso_));

    double** Amnp = Amn->pointer();
    double** Amip = Ami->pointer();
    double** Cvp = Cavir_->pointer();

    C_DGEMM('N','N', naux_ * (ULI) nso_, navir_, nso_, 1.0, Amnp[0], nso_, Cvp[0], navir_,
        0.0, Amip[0], navir_);

    Amn.reset();

    SharedMatrix Aia(new Matrix("Qab", naux_, navir_ * (ULI) navir_));
    double** Aiap = Aia->pointer();

    for (int Q = 0; Q < naux_; Q++) {
        C_DGEMM('T','N',navir_,navir_,nso_,1.0,Amip[Q],navir_,Cvp[0],navir_, 0.0, Aiap[Q], navir_);
    }

    if (debug_) {
        Cavir_->print();
        Ami->print();
        Aia->print();
    }
    // Build numpy and final matrix shape
    std::vector<int> nshape{naux_, navir_, navir_};
    Aia->set_numpy_shape(nshape);

    return Aia;
}
SharedMatrix DFTensor::Qmo()
{
    SharedMatrix Amn = Qso();
    SharedMatrix Ami(new Matrix("Qmi", naux_, nmo_ * (ULI) nso_));

    double** Amnp = Amn->pointer();
    double** Amip = Ami->pointer();
    double** Cvp = C_->pointer();

    C_DGEMM('N','N', naux_ * (ULI) nso_, nmo_, nso_, 1.0, Amnp[0], nso_, Cvp[0], nmo_,
        0.0, Amip[0], nmo_);

    Amn.reset();

    SharedMatrix Aia(new Matrix("Qmo", naux_, nmo_ * (ULI) nmo_));
    double** Aiap = Aia->pointer();

    for (int Q = 0; Q < naux_; Q++) {
        C_DGEMM('T','N',nmo_,nmo_,nso_,1.0,Amip[Q],nmo_,Cvp[0],nmo_, 0.0, Aiap[Q], nmo_);
    }

    if (debug_) {
        C_->print();
        Ami->print();
        Aia->print();
    }
    // Build numpy and final matrix shape
    std::vector<int> nshape{naux_, nmo_, nmo_};
    Aia->set_numpy_shape(nshape);

    return Aia;
}
SharedMatrix DFTensor::Imo()
{
    std::shared_ptr<MintsHelper> mints(new MintsHelper(primary_, options_, 0));
    return mints->mo_eri(C_,C_);
}
SharedMatrix DFTensor::Idfmo()
{
    SharedMatrix Amo = Qmo();
    double** Amop = Amo->pointer();

    SharedMatrix Imo(new Matrix("DF MO ERI Tensor", nmo_ * nmo_, nmo_ * nmo_));
    double** Imop = Imo->pointer();

    C_DGEMM('T','N',nmo_ * nmo_, nmo_ * nmo_, naux_, 1.0, Amop[0], nmo_ * nmo_,
        Amop[0], nmo_ * nmo_, 0.0, Imop[0], nmo_ * nmo_);

    // Build numpy and final matrix shape
    std::vector<int> nshape{nmo_, nmo_, nmo_, nmo_};
    Imo->set_numpy_shape(nshape);

    return Imo;
}

}
