/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2022 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This file is part of Psi4.
 *
 * Psi4 is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * Psi4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with Psi4; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

#include "points.h"
#include "hamiltonian.h"
#include "jk.h"
#include "v.h"

#include "psi4/libqt/qt.h"
#include "psi4/psi4-dec.h"
#include "psi4/libmints/vector.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libpsi4util/PsiOutStream.h"

#include <sstream>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace psi {

Hamiltonian::Hamiltonian(std::shared_ptr<JK> jk) : jk_(jk) { common_init(); }
Hamiltonian::Hamiltonian(std::shared_ptr<JK> jk, std::shared_ptr<VBase> v) : jk_(jk), v_(v) { common_init(); }
Hamiltonian::~Hamiltonian() {}
void Hamiltonian::common_init() {
    print_ = 1;
    debug_ = 0;
    bench_ = 0;
    exact_diagonal_ = false;
}

RHamiltonian::RHamiltonian(std::shared_ptr<JK> jk) : Hamiltonian(jk) {}
RHamiltonian::RHamiltonian(std::shared_ptr<JK> jk, std::shared_ptr<VBase> v) : Hamiltonian(jk, v) {}
RHamiltonian::~RHamiltonian() {}
SharedMatrix RHamiltonian::explicit_hamiltonian() {
    std::shared_ptr<Vector> diag = diagonal();

    auto H = std::make_shared<Matrix>("Explicit Hamiltonian", diag->nirrep(), diag->dimpi(), diag->dimpi());

    auto b = std::make_shared<Vector>(std::move(diag->clone()));
    auto s = std::make_shared<Vector>(std::move(diag->clone()));
    std::vector<std::shared_ptr<Vector> > bb;
    std::vector<std::shared_ptr<Vector> > ss;
    bb.push_back(b);
    ss.push_back(s);
    for (int h = 0; h < diag->nirrep(); h++) {
        for (int n = 0; n < diag->dimpi()[h]; n++) {
            b->zero();
            s->zero();
            b->set(h, n, 1.0);
            product(bb, ss);
            C_DCOPY(diag->dimpi()[h], s->pointer(h), 1, H->pointer(h)[n], 1);
        }
    }

    return H;
}

CPHFRHamiltonian::CPHFRHamiltonian(std::shared_ptr<JK> jk, SharedMatrix Caocc, SharedMatrix Cavir,
                                   std::shared_ptr<Vector> eps_aocc, std::shared_ptr<Vector> eps_avir,
                                   std::shared_ptr<VBase> v)
    : RHamiltonian(jk, v), Caocc_(Caocc), Cavir_(Cavir), eps_aocc_(eps_aocc), eps_avir_(eps_avir) {}
CPHFRHamiltonian::~CPHFRHamiltonian() {}
void CPHFRHamiltonian::print_header() const {
    if (print_) {
        outfile->Printf("  ==> CPHFRHamiltonian (by Rob Parrish) <== \n\n");
    }
}
std::shared_ptr<Vector> CPHFRHamiltonian::diagonal() {
    int nirrep = eps_aocc_->nirrep();
    Dimension nov(nirrep);
    for (int symm = 0; symm < nirrep; ++symm) {
        for (int h = 0; h < nirrep; ++h) {
            nov[symm] += eps_aocc_->dimpi()[h] * eps_avir_->dimpi()[symm ^ h];
        }
    }

    auto diag = std::make_shared<Vector>("CPHF Diagonal", nov);

    for (int symm = 0; symm < nirrep; ++symm) {
        long int offset = 0L;
        for (int h = 0; h < nirrep; ++h) {
            int nocc = eps_aocc_->dimpi()[h];
            int nvir = eps_avir_->dimpi()[symm ^ h];

            if (!nocc || !nvir) continue;

            double* eop = eps_aocc_->pointer(h);
            double* evp = eps_avir_->pointer(symm ^ h);
            double* dp = diag->pointer(symm);

            for (int i = 0; i < nocc; ++i) {
                for (int a = 0; a < nvir; ++a) {
                    dp[i * nvir + a + offset] = evp[a] - eop[i];
                }
            }
            offset += static_cast<long int>(nocc) * nvir;
        }
    }

    return diag;
}
std::map<std::string, SharedVector> CPHFRHamiltonian::pack(const std::map<std::string, SharedMatrix>& x) {
    int nirrep = eps_aocc_->nirrep();
    Dimension nov(nirrep);
    for (int symm = 0; symm < nirrep; ++symm) {
        for (int h = 0; h < nirrep; ++h) {
            nov[symm] += eps_aocc_->dimpi()[h] * eps_avir_->dimpi()[symm ^ h];
        }
    }

    std::map<std::string, SharedVector> X;
    for (std::map<std::string, SharedMatrix>::const_iterator it = x.begin(); it != x.end(); ++it) {
        auto v = std::make_shared<Vector>("X", nov);
        SharedMatrix x2 = (*it).second;
        int symm = x2->symmetry();
        int offset = 0;
        for (int h = 0; h < nirrep; h++) {
            int nocc = eps_aocc_->dimpi()[h];
            int nvir = eps_avir_->dimpi()[h ^ symm];

            if (!nocc || !nvir) continue;

            ::memcpy((void*)&v->pointer(symm)[offset], (void*)x2->pointer(h)[0], sizeof(double) * nocc * nvir);

            offset += nocc * nvir;
        }
        X[(*it).first] = v;
    }
    return X;
}
void CPHFRHamiltonian::product(const std::vector<std::shared_ptr<Vector> >& x,
                               std::vector<std::shared_ptr<Vector> >& b) {
    std::vector<SharedMatrix>& C_left = jk_->C_left();
    std::vector<SharedMatrix>& C_right = jk_->C_right();

    C_left.clear();
    C_right.clear();

    int nirrep = (x.size() ? x[0]->nirrep() : 0);

    for (int symm = 0; symm < nirrep; ++symm) {
        for (size_t N = 0; N < x.size(); ++N) {
            C_left.push_back(Caocc_);
            double* xp = x[N]->pointer(symm);

            std::stringstream ss;
            ss << "C_right, h = " << symm << ", N = " << N;
            auto Cr = std::make_shared<Matrix>(ss.str(), Caocc_->nirrep(), Caocc_->rowspi(), Caocc_->colspi(), symm);

            long int offset = 0L;
            for (int h = 0; h < Caocc_->nirrep(); ++h) {
                int nocc = Caocc_->colspi()[h];
                int nvir = Cavir_->colspi()[h ^ symm];
                int nso = Cavir_->rowspi()[h ^ symm];

                if (!nso || !nocc || !nvir) continue;

                double** Cvp = Cavir_->pointer(h ^ symm);
                double** Crp = Cr->pointer(h ^ symm);

                C_DGEMM('N', 'T', nso, nocc, nvir, 1.0, Cvp[0], nvir, &xp[offset], nvir, 0.0, Crp[0], nocc);

                offset += static_cast<long int>(nocc) * nvir;
            }

            C_right.push_back(Cr);
        }
    }

    jk_->compute();

    const std::vector<SharedMatrix>& J = jk_->J();
    const std::vector<SharedMatrix>& K = jk_->K();

    auto* Tp = new double[Caocc_->max_nrow() * Caocc_->max_ncol()];

    for (int symm = 0; symm < nirrep; ++symm) {
        for (size_t N = 0; N < x.size(); ++N) {
            double* bp = b[N]->pointer(symm);
            double* xp = x[N]->pointer(symm);
            long int offset = 0L;

            for (int h = 0; h < Caocc_->nirrep(); ++h) {
                int nsoocc = Caocc_->rowspi()[h];
                int nocc = Caocc_->colspi()[h];
                int nsovir = Cavir_->rowspi()[h ^ symm];
                int nvir = Cavir_->colspi()[h ^ symm];

                if (!nsoocc || !nsovir || !nocc || !nvir) continue;

                double** Cop = Caocc_->pointer(h);
                double** Cvp = Cavir_->pointer(h ^ symm);
                double* eop = eps_aocc_->pointer(h);
                double* evp = eps_avir_->pointer(h ^ symm);
                double** Jp = J[N]->pointer(h);
                double** Kp = K[N]->pointer(h);
                double** K2p = K[N]->pointer(h ^ symm);

                // 4(ia|jb)P_jb = C_im J_mn C_na
                C_DGEMM('T', 'N', nocc, nsovir, nsoocc, 1.0, Cop[0], nocc, Jp[0], nsovir, 0.0, Tp, nsovir);
                C_DGEMM('N', 'N', nocc, nvir, nsovir, 4.0, Tp, nsovir, Cvp[0], nvir, 0.0, &bp[offset], nvir);

                // -(ib|ja)P_jb = C_in K_nm C_ma
                C_DGEMM('T', 'T', nocc, nsovir, nsoocc, 1.0, Cop[0], nocc, K2p[0], nsoocc, 0.0, Tp, nsovir);
                C_DGEMM('N', 'N', nocc, nvir, nsovir, -1.0, Tp, nsovir, Cvp[0], nvir, 1.0, &bp[offset], nvir);

                // -(ij|ab)P_jb = C_im K_mn C_ra
                C_DGEMM('T', 'N', nocc, nsovir, nsoocc, 1.0, Cop[0], nocc, Kp[0], nsovir, 0.0, Tp, nsovir);
                C_DGEMM('N', 'N', nocc, nvir, nsovir, -1.0, Tp, nsovir, Cvp[0], nvir, 1.0, &bp[offset], nvir);

                for (int i = 0; i < nocc; ++i) {
                    for (int a = 0; a < nvir; ++a) {
                        bp[i * nvir + a + offset] += (evp[a] - eop[i]) * xp[i * nvir + a + offset];
                    }
                }

                offset += static_cast<long int>(nocc) * nvir;
            }
        }
    }

    delete[] Tp;

    if (debug_ > 3) {
        for (size_t N = 0; N < x.size(); N++) {
            x[N]->print();
            b[N]->print();
        }
    }
}
std::vector<SharedMatrix> CPHFRHamiltonian::unpack(const std::vector<std::shared_ptr<Vector> >& x) {
    std::vector<SharedMatrix> t1;
    int nirrep = x[0]->nirrep();
    for (size_t i = 0; i < x.size(); i++) {
        for (int symm = 0; symm < nirrep; ++symm) {
            auto t = std::make_shared<Matrix>("X", Caocc_->nirrep(), Caocc_->colspi(), Cavir_->colspi(), symm);
            long int offset = 0L;
            for (int h = 0; h < nirrep; ++h) {
                int nocc = Caocc_->colspi()[h];
                int nvir = Cavir_->colspi()[h ^ symm];

                if (!nocc || !nvir) continue;

                ::memcpy((void*)(t->pointer(h)[0]), (void*)(&x[i]->pointer(symm)[offset]),
                         sizeof(double) * nocc * nvir);

                offset += static_cast<long int>(nocc) * nvir;
            }

            t1.push_back(t);
        }
    }

    return t1;
}
}  // namespace psi
