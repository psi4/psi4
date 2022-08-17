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

UHamiltonian::UHamiltonian(std::shared_ptr<JK> jk) : Hamiltonian(jk) {}
UHamiltonian::UHamiltonian(std::shared_ptr<JK> jk, std::shared_ptr<VBase> v) : Hamiltonian(jk, v) {}
UHamiltonian::~UHamiltonian() {}

MatrixRHamiltonian::MatrixRHamiltonian(SharedMatrix M) : RHamiltonian(std::shared_ptr<JK>()), M_(M) {}
MatrixRHamiltonian::~MatrixRHamiltonian() {}
void MatrixRHamiltonian::print_header() const {
    if (print_) {
        outfile->Printf("  ==> MatrixRHamiltonian (by Rob Parrish) <== \n\n");
    }
}
std::shared_ptr<Vector> MatrixRHamiltonian::diagonal() {
    auto diag = std::make_shared<Vector>("Matrix Diagonal", M_->rowspi());
    for (int h = 0; h < M_->nirrep(); ++h) {
        int n = M_->rowspi()[h];
        if (!n) continue;
        double** Mp = M_->pointer(h);
        double* Dp = diag->pointer(h);
        for (int i = 0; i < n; ++i) {
            Dp[i] = Mp[i][i];
        }
    }
    return diag;
}
void MatrixRHamiltonian::product(const std::vector<std::shared_ptr<Vector> >& x,
                                 std::vector<std::shared_ptr<Vector> >& b) {
    for (size_t N = 0; N < x.size(); ++N) {
        for (int h = 0; h < M_->nirrep(); ++h) {
            int n = M_->rowspi()[h];
            if (!n) continue;
            double** Mp = M_->pointer(h);
            double* xp = x[N]->pointer(h);
            double* bp = b[N]->pointer(h);
            C_DGEMV('N', n, n, 1.0, Mp[0], n, xp, 1, 0.0, bp, 1);
        }
    }
}

MatrixUHamiltonian::MatrixUHamiltonian(std::pair<SharedMatrix, SharedMatrix> M)
    : UHamiltonian(std::shared_ptr<JK>()), M_(M) {}
MatrixUHamiltonian::~MatrixUHamiltonian() {}
void MatrixUHamiltonian::print_header() const {
    if (print_) {
        outfile->Printf("  ==> MatrixUHamiltonian (by Rob Parrish) <== \n\n");
    }
}
std::pair<std::shared_ptr<Vector>, std::shared_ptr<Vector> > MatrixUHamiltonian::diagonal() {
    auto diaga = std::make_shared<Vector>("Alpha Matrix Diagonal", M_.first->rowspi());
    auto diagb = std::make_shared<Vector>("Beta Matrix Diagonal", M_.first->rowspi());
    for (int h = 0; h < M_.first->nirrep(); ++h) {
        int n = M_.first->rowspi()[h];
        if (!n) continue;
        double** Map = M_.first->pointer(h);
        double* Dap = diaga->pointer(h);
        double** Mbp = M_.second->pointer(h);
        double* Dbp = diagb->pointer(h);
        for (int i = 0; i < n; ++i) {
            Dap[i] = Map[i][i];
            Dbp[i] = Mbp[i][i];
        }
    }
    return make_pair(diaga, diagb);
}
void MatrixUHamiltonian::product(const std::vector<std::pair<std::shared_ptr<Vector>, std::shared_ptr<Vector> > >& x,
                                 std::vector<std::pair<std::shared_ptr<Vector>, std::shared_ptr<Vector> > >& b) {
    for (size_t N = 0; N < x.size(); ++N) {
        for (int h = 0; h < M_.first->nirrep(); ++h) {
            int n = M_.first->rowspi()[h];
            if (!n) continue;
            double** Map = M_.first->pointer(h);
            double* xap = x[N].first->pointer(h);
            double* bap = b[N].first->pointer(h);
            C_DGEMV('N', n, n, 1.0, Map[0], n, xap, 1, 0.0, bap, 1);
            //            double** Mbp = M_.second->pointer(h);
            //            double*  xbp = x[N].second->pointer(h);
            //            double*  bbp = b[N].second->pointer(h);
            C_DGEMV('N', n, n, 1.0, Map[0], n, xap, 1, 0.0, bap, 1);
        }
    }
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

// Hamiltonian for UHF stability analysis
USTABHamiltonian::USTABHamiltonian(std::shared_ptr<JK> jk, SharedMatrix Cocca, SharedMatrix Cvira, SharedMatrix Coccb,
                                   SharedMatrix Cvirb, std::shared_ptr<Vector> eps_occa,
                                   std::shared_ptr<Vector> eps_vira, std::shared_ptr<Vector> eps_occb,
                                   std::shared_ptr<Vector> eps_virb, std::shared_ptr<VBase> v)
    : UHamiltonian(jk, v),
      Cocca_(Cocca),
      Cvira_(Cvira),
      Coccb_(Coccb),
      Cvirb_(Cvirb),
      eps_occa_(eps_occa),
      eps_vira_(eps_vira),
      eps_occb_(eps_occb),
      eps_virb_(eps_virb) {}
USTABHamiltonian::~USTABHamiltonian() {}
void USTABHamiltonian::print_header() const {
    if (print_) {
        outfile->Printf("  ==> USTABHamiltonian (by Jérôme Gonthier) <== \n");
        outfile->Printf("  ==> Inspired by R.Parrish CISRHamiltonian <== \n\n");
    }
}
std::pair<std::shared_ptr<Vector>, std::shared_ptr<Vector> > USTABHamiltonian::diagonal() {
    int nirrepa = eps_occa_->nirrep();
    Dimension nova(nirrepa);
    for (int symm = 0; symm < nirrepa; ++symm) {
        for (int h = 0; h < nirrepa; ++h) {
            nova[symm] += eps_occa_->dimpi()[h] * eps_vira_->dimpi()[symm ^ h];
        }
    }
    int nirrepb = eps_occb_->nirrep();
    Dimension novb(nirrepb);
    for (int symm = 0; symm < nirrepb; ++symm) {
        for (int h = 0; h < nirrepb; ++h) {
            novb[symm] += eps_occb_->dimpi()[h] * eps_virb_->dimpi()[symm ^ h];
        }
    }

    auto diaga = std::make_shared<Vector>("UStab Alpha Diagonal", nova);
    auto diagb = std::make_shared<Vector>("UStab Beta Diagonal", novb);

    for (int symm = 0; symm < nirrepa; ++symm) {
        long int offset = 0L;
        for (int h = 0; h < nirrepa; ++h) {
            int nocc = eps_occa_->dimpi()[h];
            int nvir = eps_vira_->dimpi()[symm ^ h];

            if (!nocc || !nvir) continue;

            double* eop = eps_occa_->pointer(h);
            double* evp = eps_vira_->pointer(symm ^ h);
            double* dp = diaga->pointer(symm);

            for (int i = 0; i < nocc; ++i) {
                for (int a = 0; a < nvir; ++a) {
                    dp[i * nvir + a + offset] = evp[a] - eop[i];
                }
            }
            offset += static_cast<long int>(nocc) * nvir;
        }
    }
    for (int symm = 0; symm < nirrepb; ++symm) {
        long int offset = 0L;
        for (int h = 0; h < nirrepb; ++h) {
            int nocc = eps_occb_->dimpi()[h];
            int nvir = eps_virb_->dimpi()[symm ^ h];

            if (!nocc || !nvir) continue;

            double* eop = eps_occb_->pointer(h);
            double* evp = eps_virb_->pointer(symm ^ h);
            double* dp = diagb->pointer(symm);

            for (int i = 0; i < nocc; ++i) {
                for (int a = 0; a < nvir; ++a) {
                    dp[i * nvir + a + offset] = evp[a] - eop[i];
                }
            }
            offset += static_cast<long int>(nocc) * nvir;
        }
    }

    if (exact_diagonal_) {
        outfile->Printf("No exact diagonal available for UStab Hamiltionan.\n\n");
        outfile->Printf("Providing orbital energy difference");
    }

    return make_pair(diaga, diagb);
}

void USTABHamiltonian::product(const std::vector<std::pair<std::shared_ptr<Vector>, std::shared_ptr<Vector> > >& x,
                               std::vector<std::pair<std::shared_ptr<Vector>, std::shared_ptr<Vector> > >& b) {
    std::vector<SharedMatrix>& C_left = jk_->C_left();
    std::vector<SharedMatrix>& C_right = jk_->C_right();

    C_left.clear();
    C_right.clear();

    int nirrepa = (x.size() ? x[0].first->nirrep() : 0);
    int nirrepb = (x.size() ? x[0].second->nirrep() : 0);

    //    Alpha orbitals are handled first

    for (int symm = 0; symm < nirrepa; ++symm) {
        for (size_t N = 0; N < x.size(); ++N) {
            double* xp = x[N].first->pointer(symm);
            long int offset = 0L;

            // Is this a delta (low-rank) vector?
            bool delta_vector = true;
            int delta_h = 0;
            int delta_i = 0;
            int delta_a = 0;
            double delta_sum = 0.0;
            for (int h = 0; h < Cocca_->nirrep(); ++h) {
                int nocc = Cocca_->colspi()[h];
                int nvir = Cvira_->colspi()[h ^ symm];
                for (int i = 0; i < nocc; i++) {
                    for (int a = 0; a < nvir; a++) {
                        double our_x = xp[i * nvir + a + offset];
                        if (our_x != 0.0 && our_x != 1.0) {
                            delta_vector = false;
                            break;
                        }
                        delta_sum += our_x;
                        if (our_x == 1.0) {
                            delta_i = i;
                            delta_a = a;
                            delta_h = h;
                        }
                    }
                }
                offset += static_cast<long int>(nocc) * nvir;
            }
            if (delta_sum != 1.0) {
                delta_vector = false;
            }

            if (delta_vector) {
                Dimension rank(Cocca_->nirrep());
                rank[delta_h] = 1;

                std::stringstream ss;
                ss << "C_left alpha, h = " << symm << ", N = " << N;
                auto Cl = std::make_shared<Matrix>(ss.str(), Cocca_->nirrep(), Cocca_->rowspi(), rank);

                double** Clp = Cl->pointer(delta_h);
                double** Cop = Cocca_->pointer(delta_h);
                C_DCOPY(Cocca_->rowspi()[delta_h], &Cop[0][delta_i], Cocca_->colspi()[delta_h], Clp[0], 1);

                std::stringstream ss2;
                ss2 << "C_right alpha, h = " << symm << ", N = " << N;
                auto Cr = std::make_shared<Matrix>(ss2.str(), Cocca_->nirrep(), Cocca_->rowspi(), rank, symm);

                double** Crp = Cr->pointer(delta_h ^ symm);
                double** Cvp = Cvira_->pointer(delta_h ^ symm);
                C_DCOPY(Cvira_->rowspi()[delta_h ^ symm], &Cvp[0][delta_a], Cvira_->colspi()[delta_h ^ symm], Crp[0],
                        1);

                C_left.push_back(Cl);
                C_right.push_back(Cr);

            } else {
                std::stringstream ss;
                ss << "C_right alpha, h = " << symm << ", N = " << N;
                auto Cr =
                    std::make_shared<Matrix>(ss.str(), Cocca_->nirrep(), Cocca_->rowspi(), Cocca_->colspi(), symm);

                offset = 0L;
                for (int h = 0; h < Cocca_->nirrep(); ++h) {
                    int nocc = Cocca_->colspi()[h];
                    int nvir = Cvira_->colspi()[h ^ symm];
                    int nso = Cvira_->rowspi()[h ^ symm];

                    if (!nso || !nocc || !nvir) continue;

                    double** Cvp = Cvira_->pointer(h ^ symm);
                    double** Crp = Cr->pointer(h ^ symm);

                    C_DGEMM('N', 'T', nso, nocc, nvir, 1.0, Cvp[0], nvir, &xp[offset], nvir, 0.0, Crp[0], nocc);

                    offset += static_cast<long int>(nocc) * nvir;
                }

                C_left.push_back(Cocca_);
                C_right.push_back(Cr);
            }
        }
    }

    //    Then beta orbitals get exactly the same treatment. Rewrite that in a smarter way.

    for (int symm = 0; symm < nirrepb; ++symm) {
        for (size_t N = 0; N < x.size(); ++N) {
            double* xp = x[N].second->pointer(symm);
            long int offset = 0L;

            // Is this a delta (low-rank) vector?
            bool delta_vector = true;
            int delta_h = 0;
            int delta_i = 0;
            int delta_a = 0;
            double delta_sum = 0.0;
            for (int h = 0; h < Coccb_->nirrep(); ++h) {
                int nocc = Coccb_->colspi()[h];
                int nvir = Cvirb_->colspi()[h ^ symm];
                for (int i = 0; i < nocc; i++) {
                    for (int a = 0; a < nvir; a++) {
                        double our_x = xp[i * nvir + a + offset];
                        if (our_x != 0.0 && our_x != 1.0) {
                            delta_vector = false;
                            break;
                        }
                        delta_sum += our_x;
                        if (our_x == 1.0) {
                            delta_i = i;
                            delta_a = a;
                            delta_h = h;
                        }
                    }
                }
                offset += static_cast<long int>(nocc) * nvir;
            }
            if (delta_sum != 1.0) {
                delta_vector = false;
            }

            if (delta_vector) {
                Dimension rank(Coccb_->nirrep());
                rank[delta_h] = 1;

                std::stringstream ss;
                ss << "C_left beta, h = " << symm << ", N = " << N;
                auto Cl = std::make_shared<Matrix>(ss.str(), Coccb_->nirrep(), Coccb_->rowspi(), rank);

                double** Clp = Cl->pointer(delta_h);
                double** Cop = Coccb_->pointer(delta_h);
                C_DCOPY(Coccb_->rowspi()[delta_h], &Cop[0][delta_i], Coccb_->colspi()[delta_h], Clp[0], 1);

                std::stringstream ss2;
                ss2 << "C_right beta, h = " << symm << ", N = " << N;
                auto Cr = std::make_shared<Matrix>(ss2.str(), Coccb_->nirrep(), Coccb_->rowspi(), rank, symm);

                double** Crp = Cr->pointer(delta_h ^ symm);
                double** Cvp = Cvirb_->pointer(delta_h ^ symm);
                C_DCOPY(Cvirb_->rowspi()[delta_h ^ symm], &Cvp[0][delta_a], Cvirb_->colspi()[delta_h ^ symm], Crp[0],
                        1);

                C_left.push_back(Cl);
                C_right.push_back(Cr);

            } else {
                std::stringstream ss;
                ss << "C_right beta, h = " << symm << ", N = " << N;
                auto Cr =
                    std::make_shared<Matrix>(ss.str(), Coccb_->nirrep(), Coccb_->rowspi(), Coccb_->colspi(), symm);

                offset = 0L;
                for (int h = 0; h < Coccb_->nirrep(); ++h) {
                    int nocc = Coccb_->colspi()[h];
                    int nvir = Cvirb_->colspi()[h ^ symm];
                    int nso = Cvirb_->rowspi()[h ^ symm];

                    if (!nso || !nocc || !nvir) continue;

                    double** Cvp = Cvirb_->pointer(h ^ symm);
                    double** Crp = Cr->pointer(h ^ symm);

                    C_DGEMM('N', 'T', nso, nocc, nvir, 1.0, Cvp[0], nvir, &xp[offset], nvir, 0.0, Crp[0], nocc);

                    offset += static_cast<long int>(nocc) * nvir;
                }

                C_left.push_back(Coccb_);
                C_right.push_back(Cr);
            }
        }
    }

    jk_->compute();

    const std::vector<SharedMatrix>& J = jk_->J();
    const std::vector<SharedMatrix>& K = jk_->K();

    //    Compute the alpha part of the b vector

    auto* Tp = new double[Cocca_->max_nrow() * Cocca_->max_ncol()];

    for (int symm = 0; symm < nirrepa; ++symm) {
        for (size_t N = 0; N < x.size(); ++N) {
            double* bp = b[N].first->pointer(symm);
            double* xp = x[N].first->pointer(symm);
            long int offset = 0L;

            for (int h = 0; h < Cocca_->nirrep(); ++h) {
                int nsoocc = Cocca_->rowspi()[h];
                int nocc = Cocca_->colspi()[h];
                int nsovir = Cvira_->rowspi()[h ^ symm];
                int nvir = Cvira_->colspi()[h ^ symm];

                if (!nsoocc || !nsovir || !nocc || !nvir) continue;

                double** Cop = Cocca_->pointer(h);
                double** Cvp = Cvira_->pointer(h ^ symm);
                double* eop = eps_occa_->pointer(h);
                double* evp = eps_vira_->pointer(h ^ symm);

                double** Jap = J[symm * x.size() + N]->pointer(h);
                double** Jbp = J[(nirrepa + symm) * x.size() + N]->pointer(h);
                double** Kp = K[symm * x.size() + N]->pointer(h);
                double** KTp = K[symm * x.size() + N]->pointer(h ^ symm);
                // We need to use h^symm representation of h because we want the
                // columns to be in h^symm so that the transpose has lines in h^symm

                // -(ij|ab)P_jb = C_im K_mn C_na
                C_DGEMM('T', 'N', nocc, nsovir, nsoocc, 1.0, Cop[0], nocc, Kp[0], nsovir, 0.0, Tp, nsovir);
                C_DGEMM('N', 'N', nocc, nvir, nsovir, -1.0, Tp, nsovir, Cvp[0], nvir, 0.0, &bp[offset], nvir);

                // -(ib|ja) P_jb = C_im (K^{T})_mn C_na
                C_DGEMM('T', 'T', nocc, nsovir, nsoocc, 1.0, Cop[0], nocc, KTp[0], nsoocc, 0.0, Tp, nsovir);
                C_DGEMM('N', 'N', nocc, nvir, nsovir, -1.0, Tp, nsovir, Cvp[0], nvir, 1.0, &bp[offset], nvir);

                // 2(ia|jb)P_jb = C_im J_mn C_na for J alpha
                C_DGEMM('T', 'N', nocc, nsovir, nsoocc, 1.0, Cop[0], nocc, Jap[0], nsovir, 0.0, Tp, nsovir);
                C_DGEMM('N', 'N', nocc, nvir, nsovir, 2.0, Tp, nsovir, Cvp[0], nvir, 1.0, &bp[offset], nvir);

                // 2(ia|jb)P_jb = C_im J_mn C_na for J beta
                C_DGEMM('T', 'N', nocc, nsovir, nsoocc, 1.0, Cop[0], nocc, Jbp[0], nsovir, 0.0, Tp, nsovir);
                C_DGEMM('N', 'N', nocc, nvir, nsovir, 2.0, Tp, nsovir, Cvp[0], nvir, 1.0, &bp[offset], nvir);

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

    //    Compute the beta part of the b vector

    Tp = new double[Coccb_->max_nrow() * Coccb_->max_ncol()];

    for (int symm = 0; symm < nirrepb; ++symm) {
        for (size_t N = 0; N < x.size(); ++N) {
            double* bp = b[N].second->pointer(symm);
            double* xp = x[N].second->pointer(symm);
            long int offset = 0L;

            for (int h = 0; h < Coccb_->nirrep(); ++h) {
                int nsoocc = Coccb_->rowspi()[h];
                int nocc = Coccb_->colspi()[h];
                int nsovir = Cvirb_->rowspi()[h ^ symm];
                int nvir = Cvirb_->colspi()[h ^ symm];

                if (!nsoocc || !nsovir || !nocc || !nvir) continue;

                double** Cop = Coccb_->pointer(h);
                double** Cvp = Cvirb_->pointer(h ^ symm);
                double* eop = eps_occb_->pointer(h);
                double* evp = eps_virb_->pointer(h ^ symm);

                double** Jap = J[symm * x.size() + N]->pointer(h);
                double** Jbp = J[(nirrepa + symm) * x.size() + N]->pointer(h);
                double** Kp = K[(nirrepa + symm) * x.size() + N]->pointer(h);
                double** KTp = K[(nirrepa + symm) * x.size() + N]->pointer(h ^ symm);

                // -(ij|ab)P_jb = C_im K_mn C_na
                C_DGEMM('T', 'N', nocc, nsovir, nsoocc, 1.0, Cop[0], nocc, Kp[0], nsovir, 0.0, Tp, nsovir);
                C_DGEMM('N', 'N', nocc, nvir, nsovir, -1.0, Tp, nsovir, Cvp[0], nvir, 0.0, &bp[offset], nvir);

                // -(ib|ja) P_jb = C_im (K^{T})_mn C_na
                C_DGEMM('T', 'T', nocc, nsovir, nsoocc, 1.0, Cop[0], nocc, KTp[0], nsoocc, 0.0, Tp, nsovir);
                C_DGEMM('N', 'N', nocc, nvir, nsovir, -1.0, Tp, nsovir, Cvp[0], nvir, 1.0, &bp[offset], nvir);

                // 2(ia|jb)P_jb = C_im J_mn C_na for J alpha
                C_DGEMM('T', 'N', nocc, nsovir, nsoocc, 1.0, Cop[0], nocc, Jap[0], nsovir, 0.0, Tp, nsovir);
                C_DGEMM('N', 'N', nocc, nvir, nsovir, 2.0, Tp, nsovir, Cvp[0], nvir, 1.0, &bp[offset], nvir);

                // 2(ia|jb)P_jb = C_im J_mn C_na for J beta
                C_DGEMM('T', 'N', nocc, nsovir, nsoocc, 1.0, Cop[0], nocc, Jbp[0], nsovir, 0.0, Tp, nsovir);
                C_DGEMM('N', 'N', nocc, nvir, nsovir, 2.0, Tp, nsovir, Cvp[0], nvir, 1.0, &bp[offset], nvir);

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
            x[N].first->print();
            b[N].first->print();
        }
        for (size_t N = 0; N < x.size(); N++) {
            x[N].second->print();
            b[N].second->print();
        }
    }
}

// No implementation here, working with pairs is annoying.
// std::vector<std::pair<SharedMatrix, SharedMatrix> > USTABHamiltonian::unpack(
//        const std::pair<std::shared_ptr<Vector>, std::shared_ptr<Vector> >& eig)
//{
//}
// New, better function for our purpose.
std::vector<std::pair<SharedMatrix, SharedMatrix> > USTABHamiltonian::unpack_paired(
    const std::shared_ptr<Vector>& eig) {
    int nirrep = eig->nirrep();
    std::vector<std::pair<SharedMatrix, SharedMatrix> > t1;
    for (int symm = 0; symm < nirrep; ++symm) {
        auto talpha = std::make_shared<Matrix>("T", Cocca_->nirrep(), Cocca_->colspi(), Cvira_->colspi(), symm);
        auto tbeta = std::make_shared<Matrix>("T", Coccb_->nirrep(), Coccb_->colspi(), Cvirb_->colspi(), symm);
        long int offseta = 0L;
        long int offsetb = 0L;
        for (int h = 0; h < nirrep; ++h) {
            int nocca = Cocca_->colspi()[h];
            int nvira = Cvira_->colspi()[h ^ symm];

            if (!nocca || !nvira) continue;

            ::memcpy((void*)(talpha->pointer(h)[0]), (void*)(&eig->pointer(symm)[offseta]),
                     sizeof(double) * nocca * nvira);

            offseta += static_cast<long int>(nocca) * nvira;
        }
        for (int h = 0; h < nirrep; ++h) {
            int noccb = Coccb_->colspi()[h];
            int nvirb = Cvirb_->colspi()[h ^ symm];

            if (!noccb || !nvirb) continue;

            ::memcpy((void*)(tbeta->pointer(h)[0]), (void*)(&eig->pointer(symm)[offseta + offsetb]),
                     sizeof(double) * noccb * nvirb);

            offsetb += static_cast<long int>(noccb) * nvirb;
        }

        t1.push_back(make_pair(talpha, tbeta));
    }

    return t1;
}
}  // namespace psi
