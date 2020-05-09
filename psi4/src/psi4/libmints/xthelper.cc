/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2019 The Psi4 Developers.
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

#include "xthelper.h"

#include <memory>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "psi4/psifiles.h"

#include "basisset.h"
#include "dimension.h"
#include "integral.h"
#include "mintshelper.h"
#include "tensor.h"

namespace psi {
template <typename Mints, typename T, typename Enable>
auto xtHelper<Mints, T, Enable>::ao_overlap_() const -> SharedMatrix_<T> {
    auto nthread_ = static_cast<const Mints *>(this)->nthread_;
    std::vector<Shared<OneBodyAOInt>> ints_vec;
    for (auto i = 0; i < nthread_; i++) {
        ints_vec.push_back(Shared<OneBodyAOInt>(static_cast<const Mints *>(this)->integral_->ao_overlap()));
    }
    auto overlap_mat = std::make_shared<Matrix_<T>>(PSIF_AO_S, static_cast<const Mints *>(this)->basisset_->nbf(),
                                                    static_cast<const Mints *>(this)->basisset_->nbf());
    one_body_ao_computer(ints_vec, overlap_mat, true);
    return overlap_mat;
}

template <typename Mints, typename T, typename Enable>
auto xtHelper<Mints, T, Enable>::ao_overlap_(Shared<BasisSet> bs1, Shared<BasisSet> bs2) const -> SharedMatrix_<T> {
    auto nthread_ = static_cast<const Mints *>(this)->nthread_;
    IntegralFactory factory(bs1, bs2, bs1, bs2);
    std::vector<Shared<OneBodyAOInt>> ints_vec;
    for (auto i = 0; i < nthread_; i++) {
        ints_vec.push_back(Shared<OneBodyAOInt>(factory.ao_overlap()));
    }
    auto overlap_mat = std::make_shared<Matrix_<T>>(PSIF_AO_S, bs1->nbf(), bs2->nbf());
    one_body_ao_computer(ints_vec, overlap_mat, false);
    return overlap_mat;
}

template <typename Mints, typename T, typename Enable>
auto xtHelper<Mints, T, Enable>::ao_kinetic_() const -> SharedMatrix_<T> {
    auto nthread_ = static_cast<const Mints *>(this)->nthread_;
    std::vector<Shared<OneBodyAOInt>> ints_vec;
    for (auto i = 0; i < nthread_; i++) {
        ints_vec.push_back(Shared<OneBodyAOInt>(static_cast<const Mints *>(this)->integral_->ao_kinetic()));
    }
    auto kinetic_mat =
        std::make_shared<Matrix_<T>>("AO-basis Kinetic Ints", static_cast<const Mints *>(this)->basisset_->nbf(),
                                     static_cast<const Mints *>(this)->basisset_->nbf());
    one_body_ao_computer(ints_vec, kinetic_mat, true);
    return kinetic_mat;
}

template <typename Mints, typename T, typename Enable>
auto xtHelper<Mints, T, Enable>::ao_kinetic_(Shared<BasisSet> bs1, Shared<BasisSet> bs2) const -> SharedMatrix_<T> {
    auto nthread_ = static_cast<const Mints *>(this)->nthread_;
    IntegralFactory factory(bs1, bs2, bs1, bs2);
    std::vector<Shared<OneBodyAOInt>> ints_vec;
    for (auto i = 0; i < nthread_; i++) {
        ints_vec.push_back(Shared<OneBodyAOInt>(factory.ao_kinetic()));
    }
    auto kinetic_mat = std::make_shared<Matrix_<T>>("AO-basis Kinetic Ints", bs1->nbf(), bs2->nbf());
    one_body_ao_computer(ints_vec, kinetic_mat, false);
    return kinetic_mat;
}

template <typename Mints, typename T, typename Enable>
auto xtHelper<Mints, T, Enable>::ao_potential_() const -> SharedMatrix_<T> {
    auto nthread_ = static_cast<const Mints *>(this)->nthread_;
    std::vector<Shared<OneBodyAOInt>> ints_vec;
    for (auto i = 0; i < nthread_; i++) {
        ints_vec.push_back(Shared<OneBodyAOInt>(static_cast<const Mints *>(this)->integral_->ao_potential()));
    }
    auto potential_mat =
        std::make_shared<Matrix_<T>>("AO-basis Potential Ints", static_cast<const Mints *>(this)->basisset_->nbf(),
                                     static_cast<const Mints *>(this)->basisset_->nbf());
    one_body_ao_computer(ints_vec, potential_mat, true);
    return potential_mat;
}

template <typename Mints, typename T, typename Enable>
auto xtHelper<Mints, T, Enable>::ao_potential_(Shared<BasisSet> bs1, Shared<BasisSet> bs2) const -> SharedMatrix_<T> {
    auto nthread_ = static_cast<const Mints *>(this)->nthread_;
    IntegralFactory factory(bs1, bs2, bs1, bs2);
    std::vector<Shared<OneBodyAOInt>> ints_vec;
    for (auto i = 0; i < nthread_; i++) {
        ints_vec.push_back(Shared<OneBodyAOInt>(factory.ao_potential()));
    }
    auto potential_mat = std::make_shared<Matrix_<T>>("AO-basis Potential Ints", bs1->nbf(), bs2->nbf());
    one_body_ao_computer(ints_vec, potential_mat, false);
    return potential_mat;
}

template <typename Mints, typename T, typename Enable>
auto xtHelper<Mints, T, Enable>::ao_eri_(Shared<IntegralFactory> input_factory) const -> SharedTensor<T, 4> {
    auto factory = input_factory ? input_factory : static_cast<const Mints *>(this)->integral_;
    return ao_helper("AO ERI Tensor", Shared<TwoBodyAOInt>(factory->eri()));
}

template <typename Mints, typename T, typename Enable>
auto xtHelper<Mints, T, Enable>::ao_eri_(Shared<BasisSet> bs1, Shared<BasisSet> bs2, Shared<BasisSet> bs3,
                                         Shared<BasisSet> bs4) const -> SharedTensor<T, 4> {
    IntegralFactory intf(bs1, bs2, bs3, bs4);
    Shared<TwoBodyAOInt> ints(intf.eri());
    return ao_helper("AO ERI Tensor", ints);
}

template <typename Mints, typename T, typename Enable>
auto xtHelper<Mints, T, Enable>::ao_eri_shell_(int M, int N, int P, int Q) -> SharedTensor<T, 4> {
    if (static_cast<const Mints *>(this)->eriInts_ == nullptr) {
        static_cast<Mints *>(this)->eriInts_ = Shared<TwoBodyAOInt>(static_cast<const Mints *>(this)->integral_->eri());
    }
    return ao_shell_getter("AO ERI Tensor", static_cast<const Mints *>(this)->eriInts_, M, N, P, Q);
}

template <typename Mints, typename T, typename Enable>
auto xtHelper<Mints, T, Enable>::one_body_ao_computer(std::vector<Shared<OneBodyAOInt>> ints, SharedMatrix_<T> out,
                                                      bool symm) const -> void {
    // Grab basis info
    Shared<BasisSet> bs1 = ints[0]->basis1();
    Shared<BasisSet> bs2 = ints[0]->basis2();

    // Limit to the number of incoming onbody ints
    auto nthread = static_cast<const Mints *>(this)->nthread_;
    if (nthread > ints.size()) {
        nthread = ints.size();
    }

    // Grab the buffers
    std::vector<const double *> ints_buff(nthread);
    for (auto thread = 0; thread < nthread; thread++) {
        ints_buff[thread] = ints[thread]->buffer();
    }

// Loop it
#pragma omp parallel for schedule(guided) num_threads(nthread)
    for (auto MU = 0; MU < bs1->nshell(); ++MU) {
        const auto num_mu = bs1->shell(MU).nfunction();
        const auto index_mu = bs1->shell(MU).function_index();

        size_t rank = 0;
#ifdef _OPENMP
        rank = omp_get_thread_num();
#endif

        if (symm) {
            // Triangular
            for (auto NU = 0; NU <= MU; ++NU) {
                const auto num_nu = bs2->shell(NU).nfunction();
                const auto index_nu = bs2->shell(NU).function_index();

                ints[rank]->compute_shell(MU, NU);

                size_t index = 0;
                for (auto mu = index_mu; mu < (index_mu + num_mu); ++mu) {
                    for (auto nu = index_nu; nu < (index_nu + num_nu); ++nu) {
                        out->set_symmetric(nu, mu, ints_buff[rank][index++]);
                    }
                }
            }  // End NU
        }      // End Symm
        else {
            // Rectangular
            for (auto NU = 0; NU < bs2->nshell(); ++NU) {
                const auto num_nu = bs2->shell(NU).nfunction();
                const auto index_nu = bs2->shell(NU).function_index();

                ints[rank]->compute_shell(MU, NU);

                size_t index = 0;
                for (auto mu = index_mu; mu < (index_mu + num_mu); ++mu) {
                    for (auto nu = index_nu; nu < (index_nu + num_nu); ++nu) {
                        // printf("%zu %zu | %zu %zu | %lf\n", MU, NU, mu, nu, ints_buff[rank][index]);
                        out->set(mu, nu, ints_buff[rank][index++]);
                    }
                }
            }  // End NU
        }      // End Rectangular
    }          // End Mu
}

template <typename Mints, typename T, typename Enable>
auto xtHelper<Mints, T, Enable>::ao_helper(const std::string &label, Shared<TwoBodyAOInt> ints) const
    -> SharedTensor<T, 4> {
    Shared<BasisSet> bs1 = ints->basis1();
    Shared<BasisSet> bs2 = ints->basis2();
    Shared<BasisSet> bs3 = ints->basis3();
    Shared<BasisSet> bs4 = ints->basis4();

    int nbf1 = bs1->nbf();
    int nbf2 = bs2->nbf();
    int nbf3 = bs3->nbf();
    int nbf4 = bs4->nbf();

    auto dims = std::array<Dimension::value_type, 4>{nbf1, nbf2, nbf3, nbf4};
    auto I = std::make_shared<Tensor<T, 4>>(label, dims);
    const double *buffer = ints->buffer();

    for (auto M = 0; M < bs1->nshell(); M++) {
        for (auto N = 0; N < bs2->nshell(); N++) {
            for (auto P = 0; P < bs3->nshell(); P++) {
                for (auto Q = 0; Q < bs4->nshell(); Q++) {
                    ints->compute_shell(M, N, P, Q);

                    for (auto m = 0, index = 0; m < bs1->shell(M).nfunction(); m++) {
                        for (auto n = 0; n < bs2->shell(N).nfunction(); n++) {
                            for (auto p = 0; p < bs3->shell(P).nfunction(); p++) {
                                for (auto q = 0; q < bs4->shell(Q).nfunction(); q++, index++) {
                                    I->set(bs1->shell(M).function_index() + m, bs2->shell(N).function_index() + n,
                                           bs3->shell(P).function_index() + p, bs4->shell(Q).function_index() + q,
                                           buffer[index]);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    return I;
}

template <typename Mints, typename T, typename Enable>
auto xtHelper<Mints, T, Enable>::ao_shell_getter(const std::string &label, Shared<TwoBodyAOInt> ints, int M, int N,
                                                 int P, int Q) const -> SharedTensor<T, 4> {
    auto mfxn = static_cast<const Mints *>(this)->basisset_->shell(M).nfunction();
    auto nfxn = static_cast<const Mints *>(this)->basisset_->shell(N).nfunction();
    auto pfxn = static_cast<const Mints *>(this)->basisset_->shell(P).nfunction();
    auto qfxn = static_cast<const Mints *>(this)->basisset_->shell(Q).nfunction();

    auto dims = std::array<Dimension::value_type, 4>{mfxn, nfxn, pfxn, qfxn};
    auto I = std::make_shared<Tensor<T, 4>>(label, dims);
    const double *buffer = ints->buffer();

    ints->compute_shell(M, N, P, Q);

    for (auto m = 0, index = 0; m < mfxn; m++) {
        for (auto n = 0; n < nfxn; n++) {
            for (auto p = 0; p < pfxn; p++) {
                for (auto q = 0; q < qfxn; q++, index++) {
                    I->set(m, n, p, q, buffer[index]);
                }
            }
        }
    }

    return I;
}

template auto xtHelper<MintsHelper, double>::ao_overlap_() const -> SharedMatrix_<double>;
template auto xtHelper<MintsHelper, double>::ao_overlap_(Shared<BasisSet>, Shared<BasisSet>) const
    -> SharedMatrix_<double>;

template auto xtHelper<MintsHelper, double>::ao_kinetic_() const -> SharedMatrix_<double>;
template auto xtHelper<MintsHelper, double>::ao_kinetic_(Shared<BasisSet>, Shared<BasisSet>) const
    -> SharedMatrix_<double>;

template auto xtHelper<MintsHelper, double>::ao_potential_() const -> SharedMatrix_<double>;
template auto xtHelper<MintsHelper, double>::ao_potential_(Shared<BasisSet>, Shared<BasisSet>) const
    -> SharedMatrix_<double>;

template auto xtHelper<MintsHelper, double>::ao_eri_(Shared<IntegralFactory> = nullptr) const
    -> SharedTensor<double, 4>;
template auto xtHelper<MintsHelper, double>::ao_eri_(Shared<BasisSet> bs1, Shared<BasisSet> bs2, Shared<BasisSet> bs3,
                                                     Shared<BasisSet> bs4) const -> SharedTensor<double, 4>;

template auto xtHelper<MintsHelper, double>::ao_eri_shell_(int M, int N, int P, int Q) -> SharedTensor<double, 4>;

// template auto xtHelper<MintsHelper, double>::so_overlap_(bool include_perturbations = true) -> SharedMatrix_<double>;
// template auto xtHelper<MintsHelper, double>::so_kinetic_(bool include_perturbations = true) -> SharedMatrix_<double>;
// template auto xtHelper<MintsHelper, double>::so_potential_(bool include_perturbations = true) ->
// SharedMatrix_<double>;

template auto xtHelper<MintsHelper, double>::one_body_ao_computer(std::vector<Shared<OneBodyAOInt>> ints,
                                                                  SharedMatrix_<double> out, bool symm) const -> void;
template auto xtHelper<MintsHelper, double>::ao_helper(const std::string &label, Shared<TwoBodyAOInt> ints) const
    -> SharedTensor<double, 4>;
template auto xtHelper<MintsHelper, double>::ao_shell_getter(const std::string &label, Shared<TwoBodyAOInt> ints, int M,
                                                             int N, int P, int Q) const -> SharedTensor<double, 4>;
}  // namespace psi
