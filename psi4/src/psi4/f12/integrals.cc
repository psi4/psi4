/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2024 The Psi4 Developers.
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

#include "mp2.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#include "psi4/libmints/basisset.h"
#include "psi4/libmints/integral.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/mintshelper.h"
#include "psi4/libmints/onebody.h"
#include "psi4/lib3index/dftensor.h"

#include "einsums.hpp"

namespace psi { namespace f12 {

void MP2F12::convert_C(einsums::Tensor<double,2> *C, OrbitalSpace bs, const int& dim1, const int& dim2)
{
    for (int p = 0; p < dim1; p++) {
        for (int q = 0; q < dim2; q++) {
            (*C)(p, q) = bs.C()->get(p, q);
        }
    }
}

void MP2F12::set_ERI(einsums::TensorView<double, 4>& ERI_Slice, einsums::Tensor<double, 4> *Slice)
{
    const auto dim1 = (*Slice).dim(0);
    const auto dim2 = (*Slice).dim(1);
    const auto dim3 = (*Slice).dim(2);
    const auto dim4 = (*Slice).dim(3);

    for (int p = 0; p < dim1; p++){
        for (int q = 0; q < dim2; q++){
            for (int r = 0; r < dim3; r++){
                for (int s = 0; s < dim4; s++){
                    ERI_Slice(p, q, r, s) = (*Slice)(p, q, r, s);
                }
            }
        }
    }
}

void MP2F12::set_ERI(einsums::TensorView<double, 3>& ERI_Slice, einsums::Tensor<double, 3> *Slice)
{
    const auto naux_= (*Slice).dim(0);
    const auto dim1 = (*Slice).dim(1);
    const auto dim2 = (*Slice).dim(2);

    for (int A = 0; A < naux_; A++){
        for (int p = 0; p < dim1; p++){
            for (int q = 0; q < dim2; q++){
                ERI_Slice(A, p, q) = (*Slice)(A, p, q);
            }
        }
    }
}

void MP2F12::two_body_ao_computer(const std::string& int_type, einsums::Tensor<double, 4> *GAO,
                                  std::shared_ptr<BasisSet> bs1, std::shared_ptr<BasisSet> bs2,
                                  std::shared_ptr<BasisSet> bs3, std::shared_ptr<BasisSet> bs4)
{
    std::shared_ptr<IntegralFactory> intf(new IntegralFactory(bs1, bs2, bs3, bs4));

    std::vector<std::shared_ptr<TwoBodyAOInt>> ints;
    if ( int_type == "F" ){
        ints.push_back(std::shared_ptr<TwoBodyAOInt>(intf->f12(cgtg_)));
    } else if ( int_type == "FG" ){
        ints.push_back(std::shared_ptr<TwoBodyAOInt>(intf->f12g12(cgtg_)));
    } else if ( int_type == "F2" ){
        ints.push_back(std::shared_ptr<TwoBodyAOInt>(intf->f12_squared(cgtg_)));
    } else if ( int_type == "Uf" ){
        ints.push_back(std::shared_ptr<TwoBodyAOInt>(intf->f12_double_commutator(cgtg_)));
    } else {
        ints.push_back(std::shared_ptr<TwoBodyAOInt>(intf->eri()));
    }

    // Make ints vector
    for (size_t thread = 1; thread < nthreads_; thread++) {
        ints.push_back(std::shared_ptr<TwoBodyAOInt>(ints[0]->clone()));
    }

    auto bs1_equiv_bs2 = (bs1 == bs2);
    auto bs3_equiv_bs4 = (bs3 == bs4);

#pragma omp parallel for collapse(2) schedule(guided) num_threads(nthreads_)
    for (size_t M = 0; M < bs1->nshell(); M++) {
        for (size_t N = 0; N < bs2->nshell(); N++) {
            if (bs1_equiv_bs2 && N < M) continue; // Only loop over unique shells

            const auto numM = bs1->shell(M).nfunction();
            const auto numN = bs2->shell(N).nfunction();
            const auto index_M = bs1->shell(M).function_index();
            const auto index_N = bs2->shell(N).function_index();

            for (size_t P = 0; P < bs3->nshell(); P++) {
                for (size_t Q = 0; Q < bs4->nshell(); Q++) {
                    if (bs3_equiv_bs4 && Q < P) continue; // Only loop over unique shells

                    size_t rank = 0;
#ifdef _OPENMP
                    rank = omp_get_thread_num();
#endif
                    const auto numP = bs3->shell(P).nfunction();
                    const auto numQ = bs4->shell(Q).nfunction();
                    const auto index_P = bs3->shell(P).function_index();
                    const auto index_Q = bs4->shell(Q).function_index();

                    ints[rank]->compute_shell(M, N, P, Q);
                    const auto *ints_buff = ints[rank]->buffers()[0];

                    if (bs1_equiv_bs2 && M != N && bs3_equiv_bs4 && P != Q) {
                        for (size_t m = 0, idx = 0; m < numM; m++) {
                            const auto fxnM = index_M + m;
                            for (size_t n = 0; n < numN; n++) {
                                const auto fxnN = index_N + n;
                                for (size_t p = 0; p < numP; p++) {
                                    const auto fxnP = index_P + p;

                                    double* targetMNPQ = &(*GAO)(fxnM, fxnN, fxnP, index_Q); // (mn|pq)
                                    double* targetNMQP = &(*GAO)(fxnN, fxnM, index_Q, fxnP); // (nm|qp)
                                    double* targetNMPQ = &(*GAO)(fxnN, fxnM, fxnP, index_Q); // (nm|pq)
                                    double* targetMNQP = &(*GAO)(fxnM, fxnN, index_Q, fxnP); // (mn|qp)

                                    for (size_t q = 0; q < numQ; q++, idx++) {
                                        // const auto fxnQ = index_Q + q;

                                        *targetMNPQ = *targetNMQP =
                                            *targetNMPQ = *targetMNQP = ints_buff[idx];

                                        targetMNPQ++;
                                        targetNMQP += (*GAO).dim(3);
                                        targetNMPQ++;
                                        targetMNQP += (*GAO).dim(3);
                                    }
                                }
                            }
                        }
                    } else if (bs1_equiv_bs2 && M != N) {
                        for (size_t m = 0, idx = 0; m < numM; m++) {
                            const auto fxnM = index_M + m;
                            for (size_t n = 0; n < numN; n++) {
                                const auto fxnN = index_N + n;
                                for (size_t p = 0; p < numP; p++) {
                                    const auto fxnP = index_P + p;

                                    double* targetMN = &(*GAO)(fxnM, fxnN, fxnP, index_Q); // (mn|pq)
                                    double* targetNM = &(*GAO)(fxnN, fxnM, fxnP, index_Q); // (nm|pq)

                                    for (size_t q = 0; q < numQ; q++, idx++) {
                                        *targetMN = *targetNM = ints_buff[idx];

                                        targetMN++;
                                        targetNM++;
                                    }
                                }
                            }
                        }
                    } else if (bs3_equiv_bs4 && P != Q) {
                        for (size_t m = 0, idx = 0; m < numM; m++) {
                            const auto fxnM = index_M + m;
                            for (size_t n = 0; n < numN; n++) {
                                const auto fxnN = index_N + n;
                                for (size_t p = 0; p < numP; p++) {
                                    const auto fxnP = index_P + p;

                                    double* targetPQ = &(*GAO)(fxnM, fxnN, fxnP, index_Q); // (mn|pq)
                                    double* targetQP = &(*GAO)(fxnM, fxnN, index_Q, fxnP); // (mn|qp)

                                    for (size_t q = 0; q < numQ; q++, idx++) {
                                        *targetPQ = *targetQP = ints_buff[idx];

                                        targetPQ++;
                                        targetQP += (*GAO).dim(3);
                                    }
                                }
                            }
                        }
                    } else {
                        for (size_t m = 0, idx = 0; m < numM; m++) {
                            const auto fxnM = index_M + m;
                            for (size_t n = 0; n < numN; n++) {
                                const auto fxnN = index_N + n;
                                for (size_t p = 0; p < numP; p++) {
                                    const auto fxnP = index_P + p;

                                    double* target = &(*GAO)(fxnM, fxnN, fxnP, index_Q); // (mn|pq)

                                    for (size_t q = 0; q < numQ; q++, idx++) {
                                        *target = ints_buff[idx];

                                        *target++;
                                    }
                                }
                            }
                        }
                    }
                } // bs4
            } // bs3
        } // bs2
    } // bs1
}

void MP2F12::three_index_ao_computer(const std::string& int_type, einsums::Tensor<double, 3> *Bpq,
                                     std::shared_ptr<BasisSet> bs1, std::shared_ptr<BasisSet> bs2)
{
    std::shared_ptr<BasisSet> zero(BasisSet::zero_ao_basis_set());
    std::shared_ptr<IntegralFactory> intf(new IntegralFactory(DFBS_, zero, bs1, bs2));

    std::vector<std::shared_ptr<TwoBodyAOInt>> ints;
    if ( int_type == "F" ){
        ints.push_back(std::shared_ptr<TwoBodyAOInt>(intf->f12(cgtg_)));
    } else if ( int_type == "FG" ){
        ints.push_back(std::shared_ptr<TwoBodyAOInt>(intf->f12g12(cgtg_)));
    } else if ( int_type == "F2" ){
        ints.push_back(std::shared_ptr<TwoBodyAOInt>(intf->f12_squared(cgtg_)));
    } else if ( int_type == "Uf" ){
        ints.push_back(std::shared_ptr<TwoBodyAOInt>(intf->f12_double_commutator(cgtg_)));
    } else {
        ints.push_back(std::shared_ptr<TwoBodyAOInt>(intf->eri()));
    }

    // Make ints vectors
    for (size_t thread = 1; thread < nthreads_; thread++) {
        ints.push_back(std::shared_ptr<TwoBodyAOInt>(ints[0]->clone()));
    }

    auto bs1_equiv_bs2 = bs1 == bs2;

#pragma omp parallel for collapse(3) schedule(guided) num_threads(nthreads_)
    for (size_t B = 0; B < DFBS_->nshell(); B++) {
        for (size_t P = 0; P < bs1->nshell(); P++) {
            for (size_t Q = 0; Q < bs2->nshell(); Q++) {
                if (bs1_equiv_bs2 && Q < P) continue; // Only loop over unique shells

                size_t rank = 0;
#ifdef _OPENMP
                rank = omp_get_thread_num();
#endif
                const auto numB = DFBS_->shell(B).nfunction();
                const auto numP = bs1->shell(P).nfunction();
                const auto numQ = bs2->shell(Q).nfunction();
                const auto index_B = DFBS_->shell(B).function_index();
                const auto index_P = bs1->shell(P).function_index();
                const auto index_Q = bs2->shell(Q).function_index();

                ints[rank]->compute_shell(B, 0, P, Q);
                const auto *ints_buff = ints[rank]->buffers()[0];

                if (bs1_equiv_bs2 && P != Q) {
                    for (size_t b = 0, idx = 0; b < numB; b++) {
                        const auto fxnB = index_B + b;
                        for (size_t p = 0; p < numP; p++) {
                            const auto fxnP = index_P + p;

                            double* targetPQ = &(*Bpq)(fxnB, fxnP, index_Q);
                            double* targetQP = &(*Bpq)(fxnB, index_Q, fxnP);

                            for (size_t q = 0; q < numQ; q++, idx++) {
                                *targetPQ = *targetQP = ints_buff[idx];

                                targetPQ++;
                                targetQP += (*Bpq).dim(2);
                            }
                        }
                    }
                } else {
                    for (size_t b = 0, index = 0; b < numB; b++) {
                        const auto B_fxns = index_B + b;
                        for (size_t p = 0; p < numP; p++) {
                            const auto P_fxns = index_P + p;

                            double* target = &(*Bpq)(B_fxns, P_fxns, index_Q);

                            for (size_t q = 0; q < numQ; q++, index++) {
                                *target = ints_buff[index];

                                target++;
                            }
                        }
                    }
                }
            }
        }
    }
}

void MP2F12::form_oeints(einsums::Tensor<double, 2> *h)
{ 
    using namespace einsums;

    outfile->Printf("   One-Electron Integrals\n");

    std::vector<char> order = {'O', 'O',
                               'O', 'C',
                               'C', 'C'};

    for (int idx = 0; idx < (order.size()/2); idx++) {
        const int i = idx * 2;
        const int o1 = ( order[i]  == 'C') ? 1 : 0;
        const int o2 = (order[i+1] == 'C') ? 1 : 0;
        const int nmo1 = (o1 == 1) ? ncabs_ : nobs_;
        const int nmo2 = (o2 == 1) ? ncabs_ : nobs_;
        const int M = (o1 == 1) ? nobs_ : 0;
        const int N = (o2 == 1) ? nobs_ : 0;

        auto mints = reference_wavefunction_->mintshelper();

        // Transform OEI AO Matrix into OEI MO Matrix
        auto t_mo = std::make_shared<Matrix>("MO-based T Integral", nmo1, nmo2);
        auto v_mo = std::make_shared<Matrix>("MO-based V Integral", nmo1, nmo2);
        {
            auto bs1 = bs_[o1].basisset();
            auto bs2 = bs_[o2].basisset();
            auto C1 = bs_[o1].C();
            auto C2 = bs_[o2].C();
            auto t_ao = mints->ao_kinetic(bs1, bs2);
            auto v_ao = mints->ao_potential(bs1, bs2);
            t_mo->transform(C1, t_ao, C2);
            v_mo->transform(C1, v_ao, C2);
            t_ao.reset();
            v_ao.reset();
        }

        // Stitch into OEI
        {
            TensorView<double, 2> h_mn{*h, Dim<2>{nmo1, nmo2}, Offset<2>{M, N}};
            for (int m = 0; m < nmo1; m++) {
                for (int n = 0; n < nmo2; n++) {
                    h_mn(m, n) = t_mo->get(m, n) + v_mo->get(m, n);
                }
            }

            if ( o1 != o2 ) {
                TensorView<double, 2> h_nm{*h, Dim<2>{nmo2, nmo1}, Offset<2>{N, M}};
                for (int n = 0; n < nmo2; n++) {
                    for (int m = 0; m < nmo1; m++) {
                        h_nm(n, m) = t_mo->get(m, n) + v_mo->get(m, n);
                    }
                }
            } // end of if statement
        }
    } // end for loop
}

void MP2F12::form_teints(const std::string& int_type, einsums::Tensor<double, 4> *ERI, std::vector<char> order)
{
    using namespace einsums;
    using namespace tensor_algebra;
    using namespace tensor_algebra::index;

    bool use_offset = true;
    if (order.size() == 4) {
        use_offset = false;
    }

    // (PQ|RS)
    for (int idx = 0; idx < (order.size()/4); idx++) {
        const int i = idx * 4;
        const int o1 = ( order[i]  == 'C') ? 1 : 0;
        const int o2 = (order[i+1] == 'C') ? 1 : 0;
        const int o3 = (order[i+2] == 'C') ? 1 : 0;
        const int o4 = (order[i+3] == 'C') ? 1 : 0;
        const auto nbf1 = bs_[o1].basisset()->nbf();
        const auto nbf2 = bs_[o2].basisset()->nbf();
        const auto nbf3 = bs_[o3].basisset()->nbf();
        const auto nbf4 = bs_[o4].basisset()->nbf();

        // Create ERI AO Tensor
        auto GAO = std::make_unique<Tensor<double, 4>>("ERI AO", nbf1, nbf2, nbf3, nbf4);
        timer_on("Two-Body AO Integrals");
        two_body_ao_computer(int_type, GAO.get(), bs_[o1].basisset(), bs_[o2].basisset(),
                                                  bs_[o3].basisset(), bs_[o4].basisset());
        timer_off("Two-Body AO Integrals");
	
        // Convert all Psi4 C Matrices to einsums Tensor<double, 2>
        const auto nmo1 = (o1) ? ncabs_ : ( order[i] == 'O' ) ? nobs_ : nocc_;
        const auto nmo2 = (o2) ? ncabs_ : (order[i+1] == 'O') ? nobs_ : nocc_;
        const auto nmo3 = (o3) ? ncabs_ : (order[i+2] == 'O') ? nobs_ : nocc_;
        const auto nmo4 = (o4) ? ncabs_ : (order[i+3] == 'O') ? nobs_ : nocc_;
	
        // Transform ERI AO Tensor to ERI MO Tensor
        auto PRQS = std::make_unique<Tensor<double, 4>>("PRQS", nmo1, nmo3, nmo2, nmo4);
        timer_on("MO Transformation");
        {
            // C1
            auto C1 = std::make_unique<Tensor<double, 2>>("C1", nbf1, nmo1);
            convert_C(C1.get(), bs_[o1], nbf1, nmo1);
            auto Pqrs = std::make_unique<Tensor<double, 4>>("Pqrs", nmo1, nbf2, nbf3, nbf4);
            einsum(Indices{P, q, r, s}, &Pqrs, Indices{p, q, r, s}, GAO, Indices{p, P}, C1);
            GAO.reset();
            C1.reset();

            // Sort
            auto rsPq = std::make_unique<Tensor<double, 4>>("rsPq", nbf3, nbf4, nmo1, nbf2);
            sort(Indices{r, s, P, q}, &rsPq, Indices{P, q, r, s}, Pqrs);
            Pqrs.reset();

            // C3
            auto C3 = std::make_unique<Tensor<double, 2>>("C3", nbf3, nmo3);
            convert_C(C3.get(), bs_[o3], nbf3, nmo3);
            auto RsPq = std::make_unique<Tensor<double, 4>>("RsPq", nmo3, nbf4, nmo1, nbf2);
            einsum(Indices{R, s, P, q}, &RsPq, Indices{r, s, P, q}, rsPq, Indices{r, R}, C3);
            rsPq.reset();
            C3.reset();

            // C2
            auto C2 = std::make_unique<Tensor<double, 2>>("C2", nbf2, nmo2);
            convert_C(C2.get(), bs_[o2], nbf2, nmo2);
            auto RsPQ = std::make_unique<Tensor<double, 4>>("RsPQ", nmo3, nbf4, nmo1, nmo2);
            einsum(Indices{R, s, P, Q}, &RsPQ, Indices{R, s, P, q}, RsPq, Indices{q, Q}, C2);
            RsPq.reset();
            C2.reset();

            // Sort
            auto PQRs = std::make_unique<Tensor<double, 4>>("PQRs", nmo1, nmo2, nmo3, nbf4);
            sort(Indices{P, Q, R, s}, &PQRs, Indices{R, s, P, Q}, RsPQ);
            RsPQ.reset();

            // C4
            auto C4 = std::make_unique<Tensor<double, 2>>("C4", nbf4, nmo4);
            convert_C(C4.get(), bs_[o4], nbf4, nmo4);
            auto PQRS = std::make_unique<Tensor<double, 4>>("PQRS", nmo1, nmo2, nmo3, nmo4);
            einsum(Indices{P, Q, R, index::S}, &PQRS, Indices{P, Q, R, s}, PQRs, Indices{s, index::S}, C4);
            PQRs.reset();
            C4.reset();

            // Switch from Chem to Phys
            sort(Indices{P, R, Q, index::S}, &PRQS, Indices{P, Q, R, index::S}, PQRS);
            PQRS.reset();
        }
        timer_off("MO Transformation");

        // Stitch into ERI Tensor
        timer_on("Set in ERI");
        {
            const auto off1 = (use_offset && o1) ? nobs_ : 0;
            const auto off2 = (use_offset && o2) ? nobs_ : 0;
            const auto off3 = (use_offset && o3) ? nobs_ : 0;
            const auto off4 = (use_offset && o4) ? nobs_ : 0;

            TensorView<double, 4> ERI_PRQS{*ERI, Dim<4>{nmo1, nmo3, nmo2, nmo4}, Offset<4>{off1, off3, off2, off4}};
            set_ERI(ERI_PRQS, PRQS.get());

            if (nbf2 != nbf1 && nbf2 != nbf3 && nbf2 != nbf4 && int_type == "F") {
                Tensor<double, 4> RPSQ{"RPSQ", nmo3, nmo1, nmo4, nmo2};
                sort(Indices{R, P, index::S, Q}, &RPSQ, Indices{P, R, Q, index::S}, PRQS);
                PRQS.reset();
                TensorView<double, 4> ERI_RPSQ{*ERI, Dim<4>{nmo3, nmo1, nmo4, nmo2}, Offset<4>{off3, off1, off4, off2}};
                set_ERI(ERI_RPSQ, &RPSQ);
            } // end of if statement

            if (nbf2 != nbf1 && nbf2 != nbf3 && nbf2 != nbf4 && int_type == "J") {
                Tensor<double, 4> QSPR{"QSPR", nmo2, nmo4, nmo1, nmo3};
                sort(Indices{Q, index::S, P, R}, &QSPR, Indices{P, R, Q, index::S}, PRQS);
                PRQS.reset();
                TensorView<double, 4> ERI_QSPR{*ERI, Dim<4>{nmo2, nmo4, nmo1, nmo3}, Offset<4>{off2, off4, off1, off3}};
                set_ERI(ERI_QSPR, &QSPR);
            } // end of if statement

            if (nbf4 != nbf1 && nbf4 != nbf2 && nbf4 != nbf3 && int_type == "K") {
                Tensor<double, 4> SQRP{"SQRP", nmo4, nmo2, nmo3, nmo1};
                sort(Indices{index::S, Q, R, P}, &SQRP, Indices{P, R, Q, index::S}, PRQS);
                PRQS.reset();
                TensorView<double, 4> ERI_SQRP{*ERI, Dim<4>{nmo4, nmo2, nmo3, nmo1}, Offset<4>{off4, off2, off3, off1}};
                set_ERI(ERI_SQRP, &SQRP);
            } // end of if statement
        }
        timer_off("Set in ERI");
    } // end of for loop
}

void MP2F12::form_metric_ints(einsums::Tensor<double, 3> *DF_ERI, bool is_fock)
{
    using namespace einsums;
    using namespace tensor_algebra;
    using namespace tensor_algebra::index;

    std::vector<char> order;

    if (is_fock) {
        order = {'O', 'O',
                 'O', 'C',
                 'C', 'O',
                 'C', 'C'};
    } else {
        order = {'o', 'O',
                 'o', 'C'};
    }
    
    for (int idx = 0; idx < (order.size()/2); idx++) {
        const int i = idx * 2;
        const int o1 = ( order[i]  == 'C') ? 1 : 0;
        const int o2 = (order[i+1] == 'C') ? 1 : 0;
        const auto nbf1 = bs_[o1].basisset()->nbf();
        const auto nbf2 = bs_[o2].basisset()->nbf();

        auto Bpq = std::make_unique<Tensor<double, 3>>("Metric AO", naux_, nbf1, nbf2);
        three_index_ao_computer("G", Bpq.get(), bs_[o1].basisset(), bs_[o2].basisset());

        const auto nmo1 = (o1) ? ncabs_ : ( order[i] == 'O' ) ? nobs_ : nocc_;
        const auto nmo2 = (o2) ? ncabs_ : (order[i+1] == 'O') ? nobs_ : nocc_;

        auto BPQ = std::make_unique<Tensor<double, 3>>("BPQ", naux_, nmo1, nmo2);
        timer_on("MO Transformation");
        {
            // C2
            auto C2 = std::make_unique<Tensor<double, 2>>("C2", nbf2, nmo2);
            convert_C(C2.get(), bs_[o2], nbf2, nmo2);
            Tensor<double, 3> BpQ{"BpQ", naux_, nbf1, nmo2};
            einsum(Indices{B, p, Q}, &BpQ, Indices{B, p, q}, Bpq, Indices{q, Q}, C2);
            C2.reset();

            // Sort
            Tensor<double, 3> BQp{"BQp", naux_, nmo2, nbf1};
            sort(Indices{B, Q, p}, &BQp, Indices{B, p, Q}, BpQ);

            // C1
            auto C1 = std::make_unique<Tensor<double, 2>>("C1", nbf1, nmo1);
            convert_C(C1.get(), bs_[o1], nbf1, nmo1);
            Tensor<double, 3> BQP{"BQP", naux_, nmo2, nmo1};
            einsum(Indices{B, Q, P}, &BQP, Indices{B, Q, p}, BQp, Indices{p, P}, C1);
            C1.reset();

            // Sort
            sort(Indices{B, P, Q}, &BPQ, Indices{B, Q, P}, BQP);
        }
        timer_off("MO Transformation");

        auto APQ = std::make_unique<Tensor<double, 3>>("APQ", naux_, nmo1, nmo2);
        {
            auto metric = std::make_shared<FittingMetric>(DFBS_, true);
            metric->form_full_eig_inverse(1.0e-12);
            SharedMatrix Jm12 = metric->get_metric();

            Tensor<double, 2> AB{"JinvAB", naux_, naux_};
            for (size_t A = 0; A < naux_; A++) {
                for (size_t B = 0; B < naux_; B++) {
                    AB(A, B) = Jm12->get(A, B);
                }
            }

            einsum(Indices{A, P, Q}, &APQ, Indices{A, B}, AB, Indices{B, P, Q}, BPQ);
        }
        BPQ.reset();

        {
            const auto R = (o1) ? nobs_ : 0;
            const auto S = (o2) ? nobs_ : 0;

            TensorView<double, 3> ERI_APQ{*DF_ERI, Dim<3>{naux_, nmo1, nmo2}, Offset<3>{0, R, S}};
            set_ERI(ERI_APQ, APQ.get());
        }
    } // end of for loop
}

void MP2F12::form_oper_ints(const std::string& int_type, einsums::Tensor<double, 3> *DF_ERI)
{
    using namespace einsums;
    using namespace tensor_algebra;
    using namespace tensor_algebra::index;

    std::vector<char> order = {'o', 'o'};
    if ( int_type != "Uf" || int_type != "FG" ) {
        order = {'o', 'O',
                 'o', 'C'};
    }

    for (int idx = 0; idx < (order.size()/2); idx++) {
        const int i = idx * 2;
        const int o1 = ( order[i]  == 'C') ? 1 : 0;
        const int o2 = (order[i+1] == 'C') ? 1 : 0;
        const auto nbf1 = bs_[o1].basisset()->nbf();
        const auto nbf2 = bs_[o2].basisset()->nbf();

        auto Bpq = std::make_unique<Tensor<double, 3>>("(B|R|pq) AO", naux_, nbf1, nbf2);
        timer_on("Three-Center AO Integrals");
        three_index_ao_computer(int_type, Bpq.get(), bs_[o1].basisset(), bs_[o2].basisset());
        timer_off("Three-Center AO Integrals");

        // Convert all Psi4 C Matrices to einsums Tensor<double, 2>
        const auto nmo1 = (o1) ? ncabs_ : ( order[i] == 'O' ) ? nobs_ : nocc_;
        const auto nmo2 = (o2) ? ncabs_ : (order[i+1] == 'O') ? nobs_ : nocc_;

        auto BPQ = std::make_unique<Tensor<double, 3>>("BPQ", naux_, nmo1, nmo2);
        timer_on("MO Transformation");
        {
            // C2
            auto C2 = std::make_unique<Tensor<double, 2>>("C2", nbf2, nmo2);
            convert_C(C2.get(), bs_[o2], nbf2, nmo2);
            Tensor<double, 3> BpQ{"BpQ", naux_, nbf1, nmo2};
            einsum(Indices{B, p, Q}, &BpQ, Indices{B, p, q}, Bpq, Indices{q, Q}, C2);
            C2.reset();

            // Sort
            Tensor<double, 3> BQp{"BQp", naux_, nmo2, nbf1};
            sort(Indices{B, Q, p}, &BQp, Indices{B, p, Q}, BpQ);

            // C1
            auto C1 = std::make_unique<Tensor<double, 2>>("C1", nbf1, nmo1);
            convert_C(C1.get(), bs_[o1], nbf1, nmo1);
            Tensor<double, 3> BQP{"BQP", naux_, nmo2, nmo1};
            einsum(Indices{B, Q, P}, &BQP, Indices{B, Q, p}, BQp, Indices{p, P}, C1);
            C1.reset();

            // Sort
            sort(Indices{B, P, Q}, &BPQ, Indices{B, Q, P}, BQP);
        }
        timer_off("MO Transformation");

        timer_on("Set in ERI");
        {
            const auto R = (o1) ? nobs_ : 0;
            const auto S = (o2) ? nobs_ : 0;

            TensorView<double, 3> ERI_BPQ{*DF_ERI, Dim<3>{naux_, nmo1, nmo2}, Offset<3>{0, R, S}};
            set_ERI(ERI_BPQ, BPQ.get());
        }
        timer_off("Set in ERI");
    } // end of for loop
}

void MP2F12::form_oper_ints(const std::string& int_type, einsums::Tensor<double, 2> *DF_ERI)
{
    using namespace einsums;
    using namespace tensor_algebra;
    using namespace tensor_algebra::index;

    std::shared_ptr<BasisSet> zero(BasisSet::zero_ao_basis_set());

    std::shared_ptr<IntegralFactory> intf(new IntegralFactory(DFBS_, zero, DFBS_, zero));

    std::vector<std::shared_ptr<TwoBodyAOInt>> ints;
    if ( int_type == "F" ){
        ints.push_back(std::shared_ptr<TwoBodyAOInt>(intf->f12(cgtg_)));
    } else if ( int_type == "FG" ){
        ints.push_back(std::shared_ptr<TwoBodyAOInt>(intf->f12g12(cgtg_)));
    } else if ( int_type == "F2" ){
        ints.push_back(std::shared_ptr<TwoBodyAOInt>(intf->f12_squared(cgtg_)));
    } else if ( int_type == "Uf" ){
        ints.push_back(std::shared_ptr<TwoBodyAOInt>(intf->f12_double_commutator(cgtg_)));
    } else {
        ints.push_back(std::shared_ptr<TwoBodyAOInt>(intf->eri()));
    }

    // Make ints vectors
    for (size_t thread = 1; thread < nthreads_; thread++) {
        ints.push_back(std::shared_ptr<TwoBodyAOInt>(ints[0]->clone()));
    }

#pragma omp parallel for collapse(2) schedule(guided) num_threads(nthreads_)
    for (size_t A = 0; A < DFBS_->nshell(); A++) {
        for (size_t B = 0; B < DFBS_->nshell(); B++) {
            size_t rank = 0;
#ifdef _OPENMP
            rank = omp_get_thread_num();
#endif
            const size_t numA = DFBS_->shell(A).nfunction();
            const size_t numB = DFBS_->shell(B).nfunction();
            const size_t index_A = DFBS_->shell(A).function_index();
            const size_t index_B = DFBS_->shell(B).function_index();

            ints[rank]->compute_shell(A, 0, B, 0);
            const auto *ints_buff = ints[rank]->buffers()[0];

            size_t index = 0;
            for (size_t a = 0; a < numA; a++) {
                for (size_t b = 0; b < numB; b++) {
                    (*DF_ERI)(index_A + a,
                          index_B + b) = ints_buff[index++];
                }
            }
        }
    }
}

void MP2F12::form_df_teints(const std::string& int_type, einsums::Tensor<double, 4> *ERI,
                         einsums::Tensor<double, 3> *J_inv_AB, std::vector<char> order)
{
    using namespace einsums;
    using namespace tensor_algebra;
    using namespace tensor_algebra::index;

    bool use_offset = true;
    if (order.size() == 4) {
        use_offset = false;
    }

    // In (PQ|RS) ordering
    for (int idx = 0; idx < (order.size()/4); idx++) {
        const int i = idx * 4;
        const auto nmo1 = ( order[ i ] == 'C' ) ? ncabs_ : ( order[ i ] == 'O' ) ? nobs_ : nocc_;
        const auto nmo2 = ( order[i+1] == 'C' ) ? ncabs_ : ( order[i+1] == 'O' ) ? nobs_ : nocc_;
        const auto nmo3 = ( order[i+2] == 'C' ) ? ncabs_ : ( order[i+2] == 'O' ) ? nobs_ : nocc_;
        const auto nmo4 = ( order[i+3] == 'C' ) ? ncabs_ : ( order[i+3] == 'O' ) ? nobs_ : nocc_;
        auto off1 = ( order[ i ] == 'C' ) ? nobs_ : 0;
        auto off2 = ( order[i+1] == 'C' ) ? nobs_ : 0;
        auto off3 = ( order[i+2] == 'C' ) ? nobs_ : 0;
        auto off4 = ( order[i+3] == 'C' ) ? nobs_ : 0;

        auto phys_robust = std::make_unique<Tensor<double, 4>>("<PR|F12|QS> MO", nmo1, nmo3, nmo2, nmo4);
        {
            Tensor<double, 4> chem_robust("(PQ|F12|RS) MO", nmo1, nmo2, nmo3, nmo4);

            timer_on("Robust DF Procedure");
            auto ARPQ = std::make_unique<Tensor<double, 3>>("(A|R|PQ) MO", naux_, nocc_, nri_);
            form_oper_ints(int_type, ARPQ.get());

            // Term 1
            Tensor left_metric  = (*J_inv_AB)(All, Range{off1, nmo1 + off1}, Range{off2, nmo2 + off2});
            Tensor right_oper = (*ARPQ)(All, Range{off3, nmo3 + off3}, Range{off4, nmo4 + off4});
            einsum(Indices{p, q, r, s}, &chem_robust, Indices{A, p, q}, left_metric, Indices{A, r, s}, right_oper);

            if ( int_type != "G" ) {
                // Term 2
                Tensor right_metric = (*J_inv_AB)(All, Range{off3, nmo3 + off3}, Range{off4, nmo4 + off4});
                {
                    Tensor left_oper  = (*ARPQ)(All, Range{off1, nmo1 + off1}, Range{off2, nmo2 + off2});
                    einsum(1.0, Indices{p, q, r, s}, &chem_robust,
                           1.0, Indices{A, p, q}, left_oper, Indices{A, r, s}, right_metric);
                }
                ARPQ.reset();

                // Term 3
                {
                    auto ARB = std::make_unique<Tensor<double, 2>>("(A|R|B) MO", naux_, naux_);
                    form_oper_ints(int_type, ARB.get());

                    Tensor<double, 3> tmp{"Temp", naux_, nmo3, nmo4};
                    einsum(Indices{A, r, s}, &tmp, Indices{A, B}, ARB, Indices{B, r, s}, right_metric);
                    ARB.reset();
                    einsum(1.0, Indices{p, q, r, s}, &chem_robust,
                           -1.0, Indices{A, p, q}, left_metric, Indices{A, r, s}, tmp);
                }
            }
            timer_off("Robust DF Procedure");

            // Switch to <PR|QS> ordering
            sort(Indices{p, r, q, s}, &phys_robust, Indices{p, q, r, s}, chem_robust);
        }

        timer_on("Set in ERI");
        {
            if (!use_offset){
                off1 = 0;
                off2 = 0;
                off3 = 0;
                off4 = 0;
            }

            TensorView<double, 4> ERI_PRQS{(*ERI), Dim<4>{nmo1, nmo3, nmo2, nmo4}, Offset<4>{off1, off3, off2, off4}};
            set_ERI(ERI_PRQS, phys_robust.get());
        }
        timer_off("Set in ERI");
    } // end of for loop
}

////////////////////////////////
//* Disk Algorithm (CONV/DF) *//
////////////////////////////////

void DiskMP2F12::set_ERI(einsums::DiskView<double, 2, 4>& ERI_Slice, einsums::TensorView<double, 2>& Slice)
{
    using namespace einsums;

    const auto dim1 = Slice.dim(0);
    const auto dim2 = Slice.dim(1);

    for (int p = 0; p < dim1; p++){
        for (int q = 0; q < dim2; q++){
            ERI_Slice(p, q) = Slice(p, q);
        }
    }
}

void DiskMP2F12::form_oeints(einsums::DiskTensor<double, 2> *h)
{
    using namespace einsums;

    outfile->Printf("   One-Electron Integrals\n");

    std::vector<char> order = {'O', 'O',
                               'O', 'C',
                               'C', 'C'};

    for (int idx = 0; idx < (order.size()/2); idx++) {
        int i = idx * 2;
        int o1 = ( order[i]  == 'C') ? 1 : 0;
        int o2 = (order[i+1] == 'C') ? 1 : 0;
        int nmo1 = (o1 == 1) ? ncabs_ : nobs_;
        int nmo2 = (o2 == 1) ? ncabs_ : nobs_;
        int M = (o1 == 1) ? nobs_ : 0;
        int N = (o2 == 1) ? nobs_ : 0;

        auto mints = reference_wavefunction_->mintshelper();

        // Transform OEI AO Matrix into OEI MO Matrix
        auto t_mo = std::make_shared<Matrix>("MO-based T Integral", nmo1, nmo2);
        auto v_mo = std::make_shared<Matrix>("MO-based V Integral", nmo1, nmo2);
        {
            auto bs1 = bs_[o1].basisset();
            auto bs2 = bs_[o2].basisset();
            auto C1 = bs_[o1].C();
            auto C2 = bs_[o2].C();
            auto t_ao = mints->ao_kinetic(bs1, bs2);
            auto v_ao = mints->ao_potential(bs1, bs2);
            t_mo->transform(C1, t_ao, C2);
            v_mo->transform(C1, v_ao, C2);
            t_ao.reset();
            v_ao.reset();
        }

        // Stitch into OEI
        {
            for (int m = 0; m < nmo1; m++) {
                auto h_view = (*h)(m + M, All);
                for (int n = 0; n < nmo2; n++) {
                    h_view(n + N) = t_mo->get(m, n) + v_mo->get(m, n);
                }
            }

            if ( o1 != o2 ) {
                for (int n = 0; n < nmo2; n++) {
                    auto h_view = (*h)(n + N, All);
                    for (int m = 0; m < nmo1; m++) {
                        h_view(m + M) = t_mo->get(m, n) + v_mo->get(m, n);
                    }
                }
            } // end of if statement
        }
    } // end for loop
}

void DiskMP2F12::form_teints(const std::string& int_type, einsums::DiskTensor<double, 4> *ERI)
{
    using namespace einsums;
    using namespace tensor_algebra;
    using namespace tensor_algebra::index;

    // In (PQ|RS) ordering
    std::vector<char> order = {'o', 'o', 'o', 'o'};
    if ( int_type == "F" ) {
        order = {'o', 'O', 'o', 'O',
                 'o', 'O', 'o', 'C',
                 'o', 'C', 'o', 'C'};
    } else if ( int_type == "F2" ) {
        order = {'o', 'o', 'o', 'O',
                 'o', 'o', 'o', 'C'};
    } else if ( int_type == "G" ) {
        order = {'o', 'O', 'o', 'O',
                 'o', 'O', 'o', 'C'};
    } else if ( int_type == "J" ) {
        order = {'O', 'O', 'o', 'o',
                 'O', 'C', 'o', 'o',
                 'C', 'C', 'o', 'o'};
    } else if ( int_type == "K" ) {
        order = {'O', 'o', 'o', 'O',
                 'O', 'o', 'o', 'C',
                 'C', 'o', 'o', 'C'};
    }

    // (PQ|RS)
    for (int idx = 0; idx < (order.size()/4); idx++) {
        const int i = idx * 4;
        const int o1 = ( order[i]  == 'C') ? 1 : 0;
        const int o2 = (order[i+1] == 'C') ? 1 : 0;
        const int o3 = (order[i+2] == 'C') ? 1 : 0;
        const int o4 = (order[i+3] == 'C') ? 1 : 0;
        const auto nbf1 = bs_[o1].basisset()->nbf();
        const auto nbf2 = bs_[o2].basisset()->nbf();
        const auto nbf3 = bs_[o3].basisset()->nbf();
        const auto nbf4 = bs_[o4].basisset()->nbf();

        // Create ERI AO Tensor
        auto GAO = std::make_unique<Tensor<double, 4>>("ERI AO", nbf1, nbf2, nbf3, nbf4);
        timer_on("Two-Body AO Integrals");
        two_body_ao_computer(int_type, GAO.get(), bs_[o1].basisset(), bs_[o2].basisset(),
                                                  bs_[o3].basisset(), bs_[o4].basisset());
        timer_off("Two-Body AO Integrals");

        // Convert all Psi4 C Matrices to einsums Tensor<double, 2>
        const auto nmo1 = (o1) ? ncabs_ : ( order[i] == 'O' ) ? nobs_ : nocc_;
        const auto nmo2 = (o2) ? ncabs_ : (order[i+1] == 'O') ? nobs_ : nocc_;
        const auto nmo3 = (o3) ? ncabs_ : (order[i+2] == 'O') ? nobs_ : nocc_;
        const auto nmo4 = (o4) ? ncabs_ : (order[i+3] == 'O') ? nobs_ : nocc_;

        // Transform ERI AO Tensor to ERI MO Tensor
        auto PRQS = std::make_unique<Tensor<double, 4>>("PRQS", nmo1, nmo3, nmo2, nmo4);
        timer_on("MO Transformation");
        {
            // C1
            auto C1 = std::make_unique<Tensor<double, 2>>("C1", nbf1, nmo1);
            convert_C(C1.get(), bs_[o1], nbf1, nmo1);
            auto Pqrs = std::make_unique<Tensor<double, 4>>("Pqrs", nmo1, nbf2, nbf3, nbf4);
            einsum(Indices{P, q, r, s}, &Pqrs, Indices{p, q, r, s}, GAO, Indices{p, P}, C1);
            GAO.reset();
            C1.reset();

            // Sort
            auto rsPq = std::make_unique<Tensor<double, 4>>("rsPq", nbf3, nbf4, nmo1, nbf2);
            sort(Indices{r, s, P, q}, &rsPq, Indices{P, q, r, s}, Pqrs);
            Pqrs.reset();

            // C3
            auto C3 = std::make_unique<Tensor<double, 2>>("C3", nbf3, nmo3);
            convert_C(C3.get(), bs_[o3], nbf3, nmo3);
            auto RsPq = std::make_unique<Tensor<double, 4>>("RsPq", nmo3, nbf4, nmo1, nbf2);
            einsum(Indices{R, s, P, q}, &RsPq, Indices{r, s, P, q}, rsPq, Indices{r, R}, C3);
            rsPq.reset();
            C3.reset();

            // C2
            auto C2 = std::make_unique<Tensor<double, 2>>("C2", nbf2, nmo2);
            convert_C(C2.get(), bs_[o2], nbf2, nmo2);
            auto RsPQ = std::make_unique<Tensor<double, 4>>("RsPQ", nmo3, nbf4, nmo1, nmo2);
            einsum(Indices{R, s, P, Q}, &RsPQ, Indices{R, s, P, q}, RsPq, Indices{q, Q}, C2);
            RsPq.reset();
            C2.reset();

            // Sort
            auto PQRs = std::make_unique<Tensor<double, 4>>("PQRs", nmo1, nmo2, nmo3, nbf4);
            sort(Indices{P, Q, R, s}, &PQRs, Indices{R, s, P, Q}, RsPQ);
            RsPQ.reset();

            // C4
            auto C4 = std::make_unique<Tensor<double, 2>>("C4", nbf4, nmo4);
            convert_C(C4.get(), bs_[o4], nbf4, nmo4);
            auto PQRS = std::make_unique<Tensor<double, 4>>("PQRS", nmo1, nmo2, nmo3, nmo4);
            einsum(Indices{P, Q, R, index::S}, &PQRS, Indices{P, Q, R, s}, PQRs, Indices{s, index::S}, C4);
            PQRs.reset();
            C4.reset();

            // Switch from Chem to Phys
            sort(Indices{P, R, Q, index::S}, &PRQS, Indices{P, Q, R, index::S}, PQRS);
            PQRS.reset();
        }
        timer_off("MO Transformation");

        // Stitch into ERI Tensor
        timer_on("Set in ERI");
        {
            const auto off1 = (o1) ? nobs_ : 0;
            const auto off2 = (o2) ? nobs_ : 0;
            const auto off3 = (o3) ? nobs_ : 0;
            const auto off4 = (o4) ? nobs_ : 0;

            for (int p = 0; p < nmo1; p++) {
                for (int r = 0; r < nmo3; r++) {
                    auto ERI_PRQS = (*ERI)(p + off1, r + off3, Range{off2, off2 + nmo2}, Range{off4, off4 + nmo4});
                    auto PRQS_view = (*PRQS)(p, r, All, All);
                    set_ERI(ERI_PRQS, PRQS_view);
                }
            }

            if (nbf4 != nbf1 && nbf4 != nbf2 && nbf4 != nbf3 && int_type == "F") {
                Tensor<double, 4> RPSQ{"RPSQ", nmo3, nmo1, nmo4, nmo2};
                sort(Indices{R, P, index::S, Q}, &RPSQ, Indices{P, R, Q, index::S}, PRQS);
                PRQS.reset();
                for (int r = 0; r < nmo3; r++) {
                    for (int p = 0; p < nmo1; p++) {
                        auto ERI_rpSQ = (*ERI)(r + off3, p + off1, Range{off4, off4 + nmo4}, Range{off2, off2 + nmo2});
                        auto rpSQ_view = RPSQ(r, p, All, All);
                        set_ERI(ERI_rpSQ, rpSQ_view);
                    }
                }
            } // end of if statement

            if (nbf2 != nbf1 && nbf2 != nbf3 && nbf2 != nbf4 && int_type == "J") {
                Tensor<double, 4> QSPR{"QSPR", nmo2, nmo4, nmo1, nmo3};
                sort(Indices{Q, index::S, P, R}, &QSPR, Indices{P, R, Q, index::S}, PRQS);
                PRQS.reset();
                for (int q = 0; q < nmo2; q++) {
                    for (int s = 0; s < nmo4; s++) {
                        auto ERI_qsPR = (*ERI)(q + off2, s + off4, Range{off1, off1 + nmo1}, Range{off3, off3 + nmo3});
                        auto qsPR_view = QSPR(q, s, All, All);
                        set_ERI(ERI_qsPR, qsPR_view);
                    }
                }
            } // end of if statement

            if (nbf4 != nbf1 && nbf4 != nbf2 && nbf4 != nbf3 && int_type == "K") {
                Tensor<double, 4> SQRP{"SQRP", nmo4, nmo2, nmo3, nmo1};
                sort(Indices{index::S, Q, R, P}, &SQRP, Indices{P, R, Q, index::S}, PRQS);
                PRQS.reset();
                for (int s = 0; s < nmo4; s++) {
                    for (int q = 0; q < nmo2; q++) {
                        auto ERI_sqRP = (*ERI)(s + off4, q + off2, Range{off3, off3 + nmo3}, Range{off1, off1 + nmo1});
                        auto sqRP_view = SQRP(s, q, All, All);
                        set_ERI(ERI_sqRP, sqRP_view);
                    }
                }
            } // end of if statement
        }
    timer_off("Set in ERI");
    } // end of for loop
}

void DiskMP2F12::form_df_teints(const std::string& int_type, einsums::DiskTensor<double, 4> *ERI,
                                einsums::Tensor<double, 3> *J_inv_AB)
{
    using namespace einsums;
    using namespace tensor_algebra;
    using namespace tensor_algebra::index;

    // In (PQ|RS) ordering
    std::vector<char> order = {'o', 'o', 'o', 'o'};
    if ( int_type == "G" ) {
        order = {'o', 'O', 'o', 'O',
                 'o', 'O', 'o', 'C'};
    } else if ( int_type == "F" ) {
        order = {'o', 'O', 'o', 'O',
                 'o', 'O', 'o', 'C',
                 'o', 'C', 'o', 'O',
                 'o', 'C', 'o', 'C',};
    } else if ( int_type == "F2" ) {
        order = {'o', 'o', 'o', 'O',
                 'o', 'o', 'o', 'C',};
    }

    // (PQ|RS)
    for (int idx = 0; idx < (order.size()/4); idx++) {
        const int i = idx * 4;
        const auto nmo1 = ( order[ i ] == 'C' ) ? ncabs_ : ( order[ i ] == 'O' ) ? nobs_ : nocc_;
        const auto nmo2 = ( order[i+1] == 'C' ) ? ncabs_ : ( order[i+1] == 'O' ) ? nobs_ : nocc_;
        const auto nmo3 = ( order[i+2] == 'C' ) ? ncabs_ : ( order[i+2] == 'O' ) ? nobs_ : nocc_;
        const auto nmo4 = ( order[i+3] == 'C' ) ? ncabs_ : ( order[i+3] == 'O' ) ? nobs_ : nocc_;
        const auto off1 = ( order[ i ] == 'C' ) ? nobs_ : 0;
        const auto off2 = ( order[i+1] == 'C' ) ? nobs_ : 0;
        const auto off3 = ( order[i+2] == 'C' ) ? nobs_ : 0;
        const auto off4 = ( order[i+3] == 'C' ) ? nobs_ : 0;

        // Dunlap's robust density-fitting
        auto phys_robust = std::make_unique<Tensor<double, 4>>("<PR|F12|QS> MO", nmo1, nmo3, nmo2, nmo4);
        {
            Tensor<double, 4> chem_robust("(PQ|F12|RS) MO", nmo1, nmo2, nmo3, nmo4);

            timer_on("Robust DF Procedure");
            auto ARPQ = std::make_unique<Tensor<double, 3>>("(A|R|PQ) MO", naux_, nocc_, nri_);
            form_oper_ints(int_type, ARPQ.get());

            // Term 1
            Tensor left_metric  = (*J_inv_AB)(All, Range{off1, nmo1 + off1}, Range{off2, nmo2 + off2});
            Tensor right_oper = (*ARPQ)(All, Range{off3, nmo3 + off3}, Range{off4, nmo4 + off4});
            einsum(Indices{p, q, r, s}, &chem_robust, Indices{A, p, q}, left_metric, Indices{A, r, s}, right_oper);

            if ( int_type != "G" ) {
                // Term 2
                Tensor right_metric = (*J_inv_AB)(All, Range{off3, nmo3 + off3}, Range{off4, nmo4 + off4});
                {
                    Tensor left_oper  = (*ARPQ)(All, Range{off1, nmo1 + off1}, Range{off2, nmo2 + off2});
                    einsum(1.0, Indices{p, q, r, s}, &chem_robust,
                           1.0, Indices{A, p, q}, left_oper, Indices{A, r, s}, right_metric);
                }
                ARPQ.reset();

                // Term 3
                {
                    auto ARB = std::make_unique<Tensor<double, 2>>("(A|R|B) MO", naux_, naux_);
                    form_oper_ints(int_type, ARB.get());

                    Tensor<double, 3> tmp{"Temp", naux_, nmo3, nmo4};
                    einsum(Indices{A, r, s}, &tmp, Indices{A, B}, ARB, Indices{B, r, s}, right_metric);
                    ARB.reset();
                    einsum(1.0, Indices{p, q, r, s}, &chem_robust,
                           -1.0, Indices{A, p, q}, left_metric, Indices{A, r, s}, tmp);
                }
            }
            timer_off("Robust DF Procedure");

            // Switch to Phys Notation
            sort(Indices{p, r, q, s}, &phys_robust, Indices{p, q, r, s}, chem_robust);
        }

        // Stitch into ERI Tensor
        timer_on("Set in ERI");
        for (int p = off1; p < (off1 + nmo1); p++) {
            for (int r = off3; r < (off3 + nmo3); r++) {
                auto ERI_prQS = (*ERI)(p, r, Range{off2, off2 + nmo2}, Range{off4, off4 + nmo4});
                auto prQS_view = (*phys_robust)(p, r, All, All);
                set_ERI(ERI_prQS, prQS_view);
            }
        }
        timer_off("Set in ERI");
    } // end of for loop
}

}} // End namespaces
