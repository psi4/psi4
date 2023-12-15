/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2023 The Psi4 Developers.
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

#include "psi4/libpsi4util/PsiOutStream.h"

#include "einsums.hpp"

namespace psi { namespace mp2f12 {

void MP2F12::form_fock(einsums::Tensor<double, 2> *f, einsums::Tensor<double, 2> *k)
{
    using namespace einsums;
    using namespace tensor_algebra;
    using namespace tensor_algebra::index;

    form_oeints(f);

    Tensor Id = create_identity_tensor("I", nocc_, nocc_);
    {
        outfile->Printf("     Forming J\n");
        auto J = std::make_unique<Tensor<double, 4>>("Coulomb", nri_, nocc_, nri_, nocc_);
        form_teints("J", J.get(), {'O', 'o', 'O', 'o',
                                   'O', 'o', 'C', 'o',
                                   'C', 'o', 'C', 'o'});

        Tensor<double, 4> J_sorted{"pqiI", nri_, nri_, nocc_, nocc_};
        sort(Indices{p, q, i, I}, &J_sorted, Indices{p, i, q, I}, J);
        J.reset();
        einsum(1.0, Indices{p, q}, &(*f), 2.0, Indices{p, q, i, I}, J_sorted, Indices{i, I}, Id);
    }

    {
        outfile->Printf("     Forming K\n");
        auto K = std::make_unique<Tensor<double, 4>>("Exhange", nri_, nocc_, nocc_, nri_);
        form_teints("K", K.get(), {'O', 'o', 'o', 'O',
                                   'O', 'o', 'o', 'C',
                                   'C', 'o', 'o', 'C'});

        Tensor<double, 4> K_sorted{"pqiI", nri_, nri_, nocc_, nocc_};
        sort(Indices{p, q, i, I}, &K_sorted, Indices{p, i, I, q}, K);
        K.reset();
        einsum(Indices{p, q}, &(*k), Indices{p, q, i, I}, K_sorted, Indices{i, I}, Id);
    }

    tensor_algebra::element([](double const &val1, double const &val2)
                            -> double { return val1 - val2; },
                            &(*f), *k);
}

void MP2F12::form_df_fock(einsums::Tensor<double, 2> *f, einsums::Tensor<double, 2> *k)
{
    using namespace einsums;
    using namespace tensor_algebra;
    using namespace tensor_algebra::index;

    form_oeints(f);

    {
        auto Metric = std::make_unique<Tensor<double, 3>>("(B|PQ) MO", naux_, nri_, nri_);
        form_metric_ints(Metric.get(), true);
        auto Oper = std::make_unique<Tensor<double, 3>>("(B|PQ) MO", naux_, nocc_, nri_);
        form_oper_ints("G", Oper.get());

        {
            outfile->Printf("     Forming J\n");
            Tensor Id = create_identity_tensor("I", nocc_, nocc_);
            Tensor J_Metric = (*Metric)(Range{0, naux_}, Range{0, nri_}, Range{0, nri_});
            Tensor J_Oper = (*Oper)(Range{0, naux_}, Range{0, nocc_}, Range{0, nocc_});

            Tensor<double, 1> tmp{"B", naux_};
            einsum(Indices{B}, &tmp, Indices{B, i, j}, J_Oper, Indices{i, j}, Id);
            einsum(1.0, Indices{P, Q}, &(*f), 2.0, Indices{B, P, Q}, J_Metric, Indices{B}, tmp);
        }

        {
            outfile->Printf("     Forming K\n");
            Tensor K_Metric = (*Metric)(Range{0, naux_}, Range{0, nri_}, Range{0, nocc_});
            Tensor K_Oper = (*Oper)(Range{0, naux_}, Range{0, nocc_}, Range{0, nri_});

            Tensor<double, 3> tmp{"", naux_, nocc_, nri_};
            sort(Indices{B, i, P}, &tmp, Indices{B, P, i}, K_Metric);
            einsum(Indices{P, Q}, &(*k), Indices{B, i, P}, tmp, Indices{B, i, Q}, K_Oper);
        }
    }

    tensor_algebra::element([](double const &val1, double const &val2)
                            -> double { return val1 - val2; },
                            &(*f), *k);
}

void MP2F12::form_V_X(einsums::Tensor<double, 4> *V, einsums::Tensor<double, 4> *X)
{
    using namespace einsums;
    using namespace tensor_algebra;
    using namespace tensor_algebra::index;

    form_teints("FG", V, {'o', 'o', 'o', 'o'});
    form_teints("F2", X, {'o', 'o', 'o', 'o'});

    {
        Tensor<double, 4> F_oooc{"<oo|F|oC>", nocc_, nocc_, nocc_, ncabs_};
        form_teints("F", &F_oooc, {'o', 'o', 'o', 'C'});

        Tensor<double, 4> tmp{"Temp", nocc_, nocc_, nocc_, nocc_};
        {
            einsum(Indices{i, j, k, l}, &tmp, Indices{i, j, m, q}, F_oooc, Indices{k, l, m, q}, F_oooc);
            sort(1.0, Indices{i, j, k, l}, &(*X), -1.0, Indices{i, j, k, l}, tmp);
            sort(1.0, Indices{i, j, k, l}, &(*X), -1.0, Indices{j, i, l, k}, tmp);
        }

        Tensor<double, 4> G_oooc{"<oo|oC>", nocc_, nocc_, nocc_, ncabs_};
        form_teints("G", &G_oooc, {'o', 'o', 'o', 'C'});

        einsum(Indices{i, j, k, l}, &tmp, Indices{i, j, m, q}, G_oooc, Indices{k, l, m, q}, F_oooc);
        sort(1.0, Indices{i, j, k, l}, &(*V), -1.0, Indices{i, j, k, l}, tmp);
        sort(1.0, Indices{i, j, k, l}, &(*V), -1.0, Indices{j, i, l, k}, tmp);
    }

    {
        Tensor<double, 4> F_oopq{"<oo|F|OO>", nocc_, nocc_, nobs_, nobs_};
        form_teints("F", &F_oopq, {'o', 'o', 'O', 'O'});
        einsum(1.0, Indices{i, j, k, l}, &(*X), -1.0, Indices{i, j, p, q}, F_oopq, Indices{k, l, p, q}, F_oopq);

        Tensor<double, 4> G_oopq{"<oo|OO>", nocc_, nocc_, nobs_, nobs_};
        form_teints("G", &G_oopq, {'o', 'o', 'O', 'O'});
        einsum(1.0, Indices{i, j, k, l}, &(*V), -1.0, Indices{i, j, p, q}, G_oopq, Indices{k, l, p, q}, F_oopq);
    }
}

void MP2F12::form_df_V_X(einsums::Tensor<double, 4> *V, einsums::Tensor<double, 4> *X,
                         einsums::Tensor<double, 3> *J_inv_AB)
{
    using namespace einsums;
    using namespace tensor_algebra;
    using namespace tensor_algebra::index;

    form_df_teints("FG", V, J_inv_AB, {'o', 'o', 'o', 'o'});
    form_df_teints("F2", X, J_inv_AB, {'o', 'o', 'o', 'o'});

    {
        Tensor<double, 4> F_oooc{"<oo|F|oC>", nocc_, nocc_, nocc_, ncabs_};
        form_df_teints("F", &F_oooc, J_inv_AB, {'o', 'o', 'o', 'C'});

        Tensor<double, 4> tmp{"Temp", nocc_, nocc_, nocc_, nocc_};
        {
            einsum(Indices{i, j, k, l}, &tmp, Indices{i, j, m, q}, F_oooc, Indices{k, l, m, q}, F_oooc);
            sort(1.0, Indices{i, j, k, l}, &(*X), -1.0, Indices{i, j, k, l}, tmp);
            sort(1.0, Indices{i, j, k, l}, &(*X), -1.0, Indices{j, i, l, k}, tmp);
        }

        Tensor<double, 4> G_oooc{"<oo|oC>", nocc_, nocc_, nocc_, ncabs_};
        form_df_teints("G", &G_oooc, J_inv_AB, {'o', 'o', 'o', 'C'});

        einsum(Indices{i, j, k, l}, &tmp, Indices{i, j, m, q}, G_oooc, Indices{k, l, m, q}, F_oooc);
        sort(1.0, Indices{i, j, k, l}, &(*V), -1.0, Indices{i, j, k, l}, tmp);
        sort(1.0, Indices{i, j, k, l}, &(*V), -1.0, Indices{j, i, l, k}, tmp);
    }

    {
        Tensor<double, 4> F_oopq{"<oo|F|OO>", nocc_, nocc_, nobs_, nobs_};
        form_df_teints("F", &F_oopq, J_inv_AB, {'o', 'O', 'o', 'O'});
        einsum(1.0, Indices{i, j, k, l}, &(*X), -1.0, Indices{i, j, p, q}, F_oopq, Indices{k, l, p, q}, F_oopq);

        Tensor<double, 4> G_oopq{"<oo|OO>", nocc_, nocc_, nobs_, nobs_};
        form_df_teints("G", &G_oopq, J_inv_AB, {'o', 'O', 'o', 'O'});
        einsum(1.0, Indices{i, j, k, l}, &(*V), -1.0, Indices{i, j, p, q}, G_oopq, Indices{k, l, p, q}, F_oopq);
    }
}

void MP2F12::form_C(einsums::Tensor<double, 4> *C, einsums::Tensor<double, 2> *f)
{
    using namespace einsums;
    using namespace tensor_algebra;
    using namespace tensor_algebra::index;

    Tensor<double, 4> F_oovc{"<oo|F|vC>", nocc_, nocc_, nvir_, ncabs_};
    {
        Tensor<double, 4> F_oopc{"<oo|F|OC>", nocc_, nocc_, nobs_, ncabs_};
        form_teints("F", &F_oopc, {'o', 'o', 'O', 'C'});
        F_oovc = F_oopc(All, All, Range{nocc_, nobs_}, All);
    }

    Tensor f_vc = (*f)(Range{nocc_, nobs_}, Range{nobs_, nri_});
    Tensor<double, 4> tmp{"Temp", nocc_, nocc_, nvir_, nvir_};

    einsum(Indices{k, l, a, b}, &tmp, Indices{k, l, a, q}, F_oovc, Indices{b, q}, f_vc);
    sort(Indices{k, l, a, b}, &(*C), Indices{k, l, a, b}, tmp);
    sort(1.0, Indices{k, l, a, b}, &(*C), 1.0, Indices{l, k, b, a}, tmp);
}

void MP2F12::form_df_C(einsums::Tensor<double, 4> *C, einsums::Tensor<double, 2> *f,
                       einsums::Tensor<double, 3> *J_inv_AB)
{
    using namespace einsums;
    using namespace tensor_algebra;
    using namespace tensor_algebra::index;

    Tensor<double, 4> F_oovc{"<oo|F|vC>", nocc_, nocc_, nvir_, ncabs_};
    {
        Tensor<double, 4> F_oopc{"<oo|F|OC>", nocc_, nocc_, nobs_, ncabs_};
        form_df_teints("F", &F_oopc, J_inv_AB, {'o', 'O', 'o', 'C'});
        F_oovc = F_oopc(All, All, Range{nocc_, nobs_}, All);
    }

    Tensor f_vc = (*f)(Range{nocc_, nobs_}, Range{nobs_, nri_});
    Tensor<double, 4> tmp{"Temp", nocc_, nocc_, nvir_, nvir_};

    einsum(Indices{k, l, a, b}, &tmp, Indices{k, l, a, q}, F_oovc, Indices{b, q}, f_vc);
    sort(Indices{k, l, a, b}, &(*C), Indices{k, l, a, b}, tmp);
    sort(1.0, Indices{k, l, a, b}, &(*C), 1.0, Indices{l, k, b, a}, tmp);
}

void MP2F12::form_B(einsums::Tensor<double, 4> *B, einsums::Tensor<double, 2> *f,
                    einsums::Tensor<double, 2> *k)
{
    using namespace einsums;
    using namespace tensor_algebra;
    using namespace tensor_algebra::index;

    Tensor<double, 4> B_nosymm{"B_klmn", nocc_, nocc_, nocc_, nocc_};
    form_teints("Uf", &B_nosymm, {'o', 'o', 'o', 'o'});

    auto tmp_1 = std::make_unique<Tensor<double, 4>>("Temp 1", nocc_, nocc_, nocc_, nocc_);
    {
        Tensor<double, 4> F2_ooo1{"<oo|F2|o1>", nocc_, nocc_, nocc_, nri_};
        form_teints("F2", &F2_ooo1, {'o', 'o', 'o', 'O',
                                     'o', 'o', 'o', 'C'});

        Tensor<double, 2> fk_o1{"Fock-Exchange Matrix", nocc_, nri_};
        {
            auto f_o1 = (*f)(Range{0, nocc_}, All);
            auto k_o1 = (*k)(Range{0, nocc_}, All);
            tensor_algebra::element([](double const &val1, double const &val2, double const &val3)
                                    -> double { return val2 + val3; },
                                    &fk_o1, f_o1, k_o1);
        }

        einsum(Indices{l, k, n, m}, &tmp_1, Indices{l, k, n, I}, F2_ooo1, Indices{m, I}, fk_o1);
        sort(1.0, Indices{k, l, m, n}, &B_nosymm, 1.0, Indices{k, l, m, n}, *tmp_1);
        sort(1.0, Indices{k, l, m, n}, &B_nosymm, 1.0, Indices{l, k, n, m}, *tmp_1);
    }
    tmp_1.reset();

    auto F = std::make_unique<Tensor<double, 4>>("<oo|F|11>", nocc_, nocc_, nri_, nri_);
    form_teints("F", F.get(), {'o', 'o', 'O', 'O',
                               'o', 'o', 'C', 'O',
                               'o', 'o', 'C', 'C'});

    auto tmp_2 = std::make_unique<Tensor<double, 4>>("Temp 2", nocc_, nocc_, nocc_, nocc_);
    {
        tmp_1 = std::make_unique<Tensor<double, 4>>("Temp 1", nocc_, nocc_, nri_, nri_);

        einsum(Indices{l, k, P, A}, &tmp_1, Indices{l, k, P, C}, *F, Indices{C, A}, *k);
        einsum(Indices{l, k, n, m}, &tmp_2, Indices{l, k, P, A}, tmp_1, Indices{n, m, P, A}, *F);
        sort(1.0, Indices{k, l, m, n}, &B_nosymm, -1.0, Indices{k, l, m, n}, *tmp_2);
        sort(1.0, Indices{k, l, m, n}, &B_nosymm, -1.0, Indices{l, k, n, m}, *tmp_2);
    }
    tmp_1.reset();
    tmp_2.reset();

    {
        Tensor F_ooo1 = (*F)(Range{0, nocc_}, Range{0, nocc_}, Range{0, nocc_}, All);
        tmp_1 = std::make_unique<Tensor<double, 4>>("Temp 1", nocc_, nocc_, nocc_, nri_);
        tmp_2 = std::make_unique<Tensor<double, 4>>("Temp 2", nocc_, nocc_, nocc_, nocc_);

        einsum(Indices{l, k, j, A}, &tmp_1, Indices{l, k, j, C}, F_ooo1, Indices{C, A}, *f);
        einsum(Indices{l, k, n, m}, &tmp_2, Indices{l, k, j, A}, tmp_1, Indices{n, m, j, A}, F_ooo1);
        sort(1.0, Indices{k, l, m, n}, &B_nosymm, -1.0, Indices{k, l, m, n}, *tmp_2);
        sort(1.0, Indices{k, l, m, n}, &B_nosymm, -1.0, Indices{l, k, n, m}, *tmp_2);
    }
    tmp_1.reset();
    tmp_2.reset();

    TensorView<double, 4> F_ooco_temp{(*F), Dim<4>{nocc_, nocc_, ncabs_, nocc_}, Offset<4>{0, 0, nobs_, 0}};
    {
        Tensor F_ooco = F_ooco_temp;
        Tensor f_oo   = (*f)(Range{0, nocc_}, Range{0, nocc_});
        tmp_1 = std::make_unique<Tensor<double, 4>>("Temp 1", nocc_, nocc_, ncabs_, nocc_);
        tmp_2 = std::make_unique<Tensor<double, 4>>("Temp 2", nocc_, nocc_, nocc_, nocc_);

        einsum(Indices{l, k, p, i}, &tmp_1, Indices{l, k, p, j}, F_ooco, Indices{j, i}, f_oo);
        einsum(Indices{l, k, n, m}, &tmp_2, Indices{l, k, p, i}, tmp_1, Indices{n, m, p, i}, F_ooco);
        sort(1.0, Indices{k, l, m, n}, &B_nosymm, 1.0, Indices{k, l, m, n}, *tmp_2);
        sort(1.0, Indices{k, l, m, n}, &B_nosymm, 1.0, Indices{l, k, n, m}, *tmp_2);
    }
    tmp_1.reset();
    tmp_2.reset();

    TensorView<double, 4> F_oovq_temp{(*F), Dim<4>{nocc_, nocc_, nvir_, nobs_}, Offset<4>{0, 0, nocc_, 0}};
    {
        Tensor F_oovq = F_oovq_temp;
        Tensor f_pq = (*f)(Range{0, nobs_}, Range{0, nobs_});
        tmp_1 = std::make_unique<Tensor<double, 4>>("Temp 1", nocc_, nocc_, nvir_, nobs_);
        tmp_2 = std::make_unique<Tensor<double, 4>>("Temp 2", nocc_, nocc_, nocc_, nocc_);

        einsum(Indices{l, k, b, p}, &tmp_1, Indices{l, k, b, r}, F_oovq, Indices{r, p}, f_pq);
        einsum(Indices{l, k, n, m}, &tmp_2, Indices{l, k, b, p}, tmp_1, Indices{n, m, b, p}, F_oovq);
        sort(1.0, Indices{k, l, m, n}, &B_nosymm, -1.0, Indices{k, l, m, n}, *tmp_2);
        sort(1.0, Indices{k, l, m, n}, &B_nosymm, -1.0, Indices{l, k, n, m}, *tmp_2);
    }
    tmp_1.reset();
    tmp_2.reset();

    {
        Tensor F_ooco = F_ooco_temp;
        Tensor F_ooc1 = (*F)(Range{0, nocc_}, Range{0, nocc_}, Range{nobs_, nri_}, All);
        Tensor f_o1   = (*f)(Range{0, nocc_}, All);
        tmp_1 = std::make_unique<Tensor<double, 4>>("Temp 1", nocc_, nocc_, ncabs_, nocc_);
        tmp_2 = std::make_unique<Tensor<double, 4>>("Temp 2", nocc_, nocc_, nocc_, nocc_);

        einsum(Indices{l, k, p, j}, &tmp_1, Indices{l, k, p, I}, F_ooc1, Indices{j, I}, f_o1);
        einsum(0.0, Indices{l, k, n, m}, &tmp_2, 2.0, Indices{l, k, p, j}, tmp_1, Indices{n, m, p, j}, F_ooco);
        sort(1.0, Indices{k, l, m, n}, &B_nosymm, -1.0, Indices{k, l, m, n}, *tmp_2);
        sort(1.0, Indices{k, l, m, n}, &B_nosymm, -1.0, Indices{l, k, n, m}, *tmp_2);
    }
    tmp_1.reset();
    tmp_2.reset();

    {
        Tensor F_oovq = F_oovq_temp;
        Tensor F_oovc = (*F)(Range{0, nocc_}, Range{0, nocc_}, Range{nocc_, nobs_}, Range{nobs_, nri_});
        Tensor f_pc   = (*f)(Range{0, nobs_}, Range{nobs_, nri_});
        tmp_1 = std::make_unique<Tensor<double, 4>>("Temp 1", nocc_, nocc_, nvir_, ncabs_);
        tmp_2 = std::make_unique<Tensor<double, 4>>("Temp 2", nocc_, nocc_, nocc_, nocc_);

        einsum(Indices{l, k, b, q}, &tmp_1, Indices{l, k, b, r}, F_oovq, Indices{r, q}, f_pc);
        einsum(0.0, Indices{l, k, n, m}, &tmp_2, 2.0, Indices{l, k, b, q}, tmp_1, Indices{n, m, b, q}, F_oovc);
        sort(1.0, Indices{k, l, m, n}, &B_nosymm, -1.0, Indices{k, l, m, n}, *tmp_2);
        sort(1.0, Indices{k, l, m, n}, &B_nosymm, -1.0, Indices{l, k, n, m}, *tmp_2);
    }
    tmp_1.reset();
    tmp_2.reset();

    (*B) = B_nosymm(All, All, All, All);
    sort(0.5, Indices{m, n, k, l}, &(*B), 0.5, Indices{k, l, m, n}, B_nosymm);
}

void MP2F12::form_df_B(einsums::Tensor<double, 4> *B, einsums::Tensor<double, 2> *f,
                       einsums::Tensor<double, 2> *k, einsums::Tensor<double, 3> *J_inv_AB)
{
    using namespace einsums;
    using namespace tensor_algebra;
    using namespace tensor_algebra::index;

    Tensor<double, 4> B_nosymm{"B_klmn", nocc_, nocc_, nocc_, nocc_};
    form_df_teints("Uf", &B_nosymm, J_inv_AB, {'o', 'o', 'o', 'o'});

    auto tmp_1 = std::make_unique<Tensor<double, 4>>("Temp 1", nocc_, nocc_, nocc_, nocc_);
    {
        Tensor<double, 4> F2_ooo1{"<oo|F2|o1>", nocc_, nocc_, nocc_, nri_};
        form_df_teints("F2", &F2_ooo1, J_inv_AB, {'o', 'o', 'o', 'O',
                                                  'o', 'o', 'o', 'C'});

        Tensor<double, 2> fk_o1{"Fock-Exchange Matrix", nocc_, nri_};
        {
            auto f_o1 = (*f)(Range{0, nocc_}, All);
            auto k_o1 = (*k)(Range{0, nocc_}, All);
            tensor_algebra::element([](double const &val1, double const &val2, double const &val3)
                                    -> double { return val2 + val3; },
                                    &fk_o1, f_o1, k_o1);
        }

        einsum(Indices{l, k, n, m}, &tmp_1, Indices{l, k, n, I}, F2_ooo1, Indices{m, I}, fk_o1);
        sort(1.0, Indices{k, l, m, n}, &B_nosymm, 1.0, Indices{k, l, m, n}, *tmp_1);
        sort(1.0, Indices{k, l, m, n}, &B_nosymm, 1.0, Indices{l, k, n, m}, *tmp_1);
    }
    tmp_1.reset();

    auto F = std::make_unique<Tensor<double, 4>>("<oo|F|11>", nocc_, nocc_, nri_, nri_);
    form_df_teints("F", F.get(), J_inv_AB, {'o', 'O', 'o', 'O',
                                            'o', 'O', 'o', 'C',
                                            'o', 'C', 'o', 'O',
                                            'o', 'C', 'o', 'C'});

    auto tmp_2 = std::make_unique<Tensor<double, 4>>("Temp 2", nocc_, nocc_, nocc_, nocc_);
    {
        tmp_1 = std::make_unique<Tensor<double, 4>>("Temp 1", nocc_, nocc_, nri_, nri_);

        einsum(Indices{l, k, P, A}, &tmp_1, Indices{l, k, P, C}, *F, Indices{C, A}, *k);
        einsum(Indices{l, k, n, m}, &tmp_2, Indices{l, k, P, A}, tmp_1, Indices{n, m, P, A}, *F);
        sort(1.0, Indices{k, l, m, n}, &B_nosymm, -1.0, Indices{k, l, m, n}, *tmp_2);
        sort(1.0, Indices{k, l, m, n}, &B_nosymm, -1.0, Indices{l, k, n, m}, *tmp_2);
    }
    tmp_1.reset();
    tmp_2.reset();

    {
        Tensor F_ooo1 = (*F)(Range{0, nocc_}, Range{0, nocc_}, Range{0, nocc_}, All);
        tmp_1 = std::make_unique<Tensor<double, 4>>("Temp 1", nocc_, nocc_, nocc_, nri_);
        tmp_2 = std::make_unique<Tensor<double, 4>>("Temp 2", nocc_, nocc_, nocc_, nocc_);

        einsum(Indices{l, k, j, A}, &tmp_1, Indices{l, k, j, C}, F_ooo1, Indices{C, A}, *f);
        einsum(Indices{l, k, n, m}, &tmp_2, Indices{l, k, j, A}, tmp_1, Indices{n, m, j, A}, F_ooo1);
        sort(1.0, Indices{k, l, m, n}, &B_nosymm, -1.0, Indices{k, l, m, n}, *tmp_2);
        sort(1.0, Indices{k, l, m, n}, &B_nosymm, -1.0, Indices{l, k, n, m}, *tmp_2);
    }
    tmp_1.reset();
    tmp_2.reset();

    TensorView<double, 4> F_ooco_temp{(*F), Dim<4>{nocc_, nocc_, ncabs_, nocc_}, Offset<4>{0, 0, nobs_, 0}};
    {
        Tensor F_ooco = F_ooco_temp;
        Tensor f_oo   = (*f)(Range{0, nocc_}, Range{0, nocc_});
        tmp_1 = std::make_unique<Tensor<double, 4>>("Temp 1", nocc_, nocc_, ncabs_, nocc_);
        tmp_2 = std::make_unique<Tensor<double, 4>>("Temp 2", nocc_, nocc_, nocc_, nocc_);

        einsum(Indices{l, k, p, i}, &tmp_1, Indices{l, k, p, j}, F_ooco, Indices{j, i}, f_oo);
        einsum(Indices{l, k, n, m}, &tmp_2, Indices{l, k, p, i}, tmp_1, Indices{n, m, p, i}, F_ooco);
        sort(1.0, Indices{k, l, m, n}, &B_nosymm, 1.0, Indices{k, l, m, n}, *tmp_2);
        sort(1.0, Indices{k, l, m, n}, &B_nosymm, 1.0, Indices{l, k, n, m}, *tmp_2);
    }
    tmp_1.reset();
    tmp_2.reset();

    TensorView<double, 4> F_oovq_temp{(*F), Dim<4>{nocc_, nocc_, nvir_, nobs_}, Offset<4>{0, 0, nocc_, 0}};
    {
        Tensor F_oovq = F_oovq_temp;
        Tensor f_pq = (*f)(Range{0, nobs_}, Range{0, nobs_});
        tmp_1 = std::make_unique<Tensor<double, 4>>("Temp 1", nocc_, nocc_, nvir_, nobs_);
        tmp_2 = std::make_unique<Tensor<double, 4>>("Temp 2", nocc_, nocc_, nocc_, nocc_);

        einsum(Indices{l, k, b, p}, &tmp_1, Indices{l, k, b, r}, F_oovq, Indices{r, p}, f_pq);
        einsum(Indices{l, k, n, m}, &tmp_2, Indices{l, k, b, p}, tmp_1, Indices{n, m, b, p}, F_oovq);
        sort(1.0, Indices{k, l, m, n}, &B_nosymm, -1.0, Indices{k, l, m, n}, *tmp_2);
        sort(1.0, Indices{k, l, m, n}, &B_nosymm, -1.0, Indices{l, k, n, m}, *tmp_2);
    }
    tmp_1.reset();
    tmp_2.reset();

    {
        Tensor F_ooco = F_ooco_temp;
        Tensor F_ooc1 = (*F)(Range{0, nocc_}, Range{0, nocc_}, Range{nobs_, nri_}, All);
        Tensor f_o1   = (*f)(Range{0, nocc_}, All);
        tmp_1 = std::make_unique<Tensor<double, 4>>("Temp 1", nocc_, nocc_, ncabs_, nocc_);
        tmp_2 = std::make_unique<Tensor<double, 4>>("Temp 2", nocc_, nocc_, nocc_, nocc_);

        einsum(Indices{l, k, p, j}, &tmp_1, Indices{l, k, p, I}, F_ooc1, Indices{j, I}, f_o1);
        einsum(0.0, Indices{l, k, n, m}, &tmp_2, 2.0, Indices{l, k, p, j}, tmp_1, Indices{n, m, p, j}, F_ooco);
        sort(1.0, Indices{k, l, m, n}, &B_nosymm, -1.0, Indices{k, l, m, n}, *tmp_2);
        sort(1.0, Indices{k, l, m, n}, &B_nosymm, -1.0, Indices{l, k, n, m}, *tmp_2);
    }
    tmp_1.reset();
    tmp_2.reset();

    {
        Tensor F_oovq = F_oovq_temp;
        Tensor F_oovc = (*F)(Range{0, nocc_}, Range{0, nocc_}, Range{nocc_, nobs_}, Range{nobs_, nri_});
        Tensor f_pc   = (*f)(Range{0, nobs_}, Range{nobs_, nri_});
        tmp_1 = std::make_unique<Tensor<double, 4>>("Temp 1", nocc_, nocc_, nvir_, ncabs_);
        tmp_2 = std::make_unique<Tensor<double, 4>>("Temp 2", nocc_, nocc_, nocc_, nocc_);

        einsum(Indices{l, k, b, q}, &tmp_1, Indices{l, k, b, r}, F_oovq, Indices{r, q}, f_pc);
        einsum(0.0, Indices{l, k, n, m}, &tmp_2, 2.0, Indices{l, k, b, q}, tmp_1, Indices{n, m, b, q}, F_oovc);
        sort(1.0, Indices{k, l, m, n}, &B_nosymm, -1.0, Indices{k, l, m, n}, *tmp_2);
        sort(1.0, Indices{k, l, m, n}, &B_nosymm, -1.0, Indices{l, k, n, m}, *tmp_2);
    }
    tmp_1.reset();
    tmp_2.reset();

    (*B) = B_nosymm(All, All, All, All);
    sort(0.5, Indices{m, n, k, l}, &(*B), 0.5, Indices{k, l, m, n}, B_nosymm);
}

////////////////////////////////
//* Disk Algorithm (CONV/DF) *//
////////////////////////////////

void DiskMP2F12::form_fock(einsums::DiskTensor<double, 2> *f, einsums::DiskTensor<double, 2> *k,
                           einsums::DiskTensor<double, 2> *fk)
{
    using namespace einsums;
    using namespace tensor_algebra;
    using namespace tensor_algebra::index;

    form_oeints(f);
    auto f_view = (*f)(All, All);

    {
        outfile->Printf("     Forming J\n");
        auto J = DiskTensor<double, 4>{state::data, "Coulomb", nri_, nocc_, nri_, nocc_};
        if (!J.existed()) form_teints("J", &J);

        for (int i = 0; i < nocc_; i++) {
            auto J_view = J(All, i, All, i);
            sort(1.0, Indices{p, q}, &f_view.get(), 2.0, Indices{p, q}, J_view.get());
        }

        auto fk_view = (*fk)(All, All);
        sort(Indices{p, q}, &fk_view.get(), Indices{p, q}, f_view.get());
    }

    {
        outfile->Printf("     Forming K\n");
        auto K = DiskTensor<double, 4>{state::data, "Exchange", nri_, nocc_, nocc_, nri_};
        if (!K.existed()) form_teints("K", &K);

        auto k_view = (*k)(All, All);
        for (int i = 0; i < nocc_; i++) {
            auto K_view = K(All, i, i, All);
            sort(1.0, Indices{p, q}, &k_view.get(), 1.0, Indices{p, q}, K_view.get());
        }

        sort(1.0, Indices{p, q}, &f_view.get(), -1.0, Indices{p, q}, k_view.get());
    }
}

void DiskMP2F12::form_df_fock(einsums::DiskTensor<double, 2> *f, einsums::DiskTensor<double, 2> *k,
                              einsums::DiskTensor<double, 2> *fk)
{
    using namespace einsums;
    using namespace tensor_algebra;
    using namespace tensor_algebra::index;

    form_oeints(f);
    auto f_view = (*f)(All, All);

    {
        auto Metric = std::make_unique<Tensor<double, 3>>("(B|PQ) MO", naux_, nri_, nri_);
        form_metric_ints(Metric.get(), true);
        auto Oper = std::make_unique<Tensor<double, 3>>("(B|PQ) MO", naux_, nocc_, nri_);
        form_oper_ints("G", Oper.get());

        {
            outfile->Printf("     Forming J\n");
            Tensor<double, 1> tmp{"B", naux_};
            {
                Tensor Id = create_identity_tensor("I", nocc_, nocc_);
                Tensor J_Oper = (*Oper)(Range{0, naux_}, Range{0, nocc_}, Range{0, nocc_});
                einsum(Indices{B}, &tmp, Indices{B, i, j}, J_Oper, Indices{i, j}, Id);
            }

            Tensor J_Metric = (*Metric)(Range{0, naux_}, Range{0, nri_}, Range{0, nri_});
            einsum(1.0, Indices{P, Q}, &f_view.get(), 2.0, Indices{B, P, Q}, J_Metric, Indices{B}, tmp);

            auto fk_view = (*fk)(All, All);
            sort(Indices{p, q}, &fk_view.get(), Indices{p, q}, f_view.get());
        }

        {
            outfile->Printf("     Forming K\n");
            Tensor<double, 3> tmp{"", naux_, nocc_, nri_};
            {
                Tensor K_Metric = (*Metric)(Range{0, naux_}, Range{0, nri_}, Range{0, nocc_});
                sort(Indices{B, i, P}, &tmp, Indices{B, P, i}, K_Metric);
            }

            auto k_view = (*k)(All, All);
            Tensor K_Oper = (*Oper)(Range{0, naux_}, Range{0, nocc_}, Range{0, nri_});
            einsum(Indices{P, Q}, &k_view.get(), Indices{B, i, P}, tmp, Indices{B, i, Q}, K_Oper);

            sort(1.0, Indices{p, q}, &f_view.get(), -1.0, Indices{p, q}, k_view.get());
        }
    }
}

void DiskMP2F12::form_V_X(einsums::DiskTensor<double, 4> *VX, einsums::DiskTensor<double, 4> *F,
                          einsums::DiskTensor<double, 4> *G_F, einsums::DiskTensor<double, 4> *FG_F2)
{
    using namespace einsums;
    using namespace tensor_algebra;
    using namespace tensor_algebra::index;

    Tensor<double, 0> tmp1{"F_oooc . G_F_oooc IJKL"};
    Tensor<double, 0> tmp2{"F_oooc . G_F_oooc JILK"};
    Tensor<double, 0> tmp3{"F_oopq . G_F_oopq IJKL"};

    for (int I = 0; I < nocc_; I++) {
        for (int J = 0; J < nocc_; J++) {
            auto VX_IJ = (*VX)(All, All, I, J);
            auto FG_F2_IJ = (*FG_F2)(I, J, All, All);
            auto F_IJ_oc = (*F)(I, J, Range{0, nocc_}, Range{nobs_, nri_});
            auto F_JI_oc = (*F)(J, I, Range{0, nocc_}, Range{nobs_, nri_});
            auto F_pq = (*F)(I, J, Range{0, nobs_}, Range{0, nobs_});

            for(int K = 0; K < nocc_; K++) {
                for (int L = 0; L < nocc_; L++) {
                    auto G_F_KL_oc = (*G_F)(K, L, Range{0, nocc_}, Range{nobs_, nri_});
                    auto G_F_LK_oc = (*G_F)(L, K, Range{0, nocc_}, Range{nobs_, nri_});
                    einsum(Indices{}, &tmp1, Indices{m, q}, F_IJ_oc.get(), Indices{m, q}, G_F_KL_oc.get());
                    einsum(Indices{}, &tmp2, Indices{m, q}, G_F_LK_oc.get(), Indices{m, q}, F_JI_oc.get());

                    auto G_F_pq = (*G_F)(K, L, Range{0, nobs_}, Range{0, nobs_});
                    einsum(Indices{}, &tmp3, Indices{p, q}, F_pq.get(), Indices{p, q}, G_F_pq.get());

                    VX_IJ(K, L) = FG_F2_IJ(K, L) - tmp1 - tmp2 - tmp3;
                }
            }
        }
    }
}

void DiskMP2F12::form_C(einsums::DiskTensor<double, 4> *C, einsums::DiskTensor<double, 4> *F,
                    einsums::DiskTensor<double, 2> *f)
{
    using namespace einsums;
    using namespace tensor_algebra;
    using namespace tensor_algebra::index;

    auto f_vc = (*f)(Range{nocc_, nobs_}, Range{nobs_, nri_});
    for(int K = 0; K < nocc_; K++) {
        for (int L = 0; L < nocc_; L++) {
            auto C_KL = (*C)(K, L, All, All);
            auto F_KL_vc = (*F)(K, L, Range{nocc_, nobs_}, Range{nobs_, nri_});
            auto F_LK_vc = (*F)(L, K, Range{nocc_, nobs_}, Range{nobs_, nri_});

            einsum(Indices{a, b}, &C_KL.get(), Indices{a, q}, F_KL_vc.get(), Indices{b, q}, f_vc.get());
            einsum(1.0, Indices{a, b}, &C_KL.get(), 1.0, Indices{a, q}, f_vc.get(), Indices{b, q}, F_LK_vc.get());
        }
    }
}

void DiskMP2F12::form_B(einsums::DiskTensor<double, 4> *B, einsums::DiskTensor<double, 4> *Uf,
                    einsums::DiskTensor<double, 4> *F2, einsums::DiskTensor<double, 4> *F,
                    einsums::DiskTensor<double, 2> *f, einsums::DiskTensor<double, 2> *fk,
                    einsums::DiskTensor<double, 2> *kk)
{
    using namespace einsums;
    using namespace tensor_algebra;
    using namespace tensor_algebra::index;

    DiskTensor<double, 4> B_nosymm{state::data, "B_klmn", nocc_, nocc_, nocc_, nocc_};

    // Term 1 and Term 2
    {
        auto fk_o1   = (*fk)(Range{0, nocc_}, All);
        for(int K = 0; K < nocc_; K++) {
            for (int L = 0; L < nocc_; L++) {
                auto B_KL = B_nosymm(K, L, All, All);

                auto Uf_KL = (*Uf)(K, L, All, All);
                sort(Indices{m, n}, &B_KL.get(), Indices{m, n}, Uf_KL.get());

                auto F2_KL_o1 = (*F2)(K, L, Range{0, nocc_}, All);
                auto F2_LK_o1 = (*F2)(L, K, Range{0, nocc_}, All);
                einsum(1.0, Indices{m, n}, &B_KL.get(), 1.0, Indices{m, I}, fk_o1.get(), Indices{n, I}, F2_LK_o1.get());
                einsum(1.0, Indices{m, n}, &B_KL.get(), 1.0, Indices{m, I}, F2_KL_o1.get(), Indices{n, I}, fk_o1.get());
            }
        }
    }

    // Term 3
    {
        Tensor<double, 2> rank2{"Contraction 1", nri_, nri_};
        Tensor<double, 0> tmp{"Contraction 2"};
        auto k_view = (*kk)(All, All);

        for(int K = 0; K < nocc_; K++) {
            for (int L = 0; L < nocc_; L++) {
                auto B_KL = B_nosymm(K, L, All, All);
                auto F_KL_11 = (*F)(K, L, All, All);
                auto F_LK_11 = (*F)(L, K, All, All);

                for(int M = 0; M < nocc_; M++) {
                    for (int N = 0; N < nocc_; N++) {
                        auto F_MN_11 = (*F)(M, N, All, All);
                        einsum(Indices{P, A}, &rank2, Indices{P, C}, F_KL_11.get(), Indices{C, A}, k_view.get());
                        einsum(Indices{}, &tmp, Indices{P, A}, rank2, Indices{P, A}, F_MN_11.get());
                        B_KL(M, N) -= tmp;

                        auto F_NM_11 = (*F)(N, M, All, All);
                        einsum(Indices{P, A}, &rank2, Indices{P, C}, F_LK_11.get(), Indices{C, A}, k_view.get());
                        einsum(Indices{}, &tmp, Indices{P, A}, rank2, Indices{P, A}, F_NM_11.get());
                        B_KL(M, N) -= tmp;
                    }
                }
            }
        }
    }

    // Term 4
    {
        Tensor<double, 2> rank2{"Contraction 1", nocc_, nri_};
        Tensor<double, 0> tmp{"Contraction 2"};
        auto f_view = (*f)(All, All);

        for(int K = 0; K < nocc_; K++) {
            for (int L = 0; L < nocc_; L++) {
                auto B_KL = B_nosymm(K, L, All, All);
                auto F_KL_o1 = (*F)(K, L, Range{0, nocc_}, All);
                auto F_LK_o1 = (*F)(L, K, Range{0, nocc_}, All);

                for(int M = 0; M < nocc_; M++) {
                    for (int N = 0; N < nocc_; N++) {
                        auto F_MN_o1 = (*F)(M, N, Range{0, nocc_}, All);
                        einsum(Indices{j, A}, &rank2, Indices{j, C}, F_KL_o1.get(), Indices{C, A}, f_view.get());
                        einsum(Indices{}, &tmp, Indices{j, A}, rank2, Indices{j, A}, F_MN_o1.get());
                        B_KL(M, N) -= tmp;

                        auto F_NM_o1 = (*F)(N, M, Range{0, nocc_}, All);
                        einsum(Indices{j, A}, &rank2, Indices{j, C}, F_LK_o1.get(), Indices{C, A}, f_view.get());
                        einsum(Indices{}, &tmp, Indices{j, A}, rank2, Indices{j, A}, F_NM_o1.get());
                        B_KL(M, N) -= tmp;
                    }
                }
            }
        }
    }

    // Term 5 and Term 7
    {
        Tensor<double, 2> rank2{"Contraction 1", ncabs_, nocc_};
        Tensor<double, 0> tmp{"Contraction 2"};
        auto f_oo = (*f)(Range{0, nocc_}, Range{0, nocc_});
        auto f_o1 = (*f)(Range{0, nocc_}, All);

        for(int K = 0; K < nocc_; K++) {
            for (int L = 0; L < nocc_; L++) {
                auto B_KL = B_nosymm(K, L, All, All);
                auto F_KL_co = (*F)(K, L, Range{nobs_, nri_}, Range{0, nocc_});
                auto F_LK_co = (*F)(L, K, Range{nobs_, nri_}, Range{0, nocc_});
                auto F_KL_c1 = (*F)(K, L, Range{nobs_, nri_}, All);
                auto F_LK_c1 = (*F)(L, K, Range{nobs_, nri_}, All);

                for(int M = 0; M < nocc_; M++) {
                    for (int N = 0; N < nocc_; N++) {
                        {
                            auto F_MN_co = (*F)(M, N, Range{nobs_, nri_}, Range{0, nocc_});
                            einsum(Indices{p, i}, &rank2, Indices{p, j}, F_KL_co.get(), Indices{i, j}, f_oo.get());
                            einsum(Indices{}, &tmp, Indices{p, i}, rank2, Indices{p, i}, F_MN_co.get());
                            B_KL(M, N) += tmp;

                            einsum(Indices{p, j}, &rank2, Indices{p, I}, F_KL_c1.get(), Indices{j, I}, f_o1.get());
                            einsum(Indices{}, &tmp, Indices{p, j}, rank2, Indices{p, j}, F_MN_co.get());
                            B_KL(M, N) -= 2.0 * tmp;
                        }

                        {
                            auto F_NM_co = (*F)(N, M, Range{nobs_, nri_}, Range{0, nocc_});
                            einsum(Indices{p, i}, &rank2, Indices{p, j}, F_LK_co.get(), Indices{i, j}, f_oo.get());
                            einsum(Indices{}, &tmp, Indices{p, i}, rank2, Indices{p, i}, F_NM_co.get());
                            B_KL(M, N) += tmp;

                            einsum(Indices{p, j}, &rank2, Indices{p, I}, F_LK_c1.get(), Indices{j, I}, f_o1.get());
                            einsum(Indices{}, &tmp, Indices{p, j}, rank2, Indices{p, j}, F_NM_co.get());
                            B_KL(M, N) -= 2.0 * tmp;
                        }
                    }
                }
            }
        }
    }

    // Term 6 and Term 8
    {
        Tensor<double, 2> rank2{"Contraction 1", nvir_, nobs_};
        Tensor<double, 0> tmp{"Contraction 2"};
        auto f_pq = (*f)(Range{0, nobs_}, Range{0, nobs_});
        auto f_pc   = (*f)(Range{0, nobs_}, Range{nobs_, nri_});

        for(int K = 0; K < nocc_; K++) {
            for (int L = 0; L < nocc_; L++) {
                auto B_KL = B_nosymm(K, L, All, All);
                auto F_KL_vq = (*F)(K, L, Range{nocc_, nobs_}, Range{0, nobs_});
                auto F_LK_vq = (*F)(L, K, Range{nocc_, nobs_}, Range{0, nobs_});
                auto F_KL_vc = (*F)(K, L, Range{nocc_, nobs_}, Range{nobs_, nri_});
                auto F_LK_vc = (*F)(L, K, Range{nocc_, nobs_}, Range{nobs_, nri_});

                for(int M = 0; M < nocc_; M++) {
                    for (int N = 0; N < nocc_; N++) {
                        {
                            auto F_MN_vq = (*F)(M, N, Range{nocc_, nobs_}, Range{0, nobs_});
                            einsum(Indices{b, q}, &rank2, Indices{b, p}, F_KL_vq.get(), Indices{p, q}, f_pq.get());
                            einsum(Indices{}, &tmp, Indices{b, q}, rank2, Indices{b, q}, F_MN_vq.get());
                            B_KL(M, N) -= tmp;

                            einsum(Indices{b, q}, &rank2, Indices{b, w}, F_KL_vc.get(), Indices{q, w}, f_pc.get());
                            einsum(Indices{}, &tmp, Indices{b, q}, rank2, Indices{b, q}, F_MN_vq.get());
                            B_KL(M, N) -= 2.0 * tmp;
                        }

                        {
                            auto F_NM_vq = (*F)(N, M, Range{nocc_, nobs_}, Range{0, nobs_});
                            einsum(Indices{b, q}, &rank2, Indices{b, p}, F_LK_vq.get(), Indices{p, q}, f_pq.get());
                            einsum(Indices{}, &tmp, Indices{b, q}, rank2, Indices{b, q}, F_NM_vq.get());
                            B_KL(M, N) -= tmp;

                            einsum(Indices{b, q}, &rank2, Indices{b, w}, F_LK_vc.get(), Indices{q, w}, f_pc.get());
                            einsum(Indices{}, &tmp, Indices{b, q}, rank2, Indices{b, q}, F_NM_vq.get());
                            B_KL(M, N) -= 2.0 * tmp;
                        }
                    }
                }
            }
        }
    }

    for(int K = 0; K < nocc_; K++) {
        for (int L = 0; L < nocc_; L++) {
            auto B_KL = (*B)(K, L, All, All);
            auto B_klmn = B_nosymm(K, L, All, All);
            auto B_mnkl = B_nosymm(All, All, K, L);
            for(int M = 0; M < nocc_; M++) {
                for (int N = 0; N < nocc_; N++) {
                    B_KL(M, N) = 0.5 * (B_klmn(M, N) + B_mnkl(M, N));
                }
            }
        }
    }
}

}} // end namespaces
