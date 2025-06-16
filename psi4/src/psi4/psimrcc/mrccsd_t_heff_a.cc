/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2025 The Psi4 Developers.
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

#include "blas.h"
#include "index_iterator.h"
#include "mrccsd_t.h"
#include "special_matrices.h"

namespace psi {
namespace psimrcc {

double MRCCSD_T::compute_A_ooo_contribution_to_Heff(int u_abs, int x_abs, int i_abs, int j_abs, int k_abs, int mu,
                                                    BlockMatrix* T3) {
    double value = 0.0;
    int i_sym = o->get_tuple_irrep(i_abs);
    int j_sym = o->get_tuple_irrep(j_abs);
    int k_sym = o->get_tuple_irrep(k_abs);
    int ijk_sym = i_sym ^ j_sym ^ k_sym;

    int x_sym = v->get_tuple_irrep(x_abs);
    int jk_sym = oo->get_tuple_irrep(j_abs, k_abs);

    size_t x_rel = v->get_tuple_rel_index(x_abs);
    size_t jk_rel = oo->get_tuple_rel_index(j_abs, k_abs);

    if (i_abs == u_abs) {
        CCIndexIterator ef(wfn_, "[vv]", ijk_sym ^ x_sym);
        for (ef.first(); !ef.end(); ef.next()) {
            int ef_sym = vv->get_tuple_irrep(ef.ind_abs<0>(), ef.ind_abs<1>());
            size_t ef_rel = vv->get_tuple_rel_index(ef.ind_abs<0>(), ef.ind_abs<1>());
            if (jk_sym == ef_sym) {
                value += 0.25 * T3->get(x_sym, x_rel, ef_rel) * V_oovv[jk_sym][jk_rel][ef_rel];
            }
        }
    }
    return value;
}

double MRCCSD_T::compute_A_ooO_contribution_to_Heff(int u_abs, int x_abs, int i_abs, int j_abs, int k_abs, int mu,
                                                    BlockMatrix* T3) {
    double value = 0.0;
    int i_sym = o->get_tuple_irrep(i_abs);
    int j_sym = o->get_tuple_irrep(j_abs);
    int k_sym = o->get_tuple_irrep(k_abs);
    int ijk_sym = i_sym ^ j_sym ^ k_sym;

    int x_sym = v->get_tuple_irrep(x_abs);
    int jk_sym = oo->get_tuple_irrep(j_abs, k_abs);

    size_t x_rel = v->get_tuple_rel_index(x_abs);
    size_t jk_rel = oo->get_tuple_rel_index(j_abs, k_abs);

    if (i_abs == u_abs) {
        CCIndexIterator ef(wfn_, "[vv]", ijk_sym ^ x_sym);
        for (ef.first(); !ef.end(); ef.next()) {
            int ef_sym = vv->get_tuple_irrep(ef.ind_abs<0>(), ef.ind_abs<1>());
            size_t ef_rel = vv->get_tuple_rel_index(ef.ind_abs<0>(), ef.ind_abs<1>());
            if (jk_sym == ef_sym) {
                value += T3->get(x_sym, x_rel, ef_rel) * V_oOvV[jk_sym][jk_rel][ef_rel];
            }
        }
    }
    return value;
}

double MRCCSD_T::compute_A_oOO_contribution_to_Heff(int u_abs, int x_abs, int i_abs, int j_abs, int k_abs, int mu,
                                                    BlockMatrix* T3) {
    double value = 0.0;
    int i_sym = o->get_tuple_irrep(i_abs);
    int j_sym = o->get_tuple_irrep(j_abs);
    int k_sym = o->get_tuple_irrep(k_abs);
    int ijk_sym = i_sym ^ j_sym ^ k_sym;

    int x_sym = v->get_tuple_irrep(x_abs);
    int jk_sym = oo->get_tuple_irrep(j_abs, k_abs);

    size_t x_rel = v->get_tuple_rel_index(x_abs);
    size_t jk_rel = oo->get_tuple_rel_index(j_abs, k_abs);

    if (i_abs == u_abs) {
        CCIndexIterator ef(wfn_, "[vv]", ijk_sym ^ x_sym);
        for (ef.first(); !ef.end(); ef.next()) {
            int ef_sym = vv->get_tuple_irrep(ef.ind_abs<0>(), ef.ind_abs<1>());
            size_t ef_rel = vv->get_tuple_rel_index(ef.ind_abs<0>(), ef.ind_abs<1>());
            if (jk_sym == ef_sym) {
                value += 0.25 * T3->get(x_sym, x_rel, ef_rel) * V_oovv[jk_sym][jk_rel][ef_rel];
            }
        }
    }
    return value;
}

}  // namespace psimrcc
}  // namespace psi
