/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2018 The Psi4 Developers.
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

#include "psi4/libmoinfo/libmoinfo.h"
#include "psi4/libpsi4util/libpsi4util.h"
#include <cstdio>

#include "blas.h"
#include "debugging.h"
#include "matrix.h"

namespace psi {

namespace psimrcc {
extern MOInfo* moinfo;

extern MemoryManager* memory_manager;

void CCBLAS::zero(const char* cstr) {
    std::string str(cstr);
    // To zero diagonals of things like "Fae[v][v]{u}"
    std::vector<std::string> names = moinfo->get_matrix_names(str);
    for (size_t n = 0; n < names.size(); ++n) {
        CCMatrix* Matrix = get_Matrix(names[n]);
        Matrix->zero_matrix();
        DEBUGGING(5, outfile->Printf("\n...setting %s to zero", names[n].c_str()););
    }
}

void CCBLAS::zero_right_four_diagonal(const char* cstr) {
    std::string str(cstr);
    // To zero diagonals of things like "Fae[v][v]{u}"
    std::vector<std::string> names = moinfo->get_matrix_names(str);
    for (size_t n = 0; n < names.size(); ++n) {
        CCMatrix* Matrix = get_Matrix(names[n]);
        Matrix->zero_right_four_diagonal();
        DEBUGGING(5, outfile->Printf("\n...setting the right diagonal terms of %s to zero", names[n].c_str()););
    }
}

void CCBLAS::zero_non_doubly_occupied(const char* cstr) {
    std::string str(cstr);
    // To zero non-doubly occupied MOs of things like "Fae[v][v]{u}"
    std::vector<std::string> names = moinfo->get_matrix_names(str);
    for (size_t n = 0; n < names.size(); ++n) {
        CCMatrix* Matrix = get_Matrix(names[n]);
        Matrix->zero_non_doubly_occupied();
        DEBUGGING(5, outfile->Printf("\n...setting the right diagonal terms of %s to zero", names[n].c_str()););
    }
}

void CCBLAS::zero_non_external(const char* cstr) {
    std::string str(cstr);
    // To zero non-external MOs of things like "Fae[v][v]{u}"
    std::vector<std::string> names = moinfo->get_matrix_names(str);
    for (size_t n = 0; n < names.size(); ++n) {
        CCMatrix* Matrix = get_Matrix(names[n]);
        Matrix->zero_non_external();
        DEBUGGING(5, outfile->Printf("\n...setting the right diagonal terms of %s to zero", names[n].c_str()););
    }
}

void CCBLAS::scale(const char* cstr, int reference, double value) {
    std::string str(cstr);
    scale(str, reference, value);
}

void CCBLAS::scale(std::string& str, int reference, double value) {
    std::string matrix_str = add_reference(str, reference);
    // Make sure that the element that we are retrieving is present
    MatrixMap::iterator iter = matrices.find(matrix_str);
    if (iter != matrices.end()) {
        load(iter->second);
        iter->second->scale(value);
        return;
    }
    throw PSIEXCEPTION("\nCCBLAS::scale() couldn't find matrix " + matrix_str);
}

void CCBLAS::reduce_spaces(const char* out, const char* in) {
    std::string in_str(in);
    std::string out_str(out);
    // To zero diagonals of things like "Fae[v][v]{u}"
    std::vector<std::string> in_names = moinfo->get_matrix_names(in_str);
    std::vector<std::string> out_names = moinfo->get_matrix_names(out_str);
    if (in_names.size() != out_names.size()) throw PSIEXCEPTION("CCBLAS::map_spaces, number of references mismatch");
    for (size_t n = 0; n < in_names.size(); ++n) {
        CCMatrix* in_Matrix = get_Matrix(in_names[n]);
        CCMatrix* out_Matrix = get_Matrix(out_names[n]);
        process_reduce_spaces(out_Matrix, in_Matrix);
    }
}

void CCBLAS::process_reduce_spaces(CCMatrix* out_Matrix, CCMatrix* in_Matrix) {
    double*** out_matrix = out_Matrix->get_matrix();
    const intvec& act_to_occ = moinfo->get_actv_to_occ();
    const intvec& act_to_vir = moinfo->get_actv_to_vir();

    std::string& out_index_label = out_Matrix->get_index_label();
    std::string& in_index_label = in_Matrix->get_index_label();

    int index_label_size = out_index_label.size();

    int** map;
    allocate2(int, map, index_label_size, moinfo->get_nmo());

    for (int k = 0; k < index_label_size; k++) {
        if (out_index_label[k] == 'a' && in_index_label[k] == 'o') {
            for (int l = 0; l < moinfo->get_nactv(); l++) {
                map[k][l] = act_to_occ[l];
            }
        } else if (out_index_label[k] == 'a' && in_index_label[k] == 'v') {
            for (int l = 0; l < moinfo->get_nactv(); l++) {
                map[k][l] = act_to_vir[l];
            }
        } else {
            for (int l = 0; l < moinfo->get_nmo(); l++) {
                map[k][l] = l;
            }
        }
    }

    if (index_label_size == 2) {
        auto* pq = new short[2];
        for (int h = 0; h < moinfo->get_nirreps(); h++) {
            for (size_t i = 0; i < out_Matrix->get_left_pairpi(h); ++i) {
                for (size_t j = 0; j < out_Matrix->get_right_pairpi(h); ++j) {
                    out_Matrix->get_two_indices(pq, h, i, j);
                    out_matrix[h][i][j] = in_Matrix->get_two_address_element(map[0][pq[0]], map[1][pq[1]]);
                }
            }
        }
        delete[] pq;
    } else if (index_label_size == 4) {
        auto* pqrs = new short[4];
        for (int h = 0; h < moinfo->get_nirreps(); h++) {
            for (size_t i = 0; i < out_Matrix->get_left_pairpi(h); ++i) {
                for (size_t j = 0; j < out_Matrix->get_right_pairpi(h); ++j) {
                    out_Matrix->get_four_indices(pqrs, h, i, j);
                    out_matrix[h][i][j] = in_Matrix->get_four_address_element(map[0][pqrs[0]], map[1][pqrs[1]],
                                                                              map[2][pqrs[2]], map[3][pqrs[3]]);
                }
            }
        }
        delete[] pqrs;
    }
    release2(map);
}

void CCBLAS::expand_spaces(const char* out, const char* in) {
    std::string in_str(in);
    std::string out_str(out);

    std::vector<std::string> in_names = moinfo->get_matrix_names(in_str);
    std::vector<std::string> out_names = moinfo->get_matrix_names(out_str);
    if (in_names.size() != out_names.size()) throw PSIEXCEPTION("CCBLAS::map_spaces, number of references mismatch");
    for (size_t n = 0; n < in_names.size(); ++n) {
        CCMatrix* in_Matrix = get_Matrix(in_names[n]);
        CCMatrix* out_Matrix = get_Matrix(out_names[n]);
        process_expand_spaces(out_Matrix, in_Matrix);
    }
}

void CCBLAS::process_expand_spaces(CCMatrix* out_Matrix, CCMatrix* in_Matrix) {
    double*** out_matrix = out_Matrix->get_matrix();
    const intvec& act_to_occ = moinfo->get_actv_to_occ();
    const intvec& act_to_vir = moinfo->get_actv_to_vir();

    std::string& out_index_label = out_Matrix->get_index_label();
    std::string& in_index_label = in_Matrix->get_index_label();

    int index_label_size = out_index_label.size();

    int** map;
    allocate2(int, map, index_label_size, moinfo->get_nmo());

    for (int k = 0; k < index_label_size; k++) {
        if (out_index_label[k] == 'a' && in_index_label[k] == 'o') {
            for (int l = 0; l < moinfo->get_nactv(); l++) {
                map[k][l] = act_to_occ[l];
            }
        } else if (out_index_label[k] == 'a' && in_index_label[k] == 'v') {
            for (int l = 0; l < moinfo->get_nactv(); l++) {
                map[k][l] = act_to_vir[l];
            }
        } else {
            for (int l = 0; l < moinfo->get_nmo(); l++) {
                map[k][l] = l;
            }
        }
    }

    if (index_label_size == 2) {
        auto* pq = new short[2];

        for (int h = 0; h < moinfo->get_nirreps(); h++) {
            for (size_t i = 0; i < out_Matrix->get_left_pairpi(h); ++i) {
                for (size_t j = 0; j < out_Matrix->get_right_pairpi(h); ++j) {
                    out_Matrix->get_two_indices(pq, h, i, j);
                    in_Matrix->set_two_address_element(map[0][pq[0]], map[1][pq[1]], out_matrix[h][i][j]);
                }
            }
        }

        delete[] pq;
    } else if (index_label_size == 4) {
        auto* pqrs = new short[4];
        for (int h = 0; h < moinfo->get_nirreps(); h++) {
            for (size_t i = 0; i < out_Matrix->get_left_pairpi(h); ++i) {
                for (size_t j = 0; j < out_Matrix->get_right_pairpi(h); ++j) {
                    out_Matrix->get_four_indices(pqrs, h, i, j);
                    in_Matrix->set_four_address_element(map[0][pqrs[0]], map[1][pqrs[1]], map[2][pqrs[2]],
                                                        map[3][pqrs[3]], out_matrix[h][i][j]);
                }
            }
        }
        delete[] pqrs;
    }
    release2(map);
}

}  // namespace psimrcc
}  // namespace psi
