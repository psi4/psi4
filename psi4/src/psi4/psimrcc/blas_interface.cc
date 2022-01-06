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

#include <cstdlib>

#include "psi4/libmoinfo/libmoinfo.h"
#include "psi4/libciomr/libciomr.h"
#include "psi4/libqt/qt.h"
#include "psi4/libpsi4util/libpsi4util.h"

#include "blas.h"
#include "index.h"
#include "matrix.h"

namespace psi {

namespace psimrcc {

void CCBLAS::add_index(const char* cstr) {
    // Make sure that the element that we are adding is not present
    std::string str(cstr);
    to_lower(str);
    if (indices.find(str) == indices.end()) {
        indices.insert(make_pair(str, new CCIndex(wfn_, str)));
    }
}

void CCBLAS::add_Matrix(const char* cstr) {
    std::string str(cstr);
    std::vector<std::string> names = wfn_->moinfo()->get_matrix_names(str);
    for (size_t n = 0; n < names.size(); ++n) add_Matrix_ref(names[n]);
}

void CCBLAS::add_Matrix(std::string str) {
    std::vector<std::string> names = wfn_->moinfo()->get_matrix_names(str);
    for (size_t n = 0; n < names.size(); ++n) add_Matrix_ref(names[n]);
}

void CCBLAS::add_Matrix_ref(std::string& str) {
    // Make sure that the element that we are adding is not present
    MatrixMap::iterator iter = matrices.find(str);
    if (iter == matrices.end()) {
        CCIndex* index_pointer[2];
        // Default: assume the [] indexing
        index_pointer[0] = get_index("[]");
        index_pointer[1] = get_index("[]");
        std::vector<std::string> index_string_vec = split_indices(str);
        for (size_t i = 0; i < index_string_vec.size(); ++i) index_pointer[i] = get_index(index_string_vec[i]);
        matrices.insert(make_pair(str, new CCMatrix(str, index_pointer[0], index_pointer[1])));
    }
}

CCIndex* CCBLAS::get_index(const char* cstr) {
    std::string str(cstr);
    to_lower(str);
    // Make sure that the element that we are retrieving is present
    IndexMap::iterator iter = indices.find(str);
    if (iter != indices.end()) {
        return (indices[str]);
    }
    throw PSIEXCEPTION("\nCCBLAS::get_index() couldn't find index " + str);
    return (nullptr);
}

CCIndex* CCBLAS::get_index(std::string& str) {
    to_lower(str);
    // Make sure that the element that we are retrieving is present
    IndexMap::iterator iter = indices.find(str);
    if (iter != indices.end()) {
        return (indices[str]);
    }
    throw PSIEXCEPTION("\nCCBLAS::get_index() couldn't find index " + str);
    return (nullptr);
}

CCMatTmp CCBLAS::get_MatTmp(std::string str, int reference, DiskOpt disk_option) {
    append_reference(str, reference);
    load(get_Matrix(str));
    return (CCMatTmp(get_Matrix(str), disk_option));
}

CCMatTmp CCBLAS::get_MatTmp(std::string str, DiskOpt disk_option) {
    load(get_Matrix(str));
    return (CCMatTmp(get_Matrix(str), disk_option));
}

CCMatTmp CCBLAS::get_MatTmp(CCMatrix* Matrix, DiskOpt disk_option) {
    load(Matrix);
    return (CCMatTmp(Matrix, disk_option));
}

CCMatIrTmp CCBLAS::get_MatIrTmp(std::string str, int reference, int irrep, DiskOpt disk_option) {
    append_reference(str, reference);
    load_irrep(get_Matrix(str), irrep);
    return (CCMatIrTmp(get_Matrix(str), irrep, disk_option));
}

CCMatIrTmp CCBLAS::get_MatIrTmp(std::string str, int irrep, DiskOpt disk_option) {
    load_irrep(get_Matrix(str), irrep);
    return (CCMatIrTmp(get_Matrix(str), irrep, disk_option));
}

CCMatIrTmp CCBLAS::get_MatIrTmp(CCMatrix* Matrix, int irrep, DiskOpt disk_option) {
    load_irrep(Matrix, irrep);
    return (CCMatIrTmp(Matrix, irrep, disk_option));
}

CCMatrix* CCBLAS::get_Matrix(const char* cstr, int reference) {
    std::string str(cstr);
    append_reference(str, reference);
    return (get_Matrix(str));
}

CCMatrix* CCBLAS::get_Matrix(const char* cstr) {
    std::string str(cstr);
    return (get_Matrix(str));
}

CCMatrix* CCBLAS::get_Matrix(std::string& str) {
    // Make sure that the element that we are retrieving is present
    MatrixMap::iterator iter = matrices.find(str);
    if (iter != matrices.end()) return (matrices[str]);
    throw PSIEXCEPTION("\nCCBLAS::get_matrix() couldn't find matrix " + str);
    return (nullptr);
}

CCMatrix* CCBLAS::get_Matrix(std::string& str, std::string& expression) {
    // Make sure that the element that we are retrieving is present
    MatrixMap::iterator iter = matrices.find(str);
    if (iter != matrices.end()) {
        return (matrices[str]);
    }
    throw PSIEXCEPTION("\n\nCCBLAS::parse() couldn't find the matrix " + str +
                       " in the CCMatrix list\n\nwhile parsing the string:\n\t " + expression + "\n\n");
    return nullptr;
}

void CCBLAS::set_scalar(const char* cstr, int reference, double value) {
    std::string str(cstr);
    set_scalar(str, reference, value);
}

void CCBLAS::set_scalar(std::string& str, int reference, double value) {
    std::string matrix_str = add_reference(str, reference);
    // Make sure that the element that we are retrieving is present
    MatrixMap::iterator iter = matrices.find(matrix_str);
    if (iter != matrices.end()) {
        load(iter->second);
        iter->second->set_scalar(value);
        return;
    }
    throw PSIEXCEPTION("\nCCBLAS::set_scalar() couldn't find matrix " + matrix_str);
}

double CCBLAS::get_scalar(const char* cstr, int reference) {
    std::string str(cstr);
    return (get_scalar(str, reference));
}

double CCBLAS::get_scalar(std::string& str, int reference) {
    std::string matrix_str(str);
    append_reference(matrix_str, reference);
    // Make sure that the element that we are retrieving is present
    MatrixMap::iterator iter = matrices.find(matrix_str);
    if (iter != matrices.end()) {
        load(iter->second);
        return (iter->second->get_scalar());
    }
    throw PSIEXCEPTION("\nCCBLAS::get_scalar() couldn't find matrix " + matrix_str);
    return (0.0);
}

double CCBLAS::get_scalar(std::string str) {
    // Make sure that the element that we are retrieving is present
    MatrixMap::iterator iter = matrices.find(str);
    if (iter != matrices.end()) {
        load(iter->second);
        return (iter->second->get_scalar());
    }
    throw PSIEXCEPTION("\nCCBLAS::get_scalar() couldn't find matrix " + str);
    return (0.0);
}

void CCBLAS::load(CCMatrix* Matrix) {
    if (!Matrix->is_allocated()) {
        // Do we have enough memory to fit the entire matrix in core?
        size_t memory_required = Matrix->get_memory2();
        make_space(memory_required);
        Matrix->load();
    }
}

void CCBLAS::load_irrep(CCMatrix* Matrix, int h) {
    if (!Matrix->is_block_allocated(h)) {
        // Do we have enough memory to fit the entire matrix in core?
        size_t memory_required = Matrix->get_memorypi2(h);
        make_space(memory_required);
        Matrix->load_irrep(h);
    }
}

void CCBLAS::make_space(size_t memory_required) {
    if (memory_required < wfn_->free_memory_)
        return;
    else {
        outfile->Printf("\nCCBLAS::make_space() not implemented yet!!!");
        // Attempt #1
    }
}

}  // namespace psimrcc
}  // namespace psi
