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

#include "psi4/libiwl/iwl.hpp"
#include "psi4/libqt/qt.h"
#include "psi4/libdpd/dpd.h"
#include "psi4/libmints/wavefunction.h"
#include "defines.h"
#include "arrays.h"
#include "dpd.h"
#include "psi4/libpsi4util/PsiOutStream.h"

#include <cmath>

namespace psi {
namespace occwave {

/********************************************************************************************/
/************************** SymBlockMatrix **************************************************/
/********************************************************************************************/
SymBlockMatrix::SymBlockMatrix() {
    matrix_ = nullptr;
    rowspi_ = nullptr;
    colspi_ = nullptr;
}

SymBlockMatrix::SymBlockMatrix(std::string name, int nirreps, int *ins_rowspi, int *ins_colspi) {
    matrix_ = nullptr;
    name_ = name;
    nirreps_ = nirreps;
    rowspi_ = new int[nirreps_];
    colspi_ = new int[nirreps_];
    for (int h = 0; h < nirreps_; h++) {
        rowspi_[h] = ins_rowspi[h];
        colspi_[h] = ins_colspi[h];
    }
    memalloc();
}  //

SymBlockMatrix::~SymBlockMatrix() {
    release();
    if (rowspi_) delete[] rowspi_;
    if (colspi_) delete[] colspi_;
}  //

void SymBlockMatrix::memalloc() {
    if (matrix_) release();
    matrix_ = (double ***)malloc(sizeof(double ***) * nirreps_);
    for (int h = 0; h < nirreps_; h++) {
        if (rowspi_[h] != 0 && colspi_[h] != 0) {
            matrix_[h] = block_matrix(rowspi_[h], colspi_[h]);
        } else
            matrix_[h] = nullptr;
    }
}  //

void SymBlockMatrix::release() {
    if (!matrix_) return;
    for (int h = 0; h < nirreps_; h++) {
        if (matrix_[h]) free_block(matrix_[h]);
    }
    matrix_ = nullptr;
}  //

void SymBlockMatrix::zero() {
    size_t size;
    for (int h = 0; h < nirreps_; h++) {
        size = rowspi_[h] * colspi_[h] * sizeof(double);
        if (size) {
            memset(&(matrix_[h][0][0]), 0, size);
        }
    }
}  //

void SymBlockMatrix::set(dpdbuf4 G) {
    for (int h = 0; h < nirreps_; ++h) {
        global_dpd_->buf4_mat_irrep_init(&G, h);
        global_dpd_->buf4_mat_irrep_rd(&G, h);
        for (int row = 0; row < G.params->rowtot[h]; ++row) {
            for (int col = 0; col < G.params->coltot[h]; ++col) {
                matrix_[h][row][col] = G.matrix[h][row][col];
            }
        }
        global_dpd_->buf4_mat_irrep_close(&G, h);
    }
}  //

double SymBlockMatrix::get(int h, int m, int n) { return matrix_[h][m][n]; }  //

/********************************************************************************************/
/************************** SymBlockVector **************************************************/
/********************************************************************************************/
SymBlockVector::SymBlockVector() {
    vector_ = nullptr;
    dimvec_ = nullptr;
}

SymBlockVector::SymBlockVector(std::string name) {
    vector_ = nullptr;
    dimvec_ = nullptr;
}

SymBlockVector::SymBlockVector(int nirreps, int *ins_dimvec) {
    vector_ = nullptr;
    nirreps_ = nirreps;
    dimvec_ = new int[nirreps_];
    for (int h = 0; h < nirreps_; h++) {
        dimvec_[h] = ins_dimvec[h];
    }
    memalloc();
}  //

SymBlockVector::SymBlockVector(std::string name, int nirreps, int *ins_dimvec) {
    vector_ = nullptr;
    nirreps_ = nirreps;
    name_ = name;
    dimvec_ = new int[nirreps_];
    for (int h = 0; h < nirreps_; h++) {
        dimvec_[h] = ins_dimvec[h];
    }
    memalloc();
}  //

SymBlockVector::~SymBlockVector() {
    release();
    if (dimvec_) delete[] dimvec_;
}  //

int *SymBlockVector::dimvec() { return dimvec_; }  //

void SymBlockVector::memalloc() {
    if (vector_) release();
    vector_ = (double **)malloc(sizeof(double **) * nirreps_);
    for (int h = 0; h < nirreps_; h++) {
        if (dimvec_[h] != 0) {
            vector_[h] = new double[dimvec_[h]];
        } else
            vector_[h] = nullptr;
    }
}  //

void SymBlockVector::release() {
    if (!vector_) return;
    for (int h = 0; h < nirreps_; h++) {
        if (vector_[h]) free(vector_[h]);
    }
    vector_ = nullptr;
}  //

void SymBlockVector::zero() {
    size_t size;
    for (int h = 0; h < nirreps_; h++) {
        size = dimvec_[h] * sizeof(double);
        if (size) {
            memset(&(vector_[h][0]), 0.0, size);
        }
    }
}  //

void SymBlockVector::copy(const SymBlockVector *Adum) {
    // Make sure that vectors are in the same size
    bool same = true;
    for (int h = 0; h < nirreps_; h++) {
        if (dimvec_[h] != Adum->dimvec_[h]) same = false;
    }

    if (same == false) {
        release();
        if (dimvec_) delete[] dimvec_;
        dimvec_ = new int[nirreps_];
        for (int i = 0; i < nirreps_; ++i) {
            dimvec_[i] = Adum->dimvec_[i];
        }
        memalloc();
    }

    // If vectors are in the same size
    for (int h = 0; h < nirreps_; h++) {
        if (dimvec_[h] != 0) {
            memcpy(&(vector_[h][0]), &(Adum->vector_[h][0]), dimvec_[h] * sizeof(double));
        }
    }
}  //

void SymBlockVector::add(const SymBlockVector *Adum) {
    double *lhs, *rhs;
    for (int h = 0; h < nirreps_; h++) {
        size_t size = dimvec_[h];
        if (size) {
            lhs = vector_[h];
            rhs = Adum->vector_[h];
            for (size_t cnt = 0; cnt < size; cnt++) {
                *lhs += *rhs;
                lhs++;
                rhs++;
            }
        }
    }
}  //

void SymBlockVector::add(int h, int i, double value) { vector_[h][i] += value; }  //

void SymBlockVector::subtract(const SymBlockVector *Adum) {
    double *lhs, *rhs;
    for (int h = 0; h < nirreps_; h++) {
        size_t size = dimvec_[h];
        if (size) {
            lhs = vector_[h];
            rhs = Adum->vector_[h];
            for (size_t cnt = 0; cnt < size; cnt++) {
                *lhs -= *rhs;
                lhs++;
                rhs++;
            }
        }
    }
}  //

void SymBlockVector::subtract(int h, int i, double value) { vector_[h][i] -= value; }  //

void SymBlockVector::scale(double a) {
    size_t size;
    for (int h = 0; h < nirreps_; h++) {
        size = dimvec_[h];
        if (size) C_DSCAL(size, a, &(vector_[h][0]), 1);
    }
}  //

double SymBlockVector::sum_of_squares() {
    double summ;
    summ = 0.0;
    for (int h = 0; h < nirreps_; h++) {
        for (int j = 0; j < dimvec_[h]; ++j) {
            summ += vector_[h][j] * vector_[h][j];
        }
    }
    return summ;
}  //

double SymBlockVector::rms() {
    double summ;
    int dim;
    summ = 0.0;
    dim = 0;

    for (int h = 0; h < nirreps_; h++) {
        if (dimvec_[h] != 0) dim += dimvec_[h];
    }

    for (int h = 0; h < nirreps_; h++) {
        for (int j = 0; j < dimvec_[h]; ++j) {
            summ += vector_[h][j] * vector_[h][j];
        }
    }
    summ = std::sqrt(summ) / dim;
    return summ;
}  //

double SymBlockVector::rms(SymBlockVector *Atemp) {
    double summ;
    int dim;
    summ = 0.0;
    dim = 0;

    for (int h = 0; h < nirreps_; h++) {
        if (dimvec_[h] != 0) dim += dimvec_[h];
    }

    for (int h = 0; h < nirreps_; h++) {
        for (int j = 0; j < dimvec_[h]; ++j) {
            summ += (vector_[h][j] * vector_[h][j]) - (Atemp->vector_[h][j] * Atemp->vector_[h][j]);
        }
    }
    summ = std::sqrt(summ) / dim;
    return summ;
}  //

double SymBlockVector::norm() {
    double summ;
    summ = 0.0;
    for (int h = 0; h < nirreps_; h++) {
        for (int j = 0; j < dimvec_[h]; ++j) {
            summ += vector_[h][j] * vector_[h][j];
        }
    }
    summ = std::sqrt(summ);
    return summ;
}  //

void SymBlockVector::set(double value) {
    size_t size;
    for (int h = 0; h < nirreps_; h++) {
        size = dimvec_[h];
        for (size_t i = 0; i < size; ++i) {
            vector_[h][i] = value;
        }
    }
}  //

void SymBlockVector::set(int h, int i, double value) { vector_[h][i] = value; }  //

void SymBlockVector::set(double *Avec) {
    int offset;
    if (Avec == nullptr) return;

    offset = 0;
    for (int h = 0; h < nirreps_; h++) {
        for (int i = 0; i < dimvec_[h]; ++i) {
            int ii = i + offset;
            vector_[h][i] = Avec[ii];
        }
        offset += dimvec_[h];
    }
}  //

double SymBlockVector::get(int h, int m) { return vector_[h][m]; }  //

double *SymBlockVector::to_vector() {
    int sizecol = 0;
    for (int h = 0; h < nirreps_; h++) {
        sizecol += dimvec_[h];
    }

    double *temp = new double[sizecol];
    int offsetcol = 0;
    for (int h = 0; h < nirreps_; h++) {
        for (int j = 0; j < dimvec_[h]; ++j) {
            temp[j + offsetcol] = vector_[h][j];
        }
        offsetcol += dimvec_[h];
    }

    return temp;
}  //

double SymBlockVector::trace() {
    double value;
    value = 0.0;
    for (int h = 0; h < nirreps_; h++) {
        if (dimvec_[h] != 0) {
            for (int j = 0; j < dimvec_[h]; ++j) value += vector_[h][j];
        }
    }
    return value;
}  //

void SymBlockVector::print(std::string out_fname) {
    auto mode = std::ostream::app;
    auto printer = out_fname == "outfile" ? outfile : std::make_shared<PsiOutStream>(out_fname, mode);
    if (name_.length()) printer->Printf("\n ## %s ##\n", name_.c_str());
    for (int h = 0; h < nirreps_; h++) {
        if (dimvec_[h] != 0) {
            printer->Printf("\n Irrep: %d\n", h + 1);
            for (int j = 0; j < dimvec_[h]; ++j) {
                printer->Printf("%20.14f \n", vector_[h][j]);
            }
        }
    }
}  //

void SymBlockVector::print() {
    if (name_.length()) outfile->Printf("\n ## %s ##\n", name_.c_str());
    for (int h = 0; h < nirreps_; h++) {
        if (dimvec_[h] != 0) {
            outfile->Printf("\n Irrep: %d\n", h + 1);
            for (int j = 0; j < dimvec_[h]; ++j) {
                outfile->Printf("%20.14f \n", vector_[h][j]);
            }
        }
    }

}  //

void SymBlockVector::set_to_unit() {
    size_t size;
    for (int h = 0; h < nirreps_; h++) {
        size = dimvec_[h] * sizeof(double);
        if (size) {
            memset(&(vector_[h][0]), 0.0, size);
            for (int i = 0; i < dimvec_[h]; i++) vector_[h][i] = 1.0;
        }
    }
}  //

void SymBlockVector::gemv(bool transa, double alpha, SymBlockMatrix *A, SymBlockVector *X, double beta) {
    char trans = transa ? 't' : 'n';

    for (int h = 0; h < nirreps_; ++h) {
        C_DGEMV(trans, A->rowspi_[h], A->colspi_[h], alpha, &(A->matrix_[h][0][0]), A->rowspi_[h], &(X->vector_[h][0]),
                1, beta, &(vector_[h][0]), 1);
    }
}  //

double SymBlockVector::dot(SymBlockVector *X) {
    double tmp = 0.0;
    for (int h = 0; h < nirreps_; ++h) {
        if (dimvec_[h] != X->dimvec_[h]) {
            printf("SymBlockVector::dot: Vectors are not of the same size.\n");
            return 0.0;
        }
        for (int i = 0; i < dimvec_[h]; ++i) {
            tmp += vector_[h][i] * X->vector_[h][i];
        }
    }
    return tmp;
}  //

/********************************************************************************************/
/********************************************************************************************/
}
}  // End Namespaces
