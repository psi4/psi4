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

/*
 *  matrix.cc
 *  matrix
 *
 *  Created by Justin Turney on 4/1/08.
 *
 */

#include "psi4/libpsi4util/exception.h"
#include "psi4/libpsi4util/libpsi4util.h"
#include "psi4/libciomr/libciomr.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libiwl/iwl.hpp"
#include "psi4/libqt/qt.h"
#include "psi4/libparallel/parallel.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/integral.h"
#include "psi4/libdpd/dpd.h"
#include "psi4/libparallel/ParallelPrinter.h"
#include "factory.h"
#include "wavefunction.h"
#include "dimension.h"
#include "molecule.h"
#include "pointgrp.h"
#include "petitelist.h"

#include "psi4/pybind11.h"

#include <cmath>
#include <cstdio>
#include <fstream>
#include <algorithm>
#include <ctype.h>
#include <sstream>
#include <string>
#include <regex>
#include <tuple>

// In molecule.cc
namespace psi {
extern int str_to_int(const std::string &s);

extern double str_to_double(const std::string &s);

#define MIN(x, y) ((x) < (y) ? (x) : (y))
#define MAX(x, y) ((x) > (y) ? (x) : (y))

Matrix::Matrix()
{
    matrix_ = NULL;
    nirrep_ = 0;
    symmetry_ = 0;
}

Matrix::Matrix(const std::string &name, int symmetry)
        : matrix_(NULL), nirrep_(0),
          name_(name), symmetry_(symmetry)
{
}

Matrix::Matrix(const Matrix &c)
        : rowspi_(c.rowspi_), colspi_(c.colspi_)
{
    matrix_ = NULL;
    nirrep_ = c.nirrep_;
    symmetry_ = c.symmetry_;
    alloc();
    copy_from(c.matrix_);
}

Matrix::Matrix(const SharedMatrix &c)
        : rowspi_(c->rowspi_), colspi_(c->colspi_)
{
    matrix_ = NULL;
    nirrep_ = c->nirrep_;
    symmetry_ = c->symmetry_;
    alloc();
    copy_from(c->matrix_);
}

Matrix::Matrix(const Matrix *c)
        : rowspi_(c->rowspi_), colspi_(c->colspi_)
{
    matrix_ = NULL;
    nirrep_ = c->nirrep_;
    symmetry_ = c->symmetry_;
    alloc();
    copy_from(c->matrix_);
}

Matrix::Matrix(int l_nirreps, const int *l_rowspi, const int *l_colspi, int symmetry)
        : rowspi_(l_nirreps), colspi_(l_nirreps)
{
    matrix_ = NULL;
    nirrep_ = l_nirreps;
    symmetry_ = symmetry;
    rowspi_ = l_rowspi;
    colspi_ = l_colspi;
    alloc();
}

Matrix::Matrix(const std::string &name, int l_nirreps, const int *l_rowspi, const int *l_colspi, int symmetry)
        : rowspi_(l_nirreps), colspi_(l_nirreps), name_(name)
{
    matrix_ = NULL;
    nirrep_ = l_nirreps;
    symmetry_ = symmetry;
    rowspi_ = l_rowspi;
    colspi_ = l_colspi;
    alloc();
}

Matrix::Matrix(const std::string &name, int rows, int cols)
        : rowspi_(1), colspi_(1), name_(name)
{
    matrix_ = NULL;
    nirrep_ = 1;
    symmetry_ = 0;
    rowspi_[0] = rows;
    colspi_[0] = cols;
    alloc();
}

Matrix::Matrix(int rows, int cols)
        : rowspi_(1), colspi_(1)
{
    matrix_ = NULL;
    nirrep_ = 1;
    symmetry_ = 0;
    rowspi_[0] = rows;
    colspi_[0] = cols;
    alloc();
}

Matrix::Matrix(int nirrep, int rows, const int *colspi)
        : rowspi_(nirrep), colspi_(nirrep)
{
    matrix_ = NULL;
    symmetry_ = 0;
    nirrep_ = nirrep;
    for (int i = 0; i < nirrep_; ++i) {
        rowspi_[i] = rows;
        colspi_[i] = colspi[i];
    }
    alloc();
}

Matrix::Matrix(int nirrep, const int *rowspi, int cols)
        : rowspi_(nirrep), colspi_(nirrep)
{
    matrix_ = NULL;
    symmetry_ = 0;
    nirrep_ = nirrep;
    for (int i = 0; i < nirrep_; ++i) {
        rowspi_[i] = rowspi[i];
        colspi_[i] = cols;
    }
    alloc();
}

Matrix::Matrix(const std::string &name, const Dimension &rows, const Dimension &cols, int symmetry)
{
    name_ = name;
    matrix_ = NULL;
    symmetry_ = symmetry;

    // This will happen in PetiteList::aotoso()
    if (rows.n() == 1) {
        nirrep_ = cols.n();
        rowspi_ = Dimension(nirrep_);
        colspi_ = Dimension(nirrep_);
        for (int i = 0; i < nirrep_; ++i) {
            rowspi_[i] = rows[0];
            colspi_[i] = cols[i];
        }
    } else {
        nirrep_ = rows.n();
        rowspi_ = Dimension(nirrep_);
        colspi_ = Dimension(nirrep_);
        for (int i = 0; i < nirrep_; ++i) {
            rowspi_[i] = rows[i];
            colspi_[i] = cols[i];
        }
    }

    alloc();
}

Matrix::Matrix(const Dimension &rows, const Dimension &cols, int symmetry)
{
    matrix_ = NULL;
    symmetry_ = symmetry;

    // This will happen in PetiteList::aotoso()
    if (rows.n() == 1) {
        nirrep_ = cols.n();
        rowspi_ = Dimension(nirrep_);
        colspi_ = Dimension(nirrep_);
        for (int i = 0; i < nirrep_; ++i) {
            rowspi_[i] = rows[0];
            colspi_[i] = cols[i];
        }
    } else {
        nirrep_ = rows.n();
        rowspi_ = Dimension(nirrep_);
        colspi_ = Dimension(nirrep_);
        for (int i = 0; i < nirrep_; ++i) {
            rowspi_[i] = rows[i];
            colspi_[i] = cols[i];
        }
    }

    alloc();
}

Matrix::Matrix(dpdfile2 *inFile)
        : rowspi_(inFile->params->nirreps), colspi_(inFile->params->nirreps), name_(inFile->label)
{
    global_dpd_->file2_mat_init(inFile);
    global_dpd_->file2_mat_rd(inFile);
    matrix_ = NULL;
    symmetry_ = inFile->my_irrep;
    nirrep_ = inFile->params->nirreps;
    for (int i = 0; i < nirrep_; ++i) {
        rowspi_[i] = inFile->params->rowtot[i];
        colspi_[i] = inFile->params->coltot[i];
    }
    alloc();
    copy_from(inFile->matrix);
    global_dpd_->file2_mat_close(inFile);
}

Matrix::~Matrix()
{
    release();
}

/// allocate a block matrix -- analogous to libciomr's block_matrix
double **Matrix::matrix(int nrow, int ncol)
{
    double **mat = (double **) malloc(sizeof(double *) * nrow);
    const size_t size = sizeof(double) * nrow * ncol;
    mat[0] = (double *) malloc(size);
    ::memset((void *) mat[0], 0, size);
    for (int r = 1; r < nrow; ++r) mat[r] = mat[r - 1] + ncol;
    return mat;
}

/// free a (block) matrix -- analogous to libciomr's free_block
void Matrix::free(double **Block)
{
    ::free(Block[0]);
    ::free(Block);
}

void Matrix::init(int l_nirreps, const int *l_rowspi, const int *l_colspi, const std::string &name, int symmetry)
{
    name_ = name;
    symmetry_ = symmetry;
    nirrep_ = l_nirreps;
    rowspi_ = Dimension(nirrep_);
    colspi_ = Dimension(nirrep_);
    for (int i = 0; i < nirrep_; ++i) {
        rowspi_[i] = l_rowspi[i];
        colspi_[i] = l_colspi[i];
    }
    alloc();
}

void Matrix::init(const Dimension &l_rowspi, const Dimension &l_colspi, const std::string &name, int symmetry)
{
    name_ = name;
    symmetry_ = symmetry;
    nirrep_ = l_rowspi.n();
    rowspi_ = Dimension(nirrep_);
    colspi_ = Dimension(nirrep_);
    for (int i = 0; i < nirrep_; ++i) {
        rowspi_[i] = l_rowspi[i];
        colspi_[i] = l_colspi[i];
    }
    alloc();
}

SharedMatrix Matrix::create(const std::string &name,
                            const Dimension &rows,
                            const Dimension &cols)
{
    return SharedMatrix(new Matrix(name, rows, cols));
}

SharedMatrix Matrix::clone() const
{
    SharedMatrix temp(new Matrix(this));
    return temp;
}

void Matrix::copy(const Matrix *cp)
{
    // Make sure we are the same size as cp
    bool same = true;
    if (nirrep_ != cp->nirrep_ || symmetry_ != cp->symmetry_) {
        same = false;
    } else {
        if (colspi_ != cp->colspi_ || rowspi_ != cp->rowspi_)
            same = false;
    }

    if (same == false) {
        release();
        nirrep_ = cp->nirrep_;
        symmetry_ = cp->symmetry_;
        rowspi_ = Dimension(nirrep_);
        colspi_ = Dimension(nirrep_);
        for (int i = 0; i < nirrep_; ++i) {
            rowspi_[i] = cp->rowspi_[i];
            colspi_[i] = cp->colspi_[i];
        }
        alloc();
    }

    // When here we are the same size
#pragma omp parallel for
    for (int h = 0; h < nirrep_; ++h) {
        if (rowspi_[h] != 0 && colspi_[h ^ symmetry_] != 0)
            memcpy(&(matrix_[h][0][0]), &(cp->matrix_[h][0][0]), rowspi_[h] * (size_t) colspi_[h ^ symmetry_] * sizeof(double));
    }
}

SharedMatrix Matrix::horzcat(const std::vector <SharedMatrix> &mats)
{
    int nirrep = mats[0]->nirrep();
    for (size_t a = 0; a < mats.size(); ++a) {
        if (nirrep != mats[a]->nirrep()) {
            throw PSIEXCEPTION("Horzcat: Matrices not of same nirrep");
        }
    }

    for (size_t a = 1; a < mats.size(); ++a) {
        for (int h = 0; h < nirrep; ++h) {
            if (mats[a]->rowspi()[h] != mats[0]->rowspi()[h]) {
                throw PSIEXCEPTION("Horzcat: Matrices must all have same row dimension");
            }
        }
    }

    Dimension colspi(nirrep);

    for (size_t a = 0; a < mats.size(); ++a) {
        colspi += mats[a]->colspi();
    }

    SharedMatrix cat(new Matrix("", nirrep, mats[0]->rowspi(), colspi));

    for (int h = 0; h < nirrep; ++h) {
        if (mats[0]->rowspi()[h] == 0 || colspi[h] == 0) continue;
        double **catp = cat->pointer(h);
        int offset = 0;
        int rows = mats[0]->rowspi()[h];
        for (size_t a = 0; a < mats.size(); ++a) {
            int cols = mats[a]->colspi()[h];
            if (cols == 0) continue;

            double **Ap = mats[a]->pointer(h);

            for (int col = 0; col < cols; ++col) {
                C_DCOPY(rows, &Ap[0][col], cols, &catp[0][col + offset], colspi[h]);
            }

            offset += cols;
        }
    }

    return cat;
}

SharedMatrix Matrix::vertcat(const std::vector <SharedMatrix> &mats)
{
    int nirrep = mats[0]->nirrep();
    for (size_t a = 0; a < mats.size(); ++a) {
        if (nirrep != mats[a]->nirrep()) {
            throw PSIEXCEPTION("Vertcat: Matrices not of same nirrep");
        }
    }

    for (size_t a = 1; a < mats.size(); ++a) {
        for (int h = 0; h < nirrep; ++h) {
            if (mats[a]->colspi()[h] != mats[0]->colspi()[h]) {
                throw PSIEXCEPTION("Vertcat: Matrices must all have same col dimension");
            }
        }
    }

    Dimension rowspi(nirrep);

    for (size_t a = 0; a < mats.size(); ++a) {
        rowspi += mats[a]->rowspi();
    }

    SharedMatrix cat(new Matrix("", nirrep, rowspi, mats[0]->colspi()));

    for (int h = 0; h < nirrep; ++h) {
        if (mats[0]->colspi()[h] == 0 || rowspi[h] == 0) continue;
        double **catp = cat->pointer(h);
        int offset = 0;
        int cols = mats[0]->colspi()[h];
        for (size_t a = 0; a < mats.size(); ++a) {
            int rows = mats[a]->rowspi()[h];
            if (rows == 0) continue;

            double **Ap = mats[a]->pointer(h);

            for (int row = 0; row < rows; ++row) {
                ::memcpy((void *) catp[row + offset], (void *) Ap[row], sizeof(double) * cols);
            }

            offset += rows;
        }
    }

    return cat;
}

SharedMatrix Matrix::matrix_3d_rotation(Vector3 axis, double phi, bool Sn)
{

    if (ncol() != 3)
        throw PSIEXCEPTION("Can only rotate matrix with 3d vectors");

    // Normalize rotation vector
    double norm = sqrt(axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2]);
    axis[0] /= norm;
    axis[1] /= norm;
    axis[2] /= norm;

    double wx, wy, wz, cp;
    wx = axis[0];
    wy = axis[1];
    wz = axis[2];
    cp = 1.0 - cos(phi);

    Matrix R("Rotation Matrix", 3, 3);
    R(0, 0) = cos(phi) + wx * wx * cp;
    R(0, 1) = -wz * sin(phi) + wx * wy * cp;
    R(0, 2) = wy * sin(phi) + wx * wz * cp;
    R(1, 0) = wz * sin(phi) + wx * wy * cp;
    R(1, 1) = cos(phi) + wy * wy * cp;
    R(1, 2) = -wx * sin(phi) + wy * wz * cp;
    R(2, 0) = -wy * sin(phi) + wx * wz * cp;
    R(2, 1) = wx * sin(phi) + wy * wz * cp;
    R(2, 2) = cos(phi) + wz * wz * cp;

    //  R * coord^t = R_coord^t or coord * R^t = R_coord
    Matrix rotated_coord(nrow(), 3);
    rotated_coord.gemm(false, true, 1.0, *this, R, 0.0);

    if (Sn) { // delta_ij - 2 a_i a_j / ||a||^2
        R.identity();
        R(0, 0) -= 2 * wx * wx;
        R(1, 1) -= 2 * wy * wy;
        R(2, 2) -= 2 * wz * wz;
        R(1, 0) = R(0, 1) = 2 * wx * wy;
        R(2, 0) = R(0, 2) = 2 * wx * wz;
        R(2, 1) = R(1, 2) = 2 * wy * wz;
        Matrix tmp(nrow(), 3);
        tmp.gemm(false, true, 1.0, rotated_coord, R, 0.0);
        rotated_coord.copy(tmp);
    }

    SharedMatrix to_return = rotated_coord.clone();
    return to_return;
}


void Matrix::copy_to_row(int h, int row, double const *const data)
{
    if (h >= nirrep_ || row >= rowspi_[h])
        throw PSIEXCEPTION("Matrix::copy_to_row: Out of bounds.");

    memcpy(matrix_[h][row], data, sizeof(double) * colspi_[h]);
}

void Matrix::copy(const Matrix &cp)
{
    copy(&cp);
}

void Matrix::copy(const SharedMatrix &cp)
{
    copy(cp.get());
}

void Matrix::alloc()
{
    if (matrix_)
        release();

    // This is probably a default constructor matrix
    if (!nirrep_) {
        matrix_ = NULL;
        return;
    }

    matrix_ = (double ***) malloc(sizeof(double ***) * nirrep_);
    for (int h = 0; h < nirrep_; ++h) {
        if (rowspi_[h] != 0 && colspi_[h ^ symmetry_] != 0)
            matrix_[h] = Matrix::matrix(rowspi_[h], colspi_[h ^ symmetry_]);
        else {
            // Force rowspi_[h] and colspi_[h^symmetry] to hard 0
            // This solves an issue where a row can have 0 dim but a col does not (or the other way).

            // This was commented out to resolve issues that people were
            // dependent on one or both containing valid dimensions and
            // not a hard zero.
            //rowspi_[h] = colspi_[h^symmetry_] = 0;
            matrix_[h] = NULL;
        }
    }
}

void Matrix::release()
{
    if (!matrix_)
        return;

    for (int h = 0; h < nirrep_; ++h) {
        if (matrix_[h])
            Matrix::free(matrix_[h]);
    }
    ::free(matrix_);
    matrix_ = NULL;
}

void Matrix::copy_from(double ***c)
{
    for (int h = 0; h < nirrep_; ++h) {
        size_t size = rowspi_[h] * colspi_[h ^ symmetry_] * sizeof(double);
        if (size)
            memcpy(&(matrix_[h][0][0]), &(c[h][0][0]), size);
    }
}

// Sets all elements of matrix to val
void Matrix::set(double val)
{
    for (int h = 0; h < nirrep_; ++h) {
        size_t size = rowspi_[h] * colspi_[h ^ symmetry_];

        for (size_t i = 0; i < size; ++i) {
            matrix_[h][0][i] = val;
        }
    }
}

void Matrix::set(const double *const tri)
{
    int h, i, j, ii, jj;
    int row_offset;

    row_offset = 0;
    for (h = 0; h < nirrep_; ++h) {
        for (i = 0; i < rowspi_[h]; ++i) {
            ii = i + row_offset;

            if (symmetry_ == 0) {
                for (j = 0; j <= i; ++j) {
                    jj = j + row_offset;
                    matrix_[h][i][j] = matrix_[h][j][i] = tri[ii * (ii + 1) / 2 + jj];
                }
            } else {
                int col_offset = 0;
                for (int g = 0; g < (h ^ symmetry_); ++g)
                    col_offset += colspi_[g];

                for (j = 0; j < colspi_[h ^ symmetry_]; ++j) {
                    jj = j + col_offset;
                    matrix_[h][i][j] = tri[ii * (ii + 1) / 2 + jj];
                    matrix_[h ^ symmetry_][j][i] = matrix_[h][i][j];
                }
            }
        }
        row_offset += rowspi_[h];
    }
}

void Matrix::set(const double *const *const sq)
{
    if (symmetry_) {
        throw PSIEXCEPTION("Matrix::set called on a non-totally symmetric matrix.");
    }

    int h, i, j, ii, jj;
    int offset;

    if (sq == NULL) {
        throw PSIEXCEPTION("Matrix::set: Set call with a NULL double** matrix");
    }
    offset = 0;
    for (h = 0; h < nirrep_; ++h) {
        for (i = 0; i < rowspi_[h]; ++i) {
            ii = i + offset;
            for (j = 0; j <= i; ++j) {
                jj = j + offset;
                matrix_[h][i][j] = sq[ii][jj];
                matrix_[h][j][i] = sq[jj][ii];
            }
        }
        offset += rowspi_[h];
    }
}

void Matrix::set(const double *const *const sq, int h)
{
    if (sq == NULL)
        throw PSIEXCEPTION("Matrix::set: Set call with a NULL double** matrix");

    for (int i = 0; i < rowspi_[h]; i++)
        for (int j = 0; j < colspi_[h]; j++)
            matrix_[h][i][j] = sq[i][j];
}

void Matrix::set_diagonal(const Vector *const vec)
{
    if (symmetry_) {
        throw PSIEXCEPTION("Matrix::set_diagonal called on a non-totally symmetric matrix.");
    }

    int h, i, size;
    zero();
    for (h = 0; h < nirrep_; ++h) {
        size = rowspi_[h];
        if (size) {
            for (i = 0; i < size; ++i)
                matrix_[h][i][i] = vec->vector_[h][i];
        }
    }
}

void Matrix::set_diagonal(const Vector &vec)
{
    if (symmetry_) {
        throw PSIEXCEPTION("Matrix::set_diagonal called on a non-totally symmetric matrix.");
    }

    int h, i, size;
    zero();
    for (h = 0; h < nirrep_; ++h) {
        size = rowspi_[h];
        if (size) {
            for (i = 0; i < size; ++i)
                matrix_[h][i][i] = vec.vector_[h][i];
        }
    }
}

void Matrix::set_diagonal(const std::shared_ptr <Vector> &vec)
{
    if (symmetry_) {
        throw PSIEXCEPTION("Matrix::set_diagonal called on a non-totally symmetric matrix.");
    }

    int h, i, size;
    zero();
    for (h = 0; h < nirrep_; ++h) {
        size = rowspi_[h];
        if (size) {
            for (i = 0; i < size; ++i)
                matrix_[h][i][i] = vec->vector_[h][i];
        }
    }
}

SharedVector Matrix::get_row(int h, int m)
{
    if (m >= rowspi_[h]) {
        throw PSIEXCEPTION("Matrix::set_row: index is out of bounds.");
    }
    SharedVector vec = SharedVector(new Vector("Row", colspi_));
    vec->zero();
    size_t size = colspi_[h];
    for (size_t i = 0; i < size; ++i) {
        vec->vector_[h][i] = matrix_[h][m][i];
    }
    return vec;
}

SharedVector Matrix::get_column(int h, int m)
{
    if (m >= colspi_[h]) {
        throw PSIEXCEPTION("Matrix::get_column: index is out of bounds.");
    }
    SharedVector vec = SharedVector(new Vector("Column", rowspi_));
    vec->zero();
    size_t size = rowspi_[h];
    for (size_t i = 0; i < size; ++i) {
        vec->vector_[h][i] = matrix_[h][i][m];
    }
    return vec;
}

void Matrix::set_row(int h, int m, SharedVector vec)
{
    if (m >= rowspi_[h]) {
        throw PSIEXCEPTION("Matrix::set_row: index is out of bounds.");
    }
    size_t size = colspi_[h];
    for (size_t i = 0; i < size; ++i) {
        matrix_[h][m][i] = vec->vector_[h][i];
    }
}

void Matrix::set_column(int h, int m, SharedVector vec)
{
    if (m >= colspi_[h]) {
        throw PSIEXCEPTION("Matrix::set_column: index is out of bounds.");
    }
    size_t size = rowspi_[h];
    for (size_t i = 0; i < size; ++i) {
        matrix_[h][i][m] = vec->vector_[h][i];
    }
}

double *Matrix::to_lower_triangle() const
{
    int sizer = 0, sizec = 0;
    for (int h = 0; h < nirrep_; ++h) {
        sizer += rowspi_[h];
        sizec += colspi_[h ^ symmetry_];
    }
    if (sizer != sizec)
        return NULL;

    double *tri = new double[ioff[sizer]];
    double **temp = to_block_matrix();
    sq_to_tri(temp, tri, sizer);
    free_block(temp);
    return tri;
}

double **Matrix::to_block_matrix() const
{
    int sizer = 0, sizec = 0;
    for (int h = 0; h < nirrep_; ++h) {
        sizer += rowspi_[h];
        sizec += colspi_[h ^ symmetry_];
    }

    int *col_offset = new int[nirrep_];
    col_offset[0] = 0;
    for (int h = 1; h < nirrep_; ++h) {
        col_offset[h] = col_offset[h - 1] + colspi_[h - 1];
    }

    double **temp = block_matrix(sizer, sizec);
    int offsetr = 0, offsetc = 0;
    for (int h = 0; h < nirrep_; ++h) {
        offsetc = col_offset[h ^ symmetry_];
        for (int i = 0; i < rowspi_[h]; ++i) {
            for (int j = 0; j < colspi_[h ^ symmetry_]; ++j) {
                temp[i + offsetr][j + offsetc] = matrix_[h][i][j];
            }
        }
        offsetr += rowspi_[h];
//        offsetc += colspi_[h^symmetry_];
    }

    delete[] col_offset;
    return temp;
}

SharedMatrix Matrix::to_block_sharedmatrix() const
{
    int sizer = 0, sizec = 0;
    for (int h = 0; h < nirrep_; ++h) {
        sizer += rowspi_[h];
        sizec += colspi_[h ^ symmetry_];
    }
    SharedMatrix ret(new Matrix(name_ + " Block Copy", sizer, sizec));
    double **temp = to_block_matrix();
    ret->set(temp, 0);
    free_block(temp);
    return ret;

}

void Matrix::print_mat(const double *const *const a, int m, int n, std::string out) const
{
    std::shared_ptr <psi::PsiOutStream> printer = (out == "outfile" ? outfile :
                                                   std::shared_ptr<OutFile>(new OutFile(out)));

    const int print_ncol = Process::environment.options.get_int("MAT_NUM_COLUMN_PRINT");
    int num_frames = int(n / print_ncol);
    int num_frames_rem = n % print_ncol; //adding one for changing 0->1 start
    int num_frame_counter = 0;
    //for each frame
    for (num_frame_counter = 0; num_frame_counter < num_frames; num_frame_counter++) {
        printer->Printf("\n");
        for (int j = print_ncol * num_frame_counter + 1; j < print_ncol * num_frame_counter + print_ncol + 1; j++) {
            if (j == print_ncol * num_frame_counter + 1) { printer->Printf("%18d", j); }
            else { printer->Printf("               %5d", j); }
        }
        printer->Printf("\n\n");

        for (int k = 1; k <= m; ++k) {
            for (int j = print_ncol * num_frame_counter + 1; j < print_ncol * num_frame_counter + print_ncol + 2; j++) {
                if (j == print_ncol * num_frame_counter + 1) { printer->Printf("%5d", k); }
                else { printer->Printf(" %20.14f", a[k - 1][j - 2]); }
            }
            printer->Printf("\n");
        }
    }

    // ALREADY DID THE FULL FRAMES BY THIS POINT
    // NEED TO TAKE CARE OF THE REMAINDER
    if (num_frames_rem != 0) {
        printer->Printf("\n");
        for (int j = print_ncol * num_frame_counter + 1; j <= n; j++) {
            if (j == print_ncol * num_frame_counter + 1) { printer->Printf("%18d", j); }
            else { printer->Printf("               %5d", j); }
        }
        printer->Printf("\n\n");

        for (int k = 1; k <= m; ++k) {
            for (int j = print_ncol * num_frame_counter + 1; j < n + 2; j++) {
                if (j == print_ncol * num_frame_counter + 1) { printer->Printf("%5d", k); }
                else { printer->Printf(" %20.14f", a[k - 1][j - 2]); }
            }
            printer->Printf("\n");
        }
    }
    printer->Printf("\n\n");
    //R.I.P. goto statements - Aug 4th 2010 - MSM
}

void Matrix::print(std::string out, const char *extra) const
{
    int h;
    std::shared_ptr <psi::PsiOutStream> printer = (out == "outfile" ? outfile :
                                                   std::shared_ptr<OutFile>(new OutFile(out)));
    if (name_.length()) {
        if (extra == NULL)
            printer->Printf("  ## %s (Symmetry %d) ##\n", name_.c_str(), symmetry_);
        else
            printer->Printf("  ## %s %s (Symmetry %d)##\n", name_.c_str(), extra, symmetry_);
    }

    for (h = 0; h < nirrep_; ++h) {
        printer->Printf("  Irrep: %d Size: %d x %d\n", h + 1, rowspi_[h], colspi_[h ^ symmetry_]);
        if (rowspi_[h] == 0 || colspi_[h ^ symmetry_] == 0)
            printer->Printf("\n\t(empty)\n");
        else
            print_mat(matrix_[h], rowspi_[h], colspi_[h ^ symmetry_], out);
        printer->Printf("\n");
    }
}

void Matrix::print_to_mathematica()
{
    if (name_.length())
        outfile->Printf("  ## %s in Mathematica form ##\n", name_.c_str());
    else
        outfile->Printf("  ## Request matrix in Mathematica form ##\n");

    outfile->Printf("{");
    for (int h = 0; h < nirrep_; ++h) {
        outfile->Printf("{");

        for (int r = 0; r < rowspi_[h]; ++r) {
            outfile->Printf("{");
            for (int c = 0; c < colspi_[h ^ symmetry_]; ++c) {
                outfile->Printf("%14.12lf", get(h, r, c));
                if (c < colspi_[h] - 1)
                    outfile->Printf(", ");
            }
            outfile->Printf("}");
            if (r < rowspi_[h] - 1)
                outfile->Printf(",\n");
        }

        outfile->Printf("}");
        if (h < nirrep_ - 1)
            outfile->Printf(",\n");
    }
    outfile->Printf("}\n");
}

void Matrix::print_atom_vector(std::string out)
{
    int i;
    std::shared_ptr <psi::PsiOutStream> printer = (out == "outfile" ? outfile :
                                                   std::shared_ptr<OutFile>(new OutFile(out)));
    if (name_.length()) {
        printer->Printf("\n  -%s:\n", name_.c_str());
    }
    printer->Printf("     Atom            X                  Y                   Z\n");
    printer->Printf("    ------   -----------------  -----------------  -----------------\n");

    for (i = 0; i < nrow(); i++) {
        printer->Printf("    %4d   ", i + 1);
        printer->Printf("  %17.12lf  %17.12lf  %17.12lf", matrix_[0][i][0], matrix_[0][i][1], matrix_[0][i][2]);
        printer->Printf("\n");
    }
    printer->Printf("\n");
}


void Matrix::eivprint(const Vector *const values, std::string out)
{
    std::shared_ptr <psi::PsiOutStream> printer = (out == "outfile" ? outfile :
                                                   std::shared_ptr<OutFile>(new OutFile(out)));
    if (symmetry_)
        throw PSIEXCEPTION("Matrix::eivprint: This print does not make sense for non-totally symmetric matrices.");

    int h;

    if (name_.length()) {
        printer->Printf("  ## %s with eigenvalues ##\n", name_.c_str());
    }

    for (h = 0; h < nirrep_; ++h) {
        printer->Printf(" Irrep: %d\n", h + 1);
        eivout(matrix_[h], values->vector_[h], rowspi_[h], colspi_[h ^ symmetry_], out);
        printer->Printf("\n");
    }
}

void Matrix::eivprint(const Vector &values, std::string out)
{
    eivprint(&values, out);
}

void Matrix::eivprint(const std::shared_ptr <Vector> &values, std::string out)
{
    eivprint(values.get(), out);
}

void Matrix::symmetrize_gradient(std::shared_ptr <Molecule> molecule)
{
    if (nirrep_ > 1 || rowspi_[0] != molecule->natom() || colspi_[0] != 3)
        throw PSIEXCEPTION("Molecule::symmetrize_gradient: Matrix cannot be symmetrized.");

    // Symmetrize the gradients to remove any noise:
    CharacterTable ct = molecule->point_group()->char_table();

    // Obtain atom mapping of atom * symm op to atom
    int **atom_map = compute_atom_map(molecule);

    SharedMatrix ret(clone());
    ret->zero();
    Matrix temp = *this;

    // Symmetrize the gradients to remove any noise
    for (int atom = 0; atom < molecule->natom(); ++atom) {
        for (int g = 0; g < ct.order(); ++g) {

            int Gatom = atom_map[atom][g];

            SymmetryOperation so = ct.symm_operation(g);

            ret->add(atom, 0, so(0, 0) * temp(Gatom, 0) / ct.order());
            ret->add(atom, 0, so(0, 1) * temp(Gatom, 1) / ct.order());
            ret->add(atom, 0, so(0, 2) * temp(Gatom, 2) / ct.order());

            ret->add(atom, 1, so(1, 0) * temp(Gatom, 0) / ct.order());
            ret->add(atom, 1, so(1, 1) * temp(Gatom, 1) / ct.order());
            ret->add(atom, 1, so(1, 2) * temp(Gatom, 2) / ct.order());

            ret->add(atom, 2, so(2, 0) * temp(Gatom, 0) / ct.order());
            ret->add(atom, 2, so(2, 1) * temp(Gatom, 1) / ct.order());
            ret->add(atom, 2, so(2, 2) * temp(Gatom, 2) / ct.order());
        }
    }
    delete_atom_map(atom_map, molecule);
    copy(ret);
    ret.reset();
}

void Matrix::identity()
{
    if (symmetry_)
        return;

    int h;
    size_t size;

    for (h = 0; h < nirrep_; ++h) {
        size = rowspi_[h] * colspi_[h] * sizeof(double);

        if (size) {
            memset(&(matrix_[h][0][0]), 0, size);
            for (int i = 0; i < MIN(rowspi_[h], colspi_[h]); ++i)
                matrix_[h][i][i] = 1.0;
        }
    }
}

void Matrix::zero()
{
    size_t size;
    int h;

    for (h = 0; h < nirrep_; ++h) {
        size = rowspi_[h] * ((size_t) colspi_[h ^ symmetry_]) * sizeof(double);

        if (size) {
            memset(&(matrix_[h][0][0]), 0, size);
        }
    }
}

void Matrix::zero_diagonal()
{
    if (symmetry_)
        return;

    int h, i;

    for (h = 0; h < nirrep_; ++h) {
        for (i = 0; i < MIN(rowspi_[h], colspi_[h]); ++i) {
            matrix_[h][i][i] = 0.0;
        }
    }
}

double Matrix::trace()
{
    if (symmetry_)
        return 0.0;

    int i, h;
    double val = (double) 0.0;

    for (h = 0; h < nirrep_; ++h) {
        for (i = 0; i < MIN(rowspi_[h], colspi_[h]); ++i) {
            val += matrix_[h][i][i];
        }
    }

    return val;
}

SharedMatrix Matrix::transpose()
{
    SharedMatrix temp(new Matrix(name_, nirrep_, colspi_, rowspi_, symmetry_));

    if (symmetry_) {

        for (int rowsym = 0; rowsym < nirrep_; ++rowsym) {
            int colsym = rowsym ^symmetry_;
            if (rowsym < colsym) continue;
            int rows = rowspi_[rowsym];
            int cols = colspi_[colsym];
            for (int row = 0; row < rows; row++) {
                for (int col = 0; col < cols; col++) {
                    temp->matrix_[colsym][col][row] = matrix_[rowsym][row][col];
                    temp->matrix_[rowsym][row][col] = matrix_[colsym][col][row];
                }
            }
        }
    } else {
        int h, i, j;
        for (h = 0; h < nirrep_; ++h) {
            for (i = 0; i < rowspi_[h]; ++i) {
                for (j = 0; j < colspi_[h]; ++j) {
                    temp->matrix_[h][j][i] = matrix_[h][i][j];
                }
            }
        }
    }

    return temp;
}

void Matrix::transpose_this()
{
    if (symmetry_) {
        for (int rowsym = 0; rowsym < nirrep_; ++rowsym) {
            int colsym = rowsym ^symmetry_;
            if (rowsym < colsym) continue;
            int rows = rowspi_[rowsym];
            int cols = colspi_[colsym];
            for (int row = 0; row < rows; row++) {
                for (int col = 0; col < cols; col++) {
                    std::swap(matrix_[colsym][col][row], matrix_[rowsym][row][col]);
                }
            }
        }
    } else {
        int h, i, j;
        if (rowspi_ == colspi_) {
            for (h = 0; h < nirrep_; ++h) {
                for (i = 0; i < rowspi_[h]; ++i) {
                    for (j = 0; j < i; ++j) {
                        std::swap(matrix_[h][i][j], matrix_[h][j][i]);
                    }
                }
            }
        } else {
            throw NOT_IMPLEMENTED_EXCEPTION();
        }
    }
}

void Matrix::add(const Matrix *const plus)
{
    for (int h = 0; h < nirrep_; ++h) {
        size_t size = rowspi_[h] * colspi_[h ^ symmetry_];
        if (size) {
            C_DAXPY(size, 1.0, plus->matrix_[h][0], 1, matrix_[h][0], 1);
        }
    }
}

void Matrix::add(const Matrix &plus)
{
    add(&plus);
}

void Matrix::add(const SharedMatrix &plus)
{
    add(plus.get());
}

void Matrix::subtract(const Matrix *const plus)
{
    for (int h = 0; h < nirrep_; ++h) {
        size_t size = rowspi_[h] * colspi_[h ^ symmetry_];
        if (size) {
            C_DAXPY(size, -1.0, plus->matrix_[h][0], 1, matrix_[h][0], 1);
        }
    }
}

void Matrix::subtract(const SharedMatrix &sub)
{
    subtract(sub.get());
}

void Matrix::apply_denominator(const Matrix *const plus)
{
    double *lhs, *rhs;
    for (int h = 0; h < nirrep_; ++h) {
        size_t size = rowspi_[h] * colspi_[h ^ symmetry_];
        if (size) {
            lhs = matrix_[h][0];
            rhs = plus->matrix_[h][0];
#pragma omp parallel for simd
            for (size_t ij = 0; ij < size; ++ij) {
                lhs[ij] /= rhs[ij];
            }
        }
    }
}

void Matrix::apply_denominator(const Matrix &plus)
{
    apply_denominator(&plus);
}

void Matrix::apply_denominator(const SharedMatrix &plus)
{
    apply_denominator(plus.get());
}

void Matrix::accumulate_product(const Matrix *const a, const Matrix *const b)
{
    gemm(false, false, 1.0, a, b, 1.0);
}

void Matrix::accumulate_product(const SharedMatrix &a,
                                const SharedMatrix &b)
{
    gemm(false, false, 1.0, a, b, 1.0);
}

void Matrix::scale(double a)
{
    int h;
    size_t size;
    for (h = 0; h < nirrep_; ++h) {
        size = rowspi_[h] * (size_t) colspi_[h ^ symmetry_];
        if (size)
            C_DSCAL(size, a, &(matrix_[h][0][0]), 1);
    }
}

void Matrix::scale_row(int h, int m, double a)
{
    C_DSCAL(colspi_[h ^ symmetry_], a, &(matrix_[h][m][0]), 1);
}

void Matrix::scale_column(int h, int n, double a)
{
    C_DSCAL(rowspi_[h], a, &(matrix_[h][0][n]), colspi_[h ^ symmetry_]);
}

double Matrix::sum_of_squares()
{
    double sum = (double) 0.0;
    for (int h = 0; h < nirrep_; ++h) {
#pragma omp parallel for reduction(+:sum)
        for (int i = 0; i < rowspi_[h]; ++i) {
            for (int j = 0; j < colspi_[h ^ symmetry_]; ++j) {
                sum += matrix_[h][i][j] * matrix_[h][i][j];
            }
        }
    }

    return sum;
}

double Matrix::rms()
{
    double sum = (double) 0.0;
    long terms = 0;
    for (int h = 0; h < nirrep_; ++h) {
#pragma omp parallel for reduction(+:sum)
        for (int i = 0; i < rowspi_[h]; ++i) {
            for (int j = 0; j < colspi_[h ^ symmetry_]; ++j) {
                sum += matrix_[h][i][j] * matrix_[h][i][j];
                terms++;
            }
        }
    }

    return sqrt(sum / terms);
}

double Matrix::absmax()
{
    double max = (double) 0.0;
    for (int h = 0; h < nirrep_; ++h) {
        for (int i = 0; i < rowspi_[h]; ++i) {
            for (int j = 0; j < colspi_[h ^ symmetry_]; ++j) {
                double absval = std::abs(matrix_[h][i][j]);
                if (absval > max) max = absval;
            }
        }
    }

    return max;
}

void Matrix::transform(const Matrix *const a, const Matrix *const transformer)
{
#ifdef PSIDEBUG
    // Check dimensions
    // 'this' should be transformer->colspi by transformer->colspi
    if (rowspi_ != transformer->colspi() || colspi_ != transformer->colspi())
        throw PSIEXCEPTION("Matrix::transformer(a, transformer): Target matrix does not have correct dimensions.");
#endif

    Matrix temp(a->rowspi(), transformer->colspi());
    temp.gemm(false, false, 1.0, a, transformer, 0.0);
    gemm(true, false, 1.0, transformer, &temp, 0.0);
}

void Matrix::transform(const SharedMatrix &a, const SharedMatrix &transformer)
{
    transform(a.get(), transformer.get());
}

void Matrix::transform(const Matrix *const transformer)
{
    Matrix temp(nirrep_, rowspi_, transformer->colspi());
    temp.gemm(false, false, 1.0, this, transformer, 0.0);

    // Might need to resize the target matrix.
    if (rowspi() != transformer->rowspi() || colspi() != transformer->colspi())
        init(transformer->colspi(), transformer->colspi(), name_, symmetry_);

    gemm(true, false, 1.0, transformer, &temp, 0.0);
}

void Matrix::transform(const SharedMatrix &transformer)
{
    transform(transformer.get());
}

void Matrix::transform(const SharedMatrix &L,
                       const SharedMatrix &F,
                       const SharedMatrix &R)
{
#ifdef PSIDEBUG
    // Check dimensions
    // 'this' should be transformer->colspi by transformer->colspi
    if (rowspi_ != L->colspi() || colspi_ != R->colspi())
        throw PSIEXCEPTION("Matrix::transformer(L, F, R): Target matrix does not have correct dimensions.");
#endif

    Matrix temp(nirrep_, F->rowspi_, R->colspi_, F->symmetry_ ^ R->symmetry_);
    temp.gemm(false, false, 1.0, F, R, 0.0);
    gemm(true, false, 1.0, L, temp, 0.0);
}

void Matrix::back_transform(const Matrix *const a, const Matrix *const transformer)
{
    Matrix temp(transformer->rowspi(), a->colspi());

    temp.gemm(false, false, 1.0, transformer, a, 0.0);
    gemm(false, true, 1.0, &temp, transformer, 0.0);
}

void Matrix::back_transform(const SharedMatrix &a, const SharedMatrix &transformer)
{
    back_transform(a.get(), transformer.get());
}

void Matrix::back_transform(const Matrix *const transformer)
{
    bool square = true;
    int h = 0;

    while (h < nirrep_ && square) {
        if (transformer->rowspi()[h] != transformer->colspi()[h]) {
            square = false;
        }
        h++;
    }

    if (square) {
        Matrix temp("", rowspi_, colspi_);
        temp.gemm(false, true, 1.0, this, transformer, 0.0);
        gemm(false, false, 1.0, transformer, &temp, 0.0);
    } else {
        Matrix temp(nirrep_, rowspi_, transformer->rowspi());
        Matrix result(nirrep_, transformer->rowspi(), transformer->rowspi());
        temp.gemm(false, true, 1.0, this, transformer, 0.0);
        result.gemm(false, false, 1.0, transformer, &temp, 0.0);
        copy(&result);
    }
}

void Matrix::back_transform(const SharedMatrix &transformer)
{
    back_transform(transformer.get());
}

void Matrix::gemm(const char &transa, const char &transb,
                  const std::vector<int> &m,
                  const std::vector<int> &n,
                  const std::vector<int> &k,
                  const double &alpha,
                  const SharedMatrix &a, const std::vector<int> &lda,
                  const SharedMatrix &b, const std::vector<int> &ldb,
                  const double &beta,
                  const std::vector<int> &ldc,
                  const std::vector<unsigned long> &offset_a,
                  const std::vector<unsigned long> &offset_b,
                  const std::vector<unsigned long> &offset_c)
{
    // For now only handle symmetric matrices right now
    if (symmetry_ || a->symmetry_ || b->symmetry_)
        throw PSIEXCEPTION("Matrix::Advanced GEMM: Can only handle totally symmetric matrices.");

    if (nirrep_ != a->nirrep_ || nirrep_ != b->nirrep_)
        throw PSIEXCEPTION("Matrix::Advanced GEMM: Number of irreps do not equal.");

    for (int h = 0; h < nirrep_; ++h) {
        int offa, offb, offc;

        offa = offset_a.size() == 0 ? 0 : offset_a[h];
        offb = offset_b.size() == 0 ? 0 : offset_b[h];
        offc = offset_c.size() == 0 ? 0 : offset_c[h];

        C_DGEMM(transa, transb, m[h], n[h], k[h],
                alpha,
                &a->matrix_[h][0][offa], lda[h],
                &b->matrix_[h][0][offb], ldb[h],
                beta,
                &matrix_[h][0][offc], ldc[h]);
    }
}

void Matrix::gemm(const char &transa, const char &transb,
                  const int &m,
                  const int &n,
                  const int &k,
                  const double &alpha,
                  const SharedMatrix &a, const int &lda,
                  const SharedMatrix &b, const int &ldb,
                  const double &beta,
                  const int &ldc,
                  const unsigned long &offset_a,
                  const unsigned long &offset_b,
                  const unsigned long &offset_c)
{
#ifdef DEBUG
    if (nirrep_ > 1)
        throw PSIEXCEPTION("Matrix::Advanced GEMM: C1 version called on symmetry objects.");
#endif

    C_DGEMM(transa, transb, m, n, k,
            alpha,
            &a->matrix_[0][0][offset_a], lda,
            &b->matrix_[0][0][offset_b], ldb,
            beta,
            &matrix_[0][0][offset_c], ldc);
}

void Matrix::gemm(bool transa, bool transb, double alpha, const Matrix *const a,
                  const Matrix *const b, double beta)
{
    if (nirrep_ != a->nirrep_ || nirrep_ != b->nirrep_)
        throw PSIEXCEPTION("Matrix::gemm error: Number of irreps do not equal.");

    // Check symmetry
    if (symmetry_ != (a->symmetry_ ^ b->symmetry_)) {
        outfile->Printf("Matrix::gemm error: Input symmetries will not result in target symmetry.\n");
        outfile->Printf(" Asym %d ^ Bsym %d != Csym %d\n", a->symmetry(), b->symmetry(), symmetry());
        outfile->Printf("Result is %d\n", a->symmetry_ ^ b->symmetry_);
        throw PSIEXCEPTION("Matrix::gemm error: Input symmetries will not result in target symmetry.");
    }

    if (transa && a->symmetry_)
        throw PSIEXCEPTION("Matrix::gemm error: a is non totally symmetric and you're trying to transpose it");
    if (transb && b->symmetry_)
        throw PSIEXCEPTION("Matrix::gemm error: b is non totally symmetric and you're trying to transpose it");

    char ta = transa ? 't' : 'n';
    char tb = transb ? 't' : 'n';
    int h, m, n, k, lda, ldb, ldc;

    for (h = 0; h < nirrep_; ++h) {
        m = rowspi_[h];
        n = colspi_[h ^ symmetry_];
        k = transa ? a->rowspi_[h] : a->colspi_[h ^ a->symmetry_];

        lda = a->colspi_[h ^ a->symmetry_];
        ldb = b->colspi_[h ^ b->symmetry_];
        ldc = colspi_[h ^ symmetry_];

        if (m && n && k) {
            C_DGEMM(ta, tb, m, n, k, alpha, &(a->matrix_[h][0][0]),
                    lda, &(b->matrix_[h ^ symmetry_ ^ b->symmetry_][0][0]), ldb, beta, &(matrix_[h][0][0]),
                    ldc);
        }
    }
}

void Matrix::gemm(bool transa, bool transb, double alpha,
                  const SharedMatrix &a, const SharedMatrix &b,
                  double beta)
{
    gemm(transa, transb, alpha, a.get(), b.get(), beta);
}

void Matrix::gemm(bool transa, bool transb, double alpha,
                  const Matrix &a, const SharedMatrix &b,
                  double beta)
{
    gemm(transa, transb, alpha, &a, b.get(), beta);
}

void Matrix::gemm(bool transa, bool transb, double alpha,
                  const SharedMatrix &a, const Matrix &b,
                  double beta)
{
    gemm(transa, transb, alpha, a.get(), &b, beta);
}

void Matrix::gemm(bool transa, bool transb, double alpha,
                  const Matrix &a, const Matrix &b, double beta)
{
    gemm(transa, transb, alpha, &a, &b, beta);
}

SharedMatrix Matrix::doublet(const SharedMatrix &A, const SharedMatrix &B, bool transA, bool transB)
{
    if (A->symmetry() || B->symmetry()) {
        throw PSIEXCEPTION("Matrix::doublet is not supported for this non-totally-symmetric thing.");
    }

    if (A->nirrep() != B->nirrep()) {
        throw PSIEXCEPTION("Matrix::doublet: Matrices do not have the same nirreps");
    }

    int nirrep = A->nirrep();
    int m[nirrep];
    int n[nirrep];
    int k[nirrep];
    int k2[nirrep];

    for (int h = 0; h < nirrep; h++) {
        m[h] = (transA ? A->colspi()[h] : A->rowspi()[h]);
        n[h] = (transB ? B->rowspi()[h] : B->colspi()[h]);
        k[h] = (!transA ? A->colspi()[h] : A->rowspi()[h]);
        k2[h] = (!transB ? B->rowspi()[h] : B->colspi()[h]);
        if (k[h] != k2[h]) {
            throw PSIEXCEPTION("Matrix::doublet: Dimension mismatch");
        }
    }

    SharedMatrix T(new Matrix("T", nirrep, m, n));

    for (int h = 0; h < nirrep; h++) {
        if (!k[h] || !m[h] || !n[h]) continue;
        C_DGEMM(
                (transA ? 'T' : 'N'),
                (transB ? 'T' : 'N'),
                m[h],
                n[h],
                k[h],
                1.0,
                A->pointer(h)[0],
                A->colspi()[h],
                B->pointer(h)[0],
                B->colspi()[h],
                0.0,
                T->pointer(h)[0],
                T->colspi()[h]);
    }

    return T;
}

SharedMatrix Matrix::triplet(const SharedMatrix &A, const SharedMatrix &B, const SharedMatrix &C, bool transA, bool transB, bool transC)
{
    SharedMatrix T = Matrix::doublet(A, B, transA, transB);
    SharedMatrix S = Matrix::doublet(T, C, false, transC);
    return S;
}

void Matrix::axpy(double a, SharedMatrix X)
{
    if (nirrep_ != X->nirrep()) {
        throw PSIEXCEPTION("Matrix::axpy: Matrices do not have the same nirreps");
    }
    for (int h = 0; h < nirrep_; h++) {
        size_t size = colspi_[h] * rowspi_[h];
        if (size != (X->rowspi()[h] * X->colspi()[h])) {
            throw PSIEXCEPTION("Matrix::axpy: Matrices sizes do not match.");
        }
        if (size) {
            double *Xp = X->pointer(h)[0];
            double *Yp = matrix_[h][0];
            C_DAXPY(size, a, Xp, 1, Yp, 1);
        }

    }
}

SharedMatrix Matrix::collapse(int dim)
{
    if (dim < 0 || dim > 1) throw PSIEXCEPTION("Matrix::collapse: dim must be 0 (row sum) or 1 (col sum)");

    if (symmetry_) {
        throw PSIEXCEPTION("Matrix::collapse is not supported for this non-totally-symmetric thing.");
    }


    Dimension ones(nirrep_);
    for (int h = 0; h < nirrep_; h++) {
        ones[h] = 1;
    }

    std::shared_ptr <Matrix> T(new Matrix("T", ((dim == 0) ? colspi_ : rowspi_), ones));

    for (int h = 0; h < nirrep_; h++) {
        int nrow = rowspi_[h];
        int ncol = colspi_[h];
        double **Mp = matrix_[h];
        double **Tp = T->pointer(h);
        if (dim == 0) {
            for (int j = 0; j < ncol; j++) {
                for (int i = 0; i < nrow; i++) {
                    Tp[j][0] += Mp[i][j];
                }
            }
        } else {
            for (int i = 0; i < nrow; i++) {
                for (int j = 0; j < ncol; j++) {
                    Tp[i][0] += Mp[i][j];
                }
            }
        }
    }

    return T;
}

namespace {

int mat_schmidt_tol(double **C, double **S, int nrow, int ncol, double tolerance, double *res)
{
    int i, j;
    int northog = 0;
    double vtmp;
    std::vector<double> v(nrow);

    if (res) *res = 1.0;

    // Orthonormalize the columsn of this wrt S.
    std::fill(v.begin(), v.end(), 0.0);

    for (int m = 0; m < ncol; ++m) {
        v[0] = C[0][m] * S[0][0];

        for (i = 1; i < nrow; ++i) {
            for (j = 0, vtmp = 0.0; j < i; j++) {
                vtmp += C[j][m] * S[i][j];
                v[j] += C[i][m] * S[i][j];
            }
            v[i] = vtmp + C[i][m] * S[i][j];
        }

        for (i = 0, vtmp = 0.0; i < nrow; ++i)
            vtmp += v[i] * C[i][m];

        if (vtmp < tolerance) continue;

        if (res && (m == 0 || vtmp < *res)) *res = vtmp;

        vtmp = 1.0 / sqrt(vtmp);

        for (i = 0; i < nrow; ++i) {
            v[i] *= vtmp;
            C[i][northog] = C[i][m] * vtmp;
        }

        for (i = m + 1, vtmp = 0.0; i < ncol; ++i) {
            for (j = 0, vtmp = 0.0; j < nrow; ++j)
                vtmp += v[j] * C[j][i];
            for (j = 0; j < nrow; ++j)
                C[j][i] -= vtmp * C[j][northog];
        }

        northog++;
    }

    return northog;
}

}

void Matrix::schmidt()
{
    for (int h = 0; h < nirrep(); ++h) {
        if (!rowspi(h) || !colspi(h)) continue;
        psi::schmidt(matrix_[h], rowspi(h), colspi(h), "STUPID");
    }
}

Dimension Matrix::schmidt_orthog_columns(SharedMatrix S, double tol, double * /*res*/)
{
    Dimension northog(nirrep());
    std::vector<double> resid(nirrep());

    for (int h = 0; h < nirrep(); ++h) {
        northog[h] = mat_schmidt_tol(matrix_[h], S->matrix_[h], rowspi(h), colspi(h), tol, &resid[h]);
    }

    return northog;
}

bool Matrix::add_and_orthogonalize_row(const SharedVector v)
{
    Vector v_copy(*v.get());
    if (v_copy.nirrep() > 1 || nirrep_ > 1)
        throw PSIEXCEPTION("Matrix::schmidt_add_and_orthogonalize: Symmetry not allowed (yet).");
    if (v_copy.dimpi()[0] != colspi_[0])
        throw PSIEXCEPTION("Matrix::schmidt_add_and_orthogonalize: Incompatible dimensions.");
    double **mat = Matrix::matrix(rowspi_[0] + 1, colspi_[0]);
    size_t n = colspi_[0] * rowspi_[0] * sizeof(double);
    if (n) {
        ::memcpy(mat[0], matrix_[0][0], n);
        Matrix::free(matrix_[0]);
    }
    matrix_[0] = mat;
    bool ret = schmidt_add_row(0, rowspi_[0], v_copy);
    rowspi_[0]++;
    return ret;
}

bool Matrix::schmidt_add_row(int h, int rows, Vector &v) throw()
{
    if (v.nirrep() > 1)
        throw PSIEXCEPTION("Matrix::schmidt_add: This function needs to be adapted to handle symmetry blocks.");

    double dotval, normval;
    int i, I;

    for (i = 0; i < rows; ++i) {
        dotval = C_DDOT(coldim(h), matrix_[h][i], 1, v.pointer(), 1);
        for (I = 0; I < coldim(h); ++I)
            v(I) -= dotval * matrix_[h][i][I];
    }

    normval = C_DDOT(coldim(h), v.pointer(), 1, v.pointer(), 1);
    normval = sqrt(normval);

    if (normval > 1.0e-5) {
        for (I = 0; I < coldim(h); ++I)
            matrix_[h][rows][I] = v(I) / normval;
        return true;
    } else
        return false;
}

bool Matrix::schmidt_add_row(int h, int rows, double *v) throw()
{
    double dotval, normval;
    int i, I;

//    outfile->Printf( "in schmidt_add\n");
//    for (i=0; i<coldim(h); ++i)
//        outfile->Printf( "%lf ", v[i]);
//    outfile->Printf( "\n");

    for (i = 0; i < rows; ++i) {
        dotval = C_DDOT(coldim(h), matrix_[h][i], 1, v, 1);
        for (I = 0; I < coldim(h); ++I)
            v[I] -= dotval * matrix_[h][i][I];
    }

    normval = C_DDOT(coldim(h), v, 1, v, 1);
    normval = sqrt(normval);

    if (normval > 1.0e-5) {
        for (I = 0; I < coldim(h); ++I)
            matrix_[h][rows][I] = v[I] / normval;

//        for (i=0; i<coldim(h); ++i)
//            outfile->Printf( "%lf ", matrix_[h][rows][i]);
//        outfile->Printf( "\n");
        return true;
    } else
        return false;
}

void Matrix::project_out(Matrix &constraints)
{
    // We're going to work through temp and add to this
    Matrix temp = *this;
    zero();

//    outfile->Printf( "in project_out:\n");

    temp.set_name("temp");
//    temp.print();

//    constraints.print();

    double *v = new double[coldim()];
//    outfile->Printf( "coldim(): %d\n", coldim());
    for (int h = 0; h < nirrep(); ++h) {
        for (int i = 0; i < rowdim(h); ++i) {
//            outfile->Printf( "i=%d, copying %d elements from temp[%d][%d] to v\n", i, coldim(h), h, i);
            memcpy(v, temp[h][i], sizeof(double) * coldim(h));

//            outfile->Printf( "temp[%d][] ", h);
//            for(int z=0; z<coldim(h); ++z)
//                outfile->Printf( "%lf ", temp[h][i][z]);
//            outfile->Printf( "\n");

//            outfile->Printf( "v[] ", h);
//            for(int z=0; z<coldim(h); ++z)
//                outfile->Printf( "%lf ", v[z]);
//            outfile->Printf( "\n");

            for (int j = 0; j < constraints.rowdim(0); ++j) {
                // hand rolled ddot
                double dotval = 0.0;
                for (int z = 0; z < coldim(h); ++z) {
                    dotval += temp[h][i][z] * constraints[0][j][z];
//                    outfile->Printf( " %lf * %lf ", temp[h][i][z], constraints[0][j][z]);
                }
//                outfile->Printf( "\n");
//                double dotval = C_DDOT(coldim(h), &(temp[h][i][0]), 1, &(constraints[0][j][0]), 1);
//                outfile->Printf( "dotval = %lf\n", dotval);
                for (int I = 0; I < coldim(h); ++I)
                    v[I] -= dotval * constraints[0][j][I];
            }

//            outfile->Printf( "after removing constraints v[] ", h);
//            for(int z=0; z<coldim(h); ++z)
//                outfile->Printf( "%lf ", v[z]);
//            outfile->Printf( "\n");

            // At this point all constraints have been projected out of "v"
            // Normalize it add Schmidt orthogonalize it against this
            double normval = C_DDOT(coldim(h), v, 1, v, 1);
            if (normval > 1.0E-10) {
                normval = sqrt(normval);
                for (int j = 0; j < coldim(h); ++j)
                    v[j] /= normval;

//                outfile->Printf( "calling schmidt_add sending i=%d\n", i);
//                for(int z=0; z<coldim(h); ++z)
//                    outfile->Printf( "%lf ", v[z]);
//                outfile->Printf( "\n");
                schmidt_add_row(h, i, v);
            }
        }
    }

    delete[] v;
}

double Matrix::vector_dot(const Matrix *const rhs)
{
    if (symmetry_ != rhs->symmetry_)
        return 0.0;

    double sum = 0.0;
    int h;
    size_t size;

    for (h = 0; h < nirrep_; ++h) {
        size = rowspi_[h] * colspi_[h ^ symmetry_];
        // Check the size of the other
        if (size != (size_t)(rhs->rowdim(h) * rhs->coldim(h ^ symmetry_)))
            throw PSIEXCEPTION("Matrix::vector_dot: Dimensions do not match!\n");

        if (size)
            sum += C_DDOT(size, (&matrix_[h][0][0]), 1, &(rhs->matrix_[h][0][0]), 1);
    }

    return sum;
}

double Matrix::vector_dot(const SharedMatrix &rhs)
{
    return vector_dot(rhs.get());
}

void Matrix::diagonalize(Matrix *eigvectors, Vector *eigvalues, diagonalize_order nMatz)
{
    if (symmetry_) {
        throw PSIEXCEPTION("Matrix::diagonalize: Matrix is non-totally symmetric.");
    }
    int h;
    for (h = 0; h < nirrep_; ++h) {
        if (rowspi_[h]) {
            sq_rsp(rowspi_[h], colspi_[h], matrix_[h], eigvalues->vector_[h], static_cast<int>(nMatz), eigvectors->matrix_[h], 1.0e-14);
        }
    }
}

void Matrix::diagonalize(SharedMatrix &eigvectors, std::shared_ptr <Vector> &eigvalues, diagonalize_order nMatz)
{
    diagonalize(eigvectors.get(), eigvalues.get(), nMatz);
}

void Matrix::diagonalize(SharedMatrix &eigvectors, Vector &eigvalues, diagonalize_order nMatz)
{
    diagonalize(eigvectors.get(), &eigvalues, nMatz);
}

void Matrix::diagonalize(SharedMatrix &metric, SharedMatrix & /*eigvectors*/, std::shared_ptr <Vector> &eigvalues, diagonalize_order /*nMatz*/)
{
    if (symmetry_) {
        throw PSIEXCEPTION("Matrix::diagonalize: Matrix non-totally symmetric.");
    }

    // this and metric are destroyed in the process, so let's make a copy
    // that we work with.
    Matrix t(*this);
    Matrix m(metric);

    int lwork = 3 * max_nrow();
    double *work = new double[lwork];

    for (int h = 0; h < nirrep_; ++h) {
        if (!rowspi_[h] && !colspi_[h])
            continue;

        int err = C_DSYGV(1, 'V', 'U',
                          rowspi_[h], t.matrix_[h][0],
                          rowspi_[h], m.matrix_[h][0],
                          rowspi_[h], eigvalues->pointer(h),
                          work, lwork);

        if (err != 0) {
            if (err < 0) {
                outfile->Printf("Matrix::diagonalize with metric: C_DSYGV: argument %d has invalid parameter.\n", -err);

                abort();
            }
            if (err > 0) {
                outfile->Printf("Matrix::diagonalize with metric: C_DSYGV: error value: %d\n", err);

                abort();
            }
        }

        // TODO: Sort the data according to eigenvalues.
    }
    delete[] work;
}

std::tuple<SharedMatrix, SharedVector, SharedMatrix> Matrix::svd_temps()
{
    Dimension rank(nirrep_);
    for (int h = 0; h < nirrep_; h++) {
        int m = rowspi_[h];
        int n = colspi_[h ^ symmetry_];
        int k = (m < n ? m : n);
        rank[h] = k;
    }
    SharedMatrix U(new Matrix("U", rowspi_, rank));
    SharedVector S(new Vector("S", rank));
    SharedMatrix V(new Matrix("V", rank, colspi_));

    return std::tuple<SharedMatrix, SharedVector, SharedMatrix>(U, S, V);
}

std::tuple<SharedMatrix, SharedVector, SharedMatrix> Matrix::svd_a_temps()
{
    Dimension rank(nirrep_);
    for (int h = 0; h < nirrep_; h++) {
        int m = rowspi_[h];
        int n = colspi_[h ^ symmetry_];
        rank[h] = std::min(m, n);
    }
    SharedMatrix U(new Matrix("U", rowspi_, rowspi_));
    SharedVector S(new Vector("S", rank));
    SharedMatrix V(new Matrix("V", colspi_, colspi_));
    return std::tuple<SharedMatrix, SharedVector, SharedMatrix>(U, S, V);
}

void Matrix::svd(SharedMatrix &U, SharedVector &S, SharedMatrix &V)
{
    // Actually, this routine takes mn + mk + nk
    for (int h = 0; h < nirrep_; h++) {
        if (!rowspi_[h] || !colspi_[h ^ symmetry_])
            continue;

        int m = rowspi_[h];
        int n = colspi_[h ^ symmetry_];
        int k = (m < n ? m : n);

        double **Ap = Matrix::matrix(m, n);
        ::memcpy((void *) Ap[0], (void *) matrix_[h][0], sizeof(double) * m * n);
        double *Sp = S->pointer(h);
        double **Up = U->pointer(h);
        double **Vp = V->pointer(h ^ symmetry_);

        int *iwork = new int[8L * k];

        // Workspace Query
        double lwork;
        int info = C_DGESDD('S', n, m, Ap[0], n, Sp, Vp[0], n, Up[0], k, &lwork, -1, iwork);

        double *work = new double[(int) lwork];

        // SVD
        info = C_DGESDD('S', n, m, Ap[0], n, Sp, Vp[0], n, Up[0], k, work, (int) lwork, iwork);

        delete[] work;
        delete[] iwork;

        if (info != 0) {
            if (info < 0) {
                outfile->Printf("Matrix::svd with metric: C_DGESDD: argument %d has invalid parameter.\n", -info);

                abort();
            }
            if (info > 0) {
                outfile->Printf("Matrix::svd with metric: C_DGESDD: error value: %d\n", info);

                abort();
            }
        }
        Matrix::free(Ap);
    }
}

void Matrix::svd_a(SharedMatrix &U, SharedVector &S, SharedMatrix &V)
{
    // Actually, this routine takes mn + mk + nk
    for (int h = 0; h < nirrep_; h++) {
        int m = rowspi_[h];
        int n = colspi_[h ^ symmetry_];
        // There is something to SVD
        if ((m != 0) && (n != 0)) {

            int k = (m < n ? m : n);

            double **Ap = Matrix::matrix(m, n);
            ::memcpy((void *) Ap[0], (void *) matrix_[h][0], sizeof(double) * m * n);
            double *Sp = S->pointer(h);
            double **Up = U->pointer(h);
            double **Vp = V->pointer(h ^ symmetry_);

            int *iwork = new int[8L * k];

            // Workspace Query
            double lwork;
            int info = C_DGESDD('A', n, m, Ap[0], n, Sp, Vp[0], n, Up[0], m, &lwork, -1, iwork);

            double *work = new double[(int) lwork];

            // SVD
            info = C_DGESDD('A', n, m, Ap[0], n, Sp, Vp[0], n, Up[0], m, work, (int) lwork, iwork);

            delete[] work;
            delete[] iwork;

            if (info != 0) {
                if (info < 0) {
                    outfile->Printf("Matrix::svd with metric: C_DGESDD: argument %d has invalid parameter.\n", -info);

                    abort();
                }
                if (info > 0) {
                    outfile->Printf("Matrix::svd with metric: C_DGESDD: error value: %d\n", info);

                    abort();
                }
            }
            Matrix::free(Ap);
        } else if ((m != 0) && (n == 0)) {
            // There is nothing to SVD, but we need set the U block to the identity matrix
            double **Up = U->pointer(h);
            for (int i = 0; i < m; ++i) {
                for (int j = 0; j < m; ++j) {
                    Up[i][j] = 0.0;
                }
                Up[i][i] = 1.0;
            }
        } else if ((m == 0) && (n != 0)) {
            // There is nothing to SVD, but we need set the V block to the identity matrix
            double **Vp = V->pointer(h ^ symmetry_);
            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < n; ++j) {
                    Vp[i][j] = 0.0;
                }
                Vp[i][i] = 1.0;
            }
        }
    }
}

SharedMatrix Matrix::pseudoinverse(double condition, bool *conditioned)
{
    std::tuple<SharedMatrix, SharedVector, SharedMatrix> svd_temp = svd_temps();
    SharedMatrix U = std::get<0>(svd_temp);
    SharedVector S = std::get<1>(svd_temp);
    SharedMatrix V = std::get<2>(svd_temp);

    svd(U, S, V);

    bool conditioned_here = false;
    for (int h = 0; h < nirrep_; h++) {
        int ncol = S->dimpi()[h];
        double *Sp = S->pointer(h);
        double S0 = (ncol ? Sp[0] : 0.0);
        for (int i = 0; i < ncol; i++) {
            if (Sp[i] > S0 * condition) {
                Sp[i] = 1.0 / Sp[i];
            } else {
                Sp[i] = 0.0;
                conditioned_here = true;
            }
        }
    }
    if (conditioned) {
        *conditioned = conditioned_here;
    }

    SharedMatrix Q(clone());

    for (int h = 0; h < nirrep_; h++) {
        int m = rowspi_[h];
        int n = colspi_[h ^ symmetry_];
        int k = S->dimpi()[h];
        if (!m || !n || !k) continue;
        double **Up = U->pointer(h);
        double *Sp = S->pointer(h);
        double **Vp = V->pointer(h ^ symmetry_);
        double **Qp = Q->pointer(h);
        for (int i = 0; i < k; i++) {
            C_DSCAL(m, Sp[i], &Up[0][i], k);
        }
        C_DGEMM('N', 'N', m, n, k, 1.0, Up[0], k, Vp[0], n, 0.0, Qp[0], n);
    }

    return Q;
}

SharedMatrix Matrix::canonical_orthogonalization(double delta, SharedMatrix eigvec)
{
    if (symmetry_) {
        throw PSIEXCEPTION("Matrix: canonical orthogonalization only works for totally symmetric matrices");
    }

    SharedMatrix U(clone());
    SharedVector a(new Vector("a", rowspi_));

    diagonalize(U, a, descending);

    if (eigvec)
        eigvec->copy(U);

    Dimension rank(nirrep_);

    for (int h = 0; h < nirrep_; h++) {
        int k = a->dimpi()[h];
        if (!k) continue;
        int sig = 0;
        double *ap = a->pointer(h);
        double a0 = ap[0];
        for (int i = 0; i < k; i++) {
            if (ap[i] > a0 * delta) {
                ap[i] = pow(ap[i], -1.0 / 2.0);
                sig++;
            } else {
                ap[i] = 0.0;
            }
        }
        rank[h] = sig;
    }

    SharedMatrix X(new Matrix("X", rowspi_, rank));

    for (int h = 0; h < nirrep_; h++) {
        int m = rowspi_[h];
        int k = rank[h];
        if (!m || !k) continue;
        double **Up = U->pointer(h);
        double **Xp = X->pointer(h);
        double *ap = a->pointer(h);
        for (int i = 0; i < k; i++) {
            C_DAXPY(m, ap[i], &Up[0][i], m, &Xp[0][i], k);
        }
    }

    return X;
}

void Matrix::swap_rows(int h, int i, int j)
{
    C_DSWAP(colspi_[h], &(matrix_[h][i][0]), 1, &(matrix_[h][j][0]), 1);
}

void Matrix::swap_columns(int h, int i, int j)
{
    C_DSWAP(rowspi_[h], &(matrix_[h][0][i]), colspi_[h], &(matrix_[h][0][j]), colspi_[h]);
}

void Matrix::cholesky_factorize()
{
    if (symmetry_) {
        throw PSIEXCEPTION("Matrix::cholesky_factorize: Matrix is non-totally symmetric.");
    }
    for (int h = 0; h < nirrep_; ++h) {
        if (rowspi_[h]) {
            int err = C_DPOTRF('L', rowspi_[h], matrix_[h][0], rowspi_[h]);
            if (err != 0) {
                if (err < 0) {
                    outfile->Printf("cholesky_factorize: C_DPOTRF: argument %d has invalid paramter.\n", -err);

                    abort();
                }
                if (err > 1) {
                    outfile->Printf("cholesky_factorize: C_DPOTRF: the leading minor of order %d is not "
                                            "positive definite, and the factorization could not be "
                                            "completed.", err);

                    abort();
                }
            }
        }
    }
}

SharedMatrix Matrix::partial_cholesky_factorize(double delta, bool throw_if_negative)
{
    if (symmetry_) {
        throw PSIEXCEPTION("Matrix::partial_cholesky_factorize: Matrix is non-totally symmetric.");
    }

    // Temporary cholesky factor (full memory)
    SharedMatrix K(new Matrix("L Temp", nirrep_, rowspi_, rowspi_));

    // Significant Cholesky columns per irrep
    int *sigpi = new int[nirrep_];
    ::memset(static_cast<void *>(sigpi), '\0', nirrep_ * sizeof(int));

    for (int h = 0; h < nirrep_; ++h) {

        if (!rowspi_[h]) continue;

        // Dimension
        int n = rowspi_[h];
        // Cholesky factor
        double **Kp = K->pointer(h);
        // Original matrix
        double **Ap = matrix_[h];

        // Diagonal (or later Schur complement diagonal)
        double *Dp = new double[n];
        for (int i = 0; i < n; i++)
            Dp[i] = Ap[i][i];

        // Vector of completed columns (absolute)
        std::vector<int> order;

        int nQ = 0;
        int Q = -1;
        while (nQ < n) {

            // Find max error on diagonal
            // (Should always be in the Schur complement)
            int imax = 0;
            for (int i = 0; i < n; i++)
                if (fabs(Dp[i]) > fabs(Dp[imax]))
                    imax = i;

            double dmax = Dp[imax];
            if (fabs(dmax) <= delta)
                break;

            if (dmax <= 0.0) {
                if (throw_if_negative)
                    throw PSIEXCEPTION("Matrix::partial_cholesky_factorize: Pivot is numerically negative or zero");
                else
                    break;
            }

            // New vector!
            nQ++;
            Q++;

            // Find the diagonal
            double diag = sqrt(dmax);

            // Update the vector
            C_DCOPY(n, &Ap[0][imax], n, &Kp[0][Q], n);
            C_DGEMV('N', n, nQ - 1, -1.0, Kp[0], n, Kp[imax], 1, 1.0, &Kp[0][Q], n);
            C_DSCAL(n, 1.0 / diag, &Kp[0][Q], n);

            // Explicitly zero out elements of the vector
            // Which are psychologically upper triangular
            for (size_t i = 0; i < order.size(); i++)
                Kp[order[i]][Q] = 0.0;

            // Place the diagonal
            Kp[imax][Q] = diag;

            // Update the Schur complement
            for (int i = 0; i < n; i++)
                Dp[i] -= Kp[i][Q] * Kp[i][Q];

            // Explicitly zero out elements of the Schur complement
            // Which are already exact, and do not really belong
            // This prevents false selection due to roundoff
            Dp[imax] = 0.0;

            // Add the diagonal index to the list of completed indices
            order.push_back(imax);
        }
        sigpi[h] = nQ;

	delete [] Dp;
    }

    // Copy out to properly sized array
    SharedMatrix L(new Matrix("Partial Cholesky Factor", nirrep_, rowspi_, sigpi));

    //K->print();
    //L->print();

    for (int h = 0; h < nirrep_; h++) {
        if (!rowspi_[h] || !sigpi[h]) continue;
        double **Kp = K->pointer(h);
        double **Lp = L->pointer(h);

        for (int i = 0; i < rowspi_[h]; i++) {
            ::memcpy(static_cast<void *>(Lp[i]), static_cast<void *>(Kp[i]), sizeof(double) * sigpi[h]);
        }
    }

    delete[] sigpi;
    return L;
}

std::pair <SharedMatrix, SharedMatrix> Matrix::partial_square_root(double delta)
{
    if (symmetry_) {
        throw PSIEXCEPTION("Matrix::partial_square_root: Matrix is non-totally symmetric.");
    }

    SharedMatrix V(new Matrix("V", colspi_, colspi_));
    SharedVector d(new Vector("d", colspi_));

    diagonalize(V, d);

    Dimension Ppi(d->nirrep());
    Dimension Npi(d->nirrep());
    for (int h = 0; h < d->nirrep(); h++) {
        for (int i = 0; i < d->dimpi()[h]; i++) {
            if (fabs(d->get(h, i)) >= delta) {
                if (d->get(h, i) >= 0) {
                    Ppi[h]++;
                } else {
                    Npi[h]++;
                }
            }
        }
    }

    SharedMatrix P(new Matrix("P", colspi_, Ppi));
    SharedMatrix N(new Matrix("N", colspi_, Npi));

    for (int h = 0; h < d->nirrep(); h++) {
        double **Vp = V->pointer(h);
        double **Pp = P->pointer(h);
        double **Np = N->pointer(h);

        int Pcounter = 0;
        int Ncounter = 0;
        for (int i = 0; i < colspi_[h]; i++) {
            if (fabs(d->get(h, i)) >= delta) {
                if (d->get(h, i) >= 0.0) {
                    // +
                    double val = sqrt(fabs(d->get(h, i)));
                    C_DAXPY(colspi_[h], val, &Vp[0][i], colspi_[h], &Pp[0][Pcounter], Ppi[h]);
                    Pcounter++;
                } else {
                    // -
                    double val = sqrt(fabs(d->get(h, i)));
                    C_DAXPY(colspi_[h], -val, &Vp[0][i], colspi_[h], &Np[0][Ncounter], Npi[h]);
                    Ncounter++;
                }
            }
        }
    }

    return std::pair<SharedMatrix, SharedMatrix>(P, N);
}

void Matrix::invert()
{
    if (symmetry_) {
        throw PSIEXCEPTION("Matrix::invert: Matrix is non-totally symmetric.");
    }

    double **work = block_matrix(max_nrow(), max_ncol());
    for (int h = 0; h < nirrep_; ++h) {
        if (rowspi_[h] && colspi_[h ^ symmetry_] && rowspi_[h] == colspi_[h ^ symmetry_]) {
            invert_matrix(matrix_[h], work, rowspi_[h], "outfile");
            memcpy(&(matrix_[h][0][0]), &(work[0][0]), sizeof(double) * rowspi_[h] * colspi_[h]);
        }
    }
    free_block(work);
}

void Matrix::general_invert()
{
    if (symmetry_) {
        throw PSIEXCEPTION("Matrix::invert: Matrix is non-totally symmetric.");
    }

    int lwork = max_nrow() * max_ncol();
    double *work = new double[lwork];
    int *ipiv = new int[max_nrow()];

    for (int h = 0; h < nirrep_; ++h) {
        if (rowspi_[h] && colspi_[h]) {
            int err = C_DGETRF(rowspi_[h], colspi_[h], matrix_[h][0], rowspi_[h], ipiv);
            if (err != 0) {
                if (err < 0) {
                    outfile->Printf("invert: C_DGETRF: argument %d has invalid paramter.\n", -err);

                    abort();
                }
                if (err > 1) {
                    outfile->Printf("invert: C_DGETRF: the (%d,%d) element of the factor U or L is "
                                            "zero, and the inverse could not be computed.\n", err, err);

                    abort();
                }
            }

            err = C_DGETRI(colspi_[h], matrix_[h][0], rowspi_[h], ipiv, work, lwork);
            if (err != 0) {
                if (err < 0) {
                    outfile->Printf("invert: C_DGETRI: argument %d has invalid paramter.\n", -err);

                    abort();
                }
                if (err > 1) {
                    outfile->Printf("invert: C_DGETRI: the (%d,%d) element of the factor U or L is "
                                            "zero, and the inverse could not be computed.\n", err, err);

                    abort();
                }
            }
        }
    }
    delete[] ipiv;
    delete[] work;
}

Dimension Matrix::power(double alpha, double cutoff)
{
    if (symmetry_) {
        throw PSIEXCEPTION("Matrix::power: Matrix is non-totally symmetric.");
    }

    Dimension remaining(nirrep_, "Number of remaining orbitals");

    for (int h = 0; h < nirrep_; ++h) {
        if (rowspi_[h] == 0) continue;

        int n = rowspi_[h];
        double **A = matrix_[h];

        double **A1 = Matrix::matrix(n, n);
        double **A2 = Matrix::matrix(n, n);
        double *a = new double[n];

        memcpy(static_cast<void *>(A1[0]), static_cast<void *>(A[0]), sizeof(double) * n * n);

        // Eigendecomposition
        double lwork;
        int stat = C_DSYEV('V', 'U', n, A1[0], n, a, &lwork, -1);
        double *work = new double[(int) lwork];
        stat = C_DSYEV('V', 'U', n, A1[0], n, a, work, (int) lwork);
        delete[] work;

        if (stat)
            throw PSIEXCEPTION("Matrix::power: C_DSYEV failed");

        memcpy(static_cast<void *>(A2[0]), static_cast<void *>(A1[0]), sizeof(double) * n * n);

        double max_a = (fabs(a[n - 1]) > fabs(a[0]) ? fabs(a[n - 1]) : fabs(a[0]));
        int remain = 0;
        for (int i = 0; i < n; i++) {

            if (alpha < 0.0 && fabs(a[i]) < cutoff * max_a)
                a[i] = 0.0;
            else {
                a[i] = pow(a[i], alpha);
                if (std::isfinite(a[i])) {
                    remain++;
                } else {
                    a[i] = 0.0;
                }
            }

            C_DSCAL(n, a[i], A2[i], 1);
        }
        remaining[h] = remain;

        C_DGEMM('T', 'N', n, n, n, 1.0, A2[0], n, A1[0], n, 0.0, A[0], n);

        delete[] a;
        Matrix::free(A1);
        Matrix::free(A2);
    }

    return remaining;
}

void Matrix::expm(int m, bool scale)
{
    if (symmetry_) {
        throw PSIEXCEPTION("Matrix::expm: Matrix is non-totally symmetric.");
    }

    // Build the Pade Table
    std::vector<double> fact;
    fact.push_back(1.0);
    for (int i = 1; i <= 2 * m; i++) {
        fact.push_back(fact[i - 1] * i);
    }
    std::vector<double> alpha;
    for (int k = 0; k <= m; k++) {
        alpha.push_back((fact[2 * m - k] * fact[m]) / (fact[2 * m] * fact[k] * fact[m - k]));
    }

    // Form the exponential
    for (int h = 0; h < nirrep_; ++h) {
        if (rowspi_[h] == 0) continue;

        int n = rowspi_[h];
        double **A = matrix_[h];

        double L;
        int S;
        if (scale) {
            // Trace reduction
            L = 0.0;
            for (int i = 0; i < n; i++) {
                L += A[i][i];
            }
            L /= (double) n;
            for (int i = 0; i < n; i++) {
                A[i][i] -= L;
            }

            // Scaling
            double norm = 0.0;
            for (int i = 0; i < n; i++) {
                double row = 0.0;
                for (int j = 0; j < n; j++) {
                    row += fabs(A[i][j]);
                }
                norm = (norm > row ? norm : row);
            }
            double mag = log(norm) / log(2.0);
            mag = (mag < 0.0 ? 0.0 : mag);
            mag = (mag > 4.0 ? 4.0 : mag);
            S = (int) (mag);
            C_DSCAL(n * (unsigned long int) n, pow(2.0, -S), A[0], 1);
        }

        double **T = Matrix::matrix(n, n);
        double **U = Matrix::matrix(n, n);
        double **X = Matrix::matrix(n, n);
        double **Y = Matrix::matrix(n, n);

        // Zero-th Order
        for (int i = 0; i < n; i++) {
            X[i][i] = 1.0;
        }

        // Build X and Y as polynomials in A
        ::memcpy((void *) T[0], (void *) A[0], sizeof(double) * n * n);
        for (int Q = 1; Q <= m; Q++) {

            if ((Q % 2) == 1)
                C_DAXPY(n * (unsigned long int) n, alpha[Q], T[0], 1, Y[0], 1);
            else
                C_DAXPY(n * (unsigned long int) n, alpha[Q], T[0], 1, X[0], 1);

            if (Q == m) break;

            // T *= A
            C_DGEMM('N', 'N', n, n, n, 1.0, T[0], n, A[0], n, 0.0, U[0], n);
            double **t = T;
            T = U;
            U = t;
        }

        // Build N and D as polynomials in A (avoids cancelation)
        double **N = T;
        double **D = U;

        ::memcpy((void *) N[0], (void *) X[0], sizeof(double) * n * n);
        ::memcpy((void *) D[0], (void *) X[0], sizeof(double) * n * n);
        C_DAXPY(n * (unsigned long int) n, -1.0, Y[0], 1, N[0], 1);
        C_DAXPY(n * (unsigned long int) n, 1.0, Y[0], 1, D[0], 1);

        //outfile->Printf("  ## N ##\n\n");
        //print_mat(N,n,n,outfile);
        //outfile->Printf("  ## D ##\n\n");
        //print_mat(D,n,n,outfile);

        // Solve exp(A) = N / D = D^{1} N = D \ N
        int *ipiv = new int[n];

        // LU = D
        int info1 = C_DGETRF(n, n, D[0], n, ipiv);
        if (info1)
            throw PSIEXCEPTION("Matrix::expm: LU factorization of D failed");

        // Transpose N before solvation (FORTRAN)
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                double temp = N[j][i];
                N[j][i] = N[i][j];
                N[i][j] = temp;
            }
        }

        //outfile->Printf("  ## LU ##\n\n");
        //print_mat(D,n,n,outfile);
        //outfile->Printf("  ## S(0) ##\n\n");
        //print_mat(N,n,n,outfile);

        // D \ N
        int info2 = C_DGETRS('N', n, n, D[0], n, ipiv, N[0], n);
        if (info2)
            throw PSIEXCEPTION("Matrix::expm: LU solution of D failed");

        delete[] ipiv;

        //outfile->Printf("  ## S ##\n\n");
        //print_mat(N,n,n,outfile);

        if (scale) {

            // Copy result back to A, transposing as you go (back to C++)
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    X[i][j] = N[j][i];
                }
            }

            // Inverse scale
            for (int i = 0; i < S; i++) {
                C_DGEMM('N', 'N', n, n, n, 1.0, X[0], n, X[0], n, 0.0, Y[0], n);
                double **t = X;
                X = Y;
                Y = t;
            }
            ::memcpy((void *) A[0], (void *) X[0], sizeof(double) * n * n);

            // Inverse trace shift
            C_DSCAL(n * (unsigned long int) n, exp(L), A[0], 1);

        } else {
            // Copy result back to A, transposing as you go (back to C++)
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    A[i][j] = N[j][i];
                }
            }
        }

        Matrix::free(X);
        Matrix::free(Y);
        Matrix::free(T);
        Matrix::free(U);
    }
}

void Matrix::zero_lower()
{
    if (symmetry_) {
        throw PSIEXCEPTION("Matrix::zero_lower: Matrix is non-totally symmetric.");
    }

    for (int h = 0; h < nirrep_; ++h) {
#pragma omp parallel for schedule(dynamic) 
        for (int m = 0; m < rowspi_[h]; ++m) {
            for (int n = 0; n < m; ++n) {
                matrix_[h][m][n] = 0.0;
            }
        }
    }

}

void Matrix::zero_upper()
{
    if (symmetry_) {
        throw PSIEXCEPTION("Matrix::zero_upper: Matrix is non-totally symmetric.");
    }

    for (int h = 0; h < nirrep_; ++h) {
#pragma omp parallel for schedule(dynamic) 
        for (int m = 0; m < rowspi_[h]; ++m) {
            for (int n = 0; n < m; ++n) {
                matrix_[h][n][m] = 0.0;
            }
        }
    }

}

void Matrix::zero_row(int h, int i)
{
    if (i >= rowspi_[h]) {
        throw PSIEXCEPTION("Matrix::zero_row: index is out of bounds.");
    }
#pragma omp parallel for simd
    for (int m = 0; m < colspi_[h]; ++m) {
        matrix_[h][i][m] = 0.0;
    }
}

void Matrix::zero_column(int h, int i)
{
    if (i >= colspi_[h]) {
        throw PSIEXCEPTION("Matrix::zero_column: index is out of bounds.");
    }
#pragma omp parallel for
    for (int m = 0; m < rowspi_[h]; ++m) {
        matrix_[h][m][i] = 0.0;
    }
}

void Matrix::copy_lower_to_upper()
{
    if (symmetry_) {
        for (int rowsym = 0; rowsym < nirrep_; ++rowsym) {
            int colsym = rowsym ^symmetry_;
            if (rowsym < colsym) continue;
            int rows = rowspi_[rowsym];
            int cols = colspi_[colsym];
            for (int row = 0; row < rows; row++) {
                for (int col = 0; col < cols; col++) {
                    matrix_[colsym][col][row] = matrix_[rowsym][row][col];
                }
            }
        }
    } else {
        for (int h = 0; h < nirrep_; ++h) {
            for (int m = 0; m < rowspi_[h]; ++m) {
                for (int n = 0; n < m; ++n) {
                    matrix_[h][n][m] = matrix_[h][m][n];
                }
            }
        }
    }
}

void Matrix::copy_upper_to_lower()
{
    if (symmetry_) {
        for (int rowsym = 0; rowsym < nirrep_; ++rowsym) {
            int colsym = rowsym ^symmetry_;
            if (rowsym > colsym) continue;
            int rows = rowspi_[rowsym];
            int cols = colspi_[colsym];
            for (int row = 0; row < rows; row++) {
                for (int col = 0; col < cols; col++) {
                    matrix_[rowsym][row][col] = matrix_[colsym][col][row];
                }
            }
        }
    } else {
        for (int h = 0; h < nirrep_; ++h) {
            for (int m = 0; m < rowspi_[h]; ++m) {
                for (int n = 0; n < m; ++n) {
                    matrix_[h][m][n] = matrix_[h][n][m];
                }
            }
        }
    }
}

void Matrix::hermitivitize()
{
    if (symmetry_) {
        throw PSIEXCEPTION("Hermitivitize: matrix is not totally symmetric");
    }

    for (int h = 0; h < nirrep_; h++) {
        if (rowspi_[h] != colspi_[h]) {
            throw PSIEXCEPTION("Hermitivitize: matrix is not square");
        }
        int n = rowspi_[h];
        if (!n) continue;
        double **M = matrix_[h];

        for (int row = 0; row < n - 1; row++) {
            for (int col = row + 1; col < n; col++) {
                M[row][col] = M[col][row] = 0.5 * (M[row][col] + M[col][row]);
            }
        }
    }
}

// Reference versions of the above functions:

void Matrix::transform(const Matrix &a, const Matrix &transformer)
{
    // Allocate adaquate size temporary matrix.
    Matrix temp(a.rowspi(), transformer.colspi());
    temp.gemm(false, false, 1.0, a, transformer, 0.0);
    gemm(true, false, 1.0, transformer, temp, 0.0);
}

void Matrix::apply_symmetry(const SharedMatrix &a, const SharedMatrix &transformer)
{
    // Check dimensions of the two matrices and symmetry
    if (a->nirrep() > 1) {
        throw PSIEXCEPTION("Matrix::apply_symmetry: first matrix must have no symmetry.\n");
    }
    if (a->nrow() != transformer->rowdim(0)
        || a->ncol() != transformer->ncol()) {
        a->print();
        transformer->print();
        throw PSIEXCEPTION("Matrix::apply_symmetry: simple to regular. Sizes are not compatible.\n");
    }

    // Create temporary matrix of proper size.
    Matrix temp(nirrep(), a->nrow(), transformer->colspi());

    char ta = 'n';
    char tb = 'n';
    int m, n, k, nca, ncb, ncc;

    // Solve F = T^ M T

    // temp = M T
    for (int h = 0; h < nirrep_; ++h) {
        m = temp.rowdim(h);
        n = temp.coldim(h);
        k = a->ncol();
        nca = k;
        ncb = n;
        ncc = n;

        if (m && n && k) {
            C_DGEMM(ta, tb, m, n, k, 1.0, &(a->matrix_[0][0][0]),
                    nca, &(transformer->matrix_[h][0][0]), ncb,
                    0.0, &(temp.matrix_[h][0][0]), ncc);
        }
    }

    // F = T^ temp
    ta = 't';
    for (int h = 0; h < nirrep_; ++h) {
        m = rowdim(h);
        n = coldim(h);
        k = transformer->rowdim(h);
        nca = m;
        ncb = n;
        ncc = n;

        if (m && n && k) {
            C_DGEMM(ta, tb, m, n, k, 1.0, &(transformer->matrix_[h][0][0]),
                    nca, &(temp.matrix_[h][0][0]), ncb, 0.0, &(matrix_[h][0][0]),
                    ncc);
        }
    }
}

void Matrix::remove_symmetry(const SharedMatrix &a, const SharedMatrix &SO2AO)
{
    // Check dimensions of the two matrices and symmetry
    if (a->nirrep() != SO2AO->nirrep()) {
        throw PSIEXCEPTION("Matrix::remove_symmetry: matrices must have same symmetry.\n");
    }
    if (nirrep() != 1) {
        throw PSIEXCEPTION("Matrix::remove_symmetry: result matrix must not have symmetry. \n");
    }
    if (ncol() != SO2AO->coldim(0) ||
        a->nrow() != SO2AO->nrow()) {
        a->print();
        SO2AO->print();
        throw PSIEXCEPTION("Matrix::remove_symmetry: Sizes are not compatible.\n");
    }

    // Ensure we're working with a clea
    zero();

    // Create temporary matrix of proper size.
    Matrix temp(SO2AO->nirrep(), SO2AO->rowspi(), SO2AO->colspi());

    char ta = 'n';
    char tb = 'n';
    int m, n, k, nca, ncb, ncc;

    // Solve F = T^ M T

    // temp = M T
    for (int h = 0; h < SO2AO->nirrep(); ++h) {
        m = temp.rowdim(h);
        n = temp.coldim(h);
        k = a->coldim(h);
        nca = k;
        ncb = n;
        ncc = n;

        if (m && n && k) {
            C_DGEMM(ta, tb, m, n, k, 1.0, &(a->matrix_[h][0][0]),
                    nca, &(SO2AO->matrix_[h][0][0]), ncb,
                    1.0, &(temp.matrix_[h][0][0]), ncc);
        }
    }

    // F = T^ temp
    ta = 't';
    for (int h = 0; h < SO2AO->nirrep(); ++h) {
        m = nrow();
        n = ncol();
        k = temp.rowdim(h);
        nca = m; //k
        ncb = n; //k
        ncc = n; //k

        if (m && n && k) {
            C_DGEMM(ta, tb, m, n, k, 1.0, &(SO2AO->matrix_[h][0][0]),
                    nca, &(temp.matrix_[h][0][0]), ncb, 1.0, &(matrix_[0][0][0]),
                    ncc);
        }
    }
}

void Matrix::transform(const Matrix &transformer)
{
    Matrix temp(this);

    temp.gemm(false, false, 1.0, *this, transformer, 0.0);
    gemm(true, false, 1.0, transformer, temp, 0.0);
}

void Matrix::back_transform(const Matrix &a, const Matrix &transformer)
{
    Matrix temp(a);

    temp.gemm(false, true, 1.0, a, transformer, 0.0);
    gemm(false, false, 1.0, transformer, temp, 0.0);
}

void Matrix::back_transform(const Matrix &transformer)
{
    Matrix temp(this);

    temp.gemm(false, true, 1.0, *this, transformer, 0.0);
    gemm(false, false, 1.0, transformer, temp, 0.0);
}

double Matrix::vector_dot(const Matrix &rhs)
{
    return vector_dot(&rhs);
}

void Matrix::diagonalize(Matrix &eigvectors, Vector &eigvalues, int nMatz)
{
    int h;
    for (h = 0; h < nirrep_; ++h) {
        if (rowspi_[h]) {
            sq_rsp(rowspi_[h], colspi_[h], matrix_[h], eigvalues.vector_[h], nMatz, eigvectors.matrix_[h], 1.0e-14);
        }
    }
}

void Matrix::write_to_dpdfile2(dpdfile2 *outFile)
{
    global_dpd_->file2_mat_init(outFile);

    if (outFile->params->nirreps != nirrep_) {
        char *str = new char[100];
        sprintf(str, "Irrep count mismatch.  Matrix class has %d irreps, but dpdfile2 has %d.",
                nirrep_, outFile->params->nirreps);
        throw SanityCheckError(str, __FILE__, __LINE__);
    }

    if (outFile->my_irrep != 0) {
        throw FeatureNotImplemented("libmints Matrix class",
                                    "Matrices whose irrep is not the symmetric one", __FILE__, __LINE__);
    }

    for (int h = 0; h < nirrep_; ++h) {
        if (outFile->params->rowtot[h] != rowspi_[h]) {
            char *str = new char[100];
            sprintf(str, "Row count mismatch for irrep %d.  Matrix class has %d rows, but dpdfile2 has %d.",
                    h, rowspi_[h], outFile->params->rowtot[h]);
            throw SanityCheckError(str, __FILE__, __LINE__);
        }
        if (outFile->params->coltot[h] != colspi_[h]) {
            char *str = new char[100];
            sprintf(str, "Column count mismatch for irrep %d.  Matrix class has %d columns, but dpdfile2 has %d.",
                    h, colspi_[h], outFile->params->coltot[h]);
            throw SanityCheckError(str, __FILE__, __LINE__);
        }

        // TODO: optimize this with memcopys
        for (int row = 0; row < rowspi_[h]; ++row) {
            for (int col = 0; col < colspi_[h]; ++col) {
                outFile->matrix[h][row][col] = matrix_[h][row][col];
            }
        }
    }

    global_dpd_->file2_mat_wrt(outFile);
    global_dpd_->file2_mat_close(outFile);
}

void Matrix::save(const std::string &filename, bool append, bool saveLowerTriangle, bool saveSubBlocks)
{
    static const char *str_block_format = "%3d %3d %3d %16.16f\n";
    static const char *str_full_format = "%3d %3d %16.16f\n";

    // We can only save lower triangle if symmetry_ if 0
    if (symmetry_ && saveLowerTriangle)
        throw PSIEXCEPTION("Matrix::save: Unable to save lower triangle for non-totally symmetric matrix.");

    std::shared_ptr <psi::PsiOutStream> printer =
            std::shared_ptr<OutFile>(
                    new OutFile(filename, (append ? APPEND : TRUNCATE)));
    printer->Printf("%s\n", name_.c_str());
    printer->Printf("symmetry %d\n", symmetry_);

    if (saveSubBlocks == false) {
        // Convert the matrix to a full matrix
        double **fullblock = to_block_matrix();

        // Need to know the size
        int sizer = 0, sizec = 0;
        for (int h = 0; h < nirrep_; ++h) {
            sizer += rowspi_[h];
            sizec += colspi_[h ^ symmetry_];
        }

        if (saveLowerTriangle) {
            // Count the number of non-zero element
            int count = 0;
            for (int i = 0; i < sizer; ++i) {
                for (int j = 0; j <= i; ++j) {
                    if (fabs(fullblock[i][j]) > 1.0e-12) {
                        count++;
                    }
                }
            }
            printer->Printf("%5d\n", count);
            for (int i = 0; i < sizer; ++i) {
                for (int j = 0; j <= i; ++j) {
                    if (fabs(fullblock[i][j]) > 1.0e-12) {
                        printer->Printf(str_full_format, i, j, fullblock[i][j]);
                    }
                }
            }
        } else {
            // Count the number of non-zero element
            int count = 0;
            for (int i = 0; i < sizer; ++i) {
                for (int j = 0; j < sizec; ++j) {
                    if (fabs(fullblock[i][j]) > 1.0e-12) {
                        count++;
                    }
                }
            }
            printer->Printf("%5d\n", count);
            for (int i = 0; i < sizer; ++i) {
                for (int j = 0; j < sizec; ++j) {
                    if (fabs(fullblock[i][j]) > 1.0e-12) {
                        printer->Printf(str_full_format, i, j, fullblock[i][j]);
                    }
                }
            }
        }
        Matrix::free(fullblock);
    } else {
        if (saveLowerTriangle) {
            // Count the number of non-zero elements
            int count = 0;
            for (int h = 0; h < nirrep_; ++h) {
                for (int i = 0; i < rowspi_[h]; ++i) {
                    for (int j = 0; j <= i; ++j) {
                        if (fabs(matrix_[h][i][j]) > 1.0e-12) {
                            count++;
                        }
                    }
                }
            }
            printer->Printf("%5d\n", count);
            for (int h = 0; h < nirrep_; ++h) {
                for (int i = 0; i < rowspi_[h]; ++i) {
                    for (int j = 0; j <= i; ++j) {
                        if (fabs(matrix_[h][i][j]) > 1.0e-12) {
                            printer->Printf(str_block_format, h, i, j, matrix_[h][i][j]);
                        }
                    }
                }
            }
        } else {
            // Count the number of non-zero elements
            int count = 0;
            for (int h = 0; h < nirrep_; ++h) {
                for (int i = 0; i < rowspi_[h]; ++i) {
                    for (int j = 0; j < colspi_[h ^ symmetry_]; ++j) {
                        if (fabs(matrix_[h][i][j]) > 1.0e-12) {
                            count++;
                        }
                    }
                }
            }
            printer->Printf("%5d\n", count);
            for (int h = 0; h < nirrep_; ++h) {
                for (int i = 0; i < rowspi_[h]; ++i) {
                    for (int j = 0; j < colspi_[h ^ symmetry_]; ++j) {
                        if (fabs(matrix_[h][i][j]) > 1.0e-12) {
                            printer->Printf(str_block_format, h, i, j, matrix_[h][i][j]);
                        }
                    }
                }
            }
        }
    }
}

void Matrix::save(std::shared_ptr <psi::PSIO> &psio,
                  unsigned int fileno,
                  SaveType savetype)
{
    save(psio.get(), fileno, savetype);
}

void Matrix::save(psi::PSIO *const psio, unsigned int fileno, SaveType st)
{
    // Check to see if the file is open
    bool already_open = false;
    if (psio->open_check(fileno)) {
        already_open = true;
    } else {
        psio->open(fileno, PSIO_OPEN_OLD);
    }

    // Need to know the size
    int sizer = 0, sizec = 0;
    for (int h = 0; h < nirrep_; ++h) {
        sizer += rowspi_[h];
        sizec += colspi_[h ^ symmetry_];
    }

    if (st == SubBlocks) {
        for (int h = 0; h < nirrep_; ++h) {
            std::string str(name_);
            str += " Symmetry " + to_string(symmetry_) + " Irrep " + to_string(h);

            // Write the sub-blocks
            if (colspi_[h ^ symmetry_] > 0 && rowspi_[h] > 0)
                psio->write_entry(fileno, const_cast<char *>(str.c_str()), (char *) matrix_[h][0],
                                  sizeof(double) * colspi_[h ^ symmetry_] * rowspi_[h]);
        }
    } else if (st == Full) {
        double **fullblock = to_block_matrix();

        // Write the full block
        if (sizer > 0 && sizec > 0)
            psio->write_entry(fileno, const_cast<char *>(name_.c_str()), (char *) fullblock[0],
                              sizeof(double) * sizer * sizec);

        Matrix::free(fullblock);
    } else if (st == LowerTriangle) {
        double *lower = to_lower_triangle();

        if (sizer > 0)
            psio->write_entry(fileno, const_cast<char *>(name_.c_str()), (char *) lower, sizeof(double) * ioff[sizer]);
        delete[] lower;
    } else {
        throw PSIEXCEPTION("Matrix::save: Unknown SaveType\n");
    }
    if (!already_open)
        psio->close(fileno, 1);     // Close and keep
}

bool Matrix::load(std::shared_ptr <psi::PSIO> &psio,
                  unsigned int fileno,
                  const std::string &tocentry,
                  int nso)
{
    return load(psio.get(), fileno, tocentry, nso);
}

bool Matrix::load(psi::PSIO *const psio, unsigned int fileno, const std::string &tocentry, int nso)
{
    if (symmetry_) {
        throw PSIEXCEPTION("Matrix::load: Matrix is non-totally symmetric.");
    }

    double *integrals = init_array(ioff[nso]);

    // If psi fails to read in the data this will abort out.
    if (!tocentry.empty())
        psi::IWL::read_one(psio, fileno, tocentry.c_str(), integrals, ioff[nso], 0, 0, "outfile");
    else
        psi::IWL::read_one(psio, fileno, name_.c_str(), integrals, ioff[nso], 0, 0, "outfile");

    set(integrals);

    ::free(integrals);

    return true;
}

void Matrix::load(psi::PSIO *const psio, unsigned int fileno, SaveType st)
{
    // The matrix must be sized correctly first
    // Check to see if the file is open
    bool already_open = false;
    if (psio->open_check(fileno)) {
        already_open = true;
    } else {
        psio->open(fileno, PSIO_OPEN_OLD);
    }

    // Need to know the size
    int sizer = 0, sizec = 0;
    for (int h = 0; h < nirrep_; ++h) {
        sizer += rowspi_[h];
        sizec += colspi_[h ^ symmetry_];
    }

    if (st == SubBlocks) {
        for (int h = 0; h < nirrep_; ++h) {
            std::string str(name_);
            str += " Symmetry " + to_string(symmetry_) + " Irrep " + to_string(h);

            // Read the sub-blocks
            if (colspi_[h] > 0 && rowspi_[h] > 0)
                psio->read_entry(fileno, str.c_str(), (char *) matrix_[h][0],
                                 sizeof(double) * colspi_[h ^ symmetry_] * rowspi_[h]);
        }
    } else if (st == Full) {
        double **fullblock = to_block_matrix();

        // Read the full block
        if (sizer > 0 && sizec > 0)
            psio->read_entry(fileno, name_.c_str(), (char *) fullblock[0], sizeof(double) * sizer * sizec);

        set(fullblock);
        Matrix::free(fullblock);
    } else if (st == LowerTriangle) {
        double *lower = to_lower_triangle();

        if (sizer > 0)
            psio->read_entry(fileno, name_.c_str(), (char *) lower, sizeof(double) * ioff[sizer]);

        set(lower);
        delete[] lower;
    } else {
        throw PSIEXCEPTION("Matrix::load: Unknown SaveType\n");
    }
    if (!already_open)
        psio->close(fileno, 1);     // Close and keep // Idempotent, win!
}

void Matrix::load(std::shared_ptr <psi::PSIO> &psio,
                  unsigned int fileno,
                  SaveType savetype)
{
    load(psio.get(), fileno, savetype);
}

//void Matrix::alloc()
//{
//    if (matrix_)
//        release();

//    matrix_ = (double***)malloc(sizeof(double***) * nirrep_);
//    for (int h=0; h<nirrep_; ++h) {
//        if (rowspi_[h] != 0 && colspi_[h^symmetry_] != 0)
//            matrix_[h] = Matrix::matrix(rowspi_[h], colspi_[h^symmetry_]);
//        else {
//            // Force rowspi_[h] and colspi_[h^symmetry] to hard 0
//            // This solves an issue where a row can have 0 dim but a col does not (or the other way).

//            // This was commented out to resolve issues that people were
//            // dependent on one or both containing valid dimensions and
//            // not a hard zero.
//            //rowspi_[h] = colspi_[h^symmetry_] = 0;
//            matrix_[h] = NULL;
//        }
//    }
//}

void Matrix::load_mpqc(const std::string &filename)
{
    // Entire file.
    std::vector <std::string> lines;
    std::string text;

    std::ifstream infile(filename.c_str());
    if (!infile)
        throw PSIEXCEPTION("Matrix::load_mpqc: Unable to open file " + filename);

    // Stupid way of reading an entire file
    while (infile.good()) {
        getline(infile, text);
        trim_spaces(text);
        if (!text.empty())
            lines.push_back(text);
    }

    release();

    // First line is label
    set_name(lines[0]);

    // Second line is symmetry information
    std::smatch match;
    int infile_symm = -1;
    std::regex symmetry_line("^\\s*symmetry\\s*(\\d+)\\s*", std::regex_constants::icase);
    if (std::regex_match(lines[1], match, symmetry_line))
        infile_symm = str_to_int(match[1]);
    else
        throw PSIEXCEPTION("Matrix::load_mpqc: Second line must be 'symmetry #'");

    symmetry_ = infile_symm;

    // Third line is number of blocks (nirrep)
    int infile_nirrep = -1;
    std::regex blocks_line("^\\s*blocks\\s*(\\d+)\\s*", std::regex_constants::icase);
    if (std::regex_match(lines[2], match, blocks_line))
        infile_nirrep = str_to_int(match[1]);
    else
        throw PSIEXCEPTION("Matrix::load_mpqc: Third line must be 'blocks #'");

    nirrep_ = infile_nirrep;

    rowspi_ = Dimension(nirrep_);
    colspi_ = Dimension(nirrep_);

    // We'll handle memory allocation here
    matrix_ = (double ***) malloc(sizeof(double ***) * nirrep_);

    // This will hold the index in lines we should be working on
    size_t nline = 3;

    // Handle each block separately
    std::regex dim_line("^\\s*(\\d+)\\s*(\\d+)\\s*", std::regex_constants::icase);
    std::regex data_line("^\\s*\\d+\\s*" NUMBER "\\s*" NUMBER "?\\s*" NUMBER "?\\s*");

    for (int h = 0; h < nirrep_; ++h) {

        // Read dimension line
        if (std::regex_match(lines[nline], match, dim_line)) {
            rowspi_[h] = str_to_int(match[1]);
            colspi_[h ^ symmetry_] = str_to_int(match[2]);
        } else
            throw PSIEXCEPTION("Matrix::load_mpqc: Expected to find dimensions.");

        // Done with dimension line
        nline++;

        // Allocate memory
        if (rowspi_[h] != 0 && colspi_[h ^ symmetry_] != 0)
            matrix_[h] = Matrix::matrix(rowspi_[h], colspi_[h ^ symmetry_]);
        else
            matrix_[h] = NULL;

        // We read no more than 3 columns at a time
        for (int col = 0; col < colspi_[h ^ symmetry_]; col += 3) {
            // skip the next line (contains column numbers, which we ignore)
            nline++;

            for (int row = 0; row < rowspi_[h]; ++row) {
                if (regex_match(lines[nline], match, data_line)) {
                    std::string s1 = match[1];
                    std::string s2 = match[2];
                    std::string s3 = match[3];

//                    outfile->Printf( "'%s' %s (%d) %s (%d) %s (%d)\n",
//                            std::string(match[0]).c_str(),
//                            s1.c_str(), s1.length(),
//                            s2.c_str(), s2.length(),
//                            s3.c_str(), s3.length());

                    if (s1.length())
                        matrix_[h][row][col] = str_to_double(s1);
                    if (s2.length())
                        matrix_[h][row][col + 1] = str_to_double(s2);
                    if (s3.length())
                        matrix_[h][row][col + 2] = str_to_double(s3);
                } else
                    throw PSIEXCEPTION("Matrix::load_mpqc: Unable to match data line:\n" + lines[nline]);

                // Move on
                nline++;
            }
        }
        // Last thing to do.
        nline++;
    }
}

void Matrix::load(const std::string &filename)
{
    // Entire file.
    std::vector <std::string> lines;
    // temp
    std::string text;

    // Stream to use
    std::ifstream infile(filename.c_str());
    if (!infile)
        throw PSIEXCEPTION("Matrix::load: Unable to open file " + filename);

    // stupid way of reading in entire file
    while (infile.good()) {
        std::getline(infile, text);
        trim_spaces(text);
        if (!text.empty())
            lines.push_back(text);
    }

    // File MUST be at least 3 lines
    if (lines.size() < 3)
        throw PSIEXCEPTION("Matrix::load: " + filename + " insufficient length.");

    // First line is label (name of the matrix)
    set_name(lines[0]);

    // Second line is symmetry information
    std::smatch match;
    int infile_symm = -1;
    std::regex symmetry_line("^\\s*symmetry\\s*(\\d+)\\s*", std::regex_constants::icase);
    if (std::regex_match(lines[1], match, symmetry_line))
        infile_symm = str_to_int(match[1]);
    if (infile_symm != symmetry_) {
        release();
        symmetry_ = infile_symm;
        alloc();
    }

    // Third line is number of nonzero elements in the matrix.
    int infile_nonzero = 0;
    std::regex nonzero_line("^\\s*(\\d+)\\s*");
    if (std::regex_match(lines[2], match, nonzero_line))
        infile_nonzero = str_to_int(match[1]);

    // Check the number of lines with the number of nonzero elements
    int nline = lines.size();
    if (nline - 3 != infile_nonzero)
        throw PSIEXCEPTION("Matrix::load: Specified number of nonzero elements does not match number of elements given.");

    // Clear out this matrix
    zero();

    // Go through the file grabbing the data.
    std::regex element_line("^\\s*(\\d+)\\s*(\\d+)\\s*(\\d+)\\s*" NUMBER);
    int h, m, n;
    double val;
    for (int elem = 0; elem < infile_nonzero; ++elem) {
        if (std::regex_match(lines[elem + 3], match, element_line)) {
            h = str_to_int(match[1]);
            m = str_to_int(match[2]);
            n = str_to_int(match[3]);
            val = str_to_double(match[4]);
        } else
            throw PSIEXCEPTION("Matrix::load: Unable to match the following line:\n" + lines[elem + 3]);

        // Check the info to make sure it is good
        if (h >= nirrep_)
            throw PSIEXCEPTION("Matrix::load: Irrep number is too large:\n" + lines[elem + 3]);
        if (m >= rowspi_[h])
            throw PSIEXCEPTION("Matrix::load: Row number is too large:\n" + lines[elem + 3]);
        if (n >= colspi_[h ^ symmetry_])
            throw PSIEXCEPTION("Matrix::load: Column number is too large:\n" + lines[elem + 3]);

        // Set the data
        set(h, m, n, val);
    }
}

void Matrix::send()
{
}

void Matrix::recv()
{
}

void Matrix::bcast(int)
{
    std::cout << "Someone is calling the Matrix bcast routine..." << std::endl;
    // Assume the user allocated the matrix to the correct size first.
    /*for (int h=0; h<nirrep_; ++h) {
        if (rowspi_[h] > 0 && colspi_[h] > 0)
            WorldComm->bcast(matrix_[h][0], rowspi_[h] * colspi_[h^symmetry_], broadcaster);
    }*/
}

void Matrix::sum()
{
    //RMR--Removed the call to here because as
    //it stood it did nothing
}

bool Matrix::equal(const Matrix &rhs)
{
    return equal(&rhs);
}

bool Matrix::equal(const SharedMatrix &rhs)
{
    return equal(rhs.get());
}

bool Matrix::equal(const Matrix *rhs)
{
    // Check dimensions
    if (rhs->nirrep() != nirrep())
        return false;

    if (symmetry_ != rhs->symmetry_)
        return false;

    for (int h = 0; h < nirrep(); ++h)
        if ((rowspi()[h] != rhs->rowspi()[h]) ||
            (colspi()[h] != rhs->colspi()[h]))
            return false;

    // Check element by element
    for (int h = 0; h < nirrep(); ++h) {
        for (int m = 0; m < rowspi()[h]; ++m) {
            for (int n = 0; n < colspi()[h ^ symmetry_]; ++n) {
                if (get(h, m, n) != rhs->get(h, m, n))
                    return false;
            }
        }
    }

    return true;
}

bool Matrix::equal_but_for_row_order(const Matrix &rhs, double TOL)
{
    return equal_but_for_row_order(&rhs, TOL);
}

bool Matrix::equal_but_for_row_order(const SharedMatrix &rhs, double TOL)
{
    return equal_but_for_row_order(rhs.get(), TOL);
}

bool Matrix::equal_but_for_row_order(const Matrix *rhs, double TOL)
{
    if (rhs->nirrep() != nirrep())
        return false;

    if (symmetry_ != rhs->symmetry_)
        return false;

    for (int h = 0; h < nirrep(); ++h)
        if ((rowspi()[h] != rhs->rowspi()[h]) || (colspi()[h] != rhs->colspi()[h]))
            return false;

    for (int h = 0; h < nirrep(); ++h) {
        for (int m = 0; m < rowspi()[h]; ++m) {
            for (int m_rhs = 0; m_rhs < rowspi()[h]; ++m_rhs) {

                int n;
                for (n = 0; n < colspi()[h ^ symmetry_]; ++n) {
                    if (fabs(get(h, m, n) - rhs->get(h, m_rhs, n)) > TOL)
                        break;
                }

                if (n == colspi()[h ^ symmetry_]) {
                    break; // whole row matched, goto next m row
                }

                if (m_rhs == rowspi()[h] - 1)
                    return false; // no matching row was found
            }
        }
    }
    return true;
}

double Matrix::pyget(const py::tuple &key)
{
    return get(key[0].cast<int>(),
               key[1].cast<int>(),
               key[2].cast<int>());
}

void Matrix::pyset(const py::tuple &key, double value)
{
    return set(key[0].cast<int>(),
               key[1].cast<int>(),
               key[2].cast<int>(),
               value);
}

void Matrix::set_by_python_list(const py::list &data)
{
    size_t rows = py::len(data);

    // Make sure nrows < rows
    if ((size_t) nrow() > rows)
        throw PSIEXCEPTION("Uh, moron!");

    for (size_t i = 0; i < rows; ++i) {
        size_t cols = py::len(data[i]);
        if ((size_t) ncol() > cols)
            throw PSIEXCEPTION("Uh, moron!");
        for (size_t j = 0; j < cols; ++j) {
            set(i, j, data[i].cast<py::list>()[j].cast<double>());
        }
    }
}

void Matrix::rotate_columns(int h, int i, int j, double theta)
{
    if (h > nirrep_)
        throw PSIEXCEPTION("In rotate columns: Invalid Irrep");
    if (!colspi_[h] || !rowspi_[h]) return;
    if (i > colspi_[h])
        throw PSIEXCEPTION("In rotate columns: Invalid column number for i");
    if (j > colspi_[h])
        throw PSIEXCEPTION("In rotate columns: Invalid column number for j");
    double costheta = cos(theta);
    double sintheta = sin(theta);
    C_DROT(rowspi_[h], &matrix_[h][0][i], colspi_[h], &matrix_[h][0][j], colspi_[h], costheta, sintheta);
}
std::vector<py::buffer_info> Matrix::array_interface(){
    std::vector<py::buffer_info> ret;

    if (numpy_shape_.size()) {
        if (nirrep_ > 1){
            throw PSIEXCEPTION("Matrix::array_interface numpy shape with more than one irrep is not valid.");
        }

        std::vector<size_t> shape(numpy_shape_.size());
        std::vector<size_t> strides(numpy_shape_.size());
        size_t current_stride = sizeof(double);

        for (size_t i = numpy_shape_.size(); i-- > 0;) {
            shape[i] = numpy_shape_[i];
            strides[i] = current_stride;
            current_stride *= numpy_shape_[i];
        }
        ret.push_back(py::buffer_info(get_pointer(0), sizeof(double),
                                      py::format_descriptor<double>::format(),
                                      numpy_shape_.size(),
                                      shape, strides));

    } else {
        for (size_t h = 0; h < nirrep_; h++) {
            ret.push_back(py::buffer_info(get_pointer(h), sizeof(double),
                          py::format_descriptor<double>::format(), 2,
                          {static_cast<size_t>(rowspi(h)), static_cast<size_t>(colspi(h))},
                          {sizeof(double) * rowspi(h), sizeof(double)}));
        }
    }
    return ret;
}

} // namespace psi
