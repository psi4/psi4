/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2016 The Psi4 Developers.
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

#include <boost/shared_ptr.hpp>
#include <exception.h>
#include <libqt/qt.h>

#include "matrix.h"
#include "view.h"

namespace psi {

//////////////////////////////////////////////////////////////////////////////
//
//  Class: View
//
//////////////////////////////////////////////////////////////////////////////

View::~View()
{
    nirrep_ = 0;
    delete[] row_offset_per_irrep_; row_offset_per_irrep_ = 0;
    delete[] col_offset_per_irrep_; col_offset_per_irrep_ = 0;
    delete[] rows_per_irrep_;       rows_per_irrep_ = 0;
    delete[] cols_per_irrep_;       cols_per_irrep_ = 0;
}

View::View(int nirrep, int *rows, int *cols)
    : nirrep_(nirrep), row_offset_per_irrep_(0), col_offset_per_irrep_(0),
      rows_per_irrep_(0), cols_per_irrep_(0)
{
    if (nirrep_ <= 0)
        throw InputException("Number of irreps is less than or equal to zero.", "nirrep", nirrep, __FILE__, __LINE__);
    if (rows == 0)
        throw InputException("Array of row sizes is 0.", "rows", 0, __FILE__, __LINE__);
    if (cols == 0)
        throw InputException("Array of column sizes is 0.", "cols", 0, __FILE__, __LINE__);

    rows_per_irrep_ = new int[nirrep_];
    cols_per_irrep_ = new int[nirrep_];
    row_offset_per_irrep_ = new int[nirrep_];
    col_offset_per_irrep_ = new int[nirrep_];

    for (int h=0; h<nirrep_; ++h) {
        rows_per_irrep_[h] = rows[h];
        cols_per_irrep_[h] = cols[h];
        row_offset_per_irrep_[h] = 0;
        col_offset_per_irrep_[h] = 0;
    }
}

View::View(int nirrep, int *rows, int *cols, int *row_offsets, int *col_offsets)
    : nirrep_(nirrep), row_offset_per_irrep_(0), col_offset_per_irrep_(0),
      rows_per_irrep_(0), cols_per_irrep_(0)
{
    if (nirrep_ <= 0)
        throw InputException("Number of irreps is less than or equal to zero.", "nirrep", nirrep, __FILE__, __LINE__);
    if (rows == 0)
        throw InputException("Array of row sizes is 0.", "rows", 0, __FILE__, __LINE__);
    if (cols == 0)
        throw InputException("Array of column sizes is 0.", "cols", 0, __FILE__, __LINE__);
    if (row_offsets == 0)
        throw InputException("Array of row offsets is 0.", "row_offsets", 0, __FILE__, __LINE__);
    if (col_offsets == 0)
        throw InputException("Array of column offsets is 0.", "col_offsets", 0, __FILE__, __LINE__);

    rows_per_irrep_ = new int[nirrep_];
    cols_per_irrep_ = new int[nirrep_];
    row_offset_per_irrep_ = new int[nirrep_];
    col_offset_per_irrep_ = new int[nirrep_];

    for (int h=0; h<nirrep_; ++h) {
        rows_per_irrep_[h] = rows[h];
        cols_per_irrep_[h] = cols[h];
        row_offset_per_irrep_[h] = row_offsets[h];
        col_offset_per_irrep_[h] = col_offsets[h];
    }
}

View::View(SharedMatrix matrix, const Dimension& rows, const Dimension& cols)
    : matrix_(matrix), nirrep_(0),
      row_offset_per_irrep_(0), col_offset_per_irrep_(0),
      rows_per_irrep_(0), cols_per_irrep_(0)

{
    nirrep_ = matrix_->nirrep();

    rows_per_irrep_ = new int[nirrep_];
    cols_per_irrep_ = new int[nirrep_];
    row_offset_per_irrep_ = new int[nirrep_];
    col_offset_per_irrep_ = new int[nirrep_];

    for (int h=0; h<nirrep_; ++h) {
        rows_per_irrep_[h] = rows[h];
        cols_per_irrep_[h] = cols[h];
        row_offset_per_irrep_[h] = 0;
        col_offset_per_irrep_[h] = 0;
    }
}

View::View(SharedMatrix matrix, const Dimension& rows, const Dimension& cols, const Dimension& row_offsets, const Dimension& col_offsets)
    : matrix_(matrix), nirrep_(0),
      row_offset_per_irrep_(0), col_offset_per_irrep_(0),
      rows_per_irrep_(0), cols_per_irrep_(0)

{
    nirrep_ = matrix_->nirrep();

    rows_per_irrep_ = new int[nirrep_];
    cols_per_irrep_ = new int[nirrep_];
    row_offset_per_irrep_ = new int[nirrep_];
    col_offset_per_irrep_ = new int[nirrep_];

    for (int h=0; h<nirrep_; ++h) {
        rows_per_irrep_[h] = rows[h];
        cols_per_irrep_[h] = cols[h];
        row_offset_per_irrep_[h] = row_offsets[h];
        col_offset_per_irrep_[h] = col_offsets[h];
    }
}

SharedMatrix View::operator ()()
{
    // Create a new matrix with needed size
    SharedMatrix matrix = SharedMatrix(new Matrix(nirrep_, rows_per_irrep_, cols_per_irrep_));

    // Copy over the data
    for (int h=0; h<nirrep_; ++h) {
        for (int i=0; i<rows_per_irrep_[h]; ++i) {
            int ioffset = i + row_offset_per_irrep_[h];
            for (int j=0; j<cols_per_irrep_[h]; ++j) {
                int joffset = j + col_offset_per_irrep_[h];

                matrix->set(h, i, j, matrix_->get(h, ioffset, joffset));
            }
        }
    }

    return matrix;
}

SharedMatrix View::view(SharedMatrix matrix)
{
    SharedMatrix old = matrix_;
    matrix_ = matrix;
    return old;
}

}