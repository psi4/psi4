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

#ifndef libmints_view_h_
#define libmints_view_h_

#include "typedefs.h"

namespace psi {

struct dpdfile2;

class PSIO;
class Matrix;
class Vector;
class SimpleVector;
class Dimension;

// Forward declarations
////////////////////////////////////////////


/*! \ingroup MINTS
 *  \class View
 *  \brief Provides a view to a region of the matrix.
 */
class View
{
protected:
    /// Matrix we are viewing.
    SharedMatrix matrix_;
    /// Number of irreps
    int nirrep_;
    /// Starting offsets in matrix_;
    int *row_offset_per_irrep_;
    int *col_offset_per_irrep_;
    /// Number of rows we are viewing
    int *rows_per_irrep_;
    int *cols_per_irrep_;

private:
    View();  // No default constructor
    View(const View& );  // No copy constructor

public:
    virtual ~View();

    /** Constructor, assumes offsets for each irrep is 0
     *  @param nirrep Number of irreps
     *  @param rows How many rows per irrep are we interested in
     *  @param cols How many cols per irrep are we interested in
     */
    View(int nirrep, int *rows, int *cols);
    /** Constructor, user provides offsets and dimensions.
     *  @param nirrep Number of irreps
     *  @param rows How many rows per irrep
     *  @param cols How many cols per irrep
     *  @param row_offsets Row offset per irrep
     *  @param col_offsets Column offset per irrep
     */
    View(int nirrep, int *rows, int *cols, int *row_offsets, int *col_offsets);
    /** Constructor, user provides a Matrix to view and desired row count
     *  @param matrix Matrix we want to view, View obtains nirrep from it
     *  @param rows How many rows per irrep
     *  @param cols How many cols per irrep
     */
    View(SharedMatrix matrix, const Dimension& rows, const Dimension& cols);
    /** Constructor, user provides a Matrix to view and desired row count
     *  @param matrix Matrix we want to view, View obtains nirrep from it
     *  @param rows How many rows per irrep
     *  @param cols How many cols per irrep
     *  @param row_offsets Row offset per irrep
     *  @param col_offsets Column offset per irrep
     */
    View(SharedMatrix matrix,
         const Dimension& rows, const Dimension& cols,
         const Dimension& row_offsets, const Dimension& col_offsets);

    /** Operator () overload. Creates a new Matrix that only contains the view.
     *  @return New Matrix containing the view.
     */
    SharedMatrix operator()();

    /** Set the Matrix that we should be viewing.
     *  @param matrix Matrix to view.
     *  @return The old Matrix we were viewing.
     */
    SharedMatrix view(SharedMatrix matrix);
};

}

#endif
