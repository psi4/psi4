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

#ifndef _psi_src_lib_libmints_factory_h_
#define _psi_src_lib_libmints_factory_h_

#include "psi4/libmints/vector.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/dimension.h"
#include "psi4/libpsio/psio.hpp"

namespace psi {

class SOBasisSet;

/*! \ingroup MINTS
 *  \class MatrixFactory
 *  \brief A class for creating Matrix and Vector objects.
 *
 * The objects this factory creates can automatically be sized based on information
 * from checkpoint.
 */
class PSI_API MatrixFactory {
    /// Number of irreps
    int nirrep_;
    /// Number of orbitals
    int nso_;
    /// Number of rows per irrep
    Dimension rowspi_;
    /// Number of columns per irrep
    Dimension colspi_;

   public:
    /// Default constructor, does nothing.
    MatrixFactory();
    /// Copy constructor.
    MatrixFactory(const MatrixFactory& copy);
    ~MatrixFactory();

    /// Manually initialize the matrix factory
    bool init_with(int nirrep, int* rowspi, int* colspi);

    /// Manually initialize the matrix factory with Dimension objects
    bool init_with(const Dimension& rows, const Dimension& cols);

    /// Manually initialize the matrix factory with SOBasisSet object
    bool init_with(const std::shared_ptr<SOBasisSet>& sobasis);

    /// Returns number of irreps
    int nirrep() const;

    /// Returns the rows per irrep array
    const Dimension& rowspi() const;

    /// Returns the number of rows in irrep h
    int nrow(int h) const;

    /// Returns the columns per irrep array
    const Dimension& colspi() const;

    /// Returns the number of columns in irrep h
    int ncol(int h) const;

    /// Returns the number of orbitals
    int norb() const;

    /// Returns a new Matrix object with default dimensions
    std::unique_ptr<Matrix> create_matrix(int symmetry = 0);

    /// Returns a new Matrix object with default dimensions
    SharedMatrix create_shared_matrix() const;

    void create_matrix(Matrix& mat, int symmetry = 0);

    /// Returns a new Matrix object named name with default dimensions
    std::unique_ptr<Matrix> create_matrix(std::string name, int symmetry = 0);

    SharedMatrix create_shared_matrix(const std::string& name) const;

    SharedMatrix create_shared_matrix(const std::string& name, int symmetry) const;

    SharedMatrix create_shared_matrix(const std::string& name, int rows, int cols) const;

    void create_matrix(Matrix& mat, std::string name, int symmetry = 0);

    /// Returns a new Vector object with default dimensions
    std::unique_ptr<Vector> create_vector();

    void create_vector(Vector& vec);

    SharedVector create_shared_vector(const std::string& name);
};

}  // namespace psi

#endif
