/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
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
 *@END LICENSE
 */

#ifndef _psi_src_lib_libmints_factory_h_
#define _psi_src_lib_libmints_factory_h_

#include <libmints/vector.h>
#include <libmints/matrix.h>
#include <libmints/dimension.h>

#include <libchkpt/chkpt.hpp>
#include <libpsio/psio.hpp>

namespace boost {
template<class T> class shared_ptr;
}

namespace psi {

class SOBasisSet;

/*! \ingroup MINTS
 *  \class MatrixFactory
 *  \brief A class for creating Matrix, SimpleMatrix, Vector, and SimpleVector objects.
 *
 * The objects this factory creates can automatically be sized based on information
 * from checkpoint.
 */
class MatrixFactory {
    /// Number of irreps
    int nirrep_;
    /// Number of orbitals
    int nso_;
    /// Number of rows per irrep
    int *rowspi_;
    /// Number of columns per irrep
    int *colspi_;

public:
    /// Default constructor, does nothing.
    MatrixFactory();
    /// Copy constructor.
    MatrixFactory(const MatrixFactory& copy);
    ~MatrixFactory();

    /// Initializes the matrix factory by creating a chkpt object with a psio reference.
    bool init_with_chkpt(boost::shared_ptr<psi::PSIO> psio);

    /// Initializes the matrix factory using the given chkpt object.
    bool init_with_chkpt(boost::shared_ptr<psi::Chkpt> chkpt);

    /// Manually initialize the matrix factory
    bool init_with(int nirrep, int *rowspi, int *colspi);

    /// Manually initialize the matrix factory with Dimension objects
    bool init_with(const Dimension& rows, const Dimension& cols);

    /// Manually initialize the matrix factory with SOBasisSet object
    bool init_with(const boost::shared_ptr<SOBasisSet>& sobasis);

    /// Returns number of irreps
    int nirrep() const;

    /// Returns the rows per irrep array
    int *rowspi() const;

    /// Returns the number of rows in irrep h
    int nrow(int h) const;

    /// Returns the columns per irrep array
    int *colspi() const;

    /// Returns the number of columns in irrep h
    int ncol(int h) const;

    /// Returns the number of orbitals
    int norb() const;

    /// Returns a new Matrix object with default dimensions
    Matrix * create_matrix(int symmetry=0);

    /// Returns a new Matrix object with default dimensions
    SharedMatrix create_shared_matrix();

    void create_matrix(Matrix& mat, int symmetry=0);

    /// Returns a new Matrix object named name with default dimensions
    Matrix * create_matrix(std::string name, int symmetry=0);

    SharedMatrix create_shared_matrix(const std::string& name);

    SharedMatrix create_shared_matrix(const std::string& name, int symmetry);

    SharedMatrix create_shared_matrix(const std::string& name, int rows, int cols);

    void create_matrix(Matrix& mat, std::string name, int symmetry=0);

    /// Returns a new Vector object with default dimensions
    Vector * create_vector();

    void create_vector(Vector& vec);

    SharedVector create_shared_vector(const std::string& name);
};

}

#endif
