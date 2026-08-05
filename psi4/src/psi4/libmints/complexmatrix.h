/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2026 The Psi4 Developers.
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

#ifndef _psi_src_lib_libmints_complexmatrix_h
#define _psi_src_lib_libmints_complexmatrix_h

#include "psi4/pragma.h"

#include <complex>
#include <memory>
#include <string>
#include <vector>
#include <array>

#include "dimension.h"

#ifdef USING_OpenOrbitalOptimizer
#ifdef USING_LAPACK_MKL
#include <mkl.h>
#define ARMA_USE_MKL
#define ARMA_USE_MKL_TYPES
#endif
#define ARMA_DONT_USE_FORTRAN_HIDDEN_ARGS
#define ARMA_DONT_USE_WRAPPER
#include <armadillo>
#else
namespace arma {
class cx_mat;
}
#endif

#ifdef USING_Einsums
#include <Einsums/Config.hpp>
#include <Einsums/Tensor.hpp>
#include <Einsums/LinearAlgebra.hpp>
#endif

namespace psi {

class PSIO;

// Certain friend free methods (e.g. `linalg::doublet`) need to access the private
// `tensor_` member to do lower-level Einsums operations. To be defined in a
// different namespace, you need to forward declare the new namespace and func.
// This requires a forward-declared ComplexMatrix class for argument types.
// This eliminates the need for adding a public getter for private members.
class ComplexMatrix;
using SharedComplexMatrix = std::shared_ptr<ComplexMatrix>;

#ifdef USING_Einsums

// Forward declare the friend functions.
namespace linalg {
template <bool, bool>
ComplexMatrix doublet(const ComplexMatrix&, const ComplexMatrix&);

template <bool, bool>
SharedComplexMatrix doublet(const SharedComplexMatrix&, const SharedComplexMatrix&);

template <bool, bool, bool>
ComplexMatrix triplet(const ComplexMatrix&, const ComplexMatrix&, const ComplexMatrix&);

template <bool, bool, bool>
SharedComplexMatrix triplet(const SharedComplexMatrix&, const SharedComplexMatrix&, const SharedComplexMatrix&);
}  // namespace linalg

/*! \ingroup MINTS
 *  \class ComplexMatrix
 *  \brief Complex blocked matrix backed by an einsums TiledTensor.
 *
 *  Wraps ``einmsums::TiledTensor<std::complex<double>, 2>`` as a private member
 *  and provides familiar API, similar to the (completely independent) ``Matrix``
 *  class. Tiles are analogous to irrep blocks in ``Matrix``, but TiledTensors
 *  allow tiles to be created on-demand, anywhere on the predefined grid.
 *  This grid is set at the constructor by rowspi and optionally colspi. Like
 *  ``Matrix``, the number of irreps in each dimension is assumed to be equal.
 *  Unlike ``Matrix``, we assume at the ctor: ``rowspi.n() == colspi.n()``.
 *
 *  TiledTensors work differently to ``Matrix``. Memory is handled by the
 *  object. When you call ``zero()`` the tiles are emptied, not actually set to
 *  zero. When you try to ``get()`` a tile, the tile is allocated and returns
 *  zeros. When you ``get()`` a specific value in an uninitialized tile, the
 *  tile may or may not be allocated, but a zero is always returned.
 *
 *  Adds operations needed by DIIS and the Python layer: clone, axpy, subtract,
 *  vector_dot, save, and load.
 */
class PSI_API ComplexMatrix {
   private:
    /// The backed einsums TiledTensor object.
    einsums::TiledTensor<std::complex<double>, 2> tensor_;

   public:
    using BlockT = einsums::Tensor<std::complex<double>, 2>;
    using TiledT = einsums::TiledTensor<std::complex<double>, 2>;
    using ValueT = std::complex<double>;

    // -- Default constructors (forward to TiledTensor) --

    ComplexMatrix() = default;
    ComplexMatrix(const ComplexMatrix&) = default;
    ComplexMatrix(ComplexMatrix&&) = default;
    ComplexMatrix& operator=(const ComplexMatrix&) = default;
    ComplexMatrix& operator=(ComplexMatrix&& other) = default;

    // -- Constructors for Dimension callers --

    ComplexMatrix(const std::string& name, const Dimension& row_sizes, const Dimension& col_sizes)
        : tensor_(name, row_sizes.blocks(), col_sizes.blocks()) {}

    ComplexMatrix(const std::string& name, const Dimension& row_sizes)
        : ComplexMatrix(name, row_sizes, row_sizes) {}

    ComplexMatrix(const Dimension& row_sizes, const Dimension& col_sizes)
        : ComplexMatrix("", row_sizes, col_sizes) {}

    ComplexMatrix(const Dimension& row_sizes)
        : ComplexMatrix("", row_sizes) {}

    /// Overload for single block of size
    ComplexMatrix(const std::string& name, int rows, int cols)
        : tensor_(name, std::vector<int>{rows}, std::vector<int>{cols}) {}

    ComplexMatrix(int rows, int cols)
        : ComplexMatrix("", rows, cols) {}

    // -- implicit conversion --

    /// Implicit conversion **from** TiledTensor
    ComplexMatrix(const TiledT& t) : tensor_(t) {}

    /// Implicit conversion **to** TiledTensor
    operator TiledT&() { return tensor_; }
    operator const TiledT&() const { return tensor_; }

    // -- ComplexMatrix-specific operations --

    /// Arithmetic operator ``+=``. Einsums asserts the dimensions at runtime.
    void operator+=(const ComplexMatrix& other) { tensor_ += other; }
    /// Arithmetic operator ``-=``. Einsums asserts the dimensions at runtime.
    void operator-=(const ComplexMatrix& other) { tensor_ -= other; }

    /// Deep copy.
    std::shared_ptr<ComplexMatrix> clone() const;

    /// In-place self += alpha * other (diagonal tiles only).
    void axpy(std::complex<double> alpha, const ComplexMatrix& other);

    /// In-place add
    void add(const ComplexMatrix& other) { tensor_ += other; }

    /// In-place subtract
    void subtract(const ComplexMatrix& other) { tensor_ -= other; }

    /// Re(Tr(self^H other)), summed over diagonal tiles.
    double vector_dot(const ComplexMatrix&) const;

    std::complex<double> product_trace(const ComplexMatrix&) const;

    /// Save diagonal tiles as raw complex sub-blocks to a PSIO file.
    void save(std::shared_ptr<PSIO>& psio, size_t fileno);

    /// Load diagonal tiles as raw complex sub-blocks from a PSIO file.
    void load(std::shared_ptr<PSIO>& psio, size_t fileno);

    /// Get the name, used for ``print()`` and Einsums runtime errors.
    const std::string& name() const { return tensor_.name(); }
    /// Set the name, used for ``print()`` and Einsums runtime errors.
    void set_name(const std::string& s) { tensor_.set_name(s); }

    /// Python-compatible printer.
    void print_out() const { print("outfile"); }
    /// Print to an ostream (delegates to the underlying TiledTensor).
    void print(std::string outfile = "outfile", const char* extra = nullptr) const;

    /// Number of possible row tiles (empty or filled)
    int nirrep() const { return tensor_.grid_size(0); }

    /// Returns dimension of row tile for given irrep.
    int rowdim(const int& h = 0) const { return tensor_.tile_size(0)[h]; }
    /// Returns dimension of column tile for given irrep.
    int coldim(const int& h = 0) const { return tensor_.tile_size(1)[h]; }

    // Unlike Matrix, ComplexMatrix does not have an internal rowspi_ function to return.
    // The lifetime of the reference in `Matrix::rowspi()` lasts as long as the
    // Matrix instance. Here, we have to create a Dimension then return it.
    // The Dimension lifetime is extended only by assigning it to a variable.
    // const Dimension would force a copy (not move) so Dimension is preferred.

    Dimension rowspi() const { return Dimension(tensor_.tile_size(0)); }
    int rowspi(const int& h) const { return rowdim(h); }
    Dimension colspi() const { return Dimension(tensor_.tile_size(1)); }
    int colspi(const int& h) const { return coldim(h); }

    /// Returns the total number of rows
    int nrow() const { return static_cast<int>(tensor_.dim(0)); }
    /// Returns the total number of columns
    int ncol() const { return static_cast<int>(tensor_.dim(1)); }

    PSI_DEPRECATED(
        "Internal einsums type used for constructors. Should be "
        "replaced with rowspi/colspi ASAP.")
    const std::array<std::vector<int>, 2> block_sizes() const { return tensor_.tile_sizes(); }

    /// Zero out all blocks. Deallocates tiles.
    void zero() { tensor_.zero(); }

    /// @{
    /// Returns the h'th irrep as ``einsums::Tensor<std::complex<double>, 2>`` type.
    BlockT& get(const int& h) { return tensor_.tile(h, h); }
    const BlockT& get(const int& h) const { return tensor_.tile(h, h); }
    BlockT& operator[](const int& h) { return tensor_.tile(h, h); }
    const BlockT& operator[](const int& h) const { return tensor_.tile(h, h); }
    BlockT& operator()(const int& h) { return tensor_.tile(h, h); }
    const BlockT& operator()(const int& h) const { return tensor_.tile(h, h); }
    /// @}

    /**
     * Returns a single element of tensor_
     *
     * @param h Subtile
     * @param m Row
     * @param n Column
     * @return value at position (h, m, n)
     */
    ValueT get(const int& h, const int& m, const int& n) const { return tensor_.tile(h, h)(m, n); }

    /// Returns a ``std::complex<double>&`` reference to element at position (h, i, j).
    ValueT& operator()(const int& h, const int& i, const int& j) { return tensor_.tile(h, h)(i, j); }
    /// Returns a ``const std::complex<double>&`` reference to element at position (h, i, j).
    const ValueT& operator()(const int& h, const int& i, const int& j) const { return tensor_.tile(h, h)(i, j); }

    /**
     * Set a single element of tensor_
     *
     * @param h Subtile
     * @param m Row
     * @param n Column
     * @param value Value
     */
    void set(const int& h, const int& m, const int& n, const ValueT& value) { tensor_.tile(h, h)(m, n) = value; }

#ifdef USING_OpenOrbitalOptimizer
    /// Returns an Armadillo complex matrix for irrep block ``h``
    arma::cx_mat to_armadillo_matrix(int h = 0);
    /// Copies data from an Armadillo complex matrix into irrep block ``h``
    void from_armadillo_matrix(const arma::cx_mat& m, int h = 0);
#endif

    template <bool, bool>
    friend ComplexMatrix linalg::doublet(const ComplexMatrix&, const ComplexMatrix&);

    template <bool, bool, bool>
    friend ComplexMatrix linalg::triplet(const ComplexMatrix&, const ComplexMatrix&, const ComplexMatrix&);

    /// Return a conjugate transpose ComplexMatrix
    SharedComplexMatrix conjT() const;
};

namespace linalg {

/**
 *  If this were placed in the source, no explicit specializations would be compiled,
 *  resulting in a linking error. Practically, this **is** a runtime error in Psi4.
 */

template <bool AdjoinA, bool AdjoinB>
ComplexMatrix doublet(const ComplexMatrix& A, const ComplexMatrix& B) {
    const Dimension C_rowspi = (AdjoinA) ? A.colspi() : A.rowspi();
    const Dimension C_colspi = (AdjoinB) ? B.rowspi() : B.colspi();

    ComplexMatrix C{"T", C_rowspi, C_colspi};

    einsums::linear_algebra::gemm<AdjoinA, AdjoinB>(std::complex<double>{1.0}, A.tensor_, B.tensor_,
                                                    std::complex<double>{0.0}, &C.tensor_);

    return C;
}

template <bool AdjoinA, bool AdjoinB>
SharedComplexMatrix doublet(const SharedComplexMatrix& A, const SharedComplexMatrix& B) {
    return std::make_shared<ComplexMatrix>(std::move(doublet<AdjoinA, AdjoinB>(*A, *B)));
}

template <bool AdjoinA, bool AdjoinB, bool AdjoinC>
ComplexMatrix triplet(const ComplexMatrix& A, const ComplexMatrix& B, const ComplexMatrix& C) {
    auto BC = doublet<AdjoinB, AdjoinC>(B, C);
    auto ABC = doublet<AdjoinA, AdjoinB>(A, BC);
    return ABC;
}

template <bool AdjoinA, bool AdjoinB, bool AdjoinC>
SharedComplexMatrix triplet(const SharedComplexMatrix& A, const SharedComplexMatrix& B, const SharedComplexMatrix& C) {
    return std::make_shared<ComplexMatrix>(std::move(triplet<AdjoinA, AdjoinB, AdjoinC>(*A, *B, *C)));
}

}  // namespace linalg

#else  // !USING_Einsums

/// Stub type so pybind can expose ComplexMatrix without Einsums.
class PSI_API ComplexMatrix {};
using SharedComplexMatrix = std::shared_ptr<ComplexMatrix>;

#endif  // USING_Einsums

}  // namespace psi

#endif
