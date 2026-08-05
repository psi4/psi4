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

#include "psi4/libmints/complexmatrix.h"

#include "psi4/libpsio/psio.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libpsi4util/exception.h"
#include "psi4/libmints/dimension.h"

// For outfile in ComplexMatrix::print()
#include "psi4/psi4-dec.h"
#include "psi4/libpsi4util/PsiOutStream.h"

#ifdef USING_Einsums
#include <Einsums/TensorAlgebra.hpp>
#include <Einsums/LinearAlgebra.hpp>
#endif

namespace psi {

#ifdef USING_Einsums

std::shared_ptr<ComplexMatrix> ComplexMatrix::clone() const { return std::make_shared<ComplexMatrix>(*this); }

#ifdef USING_OpenOrbitalOptimizer
arma::cx_mat ComplexMatrix::to_armadillo_matrix(int h) {
    int nc = coldim(h);
    int nr = rowdim(h);
    arma::cx_mat m(nr, nc);
    for (int ir = 0; ir < nr; ir++)
        for (int ic = 0; ic < nc; ic++) m(ir, ic) = get(h, ir, ic);
    return m;
}

void ComplexMatrix::from_armadillo_matrix(const arma::cx_mat& m, int h) {
    int nc = coldim(h);
    int nr = rowdim(h);
    for (int ir = 0; ir < nr; ir++)
        for (int ic = 0; ic < nc; ic++) set(h, ir, ic, m(ir, ic));
}
#endif

// self += alpha * other
//
// Implemented as a plain element loop (via operator(p, q), the same accessor used
// throughout cghf.cc/ComplexJK) rather than einsums::linear_algebra::axpy, whose
// generic TensorConcept dispatch appears to not work.
void ComplexMatrix::axpy(std::complex<double> alpha, const ComplexMatrix& other) {
    const auto& other_t = static_cast<const TiledT&>(other);
    if (tensor_.grid_size(0) != other_t.grid_size(0) || tensor_.grid_size(1) != other_t.grid_size(1)) {
        throw PSIEXCEPTION("ComplexMatrix::axpy: tile grids must match.");
    }
    for (int h = 0; h < static_cast<int>(other_t.grid_size(0)); ++h) {
        if (!other_t.has_tile(h, h) || other_t.has_zero_size(h, h)) continue;
        const auto& B = other_t.tile(h, h);
        auto& A = tensor_.tile(h, h);  // lazily allocated if missing
        const int nr = static_cast<int>(B.dim(0));
        const int nc = static_cast<int>(B.dim(1));
        for (int p = 0; p < nr; ++p) {
            for (int q = 0; q < nc; ++q) {
                A(p, q) += alpha * B(p, q);
            }
        }
    }
}

// Re(Tr(self^H other)), summed over diagonal tiles
double ComplexMatrix::vector_dot(const ComplexMatrix& other) const {
    const auto& other_t = other.tensor_;
    std::complex<double> total{0.0, 0.0};
    for (int h = 0; h < static_cast<int>(tensor_.grid_size(0)); ++h) {
        if (!tensor_.has_tile(h, h) || !other_t.has_tile(h, h)) continue;
        const auto& A = tensor_.tile(h, h);
        const auto& B = other_t.tile(h, h);
        const int nr = static_cast<int>(A.dim(0));
        const int nc = static_cast<int>(A.dim(1));
        for (int p = 0; p < nr; ++p) {
            for (int q = 0; q < nc; ++q) {
                total += std::conj(A(p, q)) * B(p, q);
            }
        }
    }
    return total.real();
}

// np.einsum("ij,ji->", this, other)
std::complex<double> ComplexMatrix::product_trace(const ComplexMatrix& other) const {
    std::complex<double> E;

    using namespace einsums;
    tensor_algebra::einsum(Indices{}, &E, Indices{index::i, index::j}, tensor_, Indices{index::j, index::i},
                           other.tensor_);

    return E;
}

// Raw per-tile complex sub-blocks to a PSIO file, mirroring Matrix::save with
// SaveType::SubBlocks (libmints/matrix.cc).
void ComplexMatrix::save(std::shared_ptr<PSIO>& psio, size_t fileno) {
    bool already_open = psio->open_check(fileno);
    if (!already_open) psio->open(fileno, PSIO_OPEN_OLD);

    for (int h = 0; h < static_cast<int>(tensor_.grid_size(0)); ++h) {
        if (!tensor_.has_tile(h, h) || tensor_.has_zero_size(h, h)) continue;
        auto& t = tensor_.tile(h, h);
        std::string entry = tensor_.name() + " Tile " + std::to_string(h);
        psio->write_entry(fileno, entry, (char*)t.data(), sizeof(std::complex<double>) * t.size());
    }

    if (!already_open) psio->close(fileno, 1);  // keep
}

// The ComplexMatrix must already have the right tile grid before loading (as with
// Matrix::load), e.g. constructed via ComplexMatrix(name, block_sizes).
void ComplexMatrix::load(std::shared_ptr<PSIO>& psio, size_t fileno) {
    bool already_open = psio->open_check(fileno);
    if (!already_open) psio->open(fileno, PSIO_OPEN_OLD);

    for (int h = 0; h < static_cast<int>(tensor_.grid_size(0)); ++h) {
        if (tensor_.tile_size(0)[h] == 0 || tensor_.tile_size(1)[h] == 0) continue;
        auto& t = tensor_.tile(h, h);  // lazily allocated to the declared size
        std::string entry = tensor_.name() + " Tile " + std::to_string(h);
        psio->read_entry(fileno, entry, (char*)t.data(), sizeof(std::complex<double>) * t.size());
    }

    if (!already_open) psio->close(fileno, 1);
}

void ComplexMatrix::print(std::string out, const char* extra) const {
    if (extra != nullptr) {
        throw PSIEXCEPTION("Not implemented");
    }

    std::shared_ptr<psi::PsiOutStream> printer = (out == "outfile" ? outfile : std::make_shared<PsiOutStream>(out));

    einsums::fprintln(*printer->stream(), tensor_);
}

SharedComplexMatrix ComplexMatrix::conjT() const {
    auto temp = std::make_shared<ComplexMatrix>(this->name() + "^H", this->colspi(), this->rowspi());

    // using keyT = const std::array<int, 2>;
    // using valueT = const einsums::Tensor<std::complex<double>, 2>;
    for (const auto& [key, value] : tensor_.tiles()) {
        const auto& [rowpi, colpi] = key;
        static_assert(std::is_same_v<decltype(rowpi), const int>);

        // necessarily: tensor_.has_tile(rowpi, colpi) == true
        // This adds a tile to temp.
        auto& ttile = temp->tensor_.tile(colpi, rowpi);

        // The template parameter <true> means take the conjugate.
#pragma message("Einsums::tensor_algebra::permute API will change as soon as v2.0")
        einsums::tensor_algebra::permute<true>(
            std::complex<double>{0.0}, einsums::Indices{einsums::index::i, einsums::index::j}, &ttile,
            std::complex<double>{1.0}, einsums::Indices{einsums::index::j, einsums::index::i}, value);
    }

    return temp;
}

#endif  // USING_Einsums

}  // namespace psi
