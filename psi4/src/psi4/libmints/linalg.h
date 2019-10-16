/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2019 The Psi4 Developers.
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

#pragma once

#include <complex>
#include <tuple>
#include <type_traits>
#include <vector>

#include <xtensor-blas/xlinalg.hpp>

#include "psi4/libpsi4util/exception.h"

#include "dimension.h"
#include "tensor.h"

namespace psi {
enum class Operation { None, Transpose, TransposeConj };

namespace detail {
/*! Whether the operation involves a transposition */
inline bool do_transpose(Operation op) {
    return ((op == Operation::Transpose || op == Operation::TransposeConj) ? true : false);
}

/*! Whether the operation involves a conjugation */
inline bool do_conjugate(Operation op) { return ((op == Operation::TransposeConj) ? true : false); }

/*! Type trait for mixable types
 *
 * These are used for SFINAE in 2-arity functions.
 * Mixable types are:
 *   - T = U (obviously)
 *   - T = float and U = double
 *   - T = double and  U = complex<double>
 * and permutations thereof.
 */
template <typename T, typename U>
struct is_2arity_mixable
    : std::integral_constant<bool, !(std::is_same<T, float>::value && std::is_same<U, std::complex<double>>::value)> {};

// NOTE for posterity this can be declared inline in C++17
template <typename T, typename U>
constexpr bool is_2arity_mixable_v = is_2arity_mixable<T, U>::value;
}  // namespace detail

/*! Return a tensor with all blocks filled with given value of same shape and value type as input
 *  \param[in] mold input tensor
 *  \param[in] fill_value the value in all blocks
 */
template <typename T, size_t Rank, typename U = T>
SharedTensor<U, Rank> full_like(const SharedTensor<T, Rank>& mold, U fill_value) noexcept {
    return std::make_shared<Tensor<U, Rank>>(mold->label(), mold->nirrep(), mold->axes_dimpi(), mold->symmetry(),
                                             static_cast<U>(fill_value));
}

/*! Return a tensor with all blocks filled with 0 of same shape and value type as input
 *  \param[in] mold input tensor
 */
template <typename T, size_t Rank, typename U = T>
SharedTensor<U, Rank> zeros_like(const SharedTensor<T, Rank>& mold) noexcept {
    return full_like(mold, static_cast<U>(0));
}

/*! Return a tensor with all blocks filled with 1 of same shape and value type as input
 *  \param[in] mold input tensor
 */
template <typename T, size_t Rank, typename U = T>
SharedTensor<U, Rank> ones_like(const SharedTensor<T, Rank>& mold) noexcept {
    return full_like(mold, static_cast<U>(1));
}

/*! Symmetry-blocking aware GEneralized Matrix Vector multiplication (GEMV)
 *  \param[in] opA preliminary operation on A (None, Transpose, TransposeConj)
 *  \param[in] alpha scaling of A*v
 *  \param[in] A matrix
 *  \param[in] x vector
 *  \param[in] beta scaling of result
 *  \param[in] y result vector, that is y := alpha * op(A) * x + beta * y
 */
template <typename T, typename U = T, typename V = decltype(std::declval<T>() * std::declval<U>()),
          typename = std::enable_if_t<detail::is_2arity_mixable_v<T, U>>>
SharedTensor<V, 1> gemv(Operation opA, double alpha, const SharedTensor<T, 2>& A, const SharedTensor<U, 1>& x,
                        double beta, SharedTensor<V, 1>& y) {
    std::string label = "GEMV: ";

    std::vector<size_t> ax_A;
    bool transA = detail::do_transpose(opA);
    bool conjA = detail::do_conjugate(opA);
    switch (opA) {
        case Operation::None:
            ax_A.push_back(1);
            label += A->label();
            break;
        case Operation::Transpose:
            ax_A.push_back(0);
            label += A->label() + "T x ";
            break;
        case Operation::TransposeConj:
            ax_A.push_back(0);
            label += A->label() + "H x ";
            break;
    }

    for (size_t hA = 0; hA < A->nirrep(); ++hA) {
        size_t hx = hA ^ (transA ? 0 : A->symmetry());
        size_t hy = hA ^ (transA ? A->symmetry() : 0);

        if (!conjA) {  // No conjugation
            y->block(hy) = xt::linalg::tensordot(alpha * A->block(hA), x->block(hx), ax_A, {0});
        } else {  // Conjugate A
            y->block(hy) = xt::linalg::tensordot(alpha * xt::conj(A->block(hA)), x->block(hx), ax_A, {0});
        }
        y->block(hy) += beta * y->block(hy);
    }
    // More descriptive label
    y->set_label(label);

    return y;
}

/*! Symmetry-blocking aware GEneralized Matrix Vector multiplication (GEMV), with internal allocation
 *  \param[in] opA preliminary operation on A (None, Transpose, TransposeConj)
 *  \param[in] alpha scaling of A*v
 *  \param[in] A matrix
 *  \param[in] x vector
 *  \param[in] beta scaling of result
 *  \param[in] y result vector, that is y := alpha * op(A) * x + beta * y
 */
template <typename T, typename U = T, typename V = decltype(std::declval<T>() * std::declval<U>()),
          typename = std::enable_if_t<detail::is_2arity_mixable_v<T, U>>>
SharedTensor<V, 1> gemv(Operation opA, double alpha, const SharedTensor<T, 2>& A, const SharedTensor<U, 1>& x,
                        double beta) {
    auto y = zeros_like<U, 1, V>(x);
    gemv(opA, alpha, A, x, beta, y);
    return y;
}

/*! Symmetry-blocking aware GEneralized Matrix Multiplication (GEMM)
 *  \param[in] opA preliminary operation on A (None, Transpose, TransposeConj)
 *  \param[in] opB preliminary operation on B (None, Transpose, TransposeConj)
 *  \param[in] alpha scaling of A*B
 *  \param[in] A first matrix
 *  \param[in] B second matrix
 *  \param[in] beta scaling of result
 *  \param[in] C result matrix, that is C := alpha * op(A) * op(B) + beta * C
 */
template <typename T, typename U = T, typename V = decltype(std::declval<T>() * std::declval<U>()),
          typename = std::enable_if_t<detail::is_2arity_mixable_v<T, U>>>
SharedTensor<V, 2> gemm(Operation opA, Operation opB, double alpha, const SharedTensor<T, 2>& A,
                        const SharedTensor<U, 2>& B, double beta, SharedTensor<V, 2>& C) {
    std::string label = "GEMM: ";

    std::vector<size_t> ax_A;
    bool transA = detail::do_transpose(opA);
    bool conjA = detail::do_conjugate(opA);
    switch (opA) {
        case Operation::None:
            ax_A.push_back(1);
            label += A->label();
            break;
        case Operation::Transpose:
            ax_A.push_back(0);
            label += A->label() + "T x ";
            break;
        case Operation::TransposeConj:
            ax_A.push_back(0);
            label += A->label() + "H x ";
            break;
    }

    std::vector<size_t> ax_B;
    bool transB = detail::do_transpose(opB);
    bool conjB = detail::do_conjugate(opB);
    switch (opB) {
        case Operation::None:
            ax_B.push_back(0);
            label += B->label();
            break;
        case Operation::Transpose:
            ax_B.push_back(1);
            label += B->label() + "T";
            break;
        case Operation::TransposeConj:
            ax_B.push_back(1);
            label += B->label() + "H";
            break;
    }

    for (size_t hA = 0; hA < A->nirrep(); ++hA) {
        size_t hB = hA ^ (transA ? 0 : A->symmetry()) ^ (transB ? B->symmetry() : 0);
        size_t hC = hA ^ (transA ? A->symmetry() : 0);

        if (!conjA && !conjB) {  // No conjugation
            C->block(hC) = xt::linalg::tensordot(alpha * A->block(hA), B->block(hB), ax_A, ax_B);
        } else if (!conjA && conjB) {  // Conjugate B
            C->block(hC) = xt::linalg::tensordot(alpha * A->block(hA), xt::conj(B->block(hB)), ax_A, ax_B);
        } else if (conjA && !conjB) {  // Conjugate A
            C->block(hC) = xt::linalg::tensordot(alpha * xt::conj(A->block(hA)), B->block(hB), ax_A, ax_B);
        } else {  // Conjugate A and B
            C->block(hC) = xt::linalg::tensordot(alpha * xt::conj(A->block(hA)), xt::conj(B->block(hB)), ax_A, ax_B);
        }
        C->block(hC) += beta * C->block(hC);
    }
    // More descriptive label
    C->set_label(label);

    return C;
}

/*! Symmetry-blocking aware GEneralized Matrix Multiplication (GEMM), with internal allocation
 *  \param[in] opA preliminary operation on A (None, Transpose, TransposeConj)
 *  \param[in] opB preliminary operation on B (None, Transpose, TransposeConj)
 *  \param[in] alpha scaling of A*B
 *  \param[in] A first matrix
 *  \param[in] B second matrix
 *  \param[in] beta scaling of result
 *  \param[in] C result matrix, that is C := alpha * op(A) * op(B) + beta * C
 */
template <typename T, typename U = T, typename V = decltype(std::declval<T>() * std::declval<U>()),
          typename = std::enable_if_t<detail::is_2arity_mixable_v<T, U>>>
SharedTensor<V, 2> gemm(Operation opA, Operation opB, double alpha, const SharedTensor<T, 2>& A,
                        const SharedTensor<U, 2>& B, double beta) {
    auto C = zeros_like<T, 2, V>(A);
    gemm(opA, opB, alpha, A, B, beta, C);
    return C;
}

/*! Simple doublet GEMM with on-the-fly allocation
 * \param[in] A The first matrix
 * \param[in] B The second matrix
 * \param[in] opA preliminary operation on A (None, Transpose, TransposeConj)
 * \param[in] opB preliminary operation on B (None, Transpose, TransposeConj)
 */
template <typename T, typename U = T, typename V = decltype(std::declval<T>() * std::declval<U>()),
          typename = std::enable_if_t<detail::is_2arity_mixable_v<T, U>>>
SharedTensor<V, 2> doublet(const SharedTensor<T, 2>& A, const SharedTensor<U, 2>& B, Operation opA = Operation::None,
                           Operation opB = Operation::None) {
    Dimension rowspi = (opA == Operation::Transpose || opA == Operation::TransposeConj) ? A->colspi() : A->rowspi();
    Dimension colspi = (opB == Operation::Transpose || opB == Operation::TransposeConj) ? B->rowspi() : B->colspi();

    auto C = std::make_shared<Tensor<V, 2>>("result", rowspi, colspi, A->symmetry() ^ B->symmetry());

    return gemm(opA, opB, 1.0, A, B, 0.0, C);
}

/*! Simple doublet GEMM with on-the-fly allocation
 * \param[in] A The first matrix
 * \param[in] B The second matrix
 * \param[in] transA Transpose the first matrix
 * \param[in] transB Transpose the second matrix
 */
template <typename T, typename U = T, typename V = decltype(std::declval<T>() * std::declval<U>()),
          typename = std::enable_if_t<detail::is_2arity_mixable_v<T, U>>>
SharedTensor<V, 2> doublet(const SharedTensor<T, 2>& A, const SharedTensor<U, 2>& B, bool transA = false,
                           bool transB = false) {
    Operation opA = transA ? Operation::Transpose : Operation::None;
    Operation opB = transB ? Operation::Transpose : Operation::None;

    return doublet(A, B, opA, opB);
}

template <typename T, size_t Rank, typename U = typename detail::is_complex<T>::real_type>
SharedTensor<U, Rank> real(const SharedTensor<T, Rank>& in) noexcept {
    auto out = zeros_like<T, Rank, U>(in);
    std::transform(in->cbegin(), in->cend(), out->begin(),
                   [](const typename Tensor<T, Rank>::block_type& blk) ->
                   typename Tensor<U, Rank>::block_type { return xt::real(blk); });
    return out;
}

template <typename T, size_t Rank, typename U = typename detail::is_complex<T>::real_type>
SharedTensor<U, Rank> imag(const SharedTensor<T, Rank>& in) noexcept {
    auto out = zeros_like<T, Rank, U>(in);
    std::transform(in->cbegin(), in->cend(), out->begin(),
                   [](const typename Tensor<T, Rank>::block_type& blk) ->
                   typename Tensor<U, Rank>::block_type { return xt::imag(blk); });
    return out;
}

template <typename T, size_t Rank>
SharedTensor<T, Rank> conj(const SharedTensor<T, Rank>& in) noexcept {
    auto out = zeros_like(in);
    std::transform(in->cbegin(), in->cend(), out->begin(),
                   [](const typename Tensor<T, Rank>::block_type& blk) ->
                   typename Tensor<T, Rank>::block_type { return xt::conj(blk); });
    return out;
}

/*! @{ Decompositions */
template <typename T>
SharedTensor<T, 2> cholesky(const SharedTensor<T, 2>& in) {
    if (in->symmetry()) {
        throw PSIEXCEPTION("Matrix is non-totally symmetric.");
    }
    auto out = zeros_like(in);
    std::transform(in->cbegin(), in->cend(), out->begin(),
                   [](const typename Tensor<T, 2>::block_type& blk) ->
                   typename Tensor<T, 2>::block_type { return xt::linalg::cholesky(blk); });
    return out;
}

template <typename T>
SharedTensor<T, 2> partial_cholesky(const SharedTensor<T, 2>& in) {
    if (in->symmetry()) {
        throw PSIEXCEPTION("Matrix is non-totally symmetric.");
    }
    auto out = zeros_like(in);
    std::transform(in->cbegin(), in->cend(), out->begin(),
                   [](const typename Tensor<T, 2>::block_type& blk) ->
                   typename Tensor<T, 2>::block_type { return xt::linalg::cholesky(blk); });
    return out;
}

template <typename T>
using QRResult = std::tuple<SharedTensor<T, 2>, SharedTensor<T, 2>>;

template <typename T>
auto qr(const SharedTensor<T, 2>& in, xt::linalg::qrmode mode = xt::linalg::qrmode::reduced) -> QRResult<T> {
    using Block = typename Tensor<T, 2>::block_type;
    using QRs = std::tuple<Block, Block>;
    auto tmp = std::vector<QRs>(in->nirrep());

    std::transform(in->cbegin(), in->cend(), tmp.begin(),
                   [mode](const Block& blk) -> QRs { return xt::linalg::qr(blk, mode); });

    auto Q = std::make_shared<Tensor<T, 2>>("Q matrix of " + in->label(), in->rowspi(), in->colspi());
    std::transform(tmp.cbegin(), tmp.cend(), Q->begin(), [](const QRs& res) -> Block { return std::get<0>(res); });

    auto R = std::make_shared<Tensor<T, 2>>("R matrix of " + in->label(), in->rowspi(), in->colspi());
    std::transform(tmp.cbegin(), tmp.cend(), R->begin(), [](const QRs& res) -> Block { return std::get<1>(res); });

    return std::make_tuple(Q, R);
}

template <typename T>
using SVDResult = std::tuple<SharedTensor<T, 2>, SharedTensor<T, 1>, SharedTensor<T, 2>>;

template <typename T>
auto svd(const SharedTensor<T, 2>& in, bool full_matrices = true, bool compute_uv = true) -> SVDResult<T> {
    using Block = typename Tensor<T, 2>::block_type;
    using SVDs = std::tuple<Block, typename Tensor<T, 1>::block_type, Block>;
    auto tmp = std::vector<SVDs>(in->nirrep());

    std::transform(in->cbegin(), in->cend(), tmp.begin(), [full_matrices, compute_uv](const Block& blk) -> SVDs {
        return xt::linalg::svd(blk, full_matrices, compute_uv);
    });

    auto Kspi = Dimension(in->nirrep());
    for (auto h = 0; h < in->nirrep(); ++h) {
        Kspi[h] = std::min(in->rowspi()[h], in->colspi()[h]);
    }

    auto U =
        std::make_shared<Tensor<T, 2>>("U matrix of " + in->label(), in->rowspi(), full_matrices ? in->rowspi() : Kspi);
    std::transform(tmp.cbegin(), tmp.cend(), U->begin(), [](const SVDs& res) -> Block { return std::get<0>(res); });

    auto S = std::make_shared<Tensor<T, 1>>("S matrix of " + in->label(), Kspi);
    std::transform(tmp.cbegin(), tmp.cend(), S->begin(),
                   [](const SVDs& res) -> typename Tensor<T, 1>::block_type { return std::get<1>(res); });

    auto Vh = std::make_shared<Tensor<T, 2>>("V^H matrix of " + in->label(), full_matrices ? in->colspi() : Kspi,
                                             in->colspi());
    std::transform(tmp.cbegin(), tmp.cend(), Vh->begin(), [](const SVDs& res) -> Block { return std::get<2>(res); });

    return std::make_tuple(U, S, Vh);
}
/*! @}*/

/*! @{ Matrix eigenvalues */
namespace eig_detail {
template <typename T>
using Eigvals = Tensor<T, 1>;

template <typename T>
using Eigvecs = Tensor<T, 2>;

template <typename T>
using Result = std::tuple<typename Eigvals<T>::block_type, typename Eigvecs<T>::block_type>;
}  // namespace eig_detail

template <typename T>
using EigenResult = std::tuple<SharedTensor<T, 1>, SharedTensor<T, 2>>;

template <typename T, typename U = std::complex<typename detail::is_complex<T>::real_type>>
auto eig(const SharedTensor<T, 2>& in) -> EigenResult<U> {
    if (in->symmetry()) {
        throw PSIEXCEPTION("Matrix is non-totally symmetric.");
    }

    auto tmp = std::vector<eig_detail::Result<U>>(in->nirrep());

    std::transform(in->cbegin(), in->cend(), tmp.begin(),
                   [](const typename eig_detail::Eigvecs<T>::block_type& blk) -> eig_detail::Result<U> {
                       return xt::linalg::eig(blk);
                   });

    auto eigvals = std::make_shared<eig_detail::Eigvals<U>>("Eigenvalues of " + in->label(), in->rowspi());
    std::transform(tmp.cbegin(), tmp.cend(), eigvals->begin(), [
    ](const eig_detail::Result<U>& res) -> typename eig_detail::Eigvals<U>::block_type { return std::get<0>(res); });

    auto eigvecs =
        std::make_shared<eig_detail::Eigvecs<U>>("Eigenvectors of " + in->label(), in->rowspi(), in->colspi());
    std::transform(tmp.cbegin(), tmp.cend(), eigvecs->begin(), [
    ](const eig_detail::Result<U>& res) -> typename eig_detail::Eigvecs<U>::block_type { return std::get<1>(res); });

    return std::make_tuple(eigvals, eigvecs);  // NOTE for posterity in C++17 simplifies to return {};
}

template <typename T, typename U = std::complex<typename detail::is_complex<T>::real_type>>
auto eigvals(const SharedTensor<T, 2>& in) -> std::shared_ptr<eig_detail::Eigvals<U>> {
    if (in->symmetry()) {
        throw PSIEXCEPTION("Matrix is non-totally symmetric.");
    }

    auto eigvals = std::make_shared<eig_detail::Eigvals<U>>("Eigenvalues of " + in->label(), in->rowspi());
    std::transform(in->cbegin(), in->cend(), eigvals->begin(),
                   [](const typename Tensor<T, 2>::block_type& blk) ->
                   typename eig_detail::Eigvals<U>::block_type { return xt::linalg::eigvals(blk); });
    return eigvals;
}

template <typename T>
auto eigh(const Tensor<T, 2>& in, char UPLO = 'L') -> std::tuple<Tensor<T, 1>, Tensor<T, 2>> {
    if (in.symmetry()) {
        throw PSIEXCEPTION("Matrix is non-totally symmetric.");
    }

    auto tmp = std::vector<eig_detail::Result<T>>(in.nirrep());

    std::transform(in.cbegin(), in.cend(), tmp.begin(),
                   [UPLO](const typename eig_detail::Eigvecs<T>::block_type& blk) -> eig_detail::Result<T> {
                       return xt::linalg::eigh(blk, UPLO);
                   });

    auto eigvals = eig_detail::Eigvals<T>("Eigenvalues of " + in.label(), in.rowspi());
    std::transform(tmp.cbegin(), tmp.cend(), eigvals.begin(), [
    ](const eig_detail::Result<T>& res) -> typename eig_detail::Eigvals<T>::block_type { return std::get<0>(res); });

    auto eigvecs = eig_detail::Eigvecs<T>("Eigenvectors of " + in.label(), in.rowspi(), in.colspi());
    std::transform(tmp.cbegin(), tmp.cend(), eigvecs.begin(), [
    ](const eig_detail::Result<T>& res) -> typename eig_detail::Eigvecs<T>::block_type { return std::get<1>(res); });

    return std::make_tuple(eigvals, eigvecs);  // NOTE for posterity in C++17 simplifies to return {};
}

template <typename T>
auto eigh(const SharedTensor<T, 2>& in, char UPLO = 'L') -> EigenResult<T> {
    if (in->symmetry()) {
        throw PSIEXCEPTION("Matrix is non-totally symmetric.");
    }

    auto tmp = std::vector<eig_detail::Result<T>>(in->nirrep());

    std::transform(in->cbegin(), in->cend(), tmp.begin(),
                   [UPLO](const typename eig_detail::Eigvecs<T>::block_type& blk) -> eig_detail::Result<T> {
                       return xt::linalg::eigh(blk, UPLO);
                   });

    auto eigvals = std::make_shared<eig_detail::Eigvals<T>>("Eigenvalues of " + in->label(), in->rowspi());
    std::transform(tmp.cbegin(), tmp.cend(), eigvals->begin(), [
    ](const eig_detail::Result<T>& res) -> typename eig_detail::Eigvals<T>::block_type { return std::get<0>(res); });

    auto eigvecs =
        std::make_shared<eig_detail::Eigvecs<T>>("Eigenvectors of " + in->label(), in->rowspi(), in->colspi());
    std::transform(tmp.cbegin(), tmp.cend(), eigvecs->begin(), [
    ](const eig_detail::Result<T>& res) -> typename eig_detail::Eigvecs<T>::block_type { return std::get<1>(res); });

    return std::make_tuple(eigvals, eigvecs);  // NOTE for posterity in C++17 simplifies to return {};
}

template <typename T>
auto eigvalsh(const SharedTensor<T, 2>& in, char UPLO = 'L') -> std::shared_ptr<eig_detail::Eigvals<T>> {
    if (in->symmetry()) {
        throw PSIEXCEPTION("Matrix is non-totally symmetric.");
    }

    auto eigvals = std::make_shared<eig_detail::Eigvals<T>>("Eigenvalues of " + in->label(), in->rowspi());
    std::transform(in->cbegin(), in->cend(), eigvals->begin(),
                   [UPLO](const typename eig_detail::Eigvecs<T>::block_type& blk) ->
                   typename eig_detail::Eigvals<T>::block_type { return xt::linalg::eigvalsh(blk, UPLO); });
    return eigvals;
}
/*! @}*/
}  // namespace psi
