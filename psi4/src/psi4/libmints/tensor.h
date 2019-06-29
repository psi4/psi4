/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2018 The Psi4 Developers.
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

#include <array>
#include <complex>
#include <memory>
#include <string>
#include <tuple>
#include <type_traits>
#include <vector>

#include <xtensor/xtensor.hpp>
#include <xtensor/xio.hpp>

#include "psi4/libpsi4util/exception.h"

#include "dimension.h"

namespace xt {
template <class T, class S>
inline auto full(S shape, T fill_value) noexcept {
    return broadcast(static_cast<T>(fill_value), std::forward<S>(shape));
}
}  // namespace xt

namespace psi {
namespace detail {
template <typename T>
struct is_complex : std::false_type {};

template <typename T>
struct is_complex<std::complex<T>> : std::true_type {};

// NOTE for posterity this can be declared inline in C++17
template <typename T>
constexpr bool is_complex_v = is_complex<T>::value;

/*! Type trait for tensorisable types
 *
 * These are used for SFINAE in the Tensor class.
 * Tensorisable types are:
 *   - int
 *   - float
 *   - double
 *   - complex of any of the above
 */
template <typename T>
struct is_tensorisable : std::integral_constant<bool, std::is_arithmetic<T>::value || is_complex_v<T>> {};

// NOTE for posterity this can be declared inline in C++17
template <typename T>
constexpr bool is_tensorisable_v = is_tensorisable<T>::value;

using rank1 = std::integral_constant<size_t, 1>;

template <size_t Rank>
struct is_rank1 : std::integral_constant<bool, std::is_same<std::integral_constant<size_t, Rank>, rank1>::value> {};

// NOTE for posterity this can be declared inline in C++17
template <size_t Rank>
constexpr bool is_rank1_v = is_rank1<Rank>::value;

using rank2 = std::integral_constant<size_t, 2>;

template <size_t Rank>
struct is_rank2 : std::integral_constant<bool, std::is_same<std::integral_constant<size_t, Rank>, rank2>::value> {};

// NOTE for posterity this can be declared inline in C++17
template <size_t Rank>
constexpr bool is_rank2_v = is_rank2<Rank>::value;

template <size_t Rank>
struct is_rankn : std::integral_constant<bool, !is_rank1_v<Rank> || !is_rank2_v<Rank>> {};

// NOTE for posterity this can be declared inline in C++17
template <size_t Rank>
constexpr bool is_rankn_v = is_rankn<Rank>::value;

template <size_t Rank>
std::string class_name() {
    return "Tensor" + std::to_string(Rank);
}

template <>
std::string class_name<1>() {
    return "Vector";
}

template <>
std::string class_name<2>() {
    return "Matrix";
}

template <typename T>
struct Type2String final {
    static std::string full() {}
    static std::string suffix() {}
};

template <>
struct Type2String<int> final {
    static std::string full() { return "int"; }
    static std::string suffix() { return "I"; }
};

template <>
struct Type2String<float> final {
    static std::string full() { return "float"; }
    static std::string suffix() { return "F"; }
};

template <>
struct Type2String<double> final {
    static std::string full() { return "double"; }
    static std::string suffix() { return "D"; }
};

template <>
struct Type2String<std::complex<double>> final {
    static std::string full() { return "complex<double>"; }
    static std::string suffix() { return "CD"; }
};

template <size_t Rank>
std::string print_shape(const std::array<size_t, Rank>& shape) {
    std::ostringstream retval;
    retval << "{";
    std::string sep;
    for (const auto& s : shape) {
        retval << sep << s;
        sep = ", ";
    }
    retval << "}";
    return retval.str();
}
}  // namespace detail

/*! Basic linear algebra storage object
 * \tparam T the underlying numerical type
 * \tparam Rank rank of the object, i.e. 1 for a vector, 2 for a matrix, etc.
 */
template <typename T, size_t Rank>
using Storage = xt::xtensor<T, Rank>;

template <typename T, size_t Rank>
class Tensor {
   public:
    /*! Access rank of Tensor as Tensor<T, Rank>::rank */
    static constexpr size_t rank = Rank;
    /*! Access arithmetic type of Tensor as Tensor<T, Rank>::value_type */
    using value_type = T;
    using Shape = std::array<size_t, Rank>;
    /*! C++ string representation of type */
    static std::string cxxClassName() {
        return detail::class_name<Rank>() + "<" + detail::Type2String<T>::full() + ">";
    }
    /*! Python string representation of type */
    static std::string pyClassName() { return detail::class_name<Rank>() + "_" + detail::Type2String<T>::suffix(); }

    /*! Labeled, blocked, symmetry-assigned, rank-n CTOR
     *  \param[in] label name of the tensor
     *  \param[in] nirrep number of irreps (a.k.a. blocks) in the tensor
     *  \param[in] axes_dimpi dimension of each axis
     *  \param[in] symmetry overall symmetry of the tensor
     *  \param[in] fill_value value of all elements in each block
     *
     *  This is the only CTOR that does actual work. All other CTORs delegate to
     *  this one and returns an tensor with all zero blocks.
     */
    template <typename T_ = T, typename = std::enable_if_t<detail::is_tensorisable_v<T_>>>
    explicit Tensor(const std::string& label, size_t nirrep, const std::array<Dimension, Rank>& axes_dimpi,
                    unsigned int symmetry, T fill_value = static_cast<T>(0))
        : symmetry_(symmetry), label_(label), axes_dimpi_(axes_dimpi), shapes_(nirrep), store_(nirrep) {
        if (Rank > 1) {
            auto ax0_n = axes_dimpi[0].n();
            auto all_axes_have_same_size = std::all_of(std::begin(axes_dimpi), std::end(axes_dimpi),
                                                       [ax0_n](const Dimension& dimpi) { return dimpi.n() == ax0_n; });
            if (!all_axes_have_same_size) {
                throw PSIEXCEPTION("In Tensor CTOR axes_dimpi do NOT have same size");
            }
        }
        for (int h = 0; h < store_.size(); ++h) {
            shapes_[h][0] = axes_dimpi[0][h];
            for (int ax = 1; ax < Rank; ++ax) shapes_[h][ax] = axes_dimpi[ax][h ^ symmetry];
            store_[h] = xt::full<T>(shapes_[h], fill_value);
        }
    }

    /*! @{ Rank-n CTORs */
    /*! Labeled, 1-irrep, rank-n CTOR
     *  \param[in] label
     *  \param[in] axes_dimpi
     *  \param[in] fill_value
     */
    template <typename T_ = T, typename = std::enable_if_t<detail::is_tensorisable_v<T_>>>
    explicit Tensor(const std::string& label, const std::array<Dimension, Rank>& axes_dimpi,
                    T fill_value = static_cast<T>(0))
        : Tensor(label, 1, axes_dimpi, 0, fill_value) {}
    /*! Unlabeled, blocked, symmetry-assigned rank-n CTOR
     *  \param[in] nirrep
     *  \param[in] axes_dimpi
     *  \param[in] symmetry
     *  \param[in] fill_value
     */
    template <typename T_ = T, typename = std::enable_if_t<detail::is_tensorisable_v<T_>>>
    explicit Tensor(size_t nirrep, const std::array<Dimension, Rank>& axes_dimpi, unsigned int symmetry,
                    T fill_value = static_cast<T>(0))
        : Tensor("", nirrep, axes_dimpi, symmetry, fill_value) {}
    /*! Unlabeled, blocked, rank-n CTOR
     *  \param[in] nirrep
     *  \param[in] axes_dimpi
     */
    template <typename T_ = T, typename = std::enable_if_t<detail::is_tensorisable_v<T_>>>
    explicit Tensor(size_t nirrep, const std::array<Dimension, Rank>& axes_dimpi, T fill_value = static_cast<T>(0))
        : Tensor("", nirrep, axes_dimpi, 0, fill_value) {}
    /*! Unlabeled, 1-irrep, rank-n CTOR
     *  \param[in] axes_dimpi
     */
    template <typename T_ = T, typename = std::enable_if_t<detail::is_tensorisable_v<T_>>>
    explicit Tensor(const std::array<Dimension, Rank>& axes_dimpi, T fill_value = static_cast<T>(0))
        : Tensor("", 1, axes_dimpi, 0, fill_value) {}
    /*! @}*/

    /*! @{ Rank-1 CTORs */
    /*! Labeled, blocked, rank-1 CTOR
     *  \param[in] label
     *  \param[in] dimpi
     */
    template <size_t Rank_ = Rank, typename = std::enable_if_t<detail::is_rank1_v<Rank_>>>
    explicit Tensor(const std::string& label, const Dimension& dimpi, T fill_value = static_cast<T>(0))
        : Tensor(label, dimpi.n(), std::array<Dimension, Rank>{dimpi}, 0, fill_value) {}
    /*! Labeled, 1-irrep, rank-1 CTOR
     *  \param[in] label
     *  \param[in] dim
     * FIXME int vs. size_t in this CTOR
     */
    template <size_t Rank_ = Rank, typename = std::enable_if_t<detail::is_rank1_v<Rank_>>>
    explicit Tensor(const std::string& label, int dim, T fill_value = static_cast<T>(0))
        : Tensor(label, 1, std::array<Dimension, Rank>{Dimension(std::vector<int>{dim})}, 0, fill_value) {}
    /*! Unlabeled, blocked, rank-1 CTOR
     *  \param[in] dimpi
     */
    template <size_t Rank_ = Rank, typename = std::enable_if_t<detail::is_rank1_v<Rank_>>>
    explicit Tensor(const Dimension& dimpi, T fill_value = static_cast<T>(0))
        : Tensor("", dimpi.n(), std::array<Dimension, Rank>{dimpi}, 0, fill_value) {}
    /*! Unlabeled, 1-irrep, rank-1 CTOR
     *  \param[in] dim
     * FIXME int vs. size_t in this CTOR
     */
    template <size_t Rank_ = Rank, typename = std::enable_if_t<detail::is_rank1_v<Rank_>>>
    explicit Tensor(int dim, T fill_value = static_cast<T>(0))
        : Tensor("", 1, std::array<Dimension, Rank>{Dimension(std::vector<int>{dim})}, 0, fill_value) {}
    /*! @}*/

    /*! @{ Rank-2 CTORs */
    /*! Labeled, blocked, symmetry-assigned, rank-2 CTOR
     *  \param[in] label
     *  \param[in] rowspi
     *  \param[in] colspi
     *  \param[in] symmetry
     */
    template <size_t Rank_ = Rank, typename = std::enable_if_t<detail::is_rank2_v<Rank_>>>
    explicit Tensor(const std::string& label, const Dimension& rowspi, const Dimension& colspi, unsigned int symmetry,
                    T fill_value = static_cast<T>(0))
        : Tensor(label, rowspi.n(), std::array<Dimension, Rank>{rowspi, colspi}, symmetry, fill_value) {}
    /*! Labeled, blocked, rank-2 CTOR
     *  \param[in] label
     *  \param[in] rowspi
     *  \param[in] colspi
     */
    template <size_t Rank_ = Rank, typename = std::enable_if_t<detail::is_rank2_v<Rank_>>>
    explicit Tensor(const std::string& label, const Dimension& rowspi, const Dimension& colspi,
                    T fill_value = static_cast<T>(0))
        : Tensor(label, rowspi.n(), std::array<Dimension, Rank>{rowspi, colspi}, 0, fill_value) {}
    /*! Labeled 1-irrep rank-2 CTOR
     *  \param[in] label
     *  \param[in] rows
     *  \param[in] cols
     * FIXME int vs. size_t in this CTOR
     */
    template <size_t Rank_ = Rank, typename = std::enable_if_t<detail::is_rank2_v<Rank_>>>
    explicit Tensor(const std::string& label, int rows, int cols, T fill_value = static_cast<T>(0))
        : Tensor(label, 1,
                 std::array<Dimension, Rank>{Dimension(std::vector<int>{rows}), Dimension(std::vector<int>{cols})}, 0,
                 fill_value) {}
    /*! Unlabeled, blocked, symmetry-assigned rank-2 CTOR
     *  \param[in] rowspi
     *  \param[in] colspi
     *  \param[in] symmetry
     */
    template <size_t Rank_ = Rank, typename = std::enable_if_t<detail::is_rank2_v<Rank_>>>
    explicit Tensor(const Dimension& rowspi, const Dimension& colspi, unsigned int symmetry,
                    T fill_value = static_cast<T>(0))
        : Tensor("", rowspi.n(), std::array<Dimension, Rank>{rowspi, colspi}, symmetry, fill_value) {}
    /*! Unlabeled blocked rank-2 CTOR
     *  \param[in] rowspi
     *  \param[in] colspi
     */
    template <size_t Rank_ = Rank, typename = std::enable_if_t<detail::is_rank2_v<Rank_>>>
    explicit Tensor(const Dimension& rowspi, const Dimension& colspi, T fill_value = static_cast<T>(0))
        : Tensor("", rowspi.n(), std::array<Dimension, Rank>{rowspi, colspi}, 0, fill_value) {}
    /*! Unlabeled 1-irrep rank-2 CTOR
     *  \param[in] rows
     *  \param[in] cols
     * FIXME int vs. size_t in this CTOR
     */
    template <size_t Rank_ = Rank, typename = std::enable_if_t<detail::is_rank2_v<Rank_>>>
    explicit Tensor(int rows, int cols, T fill_value = static_cast<T>(0))
        : Tensor("", 1,
                 std::array<Dimension, Rank>{Dimension(std::vector<int>{rows}), Dimension(std::vector<int>{cols})}, 0,
                 fill_value) {}
    /*! @}*/

    /*! Return dimension of tensor */
    size_t dim() const noexcept {
        return std::accumulate(std::begin(store_), std::end(store_), 0,
                               [](size_t s, auto& b) { return (s + b.size()); });
    }

    /*! Return number of irreducible representations */
    size_t nirrep() const noexcept { return store_.size(); }

    /*! Return shapes of blocks */
    const std::vector<Shape>& shapes() const noexcept { return shapes_; }

    /*! Return block at given irrep
     *  \param[in] h
     */
    const Storage<T, Rank>& block(size_t h = 0) const { return store_[h]; }
    /*! Return block at given irrep
     *  \param[in] h
     */
    Storage<T, Rank>& block(size_t h = 0) { return store_[h]; }
    Storage<T, Rank> operator[](size_t h) const { return store_[h]; }
    /*! Set block at given irrep
     *  \param[in] h
     *  \param[in] block
     */
    void set_block(size_t h, const Storage<T, Rank>& block) { store_[h] = block; }
    Storage<T, Rank>& operator[](size_t h) { return store_[h]; }

    std::string label() const noexcept { return label_; }
    void set_label(const std::string& label) noexcept { label_ = label; }

    unsigned int symmetry() const noexcept { return symmetry_; }
    void set_symmetry(unsigned int s) noexcept { symmetry_ = s; }

    /*! Returns Dimension objects of all axes */
    const std::array<Dimension, Rank>& axes_dimpi() const noexcept { return axes_dimpi_; }

    /*! Returns Dimension object for given axis
     * \param[in] ax
     */
    const Dimension& axes_dimpi(size_t ax) const noexcept { return axes_dimpi_.at(ax); }

    /// Returns the dimension array for a rank-1 tensor aka a vector
    template <size_t Rank_ = Rank, typename = std::enable_if_t<detail::is_rank1_v<Rank_>>>
    const Dimension& dimpi() const {
        return axes_dimpi_[0];
    }

    /// Returns the dimension array for the rows of a rank-2 tensor aka a matrix
    template <size_t Rank_ = Rank, typename = std::enable_if_t<detail::is_rank2_v<Rank_>>>
    const Dimension& rowspi() const {
        return axes_dimpi_[0];
    }
    /*! Returns the number of rows in given irrep
     *  \param[in] h
     */
    template <size_t Rank_ = Rank, typename = std::enable_if_t<detail::is_rank2_v<Rank_>>>
    size_t rows(size_t h = 0) const {
        return axes_dimpi_[0][h];
    }

    /// Returns the dimension array for the rows of a rank-2 tensor aka a matrix
    template <size_t Rank_ = Rank, typename = std::enable_if_t<detail::is_rank2_v<Rank_>>>
    const Dimension& colspi() const {
        return axes_dimpi_[1];
    }
    /*! Returns the number of columns in given irrep
     *  \param[in] h
     */
    template <size_t Rank_ = Rank, typename = std::enable_if_t<detail::is_rank2_v<Rank_>>>
    size_t cols(size_t h = 0) const {
        return axes_dimpi_[1][h];
    }

    /*! Returns pointer to given irrep
     *  \param[in] h
     */
    T* data(size_t h = 0) { return store_.at(h).data(); }
    const T* data(size_t h = 0) const { return store_.at(h).data(); }

    std::string repr() const noexcept { return cxxClassName(); }

    std::string str() const noexcept { return format(cxxClassName()); }

    /*! Tensor formatter
     *  \param[in] extra
     * TODO Generalise to any rank
     */
    std::string format(const std::string extra = "") const noexcept {
        std::ostringstream retval;
        // Title and stuff
        if (!label_.empty()) {
            retval << "  ## " << label_ << " " << extra << " (Symmetry " << symmetry_ << ") ##\n" << std::endl;
        }
        // Blocks
        for (size_t h = 0; h < store_.size(); ++h) {
            retval << "  Irrep: " << h + 1 << " Shape: " << detail::print_shape(store_[h].shape()) << std::endl;
            if (store_[h].size() == 0) {
                retval << "    (empty)" << std::endl;
            } else {
                retval << xt::print_options::line_width(120) << xt::print_options::precision(14) << store_[h]
                       << std::endl;
            }
            retval << std::endl;
        }
        return retval.str();
    }

   protected:
    unsigned int symmetry_{0};
    std::string label_;
    std::array<Dimension, Rank> axes_dimpi_{};
    std::vector<Shape> shapes_{};
    std::vector<Storage<T, Rank>> store_{};
};

template <typename T, size_t Rank>
using SharedTensor = std::shared_ptr<Tensor<T, Rank>>;
}  // namespace psi
