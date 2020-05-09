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

#include <array>
#include <complex>
#include <cstddef>
#include <initializer_list>
#include <memory>
#include <string>
#include <tuple>
#include <type_traits>
#include <vector>

#include <xtensor/xtensor.hpp>
#include <xtensor/xio.hpp>
#include <xtensor-io/xhighfive.hpp>

#include "psi4/psi4-dec.h"
#include "psi4/libpsi4util/exception.h"
#include "psi4/libpsi4util/PsiOutStream.h"

#include "dimension.h"
#include "tensor_impl.h"
#include "tensor_impl_matrix.h"
#include "tensor_impl_vector.h"

namespace psi {
template <typename T, size_t Rank, typename Enable = void>
class Tensor {};

namespace detail {
template <typename T, size_t Rank, std::ptrdiff_t... indices>
struct Accessor<T, Rank, std::integer_sequence<std::ptrdiff_t, indices...>> {
    inline auto get(size_t h, HomTs<indices, std::ptrdiff_t>... vs) const -> T {
        return static_cast<const Tensor<T, Rank>*>(this)->store_[h][{vs...}];
    }
    inline auto get(HomTs<indices, std::ptrdiff_t>... vs) const -> T { return this->get(0, vs...); }
    inline auto get(size_t h, std::initializer_list<std::ptrdiff_t> idxs) const -> T {
        return static_cast<const Tensor<T, Rank>*>(this)->store_[h][idxs];
    }
    inline auto get(std::initializer_list<std::ptrdiff_t> idxs) const -> T { return this->get(0, idxs); }

    inline auto set(size_t h, HomTs<indices, std::ptrdiff_t>... vs, T val) -> void {
        static_cast<Tensor<T, Rank>*>(this)->store_[h][{vs...}] = val;
    }
    inline auto set(HomTs<indices, std::ptrdiff_t>... vs, T val) -> void { this->set(0, vs..., val); }
    inline auto set(size_t h, std::initializer_list<std::ptrdiff_t> idxs, T val) -> void {
        static_cast<Tensor<T, Rank>*>(this)->store_[h][idxs] = val;
    }
    inline auto set(std::initializer_list<std::ptrdiff_t> idxs, T val) -> void { this->set(0, idxs, val); }
};
}  // namespace detail

template <typename T, size_t Rank>
class Tensor<T, Rank, detail::Valid<T, Rank>>
    : public detail::Accessor<T, Rank>,
      public detail::RankDependentImpl<Tensor<T, Rank, detail::Valid<T, Rank>>>,
      public std::enable_shared_from_this<Tensor<T, Rank, detail::Valid<T, Rank>>> {
    friend struct detail::Accessor<T, Rank>;
    friend struct detail::RankDependentImpl<Tensor<T, Rank, detail::Valid<T, Rank>>>;

   public:
    /*! Access rank of Tensor as Tensor<T, Rank>::rank */
    static constexpr size_t rank = Rank;

    /*! Access arithmetic type of Tensor as Tensor<T, Rank>::value_type */
    using value_type = T;

    using shape_type = std::array<size_t, Rank>;
    using index_type = std::array<std::ptrdiff_t, Rank>;
    /*! Basic linear algebra storage object for a block
     * \tparam T the underlying numerical type
     * \tparam Rank rank of the object, i.e. 1 for a vector, 2 for a matrix, etc.
     */
    using block_type = xt::xtensor<T, Rank>;
    using store_type = std::vector<block_type>;

    using accessor = detail::Accessor<T, Rank>;
    using crtp_base = detail::RankDependentImpl<Tensor<T, Rank, detail::Valid<T, Rank>>>;

    /*! @{ Main constructors */
    /*! Labeled, blocked, symmetry-assigned, rank-n CTOR with non-zero fill value
     *  \param[in] label name of the tensor
     *  \param[in] nirrep number of irreps (a.k.a. blocks) in the tensor
     *  \param[in] axes_dimpi dimension of each axis
     *  \param[in] symmetry overall symmetry of the tensor
     *  \param[in] fill_value value of all elements in each block
     *
     *  \note This CTOR initializes the storage with the given fill value (zero by default)
     */
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
            if (nirrep != ax0_n) {
                std::ostringstream err;
                err << "In Tensor CTOR axes dimensions (" << std::to_string(ax0_n) << ") and numbe of irreps ("
                    << std::to_string(nirrep) << ") do not match";
                throw PSIEXCEPTION(err.str());
            }
        }
        for (auto h = 0; h < store_.size(); ++h) {
            shapes_[h][0] = axes_dimpi[0][h];
            for (int ax = 1; ax < Rank; ++ax) shapes_[h][ax] = axes_dimpi[ax][h ^ symmetry];
            store_[h] = xt::full<T>(shapes_[h], fill_value);
        }
    }

    /*! Labeled, blocked, symmetry-assigned, rank-n CTOR
     *  \param[in] label name of the tensor
     *  \param[in] nirrep number of irreps (a.k.a. blocks) in the tensor
     *  \param[in] axes_dimpi dimension of each axis
     *  \param[in] symmetry overall symmetry of the tensor
     *
     *  \note This CTOR _does not_ initialize the storage
     */
    explicit Tensor(const std::string& label, size_t nirrep, const std::array<Dimension, Rank>& axes_dimpi,
                    unsigned int symmetry)
        : symmetry_(symmetry), label_(label), axes_dimpi_(axes_dimpi), shapes_(nirrep), store_(nirrep) {
        if (Rank > 1) {
            auto ax0_n = axes_dimpi[0].n();
            auto all_axes_have_same_size = std::all_of(std::begin(axes_dimpi), std::end(axes_dimpi),
                                                       [ax0_n](const Dimension& dimpi) { return dimpi.n() == ax0_n; });
            if (!all_axes_have_same_size) {
                throw PSIEXCEPTION("In Tensor CTOR axes_dimpi do NOT have same size");
            }
            if (nirrep != ax0_n) {
                std::ostringstream err;
                err << "In Tensor CTOR axes dimensions (" << std::to_string(ax0_n) << ") and numbe of irreps ("
                    << std::to_string(nirrep) << ") do not match";
                throw PSIEXCEPTION(err.str());
            }
        }
        for (auto h = 0; h < store_.size(); ++h) {
            shapes_[h][0] = axes_dimpi[0][h];
            for (int ax = 1; ax < Rank; ++ax) shapes_[h][ax] = axes_dimpi[ax][h ^ symmetry];
        }
    }

    /*! Empty tensor */
    explicit Tensor()
        : Tensor("empty", 1, std::array<Dimension, Rank>{{Dimension(std::vector<Dimension::value_type>{0})}}, 0,
                 static_cast<T>(0)) {}
    /*! @} */

    /*! @{ Rank-n CTORs */
    /*! Labeled, 1-irrep, rank-n CTOR
     *  \param[in] label
     *  \param[in] axes_dimpi
     *  \param[in] fill_value
     */
    explicit Tensor(const std::string& label, const std::array<Dimension, Rank>& axes_dimpi,
                    T fill_value = static_cast<T>(0))
        : Tensor(label, 1, axes_dimpi, 0, fill_value) {}
    /*! Labeled, 1-irrep, rank-n CTOR
     *  \param[in] label
     *  \param[in] axes_dims
     *  \param[in] fill_value
     */
    explicit Tensor(const std::string& label, const std::array<Dimension::value_type, Rank>& axes_dims,
                    T fill_value = static_cast<T>(0))
        : Tensor(label, 1, detail::axes_dimpi(axes_dims), 0, fill_value) {}
    /*! Unlabeled, blocked, symmetry-assigned rank-n CTOR
     *  \param[in] nirrep
     *  \param[in] axes_dimpi
     *  \param[in] symmetry
     *  \param[in] fill_value
     */
    explicit Tensor(size_t nirrep, const std::array<Dimension, Rank>& axes_dimpi, unsigned int symmetry,
                    T fill_value = static_cast<T>(0))
        : Tensor("", nirrep, axes_dimpi, symmetry, fill_value) {}
    /*! Unlabeled, blocked, rank-n CTOR
     *  \param[in] nirrep
     *  \param[in] axes_dimpi
     *  \param[in] fill_value
     */
    explicit Tensor(size_t nirrep, const std::array<Dimension, Rank>& axes_dimpi, T fill_value = static_cast<T>(0))
        : Tensor("", nirrep, axes_dimpi, 0, fill_value) {}
    /*! Unlabeled, 1-irrep, rank-n CTOR
     *  \param[in] axes_dimpi
     *  \param[in] fill_value
     */
    explicit Tensor(const std::array<Dimension, Rank>& axes_dimpi, T fill_value = static_cast<T>(0))
        : Tensor("", 1, axes_dimpi, 0, fill_value) {}
    /*! Unlabeled, 1-irrep, rank-n CTOR
     *  \param[in] axes_dims
     *  \param[in] fill_value
     */
    explicit Tensor(const std::array<Dimension::value_type, Rank>& axes_dims, T fill_value = static_cast<T>(0))
        : Tensor("", 1, detail::axes_dimpi(axes_dims), 0, fill_value) {}
    /*! @}*/

    /*! @{ Rank-1 CTORs */
    /*! Labeled, blocked, rank-1 CTOR
     *  \param[in] label
     *  \param[in] dimpi
     */
    template <size_t Rank_ = Rank, typename = std::enable_if_t<detail::is_rank1_v<Rank_>>>
    explicit Tensor(const std::string& label, const Dimension& dimpi, T fill_value = static_cast<T>(0))
        : Tensor(label, dimpi.n(), std::array<Dimension, Rank>{{dimpi}}, 0, fill_value) {}
    /*! Labeled, 1-irrep, rank-1 CTOR
     *  \param[in] label
     *  \param[in] dim
     */
    template <size_t Rank_ = Rank, typename = std::enable_if_t<detail::is_rank1_v<Rank_>>>
    explicit Tensor(const std::string& label, int dim, T fill_value = static_cast<T>(0))
        : Tensor(label, 1, std::array<Dimension, Rank>{{Dimension(std::vector<Dimension::value_type>{dim})}}, 0,
                 fill_value) {}
    /*! Unlabeled, blocked, rank-1 CTOR
     *  \param[in] dimpi
     */
    template <size_t Rank_ = Rank, typename = std::enable_if_t<detail::is_rank1_v<Rank_>>>
    explicit Tensor(const Dimension& dimpi, T fill_value = static_cast<T>(0))
        : Tensor("", dimpi.n(), std::array<Dimension, Rank>{{dimpi}}, 0, fill_value) {}
    /*! Unlabeled, 1-irrep, rank-1 CTOR
     *  \param[in] dim
     */
    template <size_t Rank_ = Rank, typename = std::enable_if_t<detail::is_rank1_v<Rank_>>>
    explicit Tensor(int dim, T fill_value = static_cast<T>(0))
        : Tensor("", 1, std::array<Dimension, Rank>{{Dimension(std::vector<Dimension::value_type>{dim})}}, 0,
                 fill_value) {}
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
        : Tensor(label, rowspi.n(), std::array<Dimension, Rank>{{rowspi, colspi}}, symmetry, fill_value) {}
    /*! Labeled, blocked, rank-2 CTOR
     *  \param[in] label
     *  \param[in] rowspi
     *  \param[in] colspi
     */
    template <size_t Rank_ = Rank, typename = std::enable_if_t<detail::is_rank2_v<Rank_>>>
    explicit Tensor(const std::string& label, const Dimension& rowspi, const Dimension& colspi,
                    T fill_value = static_cast<T>(0))
        : Tensor(label, rowspi.n(), std::array<Dimension, Rank>{{rowspi, colspi}}, 0, fill_value) {}
    /*! Labeled 1-irrep rank-2 CTOR
     *  \param[in] label
     *  \param[in] rows
     *  \param[in] cols
     */
    template <size_t Rank_ = Rank, typename = std::enable_if_t<detail::is_rank2_v<Rank_>>>
    explicit Tensor(const std::string& label, int rows, int cols, T fill_value = static_cast<T>(0))
        : Tensor(label, 1,
                 std::array<Dimension, Rank>{{Dimension(std::vector<Dimension::value_type>{rows}),
                                              Dimension(std::vector<Dimension::value_type>{cols})}},
                 0, fill_value) {}
    /*! Unlabeled, blocked, symmetry-assigned rank-2 CTOR
     *  \param[in] rowspi
     *  \param[in] colspi
     *  \param[in] symmetry
     */
    template <size_t Rank_ = Rank, typename = std::enable_if_t<detail::is_rank2_v<Rank_>>>
    explicit Tensor(const Dimension& rowspi, const Dimension& colspi, unsigned int symmetry,
                    T fill_value = static_cast<T>(0))
        : Tensor("", rowspi.n(), std::array<Dimension, Rank>{{rowspi, colspi}}, symmetry, fill_value) {}
    /*! Unlabeled blocked rank-2 CTOR
     *  \param[in] rowspi
     *  \param[in] colspi
     */
    template <size_t Rank_ = Rank, typename = std::enable_if_t<detail::is_rank2_v<Rank_>>>
    explicit Tensor(const Dimension& rowspi, const Dimension& colspi, T fill_value = static_cast<T>(0))
        : Tensor("", rowspi.n(), std::array<Dimension, Rank>{{rowspi, colspi}}, 0, fill_value) {}
    /*! Unlabeled 1-irrep rank-2 CTOR
     *  \param[in] rows
     *  \param[in] cols
     */
    template <size_t Rank_ = Rank, typename = std::enable_if_t<detail::is_rank2_v<Rank_>>>
    explicit Tensor(int rows, int cols, T fill_value = static_cast<T>(0))
        : Tensor("", 1,
                 std::array<Dimension, Rank>{{Dimension(std::vector<Dimension::value_type>{rows}),
                                              Dimension(std::vector<Dimension::value_type>{cols})}},
                 0, fill_value) {}
    /*! @}*/

    /*! @{ Iterators */
    using iterator = typename store_type::iterator;
    using const_iterator = typename store_type::const_iterator;

    iterator begin() noexcept { return store_.begin(); }
    const_iterator begin() const noexcept { return store_.begin(); }
    const_iterator cbegin() const noexcept { return store_.cbegin(); }

    iterator end() noexcept { return store_.end(); }
    const_iterator end() const noexcept { return store_.end(); }
    const_iterator cend() const noexcept { return store_.cend(); }

    using reverse_iterator = typename store_type::reverse_iterator;
    using const_reverse_iterator = typename store_type::const_reverse_iterator;

    reverse_iterator rbegin() noexcept { return store_.rbegin(); }
    const_reverse_iterator rbegin() const noexcept { return store_.rbegin(); }
    const_reverse_iterator crbegin() const noexcept { return store_.crbegin(); }

    reverse_iterator rend() noexcept { return store_.rend(); }
    const_reverse_iterator rend() const noexcept { return store_.rend(); }
    const_reverse_iterator crend() const noexcept { return store_.crend(); }
    /*! @}*/

    /*! Return dimension of tensor */
    size_t dim() const noexcept {
        return std::accumulate(std::begin(store_), std::end(store_), 0,
                               [](size_t s, auto& b) { return (s + b.size()); });
    }
    /*! Return dimension of given irrep of tensor
     *  \param[in] h
     */
    size_t dim(size_t h) const noexcept { return store_[h].size(); }

    /*! Return number of irreducible representations */
    size_t nirrep() const noexcept { return store_.size(); }

    /*! Return shapes of blocks */
    const std::vector<shape_type>& shapes() const noexcept { return shapes_; }

    /*! Return block at given irrep
     *  \param[in] h
     */
    const block_type& block(size_t h = 0) const { return store_[h]; }
    /*! Return block at given irrep
     *  \param[in] h
     */
    block_type& block(size_t h = 0) { return store_[h]; }
    block_type operator[](size_t h) const { return store_[h]; }
    /*! Set block at given irrep
     *  \param[in] h
     *  \param[in] block
     */
    void set_block(size_t h, const block_type& block) { store_[h] = block; }
    /*! Set block at totally symmetric irrep
     *  \param[in] block
     */
    void set_block(const block_type& block) { this->set_block(0, block); }
    block_type& operator[](size_t h) { return store_[h]; }

    // using crtp_base::get;  // NOTE This is to make the rank-dependent get-s are accessible
    using accessor::get;
    T get(size_t h, index_type idxs) const { return store_[h][idxs]; }
    T get(index_type idxs) const { return store_[0][idxs]; }

    // using crtp_base::set;  // NOTE This is to make the rank-dependent set-s are accessible
    using accessor::set;
    void set(size_t h, index_type idxs, T val) { store_[h][idxs] = val; }
    void set(index_type idxs, T val) { store_[0][idxs] = val; }

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

    /*! Returns pointer to given irrep
     *  \param[in] h
     */
    T* data(size_t h = 0) { return store_.at(h).data(); }
    const T* data(size_t h = 0) const { return store_.at(h).data(); }

    PSI_DEPRECATED(
        "Using `Tensor::pointer` instead of `Tensor::data` is deprecated, and in 1.4 it will "
        "stop working")
    T* pointer(size_t h = 0) { return store_.at(h).data(); }
    PSI_DEPRECATED(
        "Using `Tensor::pointer` instead of `Tensor::data` is deprecated, and in 1.4 it will "
        "stop working")
    const T* pointer(size_t h = 0) const { return store_.at(h).data(); }

    std::string repr() const noexcept { return crtp_base::cxxClassName(); }

    std::string str() const noexcept { return format(crtp_base::cxxClassName()); }

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
        for (auto h = 0; h < store_.size(); ++h) {
            retval << "  Irrep: " << h + 1 << " Shape: " << detail::print_shape(store_[h].shape()) << std::endl;
            if (store_[h].size() == 0) {
                retval << "    (empty)" << std::endl;
            } else {
                retval << xt::print_options::threshold(10000) << xt::print_options::line_width(120)
                       << xt::print_options::precision(14) << store_[h] << std::endl;
            }
            retval << std::endl;
        }
        return retval.str();
    }

    PSI_DEPRECATED(
        "Using `Tensor::print` instead of `print` is deprecated, and in 1.4 it will "
        "stop working")
    void print(std::string out = "outfile", const std::string extra = "") const noexcept {
        std::shared_ptr<psi::PsiOutStream> printer = (out == "outfile" ? outfile : std::make_shared<PsiOutStream>(out));
        printer->Printf(format(extra));
        printer->Printf("\n");
    }

    void dump_hdf5(const std::string& fname, const std::string& path) const {
        for (auto h = 0; h < store_.size(); ++h) {
            xt::dump_hdf5(fname, path + "/irrep_" + std::to_string(h), store_[h]);
        }
    }

   protected:
    unsigned int symmetry_{0};
    std::string label_;
    std::array<Dimension, Rank> axes_dimpi_{};
    std::vector<shape_type> shapes_{};
    store_type store_{};
};

template <typename T, size_t Rank>
using SharedTensor = std::shared_ptr<Tensor<T, Rank>>;

template <typename T, size_t Rank>
void print(const Tensor<T, Rank>& t, std::string out = "outfile", const std::string extra = "") noexcept {
    std::shared_ptr<psi::PsiOutStream> printer = (out == "outfile" ? outfile : std::make_shared<PsiOutStream>(out));
    printer->Printf(t.format(extra));
    printer->Printf("\n");
}

template <typename T, size_t Rank>
void print(const SharedTensor<T, Rank>& t, std::string out = "outfile", const std::string extra = "") noexcept {
    std::shared_ptr<psi::PsiOutStream> printer = (out == "outfile" ? outfile : std::make_shared<PsiOutStream>(out));
    printer->Printf(t->format(extra));
    printer->Printf("\n");
}
}  // namespace psi
