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
#include <sstream>
#include <string>
#include <tuple>
#include <type_traits>
#include <vector>

#include <xtensor/xtensor.hpp>

namespace xt {
template <class T, class S>
inline auto full(S shape, T fill_value) noexcept {
    return broadcast(static_cast<T>(fill_value), std::forward<S>(shape));
}
}  // namespace xt

namespace psi {
namespace detail {
template <typename T>
struct is_double : std::integral_constant<bool, std::is_same<double, typename std::remove_cv<T>::type>::value> {};

// NOTE for posterity this can be declared inline in C++17
template <typename T>
constexpr bool is_double_v = is_double<T>::value;

template <typename T>
struct is_complex : std::false_type {
    using real_type = T;
};

template <typename T>
struct is_complex<std::complex<T>> : std::true_type {
    using real_type = T;
};

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

using rank0 = std::integral_constant<size_t, 0>;

template <size_t Rank>
struct is_rank0 : std::integral_constant<bool, std::is_same<std::integral_constant<size_t, Rank>, rank0>::value> {};

// NOTE for posterity this can be declared inline in C++17
template <size_t Rank>
constexpr bool is_rank0_v = is_rank0<Rank>::value;

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
struct is_rankn : std::integral_constant<bool, !is_rank0_v<Rank> || !is_rank1_v<Rank> || !is_rank2_v<Rank>> {};

// NOTE for posterity this can be declared inline in C++17
template <size_t Rank>
constexpr bool is_rankn_v = is_rankn<Rank>::value;

template <typename T>
struct Type2String final {
    static std::string full() { return ""; }
    static std::string suffix() { return ""; }
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

/// Implements rank-dependent functions
template <typename T>
struct RankDependentImpl {
    // Here to satisfy the compiler
    T get() const {}
    void set() {}

    /*! C++ string representation of type */
    std::string cxxClassName() const noexcept {
        return "Tensor" + std::to_string(T::rank) + "<" + Type2String<typename T::value_type>::full() + ">";
    }
    /*! Python string representation of type */
    static std::string pyClassName() noexcept {
        return "Tensor" + std::to_string(T::rank) + "_" + Type2String<typename T::value_type>::suffix();
    }
};

/*! STL collection as a string
 *  \param[in] coll
 *  \param[in] l
 *  \param[in] r
 *  \param[in] sep
 */
template <typename T>
auto stream_collection(const T &coll, std::string l = "[", std::string r = "]", std::string sep = ", ") noexcept
    -> std::string {
    std::ostringstream os;
    bool first = true;
    os << l;
    for (auto elem : coll) {
        if (!first) os << sep;
        os << elem;
        first = false;
    }
    os << r;
    return os.str();
}

template <size_t Rank>
std::string print_shape(const std::array<size_t, Rank> &shape) noexcept {
    std::ostringstream retval;
    retval << stream_collection(shape, "{", "}");
    return retval.str();
}

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

template <typename T, size_t N>
auto operator<<(std::ostream &os, const std::array<T, N> &coll) -> std::ostream & {
    return (os << detail::stream_collection(coll));
}

template <typename T>
auto operator<<(std::ostream &os, const std::vector<T> &coll) -> std::ostream & {
    return (os << detail::stream_collection(coll));
}
}  // namespace psi
