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
#include <cmath>
#include <memory>
#include <sstream>
#include <vector>

#include <xtensor-blas/xlinalg.hpp>
#include <xtensor/xadapt.hpp>
#include <xtensor/xfixed.hpp>

namespace psi {
template <typename T>
using Vector3 = xt::xtensor_fixed<T, xt::xshape<3>>;

template <typename T>
using SharedVector3 = std::shared_ptr<Vector3<T>>;

template <typename T>
Vector3<T> from_array(const std::array<T, 3>& a) noexcept {
    return xt::adapt(a, {3});
}

template <typename T>
Vector3<T> from_vector(const std::vector<T>& v) noexcept {
    return xt::adapt(v, {3});
}

template <typename T>
Vector3<T> from_Ts(T x, T y, T z) noexcept {
    return {x, y, z};
}

template <typename T>
Vector3<T> from_Tstack(T stack[3]) noexcept {
    return xt::adapt(stack, {3});
}

template <typename T>
Vector3<T> from_Tptr(T* ptr) noexcept {
    std::vector<std::size_t> shape = {3};
    return xt::adapt(ptr, 3, xt::no_ownership(), shape);
}

template <typename T>
Vector3<T> from_Tptr(const T* ptr) noexcept {
    std::vector<std::size_t> shape = {3};
    return xt::adapt(ptr, 3, xt::no_ownership(), shape);
}

template <typename T>
inline auto dot(const Vector3<T>& u, const Vector3<T>& v) noexcept {
    return xt::linalg::vdot(u, v);
}

template <typename E1, typename E2>
inline auto dot(E1&& u, E2&& v) noexcept {
    return xt::linalg::vdot(u, v);
}

template <typename T>
inline auto norm(const Vector3<T>& u) noexcept {
    return xt::linalg::norm(u);
}

template <typename E>
inline auto norm(E&& u) noexcept {
    return xt::linalg::norm(u);
}

template <typename T>
inline auto cross(const Vector3<T>& u, const Vector3<T>& v) noexcept {
    return xt::linalg::cross(u, v);
}

template <typename E1, typename E2>
inline auto cross(E1&& u, E2&& v) noexcept {
    return xt::eval(xt::linalg::cross(u, v));
}

template <class T>
inline auto normalize(const Vector3<T>& u) noexcept {
    return (u / xt::linalg::norm(u));
}

template <class E>
inline auto normalize(E&& u) noexcept {
    return xt::eval(u / xt::linalg::norm(u));
}

template <typename T>
inline auto distance(const Vector3<T>& u, const Vector3<T>& v) noexcept {
    return xt::linalg::norm(u - v);
}

template <typename T>
inline auto angle(const Vector3<T>& u, const Vector3<T>& v) noexcept {
    Vector3<T> _u = normalize(u);
    Vector3<T> _v = normalize(v);
    return std::acos(psi::dot(_u, _v));
}

template <typename E1, typename E2>
inline auto angle(E1&& u, E2&& v) noexcept {
    auto _u = normalize(u);
    auto _v = normalize(v);
    return std::acos(xt::linalg::vdot(_u, _v));
}

template <typename T>
auto str(const Vector3<T>& u) noexcept -> std::string {
    std::ostringstream os;
    os << xt::print_options::threshold(10000) << xt::print_options::line_width(120) << xt::print_options::precision(14)
       << u << std::endl;
    return os.str();
}

template <typename T>
inline Vector3<T> perp_unit(const Vector3<T>& u, const Vector3<T>& v) noexcept {
    // try cross product
    auto result = cross(u, v);
    auto resultdotresult = dot(result, result);

    if (resultdotresult < 1.e-16) {
        // cross product is too small to normalize
        // find the largest of this and v
        double dotprodu = dot(u, u);
        double dotprodv = dot(v, v);
        auto d = from_Ts(0.0, 0.0, 0.0);
        double dotprodd;
        if (dotprodu < dotprodv) {
            d = v;
            dotprodd = dotprodv;
        } else {
            d = u;
            dotprodd = dotprodu;
        }

        // see if d is big enough
        if (dotprodd < 1.e-16) {
            // choose an arbitrary vector, since the biggest vector is small
            return from_Ts(1.0, 0.0, 0.0);
        } else {
            // choose a vector prependicular to d
            // choose it in one of the planes xy, xz, yz
            // choose the plane to be that which contains the two largest
            // components of d
            double absd[3];
            absd[0] = std::abs(d[0]);
            absd[1] = std::abs(d[1]);
            absd[2] = std::abs(d[2]);
            int axis0, axis1;
            if ((absd[1] - absd[0]) > 1.0e-12) {
                axis0 = 1;
                if ((absd[2] - absd[0]) > 1.0e-12) {
                    axis1 = 2;
                } else {
                    axis1 = 0;
                }
            } else {
                axis0 = 0;
                if ((absd[2] - absd[1]) > 1.0e-12) {
                    axis1 = 2;
                } else {
                    axis1 = 1;
                }
            }

            result[0] = 0.0;
            result[1] = 0.0;
            result[2] = 0.0;
            // do the pi/2 rotation in the plane
            result[axis0] = d[axis1];
            result[axis1] = -d[axis0];
        }
        return normalize(result);
    } else {
        // normalize the cross product and return the result
        return normalize(result);
    }
}

template <typename T>
inline Vector3<double> rotate(const Vector3<T>& in, double theta, const Vector3<T>& axis) noexcept {
    Vector3<T> unitaxis = normalize(axis);

    // split into parallel and perpendicular components along axis
    Vector3<T> parallel = unitaxis * (dot(in, unitaxis) / dot(unitaxis, unitaxis));
    Vector3<T> perpendicular = in - parallel;

    // form unit vector perpendicular to parallel and perpendicular
    Vector3<T> third_axis = perp_unit(axis, perpendicular) * xt::linalg::norm(perpendicular);

    return parallel + cos(theta) * perpendicular + sin(theta) * third_axis;
}
}  // namespace psi
