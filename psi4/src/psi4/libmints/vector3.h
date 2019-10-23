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
#include <sstream>

#include <xtensor/xfixed.hpp>
#include <xtensor-blas/xlinalg.hpp>

#include "psi4/pragma.h"

namespace psi {
template <typename T>
using Vector3_ = xt::xtensor_fixed<T, xt::xshape<3, 1>>;

template <typename T>
T dot(const Vector3_<T>& u, const Vector3_<T>& v) noexcept {
    return xt::linalg::vdot(u, v);
}

template <typename T>
double norm(const Vector3_<T>& u) noexcept {
    return xt::linalg::norm(u);
}

template <typename T>
Vector3_<T> cross(const Vector3_<T>& u, const Vector3_<T>& v) noexcept {
    return xt::linalg::cross(u, v);
}

template <class T>
Vector3_<T> normalize(const Vector3_<T>& u) noexcept {
    return (u / xt::linalg::norm(u));
}

template <typename T>
double distance(const Vector3_<T>& u, const Vector3_<T>& v) noexcept {
    return xt::linalg::norm(u - v);
}

template <typename T>
double angle(const Vector3_<T>& u, const Vector3_<T>& v) noexcept {
    Vector3_<T> _u = normalize(u);
    Vector3_<T> _v = normalize(v);
    return std::acos(psi::dot(_u, _v));
}

template <typename T>
Vector3_<T> perp_unit(const Vector3_<T>& u, const Vector3_<T>& v) noexcept {
    // try cross product
    auto result = cross(u, v);
    auto resultdotresult = dot(result, result);

    if (resultdotresult < 1.e-16) {
        // cross product is too small to normalize
        // find the largest of this and v
        double dotprodu = dot(u, u);
        double dotprodv = dot(v, v);
        auto d = Vector3_<T>({{0.0, 0.0, 0.0}});
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
            return Vector3_<T>({{1.0, 0.0, 0.0}});
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
Vector3_<T> rotate(const Vector3_<T>& in, double theta, const Vector3_<T>& axis) noexcept {
    Vector3_<T> unitaxis = normalize(axis);

    // split into parallel and perpendicular components along axis
    Vector3_<T> parallel = unitaxis * (dot(in, unitaxis) / dot(unitaxis, unitaxis));
    Vector3_<T> perpendicular = in - parallel;

    // form unit vector perpendicular to parallel and perpendicular
    Vector3_<T> third_axis = perp_unit(axis, perpendicular) * xt::linalg::norm(perpendicular);

    return parallel + cos(theta) * perpendicular + sin(theta) * third_axis;
}

class Vector3;
inline Vector3 operator*(double, const Vector3&);

class PSI_API Vector3 {
   private:
    double v_[3];

   public:
    Vector3() { v_[0] = v_[1] = v_[2] = 0.0; }
    Vector3(const double p[3]) {
        v_[0] = p[0];
        v_[1] = p[1];
        v_[2] = p[2];
    }
    Vector3(double d) { v_[0] = v_[1] = v_[2] = d; }
    Vector3(double x, double y, double z) {
        v_[0] = x;
        v_[1] = y;
        v_[2] = z;
    }
    Vector3(const Vector3& p) {
        v_[0] = p.v_[0];
        v_[1] = p.v_[1];
        v_[2] = p.v_[2];
    }
    Vector3(const std::array<double, 3> p) {
        v_[0] = p[0];
        v_[1] = p[1];
        v_[2] = p[2];
    }

    void operator=(const double* x) {
        v_[0] = x[0];
        v_[1] = x[1];
        v_[2] = x[2];
    }
    void operator=(const Vector3& x) {
        v_[0] = x.v_[0];
        v_[1] = x.v_[1];
        v_[2] = x.v_[2];
    }
    void operator=(double d) {
        v_[0] = d;
        v_[1] = d;
        v_[2] = d;
    }
    void operator+=(const Vector3& x) {
        v_[0] += x.v_[0];
        v_[1] += x.v_[1];
        v_[2] += x.v_[2];
    }
    void operator-=(const Vector3& x) {
        v_[0] -= x.v_[0];
        v_[1] -= x.v_[1];
        v_[2] -= x.v_[2];
    }
    void operator*=(double m) {
        v_[0] *= m;
        v_[1] *= m;
        v_[2] *= m;
    }
    void operator/=(double m) {
        v_[0] /= m;
        v_[1] /= m;
        v_[2] /= m;
    }

    Vector3 operator*(double d) const { return d * (*this); }
    Vector3 operator/(double d) const {
        Vector3 result;
        result[0] = v_[0] / d;
        result[1] = v_[1] / d;
        result[2] = v_[2] / d;
        return result;
    }
    Vector3 operator+(const Vector3& x) {
        Vector3 result;
        result.v_[0] = v_[0] + x.v_[0];
        result.v_[1] = v_[1] + x.v_[1];
        result.v_[2] = v_[2] + x.v_[2];
        return result;
    }
    Vector3 operator-(const Vector3& x) const {
        Vector3 result;
        result.v_[0] = v_[0] - x.v_[0];
        result.v_[1] = v_[1] - x.v_[1];
        result.v_[2] = v_[2] - x.v_[2];
        return result;
    }
    Vector3 operator*(const Vector3& x) const {
        Vector3 result;
        result.v_[0] = v_[0] * x.v_[0];
        result.v_[1] = v_[1] * x.v_[1];
        result.v_[2] = v_[2] * x.v_[2];
        return result;
    }
    Vector3 operator-() { return Vector3(-v_[0], -v_[1], -v_[2]); }

    double& operator[](int i) { return v_[i]; }
    const double& operator[](int i) const { return v_[i]; }

    /// Checks for exact equality (i.e. no tolerances)
    bool operator==(const Vector3& RHS) const {
        return (v_[0] == RHS.v_[0] && v_[1] == RHS.v_[1] && v_[2] == RHS.v_[2]);
    }

    /// True if the coordinates differ in any bit
    bool operator!=(const Vector3& RHS) const { return !(*this == RHS); }

    double get(int i) {
        if (i >= 0 && i <= 2)
            return v_[i];
        else
            return 0.0;
    }

    void set(int i, double val) {
        if (i >= 0 && i <= 2) v_[i] = val;
    }

    double dot(const Vector3& x) const { return v_[0] * x.v_[0] + v_[1] * x.v_[1] + v_[2] * x.v_[2]; }

    double distance(const Vector3& s) const {
        double x = v_[0] - s.v_[0];
        double y = v_[1] - s.v_[1];
        double z = v_[2] - s.v_[2];
        return std::sqrt(x * x + y * y + z * z);
    }
    void normalize() {
        double temp = 0.0;

        for (auto i = 0; i < 3; ++i) temp += v_[i] * v_[i];
        temp = 1.0 / std::sqrt(temp);
        for (auto i = 0; i < 3; ++i) v_[i] *= temp;
    }
    double norm() const { return sqrt(this->dot(*this)); }
    void rotate(double theta, Vector3& axis) {
        Vector3 result;
        Vector3 unitaxis = axis;
        unitaxis.normalize();

        // split into parallel and perpendicular components along axis
        Vector3 parallel = axis * (this->dot(axis) / axis.dot(axis));
        Vector3 perpendicular = (*this) - parallel;

        // form unit vector perpendicular to parallel and perpendicular
        Vector3 third_axis = axis.perp_unit(perpendicular);
        third_axis = third_axis * perpendicular.norm();

        result = parallel + cos(theta) * perpendicular + sin(theta) * third_axis;
        (*this) = result;
    }
    Vector3 perp_unit(const Vector3& v) const {
        // try cross product
        Vector3 result = cross(v);
        double resultdotresult = result.dot(result);

        if (resultdotresult < 1.e-16) {
            // cross product is too small to normalize
            // find the largest of this and v
            double dotprodt = this->dot(*this);
            double dotprodv = v.dot(v);
            const Vector3* d;
            double dotprodd;
            if (dotprodt < dotprodv) {
                d = &v;
                dotprodd = dotprodv;
            } else {
                d = this;
                dotprodd = dotprodt;
            }

            // see if d is big enough
            if (dotprodd < 1.e-16) {
                // choose an arbitrary vector, since the biggest vector is small
                result[0] = 1.0;
                result[1] = 0.0;
                result[2] = 0.0;
                return result;
            } else {
                // choose a vector prependicular to d
                // choose it in one of the planes xy, xz, yz
                // choose the plane to be that which contains the two largest
                // components of d
                double absd[3];
                absd[0] = std::fabs(d->v_[0]);
                absd[1] = std::fabs(d->v_[1]);
                absd[2] = std::fabs(d->v_[2]);
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
                result[axis0] = d->v_[axis1];
                result[axis1] = -d->v_[axis0];
            }
            result.normalize();
            return result;
        } else {
            // normalize the cross product and return the result
            result *= 1.0 / sqrt(resultdotresult);
            return result;
        }
    }

    Vector3 cross(const Vector3& x) const {
        Vector3 result(v_[1] * x.v_[2] - v_[2] * x.v_[1], v_[2] * x.v_[0] - v_[0] * x.v_[2],
                       v_[0] * x.v_[1] - v_[1] * x.v_[0]);
        return result;
    }

    std::string to_string() const {
        std::stringstream s;
        s << "[ " << v_[0] << ", " << v_[1] << ", " << v_[2] << " ]";
        return s.str();
    }
};

inline Vector3 operator*(double d, const Vector3& x) {
    Vector3 result;
    result[0] = d * x[0];
    result[1] = d * x[1];
    result[2] = d * x[2];
    return result;
}
}  // namespace psi
