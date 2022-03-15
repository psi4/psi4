/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2021 The Psi4 Developers.
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

#ifndef _psi_src_lib_libmints__vector3_h
#define _psi_src_lib_libmints__vector3_h

#include <array>
#include <cmath>
#include <sstream>

#include "psi4/pragma.h"

namespace psi {

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

    Vector3 operator*(double d) const;
    Vector3 operator/(double d) const;
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

    double distance(const Vector3&) const;
    void normalize();
    double norm() const { return sqrt(this->dot(*this)); }
    void rotate(double theta, Vector3& v);
    Vector3 perp_unit(const Vector3& v) const;

    Vector3 cross(const Vector3&) const;

    std::string to_string() const {
        std::stringstream s;
        s << "[ " << v_[0] << ", " << v_[1] << ", " << v_[2] << " ]";
        return s.str();
    }
};

Vector3 operator*(double, const Vector3&);

#include "vector3.i"

}  // namespace psi

#endif
