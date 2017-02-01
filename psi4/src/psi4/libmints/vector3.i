/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
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
 * @END LICENSE
 */

#include <cmath>

inline Vector3 operator*(double d, const Vector3& x)
{
    Vector3 result;
    result[0] = d * x[0];
    result[1] = d * x[1];
    result[2] = d * x[2];
    return result;
}

inline Vector3 Vector3::operator*(double d) const
{
    return d*(*this);
}

inline Vector3 Vector3::operator/(double d) const
{
    Vector3 result;
    result[0] = v_[0] / d;
    result[1] = v_[1] / d;
    result[2] = v_[2] / d;
    return result;
}

inline double Vector3::distance(const Vector3& s) const
{
    double x = v_[0] - s.v_[0];
    double y = v_[1] - s.v_[1];
    double z = v_[2] - s.v_[2];
    return sqrt(x*x + y*y + z*z);
}

inline void Vector3::normalize()
{
    double temp=0.0;
    int i;

    for (i=0; i<3; ++i)
        temp += v_[i] * v_[i];
    temp = 1.0 / sqrt(temp);
    for (i=0; i<3; ++i)
        v_[i] *= temp;
}

inline Vector3 Vector3::cross(const Vector3& x) const
{
    Vector3 result(v_[1] * x.v_[2] - v_[2] * x.v_[1],
                   v_[2] * x.v_[0] - v_[0] * x.v_[2],
                   v_[0] * x.v_[1] - v_[1] * x.v_[0]);
    return result;
}

inline void Vector3::rotate(double theta, Vector3& axis)
{
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

inline Vector3 Vector3::perp_unit(const Vector3& v) const
{
    // try cross product
    Vector3 result = cross(v);
    double resultdotresult = result.dot(result);

    if (resultdotresult < 1.e-16) {
        // cross product is too small to normalize
        // find the largest of this and v
        double dotprodt = this->dot(*this);
        double dotprodv = v.dot(v);
        const Vector3 *d;
        double dotprodd;
        if (dotprodt < dotprodv) {
            d = &v;
            dotprodd = dotprodv;
        }
        else {
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
        }
        else {
            // choose a vector prependicular to d
            // choose it in one of the planes xy, xz, yz
            // choose the plane to be that which contains the two largest
            // components of d
            double absd[3];
            absd[0] = fabs(d->v_[0]);
            absd[1] = fabs(d->v_[1]);
            absd[2] = fabs(d->v_[2]);
            int axis0, axis1;
            if ((absd[1] - absd[0]) > 1.0e-12) {
                axis0 = 1;
                if ((absd[2] - absd[0]) > 1.0e-12) {
                    axis1 = 2;
                }
                else {
                    axis1 = 0;
                }
            }
            else {
                axis0 = 0;
                if ((absd[2] - absd[1]) > 1.0e-12) {
                    axis1 = 2;
                }
                else {
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
    }
    else {
        // normalize the cross product and return the result
        result *= 1.0/sqrt(resultdotresult);
        return result;
    }
}
