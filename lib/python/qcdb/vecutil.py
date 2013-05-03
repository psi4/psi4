#
#@BEGIN LICENSE
#
# PSI4: an ab initio quantum chemistry software package
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
#@END LICENSE
#

r"""File for accessory procedures in the chem module.
Credit for the libmints vector3 class to Justin M. Turney and
incremental improvements by other psi4 developers.

"""
import math
import copy
from exceptions import *

ZERO = 1.0E-14


def norm(v):
    """Compute the magnitude of vector *v*."""
    return math.sqrt(sum(v[i] * v[i] for i in range(len(v))))


def add(v, u):
    """Compute sum of vectors *v* and *u*."""
    return [u[i] + v[i] for i in range(len(v))]


def sub(v, u):
    """Compute difference of vectors *v* - *u*."""
    return [v[i] - u[i] for i in range(len(v))]


def dot(v, u):
    """Compute dot product of vectors *v* and *u*."""
    return sum(u[i] * v[i] for i in range(len(v)))


def scale(v, d):
    """Compute by-element scale by *d* of vector *v*."""
    return [d * v[i] for i in range(len(v))]


def normalize(v):
    """Compute normalized vector *v*."""
    vmag = norm(v)
    return [v[i] / vmag  for i in range(len(v))]


def distance(v, u):
    """Compute the distance between points defined by vectors *v* and *u*."""
    return norm(sub(v, u))


def cross(v, u):
    """Compute cross product of length 3 vectors *v* and *u*."""
    if len(u) != 3 or len(v) != 3:
        raise ValidationError('cross() only defined for vectors of length 3\n')
    return [v[1] * u[2] - v[2] * u[1],
            v[2] * u[0] - v[0] * u[2],
            v[0] * u[1] - v[1] * u[0]]


def rotate(v, theta, axis):
    """Rotate length 3 vector *v* about *axis* by *theta* radians."""
    if len(v) != 3 or len(axis) != 3:
        raise ValidationError('rotate() only defined for vectors of length 3\n')

    unitaxis = normalize(copy.deepcopy(axis))
    # split into parallel and perpendicular components along axis
    parallel = scale(axis, dot(v, axis) / dot(axis, axis))
    perpendicular = sub(v, parallel)
    # form unit vector perpendicular to parallel and perpendicular
    third_axis = perp_unit(axis, perpendicular)
    third_axis = scale(third_axis, norm(perpendicular))

    result = add(parallel, add(scale(perpendicular, math.cos(theta)), scale(third_axis, math.sin(theta))))
    for item in range(len(result)):
        if math.fabs(result[item]) < ZERO:
            result[item] = 0.0
    return result


def perp_unit(u, v):
    """Compute unit vector perpendicular to length 3 vectors *u* and *v*."""
    if len(u) != 3 or len(v) != 3:
        raise ValidationError('perp_unit() only defined for vectors of length 3\n')

    # try cross product
    result = cross(u, v)
    resultdotresult = dot(result, result)

    if resultdotresult < 1.E-16:
        # cross product is too small to normalize
        # find the largest of this and v
        dotprodt = dot(u, u)
        dotprodv = dot(v, v)
        if dotprodt < dotprodv:
            d = copy.deepcopy(v)
            dotprodd = dotprodv
        else:
            d = copy.deepcopy(u)
            dotprodd = dotprodt

        # see if d is big enough
        if dotprodd < 1.e-16:
            # choose an arbitrary vector, since the biggest vector is small
            result = [1.0, 0.0, 0.0]
            return result
        else:
            # choose a vector perpendicular to d
            # choose it in one of the planes xy, xz, yz
            # choose the plane to be that which contains the two largest components of d
            absd = [math.fabs(d[0]), math.fabs(d[1]), math.fabs(d[2])]
            if (absd[1] - absd[0]) > 1.0e-12:
            #if absd[0] < absd[1]:
                axis0 = 1
                if (absd[2] - absd[0]) > 1.0e-12:
                #if absd[0] < absd[2]:
                    axis1 = 2
                else:
                    axis1 = 0
            else:
                axis0 = 0
                if (absd[2] - absd[1]) > 1.0e-12:
                #if absd[1] < absd[2]:
                    axis1 = 2
                else:
                    axis1 = 1
            result = [0.0, 0.0, 0.0]
            # do the pi/2 rotation in the plane
            result[axis0] = d[axis1]
            result[axis1] = -1.0 * d[axis0]
        result = normalize(result)
        return result

    else:
        # normalize the cross product and return the result
        result = scale(result, 1.0 / math.sqrt(resultdotresult))
        return result


def determinant(mat):
    """Given 3x3 matrix *mat*, compute the determinat

    """
    if len(mat) != 3 or len(mat[0]) != 3 or len(mat[1]) != 3 or len(mat[2]) != 3:
        raise ValidationError('determinant() only defined for arrays of dimension 3x3\n')

    det = mat[0][0] * mat[1][1] * mat[2][2] - mat[0][2] * mat[1][1] * mat[2][0] + \
          mat[0][1] * mat[1][2] * mat[2][0] - mat[0][1] * mat[1][0] * mat[2][2] + \
          mat[0][2] * mat[1][0] * mat[2][1] - mat[0][0] * mat[1][2] * mat[2][1]
    return det


def diagonalize3x3symmat(A):
    """Given an real symmetric 3x3 matrix A, compute the eigenvalues

    """
    if len(A) != 3 or len(A[0]) != 3 or len(A[1]) != 3 or len(A[2]) != 3:
        raise ValidationError('diagonalize3x3symmat() only defined for arrays of dimension 3x3\n')

    B = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]  # I is the identity matrix
    p = A[0][1] * A[0][1] + A[0][2] * A[0][2] + A[1][2] * A[1][2]
    if p == 0:
        # A is diagonal
        eig1 = A[0][0]
        eig2 = A[1][1]
        eig3 = A[2][2]
    else:
        q = (A[0][0] + A[1][1] + A[2][2]) / 3.0
        p = (A[0][0] - q) ** 2 + (A[1][1] - q) ** 2 + (A[2][2] - q) ** 2 + 2.0 * p
        p = math.sqrt(p / 6.0)
        for i in range(3):
            for j in range(3):
                B[i][j] = (A[i][j] - q * B[i][j]) / p
        # B = (A - q * I) / p
        r = determinant(B) / 2.0

        # In exact arithmetric for a symmetric matrix  -1 <= r <= 1
        #    but computation error can leave it slightly outside this range.
        if r <= -1:
            phi = math.pi / 3
        elif r >= 1:
            phi = 0
        else:
            phi = math.acos(r) / 3.0

        # the eigenvalues satisfy eig3 <= eig2 <= eig1
        eig1 = q + 2.0 * p * math.cos(phi)
        eig3 = q + 2.0 * p * math.cos(phi + math.pi * (2.0 / 3.0))
        eig2 = 3.0 * q - eig1 - eig3     # since trace(A) = eig1 + eig2 + eig3
    return [eig1, eig2, eig3]


def zero(m, n):
    """ Create zero matrix"""
    new_matrix = [[0 for row in range(n)] for col in range(m)]
    return new_matrix


def show(matrix):
    """ Print out matrix"""
    for col in matrix:
        print col


def mscale(matrix, d):
    """Return *matrix* scaled by scalar *d*"""
    for i in range(len(matrix)):
        for j in range(len(matrix[0])):
            matrix[i][j] *= d
    return matrix


def mult(matrix1, matrix2):
    """ Matrix multiplication"""
    if len(matrix1[0]) != len(matrix2):
        # Check matrix dimensions
        raise ValidationError('Matrices must be m*n and n*p to multiply!')

    else:
        # Multiply if correct dimensions
        new_matrix = zero(len(matrix1), len(matrix2[0]))
        for i in range(len(matrix1)):
            for j in range(len(matrix2[0])):
                for k in range(len(matrix2)):
                    new_matrix[i][j] += matrix1[i][k] * matrix2[k][j]
        return new_matrix


def transpose(matrix):
    """Return matrix transpose"""
    if len(matrix[0]) != len(matrix):
        # Check matrix dimensions
        raise ValidationError('Matrices must be square.')

    tmat = [list(i) for i in zip(*matrix)]
    return tmat
