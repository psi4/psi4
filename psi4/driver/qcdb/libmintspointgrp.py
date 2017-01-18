#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2017 The Psi4 Developers.
#
# The copyrights for code used from other parties are included in
# the corresponding files.
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
# @END LICENSE
#

from __future__ import absolute_import
from __future__ import print_function
from .exceptions import *
from .vecutil import *
import sys
if sys.version_info >= (3,0):
    basestring = str


#
# Additional modifications made by Justin Turney <jturney@ccqc.uga.edu>
# for use in PSI4.
#
# Modifications are
# Copyright (C) 1996 Limit Point Systems, Inc.
#
# Author: Edward Seidl <seidl@janed.com>
# Maintainer: LPS
#
# This file is part of the SC Toolkit.
#
# The SC Toolkit is free software; you can redistribute it and/or modify
# it under the terms of the GNU Library General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.
#
# The SC Toolkit is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Library General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
# The U.S. Government is granted a limited license as per AL 91-7.
#
#
# pointgrp.h -- definition of the point group classes
#
#      THIS SOFTWARE FITS THE DESCRIPTION IN THE U.S. COPYRIGHT ACT OF A
#      "UNITED STATES GOVERNMENT WORK".  IT WAS WRITTEN AS A PART OF THE
#      AUTHOR'S OFFICIAL DUTIES AS A GOVERNMENT EMPLOYEE.  THIS MEANS IT
#      CANNOT BE COPYRIGHTED.  THIS SOFTWARE IS FREELY AVAILABLE TO THE
#      PUBLIC FOR USE WITHOUT A COPYRIGHT NOTICE, AND THERE ARE NO
#      RESTRICTIONS ON ITS USE, NOW OR SUBSEQUENTLY.
#
#  Author:
#      E. T. Seidl
#      Bldg. 12A, Rm. 2033
#      Computer Systems Laboratory
#      Division of Computer Research and Technology
#      National Institutes of Health
#      Bethesda, Maryland 20892
#      Internet: seidl@alw.nih.gov
#      June, 1993
#


SymmOps = {'E': 0,
           'C2_z': 1,
           'C2_y': 2,
           'C2_x': 4,
           'i': 8,
           'Sigma_xy': 16,
           'Sigma_xz': 32,
           'Sigma_yz': 64,
           'ID': 128
           }

PointGroups = {
    'C1': SymmOps['E'],
    'Ci': SymmOps['E'] | SymmOps['i'],
    'C2X': SymmOps['E'] | SymmOps['C2_x'],
    'C2Y': SymmOps['E'] | SymmOps['C2_y'],
    'C2Z': SymmOps['E'] | SymmOps['C2_z'],
    'CsZ': SymmOps['E'] | SymmOps['Sigma_xy'],
    'CsY': SymmOps['E'] | SymmOps['Sigma_xz'],
    'CsX': SymmOps['E'] | SymmOps['Sigma_yz'],
    'D2': SymmOps['E'] | SymmOps['C2_x'] | SymmOps['C2_y'] | SymmOps['C2_z'],
    'C2vX': SymmOps['E'] | SymmOps['C2_x'] | SymmOps['Sigma_xy'] | SymmOps['Sigma_xz'],
    'C2vY': SymmOps['E'] | SymmOps['C2_y'] | SymmOps['Sigma_xy'] | SymmOps['Sigma_yz'],
    'C2vZ': SymmOps['E'] | SymmOps['C2_z'] | SymmOps['Sigma_xz'] | SymmOps['Sigma_yz'],
    'C2hX': SymmOps['E'] | SymmOps['C2_x'] | SymmOps['Sigma_yz'] | SymmOps['i'],
    'C2hY': SymmOps['E'] | SymmOps['C2_y'] | SymmOps['Sigma_xz'] | SymmOps['i'],
    'C2hZ': SymmOps['E'] | SymmOps['C2_z'] | SymmOps['Sigma_xy'] | SymmOps['i'],
    'D2h': SymmOps['E'] | SymmOps['C2_x'] | SymmOps['C2_y'] | SymmOps['C2_z'] | SymmOps['i'] | \
            SymmOps['Sigma_xy'] | SymmOps['Sigma_xz'] | SymmOps['Sigma_yz']
    }


# changed signature from def similar(bits, sim, cnt):
def similar(bits):
    """From *bits* of a directionalized point group, returns array of
    bits of all directions.

    """
    cs = [PointGroups['CsX'], PointGroups['CsY'], PointGroups['CsZ']]
    c2v = [PointGroups['C2vZ'], PointGroups['C2vY'], PointGroups['C2vX']]
    c2h = [PointGroups['C2hZ'], PointGroups['C2hY'], PointGroups['C2hX']]
    c2 = [PointGroups['C2Z'], PointGroups['C2Y'], PointGroups['C2X']]
    d2h = [PointGroups['D2h']]
    d2 = [PointGroups['D2']]
    ci = [PointGroups['Ci']]
    c1 = [PointGroups['C1']]

    if bits in cs:
        sim = cs
    elif bits in c2v:
        sim = c2v
    elif bits in c2h:
        sim = c2h
    elif bits in c2:
        sim = c2
    elif bits in d2h:
        sim = d2h
    elif bits in ci:
        sim = ci
    elif bits in c1:
        sim = c1
    elif bits in d2:
        sim = d2
    else:
        raise ValidationError('PointGroups::similar: Should not have reached here.')

    return sim, len(sim)


class SymmetryOperation(object):
    """The SymmetryOperation class provides a 3 by 3 matrix
    representation of a symmetry operation, such as a rotation or reflection.

    """

    def __init__(self, *args):
        """Constructor"""

        # matrix representation
        self.d = zero(3, 3)
        # bits representation
        self.bits = 0

        # Divert to constructor functions
        if len(args) == 0:
            pass
        elif len(args) == 1 and \
            isinstance(args[0], SymmetryOperation):
            self.constructor_symmetryoperation(*args)
        else:
            raise ValidationError('SymmetryOperation::constructor: Inappropriate configuration of constructor arguments')

    # <<< Methods for Construction >>>

    def constructor_symmetryoperation(self, so):
        self.bits = so.bits
        self.d = [row[:] for row in so.d]

    # <<< Simple Methods for Basic SymmetryOperation Information >>>

    def bit(self):
        """Get the bit value."""
        return self.bits

    def __getitem__(self, i, j=None):
        """Returns the (i,j)th element of the transformation matrix
        or the i'th row of the transformation matrix if *j* is None.

        """
        if j is None:
            return self.d[i]
        else:
            return self.d[i][j]

    def trace(self):
        """returns the trace of the transformation matrix"""
        return self.d[0][0] + self.d[1][1] + self.d[2][2]

    # <<< Methods for Symmetry Operations >>>

    def zero(self):
        """zero out the symop"""
        self.d = zero(3, 3)

    def unit(self):
        """Set equal to a unit matrix"""
        self.zero()
        self.d[0][0] = 1.0
        self.d[1][1] = 1.0
        self.d[2][2] = 1.0

    def E(self):
        """Set equal to E"""
        self.unit()
        self.bits = SymmOps['E']

    def i(self):
        """Set equal to an inversion"""
        self.zero()
        self.d[0][0] = -1.0
        self.d[1][1] = -1.0
        self.d[2][2] = -1.0
        self.bits = SymmOps['i']

    def sigma_xy(self):
        """Set equal to reflection in xy plane"""
        self.unit()
        self.d[2][2] = -1.0
        self.bits = SymmOps['Sigma_xy']

    def sigma_xz(self):
        """Set equal to reflection in xz plane"""
        self.unit()
        self.d[1][1] = -1.0
        self.bits = SymmOps['Sigma_xz']

    def sigma_yz(self):
        """Set equal to reflection in yz plane"""
        self.unit()
        self.d[0][0] = -1.0
        self.bits = SymmOps['Sigma_yz']

    def c2_x(self):
        """Set equal to C2 about the x axis"""
        self.i()
        self.d[0][0] = 1.0
        self.bits = SymmOps['C2_x']

    def c2_y(self):
        """Set equal to C2 about the y axis"""
        self.i()
        self.d[1][1] = 1.0
        self.bits = SymmOps['C2_y']

    def c2_z(self):
        """Set equal to C2 about the z axis"""
        self.i()
        self.d[2][2] = 1.0
        self.bits = SymmOps['C2_z']

    # <<< Methods for Operations >>>

    def analyze_d(self):
        """

        """
        temp = [self.d[0][0], self.d[1][1], self.d[2][2]]
        tol = 1.0e-5

        if all([abs(temp[idx] - val) < tol for idx, val in enumerate([1.0, 1.0, 1.0])]):
            self.bits = SymmOps['E']
        elif all([abs(temp[idx] - val) < tol for idx, val in enumerate([1.0, -1.0, -1.0])]):
            self.bits = SymmOps['C2_x']
        elif all([abs(temp[idx] - val) < tol for idx, val in enumerate([-1.0, 1.0, -1.0])]):
            self.bits = SymmOps['C2_y']
        elif all([abs(temp[idx] - val) < tol for idx, val in enumerate([-1.0, -1.0, 1.0])]):
            self.bits = SymmOps['C2_z']
        elif all([abs(temp[idx] - val) < tol for idx, val in enumerate([1.0, 1.0, -1.0])]):
            self.bits = SymmOps['Sigma_xy']
        elif all([abs(temp[idx] - val) < tol for idx, val in enumerate([1.0, -1.0, 1.0])]):
            self.bits = SymmOps['Sigma_xz']
        elif all([abs(temp[idx] - val) < tol for idx, val in enumerate([-1.0, 1.0, 1.0])]):
            self.bits = SymmOps['Sigma_yz']
        elif all([abs(temp[idx] - val) < tol for idx, val in enumerate([-1.0, -1.0, -1.0])]):
            self.bits = SymmOps['i']

    def operate(self, r):
        """This operates on this with r (i.e. return r * this)"""

        ret = SymmetryOperation()
        for i in range(3):
            for j in range(3):
                t = 0.0
                for k in range(3):
                    t += r.d[i][k] * self.d[k][j]
                ret.d[i][j] = t
        ret.analyze_d()
        return ret

    def transform(self, r):
        """This performs the transform r * this * r~"""

        # foo = r * d
        foo = SymmetryOperation()
        for i in range(3):
            for j in range(3):
                t = 0.0
                for k in range(3):
                    t += r.d[i][k] * self.d[k][j]
                foo.d[i][j] = t

        # ret = (r*d)*r~ = foo*r~
        ret = SymmetryOperation()
        for i in range(3):
            for j in range(3):
                t = 0.0
                for k in range(3):
                    t += foo.d[i][k] * r.d[j][k]
                ret.d[i][j] = t

        ret.analyze_d()
        return ret

#    SymmetryOperation & operator = (SymmetryOperation const & a); // Assignment operator

    def rotation(self, theta):
        """Set equal to a clockwise rotation by 2pi/n or theta degrees"""

        if isinstance(theta, int):
            theta = 2.0 * math.pi if theta == 0 else 2.0 * math.pi / theta
        ctheta = math.cos(theta)
        stheta = math.sin(theta)

        self.zero()
        self.d[0][0] = ctheta
        self.d[0][1] = stheta
        self.d[1][0] = -stheta
        self.d[1][1] = ctheta
        self.d[2][2] = 1.0
        self.analyze_d()

    def transpose(self):
        """Transpose matrix operation"""

        for i in range(3):
            for j in range(i):
                tmp = self.d[i][j]
                self.d[i][j] = self.d[j][i]
                self.d[j][i] = tmp
        self.analyze_d()

    # <<< Methods for Printing >>>

    def __str__(self, out=None):
        """print the matrix"""

        text = "        1          2          3\n"
        text += "  1  "
        text += "%10.7f " % (self.d[0][0])
        text += "%10.7f " % (self.d[0][1])
        text += "%10.7f \n" % (self.d[0][2])
        text += "  2  "
        text += "%10.7f " % (self.d[1][0])
        text += "%10.7f " % (self.d[1][1])
        text += "%10.7f \n" % (self.d[1][2])
        text += "  3  "
        text += "%10.7f " % (self.d[2][0])
        text += "%10.7f " % (self.d[2][1])
        text += "%10.7f \n" % (self.d[2][2])
        text += "bits_ = %d\n" % (self.bits)

        if out is None:
            return text
        else:
            with open(out, mode='w') as handle:
                handle.write(text)


class SymRep(object):
    """The SymRep class provides an n dimensional matrix representation of a
    symmetry operation, such as a rotation or reflection. The trace of a
    SymRep can be used as the character for that symmetry operation. d is
    hardwired to 5x5 since the H irrep in Ih is 5 dimensional.

    """

    def __init__(self, *args):
        """Constructor"""

        # order of representation
        self.n = 0
        # matrix representation
        self.d = zero(5, 5)

        # Divert to constructor functions
        if len(args) == 1 and \
            isinstance(args[0], int):
            self.constructor_order(*args)
        elif len(args) == 1 and \
            isinstance(args[0], SymmetryOperation):
            self.constructor_symmetryoperation(*args)
        else:
            raise ValidationError('SymRep::constructor: Inappropriate configuration of constructor arguments')

    # <<< Methods for Construction >>>

    def constructor_order(self, i):
        """Initialize order only

        """
        self.n = i
        self.zero()

    def constructor_symmetryoperation(self, so):
        """Initialize from 3x3 SymmetryOperation

        """
        self.n = 3
        self.zero()
        for i in range(3):
            for j in range(3):
                self.d[i][j] = so[i][j]

    def SymmetryOperation(self):
        """Cast SymRep to SymmetryOperation

        """
        if self.n != 3:
            raise ValidationError("SymRep::operator SymmetryOperation(): trying to cast to symop when n != 3")

        so = SymmetryOperation()
        for i in range(3):
            for j in range(3):
                so[i][j] = self.d[i][j]
        return so

    # <<< Simple Methods for Basic SymRep Information >>>

    def set_dim(self, i):
        """Set the dimension of d"""
        self.n = i

    def __getitem__(self, i, j=None):
        """Returns the (i,j)th element of the transformation matrix
        or the i'th row of the transformation matrix if *j* is None.

        """
        if j is None:
            return self.d[i]
        else:
            return self.d[i][j]

    def trace(self):
        """returns the trace of the transformation matrix

        """
        r = 0.0
        for i in range(self.n):
            r += self.d[i][i]
        return r

    # <<< Methods for Symmetry Operations >>>

    def zero(self):
        """zero out the symop"""
        self.d = zero(5, 5)

    def unit(self):
        """Set equal to a unit matrix"""
        self.zero()
        self.d[0][0] = 1.0
        self.d[1][1] = 1.0
        self.d[2][2] = 1.0
        self.d[3][3] = 1.0
        self.d[4][4] = 1.0

    def E(self):
        """Set equal to the identity"""
        self.unit()

    def i(self):
        """Set equal to an inversion"""
        self.zero()
        self.d[0][0] = -1.0
        self.d[1][1] = -1.0
        self.d[2][2] = -1.0
        self.d[3][3] = -1.0
        self.d[4][4] = -1.0

    def sigma_h(self):
        """Set equal to reflection in xy plane

        """
        self.unit()
        if self.n == 3:
            self.d[2][2] = -1.0
        elif self.n == 5:
            self.d[3][3] = -1.0
            self.d[4][4] = -1.0

    def sigma_xz(self):
        """Set equal to reflection in xz plane

        """
        self.unit()
        if self.n == 2 or self.n == 3 or self.n == 4:
            self.d[1][1] = -1.0
            if self.n == 4:
                self.d[2][2] = -1.0
        elif self.n == 5:
            self.d[2][2] = -1.0
            self.d[4][4] = -1.0

    def sigma_yz(self):
        """Set equal to reflection in yz plane

        """
        self.unit()
        if self.n == 2 or self.n == 3 or self.n == 4:
            self.d[0][0] = -1.0
            if self.n == 4:
                self.d[3][3] = -1.0
        elif self.n == 5:
            self.d[2][2] = -1.0
            self.d[3][3] = -1.0

    def c2_x(self):
        """Set equal to C2 about the x axis

        """
        self.i()
        if self.n == 2 or self.n == 3 or self.n == 4:
            self.d[0][0] = 1.0
            if self.n == 4:
                self.d[3][3] = 1.0
        elif self.n == 5:
            self.d[0][0] = 1.0
            self.d[1][1] = 1.0
            self.d[4][4] = 1.0

    def c2_y(self):
        """Set equal to C2 about the y axis

        """
        self.i()
        if self.n == 2 or self.n == 3 or self.n == 4:
            self.d[1][1] = 1.0
            if self.n == 4:
                self.d[2][2] = 1.0
        elif self.n == 5:
            self.d[0][0] = 1.0
            self.d[1][1] = 1.0
            self.d[3][3] = 1.0

    def c2_z(self):
        """Set equal to C2 about the z axis

        """
        self.i()
        if self.n == 2 or self.n == 3 or self.n == 4:
            self.d[1][1] = 1.0
            if self.n == 4:
                self.d[2][2] = 1.0
        elif self.n == 5:
            self.d[0][0] = 1.0
            self.d[1][1] = 1.0
            self.d[3][3] = 1.0

    # <<< Methods for Operations >>>

    def operate(self, r):
        """This operates on this with r (i.e. return r * this)

        """
        if r.n != self.n:
            raise ValidationError("SymRep::operate(): dimensions don't match")

        ret = SymRep(self.n)
        for i in range(self.n):
            for j in range(self.n):
                t = 0.0
                for k in range(self.n):
                    t += r[i][k] * self.d[k][j]
                ret[i][j] = t
        return ret

    def transform(self, r):
        """This performs the transform r * this * r~

        """
        if r.n != self.n:
            raise ValidationError("SymRep::operate(): dimensions don't match")

        foo = SymRep(n)
        # foo = r * d
        for i in range(self.n):
            for j in range(self.n):
                t = 0.0
                for k in range(self.n):
                    t += r[i][k] * d[k][j]
                foo[i][j] = t

        ret = SymRep(n)
        # ret = (r*d)*r~ = foo*r~
        for i in range(self.n):
            for j in range(self.n):
                t = 0.0
                for k in range(self.n):
                    t += foo[i][k] * r[j][k]
                ret[i][j] = t

        return ret

    def rotation(self, theta):
        """Set equal to a clockwise rotation by 2pi/n or theta degrees

        """
        if isinstance(theta, int):
            theta = 2.0 * math.pi if theta == 0 else 2.0 * math.pi / theta

        ctheta = math.cos(theta)
        stheta = math.sin(theta)
        c2theta = math.cos(2 * theta)
        s2theta = math.sin(2 * theta)

        self.zero()
        if self.n == 1:
            self.d[0][0] = 1.0

        elif self.n == 3:
            self.d[0][0] = ctheta
            self.d[0][1] = stheta
            self.d[1][0] = -stheta
            self.d[1][1] = ctheta
            self.d[2][2] = 1.0

        elif self.n == 2 or self.n == 4:
            self.d[0][0] = ctheta
            self.d[0][1] = stheta
            self.d[1][0] = -stheta
            self.d[1][1] = ctheta

            # this is ok since d is hardwired
            self.d[2][2] = c2theta
            self.d[2][3] = -s2theta
            self.d[3][2] = s2theta
            self.d[3][3] = c2theta

        elif self.n == 5:
            self.d[0][0] = 1.0

            self.d[1][1] = c2theta
            self.d[1][2] = s2theta
            self.d[2][1] = -s2theta
            self.d[2][2] = c2theta

            self.d[3][3] = ctheta
            self.d[3][4] = -stheta
            self.d[4][3] = stheta
            self.d[4][4] = ctheta

        else:
            raise ValidationError("SymRep::rotation(): n > 5")


class IrreducibleRepresentation(object):
    """The IrreducibleRepresentation class provides information associated
    with a particular irreducible representation of a point group. This
    includes the Mulliken symbol for the irrep, the degeneracy of the
    irrep, the characters which represent the irrep, and the number of
    translations and rotations in the irrep. The order of the point group
    is also provided (this is equal to the number of characters in an
    irrep).

    """

    def __init__(self, *args):
        """Constructor"""

        # the order of the group
        self.g = 0  # int  really self?
        # the degeneracy of the irrep
        self.degen = 0  # int  really self?
        # the number of rotations in this irrep
        self.PYnrot = 0  # int
        # the number of translations in this irrep
        self.PYntrans = 0  # int
        # true if this irrep has a complex representation
        self.PYcomplex = 0
        # mulliken symbol for this irrep
        self.symb = 0  # str
        # mulliken symbol for this irrep w/o special characters
        self.csymb = 0  # str
        # representation matrices for the symops
        self.rep = []

        # Divert to constructor functions
        if len(args) == 0:
            pass
        elif len(args) == 4 and \
            isinstance(args[0], int) and \
            isinstance(args[1], int) and \
            isinstance(args[2], basestring) and \
            isinstance(args[3], basestring):
            self.constructor_order_degen_mulliken(*args)
        else:
            raise ValidationError('IrreducibleRepresentation::constructor: Inappropriate configuration of constructor arguments')

    # <<< Methods for Construction >>>

    def constructor_order_degen_mulliken(self, order, d, lab, clab):
        """This constructor takes as arguments the *order* of the point
        group, the degeneracy *d* of the irrep, and the Mulliken symbol of
        the irrep. The Mulliken symbol is copied internally.

        """
        self.init(order, d, lab, clab)

    def init(self, order, d, lab, clab):
        """Initialize the order, degeneracy, and Mulliken symbol of the
        irrep.

        """
        self.g = order
        self.degen = d
        self.PYntrans = 0
        self.PYnrot = 0
        self.PYcomplex = 0
        self.symb = lab
        self.csymb = clab

        if order > 0:
            for i in range(order):
                self.rep.append(SymRep(d))

#    IrreducibleRepresentation(const IrreducibleRepresentation&);
#    IrreducibleRepresentation& operator=(const IrreducibleRepresentation&);

    # <<< Simple Methods for Basic IrreducibleRepresentation Information >>>

    def order(self):
        """Returns the order of the group."""
        return self.g

    def degeneracy(self):
        """Returns the degeneracy of the irrep."""
        return self.degen

    def complex(self):
        """Returns the value of complex"""
        return self.PYcomplex

    def nproj(self):
        """Returns the number of projection operators for the irrep."""
        return self.degen * self.degen

    def nrot(self):
        """Returns the number of rotations associated with the irrep."""
        return self.PYnrot

    def ntrans(self):
        """Returns the number of translations associated with the irrep."""
        return self.PYntrans

    def symbol(self):
        """Returns the Mulliken symbol for the irrep."""
        return self.symb

    def symbol_ns(self):
        """Returns the Mulliken symbol for the irrep without special
        characters.

        """
        return self.csymb if self.csymb else self.symb

    def character(self, i):
        """Returns the character for the i'th symmetry operation of the
        point group.

        """
        return 0.5 * self.rep[i].trace() if self.complex() else self.rep[i].trace()

    def p(self, x1, x2, i=None):
        """Returns the element (x1, x2) of the i'th representation matrix.
        Or Returns the character for the x1'th contribution to the x2'th
        representation matrix.

        """
        if i is None:
            dr = x1 % self.degen  # dr should be int; always seems to be
            dc = x1 / self.degen  # dc should be int; always seems to be
            #print 'need to be int', dr, dc
            i = x2
            x1 = dr
            x2 = dc

        return self.rep[i][x1][x2]

    # <<< Methods for Printing >>>

    def __str__(self, out=None):
        """This prints the irrep to the given file, or stdout if none is
        given. The second argument is an optional string of spaces to
        offset by.

        """
        if self.g == 0:
            return

        text = "  %-5s" % (self.symb)

        for i in range(self.g):
            text += " %6.3f" % (self.character(i))
        text += " | %d t, %d R\n" % (self.PYntrans, self.PYnrot)

        for d in range(self.nproj()):
            text += "       "
            for i in range(self.g):
                text += " %6.3f" % (self.p(d, i))
            text += "\n"

        if out is None:
            return text
        else:
            with open(out, mode='w') as handle:
                handle.write(text)


class CharacterTable(object):
    """The CharacterTable class provides a workable character table for
    all of the non-cubic point groups. While I have tried to match the
    ordering in Cotton's book, I don't guarantee that it is always
    followed. It shouldn't matter anyway. Also note that I don't lump
    symmetry operations of the same class together. For example, in C3v
    there are two distinct C3 rotations and 3 distinct reflections, each
    with a separate character. Thus symop has 6 elements rather than the 3
    you'll find in most published character tables.

    """

    def __init__(self, *args):
        """Constructor"""

        # order of the principal rot axis
        self.nt = 0
        # the class of the point group
        self.pg = PointGroups['C1']
        # the number of irreps in this pg
        self.PYnirrep = 0
        # an array of irreps
        self.PYgamma = 0
        # the matrices describing sym ops
        self.symop = 0
        # index of the inverse symop
        self.inv = 0
        # the Schoenflies symbol for the pg
        self.symb = 0
        # Bitwise representation of the symmetry operations
        self.PYbits = 0

        # Divert to constructor functions
        if len(args) == 0:
            pass
        elif len(args) == 1 and \
            isinstance(args[0], basestring):
            self.constructor_schoenflies(*args)
        elif len(args) == 1 and \
            isinstance(args[0], int):
            self.constructor_bits(*args)
        else:
            raise ValidationError('BasisSet::constructor: Inappropriate configuration of constructor arguments')

    # <<< Methods for Construction >>>

    def constructor_schoenflies(self, cpg):
        """This constructor takes the Schoenflies symbol of a point group
        as input.

        """
        self.symb = cpg
        # Check the symbol coming in
        self.PYbits = PointGroup.full_name_to_bits(cpg)
        if self.PYbits is None:
            raise ValidationError('CharacterTable: Invalid point group name: %s\n' % (cpg))
        self.common_init()

    def constructor_bits(self, bits):
        """This constructor takes the bitswise representation of a point
        group as input.

        """
        self.PYbits = bits
        self.symb = PointGroup.bits_to_basic_name(bits)
        self.common_init()

    def common_init(self):
        """First parse the point group symbol, this will give us the
        order of the point group(g), the type of point group (pg), the
        order of the principle rotation axis (nt), and the number of
        irreps (nirrep).

        """
        if len(self.symb) == 0:
            raise ValidationError('CharacterTable::CharacterTable: null point group')
        if self.make_table() < 0:
            raise ValidationError('CharacterTable::CharacterTable: could not make table')

#    CharacterTable(const CharacterTable&);
#CharacterTable::CharacterTable(const CharacterTable& ct)
#    : nt(0), pg(PointGroups::C1), nirrep_(0), gamma_(0), symop(0), _inv(0), symb(0),
#      bits_(0)
#{
#    *this = ct;
#}
#
#
#    CharacterTable& operator=(const CharacterTable&);
#CharacterTable&
#CharacterTable::operator=(const CharacterTable& ct)
#{
#    nt=ct.nt; pg=ct.pg; nirrep_=ct.nirrep_;
#
#    symb = ct.symb;
#
#    if (gamma_) delete[] gamma_; gamma_=0;
#    if (ct.gamma_) {
#        gamma_ = new IrreducibleRepresentation[nirrep_];
#        for (int i=0; i < nirrep_; i++) {
#            gamma_[i].init();
#            gamma_[i] = ct.gamma_[i];
#        }
#    }
#
#    if (symop)
#        delete[] symop;
#    symop=0;
#
#    if (ct.symop) {
#        symop = new SymmetryOperation[nirrep_];
#        for (int i=0; i < nirrep_; i++) {
#            symop[i] = ct.symop[i];
#        }
#    }
#
#    if (_inv)
#        delete[] _inv;
#    _inv=0;
#
#    if (ct._inv) {
#        _inv = new int[nirrep_];
#        memcpy(_inv,ct._inv,sizeof(int)* nirrep_);
#    }
#
#    return *this;
#}

    # <<< Simple Methods for Basic CharacterTable Information >>>

    def nirrep(self):
        """Returns the number of irreps."""
        return self.PYnirrep

    def order(self):
        """Returns the order of the point group"""
        return self.PYnirrep

    def symbol(self):
        """Returns the Schoenflies symbol for the point group"""
        return self.symb

    def bits(self):
        """Returns bitwise representation of symm ops"""
        return self.PYbits

    def complex(self):
        """Cn, Cnh, Sn, T, and Th point groups have complex representations.
        This function returns 1 if the point group has a complex
        representation, 0 otherwise.

        """
        return 0

    def gamma(self, i):
        """Returns the i'th irrep."""
        return self.PYgamma[i]

    def symm_operation(self, i):
        """Returns the i'th symmetry operation."""
        return self.symop[i]

    def inverse(self, i):
        """Returns the index of the symop which is the inverse of symop[i]."""
        return self.inv[i]

    def ncomp(self):
        """Returns number of compenents, including degeneracies

        """
        ret = 0
        for i in range(self.PYnirrep):
            nc = 1 if self.PYgamma[i].complex() else self.PYgamma[i].degen
            ret += nc
        return ret

    def which_irrep(self, i):
        """Returns the irrep component i belongs to.

        """
        cn = 0
        for ir in range(self.PYnirrep):
            nc = 1 if self.PYgamma[ir].complex() else self.PYgamma[ir].degen
            for c in range(nc):
                print('i =', i, 'ir =', ir, 'c =', c, 'cn =', cn, 'nc =', nc)
                if cn == i:
                    return ir
                cn += 1  # right place to increment?
        return -1

    def which_comp(self, i):
        """Returns which component i is.

        """
        cn = 0
        for ir in range(self.PYnirrep):
            nc = 1 if self.PYgamma[ir].complex() else self.PYgamma[ir].degen
            for c in range(nc):
                print('i =', i, 'ir =', ir, 'c =', c, 'cn =', cn, 'nc =', nc)
                if cn == i:
                    return c
                cn += 1  # right place to increment?
        return -1

    # <<< Methods for Operations >>>

    def make_table(self):
        """This function will generate a character table for the point
        group. This character table is in the order that symmetry
        operations are generated, not in Cotton order. If this is a
        problem, tough. Also generate the transformation matrices.
        This fills in the irrep and symop arrays

        """
        # set nt and nirrep
        if self.PYbits in [
            PointGroups['C1']]:
            self.PYnirrep = 1
            self.nt = 1

        elif self.PYbits in [
            PointGroups['CsX'],
            PointGroups['CsY'],
            PointGroups['CsZ'],
            PointGroups['Ci']]:
            self.PYnirrep = 2
            self.nt = 1

        elif self.PYbits in [
            PointGroups['C2X'],
            PointGroups['C2Y'],
            PointGroups['C2Z']]:
            self.PYnirrep = 2
            self.nt = 2

        elif self.PYbits in [
            PointGroups['C2hX'],
            PointGroups['C2hY'],
            PointGroups['C2hZ'],
            PointGroups['C2vX'],
            PointGroups['C2vY'],
            PointGroups['C2vZ'],
            PointGroups['D2']]:
            self.PYnirrep = 4
            self.nt = 2

        elif self.PYbits in [
            PointGroups['D2h']]:
            self.PYnirrep = 8
            self.nt = 2

        else:
            raise ValidationError("Should not have receached here!")

        if self.PYnirrep == 0:
            return 0

        so = SymmetryOperation()
        self.PYgamma = []
        self.symop = []
        self.inv = []
        for h in range(self.PYnirrep):
            self.PYgamma.append(IrreducibleRepresentation())
            self.symop.append(SymmetryOperation())
            self.inv.append(0)

        # this array forms a reducible representation for rotations about x,y,z
        rot = zero(self.PYnirrep, 1)

        # this array forms a reducible representation for translations along x,y,z
        trans = zero(self.PYnirrep, 1)

        # the angle to rotate about the principal axis
        theta = 2.0 * math.pi if self.nt == 0 else 2.0 * math.pi / self.nt

        # Handle irreducible representations; set PYgamma
        if self.PYbits in [
            PointGroups['C1']]:
            # no symmetry case
            self.PYgamma[0].init(1, 1, "A", "A")
            self.PYgamma[0].PYnrot = 3
            self.PYgamma[0].PYntrans = 3
            self.PYgamma[0].rep[0][0][0] = 1.0

        elif self.PYbits in [
            PointGroups['CsX'],  # reflection through the yz plane
            PointGroups['CsY'],  # reflection through the xz plane
            PointGroups['CsZ']]:  # reflection through the xy plane
            self.PYgamma[0].init(2, 1, "A'", "Ap")
            self.PYgamma[0].rep[0][0][0] = 1.0
            self.PYgamma[0].rep[1][0][0] = 1.0
            self.PYgamma[0].PYnrot = 1
            self.PYgamma[0].PYntrans = 2

            self.PYgamma[1].init(2, 1, "A\"", "App")
            self.PYgamma[1].rep[0][0][0] = 1.0
            self.PYgamma[1].rep[1][0][0] = -1.0
            self.PYgamma[1].PYnrot = 2
            self.PYgamma[1].PYntrans = 1

        elif self.PYbits in [
            PointGroups['Ci']]:
            # equivalent to S2 about the z axis
            self.PYgamma[0].init(2, 1, "Ag", "Ag")
            self.PYgamma[0].rep[0][0][0] = 1.0
            self.PYgamma[0].rep[1][0][0] = 1.0
            self.PYgamma[0].PYnrot = 3

            self.PYgamma[1].init(2, 1, "Au", "Au")
            self.PYgamma[1].rep[0][0][0] = 1.0
            self.PYgamma[1].rep[1][0][0] = -1.0
            self.PYgamma[1].PYntrans = 3

        elif self.PYbits in [
            PointGroups['C2X'],
            PointGroups['C2Y'],
            PointGroups['C2Z']]:
            self.PYgamma[0].init(2, 1, "A", "A")
            self.PYgamma[0].rep[0][0][0] = 1.0
            self.PYgamma[0].rep[1][0][0] = 1.0
            self.PYgamma[0].PYnrot = 1
            self.PYgamma[0].PYntrans = 1

            self.PYgamma[1].init(2, 1, "B", "B")
            self.PYgamma[1].rep[0][0][0] = 1.0
            self.PYgamma[1].rep[1][0][0] = -1.0
            self.PYgamma[1].PYnrot = 2
            self.PYgamma[1].PYntrans = 2

        elif self.PYbits in [
            PointGroups['C2hX'],
            PointGroups['C2hY'],
            PointGroups['C2hZ']]:
            self.PYgamma[0].init(4, 1, "Ag", "Ag")
            self.PYgamma[0].rep[0][0][0] = 1.0
            self.PYgamma[0].rep[1][0][0] = 1.0
            self.PYgamma[0].rep[2][0][0] = 1.0
            self.PYgamma[0].rep[3][0][0] = 1.0
            self.PYgamma[0].PYnrot = 1
            self.PYgamma[0].PYntrans = 0

            self.PYgamma[1].init(4, 1, "Bg", "Bg")
            self.PYgamma[1].rep[0][0][0] = 1.0
            self.PYgamma[1].rep[1][0][0] = -1.0
            self.PYgamma[1].rep[2][0][0] = 1.0
            self.PYgamma[1].rep[3][0][0] = -1.0
            self.PYgamma[1].PYnrot = 2
            self.PYgamma[1].PYntrans = 0

            self.PYgamma[2].init(4, 1, "Au", "Au")
            self.PYgamma[2].rep[0][0][0] = 1.0
            self.PYgamma[2].rep[1][0][0] = 1.0
            self.PYgamma[2].rep[2][0][0] = -1.0
            self.PYgamma[2].rep[3][0][0] = -1.0
            self.PYgamma[2].PYnrot = 0
            self.PYgamma[2].PYntrans = 1

            self.PYgamma[3].init(4, 1, "Bu", "Bu")
            self.PYgamma[3].rep[0][0][0] = 1.0
            self.PYgamma[3].rep[1][0][0] = -1.0
            self.PYgamma[3].rep[2][0][0] = -1.0
            self.PYgamma[3].rep[3][0][0] = 1.0
            self.PYgamma[3].PYnrot = 0
            self.PYgamma[3].PYntrans = 2

        elif self.PYbits in [
            PointGroups['C2vX'],
            PointGroups['C2vY'],
            PointGroups['C2vZ']]:
            self.PYgamma[0].init(4, 1, "A1", "A1")
            self.PYgamma[0].rep[0][0][0] = 1.0
            self.PYgamma[0].rep[1][0][0] = 1.0
            self.PYgamma[0].rep[2][0][0] = 1.0
            self.PYgamma[0].rep[3][0][0] = 1.0
            self.PYgamma[0].PYnrot = 0
            self.PYgamma[0].PYntrans = 1

            self.PYgamma[1].init(4, 1, "A2", "A2")
            self.PYgamma[1].rep[0][0][0] = 1.0
            self.PYgamma[1].rep[1][0][0] = 1.0
            self.PYgamma[1].rep[2][0][0] = -1.0
            self.PYgamma[1].rep[3][0][0] = -1.0
            self.PYgamma[1].PYnrot = 1
            self.PYgamma[1].PYntrans = 0

            self.PYgamma[2].init(4, 1, "B1", "B1")
            self.PYgamma[2].rep[0][0][0] = 1.0
            self.PYgamma[2].rep[1][0][0] = -1.0
            self.PYgamma[2].rep[2][0][0] = 1.0
            self.PYgamma[2].rep[3][0][0] = -1.0
            self.PYgamma[2].PYnrot = 1
            self.PYgamma[2].PYntrans = 1

            self.PYgamma[3].init(4, 1, "B2", "B2")
            self.PYgamma[3].rep[0][0][0] = 1.0
            self.PYgamma[3].rep[1][0][0] = -1.0
            self.PYgamma[3].rep[2][0][0] = -1.0
            self.PYgamma[3].rep[3][0][0] = 1.0
            self.PYgamma[3].PYnrot = 1
            self.PYgamma[3].PYntrans = 1

        elif self.PYbits in [
            PointGroups['D2']]:
            self.PYgamma[0].init(4, 1, "A", "A")
            self.PYgamma[0].rep[0][0][0] = 1.0
            self.PYgamma[0].rep[1][0][0] = 1.0
            self.PYgamma[0].rep[2][0][0] = 1.0
            self.PYgamma[0].rep[3][0][0] = 1.0
            self.PYgamma[0].PYnrot = 0
            self.PYgamma[0].PYntrans = 0

            self.PYgamma[1].init(4, 1, "B1", "B1")
            self.PYgamma[1].rep[0][0][0] = 1.0
            self.PYgamma[1].rep[1][0][0] = 1.0
            self.PYgamma[1].rep[2][0][0] = -1.0
            self.PYgamma[1].rep[3][0][0] = -1.0
            self.PYgamma[1].PYnrot = 1
            self.PYgamma[1].PYntrans = 1

            self.PYgamma[2].init(4, 1, "B2", "B2")
            self.PYgamma[2].rep[0][0][0] = 1.0
            self.PYgamma[2].rep[1][0][0] = -1.0
            self.PYgamma[2].rep[2][0][0] = 1.0
            self.PYgamma[2].rep[3][0][0] = -1.0
            self.PYgamma[2].PYnrot = 1
            self.PYgamma[2].PYntrans = 1

            self.PYgamma[3].init(4, 1, "B3", "B3")
            self.PYgamma[3].rep[0][0][0] = 1.0
            self.PYgamma[3].rep[1][0][0] = -1.0
            self.PYgamma[3].rep[2][0][0] = -1.0
            self.PYgamma[3].rep[3][0][0] = 1.0
            self.PYgamma[3].PYnrot = 1
            self.PYgamma[3].PYntrans = 1

        elif self.PYbits in [
            PointGroups['D2h']]:

            self.PYgamma[0].init(8, 1, "Ag", "Ag")
            self.PYgamma[0].rep[0][0][0] = 1.0
            self.PYgamma[0].rep[1][0][0] = 1.0
            self.PYgamma[0].rep[2][0][0] = 1.0
            self.PYgamma[0].rep[3][0][0] = 1.0
            self.PYgamma[0].rep[4][0][0] = 1.0
            self.PYgamma[0].rep[5][0][0] = 1.0
            self.PYgamma[0].rep[6][0][0] = 1.0
            self.PYgamma[0].rep[7][0][0] = 1.0
            self.PYgamma[0].PYnrot = 0
            self.PYgamma[0].PYntrans = 0

            self.PYgamma[1].init(8, 1, "B1g", "B1g")
            self.PYgamma[1].rep[0][0][0] = 1.0
            self.PYgamma[1].rep[1][0][0] = 1.0
            self.PYgamma[1].rep[2][0][0] = -1.0
            self.PYgamma[1].rep[3][0][0] = -1.0
            self.PYgamma[1].rep[4][0][0] = 1.0
            self.PYgamma[1].rep[5][0][0] = 1.0
            self.PYgamma[1].rep[6][0][0] = -1.0
            self.PYgamma[1].rep[7][0][0] = -1.0
            self.PYgamma[1].PYnrot = 1
            self.PYgamma[1].PYntrans = 0

            self.PYgamma[2].init(8, 1, "B2g", "B2g")
            self.PYgamma[2].rep[0][0][0] = 1.0
            self.PYgamma[2].rep[1][0][0] = -1.0
            self.PYgamma[2].rep[2][0][0] = 1.0
            self.PYgamma[2].rep[3][0][0] = -1.0
            self.PYgamma[2].rep[4][0][0] = 1.0
            self.PYgamma[2].rep[5][0][0] = -1.0
            self.PYgamma[2].rep[6][0][0] = 1.0
            self.PYgamma[2].rep[7][0][0] = -1.0
            self.PYgamma[2].PYnrot = 1
            self.PYgamma[2].PYntrans = 0

            self.PYgamma[3].init(8, 1, "B3g", "B3g")
            self.PYgamma[3].rep[0][0][0] = 1.0
            self.PYgamma[3].rep[1][0][0] = -1.0
            self.PYgamma[3].rep[2][0][0] = -1.0
            self.PYgamma[3].rep[3][0][0] = 1.0
            self.PYgamma[3].rep[4][0][0] = 1.0
            self.PYgamma[3].rep[5][0][0] = -1.0
            self.PYgamma[3].rep[6][0][0] = -1.0
            self.PYgamma[3].rep[7][0][0] = 1.0
            self.PYgamma[3].PYnrot = 1
            self.PYgamma[3].PYntrans = 0

            self.PYgamma[4].init(8, 1, "Au", "Au")
            self.PYgamma[4].rep[0][0][0] = 1.0
            self.PYgamma[4].rep[1][0][0] = 1.0
            self.PYgamma[4].rep[2][0][0] = 1.0
            self.PYgamma[4].rep[3][0][0] = 1.0
            self.PYgamma[4].rep[4][0][0] = -1.0
            self.PYgamma[4].rep[5][0][0] = -1.0
            self.PYgamma[4].rep[6][0][0] = -1.0
            self.PYgamma[4].rep[7][0][0] = -1.0
            self.PYgamma[4].PYnrot = 0
            self.PYgamma[4].PYntrans = 0

            self.PYgamma[5].init(8, 1, "B1u", "B1u")
            self.PYgamma[5].rep[0][0][0] = 1.0
            self.PYgamma[5].rep[1][0][0] = 1.0
            self.PYgamma[5].rep[2][0][0] = -1.0
            self.PYgamma[5].rep[3][0][0] = -1.0
            self.PYgamma[5].rep[4][0][0] = -1.0
            self.PYgamma[5].rep[5][0][0] = -1.0
            self.PYgamma[5].rep[6][0][0] = 1.0
            self.PYgamma[5].rep[7][0][0] = 1.0
            self.PYgamma[5].PYnrot = 0
            self.PYgamma[5].PYntrans = 1

            self.PYgamma[6].init(8, 1, "B2u", "B2u")
            self.PYgamma[6].rep[0][0][0] = 1.0
            self.PYgamma[6].rep[1][0][0] = -1.0
            self.PYgamma[6].rep[2][0][0] = 1.0
            self.PYgamma[6].rep[3][0][0] = -1.0
            self.PYgamma[6].rep[4][0][0] = -1.0
            self.PYgamma[6].rep[5][0][0] = 1.0
            self.PYgamma[6].rep[6][0][0] = -1.0
            self.PYgamma[6].rep[7][0][0] = 1.0
            self.PYgamma[6].PYnrot = 0
            self.PYgamma[6].PYntrans = 1

            self.PYgamma[7].init(8, 1, "B3u", "B3u")
            self.PYgamma[7].rep[0][0][0] = 1.0
            self.PYgamma[7].rep[1][0][0] = -1.0
            self.PYgamma[7].rep[2][0][0] = -1.0
            self.PYgamma[7].rep[3][0][0] = 1.0
            self.PYgamma[7].rep[4][0][0] = -1.0
            self.PYgamma[7].rep[5][0][0] = 1.0
            self.PYgamma[7].rep[6][0][0] = 1.0
            self.PYgamma[7].rep[7][0][0] = -1.0
            self.PYgamma[7].PYnrot = 0
            self.PYgamma[7].PYntrans = 1

        # Handle symmetry operations
        self.symop[0].E()

        if self.PYbits == PointGroups['C1']:
            pass

        elif self.PYbits == PointGroups['Ci']:
            self.symop[1].i()

        elif self.PYbits == PointGroups['CsX']:  # reflection through the yz plane
            self.symop[1].sigma_yz()

        elif self.PYbits == PointGroups['CsY']:  # reflection through the xz plane
            self.symop[1].sigma_xz()

        elif self.PYbits == PointGroups['CsZ']:  # reflection through the xy plane
            self.symop[1].sigma_xy()

        elif self.PYbits == PointGroups['C2X']:
            self.symop[1].c2_x()

        elif self.PYbits == PointGroups['C2Y']:
            self.symop[1].c2_y()

        elif self.PYbits == PointGroups['C2Z']:
            self.symop[1].rotation(2)

        elif self.PYbits == PointGroups['C2hX']:
            self.symop[1].c2_x()
            self.symop[2].i()
            self.symop[3].sigma_yz()

        elif self.PYbits == PointGroups['C2hY']:
            self.symop[1].c2_y()
            self.symop[2].i()
            self.symop[3].sigma_xz()

        elif self.PYbits == PointGroups['C2hZ']:
            self.symop[1].rotation(2)
            self.symop[2].i()
            self.symop[3].sigma_xy()

        elif self.PYbits == PointGroups['C2vX']:
            self.symop[1].c2_x()
            self.symop[2].sigma_xy()
            self.symop[3].sigma_xz()

        elif self.PYbits == PointGroups['C2vY']:
            self.symop[1].c2_y()
            self.symop[2].sigma_xy()
            self.symop[3].sigma_yz()

        elif self.PYbits == PointGroups['C2vZ']:
            self.symop[1].rotation(2)
            self.symop[2].sigma_xz()
            self.symop[3].sigma_yz()

        elif self.PYbits == PointGroups['D2']:
            self.symop[1].rotation(2)
            self.symop[2].c2_y()
            self.symop[3].c2_x()

        elif self.PYbits == PointGroups['D2h']:
            self.symop[1].rotation(2)
            self.symop[2].c2_y()
            self.symop[3].c2_x()
            self.symop[4].i()
            self.symop[5].sigma_xy()
            self.symop[6].sigma_xz()
            self.symop[7].sigma_yz()

        else:
            return -1

        # now find the inverse of each symop
        for gi in range(self.PYnirrep):
            for gj in range(self.PYnirrep):
                so = self.symop[gi].operate(self.symop[gj])

                # is so a unit matrix?
                if abs(1.0 - so[0][0]) < 1.0e-8 and \
                    abs(1.0 - so[1][1]) < 1.0e-8 and \
                    abs(1.0 - so[2][2]) < 1.0e-8:
                    break

            if gj == self.PYnirrep:
                # ExEnv::err0() << indent
                #      << "make_table: uh oh, can't find inverse of " << gi << endl;
                # abort();
                raise ValidationError("make_table: uh oh, can't find inverse")

            self.inv[gi] = gj

        # Check the bits of the operator make sure they make what
        #   we were given.
        sym_bits = 0
        for i in range(self.PYnirrep):
            sym_bits |= self.symop[i].bit()

        if sym_bits != self.PYbits:
            raise ValidationError("make_table: Symmetry operators did not match the point group given.")

        return 0

    # <<< Methods for Printing >>>

    def __str__(self, out=None):
        """This prints the irrep to the given file, or stdout if none is
        given.

        """
        text = ''
        if not self.PYnirrep:
            return

        text += '  point group %s\n\n' % (self.symb)
        for i in range(self.PYnirrep):
            text += self.PYgamma[i].__str__(out=None)

        text += '\n  symmetry operation matrices:\n\n'
        for i in range(self.PYnirrep):
            text += self.symop[i].__str__(out=None)

        text += '\n  inverse symmetry operation matrices:\n\n'
        for i in range(self.PYnirrep):
            text += self.symop[self.inverse(i)].__str__(out=None)

        if out is None:
            return text
        else:
            with open(out, mode='w') as handle:
                handle.write(text)


class PointGroup(object):
    """The PointGroup class is really a place holder for a CharacterTable.
    It contains a string representation of the Schoenflies symbol of a
    point group, a frame of reference for the symmetry operation
    transformation matrices, and a point of origin. The origin is not
    respected by the symmetry operations, so if you want to use a point
    group with a nonzero origin, first translate all your coordinates to
    the origin and then set the origin to zero.

    """

    def __init__(self, *args):
        """Constructor"""

        # Schoenflies symbol
        self.symb = 'c1'
        # point of origin
        self.PYorigin = [0.0, 0.0, 0.0]
        # bit representation of point group
        self.PYbits = 0

        # Divert to constructor functions
#        if len(args) == 0:
#            self.constructor_zero_ao_basis()
        if len(args) == 1 and \
            isinstance(args[0], basestring):
            self.constructor_schoenflies(*args)
        elif len(args) == 1 and \
            isinstance(args[0], int):
            self.constructor_bits(*args)
        elif len(args) == 2 and \
            isinstance(args[0], basestring) and \
            len(args[1]) == 3:
            self.constructor_schoenflies_origin(*args)
        elif len(args) == 2 and \
            isinstance(args[0], int) and \
            len(args[1]) == 3:
            self.constructor_bits_origin(*args)
        else:
            raise ValidationError('BasisSet::constructor: Inappropriate configuration of constructor arguments')

    # <<< Methods for Construction >>>

    # libmints: These 2 constructors do not work right now.
    def constructor_schoenflies(self, s):
        """This constructor takes a string containing the Schoenflies
        symbol of the point group as its only argument.

        """
        self.PYbits = self.full_name_to_bits(s)
        if self.PYbits is None:
            raise ValidationError('PointGroup: Unknown point group name provided.')
        self.symb = self.bits_to_basic_name(self.PYbits)
        self.PYorigin = [0.0, 0.0, 0.0]

    def constructor_schoenflies_origin(self, s, origin):
        """Like the above, but this constructor also takes a point of
        origin as an argument.

        """
        self.PYbits = self.full_name_to_bits(s)
        if self.PYbits is None:
            raise ValidationError('PointGroup: Unknown point group name provided.')
        self.symb = self.bits_to_basic_name(self.PYbits)
        self.PYorigin = origin

    def constructor_bits(self, bits):
        """Using the bitwise representation constructor the point group
        object.

        """
        self.PYbits = bits
        self.symb = self.bits_to_basic_name(self.PYbits)
        self.PYorigin = [0.0, 0.0, 0.0]

    def constructor_bits_origin(self, bits, origin):
        """Using the bitwise representation constructor the point group
        object.

        """
        self.PYbits = bits
        self.symb = self.bits_to_basic_name(self.PYbits)
        self.PYorigin = origin

    # <<< Simple Methods for Basic PointGroup Information >>>

    def symbol(self):
        """Returns the Schoenflies symbol for this point group."""
        return self.symb

    def set_symbol(self, sym):
        """Sets (or resets) the Schoenflies symbol."""
        self.symb = sym if (len(sym) > 0) else 'c1'

    def origin(self):
        """Returns the origin of the symmetry frame."""
        return self.PYorigin

    def bits(self):
        """Returns the bitwise representation of the point group"""
        return self.PYbits

    def char_table(self):
        """Returns the CharacterTable for this point group."""
        return CharacterTable(self.PYbits)

#    def equiv(self, grp, tol=1.0e-6):
#        """Returns 1 if the point groups *self* and *grp* are equivalent,
#        0 otherwise.
#
#        """
#        return 1 if self.symb == grp.symb else 0

#PointGroup::PointGroup(const PointGroup& pg)
#{
#    *this = pg;
#}
#
#PointGroup::PointGroup(const boost::shared_ptr<PointGroup>& pg)
#{
#    *this = *pg.get();
#}
#

#    """The PointGroup KeyVal constructor looks for three keywords:
#       symmetry, symmetry_frame, and origin. symmetry is a string
#       containing the Schoenflies symbol of the point group. origin is an
#       array of doubles which gives the x, y, and z coordinates of the
#       origin of the symmetry frame. symmetry_frame is a 3 by 3 array of
#       arrays of doubles which specify the principal axes for the
#       transformation matrices as a unitary rotation.
#
#       For example, a simple input which will use the default origin and
#       symmetry_frame ((0,0,0) and the unit matrix, respectively), might
#       look like this:
#
#       <pre>
#       pointgrp<PointGroup>: (
#         symmetry = "c2v"
#       )
#       </pre>
#
#       By default, the principal rotation axis is taken to be the z axis.
#       If you already have a set of coordinates which assume that the
#       rotation axis is the x axis, then you'll have to rotate your frame
#       of reference with symmetry_frame:
#
#       <pre>
#       pointgrp<PointGroup>: (
#         symmetry = "c2v"
#         symmetry_frame = [
#           [ 0 0 1 ]
#           [ 0 1 0 ]
#           [ 1 0 0 ]
#         ]
#       )
#       </pre>
#    """
#    // PointGroup(const Ref<KeyVal>&);
#
#    // PointGroup(StateIn&);
#    PointGroup(const PointGroup&);
#    PointGroup(const boost::shared_ptr<PointGroup>&);
#    ~PointGroup();
#
#    PointGroup& operator=(const PointGroup&);
#PointGroup& PointGroup::operator=(const PointGroup& pg)
#{
#    set_symbol(pg.symb);
#    origin_ = pg.origin_;
#    return *this;
#}
#

    # <<< Methods for Printing >>>

    def __str__(self, out=None):
        text = 'PointGroup: %s\n' % (self.symb)

        if out is None:
            return text
        else:
            with open(out, mode='w') as handle:
                handle.write(text)

    # <<< Methods for Translating Symmetry Encoding >>>

    @staticmethod
    def bits_to_full_name(bits):
        """

        """
        if bits == PointGroups['C1']:
            return "C1"
        elif bits == PointGroups['Ci']:
            return "Ci"
        elif bits == PointGroups['C2X']:
            return "C2(x)"
        elif bits == PointGroups['C2Y']:
            return "C2(y)"
        elif bits == PointGroups['C2Z']:
            return "C2(z)"
        elif bits == PointGroups['CsZ']:
            return "Cs(Z)"
        elif bits == PointGroups['CsY']:
            return "Cs(Y)"
        elif bits == PointGroups['CsX']:
            return "Cs(X)"
        elif bits == PointGroups['D2']:
            return "D2"
        elif bits == PointGroups['C2vX']:
            return "C2v(X)"
        elif bits == PointGroups['C2vY']:
            return "C2v(Y)"
        elif bits == PointGroups['C2vZ']:
            return "C2v(Z)"
        elif bits == PointGroups['C2hX']:
            return "C2h(X)"
        elif bits == PointGroups['C2hY']:
            return "C2h(Y)"
        elif bits == PointGroups['C2hZ']:
            return "C2h(Z)"
        elif bits == PointGroups['D2h']:
            return "D2h"
        else:
            raise ValidationError("Unrecognized point group bits: %d\n" % (bits))

    @staticmethod
    def bits_to_basic_name(bits):
        """From bit representation of point group, returns string of simple
        (non-directional) Schoenflies symbol.

        """
        if bits == PointGroups['C1']:
            return "c1"
        elif bits == PointGroups['Ci']:
            return "ci"
        elif bits in [PointGroups['C2X'], PointGroups['C2Y'], PointGroups['C2Z']]:
            return "c2"
        elif bits in [PointGroups['CsZ'], PointGroups['CsY'], PointGroups['CsX']]:
            return "cs"
        elif bits == PointGroups['D2']:
            return "d2"
        elif bits in [PointGroups['C2vX'], PointGroups['C2vY'], PointGroups['C2vZ']]:
            return "c2v"
        elif bits in [PointGroups['C2hX'], PointGroups['C2hY'], PointGroups['C2hZ']]:
            return "c2h"
        elif bits == PointGroups['D2h']:
            return "d2h"
        else:
            raise ValidationError('Unrecognized point group bits: %d\n' % (bits))

    @staticmethod
    def full_name_to_bits(pg):  # altered signature from (pg, bits):
        """

        """
        pgc = pg.capitalize()

        if pgc == 'C1':
            bits = PointGroups['C1']
        elif pgc == 'Ci':
            bits = PointGroups['Ci']
        elif pgc == 'C2(x)' or pgc == 'C2x' or pgc == 'C2_x':
            bits = PointGroups['C2X']
        elif pgc == 'C2(y)' or pgc == 'C2y' or pgc == 'C2_y':
            bits = PointGroups['C2Y']
        elif pgc == 'C2(z)' or pgc == 'C2z' or pgc == 'C2_z':
            bits = PointGroups['C2Z']
        elif pgc == 'Cs(x)' or pgc == 'Csx' or pgc == 'Cs_x':
            bits = PointGroups['CsX']
        elif pgc == 'Cs(y)' or pgc == 'Csy' or pgc == 'Cs_y':
            bits = PointGroups['CsY']
        elif pgc == 'Cs(z)' or pgc == 'Csz' or pgc == 'Cs_z':
            bits = PointGroups['CsZ']
        elif pgc == 'D2':
            bits = PointGroups['D2']
        elif pgc == 'C2v(x)' or pgc == 'C2vx' or pgc == 'C2v_x':  # changed from C2v(X)
            bits = PointGroups['C2vX']
        elif pgc == 'C2v(y)' or pgc == 'C2vy' or pgc == 'C2v_y':  # changed from C2v(Y)
            bits = PointGroups['C2vY']
        elif pgc == 'C2v(z)' or pgc == 'C2vz' or pgc == 'C2v_z':  # changed from C2v(Z)
            bits = PointGroups['C2vZ']
        elif pgc == 'C2h(x)' or pgc == 'C2hx' or pgc == 'C2h_x':  # changed from C2h(X)
            bits = PointGroups['C2hX']
        elif pgc == 'C2h(y)' or pgc == 'C2hy' or pgc == 'C2h_y':  # changed from C2h(Y)
            bits = PointGroups['C2hY']
        elif pgc == 'C2h(z)' or pgc == 'C2hz' or pgc == 'C2h_z':  # changed from C2h(Z)
            bits = PointGroups['C2hZ']
        elif pgc == 'D2h':
            bits = PointGroups['D2h']

        # Ok, the user gave us Cs, C2v, C2h, C2, but no directionality
        elif pgc == 'Cs':
            bits = PointGroups['CsX']
        elif pgc == 'C2v':
            bits = PointGroups['C2vZ']
        elif pgc == 'C2h':
            bits = PointGroups['C2hZ']
        elif pgc == 'C2':
            bits = PointGroups['C2Z']

        else:
            bits = None

        return bits
