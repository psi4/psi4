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

r"""Module to largely replicate in python the psi4 libmints
CoordValue and CoordEntry classes, which were developed by
Justin M. Turney, with incremental improvements by other
psi4 developers.

"""
from __future__ import absolute_import
from __future__ import print_function
import math
import copy
from .vecutil import *
from .exceptions import *
try:
    from collections import OrderedDict
except ImportError:
    from .oldpymodules import OrderedDict


class CoordValue(object):
    """An abstract class to handle storage of Cartesian coordinate values, which
    may be defined in terms of other variables through this mechanism, greatly
    simplifying Z-matrix specification, for example.

    """

    def __init__(self, fixed=False, computed=False):
        # Fixed coordinate?
        self.PYfixed = fixed
        # Whether the current value is up to date or not
        self.computed = computed

    def set_fixed(self, fixed):
        """Set whether the coordinate value is fixed or not"""
        self.PYfixed = fixed

    def fixed(self):
        """Get whether the coordinate value is fixed or not"""
        return self.PYfixed

    def invalidate(self):
        """Flag the current value as outdated"""
        self.computed = False

    def everything(self):
        print('\nCoordValue\n  Fixed = %s\n  Computed = %s\n\n' % (self.PYfixed, self.computed))


class NumberValue(CoordValue):
    """Specialization of CoordValue that is simply a number to be stored."""

    def __init__(self, value, fixed=False):
        CoordValue.__init__(self, fixed, True)
        # coordinate number value
        self.value = value

    def compute(self):
        """Computes value of coordinate from member data"""
        return self.value

    def rset(self, val):
        """Resets value of coordinate if not fixed"""
        if not self.PYfixed:
            self.value = val

    def type(self):
        """Gets specialization type of CoordValue"""
        return 'NumberType'

    def clone(self):
        """Returns new, independent NumberValue object"""
        return copy.deepcopy(self)

    def variable_to_string(self, precision):
        """Takes a CoordValue object, and returns a string for printing."""
        return "%*.*f" % (precision + 5, precision, self.compute())

    def everything(self):
        print('\nNumberValue\n  Fixed = %s\n  Computed = %s\n  Type = %s\n  Value = %f\n  FValue = %s\n\n' %
            (self.PYfixed, self.computed, self.type(), self.value, self.variable_to_string(4)))


class VariableValue(CoordValue):
    """Specialization of CoordValue, where the current value depends
    on the list of geometry values stored by the molecule.

    """
    def __init__(self, name, geometryVariables, negate=False, fixed=False):
        CoordValue.__init__(self, fixed, True)
        # Name of variable
        self.PYname = name
        # Dictionary from molecule of variable names and values
        self.geometryVariables = geometryVariables
        # Whether the coordinate value is actually the negative of the variable value
        self.negate = negate

    def compute(self):
        """Computes value of coordinate from member data"""
        vstr = self.PYname.upper()
        if vstr not in self.geometryVariables:
            raise IncompleteAtomError('Variable %s used in geometry specification has not been defined' % (vstr))
        if self.negate:
            return self.geometryVariables[vstr] * -1.0
        else:
            return self.geometryVariables[vstr]

    def negated(self):
        """Gets whether the coordinate value is actually the negative of the variable value"""
        return self.negate

    def name(self):
        """Gets the name of the variable"""
        return self.PYname

    def rset(self, val):
        """Resets value of coordinate if not fixed"""
        if not self.PYfixed:
            if self.negate:
                self.geometryVariables[self.PYname] = val * -1.0
            else:
                self.geometryVariables[self.PYname] = val

    def type(self):
        """Gets specialization type of CoordValue"""
        return 'VariableType'

    def clone(self):
        """Returns new, independent VariableValue object"""
        return copy.deepcopy(self)

    def variable_to_string(self, precision):
        """Takes a CoordValue object, and returns a string for printing."""
        if self.negate:
            return '-' + self.PYname
        else:
            return self.PYname

    def everything(self):
        print('\nVariableValue\n  Fixed = %s\n  Computed = %s\n  Type = %s\n  Value = %f\n  FValue = %s\n  Name = %s\n  Negated = %s\n  Map = %s\n\n' %
            (self.PYfixed, self.computed, self.type(), self.compute(), self.variable_to_string(4), self.name(), self.negated(), self.geometryVariables))


class CoordEntry(object):
    """Class to store all the attributes associated with an atom, not the
    larger Molecule. Specialized into CartesianEntry and ZMatrixEntry.

    """

    def __init__(self, entry_number, Z, charge, mass, symbol, label="", basis=None, shells=None):
        """Constructor"""
        # Order in full atomic list
        self.PYentry_number = entry_number
        # Whether the coordinates have been computed
        self.computed = False
        # Actual cartesian coordinates of the atom
        self.coordinates = [None, None, None]
        # Atomic number of the atom
        self.PYZ = Z
        # Charge of the atom (SAD-related)
        self.PYcharge = charge
        # Mass of the atom
        self.PYmass = mass
        # Label of the atom minus any extra info (H1 => H)
        self.PYsymbol = symbol
        # Original label from the molecule from the input file (H1)
        self.PYlabel = label
        # Is this a ghost atom?
        self.ghosted = False
        # Different types of basis sets that can be assigned to this atom.
        self.PYbasissets = basis if basis is not None else OrderedDict()
        # Hash of one-atom BasisSet attached to this atom
        self.PYshells = shells if shells is not None else OrderedDict()

    @staticmethod
    def r(a1, a2):
        """Computes the distance between two sets of coordinates"""
        if len(a1) != 3 or len(a2) != 3:
            raise ValidationError('ERROR: r() only defined for Vector3\n')
        return distance(a1, a2)

    @staticmethod
    def a(a1, a2, a3):
        """Computes the angle (in rad.) between three sets of coordinates."""
        if len(a1) != 3 or len(a2) != 3 or len(a3) != 3:
            raise ValidationError('ERROR: a() only defined for Vector3\n')
        eBA = sub(a2, a1)
        eBC = sub(a2, a3)
        eBA = normalize(eBA)
        eBC = normalize(eBC)
        costheta = dot(eBA, eBC)

        if costheta > 1.0 - 1.0E-14:
            costheta = 1.0
        if costheta < 1.0E-14 - 1.0:
            costheta = -1.0
        return math.acos(costheta)

    @staticmethod
    def d(a1, a2, a3, a4):
        """Computes the dihedral (in rad.) between four sets of coordinates."""
        if len(a1) != 3 or len(a2) != 3 or len(a3) != 3 or len(a4) != 3:
            raise ValidationError('ERROR: d() only defined for Vector3\n')
        eBA = sub(a2, a1)
        eDC = sub(a4, a3)
        eCB = sub(a3, a2)
        CBNorm = norm(eCB)
        DCxCB = cross(eDC, eCB)
        CBxBA = cross(eCB, eBA)
        return -1.0 * math.atan2(CBNorm * dot(eDC, CBxBA), dot(DCxCB, CBxBA))

    def is_computed(self):
        """Whether the current atom's coordinates are up-to-date."""
        return self.computed

    def is_equivalent_to(self, other):
        """Whether this atom has the same mass and ghost status as atom *other*.
        Also compares basis set assignment down to nbf(), has_puream() level
        with code borrowed from Robert M. Parrish's SAD guess in Psi4.

        """
        if other.PYZ != self.PYZ:
            return False
        if other.PYmass != self.PYmass:
            return False
        if other.ghosted != self.ghosted:
            return False
        if other.PYshells is not None and self.PYshells is not None:
            for bas in self.PYshells:  # do we instead care only about orbital basis?
                if bas in other.PYshells:
                    if other.PYshells[bas] != self.PYshells[bas]:
                        return False
                    #if other.PYshells[bas].nbf() != self.PYshells[bas].nbf():
                    #    return False
                    #if other.PYshells[bas].nshell() != self.PYshells[bas].nshell():
                    #    return False
                    #if other.PYshells[bas].nprimitive() != self.PYshells[bas].nprimitive():
                    #    return False
                    #if other.PYshells[bas].max_am() != self.PYshells[bas].max_am():
                    #    return False
                    #if other.PYshells[bas].max_nprimitive() != self.PYshells[bas].max_nprimitive():
                    #    return False
                    #if other.PYshells[bas].has_puream() != self.PYshells[bas].has_puream():
                    #    return False
                else:
                    raise ValidationError("""Basis set %s set for one and not other. This shouldn't happen. Investigate.""" % (bas))
        return True

    def is_ghosted(self):
        """Whether the current atom is ghosted or not."""
        return self.ghosted

    def set_ghosted(self, gh):
        """Flag the atom as either ghost or real."""
        self.ghosted = gh

    def Z(self):
        """The nuclear charge of the current atom (0 if ghosted)."""
        if self.ghosted:
            return 0.0
        else:
            return self.PYZ

    def charge(self):
        """The "atomic charge" of the current atom (for SAD purposes)."""
        return self.PYcharge

    def mass(self):
        """The atomic mass of the current atom."""
        return self.PYmass

    def symbol(self):
        """The atomic symbol."""
        return self.PYsymbol

    def label(self):
        """The atom label."""
        return self.PYlabel

    def entry_number(self):
        """The order in which this appears in the full atom list."""
        return self.PYentry_number

    def set_basisset(self, name, role='BASIS'):
        """Set the basis for this atom
        * @param type Keyword from input file, basis, ri_basis, etc.
        * @param name Value from input file

        """
        self.PYbasissets[role] = name

    def basisset(self, role='BASIS'):
        """Returns the basis name for the provided type.
        * @param type Keyword from input file.
        * @returns the value from input.

        """
        try:
            return self.PYbasissets[role]
        except ValueError:
            raise ValidationError('CoordEntry::basisset: Basisset not set for %s and type of %s' % \
                (self.PYlabel, role))

    def basissets(self):
        """Returns basisset to atom map"""
        return self.PYbasissets

    def set_shell(self, bshash, key='BASIS'):
        """Set the hash for this atom
        * @param key Keyword from input file, basis, ri_basis, etc.
        * @param bshash hash string of one-atom BasisSet

        """
        self.PYshells[key] = bshash

    def shell(self, key='BASIS'):
        """Returns the hash for the provided type.
        * @param type Keyword from input file.
        * @returns the hash string for basis.

        """
        try:
            return self.PYshells[key]
        except (ValueError, KeyError):
            raise ValidationError('CoordEntry::shells: Shells not set for %s and type of %s' % \
                (self.PYlabel, key))

    def shells(self):
        """Returns shells sets to atom map"""
        return self.PYshells

    def everything(self):
        print('\nCoordEntry\n  Entry Number = %d\n  Computed = %s\n  Z = %d\n  Charge = %f\n  Mass = %f\n  Symbol = %s\n  Label = %s\n  Ghosted = %s\n  Coordinates = %s\n  Basissets = %s\n\n  Shells = %s\n\n' %
            (self.entry_number(), self.is_computed(), self.Z(), self.charge(),
            self.mass(), self.symbol(), self.label(), self.is_ghosted(),
            self.coordinates, self.PYbasissets, self.PYshells))


class CartesianEntry(CoordEntry):
    """Class to hold all information about an atom, including its
    coordinate specification as three Cartesians.

    """

    def __init__(self, entry_number, Z, charge, mass, symbol, label, x, y, z, basis=None, shells=None):
        CoordEntry.__init__(self, entry_number, Z, charge, mass, symbol, label, basis, shells)
        self.x = x
        self.y = y
        self.z = z

    def compute(self):
        """Computes the values of the coordinates (in whichever units
        were inputted), returning them in a Vector

        """
        if self.computed:
            return self.coordinates
        self.coordinates[0] = self.x.compute()
        self.coordinates[1] = self.y.compute()
        self.coordinates[2] = self.z.compute()
        self.computed = True
        return self.coordinates

    def set_coordinates(self, x, y, z):
        """Given the current set of coordinates, updates the values of this
        atom's coordinates and any variables that may depend on it.

        """
        self.coordinates[0] = x
        self.coordinates[1] = y
        self.coordinates[2] = z

        self.x.rset(x)
        self.y.rset(y)
        self.z.rset(z)

        self.computed = True

    def type(self):
        """The type of CoordEntry specialization."""
        return 'CartesianCoord'

    def print_in_input_format(self):
        """Prints the updated geometry, in the format provided by the user."""
        xstr = self.x.variable_to_string(12)
        ystr = self.y.variable_to_string(12)
        zstr = self.z.variable_to_string(12)
        return "  %17s  %17s  %17s\n" % (xstr, ystr, zstr)
        # should go to outfile

    def print_in_input_format_cfour(self):
        """Prints the updated geometry, in the format provided by the user.
        This, for Cfour, not different from regular version.

        """
        xstr = self.x.variable_to_string(12)
        ystr = self.y.variable_to_string(12)
        zstr = self.z.variable_to_string(12)
        return " %17s %17s %17s\n" % (xstr, ystr, zstr)
        # should go to outfile

    def invalidate(self):
        """Flags the current coordinates as being outdated."""
        self.computed = False
        self.x.invalidate()
        self.y.invalidate()
        self.z.invalidate()

    def clone(self):
        """Returns new, independent CartesianEntry object"""
        return copy.deepcopy(self)

    def everything(self):
        CoordEntry.everything(self)
        print('\nCartesianEntry\n  Type = %s\n  x = %s\n  y = %s\n  z = %s\n\n' % (self.type(), self.x.variable_to_string(8), self.y.variable_to_string(8), self.z.variable_to_string(8)))


class ZMatrixEntry(CoordEntry):
    """Class to hold all information about an atom, including its
    coordinate specification as any position of ZMatrix.

    """

    def __init__(self, entry_number, Z, charge, mass, symbol, label, \
        rto=None, rval=0, ato=None, aval=0, dto=None, dval=0, basis=None, shells=None):
        """Constructor"""  # note that pos'n of basis arg changed from libmints
        CoordEntry.__init__(self, entry_number, Z, charge, mass, symbol, label, basis, shells)
        self.rto = rto
        self.rval = rval
        self.ato = ato
        self.aval = aval
        self.dto = dto
        self.dval = dval

    def invalidate(self):
        """Flags the current coordinates as being outdated"""
        self.computed = False
        if self.rval != 0:
            self.rval.invalidate()
        if self.aval != 0:
            self.aval.invalidate()
        if self.dval != 0:
            self.dval.invalidate()

    def print_in_input_format(self):
        """Prints the updated geometry, in the format provided by the user"""
        text = ""
        if self.rto == None and self.ato == None and self.dto == None:
            # The first atom
            text += "\n"
        elif self.ato == None and self.dto == None:
            # The second atom
            now_rto = self.rto.entry_number() + 1
            now_rval = self.rval.variable_to_string(10)
            text += "  %5d %11s\n" % (now_rto, now_rval)
        elif self.dto == None:
            # The third atom
            now_rto = self.rto.entry_number() + 1
            now_rval = self.rval.variable_to_string(10)
            now_ato = self.ato.entry_number() + 1
            now_aval = self.aval.variable_to_string(10)
            text += "  %5d %11s  %5d %11s\n" % (now_rto, now_rval, now_ato, now_aval)
        else:
            # Remaining atoms
            now_rto = self.rto.entry_number() + 1
            now_rval = self.rval.variable_to_string(10)
            now_ato = self.ato.entry_number() + 1
            now_aval = self.aval.variable_to_string(10)
            now_dto = self.dto.entry_number() + 1
            now_dval = self.dval.variable_to_string(10)
            text += "  %5d %11s  %5d %11s  %5d %11s\n" % \
                (now_rto, now_rval, now_ato, now_aval, now_dto, now_dval)
        return text
#        outfile

    def print_in_input_format_cfour(self):
        """Prints the updated geometry, in the format provided by the user"""
        text = ""
        if self.rto == None and self.ato == None and self.dto == None:
            # The first atom
            text += "\n"
        elif self.ato == None and self.dto == None:
            # The second atom
            now_rto = self.rto.entry_number() + 1
            now_rval = self.rval.variable_to_string(10)
            text += " %d %s\n" % (now_rto, now_rval)
        elif self.dto == None:
            # The third atom
            now_rto = self.rto.entry_number() + 1
            now_rval = self.rval.variable_to_string(10)
            now_ato = self.ato.entry_number() + 1
            now_aval = self.aval.variable_to_string(10)
            text += " %d %s %d %s\n" % (now_rto, now_rval, now_ato, now_aval)
        else:
            # Remaining atoms
            now_rto = self.rto.entry_number() + 1
            now_rval = self.rval.variable_to_string(10)
            now_ato = self.ato.entry_number() + 1
            now_aval = self.aval.variable_to_string(10)
            now_dto = self.dto.entry_number() + 1
            now_dval = self.dval.variable_to_string(10)
            text += " %d %s %d %s %d %s\n" % \
                (now_rto, now_rval, now_ato, now_aval, now_dto, now_dval)
        return text
#        outfile

    def set_coordinates(self, x, y, z):
        """Given the current set of coordinates, updates the values of this
        atom's coordinates, and any variables that may depend on it.

        """
        self.coordinates[0] = 0.0 if math.fabs(x) < 1.0E-14 else x
        self.coordinates[1] = 0.0 if math.fabs(y) < 1.0E-14 else y
        self.coordinates[2] = 0.0 if math.fabs(z) < 1.0E-14 else z

        if self.rto != None:
            if not self.rto.is_computed():
                raise ValidationError("Coordinates have been set in the wrong order")
            self.rval.rset(self.r(self.coordinates, self.rto.compute()))

        if self.ato != None:
            if not self.ato.is_computed():
                raise ValidationError("Coordinates have been set in the wrong order")
            aval = self.a(self.coordinates, self.rto.compute(), self.ato.compute())
            # Noise creeps in for linear molecules. Force linearity, if it is close enough.
            val = aval * 180.0 / math.pi
            self.aval.rset(val)

        if self.dto != None:
            if not self.dto.is_computed():
                raise ValidationError("Coordinates have been set in the wrong order")
            val = self.d(self.coordinates, self.rto.compute(), self.ato.compute(), self.dto.compute())
            # Check for NaN, and don't update if we find one
            # what is this? proper py traslation?
            if val == val:
                self.dval.rset(val * 180.0 / math.pi)

        self.computed = True

    def type(self):
        """The type of CoordEntry specialization."""
        return 'ZMatrixCoord'

    def clone(self):
        """Returns new, independent ZMatrixEntry object."""
        return copy.deepcopy(self)

    def compute(self):
        """Compute the Cartesian coordinates in Bohr of current atom's entry."""

        if self.computed:
            return self.coordinates

        # place first atom at the origin
        if self.rto == None and self.ato == None and self.dto == None:
            self.coordinates[0] = 0.0
            self.coordinates[1] = 0.0
            self.coordinates[2] = 0.0

        # place second atom directly above the first
        elif self.ato == None and self.dto == None:
            self.coordinates[0] = 0.0
            self.coordinates[1] = 0.0
            self.coordinates[2] = self.rval.compute()

        # place third atom pointing upwards
        #    this       rTo   rVal  aTo  aVal
        #      A         B           C
        elif self.dto == None:
            r = self.rval.compute()
            a = self.aval.compute() * math.pi / 180.0
            cosABC = math.cos(a)
            sinABC = math.sin(a)
            B = self.rto.compute()
            C = self.ato.compute()

            eCB = sub(B, C)
            eCB = normalize(eCB)
            eX = [0.0, 0.0, 0.0]
            eY = [0.0, 0.0, 0.0]
            if (math.fabs(1.0 - math.fabs(eCB[0])) < 1.0E-5):
                # CB is collinear with X, start by finding Y
                eY[1] = 1.0
                eX = perp_unit(eY, eCB)
                eY = perp_unit(eX, eCB)
            else:
                # CB is not collinear with X, we can safely find X first
                eX[0] = 1.0
                eY = perp_unit(eX, eCB)
                eX = perp_unit(eY, eCB)
            for xyz in range(3):
                self.coordinates[xyz] = B[xyz] + r * (eY[xyz] * sinABC - eCB[xyz] * cosABC)
                if math.fabs(self.coordinates[xyz]) < 1.E-14:
                    self.coordinates[xyz] = 0.0

        # The fourth, or subsequent, atom
        #
        # The atom specification is
        #      this       rTo   rVal  aTo  aVal   dTo   dVal
        #        A         B           C           D
        # which allows us to define the vector from C->B (eCB) as the +z axis, and eDC
        # lies in the xz plane.  Then eX, eY and eZ (=eBC) are the x, y, and z axes, respecively.
        else:
            r = self.rval.compute()
            a = self.aval.compute() * math.pi / 180.0
            d = self.dval.compute() * math.pi / 180.0
            B = self.rto.compute()
            C = self.ato.compute()
            D = self.dto.compute()

            eDC = sub(C, D)
            eCB = sub(B, C)
            eDC = normalize(eDC)
            eCB = normalize(eCB)
            cosABC = math.cos(a)
            sinABC = math.sin(a)
            cosABCD = math.cos(d)
            sinABCD = math.sin(d)
            eY = perp_unit(eDC, eCB)
            eX = perp_unit(eY, eCB)
            for xyz in range(3):
                self.coordinates[xyz] = B[xyz] + r * (eX[xyz] * sinABC * cosABCD +
                         eY[xyz] * sinABC * sinABCD - eCB[xyz] * cosABC)
                if math.fabs(self.coordinates[xyz]) < 1.E-14:
                    self.coordinates[xyz] = 0.0

        self.computed = True
        return self.coordinates

    def everything(self):
        CoordEntry.everything(self)
        print('\nZMatrixEntry\n  Type = %s\n\n' % (self.type()))
        print(self.print_in_input_format())
