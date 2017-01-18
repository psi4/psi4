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
import os
import re
import copy
import math
try:
    from collections import OrderedDict
except ImportError:
    from .oldpymodules import OrderedDict
from .periodictable import *
from .physconst import *
from .vecutil import *
from .exceptions import *
#from libmintscoordentry import *
from .libmintscoordentry import NumberValue, VariableValue, CartesianEntry, ZMatrixEntry
from .libmintspointgrp import SymmOps, similar, SymmetryOperation, PointGroup

#from libmintspointgrp import PointGroups
#print PointGroups


LINEAR_A_TOL = 1.0E-2  # When sin(a) is below this, we consider the angle to be linear
DEFAULT_SYM_TOL = 1.0E-8
FULL_PG_TOL = 1.0e-8
ZERO = 1.0E-14
NOISY_ZERO = 1.0E-8


class LibmintsMolecule(dict):
    """Class to store the elements, coordinates, fragmentation pattern,
    charge, multiplicity of a molecule. Largely replicates psi4's libmints
    Molecule class, developed by Justin M. Turney and Andy M. Simmonett
    with incremental improvements by other psi4 developers. Major
    differences from the C++ class are: no basisset handling, no symmetry,
    no pubchem, no efp, no discarding dummies. This class translated so
    that databases can function independently of psi4.

    >>> H2OH2O = qcdb.Molecule(\"\"\"
        0 1
        O1  -1.551007  -0.114520   0.000000
        H1  -1.934259   0.762503   0.000000
        H2  -0.599677   0.040712   0.000000
        --
        0 1
        X    0.000000   0.000000   0.000000
        O2   1.350625   0.111469   0.000000
        H3   1.680398  -0.373741  -0.758561
        H4   1.680398  -0.373741   0.758561
        no_com
        no_reorient
        units angstrom
        \"\"\")

    >>> H2O = qcdb.Molecule.init_with_xyz('h2o.xyz')
    """
    FullPointGroupList = ["ATOM", "C_inf_v", "D_inf_h", "C1", "Cs", "Ci", \
        "Cn", "Cnv", "Cnh", "Sn", "Dn", "Dnd", "Dnh", "Td", "Oh", "Ih"]

    def __init__(self, psi4molstr=None):
        """Initialize Molecule object from string in psi4 format"""

        # <<< Basic Molecule Information >>>

        # Molecule (or fragment) name
        self.PYname = 'default'
        # The molecular charge
        self.PYmolecular_charge = 0
        # Whether the charge was given by the user
        self.PYcharge_specified = False
        # The multiplicity (defined as 2Ms + 1)
        self.PYmultiplicity = 1
        # Whether the multiplicity was specified by the user
        self.PYmultiplicity_specified = False
        # The units used to define the geometry
        self.PYunits = 'Angstrom'
        # The conversion factor to take input units to Bohr
        self.input_units_to_au = 1.0 / psi_bohr2angstroms
        # Whether this molecule has at least one zmatrix entry
        self.zmat = False  # TODO None?

        # <<< Coordinates >>>

        # Atom info vector (no knowledge of dummy atoms)
        self.atoms = []
        # Atom info vector (includes dummy atoms)
        self.full_atoms = []
        # A list of all variables known, whether they have been set or not.
        self["all_variables"] = []
        # A listing of the variables used to define the geometries
        self["geometry_variables"] = {}

        # <<< Fragmentation >>>

        # The list of atom ranges defining each fragment from parent molecule
        self.fragments = []
        # A list describing how to handle each fragment
        self.fragment_types = []
        # The charge of each fragment
        self.fragment_charges = []
        # The multiplicity of each fragment
        self.fragment_multiplicities = []

        # <<< Frame >>>

        # Move to center of mass or not?
        self.PYmove_to_com = True
        # Reorient or not?  UNUSED
        self.PYfix_orientation = False
        # Reinterpret the coord entries or not (Default is true, except for findif)
        self.PYreinterpret_coordentries = True
        # Nilpotence boolean (flagged upon first determination of symmetry frame,
        #    reset each time a substantiative change is made)
        self.lock_frame = False

        # <<< Symmetry >>>

        # Point group to use with this molecule
        self.pg = None
        # Full point group
        self.full_pg = 'C1'
        # n of the highest rotational axis Cn
        self.PYfull_pg_n = 1
        # Symmetry string from geometry specification
        self.PYsymmetry_from_input = None
        # Number of unique atoms
        self.PYnunique = 0
        # Number of equivalent atoms per unique atom (length nunique)
        self.nequiv = 0
        # Equivalent atom mapping array (length 1st dim nunique)
        self.equiv = 0
        # Atom to unique atom mapping array (length natom)
        self.PYatom_to_unique = 0

        if psi4molstr:
            self.create_molecule_from_string(psi4molstr)

    # <<< Simple Methods for Basic Molecule Information >>>

    def name(self):
        """Get molecule name

        >>> print H2OH2O.name()
        water_dimer

        """
        return self.PYname

    def set_name(self, name):
        """Set molecule name

        >>> H2OH2O.set_name('water_dimer')

        """
        self.PYname = name

    def natom(self):
        """Number of atoms

        >>> print H2OH2O.natom()
        6

        """
        return len(self.atoms)

    def nallatom(self):
        """Number of all atoms (includes dummies)

        >>> print H2OH2O.nallatom()
        7

        """
        return len(self.full_atoms)

    def molecular_charge(self):
        """Gets the molecular charge

        >>> print H2OH2O.molecular_charge()
        -2

        """
        return self.PYmolecular_charge

    def set_molecular_charge(self, charge):
        """Sets the molecular charge

        >>> H2OH2O.set_molecular_charge(-2)

        """
        self.PYcharge_specified = True
        self.PYmolecular_charge = charge

    def charge_specified(self):
        """Whether the charge was given by the user

        >>> print H2OH2O.charge_specified()
        True

        """
        return self.PYcharge_specified

    def multiplicity(self):
        """Get the multiplicity (defined as 2Ms + 1)

        >>> print H2OH2O.multiplicity()

        """
        return self.PYmultiplicity

    def set_multiplicity(self, mult):
        """Sets the multiplicity (defined as 2Ms + 1)

        >>> H2OH2O.set_multiplicity(3)

        """
        self.PYmultiplicity_specified = True
        self.PYmultiplicity = mult

    def multiplicity_specified(self):
        """Whether the multiplicity was given by the user

        >>> print H2OH2O.multiplicity_specified()
        True

        """
        return self.PYmultiplicity_specified

    def units(self):
        """Gets the geometry units

        >>> print H2OH2O.units()
        Angstrom

        """
        return self.PYunits

    def set_units(self, units):
        """Sets the geometry units

        >>> H2OH2O.set_units('Angstom')

        """
        if units == 'Angstrom' or units == 'Bohr':
            self.PYunits = units
        else:
            raise ValidationError("""Molecule::set_units: argument must be 'Angstrom' or 'Bohr'.""")

    def has_zmatrix(self):
        """Gets the presence of any zmatrix entry

        >>> print H2OH2O.has_zmatrix()
        False

        """
        return self.zmat

    def set_has_zmatrix(self, tf):
        """Sets the presence of any zmatrix entry

        >>> H2OH2O.set_has_zmatrix(True)

        """
        self.zmat = tf

    # <<< Simple Methods for Coordinates >>>

    def Z(self, atom):
        """Nuclear charge of atom (0-indexed)

        >>> print H2OH2O.Z(4)
        1

        """
        return self.atoms[atom].Z()

    def x(self, atom):
        """x position of atom (0-indexed) in Bohr

        >>> print H2OH2O.x(4)
        3.17549201425

        """
        return self.input_units_to_au * self.atoms[atom].compute()[0]

    def y(self, atom):
        """y position of atom (0-indexed) in Bohr

        >>> print H2OH2O.y(4)
        -0.706268134631

        """
        return self.input_units_to_au * self.atoms[atom].compute()[1]

    def z(self, atom):
        """z position of atom (0-indexed) in Bohr

        >>> print H2OH2O.z(4)
        -1.43347254509

        """
        return self.input_units_to_au * self.atoms[atom].compute()[2]

    def xyz(self, atom, posn=None):
        """Returns a Vector3 with x, y, z position of atom (0-indexed)
        in Bohr or coordinate at *posn*

        >>> print H2OH2O.xyz(4)
        [3.175492014248769, -0.7062681346308132, -1.4334725450878665]

        """
        temp = scale(self.atoms[atom].compute(), self.input_units_to_au)
        if posn is None:
            return temp
        else:
            return temp[posn]

    def mass(self, atom):
        """Returns mass of atom (0-indexed)

        >>> print H2OH2O.mass(4)
        1.00782503207

        """
        if self.atoms[atom].mass() != 0.0:
            return self.atoms[atom].mass()

        if math.fabs(self.atoms[atom].Z() - int(self.atoms[atom].Z())) > 0.0:
            print("""WARNING: Obtaining masses from atom with fractional charge...may be incorrect!!!\n""")
            # TODO outfile
        return z2mass[int(self.atoms[atom].Z())]

    def symbol(self, atom):
        """Returns the cleaned up label of the atom (C2 => C, H4 = H) (0-indexed)

        >>> print H2OH2O.symbol(4)
        H

        """
        return self.atoms[atom].symbol()

    def label(self, atom):
        """Returns the original label of the atom (0-indexed) as given in the input file (C2, H4). (0-indexed)

        >>> print H2OH2O.label(4)
        H3

        """
        return self.atoms[atom].label()

    def charge(self, atom):
        """Returns charge of atom (0-indexed).
        Related to SAD guess in libmints version.

        >>> print H2OH2O.charge(4)
        1.0

        """
        return self.atoms[atom].charge()

    def fZ(self, atom):
        """Nuclear charge of atom (includes dummies)

        >>> print H2OH2O.fZ(4)
        8

        """
        return self.full_atoms[atom].Z()

    def fx(self, atom):
        """x position of atom (0-indexed, includes dummies) in Bohr

        >>> print H2OH2O.fx(4)
        2.55231135823

        """
        return self.input_units_to_au * self.full_atoms[atom].compute()[0]

    def fy(self, atom):
        """y position of atom (0-indexed, includes dummies) in Bohr

        >>> print H2OH2O.fy(4)
        0.210645882307

        """
        return self.input_units_to_au * self.full_atoms[atom].compute()[1]

    def fz(self, atom):
        """z position of atom (0-indexed, includes dummies) in Bohr

        >>> print H2OH2O.fz(4)
        0.0

        """
        return self.input_units_to_au * self.full_atoms[atom].compute()[2]

    def fxyz(self, atom):
        """Returns a Vector3 with x, y, z position of atom
        (0-indexed) in Bohr (includes dummies)

        >>> print H2OH2O.fxyz(4)
        [2.5523113582286716, 0.21064588230662976, 0.0]

        """
        return scale(self.full_atoms[atom].compute(), self.input_units_to_au)

    def fmass(self, atom):
        """Returns mass of atom (0-indexed, includes dummies)

        >>> print H2OH2O.fmass(4)
        15.9949146196

        """
        return self.full_atoms[atom].mass()

    def fsymbol(self, atom):
        """Returns the cleaned up label of the atom (C2 => C, H4 = H) (includes dummies) (0-indexed)

        >>> print H2OH2O.fsymbol(4)
        O

        """
        return self.full_atoms[atom].symbol()

    def flabel(self, atom):
        """Returns the original label of the atom (0-indexed) as given in
        the input file (C2, H4) (includes dummies)

        >>> print H2OH2O.flabel(4)
        O2

        """
        return self.full_atoms[atom].label()

    def fcharge(self, atom):
        """Returns charge of atom (0-indexed, includes dummies).
        Related to SAD guess in libmints version.

        >>> print H2OH2O.fcharge(4)
        8.0

        """
        return self.full_atoms[atom].charge()

    # <<< Simple Methods for Fragmentation >>>

    def nfragments(self):
        """The number of fragments in the molecule.

        >>> print H2OH2O.nfragments()
        2

        """
        return len(self.fragments)

    def nactive_fragments(self):
        """The number of active fragments in the molecule.

        >>> print H2OH2O.nactive_fragments()
        2

        """
        n = 0
        for fr in range(self.nfragments()):
            if self.fragment_types[fr] == 'Real':
                n += 1
        return n

    def activate_all_fragments(self):
        """Sets all fragments in the molecule to be active."""
        self.lock_frame = False
        for fr in range(self.nfragments()):
            self.fragment_types[fr] = 'Real'

    def set_active_fragment(self, fr):
        """Tags fragment index *fr* as composed of real atoms."""
        self.lock_frame = False
        self.fragment_types[fr - 1] = 'Real'

    def set_active_fragments(self, reals):
        """Tags the fragments in array *reals* as composed of real atoms."""
        self.lock_frame = False
        for fr in reals:
            self.fragment_types[fr - 1] = 'Real'

    def set_ghost_fragment(self, fr):
        """Tags fragment index *fr* as composed of ghost atoms."""
        self.lock_frame = False
        self.fragment_types[fr - 1] = 'Ghost'

    def set_ghost_fragments(self, ghosts):
        """Tags the fragments in array *ghosts* as composed of ghost atoms."""
        self.lock_frame = False
        for fr in ghosts:
            self.fragment_types[fr - 1] = 'Ghost'

    def deactivate_all_fragments(self):
        """Sets all fragments in the molecule to be inactive."""
        self.lock_frame = False
        for fr in range(self.nfragments()):
            self.fragment_types[fr] = 'Absent'

    def extract_subsets(self, reals, ghosts=[]):
        """Wrapper for :py:func:`~qcdb.molecule.extract_fragments`.
        See note there. This function can be used as long as not
        in psi4 input file. Use extract_fragments directly, then.

        >>> H2OH2O.extract_subsets(2)  # monomer B, unCP-corrected
        >>> H2OH2O.extract_subsets(2,1)  # monomer B, CP-corrected
        >>> obj.extract_subsets(1,[2,3])  # monomer A, CP-corrected if obj is tri-molecular complex

        """
        return self.extract_fragments(reals, ghosts=ghosts)

    def extract_fragments(self, reals, ghosts=[]):
        """Makes a copy of the molecule, returning a new molecule with
        only certain fragment atoms present as either ghost or real atoms
        *reals*: The list or int of fragments (1-indexed) that should be present in the molecule as real atoms.
        *ghosts*: The list or int of fragments (1-indexed) that should be present in the molecule as ghosts.
        (method name in libmints is extract_subsets. This is different
        in qcdb because the psi4 input parser tries to process lines with
        that term, giving rise to Boost:Python type conlicts.) See usage
        at :py:func:`~qcdb.molecule.extract_fragments`.

        """
        lreals = []
        try:
            for idx in reals:
                lreals.append(idx - 1)
        except TypeError:
            lreals = [reals - 1]
        lghosts = []
        try:
            for idx in ghosts:
                lghosts.append(idx - 1)
        except TypeError:
            lghosts = [ghosts - 1]
        if len(lreals) + len(lghosts) > self.nfragments():
            raise ValidationError('Molecule::extract_fragments: sum of real- and ghost-atom subsets is greater than the number of subsets')

        subset = self.clone()
        subset.deactivate_all_fragments()
        for fr in lreals:
            subset.set_active_fragment(fr + 1)  # the active fragment code subtracts 1
        for fr in lghosts:
            subset.set_ghost_fragment(fr + 1)  # the ghost fragment code subtracts 1

        subset.update_geometry()
        return subset

    # <<< Methods for Construction >>>

    def create_molecule_from_string(self, text):
        """Given a string *text* of psi4-style geometry specification
        (including newlines to separate lines), builds a new molecule.
        Called from constructor.

        """
        comment = re.compile(r'^\s*#')
        blank = re.compile(r'^\s*$')
        bohr = re.compile(r'^\s*units?[\s=]+(bohr|au|a.u.)\s*$', re.IGNORECASE)
        ang = re.compile(r'^\s*units?[\s=]+(ang|angstrom)\s*$', re.IGNORECASE)
        orient = re.compile(r'^\s*(no_reorient|noreorient)\s*$', re.IGNORECASE)
        com = re.compile(r'^\s*(no_com|nocom)\s*$', re.IGNORECASE)
        symmetry = re.compile(r'^\s*symmetry[\s=]+(\w+)\s*$', re.IGNORECASE)
        ATOM = '((([A-Z]{1,3})_\w+)|(([A-Z]{1,3})\d*))'  # match 'C', 'al', 'p88', 'p_pass' not 'Ofail', 'h99_text'  # good, but unused
        atom = re.compile(r'^(?:(?P<gh1>@)|(?P<gh2>Gh\())?(?P<label>(?P<symbol>[A-Z]{1,3})(?:(_\w+)|(\d+))?)(?(gh2)\))(?:@(?P<mass>\d+\.\d+))?$', re.IGNORECASE)
        cgmp = re.compile(r'^\s*(-?\d+)\s+(\d+)\s*$')
        frag = re.compile(r'^\s*--\s*$')
        variable = re.compile(r'^\s*(\w+)\s*=\s*(-?\d+\.\d+|-?\d+\.|-?\.\d+|-?\d+|tda)\s*$', re.IGNORECASE)
        ghost = re.compile(r'@(.*)|Gh\((.*)\)', re.IGNORECASE)

        lines = re.split('\n', text)
        glines = []
        ifrag = 0

        for line in lines:

            # handle comments
            if comment.match(line) or blank.match(line):
                pass

            # handle units
            elif ang.match(line):
                self.set_units('Angstrom')
                self.input_units_to_au = 1.0 / psi_bohr2angstroms
            elif bohr.match(line):
                self.set_units('Bohr')
                self.input_units_to_au = 1.0

            # handle no_reorient
            elif orient.match(line):
                self.fix_orientation(True)

            # handle no_com
            elif com.match(line):
                self.PYmove_to_com = False

            # handle symmetry
            elif symmetry.match(line):
                self.PYsymmetry_from_input = symmetry.match(line).group(1).lower()

            # handle variables
            elif variable.match(line):
                vname = variable.match(line).group(1).upper()
                vval = float(variable.match(line).group(2))
                tda = 360.0 * math.atan(math.sqrt(2)) / math.pi
                self.geometry_variables['%s' % vname] = tda if vname == 'TDA' else vval

            # handle charge and multiplicity
            elif cgmp.match(line):
                tempCharge = int(cgmp.match(line).group(1))
                tempMultiplicity = int(cgmp.match(line).group(2))

                if ifrag == 0:
                    self.PYcharge_specified = True
                    self.PYmultiplicity_specified = True
                    self.PYmolecular_charge = tempCharge
                    self.PYmultiplicity = tempMultiplicity
                self.fragment_charges.append(tempCharge)
                self.fragment_multiplicities.append(tempMultiplicity)

            # handle fragment markers and default fragment cgmp
            elif frag.match(line):
                try:
                    self.fragment_charges[ifrag]
                except:
                    self.fragment_charges.append(0)
                    self.fragment_multiplicities.append(1)
                ifrag += 1
                glines.append(line)

            elif atom.match(line.split()[0].strip()):
                glines.append(line)
            else:
                raise ValidationError('Molecule::create_molecule_from_string: Unidentifiable line in geometry specification: %s' % (line))

        # catch last default fragment cgmp
        try:
            self.fragment_charges[ifrag]
        except:
            self.fragment_charges.append(0)
            self.fragment_multiplicities.append(1)

        # Now go through the rest of the lines looking for fragment markers
        ifrag = 0
        iatom = 0
        tempfrag = []
        atomSym = ""
        atomLabel = ""
        zmatrix = False
        for line in glines:

            # handle fragment markers
            if frag.match(line):
                ifrag += 1
                self.fragments.append([tempfrag[0], tempfrag[-1]])
                self.fragment_types.append('Real')
                tempfrag = []

            # handle atom markers
            else:
                entries = re.split(r'\s+|\s*,\s*', line.strip())
                atomm = atom.match(line.split()[0].strip().upper())
                atomLabel = atomm.group('label')
                atomSym = atomm.group('symbol')

                # We don't know whether the @C or Gh(C) notation matched. Do a quick check.
                ghostAtom = False if (atomm.group('gh1') is None and atomm.group('gh2') is None) else True

                # Check that the atom symbol is valid
                if not atomSym in el2z:
                    raise ValidationError('Molecule::create_molecule_from_string: Illegal atom symbol in geometry specification: %s' % (atomSym))

                zVal = el2z[atomSym]
                atomMass = el2mass[atomSym] if atomm.group('mass') is None else float(atomm.group('mass'))
                charge = float(zVal)
                if ghostAtom:
                    zVal = 0
                    charge = 0.0

                # handle cartesians
                if len(entries) == 4:
                    tempfrag.append(iatom)
                    xval = self.get_coord_value(entries[1])
                    yval = self.get_coord_value(entries[2])
                    zval = self.get_coord_value(entries[3])
                    self.full_atoms.append(CartesianEntry(iatom, zVal, charge, \
                        atomMass, atomSym, atomLabel, \
                        xval, yval, zval))

                # handle first line of Zmat
                elif len(entries) == 1:
                    zmatrix = True
                    tempfrag.append(iatom)
                    self.full_atoms.append(ZMatrixEntry(iatom, zVal, charge, \
                        atomMass, atomSym, atomLabel))

                # handle second line of Zmat
                elif len(entries) == 3:
                    zmatrix = True
                    tempfrag.append(iatom)

                    rTo = self.get_anchor_atom(entries[1], line)
                    if rTo >= iatom:
                        raise ValidationError("Molecule::create_molecule_from_string: Error on geometry input line %s. Atom %s has not been defined yet.\n" % (line, entries[1]))
                    rval = self.get_coord_value(entries[2])

                    if self.full_atoms[rTo].symbol() == 'X':
                        rval.set_fixed(True)

                    self.full_atoms.append(ZMatrixEntry(iatom, zVal, charge, \
                        atomMass, atomSym, atomLabel, \
                        self.full_atoms[rTo], rval))

                # handle third line of Zmat
                elif len(entries) == 5:
                    zmatrix = True
                    tempfrag.append(iatom)

                    rTo = self.get_anchor_atom(entries[1], line)
                    if rTo >= iatom:
                        raise ValidationError("Molecule::create_molecule_from_string: Error on geometry input line %s. Atom %s has not been defined yet.\n" % (line, entries[1]))
                    aTo = self.get_anchor_atom(entries[3], line)
                    if aTo >= iatom:
                        raise ValidationError("Molecule::create_molecule_from_string: Error on geometry input line %s. Atom %s has not been defined yet.\n" % (line, entries[3]))
                    if aTo == rTo:
                        raise ValidationError("Molecule::create_molecule_from_string: Atom used multiple times on line %s." % (line))
                    rval = self.get_coord_value(entries[2])
                    aval = self.get_coord_value(entries[4])

                    if self.full_atoms[rTo].symbol() == 'X':
                        rval.set_fixed(True)
                    if self.full_atoms[aTo].symbol() == 'X':
                        aval.set_fixed(True)

                    self.full_atoms.append(ZMatrixEntry(iatom, zVal, charge, \
                        atomMass, atomSym, atomLabel, \
                        self.full_atoms[rTo], rval, \
                        self.full_atoms[aTo], aval))

                # handle fourth line of Zmat
                elif len(entries) == 7:
                    zmatrix = True
                    tempfrag.append(iatom)

                    rTo = self.get_anchor_atom(entries[1], line)
                    if rTo >= iatom:
                        raise ValidationError("Molecule::create_molecule_from_string: Error on geometry input line %s. Atom %s has not been defined yet.\n" % (line, entries[1]))
                    aTo = self.get_anchor_atom(entries[3], line)
                    if aTo >= iatom:
                        raise ValidationError("Molecule::create_molecule_from_string: Error on geometry input line %s. Atom %s has not been defined yet.\n" % (line, entries[3]))
                    dTo = self.get_anchor_atom(entries[5], line)
                    if dTo >= iatom:
                        raise ValidationError("Molecule::create_molecule_from_string: Error on geometry input line %s. Atom %s has not been defined yet.\n" % (line, entries[5]))
                    if aTo == rTo or rTo == dTo or aTo == dTo:  # for you star wars fans
                        raise ValidationError("Molecule::create_molecule_from_string: Atom used multiple times on line %s" % (line))

                    rval = self.get_coord_value(entries[2])
                    aval = self.get_coord_value(entries[4])
                    dval = self.get_coord_value(entries[6])

                    if self.full_atoms[rTo].symbol() == 'X':
                        rval.set_fixed(True)
                    if self.full_atoms[aTo].symbol() == 'X':
                        aval.set_fixed(True)
                    if self.full_atoms[dTo].symbol() == 'X':
                        dval.set_fixed(True)

                    self.full_atoms.append(ZMatrixEntry(iatom, zVal, charge, \
                        atomMass, atomSym, atomLabel, \
                        self.full_atoms[rTo], rval, \
                        self.full_atoms[aTo], aval, \
                        self.full_atoms[dTo], dval))

                else:
                    raise ValidationError('Molecule::create_molecule_from_string: Illegal geometry specification line : %s. \
                        You should provide either Z-Matrix or Cartesian input' % (line))

                iatom += 1

        self.fragments.append([tempfrag[0], tempfrag[-1]])
        self.fragment_types.append('Real')
        self.set_has_zmatrix(zmatrix)

    def init_with_checkpoint(self, chkpt):
        """ **NYI** Pull information from the *chkpt* object passed
        (method name in libmints is init_with_chkpt)

        """
        raise FeatureNotImplemented('Molecule::init_with_checkpoint')  # FINAL

    def init_with_io(self, psio):
        """ **NYI** Pull information from a chkpt object created from psio
        (method name in libmints is init_with_psio)

        """
        raise FeatureNotImplemented('Molecule::init_with_io')  # FINAL

    @classmethod
    def init_with_xyz(cls, xyzfilename):
        """Pull information from an XYZ file. No fragment or chg/mult info detected.

        >>> H2O = qcdb.Molecule.init_with_xyz('h2o.xyz')

        """
        instance = cls()
        instance.lock_frame = False

        try:
            infile = open(xyzfilename, 'r')
        except IOError:
            raise ValidationError("""Molecule::init_with_xyz: given filename '%s' does not exist.""" % (xyzfilename))
        if os.stat(xyzfilename).st_size == 0:
            raise ValidationError("""Molecule::init_with_xyz: given filename '%s' is blank.""" % (xyzfilename))
        text = infile.readlines()

        xyz1 = re.compile(r"^\s*(\d+)\s*(bohr|au)?\s*$", re.IGNORECASE)
        xyzN = re.compile(r"(?:\s*)([A-Z](?:[a-z])?)(?:\s+)(-?\d+\.\d+)(?:\s+)(-?\d+\.\d+)(?:\s+)(-?\d+\.\d+)(?:\s*)", re.IGNORECASE)

        # Try to match the first line
        if xyz1.match(text[0]):
            fileNatom = int(xyz1.match(text[0]).group(1))
            if xyz1.match(text[0]).group(2) == None:
                fileUnits = 'Angstrom'
            else:
                fileUnits = 'Bohr'
        else:
            raise ValidationError("Molecule::init_with_xyz: Malformed first line\n%s" % (text[0]))

        # Skip the second line

        # Next line begins the useful information.
        for i in range(fileNatom):
            try:
                if xyzN.match(text[2 + i]):

                    fileAtom = xyzN.match(text[2 + i]).group(1).upper()
                    fileX = float(xyzN.match(text[2 + i]).group(2))
                    fileY = float(xyzN.match(text[2 + i]).group(3))
                    fileZ = float(xyzN.match(text[2 + i]).group(4))

                    # Coordinates in Molecule must be bohr.
                    if fileUnits == 'Angstrom':
                        fileX /= psi_bohr2angstroms
                        fileY /= psi_bohr2angstroms
                        fileZ /= psi_bohr2angstroms

                    # Check that the atom symbol is valid
                    if not fileAtom in el2z:
                        raise ValidationError('Molecule::init_with_xyz: Illegal atom symbol in geometry specification: %s' % (atomSym))

                    # Add it to the molecule.
                    instance.add_atom(el2z[fileAtom], fileX, fileY, fileZ, fileAtom, el2mass[fileAtom])

                else:
                    raise ValidationError("Molecule::init_with_xyz: Malformed atom information line %d." % (i + 3))
            except IndexError:
                raise ValidationError("Molecule::init_with_xyz: Expected atom in file at line %d.\n%s" % (i + 3, text[i + 2]))

        # We need to make 1 fragment with all atoms
        instance.fragments.append([0, fileNatom - 1])
        instance.fragment_types.append('Real')
        instance.fragment_charges.append(0)
        instance.fragment_multiplicities.append(1)
        # Set the units to bohr since we did the conversion above, if needed.
        instance.PYunits = 'Bohr'
        instance.input_units_to_au = 1.0

        instance.update_geometry()
        return instance

    def clone(self):
        """Returns new, independent Molecule object.

        >>> dimer = H2OH2O.clone()

        """
        return copy.deepcopy(self)

    # <<< Methods for Printing >>>

    def print_out(self):
        """Print the molecule.
        (method name in libmints is print)

        >>> H2OH2O.print_out()
        Geometry (in Angstrom), charge = -2, multiplicity = 3:
           Center              X                  Y                   Z
        ------------   -----------------  -----------------  -----------------
               O         -1.551007000000    -0.114520000000     0.000000000000
               H         -1.934259000000     0.762503000000     0.000000000000
               H         -0.599677000000     0.040712000000     0.000000000000
               O          1.350625000000     0.111469000000     0.000000000000
               H          1.680398000000    -0.373741000000    -0.758561000000
               H          1.680398000000    -0.373741000000     0.758561000000

        """
        text = ""
        if self.natom():
            if self.pg:
                text += """    Molecular point group: %s\n""" % (self.pg.symbol())
            if self.full_pg:
                text += """    Full point group: %s\n\n""" % (self.get_full_point_group())
            text += """    Geometry (in %s), charge = %d, multiplicity = %d:\n\n""" % \
                ('Angstrom' if self.units() == 'Angstrom' else 'Bohr', self.molecular_charge(), self.multiplicity())
            text += """       Center              X                  Y                   Z       \n"""
            text += """    ------------   -----------------  -----------------  -----------------\n"""

            for i in range(self.natom()):
                geom = self.atoms[i].compute()
                text += """    %8s%4s """ % (self.symbol(i), "" if self.Z(i) else "(Gh)")
                for j in range(3):
                    text += """  %17.12f""" % (geom[j])
                text += "\n"
            # TODO if (Process::environment.options.get_int("PRINT") > 2) {
            text += "\n"
            for i in range(self.natom()):
                text += """    %8s\n""" % (self.label(i))
                for bas in self.atoms[i].basissets().keys():
                    text += """              %-15s %-20s""" % (bas,
                        self.atoms[i].basissets()[bas])
                    if bas in self.atoms[i].shells():
                        text += """%s""" % (self.atoms[i].shells()[bas])
                    text += '\n'
            text += "\n"
        else:
            text += "  No atoms in this molecule.\n"
        print(text)
        # TODO outfile

    def print_out_in_bohr(self):
        """Print the molecule in Bohr. Same as :py:func:`print_out` only in Bohr.
        (method name in libmints is print_in_bohr)

        """
        text = ""
        if self.natom():
            if self.pg:
                text += """    Molecular point group: %s\n""" % (self.pg.symbol())
            if self.full_pg:
                text += """    Full point group: %s\n\n""" % (self.get_full_point_group())
            text += """    Geometry (in %s), charge = %d, multiplicity = %d:\n\n""" % \
                ('Bohr', self.molecular_charge(), self.multiplicity())
            text += """       Center              X                  Y                   Z       \n"""
            text += """    ------------   -----------------  -----------------  -----------------\n"""

            for i in range(self.natom()):
                text += """    %8s%4s """ % (self.symbol(i), "" if self.Z(i) else "(Gh)")
                for j in range(3):
                    text += """  %17.12f""" % (self.xyz(i, j))
                text += "\n"
            text += "\n"
        else:
            text += "  No atoms in this molecule.\n"
        print(text)
        # TODO outfile

    def print_out_in_angstrom(self):
        """Print the molecule in Angstroms. Same as :py:func:`print_out` only always in Angstroms.
        (method name in libmints is print_in_angstrom)

        """
        text = ""
        if self.natom():
            if self.pg:
                text += """    Molecular point group: %s\n""" % (self.pg.symbol())
            if self.full_pg:
                text += """    Full point group: %s\n\n""" % (self.get_full_point_group())
            text += """    Geometry (in %s), charge = %d, multiplicity = %d:\n\n""" % \
                ('Angstrom', self.molecular_charge(), self.multiplicity())
            text += """       Center              X                  Y                   Z       \n"""
            text += """    ------------   -----------------  -----------------  -----------------\n"""

            for i in range(self.natom()):
                text += """    %8s%4s """ % (self.symbol(i), "" if self.Z(i) else "(Gh)")
                for j in range(3):
                    text += """  %17.12f""" % (self.xyz(i, j) * psi_bohr2angstroms)
                text += "\n"
            text += "\n"
        else:
            text += "  No atoms in this molecule.\n"
        print(text)
        # TODO outfile

    def print_full(self):
        """Print full atom list. Same as :py:func:`print_out` only displays dummy atoms.

        """
        text = ""
        if self.natom():
            if self.pg:
                text += """    Molecular point group: %s\n""" % (self.pg.symbol())
            if self.full_pg:
                text += """    Full point group: %s\n\n""" % (self.get_full_point_group())
            text += """    Geometry (in %s), charge = %d, multiplicity = %d:\n\n""" % \
                (self.units(), self.molecular_charge(), self.multiplicity())
            text += """       Center              X                  Y                   Z       \n"""
            text += """    ------------   -----------------  -----------------  -----------------\n"""

            for i in range(self.nallatom()):
                geom = self.full_atoms[i].compute()
                text += """    %8s%4s """ % (self.fsymbol(i), "" if self.fZ(i) else "(Gh)")
                for j in range(3):
                    text += """  %17.12f""" % (geom[j])
                text += "\n"
            text += "\n"
        else:
            text += "  No atoms in this molecule.\n"
        print(text)
        # TODO outfile

    def print_in_input_format(self):
        """Print the molecule in the same format that the user provided.
        """
        text = ""
        if self.nallatom():
            text += "    Geometry (in %s), charge = %d, multiplicity = %d:\n\n" % \
                    ("Angstrom" if self.units() == 'Angstrom' else "Bohr",
                    self.molecular_charge(), self.multiplicity())
            for i in range(self.nallatom()):
                if self.fZ(i) or self.fsymbol(i) == "X":
                    text += "    %-8s" % (self.fsymbol(i))
                else:
                    text += "    %-8s" % ("Gh(" + self.fsymbol(i) + ")")
                text += self.full_atoms[i].print_in_input_format()
            text += "\n"
            if len(self.geometry_variables):
                for vb, val in self.geometry_variables.items():
                    text += """    %-10s=%16.10f\n""" % (vb, val)
                text += "\n"

        print(text)
        # TODO outfile

    def everything(self):
        """Quick print of class data"""
        text = """  ==> qcdb Molecule %s <==\n\n""" % (self.name())
        text += """  Natom         %d\t\tNallatom       %d\n""" % (self.natom(), self.nallatom())
        text += """  charge        %d\t\tspecified?     %s\n""" % (self.molecular_charge(), self.charge_specified())
        text += """  multiplicity  %d\t\tspecified?     %s\n""" % (self.multiplicity(), self.multiplicity_specified())
        text += """  units         %s\tconversion     %f\n""" % (self.units(), self.input_units_to_au)
        text += """  DOcom?        %s\t\tDONTreorient?  %s\n""" % (self.PYmove_to_com, self.orientation_fixed())
        text += """  reinterpret?  %s\t\tlock_frame?    %s\n""" % (self.PYreinterpret_coordentries, self.lock_frame)
        text += """  input symm    %s\n""" % (self.symmetry_from_input())
        text += """  Nfragments    %d\t\tNactive        %d\n""" % (self.nfragments(), self.nactive_fragments())
        text += """  zmat?         %s\n""" % (self.has_zmatrix())
        print(text)

    def create_psi4_string_from_molecule(self, force_c1=False):
        """Regenerates a input file molecule specification string from the
        current state of the Molecule. Contains geometry info,
        fragmentation, charges and multiplicities, and any frame
        restriction.
        """
        text = ""
        if self.nallatom():

            # append units and any other non-default molecule keywords
            text += "    units %-s\n" % ("Angstrom" if self.units() == 'Angstrom' else "Bohr")
            if not self.PYmove_to_com:
                text += "    no_com\n"
            if self.PYfix_orientation:
                text += "    no_reorient\n"
            if force_c1:
                text += "    symmetry c1\n"

            # append atoms and coordentries and fragment separators with charge and multiplicity
            Pfr = 0
            for fr in range(self.nfragments()):
                if self.fragment_types[fr] == 'Absent' and not self.has_zmatrix():
                    continue
                text += "%s    %s%d %d\n" % (
                    "" if Pfr == 0 else "    --\n",
                    "#" if self.fragment_types[fr] == 'Ghost' or self.fragment_types[fr] == 'Absent' else "",
                    self.fragment_charges[fr], self.fragment_multiplicities[fr])
                Pfr += 1
                for at in range(self.fragments[fr][0], self.fragments[fr][1] + 1):
                    if self.fragment_types[fr] == 'Absent':
                        text += "    %-8s" % ("X")
                    elif self.fZ(at) or self.fsymbol(at) == "X":
                        text += "    %-8s" % (self.flabel(at))
                    else:
                        text += "    %-8s" % ("Gh(" + self.flabel(at) + ")")
                    text += "    %s" % (self.full_atoms[at].print_in_input_format())
            text += "\n"

            # append any coordinate variables
            if len(self.geometry_variables):
                for vb, val in self.geometry_variables.items():
                    text += """    %-10s=%16.10f\n""" % (vb, val)
                text += "\n"

        return text

    # <<< Involved Methods for Coordinates >>>

    def get_coord_value(self, vstr):
        """Attempts to interpret a string as a double, if not it assumes it's a variable.

        """
        vstr = vstr.upper()
        realNumber = re.compile(r"""[-+]?(?:(?:\d*\.\d+)|(?:\d+\.?))(?:[Ee][+-]?\d+)?""", re.VERBOSE)

        # handle number values
        if realNumber.match(vstr):
            return NumberValue(float(vstr))

        # handle variable values, whether defined or not
        else:
            if vstr == 'TDA':
                self.geometry_variables[vstr] = 360.0 * math.atan(math.sqrt(2)) / math.pi

            # handle negative variable values (ignore leading '-' and return minus the value)
            if vstr[0] == '-':
                self.all_variables.append(vstr[1:])
                return VariableValue(vstr[1:], self.geometry_variables, True)

            # handle normal variable values
            else:
                self.all_variables.append(vstr)
                return VariableValue(vstr, self.geometry_variables)

    def add_atom(self, Z, x, y, z, label="", mass=0.0, charge=0.0, lineno=-1):
        """Add an atom to the molecule
        *Z* atomic number
        *x* cartesian coordinate
        *y* cartesian coordinate
        *z* cartesian coordinate
        *symb* atomic symbol to use
        *mass* mass to use if non standard
        *charge* charge to use if non standard
        *lineno* line number when taken from a string

        """
        self.lock_frame = False

        if self.atom_at_position([x, y, z]) == -1:
            # Dummies go to full_atoms, ghosts need to go to both.
            self.full_atoms.append(CartesianEntry(self.nallatom(), Z, charge, mass, label, label, \
                NumberValue(x), NumberValue(y), NumberValue(z)))
            if label.upper() != 'X':
                self.atoms.append(self.full_atoms[-1])
        else:
            raise ValidationError("Molecule::add_atom: Adding atom on top of an existing atom.")

    def atom_entry(self, atom):
        """Returns the CoordEntry for an atom."""
        return self.atoms[atom]

    def atom_at_position(self, b, tol=0.05):
        """Tests to see of an atom is at the passed position *b* in Bohr with a tolerance *tol*.

        >>> print H2OH2O.atom_at_position([1.35*(1.0/psi_bohr2angstroms), 0.10*(1.0/psi_bohr2angstroms), 0.0*(1.0/psi_bohr2angstroms)])
        3

        """
        if len(b) != 3:
            raise ValidationError('Molecule::atom_at_position: Argument vector not of length 3\n')

        for at in range(self.natom()):
            a = self.xyz(at)
            if distance(b, a) < tol:
                return at
        return -1

    def is_variable(self, vstr):
        """Checks to see if the variable str is in the list, returns
        true if it is, and returns false if not.

        >>> H2OH2O.is_variable('R')
        False

        """
        #if self.all_variables
        #print 'vstr', vstr, 'all_variables', self.all_variables, (vstr.upper() in self.all_variables), '\n'
        return True if vstr.upper() in self.all_variables else False

    def get_variable(self, vstr):
        """Checks to see if the variable str is in the list, sets it to
        val and returns true if it is, and returns false if not.

        """
        vstr = vstr.upper()
        try:
            return self["geometry_variables"][vstr]
        except KeyError:
            raise ValidationError('Molecule::get_variable: Geometry variable %s not known.\n' % (vstr))

    def set_variable(self, vstr, val):
        """Assigns the value val to the variable labelled string in the
        list of geometry variables. Also calls update_geometry()

        """
        self.lock_frame = False
        self["geometry_variables"][vstr.upper()] = val
        print("""Setting geometry variable %s to %f""" % (vstr.upper(), val))
        try:
            self.update_geometry()
        except IncompleteAtomError:
            # Update geometry might have added some atoms, delete them to be safe.
            self.atoms = []
        # TODO outfile

    def __setattr__(self, name, value):
        """Function to overload setting attributes to allow geometry
        variable assigment as if member data.

        """
    
        if "all_variables" not in self.keys():
            self["all_variables"] = []

        if name.upper() in list(self["all_variables"]):
            self.set_variable(name, value)
        else:
            self[name] = value

    def __getattr__(self, name):
        """Function to overload accessing attribute contents to allow
        retrivial geometry variable values as if member data.

        """

        if name.upper() in list(self["all_variables"]):
            return self.get_variable(name)
        elif name in list(self):
            return self[name]
        else:
            raise AttributeError


    def get_anchor_atom(self, vstr, line):
        """Attempts to interpret a string *vstr* as an atom specifier in
        a zmatrix. Takes the current *line* for error message printing.
        Returns the atom number (adjusted to zero-based counting).

        """
        integerNumber = re.compile(r"(-?\d+)", re.IGNORECASE)
        if integerNumber.match(vstr):
            # This is just a number, return it
            return int(vstr) - 1
        else:
            # Look to see if this string is known
            for i in range(self.nallatom()):
                if self.full_atoms[i].label() == vstr:
                    return i
            raise ValidationError("Molecule::get_anchor_atom: Illegal value %s in atom specification on line %s.\n" % (vstr, line))

    def geometry(self):
        """Returns the geometry in Bohr as a N X 3 array.

        >>> print H2OH2O.geometry()
        [[-2.930978460188563, -0.21641143673806384, 0.0], [-3.655219780069251, 1.4409218455037016, 0.0], [-1.1332252981904638, 0.0769345303220403, 0.0], [2.5523113582286716, 0.21064588230662976, 0.0], [3.175492014248769, -0.7062681346308132, -1.4334725450878665], [3.175492014248769, -0.7062681346308132, 1.4334725450878665]]

        """
        geom = []
        for at in range(self.natom()):
            geom.append([self.x(at), self.y(at), self.z(at)])
        return geom

    def full_geometry(self):
        """Returns the full (dummies included) geometry in Bohr as a N X 3 array.

        >>> print H2OH2O.full_geometry()
        [[-2.930978460188563, -0.21641143673806384, 0.0], [-3.655219780069251, 1.4409218455037016, 0.0], [-1.1332252981904638, 0.0769345303220403, 0.0], [0.0, 0.0, 0.0], [2.5523113582286716, 0.21064588230662976, 0.0], [3.175492014248769, -0.7062681346308132, -1.4334725450878665], [3.175492014248769, -0.7062681346308132, 1.4334725450878665]]

        """
        geom = []
        for at in range(self.nallatom()):
            geom.append([self.fx(at), self.fy(at), self.fz(at)])
        return geom

    def set_geometry(self, geom):
        """Sets the geometry, given a N X 3 array of coordinates *geom* in Bohr.

        >>> H2OH2O.set_geometry([[1,2,3],[4,5,6],[7,8,9],[-1,-2,-3],[-4,-5,-6],[-7,-8,-9]])

        """
        self.lock_frame = False
        for at in range(self.natom()):
            self.atoms[at].set_coordinates(geom[at][0] / self.input_units_to_au,
                                           geom[at][1] / self.input_units_to_au,
                                           geom[at][2] / self.input_units_to_au)

    def set_full_geometry(self, geom):
        """Sets the full geometry (dummies included), given a N X 3 array of coordinates *geom* in Bohr.

        >>> H2OH2O.set_full geometry([[1,2,3],[4,5,6],[7,8,9],[0,0,0],[-1,-2,-3],[-4,-5,-6],[-7,-8,-9]])

        """
        self.lock_frame = False
        for at in range(self.nallatom()):
            self.full_atoms[at].set_coordinates(geom[at][0] / self.input_units_to_au,
                                                geom[at][1] / self.input_units_to_au,
                                                geom[at][2] / self.input_units_to_au)

    def distance_matrix(self):
        """Computes a matrix depicting distances between atoms. Prints
        formatted and returns array.

        >>> H2OH2O.distance_matrix()
                Interatomic Distances (Angstroms)
                          [1]         [2]         [3]         [4]         [5]         [6]
          [1]         0.00000
          [2]         0.95711     0.00000
          [3]         0.96391     1.51726     0.00000
          [4]         2.91042     3.34878     1.95159     0.00000
          [5]         3.32935     3.86422     2.43843     0.95895     0.00000
          [6]         3.32935     3.86422     2.43843     0.95895     1.51712     0.00000

        """
        distm = [[None] * self.natom() for _ in range(self.natom())]
        text = "        Interatomic Distances (Angstroms)\n\n          "
        for i in range(self.natom()):
            text += '%11s ' % ('[' + str(i + 1) + ']')
        text += "\n"
        for i in range(self.natom()):
            text += '  %-8s ' % ('[' + str(i + 1) + ']')
            for j in range(self.natom()):
                if j > i:
                    continue
                if j == i:
                    text += '%10.5f  ' % (0.0)
                    distm[i][j] = 0.0
                else:
                    eij = sub(self.xyz(j), self.xyz(i))
                    dist = norm(eij) * psi_bohr2angstroms
                    text += '%10.5f  ' % (dist)
                    distm[i][j] = dist
                    distm[j][i] = dist
            text += "\n"
        text += "\n\n"
        print(text)
        return distm
        # TODO outfile

    def print_distances(self):
        """Print the geometrical parameters (distances) of the molecule.
        suspect libmints version actually prints Bohr.

        >>> print H2OH2O.print_distances()
        Interatomic Distances (Angstroms)
        Distance 1 to 2 0.957
        Distance 1 to 3 0.964
        Distance 1 to 4 2.910
        ...

        """
        text = "        Interatomic Distances (Angstroms)\n\n"
        for i in range(self.natom()):
            for j in range(i + 1, self.natom()):
                eij = sub(self.xyz(j), self.xyz(i))
                dist = norm(eij) * psi_bohr2angstroms
                text += "        Distance %d to %d %-8.3lf\n" % (i + 1, j + 1, dist)
        text += "\n\n"
        return text
        # TODO outfile

    def print_bond_angles(self):
        """Print the geometrical parameters (bond_angles) of the molecule.

        >>> print H2OH2O.print_bond_angles()
        Bond Angles (degrees)
        Angle 2-1-3:  104.337
        Angle 2-1-4:  109.152
        Angle 2-1-5:  117.387
        ...

        """
        text = "        Bond Angles (degrees)\n\n"
        for j in range(self.natom()):
            for i in range(self.natom()):
                if j == i:
                    continue
                for k in range(i + 1, self.natom()):
                    if j == k:
                        continue
                    eji = sub(self.xyz(i), self.xyz(j))
                    eji = normalize(eji)
                    ejk = sub(self.xyz(k), self.xyz(j))
                    ejk = normalize(ejk)
                    dotproduct = dot(eji, ejk)
                    phi = 180.0 * math.acos(dotproduct) / math.pi
                    text += "        Angle %d-%d-%d: %8.3lf\n" % (i + 1, j + 1, k + 1, phi)
        text += "\n\n"
        return text
        # TODO outfile

    def print_dihedrals(self):
        """Print the geometrical parameters (dihedrals) of the molecule.

        >>> print H2OH2O.print_dihedrals()
        Dihedral Angles (Degrees)
        Dihedral 1-2-3-4:  180.000
        Dihedral 1-2-3-5:  133.511
        Dihedral 1-2-3-6:  133.511
        ...

        """
        text = "        Dihedral Angles (Degrees)\n\n"
        for i in range(self.natom()):
            for j in range(self.natom()):
                if i == j:
                    continue
                for k in range(self.natom()):
                    if i == k or j == k:
                        continue
                    for l in range(self.natom()):
                        if i == l or j == l or k == l:
                            continue
                        eij = sub(self.xyz(j), self.xyz(i))
                        eij = normalize(eij)
                        ejk = sub(self.xyz(k), self.xyz(j))
                        ejk = normalize(ejk)
                        ekl = sub(self.xyz(l), self.xyz(k))
                        ekl = normalize(ekl)
                        # Compute angle ijk
                        angleijk = math.acos(dot(scale(eij, -1.0), ejk))
                        # Compute angle jkl
                        anglejkl = math.acos(dot(scale(ejk, -1.0), ekl))
                        # compute term1 (eij x ejk)
                        term1 = cross(eij, ejk)
                        # compute term2 (ejk x ekl)
                        term2 = cross(ejk, ekl)
                        numerator = dot(term1, term2)
                        denominator = math.sin(angleijk) * math.sin(anglejkl)
                        try:
                            costau = numerator / denominator
                        except ZeroDivisionError:
                            costau = 0.0
                        if costau > 1.00 and costau < 1.000001:
                            costau = 1.00
                        if costau < -1.00 and costau > -1.000001:
                            costau = -1.00
                        tau = 180.0 * math.acos(costau) / math.pi
                        text += "        Dihedral %d-%d-%d-%d: %8.3lf\n" % (i + 1, j + 1, k + 1, l + 1, tau)
        text += "\n\n"
        return text
        # TODO outfile

    def print_out_of_planes(self):
        """Print the geometrical parameters (out_of_planes) of the molecule.

        >>> print H2OH2O.print_out_of_planes()
        Out-Of-Plane Angles (Degrees)
        Out-of-plane 1-2-3-4:    0.000
        Out-of-plane 1-2-3-5:   -7.373
        Out-of-plane 1-2-3-6:    7.373
        ...

        """
        text = "        Out-Of-Plane Angles (Degrees)\n\n"
        for i in range(self.natom()):
            for j in range(self.natom()):
                if i == j:
                    continue
                for k in range(self.natom()):
                    if i == k or j == k:
                        continue
                    for l in range(self.natom()):
                        if i == l or j == l or k == l:
                            continue
                        # Compute vectors we need first
                        elj = sub(self.xyz(j), self.xyz(l))
                        elj = normalize(elj)
                        elk = sub(self.xyz(k), self.xyz(l))
                        elk = normalize(elk)
                        eli = sub(self.xyz(i), self.xyz(l))
                        eli = normalize(eli)
                        # Denominator
                        denominator = math.sin(math.acos(dot(elj, elk)))
                        # Numerator
                        eljxelk = cross(elj, elk)
                        numerator = dot(eljxelk, eli)
                        # compute angle
                        try:
                            sinetheta = numerator / denominator
                        except ZeroDivisionError:
                            sinetheta = 0.0
                        if sinetheta > 1.00:
                            sinetheta = 1.000
                        if sinetheta < -1.00:
                            sinetheta = -1.000
                        theta = 180.0 * math.asin(sinetheta) / math.pi
                        text += "        Out-of-plane %d-%d-%d-%d: %8.3lf\n" % (i + 1, j + 1, k + 1, l + 1, theta)
        text += "\n\n"
        return text
        # TODO outfile

    def reinterpret_coordentry(self, rc):
        """Do we reinterpret coordentries during a call to update_geometry?
        (method name in libmints is set_reinterpret_coordentry)

        """
        self.PYreinterpret_coordentries = rc

    def reinterpret_coordentries(self):
        """Reinterpret the fragments for reals/ghosts and build the atom list.

        """
        self.atoms = []
        for item in self.full_atoms:
            item.invalidate()

        temp_charge = self.PYmolecular_charge
        temp_multiplicity = self.PYmultiplicity
        self.PYmolecular_charge = 0
        self.PYmultiplicity = 1

        for fr in range(self.nfragments()):
            if self.fragment_types[fr] == 'Absent':
                continue

            if self.fragment_types[fr] == 'Real':
                self.PYmolecular_charge += self.fragment_charges[fr]
                self.PYmultiplicity += self.fragment_multiplicities[fr] - 1

            for at in range(self.fragments[fr][0], self.fragments[fr][1] + 1):
                self.full_atoms[at].compute()
                self.full_atoms[at].set_ghosted(self.fragment_types[fr] == 'Ghost')
                if self.full_atoms[at].symbol() != 'X':
                    self.atoms.append(self.full_atoms[at])

        # TODO: This is a hack to ensure that set_multiplicity and set_molecular_charge
        #   work for single-fragment molecules.
        if self.nfragments() < 2:
            self.PYmolecular_charge = temp_charge
            self.PYmultiplicity = temp_multiplicity

    def update_geometry(self):
        """Updates the geometry, by (re)interpreting the string used to
        create the molecule, and the current values of the variables.
        The atoms list is cleared, and then rebuilt by this routine.
        This function must be called after first instantiation of Molecule.

        >>> H2 = qcdb.Molecule("H\\nH 1 0.74\\n")
        >>> print H2.natom()
        0
        >>> H2.update_geometry()
        >>> print H2.natom()
        2

        """
        if self.nfragments() == 0:
            raise ValidationError("Molecule::update_geometry: There are no fragments in this molecule.")

        # Idempotence condition
        if self.lock_frame:
            return

        #print("beginning update_geometry:")
        #self.print_full()
        if self.PYreinterpret_coordentries:
            self.reinterpret_coordentries()
        #print("after reinterpret_coordentries:")
        #self.print_full()

        if self.PYmove_to_com:
            self.move_to_com()
        #print "after com:"
        #self.print_full()

        # If the no_reorient command was given, don't reorient
        if not self.PYfix_orientation:
            # Now we need to rotate the geometry to its symmetry frame
            # to align the axes correctly for the point group
            # symmetry_frame looks for the highest point group so that we can align
            # the molecule according to its actual symmetry, rather than the symmetry
            # the the user might have provided.
            frame = self.symmetry_frame()
            self.rotate_full(frame)
            #print "after rotate:"
            #self.print_full()

        # Recompute point group of the molecule, so the symmetry info is updated to the new frame
        self.set_point_group(self.find_point_group())
        self.set_full_point_group()

        # Disabling symmetrize for now if orientation is fixed, as it is not
        #   correct.  We may want to fix this in the future, but in some cases of
        #   finite-differences the set geometry is not totally symmetric anyway.
        # Symmetrize the molecule to remove any noise
        self.symmetrize()
        #print "after symmetry:"
        #self.print_full()

        self.lock_frame = True

    # <<< Methods for Miscellaneous >>>

    def clear(self):
        """Zero it out."""
        self.lock_frame = False
        self.atoms = []
        self.full_atoms = []

    def nuclear_repulsion_energy(self):
        """Computes nuclear repulsion energy.

        >>> print H2OH2O.nuclear_repulsion_energy()
        36.6628478528

        """
        e = 0.0
        for at1 in range(self.natom()):
            for at2 in range(self.natom()):
                if at2 < at1:
                    Zi = self.Z(at1)
                    Zj = self.Z(at2)
                    dist = distance(self.xyz(at1), self.xyz(at2))
                    e += Zi * Zj / dist
        return e

    def nuclear_repulsion_energy_deriv1(self):
        """Computes nuclear repulsion energy derivatives

        >>> print H2OH2O.nuclear_repulsion_energy_deriv1()
        [[3.9020946901323774, 2.76201566471991, 0.0], [1.3172905807089021, -2.3486366050337293, 0.0], [-1.8107598525022435, -0.32511212499256564, 0.0], [-1.217656141385739, -2.6120090867576717, 0.0], [-1.0954846384766488, 1.2618710760320282, 2.1130743287465603], [-1.0954846384766488, 1.2618710760320282, -2.1130743287465603]]

        """
        de = []
        for i in range(self.natom()):
            entry = [0.0, 0.0, 0.0]
            for j in range(self.natom()):
                if i != j:
                    temp = distance(self.xyz(i), self.xyz(j)) ** 3.0
                    Zi = self.Z(i)
                    Zj = self.Z(j)
                    entry[0] -= (self.x(i) - self.x(j)) * Zi * Zj / temp
                    entry[1] -= (self.y(i) - self.y(j)) * Zi * Zj / temp
                    entry[2] -= (self.z(i) - self.z(j)) * Zi * Zj / temp
            de.append(entry)
        return de

    def nuclear_repulsion_energy_deriv2(self):
        """ **NYI** Computes nuclear repulsion energy second derivatives"""
        raise FeatureNotImplemented('Molecule::nuclear_repulsion_energy_deriv2')  # FINAL

    def set_basis_all_atoms(self, name, role="BASIS"):
        """Assigns basis *name* to all atoms."""
        uc = name.upper()
        if uc in ['SPECIAL', 'GENERAL', 'CUSTOM']:
            # These aren't really basis set specifications, just return.
            return
        for atom in self.full_atoms:
            atom.set_basisset(name, role)

    def set_basis_by_symbol(self, symbol, name, role="BASIS"):
        """Assigns basis *name* to all *symbol* atoms."""
        for atom in self.full_atoms:
            if symbol.upper() == atom.symbol():
                atom.set_basisset(name, role)

    def clear_basis_all_atoms(self):
        """Remove all basis information from atoms."""
        for atom in self.full_atoms:
            atom.PYbasissets = OrderedDict()

    def set_basis_by_number(self, number, name, role="BASIS"):
        """Assigns basis *name* to atom number *number* (0-indexed, excludes dummies)."""
        # change from libmints to 0-indexing and to real/ghost numbering, dummies not included (libmints >= error)
        if number >= self.natom():
            raise ValidationError("Molecule::set_basis_by_number: Basis specified for atom %d, but there are only %d atoms in this molecule." % \
                (number, self.natom()))
        self.atoms[number].set_basisset(name, role)

    def set_basis_by_label(self, label, name, role="BASIS"):
        """Assigns basis *name* to all atoms with *label*."""
        for atom in self.full_atoms:
            if label.upper() == atom.label():
                atom.set_basisset(name, role)

    def set_shell_by_number(self, number, bshash, role="BASIS"):
        """Assigns BasisSet *bshash* to atom number *number* (0-indexed, excludes dummies)."""
        self.lock_frame = False
        if number >= self.natom():
            raise ValidationError("Molecule::set_shell_by_number: Basis specified for atom %d, but there are only %d atoms in this molecule." % \
                (number, self.natom()))
        self.atoms[number].set_shell(bshash, role)

    def nfrozen_core(self, depth=False):
        """Number of frozen core for molecule given freezing state.

        >>> print H2OH2O.nfrozen_core()
        2

        """
        if depth == False or depth.upper() == 'FALSE':
            return 0

        elif depth == True or depth.upper() == 'TRUE':
            # Freeze the number of core electrons corresponding to the
            # nearest previous noble gas atom.  This means that the 4p block
            # will still have 3d electrons active.  Alkali earth atoms will
            # have one valence electron in this scheme.
            nfzc = 0
            for A in range(self.natom()):
                if self.Z(A) > 2:
                    nfzc += 1
                if self.Z(A) > 10:
                    nfzc += 4
                if self.Z(A) > 18:
                    nfzc += 4
                if self.Z(A) > 36:
                    nfzc += 9
                if self.Z(A) > 54:
                    nfzc += 9
                if self.Z(A) > 86:
                    nfzc += 16
                if self.Z(A) > 108:
                    raise ValidationError("Molecule::nfrozen_core: Invalid atomic number")
            return nfzc

        else:
            raise ValidationError("Molecule::nfrozen_core: Frozen core '%s' is not supported, options are {true, false}." % (depth))

    # <<< Involved Methods for Frame >>>

    def translate(self, r):
        """Translates molecule by r.

        >>> H2OH2O.translate([1.0, 1.0, 0.0])

        """
        temp = [None, None, None]
        for at in range(self.nallatom()):
            temp = scale(self.full_atoms[at].compute(), self.input_units_to_au)
            temp = add(temp, r)
            temp = scale(temp, 1.0 / self.input_units_to_au)
            self.full_atoms[at].set_coordinates(temp[0], temp[1], temp[2])

    def center_of_mass(self):
        """Computes center of mass of molecule (does not translate molecule).

        >>> H2OH2O.center_of_mass()
        [-0.12442647346606871, 0.00038657002584110707, 0.0]

        """
        ret = [0.0, 0.0, 0.0]
        total_m = 0.0

        for at in range(self.natom()):
            m = self.mass(at)
            ret = add(ret, scale(self.xyz(at), m))
            total_m += m

        ret = scale(ret, 1.0 / total_m)
        return ret

    def move_to_com(self):
        """Moves molecule to center of mass

        """
        com = scale(self.center_of_mass(), -1.0)
        self.translate(com)

    def set_com_fixed(self, _fix=True):
        """ **NYI** Fix the center of mass at its current frame.
        Not used in libmints so not implemented.

        """
        raise FeatureNotImplemented('Molecule::set_com_fixed')  # FINAL

    def inertia_tensor(self):
        """Compute inertia tensor.

        >>> print H2OH2O.inertia_tensor()
        [[8.704574864178731, -8.828375721817082, 0.0], [-8.828375721817082, 280.82861714077666, 0.0], [0.0, 0.0, 281.249500988553]]

        """
        tensor = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]

        for i in range(self.natom()):
            # I(alpha, alpha)
            tensor[0][0] += self.mass(i) * (self.y(i) * self.y(i) + self.z(i) * self.z(i))
            tensor[1][1] += self.mass(i) * (self.x(i) * self.x(i) + self.z(i) * self.z(i))
            tensor[2][2] += self.mass(i) * (self.x(i) * self.x(i) + self.y(i) * self.y(i))

            # I(alpha, beta)
            tensor[0][1] -= self.mass(i) * self.x(i) * self.y(i)
            tensor[0][2] -= self.mass(i) * self.x(i) * self.z(i)
            tensor[1][2] -= self.mass(i) * self.y(i) * self.z(i)

        # mirror
        tensor[1][0] = tensor[0][1]
        tensor[2][0] = tensor[0][2]
        tensor[2][1] = tensor[1][2]

        # Check the elements for zero and make them a hard zero.
        for i in range(3):
            for j in range(3):
                if math.fabs(tensor[i][j]) < ZERO:
                    tensor[i][j] = 0.0
        return tensor

    def rotational_constants(self, tol=FULL_PG_TOL):
        """Compute the rotational constants and return them in wavenumbers"""

        evals, evecs = diagonalize3x3symmat(self.inertia_tensor())
        evals = sorted(evals)
        isLinear = True if evals[0] < ZERO else False
        isOne = True if evals[1] < ZERO else False
        isAtom = True if evals[2] < ZERO else False

        im_amuA = psi_bohr2angstroms * psi_bohr2angstroms
        im_ghz = psi_h * psi_na * 1E14 / (8.0 * math.pi * math.pi * psi_bohr2angstroms * psi_bohr2angstroms)
        im_mhz = im_ghz * 1000
        im_cm = im_ghz * 1E7 / psi_c

        text = "        Moments of Inertia and Rotational Constants\n\n"
        text += '  %-12s    %3s %16.8f    %3s %16.8f    %3s %16.8f\n' % \
            ('[amu B^2]', 'I_A', evals[0], 'I_B', evals[1], 'I_C', evals[2])
        text += '  %-12s    %3s %16.8f    %3s %16.8f    %3s %16.8f\n' % \
            ('[amu A^2]', 'I_A', evals[0] * im_amuA, 'I_B', evals[1] * im_amuA, 'I_C', evals[2] * im_amuA)
        text += '  %-12s    %3s %16s    %3s %16s    %3s %16s\n' % ('[GHz]', \
            'A', '%16.8f' % (im_ghz / evals[0]) if not isLinear else '*****', \
            'B', '%16.8f' % (im_ghz / evals[1]) if not isOne else '*****', \
            'C', '%16.8f' % (im_ghz / evals[2]) if not isAtom else '*****')
        text += '  %-12s    %3s %16s    %3s %16s    %3s %16s\n' % ('[MHz]', \
            'A', '%16.8f' % (im_mhz / evals[0]) if not isLinear else '*****', \
            'B', '%16.8f' % (im_mhz / evals[1]) if not isOne else '*****', \
            'C', '%16.8f' % (im_mhz / evals[2]) if not isAtom else '*****')
        text += '  %-12s    %3s %16s    %3s %16s    %3s %16s\n' % ('[cm^-1]', \
            'A', '%16.8f' % (im_cm / evals[0]) if not isLinear else '*****', \
            'B', '%16.8f' % (im_cm / evals[1]) if not isOne else '*****', \
            'C', '%16.8f' % (im_cm / evals[2]) if not isAtom else '*****')
        print(text)
        # TODO outfile
        rot_const = []
        rot_const.append(im_cm / evals[0] if not isLinear else None)
        rot_const.append(im_cm / evals[1] if not isOne else None)
        rot_const.append(im_cm / evals[2] if not isAtom else None)
        return rot_const

    def rotor_type(self, tol=FULL_PG_TOL):
        """Returns the rotor type.

        >>> H2OH2O.rotor_type()
        RT_ASYMMETRIC_TOP

        """
        rot_const = self.rotational_constants()
        for i in range(3):
            if rot_const[i] == None:
                rot_const[i] = 0.0

        # Determine degeneracy of rotational constants.
        degen = 0
        for i in range(2):
            for j in range(i + 1, 3):
                if degen >= 2:
                    continue
                rabs = math.fabs(rot_const[i] - rot_const[j])
                tmp = rot_const[i] if rot_const[i] > rot_const[j] else rot_const[j]
                if rabs > ZERO:
                    rel = rabs / tmp
                else:
                    rel = 0.0
                if rel < tol:
                    degen += 1
        #print "\tDegeneracy is %d\n" % (degen)

        # Determine rotor type
        if self.natom() == 1:
            rotor_type = 'RT_ATOM'
        elif rot_const[0] == 0.0:
            rotor_type = 'RT_LINEAR'          # 0  <  IB == IC      inf > B == C
        elif degen == 2:
            rotor_type = 'RT_SPHERICAL_TOP'   # IA == IB == IC       A == B == C
        elif degen == 1:
            rotor_type = 'RT_SYMMETRIC_TOP'   # IA <  IB == IC       A >  B == C --or--
                                              # IA == IB <  IC       A == B >  C
        else:
            rotor_type = 'RT_ASYMMETRIC_TOP'  # IA <  IB <  IC       A >  B >  C
        return rotor_type

    def rotate(self, R):
        """Rotates the molecule using rotation matrix *R*.

        >>> H2OH2O.rotate([[0,-1,0],[-1,0,0],[0,0,1]])

        """
        new_geom = zero(3, self.natom())
        geom = self.geometry()
        new_geom = mult(geom, R)
        self.set_geometry(new_geom)

    def rotate_full(self, R):
        """Rotates the full molecule using rotation matrix *R*.

        >>> H2OH2O.rotate_full([[0,-1,0],[-1,0,0],[0,0,1]])

        """
        new_geom = zero(3, self.nallatom())
        geom = self.full_geometry()
        new_geom = mult(geom, R)
        self.set_full_geometry(new_geom)

    def orientation_fixed(self):
        """Get whether or not orientation is fixed.

        >>> H2OH2O.orientation_fixed()
        True

        """
        return self.PYfix_orientation

    def fix_orientation(self, _fix=True):
        """Fix the orientation at its current frame
        (method name in libmints is set_orientation_fixed)

        """
        if _fix:
            self.PYfix_orientation = True  # tells update_geometry() not to change orientation
            # Compute original cartesian coordinates - code coped from update_geometry()
            self.atoms = []
            for item in self.full_atoms:
                item.invalidate()

            for fr in range(self.nfragments()):
                for at in range(self.fragments[fr][0], self.fragments[fr][1] + 1):
                    self.full_atoms[at].compute()
                    self.full_atoms[at].set_ghosted(self.fragment_types[fr] == 'Ghost')
                    if self.full_atoms[at].symbol() != 'X':
                        self.atoms.append(self.full_atoms[at])
        else:  # release orientation to be free
            self.PYfix_orientation = False

    # <<< Methods for Saving >>>

    def save_string_xyz(self, save_ghosts=True):
        """Save a string for a XYZ-style file.

        >>> H2OH2O.save_string_xyz()
        6
        _
         O   -1.551007000000   -0.114520000000    0.000000000000
         H   -1.934259000000    0.762503000000    0.000000000000
         H   -0.599677000000    0.040712000000    0.000000000000
         O    1.350625000000    0.111469000000    0.000000000000
         H    1.680398000000   -0.373741000000   -0.758561000000
         H    1.680398000000   -0.373741000000    0.758561000000

        """
        factor = 1.0 if self.PYunits == 'Angstrom' else psi_bohr2angstroms

        N = self.natom()
        if not save_ghosts:
            N = 0
            for i in range(self.natom()):
                if self.Z(i):
                    N += 1
        text = "%d\n\n" % (N)

        for i in range(self.natom()):
            [x, y, z] = self.atoms[i].compute()
            if save_ghosts or self.Z(i):
                text += '%2s %17.12f %17.12f %17.12f\n' % ((self.symbol(i) if self.Z(i) else "Gh"), \
                    x * factor, y * factor, z * factor)
        return text

    def save_xyz(self, filename, save_ghosts=True):
        """Save an XYZ file.

        >>> H2OH2O.save_xyz('h2o.xyz')

        """
        outfile = open(filename, 'w')
        outfile.write(self.save_string_xyz(save_ghosts))
        outfile.close()

    def save_to_checkpoint(self, chkpt, prefix=""):
        """ **NYI** Save information to checkpoint file
        (method name in libmints is save_to_chkpt)

        """
        raise FeatureNotImplemented('Molecule::save_to_checkpoint')  # FINAL

    # <<< Methods for Symmetry >>>

    def has_symmetry_element(self, op, tol=DEFAULT_SYM_TOL):
        """ **NYI** Whether molecule satisfies the vector symmetry
        operation *op*. Not used by libmints.

        """
        raise FeatureNotImplemented('Molecule::has_symmetry_element')  # FINAL
        for i in range(self.natom()):
            result = naivemult(self.xyz(i), op)
            atom = self.atom_at_position(result, tol)

            if atom != -1:
                if not self.atoms[atom].is_equivalent_to(self.atoms[i]):
                    return False
            else:
                return False
        return True

    def point_group(self):
        """Returns the point group (object) if set"""
        if self.pg is None:
            raise ValidationError("Molecule::point_group: Molecular point group has not been set.")
        return self.pg

    def set_point_group(self, pg):
        """Set the point group to object *pg* """
        self.pg = pg
        # Call this here, the programmer will forget to call it, as I have many times.
        self.form_symmetry_information()

    def set_full_point_group(self, tol=FULL_PG_TOL):
        """Determine and set FULL point group. self.PYfull_pg_n is highest
        order n in Cn. 0 for atoms or infinity.

        """
        verbose = 1  # TODO
        # Get cartesian geometry and put COM at origin
        geom = self.geometry()
        com = self.center_of_mass()
        for at in range(self.natom()):
            geom[at][0] += -com[0]
            geom[at][1] += -com[1]
            geom[at][2] += -com[2]

        # Get rotor type
        rotor = self.rotor_type(tol)
        if verbose > 2:
            print("""  Rotor type                       : %s""" % (rotor))

        # Get the D2h point group from Jet and Ed's code: c1 ci c2 cs d2 c2v c2h d2h
        #   and ignore the user-specified subgroup in this case.
        pg = self.find_highest_point_group(tol)
        d2h_subgroup = pg.symbol()
        if verbose > 2:
            print("""  D2h_subgroup                     : %s""" % (self.point_group().symbol()))

        # Check inversion
        v3_zero = [0.0, 0.0, 0.0]
        op_i = self.has_inversion(v3_zero, tol)
        if verbose > 2:
            print("""  Inversion symmetry               : %s""" % ('yes' if op_i else 'no'))

        x_axis = [1, 0, 0]
        y_axis = [0, 1, 0]
        z_axis = [0, 0, 1]
        rot_axis = [0.0, 0.0, 0.0]

        if rotor == 'RT_ATOM':  # atoms
            self.full_pg = 'ATOM'
            self.PYfull_pg_n = 0

        elif rotor == 'RT_LINEAR':  # linear molecules
            self.full_pg = 'D_inf_h' if op_i else 'C_inf_v'
            self.PYfull_pg_n = 0

        elif rotor == 'RT_SPHERICAL_TOP':  # spherical tops
            if not op_i:  # The only spherical top without inversion is Td.
                self.full_pg = 'Td'
                self.PYfull_pg_n = 3
            else:  # Oh or Ih ?
                # Oh has a S4 and should be oriented properly already.
                test_mat = matrix_3d_rotation(geom, z_axis, math.pi / 2.0, True)
                op_symm = equal_but_for_row_order(geom, test_mat, tol)
                if verbose > 2:
                    print("""  S4z                              : %s""" % ('yes' if op_symm else 'no'))

                if op_symm:
                    self.full_pg = 'Oh'
                    self.PYfull_pg_n = 4
                else:
                    self.full_pg = 'Ih'
                    self.PYfull_pg_n = 5

        elif rotor == 'RT_ASYMMETRIC_TOP':  # asymmetric tops cannot exceed D2h, right?

            if d2h_subgroup == 'c1':
                self.full_pg = 'C1'
                self.PYfull_pg_n = 1

            elif d2h_subgroup == 'ci':
                self.full_pg = 'Ci'
                self.PYfull_pg_n = 1

            elif d2h_subgroup == 'c2':
                self.full_pg = 'Cn'
                self.PYfull_pg_n = 2

            elif d2h_subgroup == 'cs':
                self.full_pg = 'Cs'
                self.PYfull_pg_n = 1

            elif d2h_subgroup == 'd2':
                self.full_pg = 'Dn'
                self.PYfull_pg_n = 2

            elif d2h_subgroup == 'c2v':
                self.full_pg = 'Cnv'
                self.PYfull_pg_n = 2

            elif d2h_subgroup == 'c2h':
                self.full_pg = 'Cnh'
                self.PYfull_pg_n = 2

            elif d2h_subgroup == 'd2h':
                self.full_pg = 'Dnh'
                self.PYfull_pg_n = 2

            else:
                print("""  Warning: Cannot determine point group.""")

        elif rotor in ['RT_SYMMETRIC_TOP', 'RT_PROLATE_SYMMETRIC_TOP', 'RT_OBLATE_SYMMETRIC_TOP']:

            # Find principal axis that is unique and make it z-axis.
            It = self.inertia_tensor()
            I_evals, I_evecs = diagonalize3x3symmat(It)
            ev_list = list(zip(I_evals, transpose(I_evecs)))  # eigenvectors are cols of I_evecs
            ev_list.sort(key=lambda tup: tup[0], reverse=False)
            I_evals, I_evecs = zip(*ev_list)  # sorted eigenvectors are now rows of I_evecs
            if verbose > 2:
                print("""  I_evals: %15.10lf %15.10lf %15.10lf""" % (I_evals[0], I_evals[1], I_evals[2]))

            unique_axis = 1
            if abs(I_evals[0] - I_evals[1]) < tol:
                unique_axis = 2
            elif abs(I_evals[1] - I_evals[2]) < tol:
                unique_axis = 0

            # Compute angle between unique axis and the z-axis
            old_axis = I_evecs[unique_axis]

            ddot = dot(z_axis, old_axis)
            if abs(ddot - 1) < 1.0e-10:
                phi = 0.0
            elif abs(ddot + 1) < 1.0e-10:
                phi = math.pi
            else:
                phi = math.acos(ddot)

            # Rotate geometry to put unique axis on the z-axis, if it isn't already.
            if abs(phi) > 1.0e-14:
                rot_axis = cross(z_axis, old_axis)  # right order?
                test_mat = matrix_3d_rotation(geom, rot_axis, phi, False)
                if verbose > 2:
                    print("""  Rotating by %lf to get principal axis on z-axis ...""" % (phi))
                geom = [row[:] for row in test_mat]

            if verbose > 2:
                print("""  Geometry to analyze - principal axis on z-axis:""")
                for at in range(self.natom()):
                    print("""%20.15lf %20.15lf %20.15lf""" % (geom[at][0], geom[at][1], geom[at][2]))
                print('\n')

            # Determine order Cn and Sn of principal axis.
            Cn_z = matrix_3d_rotation_Cn(geom, z_axis, False, tol)
            if verbose > 2:
                print("""  Highest rotation axis (Cn_z)     : %d""" % (Cn_z))

            Sn_z = matrix_3d_rotation_Cn(geom, z_axis, True, tol)
            if verbose > 2:
                print("""  Rotation axis (Sn_z)             : %d""" % (Sn_z))

            # Check for sigma_h (xy plane).
            op_sigma_h = False
            for at in range(self.natom()):
                if abs(geom[at][2]) < tol:
                    continue  # atom is in xy plane
                else:
                    test_atom = [geom[at][0], geom[at][1], -1 * geom[at][2]]
                    if not atom_present_in_geom(geom, test_atom, tol):
                        break
            else:
                op_sigma_h = True
            if verbose > 2:
                print("""  sigma_h                          : %s""" % ('yes' if op_sigma_h else 'no'))

            # Rotate one off-axis atom to the yz plane and check for sigma_v's.
            for at in range(self.natom()):
                dist_from_z = math.sqrt(geom[at][0] * geom[at][0] + geom[at][1] * geom[at][1])
                if abs(dist_from_z) > tol:
                    pivot_atom_i = at
                    break

            if pivot_atom_i == self.natom():  # needs to be in else clause?
                raise ValidationError("Molecule::set_full_point_group: Not a linear molecule but could not find off-axis atom.")

            # Rotate around z-axis to put pivot atom in the yz plane
            xy_point = normalize([geom[pivot_atom_i][0], geom[pivot_atom_i][1], 0])
            ddot = dot(y_axis, xy_point)
            if abs(ddot - 1) < 1.0e-10:
                phi = 0.0
            elif abs(ddot + 1) < 1.0e-10:
                phi = math.pi
            else:
                phi = math.acos(ddot)

            is_D = False
            if abs(phi) > 1.0e-14:
                test_mat = matrix_3d_rotation(geom, z_axis, phi, False)
                if verbose > 2:
                    print("""  Rotating by %8.3e to get atom %d in yz-plane ...""" % (phi, pivot_atom_i + 1))
                geom = [row[:] for row in test_mat]

            # Check for sigma_v (yz plane).
            op_sigma_v = False
            for at in range(self.natom()):
                if abs(geom[at][0]) < tol:
                    continue  # atom is in yz plane
                else:
                    test_atom = [-1 * geom[at][0], geom[at][1], geom[at][2]]
                    if not atom_present_in_geom(geom, test_atom, tol):
                        break
            else:
            #if at == self.natom():
                op_sigma_v = True
            if verbose > 2:
                print("""  sigma_v                          : %s""" % ('yes' if op_sigma_v else 'no'))

                print("""  geom to analyze - one atom in yz plane:""")
                for at in range(self.natom()):
                    print("""%20.15lf %20.15lf %20.15lf""" % (geom[at][0], geom[at][1], geom[at][2]))
                print('\n')

            # Check for perpendicular C2's.
            # Loop through pairs of atoms to find c2 axis candidates.
            for i in range(self.natom()):
                A = [geom[i][0], geom[i][1], geom[i][2]]
                AdotA = dot(A, A)
                for j in range(at):
                    if self.Z(at) != self.Z(j):
                        continue  # ensure same atomic number

                    B = [geom[j][0], geom[j][1], geom[j][2]]  # ensure same distance from com
                    if abs(AdotA - dot(B, B)) > 1.0e-6:
                        continue  # loose check

                    # Use sum of atom vectors as axis if not 0.
                    axis = add(A, B)
                    if norm(axis) < 1.0e-12:
                        continue
                    axis = normalize(axis)

                    # Check if axis is perpendicular to z-axis.
                    if abs(dot(axis, z_axis)) > 1.0e-6:
                        continue

                    # Do the thorough check for C2.
                    if matrix_3d_rotation_Cn(geom, axis, False, tol, 2) == 2:
                        is_D = True
            if verbose > 2:
                print("""  perp. C2's                       : %s""" % ('yes' if is_D else 'no'))

            # Now assign point groups!  Sn first.
            if Sn_z == 2 * Cn_z and not is_D:
                self.full_pg = 'Sn'
                self.PYfull_pg_n = Sn_z
                return

            if is_D:  # has perpendicular C2's
                if op_sigma_h and op_sigma_v:  # Dnh : Cn, nC2, sigma_h, nSigma_v
                    self.full_pg = 'Dnh'
                    self.PYfull_pg_n = Cn_z

                elif Sn_z == 2 * Cn_z:  # Dnd : Cn, nC2, S2n axis coincident with Cn
                    self.full_pg = 'Dnd'
                    self.PYfull_pg_n = Cn_z

                else:  # Dn : Cn, nC2
                    self.full_pg = 'Dn'
                    self.PYfull_pg_n = Cn_z

            else:  # lacks perpendicular C2's
                if op_sigma_h and Sn_z == Cn_z:  # Cnh : Cn, sigma_h, Sn coincident with Cn
                    self.full_pg = 'Cnh'
                    self.PYfull_pg_n = Cn_z

                elif op_sigma_v:  # Cnv : Cn, nCv
                    self.full_pg = 'Cnv'
                    self.PYfull_pg_n = Cn_z

                else:  # Cn  : Cn
                    self.full_pg = 'Cn'
                    self.PYfull_pg_n = Cn_z
        return

    def has_inversion(self, origin, tol=DEFAULT_SYM_TOL):
        """Does the molecule have an inversion center at origin

        """
        for i in range(self.natom()):
            inverted = sub(origin, sub(self.xyz(i), origin))
            atom = self.atom_at_position(inverted, tol)
            if atom < 0 or not self.atoms[atom].is_equivalent_to(self.atoms[i]):
                return False
        return True

    def is_plane(self, origin, uperp, tol=DEFAULT_SYM_TOL):
        """Is a plane?

        """
        for i in range(self.natom()):
            A = sub(self.xyz(i), origin)
            Apar = scale(uperp, dot(uperp, A))
            Aperp = sub(A, Apar)
            A = add(sub(Aperp, Apar), origin)
            atom = self.atom_at_position(A, tol)
            if atom < 0 or not self.atoms[atom].is_equivalent_to(self.atoms[i]):
                return False
        return True

    def is_axis(self, origin, axis, order, tol=DEFAULT_SYM_TOL):
        """Is *axis* an axis of order *order* with respect to *origin*?

        """
        for i in range(self.natom()):
            A = sub(self.xyz(i), origin)
            for j in range(1, order):
                R = A
                R = rotate(R, j * 2.0 * math.pi / order, axis)
                R = add(R, origin)
                atom = self.atom_at_position(R, tol)
                if atom < 0 or not self.atoms[atom].is_equivalent_to(self.atoms[i]):
                    return False
        return True

    def is_linear_planar(self, tol=DEFAULT_SYM_TOL):
        """Is the molecule linear, or planar?

        >>> print H2OH2O.is_linear_planar()
        (False, False)

        """
        linear = None
        planar = None

        if self.natom() < 3:
            linear = True
            planar = True
            return linear, planar

        # find three atoms not on the same line
        A = self.xyz(0)
        B = self.xyz(1)
        BA = sub(B, A)
        BA = normalize(BA)
        CA = [None, None, None]

        min_BAdotCA = 1.0
        for i in range(2, self.natom()):
            tmp = sub(self.xyz(i), A)
            tmp = normalize(tmp)
            if math.fabs(dot(BA, tmp)) < min_BAdotCA:
                CA = copy.deepcopy(tmp)
                min_BAdotCA = math.fabs(dot(BA, tmp))
        if min_BAdotCA >= 1.0 - tol:
            linear = True
            planar = True
            return linear, planar

        linear = False
        if self.natom() < 4:
            planar = True
            return linear, planar

        # check for nontrivial planar molecules
        BAxCA = normalize(cross(BA, CA))
        for i in range(2, self.natom()):
            tmp = sub(self.xyz(i), A)
            if math.fabs(dot(tmp, BAxCA)) > tol:
                planar = False
                return linear, planar
        planar = True
        return linear, planar

    @staticmethod
    def like_world_axis(axis, worldxaxis, worldyaxis, worldzaxis):
        """Returns which worldaxis *axis* most overlaps with.
        Inverts axis when indicated.

        """
        like = None
        xlikeness = math.fabs(dot(axis, worldxaxis))
        ylikeness = math.fabs(dot(axis, worldyaxis))
        zlikeness = math.fabs(dot(axis, worldzaxis))

        if (xlikeness - ylikeness) > 1.0E-12 and (xlikeness - zlikeness) > 1.0E-12:
            like = 'XAxis'
            if dot(axis, worldxaxis) < 0:
                axis = scale(axis, -1.0)
        elif (ylikeness - zlikeness) > 1.0E-12:
            like = 'YAxis'
            if dot(axis, worldyaxis) < 0:
                axis = scale(axis, -1.0)
        else:
            like = 'ZAxis'
            if dot(axis, worldzaxis) < 0:
                axis = scale(axis, -1.0)
        return like, axis

    def find_point_group(self, tol=DEFAULT_SYM_TOL):
        """Find computational molecular point group, user can override
        this with the "symmetry" keyword. Result is highest D2h subgroup
        attendant on molecule and allowed by the user.

        """
        pg = self.find_highest_point_group(tol)  # D2h subgroup
        user = self.symmetry_from_input()

        if user is not None:
            # Need to handle the cases that the user only provides C2, C2v, C2h, Cs.
            #   These point groups need directionality.

            # Did the user provide directionality? If they did, the last letter would be x, y, or z
            #   Directionality given, assume the user is smart enough to know what they're doing.
            user_specified_direction = True if user[-1] in ['X', 'x', 'Y', 'y', 'Z', 'z'] else False

            if self.symmetry_from_input() != pg.symbol():
                user = PointGroup(self.symmetry_from_input())

                if user_specified_direction:
                    # Assume the user knows what they're doing.
                    #   Make sure user is subgroup of pg
                    if (pg.bits() & user.bits()) != user.bits():
                        raise ValidationError("Molecule::find_point_group: User specified point group (%s) is not a subgroup of the highest detected point group (%s)" % (PointGroup.bits_to_full_name(user.bits()), PointGroup.bits_to_full_name(pg.bits())))

                else:
                    similars, count = similar(user.bits())
                    found = False
                    for typ in range(count):
                        # If what the user specified and the similar type
                        #   matches the full point group we've got a match
                        if (similars[typ] & pg.bits()) == similars[typ]:
                            found = True
                            break

                    if found:
                        # Construct a point group object using the found similar
                        user = PointGroup(similars[typ])

                    else:
                        raise ValidationError("Molecule::find_point_group: User specified point group (%s) is not a subgroup of the highest detected point group (%s). If this is because the symmetry increased, try to start the calculation again from the last geometry, after checking any symmetry-dependent input, such as DOCC." % (PointGroup.bits_to_full_name(user.bits()), PointGroup.bits_to_full_name(pg.bits())))

                # If we make it here, what the user specified is good.
                pg = user

        return pg

    def reset_point_group(self, pgname):
        """Override symmetry from outside the molecule string"""
        self.PYsymmetry_from_input = pgname.lower()
        self.set_point_group(self.find_point_group())

    def find_highest_point_group(self, tol=DEFAULT_SYM_TOL):
        """Find the highest D2h point group from Jet and Ed's code: c1
        ci c2 cs d2 c2v c2h d2h. Ignore the user-specified subgroup in
        this case.

        """
        pg_bits = 0

        # The order of the next 2 arrays MUST match!
        symm_bit = [
            SymmOps['C2_z'],
            SymmOps['C2_y'],
            SymmOps['C2_x'],
            SymmOps['i'],
            SymmOps['Sigma_xy'],
            SymmOps['Sigma_xz'],
            SymmOps['Sigma_yz']]

        symm_func = [
            SymmetryOperation.c2_z,
            SymmetryOperation.c2_y,
            SymmetryOperation.c2_x,
            SymmetryOperation.i,
            SymmetryOperation.sigma_xy,
            SymmetryOperation.sigma_xz,
            SymmetryOperation.sigma_yz]

        symop = SymmetryOperation()
        matching_atom = -1

        # Only needs to detect the 8 symmetry operations
        for g in range(7):

            # Call the function pointer
            symm_func[g](symop)
            found = True

            for at in range(self.natom()):
                op = [symop[0][0], symop[1][1], symop[2][2]]
                pos = naivemult(self.xyz(at), op)

                matching_atom = self.atom_at_position(pos, tol)
                if matching_atom >= 0:
                    if not self.atoms[at].is_equivalent_to(self.atoms[matching_atom]):
                        found = False
                        break
                else:
                    found = False
                    break

            if found:
                pg_bits |= symm_bit[g]

        return PointGroup(pg_bits)

    def symmetry_frame(self, tol=DEFAULT_SYM_TOL):
        """Determine symmetry reference frame. If noreorient is not set,
        this is the rotation matrix applied to the geometry in update_geometry.

        >>> print H2OH2O.symmetry_frame()
        [[1.0, -0.0, 0.0], [0.0, 1.0, 0.0], [0.0, -0.0, 1.0]]

        """
        com = self.center_of_mass()
        worldxaxis = [1.0, 0.0, 0.0]
        worldyaxis = [0.0, 1.0, 0.0]
        worldzaxis = [0.0, 0.0, 1.0]

        sigma = [0.0, 0.0, 0.0]
        sigmav = [0.0, 0.0, 0.0]
        c2axis = [0.0, 0.0, 0.0]
        c2axisperp = [0.0, 0.0, 0.0]

        linear, planar = self.is_linear_planar(tol)
        have_inversion = self.has_inversion(com, tol)

        # check for C2 axis
        have_c2axis = False
        if self.natom() < 2:
            have_c2axis = True
            c2axis = [0.0, 0.0, 1.0]

        elif linear:
            have_c2axis = True
            c2axis = sub(self.xyz(1), self.xyz(0))
            c2axis = normalize(c2axis)

        elif planar and have_inversion:
            # there is a c2 axis that won't be found using the usual
            #   algorithm. find two noncolinear atom-atom vectors (we know
            #   that linear == 0)
            BA = sub(self.xyz(1), self.xyz(0))
            BA = normalize(BA)
            for i in range(2, self.natom()):
                CA = sub(self.xyz(i), self.xyz(0))
                CA = normalize(CA)
                BAxCA = cross(BA, CA)
                if norm(BAxCA) > tol:
                    have_c2axis = True
                    BAxCA = normalize(BAxCA)
                    c2axis = copy.deepcopy(BAxCA)
                    break

        else:
            # loop through pairs of atoms to find c2 axis candidates
            for i in range(self.natom()):
                A = sub(self.xyz(i), com)
                AdotA = dot(A, A)
                for j in range(i + 1):
                    # the atoms must be identical
                    if not self.atoms[i].is_equivalent_to(self.atoms[j]):
                        continue
                    B = sub(self.xyz(j), com)
                    # the atoms must be the same distance from the com
                    if math.fabs(AdotA - dot(B, B)) > tol:
                        continue
                    axis = add(A, B)
                    # atoms colinear with the com don't work
                    if norm(axis) < tol:
                        continue
                    axis = normalize(axis)
                    if self.is_axis(com, axis, 2, tol):
                        have_c2axis = True
                        c2axis = copy.deepcopy(axis)
                        break
                else:
                    continue
                break

        # symmframe found c2axis
        c2like = 'ZAxis'
        if have_c2axis:
            # try to make the sign of the axis correspond to one of the world axes
            c2like, c2axis = self.like_world_axis(c2axis, worldxaxis, worldyaxis, worldzaxis)

        # check for c2 axis perp to first c2 axis
        have_c2axisperp = False
        if have_c2axis:
            if self.natom() < 2:
                have_c2axisperp = True
                c2axisperp = [1.0, 0.0, 0.0]

            elif linear:
                if have_inversion:
                    have_c2axisperp = True
                    c2axisperp = perp_unit(c2axis, [0.0, 0.0, 1.0])

            else:
                # loop through paris of atoms to find c2 axis candidates
                for i in range(self.natom()):
                    A = sub(self.xyz(i), com)
                    AdotA = dot(A, A)
                    for j in range(i):
                        # the atoms must be identical
                        if not self.atoms[i].is_equivalent_to(self.atoms[j]):
                            continue
                        B = sub(self.xyz(j), com)
                        # the atoms must be the same distance from the com
                        if math.fabs(AdotA - dot(B, B)) > tol:
                            continue
                        axis = add(A, B)
                        # atoms colinear with the com don't work
                        if norm(axis) < tol:
                            continue
                        axis = normalize(axis)
                        # if axis is not perp continue
                        if math.fabs(dot(axis, c2axis)) > tol:
                            continue
                        if self.is_axis(com, axis, 2, tol):
                            have_c2axisperp = True
                            c2axisperp = copy.deepcopy(axis)
                        break
                    else:
                        continue
                    break

        # symmframe found c2axisperp
        if have_c2axisperp:
            # try to make the sign of the axis correspond to one of the world axes
            c2perplike, c2axisperp = self.like_world_axis(c2axisperp, worldxaxis, worldyaxis, worldzaxis)

            # try to make c2axis the z axis
            if c2perplike == 'ZAxis':
                tmpv = copy.deepcopy(c2axisperp)
                c2axisperp = copy.deepcopy(c2axis)
                c2axis = copy.deepcopy(tmpv)
                c2perplike = c2like
                c2like = 'ZAxis'

            if c2like != 'ZAxis':
                if c2like == 'XAxis':
                    c2axis = cross(c2axis, c2axisperp)
                else:
                    c2axis = cross(c2axisperp, c2axis)
                c2like, c2axis = self.like_world_axis(c2axis, worldxaxis, worldyaxis, worldzaxis)

            # try to make c2axisperplike the x axis
            if c2perplike == 'YAxis':
                c2axisperp = cross(c2axisperp, c2axis)
                c2perplike, c2axisperp = self.like_world_axis(c2axisperp, worldxaxis, worldyaxis, worldzaxis)

        # Check for vertical plane
        have_sigmav = False
        if have_c2axis:
            if self.natom() < 2:
                have_sigmav = True
                sigmav = copy.deepcopy(c2axisperp)

            elif linear:
                have_sigmav = True
                if have_c2axisperp:
                    sigmav = copy.deepcopy(c2axisperp)

                else:
                    sigmav = perp_unit(c2axis, [0.0, 0.0, 1.0])
            else:
                # loop through pairs of atoms to find sigma v plane candidates
                for i in range(self.natom()):
                    A = sub(self.xyz(i), com)
                    AdotA = dot(A, A)
                    # the second atom can equal i because i might be in the plane
                    for j in range(i + 1):
                        # the atoms must be identical
                        if not self.atoms[i].is_equivalent_to(self.atoms[j]):
                            continue
                        B = sub(self.xyz(j), com)
                        # the atoms must be the same distance from the com
                        if math.fabs(AdotA - dot(B, B)) > tol:
                            continue
                        inplane = add(B, A)
                        norm_inplane = norm(inplane)
                        if norm_inplane < tol:
                            continue
                        inplane = scale(inplane, 1.0 / norm_inplane)
                        perp = cross(c2axis, inplane)
                        norm_perp = norm(perp)
                        if norm_perp < tol:
                            continue
                        perp = scale(perp, 1.0 / norm_perp)
                        if self.is_plane(com, perp, tol):
                            have_sigmav = True
                            sigmav = copy.deepcopy(perp)
                            break
                    else:
                        continue
                    break

        # symmframe found sigmav
        if have_sigmav:
            # try to make the sign of the oop vec correspond to one of the world axes
            sigmavlike, sigmav = self.like_world_axis(sigmav, worldxaxis, worldyaxis, worldzaxis)

            # Choose sigmav to be the world x axis, if possible
            if c2like == 'ZAxis' and sigmavlike == 'YAxis':
                sigmav = cross(sigmav, c2axis)
            elif c2like == 'YAxis' and sigmavlike == 'ZAxis':
                sigmav = cross(c2axis, sigmav)

        # under certain conditions i need to know if there is any sigma plane
        have_sigma = False
        if not have_inversion and not have_c2axis:
            if planar:
                # find two noncolinear atom-atom vectors
                # we know that linear==0 since !have_c2axis
                BA = sub(self.xyz(1), self.xyz(0))
                BA = normalize(BA)
                for i in range(2, self.natom()):
                    CA = sub(self.xyz(i), self.xyz(0))
                    CA = normalize(CA)
                    BAxCA = cross(BA, CA)
                    if norm(BAxCA) > tol:
                        have_sigma = True
                        BAxCA = normalize(BAxCA)
                        sigma = copy.deepcopy(BAxCA)
                        break
            else:
                # loop through pairs of atoms to contruct trial planes
                for i in range(self.natom()):
                    A = sub(self.xyz(i), com)
                    AdotA = dot(A, A)
                    for j in range(i):
                        # the atoms must be identical
                        if not self.atoms[i].is_equivalent_to(self.atoms[j]):
                            continue
                        B = sub(self.xyz(j), com)
                        BdotB = dot(B, B)
                        # the atoms must be the same distance from the com
                        if math.fabs(AdotA - BdotB) > tol:
                            continue
                        perp = sub(B, A)
                        norm_perp = norm(perp)
                        if norm_perp < tol:
                            continue
                        perp = scale(perp, 1.0 / norm_perp)
                        if self.is_plane(com, perp, tol):
                            have_sigma = True
                            sigma = copy.deepcopy(perp)
                            break
                    else:
                        continue
                    break

        # foundsigma
        if have_sigma:
            # try to make the sign of the oop vec correspond to one of the world axes
            xlikeness = math.fabs(dot(sigma, worldxaxis))
            ylikeness = math.fabs(dot(sigma, worldyaxis))
            zlikeness = math.fabs(dot(sigma, worldzaxis))

            if xlikeness > ylikeness and xlikeness > zlikeness:
                if dot(sigma, worldxaxis) < 0:
                    sigma = scale(sigma, -1.0)
            elif ylikeness > zlikeness:
                if dot(sigma, worldyaxis) < 0:
                    sigma = scale(sigma, -1.0)
            else:
                if dot(sigma, worldzaxis) < 0:
                    sigma = scale(sigma, -1.0)

        # Find the three axes for the symmetry frame
        xaxis = copy.deepcopy(worldxaxis)
        zaxis = copy.deepcopy(worldzaxis)
        if have_c2axis:
            zaxis = copy.deepcopy(c2axis)
            if have_sigmav:
                xaxis = copy.deepcopy(sigmav)
            elif have_c2axisperp:
                xaxis = copy.deepcopy(c2axisperp)
            else:
                # any axis orthogonal to the zaxis will do
                xaxis = perp_unit(zaxis, zaxis)
        elif have_sigma:
            zaxis = copy.deepcopy(sigma)
            xaxis = perp_unit(zaxis, zaxis)

        # Clean up our z axis
        if math.fabs(zaxis[0]) < NOISY_ZERO:
            zaxis[0] = 0.0
        if math.fabs(zaxis[1]) < NOISY_ZERO:
            zaxis[1] = 0.0
        if math.fabs(zaxis[2]) < NOISY_ZERO:
            zaxis[2] = 0.0

        # Clean up our x axis
        if math.fabs(xaxis[0]) < NOISY_ZERO:
            xaxis[0] = 0.0
        if math.fabs(xaxis[1]) < NOISY_ZERO:
            xaxis[1] = 0.0
        if math.fabs(xaxis[2]) < NOISY_ZERO:
            xaxis[2] = 0.0

        # the y is then -x cross z
        yaxis = scale(cross(xaxis, zaxis), -1.0)

        #print "xaxis %20.14lf %20.14lf %20.14lf" % (xaxis[0], xaxis[1], xaxis[2])
        #print "yaxis %20.14lf %20.14lf %20.14lf" % (yaxis[0], yaxis[1], yaxis[2])
        #print "zaxis %20.14lf %20.14lf %20.14lf" % (zaxis[0], zaxis[1], zaxis[2])

        frame = zero(3, 3)
        for i in range(3):
            frame[i][0] = xaxis[i]
            frame[i][1] = yaxis[i]
            frame[i][2] = zaxis[i]
        return frame

    def release_symmetry_information(self):
        """Release symmetry information"""
        self.PYnunique = 0
        self.nequiv = 0
        self.PYatom_to_unique = 0
        self.equiv = 0

    def form_symmetry_information(self, tol=DEFAULT_SYM_TOL):
        """Initialize molecular specific symmetry information.
        Uses the point group object obtain by calling point_group()

        """
        if self.equiv:
            self.release_symmetry_information()

        if self.natom() == 0:
            self.PYnunique = 0
            self.nequiv = 0
            self.PYatom_to_unique = 0
            self.equiv = 0
            print("""No atoms detected, returning\n""")
            return

        self.nequiv = []
        self.PYatom_to_unique = [0] * self.natom()
        self.equiv = []

        if self.point_group().symbol() == 'c1':
            self.PYnunique = self.natom()
            for at in range(self.natom()):
                self.nequiv.append(1)
                self.PYatom_to_unique[at] = at
                self.equiv.append([at])
            return

        # The first atom is always unique
        self.PYnunique = 1
        self.nequiv.append(1)
        self.PYatom_to_unique[0] = 0
        self.equiv.append([0])

        ct = self.point_group().char_table()
        so = SymmetryOperation()
        np = [0.0, 0.0, 0.0]

        # Find the equivalent atoms
        for i in range(1, self.natom()):
            ac = self.xyz(i)
            i_is_unique = True
            i_equiv = 0

            # Apply all symmetry ops in the group to the atom
            for g in range(ct.order()):
                so = ct.symm_operation(g)
                for ii in range(3):
                    np[ii] = 0
                    for jj in range(3):
                        np[ii] += so[ii][jj] * ac[jj]

                # See if the transformed atom is equivalent to a unique atom
                for j in range(self.PYnunique):
                    unique = self.equiv[j][0]
                    aj = self.xyz(unique)
                    if distance(np, aj) < tol and \
                        self.Z(unique) == self.Z(i) and \
                        abs(self.mass(unique) - self.mass(i)) < tol:
                        i_is_unique = False
                        i_equiv = j
                        break

            if i_is_unique:
                self.nequiv.append(1)
                self.PYatom_to_unique[i] = self.PYnunique
                self.equiv.append([i])
                self.PYnunique += 1

            else:
                self.equiv[i_equiv].append(i)
                self.nequiv[i_equiv] += 1
                self.PYatom_to_unique[i] = i_equiv

        # The first atom in the equiv list is considered the primary
        #   unique atom. Just to make things look pretty, make the
        #   atom with the most zeros in its x, y, z coordinate the
        #   unique atom. Nothing else should rely on this being done.
        ztol = 1.0e-5
        for i in range(self.PYnunique):
            maxzero = 0
            jmaxzero = 0
            for j in range(self.nequiv[i]):
                nzero = 0
                for k in range(3):
                    tmp = self.equiv[i][j]
                    if abs(self.xyz(tmp, k)) < ztol:
                        nzero += 1
                if nzero > maxzero:
                    maxzero = nzero
                    jmaxzero = j
            tmp = self.equiv[i][jmaxzero]
            self.equiv[i][jmaxzero] = self.equiv[i][0]
            self.equiv[i][0] = tmp
        #print 'nunique', self.PYnunique
        #print 'nequiv', self.nequiv
        #print 'atom_to_unique', self.PYatom_to_unique
        #print 'equiv', self.equiv

    def sym_label(self):
        """Returns the symmetry label"""
        if self.pg is None:
            self.set_point_group(self.find_point_group())
        return self.pg.symbol()

    def irrep_labels(self):
        """Returns the irrep labels"""
        if self.pg is None:
            self.set_point_group(self.find_point_group())
        return [self.pg.char_table().gamma(i).symbol_ns() for i in range(self.pg.char_table().nirrep())]

    def symmetry_from_input(self):
        """Returns the symmetry specified in the input.

        >>> print H2OH2O.symmetry_from_input()
        C1

        """
        return self.PYsymmetry_from_input

    def symmetrize(self):
        """Force the molecule to have the symmetry specified in pg.
        This is to handle noise coming in from optking.

        """
        #raise FeatureNotImplemented('Molecule::symmetrize')  # FINAL SYMM
        temp = zero(self.natom(), 3)
        ct = self.point_group().char_table()

        # Obtain atom mapping of atom * symm op to atom
        atom_map = compute_atom_map(self)

        # Symmetrize the molecule to remove any noise
        for at in range(self.natom()):
            for g in range(ct.order()):
                Gatom = atom_map[at][g]
                so = ct.symm_operation(g)

                # Full so must be used if molecule is not in standard orientation
                temp[at][0] += so[0][0] * self.x(Gatom) / ct.order()
                temp[at][0] += so[0][1] * self.y(Gatom) / ct.order()
                temp[at][0] += so[0][2] * self.z(Gatom) / ct.order()
                temp[at][1] += so[1][0] * self.x(Gatom) / ct.order()
                temp[at][1] += so[1][1] * self.y(Gatom) / ct.order()
                temp[at][1] += so[1][2] * self.z(Gatom) / ct.order()
                temp[at][2] += so[2][0] * self.x(Gatom) / ct.order()
                temp[at][2] += so[2][1] * self.y(Gatom) / ct.order()
                temp[at][2] += so[2][2] * self.z(Gatom) / ct.order()

        # Set the geometry to ensure z-matrix variables get updated
        self.set_geometry(temp)

    def schoenflies_symbol(self):
        """Returns the Schoenflies symbol"""
        return self.point_group().symbol()

    def valid_atom_map(self, tol=0.01):
        """Check if current geometry fits current point group

        """
        np = [0.0, 0.0, 0.0]
        ct = self.point_group().char_table()

        # loop over all centers
        for at in range(self.natom()):
            ac = self.xyz(at)

            # For each operation in the pointgroup, transform the coordinates of
            #   center "at" and see which atom it maps into
            for g in range(ct.order()):
                so = ct.symm_operation(g)

                for ii in range(3):
                    np[ii] = 0
                    for jj in range(3):
                        np[ii] += so[ii][jj] * ac[jj]

                if self.atom_at_position(np, tol) < 0:
                    return False
        return True

    def full_point_group_with_n(self):
        """Return point group name such as Cnv or Sn."""
        return self.full_pg

    def full_pg_n(self):
        """Return n in Cnv, etc.; If there is no n (e.g. Td)
        it's the highest-order rotation axis.

        """
        return self.PYfull_pg_n

    def get_full_point_group(self):
        """Return point group name such as C3v or S8.
        (method name in libmints is full_point_group)

        """
        pg_with_n = self.full_pg
        if pg_with_n in ['D_inf_h', 'C_inf_v', 'C1', 'Cs', 'Ci', 'Td', 'Oh', 'Ih']:
            return pg_with_n  # These don't need changes - have no 'n'.
        else:
            return pg_with_n.replace('n', str(self.PYfull_pg_n), 1)

    # <<< Methods for Uniqueness >>> (assume molecular point group has been determined)

    def nunique(self):
        """Return the number of unique atoms."""
        return self.PYnunique

    def unique(self, iuniq):
        """Returns the overall number of the iuniq'th unique atom."""
        return self.equiv[iuniq][0]

    def nequivalent(self, iuniq):
        """Returns the number of atoms equivalent to iuniq."""
        return self.nequiv[iuniq]

    def equivalent(self, iuniq, j):
        """Returns the j'th atom equivalent to iuniq."""
        return self.equiv[iuniq][j]

    def atom_to_unique(self, iatom):
        """Converts an atom number to the number of its generating unique atom.
        The return value is in [0, nunique).

        """
        return self.PYatom_to_unique[iatom]

    def atom_to_unique_offset(self, iatom):
        """Converts an atom number to the offset of this atom
        in the list of generated atoms. The unique atom itself is allowed offset 0.

        """
        iuniq = self.PYatom_to_unique[iatom]
        nequiv = self.nequiv[iuniq]
        for i in range(nequiv):
            if self.equiv[iuniq][i] == iatom:
                return i
        raise ValidationError("Molecule::atom_to_unique_offset: I should've found the atom requested...but didn't.")
        return -1

    def max_nequivalent(self):
        """Returns the maximum number of equivalent atoms."""
        mmax = 0
        for i in range(self.nunique()):
            if mmax < self.nequivalent(i):
                mmax = self.nequivalent(i)
        return mmax


def atom_present_in_geom(geom, b, tol=DEFAULT_SYM_TOL):
    """Function used by set_full_point_group() to scan a given geometry
    and determine if an atom is present at a given location.

    """
    for i in range(len(geom)):
        a = [geom[i][0], geom[i][1], geom[i][2]]
        if distance(b, a) < tol:
            return True
    return False


def matrix_3d_rotation_Cn(coord, axis, reflect, tol=DEFAULT_SYM_TOL, max_Cn_to_check=-1):
    """Find maximum n in Cn around given axis, i.e., the highest-order rotation axis.
    @param coord Matrix    : points to rotate - column dim is 3
    @param axis  Vector3   : axis around which to rotate, does not need to be normalized
    @param bool  reflect   : if true, really look for Sn not Cn
    @returns n

    """
    # Check all atoms. In future, make more intelligent.
    max_possible = len(coord) if max_Cn_to_check == -1 else max_Cn_to_check

    Cn = 1  # C1 is there for sure
    for n in range(2, max_possible + 1):
        rotated_mat = matrix_3d_rotation(coord, axis, 2 * math.pi / n, reflect)
        if equal_but_for_row_order(coord, rotated_mat, tol):
            Cn = n
    return Cn


def matrix_3d_rotation(mat, axis, phi, Sn):
    """For a matrix of 3D vectors (ncol==3), rotate a set of points around an
    arbitrary axis.  Vectors are the rows of the matrix.
    @param  axis  Vector3   : axis around which to rotate (need not be normalized)
    @param  phi   double    : magnitude of rotation in rad
    @param  Sn    bool      : if true, then also reflect in plane through origin and perpendicular to rotation
    @returns SharedMatrix with rotated points (rows)

    """
    if len(mat[0]) != 3 or len(axis) != 3:
        raise ValidationError("matrix_3d_rotation: Can only rotate matrix with 3d vectors")

    # Normalize rotation vector
    [wx, wy, wz] = normalize(axis)
    cp = 1.0 - math.cos(phi)

    R = zero(3, 3)
    R[0][0] = wx * wx * cp + math.cos(phi)
    R[0][1] = wx * wy * cp + math.sin(phi) * wz * -1
    R[0][2] = wx * wz * cp + math.sin(phi) * wy
    R[1][0] = wx * wy * cp + math.sin(phi) * wz
    R[1][1] = wy * wy * cp + math.cos(phi)
    R[1][2] = wy * wz * cp + math.sin(phi) * wx * -1
    R[2][0] = wx * wz * cp + math.sin(phi) * wy * -1
    R[2][1] = wy * wz * cp + math.sin(phi) * wx
    R[2][2] = wz * wz * cp + math.cos(phi)

    # R * coord^t = R_coord^t or coord * R^t = R_coord
    #Matrix rotated_coord(nrow(),3);
    #rotated_coord.gemm(false, true, 1.0, *this, R, 0.0);
    rotated_coord = mult(mat, transpose(R))
#    print 'after C'
#    show(rotated_coord)

    if Sn:  # delta_ij - 2 a_i a_j / ||a||^2
        R = identity(3)
        #R = zero(3, 3)
        R[0][0] -= 2 * wx * wx
        R[1][1] -= 2 * wy * wy
        R[2][2] -= 2 * wz * wz
        #R[0][0] = 1 - 2 * wx * wx
        #R[1][1] = 1 - 2 * wy * wy
        #R[2][2] = 1 - 2 * wz * wz
        R[1][0] = 2 * wx * wy
        R[2][0] = 2 * wx * wz
        R[2][1] = 2 * wy * wz
        R[0][1] = 2 * wx * wy
        R[0][2] = 2 * wx * wz
        R[1][2] = 2 * wy * wz
        rotated_coord = mult(rotated_coord, transpose(R))
        #tmp = mult(rotated_coord, transpose(R))
        #Matrix tmp(nrow(),3);
        #tmp.gemm(false, true, 1.0, rotated_coord, R, 0.0);
        #rotated_coord.copy(tmp);
        #rotated_coord = [row[:] for row in tmp]

    #SharedMatrix to_return = rotated_coord.clone();
    #return to_return
    return rotated_coord


def equal_but_for_row_order(mat, rhs, tol=DEFAULT_SYM_TOL):
    """Checks matrix equality, but allows rows to be in a different order.
    @param rhs Matrix to compare to.
    @returns true if equal, otherwise false.

    """
    for m in range(len(mat)):
        for m_rhs in range(len(mat)):

            for n in range(len(mat[m])):
                if abs(mat[m][n] - rhs[m_rhs][n]) > tol:
                    break  # from n
            else:
                # whole row matched, goto next m row
                break  # from m_rhs
        else:
            # no matching row was found
            return False
    else:
        return True


def compute_atom_map(mol):
    """Computes atom mappings during symmetry operations. Useful in
    generating SO information and Cartesian displacement SALCs.
    param mol Molecule to form mapping matrix from.
    returns Integer matrix of dimension natoms X nirreps.

    """
    # create the character table for the point group
    ct = mol.point_group().char_table()

    natom = mol.natom()
    ng = ct.order()
    atom_map = [0] * natom
    for i in range(natom):
        atom_map[i] = [0] * ng

    np = [0.0, 0.0, 0.0]
    so = SymmetryOperation()

    # loop over all centers
    for i in range(natom):
        ac = mol.xyz(i)

        # then for each symop in the pointgroup, transform the coordinates of
        #   center "i" and see which atom it maps into
        for g in range(ng):
            so = ct.symm_operation(g)

            for ii in range(3):
                np[ii] = 0
                for jj in range(3):
                    np[ii] += so[ii][jj] * ac[jj]

            atom_map[i][g] = mol.atom_at_position(np, 0.05)
            if atom_map[i][g] < 0:
                print("""  Molecule:\n""")
                mol.print_out()
                print("""  attempted to find atom at\n""")
                print("""    %lf %lf %lf\n""" % (np[0], np[1], np[2]))
                raise ValidationError("ERROR: Symmetry operation %d did not map atom %d to another atom:\n" % (g, i + 1))

    return atom_map




# TODO outfile
# ignored =, +, 0, += assignment operators
# no pubchem
# TODO rename save_string_for_psi4
# TODO add no_com no_reorint in save string for psi4
