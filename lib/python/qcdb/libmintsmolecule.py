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

import os
import re
import math
import copy
from periodictable import *
from physconst import *
from vecutil import *
from exceptions import *
from libmintscoordentry import *

LINEAR_A_TOL = 1.0E-2  # When sin(a) is below this, we consider the angle to be linear
DEFAULT_SYM_TOL = 1.0E-8
FULL_PG_TOL = 1.0e-8
ZERO = 1.0E-14
NOISY_ZERO = 1.0E-8


class LibmintsMolecule(object):
    """Class to store the elements, coordinates, fragmentation pattern,
    charge, multiplicity of a molecule. Largely replicates psi4's libmints
    Molecule class, developed by Justin M. Turney with incremental
    improvements by other psi4 developers. Major differences from the C++
    class are: no basisset handling, no symmetry, no pubchem. This class
    translated so that databases can function independently of psi4.

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
        # Whether the charge was given by the user  UNUSED
        self.PYcharge_specified = False
        # The multiplicity (defined as 2Ms + 1)
        self.PYmultiplicity = 1
        # Whether the multiplicity was specified by the user  UNUSED
        self.PYmultiplicity_specified = False
        # The units used to define the geometry
        self.PYunits = 'Angstrom'
        # The conversion factor to take input units to Bohr
        self.input_units_to_au = 1.0 / psi_bohr2angstroms
        # Whether this molecule has at least one zmatrix entry
        self.zmat = False

        # <<< Coordinates >>>

        # Atom info vector (no knowledge of dummy atoms)
        self.atoms = []
        # Atom info vector (includes dummy atoms)
        self.full_atoms = []
        # A list of all variables known, whether they have been set or not.
        self.all_variables = []
        # A listing of the variables used to define the geometries
        self.geometry_variables = {}

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

        # Point group to use with this molecule UNUSED
        self.pg = None
        # Full point group UNUSED
        self.full_pg = 'PG_C1'
        # n of the highest rotational axis Cn UNUSED
        self.full_pg_n = 1
        # Symmetry string from geometry specification
        self.PYsymmetry_from_input = 'C1'
        # Number of unique atoms
        self.PYnunique = 0
        # Number of equivalent atoms per unique atom (length nunique)
        self.nequiv = 0
        # Equivalent atom mapping array
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
            raise ValidationError("""Argument to Molecule::set_units must be 'Angstrom' or 'Bohr'.""")

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
        if posn >= 0:
            return temp[posn]
        else:
            return temp

    def mass(self, atom):
        """Returns mass of atom (0-indexed)

        >>> print H2OH2O.mass(4)
        1.00782503207

        """
        if self.atoms[atom].mass() != 0.0:
            return self.atoms[atom].mass()

        if math.fabs(self.atoms[atom].Z() - int(self.atoms[atom].Z())) > 0.0:
            print "WARNING: Obtaining masses from atom with fractional charge...may be incorrect!!!\n"
            # TODO outfile
        return an2masses[int(self.atoms[atom].Z())]

    def symbol(self, atom):
        """Returns the cleaned up label of the atom (C2 => C, H4 = H) (0-indexed)

        >>> print H2OH2O.symbol(4)
        H

        """
        return self.atoms[atom].symbol()

    def label(self, atom):
        """Returns the original label of the atom (0-indexed) as given in the input file (C2, H4).

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
        """Returns the cleaned up label of the atom (C2 => C, H4 = H) (includes dummies)

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
        print 'it lives', 'activate all'
        for fr in range(self.nfragments()):
            print 'reviving', fr
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
        print 'doomed', ghosts
        for fr in ghosts:
            print 'killing', fr - 1
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
        return self.extract_fragments(reals, ghosts=[])

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
            raise ValidationError('The sum of real- and ghost-atom subsets is greater than the number of subsets')

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
        """Given a string *geom* of psi4-style geometry specification
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
        atom = re.compile(r'^\s*(@?[A-Z]{1,2})\s*', re.IGNORECASE)
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
                tempSymm = symmetry.match(line).group(1)
                temp2 = re.sub('[23456789]', 'n', tempSymm).upper()
                if temp2 in (item.upper() for item in self.FullPointGroupList):
                    self.PYsymmetry_from_input = tempSymm

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

            elif atom.match(line):
                glines.append(line)

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
            elif atom.match(line):
                entries = re.split(r'\s+|\s*,\s*', line.strip())
                atomLabel = entries[0]

                # handle ghost atoms
                ghostAtom = False
                if ghost.match(atomLabel):
                    # We don't know whether the @C or Gh(C) notation matched.  Do a quick check.
                    atomLabel = ghost.match(atomLabel).group(2) if not ghost.match(atomLabel).group(1) \
                        else ghost.match(atomLabel).group(1)
                    ghostAtom = True

                # Save the actual atom symbol (H1 => H)
                atomSym = re.split('(\d+)', atomLabel)[0].upper()

                # Check that the atom symbol is valid
                if not atomSym in el2z:
                    raise ValidationError('Illegal atom symbol in geometry specification: %s' % (atomSym))

                zVal = el2z[atomSym]
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
                        el2masses[atomSym], atomSym, atomLabel, \
                        xval, yval, zval))

                # handle first line of Zmat
                elif len(entries) == 1:
                    zmatrix = True
                    tempfrag.append(iatom)
                    self.full_atoms.append(ZMatrixEntry(iatom, zVal, charge, \
                        el2masses[atomSym], atomSym, atomLabel))

                # handle second line of Zmat
                elif len(entries) == 3:
                    zmatrix = True
                    tempfrag.append(iatom)

                    rTo = self.get_anchor_atom(entries[1], line)
                    if rTo >= iatom:
                        raise ValidationError("Error on geometry input line %s. Atom %s has not been defined yet.\n" % (line, entries[1]))
                    rval = self.get_coord_value(entries[2])

                    if self.full_atoms[rTo].symbol() == 'X':
                        rval.set_fixed(True)

                    self.full_atoms.append(ZMatrixEntry(iatom, zVal, charge, \
                        el2masses[atomSym], atomSym, atomLabel, \
                        self.full_atoms[rTo], rval))

                # handle third line of Zmat
                elif len(entries) == 5:
                    zmatrix = True
                    tempfrag.append(iatom)

                    rTo = self.get_anchor_atom(entries[1], line)
                    if rTo >= iatom:
                        raise ValidationError("Error on geometry input line %s. Atom %s has not been defined yet.\n" % (line, entries[1]))
                    aTo = self.get_anchor_atom(entries[3], line)
                    if aTo >= iatom:
                        raise ValidationError("Error on geometry input line %s. Atom %s has not been defined yet.\n" % (line, entries[3]))
                    if aTo == rTo:
                        raise ValidationError("Atom used multiple times on line %s." % (line))
                    rval = self.get_coord_value(entries[2])
                    aval = self.get_coord_value(entries[4])

                    if self.full_atoms[rTo].symbol() == 'X':
                        rval.set_fixed(True)
                    if self.full_atoms[aTo].symbol() == 'X':
                        aval.set_fixed(True)

                    self.full_atoms.append(ZMatrixEntry(iatom, zVal, charge, \
                        el2masses[atomSym], atomSym, atomLabel, \
                        self.full_atoms[rTo], rval, \
                        self.full_atoms[aTo], aval))

                # handle fourth line of Zmat
                elif len(entries) == 7:
                    zmatrix = True
                    tempfrag.append(iatom)

                    rTo = self.get_anchor_atom(entries[1], line)
                    if rTo >= iatom:
                        raise ValidationError("Error on geometry input line %s. Atom %s has not been defined yet.\n" % (line, entries[1]))
                    aTo = self.get_anchor_atom(entries[3], line)
                    if aTo >= iatom:
                        raise ValidationError("Error on geometry input line %s. Atom %s has not been defined yet.\n" % (line, entries[3]))
                    dTo = self.get_anchor_atom(entries[5], line)
                    if dTo >= iatom:
                        raise ValidationError("Error on geometry input line %s. Atom %s has not been defined yet.\n" % (line, entries[5]))
                    if aTo == rTo or rTo == dTo or aTo == dTo:  # for you star wars fans
                        raise ValidationError("Atom used multiple times on line %s" % (line))

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
                        el2masses[atomSym], atomSym, atomLabel, \
                        self.full_atoms[rTo], rval, \
                        self.full_atoms[aTo], aval, \
                        self.full_atoms[dTo], dval))

                else:
                    raise ValidationError('Illegal geometry specification line : %s. \
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
                        raise ValidationError('Illegal atom symbol in geometry specification: %s' % (atomSym))

                    # Add it to the molecule.
                    instance.add_atom(el2z[fileAtom], fileX, fileY, fileZ, fileAtom, el2masses[fileAtom])

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
#            if self.full_pg:  TODO symmetry
#                text += """    Full point group: %s\n\n""" % (self.full_point_group())  TODO symmetry
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
            text += "\n"
        else:
            text += "  No atoms in this molecule.\n"
        print text
        # TODO outfile

    def print_out_in_bohr(self):
        """Print the molecule in Bohr. Same as :py:func:`print_out` only in Bohr.
        (method name in libmints is print_in_bohr)

        """
        text = ""
        if self.natom():
            if self.pg:
                text += """    Molecular point group: %s\n""" % (self.pg.symbol())
#            if self.full_pg:  TODO symmetry
#                text += """    Full point group: %s\n\n""" % (self.full_point_group())  TODO symmetry
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
        print text
        # TODO outfile

    def print_out_in_angstrom(self):
        """Print the molecule in Angstroms. Same as :py:func:`print_out` only always in Angstroms.
        (method name in libmints is print_in_angstrom)

        """
        text = ""
        if self.natom():
            if self.pg:
                text += """    Molecular point group: %s\n""" % (self.pg.symbol())
#            if self.full_pg:  TODO symmetry
#                text += """    Full point group: %s\n\n""" % (self.full_point_group())  TODO symmetry
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
        print text
        # TODO outfile

    def print_full(self):
        """Print full atom list. Same as :py:func:`print_out` only displays dummy atoms.

        """
        text = ""
        if self.natom():
            if self.pg:
                text += """    Molecular point group: %s\n""" % (self.pg.symbol())
#            if self.full_pg:  TODO symmetry
#                text += """    Full point group: %s\n\n""" % (self.full_point_group())  TODO symmetry
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
        print text
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

        print text
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
        print text

    def create_psi4_string_from_molecule(self):
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
                        text += "    %-8s" % (self.fsymbol(at))
                    else:
                        text += "    %-8s" % ("Gh(" + self.fsymbol(at) + ")")
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
            raise ValidationError('ERROR: atom_at_position() requires as argument a vector of length 3\n')

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
        return True if vstr.upper() in self.all_variables else False

    def get_variable(self, vstr):
        """Checks to see if the variable str is in the list, sets it to
        val and returns true if it is, and returns false if not.

        """
        vstr = vstr.upper()
        try:
            return self.geometry_variables[vstr]
        except KeyError:
            raise ValidationError('ERROR: Geometry variable %s not known.\n' % (vstr))

    def set_variable(self, vstr, val):
        """Assigns the value val to the variable labelled string in the
        list of geometry variables. Also calls update_geometry()

        """
        self.__dict__['lock_frame'] = False
        self.__dict__['geometry_variables'][vstr.upper()] = val
        print "Setting geometry variable %s to %f" % (vstr.upper(), val)
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
        try:
            if name.upper() in self.__dict__['all_variables']:
                self.set_variable(name, value)
            else:
                self.__dict__[name] = value
        except KeyError:
            self.__dict__[name] = value

    def __getattr__(self, name):
        """Function to overload accessing attribute contents to allow
        retrivial geometry variable values as if member data.

        """
        if not name in self.__dict__:
            if object.__getattribute__(self, 'is_variable')(name):
                return object.__getattribute__(self, 'get_variable')(name)
            else:
                raise AttributeError
        else:
            return self.__dict__[name]

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
            raise ValidationError("Illegal value %s in atom specification on line %s.\n" % (vstr, line))

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
        print text
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

        #print "beginning update_geometry:"
        #self.print_full()
        if self.PYreinterpret_coordentries:
            self.reinterpret_coordentries()
        #print "after reinterpret_coordentries:"
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
# TODO        self.set_point_group(self.find_point_group())
# TODO        self.set_full_point_group()

        # Disabling symmetrize for now if orientation is fixed, as it is not
        #   correct.  We may want to fix this in the future, but in some cases of
        #   finite-differences the set geometry is not totally symmetric anyway.
        # Symmetrize the molecule to remove any noise
# TODO        self.symmetrize()
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

    def set_basis_all_atoms(self, name, type="BASIS"):
        """ **NYI** Assigns basis *name* to all atoms."""
        raise FeatureNotImplemented('Molecule::set_basis_all_atoms')  # FINAL

    def set_basis_by_symbol(self, symbol, name, type="BASIS"):
        """ **NYI** Assigns basis *name* to all *symbol* atoms."""
        raise FeatureNotImplemented('Molecule::set_basis_by_symbol')  # FINAL

    def set_basis_by_number(self, number, name, type="BASIS"):
        """ **NYI** Assigns basis *name* to atom number *number* (1-indexed, includes dummies)."""
        raise FeatureNotImplemented('Molecule::set_basis_by_number')  # FINAL

    def set_basis_by_label(self, label, name, type="BASIS"):
        """ **NYI** Assigns basis *name* to all atoms with *label*."""
        raise FeatureNotImplemented('Molecule::set_basis_by_label')  # FINAL

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
                    raise ValidationError("Invalid atomic number")
            return nfzc

        else:
            raise ValidationError("Frozen core '%s' is not supported, options are {true, false}." % (depth))

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
        print text
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
        elif rot_const[0] == 0.0:             # A == 0, B == C
            rotor_type = 'RT_LINEAR'
        elif degen == 2:                      # A == B == C
            rotor_type = 'RT_SPHERICAL_TOP'
        elif degen == 1:                      # A > B == C
            rotor_type = 'RT_SYMMETRIC_TOP'   # A == B > C
        else:
            rotor_type = 'RT_ASYMMETRIC_TOP'  # A != B != C
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
        -2 3 water_dimer
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
        text = "%d\n" % (N)
        text += '%d %d %s\n' % (self.molecular_charge(), self.multiplicity(), self.name())

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
        """ **NYI** Whether molecule satisfies the vector symmetry operation *op* """
        raise FeatureNotImplemented('Molecule::has_symmetry_element')  # FINAL SYMM

    def point_group(self):
        """ **NYI** Returns the point group (object) if set"""
        raise FeatureNotImplemented('Molecule::point_group')  # FINAL SYMM

    def set_point_group(self, pg):
        """ **NYI** Set the point group to object *pg* """
        raise FeatureNotImplemented('Molecule::set_point_group')  # FINAL SYMM

    def set_full_point_group(self, tol=FULL_PG_TOL):
        """ **NYI** Determine and set FULL point group"""
        raise FeatureNotImplemented('Molecule::set_full_point_group')  # FINAL SYMM

    def has_inversion(self, origin, tol=DEFAULT_SYM_TOL):
        """Does the molecule have an inversion center at origin"""
        for i in range(self.natom()):
            inverted = sub(origin, sub(self.xyz(i), origin))
            atom = self.atom_at_position(inverted, tol)
            if atom < 0 or not self.atoms[atom].is_equivalent_to(self.atoms[i]):
                return False
        return True

    def is_plane(self, origin, uperp, tol=DEFAULT_SYM_TOL):
        """Is a plane?"""
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
        """Is *axis* an axis of order *order* with respect to *origin*?"""
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
#        if xlikeness > ylikeness and xlikeness > zlikeness:
            like = 'XAxis'
            if dot(axis, worldxaxis) < 0:
                axis = scale(axis, -1.0)
        elif (ylikeness - zlikeness) > 1.0E-12:
#        elif ylikeness > zlikeness:
            like = 'YAxis'
            if dot(axis, worldyaxis) < 0:
                axis = scale(axis, -1.0)
        else:
            like = 'ZAxis'
            if dot(axis, worldzaxis) < 0:
                axis = scale(axis, -1.0)
        return like, axis

    def find_point_group(self, tol=DEFAULT_SYM_TOL):
        """ **NYI** Find computational molecular point group,
        user can override this with the "symmetry" keyword

        """
        raise FeatureNotImplemented('Molecule::find_point_group')  # FINAL SYMM

    def reset_point_group(self, pgname):
        """ **NYI** Override symmetry from outside the molecule string"""
        raise FeatureNotImplemented('Molecule::reset_point_group')  # FINAL SYMM

    def find_highest_point_group(self, tol=DEFAULT_SYM_TOL):
        """ **NYI** Find highest molecular point group"""
        raise FeatureNotImplemented('Molecule::find_highest_point_group')  # FINAL SYMM

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
        """ **NYI** Release symmetry information"""
        raise FeatureNotImplemented('Molecule::release_symmetry_information')  # FINAL SYMM

    def form_symmetry_information(self, tol=DEFAULT_SYM_TOL):
        """ **NYI** Initialize molecular specific symmetry information.
        Uses the point group object obtain by calling point_group()

        """
        raise FeatureNotImplemented('Molecule::form_symmetry_information')  # FINAL SYMM

    def sym_label(self):
        """ **NYI** Returns the symmetry label"""
        raise FeatureNotImplemented('Molecule::sym_label')  # FINAL SYMM

    def irrep_labels(self):
        """ **NYI** Returns the irrep labels"""
        raise FeatureNotImplemented('Molecule::irrep_labels')  # FINAL SYMM

    def symmetry_from_input(self):
        """Returns the symmetry specified in the input.

        >>> print H2OH2O.symmetry_from_input()
        C1

        """
        return self.PYsymmetry_from_input

    def symmetrize(self):
        """ **NYI** Force the molecule to have the symmetry specified in pg.
        This is to handle noise coming in from optking.

        """
        raise FeatureNotImplemented('Molecule::symmetrize')  # FINAL SYMM

    def schoenflies_symbol(self):
        """ **NYI** Returns the Schoenflies symbol"""
        raise FeatureNotImplemented('Molecule::schoenflies_symbol')  # FINAL SYMM

    def valid_atom_map(self, tol=0.01):
        """ **NYI** Check if current geometry fits current point group"""
        raise FeatureNotImplemented('Molecule::valid_atom_map')  # FINAL SYMM

    def full_point_group_with_n(self):
        """ **NYI** Return point group name such as Cnv or Sn."""
        #return FullPointGroupList[self.full_pg]
        raise FeatureNotImplemented('Molecule::full_point_group_n')  # FINAL SYMM

    def full_pg_n(self):
        """ **NYI** Return n in Cnv, etc.; If there is no n (e.g. Td)
        it's the highest-order rotation axis.

        """
        #return self.full_pg_n
        raise FeatureNotImplemented('Molecule::full_pg_n')  # FINAL SYMM

    def get_full_point_group(self):
        """ **NYI** Return point group name such as C3v or S8.
        (method name in libmints is full_point_group)

        """
        raise FeatureNotImplemented('Molecule::get_full_point_group')  # FINAL SYMM

    # <<< Methods for Uniqueness >>> (assume molecular point group has been determined)

    def nunique(self):
        """ **NYI** Return the number of unique atoms."""
        #w#return PYnunique
        raise FeatureNotImplemented('Molecule::nunique')  # FINAL SYMM

    def unique(self, iuniq):
        """ **NYI** Returns the overall number of the iuniq'th unique atom."""
        #w#return self.equiv[iuniq][0]
        raise FeatureNotImplemented('Molecule::unique')  # FINAL SYMM

    def nequivalent(self, iuniq):
        """ **NYI** Returns the number of atoms equivalent to iuniq."""
        #w#return self.nequiv[iuniq]
        raise FeatureNotImplemented('Molecule::nequivalent')  # FINAL SYMM

    def equivalent(self, iuniq, j):
        """ **NYI** Returns the j'th atom equivalent to iuniq."""
        #w#return self.equiv[iuniq][j]
        raise FeatureNotImplemented('Molecule::equivalent')  # FINAL SYMM

    def atom_to_unique(self, iatom):
        """ **NYI** Converts an atom number to the number of its generating unique atom.
        The return value is in [0, nunique).

        """
        #w#return PYatom_to_unique[iatom]
        raise FeatureNotImplemented('Molecule::atom_to_unique')  # FINAL SYMM

    def atom_to_unique_offset(self, iatom):
        """ **NYI** Converts an atom number to the offset of this atom
        in the list of generated atoms. The unique atom itself is allowed offset 0.

        """
        raise FeatureNotImplemented('Molecule::atom_to_unique_offset')  # FINAL SYMM

    def max_nequivalent(self):
        """  **NYI** Returns the maximum number of equivalent atoms."""
        raise FeatureNotImplemented('Molecule::max_nequivalent')  # FINAL SYMM



# TODO outfile
# ignored =, +, 0, += assignment operators
# no pubchem
# no symmetry
# TODO rename save_string_for_psi4
# TODO add no_com no_reorint in save string for psi4
