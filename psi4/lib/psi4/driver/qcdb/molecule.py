#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2024 The Psi4 Developers.
#
# The copyrights for code used from other parties are included in
# the corresponding files.
#
# This file is part of Psi4.
#
# Psi4 is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, version 3.
#
# Psi4 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License along
# with Psi4; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
# @END LICENSE
#

import os
import hashlib
import collections
from typing import Dict, List, Optional, Tuple, Union

import numpy as np

import qcelemental as qcel

import psi4
from .constants import constants
from .util import parse_dertype
from .libmintsmolecule import *
from .testing import compare_values, compare_integers, compare_molrecs
from .bfs import BFS

qcdbmol = "psi4.driver.qcdb.molecule.Molecule"

class Molecule(LibmintsMolecule):
    """Class to store the elements, coordinates, fragmentation pattern,
    charge, multiplicity of a molecule. Largely replicates psi4's libmints
    Molecule class, developed by Justin M. Turney and Andy M. Simmonett
    with incremental improvements by other psi4 developers. Major

    This class extends `qcdb.LibmintsMolecule` and occasionally
    `psi4.core.Molecule` itself.

    """

    def __init__(self,
                 molinit=None,
                 dtype=None,
                 geom=None,
                 elea=None,
                 elez=None,
                 elem=None,
                 mass=None,
                 real=None,
                 elbl=None,
                 name=None,
                 units='Angstrom',
                 input_units_to_au=None,
                 fix_com=None,
                 fix_orientation=None,
                 fix_symmetry=None,
                 fragment_separators=None,
                 fragment_charges=None,
                 fragment_multiplicities=None,
                 molecular_charge=None,
                 molecular_multiplicity=None,
                 comment=None,
                 provenance=None,
                 connectivity=None,
                 enable_qm=True,
                 enable_efp=True,
                 missing_enabled_return_qm='none',
                 missing_enabled_return_efp='none',
                 missing_enabled_return='error',
                 tooclose=0.1,
                 zero_ghost_fragments=False,
                 nonphysical=False,
                 mtol=1.e-3,
                 verbose=1):
        """Initialize Molecule object from LibmintsMolecule"""
        super(Molecule, self).__init__()

        if molinit is not None or geom is not None:
            if isinstance(molinit, dict):
                molrec = molinit

            elif isinstance(molinit, str):
                compound_molrec = qcel.molparse.from_string(
                    molstr=molinit,
                    dtype=dtype,
                    name=name,
                    fix_com=fix_com,
                    fix_orientation=fix_orientation,
                    fix_symmetry=fix_symmetry,
                    return_processed=False,
                    enable_qm=enable_qm,
                    enable_efp=enable_efp,
                    missing_enabled_return_qm=missing_enabled_return_qm,
                    missing_enabled_return_efp=missing_enabled_return_efp,
                    verbose=verbose)
                molrec = compound_molrec['qm']

            elif molinit is None and geom is not None:
                molrec = qcel.molparse.from_arrays(
                    geom=geom,
                    elea=elea,
                    elez=elez,
                    elem=elem,
                    mass=mass,
                    real=real,
                    elbl=elbl,
                    name=name,
                    units=units,
                    input_units_to_au=input_units_to_au,
                    fix_com=fix_com,
                    fix_orientation=fix_orientation,
                    fix_symmetry=fix_symmetry,
                    fragment_separators=fragment_separators,
                    fragment_charges=fragment_charges,
                    fragment_multiplicities=fragment_multiplicities,
                    molecular_charge=molecular_charge,
                    molecular_multiplicity=molecular_multiplicity,
                    comment=comment,
                    provenance=provenance,
                    connectivity=connectivity,
                    domain='qm',
                    missing_enabled_return=missing_enabled_return,
                    tooclose=tooclose,
                    zero_ghost_fragments=zero_ghost_fragments,
                    nonphysical=nonphysical,
                    mtol=mtol,
                    verbose=verbose)

            # ok, got the molrec dictionary; now build the thing
            self._internal_from_dict(molrec, verbose=verbose)

        # The comment line
        self.tagline = ""

    def __str__(self):
        text = """  ==> qcdb Molecule %s <==\n\n""" % (self.name())
        text += """   => %s <=\n\n""" % (self.tagline)
        text += self.create_psi4_string_from_molecule()
        return text

    def __setattr__(self, name, value):
        """Function to overload setting attributes to allow geometry
        variable assigment as if member data.

        """
        if 'all_variables' in self.__dict__:
            if name.upper() in self.__dict__['all_variables']:
                self.set_variable(name, value)
        super(Molecule, self).__setattr__(name, value)

    def __getattr__(self, name):
        """Function to overload accessing attribute contents to allow
        retrival of geometry variable values as if member data.

        """
        if 'all_variables' in self.__dict__ and name.upper() in self.__dict__['all_variables']:
            return self.get_variable(name)
        else:
            raise AttributeError

    @classmethod
    def init_with_xyz(cls, xyzfilename, no_com=False, no_reorient=False, contentsNotFilename=False):
        """Pull information from an XYZ file. No fragment info detected.
        Bohr/Angstrom pulled from first line if available.  Charge,
        multiplicity, tagline pulled from second line if available.  Body
        accepts atom symbol or atom charge in first column. Arguments
        *no_com* and *no_reorient* can be used to turn off shift and
        rotation. If *xyzfilename* is a string of the contents of an XYZ
        file, rather than the name of a file, set *contentsNotFilename*
        to ``True``.

        >>> H2O = qcdb.Molecule.init_with_xyz('h2o.xyz')

        """
        raise FeatureDeprecated(
            """qcdb.Molecule.init_with_xyz. Replace with: qcdb.Molecule.from_string(..., dtype='xyz+')""")

    @classmethod
    def init_with_mol2(cls, xyzfilename, no_com=False, no_reorient=False, contentsNotFilename=False):
        """Pull information from a MOl2 file. No fragment info detected.
        Bohr/Angstrom pulled from first line if available.  Charge,
        multiplicity, tagline pulled from second line if available.  Body
        accepts atom symbol or atom charge in first column. Arguments
        *no_com* and *no_reorient* can be used to turn off shift and
        rotation. If *xyzfilename* is a string of the contents of an XYZ
        file, rather than the name of a file, set *contentsNotFilename*
        to ``True``.

        NOTE: chg/mult NYI

        >>> H2O = qcdb.Molecule.init_with_mol2('h2o.mol2')

        """
        instance = cls()
        instance.lock_frame = False
        instance.PYmove_to_com = not no_com
        instance.PYfix_orientation = no_reorient

        if contentsNotFilename:
            text = xyzfilename.splitlines()
        else:
            try:
                infile = open(xyzfilename, 'r')
            except IOError:
                raise ValidationError(
                    """Molecule::init_with_mol2: given filename '%s' does not exist.""" % (xyzfilename))
            if os.stat(xyzfilename).st_size == 0:
                raise ValidationError("""Molecule::init_with_mol2: given filename '%s' is blank.""" % (xyzfilename))
            text = infile.readlines()

        # fixed-width regex ((?=[ ]*-?\d+)[ -\d]{5})
        v2000 = re.compile(r'^((?=[ ]*\d+)[ \d]{3})((?=[ ]*\d+)[ \d]{3})(.*)V2000\s*$')
        vend = re.compile(r'^\s*M\s+END\s*$')
        NUMBER = "((?:[-+]?\\d*\\.\\d+(?:[DdEe][-+]?\\d+)?)|(?:[-+]?\\d+\\.\\d*(?:[DdEe][-+]?\\d+)?))"
        xyzM = re.compile(
            r'^(?:\s*)' + NUMBER + r'(?:\s+)' + NUMBER + r'(?:\s+)' + NUMBER + r'(?:\s+)([A-Z](?:[a-z])?)(?:\s+)(.*)',
            re.IGNORECASE)

        ## now charge and multiplicity
        #   $chargem = 0  ; $multm = 1 ;
        #while (<MOL>) {
        #if (/CHARGE/) { $chargem = <MOL> ; chop($chargem) ;}
        #if (/MULTIPLICITY/) { $multm = <MOL> ; chop($multm) }
        #        } # end while charge and multiplicity

        if not text:
            raise ValidationError("Molecule::init_with_mol2: file blank")
        # Try to match header/footer
        if vend.match(text[-1]):
            pass
        else:
            raise ValidationError("Molecule::init_with_mol2: Malformed file termination\n%s" % (text[-1]))
        sysname = '_'.join(text[0].strip().split())
        comment = text[2].strip()
        if comment:
            instance.tagline = sysname + ' ' + comment
        else:
            instance.tagline = sysname
        #instance.tagline = text[0].strip() + ' ' + text[2].strip()
        fileUnits = 'Angstrom'  # defined for MOL
        #instance.set_molecular_charge(int(xyz2.match(text[1]).group(1)))
        #instance.set_multiplicity(int(xyz2.match(text[1]).group(2)))
        if v2000.match(text[3]):
            fileNatom = int(v2000.match(text[3]).group(1))
            fileNbond = int(v2000.match(text[3]).group(2))
        else:
            raise ValidationError("Molecule::init_with_mol2: Malformed fourth line\n%s" % (text[3]))
        if fileNatom < 1:
            raise ValidationError("Molecule::init_with_mol2: Malformed Natom\n%s" % (str(fileNatom)))

        # Next line begins the useful information.
        for i in range(fileNatom):
            try:
                if xyzM.match(text[4 + i]):

                    fileX = float(xyzM.match(text[4 + i]).group(1))
                    fileY = float(xyzM.match(text[4 + i]).group(2))
                    fileZ = float(xyzM.match(text[4 + i]).group(3))
                    fileAtom = xyzM.match(text[4 + i]).group(4).upper()

                    # Check that the atom symbol is valid
                    z = qcel.periodictable.to_Z(fileAtom)

                    # Add it to the molecule.
                    instance.add_atom(z, fileX, fileY, fileZ, fileAtom, qcel.periodictable.to_mass(fileAtom), z)

                else:
                    raise ValidationError("Molecule::init_with_mol2: Malformed atom information line %d." % (i + 5))

            except IndexError:
                raise ValidationError(
                    "Molecule::init_with_mol2: Expected atom in file at line %d.\n%s" % (i + 5, text[i + 4]))

        # We need to make 1 fragment with all atoms
        instance.fragments.append([0, fileNatom - 1])
        instance.fragment_types.append('Real')
        instance.fragment_charges.append(instance.molecular_charge())
        instance.fragment_multiplicities.append(instance.multiplicity())
        # Set the units properly
        instance.PYunits = fileUnits
        if fileUnits == 'Bohr':
            instance.PYinput_units_to_au = 1.0
        elif fileUnits == 'Angstrom':
            instance.PYinput_units_to_au = 1.0 / constants.bohr2angstroms

        instance.update_geometry()
        return instance

    def save_string_xyz(self, save_ghosts=True, save_natom=False):
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
        factor = 1.0 if self.PYunits == 'Angstrom' else constants.bohr2angstroms

        N = self.natom()
        if not save_ghosts:
            N = 0
            for i in range(self.natom()):
                if self.Z(i):
                    N += 1
        text = ''
        if save_natom:
            text += "%d\n" % (N)
        text += '%d %d %s\n' % (self.molecular_charge(), self.multiplicity(), self.tagline)

        for i in range(self.natom()):
            [x, y, z] = self.atoms[i].compute()
            if save_ghosts or self.Z(i):
                text += '%2s %17.12f %17.12f %17.12f\n' % ((self.symbol(i) if self.Z(i) else "Gh"),
                    x * factor, y * factor, z * factor)
        return text

    def save_xyz(self, filename, save_ghosts=True, save_natom=True):
        """Save an XYZ file.

        >>> H2OH2O.save_xyz('h2o.xyz')

        """
        outfile = open(filename, 'w')
        outfile.write(self.save_string_xyz(save_ghosts, save_natom))
        outfile.close()

    def format_molecule_for_numpy(self, npobj=True):
        """Returns a NumPy array of the non-dummy atoms of the geometry
        in Cartesian coordinates in Angstroms with element encoded as
        atomic number. If *npobj* is False, returns representation of
        NumPy array.

        """
        factor = 1.0 if self.PYunits == 'Angstrom' else constants.bohr2angstroms
        self.update_geometry()

        # TODO fn title is format_mol... but return args not compatible
        geo = []
        for i in range(self.natom()):
            [x, y, z] = self.atoms[i].compute()
            geo.append([self.Z(i), x * factor, y * factor, z * factor])

        nparr = np.array(geo)
        return nparr if npobj else np.array_repr(nparr)

    def format_molecule_for_psi4(self):
        """Returns string of molecule definition block."""
        text = 'molecule mol {\n'
        for line in self.create_psi4_string_from_molecule().splitlines():
            text += '   ' + line + '\n'
        text += '}\n'
        return text

    def format_molecule_for_qchem_old(self, mixedbas=True):
        """Returns geometry section of input file formatted for Q-Chem.
        For ghost atoms, prints **Gh** as elemental symbol, with expectation
        that element identity will be established in mixed basis section.
        For ghost atoms when *mixedbas* is False, prints @ plus element symbol.

        prints whole dimer for unCP mono when called dir (as opposed to passing thru str
        no frag markers
        """
        factor = 1.0 if self.PYunits == 'Angstrom' else constants.bohr2angstroms

        text = ""
        text += '$molecule\n'
        text += '%d %d\n' % (self.molecular_charge(), self.multiplicity())

        for i in range(self.natom()):
            [x, y, z] = self.atoms[i].compute()
            if mixedbas:
                text += '%2s ' % (self.symbol(i) if self.Z(i) else "Gh")
            else:
                text += '%-3s ' % (('' if self.Z(i) else '@') + self.symbol(i))
            text += '%17.12f %17.12f %17.12f\n' % (x * factor, y * factor, z * factor)
        text += '$end\n\n'

        # prepare molecule keywords to be set as c-side keywords
        options = collections.defaultdict(lambda: collections.defaultdict(dict))
        #options['QCHEM'['QCHEM_CHARGE']['value'] = self.molecular_charge()
        #options['QCHEM'['QCHEM_MULTIPLICITY']['value'] = self.multiplicity()
        options['QCHEM']['QCHEM_INPUT_BOHR']['value'] = False
        #options['QCHEM']['QCHEM_COORDINATES']['value'] = 'CARTESIAN'
        #SYM_IGNORE equiv to no_reorient, no_com, symmetry c1

        options['QCHEM']['QCHEM_INPUT_BOHR']['clobber'] = True

        return text, options

    def format_molecule_for_psi4_xyz(self):
        """not much examined

        """
        text = ""
        if self.nallatom():

            factor = 1.0 if self.PYunits == 'Angstrom' else constants.bohr2angstroms
            # append units and any other non-default molecule keywords
            text += "units Angstrom\n"
            #text += "    units %-s\n" % ("Angstrom" if self.units() == 'Angstrom' else "Bohr")
            if not self.PYmove_to_com:
                text += "no_com\n"
            if self.PYfix_orientation:
                text += "no_reorient\n"

            # append atoms and coordentries and fragment separators with charge and multiplicity
            Pfr = 0
            for fr in range(self.nfragments()):
                if self.fragment_types[fr] == 'Absent' and not self.has_zmatrix():
                    continue
                text += "%s%s%d %d\n" % ("" if Pfr == 0 else "--\n", "#" if self.fragment_types[fr] == 'Ghost'
                                         or self.fragment_types[fr] == 'Absent' else "", self.fragment_charges[fr],
                                         self.fragment_multiplicities[fr])
                Pfr += 1
                for at in range(self.fragments[fr][0], self.fragments[fr][1] + 1):
                    if self.fragment_types[fr] == 'Absent' or self.fsymbol(at) == "X":
                        pass
                    else:
                        if self.fZ(at):
                            text += "%-8s" % (self.flabel(at))
                        else:
                            text += "%-8s" % ("Gh(" + self.flabel(at) + ")")
                        [x, y, z] = self.full_atoms[at].compute()
                        text += '%17.12f %17.12f %17.12f\n' % \
                            (x * factor, y * factor, z * factor)
            text += "\n"

        wtext = 'molecule mol {\n'
        for line in text.splitlines():
            wtext += '   ' + line + '\n'
        wtext += '}\n'
        return wtext

    def format_molecule_for_molpro(self):
        """

        """
        factor = 1.0 if self.PYunits == 'Angstrom' else constants.bohr2angstroms
        # TODO keep fix_or?  # Jan 2015 turning off fix_or
        #self.fix_orientation(True)
        #self.PYmove_to_com = False
        self.update_geometry()

        text = ""
        text += 'angstrom\n'
        text += 'geometry={\n'
        dummy = []

        for i in range(self.natom()):
            [x, y, z] = self.atoms[i].compute()
            text += '%-2s %17.12f %17.12f %17.12f\n' % (self.symbol(i),
                x * factor, y * factor, z * factor)
            if not self.Z(i):
                dummy.append(str(i + 1))  # Molpro atom number is 1-indexed

        text += '}\n\n'
        text += 'SET,CHARGE=%d\n' % (self.molecular_charge())
        text += 'SET,SPIN=%d\n' % (self.multiplicity() - 1)  # Molpro wants (mult-1)
        if len(dummy) > 0:
            text += 'dummy,' + ','.join(dummy) + '\n'
        return text

    def format_molecule_for_cfour(self):
        """Function to print Molecule in a form readable by Cfour.

        """
        self.update_geometry()
        factor = 1.0 if self.PYunits == 'Angstrom' else constants.bohr2angstroms
        #factor = 1.0 if self.PYunits == 'Bohr' else 1.0/psi_bohr2angstroms

        text = 'auto-generated by qcdb from molecule %s\n' % (self.tagline)

        # append atoms and coordentries
        for i in range(self.natom()):
            [x, y, z] = self.atoms[i].compute()
            text += '%-2s %17.12f %17.12f %17.12f\n' % ((self.symbol(i) if self.Z(i) else "GH"),
                x * factor, y * factor, z * factor)

        #for fr in range(self.nfragments()):
        #    if self.fragment_types[fr] == 'Absent':
        #        pass
        #    else:
        #        for at in range(self.fragments[fr][0], self.fragments[fr][1] + 1):
        #            [x, y, z] = self.atoms[at].compute()
        #            text += '%-2s %17.12f %17.12f %17.12f\n' % ((self.symbol(at) if self.Z(at) else "GH"), \
        #                x * factor, y * factor, z * factor)
        text += '\n'

        # prepare molecule keywords to be set as c-side keywords
        options = collections.defaultdict(lambda: collections.defaultdict(dict))
        options['CFOUR']['CFOUR_CHARGE']['value'] = self.molecular_charge()
        options['CFOUR']['CFOUR_MULTIPLICITY']['value'] = self.multiplicity()
        options['CFOUR']['CFOUR_UNITS']['value'] = 'ANGSTROM'
        #options['CFOUR']['CFOUR_UNITS']['value'] = 'BOHR'
        options['CFOUR']['CFOUR_COORDINATES']['value'] = 'CARTESIAN'
        #options['CFOUR']['CFOUR_SUBGROUP']['value'] = self.symmetry_from_input().upper()
        #print self.inertia_tensor()
        #print self.inertial_system()

        options['CFOUR']['CFOUR_CHARGE']['clobber'] = True
        options['CFOUR']['CFOUR_MULTIPLICITY']['clobber'] = True
        options['CFOUR']['CFOUR_UNITS']['clobber'] = True
        options['CFOUR']['CFOUR_COORDINATES']['clobber'] = True

        return text, options

    def format_basis_for_cfour(self, puream):
        """Function to print the BASIS=SPECIAL block for Cfour according
        to the active atoms in Molecule. Special short basis names
        are used by Psi4 libmints GENBAS-writer in accordance with
        Cfour constraints.

        """
        text = ''
        cr = 1
        for fr in range(self.nfragments()):
            if self.fragment_types[fr] == 'Absent':
                pass
            else:
                for at in range(self.fragments[fr][0], self.fragments[fr][1] + 1):
                    text += """%s:P4_%d\n""" % (self.symbol(at).upper(), cr)
                    cr += 1
        text += '\n'

        options = collections.defaultdict(lambda: collections.defaultdict(dict))
        options['CFOUR']['CFOUR_BASIS']['value'] = 'SPECIAL'
        options['CFOUR']['CFOUR_SPHERICAL']['value'] = puream

        options['CFOUR']['CFOUR_BASIS']['clobber'] = True
        options['CFOUR']['CFOUR_SPHERICAL']['clobber'] = True

        options['CFOUR']['CFOUR_BASIS']['superclobber'] = True
        options['CFOUR']['CFOUR_SPHERICAL']['superclobber'] = True

        return text, options

    def format_molecule_for_orca(self):
        """
        Format the molecule into an orca xyz format
        """
        options = collections.defaultdict(lambda: collections.defaultdict(dict))
        self.update_geometry()
        factor = 1.0 if self.PYunits == 'Angstrom' else constants.bohr2angstroms

        text = ""
        text += '* xyz {} {}\n'.format(self.molecular_charge(), self.multiplicity())

        n_frags = self.nfragments()
        for fr in range(n_frags):
            if self.fragment_types[fr] == 'Absent':
                pass
            else:
                for at in range(self.fragments[fr][0], self.fragments[fr][1] + 1):
                    if self.fragment_types[fr] == 'Ghost':
                        # TODO: add support for ghost atoms
                        # atom += ':'
                        continue
                    x, y, z = self.atoms[at].compute()
                    atom = self.symbol(at)
                    if n_frags > 1:
                        text += '    {:2s}({:d}) {:> 17.12f} {:> 17.12f} {:> 17.12f}\n'.format(
                                atom, fr + 1, x * factor, y * factor, z * factor)
                    else:
                        text += '    {:2s} {:> 17.12f} {:> 17.12f} {:> 17.12f}\n'.format(
                                atom, x * factor, y * factor, z * factor)
        text += '*'

        return text, options

    def format_molecule_for_qchem(self, mixedbas=True):
        """Returns geometry section of input file formatted for Q-Chem.
        For ghost atoms, prints **Gh** as elemental symbol, with expectation
        that element identity will be established in mixed basis section.
        For ghost atoms when *mixedbas* is False, prints @ plus element symbol.

        candidate modeled after psi4_xyz so that absent fragments observed force xyz

        """
        text = ""
        if self.nallatom():
            factor = 1.0 if self.PYunits == 'Angstrom' else constants.bohr2angstroms
            Pfr = 0
            # any general starting notation here <<<
            text += '$molecule\n'
            text += '%d %d\n' % (self.molecular_charge(), self.multiplicity())
            # >>>
            for fr in range(self.nfragments()):
                if self.fragment_types[fr] == 'Absent' and not self.has_zmatrix():
                    continue
                # any fragment marker here <<<
                if self.nactive_fragments() > 1:
                    # this only distinguishes Real frags so Real/Ghost don't get
                    #   fragmentation. may need to change
                    text += """--\n"""
                # >>>
                # any fragment chgmult here <<<
                if self.nactive_fragments() > 1:
                    text += """{}{} {}\n""".format('!' if self.fragment_types[fr] in ['Ghost', 'Absent'] else '',
                                                   self.fragment_charges[fr], self.fragment_multiplicities[fr])
                # >>>
                Pfr += 1
                for at in range(self.fragments[fr][0], self.fragments[fr][1] + 1):
                    if self.fragment_types[fr] == 'Absent' or self.fsymbol(at) == "X":
                        pass
                    else:
                        if self.fZ(at):
                            # label for real live atom <<<
                            text += """{:>3s} """.format(self.fsymbol(at))
                            # >>>
                        else:
                            # label for ghost atom <<<
                            text += """{:>3s} """.format('Gh' if mixedbas else ('@' + self.fsymbol(at)))
                            # >>>
                        [x, y, z] = self.full_atoms[at].compute()
                        # Cartesian coordinates <<<
                        text += """{:>17.12f} {:>17.12f} {:>17.12f}\n""".format(x * factor, y * factor, z * factor)
                        # >>>
            # any general finishing notation here <<<
            text += '$end\n\n'
            # >>>

        # prepare molecule keywords to be set as c-side keywords
        options = collections.defaultdict(lambda: collections.defaultdict(dict))
        #options['QCHEM'['QCHEM_CHARGE']['value'] = self.molecular_charge()
        #options['QCHEM'['QCHEM_MULTIPLICITY']['value'] = self.multiplicity()
        options['QCHEM']['QCHEM_INPUT_BOHR']['value'] = False
        #options['QCHEM']['QCHEM_COORDINATES']['value'] = 'CARTESIAN'
        if (not self.PYmove_to_com) or self.PYfix_orientation:
            options['QCHEM']['QCHEM_SYM_IGNORE']['value'] = True
            #SYM_IGNORE equiv to no_reorient, no_com, symmetry c1

        options['QCHEM']['QCHEM_INPUT_BOHR']['clobber'] = True
        options['QCHEM']['QCHEM_SYM_IGNORE']['clobber'] = True

        return text, options

    def format_molecule_for_cfour_old(self):
        """Function to print Molecule in a form readable by Cfour. This
        version works as long as zmat is composed entirely of variables,
        not internal values, while cartesian is all internal values,
        no variables. Cutting off this line of development because,
        with getting molecules after passing through libmints Molecule,
        all zmats with dummies (Cfour's favorite kind) have already been
        converted into cartesian. Next step, if this line was pursued
        would be to shift any zmat internal values to external and any
        cartesian external values to internal.

        """

        text = ''
        text += 'auto-generated by qcdb from molecule %s\n' % (self.tagline)

        ## append units and any other non-default molecule keywords
        #text += "    units %-s\n" % ("Angstrom" if self.units() == 'Angstrom' else "Bohr")
        #if not self.PYmove_to_com:
        #    text += "    no_com\n"
        #if self.PYfix_orientation:
        #    text += "    no_reorient\n"

        # append atoms and coordentries and fragment separators with charge and multiplicity
        Pfr = 0
        isZMat = False
        isCart = False
        for fr in range(self.nfragments()):
            if self.fragment_types[fr] == 'Absent' and not self.has_zmatrix():
                continue
#            text += "%s    %s%d %d\n" % (
#                "" if Pfr == 0 else "    --\n",
#                "#" if self.fragment_types[fr] == 'Ghost' or self.fragment_types[fr] == 'Absent' else "",
#                self.fragment_charges[fr], self.fragment_multiplicities[fr])
            Pfr += 1
            for at in range(self.fragments[fr][0], self.fragments[fr][1] + 1):
                if type(self.full_atoms[at]) == ZMatrixEntry:
                    isZMat = True
                elif type(self.full_atoms[at]) == CartesianEntry:
                    isCart = True
                if self.fragment_types[fr] == 'Absent':
                    text += "%s" % ("X")
                elif self.fZ(at) or self.fsymbol(at) == "X":
                    text += "%s" % (self.fsymbol(at))
                else:
                    text += "%s" % ("GH")  # atom info is lost + self.fsymbol(at) + ")")
                text += "%s" % (self.full_atoms[at].print_in_input_format_cfour())
        text += "\n"

        # append any coordinate variables
        if len(self.geometry_variables):
            for vb, val in self.geometry_variables.items():
                text += """%s=%.10f\n""" % (vb, val)
            text += "\n"

        # prepare molecule keywords to be set as c-side keywords
        options = collections.defaultdict(lambda: collections.defaultdict(dict))
        options['CFOUR']['CFOUR_CHARGE']['value'] = self.molecular_charge()
        options['CFOUR']['CFOUR_MULTIPLICITY']['value'] = self.multiplicity()
        options['CFOUR']['CFOUR_UNITS']['value'] = self.units()
        if isZMat and not isCart:
            options['CFOUR']['CFOUR_COORDINATES']['value'] = 'INTERNAL'
        elif isCart and not isZMat:
            options['CFOUR']['CFOUR_COORDINATES']['value'] = 'CARTESIAN'
        else:
            raise ValidationError("""Strange mix of Cartesian and ZMatrixEntries in molecule unsuitable for Cfour.""")

        return text, options

    def format_molecule_for_nwchem(self):
        """

        """
        factor = 1.0 if self.PYunits == 'Angstrom' else constants.bohr2angstroms

        text = ""
        text += '%d %d %s\n' % (self.molecular_charge(), self.multiplicity(), self.tagline)

        for i in range(self.natom()):
            [x, y, z] = self.atoms[i].compute()
            text += '%4s %17.12f %17.12f %17.12f\n' % (("" if self.Z(i) else 'Bq') + self.symbol(i),
                x * factor, y * factor, z * factor)
        return text
        pass

    #    if symm   print M2OUT "nosym\nnoorient\n";
    #    print DIOUT "angstrom\ngeometry={\n";

    def inertia_tensor(self, masswt=True, zero=ZERO):
        """Compute inertia tensor.

        >>> print H2OH2O.inertia_tensor()
        [[8.704574864178731, -8.828375721817082, 0.0], [-8.828375721817082, 280.82861714077666, 0.0], [0.0, 0.0, 281.249500988553]]

        """
        return self.inertia_tensor_partial(range(self.natom()), masswt, zero)

    def inertia_tensor_partial(self, part, masswt=True, zero=ZERO):
        """Compute inertia tensor based on atoms in *part*.

        """
        tensor = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]

        for i in part:
            if masswt:
                # I(alpha, alpha)
                tensor[0][0] += self.mass(i) * (self.y(i) * self.y(i) + self.z(i) * self.z(i))
                tensor[1][1] += self.mass(i) * (self.x(i) * self.x(i) + self.z(i) * self.z(i))
                tensor[2][2] += self.mass(i) * (self.x(i) * self.x(i) + self.y(i) * self.y(i))

                # I(alpha, beta)
                tensor[0][1] -= self.mass(i) * self.x(i) * self.y(i)
                tensor[0][2] -= self.mass(i) * self.x(i) * self.z(i)
                tensor[1][2] -= self.mass(i) * self.y(i) * self.z(i)

            else:
                # I(alpha, alpha)
                tensor[0][0] += self.y(i) * self.y(i) + self.z(i) * self.z(i)
                tensor[1][1] += self.x(i) * self.x(i) + self.z(i) * self.z(i)
                tensor[2][2] += self.x(i) * self.x(i) + self.y(i) * self.y(i)

                # I(alpha, beta)
                tensor[0][1] -= self.x(i) * self.y(i)
                tensor[0][2] -= self.x(i) * self.z(i)
                tensor[1][2] -= self.y(i) * self.z(i)

        # mirror
        tensor[1][0] = tensor[0][1]
        tensor[2][0] = tensor[0][2]
        tensor[2][1] = tensor[1][2]

        # Check the elements for zero and make them a hard zero.
        for i in range(3):
            for j in range(3):
                if math.fabs(tensor[i][j]) < zero:
                    tensor[i][j] = 0.0
        return tensor

    def inertial_system_partial(self, part, masswt=True, zero=ZERO):
        """Solve inertial system based on atoms in *part*"""
        return diagonalize3x3symmat(self.inertia_tensor_partial(part, masswt, zero))

    def inertial_system(self, masswt=True, zero=ZERO):
        """Solve inertial system"""
        return diagonalize3x3symmat(self.inertia_tensor(masswt, zero))

    def print_ring_planes(self, entity1, entity2, entity3=None, entity4=None):
        """(reals only, 1-indexed)

        """
        pass
        # TODO allow handle lines
        text = ""
        summ = []

        #for entity in [entity1, entity2, entity3, entity4]:
        for item in [entity1, entity2]:

            text += """\n  ==> Entity %s <==\n\n""" % (item)

            # convert plain atoms into list and move from 1-indexed to 0-indexed
            entity = []
            try:
                for idx in item:
                    entity.append(idx - 1)
            except TypeError:
                entity = [item - 1]

            if len(entity) == 1:
                dim = 'point'
            elif len(entity) == 2:
                dim = 'line'
            else:
                dim = 'plane'

            # compute centroid
            cent = [0.0, 0.0, 0.0]
            for at in entity:
                cent = add(cent, self.xyz(at))
            cent = scale(cent, 1.0 / len(entity))
            text += '  Centroid:      %14.8f %14.8f %14.8f                  [Angstrom]\n' % \
                (cent[0] * constants.bohr2angstroms,
                 cent[1] * constants.bohr2angstroms,
                 cent[2] * constants.bohr2angstroms)
            text += '  Centroid:      %14.8f %14.8f %14.8f                  [Bohr]\n' % \
                (cent[0], cent[1], cent[2])

            if dim == 'point':
                summ.append({'dim': dim, 'geo': cent, 'cent': cent})
                # TODO: figure out if should be using mass-weighted

            self.translate(scale(cent, -1))
            evals, evecs = self.inertial_system_partial(entity, masswt=False)
            midx = evals.index(max(evals))

            text += '  Normal Vector: %14.8f %14.8f %14.8f                  [unit]\n' % \
                (evecs[0][midx], evecs[1][midx], evecs[2][midx])
            text += '  Normal Vector: %14.8f %14.8f %14.8f                  [unit]\n' % \
                (evecs[0][midx] + cent[0], evecs[1][midx] + cent[1], evecs[2][midx] + cent[2])
            xplane = [evecs[0][midx], evecs[1][midx], evecs[2][midx],
                -1.0 * (evecs[0][midx] * cent[0] + evecs[1][midx] * cent[1] + evecs[2][midx] * cent[2])]
            text += '  Eqn. of Plane: %14.8f %14.8f %14.8f %14.8f   [Ai + Bj + Ck + D = 0]\n' % \
                (xplane[0], xplane[1], xplane[2], xplane[3])
            dtemp = math.sqrt(evecs[0][midx] * evecs[0][midx] + evecs[1][midx] * evecs[1][midx] +
                              evecs[2][midx] * evecs[2][midx])
            hessplane = [evecs[0][midx] / dtemp, evecs[1][midx] / dtemp, evecs[2][midx] / dtemp, xplane[3] / dtemp]
            hessplane2 = [xplane[0] / dtemp, xplane[1] / dtemp, xplane[2] / dtemp, xplane[3] / dtemp]
            text += '  Eqn. of Plane: %14.8f %14.8f %14.8f %14.8f   [Ai + Bj + Ck + D = 0] H\n' % \
                (hessplane[0], hessplane[1], hessplane[2], hessplane[3])
            text += '  Eqn. of Plane: %14.8f %14.8f %14.8f %14.8f   [Ai + Bj + Ck + D = 0] H2\n' % \
                (hessplane2[0], hessplane2[1], hessplane2[2], hessplane2[3])

            self.translate(cent)

            if dim == 'plane':
                summ.append({'dim': dim, 'geo': xplane, 'cent': cent})

        #print summ
        text += """\n  ==> 1 (%s) vs. 2 (%s) <==\n\n""" % (summ[0]['dim'], summ[1]['dim'])

        #if summ[0]['dim'] == 'plane' and summ[1]['dim'] == 'point':
        #    cent = summ[1]['geo']
        #    plane = summ[0]['geo']
        #    print cent, plane
        #
        #    D = math.fabs(plane[0] * cent[0] + plane[1] * cent[1] + plane[2] * cent[2] + plane[3]) / \
        #        math.sqrt(plane[0] * plane[0] + plane[1] * plane[1] + plane[2] * plane[2])
        #    text += '  Pt to Plane: %14.8f [Angstrom]\n' % (D * psi_bohr2angstroms)

        #if summ[0]['dim'] == 'plane' and summ[1]['dim'] == 'plane':
        if summ[0]['dim'] == 'plane' and (summ[1]['dim'] == 'plane' or summ[1]['dim'] == 'point'):
            cent1 = summ[0]['cent']
            cent2 = summ[1]['cent']
            plane1 = summ[0]['geo']
            #plane2 = summ[1]['geo']

            distCC = distance(cent1, cent2)
            text += '  Distance from Center of %s to Center of %s:                   %14.8f   [Angstrom]\n' % \
                ('2', '1', distCC * constants.bohr2angstroms)

            distCP = math.fabs(plane1[0] * cent2[0] + plane1[1] * cent2[1] + plane1[2] * cent2[2] + plane1[3])
            # distCP expression has a denominator that's one since plane constructed from unit vector
            text += '  Distance from Center of %s to Plane of %s:                    %14.8f   [Angstrom]\n' % \
                ('2', '1', distCP * constants.bohr2angstroms)

            distCPC = math.sqrt(distCC * distCC - distCP * distCP)
            text += '  Distance from Center of %s to Center of %s along Plane of %s:  %14.8f   [Angstrom]\n' % \
                ('2', '1', '1', distCPC * constants.bohr2angstroms)

        print(text)

#        text = "        Interatomic Distances (Angstroms)\n\n"
#        for i in range(self.natom()):
#            for j in range(i + 1, self.natom()):
#                eij = sub(self.xyz(j), self.xyz(i))
#                dist = norm(eij) * psi_bohr2angstroms
#                text += "        Distance %d to %d %-8.3lf\n" % (i + 1, j + 1, dist)
#        text += "\n\n"
#        return text

    def rotor_type(self, tol=FULL_PG_TOL):
        """Returns the rotor type.

        >>> H2OH2O.rotor_type()
        RT_ASYMMETRIC_TOP

        """
        evals, evecs = diagonalize3x3symmat(self.inertia_tensor())
        evals = sorted(evals)

        rot_const = [
            1.0 / evals[0] if evals[0] > 1.0e-6 else 0.0, 1.0 / evals[1] if evals[1] > 1.0e-6 else 0.0,
            1.0 / evals[2] if evals[2] > 1.0e-6 else 0.0
        ]

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
            rotor_type = 'RT_LINEAR'  # 0  <  IB == IC      inf > B == C
        elif degen == 2:
            rotor_type = 'RT_SPHERICAL_TOP'  # IA == IB == IC       A == B == C
        elif degen == 1:
            if (rot_const[1] - rot_const[2]) < 1.0e-6:
                rotor_type = 'RT_PROLATE_SYMMETRIC_TOP'  # IA <  IB == IC       A >  B == C
            elif (rot_const[0] - rot_const[1]) < 1.0e-6:
                rotor_type = 'RT_OBLATE_SYMMETRIC_TOP'  # IA == IB <  IC       A == B >  C
        else:
            rotor_type = 'RT_ASYMMETRIC_TOP'  # IA <  IB <  IC       A  > B >  C
        return rotor_type

    def center_of_charge(self):
        """Computes center of charge of molecule (does not translate molecule).

        >>> H2OH2O.center_of_charge()
        [-0.073339893272065401, 0.002959783555632145, 0.0]

        """
        ret = [0.0, 0.0, 0.0]
        total_c = 0.0

        for at in range(self.natom()):
            c = self.charge(at)
            ret = add(ret, scale(self.xyz(at), c))
            total_c += c

        ret = scale(ret, 1.0 / total_c)
        return ret

    def move_to_coc(self):
        """Moves molecule to center of charge

        """
        coc = scale(self.center_of_charge(), -1.0)
        self.translate(coc)

    def rotational_symmetry_number(self):
        """Number of unique orientations of the rigid molecule that only interchange identical atoms.

        Notes
        -----
        Source http://cccbdb.nist.gov/thermo.asp (search "symmetry number")

        """
        pg = self.get_full_point_group()
        pg = self.full_point_group_with_n()
        if pg in ['ATOM', 'C1', 'Ci', 'Cs', 'C_inf_v']:
            sigma = 1
        elif pg == 'D_inf_h':
            sigma = 2
        elif pg in ['T', 'Td']:
            sigma = 12
        elif pg == 'Oh':
            sigma = 24
        elif pg == 'Ih':
            sigma = 60
        elif pg in ['Cn', 'Cnv', 'Cnh']:
            sigma = self.full_pg_n()
        elif pg in ['Dn', 'Dnd', 'Dnh']:
            sigma = 2 * self.full_pg_n()
        elif pg == 'Sn':
            sigma = self.full_pg_n() / 2
        else:
            raise ValidationError("Can't ID full symmetry group: " + pg)

        return sigma

    def axis_representation(self, zero=1e-8):
        """Molecule vs. laboratory frame representation (e.g., IR or IIIL).

        Parameters
        ----------
        zero : float, optional
            Screen for inertial tensor elements

        Returns
        -------
        str
            Representation code IR, IIR, IIIR, IL, IIL, IIIL. When
            molecule not in inertial frame, string is prefixed by "~".

        Notes
        -----
        Not carefully handling degenerate inertial elements.

        """
        it = self.inertia_tensor(zero=zero)
        Iidx = np.argsort(np.diagonal(it))
        if np.array_equal(Iidx, np.asarray([1, 2, 0])):
            ar = 'IR'
        elif np.array_equal(Iidx, np.asarray([2, 0, 1])):
            ar = 'IIR'
        elif np.array_equal(Iidx, np.asarray([0, 1, 2])):
            ar = 'IIIR'
        elif np.array_equal(Iidx, np.asarray([2, 1, 0])):
            ar = 'IL'
        elif np.array_equal(Iidx, np.asarray([0, 2, 1])):
            ar = 'IIL'
        elif np.array_equal(Iidx, np.asarray([1, 0, 2])):
            ar = 'IIIL'

        # if inertial tensor has non-zero off-diagonals, this whole classification is iffy
        if np.count_nonzero(it - np.diag(np.diagonal(it))):
            ar = '~' + ar

        return ar

    def to_arrays(self, dummy: bool = False, ghost_as_dummy: bool = False) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """Exports coordinate info into NumPy arrays.

        Parameters
        ----------
        dummy
            Whether or not to include dummy atoms in returned arrays.
        ghost_as_dummy
            Whether or not to treat ghost atoms as dummies.

        Returns
        -------
        geom, mass, elem, elez, uniq : numpy.ndarray, numpy.ndarray, numpy.ndarray, numpy.ndarray, numpy.ndarray
            (nat, 3) geometry [a0].
            (nat,) mass [u].
            (nat,) element symbol.
            (nat,) atomic number.
            (nat,) hash of element symbol and mass.
            Note that coordinate, orientation, and element information is
            preserved but fragmentation, chgmult, and dummy/ghost is lost.

        Usage
        -----
        geom, mass, elem, elez, uniq = molinstance.to_arrays()

        """
        self.update_geometry()

        if dummy:
            if isinstance(self, Molecule):
                # normal qcdb.Molecule
                geom = self.full_geometry(np_out=True)
            else:
                # psi4.core.Molecule
                geom = np.array(self.full_geometry())
            mass = np.asarray(
                [(0. if (ghost_as_dummy and self.fZ(at) == 0) else self.fmass(at)) for at in range(self.nallatom())])
            elem = np.asarray(
                ['X' if (ghost_as_dummy and self.fZ(at) == 0) else self.fsymbol(at) for at in range(self.nallatom())])
            elez = np.asarray(
                [0 if (ghost_as_dummy and self.fZ(at) == 0) else self.fZ(at) for at in range(self.nallatom())])
            uniq = np.asarray([
                hashlib.sha1((str(elem[at]) + str(mass[at])).encode('utf-8')).hexdigest()
                for at in range(self.nallatom())
            ])
        else:
            if isinstance(self, Molecule):
                # normal qcdb.Molecule
                geom = self.geometry(np_out=True)
            else:
                # psi4.core.Molecule
                geom = np.array(self.geometry())
            mass = np.asarray([self.mass(at) for at in range(self.natom())])
            elem = np.asarray([self.symbol(at) for at in range(self.natom())])
            elez = np.asarray([self.Z(at) for at in range(self.natom())])
            uniq = np.asarray([
                hashlib.sha1((str(elem[at]) + str(mass[at])).encode('utf-8')).hexdigest() for at in range(self.natom())
            ])

        return geom, mass, elem, elez, uniq

    @staticmethod
    def from_string(molstr,
                    dtype=None,
                    name=None,
                    fix_com=None,
                    fix_orientation=None,
                    fix_symmetry=None,
                    return_dict=False,
                    enable_qm=True,
                    enable_efp=True,
                    missing_enabled_return_qm='none',
                    missing_enabled_return_efp='none',
                    verbose=1):
        molrec = qcel.molparse.from_string(
            molstr=molstr,
            dtype=dtype,
            name=name,
            fix_com=fix_com,
            fix_orientation=fix_orientation,
            fix_symmetry=fix_symmetry,
            return_processed=False,
            enable_qm=enable_qm,
            enable_efp=enable_efp,
            missing_enabled_return_qm=missing_enabled_return_qm,
            missing_enabled_return_efp=missing_enabled_return_efp,
            verbose=verbose)
        if return_dict:
            return Molecule.from_dict(molrec['qm']), molrec
        else:
            return Molecule.from_dict(molrec['qm'])

    @staticmethod
    def from_arrays(geom=None,
                    elea=None,
                    elez=None,
                    elem=None,
                    mass=None,
                    real=None,
                    elbl=None,
                    name=None,
                    units='Angstrom',
                    input_units_to_au=None,
                    fix_com=False,
                    fix_orientation=False,
                    fix_symmetry=None,
                    fragment_separators=None,
                    fragment_charges=None,
                    fragment_multiplicities=None,
                    molecular_charge=None,
                    molecular_multiplicity=None,
                    comment=None,
                    provenance=None,
                    connectivity=None,
                    missing_enabled_return='error',
                    tooclose=0.1,
                    zero_ghost_fragments=False,
                    nonphysical=False,
                    mtol=1.e-3,
                    verbose=1,
                    return_dict=False):
        """Construct Molecule from unvalidated arrays and variables.

        Light wrapper around :py:func:`~qcelemental.molparse.from_arrays`
        that is a full-featured constructor to dictionary representa-
        tion of Molecule. This follows one step further to return
        Molecule instance.

        Parameters
        ----------
        See :py:func:`~qcelemental.molparse.from_arrays`.
        return_dict : bool, optional
            Additionally return Molecule dictionary intermediate.

        Returns
        -------
        mol : :py:class:`~qcdb.Molecule`
        molrec : dict, optional
            Dictionary representation of instance.
            Only provided if `return_dict` is True.

        """
        molrec = qcel.molparse.from_arrays(
            geom=geom,
            elea=elea,
            elez=elez,
            elem=elem,
            mass=mass,
            real=real,
            elbl=elbl,
            name=name,
            units=units,
            input_units_to_au=input_units_to_au,
            fix_com=fix_com,
            fix_orientation=fix_orientation,
            fix_symmetry=fix_symmetry,
            fragment_separators=fragment_separators,
            fragment_charges=fragment_charges,
            fragment_multiplicities=fragment_multiplicities,
            molecular_charge=molecular_charge,
            molecular_multiplicity=molecular_multiplicity,
            comment=comment,
            provenance=provenance,
            connectivity=connectivity,
            domain='qm',
            missing_enabled_return=missing_enabled_return,
            tooclose=tooclose,
            zero_ghost_fragments=zero_ghost_fragments,
            nonphysical=nonphysical,
            mtol=mtol,
            verbose=verbose)
        if return_dict:
            return Molecule.from_dict(molrec), molrec
        else:
            return Molecule.from_dict(molrec)

    def to_string(self, dtype, units=None, atom_format=None, ghost_format=None, width=17, prec=12):
        """Format a string representation of QM molecule."""

        molrec = self.to_dict(np_out=True)

        # flip zeros
        molrec['geom'][np.abs(molrec['geom']) < 5**(-(prec))] = 0

        smol = qcel.molparse.to_string(
            molrec,
            dtype=dtype,
            units=units,
            atom_format=atom_format,
            ghost_format=ghost_format,
            width=width,
            prec=prec)
        return smol

    def run_dftd3(self, func: Optional[str] = None, dashlvl: Optional[str] = None, dashparam: Optional[Dict] = None, dertype: Union[int, str, None] = None, verbose: int = 1):
        """Compute dispersion correction via Grimme's DFTD3 program.

        Parameters
        ----------
        func
            Name of functional (func only, func & disp, or disp only) for
            which to compute dispersion (e.g., blyp, BLYP-D2, blyp-d3bj,
            blyp-d3(bj), hf+d). Any or all parameters initialized
            from `dashcoeff[dashlvl][func]` can be overwritten via
            `dashparam`.
        dashlvl
            Name of dispersion correction to be applied (e.g., d, D2,
            d3(bj), das2010). Must be key in `dashcoeff` or "alias" or
            "formal" to run.
        dashparam
            Values for the same keys as `dashcoeff[dashlvl]['default']`
            used to override any or all values initialized by `func`.
            Extra parameters will error.
        dertype
            Maximum derivative level at which to run DFTD3. For large
            molecules, energy-only calculations can be significantly more
            efficient. Influences return values, see below.
        verbose
            Amount of printing.

        Returns
        -------
        energy : float
            When `dertype=0`, energy [Eh].
        gradient : ~numpy.ndarray
            When `dertype=1`, (nat, 3) gradient [Eh/a0].
        (energy, gradient) : tuple of float and ~numpy.ndarray
            When `dertype=None`, both energy [Eh] and (nat, 3) gradient [Eh/a0].

        """
        import qcengine as qcng

        if dertype is None:
            derint, derdriver = -1, 'gradient'
        else:
            derint, derdriver = parse_dertype(dertype, max_derivative=1)

        resinp = {
            'molecule': self.to_schema(dtype=2),
            'driver': derdriver,
            'model': {
                'method': func,
                'basis': '(auto)',
            },
            'keywords': {
                'verbose': verbose,
            },
        }
        if dashlvl:
            resinp['keywords']['level_hint'] = dashlvl
        if dashparam:
            resinp['keywords']['params_tweaks'] = dashparam
        jobrec = qcng.compute(resinp, 'dftd3', raise_error=True)
        jobrec = jobrec.dict()

        # hack as not checking type GRAD
        for k, qca in jobrec['extras']['qcvars'].items():
            if isinstance(qca, (list, np.ndarray)):
                jobrec['extras']['qcvars'][k] = np.array(qca).reshape(-1, 3)

        if isinstance(self, Molecule):
            pass
        else:
            from psi4 import core

            for k, qca in jobrec['extras']['qcvars'].items():
                if not isinstance(qca, (list, np.ndarray)):
                    core.set_variable(k, float(qca))

        if derint == -1:
            return (float(jobrec['extras']['qcvars']['DISPERSION CORRECTION ENERGY']),
                    jobrec['extras']['qcvars']['DISPERSION CORRECTION GRADIENT'])
        elif derint == 0:
            return float(jobrec['extras']['qcvars']['DISPERSION CORRECTION ENERGY'])
        elif derint == 1:
            return jobrec['extras']['qcvars']['DISPERSION CORRECTION GRADIENT']

    def run_dftd4(self, func: Optional[str] = None, dashlvl: Optional[str] = None, dashparam: Optional[Dict] = None, dertype: Union[int, str, None] = None, verbose: int = 1):
        """Compute dispersion correction via Grimme's DFTD4 program.

        Parameters
        ----------
        func
            Name of functional (func only, func & disp, or disp only) for
            which to compute dispersion (e.g., blyp, BLYP-D2, blyp-d3bj,
            blyp-d3(bj), hf+d). Unlike run_dftd3, ``func`` overwrites any
            parameter initialized via `dashparam`.
        dashlvl
            Name of dispersion correction to be applied (e.g., d, D2,
            d3(bj), das2010). Must be key in `dashcoeff` or "alias" or
            "formal" to run.
        dashparam
            Values for the same keys as `dashcoeff[dashlvl]['default']`
            used to provide custom values. Unlike run_dftd3, will not have
            effect if `func` given. Must provide all parameters.
            Extra parameters will error.
        dertype
            Maximum derivative level at which to run DFTD3. For large
            molecules, energy-only calculations can be significantly more
            efficient. Influences return values, see below.
        verbose
            Amount of printing.

        Returns
        -------
        energy : float
            When `dertype=0`, energy [Eh].
        gradient : ndarray
            When `dertype=1`, (nat, 3) gradient [Eh/a0].
        (energy, gradient) : tuple of float and ndarray
            When `dertype=None`, both energy [Eh] and (nat, 3) gradient [Eh/a0].

        Notes
        -----
        This function wraps the QCEngine dftd4 harness which wraps the internal DFTD4 Python API.
        As such, the upstream convention of `func` trumping `dashparam` holds, rather than the
        :py:func:`run_dftd3` behavior of `dashparam` extending or overriding `func`.

        """
        import qcengine as qcng

        if dertype is None:
            derint, derdriver = -1, 'gradient'
        else:
            derint, derdriver = parse_dertype(dertype, max_derivative=1)

        resinp = {
            'molecule': self.to_schema(dtype=2),
            'driver': derdriver,
            'model': {
                'method': func,
                'basis': '(auto)',
            },
            'keywords': {
                'verbose': verbose,
            },
        }
        if dashlvl:
            resinp['keywords']['level_hint'] = dashlvl
        if dashparam:
            resinp['keywords']['params_tweaks'] = dashparam

        jobrec = qcng.compute(resinp, 'dftd4', raise_error=True)
        jobrec = jobrec.dict()

        # hack as not checking type GRAD
        for k, qca in jobrec['extras']['qcvars'].items():
            if isinstance(qca, (list, np.ndarray)):
                jobrec['extras']['qcvars'][k] = np.array(qca).reshape(-1, 3)

        if isinstance(self, Molecule):
            pass
        else:
            from psi4 import core

            for k, qca in jobrec['extras']['qcvars'].items():
                if not isinstance(qca, (list, np.ndarray)):
                    core.set_variable(k, float(qca))

        if derint == -1:
            return (float(jobrec['extras']['qcvars']['DISPERSION CORRECTION ENERGY']),
                    jobrec['extras']['qcvars']['DISPERSION CORRECTION GRADIENT'])
        elif derint == 0:
            return float(jobrec['extras']['qcvars']['DISPERSION CORRECTION ENERGY'])
        elif derint == 1:
            return jobrec['extras']['qcvars']['DISPERSION CORRECTION GRADIENT']

    def run_gcp(self, func: Optional[str] = None, dertype: Union[int, str, None] = None, verbose: int = 1):
        """Compute geometrical BSSE correction via Grimme's GCP program.

        Function to call Grimme's GCP program
        https://www.chemie.uni-bonn.de/pctc/mulliken-center/software/gcp/gcp
        to compute an a posteriori geometrical BSSE correction to *self* for
        several HF, generic DFT, and specific HF-3c and PBEh-3c method/basis
        combinations, *func*. Returns energy if *dertype* is 0, gradient
        if *dertype* is 1, else tuple of energy and gradient if *dertype*
        unspecified. The gcp executable must be independently compiled and
        found in :envvar:`PATH` or :envvar:`PSIPATH`. *self* may be either a
        qcdb.Molecule (sensibly) or a psi4.Molecule (works b/c psi4.Molecule
        has been extended by this method py-side and only public interface
        fns used) or a string that can be instantiated into a qcdb.Molecule.

        Parameters
        ----------
        func
            Name of method/basis combination or composite method for which to compute the correction
           (e.g., HF/cc-pVDZ, DFT/def2-SVP, HF3c, PBEh3c).
        dertype
            Maximum derivative level at which to run GCP. For large
            molecules, energy-only calculations can be significantly more
            efficient. Influences return values, see below.
        verbose
            Amount of printing. Unused at present.

        Returns
        -------
        energy : float
            When `dertype=0`, energy [Eh].
        gradient : ndarray
            When `dertype=1`, (nat, 3) gradient [Eh/a0].
        (energy, gradient) : tuple of float and ndarray
            When `dertype=None`, both energy [Eh] and (nat, 3) gradient [Eh/a0].

        """
        import qcengine as qcng

        if dertype is None:
            derint, derdriver = -1, 'gradient'
        else:
            derint, derdriver = parse_dertype(dertype, max_derivative=1)

        resinp = {
            'molecule': self.to_schema(dtype=2),
            'driver': derdriver,
            'model': {
                'method': func,
                'basis': '(auto)',
            },
            'keywords': {
                'verbose': verbose,
            },
        }
        try:
            jobrec = qcng.compute(resinp, 'mctc-gcp', raise_error=True)
        except qcng.exceptions.ResourceError:
            jobrec = qcng.compute(resinp, 'gcp', raise_error=True)

        jobrec = jobrec.dict()

        # hack (instead of checking dertype GRAD) to collect `(nat, 3)` ndarray of gradient if present
        for variable_name, qcv in jobrec['extras']['qcvars'].items():
            if isinstance(qcv, (list, np.ndarray)):
                jobrec['extras']['qcvars'][variable_name] = np.array(qcv).reshape(-1, 3)

        if isinstance(self, Molecule):
            pass
        else:
            from psi4 import core

            for variable_name, qcv in jobrec['extras']['qcvars'].items():
                if not isinstance(qcv, (list, np.ndarray)):
                    core.set_variable(variable_name, float(qcv))

        if derint == -1:
            return (float(jobrec['extras']['qcvars']['GCP CORRECTION ENERGY']),
                    jobrec['extras']['qcvars']['GCP CORRECTION GRADIENT'])
        elif derint == 0:
            return float(jobrec['extras']['qcvars']['GCP CORRECTION ENERGY'])
        elif derint == 1:
            return jobrec['extras']['qcvars']['GCP CORRECTION GRADIENT']

    @staticmethod
    def from_schema(molschema, return_dict=False, verbose=1):
        """Construct Molecule from non-Psi4 schema.

        Light wrapper around :py:func:`~qcdb.Molecule.from_arrays`.

        Parameters
        ----------
        molschema : dict
            Dictionary form of Molecule following known schema.
        return_dict : bool, optional
            Additionally return Molecule dictionary intermediate.
        verbose : int, optional
            Amount of printing.

        Returns
        -------
        mol : :py:class:`~qcdb.Molecule`
        molrec : dict, optional
            Dictionary representation of instance.
            Only provided if `return_dict` is True.

        """

        molrec = qcel.molparse.from_schema(molschema, verbose=verbose)

        if return_dict:
            return Molecule.from_dict(molrec), molrec
        else:
            return Molecule.from_dict(molrec)

    def to_schema(self, dtype, units='Bohr'):
        """Serializes instance into dictionary according to schema `dtype`."""

        molrec = self.to_dict(np_out=True)
        schmol = qcel.molparse.to_schema(molrec, dtype=dtype, units=units)
        return schmol

    def to_dict(self, force_c1=False, force_units=False, np_out=True):
        """Serializes instance into Molecule dictionary."""

        self.update_geometry()
        molrec = {}

        if self.name() not in ['', 'default']:
            molrec['name'] = self.name()

        if self.comment() not in ['', 'default']:
            molrec['comment'] = self.comment()

        # qcdb does not add prov, so rely upon all qcdb.Mol creation happening in molparse for this to return valid value (not [])
        molrec['provenance'] = copy.deepcopy(self.provenance())

        if self.connectivity() != []:
            molrec['connectivity'] = copy.deepcopy(self.connectivity())

        if force_units == 'Bohr':
            molrec['units'] = 'Bohr'
        elif force_units == 'Angstrom':
            molrec['units'] = 'Angstrom'
        else:
            units = self.units()
            molrec['units'] = units
            if units == 'Angstrom' and abs(self.input_units_to_au() * constants.bohr2angstroms - 1.) > 1.e-6:
                molrec['input_units_to_au'] = self.input_units_to_au()

        molrec['fix_com'] = self.com_fixed()
        molrec['fix_orientation'] = self.orientation_fixed()
        if force_c1:
            molrec['fix_symmetry'] = 'c1'
        elif self.symmetry_from_input():
            molrec['fix_symmetry'] = self.symmetry_from_input()

        # if self.has_zmatrix:
        #     moldict['zmat'] = self.zmat
        # TODO zmat, geometry_variables

        nat = self.natom()
        geom = np.array(self.geometry())  # [a0]
        if molrec['units'] == 'Angstrom':
            geom *= constants.bohr2angstroms  #self.input_units_to_au()
        molrec['geom'] = geom.reshape((-1))

        molrec['elea'] = np.array([self.mass_number(at) for at in range(nat)])
        molrec['elez'] = np.array([qcel.periodictable.to_Z(self.symbol(at)) for at in range(nat)])
        molrec['elem'] = np.array([self.symbol(at).capitalize() for at in range(nat)])
        molrec['mass'] = np.array([self.mass(at) for at in range(nat)])
        molrec['real'] = np.array([bool(self.Z(at)) for at in range(nat)])
        molrec['elbl'] = np.array([self.label(at)[len(self.symbol(at)):].lower() for at in range(nat)])

        fragments = [x[:] for x in self.get_fragments()]
        fragment_charges = [float(f) for f in self.get_fragment_charges()]
        fragment_multiplicities = [m for m in self.get_fragment_multiplicities()]

        # do trimming not performed in Molecule class b/c fragment_* member data never directly exposed
        for ifr, fr in reversed(list(enumerate(self.get_fragment_types()))):
            if fr == 'Ghost':
                fragment_charges[ifr] = 0.
                fragment_multiplicities[ifr] = 1
            elif fr == 'Absent':
                del fragment_charges[ifr]
                del fragment_multiplicities[ifr]
                # readjust atom indices for subsequent fragments
                renum = fragments[ifr][0]
                for iffr, ffr in enumerate(fragments):
                    if iffr <= ifr:
                        continue
                    lenfr = ffr[1] - ffr[0]
                    fragments[iffr] = [renum, renum + lenfr]
                    renum += lenfr
                del fragments[ifr]

        molrec['fragment_separators'] = [int(f[0]) for f in fragments[1:]]  # np.int --> int
        molrec['fragment_charges'] = fragment_charges
        molrec['fragment_multiplicities'] = fragment_multiplicities

        molrec['molecular_charge'] = float(self.molecular_charge())
        molrec['molecular_multiplicity'] = self.multiplicity()

        # * mass number (elea) untouched by qcdb.Molecule/psi4.core.Molecule and
        #   likely to be array of -1s, so let from_arrays fill in the values and
        #   (1) don't complain about the difference and
        #   (2) return the from_arrays filled-in values
        # * from.arrays is expecting speclabel "Co_userlbl" for elbl, but we're
        #   sending "_userlbl", hence speclabel=False
        # * from.arrays sets difference provenance than input mol
        forgive = ['elea', 'provenance']

        # * from_arrays and comparison lines below are quite unnecessary to
        #   to_dict, but is included as a check. in practice, only fills in mass
        #   numbers and heals user chgmult.
        try:
            validated_molrec = qcel.molparse.from_arrays(speclabel=False, verbose=0, domain='qm', **molrec)
        except qcel.ValidationError as err:
            # * this can legitimately happen if total chg or mult has been set
            #   independently b/c fragment chg/mult not reset. so try again.
            print(
                """Following warning is harmless if you've altered chgmult through `set_molecular_change` or `set_multiplicity`. Such alterations are an expert feature. Specifying in the original molecule string is preferred. Nonphysical masses may also trigger the warning."""
            )
            molrec['fragment_charges'] = [None] * len(fragments)
            molrec['fragment_multiplicities'] = [None] * len(fragments)
            validated_molrec = qcel.molparse.from_arrays(speclabel=False, nonphysical=True, verbose=0, domain='qm', **molrec)
            forgive.append('fragment_charges')
            forgive.append('fragment_multiplicities')
        compare_molrecs(validated_molrec, molrec, 'to_dict', atol=1.e-6, forgive=forgive, verbose=0)

        # from_arrays overwrites provenance
        validated_molrec['provenance'] = copy.deepcopy(molrec['provenance'])

        if not np_out:
            validated_molrec = qcel.util.unnp(validated_molrec)

        return validated_molrec

    @classmethod
    def from_dict(cls, molrec, verbose=1):

        mol = cls()
        mol._internal_from_dict(molrec=molrec, verbose=verbose)
        return mol

    def _internal_from_dict(self, molrec, verbose=1):
        """Constructs instance from fully validated and defaulted dictionary `molrec`."""

        # Compromises for qcdb.Molecule
        # * molecular_charge is int, not float
        # * fragment_charges are int, not float

        self.lock_frame = False

        if 'name' in molrec:
            self.set_name(molrec['name'])

        if 'comment' in molrec:
            self.set_comment(molrec['comment'])

        self.set_provenance(copy.deepcopy(molrec['provenance']))

        if 'connectivity' in molrec:
            self.set_connectivity(copy.deepcopy(molrec['connectivity']))

        self.set_units(molrec['units'])
        if 'input_units_to_au' in molrec:
            self.set_input_units_to_au(molrec['input_units_to_au'])

        if 'geom_unsettled' in molrec:
            nat = len(molrec['geom_unsettled'])
            unsettled = True

            for iat in range(nat):
                entry = molrec['geom_unsettled'][iat]
                label = molrec['elem'][iat] + molrec['elbl'][iat]
                Z = molrec['elez'][iat] * int(molrec['real'][iat])
                self.add_unsettled_atom(Z, entry, molrec['elem'][iat], molrec['mass'][iat], Z, label,
                                        molrec['elea'][iat])
            for var in molrec['variables']:
                self.set_geometry_variable(var[0], var[1])

        else:
            geom = np.array(molrec['geom']).reshape((-1, 3))
            nat = geom.shape[0]
            unsettled = False

            for iat in range(nat):
                x, y, z = geom[iat]
                label = molrec['elem'][iat] + molrec['elbl'][iat]
                Z = molrec['elez'][iat] * int(molrec['real'][iat])
                self.add_atom(Z, x, y, z, molrec['elem'][iat], molrec['mass'][iat], Z, label, molrec['elea'][iat])
                # TODO charge and 2nd elez site
                # TODO real back to type Ghost?

        # apparently py- and c- sides settled on a diff convention of 2nd of pair in fragments_
        fragment_separators = np.array(molrec['fragment_separators'], dtype=int)
        fragment_separators = np.insert(fragment_separators, 0, 0)
        fragment_separators = np.append(fragment_separators, nat)
        fragments = [[fragment_separators[ifr], fr - 1] for ifr, fr in enumerate(fragment_separators[1:])]

        self.set_fragment_pattern(fragments, ['Real'] * len(fragments), [int(f) for f in molrec['fragment_charges']],
                                  molrec['fragment_multiplicities'])

        self.set_molecular_charge(int(molrec['molecular_charge']))
        self.set_multiplicity(molrec['molecular_multiplicity'])

        self.fix_com(molrec['fix_com'])
        self.fix_orientation(molrec['fix_orientation'])
        if 'fix_symmetry' in molrec:
            # Save the user-specified symmetry, but don't set it as the point group
            # That step occurs in update_geometry, after the atoms are added
            self.PYsymmetry_from_input = molrec['fix_symmetry'].lower()

        ## hack to prevent update_geometry termination upon no atoms
        #if nat == 0:
        #    self.set_lock_frame(True)

        if not unsettled:
            self.update_geometry()

    def BFS(self,
            seed_atoms: List = None,
            bond_threshold: float = 1.20,
            return_arrays: bool = False,
            return_molecules: bool = False,
            return_molecule: bool = False):
        """Detect fragments among real atoms through a breadth-first search (BFS) algorithm.

        Parameters
        ----------
        self : qcdb.Molecule or psi4.core.Molecule
        seed_atoms
            List of lists of atoms (0-indexed) belonging to independent fragments.
            Useful to prompt algorithm or to define intramolecular fragments through
            border atoms. Example: `[[1, 0], [2]]`
        bond_threshold
            Factor beyond average of covalent radii to determine bond cutoff.
        return_arrays
            If `True`, also return fragments as list of arrays.
        return_molecules
            If True, also return fragments as list of Molecules.
        return_molecule
            If True, also return one big Molecule with fragmentation encoded.

        Returns
        -------
        bfs_map : list of lists
            Array of atom indices (0-indexed) of detected fragments.
        bfs_arrays : tuple of lists of ndarray, optional
            geom, mass, elem info per-fragment.
            Only provided if `return_arrays` is True.
        bfs_molecules : list of qcdb.Molecule or psi4.core.Molecule, optional
            List of molecules, each built from one fragment. Center and
            orientation of fragments is fixed so orientation info from `self` is
            not lost. Loses chgmult and ghost/dummy info from `self` and contains
            default chgmult.
            Only provided if `return_molecules` is True.
            Returned are of same type as `self`.
        bfs_molecule : qcdb.Molecule or psi4.core.Molecule, optional
            Single molecule with same number of real atoms as `self` with atoms
            reordered into adjacent fragments and fragment markers inserted.
            Loses ghost/dummy info from `self`; keeps total charge but not total mult.
            Only provided if `return_molecule` is True.
            Returned is of same type as `self`.

        Authors
        -------
        Original code from Michael S. Marshall, linear-scaling algorithm from
        Trent M. Parker, revamped by Lori A. Burns

        Notes
        -----
        Relies upon van der Waals radii and so faulty for close (especially hydrogen-bonded) fragments. See` `seed_atoms``.
        Any existing fragmentation info/chgmult encoded in ``self`` is lost.

        """
        self.update_geometry()
        if self.natom() != self.nallatom():
            raise ValidationError("""BFS not adapted for dummy atoms""")

        cgeom, cmass, celem, celez, cuniq = self.to_arrays()
        frag_pattern = BFS(cgeom, celez, seed_atoms=seed_atoms, bond_threshold=bond_threshold)
        outputs = [frag_pattern]

        if return_arrays:
            fgeoms = [cgeom[fr] for fr in frag_pattern]
            fmasss = [cmass[fr] for fr in frag_pattern]
            felems = [celem[fr] for fr in frag_pattern]
            outputs.append((fgeoms, fmasss, felems))

        if return_molecules:
            molrecs = [
                qcel.molparse.from_arrays(
                    geom=cgeom[fr],
                    mass=cmass[fr],
                    elem=celem[fr],
                    elez=celez[fr],
                    units='Bohr',
                    fix_com=True,
                    fix_orientation=True) for fr in frag_pattern
            ]
            if isinstance(self, Molecule):
                ret_mols = [Molecule.from_dict(molrec) for molrec in molrecs]
            else:
                from psi4 import core
                ret_mols = [core.Molecule.from_dict(molrec) for molrec in molrecs]
            outputs.append(ret_mols)

        if return_molecule:
            dcontig = qcel.molparse.contiguize_from_fragment_pattern(
                frag_pattern, geom=cgeom, elez=celez, elem=celem, mass=cmass)
            molrec = qcel.molparse.from_arrays(
                geom=dcontig['geom'],
                mass=dcontig['mass'],
                elem=dcontig['elem'],
                elez=dcontig['elez'],
                units='Bohr',
                molecular_charge=self.molecular_charge(),
                # molecular_multiplicity may not be conservable upon fragmentation
                #   potentially could do two passes and try to preserve it
                fix_com=self.com_fixed(),
                fix_orientation=self.orientation_fixed(),
                fix_symmetry=(None if self.symmetry_from_input() == '' else self.symmetry_from_input()),
                fragment_separators=dcontig['fragment_separators'])
            if isinstance(self, Molecule):
                ret_mol = Molecule.from_dict(molrec)
            else:
                from psi4 import core
                ret_mol = core.Molecule.from_dict(molrec)

            outputs.append(ret_mol)

        outputs = tuple(outputs)
        return (frag_pattern, ) + outputs[1:]

    def B787(concern_mol: Union[qcdbmol, psi4.core.Molecule],
             ref_mol: Union[qcdbmol, psi4.core.Molecule],
             do_plot: bool = False,
             verbose: int = 1,
             atoms_map: bool = False,
             run_resorting: bool = False,
             mols_align: bool = False,
             run_to_completion: bool = False,
             uno_cutoff: float = 1.e-3,
             run_mirror: bool = False):
        """Finds shift, rotation, and atom reordering of `concern_mol` that best
        aligns with `ref_mol`.

        Wraps :py:func:`qcelemental.molutil.B787` for :py:class:`psi4.driver.qcdb.Molecule` or
        :py:class:`psi4.core.Molecule`. Employs the Kabsch, Hungarian, and
        Uno algorithms to exhaustively locate the best alignment for
        non-oriented, non-ordered structures.

        Parameters
        ----------
        concern_mol
            Molecule of concern, to be shifted, rotated, and reordered into
            best coincidence with `ref_mol`.
        ref_mol
            Molecule to match.
        atoms_map
            Whether atom1 of `ref_mol` corresponds to atom1 of `concern_mol`, etc.
            If true, specifying `True` can save much time.
        mols_align
            Whether `ref_mol` and `concern_mol` have identical geometries by eye
            (barring orientation or atom mapping) and expected final RMSD = 0.
            If `True`, procedure is truncated when RMSD condition met, saving time.
        do_plot
            Pops up a mpl plot showing before, after, and ref geometries.
        run_to_completion
            Run reorderings to completion (past RMSD = 0) even if unnecessary because
            `mols_align=True`. Used to test worst-case timings.
        run_resorting
            Run the resorting machinery even if unnecessary because `atoms_map=True`.
        uno_cutoff
            TODO
        run_mirror
            Run alternate geometries potentially allowing best match to `ref_mol`
            from mirror image of `concern_mol`. Only run if system confirmed to
            be nonsuperimposable upon mirror reflection.

        Returns
        -------
        float, tuple, qcdb.Molecule or psi4.core.Molecule
            First item is RMSD [A] between `ref_mol` and the optimally aligned
            geometry computed.
            Second item is a AlignmentMill namedtuple with fields
            (shift, rotation, atommap, mirror) that prescribe the transformation
            from `concern_mol` and the optimally aligned geometry.
            Third item is a crude charge-, multiplicity-, fragment-less Molecule
            at optimally aligned (and atom-ordered) geometry. Return type
            determined by `concern_mol` type.

        """
        rgeom, rmass, relem, relez, runiq = ref_mol.to_arrays()
        cgeom, cmass, celem, celez, cuniq = concern_mol.to_arrays()

        rmsd, solution = qcel.molutil.B787(
            cgeom=cgeom,
            rgeom=rgeom,
            cuniq=cuniq,
            runiq=runiq,
            do_plot=do_plot,
            verbose=verbose,
            atoms_map=atoms_map,
            run_resorting=run_resorting,
            mols_align=mols_align,
            run_to_completion=run_to_completion,
            run_mirror=run_mirror,
            uno_cutoff=uno_cutoff)

        ageom, amass, aelem, aelez, auniq = solution.align_system(cgeom, cmass, celem, celez, cuniq, reverse=False)
        adict = qcel.molparse.from_arrays(
            geom=ageom,
            mass=amass,
            elem=aelem,
            elez=aelez,
            units='Bohr',
            molecular_charge=concern_mol.molecular_charge(),
            molecular_multiplicity=concern_mol.multiplicity(),
            fix_com=True,
            fix_orientation=True)
        if isinstance(concern_mol, Molecule):
            amol = Molecule.from_dict(adict)
        else:
            from psi4 import core
            amol = core.Molecule.from_dict(adict)

        compare_values(
            concern_mol.nuclear_repulsion_energy(),
            amol.nuclear_repulsion_energy(),
            4,
            'Q: concern_mol-->returned_mol NRE uncorrupted',
            verbose=verbose - 1)
        if mols_align:
            compare_values(
                ref_mol.nuclear_repulsion_energy(),
                amol.nuclear_repulsion_energy(),
                4,
                'Q: concern_mol-->returned_mol NRE matches ref_mol',
                verbose=verbose - 1)
            compare_integers(
                True,
                np.allclose(ref_mol.geometry(), amol.geometry(), atol=4),
                'Q: concern_mol-->returned_mol geometry matches ref_mol',
                verbose=verbose - 1)

        return rmsd, solution, amol

    def scramble(ref_mol: "Molecule",
                 do_shift: Union[bool, np.ndarray, List] = True,
                 do_rotate: Union[bool, np.ndarray, List[List]] = True,
                 do_resort: Union[bool, List] = True,
                 deflection: float = 1.0,
                 do_mirror: bool = False,
                 do_plot: bool = False,
                 run_to_completion: bool = False,
                 run_resorting: bool = False,
                 verbose: int = 1):
        """Tester for B787 by shifting, rotating, and atom shuffling `ref_mol` and
        checking that the aligner returns the opposite transformation.

        Parameters
        ----------
        ref_mol
            Molecule to perturb.
        do_shift
            Whether to generate a random atom shift on interval [-3, 3) in each
            dimension (`True`) or leave at current origin. To shift by a specified
            vector, supply a 3-element list.
        do_rotate
            Whether to generate a random 3D rotation according to algorithm of Arvo.
            To rotate by a specified matrix, supply a 9-element list of lists.
        do_resort
            Whether to shuffle atoms (`True`) or leave 1st atom 1st, etc. (`False`).
            To specify shuffle, supply a nat-element list of indices.
        deflection
            If `do_rotate`, how random a rotation: 0.0 is no change, 0.1 is small
            perturbation, 1.0 is completely random.
        do_mirror
            Whether to construct the mirror image structure by inverting y-axis.
        do_plot
            Pops up a mpl plot showing before, after, and ref geometries.
        run_to_completion
            By construction, scrambled systems are fully alignable (final RMSD=0).
            Even so, `True` turns off the mechanism to stop when RMSD reaches zero
            and instead proceed to worst possible time.
        run_resorting
            Even if atoms not shuffled, test the resorting machinery.
        verbose
            Print level.

        Returns
        -------
        None

        """
        rgeom, rmass, relem, relez, runiq = ref_mol.to_arrays()
        nat = rgeom.shape[0]

        perturbation = qcel.molutil.compute_scramble(
            rgeom.shape[0],
            do_shift=do_shift,
            do_rotate=do_rotate,
            deflection=deflection,
            do_resort=do_resort,
            do_mirror=do_mirror)
        cgeom, cmass, celem, celez, cuniq = perturbation.align_system(rgeom, rmass, relem, relez, runiq, reverse=True)
        cmol = Molecule.from_arrays(
            geom=cgeom,
            mass=cmass,
            elem=celem,
            elez=celez,
            units='Bohr',
            molecular_charge=ref_mol.molecular_charge(),
            molecular_multiplicity=ref_mol.multiplicity(),
            fix_com=True,
            fix_orientation=True)

        rmsd = np.linalg.norm(cgeom - rgeom) * constants.bohr2angstroms / np.sqrt(nat)
        if verbose >= 1:
            print('Start RMSD = {:8.4f} [A]'.format(rmsd))

        rmsd, solution, amol = cmol.B787(
            ref_mol,
            do_plot=do_plot,
            atoms_map=(not do_resort),
            run_resorting=run_resorting,
            mols_align=True,
            run_to_completion=run_to_completion,
            run_mirror=do_mirror,
            verbose=verbose)

        compare_integers(
            True, np.allclose(solution.shift, perturbation.shift, atol=6), 'shifts equiv', verbose=verbose - 1)
        if not do_resort:
            compare_integers(
                True,
                np.allclose(solution.rotation.T, perturbation.rotation),
                'rotations transpose',
                verbose=verbose - 1)
        if solution.mirror:
            compare_integers(True, do_mirror, 'mirror allowed', verbose=verbose - 1)

    def set_fragment_pattern(self, frl, frt, frc, frm):
        """Set fragment member data through public method analogous to psi4.core.Molecule"""

        if not (len(frl) == len(frt) == len(frc) == len(frm)):
            raise ValidationError("""Molecule::set_fragment_pattern: fragment arguments not of same length.""")

        self.fragments = frl
        self.fragment_types = frt
        self.fragment_charges = frc
        self.fragment_multiplicities = frm


# Attach methods to qcdb.Molecule class
from .parker import xyz2mol as _parker_xyz2mol_yo
Molecule.format_molecule_for_mol = _parker_xyz2mol_yo
