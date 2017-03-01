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
#import re
#import math
#import copy
#from periodictable import *
#from physconst import *
#from vecutil import *
#from exceptions import *
#from coordentry import *
import subprocess
import socket
import shutil
import random
from collections import defaultdict
from .libmintsmolecule import *


class Molecule(LibmintsMolecule):
    """Class to store python extensions to the MoleculeLibmints class.
    Multiple classes allows separation of libmints and extension methods.

    """

    def __init__(self, psi4molstr=None):
        """Initialize Molecule object from LibmintsMolecule"""
        LibmintsMolecule.__init__(self, psi4molstr)

        # The comment line
        self.tagline = ""

    def __str__(self):
        text = """  ==> qcdb Molecule %s <==\n\n""" % (self.name())
        text += """   => %s <=\n\n""" % (self.tagline)
        text += self.create_psi4_string_from_molecule()
        return text

#    def __getstate__(self):
#        print 'im being pickled'
#        return self.__dict__

#    def __setstate__(self, d):
#        print 'im being unpickled with these values', d
#        self.__dict__ = d

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
                raise ValidationError("""Molecule::init_with_xyz: given filename '%s' does not exist.""" % (xyzfilename))
            if os.stat(xyzfilename).st_size == 0:
                raise ValidationError("""Molecule::init_with_xyz: given filename '%s' is blank.""" % (xyzfilename))
            text = infile.readlines()

        xyz1 = re.compile(r"^\s*(\d+)\s*(bohr|au)?\s*$", re.IGNORECASE)
        xyz2 = re.compile(r'^\s*(-?\d+)\s+(\d+)\s+(.*)\s*$')
        NUMBER = "((?:[-+]?\\d*\\.\\d+(?:[DdEe][-+]?\\d+)?)|(?:[-+]?\\d+\\.\\d*(?:[DdEe][-+]?\\d+)?))"
        xyzN = re.compile(r'(?:\s*)([A-Z](?:[a-z])?)(?:\s+)' +
            NUMBER + '(?:\s+)' + NUMBER + '(?:\s+)' + NUMBER + '(?:\s*)', re.IGNORECASE)
        xyzC = re.compile(r'(?:\s*)(\d+\.?\d*)(?:\s+)' +
            NUMBER + '(?:\s+)' + NUMBER + '(?:\s+)' + NUMBER + '(?:\s*)', re.IGNORECASE)

        # Try to match the first line
        if xyz1.match(text[0]):
            fileNatom = int(xyz1.match(text[0]).group(1))
            if xyz1.match(text[0]).group(2) == None:
                fileUnits = 'Angstrom'
            else:
                fileUnits = 'Bohr'
        else:
            raise ValidationError("Molecule::init_with_xyz: Malformed first line\n%s" % (text[0]))

        # Try to match the second line
        if xyz2.match(text[1]):
            instance.set_molecular_charge(int(xyz2.match(text[1]).group(1)))
            instance.set_multiplicity(int(xyz2.match(text[1]).group(2)))
            instance.tagline = xyz2.match(text[1]).group(3).strip()
        else:
            instance.tagline = text[1].strip()

        # Next line begins the useful information.
        for i in range(fileNatom):
            try:
                if xyzN.match(text[2 + i]):

                    fileAtom = xyzN.match(text[2 + i]).group(1).upper()
                    fileX = float(xyzN.match(text[2 + i]).group(2))
                    fileY = float(xyzN.match(text[2 + i]).group(3))
                    fileZ = float(xyzN.match(text[2 + i]).group(4))

                    # Check that the atom symbol is valid
                    if not fileAtom in el2z:
                        raise ValidationError('Illegal atom symbol in geometry specification: %s' % (fileAtom))

                    # Add it to the molecule.
                    instance.add_atom(el2z[fileAtom], fileX, fileY, fileZ, fileAtom, el2mass[fileAtom], el2z[fileAtom])

                elif xyzC.match(text[2 + i]):

                    fileAtom = int(float(xyzC.match(text[2 + i]).group(1)))
                    fileX = float(xyzC.match(text[2 + i]).group(2))
                    fileY = float(xyzC.match(text[2 + i]).group(3))
                    fileZ = float(xyzC.match(text[2 + i]).group(4))

                    # Check that the atomic number is valid
                    if not fileAtom in z2el:
                        raise ValidationError('Illegal atom symbol in geometry specification: %d' % (fileAtom))

                    # Add it to the molecule.
                    instance.add_atom(fileAtom, fileX, fileY, fileZ, z2el[fileAtom], z2mass[fileAtom], fileAtom)

                else:
                    raise ValidationError("Molecule::init_with_xyz: Malformed atom information line %d." % (i + 3))
            except IndexError:
                raise ValidationError("Molecule::init_with_xyz: Expected atom in file at line %d.\n%s" % (i + 3, text[i + 2]))

        # We need to make 1 fragment with all atoms
        instance.fragments.append([0, fileNatom - 1])
        instance.fragment_types.append('Real')
        instance.fragment_charges.append(instance.molecular_charge())
        instance.fragment_multiplicities.append(instance.multiplicity())
        # Set the units properly
        instance.PYunits = fileUnits
        if fileUnits == 'Bohr':
            instance.input_units_to_au = 1.0
        elif fileUnits == 'Angstrom':
            instance.input_units_to_au = 1.0 / psi_bohr2angstroms

        instance.update_geometry()
        return instance

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
                raise ValidationError("""Molecule::init_with_mol2: given filename '%s' does not exist.""" % (xyzfilename))
            if os.stat(xyzfilename).st_size == 0:
                raise ValidationError("""Molecule::init_with_mol2: given filename '%s' is blank.""" % (xyzfilename))
            text = infile.readlines()

        # fixed-width regex ((?=[ ]*-?\d+)[ -\d]{5})
        v2000 = re.compile(r'^((?=[ ]*\d+)[ \d]{3})((?=[ ]*\d+)[ \d]{3})(.*)V2000\s*$')
        vend = re.compile(r'^\s*M\s+END\s*$')
        NUMBER = "((?:[-+]?\\d*\\.\\d+(?:[DdEe][-+]?\\d+)?)|(?:[-+]?\\d+\\.\\d*(?:[DdEe][-+]?\\d+)?))"
        xyzM = re.compile(r'^(?:\s*)' + NUMBER + '(?:\s+)' + NUMBER + '(?:\s+)' + NUMBER +
                          '(?:\s+)([A-Z](?:[a-z])?)(?:\s+)(.*)', re.IGNORECASE)

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
                    if not fileAtom in el2z:
                        raise ValidationError('Illegal atom symbol in geometry specification: %s' % (fileAtom))

                    # Add it to the molecule.
                    instance.add_atom(el2z[fileAtom], fileX, fileY, fileZ, fileAtom, el2mass[fileAtom], el2z[fileAtom])

                else:
                    raise ValidationError("Molecule::init_with_mol2: Malformed atom information line %d." % (i + 5))
            except IndexError:
                raise ValidationError("Molecule::init_with_mol2: Expected atom in file at line %d.\n%s" % (i + 5, text[i + 4]))

        # We need to make 1 fragment with all atoms
        instance.fragments.append([0, fileNatom - 1])
        instance.fragment_types.append('Real')
        instance.fragment_charges.append(instance.molecular_charge())
        instance.fragment_multiplicities.append(instance.multiplicity())
        # Set the units properly
        instance.PYunits = fileUnits
        if fileUnits == 'Bohr':
            instance.input_units_to_au = 1.0
        elif fileUnits == 'Angstrom':
            instance.input_units_to_au = 1.0 / psi_bohr2angstroms

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
        factor = 1.0 if self.PYunits == 'Angstrom' else psi_bohr2angstroms

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
                text += '%2s %17.12f %17.12f %17.12f\n' % ((self.symbol(i) if self.Z(i) else "Gh"), \
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
        import numpy as np
        factor = 1.0 if self.PYunits == 'Angstrom' else psi_bohr2angstroms
        self.update_geometry()

        # TODO fn title is format_mol... but return args not compatible
        geo = []
        for i in range(self.natom()):
            [x, y, z] = self.atoms[i].compute()
            geo.append([self.Z(i), x * factor, y * factor, z * factor])

        nparr = np.array(geo)
        return nparr if npobj else np.array_repr(nparr)

#    def save_string_for_psi4(self):
#        """Returns a string of Molecule formatted for psi4.
#        Includes fragments and reorienting, if specified.
#
#        >>> print H2OH2O.save_string_for_psi4()
#        6
#        0 1
#        O         -1.55100700      -0.11452000       0.00000000
#        H         -1.93425900       0.76250300       0.00000000
#        H         -0.59967700       0.04071200       0.00000000
#        --
#        0 1
#        @X         0.00000000       0.00000000       0.00000000
#        O          1.35062500       0.11146900       0.00000000
#        H          1.68039800      -0.37374100      -0.75856100
#        H          1.68039800      -0.37374100       0.75856100
#        units Angstrom
#
#        """
#        Nfr = 0
#        text = ""
#        for fr in range(self.nfragments()):
#            if self.fragment_types[fr] == 'Absent':
#                continue
#            if Nfr != 0:
#                text += """--\n"""
#            Nfr += 1
#            text += """%d %d\n""" % (self.fragment_charges[fr], self.fragment_multiplicities[fr])
#            for at in range(self.fragments[fr][0], self.fragments[fr][1] + 1):
#                geom = self.full_atoms[at].compute()
#                text += """%-3s  %16.8f %16.8f %16.8f\n""" % \
#                    (("" if self.fZ(at) else "@") + self.full_atoms[at].symbol(), \
#                    geom[0], geom[1], geom[2])
#        text += """units %s\n""" % (self.units().lower())
#        return text

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
        factor = 1.0 if self.PYunits == 'Angstrom' else psi_bohr2angstroms

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
        options = defaultdict(lambda: defaultdict(dict))
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

            factor = 1.0 if self.PYunits == 'Angstrom' else psi_bohr2angstroms
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
                text += "%s%s%d %d\n" % (
                    "" if Pfr == 0 else "--\n",
                    "#" if self.fragment_types[fr] == 'Ghost' or self.fragment_types[fr] == 'Absent' else "",
                    self.fragment_charges[fr], self.fragment_multiplicities[fr])
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
        factor = 1.0 if self.PYunits == 'Angstrom' else psi_bohr2angstroms
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
            text += '%-2s %17.12f %17.12f %17.12f\n' % (self.symbol(i), \
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
        factor = 1.0 if self.PYunits == 'Angstrom' else psi_bohr2angstroms
        #factor = 1.0 if self.PYunits == 'Bohr' else 1.0/psi_bohr2angstroms

        text = 'auto-generated by qcdb from molecule %s\n' % (self.tagline)

        # append atoms and coordentries
        for i in range(self.natom()):
            [x, y, z] = self.atoms[i].compute()
            text += '%-2s %17.12f %17.12f %17.12f\n' % ((self.symbol(i) if self.Z(i) else "GH"), \
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
        options = defaultdict(lambda: defaultdict(dict))
        options['CFOUR']['CFOUR_CHARGE']['value'] = self.molecular_charge()
        options['CFOUR']['CFOUR_MULTIPLICITY']['value'] = self.multiplicity()
        options['CFOUR']['CFOUR_UNITS']['value'] = 'ANGSTROM'
#        options['CFOUR']['CFOUR_UNITS']['value'] = 'BOHR'
        options['CFOUR']['CFOUR_COORDINATES']['value'] = 'CARTESIAN'
#        options['CFOUR']['CFOUR_SUBGROUP']['value'] = self.symmetry_from_input().upper()
#        print self.inertia_tensor()
#        print self.inertial_system()

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

        options = defaultdict(lambda: defaultdict(dict))
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
        options = defaultdict(lambda: defaultdict(dict))
        self.update_geometry()
        factor = 1.0 if self.PYunits == 'Angstrom' else psi_bohr2angstroms

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
                        text += '    {:2s}({:d}) {:> 17.12f} {:> 17.12f} {:> 17.12f}\n'.format(\
                                atom, fr + 1, x * factor, y * factor, z * factor)
                    else:
                        text += '    {:2s} {:> 17.12f} {:> 17.12f} {:> 17.12f}\n'.format(\
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
            factor = 1.0 if self.PYunits == 'Angstrom' else psi_bohr2angstroms
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
                    # this only distiguishes Real frags so Real/Ghost don't get 
                    #   fragmentation. may need to change
                    text += """--\n"""
                                         # >>>
                # any fragment chgmult here <<<
                if self.nactive_fragments() > 1:
                    text += """{}{} {}\n""".format(
                        '!' if self.fragment_types[fr] in ['Ghost', 'Absent'] else '',
                        self.fragment_charges[fr], self.fragment_multiplicities[fr])
                                          # >>>
                Pfr += 1
                for at in range(self.fragments[fr][0], self.fragments[fr][1] + 1):
                    if self.fragment_types[fr] == 'Absent' or self.fsymbol(at) == "X":
                        pass
                    else:
                        if self.fZ(at):
                            # label for real live atom <<<
                            text += """{:>3s} """.format(self.symbol(at))
                                                     # >>>
                        else:
                            # label for ghost atom <<<
                            text += """{:>3s} """.format(
                                'Gh' if mixedbas else ('@' + self.symbol(at)))                
                                                 # >>>
                        [x, y, z] = self.full_atoms[at].compute()
                        # Cartesian coordinates <<<
                        text += """{:>17.12f} {:>17.12f} {:>17.12f}\n""".format(
                            x * factor, y * factor, z * factor)
                                              # >>>
            # any general finishing notation here <<<
            text += '$end\n\n'
                                                # >>>

        # prepare molecule keywords to be set as c-side keywords
        options = defaultdict(lambda: defaultdict(dict))
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

#        # append units and any other non-default molecule keywords
#        text += "    units %-s\n" % ("Angstrom" if self.units() == 'Angstrom' else "Bohr")
#        if not self.PYmove_to_com:
#            text += "    no_com\n"
#        if self.PYfix_orientation:
#            text += "    no_reorient\n"

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
        options = defaultdict(lambda: defaultdict(dict))
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
        factor = 1.0 if self.PYunits == 'Angstrom' else psi_bohr2angstroms

        text = ""
        text += '%d %d %s\n' % (self.molecular_charge(), self.multiplicity(), self.tagline)

        for i in range(self.natom()):
            [x, y, z] = self.atoms[i].compute()
            text += '%4s %17.12f %17.12f %17.12f\n' % (("" if self.Z(i) else 'Bq') + self.symbol(i), \
                x * factor, y * factor, z * factor)
        return text
        pass

    #    if symm   print M2OUT "nosym\nnoorient\n";
    #    print DIOUT "angstrom\ngeometry={\n";

    def auto_fragments(self):
        """Detects fragments in an unfragmented molecule using BFS
        algorithm. Returns a new Molecule in Cartesian, fixed-geom
        (no variable values), no dummy-atom format. Any non-default
        charge and multiplicity assigned to first fragment.

        """
        if self.nfragments() != 1:
            print("""Molecule already fragmented so no further action by auto_fragments().""")
            return self

        flist = self.BFS()

        # form new molecule through a string since self may contain
        #   dummies or zmatrix specs that mayn't be valid with atom shuffling
        new_geom = '\n'

        if self.PYcharge_specified or self.PYmultiplicity_specified:
            new_geom = """\n   %d %d\n""" % (self.molecular_charge(), self.multiplicity())

        for fr in range(len(flist)):
            new_geom += "" if fr == 0 else "   --\n"
            for at in flist[fr]:
                geom = self.atoms[at].compute()
                new_geom += """%-4s """ % (("" if self.Z(at) else "@") + self.symbol(at))
                for j in range(3):
                    new_geom += """  %17.12f""" % (geom[j])
                new_geom += "\n"
        new_geom += "   units %s\n" % (self.units())
        if not self.PYmove_to_com:
            new_geom += "   no_com\n"
        if self.orientation_fixed():
            new_geom += "   no_reorient\n"

        subset = Molecule(new_geom)
        subset.update_geometry()
        return subset

    def BFS(self):
        """Perform a breadth-first search (BFS) on the real atoms
        in molecule, returning an array of atom indices of fragments.
        Relies upon van der Waals radii and so faulty for close
        (esp. hydrogen-bonded) fragments. Original code from
        Michael S. Marshall.

        """
        vdW_diameter = {
            #'H':  1.001 / 1.5,  # JMol
            'HE': 1.012 / 1.5,  # JMol
            'LI': 0.825 / 1.5,  # JMol
            'BE': 1.408 / 1.5,  # JMol
            #'B':  1.485 / 1.5,  # JMol
            #'C':  1.452 / 1.5,  # JMol
            #'N':  1.397 / 1.5,  # JMol
            #'O':  1.342 / 1.5,  # JMol
            #'F':  1.287 / 1.5,  # JMol
            'NE': 1.243 / 1.5,  # JMol
            'NA': 1.144 / 1.5,  # JMol
            'MG': 1.364 / 1.5,  # JMol
            'AL': 1.639 / 1.5,  # JMol
            #'SI': 1.716 / 1.5,  # JMol
            #'P':  1.705 / 1.5,  # JMol
            #'S':  1.683 / 1.5,  # JMol
            #'CL': 1.639 / 1.5,  # JMol
            'AR': 1.595 / 1.5,  # JMol

            'H': 1.06 / 1.5,  # Bondi JPC 68 441 (1964)
            'B': 1.65 / 1.5,  # Bondi JPC 68 441 (1964)
            'C': 1.53 / 1.5,  # Bondi JPC 68 441 (1964)
            'N': 1.46 / 1.5,  # Bondi JPC 68 441 (1964)
            'O': 1.42 / 1.5,  # Bondi JPC 68 441 (1964)
            'F': 1.40 / 1.5,  # Bondi JPC 68 441 (1964)
            'SI': 1.93 / 1.5,  # Bondi JPC 68 441 (1964)
            'P': 1.86 / 1.5,  # Bondi JPC 68 441 (1964)
            'S': 1.80 / 1.5,  # Bondi JPC 68 441 (1964)
            'CL': 1.75 / 1.5,  # Bondi JPC 68 441 (1964)
            'GE': 1.98 / 1.5,  # Bondi JPC 68 441 (1964)
            'AS': 1.94 / 1.5,  # Bondi JPC 68 441 (1964)
            'SE': 1.90 / 1.5,  # Bondi JPC 68 441 (1964)
            'BR': 1.87 / 1.5,  # Bondi JPC 68 441 (1964)
            'SN': 2.16 / 1.5,  # Bondi JPC 68 441 (1964)
            'SB': 2.12 / 1.5,  # Bondi JPC 68 441 (1964)
            'TE': 2.08 / 1.5,  # Bondi JPC 68 441 (1964)
            'I': 2.04 / 1.5,  # Bondi JPC 68 441 (1964)
            'XE': 2.05 / 1.5}  # Bondi JPC 68 441 (1964)

        Queue = []
        White = range(self.natom())  # untouched
        Black = []  # touched and all edges discovered
        Fragment = []  # stores fragments

        start = 0  # starts with the first atom in the list
        Queue.append(start)
        White.remove(start)

        # Simply start with the first atom, do a BFS when done, go to any
        #   untouched atom and start again iterate until all atoms belong
        #   to a fragment group
        while len(White) > 0 or len(Queue) > 0:  # Iterates to the next fragment
            Fragment.append([])

            while len(Queue) > 0:                # BFS within a fragment
                for u in Queue:                  # find all (still white) nearest neighbors to vertex u
                    for i in White:
                        dist = distance(self.xyz(i), self.xyz(u)) * psi_bohr2angstroms
                        if dist < vdW_diameter[self.symbol(u)] + vdW_diameter[self.symbol(i)]:
                            Queue.append(i)      # if you find you, put in the queue
                            White.remove(i)      # and remove it from the untouched list
                Queue.remove(u)                  # remove focus from Queue
                Black.append(u)
                Fragment[-1].append(int(u))      # add to group (0-indexed)
                Fragment[-1].sort()              # preserve original atom ordering

            if len(White) != 0:                  # can't move White -> Queue if no more exist
                Queue.append(White[0])
                White.remove(White[0])

        return Fragment

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
                (cent[0] * psi_bohr2angstroms, \
                 cent[1] * psi_bohr2angstroms, \
                 cent[2] * psi_bohr2angstroms)
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
            xplane = [evecs[0][midx], evecs[1][midx], evecs[2][midx], \
                -1.0 * (evecs[0][midx] * cent[0] + evecs[1][midx] * cent[1] + evecs[2][midx] * cent[2])]
            text += '  Eqn. of Plane: %14.8f %14.8f %14.8f %14.8f   [Ai + Bj + Ck + D = 0]\n' % \
                (xplane[0], xplane[1], xplane[2], xplane[3])
            dtemp = math.sqrt(evecs[0][midx] * evecs[0][midx] + evecs[1][midx] * evecs[1][midx] + evecs[2][midx] * evecs[2][midx])
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

#        if summ[0]['dim'] == 'plane' and summ[1]['dim'] == 'point':
#            cent = summ[1]['geo']
#            plane = summ[0]['geo']
#            print cent, plane
#
#            D = math.fabs(plane[0] * cent[0] + plane[1] * cent[1] + plane[2] * cent[2] + plane[3]) / \
#                math.sqrt(plane[0] * plane[0] + plane[1] * plane[1] + plane[2] * plane[2])
#            text += '  Pt to Plane: %14.8f [Angstrom]\n' % (D * psi_bohr2angstroms)

        #if summ[0]['dim'] == 'plane' and summ[1]['dim'] == 'plane':
        if summ[0]['dim'] == 'plane' and (summ[1]['dim'] == 'plane' or summ[1]['dim'] == 'point'):
            cent1 = summ[0]['cent']
            cent2 = summ[1]['cent']
            plane1 = summ[0]['geo']
            #plane2 = summ[1]['geo']

            distCC = distance(cent1, cent2)
            text += '  Distance from Center of %s to Center of %s:                   %14.8f   [Angstrom]\n' % \
                ('2', '1', distCC * psi_bohr2angstroms)

            distCP = math.fabs(plane1[0] * cent2[0] + plane1[1] * cent2[1] + plane1[2] * cent2[2] + plane1[3])
            # distCP expression has a denominator that's one since plane constructed from unit vector
            text += '  Distance from Center of %s to Plane of %s:                    %14.8f   [Angstrom]\n' % \
                ('2', '1', distCP * psi_bohr2angstroms)

            distCPC = math.sqrt(distCC * distCC - distCP * distCP)
            text += '  Distance from Center of %s to Center of %s along Plane of %s:  %14.8f   [Angstrom]\n' % \
                ('2', '1', '1', distCPC * psi_bohr2angstroms)

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

        rot_const = [1.0 / evals[0] if evals[0] > 1.0e-6 else 0.0,
                     1.0 / evals[1] if evals[1] > 1.0e-6 else 0.0,
                     1.0 / evals[2] if evals[2] > 1.0e-6 else 0.0]

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
            rotor_type = 'RT_LINEAR'                     # 0  <  IB == IC      inf > B == C
        elif degen == 2:
            rotor_type = 'RT_SPHERICAL_TOP'              # IA == IB == IC       A == B == C
        elif degen == 1:
            if (rot_const[1] - rot_const[2]) < 1.0e-6:
                rotor_type = 'RT_PROLATE_SYMMETRIC_TOP'  # IA <  IB == IC       A >  B == C
            elif (rot_const[0] - rot_const[1]) < 1.0e-6:
                rotor_type = 'RT_OBLATE_SYMMETRIC_TOP'   # IA == IB <  IC       A == B >  C
        else:
            rotor_type = 'RT_ASYMMETRIC_TOP'             # IA <  IB <  IC       A  > B >  C
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

# Attach methods to qcdb.Molecule class
from .interface_dftd3 import run_dftd3 as _dftd3_qcdb_yo
Molecule.run_dftd3 = _dftd3_qcdb_yo
from .parker import xyz2mol as _parker_xyz2mol_yo
Molecule.format_molecule_for_mol2 = _parker_xyz2mol_yo
from .parker import bond_profile as _parker_bondprofile_yo
Molecule.bond_profile = _parker_bondprofile_yo
from .interface_gcp import run_gcp as _gcp_qcdb_yo
Molecule.run_gcp = _gcp_qcdb_yo
