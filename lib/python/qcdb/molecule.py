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
import collections
from libmintsmolecule import *


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
        text += self.save_string_for_psi4()
        return text

#    def __getstate__(self):
#        print 'im being pickled'
#        return self.__dict__

#    def __setstate__(self, d):
#        #print 'im being unpickled with these values', d
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
        xyzN = re.compile(r"(?:\s*)([A-Z](?:[a-z])?)(?:\s+)(-?\d+\.\d+)(?:\s+)(-?\d+\.\d+)(?:\s+)(-?\d+\.\d+)(?:\s*)", re.IGNORECASE)
        xyzC = re.compile(r"(?:\s*)(\d+\.?\d*)(?:\s+)(-?\d+\.\d+)(?:\s+)(-?\d+\.\d+)(?:\s+)(-?\d+\.\d+)(?:\s*)", re.IGNORECASE)

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
                    instance.add_atom(el2z[fileAtom], fileX, fileY, fileZ, fileAtom, el2masses[fileAtom], el2z[fileAtom])

                elif xyzC.match(text[2 + i]):

                    fileAtom = int(float(xyzC.match(text[2 + i]).group(1)))
                    fileX = float(xyzC.match(text[2 + i]).group(2))
                    fileY = float(xyzC.match(text[2 + i]).group(3))
                    fileZ = float(xyzC.match(text[2 + i]).group(4))

                    # Check that the atomic number is valid
                    if not fileAtom in z2el:
                        raise ValidationError('Illegal atom symbol in geometry specification: %d' % (fileAtom))

                    # Add it to the molecule.
                    instance.add_atom(fileAtom, fileX, fileY, fileZ, z2el[fileAtom], z2masses[fileAtom], fileAtom)

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
        text += '%d %d %s\n' % (self.molecular_charge(), self.multiplicity(), self.tagline)

        for i in range(self.natom()):
            [x, y, z] = self.atoms[i].compute()
            if save_ghosts or self.Z(i):
                text += '%2s %17.12f %17.12f %17.12f\n' % ((self.symbol(i) if self.Z(i) else "Gh"), \
                    x * factor, y * factor, z * factor)
        return text

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

    def format_molecule_for_qchem(self):
        """

        """
        factor = 1.0 if self.PYunits == 'Angstrom' else psi_bohr2angstroms

        text = ""
        text += '%d %d %s\n' % (self.molecular_charge(), self.multiplicity(), self.tagline)

        for i in range(self.natom()):
            [x, y, z] = self.atoms[i].compute()
            text += '%2s %17.12f %17.12f %17.12f\n' % ((self.symbol(i) if self.Z(i) else "Gh"), \
                x * factor, y * factor, z * factor)
        return text
        pass

    def format_molecule_for_molpro(self):
        """

        """
        factor = 1.0 if self.PYunits == 'Angstrom' else psi_bohr2angstroms
        # TODO keep fix_or?
        self.fix_orientation(True)
        self.PYmove_to_com = False
        self.update_geometry()

        text = ""
        text += 'angstrom\n'
        text += 'geometry={\n'

        for fr in range(self.nfragments()):
            if self.fragment_types[fr] == 'Absent':
                pass
            else:
                for at in range(self.fragments[fr][0], self.fragments[fr][1] + 1):
                    [x, y, z] = self.atoms[at].compute()
                    text += '%2s %17.12f %17.12f %17.12f\n' % (self.symbol(at), \
                        x * factor, y * factor, z * factor)
        text += '}\n\n'
        text += 'SET,CHARGE=%d\n' % (self.molecular_charge())
        text += 'SET,SPIN=%d\n' % (self.multiplicity() - 1)  # Molpro wants (mult-1)

        textDummy = "dummy"
        for fr in range(self.nfragments()):
            if self.fragment_types[fr] == 'Ghost':
                for at in range(self.fragments[fr][0], self.fragments[fr][1] + 1):
                    textDummy += """,%d""" % (at + 1)  # Molpro atom numbering is 1-indexed
        textDummy += '\n'
        if len(textDummy) > 6:
            text += textDummy
        return text

    def format_molecule_for_cfour(self):
        """Function to print Molecule in a form readable by Cfour.

        """
        self.update_geometry()
        factor = 1.0 if self.PYunits == 'Angstrom' else psi_bohr2angstroms
        #factor = 1.0 if self.PYunits == 'Bohr' else 1.0/psi_bohr2angstroms

        text = 'auto-generated by qcdb from molecule %s\n' % (self.tagline)

        # append atoms and coordentries
        for fr in range(self.nfragments()):
            if self.fragment_types[fr] == 'Absent':
                pass
            else:
                for at in range(self.fragments[fr][0], self.fragments[fr][1] + 1):
                    [x, y, z] = self.atoms[at].compute()
                    text += '%-2s %17.12f %17.12f %17.12f\n' % ((self.symbol(at) if self.Z(at) else "GH"), \
                        x * factor, y * factor, z * factor)
        text += '\n'

        # prepare molecule keywords to be set as c-side keywords
        options = collections.defaultdict(lambda: collections.defaultdict(dict))
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

        options = collections.defaultdict(lambda: collections.defaultdict(dict))
        options['CFOUR']['CFOUR_BASIS']['value'] = 'SPECIAL'
        options['CFOUR']['CFOUR_SPHERICAL']['value'] = puream

        options['CFOUR']['CFOUR_BASIS']['clobber'] = True
        options['CFOUR']['CFOUR_SPHERICAL']['clobber'] = True

        options['CFOUR']['CFOUR_BASIS']['superclobber'] = True
        options['CFOUR']['CFOUR_SPHERICAL']['superclobber'] = True

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
            print 'Molecule already fragmented so no further action by auto_fragments().'
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

            'H':   1.06 / 1.5,  # Bondi JPC 68 441 (1964)
            'B':   1.65 / 1.5,  # Bondi JPC 68 441 (1964)
            'C':   1.53 / 1.5,  # Bondi JPC 68 441 (1964)
            'N':   1.46 / 1.5,  # Bondi JPC 68 441 (1964)
            'O':   1.42 / 1.5,  # Bondi JPC 68 441 (1964)
            'F':   1.40 / 1.5,  # Bondi JPC 68 441 (1964)
            'SI':  1.93 / 1.5,  # Bondi JPC 68 441 (1964)
            'P':   1.86 / 1.5,  # Bondi JPC 68 441 (1964)
            'S':   1.80 / 1.5,  # Bondi JPC 68 441 (1964)
            'CL':  1.75 / 1.5,  # Bondi JPC 68 441 (1964)
            'GE':  1.98 / 1.5,  # Bondi JPC 68 441 (1964)
            'AS':  1.94 / 1.5,  # Bondi JPC 68 441 (1964)
            'SE':  1.90 / 1.5,  # Bondi JPC 68 441 (1964)
            'BR':  1.87 / 1.5,  # Bondi JPC 68 441 (1964)
            'SN':  2.16 / 1.5,  # Bondi JPC 68 441 (1964)
            'SB':  2.12 / 1.5,  # Bondi JPC 68 441 (1964)
            'TE':  2.08 / 1.5,  # Bondi JPC 68 441 (1964)
            'I':   2.04 / 1.5,  # Bondi JPC 68 441 (1964)
            'XE':  2.05 / 1.5}  # Bondi JPC 68 441 (1964)

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
            print 'denom', dtemp
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

        print text

#        text = "        Interatomic Distances (Angstroms)\n\n"
#        for i in range(self.natom()):
#            for j in range(i + 1, self.natom()):
#                eij = sub(self.xyz(j), self.xyz(i))
#                dist = norm(eij) * psi_bohr2angstroms
#                text += "        Distance %d to %d %-8.3lf\n" % (i + 1, j + 1, dist)
#        text += "\n\n"
#        return text

    def grimme_dftd3(self, func=None, dashlvl=None, dashparam=None, verbosity=1):
        """Function to call Grimme's dftd3 program (http://toc.uni-muenster.de/DFTD3/)
        to compute the -D correction of level *dashlvl* using parameters for
        the functional *func*. The dictionary *dashparam* can be used to supply
        a full set of dispersion parameters in the absense of *func* or to supply
        individual overrides in the presence of *func*. The dftd3 executable must be
        independently compiled and found in :envvar:PATH.

        """
        # Parameters from http://toc.uni-muenster.de/DFTD3/ on September 25, 2012
        #   dict keys translated from Turbomole to Psi4 functional names
        dashcoeff = {
            'd2': {
                'blyp'        : {'s6': 1.2,  'sr6': 1.1,   's8': 0.0,   'alpha6': 20.0},  # in psi4
                'bp86'        : {'s6': 1.05, 'sr6': 1.1,   's8': 0.0,   'alpha6': 20.0},  # in psi4
                'b97-d'       : {'s6': 1.25, 'sr6': 1.1,   's8': 0.0,   'alpha6': 20.0},  # in psi4
                'revpbe'      : {'s6': 1.25, 'sr6': 1.1,   's8': 0.0,   'alpha6': 20.0},
                'pbe'         : {'s6': 0.75, 'sr6': 1.1,   's8': 0.0,   'alpha6': 20.0},  # in psi4
                'tpss'        : {'s6': 1.0,  'sr6': 1.1,   's8': 0.0,   'alpha6': 20.0},
                'b3lyp'       : {'s6': 1.05, 'sr6': 1.1,   's8': 0.0,   'alpha6': 20.0},  # in psi4
                'pbe0'        : {'s6': 0.6,  'sr6': 1.1,   's8': 0.0,   'alpha6': 20.0},  # in psi4
                'pw6b95'      : {'s6': 0.5,  'sr6': 1.1,   's8': 0.0,   'alpha6': 20.0},
                'tpss0'       : {'s6': 0.85, 'sr6': 1.1,   's8': 0.0,   'alpha6': 20.0},
                'b2plyp'      : {'s6': 0.55, 'sr6': 1.1,   's8': 0.0,   'alpha6': 20.0},  # in psi4
                'b2gp-plyp'   : {'s6': 0.4,  'sr6': 1.1,   's8': 0.0,   'alpha6': 20.0},
                'dsd-blyp'    : {'s6': 0.41, 'sr6': 1.1,   's8': 0.0,   'alpha6': 60.0},  # in psi4
            },
            'd3zero': {
                'b1b95'       : {'s6': 1.0,  'sr6': 1.613, 's8': 1.868, 'alpha6': 14.0},
                'b2gpplyp'    : {'s6': 0.56, 'sr6': 1.586, 's8': 0.760, 'alpha6': 14.0},
                'b3lyp'       : {'s6': 1.0,  'sr6': 1.261, 's8': 1.703, 'alpha6': 14.0},  # in psi4
                'b97-d'       : {'s6': 1.0,  'sr6': 0.892, 's8': 0.909, 'alpha6': 14.0},  # in psi4
                'bhlyp'       : {'s6': 1.0,  'sr6': 1.370, 's8': 1.442, 'alpha6': 14.0},
                'blyp'        : {'s6': 1.0,  'sr6': 1.094, 's8': 1.682, 'alpha6': 14.0},
                'bp86'        : {'s6': 1.0,  'sr6': 1.139, 's8': 1.683, 'alpha6': 14.0},  # in psi4
                'bpbe'        : {'s6': 1.0,  'sr6': 1.087, 's8': 2.033, 'alpha6': 14.0},
                'mpwlyp'      : {'s6': 1.0,  'sr6': 1.239, 's8': 1.098, 'alpha6': 14.0},
                'pbe'         : {'s6': 1.0,  'sr6': 1.217, 's8': 0.722, 'alpha6': 14.0},  # in psi4
                'pbe0'        : {'s6': 1.0,  'sr6': 1.287, 's8': 0.928, 'alpha6': 14.0},  # in psi4
                'pw6b95'      : {'s6': 1.0,  'sr6': 1.532, 's8': 0.862, 'alpha6': 14.0},
                'pwb6k'       : {'s6': 1.0,  'sr6': 1.660, 's8': 0.550, 'alpha6': 14.0},
                'revpbe'      : {'s6': 1.0,  'sr6': 0.923, 's8': 1.010, 'alpha6': 14.0},
                'tpss'        : {'s6': 1.0,  'sr6': 1.166, 's8': 1.105, 'alpha6': 14.0},
                'tpss0'       : {'s6': 1.0,  'sr6': 1.252, 's8': 1.242, 'alpha6': 14.0},
                'tpssh'       : {'s6': 1.0,  'sr6': 1.223, 's8': 1.219, 'alpha6': 14.0},
                'bop'         : {'s6': 1.0,  'sr6': 0.929, 's8': 1.975, 'alpha6': 14.0},
                'mpw1b95'     : {'s6': 1.0,  'sr6': 1.605, 's8': 1.118, 'alpha6': 14.0},
                'mpwb1k'      : {'s6': 1.0,  'sr6': 1.671, 's8': 1.061, 'alpha6': 14.0},
                'olyp'        : {'s6': 1.0,  'sr6': 0.806, 's8': 1.764, 'alpha6': 14.0},
                'opbe'        : {'s6': 1.0,  'sr6': 0.837, 's8': 2.055, 'alpha6': 14.0},
                'otpss'       : {'s6': 1.0,  'sr6': 1.128, 's8': 1.494, 'alpha6': 14.0},
                'pbe38'       : {'s6': 1.0,  'sr6': 1.333, 's8': 0.998, 'alpha6': 14.0},
                'pbesol'      : {'s6': 1.0,  'sr6': 1.345, 's8': 0.612, 'alpha6': 14.0},
                'revssb'      : {'s6': 1.0,  'sr6': 1.221, 's8': 0.560, 'alpha6': 14.0},
                'ssb'         : {'s6': 1.0,  'sr6': 1.215, 's8': 0.663, 'alpha6': 14.0},
                'b3pw91'      : {'s6': 1.0,  'sr6': 1.176, 's8': 1.775, 'alpha6': 14.0},
                'bmk'         : {'s6': 1.0,  'sr6': 1.931, 's8': 2.168, 'alpha6': 14.0},
                'camb3lyp'    : {'s6': 1.0,  'sr6': 1.378, 's8': 1.217, 'alpha6': 14.0},
                'lcwpbe'      : {'s6': 1.0,  'sr6': 1.355, 's8': 1.279, 'alpha6': 14.0},
                'm05-2x'      : {'s6': 1.0,  'sr6': 1.417, 's8': 0.00 , 'alpha6': 14.0},  # in psi4
                'm05'         : {'s6': 1.0,  'sr6': 1.373, 's8': 0.595, 'alpha6': 14.0},  # in psi4
                'm062x'       : {'s6': 1.0,  'sr6': 1.619, 's8': 0.00 , 'alpha6': 14.0},
                'm06hf'       : {'s6': 1.0,  'sr6': 1.446, 's8': 0.00 , 'alpha6': 14.0},
                'm06l'        : {'s6': 1.0,  'sr6': 1.581, 's8': 0.00 , 'alpha6': 14.0},
                'm06'         : {'s6': 1.0,  'sr6': 1.325, 's8': 0.00 , 'alpha6': 14.0},
                'hcth120'     : {'s6': 1.0,  'sr6': 1.221, 's8': 1.206, 'alpha6': 14.0},  # in psi4
                'b2plyp'      : {'s6': 0.64, 'sr6': 1.427, 's8': 1.022, 'alpha6': 14.0},  # in psi4
                'dsd-blyp'    : {'s6': 0.50, 'sr6': 1.569, 's8': 0.705, 'alpha6': 14.0},  # in psi4
                'ptpss'       : {'s6': 0.75, 'sr6': 1.541, 's8': 0.879, 'alpha6': 14.0},
                'pwpb95'      : {'s6': 0.82, 'sr6': 1.557, 's8': 0.705, 'alpha6': 14.0},
                'revpbe0'     : {'s6': 1.0,  'sr6': 0.949, 's8': 0.792, 'alpha6': 14.0},
                'revpbe38'    : {'s6': 1.0,  'sr6': 1.021, 's8': 0.862, 'alpha6': 14.0},
                'rpw86pbe'    : {'s6': 1.0,  'sr6': 1.224, 's8': 0.901, 'alpha6': 14.0},
            },
            'd3bj': {
                'b1b95'       : {'s6': 1.000, 'a1':  0.2092, 's8':  1.4507, 'a2': 5.5545},
                'b2gpplyp'    : {'s6': 0.560, 'a1':  0.0000, 's8':  0.2597, 'a2': 6.3332},
                'b3pw91'      : {'s6': 1.000, 'a1':  0.4312, 's8':  2.8524, 'a2': 4.4693},
                'bhlyp'       : {'s6': 1.000, 'a1':  0.2793, 's8':  1.0354, 'a2': 4.9615},
                'bmk'         : {'s6': 1.000, 'a1':  0.1940, 's8':  2.0860, 'a2': 5.9197},
                'bop'         : {'s6': 1.000, 'a1':  0.4870, 's8':  3.295,  'a2': 3.5043},
                'bpbe'        : {'s6': 1.000, 'a1':  0.4567, 's8':  4.0728, 'a2': 4.3908},
                'camb3lyp'    : {'s6': 1.000, 'a1':  0.3708, 's8':  2.0674, 'a2': 5.4743},
                'lcwpbe'      : {'s6': 1.000, 'a1':  0.3919, 's8':  1.8541, 'a2': 5.0897},
                'mpw1b95'     : {'s6': 1.000, 'a1':  0.1955, 's8':  1.0508, 'a2': 6.4177},
                'mpwb1k'      : {'s6': 1.000, 'a1':  0.1474, 's8':  0.9499, 'a2': 6.6223},
                'mpwlyp'      : {'s6': 1.000, 'a1':  0.4831, 's8':  2.0077, 'a2': 4.5323},
                'olyp'        : {'s6': 1.000, 'a1':  0.5299, 's8':  2.6205, 'a2': 2.8065},
                'opbe'        : {'s6': 1.000, 'a1':  0.5512, 's8':  3.3816, 'a2': 2.9444},
                'otpss'       : {'s6': 1.000, 'a1':  0.4634, 's8':  2.7495, 'a2': 4.3153},
                'pbe38'       : {'s6': 1.000, 'a1':  0.3995, 's8':  1.4623, 'a2': 5.1405},
                'pbesol'      : {'s6': 1.000, 'a1':  0.4466, 's8':  2.9491, 'a2': 6.1742},
                'ptpss'       : {'s6': 0.750, 'a1':  0.000,  's8':  0.2804, 'a2': 6.5745},
                'pwb6k'       : {'s6': 1.000, 'a1':  0.1805, 's8':  0.9383, 'a2': 7.7627},
                'revssb'      : {'s6': 1.000, 'a1':  0.4720, 's8':  0.4389, 'a2': 4.0986},
                'ssb'         : {'s6': 1.000, 'a1': -0.0952, 's8': -0.1744, 'a2': 5.2170},
                'tpssh'       : {'s6': 1.000, 'a1':  0.0000, 's8':  0.4243, 'a2': 5.5253},
                'hcth120'     : {'s6': 1.000, 'a1':  0.3563, 's8':  1.0821, 'a2': 4.3359},  # in psi4
                'b2plyp'      : {'s6': 0.640, 'a1':  0.3065, 's8':  0.9147, 'a2': 5.0570},  # in psi4
                'b3lyp'       : {'s6': 1.000, 'a1':  0.3981, 's8':  1.9889, 'a2': 4.4211},  # in psi4
                'b97-d'       : {'s6': 1.000, 'a1':  0.5545, 's8':  2.2609, 'a2': 3.2297},  # in psi4
                'blyp'        : {'s6': 1.000, 'a1':  0.4298, 's8':  2.6996, 'a2': 4.2359},  # in psi4
                'bp86'        : {'s6': 1.000, 'a1':  0.3946, 's8':  3.2822, 'a2': 4.8516},  # in psi4
                'dsd-blyp'    : {'s6': 0.500, 'a1':  0.000,  's8':  0.2130, 'a2': 6.0519},  # in psi4
                'pbe0'        : {'s6': 1.000, 'a1':  0.4145, 's8':  1.2177, 'a2': 4.8593},  # in psi4
                'pbe'         : {'s6': 1.000, 'a1':  0.4289, 's8':  0.7875, 'a2': 4.4407},  # in psi4
                'pw6b95'      : {'s6': 1.000, 'a1':  0.2076, 's8':  0.7257, 'a2': 6.3750},
                'pwpb95'      : {'s6': 0.820, 'a1':  0.0000, 's8':  0.2904, 'a2': 7.3141},
                'revpbe0'     : {'s6': 1.000, 'a1':  0.4679, 's8':  1.7588, 'a2': 3.7619},
                'revpbe38'    : {'s6': 1.000, 'a1':  0.4309, 's8':  1.4760, 'a2': 3.9446},
                'revpbe'      : {'s6': 1.000, 'a1':  0.5238, 's8':  2.3550, 'a2': 3.5016},
                'rpw86pbe'    : {'s6': 1.000, 'a1':  0.4613, 's8':  1.3845, 'a2': 4.5062},
                'tpss0'       : {'s6': 1.000, 'a1':  0.3768, 's8':  1.2576, 'a2': 4.5865},
                'tpss'        : {'s6': 1.000, 'a1':  0.4535, 's8':  1.9435, 'a2': 4.4752},
            }
        }

        # Validate arguments
        dashlvl = dashlvl.lower()
        if dashlvl not in dashcoeff.keys():
            raise ValidationError("""-D correction level %s is not available. Choose among %s.""" % (dashlvl, dashcoeff.keys()))

        if func is None:
            if dashparam is None:
                # defunct case
                raise ValidationError("""Parameters for -D correction missing. Provide a func or a dashparam kwarg.""")
            else:
                # case where all param read from dashparam dict (which must have all correct keys)
                func = 'custom'
                dashcoeff[dashlvl][func] = {}
                dashparam = dict((k.lower(), v) for k, v in dashparam.iteritems())
                for key in dashcoeff[dashlvl]['b3lyp'].keys():
                    if key in dashparam.keys():
                        dashcoeff[dashlvl][func][key] = dashparam[key]
                    else:
                        raise ValidationError("""Parameter %s is missing from dashparam dict %s.""" % (key, dashparam))
        else:
            func = func.lower()
            if func not in dashcoeff[dashlvl].keys():
                raise ValidationError("""Functional %s is not available for -D level %s.""" % (func, dashlvl))
            if dashparam is None:
                # (normal) case where all param taken from dashcoeff above
                pass
            else:
                # case where items in dashparam dict can override param taken from dashcoeff above
                dashparam = dict((k.lower(), v) for k, v in dashparam.iteritems())
                for key in dashcoeff[dashlvl]['b3lyp'].keys():
                    if key in dashparam.keys():
                        dashcoeff[dashlvl][func][key] = dashparam[key]

        # Move ~/.dftd3par.<hostname> out of the way so it won't interfere
        defaultfile = os.path.expanduser('~') + '/.dftd3par.' + socket.gethostname()
        defmoved = False
        if os.path.isfile(defaultfile):
            os.rename(defaultfile, defaultfile + '_hide')
            defmoved = True

        # Setup unique scratch directory and move in
        current_directory = os.getcwd()
        dftd3_tmpdir = 'dftd3_' + str(random.randint(0, 99999))
        if os.path.exists(dftd3_tmpdir) is False:
            os.mkdir(dftd3_tmpdir)
        os.chdir(dftd3_tmpdir)

        # Write dftd3_parameters file that governs dispersion calc
        paramfile = './dftd3_parameters'
        pfile = open(paramfile, 'w')

        if dashlvl == 'd2':
            # d2:      s6 sr6 s8 a2=None alpha6 version=2
            pfile.write('%12.6f %12.6f %12.6f %12.6f %12.6f %6d\n' %
                (dashcoeff[dashlvl][func]['s6'], dashcoeff[dashlvl][func]['sr6'], dashcoeff[dashlvl][func]['s8'],
                0.0, dashcoeff[dashlvl][func]['alpha6'], 2))
        elif dashlvl == 'd3zero':
            # d3zero:  s6 sr6 s8 a2=None alpha6 version=3
            pfile.write('%12.6f %12.6f %12.6f %12.6f %12.6f %6d\n' %
                (dashcoeff[dashlvl][func]['s6'], dashcoeff[dashlvl][func]['sr6'], dashcoeff[dashlvl][func]['s8'],
                1.0, dashcoeff[dashlvl][func]['alpha6'], 3))
        elif dashlvl == 'd3bj':
            # d3bj:    s6 a1 s8 a2 alpha6=None version=4
            pfile.write('%12.6f %12.6f %12.6f %12.6f %12.6f %6d\n' %
                (dashcoeff[dashlvl][func]['s6'], dashcoeff[dashlvl][func]['a1'], dashcoeff[dashlvl][func]['s8'],
                dashcoeff[dashlvl][func]['a2'], 0.0, 4))
        pfile.close()

        # Write dftd3_geometry file that supplies geometry to dispersion calc
        geomfile = './dftd3_geometry.xyz'
        gfile = open(geomfile, 'w')
        gfile.write(self.save_string_xyz())
        gfile.close()

        # Call dftd3 program
        try:
            dashout = subprocess.Popen(['dftd3', geomfile, '-grad'], stdout=subprocess.PIPE)
        except OSError:
            raise ValidationError('Program dftd3 not found in path.')
        out, err = dashout.communicate()
        if verbosity >= 3:
            print out

        # Parse output (could go further and break into E6, E8, E10 and Cn coeff)
        success = False
        for line in out.splitlines():
            if re.match(' Edisp /kcal,au', line):
                sline = line.split()
                dashd = float(sline[3])
            if re.match(' normal termination of dftd3', line):
                success = True

        if not success:
            raise ValidationError('Program dftd3 did not complete successfully.')

        # Parse grad output
        derivfile = './dftd3_gradient'
        dfile = open(derivfile, 'r')
        dashdderiv = []
        for at in dfile.readlines():
            dashdderiv.append([float(x.replace('D', 'E')) for x in at.split()])
        dfile.close()
        if len(dashdderiv) != self.natom():
            raise ValidationError('Program dftd3 gradient file has %d atoms- %d expected.' % \
                (len(dashdderiv), self.natom()))

        # Clean up files and remove scratch directory
#        os.unlink(paramfile)
#        os.unlink(geomfile)
#        os.unlink(derivfile)
        if defmoved is True:
            os.rename(defaultfile + '_hide', defaultfile)

        os.chdir('..')
#        try:
#            shutil.rmtree(dftd3_tmpdir)
#        except OSError as e:
#            ValidationError('Unable to remove dftd3 temporary directory %s' % e, file=sys.stderr)
        os.chdir(current_directory)

        # return -D & d(-D)/dx
        return dashd, dashdderiv

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
