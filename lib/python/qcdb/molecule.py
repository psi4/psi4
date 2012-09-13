#import os
#import re
#import math
#import copy
#from periodictable import *
#from physconst import *
#from vecutil import *
#from exceptions import *
#from coordentry import *
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

    @classmethod
    def init_with_xyz(cls, xyzfilename, no_com=False, no_reorient=False):
        """Pull information from an XYZ file. No fragment info detected.
        Charge, multiplicity, tagline pulled from second line if available.

        >>> H2O = qcdb.Molecule.init_with_xyz('h2o.xyz')

        """
        instance = cls()
        instance.lock_frame = False
        instance.PYmove_to_com = not no_com
        instance.PYfix_orientation = no_reorient

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

    def save_string_for_psi4(self):
        """Returns a string of Molecule formatted for psi4.
        Includes fragments and reorienting, if specified.

        >>> print H2OH2O.save_string_for_psi4()
        6
        0 1
        O         -1.55100700      -0.11452000       0.00000000
        H         -1.93425900       0.76250300       0.00000000
        H         -0.59967700       0.04071200       0.00000000
        --
        0 1
        @X         0.00000000       0.00000000       0.00000000
        O          1.35062500       0.11146900       0.00000000
        H          1.68039800      -0.37374100      -0.75856100
        H          1.68039800      -0.37374100       0.75856100
        units Angstrom

        """
        Nfr = 0
        text = ""
        for fr in range(self.nfragments()):
            if self.fragment_types[fr] == 'Absent':
                continue
            if Nfr != 0:
                text += """--\n"""
            Nfr += 1
            text += """%d %d\n""" % (self.fragment_charges[fr], self.fragment_multiplicities[fr])
            for at in range(self.fragments[fr][0], self.fragments[fr][1] + 1):
                geom = self.full_atoms[at].compute()
                text += """%-3s  %16.8f %16.8f %16.8f\n""" % \
                    (("" if self.fZ(at) else "@") + self.full_atoms[at].symbol(), \
                    geom[0], geom[1], geom[2])
        text += """units %s\n""" % (self.units().lower())
        return text

    def format_string_for_qchem(self):
        pass

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
            'H':  1.001 / 1.5,
            'HE': 1.012 / 1.5,
            'LI': 0.825 / 1.5,
            'BE': 1.408 / 1.5,
            'B':  1.485 / 1.5,
            'C':  1.452 / 1.5,
            'N':  1.397 / 1.5,
            'O':  1.342 / 1.5,
            'F':  1.287 / 1.5,
            'NE': 1.243 / 1.5,
            'NA': 1.144 / 1.5,
            'MG': 1.364 / 1.5,
            'AL': 1.639 / 1.5,
            'SI': 1.716 / 1.5,
            'P':  1.705 / 1.5,
            'S':  1.683 / 1.5,
            'CL': 1.639 / 1.5,
            'AR': 1.595 / 1.5}

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
