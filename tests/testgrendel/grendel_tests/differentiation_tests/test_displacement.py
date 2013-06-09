
import unittest
import sys
import os

# Add the directory containing 'grendel' to the sys.path
sys.path.append(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir, os.pardir))

# Add the directory containing the 'grendel_tests' package to sys.path
sys.path.append(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir))

from grendel_tests import long_test, profile
from grendel import *
from grendel.differentiation.displacement import Displacement

class DisplacementTest(unittest.TestCase):

    def test_displace_bond_length(self):
        mol = Molecule.from_z_matrix("""
                H
                O 1 1.0
            """, create_representation=True
        )
        rep = mol.internal_representation

        d = Displacement.from_increments(Vector(1), rep, Vector(0.1))
        self.assertAlmostEqual(d.displaced_representation[0].value, 1.1)

    def test_displace_bond_angle(self):
        deg_to_def = Degrees.to(AngularUnit.default)

        mol = Molecule.from_z_matrix("""
                H
                O 1 1.0
                H 2 1.0 1 90
            """, create_representation=True
        )
        rep = mol.internal_representation

        d = Displacement.from_increments(Vector(1, 0, 0), rep, Vector([0.1]*3))
        self.assertAlmostEqual(d.displaced_representation[0].value, 1.1)
        self.assertAlmostEqual(d.displaced_representation[2].value, 90.*deg_to_def)

    def test_displace_torsion(self):
        mol = Molecule.from_z_matrix("""
                H
                O 1 1.0
                O 2 1.0 1 90
                H 3 1.0 2 90 1 90
            """, create_representation=True
        )
        rep = mol.internal_representation
        d = Displacement.from_increments(Vector(0,0,0,0,0,1), rep, Vector(0,0,0,0,0,0.1))
        deg_to_def = Degrees.to(AngularUnit.default)
        self.assertAlmostEqual(d.displaced_representation[0].value, 1.0)
        self.assertAlmostEqual(d.displaced_representation[5].value, 90.1 * deg_to_def)

    def test_displace_all_zmat(self):
        mol = Molecule.from_z_matrix("""
                H
                O 1 1.0
                O 2 1.0 1 90
                H 3 1.0 2 90 1 -90
            """, create_representation=True
        )
        rep = mol.internal_representation
        d = Displacement.from_increments(Vector(1,-1,1,-1,1,-1), rep, Vector([0.1]*6))
        deg_to_def = Degrees.to(AngularUnit.default)
        self.assertAlmostEqual(d.displaced_representation[0].value, 1.1)
        self.assertAlmostEqual(d.displaced_representation[1].value, 0.9)
        self.assertAlmostEqual(d.displaced_representation[2].value, 90.1 * deg_to_def)
        self.assertAlmostEqual(d.displaced_representation[3].value, 0.9)
        self.assertAlmostEqual(d.displaced_representation[4].value, 90.1 * deg_to_def)
        self.assertAlmostEqual(d.displaced_representation[5].value, -90.1 * deg_to_def)

    def test_displace_all_hooh(self):
        mol = Molecule.from_z_matrix("""
                H
                O 1 0.963242
                O 2 1.449863 1 100.120071
                H 3 0.963242 2 100.120071 1 -67.344079
            """, create_representation=True
        )
        rep = mol.internal_representation
        deg_to_def = Degrees.to(AngularUnit.default)
        self.assertAlmostEqual(rep[5].value, -67.344079 * deg_to_def)
        d = Displacement.from_increments(Vector(0,-1,0,0,0,0), rep, Vector([0.1]*6))
        self.assertAlmostEqual(0.963242, d.displaced_representation[0].value)
        self.assertAlmostEqual(1.349863, d.displaced_representation[1].value)
        self.assertAlmostEqual(100.120071 * deg_to_def, d.displaced_representation[2].value_with_units.in_units(AngularUnit.default))
        self.assertAlmostEqual(0.963242, d.displaced_representation[3].value)
        self.assertAlmostEqual(100.120071 * deg_to_def, d.displaced_representation[4].value_with_units.in_units(AngularUnit.default))
        self.assertAlmostEqual(-67.344079 * deg_to_def, d.displaced_representation[5].value_with_units.in_units(AngularUnit.default))

    def test_displace_hooh_cart(self):
        mol = Molecule("""
            H         -0.8026739226       -0.8820752737        0.4535861764
            O         -0.0110284705       -0.7248476524       -0.0721631764
            O          0.0110284705        0.7248476524       -0.0721631764
            H          0.8026739226        0.8820752737        0.4535861764
        """
        )
        mol.convert_units(Bohr)
        mol.internal_representation = rep = InternalRepresentation(mol, """
            bond 2 1
            bond 3 2
            bend 3 2 1
            bond 4 3
            bend 4 3 2
            tors 4 3 2 1
        """
        )
        d = Displacement.from_increments((0,-1,0,0,0,0), rep,
            Vector([0.1*Angstrom, 0.1*Angstrom, 0.2*Degrees, 0.1*Angstrom, 0.2*Degrees, 0.2*Degrees])
        )
        deg_to_def = Degrees.to(AngularUnit.default)
        self.assertAlmostEqual(
            rep[1].in_units(Angstrom) - 0.1*Angstrom,
            d.displaced_representation[1].in_units(Angstrom)
        )

