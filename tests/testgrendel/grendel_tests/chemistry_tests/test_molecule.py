from copy import copy
import unittest
import sys
import os

from numpy.testing import assert_array_almost_equal

# Add the directory containing 'grendel' to the sys.path
sys.path.append(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir, os.pardir))

# Add the directory containing the 'grendel_tests' package to sys.path
sys.path.append(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir))

from grendel import *
import grendel.chemistry
from grendel.chemistry.atom import Atom
from grendel.chemistry.molecule import InvalidXYZFormatError, InvalidZMatrixException, Molecule
from grendel.chemistry import SampleMolecules
from grendel_tests import skip

class moleculeTest(unittest.TestCase):

    ##############
    # Assertions #
    ##############

    def assertAlmostEqualMolecules(self, act, exp, msg=None, places=7):
        a_atoms = copy(exp.atoms)
        b_atoms = copy(act.atoms)
        for atom in a_atoms:
            atom.position = Vector([round(num, places) for num in atom.position])
        for atom in b_atoms:
            atom.position = Vector([round(num, places) for num in atom.position])
        round_a = Molecule(a_atoms)
        round_b = Molecule(b_atoms)
        if not round_a == round_b:
            seq_a = []
            for a in round_a: seq_a.append((a.symbol, a.position))
            seq_b = []
            for b in round_b: seq_b.append((b.symbol, b.position))
            seq_msg=None
            try:
                self.assertSequenceEqual(seq_a, seq_b)
            except AssertionError as err:
                seq_msg = err.message
            if not seq_msg:
                raise RuntimeError("Molecule almost equal assertion inconclusive: Rounded sequences are equal but rounded molecules are not.")
            if msg is None:

                msg = "\nExpected:\n{exp}\nGot:\n{act}\n\nAs sequences:\n".format(
                    exp=exp, act=act, ndigits=places
                ) + seq_msg
            raise unittest.TestCase.failureException(msg)

    ##################
    # Setup/Teardown #
    ##################

    def setUp(self):
        grendel.chemistry.init_sample_molecules()
        self.addTypeEqualityFunc(Molecule, self.assertAlmostEqualMolecules)

    def tearDown(self):
        pass

    #########
    # Tests #
    #########

    def test_init(self):
        self.assertEqual(Molecule("H 1 -1 0"), Molecule([Atom('H', [1, -1, 0])]))

    def test_rotational_constants(self):
        mol = SampleMolecules['Benzene']
        self.assertAlmostEqual(mol.A_e.value, 0.19131066301)
        self.assertAlmostEqual(mol.B_e.value, 0.191310586465)
        self.assertAlmostEqual(mol.C_e.value, 0.0956553123686)

    def test_errors_1(self):
        # Stupid things to pass to _init_xyz_string
        self.assertRaises(InvalidXYZFormatError, Molecule,
            """
            4

            H 1.00 1.00 7.5
            C 3.2 1.0 1.0
            O 2.2 2.2 2.2

            """
        )

    def test_errors_2(self):
        self.assertRaises(InvalidXYZFormatError, Molecule,
            """
            3

            H 1.00 1.00 7.5
            C 3.2a 1.0 1.0
            O 2.2 2.2 2.2

            """
        )

    def test_errors_3(self):
        # Invalid z-matrices
        self.assertRaises(InvalidZMatrixException,
            Molecule.from_z_matrix,
            """
            O 1 1.0
            """
        )

    def test_errors_4(self):
        self.assertRaises(InvalidZMatrixException,
            Molecule.from_z_matrix,
            """
            O
            H 1
            """
        )

    def test_errors_5(self):
        self.assertRaises(InvalidZMatrixException,
            Molecule.from_z_matrix,
            """
            O
            H 1 1.0 2
            """
        )

    def test_errors_6(self):
        self.assertRaisesRegexp(IndexError, r'z-matrix index',
            Molecule.from_z_matrix,
            """
            O
            H 2 1.0
            """
        )

    def test_errors_7(self):
        self.assertRaises(InvalidZMatrixException,
            Molecule.from_z_matrix,
            """
            O
            H 1 1.0
            H 2 1.0 1
            """
        )

    def test_errors_8(self):
        self.assertRaises(InvalidZMatrixException,
            Molecule.from_z_matrix,
            """
            O
            H 1 1.0
            H 2 1.0 1 3.7 2
            """
        )

    def test_errors_9(self):
        self.assertRaises(InvalidZMatrixException,
            Molecule.from_z_matrix,
            """
            O
            H 1 1.0
            H 2 1.0 1 3.7
            O 1 1.1 2 1.2 3
            """
        )

    def test_errors_10(self):
        self.assertRaises(NotImplementedError,
            Molecule.from_z_matrix,
            """
            O
            H 1 1.0
            H 2 1.0 1 3.7
            O 1 1.1 2 1.2 3 3.5 2
            """
        )

    def test_errors_11(self):
        with self.assertRaisesRegexp(ValueError, r'Ambiguous atom identifier in z-matrix'):
            Molecule.from_z_matrix(
                """
                O
                H O 1.0
                H O -1.0 H 37
                O2 H 1.5 O 150 H 25
                """
            )

    def test_errors_12(self):
        self.assertRaisesRegexp(ValueError, r'Unknown atom identifier in z-matrix',
            Molecule.from_z_matrix,
            """
            O
            H1 O 1.0
            H2 H1 1.0 O1 3.7
            """
        )

    def test_errors_13(self):
        self.assertRaisesRegexp(ValueError, r'Invalid atom identifier in z-matrix: hello',
            Molecule.from_z_matrix,
            """
            O
            H1 hello 1.0
            H2 H1 1.0 O1 3.7
            """
        )

    def test_errors_14(self):
        with self.assertRaises(InvalidXYZFormatError):
            Molecule(
                """
                    H 1.00 1.00 7.5
                    C 3.2a 1.0 1.0
                    O 2.2 2.2 2.2
                """
            )

    def test_errors_15(self):
        # Stupid things to do to _init_atom_list
        self.assertRaises(TypeError, Molecule, [Atom("O", 1.0, 1.0, 1.0), [7.5, 5.0, 3.5]])

    def test_errors_16(self):
        # Stupid things to do to _init_atoms_matrix
        self.assertRaises(ValueError, Molecule, ["O","H","H"], Matrix((1,2,3),(4,5,6)))
        self.assertRaises(ValueError, Molecule, ["O","H","H"], Matrix(1,2,3,4,5,6,7,8,9))
        self.assertRaises(ValueError, Molecule, ["O","H","H"], Matrix((1,2,3,4),(4,5,6,7),(8,9,10,11)))

    def test_errors_17(self):
        # Z-Mat with the wrong number or type of arguments
        self.assertRaisesRegexp(TypeError, r'takes at least ', Molecule.from_z_matrix)
        self.assertRaises(TypeError, Molecule.from_z_matrix, 42)
        self.assertRaises(TypeError, Molecule.from_z_matrix, [42,42])
        self.assertRaisesRegexp(ValueError, r'Invalid atom label in z-matrix', Molecule.from_z_matrix, [[42],[42]])

    def test_errors_18(self):
        with self.assertRaisesRegexp(IndexError, r"z-matrix index '0' is out of range"):
            Molecule.from_z_matrix('''
                O
                H 1 0.93
                C 2 1.23 1 104.5
                B 2 1.35 1  89.5 0 84.3
            ''', create_representation=True)

    def test_errors_19(self):
        with self.assertRaisesRegexp(InvalidZMatrixException, r'atom indices in z-matrix line'):
            Molecule.from_z_matrix('''
                O
                H 1 0.93
                C 2 1.23 1 104.5
                B 2 1.35 1  89.5 1 84.3
            ''', create_representation=True)

    @skip("Not yet implemented")
    def test_future_errors(self):
        with self.assertRaisesRegexp(ValueError, r'Ambiguous atom identifier in z-matrix'):
            Molecule.from_z_matrix(
                """
                O
                H O 1.0
                H H 1.0 O 3.7
                """
            )

    def test_misc(self):
        # getitem
        self.assertEqual(Molecule.from_z_matrix("O")[0].symbol, "O")

        # contains
        mol = Molecule.from_z_matrix("O")
        self.assertIn(mol[0], mol)
        self.assertNotIn(Atom("H", 0, 0, 0), mol)

        # less than
        self.assertLess(mol, Molecule("O 0.0 0.0 1.0"))

    def test_zmat(self):
        self.assertEqual(
            Molecule([Atom("O", [0.0, 0.0, 0.0])]),
            Molecule.from_z_matrix("""
                O
            """
            )
        )
        self.assertEqual(
            Molecule([
                Atom("O", [0.0, 0.0, 0.0]),
                Atom("O", [0.0, 0.0, 1.0])
            ]),
            Molecule.from_z_matrix("""
                O1
                O2 O1 1.0
            """
            )
        )
        self.assertEqual(
            Molecule([
                Atom("O", [0.0,  0.0, 0.0]),
                Atom("H", [0.0,  0.0, 1.0]),
                Atom("H", [0.0, -1.0, 0.0])
            ]),
            Molecule.from_z_matrix("""
                O
                H1 O 1.0
                H2 O 1.0 H1 90
            """
            )
        )

    def test_displace(self):
        mol = Molecule.from_z_matrix('H')
        mol.displace(Vector(1., -1., 0.))
        self.assertEqual(mol, Molecule("H 1 -1 0"))

        mol = Molecule.from_z_matrix('H\nH 1 0.9')
        mol.displace(Vector(1., 0., 0., 1., 0., 0.))
        self.assertEqual(mol, Molecule("H 1 0 0.0\nH 1 0 0.9 "))

        with self.assertRaisesRegexp(ValueError, 'dimension mismatch'):
            mol.displace(Vector(1., 0., 0.))

        with self.assertRaises(TypeError):
            mol.displace(Matrix(1., 0., 0.))

    def test_create_representation_bond(self):
        with self.assertRaisesRegexp(ValueError, r'atomics.*only representable as.*Cannot create'):
            Molecule.from_z_matrix('O', create_representation=True)

        mol = Molecule.from_z_matrix('H\nH 1 0.9', create_representation=True)
        self.assertEqual(len(mol.internal_representations), 1)
        self.assertEqual(len(mol.internal_representations[0]), 1)

    def test_create_representation_ang(self):
        mol = Molecule.from_z_matrix("""
                O
                H1 O 1.0
                H2 O 1.0 H1 90
            """, create_representation=True
            )
        self.assertEqual(len(mol.internal_representations), 1)
        self.assertEqual(len(mol.internal_representations[0]), 3)
        self.assertEqual(
            mol.internal_representations[0].coords[2].value_with_units.in_units(AngularUnit.default),
            (90.*Degrees).in_units(AngularUnit.default)
        )

    def test_create_representation_tors(self):
        mol = Molecule.from_z_matrix("""
                H
                O 1 1.0
                O 2 1.0 1 90
                H 3 1.0 2 90 1 90
            """, create_representation=True
        )
        self.assertEqual(len(mol.internal_representations), 1)
        self.assertEqual(len(mol.internal_representations[0]), 6)
        self.assertEqual(mol.internal_representations[0].coords[2].value_with_units.in_units(AngularUnit.default),
            (90.*Degrees).in_units(AngularUnit.default))
        self.assertEqual(mol.internal_representations[0].coords[5].value_with_units.in_units(AngularUnit.default),
            (90.*Degrees).in_units(AngularUnit.default))

    def test_create_representation_tors_2(self):
        mol = Molecule.from_z_matrix("""
                H
                O 1 1.0
                O 2 1.0 1 90
                H 3 1.0 2 90 1 -90
            """, create_representation=True
        )
        self.assertEqual(mol.internal_representations[0].coords[5].value_with_units.in_units(AngularUnit.default),
            (-90.*Degrees).in_units(AngularUnit.default))

    def test_create_rep_hooh(self):
        mol = Molecule.from_z_matrix("""
                H
                O 1 0.963242
                O 2 1.449863 1 100.120071
                H 3 0.963242 2 100.120071 1 -67.344079
            """, create_representation=True
        )
        rep = mol.internal_representation
        deg_to_def = Degrees.to(AngularUnit.default)
        self.assertAlmostEqual(
            rep[5].value_with_units.in_units(AngularUnit.default),
            -67.344079 * deg_to_def)

    def test_reorient_1(self):
        m = Molecule('H 0 0 0\nH 1 0 0')
        self.assertEqual(m.reoriented('II'), Molecule('H 0 -0.5 0\nH 0  0.5 0'))

    def test_reorient_2(self):
        m = Molecule('H 0 0 0\nH 0 1 0')
        self.assertEqual(m.reoriented('II'), Molecule('H 0 -0.5 0\nH 0 0.5 0'))

    #--------------------------------------------------------------------------------#

    def test_rotate_1(self):
        m = SampleMolecules['water']
        axis = Vector([1,2,3]).normalized()
        m2 = m.rotated(angle=180.*Degrees, axis=axis)
        m3 = m2.rotated(angle=180.*Degrees, axis=axis)
        assert_array_almost_equal(m.position, m3.position)

    def test_rotate_2(self):
        m = SampleMolecules['water']
        axis = Vector([1,2,3]).normalized()
        m2 = m.rotated(angle=120.*Degrees, axis=axis)
        m3 = m2.rotated(angle=120.*Degrees, axis=axis)
        m4 = m3.rotated(angle=120.*Degrees, axis=axis)
        assert_array_almost_equal(m.position, m4.position)

    #--------------------------------------------------------------------------------#

    def test_groups_1(self):
        m = Molecule([
                Atom('C', [  0.00000000,  0.00000000,  0.00000000 ] ),
                Atom('C', [  0.00000000,  0.00000000,  4.76992933 ] ),
                Atom('O', [  0.00000000, -1.04316184,  0.61707065 ] ),
                Atom('O', [  0.01905095,  1.04298787,  4.15285868 ] ),
                Atom('C', [ -0.11039651,  1.34908096,  0.68132447 ] ),
                Atom('C', [ -0.13501595, -1.34683982,  4.08860486 ] ),
                Atom('C', [  0.10780157,  0.01502933, -1.51597276 ] ),
                Atom('C', [  0.10750912, -0.01699557,  6.28590209 ] ),
                Atom('H', [ -0.08557151,  1.24276213,  1.76717696 ] ),
                Atom('H', [ -0.10825342, -1.24099210,  3.00275237 ] ),
                Atom('H', [  0.69789248,  2.01081145,  0.34934100 ] ),
                Atom('H', [  0.66105324, -2.02322149,  4.42058833 ] ),
                Atom('H', [ -1.04824273,  1.83250625,  0.38051647 ] ),
                Atom('H', [ -1.08153441, -1.81305690,  4.38941286 ] ),
                Atom('H', [  0.11566492, -1.00528185, -1.90094854 ] ),
                Atom('H', [  0.13400478,  1.00300183,  6.67087787 ] ),
                Atom('H', [ -0.72590461,  0.57377279, -1.95554705 ] ),
                Atom('H', [ -0.73626218, -0.56042012,  6.72547638 ] ),
                Atom('H', [  1.02679855,  0.52908924, -1.82068069 ] ),
                Atom('H', [  1.01696471, -0.54775311,  6.59061002 ] )
            ], units=Angstroms
        )

        groups = m.geometric_subgroups()
        self.assertEqual(len(groups), 2)


