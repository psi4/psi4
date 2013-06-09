from copy import copy
from functools import partial
from itertools import product
import re
import unittest
import sys
import os
import math

import numpy as np
from numpy.testing import assert_array_almost_equal

from numpy.testing.utils import assert_allclose

# Add the directory containing 'grendel' to the sys.path
sys.path.append(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir, os.pardir))

# Add the directory containing the 'grendel_tests' package to sys.path
sys.path.append(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir))

from grendel import *
from grendel_tests import allow_generators, long_generator, generator, long_test, profile, skip, expected_failure
from grendel.chemistry.molecule import Molecule
from grendel.chemistry import init_sample_molecules, SampleMolecules
from grendel.differentiation.finite_difference import FiniteDifferenceDerivative
from grendel.gmath.tensor import Tensor
from grendel.util.iteration import flattened

@allow_generators
class FirstOrderBTensorTest(unittest.TestCase):

    def assertSequenceAlmostEqual(self, seq1, seq2, msg=None, seq_type=None, places=7):
        rseq1 = map(lambda x: round(x, places), seq1)
        rseq2 = map(lambda x: round(x, places), seq2)
        self.assertSequenceEqual(rseq1, rseq2, msg, seq_type)

    def assertHasValidBVector(self, coord, robustness=8, places=9):
        fdiffs = [FiniteDifferenceDerivative(coord, c, robustness=robustness) for c in coord.variables]
        values = [float(f.value) for f in fdiffs]
        self.assertSequenceAlmostEqual(values, coord.b_vector, places=places)

    def test_bond_length(self):
        mol = Molecule.from_z_matrix('H\nO 1 1.0', create_representation=True)
        coord = mol.internal_representation[0]
        self.assertHasValidBVector(coord)

    def test_bond_angle(self):
        mol = Molecule.from_z_matrix('H\nO 1 1.0\nH 2 1.0 1 90', create_representation=True)
        coord = mol.internal_representation[2]
        self.assertHasValidBVector(coord)

    @long_generator
    def bond_length_test_generator(cls):
        npoints = 10
        minval = 0.7
        maxval = 14.0
        def _generated(self, bval):
            mol = Molecule.from_z_matrix('H\nO 1 '+str(bval), create_representation=True)
            coord = mol.internal_representation[0]
            self.assertHasValidBVector(coord)
        diff = maxval - minval
        step = diff / (npoints-1)
        for i in range(npoints):
            currval = minval + (i * step)
            _generated.__name__ = 'test_bond_length_' + re.sub('\.', '_', str(currval))
            yield _generated, currval

    @long_generator
    def bond_angle_test_generator(cls):
        npoints = 30
        minval = 28
        maxval = 174
        def _generated(self, aval):
            mol = Molecule.from_z_matrix('H\nO 1 1.0\nH 2 1.0 1 ' + str(aval), create_representation=True)
            coord = mol.internal_representation[2]
            self.assertHasValidBVector(coord)
        diff = maxval - minval
        step = diff / (npoints-1)
        for i in range(npoints):
            currval = minval + (i * step)
            _generated.__name__ = 'test_bond_angle_' + re.sub('[\. ]', '_', str(currval))
            yield _generated, currval


#--------------------------------------------------------------------------------#
# Helper function that checks the B tensor of a given order by finite difference #
#   of B tensors of order one less.
#--------------------------------------------------------------------------------#

def assert_has_valid_b_tensor(
        coord,
        order,
        robustness=4,
        places=7,
        print_precision=4,
        forward=False,
        highlight_char='*',
        show_vals=False,
        dry_run=False,
        use_coordinate_indices=False,
        print_line_width=120):
    #--------------------------------------------------------------------------------#
    # output helper functions
    def numbered(instr):
        width=len(str(len(instr.splitlines())))
        return ''.join(('#%'+str(width)+'d') % (i+1) + "| " + line + '\n' for i, line in enumerate(instr.splitlines()))
    #----------------------------------------#
    def difference_string(tens, diff_from):
        ret_val = ''
        tstr = numbered(str(tens))
        dstr = numbered(str(diff_from))
        for line, should_be in zip(tstr.splitlines(), dstr.splitlines()):
            ret_val += line + "\n"
            num_re = r'-?\d+\.\d*|inf|nan'
            diff_line = ' ' * print_line_width
            for m, sbm in zip(re.finditer(num_re, line), re.finditer(num_re, should_be)):
                num, num_sb = float(m.group(0)), float(sbm.group(0))
                if abs(num-num_sb) > 10.**(-float(places)):
                    width = m.end() - m.start()
                    if show_vals:
                        corr_val = ('{:'+str(width)+'s}').format(str(num_sb))
                        diff_line = diff_line[:m.start()-1] + highlight_char + corr_val + highlight_char \
                                        + diff_line[m.end()+1:]
                    else:
                        diff_line = diff_line[:m.start()] + highlight_char * width + diff_line[m.end():]

            diff_line = diff_line.rstrip()
            if any(c != ' ' for c in diff_line):
                ret_val += diff_line + "\n"
        return ret_val
    #----------------------------------------#
    np.set_printoptions(
        precision=print_precision,
        edgeitems=300,
        linewidth=print_line_width,
        suppress=True)
    #--------------------------------------------------------------------------------#
    if not use_coordinate_indices:
        btens = coord.parent_representation.b_tensor()
    else:
        btens = coord.get_b_tensor(order)
    if not dry_run:
        fdiff_b_tens = coord.finite_difference_b_tensor(
            order,
            robustness=robustness,
            forward=forward,
            use_parent_indices=not use_coordinate_indices
            )
    natoms = coord.molecule.natoms if not use_coordinate_indices else len(coord.atoms)
    # TODO use the coordinate's internal indexing rather than the parent molecule's
    fdiff_tensor = Tensor(shape=(3*natoms,)*order)
    analytic_tensor = Tensor(shape=(3*natoms,)*order)
    for cartcoords in product(coord.variables, repeat=order):
        if not use_coordinate_indices:
            indices = tuple(c.index for c in cartcoords)
        else:
            indices = coord.internal_indices_for([c.index for c in cartcoords])
        if not use_coordinate_indices:
            analytic_tensor[indices] = btens[(coord,) + cartcoords]
        else:
            analytic_tensor[indices] = btens[indices]
        if not dry_run:
            fdiff_tensor[indices] = fdiff_b_tens[indices]
    #--------------------------------------------------------------------------------#
    #noinspection PyTypeChecker
    assert_array_almost_equal(fdiff_tensor, analytic_tensor, decimal=places, err_msg=""
        "\nArrays not equal to {} decimals"
        "\n"
        "\nfinite difference version:\n{}"
        "\n"
        "\nanalytic version:\n{}"
        "\n".format(
            places,
            numbered(str(fdiff_tensor)),
            difference_string(analytic_tensor, fdiff_tensor)
        )
    )

#--------------------------------------------------------------------------------#
#                              Bond Length Tests                                 #
#--------------------------------------------------------------------------------#

class BTensorBondLengthTest(unittest.TestCase):

    def test_bond_length_1(self):
        mol = Molecule.from_z_matrix('H\nO 1 0.9', create_representation=True)
        mol.recenter()
        # 'randomize' the position in xyz space to avoid miniscule displacement amounts
        mol.rotate(Vector(1,2,1), 32.*Degrees)
        coord = mol.internal_representation[0]
        assert_has_valid_b_tensor(coord, 1)

    def test_bond_length_2(self):
        mol = Molecule.from_z_matrix('H\nO 1 0.9', create_representation=True)
        mol.recenter()
        # 'randomize' the position in xyz space to avoid miniscule displacement amounts
        mol.rotate(Vector(1,2,1), 32.*Degrees)
        coord = mol.internal_representation[0]
        assert_has_valid_b_tensor(coord, 2)

    def test_bond_length_3(self):
        mol = Molecule.from_z_matrix('H\nO 1 0.9', create_representation=True)
        mol.recenter()
        coord = mol.internal_representation[0]
        assert_has_valid_b_tensor(coord, 3)

    def test_bond_length_4(self):
        mol = Molecule.from_z_matrix('H\nO 1 0.9', create_representation=True)
        mol.recenter()
        coord = mol.internal_representation[0]
        assert_has_valid_b_tensor(coord, 4)

    @long_test
    def test_bond_length_5(self):
        mol = Molecule.from_z_matrix('H\nO 1 0.9', create_representation=True)
        mol.recenter()
        coord = mol.internal_representation[0]
        assert_has_valid_b_tensor(coord, 5)

    @long_test
    def test_bond_length_6(self):
        mol = Molecule.from_z_matrix('H\nO 1 0.9', create_representation=True)
        mol.recenter()
        coord = mol.internal_representation[0]
        assert_has_valid_b_tensor(coord, 6)

    @profile("bond_length_4", "bond_length", "b_tensor_short", "b_tensor_all")
    def test_bond_length_4_dry_run(self):
        self.setUp()
        order = 4
        mol = Molecule.from_z_matrix('H\nO 1 0.9', create_representation=True)
        mol.recenter()
        coord = mol.internal_representation[0]
        natoms = mol.natoms
        btens = coord.parent_representation.b_tensor()
        analytic_tensor = Tensor(shape=(3*natoms,)*order)
        for cartcoords in product(coord.variables, repeat=order):
            analytic_tensor[tuple(c.index for c in cartcoords)] = btens[(coord,) + cartcoords]

    @long_test
    @profile("bond_length_5", "bond_length", "b_tensor_short", "b_tensor_all")
    def test_bond_length_5_dry_run(self):
        self.setUp()
        order = 5
        mol = Molecule.from_z_matrix('H\nO 1 0.9', create_representation=True)
        mol.recenter()
        coord = mol.internal_representation[0]
        natoms = mol.natoms
        btens = coord.parent_representation.b_tensor()
        analytic_tensor = Tensor(shape=(3*natoms,)*order)
        for cartcoords in product(coord.variables, repeat=order):
            analytic_tensor[tuple(c.index for c in cartcoords)] = btens[(coord,) + cartcoords]

    @long_test
    @profile("bond_length_6", "bond_length", "b_tensor", "b_tensor_all")
    def test_bond_length_6_dry_run(self):
        self.setUp()
        order = 6
        mol = Molecule.from_z_matrix('H\nO 1 0.9', create_representation=True)
        mol.recenter()
        coord = mol.internal_representation[0]
        natoms = mol.natoms
        btens = coord.parent_representation.b_tensor()
        analytic_tensor = Tensor(shape=(3*natoms,)*order)
        for cartcoords in product(coord.variables, repeat=order):
            analytic_tensor[tuple(c.index for c in cartcoords)] = btens[(coord,) + cartcoords]

    @long_test
    @profile("bond_length_7", "bond_length", "b_tensor", "b_tensor_all")
    def test_bond_length_7_dry_run(self):
        self.setUp()
        order = 7
        mol = Molecule.from_z_matrix('H\nO 1 0.9', create_representation=True)
        mol.recenter()
        coord = mol.internal_representation[0]
        natoms = mol.natoms
        btens = coord.parent_representation.b_tensor()
        analytic_tensor = Tensor(shape=(3*natoms,)*order)
        for cartcoords in product(coord.variables, repeat=order):
            analytic_tensor[tuple(c.index for c in cartcoords)] = btens[(coord,) + cartcoords]

#--------------------------------------------------------------------------------#
#                              Bond Angle Tests                                  #
#--------------------------------------------------------------------------------#

class BTensorBondAngleTest(unittest.TestCase):

    def setUp(self):
        init_sample_molecules()
        self.mol = SampleMolecules['water']
        self.mol.recenter()
        rep = InternalRepresentation(self.mol,
        """
            STRE 1 2
            STRE 1 3
            BEND 2 1 3
        """
        )
        # 'randomize' the position in xyz space to avoid miniscule displacement amounts
        self.mol.rotate(Vector(1,2,1), 32.*Degrees)
        self.coord = rep[2]

    def test_bond_angle_1(self):
        assert_has_valid_b_tensor(self.coord, 1)

    def test_bond_angle_2(self):
        assert_has_valid_b_tensor(self.coord, 2)

    def test_bond_length_offset_3(self):
        assert_has_valid_b_tensor(self.mol.internal_representation[1], 3)

    def test_bond_angle_3(self):
        assert_has_valid_b_tensor(self.coord, 3, show_vals=True, print_precision=6)

    @profile("bond_angle_3", "bond_angle", "b_tensor_all")
    def test_bond_angle_3_dry_run(self):
        self.setUp()
        order = 3
        coord = self.coord
        natoms = self.mol.natoms
        btens = coord.parent_representation.b_tensor()
        analytic_tensor = Tensor(shape=(3*natoms,)*order)
        for cartcoords in product(coord.variables, repeat=order):
            analytic_tensor[tuple(c.index for c in cartcoords)] = btens[(coord,) + cartcoords]

    @long_test
    def test_bond_angle_4(self):
        assert_has_valid_b_tensor(self.coord, 4)

    @long_test
    @profile("bond_angle_4", "bond_angle", "b_tensor_short", "b_tensor_all")
    def test_bond_angle_4_dry_run(self):
        self.setUp()
        order = 4
        coord = self.coord
        natoms = self.mol.natoms
        btens = coord.parent_representation.b_tensor()
        analytic_tensor = Tensor(shape=(3*natoms,)*order)
        for cartcoords in product(coord.variables, repeat=order):
            analytic_tensor[tuple(c.index for c in cartcoords)] = btens[(coord,) + cartcoords]

    @long_test
    def test_bond_angle_5(self):
        assert_has_valid_b_tensor(self.coord, 5)

    @long_test
    @profile("bond_angle_5", "bond_angle", "b_tensor_all")
    def test_bond_angle_5_dry_run(self):
        self.setUp()
        order = 5
        coord = self.coord
        natoms = self.mol.natoms
        btens = coord.parent_representation.b_tensor()
        analytic_tensor = Tensor(shape=(3*natoms,)*order)
        for cartcoords in product(coord.variables, repeat=order):
            analytic_tensor[tuple(c.index for c in cartcoords)] = btens[(coord,) + cartcoords]

    @long_test
    @profile("bond_angle_6", "bond_angle", "b_tensor_all")
    def test_bond_angle_6_dry_run(self):
        self.setUp()
        order = 6
        coord = self.coord
        natoms = self.mol.natoms
        btens = coord.parent_representation.b_tensor()
        analytic_tensor = Tensor(shape=(3*natoms,)*order)
        for cartcoords in product(coord.variables, repeat=order):
            analytic_tensor[tuple(c.index for c in cartcoords)] = btens[(coord,) + cartcoords]


@allow_generators
class BTensorTorsionTest(unittest.TestCase):

    def setUpForValue(self, val, randomize_cartesian_space=True):
        self.mol = Molecule.from_z_matrix('''
            O
            H 1 0.93
            C 2 1.23 1 104.5
            B 3 1.35 2  89.5 1 ''' + str(val), create_representation=True)
        self.mol.recenter()
        if randomize_cartesian_space:
            # 'randomize' the position in xyz space to avoid miniscule displacement amounts
            self.mol.rotate(Vector(0,0,1), 32 * Degrees.to(Radians))
            self.mol.rotate(Vector(0,1,0), 32 * Degrees.to(Radians))
            self.mol.rotate(Vector(1,0,0), 32 * Degrees.to(Radians))
        self.coord = self.mol.internal_representation[5]


    #--------------------------------------------------------------------------------#
    # Higher order B tensor tests                                                    #
    #--------------------------------------------------------------------------------#

    def test_torsion_2(self):
        self.setUpForValue(84.3)
        assert_has_valid_b_tensor(self.coord, 2, print_precision=4, print_line_width=120)

    @profile("torsion_2", "torsion", "b_tensor_short", "b_tensor_all")
    def test_torsion_2_dry_run(self):
        self.setUpForValue(84.3)
        order = 2
        coord = self.coord
        natoms = self.mol.natoms
        btens = coord.parent_representation.b_tensor()
        analytic_tensor = Tensor(shape=(3*natoms,)*order)
        for cartcoords in product(coord.variables, repeat=order):
            analytic_tensor[tuple(c.index for c in cartcoords)] = btens[(coord,) + cartcoords]

    @long_test
    def test_torsion_3(self):
        self.setUpForValue(84.3)
        assert_has_valid_b_tensor(self.coord, 3,
            robustness=4,
            print_precision=4,
            places=5,
            print_line_width=120)

    @long_test
    @profile("torsion_3", "torsion", "b_tensor_short", "b_tensor_all")
    def test_torsion_3_dry_run(self):
        self.setUpForValue(84.3)
        order = 3
        coord = self.coord
        natoms = self.mol.natoms
        btens = coord.parent_representation.b_tensor()
        analytic_tensor = Tensor(shape=(3*natoms,)*order)
        for cartcoords in product(coord.variables, repeat=order):
            analytic_tensor[tuple(c.index for c in cartcoords)] = btens[(coord,) + cartcoords]

    @long_test
    def test_torsion_4(self):
        self.setUpForValue(84.3)
        assert_has_valid_b_tensor(self.coord, 4,
            robustness=4,
            print_precision=3,
            places=4,
            print_line_width=120,
            show_vals=True)

    @long_test
    @profile("torsion_4", "torsion", "b_tensor_all")
    def test_torsion_4_dry_run(self):
        self.setUpForValue(84.3)
        order = 4
        coord = self.coord
        natoms = self.mol.natoms
        btens = coord.parent_representation.b_tensor()
        analytic_tensor = Tensor(shape=(3*natoms,)*order)
        for cartcoords in product(coord.variables, repeat=order):
            analytic_tensor[tuple(c.index for c in cartcoords)] = btens[(coord,) + cartcoords]

    @long_test
    @profile("multi_torsion_4")
    def test_2_torsions_4_dry_run(self):
        mol = Molecule.from_z_matrix('''
            O
            H 1 0.93
            C 2 1.23 1 104.5
            B 3 1.35 2  89.5 1 83.2
            O 4 1.25 3  89.5 2 -23.2
        ''', create_representation=True)
        mol.recenter()
        order = 4
        coord = mol.internal_representation[5]
        natoms = mol.natoms
        btens = coord.parent_representation.b_tensor()
        analytic_tensor = Tensor(shape=(3*natoms,)*order)
        for cartcoords in product(coord.variables, repeat=order):
            analytic_tensor[tuple(c.index for c in cartcoords)] = btens[(coord,) + cartcoords]
        coord = mol.internal_representation[8]
        natoms = mol.natoms
        btens = coord.parent_representation.b_tensor()
        analytic_tensor = Tensor(shape=(3*natoms,)*order)
        for cartcoords in product(coord.variables, repeat=order):
            analytic_tensor[tuple(c.index for c in cartcoords)] = btens[(coord,) + cartcoords]

    @long_test
    def test_torsion_5(self):
        self.setUpForValue(84.3)
        assert_has_valid_b_tensor(self.coord, 5)

    @long_test
    @profile("torsion_5", "torsion", "b_tensor_all")
    def test_torsion_4_dry_run(self):
        self.setUpForValue(84.3)
        order = 5
        coord = self.coord
        natoms = self.mol.natoms
        btens = coord.parent_representation.b_tensor()
        analytic_tensor = Tensor(shape=(3*natoms,)*order)
        for cartcoords in product(coord.variables, repeat=order):
            analytic_tensor[tuple(c.index for c in cartcoords)] = btens[(coord,) + cartcoords]

    #--------------------------------------------------------------------------------#
    # B Vector tests                                                                 #
    #--------------------------------------------------------------------------------#

    def test_torsion_1(self):
        self.setUpForValue(84.3)
        assert_has_valid_b_tensor(self.coord, 1)

    def test_torsion_1_90(self):
        self.setUpForValue(90)
        assert_has_valid_b_tensor(self.coord, 1, places=6)

    def test_torsion_1_neg_89(self):
        self.setUpForValue(-89)
        assert_has_valid_b_tensor(self.coord, 1, robustness=15)

    def test_torsion_1_neg_90(self):
        self.setUpForValue(-90)
        assert_has_valid_b_tensor(self.coord, 1, robustness=15)

    def test_torsion_1_neg_45(self):
        self.setUpForValue(-45)
        assert_has_valid_b_tensor(self.coord, 1, robustness=15, use_coordinate_indices=True)

    @long_generator
    def torsion_1_test_generator(cls):
        npoints = 75
        minval = -185
        maxval = 185
        def _generated(self, tval):
            self.setUpForValue(tval)
            assert_has_valid_b_tensor(self.coord, 1, places=8, robustness=4)
        diff = maxval - minval
        step = diff / (npoints-1)
        for i in range(npoints):
            currval = minval + (i * step)
            _generated.__name__ = 'test_gen_bvect_tors_' + re.sub('-', 'neg_', re.sub('[\. ]', '_', str(currval)))
            yield _generated, currval

class FunctionCoordinateBTest(unittest.TestCase):

    def setUp(self):
        init_sample_molecules()
        self.mol = SampleMolecules['water']
        self.mol.recenter()
        rep = InternalRepresentation(self.mol,
            """
                STRE 1 2
                STRE 1 3
                BEND 2 1 3
            """
        )
        # 'randomize' the position in xyz space to avoid miniscule displacement amounts
        self.mol.rotate(Vector(1,2,1), 32.*Degrees)

    #--------------------------------------------------------------------------------#

    def test_value(self):
        coord = SumOfInternalCoordinates(subcoordinates=[
            self.mol.internal_representation[0],
            self.mol.internal_representation[1]
        ], coefficients = [1.0, -1.0])
        self.assertAlmostEqual(coord.value, 0)


    def test_sum_1(self):
        coord = SumOfInternalCoordinates(subcoordinates=[
            self.mol.internal_representation[0],
            self.mol.internal_representation[1]
        ])
        assert_has_valid_b_tensor(coord, 1, use_coordinate_indices=True)

    def test_sum_1_sub(self):
        coord = SumOfInternalCoordinates(subcoordinates=[
            self.mol.internal_representation[0],
            self.mol.internal_representation[1]
        ], coefficients = [1.0, -0.5])
        assert_has_valid_b_tensor(coord, 1, use_coordinate_indices=True)

    def test_sum_2(self):
        coord = SumOfInternalCoordinates(subcoordinates=[
            self.mol.internal_representation[0],
            self.mol.internal_representation[1]
        ])
        assert_has_valid_b_tensor(coord, 2, use_coordinate_indices=True)

    def test_sum_2_sub(self):
        coord = SumOfInternalCoordinates(subcoordinates=[
            self.mol.internal_representation[0],
            self.mol.internal_representation[1]
        ], coefficients = [1/math.sqrt(2), -1/math.sqrt(2)])
        assert_has_valid_b_tensor(coord, 2, use_coordinate_indices=True)

    def test_sum_3(self):
        coord = SumOfInternalCoordinates(subcoordinates=[
            self.mol.internal_representation[0],
            self.mol.internal_representation[1]
        ])
        assert_has_valid_b_tensor(coord, 3, use_coordinate_indices=True)

    def test_sum_3_sub(self):
        coord = SumOfInternalCoordinates(subcoordinates=[
            self.mol.internal_representation[0],
            self.mol.internal_representation[1]
        ], coefficients = [37, -4.37])
        assert_has_valid_b_tensor(coord, 2, use_coordinate_indices=True)

    def test_sum_4(self):
        coord = SumOfInternalCoordinates(subcoordinates=[
            self.mol.internal_representation[0],
            self.mol.internal_representation[1]
        ])
        assert_has_valid_b_tensor(coord, 4, use_coordinate_indices=True)

class InternalCartesianBTest(unittest.TestCase):

    def setUp(self):
        self.mol = Molecule("""
        4
        H2O2
        H    -1.21670000  -0.75630000   0.00000000
        O    -0.73020000   0.07940000  -0.00000000
        O     0.73020000  -0.07940000  -0.00000000
        H     1.21670000   0.75630000   0.00000000
        """)
        self.rep = InternalRepresentation(self.mol, [
                BondLength(1, 2, self.mol),
                BondLength(1, 3, self.mol),
                BondAngle(2, 1, 3, self.mol),
                InternalCartesianX(1, 2, 3, 4, self.mol),
                InternalCartesianY(1, 2, 3, 4, self.mol),
                InternalCartesianZ(1, 2, 3, 4, self.mol)
            ]
        )

    def test_intcart_1(self):
        assert_has_valid_b_tensor(self.rep[3], 1)
        assert_has_valid_b_tensor(self.rep[4], 1)
        assert_has_valid_b_tensor(self.rep[5], 1)

#class SinPhiBTest(unittest.TestCase):
#
#    def setUp(self):
#        self.mol = Molecule.from_z_matrix('''
#            O
#            H 1 0.93
#            C 2 1.23 1 104.5
#        ''', create_representation=True)
#        coords = [
#            BondLength(self.mol[0], self.mol[1]),
#            BondLength(self.mol[1], self.mol[2]),
#            SinBondAngleCoordinate(self.mol.internal_representation[2], delta=0.001)
#        ]
#        rep = InternalRepresentation(self.mol, coords)
#        self.mol.recenter()
#        self.coord = self.mol.internal_representation[2]
#
#    def test_sinphi_2(self):
#        assert_has_valid_b_tensor(self.coord, 2, print_precision=4, print_line_width=120)
#
#    @expected_failure
#    def test_sinphi_3(self):
#        assert_has_valid_b_tensor(self.coord, 3, print_precision=4, print_line_width=120)
#
#    def test_sinphi_4(self):
#        assert_has_valid_b_tensor(self.coord, 4, print_precision=4, print_line_width=120)

