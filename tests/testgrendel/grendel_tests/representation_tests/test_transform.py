import sys
import os
import unittest
import numpy as np
from copy import copy, deepcopy
from itertools import permutations

from numpy.testing import assert_array_almost_equal

# Add the directory containing 'grendel' to the sys.path
sys.path.append(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir, os.pardir))

# Add the directory containing the 'grendel_tests' package to sys.path
sys.path.append(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir))

from grendel import *
from grendel.chemistry import SampleMolecules, init_sample_molecules
from grendel_tests import skip, profile, long_test

#noinspection PyUnresolvedReferences
class TransformationTest(unittest.TestCase):

    def setUp(self):
        init_sample_molecules()

    @skip
    def test_cartesian_180(self):
        tmp = Vector([0,0,1,0,0,-1])
        axis = Vector([1,1,0]).normalized()
        mol1 = Molecule("""
            O 0 1 0
            H 0 0 0
        """)
        mol1.recenter()
        rep1 = mol1.cartesian_representation
        rtens = RepresentationDependentTensor(deepcopy(tmp), representation=rep1)
        mol2 = mol1.rotated(angle=180.*Degrees, axis=axis)
        rep2 = mol2.cartesian_representation
        rtens = rtens.in_representation(rep2)
        assert_array_almost_equal(chopped(rtens.view(np.ndarray)), chopped(-tmp))

    def test_simplest_cartesian(self):
        axis = Vector([1,2,3]).normalized()
        mol1 = SampleMolecules['benzene']
        tmp = np.arange(1,3*mol1.natoms + 1)
        mol1.recenter()
        rep1 = mol1.cartesian_representation
        rtens = RepresentationDependentTensor(deepcopy(tmp), representation=rep1)
        mol2 = mol1.rotated(angle=360. * Degrees, axis=axis)
        rep2 = mol2.cartesian_representation
        rtens = rtens.in_representation(rep2)
        assert_array_almost_equal(rtens.view(np.ndarray), tmp)

    def test_simpler_cartesian(self):
        axis = Vector([1,2,3]).normalized()
        #----------------------------------------#
        mol1 = SampleMolecules['water']
        tmp = np.arange(1,3*mol1.natoms + 1)
        mol1.recenter()
        rep1 = mol1.cartesian_representation
        rtens = RepresentationDependentTensor(deepcopy(tmp), representation=rep1)
        #----------------------------------------#
        mol2 = mol1.rotated(angle=180. * Degrees, axis=axis)
        rep2 = mol2.cartesian_representation
        rtens = rtens.in_representation(rep2)
        #----------------------------------------#
        mol3 = mol2.rotated(angle=180. * Degrees, axis=axis)
        rep3 = mol3.cartesian_representation
        rtens = rtens.in_representation(rep3)
        #----------------------------------------#
        assert_array_almost_equal(rtens.view(np.ndarray), tmp)

    def test_cartesian_norms(self):
        axis = Vector([1,2,3]).normalized()
        mol1 = Molecule.get("glycine") #SampleMolecules['benzene']
        tmp = np.arange(1,3*mol1.natoms + 1)
        mol1.recenter()
        #----------------------------------------#
        rep1 = mol1.cartesian_representation
        rtens = RepresentationDependentTensor(deepcopy(tmp), representation=rep1)
        #----------------------------------------#
        mol2 = mol1.rotated(angle=120.*Degrees, axis=axis)
        rep2 = mol2.cartesian_representation
        rtens = rtens.in_representation(rep2)
        #----------------------------------------#
        for orig, trans in zip(grouper(3, tmp), grouper(3, rtens.view(np.ndarray))):
            self.assertAlmostEqual(Vector(orig).magnitude(), Vector(trans).magnitude())


    def test_simple_cartesian(self):
        axis = Vector([-1,-2,-3]).normalized()
        mol1 = SampleMolecules['water']
        tmp = np.arange(1,3*len(mol1.atoms) + 1)
        mol1.recenter()
        #----------------------------------------#
        rep1 = mol1.cartesian_representation
        rtens = RepresentationDependentTensor(deepcopy(tmp), representation=rep1)
        #----------------------------------------#
        mol2 = mol1.rotated(angle=120.*Degrees, axis=axis)
        rep2 = mol2.cartesian_representation
        rtens = rtens.in_representation(rep2)
        #----------------------------------------#
        mol3 = mol2.rotated(angle=120.*Degrees, axis=axis)
        rep3 = mol3.cartesian_representation
        rtens = rtens.in_representation(rep3)
        #----------------------------------------#
        mol4 = mol3.rotated(angle=120.*Degrees, axis=axis)
        rep4 = mol4.cartesian_representation
        rtens = rtens.in_representation(rep4)
        #----------------------------------------#
        assert_array_almost_equal(rtens.view(np.ndarray), tmp)

    @skip
    def test_diatomic_cartesian_to_internal(self):
        tmp = Vector([0,0,-1,0,0,1])
        axis = Vector([1,1,0]).normalized()
        #axis = Vector([-1,-2,-3]).normalized()
        #----------------------------------------#
        mol1 = Molecule("""
            H 0 0 1
            H 0 0 0
        """)
        mol1.recenter()
        rep1 = InternalRepresentation(mol1, """
            STRE 1 2
        """)
        rtens1 = RepresentationDependentTensor(
            deepcopy(tmp),
            representation=mol1.cartesian_representation
        )
        print(rtens1)
        print(rep1.b_matrix)
        print(chopped(np.linalg.pinv(rep1.b_matrix)))
        #----------------------------------------#
        mol2 = mol1.rotated(axis=axis, angle=180.*Degrees)
        rep2 = InternalRepresentation(mol2, """
            STRE 1 2
        """)
        rtens2 = rtens1.in_representation(mol2.cartesian_representation)
        print(rtens2)
        print(chopped(rep2.b_matrix))
        print(chopped(np.linalg.pinv(rep2.b_matrix)))
        #----------------------------------------#
        assert_array_almost_equal(
            rtens1.in_representation(rep1).view(np.ndarray),
            rtens2.in_representation(rep2).view(np.ndarray)
        )

    def test_cartesian_to_internal(self):
        tmp = Vector([0,2,3,0,-2,3,0,0,1])
        axis = Vector([1,0,0]).normalized()
        #axis = Vector([-1,-2,-3]).normalized()
        #----------------------------------------#
        mol1 = SampleMolecules['water']
        mol1.recenter()
        rep1 = InternalRepresentation(mol1, """
            STRE 1 2
            STRE 1 3
            BEND 2 1 3
        """)
        rtens1 = RepresentationDependentTensor(
            deepcopy(tmp),
            representation=mol1.cartesian_representation
        )
        print(rtens1)
        #----------------------------------------#
        mol2 = mol1.rotated(axis=axis, angle=45.*Degrees)
        rep2 = InternalRepresentation(mol2, """
            STRE 1 2
            STRE 1 3
            BEND 2 1 3
        """)
        rtens2 = rtens1.in_representation(mol2.cartesian_representation)
        print(rtens2)
        #----------------------------------------#
        assert_array_almost_equal(
            rtens1.in_representation(rep1).view(np.ndarray),
            rtens2.in_representation(rep2).view(np.ndarray)
        )

    def test_i2c2i_2(self):
        data = open(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir, "test_files", "hooh_ff.dat")).read()
        sections = re.split(r'\s+0\s+', data, flags=re.M)
        molstr = ''.join(sum(map(list, zip(['H ', 'O', 'O', 'H'], list(sections[1].splitlines(True)))), []))
        mol = Molecule(molstr)
        rep = InternalRepresentation(mol, sections[0])
        ff = ForceField(max_order=2, representation=rep, property=Energy)
        ff.for_order(1)[...] = 0.0
        #for i, section in enumerate(sections[2]):
        for line in sections[2].splitlines():
            parts = line.split()
            parts = [p for p in parts if p != '']
            if len(parts) == 3:
                idxs = map(int, parts[:-1])
                idxs = tuple(i-1 for i in idxs)
                val = float(parts[-1])
                for perm in permutations(idxs):
                    ff[perm] = val
        assert_array_almost_equal(
            ff.for_order(2).view(Tensor),
            ff.in_representation(mol.cartesian_representation).in_representation(rep).for_order(2).view(Tensor)
        )

    @profile('transform')
    def test_i2c_ff_3_dry_run(self):
        data = open(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir, "test_files", "hooh_ff.dat")).read()
        sections = re.split(r'\s+0\s+', data, flags=re.M)
        molstr = ''.join(sum(map(list, zip(['H ', 'O', 'O', 'H'], list(sections[1].splitlines(True)))), []))
        mol = Molecule(molstr)
        rep = InternalRepresentation(mol, sections[0])
        ff = ForceField(max_order=3, representation=rep, property=Energy)
        ff.for_order(1)[...] = 0.0
        for i, section in enumerate(sections[2:4]):
            for line in section.splitlines():
                parts = line.split()
                parts = [p for p in parts if p != '']
                if len(parts) == 3 + i:
                    idxs = map(int, parts[:-1])
                    idxs = tuple(i-1 for i in idxs)
                    val = float(parts[-1])
                    for perm in permutations(idxs):
                        ff[perm] = val
        ff.in_representation(mol.cartesian_representation)

    def test_i2c2i_ff_3(self):
        data = open(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir, "test_files", "hooh_ff.dat")).read()
        sections = re.split(r'\s+0\s+', data, flags=re.M)
        molstr = ''.join(sum(map(list, zip(['H ', 'O', 'O', 'H'], list(sections[1].splitlines(True)))), []))
        mol = Molecule(molstr)
        rep = InternalRepresentation(mol, sections[0])
        ff = ForceField(max_order=3, representation=rep, property=Energy)
        ff.for_order(1)[...] = 0.0
        for i, section in enumerate(sections[2:4]):
            for line in section.splitlines():
                parts = line.split()
                parts = [p for p in parts if p != '']
                if len(parts) == 3 + i:
                    idxs = map(int, parts[:-1])
                    idxs = tuple(i-1 for i in idxs)
                    val = float(parts[-1])
                    for perm in permutations(idxs):
                        ff[perm] = val
        assert_array_almost_equal(
            ff.for_order(3).view(Tensor),
            ff.in_representation(mol.cartesian_representation).in_representation(rep).for_order(3).view(Tensor)
        )

    def test_i2c2i_ff_3_random(self):
        mol = Molecule(
            """
                O  0.000000  0.000000  0.118154
                H  0.000000  0.758734 -0.472614
                H  0.000000 -0.758734 -0.472614
            """,
            units=Angstroms
        )
        rep = InternalRepresentation(mol, """
            STRE 1 2
            STRE 1 3
            BEND 2 1 3
        """)
        ff = ForceField(max_order=3, representation=rep, property=Energy)
        for o in xrange(1, ff.max_order+1):
            ff.for_order(o)[...] = np.random.random((len(rep),)*o)
        assert_array_almost_equal(
            ff.for_order(3).view(Tensor),
            ff.in_representation(mol.cartesian_representation).in_representation(rep).for_order(3).view(Tensor)
        )

    def test_i2c2i_ff_4_diatomic(self):
        mol = Molecule(
            """
                O  0.000000  0.000000  0.118154
                H  0.000000  0.758734 -0.472614
            """,
            units=Angstroms
        )
        rep = InternalRepresentation(mol, """
            STRE 1 2
        """)
        ff = ForceField(max_order=4, representation=rep, property=Energy)
        for o in xrange(1, ff.max_order+1):
            ff.for_order(o)[...] = np.random.random((len(rep),)*o)
        assert_array_almost_equal(
            ff.for_order(4).view(Tensor),
            ff.in_representation(mol.cartesian_representation).in_representation(rep).for_order(4).view(Tensor)
        )

    @long_test
    def test_i2c2i_ff_4_random(self):
        mol = Molecule(
            """
                O  0.000000  0.000000  0.118154
                H  0.000000  0.758734 -0.472614
                H  0.000000 -0.758734 -0.472614
            """,
            units=Angstroms
        )
        rep = InternalRepresentation(mol, """
            STRE 1 2
            STRE 1 3
            BEND 2 1 3
        """)
        ff = ForceField(max_order=4, representation=rep, property=Energy)
        for o in xrange(1, ff.max_order+1):
            ff.for_order(o)[...] = np.random.random((len(rep),)*o)
        assert_array_almost_equal(
            ff.for_order(4).view(Tensor),
            ff.in_representation(mol.cartesian_representation).in_representation(rep).for_order(4).view(Tensor)
        )

    @long_test
    def test_i2c2i_ff_5_random(self):
        mol = Molecule(
            """
                O  0.000000  0.000000  0.118154
                H  0.000000  0.758734 -0.472614
                H  0.000000 -0.758734 -0.472614
            """,
            units=Angstroms
        )
        rep = InternalRepresentation(mol, """
            STRE 1 2
            STRE 1 3
            BEND 2 1 3
        """)
        ff = ForceField(max_order=5, representation=rep, property=Energy)
        for o in xrange(1, ff.max_order+1):
            ff.for_order(o)[...] = np.random.random((len(rep),)*o)
        assert_array_almost_equal(
            ff.for_order(5).view(Tensor),
            ff.in_representation(mol.cartesian_representation).in_representation(rep).for_order(5).view(Tensor)
        )

    @long_test
    def test_i2c2i_ff_6_random_diatomic(self):
        mol = Molecule(
            """
                O  0.000000  0.000000  0.118154
                H  0.000000  0.758734 -0.472614
            """,
            units=Angstroms
        )
        rep = InternalRepresentation(mol, """
            STRE 1 2
        """)
        ff = ForceField(max_order=6, representation=rep, property=Energy)
        for o in xrange(1, ff.max_order+1):
            ff.for_order(o)[...] = np.random.random((len(rep),)*o)
        assert_array_almost_equal(
            ff.for_order(6).view(Tensor),
            ff.in_representation(mol.cartesian_representation).in_representation(rep).for_order(6).view(Tensor)
        )

    @long_test
    def test_i2c2i_ff_6_random(self):
        mol = Molecule(
            """
                O  0.000000  0.000000  0.118154
                H  0.000000  0.758734 -0.472614
                H  0.000000 -0.758734 -0.472614
            """,
            units=Angstroms
        )
        rep = InternalRepresentation(mol, """
            STRE 1 2
            STRE 1 3
            BEND 2 1 3
        """)
        ff = ForceField(max_order=6, representation=rep, property=Energy)
        for o in xrange(1, ff.max_order+1):
            ff.for_order(o)[...] = np.random.random((len(rep),)*o)
        assert_array_almost_equal(
            ff.for_order(6).view(Tensor),
            ff.in_representation(mol.cartesian_representation).in_representation(rep).for_order(6).view(Tensor)
        )

