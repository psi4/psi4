from copy import copy, deepcopy
import unittest
import sys
import os

# Add the directory containing 'grendel' to the sys.path
sys.path.append(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir, os.pardir))

# Add the directory containing the 'grendel_tests' package to sys.path
sys.path.append(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir))

from grendel import *
from grendel_tests import allow_generators
from grendel_tests.gmath_tests.test_einsum import range_tensor, index_range_tensor

@allow_generators
class TensorTest(unittest.TestCase):

    sample_tensors = []

    initial_sample_tensors = [
        Tensor([0]),

        # Should initialize to a 1x1x1 Tensor with the element 1.0
        Tensor(
            [ [ [ 1 ] ] ]
        ),

        # Should initialize to a 1x1x1x1 Tensor with the element 1.0
        Tensor(
            [ [ [ [ 1 ] ] ] ]
        ),

        # Should be a 3x4x5 Tensor with all elements 3.14
        Tensor(shape=(3,4,5), default_val=3.14),

        # Should be a 5x4x3 Tensor with all elements 3.14
        Tensor(shape=(5,4,3), default_val=3.14),

        # Should be a 4x4 Tensor with all elements 1.0-16.0
        Tensor(
            [ 1,  2,  3,  4],
            [ 5,  6,  7,  8],
            [ 9, 10, 11, 12],
            [13, 14, 15, 16]
        ),

        # Should be a 4x2x2 Tensor
        Tensor(
            [[ 1,  2], [ 3,  4]],
            [[ 5,  6], [ 7,  8]],
            [[ 9, 10], [11, 12]],
            [[13, 14], [15, 16]]
        ),

        # Should be a 2x2x2x2 Tensor
        Tensor(
            [[[ 1,  2], [ 3,  4]],
                [[ 5,  6], [ 7,  8]]],
            [[[ 9, 10], [11, 12]],
                [[13, 14], [15, 16]]]
        ),
    ]

    def setUp(self):
        self.sample_tensors = deepcopy(self.initial_sample_tensors)

    def tearDown(self):
        self.sample_tensors = []

    @classmethod
    def eq_test_generator(cls):
        def _generated(self, tens1, tens2, expectedEq):
            if expectedEq:
                self.assertEqual(tens1, tens2)
            else:
                self.assertNotEqual(tens1, tens2)

        for i, itens in enumerate(cls.initial_sample_tensors):
            for j, jtens in enumerate(cls.initial_sample_tensors):
                _generated.__name__ = "test_eq_" + str(i) + "_" + str(j)
                yield _generated, itens, jtens, i==j

    def test_zero_1(self):
        self.assertTrue(self.sample_tensors[0].is_zero())

    def test_zero_2(self):
        self.assertFalse(self.sample_tensors[4].is_zero())

    def test_chop(self):
        t = range_tensor((3,3,3)) + 1e-15
        self.assertSequenceEqual(chopped(t), range_tensor((3,3,3)))
        chop(t)
        self.assertIsNot(t, chopped(t))
        self.assertSequenceEqual(t, range_tensor((3,3,3)))


class ComputedTensorTest(unittest.TestCase):

    def test_get_item_1(self):
        def idx_sum(tens, indices):
            return sum(indices)
        t = ComputableTensor(
            shape=(3,4,5),
            uncomputed=True,
            compute_function=idx_sum)
        for i in xrange(3):
            for j in xrange(4):
                for k in xrange(5):
                    self.assertEqual(t[i,j,k], i + j + k)

    def test_get_item_2(self):
        def idx_sum(tens, indices):
            return sum(indices)
        t = ComputableTensor(
            shape=(3,4,5),
            uncomputed=True,
            compute_function=idx_sum)
        for i in xrange(3):
            for j in xrange(4):
                for k in xrange(5):
                    self.assertEqual(t[i][j][k], i + j + k)

    def test_get_item_3(self):
        global callcount
        callcount = 0
        def idx_sum(tens, indices):
            tens[indices[0]] = indices[0]
            global callcount
            callcount += 1
            return indices[0]
        t = ComputableTensor(
            shape=(3,4,5),
            uncomputed=True,
            compute_function=idx_sum)
        for i in xrange(3):
            for j in xrange(4):
                for k in xrange(5):
                    self.assertEqual(t[i][j][k], i)
        self.assertEqual(callcount, 3)


