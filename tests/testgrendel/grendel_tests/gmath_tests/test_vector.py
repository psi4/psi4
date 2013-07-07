
import unittest
import sys
import os

# Add the directory containing 'grendel' to the sys.path
sys.path.append(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir, os.pardir))

# Add the directory containing the 'grendel_tests' package to sys.path
sys.path.append(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir))

from grendel import *
import numpy as np

class VectorTest(unittest.TestCase):

    def test_cartesian(self):
        v2 = Vector(1,2)
        v3 = Vector(1,2,3)
        v4 = Vector(1,2,3,4)
        with self.assertRaises(IndexError):
            v2.z
        with self.assertRaises(IndexError):
            v2.z = 1.5
        with self.assertRaises(IndexError):
            v4.x
        with self.assertRaises(IndexError):
            v4.x = 1.5
        with self.assertRaises(IndexError):
            v4.y
        with self.assertRaises(IndexError):
            v4.y = 1.5
        with self.assertRaises(IndexError):
            v4.z
        with self.assertRaises(IndexError):
            v4.z = 1.5
        v3.x, v3.y, v3.z = 3, 2, 1
        self.assertEqual(v3, Vector(3,2,1))
        self.assertEqual(v3.x, 3)
        self.assertEqual(v3.y, 2)
        self.assertEqual(v3.z, 1)

    def test_misc(self):
        # Reshape
        self.assertEqual(Vector(range(1,9)).reshape((2,2,2)), Tensor([[1,2],[3,4]],[[5,6],[7,8]]))

        # normalized, test that a copy is made
        v = Vector(1,2,3)
        v2 = v.normalized()
        v2[2] = 5
        self.assertNotEqual(v2[2], v[2])
        v.normalize()
        self.assertNotEqual(v2, v)
        v2 = Vector(1,2,3).normalized()
        self.assertEqual(v2, v)

        # cross
        v = Vector(range(3))
        self.assertSequenceEqual(v.cross(v+1), Vector([-1, 2, -1]))
        self.assertIs(type(cross(Tensor([1,2,3]), Tensor([4,5,6]))), Tensor)
        self.assertIs(type(cross(np.array([1,2,3]), np.array([4,5,6]))), np.ndarray)

        # magnitude
        self.assertEqual(magnitude(Vector(3,4)), 5.0)


    def test_errors(self):
        with self.assertRaises(ZeroDivisionError):
            Vector(0,0,0).normalize()

    def test_multiply(self):
        m = Matrix(range(9), shape=(3,3))
        v = Vector(range(3))
        self.assertSequenceEqual(m*v, Vector(5, 14, 23))

