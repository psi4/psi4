""" A module containing the abstract superclass for elements of point groups.
"""
from abc import ABCMeta, abstractproperty
from functools import total_ordering
import math
import numpy
from grendel.gmath.geometry import angle_between_vectors
from grendel.gmath.vector import magnitude, Vector
import grendel as g
from grendel.util.strings import classname
from grendel.util.units import *

__all__ = [
    'SymmetryOperation'
]

@total_ordering
class SymmetryOperation(object):
    """ The abstract superclass for elements of point groups.


    Attributes
    ----------
    name : str
        An attempt at creating a human-readable name for the operation.


    """
    __metaclass__ = ABCMeta

    ####################
    # Class Attributes #
    ####################

    same_axis_tolerance = 1e-4
    same_operation_tolerance = 1e-4
    zero_vector_cutoff = 1e-5
    perpendicular_tolerance = 1.0 * Degrees.to(Radians)  #TODO compute this from symmetry tolerance


    ##############
    # Attributes #
    ##############

    point_group = None
    conjugacy_class = None
    inverse = None
    name = None


    ######################
    # Private Attributes #
    ######################

    _matrix = None


    #######################
    # Abstract Properties #
    #######################

    @abstractproperty
    def matrix(self):
        """ The Matrix that transforms the molecule according to the symmetry operation `self`
        """
        return NotImplemented


    ##############
    # Properties #
    ##############

    @property
    def molecule(self):
        """ The Molecule that the operation acts on
        """
        return self.point_group.molecule

    ###################
    # Special Methods #
    ###################

    def __eq__(self, other):
        if not isinstance(other, SymmetryOperation):
            raise TypeError("Equality comparison of SymmetryOperation to " + classname(other) + " is not allowed.")
        diff = self.matrix - other.matrix
        return diff.norm() < self.same_operation_tolerance

    # TODO prefer ordering by class containment if self.point_group.classes is populated
    def __lt__(self, other):
        """ Determines a sort order based on Cotton ordering (or as close as possible to it...)
        (i.e. the ordering used by Cotton in his character tables)
        Note:  Cotton isn't really that consistent with his ordering, particularly for larger point groups.  This
        represents my best guess, basically just to get something consistant out on the table.
        """
        if self is other:
            return False
        if not isinstance(other, SymmetryOperation):
            raise TypeError("Comparison of SymmetryOperation with " + classname(other) + " is not allowed.")
        elif not self.point_group is other.point_group:
            raise ValueError("Cannot compare operations from different point groups.")
        if isinstance(self, g.IdentityOperation):
            return True
        if isinstance(other, g.IdentityOperation):
            return False
        elif isinstance(self, g.Rotation):
            if not isinstance(other, g.Rotation):
                return True
            elif (self.n, self.n - self.exponent) != (other.n, other.n - other.exponent):
                # Reverse ordering by order, forward ordering by exponent
                return (self.n, self.n - self.exponent) > (other.n, other.n - other.exponent)
            # Then put in the order z, y, x
            elif SymmetryOperation.is_same_axis(self.axis, Vector(0,0,1)):
                return True
            elif SymmetryOperation.is_same_axis(other.axis, Vector(0,0,1)):
                return False
            elif SymmetryOperation.is_same_axis(self.axis, Vector(0,1,0)):
                return True
            elif SymmetryOperation.is_same_axis(other.axis, Vector(0,1,0)):
                return False
            elif SymmetryOperation.is_same_axis(self.axis, Vector(1,0,0)):
                return True
            elif SymmetryOperation.is_same_axis(other.axis, Vector(1,0,0)):
                return False
            # Finally, order by proximity of axis to the z-axis (not an official part of Cotton ordering, but doesn't come up too often
            else:
                return (Vector(0,0,1) - self.axis).magnitude() < (Vector(0,0,1) - other.axis).magnitude()
        elif isinstance(other, g.Rotation):
            return False
        elif isinstance(self, g.Inversion):
            return True
        elif isinstance(other, g.Inversion):
            return False
        elif isinstance(self, g.ImproperRotation):
            if not isinstance(other, g.ImproperRotation):
                return True
            elif (self.n, self.n - self.exponent) != (other.n, other.n - other.exponent):
                # Reverse ordering by order, forward ordering by exponent
                return (self.n, self.n - self.exponent) > (other.n, other.n - other.exponent)
            # Then put in the order z, y, x
            elif SymmetryOperation.is_same_axis(self.axis, Vector(0,0,1)):
                return True
            elif SymmetryOperation.is_same_axis(other.axis, Vector(0,0,1)):
                return False
            elif SymmetryOperation.is_same_axis(self.axis, Vector(0,1,0)):
                return True
            elif SymmetryOperation.is_same_axis(other.axis, Vector(0,1,0)):
                return False
            elif SymmetryOperation.is_same_axis(self.axis, Vector(1,0,0)):
                return True
            elif SymmetryOperation.is_same_axis(other.axis, Vector(1,0,0)):
                return False
            # Finally, order by proximity of axis to the z-axis (not an official part of Cotton ordering, but doesn't come up too often
            else:
                return (Vector(0,0,1) - self.axis).magnitude() < (Vector(0,0,1) - other.axis).magnitude()
        elif isinstance(other, g.ImproperRotation):
            return False
        else:  # Both are Reflections
            if self.is_principal_reflection():
                return True
            elif other.is_principal_reflection():
                return False
            elif self.is_dihedral():
                if not other.is_dihedral():
                    return True
                # Otherwise same axis ordering as rotations
                elif SymmetryOperation.is_same_axis(self.axis, Vector(0,0,1)):
                    return True
                elif SymmetryOperation.is_same_axis(other.axis, Vector(0,0,1)):
                    return False
                elif SymmetryOperation.is_same_axis(self.axis, Vector(0,1,0)):
                    return True
                elif SymmetryOperation.is_same_axis(other.axis, Vector(0,1,0)):
                    return False
                elif SymmetryOperation.is_same_axis(self.axis, Vector(1,0,0)):
                    return True
                elif SymmetryOperation.is_same_axis(other.axis, Vector(1,0,0)):
                    return False
                else:
                    return (Vector(0,0,1) - self.axis).magnitude() < (Vector(0,0,1) - other.axis).magnitude()
            elif other.is_dihedral():
                return False
            else: # Both are sigma_v's
                # Same axis ordering as rotations
                if SymmetryOperation.is_same_axis(self.axis, Vector(0,0,1)):
                    return True
                elif SymmetryOperation.is_same_axis(other.axis, Vector(0,0,1)):
                    return False
                elif SymmetryOperation.is_same_axis(self.axis, Vector(0,1,0)):
                    return True
                elif SymmetryOperation.is_same_axis(other.axis, Vector(0,1,0)):
                    return False
                elif SymmetryOperation.is_same_axis(self.axis, Vector(1,0,0)):
                    return True
                elif SymmetryOperation.is_same_axis(other.axis, Vector(1,0,0)):
                    return False
                else:
                    return (Vector(0,0,1) - self.axis).magnitude() < (Vector(0,0,1) - other.axis).magnitude()

    # TODO Cache the multiplication table in the parent PointGroup object
    def __mul__(self, other):
        return self.point_group.element_with_matrix(self.matrix * other.matrix)

    #################
    # Class Methods #
    #################

    @classmethod
    def is_same_axis(cls, a1, a2, parallel_only=False):
        """ True if the axes `a1` and `a2` are parallel or antiparallel.
        Note that normalized versions of `a1`  and `a2` are used, so normalizing before
        passing in will just slow things down.  If `parallel_only` is `True` (it is `False` by default), this method
        only returns `True` if the two axes are parallel, not anti_parallel.
        """
        n1 = a1.normalized()
        n2 = a2.normalized()
        if parallel_only:
            return magnitude(n1-n2) < cls.same_axis_tolerance
        else:
            return magnitude(n1-n2) < cls.same_axis_tolerance or magnitude(n1+n2) < cls.same_axis_tolerance

    @classmethod
    def is_same_matrix(cls, m1, m2):
        """ True if `(a1-a2).norm() < SymmetryOperation.same_operation_tolerance`
        (See `~pyobj:grendel.gmath.tensor.Tensor.norm` for description of what is meant by norm here)
        """
        return (m1-m2).norm() < cls.same_operation_tolerance

    @classmethod
    def is_perpendicular(cls, v1, v2):
        return abs(math.pi / 2.0 - angle_between_vectors(v1, v2)) < cls.perpendicular_tolerance

    ###########
    # Methods #
    ###########



