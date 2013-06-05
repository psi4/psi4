"""
"""
from grendel.gmath import divisors
from grendel.gmath.matrix import Matrix
from grendel.gmath.vector import Vector, cross
from grendel.util.decorators import CachedMethod
from grendel.util.strings import classname
from symmetry_operation import SymmetryOperation

__all__ = [
    'Reflection'
]

class Reflection(SymmetryOperation):
    """ A reflection (sigma) element of a point group.
    """

    ##############
    # Attributes #
    ##############

    axis = None
    principal_reflection = False # True if the reflection is a sigma_h

    ##################
    # Initialization #
    ##################

    def __init__(self, axis, point_group = None):
        """
        """
        self.axis = axis
        self.axis.normalize()
        self.point_group = point_group


    ##############
    # Properties #
    ##############

    @property
    def matrix(self):
        if not self._matrix is None:
            return self._matrix
        x = self.axis[0]
        y = self.axis[1]
        z = self.axis[2]
        self._matrix = Matrix([
            [1 - 2*x*x, -2*x*y, -2*x*z],
            [-2*x*y, 1 - 2*y*y, -2*y*z],
            [-2*x*z, -2*y*z, 1 - 2*z*z]
        ])
        return self._matrix


    ###################
    # Special Methods #
    ###################

    def __str__(self):
        if not self.name is None:
            return self.name
        else:
            return "sigma(axis=" + str(self.axis) + ")"

    def __repr__(self):
        return "Reflection(" + repr(self.axis) + ")"

    def __mul__(self, other):
        if not isinstance(other, SymmetryOperation):
            raise TypeError("Invalid multiplication of " + classname(self) + " by " + classname(other))
        return super(Reflection, self).__mul__(other)


    #################
    # Class Methods #
    #################

    @classmethod
    def with_normal(cls, molecule, axis):
        op = Reflection(axis)
        if molecule.has_symmetry(op):
            return op
        else:
            return None


    ###########
    # Methods #
    ###########

    def is_principal_reflection(self):
        """
        """
        return self.principal_reflection

    def is_dihedral(self):
        """ Whether or not the reflection is a sigma_d
        """
        if self.principal_reflection:
            return False
        for atom in self.molecule:
            # If there's an atom in the plane, it's a sigma_v.  If there's an atom in the direction of the normal (or it's negative),
            # it's sigma_v.  Otherwise, it's a sigma_d.  If the atom is at the origin, ignore it.
            if atom.pos.is_zero():
                continue
            elif SymmetryOperation.is_same_axis(self.axis, atom.pos):
                return False
            elif SymmetryOperation.is_perpendicular(self.axis, atom.pos):
                return False
        return True

