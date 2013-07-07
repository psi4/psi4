"""
"""
from fractions import gcd
import math
from grendel.gmath import divisors
from grendel.gmath.matrix import Matrix
from grendel.util.strings import classname, superscript, subscript
from symmetry_operation import SymmetryOperation

__all__ = [
    'Rotation'
]


class Rotation(SymmetryOperation):
    """ A rotation (C_n^m) element of a point group.

    :Attributes:

    n : int
        The order of the rotation (e.g. for a C_3 operation, `n`=3)
    axis : Vector
        The axis about which the operation occurs
    exponent : int
        The number of times the operation is applied.  (e.g. for C_3, `exponent` = 1, for C_3^2, `exponent` = 2, etc.)

    """

    ####################
    # Class Attributes #
    ####################


    ##############
    # Attributes #
    ##############

    n = None
    axis = None
    exponent = None


    ##################
    # Initialization #
    ##################

    def __init__(self, n, axis, exponent = 1, point_group = None):
        """ Create a new ImproperRotation about `axis` with order `n` and exponent `exponent`

        :Parameters:

        n : int
            The order of the rotation (e.g. for a C_3 operation, `n`=3)
        axis : Vector
            The axis about which the operation occurs
        exponent : int, optional
            The number of times the operation is applied.  (e.g. for C_3, `exponent` = 1, for C_3^2, `exponent` = 2, etc.)

        """
        self.point_group = point_group
        self.n = n
        self.axis = axis.normalized()
        self.exponent = exponent

    ##############
    # Properties #
    ##############

    @property
    def theta(self):
        """ The angle of rotation, in Radians
        """
        return 2 * self.exponent * math.pi / self.n

    @property
    def matrix(self):
        if self._matrix is not None:
            return self._matrix
        u = self.axis
        cos = math.cos(self.theta)
        sin = math.sin(self.theta)
        self._matrix = Matrix(
            [
                [ cos + u.x**2*(1-cos), u.x*u.y*(1-cos) - u.z*sin, u.x*u.z*(1-cos) + u.y*sin ],
                [ u.y*u.x*(1-cos) + u.z*sin, cos + u.y**2*(1-cos), u.y*u.z*(1-cos) - u.x*sin ],
                [ u.z*u.x*(1-cos) - u.y*sin, u.z*u.y*(1-cos) + u.x*sin, cos + u.z**2*(1-cos) ]
            ]
        )
        return self._matrix


    ###################
    # Special Methods #
    ###################

    def __str__(self):
        if not self.name is None:
            return self.name
        else:
            return 'C' + str(self.n) + (("^" + str(self.exponent)) if not self.exponent == 1 else '')

    def __unicode__(self):
        return u'C' + subscript(self.n) + ((superscript(self.exponent)) if not self.exponent == 1 else '')

    def __repr__(self):
        return 'Rotation(' + str(self.n) + ', ' + repr(self.axis)\
               + ((', ' + str(self.exponent)) if not self.exponent == 1 else '') + ')'

    def __latex__(self):
        return '\ensuremath{\hat{C}_{' + str(self.n) + "}" + (("^{" + str(self.exponent) + "}") if not self.exponent == 1 else '') + '}'

    #def __mul__(self, other):
    #    if not isinstance(other, SymmetryOperation):
    #        raise TypeError("Invalid multiplication of " + classname(self) + " by " + classname(other))
    #    if type(other) is Rotation:
    #        if self.is_same_axis(self.axis, other.axis):
    #            new_n = gcd(self.n, other.n)
    #            new_exp = ((new_n/self.n)*self.exponent + (new_n/other.n)*other.exponent) % new_n
    #            return self.point_group.element_for(Rotation(new_n, self.axis, new_exp))
    #        else:
    #            return super(Rotation, self).__mul__(other)
    #    else:
    #        return super(Rotation, self).__mul__(other)


    #################
    # Class Methods #
    #################

    @classmethod
    def from_matrix(cls, matrix):
        return NotImplemented

    @classmethod
    def about_axis(cls, molecule, axis, max_order):
        ret_val = []
        # Find the largest n...
        max_n = 1
        for n in xrange(2, max_order+1):
            if molecule.has_symmetry(Rotation(n, axis)):
                max_n = n
        for i in range(1, max_n):
            ret_val.append(Rotation(max_n, axis, i))
        return ret_val





