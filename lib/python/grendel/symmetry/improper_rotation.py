"""
"""
from copy import copy
from grendel.gmath import divisors
from grendel.symmetry.reflection import Reflection
from grendel.symmetry.rotation import Rotation
from grendel.util.strings import superscript, subscript
from symmetry_operation import SymmetryOperation

__all__ = [
    'ImproperRotation'
]

class ImproperRotation(SymmetryOperation):
    """ An improper rotation (S_n^m) element of a point group.

    :Attributes:

    n : int
        The order of the rotation (e.g. for a S_3 operation, n=3)
    axis : Vector
        The axis about which the operation occurs
    exponent : int
        The number of times the operation is applied.  (e.g. for S_5, exponent = 1, for S_5^4, exponent = 4, etc.)

    """

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
        """ Create a new ImproperRotation about axis with order n and exponent exponent

        :Parameters:

        n : int
            The order of the rotation (e.g. for a S_3 operation, n=3)
        axis : Vector
            The axis about which the operation occurs
        exponent : int, optional
            The number of times the operation is applied.  (e.g. for S_3, exponent = 1, for S_3^2, exponent = 2, etc.)

        """
        self.point_group = point_group
        self.n = n
        self.axis = axis
        self.exponent = exponent


    ##############
    # Properties #
    ##############

    @property
    def matrix(self):
        if not self._matrix is None:
            return self._matrix
        original_matrix = Rotation(self.n, self.axis).matrix * Reflection(self.axis).matrix
        self._matrix = copy(original_matrix)
        for i in xrange(1, self.exponent):
            self._matrix = self._matrix * original_matrix

        return self._matrix


    ###################
    # Special Methods #
    ###################

    def __str__(self):
        if not self.name is None:
            return self.name
        else:
            return 'S' + str(self.n) + (("^" + str(self.exponent)) if not self.exponent == 1 else '')

    def __unicode__(self):
        return u'S' + subscript(self.n) + ((superscript(self.exponent)) if not self.exponent == 1 else '')

    def __repr__(self):
        return 'ImproperRotation(' + str(self.n) + ', ' + repr(self.axis) \
            + ((', ' + str(self.exponent)) if not self.exponent == 1 else '') + ')'

    def __latex__(self):
        return '$\hat{S}_{' + str(self.n) + "}" + (("^{" + str(self.exponent) + "}") if not self.exponent == 1 else '') + '$'


    #################
    # Class Methods #
    #################

    @classmethod
    def about_axis(cls, molecule, axis, max_ring):
        ret_val = []
        #TODO be more efficient/intelligent about this
        # Find the largest n...
        max_n = 1
        valid_ns = []
        for n in xrange(3, 2 * max_ring+1):
            if molecule.has_symmetry(ImproperRotation(n, axis)):
                valid_ns.append(n)
        for max_n in valid_ns:
            if max_n % 2 == 0:
                for i in xrange(1, max_n):
                    ret_val.append(ImproperRotation(max_n, axis, i))
            else:
                for i in xrange(1, 2*max_n):
                    ret_val.append(ImproperRotation(max_n, axis, i))
        return ret_val


