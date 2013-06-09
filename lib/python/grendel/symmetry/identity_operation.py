"""
"""
from grendel.gmath.matrix import Matrix
from grendel.util.strings import classname
from symmetry_operation import SymmetryOperation

__all__ = [
    'IdentityOperation'
]

class IdentityOperation(SymmetryOperation):
    """ The identity (E) element of a point group.
    """

    ##############
    # Attributes #
    ##############


    ##################
    # Initialization #
    ##################

    def __init__(self, point_group = None):
        """ Construct the identity element for the point group `point_group`.
        Raises a `GroupTheoryError` if the point group already has an identity operation.
        """
        self.point_group = point_group


    ##############
    # Properties #
    ##############

    @property
    def matrix(self):
        if self._matrix is not None:
            return self._matrix
        self._matrix = Matrix([
            [ 1., 0., 0.],
            [ 0., 1., 0.],
            [ 0., 0., 1.]
        ])
        return self._matrix


    ###################
    # Special Methods #
    ###################

    def __str__(self):
        return "E"

    def __repr__(self):
        return "IdentityOperation()"

    def __unicode__(self):
        return u'E'

    def __latex__(self):
        return "\ensuremath{\hat{E}}"

    def __mul__(self, other):
        if not isinstance(other, SymmetryOperation):
            raise TypeError("Invalid multiplication of " + classname(self) + " by " + classname(other))
        return other

    def __rmul__(self, other):
        if not isinstance(other, SymmetryOperation):
            raise TypeError("Invalid multiplication of " + classname(self) + " by " + classname(other))
        return other




