"""
"""
from grendel.gmath.matrix import Matrix
from symmetry_operation import SymmetryOperation

__all__ = [
    'Inversion'
]

class Inversion(SymmetryOperation):
    """ An inversion (i) element of a point group.
    """

    ##############
    # Attributes #
    ##############


    ##################
    # Initialization #
    ##################

    def __init__(self, point_group = None):
        """
        """
        self.point_group = point_group


    ##############
    # Properties #
    ##############

    @property
    def matrix(self):
        """
        """
        if not self._matrix is None:
            return self._matrix
        self._matrix = Matrix([[-1,0,0],[0,-1,0],[0,0,-1]])
        return self._matrix


    ###################
    # Special Methods #
    ###################

    def __str__(self):
        return 'i'

    def __repr__(self):
        return 'Inversion()'


    #################
    # Class Methods #
    #################

    @classmethod
    def exists_for_molecule(cls, molecule):
        return molecule.has_symmetry(Inversion())

