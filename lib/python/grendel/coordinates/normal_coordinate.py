from __future__ import division
from grendel.coordinates.coordinate import Coordinate
from grendel.gmath import Vector
from grendel.util.decorators import typechecked
from grendel.util.units import *

__all__ = ["NormalCoordinate"]

class NormalCoordinate(Coordinate):
    """
    """

    ####################
    # Class Attributes #
    ####################

    # TODO think through units...I'm not sure this is right
    default_delta = 0.001 * (Angstroms / AMU**(1/2))
    coordinate_symbol = 'q'

    ##############
    # Attributes #
    ##############

    b_vector = None
    frequency = None
    atoms = None
    index = None

    ##################
    # Initialization #
    ##################

    @typechecked(
        b_vector=Vector,
        frequency=(float, complex),
        parent='NormalRepresentation'
    )
    def __init__(self, b_vector, frequency, parent, index, **kwargs):
        self.frequency = frequency
        self.b_vector = b_vector
        self.index = index
        self.atoms = parent.molecule.atoms[:]
        self._init(
            parent=parent,
            **kwargs
        )

    ###################
    # Special Methods #
    ###################

    def __copy_kwargs__(self):
        kw = super(NormalCoordinate, self).__copy_kwargs__()
        kw.update(
            frequency=self.frequency,
            b_vector=self.b_vector
        )

    ###########
    # Methods #
    ###########

    #--------------------------------#
    # Methods abstract in Coordinate #
    #--------------------------------#

    def copy_for_representation(self, rep, **kwargs):
        raise ValueError("NormalCoordinate instances can't be copied to different representations.")

    def value_for_molecule_matrix(self, mat):
        # TODO implement this
        raise NotImplementedError

    @classmethod
    def value_for_positions(cls, *pos):
        raise ValueError("orphaned NormalCoordinates don't have values")



#####################
# Dependent Imports #
#####################

from grendel.representations.normal_representation import NormalRepresentation

