""" A module containing a class for conjugacy classes of point group elements.
"""
import itertools
from grendel.util.decorators import CachedProperty
from symmetry_operation import SymmetryOperation

__all__ = [
    'ConjugacyClass'
]

class ConjugacyClass(object):
    """
    """

    ##############
    # Attributes #
    ##############

    elements = []

    ##################
    # Initialization #
    ##################

    def __init__(self, *initial_elements):
        """
        """
        self.elements = []
        for el in initial_elements:
            if not isinstance(el, SymmetryOperation):
                raise TypeError("Only SymmetryOperations can be elements of a ConjugacyClass")
            self.elements.append(el)


    ##############
    # Properties #
    ##############

    @CachedProperty
    def first_element(self):
        self.elements.sort()
        return self.elements[0]

    @property
    def order(self):
        return len(self)

    ###################
    # Special Methods #
    ###################

    def __len__(self):
        return len(self.elements)

    def __iter__(self):
        return itertools.chain(self.elements)

    def __str__(self):
        return "{" + ((str(len(self)) + " ") if len(self) > 1 else '')  + str(self.first_element) + "}"
    __repr__ = __str__  # for now, at least...


    ###########
    # Methods #
    ###########

    def add_element(self, el):
        """ Adds an element to the conjugacy class (if it has not already been added)
        """
        if el not in self:
            self.elements.append(el)



