"""
"""
from copy import copy
from grendel.interface.computation import Computation


class Optimization(object):
    """
    """


    ##############
    # Attributes #
    ##############

    steps = None
    starting_molecule = None

    pass


class OptimizationComputation(Optimization, Computation):


    def __new__(cls, **kwargs):
        if 'details' in kwargs:
            kwargs['details'] = copy(kwargs['details'])
            kwargs['details']['optimization'] = True
        # if details is not in kwargs, an error will
        #   be raise in Computation
        return Computation.__new__(cls, **kwargs)
