from abc import abstractproperty, abstractmethod
from grendel.coordinates.simple_internal_coordinate import SimpleInternalCoordinate

class PeriodicCoordinate(SimpleInternalCoordinate):
    """
    Base class for all coordinates with periodicity.
    """

    ####################
    # Class Attributes #
    ####################

    preferred_range = abstractproperty()
    #singularity = abstractproperty()

    ##############
    # Properties #
    ##############

    @property
    def period(self):
        return self.preferred_range[1] - self.preferred_range[0]

    ##########################
    # Abstract Class Methods #
    ##########################

    @abstractmethod
    def noncanonical_value_for_xyz(cls, xyz):
        return NotImplemented

    #################
    # Class Methods #
    #################

    @classmethod
    def value_for_xyz(cls, xyz):
        """
        Give the canonical value in the range described by `preferred_range` attribute.
        """
        ret_val = cls.noncanonical_value_for_xyz(xyz)
        return cls.canonicalize_value(ret_val)

    @classmethod
    def possible_values_for_xyz(cls, xyz):
        """
        Give the canonical value as well as the value plus one period and the value minus one period
        (returns a tuple of values in this order)
        """
        period = cls.preferred_range[1] - cls.preferred_range[0]
        val = cls.value_for_xyz(xyz)
        return val, val + period, val - period

    @classmethod
    def canonicalize_value(cls, val):
        period = cls.preferred_range[1] - cls.preferred_range[0]
        if val > cls.preferred_range[1] or val < cls.preferred_range[0]:
            return ((val - cls.preferred_range[0]) % period) + cls.preferred_range[0]
        else:
            return val

