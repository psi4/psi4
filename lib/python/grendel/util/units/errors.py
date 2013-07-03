from grendel.util.strings import classname

__all__ = [
    'IncompatibleUnitsError',
    'UnknownUnitError',
    'UnitizedObjectError'
]

class IncompatibleUnitsError(Exception):
    """ Exception for attempted incompatible unit conversions
    """

    def __init__(self, unit1, unit2):
        """
        """
        name1 = unit1.name if isinstance(unit1, CompositeUnit) else '<unitless>' if unit1 is None else classname(unit1)
        name2 = unit2.name if isinstance(unit2, CompositeUnit) else '<unitless>' if unit2 is None else classname(unit2)
        message = "Cannot convert from units " + name1 + " to units " + name2 + "."
        super(IncompatibleUnitsError, self).__init__(message)

class UnknownUnitError(Exception):
    """ Exception for attempted conversion of something that is not a subclass of Unit
    """

    def __init__(self, unit1):
        """
        """
        name1 = unit1.name if isinstance(unit1, CompositeUnit) else classname(unit1)
        message = "Unit " + name1 + " is not a Unit that Grendel knows about."
        super(UnknownUnitError, self).__init__(message)

class UnitizedObjectError(ValueError):
    """ For errors encountered in handling unitized objects
    """
    pass

#####################
# Dependent Imports #
#####################

from grendel.util.units.composite import CompositeUnit

